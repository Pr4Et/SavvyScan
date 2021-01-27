/*
 *  gpubp.cu -- Kernel and C code for CUDA-based backprojection, reprojection
 *               and Fourier filtering
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 2009 by Boulder Laboratory for 3-Dimensional Electron
 *  Microscopy of Cells ("BL3DEMC") and the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 *  $Id$
 */
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cuda_runtime_api.h"
#include "cuda.h"
#include "cufft.h"
#include "b3dutil.h"
#include "imodconfig.h"

#ifdef F77FUNCAP
#define gpuavailable GPUAVAILABLE
#define gpuallocarrays GPUALLOCARRAYS
#define gpuloadproj GPULOADPROJ
#define gpushiftproj GPUSHIFTPROJ
#define gpubpxtilt GPUBPXTILT
#define gpubpnox GPUBPNOX
#define gpubplocal GPUBPLOCAL
#define gpuloadlocals GPULOADLOCALS
#define gpuloadfilter GPULOADFILTER
#define gpufilterlines GPUFILTERLINES
#define gpureproject GPUREPROJECT
#define gpureprojlocal GPUREPROJLOCAL
#define gpureprojoneslice GPUREPROJONESLICE
#define gpudone GPUDONE
#else
#define gpuavailable gpuavailable_
#define gpuallocarrays gpuallocarrays_
#define gpuloadproj gpuloadproj_
#define gpushiftproj gpushiftproj_
#define gpubpxtilt gpubpxtilt_
#define gpubpnox gpubpnox_
#define gpubplocal gpubplocal_
#define gpuloadlocals gpuloadlocals_
#define gpuloadfilter gpuloadfilter_
#define gpufilterlines gpufilterlines_
#define gpureproject gpureproject_
#define gpureprojlocal gpureprojlocal_
#define gpureprojoneslice gpureprojoneslice_
#define gpudone gpudone_
#endif

#if CUDA_VERSION >= 4000
#define HAS_LAYERS
#else
#define cudaArrayLayered 0
#endif

#ifdef __cplusplus
extern "C" {
  int gpuavailable(int *nGPU, float *memory, int *maxTex2D, int *maxTexLayer, 
                   int *maxTex3D, int *debug);
  int gpuallocarrays(int *width, int *nyout, int *nxProjPad, int *nyProj,
                     int *nplanes, int *nviews, int *numWarps, int *numDelz,
                     int *nfilt, int *nreproj, int *firstNpl, int *lastNpl, int *use3D);
  int gpubpnox(float *slice, float *lines, float *sinBeta, float *cosBeta,
               int *nxprj,
               float *xcenIn, float *xcen, float *ycen, float *edgefill);
  int gpushiftproj(int *numPlanes, int *lsliceStart, int *loadStart);
  int gpuloadproj(float *lines, int *numPlanes, int *lsliceStart, 
                  int *loadStart);
  int gpubpxtilt(float *slice, float *sinBeta, float *cosBeta, float *sinAlpha,
                 float *cosAlpha, float *xzfac, float *yzfac, int *nxprj,
                 int *nyProj, float *xcenIn, float *xcen, float *ycen,
                 int *lslice, float *centerSlice, float *edgefill);
  int gpubplocal(float *slice, int *lslice, int *nxwarp, int *nywarp,
                 int *ixswarp, int *iyswarp, int *idxwarp, int *idywarp,
                 int *nxprj, float *xcen, float *xcenIn, float *delxx,
                 float *ycen, float *centerSlice, float *edgefill);
  int gpuloadfilter(float *lines);
  int gpuloadlocals(float *packed, int *numWarps);
  int gpufilterlines(float *lines, int *lslice, int *filterSet);
  int gpureproject(float *lines, float *sinBeta, float *cosBeta, float *sinAlpha, 
                   float *cosAlpha, float *xzfac, float *yzfac, float *delz,
                   int *lsStart, int *lsEnd, int *ithick,
                   float *xcen, float *xcenPdelxx, int *minXreproj, 
                   float *xprjOffset, float *ycen, int *minYreproj,
                   float *yprjOffset, float *centerSlice, int *ifalpha, 
                   float *pmean);
  int gpureprojoneslice(float *slice, float *lines, float *sinBeta, float *cosBeta,
                        float *ycen, int *numproj, float *pmean);
  int gpureprojlocal
  (float *lines, float *sinBeta, float *cosBeta, float *sinAlpha, float *cosAlpha,
   float *xzfac, float *yzfac, int *nxwarp, int *nywarp, int *ixswarp, 
   int *iyswarp, int *idxwarp, int *idywarp, float *warpDelz, int *nWarpDelz, 
   float *dxWarpDelz,float *xprojMin,float *xprojMax, int *lsStart, int *lsEnd,
   int *ithick, int *iview, float *xcen, float *xcenIn, float *delxx, 
   int *minXload, float *xprjOffset, float *ycenAdj, float *yprjOffset,
   float *centerSlice, float *pmean);
  void gpudone();
}
#endif

static int checkProjLoad(int *numPlanes, int *lsliceStart, int startm1);
static int testReportErr(const char *mess);
static int loadBetaInvertCos(float *cosBeta, float *sinBeta, float *costmp,
                             int num);
static int synchronizeCopySlice(float *devslc, int pitch, float *slice,
                                int width, int numLines);
static void pflush(const char *format, ...);
static void pflerr(const char *format, ...);
static void allocerr(const char *mess, int *nplanes, int *firstNpl,
                     int *lastNpl, int ifcuda);



// Offsets to positions in constant array
// For some reason 6 separate arrays did not work for xtilt case
// 7 arrays in 65536 bytes would allow 2340
#define DELTA_OFS  2200
#define MAX_TABLE (6 * DELTA_OFS)
__constant__ float tables[MAX_TABLE];
__constant__ int rpNumz[DELTA_OFS];

#define COSOFS 0
#define SINOFS (1 * DELTA_OFS)
#define CALOFS (2 * DELTA_OFS)
#define SALOFS (3 * DELTA_OFS)
#define XZFOFS (4 * DELTA_OFS)
#define YZFOFS (5 * DELTA_OFS)
#define INVOFS (2 * DELTA_OFS)
#define SINVOFS (3 * DELTA_OFS)

// Definitions for accessing the local alignments arrays with texture calls
#define F11IND 0.f
#define F21IND 1.f
#define F12IND 2.f
#define F22IND 3.f
#define F13IND 4.f
#define F23IND 5.f
#define CAIND 6.f
#define SAIND 7.f
#define CBIND 8.f
#define SBIND 9.f
#define XZFIND 10.f
#define YZFIND 11.f


// declare texture reference for 2D float textures
texture<float, 2, cudaReadModeElementType> projtex2D;
texture<float, 3, cudaReadModeElementType> projtex3D;
#ifdef HAS_LAYERS
texture<float, cudaTextureType2DLayered> projtexLayer;
#endif
texture<float, 2, cudaReadModeElementType> localtex;
texture<float, 2, cudaReadModeElementType> rpSlicetex;
texture<float, 2, cudaReadModeElementType> pfactex;
texture<float, 2, cudaReadModeElementType> delztex;

// Static variables for device arrays
static float *devSlice = NULL;
static cudaArray* devProj = NULL;
static float *devXprojFix = NULL;
static float *devXprojZ = NULL;
static float *devYprojFix = NULL;
static float *devYprojZ = NULL;
static cudaArray *devLocalData = NULL;
static cudaArray *devLocalPfac = NULL;
static cudaArray *devDelz = NULL;
static float *devRadialFilt = NULL;
static float *devFFT = NULL;
static cudaArray *devRpSlice = NULL;
static float *devReproj = NULL;

// Other static variables
static cufftHandle sForwardPlan = 0, sInversePlan = 0;
static int sMaxGflopsDevice = -1;
static int sDeviceSelected = 0;
static size_t sSlicePitch;
static size_t sReprojPitch;
static size_t sLocalPitch;
static int sSliceThick, sSliceWidth, sNumViews, sNumProjPlanes;
static int sLsliceFirst, sNumLoadedPlanes, sNxPlane, sNyPlane, sNumFilts;
static int sCopyFilteredOK = 0;
static int *sPlaneLoaded;
static int sUse3dTexture;

/*
 *  SETUP/SHUTDOWN ROUTINES
 */

/*
 * Test whether a GPU is available, either a GPU of the given number if nGPU is
 * > 0, or the one with the best processing rate if nGPU is 0, and return the
 * memory in bytes.  Return value is 1 for success, 0 for failure.
 */
int gpuavailable(int *nGPU, float *memory, int *maxTex2D, int *maxTexLayer, int *maxTex3D,
                 int *debug)
{
  int current_device = 0;
  int device_count = 0;
  int totalCores, ind;
  float gflops;
  struct cudaDeviceProp device_properties, best_properties;

  // The Mac mini comes through with a clock rate of 0 so allow a 0 product
  float max_gflops = -1.;
  *memory = 0;
  cudaGetDeviceCount( &device_count );
  if (*debug) {
#if CUDA_VERSION >= 3000
    int version, version2;
    cudaRuntimeGetVersion(&version2);
    cudaDriverGetVersion(&version);
    pflush("CUDA version - driver: %d.%02d  runtime: %d.%02d\n", version / 1000,
           version % 1000, version2 / 1000, version2 % 1000);
#endif
    pflush("Device count = %d\n", device_count);
  }
  if (*nGPU != 0) {
    if (*nGPU < 0 || *nGPU > device_count) {
      pflush("The requested GPU number, %d, is out of range; there are only %d devices\n",
             *nGPU, device_count);
      return 0;
    }
    current_device = *nGPU - 1;
    device_count = *nGPU;
  }
  for (; current_device < device_count; current_device++) {
    if (cudaGetDeviceProperties( &device_properties, current_device)
        != cudaSuccess) {
      pflerr("Error returned from trying to get properties of GPU device %d",
               current_device);
      return 0;
    }
    totalCores = totalCudaCores(device_properties.major, device_properties.minor,
                                device_properties.multiProcessorCount);
    if (*debug)
      pflush("Device %d (%s): cores %d  cr %d  major %d minor %d  mem %.0f\n",
             current_device, device_properties.name, 
             totalCores, device_properties.clockRate,
             device_properties.major, device_properties.minor,
             (float)device_properties.totalGlobalMem);
    gflops = totalCores * (float)device_properties.clockRate;

    // This is the best place to catch that the GPU is too old for the version
#ifdef HAS_LAYERS
    if (device_properties.major < 2) {
      if (*debug)
        pflush("The compute level of this GPU is only %d.%d and it cannot be used\n"
               "   with an IMOD package built with CUDA 4 or higher\n", 
               device_properties.major, device_properties.minor);
      gflops = -999.;
    }
#endif

    // Exclude emulation mode (?) which shows up on the Mac
    if( gflops > max_gflops && device_properties.major != 9999) {
      max_gflops = gflops;
      sMaxGflopsDevice = current_device;
      best_properties = device_properties;
    }
  }
    
  if (sMaxGflopsDevice >= 0) {
    *memory = best_properties.totalGlobalMem;
    for (ind = 0; ind < 3; ind++) {
#ifdef HAS_LAYERS
      maxTexLayer[ind] = best_properties.maxTexture2DLayered[ind];
#else
      maxTexLayer[ind] = 0;
#endif
      maxTex3D[ind] = best_properties.maxTexture3D[ind];
      if (ind < 2)
        maxTex2D[ind] = best_properties.maxTexture2D[ind];
    }
    return 1;
  }
  return 0;
}

/*
 * Allocate all needed arrays on the GPU.  Allocate a reconstructed slice or
 * reprojected line array of width x nyout, an array for nplanes of input data 
 * each with nyProj lines of length nxProjPad, and local alignment arrays if
 * numWarps > 0.  If numDelz > 0, this indicates reprojection with local
 * alignments and causes local projection factor arrays to be allocated for
 * nplanes lines, allocation of a CUDA array for those factors too, and 
 * allocation of an array of numDelz x nplanes for warpDelz values.  If 
 * nfilt > 0, also allocate arrays for line filtering with nfilt sets of
 * filters.  If nreproj > 0, allocate separate arrays for reprojecting one
 * slice while still doing regular backprojection arrays.
 */
int gpuallocarrays(int *width, int *nyout, int *nxProjPad, int *nyProj,
                   int *nplanes, int *nviews, int *numWarps, int *numDelz,
                   int *nfilt, int *nreproj, int *firstNpl, int *lastNpl, int *use3D)
{
  size_t pitch1, pitch2, pitch3, memTot;
  int nlines;
  cudaError_t err;

  sUse3dTexture = *use3D;
  sSliceWidth = *width;
  sSliceThick = *nyout;    // Only good for backprojection!
  sNumViews = *nviews;
  sNumProjPlanes = *nplanes;
  sNxPlane = *nxProjPad;
  sNyPlane = *nyProj;
  sNumFilts = *nfilt;

  if (sMaxGflopsDevice < 0)
    return 1;
  if (!sDeviceSelected && cudaSetDevice(sMaxGflopsDevice) != cudaSuccess) {
    allocerr("Error selecting GPU device", nplanes, firstNpl, lastNpl, 1);
    return 1;
  }
  sDeviceSelected = 1;

  if (sNumViews > DELTA_OFS) {
    allocerr("Too many views for the constant memory available on the GPU\n",
             nplanes, firstNpl, lastNpl, 0);
    return 1;
  }

  if (sUse3dTexture < 0 && cudaArrayLayered == 0) {
    pflerr("Bad parameter: texture type -1 was specified, but no layered textures\n"
           "   are available with this version of IMOD\n");
    return 1;
  }

  // Allocate memory for slice or reprojected lines on device
  size_t sizetmp = sSliceWidth * sizeof(float);
  if (cudaMallocPitch((void **)&devSlice, &sSlicePitch, sizetmp, sSliceThick) != 
      cudaSuccess) {
    allocerr("Failed to allocate slice array on GPU device", nplanes, 
             firstNpl, lastNpl, 1);
    return 1;
  }
  //pflush("reproj array size %d %d\n", sSliceWidth, sSliceThick);

  // Allocate memory for projection lines or slices to reproject
  cudaChannelFormatDesc projDesc = cudaCreateChannelDesc
    (32, 0, 0, 0, cudaChannelFormatKindFloat);
  if (!sUse3dTexture)
    err = cudaMallocArray(&devProj, &projDesc, sNxPlane, sNyPlane * sNumProjPlanes);
  else
    err = cudaMalloc3DArray(&devProj, &projDesc, make_cudaExtent 
                            (sNxPlane, sNyPlane, sNumProjPlanes)
#ifdef HAS_LAYERS
                            , sUse3dTexture > 0 ? 0 : cudaArrayLayered
#endif
                            );
  if (err != cudaSuccess) {
    pflush("malloc %d %d %d %d\n", sNxPlane, sNyPlane, sNumProjPlanes,
           sNyPlane * sNumProjPlanes);
    allocerr("Failed to allocate projection array on GPU device", nplanes, 
             firstNpl, lastNpl, 1);
    return 1;
  }
  memTot = sizetmp * sSliceThick + 4 * sNxPlane * sNyPlane * sNumProjPlanes;
  //pflush("input slice array size %d %d %d\n", sNxPlane, sNyPlane, sNumProjPlanes);

  // set texture parameters

  
  // Bind the array to the texture
  if (sUse3dTexture > 0) {
    projtex3D.addressMode[0] = cudaAddressModeClamp;
    projtex3D.addressMode[1] = cudaAddressModeClamp;
    projtex3D.filterMode = cudaFilterModeLinear;
    projtex3D.normalized = false;
    err = cudaBindTextureToArray(projtex3D, devProj, projDesc);
  } else if (sUse3dTexture < 0) {
#ifdef HAS_LAYERS
    projtexLayer.addressMode[0] = cudaAddressModeClamp;
    projtexLayer.addressMode[1] = cudaAddressModeClamp;
    projtexLayer.filterMode = cudaFilterModeLinear;
    projtexLayer.normalized = false;
    err = cudaBindTextureToArray(projtexLayer, devProj, projDesc);
#endif
  } else {
    projtex2D.addressMode[0] = cudaAddressModeClamp;
    projtex2D.addressMode[1] = cudaAddressModeClamp;
    projtex2D.filterMode = cudaFilterModeLinear;
    projtex2D.normalized = false;
    err = cudaBindTextureToArray(projtex2D, devProj, projDesc);
  }
  if (err != cudaSuccess) {
    allocerr("Failed to bind projection array to texture", nplanes, firstNpl, lastNpl, 1);
    return 1;
  }

  if (sNumProjPlanes > 1) {
    sPlaneLoaded = (int *)malloc(sNumProjPlanes * sizeof(int));
    if (!sPlaneLoaded) {
      allocerr("Failed to malloc little array sPlaneLoaded\n", nplanes,
               firstNpl, lastNpl, 0);
      return 1;
    }
  }

  // Get arrays for reprojection of one slice
  if (*nreproj) {
    if (cudaMallocArray(&devRpSlice, &projDesc, sSliceWidth, sSliceThick) !=
        cudaSuccess) {
      allocerr("Failed to allocate slice array for reprojection on GPU device",
               nplanes, firstNpl, lastNpl, 1);
      return 1;
    }
    if (cudaBindTextureToArray(rpSlicetex, devRpSlice, projDesc) != cudaSuccess) {
      allocerr("Failed to bind reprojection slice array to texture", nplanes, 
               firstNpl, lastNpl, 1);
      return 1;
    }
    rpSlicetex.addressMode[0] = cudaAddressModeClamp;
    rpSlicetex.addressMode[1] = cudaAddressModeClamp;
    rpSlicetex.filterMode = cudaFilterModeLinear;
    rpSlicetex.normalized = false;
  
    if (cudaMallocPitch((void **)&devReproj, &sReprojPitch, 
                        (size_t)(sNxPlane * sizeof(float)), *nreproj) != cudaSuccess) {
      allocerr("Failed to allocate reprojected line array on GPU device", 
               nplanes, firstNpl, lastNpl, 1);
      return 1;
    }
    memTot += 4 * sSliceWidth * sSliceThick + sNxPlane * *nreproj;
  }

  // Get arrays for local proj factors
  if (*numWarps > 0) {
    nlines = sNyPlane;

    // Adjust and allocate for reprojection
    if (*numDelz) {
      nlines = sNumProjPlanes;
      sizetmp = sNxPlane * sizeof(float);

      if (cudaMallocArray(&devLocalPfac, &projDesc, sNxPlane, 4 * nlines) != cudaSuccess)
        {
          allocerr("Failed to allocate local factor texture array on GPU device",
                   nplanes, firstNpl, lastNpl, 1);
          return 1;
        }
      //pflush("local factor texture  %d %d\n", sNxPlane, 4 * nlines);
      
      pfactex.filterMode = cudaFilterModePoint;
      pfactex.normalized = false;
      if (cudaBindTextureToArray(pfactex, devLocalPfac, projDesc) != cudaSuccess) {
        allocerr("Failed to bind local factor arrays to texture", nplanes, 
                 firstNpl, lastNpl, 1);
        return 1;
      }
      if (cudaMallocArray(&devDelz, &projDesc, *numDelz, nlines) != cudaSuccess) {
        allocerr("Failed to allocate warpDelz texture array on GPU device",
                 nplanes, firstNpl, lastNpl, 1);
        return 1;
      }
      //pflush("warpdelz texture  %d %d\n", *numDelz, nlines);
      delztex.filterMode = cudaFilterModePoint;
      delztex.normalized = false;
      if (cudaBindTextureToArray(delztex, devDelz, projDesc) != cudaSuccess) {
        allocerr("Failed to bind warpDelz array to texture", nplanes, 
                 firstNpl, lastNpl, 1);
        return 1;
      }
      memTot += 4 * nlines * (4 * sNxPlane + *numDelz);
    }

    // Allocate the arrays always used for local data
    if (cudaMallocPitch((void **)&devXprojFix, &pitch1, sizetmp, nlines) != cudaSuccess ||
        cudaMallocPitch((void **)&devXprojZ, &pitch2, sizetmp, nlines) != cudaSuccess ||
        cudaMallocPitch((void **)&devYprojFix, &pitch3, sizetmp, nlines) != cudaSuccess ||
        cudaMallocPitch((void **)&devYprojZ, &sLocalPitch, sizetmp, nlines) != 
        cudaSuccess  || cudaMallocArray(&devLocalData, &projDesc, *numWarps * sNumViews, 
                                        12) != cudaSuccess) {
      allocerr("Failed to allocate local factor arrays on GPU device", nplanes,
               firstNpl, lastNpl, 1);
      return 1;
    }
    /* pflush("xdevYprojFix pitches  %d %d    localdata %d\n", sNxPlane, nlines,
     *numWarps * sNumViews); */
    if (pitch2 != pitch1 || pitch3 != pitch1 || sLocalPitch != pitch1) {
      allocerr("Array pitches for local GPU arrays do NOT match\n", nplanes,
               firstNpl, lastNpl, 0);
      return 1;
    }

    localtex.filterMode = cudaFilterModePoint;
    localtex.normalized = false;
    if (cudaBindTextureToArray(localtex, devLocalData, projDesc) != cudaSuccess) {
      allocerr("Failed to bind local factor arrays to texture", nplanes, 
               firstNpl, lastNpl, 1);
      return 1;
    }
    memTot += 4 * sizetmp * nlines + 48 * *numWarps * sNumViews;
  }

  // Get arrays for radial filtering
  if (sNumFilts > 0 || sNumProjPlanes > 1) {
    sizetmp = sNxPlane * sNyPlane * sizeof(float);
    if (cudaMalloc((void **)&devFFT, sizetmp)  != cudaSuccess ||
        (sNumFilts > 0 && cudaMalloc((void **)&devRadialFilt, sizetmp * sNumFilts) 
         != cudaSuccess)) {
      allocerr("Failed to allocate GPU arrays for radial filtering", nplanes,
               firstNpl, lastNpl, 1);
      return 1;
    }
    memTot += (1 + sNumFilts) * sizetmp;
  }

  pflush("Allocated %4d MB for arrays (including %d input planes) on the GPU\n"
         , (memTot + 512*1024)/(1024*1024), sNumProjPlanes);
  return 0;
}

// Routine to free all allocated resources
void gpudone()
{
  cudaFree(devSlice);
  cudaFreeArray(devProj);
  cudaFree(devXprojFix);
  cudaFree(devXprojZ);
  cudaFree(devYprojFix);
  cudaFree(devYprojZ);
  cudaFreeArray(devLocalData);
  cudaFreeArray(devLocalPfac);
  cudaFreeArray(devDelz);
  cudaFree(devFFT);
  cudaFree(devRadialFilt);
  cudaFree(devReproj);
  cudaFreeArray(devRpSlice);
  if (sForwardPlan)
    cufftDestroy(sForwardPlan);
  if (sInversePlan)
    cufftDestroy(sInversePlan);
  devSlice = NULL;
  devProj = NULL;
  devXprojFix = NULL;
  devXprojZ = NULL;
  devYprojFix = NULL;
  devYprojZ = NULL;
  devLocalData = NULL;
  devLocalPfac = NULL;
  devDelz = NULL;
  devFFT = NULL;
  devRadialFilt = NULL;
  devReproj = NULL;
  devRpSlice = NULL;
  sForwardPlan = 0;
  sInversePlan = 0;
}

/*
 * ROUTINES FOR LOADING/MAINTAINING STACK OF PLANES ON GPU
 */ 

// Function to shift existing data in preparation for loading new data starting
// in position loadStart (numbered from 1) and with starting slice number
// lsliceStart
int gpushiftproj(int *numPlanes, int *lsliceStart, int *loadStart)
{
  int startm1 = *loadStart - 1;
  int shift, shiftStart, numToShift, todo, dstY, srcY;
  size_t sizetmp = sNxPlane * sizeof(float);
  cudaMemcpy3DParms cpyParms = {0};
  cudaMemcpy3DParms tmpParms = {0};
  if (startm1 > 0) {
    if (checkProjLoad(numPlanes, lsliceStart, startm1))
      return 1;

    // Copy data down without overlap if it goes into occupied planes
    if (startm1 < sNumLoadedPlanes) {
      shift = sNumLoadedPlanes - startm1;
      numToShift = startm1;
      shiftStart = 0;
      if (sUse3dTexture) {
          tmpParms.dstPos = make_cudaPos(0, 0, 0);
          tmpParms.dstPtr = make_cudaPitchedPtr(devFFT, sNxPlane * sizeof(float), 
                                                sNxPlane, sNyPlane);
          tmpParms.srcArray = devProj;
          tmpParms.extent = make_cudaExtent(sNxPlane, sNyPlane, 1);
          tmpParms.kind = cudaMemcpyDeviceToDevice;
          cpyParms.srcPos = make_cudaPos(0, 0, 0);
          cpyParms.srcPtr = make_cudaPitchedPtr(devFFT, sNxPlane * sizeof(float),
                                                sNxPlane, sNyPlane);
          cpyParms.dstArray = devProj;
          cpyParms.extent = make_cudaExtent(sNxPlane, sNyPlane, 1);
          cpyParms.kind = cudaMemcpyDeviceToDevice;
      }

      // Loop on the planes or sets of planes to copy
      while (numToShift > 0) {
        if (!sUse3dTexture) {
          todo = shift;
          if (todo > numToShift)
            todo = numToShift;
          dstY = shiftStart * sNyPlane;
          srcY = dstY + shift * sNyPlane;
          //pflush("Copying down %d\n", todo);
          if (cudaMemcpy2DArrayToArray(devProj, 0, dstY, devProj, 0, srcY,
                                       sizetmp, todo * sNyPlane,
                                       cudaMemcpyDeviceToDevice) != cudaSuccess){
            pflerr("Error copying segment of projection array down");
            sNumLoadedPlanes = 0;
            return 1;
          }
        } else {

          // Sadly this can only copy one plane at a time at least for the layered array
          todo = 1;
          tmpParms.srcPos = make_cudaPos(0, 0, shiftStart + shift);
          if (cudaMemcpy3D(&tmpParms) != cudaSuccess) {
            pflerr("Error copying plane %d of projection array to devFFT", 
                   shiftStart + shift);
            sNumLoadedPlanes = 0;
            return 1;
          }
          cpyParms.dstPos = make_cudaPos(0, 0, shiftStart);
          if (cudaMemcpy3D(&cpyParms) != cudaSuccess) {
            pflerr("Error copying devFFT to plane %d of projection array", 
                   shiftStart);
            sNumLoadedPlanes = 0;
            return 1;
          }
        }
        numToShift -= todo;
        shiftStart += todo;
      }
    }
  }
  sNumLoadedPlanes = startm1;
  sLsliceFirst = *lsliceStart - startm1;

  /*pflush("Initializing array num %d  first %d  loaded %d\n", sNumProjPlanes, 
    sLsliceFirst, sNumLoadedPlanes); */
  // Initialize array for keeping track of copied planes, and enable copying
  for (todo = 0; todo < sNumProjPlanes; todo++)
    sPlaneLoaded[todo] = todo < sNumLoadedPlanes ? 1 : 0;
  sCopyFilteredOK = 1;
  return 0;
}

// Function to load numPlanes planes of input data, starting in position
// loadStart (numbered from 1) and with starting slice number lsliceStart
int gpuloadproj(float *lines, int *numPlanes, int *lsliceStart, int *loadStart)
{
  int startm1 = *loadStart - 1;
  int todo, dstY, numCopy = 0;
  cudaMemcpy3DParms cpyParms = {0};
  cudaError_t err;

  if (startm1 > 0 && checkProjLoad(numPlanes, lsliceStart, startm1)) {
    sCopyFilteredOK = 0;
    return 1;
  }

  // Check for valid load
  if (startm1 + *numPlanes > sNumProjPlanes) {
    pflush("Trying to load past end of projection array\n");
    sCopyFilteredOK = 0;
    sNumLoadedPlanes = 0;
    return 1;
  }
  
  // Find the number to copy by the last plane not already loaded
  if (sCopyFilteredOK) {
    for (todo = startm1; todo < startm1 + *numPlanes; todo++)
      if (!sPlaneLoaded[todo])
        numCopy = todo + 1 - startm1;
  }
  sCopyFilteredOK = 0;

  // Finally do the load
  //if (numCopy) pflush("Loading %d planes\n", numCopy);
  if (numCopy) {
    if (sUse3dTexture) {
    
      cpyParms.srcPos = make_cudaPos(0, 0, 0);
      cpyParms.dstPos = make_cudaPos(0, 0, startm1);
      cpyParms.srcPtr = make_cudaPitchedPtr(lines, sNxPlane * sizeof(float), sNxPlane, 
                                            sNyPlane);
      cpyParms.dstArray = devProj;
      cpyParms.extent = make_cudaExtent(sNxPlane, sNyPlane, numCopy);
      cpyParms.kind = cudaMemcpyHostToDevice;
      err = cudaMemcpy3D(&cpyParms);
    } else {
      dstY = startm1 * sNyPlane;
      todo = numCopy * sNyPlane * sNxPlane * 4;
      err = cudaMemcpyToArray(devProj, 0, dstY, lines, todo, cudaMemcpyHostToDevice);
    }
    if (err != cudaSuccess) {
      pflerr("Failed to copy projection array to device");
      sNumLoadedPlanes = 0;
      return 1;
    }
  }
  sNumLoadedPlanes = startm1 + *numPlanes;
  sLsliceFirst = *lsliceStart - startm1;
  return 0;
}

// Function to do initial check on parameters in load/shift calls
static int checkProjLoad(int *numPlanes, int *lsliceStart, int startm1)
{
  if (!sNumLoadedPlanes) {
    pflush("Trying to load into higher planes when none are loaded\n");
    return 1;
  }
  if (sLsliceFirst + sNumLoadedPlanes != *lsliceStart) {
    pflush("Starting slice %d does not match first slice %d + num loaded %d"
            "\n", *lsliceStart, sLsliceFirst, sNumLoadedPlanes);
    sNumLoadedPlanes = 0;
    return 1;
  }
  if (startm1 > sNumLoadedPlanes) {
    pflush("Starting plane %d is past number loaded %d\n", startm1+1, 
           sNumLoadedPlanes);
    sNumLoadedPlanes = 0;
    return 1;
  }
  return 0;
}

/*
 * ROUTINES FOR RADIAL FILTERING OF INPUT LINES
 */

// Kernel to multiply the FFT by the filter
__global__ void filterFFT(float *FFT, float *filter, int nxProjPad, int nviews, 
                          float scale)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  if (i < nviews && j < nxProjPad)
    FFT[i * nxProjPad + j] *= filter[i * nxProjPad + j] * scale;
}

// Function to load the filter lines into the array and generate plans
int gpuloadfilter(float *lines)
{
  size_t sizetmp = sNxPlane * sNumViews * sNumFilts * sizeof(float);
  if (cudaMemcpy(devRadialFilt, lines, sizetmp, cudaMemcpyHostToDevice) !=
      cudaSuccess) {
    pflerr("Failed to copy radial filters to GPU array");
    gpudone();
    return 1;
  }
  if (cufftPlan1d(&sForwardPlan, sNxPlane - 2, CUFFT_R2C, sNumViews) != 
      CUFFT_SUCCESS || cufftPlan1d(&sInversePlan, sNxPlane - 2, CUFFT_C2R, 
                                   sNumViews) != CUFFT_SUCCESS) {
    pflush("Failed to generate a plan for CUFFT\n");
    gpudone();
    return 1;
  }
  return 0;
}

// Function to filter the set of input lines
int gpufilterlines(float *lines, int *lslice, int *filterSet)
{
  int ind, blockX = 16;
  size_t sizetmp = sNxPlane * sNumViews * sizeof(float);
  float scale = 1.f / (sNxPlane - 2);
  cudaMemcpy3DParms cpyParms = {0};
  cudaError_t err;
  if (cudaMemcpy(devFFT, lines, sizetmp, cudaMemcpyHostToDevice) !=
      cudaSuccess) {
    pflerr("Failed to copy lines to GPU array for radial filtering");
    return 1;
  }
  if (cufftExecR2C(sForwardPlan, devFFT, (cufftComplex *)devFFT) != 
      CUFFT_SUCCESS) {
    pflush("Failure in forward FFT on GPU\n");
    return 1;
  }
  
  // Filter!!!
  dim3 blockSize(blockX, 16, 1);
  dim3 gridSize((sNxPlane + blockSize.x - 1) / blockSize.x, 
                (sNumViews + blockSize.y - 1) / blockSize.y, 1);

  filterFFT<<<gridSize, blockSize>>>
    (devFFT, devRadialFilt + (*filterSet - 1) * sNxPlane * sNumViews, sNxPlane, 
     sNumViews, scale);
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    pflerr("Error executing threads for filtering"); 
    return 1;
  }
  if (cudaThreadSynchronize() != cudaSuccess) {
    pflerr("Error return from synchronizing after filtering");
    return 1;
  }

  if (cufftExecC2R(sInversePlan, (cufftComplex *)devFFT, devFFT) != 
      CUFFT_SUCCESS) {
    pflush("Failure in inverse FFT on GPU\n");
    return 1;
  }
  if (cudaMemcpy(lines, devFFT, sizetmp, cudaMemcpyDeviceToHost) !=
      cudaSuccess) {
    pflerr("Failed to copy radial filtered lines back from GPU array");
    return 1;
  }
  
  // If copying is OK and it is a slice in needed range, copy it to proj
  if (sCopyFilteredOK) {
    ind = *lslice - sLsliceFirst;
    if (ind >= 0 && ind < sNumProjPlanes) {
      //pflush("Copying %d to plane %d\n", *lslice,ind);
      if (sUse3dTexture) {
        cpyParms.srcPos = make_cudaPos(0, 0, 0);
        cpyParms.dstPos = make_cudaPos(0, 0, ind);
        cpyParms.srcPtr = make_cudaPitchedPtr(devFFT, sNxPlane * sizeof(float), sNxPlane,
                                              sNyPlane);
        cpyParms.dstArray = devProj;
        cpyParms.extent = make_cudaExtent(sNxPlane, sNyPlane, 1);
        cpyParms.kind = cudaMemcpyDeviceToDevice;
        if (cudaMemcpy3D(&cpyParms) == cudaSuccess)
          sPlaneLoaded[ind] = 1;
      } else {
        if (cudaMemcpyToArray(devProj, 0, ind * sNumViews, devFFT, sizetmp,
                            cudaMemcpyDeviceToDevice) == cudaSuccess)
          sPlaneLoaded[ind] = 1;
      }
    }
  }
  return 0;
}

/*
 * ROUTINES FOR SIMPLE BACK-PROJECTION (NO X-AXIS TILT, ETC)
 */

// Kernel for simple back-projection with testing at ends of lines
__global__ void bpNoXtTest(float *slice, int pitch, int jbase, int iwide,
                             int nxProj, int ithick, int nviews, 
                             float xcenIn, float xcenOut, float ycenOut, 
                             float edgefill)
{
  float cosBeta, sinBeta, zpart, kproj, xp;
  float sum = 0.;
  int iv;
  int j = blockIdx.x * blockDim.x + threadIdx.x + jbase;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  if (j < iwide && i < ithick) {
    for (iv = 0; iv < nviews; iv++) {
      cosBeta = tables[iv+COSOFS];
      sinBeta = tables[iv+SINOFS];
      zpart = (i + 1 - ycenOut) * sinBeta + xcenIn;
      xp =  zpart + (j + 1 - xcenOut) * cosBeta - 0.5f;
      //if (i == 150) printf("%d %d %d  %d  %d  %.2f  %.2f\n", iv, j, i, jlft, jrt, zpart + (1 - xcenOut) * cosBeta - 0.5f, zpart + (nxProj - xcenOut) * cosBeta - 0.5f);
      if (xp >= 0.5 && xp <= nxProj - 0.5) {
        kproj = iv + 0.5f;
        sum += tex2D(projtex2D, xp, kproj);
      } else {
        sum += edgefill;
      }
    }
    slice[i * pitch + j] = sum;
  }
}

// Kernel for simple back-projection with no testing
__global__ void bpNoXtFast(float *slice, int pitch, int jbase, int iwide,
                             int ithick, int nviews, 
                             float xcenIn, float xcenOut, float ycenOut)
{
  float cosBeta, sinBeta, zpart, kproj, xp;
  float sum = 0.;
  int iv;
  int j = blockIdx.x * blockDim.x + threadIdx.x + jbase;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  if (i < ithick) {
    for (iv = 0; iv < nviews; iv++) {
      cosBeta = tables[iv+COSOFS];
      sinBeta = tables[iv+SINOFS];
      zpart = (i + 1 - ycenOut) * sinBeta + xcenIn;
      kproj = iv + 0.5f;
      xp =  zpart + (j + 1 - xcenOut) * cosBeta - 0.5f;
      sum += tex2D(projtex2D, xp, kproj);
    }
    slice[i * pitch + j] = sum;
  }
}

// Function to run simple backprojection
int gpubpnox(float *slice, float *lines, float *sinBeta, float *cosBeta,
             int *nxProj, float *xcenIn, float *xcenOut, float *ycenOut,
             float *edgefill)
{
  size_t sizetmp = sizeof(float) * sNxPlane * sNumViews;
  float cosinv[DELTA_OFS];
  int iv, jlft, jrt, jlftmax, jrtmin, gridLeft, gridFast, gridRight;
  float zpart, xlft, xrt, xlfttmp, xrttmp;
  int i, blockX = 16;

  if (sUse3dTexture) {
    pflush("GPU problem: Program called gpubpnox after setting up to use texture type "
           "%d\n", sUse3dTexture);
    return 1;
  }
  if (loadBetaInvertCos(cosBeta, sinBeta, cosinv, sNumViews))
    return 1;

  // Copy projections
  if (cudaMemcpyToArray(devProj, 0, 0, lines, sizetmp, cudaMemcpyHostToDevice)
      != cudaSuccess) {
    pflerr("Failed to copy projection array to device");
    return 1;
  }

  // Find limits of region that needs no testing
  jlftmax = 1;
  jrtmin = sSliceWidth;
  for (iv = 0; iv < sNumViews; iv++) {
    for (i = 0; i <= sSliceThick - 1; i += sSliceThick - 1) {
      zpart = (i + 1 - *ycenOut) * sinBeta[iv] + *xcenIn;
      xlfttmp = (1. - zpart) * cosinv[iv] + *xcenOut;
      xrttmp = (*nxProj - zpart) * cosinv[iv] + *xcenOut;
      xlft = fmin(xlfttmp, xrttmp);
      xrt = fmax(xlfttmp, xrttmp);
      jlft = (int)ceilf(xlft);
      jrt = (int)ceilf(xrt) - 1;
      jlftmax = max(jlftmax, jlft);
      jrtmin = min(jrtmin, jrt);
      //printf("%d %f %d %.2f %d  %d  %.2f  %.2f\n", iv, cbet, i, zpart, jlft, jrt, zpart + (1 - *xcenOut) * cbet - 0.5f, zpart + (*nxProj - *xcenOut) * cbet - 0.5f);
    }
  }

  // Figure out grid sizes for left test, fast, and right test regions
  dim3 blockSize(blockX, 16, 1);
  dim3 gridSize((sSliceWidth + blockSize.x - 1) / blockSize.x, 
                (sSliceThick + blockSize.y - 1) / blockSize.y, 1);

  gridLeft = (jlftmax - 1 + blockX - 1) / blockX;
  gridFast = jrtmin / blockX - gridLeft;
  if (gridFast <= 0) {
    gridLeft = gridSize.x;
    gridRight = 0;
  } else
    gridRight = gridSize.x - (gridFast + gridLeft);

  if (gridLeft > 0) {
    gridSize.x = gridLeft;
    bpNoXtTest<<<gridSize, blockSize>>>
      (devSlice, sSlicePitch / 4, 0, sSliceWidth, *nxProj, sSliceThick, 
       sNumViews, *xcenIn, *xcenOut, *ycenOut, *edgefill);
    if (testReportErr("in left test region of backprojection"))
      return 1;
  }

  if (gridFast > 0) {
    gridSize.x = gridFast;
    bpNoXtFast<<<gridSize, blockSize>>>
      (devSlice, sSlicePitch / 4, blockX * gridLeft, sSliceWidth,
       sSliceThick, sNumViews, *xcenIn, *xcenOut, *ycenOut);
    if (testReportErr("in no-test region of backprojection"))
      return 1;
  }

  if (gridRight > 0) {
    gridSize.x = gridRight;
    bpNoXtTest<<<gridSize, blockSize>>>
      (devSlice, sSlicePitch / 4, blockX * (gridLeft + gridFast), sSliceWidth, 
       *nxProj, sSliceThick, sNumViews, *xcenIn, *xcenOut, *ycenOut, *edgefill);
    if (testReportErr("in right test region of backprojection"))
      return 1;
  }

  return (synchronizeCopySlice(devSlice, sSlicePitch, slice, sSliceWidth, 
                               sSliceThick));
    
}

/*
 * ROUTINES FOR BACK-PROJECTION WITH X AXIS TILT AND/OR Z FACTORS
 */

// Kernel for BP with X-axis tilt/Z-factors and testing at ends of lines
#define BPXTTEST_START \
    for (iv = 0; iv < nviews; iv++) {    \
      cosBeta = tables[iv+COSOFS];         \
      sinBeta = tables[iv+SINOFS];         \
      sinAlpha = tables[iv+SALOFS];         \
      cosAlpha = tables[iv+CALOFS];         \
      zpart = yy * sinAlpha * sinBeta +          \
        zz * (cosAlpha * sinBeta + tables[iv+XZFOFS]) + xcenIn;         \
      yproj = yy * cosAlpha - zz * (sinAlpha - tables[iv+YZFOFS]) + centerSlice;   \
      xp =  zpart + xx * cosBeta - 0.5f;         \
      if (yproj >= 1. - ytol && yproj <= nyProj + ytol && xp >= 0.5 &&          \
          xp < nxProj - 0.5) {         \
        yproj = fmax(1.f, fmin((float)nyProj, yproj));         \
        jproj = min((int)yproj, nyProj - 1);         \
        fj = yproj - jproj;

#define BPXTTEST_END           \
      } else {         \
        sum += edgefill;         \
      }         \
    }

__global__ void bpXtiltTest(float *slice, int pitch, int jbase, int iwide,
                            int nxProj, int nyProj, int ithick, int nviews, 
                            float xcenIn, float xcenOut, float ycenOut, float yy,
                            float centerSlice, int lsliceBase, float edgefill, int use3D)
{
  float cosBeta, sinBeta, zpart, kproj, xp, zz, cosAlpha, sinAlpha, fj, yproj, xx;
  float sum = 0.;
  int iv, jproj;
#ifdef HAS_LAYERS
  int jslice;
#endif
  float ytol = 3.05f;
  int j = blockIdx.x * blockDim.x + threadIdx.x + jbase;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  if (j < iwide && i < ithick) {
    zz = (i + 1 - ycenOut);
    xx = (j + 1 - xcenOut);
#ifdef HAS_LAYERS
    if (use3D == 0) {
#endif
      BPXTTEST_START;
      kproj = (jproj - lsliceBase) * nviews + iv + 0.5f;
      sum += (1.f - fj) * tex2D(projtex2D, xp, kproj) + 
        fj * tex2D(projtex2D, xp, kproj + nviews);
      BPXTTEST_END;
#ifdef HAS_LAYERS
    } else {
      BPXTTEST_START;
      jslice = jproj - lsliceBase;
      kproj = iv + 0.5f;
      sum += (1.f - fj) * tex2DLayered(projtexLayer, xp, kproj, jslice) + 
        fj * tex2DLayered(projtexLayer, xp, kproj, jslice + 1);
      BPXTTEST_END;
    }
#endif
    slice[i * pitch + j] = sum;
  }
}

// Kernel for BP with X-axis tilt/Z-factors and no testing 
#define BPXTFAST_ALL \
      cosBeta = tables[iv+COSOFS];         \
      sinBeta = tables[iv+SINOFS];         \
      sinAlpha = tables[iv+SALOFS];         \
      cosAlpha = tables[iv+CALOFS];         \
      zpart = yy * sinAlpha * sinBeta +          \
        zz * (cosAlpha * sinBeta + tables[iv+XZFOFS]) + xcenIn;         \
      yproj = yy * cosAlpha - zz * (sinAlpha - tables[iv+YZFOFS]) + centerSlice;     \
      jproj = (int)yproj;         \
      fj = yproj - jproj;         \
      xp =  zpart + xx * cosBeta - 0.5f;

__global__ void bpXtiltFast(float *slice, int pitch, int jbase, int iwide, int ithick,
                            int nviews, float xcenIn, float xcenOut, float ycenOut,
                            float yy, float centerSlice, int lsliceBase, int use3D)
{
  float cosBeta, sinBeta, zpart, kproj, xp, zz, cosAlpha, sinAlpha, fj, yproj, xx;
  float sum = 0.;
  int iv, jproj;
#ifdef HAS_LAYERS
  int jslice;
#endif
  int j = blockIdx.x * blockDim.x + threadIdx.x + jbase;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  if (i < ithick) {
    zz = (i + 1 - ycenOut);
    xx = (j + 1 - xcenOut);
#ifdef HAS_LAYERS
    if (use3D == 0) {
#endif
      for (iv = 0; iv < nviews; iv++) {
        BPXTFAST_ALL;
        kproj = (jproj - lsliceBase) * nviews + iv + 0.5f;
        sum += (1.f - fj) * tex2D(projtex2D, xp, kproj) + 
          fj * tex2D(projtex2D, xp, kproj + nviews);
      }
#ifdef HAS_LAYERS
    } else {
      for (iv = 0; iv < nviews; iv++) {
        BPXTFAST_ALL;
        jslice = jproj - lsliceBase;
        kproj = iv + 0.5f;
        sum += (1.f - fj) * tex2DLayered(projtexLayer, xp, kproj, jslice) + 
          fj * tex2DLayered(projtexLayer, xp, kproj, jslice + 1);
      }
    }
#endif
    slice[i * pitch + j] = sum;
  }
}

// Function to run back-projection with X-axis tilt/Z-factors
int gpubpxtilt(float *slice, float *sinBeta, float *cosBeta, 
               float *sinAlpha, float *cosAlpha, float *xzfac, float *yzfac,
               int *nxProj, int *nyProj, float *xcenIn, float *xcenOut, float *ycenOut,
               int *lslice, float *centerSlice, float *edgefill)
{
  int iv, jlft, jrt, jlftmax, jrtmin, gridLeft, gridFast, gridRight;
  float zpart, xlft, xrt, xlfttmp, xrttmp, yy, zz, yproj;
  int i, ytest, blockX = 16;
  float cosinv[DELTA_OFS];

  if (sUse3dTexture > 0) {
    pflush("GPU problem: Program called gpubpxtilt after setting up to use 3D textures"
           "\n"); 
    return 1;
  }

  if (loadBetaInvertCos(cosBeta, sinBeta, cosinv, sNumViews))
    return 1;

  // Copy alphas and z factors
  iv = sNumViews * sizeof(float);
  if (cudaMemcpyToSymbol(tables, cosAlpha, iv, CALOFS*4, cudaMemcpyHostToDevice)
      || cudaMemcpyToSymbol(tables, sinAlpha, iv, SALOFS*4,
                            cudaMemcpyHostToDevice) ||
      cudaMemcpyToSymbol(tables, xzfac, iv, XZFOFS*4, cudaMemcpyHostToDevice)
      || cudaMemcpyToSymbol(tables, yzfac, iv, YZFOFS*4,
                            cudaMemcpyHostToDevice)) {
    pflerr("Failed to copy constant data to GPU");
    return 1;
  }

  // Find limits of region that needs no testing.  Test every angle top & bot
  jlftmax = 1;
  jrtmin = sSliceWidth;
  yy = *lslice - *centerSlice;
  ytest = 0;
  for (iv = 0; iv < sNumViews; iv++) {
    for (i = 0; i <= sSliceThick - 1; i += sSliceThick - 1) {
      zz = (i + 1 - *ycenOut);
      zpart = yy * sinAlpha[iv] * sinBeta[iv] + zz * (cosAlpha[iv] * sinBeta[iv] +
                                                  xzfac[iv]) + *xcenIn;
      yproj = yy * cosAlpha[iv] - zz * (sinAlpha[iv] - yzfac[iv]) + *centerSlice;
      if (yproj < 1 || yproj > *nyProj - 1)
        ytest = 1;
      xlfttmp = (1. - zpart) * cosinv[iv] + *xcenOut;
      xrttmp = (*nxProj - zpart) * cosinv[iv] + *xcenOut;
      xlft = fmin(xlfttmp, xrttmp);
      xrt = fmax(xlfttmp, xrttmp);
      jlft = (int)ceilf(xlft);
      jrt = (int)ceilf(xrt) - 1;
      jlftmax = max(jlftmax, jlft);
      jrtmin = min(jrtmin, jrt);
      //printf("%d %f %d %.2f %d  %d  %.2f  %.2f\n", iv, cbet, i, zpart, jlft, jrt, zpart + (1 - *xcenOut) * cbet - 0.5f, zpart + (*nxProj - *xcenOut) * cbet - 0.5f);
    }
  }

  // Figure out grid sizes for left test, fast, and right test regions
  dim3 blockSize(blockX, 16, 1);
  dim3 gridSize((sSliceWidth + blockSize.x - 1) / blockSize.x, 
                (sSliceThick + blockSize.y - 1) / blockSize.y, 1);

  gridLeft = (jlftmax - 1 + blockX - 1) / blockX;
  gridFast = jrtmin / blockX - gridLeft;
  if (gridFast <= 0 || ytest) {
    gridLeft = gridSize.x;
    gridRight = 0;
    gridFast = 0;
  } else
    gridRight = gridSize.x - (gridFast + gridLeft);

  if (gridLeft > 0) {
    gridSize.x = gridLeft;
    bpXtiltTest<<<gridSize, blockSize>>>
      (devSlice, sSlicePitch / 4, 0, sSliceWidth, *nxProj, *nyProj, sSliceThick, 
       sNumViews, *xcenIn, *xcenOut, *ycenOut, yy, *centerSlice, sLsliceFirst, *edgefill,
       sUse3dTexture);
    if (testReportErr("in left test region of backprojection"))
      return 1;
  }

  if (gridFast > 0) {
    gridSize.x = gridFast;
    bpXtiltFast<<<gridSize, blockSize>>>
      (devSlice, sSlicePitch / 4, blockX * gridLeft, sSliceWidth, sSliceThick, sNumViews,
       *xcenIn, *xcenOut, *ycenOut, yy, *centerSlice, sLsliceFirst, sUse3dTexture);
    if (testReportErr("in no-test region of backprojection"))
      return 1;
  }

  if (gridRight > 0) {
    gridSize.x = gridRight;
    bpXtiltTest<<<gridSize, blockSize>>>
      (devSlice, sSlicePitch / 4, blockX * (gridLeft + gridFast), sSliceWidth, *nxProj,
       *nyProj, sSliceThick, sNumViews, *xcenIn, *xcenOut, *ycenOut, yy, *centerSlice, 
       sLsliceFirst, *edgefill, sUse3dTexture);
    if (testReportErr("in right test region of backprojection"))
      return 1;
  }

  return (synchronizeCopySlice(devSlice, sSlicePitch, slice, sSliceWidth,
                               sSliceThick));
}

/*
 * ROUTINES FOR BACK-PROJECTION WITH LOCAL ALIGNMENTS
 */

// Kernel for back-projection using local projection factors, testing as needed
#define BPLOCAL_START         \
      for (iv = 0; iv < nviews; iv++) {         \
        ind = iv * sLocalPitch + j;         \
        xp = xprojf[ind] + zz * xprojz[ind] - 0.5f;         \
        yproj = yprojf[ind] + zz * yprojz[ind];         \
        if (yproj >= lsliceBase - ytol && yproj <= lsliceLast + ytol &&          \
            xp >= 0.5f && xp < nxProj - 0.5f) {         \
          yproj = fmax((float)lsliceBase, fmin((float)lsliceLast, yproj));         \

#define BPLOCAL_END         \
      } else {         \
        sum += edgeFill;         \
      }         \
    }

__global__ void bpLocalTest(float *slice, int slPitch, float *xprojf, 
                            float *xprojz, float *yprojf, float *yprojz, 
                            int sLocalPitch, int iwide,
                            int nxProj, int lsliceLast, int ithick, int nviews,
                            float ycenOut, int lsliceBase, float edgeFill, int use3D)
{
  float kproj, xp, zz, fj, yproj, baseAdj;
  float sum = 0.;
  float ytol = 3.05f;
  int iv, jproj, ind;
#ifdef HAS_LAYERS
  int jslice;
#endif
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  if (i < ithick && j < iwide) {
    zz = (i + 1 - ycenOut);
    if (use3D > 0) {
      baseAdj = (float)lsliceBase - 0.5f;
      BPLOCAL_START;
      kproj = iv + 0.5f;
      sum += tex3D(projtex3D, xp, kproj, yproj - baseAdj);
      BPLOCAL_END;
#ifdef HAS_LAYERS
    } else if (use3D < 0) {
      BPLOCAL_START;
      kproj = iv + 0.5f;
      jproj = min((int)yproj, lsliceLast - 1);
      fj = yproj - jproj;
      jslice = jproj - lsliceBase;
      sum += (1.f - fj) * tex2DLayered(projtexLayer, xp, kproj, jslice) + 
        fj * tex2DLayered(projtexLayer, xp, kproj, jslice + 1);
      BPLOCAL_END;
    } else {
#endif
      BPLOCAL_START;
      jproj = min((int)yproj, lsliceLast - 1);
      fj = yproj - jproj;
      kproj = (jproj - lsliceBase) * nviews + iv + 0.5f;
      sum += (1.f - fj) * tex2D(projtex2D, xp, kproj) + 
        fj * tex2D(projtex2D, xp, kproj + nviews);
      BPLOCAL_END;
      }       
    slice[i * slPitch + j] = sum;
  }
}

// Kernel for computing the local projection factors from warping data
__global__ void localProjFactors
(float *xprjf, float *xprjz, float *yprjf, float *yprjz, int pitch, int iv, 
 int nviews, int iwide, int minX, int lslice, int nlines, int nxWarp, int nyWarp,
 int ixStartWarp, int iyStartWarp, int iDelXwarp, int iDelYwarp, float xcenOut,
 float xcenIn, float xcenPaxisOfs, float centerSlice)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int line = blockIdx.y * blockDim.y + threadIdx.y;
  int ind1, ind2, ind3, ind4, ixc, ixt, ixpos, iyt, iypos;
  float fnd1, fnd2, fnd3, fnd4, yzf1, yzf2, yzf3, yzf4, jpos;
  float f1, f2, f3, f4, xx, yy, fx, fy;
  float calf, salf, a11, a12, a21, a22, xadd, yadd, xalladd, yalladd;
  float calf2, salf2, a112, a122, a212, a222, xadd2, yadd2;
  float calf3, salf3, a113, a123, a213, a223, xadd3, yadd3;
  float calf4, salf4, a114, a124, a214, a224, xadd4, yadd4;
  float f1x, f2x, f3x, f4x, f1xy, f2xy, f3xy, f4xy;
  float f1y, f2y, f3y, f4y, f1yy, f2yy, f3yy, f4yy;
  float xp1f, xp1z, yp1f, xp2f, xp2z, yp2f, xp3f, xp3z, yp3f, xp4f, xp4z, yp4f;
  float cosBeta, sinBeta, cosBeta2, sinBeta2, cosBeta3, sinBeta3, cosBeta4, sinBeta4;

  if (j >= iwide || line >= nlines)
    return;
  if (iv < 0)
    iv = line;
  else
    lslice += line;

  // Need to add 1 to j when it is used as a position
  jpos = j + minX + 1;
  ixc = (int)floor(jpos - xcenOut + xcenPaxisOfs + 0.5f);
  ixt = min(max(ixc - ixStartWarp, 0), (nxWarp - 1) * iDelXwarp);
  ixpos = min(ixt / iDelXwarp + 1, nxWarp - 1);
  fx = ((float)(ixt - (ixpos - 1) * iDelXwarp)) / iDelXwarp;
  iyt = min(max(lslice - iyStartWarp, 0), (nyWarp - 1) * iDelYwarp);
  iypos = min(iyt / iDelYwarp + 1, nyWarp - 1);
  fy = ((float)(iyt - (iypos - 1) * iDelYwarp)) / iDelYwarp;

  ind1 = (nxWarp * (iypos - 1) + ixpos - 1) * nviews + iv;
  ind2 = ind1 + nviews;
  ind3 = ind1 + nxWarp * nviews;
  ind4 = ind3 + nviews;
  f1 = (1. - fy) * (1. - fx);
  f2 = (1. - fy) * fx;
  f3 = fy * (1. - fx);
  f4 = fy * fx;
  fnd1 = ind1;
  fnd2 = ind2;
  fnd3 = ind3;
  fnd4 = ind4;
  
  cosBeta = tex2D(localtex, fnd1, CBIND);
  sinBeta = tex2D(localtex, fnd1, SBIND);
  calf = tex2D(localtex, fnd1, CAIND);
  salf = tex2D(localtex, fnd1, SAIND);
  a11 = tex2D(localtex, fnd1, F11IND);
  a12 = tex2D(localtex, fnd1, F12IND);
  a21 = tex2D(localtex, fnd1, F21IND);
  a22 = tex2D(localtex, fnd1, F22IND);
  xadd = tex2D(localtex, fnd1, F13IND) + xcenIn - xcenIn * a11 - centerSlice * a12;
  yadd = tex2D(localtex, fnd1, F23IND) + centerSlice - xcenIn * a21 - centerSlice * a22;

  cosBeta2 = tex2D(localtex, fnd2, CBIND);
  sinBeta2 = tex2D(localtex, fnd2, SBIND);
  calf2 = tex2D(localtex, fnd2, CAIND);
  salf2 = tex2D(localtex, fnd2, SAIND);
  a112 = tex2D(localtex, fnd2, F11IND);
  a122 = tex2D(localtex, fnd2, F12IND);
  a212 = tex2D(localtex, fnd2, F21IND);
  a222 = tex2D(localtex, fnd2, F22IND);
  xadd2 = tex2D(localtex, fnd2, F13IND) + xcenIn - xcenIn * a112 - centerSlice * a122;
  yadd2 = tex2D(localtex, fnd2, F23IND) + centerSlice - xcenIn * a212 - 
    centerSlice * a222;

  cosBeta3 = tex2D(localtex, fnd3, CBIND);
  sinBeta3 = tex2D(localtex, fnd3, SBIND);
  calf3 = tex2D(localtex, fnd3, CAIND);
  salf3 = tex2D(localtex, fnd3, SAIND);
  a113 = tex2D(localtex, fnd3, F11IND);
  a123 = tex2D(localtex, fnd3, F12IND);
  a213 = tex2D(localtex, fnd3, F21IND);
  a223 = tex2D(localtex, fnd3, F22IND);
  xadd3 = tex2D(localtex, fnd3, F13IND) + xcenIn - xcenIn * a113 - centerSlice * a123;
  yadd3 = tex2D(localtex, fnd3, F23IND) + centerSlice - xcenIn * a213 - 
    centerSlice * a223;

  cosBeta4 = tex2D(localtex, fnd4, CBIND);
  sinBeta4 = tex2D(localtex, fnd4, SBIND);
  calf4 = tex2D(localtex, fnd4, CAIND);
  salf4 = tex2D(localtex, fnd4, SAIND);
  a114 = tex2D(localtex, fnd4, F11IND);
  a124 = tex2D(localtex, fnd4, F12IND);
  a214 = tex2D(localtex, fnd4, F21IND);
  a224 = tex2D(localtex, fnd4, F22IND);
  xadd4 = tex2D(localtex, fnd4, F13IND) + xcenIn - xcenIn * a114 - centerSlice * a124;
  yadd4 = tex2D(localtex, fnd4, F23IND) + centerSlice - xcenIn * a214 - 
    centerSlice * a224;
       
  f1x = f1 * a11;
  f2x = f2 * a112;
  f3x = f3 * a113;
  f4x = f4 * a114;
  f1xy = f1 * a12;
  f2xy = f2 * a122;
  f3xy = f3 * a123;
  f4xy = f4 * a124;

  f1y = f1 * a21;
  f2y = f2 * a212;
  f3y = f3 * a213;
  f4y = f4 * a214;
  f1yy = f1 * a22;
  f2yy = f2 * a222;
  f3yy = f3 * a223;
  f4yy = f4 * a224;

  xalladd = f1 * xadd + f2 * xadd2 + f3 * xadd3 + f4 * xadd4;
  yalladd = f1 * yadd + f2 * yadd2 + f3 * yadd3 + f4 * yadd4;
       
  // Each projection position is a sum of a fixed factor ("..f")
  // and a factor that multiplies z ("..z")
   
  xx = jpos - xcenOut;
  yy = lslice - centerSlice;
  xp1f = xx * cosBeta + yy * salf * sinBeta + xcenPaxisOfs;
  xp1z = calf * sinBeta + tex2D(localtex, fnd1, XZFIND);
  xp2f = xx * cosBeta2 + yy * salf2 * sinBeta2 + xcenPaxisOfs;
  xp2z = calf2 * sinBeta2 + tex2D(localtex, fnd2, XZFIND);
  xp3f = xx * cosBeta3 + yy * salf3 * sinBeta3 + xcenPaxisOfs;
  xp3z = calf3 * sinBeta3 + tex2D(localtex, fnd3, XZFIND);
  xp4f = xx * cosBeta4 + yy * salf4 * sinBeta4 + xcenPaxisOfs;
  xp4z = calf4 * sinBeta4 + tex2D(localtex, fnd4, XZFIND);

  yp1f = yy * calf + centerSlice;
  yp2f = yy * calf2 + centerSlice;
  yp3f = yy * calf3 + centerSlice;
  yp4f = yy * calf4 + centerSlice;

  // store the fixed and z - dependent component of the
  // projection coordinates
  yzf1 = tex2D(localtex, fnd1, YZFIND);
  yzf2 = tex2D(localtex, fnd2, YZFIND);
  yzf3 = tex2D(localtex, fnd3, YZFIND);
  yzf4 = tex2D(localtex, fnd4, YZFIND);
  ind1 = pitch * line + j;
  xprjf[ind1] = f1x * xp1f + f2x * xp2f + f3x * xp3f + f4x * xp4f + 
    f1xy * yp1f + f2xy * yp2f + f3xy * yp3f + f4xy * yp4f + xalladd;
  xprjz[ind1] = f1x * xp1z + f2x * xp2z + f3x * xp3z + f4x * xp4z - 
    (f1xy * (salf - yzf1) + f2xy * (salf2 - yzf2) + f3xy * (salf3 - yzf3) + 
     f4xy * (salf4 - yzf4));
  yprjf[ind1] = f1y * xp1f + f2y * xp2f + f3y * xp3f + f4y * xp4f + 
    f1yy * yp1f + f2yy * yp2f + f3yy * yp3f + f4yy * yp4f + yalladd;
  yprjz[ind1] = f1y * xp1z + f2y * xp2z + f3y * xp3z + f4y * xp4z - 
    (f1yy * (salf - yzf1) + f2yy * (salf2 - yzf2) + f3yy * (salf3 - yzf3) + 
     f4yy * (salf4 - yzf4));
}

// Function to load the local alignment data
int gpuloadlocals(float *packed, int *numWarps)
{
  size_t sizetmp = sizeof(float) * *numWarps * sNumViews * 12;
  if (cudaMemcpyToArray(devLocalData, 0, 0, packed, sizetmp,
                        cudaMemcpyHostToDevice) != cudaSuccess) {
    pflerr("Failed to copy local data to GPU array");
    gpudone();
    return 1;
  }
  return 0;
}

// Function to run back-projection with local alignments, first computing the
// the projection factors for all positions and views, then running the 
// back projection kernel
int gpubplocal(float *slice, int *lslice, int *nxWarp, int *nyWarp,
               int *ixStartWarp, int *iyStartWarp, int *iDelXwarp, int *iDelYwarp,
               int *nxProj, float *xcenOut, float *xcenIn, float *axisXoffset,
               float *ycenOut, float *centerSlice, float *edgefill)
{
  int blockX = 16;

  // Compute the local projection factors
  dim3 blockFac(blockX, 16, 1);
  dim3 gridFac((sSliceWidth + blockFac.x - 1) / blockFac.x, 
                (sNumViews + blockFac.y - 1) / blockFac.y, 1);
  localProjFactors<<<gridFac, blockFac>>>
    (devXprojFix, devXprojZ, devYprojFix, devYprojZ, sLocalPitch / 4, -1, sNumViews, 
     sSliceWidth, 0, *lslice, sNumViews, *nxWarp, *nyWarp, *ixStartWarp, *iyStartWarp, 
     *iDelXwarp, *iDelYwarp, *xcenOut, *xcenIn, *xcenIn + *axisXoffset, *centerSlice);
  if (testReportErr("computing localProjFactors"))
      return 1;

  if (cudaThreadSynchronize() != cudaSuccess) {
    pflerr("Error return from synchronizing after computing local factors");
    return 1;
  }

  // Do the backprojection
  dim3 blockSize(blockX, 16, 1);
  dim3 gridSize((sSliceWidth + blockSize.x - 1) / blockSize.x, 
                (sSliceThick + blockSize.y - 1) / blockSize.y, 1);

  bpLocalTest<<<gridSize, blockSize>>>
    (devSlice, sSlicePitch / 4, devXprojFix, devXprojZ, devYprojFix, devYprojZ, 
     sLocalPitch / 4, sSliceWidth, *nxProj, sLsliceFirst + sNumLoadedPlanes - 1, 
     sSliceThick, sNumViews, *ycenOut, sLsliceFirst, *edgefill, sUse3dTexture);
  if (testReportErr("for local backprojection"))
      return 1;

  return (synchronizeCopySlice(devSlice, sSlicePitch, slice, sSliceWidth, 
                               sSliceThick));
}

/*
 * ROUTINES FOR REPROJECTION
 */

// Kernel to do simple reprojection (no X axis tilt or Z factors)
__global__ void reprojNox(float *lines, int pitch, int iwide, int ithick, 
                          int lsliceStart, int lsliceEnd, int lsliceBase, 
                          float xxlim, float xcenAdj, float xcenPaxisOfs,
                          float xProjOffset, float ycenAdj, float sinBeta,
                          float cbetinv, float delz, int numz, float pmean)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  int line, kz;
  float zz, sum, frac, zslice, xproj, xx;
  line = i + lsliceStart;
  if (j >= iwide || line > lsliceEnd)
    return;
  sum = 0.;
  xproj = j + 1 + xProjOffset;
  for (kz = 0; kz < numz; kz++) {
    zz = 1 + kz * delz;
    frac = 1.;
    if (zz > ithick) {
      frac = 1. - (zz - (int)zz);
      zz = ithick;
    }
    zslice = zz - 0.5f;
    zz -= ycenAdj;

    // the usual -0.5 is incorporated into xcenAdj
    xx = (xproj - (zz  * sinBeta + xcenPaxisOfs)) * cbetinv + xcenAdj;
    if (xx < 0.5f || xx > xxlim) {
      sum += frac * pmean;
    } else {
      zslice += (line - lsliceBase) * ithick;
      sum += frac * tex2D(projtex2D, xx, zslice);
    }
  }
  lines[pitch * i + j] = sum;
}

// Kernel to do reprojection with X axis tilt and/or Z factors
__global__ void reprojXtilt
(float *lines, int pitch, int iwide, int ithick, int lsliceStart, int lsliceEnd, 
 int lsliceBase, int lsliceLast, float xxlim, float xcenAdj, float xcenPaxisOfs, 
 float xProjOffset, float centerSlice, float yProjOffset, float ycenAdj, float cbetinv, 
 float calfinv, float salfmyz, float salfsbet, float calsbetpxz, float delz, int numz, 
 float pmean, int use3D)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  int line, kz, iys;
#ifdef HAS_LAYERS
  int lslice;
#endif
  float zz, sum, frac, zslice, yproj, yy, yslice, xproj, xx, fy, ysliceUse, baseAdj;
  float ytol = 3.05f;
  line = i + lsliceStart;
  if (j >= iwide || line > lsliceEnd)
    return;
  sum = 0.;
  xproj = j + 1 + xProjOffset;
  yproj = line + yProjOffset;
  baseAdj = (float)lsliceBase - 0.5f;
  for (kz = 0; kz < numz; kz++) {
    zz = 1 + kz * delz;
    frac = 1.;
    if (zz > ithick) {
      frac = 1. - (zz - (int)zz);
      zz = ithick;
    }
    zslice = zz - 0.5f;
    zz -= ycenAdj;
    yy = (yproj + zz * salfmyz - centerSlice) * calfinv;
    yslice = yy + centerSlice - yProjOffset;

    // the usual -0.5 is incorporated into xcenAdj
    xx = (xproj - (yy * salfsbet + zz * calsbetpxz + xcenPaxisOfs)) * cbetinv + xcenAdj;
    if (xx < 0.5f || xx > xxlim || yslice < lsliceBase - ytol ||
        yslice > lsliceLast + ytol) {
      sum += frac * pmean;
    } else if (use3D > 0) {
      ysliceUse = fmin((float)lsliceLast, fmax((float)lsliceBase, yslice)) - baseAdj;
      sum += frac * tex3D(projtex3D, xx, zslice, ysliceUse);
    } else {
      iys = (int)yslice;
      if (iys < lsliceBase) {
        iys = lsliceBase;
        fy = 0.;
      } else if (iys >= lsliceLast) {
        iys = lsliceLast - 1;
        fy = 1.;
      } else {
        fy = yslice - iys;
      }
#ifdef HAS_LAYERS
      if (use3D < 0) {
        lslice = iys - lsliceBase;
        sum += frac * ((1. - fy) * tex2DLayered(projtexLayer, xx, zslice, lslice) + 
                       fy * tex2DLayered(projtexLayer, xx, zslice, lslice + 1));
      } else {
#endif
        zslice += (iys - lsliceBase) * ithick;
        sum += frac * ((1. - fy) * tex2D(projtex2D, xx, zslice) + 
                       fy * tex2D(projtex2D, xx, zslice + ithick));
#ifdef HAS_LAYERS
      }
#endif
    }
  }
  lines[pitch * i + j] = sum;
}

// Kernel to do simple reprojection at high angles (no X axis tilt or Z factors) 
__global__ void reprojNoxHigh
(float *lines, int pitch, int iwide, int ithick, int lsliceStart, int lsliceEnd, 
 int lsliceBase, float zzlim, float xcenAdj, float xcenPaxisOfs, float xProjOffset, 
 float ycenAdj, float cosBeta, float denomInv, float delx, int numx, float pmean)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  int line, kx;
  float zz, sum, frac, zslice, xproj, xx;
  line = i + lsliceStart;
  if (j >= iwide || line > lsliceEnd)
    return;
  sum = 0.;
  xproj = j + 1 + xProjOffset;
  for (kx = 0; kx < numx; kx++) {
    xx = 1.f + kx * delx;
    frac = 1.f;
    if (xx > iwide) {
      frac = 1.f - (xx - (int)xx);
      xx = iwide;
    }
    
    zz = (xproj - xcenPaxisOfs - (xx - xcenAdj) * cosBeta) * denomInv;
    zslice = zz + ycenAdj;

    if (zslice < 0.5f || zslice > zzlim) {
      sum += frac * pmean;
    } else {
      zslice += (line - lsliceBase) * ithick;
      sum += frac * tex2D(projtex2D, xx - 0.5f, zslice);
    }
  }
  lines[pitch * i + j] = sum;
}

// Kernel to do reprojection at high angles with X axis tilt and/or Z factors
__global__ void reprojXtiltHigh
(float *lines, int pitch, int iwide, int ithick, int lsliceStart, int lsliceEnd, 
 int lsliceBase, int lsliceLast, float zzlim, float xcenAdj, float xcenPaxisOfs, 
 float xProjOffset, float centerSlice, float yProjOffset, float ycenAdj, float cosBeta, 
 float calfinv, float salfmyz, float salfsbetdcal, float denomInv, float delx, int numx, 
 float pmean, int use3D)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  int line, kx, iys;
#ifdef HAS_LAYERS
  int lslice;
#endif
  float zz, sum, frac, zslice, yproj, yy, yslice, xproj, xx, fy, ysliceUse, baseAdj;
  float ytol = 3.05f;
  line = i + lsliceStart;
  if (j >= iwide || line > lsliceEnd)
    return;
  sum = 0.;
  xproj = j + 1 + xProjOffset;
  yproj = line + yProjOffset;
  baseAdj = (float)lsliceBase - 0.5f;
  for (kx = 0; kx < numx; kx++) {
    xx = 1.f + kx * delx;
    frac = 1.f;
    if (xx > iwide) {
      frac = 1.f - (xx - (int)xx);
      xx = iwide;
    }
    
    // Here xcenAdj does not have the -0.5 and ycenAdj does
    zz = (xproj - (yproj - centerSlice) * salfsbetdcal - xcenPaxisOfs - 
          (xx - xcenAdj) * cosBeta) * denomInv;
    yy = (yproj + zz * salfmyz - centerSlice) * calfinv;
    yslice = yy + centerSlice - yProjOffset;
    zslice = zz + ycenAdj;

    if (zslice < 0.5f || zslice > zzlim || yslice < lsliceBase - ytol ||
        yslice > lsliceLast + ytol) {
      sum += frac * pmean;
    } else {
      xx -= 0.5f;
      if (use3D > 0) {
        ysliceUse = fmin((float)lsliceLast, fmax((float)lsliceBase, yslice)) - baseAdj;
        sum += frac * tex3D(projtex3D, xx, zslice, ysliceUse);
      } else {
        iys = (int)yslice;
        if (iys < lsliceBase) {
          iys = lsliceBase;
          fy = 0.;
        } else if (iys >= lsliceLast) {
          iys = lsliceLast - 1;
          fy = 1.;
        } else {
          fy = yslice - iys;
        }
#ifdef HAS_LAYERS
        if (use3D < 0) {
          lslice = iys - lsliceBase;
          sum += frac * ((1. - fy) * tex2DLayered(projtexLayer, xx, zslice, lslice) + 
                         fy * tex2DLayered(projtexLayer, xx, zslice, lslice + 1));
        } else {
#endif
          zslice += (iys - lsliceBase) * ithick;
          sum += frac * ((1. - fy) * tex2D(projtex2D, xx, zslice) + 
                         fy * tex2D(projtex2D, xx, zslice + ithick));
#ifdef HAS_LAYERS
        }
#endif
      }
    }
  }
  lines[pitch * i + j] = sum;
}

// Function to run reprojection for all cases except local alignments
int gpureproject(float *lines, float *sinBeta, float *cosBeta, float *sinAlpha, 
                 float *cosAlpha, float *xzfac, float *yzfac, float *delz,
                 int *lsliceStart, int *lsliceEnd, int *ithick,
                 float *xcenOut, float *xcenPaxisOfs, int *minXreproj, 
                 float *xProjOffset, float *ycenOut, int *minYreproj,
                 float *yProjOffset, float *centerSlice, int *ifAlpha, float *pmean)
{ 
  int blockX = 16;
  int numz, numx, numLines = *lsliceEnd + 1 - *lsliceStart;
  int lastSlice = sLsliceFirst + sNumLoadedPlanes - 1;
  float znum, xcenAdj, salfsbet, calsbetpxz, ycenAdj, salfmyz, cbetinv,calfinv;
  float delx, xnum, salsbetdcal, denomInv;

  dim3 blockSize(blockX, 16, 1);
  dim3 gridSize((sSliceWidth + blockSize.x - 1) / blockSize.x, 
                (numLines + blockSize.y - 1) / blockSize.y, 1);

  if (*ifAlpha == 0 && sUse3dTexture != 0) {
    pflush("GPU problem: Program called gpureproject after setting up to use texture "
           "type %d\n", sUse3dTexture);
    return 1;
  }

  // Common items
  xcenAdj = *xcenOut - (*minXreproj-1) - 0.5;
  ycenAdj = *ycenOut + 1 - *minYreproj;
  salfmyz = *sinAlpha - *yzfac;
  calfinv = 1. / *cosAlpha;
  calsbetpxz = *cosAlpha * *sinBeta + *xzfac;
  
  if (fabs(*sinBeta * *ithick) <= fabs(*cosBeta * sSliceWidth)) {

    // Regular low-angle lines
    znum = 1. + (*ithick - 1) / *delz;
    numz = (int)znum;
    if (znum - numz > 0.1)
      numz++;
    salfsbet = *sinAlpha * *sinBeta;
    cbetinv = 1. / *cosBeta;

    if (*ifAlpha) {
      reprojXtilt<<<gridSize, blockSize>>>
        (devSlice, sSlicePitch / 4, sSliceWidth, *ithick, *lsliceStart, *lsliceEnd, 
         sLsliceFirst, lastSlice, sNxPlane - 0.5, xcenAdj, *xcenPaxisOfs,
         *xProjOffset, *centerSlice, *yProjOffset, ycenAdj, cbetinv, calfinv, salfmyz,
         salfsbet, calsbetpxz, *delz, numz, *pmean, sUse3dTexture);
    } else {
      reprojNox<<<gridSize, blockSize>>>
        (devSlice, sSlicePitch / 4, sSliceWidth, *ithick, *lsliceStart, *lsliceEnd, 
         sLsliceFirst, sNxPlane - 0.5, xcenAdj, *xcenPaxisOfs, *xProjOffset,
         ycenAdj, *sinBeta, cbetinv, *delz, numz, *pmean);
    }
  } else {

    // High angle vertical lines
    delx = (float)fabs(*sinBeta);
    xnum = 1. + (sSliceWidth - 1) / delx;
    numx = (int)xnum;
    if (xnum - numx > 0.1)
      numx++;
    salsbetdcal = *sinAlpha * *sinBeta / *cosAlpha;
    denomInv = 1. / (salfmyz * salsbetdcal + calsbetpxz);
    if (*ifAlpha) {
      reprojXtiltHigh<<<gridSize, blockSize>>>
        (devSlice, sSlicePitch / 4, sSliceWidth, *ithick, *lsliceStart, *lsliceEnd, 
         sLsliceFirst, lastSlice, *ithick - 0.5, xcenAdj + 0.5, *xcenPaxisOfs,
         *xProjOffset, *centerSlice, *yProjOffset, ycenAdj - 0.5, *cosBeta, calfinv, 
         salfmyz, salsbetdcal, denomInv, delx, numx, *pmean, sUse3dTexture);
    } else {
      reprojNoxHigh<<<gridSize, blockSize>>>
        (devSlice, sSlicePitch / 4, sSliceWidth, *ithick, *lsliceStart, *lsliceEnd, 
         sLsliceFirst, *ithick - 0.5, xcenAdj + 0.5, *xcenPaxisOfs, *xProjOffset,
         ycenAdj - 0.5, *cosBeta, denomInv, delx, numx, *pmean);
    }
  }
  if (testReportErr("for reprojection"))
    return 1;
  return (synchronizeCopySlice(devSlice, sSlicePitch, lines, sSliceWidth,
                               numLines));
}

/*
 * ROUTINES TO REPROJECT A SINGLE SLICE
 */

// Kernel to reproject one slice
__global__ void reprojOneSlice(float *lines, int pitch, int iwide, int ithick, 
                               float ycen, int numproj, float pmean)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  int kz;
  float zz, sum, frac, xcenAdj, xx;
  sum = 0.;
  if (j >= iwide || i >= numproj )
    return;
  for (kz = 0; kz < rpNumz[i]; kz++) {
    zz = 1 + kz * tables[COSOFS + i];
    frac = 1.;
    if (zz > ithick) {
      frac = 1. - (zz - (int)zz);
      zz = ithick;
    }
    xcenAdj = iwide / 2;

    // Invert what is multipled by sine because these sines were never inverted
    // inside tilt.f, unlike the signs for regular reproj
    // The usual 0.5 is incorporated into xcenAdj
    xx = (j + 1 - ((ycen - zz)  * tables[SINOFS+i] + xcenAdj + 0.5f)) * 
      tables[INVOFS+i] + xcenAdj;
    if (xx < 0.5f || xx > iwide - 0.5f) {
      sum += frac * pmean;
    } else {
      sum += frac * tex2D(rpSlicetex, xx, zz - 0.5f);
    }
  }
  lines[pitch * i + j] = sum;
}

// Kernel to reproject one slice at high angles
__global__ void reprojOneHighSlice(float *lines, int pitch, int iwide, int ithick, 
                                   float ycen, int numproj, float pmean)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  int kz;
  float zz, sum, frac, xcenAdj, xx, delta;
  if (j >= iwide || i >= numproj )
    return;
  sum = 0.f;
  if (rpNumz[i] >= 0) {
    xcenAdj = iwide / 2;
    delta = fabsf(tables[COSOFS + i]);
    for (kz = 0; kz < rpNumz[i]; kz++) {
      zz = 1 + kz * delta;
      frac = 1.f;
      if (zz > ithick) {
        frac = 1.f - (zz - (int)zz);
        zz = ithick;
      }
      
    // Invert what is multipled by sine because these sines were never inverted
    // inside tilt.f, unlike the signs for regular reproj
    // The usual -0.5 is incorporated into xcenAdj
      xx = (j + 1 - ((ycen - zz)  * tables[SINOFS+i] + xcenAdj + 0.5f)) * 
        tables[INVOFS+i] + xcenAdj;
      if (xx < 0.5f || xx > iwide - 0.5f) {
        sum += frac * pmean;
      } else {
        sum += frac * tex2D(rpSlicetex, xx, zz - 0.5f);
      }
    }
  } else {

    // Going across in X.  Here ycen incorporates the -0.5 and xcenAdj does not
    xcenAdj = iwide / 2 + 0.5f;
    ycen -= 0.5;
    delta = fabsf(tables[SINOFS + i]);
    for (kz = 0; kz < -rpNumz[i]; kz++) {
      xx = 1 + kz * delta;
      frac = 1.f;
      if (xx > iwide) {
        frac = 1.f - (xx - (int)xx);
        xx = iwide;
      }
      
      zz = ((xx - xcenAdj) * tables[COSOFS + i] - j - 1 + xcenAdj) * tables[SINVOFS + i] +
        ycen;
      if (zz < 0.5f || zz > ithick - 0.5f) {
        sum += frac * pmean;
      } else {
        sum += frac * tex2D(rpSlicetex, xx - 0.5f, zz);
      }
    }
  }
  lines[pitch * i + j] = sum;
}

// Function to reproject one slice
int gpureprojoneslice(float *slice, float *lines, float *sinBeta, float *cosBeta,
                      float *ycen, int *numproj, float *pmean)
{
  float znum, cosinv[DELTA_OFS], sininv[DELTA_OFS];
  int numz[DELTA_OFS];
  int blockX = 16;
  int iv, high = 0;

  // Get limited inverse cosines and number of points to do in Z
  loadBetaInvertCos(cosBeta, sinBeta, cosinv, *numproj);

  for (iv = 0; iv < *numproj; iv++) {
    if (fabs(sinBeta[iv] * sSliceThick) <= fabs(cosBeta[iv] * sSliceWidth)) {
      znum = 1. + (sSliceThick - 1) * fabs(cosinv[iv]);
      numz[iv] = (int)znum;
      if (znum - numz[iv] > 0.1)
        numz[iv]++;
      sininv[iv] = 0.;
    } else {

      // For high angle slice, get the number of columns in X, save as a negative
      high = 1;
      sininv[iv] = 1. / sinBeta[iv];
      znum = 1. + (sSliceWidth - 1) * fabs(sininv[iv]);
      numz[iv] = (int)znum;
      if (znum - numz[iv] > 0.1)
        numz[iv]++;
      numz[iv] = -numz[iv];
    }
  }

  // Load constant data
  iv = *numproj * sizeof(float);
  if (cudaMemcpyToSymbol(tables, cosinv, iv, INVOFS*4, cudaMemcpyHostToDevice)
      || cudaMemcpyToSymbol(rpNumz, numz, iv, 0, cudaMemcpyHostToDevice) ||
      (high && cudaMemcpyToSymbol(tables, sininv, iv, SINVOFS*4, 
                                  cudaMemcpyHostToDevice))) {
    pflerr("Failed to copy constant data to GPU");
    return 1;
  }
  
  // Copy slice
  iv = sizeof(float) * sSliceWidth * sSliceThick;
  if (cudaMemcpyToArray(devRpSlice, 0, 0, slice, iv, cudaMemcpyHostToDevice)
      != cudaSuccess) {
    pflerr("Failed to copy slice array to device");
    return 1;
  }
  dim3 blockSize(blockX, 16, 1);
  dim3 gridSize((sSliceWidth + blockSize.x - 1) / blockSize.x, 
                (*numproj + blockSize.y - 1) / blockSize.y, 1);
  if (high)
    reprojOneHighSlice<<<gridSize, blockSize>>>
      (devReproj, sReprojPitch / 4, sSliceWidth, sSliceThick, *ycen, *numproj, *pmean);
  else
    reprojOneSlice<<<gridSize, blockSize>>>
      (devReproj, sReprojPitch / 4, sSliceWidth, sSliceThick, *ycen, *numproj, *pmean);

  if (testReportErr("for reprojection"))
    return 1;

  return (synchronizeCopySlice(devReproj, sReprojPitch, lines, sNxPlane, *numproj));
}

/*
 * ROUTINES FOR REPROJECTION WITH LOCAL ALIGNMENTS
 */

/*
  Finds loaded point that projects to xproj, yproj at centered Z value
  zz, using stored values for [xy]zfac[fv].  Takes starting value in xx,yy
  and returns found value.
  Xproj, yproj are coordinates in original aligned stack.
  XX coordinate is in terms of the loaded data in X
  YY coordinate is in yterms of slices of reconstruction
*/
__device__ void loadedProjectingPoint
(float xproj, float yproj, float zz, float ofsxpz, float ofsypf, float ofsypz, 
 int nxload, int lsliceBase, int lsliceLast, float *xx, float *yy)
{
  int iter, ix, iy, ifout;
  float xp11, yp11, xp12, yp12, xp21, yp21, xerr, yerr, dypx, dxpy,dxpx;
  float dypy, den, fx, fy, findx1, findx2, findy1, findy2;

  for (iter = 0; iter < 5; iter++) {
    ix = (int)floor(*xx);
    iy = (int)floor(*yy);
    ifout = 0;
    if (ix < 1 || ix >= nxload || iy < lsliceBase || iy >= lsliceLast) {
      ifout = 1;
      ix = min(nxload - 1, max(1, ix));
      iy = min(lsliceLast - 1, max(lsliceBase, iy));
    }

    findx1 = ix - 1;
    findx2 = findx1 + 1.;
    findy1 = iy - lsliceBase;
    findy2 = findy1 + 1;
    //*yy = tex2D(pfactex, findx1, findy1 + ofsypf); return;
    xp11 = tex2D(pfactex, findx1, findy1) + 
      tex2D(pfactex, findx1, findy1 + ofsxpz) * zz;
    yp11 = tex2D(pfactex, findx1, findy1 + ofsypf) + 
      tex2D(pfactex, findx1, findy1 + ofsypz) * zz;
    xp21 = tex2D(pfactex, findx2, findy1) + 
      tex2D(pfactex, findx2, findy1 + ofsxpz) * zz;
    yp21 = tex2D(pfactex, findx2, findy1 + ofsypf) + 
      tex2D(pfactex, findx2, findy1 + ofsypz) * zz;
    xp12 = tex2D(pfactex, findx1, findy2) + 
      tex2D(pfactex, findx1, findy2 + ofsxpz) * zz;
    yp12 = tex2D(pfactex, findx1, findy2 + ofsypf) + 
      tex2D(pfactex, findx1, findy2 + ofsypz) * zz;
 
    xerr = xproj - xp11;
    yerr = yproj - yp11;
    dxpx = xp21 - xp11;
    dxpy = xp12 - xp11;
    dypx = yp21 - yp11;
    dypy = yp12 - yp11;
    den = dxpx * dypy - dxpy * dypx;
    fx = (xerr * dypy - yerr * dxpy) / den;
    fy = (dxpx * yerr - dypx * xerr) / den;
    *xx = ix + fx;
    *yy = iy + fy;
    if (fx > -0.1 & fx < 1.1 && fy > -0.1 && fy < 1.1) 
      return;
    if (ifout && (iter > 0 ||  *xx < 0. || *xx > nxload + 1 || 
                  *yy < lsliceBase - 1. || *yy > lsliceLast + 1.))
      return;
  }
}

// Kernel for reprojection with local alignments
__global__ void reprojLocal
(float *lines, int pitch, int nWarpDelz, float dxWarpDelz, int nxload,
 int iwide, int ithick, int lsliceStart, int lsliceEnd, int lsliceBase, int lsliceLast,
 float xprojMin, float xprojMax, float xcenAdj, float xcenPaxisOfs,
 float xProjOffset, float centerSlice, float yProjOffset, float ycenAdj, float cosBeta,
 float sinBeta, float cbetinv, float calfinv, float salfmyz, float salfsbet,
 float calsbetpxz, float pmean, int use3D)
{
  int j = blockIdx.x * blockDim.x + threadIdx.x;
  int i = blockIdx.y * blockDim.y + threadIdx.y;
  int line, lastZdone, iy;
#ifdef HAS_LAYERS
  int lslice;
#endif
  float zz, sum, frac, zslice, yproj, yy, xproj, xx, fy, zind, fline, ofsypz;
  float xxtex, ofsxpz, ofsypf, baseAdj;
  float ytol = 3.05f;
  float zzlim, lbaseMtol, llastPtol, dxWarpInv;
  //  int skip =390;

  line = i + lsliceStart;
  sum = 0.;
  if (j >= iwide || line > lsliceEnd)
    return;

  ofsxpz = lsliceLast + 1 - lsliceBase;
  ofsypf = ofsxpz + ofsxpz;
  ofsypz = ofsypf + ofsxpz;
  fline = i;
  yproj = line + yProjOffset;
  baseAdj = (float)lsliceBase - 0.5f;

  /* Get x projection coord, starting centered Z coordinate, and
     approximate x and y coordinates 
     X coordinate needs to be a loaded X index
     Y coordinate is in slices of reconstruction */

  // ycenAdj needs to be ycen - (minYreproj - 1)
  // xcenAdj = xcen - (minXload - 1)
  xproj = j + 1 + xProjOffset;
  zz = 1. - ycenAdj;
  yy = (yproj + zz * salfmyz - centerSlice) * calfinv + centerSlice;
  xx = (xproj - (yy*salfsbet + zz * calsbetpxz + xcenPaxisOfs)) * cbetinv +
    xcenAdj;
  yy -= yProjOffset;
  //lines[pitch * i + j] = yy; return;

  // Precalculate some items, doesn't help
  zzlim = ithick + 1 - ycenAdj;
  lbaseMtol = lsliceBase - ytol;
  llastPtol = lsliceLast + ytol;
  dxWarpInv = 1. / dxWarpDelz;

  // Move on ray up in Z
  lastZdone = 0;
              
  while (zz < zzlim && !lastZdone) {

    // xprojMin/Max already adjusted by 5
    if (xproj < xprojMin || xproj > xprojMax) {
      sum = sum + pmean;
      //if (zz + ycenAdj > ithick - skip) {lines[pitch * i + j] = 0; return;}
    } else {
      loadedProjectingPoint(xproj, yproj, zz, ofsxpz, ofsypf, ofsypz,
                            nxload, lsliceBase, lsliceLast, &xx, &yy);
      //if (zz + ycenAdj > ithick - skip) {lines[pitch * i + j] = yy; return;}

      // If X or Y is out of bounds, fill with mean
      if (yy < lbaseMtol || yy > llastPtol || xx < 1. || xx >= nxload) {
        sum = sum + pmean;
      } else {

        // otherwise, get x, y, z indexes, clamp y to limits, allow
        // a fractional Z pixel at top of volume
        xxtex = xx - 0.5f;
        yy = max((float)lsliceBase, min(lsliceLast - 0.01, yy));
        zslice = zz + ycenAdj;
        frac = 1.;
        if (zslice > ithick) {
          frac = 1. - (zslice - (int)zslice);
          zslice = ithick - 0.5f;
          lastZdone = 1;
        } else
          zslice -= 0.5f;
                     
        // Do the interpolation
        if (use3D > 0) {
          sum += frac * tex3D(projtex3D, xxtex, zslice, yy - baseAdj);
#ifdef HAS_LAYERS
        } else if (use3D < 0) {
          iy = yy;
          fy = yy - iy;
          lslice = iy - lsliceBase;
          sum += frac * ((1. - fy) * tex2DLayered(projtexLayer, xxtex, zslice, lslice) +
                         fy * tex2DLayered(projtexLayer, xxtex, zslice, lslice + 1));
#endif
        } else {
          iy = yy;
          fy = yy - iy;
          zslice += (iy - lsliceBase) * ithick;
          sum += frac * ((1. - fy) * tex2D(projtex2D, xxtex, zslice) +
                         fy * tex2D(projtex2D, xxtex, zslice + ithick));
        }

        // ELIMINATED JUMPING, IT TAKES 50% LONGER
      }
    }
                 
    // Adjust Z by local factor, move X approximately for next pixel
    zind = max(0., min(nWarpDelz - 1., xx * dxWarpInv));
    zz = zz + tex2D(delztex, zind, fline);
    xx = xx + sinBeta;
  }
  lines[pitch * i + j] = sum;
}

// Function to do reprojection with local alignments
int gpureprojlocal
(float *lines, float *sinBeta, float *cosBeta, float *sinAlpha, float *cosAlpha,
 float *xzfac, float *yzfac, int *nxWarp, int *nyWarp, int *ixStartWarp, 
 int *iyStartWarp, int *iDelXwarp, int *iDelYwarp, float *warpDelz, int *nWarpDelz, 
 float *dxWarpDelz, float *xprojMin, float *xprojMax, int *lsliceStart, int *lsliceEnd,
 int *ithick, int *iview, float *xcenOut, float *xcenIn, float *axisXoffset, 
 int *minXload, float *xProjOffset, float *ycenAdj, float *yProjOffset,
 float *centerSlice, float *pmean)
{
  int blockX = 16;
  int numLines = *lsliceEnd + 1 - *lsliceStart;
  int lastSlice = sLsliceFirst + sNumLoadedPlanes - 1;
  int nbd, nbp;
  float xcenAdj, salfsbet, calsbetpxz, salfmyz, cbetinv,calfinv;

  xcenAdj = *xcenOut - (*minXload-1);
  salfsbet = *sinAlpha * *sinBeta;
  calsbetpxz = *cosAlpha * *sinBeta + *xzfac;
  salfmyz = *sinAlpha - *yzfac;
  cbetinv = 1. / *cosBeta;
  calfinv = 1. / *cosAlpha;
  nbd = (int)floor(*yProjOffset + 0.5);

  // Compute the local projection factors
  dim3 blockFac(blockX, 16, 1);
  dim3 gridFac((sNxPlane + blockFac.x - 1) / blockFac.x, 
                (sNumLoadedPlanes + blockFac.y - 1) / blockFac.y, 1);
  localProjFactors<<<gridFac, blockFac>>>
    (devXprojFix, devXprojZ, devYprojFix, devYprojZ, sLocalPitch / 4, *iview - 1, 
     sNumViews, sNxPlane, *minXload - 1, sLsliceFirst + nbd, sNumLoadedPlanes, *nxWarp,
     *nyWarp, *ixStartWarp, *iyStartWarp, *iDelXwarp, *iDelYwarp, *xcenOut, *xcenIn, 
     *xcenIn+*axisXoffset, *centerSlice);
  if (testReportErr("computing localProjFactors"))
      return 1;
  /* return (synchronizeCopySlice(devYprojFix, sLocalPitch, lines, sSliceWidth,
     numLines)); */

  if (cudaThreadSynchronize() != cudaSuccess) {
    pflerr("Error return from synchronizing after computing local factors");
    return 1;
  }

  // Load the texture arrays
  nbd = sizeof(float) * *nWarpDelz * numLines;
  nbp = sizeof(float) * sNxPlane;
  if (cudaMemcpyToArray(devDelz, 0, 0, warpDelz, nbd, cudaMemcpyHostToDevice)
      != cudaSuccess ||
      cudaMemcpy2DToArray(devLocalPfac, 0, 0, devXprojFix, sLocalPitch, nbp, 
                          sNumLoadedPlanes, cudaMemcpyDeviceToDevice) 
      != cudaSuccess ||
      cudaMemcpy2DToArray(devLocalPfac, 0, sNumLoadedPlanes, devXprojZ, sLocalPitch,
                          nbp, sNumLoadedPlanes, cudaMemcpyDeviceToDevice) 
      != cudaSuccess ||
      cudaMemcpy2DToArray(devLocalPfac, 0, 2*sNumLoadedPlanes, devYprojFix, sLocalPitch,
                          nbp, sNumLoadedPlanes, cudaMemcpyDeviceToDevice) 
      != cudaSuccess ||
      cudaMemcpy2DToArray(devLocalPfac, 0, 3*sNumLoadedPlanes, devYprojZ, sLocalPitch,
                          nbp, sNumLoadedPlanes, cudaMemcpyDeviceToDevice) 
      != cudaSuccess) {
    pflerr("Failed to copy local proj factors to texture array");
    return 1;
  }

  // Do the reprojection
  dim3 blockSize(blockX, 16, 1);
  dim3 gridSize((sSliceWidth + blockSize.x - 1) / blockSize.x, 
                (numLines + blockSize.y - 1) / blockSize.y, 1);
  reprojLocal<<<gridSize, blockSize>>>
    (devSlice, sSlicePitch / 4, *nWarpDelz, *dxWarpDelz, sNxPlane, sSliceWidth,
     *ithick, *lsliceStart, *lsliceEnd, sLsliceFirst, lastSlice, *xprojMin, *xprojMax,
     xcenAdj, *xcenIn + *axisXoffset, *xProjOffset, *centerSlice, *yProjOffset, *ycenAdj,
     *cosBeta, *sinBeta, cbetinv, calfinv, salfmyz, salfsbet, calsbetpxz, *pmean, 
     sUse3dTexture);
  if (testReportErr("for local reprojection"))
      return 1;
  return (synchronizeCopySlice(devSlice, sSlicePitch, lines, sSliceWidth,
                               numLines));
}

/*
 * UTILITY ROUTINES
 */
   
// Load cosine and sine beta into constant array and compute inverse cosine
static int loadBetaInvertCos(float *cosBeta, float *sinBeta, float *cosinv, 
                             int num)
{
  int i, iv;
  float yy;

  // Invert cosines with limit
  for (i = 0; i < num; i++) {
    yy = cosBeta[i];
    if (fabs(yy) < 0.001f)
      yy = yy >= 0 ? 0.001f : -0.001f;
    cosinv[i] = 1.f / yy;
  }

  // Copy sines/cosines
  iv = num * sizeof(float);
  if (cudaMemcpyToSymbol(tables, cosBeta, iv, 0, cudaMemcpyHostToDevice) ||
      cudaMemcpyToSymbol(tables, sinBeta, iv, SINOFS*4,
                            cudaMemcpyHostToDevice)) {
    pflerr("Failed to copy constant data to GPU");
    return 1;
  }
  return 0;
}

// Synchronize the threads and copy computed data back to caller's array
static int synchronizeCopySlice(float *devslc, int pitch, float *slice,
                                int width, int numLines)
{
  int sizetmp;
  if (cudaThreadSynchronize() != cudaSuccess) {
    pflerr("Error return from synchronizing after backprojection");
    return 1;
  }

  // Get slice back
  sizetmp = sizeof(float) * width;
  if (cudaMemcpy2D(slice, sizetmp, devslc, pitch, sizetmp, numLines, 
                   cudaMemcpyDeviceToHost) != cudaSuccess) {
    pflerr("Error copying slice back to host");
    return 1;
  }
  return 0;
}

// Test for and report error after executing threads           
static int testReportErr(const char *mess)
{
  cudaError_t err;
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    pflush("Error executing threads %s: %s\n", mess,
           cudaGetErrorString(err));
    return 1;
  }
  return 0;
}

// Print a message with flushes to get it out before fortran output
static void pflush(const char *format, ...)
{
  char errorMess[512];
  va_list args;
  va_start(args, format);
  vsprintf(errorMess, format, args);
  printf("%s", errorMess);
  fflush(stdout);  
  fflush(stdout);
  va_end(args);
}

// In case of error, find the error string and print it with message
static void pflerr(const char *format, ...)
{
  cudaError_t err;
  char errorMess[512];
  va_list args;
  va_start(args, format);
  vsprintf(errorMess, format, args);
  printf("%s", errorMess);
  err = cudaGetLastError();
  pflush(": %s\n", cudaGetErrorString(err));
  fflush(stdout);  
  fflush(stdout);
  va_end(args);
}

// Print appropriate error from allocation and free all arrays
static void allocerr(const char *mess, int *nplanes, int *firstNpl,
                     int *lastNpl, int ifcuda)
{
  const char *whichText[3] = {"first", "last", "only"};
  int which = 2;
  gpudone();
  if (*firstNpl != *lastNpl) {
    if (*nplanes == *firstNpl)
      which = 0;
    else if (*nplanes == *lastNpl)
      which = 1;
    else
      return;
  }
  if (ifcuda)
    pflerr("On %s try (for %d planes), %s", whichText[which], *nplanes, mess);
  else
    pflush("On %s try (for %d planes), %s", whichText[which], *nplanes, mess);
}


