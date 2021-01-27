/*
 *  gpuctf.cu  -  GPU module for CTF correction by ctfphaseflip
 *
 *  Author: David Mastronarde
 *
 *  Copyright (C) 2018 by  the Regents of the University of
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 *  $Id: ctfphaseflip.cpp,v 4d4a2804e2b6 2018/09/25 02:59:33 mast $
 */
#include <stdarg.h>
#include "cuda_runtime_api.h"
#include "cuda.h"
#include "cufft.h"
#include "b3dutil.h"
#include "gpuctf.h"
#include "frameutil.h"
#include "cppdefs.h"

// Static functions
//#define USE_DUMPS
#ifdef USE_DUMPS
static void dumpFFT(float *fft, int nxPad, int nyPad, const char *descrip, int doReal);
static void dumpImage(float *image, int nxDim, int nxPad, int nyPad, const char *descrip);
#endif
static int testReportErr(const char *mess);
static void pflerr(const char *format, ...);
static void freeCudaArray(float **array);
static void destroyPlan(cufftHandle &plan);

/*
 * Static variables in CPU or device
 */
static cufftHandle sForwardPlan = 0;
static cufftHandle sInversePlan = 0;
static float *sFullSlice = NULL;
static float *sFullOnDev = NULL;
static float *sFullPadImage = NULL;
static float *sFullCopy1 = NULL;
static float *sFullCopy2 = NULL;
static float *sStripArray1 = NULL;
static float *sStripArray2 = NULL;
static float *sOutputArray = NULL;
static bool sDoFullImages;
static int sNxSlice;
static int sNySlice;
static bool sNeedFullXform = false;
static int sStripXdim = 0;
static int sNxPad;
static int sNyPad;
static int sDebug = 0;
static double sWallCopy = 0.;
static double sWallInterp = 0.;
static double sWallPrep = 0.;
static double sWallFFT = 0.;
static double sWallCorrect = 0.;
static double sWallStart, sWallNow;

// Macros to get time only when needed
#define START_TIMER  if (sDebug) sWallStart = wallTime();
#define ADD_TIME(a) if (sDebug)                                         \
  {sWallNow = wallTime(); a += sWallNow - sWallStart; sWallStart = sWallNow;}

#if CUDA_VERSION < 4000
#define cudaDeviceSynchronize cudaThreadSynchronize
#endif

/*
 * KERNELS
 */

/*
 * Kernel to taper and pad an image
 */
__global__ void taperInPadKernel(float *fullImage, int nxDimIn, int ixStart, int nxBox,
                                 int iyStart, int nyBox, float *outArr, int nxDimOut,
                                 int nx, int ny, int nxTaper, int nyTaper, float dmean)
{
  int padInX, padInY, outInd, imageX, imageY, fullInd, xInTaper, yInTaper;
  float value, atten;
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix >= nx || iy >= ny)
    return;
  padInX = nx / 2 - nxBox / 2;
  padInY = ny / 2 - nyBox / 2;
  outInd = ix + iy * nxDimOut;
  imageX = ix - padInX;
  imageY = iy - padInY;
  
  if (imageX < 0 || imageX >= nxBox || imageY < 0 || imageY >= nyBox) {
    outArr[outInd] = dmean;
  } else {
    fullInd = imageX + ixStart + nxDimIn * (imageY + iyStart);
    value = fullImage[fullInd];
    xInTaper = min(imageX, (nxBox - 1) - imageX);
    yInTaper = min(imageY, (nyBox - 1) - imageY);
    if (xInTaper < nxTaper || yInTaper < nyTaper) {
      atten = min((xInTaper + 1.f) / (nxTaper + 1.f), (yInTaper + 1.f) / (nyTaper + 1.f));
      value = atten * (value - dmean) + dmean;
    }
    outArr[outInd] = value;
  }
}


/*
 * Kernel to correct the CTF
 */
__global__ void correctCtfKernel
(float *curStrip, int stripXdim, int ny, float freq_scalex, float freq_scaley,
 float pointDefocus, float cosAstig, float sinAstig, float focusSum, float focusDiff, 
 float cutonAngstroms, float phaseFracFactor, float phaseShift, float ampAngle, float C1,
 float C2, float scaleByPower, int powerIsHalf, int generalPower, float firstZeroFreqSq,
 float attenStartFrac, float minAttenFreqSq)
{
  float gx, gy, f2, denom, cosSum, phaseFrac = 1., waveAberration, attenFrac, ctf;
  int fx = blockIdx.x * blockDim.x + threadIdx.x;
  int fy = blockIdx.y * blockDim.y + threadIdx.y;
  int index;
  if (fx >= stripXdim / 2 || fy >= ny)
    return;
  index = fy * stripXdim + 2 * fx;
  if (fy > ny / 2)
    fy -= ny;
  gy = fy * freq_scaley;
  gx = fx * freq_scalex;
  f2 = gx * gx + gy * gy;
  if (focusDiff && (fx || fy)) {
    denom = sqrtf(f2);
    cosSum = (cosAstig * gx + sinAstig * gy) / denom;
    pointDefocus = focusSum + focusDiff * (2. * cosSum * cosSum - 1.);
  }
  if (cutonAngstroms > 0.)
    phaseFrac = phaseFracFactor * (1. - exp(-sqrtf(f2) / cutonAngstroms));
  waveAberration = (C2 * f2 - C1 * pointDefocus) * f2 - phaseFrac * phaseShift;

  // Produce a positive ctf for consistency and so it can be used for scaling
  // (Here is the formal equation before simplifying)
  /*ctf = -(sqrt(1 - ampContrast * ampContrast)) * sin(waveAberration)
    + ampContrast * cos(waveAberration);*/
  ctf = -sin(waveAberration - ampAngle);
  if (scaleByPower > 0. && f2 > minAttenFreqSq) {
    if (powerIsHalf)
      ctf = sqrt(fabs(ctf)) * (ctf >= 0. ? 1. : -1.);
    else if (generalPower) 
      ctf = pow(fabs(ctf), scaleByPower) * (ctf >= 0. ? 1. : -1.);
    if (f2 < firstZeroFreqSq) {
      attenFrac = (firstZeroFreqSq - f2) / 
        ((1. - attenStartFrac) * firstZeroFreqSq);
      ctf = attenFrac + (1. - attenFrac) * ctf;
    }

    // *= did not work here
    curStrip[index] = curStrip[index] * ctf;
    curStrip[index + 1] = curStrip[index + 1] * ctf;
  } else if (ctf < 0) {
    curStrip[index] = -curStrip[index];
    curStrip[index + 1] = -curStrip[index + 1];
  }
}

/*
 * Kernel to copy columns
 */
__global__ void copyColumnsKernel(float *curStrip, float *restoredArray, int nxDim,
                                  int nyFile, int stripXdim, int xoff, int yoff, 
                                  int startCol, int endCol)
{
  int column = blockIdx.x * blockDim.x + threadIdx.x + startCol;
  int row = blockIdx.y * blockDim.y + threadIdx.y;
  if (column <= endCol && row < nyFile)
    restoredArray[row * nxDim + column] =
      curStrip[(row + yoff) * stripXdim + column + xoff];
}

/*
 * Kernel to interpolate between two strips
 */
__global__ void interpolateKernel(float *curStrip, float *lastStrip, float *restoredArray,
                                  int nxDim, int nyFile, int stripXdim, int yoff, 
                                  int stripStride, int stripMid, int halfStrip,
                                  int curOffset, int lastOffset)
{
  int stripDist0, stripDist1;
  float curFrac, lastFrac;
  int startCol = stripMid - stripStride + 1;
  int column = blockIdx.x * blockDim.x + threadIdx.x + startCol;
  int row = blockIdx.y * blockDim.y + threadIdx.y;
  if (column <= stripMid && row < nyFile) {
    stripDist0 = column - stripMid + stripStride - 1;
    stripDist1 = stripMid + 1 - column;
    curFrac = stripDist0 / (float)stripStride;
    lastFrac = stripDist1 / (float)stripStride;
    restoredArray[row * nxDim + column] =
      curFrac * curStrip[(row + yoff) * stripXdim + curOffset + column] +
      lastFrac * lastStrip[(row + yoff) * stripXdim + lastOffset + column];
  }
}

/*
 * Kernel to copy along diagonals from a full image
 */
__global__ void copyDiagonalKernel
(float *curStrip, float *restoredArray, int nx, int nyFile, int stripXdim, int xoff,
 int yoff, float sinViewAxis, float cosViewAxis, float lowLim, float highLim)
{
  float axisDist;
  float xPixCenter = nx / 2.f - 0.5f;
  float yPixCenter = nyFile / 2.f - 0.5f;
  int column = blockIdx.x * blockDim.x + threadIdx.x;
  int row = blockIdx.y * blockDim.y + threadIdx.y;
  if (column < nx && row < nyFile) {
    axisDist = -sinViewAxis * (column - xPixCenter) + cosViewAxis * (row - yPixCenter);
    if (axisDist >= lowLim && axisDist <= highLim) {
      restoredArray[row * nx + column] =
        curStrip[(row + yoff) * stripXdim + column + xoff];
    }
  }
}

/*
 * Kernel to interpolate along diagonals between two full images
 */
__global__ void interpDiagonalsKernel
(float *curStrip, float *lastStrip, float *restoredArray, int nx, int nyFile,
 int stripXdim, int xoff, int yoff, float stripStride, float sinViewAxis,
 float cosViewAxis, float lastAxisDist, float curAxisDist)
{
  float axisDist, curAxFrac;
  float xPixCenter = nx / 2.f - 0.5f;
  float yPixCenter = nyFile / 2.f - 0.5f;
  int column = blockIdx.x * blockDim.x + threadIdx.x;
  int row = blockIdx.y * blockDim.y + threadIdx.y;
  if (column < nx && row < nyFile) {
    axisDist = -sinViewAxis * (column - xPixCenter) + cosViewAxis * (row - yPixCenter);
    if (axisDist >= lastAxisDist + 0.5f && axisDist <= curAxisDist + 0.5f) {
      axisDist = max(lastAxisDist, min(curAxisDist, axisDist));
      curAxFrac = (min(axisDist, curAxisDist) - lastAxisDist) / stripStride;
      restoredArray[row * nx + column] = 
        curAxFrac * curStrip[(row + yoff) * stripXdim + xoff + column] +
        (1. - curAxFrac) * lastStrip[(row + yoff) * stripXdim + xoff + column];
    }
  }
}

/*
 * EXTERNALLY CALLED FUNCTTIONS
 */

/*
 * Test whether a GPU is available, either a GPU of the given number if nGPU is
 * > 0, or the one with the best processing rate if nGPU is 0, and return the
 * memory in bytes.  Return value is 1 for success, 0 for failure.
 */
int gpuAvailable(int nGPU, float *memory, int debug)
{
  int current_device = 0;
  int device_count = 0;
  int totalCores, max_gflops_device;
  float gflops;
  struct cudaDeviceProp device_properties, best_properties;

  // The Mac mini comes through with a clock rate of 0 so allow a 0 product
  float max_gflops = -1.;
  sDebug = debug;
  *memory = 0;
  cudaGetDeviceCount( &device_count );
  if (debug) {
#if CUDA_VERSION >= 3000
    int version, version2;
    cudaRuntimeGetVersion(&version2);
    cudaDriverGetVersion(&version);
    utilPrint("CUDA version - driver: %d.%02d  runtime: %d.%02d\n", version / 1000,
              version % 1000, version2 / 1000, version2 % 1000);
#endif
    utilPrint("Device count = %d\n", device_count);
  }
  if (nGPU != 0) {
    if (nGPU < 0 || nGPU > device_count) {
      utilPrint("The requested GPU number, %d, is out of range; there are only %d "
                "devices\n", nGPU, device_count);
      return 0;
    }
    current_device = nGPU - 1;
    device_count = nGPU;
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
    if (debug) {
      utilPrint("Device %d (%s): mp %d cores %d  cr %d  major %d minor %d  mem %.0f",
                current_device, device_properties.name,
                device_properties.multiProcessorCount, totalCores,
                device_properties.clockRate, device_properties.major,
                device_properties.minor, (float)device_properties.totalGlobalMem);
#if CUDA_VERSION >= 4000
      utilPrint("  tex1d %d", device_properties.maxTexture1DLinear);
#endif
      utilPrint("\n");
    }
    gflops = totalCores * (float)device_properties.clockRate;

    // Exclude emulation mode (?) which shows up on the Mac
    if( gflops > max_gflops && device_properties.major != 9999) {
      max_gflops = gflops;
      max_gflops_device = current_device;
      best_properties = device_properties;
    }
  }
    
  if (max_gflops_device >= 0) {
    *memory = best_properties.totalGlobalMem;
    if (cudaSetDevice(max_gflops_device) != cudaSuccess) {
      pflerr("Error selecting GPU device %d", max_gflops_device + 1);
      return 0;
    }
    return 1;
  }
  return 0;
}

/*
 * Initialize operations for one slice, doing some one-time operations in the first call
 * i.e., allocate arrays and set up FFT plans
 */
int gpuInitializeSlice(float *sliceData, int nxFile, int nyFile, int stripXdim, int nxPad, 
                       int nyPad, bool doFullImages)
{
  int error;
  int sliceBytes = nxFile * nyFile * sizeof(float);
  int padBytes = stripXdim * nyPad * sizeof(float);
  sFullSlice = sliceData;
  sDoFullImages = doFullImages;
  sNxSlice = nxFile;
  sNySlice = nyFile;
  
  if (!sFullOnDev) {
    // First slice: allocate array
    if (cudaMalloc((void **)&sFullOnDev, sliceBytes) != cudaSuccess || 
        cudaMalloc((void **)&sOutputArray, sliceBytes) != cudaSuccess) {
      pflerr("Failed to allocate full slice array on GPU");
      return 1;
    }
    
    // If doing full slice, also allocate the padded array and its two copies here
    if (doFullImages) {
      if (cudaMalloc((void **)&sFullPadImage, padBytes) != cudaSuccess || 
          cudaMalloc((void **)&sFullCopy1, padBytes) != cudaSuccess || 
          cudaMalloc((void **)&sFullCopy2, padBytes) != cudaSuccess) {
        pflerr("Failed to allocate padded full slice arrays on GPU");
        return 1;
      }
    }
  }

  if (!doFullImages && stripXdim != sStripXdim) {
    freeCudaArray(&sStripArray1);
    freeCudaArray(&sStripArray2);
    if (cudaMalloc((void **)&sStripArray1, padBytes) != cudaSuccess || 
        cudaMalloc((void **)&sStripArray2, padBytes) != cudaSuccess) {
      pflerr("Failed to allocate strip arrays on GPU");
      return 1;
    }
  }

  // Make FFT plans once or whenever strip size changes
  if (stripXdim != sStripXdim) {
    destroyPlan(sForwardPlan);
    destroyPlan(sInversePlan);
    error = cufftPlan2d(&sForwardPlan, nyPad, nxPad, CUFFT_R2C);

    if (error == CUFFT_SUCCESS)
      error = cufftPlan2d(&sInversePlan, nyPad, nxPad, CUFFT_C2R);
    if (error != CUFFT_SUCCESS) {
      utilPrint("Failed to make plan for FFTs (error %d)\n", error);
      return 1;
    }
  }    

  // Copy full array to the GPU
  START_TIMER;
  if (cudaMemcpy(sFullOnDev, sFullSlice, sNxSlice * sNySlice * sizeof(float), 
                 cudaMemcpyHostToDevice) != cudaSuccess) {
    pflerr("Failed to copy slice image to GPU array");
    return 1;
  }
  ADD_TIME(sWallCopy);

  sNeedFullXform = sDoFullImages;
  sStripXdim = stripXdim;
  sNxPad = nxPad;
  sNyPad = nyPad;
  return 0;
}

/*
 * Extract a strip between stripBegin and stripEnd with the given tapering and padded
 * as specified previously, and get the transform.  If doing full images, do this operation
 * only on the first strip and after that, copy the transform to the appropriate buffer
 */
int gpuExtractAndTransform(int stripInd, int stripBegin, int stripEnd, int nxTaper,
                           int nyTaper)
{
  float *curStrip;
  float dmean;
  int err;
  int blockX = 32;
  if (sDoFullImages)
    curStrip = sFullPadImage;
  else
    curStrip = (stripInd % 2) ? sStripArray2 : sStripArray1;
    
  // Need to taper-pad into curStrip for the strips or the full padded image
  if (!sDoFullImages || sNeedFullXform) {
    if (sDoFullImages) {
      stripBegin = 0;
      stripEnd = sNxSlice - 1;
    }
    dmean = sliceEdgeMean(sFullSlice, sNxSlice, stripBegin, stripEnd, 0, sNySlice - 1);
    START_TIMER;
    dim3 blockSize(blockX, 8, 1);
    dim3 gridSize((sNxPad + blockSize.x - 1) / blockSize.x,
                  (sNyPad + blockSize.y - 1) / blockSize.y, 1);
    taperInPadKernel<<<gridSize, blockSize>>>
      (sFullOnDev, sNxSlice, stripBegin, stripEnd + 1 - stripBegin, 0, sNySlice, curStrip,
       sStripXdim, sNxPad, sNyPad, nxTaper, nyTaper, dmean);
    if (testReportErr("to taper-pad image on GPU"))
      return 1;
    ADD_TIME(sWallPrep);
    if (stripInd < 5) {
      //dumpImage(curStrip, sStripXdim, sNxPad, sNyPad, "taperpad");
    }
    // And then take the FFT
    err = cufftExecR2C(sForwardPlan, curStrip, (cufftComplex *)curStrip);
    cudaDeviceSynchronize();
    if (err != CUFFT_SUCCESS) {
      utilPrint("Failure in forward FFT on GPU (CUFFT error %d)\n", err);
      return 1;
    }
    ADD_TIME(sWallFFT);
    sNeedFullXform = false;
    if (stripInd < 5) {
      //dumpFFT(curStrip, sNxPad, sNyPad, "taperpad-fft", 0);
    }
  }
  
  // If doing full image, now copy to current copy array
  if (sDoFullImages) {
    START_TIMER;
    if (cudaMemcpy((stripInd % 2) ? sFullCopy2 : sFullCopy1, sFullPadImage,
                   sStripXdim * sNyPad * sizeof(float), cudaMemcpyDeviceToDevice) != 
        cudaSuccess) {
      pflerr("Failed to copy full padded FFT to GPU current copy array");
      return 1;
    }
    ADD_TIME(sWallCopy);
  }
  
  return 0;
}

/*
 * Correct the CTF in the current strip and back-transform it
 */
int gpuCorrectCTF(int stripInd, float freq_scalex, float freq_scaley, float pointDefocus,
                  float cosAstig, float sinAstig, float focusSum, float focusDiff,
                  float cutonAngstroms, float phaseFracFactor, float phaseShift, 
                  float ampAngle, float C1, float C2, float scaleByPower, bool powerIsHalf,
                  bool generalPower, float firstZeroFreqSq, float attenStartFrac,
                  float minAttenFreqSq)
{
  float *curStrip;
  int err;
  int blockX = 32;
  if (sDoFullImages)
    curStrip = (stripInd % 2) ? sFullCopy2 : sFullCopy1;
  else
    curStrip = (stripInd % 2) ? sStripArray2 : sStripArray1;

  // Correct the CTF in the current strip
  START_TIMER;
  dim3 blockSize(blockX, 8, 1);
  dim3 gridSize((sStripXdim / 2 + blockSize.x - 1) / blockSize.x,
                (sNyPad + blockSize.y - 1) / blockSize.y, 1);
  correctCtfKernel<<<gridSize, blockSize>>>
    (curStrip, sStripXdim, sNyPad, freq_scalex, freq_scaley, pointDefocus, cosAstig,
     sinAstig, focusSum, focusDiff, cutonAngstroms, phaseFracFactor, phaseShift, 
     ampAngle, C1, C2, scaleByPower, powerIsHalf ? 1 : 0, generalPower ? 1 : 0,
     firstZeroFreqSq, attenStartFrac, minAttenFreqSq);
  if (testReportErr("to apply CTF corrections on GPU"))
    return 1;
  ADD_TIME(sWallCorrect);
  if (stripInd < 5) {
    //dumpFFT(curStrip, sNxPad, sNyPad, "corr-fft", 0);
  }
    
  // And then take the inverse FFT
  err = cufftExecC2R(sInversePlan, (cufftComplex *)curStrip, curStrip);
  cudaDeviceSynchronize();
  if (err != CUFFT_SUCCESS) {
    utilPrint("Failure in inverse FFT on GPU (CUFFT error %d)\n", err);
    return 1;
  }
  ADD_TIME(sWallFFT);
  //if (stripInd < 5) {
  //dumpImage(curStrip, sStripXdim, sNxPad, sNyPad, "corr-img");
    //}

  return 0;
}

/*
 * Interpolate columns between the last and current strip from stripMid - stripStride + 1 
 * to stripMid, where the offsets are ADDED when getting the strip coordinates 
 */
int gpuInterpolateColumns(int stripInd, int yoff, int stripStride, int stripMid,
                          int halfStrip, int curOffset, int lastOffset)
{
  float *curStrip, *lastStrip;
  int blockX = 32;
  if (sDoFullImages) {
    curStrip = (stripInd % 2) ? sFullCopy2 : sFullCopy1;
    lastStrip = (stripInd % 2) ? sFullCopy1 : sFullCopy2;
  } else {
    curStrip = (stripInd % 2) ? sStripArray2 : sStripArray1;
    lastStrip = (stripInd % 2) ? sStripArray1 : sStripArray2;
  }
  START_TIMER;
  dim3 blockSize(blockX, 8, 1);
  dim3 gridSize((stripStride + blockSize.x - 1) / blockSize.x,
                (sNySlice + blockSize.y - 1) / blockSize.y, 1);
  interpolateKernel<<<gridSize, blockSize>>>
    (curStrip, lastStrip, sOutputArray, sNxSlice, sNySlice, sStripXdim, yoff,
     stripStride, stripMid, halfStrip, curOffset, lastOffset);
  ADD_TIME(sWallInterp);
  return testReportErr("to interpolate columns on GPU");
}

/*
 * Copy columns from the current strip between startCol and endCol to the output, where
 * the offsets are ADDED to output coordinates to get the strip coordinates 
 */
int gpuCopyColumns(int stripInd, int xoff, int yoff, int startCol, int endCol)
{
  float *curStrip;
  int blockX = 32;
  if (sDoFullImages)
    curStrip = (stripInd % 2) ? sFullCopy2 : sFullCopy1;
  else
    curStrip = (stripInd % 2) ? sStripArray2 : sStripArray1;

  START_TIMER;
  dim3 blockSize(blockX, 8, 1);
  dim3 gridSize(((endCol + 1 - startCol) + blockSize.x - 1) / blockSize.x,
                (sNySlice + blockSize.y - 1) / blockSize.y, 1);
  copyColumnsKernel<<<gridSize, blockSize>>>
    (curStrip, sOutputArray, sNxSlice, sNySlice, sStripXdim, xoff, yoff, startCol,
     endCol);
  ADD_TIME(sWallInterp);
  return testReportErr("to copy columns on GPU");
}

/*
 * Interpolate along diagonals
 */
int gpuInterpDiagonals(int stripInd, int xoff, int yoff, int stripStride, 
                       float sinViewAxis, float cosViewAxis, float lastAxisDist,
                       float curAxisDist)
{
  float *curStrip, *lastStrip;
  int blockX = 32;
  curStrip = (stripInd % 2) ? sFullCopy2 : sFullCopy1;
  lastStrip = (stripInd % 2) ? sFullCopy1 : sFullCopy2;
  START_TIMER;
  dim3 blockSize(blockX, 8, 1);
  dim3 gridSize((sNxSlice + blockSize.x - 1) / blockSize.x,
                (sNySlice + blockSize.y - 1) / blockSize.y, 1);
  interpDiagonalsKernel<<<gridSize, blockSize>>>
    (curStrip, lastStrip, sOutputArray, sNxSlice, sNySlice, sStripXdim, xoff, yoff, 
     (float)stripStride, sinViewAxis, cosViewAxis, lastAxisDist, curAxisDist);
  ADD_TIME(sWallInterp);
  return testReportErr("to interpolate diagonals on GPU");
}

/*
 * Copy along diagonals
 */
int gpuCopyDiagonals(int stripInd, int xoff, int yoff, float sinViewAxis, 
                     float cosViewAxis, float lowLim, float highLim)
{
  float *curStrip;
  int blockX = 32;
  curStrip = (stripInd % 2) ? sFullCopy2 : sFullCopy1;
  START_TIMER;
  dim3 blockSize(blockX, 8, 1);
  dim3 gridSize((sNxSlice + blockSize.x - 1) / blockSize.x,
                (sNySlice + blockSize.y - 1) / blockSize.y, 1);
  copyDiagonalKernel<<<gridSize, blockSize>>>
    (curStrip, sOutputArray, sNxSlice, sNySlice, sStripXdim, xoff, yoff, sinViewAxis,
     cosViewAxis, lowLim, highLim);
  ADD_TIME(sWallInterp);
  return testReportErr("to copy diagonals on GPU");
}

/*
 * Return the corrected image
 */
int gpuReturnImage(float *finalImage)
{
  START_TIMER;
  if (cudaMemcpy(finalImage, sOutputArray, sNxSlice * sNySlice * sizeof(float),
                 cudaMemcpyDeviceToHost) != cudaSuccess) {
    pflerr("Failed to copy restored image back from GPU");
    return 1;
  }
  ADD_TIME(sWallCopy);
  for (int i = 0; i < sNxSlice * sNySlice; i++)
    finalImage[i] /= sNxPad * sNyPad;
  return 0;
}

/*
 * Return the diagnostic times
 */
void gpuGetTimes(double &copy, double &prep, double &FFT, double &correct, double &interp)
{
  copy = sWallCopy;
  prep = sWallPrep;
  FFT = sWallFFT;
  correct = sWallCorrect;
  interp = sWallInterp;
}

/*
 * UTILITIES
 */

/*
 * Cleanup functions when changing sizes
 */
static void freeCudaArray(float **array)
{
  if (*array)
    cudaFree(*array);
  *array = NULL;
}

static void destroyPlan(cufftHandle &plan)
{
  if (plan)
    cufftDestroy(plan);
  plan = 0;
}

/*
 * Test for an error launching threads then synchronize, which flushes out other errors
 */
static int testReportErr(const char *mess)
{
  cudaError_t err;
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    utilPrint("Error executing threads %s: %s\n", mess,
              cudaGetErrorString(err));
    return 1;
  }
  if (cudaDeviceSynchronize() != cudaSuccess) {
    pflerr("Error return from synchronizing after running kernel %s", mess);
    return 1;
  }
  return 0;
}

/*
 * In case of error, find the error string and print it with message
 */
static void pflerr(const char *format, ...)
{
  cudaError_t err;
  char errorMess[512];
  va_list args;
  va_start(args, format);
  vsprintf(errorMess, format, args);
  err = cudaGetLastError();
  utilPrint("%s: %s\n", errorMess, cudaGetErrorString(err));
  va_end(args);
}

/*
 * Wrappers to the utility dump functions
 */
#ifdef USE_DUMPS
static void dumpFFT(float *fft, int nxPad, int nyPad, const char *descrip, int doReal)
{
  int sizeTmp = (nxPad + 2) * nyPad;
  float *temp = B3DMALLOC(float, sizeTmp);
  if (!temp)
    return;
  if (cudaMemcpy(temp, fft, sizeTmp * sizeof(float), cudaMemcpyDeviceToHost) ==
      cudaSuccess)
    utilDumpFFT(temp, nxPad, nyPad, descrip, doReal, 0, 1);
  free(temp);
}

static void dumpImage(float *image, int nxDim, int nxPad, int nyPad, const char *descrip)
{
  int sizeTmp = nxDim * nyPad;
  float *temp = B3DMALLOC(float, sizeTmp);
  if (!temp)
    return;
  if (cudaMemcpy(temp, image, sizeTmp * sizeof(float), cudaMemcpyDeviceToHost) == 
      cudaSuccess)
    utilDumpImage(temp, nxDim, nxPad, nyPad, 0, descrip);
  free(temp);
}
#endif

