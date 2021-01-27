
/*
 *  gpuframe.cu -- Kernel and supporting code for frame summing on GPU
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 2016 by the Regents of the University of 
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
#include "gpuframe.h"
#include "framealign.h"
#include "frameutil.h"

#define PI 3.141593
#define FULL_INDENT 8

// Surfaces save ~1% in pre-processing and 7-10% in reduction, as well as saving memory
#if CUDA_VERSION < 4000
#define cudaDeviceSynchronize cudaThreadSynchronize
#define NO_SURFACES
#endif

static void destroyPlan(cufftHandle &plan);
static void free2Darray(cudaArray **array);
static cudaError_t malloc2Darray(cudaArray **arrayPtr, int type, int sizeX, int sizeY);
static cudaError_t bindUnpadArray(cudaArray *arrayPtr, int type);

static FrameGPU sFGPU;

// Use almost all of constant memory for the sines and cosines
#define MAX_KERNEL 100
#define INT_KERNEL_SCALE 16000.f
#define MAX_TABLE (65536 / 4 - MAX_KERNEL)
__constant__ float trigTable[MAX_TABLE];
__constant__ int intKernel[MAX_KERNEL];

/*
 * Static CUDA variables (to keep them out of the .h file)
 */
static cufftHandle sFullForwardPlan = 0;
static cufftHandle sSumInversePlan = 0;
static cufftHandle sAlignForwardPlan = 0;
static cufftHandle sAlignInversePlan = 0;
static texture<float, 1, cudaReadModeElementType> sSumTex;
static texture<float, 1, cudaReadModeElementType> sOddTex;
static texture<float, 1, cudaReadModeElementType> sFullTex;
static texture<float, 1, cudaReadModeElementType> sMaskTex;
static texture<float, 1, cudaReadModeElementType> sDWFilterTex;
static texture<float, 2, cudaReadModeElementType> sUnpadFloatTex;
static texture<unsigned char, 2, cudaReadModeElementType> sUnpadByteTex;
static texture<short int, 2, cudaReadModeElementType> sUnpadShortTex;
static texture<unsigned short int, 2, cudaReadModeElementType> sUnpadUShortTex;
static texture<float, 2, cudaReadModeElementType> sGainRefTex;
static texture<unsigned char, 2, cudaReadModeElementType> sDefectTex;
static cudaChannelFormatDesc sChanDesc;
static cudaChannelFormatDesc sByteChanDesc;
static cudaChannelFormatDesc sShortChanDesc;
static cudaChannelFormatDesc sUShortChanDesc;
static cudaArray *sGainRefArray = NULL;
static cudaArray *sDefectMapArray = NULL;
static cudaArray *sTempRawArray = NULL;
static cudaArray *sTempFloatArray = NULL;
#ifndef NO_SURFACES
static surface<void, 2> sTempSurfRef;
#endif
std::vector<cudaArray *>sSavedUnpadded;
static cudaArray *sReducedInXarray = NULL;

#define START_TIMER  if (mTrackTime) mWallStart = wallTime();
#define ADD_TIME(a) if (mTrackTime) a += wallTime() - mWallStart;

#define SETUP_TEXTURE(tex)                    \
  tex.addressMode[0] = cudaAddressModeClamp;    \
  tex.addressMode[1] = cudaAddressModeClamp;    \
  tex.filterMode = cudaFilterModePoint;         \
  tex.normalized = false;

/*
 * KERNELS
 */

// Shift an FFT array and add into a sum array, possibly with size reduction
__global__ void shiftAndAddToSum(float *sumArr, float *nonDWsum, int ixStart,
                                 int iyInStart, int iyOutStart, int numXdo, int numYdo,
                                 int nxInDim, int nxOutDim, int tableYoffset, 
                                 int numDoseFilt, float doseFiltDelta, float freqDelX, 
                                 float freqDelY)
{
  int indIn, indOut, trigXind, trigYind, ixFull, iyFull;
  float real, imag, xcos, xsin, ycos, ysin, phre, phim, sumReal, sumImag, xfreq, yfreq;
  float freq, newReal, newImag, atten = 1.f;
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;

  // textures only make a tiny difference here
  if (ix < numXdo && iy < numYdo) {
    ixFull = ix + ixStart;
    iyFull = iy + iyInStart;
    if (numDoseFilt > 0) {

      // Determine attenuation factor based on frequency, or return if it will be 0
      // Multiple if statements testing against a maximum frequency do not save any time
      // Ie., testing "if (freq > (numDoseFilt - 1) * doseFiltDelta"
      xfreq = ixFull * freqDelX;
      yfreq = iyFull * freqDelY;
      if (yfreq > 0.5f)
        yfreq = 1.f - yfreq;
      freq = sqrt(xfreq * xfreq + yfreq * yfreq);
      atten = tex1Dfetch(sDWFilterTex, freq / doseFiltDelta + 0.5f);
    }
    indIn = 2 * (ixFull + nxInDim * iyFull);
    indOut = 2 * (ix + ixStart + nxOutDim * (iy + iyOutStart));
    real = tex1Dfetch(sFullTex, indIn);
    imag = tex1Dfetch(sFullTex, indIn + 1);
    //real = fullArr[indIn];
    //imag = fullArr[indIn + 1];
    trigXind = 2 * ix;
    trigYind = tableYoffset + 2 * iy;
    xcos = trigTable[trigXind];
    xsin = trigTable[trigXind + 1];
    ycos = trigTable[trigYind];
    ysin = trigTable[trigYind + 1];
    phre = xcos * ycos - xsin * ysin;
    phim = xsin * ycos + xcos * ysin;
    newReal = phre * real - phim * imag;
    newImag = phim * real + phre * imag;
    sumReal = tex1Dfetch(sSumTex, indOut);
    sumImag = tex1Dfetch(sSumTex, indOut + 1);
    sumArr[indOut] = sumReal + atten * newReal;
    sumArr[indOut + 1] = sumImag + atten * newImag;
    if (nonDWsum) {
      nonDWsum[indOut] += newReal;
      nonDWsum[indOut + 1] += newImag;
    }
  }
}

// Shift an array and store it back shifted, and add it to a sum, with no reduction
__global__ void shiftInPlaceAddToSum(float *fullArr, float *sumArr, int ixStart,
                                     int iyStart, int numXdo, int numYdo, int nxDim,
                                     int tableYoffset)
{
  int indIn, trigXind, trigYind;
  float real, imag, xcos, xsin, ycos, ysin, phre, phim, sumReal, sumImag;
  float shiftReal, shiftImag;
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix < numXdo && iy < numYdo) {
    indIn = 2 * (ix + ixStart + nxDim * (iy + iyStart));
    real = tex1Dfetch(sFullTex, indIn);
    imag = tex1Dfetch(sFullTex, indIn + 1);
    //real = fullArr[indIn];
    //imag = fullArr[indIn + 1];
    trigXind = 2 * ix;
    trigYind = tableYoffset + 2 * iy;
    xcos = trigTable[trigXind];
    xsin = trigTable[trigXind + 1];
    ycos = trigTable[trigYind];
    ysin = trigTable[trigYind + 1];
    phre = xcos * ycos - xsin * ysin;
    phim = xsin * ycos + xcos * ysin;
    sumReal = tex1Dfetch(sSumTex, indIn);
    sumImag = tex1Dfetch(sSumTex, indIn + 1);
    shiftReal = phre * real - phim * imag;
    shiftImag = phim * real + phre * imag;
    sumArr[indIn] = sumReal + shiftReal;
    sumArr[indIn + 1] = sumImag + shiftImag;
    //sumArr[indIn] += phre * real - phim * imag;
    //sumArr[indIn + 1] += phim * real + phre * imag;
    fullArr[indIn] = shiftReal;
    fullArr[indIn + 1] = shiftImag;
  }
}

// Add the odd sum to the even sum
__global__ void addOddToEvenSum(float *evenArr, int sumXplus, int sumYpad)
{

  // textures did not speed this up
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix < sumXplus && iy < sumYpad) {
    int ind = ix + sumXplus * iy;
    evenArr[ind] = tex1Dfetch(sSumTex, ind) + tex1Dfetch(sOddTex, ind);
  }
}

// Subtract one array from another, multiple by a filter mask, and store in sumArr
__global__ void subtractFilterSum(float *subArr, float *sumArr, int sumXplus, int sumYpad)
{
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix < sumXplus && iy < sumYpad) {
    int ind = ix + sumXplus * iy;
    sumArr[ind] = (tex1Dfetch(sSumTex, ind) - subArr[ind]) * tex1Dfetch(sMaskTex, ind);
  }
}

// Apply filter mask to an FTT
__global__ void filterAlignFFT(float *fft, int alignXplus, int alignYpad)
{
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix < alignXplus && iy < alignYpad) {
    int ind = ix + alignXplus * iy;
    fft[ind] *= tex1Dfetch(sMaskTex, ind);
  }
}

// Compute conjugate product of two FFTs and store in prod array
__global__ void conjugateProduct(float *array, float *brray, float *prod, int nxFFT,
                                 int alignYpad)
{
  int jx, jp1;
  float a, b, c, d;
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix < nxFFT && iy < alignYpad) {
    jx = 2 * (ix + nxFFT * iy);
    jp1 = jx + 1;
    a = array[jx];
    b = array[jp1];
    c = brray[jx];
    d = brray[jp1];
    prod[jx] = a * c + b * d;
    prod[jp1] = b * c - a * d;
  }
}

// Extract the corners of a cross-correlation image into small array with origin at center
__global__ void wrapCorners(float *fromArr, float *toArr, int nxFrom, int nyFrom, 
                            int nxTo, int nyTo, int ixFrom, int iyFrom)
{
  int fromX, fromY;
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix < nxFrom && iy < nyFrom) {
    fromX = (ix + ixFrom) % nxFrom;
    fromY = (iy + iyFrom) % nyFrom;
    toArr[ix + iy * nxTo] = tex1Dfetch(sFullTex, fromX + fromY * nxFrom);
  }
}

// Sum different numbers of FFTs into groups
__global__ void sum2IntoGroup(float *arr1, float *arr2, float *groupArr, int alignXplus,
                              int alignYpad)
{
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix < alignXplus && iy < alignYpad) {
    int ind = ix + alignXplus * iy;
    groupArr[ind] = arr1[ind] + arr2[ind];
  }
}

__global__ void sum3IntoGroup(float *arr1, float *arr2, float * arr3, float *groupArr,
                              int alignXplus, int alignYpad)
{
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix < alignXplus && iy < alignYpad) {
    int ind = ix + alignXplus * iy;
    groupArr[ind] = arr1[ind] + arr2[ind] + arr3[ind];
  }
}

__global__ void sum4IntoGroup(float *arr1, float *arr2, float *arr3, float *arr4, 
                              float *groupArr, int alignXplus, int alignYpad)
{
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix < alignXplus && iy < alignYpad) {
    int ind = ix + alignXplus * iy;
    groupArr[ind] = arr1[ind] + arr2[ind] + arr3[ind] + arr4[ind];
  }
}

__global__ void sum5IntoGroup(float *arr1, float *arr2, float *arr3, float *arr4, 
                              float *arr5, float *groupArr, int alignXplus, int alignYpad)
{
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix < alignXplus && iy < alignYpad) {
    int ind = ix + alignXplus * iy;
    groupArr[ind] = arr1[ind] + arr2[ind] + arr3[ind] + arr4[ind] + arr5[ind];
  }
}

// Macro for preprocessing a specific type with its texture
// Multiply by gain if doing that
// If pixel is above truncation limit or marked as a defect, average surrounding 
// box of 9x9 pixels excluding central 9, picking only pixels that pass the same test
#define PREPROC_PIXEL(mtyp, tex)                                        \
  case mtyp:                                                            \
  pixVal = tex2D(tex, ix, iy);                                          \
  if (doGain)                                                           \
    pixVal *= tex2D(sGainRefTex, ix, iy);                               \
  if ((truncLimit > 0 && pixVal > truncLimit) ||                        \
      (doDefects > 0 && tex2D(sDefectTex, ix, iy) > 0)) {               \
    for (delY = -4; delY <= 4; delY++) {                                \
      for (delX = -4; delX <= 4; delX++) {                              \
        if (delX >= -1 && delX <= 1 && delY >= -1 && delY <= 1)         \
          continue;                                                     \
        inX = ix + delX;                                                \
        inY = iy + delY;                                                \
        pixVal = tex2D(tex, inX, inY);                                  \
        if (doGain)                                                     \
          pixVal *= tex2D(sGainRefTex, inX, inY);                       \
        if ((truncLimit > 0 && pixVal > truncLimit) ||                  \
            (doDefects > 0 && tex2D(sDefectTex, inX, inY) > 0))         \
          continue;                                                     \
        sum += pixVal;                                                  \
        nsum++;                                                         \
      }                                                                 \
    }                                                                   \
    pixVal = sum / nsum;                                                \
  }                                                                     \
  break;

// Kernel to pre-process and writ e to the surface reference
__global__ void preprocessFrame(float *linearArr, int type, int nx, int ny, int doGain, 
                                float truncLimit, int doDefects)
{
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  float pixVal, sum = 0.;
  int nsum = 0, inX, inY, delX, delY;
  if (ix >= nx || iy >= ny)
    return;
  switch (type) {
    PREPROC_PIXEL(MRC_MODE_BYTE, sUnpadByteTex);
    PREPROC_PIXEL(MRC_MODE_SHORT, sUnpadShortTex);
    PREPROC_PIXEL(MRC_MODE_USHORT, sUnpadUShortTex);
    PREPROC_PIXEL(MRC_MODE_FLOAT, sUnpadFloatTex);
  }
#ifdef NO_SURFACES
  linearArr[ix + iy * nx] = pixVal;
#else
  surf2Dwrite(pixVal, sTempSurfRef, ix * sizeof(float), iy);
#endif
}

// Pad a full image with noise that tapers down to the mean
__global__ void noiseTaperPad(int type, int nxBox, int nyBox, float *outArr,
                              int nxDimOut, int nx, int ny, int noiseLength, 
                              int noiseRows, int cornerSize, float dmean, int seed)
{
  int padInX, padInY, outInd, imageX, imageY, outInX = 0, xNear, outInY = 0, yNear;
  float xAtten, yAtten, atten;
  int nxNoise, nyNoise, ixNoise, iyNoise, ixHigh, iyHigh, inX, inY;
  int indInBox, pseudo, ixInBox, iyInBox;
  int xDirToData = 1, yDirToData = 1;
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix >= nx || iy >= ny)
    return;
  padInX = (nx - nxBox) / 2;
  padInY = (ny - nyBox) / 2;
  outInd = ix + iy * nxDimOut;
  imageX = ix - padInX;
  imageY = iy - padInY;
  xDirToData = 1;
  
  // Determine if it is outside the image, the direction toward data, and nearest pixel 
  // in each direction, as well as attenuation for the taper
  if (imageX < 0) {
    outInX = 1;
    xNear = 0;
    xAtten = (float)ix / (padInX + 1.f);
  }
  if (imageX >= nxBox) {
    outInX = 1;
    xDirToData = -1;
    xNear = nxBox - 1;
    xAtten = (float)(nx - ix - 1) / (padInX + 1.f);
  }
  if (imageY < 0)  {
    outInY = 1;
    yNear = 0;
    yAtten = (float)iy / (padInY + 1.);
  }
  if (imageY >= nyBox) {
    outInY = 1;
    yDirToData = -1;
    yNear = nyBox - 1;
    yAtten = (float)(ny - iy - 1) / (padInY + 1.f);
  }

  // If inside, we're done, copy the pixel
  if (!outInX && !outInY) {
    switch (type) {
    case MRC_MODE_BYTE:
      outArr[outInd] = tex2D(sUnpadByteTex, imageX, imageY);
      break;
    case MRC_MODE_SHORT:
      outArr[outInd] = tex2D(sUnpadShortTex, imageX, imageY);
      break;
    case MRC_MODE_USHORT:
      outArr[outInd] = tex2D(sUnpadUShortTex, imageX, imageY);
      break;
    default:
      outArr[outInd] = tex2D(sUnpadFloatTex, imageX, imageY);
      break;
    }
  } else {

    // Otherwise is in a corner
    if (outInX && outInY) {
      nxNoise = cornerSize;
      nyNoise = cornerSize;
      ixNoise = xNear;
      iyNoise = yNear;
      atten = max(xAtten, yAtten);

      // Or in left or right border, position box along edge
    } else if (outInX) {
      ixNoise = xNear;
      nxNoise = noiseRows;
      iyNoise = imageY - noiseLength / 2;
      
      // Expand box in other direction if it is near corner
      if (iyNoise < 0) {
        nxNoise = min(cornerSize, nxNoise - iyNoise);
        iyNoise = 0;
      }
      iyHigh = imageY + noiseLength / 2;  // One past
      if (iyHigh > nyBox) {
        nxNoise = min(cornerSize, nxNoise + (iyHigh - nyBox));
        iyHigh = nyBox;
      }
      nyNoise = iyHigh - iyNoise;
      atten = xAtten;

      // Or in top or bottom border, do same tests
    } else {
      iyNoise = yNear;
      nyNoise = noiseRows;
      ixNoise = imageX - noiseLength / 2;
      if (ixNoise < 0) {
        nyNoise = min(cornerSize, nyNoise - ixNoise);
        ixNoise = 0;
      }
      ixHigh = imageX + noiseLength / 2;  // One past
      if (ixHigh > nxBox) {
        nyNoise = min(cornerSize, nyNoise + (ixHigh - nxBox));
        ixHigh = nxBox;
      }
      nxNoise = ixHigh - ixNoise;
      atten = yAtten;
    } 

    // Get the random position, convert to an index in box and coordinates in box
    pseudo = (seed + 197 * (ix + 157) * (iy + 179)) & 0xFFFFF;
    pseudo = (197 * (pseudo + 1)) & 0xFFFFF;
    indInBox = (pseudo >> 4) % (nxNoise * nyNoise);
    ixInBox = indInBox % nxNoise;
    iyInBox = indInBox / nxNoise;
    inX = ixNoise + ixInBox * xDirToData;
    inY = iyNoise + iyInBox * yDirToData;
    switch (type) {
    case MRC_MODE_BYTE:
      outArr[outInd] = atten * (tex2D(sUnpadByteTex, inX, inY) - dmean) + dmean;
      break;
    case MRC_MODE_SHORT:
      outArr[outInd] = atten * (tex2D(sUnpadShortTex, inX, inY) - dmean) + dmean;
      break;
    case MRC_MODE_USHORT:
      outArr[outInd] = atten * (tex2D(sUnpadUShortTex, inX, inY) - dmean) + dmean;
      break;
    default:
      outArr[outInd] = atten * (tex2D(sUnpadFloatTex, inX, inY) - dmean) + dmean;
      break;
    }
  }
}


/*
 * Image reduction: column reduction macros and kernel, followed by row reduction kernel
 * This is an odd set of macros to minimize the number of lines in the column reduction 
 * kernel.  It saves about 25% to specify the unroll
 */
#define REDUCE_COL_BYTE                                                 \
  for (inX = 0; inX < numLoop; inX++)                                   \
    isum += intKernel[inX] * tex2D(sUnpadByteTex, inX + addXoffset, iySrc); \
  sum = isum / INT_KERNEL_SCALE;                                        \
  break;                                                                \
 case MRC_MODE_SHORT:

#define REDUCE_COL_SHORT                                                \
  for (inX = 0; inX < numLoop; inX++)                                   \
    isum += intKernel[inX] * tex2D(sUnpadShortTex, inX + addXoffset, iySrc); \
  sum = isum / INT_KERNEL_SCALE;                                        \
  break;                                                                \
 case MRC_MODE_USHORT:

#define REDUCE_COL_USHORT                                               \
  for (inX = 0; inX < numLoop; inX++)                                   \
    isum += intKernel[inX] * tex2D(sUnpadUShortTex, inX + addXoffset, iySrc); \
  sum = isum / INT_KERNEL_SCALE;                                        \
  break;                                                                \
 case MRC_MODE_FLOAT:

#define REDUCE_COL_FLOAT                                                \
  for (inX = 0; inX < numLoop; inX++)                                   \
    sum += trigTable[inX] * tex2D(sUnpadFloatTex, inX + addXoffset, iySrc); \
  break;

// Kernel to reduce portion of full image in X from different types of data, outputting
// floats through a surface reference
__global__ void reduceColumns(int binning, float *linearArr, int type, int nxOut,
                              int nyOut, int iyStart, int delXstart, int delXend) 
{
  int ixSrc, iySrc, isum = 0;
  int addXoffset, inX, numLoop = delXend + 1 - delXstart;
  float sum = 0.;
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix >= nxOut || iy >= nyOut)
    return;
  
  // ixSrc is a source coordinate in full image relative to the mXstart and delX needs
  // to be added to get to the absolute source coord
  // iySrc is a absolute source coordinate in full image, and iyStart is mYstart minus 
  // half the support (kernel radius)
  ixSrc = ix * binning;
  iySrc = iy + iyStart;
  
  addXoffset = ixSrc + delXstart;
  if (numLoop >= 40) {
    switch (type) {
    case MRC_MODE_BYTE:
#pragma unroll 40
      REDUCE_COL_BYTE;
#pragma unroll 40
      REDUCE_COL_SHORT;
#pragma unroll 40
      REDUCE_COL_USHORT;
#pragma unroll 40
      REDUCE_COL_FLOAT;
    }
  } else if (numLoop >= 32) {
    switch (type) {
    case MRC_MODE_BYTE:
#pragma unroll 32
      REDUCE_COL_BYTE;
#pragma unroll 32
      REDUCE_COL_SHORT;
#pragma unroll 32
      REDUCE_COL_USHORT;
#pragma unroll 32
      REDUCE_COL_FLOAT;
    }
  } else if (numLoop >= 24) {
    switch (type) {
    case MRC_MODE_BYTE:
#pragma unroll 24
      REDUCE_COL_BYTE;
#pragma unroll 24
      REDUCE_COL_SHORT;
#pragma unroll 24
      REDUCE_COL_USHORT;
#pragma unroll 24
      REDUCE_COL_FLOAT;
    }
  } else if (numLoop >= 16) {
    switch (type) {
    case MRC_MODE_BYTE:
#pragma unroll 16
      REDUCE_COL_BYTE;
#pragma unroll 16
      REDUCE_COL_SHORT;
#pragma unroll 16
      REDUCE_COL_USHORT;
#pragma unroll 16
      REDUCE_COL_FLOAT;
    }
  } else if (numLoop >= 8) {
    switch (type) {
    case MRC_MODE_BYTE:
#pragma unroll 8
      REDUCE_COL_BYTE;
#pragma unroll 8
      REDUCE_COL_SHORT;
#pragma unroll 8
      REDUCE_COL_USHORT;
#pragma unroll 8
      REDUCE_COL_FLOAT;
    }
  } else {
    switch (type) {
    case MRC_MODE_BYTE:
      REDUCE_COL_BYTE;
      REDUCE_COL_SHORT;
      REDUCE_COL_USHORT;
      REDUCE_COL_FLOAT;
    }
  }
  
#ifdef NO_SURFACES
  linearArr[ix + iy * nxOut] = sum;
#else
  surf2Dwrite(sum, sTempSurfRef, ix * sizeof(float), iy);
#endif
}

// Second kernel to reduce the X-reduced image in Y.  This is much quicker (36 ms vs 
// 221 ms for X reduction) and the unroll actuaaly cost a bit of time
__global__ void reduceRowsTaperPad(float *outArr, int binning, int inYstart,
                                   int ixLow, int ixHigh, int iyLow, int iyHigh,
                                   int nxDimOut, int alignXpad, int alignYpad, 
                                   int nxTaper, int nyTaper,
                                   float dmean, int delYstart, int delYend)
{
  int outInd, ixSrc, iySrc;
  int addYoffset, inY, numLoop = delYend + 1 - delYstart;
  float fracX, fracY, fmin, sum;
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix >= alignXpad || iy >= alignYpad)
    return;
  outInd = ix + nxDimOut * iy;
  if (ix < ixLow || ix >= ixHigh || iy < iyLow || iy >= iyHigh) {
    outArr[outInd] = dmean;
    return;
  }

  // ixSrc is absolute coordinate in X-reduced image
  // iySrc is location in full image relative to mYstart there, adjusted downward by
  // start of X-reduced image in Y (mYstart - half-support), so when delY is added it
  // will be location in X-reduced image
  ixSrc = ix - ixLow;
  iySrc = (iy - iyLow) * binning - inYstart;

  addYoffset = iySrc + delYstart;
  sum = 0;
  for (inY = 0; inY < numLoop; inY++)
    sum += trigTable[inY] * tex2D(sUnpadFloatTex, ixSrc, inY + addYoffset);
  fracX = 1.;
  fracY = 1.;
  if (ix < nxTaper + ixLow)
    fracX = (ix + 1.f - ixLow) / (nxTaper + 1.);
  else if (ix >= ixHigh - nxTaper)
    fracX = (ixHigh - ix) / (nxTaper + 1.);
  if (iy < nyTaper + iyLow)
    fracY = (iy + 1.f - iyLow) / (nyTaper + 1.);
  else if (iy >= iyHigh - nyTaper)
    fracY = (iyHigh - iy) / (nyTaper + 1.);
  if (fracX < 1 || fracY < 1.) {
    fmin = min(fracX, fracY);
    outArr[outInd] = fmin * (sum - dmean) + dmean;
  } else {
    outArr[outInd] = sum;
  }
}

// Kernel to simply trim, pad, and taper inside when there is no reduction
__global__ void trimTaperPad(float *outArr, int type, int ixLow, int ixHigh, int iyLow,
                             int iyHigh, int nxDimOut, int alignXpad, int alignYpad,
                             int nxTaper, int nyTaper, float dmean)
{
  int outInd, ixSrc, iySrc;
  float fracX, fracY, fmin, sum;
  int ix = blockIdx.x * blockDim.x + threadIdx.x;
  int iy = blockIdx.y * blockDim.y + threadIdx.y;
  if (ix >= alignXpad || iy >= alignYpad)
    return;
  outInd = ix + nxDimOut * iy;
  if (ix < ixLow || ix >= ixHigh || iy < iyLow || iy >= iyHigh) {
    outArr[outInd] = dmean;
    return;
  }
  ixSrc = (ix + ixLow);
  iySrc = (iy + iyLow);

  // Does reduction on CPU reach outside the subarea being computed?
  switch (type) {
  case MRC_MODE_BYTE:
    sum = tex2D(sUnpadByteTex, ixSrc, iySrc);
    break;
  case MRC_MODE_SHORT:
    sum = tex2D(sUnpadShortTex, ixSrc, iySrc);
    break;
  case MRC_MODE_USHORT:
    sum = tex2D(sUnpadUShortTex, ixSrc, iySrc);
    break;
  default:
    sum = tex2D(sUnpadFloatTex, ixSrc, iySrc);
    break;
  }

  fracX = 1.;
  fracY = 1.;
  if (ix < nxTaper + ixLow)
    fracX = (ix + 1.f - ixLow) / (nxTaper + 1.);
  else if (ix >= ixHigh - nxTaper)
    fracX = (ixHigh - ix) / (nxTaper + 1.);
  if (iy < nyTaper + iyLow)
    fracY = (iy + 1.f - iyLow) / (nyTaper + 1.);
  else if (iy >= iyHigh - nyTaper)
    fracY = (iyHigh - iy) / (nyTaper + 1.);
  if (fracX < 1 || fracY < 1.) {
    fmin = min(fracX, fracY);
    outArr[outInd] = fmin * (sum - dmean) + dmean;
  } else {
    outArr[outInd] = sum;
  }
}


/////////////////////////////////////////////////////////////////////////////
// THE FrameGPU CLASS
/////////////////////////////////////////////////////////////////////////////

FrameGPU::FrameGPU()
{
  sChanDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
  sByteChanDesc = cudaCreateChannelDesc(8, 0, 0, 0, cudaChannelFormatKindUnsigned);
  sShortChanDesc = cudaCreateChannelDesc(16, 0, 0, 0, cudaChannelFormatKindSigned);
  sUShortChanDesc = cudaCreateChannelDesc(16, 0, 0, 0, cudaChannelFormatKindUnsigned);
  clearAllItems();
  // framealign calls zeroTimers which is implemented in .h file
}

/*
 * complete cleanup of both functionalities
 */
void FrameGPU::cleanup()
{
  B3DFREE(mXshiftTrig);
  B3DFREE(mYshiftTrig);
  cleanAlignItems();
  cleanSumItems();
  cleanPreProc();
  clearAllItems();
}

/*
 * Clean up and free resources for alignment
 */
void FrameGPU::cleanAlignItems()
{
  unbindVariableBindings();
  destroyPlan(sAlignForwardPlan);
  destroyPlan(sAlignInversePlan);
  freeCudaStack(mSavedBinPad);
  freeCudaStack(mSavedGroups);
  freeCudaArray(&mWorkBinPad);
  freeCudaArray(&mCorrBinPad);
  freeCudaArray(&mAlignSum);
  freeCudaArray(&mSubareaCorr);
  freeCudaArray(&mRealCorr);
  if (mFiltMask) {
    freeCudaArray(&mFiltMask);
    cudaUnbindTexture(sMaskTex);
  }
  B3DFREE(mHostSubarea);
  free2Darray(&sReducedInXarray);
  freeCudaArray(&mReducedInXlinear);
  mRedColX = mRedColY = 0;
  mAlignXpad = mAlignYpad = 0;
  mAliFiltSize = 0;
  mUnpaddedBytes = 0;
  mDoAlignBinPad = 0;
}

/*
 * Clean up and free resources for summing
 */
void FrameGPU::cleanSumItems()
{
  unbindVariableBindings();
  freeCudaArray(&mWorkFullSize);
  freeCudaArray(&mEvenSum);
  freeCudaArray(&mOddSum);
  freeCudaArray(&mNonDWsum);
  if (mDoseWgtFilter)
    cudaUnbindTexture(sDWFilterTex);
  for (int ind = 0; ind < (int)sSavedUnpadded.size(); ind++)
    free2Darray(&sSavedUnpadded[ind]);
  sSavedUnpadded.clear();
  mSavedFrameNums.clear();
  mUnpadEdgeMeans.clear();
  mNumOnUnpadStack = 0;
  freeCudaArray(&mDoseWgtFilter);
  destroyPlan(sFullForwardPlan);
  destroyPlan(sSumInversePlan);
  mFullXpad = mFullYpad = 0;
  mSumXpad = mSumYpad = 0;
  mDoNoiseTaper = 0;
}

/*
 * Free up the two bindings that can be set to different arrays
 */
void FrameGPU::unbindVariableBindings()
{
  if (mBoundToFull)
    cudaUnbindTexture(sFullTex);
  if (mBoundSum >= 0)
    cudaUnbindTexture(sSumTex);
  mBoundToFull = NULL;
  mBoundSum = -1;
}

/*
 * Cleans up pre-processing arrays
 */
void FrameGPU::cleanPreProc()
{
  if (sGainRefArray) {
    cudaUnbindTexture(sGainRefTex);
    free2Darray(&sGainRefArray);
    sGainRefArray = NULL;
  }
  if (sDefectMapArray) {
    cudaUnbindTexture(sDefectTex);
    free2Darray(&sDefectMapArray);
  }
  free2Darray(&sTempFloatArray);
  freeCudaArray(&mProcessedLinear);
  mNxGain = 0;
  mNyGain = 0;
  mCamSizeX = 0;
  mCamSizeY = 0;
  mTruncLimit = 0.;
  mTempFloatSizeX = 0;
  mTempFloatSizeY = 0;
}

/*
 * Initializes all pointers and tested member variables on cleanup or construction
 */
void FrameGPU::clearAllItems()
{
  mWorkFullSize = NULL;
  mEvenSum = NULL;
  mOddSum = NULL;
  mNonDWsum = NULL;
  mDoseWgtFilter = NULL;
  mXshiftTrig = NULL;
  mYshiftTrig = NULL;
  mWorkBinPad = NULL;
  mCorrBinPad = NULL;
  mFiltMask = NULL;
  mAlignSum = NULL;
  mSubareaCorr = NULL;
  mRealCorr = NULL;
  mBoundToFull = NULL;
  mHostSubarea = NULL;
  mReducedInXlinear = NULL;
  mProcessedLinear = NULL;
#ifdef NO_SURFACES
  mNoSurfaces = 1;
#else
  mNoSurfaces = 0;
#endif
  mAliFiltSize = 0;
  mFullXpad = mFullYpad = 0;
  mSumXpad = mSumYpad = 0;
  mAlignXpad = mAlignYpad = 0;
  mBoundSum = -1;
  mXtrigSize = 0;
  mYtrigSize = 0;
  mGroupSize = 1;
  mDWFilterDelta = 0.;
  mDWFilterSize = 0;
  mNxGain = 0;
  mNyGain = 0;
  mCamSizeX = 0;
  mCamSizeY = 0;
  mTruncLimit = 0.;
  mTempFloatSizeX = 0;
  mTempFloatSizeY = 0;
  free2Darray(&sTempRawArray);
  mNxRaw = 0;
  mNyRaw = 0;
  mNumOnUnpadStack = 0;
  mRedColX = mRedColY = 0;
}

// Free one regular array, set its pointer to NULL
void FrameGPU::freeCudaArray(float **array)
{
  if (*array)
    cudaFree(*array);
  *array = NULL;
}

/* convenience function to free a stack of regular arrays */
void FrameGPU::freeCudaStack(std::vector<float *> &saved)
{
  for (int ind = 0; ind < saved.size(); ind++)
    freeCudaArray(&saved[ind]);
  saved.clear();
}

// Free a 2D array and set its pointer to null
static void free2Darray(cudaArray **array)
{
  if (*array)
    cudaFreeArray(*array);
  *array = NULL;
}

// Destroy an FFT plan and set to 0
static void destroyPlan(cufftHandle &plan)
{
  if (plan)
    cufftDestroy(plan);
  plan = 0;
}

/*
 * Test whether a GPU is available, either a GPU of the given number if nGPU is
 * > 0, or the one with the best processing rate if nGPU is 0, and return the
 * memory in bytes.  Return value is 1 for success, 0 for failure.
 */
int FrameGPU::gpuAvailable(int nGPU, float *memory, int debug)
{
  int current_device = 0;
  int device_count = 0;
  int totalCores;
  float gflops;
  struct cudaDeviceProp device_properties, best_properties;

  // The Mac mini comes through with a clock rate of 0 so allow a 0 product
  float max_gflops = -1.;
  mDebug = debug;
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
      mMax_gflops_device = current_device;
      best_properties = device_properties;
    }
  }
    
  if (mMax_gflops_device >= 0) {
    *memory = best_properties.totalGlobalMem;
    if (cudaSetDevice(mMax_gflops_device) != cudaSuccess) {
      pflerr("Error selecting GPU device %d", mMax_gflops_device + 1);
      return 0;
    }
    mDeviceSelected = 1;
    return 1;
  }
  return 0;

}

/*
 * Set parameters and arrays for pre-processing if any.
 * This needs to be called (only) on first frame of a set since processing parameters
 * are not available in framealign::initialize
 */
int FrameGPU::setPreProcParams(float *gainRef, int nxGain, int nyGain, float truncLimit,
                               unsigned char *defectMap, int camSizeX, int camSizeY)
{
  int err = 0;

  // Set flags for actions to do based on what has been passed in
  mCorrectDefects = (defectMap != NULL && camSizeX > 0)  ? 1 : 0;
  mDoGainNorm = (gainRef != NULL && nxGain > 0) ? 1 : 0;
  mDoPreprocess = (mCorrectDefects || mDoGainNorm || truncLimit > 0) ? 1 : 0;

  // Cleanup if anything has changed
  if (!mDoPreprocess || nxGain != mNxGain || nyGain != mNyGain || 
      (mNxGain > 0 && !mDoGainNorm) || camSizeX != mCamSizeX || camSizeY != mCamSizeY || 
      (mCamSizeX > 0 && !mCorrectDefects))
    cleanPreProc();

  // Setup gain normalization
  if (mDoGainNorm && !mNxGain) {
    if (testErrCode(cudaMallocArray(&sGainRefArray, &sChanDesc, nxGain, nyGain),
                    "allocate gain reference array on GPU", 0)) {
      err = 1;
    } else {
    
      // Bind texture
      SETUP_TEXTURE(sGainRefTex);
      if (testErrCode(cudaBindTextureToArray(sGainRefTex, sGainRefArray, sChanDesc),
                      "bind gain reference array to texture", 0)) {
        err = 1;
      } else {
        mNxGain = nxGain;
        mNyGain = nyGain;

        // Copy to array
        if (testErrCode(cudaMemcpyToArray(sGainRefArray, 0, 0, gainRef, nxGain * nyGain * 
                                          sizeof(float), cudaMemcpyHostToDevice),
                        "copy gain reference to GPU", 0))
          err = 1;
      }
    }
    //dumpUnpadArray(sGainRefArray, nxGain, nyGain, MRC_MODE_FLOAT, "gain reference");
  }

  // Set up defect correction
  if (!err && mCorrectDefects && camSizeX != mCamSizeX) {
    if (testErrCode(cudaMallocArray(&sDefectMapArray, &sByteChanDesc, mUnpaddedX, 
                                    mUnpaddedY) , "allocate defect map array on GPU", 0)){
      err = 1;
    } else {

      // Bind texture
      SETUP_TEXTURE(sDefectTex);
      if (testErrCode(cudaBindTextureToArray(sDefectTex, sDefectMapArray, sByteChanDesc),
                      "bind defect map array to texture", 0)) {
        err = 1;
      } else {
        mCamSizeX = camSizeX;
        mCamSizeY = camSizeY;
        
        // Copy to array
        if (testErrCode(cudaMemcpyToArray(sDefectMapArray, 0, 0, defectMap, mUnpaddedX * 
                                          mUnpaddedY, cudaMemcpyHostToDevice),
                        "copy defect map to GPU", 0))
          err = 1;
      }
    }
  }

  // Need a 2D temp float array if processing at all
  if (!err && mDoPreprocess && 
      (!sTempFloatArray || (mNoSurfaces && !mProcessedLinear) || 
       mUnpaddedX != mTempFloatSizeX || mUnpaddedY != mTempFloatSizeY)) {
    free2Darray(&sTempFloatArray);
    freeCudaArray(&mProcessedLinear);
#ifdef NO_SURFACES
    if (cudaMallocArray(&sTempFloatArray, &sChanDesc, mUnpaddedX, mUnpaddedY) != 
        cudaSuccess || cudaMalloc(&mProcessedLinear, mUnpaddedX * mUnpaddedY * 
                                  sizeof(float)) != cudaSuccess) {
      pflerr("allocate temp float 2D array or linear array on GPU");
      err = 1;
#else
    if (testErrCode(cudaMallocArray(&sTempFloatArray, &sChanDesc, mUnpaddedX, mUnpaddedY,
                                    cudaArraySurfaceLoadStore),
                    "allocate temp float 2D array on GPU", 0)) {
      err = 1;
#endif
    } else {
      mTempFloatSizeX = mUnpaddedX;
      mTempFloatSizeY = mUnpaddedY;
    }
  }
  mTruncLimit = truncLimit;

  // An error at this stack can be handled ny simply canceling all the initial GPU 
  // operations
  if (err) {
    cleanPreProc();
    mDoNoiseTaper = mDoAlignBinPad = mDoPreprocess = mStackUnpadded = 0;
  }
  return err;
}

/*
 * Allocate a 2D array of the given type and size
 */
static cudaError_t malloc2Darray(cudaArray **arrayPtr, int type, int sizeX, int sizeY)
{
  switch (type) {
  case MRC_MODE_BYTE:
    return(cudaMallocArray(arrayPtr, &sByteChanDesc, sizeX, sizeY));
  case MRC_MODE_SHORT:
    return(cudaMallocArray(arrayPtr, &sShortChanDesc, sizeX, sizeY));
  case MRC_MODE_USHORT:
    return(cudaMallocArray(arrayPtr, &sUShortChanDesc, sizeX, sizeY));
  case MRC_MODE_FLOAT:
    return(cudaMallocArray(arrayPtr, &sChanDesc, sizeX, sizeY));
  }
  return cudaErrorUnknown;
}

/*
 * Bind one of the unpadded 2D arrays to a texture of the right type
 */
static cudaError_t bindUnpadArray(cudaArray *arrayPtr, int type)
{
  switch (type) {
  case MRC_MODE_BYTE:
    SETUP_TEXTURE(sUnpadByteTex);
    return(cudaBindTextureToArray(sUnpadByteTex, arrayPtr));
  case MRC_MODE_SHORT:
    SETUP_TEXTURE(sUnpadShortTex);
    return(cudaBindTextureToArray(sUnpadShortTex, arrayPtr));
  case MRC_MODE_USHORT:
    SETUP_TEXTURE(sUnpadUShortTex);
    return(cudaBindTextureToArray(sUnpadUShortTex, arrayPtr));
  case MRC_MODE_FLOAT:
    SETUP_TEXTURE(sUnpadFloatTex);
    return(cudaBindTextureToArray(sUnpadFloatTex, arrayPtr));
  }
  return cudaErrorUnknown;
}

/*
 * Unbind an unpadded 2D array of the given type
 */
void FrameGPU::unbindUnpadArray(int type2d)
{
  switch (type2d) {
  case MRC_MODE_BYTE:
    cudaUnbindTexture(sUnpadByteTex);
    break;
  case MRC_MODE_SHORT:
    cudaUnbindTexture(sUnpadShortTex);
    break;
  case MRC_MODE_USHORT:
    cudaUnbindTexture(sUnpadUShortTex);
    break;
  case MRC_MODE_FLOAT:
    cudaUnbindTexture(sUnpadFloatTex);
    break;
  }
}

/*
 * Set parameters for the reduction tapering and the length parameter for noise tapering
 * This can be called on all frames of a set and should always be called, regardless of
 * whether doing binpad
 */
void FrameGPU::setBinPadParams(int xstart, int xend, int ystart, int yend, int binning,
                               int nxTaper, int nyTaper, int type, int filtType,
                               int noiseLen)
{
  mXstart = xstart;
  mXend = xend;
  mYstart = ystart;
  mYend = yend;
  mAliBinning = binning;
  mNxTaper = nxTaper;
  mNyTaper = nyTaper;
  mStackType = type;
  mNoiseLength = noiseLen;
  mAntiFiltType = filtType;
}

/*
 * Allocate a temporary 2D array for raw data of the current "stack" type when there
 * is no stack or the stack is not big enough; and also an X-reduced array when
 * doing reduction - padding,
 * This has to get called by summing and aligning as it is needed independently by
 * those two operations
 */
int FrameGPU::manageRawTempArray(int aligning)
{
  int support, redColY;
  int redColX = (mXend + 1 - mXstart) / mAliBinning;

  // Need a "raw" temporary array if not stacking 
  if ((mDoAlignBinPad || mDoNoiseTaper) && (!mStackUnpadded || mStackIsLimited)) {
    if (!sTempRawArray || mNxRaw != mUnpaddedX || mNyRaw != mUnpaddedY) {
      free2Darray(&sTempRawArray);
      if (testErrCode(malloc2Darray(&sTempRawArray, mStackType, mUnpaddedX, mUnpaddedY), 
                      "allocate raw 2D temp array on GPU", 1))
        return 1;
      mNxRaw = mUnpaddedX;
      mNyRaw = mUnpaddedY;
    }
  } else {
    free2Darray(&sTempRawArray);
    mNxRaw = mNyRaw = 0;
  }

  // Manage column-reduced array(s) for reducing
  if (mDoAlignBinPad && mAliBinning > 1) {
    if (aligning) {
      selectZoomFilter(mAntiFiltType, 1. / mAliBinning, &support);
      redColY = (mYend + 1 - mYstart) + 2 * (support / 2 + 1);
      if (!sReducedInXarray || mRedColX != redColX || mRedColY != redColY ||
          (mNoSurfaces && !mReducedInXlinear)) {
        free2Darray(&sReducedInXarray);
        freeCudaArray(&mReducedInXlinear);
#ifdef NO_SURFACES
        if (testErrCode(cudaMallocArray(&sReducedInXarray, &sChanDesc, redColX, redColY),
                        "allocate row reduction 2D array on GPU", 1))
          return 1;
        if (testErrCode(cudaMalloc(&mReducedInXlinear, redColX * redColY * sizeof(float)),
                        "allocate row reduction output array on GPU", 1))
          return 1;
#else
        if (testErrCode(cudaMallocArray(&sReducedInXarray, &sChanDesc, redColX, redColY,
                                        cudaArraySurfaceLoadStore),
                        "allocate row reduction 2D array on GPU", 1))
          return 1;
#endif
        mRedColX = redColX;
        mRedColY = redColY;
      }
    }
  } else {
    free2Darray(&sReducedInXarray);
    mRedColX = mRedColY = 0;
    freeCudaArray(&mReducedInXlinear);
  }
  return 0;
}

/*
 * Set some basic parameters for the initial processing steps
 * This needs to be called before setup of align OR summing; summiong setup is
 * done in two places in framealign
 */
void FrameGPU::setUnpaddedSize(int unpadX, int unpadY, int flags, int debug)
{
  mUnpaddedX = unpadX;
  mUnpaddedY = unpadY;
  mDoNoiseTaper = (flags & GPU_DO_NOISE_TAPER) ? 1 : 0;
  mDoAlignBinPad = (flags & GPU_DO_BIN_PAD) ? 1 : 0;
  mStackUnpadded = (flags & STACK_FULL_ON_GPU) ? 1 : 0;
  mStackIsLimited = (flags & GPU_STACK_LIMITED) ? 1 : 0;
  mDebug = debug % 10;
  mTrackTime = (debug / 10) % 10;
}

/*
 * Set parameters and initialize for summing aligned images into one or two arrays,
 * allocating new arrays as needed
 */
int FrameGPU::setupSumming(int fullXpad, int fullYpad, int sumXpad,
                           int sumYpad, int evenOdd)
{
  //size_t workSize;
  int error;
  int nonDW = evenOdd & 2;
  bool sumChanged = sumXpad != mSumXpad || sumYpad != mSumYpad;
  evenOdd &= 1;
  if (!mDeviceSelected)
    return -1;

  // Do not clean up on any failures, so that deferred summing can recover the stack
  // proceed on CPU.  Initial framealign setupSumming does full cleanup

  // Manage work array if size has changed
  mFullBytes = (fullXpad + 2) * fullYpad * sizeof(float);
  if (fullXpad != mFullXpad || fullYpad != mFullYpad) {
    freeCudaArray(&mWorkFullSize);
    if (testErrCode(cudaMalloc((void **)&mWorkFullSize, mFullBytes),
                    "allocate full work array on GPU", 0))
      return 1;
  }

  // Manage even sum array if size has changed
  mSumBytes = (sumXpad + 2) * sumYpad * sizeof(float);
  if (sumChanged) {
    freeCudaArray(&mEvenSum);
    if (testErrCode(cudaMalloc((void **)&mEvenSum, mSumBytes),
                    "allocate array for main sum on GPU", 0))
      return 1;

    SETUP_TEXTURE(sSumTex);
    if (bindSumArray(0)) {
      return 1;
    }
  }

  // Manage odd sum array if size has changed or need for it has changed
  if (sumChanged || evenOdd != mDoEvenOdd) {
    freeCudaArray(&mOddSum);
    if (evenOdd && testErrCode(cudaMalloc((void **)&mOddSum, mSumBytes),
                               "allocate array for odd sum on GPU", 0))
      return 1;
  }

  // Manage non-dose weight sum array if size or need for it has changed
  if (sumChanged || nonDW != mDoUnDWsum) {
    freeCudaArray(&mNonDWsum);
    if (nonDW && testErrCode(cudaMalloc((void **)&mNonDWsum, mSumBytes),
                             "allocate array for unweighted sum on GPU", 0))
      return 1;
  }

  // Manage FFT plans;
  destroyPlan(sFullForwardPlan);
  destroyPlan(sSumInversePlan);
  //START_TIMER;
  /* Not available in CUDA 4:
  size_t workSize;
  cufftEstimate2d(fullYpad, fullXpad, CUFFT_R2C, &workSize);
  utilPrint("Work size estimate for forward  %u\n", workSize); */
  error = cufftPlan2d(&sFullForwardPlan, fullYpad, fullXpad, CUFFT_R2C);
  if (error != CUFFT_SUCCESS) {
    utilPrint("Failed to make plan for full forward FFT (error %d)\n", error);
    return 1;
  }
  //if (mDebug)
  //PRINT2("plan time: ", wallTime() - mWallStart);
  /*cufftGetSize2d(sFullForwardPlan, fullYpad, fullXpad, CUFFT_R2C, &workSize);
    utilPrint("Work size for forward plan %u\n", workSize);*/

  if (manageShiftTrigs(sumXpad, sumYpad)) {
    return 1;
  }

  // Clear sum arrays
  if (cudaMemset(mEvenSum, 0, mSumBytes) != cudaSuccess || 
      (evenOdd && cudaMemset(mOddSum, 0, mSumBytes) != cudaSuccess) || 
      (nonDW && cudaMemset(mNonDWsum, 0, mSumBytes) != cudaSuccess)) {
    pflerr("Failed to zero out sum array on GPU");
    return 1;
  }
  
  mFullXpad = fullXpad;
  mFullYpad = fullYpad;
  mSumXpad = sumXpad;
  mSumYpad = sumYpad;
  mDoEvenOdd = evenOdd;
  mDoUnDWsum = nonDW;
  mNumFramesSummed = 0;
  mNumAlignedFrames = 0;
  mDWFilterSize = 0;
  return 0;
}

/*
 * Set parameters and initialize for aligning by a particular strategy, allocating the
 * appropriate arrays
 */
int FrameGPU::setupAligning(int alignXpad, int alignYpad, int sumXpad, int sumYpad,
                            float *alignMask, int aliFiltSize, int groupSize, 
                            int expectStackSize, int doAlignSum)
{
  int error;
  int xpad = sumXpad;
  int ypad = sumYpad;

  if (doAlignSum || !expectStackSize) {
    xpad = B3DMAX(sumXpad, alignXpad);
    ypad = B3DMAX(sumYpad, alignYpad);
  }
  if (xpad && manageShiftTrigs(xpad, ypad)) {
    cleanup();
    return 1;
  }

  // Manage align arrays etc if size has changed
  if (alignXpad != mAlignXpad || alignYpad != mAlignYpad || groupSize != mGroupSize ||
      expectStackSize != mExpectStackSize || doAlignSum != mDoAlignSum) {

    // Cleanup existing
    error = mDoAlignBinPad;
    cleanAlignItems();
    mDoAlignBinPad = error;
    
    // Allocate arrays that are needed
    mAlignBytes = (alignXpad + 2) * alignYpad * sizeof(float);
    if (cudaMalloc((void **)&mWorkBinPad, mAlignBytes) != cudaSuccess ||
        cudaMalloc((void **)&mFiltMask, mAlignBytes) != cudaSuccess ||
        ((doAlignSum || !expectStackSize) && 
         (cudaMalloc((void **)&mAlignSum, mAlignBytes) != cudaSuccess ||
          cudaMalloc((void **)&mCorrBinPad, mAlignBytes) != cudaSuccess)) ||
        cudaMalloc((void **)&mRealCorr, alignXpad * alignYpad * sizeof(float)) !=
        cudaSuccess) {
      pflerr("Failed to allocate arrays for aligning on GPU");
      cleanup();
      return 1;
    }

    // Make plans
    error = cufftPlan2d(&sAlignForwardPlan, alignYpad, alignXpad, CUFFT_R2C);
    if (error == CUFFT_SUCCESS)
      error = cufftPlan2d(&sAlignInversePlan, alignYpad, alignXpad, CUFFT_C2R);
    if (error != CUFFT_SUCCESS) {
      utilPrint("Failed to make plan for align FFTs (error %d)\n", error);
      cleanup();
      return 1;
    }
    
    // Copy mask array
    START_TIMER;
    if (testErrCode(cudaMemcpy(mFiltMask, alignMask, mAlignBytes, cudaMemcpyHostToDevice),
                    "copy filter mask to GPU array", 1))
      return 1;
    ADD_TIME(mWallCopy);
    
    // Bind filter mask
    sMaskTex.filterMode = cudaFilterModePoint;
    sMaskTex.normalized = false;
    if (testErrCode(cudaBindTexture(NULL, sMaskTex, mFiltMask, sChanDesc, mAlignBytes),
                    "bind filter mask array to texture", 1))
      return 1;
  }
    
  // Get arrays for extracting subarea from correlation, both on device and host make
  // it a multiple of a nice size
  if (!mSubareaCorr || !mHostSubarea || aliFiltSize != mAliFiltSize) {
    mBigSubareaSize = NICE_GPU_DIVISOR * (aliFiltSize / NICE_GPU_DIVISOR + 1);
    freeCudaArray(&mSubareaCorr);
    B3DFREE(mHostSubarea);
    if (testErrCode(cudaMalloc((void **)&mSubareaCorr, 
                               mBigSubareaSize * mBigSubareaSize * sizeof(float)),
                    "allocate array for subarea of correlation on GPU", 1))
      return 1;
    mHostSubarea = B3DMALLOC(float, mBigSubareaSize * mBigSubareaSize);
    if (!mHostSubarea) {
      utilPrint("Failed to allocate memory for oversized subarea\n");
      cleanup();
      return 1;
    }
  }
  
  // Clear sum arrays
  if ((doAlignSum || !expectStackSize) && clearAlignSum()) {
    cleanup();
    return 1;
  }
  
  mAlignXpad = alignXpad;
  mAlignYpad = alignYpad;
  mExpectStackSize = expectStackSize;
  mDoAlignSum = doAlignSum;
  mAliFiltSize = aliFiltSize;
  mGroupSize = groupSize;
  mNumOnUnpadStack = 0;
  mNumFramesSummed = 0;
  mNumAlignedFrames = 0;
  return 0;
}

/*
 * Set the dose weighting filter to use on the next final summing operation.  This can
 * be called before the addToSums call even when not dose weighting, although the filter
 * size is initialized to 0 when setting up
 */
int FrameGPU::setupDoseWeighting(float *filter, int filtSize, float delta)
{
  std::vector<float> filtSubset;
  int numBytes;
  mDWFilterSize = filtSize;
  if (!filtSize) {
    mDWFilterDelta = 0.;
    return 0;
  }

  mDWFilterDelta = delta;
  numBytes = mDWFilterSize * sizeof(float);
  if (!mDoseWgtFilter) {
    if (testErrCode(cudaMalloc((void **)&mDoseWgtFilter, numBytes),
                    "allocate array for dose weight filter on GPU", 0))
      return 1;

    sDWFilterTex.addressMode[0] = cudaAddressModeClamp;
    sDWFilterTex.addressMode[1] = cudaAddressModeClamp;
    sDWFilterTex.filterMode = cudaFilterModeLinear;
    sDWFilterTex.normalized = false;
    if (testErrCode(cudaBindTexture(NULL, sDWFilterTex, mDoseWgtFilter, sChanDesc,
                                    numBytes),
                    "bind dose weight filter array to texture", 0))
      return 1;
  }

  // Copy the filter
  if (testErrCode(cudaMemcpy(mDoseWgtFilter, filter, numBytes, cudaMemcpyHostToDevice),
                  "copy dose weight filter to GPU array", 0))
    return 1;
  return 0;
}

/*
 * Zero out the sum for adding up an alignment reference
 */
int FrameGPU::clearAlignSum()
{
  if (testErrCode(cudaMemset(mAlignSum, 0, mAlignBytes),
                  "zero out align sum array on GPU", 0))
    return 1;
  return 0;
}

/*
 * Take the FFT of a full image and shift add it to the appropriate sum
 */
int FrameGPU::addToFullSum(float *fullArr, float shiftX, float shiftY)
{
  float *sumArr = mEvenSum;
  int err, dataSize, ind, needBound = 0;
  cudaArray *dev2dArr;
  int type2d = mStackType;
  float edgeMean;
  int blockX = 16;
  static int seed = 123456;

  if (manageRawTempArray(0))
    return 1;
  dev2dArr = sTempRawArray;

  // Select array
  if (mDoEvenOdd && (mNumFramesSummed % 2) != 0) {
    sumArr = mOddSum;
    needBound = 1;
  }

  /* Data flow for doing noise pad on GPU:
         No stack                         Stack
     load to sTempRawArray       get from raw stack array
     
                         Need preproc
     proc to sTempFloatArray     proc to sTempFloatArray  
  */
  
  if (!fullArr && (!mNumOnUnpadStack || mSavedFrameNums[0] != mNumFramesSummed)) {
    utilPrint("A NULL array was passed to addToFullSum but frame %d instead of %d is "
              "first on GPU stack\n", mNumOnUnpadStack ? mSavedFrameNums[0] : -1,
              mNumFramesSummed);
    return 1;
  }
  if (mDoNoiseTaper) {

    // If stacking, get array from stack
    if (!fullArr) {
      dev2dArr = sSavedUnpadded[0];
      edgeMean = mUnpadEdgeMeans[0];
      
    } else {

      // Otherwise get an edge mean and copy to raw temp array
      // Should these edge means be indented some?
      edgeMean = frameEdgeMean(fullArr, type2d, mUnpaddedX, FULL_INDENT, 
                               mUnpaddedX - FULL_INDENT - 1, mUnpaddedY - FULL_INDENT,
                               mUnpaddedY - FULL_INDENT - 1);
      dataSizeForMode(mStackType, &dataSize, &ind);
      START_TIMER;
      if (testErrCode(cudaMemcpyToArray(dev2dArr, 0, 0, fullArr, mUnpaddedX * mUnpaddedY *
                                        dataSize, cudaMemcpyHostToDevice),
                      "copy unpadded image to GPU array for noise/pad", 0))
        return 1;
      ADD_TIME(mWallCopy);
    }

    // Preprocess and change to output array and type
    if (mDoPreprocess) {
      if (runPreprocess(dev2dArr, mStackType, mNumFramesSummed))
        return 1;
      dev2dArr = sTempFloatArray;
      type2d = MRC_MODE_FLOAT;
    }

    if (testErrCode(bindUnpadArray(dev2dArr, type2d), "bind unpadded array to texture", 
                    0))
      return 1;
      
    // Do the noise pad
    START_TIMER;
    seed = (197 * (seed + 1)) & 0xFFFFF;
    dim3 blockSize(blockX, 8, 1);
    dim3 gridSize((mFullXpad + blockSize.x - 1) / blockSize.x,
                  (mFullYpad + blockSize.y - 1) / blockSize.y, 1);
    noiseTaperPad<<<gridSize, blockSize>>>
      (type2d, mUnpaddedX, mUnpaddedY, mWorkFullSize, mFullXpad + 2, mFullXpad,
       mFullYpad, mNoiseLength, 8, 20, edgeMean, seed);
    if (testReportErr("to noise pad full image"))
      return 1;
    
    if (cudaDeviceSynchronize() != cudaSuccess) {
      pflerr("Error return from synchronizing after noise padding full image");
      return 1;
    }
    ADD_TIME(mWallNoise);
    unbindUnpadArray(type2d);
    
  } else {
    
    // Copy to device
    START_TIMER;
    if (testErrCode(cudaMemcpy(mWorkFullSize, fullArr, mFullBytes, cudaMemcpyHostToDevice)
                    , "copy full padded image to GPU array", 0))
      return 1;
    ADD_TIME(mWallCopy);
  }    
  //dumpImage(mWorkFullSize, mFullXpad + 2, mFullXpad, mFullYpad, 0, "noise pad");
  
  // take FFT
  START_TIMER;
  err = cufftExecR2C(sFullForwardPlan, mWorkFullSize, (cufftComplex *)mWorkFullSize);
  if (err != CUFFT_SUCCESS) {
    utilPrint("Failure in forward full FFT on GPU (CUFFT error %d)\n", err);
    return 1;
  }
  if (mTrackTime && cudaDeviceSynchronize() != cudaSuccess) {
    pflerr("Error return from synchronizing after forward full FFT on GPU");
    return 1;
  }
  ADD_TIME(mWallFFT);
  //dumpFFT(mWorkFullSize, mFullXpad, mFullYpad, "full fft", 1);

  //dumpFFT(sumArr, mSumXpad, mSumYpad, "reduced sum", 1);
  if (shiftAddCommon(mWorkFullSize, sumArr, needBound, mFullXpad, mFullYpad,
                     mSumXpad, mSumYpad, shiftX, shiftY, 0, true))
    return 1;
                         
  mNumFramesSummed++;

  // Roll stack and reduce number on it
  if (!fullArr) {
    dev2dArr = sSavedUnpadded[0];
    for (ind = 0; ind < mNumOnUnpadStack - 1; ind++) {
      sSavedUnpadded[ind] = sSavedUnpadded[ind + 1];
      mSavedFrameNums[ind] = mSavedFrameNums[ind + 1];
      mUnpadEdgeMeans[ind] = mUnpadEdgeMeans[ind + 1];
    }
    sSavedUnpadded[ind] = dev2dArr;
    mNumOnUnpadStack--;
  }
      
  return 0;
}

/*
 * Common routine for shifting an FFT and adding it to a sum array
 */
int FrameGPU::shiftAddCommon(float *fullArr, float *sumArr, int needBound, int fullXpad,
                             int fullYpad, int sumXpad, int sumYpad, float shiftX,
                             float shiftY, int shiftSource, bool applyDoseFilt)
{
  int blockX = 32;
  int nxFFT = sumXpad / 2 + 1;
  float redFac = (float) fullXpad / (float)sumXpad;
  float dxy = -(redFac - 1) / (2. * redFac);
  float dxOut = shiftX / redFac + dxy;
  float dyOut = shiftY / redFac + dxy;
  bool reducing = sumXpad < fullXpad;
  double freq, arg;
  int ix, iy, ixSuper, iySuper, ixStart, iyStart, tableYoffset, iyInStart, ind;
  int numYdo, numXdo, byteScale, trigSize, numXsuperBlocks, numYsuperBlocks;

  // Setup all sines and cosines in host memory arrays
  START_TIMER;
  for (ix = 0; ix < nxFFT; ix++) {
    ind = 2 * ix;
    freq = 0.5 * ix / (nxFFT - 1.);
    arg = -2. * PI * freq * dxOut;
    mXshiftTrig[ind] = (float)cos(arg);
    mXshiftTrig[ind + 1] = (float)sin(arg);
  }

  for (iy = 0; iy < sumYpad; iy++) {
    ind = 2 * iy;
    freq = iy / (float)sumYpad;
    if (freq > 0.5)
      freq = freq - 1.;
    arg = -2. * PI * freq * dyOut;
    mYshiftTrig[ind] = (float)cos(arg);
    mYshiftTrig[ind + 1] = (float)sin(arg);
  }    

  // Take care of binding
  if (bindSumArray(needBound))
    return 1;
  
  if (bindFullOrCorrArray(fullArr, (fullXpad + 2 ) * fullYpad * sizeof(float)))
    return 1;

  // Determine blocking of the shift/reduction
  trigSize = 2 * (nxFFT + 1 + sumYpad);
  numXsuperBlocks = (trigSize + MAX_TABLE - 1) / MAX_TABLE;
  numYsuperBlocks = numXsuperBlocks;
  if (reducing)
    numYsuperBlocks = 2 * ((numYsuperBlocks + 1) / 2);
  iyStart = 0;
  byteScale = 2 * sizeof(float);
  tableYoffset = ((nxFFT + 1) / numXsuperBlocks);

  // Loop on the blocks
  for (iySuper = 0; iySuper < numYsuperBlocks; iySuper++) {
    if (reducing && iySuper == numYsuperBlocks / 2 - 1)
      numYdo = B3DMIN(sumYpad / 2 - iyStart, sumYpad / numYsuperBlocks);
    else
      numYdo = B3DMIN(sumYpad - iyStart, sumYpad / numYsuperBlocks);
    
    if (testErrCode(cudaMemcpyToSymbol(trigTable, &mYshiftTrig[2 * iyStart],
                                       numYdo * byteScale, byteScale * tableYoffset,
                                       cudaMemcpyHostToDevice), 
                    "copy constant data to GPU", 0))
      return 1;

    ixStart = 0;
    for (ixSuper = 0; ixSuper < numXsuperBlocks; ixSuper++) {
      numXdo = B3DMIN(nxFFT - ixStart, (nxFFT + 1) / numXsuperBlocks);
      if (testErrCode(cudaMemcpyToSymbol(trigTable, &mXshiftTrig[2 * ixStart],
                                         numXdo * byteScale, 0, cudaMemcpyHostToDevice),
                      "copy constant data to GPU", 0))
        return 1;

      // Do one superblock
      dim3 blockSize(blockX, 8, 1);
      dim3 gridSize((nxFFT + blockSize.x - 1) / blockSize.x, 
                    (sumYpad + blockSize.y - 1) / blockSize.y, 1);
      iyInStart = iyStart;
      if (reducing && iySuper >= numYsuperBlocks / 2)
        iyInStart += fullYpad - sumYpad;

      if (shiftSource)
        shiftInPlaceAddToSum<<<gridSize, blockSize>>>
          (fullArr, sumArr, ixStart, iyStart, numXdo, numYdo, nxFFT, 2 * tableYoffset);
      else
        shiftAndAddToSum<<<gridSize, blockSize>>>
          (sumArr, mDoUnDWsum ? mNonDWsum : NULL, ixStart, iyInStart, iyStart, numXdo,
           numYdo, fullXpad / 2 + 1, nxFFT, 2 * tableYoffset, 
           applyDoseFilt ? mDWFilterSize : 0, mDWFilterDelta, 1.f / fullXpad,
           1.f / fullYpad); 
      if (testReportErr("to add to sum")) {
        return 2;
      }
      if (cudaDeviceSynchronize() != cudaSuccess) {
        pflerr("Error return from synchronizing after shift and add block %d %d", ixSuper,
               iySuper);
        return 1;
      }
      ixStart += numXdo;
    }
    iyStart += numYdo;
  }
  //dumpFFT(sumArr, sumXpad, sumYpad, "fft of sum", 0);
  ADD_TIME(mWallShift);
  return 0;
}

/*
 * Return summed FFTs, add them and return real sum
 */
int FrameGPU::returnSums(float *sumArr, float *evenArr, float *oddArr, int evenOddOnly)
{
  int blockX = 32;
  int err, error = 0;
  int sumXplus = mSumXpad + 2;
  // FFT scaling consists of the standard forward fft scaling, divided by the reduction 
  // factor as in fourierReduceImage
  // Image scaling then includes the standard FFT scaling for the inverse FFT
  float fftScale = 1. / (sqrt((double)mFullXpad * mFullYpad) * (float)mFullYpad/mSumYpad);
  float imScale = (fftScale / sqrt((double)mSumXpad * mSumYpad));

  // Copy the even/odd arrays back if they exist
  // Not sure this is right if emergency call is made and there are no odd sums
  START_TIMER;
  if (mDoEvenOdd && evenArr && oddArr) {
    if (cudaMemcpy(evenArr, mEvenSum, mSumBytes, cudaMemcpyDeviceToHost) !=cudaSuccess ||
        cudaMemcpy(oddArr, mOddSum, mSumBytes, cudaMemcpyDeviceToHost) != cudaSuccess) {
      pflerr("Failure to copy even/odd sums back from GPU");
      error = 1;
    }
    normalize(evenArr, fftScale, sumXplus * mSumYpad);
    normalize(oddArr, fftScale, sumXplus * mSumYpad);
  }
  ADD_TIME(mWallCopy);
  if (evenOddOnly)
    return error;

  // If there are even and odd, add them into even
  if (mDoEvenOdd) {
    if (bindSumArray(0)) {
      error += 2;
    } else {
      sOddTex.filterMode = cudaFilterModePoint;
      sOddTex.normalized = false;
      if (cudaBindTexture(NULL, sOddTex, mOddSum, sChanDesc, mSumBytes) != cudaSuccess) {
        error += 2;
      }
    }

    if (error < 2) {
      START_TIMER;
      dim3 blockSize(blockX, 8, 1);
      dim3 gridSize((sumXplus + blockSize.x - 1) / blockSize.x, 
                    (mSumYpad + blockSize.y - 1) / blockSize.y, 1);
      addOddToEvenSum<<<gridSize, blockSize>>>
        (mEvenSum, sumXplus, mSumYpad);
      if (testReportErr("to add odd and even sums"))
        error += 2;
      if (error < 2 && cudaDeviceSynchronize() != cudaSuccess) {
        pflerr("Error return from synchronizing after adding odd and even sums");
        error += 2;
      }
      ADD_TIME(mWallAddEO);
    }
    cudaUnbindTexture(sOddTex);
  }

  // Inverse FFT after destroying forward plan
  destroyPlan(sFullForwardPlan);
  err = cufftPlan2d(&sSumInversePlan, mSumYpad, mSumXpad, CUFFT_C2R);
  if (err != CUFFT_SUCCESS) {
    utilPrint("Failed to make plan for inverse sum FFT (size %d %d, error %d)\n",
              mSumXpad, mSumYpad, err);
    error += 2;
  }

  START_TIMER;
  if (error < 2) {
    err = cufftExecC2R(sSumInversePlan, (cufftComplex *)mEvenSum, mEvenSum);
    if (err != CUFFT_SUCCESS) {
      utilPrint("Failure in sum inverse FFT on GPU (error %d)\n", err);
      error += 2;
    }
  }
  if (error < 2 && mTrackTime && cudaDeviceSynchronize() != cudaSuccess) {
    pflerr("Error return from synchronizing after inverse FFT on GPU");
    return error + 2;
  }
  ADD_TIME(mWallFFT);
  destroyPlan(sSumInversePlan);

  // Get result back
  START_TIMER;
  if (error < 2 && cudaMemcpy(sumArr, mEvenSum, mSumBytes, cudaMemcpyDeviceToHost) != 
      cudaSuccess) {
    pflerr("Failure to copy real sum from GPU");
    error += 2;
  }
  ADD_TIME(mWallCopy);

  // If there is a failure in result and even only not copied yet, copy it
  if (error == 2 && !mDoEvenOdd) {
    if (cudaMemcpy(evenArr, mEvenSum, mSumBytes, cudaMemcpyDeviceToHost) != cudaSuccess){
      error++;
    } else if (mOddSum) {
      memset(mOddSum, 0, mSumBytes);
      normalize(evenArr, fftScale, sumXplus * mSumYpad);
    }
  }
  if (error < 2)
    normalize(sumArr, imScale, sumXplus * mSumYpad);

  return error;
}

/*
 * Inverse transform and return the non-dose weighted sum
 */
int FrameGPU::returnUnweightedSum(float *sumArr)
{
  float fftScale = 1. / (sqrt((double)mFullXpad * mFullYpad) * (float)mFullYpad/mSumYpad);
  float imScale = (fftScale / sqrt((double)mSumXpad * mSumYpad));
  int err;
  if (!mDoUnDWsum)
    return 1;
  err = cufftPlan2d(&sSumInversePlan, mSumYpad, mSumXpad, CUFFT_C2R);
  if (err != CUFFT_SUCCESS) {
    utilPrint("Failed to make plan for inverse nonDW FFT (size %d %d, error %d)\n",
              mSumXpad, mSumYpad, err);
    return 1;
  }

  START_TIMER;
  err = cufftExecC2R(sSumInversePlan, (cufftComplex *)mNonDWsum, mNonDWsum);
  if (err != CUFFT_SUCCESS) {
    utilPrint("Failure in nonDW sum inverse FFT on GPU (error %d)\n", err);
    return 1;
  }
  if (mTrackTime && cudaDeviceSynchronize() != cudaSuccess) {
    pflerr("Error return from synchronizing after inverse FFT on GPU");
    return 1;
  }
  ADD_TIME(mWallFFT);
  destroyPlan(sSumInversePlan);

  // Get result back. The supplied array must be big enough: (x + 2) * y
  START_TIMER;
  if (testErrCode(cudaMemcpy(sumArr, mNonDWsum, mSumBytes, cudaMemcpyDeviceToHost), 
                  "copy nonDW sum from GPU", 0))
    return 1;
  ADD_TIME(mWallCopy);
  normalize(sumArr, imScale, (mSumXpad + 2) * mSumYpad);
  return 0;
}

/*
 * Take in the next align image, copy to appropriate place, transform and filter
 * Errors in bin/pad and stacking operations return a 2 so that only initial operations
 * can be moved to CPU, otherwise return a 1
 */
int FrameGPU::processAlignImage(float *binArr, int stackInd, int groupInd, int stackOnGpu)
{  int err, type2d, dataSize, frmStkInd, blockX = 32;
  int ind, nxBox, nyBox, ixLow, ixHigh, iyLow, iyHigh, iyStart;
  float *devArr = mWorkBinPad;
  float *groupArr;
  float edgeMean;
  bool makingGroup = mGroupSize > 1 && groupInd >= 0;
  FloatVec weights;
  IntVec iWeights;
  int support, delXstart, delXend, delYstart, delYend;
  cudaArray *dev2dArr;

  if (manageRawTempArray(1))
    return 2;

  // Set array to copy into, make new if needed and push on stack
  if (stackInd >= 0) {
    if (stackInd < mSavedBinPad.size()) {
      devArr = mSavedBinPad[stackInd];
    } else if (stackInd == mSavedBinPad.size()) {
      if (testErrCode(cudaMalloc((void **)&devArr, mAlignBytes), 
                      "allocate new array for align image on GPU", 0))
        return 1;
      mSavedBinPad.push_back(devArr);
    } else {
      utilPrint("Stack index in call to processAlignImage is %d, stack size is only %d\n",
                stackInd, mSavedBinPad.size());
      return 1;
    }
  }

  // If grouping, set array for that too
  if (makingGroup) {
    if (groupInd > (int)mSavedGroups.size()) {
      utilPrint("Group index in call to processAlignImage is %d, group stack size is only"
                " %d\n", groupInd, mSavedGroups.size());
      return 1;
    }
    if (groupInd < mSavedGroups.size()) {
      groupArr = mSavedGroups[groupInd];
    } else {
      if (testErrCode(cudaMalloc((void **)&groupArr, mAlignBytes),
                      "allocate new array for group align image on GPU", 0))
        return 1;
      mSavedGroups.push_back(groupArr);
    }
  }

  /* Data flow for doing bin/pad on GPU:
         No stack                         Stack
     load to sTempRawArray       load to raw stack array

                       Need preproc
     proc to sTempFloatArray     proc to sTempFloatArray  
   */

  if (mDoAlignBinPad) {

    // Set where to load this image to
    dev2dArr = sTempRawArray;
    type2d = mStackType;
    if (stackOnGpu) {

      // When stacking, we need to get the edge mean now while we have the data
      edgeMean = frameEdgeMean(binArr, type2d, mUnpaddedX, FULL_INDENT, 
                               mUnpaddedX - FULL_INDENT - 1, mUnpaddedY - FULL_INDENT,
                               mUnpaddedY - FULL_INDENT - 1);

      // Put in existing spot at end of stack or make a new spot on stack
      frmStkInd = mNumOnUnpadStack;
      if (frmStkInd < sSavedUnpadded.size()) {
        dev2dArr = sSavedUnpadded[frmStkInd];
        mUnpadEdgeMeans[frmStkInd] = edgeMean;
        mSavedFrameNums[frmStkInd] = mNumAlignedFrames;
        mNumOnUnpadStack++;
      } else {
        if (testErrCode(malloc2Darray(&dev2dArr, mStackType, mUnpaddedX, mUnpaddedY),
                        "allocate 2D array for full stack on GPU", 0))
          return 2;
        sSavedUnpadded.push_back(dev2dArr);
        mUnpadEdgeMeans.push_back(edgeMean);
        mSavedFrameNums.push_back(mNumAlignedFrames);
        mNumOnUnpadStack++;
      }
    }

    // copy to array
    dataSizeForMode(mStackType, &dataSize, &err);
    START_TIMER;
    if (testErrCode(cudaMemcpyToArray(dev2dArr, 0, 0, binArr, mUnpaddedX * mUnpaddedY *
                                      dataSize, cudaMemcpyHostToDevice),
                    "copy unpadded image to GPU array for bin/pad", 0))
      return 2;
    ADD_TIME(mWallCopy);
    //dumpUnpadArray(dev2dArr, mUnpaddedX, mUnpaddedY, type2d, "input array");
    
    if (mDoPreprocess) {
      if (runPreprocess(dev2dArr, type2d, mNumAlignedFrames))
        return 2;
      
      dev2dArr = sTempFloatArray;
      type2d = MRC_MODE_FLOAT;
    }
    //dumpUnpadArray(dev2dArr, mUnpaddedX, mUnpaddedY, type2d, "preproc array");

    if (testErrCode(bindUnpadArray(dev2dArr, type2d), 
                    "bind input array to texture for reduction", 0))
      return 2;

    nxBox = (mXend + 1 - mXstart) / mAliBinning;
    nyBox = (mYend + 1 - mYstart) / mAliBinning;
    ixLow = mAlignXpad / 2 - nxBox / 2;
    ixHigh = ixLow + nxBox;
    iyLow = mAlignYpad / 2 - nyBox / 2;
    iyHigh = iyLow + nyBox;
    
    edgeMean = frameEdgeMean(binArr, mStackType, mUnpaddedX, mXstart, mXend, mYstart,
                             mYend);

    if (mAliBinning > 1) {

      // Set up the reduction filter for X and copy to constant memory
      selectZoomFilter(mAntiFiltType, 1. / mAliBinning, &support);
      makeReductionWeights(mXstart, support, weights, delXstart, delXend);
      if (type2d == MRC_MODE_FLOAT) {
        if (testErrCode(cudaMemcpyToSymbol(trigTable, &weights[0], weights.size() * 
                                           sizeof(float), 0, cudaMemcpyHostToDevice),
                        "copy constant X weight data to GPU", 0))
          return 2;
      } else {
        for (ind = 0; ind < (int)weights.size(); ind++)
          iWeights.push_back(B3DNINT(INT_KERNEL_SCALE * weights[ind]));
        if (testErrCode(cudaMemcpyToSymbol(intKernel, &iWeights[0], iWeights.size() * 
                                           sizeof(int), 0, cudaMemcpyHostToDevice),
                        "copy integer constant weight data to GPU", 0))
          return 2;
      }
      iyStart = mYstart - (support / 2 + 1);

#ifndef NO_SURFACES
      if (testErrCode(cudaBindSurfaceToArray(sTempSurfRef, sReducedInXarray),
                      "bind X-reduced array to surface for reduction", 0))
        return 2;
#endif

      START_TIMER;
      dim3 blockSize(blockX, 8, 1);
      dim3 gridSize((mRedColX + blockSize.x - 1) / blockSize.x, 
                    (mRedColY + blockSize.y - 1) / blockSize.y, 1);
      reduceColumns<<<gridSize, blockSize>>>
        (mAliBinning, mReducedInXlinear, type2d, mRedColX, mRedColY, iyStart, delXstart,
         delXend);
      if (testReportErr("to do reduction in X"))
        return 2;

      // Do not synchronize here  - saves 3% of time
      /*if (cudaDeviceSynchronize() != cudaSuccess) {
        pflerr("Error return from synchronizing after reducing/padding image");
        return 2;
        } */ 
      unbindUnpadArray(type2d);
      //dumpUnpadArray(sReducedInXarray, mRedColX, mRedColY, MRC_MODE_FLOAT, 
      //"x-reduced array");

      if (mNoSurfaces && testErrCode
          (cudaMemcpyToArray(sReducedInXarray, 0, 0, mReducedInXlinear, mRedColX * 
                             mRedColY * sizeof(float), cudaMemcpyDeviceToDevice),
           "copy X-reduced output array to 2D array", 0))
        return 1;
    
      if (testErrCode(bindUnpadArray(sReducedInXarray, MRC_MODE_FLOAT), 
                      "bind X-reduced array to texture for Y reduction", 0))
        return 2;

      weights.clear();
      makeReductionWeights(mYstart, support, weights, delYstart, delYend);
      if (testErrCode(cudaMemcpyToSymbol(trigTable, &weights[0], weights.size() * 
                                         sizeof(float), 0, cudaMemcpyHostToDevice),
                      "copy constant Y weight data to GPU", 0))
        return 2;

      //dim3 blockSize(blockX, 8, 1);
      dim3 gridSize2((mUnpaddedX + blockSize.x - 1) / blockSize.x, 
                     (mUnpaddedY + blockSize.y - 1) / blockSize.y, 1);
      reduceRowsTaperPad<<<gridSize2, blockSize>>>
        (devArr, mAliBinning, iyStart, ixLow, ixHigh, iyLow, 
         iyHigh, mAlignXpad + 2, mAlignXpad, mAlignYpad, mNxTaper, mNyTaper, edgeMean, 
         delYstart, delYend);
      if (testReportErr("to reduce, taper and pad image"))
        return 2;
    } else {

      // binning 1 is a simple copy of trimmed area with padding and taper inside
      START_TIMER;
      dim3 blockSize(blockX, 8, 1);
      dim3 gridSize((mUnpaddedX + blockSize.x - 1) / blockSize.x, 
                    (mUnpaddedY + blockSize.y - 1) / blockSize.y, 1);
      trimTaperPad<<<gridSize, blockSize>>>
        (devArr, type2d, ixLow, ixHigh, iyLow, iyHigh, mAlignXpad + 2, mAlignXpad,
         mAlignYpad, mNxTaper, mNyTaper, edgeMean);
      if (testReportErr("to trim, taper and pad image"))
        return 2;
    }
    
    if (cudaDeviceSynchronize() != cudaSuccess) {
      pflerr("Error return from synchronizing after reducing/padding image");
      return 2;
    }
    ADD_TIME(mWallRedPad);
    unbindUnpadArray(type2d);

  } else {

    // NO BIN/PAD HERE, copy taper-padded image to device
    START_TIMER;
    if (testErrCode(cudaMemcpy(devArr, binArr, mAlignBytes, cudaMemcpyHostToDevice),
                    "copy align image to GPU array", 0))
      return 1;
    ADD_TIME(mWallCopy);
  }
  //dumpImage(devArr, mAlignXpad + 2, mAlignXpad, mAlignYpad, 0, "reduction", mNumAlignedFrames);
  
  // If grouping and not refining, just sum real-space stack into group array
  if (makingGroup && !mDoAlignSum) {
    if (sumIntoGroup(stackInd, groupInd))
      return 1;
    devArr = groupArr;
  }
  
  if (mGroupSize == 1 || mDoAlignSum || groupInd >= 0) {
    
    // Take the FFT
    START_TIMER;
    err = cufftExecR2C(sAlignForwardPlan, devArr, (cufftComplex *)devArr);
    if (err != CUFFT_SUCCESS) {
      utilPrint("Failure in forward FFT of align image on GPU (error %d)\n", err);
      return 1;
    }
    if (mTrackTime && cudaDeviceSynchronize() != cudaSuccess) {
      pflerr("Error return from synchronizing after forward FFT of align image on GPU");
      return 1;
    }
    ADD_TIME(mWallFFT);
    
    // Filter the array
    START_TIMER;
    dim3 blockSize(blockX, 8, 1);
    dim3 gridSize((mAlignXpad + 2 + blockSize.x - 1) / blockSize.x, 
                  (mAlignYpad + blockSize.y - 1) / blockSize.y, 1);
    filterAlignFFT<<<gridSize, blockSize>>>
      (devArr, mAlignXpad + 2, mAlignYpad);
    if (testReportErr("to filter align image FFT"))
      return 1;
    
    if (cudaDeviceSynchronize() != cudaSuccess) {
      pflerr("Error return from synchronizing after filtering align image FFT");
      return 1;
    }
    ADD_TIME(mWallFilt);
  }

  // Or if grouping and refining, now sum the filtered FFTs
  if (makingGroup && mDoAlignSum && sumIntoGroup(stackInd, groupInd))
    return 1;
  
  mNumAlignedFrames++;
  return 0;
}

/* 
 * Do the preprocessing: gain normalization, truncation, defect removal
 */
 int FrameGPU::runPreprocess(void *dev2dArr, int type2d, int frame)
{
  int blockX = 16;

  //dumpUnpadArray(dev2dArr, mUnpaddedX, mUnpaddedY, type2d, "raw image");
  
  // Bind input array to texture and output array to surface
  if (testErrCode(bindUnpadArray((cudaArray *)dev2dArr, type2d), 
                  "bind array to texture for preprocessing", 0))
    return 1;
#ifndef NO_SURFACES
  if (testErrCode(cudaBindSurfaceToArray(sTempSurfRef, sTempFloatArray),
                  "bind output array to surface for preprocessing", 0))
    return 1;
#endif
  
  START_TIMER;
  dim3 blockSize(blockX, 8, 1);
  dim3 gridSize((mUnpaddedX + 2 + blockSize.x - 1) / blockSize.x, 
                (mUnpaddedY + blockSize.y - 1) / blockSize.y, 1);
  preprocessFrame<<<gridSize, blockSize>>>
    (mProcessedLinear, type2d, mUnpaddedX, mUnpaddedY, mDoGainNorm, mTruncLimit,
        mCorrectDefects);
  if (testReportErr("to preprocess image"))
    return 1;
      
  if (cudaDeviceSynchronize() != cudaSuccess) {
    pflerr("Error return from synchronizing after preprocessing image");
    return 1;
  }
  if (mNoSurfaces && testErrCode
      (cudaMemcpyToArray(sTempFloatArray, 0, 0, mProcessedLinear, mUnpaddedX *
                         mUnpaddedY * sizeof(float), cudaMemcpyDeviceToDevice),
       "copy pre-processed output array to 2D array", 0))
    return 1;
  ADD_TIME(mWallPreproc);

  // Unbind input texture: there is no surface unbinding as the next bind unbinds it
  unbindUnpadArray(type2d);
  //dumpUnpadArray(sTempFloatArray, mUnpaddedX, mUnpaddedY, MRC_MODE_FLOAT,
  //             "gain normed", frame);
  return 0;
}

/*
 * Get weights for image reduction in one dimension and offsets that apply to given 
 * starting coordinate in the image
 */
void FrameGPU::makeReductionWeights(int startCoord, int support, FloatVec &weights,
                                    int &delStart, int &delEnd)
{
  // Replicate variables and equations in zoomFiltInterp routine when passed startCoord
  float cen = startCoord + 0.5 * mAliBinning;
  int ind, i0 = floor(cen + 0.5 - support);
  int i1 = ceil(cen + 0.5 + support);
  float filtVal;
  delStart = -999;

  // Loop on support and keep track of first and last non-zero value
  // The function returns normalized values
  for (ind = i0; ind <= i1; ind++) {
    filtVal = zoomFiltValue(ind + 0.5 - cen);
    if (filtVal) {
      weights.push_back(filtVal);
      delEnd = ind;
      if (delStart == -999)
        delStart = ind;
    }
  }
}

/* 
 * Sum from 2 to 5 images or FFTs into a group array
 */
int FrameGPU::sumIntoGroup(int stackInd, int groupInd)
{
  int blockX = 32;
  float *groupArr;
  if (stackInd + 1 - mGroupSize < 0 || stackInd >= (int)mSavedBinPad.size() || 
      groupInd < 0 || groupInd >= (int)mSavedGroups.size()) {
    utilPrint("Index in call to sumIntoGroup is out of range (%d and %d - sizes %d and "
              "%d\n", stackInd, groupInd, mSavedBinPad.size(), mSavedGroups.size());
    return 1;
  }
  groupArr = mSavedGroups[groupInd];

  START_TIMER;
  dim3 blockSize(blockX, 8, 1);
  dim3 gridSize((mAlignXpad + 2 + blockSize.x - 1) / blockSize.x, 
                (mAlignYpad + blockSize.y - 1) / blockSize.y, 1);
  if (mGroupSize == 2)
    sum2IntoGroup<<<gridSize, blockSize>>>
      (mSavedBinPad[stackInd - 1], mSavedBinPad[stackInd], 
       groupArr, mAlignXpad + 2, mAlignYpad);
  else if (mGroupSize == 3)
    sum3IntoGroup<<<gridSize, blockSize>>>
      (mSavedBinPad[stackInd - 2], mSavedBinPad[stackInd - 1], mSavedBinPad[stackInd], 
       groupArr, mAlignXpad + 2, mAlignYpad);
  else if (mGroupSize == 4)
    sum4IntoGroup<<<gridSize, blockSize>>>
      (mSavedBinPad[stackInd - 3], mSavedBinPad[stackInd - 2], mSavedBinPad[stackInd - 1],
       mSavedBinPad[stackInd], groupArr, mAlignXpad + 2, mAlignYpad);
  else
    sum5IntoGroup<<<gridSize, blockSize>>>
      (mSavedBinPad[stackInd - 4], mSavedBinPad[stackInd - 3], mSavedBinPad[stackInd - 2],
       mSavedBinPad[stackInd - 1], mSavedBinPad[stackInd], groupArr, mAlignXpad + 2,
       mAlignYpad);
  if (testReportErr("to sum group on GPU"))
    return 1;
  
  if (cudaDeviceSynchronize() != cudaSuccess) {
    pflerr("Error return from synchronizing after summing group on GPU");
    return 1;
  }
  ADD_TIME(mWallGroup);
  return 0;
}

/*
 * Correlate two images that have been stored somewhere
 * aliInd is a stack index or -1 to align mWorkBinPad; refInd is a stack index or -1 to
 * align to mAlignSum or -2 to align to mWorkBinPad
 */
int FrameGPU::crossCorrelate(int aliInd, int refInd, float *subarea, int subXoffset,
                             int subYoffset)
{
  float *aliArr = mWorkBinPad;
  float *refArr = mAlignSum;
  float *corrArr = mCorrBinPad;
  int ixFrom0[4], ixTo0[4], iyFrom0[4], iyTo0[4], ixFrom1[4], ixTo1[4], iyFrom1[4];
  int iyTo1[4], roundXoffset, roundYoffset, ix0, iy0, bigSubBytes, err;
  int blockX = 32;
  std::vector<float *> *savedArr = mGroupSize > 1 ? &mSavedGroups : &mSavedBinPad;

  // Check the indexes and set the source and destination arrays 
  if (aliInd >= (int)savedArr->size() || refInd >= (int)savedArr->size()) {
    utilPrint("Index for align image (%d) or reference (%d) is too big (stack size %d)",
              aliInd, refInd, savedArr->size());
    return 1;
  }
  if (aliInd >= 0)
    aliArr = (*savedArr)[aliInd];
  if (refInd >= 0)
    refArr = (*savedArr)[refInd];
  if (refInd < -1)
    refArr = mWorkBinPad;
  if (aliInd >= 0 && refInd >= 0)
    corrArr = mWorkBinPad;

  // Take the conjugate product
  START_TIMER;
  dim3 blockSize(blockX, 8, 1);
  dim3 gridSize((mAlignXpad + 2 + blockSize.x - 1) / blockSize.x, 
                (mAlignYpad + blockSize.y - 1) / blockSize.y, 1);
  conjugateProduct<<<gridSize, blockSize>>>
    (aliArr, refArr, corrArr, mAlignXpad / 2 + 1, mAlignYpad);
  if (testReportErr("to take conjugate product"))
    return 1;

  if (cudaDeviceSynchronize() != cudaSuccess) {
    pflerr("Error return from synchronizing after taking conjugate product");
    return 1;
  }
  ADD_TIME(mWallConj);

  START_TIMER;
  err = cufftExecC2R(sAlignInversePlan, (cufftComplex *)corrArr, mRealCorr);
  if (err != CUFFT_SUCCESS) {
    utilPrint("Failure in inverse FFT of align image on GPU (error %d)\n", err);
    return 1;
  }
  if (mTrackTime && cudaDeviceSynchronize() != cudaSuccess) {
    pflerr("Error return from synchronizing after inverse FFT of align image on GPU");
    return 1;
  }
  ADD_TIME(mWallFFT);

  // extract the subarea
  // This is an intrinsically slow operation, but copying lines direct from device
  // to host was almost twice as slow
  roundXoffset = NICE_GPU_DIVISOR * B3DNINT((float)subXoffset / NICE_GPU_DIVISOR);
  roundYoffset = NICE_GPU_DIVISOR * B3DNINT((float)subYoffset / NICE_GPU_DIVISOR);
  utilCoordsForWrap(mAlignXpad, mAlignYpad, mBigSubareaSize, mBigSubareaSize,
                    roundXoffset, roundYoffset, ixFrom0, ixTo0, iyFrom0, iyTo0,
                    ixFrom1, ixTo1, iyFrom1, iyTo1);
  
  if (bindFullOrCorrArray(mRealCorr, mAlignBytes))
    return 1;
  
  dim3 gridSize2((mBigSubareaSize + blockSize.x - 1) / blockSize.x, 
                 (mBigSubareaSize + blockSize.y - 1) / blockSize.y, 1);
  wrapCorners<<<gridSize2, blockSize>>>
    (corrArr, mSubareaCorr, mAlignXpad, mAlignYpad, mBigSubareaSize, mBigSubareaSize, 
     ixFrom0[2], iyFrom0[2]);
  if (testReportErr("to wrap subarea corr in GPU"))
    return 1;
  
  if (cudaDeviceSynchronize() != cudaSuccess) {
    pflerr("Error return from synchronizing after wrapping subarea corr in GPU");
    return 1;
  }
  ADD_TIME(mWallExtract);
  //dumpImage(mSubareaCorr, mBigSubareaSize, mBigSubareaSize, mBigSubareaSize, 0,
  //"big subarea dev");
  
  // Unload the subarea
  mWallStart= wallTime();
  bigSubBytes = mBigSubareaSize * mBigSubareaSize * sizeof(float);
  if (cudaMemcpy(mHostSubarea, mSubareaCorr, bigSubBytes, cudaMemcpyDeviceToHost) != 
      cudaSuccess) {
    pflerr("Copying subarea correlation to host");
    return 1;
  }
  
  // Extract proper area into subarea array
  ix0 = (mBigSubareaSize - mAliFiltSize) / 2 + (subXoffset - roundXoffset);
  iy0 = (mBigSubareaSize - mAliFiltSize) / 2 + (subYoffset - roundYoffset);
  sliceTaperInPad(mHostSubarea, SLICE_MODE_FLOAT, mBigSubareaSize, ix0, 
                  ix0 + mAliFiltSize - 1, iy0, iy0 + mAliFiltSize - 1, subarea,
                  mAliFiltSize + 2, mAliFiltSize, mAliFiltSize, 0, 0);
  ADD_TIME(mWallCopy);
  return 0;
}

/*
 * Shift one alignment image and add it into align sum, possibly in place
 */
int FrameGPU::shiftAddToAlignSum(int stackInd, float shiftX, float shiftY,
                                 int shiftSource)
{
  float *binArr = mWorkBinPad;
  if (stackInd >= (int)mSavedBinPad.size()) {
    utilPrint("Stack index in call to shiftAddToalignSum is %d, stack size is "
              "only %d\n", stackInd, mSavedBinPad.size());
    return 1;
  }
  if (stackInd >= 0)
    binArr = mSavedBinPad[stackInd];
  //dumpFFT(binArr, mAlignXpad, mAlignYpad, "fft of shift/add", 1);
  int err = shiftAddCommon(binArr, mAlignSum, 2, mAlignXpad, mAlignYpad, mAlignXpad, 
                           mAlignYpad, shiftX, shiftY, shiftSource, false);
  //dumpFFT(binArr, mAlignXpad, mAlignYpad, "fft after shifing", 1);
  //dumpFFT(mAlignSum, mAlignXpad, mAlignYpad, "Align sum", 1);
  return err;
}

/*
 * Upload a new filter mask to the the device for high-frequency filtering in refinement
 */
int FrameGPU::newFilterMask(float *alignMask)
{
  START_TIMER;
  if (testErrCode(cudaMemcpy(mFiltMask, alignMask, mAlignBytes, cudaMemcpyHostToDevice),
                  "copy filter mask to GPU array", 0))
    return 1;
  ADD_TIME(mWallCopy);
  return 0;
}

/*
 * Subtract one image on stack from the align sum and apply the filter mask to the
 * leave-one-out sum
 */
int FrameGPU::subtractAndFilterAlignSum(int stackInd, int groupRefine)
{
  int blockX = 32;
  std::vector<float *> *savedArr = groupRefine ? &mSavedGroups : &mSavedBinPad;

  if (stackInd >= (int)savedArr->size() || stackInd < 0) {
    utilPrint("Stack index in call to addOrSubtractAlignSum is %d, stack size is "
              "only %d\n", stackInd, savedArr->size());
    return 1;
  }

  START_TIMER;
  bindSumArray(2);
  dim3 blockSize(blockX, 8, 1);
  dim3 gridSize((mAlignXpad + 2 + blockSize.x - 1) / blockSize.x, 
                (mAlignYpad + blockSize.y - 1) / blockSize.y, 1);
  subtractFilterSum<<<gridSize, blockSize>>>
    ((*savedArr)[stackInd], mWorkBinPad, mAlignXpad + 2, mAlignYpad);
  if (testReportErr("to subtract from aligned sum"))
    return 1;

  if (cudaDeviceSynchronize() != cudaSuccess) {
    pflerr("Error return from synchronizing after subtracting from align sum");
    return 1;
  }
  ADD_TIME(mWallSubtract);
  return 0;
}

/*
 * Roll the stack of alignment images
 */
void FrameGPU::rollAlignStack()
{
  if (mSavedBinPad.size() == 0)
    return;
  utilRollSavedFrames(mSavedBinPad, mSavedBinPad.size());
}

/*
 * Roll the stack of group images
 */
void FrameGPU::rollGroupStack()
{
  if (mGroupSize > 1 && mSavedGroups.size() > 0)
    utilRollSavedFrames(mSavedGroups, mSavedGroups.size());
}

/*
 * Return the number of arrays saved in savedBinPad and savedGroups so that allocations
 * can be made for fetching them back after an error
 */
void FrameGPU::numberOfAlignFFTs(int *numBinPad, int *numGroups)
{
  *numBinPad = (int)mSavedBinPad.size();
  *numGroups = (int)mSavedGroups.size();
}

/*
 * Return all the data from GPU being used for alignment to try to recover from error
 */
int FrameGPU::returnAlignFFTs(float **saved, float **groups, float *alignSum,
                              float *workArr)
{
  int ind, loop, alignPix = (mAlignXpad + 2) * mAlignYpad;
  float *array;
  float fftScale = 1. / sqrt((double)mAlignXpad * mAlignYpad);
  std::vector<float *> *devSaved = &mSavedBinPad;
  float **hostSaved = saved;
  for (loop = 0; loop < 2; loop ++) {
    for (ind = 0; ind < devSaved->size(); ind++) {
      if (!hostSaved)
        return 1;
      array = hostSaved[ind];
      if (cudaMemcpy(array, devSaved->at(ind), mAlignBytes, cudaMemcpyDeviceToHost) !=
          cudaSuccess) {
        pflerr("Error copying align FFT back from GPU");
        return 1;
      }

      // Normalize unless it is a real-space image kept to make group sum in mSavedBinPad
      if (!(loop == 0 && mGroupSize > 1 && !mDoAlignSum))
        normalize(array, fftScale, alignPix);
    }
    devSaved = &mSavedGroups;
    hostSaved = groups;
  }
  if (alignSum) {
    if (cudaMemcpy(alignSum, mAlignSum, mAlignBytes, cudaMemcpyDeviceToHost) != 
        cudaSuccess) {
      pflerr("Error copying align sum FFT back from GPU");
      return 1;
    }
    normalize(array, fftScale, alignPix);
  }
  if (workArr) {
    if (cudaMemcpy(workArr, mWorkBinPad, mAlignBytes, cudaMemcpyDeviceToHost) != 
        cudaSuccess) {
      pflerr("Error copying align work FFT back from GPU");
      return 1;
    }
    normalize(workArr, fftScale, alignPix);
  }
  return 0;
}

/*
 * Return a stacked full frame from the first position on that stack, and its frame number
 */
int FrameGPU::returnStackedFrame(float *array, int *frameNum)
{
  cudaArray *dev2dArr;
  int err, dataSize;
  dataSizeForMode(mStackType, &dataSize, &err);
  if (!mNumOnUnpadStack) {
    utilPrint("Program error: returnStackedFrame called with no frames left on stack\n");
    return 1;
  }
  dev2dArr = sSavedUnpadded[0];
  *frameNum = mSavedFrameNums[0];
  if (cudaMemcpyFromArray(array, dev2dArr, 0, 0, dataSize * mUnpaddedX * mUnpaddedY,
                      cudaMemcpyDeviceToHost) != cudaSuccess) {
    pflerr("Error copying stacked full frame back from GPU");
    return 1;
  }
  free2Darray(&dev2dArr);

  // Remove from stack
  sSavedUnpadded.erase(sSavedUnpadded.begin());
  mSavedFrameNums.erase(mSavedFrameNums.begin());
  mUnpadEdgeMeans.erase(mUnpadEdgeMeans.begin());
  mNumOnUnpadStack--;
  return 0;
}

/*
 * Bind whichever sum array is needed
 */
int FrameGPU::bindSumArray(int needBound) 
{
  float *array = mEvenSum;
  int sizeArr = mSumBytes;
  if (needBound == mBoundSum)
    return 0;
  if (needBound == 1) {
    array = mOddSum;
  } else if (needBound == 2) {
    array = mAlignSum;
    sizeArr = mAlignBytes;
  }

  if (mBoundSum >= 0)
    cudaUnbindTexture(sSumTex);
  if (testErrCode(cudaBindTexture(NULL, sSumTex, array, sChanDesc, sizeArr),
                  "bind array for sum to texture", 0))
    return 1;
  mBoundSum = needBound;
  return 0;
}

/*
 * And bind whichever full or correlation array is needed to this texture
 */
int FrameGPU::bindFullOrCorrArray(float *fullArr, size_t sizeTmp)
{
  if (mBoundToFull != fullArr) {
    sFullTex.filterMode = cudaFilterModePoint;
    sFullTex.normalized = false;
    if (mBoundToFull)
      cudaUnbindTexture(sFullTex);
    mBoundToFull = NULL;
    if (testErrCode(cudaBindTexture(NULL, sFullTex, fullArr, sChanDesc, sizeTmp), 
                    "bind source array for summing to texture", 0))
      return 1;
    mBoundToFull = fullArr;
  }
  return 0;
}

/*
 * Normalize an array after FFT
 */
void FrameGPU::normalize(float *data, float scale, int numPix)
{
  for (int ind = 0; ind < numPix; ind++)
    data[ind] *= scale;
}

/*
 * Make sure the arrays for storing shift sines and cosines for copy to constant
 * memory are big enough
 */
int FrameGPU::manageShiftTrigs(int xpad, int ypad)
{
  if (xpad + 2 > mXtrigSize) {
    B3DFREE(mXshiftTrig);
    mXshiftTrig = B3DMALLOC(float, xpad + 2);
    mXtrigSize = xpad + 2;
  }

  if (ypad * 2 > mYtrigSize) {
    B3DFREE(mYshiftTrig);
    mYshiftTrig = B3DMALLOC(float, ypad * 2);
    mYtrigSize = ypad * 2;
  }
  if (!mXshiftTrig || !mYshiftTrig) {
    utilPrint("Failed to allocate arrays for sines/cosines\n");
    return 1;
  }
  return 0;
}

void FrameGPU::printTimers()
{
  utilPrint("GPU: copy %.4f  pre %.4f  FFT %.4f  shift %.4f  e/o add %.4f  filt %.4f\n "
            "   red %.4f  noise %.4f  conj %.4f  extr %.4f  subtr %.4f "
            " group %.4f\n", mWallCopy, mWallPreproc,
            mWallFFT, mWallShift, mWallAddEO, mWallFilt,  mWallRedPad,
            mWallNoise, mWallConj, mWallExtract, mWallSubtract, mWallGroup);
}

// Test for and report error after executing threads           
int FrameGPU::testReportErr(const char *mess)
{
  cudaError_t err;
  err = cudaGetLastError();
  if (err != cudaSuccess) {
    utilPrint("Error executing threads %s: %s\n", mess,
              cudaGetErrorString(err));
    return 1;
  }
  return 0;
}

// Test error code from cuda calls and output messages
int FrameGPU::testErrCode(int errCode, const char *mess, int cleanAll)
{
  if (errCode != cudaSuccess) {
    pflerr("Failed to %s", mess);
    if (cleanAll)
      cleanup();
    return 1;
  }
  return 0;
}

// In case of error, find the error string and print it with message
void FrameGPU::pflerr(const char *format, ...)
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
 * Routines to output images afer unloading from GPU
 */
void FrameGPU::dumpFFT(float *fft, int nxPad, int nyPad, const char *descrip, int doReal)
{
  int sizeTmp = (nxPad + 2) * nyPad;
  float *temp = B3DMALLOC(float, sizeTmp);
  if (!temp)
    return;
  if (cudaMemcpy(temp, fft, sizeTmp * sizeof(float), cudaMemcpyDeviceToHost) ==
      cudaSuccess)
    utilDumpFFT(temp, nxPad, nyPad, descrip, doReal);
  free(temp);
}

void FrameGPU::dumpImage(float *image, int nxDim, int nxPad, int nyPad, int isCorr,
                         const char *descrip, int frame)
{
  int sizeTmp = nxDim * nyPad;
  float *temp = B3DMALLOC(float, sizeTmp);
  if (!temp)
    return;
  if (cudaMemcpy(temp, image, sizeTmp * sizeof(float), cudaMemcpyDeviceToHost) == 
      cudaSuccess)
    utilDumpImage(temp, nxDim, nxPad, nyPad, isCorr, descrip);
  free(temp);
}

void FrameGPU::dumpUnpadArray(void *array, int sizeX, int sizeY, int type,
                              const char *descrip, int frame)
{
  int err, dataSize, sizeTmp;
  cudaArray *dev2dArr = (cudaArray *)array;
  dataSizeForMode(type, &dataSize, &err);
  sizeTmp = sizeX * sizeY * dataSize;
  unsigned char *temp = B3DMALLOC(unsigned char, sizeTmp);
  if (!temp)
    return;
  float *ftemp = B3DMALLOC(float, sizeTmp);
  if (!ftemp) {
    free(temp);
    return;
  }
  if (cudaMemcpyFromArray(temp, dev2dArr, 0, 0, sizeTmp, cudaMemcpyDeviceToHost) == 
      cudaSuccess) {
    sliceTaperOutPad(temp, type, sizeX, sizeY, ftemp, sizeX,sizeX, sizeY , 0, 0.);
    utilDumpImage(ftemp, sizeX, sizeX, sizeY, 0, descrip, frame);
  }

  free(temp);
  free(ftemp);
}

/*
 * Possible addition to libcfshr/filtxcorr, takes an edge mean for any data type
 */
float FrameGPU::frameEdgeMean(void *array, int type, int nxdim, int ixlo, int ixhi,
                              int iylo, int iyhi)
{
  double sum = 0.;
  float dmean;
  float *fdata = (float *)array;
  unsigned char *bdata = (unsigned char *)array;
  short *sdata = (short *)array;
  unsigned short *usdata = (unsigned short *)array;
  int ix, iy;

  switch (type) {
  case MRC_MODE_BYTE:
    for (ix = ixlo; ix <= ixhi; ix++) {
      sum += bdata[ix + iylo * nxdim];
      sum += bdata[ix + iyhi * nxdim];
    }
    for (iy = iylo + 1; iy < iyhi; iy++) {
      sum += bdata[ixlo + iy * nxdim];
      sum += bdata[ixhi + iy * nxdim];
    }
    break;
  case MRC_MODE_SHORT:
    for (ix = ixlo; ix <= ixhi; ix++)
      sum += sdata[ix + iylo * nxdim] + sdata[ix + iyhi * nxdim];
    for (iy = iylo + 1; iy < iyhi; iy++)
      sum += sdata[ixlo + iy * nxdim] + sdata[ixhi + iy * nxdim];
    break;
  case MRC_MODE_USHORT:
    for (ix = ixlo; ix <= ixhi; ix++)
      sum += usdata[ix + iylo * nxdim] + usdata[ix + iyhi * nxdim];
    for (iy = iylo + 1; iy < iyhi; iy++)
      sum += usdata[ixlo + iy * nxdim] + usdata[ixhi + iy * nxdim];
    break;
  case MRC_MODE_FLOAT:
    for (ix = ixlo; ix <= ixhi; ix++)
      sum += fdata[ix + iylo * nxdim] + fdata[ix + iyhi * nxdim];
    for (iy = iylo + 1; iy < iyhi; iy++)
      sum += fdata[ixlo + iy * nxdim] + fdata[ixhi + iy * nxdim];
    break;
  }
  dmean = sum / (2 * (ixhi - ixlo + iyhi - iylo));
  return dmean;
}


////////////////////////////////////////////////////////////////////////////
// WRAPPER FUNCTIONS
///////////////////////////////////////////////////////////////////////////

DLL_EX_IM int fgpuGpuAvailable(int nGPU, float *memory, int debug)
{
  return sFGPU.gpuAvailable(nGPU, memory, debug);
}

DLL_EX_IM void fgpuSetUnpaddedSize(int unpadX, int unpadY, int flags, int debug)
{
  sFGPU.setUnpaddedSize(unpadX, unpadY, flags, debug);
}

DLL_EX_IM int fgpuSetPreProcParams(float *gainRef, int nxGain, int nyGain,
                                   float truncLimit, unsigned char *defectMap,
                                   int camSizeX, int camSizeY)
{
  return sFGPU.setPreProcParams(gainRef, nxGain, nyGain, truncLimit, defectMap,
                               camSizeX, camSizeY);
}

DLL_EX_IM void fgpuSetBinPadParams(int xstart, int xend, int ystart, int yend,
                                   int binning, int nxTaper, int nyTaper, int type,
                                   int filtType, int noiseLen)
{
  sFGPU.setBinPadParams(xstart, xend, ystart, yend, binning, nxTaper, nyTaper,
                        type, filtType, noiseLen);
}

DLL_EX_IM int fgpuSetupSumming(int fullXpad, int fullYpad, int sumXpad, int sumYpad,
                               int evenOdd)
{
  return sFGPU.setupSumming(fullXpad, fullYpad, sumXpad, sumYpad, evenOdd);
}

DLL_EX_IM int fgpuSetupAligning(int alignXpad, int alignYpad, int sumXpad, int sumYpad,
                                float *alignMask, int aliFiltSize, int groupSize,
                                int expectStackSize, int doAlignSum)
{
  return sFGPU.setupAligning(alignXpad, alignYpad, sumXpad, sumYpad, alignMask,
                             aliFiltSize, groupSize, expectStackSize, doAlignSum);
}

DLL_EX_IM int fgpuSetupDoseWeighting(float *filter, int filtSize, float delta)
{
  return sFGPU.setupDoseWeighting(filter, filtSize, delta);
}

DLL_EX_IM int fgpuAddToFullSum(float *fullArr, float shiftX, float shiftY)
{
  return sFGPU.addToFullSum(fullArr, shiftX, shiftY);
}

DLL_EX_IM int fgpuReturnSums(float *sumArr, float *evenArr, float *oddArr,
                             int evenOddOnly)
{
  return sFGPU.returnSums(sumArr, evenArr, oddArr, evenOddOnly);
}

DLL_EX_IM int fgpuReturnUnweightedSum(float *sumArr)
{
  return sFGPU.returnUnweightedSum(sumArr);
}

DLL_EX_IM void fgpuCleanup()
{
  sFGPU.cleanup();
}

DLL_EX_IM void fgpuRollAlignStack()
{
  sFGPU.rollAlignStack();
}

DLL_EX_IM void fgpuRollGroupStack()
{
  sFGPU.rollGroupStack();
}

DLL_EX_IM int fgpuSubtractAndFilterAlignSum(int stackInd, int groupRefine)
{
  return sFGPU.subtractAndFilterAlignSum(stackInd, groupRefine);
}

DLL_EX_IM int fgpuNewFilterMask(float *alignMask)
{
  return sFGPU.newFilterMask(alignMask);
}

DLL_EX_IM int fgpuShiftAddToAlignSum(int stackInd, float shiftX, float shiftY,
                                     int shiftSource)
{
  return sFGPU.shiftAddToAlignSum(stackInd, shiftX, shiftY, shiftSource);
}

DLL_EX_IM int fgpuCrossCorrelate(int aliInd, int refInd, float *subarea, int subXoffset,
                                 int subYoffset)
{
  return sFGPU.crossCorrelate(aliInd, refInd, subarea, subXoffset, subYoffset);
}

DLL_EX_IM int fgpuProcessAlignImage(float *binArr, int stackInd, int groupInd, 
                                    int stackOnGpu)
{
  return sFGPU.processAlignImage(binArr, stackInd, groupInd, stackOnGpu);
}

DLL_EX_IM void fgpuNumberOfAlignFFTs(int *numBinPad, int *numGroups)
{
  sFGPU.numberOfAlignFFTs(numBinPad, numGroups);
}

DLL_EX_IM int fgpuReturnAlignFFTs(float **saved, float **groups, float *alignSum, 
                                  float *workArr)
{
  return sFGPU.returnAlignFFTs(saved, groups, alignSum, workArr);
}

DLL_EX_IM int fgpuReturnStackedFrame(float *array, int *frameNum)
{
  return sFGPU.returnStackedFrame(array, frameNum);
}

DLL_EX_IM void fgpuCleanSumItems()
{
  sFGPU.cleanSumItems();
}

DLL_EX_IM void fgpuCleanAlignItems()
{
  sFGPU.cleanAlignItems();
}

DLL_EX_IM void fgpuZeroTimers()
{
  sFGPU.zeroTimers();
}

DLL_EX_IM void fgpuPrintTimers()
{
  sFGPU.printTimers();
}

DLL_EX_IM int fgpuClearAlignSum()
{
  return sFGPU.clearAlignSum();
}

DLL_EX_IM int fgpuSumIntoGroup(int stackInd, int groupInd)
{
  return sFGPU.sumIntoGroup(stackInd, groupInd);
}

DLL_EX_IM void fgpuSetGroupSize(int inVal)
{
  sFGPU.setGroupSize(inVal);
}

// These two will be called only if this is a dll
DLL_EX_IM void fgpuSetPrintFunc(CharArgType func)
{
  utilSetPrintFunc(func);
}

DLL_EX_IM int fgpuGetVersion(void)
{
  return GPUFRAME_VERSION;
}

#if defined(_WIN32) && defined(DELAY_LOAD_FGPU)
extern "C"
BOOL WINAPI DllMain(HINSTANCE hInstance, DWORD dwReason, LPVOID /*lpReserved*/)
{
  if (dwReason == DLL_PROCESS_ATTACH) {

    // This disables notifications about threads
    DisableThreadLibraryCalls(hInstance);
  } else if (dwReason == DLL_PROCESS_DETACH) {
  }
  return TRUE;  
}
#endif
