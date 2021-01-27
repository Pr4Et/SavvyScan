/*  multibinstat.c - Measure local mean/SD in array of boxes at multiple binnings
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 2014 by the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 * $Id$
 */
#include "hvemtypes.h"
#include "b3dutil.h"

#ifdef F77FUNCAP
#define multibinsetup MULTIBINSETUP
#define multibinstats MULTIBINSTATS
#else
#define multibinsetup multibinsetup_
#define multibinstats multibinstats_
#endif


/*!
 * Does initial setup for determining the local mean and standard deviation in an array of
 * boxes within an image volume at a series of binnings (scales).  For all of the 
 * array variables showing a dimension of 3, the index in the dimension is 0 for the X,
 * 1 for the Y, and 2 for the Z value.  ^
 * Inputs:  ^
 * [binning]         - Binning in X, Y, and Z for each scale  ^
 * [boxSize]         - Size of box to analyze in X, Y, and Z for each scale, in binned 
 * pixels   ^
 * [boxSpacing]      - Spacing between adjacent boxes in X, Y, and Z for each scale, in 
 * binned pixels  ^
 * [numScales]       - Number of binnings to analyze  ^
 * [startCoord]      - Starting unbinned coordinate of volume to analyze in X, Y, and Z  ^
 * [endCoord]        - Ending coordinate of volume to analyze in X, Y, and Z (inclusive) ^
 * Outputs:  ^
 * [boxStart]        - Starting coordinate in binned slices of first box in X, Y, and Z 
 * for each scale (equals 0 at start of slice, not [startCoord])  ^
 * [numBoxes]        - Number of boxes i nX, Y, and Z for each scale  ^
 * [bufferStartInds] - Array to receive starting index of the binned slice for each scale;
 * must be dimensioned to at least [numScales] + 1.  The size of a floating point array to
 * allocate for the loaded and binned slices is returned in the extra element.  ^
 * [statStartInds]   - Array to receive starting index of the box mean and SD values for
 * each scale; must be dimensioned to at least [numScales] + 1.  The size of floating 
 * point arrays to allocate for the means and SDs is returned in the extra element.  ^
 * Returns 1 for inappropriate coordinate limits, 2 for more than the number of scales
 * defined by MAX_MBS_SCALES, 3 for box size bigger than the binned size of the volume at 
 * some scale, or 4 for [boxSpacing] less than 1.
 */
int multiBinSetup(int binning[][3], int boxSize[][3], int boxSpacing[][3], int numScales,
                  int startCoord[3], int endCoord[3], int boxStart[][3], 
                  int numBoxes[][3], int *bufferStartInds, int *statStartInds)
{
  int dim, ind, numUBpixels[3], numBinnedPix[3];

  /* Get the number of unbinned pixels and check ranges */
  for (dim = 0; dim < 3; dim++) {
    numUBpixels[dim] = endCoord[dim] + 1 - startCoord[dim];
    if (numUBpixels[dim] < 4)
      return 1;
  }
  if (numScales > MAX_MBS_SCALES)
    return 2;

  /* Go through the scales, figure out how many binned pixels, and the number and
   * starting coordinates of the boxes in each dimension. 
   * Allow space in bottom of buffer for an array of flags for each slice */
  bufferStartInds[0] = numUBpixels[b3dX] * numUBpixels[b3dY] + numUBpixels[b3dZ] / 4 + 1;
  statStartInds[0] = 0;
  for (ind = 0; ind < numScales; ind++) {
    for (dim = 0; dim < 3; dim++) {
      numBinnedPix[dim] = numUBpixels[dim] / binning[ind][dim];
      if (numBinnedPix[dim] < boxSize[ind][dim])
        return 3;
      if (boxSpacing[ind][dim] < 1)
        return 4;
      numBoxes[ind][dim] = (numBinnedPix[dim] - boxSize[ind][dim]) / boxSpacing[ind][dim]
        + 1;
      boxStart[ind][dim] = (numBinnedPix[dim] - boxSize[ind][dim] - 
                            (numBoxes[ind][dim] - 1) * boxSpacing[ind][dim]) / 2;
    }

    /* Get the starting indexes for binned buffers and the means and SDs */
    bufferStartInds[ind + 1] = bufferStartInds[ind] + 
      numBinnedPix[b3dX] * numBinnedPix[b3dY];
    statStartInds[ind + 1] = statStartInds[ind] + 
      numBoxes[ind][b3dX] * numBoxes[ind][b3dY] * numBoxes[ind][b3dZ];
  }
  return 0;
}

/*!
 * Fortran wrapper for @@multiBinSetup@.  In Fortran, the 2-D arrays would all be 
 * dimensioned (3,*), with the X/Y/Z index in the first dimension.
 */
int multibinsetup(int binning[][3], int boxSize[][3], int boxSpacing[][3], int *numScales,
                  int *startCoord, int *endCoord, int boxStart[][3], 
                  int numBoxes[][3], int *bufferStartInds, int *statStartInds)
{
  return multiBinSetup(binning, boxSize, boxSpacing, *numScales, startCoord, endCoord,
                       boxStart, numBoxes, bufferStartInds, statStartInds);
}

/*!
 * Computes the local mean and standard deviation in an array of
 * boxes within an image volume at a series of binnings.  [binning], [boxSize], 
 * [boxSpacing], [numScales], [startCoord], [endCoord], [boxStart], [numBoxes],
 * [bufferStartInds], and [statStartInds] must be the same variables supplied to 
 * @@multiBinSetup@.  [buffer] is an array allocated to the size given by the value after 
 * [numScales] values in [bufferStartInds].  Mean and standard deviation values are
 * returned in [means] and [SDs], which must be allocated to the value after [numScales]
 * values in [statStartInds].  The callback function [getSliceFunc] is called to get the
 * one slice of data within the given coordinate range by
 * getSliceFunc(&iz, funcData, buffer), where [iz] is the Z coordinate and [funcData] is
 * filled with whatever the loading function requires (i.e., if it needs the starting and
 * ending X and Y coordinates, they must be passed in [funcData]).  Returns 
 * an error code from [getSliceFunc].
 */
int multiBinStats(int binning[][3], int boxSize[][3], int boxSpacing[][3], int numScales,
                  int startCoord[3], int endCoord[3], int boxStart[][3], 
                  int numBoxes[][3], int *bufferStartInds, int *statStartInds, 
                  float *buffer, float *means, float *SDs, int *funcData,
                  int (*getSliceFunc)(int *, int *, float *))
{
  int dim, ind, numUBpixels[3], numBoxPix[MAX_MBS_SCALES], slicesAdded[MAX_MBS_SCALES];
  int nxBin[MAX_MBS_SCALES], nyBin[MAX_MBS_SCALES];
  int scl, iz, ix, iy, ixBox, iyBox, izBox, err, ixy, zStart, zCur, zSpacing;
  int ixStart, iyStart, boxInd, bufBase, lastZused;
  float pixel, oneMean, oneSD;
  float *readBuf;
  unsigned char *needSlice;

  /* Initialize: set up the binned sizes etc. */
  for (dim = 0; dim < 3; dim++)
    numUBpixels[dim] = endCoord[dim] + 1 - startCoord[dim];
  for (ind = 0; ind < statStartInds[numScales]; ind++)
    means[ind] = SDs[ind] = 0.;
  for (scl = 0; scl < numScales; scl++) {
    nxBin[scl] = numUBpixels[b3dX] / binning[scl][b3dX];
    nyBin[scl] = numUBpixels[b3dY] / binning[scl][b3dY];
    numBoxPix[scl] = boxSize[scl][b3dX] * boxSize[scl][b3dY] * boxSize[scl][b3dZ];
  }

  /* Determine which slices are needed */
  needSlice = (unsigned char *)buffer;
  memset(needSlice, 0, numUBpixels[b3dZ]);
  for (scl = 0; scl < numScales; scl++) {
    for (izBox = 0; izBox < numBoxes[scl][b3dZ]; izBox++) {
      zStart = (boxStart[scl][b3dZ] + izBox * boxSpacing[scl][b3dZ]) * binning[scl][b3dZ];
      for (iz = zStart; iz < zStart + boxSize[scl][b3dZ] * binning[scl][b3dZ]; iz++)
        needSlice[iz] = 1;
    }
  }
  lastZused = -999;
  readBuf = buffer + numUBpixels[b3dZ] / 4 + 1;
    
  /* Loop on slices */
  for (iz = startCoord[b3dZ]; iz <= endCoord[b3dZ]; iz++) {

    /* Skip if slice is unneeded */
    if (!needSlice[iz - startCoord[b3dZ]])
      continue;

    /* If starting a batch of slices, clear out the binned buffers */
    if (iz != lastZused + 1) {
      for (scl = 0; scl < numScales; scl++) {
        slicesAdded[scl] = 0;
        for (ixy = 0; ixy < nxBin[scl] * nyBin[scl]; ixy++)
          buffer[bufferStartInds[scl] + ixy] = 0.;
      }
    }
    lastZused = iz;
    
    /* Get the next slice */
    err = getSliceFunc(&iz, funcData, readBuf);
    if (err) {
      return err;
    }
    
    /* Add it into each binned buffer */
    for (scl = 0; scl < numScales; scl++) {
      binIntoSlice(readBuf, numUBpixels[b3dX], &buffer[bufferStartInds[scl]], nxBin[scl],
                   nyBin[scl], binning[scl][b3dX], binning[scl][b3dY], 
                   1. / binning[scl][b3dZ]);
      slicesAdded[scl]++;
      if (slicesAdded[scl] >= binning[scl][b3dZ]) {
        
        /* When binned slice is done, add it into the needed boxes */
        zStart = boxStart[scl][b3dZ];
        zSpacing = boxSpacing[scl][b3dZ];
        zCur = (iz - startCoord[b3dZ]) / binning[scl][b3dZ];
        izBox = (zCur - zStart) / zSpacing;
        izBox = B3DMIN(izBox, numBoxes[scl][b3dZ] - 1);

        /* Loop backwards in Z until the slice is past the box */
        while (izBox >= 0 && zCur >= zStart && zStart + zSpacing * izBox + 
               boxSize[scl][b3dZ] - 1 >= zCur) {
          
          /* In each box, add in all the pixels in the box */
          /* TODO: find out how effective this is! */
#pragma omp parallel for  num_threads(numOMPthreads(8))                 \
  shared(numBoxes, boxSpacing, scl, boxStart, statStartInds, izBox, boxSize, means, SDs) \
  private(iyBox, iyStart, ixBox, ixStart, boxInd, iy, bufBase, ix, pixel)
          for (iyBox = 0; iyBox < numBoxes[scl][b3dY]; iyBox++) {
            iyStart = boxStart[scl][b3dY] + boxSpacing[scl][b3dY] * iyBox;
            for (ixBox = 0; ixBox < numBoxes[scl][b3dX]; ixBox++) {
              ixStart = boxStart[scl][b3dX] + boxSpacing[scl][b3dX] * ixBox;
              boxInd = statStartInds[scl] + (izBox * numBoxes[scl][b3dY] + iyBox) * 
                numBoxes[scl][b3dX] + ixBox;
              for (iy = iyStart; iy < iyStart + boxSize[scl][b3dY]; iy++) {
                bufBase = bufferStartInds[scl] + nxBin[scl] * iy;
                for (ix = ixStart; ix < ixStart + boxSize[scl][b3dX]; ix++) {
                  pixel = buffer[bufBase + ix];
                  means[boxInd] += pixel;
                  SDs[boxInd] += pixel * pixel;
                }
              }
            }
          }

          izBox--;
        }

        /* Reset for adding more slices, and increment Z */
        slicesAdded[scl] = 0;
        for (ixy = 0; ixy < nxBin[scl] * nyBin[scl]; ixy++)
          buffer[bufferStartInds[scl] + ixy] = 0.;
      }
    }
  }
   
  /* Compute the means and SDs */
  for (scl = 0; scl < numScales; scl++) {
    for (izBox = 0; izBox < numBoxes[scl][b3dZ]; izBox++) {
      for (iyBox = 0; iyBox < numBoxes[scl][b3dY]; iyBox++) {
        for (ixBox = 0; ixBox < numBoxes[scl][b3dX]; ixBox++) {
          boxInd = statStartInds[scl] + (izBox * numBoxes[scl][b3dY] + iyBox) * 
            numBoxes[scl][b3dX] + ixBox;
          sumsToAvgSD(means[boxInd], SDs[boxInd], numBoxPix[scl], &oneMean, &oneSD);
          means[boxInd] = oneMean;
          SDs[boxInd] = oneSD;
        }
      }
    }
  }
  return 0;
}

/*! Fortran wrapper to @getSliceFunc */
int multibinstats(int binning[][3], int boxSize[][3], int boxSpacing[][3], int *numScales,
                  int *startCoord, int *endCoord, int boxStart[][3], int numBoxes[][3],
                  int *bufferStartInds, int *statStartInds, float *buffer, float *means,
                  float *SDs, int *funcData, int (*getSliceFunc)(int *, int *, float *))
{
  return multiBinStats(binning, boxSize, boxSpacing, *numScales, startCoord, endCoord,
                       boxStart, numBoxes, bufferStartInds, statStartInds, buffer,
                       means, SDs, funcData, getSliceFunc);
}
