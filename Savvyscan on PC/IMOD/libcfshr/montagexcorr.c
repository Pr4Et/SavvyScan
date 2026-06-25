/*  montagexcorr.c : Functions for correlation of overlap zones between montage pieces
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 2014 by the Regents of the University of Colorado.  
 *  See dist/COPYRIGHT for full copyright notice.
 *
 * $Id$
 */

#include "imodconfig.h"
#include "b3dutil.h"

#ifdef F77FUNCAP
#define montxcbasicsizes MONTXCBASICSIZES
#define montxcindsandctf MONTXCINDSANDCTF
#define montxcfindbinning MONTXCFINDBINNING
#define montxcorredge MONTXCORREDGE
#else
#define montxcbasicsizes montxcbasicsizes_
#define montxcindsandctf montxcindsandctf_
#define montxcfindbinning montxcfindbinning_
#define montxcorredge montxcorredge_
#endif

/*!
 * Sets up most of the sizes for the overlap zone correlations.  Inputs:  ^
 * [ixy]           - 0 for X edge, 1 for Y edge ^
 * [nbin]          - Binning ^
 * [indentXC]      - Border from edge of image  ^
 * [nxyPiece]      - Size of pieces in X and Y ^
 * [nxyOverlap]    - Overlap between pieces in X and Y ^
 * [aspectMax]     - Maximum aspect ratio to correlate ^
 * [extraWidth]    - Fraction of overlap zone width to add for extra extent ^
 * [padFrac]       - Fraction of box size to pad on each edge ^
 * [niceLimit]     - Largest prime factor for FFT sizes ^
 * Outputs:  ^
 * [indentUse]     - Actual border being used ^
 * [nxyBox]        - Binned size  in X and Y of boxes to be extracted from images^
 * [numExtra]      - Number of extra pixels in X and Y ^
 * [nxPad]         - Size of padded image in X  ^
 * [nyPad]         - Size of padded image in Y  ^
 * [maxLongShift]  - Maximum shift along the long axis of the overlap zone ^
 */
void montXCBasicSizes(int ixy, int nbin,  int indentXC, int *nxyPiece, int *nxyOverlap, 
                      float aspectMax, float extraWidth, float padFrac, int niceLimit,
                      int *indentUse, int *nxyBox, int *numExtra, int *nxPad, int *nyPad,
                      int *maxLongShift)
{
  int iyx = 1 - ixy, nxyBorder[2];
  *indentUse = B3DMIN(indentXC, (nxyOverlap[ixy] - 8) / 2);
  nxyBox[ixy] = (nxyOverlap[ixy] - *indentUse * 2) / nbin;
  nxyBox[iyx] = B3DMIN((nxyPiece[iyx] - B3DMAX(2 * nbin, nxyPiece[ixy] / 20)), 
                       (int)(aspectMax * nxyOverlap[ixy])) / nbin;
  numExtra[iyx] = 0;
  numExtra[ixy] = B3DMIN(2 * (B3DNINT(extraWidth * nxyBox[ixy]) / 2), 
                         (nxyPiece[ixy] - b3dIMax(3, nbin, *indentUse, nxyPiece[ixy] / 20)
                          * 2) / nbin - nxyBox[ixy]);
  nxyBox[ixy] = nxyBox[ixy] + numExtra[ixy];
  *maxLongShift = B3DNINT(B3DMAX(1.9 * nxyOverlap[ixy] / nbin, 1.5 * nxyBox[ixy]));
  
  /* get the padded size */
  /* Limit the long dimension padding to that needed for the maximum shift */
  nxyBorder[ixy] = B3DMAX(5, B3DNINT(padFrac * nxyBox[ixy]));
  nxyBorder[iyx] = B3DMIN(B3DMAX(5, B3DNINT(padFrac * nxyBox[iyx])),
                          B3DMAX(5, B3DNINT(0.45 * *maxLongShift)));
  *nxPad = niceFrame(nxyBox[0] + 2 * nxyBorder[0], 2, niceLimit);
  *nyPad = niceFrame(nxyBox[1] + 2 * nxyBorder[1], 2, niceLimit);
}

/*!
 * Fortran wrapper for @@montXCBasicSizes@.  [ixy] should be 1 or 2.
 */
void montxcbasicsizes(int *ixy, int *nbin,  int *indentXC, int *nxyPiece, int *nxyOverlap,
                    float *aspectMax, float *extraWidth, float *padFrac, int *niceLimit,
                    int *indentUse, int *nxyBox, int *numExtra, int *nxPad, int *nyPad,
                    int *maxLongShift)
{
  montXCBasicSizes(*ixy - 1, *nbin, *indentXC, nxyPiece, nxyOverlap, *aspectMax,
                   *extraWidth, *padFrac, *niceLimit, indentUse, nxyBox, numExtra, nxPad,
                   nyPad, maxLongShift);
}

/*!
 * Sets up indices for extracting boxes, and the filter function.  Inputs:
 * [ixy]           - 0 for X edge, 1 for Y edge ^
 * [nxyPiece]      - Size of pieces in X and Y ^
 * [nxyOverlap]    - Overlap between pieces in X and Y ^
 * [nxyBox]        - Binned size  in X and Y of boxes to be extracted from images^
 * [nbin]          - Binning ^
 * [indentUse]     - Border from edge of image  ^
 * [numExtra]      - Number of extra pixels in X and Y ^
 * [nxPad]         - Size of padded image in X  ^
 * [nyPad]         - Size of padded image in Y  ^
 * [numSmooth]     - Number of pixels of smoothing if possible  ^
 * [sigma1]        - Sigma for inverted gaussian at low frequency ^
 * [sigma2]        - Sigma for high frequency filter ^
 * [radius1]       - Radius out to which low frequency filter is 0 ^
 * [radius2]       - Radius for high frequency filter ^
 * [evalCCC]       - 1 to evaluate correlation coefficients ^
 * Outputs:  ^
 * [ind0Lower]     - Starting index to extract in X and Y for lower piece ^
 * [ind1Lower]     - Ending index to extract in X and Y for lower piece ^
 * [ind0Upper]     - Starting index to extract in the [ixy] dimension for uppoer piece ^
 * [ind1Upper]     - Ending index to extract in the [ixy] dimension for uppoer piece ^
 * [nxSmooth]      - Size of edge-smoothed image in X ^
 * [nySmooth]      - Size of edge-smoothed image in Y ^
 * [ctf]           - Array of filter values; must be dimensioned at least 8193 ^
 * [delta]         - Change in radius per index in the [ctf] array ^
 */
void montXCIndsAndCTF(int ixy, int *nxyPiece, int *nxyOverlap, int *nxyBox, int nbin, 
                      int indentUse, int *numExtra, int nxPad, int nyPad, int numSmooth, 
                      float sigma1, float sigma2, float radius1, float radius2, 
                      int evalCCC, int *ind0Lower, int *ind1Lower, int *ind0Upper,
                      int *ind1Upper, int *nxSmooth, int *nySmooth, float *ctf,
                      float *delta)
{
  int iyx = 1 - ixy;
  ind0Lower[iyx] = nxyPiece[iyx] / 2 - (nbin * nxyBox[iyx]) / 2;
  ind1Lower[iyx] = ind0Lower[iyx] + nbin * nxyBox[iyx] - 1;
  ind0Lower[ixy] = nxyPiece[ixy] - nxyOverlap[ixy] + indentUse - nbin * numExtra[ixy];
  ind1Lower[ixy] = ind0Lower[ixy] + nbin * nxyBox[ixy] - 1;
  *ind0Upper = indentUse;
  *ind1Upper = indentUse + nbin * nxyBox[ixy] - 1;
      
  /* Set up smoothing over some pixels, but no more than half of the pad */
  *nxSmooth = nxyBox[0] + B3DMIN(2 * numSmooth, (nxPad - nxyBox[0]) / 2);
  *nySmooth = nxyBox[1] + B3DMIN(2 * numSmooth, (nyPad - nxyBox[1]) / 2);

  /* Multiply high-frequency filtering parameters by the binning so they are equivalent
     to frequencies in unbinned images */
  XCorrSetCTF(sigma1, nbin * sigma2, radius1, nbin * radius2, ctf, nxPad, nyPad, delta);
  if (evalCCC)
    for (iyx = 0; iyx < 8193; iyx++)
      ctf[iyx] = (float)sqrt((double)ctf[iyx]);
}

/*!
 * Fortran wrapper for @@montXCIndsAndCTF@.  [ixy] should be 1 or 2.
 */
void montxcindsandctf(int *ixy, int *nxyPiece, int *nxyOverlap, int *nxyBox, int *nbin, 
                    int *indentUse, int *numExtra, int *nxPad, int *nyPad, int *numSmooth,
                    float *sigma1, float *sigma2, float *radius1, float *radius2,
                    int *evalCCC, int *ind0Lower, int *ind1Lower, int *ind0Upper,
                    int *ind1Upper, int *nxSmooth, int *nySmooth, float *ctf,
                    float *delta)
{
  montXCIndsAndCTF(*ixy - 1, nxyPiece, nxyOverlap, nxyBox, *nbin, *indentUse, numExtra,
                   *nxPad, *nyPad, *numSmooth, *sigma1, *sigma2, *radius1, *radius2,
                   *evalCCC, ind0Lower, ind1Lower, ind0Upper, ind1Upper,
                   nxSmooth, nySmooth, ctf, delta);
}

/*!
 * Finds binning needed to keep boxed out area smaller than a target size.  Inputs:
 * [maxBin]        - Maximum binning to select ^
 * [targetSize]    - Target size to make the number of boxed pixels be less than the 
 * square of ^
 * [indentXC]      - Border from edge of image to assume ^
 * [nxyPiece]      - Size of pieces in X and Y ^
 * [nxyOverlap]    - Overlap between pieces in X and Y ^
 * [aspectMax]     - Maximum aspect ratio to correlate ^
 * [extraWidth]    - Fraction of overlap zone width to add for extra extent ^
 * [padFrac]       - Fraction of box size to pad on each edge ^
 * [niceLimit]     - Largest prime factor for FFT sizes ^
 * Outputs:  ^
 * [numPaddedPix]  - Total number of pixels for padded images, with extra allowance ^
 * [numBoxedPix]   - Total number of pixels for boxed images, with extra allowance ^
 * The return value is the binning.
 */
int montXCFindBinning(int maxBin, int targetSize, int indentXC, int *nxyPiece,
                      int *nxyOverlap, float aspectMax, float extraWidth,
                      float padFrac, int niceLimit, int *numPaddedPix, int *numBoxedPix)
{
  int nxPad, nyPad, indentUse, nxyBox[2], numExtra[2], maxLongShift;
  int ixy, nbin;
  for (nbin = 1; nbin <= maxBin; nbin++) {
    *numPaddedPix = 0;
    *numBoxedPix = 0;
    for (ixy = 0; ixy < 2; ixy++) {
      montXCBasicSizes(ixy, nbin, indentXC, nxyPiece, nxyOverlap, aspectMax, extraWidth,
                       padFrac, niceLimit, &indentUse, &nxyBox[0], &numExtra[0], &nxPad,
                       &nyPad, &maxLongShift);
      *numPaddedPix = B3DMAX(*numPaddedPix, (nxPad + 8) * (nyPad + 8));
      *numBoxedPix = B3DMAX(*numBoxedPix, (nxyBox[0] + 4) * (nxyBox[1] + 4));
    }
    if (*numBoxedPix <= targetSize * targetSize)
      return nbin;
  }
  return maxBin;
}

/*!
 * Fortran wrapper for @@montXCFindBinning@.
 */
int montxcfindbinning(int *maxBin, int *targetSize, int *indentXC, int *nxyPiece,
                      int *nxyOverlap, float *aspectMax, float *extraWidth,
                      float *padFrac, int *niceLimit, int *numPaddedPix, int *numBoxedPix)
{
  return montXCFindBinning(*maxBin, *targetSize, *indentXC, nxyPiece, nxyOverlap,
                           *aspectMax, *extraWidth, *padFrac, *niceLimit, numPaddedPix,
                           numBoxedPix);
}

/*!
 * Performs Fourier correlations and evaluation of real-space correlations to find 
 * displacement on an edge.  Whether the CCC's are evaluated depends on the value of
 * [numXcorrPeaks].  If it is 1, they are not computed, but 16 peaks are found and the 
 * first peak is taken that does not exceed the lateral shift in [maxLongShift].  If
 * it is greater than 1, at least 16 peaks are still found and checked against 
 * [maxLongShift], then the CC is evaluated at the given number of peaks. Inputs:
 * [lowerIn]       - Boxed, possibly binned image extracted from lower piece ^
 * [upperIn]       - Boxed, possibly binned image extracted from upper piece ^
 * [nxyBox]        - Size of boxed input images in X and Y ^
 * [nxyPiece]      - Size of full pieces in X and Y ^
 * [nxyOverlap]    - Overlap between pieces in X and Y ^
 * [nxSmooth]      - Size of edge-smoothed image in X ^
 * [nySmooth]      - Size of edge-smoothed image in Y ^
 * [nxPad]         - Size of padded image in X  ^
 * [nyPad]         - Size of padded image in Y  ^
 * [lowerPad]      - Temporary array for lower padded piece ^
 * [upperPad]      - Temporary array for upper padded piece ^
 * [lowerCopy]     - Temporary array for copy of lower padded piece, needed only if
 * evaluating CCC's at multiple peaks ^
 * [numXcorrPeaks] - Number of correlation peaks to examine ^
 * [legacy]        - Non-zero to do legacy correlations ^
 * [ctf]           - Array with filter function ^
 * [delta]         - Change in radius per index in the [ctf] array ^
 * [numExtra]      - Number of extra pixels in X and Y ^
 * [nbin]          - Binning ^
 * [ixy]           - 0 for X edge, 1 for Y edge ^
 * [maxLongShift]  - Maximum shift along the long axis of the overlap zone ^
 * [weightCCC]     - Non-zero to weight CCC by power function of overlap area ^
 * Outputs, functions, debugging variables: ^
 * [xDisplace]     - Displacement in X: amount to shift upper piece to align to lower ^
 * [yDisplace]     - Displacement in Y ^
 * [CCC]           - Value of CCC at the chosen peak, if computed ^
 * [twoDfft]       - 2-D FFT function to receive array, [nxPad], [nyPad], and 0 for 
 * forward or 1 for inverse ^
 * [dumpEdge]      - If non-null, a function for writing images to file.  Its arguments
 * are the image array, X dimension, X size, Y size, [ixy] + 1,  and 0 for an image or
 * 1 for a correlation  ^
 * [debugStr]      - String to receive debugging text.  To hold all text with full output,
 * it should be [numXcorrPeaks] + 1 times MONTXC_MAX_DEBUG_LINE ^
 * [debugLen]      - Length of [debugLen] ^
 * [debugLevel]    - 1 for summary of selected peak and the displacement, 2 for full 
 * output about each peak, 3 to include eliminated peaks ^
 */
void montXCorrEdge(float *lowerIn, float *upperIn, int *nxyBox, int *nxyPiece, 
                   int *nxyOverlap, int nxSmooth, int nySmooth, int nxPad, int nyPad,
                   float *lowerPad, float *upperPad, float *lowerCopy, int numXcorrPeaks,
                   int legacy, float *ctf, float delta, int *numExtra, int nbin, int ixy,
                   int maxLongShift, int weightCCC, float *xDisplace, float *yDisplace,
                   float *CCC, void (*twoDfft)(float *, int *, int *, int *),
                   void (*dumpEdge)(float *, int *, int *, int *, int *, int *), 
                   char *debugStr, int debugLen, int debugLevel)
{
  int ind, i, nxTrim, nyTrim, numPixel, indPeak, curDebugLen = 0;
  int nxPadDim = nxPad + 2;
  float *arrayIn = lowerIn;
  float *arrayOut = lowerPad;
  float xpeak[MONTXC_MAX_PEAKS], ypeak[MONTXC_MAX_PEAKS], peak[MONTXC_MAX_PEAKS];
  int zero = 0;
  int one = 1;
  int ixyP1 = ixy + 1;
  int evalCCC = (numXcorrPeaks > 1 && !legacy) ? 1 : 0;
  double ccc, cccMax, wgtCCC, fracArea;
  double overlapPow = 0.166667;

  for (ind = 0; ind < 2; ind++) {
    if (nxSmooth > nxyBox[0] && nySmooth > nxyBox[1]) {
      sliceSmoothOutPad(arrayIn, SLICE_MODE_FLOAT, nxyBox[0], nxyBox[1], arrayOut,
                        nxSmooth, nxSmooth, nySmooth);
      sliceTaperOutPad(arrayOut, SLICE_MODE_FLOAT, nxSmooth, nySmooth, arrayOut, 
                            nxPadDim, nxPad, nyPad, 0, 0.);
    } else {
      sliceTaperOutPad(arrayIn, SLICE_MODE_FLOAT, nxyBox[0], nxyBox[1], arrayOut,
                            nxPadDim, nxPad, nyPad, 0, 0.);
    }
    XCorrMeanZero(arrayOut, nxPadDim, nxPad, nyPad);
    if (dumpEdge)
      dumpEdge(arrayOut, &nxPadDim, &nxPad, &nyPad, &ixyP1, &zero);
    twoDfft(arrayOut, &nxPad, &nyPad, &zero);

    /* If filtering, apply to lower, and to upper as well if evaluating CCC's */
    if (delta > 0. && (!ind || evalCCC))
      XCorrFilterPart(arrayOut, arrayOut, nxPad, nyPad, ctf, delta);
    arrayIn = upperIn;
    arrayOut = upperPad;
  }
  if (delta > 0. && evalCCC)
    memcpy(lowerCopy, lowerPad, nxPadDim * nyPad * sizeof(float));

  /* multiply lower by complex conjugate of upper, put back in lower */
  conjugateProduct(lowerPad, upperPad, nxPad, nyPad);
  twoDfft(lowerPad, &nxPad, &nyPad, &one);
  XCorrPeakFind(lowerPad, nxPadDim, nyPad, xpeak, ypeak, peak, B3DMAX(16, numXcorrPeaks));
  
  /* Eliminate any peaks that shift beyond maximum along edge */
  /* leave indPeak pointing to first good peak */
  indPeak = -1;
  for (i = 0; i < B3DMAX(16, numXcorrPeaks); i++) {
    if ((ixy == 0 && fabs((double)ypeak[i]) > maxLongShift) ||
        (ixy == 1 && fabs((double)xpeak[i]) > maxLongShift)) {
      peak[i] = -1.e30;
      if (debugLevel > 2 && debugLen > curDebugLen + MONTXC_MAX_DEBUG_LINE) {
        sprintf(&debugStr[curDebugLen], "Eliminated peak %d at %.1f %.1f\n", i, xpeak[i],
                ypeak[i]);
        curDebugLen += strlen(&debugStr[curDebugLen]);
      }
    } else if (indPeak == -1 && peak[i] > -1.e29) {
      indPeak = i;
    }
  }
  
  /* But if no peak was legal, zero out the shift */
  if (indPeak == -1) {
    indPeak = 0;
    xpeak[0] = numExtra[0];
    ypeak[0] = numExtra[1];
  }
  *CCC = -1.5f;
  if (evalCCC) {
    
    /* If there was no filtering, simply pad images again */
    if (delta == 0) {
      sliceTaperOutPad(lowerIn, SLICE_MODE_FLOAT, nxyBox[0], nxyBox[1], lowerCopy,
                            nxPadDim, nxPad, nyPad, 0, 0.);
      sliceTaperOutPad(upperIn, SLICE_MODE_FLOAT, nxyBox[0], nxyBox[1], upperPad,
                            nxPadDim, nxPad, nyPad, 0, 0.);
    } else {
      
      /* Otherwise, back-transform the filtered images */
      twoDfft(lowerCopy, &nxPad, &nyPad, &one);
      twoDfft(upperPad, &nxPad, &nyPad, &one);
      if (dumpEdge) {
        dumpEdge(lowerCopy, &nxPadDim, &nxPad, &nyPad, &ixyP1, &zero);
        dumpEdge(upperPad, &nxPadDim, &nxPad, &nyPad, &ixyP1, &zero);
      }
    }
    cccMax = -1.5;
    for (i = 0; i < numXcorrPeaks; i++) {
      if (peak[i] > -1.e29) {
        
        /* Reject peak at no zero image offset from fixed pattern noise */
        if (!(ixy == 0 && fabs(nxyPiece[0] + nbin * (xpeak[i] - numExtra[0]) - 
                               nxyOverlap[0]) <= 3.) ||
            (ixy == 1 && fabs(nxyPiece[1] + nbin * (ypeak[i] - numExtra[1]) - 
                              nxyOverlap[1]) <= 3.)) {
          nxTrim = B3DMIN(4, nxyBox[0] / 8) + (nxPad - nxyBox[0]) / 2;
          nyTrim = B3DMIN(4, nxyBox[1] / 8) + (nyPad - nxyBox[1]) / 2;
          ccc = XCorrCCCoefficient(lowerCopy, upperPad, nxPadDim, nxPad, nyPad, 
                                   xpeak[i], ypeak[i], nxTrim, nyTrim, &numPixel);
          fracArea = (double)numPixel / ((nxPad - 2 * nxTrim) * (nyPad - 2 * nyTrim));
          wgtCCC = ccc * pow(fracArea, overlapPow);

          /* Use weighted peak if indicated, or just reject peaks with < 1/8 area */
          if (weightCCC && wgtCCC > cccMax) {
            cccMax = wgtCCC;
            indPeak = i;
          } else if (!weightCCC && ccc > cccMax && (!i || fracArea > 0.125)) {
            cccMax = ccc;
            indPeak = i;
          }
          if (debugLevel > 1 && debugLen > curDebugLen + MONTXC_MAX_DEBUG_LINE) {
            sprintf(&debugStr[curDebugLen], "%2d: at %7.1f %7.1f peak %14.7e  frac "
                    "%.3f CCC %.5f%s wgt %.5f%s\n", i, xpeak[i], ypeak[i], peak[i], 
                    fracArea, ccc, !weightCCC && indPeak == i ? "*" : " ", wgtCCC,
                    weightCCC && indPeak == i ? "*" : " ");
            curDebugLen += strlen(&debugStr[curDebugLen]);
          }
        }
      }
    }
    i = indPeak;
    if (debugLevel == 1 && debugLen > curDebugLen + MONTXC_MAX_DEBUG_LINE) {
      sprintf(&debugStr[curDebugLen], "Peak %d at %7.1f %7.1f  peak = %14.7g  CCC = "
              "%.5f\n", i, xpeak[i], ypeak[i], peak[i], cccMax);
      curDebugLen += strlen(&debugStr[curDebugLen]);
    }      
    *CCC = (float)cccMax;
  }
  if (dumpEdge)
    dumpEdge(lowerPad, &nxPadDim, &nxPad, &nyPad, &ixyP1, &one);
  
  /* return the amount to shift upper to align it to lower (verified) */
  *xDisplace = nbin * (xpeak[indPeak] - numExtra[0]);
  *yDisplace = nbin * (ypeak[indPeak] - numExtra[1]);
  if (debugLevel && debugLen > curDebugLen + MONTXC_MAX_DEBUG_LINE)
    sprintf(&debugStr[curDebugLen], "Peak at %8.2f %8.2f  Displacement %8.2f %8.2f\n",
            xpeak[indPeak], ypeak[indPeak], *xDisplace, *yDisplace);
}

/*!
 * Fortran wrapper for @@montXCorrEdge@.  The wrapper will collect and print the debug
 * strings if [debugLevel] is non-zero.
 */
void montxcorredge(float *lowerIn, float *upperIn, int *nxyBox, int *nxyPiece, 
                   int *nxyOverlap, int *nxSmooth, int *nySmooth, int *nxPad, int *nyPad,
                   float *lowerPad, float *upperPad, float *lowerCopy, int *numXcorrPeaks,
                   int *legacy, float *ctf, float *delta, int *numExtra, int *nbin, 
                   int *ixy, int *maxLongShift, int *weightCCC, float *xDisplace,
                   float *yDisplace, float *CCC, 
                   void (*twoDfft)(float *, int *, int *, int *),
                   void (*dumpEdge)(float *, int *, int *, int *, int *, int *), 
                   int *debugLevel)
{
  int debugLen = MONTXC_MAX_PEAKS * MONTXC_MAX_DEBUG_LINE;
  char debugStr[MONTXC_MAX_PEAKS * MONTXC_MAX_DEBUG_LINE];
  char *curDebug = &debugStr[0];
  char *lineEnd;
  montXCorrEdge(lowerIn, upperIn, nxyBox, nxyPiece, 
                nxyOverlap, *nxSmooth, *nySmooth, *nxPad, *nyPad,
                lowerPad, upperPad, lowerCopy, *numXcorrPeaks,
                *legacy, ctf, *delta, numExtra, *nbin, *ixy - 1,
                *maxLongShift, *weightCCC, xDisplace, yDisplace, CCC,
                twoDfft, dumpEdge, debugStr, debugLen, *debugLevel);
  if (*debugLevel) {
    while ((lineEnd = strchr(curDebug, '\n')) != NULL) {
      *lineEnd = 0x00;
      printf("%s\n", curDebug);
      curDebug = lineEnd + 1;
    }
    fflush(stdout);
  }
}
