/*
 *  findsection.cpp - Find section boundaries for various purposes
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 2014 by the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 *  $Id$
 */
#include <stdlib.h>
#include <math.h>
#include "cppdefs.h"
#include "iiunit.h"
#include "parse_params.h"
#include "hvemtypes.h"
#include "imodel.h"
#include "findsection.h"

#define MAX_DEFAULT_SCALES 12
static int getSlice(int *iz, int *fdata, float *buffer);
static void amoebaFunc(float *yvec, float *thickness);
static FindSect sFindSect;

/* 
 * Main entry
 */
int main( int argc, char *argv[])
{
  sFindSect.main(argc, argv);
  exit(0);
}

/*
 * Class constructor
 */
FindSect::FindSect()
{
  // Initializations
  mDebugOutput = 0;
  mNumThicknesses = 0;
  mFitPitchSeparately = 0;
  mScanBlockSize = -1;
  mColMedians = NULL;
  mColSlice = NULL;
  mBeadDiameter = 5;
  mBoostHighSDThickness = 0.;
  mMinNumThickForCheck = 7;   // Minimum number of thickness to do checking with
  
  // Criterion MADNs below median thickness for eliminating both points if other analysis
  // does not do either one
  mCritThickMADN = 8.;        

  // Parameters
  // 1: Minimum # of points for using robust fit to get pitch line on one surface
  mMinForRobustPitch = 6;

  // findColumnMidpoint parameters
  // 3: Number of edge MADN's above edge median that maximum value must be to proceed
  mColMaxEdgeDiffCrit = 2.;
  // 4: Fraction of maximum - edge difference to achieve
  mFracColMaxEdgeDiff = 0.5;
  // 5: Number of edge MADNs above edge to achieve as well
  mCritEdgeMADN = 3.;
  // 6: Number of box medians that need to be above those criteria
  mNumHighInsideCrit = 3;

  // fitColumnBoundaries parameters
  // 7: Number of center MADN's below the center median for inside median to be too low
  mColumnToCenMADNCrit = 5.;
  // 8: Fraction of inside - edge median difference that it must fall toward edge median
  mMaxFalloffFrac = 0.3f;
  // 9, 10: Low and high limits of range of fractions of inside - edge median difference 
  // to fit
  mLowFitFrac = 0.2f;
  mHighFitFrac = 0.8f;
  // 11: Fraction of inside - edge median difference at which to save boundary
  mBoundaryFrac = 0.5f;
  // 12: Fraction of difference at which to estimate extra boundary distance for pitch 
  // output
  mPitchBoundaryFrac = 0.25f;
  // 13: Minimum fraction of boxes in column that must yield boundaries
  mMinFracBoundsInCol = 0.5f;

  // checkBlockThicknesses parameters
  // 14: Criterion fraction of median thickness for considering block too thin
  mTooThinCrit = 0.5f;
  // 15: Drop a boundary if it is this much farther from local mean than other boundary is
  mFartherFromMeanCrit = 2.;
  // 16: Drop a boundary if its difference from the mean is this fraction of median 
  // thickness
  mMeanDiffThickFrac = 0.35f;

  // Robust fitting parameters
  // 17: K-factor for the weighting function
  mKfactor = 4.68;
  // 18: Maximum change in weights for terminatiom
  mMaxChange = 0.02;
  // 19: Maximum change in weights for terminating on an oscillation
  mMaxOscill = 0.05;
  // 20: Maximum iterations
  mMaxIter = 30;

  // 28: Fraction of beads to exclude from spread measurement
  mBeadExclPctlSpread = 0.04f;
  // 29: Basic amount to weight bead separation in the spread measure
  mBeadWeightFac = 0.33f;
  // 31: Type of data to use for layer projections: 1=mean, 2=median, 3=75th %ile,
  //     4 = fraction of boxes above criterion
  mProjLayerStatType = 0;
  // 32: Fraction of way from baseline to peak for finding rise point of layer projections
  mProjEdgeCrit = 0.1f;
  // 33: Use extrapolation to baseline rather than point where projection crosses crit
  mProjUseExtrap = 0;
  // 34: 1 to use 2nd moment of projection values, 2 to use 4th moment as spread measure
  mUseProjForSpread = 0;
  // 35: Fraction of beads to exclude from determining low and high boundaries
  mBeadExclPctlThick = 0.01f;
}

/*
 * Class main entry
 */
void FindSect::main( int argc, char *argv[])
{
  char *progname = imodProgName(argv[0]);
  int boxSize[MAX_MBS_SCALES][3];
  int boxStart[MAX_MBS_SCALES][3], startCoord[3], endCoord[3], mxyz[3];
  int numOptArgs, numNonOptArgs, numTomos, ind, ifFlip, ixyz, ierr, iz, mode, scl;
  int numBinnings, numSpacings, numSizes, tomo, funcData[5], inUnitBase = 3, nxyzTmp[3];
  int edgeExtent[3], centerExtent[3], starts[5], ends[5], ixyCen[2];
  int loop, numStat, ixBox, iyBox, size, edgeBoxes, cenBoxes, botTop, numBound, numGood;
  int yInd, numXblocks, numYblocks, iyBlock, ixBlock, major, minor, maxZeroWgt, numFit;
  int cenStarts[MAX_MBS_SCALES][3], cenEnds[MAX_MBS_SCALES][3], patchLim[2];
  int maxEdgeBoxes, maxCenBoxes, maxColumnBuf, zRange, rem, boxNum, izInside[2], izMid;
  int numIter, wgtColIn, wgtColOut, numSamples, maxPitchFit, samBlockSpace, tomoSeq;
  int maxBlocks, numTomoOpt, numBinEntries, minBoxes, maxPixels, tomoArrSize;
  int maxColBoxes, maxColSlice, boundArrSize, lowSDerror[MAX_MBS_SCALES];
  float dmin, dmax, dmean, tmin, tmax, tmean, insideMed, lastRatio, ratioDiff, lastDiff;
  float pixelDelta[3], newCell[6], origin[3], newOrigin[3], curTilt[3];
  float ratio, fitConst[2], xx, yy, xySpacing, madnFac, numBoxPerBlock;
  float fitSD[MAX_FIT_COL], fitMean[MAX_FIT_COL], fitSolution[MAX_FIT_COL];
  int defaultScales[MAX_DEFAULT_SCALES] = {1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64};
  float axisRotation, tiltXvert[4], tiltYvert[4], rotMat[3][2];
  int ifReconArea, nxSeries, nySeries, sign, surfaceLim[2], combineLim[2], colXyz[3];
  double meanSum;
  float edgeDenMeans[MAX_MBS_SCALES], cenDenMeans[MAX_MBS_SCALES];
  float fracAboveEdge[MAX_MBS_SCALES];
  int minRange[3] = {0, 0, 0};
  char *linearFormat = " %2d %8s  %7.0f  %6.1f %5.1f   %7.0f  %6.1f %5.1f  %8.2f\n";
  char *sqrtFormat = " %2d %8s  %7.1f  %6.2f %5.2f   %7.1f  %6.2f %5.2f  %8.3f\n";
  char *tableFormat = linearFormat;
  char *filename;
  char *binStr;
  bool anySamples = false;
  float *allMeans, *allSDs, *allBoundaries;
  char *volRoot = NULL;
  char *pointRoot = NULL;
  char *pitchName = NULL;
  char *surfaceName = NULL;
  char *beadFile = NULL;
  Imod *pitchModel = NULL;
  Imod *beadModel = NULL;
  IrefImage useRef, *modRefp;
  IntVec tomoInds;
  Ipoint unitPt = {1., 1., 1.};
  int sampleExtent = 0;
  int numPitchPairs = 0;
  float highSDcrit = 0.;
  int lowestSDforEdges = 0;
  int bufferStartInds[MAX_MBS_SCALES];
  MrcHeader *inHeader;
  float *writeSrc, *writeBuf, *smoothBound;
  char *lowSDerrStrings[] = 
    {"No peak in median SD value found except at top or bottom",
     "No minimum in median SD value found on one side of peak",
     "Minimum values of median SD do not occur on both sides of middle",
     "Median SD does not rise far enough above its minimum value"};

  // values used to set default center and edge areas
  float borderFrac = 0.05f;
  float edgeThickFrac = 0.025f;
  float edgeAreaFrac = 0.5f;
  float cenAreaFrac = 0.33f;

  // Parameters
  // 21: Fraction that the difference between distinguishability of center from edge
  // points must improve to adopt a higher scaling for analysis
  float cenEdgeRatioDiffCrit = 0.33f;
  // 22: Threshold weight from robust fit for including a point in the final smoothing fit
  float wgtThresh = 0.2f;
  // 23: Threshold weight from robust fit for counting a point as "good"
  float goodThresh = 0.6f;
  // 24: Fraction for percentile of positions included in auto-combine Z limits
  float combineHighPctl = 0.1f;
  // 25: Percentile of positions less than combineOutsidePix outside the Z limits
  float combineLowPctl = 0.01f;
  // 26: Number of pixels to back off from the low-percentile of positions
  int combineOutsidePix = 20;
  // 27: Fraction of depth extent to use for center samples (.1, or .4 if highSD)
  float cenThickFrac = 0.;
  // 30: Take square root of SD values
  int takeSqrt = -1;

  // Fallbacks from    ../manpages/autodoc2man 2 1 findsection
  int numOptions = 28;
  const char *options[] = {
    "tomo:TomogramFile:FNM:", "surface:SurfaceModel:FN:", "pitch:TomoPitchModel:FN:",
    "separate:SeparatePitchLineFits:B:", "samples:NumberOfSamples:I:",
    "extent:SampleExtentInY:I:", "high:HighSDboxCriterion:F:", "bead:BeadModelFile:FN:",
    "diameter:BeadDiameter:F:", "boost:BoostHighSDThickness:F:",
    "lowest:LowestSDforEdges:B:", "scales:NumberOfDefaultScales:I:",
    "binning:BinningInXYZ:ITM:", "size:SizeOfBoxesInXYZ:ITM:",
    "spacing:SpacingInXYZ:ITM:", "block:BlockSize:I:", "xminmax:XMinAndMax:IP:",
    "yminmax:YMinAndMax:IP:", "zminmax:ZMinAndMax:IP:", "flipped:ThickDimensionIsY:I:",
    "axis:AxisRotationAngle:F:", "tilt:TiltSeriesSizeXY:IP:",
    "edge:EdgeExtentInXYZ:IT:", "center:CenterExtentInXYZ:IT:",
    "control:ControlValue:FPM:", "volume:VolumeRootname:CH:", "point:PointRootname:CH:",
    "debug:DebugOutput:I:"};
  
  // Startup with fallback
  PipReadOrParseOptions(argc, argv, options, numOptions, progname, 
                        2, 1, 1, &numOptArgs, &numNonOptArgs, imodUsageHeader);

  // Get input files and make sure multiple files are same size
  ierr = PipNumberOfEntries("TomogramFile", &numTomoOpt);
  numTomos = numTomoOpt + numNonOptArgs;
  if (!numTomos)
    exitError("No input tomogram file(s) specified");

  for (ind = 0; ind < numTomos; ind++) {
    if (ind < numTomoOpt)
      PipGetString("TomogramFile", &filename);
    else
      PipGetNonOptionArg(ind - numTomoOpt, &filename);
    if (iiuOpen(inUnitBase + ind, filename, "RO"))
      exit(1);
    free(filename);
    iiuRetBasicHead(inUnitBase + ind, nxyzTmp, mxyz, &mode, &dmin, &dmax, &dmean);
    if (ind) {
      if (nxyzTmp[0] != mNxyz[0] || nxyzTmp[1] != mNxyz[1] || nxyzTmp[2] != mNxyz[2])
        exitError("All sample volumes must be the same size in X, Y, and Z");
    } else {
      mNxyz[0] = nxyzTmp[0];
      mNxyz[1] = nxyzTmp[1];
      mNxyz[2] = nxyzTmp[2];
    }
  }

  // Set the index for the thickness axis
  mThickInd = b3dZ;
  ifFlip = -1;
  PipGetInteger("ThickDimensionIsY", &ifFlip);
  if (numTomos > 1 && !ifFlip)
    exitError("Multiple tomograms must have thickness in Y dimension");
  if (ifFlip > 0 || numTomos > 1 || (ifFlip < 0 && mNxyz[b3dY] < mNxyz[b3dZ]))
    mThickInd = b3dY;
  yInd = 3 - mThickInd;

  // Get whether high SD is being done, adjust default thicknesses
  PipGetFloat("HighSDboxCriterion", &highSDcrit);
  if (highSDcrit > 0 && numTomos > 1)
    exitError("You cannot analyze for high SD with multiple tomograms");
  if (cenThickFrac <= 0.)
    cenThickFrac = highSDcrit > 0. ? 0.4f : 0.1f;
  edgeThickFrac = highSDcrit > 0. ? 0.075f : 0.025f;
  PipGetBoolean("LowestSDforEdges", &lowestSDforEdges);
  PipGetFloat("BoostHighSDThickness", &mBoostHighSDThickness);

  // Allow any parameter to be set
  PipNumberOfEntries("ControlValue", &ixyz);
  for (ind = 0; ind < ixyz; ind++) {
    PipGetTwoFloats("ControlValue", &xx, &yy);
    switch (B3DNINT(xx)) {
      SET_CONTROL_INT(1, mMinForRobustPitch);
      SET_CONTROL_FLOAT(3, mColMaxEdgeDiffCrit);
      SET_CONTROL_FLOAT(4, mFracColMaxEdgeDiff);
      SET_CONTROL_FLOAT(5, mCritEdgeMADN);
      SET_CONTROL_INT(6, mNumHighInsideCrit);
      SET_CONTROL_FLOAT(7, mColumnToCenMADNCrit);
      SET_CONTROL_FLOAT(8, mMaxFalloffFrac);
      SET_CONTROL_FLOAT(9, mLowFitFrac);
      SET_CONTROL_FLOAT(10, mHighFitFrac);
      SET_CONTROL_FLOAT(11, mBoundaryFrac);
      SET_CONTROL_FLOAT(12, mPitchBoundaryFrac);
      SET_CONTROL_FLOAT(13, mMinFracBoundsInCol);
      SET_CONTROL_FLOAT(14, mTooThinCrit);
      SET_CONTROL_FLOAT(15, mFartherFromMeanCrit);
      SET_CONTROL_FLOAT(16, mMeanDiffThickFrac);
      SET_CONTROL_FLOAT(17, mKfactor);
      SET_CONTROL_FLOAT(18, mMaxChange);
      SET_CONTROL_FLOAT(19, mMaxOscill);
      SET_CONTROL_INT(20, mMaxIter);
      SET_CONTROL_FLOAT(21, cenEdgeRatioDiffCrit);
      SET_CONTROL_FLOAT(22, wgtThresh);
      SET_CONTROL_FLOAT(23, goodThresh);
      SET_CONTROL_FLOAT(24, combineHighPctl);
      SET_CONTROL_FLOAT(25, combineLowPctl);
      SET_CONTROL_INT(26, combineOutsidePix);
      SET_CONTROL_FLOAT(27, cenThickFrac);
      SET_CONTROL_FLOAT(28, mBeadExclPctlSpread);
      SET_CONTROL_FLOAT(29, mBeadWeightFac);
      SET_CONTROL_INT(30, takeSqrt);
      SET_CONTROL_INT(31, mProjLayerStatType);
      SET_CONTROL_FLOAT(32, mProjEdgeCrit);
      SET_CONTROL_INT(33, mProjUseExtrap);
      SET_CONTROL_INT(34, mUseProjForSpread);
      SET_CONTROL_FLOAT(35, mBeadExclPctlThick);
    }
  }
  PipGetInteger("DebugOutput", &mDebugOutput);

  // Get the numbers of binnings, sizes, and spacings
  numBinnings = 1;
  ierr = PipGetInteger("NumberOfDefaultScales", &numBinnings);
  PipNumberOfEntries("BinningInXYZ", &numBinEntries);
  if (!ierr && numBinEntries)
    exitError("You cannot enter both -scalings and -binning");
  if (numBinnings > MAX_DEFAULT_SCALES)
    exitError("There are only %d default scalings available", MAX_DEFAULT_SCALES);
  numBinnings = B3DMAX(numBinnings, numBinEntries);
  if (numBinnings > MAX_MBS_SCALES)
    exitError("Too many binnings for array sizes; only %d allowed", MAX_MBS_SCALES);
  
  PipNumberOfEntries("SpacingInXYZ", &numSpacings);
  if (PipNumberOfEntries("SizeOfBoxesInXYZ", &numSizes) || !numSizes)
    exitError("Size of boxes must be entered for at least the first binning");
  if (numSizes > numBinnings || numSpacings > numBinnings)
    exitError("It makes no sense to enter more sizes or spacings than binnings");

  // Get the entries for each binning.  When there is no size entry, scale it down from
  // the last size entered.  When there is no spacing, scale it down from the last spacing
  // entered.  When there are no spacings at all, assign it the box size divided by 2
  for (scl = 0; scl < numBinnings; scl++) {

    // Binning, read in the 3 values or take the default isotropic binning
    if (scl < numBinEntries) {
      PipGetThreeIntegers("BinningInXYZ", &mBinning[scl][0], &mBinning[scl][1], 
                          &mBinning[scl][2]); 
    } else {
      for (ixyz = 0; ixyz < 3; ixyz++)
        mBinning[scl][ixyz] = defaultScales[scl];
    }

    // Size : read in the size or scale down the last size
    if (scl < numSizes) {
      PipGetThreeIntegers("SizeOfBoxesInXYZ", &boxSize[scl][0], &boxSize[scl][1],
                          &boxSize[scl][2]);
    } else {
      for (ixyz = 0; ixyz < 3; ixyz++) {
        boxSize[scl][ixyz] = B3DNINT(((float)boxSize[numSizes - 1][ixyz] *
                                      mBinning[numSizes - 1][ixyz]) /mBinning[scl][ixyz]);
        B3DCLAMP(boxSize[scl][ixyz], 1, mNxyz[ixyz] / mBinning[scl][ixyz]);
      }
      if (mDebugOutput)
        printf("Scale %d  box size %d %d %d\n", scl + 1, boxSize[scl][0], boxSize[scl][1],
               boxSize[scl][2]);
    }
    for (ixyz = 0; ixyz < 3; ixyz++)
      minRange[ixyz] = B3DMAX(minRange[ixyz], boxSize[scl][ixyz] * mBinning[scl][ixyz]);

    // Spacing: read it in, or set up initial one, or scale down the last one
    if (scl < numSpacings) {
      PipGetThreeIntegers("SpacingInXYZ", &mBoxSpacing[scl][0], &mBoxSpacing[scl][1],
                          &mBoxSpacing[scl][2]);
    } else if (!scl) {
      for (ixyz = 0; ixyz < 3; ixyz++)
        mBoxSpacing[scl][ixyz] = B3DMAX(1, boxSize[scl][ixyz] / 2);
      numSpacings = 1;
    } else {  
      for (ixyz = 0; ixyz < 3; ixyz++)
        mBoxSpacing[scl][ixyz] = 
          B3DMAX(1, B3DNINT(((float)mBoxSpacing[numSpacings - 1][ixyz] *
                             mBinning[numSpacings - 1][ixyz]) / mBinning[scl][ixyz]));
      if (mDebugOutput)
        printf("Scale %d  box spacing %d %d %d\n", scl + 1, mBoxSpacing[scl][0],
               mBoxSpacing[scl][1], mBoxSpacing[scl][2]);
    }
  }
  
  // Get output file names and check for validity
  PipGetString("VolumeRootname", &volRoot);
  PipGetString("PointRootname", &pointRoot);
  PipGetString("SurfaceModel", &surfaceName);
  PipGetBoolean("SeparatePitchLineFits", &mFitPitchSeparately);
  PipGetString("TomoPitchModel", &pitchName);
  if (takeSqrt < 0)
    takeSqrt = highSDcrit > 0 ? 1 : 0;
  PipGetString("BeadModelFile", &beadFile);
  if (surfaceName && numTomos > 1)
    exitError("You cannot output a surface model with multiple input tomograms");
  ierr = PipGetInteger("NumberOfSamples", &numSamples);
  if (!ierr && (numTomos > 1 || highSDcrit > 0.))
    exitError("You cannot specify sampling with multiple input tomograms or analysis "
              "for high SD");
  if (ierr && pitchName && highSDcrit <= 0. && numTomos == 1)
    exitError("You must specify the number of samples for a boundary model from a "
              "single tomogram");
  ierr = PipGetInteger("SampleExtentInY", &sampleExtent);
  if (!ierr && (numTomos > 1 || highSDcrit > 0.))
    exitError("You cannot specify sample extent with multiple tomograms or analysis "
              "for high SD"); 
  if (!ierr && !pitchName)
    exitError("You cannot specify sample extent unless outputting a model for "
              "tomopitch"); 

  // Get bead model and transform it to this volume
  if (beadFile) {
    beadModel = imodRead(beadFile);
    if (!beadModel)
      exitError("Reading in bead model file %s", beadModel);
    modRefp = beadModel->refImage;
    if (!modRefp || !(beadModel->flags & IMODF_OTRANS_ORIGIN))
      exitError("Bead model does not contain coordinate information needed to transform "
                "to volume being analyzed");
    PipGetFloat("BeadDiameter", &mBeadDiameter);
    iiuRetOrigin(inUnitBase, &useRef.ctrans.x, &useRef.ctrans.y, &useRef.ctrans.z);
    iiuRetDelta(inUnitBase, pixelDelta);
    iiuRetTilt(inUnitBase, curTilt);
    useRef.crot.x = curTilt[0];
    useRef.crot.y = curTilt[1];
    useRef.crot.z = curTilt[2];
    useRef.cscale.x = pixelDelta[0];
    useRef.cscale.y = pixelDelta[1];
    useRef.cscale.z = pixelDelta[2];
    useRef.otrans = modRefp->ctrans;
    useRef.orot = modRefp->crot;
    useRef.oscale = modRefp->cscale;
    imodTransFromRefImage(beadModel, &useRef, unitPt);
  }

  borderFrac = surfaceName ? 0.025f : 0.05f;

  // Initialize coordinate limits with a border on the other two axes; get limits
  for (ind = 0; ind < 3; ind++) {
    startCoord[ind] = ind == mThickInd ? 0 : B3DNINT(borderFrac * mNxyz[ind]);
    if (mNxyz[ind] - 2 * startCoord[ind] < minRange[ind])
      startCoord[ind] = (mNxyz[ind] - minRange[ind]) / 2;
    endCoord[ind] = mNxyz[ind] - 1 - startCoord[ind];
  }
  PipGetTwoIntegers("XMinAndMax", &startCoord[0], &endCoord[0]);
  PipGetTwoIntegers("YMinAndMax", &startCoord[1], &endCoord[1]);
  PipGetTwoIntegers("ZMinAndMax", &startCoord[2], &endCoord[2]);
  for (ind = 0; ind < 3; ind++)
    if (startCoord[ind] < 0 || endCoord[ind] >= mNxyz[ind])
      exitError("Starting or ending %c coordinate outside of range for volume%s",
                'X' + ind, mBinning[0][ind] > 1 ?
                "; enter box size in binned coordinates" : "");

  if (mDebugOutput)
    printf("Analyzing X %d %d  Y %d %d  Z %d %d\n", startCoord[0], endCoord[0],
           startCoord[1], endCoord[1], startCoord[2], endCoord[2]);

  // Get rotation angle and stack size and make a boundary of reconstructable region
  nxSeries = mNxyz[b3dX];
  nySeries = mNxyz[yInd];
  ierr = PipGetTwoIntegers("TiltSeriesSizeXY", &nxSeries, &nySeries);
  ifReconArea = 1 - PipGetFloat("AxisRotationAngle", &axisRotation);
  if (!ierr && !ifReconArea)
    exitError("Axis rotation angle must also be entered with tilt series size");
  if (ifReconArea) {
    tiltXvert[0] = mNxyz[b3dX] / 2. - nxSeries / 2.;
    tiltYvert[0] = mNxyz[yInd] / 2. - nySeries / 2.;
    tiltXvert[1] = mNxyz[b3dX] / 2. + nxSeries / 2.;
    tiltYvert[1] = tiltYvert[0];
    tiltXvert[2] = tiltXvert[1];
    tiltYvert[2] = mNxyz[yInd] / 2. + nySeries / 2.;
    tiltXvert[3] = tiltXvert[0];
    tiltYvert[3] = tiltYvert[2];
    rotMat[0][0] = cos(axisRotation * RADIANS_PER_DEGREE);
    rotMat[1][0] = sin(axisRotation * RADIANS_PER_DEGREE);
    rotMat[0][1] = -rotMat[1][0];
    rotMat[1][1] = rotMat[0][0];
    rotMat[2][0] = rotMat[2][1] = 0.;
    for (ind = 0; ind < 4; ind++) 
      xfApply(&rotMat[0][0], mNxyz[b3dX] / 2., mNxyz[yInd] / 2., tiltXvert[ind],
              tiltYvert[ind], &tiltXvert[ind], &tiltYvert[ind], 2);
  }

  // set up and get the extent of edge and center analysis
  for (ixyz = 0; ixyz < 3; ixyz++) {
    size = endCoord[ixyz] + 1 - startCoord[ixyz];
    if (ixyz == mThickInd) {
      edgeExtent[ixyz] = B3DMAX(1, B3DNINT(edgeThickFrac * size));
      centerExtent[ixyz] = B3DMAX(1, B3DNINT(cenThickFrac * size));
    } else {
      edgeExtent[ixyz] = B3DMAX(1, B3DNINT(edgeAreaFrac * size));
      centerExtent[ixyz] = B3DMAX(1, B3DNINT(cenAreaFrac * size));
    }
  }
  PipGetThreeIntegers("EdgeExtentInXYZ", &edgeExtent[0], &edgeExtent[1], &edgeExtent[2]);
  PipGetThreeIntegers("CenterExtentInXYZ", &centerExtent[0], &centerExtent[1],
                      &centerExtent[2]);
  if (mDebugOutput)
    printf("Extents edge %d %d %d  center %d %d %d\n", edgeExtent[0], edgeExtent[1],
           edgeExtent[2], centerExtent[0], centerExtent[1], centerExtent[2]);

  // Set up the computation
  ierr = multiBinSetup(mBinning, boxSize, mBoxSpacing, numBinnings, startCoord, endCoord,
                       boxStart, mNumBoxes, bufferStartInds, mStatStartInds);
  if (ierr)
    exitError("%s", ierr == 1 ? "Coordinate limits are not usable" : 
              "Box size is too large for binned volume size");

  // Set default scanning size if not entered, then make bigger if # of boxes is limited 
  // in one direction, to give equivalent number of boxes in block
  PipGetInteger("BlockSize", &mScanBlockSize);
  if (mScanBlockSize < 0)
    mScanBlockSize = surfaceName ? 200. : 100.;
  minBoxes = B3DMIN(mNumBoxes[0][b3dX], mNumBoxes[0][yInd]);
  maxPixels = B3DMAX(boxSize[0][b3dX] * mBinning[0][b3dX], 
                     boxSize[0][yInd] * mBinning[0][yInd]);
  numBoxPerBlock = (float)mScanBlockSize / maxPixels;
  if (minBoxes < numBoxPerBlock)
    mScanBlockSize = maxPixels * B3DNINT(numBoxPerBlock * numBoxPerBlock / minBoxes);

  // Get sizes for center and edge samples and make sure array will be big enough
  maxEdgeBoxes = 0;
  maxCenBoxes = 0;
  maxPitchFit = 0;
  maxBlocks = 0;
  maxColBoxes = 0;
  maxColSlice = 0;
  for (scl = 0; scl < numBinnings; scl++) {
    edgeBoxes = 2;
    cenBoxes = 1;
    for (ixyz = 0; ixyz < 3; ixyz++) {
      edgeBoxes *= B3DMAX(1, B3DNINT(((float)edgeExtent[ixyz] / mBinning[scl][ixyz] - 
                                      boxSize[scl][ixyz]) / mBoxSpacing[scl][ixyz] + 1.));
      numStat = B3DNINT(((float)centerExtent[ixyz] / mBinning[scl][ixyz] - 
                         boxSize[scl][ixyz]) / mBoxSpacing[scl][ixyz] + 1.);
      B3DCLAMP(numStat, 1, mNumBoxes[scl][ixyz]);
      cenStarts[scl][ixyz] = (mNumBoxes[scl][ixyz] - numStat) / 2;
      cenEnds[scl][ixyz] = cenStarts[scl][ixyz] + numStat - 1;
      cenBoxes *= numStat;
    }
    invertYifFlipped(cenStarts[scl][1], cenEnds[scl][1], &mNumBoxes[scl][0]);
    maxEdgeBoxes = B3DMAX(maxEdgeBoxes, edgeBoxes);
    maxCenBoxes = B3DMAX(maxCenBoxes, cenBoxes);
    if (mDebugOutput)
      printf("%d cen s-e x %d %d  y %d %d  z %d %d\n", scl, cenStarts[scl][0],
             cenEnds[scl][0], cenStarts[scl][1], cenEnds[scl][1], cenStarts[scl][2],
             cenEnds[scl][2]);

    setupBlocks(mNumBoxes[scl][b3dX], scl, b3dX);
    setupBlocks(mNumBoxes[scl][yInd], scl, yInd);
    maxBlocks = B3DMAX(maxBlocks, (mNumBlocks[scl][b3dX] + 2) * 
                       (mNumBlocks[scl][yInd] + 2));

    // fitColumnBoundaries needs space for a block of values on each surface and for
    // two fitting arrays that could be an unknown fraction of the Z extent
    maxColumnBuf = B3DMAX(maxColumnBuf, 2 * mNumBoxes[scl][mThickInd] + 4 *
                          (mNumInBlock[scl][b3dX] + 3) * (mNumInBlock[scl][yInd] + 3));

    // Make sure there is enough space for fitting lines for pitch model
    if (pitchName) {
      if (numTomos > 1) {
        size = mNumBlocks[scl][yInd];
      } else if (!sampleExtent) {
        size = 1;
      } else {
        size = mBinning[scl][yInd] * (mNumInBlock[scl][yInd] * mBoxSpacing[scl][yInd] +
                                      boxSize[scl][yInd]);
        size = B3DMAX(1, B3DNINT((float)sampleExtent / size));
      }
      maxPitchFit = B3DMAX(maxPitchFit, size * mNumBlocks[scl][b3dX]);
    }

    // Keep track of biggest thickness and biggest box slice in X-Z for column outputs
    maxColBoxes = B3DMAX(maxColBoxes, mNumBoxes[scl][mThickInd]);
    maxColSlice = B3DMAX(maxColSlice, mNumBoxes[scl][b3dX]*mNumBoxes[scl][mThickInd]);
  }

  size = b3dIMax(6, bufferStartInds[numBinnings], 2 * maxEdgeBoxes, 2 * maxCenBoxes, 
                 B3DMAX(maxCenBoxes, maxColumnBuf) + maxColumnBuf, 3 * maxPitchFit, 
                 maxBlocks);
  mBuffer = B3DMALLOC(float, size * numTomos);
  tomoArrSize = mStatStartInds[numBinnings];
  allMeans = mMeans = B3DMALLOC(float, tomoArrSize * numTomos);
  allSDs = mSDs = B3DMALLOC(float, tomoArrSize * numTomos);
  mThicknesses = B3DMALLOC(float, maxBlocks * numTomos);
  if (volRoot) {
    mColSlice = B3DMALLOC(float, maxColSlice + maxColBoxes);
    mColMedians = mColSlice + maxColSlice;
  }
  if (!mBuffer || !mMeans || !mSDs || !mThicknesses || (volRoot && !mColSlice))
    exitError("Allocating arrays for image data and statistics");

  // Set up a tomopitch model
  if (pitchName) {
    pitchModel = imodNew();
    if (!pitchModel || imodNewObject(pitchModel))
      exitError("Creating model for tomopitch output");
    pitchModel->obj[0].extra[IOBJ_EX_PNT_LIMIT] = 2;
    pitchModel->obj[0].flags = IMOD_OBJFLAG_OPEN | IMOD_OBJFLAG_PLANAR;
    if (numTomos > 1)
      pitchModel->obj[0].flags |= IMOD_OBJFLAG_TIME;
    pitchModel->xmax = mNxyz[0];
    pitchModel->ymax = mNxyz[1];
    pitchModel->zmax = mNxyz[2];
  }

  // Set up sequences indexes to do tomograms from center outward
  for (tomoSeq = 0; tomoSeq < numTomos; tomoSeq++) {
    tomo = numTomos / 2;
    if (tomoSeq % 2)
      tomo -= (tomoSeq + 1) / 2;
    else
      tomo += tomoSeq / 2;
    tomoInds.push_back(tomo);
  }

  // START LOOP ON TOMOGRAMS
  mNumExtraSum = 0;
  mExtraPitchSum = 0.;
  for (tomoSeq = 0; tomoSeq < numTomos; tomoSeq++) {
    tomo = tomoInds[tomoSeq];
    mSDs = &allSDs[tomo * tomoArrSize];
    mMeans = &allMeans[tomo * tomoArrSize];

    funcData[0] = inUnitBase + tomo;
    funcData[1] = startCoord[0];
    funcData[2] = endCoord[0];
    funcData[3] = startCoord[1];
    funcData[4] = endCoord[1];

    ierr = multiBinStats(mBinning, boxSize, mBoxSpacing, numBinnings, startCoord,
                         endCoord, boxStart, mNumBoxes, bufferStartInds, mStatStartInds,
                         mBuffer, mMeans, mSDs, funcData, getSlice);
    if (ierr)
      exitError("Reading data from file");

    if (takeSqrt) {
      tableFormat = sqrtFormat;
      for (scl = 0; scl < numBinnings; scl++) {
        for (iz = 0; iz < mNumBoxes[scl][b3dZ]; iz++) {
          writeBuf = mSDs + mStatStartInds[scl] + iz *  mNumBoxes[scl][b3dY] * 
            mNumBoxes[scl][b3dX];
          for (iyBox = 0; iyBox < mNumBoxes[scl][b3dY]; iyBox++)  {
            for (ixBox = 0; ixBox < mNumBoxes[scl][b3dX]; ixBox++) {
              ind = ixBox + iyBox * mNumBoxes[scl][b3dX];
              writeBuf[ind] = sqrt(writeBuf[ind]);
            }
          }
        }
      }
    }

    // Write data
    if (volRoot) {
      writeSrc = mMeans;
      iiuRetDelta(inUnitBase + tomo, pixelDelta);
      iiuRetOrigin(inUnitBase + tomo, &origin[0], &origin[1], &origin[2]);
      iiuRetTilt(inUnitBase + tomo, curTilt);
      for (loop = 0; loop < 2; loop++) {
        for (scl = 0; scl < numBinnings; scl++) {

          // Adjust origin and pixel size
          for (ixyz = 0; ixyz < 3; ixyz++) {
            newCell[ixyz] = mNumBoxes[scl][ixyz] * pixelDelta[ixyz] * 
              mBinning[scl][ixyz] * mBoxSpacing[scl][ixyz]; 
            newOrigin[ixyz] = origin[ixyz] - pixelDelta[ixyz] * 
              (startCoord[ixyz] + boxStart[scl][ixyz] * mBinning[scl][ixyz]);
            newCell[ixyz + 3] = 90.;
          }
          sprintf((char *)mBuffer, "%s%d-scale%d.%s", volRoot, tomo, scl,
                  loop ? "SDs" : "means");
          iiuOpen(1, (char *)mBuffer, "NEW");
          iiuCreateHeader(1, &mNumBoxes[scl][0], &mNumBoxes[scl][0], 2, NULL, 0);
          iiuAltCell(1, newCell);
          iiuAltOrigin(1, newOrigin[0], newOrigin[1], newOrigin[2]);
          iiuAltTilt(1, curTilt);
          dmin = 1.e37;
          dmax = -dmin;
          dmean = 0;
          for (iz = 0; iz < mNumBoxes[scl][b3dZ]; iz++) {
            writeBuf = writeSrc + mStatStartInds[scl] + iz *  mNumBoxes[scl][b3dY] * 
              mNumBoxes[scl][b3dX];

            iiuWriteSection(1, (char *)writeBuf);
            arrayMinMaxMean(writeBuf, mNumBoxes[scl][b3dX], mNumBoxes[scl][b3dY], 0,
                            mNumBoxes[scl][b3dX] - 1, 0 , mNumBoxes[scl][b3dY] - 1, &tmin,
                            &tmax, &tmean);
            dmin = B3DMIN(dmin, tmin);
            dmax = B3DMAX(dmax, tmax);
            dmean += tmean / mNumBoxes[scl][b3dZ];
          }
          iiuWriteHeaderStr(1, "FINDSECTION", 0, dmin, dmax, dmean);
          iiuClose(1);
        }          
        writeSrc = mSDs;
      } 

      // Set up output file for column medians
      for (scl = 0; scl < numBinnings; scl++) {
        sprintf((char *)mBuffer, "%s%d-scale%d.colmed", volRoot, tomo, scl);
        iiuOpen(1, (char *)mBuffer, "NEW");
        setupBlocks(mNumBoxes[scl][b3dX], scl, b3dX);
        setupBlocks(mNumBoxes[scl][yInd], scl, yInd);
        printf("%d: %d x blocks of %d+   %d y blocks of %d+\n", scl,mNumBlocks[scl][b3dX],
               mNumInBlock[scl][0], mNumBlocks[scl][yInd], mNumInBlock[scl][yInd]);
        colXyz[0] = mNumBoxes[scl][0];
        colXyz[1] = mNumBoxes[scl][mThickInd];
        colXyz[2] = mNumBlocks[scl][yInd];
        iiuCreateHeader(1, &colXyz[0], &colXyz[0], 2, NULL, 0);
        dmin = 1.e37;
        dmax = -dmin;
        dmean = 0;
        
        for (iyBox = 0; iyBox < mNumBlocks[scl][yInd]; iyBox++) {
          for (ixBox = 0; ixBox < mNumBlocks[scl][b3dX]; ixBox++) {
            
            // Get the range of boxes in the block
            for (ixyz = 0; ixyz <= yInd; ixyz += yInd)
              balancedGroupLimits(mNumBoxes[scl][ixyz], mNumBlocks[scl][ixyz], 
                                  ixyz ? iyBox : ixBox, &starts[ixyz], &ends[ixyz]);
            
            // Get the medians
            findColumnMidpoint(starts[0], ends[0], starts[yInd], ends[yInd], scl,
                               &mBuffer[maxCenBoxes], izInside, izMid, insideMed);
            
            // Replicate medians into the boxes within the slice
            for (iz = 0; iz < mNumBoxes[scl][mThickInd]; iz++)
              for (ind = starts[0]; ind <= ends[0]; ind++)
                mColSlice[iz * mNumBoxes[scl][0] + ind] = mColMedians[iz];
            
          }
          
          // Write slice for this Y
          iiuWriteSection(1, (char *)mColSlice);
          arrayMinMaxMean(mColSlice, colXyz[0], colXyz[1], 0, colXyz[0] - 1, 0, 
                          colXyz[1] - 1, &tmin, &tmax, &tmean);
          dmin = B3DMIN(dmin, tmin);
          dmax = B3DMAX(dmax, tmax);
          dmean += tmean / mNumBlocks[scl][yInd];
        }
        iiuWriteHeaderStr(1, "FINDSECTION", 0, dmin, dmax, dmean);
        iiuClose(1);
      }
    }
  }

  // Evaluate statistics of edge first
  for (scl = 0; scl < numBinnings; scl++) {
    for (ixyz = 0; ixyz < 3; ixyz++) {
      numStat = B3DNINT(((float)edgeExtent[ixyz] / mBinning[scl][ixyz] - 
                         boxSize[scl][ixyz]) / mBoxSpacing[scl][ixyz] + 1.);
      B3DCLAMP(numStat, 1, mNumBoxes[scl][ixyz]);
      
      // For thickness, get inclusive limits of excluded region
      if (ixyz == mThickInd)
        numStat = mNumBoxes[scl][ixyz] - numStat;
      starts[ixyz] = (mNumBoxes[scl][ixyz] - numStat) / 2;
      ends[ixyz] = starts[ixyz] + numStat - 1;
    }
    invertYifFlipped(starts[1], ends[1], &mNumBoxes[scl][0]);

    // Save excluded region in the thickness dimension in s/e[3] and the absolute
    // limits in s/e[4], set fallbacks for edge limits
    starts[3] = starts[mThickInd];
    ends[3] = ends[mThickInd];
    starts[4] = 0;
    ends[4] = mNumBoxes[scl][mThickInd] - 1;
    mBestLowEdge[scl] = 0;
    mBestHighEdge[scl] = ends[4];

    if (highSDcrit > 0. && lowestSDforEdges) {
      lowSDerror[scl] = findLowestSDEdges(starts, ends, scl);
      if (lowSDerror[scl] && mDebugOutput)
        printf("Error %d finding best limits for edge at scale %d\n", ierr, scl + 1);
    }
    
    numStat = 0;
    meanSum = 0.;
    if (mDebugOutput)
      printf("edge boxes %d %d %d %d %d %d %d %d\n", starts[0], ends[0], starts[1], 
             ends[1], starts[2], ends[2], starts[3], ends[3]);
    
    // Now add each tomogram into edge sample
    for (tomoSeq = 0; tomoSeq < numTomos; tomoSeq++) {
      tomo = tomoInds[tomoSeq];
      mSDs = &allSDs[tomo * tomoArrSize];
      mMeans = &allMeans[tomo * tomoArrSize];
      
      // Set start and end of lower sample from absolute start and one below excluded
      starts[mThickInd] = starts[4];
      ends[mThickInd] = starts[3] - 1;
      addBoxesToSample(starts[0], ends[0], starts[1], ends[1], starts[2], ends[2], scl,
                       numStat, meanSum);
      
      starts[mThickInd] = ends[3] + 1;
      ends[mThickInd] = ends[4];
      if (mDebugOutput)
        printf("AND  boxes %d %d %d %d %d %d\n", starts[0], ends[0], starts[1], ends[1],
               starts[2], ends[2]);
      addBoxesToSample(starts[0], ends[0], starts[1], ends[1], starts[2], ends[2], scl,
                       numStat, meanSum);
    }

    edgeDenMeans[scl] = meanSum / numStat;
    rsFastMedian(mBuffer, numStat, &mBuffer[numStat], &mEdgeMedians[scl]);
    rsFastMADN(mBuffer, numStat, mEdgeMedians[scl], &mBuffer[numStat], 
               &mEdgeMADNs[scl]);
    /*printf("cen: %6.1f %5.1f  %6.1f %5.1f   edge: %6.1f %5.1f  %6.1f %5.1f", 
      edgCenMean[0], edgCenSD[0], medians[scl][0], MADNs[scl][0],
      edgCenMean[1], edgCenSD[1], medians[scl][1], MADNs[scl][1]); */
    if (mDebugOutput)
      printf("edge: %d  %8.3f %8.3f\n", numStat, mEdgeMedians[scl], mEdgeMADNs[scl]);
    /*  for (loop = 0; loop < 4; loop++)
        printf("  %.3f",(medians[scl][0] - (1. + 0.5 * loop) * MADNs[scl][0] - 
        medians[scl][1]) / MADNs[scl][1]);*/
    
    //printf("   %.3f", (edgCenMean[0] - 2. * edgCenSD[0] - edgCenMean[1])/edgCenSD[1]);
  }
  if (mDebugOutput)
    printf("\n");

  printf("               center   center SDs      edge     edge SDs      distinct-\n"
         "scale  binning  mean   median  MADN     mean   median  MADN      ness\n");
  for (scl = 0; scl < numBinnings; scl++) {
      
    // Get the block size and number of blocks for center sampling
    setupBlocks(cenEnds[scl][b3dX] + 1 - cenStarts[scl][b3dX], scl, b3dX);
    setupBlocks(cenEnds[scl][yInd] + 1 - cenStarts[scl][yInd], scl, yInd);
    zRange = cenEnds[scl][mThickInd] + 1 - cenStarts[scl][mThickInd];

    mCenMedians[scl] = 0.;    
    numStat = 0;
    meanSum = 0.;
    /*printf("%d: %d x blocks of %d+   %d y blocks of %d+\n", scl,mNumBlocks[scl][b3dX],
      mNumInBlock[scl][0], mNumBlocks[scl][yInd], mNumInBlock[scl][yInd]);*/
    
    // Do each tomogram
    for (tomoSeq = 0; tomoSeq < numTomos; tomoSeq++) {
      tomo = tomoInds[tomoSeq];
      mSDs = &allSDs[tomo * tomoArrSize];
      mMeans = &allMeans[tomo * tomoArrSize];

      // Loop on the blocks
      for (iyBox = 0; iyBox < mNumBlocks[scl][yInd]; iyBox++) {
        for (ixBox = 0; ixBox < mNumBlocks[scl][b3dX]; ixBox++) {
          
          // Get the range of boxes in the block
          for (ixyz = 0; ixyz <= yInd; ixyz += yInd) {
            boxNum = ixyz ? iyBox : ixBox;
            rem = (cenEnds[scl][ixyz] + 1 - cenStarts[scl][ixyz]) % mNumBlocks[scl][ixyz];
            starts[ixyz] = cenStarts[scl][ixyz] + boxNum * mNumInBlock[scl][ixyz] + 
              (boxNum <= rem ? boxNum : rem);
            ends[ixyz] = starts[ixyz] + mNumInBlock[scl][ixyz] + (boxNum < rem ? 0 : -1);
          }
          /*printf("block %d %d  x %d %d y %d %d\n", ixBox, iyBox, starts[0], ends[0], 
            starts[yInd], ends[yInd]);*/
          
          if (highSDcrit > 0.) {
            
            // Shift the range to be centered on the best edge limits
            if (lowestSDforEdges) {
              starts[mThickInd] = B3DMAX((mBestLowEdge[scl] + mBestHighEdge[scl] - zRange)
                                         / 2, 0);
              ends[mThickInd] = B3DMIN(starts[mThickInd] + zRange - 1, 
                                       mNumBoxes[scl][mThickInd]);
            }
            addBoxesToSample(starts[0], ends[0], starts[1], ends[1], starts[2], ends[2],
                             scl, numStat, meanSum);
              
            // Get the midpoint and ranges inside and add boxes
          } else if (!findColumnMidpoint(starts[0], ends[0], starts[yInd], ends[yInd], 
                                         scl, &mBuffer[maxCenBoxes], izInside, izMid,
                                         insideMed)) {
            starts[mThickInd] = B3DMAX(izInside[0], izMid - zRange / 2);
            ends[mThickInd] = B3DMIN(izInside[1], starts[mThickInd] + zRange - 1);
            addBoxesToSample(starts[0], ends[0], starts[1], ends[1], starts[2], ends[2],
                             scl, numStat, meanSum);
          }
        }
      }
    }
    if (numStat < 10)
        continue;

    cenDenMeans[scl] = meanSum / numStat;
    anySamples = true;

    rsFastMedian(mBuffer, numStat, &mBuffer[numStat], &mCenMedians[scl]);
    rsFastMADN(mBuffer, numStat, mCenMedians[scl], &mBuffer[numStat], 
               &mCenMADNs[scl]);
    madnFac = (mCenMedians[0] - mEdgeMedians[0]) / mCenMADNs[0];
    B3DCLAMP(madnFac, 1., 2.);
    fracAboveEdge[scl] = (mCenMedians[scl] - madnFac * mCenMADNs[scl] - 
                          mEdgeMedians[scl]) / mEdgeMADNs[scl];

    // If doing high SD, find fraction of center boxes distinct from the edge
    if (highSDcrit > 0.) {
      numGood = 0;
      for (ind = 0; ind < numStat; ind++)
        if ((mBuffer[ind] - mEdgeMedians[scl]) / mEdgeMADNs[scl] > madnFac)
          numGood++;
      fracAboveEdge[scl] = ((float)numGood) / numStat;
    }
    binStr = (char *)&mBuffer[0];
    sprintf(binStr, "%d,%d,%d", mBinning[scl][0], mBinning[scl][1], 
            mBinning[scl][1]);
    printf(tableFormat, scl + 1, binStr,
           cenDenMeans[scl], mCenMedians[scl], mCenMADNs[scl], edgeDenMeans[scl], 
           mEdgeMedians[scl], mEdgeMADNs[scl], fracAboveEdge[scl]);
  }

  if (!anySamples)
    exitError("%soo few boxes in center sample where boundaries of section could "
              "be detected", numBinnings > 1 ? "For all scalings, t" : "T");
  
  // Next, pick the best binning
  mBestScale = 0;
  lastRatio = 0.;
  ratioDiff = 0.;
  for (scl = 0; scl < numBinnings; scl++) {
    ratio = fracAboveEdge[scl];
    if (mCenMedians[scl] == 0.) {
      printf("WARNING: Too few boxes in center sample where boundaries of section could"
             " be detected for scaling # %d\n", scl + 1);
      continue;
    }

    // Require SOME distinction between edge and center
    if (scl && (highSDcrit > 0 || 
                (mCenMedians[scl] - mEdgeMedians[scl]) / mEdgeMADNs[scl] > 0.05)) {

      // Compute ratio change per step since the last best scale
      ratioDiff = (ratio - lastRatio) / (scl - mBestScale);
      if (ratioDiff < cenEdgeRatioDiffCrit * lastDiff)
        continue;
      mBestScale = scl;
    } 
    lastRatio = ratio;
    lastDiff = ratioDiff;
  }
  if (numBinnings > 1)
    printf("Selected scaling # %d as the best one for analysis\n", mBestScale + 1);

  // report fate of minimum SD edge picking for that scaling
  if (highSDcrit > 0. && lowestSDforEdges) {
    if (lowSDerror[mBestScale])
      printf("For this scaling, an error occurred finding slices with a\n"
             "   value of median SD to use for edge samples with this scaling:\n   %s\n"
             , lowSDerrStrings[lowSDerror[mBestScale]]);
    else
      printf("The edge statistics for this scaling were computed from boxes around\n"
             "   %d and %d (out of %d) instead of from top and bottom boxes\n", 
             mBestLowEdge[mBestScale], mBestHighEdge[mBestScale], 
             mNumBoxes[mBestScale][mThickInd]);
  }

  // Do everything for cryo (highSD) analysis and exit
  if (highSDcrit > 0) {
    analyzeHighSD(highSDcrit, beadModel, &boxSize[mBestScale][0], startCoord,
                  endCoord, pitchModel);
    if (pointRoot || pitchName)
      inHeader = iiuMrcHeader(inUnitBase + tomo, "findsection", 1, 0);
    if (pointRoot) {
      sprintf((char *)mBuffer, "%s%d-highSDbound.mod", pointRoot, tomo);
      dumpPointModel(mBoundaries, mNumHighSD, (char *)mBuffer, inHeader);
    }
    if (pitchName) {
      imodSetRefImage(pitchModel, inHeader);
      writeModel(pitchName, pitchModel);
    }
    exit(0);
  }
    
  // Set up the blocking from scratch and allocate more arrays
  scl = mBestScale;
  setupBlocks(mNumBoxes[scl][b3dX], scl, b3dX);
  setupBlocks(mNumBoxes[scl][yInd], scl, yInd);
  numXblocks = mNumBlocks[scl][b3dX];
  numYblocks = mNumBlocks[scl][yInd];
  boundArrSize = 2 * numXblocks * numYblocks;
  allBoundaries = mBoundaries = B3DMALLOC(float, boundArrSize * numTomos);
  mBlockCenters = B3DMALLOC(float, boundArrSize);
  smoothBound = B3DMALLOC(float, boundArrSize);
  if (!mBoundaries || !mBlockCenters || !smoothBound)
    exitError("Allocating array for boundary positions");

  // Loop on tomos to get boundaries and all thickness
  for (tomoSeq = 0; tomoSeq < numTomos; tomoSeq++) {
    tomo = tomoInds[tomoSeq];
    mSDs = &allSDs[tomo * tomoArrSize];
    mMeans = &allMeans[tomo * tomoArrSize];
    mBoundaries = &allBoundaries[tomo * boundArrSize];
    
    // Now find a boundary position in each block; loop on the blocks
    for (iyBlock = 0; iyBlock < numYblocks; iyBlock++) {
      for (ixBlock = 0; ixBlock < numXblocks; ixBlock++) {
        ind = 2 * (iyBlock * numXblocks + ixBlock);

        // Get the range of boxes in the block
        for (ixyz = 0; ixyz <= yInd; ixyz += yInd) {
          boxNum = ixyz ? iyBlock : ixBlock;
          rem = mNumBoxes[scl][ixyz] % mNumBlocks[scl][ixyz];
          starts[ixyz] = boxNum * mNumInBlock[scl][ixyz] + (boxNum <= rem ? boxNum : rem);
          ends[ixyz] = starts[ixyz] + mNumInBlock[scl][ixyz] + (boxNum < rem ? 0 : -1);

          // And get the unbinned center coordinate of the block
          mBlockCenters[ind + B3DMIN(1, ixyz)] = startCoord[ixyz] + mBinning[scl][ixyz] *
            (((starts[ixyz] + ends[ixyz]) * mBoxSpacing[scl][ixyz] + boxSize[scl][ixyz]) /
             2. + boxStart[scl][ixyz]);
        }
        
        // Get the boundaries and scale to unbinned coordinates there too
        if (ifReconArea && !InsideContour(tiltXvert, tiltYvert, 4, mBlockCenters[ind], 
                                          mBlockCenters[ind + 1])) {
          mBoundaries[ind] = mBoundaries[ind + 1] = -1.;
        } else {
          fitColumnBoundaries(starts[0], ends[0], starts[yInd], ends[yInd],
                              &mBoundaries[ind]);
        }
        for (loop = 0; loop < 2; loop++)
          if (mBoundaries[ind + loop] > 0.)
            mBoundaries[ind + loop] = startCoord[mThickInd] + mBinning[scl][mThickInd] *
              (mBoundaries[ind + loop] * mBoxSpacing[scl][mThickInd] + 
               0.5 * boxSize[scl][mThickInd] + boxStart[scl][mThickInd]);
        
        // Add thicknesses to a collection
        if (mBoundaries[ind] >= 0. && mBoundaries[ind + 1] >= 0.)
          mThicknesses[mNumThicknesses++] = mBoundaries[ind + 1] - mBoundaries[ind];
      }
    }
  }
  if (mNumExtraSum)
    mExtraForPitch = mExtraPitchSum / mNumExtraSum;
  if (mDebugOutput)
    printf("Mean extra distance for tomopitch lines = %.1f\n", mExtraForPitch);
  
  if (mNumThicknesses >= mMinNumThickForCheck) {
    rsFastMedian(mThicknesses, mNumThicknesses, mBuffer, &mThickMedian);
    rsFastMADN(mThicknesses, mNumThicknesses, mThickMedian, mBuffer, &mThickMADN);
    if (mDebugOutput)
      PRINT3(mThickMedian, mThickMADN, mNumThicknesses);
  }


  // START OF BIG REMAINING LOOP ON TOMOGRAMS
  for (tomoSeq = 0; tomoSeq < numTomos; tomoSeq++) {
    tomo = tomoInds[tomoSeq];
    mSDs = &allSDs[tomo * tomoArrSize];
    mMeans = &allMeans[tomo * tomoArrSize];
    mBoundaries = &allBoundaries[tomo * boundArrSize];
    
    // Eliminate points that give a bad thickness if possible
    checkBlockThicknesses();

    inHeader = iiuMrcHeader(inUnitBase + tomo, "findsection", 1, 0);
    if (pointRoot) {
      sprintf((char *)mBuffer, "%s%d-colbound.mod", pointRoot, tomo);
      dumpPointModel(mBoundaries, numXblocks * numYblocks, (char *)mBuffer, inHeader);
    }

    // Now for fitting/smoothing, loop on every block, both surfaces
    for (botTop = 0; botTop < 2; botTop++) {
      for (iyBlock = 0; iyBlock < numYblocks; iyBlock++) {
        for (ixBlock = 0; ixBlock < numXblocks; ixBlock++) {
          ixyCen[0] = ixBlock;
          ixyCen[1] = iyBlock;
          ind = 2 * (iyBlock * numXblocks + ixBlock) + botTop;
          smoothBound[ind] = -1.;
          if (mBoundaries[ind] < 0.)
            continue;

          // Try for a square 5 x 5 region first, but if there are not 2 good rows of
          // data at least, drop back to 7 x 2 region
          getFittingRegion(ixyCen, 5, 5, botTop, starts, ends, numBound, xySpacing);
          if (mGoodRowCol[0] < 2 || mGoodRowCol[1] < 2)
            getFittingRegion(ixyCen, 7, 2, botTop, starts, ends, numBound, xySpacing);

          // Set up major and minor axis 
          wgtColIn = -1;
          major = 0;
          if (ends[0] - starts[0] < ends[1] - starts[1])
            major = 1;
          minor = 1 - major;

          // If there are enough points for robust fit, set the variable list
          if (numBound >= 6) {
            buildVariableList(major, minor, numBound, numBound);
            maxZeroWgt = numBound / 6;

            // Load matrix and do the fit; if it works, record weight column
            loadFittingMatrix(ixyCen, starts, ends, botTop, xySpacing, 0., -1, -1,
                              numFit);
            ierr = robustRegress(&mFitMat[0][0], MAX_FIT_DATA, 0, mNumVars, numFit, 1, 
                                 fitSolution, MAX_FIT_COL, &fitConst[0], fitMean, fitSD, 
                                 mFitWork, mKfactor, &numIter, mMaxIter, maxZeroWgt, 
                                 mMaxChange, mMaxOscill);
            if (mDebugOutput > 1)
              printf("block %d %d bt %d fr %d %d %d %d err %d iter %d  c %.1f\n", ixBlock,
                     iyBlock, botTop, starts[0], ends[0], starts[1], ends[1],
                     ierr, numIter, fitConst[0]);
            if (!ierr)
              wgtColIn = mNumVars + 1;
          }

          numFit = numGood = numBound;
          
          if (wgtColIn > 0) {
            
            // If there are weights, count up the number of points still included in
            // numBound, and # with weights above a higher threshold in numGood
            numBound = numGood = 0;
            for (ixyz = 0; ixyz < numFit; ixyz++) {
              if (mFitMat[wgtColIn][ixyz] >= wgtThresh)
                numBound++;
              if (mFitMat[wgtColIn][ixyz] >= goodThresh)
                numGood++;
            }
            if (mDebugOutput > 1)
              printf("total %d  retain %d  good %d\n", numFit, numBound, numGood);
          }

          // Set up variables for final fit with possible weighting, with requirements on
          // number to fit as well as number with higher weighting; load and fit
          smoothBound[ind] = mBoundaries[ind];
          if (numBound >= 5) {
            buildVariableList(major, minor, numBound, numGood);
            wgtColOut = mNumVars + 1;
            loadFittingMatrix(ixyCen, starts, ends, botTop, xySpacing, wgtThresh,
                              wgtColIn, wgtColOut, numFit);
            ierr = multRegress(&mFitMat[0][0], MAX_FIT_DATA, 0, mNumVars, numFit, 1,
                               wgtColIn, fitSolution, MAX_FIT_COL, &fitConst[0], fitMean,
                               fitSD, mFitWork);
            if (mDebugOutput > 1)
              printf("block %d %d bt %d fr %d %d %d %d err %d xm %.1f %.1f %.1f b %f %f "
                     "c %.1f\n", ixBlock, iyBlock, botTop, starts[0], ends[0], starts[1],
                     ends[1], ierr, fitMean[0], fitMean[1], fitMean[2], fitSolution[0], 
                     fitSolution[1], fitConst[0]);
            if (!ierr)
              smoothBound[ind] = fitConst[0];
          }

        }
      }
    }

    if (pointRoot) {
      sprintf((char *)mBuffer, "%s%d-smooth.mod", pointRoot, tomo);
      dumpPointModel(smoothBound, numXblocks * numYblocks, (char *)mBuffer, inHeader);
    }

    if (surfaceName)
      makeSurfaceModel(smoothBound, surfaceName, inHeader, mNxyz);

    // For outputting tomopitch model, first assign the header for first tomo
    if (pitchName) {
      if (!tomo)
        imodSetRefImage(pitchModel, inHeader);

      // Add to model for sample tomograms
      if (numTomos > 1) {
        if (!addToPitchModel(pitchModel, smoothBound, 0, numYblocks - 1, tomo + 1))
          numPitchPairs++;
      } else {

        // Otherwise figure out number of blocks in sample extent
        if (sampleExtent) {
          size = (mBlockCenters[2 * numXblocks * (numYblocks - 1) + 1] - mBlockCenters[1])
            / (float)B3DMAX(numYblocks - 1, 1);
          size = B3DMAX(1, B3DNINT((float)sampleExtent / size));
        } else {
          size = 1;
        }

        // Limit number of samples if needed
        if (numYblocks / size < numSamples) {
          numSamples = numYblocks / size;
          printf("WARNING: With a block size of %d, there can be only %d samples\n",
                 mScanBlockSize, numSamples);
        }

        // Indent by one block if there is enough extra stuff
        ind = 0;
        if (size * numSamples < (numYblocks - 2 * size) / 2)
          ind = size;

        // Get spacing between samples and # requiring spacing + 1
        samBlockSpace = (numYblocks - 2 * ind - size) / B3DMAX(1, numSamples - 1);
        rem = (numYblocks - 2 * ind - size) % B3DMAX(1, numSamples - 1);

        // Loop on the samples
        for (loop = 0; loop < numSamples; loop++) {
          if (!addToPitchModel(pitchModel, smoothBound, ind, ind + size - 1, 0))
            numPitchPairs++;
          ind += samBlockSpace + (loop < rem ? 1 : 0);
        }
      }
    }

    iiuClose(inUnitBase + tomo);
  }

  // Finish tomopitch model
  if (pitchName) {
    
    writeModel(pitchName, pitchModel);
    if (numPitchPairs < 2)
      exitError("Only one pair of lines was placed in the model for tomopitch");
    if (numPitchPairs < numTomos || (numTomos == 1 && numPitchPairs < numSamples))
      printf("WARNING: %s - A pair of lines was found for only %d of the %d samples %s\n",
             progname, numPitchPairs, numTomos > 1 ? numTomos : numSamples,
             numTomos > 1 ? "tomograms" : "positions");
  }
  
  // For single tomogram, get the median on each surface and output Z limits
  if (numTomos == 1) {
    for (loop = 0; loop < 2; loop++) {
      patchLim[loop] = -2;
      numStat = 0;
      sign = 2 * loop - 1;
      for (ind = 0; ind < numXblocks * numYblocks; ind++) {
        if (smoothBound[2 * ind + loop] >= 0.)
          mBuffer[numStat++] = smoothBound[2 * ind + loop];
      }
      if (numStat > 3) {
        rsSortFloats(mBuffer, numStat);
        rsMedianOfSorted(mBuffer, numStat, &insideMed);
        
        // Determine the integer slice that the median plus extra amount occurs in
        patchLim[loop] = B3DNINT(insideMed + sign * mExtraForPitch - 0.5);
        B3DCLAMP(patchLim[loop], 0, mNxyz[mThickInd] - 1);

        // Determine the same for absolute limits of surfaces
        if (loop)
          surfaceLim[loop] = B3DNINT(mBuffer[numStat - 1] + mExtraForPitch - 0.5);
        else
          surfaceLim[loop] = B3DNINT(mBuffer[0] - mExtraForPitch - 0.5);
        B3DCLAMP(surfaceLim[loop], 0, mNxyz[mThickInd] - 1);

        // Do combination of the more extreme of a "high" percentile limit, and a
        // very low percentile limit with some number pixels allowed to be outside
        ratio = combineLowPctl * (1 - loop) + loop * (1. - combineLowPctl);
        rsPercentileOfSorted(mBuffer, numStat, ratio, &insideMed);
        ratio = combineHighPctl * (1 - loop) + loop * (1. - combineHighPctl);
        rsPercentileOfSorted(mBuffer, numStat, ratio, &yy);
        if (loop)
          combineLim[loop] = B3DNINT(B3DMAX(yy, insideMed - combineOutsidePix) +
                                     mExtraForPitch - 0.5);
        else
          combineLim[loop] = B3DNINT(B3DMIN(yy, insideMed + combineOutsidePix) -
                                     mExtraForPitch - 0.5);
      }
    }

    // Output results if any
    if (patchLim[0] >= 0  && patchLim[1] >= 0) {
      invertYifFlipped(patchLim[0], patchLim[1], mNxyz);
      printf("Median Z values of surfaces, numbered from 1, are: %d  %d\n", 
             patchLim[0] + 1, patchLim[1] + 1);
      invertYifFlipped(combineLim[0], combineLim[1], mNxyz);
      printf("Z limits for autopatchfit combine, numbered from 1, are: %d  %d\n", 
             combineLim[0] + 1, combineLim[1] + 1);
      invertYifFlipped(surfaceLim[0], surfaceLim[1], mNxyz);
      printf("Absolute limits of surfaces, numbered from 1, are: %d  %d\n", 
             surfaceLim[0] + 1, surfaceLim[1] + 1);
    } else
      printf("Too few surface points to determine summary Z values for surface\n");
  }

  exit(0);
}

/*
 * Write out a model of points, either the raw column boundaries or the smoothed 
 * boundaries
 */
void FindSect::dumpPointModel(float *boundaries, int numPts,
                              const char *filename, MrcHeader *inHeader)
{
  Imod *imod = imodNew();
  if (!imod || imodNewObject(imod) || imodNewContour(imod) || imodNewObject(imod) || 
      imodNewContour(imod))
    exitError("Setting up model");
  bool flipped = mThickInd == 1;
  int pt, ind;
  imodSetRefImage(imod, inHeader);
  for (pt = 0; pt < numPts; pt++) {
    for (ind = 0; ind < 2; ind++)
      if (boundaries[2 * pt + ind] >= 0 && !imodPointAppendXYZ
          (imod->obj[ind].cont, mBlockCenters[2 * pt],
           flipped ? boundaries[2 * pt + ind] : mBlockCenters[2 * pt + 1],
           flipped ? mBlockCenters[2 * pt + 1] : boundaries[2 * pt + ind]))
        exitError("Adding point to model");
  }
  imod->obj[0].flags = IMOD_OBJFLAG_SCAT | IMOD_OBJFLAG_PNT_ON_SEC;
  imod->obj[0].pdrawsize = 4;
  imod->obj[1].flags = IMOD_OBJFLAG_SCAT | IMOD_OBJFLAG_PNT_ON_SEC;
  imod->obj[1].pdrawsize = 4;
  writeModel(filename, imod);
  imodDelete(imod);
}

/*
 * Put out a model of the smoothed boundaries as open contours along the surface,
 * usable for flattenwarp
 */
void FindSect::makeSurfaceModel(float *boundaries, const char *filename,
                                MrcHeader *inHeader, int *nxyz)
{
  Imod *imod = imodNew();
  Icont *cont;
  if (!imod || imodNewObject(imod))
    exitError("Setting up model");
  bool flipped = mThickInd == 1;
  int numXblocks = mNumBlocks[mBestScale][b3dX];
  int numYblocks = mNumBlocks[mBestScale][3 - mThickInd];
  int pt, ind, iy, ix;
  float xscale, yscale, zscale;
  float zRound;
  imodSetRefImage(imod, inHeader);
  mrc_get_scale(inHeader, &xscale, &yscale, &zscale);
  imod->xmax = mNxyz[0];
  imod->ymax = mNxyz[1];
  imod->zmax = mNxyz[2];
  
  // Loop on levels in Y, and on bottom and top surfaces
  for (iy = 0; iy < numYblocks; iy++) {
    for (ind = 0; ind < 2; ind++) {
      cont = NULL;

      // Loop across
      for (ix = 0; ix < numXblocks; ix++) {
        pt = iy * numXblocks + ix;
        if (boundaries[2 * pt + ind] >= 0) {

          // Add contour only when the first point is found
          if (!cont) {
            if (imodNewContour(imod) || (cont = imodContourGet(imod)) == NULL)
              exitError("Adding contour to model");
          }
          zRound = B3DNINT(mBlockCenters[2 * pt + 1]);
          if (!imodPointAppendXYZ(cont, mBlockCenters[2 * pt],
                                  flipped ? boundaries[2 * pt + ind] : zRound,
                                  flipped ? zRound : boundaries[2 * pt + ind]))
            exitError("Adding point to model");
        }
      }
    }
  }
  imod->obj[0].flags = IMOD_OBJFLAG_OPEN;
  imod->obj[0].symsize = 7;
  imod->obj[0].symbol = IOBJ_SYM_CIRCLE;
  if (xscale && xscale != 1.0f) {
    imod->pixsize = xscale / 10.f;
    imod->units = IMOD_UNIT_NM;
  }
  writeModel(filename, imod);
  imodDelete(imod);
}

void FindSect::writeModel(const char *filename, Imod *imod)
{
  imodBackupFile(filename);
  FILE *fp = fopen(filename, "wb");
  if (!fp || imodWrite(imod, fp))
    exitError("Opening or writing model %s", filename);
  fclose(fp);
}

/*
 * Fit two lines separately or together with the same slope to the points on bottom and
 * top surfaces in the given range of blocks in Y.  Use robust fitting if there are
 * enough points.  Returns 1 if there are inadequate points or a fitting error.
 */
int FindSect::addToPitchModel(Imod *imod, float *boundaries, int yStart, int yEnd, 
                               int time)
{
  Icont *cont;
  float fitSD[MAX_FIT_COL], fitMean[MAX_FIT_COL], fitSolution[MAX_FIT_COL];
  bool flipped = mThickInd == 1;
  int numXblocks = mNumBlocks[mBestScale][b3dX];
  int pt, ind, iy, ix, err, numPts = 0, numBot, numFit, maxZeroWgt, indFit, numIter;
  int maxPts = 2 * numXblocks * (yEnd + 1 - yStart);
  float zRound, yVal, yShift[2];
  float *xfit = mBuffer;
  float *yfit = &mBuffer[maxPts];
  float *zfit = &mBuffer[2 * maxPts];
  float slopes[2], intercepts[2], xmin[2], xmax[2];
  int minOnSide = mFitPitchSeparately ? 3 : 2;
  
  // Loop on bottom and top surfaces
  for (ind = 0; ind < 2; ind++) {
    xmin[ind] = 1.e37;
    xmax[ind] = -1.e37;
    numBot = numPts;

    // Load the points into arrays
    for (iy = yStart; iy <= yEnd; iy++) {
      for (ix = 0; ix < numXblocks; ix++) {
        pt = iy * numXblocks + ix;
        if (boundaries[2 * pt + ind] >= 0) {
          xfit[numPts] = mBlockCenters[2 * pt];
          xmin[ind] = B3DMIN(xmin[ind], xfit[numPts]);
          xmax[ind] = B3DMAX(xmax[ind], xfit[numPts]);
          yfit[numPts] = boundaries[2 * pt + ind];
          zfit[numPts] = ind;
          numPts++;
        }
      }
    }
  }
  if ((mFitPitchSeparately && numPts < 5) || numBot < minOnSide || 
      numPts - numBot < minOnSide)
    return 1;

  // Do fit to each surface separately or to both at once, set up slope/intercept for
  // top in the latter case
  err = 1;
  if (mFitPitchSeparately) {
    numFit = numBot;
    indFit = 0;
    for (ind = 0; ind < 2; ind++) {

      // Use robust fit if there are enough points
      if (numFit >= mMinForRobustPitch && numFit <= MAX_FIT_DATA) {
        for (ix = 0; ix < numFit; ix++) {
          mFitMat[0][ix] = xfit[ix + indFit];
          mFitMat[1][ix] = yfit[ix + indFit];
        }
        maxZeroWgt = numFit / 6;
        err = robustRegress(&mFitMat[0][0], MAX_FIT_DATA, 0, 1, numFit, 1, 
                            fitSolution, MAX_FIT_COL, &intercepts[ind], fitMean, fitSD, 
                            mFitWork, mKfactor, &numIter, mMaxIter, maxZeroWgt, 
                            mMaxChange, mMaxOscill);
        slopes[ind] = fitSolution[0];
      } 

      // Otherwise, or if there was an error in the robust fit, do standard fit
      if (err) {
        lsFit(&xfit[indFit], &yfit[indFit], numFit, &slopes[ind], &intercepts[ind],
              &zRound);
      }
      numFit = numPts - numBot;
      indFit = numBot;
    }
  } else {

    // Require 2 more points for robust fit
    if (numPts >= mMinForRobustPitch + 2 && numPts <= MAX_FIT_DATA) {
      for (ix = 0; ix < numPts; ix++) {
        mFitMat[0][ix] = xfit[ix];
        mFitMat[1][ix] = zfit[ix];
        mFitMat[2][ix] = yfit[ix];
      }
      maxZeroWgt = numPts / 6;
      err = robustRegress(&mFitMat[0][0], MAX_FIT_DATA, 0, 2, numPts, 1, 
                         fitSolution, MAX_FIT_COL, &intercepts[0], fitMean, fitSD, 
                         mFitWork, mKfactor, &numIter, mMaxIter, maxZeroWgt, 
                         mMaxChange, mMaxOscill);
      slopes[0] = fitSolution[0];
      zRound = fitSolution[1];
    }

    // Do regular fit if robust not tried or failed
    if (err) {
      lsFit2(xfit, zfit, yfit, numPts, &slopes[0], &zRound, &intercepts[0]);
    }
    slopes[1] = slopes[0];
    intercepts[1] = intercepts[0] + zRound;
  }

  // Find a shift equal to maximum residual on each side
  yShift[0] = yShift[1] = 0.;
  for (pt = 0; pt < numPts; pt++) {
    ind = pt < numBot ? 0 : 1;
    yVal = yfit[pt] - (xfit[pt] * slopes[ind] + intercepts[ind]);
    if ((ind && yVal > yShift[ind])|| (!ind && yVal < yShift[ind]))
      yShift[ind] = yVal;
  }

  // Make the lines
  for (ind = 0; ind < 2; ind++) {
    yShift[ind] += (2 * ind - 1) * mExtraForPitch;
    if (imodNewContour(imod) || (cont = imodContourGet(imod)) == NULL)
      exitError("Adding contour to model");
    zRound = B3DNINT((mBlockCenters[2 * yStart * numXblocks + 1] + 
                      mBlockCenters[2 * yEnd * numXblocks + 1]) / 2.);
    yVal = xmin[ind] * slopes[ind] + intercepts[ind] + yShift[ind];
    if (!imodPointAppendXYZ(cont, xmin[ind], flipped ? yVal : zRound,
                            flipped ? zRound : yVal))
      exitError("Adding point to model");
    yVal = xmax[ind] * slopes[ind] + intercepts[ind] + yShift[ind];
    if (!imodPointAppendXYZ(cont, xmax[ind], flipped ? yVal : zRound,
                            flipped ? zRound : yVal))
      exitError("Adding point to model");
    cont->time = time;
  }
  return 0;
}

/*
 * Determine number of non-overlapping blocks and number of analyzed points in block for
 * a given scale index and axis
 */
void FindSect::setupBlocks(int numBoxes, int sclInd, int ixyz)
{
  mNumInBlock[sclInd][ixyz] = mScanBlockSize / 
    (mBoxSpacing[sclInd][ixyz] * mBinning[sclInd][ixyz]);
  B3DCLAMP(mNumInBlock[sclInd][ixyz], 1, numBoxes);
  mNumBlocks[sclInd][ixyz] = numBoxes / mNumInBlock[sclInd][ixyz];
  mNumInBlock[sclInd][ixyz] = numBoxes / mNumBlocks[sclInd][ixyz];
}

/*
 * For the range of boxes indicated, and the scaling binInd, add the SD values to 
 * the list in mBuffer and accumulate the sum of means in meanSum.
 */
void FindSect::addBoxesToSample(int startX, int endX, int startY, int endY, int startZ,
                                int endZ, int binInd, int &numStat, double &meanSum)
{
  int boxInd;
  for (int izBox = startZ; izBox <= endZ; izBox++) {
    for (int iyBox = startY; iyBox <= endY; iyBox++) {
      for (int ixBox = startX; ixBox <= endX; ixBox++) {
        boxInd = mStatStartInds[binInd] + (izBox * mNumBoxes[binInd][b3dY] + iyBox) * 
          mNumBoxes[binInd][b3dX] + ixBox;
        mBuffer[numStat++] = mSDs[boxInd];
        meanSum += mMeans[boxInd];
      }
    }
  }  
}

/*
 * Adjust a range in Y to apply to a rotated volume if the data are flipped; this
 * gives consistency between results from original volumes in flipped orientation and
 * volume rotated with rotx.
 */
void FindSect::invertYifFlipped(int &start, int &end, int *dims)
{
  if (mThickInd == b3dY) {
    int tmp = start;
    start = dims[1] - 1 - end;
    end = dims[1] - 1 - tmp;
  }
}

/*
 * For a column of boxes through the thickness dimension with the given box extents in X 
 * and Y, find the midpoint of the region with structure by looking from each edge for the
 * Z value where the SD goes above criterion (returned in izInside).  The middle of those
 * values is returned in izMid, and the median of SD values between that range is 
 * returned in insideMedian.  Return value is 1 if the maximum median in the column is
 * not different enough from edge; 2 if there are too few median SD's above a criterion;
 * or 3 if the inside positions end up being crossed.
 */
int FindSect::findColumnMidpoint(int startX, int endX, int startY, int endY, int binInd,
                                 float *buffer, int *izInside, int &izMid,
                                 float &insideMedian)
{
  float medSD, edgeDiff, edgeCrit, sdMax = -1.;
  int numXYbox = (endY + 1 - startY) * (endX + 1 - startX);
  int boxInd, loop, dir, yStride = 1, zStride = 1;
  int ixBox, iyBox, izBox, zStart, zEnd, numAbove;

  if (mThickInd == b3dZ)
    zStride = mNumBoxes[binInd][b3dY];
  else
    yStride = mNumBoxes[binInd][b3dY];

  // First find a maximum in the column, taking the median of values across each plane
  for (izBox = 0; izBox < mNumBoxes[binInd][mThickInd]; izBox++) {
    dir = 0;
    for (iyBox = startY; iyBox <= endY; iyBox++) {
      for (ixBox = startX; ixBox <= endX; ixBox++) {
        boxInd = mStatStartInds[binInd] + (izBox * zStride + iyBox * yStride) * 
          mNumBoxes[binInd][b3dX] + ixBox;
        buffer[dir++] = mSDs[boxInd];
      }
    }
    rsFastMedianInPlace(buffer, numXYbox, &medSD);
    sdMax = B3DMAX(sdMax, medSD);
    buffer[numXYbox + izBox] = medSD;
    if (mColMedians)
      mColMedians[izBox] = medSD;
  }

  // If difference from edge is below criterion, return error
  edgeDiff = sdMax - mEdgeMedians[binInd];
  if (edgeDiff < mColMaxEdgeDiffCrit * mEdgeMADNs[binInd])
    return 1;
  
  // Criterion is maximum of a number of MADNs above edge and a fraction of max-edge diff
  edgeCrit = mEdgeMedians[binInd] + B3DMAX(edgeDiff * mFracColMaxEdgeDiff, 
                                           mCritEdgeMADN * mEdgeMADNs[binInd]);

  // From each direction, find an edge above the criterion
  zStart = 0;
  zEnd = mNumBoxes[binInd][mThickInd];
  for (loop = 0; loop < 2; loop++) {
    dir = 1 - 2 * loop;
    numAbove = 0;
    for (izBox = zStart; izBox != zEnd; izBox += dir) {
      if (buffer[numXYbox + izBox] >= edgeCrit) {
        numAbove++;
        if (numAbove >= mNumHighInsideCrit) {
          izInside[loop] = izBox;
          break;
        }
      }
    }
    if (numAbove < mNumHighInsideCrit)
      return 2;
    zStart = mNumBoxes[binInd][mThickInd] - 1;
    zEnd = -1;
  }

  //printf("sdMax  %f  edgeCrit %f  %d  %d\n", sdMax, edgeCrit, izInside[0], izInside[1]);
  if (izInside[0] > izInside[1])
    return 3;
  izMid = (izInside[0] + izInside[1]) / 2;
  rsFastMedianInPlace(&buffer[numXYbox + izInside[0]], izInside[1] + 1 - izInside[0], 
                      &insideMedian);
  return 0;
}

/*
 * Finds boundaries of a column defined by the range of boxes in X and Y, looking outward
 * from middle, and fits a line to the falling phase to find the boundary at a set level.
 */
void FindSect::fitColumnBoundaries(int startX, int endX, int startY, int endY,
                                    float *boundary)
{
  int boxInd, loop, ind, dir, yStride = 1, zStride = 1, numExtra = 0;
  int izMid, izInside[2], numInCol[2];
  float insideMed, fallCrit, lowFitCrit, highFitCrit, boundaryLevel, pitchLevel;
  float lastSD, slope, intercept, ro, extraMed;
  int ixBox, iyBox, izBox, iz, izAdd, fitStart, numFit;
  int scl = mBestScale;
  int zRange = mNumBoxes[scl][mThickInd];
  int maxInCol = (endX + 1 - startX) * (endY + 1 - startY);
  float *xfit = mBuffer + 4 * maxInCol;
  float *yfit = xfit + zRange;

  if (mThickInd == b3dZ)
    zStride = mNumBoxes[scl][b3dY];
  else
    yStride = mNumBoxes[scl][b3dY];
  boundary[0] = boundary[1] = -1.;

  // Find midpoint and inside median, reject if it is too low
  if (findColumnMidpoint(startX, endX, startY, endY, scl, mBuffer, izInside, izMid,
                         insideMed))
    return;
  if (insideMed < mCenMedians[scl] - mColumnToCenMADNCrit * mCenMADNs[scl]) {
    if (mDebugOutput > 1)
      printf("skipping %d %d low median %f\n", startX, startY, insideMed);
    return;
  }

  // Set various criteria
  fallCrit = mEdgeMedians[scl] + mMaxFalloffFrac * (insideMed - mEdgeMedians[scl]);
  lowFitCrit = mEdgeMedians[scl] + mLowFitFrac * (insideMed - mEdgeMedians[scl]);
  highFitCrit = mEdgeMedians[scl] + mHighFitFrac * (insideMed - mEdgeMedians[scl]);
  boundaryLevel = mEdgeMedians[scl] + mBoundaryFrac * 
    (insideMed - mEdgeMedians[scl]);
  pitchLevel = mEdgeMedians[scl] + mPitchBoundaryFrac *
    (insideMed - mEdgeMedians[scl]);

  numInCol[0] = numInCol[1] = 0;
  for (iyBox = startY; iyBox <= endY; iyBox++) {
    for (ixBox = startX; ixBox <= endX; ixBox++) {
      for (loop = 0; loop < 2; loop++) {
        dir = 2 * loop - 1;
        numFit = 0;
        fitStart = -1;
        for (izBox = izMid; izBox >= 0 && izBox < zRange; izBox += dir) {
          
          boxInd = mStatStartInds[scl] + (izBox * zStride + iyBox * yStride) * 
            mNumBoxes[scl][b3dX] + ixBox;
          if (mSDs[boxInd] < fallCrit) {

            // Going below the fall criterion triggers various checks
            // First, if the starting point is below it, skip this box column
            if (izBox == izMid)
              break;

            // If this point is below the low fit criterion or is higher than the
            // the last, start fit on previous point, unless it is above the high
            // criterion; in which case start on this one
            if (mSDs[boxInd] < lowFitCrit || mSDs[boxInd] > lastSD) {
              fitStart = izBox - dir;
              if (lastSD > highFitCrit)
                fitStart = izBox;
              
              // But if we are at end of range, start fit on this point
            } else if (izBox == 0 || izBox == zRange - 1) {
              fitStart = izBox;
            }
            
            if (fitStart >= 0) {
              
              // Add points to fit arrays until it goes above the high crit
              for (iz = fitStart; iz != izMid; iz -= dir) {
                ind = mStatStartInds[scl] + (iz * zStride + iyBox * yStride) *
                  mNumBoxes[scl][b3dX] + ixBox;
                if (mSDs[ind] > highFitCrit)
                  break;
                xfit[numFit] = iz;
                yfit[numFit++] = mSDs[ind];
              }
              
              // If there is only one point, add one that is out of the fit range
              // in the other direction if possible
              if (numFit == 1) {
                iz = B3DNINT(xfit[0]);
                izAdd = -1;
                if (yfit[0] < 0.5 * (lowFitCrit + highFitCrit)) {
                  if (iz - dir != izMid)
                    izAdd = iz - dir;
                } else {
                  if (iz + dir >= 0 && iz + dir < zRange)
                    izAdd = iz + dir;
                }
                if (izAdd >= 0) {
                  ind = mStatStartInds[scl] + 
                    (izAdd * zStride + iyBox * yStride) * mNumBoxes[scl][b3dX] + ixBox;
                  xfit[numFit] = izAdd;
                  yfit[numFit++] = mSDs[ind];
                }
              }
              
              // Do the fit if at least 2 points and get Z value at boundary level
              if (numFit >= 2) {
                lsFit(xfit, yfit, numFit, &slope, &intercept, &ro);
                mBuffer[loop * maxInCol + numInCol[loop]++] = 
                  (boundaryLevel - intercept) / slope;
                mBuffer[2 * maxInCol + numExtra++] = 
                  fabs((double)(boundaryLevel - pitchLevel) / slope);
              }
              break;
            }    
          }
          lastSD = mSDs[boxInd];
        }
      }
    }
  }

  // Get median of boundary values if there are enough of them
  for (loop = 0; loop < 2; loop++)
    if (numInCol[loop] >= B3DNINT(B3DMAX(1., mMinFracBoundsInCol * maxInCol)))
      rsFastMedianInPlace(&mBuffer[loop * maxInCol], numInCol[loop], &boundary[loop]);

  // Get the median of extra Z values if there are enough, and add to sum
  if (numExtra >= B3DNINT(B3DMAX(1., mMinFracBoundsInCol * 2. * maxInCol))) {
    mNumExtraSum++;
    rsFastMedianInPlace(&mBuffer[2 * maxInCol], numExtra, &extraMed);
    mExtraPitchSum += extraMed * mBinning[scl][mThickInd];
  }
}

/*
 * Tries to identify blocks that have insufficient thickness and eliminate one or
 * both boundaries
 */
void FindSect::checkBlockThicknesses()
{
  int scl = mBestScale;
  int yInd = 3 - mThickInd;
  int numXblocks = mNumBlocks[scl][b3dX];
  int numYblocks = mNumBlocks[scl][yInd];
  int ixBlock, iyBlock, ind, ixyCen[2], starts[2], ends[2], botTop, numBound, ix, iy, jj;
  float thickMedian, thickMADN, xySpacing, meanBound, meanDiff[2], tmp1, tmp2, numMADNs;
  int minBound = 10000;

  if (mNumThicknesses <  mMinNumThickForCheck)
    return;

  // Look for outlier thicknesses
  for (iyBlock = 0; iyBlock < numYblocks; iyBlock++) {
    for (ixBlock = 0; ixBlock < numXblocks; ixBlock++) {
      ind = 2 * (iyBlock * numXblocks + ixBlock);
      if (mBoundaries[ind] < 0. || mBoundaries[ind + 1] < 0. ||
          (mBoundaries[ind + 1] - mBoundaries[ind]) / mThickMedian > mTooThinCrit)
        continue;

      if (mDebugOutput)
        printf("block too thin %d %d  %.1f\n", ixBlock, iyBlock, 
               mBoundaries[ind + 1] - mBoundaries[ind]);
      ixyCen[0] = ixBlock;
      ixyCen[1] = iyBlock;

      // Find a sampling region and extract boundaries for ones not involved in an
      // outlier thickness
      for (botTop = 0; botTop < 2; botTop++) {
        getFittingRegion(ixyCen, 5, 5, botTop, starts, ends, numBound, xySpacing);
        if (mGoodRowCol[0] < 2 || mGoodRowCol[1] < 2)
          getFittingRegion(ixyCen, 7, 2, botTop, starts, ends, numBound, xySpacing);
        numBound = 0;
        for (ix = starts[0]; ix <= ends[0]; ix++) {
          for (iy = starts[1]; iy <= ends[1]; iy++) {
            jj = 2 * (iy * numXblocks + ix);
            if (mBoundaries[jj + botTop] >= 0. && 
                (mBoundaries[jj + 1 - botTop] < 0. ||
                 (mBoundaries[jj + 1] - mBoundaries[jj]) / mThickMedian > mTooThinCrit))
              mBuffer[numBound++] = mBoundaries[jj + botTop];
          }
        }
        ACCUM_MIN(minBound, numBound);
        meanDiff[botTop] = -1.;
        if (numBound > 2) {
          if (numBound > 5)
            rsFastMedianInPlace(mBuffer, numBound, &meanBound);
          else
            avgSD(mBuffer, numBound, &meanBound, &tmp1, &tmp2);
          meanDiff[botTop] = fabs(mBoundaries[ind + botTop] - meanBound);
        }
      }

      // Evaluate the difference from the local mean for each boundary
      // Drop it if it is sufficiently larger than the difference for the other boundary
      // or if it is bigger than a fraction of the median thickness
      if (meanDiff[0] >= 0. && meanDiff[1] >= 0.) {
        if (meanDiff[0] > mFartherFromMeanCrit * meanDiff[1] ||
            meanDiff[0] > mMeanDiffThickFrac * mThickMedian) {
          if (mDebugOutput)
            printf("Dropping bottom %.1f diffs %.1f %.1f\n", mBoundaries[ind],
                   meanDiff[0], meanDiff[1]);
          mBoundaries[ind] = -1.;
        } 
        if (meanDiff[1] > mFartherFromMeanCrit * meanDiff[0] ||
            meanDiff[1] > mMeanDiffThickFrac * mThickMedian) {
          if (mDebugOutput)
            printf("Dropping top %.1f diffs %.1f %.1f\n", mBoundaries[ind + 1],
                   meanDiff[0], meanDiff[1]);
          mBoundaries[ind + 1] = -1.;
        }

      }

      // If that didn't eliminate either, apply an outlier criterion to the deviation
      // from median thickness as long as there are at least two boundary points left
      numMADNs = (mThickMedian - (mBoundaries[ind + 1] - mBoundaries[ind])) / mThickMADN;
      if (mBoundaries[ind] >= 0. && mBoundaries[ind + 1] >= 0. && minBound >= 2 &&
          numMADNs > mCritThickMADN) {
        if (mDebugOutput)
          printf("Dropping both: %.1f MADNs from median\n", numMADNs);
        mBoundaries[ind] = -1.;
        mBoundaries[ind + 1] = -1.;
      }
    }
  }
}
        
/*
 * For a block whose X, Y center is in blockCen, finds a region to fit with desired size 
 * extentNum in the major direction, and maximum extent maxDepth in the other, for the
 * surface in botTop.   Returns starting and ending blocks in blockStart and blockEnd,
 * num of positions with data in numPts, and the spacing between blocks in xySpacing.
 */
void FindSect::getFittingRegion(int blockCen[], int extentNum, int maxDepth, int botTop,
                                int *blockStart, int *blockEnd, int &numPts, 
                                float &xySpacing)
{
  int scl = mBestScale;
  int yInd = 3 - mThickInd;
  int numXblocks = mNumBlocks[scl][b3dX];
  int numYblocks = mNumBlocks[scl][yInd];
  int stride[2], numBlk[2];
  int extentAxis = -1;
  float distance, spacing[2] = {0., 0.};
  bool shifted;
  int ixy, ind, other, inner, outer, end, start, depAxis, numGood;
  stride[0] = 2;
  stride[1] = 2 * numXblocks;

  // Find direction that defines the extent
  for (ixy = 0; ixy < 2; ixy++) {
    start = B3DMAX(0, blockCen[ixy] - 1);
    end = B3DMIN(mNumBlocks[scl][ixy * yInd] - 1, start + 2);
    start = B3DMAX(0, end - 2);
    if (start == end) {
      extentAxis = 1 - ixy;
      continue;
    }
    ind = blockCen[1 - ixy] * stride[1 - ixy] + ixy;
    spacing[ixy] = (mBlockCenters[end * stride[ixy] + ind] - 
                    mBlockCenters[start * stride[ixy] + ind]) / (end - start);
  }

  // If either axis is feasible, take the one with the bigger number of blocks if either
  // is limited, or the one with smaller spacing 
  if (extentAxis < 0) {
    if (numXblocks < extentNum || numYblocks < extentNum)
      extentAxis = numXblocks < numYblocks ? 1 : 0;
    else
      extentAxis = B3DCHOICE(spacing[1] < 0.9 * spacing[0], 1, 0);
  }
  depAxis = 1 - extentAxis;
  numBlk[extentAxis] = extentNum;
  xySpacing = B3DMAX(1., spacing[extentAxis]);

  // Get the number in the other direction: shoot for equal extent, limit by maximum
  if (!spacing[depAxis]) {
    numBlk[depAxis] = 1;
  } else {
    distance = spacing[extentAxis] * (extentNum - 1);
    numBlk[depAxis] = B3DNINT(1. + distance / spacing[depAxis]);
    numBlk[depAxis] = B3DMIN(numBlk[depAxis], maxDepth);
  }

  // Get the start and end in each direction
  for (ixy = 0; ixy < 2; ixy++) {
    blockStart[ixy] = B3DMAX(0, blockCen[ixy] - (numBlk[ixy] / 2));
    blockEnd[ixy] = B3DMIN(mNumBlocks[scl][ixy * yInd] - 1,
                           blockStart[ixy] + numBlk[ixy] - 1);
    blockStart[ixy] = B3DMAX(0, blockEnd[ixy] + 1 - numBlk[ixy]);
    numBlk[ixy] = blockEnd[ixy] + 1 - blockStart[ixy];
  }

  // Try to slide region if there are empty rows on one side
  for (ixy = 0; ixy < 2; ixy++) {
    shifted = false;
    other = 1 - ixy;
    while (!hasBoundaries(ixy, blockStart[ixy], blockStart[other], blockEnd[other], 
                          botTop) && 
           hasBoundaries(ixy, blockEnd[ixy] + 1, blockStart[other], blockEnd[other],
                         botTop) && blockStart[ixy] < blockCen[ixy]) {
      shifted = true;
      blockStart[ixy]++;
      blockEnd[ixy]++;
      if (mDebugOutput > 1)
        printf("Shifted + on axis %d\n", ixy);
    }
    while (!shifted && !hasBoundaries(ixy, blockEnd[ixy], blockStart[other],
                                      blockEnd[other], botTop) &&
           hasBoundaries(ixy, blockStart[ixy] - 1, blockStart[other], blockEnd[other],
                         botTop) && blockEnd[ixy] > blockCen[ixy]) {
      blockStart[ixy]--;
      blockEnd[ixy]--;
      if (mDebugOutput > 1)
        printf("Shifted - on axis %d\n", ixy);
    }
  }

  // Count the good rows and columns and total points
  for (ixy = 0; ixy < 2; ixy++) {
    other = 1 - ixy;

    // For counting in a direction, this determines number good in other direction
    // Loop on other direction
    mGoodRowCol[other] = 0;
    numPts = 0;
    for (outer = blockStart[other]; outer <= blockEnd[other]; outer++) {
      numGood = 0;

      // Loop on main direction and count ones with data
      for (inner = blockStart[ixy]; inner <= blockEnd[ixy]; inner++)
        if (mBoundaries[inner * stride[ixy] + outer * stride[other] + botTop] >= 0)
          numGood++;

      // The line is good if it has at least half of the blocks with data
      if (numGood >= B3DMAX(0.99, numBlk[ixy] / 2. - 0.1))
        mGoodRowCol[other]++;
      numPts += numGood;
    }
  }
}

/*
 * Tests whether there are any boundaries along one axis, at the given row or column,
 * between start and end on the other axis, on surface given by botTop
 */
bool FindSect::hasBoundaries(int axis, int rowCol, int start, int end, int botTop)
{
  int stride[2] = {2, 2};
  stride[1] = 2 * mNumBlocks[mBestScale][b3dX];
  if (rowCol < 0 || rowCol >= mNumBlocks[mBestScale][axis * (3 - mThickInd)])
    return false;
  for (int ind = start; ind <= end; ind++)
    if (mBoundaries[rowCol * stride[axis] + ind * stride[1 - axis] + botTop] >= 0)
      return true;
  return false;
}

/*
 *
 */
void FindSect::loadFittingMatrix(int *ixyCen, int *starts, int *ends, int botTop,
                                 float xySpacing, float wgtThresh, int wgtColIn,
                                 int wgtColOut, int &numFit)
{
  int ixBlock, iyBlock, cenInd, blkInd, ixy, var, indWgt = -1;
  numFit = 0;
  cenInd = 2 * (ixyCen[0] + ixyCen[1] * mNumBlocks[mBestScale][b3dX]);
  for (iyBlock = starts[1]; iyBlock <= ends[1]; iyBlock++) {
    for (ixBlock = starts[0]; ixBlock <= ends[0]; ixBlock++) {
      blkInd = 2 * (ixBlock + iyBlock * mNumBlocks[mBestScale][b3dX]);
      if (mBoundaries[blkInd + botTop] >= 0) {
        indWgt++;
        if (wgtColIn > 0 && mFitMat[wgtColIn][indWgt] < wgtThresh)
          continue;
        if (wgtColIn > 0 && wgtColOut > 0)
          mFitMat[wgtColOut][numFit] = mFitMat[wgtColIn][indWgt];
        for (var = 0; var < mNumVars; var++) {
          ixy = mVarList[var][0];
          mFitMat[var][numFit] = (mBlockCenters[blkInd + ixy] - 
                                  mBlockCenters[cenInd + ixy]) / xySpacing;
          ixy = mVarList[var][1];
          if (ixy >= 0)
            mFitMat[var][numFit] *= (mBlockCenters[blkInd + ixy] - 
                                         mBlockCenters[cenInd + ixy]) / xySpacing;
        }
        mFitMat[mNumVars][numFit] = mBoundaries[blkInd + botTop];
        numFit++;
      }
    }
  }
}

/*
 * Sets up indices for a polynomial fit appropriate to the number of points and
 * number of "good" points with higher weights.  A term can be a single variable (if
 * varlist[][1] is -1 or a product of two variables in varList[][0] and varList[1].
 */
void FindSect::buildVariableList(int major, int minor, int numBound, int numGood)
{
  mNumVars = 1;
  mVarList[0][0] = major;
  mVarList[0][1] = -1;
  if (numBound >= 10 && numGood > 8 && mGoodRowCol[minor] >= 2) {
    mVarList[1][0] = minor;
    mVarList[1][1] = -1;
    mNumVars = 2;
  }
  if (numBound >= 14 && numGood > 11 && mGoodRowCol[minor] >= 2) {
    mVarList[2][0] = major;
    mVarList[2][1] = major;
    mNumVars = 3;
  }
  if (numBound >= 20 && numGood > 17 && mGoodRowCol[minor] >= 4) {
    mVarList[3][0] = major;
    mVarList[3][1] = minor;
    mVarList[4][0] = minor;
    mVarList[4][1] = minor;
    mNumVars = 5;
  }
}

/*
 * Does samples through the entire thickness extent with the size of the edge samples and 
 * searches for a peak in the middle and dips on either side.
 */
int FindSect::findLowestSDEdges(int *starts, int *ends, int scl)
{
  int numStat, numSamples, sampleThick, dividedExtent, sampSpacing, samp, ind;
  double meanSum;
  float peak, dip;
  int dipInd[2], peakInd, direc, numAbove;
  float aboveDipMADNcrit = 2.;
  int aboveDipNumCrit = 3;
  vector<float> edgeMedians, edgeMADNs;

  // Get thickness, number of samples and spacing
  sampleThick =  starts[3];
  starts[mThickInd] = 0;
  ends[mThickInd] = starts[3] - 1;
  sampSpacing = B3DMAX(1, sampleThick / 2);
  dividedExtent = mNumBoxes[scl][mThickInd] - (sampleThick - sampSpacing);
  numSamples = (dividedExtent + sampSpacing - 1) / sampSpacing;
  edgeMedians.resize(numSamples);
  edgeMADNs.resize(numSamples);
  for (samp = 0; samp < numSamples; samp++) {
    numStat = 0;
    meanSum = 0.;
    ends[mThickInd] = B3DMIN(samp * sampSpacing + sampleThick - 1, ends[4]);
    starts[mThickInd] = ends[mThickInd] + 1 - sampleThick;
    addBoxesToSample(starts[0], ends[0], starts[1], ends[1], starts[2], ends[2], 
                     scl, numStat, meanSum);
    rsFastMedian(mBuffer, numStat, &mBuffer[numStat], &edgeMedians[samp]);
    rsFastMADN(mBuffer, numStat, edgeMedians[samp], &mBuffer[numStat], &edgeMADNs[samp]);
    if (mDebugOutput > 1)
      printf("%d %d %8.3f %8.3f\n", starts[mThickInd], ends[mThickInd], edgeMedians[samp],
             edgeMADNs[samp]);
  }

  // Find peak, disqualifying a peak at the end of the range
  peak = -1.;
  peakInd = -1;
  for (samp = 1; samp < numSamples - 1; samp++) {
    if (edgeMedians[samp] > edgeMedians[samp - 1] && 
        edgeMedians[samp] > edgeMedians[samp + 1] && edgeMedians[samp] > peak) {
      peak = edgeMedians[samp];
      peakInd = samp;
    }
  }
  if (peakInd < 0)
    return 1;

  // Find minima on each side of peak
  for (direc = -1; direc <= 1; direc += 2) {
    ind = (direc + 1) / 2;
    samp = peakInd + direc;
    dip = peak + 10.;
    while (samp >= 0 && samp < numSamples) {
      if (edgeMedians[samp] < dip) {
        dip = edgeMedians[samp];
        dipInd[ind] = samp;
      }
      samp += direc;
    }
    if (dip >= peak)
      return 2;
  }
  
  // Sanity check: are the two dips on their respective sides of the middle?
  if (dipInd[0] >= numSamples / 2 || dipInd[1] <= numSamples / 2)
    return 3;

  // And structure test: how many samples are above n MADNs from the dip?
  numAbove = 0;
  ind = dipInd[0];
  if (edgeMedians[ind] > edgeMedians[dipInd[1]])
    ind = dipInd[1];
  for (samp = dipInd[0] + 1; samp < dipInd[1]; samp++)
    if ((edgeMedians[samp] - edgeMedians[ind]) / edgeMADNs[ind] > aboveDipMADNcrit)
      numAbove++;
  if (numAbove < aboveDipNumCrit)
    return 4;

  // Set the best edges to middle of dip samples
  if (dipInd[0] > 0)
    mBestLowEdge[scl] = (dipInd[0] + 1) * sampSpacing;
  if (dipInd[1] < ends[4])
    mBestHighEdge[scl] = (dipInd[1] + 1) * sampSpacing;

  // Set the limits for the edge samples; i.e. set absolute limits in s/e[4] and 
  // make up excluded region in s/e[3]
  starts[4] = mBestLowEdge[scl];
  ends[4] = mBestHighEdge[scl];
  starts[3] = starts[4] + sampleThick;
  ends[3] = ends[4] - sampleThick;
  if (mDebugOutput)
    printf("Set new limits for edge at %d %d\n", mBestLowEdge[scl], mBestHighEdge[scl]);
  return 0;
}

#define MAX_AMVAR 2
#define MAX_MULT 16
/*
 * Separate analysis of thickness based on distribution of high SD values and on optional
 * bead model
 */
int FindSect::analyzeHighSD(float highSDcrit, Imod *beadModel, int *boxSize, 
                            int *startCoord, int *endCoord, Imod *pitchModel) 
{
  int boxInd, ixBox, iyBox, izBox, indObj, indCont, indPt, ilong, thick, ix, iy;
  int ind, midInd, indMax, mult, binWidth, numBins, bestMult, dir, midAtRise, numLow;
  int redo, numFit, numAbove;
  float crossCrit[2], crossBase[2], peakSlope[2], projSlope[2], projIntcp[2];
  float baseAtCrit[2];
  float *projStat, *projTemp;
  float maxBin, sem, crit, beadY, projMin, projPeak, extrap, base, aboveCrit, ro;
  float cumul, baseAvg, baseSD, zCol, zfrac, xCoeff, yCoeff, zCoeff, xyConst;
  float zSlice, xCol, yCol, slope, intcp, slopeDiff, peakMean, edgeMean;
  int numDips[MAX_MULT];
  Icont *cont;
  Ipoint *pts;
  Ipoint levelPt, pitchPt;
  float *beadCenters, *beadBounds;
  int scl = mBestScale;
  int loop, minBase;
  int yInd = 3 - mThickInd;
  bool flipped = mThickInd != 2;
  float da[MAX_AMVAR] = {2., 2.};
  float a[MAX_AMVAR] = {0., 0.};
  float yy[MAX_AMVAR + 1];
  float ptolFacs[MAX_AMVAR] = {0.05f, .003f}, ftolFacs[MAX_AMVAR] = {5.e-4f, 1.e-5f};
  float delfac = 2.;
  int iter, nvar;
  float alpha, beta, thickness, colMin, colMax, colVal, spread, cosAlpha, sinAlpha;
  float *bins, *comboBins;
  float *projSlice, *projMeans, *projMedians;
  float *projPct75, *projPct90;
  int numInSlice, iz;
  float monoCrit = 0.1;
  float minBaseFrac = 0.002;
  float allRisingCrit = 3.;
  float riseBackoffCrit = 2.;
  int maxRiseBackoff = 3;
  float maxLowFracPastRise = 0.15;
  float maxBumpLowFrac = 0.1;
  float bumpSum, bumpCrit;
  float bumpSumCritFac = 20.;
  bool goodBump;
  float xcen = mNxyz[0] / 2.;
  float ycen = mNxyz[yInd] / 2.;
  float zcen = mNxyz[mThickInd] / 2.;
  int *binnings = &mBinning[scl][0];
  int *spacings = &mBoxSpacing[scl][0];
  int *numBoxes = &mNumBoxes[scl][0];
  Imat *mat = imodMatNew(3);
  nvar = 2;
  
  mBestBoxSize = boxSize;
  mStartCoord = startCoord;
  mNumBeads = 0;
  if (beadModel)  {
    for (indObj = 0; indObj < beadModel->objsize; indObj++)
      for (indCont = 0; indCont < beadModel->obj[indObj].contsize; indCont++)
        mNumBeads += beadModel->obj[indObj].cont[indCont].psize;
  }

  // Loop on columns twice looking for points above threshold
  for (loop = 0; loop < 2; loop++) {
    mNumHighSD = 0;
    for (ilong = 0; ilong < numBoxes[yInd]; ilong++) {
      for (ixBox = 0; ixBox < numBoxes[b3dX]; ixBox++) {
        colMin = 1.e10;
        colMax = -1.e10;
        for (thick = mBestLowEdge[scl]; thick <= mBestHighEdge[scl]; thick++) {
          izBox = flipped ? ilong : thick;
          iyBox = flipped ? thick : ilong;
          boxInd = mStatStartInds[scl] + (izBox * numBoxes[b3dY] + iyBox) * 
            numBoxes[b3dX] + ixBox;
          if ((mSDs[boxInd] - mEdgeMedians[scl]) / mEdgeMADNs[scl] > highSDcrit) {
            if (!loop) {
              
              // First time, just counting columns with such points
              mNumHighSD++;
              break;
            } else {

              // Second time, get Z value and keep track of min/max
              colVal = startCoord[mThickInd] + binnings[mThickInd] *
                (thick * spacings[mThickInd] + boxSize[mThickInd] / 2.);
              colMin = B3DMIN(colMin, colVal);
              colMax = B3DMAX(colMax, colVal);
            }
          }
        }

        // End of column, add the min and max points, 
        if (loop && colMax >= colMin) {
          mBlockCenters[2 * mNumHighSD] = startCoord[b3dX] + binnings[b3dX] *
            (ixBox * spacings[b3dX] + boxSize[b3dX] / 2.);
          mBlockCenters[2 * mNumHighSD + 1] = startCoord[yInd] + binnings[yInd] *
            (ilong * spacings[yInd] + boxSize[yInd] / 2.);
          mBoundaries[2 * mNumHighSD] = colMin;
          mBoundaries[2 * mNumHighSD + 1] = colMax;
          mNumHighSD += 1;
        }
      }
    }

    // First time allocate
    if (!loop) {
      if (!mNumHighSD)
        exitError("There are no boxes with SD above the criterion value");
      if (mDebugOutput)
        printf("numhigh %d\n", mNumHighSD);
      mBlockCenters = B3DMALLOC(float, 2 * mNumHighSD + 2 * mNumBeads);
      mBoundaries = B3DMALLOC(float, 2 * mNumHighSD + mNumBeads);
      mBoundRot = B3DMALLOC(float, 3 * mNumHighSD + mNumBeads);
      bins = B3DMALLOC(float, mNxyz[mThickInd] / spacings[mThickInd] + 5);
      comboBins = B3DMALLOC(float, mNxyz[mThickInd] / spacings[mThickInd] + 5);
      projSlice =  B3DMALLOC(float, numBoxes[yInd] * numBoxes[b3dX]);
      ind = numBoxes[mThickInd];
      projMeans = B3DMALLOC(float, ind);
      projMedians = B3DMALLOC(float, ind);
      projPct75 = B3DMALLOC(float, ind);
      projPct90 = B3DMALLOC(float, ind);
      if (!mBlockCenters || !mBoundaries || !mBoundRot || !bins || !comboBins ||
          !projSlice || !projMeans || ! projMedians || !projPct75 || !projPct90)
        exitError("Allocating arrays for high-SD extreme values");
      mProjSlice = projSlice;
      mProjMeans = projMeans;
    }
  }
  beadCenters = &mBlockCenters[2 * mNumHighSD];
  beadBounds = &mBoundaries[2 * mNumHighSD];

  // Put bead positions on the ends of these arrays
  if (beadModel) {
    mNumBeads = 0;
    for (indObj = 0; indObj < beadModel->objsize; indObj++) {
      for (indCont = 0; indCont < beadModel->obj[indObj].contsize; indCont++) {
        cont = &beadModel->obj[indObj].cont[indCont];
        pts = cont->pts;
        for (indPt = 0; indPt < cont->psize; indPt++) {
          beadY = flipped ? pts[indPt].z : pts[indPt].y;
          if (pts[indPt].x >= startCoord[b3dX] && pts[indPt].x <= endCoord[b3dX] &&
              beadY >= startCoord[yInd] && beadY <= endCoord[yInd]) {
            beadCenters[mNumBeads * 2] = pts[indPt].x;
            beadCenters[mNumBeads * 2 + 1] = beadY;
            beadBounds[mNumBeads++] = flipped ? pts[indPt].y : pts[indPt].z;
          }
        }
      }
    }
    if (mNumBeads) {
      ind = 2 * B3DMAX(2 * mNumHighSD, mNumBeads);
      mConvexXtmp = B3DMALLOC(float, ind);
      mConvexYtmp = B3DMALLOC(float, ind);
      mConvexCont = imodContourNew();
      if (mConvexCont)
        mConvexCont->pts = B3DMALLOC(Ipoint, ind);
      if (!mConvexXtmp || !mConvexYtmp || !mConvexCont || !mConvexCont->pts)
        exitError("Allocating arrays for convex hulls");
    }
  }

  // Run amoeba fit twice and get solution
  mSimplexIter = 0;
  mAfterSimplex = false;
  dualAmoeba(yy, nvar, delfac, ptolFacs, ftolFacs, a, da, amoebaFunc, &iter);
  mAfterSimplex = true;
  amoebaFunc(a, &spread);

  // Get final rotation matrix first for forward then for back rotation
  alpha = a[0];
  beta = a[1];
  imodMatRot(mat, -beta, b3dY);
  imodMatRot(mat, -alpha, b3dX);
  cosAlpha = cos(alpha * RADIANS_PER_DEGREE);
  sinAlpha = sin(alpha * RADIANS_PER_DEGREE);
  xCoeff = mat->data[2];
  yCoeff = mat->data[6];
  zCoeff = mat->data[10];
  imodMatId(mat);
  imodMatRot(mat, alpha, b3dX);
  imodMatRot(mat, beta, b3dY);

  if (mProjLayerStatType > 0) {

    // If using projection layers, fetch data from slices at this rotation and average them
    for (thick = 0; thick < numBoxes[mThickInd]; thick++) {
      numInSlice = 0;
      numAbove = 0;
      zSlice = startCoord[mThickInd] + binnings[mThickInd] *
        (thick * spacings[mThickInd] + boxSize[mThickInd] / 2.);
      for (ilong = 0; ilong < numBoxes[yInd]; ilong++) {
        for (ixBox = 0; ixBox < numBoxes[b3dX]; ixBox++) {
          xCol = startCoord[b3dX] + binnings[b3dX] *
            (ixBox * spacings[b3dX] + boxSize[b3dX] / 2.);
          yCol = startCoord[yInd] + binnings[yInd] *
            (ilong * spacings[yInd] + boxSize[yInd] / 2.);
          xyConst = xCoeff * (xCol - xcen) + yCoeff * (yCol - ycen);
          zCol = ((zSlice - zcen - xyConst) / zCoeff + zcen) / binnings[mThickInd];
          if (zCol >= 0 && zCol <= numBoxes[mThickInd]) {
            iz = zCol;
            zfrac = zCol - iz;
            izBox = flipped ? ilong : iz;
            iyBox = flipped ? iz : ilong;
            boxInd = mStatStartInds[scl] + (izBox * numBoxes[b3dY] + iyBox) * 
              numBoxes[b3dX] + ixBox;
            projSlice[numInSlice] = (1. - zfrac) * mSDs[boxInd];

            izBox = flipped ? ilong : iz + 1;
            iyBox = flipped ? iz + 1 : ilong;
            boxInd = mStatStartInds[scl] + (izBox * numBoxes[b3dY] + iyBox) * 
              numBoxes[b3dX] + ixBox;
            projSlice[numInSlice] = (projSlice[numInSlice] + zfrac * mSDs[boxInd] -
                                     mEdgeMedians[scl]) / mEdgeMADNs[scl];
            if (projSlice[numInSlice] > highSDcrit)
              numAbove++;
            numInSlice++;
          }
        }
      }
      avgSD(projSlice, numInSlice, &projMeans[thick], &colMin, &colMax);
      rsFastMedianInPlace(projSlice, numInSlice, &projMedians[thick]);
      projPct75[thick] = percentileFloat(B3DNINT(0.75 * numInSlice), projSlice, 
                                         numInSlice);
      projPct90[thick] = (float)numAbove / (float)numInSlice;
      if (mDebugOutput > 1)
        printf("%d  %.3f  %.3f %.3f %.3f\n", B3DNINT(zSlice), projMeans[thick], 
               projMedians[thick], projPct75[thick], projPct90[thick]);
    }

    // Assign the chosen array
    projTemp = projPct90;
    if (mProjLayerStatType == 1)
      projStat = projMeans;
    else if (mProjLayerStatType == 2)
      projStat = projMedians;
    else if (mProjLayerStatType == 3)
      projStat = projPct75;
    else {
      projStat = projPct90;
      projTemp = projPct75;
    }

    // Find peak of projected values
    indMax = -1;
    for (ind = 0; ind <  numBoxes[mThickInd]; ind++) {
      if (indMax < 0 || projStat[ind] > projPeak) {
        projPeak =  projStat[ind];
        indMax = ind;
      }
    }
    peakMean = (projStat[indMax] + projStat[B3DMAX(0, indMax - 1)] + 
                projStat[B3DMAX(0, indMax + 1)]) / 3.;
    ind = numBoxes[mThickInd] - 3;
    edgeMean = (projStat[0] + projStat[1] + projStat[2] + projStat[ind] + projStat[ind + 1]
                + projStat[ind + 2]) / 6.;
    if (mDebugOutput)
      printf("Peak %.3f  edge %.3f  difference %.3f\n", peakMean, edgeMean, 
             peakMean - edgeMean);

    // Fit to the bottom and the top at the edges
    numFit = B3DMAX(10, mBestLowEdge[scl] + 1);
    for (ind = 0; ind < numBoxes[mThickInd]; ind++)
      projTemp[ind] = ind;
    lsFit(projTemp, projStat, numFit, &projSlope[0], &projIntcp[0], &ro);
    ind = numFit;
    numFit = B3DMAX(10, numBoxes[mThickInd] - mBestHighEdge[scl]);
    thick = numBoxes[mThickInd] - numFit;
    lsFit(&projTemp[thick], &projStat[thick], numFit, &projSlope[1], &projIntcp[1], &ro);

    // For each side, walk in from fit looking for point that is up by the criterion above
    // the baseline, where the base is the maximum of the minimum value seen and the
    // current extrapolation from the baseline
    for (loop = 0; loop < 2; loop++) {
      projMin = projStat[ind];
      dir = 1 - 2 * loop;
      crossCrit[loop] = -1.;
      while (dir * (indMax - ind) > 0) {
        extrap = projSlope[loop] * ind + projIntcp[loop];
        projMin = B3DMIN(projMin, projStat[ind]);
        base = B3DMAX(projMin, extrap);
        if (mDebugOutput > 1)
          PRINT4(extrap, projMin, base, projStat[ind]);
        aboveCrit = projStat[ind] - (base + mProjEdgeCrit * (projPeak - base));

        // Get position where it crosses criterion, and extrapolate back to baseline
        // from immediate slope.  Also get slope from here to peak
        if (aboveCrit >= 0) {
          peakSlope[loop] = (projPeak - projStat[ind]) / fabs((double)ind - indMax);
          crossCrit[loop] = ind - dir * aboveCrit / (projStat[ind] - projStat[ind + dir]);
          lsFit(&projTemp[ind - 2 * loop], &projStat[ind - 2 * loop], 3, &slope, &intcp, 
                &ro);
          crossBase[loop] = (base - intcp) / slope;
          baseAtCrit[loop] = base;
          break;
        }
        ind += dir;
      }

      // Starting point for other direction
      ind = thick - 1;
    }

    // If only one side exists or if other side has steep slope, redo that side with
    // the base from the good side
    if (crossCrit[0] < 0 && crossCrit[1] < 0)
      return 1;
    slopeDiff = fabs(fabs(projSlope[0]) - fabs(projSlope[1]));
    redo = -1;
    if (crossCrit[0] < 0 || peakSlope[0] - projSlope[0] < slopeDiff)
      redo = 0;
    if (crossCrit[1] < 0 || peakSlope[1] - projSlope[1] < slopeDiff)
      redo = 1;
    if (redo >= 0) {
      crossCrit[redo] = -1;
      base = baseAtCrit[1 - redo];
      ind = redo > 0 ?  numBoxes[mThickInd] - 1 : 0;
      dir = 1 - 2 * redo;
      while (dir * (indMax - ind) > 0) {
        aboveCrit = projStat[ind] - (base + mProjEdgeCrit * (projPeak - base));
        if (aboveCrit >= 0) {
          crossCrit[redo] = ind - dir * aboveCrit / (projStat[ind] - projStat[ind + dir]);
          lsFit(&projTemp[ind - 2 * redo], &projStat[ind - 2 * redo], 3, &slope, &intcp, 
                &ro);
          crossBase[redo] = (base - intcp) / slope;
          break;
        }
        ind += dir;
      }
      if (crossCrit[redo] < 0)
        return 1;
    }

    // scale both sets of numbers
    for (ind = 0; ind < 2; ind++) {
      crossCrit[ind] = startCoord[mThickInd] + binnings[mThickInd] *
        (crossCrit[ind] * spacings[mThickInd] + boxSize[mThickInd] / 2.);
      crossBase[ind] = startCoord[mThickInd] + binnings[mThickInd] *
        (crossBase[ind] * spacings[mThickInd] + boxSize[mThickInd] / 2.);
    }
    printf("Surfaces at crossing point %.1f  %.1f  thick %1.f\n", crossCrit[0], 
           crossCrit[1], crossCrit[1] - crossCrit[0]);
    printf("Surfaces extrapolated to base %.1f  %.1f  thick %1.f\n", crossBase[0], 
           crossBase[1], crossBase[1] - crossBase[0]);
    if (mProjUseExtrap) {
      mBoundLow = crossBase[0];
      mBoundHigh = crossBase[1];
    } else {
      mBoundLow = crossCrit[0];
      mBoundHigh = crossCrit[1];
    }

  } else {

    // Otherwise (default) work with distributions of column boundaries

    // Analyze slightly smoothed histograms for minimal bumpiness above 0.1 of peak
    for (loop = 0; loop < 2; loop++) {
      for (mult = 1; mult < MAX_MULT; mult++) {
        makeCombinedBins(&mBoundRot[loop * mNumHighSD], mNumHighSD, bins, comboBins, mult,
                         numBins, binWidth, maxBin, indMax);
        if (mDebugOutput)
          PRINT4(numBins, binWidth, maxBin, indMax);

        // Count non-monotonic points above criterion level
        numDips[mult] = 0;
        for (ind = 1; ind < numBins - 1; ind++) {
          if (comboBins[ind] >= monoCrit * maxBin && comboBins[ind] < comboBins[ind - 1] &&
              comboBins[ind] < comboBins[ind + 1])
            numDips[mult]++;
        }
        if (mDebugOutput)
          PRINT2(mult, numDips[mult]);
        if (!numDips[mult])
          break;
      }
      bestMult = 1;
      for (mult = 1; mult < MAX_MULT; mult++) {
      
        // Done if 0 dips
        if (!numDips[mult]) {
          bestMult = mult;
          break;
        }
      
        // New minimum, take it.  Otherwise stop when no new minimum and 1 or 2 dips,
        // or 2nd time after no new minimum with 3 dips
        if (numDips[mult] < numDips[bestMult]) {
          bestMult = mult;
        } else if ((numDips[bestMult] < 3 && mult - bestMult > 0) || 
                   (numDips[bestMult] == 3 && mult - bestMult >= 2)) {
          break;
        }
      }
      if (bestMult == 16) {
        return 1;
      }
      if (mDebugOutput)
        PRINT1(bestMult);
      makeCombinedBins(&mBoundRot[loop * mNumHighSD], mNumHighSD, bins, comboBins, 
                       bestMult, numBins, binWidth, maxBin, indMax);
      if (indMax <= 3 || indMax >= numBins - 4) {
        return 1;
      }
      if (mDebugOutput > 1) {
        for (ind = 0; ind < numBins; ind++)
          printf("%d  %.1f\n", ind, comboBins[ind]);
      }

      // find minimum baseline
      minBase = B3DMAX(2, minBaseFrac * mNumHighSD);
      if (loop) {
        dir = -1;
        ind = numBins - 1;
      } else {
        dir = 1;
        ind = 0;
      }

      // Advance past minimum but make sure not to get too close to peak
      cumul = 0;
      for (; cumul < minBase; ind += dir)
        cumul += bins[ind];
      if (loop)
        B3DCLAMP(ind, indMax + 4, numBins - 3);
      else
        B3DCLAMP(ind, 3, indMax - 4);
      midInd = ind;

      // Look for point where it rises; get baseline mean and SD
      midAtRise = -1;
      for (; dir * (indMax - 3 - midInd) > 0; midInd += dir) {
        avgSD(&comboBins[loop ? midInd : 0], loop ? (numBins - 1 - midInd) : midInd, 
              &baseAvg, &baseSD, &sem);
        crit = baseAvg + allRisingCrit * baseSD;

        // If next three bins are all above triggering criterion, make sure bins stay 
        // high to peak, or that a contiguous "bump" contains a lot of counts above
        // the criterion
        if (baseSD > 0 && comboBins[midInd + dir] > crit && 
            comboBins[midInd + 2 * dir] > crit && comboBins[midInd + 3 * dir] > crit) {
          bumpCrit = bumpSumCritFac * crit;
          crit = baseAvg + riseBackoffCrit * baseSD;
          numLow = 0;
          bumpSum = 0.;
          goodBump = false;
          for (ind = 1; ind < dir * (indMax - midInd) ; ind++) {
            if (comboBins[midInd + ind * dir] <= crit)
              numLow++;
            bumpSum += comboBins[midInd + ind * dir];
            if (bumpSum > bumpCrit && 
                numLow <= B3DNINT(maxBumpLowFrac * ind))
              goodBump = true;
          }
          if (!goodBump && 
              numLow > B3DNINT(maxLowFracPastRise * (dir * (indMax - midInd) - 4)))
            continue;

          // If it passes that, then try to back off from higher criterion
          for (ind = 0; ind < maxRiseBackoff; ind++) {
            if (comboBins[midInd] <= crit)
              break;
            midInd -= dir;
          }
          midAtRise = midInd;
          break;
        }
      }

      // Error if no good spot found
      if (ind < 0) {
        return 1;
      }
      if (mDebugOutput)
        printf("%s rise after bin %d\n", loop ? "High:" : "Low: ", midAtRise);
    
      if (loop)
        mBoundHigh = (midAtRise) * binWidth;
      else
        mBoundLow = (midAtRise + 1) * binWidth;
    }
  }

  if (mBoostHighSDThickness > 0.) {
    zfrac = 0.5 * mBoostHighSDThickness * (mBoundHigh - mBoundLow);
    mBoundLow = B3DMAX(0, mBoundLow - zfrac);
    mBoundHigh = B3DMIN(mNxyz[mThickInd] - 1., mBoundHigh + zfrac);
  }

  // Now combine results from beads
  if (mNumBeads) {
    computeBeadLimits(mBeadExclPctlThick);
    mBoundLow = B3DMIN(mBoundLow, mBeadLow + zcen);
    mBoundHigh = B3DMAX(mBoundHigh, mBeadHigh + zcen);
  }

  thickness = mBoundHigh - mBoundLow;

  if (mDebugOutput > 0)
    printf("Iterations %d  spread measure %.2f\n", mSimplexIter, spread);
  ix = B3DNINT(mBoundLow);
  iy = B3DNINT(mBoundHigh);
  invertYifFlipped(ix, iy, mNxyz);
  printf("Thickness = %.1f  boundaries (from 1) = %d %d  alpha = %.2f  "
         "beta = %.2f\n", thickness, ix + 1, iy + 1, alpha, beta);
  if (!pitchModel)
    return 0;

  // Make the tomopitch model, end up with lines in the same slice
  for (iy = -1; iy <= 1; iy++) {
    for (loop = 0; loop < 2; loop++) {
      if (imodNewContour(pitchModel) || (cont = imodContourGet(pitchModel)) == NULL)
        exitError("Adding contour to model");
      pitchPt.y = iy * mNxyz[yInd] / 3.;
      levelPt.z = (loop ? mBoundHigh : mBoundLow) - zcen;
      levelPt.y = (pitchPt.y + levelPt.z * sinAlpha) / cosAlpha;
      for (ix = -1; ix <= 1; ix += 2) {
        levelPt.x = ix * mNxyz[0] * 0.84 / 2.;
        imodMatTransform(mat, &levelPt, &pitchPt); 
        if (!imodPointAppendXYZ(cont, pitchPt.x + xcen, 
                                flipped ? pitchPt.z + zcen : pitchPt.y + ycen, 
                                flipped ? pitchPt.y + ycen : pitchPt.z + zcen))
          exitError("Adding point to model");
      }
    }
  }
  return 0;
}

/* 
 * Make a simple histogram of the given depth values with the given multiplier of the 
 * basic Z spacing and then make simple weighted * combination of 3 adjacent bins in 
 * comboBins
 */
void FindSect::makeCombinedBins(float *values, int numVals, float *bins, float *comboBins,
                                int mult, int &numBins, int &binWidth, float &maxBin, 
                                int &indMax)
{
  float lastVal, firstVal;
  int ind;
  binWidth = mult * mBoxSpacing[mBestScale][mThickInd];
  numBins = mNxyz[mThickInd] / binWidth;
  firstVal = -0.5 * mNxyz[mThickInd];
  lastVal = firstVal + binWidth * numBins;
  kernelHistogram(values, numVals, bins, numBins, firstVal, lastVal, 0, 0);
      
  // Make combined bins and find peak
  indMax = -1;
  for (ind = 0; ind < numBins; ind++) {
    comboBins[ind] = 0.25 * bins[B3DMAX(0, ind - 1)] + 0.5 * bins[ind] + 
      0.25 * bins[B3DMIN(numBins - 1, ind + 1)];
    if (indMax < 0 || maxBin < comboBins[ind]) {
      indMax = ind;
      maxBin = comboBins[ind];
    }
  }
}

/*
 * The actual function for computing a normalized measure of spreading
 */
void FindSect::spreadFunc(float *yvec, float *spread)
{
  float xcen = mNxyz[0] / 2.;
  float ycen = mNxyz[3 - mThickInd] / 2.;
  float zcen = mNxyz[mThickInd] / 2.;
  Imat *mat = imodMatNew(3);
  float xCoeff, yCoeff, zCoeff, xyConst, weight, zeroDist;
  float zeroBoostLim = 0.05;
  float *beadBounds = &mBoundRot[2 * mNumHighSD];
  float *highBounds = &mBoundRot[mNumHighSD];
  float *temp = &mBoundRot[2 * mNumHighSD + mNumBeads];
  float medians[2], MADNs[2];
  int ind, lohi;
  int thick, numInSlice, ilong, ixBox, iyBox, izBox, boxInd, iz;
  int scl = mBestScale;
  int yInd = 3 - mThickInd;
  bool flipped = mThickInd != 2;
  int *binnings = &mBinning[scl][0];
  int *spacings = &mBoxSpacing[scl][0];
  int *numBoxes = &mNumBoxes[scl][0];
  float zSlice, xCol, yCol, zCol, zfrac, dum1, dum2, projMin;
  //bool doDebug = mAfterSimplex && mDebugOutput > 0;

  imodMatRot(mat, -yvec[1], b3dY);
  imodMatRot(mat, -yvec[0], b3dX);
  xCoeff = mat->data[2];
  yCoeff = mat->data[6];
  zCoeff = mat->data[10];

  // Compute the rotated boundaries, rotate around center to retain precision
  for (ind = 0; ind < mNumHighSD; ind++) {
    xyConst = xCoeff * (mBlockCenters[2 * ind] - xcen) +
      yCoeff * (mBlockCenters[2 * ind + 1] - ycen);
    mBoundRot[ind] = xyConst + zCoeff * (mBoundaries[2 * ind] - zcen);
    highBounds[ind] = xyConst + zCoeff * (mBoundaries[2 * ind + 1] - zcen);
  }

  // Rotate the beads if any
  for (ind = 0; ind < mNumBeads; ind++)
    beadBounds[ind] = xCoeff * (mBlockCenters[2 * (ind + mNumHighSD)] - xcen) + 
      yCoeff * (mBlockCenters[2 * (ind + mNumHighSD) + 1] - ycen) +
      zCoeff * (mBoundaries[ind + 2 * mNumHighSD] - zcen);

  // Compute median and MADN of low and high boundaries
  for (lohi = 0; lohi < 2; lohi++) {
    rsFastMedianInPlace(&mBoundRot[lohi * mNumHighSD], mNumHighSD, &medians[lohi]);
    rsFastMADN(&mBoundRot[lohi * mNumHighSD], mNumHighSD, medians[lohi], temp,
               &MADNs[lohi]);
  }

  // Fetch data from slices at this rotation and average them
  if (mUseProjForSpread) {
    projMin = 1.e20;
    for (thick = 0; thick < numBoxes[mThickInd]; thick++) {
      numInSlice = 0;
      zSlice = mStartCoord[mThickInd] + binnings[mThickInd] *
        (thick * spacings[mThickInd] + mBestBoxSize[mThickInd] / 2.);
      for (ilong = 0; ilong < numBoxes[yInd]; ilong++) {
        for (ixBox = 0; ixBox < numBoxes[b3dX]; ixBox++) {
          xCol = mStartCoord[b3dX] + binnings[b3dX] *
            (ixBox * spacings[b3dX] + mBestBoxSize[b3dX] / 2.);
          yCol = mStartCoord[yInd] + binnings[yInd] *
            (ilong * spacings[yInd] + mBestBoxSize[yInd] / 2.);
          xyConst = xCoeff * (xCol - xcen) + yCoeff * (yCol - ycen);
          zCol = ((zSlice - zcen - xyConst) / zCoeff + zcen) / binnings[mThickInd];
          if (zCol >= 0 && zCol <= numBoxes[mThickInd]) {
            iz = zCol;
            zfrac = zCol - iz;
            izBox = flipped ? ilong : iz;
            iyBox = flipped ? iz : ilong;
            boxInd = mStatStartInds[scl] + (izBox * numBoxes[b3dY] + iyBox) * 
              numBoxes[b3dX] + ixBox;
            mProjSlice[numInSlice] = (1. - zfrac) * mSDs[boxInd];

            izBox = flipped ? ilong : iz + 1;
            iyBox = flipped ? iz + 1 : ilong;
            boxInd = mStatStartInds[scl] + (izBox * numBoxes[b3dY] + iyBox) * 
              numBoxes[b3dX] + ixBox;
            mProjSlice[numInSlice] = (mProjSlice[numInSlice] + zfrac * mSDs[boxInd] -
                                      mEdgeMedians[scl]) / mEdgeMADNs[scl];
            numInSlice++;
          }
        }
      }
      mProjMeans[thick] = 0.;
      if (numInSlice)
        avgSD(mProjSlice, numInSlice, &mProjMeans[thick], &dum1, &dum2);
      projMin = B3DMIN(projMin, mProjMeans[thick]);
    }

    // Need centroid and moments of distribution
    double sum, tot, centroid, sum4;
    sum = 0.;
    tot = 0.;
    for (ind = 0; ind < numBoxes[mThickInd]; ind++) {
      sum += ind * (mProjMeans[ind] - projMin);
      tot += (mProjMeans[ind] - projMin);
    }
    centroid = sum / tot;
    sum = 0.;
    sum4 = 0.;
    for (ind = 0; ind < numBoxes[mThickInd]; ind++) {
      sum += pow(ind - centroid, 2) * (mProjMeans[ind] - projMin);
      sum4 += pow(ind - centroid, 4) * (mProjMeans[ind] - projMin);
    }
    sum = sqrt(sum / tot);
    sum4 = pow(sum4/ tot, .25);
    //PRINT3(centroid, sum, sum4);
    if (mUseProjForSpread > 1)
      *spread = sum4;
    else
      *spread = sum;
  }

  // Now do beads.  If the percentile is 0, just compute min/max
  if (mNumBeads) 
    computeBeadLimits(mBeadExclPctlSpread);

  // First time, initialize the normalization factors to the current values
  if (mSimplexIter == 0) {
    mMinSpread = 1.e10;
    mSpreadNorms[0] = MADNs[0];
    mSpreadNorms[1] = MADNs[1];
    if (mDebugOutput)
      PRINT2(mSpreadNorms[0], mSpreadNorms[1]);
    if (mNumBeads) {
      mSpreadNorms[2] = mBeadHigh - mBeadLow;

      // Measure area covered by structure and beads to base a weighting on that
      mStructArea = convexAreaCovered(mNumHighSD, 0);
      mBeadArea = convexAreaCovered(mNumBeads, mNumHighSD);
      if (mDebugOutput)
        PRINT3(mSpreadNorms[2], mStructArea, mBeadArea);
    }
  }

  // Compute spread as mean of normalized values, then take weighted average with
  // normalized thickness based on beads
  if (!mUseProjForSpread)
    *spread = (MADNs[0] / mSpreadNorms[0] + MADNs[1] / mSpreadNorms[1]) / 2.;
  if (mNumBeads) {
    weight = mBeadWeightFac * B3DMIN(1., mBeadArea / mStructArea);
    B3DCLAMP(weight, 0., 1.);
    *spread = *spread * (1. - weight) + weight * (mBeadHigh - mBeadLow) / mSpreadNorms[2];
  }
  zeroDist = sqrt(yvec[0] * yvec[0] + yvec[1] * yvec[1]);
  if (zeroDist < zeroBoostLim)
    *spread *= 1. + (zeroBoostLim - zeroDist) / zeroBoostLim;
  imodMatDelete(mat);
  mSimplexIter++;
  if (mDebugOutput > 1) {
    printf("%d %f %f %f %s", mSimplexIter, yvec[0], yvec[1], *spread,
           (*spread < mMinSpread) ? "*\n" : "\n");
  }
  mMinSpread = B3DMIN(mMinSpread, *spread);
}

/*
 * Compute the area occupied by either beads or structure from center positions
 */
float FindSect::convexAreaCovered(int numPts, int offset)
{
  int ind, numVert;
  float cbXcen, cbYcen;
  for (ind = 0; ind < numPts; ind++) {
    mConvexXtmp[ind] = mBlockCenters[2 * (offset + ind)];
    mConvexYtmp[ind] = mBlockCenters[2 * (offset + ind) + 1];
  }
  convexBound(mConvexXtmp, mConvexYtmp, numPts, 0.01, 0., &mConvexXtmp[numPts],
              &mConvexYtmp[numPts], &numVert, &cbXcen, &cbYcen, numPts);
  for (ind = 0; ind < numVert; ind++) {
    mConvexCont->pts[ind].x = mConvexXtmp[numPts + ind];
    mConvexCont->pts[ind].y = mConvexYtmp[numPts + ind];
    mConvexCont->pts[ind].z = 0.;
  }
  mConvexCont->psize = numVert;
  return imodContourArea(mConvexCont);
}

/*
 * Get mBeadLow and mBeadHigh limits with the given fraction excluded
 */
void FindSect::computeBeadLimits(float excludePctl)
{
  float *beadBounds = &mBoundRot[2 * mNumHighSD];
  int ind;

  // If no exclusion, find min/max
  if (excludePctl <= 0.) {
    mBeadLow = 1.e10;
    mBeadHigh = -1.e10;
    for (ind = 0; ind < mNumBeads; ind++) {
      mBeadLow = B3DMIN(mBeadLow, beadBounds[ind]);
      mBeadHigh = B3DMAX(mBeadHigh, beadBounds[ind]);
    }
    mBeadLow -= mBeadDiameter / 2.;
    mBeadHigh += mBeadDiameter / 2.;
  } else {

    // Otherwise compute the percentiles
    ind = B3DNINT(mNumBeads * excludePctl);
    B3DCLAMP(ind, 1, mNumBeads);
    mBeadLow = percentileFloat(ind, beadBounds, mNumBeads) - mBeadDiameter / 2.;
    ind = B3DNINT(mNumBeads * (1. - excludePctl));
    B3DCLAMP(ind, 1, mNumBeads);
    mBeadHigh = percentileFloat(ind, beadBounds, mNumBeads) + mBeadDiameter / 2.;
  }
}

/*
 * Function to call from amoeba
 */
static void amoebaFunc(float *yvec, float *spread)
{
  sFindSect.spreadFunc(yvec, spread);
}

/*
 * Callback function to load a needed slice for multiBinStats
 */
static int getSlice(int *iz, int *fdata, float *buffer)
{
  iiuSetPosition(fdata[0], *iz, 0);
  return iiuReadSecPart(fdata[0], (char *)buffer, fdata[2] + 1 - fdata[1], fdata[1],
                        fdata[2], fdata[3], fdata[4]);
}

