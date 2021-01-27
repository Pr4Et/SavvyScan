/*
 *  ctfphaseflip.cpp  -  CTF correction of tilted images
 *
 *  Authors: Quanren Xiong and David Mastronarde
 *
 *  Copyright (C) 2007-2018 by the Regents of the University of
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 *  $Id$
 */

#include <limits>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "b3dutil.h"
#include "iimage.h"
#include "sliceproc.h"
#include "ilist.h"
#include "ctfutils.h"
#include "cfft.h"
#include "parse_params.h"
#include "gpuctf.h"
#include "frameutil.h"
#include "cppdefs.h"

using namespace std;

#define MIN_ANGLE 0.001  //tilt Angle less than this in degrees is treated as 0.0;
#define MY_PI 3.1415926
#define UNUSED_DEFOCUS -2000000.0
#define MAXLINE 160

#define START_TIMER  if (debugMode) wallStart = wallTime();
#define ADD_TIME(a) if (debugMode)                                         \
  {curTime = wallTime(); a += curTime - wallStart; wallStart = curTime;}


static double getAnyOneZero(double defocus, double phase, int zeroNum, double ampContrast,
                            double Cs, double pixelSize, double voltage);
static double firstZeroShift(double defocus, double phase, double angle, int extent,
                             double ampContrast, double Cs, double pixelSize,
                             double voltage);
static void interpolateTable(float *defocus, int nz, bool ifAngles);

// Macro for copying pixels along a diagonal
#define COPY_DIAGONALS                          \
  if (useGPU < 0) {                                         \
    for (row = 0; row < nyfile; row++) {                    \
      for (column = 0; column < nx; column++) {                           \
        axisDist = -sinViewAxis * (column - xPixCenter) + cosViewAxis * \
          (row - yPixCenter);                                                 \
        if (axisDist >= lowLim && axisDist <= highLim)                  \
          restoredArray[row * nx + column] =                            \
            curStrip[(row + yoff) * stripXdim + column + xoff];         \
      }                                                                   \
    }                                                                   \
  } else {                                                            \
    if (gpuCopyDiagonals(stripInd, xoff, yoff, sinViewAxis, cosViewAxis, \
                         lowLim, highLim))                          \
      exitError("Calling gpuCopyDiagonals for view %d", view);      \
  }


/* Notes on the angle of axis of constant Z with x-tilt
   a = x-axis tilt in IMOD, originally around Y axis in microscope
   b = tilt angle around Y in IMOD, but originally (in microscope) around X axis
   In original microscope space: 
   Points on a tilted plane that is level in the Y dimension (tilted only around Y)
   are defined by:
   (1)  (x, y, x tan(a))
   A point at (x, y, z) at 0 deg is tilted to (x', y', z'):
   (2)  (x, y cos(b) - z sin(b), y sin(b) + z cos(b))
   Points on the tilted plane are defined by substituting x tan(a) for z in (2):
   (3)  (x, y cos(b) - x tan(a) sin(b), y sin(b) + x tan(a) cos(b))
   The axis of equal z-height is the set of points with z' = 0, which in (3) satisfies:
        0 = y sin(b) + x tan(a) cos(b)
   (4)  y = -x tan(a) / tan(b)
   Substituting this into 3 gives:
        (x, -x (tan(a) / tan(b)) cos(b) - x tan(a) sin(b), 0)
        (x, -x tan(a) (cos(b)^2 / sin(b) + sin(b)), 0)
        (x, -x (tan(a) / sin(b)) (cos(b)^2 + sin(b)^2), 0)
   (5)  (x, -x tan(a) / sin(b), 0)
   This is a line of slope -tan(a) / sin(b)
   To rotate by 90 to bring the tilt axis vertical, apply 90 deg rotation to (5):
   (6)  (x tan(a) / sin(b), x, 0)
   This is now a line of slope sin(b) / tan(a)  and the angle g of the line of equal 
   z-height is
   (7)  g = atan(sin(b) / tan(a))
   However, the atan will put the angle in the first and fourth quadrants.  Since this 
   angle should converge toward the 90 degree axis at high tilt, the angle must be kept
   in the first and second quadrants instead, or the axis direction will be wrong.

   The same result can be obtained by using the rotated coordinates in the first place,
   as is done in IMOD, where the plane is 
   (1')   (0, 0, y tan(a)) 
   and the result of tilting (x, y, z) is the less intuitive
   (2')   (x cos(b) + z sin(b), y, -x sin(b) + z cos(b))
   and the tilted plane consists of points
   (3')   (x cos(b) + y tan(a) sin(b), y, -x sin(b) + y tan(a) cos(b))
   
   The distance of a point (x', y') in the image to this axis is its Y coordinate after 
   rotating to bring the axis to the X axis, 
       -x' sin(g) + y' cos(g)

   Note that if there is no X-axis tilt and the tilt axis is vertical (at 90 degrees), 
   the distance from the axis for points to the right of the tilt axis is -x'

   How does Dz, the change in z, depend on distance from the center of the image 
   perpendicular to this axis?  The line of points with z' = Dz satisfies:
      Dz = -x sin(b) + y tan(a) cos(b)    in original x, y, z IMOD space
   and we can pick the point originally at y = 0 and x = -Dz / sin(b), which
   substituted into 3' rotates to
   (-Dz cos(b) / sin(b), 0, Dz).  This point's distance from the axis at angle g is
        d = (Dz cos(b) / sin(b)) sin(g) = Dz sin(g) / tan(b)
   (8)  Dz = d tan(b) / sin(g)
   where the case of b = 0 and g = 0 has to be handled separately.
   But the change in defocus is the negative of this Dz

   Finally, if the tilt axis is rotated by g' from the vertical, i.e. in unaligned images,
   then as long as the X-axis tilt is defined as being around the axis perpendicular to 
   the tilt axis, i.e., around the microscope Y axis, the coordinate system 
   can be rotated by -g' to make the equations above apply.  The angle of the axis of
   equal z-height is thus at g + g'.  d must be computed from g + g' but the relation 
   between Dz and d depends on g alone.
*/

/*
 * One big main
 */
int main(int argc, char *argv[])
{
  int numOptArgs, numNonOptArgs;

  // Fallbacks from   ../manpages/autodoc2man 2 1 ctfphaseflip
  int numOptions = 30;
  const char *options[] = {
    "input:InputStack:FN:", "output:OutputFileName:FN:", "angleFn:AngleFile:FN:",
    "invert:InvertTiltAngles:B:", "axis:AxisAngle:F:", "xtilt:XAxisTilt:F:",
    "defFn:DefocusFile:FN:", "zoff:OffsetInZ:F:", "xform:TransformFile:FN:",
    "defTol:DefocusTol:I:", "maxWidth:MaximumStripWidth:I:",
    "iWidth:InterpolationWidth:I:", "zero:MinimumZeroSpacing:I:",
    "pixelSize:PixelSize:F:", "volt:Voltage:I:", "cs:SphericalAberration:F:",
    "ampContrast:AmplitudeContrast:F:", "phase:PhasePlateShift:F:",
    "degPhase:PhaseShiftInDegrees:F:", "cuton:CutOnFrequency:F:",
    "scale:ScaleByCtfPower:F:", "gpu:UseGPU:I:", "action:ActionIfGPUFails:IP:",
    "views:StartingEndingViews:IP:", "totalViews:TotalViews:IP:",
    "boundary:BoundaryInfoFile:FN:", "skipFlag:SkipCorrectedFlag:I:",
    "debug:DebugOutput:I:", "param:Parameter:PF:", "help:usage:B:"};

  char *stackFn, *angleFn, *outFn, *defFn, *xformFn = NULL;
  char *boundFn = NULL;
  char *gpuEnv;
  int volt, iWidth, defocusTol, ii, ierr, defFlags, maxWidth, minWidth;
  float pixelSize, cs, ampContrast, stripDefocus, offsetInZ = 0.;
  int startingView, endingView, startingTotal, endingTotal, defVersion;
  bool isSingleRun = false;
  int invertAngles = 0;
  int maxStripWidth = 0;
  int minStripWidth = 128;
  int minWidthToScaleInterp = 256;  // The former maximum strip width
  float minZeroShift = 0.6f;
  float freqForInterZero = 0.8f;    // Nyquist units at which to assess interzero distance
  float minInterZeroPixels = 8.;
  float attenStartFrac = 0.81f;
  float tiltAxisAngle = 0.;
  float xAxisTilt = 0.;
  int useGPU = -1;
  int skipDoneFlag = 0;
  bool ifGpuByEnv = false;
  int actGpuFailOption = 0, actGpuFailEnviron = 0;
  int debugMode = 0;
  int hushDefoci = 0;
  double angleSign, zeroShift;
  bool maxWasEntered, scaleInterp, doFullImages, doDiagonals, angleIsZero;
  float minAngle, maxAngle, scaleByPower = 0., phasePlateShift = 0., cutOnEntered = 0.;
  float phaseDeg, gpuMemory;
  float *tiltAngles = NULL;
  float *xShifts = NULL;
  float *rotations = NULL;
  Ilist *defocusList;
  char line[MAXLINE];
  char *progname = imodProgName(argv[0]);

  PipReadOrParseOptions(argc, argv, options, numOptions, progname,
                        1, 0, 0, &numOptArgs, &numNonOptArgs, NULL);

  if (!PipGetBoolean("usage", &ierr)) {
    PipPrintHelp(progname, 0, 0, 0);
    exit(0);
  }
  PipGetInteger("DebugOutput", &debugMode);
  if (debugMode >= 10) {
    hushDefoci = 1;
    debugMode -= 10;
  }
  if (PipGetString("InputStack", &stackFn))
    exitError("No stack specified");
  if (PipGetString("AngleFile", &angleFn)) {
    angleFn = NULL;
    printf("No angle file is specified, tilt angle is assumed to be 0.0\n");
  }
  if (PipGetString("DefocusFile", &defFn))
    exitError("No defocus file is specified");

  // Get axis angle and X axis tilt and decide if doing diagonals
  PipGetFloat("AxisAngle", &tiltAxisAngle);
  PipGetFloat("XAxisTilt", &xAxisTilt);
  doDiagonals = (fabs(tiltAxisAngle) > MIN_ANGLE && angleFn) || fabs(xAxisTilt) > 0.1;
  tiltAxisAngle *= RADIANS_PER_DEGREE;
  xAxisTilt *= RADIANS_PER_DEGREE;

  if (PipGetInteger("DefocusTol", &defocusTol))
    exitError("No DefocusTol specified");
  if (PipGetInteger("InterpolationWidth", &iWidth) || !iWidth)
    exitError("No InterpolationWidth specified, or 0 value entered");
  scaleInterp = iWidth > 0;
  if (!scaleInterp)
    iWidth = -iWidth;
  if (!PipGetFloat("MinimumZeroSpacing", &minInterZeroPixels))
    if (minInterZeroPixels < 3 || minInterZeroPixels > 30)
      exitError("MinimumZeroSpacing must be between 3 and 30");
  if (PipGetFloat("PixelSize", &pixelSize))
    exitError("No PixelSize specified");
  if (PipGetInteger("Voltage", &volt))
    exitError("Voltage is not specified");
  if (PipGetFloat("SphericalAberration", &cs))
    exitError("SphericalAberration is not specified");
  cs = B3DMAX(0.01, cs);
  if (PipGetFloat("AmplitudeContrast", &ampContrast))
    exitError("No AmplitudeContrast is specified");
  ierr = PipGetFloat("PhasePlateShift", &phasePlateShift);
  if (!PipGetFloat("PhaseShiftInDegrees", &phaseDeg)) {
    if (!ierr)
      exitError("You cannot enter phase shift in both degrees and radians");
    phasePlateShift = phaseDeg * RADIANS_PER_DEGREE;
  }
  PipGetFloat("CutOnFrequency", &cutOnEntered);
  PipGetFloat("ScaleByCtfPower", &scaleByPower);
  PipGetFloat("OffsetInZ", &offsetInZ);
  if (PipGetTwoIntegers("TotalViews", &startingTotal, &endingTotal))
    isSingleRun = true; // TotalViews is not specified;
  if (PipGetString("OutputFileName", &outFn))
    exitError("OutputFileName is not specified");
  PipGetString("BoundaryInfoFile", &boundFn);
  PipGetInteger("SkipCorrectedFlag", &skipDoneFlag);
  PipGetString("TransformFile", &xformFn);
  PipGetBoolean("InvertTiltAngles", &invertAngles);
  angleSign = invertAngles ? -1. : 1.;

  // Get GPU specification and flag its source
  gpuEnv = getenv("IMOD_USE_GPU");
  if (gpuEnv) {
    ifGpuByEnv = true;
    useGPU = atoi(gpuEnv);
  }
  if (!PipGetInteger("UseGPU", &useGPU))
    ifGpuByEnv = false;
  gpuEnv = getenv("IMOD_USE_GPU2");
  if (gpuEnv) {
    ifGpuByEnv = true;
    useGPU = atoi(gpuEnv);
  }
  PipGetTwoIntegers("ActionIfGPUFails", &actGpuFailOption, &actGpuFailEnviron);

  printf("stackFn = %s, angleFn=%s,  invertAngles=%d\n", stackFn, angleFn, invertAngles);
  printf("volt=%d Kv, interpolationWidth=%d pixels, defocusTol=%d nm \n",
         volt, iWidth, defocusTol);
  printf("pixelSize=%f nm, cs=%f mm, ampContrast=%f \n", pixelSize, cs, ampContrast);

  FILE *fpStack;
  if ((fpStack = iiFOpen(stackFn, "rb")) == 0)
    exitError("Could not open input file %s", stackFn);
  defocusList = readDefocusFile(defFn, defVersion, defFlags);
  if (!ilistSize(defocusList))
    exitError("The defocus file %s is non-existent or empty - did you save in ctfplotter?"
              , defFn);

  FILE *foutput;
  MrcHeader header;
  MrcHeader outHeader;
  int sliceMode;
  ImodImageFile *iiFile;
  bool parallelHDF = false;
  
  /* read header */
  if (mrc_head_read(fpStack, &header))
    exitError("Reading header of input file %s",  stackFn);
  if (!(skipDoneFlag & 1) && (header.imodFlags & MRC_FLAGS_CTF_CORRECTED))
    exitError("The input file stack header has the flag set that it has been CTF "
              "corrected");
  outHeader = header;
  maxWasEntered = PipGetInteger("MaximumStripWidth", &maxStripWidth) == 0;
  if (!maxWasEntered)
    maxStripWidth = header.nx / 2;
  doFullImages = maxStripWidth >= header.nx;
  if (doDiagonals && maxWasEntered && !doFullImages)
    exitError("You cannot enter a maximum strip width less than the image size with a "
              "non-zero X-axis tilt or tilt axis angle");
  if (doDiagonals)
    doFullImages = true;
  if (!doFullImages)
    maxStripWidth = 2 * (maxStripWidth / 2);
  if (useGPU >= 0 && !gpuAvailable(useGPU, &gpuMemory, debugMode)) {
    ierr = ifGpuByEnv ? actGpuFailEnviron : actGpuFailOption;
    if (ierr > 1)
      exitError("Use of a GPU was requested but none is available");
    printf("%sUse of a GPU was requested but none is available; using the CPU\n",
           ierr ? "MESSAGE: Ctfphaseflip - " : "");
    useGPU = -1;
  }

  // Allocate arrays for shifts and rotations and set to 0
  xShifts = B3DMALLOC(float, header.nz);
  rotations = B3DMALLOC(float, header.nz);
  if (!xShifts || !rotations)
    exitError("Allocating array for X shifts or rotations");
  for (ii = 0; ii < header.nz; ii++)
    xShifts[ii] = rotations[ii] = 0.;

  // If transforms, read them and store shifts and rotations
  if (xformFn) {
    float dum1, dum2, dum3, dum4, dum5, dum6, dum7;
    FILE *fpXF = fopen(xformFn, "r");
    if (!fpXF)
      exitError("Opening transform file %s", xformFn);
    for (ii = 0; ii < header.nz; ii++) {
      ierr = fgetline(fpXF, line, MAXLINE);
      if (ierr == -1)
        exitError("Reading transform file %s", xformFn);
      if (ierr == -2)
        exitError("End of file after reading %d transforms; not enough transforms in "
                  "file", ii);
      ierr = sscanf(line, "%f %f %f %f %f", &dum1, &dum2, &dum3, &dum4, &xShifts[ii]);
      if (ierr <= 0) {
        ii--;
        continue;
      }
      amatToRotmagstr(dum1, dum2, dum3, dum4, &rotations[ii], &dum5, &dum6, &dum7);
    }
  }

  startingView = 1;
  endingView = header.nz;
  PipGetTwoIntegers("StartingEndingViews", &startingView, &endingView);
  sliceMode = sliceModeIfReal(header.mode);
  if (sliceMode < 0)
    exitError("File mode is %d; only byte, short integer, or real allowed", header.mode);

  //The number of slices this run deals with;
  int currNz = endingView - startingView + 1;
  if (isSingleRun) {
    outHeader.nz = currNz;
    outHeader.mz = currNz;
  } else {
    outHeader.nz = endingTotal - startingTotal + 1;
    outHeader.mz = endingTotal - startingTotal + 1;

    // Determine if doing parallel HDF by output type setting for setup run, or
    // by simple test of file for real runs
    if (startingView == -1 && endingView == -1)
      parallelHDF =  b3dOutputFileType() == OUTPUT_TYPE_HDF;
    else
      parallelHDF = boundFn && iiTestIfHDF(outFn) > 0;
  }
  outHeader.zlen = header.zlen * outHeader.mz / header.mz;
  mrcInitOutputHeader(&outHeader);
  mrc_head_label(&outHeader, "ctfPhaseFlip: CTF correction with phase flipping only");
  setOrClearFlags((b3dUInt32 *)&outHeader.imodFlags, MRC_FLAGS_CTF_CORRECTED,
                  (skipDoneFlag & 2) ? 0 : 1);

  if ((startingView == -1 && endingView == -1) || isSingleRun) {
    if (!isSingleRun && b3dOutputFileType() == OUTPUT_TYPE_TIFF)
      exitError("Cannot do parallel writing to a TIFF output file");
    imodBackupFile(outFn);
    foutput = iiFOpen(outFn, "wb");
  } else if (!parallelHDF) {
    foutput = iiFOpen(outFn, "r+b");
    if (!foutput)
      exitError("fopen() failed to open %s", outFn);
  }


  // Starting run of parallel run: write header and exit
  if (startingView == -1 && endingView == -1 && !isSingleRun) {
    if (parallelHDF) {
      iiFile = iiLookupFileFromFP(foutput);
      iiSyncFromMrcHeader(iiFile, &outHeader);
      if (!iiFile || iiFile->file != IIFILE_HDF)
        exitError("Expected an HDF file but new file %s is not HDF", outFn);
      if (hdfWriteDummySection(iiFile, boundFn, 0))
        exitError("Writing dummy dataset for first section to HDF file");
    }
    if (mrc_head_write(foutput, &outHeader))
      exitError("Error when write out header");
    iiFClose(fpStack);
    iiFClose(foutput);
    return 0;
  }

  int err = parWrtInitialize(boundFn, header.nx * (parallelHDF ? -1 : 1), header.ny);
  if (err)
    exitError("Initializing parallel writing with boundary info file %s "
              "(error %d)", boundFn, err);
  if (parallelHDF) {
    parWrtProperties(&ii, &ierr, &ierr);
    if (b3dLockFile(ii))
      exitError("Failed to obtain initial lock on HDF file");
    foutput = iiFOpen(outFn, "r+b");
    if (!foutput) {
      b3dUnlockFile(ii);
      exitError("iiFOpen() failed to open %s", outFn);
    }
    iiFile = iiLookupFileFromFP(foutput);
    if (!iiFile || iiFile->file != IIFILE_HDF)
      exitError("Expected an HDF file but existing file %s is not HDF", outFn);
    if (parWrtRecloseHDF(iiFile, NULL))
      exitError("Closing or unlocking HDF file\n");
    foutput = (FILE *)iiFile;
  }

  int nx = header.nx;
  int nyfile = header.ny;
  int nz = header.nz;
  float currAngle;
  Islice *currSlice;
  int stripPixelNum, interPixelNum;
  int k, view, row, column, i, nxPad, fullXdim, widthForScaling, maxCopyPixels;
  int stripInd, fy, fyy,  fx, ny, yoff, xoff, curOffset, lastOffset;
  float WL, C1, C2, f2, ctf, freq_scalex, freq_scaley, fyComponent, freqScaleXsq;
  float waveAberration, ampAngle, attenFrac, firstZeroFreqSq, minAttenFreqSq;
  float phaseShift, denom, cosSum, gx, gy, pointDefocus, focusSum, focusDiff, cutonToUse;
  float sinAstig, cosAstig, stripFocus2, phaseFrac, phaseFracFactor, cutonAngstroms;
  float interZero, lastZero, curZero;
  bool doScaleByPower = scaleByPower > 0.;
  bool powerIsHalf = fabs(scaleByPower - 0.5) < 1.e-5;
  bool generalPower = !powerIsHalf && fabs(scaleByPower - 1.) > 1.e-5;
  bool haveAstig = (defFlags & DEF_FILE_HAS_ASTIG) != 0;
  bool havePhase = (defFlags & DEF_FILE_HAS_PHASE) != 0;
  bool haveCuton = (defFlags & DEF_FILE_HAS_CUT_ON) != 0;

  //Get the tilt angles and detected defocus for each slice;
  SavedDefocus *item;
  float *defocus = B3DMALLOC(float, nz);
  float *defocus2 = B3DMALLOC(float, nz);
  float *astigAngle = B3DMALLOC(float, nz);
  float *platePhase = B3DMALLOC(float, nz);
  float *cutonFreq = B3DMALLOC(float, nz);
  float *curFrac = B3DMALLOC(float, nx);
  float *lastFrac = B3DMALLOC(float, nx);
  double wallPrep = 0., wallFFT = 0., wallCorr = 0., wallInterp = 0., wallStart, curTime;
  if (!defocus || !defocus2 || !astigAngle || !platePhase || !cutonFreq || !curFrac ||
      !lastFrac)
    exitError("Allocating memory for defocus arrays");

  if (angleFn)
    tiltAngles = readTiltAngles(angleFn, nz, angleSign, minAngle, maxAngle);

  // Check the defocus list if there is more than one value
  if (ilistSize(defocusList) > 1 && 
      checkAndFixDefocusList(defocusList, tiltAngles, nz, defVersion))
    printf("WARNING: ctfphaseflip - View numbers in defocus file are not all "
           "consistent with the angular ranges\n");

  //sets to UNUSED_DEFOCUS since ctfplotter saves defocus with that value
  // when defocus is not computed.
  // Also detect if all angles are 0 and 
  ierr = 1;
  for (k = 0; k < nz; k++) {
    defocus[k] = defocus2[k] = astigAngle[k] = platePhase[k] = cutonFreq[k] =
      UNUSED_DEFOCUS;
    if (angleFn && fabs(tiltAngles[k]) > MIN_ANGLE)
      ierr = 0;
  }
  if (ierr) {
    doFullImages = true;
    maxStripWidth = nx;
    printf("All tilt angles are 0, doing single correction of each full image\n");
  }

  // Process the defocus values; slice numbers are now numbered from 0
  for (i = 0; i < ilistSize(defocusList); i++) {
    item = (SavedDefocus *)ilistItem(defocusList, i);
    k = (item->startingSlice + item->endingSlice) / 2;
    if (k < 0 || k >= nz)
      exitError("View numbers in defocus file are out of range");
    // They are already in microns
    defocus[k] = item->defocus;
    if (item->defocus2 != 0.) {
      defocus2[k] = item->defocus2;
      astigAngle[k] = item->astigAngle;
    }
    if (item->platePhase != 0.) {
      platePhase[k] = fabs(item->platePhase);
    }
    if (item->cutOnFreq > 0.) {
      cutonFreq[k] = item->cutOnFreq;
    }
    //printf("beginNum=%d endNum=%d k=%d defocus=%f\n", beginNum, endNum, k,
    //defocus[k]);
  }

  // interpolation for defocus, astig angle, and phase shifts can use the function
  interpolateTable(defocus, nz, false);
  if (haveAstig)
    interpolateTable(astigAngle, nz, true);
  if (havePhase)
    interpolateTable(platePhase, nz, false);
  if (haveCuton)
    interpolateTable(cutonFreq, nz, false);

  // But if there are any astigmatism entries, we need to fill in defocus values to 
  // interpolate the astigmatism amplitude between existing information for that
  if (haveAstig) {
    int first = -1, second = 0;
    for (k = 0; k < nz; k++) {
      if (astigAngle[k] != UNUSED_DEFOCUS) {
        astigAngle[k] += rotations[k];
        if (astigAngle[k] < -90)
          astigAngle[k] += 180;
        else if (astigAngle[k] > 90)
          astigAngle[k] -= 180.;
      }
      if (defocus2[k] == UNUSED_DEFOCUS)
        continue;
      second = k;
      if (first == -1) {
        for (row = 0; row < second; row++)
          defocus2[row] = defocus[row] * (defocus2[second] / defocus[second]);
      } else {
        for (row = first + 1; row < second; row++)
          defocus2[row] = defocus[row] * 
            ((row - first) * (defocus2[second] / defocus[second]) + (second - row) *
             (defocus2[first] / defocus[first])) / (float)(second - first);
      }
      first = second;
    }
    
    for (k = nz - 1; k >= 0; k--)
      if (defocus2[k] == UNUSED_DEFOCUS)
        defocus2[k] = defocus[k] * (defocus2[second] / defocus[second]);
      else
        break;
  }

  // Report the angles, after applying the offset if any
  offsetInZ *= 0.001 * pixelSize;
  focusDiff = offsetInZ;
  for (k = 0; k < nz; k++) {
    if (offsetInZ) {

      // The rationale for the sign is that the right side of an aligned image is both
      // lower in the scope and at a higher underfocus according to the man page on
      // invertAngles, and lower in the scope should be lower in tomogram.  Thus
      // negative Z offset in the tomogram should be higher underfocus
      // The vertical (z) distance between tilted planes is the hypotenuse of a triangle
      // with the side perpendicular to the planes (the entered offset) adjacent to 
      // the angle, hence the cosine here.
      if (tiltAngles)
        focusDiff = offsetInZ / cos(RADIANS_PER_DEGREE * tiltAngles[k]);
      defocus[k] -= focusDiff;
      if (haveAstig && astigAngle[k] != UNUSED_DEFOCUS)
        defocus2[k] -= focusDiff;
      if (!hushDefoci)
        printf("adjusted defocus[%d] = %f microns", k , defocus[k]);
    } else if (!hushDefoci)
      printf("defocus[%d] = %f microns", k , defocus[k]);
    if (haveAstig && !hushDefoci && astigAngle[k] != UNUSED_DEFOCUS)
      printf("   defocus2 = %f   angle = %.2f", defocus2[k], astigAngle[k]);
    if (havePhase && !hushDefoci)
      printf("   phase shift = %f", platePhase[k]);
    if (!hushDefoci)
      printf("\n");
  }

  WL = 12.41 / sqrt(volt * (volt + 1022.0)); //wavelength;
  C1 = MY_PI * WL;
  C2 = C1 * cs * 1000000.0 * WL * WL / 2.0;


  int stripDist[2];
  float *restoredArray, *strip;
  float *fullImage = NULL, *fullCopy1 = NULL, *fullCopy2 = NULL;
  double meanSum = 0.0;
  double amin = 0.1 * numeric_limits<double>::max();
  double amax = -amin;
  int niceLimit = niceFFTlimit();
  // if (useGPU >= 0)
    niceLimit = 5;

  Islice *outSlice;

  // Pad the extent in Y if necessary
  ny = niceFrame(nyfile, 2, niceLimit);
  yoff = (ny - nyfile) / 2;
  xoff = 0;
  int currK = 0;
  if (!isSingleRun)
    currK = startingView - startingTotal;
  if (doFullImages) {
    nxPad = niceFrame(nx, 2, niceLimit);
    fullXdim = nxPad + 2;
    xoff = (nxPad - nx) / 2;
  }
  fflush(stdout);

  for (view = startingView; view <= endingView; view++) {
    if (tiltAngles) {
      currAngle = tiltAngles[view - 1];
      printf("Slice %d, tilt angle is %.2f degrees. \n", view, currAngle);
    } else {
      currAngle = 0.0;
      printf("Slice %d, no angle is specified, set to 0.0\n", view);
    }

    if (defocus[view - 1] == UNUSED_DEFOCUS)
      exitError("specified defocus is wrong for slice %d", view);
    ampAngle = atan(ampContrast / sqrt(1. - ampContrast * ampContrast));

    // Assign phase shift from table or entered value
    if (havePhase)
      phaseShift = platePhase[view - 1];
    else
      phaseShift = phasePlateShift;
    if (phaseShift == UNUSED_DEFOCUS)
      exitError("specified phase shift is wrong for slice %d", view);

    // Assign cut-on frequency from table or entered value
    if (haveCuton)
      cutonToUse = cutonFreq[view - 1];
    else
      cutonToUse = cutOnEntered;
    if (cutonToUse == UNUSED_DEFOCUS)
      exitError("specified cut-on frequency is wrong for slice %d", view);
    phaseFrac = 1.;
    if (cutonToUse > 0.) {
      cutonAngstroms =  cutonToUse / 10.;
      phaseFracFactor = 1. / (1. - exp(-FREQ_FOR_PHASE / cutonToUse));
      printf("phaseFracFactor %f  coxpx2 %f\n", phaseFracFactor, cutonAngstroms);
    }

    currSlice = sliceCreate(nx, nyfile, sliceMode);
    outSlice = sliceCreate(nx, nyfile, SLICE_MODE_FLOAT);
    if (!currSlice || !outSlice)
      exitError("creating outslice or currSlice");
    restoredArray = outSlice->data.f;

    //startingView starts at 1, the API starts 0;
    if (mrc_read_slice(currSlice->data.b, fpStack, &header, view - 1, 'Z'))
      exitError("reading slice");

    //convert slice to floats
    if (sliceMode != SLICE_MODE_FLOAT)
      if (sliceNewMode(currSlice, SLICE_MODE_FLOAT) < 0)
        exitError("converting slice to float");

    angleIsZero = fabs(currAngle) <= MIN_ANGLE && (!maxWasEntered || doFullImages);
    currAngle = currAngle * MY_PI / 180.0;
    if (!angleIsZero) {
      stripPixelNum = 2 * ((int)(fabs(defocusTol / tan(currAngle)) / pixelSize) / 2);
      maxCopyPixels = stripPixelNum / 2;
    } else {
      maxCopyPixels = stripPixelNum = nx;
    }
    if (stripPixelNum > nx)
      stripPixelNum = nx;
    
    // Evaluate a minimum strip width that allows some criterion pixels between zeros
    for (i = 1; i < 100; i++) {
      curZero = getAnyOneZero(defocus[view - 1], phaseShift, i, ampContrast, cs, 
                              pixelSize, volt);
      if (i > 1) {
        interZero = curZero - lastZero;
        if (curZero > freqForInterZero)
          break;
      }
      lastZero = curZero;
    }
    minWidth = minInterZeroPixels * 2. / interZero;
    ACCUM_MAX(minWidth, minStripWidth);

    // If the defocus tolerance allows a width bigger than the basic amount, and no 
    // specific maximum width was entered, find the first maximum width that gives
    // sufficient zero shift across the whole width of a strip to eliminate rings at the
    // zeros
    maxWidth = maxStripWidth;
    if (angleIsZero)
      maxWidth = nx;
    else if (!maxWasEntered && stripPixelNum > minWidthToScaleInterp) {
      maxWidth = minWidthToScaleInterp;
      zeroShift = 0.5 * maxWidth * 
        firstZeroShift(defocus[view - 1], phaseShift, currAngle, maxWidth, ampContrast,
                       cs, pixelSize, volt);
      while (maxWidth < maxStripWidth && 
             zeroShift < (1. - fabs(tan(currAngle))) * minZeroShift) {
        maxWidth += 2;
        zeroShift = 0.5 * maxWidth * 
          firstZeroShift(defocus[view - 1], phaseShift, currAngle, maxWidth, ampContrast,
                         cs, pixelSize, volt);
      }
    }

    // Scale the spacing between strips if this makes width be bigger than basic amount
    widthForScaling = B3DMIN(stripPixelNum, maxWidth);
    interPixelNum = iWidth;
    if (angleIsZero)
      interPixelNum = nx;
    else if (iWidth > 1 && scaleInterp && widthForScaling > minWidthToScaleInterp)
      interPixelNum = B3DMIN((widthForScaling * iWidth) / minWidthToScaleInterp,
                             stripPixelNum / 2 - 1);

    // Now adjust either the min if max entered, or the max if it was not
    if (maxWasEntered)
      minWidth = B3DMIN(minWidth, maxWidth);
    else
      maxWidth = B3DMAX(minWidth, maxWidth);

    // Limit the strip width then get it to a nice size
    B3DCLAMP(stripPixelNum, minWidth, maxWidth);
    stripPixelNum = niceFrame(stripPixelNum, 2 , niceLimit);
    zeroShift = 0.5 * stripPixelNum * firstZeroShift
      (defocus[view - 1], phaseShift, currAngle, stripPixelNum, ampContrast, cs,
       pixelSize, volt);

    //interPixelNum must be less than stripPixelNum/2;
    if (interPixelNum >= stripPixelNum / 2 && !angleIsZero)
      exitError("Interpolation width is too high, must be less than %d", 
                minStripWidth / 2);

    // All of the above was needed to get dynamic interpolation width, now set for full
    if (doFullImages)
      stripPixelNum = nxPad;
    printf("stripPixelNum=%d interPixelNum=%d zeroShift=%f inter-strip shift=%f\n", 
           stripPixelNum, interPixelNum, zeroShift, firstZeroShift
           (defocus[view - 1], phaseShift, currAngle, interPixelNum, ampContrast, cs,
            pixelSize, volt));

    int stripXdim = stripPixelNum + 2;
    bool finished = false, doRightShiftedMid = false, inRightShiftedMid = false;
    bool viewHasAstig = haveAstig && astigAngle[view - 1] != UNUSED_DEFOCUS;
    float *curStrip, *lastStrip;
    int stripBegin, fullFirstMid, intervals, numLeftShiftedMid;
    int stripEnd = 0;
    int stripStride, halfStrip, stripMid;
    int effectiveNx;
    float effectiveXcen, xPixCenter, yPixCenter, viewAxisAngle, cornerDist1, cornerDist2;
    float sinViewAxis, cosViewAxis, lowLim, highLim, axisDist, curAxisDist, lastAxisDist;
    float curAxFrac, deltaZ, constantZaxis;

    effectiveNx = nx;
    effectiveXcen = nx / 2;

    // Do initial operations on GPU so that if this one fails, it can fall back
    if (useGPU >= 0 && gpuInitializeSlice(currSlice->data.f, nx, nyfile, stripXdim, 
                                          stripPixelNum, ny, doFullImages)) {
      ierr = ifGpuByEnv ? actGpuFailEnviron : actGpuFailOption;
      if (ierr > 1)
        exitError("Failure in initial call to GPU for view %d", view);
      printf("%sFailure in initial call to GPU for view %d; falling back to CPU\n",
             ierr ? "MESSAGE: Ctfphaseflip - " : "", view);
      useGPU = -1;
    }

    //Allocate 2 strips, even and odd strip;
    if (!doFullImages) {
      xoff = 0.;
      if (useGPU < 0)  {
        strip = B3DMALLOC(float, 2 * ny * stripXdim);
        if (!strip)
          exitError("Allocating arrays for strips");
      }
      if (!angleIsZero) {
        numLeftShiftedMid = ((stripPixelNum / 2 - maxCopyPixels) + interPixelNum - 1) /
          interPixelNum;
        ACCUM_MAX(numLeftShiftedMid, 0);
        while (stripPixelNum / 2 - numLeftShiftedMid * interPixelNum < 2)
          numLeftShiftedMid--;
      }
    } else {

      // Or taper and pad full image and take its FFT
      if (useGPU < 0)  {

        // First time this is encountered, allocate the arrays - do it this way
        // in case started on GPU and fell back
        if (!fullImage) {
          fullImage = B3DMALLOC(float, ny * fullXdim);
          fullCopy1 = B3DMALLOC(float, ny * fullXdim);
          fullCopy2 = B3DMALLOC(float, ny * fullXdim);
          if (!fullImage || !fullCopy1 || !fullCopy2)
            exitError("Allocating arrays for full images/FFTs");
        }
        START_TIMER;
        sliceTaperInPad(currSlice->data.f, SLICE_MODE_FLOAT, nx, 0, nx - 1, 0, nyfile - 1,
                        fullImage, fullXdim, nxPad, ny, 9, 9);
        ADD_TIME(wallPrep);
        todfftc(fullImage, nxPad, ny, 0);
        //utilDumpFFT(fullImage, nxPad, ny, "cpu-pad-fft", 0);
        ADD_TIME(wallFFT);
      }

      // Determine how to do full image along diagonals
      // Angles are already radians
      if (doDiagonals) {
        xPixCenter = nx / 2. - 0.5;
        yPixCenter = nyfile / 2. - 0.5;
        
        // Get the axis of constant Z before rotation by the tilt axis angle
        // Keep it in the 1st and 2nd quadrant, on both sides of tilt axis
        if (xAxisTilt != 0.) {
          constantZaxis = atan(sin(currAngle) / tan(xAxisTilt));
          if (constantZaxis < 0)
            constantZaxis += 180. * RADIANS_PER_DEGREE;
        } else
          constantZaxis = 90. * RADIANS_PER_DEGREE;

        // Add tilt axis angle to get actual axis in image
        viewAxisAngle = constantZaxis + tiltAxisAngle;
        sinViewAxis = sin(viewAxisAngle);
        cosViewAxis = cos(viewAxisAngle);
        cornerDist1 = fabs(-sinViewAxis * nx / 2. - cosViewAxis * nyfile / 2.);
        cornerDist2 = fabs(-sinViewAxis * nx / 2. + cosViewAxis * nyfile / 2.);
        effectiveNx = 2 * (int)ceil(B3DMAX(cornerDist1, cornerDist2));
        effectiveXcen = effectiveNx / 2. - 0.5;
        printf("Axis of constant Z: %.2f\n", viewAxisAngle / RADIANS_PER_DEGREE);
      }

      intervals = (effectiveNx - 4) / interPixelNum; 
      fullFirstMid = (effectiveNx - intervals * interPixelNum) / 2;
      if (angleIsZero && !doDiagonals)
        fullFirstMid = nx / 2;
      //printf("%d %d %d\n", intervals, fullFirstMid, interPixelNum);
    }

    // convert pixelSize to Angstroms from nm to scale frequency to 1/A
    freq_scalex = 1.0 / (pixelSize * 10.0 * stripPixelNum);
    freqScaleXsq = freq_scalex * freq_scalex;
    freq_scaley = 1.0 / (pixelSize * 10.0 * ny);
    stripInd = 0;
    while (!finished) {

      /*
       * Set up the strip variables based on the current index
       *
       stripBegin, stripEnd, stripMid are in coordinates of the original image
       stripMid is the column whose defocus will be corrected for, not necesary in the
       middle.  
       The offsets are added to the computation of coordinates in the strip 
       The strip is never padded in X so it has an offset of 0 initially, but has
       a negative offset for the interpolations on the right
       A full image IS padded so offset is positive
       halfStrip gets set here to make the offsets work out right so it is not always
       half the strip size
       When not doing full images one way or another, the full pattern is to copy
       pixels on the left up to maxCopyPixels, setting stripMid the end of that copy,
       then do multiple strips starting at 0 and with stripMid stepped up to half the
       strip width, 
       then switch to advancing the strip location each time with stripMid in middle, 
       then do multiple strips ending at the end if needed to advance stripMid to
       within maxCopyPixels of the end, then copy up to that many pixels
      */

      if (doFullImages) {

        // Full image, shifting stripMid
        stripMid = stripInd * interPixelNum + fullFirstMid; 
        finished = stripMid + interPixelNum > effectiveNx - 4;
        curStrip = (stripInd % 2) ? fullCopy2 : fullCopy1;
        stripStride = interPixelNum;
        halfStrip = stripMid + (stripInd == 0 ? 1 : 0);
      } else {

          // Also full image, just one operation
        if (angleIsZero) {
          stripBegin = 0;
          stripStride = interPixelNum;
          stripEnd = nx - 1;
          finished = true;
          xoff = (stripPixelNum - nx) / 2;
          stripMid = (stripBegin + stripEnd) / 2;
          halfStrip = stripPixelNum / 2;

          // Left side matching strips with different middles, offset is zero
        } else if (stripInd <= numLeftShiftedMid) {
          stripBegin = 0;
          stripEnd = stripPixelNum - 1;
          stripStride = interPixelNum;
          stripMid = stripPixelNum / 2 - (numLeftShiftedMid - stripInd) * interPixelNum;
          halfStrip = stripMid + (stripInd == 0 ? 1 : 0);

          // Classic region with stripMid in middle of strip that is shifted from last
        } else if ((stripInd - numLeftShiftedMid) * interPixelNum + stripPixelNum - 1 < 
                   nx) {
          stripBegin = (stripInd - numLeftShiftedMid) * interPixelNum;
          stripEnd = stripBegin + stripPixelNum - 1;
          stripStride = interPixelNum;
          stripMid = (stripBegin + stripEnd) / 2;
          halfStrip = stripPixelNum / 2;
          if (stripEnd == nx - 1) {
            doRightShiftedMid = true;
            finished = nx - stripMid <= maxCopyPixels || stripMid + interPixelNum > nx -2;
          }

          // In right side matching strips with different middles
        } else if (doRightShiftedMid) {
          stripStride = interPixelNum;
          stripMid += interPixelNum;
          halfStrip = stripMid;
          inRightShiftedMid = true;
          xoff = -stripBegin;
          finished = nx - stripMid <= maxCopyPixels || stripMid + interPixelNum > nx - 2;

          // Or set up last possible strip with unique stride from last
        } else {
          stripStride = (nx - stripPixelNum) - stripBegin;
          stripBegin = nx - stripPixelNum;
          stripEnd = nx - 1;
          stripMid = (stripBegin + stripEnd) / 2;
          halfStrip = stripPixelNum / 2;
          doRightShiftedMid = true;
          finished = nx - stripMid <= maxCopyPixels || stripMid + interPixelNum > nx - 2;
        }
        curStrip = strip + (stripInd % 2) * ny * stripXdim;
        if (debugMode > 1)
          printf("stripInd=%d stripBegin=%d stripEnd=%d   ", stripInd, stripBegin, 
                 stripEnd);
      }
      
      // Get a change in Z height and a defocus in nm;
      if (doDiagonals) {

        // For diagonals, the axis distance is based on the strip coordinate;
        // handle the special case where the constant Z axis is at 0, and apply
        // formula for axis distance to get change in Z height
        axisDist = (stripMid - (effectiveXcen + xShifts[view - 1])) * pixelSize;
        if (xAxisTilt != 0.  && (fabs(currAngle) <= MIN_ANGLE * RADIANS_PER_DEGREE ||
                                 fabs(constantZaxis) < 1.e-6)) {
          deltaZ = axisDist * tan(xAxisTilt);
        } else {
          deltaZ = axisDist * tan(currAngle) / sin(constantZaxis);
        }
      } else {

        // For regular tilt, the axis distance is negative x at positive tilt
        deltaZ = (nx / 2 + xShifts[view - 1] - stripMid) * pixelSize * tan(currAngle);
      }

      // The defocus is higher with negative deltaZ, so subtract it
      stripDefocus = defocus[view - 1] * 1000.0 - deltaZ;
      //printf("defocus is %6.1f\n", stripDefocus);
      
      //convert defocus to Angstroms
      stripDefocus *= 10.;
      pointDefocus = stripDefocus;

      // For astigmatism, get the other defocus and the sine and cosine of angle
      cosAstig = sinAstig = focusSum = focusDiff = 0.;
      if (viewHasAstig) {
        stripFocus2 = (defocus2[view - 1] * 1000.0 - deltaZ) * 10.;
        sinAstig = sin(astigAngle[view - 1] * RADIANS_PER_DEGREE);
        cosAstig = cos(astigAngle[view - 1] * RADIANS_PER_DEGREE);
        focusSum = 0.5 * (stripDefocus + stripFocus2);
        focusDiff = 0.5 * (stripDefocus - stripFocus2);
      }

      // Get the first zero and frequency to start attenuation at (if any) as a square
      firstZeroFreqSq = pow(getAnyOneZero(stripDefocus / 10000., phaseShift, 1,
                                          ampContrast, cs, pixelSize, volt) / 
                            (20. * pixelSize), 2.);
      minAttenFreqSq = attenStartFrac * firstZeroFreqSq;

      if (useGPU < 0) {

        // Copy full image transform or taper-pad strip to curStrip and transform it
        if (doFullImages) {
          memcpy(curStrip, fullImage, fullXdim * ny * sizeof(float));
        } else {
          START_TIMER;
          sliceTaperInPad(currSlice->data.f, SLICE_MODE_FLOAT, nx, stripBegin,
                          stripEnd, 0, nyfile - 1, curStrip,
                          stripPixelNum + 2, stripPixelNum, ny, 9, 9);
          ADD_TIME(wallPrep);
          todfftc(curStrip, stripPixelNum, ny, 0);
          ADD_TIME(wallFFT);
        }

        // Correct the CTF
        for (fy = 0; fy < ny; fy++) {
          fyy = fy;
          if (fy > ny / 2)
            fyy -= ny;
          gy = fyy * freq_scaley;
          fyComponent = gy * gy;
          for (fx = 0; fx < stripXdim / 2; fx++) {
            f2 = fx * fx * freqScaleXsq + fyComponent;
            if (viewHasAstig && (fx || fy)) {

              // The equation here is from Rohou and Grigorieff, 2015
              // 0.5 * (df1 + df2 + (df1 - df2) * cos(2 * (angle - axis)))
              // cos(2 * (angle - axis)) = 2 * cos(angle - axis)**2 - 1
              // cos(angle - axis) = cos(angle) * cos(axis) + sin(angle) * sin(axis)
              // cos(angle) = gx / sqrt(gx**2 + gy**2)
              // sin(angle) = gy / sqrt(gx**2 + gy**2)
              gx = fx * freq_scalex;
              denom = sqrt(gx * gx + fyComponent);
              cosSum = (cosAstig * gx + sinAstig * gy) / denom;
              pointDefocus = focusSum + focusDiff * (2. * cosSum * cosSum - 1.);
            }
            if (cutonToUse > 0.)
              phaseFrac = phaseFracFactor * (1. - exp(-sqrtf(f2) / cutonAngstroms));
            waveAberration = (C2 * f2 - C1 * pointDefocus) * f2 - phaseFrac * phaseShift;

            // Produce a positive ctf for consistency and so it can be used for scaling
            // (Here is the formal equation before simplifying)
            /*ctf = -(sqrt(1 - ampContrast * ampContrast)) * sin(waveAberration)
              + ampContrast * cos(waveAberration);*/
            ctf = -sin(waveAberration - ampAngle);
            if (doScaleByPower && f2 > minAttenFreqSq) {
              if (powerIsHalf)
                ctf = sqrt(fabs(ctf)) * (ctf >= 0. ? 1. : -1.);
              else if (generalPower) 
                ctf = pow(fabs(ctf), scaleByPower) * (ctf >= 0. ? 1. : -1.);
              if (f2 < firstZeroFreqSq) {
                attenFrac = (firstZeroFreqSq - f2) / 
                  ((1. - attenStartFrac) * firstZeroFreqSq);
                ctf = attenFrac + (1. - attenFrac) * ctf;
              }
              curStrip[fy * stripXdim + 2 * fx] *= ctf;
              curStrip[fy * stripXdim + 2 * fx + 1] *= ctf;
            } else if (ctf < 0) {
              curStrip[fy * stripXdim + 2 * fx] *= -1.;
              curStrip[fy * stripXdim + 2 * fx + 1] *= -1.;
            }
          }
        }
        ADD_TIME(wallCorr);
        //utilDumpFFT(curStrip, stripPixelNum, ny, "cpu-corr-fft", 0);

        //inverse FFT;
        todfftc(curStrip, stripPixelNum, ny, 1);
        ADD_TIME(wallFFT);

      } else {

        // On GPU, do extraction and transform or full copy, then correct the CTF
        if (gpuExtractAndTransform(stripInd, stripBegin, stripEnd, 9, 9))
          exitError("Calling gpuExtractAndTransform for view %d", view);
        if (gpuCorrectCTF(stripInd, freq_scalex, freq_scaley, pointDefocus, cosAstig,
                        sinAstig, focusSum, focusDiff, cutonAngstroms, phaseFracFactor,
                        phaseShift, ampAngle, C1, C2, scaleByPower, powerIsHalf,
                        generalPower, firstZeroFreqSq, attenStartFrac, minAttenFreqSq))
          exitError("Calling gpuApplyCTF for view %d", view);
      }

      /*
       * Put corrected data column(s) into the restored array
       */
      // The starting strip requires a copy of columns
      if (stripInd == 0) { 
        if (doDiagonals) {

          // Diagonals
          lowLim = -0.5 - effectiveXcen;
          highLim = halfStrip - 0.5 - effectiveXcen;
          COPY_DIAGONALS;
        } else {

          // Regular columns
          if (useGPU < 0) {
            for (row = 0; row < nyfile; row++)
              for (column = 0; column < halfStrip; column++)
                restoredArray[row * nx + column] =
                  curStrip[(row + yoff) * stripXdim + column + xoff];
          } else {
            if (gpuCopyColumns(stripInd, xoff, yoff, 0, halfStrip - 1))
              exitError("Calling gpuCopyColumns for columns %d to %d of view %d", 
                        0, halfStrip - 1, view);
          }
        }
        if (debugMode > 1)
          printf("column=1 ... %d stripMid %d  xoff %d\n", halfStrip, stripMid, xoff);
        finished = finished || (stripEnd == effectiveNx - 1 && !doRightShiftedMid);
      } else {

        // Single pixel from interpolation width 1 is just a copy
        if (interPixelNum == 1) {
          if (doDiagonals) {

            // Diagonals
            lowLim = stripMid - 0.5 - effectiveXcen;
            highLim = stripMid + 0.5 - effectiveXcen;
            COPY_DIAGONALS;
          } else {

            // Regular columns
            if (useGPU < 0) {
              for (row = 0; row < nyfile; row++) {
                restoredArray[row * nx + stripMid] = 
                  lastStrip[(row + yoff) * stripXdim + halfStrip + xoff];
              }
            } else {
              if (doFullImages || stripInd <= numLeftShiftedMid || inRightShiftedMid)
                curOffset = xoff;
              else
                curOffset = halfStrip - stripMid;
              if (gpuCopyColumns(stripInd + 1, curOffset, yoff, stripMid, stripMid))
                exitError("Calling gpuCopyColumns for column %d of view %d", stripMid, k);
            }
          }
          
          // Otherwise do the interpolation
        } else {
          if (doDiagonals) {

            // Diagonals
            lastAxisDist = (stripMid - stripStride) - effectiveXcen;
            curAxisDist = stripMid - effectiveXcen;
            if (useGPU < 0) {
              lowLim = lastAxisDist + 0.5;
              highLim = curAxisDist + 0.5;
              for (row = 0; row < nyfile; row++) {
                for (column = 0; column < nx; column++) {
                  axisDist = -sinViewAxis * (column - xPixCenter) + cosViewAxis * 
                    (row - yPixCenter);
                  if (axisDist >= lowLim && axisDist <= highLim) {
                    B3DCLAMP(axisDist, lastAxisDist, curAxisDist);
                    curAxFrac = (B3DMIN(axisDist, curAxisDist) - lastAxisDist) / 
                      stripStride;
                    restoredArray[row * nx + column] = 
                      curAxFrac * curStrip[(row + yoff) * stripXdim + xoff + column] +
                      (1. - curAxFrac) * lastStrip[(row + yoff) * stripXdim + xoff +
                                                   column];
                  }
                }
              }
            } else {
              if (gpuInterpDiagonals(stripInd, xoff, yoff, stripStride, sinViewAxis,
                                     cosViewAxis, lastAxisDist, curAxisDist))
                exitError("Calling gpuInterpDiagonals for view %d", view);
            }
          } else {
            
            // Regular columns: set up arrays of fractions and set offsets
            for (column = stripMid - stripStride + 1; column < stripMid + 1; column++) {
              stripDist[0] = column - stripMid + stripStride - 1;
              stripDist[1] = stripMid + 1 - column;
              curFrac[column] = stripDist[0] / (float)stripStride;
              lastFrac[column] = stripDist[1] / (float)stripStride;
            }
            if (doFullImages || stripInd <= numLeftShiftedMid || inRightShiftedMid) {
              curOffset = lastOffset = xoff;
            } else {
              curOffset = halfStrip - (stripMid + 1);
              lastOffset = halfStrip + stripStride - (stripMid + 1);
            }
          
            // Do the interpolation
            if (useGPU < 0) {
              for (row = 0; row < nyfile; row++) {
                for (column = stripMid - stripStride + 1; column < stripMid + 1;
                     column++) {
                  restoredArray[row * nx + column] = 
                    curFrac[column] * curStrip[(row + yoff) * stripXdim + curOffset +
                                               column]
                    + lastFrac[column] * lastStrip[(row + yoff) * stripXdim + lastOffset +
                                                   column];
                }
              }
            } else {
              if (gpuInterpolateColumns(stripInd, yoff, stripStride, stripMid, halfStrip,
                                        curOffset, lastOffset))
                exitError("Calling gpuInterpolateColumns for columns %d to %d of view %d",
                          stripMid - stripStride + 1, stripMid, k);
            }
          }
        }
        if (debugMode > 1) {
          printf("column=%d ... %d", stripMid- stripStride+1+1, stripMid+1);
          if (interPixelNum == 1)
            printf("  offset %d\n", useGPU < 1 ? halfStrip + xoff : curOffset);
          else
            printf(" curOff %d  lastOff %d\n", curOffset, lastOffset);
        }
        if (useGPU < 0)
          ADD_TIME(wallInterp);
      }

      // Do a copy from the last strip when finished is set
      if (finished) {
        if (doDiagonals) {

          // Diagonals
          lowLim = (stripMid + 1) - 0.5 - effectiveXcen;
          highLim = effectiveNx - 0.5 - effectiveXcen;
          COPY_DIAGONALS;
        } else {

          // Regular columns
          if (doFullImages || inRightShiftedMid)
            curOffset = xoff;
          else
            curOffset = halfStrip - stripMid - 1;
          if (useGPU < 0) {
            for (row = 0; row < nyfile; row++)
              for (column = stripMid + 1; column < nx; column++)
                restoredArray[row * nx + column]  =
                  curStrip[(row + yoff) * stripXdim + curOffset + column];
          } else {
            if (gpuCopyColumns(stripInd, curOffset, yoff, stripMid + 1, nx - 1))
              exitError("Calling gpuCopyColumns for columns %d to %d of view %d", 
                        stripMid + 1, nx - 1, k);
          }
        }
        if (debugMode > 1)
          printf("column=%d ... %d  curOff %d\n", stripMid+1+1, nx, curOffset);
      }

      //printf("stripInd=%d\n", stripInd);
      stripInd++;
      lastStrip = curStrip;
    }//while strip loop

    // Get result from GPU if any, process and write
    if (useGPU >= 0 && gpuReturnImage(restoredArray))
      exitError("Calling gpuReturnImage for view %d", k);

    if (!doFullImages && useGPU < 0)
      free(strip);
    sliceMMM(outSlice);
    if (outSlice->min < amin)
      amin = outSlice->min;
    if (outSlice->max > amax)
      amax = outSlice->max;
    meanSum += outSlice->mean;
    if (sliceMode != SLICE_MODE_FLOAT)
      if (sliceNewMode(outSlice, sliceMode) < 0)
        exitError("converting slice to original mode");

    if (parallelWriteSlice(outSlice->data.b, foutput, &outHeader, currK))
      exitError("Writing slice %d error", currK);
    currK++;
    sliceFree(currSlice);
    sliceFree(outSlice);
    fflush(stdout);
  }//k slice

  if (isSingleRun) {
    outHeader.amin = amin;
    outHeader.amax = amax;
    outHeader.amean = meanSum / (double)currNz;
    if (mrc_head_write(foutput, &outHeader))
      exitError("Writing slice header error");
  } else {//for collectmmm
    printf("min, max, mean, # pixels= %f  %f  %f %d \n",
           amin, amax, meanSum / (double)currNz, nx * ny * currNz);
  }
  if (debugMode) {
    if (useGPU >= 0) {
      gpuGetTimes(curTime, wallPrep, wallFFT, wallCorr, wallInterp);
      printf("Copy %.3f  ", curTime);
    }
    printf("Prep %.3f  FFT %.3f  correct %.3f  Interp %.3f\n", wallPrep, wallFFT, 
           wallCorr, wallInterp);
  }
  if (parallelHDF) {
    if (parWrtFlushBuffers(iiFile, &outHeader))
      exitError("Doing final write of buffer to output HDF file");
    iiDelete(iiFile);
    parWrtClose();
  } else
    iiFClose(foutput);
  iiFClose(fpStack);
  free(defocus);
  free(defocus2);
  free(astigAngle);
  free(platePhase);
  free(cutonFreq);
  free(curFrac);
  free(lastFrac);
  B3DFREE(fullImage);
  B3DFREE(fullCopy1);
  B3DFREE(fullCopy2);
}

/*
 * Computes the shift in the first zero across an extent in pixels, given focus in 
 * microns, tilt angle in radians and the basic parameters of the CTF.  
 * The shift is in Nyquist units
 */
static double firstZeroShift(double defocus, double phase, double angle, int extent,
                             double ampContrast, double Cs, double pixelSize,
                             double voltage)
{
  double firstZero, focus;
  focus = defocus - 0.5 * extent * tan(angle) * pixelSize / 1000.;
  firstZero = getAnyOneZero(focus, phase, 1, ampContrast, Cs, pixelSize, voltage);
  focus = defocus + 0.5 * extent * tan(angle) * pixelSize / 1000.;
  return fabs(firstZero - getAnyOneZero(focus, phase, 1, ampContrast, Cs, pixelSize,
                                        voltage));
}

/*
 * Computes the given zero position in Nyquist units given focus in microns and phase 
 * shift in radians, and the basic parameters of the CTF
 */
static double getAnyOneZero(double defocus, double phase, int zeroNum, double ampContrast,
                            double Cs, double pixelSize, double voltage)
{
  double delz, theta;
  double wavelength = 1.241 / sqrt(voltage * (voltage +  1022.0));
  double CsOne = sqrt(Cs * wavelength); // deltaZ=-deltaZ'/mCs1;  In microns
  double CsTwo = sqrt(sqrt(1000000.0 * Cs / wavelength)); //theta=theta'*mCs2;
  double ampAngle = 2. * atan(ampContrast / sqrt(1. - ampContrast * ampContrast)) / MY_PI;
  delz = defocus / CsOne;
  theta = sqrt(delz - sqrt(B3DMAX(0., delz * delz + ampAngle + 2. * phase / MY_PI - 
                                  2. * zeroNum)));
  return theta * pixelSize * 2.0 / (wavelength * CsTwo);
}


/*
 * Interpolates or extend values from for one item in the table to fill a whole array
 */
static void interpolateTable(float *defocus, int nz, bool ifAngles)
{
  int k, row;
  int first = -1, second = 0;
  float diff, frac;

  // Skip to the next measured value
  for (k = 0; k < nz; k++) {
    if (defocus[k] == UNUSED_DEFOCUS)
      continue;
    second = k;

    // Then copy into start of table or interpolate between last and this one
    if (first == -1) {
      for (row = 0; row < second; row++)
        defocus[row] = defocus[second];
    } else {
      for (row = first + 1; row < second; row++) {
        if (ifAngles) {
          diff = angleWithinLimits(defocus[second] - defocus[first], -90., 90.);
          frac = (row - first) / (float)(second - first);
          defocus[row] = angleWithinLimits(defocus[first] + frac * diff, -90., 90.);
        } else {
          defocus[row] = ((row - first) * defocus[second] + (second - row) *
                          defocus[first]) / (float)(second - first);
        }
      }
    }
    first = second;
  }

  // Then copy last into end of table
  for (k = nz - 1; k >= 0; k--)
    if (defocus[k] == UNUSED_DEFOCUS)
      defocus[k] = defocus[second];
    else
      break;
}
