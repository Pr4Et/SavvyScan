/*
 *  alignframes - align movie frames and stack multiple frame files
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 2016-2019 by the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 *  $Id$
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <string>
#include "b3dutil.h"
#include "autodoc.h"
#include "cppdefs.h"
#include "iimage.h"
#include "parse_params.h"
#include "framealign.h"
#include "frameutil.h"
#ifdef TEST_SHRMEM
#include "ShrMemClient.h"
#endif

#define MAX_BINNINGS 6 
#define SIG2_ROUND_FAC 10000.
#define MAX_LINE 600
#define MAX_READ_THREADS 16
#define FRAME_DOSE_KEY "FrameDosesAndNumber"

class AliFrame
{
public:
  AliFrame() {};

  void main( int argc, char *argv[]);
 private:

  // Methods
  FILE *openAndReadHeader(const char *filename, MrcHeader *head, 
                          const char *descrip, bool testMode = false);
  void checkInputFile(const char *filename, MrcHeader *head, int nx, int ny, 
                      int combine);
  void addToSumBuffer(void *readBuf, int inMode, void *sumBuf, int &useMode, 
                      int nxy);
  void extractFileTail(const char *filename, std::string &sstr);
  int adjustTitleBinning(const char *title, std::string &sstr, float relBinning,
                         bool printChange);
  void minMaxSetSize(int basicSize, int numFrames, int &minSize, int &maxSize);
  int analyzeExtraHeader(FILE *inFP, MrcHeader *head, int startFrame, 
                         int endFrame, bool axisPixOnly, float **buffer,
                         int &bufSize, FloatVec &tilts, IntVec &setStarts,
                         int &minSet, int &maxSet, float &axisAngle, float &pixSize);
  void rotateFlipGainReference(float *ref, int &nxGain, int &nyGain, int rotateFlip);
  int openMdocFile(const char *filename, int &numSect, int &adocType);
  void expandFrameDosesNumbers(std::string &line, float *tempArr, int maxFrames,
                              int &totalFrames, float &totalDose, float *frameDoses);
  void frameGroupLimits(int numFetch, int nzAlign, int group, int &groupStart, 
                      int &groupEnd, int zStart, int zDir, int &izLow, int &izHigh);
};

#ifdef TEST_SHRMEM
static void defectFileToString(char *name, std::string &outStr);
static ShrMemClient sSMC;
#endif

// Static class instance of FrameAlign object
static FrameAlign sFA;

int main( int argc, char *argv[])
{
  AliFrame ali;
  ali.main(argc, argv);
  exit(0);
}

// A BIG MAIN METHOD
void AliFrame::main( int argc, char *argv[])
{
  const char *progname = imodProgName(argv[0]);
  char *filename, *tempname, *xfExt = NULL;
  std::vector<char *> inFiles;
  std::string xfName, sstr;
  std::vector<float> tiltAngles;
  int binsToTest[MAX_BINNINGS], numFiltTests[MAX_BINNINGS];
  float varyRadius2[MAX_BINNINGS][MAX_FILTERS];
  float varySigma2[MAX_BINNINGS][MAX_FILTERS];
  int numTimesBest[MAX_BINNINGS][MAX_FILTERS];
  char title[MRC_LABEL_SIZE + 1];
  char inLine[MAX_LINE + 1];
  ImodImageFile *fileCopies[MAX_READ_THREADS];
  
  char *framePath = NULL;
  char *extraName, *listName = NULL;
  char *stackName = NULL, *mdocName = NULL, *tiltName = NULL, *gainName = NULL;
  char *savedListName = NULL;
  char *outNames[2] = {NULL, NULL};
  char *openTsNames[2] = {NULL, NULL};
  FILE *inFP;
  FILE *outFPs[2];
  FILE *extraFP, *stackFP, *frcFP = NULL, *plotFP = NULL, *fileListFP = NULL;
  unsigned char *readBuf, *sumBuf, *useBuf;
  int needAlloc, bufAllocSize = 0;
  MrcHeader inHead, mainHead, gainHead, darkHead, stackHead, unwgtHead;
  MrcHeader *outHeads[2] = {&mainHead, &unwgtHead}, *headPtr;
  Islice *gainSlice = NULL;
  Islice *darkSlice = NULL;
  float *summed, *rotSum = NULL, *useSum, *unwgtSum = NULL;
  CameraDefects defects;
  int camSizeX = 0, camSizeY = 0;  // Camera size in X and Y represented in table
  int ignoreZvalue = 0;
  int testMode = 0;
  int defaultBinnings[] = {2, 3, 4, 6, 8};
  int targetAliSize = 1250;
  IloadInfo li;
  bool doRobust, sumInOnePass, copyShifts, wasGoodEnough, namesFromMdoc = false;
  bool gettingFRC, sumWithAlign, useBlockGroup = false, endReached = false;
  float taperFrac = 0.1;
  float scale = 1., totalScale = 0., defaultByteScale = 30.;
  float kFactor = 4.5;
  int shiftLimit = 20;
  int antiFiltType = 4;
  int hybridShifts = 0;
  int sumBin = 1;
  float sizeDiffCrit = 0.1;
  float maxMaxWeight = 0.1;
  float goodEnough = 0.;
  int refineAtEnd = 0;
  int rotationFlip = 0;
  int sumRotationFlip = 0;
  int groupSize = 1;
  int combineFiles = 0;
  int breakSetSize = 0;
  float refRadius2 = 0., refSigma2;
  int skipChecks = 0;
  int adjustMdoc = 0;
  int splineSmooth = 1;
  int minNumForSpline = 20;
  int trimCrit = 10;
  float iterCrit = 0.1;
  int groupRefine = 0;
  float truncLimit = 0.;
  int useGPU = -1;
  int gpuFlags = 0;
  int deferSum = 0;
  float gpuFracMem = 0.85;
  int numOutFiles = 1;
  int ifDebug = 0, corDefBinning = 1, scaleDefects = 0, maxDataSize = 0;
  int nx, ny, nz, alignBin, ind, nxSum, nySum, ix, iy, itest, faBestFilt, doSpline;
  int nxStack, nyStack, stackMode, numSect, adocInd, relXbin, relYbin, numAVAuse;
  int numSingleFiles, startCombine, endCombine, adocType, alignBinIn, dataSize;
  size_t tind;
  float xScale, yScale, zScale, relBinning, error, minError, minMean, halfCross;
  float quartCross, eighthCross, halfNyq, gpuMemory, needed, needForAli, gpuUsableMem;
  float needForPreOps;
  float memLimits[2], gpuMemLimit = 0., memoryLimit = 12.;
  float fullTaperFrac = 0.02;

  // framealign uses fullTaperFrac as the padding fraction so a default trimming by the
  // same amount will keep the padded align within the original size, good if it is 4K
  float trimFrac = fullTaperFrac;   
  double diff, minDiff, wallStart, wallRead = 0.;
  float imagesBinned = -1.;
  int numAVAinput = 7;
  int minFractionalAVA = 7;
  int reverse = 0;
  int startFrame = -1, endFrame = -1;
  bool warnedTwoPass = false;
  float frcDeltaR = 0.005;
  float ringCorrs[510];
  float radius1 = 0., radius2 = 0.06, sigma1 = 0.03, sigma2 = 0.0086;
  float *xShifts, *yShifts, *bestXshifts, *bestYshifts, *rawXshifts, *rawYshifts;
  float *bestXraw, *bestYraw;
  int alsumBinEntered, minBinningToTest = 100;
  int ierr, iz, zStart, zEnd, zDir, numOptArgs, numNonOptArgs, outMode, fullDataSize;
  int nxOut, nyOut, numInByOpt, numInFiles, ifile, maxNumZ = 0;
  int numBinTests, nxGain, nyGain, numVaries, indBestBin, indBestFilt, slideGrpSize;
  int useInd, summingMode, numTestLoops, useStart, useEnd, indBinUse, indFiltUse;
  int numFiltUse, numDone, numFetch, filt, inDataSize, stackBin, numHoldFull, outNum;
  float sumPadSize, alignPadSize, fullPadSize, needForGPUsum, pixTemp, physMem;
  int nzAlign, groupEnd, groupStart, izLow, izHigh, useMode, group, blockGrpSize;
  int startAssess = -1, endAssess = -1;
  int numFrameUse, numSets, minSet, maxSet, fileHasTilts, maxReadThreads = 1;
  int extraHasTilts = 0, extraHasGainRef = 0, extraHasAxisPix = 0, numReadThreads;
  int numAllSets, minSetSize, maxSetSize, superFile, setInFile, izRead, numAllVsAll;
  float fileAxis, extraAxis, filePix, extraPixSize, optionPixSize = 0.;
  bool hasExtra, enteredScale, enteredMode, needPreprocess, parallelRead = false;
  int useShrMem = 0;
  float stackBinX, stackBinY;
  FloatVec extraTilts;
  IntVec setStarts, savedFrames;
  float *extraBuf = NULL;
  int extraBufSize = 0;
  float resMean[MAX_FILTERS + 1], resSD[MAX_FILTERS + 1], meanResMax[MAX_FILTERS + 1];
  float maxResMax[MAX_FILTERS + 1], meanRawMax[MAX_FILTERS + 1];
  float maxRawMax[MAX_FILTERS + 1], smoothDist[MAX_FILTERS + 1], rawDist[MAX_FILTERS + 1];
  float  tmin, tmax, tmean, alreadyScaledBy;

  // Dose weighting variables
  float totalDose = 0., initialDose = 0., doseAfac = 0., doseBfac = 0., doseCfac = 0.;
  float doseScaling = 1., critDoseScale200KV = 0.8;
  float priorTemp, doseTemp, sumOfDoses;
  int doseFileType = -1, voltage = 300, numBidir = 0, doseAccumulates = -1;
  int numMdocSect, sumOfFrames, maxFrameDoses = 0;
  char *doseName = NULL;
  FloatVec doseFromMdoc, priorFromMdoc, tempVal1, totalDoseVec, priorDoseVec, frameDoses;
  std::vector<std::string> frameDoseLines, mdocFileTails;
  std::vector<char *> fullFramePaths;
  std::string fixedFrameDoses;
  FloatVec reweightOnes;
  float *reweightFilt = NULL;
  std::string defectString;
  
  IntVec izPiece;

  // Fallbacks from    ../manpages/autodoc2man 2 1 alignframes
  int numOptions = 72;
  const char *options[] = {
    "input:InputFile:FNM:", "output:OutputImageFile:FN:", "list:ListOfInputFiles:FN:",
    "break:BreakFramesIntoSets:I:", "saved:SavedFrameListFile:FN:",
    "skip:SkipFileChecks:B:", "stack:CorrespondingStack:FN:", "mdoc:MetadataFile:FN:",
    "path:PathToFramesInMdoc:CH:", "ignore:IgnoreZvaluesInMdoc:B:",
    "adjust:AdjustAndWriteMdoc:B:", "pixel:PixelSize:F:",
    "binning:AlignAndSumBinning:IP:", "frames:StartingEndingFrames:IP:",
    "mode:ModeToOutput:I:", "scale:ScalingOfSum:F:", "total:TotalScalingOfData:F:",
    "rfsum:SumRotationAndFlip:I:", "tilt:TiltAngleFile:FN:",
    "xfext:TransformExtension:CH:", "frc:FRCOutputFile:FN:",
    "ring:RingSpacingForFRC:F:", "plottable:PlottableShiftFile:FN:",
    "nosum:NoSumsOutput:B:", "gain:GainReferenceFile:FN:",
    "rotation:RotationAndFlip:I:", "dark:DarkReferenceFile:FN:",
    "defect:CameraDefectFile:FN:", "double:DoubleDefectCoords:B:",
    "imagebinned:ImagesAreBinned:F:", "truncate:TruncateAbove:FP:",
    "pair:PairwiseFrames:I:", "reverse:ReverseOrder:B:", "shift:ShiftLimit:I:",
    "group:GroupSize:I:", "radius2:FilterRadius2:F:", "vary:VaryFilter:FAM:",
    "hybrid:UseHybridShifts:B:", "refine:RefineAlignment:I:",
    "rgroup:RefineWithGroupSums:B:", "stop:StopIterationsAtShift:F:",
    "rrad2:RefineRadius2:F:", "smooth:MinForSplineSmoothing:I:", "gpu:UseGPU:I:",
    "memory:MemoryLimitGB:FA:", "dtype:TypeOfDoseFile:I:",
    "dfile:DoseWeightingFile:FN:", "dtotal:FixedTotalDose:F:",
    "dframe:FixedFrameDoses:F:", "dprior:InitialPriorDose:F:",
    "bidir:BidirectionalNumViews:I:", "accum:DoseAccumulates:I:",
    "normalize:NormalizeDoseWeighting:B:", "volt:Voltage:I:",
    "optimal:OptimalDoseScaling:F:", "critical:CriticalDoseFactors:FT:",
    "unweight:UnweightedOutputFile:FN:", "test:TestBinnings:IA:",
    "assess:AssessWithFrames:IP:", "good:GoodEnoughError:F:",
    "weight:MaxResidualWeight:F:", "trim:TrimFraction:F:", "taper:TaperFraction:F:",
    "antialias:AntialiasFilter:I:", "radius1:FilterRadius1:F:",
    "sigma1:FilterSigma1:F:", "sigma2:FilterSigma2:F:", "kfactor:KFactorForFits:F:",
    "debug:DebugOutput:I:", "flags:FlagsForGPU:I:", "shrmem:ShrMemTest:B:",
    "help:usage:B:"};

  // Startup with fallback
  PipReadOrParseOptions(argc, argv, options, numOptions, progname, 
                        2, 1, 1, &numOptArgs, &numNonOptArgs, imodUsageHeader);

  // Get output file and number of input files
  PipGetBoolean("NoSumsOutput", &testMode);
  gettingFRC = !testMode;
  PipNumberOfEntries("InputFile", &numInByOpt);
  numInFiles = numInByOpt + numNonOptArgs;
  if (PipGetString("OutputImageFile", &outNames[0])) {
    if (!testMode) {
      if (!numNonOptArgs)
        exitError("No output file specified");
      numInFiles--;
      PipGetNonOptionArg(numNonOptArgs - 1, &outNames[0]);
    }
  } else if (testMode) {
    exitError("No output file should be specified when not making sums");
  }

  if (PipGetString("ListOfInputFiles", &listName) == 0) {
    if (numInFiles)
      exitError("You cannot enter input files as arguments with the -list option");
    fileListFP = fopen(listName, "r");
    if (!fileListFP)
      exitError("Could not open list of input files, %s", listName);
    numInFiles = 2000000000;
  }
  PipGetInteger("BreakFramesIntoSets", &breakSetSize);
  if (breakSetSize != 0 && breakSetSize < 2)
      exitError("The entry for -break must be at least 2");
  PipGetBoolean("SkipFileChecks", &skipChecks);
  if (skipChecks)
    iiAssumeDMfileMatches(1);

  // Get a list of frames saved from SEMCCD in a single exposure tilt series 
  PipGetString("SavedFrameListFile", &savedListName);
  if (savedListName) {
    if (numInFiles != 1)
      exitError("There must be only a single input file with a saved frame list");
    if (breakSetSize)
      exitError("You cannot use the -break option with a saved frame list");
    inFP = fopen(savedListName, "r");
    if (!inFP)
      exitError("Opening frame list file %s", savedListName);
    while (1) {
      ierr = fgetline(inFP, title, MRC_LABEL_SIZE);
      if (!ierr)
        continue;
      if (ierr == -2)
        break;
      if (ierr == -1)
        exitError("Reading frame list file %s", tiltName);
      savedFrames.push_back(atoi(title));
      if (ierr < 0)
        break;
    }
    if (savedFrames.size() < 10)
      exitError("There are only %d numbers in the frame list file %s", 
                savedFrames.size(), savedListName);

    // Analyze it now so the number is known if tilt angles come in
    // Add non-existent set at end for ease of use
    setStarts.push_back(0);
    for (ind = 1; ind < savedFrames.size(); ind++)
      if (savedFrames[ind] != savedFrames[ind - 1] + 1)
        setStarts.push_back(ind);
    setStarts.push_back(ind);
  }
  
  // Find out what auxiliary files are being used
  PipGetString("MetadataFile", &mdocName);
  PipGetString("CorrespondingStack", &stackName);
  if (numInFiles > 0 && mdocName && stackName)
    exitError("You cannot enter -mdoc with -stack; the mdoc would not be used");
  if (!numInFiles && !mdocName)
    exitError("Input file(s) must be specified with arguments, an mdoc file, or a "
              "list file");

  // Open the mdoc now and set flag to get names from it if necessary
  if (mdocName) {
    adocInd = openMdocFile(mdocName, numSect, adocType);
    if (AdocGetInteger(ADOC_GLOBAL_NAME, 0, "DataMode", &stackMode) || 
        AdocGetTwoIntegers(ADOC_GLOBAL_NAME, 0, "ImageSize", &nxStack, &nyStack))
      exitError("Getting data mode or image size from mdoc file");

    namesFromMdoc = !numInFiles;
    if (namesFromMdoc) {
      numInFiles = numSect;
      PipGetString("PathToFramesInMdoc", &framePath);
    }
    PipGetBoolean("IgnoreZvaluesInMdoc", &ignoreZvalue);
    if (!testMode)
      PipGetBoolean("AdjustAndWriteMdoc", &adjustMdoc);
  }
  
  if (PipGetTwoIntegers("StartingEndingFrames", &startFrame, &endFrame) == 0) {
    if (startFrame <= 0 || endFrame < startFrame)
      exitError("Values for starting and ending frames are out of range");
    if (savedListName)
      exitError("You cannot use -frame with a saved frame list file");
  }
  if (PipGetTwoIntegers("AssessWithFrames", &startAssess, &endAssess) == 0 &&
      (startAssess <= 0 || endAssess < startAssess))
    exitError("Values for starting and ending frames to use for assessment are out of "
              "range");
  if ((breakSetSize > 0 || savedListName) && startAssess >= 0)
    exitError("You cannot enter -assess when breaking frames into sets to sum");
  if (breakSetSize > 0 && namesFromMdoc)
    exitError("You cannot break frames into sets when input filenames come from an mdoc "
              "file");

  PipGetFloat("TotalScalingOfData", &totalScale);
  enteredScale = PipGetFloat("ScalingOfSum", &scale) == 0;
  if (enteredScale && totalScale > 0.)
    exitError("You cannot enter both -scale and -total");
  alreadyScaledBy = 1.;
  enteredMode = PipGetInteger("ModeToOutput", &outMode) == 0;
  if (enteredMode && sliceModeIfReal(outMode) < 0)
    exitError("Output mode of %d is not allowed", outMode);
  PipGetString("GainReferenceFile", &gainName);
  PipGetInteger("RotationAndFlip", &rotationFlip);
  PipGetString("TiltAngleFile", &tiltName);
  PipGetString("CorrespondingStack", &stackName);
  PipGetInteger("SumRotationAndFlip", &sumRotationFlip);

  // Get dose-weighting related options
  filename = NULL;
  if (PipGetFloat("FixedTotalDose", &totalDose) + 
      PipGetInteger("TypeOfDoseFile", &doseFileType) + PipGetString("FixedFrameDoses",
                                                                    &filename) < 2)
    exitError("You can enter only one of the dose weighting options -dtype, -dtotal,"
              " or -dframe");
  
  if (doseFileType > 0 || totalDose > 0 || filename) {
    if (testMode)
      exitError("You cannot enter dose weighting options with no summing");
    PipGetFloat("InitialPriorDose", &initialDose);
    PipGetInteger("Voltage", &voltage);
    PipGetInteger("BidirectionalNumViews", &numBidir);
    PipGetFloat("OptimalDoseScaling", &doseScaling);
    PipGetThreeFloats("CriticalDoseFactors", &doseAfac, &doseBfac, &doseCfac);
    PipGetInteger("DoseAccumulates", &doseAccumulates);
    ind = 0;
    PipGetBoolean("NormalizeDoseWeighting", &ind);
    if (ind) {
      reweightOnes.resize(9000, 1.);
      reweightFilt = &reweightOnes[0];
    }

    // Incorporate voltage info into scaling factor
    if (voltage != 300) {
      if (voltage != 200)
        exitError("Voltage must be either 200 or 300");
      doseScaling *= critDoseScale200KV;
    }
    if (doseAccumulates < 0)
      doseAccumulates = (stackName || mdocName || savedListName) ? 1 : 0;
    if (filename) {
      fixedFrameDoses = filename;
      free(filename);
    }

    // Get dose file of various kinds
    if (doseFileType > 0) {
      PipGetString("DoseWeightingFile", &doseName);
      if (doseFileType == 4) {

        // Mdoc file, either existing one or specified one
        if (doseName && mdocName)
          exitError("You cannot enter a dose weighting mdoc file name if you also enter "
                    "the -mdoc option");
        if (mdocName)
          numMdocSect = numSect;
        else
          adocInd = openMdocFile(doseName, numMdocSect, adocType);

        // Set up vectors to call for doses from it
        doseFromMdoc.resize(numMdocSect);
        priorFromMdoc.resize(numMdocSect);
        izPiece.resize(numMdocSect);
        frameDoseLines.resize(numInFiles);
        for (iz = 0; iz < numMdocSect; iz++)
          izPiece[iz] = iz;
        ierr = getMetadataWeightingDoses(adocInd, adocType, numMdocSect, &izPiece[0],
                                         numBidir, &priorFromMdoc[0], &doseFromMdoc[0]);
        if (ierr == 1)
          exitError("Problems occurred accessing data in the mdoc file");
        if (ierr == 2)
          exitError("The dose information in the mdoc file was not usable");
      } else {

        // Other text files, the name must be provided; open the file and read lines
        if (!doseName)
          exitError("You must enter a dose weighting file also");
        inFP = fopen(doseName, "r");
        if (!inFP)
          exitError("Opening dose file %s", doseName);
        while (1) {
          ierr = fgetline(inFP, inLine, MAX_LINE);
          if (!ierr)
            continue;
          if (ierr == -2)
            break;
          if (ierr == -1)
            exitError("Reading dose file %s", doseName);

          // A single dose value for type 1, a line for type 4, or prior and another value
          if (doseFileType == 1) {
            totalDoseVec.push_back(atof(inLine));
            priorDoseVec.push_back(0.);
          } else if (doseFileType > 4) {
            frameDoseLines.push_back(std::string(inLine));
          } else {
            sscanf(inLine, "%f %f", &priorTemp, &doseTemp);
            if (doseFileType == 3)
              doseTemp -= priorTemp;
            totalDoseVec.push_back(doseTemp);
            priorDoseVec.push_back(priorTemp);
          }
          if (ierr < 0)
            break;
        }
      }
      B3DFREE(doseName);
    }

    // Get option for unweighted output also
    if (!PipGetString("UnweightedOutputFile", &outNames[1]))
      numOutFiles = 2;
  }

  // Open and read header of every input file before starting
  numAllSets = 0;
  minSetSize= 0;
  for (ind = 0; ind < numInFiles; ind++) {
    if (namesFromMdoc) {

      // Get the name from the mdoc file and extract it from the path
      if (ignoreZvalue)
        iz = ind;
      else
        iz = AdocLookupByNameValue(ADOC_ZVALUE_NAME, ind);
      if (iz < 0)
        exitError("Looking up section with Z value %d in mdoc file", ind);
      if (AdocGetString(ADOC_ZVALUE_NAME, iz, "SubFramePath", &filename))
        exitError("Getting SubFramePath for %s %d in mdoc file", ind, 
                  ignoreZvalue ? "section" : "Z value");
      extractFileTail(filename, sstr);
      free(filename);
      if (framePath) {
        sstr.insert(0, "/");
        sstr.insert(0, framePath);
      }
      filename = strdup(sstr.c_str());
      if (!filename)
        exitError("Duplicating filename from mdoc file");

      // Try to get a frame dose line from this section, store the dose data
      if (doseFileType == 4) {
        ierr = AdocGetString(ADOC_ZVALUE_NAME, iz, FRAME_DOSE_KEY, &tempname);
        if (ierr < 0)
          exitError("Trying to access "FRAME_DOSE_KEY" for %s %d in mdoc file", ind,
                    ignoreZvalue ? "section" : "Z value");
        if (!ierr) {
          frameDoseLines[ind] = tempname;
          free(tempname);
        }
        totalDoseVec.push_back(doseFromMdoc[iz]);
        priorDoseVec.push_back(doseAccumulates > 0 ? priorFromMdoc[iz] : 0.);
      }
          
    } else if (fileListFP) {

      // Or get name from file
      if (endReached)
        break;
      iz = fgetline(fileListFP, inLine, MAX_LINE);
      if (!iz)
        continue;
      if (iz == -2)
        break;
      if (iz == -1)
        exitError("Reading line %d of list of input files %s", ind + 1, listName);
      if (iz < 0)
        endReached = true;
      filename = strdup(inLine);
      if (!filename)
        exitError("Duplicating filename from list of input files");

    } else if (ind < numInByOpt) {

      // Or get arguments
      PipGetString("InputFile", &filename);
    } else {
      PipGetNonOptionArg(ind - numInByOpt, &filename);
    }

    // Save filename and check file
    inFiles.push_back(filename);
    if (!ind || !skipChecks) {
      inFP = openAndReadHeader(filename, &inHead, "input image", 
                             testMode && ind == numInFiles - 1);
      checkInputFile(filename, &inHead, ind ? nx : 0, ny, combineFiles);
      dataSizeForMode(inHead.mode, &dataSize, &iz);
      ACCUM_MAX(maxDataSize, dataSize);
    }
    
    // Set size and mode from first file.  Default to not do bytes as output
    if (!ind) {
      if (!enteredMode)
        outMode = inHead.mode == MRC_MODE_BYTE ? MRC_MODE_SHORT : inHead.mode;
      nx = inHead.nx;
      ny = inHead.ny;
      if (savedListName && inHead.nz != savedFrames.size())
        exitError("The number of frames (%d) does not match the number of entries in "
                  "the saved frame list (%d)", inHead.nz, savedFrames.size());

      mainHead = inHead;
      unwgtHead = inHead;
      mrc_get_scale(&inHead, &xScale, &yScale, &zScale);

      // Also get the rotation/flip if needed
      if (rotationFlip < 0 || totalScale > 0.) {
        for (ix = 0; ix < inHead.nlabl; ix++) {
          if (rotationFlip < 0) {
            extraName = strstr(inHead.labels[ix], " r/f ");
            if (extraName)
              rotationFlip = atoi(extraName + 4);
          }
          extraName = strstr(inHead.labels[ix], ", scaled by");
          if (extraName)
            alreadyScaledBy = atof(extraName + 11);
        }
        if (rotationFlip < 0)
          exitError("Cannot find r/f entry in header of first input file");
      }

      // And commit to combining files if one frame and breaking into sets, and disallow
      // frame subsets
      if (breakSetSize > 0 && inHead.nz == 1)
        combineFiles = breakSetSize;
      if (combineFiles > 0 && (startFrame >= 0 || startAssess >= 0))
        exitError("You cannot enter -frames when combining single-frame files");
      if (!combineFiles && skipChecks)
        exitError("You cannot skip file checks unless combining single-frame files");

      // Determine if reading TIFF and if so, get number of threads to use
      if (!combineFiles) {
        fileCopies[0] = iiLookupFileFromFP(inFP);
        if (!fileCopies[0])
          printf("WARNING: %s - Could not find iiFile from file pointer to assess " \
                 "whether to read a TIFF file in parallel\n", progname);
        if (fileCopies[0] && fileCopies[0]->file == IIFILE_TIFF) 
          maxReadThreads = tiffNumReadThreads(nx, ny, fileCopies[0]->tiffCompression,
                                              MAX_READ_THREADS);
      }
    }

    // Get frames to use from file and make sure it is legal
    numFrameUse = inHead.nz;
    if (startFrame > 0)
      numFrameUse = B3DMIN(inHead.nz, endFrame) + 1 - startFrame;
    if (numFrameUse < 1)
      exitError("No frames would be included for file %s which has only %d frames",
                inFiles[ind], inHead.nz);
    if (!combineFiles && numFrameUse < breakSetSize)
      exitError("The available frames for file %s is %d, less than the set size of %d",
                inFiles[ind], numFrameUse, breakSetSize);

    // Check for gain reference if one not entered
    hasExtra = inHead.next && !extraIsNbytesAndFlags(inHead.nint, inHead.nreal);
    iz = 0;
    if (hasExtra && !gainName && 
        inHead.next >= inHead.nz * 4 * (inHead.nint + inHead.nreal) + 4 * nx * ny)
      iz = 1;
    if (!ind)
      extraHasGainRef = iz;
    else if (extraHasGainRef != iz)
      exitError("All files must have gain references in their extended header if any do;"
                " it is missing in %s", inFiles[ind]);

    // Get the min and max set sizes if breaking, or if extra header has tilt angles,
    // or for a frame list file
    minSet = 0;
    if (breakSetSize > 0 && !combineFiles) {
      minMaxSetSize(breakSetSize, numFrameUse, minSet, maxSet);
      numSets = numFrameUse / breakSetSize;
    } else if (savedListName) {
      numSets = setStarts.size() - 1;
      minSet = inHead.nz;
      maxSet = 0;
      for (ind = 0; ind < numSets; ind++) {
        ACCUM_MIN(minSet, setStarts[ind + 1] - setStarts[ind]);
        ACCUM_MAX(maxSet, setStarts[ind + 1] - setStarts[ind]);
      }

    } else if (!combineFiles && hasExtra) {
      ierr = analyzeExtraHeader(inFP, &inHead, startFrame, endFrame, breakSetSize > 0, 
                                &extraBuf, extraBufSize, extraTilts, setStarts, minSet,
                                maxSet, fileAxis, filePix);

      // Save the axis rotation and pixel size if any, and make sure all files are 
      // consistent
      fileHasTilts = (namesFromMdoc || breakSetSize > 0) ? 0 : ierr / 2;
      if (!ind) {
        extraHasAxisPix = ierr % 2;
        if (ierr % 2) {
          extraAxis = fileAxis;
          extraPixSize = filePix;
        }
        extraHasTilts = fileHasTilts;
      }

      if (extraHasAxisPix != ierr % 2 || 
          (extraHasAxisPix && (fabs(extraAxis - fileAxis) > 0.01 || 
                               fabs(extraPixSize - filePix) > 0.01)))
        exitError("All files must have the same axis rotation angles and pixel sizes in "
                  "the extended header if any do; they differ in %s", inFiles[ind]);
      if (extraHasTilts != fileHasTilts)
        exitError("All files must have valid tilt angles in extended header if any do "
                  "and if -break is not entered; they are invalid in %s", inFiles[ind]);
      if (extraHasTilts) {
        numSets = extraTilts.size();
        if (!tiltName && !stackName) {
          tiltAngles.reserve(tiltAngles.size() + extraTilts.size());
          tiltAngles.insert(tiltAngles.end(), extraTilts.begin(), extraTilts.end());
        }
      }
    }

    // Keep track of minimum and maximum set size if sets came out either way
    if (minSet) {
      if (!minSetSize) {
        minSetSize = minSet;
        maxSetSize = maxSet;
      } else {
        ACCUM_MIN(minSetSize, minSet);
        ACCUM_MAX(maxSetSize, maxSet);
      }
      numAllSets += numSets;
    }
    if (!ind || !skipChecks)
      iiFClose(inFP);

    ACCUM_MAX(maxNumZ, numFrameUse);
    ACCUM_MAX(maxFrameDoses, inHead.nz);
  }

  // Finish up with processing a list of input files: just fix the # of files
  if (fileListFP) {
    fclose(fileListFP);
    numInFiles = (int)inFiles.size();
    if (!numInFiles)
      exitError("There were no input files in the list file %s", listName);
    free(listName);
    printf("%d files in input file list\n", numInFiles);
  }

  // Handle combination of single-frame files or breaking frames into sets
  if (combineFiles) {
    if (numInFiles < combineFiles)
      exitError("The break entry, %d, is bigger than the number of single-frame input "
                "files, %d", combineFiles, numInFiles);
    numSingleFiles = numInFiles;
    minMaxSetSize(combineFiles, numInFiles, ierr, maxNumZ);
    maxFrameDoses = maxNumZ;
    numInFiles /= combineFiles;
    printf("%d files will be combined into %d summed images\n", numSingleFiles,
           numInFiles);
  } else if (breakSetSize > 0 || extraHasTilts || savedListName) {
    numSingleFiles = numInFiles;
    numInFiles = numAllSets;
    maxNumZ = maxSetSize;
    maxFrameDoses = maxSetSize;
    if (extraHasTilts)
      printf("Tilt angles from extended header will be used to break frames into sets\n");
    printf("Frames from %d files will be broken into %d summed images\n", numSingleFiles,
           numInFiles);
  }
  if (extraHasGainRef)
    printf("Gain reference from extended header will be applied to frames\n");
  if ((gainName || extraHasGainRef) && inHead.mode == MRC_MODE_BYTE && !enteredScale &&
      !totalScale && outMode != MRC_MODE_FLOAT) {
    printf("Applying default total scaling of %g because byte values are being "
           "gain-normalized\n", defaultByteScale);
    totalScale = defaultByteScale;
  }

  // Now that number of "files" is known, and maximum number of sections, return to
  // dealing with dose-weighting
  if (totalDose > 0. || doseFileType > 0 || !fixedFrameDoses.empty()) {
    frameDoses.resize(maxFrameDoses);
    tempVal1.resize(2 * maxFrameDoses);
  }
  if (doseFileType > 0) {
    ind = doseFileType < 4 ? totalDoseVec.size() : frameDoseLines.size();
    if (doseFileType != 4 && ind < numInFiles)
      exitError("The dose file has fewer lines (%d) than frame sets to be aligned %d",
                ind, numInFiles);
    if ((combineFiles > 0 || breakSetSize > 0) && doseFileType == 4)
      exitError("You cannot use an mdoc for a dose file when combining files or "
                "breaking frames into sets");
    if (doseFileType == 4 && !namesFromMdoc) {
      
      // Need to match frame paths in mdoc to actual files being aligned
      fullFramePaths.resize(numMdocSect, NULL);
      tempVal1.resize(numMdocSect, 0);
      if (getMetadataByKey(adocInd, adocType, numMdocSect, "SubFramePath", 0,
                           &tempVal1[0], NULL, NULL, &fullFramePaths[0], &ind, &iz,
                           numMdocSect, &izPiece[0]) || !iz)
        exitError("Getting all frame paths from mdoc file");

        // Got some names: reduce them all to filename only
        mdocFileTails.resize(numMdocSect);
        for (ind = 0; ind < numMdocSect; ind++) {
          if (fullFramePaths[ind]) {
            extractFileTail(fullFramePaths[ind], sstr);
            mdocFileTails[ind] = sstr;
            B3DFREE(fullFramePaths[ind]);
          }
        }

        // Loop on the filenames
        for (ifile = 0; ifile < numInFiles; ifile++) {
          extractFileTail(inFiles[ifile], sstr);
          for (ind = 0; ind < numMdocSect; ind++) {
            if (!mdocFileTails[ind].compare(sstr)) {
              totalDoseVec.push_back(doseFromMdoc[ind]);
              priorDoseVec.push_back
                (B3DCHOICE((mdocName && doseAccumulates > 0) || doseAccumulates > 1,
                           priorFromMdoc[ind], 0.));
              mdocFileTails[ind] = "";
              break;
            }
          }
          if (ind >= numMdocSect)
            exitError("No section was found in the mdoc file with a filename matching "
                      "input file %s", inFiles[ifile]);
          ierr = AdocGetString(ADOC_ZVALUE_NAME, ind, FRAME_DOSE_KEY, &tempname);
          if (ierr < 0)
            exitError("Trying to access "FRAME_DOSE_KEY" from section %d in mdoc file",
                      ind);
          if (!ierr) {
            frameDoseLines[ifile] = tempname;
            free(tempname);
          } 
        }
      }

    // Now want to check frameDoses for reasonableness and get a total dose for type 5
    if (doseFileType > 3) {
      for (ifile = 0; ifile < numInFiles; ifile++) {
        expandFrameDosesNumbers(frameDoseLines[ifile], &tempVal1[0], maxFrameDoses, ind, 
                               doseTemp, NULL);
        if (doseFileType > 4) {
          totalDoseVec.push_back(doseTemp);
          priorDoseVec.push_back(0.);
        }
      }
    }
  }

  // Single frame dose entry: check it, fill arrays, and pretend it is type 5
  if (!fixedFrameDoses.empty()) {
    expandFrameDosesNumbers(fixedFrameDoses, &tempVal1[0], maxFrameDoses,
                            ind, doseTemp, NULL);
    for (ifile = 0; ifile < numInFiles; ifile++) {
      frameDoseLines.push_back(fixedFrameDoses);
      totalDoseVec.push_back(doseTemp);
      priorDoseVec.push_back(0.);
    }
    doseFileType = 5;
  }

  // Fixed dose, fill the arrays
  if (totalDose > 0.) {
    totalDoseVec.resize(numInFiles, totalDose);
    priorDoseVec.resize(numInFiles, 0.);
  }
    
  // If assuming a tilt series, assign prior doses when none available
  if ((doseAccumulates > 0 && (doseFileType == 1 || doseFileType > 4 || totalDose > 0.))
    || (doseAccumulates == 1 && doseFileType == 4 && !mdocName))
    priorDosesFromImageDoses(&totalDoseVec[0], (int)totalDoseVec.size(), numBidir,
                             &priorDoseVec[0]);

  // Adjust default memory if physical memory is available, and get 
  physMem = b3dPhysicalMemory() / (1024. * 1024. * 1024.);
  if (physMem > 0.) {
    if (physMem < 16.)
      memoryLimit = 0.75 * physMem;
    if (physMem > 24.)
      memoryLimit = 0.5 * physMem;
  }
  ind = 0;
  if (PipGetFloatArray("MemoryLimitGB", memLimits, &ind, 2) == 0) {
    memoryLimit = memLimits[0];
    if (memLimits[0] < -0.95 || (ind > 1 && memLimits[1] < -0.95) ||
        fabs(memLimits[0]) < 0.05 || (ind > 1 && fabs(memLimits[1]) < 0.05))
      exitError("You cannot enter a memory limit below -0.95 or between -0.05 and 0.05");
    if (memLimits[0] < 0) {
      if (!physMem)
        exitError("You cannot enter a negative CPU memory limit: system memory "
                  "not available");
      memoryLimit *= -physMem;
    }
    if (ind > 1)
      gpuMemLimit = memLimits[1];
  }

  // Get lots more options
  if (totalScale > 0.)
    scale = totalScale / alreadyScaledBy;
  xShifts = B3DMALLOC(float, maxNumZ);
  yShifts = B3DMALLOC(float, maxNumZ);
  bestXshifts = B3DMALLOC(float, maxNumZ);
  bestYshifts = B3DMALLOC(float, maxNumZ);
  rawXshifts = B3DMALLOC(float, maxNumZ);
  rawYshifts = B3DMALLOC(float, maxNumZ);
  bestXraw = B3DMALLOC(float, maxNumZ);
  bestYraw = B3DMALLOC(float, maxNumZ);
  if (!xShifts || !yShifts || !bestXshifts || !bestYshifts || !bestXraw || !bestYraw ||
      !rawXshifts || !rawYshifts)
    exitError("Allocating arrays for shifts");
  PipGetInteger("PairwiseFrames", &numAVAinput);
  PipGetFloat("TaperFraction", &taperFrac);
  if (!taperFrac)
    fullTaperFrac = 0.05;
  PipGetFloat("TrimFraction", &trimFrac);
  PipGetBoolean("ReverseOrder",  &reverse);
  PipGetInteger("ShiftLimit", &shiftLimit);
  PipGetFloat("TruncateAbove", &truncLimit);
  PipGetString("TransformExtension", &xfExt);
  PipGetInteger("DebugOutput", &ifDebug);
#ifdef TEST_SHRMEM
  PipGetInteger("ShrMem", &useShrMem);
#endif
  PipGetFloat("KFactorForFits", &kFactor);
  PipGetFloat("MaxResidualWeight", &maxMaxWeight);
  PipGetFloat("GoodEnoughError", &goodEnough);
  PipGetBoolean("UseHybridShifts", &hybridShifts);
  PipGetFloat("FilterSigma1", &sigma1);
  PipGetFloat("FilterSigma2", &sigma2);
  PipGetFloat("FilterRadius1", &radius1);
  PipGetFloat("FilterRadius2", &radius2);
  PipGetInteger("RefineAlignment", &refineAtEnd);
  PipGetFloat("RefineRadius2", &refRadius2);
  PipGetInteger("AntialiasFilter", &antiFiltType);
  PipGetFloat("RingSpacingForFRC", &frcDeltaR);
  PipGetInteger("UseGPU", &useGPU);
  PipGetInteger("GroupSize", &groupSize);
  PipGetInteger("RefineWithGroupSums", &groupRefine);
  PipGetFloat("StopIterationsAtShift", &iterCrit);
  PipGetFloat("PixelSize", &optionPixSize);
  splineSmooth = PipGetInteger("MinForSplineSmoothing", &minNumForSpline);
  if (minNumForSpline < 8)
    splineSmooth = 0;
  B3DCLAMP(antiFiltType, 1, 6);
  if (refineAtEnd < 0)
    exitError("Entry for -refine cannot be negative");
  numAllVsAll = numAVAinput;
  if (numAllVsAll < 0) {
    if (numAVAinput < -4)
      exitError("The value for the -pair option cannot be more negative than -4");
    if (numAVAinput == -1)
      numAllVsAll = B3DMIN(MAX_ALL_VS_ALL, maxNumZ + 4);
    else
      numAllVsAll = B3DMAX(minFractionalAVA, 
                           (maxNumZ - 1 - numAVAinput) / (-numAVAinput));
  }
  if (startAssess > 0 && !numAllVsAll)
    exitError("You cannot set frames for assessing fits with cumulative correlations");
  if (reverse && (breakSetSize > 0 || extraHasTilts))
    exitError("You cannot process in reverse when breaking frames into sets");

  // Get default binning for size
  minDiff = 1.e20;
  for (ind = 0; ind < (int)(sizeof(defaultBinnings) / sizeof(int)); ind++) {
    diff = fabs(targetAliSize - sqrt((double)nx * ny) / defaultBinnings[ind]);
    if (diff < minDiff) {
      minDiff = diff;
      alignBin = defaultBinnings[ind];
    }
  }

  // Set output size based on this binning of frames
  alignBinIn = alignBin;
  alsumBinEntered = 1 - PipGetTwoIntegers("AlignAndSumBinning", &alignBinIn, &sumBin);
  if (alignBinIn == 0 || sumBin < 1 || alignBinIn > 16 || sumBin > 16)
    exitError("Binning value is out of allowed range");
  if (alignBinIn > 0)
    alignBin = alignBinIn;
  if (sumRotationFlip < 0 || sumRotationFlip > 7)
    exitError("Inappropriate value of rotation and flip for sum entered");
  nxOut = nxSum = B3DCHOICE(sumRotationFlip % 2, ny, nx) / sumBin;
  nyOut = nySum = B3DCHOICE(sumRotationFlip % 2, nx, ny) / sumBin;

  // Get tilt angles from file
  if (tiltName) {
    inFP = fopen(tiltName, "r");
    if (!inFP)
      exitError("Opening tilt angle file %s", tiltName);
    while (1) {
      ierr = fgetline(inFP, title, MRC_LABEL_SIZE);
      if (ierr == -2)
        break;
      if (ierr == -1)
        exitError("Reading tilt angle file %s", tiltName);
      tiltAngles.push_back((float)atof(title));
      if (ierr < 0)
        break;
    }
    if ((int)tiltAngles.size() < numInFiles)
      exitError("There are fewer tilt angles in the file (%d) than frame files or sets"
                " (%d)", tiltAngles.size(), numInFiles);
    if ((int)tiltAngles.size() > numInFiles)
      printf("WARNING: %s - There are fewer frame files or sets (%d) than tilt angles "
             "in the file (%d)\n", progname, numInFiles, (int)tiltAngles.size());
  }

  // Collect information from stack or mdoc file
  if (stackName) {
    stackFP = openAndReadHeader(stackName, &stackHead, "stack");
    stackMode = stackHead.mode;
    nxStack = stackHead.nx;
    nyStack = stackHead.ny;
    mainHead = stackHead;
    unwgtHead = stackHead;
    numSect = stackHead.nz;
  } 

  if (mdocName) {
    if (stackName)
      exitError("You cannot enter both a corresponding stack and an mdoc file");
    
    // Get tilt angles if not already got them
    if (!tiltName) {
      tiltAngles.resize(numSect);
      for (ind = 0; ind < numSect; ind++) {
        if (ignoreZvalue) 
          iz = ind;
        else
          iz = AdocLookupByNameValue(ADOC_ZVALUE_NAME, ind);
        if (iz < 0)
          exitError("Looking up section with Z value %d in mdoc file", ind);
        if (AdocGetFloat(ADOC_ZVALUE_NAME, iz, "TiltAngle", &tiltAngles[ind]))
          exitError("Getting tilt angle for %s %d in mdoc file", ind,
                    ignoreZvalue ? "section" : "Z value");
      }
    }
    
    // Get titles
    ind = AdocGetNumberOfSections("T");
    if (ind < 0)
      exitError("Looking up titles in mdoc file");
    for (outNum = 0; outNum < numOutFiles; outNum++) {
      for (iz = 0; iz < ind; iz++) {
        if (AdocGetSectionName("T", iz, &filename))
          exitError("Getting title from mdoc file");
        sstr = filename;
        if ((sstr.find("Tilt axis angle") != string::npos) && (sstr[0] != ' '))
          sstr = "    " + sstr;
        strncpy(&outHeads[outNum]->labels[iz][0], sstr.c_str(), MRC_LABEL_SIZE);
        fixTitlePadding(&outHeads[outNum]->labels[iz][0]);
        free(filename);
      }
      outHeads[outNum]->nlabl = iz;
    }
  }

  // Make sure things work out
  if (stackName || mdocName) {
    if (stackMode != MRC_MODE_BYTE && !enteredMode)
      outMode = stackMode;
    if (numSect < numInFiles && !tiltName)
      exitError("There are fewer sections in the %s (%d) than frame files or sets (%d)", 
                stackName ? "stack" : "mdoc file", numSect, numInFiles);
    if (numSect > numInFiles && !tiltName)
      printf("WARNING: %s - There are fewer frame sets or files (%d) than sections in "
             "the %s (%d)\n", progname, numInFiles, stackName ? "stack" : "mdoc file",
             numSect);
    
    // Figure out if size works, requires trimming, or implies a binning relative to stack
    if ((nxStack <= nxOut && nyStack > nyOut) || (nxStack > nxOut && nyStack <= nyOut))
      exitError("The image size from the %s is bigger in one dimension than the frame "
                "size", stackName ? "stack" : "mdoc file");
    ierr = 0;

    // Look for integer binning difference that matches in each direction
    // And make sure each size is close enough after scaling
    if (nxStack <= nxOut) {
      relXbin = B3DNINT((double)nxOut / nxStack);
      relYbin = B3DNINT((double)nyOut / nyStack);
      if (fabs((double)relXbin * nxStack - nxOut) > sizeDiffCrit * nxOut || 
          fabs((double)relYbin * nyStack - nyOut) > sizeDiffCrit * nyOut || 
          relXbin != relYbin)
        ierr = 1;
      relBinning = 1. / relXbin;
    } else {
      relXbin = B3DNINT((double)nxStack / nxOut);
      relYbin = B3DNINT((double)nyStack / nyOut);
      if (fabs((double)relXbin * nxOut - nxStack) > sizeDiffCrit * nxStack || 
          fabs((double)relYbin * nyOut - nyStack) > sizeDiffCrit * nyStack || 
          relXbin != relYbin)
        ierr = 1;
      relBinning = relXbin;
    }
    if (ierr) 
      exitError("The image size does not correspond well enough between the %s and the "
                "frames to deduce their relationship", stackName ? "stack" : "mdoc file");

    // For same binning, trim if there is a small difference
    if (relXbin == 1) {
      if (nxOut - nxStack > trimCrit || nyOut - nyStack > trimCrit) {
        printf("The image size from the %s is significantly smaller and frames will not "
               "be trimmed to that size\n", stackName ? "stack" : "mdoc file");
      } else if (nxStack < nxOut || nyStack < nyOut) {
        printf("The image size from the %s is slightly smaller and frames will be trimmed"
               " to that size\n", stackName ? "stack" : "mdoc file");
        nxOut = nxStack;
        nyOut = nyStack;
      }
    } else {

      // Different binnings: look at labels to try to adjust it there
      printf("The %s is at a different binning from the frames; frame sizes will not be"
             " adjusted\n", stackName ? "stack" : "mdoc file"); 
      for (outNum = 0; outNum < numOutFiles; outNum++) {
        for (ind = 0; ind < outHeads[outNum]->nlabl; ind++) {
          strncpy(title, (const char *)&outHeads[outNum]->labels[ind], MRC_LABEL_SIZE);
          title[MRC_LABEL_SIZE] = 0x00;
          if (adjustTitleBinning(title, sstr, relBinning, true)) {
            strncpy((char *)&outHeads[outNum]->labels[ind][0], &sstr[0], MRC_LABEL_SIZE);
            fixTitlePadding(&outHeads[outNum]->labels[ind][0]);
            break;
          }
        }
      }

      // Find and fix title in mdoc too
      if (adjustMdoc) {
        nz = AdocGetNumberOfSections("T");
        for (ind = 0; ind < nz; ind++) {
          if (AdocGetSectionName("T", ind, &filename) >= 0) {
            if (adjustTitleBinning(filename, sstr, relBinning, false)) {
              if (AdocChangeSectionName("T", ind, sstr.c_str()))
                exitError("Adjusting title with binning in mdoc file");
              free(filename);
              break;
            }
            free(filename);
          }
        }
      }
    }      
  }

  groupSize = B3DMAX(1, groupSize);
  if (groupSize > 1) {
    ierr = B3DMIN(maxNumZ, numAllVsAll + groupSize - 1) + 1 - groupSize;
    useBlockGroup = ((ierr + 1 - groupSize) * (ierr - groupSize)) / 2 < ierr;
    if (useBlockGroup)
      printf("Using block grouping; too few frames being fit for slide grouping\n");
    else if (numAllVsAll <= MAX_ALL_VS_ALL + 1 - groupSize)
      numAllVsAll += groupSize - 1;
  } else
    groupRefine = 0;

  // Set up the binnings to test
  numBinTests = 1;
  numFiltTests[0] = 1;
  numVaries = 0;
  varyRadius2[0][0] = radius2;
  varySigma2[0][0] = sigma2;
  ierr = 0;
  if (PipGetIntegerArray("TestBinnings", binsToTest, &ierr, MAX_BINNINGS) == 0) {
    numBinTests = ierr;
    for (ind = 0; ind < numBinTests; ind++) {
      if (binsToTest[ind] < 1 || binsToTest[ind] > 16)
        exitError("Binning value %d not allowed", binsToTest[ind]);
      ACCUM_MIN(minBinningToTest, binsToTest[ind]);
    }
  } else {
    if (!alsumBinEntered || alignBinIn < 0)
      printf("Selected the default binning of %d for this image size\n", alignBin);
    binsToTest[0] = alignBin;
    minBinningToTest = alignBin;
  }

  // Get the filters to test
  PipNumberOfEntries("VaryFilter", &numVaries);
  
  for (ind = 0; ind < numBinTests; ind++) {
    if (ind < numVaries) {
      numFiltTests[ind] = 0;
      PipGetFloatArray("VaryFilter", &varyRadius2[ind][0], &numFiltTests[ind],
                       MAX_FILTERS);
      rsSortFloats(&varyRadius2[ind][0], numFiltTests[ind]);
      useInd = ind;
    } else {
      useInd = B3DMAX(0, numVaries - 1);
      numFiltTests[ind] = numFiltTests[useInd];
    }
    for (iz = 0; iz < numFiltTests[ind]; iz++) {
      varyRadius2[ind][iz] = varyRadius2[useInd][iz];
      varySigma2[ind][iz] = B3DNINT(SIG2_ROUND_FAC * sigma2 * varyRadius2[ind][iz] /
                                    radius2) / SIG2_ROUND_FAC;
      numTimesBest[ind][iz] = 0;
    }
  }

  if (!refRadius2)
    refRadius2 = varyRadius2[0][0];
  refSigma2 = B3DNINT(SIG2_ROUND_FAC * sigma2 * refRadius2 / radius2) / SIG2_ROUND_FAC;

  // Gain reference
  gainHead.nx = gainHead.ny = 0;
  if (gainName) {
    extraFP = openAndReadHeader(gainName, &gainHead,"gain reference");
    free(gainName);
    if (gainHead.mode != MRC_MODE_FLOAT)
      exitError("Gain reference must be floating point");
    gainSlice = sliceReadMRC(&gainHead, 0, 'Z');
    if (!gainSlice)
      exitError("Reading gain reference file");
    iiFClose(extraFP);
    nxGain = gainHead.nx;
    nyGain = gainHead.ny;
    if (rotationFlip)
      rotateFlipGainReference(gainSlice->data.f, nxGain, nyGain, rotationFlip);
    if (nxGain < nx || nyGain < ny)
      exitError("Gain reference is smaller than image in %s%s%s", nxGain < nx ? "X" : "",
                nxGain < nx && nyGain < ny ? " and " : "", nyGain < ny ? "Y" : "");
  } else if (extraHasGainRef) {
    gainSlice = sliceCreate(nx, ny, SLICE_MODE_FLOAT);
    if (!gainSlice)
      exitError("Allocating memory for gain reference");
  }
      
  // Dark reference
  if (PipGetString("DarkReferenceFile", &extraName) == 0) {
    extraFP = openAndReadHeader(extraName, &darkHead, "dark reference");
    free(extraName);
    if (darkHead.mode != MRC_MODE_SHORT && darkHead.mode != MRC_MODE_USHORT)
      exitError("Dark reference must be signed or unsigned short integers");
    if (darkHead.nx != nx || darkHead.ny != ny)
      exitError("Dark reference is not the same size as the image");
    darkSlice = sliceReadMRC(&darkHead, 0, 'Z');
    if (!darkSlice)
      exitError("Reading dark reference file");
    iiFClose(extraFP);
  }
      
  // Defect file
  if (PipGetString("CameraDefectFile", &extraName) == 0) {
#ifdef TEST_SHRMEM
    if (useShrMem)
      defectFileToString(extraName, defectString);
#endif
    
    ierr = CorDefParseDefects(extraName, 0, defects, camSizeX, camSizeY);
    if (ierr)
      exitError("%s defect file %s\n", (ierr == 1) ? "Opening" : 
                  "Reading or parsing lines in", extraName);
    free(extraName);
    if (!camSizeX || !camSizeY)
      exitError("Defect list file must have CameraSizeX and CameraSizeY entries");
    PipGetFloat("ImagesAreBinned", &imagesBinned);
    PipGetInteger("DoubleDefectCoords", &scaleDefects);
    CorDefFlipDefectsInY(&defects, camSizeX, camSizeY, 0);
    CorDefFindTouchingPixels(defects, camSizeX, camSizeY, 0);
    if (CorDefSetupToCorrect(inHead.nx, inHead.ny, defects, camSizeX, camSizeY, 
                             scaleDefects, imagesBinned, corDefBinning, "-imagebinned"))
      exitError("Image size is more than twice the size stored in the camera defect "
                "list");
  }

  // Get sizes
  sFA.getPadSizesBytes(nx, ny, fullTaperFrac, sumBin, minBinningToTest, fullPadSize, 
                       sumPadSize, alignPadSize);

  // See about GPU
  gpuFlags = 0;
  needPreprocess = truncLimit > 0. || camSizeX > 0 || gainSlice != NULL;
  fullDataSize = sizeof(float);
  if (useGPU && taperFrac <= 0. && trimFrac <= 0.) {
    printf("The GPU cannot be used when the taper fraction is set to 0\n");
      useGPU = -1;
  }

  if (useGPU >= 0) {
    if (useShrMem) {
#ifdef TEST_SHRMEM
      if (!sSMC.gpuAvailable(useGPU, &gpuMemory, ifDebug % 10))
        useGPU = -1;
#endif
    } else {
    if (!sFA.gpuAvailable(useGPU, &gpuMemory, ifDebug % 10))
      useGPU = -1;
    }
  }
  if (useGPU >= 0) {
    if (sFA.gpuAvailable(useGPU, &gpuMemory, ifDebug % 10)) {
      gpuUsableMem = gpuMemory * gpuFracMem;
      if (gpuMemLimit > 0)
        gpuUsableMem = 1024. * 1024. * 1024. * gpuMemLimit;
      else if (gpuMemLimit < 0)
        gpuUsableMem = -gpuMemory * gpuMemLimit;
      needed = 0.;
      nzAlign = maxNumZ;
      if (useBlockGroup)
        nzAlign = B3DMAX(1, nzAlign / groupSize);
      sFA.gpuMemoryNeeds(fullPadSize, sumPadSize, alignPadSize, numAllVsAll, nzAlign,
                         refineAtEnd, groupSize, needForGPUsum, needForAli);
      if (numOutFiles > 1) 
        needForGPUsum += sumPadSize;

      // First see if summing can be done
      if (!testMode) {
        if (needForGPUsum > gpuUsableMem) {
          printf("Insufficient memory on GPU to use it for summing (%.0f MB needed of "
                 "%.0f MB total)\n", needForGPUsum / 1.048e6, gpuMemory / 1.048e6);
          useGPU = -1;
        } else {
          gpuFlags = GPU_FOR_SUMMING + (numOutFiles > 1 ? GPU_DO_UNWGT_SUM : 0);
          needed = needForGPUsum;
        }
      }

      // Summing will be done with alignment if only one binning, only one filter or
      // using hybrid shift, and not assessing or doing spline or refining at end
      doSpline = (splineSmooth && maxNumZ >= minNumForSpline) ? 1 : 0;
      sumWithAlign = ((numBinTests == 1 && (hybridShifts || numFiltTests[0] == 1)) ||
                      startAssess >= 0) && !testMode && !doSpline && !refineAtEnd;

      // If that is the case, and alignment alone would fit but both would not,
      // then see if deferring the sum will require less memory than the limit and if
      // so, then defer the summing
      // Call this first unconditionally to get numHoldFull set properly
      sFA.totalMemoryNeeds(fullPadSize, 4, sumPadSize, alignPadSize, numAllVsAll, nzAlign,
                           refineAtEnd, numBinTests, numFiltTests[0], hybridShifts, 
                           groupSize, doSpline, GPU_FOR_ALIGNING, 0, testMode,
                           startAssess, sumInOnePass, numHoldFull);
      if (sumWithAlign && needForAli < gpuUsableMem && needForAli + needed > 
          gpuUsableMem) {
        tmean = sFA.totalMemoryNeeds
          (fullPadSize, 4, sumPadSize, alignPadSize, numAllVsAll, nzAlign, refineAtEnd,
           numBinTests, numFiltTests[0], hybridShifts, groupSize, doSpline,
           GPU_FOR_ALIGNING, 1, testMode, startAssess, sumInOnePass, numHoldFull);
        if (tmean < memoryLimit) {
          sumWithAlign = false;
          deferSum = 1;
        }
      }

      // If summing is done with aligning add the usage for summing if any
      if (sumWithAlign)
        needForAli += needed;
      
      // Decide if alignment can be done there
      if (needForAli > gpuUsableMem) {
        printf("Insufficient memory on GPU to do alignment (%.0f MB needed"
               " of %.0f MB total)\n", needForAli / 1.048e6, gpuMemory / 1.048e6);
      } else {
        if (sumWithAlign)
          needed = needForAli;
        gpuFlags |= GPU_FOR_ALIGNING;
      }

      // Now if summing, see if there is room for even/odd
      if (gpuFlags & GPU_FOR_SUMMING) {
        needed += sumPadSize;
        if (needed > gpuUsableMem) {
          printf("Insufficient memory on GPU to get even/odd sums for FRC (%.0f MB needed"
                 " of %.0f MB total)\n", needed / 1.048e6, gpuMemory / 1.048e6);
          gettingFRC = false;
          needed -= sumPadSize;
        } else {
          gpuFlags |= GPU_DO_EVEN_ODD;
        }
      }

      // Now see about adding bin, noise, and preproc to GPU
      if (gpuFlags) {

        // If option is entered, just use it to set the flags
        ind = 0;
        if (PipGetInteger("FlagsForGPU", &ind) == 0) {
          if (ind % 10)
            gpuFlags |= GPU_DO_NOISE_TAPER;
          if ((ind / 10) % 10)
            gpuFlags |= GPU_DO_BIN_PAD;
          if (gpuFlags & (GPU_DO_NOISE_TAPER | GPU_DO_BIN_PAD)) {
            if ((ind / 100) % 10) {
              gpuFlags |= STACK_FULL_ON_GPU;
              if (ind / 10000)
                gpuFlags |= GPU_STACK_LIMITED + ((ind / 10000) << GPU_STACK_LIM_SHIFT);
            }
            if ((ind / 1000 % 10) && needPreprocess) {
              gpuFlags |= GPU_DO_PREPROCESS;
              if (gainSlice)
                gpuFlags |= GPU_DO_GAIN_NORM;
              if (camSizeX > 0)
                gpuFlags |= GPU_CORRECT_DEFECTS;
            }
          }
        } else {

          // Otherwise, need to analyze which preprocessing/bin/pad operations to perform
          needForPreOps = sFA.findPreprocPadGpuFlags
            (nx, ny, darkSlice ? sizeof(float) : maxDataSize, minBinningToTest,
             !darkSlice && gainSlice != NULL, !darkSlice && camSizeX > 0,
             !darkSlice && truncLimit > 0., B3DMAX(numHoldFull, 1), 
             gpuUsableMem - needed, 0.5 * (1. - gpuFracMem) * gpuUsableMem, gpuFlags, 
             gpuFlags);
          needed += needForPreOps;
        }
        ind = (gpuFlags >> GPU_STACK_LIM_SHIFT) & GPU_STACK_LIM_MASK;
        if (ifDebug && (gpuFlags & (GPU_DO_NOISE_TAPER | GPU_DO_BIN_PAD)))
          printf("%s  %s  %s  %s %s %d on GPU\n", (gpuFlags & GPU_DO_NOISE_TAPER) ? 
                 "noise-pad" : "", (gpuFlags & GPU_DO_BIN_PAD) ? "bin-pad" : "",
                 (gpuFlags & GPU_DO_PREPROCESS) ? "preprocess" : "", 
                 (gpuFlags & STACK_FULL_ON_GPU) ? "stack" : "", 
                 (gpuFlags & GPU_STACK_LIMITED) ? "limit" : "   ",
                 (gpuFlags & GPU_STACK_LIMITED) ? ind : 0);

        // Can stack smaller size if doing either operation on GPU and either there is
        // no preprocess or preprocessing is on GPU
        if ((gpuFlags & (GPU_DO_BIN_PAD | GPU_DO_NOISE_TAPER)) && 
            (!needPreprocess || (gpuFlags & GPU_DO_PREPROCESS)))
          fullDataSize = maxDataSize;
      }
    }
  }


  // FRC output file
  if (PipGetString("FRCOutputFile", &extraName) == 0) {
    if (testMode)
      exitError("There is no FRC output available when not making sums");
    imodBackupFile(extraName);
    frcFP = fopen(extraName, "w");
    if (!frcFP)
      exitError("Opening file for FRC curves, %s", extraName);
    free(extraName);
  }

  // Shift output file
  if (PipGetString("PlottableShiftFile", &extraName) == 0) {
    imodBackupFile(extraName);
    plotFP = fopen(extraName, "w");
    if (!plotFP)
      exitError("Opening file for shift curves, %s", extraName);
    free(extraName);
  }
  PipDone();

  if (!testMode) {

    // Set up output header(s)
    for (outNum = 0; outNum < numOutFiles; outNum++) {
      headPtr = outHeads[outNum];
      headPtr->nz = numInFiles;
      headPtr->mz = numInFiles;
      headPtr->nx = nxOut;
      headPtr->ny = nyOut;
      headPtr->mx = headPtr->nx;
      headPtr->my = headPtr->ny;
      headPtr->mode = outMode;
      headPtr->amax = -1.e30;
      headPtr->amin = 1.e30;
      headPtr->amean = 0.;
      if (xScale == 1.0 && extraHasAxisPix)
        xScale = yScale = zScale = extraPixSize;
      else if (xScale == 1.0 && stackName) {
        stackBinX = nxOut * sumBin / stackHead.nx;
        stackBinY = nyOut * sumBin / stackHead.ny;
        stackBin = B3DNINT(stackBinX);
        if (fabs(B3DNINT(stackBinX) - stackBinX) < 0.05 && 
            fabs(B3DNINT(stackBinY) - stackBinY) < 0.05 && stackBin == B3DNINT(stackBinY))
          mrc_get_scale(&stackHead, &xScale, &yScale, &zScale);
      }
      if (optionPixSize > 0.)
        xScale = yScale = zScale = 10. * optionPixSize;
      mrc_set_scale(headPtr, sumBin * xScale, sumBin * yScale, sumBin * zScale);
      if (extraHasAxisPix) {
        sprintf(title, "    Tilt axis angle = %.2f", extraAxis);
        ind = B3DMIN(headPtr->nlabl, 8);
        strncpy((char *)&headPtr->labels[ind][0], title, MRC_LABEL_SIZE);
        fixTitlePadding(&headPtr->labels[ind][0]);
        headPtr->nlabl = B3DMAX(headPtr->nlabl, ind + 1);
      }
      if (sumBin > 1)
        sprintf(title, "alignframes: summed frames scaled by %g, reduced by %d", scale, 
                sumBin);
      else
        sprintf(title, "alignframes: summed frames scaled by %g", scale);
      mrc_head_label(headPtr, title);
      mrcInitOutputHeader(headPtr);
    
      // Set up output file(s) and li for writing
      imodBackupFile(outNames[outNum]);
      outFPs[outNum] = iiFOpen(outNames[outNum], "wb");
      if (!outFPs[outNum])
        exitError("Opening output file %s", outNames[outNum]);

      // Also create an .openTS file if it seems to be a tilt series
      if (tiltAngles.size() > 0 || stackName || savedListName) {
        openTsNames[outNum] = B3DMALLOC(char, strlen(outNames[outNum]) + 10);
        if (openTsNames[outNum]) {
          sprintf(openTsNames[outNum], "%s.openTS", outNames[outNum]);
          extraFP = fopen(openTsNames[outNum], "w");
          if (extraFP)
            fclose(extraFP);
        }
      }
      
      mrc_init_li(&li, NULL);
      mrc_init_li(&li, headPtr);
      headPtr->fp = outFPs[outNum];
  
      // Need to test for output file type 
      if (b3dOutputFileType() == OUTPUT_TYPE_MRC) {
        if (stackName && stackHead.next && !tiltAngles.size()) {
          headPtr->next = stackHead.next;
          headPtr->nint = stackHead.nint;
          headPtr->nreal = stackHead.nreal;
          headPtr->headerSize = stackHead.headerSize;
          if (mrcCopyExtraHeader(&stackHead, headPtr))
            exitError("Copying extended header from stack to output file");
          iiFClose(stackFP);
        } else if (tiltAngles.size() > 0) {
          numSect = tiltAngles.size();
          headPtr->next = 4 * numSect;
          headPtr->nint = 0;
          headPtr->nreal = 1;
          headPtr->headerSize += 4 * numSect;
          if (fseek(outFPs[outNum], 1024, SEEK_SET) ||
              (int)fwrite(&tiltAngles[0], 4, numSect, outFPs[outNum]) != numSect)
            exitError("Writing tilt angles to extended header of output file");
        }
      }

      // Do modifications to the mdoc
      if (!outNum && adjustMdoc && 
          (AdocSetKeyValue(ADOC_GLOBAL_NAME, 0, "ImageFile", outNames[0]) ||
           AdocSetTwoIntegers(ADOC_GLOBAL_NAME, 0, "ImageSize", nxOut, nyOut) ||
           AdocSetInteger(ADOC_GLOBAL_NAME, 0, "DataMode", outMode) ||
           AdocGetFloat(ADOC_GLOBAL_NAME, 0, "PixelSpacing", &pixTemp) ||
           AdocSetFloat(ADOC_GLOBAL_NAME, 0, "PixelSpacing", pixTemp * relBinning)))
        exitError("Adjusting global values in mdoc");
    }
  }

  // Get data
  summed = B3DMALLOC(float, (nx / sumBin) * (ny / sumBin));
  if (numOutFiles > 1)
    unwgtSum = B3DMALLOC(float, (nx / sumBin) * (ny / sumBin));
  if (!summed || (numOutFiles > 1 && !unwgtSum))
    exitError("Allocating memory for summed slice");
  if (sumRotationFlip) {
    rotSum = B3DMALLOC(float, (nx / sumBin) * (ny / sumBin));
    if (!rotSum)
      exitError("Allocating memory for rotating summed slice");
  }

  // Loop on frame files or frame sets
  superFile = 0;
  setInFile = 0;
  numSets = 0;
  for (ifile = 0; ifile < numInFiles; ifile++) {
    minError = 1.e30;
    if (combineFiles) {
      balancedGroupLimits(numSingleFiles, numInFiles, ifile, &startCombine, 
                          &endCombine);
      filename = inFiles[startCombine];
    } if (breakSetSize > 0 || extraHasTilts || savedListName) {
      filename = inFiles[superFile];
    } else {
      filename = inFiles[ifile];
    }

    // Open file if it is time to do so
    if (!setInFile) {
      inFP = openAndReadHeader(filename, &inHead, "input image");
      checkInputFile(filename, &inHead, nx, ny, combineFiles);
      nz = inHead.nz;
      parallelRead = false;

      // If not combining files, see if the file is a TIFF for parallel reading
      if (!combineFiles && maxReadThreads > 1) {
        fileCopies[0] = iiLookupFileFromFP(inFP);
        if (fileCopies[0] && fileCopies[0]->file == IIFILE_TIFF) {
          numReadThreads = iiOpenCopiesForThreads(fileCopies, maxReadThreads);
          parallelRead = numReadThreads > 1;
        }
      }

      // Get gain reference if it is in there
      if (extraHasGainRef) {
        if (b3dFseek(inFP, MRC_HEADER_SIZE + 4 * (inHead.nint + inHead.nreal) * nz, 
                     SEEK_SET) || 
            (int)b3dFread(gainSlice->data.f, 4, nx * ny, inFP) != nx *ny)
          exitError("Reading gain reference from extended header");
        nxGain = nx;
        nyGain = ny;
        if (rotationFlip > 0)
          rotateFlipGainReference(gainSlice->data.f, nxGain, nyGain, rotationFlip);
      }

      // Set or get group limits
      if (breakSetSize > 0 && !combineFiles) {
        numFrameUse = nz;
        if (startFrame > 0)
          numFrameUse = B3DMIN(nz, endFrame) + 1 - startFrame;
        numSets = numFrameUse / breakSetSize;

        // Get the set limits and then adjust by the start frame
        setStarts.clear();
        for (ind = 0; ind < numSets; ind++) {
          balancedGroupLimits(numFrameUse, numSets, ind, &startCombine, &endCombine);
          setStarts.push_back(startCombine);
        }
        setStarts.push_back(endCombine + 1);
        if (startFrame > 1)
          for (ind = 0; ind <= numSets; ind++)
            setStarts[ind] += startFrame - 1;
        
      } else if (savedListName) {
        numSets = setStarts.size() - 1;

      } else if (extraHasTilts) {
        ierr = analyzeExtraHeader(inFP, &inHead, startFrame, endFrame, false,
                                  &extraBuf, extraBufSize, extraTilts, setStarts, minSet,
                                  maxSet, fileAxis, filePix);
        numSets = extraTilts.size();
      }
    }

    // Now proceed with the current file/set of frames
    if (combineFiles) {
      nz = endCombine + 1 - startCombine;
      fclose(inFP);
    } else if (numSets) {
      startCombine = setStarts[setInFile];
      endCombine = setStarts[setInFile + 1] - 1;
      nz = endCombine + 1 - startCombine;
    }
    extractFileTail(filename, sstr);

    nzAlign = nz;
    if (startFrame > 0 && !numSets)
      nzAlign = B3DMIN(nz, endFrame) + 1 - startFrame;
    numAVAuse = numAVAinput;
    if (numAVAuse == -1)
      numAVAuse = nzAlign;
    else if (numAVAuse < 0)
      numAVAuse = B3DMAX(minFractionalAVA, (nzAlign - 1 - numAVAinput) / (-numAVAinput));

    // Make per-file decision on grouping
    blockGrpSize = slideGrpSize = 1;
    if (groupSize > 1) {
      ierr = nzAlign + 1 - groupSize;
      if (useBlockGroup || 
          ((ierr + 1 - groupSize) * (ierr - groupSize)) / 2 < ierr) {
        blockGrpSize = groupSize;
        if (!useBlockGroup)
          printf("Using block grouping instead of sliding grouping for file # %d\n",
                 ifile + 1);
      } else {
        slideGrpSize = groupSize;
        if (numAVAuse <= MAX_ALL_VS_ALL + 1 - groupSize)
          numAVAuse += groupSize - 1;
      }
    }

    dataSizeForMode(inHead.mode, &inDataSize, &ind);
    needAlloc = inDataSize * nx * ny;
    if (needAlloc > bufAllocSize) {
      readBuf = B3DMALLOC(unsigned char, needAlloc);
      if (blockGrpSize > 1)
        sumBuf = B3DMALLOC(unsigned char, needAlloc);
      if (!readBuf || (blockGrpSize > 1 && !sumBuf))
        exitError("Allocating memory for reading and/or grouping frames");
      bufAllocSize = needAlloc;
    }

    nzAlign = B3DMAX(1, nzAlign / blockGrpSize);
    
    // Set up for spline scaling if criteria met
    doSpline = 0;
    if (splineSmooth > 0 && nzAlign >= minNumForSpline)
      doSpline = 1;
                 
    // Estimate memory usage
    // Does summing in two passes work with cumulative alignment?
    tmean = sFA.totalMemoryNeeds(fullPadSize, fullDataSize, sumPadSize, alignPadSize,
                             numAVAuse, nzAlign, refineAtEnd, numBinTests,
                             numFiltTests[0], hybridShifts, groupSize, doSpline, gpuFlags,
                             deferSum, testMode, startAssess, sumInOnePass, numHoldFull);
    if (tmean > memoryLimit) {
      if (startAssess >= 0) 
        exitError("The memory limit is too low to allow initial assessment and summing"
                  " of %d frames in one pass %s", nzAlign, refineAtEnd || doSpline ? 
                  "with refinement or smoothing at the end" : 
                  "with this many pairwise comparisons");
      if (sumInOnePass) {
        sumInOnePass = false;
        if (!warnedTwoPass)
          printf("Using two passes: a single pass with %d frames requires %.1f GB, "
                 "above the limit of %.1f GB\n", nzAlign, tmean, memoryLimit);
        warnedTwoPass = true;
        fflush(stdout);
      }
    }
    numTestLoops = numBinTests + B3DCHOICE(sumInOnePass || testMode, 0, 1);

    // Loop on conditions
    wasGoodEnough = false;
    for (itest = 0; itest < numTestLoops; itest++) {

      // Set summing flag for summing with alignment if it can be done, otherwise set
      // for sum only on final loop unless assessing from subset, otherwise skip the sum
      if (sumInOnePass)
        summingMode = 0;
      else if (itest == numBinTests)
        summingMode = startAssess >= 0 ? 0 : -1;
      else
        summingMode = 1;

      // Set limits and indices for the extra loop; if it has to compute alignment because
      // it was assessed on a subset, then it needs either best filter or the whole set
      // to do hybrid
      useStart = numSets ? 0 : startFrame;
      useEnd = endFrame;
      if (itest == numBinTests) {
        indBinUse = indBestBin;
        if (hybridShifts) {
          indFiltUse = 0;
          numFiltUse = numFiltTests[itest];
        } else {
          indFiltUse = indBestFilt;
          numFiltUse = 1;
        }
      } else {

        // Otherwise set up for this round
        indBinUse = itest;
        indFiltUse = 0;
        numFiltUse = numFiltTests[itest];
        if (startAssess >= 0) {
          useStart = startAssess;
          useEnd = endAssess;
        }
      }
      if (ifDebug % 10)
        PRINT4(itest, summingMode, indBinUse, indFiltUse);

      // Set up actual frame limits for this file
      zDir = reverse ? -1 : 1;
      if (useStart > 0) {
        zStart = B3DCHOICE(reverse, B3DMIN(useEnd, nz) - 1, useStart - 1);
        zEnd = B3DCHOICE(reverse, useStart - 1, B3DMIN(useEnd, nz) - 1);
      } else {
        zStart = reverse ? nz - 1 : 0;
        zEnd = reverse ? 0 : nz - 1;
      }
      numFetch = zDir * (zEnd - zStart) + 1;
      nzAlign = B3DMAX(1, numFetch / blockGrpSize);

      // Get file, initialize
      if (useShrMem) {
#ifdef TEST_SHRMEM        
        ierr = sSMC.initialize
          (sumBin, binsToTest[indBinUse], trimFrac, numAVAuse, 
           refineAtEnd, hybridShifts, 
           (summingMode == 0 && (deferSum || doSpline)) ? 1 : 0,
           slideGrpSize, nx, ny,
           fullTaperFrac, taperFrac, antiFiltType - 1, radius1, 
           &varyRadius2[indBinUse][indFiltUse], sigma1, 
           &varySigma2[indBinUse][indFiltUse], numFiltUse, shiftLimit,
           kFactor, maxMaxWeight, summingMode, nzAlign, numOutFiles > 1,
           gpuFlags, ifDebug);
#endif        
      } else {
        ierr = sFA.initialize
          (sumBin, binsToTest[indBinUse], trimFrac, numAVAuse, 
           refineAtEnd, hybridShifts, 
           (summingMode == 0 && (deferSum || doSpline)) ? 1 : 0,
           slideGrpSize, nx, ny,
           fullTaperFrac, taperFrac, antiFiltType - 1, radius1, 
           &varyRadius2[indBinUse][indFiltUse], sigma1, 
           &varySigma2[indBinUse][indFiltUse], numFiltUse, shiftLimit,
           kFactor, maxMaxWeight, summingMode, nzAlign, numOutFiles > 1,
           gpuFlags, ifDebug);
      }
      if (ierr)
        exitError("Error %d initializing frame summing for file # %d", ierr, ifile + 1);

      // Initialize dose weighting: start by getting doses for all underlying frames
      if ((totalDose > 0. || doseFileType > 0) && !testMode) {
        sumOfFrames = 0;
        if (doseFileType > 3) {
          expandFrameDosesNumbers(frameDoseLines[ifile], &tempVal1[0], nz, sumOfFrames,
                                 sumOfDoses, &frameDoses[0]);
          if (sumOfFrames > 0 && sumOfFrames < nz)
            printf("WARNING: %s - The frame doses and numbers for file # %d include too "
                   "few frames (%d vs %d)", progname, ifile + 1, sumOfFrames, nz);
        }

        // Fall back to equal division
        if (sumOfFrames < nz)
          for (iz = 0; iz < nz; iz++)
            frameDoses[iz] = totalDoseVec[ifile] / nz;

        // Combine doses for grouped frames and frames actually being used
        tempVal1.resize(nzAlign);
        priorTemp = initialDose + priorDoseVec[ifile];
        if (useStart > 0)
          for (iz = 0; iz < useStart; iz++)
            priorTemp += frameDoses[iz];
        for (group = 0; group < nzAlign; group++) {
          frameGroupLimits(numFetch, nzAlign, group, groupStart, groupEnd, zStart, zDir,
                           izLow, izHigh);
          tempVal1[group] = 0.;
          for (iz = izLow; zDir * (iz - izHigh) <= 0; iz += zDir)
            tempVal1[group] += frameDoses[iz];
        }
        if (ifDebug % 10) {
          printf("Prior dose %.3f   total dose %.3f  frame doses:\n", priorTemp,
                 totalDoseVec[ifile]);
          for (iz = 0; iz < nzAlign; iz++) {
            printf(" %.3f", tempVal1[iz]);
            if ((iz + 1) % 12 == 0 || iz == nzAlign -1)
              printf("\n");
          }
        }
        ierr = sFA.setupDoseWeighting(priorTemp, &tempVal1[0], xScale, doseScaling,
                                      doseAfac, doseBfac, doseCfac, reweightFilt, iz);
        if (ierr)
          exitError("Error %d setting up dose weighting for file # %d", ierr, ifile + 1);
      }
    
      // Loop on frames in selected order, but do groups backwards to put largest at end
      for (group = 0; group < nzAlign; group++) {
        frameGroupLimits(numFetch, nzAlign, group, groupStart, groupEnd, zStart, zDir,
                       izLow, izHigh);
        useBuf = readBuf;
        useMode = inHead.mode;
        for (iz = izLow; zDir * (iz - izHigh) <= 0; iz += zDir) {
          if (combineFiles) {
            inFP = openAndReadHeader(inFiles[startCombine + iz], &inHead, "input image");
            checkInputFile(inFiles[startCombine + iz], &inHead, nx, ny, combineFiles);
            if (mrc_read_slice(readBuf, inFP, &inHead, 0, 'Z'))
              exitError("Reading from file # %d: %s", startCombine + iz, 
                        inFiles[startCombine + iz]);
            fclose(inFP);
          } else {
            izRead = iz;
            if (numSets)
              izRead += startCombine;

            // Read in data in parallel for TIFF file or with regular call
            wallStart = wallTime();
            if (parallelRead) {
              if (tiffParallelRead(fileCopies, numReadThreads, 0, nx - 1, 0, ny - 1,
                                   inDataSize, (char *)readBuf, izRead, MRSA_NOPROC))
                exitError("Reading frame %d from file # %d: %s", izRead, ifile + 1,
                          b3dGetError());
              
            } else {
              if (mrc_read_slice(readBuf, inFP, &inHead, izRead, 'Z'))
                exitError("Reading frame %d from file # %d", izRead, ifile + 1);
            }
            wallRead += wallTime() - wallStart;
          }
          if (izLow != izHigh) {
            useBuf = sumBuf;
            if (iz == izLow)
              memset(sumBuf, 0, needAlloc);
            addToSumBuffer(readBuf, inHead.mode, sumBuf, useMode, nx * ny);
          }
        }
        if (useShrMem) {
#ifdef TEST_SHRMEM        
          ierr = sSMC.nextFrame(useBuf, useMode, gainSlice ? gainSlice->data.f : NULL,
                                nxGain, nyGain, 
                                truncLimit, defectString, camSizeX, camSizeY, 
                                -1, bestXshifts[group], bestYshifts[group]);
#endif
        } else {
          ierr = sFA.nextFrame(useBuf, useMode, gainSlice ? gainSlice->data.f : NULL,
                               nxGain, nyGain, darkSlice ? darkSlice->data.f : NULL,
                               truncLimit, &defects, camSizeX, camSizeY, 
                               corDefBinning, bestXshifts[group], bestYshifts[group]);
        }
        if (ierr)
          exitError("Error %d processing frame/group %d from %s %d", ierr, group,
                    numSets ? "set" : "file", ifile + 1);
      }
      numDone = nzAlign;
      ierr = B3DMIN(numAVAuse, numDone) + 1 - groupSize;
      doRobust = ((ierr + 1 - groupSize) * (ierr - groupSize)) / 2 >= 2 * ierr &&
        kFactor > 0;

      // Finish up and get results; 
      if (useShrMem) {
#ifdef TEST_SHRMEM        
        float *alisum;
        ierr = sSMC.finishAlignAndSum
          (nx / sumBin, ny / sumBin, MAX_FILTERS + 1, 510, refRadius2, refSigma2,
           iterCrit, groupRefine, doSpline, &alisum, xShifts, yShifts,
                                     rawXshifts, rawYshifts, ringCorrs, frcDeltaR,
                                     faBestFilt, smoothDist, rawDist, resMean,
                                     resSD, meanResMax, maxResMax, meanRawMax, maxRawMax);
        memcpy(summed, alisum, 4 * (nx / sumBin) * (ny / sumBin));
#endif
      } else {
        ierr = sFA.finishAlignAndSum(refRadius2, refSigma2, iterCrit, groupRefine,
                                     doSpline, summed, xShifts, yShifts,
                                     rawXshifts, rawYshifts, ringCorrs, frcDeltaR,
                                     faBestFilt, smoothDist, rawDist, resMean,
                                     resSD, meanResMax, maxResMax, meanRawMax, maxRawMax);
      }
      if (ierr == 3)
        exitError("An unrecoverable error in GPU processing occurred for %s # %d",
                  numSets ? "set" : "file", ifile);
      else if (ierr)
        exitError("No frames were aligned for file # %d", ifile);
      if (itest < numBinTests) {
        numFetch = 1;
        if (numAVAuse && numFiltUse > 1)
          numFetch = numFiltUse + 1;
        if (!itest)
          printf("%s %d (%s): %d frames%s", numSets ? "Set" : "File", ifile + 1,
                 sstr.c_str(), zDir * (zEnd - zStart) + 1, numSets ? "" : "\n");
        if (!itest && numSets)
          printf(" from %d to %d\n", startCombine + 1, endCombine + 1);
        
        for (filt = 0; filt < numFetch; filt++) {
          if (B3DMIN(numAVAuse, numDone) >= 3) {
            if (numBinTests * numFiltTests[0] > 1) {
              if (filt < numFiltTests[itest] || numFiltTests[itest] == 1)
                printf("Results with bin = %d  rad2 = %.3f  sig2 = %.4f\n", 
                       binsToTest[itest], varyRadius2[itest][filt],
                       varySigma2[itest][filt]);
              else
                printf("Hybrid results,  bin = %d\n", binsToTest[itest]);
            } 
            printf("  %sesidual mean = %.3f, SD = %.3f, mean max = %.2f, max max "
                   "= %.2f\n", doRobust ? "Weighted r" : "R", resMean[filt],
                   resSD[filt], meanResMax[filt], maxResMax[filt]);
          } 
          if (doRobust && numAVAuse)
            printf("  Max unweighted resid mean = %.2f, max = %.2f ", 
                   meanRawMax[filt], maxRawMax[filt]);
          else
            printf("                                               ");
          printf("  Dist = %.2f, smoothed = %.2f\n", rawDist[filt], smoothDist[filt]);
          fflush(stdout);
        }
      }
        
      copyShifts = itest == numBinTests && summingMode == 0;
      error = resMean[faBestFilt] * (1. - maxMaxWeight) + 
        maxResMax[faBestFilt] * maxMaxWeight;
      if (error < minError && itest < numBinTests) {
        indBestBin = itest;
        indBestFilt = faBestFilt;
        minError = error;
        minMean = resMean[faBestFilt];
        copyShifts = true;
      }
      if (copyShifts) {
        if (ifDebug % 10)
          printf("Copy shifts test %d\n", itest);
        for (ix = 0; ix < numDone; ix++) {
          bestXshifts[ix] = xShifts[ix];
          bestYshifts[ix] = yShifts[ix];
          bestXraw[ix] = rawXshifts[ix];
          bestYraw[ix] = rawYshifts[ix];
          if (ifDebug % 10)
            printf("%.2f  %.2f\n", xShifts[ix], yShifts[ix]);
        }
      }

      // If the error is now good enough, advance to the end of the test runs
      if (minError <= goodEnough && itest < numBinTests - 1) {
        itest = numBinTests - 1;
        wasGoodEnough = true;
      }

    }  // End of test loop

    // Pick up the unweighted sum
    if (numOutFiles > 1) {
      if (sFA.getUnweightedSum(unwgtSum))
        exitError("getting non-dose-weighted sum");
    }

    // Advance set number and wrap it back to 0 after last file; close file when needed
    if (numSets) {
      setInFile++;
      if (setInFile == numSets) {
        superFile++;
        setInFile = 0;
      }
    }
    if (!combineFiles && !setInFile) {
      if (parallelRead)
        for (ind = 1; ind < numReadThreads; ind++)
          iiDelete(fileCopies[ind]);
      iiFClose(inFP);
      if (ifDebug && parallelRead)
        PRINT2(numReadThreads, wallRead);
    }

    // Write out data at end of loop
    if (summingMode <= 0) {
      if ((ifDebug % 10) > 1 && gettingFRC)
        for (ix = 0; ix < (int)floor(0.5 / frcDeltaR); ix++)
          printf("%.4f  %.5f\n", (ix + 0.5) * frcDeltaR, ringCorrs[ix]);

      for (outNum = 0; outNum < numOutFiles; outNum++) {
        headPtr = outHeads[outNum];
        useSum = outNum ? unwgtSum : summed;
        if (sumRotationFlip) {
          rotateFlipImage(useSum, MRC_MODE_FLOAT, nx / sumBin, ny / sumBin, 
                          sumRotationFlip, 0, 0, 0, rotSum, &ix, &iy, 0);
          useSum = rotSum;
        }

        // Trim if size is over
        ix = (nxSum - nxOut) / 2;
        iy = (nySum - nyOut) / 2;
        if (ix || iy)
          extractWithBinning(useSum, SLICE_MODE_FLOAT, nxSum, ix, ix + nxOut - 1, iy, 
                             iy + nyOut - 1, 1, useSum, 0, &iz, &ierr);
      
        // Scale the data and/or apply scale for binning
        if (scale != 1. || sumBin > 1)
          for (iz = 0; iz < headPtr->nx * headPtr->ny; iz++)
            useSum[iz] *= scale * sumBin * sumBin;
      
        // Manage header mmm and write the data
        arrayMinMaxMean(useSum, headPtr->nx, headPtr->ny, 0, headPtr->nx - 1, 0, 
                        headPtr->ny - 1, &tmin, &tmax, &tmean);
        ACCUM_MIN(headPtr->amin, tmin);
        ACCUM_MAX(headPtr->amax, tmax);
        headPtr->amean += tmean / numInFiles;
        ierr = mrcWriteZFloat(headPtr, &li, useSum, ifile);
        if (ierr)
          exitError("Writing summed data to file for input file # %d (error # %d)", 
                    ifile + 1, ierr);
        if (adjustMdoc && !outNum) {
          if (AdocGetFloat(ADOC_ZVALUE_NAME, ifile, "PixelSpacing", &pixTemp))
            exitError("Getting pixel spacing in mdoc for output image");
          if (AdocSetFloat(ADOC_ZVALUE_NAME, ifile, "PixelSpacing", pixTemp * relBinning))
            exitError("Adjusting pixel spacing in mdoc for output image");
          if (AdocSetThreeFloats(ADOC_ZVALUE_NAME, ifile, "MinMaxMean", tmin, tmax,tmean))
            exitError("Adjusting pixel spacing in mdoc for output image");

          // Binning was a later addition so allow it not to exist
          ierr = AdocGetFloat(ADOC_ZVALUE_NAME, ifile, "Binning", &pixTemp);
          if (ierr < 0)
            exitError("Getting binning in mdoc for output image");
          if (!ierr) {
            pixTemp *= relBinning;
            if (pixTemp > 0.55)
              pixTemp = B3DNINT(pixTemp);
            if (AdocSetFloat(ADOC_ZVALUE_NAME, ifile, "Binning", pixTemp))
              exitError("Adjusting binning in mdoc for output image");
          }
          if (ignoreZvalue) {
            sprintf(inLine, "%d", ifile);
            if (AdocChangeSectionName(ADOC_ZVALUE_NAME, ifile, inLine))
              exitError("Changing section name to new Z value");
          }
        }
      }
    }

    if (numBinTests * numFiltTests[0] > 1) {
        printf("File %d: %s at bin = %d  rad2 = %.3f  sig2 = %.4f  mean res = %.3f\n",
               ifile + 1, wasGoodEnough ? "Good enough" : "Best", binsToTest[indBestBin],
               varyRadius2[indBestBin][indBestFilt],
               varySigma2[indBestBin][indBestFilt], minMean);
        numTimesBest[indBestBin][indBestFilt]++;
    }

    // Find FRC crossing and mean near half-nyquist
    if (!testMode && gettingFRC) {
      if (useShrMem) {
#ifdef TEST_SHRMEM        
        sSMC.analyzeFRCcrossings(510, ringCorrs, frcDeltaR, halfCross, quartCross,
                                 eighthCross, halfNyq);
#endif
      } else {
        sFA.analyzeFRCcrossings(ringCorrs, frcDeltaR, halfCross, quartCross, eighthCross,
                                halfNyq);
      }
      printf(" FRC crossings 0.5: %.4f  0.25: %.4f  0.125: %.4f  is %.4f at "
             "0.25/pix\n", halfCross, quartCross, eighthCross, halfNyq);
    }
    if (ifile < numInFiles - 1 && numBinTests * numFiltTests[0] > 1)
      printf("\n");

    // Output the transforms
    if (xfExt) {
      xfName = filename;
      tind = xfName.find_last_of('.');
      if (tind >= xfName.length() - 5 && tind > 1 && tind != string::npos)
        xfName.resize(tind + 1);
      xfName += xfExt;
      imodBackupFile(xfName.c_str());
      extraFP = fopen(xfName.c_str(), "w");
      if (!extraFP)
        exitError("Opening output file %s for transforms", xfName.c_str());
      /*if (startFrame > 0)
        for (iz = 1; iz < startFrame; iz++)
          fprintf(extraFP, " 1.00000    0.00000    0.00000   1.00000     0.000   "
          " 0.000\n"); */
      /*for (iz = zStart; zDir * (zEnd - iz) >= 0; iz += zDir) {
        ind = startFrame > 0 ? iz + 1 - startFrame : iz;
        fprintf(extraFP, " 1.00000    0.00000    0.00000   1.00000  %8.3f %8.3f\n",
                bestXshifts[ind], bestYshifts[ind]);
                }*/
      numFetch = zDir * (zEnd - zStart) + 1;
      for (group = 0; group < nzAlign; group++) {
        if (reverse) {
          balancedGroupLimits(numFetch, nzAlign, group, &groupStart, 
                              &groupEnd);
          ind = nzAlign - 1 - group;
        } else {
          balancedGroupLimits(numFetch, nzAlign, nzAlign - 1 - group, &groupStart, 
                              &groupEnd);
          ind = group;
        }
        if (startFrame > 0 && !group)
          groupEnd += startFrame - 1;
        for (iz = groupStart; iz <= groupEnd; iz++)
          fprintf(extraFP, " 1.00000    0.00000    0.00000   1.00000  %8.3f %8.3f\n",
                  bestXshifts[ind], bestYshifts[ind]);
      }
      fclose(extraFP);
    }
    
    // Output plottable shifts
    if (plotFP) {
      for (ind = 0; ind < numDone; ind++)
        fprintf(plotFP, "%3d  %.3f  %.3f\n", 10 * ifile + 10, bestXraw[ind],
                bestYraw[ind]);
      if (doSpline)
        for (ind = 0; ind < numDone; ind++)
          fprintf(plotFP, "%3d  %.3f  %.3f\n", 10 * ifile + 11, bestXshifts[ind], 
                  bestYshifts[ind]);
    }

    // Output the FRC
    if (frcFP)
      for (ind = 0; ind < (int)floor(0.5 / frcDeltaR); ind++)
        fprintf(frcFP, "%2d  %.4f %10.6f\n", ifile + 1, (ind + 0.5) * frcDeltaR,
                ringCorrs[ind]);
  }

  if (numInFiles > 1 && numBinTests * numFiltTests[0] > 1) {
    printf("\nNumber of times each condition is best  (rad2 in parentheses):\n");
    for (itest = 0; itest < numBinTests; itest++) {
      printf("bin = %d  ", binsToTest[itest]);
      for (filt = 0; filt < numFiltTests[itest]; filt++)
        printf("  %3d (%.3f)", numTimesBest[itest][filt], varyRadius2[itest][filt]);
      printf("\n");
    }
  }

  if (adjustMdoc) {
    sstr = outNames[0];
    sstr += ".mdoc";
    if (AdocWrite(sstr.c_str()))
      exitError("Writing adjusted mdoc file");
  }

  // Finish up
  if (useShrMem) {
#ifdef TEST_SHRMEM        
    sSMC.cleanup();
    sSMC.Disconnect();
#endif
  } else {
    sFA.cleanup();
  }
  free(summed);
  B3DFREE(rotSum);
  B3DFREE(unwgtSum);
  free(savedListName);
  if (!testMode) {
    for (outNum = 0; outNum < numOutFiles; outNum++) {
      if (mrc_head_write(outFPs[outNum], outHeads[outNum]))
        exitError("Writing header to output file");
      iiFClose(outFPs[outNum]);
      if (openTsNames[outNum])
        remove(openTsNames[outNum]);
      free(outNames[outNum]);
      free(openTsNames[outNum]);
    }
  }
  if (frcFP)
    fclose(frcFP);
  exit(0);
}

// open an MRC file and read its header
FILE *AliFrame::openAndReadHeader(const char *filename, MrcHeader *head,
                                  const char *descrip, bool testMode)
{
  FILE *inFP = iiFOpen(filename, "rb");
  if (!inFP)
    exitError("Opening %s file %s%s", descrip, filename, 
              testMode ? "; do not specify an output file when not making sums" : "");
  if (mrc_head_read(inFP, head))
    exitError("Reading header of %s file %s", descrip, filename);
  return inFP;
}

void AliFrame::checkInputFile(const char *filename, MrcHeader *head, int nx, int ny, 
                              int combine)
{
  if (sliceModeIfReal(head->mode) < 0)
    exitError("File mode for %s is %d; only byte, short, float allowed", filename, 
              head->mode);
  if (nx > 0 && (nx != head->nx || ny != head->ny))
    exitError("File %s has a different size (%d x %d) from previous files (%d x %d)",
              filename, head->nx, head->ny, nx, ny);
  if (combine && head->nz > 1)
    exitError("File %s has more than one slice and cannot be used with -combine",
              filename);
}

// Add a read-in image to a sum buffer of the proper type
void AliFrame::addToSumBuffer(void *readBuf, int inMode, void *sumBuf, int &useMode, 
                              int nxy)
{
  unsigned char *bData = (unsigned char *)readBuf;
  short *sData = (short *)readBuf;
  short *sBuf = (short *)sumBuf;
  unsigned short *usBuf = (unsigned short *)sumBuf;
  unsigned short *usData = (unsigned short *)readBuf;
  float *fBuf = (float *)sumBuf;
  float *fData = (float *)readBuf;
  int ix;
  switch (inMode) {
  case MRC_MODE_BYTE:
    useMode = SLICE_MODE_SHORT;
    for (ix = 0; ix < nxy; ix++)
      *sBuf++ += *bData++;
    break;
  case MRC_MODE_SHORT:
    useMode = SLICE_MODE_SHORT;
    for (ix = 0; ix < nxy; ix++)
      *sBuf++ += *sData++;
    break;
  case MRC_MODE_USHORT:
    useMode = SLICE_MODE_USHORT;
    for (ix = 0; ix < nxy; ix++)
      *usBuf++ += *usData++;
    break;
  case MRC_MODE_FLOAT:
    useMode = SLICE_MODE_FLOAT;
    for (ix = 0; ix < nxy; ix++)
      *fBuf++ += *fData++;
    break;
  }
}

void AliFrame::extractFileTail(const char *filename, std::string &sstr)
{
  sstr = filename;
  size_t iz = sstr.find_last_of("/\\");
  if (iz == string::npos)
    iz = 0;
  else
    iz++;
  sstr = sstr.substr(iz);
}

// Checks if the given title has the binning value in it, adjusts by relBinning and places
// into sstr, and gives message if printChange true.  Returns 1 if binning found
int AliFrame::adjustTitleBinning(const char *title, std::string &sstr, float relBinning,
                                 bool printChange)
{
  int ierr, stackBin;
  char binText[32];
  size_t iz;
  sstr = title;

  // Find binning = string
  iz = sstr.find("binning =");
  if (iz == 0 || iz == string::npos)
    return 0;
  iz += 9;
    
  // Find beginning and end of the binning number
  while (title[iz] == ' ')
    iz++;
  if (title[iz] == 0x00)
    return 0;
  ierr = iz;
  while (title[iz] != ' ' && title[iz] != 0x00)
    iz++;
  stackBin = atoi(&title[ierr]);

  // Make sure it reads, then adjust and make sure that is a sensible value
  if (!stackBin)
    return 0;
  relBinning *= stackBin;
  if (relBinning < 0.46 || (relBinning > 0.54 && relBinning < 0.96) || 
      (relBinning > 0.96 && fabs(B3DNINT(relBinning) - relBinning) > 0.04))
    return 0;

  // Replace with new value
  if (relBinning < 0.54) {
    sprintf(binText, "%.1f", relBinning);
    if (printChange)
      printf("Adjusted binning in header title to %.1f\n", relBinning);
  } else {
    sprintf(binText, "%d", B3DNINT(relBinning));
    if (printChange)
      printf("Adjusted binning in header title to %d\n", B3DNINT(relBinning));
  }
  sstr.replace(ierr, iz - ierr, binText);
  if (sstr[0] != ' ')
    sstr.insert(0, "    ");
  return 1;
}

/*
 * Determine actual minimum and maximum size for sets or groups given the nominal size and
 * the total to be divided into the groups.
 */
void AliFrame::minMaxSetSize(int basicSize, int numFrames, int &minSize, int &maxSize)
{
  int numSets = numFrames / basicSize;
  int remainder = numFrames % basicSize;
  minSize = basicSize + remainder / numSets;
  maxSize = minSize + (remainder % numSets > 0 ? 1 : 0);
}

/*
 * Analyze the extra header from a UCSFtomo file for pixel size, rotation angle, and
 * tilt angles and determine division into groups by tilt angle
 */
int AliFrame::analyzeExtraHeader(FILE *inFP, MrcHeader *head, int startFrame, 
                                 int endFrame, bool axisPixOnly, float **buffer,
                                 int &bufSize, FloatVec &tilts, IntVec &setStarts,
                                 int &minSet, int &maxSet, float &axisAngle,
                                 float &pixSize)
{
  int numIntReal = head->nint + head->nreal;
  int num = numIntReal * head->nz;
  int iz, set, izStart, izEnd, size, lastSize = 0, numAtSize, retVal = 0, gotDouble = 0;
  float temp, lastTilt;
  if (!num)
    return 0;
  if (bufSize < num) {
    B3DFREE(*buffer);
    *buffer = B3DMALLOC(float, num);
    if (!*buffer)
      exitError("Allocating buffer for extended header data");
    bufSize = num;
  }
  if (b3dFseek(inFP, MRC_HEADER_SIZE, SEEK_SET) || 
      (int)b3dFread(*buffer, 4, num, inFP) != num)
    exitError("Reading extended header data");

  // Look for a legal pixel size and rotation angle
  if (head->nreal >= 12) {
    pixSize = 0.;
    temp = (*buffer)[head->nint + 11];

    // UCSF tomo puts out angstroms, but it could still be meters.  These are limits for
    // angstroms in flib/image/header.f90
    if (temp > 0.05 && temp < 100000.) {
      pixSize = temp;
    } else {
      temp *= 1.e10;
      if (temp > 0.05 && temp < 100000.)
        pixSize = temp;
    }
    temp = (*buffer)[head->nint + 10];

    // set return to 1 if both are legal
    if (pixSize > 0. && temp >= -360. && temp <= 360.) {
      retVal = 1;
      if (temp < -180.)
        temp += 360.;
      if (temp > 180.)
        temp -= 360.;
      axisAngle = temp;
    }
  }
    
  if (axisPixOnly)
    return retVal;

  // Now analyze tilt angles in slot 1.  Loop on the subset of frames if any
  izStart = 0;
  izEnd = head->nz - 1;
  if (startFrame > 0) {
    izStart = startFrame - 1;
    izEnd = endFrame - 1;
  }
  tilts.clear();
  setStarts.clear();

  // Look for each place where tilt angle changes and save the angle and set start
  for (iz = izStart; iz <= izEnd; iz++) {
    temp = (*buffer)[head->nint + iz * numIntReal];
    if (temp < -180. || temp > 180.)
      return retVal;
    if (iz == izStart || fabs(temp - lastTilt) > 0.01) {
      if (iz > izStart) {

        // Keep track of size of sets and number of ones at the last size
        // If there have been at least 5 in a row at a size and there is one at twice 
        // the size, it must be the repeated one at the starting angle
        size = iz - setStarts.back();
        if (numAtSize > 5 && size == 2 * lastSize && !gotDouble) {
          setStarts.push_back(iz - lastSize);
          tilts.push_back(temp);
          size = lastSize;
          gotDouble = 1;
        }
        if (size == lastSize) {
          numAtSize++;
        } else {
          numAtSize = 1;
          lastSize = size;
        }
      }
      setStarts.push_back(iz);
      tilts.push_back(temp);
      lastTilt = temp;
    }
  }
  setStarts.push_back(iz);

  // Get the min and max set size
  for (iz = 0; iz < (int)tilts.size(); iz++) {
    set = setStarts[iz + 1] - setStarts[iz];
    if (!iz) 
      minSet = maxSet = set;
    ACCUM_MIN(minSet, set);
    ACCUM_MAX(maxSet, set);
  }
  return 2 + retVal;
}

/*
 * Apply rotation and flip operation to gain reference: now it's simple
 */
void AliFrame::rotateFlipGainReference(float *ref, int &nxGain, int &nyGain, 
                                       int rotationFlip)
{
  int nxIn = nxGain, nyIn = nyGain;
  float *summed = B3DMALLOC(float, nxGain * nyGain);
  if (!summed)
    exitError("Allocating array for rotating gain reference");
  if (rotateFlipImage(ref, MRC_MODE_FLOAT, nxIn, nyIn, rotationFlip, 0, 0, 0, summed,
                      &nxGain, &nyGain, 0))
    exitError("Inappropriate rotation/flip value %d entered", rotationFlip);
  memcpy(ref, summed, nxGain * nyGain * sizeof(float));
  free(summed);
}

/*
 * Common operations when opening mdoc for either option
 */
int AliFrame::openMdocFile(const char *filename, int &numSect, int &adocType)
{
  int adocInd, ind;
  adocInd = AdocOpenImageMetadata(filename, 0, &ind, &numSect, &adocType);
  if (adocInd == -1)
    exitError("Opening or reading mdoc file %s", filename);
  if (adocInd == -2)
    exitError("Metadata file %s does not exist", filename);
  if (adocInd == -3)
    exitError("Metadata file %s does not have image stack information", filename);
  return adocInd;
}

/*
 * Convert a test line for frame doses and number into 
 */
void AliFrame::expandFrameDosesNumbers(std::string &line, float *tempArr, int maxFrames,
                                      int &totalFrames, float &totalDose,
                                      float *frameDoses)
{
  int ind, add, numSame, numToGet = 0;
  totalFrames = 0;
  totalDose = 0.;
  if (line.empty())
    return;
  if (PipGetLineOfValues(FRAME_DOSE_KEY, line.c_str(), tempArr, PIP_FLOAT, 
                         &numToGet, 2 * maxFrames))
    exitError("Processing an entry for frame doses and numbers: %s", line.c_str());
  if (numToGet % 2)
    exitError("Odd number of numbers in entry for frame doses and numbers: %s", 
              line.c_str());
  for (ind = 0; ind < numToGet; ind += 2) {
    numSame = B3DNINT(tempArr[ind + 1]);
    if (fabs(numSame - tempArr[ind + 1]) > 1.e-3)
      exitError("Non-integer value for count in entry for frame doses and numbers: %s",
                line.c_str());
    if (totalFrames + numSame > maxFrames)
      exitError("Frame numbers add up to more than maximum expected number of frames "
                "in: %s", line.c_str());
    totalDose += numSame * tempArr[ind];
    if (frameDoses)
      for (add = 0; add < numSame; add++)
        frameDoses[totalFrames++] = tempArr[ind];
    else
      totalFrames += numSame;
  }
}

/*
 * Compute limits for looping over frame groups.  Yes, class members would be easier
 */
void AliFrame::frameGroupLimits(int numFetch, int nzAlign, int group, int &groupStart, 
                              int &groupEnd, int zStart, int zDir, int &izLow, 
                              int &izHigh)
{
  int iz;
  balancedGroupLimits(numFetch, nzAlign, nzAlign - 1 - group, &groupStart, &iz);
  groupEnd = numFetch - 1 - groupStart;
  groupStart = numFetch - 1 - iz;
  izLow = zStart + zDir * groupStart;
  izHigh = zStart + zDir * groupEnd;
}

#ifdef TEST_SHRMEM

// Function to get a string from the defect file that can be passed and processed
static void defectFileToString(char *name, std::string &outStr)
{
  int len;
  bool endReached = false;
  FILE *fp = fopen(name, "r");
  char buf[MAX_LINE];
  outStr = "";
  if (!fp)
    exitError("Failed to open defect file");
  while (true) {
    if (endReached)
      break;
    len = fgetline(fp, buf, MAX_LINE);

    // interpret the len value
    if (!len)
      continue;
    if (len == -2)
      break;
    if (len == -1)
      exitError("Error reading defect file");
    if (len < 0) {
      endReached = true;
      len = -len - 2;
    }
    outStr += buf;
    outStr += "\n";
    if (endReached)
      break;
  }
  fclose(fp);
}
#endif
