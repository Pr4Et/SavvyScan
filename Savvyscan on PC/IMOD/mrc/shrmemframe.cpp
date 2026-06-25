/*
 *  shrmemframe - Frame alignment program with interface through Windows shared memory
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 2019 by the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 *  $Id$
 *
 * To build:
 * setup -no_tiff -no_shared -mkl
 * make -j4 clibs
 * cd mrc
 * make -j4 shrmemframe.exe
 * addmanifest shrmemframe.exe
 *
 * Define SHRMEMFRAME_LOGDIR for logging
 * 
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Windows.h>
#include "b3dutil.h"
#include "cppdefs.h"
#include "shrmemframe.h"
#include "framealign.h"
#include "frameutil.h"
#include <time.h>

#ifndef _delayimp_h
extern "C" IMAGE_DOS_HEADER __ImageBase;
#endif

#define INT_SIZE sizeof(int)

// Pointer/counter for storing messages to pass back from call
static char *sMessageBuf;
static int sCumMessLength;

static FILE *sLogFP = NULL;

// Save a message in the buffer
static void framePrintFunc(const char *strMessage)
{
  int messLen;
  if (!sMessageBuf)
    return;
  messLen = (int)strlen(strMessage);
  if (sCumMessLength + messLen >= MESS_BUF_SIZE - 1)
    return;
  strcpy(sMessageBuf + sCumMessLength, strMessage);
  sCumMessLength += messLen;
}

// Tiny main; instantiate the class
int main( int argc, char *argv[])
{
  ShrMemFrame sharer;
  sharer.main(argc, argv);
  exit(0);
}

ShrMemFrame::ShrMemFrame()
{
  mMapFile = NULL;
  mMappedBuf = NULL;
  mActionSignal = NULL;
  mDoneSignal = NULL;
}

/*
 * Main method
 */
void ShrMemFrame::main(int argc, char *argv[])
{
  int reply, funcCode, idVal, ind;
  int sleepTime = 100;
  char buffer[MAX_PATH];
  TCHAR dllPath[MAX_PATH];
  HMODULE thisModule = reinterpret_cast<HMODULE>(&__ImageBase);
  std::string hereStr;

  // Get the ID value, needed to prevent orphaned copies from responding
  mIdValue = -999;
  if (argc > 1)
    mIdValue = atoi(argv[1]);

  // Set up logging
  if (getenv("SHRMEMFRAME_LOGDIR")) {
    sprintf(buffer, "%s/shrmemframe%d.log", getenv("SHRMEMFRAME_LOGDIR"), mIdValue);
    sLogFP = fopen(buffer, "w");
    if (!sLogFP)
      fprintf(sLogFP, "starting log\n"); fflush(sLogFP);
  }

  // Add the upper directory to the DLL search path so FrameGPU can be found
  if (GetModuleFileName(thisModule, &dllPath[0], MAX_PATH)) {
    hereStr = dllPath;
    ind = hereStr.find_last_of("/\\");
    if (ind > 0) {
      hereStr.resize(ind);
      ind = hereStr.find_last_of("/\\");
      if (ind > 0) {
        hereStr.resize(ind);
        SetDllDirectory(hereStr.c_str());
      }
    }
  }
  
  // Set size of shared memory: client uses same formula
  mShMemSize = B3DMAX(MAX_FOR_RAW_BYTE, MAX_FOR_RAW_SHORT) + ALL_ADDED_SPACE;

  // Create the file mapping
  mMapFile = CreateFileMapping(INVALID_HANDLE_VALUE,    // use paging file
                               NULL,                    // default security 
                               PAGE_READWRITE,          // read/write access
                               0,                       // max. object size 
                               mShMemSize,              // buffer size  
                               MAPPING_NAME);           // name of mapping object
  if (mMapFile == NULL) {
    if (sLogFP)
      fprintf(sLogFP, "Could not create file mapping object (%d).\n", GetLastError());
    ExitProcess(1);
  }
   
  // Get the buffer address
  mMappedBuf = (char *)MapViewOfFile(mMapFile,   // handle to map object
                                     FILE_MAP_ALL_ACCESS, // read/write permission
                                     0,                   
                                     0,                   
                                     mShMemSize);
  if (!mMappedBuf) {
    if (sLogFP)
      fprintf(sLogFP, "Could not map view of memory object (%d).\n", GetLastError());
    closeAndExit(2);
  }

  sprintf(buffer, "%s%d", ACTION_BASE_NAME, mIdValue);
  mActionSignal =  CreateEventA(NULL, FALSE, FALSE, (LPCSTR)buffer);
  if (!mActionSignal) {
    if (sLogFP)
      fprintf(sLogFP, "Could create action event (%d).\n", GetLastError());
    closeAndExit(3);
  }
  sprintf(buffer, "%s%d", DONE_BASE_NAME, mIdValue);
  mDoneSignal =  CreateEventA(NULL, FALSE, FALSE, (LPCSTR)buffer);
  if (!mDoneSignal) {
    if (sLogFP)
      fprintf(sLogFP, "Could create done event (%d).\n", GetLastError());
    closeAndExit(4);
  }

  // Send the event that we are ready to go
  SetEvent(mDoneSignal);
  mFA.setPrintFunc(framePrintFunc);

  // Start infinite loop waiting for signals
  int lastFunc = 0;
  for (;;) {

    // Every call includes ID value, function/action code, and possible arguments
    if (WaitForSingleObject(mActionSignal, sleepTime) != WAIT_OBJECT_0) {
      memcpy(&idVal, mMappedBuf, INT_SIZE);

      // If you see someone else's ID in the shared buffer, close up and quit
      if (idVal != mIdValue)
        closeAndExit(0);
      continue;
    }

    memcpy(&funcCode, mMappedBuf + INT_SIZE, INT_SIZE);

    // Respond to the code
    switch (funcCode) {
    case SHMFR_DONE:
      sendReply(SHMFR_DONE);    // Avoid deadlocks by signalling done on nonsense
      break;

    case SHMFR_VERSION:
      reply = SHRMEMFRAME_VERSION;
      memcpy(mMappedBuf + 2 * INT_SIZE, &reply, INT_SIZE);
      sendReply(SHMFR_DONE);
      break;
      
    case SHMFR_EXIT:
      sendReply(SHMFR_DONE);
      closeAndExit(0);
      break;

    case SHMFR_INITIALIZE:
      initialize();
      break;

    case SHMFR_NEXTFRAME:
      nextFrame();
      break;

    case SHMFR_FINISH:
      finishSum();
      break;

    case SHMFR_FRC_CROSS:
      analyzeFRC();
      break;

    case SHMFR_GPU_AVAIL:
      gpuAvailable();
      break;

    case SHMFR_CLEANUP:
      sMessageBuf = NULL;
      mFA.cleanup();
      sendReply(SHMFR_DONE);
      break;

    default:
      sendReply(SHMFR_DONE);
      break;
    }
  }

  // This should not be reached
  closeAndExit(0);
}

// Boilerplate items for the functions
#define GET_PARAMS(type)                                    \
  type params;                                              \
  memcpy(&params, mMappedBuf + 2 * INT_SIZE, sizeof(type));

#define INIT_MESS_BUF                                    \
  sMessageBuf = &params.messages[0];                     \
  sCumMessLength = 0;                                    \
  params.messages[0] = 0x00;

#define PUT_PARAMS(type)                               \
  memcpy(mMappedBuf + 2 * INT_SIZE, &params, sizeof(type)); \
  sendReply(SHMFR_DONE);

/*
 * Initialize
 */
void ShrMemFrame::initialize()
{
  GET_PARAMS(InitializeParams);
  INIT_MESS_BUF;
  params.retVal = mFA.initialize
    (params.binSum, params.binAlign, params.trimFrac, params.numAllVsAll,
     params.cumAlignAtEnd, params.useHybrid, params.deferSum,
     params.groupSize, params.nx, params.ny, params.padFrac, params.taperFrac,
     params.antiFiltType, params.radius1, &params.radius2[0], params.sigma1,
     &params.sigma2[0], params.numFilters, params.maxShift, params.kFactor,
     params.maxMaxWeight, params.summingMode, params.expectedZ,
     params.makeUnwgtSum, params.gpuFlags, params.debug);
  PUT_PARAMS(InitializeParams);
  mNx = params.nx;
  mNy = params.ny;
  mNumFrames = 0;
  mCamSizeX = mCamSizeY = 0;
}

/*
 * Next frame
 */
void ShrMemFrame::nextFrame()
{
  void *frame;
  float *gainRef = NULL;
  char *defectStr;
  int temp;
  GET_PARAMS(NextFrameParams);
  INIT_MESS_BUF;
  params.retVal = 0;

  // Set pointers to data and gain ref
  frame = mMappedBuf + params.frameOffset;
  if (params.gainRefOffset)
    gainRef = (float *)(mMappedBuf + params.gainRefOffset);

  // Get the defect information on the first frame only
  if (params.defectsOffset && !mNumFrames) {
    defectStr = mMappedBuf + params.defectsOffset;

    // Process the defect string if it has changed, just as in plugin
    if (!mDefectString.length() || mDefectString.compare(defectStr)) {
      try {
        mDefectString = defectStr;
        mDefCorBinning = params.defBin;
        if (CorDefParseDefects(mDefectString.c_str(), 1, mCamDefects,
                               mCamSizeX, mCamSizeY)) {
          params.retVal = 5;
        } else {
          if (params.defBin < 0) {

            // Test call from alignframes: process as there
            mDefCorBinning = 1;
            CorDefFlipDefectsInY(&mCamDefects, mCamSizeX, mCamSizeY, 0);
            CorDefFindTouchingPixels(mCamDefects, mCamSizeX, mCamSizeY, 0);
            CorDefSetupToCorrect(mNx, mNy, mCamDefects, mCamSizeX, mCamSizeY, 
                                 0, -1., mDefCorBinning, NULL);
          } else {

            // Otherwise process as in plugin
            // The plugin sets binning to 1 for counting and 2 for super-res
            // because it scales up to SR size regardless, so use the value directly
            mDefCorBinning = params.defBin;
            CorDefFindTouchingPixels(mCamDefects, mCamSizeX, mCamSizeY, 0);

            if (!mCamDefects.wasScaled)
              CorDefScaleDefectsForK2(&mCamDefects, false);
            if (mCamDefects.rotationFlip % 2) {
              temp = mCamSizeX;
              mCamSizeX = mCamSizeY;
              mCamSizeY = temp;
            }
            CorDefRotateFlipDefects(mCamDefects, 0, mCamSizeX, mCamSizeY);
          }
        }
      }
      catch (...) {
        params.retVal = 4;
      }
    }

    // return if error
    if (params.retVal) {
      PUT_PARAMS(NextFrameParams);
      return;
    }
  }

  // Make the call
  params.retVal = mFA.nextFrame(frame, params.type, gainRef, params.nxGain, params.nyGain,
                                NULL, params.truncLimit, &mCamDefects, mCamSizeX,
                                mCamSizeY, mDefCorBinning, params.shiftX,
                                params.shiftY);
  PUT_PARAMS(NextFrameParams);
  mNumFrames++;
}

/*
 * Finish align and sum
 */
void ShrMemFrame::finishSum()
{
  float *xShifts, *yShifts, *rawXshifts, *rawYshifts, *alisum;
  GET_PARAMS(FinishAlignAndSumParams);
  INIT_MESS_BUF;

  // Set pointers to return image and shifts
  alisum = (float *)(mMappedBuf + params.alisumOffset);
  xShifts = (float *)(mMappedBuf + params.xShiftsOffset);
  yShifts = (float *)(mMappedBuf + params.yShiftsOffset);
  rawXshifts = (float *)(mMappedBuf + params.rawXshiftsOffset);
  rawYshifts = (float *)(mMappedBuf + params.rawYshiftsOffset);

  // Make the call
  params.retVal = mFA.finishAlignAndSum
    (params.refineRadius2, params.refineSigma2, params.iterCrit, params.groupRefine,
     params.doSpline, alisum, xShifts, yShifts, rawXshifts, rawYshifts,
     &params.ringCorrs[0], params.deltaR, params.bestFilt,
     &params.smoothDist[0], &params.rawDist[0], &params.resMean[0], &params.resSD[0],
     &params.meanResMax[0], &params.maxResMax[0],
     &params.meanRawMax[0], &params.maxRawMax[0]);
  PUT_PARAMS(FinishAlignAndSumParams);
}

/*
 * Analyze FRC crossings
 */
void ShrMemFrame::analyzeFRC()
{
  GET_PARAMS(AnalyzeFRCcrossingsParams);
  sMessageBuf = NULL;
  mFA.analyzeFRCcrossings(&params.ringCorrs[0], params.frcDelta, params.halfCross, 
                          params.quartCross, params.eighthCross, params.halfNyq);
  PUT_PARAMS(AnalyzeFRCcrossingsParams);
}

/*
 * Is GPU available here
 */
void ShrMemFrame::gpuAvailable()
{
  GET_PARAMS(GpuAvailableParams);
  INIT_MESS_BUF;
  params.retVal = mFA.gpuAvailable(params.nGPU, &params.memory, params.debug);
  PUT_PARAMS(GpuAvailableParams);
}

/*
 * Put out the ID and the code for every reply
 */
void ShrMemFrame::sendReply(int reply)
{
  if (sLogFP)
    fprintf(sLogFP,  "Sending reply %d  %d\n", mIdValue, reply);fflush(sLogFP);
  memcpy(mMappedBuf, &mIdValue, INT_SIZE);
  memcpy(mMappedBuf + INT_SIZE, &reply, INT_SIZE);
  SetEvent(mDoneSignal);
}

/*
 * Close the mapping and exit with the code
 */
void ShrMemFrame::closeAndExit(int exitCode)
{
  if (mMappedBuf)
    UnmapViewOfFile(mMappedBuf);
  if (mMapFile)
    CloseHandle(mMapFile);
  if (mActionSignal)
    CloseHandle(mActionSignal);
  if (mDoneSignal)
    CloseHandle(mDoneSignal);
  ExitProcess(exitCode);
}
