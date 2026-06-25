/*
 *  ShrMemClient - interface for a client (alignframes or SEMCCD) to do frame 
 *  alignment through Windows shared memory to shrmemframe
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 2019 by the Regents of the University of 
 *  Colorado.  See IMOD/dist/COPYRIGHT for full copyright notice.
 *
 */
#include "stdafx.h"
#include <sys/stat.h>
#include <time.h>
#include "ShrMemClient.h"
#ifdef TEST_SHRMEM
#include "shrmemframe.h"
#else
#include "Shared/shrmemframe.h"
#include "TemplatePlugIn.h"
#endif

#ifndef _delayimp_h
extern "C" IMAGE_DOS_HEADER __ImageBase;
#endif

#define INT_SIZE sizeof(int)

#ifdef TEST_SHRMEM

// Definitions in the alignframes test case
double TickInterval(double start)
{
  double interval = GetTickCount() - start;
  if (interval < 0)
    interval += 4294967296.;
  return interval;
}
#define DEB_PRINT printf
#define PRN_FMT "%s",
#define PREFIX
#define ERR_PRINT printf
#else

// Definitions in the plugin case
extern PlugInWrapper gPlugInWrapper;
#define DEB_PRINT gPlugInWrapper.DebugToResult
#define PRN_FMT 
#define PREFIX ,"SerialCCD ShrMemClient: "
#define ERR_PRINT gPlugInWrapper.ErrorToResult
#endif

// Constructor: initialize items that may be tested for closing
ShrMemClient::ShrMemClient(void)
{
  mShMemSize = B3DMAX(MAX_FOR_RAW_BYTE, MAX_FOR_RAW_SHORT) + ALL_ADDED_SPACE;
  mMappedBuf = NULL;
  mMapFile = NULL;
  mProcessHandle = NULL;
  mActionSignal = NULL;
  mDoneSignal = NULL;
}

// Destructor: try to close the server
ShrMemClient::~ShrMemClient(void)
{
  Disconnect(50);
}

// Call to check if connection exists and make it if not
int ShrMemClient::ConnectIfNeeded(void)
{
  char buffer[160];
  STARTUPINFO sinfo;
  PROCESS_INFORMATION pinfo;
  std::string searchDirs[5], shrMemPath;
  struct _stat statbuf;
  int ind, version;
  char *argString;
  TCHAR dllPath[MAX_PATH];
  int initialTimeout = 5000, initTryInterval = 50;
  double startTime;
  HMODULE thisModule = reinterpret_cast<HMODULE>(&__ImageBase);
  
  if (mMappedBuf && !CheckIfProcessDied(
    "shrmemframe is no longer running; trying to reconnect\n"))
    return 0;

  // Get a new ID value
  mServerID = time(NULL) % 10000000;
  
  // Search for the executable in env var, plugin location (we hope) and subdir
  if (getenv("SHRMEMFRAME_DIR"))
    searchDirs[0] = getenv("SHRMEMFRAME_DIR");
  if (GetModuleFileName(thisModule, &dllPath[0], MAX_PATH)) {
    searchDirs[1] = dllPath;
    ind = (int)searchDirs[1].find_last_of("/\\");
    if (ind > 0) {
      searchDirs[1].resize(ind);
      searchDirs[2] = searchDirs[1] + "\\Framealign";
    } else
      searchDirs[1] = "";
  }
    
  searchDirs[3] = "C:\\ProgramData\\Gatan\\Plugins";
  searchDirs[4] = searchDirs[1] + "\\Framealign";
  for (ind = 0; ind < 4; ind++) {
    if (searchDirs[ind].length()) {
      shrMemPath = searchDirs[ind] + "\\shrmemframe.exe";
      if (!_stat(shrMemPath.c_str(), &statbuf))
        break;
    }
  }
  if (ind >= 4) {
    ERR_PRINT(PRN_FMT "Could not find shrmemframe.exe for frame alignment\n" PREFIX);
    return 1;
  } 
  sprintf(buffer, "\" %d", mServerID);
  shrMemPath = "\"" + shrMemPath + buffer; 
  argString = strdup(shrMemPath.c_str());
  if (!argString) {
    ERR_PRINT(PRN_FMT "Failure to duplicate command string for starting shrmemframe\n" 
      PREFIX);
    free(argString);
    return 1;
  } 
  
  // Try to start the process
  ZeroMemory( &sinfo, sizeof(sinfo) );
  sinfo.cb = sizeof(sinfo);
  ZeroMemory( &pinfo, sizeof(pinfo) );
  if (!CreateProcess(NULL, argString, NULL, NULL, FALSE,
                     CREATE_NO_WINDOW/*DETACHED_PROCESS*/,  // It lingered with DETACHED
                     NULL, NULL,  // Use parent's directory
                     &sinfo, &pinfo)) {
    sprintf(buffer, "Could not run %s with CreateProcess (error %d).\n", 
            shrMemPath.c_str(), GetLastError());
    ERR_PRINT(PRN_FMT buffer PREFIX);
    free(argString);
    return 1;
  } 
  free(argString);
  mProcessHandle = pinfo.hProcess;
  mThreadHandle = pinfo.hThread;

  sprintf(buffer, "%s%d", ACTION_BASE_NAME, mServerID);
  mActionSignal =  CreateEventA(NULL, FALSE, FALSE, (LPCSTR)buffer);
  sprintf(buffer, "%s%d", DONE_BASE_NAME, mServerID);
  mDoneSignal =  CreateEventA(NULL, FALSE, FALSE, (LPCSTR)buffer);
  if (!mActionSignal) {
    ERR_PRINT(PRN_FMT "Could not create one of the events for signaling with "
              "shrmemframe\n" PREFIX);
    Disconnect(-1);
    return 1;
  }

  startTime = GetTickCount();
  ind = 0;
  while (TickInterval(startTime) < initialTimeout) {
    if (WaitForSingleObject(mDoneSignal, initTryInterval) == WAIT_OBJECT_0) {
      ind = 1;
      break;
    }
    if (CheckIfProcessDied("shrmemframe died while waiting for a initial signal\n"))
      return 3;
  }
  if (!ind) {
    ERR_PRINT(PRN_FMT "Timeout waiting for initial signal from shrmemframe\n" PREFIX); 
    Disconnect(-1);
    return 3;
  }

  mMapFile = OpenFileMapping(FILE_MAP_ALL_ACCESS,   // read/write access
                             FALSE,                 // do not inherit the name
                             MAPPING_NAME);         // name of mapping object 
  if (mMapFile == NULL) {
    Disconnect(-1);
    sprintf(buffer, "Could not open file mapping object (error %d).\n", GetLastError());
    ERR_PRINT(PRN_FMT buffer PREFIX);
    return 1;
  } 

  mMappedBuf = (char *)MapViewOfFile(mMapFile,   // handle to map object
    FILE_MAP_ALL_ACCESS, // read/write permission
    0, 0, mShMemSize);
  if (!mMappedBuf) {
    CloseHandle(mMapFile);
    mMapFile = NULL;
    Disconnect(-1);
    sprintf(buffer, "Could not map view of file mapping object (error %d).\n", 
      GetLastError());
    ERR_PRINT(PRN_FMT buffer PREFIX);
    return 2;
  }
  mParamBuf = mMappedBuf + 2 * INT_SIZE;

  double wallStart = wallTime();
  if (SetCodeAndWaitForReply(SHMFR_VERSION, 2000))
    return 1;
  printf("Version obtained in %.3f ms\n", 1000. * (wallTime() - wallStart));
  memcpy(&version, mParamBuf, INT_SIZE);
  if (version != SHRMEMFRAME_VERSION) {
    sprintf(buffer, "The plugin is expecting version %d and shrmemframe is version %d; "
      "%s needs to be upgraded to match\n", SHRMEMFRAME_VERSION, version, 
      version < SHRMEMFRAME_VERSION ? "shrmemframe.exe" : "the SEMCCD plugin");
    ERR_PRINT(PRN_FMT buffer PREFIX);
    Disconnect();
    return 1;
  }
  return 0;
}

// Close up everything; wait defaults to 2000, if <= 0, it skips sending the exit code
// and terminates the process if it is < 0
void ShrMemClient::Disconnect(int wait)
{
  int funcCode = SHMFR_EXIT;
  if (mMappedBuf) {

    // All error returns from this call back into this function
    if (wait > 0 && SetCodeAndWaitForReply(SHMFR_EXIT, wait))
      return;
  } else
    wait = -1;

  ClearProcessVars(wait < 0);
  UnmapViewOfFile(mMappedBuf);
  CloseHandle(mMapFile);
  mMapFile = NULL;
  mMappedBuf = NULL;
  if (mActionSignal)
    CloseHandle(mActionSignal);
  mActionSignal = NULL;
  if (mDoneSignal)
    CloseHandle(mDoneSignal);
  mDoneSignal = NULL;
}

// The initialize call
int ShrMemClient::initialize(int binSum, int binAlign, float trimFrac, int numAllVsAll,
                 int cumAlignAtEnd, int useHybrid, int deferSum, int groupSize,
                 int nx, int ny, float padFrac, float taperFrac, int antiFiltType,
                 float radius1, float *radius2, float sigma1, float *sigma2, 
                 int numFilters, int maxShift, float kFactor, float maxMaxWeight,
                 int summingMode, int expectedZ, int makeUnwgtSum, int gpuFlags,
                 int debug)
{
  InitializeParams params;
  if (ConnectIfNeeded())
    return -1;
  params.binSum = binSum;
  params.binAlign = binAlign;
  params.trimFrac = trimFrac;
  params.numAllVsAll = numAllVsAll;
  params.cumAlignAtEnd = cumAlignAtEnd;
  params.useHybrid = useHybrid;
  params.deferSum = deferSum;
  params.groupSize = groupSize;
  params.nx = nx;
  params.ny = ny;
  params.padFrac = padFrac;
  params.taperFrac = taperFrac;
  params.antiFiltType = antiFiltType;
  params.radius1 = radius1;
  memcpy(&params.radius2[0], radius2, numFilters * sizeof(float));
  params.sigma1 = sigma1;
  memcpy(&params.sigma2[0], sigma2, numFilters * sizeof(float));
  params.numFilters = numFilters;
  params.maxShift = maxShift;
  params.kFactor = kFactor;
  params.maxMaxWeight = maxMaxWeight;
  params.summingMode = summingMode;
  params.expectedZ = expectedZ;
  params.makeUnwgtSum = makeUnwgtSum;
  params.gpuFlags = gpuFlags;
  params.debug = debug;
  memcpy(mParamBuf, &params, sizeof(InitializeParams));
  if (SetCodeAndWaitForReply(SHMFR_INITIALIZE, 5000))
    return 1;
  memcpy(&params, mParamBuf, sizeof(InitializeParams));
  RelayMessages(params.messages);
  mNumFrames = 0;
  mNx = nx;
  mNy = ny;
  return params.retVal;
}

// Test for GPU available
int ShrMemClient::gpuAvailable(int nGPU, float *memory, int debug)
{
  GpuAvailableParams params;
  if (ConnectIfNeeded())
    return -1;
  params.nGPU = nGPU;
  params.debug = debug;
  memcpy(mParamBuf, &params, sizeof(GpuAvailableParams));
  if (SetCodeAndWaitForReply(SHMFR_GPU_AVAIL, 5000))
    return 0;
  memcpy(&params, mParamBuf, sizeof(GpuAvailableParams));
  RelayMessages(params.messages);
  *memory = params.memory;
  return params.retVal;
}

// Pass the next frame
int ShrMemClient::nextFrame(void *frame, int type, float *gainRef, int nxGain, int nyGain,
                float truncLimit,
                std::string &defectStr, int camSizeX, int camSizeY, int defBin,
                float shiftX, float shiftY)
{
  NextFrameParams params;
  int dataSize, dum;
  unsigned int frameOffset, gainRefOffset = 0, defectsOffset = 0, cumOffset;
  if (CheckIfProcessDied(mNumFrames ? "shrmemframe died after processing last frame\n" : 
    "shrmemframe died after initializing frame alignment\n"))
    return 1;
  dataSizeForMode(type, &dataSize, &dum);
  frameOffset = 6 * INT_SIZE + sizeof(NextFrameParams);
  cumOffset = frameOffset + dataSize * mNx * mNy;
  if (gainRef) {
    gainRefOffset = cumOffset;
    cumOffset += sizeof(float) * nxGain * nyGain;
  }
  if (defectStr.length() && camSizeX) {
    defectsOffset = cumOffset;
    cumOffset += (unsigned int)defectStr.length() + 4;
  }
  if (cumOffset >= mShMemSize) {
    ERR_PRINT(PRN_FMT "Data are too large for the shared memory allocated\n" PREFIX);
    return 1;
  }
  params.frameOffset = frameOffset;
  params.type = type;
  params.gainRefOffset = gainRefOffset;
  params.nxGain = nxGain;
  params.nyGain = nyGain;
  params.truncLimit = truncLimit;
  params.defectsOffset = defectsOffset;
  params.camSizeX = camSizeX;
  params.camSizeY = camSizeY;
  params.defBin = defBin;
  params.shiftX = shiftX;
  params.shiftY = shiftY;

  // Copy params and frame
  memcpy(mParamBuf, &params, sizeof(NextFrameParams));
  if (frame != mMappedBuf + frameOffset)
    memcpy(mMappedBuf + frameOffset, frame, dataSize * mNx * mNy);

  // Copy gain ref and defects the first time
  if (!mNumFrames) {
    if (gainRef)
      memcpy(mMappedBuf + gainRefOffset, gainRef, sizeof(float) * nxGain * nyGain);
    if (defectsOffset) {
      strcpy(mMappedBuf + defectsOffset, defectStr.c_str());
    }
  }
  if (SetCodeAndWaitForReply(SHMFR_NEXTFRAME, 20000))
    return 1;
  memcpy(&params, mParamBuf, sizeof(NextFrameParams));
  RelayMessages(params.messages);
  mNumFrames++;
  return params.retVal;
}

// return address of part of buffer for frame so it can be put there directly
void *ShrMemClient::getFrameBuffer(void)
{
  return mMappedBuf + 6 * INT_SIZE + sizeof(NextFrameParams);
} 


// Finish the alignment and sum
int ShrMemClient::finishAlignAndSum(int nxOut, int nyOut, int numFilt, int maxRing, 
    float refineRadius2, float refineSigma2, float iterCrit, int groupRefine, 
    int doSpline, float **alisum, float *xShifts, float *yShifts, float *rawXshifts,
    float *rawYshifts, float *ringCorrs, float deltaR, int &bestFilt, float *smoothDist,
    float *rawDist, float *resMean, float *resSD, float *meanResMax, float *maxResMax,
    float *meanRawMax, float *maxRawMax)
{
  int ind, shiftBytes = mNumFrames * sizeof(float);
  unsigned int shiftOffset;
  FinishAlignAndSumParams params;
  if (CheckIfProcessDied("shrmemframe died after last frame was sent\n"))
    return 1;
  params.alisumOffset = 6 * INT_SIZE + sizeof(FinishAlignAndSumParams);
  params.refineRadius2 = refineRadius2;
  params.refineSigma2 = refineSigma2;
  params.iterCrit = iterCrit;
  params.groupRefine = groupRefine;
  params.doSpline = doSpline;
  params.deltaR = deltaR;
  shiftOffset = params.alisumOffset + nxOut * nyOut * sizeof(float);;
  params.xShiftsOffset = shiftOffset;
  params.yShiftsOffset = shiftOffset + shiftBytes;
  params.rawXshiftsOffset = shiftOffset + 2 * shiftBytes;
  params.rawYshiftsOffset = shiftOffset + 3 * shiftBytes;
  memcpy(mParamBuf, &params, sizeof(FinishAlignAndSumParams));
  if (SetCodeAndWaitForReply(SHMFR_FINISH, 120000))
    return 1;
  memcpy(&params, mParamBuf, sizeof(FinishAlignAndSumParams));
  RelayMessages(params.messages);
  *alisum = (float *)(mMappedBuf + params.alisumOffset);
  bestFilt = params.bestFilt;
  memcpy(xShifts, mMappedBuf + params.xShiftsOffset, shiftBytes);
  memcpy(yShifts, mMappedBuf + params.yShiftsOffset, shiftBytes);
  memcpy(rawXshifts, mMappedBuf + params.rawXshiftsOffset, shiftBytes);
  memcpy(rawYshifts, mMappedBuf + params.rawYshiftsOffset, shiftBytes);
  if (ringCorrs)
    for (ind = 0; ind < maxRing; ind++)
      ringCorrs[ind] = params.ringCorrs[ind];
  for (ind = 0; ind < numFilt; ind++) {
    smoothDist[ind] = params.smoothDist[ind];
    rawDist[ind] = params.rawDist[ind];
    resMean[ind] = params.resMean[ind];
    resSD[ind] = params.resSD[ind];
    meanResMax[ind] = params.meanResMax[ind];
    maxResMax[ind] = params.maxResMax[ind];
    meanRawMax[ind] = params.meanRawMax[ind];
    maxRawMax[ind] = params.maxRawMax[ind];
  }
  return params.retVal;
}

// Cleanup call
void ShrMemClient::cleanup(void)
{
  if (CheckIfProcessDied("shrmemframe died before call to clean up\n"))
    return;
  SetCodeAndWaitForReply(SHMFR_CLEANUP, 5000);
}

// Analyze the FRC crossings in the curve
void ShrMemClient::analyzeFRCcrossings(int maxRing, float *ringCorrs, float frcDelta, 
  float &halfCross, float &quartCross, float &eighthCross, float &halfNyq)
{
  AnalyzeFRCcrossingsParams params;
  params.frcDelta = frcDelta;
  for (int ind = 0; ind < maxRing; ind++)
    params.ringCorrs[ind] = ringCorrs[ind];
  if (CheckIfProcessDied("shrmemframe died before call to analyze FRC crossings\n"))
    return;
  memcpy(mParamBuf, &params, sizeof(AnalyzeFRCcrossingsParams));
  SetCodeAndWaitForReply(SHMFR_FRC_CROSS, 2000);
  memcpy(&params, mParamBuf, sizeof(AnalyzeFRCcrossingsParams));
  halfCross = params.halfCross;
  quartCross = params.quartCross;
  eighthCross = params.eighthCross;
  halfNyq = params.halfNyq;
}

// Checks if the process is still alive.  If it is not, closes process handles, unmaps the
// buffer and returns 1
int ShrMemClient::CheckIfProcessDied(const char *message)
{
  if (!mProcessHandle)
    return 1;
  if (WaitForSingleObject(mProcessHandle, 2) == WAIT_TIMEOUT)
    return 0;
  ERR_PRINT(PRN_FMT message PREFIX);
  Disconnect(0);
  return 1;
}

// Closes up variables related to process, terminating it if called for
void ShrMemClient::ClearProcessVars(bool terminate)
{
  if (mProcessHandle) {
    if (terminate)
      TerminateProcess(mProcessHandle, 1);
    CloseHandle(mProcessHandle);
    CloseHandle(mThreadHandle);
  }
  mProcessHandle = NULL;
}

// Set the action code and ID and set the signal, then wait for the done event
int ShrMemClient::SetCodeAndWaitForReply(int funcCode, int timeout)
{
  int reply, idVal;
  int replyWait = B3DMAX(timeout / 10, B3DMIN(50, timeout));
  replyWait = B3DMIN(replyWait, B3DMIN(1000, timeout));
  double startTime = GetTickCount();
  char buffer[100];
  sprintf(buffer, "shrmemframe died while waiting for a reply to code %d\n", funcCode);
  
  memcpy(mMappedBuf + INT_SIZE, &funcCode, INT_SIZE);
  memcpy(mMappedBuf, &mServerID, INT_SIZE);
  SetEvent(mActionSignal);

  // Check periodically for whether the program is still alive, and make sure
  // the codes are consistent if there is a signal
  while (TickInterval(startTime) < timeout) {
    if (WaitForSingleObject(mDoneSignal, replyWait) == WAIT_OBJECT_0) {
      memcpy(&idVal, mMappedBuf, INT_SIZE);
      memcpy(&reply, mMappedBuf + INT_SIZE, INT_SIZE);
      if (reply == SHMFR_DONE && idVal == mServerID)
        return 0;
      ERR_PRINT(PRN_FMT "shrmemframe gave the wrong respond on last function call; "
                "killing it\n" PREFIX); 
      Disconnect(-1);
      return 1;
    }
    if (CheckIfProcessDied(buffer))
      return 1;
  }
  ERR_PRINT(PRN_FMT "timed out on last function call to shrmemframe; killing it\n"
            PREFIX); 
  Disconnect(-1);
  return 1;
}

// Pass the messages on through whatever mechanism prints them
void ShrMemClient::RelayMessages(const char *messageBuf)
{
  if (messageBuf[0])
    DEB_PRINT(PRN_FMT messageBuf PREFIX);
}
