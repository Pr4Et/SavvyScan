/*
 * ShrMemClient.h - header file for ShrMemClient.cpp
 */
#pragma once
#include <string>

#ifdef TEST_SHRMEM
#include <Windows.h>
#endif

class ShrMemClient
{
public:
  ShrMemClient(void);
  ~ShrMemClient(void);
  int ConnectIfNeeded(void);
  void Disconnect(int wait = 2000);
  int initialize(int binSum, int binAlign, float trimFrac, int numAllVsAll,
                 int cumAlignAtEnd, int useHybrid, int deferSum, int groupSize,
                 int nx, int ny, float padFrac, float taperFrac, int antiFiltType,
                 float radius1, float *radius2, float sigma1, float *sigma2, 
                 int numFilters, int maxShift, float kFactor, float maxMaxWeight,
                 int summingMode, int expectedZ, int makeUnwgtSum, int gpuFlags,
                 int debug);
  int gpuAvailable(int nGPU, float *memory, int debug);
  int nextFrame(void *frame, int type, float *gainRef, int nxGain, int nyGain,
                float truncLimit,
                std::string &defectStr, int camSizeX, int camSizeY, int defBin,
                float shiftX, float shiftY);
  int finishAlignAndSum(int nxOut, int nyOut, int numFilt, int maxRing, 
    float refineRadius2, float refineSigma2, float iterCrit, int groupRefine, 
    int doSpline, float **alisum, float *xShifts, float *yShifts, float *rawXshifts,
    float *rawYshifts, float *ringCorrs, float deltaR, int &bestFilt, float *smoothDist,
    float *rawDist, float *resMean, float *resSD, float *meanResMax, float *maxResMax,
    float *meanRawMax, float *maxRawMax);
  void cleanup(void);
  void analyzeFRCcrossings(int maxRing, float *ringCorrs, float frcDelta, 
    float &halfCross, float &quartCross, float &eighthCross, float &halfNyq);
  void *getFrameBuffer(void);
private:
  HANDLE mMapFile;            // Handle to file mapping
  char *mMappedBuf;           // The mapped buffer
  char *mParamBuf;            // Convenience pointer to where the param starts
  unsigned int mShMemSize;    // Size of memory
  HANDLE mProcessHandle;      // Handle to process created
  HANDLE mThreadHandle;       // Handle to the "thread" on process creation
  HANDLE mActionSignal;       // Event for notifying of new action
  HANDLE mDoneSignal;         // Event for notifying action is done
  int mNumFrames;             // Keep track of number of frames so gain ref is copied once
  int mNx, mNy;               // Keep track of image size
  int mServerID;              // The unique ID for this run, send as an argument
  
public:
  int CheckIfProcessDied(const char *message);
  void ClearProcessVars(bool terminate);
  int SetCodeAndWaitForReply(int funcCode, int timeout);
  void RelayMessages(const char * messageBuf);
};
