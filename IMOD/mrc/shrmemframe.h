/*
 *  shrmemframe.h -- header file for shrmemframe
 */
#ifndef SHRMEMFRAME_H
#define SHRMEMFRAME_H

#include <string>
#include "framealign.h"
//#include "CorrectDefects.h"

#define SHRMEMFRAME_VERSION 101

// Set some sizes, and components for the shared memory size
#define MAX_RING_CORRS      510
#define MESS_BUF_SIZE       4096
#define MAX_FOR_RAW_BYTE    (5 * (11520 * 8184))
#define MAX_FOR_RAW_SHORT   (6 * (7420 * 7676))
#define MAX_STRUCT_SIZE     (4 * (40 + 8 * MAX_FILTERS + MAX_RING_CORRS) + MESS_BUF_SIZE)
#define ADD_FOR_DEFECTS     100000
#define ADD_FOR_MARGIN      (MAX_FOR_RAW_BYTE / 20)
#define ALL_ADDED_SPACE     (MAX_STRUCT_SIZE + ADD_FOR_DEFECTS + ADD_FOR_MARGIN);
#define MAPPING_NAME        "FRAMEALIGN_MAPPING"
#define ACTION_BASE_NAME    "FRAMEALIGN_ACTION_EVENT"
#define DONE_BASE_NAME      "FRAMEALIGN_DONE_EVENT"

enum {SHMFR_DONE = 0, SHMFR_INITIALIZE, SHMFR_NEXTFRAME, SHMFR_FINISH, SHMFR_FRC_CROSS,
      SHMFR_GPU_AVAIL, SHMFR_CLEANUP, SHMFR_VERSION, SHMFR_EXIT};

// The shared memory size will be:
//    max(MAX_FOR_RAW_BYTE, MAX_FOR_RAW_SHORT) + ALL_ADDED_SPACE

// Structures for each process contain the relevant variables, plust a return value
// at the start and possible message buffer at the end
struct InitializeParams {
  int retVal;
  int binSum;
  int binAlign;
  float trimFrac;
  int numAllVsAll;
  int cumAlignAtEnd;
  int useHybrid;
  int deferSum;
  int groupSize;
  int nx;
  int ny;
  float padFrac;
  float taperFrac;
  int antiFiltType;
  float radius1;
  float radius2[MAX_FILTERS];
  float sigma1;
  float sigma2[MAX_FILTERS];
  int numFilters;
  int maxShift;
  float kFactor;
  float maxMaxWeight;
  int summingMode;
  int expectedZ;
  int makeUnwgtSum;
  int gpuFlags;
  int debug;
  char messages[MESS_BUF_SIZE];
};

struct NextFrameParams {
  int retVal;
  unsigned int frameOffset;
  int type;
  unsigned int gainRefOffset;
  int nxGain;
  int nyGain;
  float truncLimit;
  unsigned int defectsOffset;
  int camSizeX;
  int camSizeY;
  int defBin;
  float shiftX;
  float shiftY;
  char messages[MESS_BUF_SIZE];
};

struct FinishAlignAndSumParams {
  int retVal;
  float refineRadius2;
  float refineSigma2;
  float iterCrit;
  int groupRefine;
  int doSpline;
  unsigned int alisumOffset;
  unsigned int xShiftsOffset;
  unsigned int yShiftsOffset;
  unsigned int rawXshiftsOffset;
  unsigned int rawYshiftsOffset;
  float ringCorrs[MAX_RING_CORRS];
  float deltaR;
  int bestFilt;
  float smoothDist[MAX_FILTERS + 1];
  float rawDist[MAX_FILTERS + 1];
  float resMean[MAX_FILTERS + 1];
  float resSD[MAX_FILTERS + 1];
  float meanResMax[MAX_FILTERS + 1];
  float maxResMax[MAX_FILTERS + 1];
  float meanRawMax[MAX_FILTERS + 1];
  float maxRawMax[MAX_FILTERS + 1];
  char messages[MESS_BUF_SIZE];
};

struct AnalyzeFRCcrossingsParams {
  float ringCorrs[MAX_RING_CORRS];
  float frcDelta;
  float halfCross;
  float quartCross;
  float eighthCross;
  float halfNyq;
};

struct GpuAvailableParams {
  int retVal;
  int nGPU;
  float memory;
  int debug;
  char messages[MESS_BUF_SIZE];
};

class ShrMemFrame
{
public:
  ShrMemFrame();
  
  void main( int argc, char *argv[]);
 private:
  void initialize();
  void nextFrame();
  void finishSum();
  void analyzeFRC();
  void gpuAvailable();
  void sendReply(int reply);
  void closeAndExit(int exitCode);
  
  FrameAlign mFA;            // Instance of framealign
  HANDLE mMapFile;           // handle to mapping
  char *mMappedBuf;          // The shared memory buffer
  unsigned int mShMemSize;   // Its size
  HANDLE mActionSignal;      // Handles for the signals
  HANDLE mDoneSignal;
  int mIdValue;              // Unique ID value for this run
  int mNx, mNy;              // Input image size from initialize
  int mNumFrames;            // Have to keep track of frames passed here too!
  CameraDefects mCamDefects; // Parsed, process defects ready for every call
  int mCamSizeX, mCamSizeY;  // And the processed sizes for call
  int mDefCorBinning;        // And the binning is a computed value too!
  std::string mDefectString; // Defects converted to string
};
#endif
