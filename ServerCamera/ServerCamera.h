#pragma once
#include <atlbase.h>
extern CComModule _Module;
#include <atlcom.h>

#pragma warning (disable: 4278)
#import <mscorlib.tlb> raw_interfaces_only


// These are defined here because everyone reads it; they are provided by top level module
void ErrorToLog(const char *message);
void DebugToLog(const char *message);
void ErrorFmt(char *fmt, ...);
void DebugFmt(char *fmt, ...);
double TickInterval(double start);

// Needed on both sides of socket
enum {JS_GetErrorString = 1, JS_GetNumberOfCameras, JS_InitializeCamera, 
      JS_UninitializeCameras, JS_SelectCamera, JS_AcquireCCDImage, JS_SetDebugMode,
      JS_StopContinuous, JS_ChunkHandshake, JS_AcquireSTEMImage, JS_GetNextChannel, 
      JS_GetSTEMProperties};

#define MAX_CAM_BINNINGS 16
#define SINKID_COUNTEREVENTS 0
#define MAX_SCAN_CHANNELS  4


  class _FrameImageCom  //shahar
	{public:
	 bool flag;
	 _FrameImageCom (bool _flag) {flag=_flag;}
  };
  class _CameraComInterfaceClass  //shahar
	{public:
	 bool flag;
	 _CameraComInterfaceClass (bool _flag) {flag=_flag;}
	 void raw_TemUnInitialize(){};
  };
  class _EventReceiver  //shahar
	{public:
	 bool flag;
	 _EventReceiver (bool _flag) {flag=_flag;}
  };
  class _CameraSettingCom  //shahar
	{public:
	 bool flag;
	 _CameraSettingCom (bool _flag) {flag=_flag;}
  };

 //shahar
   typedef struct FRAMESIZE 
  {
	public:
		int Width;
		int Height;
  } FRAMESIZE;



class CJediCamera
{
 public:
  CJediCamera(void);
  ~CJediCamera(void);

 private:


 public:

  int InitializeObjects();
  int UninitializeObjects();
  int GetNumberOfCameras();
  int InitializeCamera(int camNum);
  int UninitializeCameras(void);
  int SelectCamera(int camNum);
  int AcquireCCDImage(short int *array, long *arrSize, long *sizeX, long *sizeY,
                             int top, int left, int binning, int processing,
                             double exposure, int flags);
  int AcquireSTEMImage(short **arrays, long *arrSize, long *sizeX, long *sizeY, int top, 
    int left, int binning, double exposure, double rotation, int integration, int numChan,
    int *channels,  int flags, int *numAcquired);
  int SetDebugMode(int debug);
  int StopContinuous();
  int ReturnLiveFrame(short *array, bool firstTime, long *sizeX, long *sizeY);
  int GetSTEMProperties(int useSizeX, int useSizeY, double *minPixel, double *maxPixel, 
    double *pixelInc, double *rotationInc, double *ddum, int *maxIntegration, int *idum);
  
   _EventReceiver *m_pEventReceiver;
  _CameraComInterfaceClass *m_pCameraComInterfaceClass;
  _CameraSettingCom *m_pCameraSettingCom;
  _FrameImageCom *m_pFrameImageCom;
	
   

private:
  void RestoreUserSettings();
  //void DumpSettings(_ScanSettingCom *pScanSettingCom, const char *when);

  int mNumCameras;
  int mNumBinnings;
  int mBinnings[MAX_CAM_BINNINGS];
  int mCamSizeX, mCamSizeY;
  int mMinCamSizeX, mMinCamSizeY;
  int mScanSizeX, mScanSizeY;
  int mMinScanSizeX, mMinScanSizeY;
  int mAIPnum;
  double mRotationStep;
  int mMinFrameIntegration;
  int mMaxFrameIntegration;
  int mMinExpTimeIndex;
  int mMaxExpTimeIndex;
  double mMinPixelTime;
  double mMaxPixelTime;
  double mPixelTimeInc;
  BSTR mDetectorNames[MAX_SCAN_CHANNELS + 1];
  int mAssignedIndex[MAX_SCAN_CHANNELS];
  bool mHasDetectorNames;
  int mUserSettingIndex;

 };
