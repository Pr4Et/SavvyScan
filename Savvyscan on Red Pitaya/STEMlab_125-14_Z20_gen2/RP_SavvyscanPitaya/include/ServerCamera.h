
#pragma once
#include <string>
#include <cstdint>

void ErrorToLog(const char *message);
void DebugToLog(const char *message);
void ErrorFmt(char *fmt, ...);
void DebugFmt(char *fmt, ...);

double TickInterval(double start);

// Protocol enums stay the same
enum { JS_GetErrorString = 1, JS_GetNumberOfCameras, JS_InitializeCamera,
       JS_UninitializeCameras, JS_SelectCamera, JS_AcquireCCDImage, JS_SetDebugMode,
       JS_StopContinuous, JS_ChunkHandshake, JS_AcquireSTEMImage, JS_GetNextChannel,
       JS_GetSTEMProperties };

#define MAX_CAM_BINNINGS 16
#define MAX_SCAN_CHANNELS 4

// Plain C++ placeholder classes (no COM)
class _FrameImageCom { public: explicit _FrameImageCom(bool f): flag(f){} bool flag; };
class _CameraComInterfaceClass { public: explicit _CameraComInterfaceClass(bool f): flag(f){} bool flag; void raw_TemUnInitialize(){} };
class _EventReceiver { public: explicit _EventReceiver(bool f): flag(f){} bool flag; };
class _CameraSettingCom { public: explicit _CameraSettingCom(bool f): flag(f){} bool flag; };

typedef struct FRAMESIZE { int Width; int Height; } FRAMESIZE;

class CJediCamera {
public:
    CJediCamera(void);
    ~CJediCamera(void);

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
                         int *channels, int flags, int *numAcquired);

    int SetDebugMode(int debug);
    int StopContinuous();
    int ReturnLiveFrame(short *array, bool firstTime, long *sizeX, long *sizeY);

    int GetSTEMProperties(int useSizeX, int useSizeY, double *minPixel, double *maxPixel,
                          double *pixelInc, double *rotationInc, double *ddum, int *maxIntegration, int *idum);

private:
    void RestoreUserSettings();

    int mNumCameras{};
    int mNumBinnings{};
    int mBinnings[MAX_CAM_BINNINGS]{};
    int mCamSizeX{}, mCamSizeY{};
    int mMinCamSizeX{}, mMinCamSizeY{};
    int mScanSizeX{}, mScanSizeY{};
    int mMinScanSizeX{}, mMinScanSizeY{};
    int mAIPnum{};
    double mRotationStep{};
    int mMinFrameIntegration{};
    int mMaxFrameIntegration{};
    int mMinExpTimeIndex{};
    int mMaxExpTimeIndex{};
    double mMinPixelTime{};
    double mMaxPixelTime{};
    double mPixelTimeInc{};

    std::string mDetectorNames[MAX_SCAN_CHANNELS + 1];
    int mAssignedIndex[MAX_SCAN_CHANNELS]{};
    bool mHasDetectorNames{};
    int mUserSettingIndex{};

    _EventReceiver *m_pEventReceiver{};
    _CameraComInterfaceClass *m_pCameraComInterfaceClass{};
    _CameraSettingCom *m_pCameraSettingCom{};
    _FrameImageCom *m_pFrameImageCom{};
};
