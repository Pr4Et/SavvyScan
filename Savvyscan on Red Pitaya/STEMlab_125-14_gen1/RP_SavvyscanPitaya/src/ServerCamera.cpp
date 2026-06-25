
// Linux-adapted implementation of ServerCamera.cpp (no Win32/COM)
// - Uses std::chrono for timing
// - Uses std::mutex/std::atomic for synchronization
// - Provides a STEM simulation path that loads ./test.bmp (8-bit grayscale)
//   and returns it as 16-bit short image data.
//
// Build: GNU C++20 (requires ServerCamera.h, Logging.h)

#include "ServerCamera.h"
#include "feeder.h"   
#include "Logging.h"

#include <cmath>
#include <chrono>
#include <thread>
#include <mutex>
#include <atomic>
#include <cstring>
#include <cstdio>
#include <cstdarg>
#include <vector>
#include <string>
#include <new>
#include <algorithm>

#define MSG_STR_SIZE        320
#define IMAGE_MUTEX_WAIT    1000     // not used with std::mutex; kept for parity
#define CONTINUOUS_TIMEOUT  2000
#define SETTING_INDEX       9

// --------------------- Module / static state (Linux) ---------------------
static bool   sObjInitialized    = false;
static bool   sCamObjExists      = false;
static bool   sScanObjExists     = false;
static bool   sCamInitialized    = false;
static int    sScanInitialized   = 0;
static char   sMessageStr[MSG_STR_SIZE];
static char   sErrorStr[MSG_STR_SIZE];
static int    sCurrentCam        = 0;

static short* sLiveBuffer        = nullptr;
static int    sLiveWidth         = 0;
static int    sLiveHeight        = 0;
static int    sArrSize           = 0;
static int    sDivideBy2         = 0;
static int    sLiveError         = 0;
static std::atomic<int> sFrameReady{0};
static std::atomic<int> sWaitingForFrame{0};
static std::mutex sImageMutex;
static int sDebug = 0;

static double flyback_us = 160.0;    // retained from original notes
static const double PI    = 3.141592653589793;

// --------------------- Time helpers ---------------------
static inline uint64_t get_millis()
{
    using namespace std::chrono;
    return duration_cast<milliseconds>(steady_clock::now().time_since_epoch()).count();
}

double TickInterval(double start_ms)
{
    double now = static_cast<double>(get_millis());
    double interval = now - start_ms;
    if (interval < 0) interval = 0;
    return interval;
}

// --------------------- Logging helpers ---------------------
void ErrorFmt(char *fmt, ...)
{
    va_list args; va_start(args, fmt);
    std::vsnprintf(sErrorStr, MSG_STR_SIZE, fmt, args);
    va_end(args);
    ErrorToLog(sErrorStr);
}

void DebugFmt(char *fmt, ...)
{
    va_list args; va_start(args, fmt);
    std::vsnprintf(sMessageStr, MSG_STR_SIZE, fmt, args);
    va_end(args);
    DebugToLog(sMessageStr);
}

// --------------------- CJediCamera implementation ---------------------
CJediCamera::CJediCamera(void)
: m_pCameraComInterfaceClass(nullptr)
, m_pCameraSettingCom(nullptr)
, m_pFrameImageCom(nullptr)
, m_pEventReceiver(nullptr)
, mNumCameras(0)
, mUserSettingIndex(-1)
{
    // No Win32 CreateMutex needed; std::mutex is static above.
}

CJediCamera::~CJediCamera(void)
{
}

// Initialize non-hardware objects (pure software state).
int CJediCamera::InitializeObjects()
{
    if (sObjInitialized)
        return 0;

    // In the Windows version, COM _Module.Init would be here. Not needed on Linux.
    // If your project needs any Linux-side init, do it here.

    // Advertise two logical cameras like original (0: CCD, 1: STEM)
    mNumCameras     = 2;
    sObjInitialized = true;
    sCamObjExists   = true;
    sScanObjExists  = true;

    // You can set reasonable defaults for STEM properties, if desired:
    //mScanSizeX = 2048; mScanSizeY = 2048;
    //mMinScanSizeX = 8; mMinScanSizeY = 8;
    //mRotationStep = 0.1;
    //mMinFrameIntegration = 1;
    //mMaxFrameIntegration = 16;
    //mMinPixelTime = 1.0;     // microseconds
    //mMaxPixelTime = 100000;  // microseconds
    //mPixelTimeInc = 1.0;

    return 0;
}

int CJediCamera::UninitializeObjects()
{
    sObjInitialized = false;
    // No COM _Module.Term on Linux
    return 0;
}

int CJediCamera::GetNumberOfCameras()
{
    return mNumCameras;
}

int CJediCamera::InitializeCamera(int camNum)
{
    if (!sObjInitialized || camNum < 0 || camNum > MAX_SCAN_CHANNELS)
        return 1;

    if (!camNum) {
        // Index 0 is CCD
        if (sCamInitialized)
            return 0;
        if (!sCamObjExists)
            return 1;
        sScanInitialized = camNum;   // retain original semantics
    } else {
        // STEM
        sScanInitialized = 1;
    }
	monitor_init();//Start to monitor serialem_in.txt, run in feeder.cpp 
    return 0;
}

int CJediCamera::GetSTEMProperties(int useSizeX, int useSizeY, double *minPixel,
                                   double *maxPixel, double *pixelInc,
                                   double *rotationInc, double *ddum, int *maxIntegration,
                                   int *idum)
{
    if (!sScanInitialized)
        return 1;

    if (useSizeX > 0 && useSizeY > 0) {
        mScanSizeX = useSizeX;
        mScanSizeY = useSizeY;
    }
    *minPixel       = mMinPixelTime;
    *maxPixel       = mMaxPixelTime;
    *pixelInc       = mPixelTimeInc;
    *rotationInc    = mRotationStep;
    *ddum           = 0.0;
    *maxIntegration = mMaxFrameIntegration;
    *idum           = 0;
    return 0;
}

int CJediCamera::UninitializeCameras(void)
{
    if (sCamInitialized) {
        sCamInitialized = false;
        if (m_pCameraComInterfaceClass)
            m_pCameraComInterfaceClass->raw_TemUnInitialize();
    }
    if (sScanInitialized) {
        sScanInitialized = 0;
    }
	//monitor_shutdown();//related to monitor serialem_in.txt
    return 0;
}

int CJediCamera::SelectCamera(int camNum)
{
    // Original code: return error if camNum / 2 (i.e., camNum >= 2)
    if (camNum / 2)
        return 1;
    sCurrentCam = camNum;
    return 0;
}

int CJediCamera::SetDebugMode(int value)
{
    sDebug = value;
    return 0;
}

int CJediCamera::AcquireCCDImage(short int *array, long *arrSize, long *sizeX,
                                 long *sizeY, int top, int left, int binning,
                                 int processing, double exposure, int flags)
{
    // Not implemented in this Linux port (server currently focused on STEM)
    (void)array; (void)arrSize; (void)sizeX; (void)sizeY;
    (void)top; (void)left; (void)binning; (void)processing; (void)exposure; (void)flags;
    return 1;
}

// Helper: Load an 8-bit grayscale BMP from disk into vector<unsigned char> and get width/height.
// Returns 0 on success, nonzero on error.
static int LoadBMP8uGray(const char* filename, int &width, int &height, std::vector<unsigned char> &data)
{
    FILE* f = std::fopen(filename, "rb");
    if (!f) return 1;

    // BMP header (minimal)
    unsigned char header[54];
    if (std::fread(header, 1, 54, f) != 54) { std::fclose(f); return 2; }
    if (header[0] != 'B' || header[1] != 'M') { std::fclose(f); return 3; }

    // Extract fields
    const uint32_t dataOffset = *reinterpret_cast<uint32_t*>(&header[10]);
    const uint32_t dibSize    = *reinterpret_cast<uint32_t*>(&header[14]);
    (void)dibSize;
    width  = *reinterpret_cast<int32_t*>(&header[18]);
    height = *reinterpret_cast<int32_t*>(&header[22]);
    const uint16_t planes = *reinterpret_cast<uint16_t*>(&header[26]);
    const uint16_t bpp    = *reinterpret_cast<uint16_t*>(&header[28]);
    const uint32_t comp   = *reinterpret_cast<uint32_t*>(&header[30]);
    const uint32_t imgSz  = *reinterpret_cast<uint32_t*>(&header[34]);

    if (planes != 1 || comp != 0) { std::fclose(f); return 4; }
    if (bpp != 8) { std::fclose(f); return 5; } // Expect 8-bit grayscale (with palette)

    // Palette size is typically 1024 bytes for 8-bit BMP
    // We will seek to dataOffset and read the image payload
    if (std::fseek(f, dataOffset, SEEK_SET) != 0) { std::fclose(f); return 6; }

    // BMP rows are padded to 4-byte boundary
    const int rowStride = ((width + 3) / 4) * 4;
    const size_t payload = (size_t)rowStride * (size_t)height;
    std::vector<unsigned char> raw(payload);
    if (std::fread(raw.data(), 1, payload, f) != payload) { std::fclose(f); return 7; }
    std::fclose(f);

    // Convert from padded rows and (likely) bottom-up to a contiguous, top-down width*height buffer
    data.resize((size_t)width * (size_t)height);
    for (int y = 0; y < height; ++y) {
        const int srcY = height - 1 - y; // BMP is bottom-up by default
        std::memcpy(&data[(size_t)y * (size_t)width], &raw[(size_t)srcY * (size_t)rowStride], (size_t)width);
    }
    return 0;
}

// Acquire STEM image given binned coordinates/size and binning.
// Simulation mode: reads ./test.bmp and returns it as a 16-bit image (short per pixel).

int CJediCamera::AcquireSTEMImage(
    short **arrays, long *arrSize, long *sizeX, long *sizeY,
    int top, int left, int binning, double exposure,
    double rotation, int integration, int numChan,
    int *channels, int flags, int *numAcquired)
{
	//checkworking() 
    // Basic arg handling / early-outs
    int width  = static_cast<int>(*sizeX);
    int height = static_cast<int>(*sizeY);
    sDivideBy2 = flags & 1;
    const bool doContinuous = (flags & 2) != 0;

    if (!integration) integration = 1;
    if (!sScanInitialized) {
        ErrorToLog("You need to restart SerialEM to reinitialize the Savvyscan access");
        return 1;
    }

    //if (!doContinuous && sLiveBuffer)
    //    StopContinuous();

    //if (doContinuous && sLiveBuffer) {
    //    *numAcquired = 1;
    //    return ReturnLiveFrame(arrays[0], false, sizeX, sizeY);
    //}

    // ----- Scan geometry + feeder parameters (software mode) -----
    const int  ScanMode           = 1;
    const int  scTop              = top;
    const int  scLeft             = left;
    int        scXsize            = static_cast<int>(width * (1.0 + 0.414));  // include flyback strip
    int        scYsize            = height;
    int        oversampling       = 1;                                        // 1 point per pixel for now
    int        fsize              = height * scXsize + height * static_cast<int>(flyback_us / std::max(exposure, 1e-6));
    double     sampling_exposure  = exposure * (double)width * height / (double)fsize; //concerned with the sampling rate                                 // AWG clock ~ pixel time
    const double AspectRatio      = 1.0;

    // Ensure output buffer exists for channel 0
    const int pixels = width * height;
    if (!arrays || !arrays[0]) {
        // Actually done with sChanArrays in ServerMain: arrays[0] = new (std::nothrow) short[pixels]; 
        ErrorToLog("arrays[0] is null: caller must allocate the image buffer");
		std::fprintf(stderr,"arrays[0] is null: caller must allocate the image buffer");
        return 1;
    }

    // ----- Prepare scan, (optionally) start it, and acquire image -----
    prepare_AWG(binning, scXsize, scYsize, AspectRatio, rotation, oversampling,
                exposure, ScanMode, flyback_us, width, height,
                sampling_exposure, fsize);

    // In software-only builds this is a stub that returns 0
	int number_of_Arina_triggers=width * height;
    const int fail = activate_scan(sampling_exposure, number_of_Arina_triggers);
    if (fail) {
        // RestoreUserSettings(); // if needed in your flow
        return fail;
    }

    // Produce the image into arrays[0]
    if (!acquire_image(width, height, arrays[0])) {
        ErrorToLog("acquire_image failed");
        return 1;
    }

    const int nToAcquire = std::min(numChan, 1);
    *numAcquired = nToAcquire;
    *sizeX = width;
    *sizeY = height;
    *arrSize = static_cast<long>(pixels);
    return 0;
}


void CJediCamera::RestoreUserSettings()
{
    // No-op in Linux port
}

int CJediCamera::ReturnLiveFrame(short *array, bool firstTime, long *sizeX, long *sizeY)
{
    int error = 0;
    double startTime = static_cast<double>(get_millis());
    sWaitingForFrame = 1;

    for (;;) {
        if (sFrameReady || (!firstTime && TickInterval(startTime) > CONTINUOUS_TIMEOUT)) {
            std::lock_guard<std::mutex> lk(sImageMutex);
            if (sLiveError) {
                StopContinuous();
                error = sLiveError;
            }
            if (!error && !sLiveBuffer)
                error = 1;

            if (!error) {
                std::memcpy(array, sLiveBuffer, (size_t)sLiveWidth * (size_t)sLiveHeight * sizeof(short));
                *sizeX = sLiveWidth;
                *sizeY = sLiveHeight;
            }
            sFrameReady = 0;
            sWaitingForFrame = 0;
            return error;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
}

int CJediCamera::StopContinuous()
{
    int retVal = 0;
    {
        std::lock_guard<std::mutex> lk(sImageMutex);
        delete [] sLiveBuffer;
        sLiveBuffer = nullptr;
    }
    if (sCurrentCam > 0)
        RestoreUserSettings();
    return retVal;
}
