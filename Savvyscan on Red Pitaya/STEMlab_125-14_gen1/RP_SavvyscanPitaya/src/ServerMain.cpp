
// Linux main for Savvyscan camera server (port of SERVER-SEMCamServer.cpp)
// Removes WinSock/COM/MessageBox and uses POSIX + std::thread

#include "ServerCamera.h"
#include "BaseServer.h"
#include "Logging.h"

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <csignal>
#include <string>

// Subclass
class CServer : public CBaseServer {
public:
    CServer() { mHandshakeCode = JS_ChunkHandshake; }
    ~CServer() { };
    static int ProcessCommand(int sockInd, int numBytes);
    int StartSocket(int &osError);
    static void CleanupChanArrays(int indStart, int indEnd);
};

// Function prototypes required by BaseServer
int DoFinishStartup(int sockInd);
int DoProcessCommand(int sockInd, int numExpected);
int GetDebugVal();

// Statics
static CJediCamera sJediCam;
static CServer     sServer;
static bool        sExiting = false;
static FILE       *sFPdebug = nullptr;
static int         sDebug   = 0;
static short      *sChanArrays[MAX_SCAN_CHANNELS];
static int         sNumChanAcquired = 0;
static int         sNextChanIndex   = 0;
static long        sArrSize         = 0;

// --- Logging impls ---
void EitherToLog(const char *prefix, const char *message, bool saveErr)
{
    (void)saveErr; // not used in Linux port; kept for signature compatibility
    std::fprintf(stdout, "%s%s\n", prefix ? prefix : "", message ? message : "");
    std::fflush(stdout);
    if (sFPdebug) { std::fprintf(sFPdebug, "%s%s\n", prefix ? prefix : "", message ? message : ""); std::fflush(sFPdebug); }
}

void ErrorToLog(const char *message)
{
    EitherToLog("Error: ", message);
    // copy to error buffer for JS_GetErrorString
    std::snprintf(CServer::mErrorBuf, MESS_ERR_BUFF_SIZE, "%s", message ? message : "");
}

void DebugToLog(const char *message)
{
    if (sDebug) EitherToLog("", message ? message : "");
}

// --- Required globals ---
int DoFinishStartup(int sockInd)
{
    (void)sockInd;
    EitherToLog("This is a camera server for SerialEM\n",
                "Closing this program will end SerialEM's ability to connect to the Savvyscan camera\n");
    return 0;
}

int DoProcessCommand(int sockInd, int numExpected)
{
    return CServer::ProcessCommand(sockInd, numExpected);
}

int GetDebugVal()
{
    return sDebug;
}

// --- Signal handler ---
static void CtrlHandler(int signo)
{
    (void)signo;
    sExiting = true;
    CServer::ShutdownSocket();
    sJediCam.UninitializeCameras();
    sJediCam.UninitializeObjects();
    DebugToLog("Termination signal received\n");
}

// --- Main ---
int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    // Install basic signal handlers for Ctrl-C/termination
    std::signal(SIGINT,  CtrlHandler);
    std::signal(SIGTERM, CtrlHandler);

    int osError = 0;
    if (sJediCam.InitializeObjects()) {
        ErrorToLog("Failed to connect properly to the Savvyscan camera");
        return 1;
    }

    int ret = sServer.StartSocket(osError);

    sJediCam.UninitializeCameras();
    sJediCam.UninitializeObjects();
    return ret;
}

// --- Implementation ---
int CServer::StartSocket(int &osError)
{
    osError = 0;
    int sockInd = 0;

    // Default port; allow override via environment variable
    mPort[0] = 48901;
    if (const char *portStr = std::getenv("SAVVYSCAN_SEM_SERVER_PORT")) {
        int iPort = std::atoi(portStr);
        if (iPort <= 1024 || iPort > 65535) return 11; // invalid port
        mPort[0] = static_cast<unsigned short>(iPort);
    }

    if (const char *dbgStr = std::getenv("SAVVYSCAN_SEM_SERVER_DEBUG")) {
        sDebug = std::atoi(dbgStr);
    }
    if (sDebug & 4) {
        sFPdebug = std::fopen("/tmp/SavvySemServerDebug.txt", "w");
    }
    sDebug = sDebug % 4;
    sJediCam.SetDebugMode(sDebug);

    mInitialized[0] = true;
	//std::snprintf(CServer::mMessageBuf[0], MESS_ERR_BUFF_SIZE,"Starting SocketProc on port %d\n", (int)mPort[0]);
	//DebugToLog(CServer::mMessageBuf[0]);

    // Run socket loop in this thread (like original StartSocket did)
    SocketProc(sockInd);
    osError = mLastWSAerror[0];  // POSIX errno mapped in BaseServer
    return mStartupError[0];
}

void CServer::CleanupChanArrays(int indStart, int indPastEnd)
{
    for (int ind = indStart; ind < indPastEnd; ++ind) {
        delete [] sChanArrays[ind];
        sChanArrays[ind] = nullptr;
    }
}

// Table of functions
static ArgDescriptor sFuncTable[] = {
    {JS_InitializeCamera, 1, 0, 0, 0, 0, 0, 0, "InitializeCamera"},
    {JS_UninitializeCameras, 0, 0, 0, 0, 0, 0, 0, "UninitializeCameras"},
    {JS_GetNumberOfCameras, 0, 0, 0, 1, 0, 0, 0, "GetNumberOfCameras"},
    {JS_SelectCamera, 1, 0, 0, 0, 0, 0, 0, "SelectCamera"},
    {JS_SetDebugMode, 1, 0, 0, 0, 0, 0, 0, "SetDebugMode"},
    {JS_AcquireCCDImage, 8, 0, 1, 4, 0, 0, 0, "AcquireCCDImage"},
    {JS_AcquireSTEMImage, 9, 0, 2, 4, 0, 0, 1, "AcquireSTEMImage"},
    {JS_GetSTEMProperties, 2, 0, 0, 2, 0, 5, 0, "GetSTEMProperties"},
    {JS_StopContinuous, 0, 0, 0, 0, 0, 0, 0, "StopContinuous"},
    {JS_GetErrorString, 0, 0, 0, 1, 0, 0, 2, "GetErrorString"},
    {JS_GetNextChannel, 0, 0, 0, 4, 0, 0, 0, "GetNextChannel"},
    {-1, 0,0,0,0,0,0,0, nullptr}
};

int CServer::ProcessCommand(int sockInd, int numBytes)
{
    (void)sockInd;
    int ind, numChan, retSend = 0;
    int32_t *longArgs = mLongArgs[0];
    short *imArray = nullptr;

    if (PrepareCommand(0, numBytes, sFuncTable, "YOU PROBABLY NEED TO UPGRADE THIS SERVER TO MATCH THE CURRENT SERIALEM VERSION", ind))
        return 1;

    // If there are still channels waiting and this is not a GetNextChannel, clean up
    if (sNumChanAcquired > 0 && mLongArgs[0][0] != JS_GetNextChannel) {
        CleanupChanArrays(0, sNumChanAcquired);
        sNumChanAcquired = 0;
    }

    try {
        switch (mLongArgs[0][0]) {
        case JS_SelectCamera:
            SendArgsBack(sJediCam.SelectCamera(longArgs[1]));
            break;
        case JS_InitializeCamera:
            SendArgsBack(sJediCam.InitializeCamera(longArgs[1]));
            break;
        case JS_UninitializeCameras:
            SendArgsBack(sJediCam.UninitializeCameras());
            break;
        case JS_GetNumberOfCameras:
            mLongArgs[0][1] = sJediCam.GetNumberOfCameras();
            SendArgsBack(mLongArgs[0][1] <= 0 ? 1 : 0);
            break;
 
		case JS_AcquireCCDImage:{
			// Allocate image buffer based on requested size (original behavior)
			short *imArray = new short[mLongArgs[0][1]];

			// Use 64-bit 'long' temporaries for the camera API, then copy back to int32_t
			long arrSz = static_cast<long>(mLongArgs[0][1]);
			long szX   = static_cast<long>(mLongArgs[0][2]);
			long szY   = static_cast<long>(mLongArgs[0][3]);

			int ret = sJediCam.AcquireCCDImage(
				imArray,
				&arrSz, &szX, &szY,
				mLongArgs[0][4],           // top
				mLongArgs[0][5],           // left
				mLongArgs[0][6],           // binning
				mLongArgs[0][7],           // processing
				mDoubleArgs[0][0],         // exposure (double)
				mLongArgs[0][8]            // flags
			);

			// Copy results back to 32-bit protocol args
			mLongArgs[0][1] = static_cast<int32_t>(arrSz);
			mLongArgs[0][2] = static_cast<int32_t>(szX);
			mLongArgs[0][3] = static_cast<int32_t>(szY);

			// Send image (2 bytes per pixel for 'short')
			SendImageBack(ret, imArray, 2);
			break;	}
        case JS_AcquireSTEMImage: {
            numChan = longArgs[9];
            int32_t width  = longArgs[2];
            int32_t height = longArgs[3];
            sArrSize = static_cast<long>(width) * static_cast<long>(height);
            for (ind = 0; ind < numChan; ++ind) sChanArrays[ind] = new short[sArrSize];
            sNumChanAcquired = 0;

            long lx = width, ly = height, arrSz = sArrSize;
            retSend = sJediCam.AcquireSTEMImage(sChanArrays, &arrSz, &lx, &ly,
                                                longArgs[4], longArgs[5], longArgs[6],
                                                mDoubleArgs[0][0], mDoubleArgs[0][1],
                                                longArgs[7], numChan, (int*)mLongArray[0], longArgs[8], &sNumChanAcquired);
            sArrSize = arrSz;
            longArgs[1] = static_cast<int32_t>(sArrSize);
            longArgs[2] = static_cast<int32_t>(lx) + (static_cast<int32_t>(ly) << 16);
            longArgs[3] = static_cast<int32_t>(sNumChanAcquired);

            SendImageBack(retSend, sChanArrays[0], 2);
			usleep(100000);
            sChanArrays[0] = nullptr;
            CleanupChanArrays(sNumChanAcquired, numChan);
            sNextChanIndex = 1;
            if (sNumChanAcquired == 1) sNumChanAcquired = 0;
			usleep(100000);
            break; }
        case JS_GetNextChannel:
            if (!sNumChanAcquired || sNextChanIndex >= sNumChanAcquired || sNextChanIndex < 1 ||
                sNextChanIndex >= MAX_SCAN_CHANNELS || !sChanArrays[sNextChanIndex]) {
                ErrorToLog("GetNextChannel called inappropriately: no available channels\n");
                retSend = 1;
                break;
            }
            longArgs[1] = static_cast<int32_t>(sArrSize);
            longArgs[2] = static_cast<int32_t>(sNumChanAcquired);
            longArgs[3] = 0;
            SendImageBack(retSend, sChanArrays[sNextChanIndex], 2);
            sChanArrays[sNextChanIndex++] = nullptr;
            if (sNextChanIndex >= sNumChanAcquired) sNumChanAcquired = 0;
            return 0;
        case JS_GetSTEMProperties:
            SendArgsBack(sJediCam.GetSTEMProperties(longArgs[1], longArgs[2],
                                                    &mDoubleArgs[0][0], &mDoubleArgs[0][1], &mDoubleArgs[0][2], &mDoubleArgs[0][3],
                                                    &mDoubleArgs[0][4], (int*)(&longArgs[1]), (int*)(&longArgs[2])));
            break;
        case JS_StopContinuous:
            SendArgsBack(sJediCam.StopContinuous());
            break;
        case JS_SetDebugMode:
            sDebug = mLongArgs[0][1];
            SendArgsBack(sJediCam.SetDebugMode(mLongArgs[0][1]));
            break;
        case JS_GetErrorString:
            mLongArgs[0][1] = (int)((std::strlen(mErrorBuf) + 4) / 4);
            mLongArray[0] = (int32_t*)mErrorBuf;
            SendArgsBack(0);
            break;
        default:
            SendArgsBack(-1);
            break;
        }
    } catch (...) {
        SendArgsBack(-2);
    }
    return 0;
}
