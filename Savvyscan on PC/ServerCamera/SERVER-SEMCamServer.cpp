// Based on JEOL-SEMCamServer.cpp from another project
//

#include "stdafxServer.h"
#include <stdio.h>
#include <winsock.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include "ServerCamera.h"
#include "BaseServer.h"

// Insist on at least Winsock v1.1
const unsigned int VERSION_MAJOR = 1;
const unsigned int VERSION_MINOR = 1;

// Subclass
class CServer : public CBaseServer
{
public:
  CServer() {mHandshakeCode = JS_ChunkHandshake;};
  ~CServer() {};
  static int ProcessCommand(int sockInd, int numBytes);
  int StartSocket(int &wsaError);
  static void CleanupChanArrays(int indStart, int indEnd);
};

static BOOL CtrlHandler( DWORD fdwCtrlType );

// Table of functions, with # of incoming long, bool, and double, # of outgoing
// long, bool and double, and whether there is a long array at the end, whose size is
// in the last long argument. Include that size in the number of incoming or outgoing
// long's.
static ArgDescriptor sFuncTable[] = {
  {JS_InitializeCamera,     1, 0, 0,   0, 0, 0,   0, "InitializeCamera"},
  {JS_UninitializeCameras,  0, 0, 0,   0, 0, 0,   0, "UninitializeCameras"},
  {JS_GetNumberOfCameras,   0, 0, 0,   1, 0, 0,   0, "GetNumberOfCameras"},
  {JS_SelectCamera,         1, 0, 0,   0, 0, 0,   0, "SelectCamera"},
  {JS_SetDebugMode,         1, 0, 0,   0, 0, 0,   0, "SetDebugMode"},
  {JS_AcquireCCDImage,      8, 0, 1,   4, 0, 0,   0, "AcquireCCDImage"},
  {JS_AcquireSTEMImage,     9, 0, 2,   4, 0, 0,   1, "AcquireSTEMImage"},
  {JS_GetSTEMProperties,    2, 0, 0,   2, 0, 5,   0, "GetSTEMProperties"},
  {JS_StopContinuous,       0, 0, 0,   0, 0, 0,   0, "StopContinuous"},
  {JS_GetErrorString,       0, 0, 0,   1, 0, 0,   2, "GetErrorString"},
  {JS_GetNextChannel,       0, 0, 0,   4, 0, 0,   0, "GetNextChannel"},
  {-1, 0,0,0,0,0,0,0, NULL}
};

// Static class instances
static CJediCamera sJediCam;
static CServer sServer;
static bool sComInitialized = false;
static short *sChanArrays[MAX_SCAN_CHANNELS];
static int sNumChanAcquired = 0;
static int sNextChanIndex = 0;
static long sArrSize;
static bool sExiting = false;
static FILE *sFPdebug = NULL;
static int sDebug = 0;


//int _tmain(int argc, _TCHAR* argv[])
int main(int argc, _TCHAR* argv[])
{
  int error = 0;
  if (CoInitializeEx(NULL, COINIT_MULTITHREADED) != S_OK)
    return 1;
  if (sJediCam.InitializeObjects()) {
    MessageBox(NULL, _T("Failed to connect properly to the Savvyscan camera"),
      _T("SERVER-SEMCamServer Error"), MB_OK | MB_TOPMOST);
  } else {
    if (sServer.StartSocket(error)) {
     if (!sExiting)
       MessageBox(NULL, _T("Failed to start the socket that allows SerialEM to connect to this server"), _T("SERVER-SEMCamServer Error"), MB_OK | MB_TOPMOST);
     else
       Sleep(5000);
    }
  }
  sJediCam.UninitializeCameras();
  sJediCam.UninitializeObjects();
  return error;
}
// Required globals and statics
int DoFinishStartup(int SockInd)
{
  EitherToLog("This is a camera server for SerialEM\n",
    "Closing this window will end SerialEM's ability to connect to the Savvyscan camera\n");
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

void ErrorToLog(const char *message)
{
  EitherToLog("Error: ", message);
  strncpy_s(CServer::mErrorBuf, MESS_ERR_BUFF_SIZE, message, _TRUNCATE);
}

void DebugToLog(const char *message)
{
  if (sDebug)
    EitherToLog("", message);
}

void EitherToLog(const char *prefix, const char *message, bool saveErr)
{
  printf("%s%s\n", prefix, message);
  fflush(stdout);
  if (sFPdebug) {
    fprintf(sFPdebug, "%s%s\n", prefix, message);
    fflush(sFPdebug);
  }
}

// On Windows 7 there is no choice about exiting, it will exit.
static BOOL CtrlHandler(DWORD dwCtrlType)
{
  sExiting = true;
  CServer::ShutdownSocket();
  sJediCam.UninitializeCameras();
  sJediCam.UninitializeObjects();
  DebugToLog("Termination event received\n");
  if (dwCtrlType == CTRL_C_EVENT || dwCtrlType == CTRL_CLOSE_EVENT) {
    MessageBox(NULL, _T("The server that allows SerialEM to connect\nto the camera is exiting.\n\nYou will need to restart both the server\nand SerialEM to reconnect properly"), _T("Exiting Camera Server to SerialEM"), MB_OK | MB_TOPMOST);
  }
  return TRUE;
}

// Call from main or a program component to test for port variable and optionally, 
// start a thread to manage sockets (which is needed for a program that needs to do other
// things)
int CServer::StartSocket(int &wsaError)
{
  //DWORD threadID;
  WSADATA WSData;
  char *portStr;
  int sockInd = 0;
  mPort[0] = 48901;

  // For getting port from environment variable
  _dupenv_s(&portStr, NULL, "SAVVYSCAN_SEM_SERVER_PORT");
  wsaError = 0;
  if (portStr) {
    int iPort = atoi(portStr);
    if (iPort <= 1024 || iPort > 65535)
      return 11;
    mPort[0] = (unsigned short)iPort;
    free(portStr);
  } /* else    // Uncomment if variable is required
    return 0; */

  // Set a debug level and optional log file with variable
  _dupenv_s(&portStr, NULL, "SAVVYSCAN_SEM_SERVER_DEBUG");
  if (portStr) {
    sDebug = atoi(portStr);
    free(portStr);
  }
  if (sDebug & 4) {
    if (fopen_s(&sFPdebug, "C:\\Program Files\\SerialEM\\SavvySemServerDebug.txt", "w"))
      sFPdebug = NULL;
  }
  sDebug = sDebug % 4;
  sJediCam.SetDebugMode(sDebug);

  if (!SetConsoleCtrlHandler((PHANDLER_ROUTINE)CtrlHandler, TRUE))
    ErrorToLog("Failed to register console control handler\n");

  if (WSAStartup(MAKEWORD(VERSION_MAJOR, VERSION_MINOR), &WSData)) {
    wsaError = WSAGetLastError();
    return 1;
  }
  mInitialized[0] = true;

  SocketProc(&sockInd);
  wsaError = mLastWSAerror[0];
  return mStartupError[0];
}

void CServer::CleanupChanArrays(int indStart, int indPastEnd)
{
  for (int ind = indStart; ind < indPastEnd; ind++) {
    delete sChanArrays[ind];
    sChanArrays[ind] = NULL;
  }
}

// Process a received message
int CServer::ProcessCommand(int sockInd, int numBytes)
{
  int ind, numChan, retSend = 0;
  long sizeX, sizeY;
  long *longArgs = mLongArgs[0];
  short *imArray;

  if (PrepareCommand(0, numBytes, sFuncTable, "YOU PROBABLY NEED TO UPGRADE THIS"
    " SERVER TO MATCH THE CURRENT SERIALEM VERSION", ind))
    return 1;

  // If there are still channels waiting and this is not a GetNextChannel, clean up
  if (sNumChanAcquired > 0 && mLongArgs[0][0] != JS_GetNextChannel) {
    CleanupChanArrays(0, sNumChanAcquired);
    sNumChanAcquired = 0;
  }

  // THE FUNCTION CALLS
  // A main point about the socket interface: the function code is received in
  // mLongArgs[0] and the return value is passed back in mLongArgs[0], so the first
  // integer argument passed in either direction is in mLongArgs[1]
  try {
    switch (mLongArgs[0][0]) {
    case JS_SelectCamera:
      SendArgsBack(sJediCam.SelectCamera(mLongArgs[0][1]));
      break;

    case JS_InitializeCamera:
      SendArgsBack(sJediCam.InitializeCamera(mLongArgs[0][1]));
      break;

    case JS_UninitializeCameras:
      SendArgsBack(sJediCam.UninitializeCameras());
      break;

    case JS_GetNumberOfCameras:
      mLongArgs[0][1] = sJediCam.GetNumberOfCameras();
      SendArgsBack(mLongArgs[0][1] <= 0 ? 1 : 0);
      break;

    case JS_AcquireCCDImage:
      imArray = new short[mLongArgs[0][1]];
      SendImageBack(sJediCam.AcquireCCDImage(imArray, &mLongArgs[0][1], &mLongArgs[0][2], 
        &mLongArgs[0][3], mLongArgs[0][4], mLongArgs[0][5], mLongArgs[0][6], 
        mLongArgs[0][7], mDoubleArgs[0][0], mLongArgs[0][8]), imArray, 2);
      break;

    
	case JS_AcquireSTEMImage:
      numChan = longArgs[9];
      sizeX = longArgs[2];
      sizeY = longArgs[3];
      sArrSize = longArgs[2] * longArgs[3];
      for (ind = 0; ind < numChan; ind++)
        sChanArrays[ind] = new short[sArrSize];
      sNumChanAcquired = 0;
      retSend = sJediCam.AcquireSTEMImage(sChanArrays, &sArrSize, &sizeX, &sizeY, 
        longArgs[4], longArgs[5], longArgs[6], mDoubleArgs[0][0], mDoubleArgs[0][1],
        longArgs[7], longArgs[9], (int *)mLongArray[0], longArgs[8], &sNumChanAcquired);
      longArgs[1] = sArrSize;
      longArgs[2] = sizeX + (sizeY << 16);
      longArgs[3] = sNumChanAcquired;
      SendImageBack(retSend, sChanArrays[0], 2);
      sChanArrays[0] = NULL;
      CleanupChanArrays(sNumChanAcquired, numChan);
      sNextChanIndex = 1;
      if (sNumChanAcquired == 1)
        sNumChanAcquired = 0;
      break;

    case JS_GetNextChannel:
      if (!sNumChanAcquired || sNextChanIndex >= sNumChanAcquired || sNextChanIndex < 1 ||
        sNextChanIndex >= MAX_SCAN_CHANNELS || !sChanArrays[sNextChanIndex]) {
          ErrorToLog("GetNextChannel called inappropriately: no available channels\n");
          retSend = 1;
          break;
        }
      longArgs[1] = sArrSize;
      longArgs[2] = sNumChanAcquired;
      longArgs[3] = 0;
      SendImageBack(retSend, sChanArrays[sNextChanIndex], 2);
      sChanArrays[sNextChanIndex++] = NULL;
      if (sNextChanIndex >= sNumChanAcquired)
        sNumChanAcquired = 0;
      return 0;

    case JS_GetSTEMProperties:
      SendArgsBack(sJediCam.GetSTEMProperties(longArgs[1], longArgs[2], 
        &mDoubleArgs[0][0], &mDoubleArgs[0][1], &mDoubleArgs[0][2], &mDoubleArgs[0][3],
        &mDoubleArgs[0][4], (int *)(&longArgs[1]), (int *)(&longArgs[2])));
      break;

    case JS_StopContinuous:
      SendArgsBack(sJediCam.StopContinuous());
      break;

    case JS_SetDebugMode:
      sDebug = mLongArgs[0][1];
      SendArgsBack(sJediCam.SetDebugMode(mLongArgs[0][1]));
      break;

    case JS_GetErrorString:
      mLongArgs[0][1] = (int)(strlen(mErrorBuf) + 4) / 4;
      mLongArray[0] = (long *)mErrorBuf;
      SendArgsBack(0);
      break;

    default:
      SendArgsBack(-1);  // Incorrect command
        break;
    }
  }
  catch (...) {
    SendArgsBack(-2);  // Memory allocation or other exception
  }
  return 0;
}
