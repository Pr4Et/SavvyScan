// BaseServer.cpp - Base socket class that can be used in any socket server
//
// Copyright (C) 2013-2016 by the Regents of the University of
// Colorado.  See Copyright.txt for full notice of copyright and limitations.
//
// Author: David Mastronarde
//
//Converted to LINUX by Shahar Seifer & Microsoft co-pilot


// Linux port of BaseServer: replaces WinSock & Win32 with POSIX sockets & std::thread
#include "BaseServer.h"
#include "feeder.h"   //to use void monitor_pump_once()
#include "Logging.h"

#include <cassert>
#include <cerrno>
#include <cstring>
#include <cstdio>
#include <string>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <thread>

#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>

// Static member definitions
int CBaseServer::mNumLongSend[MAX_SOCK_CHAN];
int CBaseServer::mNumBoolSend[MAX_SOCK_CHAN];
int CBaseServer::mNumDblSend[MAX_SOCK_CHAN];
int CBaseServer::mNumLongRecv[MAX_SOCK_CHAN];
int CBaseServer::mNumBoolRecv[MAX_SOCK_CHAN];
int CBaseServer::mNumDblRecv[MAX_SOCK_CHAN];

int32_t *CBaseServer::mLongArray[MAX_SOCK_CHAN];
int32_t  CBaseServer::mLongArgs[MAX_SOCK_CHAN][MAX_LONG_ARGS];  // 4-byte words on wire
double   CBaseServer::mDoubleArgs[MAX_SOCK_CHAN][MAX_DBL_ARGS];
BOOL32   CBaseServer::mBoolArgs[MAX_SOCK_CHAN][MAX_BOOL_ARGS];  // 4-byte words on wire

char *CBaseServer::mArgsBuffer[MAX_SOCK_CHAN];
int  CBaseServer::mNumBytesSend[MAX_SOCK_CHAN];
int  CBaseServer::mArgBufSize[MAX_SOCK_CHAN];
BOOL32 CBaseServer::mSendLongArray[MAX_SOCK_CHAN];

bool CBaseServer::mInitialized[MAX_SOCK_CHAN];
unsigned short CBaseServer::mPort[MAX_SOCK_CHAN];
std::thread *CBaseServer::mSocketThread[MAX_SOCK_CHAN];

int  CBaseServer::mStartupError[MAX_SOCK_CHAN];
int  CBaseServer::mLastWSAerror[MAX_SOCK_CHAN];   // keep name for log compatibility
bool CBaseServer::mCloseForExit[MAX_SOCK_CHAN];

char CBaseServer::mMessageBuf[MAX_SOCK_CHAN][MESS_ERR_BUFF_SIZE];
char CBaseServer::mErrorBuf[MESS_ERR_BUFF_SIZE] = {0x00};

// Chunking/handshake (keep values compatible with original)
int CBaseServer::mChunkSize = 16810000;        // ~16 MB chunks
int CBaseServer::mHandshakeCode = 0;
int CBaseServer::mSuperChunkSize = 33620000;   // ~32 MB superchunks
SOCKET CBaseServer::mHListener[MAX_SOCK_CHAN];
SOCKET CBaseServer::mHClient[MAX_SOCK_CHAN];
bool CBaseServer::mProcessingCommand = false;

static int bytesAccum = 0;
static int expectedSize = 0;


CBaseServer::CBaseServer()
{
    for (int i = 0; i < MAX_SOCK_CHAN; i++) {
        mInitialized[i] = false;
        mSocketThread[i] = nullptr;
        mStartupError[i] = -1;
        mLastWSAerror[i] = 0;
        mCloseForExit[i] = false;
        mArgsBuffer[i] = nullptr;
        mArgBufSize[i] = 0;
        mMessageBuf[i][0] = 0x00;
        mHClient[i] = INVALID_SOCKET;
        mHListener[i] = INVALID_SOCKET;
    }
    // parity with old Sleep(1)
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
}

void CBaseServer::ShutdownSocket(int sockInd)
{
    if (!mInitialized[sockInd]) return;

    mCloseForExit[sockInd] = true;

    // Close listener so select returns
    if (mHListener[sockInd] != INVALID_SOCKET) {
        ::shutdown(mHListener[sockInd], SHUT_RDWR);
        ::close(mHListener[sockInd]);
        mHListener[sockInd] = INVALID_SOCKET;
    }
    // Close client
    CloseClient(sockInd);

    // Join thread (if we ever use a separate thread)
    if (mSocketThread[sockInd]) {
        if (mSocketThread[sockInd]->joinable()) mSocketThread[sockInd]->join();
        delete mSocketThread[sockInd];
        mSocketThread[sockInd] = nullptr;
    }
    Cleanup(sockInd);
}

// POSIX thread entry / main loop
void CBaseServer::SocketProc(int sockInd)
{
    SOCKET hListener;
    sockaddr_in sockaddr{};
    timeval tv{};
    int yes = 1;
    int numBytes, err, numExpected;
    fd_set readFds;

    mArgsBuffer[sockInd] = (char*)std::malloc(ARGS_BUFFER_CHUNK);
    if (!mArgsBuffer[sockInd]) {
        mStartupError[sockInd] = 8;
        return;
    }
    mArgBufSize[sockInd] = ARGS_BUFFER_CHUNK;

    // Listener socket
    hListener = ::socket(PF_INET, SOCK_STREAM, 0);
    if (hListener == INVALID_SOCKET) {
        mLastWSAerror[sockInd] = errno;
        mStartupError[sockInd] = 4;
        return;
    }

    if (::setsockopt(hListener, SOL_SOCKET, SO_REUSEADDR, (char*)&yes, sizeof(int))) {
        mLastWSAerror[sockInd] = errno;
        mStartupError[sockInd] = 5;
        return;
    }

    sockaddr.sin_family = AF_INET;
    sockaddr.sin_port   = htons(mPort[sockInd]);
    sockaddr.sin_addr.s_addr = INADDR_ANY;

    if (::bind(hListener, (struct sockaddr*)(&sockaddr), sizeof(sockaddr))) {
        mLastWSAerror[sockInd] = errno;
        mStartupError[sockInd] = 6;
        return;
    }

    tv.tv_sec  = 0;
    tv.tv_usec = 1000 * SELECT_TIMEOUT;

    if (::listen(hListener, 2)) {
        mLastWSAerror[sockInd] = errno;
        mStartupError[sockInd] = 7;
        return;
    }
    mHListener[sockInd] = hListener;

    // Finish startup hook
    if (DoFinishStartup(sockInd))
        return;

    mStartupError[sockInd] = 0;
    std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE,
                  "Listening for connections on socket %d port %d\n", mHListener[sockInd],
                  (int)mPort[sockInd]);
    DebugToLog(mMessageBuf[sockInd]);

    // Main loop
    for (;;) {
        FD_ZERO(&readFds);
        FD_SET(hListener, &readFds);
        int maxfd = hListener;
        if (mHClient[sockInd] != INVALID_SOCKET) {
            FD_SET(mHClient[sockInd], &readFds);
            maxfd = std::max(maxfd, mHClient[sockInd]);
        }
        err = ::select(maxfd + 1, &readFds, NULL, NULL, &tv);
        if (err < 0 || mCloseForExit[sockInd]) {
            DebugToLog("Closing socket\n");
            CloseClient(sockInd);
            if (hListener != INVALID_SOCKET) ::close(hListener);
            if (err < 0 && !mCloseForExit[sockInd]) {
                mLastWSAerror[sockInd] = errno;
                mStartupError[sockInd] = 7;
                std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE,
                              "POSIX Error %d on select command\n", mLastWSAerror[sockInd]);
                ErrorToLog(mMessageBuf[sockInd]);
            }
            DebugToLog("returning\n");
            return;
        }
		
		monitor_pump_once();  //Run periodically check if there is update from SerialEM script, the file serialEM_in.txt.(run in feeder.cpp).
		
		
        if (GetDebugVal() > 1) {
            std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE,
                          "Select ready: listener %d client %d\n",
                          FD_ISSET(hListener, &readFds),
                          (mHClient[sockInd] != INVALID_SOCKET && FD_ISSET(mHClient[sockInd], &readFds)) ? 1 : 0);
            DebugToLog(mMessageBuf[sockInd]);
        }

		// Client data available?
		if (mHClient[sockInd] != INVALID_SOCKET && FD_ISSET(mHClient[sockInd], &readFds)) {
			int ret = ::recv(mHClient[sockInd],
							 mArgsBuffer[sockInd] + bytesAccum,
							 mArgBufSize[sockInd] - bytesAccum,
							 0);
			if (ret <= 0) {
				ReportErrorAndClose(sockInd, ret, "recv");
				bytesAccum = 0;
				expectedSize = 0;
				continue;
			}
			bytesAccum += ret;
			// Step 1: do we know message size?
			if (expectedSize == 0 && bytesAccum >= (int)sizeof(int)) {
				std::memcpy(&expectedSize, mArgsBuffer[sockInd], sizeof(int));
				// resize buffer if needed
				if (expectedSize > mArgBufSize[sockInd]) {
					mArgBufSize[sockInd] =
						((expectedSize + ARGS_BUFFER_CHUNK - 1) / ARGS_BUFFER_CHUNK) * ARGS_BUFFER_CHUNK;

					mArgsBuffer[sockInd] =
						(char*)std::realloc(mArgsBuffer[sockInd], mArgBufSize[sockInd]);

					if (!mArgsBuffer[sockInd]) {
						ErrorToLog("Failed to reallocate buffer");
						return;
					}
				}
			}
			// Step 2: do we have full message?
			if (expectedSize > 0 && bytesAccum >= expectedSize) {
				if (GetDebugVal() > 1) {
					std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE,
								  "Got full message: %d bytes\n", expectedSize);
					DebugToLog(mMessageBuf[sockInd]);
				}
				mProcessingCommand = true;
				DoProcessCommand(sockInd, expectedSize);
				mProcessingCommand = false;
				// Reset for next message
				bytesAccum = 0;
				expectedSize = 0;
			}
		}


        // New connection?
        if (FD_ISSET(hListener, &readFds)) {
            CloseClient(sockInd);
			//if (mHClient[sockInd] != INVALID_SOCKET) {
			//	EitherToLog("", "Accepted connection from client program\n");
			//}

            mHClient[sockInd] = ::accept(hListener, NULL, NULL);
            if (mHClient[sockInd] == INVALID_SOCKET)
                ReportErrorAndClose(sockInd, SOCKET_ERROR, "accept connection from ready client\n");
            else
                EitherToLog("", "Accepted connection from client program\n");
        }
    }
}

void CBaseServer::CloseClient(int sockInd)
{
    if (mHClient[sockInd] == INVALID_SOCKET) return;
    std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE,
                  "Closing connection to client via socket %d\n", mHClient[sockInd]);
    if (!mCloseForExit[sockInd])
        EitherToLog("", mMessageBuf[sockInd]);
    ::shutdown(mHClient[sockInd], SHUT_RDWR);
    ::close(mHClient[sockInd]);
    mHClient[sockInd] = INVALID_SOCKET;
}

void CBaseServer::Cleanup(int sockInd)
{
    // No WSACleanup on Linux
    mInitialized[sockInd] = false;
}

int CBaseServer::FinishGettingBuffer(int sockInd, int numReceived, int numExpected)
{
    int numNew, ind;
    while (numReceived < numExpected) {
        ind = numReceived;
        if (numExpected > mArgBufSize[sockInd]) ind = 0; // throw away start if oversized
        numNew = ::recv(mHClient[sockInd], &mArgsBuffer[sockInd][ind], mArgBufSize[sockInd] - ind, 0);
        if (numNew <= 0) {
            ReportErrorAndClose(sockInd, numNew, "recv to get expected number of bytes\n");
            return 1;
        }
        numReceived += numNew;
    }
    return 0;
}

int CBaseServer::PrepareCommand(int sockInd, int numBytes, ArgDescriptor *funcTable,
                                const char *upgradeMess, int &ind)
{
    int funcCode, needed, needAdd = 0;
    if (numBytes < 8 || numBytes > mArgBufSize[sockInd]) { SendArgsBack(sockInd, numBytes < 8 ? -4 : -5); return 1; }

    std::memcpy(&funcCode, &mArgsBuffer[sockInd][sizeof(int)], sizeof(int32_t));

    ind = 0;
    while (funcTable[ind].funcCode >= 0 && funcTable[ind].funcCode != funcCode) ind++;
    if (funcTable[ind].funcCode < 0) {
        std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, "Function code not found: %d\n%s",
                      funcCode, upgradeMess);
        ErrorToLog(mMessageBuf[sockInd]);
        SendArgsBack(sockInd, -1);
        return 1;
    }

    mNumLongRecv[sockInd] = funcTable[ind].numLongRecv + 1; // +1 for function code
    mNumBoolRecv[sockInd] = funcTable[ind].numBoolRecv;
    mNumDblRecv[sockInd]  = funcTable[ind].numDblRecv;
    mNumLongSend[sockInd] = funcTable[ind].numLongSend + 1; // +1 for retval
    mNumBoolSend[sockInd] = funcTable[ind].numBoolSend;
    mNumDblSend[sockInd]  = funcTable[ind].numDblSend;
    mSendLongArray[sockInd] = (funcTable[ind].hasLongArray & 2) > 0;

    needed = sizeof(int) + mNumLongRecv[sockInd] * sizeof(int32_t) +
             mNumBoolRecv[sockInd] * sizeof(BOOL32) + mNumDblRecv[sockInd] * sizeof(double);

    if (needed > numBytes) {
        std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE,
            "Command %d %s not long enough: needed = %d (4 + %d x 4 + %d x 4 + %d x 8) numBytes = %d\n",
            funcCode, funcTable[ind].label, needed, mNumLongRecv[sockInd], mNumBoolRecv[sockInd], mNumDblRecv[sockInd], numBytes);
        ErrorToLog(mMessageBuf[sockInd]);
        SendArgsBack(sockInd, -4);
        return 1;
    }

    if (UnpackReceivedData(sockInd)) { SendArgsBack(sockInd, -5); return 1; }

    if (funcTable[ind].hasLongArray & 1)
        needAdd = mLongArgs[sockInd][mNumLongRecv[sockInd] - 1] * sizeof(int32_t);
    needed += needAdd;

    if (needed != numBytes) {
        std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE,
            "Command %d %s wrong length: needed = %d (4 + %d x 4 + %d x 4 + %d x 8 + %d) numBytes = %d\n",
            funcCode, funcTable[ind].label, needed, mNumLongRecv[sockInd], mNumBoolRecv[sockInd], mNumDblRecv[sockInd], needAdd, numBytes);
        ErrorToLog(mMessageBuf[sockInd]);
        SendArgsBack(sockInd, -6);
        return 1;
    }

    if (GetDebugVal() > 1) {
        std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, "Processing command %d %s\n",
                      funcCode, funcTable[ind].label);
        DebugToLog(mMessageBuf[sockInd]);
    }
    return 0;
}

int CBaseServer::SendBuffer(int sockInd, char *buffer, int numBytes)
{
    int numTotalSent = 0;
    int numToSend, numSent;
    if (GetDebugVal() > 1) {
        std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE,
                      "In SendBuffer, socket %d, sending %d bytes\n", mHClient[sockInd], numBytes);
        DebugToLog(mMessageBuf[sockInd]);
    }
    while (numTotalSent < numBytes) {
        numToSend = numBytes - numTotalSent;
        if (numToSend > mChunkSize) numToSend = mChunkSize;
        if (GetDebugVal() > 1) {
            std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE,
                          "Going to send %d bytes to socket %d\n", numToSend, mHClient[sockInd]);
            DebugToLog(mMessageBuf[sockInd]);
        }
        numSent = ::send(mHClient[sockInd], &buffer[numTotalSent], numToSend, 0);
		if (numSent < 0) {
			if (errno == EWOULDBLOCK || errno == EAGAIN) {
				usleep(1000); //1ms
				continue;
			}
			ReportErrorAndClose(sockInd, numSent, "send failed");
			return 1;
		}
		if (numSent == 0) {
			// socket buffer full → wait briefly and retry
			usleep(1000); //1ms
			continue;
		}

        numTotalSent += numSent;
    }
    return 0;
}

void CBaseServer::ReportErrorAndClose(int sockInd, int retval, const char *message)
{
    if (retval == SOCKET_ERROR) {
        mLastWSAerror[sockInd] = errno;
        std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, "POSIX Error %d on call to %s\n",
                      mLastWSAerror[sockInd], message);
        ErrorToLog(mMessageBuf[sockInd]);
    } else {
        std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, "retval %d on call to %s",
                      retval, message);
        ErrorToLog(mMessageBuf[sockInd]);
    }
    CloseClient(sockInd);
}

void CBaseServer::CloseOnExitOrSelectError(int sockInd, int err)
{
    DebugToLog("Closing socket\n");
    CloseClient(sockInd);
    if (mHListener[sockInd] != INVALID_SOCKET) ::close(mHListener[sockInd]);
    if (err < 0) {
        mLastWSAerror[sockInd] = errno;
        mStartupError[sockInd] = 7;
        std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, "POSIX Error %d on select command\n",
                      mLastWSAerror[sockInd]);
        ErrorToLog(mMessageBuf[sockInd]);
    }
}

//not used anymore, note that numBytes = ::recv(mHClient ... is modifed in another place but not here.
int CBaseServer::ListenForHandshake(int sockInd, int superChunk)
{
    timeval tv{};
    int numBytes, err, numExpected, command;
    fd_set readFds;

    tv.tv_sec = 0;
    tv.tv_usec = superChunk / 5; // 5 MB/sec pacing (same heuristic)

    FD_ZERO(&readFds);
    FD_SET(mHClient[sockInd], &readFds);
    err = ::select(mHClient[sockInd] + 1, &readFds, NULL, NULL, &tv);
    if (err < 0 || mCloseForExit[sockInd]) { CloseOnExitOrSelectError(sockInd, err); return mStartupError[sockInd]; }

    // timeout: client fails
    if (!err) { ReportErrorAndClose(sockInd, 0, "timeout on handshake from client"); return 1; }

    numBytes = ::recv(mHClient[sockInd], mArgsBuffer[sockInd], mArgBufSize[sockInd], 0);
    std::memcpy(&numExpected, &mArgsBuffer[sockInd][0], sizeof(int));
    std::memcpy(&command,     &mArgsBuffer[sockInd][4], sizeof(int));

    if (command != mHandshakeCode || numExpected != 8 || numBytes != 8) {
        ReportErrorAndClose(sockInd, numBytes, "recv handshake from ready client");
        return 1;
    }
    return 0;
}

int CBaseServer::SendArgsBack(int sockInd, int retval)
{
    mLongArgs[sockInd][0] = retval;
    if (retval < 0) {
        mNumLongSend[sockInd] = 1;
        mNumBoolSend[sockInd] = 0;
        mNumDblSend[sockInd]  = 0;
        mSendLongArray[sockInd] = false;
    }
    if (PackDataToSend(sockInd)) {
        ErrorToLog("DATA BUFFER NOT BIG ENOUGH TO SEND REPLY TO SERIALEM");
        SendArgsBack(sockInd, -3);
        return 1;
    }
    return SendBuffer(sockInd, mArgsBuffer[sockInd], mNumBytesSend[sockInd]);
}

void CBaseServer::SendImageBack(int sockInd, int retval, short *imArray, int bytesPerPixel)
{
    int numChunks, chunkSize, numToSend, numLeft, err, imSize, totalSent = 0;
    imSize = mLongArgs[sockInd][1] * bytesPerPixel;
    numChunks = (imSize + mSuperChunkSize - 1) / mSuperChunkSize;
    mLongArgs[sockInd][4] = numChunks;
    err = SendArgsBack(sockInd, retval);
    std::snprintf(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE,
                  "retval = %d, err sending args %d, sending image %d in %d chunks\n",
                  retval, err, imSize, numChunks);
    DebugToLog(mMessageBuf[sockInd]);

    if (!err && !retval) {
        numLeft = imSize;
        chunkSize = (imSize + numChunks - 1) / numChunks;
        while (totalSent < imSize) {
            numToSend = chunkSize;
            if (chunkSize > imSize - totalSent) numToSend = imSize - totalSent;
            if (SendBuffer(sockInd, (char*)imArray + totalSent, numToSend)) break;
            totalSent += numToSend;
            //if (totalSent < imSize && ListenForHandshake(sockInd, numToSend)) break;
        }
    }
    delete [] imArray;
}

int CBaseServer::UnpackReceivedData(int sockInd)
{
    int numBytes, numUnpacked = sizeof(int);
    if (mNumLongRecv[sockInd] > MAX_LONG_ARGS || mNumBoolRecv[sockInd] > MAX_BOOL_ARGS || mNumDblRecv[sockInd] > MAX_DBL_ARGS)
        return 1;

    numBytes = mNumLongRecv[sockInd] * sizeof(int32_t);
    std::memcpy(mLongArgs[sockInd], &mArgsBuffer[sockInd][numUnpacked], numBytes);
    numUnpacked += numBytes;

    numBytes = mNumBoolRecv[sockInd] * sizeof(BOOL32);
    if (numBytes) std::memcpy(mBoolArgs[sockInd], &mArgsBuffer[sockInd][numUnpacked], numBytes);
    numUnpacked += numBytes;

    numBytes = mNumDblRecv[sockInd] * sizeof(double);
    if (numBytes) std::memcpy(mDoubleArgs[sockInd], &mArgsBuffer[sockInd][numUnpacked], numBytes);
    numUnpacked += numBytes;

    mLongArray[sockInd] = (int32_t*)(&mArgsBuffer[sockInd][numUnpacked]);
    return 0;
}

#define ADD_TO_BUFFER(a) \
    if (numAdd + mNumBytesSend[sockInd] > mArgBufSize[sockInd]) return 1; \
    std::memcpy(&mArgsBuffer[sockInd][mNumBytesSend[sockInd]], a, numAdd); \
    mNumBytesSend[sockInd] += numAdd;

int CBaseServer::PackDataToSend(int sockInd)
{
    int numAdd;
    mNumBytesSend[sockInd] = sizeof(int);

    if (mNumLongSend[sockInd]) { numAdd = mNumLongSend[sockInd] * sizeof(int32_t); ADD_TO_BUFFER(mLongArgs[sockInd]); }
    if (mNumBoolSend[sockInd]) { numAdd = mNumBoolSend[sockInd] * sizeof(BOOL32);  ADD_TO_BUFFER(&mBoolArgs[sockInd]); }
    if (mNumDblSend[sockInd])  { numAdd = mNumDblSend[sockInd]  * sizeof(double);  ADD_TO_BUFFER(&mDoubleArgs[sockInd]); }
    if (mSendLongArray[sockInd]) {
        numAdd = mLongArgs[sockInd][mNumLongSend[sockInd] - 1] * sizeof(int32_t);
        ADD_TO_BUFFER(mLongArray[sockInd]);
    }

    std::memcpy(&mArgsBuffer[sockInd][0], &mNumBytesSend[sockInd], sizeof(int));
    return 0;
}
