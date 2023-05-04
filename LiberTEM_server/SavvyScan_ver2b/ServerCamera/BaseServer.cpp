// BaseServer.cpp - Base socket class that can be used in any socket server
//
// Copyright (C) 2013-2016 by the Regents of the University of
// Colorado.  See Copyright.txt for full notice of copyright and limitations.
//
// Author: David Mastronarde
//


#include "BaseServer.h"
#include <stdio.h>
#include <assert.h>
#include <string>
#include <fstream>
#include <zmq.h> 

int CBaseServer::mNumLongSend[MAX_SOCK_CHAN];
int CBaseServer::mNumBoolSend[MAX_SOCK_CHAN];
int CBaseServer::mNumDblSend[MAX_SOCK_CHAN];
int CBaseServer::mNumLongRecv[MAX_SOCK_CHAN];
int CBaseServer::mNumBoolRecv[MAX_SOCK_CHAN];
int CBaseServer::mNumDblRecv[MAX_SOCK_CHAN];
long *CBaseServer::mLongArray[MAX_SOCK_CHAN];
long CBaseServer::mLongArgs[MAX_SOCK_CHAN][MAX_LONG_ARGS];   // Max is 16
double CBaseServer::mDoubleArgs[MAX_SOCK_CHAN][MAX_DBL_ARGS];  // Max is 8
BOOL CBaseServer::mBoolArgs[MAX_SOCK_CHAN][MAX_BOOL_ARGS];   // Max is 8
char *CBaseServer::mArgsBuffer[MAX_SOCK_CHAN];
int CBaseServer::mNumBytesSend[MAX_SOCK_CHAN];
int CBaseServer::mArgBufSize[MAX_SOCK_CHAN];
BOOL CBaseServer::mSendLongArray[MAX_SOCK_CHAN];

bool CBaseServer::mInitialized[MAX_SOCK_CHAN];
unsigned short CBaseServer::mPort[MAX_SOCK_CHAN]; 
HANDLE CBaseServer::mHSocketThread[MAX_SOCK_CHAN];
int CBaseServer::mStartupError[MAX_SOCK_CHAN];
int CBaseServer::mLastWSAerror[MAX_SOCK_CHAN];
bool CBaseServer::mCloseForExit[MAX_SOCK_CHAN];
char CBaseServer::mMessageBuf[MAX_SOCK_CHAN][MESS_ERR_BUFF_SIZE];
char CBaseServer::mErrorBuf[MESS_ERR_BUFF_SIZE] = {0x00};
int CBaseServer::mChunkSize = 16810000;     // Tests indicated 16777216 was optimal size
int CBaseServer::mHandshakeCode;
//shahar, variables for zeromq communication
 void *context;
 void *responder;
 void *requester;
 int rc;
 void *context2;
 void *responder2;
 void *requester2;
 int rc2;
 void* context3;
 void* responder3;
 void* requester3;
 int rc3;
 int returnsize, zmq_msg_pos;
 char zmq_msg_buffer [250];
 char zmq_msg_buffer2 [250];
 char zmq_msg_buffer3[200];
 char zmq_sub_buffer [250];
 std::string zmq_msg;
 std::string zmq_cmd_msg;
//shahar, global vairbales from GUI to RecRepMulti.cpp
 char GUI_LiberTEMcommand[200]="OFF";
 char GUI_ArinaFileName[200];
 char GUI_foldername[239];
 char GUI_ImmediateDir[239]="A:\\SavvyscanData";
 char GUI_temp[239];
 char GUI_wtext[239];//was wchar_t
 char GUI_wtext2[239];
 int GUI_scanmode=1;
 double GUI_scanargument=0;
 double GUI_aspectratio=1;
 int GUI_InputAmpmV=1000;
 int GUI_OutputAmpmV=3000;
 int GUI_BiasOutputP = 0;
 int GUI_ch1_avg;
 int GUI_ch2_avg;
 int GUI_ch3_avg;
 int GUI_ch4_avg;
 int GUI_ch5_avg;
 bool GUI_ch_update_ready=false;
 bool GUI_tomography=false;
 int GUI_tomoIndex=0;
 int GUI_numberOfSlices=1;
 int GUI_chosenCH=7;
 int GUI_cella=61695;
 int GUI_LRcella = 61695;
 int GUI_defocus;
 int PYTHON_ask_defocus=0; //idle
 int GUI_Lothar=0;//idle, 1=active
 double GUI_tiltangle=0;
 bool GUI_saveMAT = false;
 bool GUI_saveMATKEY = false;
 bool GUI_align = false;
 int GUI_CalibrateBias = 0;
 int GUI_SerialEMAligned = 1;
 int GUI_LowResTimeS = 8;//maximum low res scan time in seconds
 bool GUI_ArinaON = false;
 bool save4python_latch = false;
 bool GUI_BF_chs_blocked = false;
 bool using_LiberTEM = false;
 bool LiberTEM_comp_connected = false;
 extern FILE** fpv;//connect to the RecRepMulti defined mrc files to close them


// 33554432 is exactly 4K x 2K x 2, failures occurred above twice this size
// This value fits the F416 4200x4200 buffer size
int CBaseServer::mSuperChunkSize = 33620000;  
SOCKET CBaseServer::mHListener[MAX_SOCK_CHAN];
SOCKET CBaseServer::mHClient[MAX_SOCK_CHAN]; 
bool CBaseServer::mProcessingCommand = false;

CBaseServer::CBaseServer()
{
  for (int i = 0; i < MAX_SOCK_CHAN; i++) {
    mInitialized[i] = false;
    mHSocketThread[i] = 0;
    mStartupError[i] = -1;
    mLastWSAerror[i] = 0;
    mCloseForExit[i] = false;
    mArgsBuffer[i] = NULL;
    mArgBufSize[i] = 0;
    mMessageBuf[i][0] = 0x00;
    mHClient[i] = INVALID_SOCKET;
  }

    
  //shahar: initiate communication with zeromq client
    //Looking for Shadow GUI interface but without wait (will catch up later)
    //start and test the contorl socket with GUI
	context = zmq_ctx_new ();
    responder = zmq_socket (context, ZMQ_REP);//I am the server (responder)
    rc = zmq_bind (responder, "tcp://*:5555"); //bind to client with non-fixed address
	assert (rc == 0);
    returnsize=zmq_recv (responder, zmq_msg_buffer, 250, ZMQ_DONTWAIT);
    if (returnsize>0)
		printf ("Received %s from GUI.\n",zmq_msg_buffer);
	else
		printf("No GUI interface\n");
    Sleep(1);          //  Do some 'work'
    zmq_send (responder, "World", 5, ZMQ_DONTWAIT);
	//printf(zmq_strerror (zmq_errno()) );
	//start and test the data channel socket with GUI
    //requester = zmq_socket (context, ZMQ_REQ);//I am the requester on socket 5556, but still the server
    //rc2 = zmq_bind (requester, "tcp://*:5556"); //bind to client with non-fixed address
	//assert (rc2 == 0);
	//zmq_send (requester, "test=1", 6, ZMQ_DONTWAIT);   //
	//returnsize=zmq_recv (responder, zmq_msg_buffer, 250, ZMQ_DONTWAIT); //must receive the reply

	//Initiate communiation with Lothar's python program via zeromq
	context2 = zmq_ctx_new ();
	requester2 = zmq_socket (context2, ZMQ_REQ);//I am the requester on socket 5557
	rc2=zmq_connect(requester2,"tcp://localhost:5557"); //bind to server to get answers (defocus number is the answer)
   

}

void CBaseServer::ShutdownSocket(int sockInd)
{
  DWORD code;
  if (!mInitialized[sockInd])
    return;
  if (mHSocketThread[sockInd]) {
    mCloseForExit[sockInd] = true;
    WaitForSingleObject(mHSocketThread[sockInd], 3 * SELECT_TIMEOUT);
    GetExitCodeThread(mHSocketThread[sockInd], &code);
    if (code == STILL_ACTIVE) {
      CloseClient(sockInd);
      closesocket(mHListener[sockInd]);

      // Terminating is cleaner with the GatanSocket case, could try suspend
      TerminateThread(mHSocketThread[sockInd], 1);
    }
    CloseHandle(mHSocketThread[sockInd]);
  } else {
    CloseClient(sockInd);
    closesocket(mHListener[sockInd]);
  }
  Cleanup(sockInd);
}


// The main socket thread routine, starts a listener, loops on getting connections and
// commands
DWORD WINAPI CBaseServer::SocketProc(LPVOID pParam)
{
  SOCKET hListener;
  SOCKADDR_IN sockaddr;
  struct timeval tv;
  BOOL yes = TRUE;
  int sockInd = *((int *)pParam);
  int numBytes, err, numExpected;
  fd_set readFds;      // file descriptor list for select()




  mArgsBuffer[sockInd] = (char *)malloc(ARGS_BUFFER_CHUNK);
  if (!mArgsBuffer[sockInd]) {
    mStartupError[sockInd] = 8;
    return mStartupError[sockInd];
  }
  mArgBufSize[sockInd] = ARGS_BUFFER_CHUNK;

  // Get the listener socket
  hListener = socket(PF_INET, SOCK_STREAM, 0);
  if (hListener == INVALID_SOCKET) {
    mLastWSAerror[sockInd] = WSAGetLastError();
    mStartupError[sockInd] = 4;
    return mStartupError[sockInd];
  }

  // Avoid "Address already in use" error message
  if (setsockopt(hListener, SOL_SOCKET, SO_REUSEADDR, (char *)&yes, sizeof(BOOL))) {
    mLastWSAerror[sockInd] = WSAGetLastError();
    mStartupError[sockInd] = 5;
    return mStartupError[sockInd];
  }

  // Get socket address for listener on the port
  sockaddr.sin_family = AF_INET;
  sockaddr.sin_port = htons(mPort[sockInd]);     // short, network byte order
  sockaddr.sin_addr.s_addr = INADDR_ANY;
  memset(sockaddr.sin_zero, '\0', sizeof(sockaddr.sin_zero));

  // Bind the listener socket to the port
  if (bind(hListener, (struct sockaddr *)(&sockaddr), sizeof(sockaddr))) {
    mLastWSAerror[sockInd] = WSAGetLastError();
    mStartupError[sockInd] = 6;
    return mStartupError[sockInd];
  }
  
  tv.tv_sec = 0;
  tv.tv_usec = 1000 * SELECT_TIMEOUT;

  // Listen
  if (listen(hListener, 2)) {
    mLastWSAerror[sockInd] = WSAGetLastError();
    mStartupError[sockInd] = 7;
    return mStartupError[sockInd];
  }
  mHListener[sockInd] = hListener;

  // Call out for special tasks
  if (DoFinishStartup(sockInd))
    return mStartupError[sockInd];
  
  // Set the value to indicate we are through all the startup
  mStartupError[sockInd] = 0;

  sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, 
    "Listening for connections on socket %d  port %d\n", mHListener[sockInd], 
    (int)mPort[sockInd]);
  DebugToLog(mMessageBuf[sockInd]);
  
 

  // Loop on listening for connections and accepting them or receiving commands
  for (;;) {
    FD_ZERO(&readFds);
    FD_SET(hListener, &readFds);
    if (mHClient[sockInd] != INVALID_SOCKET)
      FD_SET(mHClient[sockInd], &readFds);
    err = select(2, &readFds, NULL, NULL, &tv);
    if (err < 0 || mCloseForExit[sockInd]) {

      // Close up on error or signal from plugin
      DebugToLog("Closing socket\n");
      CloseClient(sockInd);
      closesocket(hListener);
      if (err < 0 && !mCloseForExit[sockInd]) {
        mLastWSAerror[sockInd] = WSAGetLastError();
        mStartupError[sockInd] = 7;
        sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, 
          "WSA Error %d on select command\n", mLastWSAerror[sockInd]);
        ErrorToLog(mMessageBuf[sockInd]);
      }
      DebugToLog("returning\n");
      return mStartupError[sockInd];
    }

	  
	//Shahar: transmit or receive message from Lothar's Python program
	if (PYTHON_ask_defocus==2)
	{
		returnsize=zmq_recv (requester2, zmq_msg_buffer2, 250, ZMQ_DONTWAIT | ZMQ_RCVMORE);
		if (returnsize>=0)
		{
			sscanf(zmq_msg_buffer2,"%d",& GUI_defocus); 
			PYTHON_ask_defocus=0;//should be =3
		}
	}
	if (PYTHON_ask_defocus==1)
	{
		zmq_send(requester2,"Defocus?",8, ZMQ_DONTWAIT);
		PYTHON_ask_defocus=2;
	}
	
	//Shahar: read messages from Shadow GUI (via zeroMQ)
    float tempfnumber;
    memset(zmq_msg_buffer,0,sizeof(zmq_msg_buffer));
	memset(zmq_sub_buffer,0,sizeof(zmq_sub_buffer));
	returnsize=zmq_recv (responder, zmq_msg_buffer, 250, ZMQ_DONTWAIT | ZMQ_RCVMORE);
	if (returnsize>=0)
	{
		zmq_msg=(std::string)zmq_msg_buffer;
		if (zmq_msg=="Hello")    
		{	zmq_send (responder, "World", 5, ZMQ_DONTWAIT);
			printf("GUI interface started");
		}
		else if (zmq_msg=="data?")
		{	
			if (GUI_ch_update_ready)
			{
                if (GUI_align && GUI_SerialEMAligned == 0)
                {
                    zmq_send(responder, "Use for alignment", 17, ZMQ_DONTWAIT);
                }
                if (GUI_align && GUI_SerialEMAligned == 1) //image for intenal use of SerialEM
                {
                    zmq_send(responder, "Use for @SerialEM", 17, ZMQ_DONTWAIT);
                }
                if (!GUI_align)
                { 
                    if (PYTHON_ask_defocus==3)
				    {
					    sprintf(zmq_sub_buffer,"ch1=%d,ch2=%d,ch3=%d,ch4=%d,ch5=%d,defocus=%d",GUI_ch1_avg,GUI_ch2_avg,GUI_ch3_avg,GUI_ch4_avg, GUI_ch5_avg, GUI_defocus);
					    PYTHON_ask_defocus=0;
                        save4python_latch = true;//second time onward the defocus value will be correct based on updated images
				    }
                    else
                    {
                        sprintf(zmq_sub_buffer, "ch1=%d,ch2=%d,ch3=%d,ch4=%d", GUI_ch1_avg, GUI_ch2_avg, GUI_ch3_avg, GUI_ch4_avg, GUI_ch5_avg);
                    }
                    zmq_send(responder, zmq_sub_buffer, sizeof(zmq_sub_buffer), ZMQ_DONTWAIT);
                }
                GUI_align = false;
                GUI_ch_update_ready=false;
			}
			else
				if (PYTHON_ask_defocus==3)
				{
					sprintf(zmq_sub_buffer,"defocus=%d",GUI_defocus);
					PYTHON_ask_defocus=0;
                    save4python_latch = true;//second time onward the defocus value will be correct based on updated images
				}
				else
					zmq_send (responder, "none", 4, 0);
		}
        else
		{			
			printf ("GUI: %s\n",zmq_msg_buffer);
			zmq_send (responder, "OK", 2, 0);//In REQ-REP mode I must reply to be able to send again
		}
		zmq_msg_pos=zmq_msg.find("=");
		if (zmq_msg_pos>=0)
		{
			zmq_cmd_msg=zmq_msg.substr(0,zmq_msg_pos); //the command in string variable: like "scanmode" in "scanmode=1"
			strcpy(zmq_sub_buffer,zmq_msg_buffer+zmq_msg_pos+1); //the argument in char array: like "1" in "scanmode=1"
			if (zmq_cmd_msg=="save") 
			{
				strcpy(GUI_foldername,zmq_sub_buffer);//save=foldername will copy existing files to the foldername
				//std::filesystem::copy; not in this c++ version
				//mbstowcs(GUI_wtext,GUI_foldername,strlen(GUI_foldername)+1);//GUI_wtext= GUI_foldername in wide char format, then convert to LPCSTR in win32
                bool foldercreated = CreateDirectory((LPCTSTR)GUI_foldername, NULL); //create if no such directory
                bool copied=false;
				for (int i=0; i<=8;i++) //copy all files from temporary directory to newly created directory
				{
					sprintf(GUI_wtext,"%s\\CH%d.mrc",GUI_ImmediateDir,i);
                    sprintf(GUI_wtext2,"%s\\CH%d.mrc",GUI_foldername,i);
                    if (i == 7 || (i == 8 && GUI_chosenCH == 8) || (i >= 0 && i <= 6 && !GUI_BF_chs_blocked) || (i == 8 && using_LiberTEM))
                    {
                        copied = CopyFile((LPCTSTR)GUI_wtext, (LPCTSTR)GUI_wtext2, true);//copy all files between directories in windows
                    }
                    //sprintf(GUI_wtext, "%s\\CH%d.dat", GUI_ImmediateDir, i);
                    //sprintf(GUI_wtext2, "%s\\CH%d.dat", GUI_foldername, i);
                    //copied = CopyFile((LPCTSTR)GUI_wtext, (LPCTSTR)GUI_wtext2, true);//copy all files between directories in windows
                    if (GUI_saveMAT && (i==7 || (i>=0 && i<=6 && !GUI_BF_chs_blocked)))
                    {
                        sprintf(GUI_wtext, "%s\\CH%d.mat", GUI_ImmediateDir, i);
                        sprintf(GUI_wtext2, "%s\\CH%d.mat", GUI_foldername, i);
                        copied = CopyFile((LPCTSTR)GUI_wtext, (LPCTSTR)GUI_wtext2, true);//copy all files between directories in windows
                    }
 				}
                if (copied && foldercreated)	printf("Saved\n");
                else  printf("Problem with creating directory or file copying !\n");
			}
            if (zmq_cmd_msg == "ArinaFile")
            {
                strcpy(GUI_ArinaFileName, zmq_sub_buffer);
            }
            if (zmq_cmd_msg=="scanmode") sscanf(zmq_sub_buffer,"%d",& GUI_scanmode); //write the number in text to variable GUI_...
            if (zmq_cmd_msg == "scanargument")
            {
                sscanf(zmq_sub_buffer, "%f", &tempfnumber);
                GUI_scanargument = tempfnumber;
            }
             if (zmq_cmd_msg == "aspectratio")
            {
                sscanf(zmq_sub_buffer, "%s", &tempfnumber);
                GUI_aspectratio = tempfnumber;
            }
			if (zmq_cmd_msg=="inputampmv") sscanf(zmq_sub_buffer,"%d",& GUI_InputAmpmV); 
			if (zmq_cmd_msg=="Outputampmv")  sscanf(zmq_sub_buffer,"%d",& GUI_OutputAmpmV);
            if (zmq_cmd_msg == "BiasOutputP")  sscanf(zmq_sub_buffer, "%d", &GUI_BiasOutputP);
            if (zmq_cmd_msg=="numberOfSlices")  sscanf(zmq_sub_buffer,"%d",& GUI_numberOfSlices);
			if (zmq_cmd_msg=="ch")  sscanf(zmq_sub_buffer,"%d",& GUI_chosenCH);
			if (zmq_cmd_msg=="cella")  sscanf(zmq_sub_buffer,"%d",& GUI_cella);
            if (zmq_cmd_msg == "LRcella")  sscanf(zmq_sub_buffer, "%d", &GUI_LRcella);
            if (zmq_cmd_msg == "CalibrateBias")  sscanf(zmq_sub_buffer, "%d", &GUI_CalibrateBias);
            if (zmq_cmd_msg == "LowResTimeS")  sscanf(zmq_sub_buffer, "%d", &GUI_LowResTimeS);
            if (zmq_cmd_msg == "ArinaON")  sscanf(zmq_sub_buffer, "%d", &GUI_ArinaON);
            if (zmq_cmd_msg == "BFchsBlocked")  sscanf(zmq_sub_buffer, "%d", &GUI_BF_chs_blocked);
            if (zmq_cmd_msg == "tiltangle")
            {
                sscanf(zmq_sub_buffer, "%f", &tempfnumber);
                GUI_tiltangle = tempfnumber;
            }
            if (zmq_cmd_msg == "LiberTEM")
            {
                strcpy(GUI_LiberTEMcommand, zmq_sub_buffer);
                if (strncmp(GUI_LiberTEMcommand, "CONNECT", 7) == 0)
                {
 
                    //Connect to LiberTEM computer by ZMQ while it continuously running our Pyhton code 
                    if (!LiberTEM_comp_connected)
                    {
                        context3 = zmq_ctx_new();
                        requester3 = zmq_socket(context3, ZMQ_REQ);//I am the requester on socket 5559
                        rc3 = zmq_connect(requester3, "tcp://192.168.100.80:5559"); //bind to LiberTEM server to send commands and get images
                        assert(rc3 == 0);
                        zmq_send(requester3, "hello", 5, ZMQ_DONTWAIT);   //
                        Sleep(200);
                        returnsize = zmq_recv(requester3, zmq_msg_buffer3, 2, ZMQ_DONTWAIT);
                        if (returnsize > 0)
                        {
                            printf("LiberTEM: %s \n", zmq_msg_buffer3);
                        }
                        strcpy(GUI_LiberTEMcommand, "OFF");
                        LiberTEM_comp_connected = true;
                    }
                    else
                    {
                        printf("Connected, if you run the python script on LiberTEM computer\n");
                        strcpy(GUI_LiberTEMcommand, "Hello");
                    }

                }
                //send command to LiberTEM , command START is sent from RecRepMulti unless last command is OFF
                zmq_send(requester3, GUI_LiberTEMcommand, strlen(GUI_LiberTEMcommand),0);
                Sleep(500);
                returnsize = zmq_recv(requester3, zmq_msg_buffer3, 2, 0); //should get OK
                if (returnsize > 0)
                {
                    printf("LiberTEM: %s \n", GUI_LiberTEMcommand);
                }
            }
            if (zmq_cmd_msg == "saveMAT")  sscanf(zmq_sub_buffer, "%d", &GUI_saveMAT);
            if (zmq_cmd_msg == "saveKEY")  sscanf(zmq_sub_buffer, "%d", &GUI_saveMATKEY);
            if (zmq_cmd_msg=="Lothar")
			{
				sscanf(zmq_sub_buffer,"%d",& GUI_Lothar);
				PYTHON_ask_defocus=0;
			}
			if (zmq_cmd_msg=="tomography")  
            {	
                sscanf(zmq_sub_buffer,"%d",& GUI_tomoIndex);
                if (GUI_tomoIndex > -1)
                {
                    if (!GUI_tomography)
                    {
                        GUI_tomography = true;
                        for (int chn = 0; chn < 8; chn++)
                        {
                            fpv[chn] = NULL;
                        }
                    }
                }
			    else
			    {	GUI_tomography=false;
				    GUI_tomoIndex=0;
				    GUI_numberOfSlices=1;
                    //for (int chn = 0; chn < 8; chn++)
                    //{
                    //    if (fpv[chn] != NULL) fclose(fpv[chn]);//close all mrc stacks
                    //}
			    }
			}
            if (zmq_cmd_msg == "SerialEMAlign")
            {
                sscanf(zmq_sub_buffer, "%d", & GUI_SerialEMAligned);
            }

		}
	}
	// Just a timeout - continue the loop
    if (!err)
      continue;

    if (GetDebugVal() > 1) {
      sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, 
        "Select returned with Ready channel: listener %d client %d\n", 
        FD_ISSET(hListener, &readFds), (mHClient[sockInd] != INVALID_SOCKET && 
        FD_ISSET(mHClient[sockInd], &readFds)) ? 1:0);
      DebugToLog(mMessageBuf[sockInd]);
    }

    // There is something to do.  Check the client first (Does ISSET Work?)
    if (mHClient[sockInd] != INVALID_SOCKET && FD_ISSET(mHClient[sockInd], &readFds)) {
      numBytes = recv(mHClient[sockInd], mArgsBuffer[sockInd], mArgBufSize[sockInd], 0);
 
      // Close client on error or disconnect, but allow new connect
      if (numBytes <= 0) {
        ReportErrorAndClose(sockInd, numBytes, "recv from ready client");
      } else {
        memcpy(&numExpected, &mArgsBuffer[sockInd][0], sizeof(int));

        // Reallocate buffer if necessary
        if (numExpected > mArgBufSize[sockInd] - 4) {
          mArgBufSize[sockInd] = ((numExpected + ARGS_BUFFER_CHUNK - 1) /
                                  ARGS_BUFFER_CHUNK) * ARGS_BUFFER_CHUNK;
          mArgsBuffer[sockInd] = (char *)realloc(mArgsBuffer[sockInd],
                                                 mArgBufSize[sockInd]);
          if (!mArgsBuffer[sockInd]) {
            mStartupError[sockInd] = 8;
            ErrorToLog("Failed to reallocate buffer for receiving data");
            return mStartupError[sockInd];
          }
          DebugToLog("Reallocated the argument buffer\n");
        }

        if (!FinishGettingBuffer(sockInd, numBytes, numExpected)) {
          if (GetDebugVal() > 1) {
            sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, 
              "Got %d bytes via recv on socket %d\n", numExpected, mHClient[sockInd]);
            DebugToLog(mMessageBuf[sockInd]);
          }
          mProcessingCommand = true;
          DoProcessCommand(sockInd, numExpected);
          mProcessingCommand = false;
        }
      }
    }

    // Now check the listener, close an existing client if any, get new client
    if (FD_ISSET(hListener, &readFds)) {
      CloseClient(sockInd);
      mHClient[sockInd] = accept(hListener, NULL, NULL);
      if (mHClient[sockInd] == INVALID_SOCKET)
        ReportErrorAndClose(sockInd, SOCKET_ERROR, 
        "accept connection from ready client\n");
      else {
        EitherToLog("", "Accepted connection from client program\n");
      }
    }
  }
 
  return 0;
}

// Close the socket and mark as invalid
void CBaseServer::CloseClient(int sockInd)
{
  if (mHClient[sockInd] == INVALID_SOCKET)
    return;
  sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, 
    "Closing connection to client via socket %d\n", mHClient[sockInd]);
  if (!mCloseForExit[sockInd])
    EitherToLog("", mMessageBuf[sockInd]);
  closesocket(mHClient[sockInd]);
  mHClient[sockInd] = INVALID_SOCKET;
}

// Call the Winsock cleanup function on last uninitialization
void CBaseServer::Cleanup(int sockInd)
{
  int numInit = 0;
  for (int i = 0; i < MAX_SOCK_CHAN; i++)
    if (mInitialized[i])
      numInit++;
  if (numInit == 1)
    WSACleanup();
  mInitialized[sockInd] = false;
}

// Get the rest of the message into the buffer if it is not there yet
int CBaseServer::FinishGettingBuffer(int sockInd, int numReceived, int numExpected)
{
  int numNew, ind;
  while (numReceived < numExpected) {

    // If message is too big for buffer, just get it all and throw away the start
    ind = numReceived;
    if (numExpected > mArgBufSize[sockInd])
      ind = 0;
    numNew = recv(mHClient[sockInd], &mArgsBuffer[sockInd][ind],
                  mArgBufSize[sockInd] - ind, 0);
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

  // Get the function code as the second element of the buffer
  if (numBytes < 8 || numBytes > mArgBufSize[sockInd]) {
    SendArgsBack(sockInd, numBytes < 8 ? -4 : -5);  // Inadequate length or too big
    return 1;
  }
  memcpy(&funcCode, &mArgsBuffer[sockInd][sizeof(int)], sizeof(long));

  // Look up the function code in the table
  ind = 0;
  while (funcTable[ind].funcCode >= 0 && funcTable[ind].funcCode != funcCode)
    ind++;
  if (funcTable[ind].funcCode < 0) {
    sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, "Function code not found: %d\n%s",
      funcCode, upgradeMess);
    ErrorToLog(mMessageBuf[sockInd]);
    SendArgsBack(sockInd, -1);
    return 1;
  }
  
  // Set the variables for receiving and sending arguments.  Add 1 to the longs for
  // the function code coming in and the return value going out
  mNumLongRecv[sockInd] = funcTable[ind].numLongRecv + 1;
  mNumBoolRecv[sockInd] = funcTable[ind].numBoolRecv;
  mNumDblRecv[sockInd] = funcTable[ind].numDblRecv;
  mNumLongSend[sockInd] = funcTable[ind].numLongSend + 1;
  mNumBoolSend[sockInd] = funcTable[ind].numBoolSend;
  mNumDblSend[sockInd] = funcTable[ind].numDblSend;
  mSendLongArray[sockInd] = (funcTable[ind].hasLongArray & 2) > 0;
  needed = sizeof(int) + mNumLongRecv[sockInd] * sizeof(long) +
    mNumBoolRecv[sockInd] * sizeof(BOOL) + mNumDblRecv[sockInd] * sizeof(double);
  if (needed > numBytes) {
    sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, "Command %d  %s not long enough:"
      " needed = %d (4 + %d x 4 + %d x 4 + %d x 8)  numBytes = %d\n", funcCode, 
      funcTable[ind].label, needed, mNumLongRecv[sockInd], mNumBoolRecv[sockInd], 
      mNumDblRecv[sockInd], numBytes);
    ErrorToLog(mMessageBuf[sockInd]);

    SendArgsBack(sockInd, -4);  // Inadequate length, don't even unpack it
    return 1;
  }
  if (UnpackReceivedData(sockInd)) {
    SendArgsBack(sockInd, -5);  // Message too big
    return 1;
  }
  if (funcTable[ind].hasLongArray & 1)
    needAdd = mLongArgs[sockInd][mNumLongRecv[sockInd] - 1] * sizeof(long);
  needed += needAdd;
  if (needed != numBytes) {
    sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, "Command %d  %s wrong length:"
      " needed = %d (4 + %d x 4 + %d x 4 + %d " "x 8 + %d) numBytes = %d\n", funcCode, 
      funcTable[ind].label, needed, mNumLongRecv[sockInd], mNumBoolRecv[sockInd], 
      mNumDblRecv[sockInd], needAdd, numBytes);
    ErrorToLog(mMessageBuf[sockInd]);
    SendArgsBack(sockInd, -6);   // Wrong length
    return 1;
  }
  if (GetDebugVal() > 1) {
    sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, "Processing command %d  %s\n", 
      funcCode, funcTable[ind].label);
    DebugToLog(mMessageBuf[sockInd]);
  }
  return 0;
}


// Send a buffer back, in chunks if necessary
int CBaseServer::SendBuffer(int sockInd, char *buffer, int numBytes)
{
  int numTotalSent = 0;
  int numToSend, numSent;
  if (GetDebugVal() > 1) {
    sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, "In SendBuffer, socket %d, "
      "sending %d bytes\n",  mHClient[sockInd], numBytes);
    DebugToLog(mMessageBuf[sockInd]);
  }
  while (numTotalSent < numBytes) {
    numToSend = numBytes - numTotalSent;
    if (numToSend > mChunkSize)
      numToSend = mChunkSize;
    if (GetDebugVal() > 1) {
      sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, "Going to send %d bytes to "
      "socket %d\n", numToSend, mHClient[sockInd]);
      DebugToLog(mMessageBuf[sockInd]);
    }
    numSent = send(mHClient[sockInd], &buffer[numTotalSent], numToSend, 0);
    if (numSent < 0) {
      ReportErrorAndClose(sockInd, numSent, "send a chunk of bytes back\n");
      return 1;
    }
    numTotalSent += numSent;
  }
  return 0;
}

// Close the connection upon error; report it unless it is clearly a SerialEM disconnect
void CBaseServer::ReportErrorAndClose(int sockInd, int retval, const char *message)
{
  if (retval == SOCKET_ERROR) {
    mLastWSAerror[sockInd] = WSAGetLastError();
    sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, "WSA Error %d on call to %s\n", 
      mLastWSAerror[sockInd], message);
    if (mLastWSAerror[sockInd] == WSAECONNRESET)
      DebugToLog(mMessageBuf[sockInd]);
    else
      ErrorToLog(mMessageBuf[sockInd]);
  } else {
    sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, "retval %d on call to %s", 
      retval, message);
    ErrorToLog(mMessageBuf[sockInd]);
  }
  CloseClient(sockInd);
}

// Close up on error or signal from plugin
void CBaseServer::CloseOnExitOrSelectError(int sockInd, int err)
{
  DebugToLog("Closing socket\n");
  CloseClient(sockInd);
  closesocket(mHListener[sockInd]);
  //if (mCloseForExit)
    //Cleanup();
  if (err < 0) {
    mLastWSAerror[sockInd] = WSAGetLastError();
    mStartupError[sockInd] = 7;
    sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, 
      "WSA Error %d on select command\n");
    ErrorToLog(mMessageBuf[sockInd]);
  }
}

// Wait for the client to acknowledge receipt of a superchunk of image
int CBaseServer::ListenForHandshake(int sockInd, int superChunk)
{
  struct timeval tv;
  int numBytes, err, numExpected, command;
  fd_set readFds;      // file descriptor list for select()
  tv.tv_sec = 0;
  tv.tv_usec = superChunk / 5;    // This is 5 MB /sec

  FD_ZERO(&readFds);
  FD_SET(mHClient[sockInd], &readFds);
  err = select(1, &readFds, NULL, NULL, &tv);
  if (err < 0 || mCloseForExit[sockInd]) {
    CloseOnExitOrSelectError(sockInd, err);
    return mStartupError[sockInd];
  }

  // A timeout - close client for this so client fails
  if (!err) {
    ReportErrorAndClose(sockInd, 0, "timeout on handshake from client");
    return 1;
  }

  numBytes = recv(mHClient[sockInd], mArgsBuffer[sockInd], mArgBufSize[sockInd], 0);

  // Close client on error or disconnect or too few bytes or anything wrong
  memcpy(&numExpected, &mArgsBuffer[sockInd][0], sizeof(int));
  memcpy(&command, &mArgsBuffer[sockInd][4], sizeof(int));
  if (command != mHandshakeCode || numExpected != 8 || numBytes != 8) {
    ReportErrorAndClose(sockInd, numBytes, "recv handshake from ready client");
    return 1;
  }
  return 0;
}

// Send the arguments back, packing the return value in the first long
int CBaseServer::SendArgsBack(int sockInd, int retval)
{
  mLongArgs[sockInd][0] = retval;

  // Emergency error code is negative, just send one word back
  if (retval < 0) {
    mNumLongSend[sockInd] = 1;
    mNumBoolSend[sockInd] = 0;
    mNumDblSend[sockInd] = 0;
    mSendLongArray[sockInd] = false;
  }
  if (PackDataToSend(sockInd)) {
    ErrorToLog("DATA BUFFER NOT BIG ENOUGH TO SEND REPLY TO SERIALEM");
    SendArgsBack(sockInd, -3);
    return 1;
  }
  return SendBuffer(sockInd, mArgsBuffer[sockInd], mNumBytesSend[sockInd]);
}

// Send the arguments from an image acquisition back then send the image if there is no
// error
void CBaseServer::SendImageBack(int sockInd, int retval, short *imArray, int bytesPerPixel)
{
  int numChunks, chunkSize, numToSend, numLeft, err, imSize, totalSent = 0;

  // determine number of superchunks and send that back as fourth long
  imSize = mLongArgs[sockInd][1] * bytesPerPixel;
  numChunks = (imSize + mSuperChunkSize - 1) / mSuperChunkSize;
  mLongArgs[sockInd][4] = numChunks;
  err = SendArgsBack(sockInd, retval);
  sprintf_s(mMessageBuf[sockInd], MESS_ERR_BUFF_SIZE, "retval = %d, err sending args %d,"
    " sending image %d in %d chunks\n", retval, err, imSize, numChunks);
  DebugToLog(mMessageBuf[sockInd]);
  if (!err && !retval) {

    // Loop on the chunks until done, getting acknowledgement after each
    numLeft = imSize;
    chunkSize = (imSize + numChunks - 1) / numChunks;
    while (totalSent < imSize) {
      numToSend = chunkSize;
      if (chunkSize > imSize - totalSent)
        numToSend = imSize - totalSent;
      if (SendBuffer(sockInd, (char *)imArray + totalSent, numToSend))
        break;
      totalSent += numToSend;
      if (totalSent < imSize && ListenForHandshake(sockInd, numToSend))
        break;
    }
  }
  delete [] imArray;
}

// Unpack the received argument buffer, skipping over the first word, the byte count
int CBaseServer::UnpackReceivedData(int sockInd)
{
  int numBytes, numUnpacked = sizeof(int);
  if (mNumLongRecv[sockInd] > MAX_LONG_ARGS || mNumBoolRecv[sockInd] > MAX_BOOL_ARGS || 
    mNumDblRecv[sockInd] > MAX_DBL_ARGS)
    return 1;
  numBytes = mNumLongRecv[sockInd] * sizeof(long);
  memcpy(mLongArgs[sockInd], &mArgsBuffer[sockInd][numUnpacked], numBytes);
  numUnpacked += numBytes;
  numBytes = mNumBoolRecv[sockInd] * sizeof(BOOL);
  if (numBytes)
    memcpy(mBoolArgs[sockInd], &mArgsBuffer[sockInd][numUnpacked], numBytes);
  numUnpacked += numBytes;
  numBytes = mNumDblRecv[sockInd] * sizeof(double);
  if (numBytes)
    memcpy(mDoubleArgs[sockInd], &mArgsBuffer[sockInd][numUnpacked], numBytes);
  numUnpacked += numBytes;

  // Here is the starting address of whatever comes next for the few expecting it
  mLongArray[sockInd] = (long *)(&mArgsBuffer[sockInd][numUnpacked]);
  return 0;
}


#define ADD_TO_BUFFER(a) \
  if (numAdd + mNumBytesSend[sockInd] > mArgBufSize[sockInd]) \
    return 1; \
  memcpy(&mArgsBuffer[sockInd][mNumBytesSend[sockInd]], a, numAdd); \
  mNumBytesSend[sockInd] += numAdd;

// Pack the data into the argument buffer as longs, BOOLS, doubles
int CBaseServer::PackDataToSend(int sockInd)
{
  int numAdd;
  mNumBytesSend[sockInd] = sizeof(int);
  if (mNumLongSend[sockInd]) {
    numAdd = mNumLongSend[sockInd] * sizeof(long);
    ADD_TO_BUFFER(mLongArgs[sockInd]);
  }
  if (mNumBoolSend[sockInd]) {
    numAdd = mNumBoolSend[sockInd] * sizeof(BOOL);
    ADD_TO_BUFFER(&mBoolArgs[sockInd]);
  }
  if (mNumDblSend[sockInd]) {
    numAdd = mNumDblSend[sockInd] * sizeof(double);
    ADD_TO_BUFFER(&mDoubleArgs[sockInd]);
  }
  if (mSendLongArray[sockInd]) {
    numAdd = mLongArgs[sockInd][mNumLongSend[sockInd] - 1] * sizeof(long);
    ADD_TO_BUFFER(mLongArray[sockInd]);
  }

  // Put the number of bytes at the beginning of the message
  memcpy(&mArgsBuffer[sockInd][0], &mNumBytesSend[sockInd], sizeof(int));
  return 0;
}
