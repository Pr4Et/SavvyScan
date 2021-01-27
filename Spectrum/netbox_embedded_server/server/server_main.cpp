/*
*******************************************************************************

server_main.cpp                                               (c) Spectrum GmbH

*******************************************************************************

Documentation for the API as well as a detailed description of the hardware
can be found in the manual for each device which can be found on our website:
https://www.spectrum-instrumentation.com/en/downloads

Further information can be found online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/knowledge-base-overview

*******************************************************************************
*/

#ifdef WIN32
#   include <winsock2.h> // before everything else!
#   include <Ws2tcpip.h>



#   undef min
#   undef max

typedef int ssize_t;
#else
#   include <sys/socket.h>
#   include <netinet/in.h>
#   include <arpa/inet.h>
#   include <unistd.h>
#   include <errno.h>
#   include <sys/types.h>
#   include <netdb.h>

typedef int SOCKET;
#define SOCKET_ERROR -1
#endif

#include <cstring>
#include <ctime>

#include <iostream>
#include <iomanip>
#include <algorithm>

#include <cstdlib>

// ----- include standard driver header from library -----
#include "../../c_header/dlltyp.h"
#include "../../c_header/regs.h"
#include "../../c_header/spcerr.h"
#include "../../c_header/spcm_drv.h"

#include "../../common/ostools/spcm_oswrap.h"
#include "../../common/ostools/spcm_ostools.h"
#include "../../common/ostools/spcm_network_winLin.h"


#define SERVER_PORT     22927
#define START 0
#define STOP  1

bool bDataCanBeReadFromSocket (SOCKET lSocket)
    {
    fd_set          fdset;
    struct timeval  stNoWait;

    FD_ZERO (&fdset);
    FD_SET (lSocket, &fdset);
    memset (&stNoWait, 0, sizeof (stNoWait));

    // returns immediately if data can be read from socket without blocking
    lSelect (lSocket + 1, &fdset, NULL, NULL, &stNoWait);
    if (FD_ISSET (lSocket, &fdset))
        return true;
    else
        return false;
    }


bool bSetupCard ()
    {

    return true;
    }

int main ()
    {
    dwInitNetwork ();

    SOCKET listenfd = lSocket (AF_INET, SOCK_STREAM, 0);

    struct sockaddr_in serv_addr;
    memset(&serv_addr, '0', sizeof(serv_addr));
    serv_addr.sin_family = AF_INET;
    serv_addr.sin_addr.s_addr = htonl(INADDR_ANY);
    serv_addr.sin_port = htons(SERVER_PORT);

    lBind (listenfd, (struct sockaddr*)&serv_addr, sizeof(serv_addr));

    lListen (listenfd, 10);

    SOCKET lAcceptedSocket = INVALID_SOCKET;
    fd_set rfds;
    SOCKET lMaxSockFd = listenfd;
    char acBuffer[1024];
    while (true)
        {
        FD_ZERO (&rfds);
        lMaxSockFd = -1;
        if (lAcceptedSocket == INVALID_SOCKET)
            {
            std::cout << "Waiting for client..." << std::endl;
            FD_SET (listenfd, &rfds);
            lMaxSockFd = listenfd;
            }
        else
            {
            FD_SET (lAcceptedSocket, &rfds);
            lMaxSockFd = std::max (lMaxSockFd, lAcceptedSocket);
            }

        select ((int)(lMaxSockFd + 1), &rfds, NULL, NULL, NULL);

        if (FD_ISSET(listenfd, &rfds))
            {
            std::cout << "Accepted client..." << std::endl;
            lAcceptedSocket = accept (listenfd, (struct sockaddr*)NULL, NULL);
            }

        // ----- read socket and print received data -----
        if (FD_ISSET (lAcceptedSocket, &rfds))
            {
            // ----- read data from sender -----
            unsigned dwReadBytes = dwRead (lAcceptedSocket, acBuffer, sizeof (acBuffer) - 1);
            if (dwReadBytes == 0)
                {
                // connection closed
                lClose (lAcceptedSocket);
                lAcceptedSocket = INVALID_SOCKET;
                }
            else if (dwReadBytes == SOCKET_ERROR)
                {
                // error occured
                std::cerr << "NETWORK ERROR" << std::endl;
                lClose (lAcceptedSocket);
                lAcceptedSocket = INVALID_SOCKET;
                }
            else
                {
                uint32 dwMagic = *reinterpret_cast < uint32* > (acBuffer);
                switch (dwMagic)
                    {
                    case START:
                        {
                        drv_handle  hCard;
                        int32       lCardType, lSerialNumber, lFncType;
                        int16*      pnData;
                        char        szErrorTextBuffer[ERRORTEXTLEN];
                        uint32      dwError;
                        int32       lStatus;
                        int64       llAvailUser, llPCPos;
                        uint64      qwTotalMem = 0;
                        uint64      qwToTransfer = MEGA_B(64);

                        // settings for the FIFO mode buffer handling
                        int64       llBufferSize =  MEGA_B(4);
                        int32       lNotifySize =   KILO_B(16);


                        // open card
                        // uncomment the second line and replace the IP address to use remote
                        // cards like in a digitizerNETBOX
                        hCard = spcm_hOpen ("/dev/spcm0");
                        if (!hCard)
                            {
                            std::cerr << "no card found..." << std::endl;
                            return EXIT_FAILURE;
                            }

#ifndef WIN32
                        system ("netbox_led_client conngreen=1");
#endif

                        // read type, function and sn and check for A/D card
                        spcm_dwGetParam_i32 (hCard, SPC_PCITYP,         &lCardType);
                        spcm_dwGetParam_i32 (hCard, SPC_PCISERIALNO,    &lSerialNumber);
                        spcm_dwGetParam_i32 (hCard, SPC_FNCTYPE,        &lFncType);

                        std::cout << "Found: ";
                        switch (lCardType & TYP_SERIESMASK)
                            {
                            case TYP_M2ISERIES:     std::cout << "M2i." << std::hex << (lCardType & TYP_VERSIONMASK);            break;
                            case TYP_M2IEXPSERIES:  std::cout << "M2i." << std::hex << (lCardType & TYP_VERSIONMASK) << "-Exp";  break;
                            case TYP_M3ISERIES:     std::cout << "M3i." << std::hex << (lCardType & TYP_VERSIONMASK);            break;
                            case TYP_M3IEXPSERIES:  std::cout << "M3i." << std::hex << (lCardType & TYP_VERSIONMASK) << "-Exp";  break;
                            case TYP_M4IEXPSERIES:  std::cout << "M4i." << std::hex << (lCardType & TYP_VERSIONMASK) << "-x8";   break;
                            case TYP_M4XEXPSERIES:  std::cout << "M4x." << std::hex << (lCardType & TYP_VERSIONMASK) << "-x4";   break;
                            default:                std::cout << "unknown type";                                               break;
                            }
                        std::cout << " sn " << std::setw (5) << lSerialNumber << std::endl;

                        if (lFncType != SPCM_TYPE_AI)
                            {
                            std::cerr << "Card not supported by example" << std::endl;
                            return EXIT_FAILURE;
                            }


                        // do a simple standard setup
                        spcm_dwSetParam_i32 (hCard, SPC_CHENABLE,        CHANNEL0);             // just 1 channel enabled
                        spcm_dwSetParam_i32 (hCard, SPC_PATH0,           0);
                        spcm_dwSetParam_i32 (hCard, SPC_AMP0,            1000);
                        spcm_dwSetParam_i32 (hCard, SPC_CARDMODE,        SPC_REC_FIFO_MULTI);   // single FIFO mode
                        spcm_dwSetParam_i32 (hCard, SPC_SEGMENTSIZE,     lNotifySize);
                        spcm_dwSetParam_i32 (hCard, SPC_POSTTRIGGER,     lNotifySize - 64);
                        spcm_dwSetParam_i32 (hCard, SPC_PRETRIGGER,      1024);                 // 1k of pretrigger data at start of FIFO mode
                        spcm_dwSetParam_i32 (hCard, SPC_LOOPS,           0);                    // endless FIFO
                        spcm_dwSetParam_i32 (hCard, SPC_TIMEOUT,         5000);                 // timeout 5 s
                        spcm_dwSetParam_i32 (hCard, SPC_TRIG_ORMASK,     SPC_TMASK_NONE);       // trigger set to channel
                        spcm_dwSetParam_i32 (hCard, SPC_TRIG_CH_ORMASK0, SPC_TMASK0_CH0);       // trigger set to channel0
                        spcm_dwSetParam_i32 (hCard, SPC_TRIG_CH0_MODE,   SPC_TM_POS);
                        spcm_dwSetParam_i32 (hCard, SPC_TRIG_CH0_LEVEL0, 20);
                        spcm_dwSetParam_i32 (hCard, SPC_TRIG_ANDMASK,    0);                    // ...
                        spcm_dwSetParam_i32 (hCard, SPC_CLOCKMODE,       SPC_CM_INTPLL);        // clock mode internal PLL

                        // we try to set the samplerate to 100 kHz (M2i) or 20 MHz (M3i and M4i) on internal PLL, no clock output
                        if (((lCardType & TYP_SERIESMASK) == TYP_M2ISERIES) || ((lCardType & TYP_SERIESMASK) == TYP_M2IEXPSERIES))
                            spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE, KILO(100));
                        else
                            spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE, MEGA(20));

                        spcm_dwSetParam_i32 (hCard, SPC_CLOCKOUT,       0);                     // no clock output


                        // define the data buffer
                        pnData = (int16*) pvAllocMemPageAligned ((uint64) llBufferSize);
                        if (!pnData)
                            {
                            std::cerr << "memory allocation failed" << std::endl;
                            spcm_vClose (hCard);
                            return false;
                            }

                        spcm_dwDefTransfer_i64 (hCard, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, lNotifySize, pnData, 0, llBufferSize);


                        // start everything
                        dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_DATA_STARTDMA);

                        // check for error
                        if (dwError != ERR_OK)
                            {
                            spcm_dwGetErrorInfo_i32 (hCard, NULL, NULL, szErrorTextBuffer);
                            std::cout << szErrorTextBuffer << std::endl;
                            vFreeMemPageAligned (pnData, (uint64) llBufferSize);
                            spcm_vClose (hCard);
                            return 0;
                            }

                        // ----- start data processing loop -----
                        while (qwTotalMem < qwToTransfer)
                            {
                            // ----- check if client has been clocked or if it has sent an abort command -----
                            bool bKeepRunning = true;
                            if (bDataCanBeReadFromSocket (lAcceptedSocket))
                                {
                                unsigned dwReadBytes = dwRead (lAcceptedSocket, acBuffer, sizeof (acBuffer) - 1);
                                if (dwReadBytes == 0)
                                    {
                                    // connection closed
                                    lClose (lAcceptedSocket);
                                    lAcceptedSocket = INVALID_SOCKET;
                                    bKeepRunning = false;
                                    }
                                else if (dwReadBytes == SOCKET_ERROR)
                                    {
                                    // error occured
                                    lClose (lAcceptedSocket);
                                    lAcceptedSocket = INVALID_SOCKET;
                                    bKeepRunning = false;
                                    }
                                else
                                    {
                                    dwMagic = *reinterpret_cast < uint32* > (acBuffer);
                                    if (dwMagic == STOP)
                                        {
                                        bKeepRunning = false;
                                        break;
                                        }
                                    }
                                }
                            if (!bKeepRunning)
                                {
                                // send the stop command
                                dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_STOP | M2CMD_DATA_STOPDMA);

                                spcm_vClose (hCard);
                                break;
                                }

                            if ((dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_DATA_WAITDMA)) != ERR_OK)
                                {
                                if (dwError == ERR_TIMEOUT)
                                    std::cout << "... Timeout" << std::endl;
                                else
                                    std::cerr << "... Error: " << dwError << std::endl;
                                break;
                                }

                            else
                                {
                                spcm_dwGetParam_i32 (hCard, SPC_M2STATUS,             &lStatus);
                                spcm_dwGetParam_i64 (hCard, SPC_DATA_AVAIL_USER_LEN,  &llAvailUser);
                                spcm_dwGetParam_i64 (hCard, SPC_DATA_AVAIL_USER_POS,  &llPCPos);

                                if (llAvailUser >= lNotifySize)
                                    {
#ifndef WIN32
                                    system ("netbox_led_client armgreen=1");
#endif
                                    uint32 dwNumSamples = lNotifySize / sizeof (int16);
                                    qwTotalMem += dwNumSamples;

                                    char* pcStartPos = (char*)pnData + llPCPos;
                                    char* pcEndPos   = pcStartPos + lNotifySize;

                                    // this is the point to do anything with the data
                                    // since notify size is multiple of buffer size, we don't care about buffer wrap-around
                                    int16 nMax = *std::max_element ((int16*)pcStartPos, (int16*)pcEndPos);
                                    int16 nMin = *std::min_element ((int16*)pcStartPos, (int16*)pcEndPos);
                                    int32 lAvg = 0;
                                    for (unsigned dwIdx = 0; dwIdx < dwNumSamples; dwIdx++)
                                        {
                                        lAvg += *(((int16*)pcStartPos) + dwIdx);
                                        }
                                    lAvg /= dwNumSamples;

                                    char acText[128];
                                    snprintf (acText, 128, "Maximum: %hd  Minimum: %hd  Average: %d\n", nMax, nMin, lAvg);

                                    ssize_t dwWritten_bytes = dwWrite (lAcceptedSocket, acText, strlen (acText));

                                    spcm_dwSetParam_i32 (hCard, SPC_DATA_AVAIL_CARD_LEN,  lNotifySize);
#ifndef WIN32
                                    system ("netbox_led_client armgreen=0");
#endif
                                    }
                                }
                            }
                        break;
                        }
                    default:
                        // unknown magic number. Should not reach here
                        break;
                    }
                }
            }
        }

    return EXIT_SUCCESS;
    }
