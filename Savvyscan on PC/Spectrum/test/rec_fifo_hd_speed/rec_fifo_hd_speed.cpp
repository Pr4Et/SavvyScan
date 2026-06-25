/*
**************************************************************************

rec_fifo_hd_speed.cpp                                    (c) Spectrum GmbH

**************************************************************************

this example supports all acquisition cards

Does FIFO acquisition to hard disk to test the maximum writing performance
of the hard disk

This program only runs under Windows as it uses some windows specific API
calls for data writing, time measurement and key checking

Documentation for the API as well as a detailed description of the hardware
can be found in the manual for each device which can be found on our website:
https://www.spectrum-instrumentation.com/en/downloads

Further information can be found online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/knowledge-base-overview

**************************************************************************
*/



// ----- standard c include files -----
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <conio.h>

// ----- include of common example librarys -----
#include "../common/spcm_lib_card.h"
#include "../common/spcm_lib_data.h"
#include "../common/spcm_lib_thread.h"


// ----- global setup for the run (can be changed interactively) -----
int64   g_llSamplingRate =  MEGA(20);
int32   g_lNotifySize =     KILO_B(1024);
int64   g_llBufferSize =     MEGA_B(128);
bool    g_bThread =         false;
uint64  g_qwChannelEnable = 1;
uint32  g_dwUpdateBuffers = 1;
uint32  g_dwUpdateCount =   0;
enum    {eStandard, eHDSpeedTest, eSpeedTest} g_eMode = eStandard;

#define FILENAME "stream_test"



/*
**************************************************************************
bDoCardSetup: setup matching the calculation routine
**************************************************************************
*/

bool bDoCardSetup (ST_SPCM_CARDINFO *pstCard)
    {


    // FIFO mode setup, we run continuously and have 16 samples of pre data before trigger event
    // all available channels are activated
    bSpcMSetupModeRecFIFOSingle (pstCard, g_qwChannelEnable, 16);

    // we try to set the samplerate on internal PLL, no clock output
    if (g_llSamplingRate > pstCard->llMaxSamplerate)
        g_llSamplingRate = pstCard->llMaxSamplerate;

    // for M4i.44xx series we activate the fine granularity clock setup mode
    if (pstCard->bM4i && ((pstCard->lCardType & TYP_FAMILYMASK) == 0x4400))
        spcm_dwSetParam_i64 (pstCard->hDrv, SPC_SPECIALCLOCK, 1);

    bSpcMSetupClockPLL (pstCard, g_llSamplingRate, false);
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstCard->llSetSamplerate / MEGA(1));

    // we set software trigger, no trigger output
    bSpcMSetupTrigSoftware (pstCard, false);

    return pstCard->bSetError;
    }



/*
**************************************************************************
Working routine data
**************************************************************************
*/


struct ST_WORKDATA
    {
    int64           llWritten;
    HANDLE          hFile;
    char            szFileName[100];
    LARGE_INTEGER   uStartTime;
    LARGE_INTEGER   uLastTime;
    LARGE_INTEGER   uHighResFreq;
    };



/*
**************************************************************************
Setup working routine
**************************************************************************
*/

bool bWorkInit (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA* pstWorkData = (ST_WORKDATA*) pvWorkData;

    // setup for the transfer, to avoid overrun we use quite large blocks as this has a better throughput to hard disk
    pstBufferData->dwDataBufLen = g_llBufferSize;
    pstBufferData->dwDataNotify = g_lNotifySize;

    // setup for the work
    pstWorkData->llWritten = 0;

    sprintf (pstWorkData->szFileName, "%s.bin", FILENAME);

    printf ("\n");
    printf ("Written      HW-Buf      SW-Buf   Average   Current\n-----------------------------------------------------\n");
    if ((g_eMode == eStandard) || (g_eMode == eHDSpeedTest))
        pstWorkData->hFile = CreateFile (pstWorkData->szFileName, 
            GENERIC_WRITE, 
            FILE_SHARE_READ | FILE_SHARE_WRITE, 
            NULL, 
            CREATE_ALWAYS, 
            FILE_ATTRIBUTE_NORMAL | FILE_FLAG_NO_BUFFERING, 
            NULL);

    QueryPerformanceFrequency (&pstWorkData->uHighResFreq);
    pstWorkData->uStartTime.QuadPart = 0;

    return ((pstWorkData->hFile != NULL) || (g_eMode == eSpeedTest));
    }



/*
**************************************************************************
bWorkDo: stores data to hard disk
**************************************************************************
*/

bool bWorkDo (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA*    pstWorkData = (ST_WORKDATA*) pvWorkData;
    uint32          dwWritten;
    int64           llBufferFillPromille;
    LARGE_INTEGER   uTime;
    double          dAverageTime = 0;
    double          dLastTime = 0;
    double          dAverageSpeed = 0;
    double          dLastSpeed = 0;
    


    // first call will initialize all time measurings
    QueryPerformanceCounter (&uTime);
    if (pstWorkData->uStartTime.QuadPart == 0)
        {
        pstWorkData->uStartTime.QuadPart = uTime.QuadPart;
        pstWorkData->uLastTime.QuadPart =  uTime.QuadPart;
        }

    // calc transfer speed
    else
        {
        dAverageTime =  (double) (uTime.QuadPart - pstWorkData->uStartTime.QuadPart) / pstWorkData->uHighResFreq.QuadPart;
        dLastTime =     (double) (uTime.QuadPart - pstWorkData->uLastTime.QuadPart) / pstWorkData->uHighResFreq.QuadPart;
        dAverageSpeed = (double) pstWorkData->llWritten / dAverageTime / MEGA_B(1);
        dLastSpeed =    (double) pstBufferData->dwDataNotify / dLastTime / MEGA_B(1);
        pstWorkData->uLastTime.QuadPart = uTime.QuadPart;
        }


    // write the data and count the samples
    if (g_eMode == eSpeedTest)
        dwWritten = pstBufferData->dwDataNotify;
    else
        WriteFile (pstWorkData->hFile, pstBufferData->pvDataCurrentBuf, pstBufferData->dwDataNotify, &dwWritten, NULL);

    pstWorkData->llWritten += dwWritten;
    if (dwWritten != pstBufferData->dwDataNotify)
        {
        printf ("\nData Write error\n");
        return false;
        }

    // current status
    if (--g_dwUpdateCount == 0)
        {
        g_dwUpdateCount = g_dwUpdateBuffers;
        spcm_dwGetParam_i64 (pstBufferData->pstCard->hDrv, SPC_FILLSIZEPROMILLE, &llBufferFillPromille);

        printf ("\r");
        if (pstBufferData->llDataTransferred > GIGA_B(1))
            printf ("%7.2lf GB", (double) pstBufferData->llDataTransferred / GIGA_B(1));
        else
            printf ("%7.2lf MB", (double) pstBufferData->llDataTransferred / MEGA_B(1));

        printf (" %6.1lf %%", (double) llBufferFillPromille / 10.0);

        printf ("    %6.1lf %%", 100.0 * (double) pstBufferData->dwDataAvailBytes / pstBufferData->dwDataBufLen);

        // print transfer speed
        printf ("   %6.2lf MB/s", dAverageSpeed);
        printf ("   %6.2lf MB/s", dLastSpeed);
        }

    pstBufferData->dwDataAvailBytes = pstBufferData->dwDataNotify;

    return true;
    }




/*
**************************************************************************
vWorkClose: Close the work and clean up
            For speed reason is the bKeyAbortCheck function (with 
            kbhit() inside) not used!
**************************************************************************
*/

void vWorkClose (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA* pstWorkData = (ST_WORKDATA*) pvWorkData;

    if (pstWorkData->hFile && (g_eMode != eSpeedTest))
        CloseHandle (pstWorkData->hFile);
    }


/*
**************************************************************************
bKeyCheckAsync
**************************************************************************
*/
int g_nKeyPress;
bool bKeyCheckAsync (void* , ST_BUFFERDATA*)
    {
    return (g_nKeyPress != GetAsyncKeyState(VK_ESCAPE));
    }



/*
**************************************************************************
bSetup: returns true if start, false if abort
**************************************************************************
*/

bool bSetup (ST_SPCM_CARDINFO* pstCard)
    {
    double dTmp;
    uint32 dwTmp;
    int32  lChannels;
    char   szErrorText[ERRORTEXTLEN], szNameBuffer[100];
    uint64 qwContBufLen;

    // read out cont buf len and set default buffer size to it
    void* pvTmp;
    spcm_dwSetParam_i64 (pstCard->hDrv, SPC_M2CMD, M2CMD_CARD_RESET);
    spcm_dwGetContBuf_i64 (pstCard->hDrv, SPCM_BUF_DATA, &pvTmp, &qwContBufLen);
    if (qwContBufLen > 0)
        g_llBufferSize = (int32) qwContBufLen;

    while (1)
        {
        printf ("\n\n");
        printf ("Current Setup\n-------------\n");
        printf ("          Card Selection:   %s", pszSpcMPrintCardInfo (pstCard, szNameBuffer, sizeof (szNameBuffer), false)); 
        printf ("I ....... Interface speed:  ");
        switch (g_eMode)
            {
            case eStandard:    printf ("Normal FIFO mode to HD\n"); break;
            case eHDSpeedTest: printf ("Max PCI/PCIe interface speed to HD\n"); break;
            case eSpeedTest:   printf ("Max PCI/PCIe interface speed only\n"); break;
            }
        printf ("B ....... Buffer Size:      %.2lf MByte (Continuous Buffer: %d MByte)\n", (double) g_llBufferSize / MEGA_B(1), (int32) (qwContBufLen / MEGA_B(1)));
        printf ("N ....... Notify Size:      %d kByte\n", g_lNotifySize / KILO_B(1));
        if (g_eMode == eStandard)
            {
            printf ("S ....... Sampling Rate:    %.2lf MS/s\n", (double) g_llSamplingRate / MEGA(1));
            printf ("T ....... Thread Mode:      %s\n", g_bThread ? "on" : "off");
            printf ("C ....... Channel Enable:   %x\n", g_qwChannelEnable);
            }
        printf ("Enter ... Start Test\n");
        printf ("Q ....... Quit\n");

        // for M4i.44xx series we activate the fine granularity clock setup mode
        if (pstCard->bM4i && ((pstCard->lCardType & TYP_FAMILYMASK) == 0x4400))
            spcm_dwSetParam_i64 (pstCard->hDrv, SPC_SPECIALCLOCK, 1);

        spcm_dwSetParam_i64 (pstCard->hDrv, SPC_CHENABLE, g_qwChannelEnable);
        spcm_dwSetParam_i64 (pstCard->hDrv, SPC_SAMPLERATE, g_llSamplingRate);
        spcm_dwGetParam_i64 (pstCard->hDrv, SPC_SAMPLERATE, &g_llSamplingRate);
        spcm_dwGetParam_i32 (pstCard->hDrv, SPC_CHCOUNT, &lChannels);
        spcm_dwSetParam_i32 (pstCard->hDrv, SPC_TEST_FIFOSPEED, (g_eMode != eStandard) ? 1 : 0);



        if (spcm_dwGetErrorInfo_i32 (pstCard->hDrv, NULL, NULL, szErrorText) != ERR_OK)
            printf ("\nSetup Error:\n------------\n%s\n\n", szErrorText);
        else
            {
            double dTransferSpeed;
                if (pstCard->eCardFunction == AnalogIn)
                    dTransferSpeed = (double)  g_llSamplingRate * lChannels * pstCard->lBytesPerSample;
                else
                    dTransferSpeed = (double)  g_llSamplingRate * lChannels / 8;

            if (g_eMode == eStandard)
                {
                printf ("          Sampling Rate adjusted to: %.2lf MS/s\n", (double) g_llSamplingRate / MEGA(1));
                printf ("          Transfer Speed: %.2lf MByte/s\n", dTransferSpeed / MEGA_B(1));
                }
            else
                printf ("          Transfer Speed: max\n");

            // calc the display update rate in buffers to x/second to keep display overhead small
            g_dwUpdateBuffers = (uint32) (dTransferSpeed / g_lNotifySize / 4);
            if (g_dwUpdateBuffers < 1)
                g_dwUpdateBuffers = 1;
            g_dwUpdateCount = g_dwUpdateBuffers;
            }
        printf ("\n");


        switch (_getch())
            {
            case 'q': 
            case 'Q':
                return false;

            case 13: 
                if (g_llBufferSize <= (int32) qwContBufLen)
                    printf ("\n***** Continuous Buffer from Kernel Driver used *****\n\n");
                return true;

            case 'i':
            case 'I':
                switch (g_eMode)
                    {
                    case eStandard:     g_eMode = eHDSpeedTest; break;
                    case eHDSpeedTest:  g_eMode = eSpeedTest; break;
                    case eSpeedTest:    g_eMode = eStandard; break;
                    }
                break;

            case 't':
            case 'T':
                g_bThread = !g_bThread;
                break;

            case 's':
            case 'S':
                printf ("Sampling Rate (MS/s): ");
                scanf ("%lf", &dTmp);
                g_llSamplingRate = (int32) (dTmp * MEGA(1));
                break;

            case 'b':
            case 'B':
                printf ("Buffer Size (MByte): ");
                scanf ("%lf", &dTmp);
                g_llBufferSize = (int32) (dTmp * MEGA_B(1));
                break;

            case 'n':
            case 'N':
                printf ("Notify Size (kByte): ");
                scanf ("%lf", &dTmp);
                g_lNotifySize = (int32) (dTmp * KILO_B(1));
                break;

            case 'c':
            case 'C':
                printf ("Channel Enable Mask (hex): ");
                scanf ("%x", &dwTmp);
                g_qwChannelEnable = dwTmp;
                break;

            }
        }
    }



/*
**************************************************************************
main 
**************************************************************************
*/

int main ()
    {
    char                szBuffer[1024];     // a character buffer for any messages
    ST_SPCM_CARDINFO    astCard[MAXBRD];    // info structure of my card
    ST_BUFFERDATA       stBufferData;       // buffer and transfer definitions
    ST_WORKDATA         stWorkData;         // work data for the working functions
    int32               lCardIdx =      0;
    int32               lCardCount =    0;

    // ------------------------------------------------------------------------
    // init cards, get some information and print it
    for (lCardIdx = 0; lCardIdx < MAXBRD; lCardIdx++)
        {
        if (bSpcMInitCardByIdx (&astCard[lCardCount], lCardIdx))
            {
            printf (pszSpcMPrintCardInfo (&astCard[lCardCount], szBuffer, sizeof (szBuffer)));
            printf ("\n");
            lCardCount++;
            }
        }
    if (lCardCount == 0)
        {
        printf ("No Spectrum card found...\n");
        return 0;
        }

    // if we have more than one card we make the selection now
    if (lCardCount > 1)
        {
        do
            {
            printf ("\n");
            printf ("Please select the card to test:\n");
            printf ("-------------------------------\n");
            for (lCardIdx = 0; lCardIdx < lCardCount; lCardIdx++)
                printf ("%d ..... M2i.%04x sn %05d\n", lCardIdx, astCard[lCardIdx].lCardType & TYP_VERSIONMASK, astCard[lCardIdx].lSerialNumber);

            int16 nSelection = _getch();
            if (nSelection == 27)
                return 1;
            if ((nSelection >= '0') && (nSelection < ('0' + lCardIdx)))
                lCardIdx = (nSelection - '0');
            }
        while (lCardIdx == lCardCount);

        // close all the other cards allowing a second instance of the program to run
        for (int32 lCloseIdx = 0; lCloseIdx < lCardCount; lCloseIdx++)
            if (lCloseIdx != lCardIdx)
                vSpcMCloseCard (&astCard[lCloseIdx]);
        }
    else
        lCardIdx = 0;



    // check whether we support this card type in the example
    if ((astCard[lCardIdx].eCardFunction != AnalogIn) && (astCard[lCardIdx].eCardFunction != DigitalIn) && (astCard[lCardIdx].eCardFunction != DigitalIO))
        return nSpcMErrorMessageStdOut (&astCard[lCardIdx], "Error: Card function not supported by this example\n", false);


    // we start with 16 bit acquisition as this is supported by all cards
    switch (astCard[lCardIdx].eCardFunction)
        {
        case AnalogIn:
            switch (astCard[lCardIdx].lBytesPerSample)
                {
                case 1: g_qwChannelEnable = CHANNEL1 | CHANNEL0; break;
                case 2: g_qwChannelEnable =            CHANNEL0; break;
                }
            break;

        case DigitalIn:
        case DigitalIO:
            g_qwChannelEnable = 0xffff;
            break;
        }



    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
    while (bSetup(&astCard[lCardIdx]))
        {
        if (!astCard[lCardIdx].bSetError)
            bDoCardSetup (&astCard[lCardIdx]);


        // ------------------------------------------------------------------------
        // setup the data transfer thread and start it, we use a timeout of 5 s in the example
        memset (&stBufferData, 0, sizeof(stBufferData));
        stBufferData.pstCard =      &astCard[lCardIdx];
        stBufferData.bStartCard =   true;
        stBufferData.bStartData =   true;
        stBufferData.lTimeout =     5000;

        // setup for async esc check
        g_nKeyPress = GetAsyncKeyState(VK_ESCAPE);

        // start the threaded version if g_bThread is defined
        if (!astCard[lCardIdx].bSetError && g_bThread)
            vDoThreadMainLoop (&stBufferData, &stWorkData, bWorkInit, bWorkDo, vWorkClose, bKeyCheckAsync);

        // start the unthreaded version with a smaller timeout of 100 ms to gain control about the FIFO loop
        stBufferData.lTimeout =     100;
        if (!astCard[lCardIdx].bSetError && !g_bThread)
            vDoMainLoop (&stBufferData, &stWorkData, bWorkInit, bWorkDo, vWorkClose, bKeyCheckAsync);


        // ------------------------------------------------------------------------
        // print error information if an error occured
        if (astCard[lCardIdx].bSetError)
            return nSpcMErrorMessageStdOut (&astCard[lCardIdx], "An error occured while programming the card:\n", true);

        } // if (bStart)


    // clean up and close the driver
    vSpcMCloseCard (&astCard[lCardIdx]);

    return EXIT_SUCCESS;
    }

