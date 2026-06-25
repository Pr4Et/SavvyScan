/*
**************************************************************************

rep_fifo_hd_speed.cpp                                    (c) Spectrum GmbH

**************************************************************************

this example supports all generator cards

Does FIFO replay from hard disk to test the maximum writing performance
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
#include <math.h>
#include <conio.h>

// ----- include of common example librarys -----
#include "../common/spcm_lib_card.h"
#include "../common/spcm_lib_data.h"
#include "../common/ostools/spcm_oswrap.h"
#include "../common/ostools/spcm_ostools.h"


// ----- global setup for the run (can be changed interactively) -----
int32   g_lSamplingRate =   MEGA(10);
int32   g_lNotifySize =     KILO_B(1024);
int32   g_lSWBufferSize =   MEGA_B(128);
int32   g_lHWBufferSize =   MEGA_B(64);
int64   g_llFileSize =      GIGA_B (4);
uint64  g_qwChannelEnable = 1;
HANDLE  g_hFile =           NULL;

#define FILENAME "stream_test.bin"



/*
**************************************************************************
bDoCardSetup: setup matching the calculation routine
**************************************************************************
*/

bool bDoCardSetup (ST_SPCM_CARDINFO *pstCard)
    {

    bSpcMSetupModeRepFIFOSingle (pstCard, g_qwChannelEnable);

    // we try to set the samplerate to 1 MHz on internal PLL, no clock output
    bSpcMSetupClockPLL (pstCard, g_lSamplingRate, false);
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstCard->llSetSamplerate / MEGA(1));

    // we set software trigger, no trigger output
    bSpcMSetupTrigSoftware (pstCard, false);

    return pstCard->bSetError;
    }



/*
**************************************************************************
DoDataCalculation: calculates the output data. Reads data from file
**************************************************************************
*/

bool bDoDataCalculation (ST_SPCM_CARDINFO *pstCard, void* pvBuffer, int64 llBytesToCalculate)
    {
    uint32 dwRead;
    ReadFile (g_hFile, pvBuffer, (uint32) llBytesToCalculate, &dwRead, NULL);
    if (dwRead != (uint32) llBytesToCalculate)
        printf ("\nFile Read Error (completed?)\n");
    return (dwRead == (uint32) llBytesToCalculate);
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
    char   szErrorText[ERRORTEXTLEN];
    bool   bStart = false;

    while (!bStart)
        {
        printf ("\n\n");
        printf ("Current Setup\n-------------\n");
        printf ("S ....... Sampling Rate:    %.2lf MS/s\n", (double) g_lSamplingRate / MEGA(1));
        printf ("B ....... SW Buffer Size:   %.2lf MByte\n", (double) g_lSWBufferSize / MEGA_B(1));
        printf ("H ....... HW Buffer Size:   %.2lf MByte\n", (double) g_lHWBufferSize / MEGA_B(1));
        printf ("N ....... Notify Size:      %d kByte\n", g_lNotifySize / KILO_B(1));
        printf ("C ....... Channel Enable:   %x\n", g_qwChannelEnable);
        printf ("F ....... File Size:        %.2lf GByte\n", (double) g_llFileSize / GIGA_B(1));
        printf ("Enter ... Start Test\n");
        printf ("Esc ..... Abort\n");

        spcm_dwSetParam_i64 (pstCard->hDrv, SPC_CHENABLE, g_qwChannelEnable);
        spcm_dwSetParam_i64 (pstCard->hDrv, SPC_SAMPLERATE, g_lSamplingRate);
        spcm_dwGetParam_i32 (pstCard->hDrv, SPC_CHCOUNT, &lChannels);
        if (spcm_dwGetErrorInfo_i32 (pstCard->hDrv, NULL, NULL, szErrorText) != ERR_OK)
            printf ("\nSetup Error:\n------------\n%s\n\n", szErrorText);
        else
            {
            if (pstCard->eCardFunction == AnalogOut)
                printf ("          Transfer Speed:   %.2lf MByte/s\n", (double)  g_lSamplingRate * lChannels * pstCard->lBytesPerSample / MEGA_B(1));
            else
                printf ("          Transfer Speed:   %.2lf MByte/s\n", (double)  g_lSamplingRate * lChannels / 8 / MEGA_B(1));
            }
        printf ("\n");


        switch (getch())
            {
            case 27: return false;
            case 13: bStart = true; break;

            case 's':
            case 'S':
                printf ("Sampling Rate (MS/s): ");
                scanf ("%lf", &dTmp);
                g_lSamplingRate = (int32) (dTmp * MEGA(1));
                break;

            case 'f':
            case 'F':
                printf ("File Size (GByte): ");
                scanf ("%lf", &dTmp);
                g_llFileSize = (int64) (dTmp * GIGA_B(1));

                // must be multiple of page size (4k)
                g_llFileSize = ((g_llFileSize >> 12) << 12);
                break;

            case 'b':
            case 'B':
                printf ("SW Buffer Size (MByte): ");
                scanf ("%lf", &dTmp);
                g_lSWBufferSize = (int32) (dTmp * MEGA_B(1));
                break;

            case 'h':
            case 'H':
                printf ("HW Buffer Size (MByte): ");
                scanf ("%lf", &dTmp);
                g_lHWBufferSize = (int32) (dTmp * MEGA_B(1));
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

    // start was called -> generate the file or open an existing one
    while (1)
        {
        g_hFile = CreateFile (FILENAME, 
            GENERIC_READ,  
            FILE_SHARE_READ | FILE_SHARE_WRITE, 
            NULL, 
            OPEN_EXISTING, 
            FILE_ATTRIBUTE_NORMAL | FILE_FLAG_NO_BUFFERING, 
            NULL);

        // file exists, check size
        if (g_hFile != INVALID_HANDLE_VALUE)
            {
            LARGE_INTEGER liFileSize;
            GetFileSizeEx (g_hFile, &liFileSize);

            if (liFileSize.QuadPart >= g_llFileSize)
                return true;

            // too small -> delete
            CloseHandle (g_hFile);
            DeleteFile (FILENAME);
            }

        // make a new file
        g_hFile = CreateFile (FILENAME, 
            GENERIC_WRITE, 
            FILE_SHARE_READ | FILE_SHARE_WRITE, 
            NULL, 
            CREATE_ALWAYS, 
            FILE_ATTRIBUTE_NORMAL | FILE_FLAG_NO_BUFFERING, 
            NULL);

        int64  llFileStep = MEGA_B(1);
        void*  pvData = pvAllocMemPageAligned ((uint32) llFileStep);
        int64  llIdx, llRemain;
        uint32 dwWritten;

        printf ("\nWriting file with output data...\n");
        memset (pvData, 0x5A, (uint32) llFileStep);
        for (llIdx = 0; llIdx < g_llFileSize; llIdx += llFileStep)
            {
            llRemain = (g_llFileSize - llIdx) > llFileStep ? llFileStep : g_llFileSize - llIdx;

            printf ("\rFileSize: %.3lf GByte of %.3lf GByte", (double) llIdx / GIGA_B(1), (double) g_llFileSize / GIGA_B(1));
            WriteFile (g_hFile, pvData, (uint32) llRemain, &dwWritten, NULL);
            if (dwWritten != (uint32) llRemain)
                {
                printf ("File Write Error!\n");
                CloseHandle (g_hFile);
                return false;
                }
            if (kbhit())
                if (getch() == 27)
                    {
                    CloseHandle (g_hFile);
                    return false;
                    }

            }
        vFreeMemPageAligned (pvData, (uint32) llFileStep);
        CloseHandle (g_hFile);
        }

    return false;
    }



/*
**************************************************************************
main 
**************************************************************************
*/

int main ()
    {
    char                szBuffer[1024];     // a character buffer for any messages
    ST_SPCM_CARDINFO    stCard;             // info structure of my card
    void*               pvBuffer = NULL;
    uint32              dwErr;


    // some example checks
    if (g_lSWBufferSize % g_lNotifySize)
        {
        printf ("In our example we can only handle sw buffers that are a whole numbered multiple of the notify size\n");
        return 1;
        }


    // ------------------------------------------------------------------------
    // init card number 0 (the first card in the system), get some information and print it
    if (bSpcMInitCardByIdx (&stCard, 0))
        {
        printf (pszSpcMPrintDocumentationLink (&stCard, szBuffer, sizeof (szBuffer)));
        printf (pszSpcMPrintCardInfo (&stCard, szBuffer, sizeof (szBuffer)));
        }
    else
        return nSpcMErrorMessageStdOut (&stCard, "Error: Could not open card\n", true);


    // check whether we support this card type in the example
    if ((stCard.eCardFunction != AnalogOut) && (stCard.eCardFunction != DigitalOut) && (stCard.eCardFunction != DigitalIO))
        return nSpcMErrorMessageStdOut (&stCard, "Error: Card function not supported by this example\n", false);


        // we start with 16 bit acquisition as this is supported by all cards
    switch (stCard.eCardFunction)
        {
        case AnalogOut:
            switch (stCard.lBytesPerSample)
                {
                case 1: g_qwChannelEnable = CHANNEL1 | CHANNEL0; break;
                case 2: g_qwChannelEnable =            CHANNEL0; break;
                }
            break;

        default:
            g_qwChannelEnable = 0xffff;
            break;
        }



    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
    while (bSetup(&stCard))
        {
        if (!stCard.bSetError)
            bDoCardSetup (&stCard);

        // ------------------------------------------------------------------------
        // allocate and setup the fifo buffer and fill it once with data
        pvBuffer = pvAllocMemPageAligned ((uint32) g_lSWBufferSize);
        if (!pvBuffer)
            return nSpcMErrorMessageStdOut (&stCard, "Memory allocation error\n", false);

        // starting with firmware version V9 we can program the hardware buffer size to reduce the latency
        if (stCard.lCtrlFwVersion >= 9)
            {
            spcm_dwSetParam_i64 (stCard.hDrv, SPC_DATA_OUTBUFSIZE, g_lHWBufferSize);
            spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_WRITESETUP);
            }

        spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, (uint32) g_lNotifySize, pvBuffer, 0, g_lSWBufferSize);


        // do initial calculation
        printf ("Fill SW-Buffer from file...\n");
        for (int32 lPos = 0; lPos < g_lSWBufferSize; lPos += g_lNotifySize)
            bDoDataCalculation (&stCard, (void*) ((int32) pvBuffer + lPos), g_lNotifySize);


        spcm_dwSetParam_i64 (stCard.hDrv, SPC_DATA_AVAIL_CARD_LEN, g_lSWBufferSize);

        // now buffer is full of data and we start the transfer (output is not started yet), timeout is 1 second
        spcm_dwSetParam_i32 (stCard.hDrv, SPC_TIMEOUT, 1000);
        dwErr = spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);




        // ------------------------------------------------------------------------
        // this is our main output and calculation loop
        int64 llTransferredBytes = 0;
        int64 llAvailUser;
        int64 llBufferFillPromille;
        int64 llUserPos;
        bool  bStarted = false;

        while (!dwErr)
            {
            spcm_dwGetParam_i64 (stCard.hDrv, SPC_DATA_AVAIL_USER_LEN, &llAvailUser);
            spcm_dwGetParam_i64 (stCard.hDrv, SPC_FILLSIZEPROMILLE, &llBufferFillPromille);
            if (llTransferredBytes > GIGA_B(1))
                printf ("\rSW-Buffer: %3.0f%%   HW-Buffer:%3d%%, Total Bytes so far: %6.2f GB    ", (float) 100.0 * (g_lSWBufferSize - llAvailUser) / g_lSWBufferSize, (uint32) llBufferFillPromille / 10, (float) llTransferredBytes / GIGA_B(1));
            else
                printf ("\rSW-Buffer: %3.0f%%   HW-Buffer:%3d%%, Total Bytes so far: %6.2f MB", (float) 100.0 * (g_lSWBufferSize - llAvailUser) / g_lSWBufferSize, (uint32) llBufferFillPromille / 10, (float) llTransferredBytes / MEGA_B(1));


            // we recalculate the amount of data that is free and set this part available for card again
            // inhere we only take pieces of notify size
            if (llAvailUser >= g_lNotifySize)
                {
                llTransferredBytes += g_lNotifySize;
                spcm_dwGetParam_i64 (stCard.hDrv, SPC_DATA_AVAIL_USER_POS, &llUserPos);
                if (bDoDataCalculation (&stCard, ((int8*) pvBuffer) + llUserPos, g_lNotifySize))
                    dwErr = spcm_dwSetParam_i64 (stCard.hDrv, SPC_DATA_AVAIL_CARD_LEN, g_lNotifySize);
                else
                    dwErr = ERR_ABORT;
                }

            // we start the output as soon as we have a sufficient amount of data on card 
            // inhere we start if the hardware buffer is completely full
            if (!bStarted && !dwErr && (llBufferFillPromille == 1000))
                {
                printf ("\nStart the output\n");
                dwErr = spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER);
                bStarted = true;
                }

            // wait for the next buffer to be free
            if (!dwErr)
                dwErr = spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_DATA_WAITDMA);


            // check for esc=abort
            if (!dwErr)
                if (bKbhit())
                    if (cGetch() == 27)
                        {
                        printf ("\nOutput stopped\n");
                        spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_STOP);
                        dwErr = ERR_ABORT;
                        }
            }

        // show runtime errors
        if (dwErr && !stCard.bSetError)
            printf ("\nEnd with Runtime Error Code:%d\n-> %s\n\n", dwErr, pszSpcMTranslateRuntimeError (dwErr, szBuffer));

        // ------------------------------------------------------------------------
        // print error information if an error occured
        if (stCard.bSetError)
            return nSpcMErrorMessageStdOut (&stCard, "An error occured while programming the card:\n", true);

        // stop the card 
        spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_RESET);
        if (g_hFile)
            CloseHandle (g_hFile);

        vFreeMemPageAligned (pvBuffer, (uint32) g_lSWBufferSize);
        }

    // clean up and close the driver
    vSpcMCloseCard (&stCard);

    return EXIT_SUCCESS;
    }

