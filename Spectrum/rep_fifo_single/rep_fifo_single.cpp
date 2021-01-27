/*
**************************************************************************

rep_fifo_single.cpp                                       (c) Spectrum GmbH

**************************************************************************

Example for all M2i, M4i, M4x, M2p analog and digital generator cards. 
Shows FIFO replay mode as single shot.

To test output latency the example allows to change the output signal
frequency by pressing the space key.
  
Feel free to use this source for own projects and modify it in any kind.

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



// ----- include of common example librarys -----
#include "../common/spcm_lib_card.h"
#include "../common/spcm_lib_data.h"
#include "../common/ostools/spcm_oswrap.h"
#include "../common/ostools/spcm_ostools.h"
#include "../sb5_file/sb5_file.h"



/*
**************************************************************************
bDoCardSetup: setup matching the calculation routine
**************************************************************************
*/

bool bDoCardSetup (ST_SPCM_CARDINFO *pstCard)
    {
    int i;
    int lFilter = 0;

    // we try to set the samplerate to 1 MHz (M2i) or 50 MHz (M4i) on internal PLL, no clock output
    if (pstCard->bM4i)
        {
        bSpcMSetupClockPLL (pstCard, MEGA(50), false);
        lFilter = 1; // the only available filter
        }
    else
        {
        bSpcMSetupClockPLL (pstCard, MEGA(1), false);
        lFilter = 3; // highest cut-off frequency
        }
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstCard->llSetSamplerate / MEGA(1));

    // we set software trigger, no trigger output
    bSpcMSetupTrigSoftware (pstCard, false);

    // type dependent card setup
    switch (pstCard->eCardFunction)
        {

        // analog generator card setup
        case AnalogOut:

            // program all output channels to +/-1 V, zero offset and filter
            for (i=0; i < pstCard->lMaxChannels; i++)
                bSpcMSetupAnalogOutputChannel (pstCard, i, 1000, 0, lFilter);

            // FIFO mode setup, we run continuously 
            // only one chanel is activated for analog output to keep example simple
            bSpcMSetupModeRepFIFOSingle (pstCard, CHANNEL0);
            break;



        // digital generator card setup
        case DigitalIn:
        case DigitalIO:

            // FIFO mode setup, we run continuously 
            // 8 channels (1 byte) is activated
            bSpcMSetupModeRepFIFOSingle (pstCard, 0xff);
            break;
        }

    return pstCard->bSetError;
    }



/*
**************************************************************************
DoDataCalculation: calculates the output data. The calculation routine
                   is quite simple as we either have one analog channel
                   or 8 digital channels to calculate
**************************************************************************
*/
static int64 g_llOffset =       0;
static int64 g_llXDiv =         KILO_B(100);

bool bDoDataCalculation (ST_SPCM_CARDINFO *pstCard, void* pvBuffer, int64 llBytesToCalculate)
    {
    int64 i;
    int16* pnData = (int16*) pvBuffer;
    int8*  pbyData = (int8*) pvBuffer;
    double dSineXScale = 2.0 * 3.14159 / g_llXDiv;

    switch (pstCard->eCardFunction)
        {
        // analog generator card setup: 1 channel slow sine signal
        case AnalogOut:
            {
            if (pstCard->lBytesPerSample == 1)
                {
                for (i = 0; i < llBytesToCalculate; i++)
                    pbyData[i] = (int8) (pstCard->uCfg.stAO.lMaxDACValue * sin (dSineXScale * (g_llOffset + i)));
                }
            else
                {
                for (i = 0; i < llBytesToCalculate/2; i++)
                    pnData[i] = (int16) (pstCard->uCfg.stAO.lMaxDACValue * sin (dSineXScale * (g_llOffset + i)));
                }

            g_llOffset += (llBytesToCalculate / pstCard->lBytesPerSample);
            break;
            }


        // digital generator card setup: simple ramp
        case DigitalIn:
        case DigitalIO:
            {
            for (i=0; i<llBytesToCalculate; i++)
                pbyData[i] = (int8) i;
            break;
            }
        }

    return true;
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
    int                 nKeyCheck = 0;      // key check counter to avoid to much key polling

    // setup for the FIFO mode (HW buffer size can be programmed starting with firmware V9)
    int64        llHWBufSize =      KILO_B(64);
    int64        llSWBufSize =      KILO_B(128);
    int64        llNotifySize =     KILO_B(8);



    // some example checks
    if (llSWBufSize % llNotifySize)
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




    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
    if (!stCard.bSetError)
        {
        bool bError = bDoCardSetup (&stCard);
        if (bError)
            {
            printf ("Error: An error occured in card setup\n");
            return 1;
            }
        }


    // starting with firmware version V9 we can program the hardware buffer size to reduce the latency
    if (stCard.lCtrlFwVersion >= 9)
        {
        spcm_dwSetParam_i64 (stCard.hDrv, SPC_DATA_OUTBUFSIZE, llHWBufSize);
        spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_WRITESETUP);
        }

    // ------------------------------------------------------------------------
    // allocate and setup the fifo buffer and fill it once with data
    pvBuffer = pvAllocMemPageAligned ((uint32) llSWBufSize);
    if (!pvBuffer)
        return nSpcMErrorMessageStdOut (&stCard, "Memory allocation error\n", false);
    spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, (uint32) llNotifySize, pvBuffer, 0, llSWBufSize);
    bDoDataCalculation (&stCard, pvBuffer, llSWBufSize);
    spcm_dwSetParam_i64 (stCard.hDrv, SPC_DATA_AVAIL_CARD_LEN, llSWBufSize);

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
        printf ("\rSW-Buffer: %3.0f%%   HW-Buffer:%3d%%, Total Bytes so far: %6.2f MB", (float) 100.0 * (llSWBufSize - llAvailUser) / llSWBufSize, (uint32) llBufferFillPromille / 10, (float) llTransferredBytes / MEGA_B(1));


        // we recalculate the amount of data that is free and set this part available for card again
        // inhere we only take pieces of notify size
        if (llAvailUser >= llNotifySize)
            {
            llTransferredBytes += llNotifySize;
            spcm_dwGetParam_i64 (stCard.hDrv, SPC_DATA_AVAIL_USER_POS, &llUserPos);
            bDoDataCalculation (&stCard, ((int8*) pvBuffer) + llUserPos, llNotifySize);
            dwErr = spcm_dwSetParam_i64 (stCard.hDrv, SPC_DATA_AVAIL_CARD_LEN, llNotifySize);
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
            if (nKeyCheck-- == 0)
                {
                nKeyCheck = 64;
                if (bKbhit())
                    switch (cGetch())
                        {
                        case 27:
                            printf ("\nOutput stopped\n");
                            spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_STOP);
                            dwErr = ERR_ABORT;
                            break;

                        // blank changes the signal to test the output latency
                        case ' ':
                            printf ("\nChange Signal\n");
                            g_llXDiv /= 2;
                            break;

                        }
                }
        }

    // show runtime errors
    if (dwErr && !stCard.bSetError)
        printf ("\nEnd with Runtime Error Code:%d\n-> %s\n\n", dwErr, pszSpcMTranslateRuntimeError (dwErr, szBuffer));

    // ------------------------------------------------------------------------
    // print error information if an error occured
    if (stCard.bSetError)
        return nSpcMErrorMessageStdOut (&stCard, "An error occured while programming the card:\n", true);

    vFreeMemPageAligned (pvBuffer, (uint32) llSWBufSize);

    // clean up and close the driver
    vSpcMCloseCard (&stCard);

    return EXIT_SUCCESS;
    }

