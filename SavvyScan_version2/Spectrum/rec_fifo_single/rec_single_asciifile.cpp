/*
**************************************************************************

rec_single_ascii.cpp                                     (c) Spectrum GmbH

**************************************************************************

Does a continous FIFO transfer and writes data as hex to an ascii file

Change the global flag g_bThread to use the threaded version or the plain
and more simple loop.

Feel free to use this source for own projects and modify it in any kind

Documentation for the API as well as a detailed description of the hardware
can be found in the manual for each device which can be found on our website:
https://www.spectrum-instrumentation.com/en/downloads

Further information can be found online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/knowledge-base-overview

**************************************************************************
*/


// ----- include standard driver header from library -----
#include "../c_header/dlltyp.h"
#include "../c_header/regs.h"
#include "../c_header/spcerr.h"
#include "../c_header/spcm_drv.h"

// ----- include of common example librarys -----
#include "../common/spcm_lib_card.h"
#include "../common/spcm_lib_data.h"
#include "../common/spcm_lib_thread.h"

// ----- standard c include files -----
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


// ----- this is the global thread flag that defines whether we use the thread or non-thread loop -----
bool g_bThread = false;


/*
**************************************************************************
bDoCardSetup: setup matching the calculation routine
**************************************************************************
*/

bool bDoCardSetup (ST_SPCM_CARDINFO *pstCard)
    {
    int     i;
    int64 llChannelMask;

    // set mask for maximal channels
    if (pstCard->lMaxChannels >= 64)
        llChannelMask = -1; // -1 is all bits set to 1 = 0xffffffffffffffff
    else
        llChannelMask = ((int64) 1 << pstCard->lMaxChannels) - 1;

    // FIFO mode setup, we run continuously and have 32 samples of pre data before trigger event
    // all available channels are activated
    bSpcMSetupModeRecFIFOSingle (pstCard, llChannelMask, 32);

    // we try to set the samplerate to 1 MHz (M2i) or 20 MHz (M3i, M4i) on internal PLL, no clock output
    // increase this to test the read-out-after-overrun
	if (pstCard->bM2i)
        bSpcMSetupClockPLL (pstCard, MEGA(1), false);
	else if (pstCard->bM2p)
        bSpcMSetupClockPLL (pstCard, MEGA(10), false);
	else if (pstCard->bM3i || pstCard->bM4i)
        bSpcMSetupClockPLL (pstCard, MEGA(20), false);
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstCard->llSetSamplerate / 1000000);

    // we set software trigger, no trigger output
    bSpcMSetupTrigSoftware (pstCard, false);

    // type dependent card setup
    switch (pstCard->eCardFunction)
        {

        // analog acquisition card setup
        case AnalogIn:

            // program all input channels to +/-1 V and 50 ohm termination (if it's available)
            if (pstCard->bM2i || pstCard->bM2p)
                {
                for (i=0; i < pstCard->lMaxChannels; i++)
                    bSpcMSetupInputChannel (pstCard, i, 1000, true);
                }
            else
                {
                bool bTerm           = true;
                bool bACCoupling     = false;
                bool bBandwidthLimit = false;
                for (i=0; i < pstCard->lMaxChannels; i++)
                    bSpcMSetupPathInputCh (pstCard, i, 0, 1000, 0, bTerm, bACCoupling, bBandwidthLimit);
                }
            break;

        // digital acquisition card setup
        case DigitalIn:
        case DigitalIO:

            // set all input channel groups to 110 ohm termination (if it's available)
            for (i=0; i < pstCard->uCfg.stDIO.lGroups; i++)
                bSpcMSetupDigitalInput (pstCard, i, true);
            break;
        }

    return pstCard->bSetError;
    }



#define FILENAME "FIFO_ascii_test.txt"

/*
**************************************************************************
Working routine data
**************************************************************************
*/

struct ST_WORKDATA
    {
    FILE*   hFile;
    int64   llFileSize;
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
    pstBufferData->dwDataBufLen = MEGA_B(64);
    pstBufferData->dwDataNotify = MEGA_B(1);

    // setup for the work
    pstWorkData->hFile = fopen (FILENAME, "wt");
    pstWorkData->llFileSize = 0;

    return (pstWorkData->hFile != NULL);
    }



/*
**************************************************************************
bWorkDo: stores data to hard disk
**************************************************************************
*/

bool bWorkDo (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA* pstWorkData = (ST_WORKDATA*) pvWorkData;
    int32 lChannels =   pstBufferData->pstCard->lSetChannels;
    int32 lSamples =    pstBufferData->dwDataNotify / pstBufferData->pstCard->lBytesPerSample;
    FILE* f =           pstWorkData->hFile;
    int32 lIndex;

    
    // we only care for blocks of notify size
    if (pstBufferData->llDataAvailBytes < pstBufferData->dwDataNotify)
        return true;

    pstBufferData->llDataAvailBytes = pstBufferData->dwDataNotify;

    // we don't use the split data function as this takes too much time, Instead we go through data manually and write
    // as this is an universal example we need two loops for 8 bit wide analog data and for 12, 14, 16 bit analog data

    // 1 byte data (8 bit resolution)
    if (pstBufferData->pstCard->lBytesPerSample == 1)
        {
        int8*  pbyData = (int8*) pstBufferData->pvDataCurrentBuf;
        for (lIndex = 0; lIndex < lSamples; lIndex ++)
            {
            fprintf (f, "0x%02x ", *pbyData++ & 0xff);
            if ((lIndex % lChannels) == (lChannels -1))
                fprintf (f, "\n");
            }
        pstWorkData->llFileSize += 5 * lSamples + 2 * (lSamples / lChannels);
        }

    // 2 byte data (12, 14 and 16 bit resolution)
    else
        {
        int16* pnData = (int16*) pstBufferData->pvDataCurrentBuf;
        for (lIndex = 0; lIndex < lSamples; lIndex ++)
            {
            fprintf (f, "0x%04x ", *pnData++ & 0xffff);
            if ((lIndex % lChannels) == (lChannels -1))
                fprintf (f, "\n");
            }
        pstWorkData->llFileSize += 7 * lSamples + 2 * (lSamples / lChannels);
        }

    // announce the number of data that has been written
    printf ("\r%.2f MSamples (sum) written to %s (file size: %.2f MB)", 
        (double) pstBufferData->qwDataTransferred / pstBufferData->pstCard->lBytesPerSample / 1024.0 / 1024.0, 
        FILENAME,
        (double) pstWorkData->llFileSize / 1024 / 1024);


    return true;
    }




/*
**************************************************************************
vWorkClose: Close the work and clean up
**************************************************************************
*/

void vWorkClose (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA* pstWorkData = (ST_WORKDATA*) pvWorkData;

    if (pstWorkData->hFile)
        fclose (pstWorkData->hFile);
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
    ST_BUFFERDATA       stBufferData;       // buffer and transfer definitions
    ST_WORKDATA         stWorkData;         // work data for the working functions


    // ------------------------------------------------------------------------
    // init card number 0 (the first card in the system), get some information and print it
    // uncomment the second line and replace the IP address to use remote
    // cards like in a digitizerNETBOX
    if (bSpcMInitCardByIdx (&stCard, 0))
    //if (bSpcMInitCardByIdx (&stCard, "192.168.1.10", 0))
        {
        printf (pszSpcMPrintDocumentationLink (&stCard, szBuffer, sizeof (szBuffer)));
        printf (pszSpcMPrintCardInfo (&stCard, szBuffer, sizeof (szBuffer)));
        }
    else
        return nSpcMErrorMessageStdOut (&stCard, "Error: Could not open card\n", true);


    // check whether we support this card type in the example
    if ((stCard.eCardFunction != AnalogIn) && (stCard.eCardFunction != DigitalIn) && (stCard.eCardFunction != DigitalIO))
        return nSpcMErrorMessageStdOut (&stCard, "Error: Card function not supported by this example\n", false);


    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
    if (!stCard.bSetError)
        bDoCardSetup (&stCard);


    // ------------------------------------------------------------------------
    // setup the data transfer thread and start it, we use atimeout of 5 s in the example
    memset (&stBufferData, 0, sizeof(stBufferData));
    stBufferData.pstCard =      &stCard;
    stBufferData.bStartCard =   true;
    stBufferData.bStartData =   true;
    stBufferData.lTimeout =     5000;

    // start the threaded version if g_bThread is defined
    if (!stCard.bSetError && g_bThread)
        vDoThreadMainLoop (&stBufferData, &stWorkData, bWorkInit, bWorkDo, vWorkClose, bKeyAbortCheck);

    // start the unthreaded version with a smaller timeout of 100 ms to gain control about the FIFO loop
    stBufferData.lTimeout =     100;
    if (!stCard.bSetError && !g_bThread)
        vDoMainLoop (&stBufferData, &stWorkData, bWorkInit, bWorkDo, vWorkClose, bKeyAbortCheck);


    // ------------------------------------------------------------------------
    // print error information if an error occured
    if (stCard.bSetError)
        return nSpcMErrorMessageStdOut (&stCard, "An error occured while programming the card:\n", true);


    // clean up and close the driver
    vSpcMCloseCard (&stCard);

    return EXIT_SUCCESS;
    }

