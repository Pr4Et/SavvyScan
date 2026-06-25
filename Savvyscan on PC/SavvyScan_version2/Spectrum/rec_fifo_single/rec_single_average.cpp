/*
**************************************************************************

fifo_single_average.cpp                                  (c) Spectrum GmbH

**************************************************************************

Does a simple averaging of continuous FIFO data and shows the result
in one print line

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
bool g_bThread = true;


/*
**************************************************************************
bDoCardSetuo: setup matching the calculation routine
**************************************************************************
*/

bool bDoCardSetup (ST_SPCM_CARDINFO *pstCard)
    {
    int     i;

    // FIFO mode setup, we run continuously and have 32 samples of pre data before trigger event
    // all available channels are activated
    bSpcMSetupModeRecFIFOSingle (pstCard, (1 << pstCard->lMaxChannels) - 1, 32);

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
        }

    return pstCard->bSetError;
    }



/*
**************************************************************************
Working routine data
**************************************************************************
*/

struct ST_WORKDATA
    {
    bool    bFirst;
    double  dMinInBlock;
    double  dMaxInBlock;
    int16*  ppnChannelData[SPCM_MAX_AICHANNEL];
    int32   lNotifySamplesPerChannel;
    };



/*
**************************************************************************
Setup working routine
**************************************************************************
*/

bool bWorkInit (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    bool bOk;

    ST_WORKDATA* pstWorkData = (ST_WORKDATA*) pvWorkData;

    // setup for the transfer
    pstBufferData->dwDataBufLen = MEGA_B(32);
    pstBufferData->dwDataNotify = KILO_B(64);

    // setup for the work
    pstWorkData->bFirst = true;
    pstWorkData->lNotifySamplesPerChannel = pstBufferData->dwDataNotify / pstBufferData->pstCard->lSetChannels / pstBufferData->pstCard->lBytesPerSample;

    // allocate arrays for channel data of notify size
    bOk = true;
    for (int i = 0; i < pstBufferData->pstCard->lSetChannels; i++)
        {
        pstWorkData->ppnChannelData[i] = (int16*) pvAllocMemPageAligned (pstWorkData->lNotifySamplesPerChannel * 2);
        if (!pstWorkData->ppnChannelData[i])
            bOk = false;
        }

    return bOk;
    }



/*
**************************************************************************
bWorkDo: does the averageing and prints the results in one line

This function is absolutely not optimized for high data throughput. It's
just an example to show the handling of data. To have highest trhoughput
one should not use the easy DemuxAnalogToInt16 function but access the
multiplexed data directly.
**************************************************************************
*/

bool bWorkDo (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA* pstWorkData = (ST_WORKDATA*) pvWorkData;
    double       dAverage;

    // we only care for blocks of notify size
    if (pstBufferData->llDataAvailBytes < pstBufferData->dwDataNotify)
        return true;

    pstBufferData->llDataAvailBytes = pstBufferData->dwDataNotify;

    // now let's split up the data
    if (pstBufferData->pstCard->lBytesPerSample == 1)
        bSpcMDemuxAnalogData (pstBufferData->pstCard, static_cast < int8* > (pstBufferData->pvDataCurrentBuf), pstWorkData->lNotifySamplesPerChannel, pstWorkData->ppnChannelData);
    else
        bSpcMDemuxAnalogData (pstBufferData->pstCard, static_cast < int16* > (pstBufferData->pvDataCurrentBuf), pstWorkData->lNotifySamplesPerChannel, pstWorkData->ppnChannelData);
        

    // calculate average of first channel
    dAverage =  dSpcMCalcAverage (pstWorkData->ppnChannelData[0], pstWorkData->lNotifySamplesPerChannel);
    dAverage =  dSpcMIntToVoltage (pstBufferData->pstCard, 0, dAverage);

    // store min/max if changed
    if (pstWorkData->bFirst || (pstWorkData->dMinInBlock > dAverage))
        pstWorkData->dMinInBlock = dAverage;
    if (pstWorkData->bFirst || (pstWorkData->dMaxInBlock < dAverage))
        pstWorkData->dMaxInBlock = dAverage;

    if (pstWorkData->bFirst)
        printf ("\n%11s %12s %12s %12s\n", "Transferred", "Average", "Min", "Max"); 

    pstWorkData->bFirst = false;

    // print some details 
    printf ("\r%8.2f MB %9.2f mV %9.2f mV %9.2f mV", 
        (double) pstBufferData->qwDataTransferred / 1024 / 1024,
        1000.0 * dAverage,
        1000.0 * pstWorkData->dMinInBlock,
        1000.0 * pstWorkData->dMaxInBlock);

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

    for (int i = 0; i < pstBufferData->pstCard->lSetChannels; i++)
        if (pstWorkData->ppnChannelData[i])
            vFreeMemPageAligned ((void*) pstWorkData->ppnChannelData[i], pstWorkData->lNotifySamplesPerChannel * 2);
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
    if (stCard.eCardFunction != AnalogIn)
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

