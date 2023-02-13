/*
**************************************************************************

rec_fifo_gate.cpp                                       (c) Spectrum GmbH

**************************************************************************

this example supports all acquisition cards with the option gated sampling
installed. If timestamp is installed the timestamps are also
read.

Does a continous FIFO transfer and examines the gates length. Data itself 
isnot touched.

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



// ----- standard c include files -----
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// ----- include of common example librarys -----
#include "../common/spcm_lib_card.h"
#include "../common/spcm_lib_data.h"
#include "../common/spcm_lib_thread.h"
#include "../sb5_file/sb5_file.h"


// ----- this is the global thread flag that defines whether we use the thread or non-thread loop -----
bool g_bThread = true;



/*
**************************************************************************
Working routine data
**************************************************************************
*/

struct ST_WORKDATA
    {
    int32       lPreSamples;
    int32       lPostSamples;
    int64       llGatesRead;
    };



/*
**************************************************************************
bDoCardSetup: setup matching the calculation routine
**************************************************************************
*/

bool bDoCardSetup (ST_WORKDATA* pstWorkData, ST_BUFFERDATA* pstBufferData)
    {
    int     i;

    pstWorkData->lPreSamples  = 32;
    pstWorkData->lPostSamples = 32;

    // we try to set the samplerate to 1 MHz on internal PLL, no clock output
    // we try to set the samplerate to 1 MHz (M2i) or 20 MHz (M3i) on internal PLL, no clock output
    // increase this to test the read-out-after-overrun
    if (pstBufferData->pstCard->bM2i)
        bSpcMSetupClockPLL (pstBufferData->pstCard, MEGA(1), false);
    else if (pstBufferData->pstCard->bM2p)
        bSpcMSetupClockPLL (pstBufferData->pstCard, MEGA(10), false);
    else if (pstBufferData->pstCard->bM3i || pstBufferData->pstCard->bM4i)
        bSpcMSetupClockPLL (pstBufferData->pstCard, MEGA(20), false);
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstBufferData->pstCard->llSetSamplerate / MEGA(1));

    // we set external trigger high gate for gated sampling
    bSpcMSetupTrigExternal (pstBufferData->pstCard, SPC_TM_POS, false, 0);

    // type dependent card setup
    switch (pstBufferData->pstCard->eCardFunction)
        {

        // analog acquisition card setup
        case AnalogIn:
            
            // we only enable 1 channel for the example
            bSpcMSetupModeRecFIFOGate (pstBufferData->pstCard, CHANNEL0, pstWorkData->lPreSamples, pstWorkData->lPostSamples, 0);

            // program all input channels to +/-1 V and 50 ohm termination (if it's available)
            for (i=0; i < pstBufferData->pstCard->lMaxChannels; i++)
                bSpcMSetupInputChannel (pstBufferData->pstCard, i, 1000, true);
            break;

        // digital acquisition card setup
        case DigitalIn:
        case DigitalIO:
            
            // we enable 16 channels
            bSpcMSetupModeRecFIFOGate (pstBufferData->pstCard, 0xffff, pstWorkData->lPreSamples, pstWorkData->lPostSamples, 0);

            // set all input channel groups to 110 ohm termination (if it's available)
            for (i=0; i < pstBufferData->pstCard->uCfg.stDIO.lGroups; i++)
                bSpcMSetupDigitalInput (pstBufferData->pstCard, i, true);
            break;
        }

    // set up the timestamp mode to standard if timestamp is installed
    if (pstBufferData->bStartTimestamp)
        bSpcMSetupTimestamp (pstBufferData->pstCard, SPC_TSMODE_STANDARD | SPC_TSCNT_INTERNAL, 0);

    return pstBufferData->pstCard->bSetError;
    }



/*
**************************************************************************
Setup working routine
**************************************************************************
*/

bool bWorkSetup (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA* pstWorkData = (ST_WORKDATA*) pvWorkData;

    // setup for the transfer, to avoid overrun one may use even larger data blocks
    pstBufferData->dwDataBufLen =   MEGA_B(8);
    pstBufferData->dwDataNotify =   KILO_B(512);
    pstBufferData->dwTSBufLen =     MEGA_B(1);
    pstBufferData->dwTSNotify =     KILO_B(4);

    // setup for the work
    pstWorkData->llGatesRead = 0;

    return true;
    }



/*
**************************************************************************
bWorkDo: get data and examine the gate intervals
**************************************************************************
*/

bool bWorkDo (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA*    pstWorkData = (ST_WORKDATA*) pvWorkData;

    // read the timestamps, count them and calculate the gate length
    uint64 qwTSAvail = pstBufferData->llTSAvailBytes;
    uint64 qwTimestampAverageDiff = 0;
    uint64 qwGates = 0;
    uint32 dwTimestampBytes;

    // ----- M2i+M3i timestamps are 64 bit, M4i/M4x and M2p 128 bit -----
    dwTimestampBytes = 0;
    if (pstBufferData->pstCard->bM2i || pstBufferData->pstCard->bM3i)
        dwTimestampBytes = 8;
    else if (pstBufferData->pstCard->bM4i || pstBufferData->pstCard->bM2p)
        dwTimestampBytes = 16;

    // we measure the timestamp diff (gate length) and calc an average
    while (qwTSAvail >= (2 * dwTimestampBytes))
        {
        uint64* pqwTimestamp1 = (uint64*)  pstBufferData->pvTSCurrentBuf;
        uint64* pqwTimestamp0 = ((uint64*) (((char*)pstBufferData->pvTSCurrentBuf) + dwTimestampBytes));
        uint64  qwDiff = (*pqwTimestamp1) - (*pqwTimestamp0);

        // the timestamps only record the "real world" trigger, we have to add pre- and post manually
        qwDiff += pstWorkData->lPreSamples + pstWorkData->lPostSamples;

        // next gate timestamp pair
        pstWorkData->llGatesRead++;
        qwTimestampAverageDiff += qwDiff;
        qwTSAvail -= (2 * dwTimestampBytes);
        qwGates++;
        pstBufferData->pvTSCurrentBuf  =    (void*) (((char*) pstBufferData->pvTSCurrentBuf) + (2 * dwTimestampBytes));
        }

    // subtract the bytes that we've not used so far as we still need the rest next time
    pstBufferData->llTSAvailBytes -=    qwTSAvail;    

    // this would also be the position to do something with the data
    // TO DO

    // print the status
    if (qwGates)
        {
        qwTimestampAverageDiff /= qwGates;
        printf ("\r%.2f MSamples transferred %d Gate Segments, Average Gate Len = %.3f ms\n", (double) pstBufferData->qwDataTransferred / pstBufferData->pstCard->lBytesPerSample / 1024.0 / 1024.0,  (uint32) pstWorkData->llGatesRead, 1000.0 * ((uint32) qwTimestampAverageDiff) / pstBufferData->pstCard->llSetSamplerate);
        }
    else
        printf ("\r%.2f MSamples transferred %d Gate Segments\n", (double) pstBufferData->qwDataTransferred / pstBufferData->pstCard->lBytesPerSample / 1024.0 / 1024.0,  (uint32) pstWorkData->llGatesRead);
    
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
    ST_WORKDATA         stWorkData;
    ST_BUFFERDATA       stBufferData;       // buffer and start information

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

    memset (&stBufferData, 0, sizeof(stBufferData));
    stBufferData.pstCard = &stCard;


    // check whether we support this card type in the example
    if ((stCard.eCardFunction != AnalogIn) && (stCard.eCardFunction != DigitalIn) && (stCard.eCardFunction != DigitalIO))
        return nSpcMErrorMessageStdOut (&stCard, "Error: Card function not supported by this example\n", false);
    if ((stCard.lFeatureMap & SPCM_FEAT_GATE) == 0)
        return nSpcMErrorMessageStdOut (&stCard, "Error: Gated Sampling Option not installed. Examples was done especially for this option!\n", false);

    // set a flag if timestamp is installed
    stBufferData.bStartTimestamp = ((stCard.lFeatureMap & SPCM_FEAT_TIMESTAMP) != 0);



    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
    if (!stCard.bSetError)
        bDoCardSetup (&stWorkData, &stBufferData);



    // ------------------------------------------------------------------------
    // setup the data transfer and start it
    stBufferData.bStartCard =       true;
    stBufferData.bStartData =       true;
    stBufferData.bStartTimestamp =  ((stCard.lFeatureMap & SPCM_FEAT_TIMESTAMP) != 0);

    // start the threaded loop
    stBufferData.lTimeout =         5000;
    if (!stCard.bSetError && g_bThread)
        vDoThreadMainLoop (&stBufferData, &stWorkData, bWorkSetup, bWorkDo, NULL, bKeyAbortCheck);

    // this is the non threaded loop, we use a small timeout of 100 ms here as we otherwise won't return from the loop if no dat is coming
    stBufferData.lTimeout =         100;
    if (!stCard.bSetError && !g_bThread)
        vDoMainLoop (&stBufferData, &stWorkData, bWorkSetup, bWorkDo, NULL, bKeyAbortCheck);



    // ------------------------------------------------------------------------
    // print error information if an error occured
    if (stCard.bSetError)
        return nSpcMErrorMessageStdOut (&stCard, "An error occured while programming the card:\n", true);


    // clean up and close the driver
    vSpcMCloseCard (&stCard);

    return EXIT_SUCCESS;
    }

