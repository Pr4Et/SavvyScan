/*
**************************************************************************

rec_fifo_aba_thread.cpp                                  (c) Spectrum GmbH

**************************************************************************

this example supports all acquisition cards with the Multiple Recording,
ABA and Timestamp option installed. 

Does a continous FIFO transfer in Multiple Recording mode and uses three
threads to continuously read all three kinds of data:
1) segment data is read and transferred data is simply counted
2) ABA data is continuously read and averaged per block
3) timestamps are continuously read and difference is calculated

The main thread gets the current information, checks for an abort and
displays the information of the three working threads

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




/*
**************************************************************************
Working routine data
**************************************************************************
*/

struct ST_WORKDATA
    {
    SPCM_MUTEX_HANDLE   hMutexAccess;       // mutex to save the access to the shared display data

    bool                bFirstDisplay;      // first display generates headlines
    int32               lSegmentsize;       // segment size of Multiple Recording

    int64               llDataTotal;        // total bytes that data thread has performed

    int64               llTSTotal;          // total bytes that timestamp thread has performed
    double              dTSDistance;        // timestamp distance (frequency of the signal)

    int64               llABATotal;         // total bytes that ABA thread has performed
    int16*              ppnChData[SPCM_MAX_AICHANNEL]; // channel data for splitted data
    double              dABAAverage;        // current ABA average value
    };



/*
**************************************************************************
bDoCardSetup: setup matching the calculation routine
**************************************************************************
*/

bool bDoCardSetup (ST_WORKDATA* pstWorkData, ST_BUFFERDATA* pstBufferData)
    {
    int     i;

    // we first setup the working data
    memset (pstWorkData, 0, sizeof(ST_WORKDATA));
    spcm_bCreateMutex (&pstWorkData->hMutexAccess);
    pstWorkData->lSegmentsize = KILO_B(4);              // segment size
    pstWorkData->llDataTotal =  0;
    pstWorkData->llTSTotal =    0;
    pstWorkData->llABATotal =   0;
    pstWorkData->dTSDistance =  0.0;
    pstWorkData->dABAAverage =  0.0;
    pstWorkData->bFirstDisplay =true;


    // define the buffer sizes and the notify size that we want
    pstBufferData->dwTSBufLen =     MEGA_B(1);
    pstBufferData->dwTSNotify =     KILO_B(4);
    pstBufferData->dwDataBufLen =   MEGA_B(8);
    pstBufferData->dwDataNotify =   pstWorkData->lSegmentsize * pstBufferData->pstCard->lSetChannels * pstBufferData->pstCard->lBytesPerSample;
    pstBufferData->dwABABufLen =    MEGA_B(1);
    pstBufferData->dwABANotify =    KILO_B(4);


    // allocate memory for our sorted channel date for the ABA averaging
    for (i=0; i<pstBufferData->pstCard->lSetChannels; i++)
        pstWorkData->ppnChData[i] = new int16[pstBufferData->dwABANotify / pstBufferData->pstCard->lSetChannels];


    // we try to set the samplerate to 1 MHz (M2i) or 20 MHz (M3i) on internal PLL, no clock output
    // increase this to test the read-out-after-overrun
	if (pstBufferData->pstCard->bM2i)
        bSpcMSetupClockPLL (pstBufferData->pstCard, MEGA(1), false);
	else if (pstBufferData->pstCard->bM2p)
        bSpcMSetupClockPLL (pstBufferData->pstCard, MEGA(10), false);
	else if (pstBufferData->pstCard->bM3i || pstBufferData->pstCard->bM4i)
        bSpcMSetupClockPLL (pstBufferData->pstCard, MEGA(20), false);
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstBufferData->pstCard->llSetSamplerate / MEGA(1));

    // we only enable 1 channel for the example and set the ABA divider to 128
    bSpcMSetupModeRecFIFOABA (pstBufferData->pstCard, CHANNEL0, pstWorkData->lSegmentsize, pstWorkData->lSegmentsize - 128, 128);


    // we set external trigger for multiple recording
    if (pstBufferData->pstCard->bM2i)
        bSpcMSetupTrigExternal (pstBufferData->pstCard, SPC_TM_POS, false, 0);
    else // M3i, M4i, M4x, M2p
        bSpcMSetupTrigExternalLevel (pstBufferData->pstCard, SPC_TM_POS, 1500, 800, false);

    // type dependent card setup
    switch (pstBufferData->pstCard->eCardFunction)
        {

        // analog acquisition card setup
        case AnalogIn:

            // program all input channels to +/-1 V and 50 ohm termination (if it's available)
            if (pstBufferData->pstCard->bM2i || pstBufferData->pstCard->bM2p)
                {
                for (i=0; i < pstBufferData->pstCard->lMaxChannels; i++)
                    bSpcMSetupInputChannel (pstBufferData->pstCard, i, 1000, true);
                }
            else
                {
                bool bTerm           = true;
                bool bACCoupling     = false;
                bool bBandwidthLimit = false;
                for (i=0; i < pstBufferData->pstCard->lMaxChannels; i++)
                    bSpcMSetupPathInputCh (pstBufferData->pstCard, i, 0, 1000, 0, bTerm, bACCoupling, bBandwidthLimit);
                }
            break;

        // digital acquisition card setup
        case DigitalIn:
        case DigitalIO:
            printf ("Not yet implemented\n");
            return false;
            break;
        }

    // set up the timestamp mode to standard if timestamp is installed
    bSpcMSetupTimestamp (pstBufferData->pstCard, SPC_TSMODE_STANDARD | SPC_TSCNT_INTERNAL, 0);

    return pstBufferData->pstCard->bSetError;
    }



/*
**************************************************************************
bWorkDo: cares for the sample data: just counts them in our example
**************************************************************************
*/

bool bWorkDo (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA*    pstWorkData = (ST_WORKDATA*) pvWorkData;

    // count the data that has been transferred and set the calculated average
    spcm_vGetMutex (&pstWorkData->hMutexAccess);
    pstWorkData->llDataTotal += pstBufferData->llDataAvailBytes;
    spcm_vReleaseMutex (&pstWorkData->hMutexAccess);

    return true;
    }



/*
**************************************************************************
bTSWorkDo: as we get several timestamps at once we add the differences and 
           calc an average of these
**************************************************************************
*/

bool bTSWorkDo (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA*    pstWorkData = (ST_WORKDATA*) pvWorkData;
    double          dDiff = 0.0;
    uint64*         pqwTimestamps = (uint64*) pstBufferData->pvTSCurrentBuf;
    uint32          dwTimestampSize = 1; // size of timestamp in uint64
    int32           lOversampling = 1;

    // M2i and M3i work with 1 x 64 bit timestamps and oversampling while M4i works with 2 x 64 bit and without oversampling
    if (pstBufferData->pstCard->bM2i || pstBufferData->pstCard->bM3i)
        {
        dwTimestampSize = 1;
        lOversampling = pstBufferData->pstCard->lOversampling;
        }
    else if (pstBufferData->pstCard->bM4i || pstBufferData->pstCard->bM2p)
        {
        dwTimestampSize = 2;
        lOversampling = 1;
        }

    // work only with complete timestamps = 8 bytes (the possible half one
    // is automatically used at the next time the function is called)
    pstBufferData->llTSAvailBytes &= ~0x07;

    // each timestamp is 64/128 bit = 8/16 bytes, depending on card series (see above)
    if ((pstBufferData->llTSAvailBytes / (dwTimestampSize * 8)) > 1)
        {
        for (uint32 i=dwTimestampSize; i < (pstBufferData->llTSAvailBytes / (dwTimestampSize * (dwTimestampSize * 8))); i++)
            dDiff += (double) ((int64) (pqwTimestamps[i] - pqwTimestamps[i - dwTimestampSize]));
        dDiff /= pstBufferData->pstCard->llSetSamplerate / lOversampling;
        dDiff /= (pstBufferData->llTSAvailBytes / (dwTimestampSize * 8) - 1);
        }

    // count the data that has been transferred and set the calculated average
    spcm_vGetMutex (&pstWorkData->hMutexAccess);
    pstWorkData->llTSTotal += pstBufferData->llTSAvailBytes;
    pstWorkData->dTSDistance = dDiff;
    spcm_vReleaseMutex (&pstWorkData->hMutexAccess);

    return true;
    }



/*
**************************************************************************
bABAWorkDo: we calc the average of channel 0. 
**************************************************************************
*/

bool bABAWorkDo (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA*    pstWorkData = (ST_WORKDATA*) pvWorkData;
    double          dAverage = 0.0;

    // we inhere limit the length if data to proceed to our defined notify size as we allocated the channel data buffer in this size
    if (pstBufferData->llABAAvailBytes > pstBufferData->dwABANotify)
        pstBufferData->llABAAvailBytes = pstBufferData->dwABANotify;

    // split the data and calculate the average of channel 0
    if (pstBufferData->pstCard->lBytesPerSample == 1)
        bSpcMDemuxAnalogData (pstBufferData->pstCard, (int8*)pstBufferData->pvABACurrentBuf, pstBufferData->llABAAvailBytes / pstBufferData->pstCard->lBytesPerSample, pstWorkData->ppnChData);
    else
        bSpcMDemuxAnalogData (pstBufferData->pstCard, (int16*)pstBufferData->pvABACurrentBuf, pstBufferData->llABAAvailBytes / pstBufferData->pstCard->lBytesPerSample, pstWorkData->ppnChData);
    dAverage = dSpcMIntToVoltage (pstBufferData->pstCard, 0, dSpcMCalcAverage (pstWorkData->ppnChData[0], pstBufferData->llABAAvailBytes / pstBufferData->pstCard->lBytesPerSample / pstBufferData->pstCard->lSetChannels));

    // count the data that has been transferred and set the calculated average
    spcm_vGetMutex (&pstWorkData->hMutexAccess);
    pstWorkData->llABATotal += pstBufferData->llABAAvailBytes;
    pstWorkData->dABAAverage = dAverage;
    spcm_vReleaseMutex (&pstWorkData->hMutexAccess);

    return true;
    }



/*
**************************************************************************
vABAClean: we delete our channel data object
**************************************************************************
*/

void vABAClean (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA*    pstWorkData = (ST_WORKDATA*) pvWorkData;
    int             i;

    for (i=0; i<pstBufferData->pstCard->lSetChannels; i++)
        delete (pstWorkData->ppnChData[i]);
    }


/*
**************************************************************************
Abort function: is called on every update, checks for an abort and display
the results of all threads
**************************************************************************
*/

bool bAbortCheckAndDisplay (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA*    pstWorkData = (ST_WORKDATA*) pvWorkData;

    // print the headline if it's the first display
    if (pstWorkData->bFirstDisplay)
        {
        pstWorkData->bFirstDisplay = false;

        printf ("\n\n");
        printf ("******************************************************************\n");
        printf ("Data: transferred sample data in MByte\n");
        printf ("ABA:  transferred slow ABA data in MByte (current average of ch0)\n");
        printf ("TS:   transferred timestamps in MByte (average distance in block)\n");
        printf ("******************************************************************\n");
        }

    // do the display, secured by mutex as it's shared data between the three working threads and the main thread
    spcm_vGetMutex (&pstWorkData->hMutexAccess);
    printf ("\rData:%.3fM   ABA:%.3fM (Av:%.3f mV)   TS: %.3fM (%.3f ms)", 
        (double) pstWorkData->llDataTotal / MEGA_B(1),
        (double) pstWorkData->llABATotal / MEGA_B(1),
        pstWorkData->dABAAverage * KILO(1),
        (double) pstWorkData->llTSTotal / MEGA_B(1),
        pstWorkData->dTSDistance * KILO(1));
    spcm_vReleaseMutex (&pstWorkData->hMutexAccess);

    // check for abort
    if (bKbhit())
        if (cGetch() == 27)
            return true;

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
    if ((stCard.lFeatureMap & SPCM_FEAT_MULTI) == 0)
        return nSpcMErrorMessageStdOut (&stCard, "Error: Multiple Recording Option not installed. This example needs the option\n", false);
    if ((stCard.lFeatureMap & SPCM_FEAT_ABA) == 0)
        return nSpcMErrorMessageStdOut (&stCard, "Error: ABA Option not installed. This example needs the option\n", false);
    if ((stCard.lFeatureMap & SPCM_FEAT_TIMESTAMP) == 0)
        return nSpcMErrorMessageStdOut (&stCard, "Error: Timestamp Option not installed. This example needs the option\n", false);


    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
    if (!stCard.bSetError)
        bDoCardSetup (&stWorkData, &stBufferData);



    // ------------------------------------------------------------------------
    // some additional information on the acquisition
    if (!stCard.bSetError)
        {
        printf ("\nData information:\n=================\n");
        printf ("Each segment is %.3f ms long\n", 1000.0 * stWorkData.lSegmentsize / stCard.llSetSamplerate);
        printf ("Maximum pulse repetition frequency to reach with this setting is %.2f Hz\n", (double) stCard.llSetSamplerate / stWorkData.lSegmentsize);
        }



    // ------------------------------------------------------------------------
    // setup the data transfer and start it
    stBufferData.bStartCard =       true;
    stBufferData.bStartData =       true;
    stBufferData.bStartABA =        true;
    stBufferData.bStartTimestamp =  true;

    // start the threaded loop (10.000 ms = 10 s timeout)
    stBufferData.lTimeout =         KILO(10);
    if (!stCard.bSetError)
        vDoThreadMainLoop (
            &stBufferData,          // related buffer data
            &stWorkData,            // related work data
            NULL,                   // data thread setup routine
            bWorkDo,                // data thread work routine
            NULL,                   // data thread close routine
            bAbortCheckAndDisplay,  // abort routine (and display)
            true,                   // use seperate threads for all transfers
            NULL,                   // timestamp setup routine
            bTSWorkDo,              // timestamp work routine
            NULL,                   // timestamp close routine
            NULL,                   // ABA setup routine
            bABAWorkDo,             // ABA work routine
            vABAClean);             // AVA close routine


    // ------------------------------------------------------------------------
    // print error information if an error occured
    if (stCard.bSetError)
        return nSpcMErrorMessageStdOut (&stCard, "An error occured while programming the card:\n", true);


    // clean up and close the driver
    vSpcMCloseCard (&stCard);

    return EXIT_SUCCESS;
    }

