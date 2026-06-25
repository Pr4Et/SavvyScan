/*
**************************************************************************

rec_std_gate.cpp                                         (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based acquisition cards with the option 
Gated Sampling and Timestamp installed. 

Information about the different products and their drivers can be found
online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/platform-driver-and-series-differences

Shows Gated Sampling data acquisition using single mode (one shot). Start
position and length of the gates are recorded using the timestamp feature.

This example is done quite simple and defines the number of gates that 
one expects globally. However if you have no idea how many gates are 
expected and you wish to read out all the timestamps it is necessary to
implement the timestamp readout dynamically. Please have a look at the
Gated Sampling FIFO example to see how this may look like.

**************************************************************************

Feel free to use this source for own projects and modify it in any kind
This example is provided "as is" without any guarantee.

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

#include "../common/ostools/spcm_oswrap.h"
#include "../common/ostools/spcm_ostools.h"

// ----- standard c include files -----
#include <stdio.h>
#include <stdlib.h>


// to keep it simple we define globally how many gate we will handle
// in this example. Each gate produces two timestamps
#define MAXGATES 1000


// maximum output lines to keep display simple
#define MAXLINES 20



/*
**************************************************************************
vDoCardSetup
**************************************************************************
*/

void vDoCardSetup (ST_SPCM_CARDINFO *pstCard, int64 llMemsize, int32 lPreSamples, int32 lPostSamples)
    {
    int   i;
    int64 llChannelMask;

    // gated sampling with 256k memory, each gate segment has 16 samples of pre data and 16 samples of post data
    if (pstCard->lMaxChannels >= 64)
        llChannelMask = -1; // -1 is all bits set to 1 = 0xffffffffffffffff
    else
        llChannelMask = ((int64) 1 << pstCard->lMaxChannels) - 1;

    bSpcMSetupModeRecStdGate (pstCard, llChannelMask, llMemsize, lPreSamples, lPostSamples);

    // we try to set the samplerate to 10 MHz on internal PLL, no clock output, 
    // if the card can't run 10 MHz it is set to maximum sampling rate internally
    bSpcMSetupClockPLL (pstCard, MEGA(10), false);
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstCard->llSetSamplerate / 1000000);

    // we set trigger to external high gate, please connect the trigger line!
    printf ("\n!!! Using external trigger - please connect a signal to the trigger input !!!\n\n");
    bSpcMSetupTrigExternal (pstCard, SPC_TM_POS, false, 0);

    // type dependent card setup
    switch (pstCard->eCardFunction)
        {

        // analog acquisition card setup
        case AnalogIn:

            // program all input channels to +/-1 V and 50 ohm termination (if it's available)
            for (i=0; i < pstCard->lMaxChannels; i++)
                bSpcMSetupInputChannel (pstCard, i, 1000, true);
            break;

        // digital acquisition card setup
        case DigitalIn:
        case DigitalIO:

            // set all input channel groups to 110 ohm termination (if it's available)
            for (i=0; i < pstCard->uCfg.stDIO.lGroups; i++)
                bSpcMSetupDigitalInput (pstCard, i, true);
            break;
        }

    // set up the timestamp mode to standard 
    bSpcMSetupTimestamp (pstCard, SPC_TSMODE_STARTRESET | SPC_TSCNT_INTERNAL, 0);
    }




/*
**************************************************************************
vShowAnalogData: displays the recorded data and timstamps. Does an average
over each gate segment and shows start position and length
**************************************************************************
*/

int16 nShowAnalogData (ST_SPCM_CARDINFO *pstCard, int32 lPreSamples, int32 lPostSamples, void* pvBuffer, uint64* pqwTimestamps, int64 llAvailTimestamps)
    {
    int i;
    int16*  ppnChannelData[SPCM_MAX_AICHANNEL];
    int32   lGateStartPos, lGateSegmentLen;
    int32   lGateIdx;
    int64   llGateCount;

    // allocate channel data
    for (i=0; i<pstCard->lSetChannels; i++)
        {
        ppnChannelData[i] = new int16[(int32) pstCard->llSetMemsize];
        if (!ppnChannelData[i])
            return nSpcMErrorMessageStdOut (pstCard, "Memory allocation error\n", false);
        }

    // split data function
    if (pstCard->lBytesPerSample == 1)
        bSpcMDemuxAnalogData (pstCard, (int8*)pvBuffer, (int32) pstCard->llSetMemsize, ppnChannelData);
    else
        bSpcMDemuxAnalogData (pstCard, (int16*)pvBuffer, (int32) pstCard->llSetMemsize, ppnChannelData);

    // some additional information on the acquisition
    printf ("\nData information:\n=================\n");
    printf ("We've read %d timestamps, making %d gate segments\n", (int)llAvailTimestamps, (int)(llAvailTimestamps / 2));

    // as the last gate is for sure not ending within the memory we only get one timestamp for this one
    llGateCount = ((llAvailTimestamps + 1) / 2);

    // loop across all acquired gate segments and print some information
    // keep in mind that each segment produces 2 timestamps
    lGateStartPos = 0;
    printf ("%6s %8s %8s %12s %15s %15s\n", "Index", "Pos", "Length", "Average Ch0", "Timestamp", "Timediff");
    for (lGateIdx = 0; (lGateIdx < llGateCount) && (lGateIdx < MAXLINES); lGateIdx ++)
        {

        // M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
        if (pstCard->bM2i || pstCard->bM3i)
            lGateSegmentLen = (int32) (pqwTimestamps[2 * lGateIdx + 1] - pqwTimestamps[2 * lGateIdx]);
        else if (pstCard->bM4i || pstCard->bM2p)
            lGateSegmentLen = (int32) (pqwTimestamps[4 * lGateIdx + 2] - pqwTimestamps[4 * lGateIdx]);

        // the timestamps only calculate the real world gate signal, we have to add our pre and postsample
        lGateSegmentLen += (lPreSamples + lPostSamples);

        // as the last gate didn't had a second timestamp we calc this length by the available data
        if (lGateIdx == (llGateCount - 1))
            lGateSegmentLen = (int32) (pstCard->llSetMemsize - lGateStartPos);

        printf ("%6d %8d %8d %9.2f mV", 
            lGateIdx, 
            lGateStartPos, 
            lGateSegmentLen, 
            1000.0 * dSpcMIntToVoltage (pstCard, 0, dSpcMCalcAverage (&ppnChannelData[0][lGateStartPos], lGateSegmentLen)));

        // M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
        if (pstCard->bM2i || pstCard->bM3i)
            printf (" %12.6f ms ", 1000.0 * (double) pqwTimestamps[2 * lGateIdx + 0] / pstCard->llSetSamplerate / pstCard->lOversampling);
        else if (pstCard->bM4i || pstCard->bM2p)
            printf (" %12.6f ms ", 1000.0 * (double) pqwTimestamps[4 * lGateIdx + 0] / pstCard->llSetSamplerate);

        // calc the difference to the last segment
        if ((lGateIdx > 0) && (pstCard->bM2i || pstCard->bM3i))
            printf ("%12.6f ms", 1000.0 * (double) (pqwTimestamps[2 * lGateIdx] - pqwTimestamps[2 * (lGateIdx - 1)]) / pstCard->llSetSamplerate / pstCard->lOversampling);
        else if ((lGateIdx > 0) && (pstCard->bM4i || pstCard->bM2p))
            printf ("%12.6f ms", 1000.0 * (double) (pqwTimestamps[4 * lGateIdx] - pqwTimestamps[4 * (lGateIdx - 1)]) / pstCard->llSetSamplerate);

        printf ("\n");
        lGateStartPos += lGateSegmentLen;
        }

    // clean up
    for (i=0; i<pstCard->lSetChannels; i++)
        delete [] (ppnChannelData[i]);

    return 0;
    }

/*
**************************************************************************
nShowDigitalData 
**************************************************************************
*/

int16 nShowDigitalData (ST_SPCM_CARDINFO *pstCard, void* pvBuffer, uint32 dwNrOfSamplesToShow)
    {
    int32 i, lChIdx;
    uint32 dwSampleIdx;
    int8* ppbyChannelData[SPCM_MAX_DIOCHANNEL];

    // allocate channel data
    for (i=0; i<pstCard->lSetChannels; i++)
        {
        ppbyChannelData[i] = new int8[(int32) pstCard->llSetMemsize];
        if (!ppbyChannelData[i])
            return nSpcMErrorMessageStdOut (pstCard, "Memory allocation error\n", false);
        }

    // split data function
    bSpcMDemuxDigitalDataToInt8 (pstCard, pvBuffer, (int32) (pstCard->llSetMemsize), ppbyChannelData);

    // print samples for each channel
    printf ("\nPrint first %u samples of all channels:\n\n", dwNrOfSamplesToShow);

    for (dwSampleIdx=0; dwSampleIdx<dwNrOfSamplesToShow; dwSampleIdx++)
        {
            printf ("SampleNr.%u\n", dwSampleIdx);

            if (pstCard->lSetChannels < 16)
                {

                // less than 16 channels
                lChIdx = pstCard->lSetChannels-1;

                while (lChIdx >= 0)
                    {
                    printf ("[D%d] = %d ", lChIdx, ppbyChannelData[lChIdx][dwSampleIdx]);
                    lChIdx--;
                    }

                printf ("\n\n");
                }
            else
                {

                // 16 channels and more
                lChIdx = pstCard->lSetChannels-1;

                while (lChIdx > 0)
                    {
                    printf ("[D%d ......... D%d] ", lChIdx, lChIdx - 15);
                    lChIdx -= 16;
                    }

                printf ("\n");

                lChIdx = pstCard->lSetChannels-1;
                while (lChIdx >= 0)
                    {
                    printf ("%d", ppbyChannelData[lChIdx][dwSampleIdx]);

                    if (!(lChIdx%16))
                        printf (" ");
                    else
                        if (!(lChIdx%4))
                            printf (".");

                    lChIdx--;
                    }

                printf ("\n\n");
                }	
        }

    for (i=0; i<pstCard->lSetChannels; i++)
        delete [] (ppbyChannelData[i]);

    return 0;
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
    int64               llMemsize;
    int32               lPreSamples, lPostSamples;
    int32               lOversampling = 1;
    int64               llAvailTimestamps, llAvailTimestampBytes;
    uint64              qwMemInBytes;
    uint64              qwTSBufferLen_bytes = 0;
    uint64*             pqwTimestamps;
    void*               pvBuffer;


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
    if ((stCard.lFeatureMap & SPCM_FEAT_GATE) == 0)
        return nSpcMErrorMessageStdOut (&stCard, "Error: Gated Sampling Option not installed. Examples was done especially for this option!\n", false);
    if ((stCard.lFeatureMap & SPCM_FEAT_TIMESTAMP) == 0)
        return nSpcMErrorMessageStdOut (&stCard, "Error: Timestamp Option not installed. Examples was done especially for this option!\n", false);


    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
    llMemsize =     KILO_B(64);
    lPreSamples =   32;
    lPostSamples =  32;
    if (!stCard.bSetError)
        vDoCardSetup (&stCard, llMemsize, lPreSamples, lPostSamples);

    // ------------------------------------------------------------------------
    // calculate the amount of data we need and allocate memory buffer
    if (!stCard.bSetError)
        {
        switch (stCard.eCardFunction)
            {
            case AnalogIn:
                qwMemInBytes = stCard.llSetMemsize * stCard.lBytesPerSample * stCard.lSetChannels;     
                break;

            case DigitalIn:
            case DigitalIO:
                qwMemInBytes = stCard.llSetMemsize;
                break;
            }


        // M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
        if (stCard.bM2i || stCard.bM3i)
            qwTSBufferLen_bytes = sizeof (uint64) * MAXGATES * 2;
        else if (stCard.bM4i || stCard.bM2p)
            qwTSBufferLen_bytes = 2ULL * sizeof (uint64) * MAXGATES * 2;

        pqwTimestamps = (uint64*)pvAllocMemPageAligned (qwTSBufferLen_bytes);

        pvBuffer = pvAllocMemPageAligned (qwMemInBytes);
        if (!pvBuffer || !pqwTimestamps)
            return nSpcMErrorMessageStdOut (&stCard, "Memory allocation error\n", false);
        }

    // ------------------------------------------------------------------------
    // make acquisition and get data
    if (!stCard.bSetError)
        {
        // We'll start and wait untill the card has finished or until a timeout occurs

        // M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
        if (stCard.bM2i || stCard.bM3i)
            spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_TIMESTAMP, SPCM_DIR_CARDTOPC, 0, (void*) pqwTimestamps, 0, MAXGATES * 2 * sizeof (uint64));
        else if (stCard.bM4i || stCard.bM2p)
            spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_TIMESTAMP, SPCM_DIR_CARDTOPC, 0, (void*) pqwTimestamps, 0, MAXGATES * 4 * sizeof (uint64));

        spcm_dwSetParam_i32 (stCard.hDrv, SPC_TIMEOUT, 5000);
        printf ("Starting the card together with timestamp transfer and waiting for ready interrupt\n");
        if (spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_START | M2CMD_EXTRA_STARTDMA | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY) == ERR_TIMEOUT)
            {
            vFreeMemPageAligned (pvBuffer, qwMemInBytes);
            return nSpcMErrorMessageStdOut (&stCard, "... Timeout", false);
            }
        else
            {

            // we define the buffer for transfer and start the DMA transfer
            printf ("Starting the DMA transfer and waiting until data is in PC memory\n");
            spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, 0, pvBuffer, 0, qwMemInBytes);
            spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);

            // check for error code
            if (spcm_dwGetErrorInfo_i32 (stCard.hDrv, NULL, NULL, szBuffer))
                {
                vFreeMemPageAligned (pvBuffer, qwMemInBytes);
                return nSpcMErrorMessageStdOut (&stCard, szBuffer, false);
                }
            printf ("... acquisition ended, data has been transferred to PC memory\n");

            // we read out timestamps (should be already done as we started timestamp read together with card
            spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_EXTRA_WAITDMA);
            spcm_dwGetParam_i64 (stCard.hDrv, SPC_TS_AVAIL_USER_LEN, &llAvailTimestampBytes);
            llAvailTimestamps = llAvailTimestampBytes / 8;
            printf ("... timestamps have been transferred to PC memory\n");
            }
        }



    // ------------------------------------------------------------------------
    // we do something with the acquired data
    if (!stCard.bSetError)
        switch (stCard.eCardFunction)
            {
            case AnalogIn:
                nShowAnalogData (&stCard, lPreSamples, lPostSamples, pvBuffer, pqwTimestamps, llAvailTimestamps);
                break;

            case DigitalIn:
            case DigitalIO:

                // show first 10 samples for each channel
                nShowDigitalData (&stCard, pvBuffer, 10);
                break;
            }
  


    // ------------------------------------------------------------------------
    // print error information if an error occured
    if (stCard.bSetError)
        return nSpcMErrorMessageStdOut (&stCard, "An error occured while programming the card:\n", true);


    // clean up and close the driver
    vSpcMCloseCard (&stCard);
    vFreeMemPageAligned (pvBuffer, qwMemInBytes);
    vFreeMemPageAligned (pqwTimestamps, qwTSBufferLen_bytes);

    return EXIT_SUCCESS;
    }
