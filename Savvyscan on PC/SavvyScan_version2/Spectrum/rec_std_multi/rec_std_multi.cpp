/*
**************************************************************************

rec_std_multi.cpp                                        (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based acquisition cards with the option 
Multiple Recording installed. 

Information about the different products and their drivers can be found
online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/platform-driver-and-series-differences

Shows Multiple Recoridng data acquisition using single mode (one shot). If
timestamp is installed the corresponding timestamp values are also read
out and displayed.

If Timestamp and BaseXIO are installed the BaseXIO lines are set to the
timestamp acquisition mode and are sampled on every trigger event. The
samples BaseXIO lines are also shown.
  
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

#include "../common/ostools/spcm_oswrap.h"
#include "../common/ostools/spcm_ostools.h"

// ----- standard c include files -----
#include <stdio.h>
#include <stdlib.h>



/*
**************************************************************************
vDoCardSetup
**************************************************************************
*/

void vDoCardSetup (ST_SPCM_CARDINFO *pstCard, int64 llMemsize, int32 lSegmentsize, int32 lPosttrigger)
    {
    int     i;
    int64   llChannelMask;

    // Multiple Recording setup
    if (pstCard->lMaxChannels >= 64)
        llChannelMask = -1; // -1 is all bits set to 1 = 0xffffffffffffffff
    else
        llChannelMask = ((int64) 1 << pstCard->lMaxChannels) - 1;

    bSpcMSetupModeRecStdMulti (pstCard, llChannelMask, llMemsize, lSegmentsize, lPosttrigger);

    // we try to set the samplerate to 1/4 of maximum samplerate on internal PLL, no clock output, 
    // if the card can't run 10 MHz it is set to maximum sampling rate internally
    bSpcMSetupClockPLL (pstCard, pstCard->llMaxSamplerate / 4, false);
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstCard->llSetSamplerate / 1000000);

    // we set trigger to external positive edge, please connect the trigger line!
    printf ("\n!!! Using external trigger - please connect a signal to the trigger input !!!\n\n");
    bSpcMSetupTrigExternal (pstCard, SPC_TM_POS, false, 0);

    // type dependent card setup
    switch (pstCard->eCardFunction)
        {

        // analog acquisition card setup
        case AnalogIn:

            // program all input channels to +/-1 V and 50 ohm termination (if it's available)
            for (i=0; i < pstCard->lMaxChannels; i++)
                if (pstCard->bM2i || pstCard->bM2p)
                    bSpcMSetupInputChannel (pstCard, i, 1000, true);
                else if (pstCard->bM3i || pstCard->bM4i)
                    bSpcMSetupPathInputCh (pstCard, i, 0, 1000, 0, false, true, true);
            break;

        // digital acquisition card setup
        case DigitalIn:
        case DigitalIO:

            // set all input channel groups to 110 ohm termination (if it's available)
            for (i=0; i < pstCard->uCfg.stDIO.lGroups; i++)
                bSpcMSetupDigitalInput (pstCard, i, true);
            break;
        }


    // if timestamp and basexio are installed we acquire the asynchronous inputs with the timestamps
    if ((pstCard->lFeatureMap & SPCM_FEAT_TIMESTAMP) && (pstCard->lFeatureMap & SPCM_FEAT_BASEXIO))
        bSpcMSetupTimestamp (pstCard, SPC_TSMODE_STANDARD | SPC_TSCNT_INTERNAL | SPC_TSXIOACQ_ENABLE, 0);

    // set up the timestamp mode to standard if timestamp is installed
    else if ((pstCard->lFeatureMap & SPCM_FEAT_TIMESTAMP) != 0)
        {
        int32 lAvailTSModes = 0;
        spcm_dwGetParam_i32 (pstCard->hDrv, SPC_TIMESTAMP_AVAILMODES, &lAvailTSModes);

        int32 lTSCmd = SPC_TSMODE_STANDARD | SPC_TSCNT_INTERNAL;
        if (lAvailTSModes & SPC_TSFEAT_TRGSRC)
            lTSCmd |= SPC_TSFEAT_TRGSRC; // acquire trigger sources
        if (lAvailTSModes & SPC_TSXIOACQ_ENABLE)
            lTSCmd |= SPC_TSXIOACQ_ENABLE; // acquire bits on multipurpose lines
        bSpcMSetupTimestamp (pstCard, lTSCmd, 0);
        }
    }



/*
**************************************************************************
nShowAnalogData
**************************************************************************
*/

template < typename T >
int16 nShowAnalogData (ST_SPCM_CARDINFO *pstCard, int32 lSegmentsize, int32 lPosttrigger, void* pvBuffer, bool bTimestampInstalled, bool bBaseXIO, bool bTrigSrc, uint64* pqwTimestamps)
    {
    T*  ppTChannelData[SPCM_MAX_AICHANNEL];
    int     i;
    int32   lSegmentIdx;
    uint32  dwBaseXIOLines;

    for (i = 0; i < pstCard->lSetChannels; i++)
        {
        ppTChannelData[i] = new T[lSegmentsize];
        if (!ppTChannelData[i])
            return nSpcMErrorMessageStdOut (pstCard, "Memory allocation error\n", false);
        }


    // loop across all acquired segments
    if (bTimestampInstalled && bBaseXIO)
        printf ("\n%8s %12s %12s %12s %4s %15s %10s\n", "Segment", "Min", "Max", "Average", "BXIO", "Timestamp", "Timediff");
    else if (bTimestampInstalled && bTrigSrc)
        printf ("\n%8s %12s %12s %12s %15s %10s %16s\n", "Segment", "Min", "Max", "Average", "Timestamp", "Timediff", "Trigger Source");
    else if (bTimestampInstalled)
        printf ("\n%8s %12s %12s %12s %15s %10s\n", "Segment", "Min", "Max", "Average", "Timestamp", "Timediff");
    else
        printf ("\n%8s %12s %12s %12s\n", "Segment", "Min", "Max", "Average");
    for (lSegmentIdx = 0; lSegmentIdx < (pstCard->llSetMemsize / lSegmentsize); lSegmentIdx++)
        {

        // split data function
        bSpcMDemuxAnalogData (pstCard, (T*)pvGetSegmentDataPointer (pstCard, pvBuffer, lSegmentsize, lSegmentIdx, pstCard->lBytesPerSample), lSegmentsize, ppTChannelData);

        // we just look at channel 0 to keep output simple independant of the number of channels
        printf ("%8d %9.2f mV %9.2f mV %9.2f mV ", lSegmentIdx, 
            1000.0 * dSpcMIntToVoltage (pstCard, 0, TSpcMCalcMin (ppTChannelData[0], (uint32) lSegmentsize)), 
            1000.0 * dSpcMIntToVoltage (pstCard, 0, TSpcMCalcMax (ppTChannelData[0], (uint32) lSegmentsize)), 
            1000.0 * dSpcMIntToVoltage (pstCard, 0, dSpcMCalcAverage (ppTChannelData[0], (uint32) lSegmentsize)));


        // if basexio is installed we split the 8 asynchronous lines and clean up the timestamp
        if (bBaseXIO)
            {
            dwBaseXIOLines = (uint32) (pqwTimestamps[lSegmentIdx] >> 56);

            // unmask the BaseXIOLines from the timestamps
            pqwTimestamps[lSegmentIdx] &= ~(((uint64) 0xff000000) << 32);
            printf ("  %02x ", dwBaseXIOLines);
            }

        // print the timestamps and the difference, keeping track of the oversampling factor
        if (bTimestampInstalled)
            {

            // M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
            if (pstCard->bM2i || pstCard->bM3i)
                printf ("%12.6f ms ", 1000.0 * (double) ((int64) pqwTimestamps[lSegmentIdx]) / pstCard->llSetSamplerate / pstCard->lOversampling);
            else if (pstCard->bM4i || pstCard->bM2p)
                printf ("%12.6f ms ", 1000.0 * (double) ((int64) pqwTimestamps[2 * lSegmentIdx]) / pstCard->llSetSamplerate);

            // print the difference, starting with segment 1
            if (lSegmentIdx > 0)
                {

                // M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
                if (pstCard->bM2i || pstCard->bM3i)
                    printf ("%7.2f ms", 1000.0 * (double) (int64) (pqwTimestamps[lSegmentIdx] - pqwTimestamps[lSegmentIdx - 1]) / pstCard->llSetSamplerate / pstCard->lOversampling);
                else if (pstCard->bM4i || pstCard->bM2p)
                    printf ("%7.2f ms", 1000.0 * (double) (int64) (pqwTimestamps[2*lSegmentIdx] - pqwTimestamps[2*lSegmentIdx - 2]) / pstCard->llSetSamplerate);
                }
            else
                printf ("          "); // for trigger source alignment

            // print trigger sources for each segment
            if (bTrigSrc)
                {
                if ((pqwTimestamps[2*lSegmentIdx + 1]) & SPC_TRGSRC_MASK_CH0)
                    printf (" Ch0");
                if ((pqwTimestamps[2*lSegmentIdx + 1]) & SPC_TRGSRC_MASK_CH1)
                    printf (" Ch1");
                if ((pqwTimestamps[2*lSegmentIdx + 1]) & SPC_TRGSRC_MASK_CH2)
                    printf (" Ch2");
                if ((pqwTimestamps[2*lSegmentIdx + 1]) & SPC_TRGSRC_MASK_CH3)
                    printf (" Ch3");
                if ((pqwTimestamps[2*lSegmentIdx + 1]) & SPC_TRGSRC_MASK_CH4)
                    printf (" Ch4");
                if ((pqwTimestamps[2*lSegmentIdx + 1]) & SPC_TRGSRC_MASK_CH5)
                    printf (" Ch5");
                if ((pqwTimestamps[2*lSegmentIdx + 1]) & SPC_TRGSRC_MASK_CH6)
                    printf (" Ch6");
                if ((pqwTimestamps[2*lSegmentIdx + 1]) & SPC_TRGSRC_MASK_CH7)
                    printf (" Ch7");
                if ((pqwTimestamps[2*lSegmentIdx + 1]) & SPC_TRGSRC_MASK_EXT0)
                    printf (" EXT0");
                if ((pqwTimestamps[2*lSegmentIdx + 1]) & SPC_TRGSRC_MASK_EXT1)
                    printf (" EXT1");
                if ((pqwTimestamps[2*lSegmentIdx + 1]) & SPC_TRGSRC_MASK_FORCE)
                    printf (" FORCE");
                }
            }


        printf ("\n");
        }

    // clean up
    for (i = 0; i < pstCard->lSetChannels; i++)
        delete [] (ppTChannelData[i]);

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
        delete (ppbyChannelData[i]);

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
    int32               lSegmentsize, lPosttrigger;
    int32               lOversampling = 1;
    bool                bTimestampInstalled, bBaseXIOInstalled;
    bool                bTrigSrcAvailable = false;
    uint64              qwMemInBytes;
    uint64              qwTSBufferLen_bytes = 0;
    uint64*             pqwTimestamps = NULL;
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
    if ((stCard.lFeatureMap & SPCM_FEAT_MULTI) == 0)
        return nSpcMErrorMessageStdOut (&stCard, "Error: Multiple Recording Option not installed. Examples was done especially for this option!\n", false);

    // if timestamp is installed we set a flag to support this mode in the example also
    bTimestampInstalled = ((stCard.lFeatureMap & SPCM_FEAT_TIMESTAMP) != 0);
    bBaseXIOInstalled = ((stCard.lFeatureMap & SPCM_FEAT_BASEXIO) != 0);
    if (stCard.bM4i || stCard.bM2p)
        bTrigSrcAvailable = true;



    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
    llMemsize =     KILO_B(256);
    lSegmentsize =  KILO_B(16);
    lPosttrigger =  KILO_B(15);
    if (!stCard.bSetError)
        vDoCardSetup (&stCard, llMemsize, lSegmentsize, lPosttrigger);



    // ------------------------------------------------------------------------
    // calculate the amount of data we need and allocate memory buffer
    if (!stCard.bSetError)
        {
        qwMemInBytes = stCard.llSetMemsize * stCard.lBytesPerSample * stCard.lSetChannels;
        pvBuffer = pvAllocMemPageAligned (qwMemInBytes);
 
        if (bTimestampInstalled)
            {

            // M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
            if (stCard.bM2i || stCard.bM3i)
                qwTSBufferLen_bytes = sizeof (uint64) * llMemsize / lSegmentsize;
            else if (stCard.bM4i || stCard.bM2p)
                qwTSBufferLen_bytes = 2ULL*sizeof (uint64) * llMemsize / lSegmentsize;

            pqwTimestamps = (uint64*)pvAllocMemPageAligned (qwTSBufferLen_bytes);
            }

        if (!pvBuffer || (bTimestampInstalled &&!pqwTimestamps))
            return nSpcMErrorMessageStdOut (&stCard, "Memory allocation error\n", false);
        }



    // ------------------------------------------------------------------------
    // make acquisition and get data
    if (!stCard.bSetError)
        {
        printf ("\n");

        // if using timestamps we need to start the transfer before the card start to avoid an overrun of the timestamp memory
        if (bTimestampInstalled)
            {
            printf ("Starting the timestamp DMA transfer\n");

            // M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
            if (stCard.bM2i || stCard.bM3i)
                spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_TIMESTAMP, SPCM_DIR_CARDTOPC, 0, (void*) pqwTimestamps, 0, llMemsize / lSegmentsize * sizeof (uint64));
            else if (stCard.bM4i || stCard.bM2p)
                spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_TIMESTAMP, SPCM_DIR_CARDTOPC, 0, (void*) pqwTimestamps, 0, 2 * llMemsize / lSegmentsize * sizeof (uint64));

            spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_EXTRA_STARTDMA);
            }

        // We'll start and wait untill the card has finished or until a timeout occurs
        spcm_dwSetParam_i32 (stCard.hDrv, SPC_TIMEOUT, 5000);
        printf ("Starting the card and waiting for ready interrupt\n");
        if (spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY) == ERR_TIMEOUT)
            {
            vFreeMemPageAligned (pvBuffer, qwMemInBytes);
            return nSpcMErrorMessageStdOut (&stCard, "... Timeout\n", false);
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

            // wait for the timestamps (should be already done as we started the transfer before the card start)
            if (bTimestampInstalled)
                {
                spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_EXTRA_WAITDMA);
                printf ("... timestamps have been transferred to PC memory\n");
                }
            }
        }



    // ------------------------------------------------------------------------
    // we go through the segments, split the data in separate channels and show some results
    if (!stCard.bSetError)
        {

        // some additional information on the acquisition
        printf ("\nData information:\n=================\n");
        printf ("Each segment is %.3f ms long\n", 1000.0 * lSegmentsize / stCard.llSetSamplerate);
        printf ("Maximum pulse repetition frequency to reach with this setting is %.2f Hz\n", (double) stCard.llSetSamplerate / lSegmentsize);

        switch (stCard.eCardFunction)
            {
            case AnalogIn:
                if (stCard.lBytesPerSample == 2)
                    nShowAnalogData < int16 >(&stCard, lSegmentsize, lPosttrigger, pvBuffer, bTimestampInstalled, bBaseXIOInstalled, bTrigSrcAvailable, pqwTimestamps);
                else
                    nShowAnalogData < int8 >(&stCard, lSegmentsize, lPosttrigger, pvBuffer, bTimestampInstalled, bBaseXIOInstalled, bTrigSrcAvailable, pqwTimestamps);
                break;

            case DigitalIn:
            case DigitalIO:

                // show first 10 samples for each channel
                nShowDigitalData (&stCard, pvBuffer, 10);
                break;
            }

        if (bTimestampInstalled)
            vFreeMemPageAligned (pqwTimestamps, qwTSBufferLen_bytes);
        }
    


    // ------------------------------------------------------------------------
    // print error information if an error occured
    if (stCard.bSetError)
        return nSpcMErrorMessageStdOut (&stCard, "An error occured while programming the card:\n", true);


    // clean up and close the driver
    vSpcMCloseCard (&stCard);
    vFreeMemPageAligned (pvBuffer, qwMemInBytes);

    return EXIT_SUCCESS;
    }

