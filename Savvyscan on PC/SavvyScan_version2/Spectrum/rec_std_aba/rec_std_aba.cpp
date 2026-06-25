/*
**************************************************************************

rec_std_aba.cpp                                        (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based acquisition cards with the ABA option 
installed. 

Information about the different products and their drivers can be found
online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/platform-driver-and-series-differences

Shows ABA data acquisition using standard ABA mode. Slow data is
acquired until all segments of the multiple recording are taken.
If timestamp is installed the corresponding timestamp values are also read
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
    int   i;

    // ABA Recording Setup
    bSpcMSetupModeRecStdABA (pstCard, ((int64) 1 << pstCard->lMaxChannels) - 1, llMemsize, lSegmentsize, lPosttrigger, 512);

    // we try to set the samplerate to 1/4 of maximum samplerate on internal PLL, no clock output, 
    // if the card can't run 10 MHz it is set to maximum sampling rate internally
    bSpcMSetupClockPLL (pstCard, pstCard->llMaxSamplerate / 4, false);
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstCard->llSetSamplerate / 1000000);

    // we set trigger to external positive edge, please connect the trigger line!
    printf ("\n!!! Using external trigger - please connect a signal to the trigger input !!!\n");
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
            printf ("Not yet implemented\n");
            break;
        }


    // if timestamp and basexio are installed we acquire the asynchronous inputs with the timestamps
    if ((pstCard->lFeatureMap & SPCM_FEAT_TIMESTAMP) && (pstCard->lFeatureMap & SPCM_FEAT_BASEXIO))
        bSpcMSetupTimestamp (pstCard, SPC_TSMODE_STANDARD | SPC_TSCNT_INTERNAL | SPC_TSXIOACQ_ENABLE, 0);

    // set up the timestamp mode to standard if timestamp is installed
    else if ((pstCard->lFeatureMap & SPCM_FEAT_TIMESTAMP) != 0)
        bSpcMSetupTimestamp (pstCard, SPC_TSMODE_STANDARD | SPC_TSCNT_INTERNAL, 0);

    }



/*
**************************************************************************
nShowAnalogData
**************************************************************************
*/

int16 nShowAnalogData (ST_SPCM_CARDINFO *pstCard, int32 lSegmentsize, int32 lPosttrigger, void* pvBuffer, bool bTimestampInstalled, bool bBaseXIO, uint64* pqwTimestamps)
    {
    int16*  ppnChannelData[SPCM_MAX_AICHANNEL];
    int     i;
    int32   lSegmentIdx;
    uint32  dwBaseXIOLines;

    for (i=0; i<pstCard->lSetChannels; i++)
        {
        ppnChannelData[i] = new int16[lSegmentsize];
        if (!ppnChannelData[i])
            return nSpcMErrorMessageStdOut (pstCard, "Memory allocation error\n", false);
        }


    // loop across all acquired segments

    if (bTimestampInstalled && bBaseXIO)
        printf ("\n%8s %12s %12s %12s %4s %15s %10s\n", "Segment", "Min", "Max", "Average", "BXIO", "Timestamp", "Timediff");
    else if (bTimestampInstalled)
        printf ("\n%8s %12s %12s %12s %15s %10s\n", "Segment", "Min", "Max", "Average", "Timestamp", "Timediff");
    else
        printf ("\n%8s %12s %12s %12s\n", "Segment", "Min", "Max", "Average");

    for (lSegmentIdx = 0; lSegmentIdx < (pstCard->llSetMemsize / lSegmentsize); lSegmentIdx++)
        {

        // split data function
        if (pstCard->lBytesPerSample == 1)
            bSpcMDemuxAnalogData (pstCard, (int8*)pvGetSegmentDataPointer (pstCard, pvBuffer, lSegmentsize, lSegmentIdx, pstCard->lBytesPerSample), lSegmentsize, ppnChannelData);
        else
            bSpcMDemuxAnalogData (pstCard, (int16*)pvGetSegmentDataPointer (pstCard, pvBuffer, lSegmentsize, lSegmentIdx, pstCard->lBytesPerSample), lSegmentsize, ppnChannelData);

        // we just look at channel 0 to keep output simple independant of the number of channels
        printf ("%8d %9.2f mV %9.2f mV %9.2f mV ", lSegmentIdx, 
            1000.0 * dSpcMIntToVoltage (pstCard, 0, TSpcMCalcMin (ppnChannelData[0], (uint32) lSegmentsize)), 
            1000.0 * dSpcMIntToVoltage (pstCard, 0, TSpcMCalcMax (ppnChannelData[0], (uint32) lSegmentsize)), 
            1000.0 * dSpcMIntToVoltage (pstCard, 0, dSpcMCalcAverage (ppnChannelData[0], (uint32) lSegmentsize)));

        // if basexio is installed we split the 8 asynchronous lines and clean up the timestamp
        if (bBaseXIO)
            {
            dwBaseXIOLines = (uint32) (pqwTimestamps[lSegmentIdx] >> 56);

            // unmask the BaseXIOLines from the timestamps
            pqwTimestamps[lSegmentIdx] &= ~(((uint64) 0xff000000) << 32);
            printf ("  %02x ", dwBaseXIOLines);
            }

        // print the timestamps and the difference, keeping track of the oversampling factor

        // M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
        if (bTimestampInstalled)
            {

            // M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
            if (pstCard->bM2i || pstCard->bM3i)
                printf ("%12.6f ms ", 1000.0 * (double) ((int64) pqwTimestamps[lSegmentIdx]) / pstCard->llSetSamplerate / pstCard->lOversampling);
            else if (pstCard->bM4i || pstCard->bM2p)
                printf ("%12.6f ms ", 1000.0 * (double) ((int64) pqwTimestamps[2 * lSegmentIdx]) / pstCard->llSetSamplerate);
            }

        // print the difference, starting with segment 1
        if (bTimestampInstalled && (lSegmentIdx > 0))
            {

            // M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
            if (pstCard->bM2i || pstCard->bM3i)
                printf ("%7.2f ms", 1000.0 * (double) (int64) (pqwTimestamps[lSegmentIdx] - pqwTimestamps[lSegmentIdx - 1]) / pstCard->llSetSamplerate / pstCard->lOversampling);
            else if (pstCard->bM4i || pstCard->bM2p)
                printf ("%7.2f ms", 1000.0 * (double) (int64) (pqwTimestamps[2*lSegmentIdx] - pqwTimestamps[2*lSegmentIdx - 2]) / pstCard->llSetSamplerate);
            }


        printf ("\n");
        }

    // clean up
    for (i=0; i<pstCard->lSetChannels; i++)
        delete [] (ppnChannelData[i]);

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
    int32               lStatus = 0;
    int64               llAvailABABytes, llByteABAPos;
    int64               llABATransfer = 0;
    uint32              dwABANotifySize = 4096;
    bool                bTimestampInstalled = 0, bBaseXIOInstalled = 0;
    uint64              qwMemInBytes;
    uint64              qwTSBufferLen_bytes = 0;
    uint64              qwABABufferLen_bytes = 0;
    uint64*             pqwTimestamps = NULL;
    void*               pvBuffer;
    void*               pvABABuffer;

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
    if ((stCard.lFeatureMap & SPCM_FEAT_ABA) == 0)
        return nSpcMErrorMessageStdOut (&stCard, "Error: ABA Option not installed. Examples was done especially for this option!\n", false);

    // if timestamp is installed we set a flag to support this mode in the example also
    bTimestampInstalled = ((stCard.lFeatureMap & SPCM_FEAT_TIMESTAMP) != 0);
    bBaseXIOInstalled = ((stCard.lFeatureMap & SPCM_FEAT_BASEXIO) != 0);

    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
    llMemsize =     KILO_B(64);
    lSegmentsize =  KILO_B(16);
    lPosttrigger =  16128;
    if (!stCard.bSetError)
        vDoCardSetup (&stCard, llMemsize, lSegmentsize, lPosttrigger);


    // ------------------------------------------------------------------------
    // calculate the amount of data we need and allocate memory buffer
    if (!stCard.bSetError)
        {
        qwMemInBytes = stCard.llSetMemsize * stCard.lBytesPerSample * stCard.lSetChannels;
        pvBuffer = pvAllocMemPageAligned (qwMemInBytes);
        
        qwABABufferLen_bytes = 20*dwABANotifySize;
        pvABABuffer = pvAllocMemPageAligned (qwABABufferLen_bytes);
        }

    // M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
    if (stCard.bM2i || stCard.bM3i)
        qwTSBufferLen_bytes = sizeof (uint64) * llMemsize / lSegmentsize;
    else if (stCard.bM4i || stCard.bM2p)
        qwTSBufferLen_bytes = 2ULL * sizeof (uint64) * llMemsize / lSegmentsize;

    pqwTimestamps = (uint64*)pvAllocMemPageAligned (qwTSBufferLen_bytes);

    if (!pvBuffer || (bTimestampInstalled &&!pqwTimestamps))
        return nSpcMErrorMessageStdOut (&stCard, "Memory allocation error\n", false);

    // ------------------------------------------------------------------------
    // make acquisition and get data
    if (!stCard.bSetError)
        {
        printf ("\n");

        // if using timestamps we need to start the transfer before the card start to avoid an overrun of the timestamp memory
        if (bTimestampInstalled)
            {
            printf ("Defining timestamp transfer\n");

            // M2i and M3i use 64 bit timestamps and M4i and M2p 128 bit
            if (stCard.bM2i || stCard.bM3i)
                spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_TIMESTAMP, SPCM_DIR_CARDTOPC, 0, (void*) pqwTimestamps, 0, llMemsize / lSegmentsize * sizeof (uint64));
            else if (stCard.bM4i || stCard.bM2i)
                spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_TIMESTAMP, SPCM_DIR_CARDTOPC, 0, (void*) pqwTimestamps, 0, 2 * llMemsize / lSegmentsize * sizeof (uint64));
            }
        
        // starting the ABA DMA transfer
        printf ("Defining ABA transfer\n");
        spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_ABA, SPCM_DIR_CARDTOPC, dwABANotifySize, pvABABuffer, 0, 20*dwABANotifySize);
        printf ("Starting the extra DMA transfer\n");
        spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_EXTRA_STARTDMA);

        // We'll start and wait untill the card has finished
        printf ("Starting the card\n\n");
        spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER);

        while(!stCard.bSetError && (lStatus & M2STAT_CARD_READY) != M2STAT_CARD_READY)
            {

            spcm_dwGetParam_i64 (stCard.hDrv, SPC_ABA_AVAIL_USER_LEN, &llAvailABABytes);
            spcm_dwGetParam_i64 (stCard.hDrv, SPC_ABA_AVAIL_USER_POS, &llByteABAPos);

            if ((llByteABAPos + llAvailABABytes) >= 20*dwABANotifySize)
                {
                llAvailABABytes = 20*dwABANotifySize - llByteABAPos;
                }

            spcm_dwSetParam_i64 (stCard.hDrv, SPC_ABA_AVAIL_CARD_LEN, llAvailABABytes);

            spcm_dwGetParam_i32 (stCard.hDrv, SPC_M2STATUS, &lStatus );

            if (llAvailABABytes > 0)
                {
                llABATransfer = llABATransfer + llAvailABABytes;
                printf("\rABA data transfered so far: %lld; ", llABATransfer);
                printf("Actually available ABA data: %lld", llAvailABABytes);
                }

            }

        // we define the buffer for transfer and start the DMA transfer
        printf ("\n\nStarting the DMA transfer and waiting until data is in PC memory\n");
        spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, 0, pvBuffer, 0, qwMemInBytes);
        spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);

        // check for error code
        if (spcm_dwGetErrorInfo_i32 (stCard.hDrv, NULL, NULL, szBuffer))
            {
            vFreeMemPageAligned (pvBuffer, qwMemInBytes);
            vFreeMemPageAligned (pvABABuffer, qwABABufferLen_bytes);
            return nSpcMErrorMessageStdOut (&stCard, szBuffer, false);
            }
        printf ("... acquisition ended, data has been transferred to PC memory\n");
        }
            // wait for the timestamps (should be already done as we started the transfer before the card start)
            if (bTimestampInstalled)
                {
                spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_EXTRA_WAITDMA);
                printf ("... timestamps have been transferred to PC memory\n");
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
                nShowAnalogData (&stCard, lSegmentsize, lPosttrigger, pvBuffer, bTimestampInstalled, bBaseXIOInstalled, pqwTimestamps);
                break;

            case DigitalIn:
            case DigitalIO:
                break;
            }

        if (bTimestampInstalled)
            vFreeMemPageAligned (pqwTimestamps, qwTSBufferLen_bytes);
        }

    spcm_dwInvalidateBuf (stCard.hDrv, SPCM_BUF_DATA);
    spcm_dwInvalidateBuf (stCard.hDrv, SPCM_BUF_ABA);
    spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_STOP | M2CMD_DATA_STOPDMA );

    vFreeMemPageAligned (pvBuffer, qwMemInBytes);
    vFreeMemPageAligned (pvABABuffer, qwABABufferLen_bytes);

    // ------------------------------------------------------------------------
    // print error information if an error occured
    if (stCard.bSetError)
        return nSpcMErrorMessageStdOut (&stCard, "An error occured while programming the card:\n", true);


    // clean up and close the driver
    vSpcMCloseCard (&stCard);

    return EXIT_SUCCESS;
    }

