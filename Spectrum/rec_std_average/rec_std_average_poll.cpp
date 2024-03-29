/*
**************************************************************************

rec_std_average_poll.cpp                                   (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based (M2i) acquisition cards with the option 
Multiple Recording installed. 

Shows Multiple Recoridng data acquisition using single mode (one shot). If
timestamp is installed the corresponding timestamp values are also read
out and displayed.

If Timestamp and BaseXIO are installed the BaseXIO lines are set to the
timestamp acquisition mode and are sampled on every trigger event. The
samples BaseXIO lines are also shown.
  
Feel free to use this source for own projects and modify it in any kind

This example shows the polling of the status, no wait function is used

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

    // Multiple Recording setup
    bSpcMSetupModeRecStdMulti (pstCard, (1 << pstCard->lMaxChannels) - 1, llMemsize, lSegmentsize, lPosttrigger);

    // we try to set the samplerate to 10 MHz on internal PLL, no clock output, 
    // if the card can't run 10 MHz it is set to maximum sampling rate internally
    bSpcMSetupClockPLL (pstCard, MEGA(10), false);
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
                bSpcMSetupInputChannel (pstCard, i, 1000, true);
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
        printf ("\n%8s %12s %12s %12s %15s %10s\n", "Segment", "Min", "Max", "Average");
    for (lSegmentIdx = 0; lSegmentIdx < (pstCard->llSetMemsize / lSegmentsize); lSegmentIdx++)
        {

        // split data function
        bSpcMDemuxAnalogData (pstCard, pvGetSegmentDataPointer (pstCard, pvBuffer, lSegmentsize, lSegmentIdx, pstCard->lBytesPerSample), lSegmentsize, ppnChannelData);

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
        if (bTimestampInstalled)
			{

			// M2i and M3i use 64 bit timestamps and M4i/M4x/M2p 128 bit
			if (pstCard->bM2i || pstCard->bM3i)
				printf ("%12.6f ms ", 1000.0 * (double) ((int64) pqwTimestamps[lSegmentIdx]) / pstCard->llSetSamplerate / pstCard->lOversampling);
			else if (pstCard->bM4i || pstCard->bM2p)
				printf ("%12.6f ms ", 1000.0 * (double) ((int64) pqwTimestamps[2 * lSegmentIdx]) / pstCard->llSetSamplerate);
			}

        // print the difference, starting with segment 1
        if (bTimestampInstalled && (lSegmentIdx > 0))
			{

			// M2i and M3i use 64 bit timestamps and M4i/M4x/M2p 128 bit
			if (pstCard->bM2i || pstCard->bM3i)
	            printf ("%7.2f ms", 1000.0 * (double) (int64) (pqwTimestamps[lSegmentIdx] - pqwTimestamps[lSegmentIdx - 1]) / pstCard->llSetSamplerate / pstCard->lOversampling);
			else if (pstCard->bM4i || pstCard->bM2p)
	            printf ("%7.2f ms", 1000.0 * (double) (int64) (pqwTimestamps[2*lSegmentIdx] - pqwTimestamps[2*lSegmentIdx - 2]) / pstCard->llSetSamplerate);
			}


        printf ("\n");
        }

    // clean up
    for (i=0; i<pstCard->lSetChannels; i++)
        delete (ppnChannelData[i]);

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
    bool                bBaseXIOInstalled;
    uint64              qwMemInBytes;
    void*               pvBuffer;
    uint32              dwErr;


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
    if (stCard.eCardFunction != AnalogIn)
        return nSpcMErrorMessageStdOut (&stCard, "Error: Card function not supported by this example\n", false);
    if ((stCard.lFeatureMap & SPCM_FEAT_MULTI) == 0)
        return nSpcMErrorMessageStdOut (&stCard, "Error: Multiple Recording Option not installed. Examples was done especially for this option!\n", false);

    // if timestamp is installed we set a flag to support this mode in the example also
    bBaseXIOInstalled = ((stCard.lFeatureMap & SPCM_FEAT_BASEXIO) != 0);



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

        if (!pvBuffer)
            return nSpcMErrorMessageStdOut (&stCard, "Memory allocation error\n", false);
        }



    // ------------------------------------------------------------------------
    // make acquisition and get data (program a number of loops if wanted
    int32 lLoops = 3;
    while (lLoops--)
        {
        if (!stCard.bSetError)
            {
            printf ("\n");

            // We'll define the buffer for data to start everything together
            spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, 0, pvBuffer, 0, qwMemInBytes);

            printf ("Starting the card, the data transfer and poll\n");
            spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_DATA_STARTDMA);


            // this is our status polling loop
            bool bReady =           false;
            bool bDataReady =       false;
            int32 lOldStatus = 0, lStatus;
            do
                {
                dwErr = spcm_dwGetParam_i32 (stCard.hDrv, SPC_M2STATUS, &lStatus);

                // check for changed status and set flags
                if (lOldStatus != lStatus)
                    {
                    if (lStatus & M2STAT_CARD_PRETRIGGER)
                        printf ("Armed ");

                    if (lStatus & M2STAT_CARD_TRIGGER)
                        printf ("1stTrigger ");

                    if (lStatus & M2STAT_CARD_READY)
                        {
                        printf ("CardReady ");
                        bReady = true;
                        }

                    if (lStatus & M2STAT_DATA_END)
                        {
                        printf ("DataTransferred ");
                        bDataReady = true;
                        }

                    printf ("\n");
                    }

                lOldStatus = lStatus;
                }
            while (!dwErr && (!bReady || !bDataReady));

            // check for error code
            if (spcm_dwGetErrorInfo_i32 (stCard.hDrv, NULL, NULL, szBuffer))
                {
                vFreeMemPageAligned (pvBuffer, qwMemInBytes);
                return nSpcMErrorMessageStdOut (&stCard, szBuffer, false);
                }

            printf ("... acquisition ended, data has been transferred to PC memory\n");
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
                    nShowAnalogData (&stCard, lSegmentsize, lPosttrigger, pvBuffer, false, bBaseXIOInstalled, NULL);
                    break;

                case DigitalIn:
                case DigitalIO:
                    break;
                }
            }
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

