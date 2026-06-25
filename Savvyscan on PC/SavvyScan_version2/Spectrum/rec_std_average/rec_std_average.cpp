/*
**************************************************************************

rec_std_average.cpp                                      (c) Spectrum GmbH

**************************************************************************

Example for all acquisition cards with the option Average installed. 

Information about the different products and their drivers can be found
online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/platform-driver-and-series-differences

Shows Average data acquisition using single mode (one shot). If
timestamp is installed the corresponding timestamp values are also read
out and displayed.
  
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

void vDoCardSetup (ST_SPCM_CARDINFO *pstCard, int64 llMemsize, int32 lSegmentsize, int32 lPosttrigger, int32 lAverages)
    {
    int     i;
    int64   llChannelMask;

    // Channel enable setup
    llChannelMask = ((int64) 1 << pstCard->lMaxChannels) - 1;

    bSpcMSetupModeRecStdAverage(pstCard, llChannelMask, llMemsize, lSegmentsize, lPosttrigger, lAverages);

    // we try to set the samplerate to maximum samplerate on internal PLL, no clock output, 
    bSpcMSetupClockPLL (pstCard, pstCard->llMaxSamplerate , false);
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstCard->llSetSamplerate / 1000000);

    // we set trigger to external positive edge, please connect the trigger line!
    printf ("\n!!! Using external trigger - please connect a signal to the trigger input !!!\n\n");    
    bSpcMSetupTrigExternal (pstCard, SPC_TM_POS, false, 0);

    // program all input channels to +/-1 V and 50 ohm termination (if it's available)
    for (i=0; i < pstCard->lMaxChannels; i++)
        bSpcMSetupPathInputCh (pstCard, i, 0, 1000, 0);
    }



/*
**************************************************************************
nShowAnalogData
**************************************************************************
*/

int16 nShowAnalogData (ST_SPCM_CARDINFO *pstCard, int32 lSegmentsize, int32 lPosttrigger, void* pvBuffer, bool bAcquireTimestamps, uint64* pqwTimestamps, int32 lAverages)
    {
    int32*  pplChannelData[SPCM_MAX_AICHANNEL];
    int     i;
    int32   lSegmentIdx;

    for (i=0; i<pstCard->lSetChannels; i++)
        {
        pplChannelData[i] = new int32[lSegmentsize];
        if (!pplChannelData[i])
            return nSpcMErrorMessageStdOut (pstCard, "Memory allocation error\n", false);
        }

    // loop across all acquired segments
    printf ("\n%8s %12s %12s %12s\n", "Segment", "Min", "Max", "Average");

    // Display some information per segment
    for (lSegmentIdx = 0; lSegmentIdx < (pstCard->llSetMemsize / lSegmentsize); lSegmentIdx++)
        {
        // split data function to have easier per channel access to the samples
        bSpcMDemuxAnalogData (pstCard, (int32*)pvGetSegmentDataPointer (pstCard, pvBuffer, lSegmentsize, lSegmentIdx, 4), lSegmentsize, pplChannelData);
        
        // we just look at one certain to keep output simple independant of the number of activated channels
        uint8 byChannel = 0;
        printf ("%8d %9.2f mV %9.2f mV %9.2f mV", lSegmentIdx, 
            1000.0 * dSpcMIntToVoltage (pstCard, 0, TSpcMCalcMin     (pplChannelData[byChannel], (uint32) lSegmentsize)), 
            1000.0 * dSpcMIntToVoltage (pstCard, 0, TSpcMCalcMax     (pplChannelData[byChannel], (uint32) lSegmentsize)), 
            1000.0 * dSpcMIntToVoltage (pstCard, 0, dSpcMCalcAverage (pplChannelData[byChannel], (uint32) lSegmentsize)));

        printf ("\n");
        }

    // clean up
    for (i=0; i<pstCard->lSetChannels; i++)
        delete [] (pplChannelData[i]);

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
    int32               lSegmentsize, lPosttrigger, lAveragesPerSegment, lNumberOfSegments;
    int32               lOversampling = 1;
    uint64              qwMemInBytes;
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
    if (!stCard.bM4i)
        return nSpcMErrorMessageStdOut (&stCard, "Error. Averaging only works on M4i cards\n", false);
    if (stCard.eCardFunction != AnalogIn)
        return nSpcMErrorMessageStdOut (&stCard, "Error: Card function not supported by this example\n", false);
    if ((stCard.lExtFeatureMap & SPCM_FEAT_EXTFW_SEGAVERAGE) == 0)
        return nSpcMErrorMessageStdOut (&stCard, "Error: Average Option not installed. Examples was done especially for this option!\n", false);


    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values  
    lNumberOfSegments =     16;
    lSegmentsize =          KILO_B(8);
    lPosttrigger =          lSegmentsize - KILO_B(1);
    lAveragesPerSegment =   16;

    llMemsize = (int64) lSegmentsize * (int64) lNumberOfSegments;

    if (!stCard.bSetError)
        vDoCardSetup (&stCard, llMemsize, lSegmentsize, lPosttrigger, lAveragesPerSegment);

    // ------------------------------------------------------------------------
    // calculate the amount of data we need and allocate memory buffer
    if (!stCard.bSetError)
        {
        // for averaging the number of bytes per sample is fixed to 4 (32 bit samples)
        qwMemInBytes = (uint64) stCard.llSetMemsize * (uint64) sizeof(int32) * (uint64) stCard.lSetChannels;

        pvBuffer = pvAllocMemPageAligned (qwMemInBytes);

        if (!pvBuffer)
            return nSpcMErrorMessageStdOut (&stCard, "Memory allocation error\n", false);
        }

    // ------------------------------------------------------------------------
    // make acquisition and get data
    if (!stCard.bSetError)
        {
        printf ("\n");

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

        nShowAnalogData (&stCard, lSegmentsize, lPosttrigger, pvBuffer, false, NULL, lAveragesPerSegment);
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

