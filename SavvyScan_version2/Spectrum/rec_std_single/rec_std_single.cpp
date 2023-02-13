/*
**************************************************************************

rec_std_single.cpp                                       (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based acquisition cards. 

Information about the different products and their drivers can be found
online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/platform-driver-and-series-differences

Shows standard data acquisition using single mode (one shot)
  
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

void vDoCardSetup (ST_SPCM_CARDINFO *pstCard)
    {
    int i;
    int64 llChannelMask;

    // set mask for maximal channels
    if (pstCard->lMaxChannels >= 64)
        llChannelMask = -1; // -1 is all bits set to 1 = 0xffffffffffffffff
    else
        llChannelMask = ((int64) 1 << pstCard->lMaxChannels) - 1;

    // standard single, all channels, memsize=16k, posttrigge=8k -> pretrigger=8k
    bSpcMSetupModeRecStdSingle (pstCard, llChannelMask, KILO_B(16), KILO_B(8));

    // we try to set the samplerate to 10 max/4 internal PLL, no clock output
    bSpcMSetupClockPLL (pstCard, pstCard->llMaxSamplerate / 4, false);
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstCard->llSetSamplerate / 1000000);

    // we set software trigger, no trigger output
    bSpcMSetupTrigSoftware (pstCard, false);


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

            // activate digital inputs if available
            if (pstCard->bM2i)
                {
                if ((pstCard->lBytesPerSample > 1) && (pstCard->lFeatureMap & SPCM_FEAT_DIGITAL))
                    spcm_dwSetParam_i32 (pstCard->hDrv, SPC_READDIGITAL, 1); 
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
    }




/*
**************************************************************************
nShowAnalogData 
**************************************************************************
*/

int16 nShowAnalogData (ST_SPCM_CARDINFO *pstCard, void* pvBuffer)
    {
    int     i;
    double  dAverage;
    int16   nMin;
    int16   nMax;
    int16*  ppnChannelData[SPCM_MAX_AICHANNEL];

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

    // calc and print some figures
    printf ("\nSome details of the acquired data:\n");
    for (i=0; i<pstCard->lSetChannels; i++)
        {
        dAverage =  dSpcMCalcAverage (ppnChannelData[i], (uint32) pstCard->llSetMemsize);
        nMin =      TSpcMCalcMin (ppnChannelData[i], (uint32) pstCard->llSetMemsize);
        nMax =      TSpcMCalcMax (ppnChannelData[i], (uint32) pstCard->llSetMemsize);
        printf ("Ch %2d\n  Min: %7d (%9.2f mV)\n  Max: %7d (%9.2f mV)\n  Av:  %7.2f (%9.2f mV)\n\n", i, 
            nMin, 1000.0 * dSpcMIntToVoltage (pstCard, i, nMin), 
            nMax, 1000.0 * dSpcMIntToVoltage (pstCard, i, nMax),
            dAverage, 1000.0 * dSpcMIntToVoltage (pstCard, i, dAverage));
        }

    for (i=0; i<pstCard->lSetChannels; i++)
        delete [] (ppnChannelData[i]);

    return 0;
    }

/*
**************************************************************************
nShowDigitalInputData 
**************************************************************************
*/

int16 nShowDigitalInputData (ST_SPCM_CARDINFO *pstCard, void* pvBuffer, uint32 dwNrOfSamplesToShow)
    {
    int32 lNrOfDigCh, lCh, lBitShift, i;
    uint32 dwSampleIdx;
    uint8 byDigMask, byBitVal;

    uint8* ppbyChannelData[SPCM_MAX_AICHANNEL];

    if (((pstCard->lCardType & TYP_SERIESMASK) == TYP_M4IEXPSERIES) || ((pstCard->lCardType & TYP_SERIESMASK) == TYP_M4XEXPSERIES))
        {
        lNrOfDigCh = 1;
        }
    else
        {
        lNrOfDigCh = 16 - pstCard->uCfg.stAI.lResolution;
        }

    // allocate channel data
    for (i = 0; i < pstCard->lSetChannels; i++)
        {
        ppbyChannelData[i] = new uint8[(int32)pstCard->llSetMemsize];
        if (!ppbyChannelData[i])
            return nSpcMErrorMessageStdOut (pstCard, "Memory allocation error\n", false);
        }

    // demux digital input data
    bSpcMDemuxDigitalInputDataToUInt8 (pstCard, pvBuffer, (int32)pstCard->llSetMemsize, ppbyChannelData);

    // print digital inputs
    printf ("\nDigital inputs:\n");

    for (dwSampleIdx = 0; dwSampleIdx < dwNrOfSamplesToShow; dwSampleIdx++)
        {
        printf ("[D0 - D%d] : ", pstCard->lSetChannels*lNrOfDigCh - 1);

        for (lCh = 0; lCh < pstCard->lSetChannels; lCh++)
            {
            byDigMask = ppbyChannelData[lCh][dwSampleIdx]; 	
            lBitShift = 0;

            do 
                {
                byBitVal  = byDigMask & 0x01;

                printf ("%u", byBitVal);

                byDigMask = byDigMask >> 1;
                lBitShift++;
                }while (lBitShift < lNrOfDigCh);

            printf (" ");
            }
            printf ("\n");
        }

    for (i = 0; i < pstCard->lSetChannels; i++)
        delete [] (ppbyChannelData[i]);

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
    uint64              qwMemInBytes;
    void*               pvBuffer;
    int32               lValue;

    // ------------------------------------------------------------------------
    // init card number 0 (the first card in the system), get some information and print it
    // uncomment the second line and replace the IP address to use remote
    // cards like in a digitizerNETBOX
    if (bSpcMInitCardByIdx (&stCard, 0))
    //if (bSpcMInitCardByIdx (&stCard, "192.168.169.42", 0))
        {
        printf (pszSpcMPrintDocumentationLink (&stCard, szBuffer, sizeof (szBuffer)));
        printf (pszSpcMPrintCardInfo (&stCard, szBuffer, sizeof (szBuffer)));
        }
    else
        return nSpcMErrorMessageStdOut (&stCard, "Error: Could not open card\n", true);

    // check whether we support this card type in the example
    if (stCard.eCardFunction != AnalogIn && stCard.eCardFunction != DigitalIn && stCard.eCardFunction != DigitalIO)
        return nSpcMErrorMessageStdOut (&stCard, "Error: Card function not supported by this example\n", false);



    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
    if (!stCard.bSetError)
        vDoCardSetup (&stCard);


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
                qwMemInBytes = stCard.lSetChannels / 8 * stCard.llSetMemsize;     
                break;
            }

        pvBuffer = pvAllocMemPageAligned (qwMemInBytes);
        if (!pvBuffer)
            return nSpcMErrorMessageStdOut (&stCard, "Memory allocation error\n", false);
        }



    // ------------------------------------------------------------------------
    // make acquisition and get data
    if (!stCard.bSetError)
        {

        // We'll start and wait untill the card has finished or until a timeout occurs
        spcm_dwSetParam_i32 (stCard.hDrv, SPC_TIMEOUT, 5000);
        printf ("\nStarting the card and waiting for ready interrupt\n");
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
    // we split the data in separate channels and show some details
    if (!stCard.bSetError)
        switch (stCard.eCardFunction)
            {
            case AnalogIn:

                spcm_dwGetParam_i32 (stCard.hDrv, SPC_READDIGITAL, &lValue);

                // check if digital inputs are set
                if (lValue)
                    {
                    // memory allocation for analog and digital data
                    void* pvAnalogData  = (void*) new uint8[(int) qwMemInBytes]; 
                    void* pvDigitalData = (void*) new uint8[(int) (qwMemInBytes / stCard.lBytesPerSample)];
                    
                    // split data in analog and digital part
                    bSpcMSplitAnalogAndDigitalData (&stCard, pvBuffer, (uint32)(qwMemInBytes / stCard.lBytesPerSample), pvAnalogData, pvDigitalData);
                    
                    // show digital data
                    nShowDigitalInputData (&stCard, pvDigitalData, 10);

                    // show analog data
                    nShowAnalogData (&stCard, pvAnalogData);
                    }
                else
                    nShowAnalogData (&stCard, pvBuffer);
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

    return EXIT_SUCCESS;
    }

