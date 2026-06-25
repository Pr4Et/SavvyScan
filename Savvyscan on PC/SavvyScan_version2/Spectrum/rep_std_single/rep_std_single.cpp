/*
**************************************************************************

rep_std_single.cpp                                       (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based analog and digital generator cards. 

Information about the different products and their drivers can be found
online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/platform-driver-and-series-differences

Shows standard replay mode as single shot, continous or single restart
  
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

// ----- operating system dependent functions for thead, event, keyboard and mutex handling -----
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

void vDoCardSetup (ST_SPCM_CARDINFO *pstCard, int32 lReplayMode, int64 llLoops = 0)
    {
    int i;
    int64 llChannelMask;


    // set mask for maximal channels
    if (pstCard->lMaxChannels >= 64)
        llChannelMask = -1; // -1 is all bits set to 1 = 0xffffffffffffffff
    else
        llChannelMask = ((int64) 1 << pstCard->lMaxChannels) - 1;


    // we try to set the samplerate to 1 MHz (M2i) or 50 MHz (M4i) on internal PLL, no clock output
    if (pstCard->bM4i)
        bSpcMSetupClockPLL (pstCard, MEGA(50), false);
    else
        bSpcMSetupClockPLL (pstCard, MEGA(1), false);
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstCard->llSetSamplerate / 1000000);


    // setup the replay mode and the trigger
    switch (lReplayMode)
        {

        // with loops == 1: singleshot replay with software trigger
        // with loops == 0: endless continuous mode with software trigger
        case SPC_REP_STD_SINGLE:
            bSpcMSetupModeRepStdLoops  (pstCard, llChannelMask, KILO_B(64), llLoops);
            bSpcMSetupTrigSoftware (pstCard, true);

            // on M2i starting with build 1604 we can use the trigger output as a marker for each loop start
            // be sure to have the trigger output enabled for this
            if (pstCard->bM2i)
                {
                if ((pstCard->lLibVersion & 0xFFFF) >= 1604)
                    spcm_dwSetParam_i32 (pstCard->hDrv, SPC_CONTOUTMARK, 1);
                }
            break;

        // single restart (one signal on every trigger edge) with ext trigger positive edge
        case SPC_REP_STD_SINGLERESTART:
            bSpcMSetupModeRepStdSingleRestart (pstCard, llChannelMask, KILO_B(64), 0);
            bSpcMSetupTrigExternal (pstCard, SPC_TM_POS, false, 0);
            break;
        }



    // type dependent card setup
    switch (pstCard->eCardFunction)
        {

        // analog generator card setup
        case AnalogOut:

            // program all output channels to +/- 1 V with no offset
            for (i=0; i < pstCard->lMaxChannels; i++)
                bSpcMSetupAnalogOutputChannel (pstCard, i, 1000, 0, 0);
            break;

        // digital generator card setup
        case DigitalOut:
        case DigitalIO:
            for (i=0; i < pstCard->uCfg.stDIO.lGroups; i++)
                bSpcMSetupDigitalOutput (pstCard, i, SPCM_STOPLVL_LOW, 0, 3300);
            break;
        }
    }



/*
**************************************************************************
DoDataCalculation: calculates the output data
**************************************************************************
*/

bool bDoDataCalculation (ST_SPCM_CARDINFO *pstCard, void* pvBuffer)
    {
    void*   ppvChannelData[SPCM_MAX_AOCHANNEL];
    int     i;

    printf ("Calculation of output data\n");

    switch (pstCard->eCardFunction)
        {

        // analog waveform generator card, each channel gets a different waveform
        case AnalogOut:
            {
            // allocate buffers for each channel
            for (i=0; i < pstCard->lMaxChannels; i++)
                {
                if (pstCard->lBytesPerSample == 1)
                    ppvChannelData[i] = new int8[(unsigned) pstCard->llSetMemsize];
                else
                    ppvChannelData[i] = new int16[(unsigned) pstCard->llSetMemsize];
                if (!ppvChannelData[i])
                    return (nSpcMErrorMessageStdOut (pstCard, "Memory allocation error\n", false) == 0);
                }

            // calculate channel data
            for (i=0; i < pstCard->lMaxChannels; i++)
                {
                switch (i)
                    {
                    case 0: bSpcMCalcSignal (pstCard, ppvChannelData[0], (uint32) pstCard->llSetMemsize, 0, eSine);              break;
                    case 1: bSpcMCalcSignal (pstCard, ppvChannelData[1], (uint32) pstCard->llSetMemsize, 0, eTriangle);          break;
                    case 2: bSpcMCalcSignal (pstCard, ppvChannelData[2], (uint32) pstCard->llSetMemsize, 0, eSawtooth);          break;
                    case 3: bSpcMCalcSignal (pstCard, ppvChannelData[3], (uint32) pstCard->llSetMemsize, 0, eRectangle);         break;
                    case 4: bSpcMCalcSignal (pstCard, ppvChannelData[4], (uint32) pstCard->llSetMemsize, 0, eInvertedSine);      break;
                    case 5: bSpcMCalcSignal (pstCard, ppvChannelData[5], (uint32) pstCard->llSetMemsize, 0, eInvertedTriangle);  break;
                    case 6: bSpcMCalcSignal (pstCard, ppvChannelData[6], (uint32) pstCard->llSetMemsize, 0, eInvertedSawtooth);  break;
                    case 7: bSpcMCalcSignal (pstCard, ppvChannelData[7], (uint32) pstCard->llSetMemsize, 0, eInvertedRectangle); break;
                    }
                }

            // mux it into the output buffer
            bSpcMMuxData (pstCard, pvBuffer, (uint32) pstCard->llSetMemsize, ppvChannelData);

            // clean up channel buffers
            for (i=0; i < pstCard->lMaxChannels; i++)
                {
                if (pstCard->lBytesPerSample == 1)
                    delete [] (int8*)ppvChannelData[i];
                else
                    delete [] (int16*)ppvChannelData[i];
                }

            break;
            }

        // digital generator card: sine over all channels
        case DigitalOut:
        case DigitalIO:
            {
            // we need to tell the calc function the number of bytes for one complete word -> [channels/8]
            bSpcMCalcSignal (pstCard, pvBuffer, (uint32) pstCard->llSetMemsize, pstCard->lSetChannels / 8, eSine);
            break;
            }
        }


    return true;
    }

/*
**************************************************************************
bSetMultiPurposeDigOut
**************************************************************************
*/

bool bSetMultiPurposeDigOut (ST_SPCM_CARDINFO *pstCard, void* pvBuffer)
    {
    if (pstCard->uCfg.stAO.lResolution != 16)
        return false;

    int16* pnData = (int16*)pvBuffer;
    int64 llDataLength = pstCard->lSetChannels * pstCard->llSetMemsize;

    uint32 dwXMode_X0 = 0;
    uint32 dwXMode_X1 = 0;
    uint32 dwXMode_X2 = 0;
    uint32 dwXMode_X3 = 0;

    uint32 dwNumOfDigOutChannels = 0;

    // set all available and active XIO channels to digital output mode
    if (pstCard->lSetChannels >= 1 && pstCard->qwSetChEnableMap & 0x01)
        {
        dwXMode_X0 = SPCM_XMODE_DIGOUT | SPCM_XMODE_DIGOUTSRC_CH0 | SPCM_XMODE_DIGOUTSRC_BIT15;
        if (spcm_dwSetParam_i32 (pstCard->hDrv, SPCM_X0_MODE, dwXMode_X0))
            return false;
        dwNumOfDigOutChannels = 1;
        }

    if (pstCard->lSetChannels >= 2 && pstCard->qwSetChEnableMap & 0x03)
        {
        dwXMode_X1 = SPCM_XMODE_DIGOUT | SPCM_XMODE_DIGOUTSRC_CH1 | SPCM_XMODE_DIGOUTSRC_BIT15;
        if (spcm_dwSetParam_i32 (pstCard->hDrv, SPCM_X1_MODE, dwXMode_X1))
            return false;
        dwNumOfDigOutChannels = 2;
        }

    if (pstCard->lSetChannels >= 3 && pstCard->qwSetChEnableMap & 0x07)
        {
        dwXMode_X2 = SPCM_XMODE_DIGOUT | SPCM_XMODE_DIGOUTSRC_CH2 | SPCM_XMODE_DIGOUTSRC_BIT15;
        if (spcm_dwSetParam_i32 (pstCard->hDrv, SPCM_X2_MODE, dwXMode_X2))
            return false;
        dwNumOfDigOutChannels = 3;
        }

    if (pstCard->bM2p)
        {
        if (pstCard->lSetChannels >= 4 && pstCard->qwSetChEnableMap & 0x0F)
            {
            dwXMode_X3 = SPCM_XMODE_DIGOUT | SPCM_XMODE_DIGOUTSRC_CH3 | SPCM_XMODE_DIGOUTSRC_BIT15;
            if (spcm_dwSetParam_i32 (pstCard->hDrv, SPCM_X3_MODE, dwXMode_X3))
                return false;
            dwNumOfDigOutChannels = 4;
            }
        }

    bool bDigValHigh;

    for (int64 llDataIdx = 0; llDataIdx < llDataLength; llDataIdx += pstCard->lSetChannels)
        {
        if (llDataIdx < llDataLength / 2)
            bDigValHigh = true;
        else
            bDigValHigh = false;
        
        // digital output value is set by adjusting highest bit in sample
        for (uint32 dwDigChIdx = 0; dwDigChIdx < dwNumOfDigOutChannels; dwDigChIdx++)
            {
            // shift sample by one bit
            pnData[llDataIdx + dwDigChIdx] >>= 1;

            // adjust bit 15 in sample to desired digital output value 
            if (bDigValHigh)
                pnData[llDataIdx + dwDigChIdx] |= 0x8000;
            else
                pnData[llDataIdx + dwDigChIdx] &= 0x7FFF;
            }
        }

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
    uint64              qwMemInBytes;
    void*               pvBuffer = NULL;

    bool bMultiPurposeDigOut = false;

    // ------------------------------------------------------------------------
    // init card number 0 (the first card in the system), get some information and print it
    // uncomment the second line and replace the IP address to use remote
    // cards like in a generatorNETBOX
    if (bSpcMInitCardByIdx (&stCard, 0))
    //if (bSpcMInitCardByIdx (&stCard, "192.168.1.10", 0))
        {
        printf (pszSpcMPrintDocumentationLink (&stCard, szBuffer, sizeof (szBuffer)));
        printf (pszSpcMPrintCardInfo (&stCard, szBuffer, sizeof (szBuffer)));
        }
    else
        return nSpcMErrorMessageStdOut (&stCard, "Error: Could not open card\n", true);


    // check whether we support this card type in the example
    if ((stCard.eCardFunction != AnalogOut) && (stCard.eCardFunction != DigitalOut) && (stCard.eCardFunction != DigitalIO))
        return nSpcMErrorMessageStdOut (&stCard, "Error: Card function not supported by this example\n", false);


    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
    if (!stCard.bSetError)
        {
        printf ("\nPlease select the output mode:\n\n(S)ingleshot\n(C)ontinuous");
        if (stCard.eCardFunction == AnalogOut && !stCard.bM2i)
            printf (" + (D)igital Output (Multi Purpose IO Lines)");
        printf ("\nSingle (R)estart\n");

        switch (cGetch())
            {
            default:
            case 's':
            case 'S':
                vDoCardSetup (&stCard, SPC_REP_STD_SINGLE, 1/*just once*/);
                break;

            case 'c':
            case 'C':
                vDoCardSetup (&stCard, SPC_REP_STD_SINGLE, 0/*forever*/);
                break;

            case 'r':
            case 'R':
                vDoCardSetup (&stCard, SPC_REP_STD_SINGLERESTART);
                break;

            case 'd':
            case 'D':
                vDoCardSetup (&stCard, SPC_REP_STD_SINGLE, 0/*forever*/);
                if (stCard.eCardFunction == AnalogOut && !stCard.bM2i)
                    bMultiPurposeDigOut = true;
                break;
            }
        }


    // ------------------------------------------------------------------------
    // calculate the amount of data we need and allocate memory buffer
    if (!stCard.bSetError)
        {
        if (stCard.eCardFunction == DigitalOut || stCard.eCardFunction == DigitalIO)
            qwMemInBytes = stCard.llSetMemsize * stCard.lSetChannels / 8;
        else
            qwMemInBytes = stCard.llSetMemsize * stCard.lBytesPerSample * stCard.lSetChannels;

        // buffer for data transfer, containing multiplexed data later on
        pvBuffer = pvAllocMemPageAligned (qwMemInBytes);
        if (!pvBuffer)
            return nSpcMErrorMessageStdOut (&stCard, "Memory allocation error\n", false);

        // calculate the data
        if (!bDoDataCalculation (&stCard, pvBuffer))
            return nSpcMErrorMessageStdOut (&stCard, "Data calculation failed\n", false);

        if (bMultiPurposeDigOut)
            {
            if (!bSetMultiPurposeDigOut (&stCard, pvBuffer))
                return nSpcMErrorMessageStdOut (&stCard, "Multi Pupose IO Lines setup failed\n", false);
            }
        }

    // ------------------------------------------------------------------------
    // start the generation
    if (!stCard.bSetError)
        {

        // we define the buffer for transfer and start the DMA transfer
        printf ("Starting the DMA transfer and waiting until data is in board memory\n");
        spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, pvBuffer, 0, qwMemInBytes);
        spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);

        // check for error code
        if (spcm_dwGetErrorInfo_i32 (stCard.hDrv, NULL, NULL, szBuffer))
            {
            vFreeMemPageAligned (pvBuffer, qwMemInBytes);
            return nSpcMErrorMessageStdOut (&stCard, szBuffer, false);
            }
        printf ("... data has been transferred to board memory\n");

        // We'll start and wait untill the card has finished or until a timeout occurs
        spcm_dwSetParam_i32 (stCard.hDrv, SPC_TIMEOUT, 5000);
        printf ("\nStarting the card and waiting for ready interrupt\n(continuous and single restart will have timeout)\n");
        if (spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY) == ERR_TIMEOUT)
            {
            vFreeMemPageAligned (pvBuffer, qwMemInBytes);
            spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_STOP);
            return nSpcMErrorMessageStdOut (&stCard, "... Timeout", false);
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
