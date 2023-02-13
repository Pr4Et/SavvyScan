/*
**************************************************************************

simple_sync_rep_fifo.cpp                                 (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based analog generation cards. 

Information about the different products and their drivers can be found
online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/platform-driver-and-series-differences

Shows a synchronized FIFO mode example using only the few necessary commands
  
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

#include "../common/ostools/spcm_oswrap.h"
#include "../common/ostools/spcm_ostools.h"

// ----- standard c include files -----
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
**************************************************************************
szTypeToName: doing name translation
**************************************************************************
*/

char* szTypeToName (int32 lCardType)
    {
    static char szName[50];
    switch (lCardType & TYP_SERIESMASK)
        {
        case TYP_M2ISERIES:     sprintf (szName, "M2i.%04x",     lCardType & TYP_VERSIONMASK);  break;
        case TYP_M2IEXPSERIES:  sprintf (szName, "M2i.%04x-Exp", lCardType & TYP_VERSIONMASK);  break;
        case TYP_M3ISERIES:     sprintf (szName, "M3i.%04x",     lCardType & TYP_VERSIONMASK);  break;
        case TYP_M3IEXPSERIES:  sprintf (szName, "M3i.%04x-Exp", lCardType & TYP_VERSIONMASK);  break;
        case TYP_M4IEXPSERIES:  sprintf (szName, "M4i.%04x-x8",  lCardType & TYP_VERSIONMASK);  break;
        case TYP_M4XEXPSERIES:  sprintf (szName, "M4x.%04x-x4",  lCardType & TYP_VERSIONMASK);  break;
        case TYP_M2PEXPSERIES:  sprintf (szName, "M2p.%04x-x4",  lCardType & TYP_VERSIONMASK);  break;
        default:                sprintf (szName, "unknown type");                               break;
        }
    return szName;
    }



/*
**************************************************************************
main 
**************************************************************************
*/

int main ()
    {
    drv_handle  ahCard[16];
    drv_handle  hSync = 0;
    bool        bError = false;
    int32       lIdx, lCardCount;
    int32       lCardType, lSerialNumber, lFncType, lFeatures, lStarHubCarrierIdx;
    int16*      apnData[16];
    char        szErrorTextBuffer[ERRORTEXTLEN], szName[50];

    // settings for the FIFO mode buffer handling
    int64       llSamplerate =  MEGA(1);

    // setup for the FIFO mode (HW buffer size can be programmed starting with firmware V9)
    int64       llBufferSize =   MEGA_B(16);
    int64       llHWBufSize =   MEGA_B(32);
    uint32      dwNotifySize =  MEGA_B(1);

    uint32      dwNumPeriodsInBuffer = 128;

    // some example checks
    if (llBufferSize % dwNotifySize)
        {
        printf ("In our example we can only handle sw buffers that are a whole numbered multiple of the notify size\n");
        return 1;
        }

    // ------------------------------------------------------------------------
    // we try to open all cards and printout some information on them
    for (lCardCount = 0; !bError && (lCardCount < 16); lCardCount++)
        {
        // uncomment the second line and replace the IP address to use remote
        // cards like in a generatorNETBOX
        sprintf (szName, "/dev/spcm%d", lCardCount);
        // sprintf (szName, "TCPIP::192.168.1.10::inst%d::INSTR", lCardCount);
        ahCard[lCardCount] = spcm_hOpen (szName);
        apnData[lCardCount] = NULL;

        // not one card found
        if (!lCardCount && !ahCard[lCardCount])
            {
            printf ("no card found...\n");
            return 0;
            }

        // no more cards found in system
        if (!ahCard[lCardCount])
            break;

        // read out some info and print it
        spcm_dwGetParam_i32 (ahCard[lCardCount], SPC_PCITYP,            &lCardType);
        spcm_dwGetParam_i32 (ahCard[lCardCount], SPC_PCISERIALNO,       &lSerialNumber);
        spcm_dwGetParam_i32 (ahCard[lCardCount], SPC_FNCTYPE,           &lFncType);

        // we check which card carries the StarHub
        spcm_dwGetParam_i32 (ahCard[lCardCount], SPC_PCIFEATURES,       &lFeatures);
        if (lFeatures & (SPCM_FEAT_STARHUB4 | SPCM_FEAT_STARHUB16))
            lStarHubCarrierIdx = lCardCount;

        switch (lFncType)
            {
            case SPCM_TYPE_AO:  
                printf ("Found: %s sn %05d\n", szTypeToName (lCardType), lSerialNumber);
                break;

            // only D/A cards allowed in example
            default:
                printf ("Card: %s sn %05d not supported by examplee\n", szTypeToName (lCardType), lSerialNumber);            
                bError = true;
                break;
            }
        }

    int32 lMaxADCValue = 0;
    spcm_dwGetParam_i32 (ahCard[0], SPC_MIINST_MAXADCVALUE, &lMaxADCValue);

    // ------------------------------------------------------------------------
    // try to open the star-hub
    if (!bError)
        {
        hSync = spcm_hOpen ("sync0");
        if (!hSync)
            {
            printf ("no star-hub found. This is essential for the example ...\n");
            bError = true;
            }
        else
            printf ("Found Star-Hub ...\n");
        }


    // ------------------------------------------------------------------------
    // do a simple standard setup for all cards
    for (lIdx = 0; !bError && (lIdx < lCardCount); lIdx++)
        {
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_CHENABLE,        CHANNEL0);              // just 1 channel enabled
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_CARDMODE,        SPC_REP_FIFO_SINGLE);   // single FIFO mode
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_TIMEOUT,         5000);                  // timeout 5 s
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_CLOCKMODE,       SPC_CM_INTPLL);         // clock mode internal PLL
        spcm_dwSetParam_i64 (ahCard[lIdx], SPC_SAMPLERATE,      llSamplerate);          // sampling clock (100k or smallest possible)
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_CLOCKOUT,        0);                     // no clock output
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_AMP0,            1000);                  // 1V output amplitude
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_LOOPS,           0);                     // loop continuous

        // trigger mode definition (in our example it is software trigger)
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_TRIG_ORMASK,     0);                     
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_TRIG_ANDMASK,    0);                     
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_TRIG_ORMASK,     SPC_TMASK_SOFTWARE);
        }

    // starting with firmware version V9 we can program the hardware buffer size to reduce the latency
    for (lIdx = 0; !bError && (lIdx < lCardCount); lIdx++)
        {
        int32 lFWVersion = 0;
        spcm_dwGetParam_i32 (ahCard[lIdx], SPC_PCIVERSION, &lFWVersion);
        if (lFWVersion >= 9)
            {
            spcm_dwSetParam_i64 (ahCard[lIdx], SPC_DATA_OUTBUFSIZE, llHWBufSize);
            spcm_dwSetParam_i32 (ahCard[lIdx], SPC_M2CMD, M2CMD_CARD_WRITESETUP);
            }
        }

    // ------------------------------------------------------------------------
    // define and allocate the data buffers
    for (lIdx = 0; !bError && (lIdx < lCardCount); lIdx++)
        {
        apnData[lIdx] = (int16*) pvAllocMemPageAligned ((uint64) llBufferSize);
        if (!apnData[lIdx])
            {
            printf ("memory allocation failed\n");
            bError = true;
            }
        else
            spcm_dwDefTransfer_i64 (ahCard[lIdx], SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, dwNotifySize, apnData[lIdx], 0, llBufferSize);
        }

    // ------------------------------------------------------------------------
    // calculate data for output (a sine in this case)
    for (lIdx = 0; !bError && (lIdx < lCardCount); lIdx++)
        {
        for (int64 llPos = 0; llPos < llBufferSize; ++llPos)
            {
            apnData[lIdx][llPos / sizeof (int16)] = static_cast < int16 > (lMaxADCValue * sin (((2.0 * 3.141 * llPos) / llBufferSize) * dwNumPeriodsInBuffer));
            }
        }

    // ------------------------------------------------------------------------
    // do the sync setup, start and check for error
    if (!bError)
        {
        spcm_dwSetParam_i32 (hSync, SPC_SYNC_ENABLEMASK, (1 << lCardCount) - 1);

        // clock mask only used for M2i series
        spcm_dwSetParam_i32 (hSync, SPC_SYNC_CLKMASK, (1 << lStarHubCarrierIdx));

        // transfer data in buffer to card memory
        for (lIdx = 0; !bError && (lIdx < lCardCount); lIdx++)
            {
            spcm_dwSetParam_i64 (ahCard[lIdx], SPC_DATA_AVAIL_CARD_LEN, llBufferSize);

            spcm_dwSetParam_i32 (ahCard[lIdx], SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);
            }
        }


    // ------------------------------------------------------------------------
    // this is the main data loop: we wait for all cards to have space available
    // in a loop. This may be put into threads to have one thread controlling data of 
    // one card. In our example we simply put it in the main loop as each card is 
    // configured similar and will have space available at the same time
    if (!bError)
        {
        int64 llTransferredBytes[16] = { 0 };
        int64 llAvailUser          = 0;
        int64 llBufferFillPromille = 0;
        int64 llUserPos            = 0;
        bool  bStarted = false;
        bool  bAbort   = false;

        // run the FIFO loop
        uint32 dwNumCardsFilled = 0; // count number of cards whose memory has been filled (only pre-start)
        while (!bError && !bAbort)
            {
            printf ("\n");
            for (lIdx = 0; !bError && !bAbort && (lIdx < lCardCount); lIdx++)
                {
                // get available space in DMA buffer
                spcm_dwGetParam_i64 (ahCard[lIdx], SPC_DATA_AVAIL_USER_LEN, &llAvailUser);
                spcm_dwGetParam_i64 (ahCard[lIdx], SPC_FILLSIZEPROMILLE,    &llBufferFillPromille);
                printf ("Card: %d   SW-Buffer: %3.0f%%   HW-Buffer:%3d%%, Total Bytes so far: %6.2f MB\n", lIdx, (float) 100.0 * (llBufferSize - llAvailUser) / llBufferSize, (uint32) llBufferFillPromille / 10, (float) llTransferredBytes[lIdx] / MEGA_B(1));

                if (llAvailUser >= dwNotifySize)
                    {
                    // get position of free space in DMA buffer
                    spcm_dwGetParam_i64 (ahCard[lIdx], SPC_DATA_AVAIL_USER_POS, &llUserPos);

                    // avoid buffer wrap-around
                    int64 llDataToWrite = dwNotifySize;
                    if (llUserPos + dwNotifySize > llBufferSize)
                        llDataToWrite = llBufferSize - llUserPos;

                    // calculate new data
                    for (int64 llPos = llUserPos; llPos < llUserPos + llDataToWrite; ++llPos)
                        {
                        apnData[lIdx][llPos / sizeof (int16)] = static_cast < int16 > (lMaxADCValue * sin (((2.0 * 3.141 * llPos) / llBufferSize) * dwNumPeriodsInBuffer));
                        }

                    // set data available for transfer
                    spcm_dwSetParam_i64 (ahCard[lIdx], SPC_DATA_AVAIL_CARD_LEN, llDataToWrite);
                    llTransferredBytes[lIdx] += llDataToWrite;
                    }

                // we start the output as soon as we have a sufficient amount of data on card 
                // inhere we start if the hardware buffer is completely full
                if (!bStarted && !bError && (llBufferFillPromille == 1000))
                    {
                    dwNumCardsFilled++;

                    // when the memory of all cards has been filled, start the cards
                    if (dwNumCardsFilled == lCardCount)
                        {
                        printf ("\nStart the output\n");

                        // start command including enable of trigger engine -> error check
                        if (spcm_dwSetParam_i32 (hSync, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER) != ERR_OK)
                            {
                            spcm_dwGetErrorInfo_i32 (hSync, NULL, NULL, szErrorTextBuffer);
                            printf ("%s\n", szErrorTextBuffer);
                            bError = true;
                            }
                        bStarted = true;
                        }
                    }

                // wait for the next buffer to be free
                switch (spcm_dwSetParam_i32 (ahCard[lIdx], SPC_M2CMD, M2CMD_DATA_WAITDMA))
                    {
                    case ERR_TIMEOUT:
                        printf ("... timeout\n");
                        bError = true;
                        break;

                    case ERR_FIFOHWOVERRUN:
                    case ERR_FIFOBUFOVERRUN:
                        printf ("... buffer underrun\n");
                        bError = true;
                        break;

                    case ERR_OK:
                        break;

                    default:
                        bError = true;
                        break;
                    }

                // check for esape = abort
                if (bKbhit ())
                    {
                    if (cGetch () == 27)
                        bAbort = true;
                    }
                }
            }
        }


    // ------------------------------------------------------------------------
    // check for errors
    if (bError)
        {
        for (lIdx = 0; lIdx < lCardCount; lIdx++)
            spcm_dwGetErrorInfo_i32 (ahCard[lIdx], NULL, NULL, szErrorTextBuffer);
        printf ("Error card %d: %s\n", lIdx, szErrorTextBuffer);
        }

    // no error: we managed to run through completely
    else
        printf ("\nfinished ...\n");

    // send the stop command
    for (lIdx = 0; lIdx < lCardCount; lIdx++)
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_M2CMD, M2CMD_CARD_STOP | M2CMD_DATA_STOPDMA);

    // clean up
    for (lIdx = 0; lIdx < lCardCount; lIdx++)
        if (apnData[lIdx])
            vFreeMemPageAligned (apnData[lIdx], (uint64) llBufferSize);

    if (hSync)
        spcm_vClose (hSync);
    for (lIdx = 0; lIdx < lCardCount; lIdx++)
        spcm_vClose (ahCard[lIdx]);

    return EXIT_SUCCESS;
    }

