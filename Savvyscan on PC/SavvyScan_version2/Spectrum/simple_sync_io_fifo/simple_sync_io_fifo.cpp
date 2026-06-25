/*
**************************************************************************

simple_sync_io_fifo.cpp                                 (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based cards. 

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

#define MAXCARDS 16

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

typedef struct
    {
    drv_handle hCard;
    bool       bDMAStarted;
    int32      lFncType;
    int32      lMaxADCValue;
    int16*     pnData;
    } ST_CARD;

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
bSetupAICard: doing simple setup for AD cards
**************************************************************************
*/

bool bSetupAICard (drv_handle hCard)
    {
    int32 lCardType;

    uint32 dwError = 0;

    dwError += spcm_dwSetParam_i32 (hCard, SPC_CHENABLE,     CHANNEL0);              // just 1 channel enabled
    dwError += spcm_dwSetParam_i32 (hCard, SPC_CARDMODE,     SPC_REC_FIFO_SINGLE);   // single FIFO mode
    dwError += spcm_dwSetParam_i32 (hCard, SPC_PRETRIGGER,   1024);                  // 1k of pretrigger data at start of FIFO mode
    dwError += spcm_dwSetParam_i32 (hCard, SPC_TIMEOUT,      5000);                  // timeout 5 s
    dwError += spcm_dwSetParam_i32 (hCard, SPC_AMP0,         1000);                  // 1V input range 
    dwError += spcm_dwSetParam_i32 (hCard, SPC_CLOCKMODE,    SPC_CM_INTPLL);         // clock mode internal PLL
    dwError += spcm_dwSetParam_i32 (hCard, SPC_CLOCKOUT,     0);                     // no clock output
    dwError += spcm_dwSetParam_i32 (hCard, SPC_TRIG_ORMASK,  0);                     
    dwError += spcm_dwSetParam_i32 (hCard, SPC_TRIG_ANDMASK, 0);                     
    dwError += spcm_dwSetParam_i32 (hCard, SPC_TRIG_ORMASK,  SPC_TMASK_SOFTWARE);

    dwError += spcm_dwGetParam_i32 (hCard, SPC_PCITYP, &lCardType);

    // we try to set the samplerate to 100 kHz (M2i) or 10 MHz
    if (((lCardType & TYP_SERIESMASK) == TYP_M2ISERIES) || ((lCardType & TYP_SERIESMASK) == TYP_M2IEXPSERIES))
        dwError += spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE, KILO(100));
    else
        dwError += spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE, MEGA(10));

    if (!dwError)
        return false;
    else
        return true;
    }

/*
**************************************************************************
bSetupAOCard: doing simple setup for DA cards
**************************************************************************
*/

bool bSetupAOCard (drv_handle hCard)
    {
    int32 lCardType;

    uint32 dwError = 0;

    dwError += spcm_dwSetParam_i32 (hCard, SPC_CHENABLE,     CHANNEL0);              // just 1 channel enabled
    dwError += spcm_dwSetParam_i32 (hCard, SPC_CARDMODE,     SPC_REP_FIFO_SINGLE);   // single FIFO mode
    dwError += spcm_dwSetParam_i32 (hCard, SPC_TIMEOUT,      5000);                  // timeout 5 s
    dwError += spcm_dwSetParam_i32 (hCard, SPC_CLOCKMODE,    SPC_CM_INTPLL);         // clock mode internal PLL
    dwError += spcm_dwSetParam_i32 (hCard, SPC_CLOCKOUT,     0);                     // no clock output
    dwError += spcm_dwSetParam_i32 (hCard, SPC_AMP0,         1000);                  // 1V output amplitude
    dwError += spcm_dwSetParam_i32 (hCard, SPC_LOOPS,        0);                     // loop continuous
    dwError += spcm_dwSetParam_i32 (hCard, SPC_ENABLEOUT0,   1);
    dwError += spcm_dwSetParam_i32 (hCard, SPC_TRIG_ORMASK,  0);                     
    dwError += spcm_dwSetParam_i32 (hCard, SPC_TRIG_ANDMASK, 0);                     
    dwError += spcm_dwSetParam_i32 (hCard, SPC_TRIG_ORMASK,  SPC_TMASK_SOFTWARE);

    dwError += spcm_dwGetParam_i32 (hCard, SPC_PCITYP, &lCardType);

    // we try to set the samplerate to 100 kHz (M2i) or 10 MHz
    if (((lCardType & TYP_SERIESMASK) == TYP_M2ISERIES) || ((lCardType & TYP_SERIESMASK) == TYP_M2IEXPSERIES))
        dwError += spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE, KILO(100));
    else
        dwError += spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE, MEGA(10));

    if (!dwError)
        return false;
    else
        return true;
    }

/*
**************************************************************************
main 
**************************************************************************
*/

int main ()
    {
    ST_CARD astCards[MAXCARDS];
    drv_handle  hSync = 0;
    bool        bError = false;
    int32       lCardType, lSerialNumber, lFncType, lFeatures, lMaxADCValue, lStatus;
    int64       llAvailUser, llPCPos;
    char        szErrorTextBuffer[ERRORTEXTLEN], szName[50];
    
    int32  lCardsCount = 0;
    int32  lFirstAICardIdx = -1;
    bool   bStarHubAvailable = false;

    uint64 qwBytesToTransfer = MEGA_B(256);
    int32  lFillSizePercentStart = 50;
    uint32 dwNumPeriodsInBuffer = 128;

    // settings for the FIFO mode buffer handling
    int64 llBufferSize = MEGA_B(16);
    int32 lNotifySize = KILO_B(128);

    // clear card structs
    for (int32 lIdx = 0; lIdx < MAXCARDS; lIdx++)
        memset (&astCards[lIdx], 0, sizeof (ST_CARD));
        
    // we try to open all cards and printout some information on them
    for (int32 lCardIdx = 0; !bError && (lCardIdx < MAXCARDS); lCardIdx++)
        {
        // uncomment the second line and replace the IP address to use remote cards
        sprintf (szName, "/dev/spcm%d", lCardIdx);
        // sprintf (szName, "TCPIP::192.168.1.10::inst%d::INSTR", lCardIdx);
        drv_handle hCard = spcm_hOpen (szName);

        // not one card found
        if (!lCardIdx && !hCard)
            {
            printf ("no card found...\n");
            return 0;
            }

        // no more cards found in system
        if (!hCard)
            break;

        // read out some info and print it
        spcm_dwGetParam_i32 (hCard, SPC_PCITYP,      &lCardType);
        spcm_dwGetParam_i32 (hCard, SPC_PCISERIALNO, &lSerialNumber);
        spcm_dwGetParam_i32 (hCard, SPC_FNCTYPE,     &lFncType);
        spcm_dwGetParam_i32 (hCard, SPC_MIINST_MAXADCVALUE, &lMaxADCValue);

        // we check if StarHub is available
        spcm_dwGetParam_i32 (hCard, SPC_PCIFEATURES, &lFeatures);
        if (lFeatures & (SPCM_FEAT_STARHUB4 | SPCM_FEAT_STARHUB16))
            bStarHubAvailable = true;
        
        if ((lFncType != SPCM_TYPE_AI) && (lFncType != SPCM_TYPE_AO))
            {
            printf ("Card: %s sn %05d not supported by examplee\n", szTypeToName (lCardType), lSerialNumber);            
            bError = true;
            }
        else
            {
            printf ("Found: %s sn %05d\n", szTypeToName (lCardType), lSerialNumber);

            ST_CARD stCard;
            stCard.hCard = hCard;
            stCard.lFncType = lFncType;
            stCard.lMaxADCValue = lMaxADCValue;

            if ((lFirstAICardIdx < 0) && (lFncType == SPCM_TYPE_AI))
                lFirstAICardIdx = lCardsCount;

            astCards[lCardsCount++] = stCard;
            }
        }

    // try to open the star-hub
    if (!bError && bStarHubAvailable)
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

    // define the data buffer
    for (int32 lCardIdx = 0; !bError && (lCardIdx < lCardsCount); lCardIdx++)
        {
        astCards[lCardIdx].pnData = (int16*) pvAllocMemPageAligned ((uint64) llBufferSize);
        if (!astCards[lCardIdx].pnData)
            {
            printf ("memory allocation failed\n");
            bError = true;
            }

        if (!bError)
            {
            switch (astCards[lCardIdx].lFncType)
                {
                case SPCM_TYPE_AI:
                    bError = bSetupAICard (astCards[lCardIdx].hCard);
                    spcm_dwDefTransfer_i64 (astCards[lCardIdx].hCard, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, lNotifySize, astCards[lCardIdx].pnData, 0, llBufferSize);
                    break;

                case SPCM_TYPE_AO:
                    bError = bSetupAOCard (astCards[lCardIdx].hCard);
                    spcm_dwDefTransfer_i64 (astCards[lCardIdx].hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, lNotifySize, astCards[lCardIdx].pnData, 0, llBufferSize);

                    // calculate data for output (a sine in this case)
                    for (int64 llPos = 0; llPos < llBufferSize; ++llPos)
                        astCards[lCardIdx].pnData[llPos / sizeof (int16)] = static_cast < int16 > (astCards[lCardIdx].lMaxADCValue * sin (((2.0 * 3.141 * llPos) / llBufferSize) * dwNumPeriodsInBuffer));
                    break;
                }
            }
        }

    // do the sync setup
    if (!bError)
        {
        spcm_dwSetParam_i32 (hSync, SPC_SYNC_ENABLEMASK, (1 << lCardsCount) - 1);

        // transfer data in buffer to DA cards memory
        for (int32 lCardIdx = 0; !bError && (lCardIdx < lCardsCount); lCardIdx++)
            {
            astCards[lCardIdx].bDMAStarted = false;

            if (astCards[lCardIdx].lFncType == SPCM_TYPE_AO)
                {
                spcm_dwSetParam_i64 (astCards[lCardIdx].hCard, SPC_DATA_AVAIL_CARD_LEN, llBufferSize);
                spcm_dwSetParam_i32 (astCards[lCardIdx].hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);
                astCards[lCardIdx].bDMAStarted = true;
                }
            }
        }

    // ------------------------------------------------------------------------
    // this is the main data loop: we wait for all cards having data available
    // in a loop. This may be put into threads to have one thread controlling data of 
    // one card.
    uint64 qwTotalMem = 0;
    
    int64  llBufferFillPromille = 0;
    
    bool bStartCards = false;
    bool bStarted = false;

    if (!bError)
        {
        printf ("\n***** Before all cards are started the DA cards buffers are filled to %d%% *****\n", lFillSizePercentStart);

        // run the FIFO loop
        bool bAbort = false;
        while ((qwTotalMem < qwBytesToTransfer) && !bError && !bAbort)
            {
            if (bStartCards)
                {
                printf ("\n\n***** Start Cards *****\n");

                // start DMA transfer for all AD cards
                for (int32 lCardIdx = 0; lCardIdx < lCardsCount; lCardIdx++)
                    {
                    if (astCards[lCardIdx].lFncType == SPCM_TYPE_AI)
                        {
                        spcm_dwSetParam_i32 (astCards[lCardIdx].hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA);
                        astCards[lCardIdx].bDMAStarted = true;
                        }
                    }

                // start command including enable of trigger engine -> error check
                if (spcm_dwSetParam_i32 (hSync, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER) != ERR_OK)
                    {
                    spcm_dwGetErrorInfo_i32 (hSync, NULL, NULL, szErrorTextBuffer);
                    printf ("%s\n", szErrorTextBuffer);
                    bError = true;
                    }
                bStarted = true;
                bStartCards = false;
                }

           
            for (int32 lCardIdx = 0; !bError && (lCardIdx < lCardsCount); lCardIdx++)
                {
                // wait for new data to be available and check return value
                if (astCards[lCardIdx].bDMAStarted)
                    {
                    switch (spcm_dwSetParam_i32 (astCards[lCardIdx].hCard, SPC_M2CMD, M2CMD_DATA_WAITDMA))
                        {
                        case ERR_TIMEOUT:
                            printf ("... timeout\n");
                            bError = true;
                            break;

                        case ERR_FIFOHWOVERRUN:
                        case ERR_FIFOBUFOVERRUN:
                            printf ("... buffer overrun\n");
                            bError = true;
                            break;

                        case ERR_OK:
                            break;

                        default:
                            bError = true;
                            break;
                        }
                    }

                // ***** AD Card *****
                if (bStarted && (astCards[lCardIdx].lFncType == SPCM_TYPE_AI))
                    {
                    // read out status information, available data and do something with the data
                    if (!bError)
                        {
                        spcm_dwGetParam_i32 (astCards[lCardIdx].hCard, SPC_M2STATUS,             &lStatus);
                        spcm_dwGetParam_i64 (astCards[lCardIdx].hCard, SPC_DATA_AVAIL_USER_LEN,  &llAvailUser);
                        spcm_dwGetParam_i64 (astCards[lCardIdx].hCard, SPC_DATA_AVAIL_USER_POS,  &llPCPos);

                        if (llAvailUser >= lNotifySize)
                            {
                            // we count data only on first card
                            if (lCardIdx == lFirstAICardIdx)
                                qwTotalMem += lNotifySize;

                            printf ("AD: Card:%d: Stat:%08x Pos:%016llx Avail:%016llx Total:%.2fMB\r", lCardIdx, lStatus, llPCPos, llAvailUser, (double) (int64) qwTotalMem / MEGA_B(1));

                            // !!!!!!!!!! this is the point to do anything with the data !!!!!!!!!!
                            // in the example we set all data to zero to show how to access the buffer
                            // Please note that SPC_DATA_AVAIL_USER_POS and SPC_DATA_AVAIL_USER_LEN give position 
                            // and length in bytes. To acces int16 array we therefore use lPCPos/2
                            memset ((void*) &astCards[lCardIdx].pnData[llPCPos/2], 0, lNotifySize); 

                            // buffer is free for DMA transfer again
                            spcm_dwSetParam_i32 (astCards[lCardIdx].hCard, SPC_DATA_AVAIL_CARD_LEN,  lNotifySize);
                            }
                        }
                    }

                // ***** DA Card *****
                if (astCards[lCardIdx].lFncType == SPCM_TYPE_AO)
                    {
                    // get available space in DMA buffer
                    spcm_dwGetParam_i64 (astCards[lCardIdx].hCard, SPC_DATA_AVAIL_USER_LEN, &llAvailUser);
                    spcm_dwGetParam_i64 (astCards[lCardIdx].hCard, SPC_FILLSIZEPROMILLE,    &llBufferFillPromille);
                    
                    if (!bStarted)
                        printf ("DA: Card: %d   SW-Buffer: %3.0f%%   HW-Buffer:%3d%%\r", lCardIdx, (float) 100.0 * (llBufferSize - llAvailUser) / llBufferSize, (uint32) llBufferFillPromille / 10);

                    if (llAvailUser >= lNotifySize)
                        {
                        // get position of free space in DMA buffer
                        spcm_dwGetParam_i64 (astCards[lCardIdx].hCard, SPC_DATA_AVAIL_USER_POS, &llPCPos);

                        // avoid buffer wrap-around
                        int64 llDataToWrite = lNotifySize;
                        if (llPCPos + lNotifySize > llBufferSize)
                            llDataToWrite = llBufferSize - llPCPos;

                        // calculate new data
                        for (int64 llPos = llPCPos; llPos < llPCPos + llDataToWrite; ++llPos)
                            astCards[lCardIdx].pnData[llPos / sizeof (int16)] = static_cast < int16 > (lMaxADCValue * sin (((2.0 * 3.141 * llPos) / llBufferSize) * dwNumPeriodsInBuffer));
                                
                        // set data available for transfer
                        spcm_dwSetParam_i64 (astCards[lCardIdx].hCard, SPC_DATA_AVAIL_CARD_LEN, llDataToWrite);
                        }

                    // we start the output as soon as we have a sufficient amount of data on card 
                    // inhere we start if the hardware buffer is 50% filled
                    if (!bStarted && !bError && (llBufferFillPromille >= 10*lFillSizePercentStart)) 
                        bStartCards = true;
                    }
                }
                    
            // check for esape = abort
            if (bKbhit())
                {
                if (cGetch() == 27)
                    {
                    bAbort = true;
                    break;
                    }
                }
            }

        printf ("\n\n");
        }

    // stop card
    spcm_dwSetParam_i32 (hSync, SPC_M2CMD, M2CMD_CARD_STOP);

    // check for errors
    if (bError)
        {
        for (int32 lCardIdx = 0; lCardIdx < lCardsCount; lCardIdx++)
            {
            spcm_dwGetErrorInfo_i32 (astCards[lCardIdx].hCard, NULL, NULL, szErrorTextBuffer);
            printf ("Error card %d: %s\n", lCardIdx, szErrorTextBuffer);
            }
        }
    
    for (int32 lCardIdx = 0; lCardIdx < lCardsCount; lCardIdx++)
        {
        // stop DMA
        spcm_dwSetParam_i32 (astCards[lCardIdx].hCard, SPC_M2CMD, M2CMD_DATA_STOPDMA);

        // clean up
        if (astCards[lCardIdx].pnData)
            vFreeMemPageAligned (astCards[lCardIdx].pnData, (uint64) llBufferSize);
        }

    if (hSync)
        spcm_vClose (hSync);

     for (int32 lCardIdx = 0; lCardIdx < lCardsCount; lCardIdx++)
        spcm_vClose (astCards[lCardIdx].hCard);

    return EXIT_SUCCESS;
    }

