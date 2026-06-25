/*
**************************************************************************

simple_sync_rec_fifo.cpp                                 (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based analog acquisition cards. 

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
#include <string.h>
#include <stdlib.h>


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
        case TYP_M2ISERIES:     sprintf (szName, "M2i.%04x", lCardType & TYP_VERSIONMASK);      break;
        case TYP_M2IEXPSERIES:  sprintf (szName, "M2i.%04x-Exp", lCardType & TYP_VERSIONMASK);  break;
        case TYP_M3ISERIES:     sprintf (szName, "M3i.%04x", lCardType & TYP_VERSIONMASK);      break;
        case TYP_M3IEXPSERIES:  sprintf (szName, "M3i.%04x-Exp", lCardType & TYP_VERSIONMASK);  break;
        case TYP_M4IEXPSERIES:  sprintf (szName, "M4i.%04x-x8", lCardType & TYP_VERSIONMASK);   break;
        case TYP_M4XEXPSERIES:  sprintf (szName, "M4x.%04x-x4", lCardType & TYP_VERSIONMASK);   break;
        case TYP_M2PEXPSERIES:  sprintf (szName, "M2p.%04x-x4", lCardType & TYP_VERSIONMASK);   break;
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
    int32       lCardType, lSerialNumber, lFncType, lFeatures, lStarHubCarrierIdx = 0;
    int64       llTmp;
    int16*      apnData[16];
    char        szErrorTextBuffer[ERRORTEXTLEN], szName[50];
    int32       lStatus;
    int64       llAvailUser, llPCPos;
    uint64      qwTotalMem = 0;
    uint64      qwToTransfer = MEGA_B(512);

    // settings for the FIFO mode buffer handling
    int64       llBufferSize =   MEGA_B(16);
    int32       lNotifySize =   KILO_B(128);
    int64       llSamplerate =  KILO(100);

    // ------------------------------------------------------------------------
    // we try to open all cards and printout some information on them
    for (lCardCount = 0; !bError && (lCardCount < 16); lCardCount++)
        {
        // uncomment the second line and replace the IP address to use remote
        // cards like in a digitizerNETBOX
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

        // we check for minimum samplerate (on first series of M3i sync only works with min samplerate * 2)
        spcm_dwGetParam_i64 (ahCard[lCardCount], SPC_MIINST_MINADCLOCK, &llTmp);
        if ((2*llTmp) > llSamplerate)
            llSamplerate = (2*llTmp);

        // we check which card carries the StarHub
        spcm_dwGetParam_i32 (ahCard[lCardCount], SPC_PCIFEATURES,       &lFeatures);
        if (lFeatures & (SPCM_FEAT_STARHUB4 | SPCM_FEAT_STARHUB16))
            lStarHubCarrierIdx = lCardCount;
        
        switch (lFncType)
            {
            case SPCM_TYPE_AI:  
                printf ("Found: %s sn %05d\n", szTypeToName (lCardType), lSerialNumber);
                break;

            // only A/D cards allowed in example
            default:
                printf ("Card: %s sn %05d not supported by examplee\n", szTypeToName (lCardType), lSerialNumber);            
                bError = true;
                break;
            }
        }


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
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_PRETRIGGER,      1024);                  // 1k of pretrigger data at start of FIFO mode
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_CARDMODE,        SPC_REC_FIFO_SINGLE);   // single FIFO mode
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_TIMEOUT,         5000);                  // timeout 5 s
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_CLOCKMODE,       SPC_CM_INTPLL);         // clock mode internal PLL
        spcm_dwSetParam_i64 (ahCard[lIdx], SPC_SAMPLERATE,      llSamplerate);          // sampling clock (100k or smallest possible)
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_CLOCKOUT,        0);                     // no clock output

        // trigger mode definition (in our example it is software trigger)
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_TRIG_ORMASK,     0);                     
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_TRIG_ANDMASK,    0);                     
        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_TRIG_ORMASK,     SPC_TMASK_SOFTWARE);       
        }

    // ------------------------------------------------------------------------
    // define and allocate the data buffers
    for (lIdx = 0; !bError && (lIdx < lCardCount); lIdx++)
        {
        void* pvContBuf     = NULL;
        uint64 qwContBufLen = 0;
        spcm_dwGetContBuf_i64 (ahCard[lIdx], SPCM_BUF_DATA, &pvContBuf, &qwContBufLen);
        if (qwContBufLen >= llBufferSize)
            { 
            printf ("using continuous memory\n");
            apnData[lIdx] = (int16*)pvContBuf;
            }
        else
            apnData[lIdx] = (int16*) pvAllocMemPageAligned ((uint64) llBufferSize);

        if (!apnData[lIdx])
            {
            printf ("memory allocation failed\n");
            bError = true;
            }
        else
            spcm_dwDefTransfer_i64 (ahCard[lIdx], SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, lNotifySize, apnData[lIdx], 0, llBufferSize);
        }


    // ------------------------------------------------------------------------
    // do the sync setup, start and check for error
    if (!bError)
        {
        spcm_dwSetParam_i32 (hSync, SPC_SYNC_ENABLEMASK, (1 << lCardCount) - 1);

        // clock mask only used for M2i series. Ignored on M3i series
        spcm_dwSetParam_i32 (hSync, SPC_SYNC_CLKMASK, (1 << lStarHubCarrierIdx));

        // start DMA on each card (prior to card start to have DMA armed if first data is coming)
        for (lIdx = 0; !bError && (lIdx < lCardCount); lIdx++)
            spcm_dwSetParam_i32 (ahCard[lIdx], SPC_M2CMD, M2CMD_DATA_STARTDMA);

        // start command including enable of trigger engine -> error check
        if (spcm_dwSetParam_i32 (hSync, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER) != ERR_OK)
            {
            spcm_dwGetErrorInfo_i32 (hSync, NULL, NULL, szErrorTextBuffer);
            printf ("%s\n", szErrorTextBuffer);
            bError = true;
            }
        }


    // ------------------------------------------------------------------------
    // this is the main data loop: we wait for all cards having data available
    // in a loop. This may be put into threads to have one thread controlling data of 
    // one card. In our example we simply put it in the main loop as each card is 
    // configured similar and will have data available at the same time
    if (!bError)
        {

        // run the FIFO loop
        bool bAbort = false;
        while ((qwTotalMem < qwToTransfer) && !bError && !bAbort)
            {
            for (lIdx = 0; !bError && (lIdx < lCardCount); lIdx++)
                {

                // wait for new data to be available and check return value
                switch (spcm_dwSetParam_i32 (ahCard[lIdx], SPC_M2CMD, M2CMD_DATA_WAITDMA))
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

                // read out status information, available data and do something with the data
                if (!bError)
                    {
                    spcm_dwGetParam_i32 (ahCard[lIdx], SPC_M2STATUS,             &lStatus);
                    spcm_dwGetParam_i64 (ahCard[lIdx], SPC_DATA_AVAIL_USER_LEN,  &llAvailUser);
                    spcm_dwGetParam_i64 (ahCard[lIdx], SPC_DATA_AVAIL_USER_POS,  &llPCPos);

                    if (llAvailUser >= lNotifySize)
                        {
                        // we count data only on first card
                        if (lIdx == 0)
                            qwTotalMem += lNotifySize;

                        printf ("Card:%d: Stat:%08x Pos:%016llx Avail:%016llx Total:%.2fMB\n", lIdx, lStatus, llPCPos, llAvailUser, (double) (int64) qwTotalMem / MEGA_B(1));

                        // !!!!!!!!!! this is the point to do anything with the data !!!!!!!!!!
                        // in the example we set all data to zero to show how to access the buffer
                        // Please note that SPC_DATA_AVAIL_USER_POS and SPC_DATA_AVAIL_USER_LEN give position 
                        // and length in bytes. To acces int16 array we therefore use lPCPos/2
                        memset ((void*) &apnData[lIdx][llPCPos/2], 0, lNotifySize); 

                        // buffer is free for DMA transfer again
                        spcm_dwSetParam_i32 (ahCard[lIdx], SPC_DATA_AVAIL_CARD_LEN,  lNotifySize);
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

