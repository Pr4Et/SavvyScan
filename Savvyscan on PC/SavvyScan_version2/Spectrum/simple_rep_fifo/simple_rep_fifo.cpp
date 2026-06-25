/*
**************************************************************************

simple_rep_fifo.cpp                                      (c) Spectrum GmbH

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
    drv_handle  hCard = NULL_HANDLE;
    bool        bError = false;
    int32       lCardType, lSerialNumber, lFncType, lFeatures;
    int16*      anData = NULL;
    char        szErrorTextBuffer[ERRORTEXTLEN], szName[50];

    // settings for the FIFO mode buffer handling
    int64       llSamplerate =  MEGA(50);

    // setup for the FIFO mode (HW buffer size can be programmed starting with firmware V9)
    int64       llBufferSize =  MEGA_B(16);
    int64       llHWBufSize =   MEGA_B(256);
    uint32      dwNotifySize =  MEGA_B(1);

    uint32      dwNumPeriodsInBuffer = 128;

    // some example checks
    if (llBufferSize % dwNotifySize)
        {
        printf ("In our example we can only handle sw buffers that are a whole numbered multiple of the notify size\n");
        return 1;
        }

    // ------------------------------------------------------------------------
    // uncomment the second line and replace the IP address to use remote
    // cards like in a generatorNETBOX
    int32 lCardIdx = 0;
    sprintf (szName, "/dev/spcm%d", lCardIdx);
    // sprintf (szName, "TCPIP::192.168.1.10::inst%d::INSTR", lCardCount);
    hCard = spcm_hOpen (szName);

    // no card found
    if (!hCard)
    {
    printf ("no card found...\n");
    return 0;
    }

    // read out some info and print it
    spcm_dwGetParam_i32 (hCard, SPC_PCITYP,            &lCardType);
    spcm_dwGetParam_i32 (hCard, SPC_PCISERIALNO,       &lSerialNumber);
    spcm_dwGetParam_i32 (hCard, SPC_FNCTYPE,           &lFncType);

    switch (lFncType)
    {
    case SPCM_TYPE_AO:  
        printf ("Found: %s sn %05d\n", szTypeToName (lCardType), lSerialNumber);
        break;

    // only D/A cards allowed in example
    default:
        printf ("Card: %s sn %05d not supported by example\n", szTypeToName (lCardType), lSerialNumber);            
        bError = true;
        break;
    }

    int32 lMaxADCValue = 0;
    spcm_dwGetParam_i32 (hCard, SPC_MIINST_MAXADCVALUE, &lMaxADCValue);

    // ------------------------------------------------------------------------
    // do a simple standard setup for all cards
    spcm_dwSetParam_i32 (hCard, SPC_CHENABLE,        CHANNEL0);              // just 1 channel enabled
    spcm_dwSetParam_i32 (hCard, SPC_CARDMODE,        SPC_REP_FIFO_SINGLE);   // single FIFO mode
    spcm_dwSetParam_i32 (hCard, SPC_TIMEOUT,         5000);                  // timeout 5 s
    spcm_dwSetParam_i32 (hCard, SPC_CLOCKMODE,       SPC_CM_INTPLL);         // clock mode internal PLL
    spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE,      llSamplerate);          // sampling clock (100k or smallest possible)
    spcm_dwSetParam_i32 (hCard, SPC_CLOCKOUT,        0);                     // no clock output
    spcm_dwSetParam_i32 (hCard, SPC_AMP0,            1000);                  // 1V output amplitude
    spcm_dwSetParam_i32 (hCard, SPC_LOOPS,           0);                     // loop continuous

    // trigger mode definition (in our example it is software trigger)
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_ORMASK,     0);                     
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_ANDMASK,    0);                     
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_ORMASK,     SPC_TMASK_SOFTWARE);

    // starting with firmware version V9 we can program the hardware buffer size to reduce the latency
    int32 lFWVersion = 0;
    spcm_dwGetParam_i32 (hCard, SPC_PCIVERSION, &lFWVersion);
    if (lFWVersion >= 9)
        {
        spcm_dwSetParam_i64 (hCard, SPC_DATA_OUTBUFSIZE, llHWBufSize);
        spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_WRITESETUP);
        }

    // ------------------------------------------------------------------------
    // define and allocate the data buffers
    anData = (int16*) pvAllocMemPageAligned ((uint64) llBufferSize);
    if (!anData)
        {
        printf ("memory allocation failed\n");
        bError = true;
        }
    else
        {
        spcm_dwDefTransfer_i64 (hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, dwNotifySize, anData, 0, llBufferSize);

        // ------------------------------------------------------------------------
        // calculate data for output (a sine in this case)
        for (int64 llPos = 0; llPos < llBufferSize; ++llPos)
            {
            anData[llPos / sizeof (int16)] = static_cast < int16 > (lMaxADCValue * sin (((2.0 * 3.141 * llPos) / llBufferSize) * dwNumPeriodsInBuffer));
            }
        }

    // ------------------------------------------------------------------------
    // do start and check for error
    if (!bError)
        {
        // transfer data in buffer to card memory
        spcm_dwSetParam_i64 (hCard, SPC_DATA_AVAIL_CARD_LEN, llBufferSize);

        spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);
        }


    // ------------------------------------------------------------------------
    // this is the main data loop: we wait for all cards to have space available
    // in a loop. This may be put into threads to have one thread controlling data of 
    // one card. In our example we simply put it in the main loop as each card is 
    // configured similar and will have space available at the same time
    if (!bError)
        {
        int64 llTransferredBytes   = 0;
        int64 llAvailUser          = 0;
        int64 llBufferFillPromille = 0;
        int64 llUserPos            = 0;
        bool  bStarted = false;
        bool  bAbort   = false;

        // run the FIFO loop
        uint32 dwNumCardsFilled = 0; // count number of cards whose memory has been filled (only pre-start)
        while (!bError && !bAbort)
            {
            // get available space in DMA buffer
            spcm_dwGetParam_i64 (hCard, SPC_DATA_AVAIL_USER_LEN, &llAvailUser);
            spcm_dwGetParam_i64 (hCard, SPC_FILLSIZEPROMILLE,    &llBufferFillPromille);
            printf ("\rSW-Buffer: %3.0f%%   HW-Buffer:%3d%%, Total Bytes so far: %6.2f MB", (float) 100.0 * (llBufferSize - llAvailUser) / llBufferSize, (uint32) llBufferFillPromille / 10, (float) llTransferredBytes / MEGA_B(1));

            if (llAvailUser >= dwNotifySize)
                {
                // get position of free space in DMA buffer
                spcm_dwGetParam_i64 (hCard, SPC_DATA_AVAIL_USER_POS, &llUserPos);

                // avoid buffer wrap-around
                int64 llDataToWrite = dwNotifySize;
                if (llUserPos + dwNotifySize > llBufferSize)
                    llDataToWrite = llBufferSize - llUserPos;

                // calculate new data
                for (int64 llPos = llUserPos; llPos < llUserPos + llDataToWrite; ++llPos)
                    {
                    anData[llPos / sizeof (int16)] = static_cast < int16 > (lMaxADCValue * sin (((2.0 * 3.141 * llPos) / llBufferSize) * dwNumPeriodsInBuffer));
                    }

                // set data available for transfer
                spcm_dwSetParam_i64 (hCard, SPC_DATA_AVAIL_CARD_LEN, llDataToWrite);
                llTransferredBytes += llDataToWrite;
                }

            // we start the output as soon as we have a sufficient amount of data on card 
            // inhere we start if the hardware buffer is completely full
            if (!bStarted && !bError && (llBufferFillPromille == 1000))
                {
                printf ("\nStart the output\n");

                // start command including enable of trigger engine -> error check
                if (spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER) != ERR_OK)
                    {
                    spcm_dwGetErrorInfo_i32 (hCard, NULL, NULL, szErrorTextBuffer);
                    printf ("%s\n", szErrorTextBuffer);
                    bError = true;
                    }
                bStarted = true;
                }

            // wait for the next buffer to be free
            switch (spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_DATA_WAITDMA))
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

            /*
            // check for escape = abort
            if (bKbhit ())
                {
                if (cGetch () == 27)
                    bAbort = true;
                }
            */
            }
        }


    // ------------------------------------------------------------------------
    // check for errors
    if (bError)
        {
        spcm_dwGetErrorInfo_i32 (hCard, NULL, NULL, szErrorTextBuffer);
        printf ("Error: %s\n", szErrorTextBuffer);
        }

    // no error: we managed to run through completely
    else
        printf ("\nfinished ...\n");

    // send the stop command
    spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_STOP | M2CMD_DATA_STOPDMA);

    // clean up
    if (anData)
        vFreeMemPageAligned (anData, (uint64) llBufferSize);

    spcm_vClose (hCard);

    return EXIT_SUCCESS;
    }

