/*
**************************************************************************

simple_rec_fifo.cpp                                      (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based analog acquisition cards. 

Information about the different products and their drivers can be found
online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/platform-driver-and-series-differences

Shows a simple FIFO mode example using only the few necessary commands
  
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
    drv_handle  hCard;
    int32       lCardType, lSerialNumber, lFncType;
    int16*      pnData;
    char        szErrorTextBuffer[ERRORTEXTLEN];
    uint32      dwError;
    int32       lStatus;
    int64       llAvailUser, llPCPos;
    uint64      qwTotalMem = 0;
    uint64      qwToTransfer = MEGA_B(64);

    // settings for the FIFO mode buffer handling
    int64       llBufferSize =   MEGA_B(4);
    int32       lNotifySize =   KILO_B(16);


    // open card
    // uncomment the second line and replace the IP address to use remote
    // cards like in a digitizerNETBOX
    hCard = spcm_hOpen ("/dev/spcm0");
    // hCard = spcm_hOpen ("TCPIP::192.168.1.10::inst0::INSTR");
    if (!hCard)
        {
        printf ("no card found...\n");
        return 0;
        }


    // read type, function and sn and check for A/D card
    spcm_dwGetParam_i32 (hCard, SPC_PCITYP,         &lCardType);
    spcm_dwGetParam_i32 (hCard, SPC_PCISERIALNO,    &lSerialNumber);
    spcm_dwGetParam_i32 (hCard, SPC_FNCTYPE,        &lFncType);

    switch (lFncType)
        {
        case SPCM_TYPE_AI:  
            printf ("Found: %s sn %05d\n", szTypeToName (lCardType), lSerialNumber);
            break;

        default:
            printf ("Card: %s sn %05d not supported by example\n", szTypeToName (lCardType), lSerialNumber);            
            return 0;
        }


    // do a simple standard setup
    spcm_dwSetParam_i32 (hCard, SPC_CHENABLE,       CHANNEL0);              // just 1 channel enabled
    spcm_dwSetParam_i32 (hCard, SPC_PRETRIGGER,     1024);                  // 1k of pretrigger data at start of FIFO mode
    spcm_dwSetParam_i32 (hCard, SPC_CARDMODE,       SPC_REC_FIFO_SINGLE);   // single FIFO mode
    spcm_dwSetParam_i32 (hCard, SPC_TIMEOUT,        5000);                  // timeout 5 s
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_ORMASK,    SPC_TMASK_SOFTWARE);    // trigger set to software
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_ANDMASK,   0);                     // ...
    spcm_dwSetParam_i32 (hCard, SPC_CLOCKMODE,      SPC_CM_INTPLL);         // clock mode internal PLL

    // we try to set the samplerate to 100 kHz (M2i) or 20 MHz (M3i and M4i) on internal PLL, no clock output
    if (((lCardType & TYP_SERIESMASK) == TYP_M2ISERIES) || ((lCardType & TYP_SERIESMASK) == TYP_M2IEXPSERIES))
        spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE, KILO(100));
    if ((lCardType & TYP_SERIESMASK) == TYP_M2PEXPSERIES)
        spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE, MEGA(10));
    else
        spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE, MEGA(20));

    spcm_dwSetParam_i32 (hCard, SPC_CLOCKOUT,       0);                     // no clock output


    // define the data buffer
    pnData = (int16*) pvAllocMemPageAligned ((uint64) llBufferSize);
    if (!pnData)
        {
        printf ("memory allocation failed\n");
        spcm_vClose (hCard);
        return 0;
        }

    spcm_dwDefTransfer_i64 (hCard, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, lNotifySize, pnData, 0, llBufferSize);

    // start everything
    dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_DATA_STARTDMA);


    // check for error
    if (dwError != ERR_OK)
        {
        spcm_dwGetErrorInfo_i32 (hCard, NULL, NULL, szErrorTextBuffer);
        printf ("%s\n", szErrorTextBuffer);
        vFreeMemPageAligned (pnData, (uint64) llBufferSize);
        spcm_vClose (hCard);
        return 0;
        }


    // run the FIFO mode and loop through the data
    else
        {
        while (qwTotalMem < qwToTransfer)
            {
            if ((dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_DATA_WAITDMA)) != ERR_OK)
                {
                if (dwError == ERR_TIMEOUT)
                    printf ("... Timeout\n");
                else
                    printf ("... Error: %d\n", dwError);
                break;
                }

            else
                {
                spcm_dwGetParam_i32 (hCard, SPC_M2STATUS,             &lStatus);
                spcm_dwGetParam_i64 (hCard, SPC_DATA_AVAIL_USER_LEN,  &llAvailUser);
                spcm_dwGetParam_i64 (hCard, SPC_DATA_AVAIL_USER_POS,  &llPCPos);

                if (llAvailUser >= lNotifySize)
                    {
                    qwTotalMem += lNotifySize;
                    printf ("Stat:%08x Pos:%016llx Avail:%016llx Total:%.2fMB\n", lStatus, llPCPos, llAvailUser, (double) (int64) qwTotalMem / MEGA_B(1));

                    // this is the point to do anything with the data

                    spcm_dwSetParam_i32 (hCard, SPC_DATA_AVAIL_CARD_LEN,  lNotifySize);
                    }

                // check for esape = abort
                if (bKbhit ())
                    if (cGetch () == 27)
                        break;
                }
            }
        }


    // send the stop command
    dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_STOP | M2CMD_DATA_STOPDMA);

    // clean up
    printf ("Finished...\n");
    vFreeMemPageAligned (pnData, (uint64) llBufferSize);
    spcm_vClose (hCard);

    return EXIT_SUCCESS;
    }

