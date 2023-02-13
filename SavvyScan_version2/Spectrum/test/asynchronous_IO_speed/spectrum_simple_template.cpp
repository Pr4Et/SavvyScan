/*
**************************************************************************

spectrum_simple_template.cpp                             (c) Spectrum GmbH

**************************************************************************

Multi-Purpose I/O loop speed test
The tool measures the write and read time of asynchronous I/O commands
on the Multi-Purpose I/O lines of M2p and M4i cards.

Documentation for the API as well as a detailed description of the hardware
can be found in the manual for each device which can be found on our website:
https://www.spectrum-instrumentation.com/en/downloads

Further information can be found online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/knowledge-base-overview

**************************************************************************
*/


// ----- include standard driver header from library -----
#include "../../c_header/dlltyp.h"
#include "../../c_header/regs.h"
#include "../../c_header/spcerr.h"
#include "../../c_header/spcm_drv.h"

#include "../../common/ostools/spcm_oswrap.h"
#include "../../common/ostools/spcm_ostools.h"

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
    drv_handle  hCard;
    int32       lCardType, lSerialNumber;
    char        szName[50];

    // ------------------------------------------------------------------------
    // we try to open one card and printout some information on it

    // uncomment the second line and replace the IP address to use remote
    // cards like in a digitizerNETBOX
    sprintf (szName, "/dev/spcm0");
    // sprintf (szName, "TCPIP::192.168.1.10::inst0::INSTR");
    hCard = spcm_hOpen (szName);

    // not one card found
    if (!hCard)
        {
        printf ("no card found...\n");
        return 0;
        }

    // read out some info and print it
    spcm_dwGetParam_i32 (hCard, SPC_PCITYP,            &lCardType);
    spcm_dwGetParam_i32 (hCard, SPC_PCISERIALNO,       &lSerialNumber);
    printf ("Found: %s sn %05d\n", szTypeToName (lCardType), lSerialNumber);


    // check for cards that don't have multi purpose I/O lines
    switch (lCardType & TYP_SERIESMASK)
        {
        case TYP_M2ISERIES:
        case TYP_M2IEXPSERIES:
        case TYP_M3ISERIES:
        case TYP_M3IEXPSERIES:
            printf("test program needs multi-purpose I/O lines what your card type doesn't have\n");
            exit(0);
        }


    // ------------------------------------------------------------------------
    // Test Loops
    // ------------------------------------------------------------------------


    int32 i, lTmp, lMaxIO = 100000;
    LARGE_INTEGER uStartTimeWrite, uEndTimeWrite, uStartTimeRead, uEndTimeRead, uHighResFreq;
    double dRunTime, dAverageTime;
    uint32 dwErr;

    printf("Start Loop\n");

    QueryPerformanceFrequency(&uHighResFreq);

    // write loop
    dwErr = spcm_dwSetParam_i32(hCard, SPCM_X1_MODE, SPCM_XMODE_ASYNCOUT);
    if (!dwErr)
        {
        QueryPerformanceCounter(&uStartTimeWrite);
        for (i = 0; i < lMaxIO; i++)
            dwErr = spcm_dwSetParam_i32(hCard, SPCM_XX_ASYNCIO, i & 2);
        QueryPerformanceCounter(&uEndTimeWrite);
        }

    // read loop
    dwErr = spcm_dwSetParam_i32(hCard, SPCM_X1_MODE, SPCM_XMODE_ASYNCIN);
    if (!dwErr)
        {
        QueryPerformanceCounter(&uStartTimeRead);
        for (i = 0; i < lMaxIO; i++)
            dwErr = spcm_dwGetParam_i32(hCard, SPCM_XX_ASYNCIO, &lTmp);
        QueryPerformanceCounter(&uEndTimeRead);
        }

    // error check
    if (dwErr)
        {
        uint32 dwErrorReg;
        int32 lErrorValue;
        char szErrorText[200];
        spcm_dwGetErrorInfo_i32(hCard, &dwErrorReg, &lErrorValue, szErrorText);
        printf("Error:\n%s\n\n", szErrorText);
        }

    // print result
    else
        {
        printf("Finished: %d AsyncIO Write and %d Read Commands\n", lMaxIO, lMaxIO);

        dRunTime = (double)(uEndTimeWrite.QuadPart - uStartTimeWrite.QuadPart) / uHighResFreq.QuadPart;
        dAverageTime = (double)dRunTime / (double)lMaxIO;
        printf("Write Run Time:          %.8lf s\n", dRunTime);
        printf("Write Average Time:      %.8lf s\n", dAverageTime);
        printf("Write Average Frequency: %.2lf kHz\n", 1.0 / 1000.0 / dAverageTime);
        printf("\n");

        dRunTime = (double)(uEndTimeRead.QuadPart - uStartTimeRead.QuadPart) / uHighResFreq.QuadPart;
        dAverageTime = (double)dRunTime / (double)lMaxIO;
        printf("Read Run Time:           %.8lf s\n", dRunTime);
        printf("Read Average Time:       %.8lf s\n", dAverageTime);
        printf("Read Average Frequency:  %.2lf kHz\n", 1.0 / 1000.0 / dAverageTime);
        }

    spcm_vClose (hCard);
    
    return EXIT_SUCCESS;
    }

