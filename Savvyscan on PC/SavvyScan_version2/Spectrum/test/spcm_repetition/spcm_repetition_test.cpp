/*
**************************************************************************

spcm_repetion_test.cpp                                   (c) Spectrum GmbH

**************************************************************************

This example supports all M2i cards. It measures the repetion rate for 
different block sizes in either direction

Documentation for the API as well as a detailed description of the hardware
can be found in the manual for each device which can be found on our website:
https://www.spectrum-instrumentation.com/en/downloads

Further information can be found online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/knowledge-base-overview

**************************************************************************
*/



// ----- standard c include files -----
#include <stdio.h>
#include <stdlib.h>
#include <conio.h>

// ----- include driver libraries -----
#include "../../c_header/dlltyp.h"
#include "../../c_header/regs.h"
#include "../../c_header/spcerr.h"
#include "../../c_header/spcm_drv.h"

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

int main ()
    {
    drv_handle  hDrv;
    char        szErrorText[ERRORTEXTLEN];
    int32       lKernelVersion, lDLLVersion, lCardType, lSN, lBytesPerSample, lFunctionType;
    int64       llMaxSamplerate;


    // try to open the card
    printf ("Trying to open the driver ...\n");
    hDrv = spcm_hOpen ("spcm0");

    if (!hDrv)
        {
        printf ("open of device failed\n");
        spcm_dwGetErrorInfo_i32 (NULL, NULL, NULL, szErrorText);
        printf ("... %s\n", szErrorText);
        _getch();
        return EXIT_FAILURE;
        }

    // get some card information and print it
    spcm_dwGetParam_i32 (hDrv, SPC_GETKERNELVERSION,        &lKernelVersion);
    spcm_dwGetParam_i32 (hDrv, SPC_GETDRVVERSION,           &lDLLVersion);
    spcm_dwGetParam_i32 (hDrv, SPC_PCITYP,                  &lCardType);
    spcm_dwGetParam_i32 (hDrv, SPC_FNCTYPE,                 &lFunctionType);    
    spcm_dwGetParam_i32 (hDrv, SPC_PCISERIALNO,             &lSN);
    spcm_dwGetParam_i64 (hDrv, SPC_PCISAMPLERATE,           &llMaxSamplerate);
    spcm_dwGetParam_i32 (hDrv, SPC_MIINST_BYTESPERSAMPLE,   &lBytesPerSample);


    printf ("---------------------------------------\n");
#ifdef WIN32
    printf ("PC: %s\n", getenv ("COMPUTERNAME"));
#endif // WIN32
    printf ("%s sn %05d\n", szTypeToName (lCardType), lSN);
    printf ("Max. Sampling rate: %lld\n", llMaxSamplerate);
    printf ("Kernel Version:     V %d.%d build %d\n", (lKernelVersion >> 24) & 0xff, (lKernelVersion >> 16) & 0xff, lKernelVersion & 0xffff);
    printf ("Library Version:    V %d.%d build %d\n", (lDLLVersion >> 24) & 0xff, (lDLLVersion >> 16) & 0xff, lDLLVersion & 0xffff);

    switch (lCardType & TYP_SERIESMASK)
        {
        case TYP_M2IEXPSERIES:
        case TYP_M3IEXPSERIES:
        case TYP_M4IEXPSERIES:
        case TYP_M4XEXPSERIES:
        case TYP_M2PEXPSERIES:
            {
            int32 lPCIeGen = 0;
            int32 lPCIeLanes = 0;
            int32 lPCIePayload = 0;
            spcm_dwGetParam_i32 (hDrv, SPC_PCIEXPGENERATION, &lPCIeGen);
            spcm_dwGetParam_i32 (hDrv, SPC_PCIEXPLANES,      &lPCIeLanes);
            spcm_dwGetParam_i32 (hDrv, SPC_PCIEXPPAYLOAD,    &lPCIePayload);
            printf ("PCIe Gen%dx%d, Payload: %d Bytes\n", lPCIeGen, lPCIeLanes, lPCIePayload);
            break;
            }
        default:
            // keine Businfos für PCI Karten
            break;
        }
    printf ("---------------------------------------\n");

    // setup the card for test
    int16           *pnBuffer, *pnContMem;
    int32           lMaxBufsize = MEGA_B(64);
    int32           lMaxLoop;
    int32           lBufsize, lLoop;
    uint64          qwContBufLen;
    LARGE_INTEGER   uStart, uEnd, uLoopStart, uLoopEnd, uMin, uMax;
    double          dAvTime, dMinTime, dMaxTime;
    LARGE_INTEGER   uHighResFreq;

    spcm_dwSetParam_i32 (hDrv, SPC_M2CMD,           M2CMD_CARD_RESET);
    spcm_dwSetParam_i32 (hDrv, SPC_TRIG_ANDMASK,    0);
    spcm_dwSetParam_i32 (hDrv, SPC_TRIG_ORMASK,     SPC_TMASK_SOFTWARE);
    spcm_dwSetParam_i32 (hDrv, SPC_TRIG_EXT0_MODE,  SPC_TM_POS);
    spcm_dwSetParam_i32 (hDrv, SPC_CLOCKMODE,       SPC_CM_INTPLL);
    spcm_dwSetParam_i64 (hDrv, SPC_SAMPLERATE,      llMaxSamplerate);
    spcm_dwSetParam_i32 (hDrv, SPC_TIMEOUT,         0);
    spcm_dwSetParam_i32 (hDrv, SPC_CARDMODE,        SPC_REC_STD_SINGLE);
    switch (lFunctionType)
        {
        case SPCM_TYPE_AI:
            spcm_dwSetParam_i32 (hDrv, SPC_CHENABLE,        CHANNEL0);
            spcm_dwSetParam_i32 (hDrv, SPC_CARDMODE,        SPC_REC_STD_SINGLE);
            break;

        case SPCM_TYPE_AO:
        case SPCM_TYPE_DO:
            spcm_vClose (hDrv);
            printf ("\ncard type not supported yet (press key)\n");
            _getch();
            break;

        case SPCM_TYPE_DI:
        case SPCM_TYPE_DIO:
            spcm_dwSetParam_i32 (hDrv, SPC_CHENABLE,        0xff);
            spcm_dwSetParam_i32 (hDrv, SPC_CARDMODE,        SPC_REC_STD_SINGLE);
            break;
        }

    // ----- Buffer Handling -----
    pnBuffer = (int16*) VirtualAlloc (NULL, lMaxBufsize, MEM_COMMIT, PAGE_READWRITE);
    spcm_dwGetContBuf_i64 (hDrv, SPCM_BUF_DATA, (void**) &pnContMem, &qwContBufLen);
    printf ("Cont Mem: %.2f MByte   \n", (double) qwContBufLen / 1024.0 / 1024.0);
    printf ("\n");

    QueryPerformanceFrequency (&uHighResFreq);
    for (lBufsize = 1024; lBufsize <= lMaxBufsize; lBufsize *= 2)
        {
        if (lBufsize < (1024*1024))
            printf ("%4d kByte ", lBufsize / 1024);
        else
            printf ("%4d MByte ", lBufsize / 1024 / 1024);

        // calc max loops depending on the buffersize
        lMaxLoop = lMaxBufsize / lBufsize;
        if (lMaxLoop < 3) lMaxLoop = 3;
        if (lMaxLoop > 1000) lMaxLoop = 1000;

        spcm_dwSetParam_i32 (hDrv, SPC_MEMSIZE,     lBufsize/lBytesPerSample);
        spcm_dwSetParam_i32 (hDrv, SPC_POSTTRIGGER, lBufsize/lBytesPerSample - lBufsize/8);
        spcm_dwSetParam_i32 (hDrv, SPC_M2CMD,       M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY);

        uMin.QuadPart = 1000000000000;
        uMax.QuadPart = 0;

        // now the speed measuring loop starts
        QueryPerformanceCounter (&uStart);
        for (lLoop = 1; lLoop <= lMaxLoop; lLoop++)
            {
            QueryPerformanceCounter (&uLoopStart);

            // use cont mem buf if available and large enough
            if (lBufsize <= (int32) qwContBufLen)
                spcm_dwDefTransfer_i64 (hDrv, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, 0, (void*) pnContMem, 0, lBufsize);
            else
                spcm_dwDefTransfer_i64 (hDrv, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, 0, (void*) pnBuffer, 0, lBufsize);

            spcm_dwSetParam_i32 (hDrv, SPC_M2CMD,     M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);

            // check for min and max repitition time
            QueryPerformanceCounter (&uLoopEnd);
            if ((uLoopEnd.QuadPart - uLoopStart.QuadPart) > uMax.QuadPart)
                uMax.QuadPart = (uLoopEnd.QuadPart - uLoopStart.QuadPart);
            if ((uLoopEnd.QuadPart - uLoopStart.QuadPart) < uMin.QuadPart)
                uMin.QuadPart = (uLoopEnd.QuadPart - uLoopStart.QuadPart);
            }
        QueryPerformanceCounter (&uEnd);


        // calc and print the results for this buffer size
        dAvTime =  (double) (uEnd.QuadPart - uStart.QuadPart) / lMaxLoop / uHighResFreq.QuadPart;
        dMinTime = (double) (uMin.QuadPart) / uHighResFreq.QuadPart;
        dMaxTime = (double) (uMax.QuadPart) / uHighResFreq.QuadPart;

        printf (" %5.1f MB/s Rep: %6.1f Hz Time(ms) Min=%6.2f Max=%6.2f Av=%6.2f\n", 
            (double) (dAvTime == 0 ? 0.0 : lBufsize / dAvTime / 1024.0 / 1024.0),
            (double) (1.0 / dAvTime),
            1000.0 * dMinTime,
            1000.0 * dMaxTime,
            1000.0 * dAvTime);
        }

    VirtualFree (pnBuffer, 0, MEM_RELEASE);
    spcm_vClose (hDrv);
    printf ("\n...end (press key)\n");
    _getch();

    return EXIT_SUCCESS;
    }

