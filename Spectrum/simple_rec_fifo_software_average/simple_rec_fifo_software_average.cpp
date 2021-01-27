/*
**************************************************************************

simple_rec_fifo_software_average                         (c) Spectrum GmbH

**************************************************************************

Example for Spectrum M4i analog acquisition cards to show software based
block average. The example is build for a M4i.2230-x8, a 1 channel
5 GS/s 8 Bit digitizer.

The test parameter section allows to define the different test settings
and handles multi-thread averaging as well as single-thread averaging
For running this example an external trigger source in the region of
1 to 2 kHz is needed.

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
#include <string.h>


/*
**************************************************************************
Test Parameters
**************************************************************************
*/

const bool        bThreads =      true;               // test the thread based version
const int32       lThreads =      4;                  // number of threads to split the average to
const int32       lSegmentsize =  KILO_B(1024);       // segment size per trigger
const int32       lNumSegments =  1;                  // number of segments per interrupt
const int32       lAverageLoop =  1000;               // number of averages (summations) 




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
bKeyCheckAsync: faster than kbhit()
**************************************************************************
*/
int g_nKeyPress;
bool bKeyCheckAsync ()
    {
#ifdef WIN32
    return (g_nKeyPress != GetAsyncKeyState(VK_ESCAPE));
#endif
    }


/*
**************************************************************************
thread structúre and thread for average
**************************************************************************
*/

typedef struct
    {
    bool                bRunNotQuit;        // controls whether the thread is running or should stop
    int8*               pbyData;            // pointer to the data buffer
    int64               llCurDataPos;       // current position of next available data inside the buffer
    int32*              plAverageData;      // pointer to average (summation) buffer
    int32               lSegmentsize;       // size of complete average segment
    int32               lStartOffset;       // start position inside segment that the thread is handling
    int32               lAveragesize;       // size of average data that is handled by this thread
    int32               lNumSegments;       // number of sequencing segments that should be handled in one run
    SPCM_EVENT_HANDLE   hStart;             // received event to start calculation
    SPCM_EVENT_HANDLE   hEnd;               // transmitted event after end of calculation
    SPCM_THREAD_HANDLE  hThread;            // thread handle
    } SPCM_AVERAGE_DATA;

// ***********************************************************************

SPCM_THREAD_RETURN SPCM_THREAD_CALLTYPE pvAverageSegmentPart (void* pvArguments)
    {
    SPCM_AVERAGE_DATA* pstData = (SPCM_AVERAGE_DATA*) pvArguments;
    int32  i, j;
    int32  lStart;
    int32  lEnd;
    int32  lNumSegments =    pstData->lNumSegments;
    int32* plAverageData =   pstData->plAverageData;
    int8*  pbyData =         pstData->pbyData;

    while (pstData->bRunNotQuit)
        {
        spcm_vWaitEvent (&pstData->hStart);

        for (j=0; j<lNumSegments; j++)
            {
            lStart =    j * pstData->lSegmentsize + pstData->lStartOffset;
            lEnd =      pstData->lStartOffset + pstData->lAveragesize;

            for (i=lStart; i < lEnd; i++) 
                plAverageData[i] += (int32) pbyData[i];
            }

        spcm_vSignalEvent (&pstData->hEnd);
        }

    return 0;
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
    int8*       pbyData;
    int32*      plAverageData;
    int32*      plStorageData;
    char        szErrorTextBuffer[ERRORTEXTLEN];
    uint32      dwError;
    int32       lStatus;
    int64       llAvailUser, llPCPos;
    int32       lFillsize;
    uint64      qwTotalMem = 0;
    int32       i, j;

    // FIFO mode buffer handling
    int32       lNotifySize =   lSegmentsize * lNumSegments;       
    int64       llBufferSize =   lNotifySize * 16;   // software (DMA) buffer size
    int32       lSegmentCount = 0;                  // number of segments acquired so far
    int32       lAverageCount = 0;                  // current average 

    // settings for the average threads
    SPCM_AVERAGE_DATA   stAverageData[lThreads];

    // -------------------------------------------------
    // open card
    hCard = spcm_hOpen ("/dev/spcm0");
    if (!hCard)
        {
        printf ("no card found...\n");
        return 0;
        }

    // read type, function and sn and check for correct /D card
    spcm_dwGetParam_i32 (hCard, SPC_PCITYP,         &lCardType);
    spcm_dwGetParam_i32 (hCard, SPC_PCISERIALNO,    &lSerialNumber);
    spcm_dwGetParam_i32 (hCard, SPC_FNCTYPE,        &lFncType);

    switch (lFncType)
        {
        case SPCM_TYPE_AI:  
            printf ("Found: %s sn %05d\n", szTypeToName (lCardType), lSerialNumber);
            switch (lCardType & TYP_VERSIONMASK)
                {
                case 0x2230:
                case 0x2233:
                case 0x2234:
                    break;

                default:
                    printf ("The example is set up for the M4i.223x cards and runs with 1 channel 5 GS/s 8 bit. For other digitizers the example settings have to be adopted\n");
                    return 0;
                }
            break;

        default:
            printf ("Card: %s sn %05d not supported by examplee\n", szTypeToName (lCardType), lSerialNumber);            
            return 0;
        }


    // -------------------------------------------------
    // define the data buffers
    pbyData =       (int8*)  pvAllocMemPageAligned ((uint64) llBufferSize);
    plAverageData = (int32*) pvAllocMemPageAligned ((uint64) lSegmentsize * sizeof(int32));
    plStorageData = (int32*) pvAllocMemPageAligned ((uint64) lSegmentsize * sizeof(int32));
    if (!pbyData || !plAverageData || !plStorageData)
        {
        printf ("memory allocation failed\n");
        spcm_vClose (hCard);
        return 0;
        }
    memset ((void*) plAverageData, 0, lSegmentsize * sizeof(int32));

    // -------------------------------------------------
    // set up the average threads
    for (i=0; i<lThreads; i++)
        {
        stAverageData[i].bRunNotQuit =         true;
        stAverageData[i].lNumSegments =        lNumSegments;
        stAverageData[i].lAveragesize =        lSegmentsize / lThreads;
        stAverageData[i].lSegmentsize =        lSegmentsize;
        stAverageData[i].lStartOffset =        i * stAverageData[i].lAveragesize;
        stAverageData[i].pbyData =             pbyData;
        stAverageData[i].plAverageData =       plAverageData;
        spcm_bCreateEvent (&stAverageData[i].hStart);
        spcm_bCreateEvent (&stAverageData[i].hEnd);
        spcm_bCreateThread (pvAverageSegmentPart, &stAverageData[i].hThread, (void*) &stAverageData[i]);
        }


    // -------------------------------------------------
    // do a simple standard setup
    spcm_dwSetParam_i32 (hCard, SPC_CHENABLE,           CHANNEL0);              // just 1 channel enabled

    spcm_dwSetParam_i32 (hCard, SPC_CARDMODE,           SPC_REC_FIFO_MULTI);    // single FIFO mode
    spcm_dwSetParam_i32 (hCard, SPC_LOOPS,              0);                     // endless
    spcm_dwSetParam_i32 (hCard, SPC_SEGMENTSIZE,        lSegmentsize);          // 1k of pretrigger data at start of FIFO mode
    spcm_dwSetParam_i32 (hCard, SPC_POSTTRIGGER,        lSegmentsize - 32);     // 32 samples pretrigger data for each segment
    spcm_dwSetParam_i32 (hCard, SPC_TIMEOUT,            5000);                  // timeout 5 s

    spcm_dwSetParam_i32 (hCard, SPC_TRIG_ORMASK,        SPC_TMASK_EXT0);        // trigger set to external input Ext0
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_ANDMASK,       0);                     // ...
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_EXT0_MODE,     SPC_TM_POS);            // ...
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_EXT0_LEVEL0,   1000);                  // level set to 1000 mV

    spcm_dwSetParam_i32 (hCard, SPC_CLOCKMODE,          SPC_CM_INTPLL);         // clock mode internal PLL
    spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE,         MEGA (2500));           // 5 GS/s
    spcm_dwSetParam_i32 (hCard, SPC_CLOCKOUT,           0);                     // no clock output

    spcm_dwDefTransfer_i64 (hCard, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, lNotifySize, pbyData, 0, llBufferSize);

    // start everything
    dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_DATA_STARTDMA);


    // check for error
    if (dwError != ERR_OK)
        {
        spcm_dwGetErrorInfo_i32 (hCard, NULL, NULL, szErrorTextBuffer);
        printf ("%s\n", szErrorTextBuffer);
        vFreeMemPageAligned (pbyData, (uint64) llBufferSize);
        spcm_vClose (hCard);
        return 0;
        }


    // -------------------------------------------------
    // run the FIFO mode and loop through the data
    while (1)
        {

        // wait for interrupt and check on error
        if ((dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_DATA_WAITDMA)) != ERR_OK)
            {
            if (dwError == ERR_TIMEOUT)
                printf ("... Timeout\n");
            else
                printf ("... Error: %d\n", dwError);
            break;
            }

        // --------------------------------------------------------------------------
        // we have new data available
        else
            {
            spcm_dwGetParam_i32 (hCard, SPC_M2STATUS,             &lStatus);
            spcm_dwGetParam_i64 (hCard, SPC_DATA_AVAIL_USER_LEN,  &llAvailUser);
            spcm_dwGetParam_i64 (hCard, SPC_DATA_AVAIL_USER_POS,  &llPCPos);
            if (lStatus & M2STAT_DATA_OVERRUN)
                {
                printf ("\nOverrun!!!\n");
                break;
                }

            if (llAvailUser >= lNotifySize)
                {
                qwTotalMem += lNotifySize;
                lSegmentCount += lNumSegments;
                lAverageCount += lNumSegments;


                // --------------------------------------------------------------------------
                // summation loop inline
                if (!bThreads)
                    {
                    for (j=0; j<lNumSegments; j++)
                        for (i=0; i<lSegmentsize;i++) 
                            plAverageData[i] += (int32) pbyData[llPCPos + i];
                    }

                // summation loop thread based
                else
                    {
                    for (i=0; i<lThreads; i++)
                        {
                        stAverageData[i].llCurDataPos = llPCPos;
                        spcm_vSignalEvent (&stAverageData[i].hStart);
                        }
                    for (i=0; i<lThreads; i++)
                        spcm_vWaitEvent (&stAverageData[i].hEnd);
                    }

                // --------------------------------------------------------------------------
                // avarage loop reached loop: store data
                if (lAverageCount >= lAverageLoop)
                    {

                    // copy average buffer and clear it for next loop
                    lAverageCount = 0;
                    memcpy ((void*) plStorageData, (void*) plAverageData, lSegmentsize * sizeof(int32));
                    memset ((void*) plAverageData, 0, lSegmentsize * sizeof(int32));

                    // read out buffer fillsize and print it
                    spcm_dwGetParam_i32 (hCard, SPC_FILLSIZEPROMILLE, &lFillsize);
                    printf ("Stat:%08x Segments:%d Fillesize = %4d%%%% Total:%.2fMB\n", lStatus, lSegmentCount, lFillsize, (double) (int64) qwTotalMem / MEGA_B(1));
                    }

                // free the buffer
                spcm_dwSetParam_i32 (hCard, SPC_DATA_AVAIL_CARD_LEN,  lNotifySize);
                }

            // check for esape = abort
#ifdef WIN32
            if (bKeyCheckAsync())
#else
            if (bKbhit ())
#endif
                if (cGetch () == 27)
                    break;
            }
        }

    // send the stop command
    dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_STOP | M2CMD_DATA_STOPDMA);

    // clean up the average threads
    for (i=0; i<lThreads; i++)
        {
        stAverageData[i].bRunNotQuit = false;
        spcm_vSignalEvent (&stAverageData[i].hStart);
        spcm_vWaitEvent (&stAverageData[i].hEnd);
        spcm_vJoinThread (&stAverageData[i].hThread, 100);
        spcm_vCloseThread (&stAverageData[i].hThread);
        }

    // clean up memory
    printf ("Finished...\n");
    vFreeMemPageAligned (pbyData,       (uint64) llBufferSize);
    vFreeMemPageAligned (plAverageData, (uint64) lSegmentsize * sizeof(int32));
    vFreeMemPageAligned (plStorageData, (uint64) lSegmentsize * sizeof(int32));
    spcm_vClose (hCard);

    return EXIT_SUCCESS;
    }

