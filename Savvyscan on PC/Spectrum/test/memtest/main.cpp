/*
**************************************************************************

spcm_memtest.cpp                                         (c) Spectrum GmbH

**************************************************************************

Memory test like in Spectrum Control Center for use in custom applications.

Documentation for the API as well as a detailed description of the hardware
can be found in the manual for each device which can be found on our website:
https://www.spectrum-instrumentation.com/en/downloads

Further information can be found online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/knowledge-base-overview

**************************************************************************
*/

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <iomanip>

#include "../../c_header/dlltyp.h"
#include "../../c_header/regs.h"
#include "../../c_header/spcerr.h"
#include "../../c_header/spcm_drv.h"

#include "../../common/ostools/spcm_oswrap.h"
#include "../../common/ostools/spcm_ostools.h"

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


// ****************************************************************************
// ***** command line Memory test for SPCM cards
// ****************************************************************************

struct ST_SETUP
    {
    uint32 dwNumLoops;
    uint32 dwNumReadsAfterWrite;
    int64  llReducedMemsize;
    int64  llOffset;
    bool   bUseContinuousMemory; // use ContMem if available and large enough
    };

uint32 dwRunMemtest (drv_handle hCard, const struct ST_SETUP* pstSetup);

int main (int argc, char* argv[])
    {
    struct ST_SETUP stSetup = { 1, 1, 0, 0, true };
    char szCard[50];

    // uncomment the second line and replace the IP address to use remote
    // cards like in a digitizerNETBOX
    sprintf (szCard, "/dev/spcm0");
    // sprintf (szCard, "TCPIP::192.168.1.10::inst0::INSTR");


    for (int lArg = 1; lArg < argc; ++lArg)
        {
        if (strcmp (argv[lArg], "--card") == 0)
            {
            // uncomment the second line and replace the IP address to use remote
            // cards like in a digitizerNETBOX
            sprintf (szCard, "/dev/spcm%s", argv[lArg + 1]);
            // sprintf (szCard, "TCPIP::192.168.1.10::inst%s::INSTR", argv[lArg + 1]);
            lArg++;
            }            
        else if (strcmp (argv[lArg], "--loop") == 0)
            {
            int lLoops = atoi (argv[lArg + 1]);
            if (lLoops == -1)
                stSetup.dwNumLoops = 0xFFFFFFFF; // used as "forever"
            else
                stSetup.dwNumLoops = lLoops;
            lArg++;
            }
        else if (strcmp (argv[lArg], "--read") == 0)
            {
            stSetup.dwNumReadsAfterWrite = atoi (argv[lArg + 1]);
            lArg++;
            }
        else if (strcmp (argv[lArg], "--mem") == 0)
            {
            stSetup.llReducedMemsize = atoll (argv[lArg + 1]);
            lArg++;
            }
        else if (strcmp (argv[lArg], "--offset") == 0)
            {
            stSetup.llOffset = atoll (argv[lArg + 1]);
            lArg++;
            }
        else if ((strcmp (argv[lArg], "--help") == 0)
                || (strcmp (argv[lArg], "-h") == 0)
                 || (strcmp (argv[lArg], "/?") == 0))
            {
            std::cout << "Commandline Memtest for SPCM-based cards" << std::endl;
            std::cout << std::endl;
            std::cout << "Syntax: " << argv[0] << " [options]" << std::endl;
            std::cout << "Options:" << std::endl;
            std::cout << "--card <n>:    index of card to be used. Default: 0" << std::endl;
            std::cout << "--loop <n>:    number of times the memtest should run. Default: 1. Use -1 to run test in endless loop." << std::endl;
            std::cout << "--read <n>:    number of times the data should be read from the card. Default: 1" << std::endl;
            std::cout << "--mem <n>:     number of bytes to test. Default: complete onboard memory" << std::endl;
            std::cout << "--offset <n>:  offset in the onboard memory in bytes. Default: 0" << std::endl;
            std::cout << "--help:        print this text" << std::endl;
            std::cout << std::endl;
            return EXIT_SUCCESS;
            }
        else
            {
            std::cerr << "Unknown parameter \"" << argv[lArg] << "\"." << std::endl;
            std::cout << "Use --help to get help." << std::endl;
            return EXIT_FAILURE;
            }
        }

    drv_handle hCard = spcm_hOpen (szCard);
    if (hCard == NULL)
        {
        std::cerr << "Could not open card" << std::endl;
        return EXIT_FAILURE;
        }

    dwRunMemtest (hCard, &stSetup);

    spcm_vClose (hCard);

    return EXIT_SUCCESS;
    }

enum MEMTEST_STATE { MEMTEST_WRITE, MEMTEST_READ };
struct ST_MEMTESTRESULTS
    {
    MEMTEST_STATE eState;
    bool   bOk;
    uint32 dwLoopCount;
    uint32 dwReadCount; // if we read multiple times in one loop
    };

// ----------------------------------------------------------------------------
// ----- Memtest
// ----------------------------------------------------------------------------

// Two reasons to use our own modified implementation of the Microsoft random generator:
// 1. He generates only 15 Bit.
// 2. Slow because of thread safety

uint32 g_dwOwnRand_Hold = 1;

void vOwnRand_Seed (uint32 dwSeed)
    {
    g_dwOwnRand_Hold = dwSeed;
    }

inline uint32 dwOwnRand ()  
    {
    // Microsoft rand.c is a Linear Congruential Generator
    // return(((holdrand = holdrand * 214013L + 2531011L) >> 16) & 0x7fff);

    g_dwOwnRand_Hold = g_dwOwnRand_Hold * 214013L + 2531011L;

    return (g_dwOwnRand_Hold >> 15) & 0xffff;
    }


/*
**************************************************************************
the memtest function
**************************************************************************
*/

#define TEST_PATTERN     ((dwOwnRand () << 16) | dwOwnRand ());

uint32 dwRunMemtest (drv_handle hCard, const struct ST_SETUP* pstSetup)
    {
    uint32 dwReturn =         ERR_OK;
    uint32 dwError =          ERR_OK;
    int64 llBlockSize =       1 * MEGA_B(1);
    uint64 qwContBufLen =     0;
    uint32* pdwBuffer =       NULL;
    uint32 dwPercentageNew =  0;
    uint32 dwPercentageOld =  0;
    uint32 dwBufIdx =         0;
    uint32* pdwBuf =          NULL;

    // ----- no memory test for demo cards -----
    int32 lDemoCard = 0;
    spcm_dwGetParam_i32 (hCard, SPC_MIINST_ISDEMOCARD, &lDemoCard);
    if (lDemoCard != 0)
        {
        std::cerr << "Memtest for Demo cards is not supported." << std::endl;
        return ERR_FNCNOTSUPPORTED;
        }

    // ----- get and print some basic information about the used card -----
    int32 lCardType = 0;
    spcm_dwGetParam_i32 (hCard, SPC_PCITYP, &lCardType);

    int32 lSN = 0;
    spcm_dwGetParam_i32 (hCard, SPC_PCISERIALNR, &lSN);

    int64 llMemsize = 0;
    spcm_dwGetParam_i64 (hCard, SPC_PCIMEMSIZE, &llMemsize);

    std::cout << "Found " << szTypeToName (lCardType) << " sn " << std::setw (5) << std::setfill ('0') << lSN << std::setfill (' ') << " with " << llMemsize / MEGA_B(1) << "MB memory" << std::endl;

    // ----- try to use a continuous buffer for data transfer or allocate a buffer in case there’s none or it's too small. -----
    dwError = spcm_dwGetContBuf_i64 (hCard, SPCM_BUF_DATA, (void**) &pdwBuffer, &qwContBufLen);
    if (!pstSetup->bUseContinuousMemory || (qwContBufLen < (2 * (uint64) llBlockSize)))
        {
        // ContMem not used.
        qwContBufLen = 0;

        pdwBuffer = (uint32*) pvAllocMemPageAligned (2 * llBlockSize);

        // check for mem alloc error
        if (!pdwBuffer)
            return ERR_MEMALLOC;
        }

    // Use a reduced memsize?
    if (pstSetup->llReducedMemsize)
        llMemsize = pstSetup->llReducedMemsize;

    // this struct will hold the results
    ST_MEMTESTRESULTS stMemTestResult;
    memset (&stMemTestResult, 0, sizeof (stMemTestResult));

    // init the random number generator
    srand((uint32) time(NULL));
    vOwnRand_Seed ((uint32) time(NULL));

    // ----- main loop -----
    uint32 dwLoopCnt = 0;
    uint32 dwLoops = pstSetup->dwNumLoops;
    do
        {
        if (pstSetup->dwNumLoops != 0xFFFFFFFF)
            std::cout << "sn " << lSN << ": Loop " << pstSetup->dwNumLoops - dwLoops + 1 << "/" << pstSetup->dwNumLoops << std::endl;
        else
            std::cout << "sn " << lSN << ": Loop " << dwLoopCnt << std::endl;

        // initialization with random start pattern
        uint32 dwStart = rand() * rand() + rand();

        stMemTestResult.eState = MEMTEST_WRITE;

        // basic setup
        dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_RESET);
        if (dwError == ERR_OK)
            dwError = spcm_dwSetParam_i32 (hCard, SPC_MEMTEST, 1);
        if (dwError == ERR_OK)
            dwError = spcm_dwSetParam_i32 (hCard, SPC_TIMEOUT, 1000);

        dwPercentageOld = 0;
        dwBufIdx =        0;
        uint32 dwSRand = rand ();
        srand (dwSRand);

        // ----- write loop -----
        vOwnRand_Seed (dwSRand);

        if (dwError == ERR_OK)
            dwError = spcm_dwDefTransfer_i64 (hCard, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, (uint32) llBlockSize, pdwBuffer, pstSetup->llOffset, 2 * llBlockSize);
        if (dwError == ERR_OK)
            dwError = spcm_dwSetParam_i32 (hCard, SPC_DATA_AVAIL_CARD_LEN, 0);

        uint32 dwCmd = M2CMD_DATA_STARTDMA;
        for (int64 llBlock = 0; (llBlock < llMemsize) && !dwError; llBlock += llBlockSize)
            {
            std::cout << "\rWriting data... " << std::setw (3) << std::fixed << std::setprecision (0) << std::right << (100. * llBlock) / llMemsize << "%";

            // Use the buffers in ping-pong mode.
            if (dwBufIdx == 0)
                {
                pdwBuf = pdwBuffer;
                dwBufIdx = 1; // next loop use other buffer
                }
            else
                {
                pdwBuf = pdwBuffer + llBlockSize / sizeof (uint32);
                dwBufIdx = 0; // next loop use other buffer
                }

            // ----- fill the data block with test pattern -----
            for (int32 i = (uint32) llBlockSize / sizeof (uint32); i; --i)
                {
                *pdwBuf = TEST_PATTERN;
                ++pdwBuf;
                }

            // New data for transfer available.
            if (dwError == ERR_OK)
                dwError = spcm_dwSetParam_i64 (hCard, SPC_DATA_AVAIL_CARD_LEN, llBlockSize);

            // Every 1 % update the progress bar.
            dwPercentageNew = (uint32) (50.0 * llBlock / llMemsize + 0.5);
            if (dwPercentageOld < dwPercentageNew)
                {
                dwPercentageOld = dwPercentageNew;
                }

            // wait until data has been transfered to card and we can fill a new block
            if (dwError == ERR_OK)
                dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, dwCmd);  // first value StartDMA
            dwCmd = M2CMD_DATA_WAITDMA;
            }

        // ----- wait until all data has been transfered to the card -----
        int32 lAvailBytes = 0;
        uint32 dwTime_ms, dwStartTime_ms = SPCM_NAMESPACE::dwGetTickCount ();
        do
            {
            if (dwError == ERR_OK)
                dwError = spcm_dwGetParam_i32 (hCard, SPC_DATA_AVAIL_USER_LEN, &lAvailBytes);

            dwTime_ms = SPCM_NAMESPACE::dwGetTickCount ();

            // overflow, we simply start at 0 again and allow one addtional timeout time for this loop
            if (dwTime_ms < dwStartTime_ms)
                dwStartTime_ms = 0;

            // check for timeout
            if ((dwTime_ms - dwStartTime_ms) > 1000) // random timeout value
                dwError = ERR_TIMEOUT;
            }
        while (!dwError && (lAvailBytes < 2 * llBlockSize));

        if (dwError == ERR_OK)
            dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_DATA_STOPDMA);

        std::cout << std::endl;

        // ----- read + compare loop
        for (uint32 dwReadLoop = 0; dwReadLoop < pstSetup->dwNumReadsAfterWrite && !dwError; ++dwReadLoop)
            {
            srand (dwSRand);
            stMemTestResult.eState = MEMTEST_READ;
            stMemTestResult.dwReadCount = dwReadLoop;

            // clear DMA buffer
            memset (pdwBuffer, 0, (size_t) (2 * llBlockSize));

            vOwnRand_Seed (dwSRand);
            dwPercentageOld = 0;
            dwBufIdx = 0;
            stMemTestResult.bOk = true;

            // Start read transfer and wait for one buffer.
            if (dwError == ERR_OK)
                dwError = spcm_dwDefTransfer_i64 (hCard, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, (uint32) llBlockSize, pdwBuffer, pstSetup->llOffset, 2 * llBlockSize);
            if (dwError == ERR_OK)
                dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);

            for (int64 llBlock = 0; (llBlock < llMemsize) && !dwError; llBlock += llBlockSize)
                {
                std::cout << "\rReading data " << "(" << dwReadLoop + 1 << "/" << pstSetup->dwNumReadsAfterWrite << ") ... " << std::setw (3) << std::fixed << std::setprecision (0) << std::right << (100. * llBlock) / llMemsize << "%";
                uint32 dwErrorCounter = 0;

                // use the other buffer
                if (dwBufIdx == 0)
                    {
                    pdwBuf = pdwBuffer;
                    dwBufIdx = 1; // next loop use other buffer
                    }
                else
                    {
                    pdwBuf = pdwBuffer + llBlockSize / 4;
                    dwBufIdx = 0; // next loop use other buffer
                    }

                uint32 dwNomVal = 0; // holds the expected value
                for (int32 i = (uint32) (llBlockSize / sizeof (uint32)); i; --i, ++pdwBuf)
                    {
                    dwNomVal = TEST_PATTERN;

                    // check if expected value and value that we read back from card match
                    if (*pdwBuf != dwNomVal)
                        {
                        dwErrorCounter++;
                        }
                    }

                // mark memory as available for transfer again
                if (dwError == ERR_OK)
                    dwError = spcm_dwSetParam_i64 (hCard, SPC_DATA_AVAIL_CARD_LEN, llBlockSize);

                // print only a line with the error count if new errors occured
                if (dwErrorCounter)
                    {
                    stMemTestResult.bOk = false;
                    }

                // wait for new data from the board
                if (dwError == ERR_OK)
                    dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_DATA_WAITDMA);
                }   // for: block

            std::cout << " " << (stMemTestResult.bOk? "OK" : "FAIL") << std::endl;

            if (dwError == ERR_OK)
                dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_DATA_STOPDMA);
            }   // for: lReadLoops

        // clear error in driver
        if (dwError != ERR_OK)
            spcm_dwGetErrorInfo_i32 (hCard, NULL, NULL, NULL);

        // Every error should be displayed.
        if (dwError)
            dwReturn = dwError;

        // decrease number of remaining loops if we don't loop forever
        if (pstSetup->dwNumLoops != 0xFFFFFFFF)
            dwLoops--;

        dwLoopCnt++;
        }
    while (dwLoops && (dwReturn == ERR_OK));

    // free DMA memory if we allocated it ourself
    if (qwContBufLen == 0)
        vFreeMemPageAligned (pdwBuffer, 2 * llBlockSize);

    return dwReturn;
    }
