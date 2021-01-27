/*
**************************************************************************

closed_loop_ad_da.cpp                                    (c) Spectrum GmbH

**************************************************************************

This program acquires data using an analog acquisition card, and replays it
immediately using an AWG card.
  
Feel free to use this source for own projects and modify it in any kind

Documentation for the API as well as a detailed description of the hardware
can be found in the manual for each device which can be found on our website:
https://www.spectrum-instrumentation.com/en/downloads

Further information can be found online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/knowledge-base-overview

**************************************************************************
*/

/*
**************************************************************************
This example needs one A/D and one D/A card with the following connection:
A/D X0 (trigger out) to D/A trigger input
A/D clock-out to D/A clock input

Best performance reached:
- compile the release version
- close the developemnt software
- close all other programs running
- directly start the release version
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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <time.h>


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
    uint32      dwError = ERR_OK;
    uint64      qwToTransfer = GIGA_B(1024);

    uint64      qwTotalMem = 0;
    int32       lStatusAD = 0;
    int32       lStatusDA = 0;
    int64       llAvailUserAD = 0;
    int64       llPCPosAD = 0;
    int64       llAvailUserDA = 0;
    int64       llPCPosDA = 0;
    int64       llBytesToProcess = 0;
    int32       lFillsizeAD = 0;
    int32       lFillsizeDA = 0;
	int32		lRefClockOutput;
	time_t		stTime;

    // settings for the example
    int32       lIR_mV = 1000;

	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // these are the main settings that control stability and latency
	int32       lNumCh =			1;				// number of channels active on both cards
	int64       llSampleRate =      KILO(500000);   // sampling rate in kHz
    int32       lNotifySize =       KILO_B(256);	// size of one notify size (data transfererd until interrupt is issued)
    int64       llBufferSizeDA =    4 * lNotifySize;// total size of software buffer
	// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
													
	int64       llBufferSizeAD =    llBufferSizeDA;
	
    // ----- calculate delay -----
    printf ("\nCalculated delay between A/D and D/A:\n");
    printf ("D/A buffer: %0.2lf ms\n\n", 1000.0 * (double)llBufferSizeDA / 2.0 / (double)llSampleRate / (double) lNumCh);

    // ----- open cards -----
    // ----- this program will use the first AD and first DA card it finds -----
    drv_handle hCardAD = NULL_HANDLE;
    drv_handle hCardDA = NULL_HANDLE;
    int lNumCards = 0;
    while ((hCardAD == NULL_HANDLE) || (hCardDA == NULL_HANDLE))
        {
        // open card
        char szDeviceName[12] = { '\0' };
        sprintf (szDeviceName, "/dev/spcm%d", lNumCards);
        drv_handle hCard = spcm_hOpen (szDeviceName);
        if (!hCard)
            {
            printf ("no card found...\n");
            return EXIT_FAILURE;
            }

        // check if the card is AD, DA or something else
        int32 lFncType = 0;
        spcm_dwGetParam_i32 (hCard, SPC_FNCTYPE, &lFncType);
        if (lFncType == SPCM_TYPE_AI)
            hCardAD = hCard;
        else if (lFncType == SPCM_TYPE_AO)
            hCardDA = hCard;
        else
            spcm_vClose (hCard);
        lNumCards++;
        }            
        

    // ----- print some info about cards that are going to be used -----
    int32 lCardType =     0;
    int32 lSerialNumber = 0;
    spcm_dwGetParam_i32 (hCardAD, SPC_PCITYP,         &lCardType);
    spcm_dwGetParam_i32 (hCardAD, SPC_PCISERIALNO,    &lSerialNumber);
    printf ("AD card: %s sn %05d\n", szTypeToName (lCardType), lSerialNumber);

    spcm_dwGetParam_i32 (hCardDA, SPC_PCITYP,         &lCardType);
    spcm_dwGetParam_i32 (hCardDA, SPC_PCISERIALNO,    &lSerialNumber);
    printf ("DA card: %s sn %05d\n", szTypeToName (lCardType), lSerialNumber);

    // calculate channel mask from active channel count
    int32 lChMask = (1 << lNumCh) - 1;

    // ----- setup AD card -----
    spcm_dwSetParam_i32 (hCardAD, SPC_CHENABLE,				lChMask);
    spcm_dwSetParam_i32 (hCardAD, SPC_PRETRIGGER,			32);
    spcm_dwSetParam_i32 (hCardAD, SPC_CARDMODE,				SPC_REC_FIFO_SINGLE);   // single FIFO mode
    spcm_dwSetParam_i32 (hCardAD, SPC_TIMEOUT,				5000);                  // timeout 5 s
    spcm_dwSetParam_i32 (hCardAD, SPC_TRIG_ORMASK,			SPC_TMASK_EXT0);        // waiting for a trigger that isn't connected to allow force trigger
    spcm_dwSetParam_i32 (hCardAD, SPC_TRIG_ANDMASK,			0);                     // ...
    spcm_dwSetParam_i32 (hCardAD, SPCM_X0_MODE,				SPCM_XMODE_TRIGOUT);    // trigger output to D/A card
    spcm_dwSetParam_i32 (hCardAD, SPC_CLOCKMODE,			SPC_CM_INTPLL);         // clock mode internal PLL
    spcm_dwSetParam_i32 (hCardAD, SPC_CLOCKOUT,				1);                     // clock output
    spcm_dwSetParam_i32 (hCardAD, SPC_SPECIALCLOCK,			1);                     // use special clock to be able to use odd sampling rates on M4i.44xx
    spcm_dwSetParam_i64 (hCardAD, SPC_SAMPLERATE,			llSampleRate);
    for (int i = 0; i < lNumCh; ++i)
        {
        spcm_dwSetParam_i64 (hCardAD, SPC_PATH0 + i*100,  0);						// path 0
        spcm_dwSetParam_i64 (hCardAD, SPC_AMP0  + i*100,  lIR_mV);					// input range as defined above
        spcm_dwSetParam_i64 (hCardAD, SPC_ACDC0 + i*100,  COUPLING_DC);				// use DC coupling
        }
    spcm_dwSetParam_i32 (hCardAD, SPC_M2CMD,				M2CMD_CARD_WRITESETUP); // write setup to avoid later settling wait time and to activate clock output
	spcm_dwGetParam_i32 (hCardAD, SPC_CLOCKOUTFREQUENCY,	&lRefClockOutput);		// read the refrence clock output frequency to write this to the D/A card

    // ----- setup DA card -----
    spcm_dwSetParam_i32 (hCardDA, SPC_CHENABLE,				lChMask);
    spcm_dwSetParam_i32 (hCardDA, SPC_CARDMODE,				SPC_REP_FIFO_SINGLE);   // single FIFO mode
    spcm_dwSetParam_i32 (hCardDA, SPC_TIMEOUT,				5000);                  // timeout 5 s
    spcm_dwSetParam_i32 (hCardDA, SPC_TRIG_ORMASK,			SPC_TMASK_EXT0);        // wait fro trigger from A/D card
    spcm_dwSetParam_i32 (hCardDA, SPC_TRIG_ANDMASK,			0);                     // ...
    spcm_dwSetParam_i32 (hCardDA, SPC_TRIG_EXT0_MODE,		SPC_TM_POS);            // positive edge
    spcm_dwSetParam_i32 (hCardDA, SPC_TRIG_TERM,			0);                     // high impedance input
    spcm_dwSetParam_i32 (hCardDA, SPC_TRIG_EXT0_LEVEL0,		1500);                  // 1.5 V trigger level
    spcm_dwSetParam_i32 (hCardDA, SPC_CLOCKMODE,			SPC_CM_EXTREFCLOCK);    // clock mode external refrence clock (from A/D card)
	spcm_dwSetParam_i32 (hCardDA, SPC_REFERENCECLOCK,		lRefClockOutput);		// reference clock from A/D card
    spcm_dwSetParam_i32 (hCardDA, SPC_CLOCKOUT,				0);                     // no clock output
    spcm_dwSetParam_i64 (hCardDA, SPC_SAMPLERATE,			llSampleRate);
    for (int i = 0; i < lNumCh; ++i)
        {
        spcm_dwSetParam_i64 (hCardDA, SPC_AMP0       + i*100,  lIR_mV);				// use same voltage level as the AD card (also see below for DATACONVERSION)
        spcm_dwSetParam_i64 (hCardDA, SPC_ENABLEOUT0 + i*100,  1);					// enable the output
        }
    spcm_dwSetParam_i32 (hCardDA, SPC_M2CMD,				M2CMD_CARD_WRITESETUP); // write setup to avoid later settling wait time

    // ----- if resolution of AD card and DA card do not match we use the data conversion feature of the DA card to get similar output levels -----
    int32 lResolutionAD = 0;
    int32 lResolutionDA = 0;
    spcm_dwGetParam_i32 (hCardAD, SPC_MIINST_BITSPERSAMPLE, &lResolutionAD);
    spcm_dwGetParam_i32 (hCardDA, SPC_MIINST_BITSPERSAMPLE, &lResolutionDA);
    switch (lResolutionDA)
        {
        case 16:
            {
            if (lResolutionAD == 14)
                spcm_dwSetParam_i32 (hCardDA, SPC_DATACONVERSION, SPCM_DC_14BIT_TO_16BIT);
            else if (lResolutionAD == 12)
                spcm_dwSetParam_i32 (hCardDA, SPC_DATACONVERSION, SPCM_DC_12BIT_TO_16BIT);
            break;
            }
        case 14:
            {
            if (lResolutionAD == 16)
                spcm_dwSetParam_i32 (hCardDA, SPC_DATACONVERSION, SPCM_DC_16BIT_TO_14BIT);
            else if (lResolutionAD == 12)
                spcm_dwSetParam_i32 (hCardDA, SPC_DATACONVERSION, SPCM_DC_12BIT_TO_14BIT);
            break;
            }
        }

    // ----- print used sample rates of both cards -----
    int64 llUsedSampleRate = 0;
    spcm_dwGetParam_i64 (hCardAD, SPC_SAMPLERATE,     &llUsedSampleRate);
    printf ("Used Samplerate AD: %lld\n", llUsedSampleRate);
    spcm_dwGetParam_i64 (hCardDA, SPC_SAMPLERATE,     &llUsedSampleRate);
    printf ("Used Samplerate DA: %lld\n", llUsedSampleRate);

    // ----- define the data buffers -----
    int16* pnDataAD = (int16*) pvAllocMemPageAligned ((uint64) llBufferSizeAD);
    if (!pnDataAD)
        {
        printf ("memory allocation failed\n");
        spcm_vClose (hCardAD);
        spcm_vClose (hCardDA);
        return EXIT_FAILURE;
        }

    int16* pnDataDA = (int16*) pvAllocMemPageAligned ((uint64) llBufferSizeDA);
    if (!pnDataDA)
        {
        printf ("memory allocation failed\n");

        vFreeMemPageAligned (pnDataAD, (uint64) llBufferSizeAD);
        spcm_vClose (hCardAD);
        spcm_vClose (hCardDA);
        return EXIT_FAILURE;
        }

    // define the buffers for transfer
    dwError = spcm_dwDefTransfer_i64 (hCardAD, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, lNotifySize, pnDataAD, 0, llBufferSizeAD);
    dwError = spcm_dwDefTransfer_i64 (hCardDA, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, lNotifySize, pnDataDA, 0, llBufferSizeDA);

    printf ("\nClosed Loop starts now. Press CTRL-C to abort\n... no output during loop for best performance\n");
	stTime = time(NULL);
	printf("\nStart Time: %s\n", ctime(&stTime));

    // -------------------------------------------------------
    // -------------------------------------------------------
    // ----- run the FIFO mode and loop through the data -----

    // ----- transfer empty buffer data to the D/A card -----
    // this defines the delay between input and output of the loop
    if (!dwError)
        {
        dwError = spcm_dwGetParam_i64(hCardDA, SPC_DATA_AVAIL_USER_LEN, &llAvailUserDA);
        dwError = spcm_dwGetParam_i64(hCardDA, SPC_DATA_AVAIL_USER_POS, &llPCPosDA);
        dwError = spcm_dwSetParam_i64(hCardDA, SPC_DATA_AVAIL_CARD_LEN, llBufferSizeDA);
        dwError = spcm_dwSetParam_i32(hCardDA, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);
        }

    // ----- start D/A and wait until it is armed -----
    if (!dwError)
        {
        dwError = spcm_dwSetParam_i32(hCardDA, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER);
        do
            {
            dwError = spcm_dwGetParam_i32(hCardDA, SPC_M2STATUS, &lStatusDA);
            } while ((lStatusDA & M2STAT_CARD_PRETRIGGER) == 0);
        }

    // ----- start AD card -----
    // this will automatically trigger the D/A card and start the output of the empty buffer
    if (!dwError)
        {
        dwError = spcm_dwSetParam_i32(hCardAD, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER);
        dwError = spcm_dwSetParam_i32(hCardAD, SPC_M2CMD, M2CMD_DATA_STARTDMA);
        dwError = spcm_dwSetParam_i32(hCardAD, SPC_M2CMD, M2CMD_CARD_FORCETRIGGER); // starts the acquisition and via trigger output the generation also
        spcm_dwSetParam_i32(hCardAD, SPCM_XX_ASYNCIO, 7);
        }

    // ----- LOOP LOOP LOOP -------------------------------------------------------------------------------
    // in the loop each buffer of data received from the A/D card is copied and transferred to the D/A card
    while (!dwError && (qwTotalMem < qwToTransfer))
        {

        // wait for next A/D buffer to be ready
        dwError = spcm_dwSetParam_i32 (hCardAD, SPC_M2CMD, M2CMD_DATA_WAITDMA);

        // process next block of data
        if (!dwError)
            {
            dwError = spcm_dwGetParam_i32(hCardAD, SPC_M2STATUS, &lStatusAD);
            spcm_dwGetParam_i64(hCardAD, SPC_DATA_AVAIL_USER_LEN, &llAvailUserAD);
            spcm_dwGetParam_i64(hCardAD, SPC_DATA_AVAIL_USER_POS, &llPCPosAD);
            if (lStatusAD & M2STAT_CARD_READY)
                {
                printf("\n\nA/D card stopped acquisition\n");
                dwError = ERR_ABORT;
                }
            }

        if (!dwError)
            {
            dwError = spcm_dwGetParam_i32(hCardDA, SPC_M2STATUS, &lStatusDA);
            if (lStatusDA & M2STAT_CARD_READY)
                {
                printf("\n\nD/A card stopped replay\n");
                dwError = ERR_ABORT;
                }
            }

        // ----- process data if at least one notify size is available (this should always be the case after WAITDMA) -----
        if (!dwError && (llAvailUserAD >= lNotifySize))
            {
            spcm_dwGetParam_i64 (hCardDA, SPC_DATA_AVAIL_USER_LEN,  &llAvailUserDA);
            spcm_dwGetParam_i64 (hCardDA, SPC_DATA_AVAIL_USER_POS,  &llPCPosDA);

            // we will handle each block serparately to get best performance
            llBytesToProcess = (int64) lNotifySize;

            // this is the point to do something with the data
            // we will simply copy it unmodified to the buffer of the DA card
            memcpy (pnDataDA + (llPCPosDA / sizeof (int16)), pnDataAD + (llPCPosAD / sizeof (int16)), (size_t) llBytesToProcess);

            // mark data bytes as processed (=free) for AD card
            spcm_dwSetParam_i64 (hCardAD, SPC_DATA_AVAIL_CARD_LEN,  llBytesToProcess);

            // mark data bytes as available for DA card and wait for end of data transfer
            spcm_dwSetParam_i64 (hCardDA, SPC_DATA_AVAIL_CARD_LEN,  llBytesToProcess);
            dwError = spcm_dwSetParam_i32(hCardDA, SPC_M2CMD, M2CMD_DATA_WAITDMA);

            qwTotalMem += llBytesToProcess;

            // ----- print some status info -----
            // to get best performance no printf function should be used in the loop
            spcm_dwGetParam_i32 (hCardAD, SPC_FILLSIZEPROMILLE, &lFillsizeAD);
            spcm_dwGetParam_i32 (hCardDA, SPC_FILLSIZEPROMILLE, &lFillsizeDA);

            // printf decreases performance a lot!
            //printf ("\rStatAD:%08x StatDA:%08x FillAD: %3d FillDA: %3d Total:%.2fMB", lStatusAD, lStatusDA, lFillsizeAD, lFillsizeDA,  (double) (int64) qwTotalMem / MEGA_B(1));
            }
        }
    // -------------------------------------------------------
    // -------------------------------------------------------
    // -------------------------------------------------------

    printf("\n\n");

    // check for timeout
    if (dwError == ERR_TIMEOUT)
        printf("Timeout occurred, reduce speed or enlarge buffers\n");

    // error checking
    uint32 dwADError = 0;
    uint32 dwADErrorReg = 0;
    int32  lADErrorValue = 0;
    char szADErrorText[ERRORTEXTLEN];

    uint32 dwDAError = 0;
    uint32 dwDAErrorReg = 0;
    int32  lDAErrorValue = 0;
    char szDAErrorText[ERRORTEXTLEN];

    dwADError = spcm_dwGetErrorInfo_i32 (hCardAD, &dwADErrorReg, &lADErrorValue, szADErrorText);
    if (dwADError)
        printf("AD Error: %d\n%s\n", dwADError, szADErrorText);
    dwDAError = spcm_dwGetErrorInfo_i32 (hCardDA, &dwDAErrorReg, &lDAErrorValue, szDAErrorText);
    if (dwDAError)
        printf("DA Error: %d\n%s\n", dwDAError, szDAErrorText);

    // send the stop command
    dwError = spcm_dwSetParam_i32 (hCardDA, SPC_M2CMD, M2CMD_CARD_STOP | M2CMD_DATA_STOPDMA);
    dwError = spcm_dwSetParam_i32 (hCardAD, SPC_M2CMD, M2CMD_CARD_STOP | M2CMD_DATA_STOPDMA);

    // clean up
    printf ("Finished...\n");
	stTime = time(NULL);
	printf("End Time: %s\n", ctime(&stTime));

    vFreeMemPageAligned (pnDataAD, (uint64) llBufferSizeAD);
    vFreeMemPageAligned (pnDataDA, (uint64) llBufferSizeDA);	stTime = time(NULL);

    spcm_vClose (hCardAD);
    spcm_vClose (hCardDA);

    cGetch();

    if (dwError)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
    }
