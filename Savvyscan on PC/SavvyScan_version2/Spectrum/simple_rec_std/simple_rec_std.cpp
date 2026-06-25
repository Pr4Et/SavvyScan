/*
**************************************************************************

simple_rec_std.cpp                                      (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based analog acquisition cards. 

Information about the different products and their drivers can be found
online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/platform-driver-and-series-difference

Shows a simple standard mode example using only the few necessary commands
  
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
    int16*      pnData = NULL;
	int8*       pbyData = NULL;
    char        szErrorTextBuffer[ERRORTEXTLEN];
    uint32      dwError;
    
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

	int64 llMemsize = 16384; // set memsize to 16 kS

	// do a simple standard setup
	spcm_dwSetParam_i32 (hCard, SPC_CARDMODE,        SPC_REC_STD_SINGLE);    // single standard mode
    spcm_dwSetParam_i32 (hCard, SPC_CHENABLE,        CHANNEL0);              // just 1 channel enabled
	spcm_dwSetParam_i64 (hCard, SPC_MEMSIZE,         llMemsize);             // acquire 16 kS in total
	spcm_dwSetParam_i64 (hCard, SPC_POSTTRIGGER,     llMemsize / 2);         // half of the total number of samples after trigger event
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_ORMASK,     SPC_TMASK_SOFTWARE);    // trigger set to software
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_ANDMASK,    0);                     // ...
	spcm_dwSetParam_i32 (hCard, SPC_CLOCKMODE,       SPC_CM_INTPLL);         // clock mode internal PLL
	spcm_dwSetParam_i32 (hCard, SPC_CLOCKOUT,        0);                     // no clock output
    spcm_dwSetParam_i32 (hCard, SPC_TIMEOUT,         5000);                  // timeout 5 s
    
    // set card to maximum sampling rate
	int64 llMaxSamplingrate = 0;
	spcm_dwGetParam_i64 (hCard, SPC_PCISAMPLERATE, &llMaxSamplingrate);
	spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE, llMaxSamplingrate);

	
	// start card and wait for card ready
    dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY);

	// check for error
    if (dwError != ERR_OK)
        {
        spcm_dwGetErrorInfo_i32 (hCard, NULL, NULL, szErrorTextBuffer);
        printf ("%s\n", szErrorTextBuffer);
        spcm_vClose (hCard);
        return 0;
        }


	// setup transfer buffer and start DMA to transfer data from card to pc memory

	// read bytes per sample value (8 bit cards = 1 bytes, 12, 14, 16 bit cards = 2 bytes)
	int32 lBytesPerSample = 0;
	spcm_dwGetParam_i32 (hCard, SPC_MIINST_BYTESPERSAMPLE, &lBytesPerSample);
	
	// define data buffer
	int64 llBufferSize = llMemsize * lBytesPerSample;

	switch (lBytesPerSample)
		{
		case 1:
			pbyData = (int8*) pvAllocMemPageAligned (llBufferSize);
			if (pbyData)
				spcm_dwDefTransfer_i64 (hCard, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, 0, pbyData, 0, llBufferSize);
			break;

		case 2:
			pnData = (int16*) pvAllocMemPageAligned (llBufferSize);
			if (pnData)
				spcm_dwDefTransfer_i64 (hCard, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, 0, pnData, 0, llBufferSize);
			break;
		}

    if (!pbyData && !pnData)
        {
        printf ("memory allocation failed\n");
        spcm_vClose (hCard);
        return 0;
        }

	// start DMA and wait for DMA transfer ready state
    dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);
	

	// get minimum and maximum values from data
	int32 lMin = 32767;
	int32 lMax = -32767;

	for (int64 llDataIdx = 0; llDataIdx < llMemsize; llDataIdx++)
		{
		switch (lBytesPerSample)
			{
			case 1:
				if (pbyData[llDataIdx] < lMin) lMin = pbyData[llDataIdx];
				if (pbyData[llDataIdx] > lMax) lMax = pbyData[llDataIdx];
				break;
			case 2:
				if (pnData[llDataIdx] < lMin) lMin = pnData[llDataIdx];
				if (pnData[llDataIdx] > lMax) lMax = pnData[llDataIdx];
				break;
			}
		}

	printf ("\nMinimum: %d\n", lMin);
	printf ("Maximum: %d\n\n", lMax);

	// close card
	spcm_vClose (hCard);

	switch (lBytesPerSample)
		{
		case 1: vFreeMemPageAligned (pbyData, llBufferSize); break;
        case 2: vFreeMemPageAligned (pnData,  llBufferSize); break;
        }

	return EXIT_SUCCESS;
    }

