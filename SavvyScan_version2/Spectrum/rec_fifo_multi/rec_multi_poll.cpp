/*
**************************************************************************

rec_multi_poll.cpp                                       (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based acquisition cards with the option 
Multiple Recording installed. 

Information about the different products and their drivers can be found
online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/platform-driver-and-series-differences

Shows Multiple Recording data acquisition using fifo mode. If timestamp 
is installed the corresponding timestamp values are also read
out and displayed.

If Timestamp and BaseXIO are installed the BaseXIO lines are set to the
timestamp acquisition mode and are sampled on every trigger event. The
samples BaseXIO lines are also shown.
  
Feel free to use this source for own projects and modify it in any kind.

This example shows the polling of the status, no wait function is used

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

// ----- include of common example librarys -----
#include "../common/spcm_lib_card.h"
#include "../common/spcm_lib_data.h"

// ----- operating system dependent functions for thead, event and mutex handling -----
#include "../common/ostools/spcm_oswrap.h"
#include "../common/ostools/spcm_ostools.h"


// ----- standard c include files -----
#include <stdio.h>
#include <stdlib.h>



/*
**************************************************************************
vDoCardSetup
**************************************************************************
*/

void vDoCardSetup (ST_SPCM_CARDINFO *pstCard, int32 lSegmentsize, int32 lPosttrigger, int32 lLoops)
    {
    int   i;

    // Multiple Recording setup (1 channel to keep example simple)
    bSpcMSetupModeRecFIFOMulti (pstCard, CHANNEL0, lSegmentsize, lPosttrigger, lLoops);

    // we try to set the samplerate to 1 MHz (M2i) or 20 MHz (M3i) on internal PLL, no clock output
    // increase this to test the read-out-after-overrun
	if (pstCard->bM2i)
        bSpcMSetupClockPLL (pstCard, MEGA(1), false);
	else if (pstCard->bM2p)
        bSpcMSetupClockPLL (pstCard, MEGA(10), false);
	else if (pstCard->bM3i || pstCard->bM4i)
        bSpcMSetupClockPLL (pstCard, MEGA(20), false);
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstCard->llSetSamplerate / 1000000);

    // we set trigger to external positive edge, please connect the trigger line!
    printf ("\n");
    printf ("*************************************************************************\n");
    printf ("* Using external trigger - please connect a signal to the trigger input *\n");
    printf ("* Example is best working with a 10 kHz Trigger signal                  *\n");
    printf ("*************************************************************************\n\n");
    bSpcMSetupTrigExternal (pstCard, SPC_TM_POS, false, 0);

    // type dependent card setup
    switch (pstCard->eCardFunction)
        {

        // analog acquisition card setup
        case AnalogIn:

            // program all input channels to +/-1 V and 50 ohm termination (if it's available)
            if (pstCard->bM2i || pstCard->bM2p)
                {
                for (i=0; i < pstCard->lMaxChannels; i++)
                    bSpcMSetupInputChannel (pstCard, i, 1000, true);
                }
            else
                {
                bool bTerm           = true;
                bool bACCoupling     = false;
                bool bBandwidthLimit = false;
                for (i=0; i < pstCard->lMaxChannels; i++)
                    bSpcMSetupPathInputCh (pstCard, i, 0, 1000, 0, bTerm, bACCoupling, bBandwidthLimit);
                }
            break;
        }


    // if timestamp and basexio are installed we acquire the asynchronous inputs with the timestamps
    if ((pstCard->lFeatureMap & SPCM_FEAT_TIMESTAMP) && (pstCard->lFeatureMap & SPCM_FEAT_BASEXIO))
        bSpcMSetupTimestamp (pstCard, SPC_TSMODE_STANDARD | SPC_TSCNT_INTERNAL | SPC_TSXIOACQ_ENABLE, 0);

    // set up the timestamp mode to standard if timestamp is installed
    else if ((pstCard->lFeatureMap & SPCM_FEAT_TIMESTAMP) != 0)
        bSpcMSetupTimestamp (pstCard, SPC_TSMODE_STANDARD | SPC_TSCNT_INTERNAL, 0);

    }



/*
**************************************************************************
nShowData
**************************************************************************
*/

int16 nShowData (ST_SPCM_CARDINFO *pstCard, int32 lDataAvailBytes, void* pvDataCurrentBuf)
	{

	// add your data processing here
    printf ("Data Buffer Ready: %d Bytes\n", lDataAvailBytes);

	return ERR_OK;
	}



/*
**************************************************************************
nShowTimestamps 
**************************************************************************
*/

int16 nShowTimestamps (ST_SPCM_CARDINFO *pstCard, int32 lTSAvailBytes, uint64* pqwTSCurrentBuf)
	{

	// add your timestamp processing here
    printf ("Timestamp data ready: %d Bytes\n", lTSAvailBytes);

	return ERR_OK;
	}



/*
**************************************************************************
main 
**************************************************************************
*/

int main ()
    {
    char                szBuffer[1024];     // a character buffer for any messages
    ST_SPCM_CARDINFO    stCard;             // info structure of my card
	int32				lDataBufLen, lDataNotify, lTSBufLen, lTSNotify; 
	int32				lAvailPos, lDataAvailBytes = 0, lTSAvailBytes = 0;
    int32               lSegmentsize, lPosttrigger, lLoops;
    int32               lOversampling = 1;
    bool                bTimestampInstalled, bBaseXIOInstalled;
    uint64*             pqwTSBuffer = NULL;
	uint64*				pqwTSCurrentBuf;
    void*               pvDataBuffer;
	void*				pvDataCurrentBuf;
	int64				llDataTransferred = 0, llTSTransferred = 0;	
    uint32              dwErr;
	bool				bDataReady, bTimestampsReady,  bReady = false;



    // ------------------------------------------------------------------------
    // init card number 0 (the first card in the system), get some information and print it
    // uncomment the second line and replace the IP address to use remote
    // cards like in a digitizerNETBOX
    if (bSpcMInitCardByIdx (&stCard, 0))
    //if (bSpcMInitCardByIdx (&stCard, "192.168.1.10", 0))
        {
        printf (pszSpcMPrintDocumentationLink (&stCard, szBuffer, sizeof (szBuffer)));
        printf (pszSpcMPrintCardInfo (&stCard, szBuffer, sizeof (szBuffer)));
        }
    else
        return nSpcMErrorMessageStdOut (&stCard, "Error: Could not open card\n", true);


    // check whether we support this card type in the example
    if (stCard.eCardFunction != AnalogIn)
        return nSpcMErrorMessageStdOut (&stCard, "Error: Card function not supported by this example\n", false);
    if ((stCard.lFeatureMap & SPCM_FEAT_MULTI) == 0)
        return nSpcMErrorMessageStdOut (&stCard, "Error: Multiple Recording Option not installed. Examples was done especially for this option!\n", false);

    // if timestamp is installed we set a flag to support this mode in the example also
    bTimestampInstalled = ((stCard.lFeatureMap & SPCM_FEAT_TIMESTAMP) != 0);
    bBaseXIOInstalled = ((stCard.lFeatureMap & SPCM_FEAT_BASEXIO) != 0);



    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
	lDataBufLen  = MEGA_B(1);
	lDataNotify	 = KILO_B(4);
	lTSBufLen    = KILO_B(128);
	lTSNotify    = KILO_B(4);
    lSegmentsize = KILO_B(4);
    lPosttrigger = lSegmentsize / 2;

	lLoops = 0; // recording infinite number of segments

    if (!stCard.bSetError)
        vDoCardSetup (&stCard, lSegmentsize, lPosttrigger, lLoops);



    // ------------------------------------------------------------------------
    // allocate memory buffer
    if (!stCard.bSetError)
        {
        pvDataBuffer = (void*) new uint8[(int) lDataBufLen];
        if (bTimestampInstalled)
            pqwTSBuffer = new uint64[(int) lTSBufLen];
        if (!pvDataBuffer || (bTimestampInstalled &&!pqwTSBuffer))
            return nSpcMErrorMessageStdOut (&stCard, "Memory allocation error\n", false);
        }



    // ------------------------------------------------------------------------
    // start the acquisition

    // if using timestamps we need to start the transfer before the card start to avoid an overrun of the timestamp memory
    if (bTimestampInstalled)
        {
        printf ("Starting the timestamp DMA transfer\n");
        spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_TIMESTAMP, SPCM_DIR_CARDTOPC, lTSNotify, (void*) pqwTSBuffer, 0, lTSBufLen);
        spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_EXTRA_STARTDMA);
        }

    // We'll define the buffer for data to start everything together
    spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, lDataNotify, pvDataBuffer, 0, lDataBufLen);

    printf ("Starting the card, the data transfer and poll\n");
    spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_DATA_STARTDMA);


    while (!bReady)
        {
        if (!stCard.bSetError)
            {
			printf ("\n");

            // this is our status polling loop
            int32 lOldStatus = 0, lStatus = 0;
			int32 lOldTrigCount = 0, lTrigCount = 0;

            bTimestampsReady = false;
            bDataReady =       false;
            do
                {
                dwErr = spcm_dwGetParam_i32 (stCard.hDrv, SPC_M2STATUS, &lStatus);

                // check for changed status and set flags
                if (lOldStatus != lStatus)
                    {
                    printf ("Status: ");
                    if (lStatus & M2STAT_CARD_PRETRIGGER)
                        printf ("Armed ");

                    if (lStatus & M2STAT_CARD_TRIGGER)
                        printf ("1stTrigger ");

                    if (lStatus & M2STAT_CARD_READY)
                        {
                        printf ("CardReady ");
                        bReady = true;
                        }

                    if (lStatus & M2STAT_DATA_BLOCKREADY)
                        {
                        printf ("DataReady ");
                        bDataReady = true;
                        }

                    if (lStatus & M2STAT_EXTRA_BLOCKREADY)
                        {
                        printf ("TimestampReady ");
                        bTimestampsReady = true;
                        }
                    printf ("\n");
                    }

                lOldStatus = lStatus;

				// check trigger count, available with firmware V10 (M2i) or V20 (M3i)
/*
                if (!dwErr)
                    {
                    dwErr = spcm_dwGetParam_i32 (stCard.hDrv, SPC_TRIGGERCOUNTER, &lTrigCount);
                    if (lTrigCount != lOldTrigCount)
                        printf ("TrigCount = %d\n", lTrigCount);
                    lOldTrigCount = lTrigCount;
                    }
*/
                }
            while (!dwErr && !(bReady || bTimestampsReady || bDataReady));

            // check for error code
            if (spcm_dwGetErrorInfo_i32 (stCard.hDrv, NULL, NULL, szBuffer))
                {
                delete [] ((uint8*) pvDataBuffer);
				if (bTimestampInstalled)
					delete [] (pqwTSBuffer);
                return nSpcMErrorMessageStdOut (&stCard, szBuffer, false);
                }

            // process data and set dat aavailable again
			if (bDataReady)
				{
				spcm_dwGetParam_i32 (stCard.hDrv, SPC_DATA_AVAIL_USER_LEN, &lDataAvailBytes);
				spcm_dwGetParam_i32 (stCard.hDrv, SPC_DATA_AVAIL_USER_POS, &lAvailPos);
				if ((lAvailPos + lDataAvailBytes) >= lDataBufLen)
					lDataAvailBytes = (uint32) (lDataBufLen - lAvailPos);
				pvDataCurrentBuf = (void*) (((char*) pvDataBuffer) + lAvailPos);
				llDataTransferred += lDataAvailBytes;

                nShowData (&stCard, lDataAvailBytes, pvDataCurrentBuf);

				dwErr = spcm_dwSetParam_i32 (stCard.hDrv, SPC_DATA_AVAIL_CARD_LEN, lDataAvailBytes);
				}

            // process timestamps and set data aavailable again
			if (bTimestampsReady && bTimestampInstalled)
				{
				spcm_dwGetParam_i32 (stCard.hDrv, SPC_TS_AVAIL_USER_LEN, &lTSAvailBytes);
				spcm_dwGetParam_i32 (stCard.hDrv, SPC_TS_AVAIL_USER_POS, &lAvailPos);

                // align available bytes to full timestamps
                if (stCard.bM4i || stCard.bM2p)
                    lTSAvailBytes &= 0xFFFFFFF0; // 16 byte per stamp
                else
                    lTSAvailBytes &= 0xFFFFFFF8; //  8 byte per stamp

				if ((lAvailPos + lTSAvailBytes) >= lTSBufLen)
					lTSAvailBytes = (uint32) (lTSBufLen - lAvailPos);
				pqwTSCurrentBuf = (uint64*) (((char*) pqwTSBuffer) + lAvailPos);
				llTSTransferred +=   lTSAvailBytes;

                nShowTimestamps (&stCard, lTSAvailBytes, pqwTSCurrentBuf);

				dwErr = spcm_dwSetParam_i32 (stCard.hDrv, SPC_TS_AVAIL_CARD_LEN, lTSAvailBytes);
    			}
            }

		if (bKbhit ())
			{
			spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_STOP);
			cGetch();
			}
        }    


    // ------------------------------------------------------------------------
    // print error information if an error occured
    if (stCard.bSetError)
        return nSpcMErrorMessageStdOut (&stCard, "An error occured while programming the card:\n", true);


    // clean up and close the driver
    vSpcMCloseCard (&stCard);
    delete [] ((uint8*) pvDataBuffer);
    if (bTimestampInstalled)
        delete [] (pqwTSBuffer);

    return EXIT_SUCCESS;
    }

