/*
**************************************************************************

rec_fifo_single_poll.cpp                                 (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based analog and digital acquisition cards. 

Information about the different products and their drivers can be found
online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/platform-driver-and-series-differences

Shows FIFO acquisitition.

This example doesn't use interrupt mode but shows status polling. Be 
aware that this of course uses complete CPU time just for polling data.

This example also shows how to read out the remaining data from the on-board 
memory in case that an overung has occurred.

Feel free to use this source for own projects and modify it in any kind

Documentation for the API as well as a detailed description of the hardware
can be found in the manual for each device which can be found on our website:
https://www.spectrum-instrumentation.com/en/downloads

Further information can be found online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/knowledge-base-overview

**************************************************************************
*/



// ----- standard c include files -----
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>



// ----- include of common example librarys -----
#include "../common/spcm_lib_card.h"
#include "../common/spcm_lib_data.h"
#include "../common/ostools/spcm_oswrap.h"
#include "../common/ostools/spcm_ostools.h"



/*
**************************************************************************
bDoCardSetup: setup of acquisition card
**************************************************************************
*/

bool bDoCardSetup (ST_SPCM_CARDINFO *pstCard)
    {
    int     i;

    // we try to set the samplerate to 1 MHz (M2i) or 20 MHz (M3i, M4i) on internal PLL, no clock output
    // increase this to test the read-out-after-overrun
	if (pstCard->bM2i)
        bSpcMSetupClockPLL (pstCard, MEGA(1), false);
	else if (pstCard->bM2p)
        bSpcMSetupClockPLL (pstCard, MEGA(10), false);
	else if (pstCard->bM3i || pstCard->bM4i)
        bSpcMSetupClockPLL (pstCard, MEGA(20), false);

    printf ("Sampling rate set to %.1lf MHz\n", (double) pstCard->llSetSamplerate / 1000000);

    // we set software trigger, no trigger output
    bSpcMSetupTrigSoftware (pstCard, false);

    // type dependent card setup
    switch (pstCard->eCardFunction)
        {
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

            // FIFO mode setup, we run continuously 
            // all channels are activated
            bSpcMSetupModeRecFIFOSingle (pstCard, (1 << pstCard->lMaxChannels) - 1, 32);
            break;

        case DigitalIn:
        case DigitalIO:

            // FIFO mode setup, we run continuously, 16 channels activated
            bSpcMSetupModeRecFIFOSingle (pstCard, 0xffff, 4);
            break;
        }

    return pstCard->bSetError;
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
    void*               pvBuffer = NULL;
    uint32              dwErr;

    // setup for the FIFO mode
    int64        llSWBufSize =      KILO_B(512);
    int64        llNotifySize =     KILO_B(64);



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
    if ((stCard.eCardFunction != AnalogIn) && (stCard.eCardFunction != DigitalIn) && (stCard.eCardFunction != DigitalIO))
        return nSpcMErrorMessageStdOut (&stCard, "Error: Card function not supported by this example\n", false);



    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
    if (!stCard.bSetError)
        bDoCardSetup (&stCard);



    // ------------------------------------------------------------------------
    // allocate and setup the fifo buffer
    pvBuffer = pvAllocMemPageAligned ((uint32) llSWBufSize);
    if (!pvBuffer)
        return nSpcMErrorMessageStdOut (&stCard, "Memory allocation error\n", false);
    spcm_dwDefTransfer_i64 (stCard.hDrv, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, (uint32) llNotifySize, pvBuffer, 0, llSWBufSize);

    // we now start everything
    if (!stCard.bSetError)
        dwErr = spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_DATA_STARTDMA);



    // ------------------------------------------------------------------------
    // this is our main polling loop
    printf ("\n\n");
    if (!stCard.bSetError)
        {
        int64 llTransferredBytes = 0;
        int64 llAvailUser;
        int64 llBufferFillPromille;
        int32 lStatus;

        // we loop as long as everthing is ok or an overrun has occurred
        while ((dwErr == ERR_OK) || (dwErr == ERR_FIFOHWOVERRUN))
            {

            // check status and fill size of buffers
            if (!dwErr) dwErr = spcm_dwGetParam_i32 (stCard.hDrv, SPC_M2STATUS, &lStatus);
            if (!dwErr) dwErr = spcm_dwGetParam_i64 (stCard.hDrv, SPC_DATA_AVAIL_USER_LEN, &llAvailUser);
            if (!dwErr) dwErr = spcm_dwGetParam_i64 (stCard.hDrv, SPC_FILLSIZEPROMILLE, &llBufferFillPromille);

            // give some status messages
            printf("\r[");
            if (lStatus & M2STAT_CARD_PRETRIGGER)   printf (" Pre");  else printf ("    ");
            if (lStatus & M2STAT_CARD_TRIGGER)      printf (" Trg");  else printf ("    ");
            if (lStatus & M2STAT_CARD_TRIGGER)      printf (" Rdy");  else printf ("    ");
            if (lStatus & M2STAT_DATA_BLOCKREADY)   printf (" DBlk"); else printf ("     ");
            if (lStatus & M2STAT_DATA_END)          printf (" DEnd"); else printf ("     ");
            if (lStatus & M2STAT_DATA_OVERRUN)      printf (" DOvr"); else printf ("     ");
            printf ("] ");

            // some data has been transferred
            if (llAvailUser > 0)
                {
                llTransferredBytes += llAvailUser;

                // this is the point to do something with the data, we simple give a message here
                printf ("Avail:%6dkB  SW:%3.0f%%  HW:%3d%%  Total:%6.2fMB", (int32) (llAvailUser / KILO_B(1)), (float) 100.0 * llAvailUser / llSWBufSize, (uint32) llBufferFillPromille / 10, (float) llTransferredBytes / MEGA_B(1));

                // set the data free for FIFO mode again
                dwErr = spcm_dwSetParam_i64 (stCard.hDrv, SPC_DATA_AVAIL_CARD_LEN, llAvailUser);
                }

            // check for data end, being emerged after an overrun and all the remaining data has been read out
            // the loop continues after the overrun to get all data that is still in the on-board memory
            if (lStatus & M2STAT_DATA_END)
                {
                dwErr = ERR_ABORT;
                printf ("\nOverrun and all remaining data has been read out, we quit now!\n");
                }

            // check for esc=abort
            if (!dwErr)
                if (bKbhit())
                    if (cGetch() == 27)
                        {
                        spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_STOP);
                        dwErr = ERR_ABORT;
                        printf ("\n\n");
                        }
            }
        }

    

    // ------------------------------------------------------------------------
    // print error information if an error occured
    if (stCard.bSetError)
        return nSpcMErrorMessageStdOut (&stCard, "An error occured while programming the card:\n", true);


    // clean up and close the driver
    vSpcMCloseCard (&stCard);

    printf ("\n\n");
    return EXIT_SUCCESS;
    }

