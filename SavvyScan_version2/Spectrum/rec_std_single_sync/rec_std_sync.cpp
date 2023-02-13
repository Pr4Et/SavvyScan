/*
**************************************************************************

rec_std_sync.cpp                                         (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based acquisition cards. 

Information about the different products and their drivers can be found
online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/platform-driver-and-series-differences

Shows the synchronization of one or more cards with acquisition mode

The Star-Hub is accessed directly as this is very simple
  
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

// ----- include of common example librarys -----
#include "../common/spcm_lib_card.h"
#include "../common/spcm_lib_data.h"

// ----- standard c include files -----
#include <stdio.h>
#include <stdlib.h>



/*
**************************************************************************
vDoCardSetup
**************************************************************************
*/

void vDoCardSetup (ST_SPCM_CARDINFO *pstCard)
    {
    int i;

    // standard single, all channels, memsize=16k, posttrigge=8k -> pretrigger=8k
    bSpcMSetupModeRecStdSingle (pstCard, pstCard->lMaxChannels == 64 ? -1 : (1 << pstCard->lMaxChannels) - 1, 16 * 1024, 8 * 1024);

    // we try to set the samplerate to 1 MHz (M2i) or 20 MHz (M3i and M4i) on internal PLL, no clock output
    if (pstCard->bM2i)
        bSpcMSetupClockPLL (pstCard, MEGA(1), false);
    else if (pstCard->bM2p)
        bSpcMSetupClockPLL (pstCard, MEGA(10), false);
    else
        bSpcMSetupClockPLL (pstCard, MEGA(20), false);
    printf ("Sampling rate card %d set to %.1lf MHz\n", pstCard->lCardIdx, (double) pstCard->llSetSamplerate / MEGA (1));

    // type dependent card setup
    switch (pstCard->eCardFunction)
        {

        // analog acquisition card setup
        case AnalogIn:

            // program all input channels to +/-1 V and 50 ohm termination (if it's available)
            for (i=0; i < pstCard->lMaxChannels; i++)
                if (pstCard->bM2i || pstCard->bM2p)
                    bSpcMSetupInputChannel (pstCard, i, 1000, true);
                else if (pstCard->bM3i || pstCard->bM4i)
                    bSpcMSetupPathInputCh (pstCard, i, 0, 1000, 0, false, true, true); 
            break;

        // digital acquisition card setup
        case DigitalIn:
        case DigitalIO:
            printf ("Not yet implemented\n");
            break;
        }
    }




/*
**************************************************************************
main 
**************************************************************************
*/

int main ()
    {
    char                szBuffer[1024];     // a character buffer for any messages
    ST_SPCM_CARDINFO*   pstCard[MAXBRD];    // info structure of my card
    int16               nCardCount, nIdx;
    int16               nStarHubCarrierIdx = 0;
    bool                bStarhubFound = false;
    drv_handle          hSync;
    bool                bOk = true;
    int32               lSyncdCards;



    // ------------------------------------------------------------------------
    // init all cards in the system
    for (nCardCount = 0; nCardCount < MAXBRD; nCardCount++)
        {
        pstCard[nCardCount] = new ST_SPCM_CARDINFO;

        // uncomment the second line and replace the IP address to use remote
        // cards like in a digitizerNETBOX
        if (!bSpcMInitCardByIdx (pstCard[nCardCount], nCardCount))
        //if (!bSpcMInitCardByIdx (pstCard[nCardCount], "192.168.1.10", nCardCount))
            {

            // no card found
            if (!nCardCount)
                return nSpcMErrorMessageStdOut (pstCard[nCardCount], "Error: Could not open card\n", true);

            // clean up
            delete pstCard[nCardCount];
            break;
            }

        // check for star hub on the card
        else
            {
            if ((pstCard[nCardCount]->lFeatureMap & SPCM_FEAT_STARHUB5) || (pstCard[nCardCount]->lFeatureMap & SPCM_FEAT_STARHUB16))
                {
                bStarhubFound = true;
                nStarHubCarrierIdx = nCardCount;
                }
            }
        }

    printf (pszSpcMPrintDocumentationLink (pstCard[0], szBuffer, sizeof (szBuffer)));

    // print info of all cards
    for (nIdx = 0; nIdx < nCardCount; nIdx++)
        printf ("Card %d\n%s\n\n", nIdx, pszSpcMPrintCardInfo (pstCard[nIdx], szBuffer, sizeof (szBuffer)));

    //  not our example if there's no starhub
    if (!bStarhubFound)
        {
        printf ("\nThere's no starhub in the system, this example can't run\n");
        bOk = false;
        }

    // the star hub is accessed under it's own handle
    if (bOk)
        {
        hSync = spcm_hOpen ("sync0");
        if (!hSync)
            {
            printf ("\nCan't open starhub handle\n");
            bOk = false;
            }

        spcm_dwGetParam_i32 (hSync, SPC_SYNC_READ_SYNCCOUNT, &lSyncdCards);
        }



    // ------------------------------------------------------------------------
    // show cable connection info
    if (bOk)
        {
        printf ("Star-Hub information:\n");
        printf ("star-hub is connected with %d cards\n", lSyncdCards);
        for (nIdx = 0; nIdx < lSyncdCards; nIdx++)
            {
            int32 lCable;
            spcm_dwGetParam_i32 (hSync, SPC_SYNC_READ_CABLECON0 + nIdx, &lCable);
            printf ("   Card Idx %d (sn %05d) is", nIdx, pstCard[nIdx]->lSerialNumber);
            if (lCable != -1) 
                printf (" connected on cable %d\n", lCable);
            else
                printf (" not connected with the star-hub\n");
            }
        printf ("\n");
        }



    // ------------------------------------------------------------------------
    // setup
    if (bOk)
        {

        // all cards got a similar setup, trigger sources disabled as default
        for (nIdx = 0; (nIdx < nCardCount) && bOk; nIdx++)
            {
            vDoCardSetup (pstCard[nIdx]);
            bOk = bSpcMSetupTrigMask (pstCard[nIdx], 0, 0, 0, 0, 0, 0);
            }

        // 1st card is used as trigger master (un-comment the second line to have external trigger on card 0)
        bOk =bSpcMSetupTrigSoftware (pstCard[0], false);
        //bOk =bSpcMSetupTrigExternal (pstCard[0], SPC_TM_POS, false, 0);
        }

    // and run
    if (bOk)
        {
        uint32 dwError = ERR_OK;

        // sync setup, all card activated, last card is clock master
        if (!dwError) dwError = spcm_dwSetParam_i32 (hSync, SPC_SYNC_ENABLEMASK, (1 << nCardCount) - 1);
        if (!dwError) dwError = spcm_dwSetParam_i32 (hSync, SPC_SYNC_CLKMASK, (1 << nStarHubCarrierIdx));

        // start the cards and wait for ready with a timeout of 5 seconds (5000 ms)
        if (!dwError) dwError = spcm_dwSetParam_i32 (hSync, SPC_TIMEOUT, KILO(5));
        if (!dwError)
            {
            printf ("Acquisition started for all cards\n");
            dwError = spcm_dwSetParam_i32 (hSync, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY);
            if (dwError == ERR_TIMEOUT)
                printf ("... Timeout\n");

            else if (!dwError)
                {
                printf ("... Sucessfully finished\n");

                // ----- read data from the cards -----
                void* apvBuffer[MAXBRD] = { NULL };
                for (nIdx = 0; (nIdx < nCardCount); nIdx++)
                    {
                    uint64 qwBufferLen_samples = pstCard[nIdx]->llSetMemsize * pstCard[nIdx]->lSetChannels;
                    uint64 qwBufferLen_bytes   = qwBufferLen_samples * pstCard[nIdx]->lBytesPerSample;
                    if (pstCard[nIdx]->lBytesPerSample == 1)
                        apvBuffer[nIdx] = new int8[(size_t)qwBufferLen_samples];
                    else
                        apvBuffer[nIdx] = new int16[(size_t)qwBufferLen_samples];
                    spcm_dwDefTransfer_i64 (pstCard[nIdx]->hDrv, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, 0, apvBuffer[nIdx], 0, qwBufferLen_bytes);

                    // start DMA transfer and wait until it is finished
                    spcm_dwSetParam_i32 (hSync, SPC_TIMEOUT, 0);
                    spcm_dwSetParam_i32(pstCard[nIdx]->hDrv, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);
                    }

                // now all data from each card is in the corresponding buffer
                // do something with the data
                // ...

                // after we have processed the data, we delete the buffers again
                for (nIdx = 0; (nIdx < nCardCount); nIdx++)
                    {
                    if (pstCard[nIdx]->lBytesPerSample == 1)
                        delete [] static_cast < int8* > (apvBuffer[nIdx]);
                    else
                        delete [] static_cast < int16* > (apvBuffer[nIdx]);
                    }
                }
            }

        // error message if something went wrong
        if (dwError && (dwError != ERR_TIMEOUT))
            {
            char szErrorMsg[ERRORTEXTLEN];
            spcm_dwGetErrorInfo_i32 (hSync, NULL, NULL, szErrorMsg);
            printf ("\nError:\n%s\n", szErrorMsg);
            }
        }




    // ------------------------------------------------------------------------
    // clean up and close the driver
    if (bOk && hSync)
        spcm_vClose (hSync);

    for (nIdx = 0; nIdx < nCardCount; nIdx++)
        {
        vSpcMCloseCard (pstCard[nIdx]);
        delete (pstCard[nIdx]);
        }

    return EXIT_SUCCESS;
    }

