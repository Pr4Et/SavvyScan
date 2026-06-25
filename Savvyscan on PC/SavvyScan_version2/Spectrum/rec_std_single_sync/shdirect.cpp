/*
**************************************************************************

shdirect.cpp                                    (c) Spectrum GmbH, 10/2018

**************************************************************************

Example for all SpcMDrv based (M2i, M2i-Express and M2p) cards. The example works 
with two cards from the same series, one of them must have a star-hub installed.

Shows the synchronization of two cards using the special SH-Direct clock mode.
In that mode one card (or one set of cards) is running in a normal synchronisation
and another card (or several cards) just take the clock from the star-hub as an
external clock. Besides using this clock the cards run totally independent. 

This mode can for example be used to run a waveform or pattern generator 
continuously and run a data acquisition card with the same speed but start
it several times with different setup all while the generator card is still
running.

This example expects the star hub to be present on a generator card!
  
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

void vDoGeneratorCardSetup (ST_SPCM_CARDINFO *pstCard)
    {
    int i;

    // we try to set the samplerate to 1 MHz on internal PLL with enabled clock output
    bSpcMSetupClockPLL (pstCard, MEGA (1), true);

    // type dependent card setup
    switch (pstCard->eCardFunction)
        {

        // analog acquisition card setup
        case AnalogOut:

            // continuous mode with one channel
            bSpcMSetupModeRepStdLoops (pstCard, CHANNEL0, KILO_B(16), 0);

            // program all analog outputs to +/-1 V, no offset, no filter
            for (i=0; i < pstCard->lMaxChannels; i++)
                bSpcMSetupAnalogOutputChannel (pstCard, i, 1000);
            break;

        // digital output card setup
        case DigitalOut:
        case DigitalIO:

            // continuous mode with 16 channels
            bSpcMSetupModeRepStdLoops (pstCard, 0xffff, KILO_B(16), 0);
            bSpcMSetupDigitalOutput (pstCard, 0);
            break;
        }

    // software trigger for start
    bSpcMSetupTrigSoftware (pstCard);
    }

// ***********************************************************************

void vDoAcquisitionCardSetup (ST_SPCM_CARDINFO *pstCard)
    {
    int i;

    // the card is running with the direct star-hub clock here! (same clock or diveded from master)
    spcm_dwSetParam_i32 (pstCard->hDrv, SPC_CLOCKMODE, SPC_CM_SHDIRECT);
    spcm_dwSetParam_i64 (pstCard->hDrv, SPC_SAMPLERATE, KILO(500));
    spcm_dwSetParam_i32 (pstCard->hDrv, SPC_CLOCKOUT, 1);

    // type dependent card setup
    switch (pstCard->eCardFunction)
        {

        // analog acquisition card setup
        case AnalogIn:

            // single shot mode with one channel
            bSpcMSetupModeRecStdSingle (pstCard, CHANNEL0, KILO_B(16), KILO_B(8));

            // program all analog inputs to +/-1 V, 50 ohm termination
            for (i=0; i < pstCard->lMaxChannels; i++)
                bSpcMSetupInputChannel (pstCard, i, 1000, true);
            break;

        // digital input card setup
        case DigitalIn:
        case DigitalIO:

            // single shot mode with 16 channels
            bSpcMSetupModeRecStdSingle (pstCard, 0xffff, KILO_B(16), KILO_B(8));

            // termination active
            bSpcMSetupDigitalInput (pstCard, 0, true);
            break;
        }

    // software trigger for start
    bSpcMSetupTrigSoftware (pstCard);
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
    int16               nCardCount, nIdx, nStarHubCard = 0, nAcquisitionCard = 0;
    bool                bStarhubFound = false;
    drv_handle          hSync = NULL;
    bool                bOk = true;
    int32               lSyncdCards;
    uint32              dwError = ERR_OK;



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

        // check for star hub on the card and store the card index
        else
            if ((pstCard[nCardCount]->lFeatureMap & SPCM_FEAT_STARHUB5) || (pstCard[nCardCount]->lFeatureMap & SPCM_FEAT_STARHUB16))
                {
                bStarhubFound = true;
                nStarHubCard = nCardCount;
                }
            else
                nAcquisitionCard = nCardCount;
        }

    printf (pszSpcMPrintDocumentationLink (pstCard[0], szBuffer, sizeof (szBuffer)));

    // print info of all cards
    for (nIdx = 0; nIdx < nCardCount; nIdx++)
        printf ("Card %d   %s\n%s\n\n", nIdx, (nIdx == nStarHubCard? "(SH carrier)" : ""), pszSpcMPrintCardInfo (pstCard[nIdx], szBuffer, sizeof (szBuffer)));

    //  not our example if there's no starhub
    if (!bStarhubFound)
        {
        printf ("\nThere's no starhub in the system, this example can't run\n");
        bOk = false;
        }

    // the star-hub card must be a generator card for the example
    if (bOk)
        switch (pstCard[nStarHubCard]->eCardFunction)
            {
            case AnalogIn:
            case DigitalIn:
                printf ("\nThe example expects a generator card as the star-hub holding card\n");
                bOk = false;
                break;
            }

    // all other cards must be acquisition cards
    if (bOk)
        for (nIdx = 0; (nIdx < nCardCount) && bOk; nIdx++)
            if (nIdx != nStarHubCard)
                switch (pstCard[nIdx]->eCardFunction)
                    {
                    case AnalogOut:
                    case DigitalOut:
                        printf ("\nThe example expects acquisition cards for all cards that do not hold the star-hub\n");
                        bOk = false;
                        break;
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
    // setup and run all
    if (bOk)
        {

        // setup
        vDoGeneratorCardSetup (pstCard[nStarHubCard]);
        vDoAcquisitionCardSetup (pstCard[nAcquisitionCard]);

        // sync setup, only the star-hub card runs synchronised with itself
        if (!dwError) dwError = spcm_dwSetParam_i32 (hSync, SPC_SYNC_ENABLEMASK, (1 << nStarHubCard));
        if (!dwError) dwError = spcm_dwSetParam_i32 (hSync, SPC_SYNC_CLKMASK, (1 << nStarHubCard));

        // the star-hub card must be started first as it provides clock information for
        // the other cards!

        // start the generator card
        printf ("... start generator card\n");
        if (!dwError) dwError = spcm_dwSetParam_i32 (hSync, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER);

        // read out the sampling rates and show it
        int64 llTest;
        spcm_dwGetParam_i64 (pstCard[nStarHubCard]->hDrv, SPC_SAMPLERATE, &llTest);
        printf ("Sampling rate generator card:   %.2lf MS/s\n", (double) llTest / MEGA(1));
        spcm_dwGetParam_i64 (pstCard[nAcquisitionCard]->hDrv, SPC_SAMPLERATE, &llTest);
        printf ("Sampling rate acquisition card: %.2lf MS/s\n", (double) llTest / MEGA(1));
        printf ("\n");

        // we now start the acquisition card a couple of times, it runs with the clock
        // taken from the star-hub. To keep the example simple we only use one acquisition
        // card here
        for (int i=0; (i<10) && (dwError == ERR_OK); i++)
            {
            printf ("... start acquisition card and wait for ready\n");
            dwError = spcm_dwSetParam_i32 (pstCard[nAcquisitionCard]->hDrv, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY);
            if (dwError == ERR_OK)
                printf ("... ready\n");
            }
        }

    // error message if something went wrong
    if (dwError && (dwError != ERR_TIMEOUT))
        {
        char szErrorMsg[ERRORTEXTLEN];
        if (!spcm_dwGetErrorInfo_i32 (hSync, NULL, NULL, szErrorMsg))
            spcm_dwGetErrorInfo_i32 (pstCard[nAcquisitionCard]->hDrv, NULL, NULL, szErrorMsg);
        printf ("\nError:\n%s\n", szErrorMsg);
        }


    // ------------------------------------------------------------------------
    // clean up and close the driver
    if (hSync)
        spcm_vClose (hSync);

    for (nIdx = 0; nIdx < nCardCount; nIdx++)
        {
        vSpcMCloseCard (pstCard[nIdx]);
        delete (pstCard[nIdx]);
        }

    return EXIT_SUCCESS;
    }

