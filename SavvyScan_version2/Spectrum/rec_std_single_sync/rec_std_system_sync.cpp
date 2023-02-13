/*
**************************************************************************

rec_std_system_sync.cpp                                  (c) Spectrum GmbH

**************************************************************************

Example for all M2i acquisition cards. 

Shows the synchronization of one or more cards on each system
starhub with standard single acquisition mode.

For simplicity this example assumes that at least one
"system starhub master" and one "system starhub slave" are both
installed in the same PC system, to gain easy software
access to both devices, without the need for inter-system
software communication.

Such a setup is rather unlikely for real-world use, because
such setup would render the usage of a system starhub over a
standard starhub rather useless.
  
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
#define NO_OF_STARHUBS 2

int main ()
    {
    char                szBuffer[1024];     // a character buffer for any messages
    ST_SPCM_CARDINFO*   pstCard[MAXBRD];    // info structure of my card
    int16               nCardCount;
    int16               nStarhubCount;
    int16               nIdx;
    drv_handle          hSync[NO_OF_STARHUBS]; 
    bool                bOk = true;
    bool                bStarhubFound = false;
    bool                bSystemMasterStarhubFound = false;
    bool                bSystemSlaveStarhubFound = false;
    int32               plSyncedCards[NO_OF_STARHUBS];
    char                szSyncName[20];



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
        else
            {
            // check for star hub on the card
            if ((pstCard[nCardCount]->lFeatureMap & SPCM_FEAT_STARHUBSYSMASTER) || (pstCard[nCardCount]->lFeatureMap & SPCM_FEAT_STARHUB16))
                bStarhubFound = true;

            if (pstCard[nCardCount]->lFeatureMap & SPCM_FEAT_STARHUBSYSMASTER)
                bSystemMasterStarhubFound = true;

            if (pstCard[nCardCount]->lFeatureMap & SPCM_FEAT_STARHUBSYSSLAVE)
                bSystemSlaveStarhubFound = true;
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

    //  not our example if there's no system starhub master
    if (!bSystemMasterStarhubFound)
        {
        printf ("\nThere's no system starhub master in the system, this example can't run\n");
        bOk = false;
        }

    //  not our example if there's no system starhub slave
    if (!bSystemSlaveStarhubFound)
        {
        printf ("\nThere's no system starhub slave in the system, this example can't run\n");
        bOk = false;
        }


    // the star hubs are accessed by their own handles
    if (bOk)
        {        
        for (nStarhubCount = 0; nStarhubCount < NO_OF_STARHUBS; nStarhubCount++)
            {
            sprintf (szSyncName, "/dev/sync%d", nStarhubCount);
            hSync[nStarhubCount] = spcm_hOpen (szSyncName);
            if (!hSync[nStarhubCount])
                {
                printf ("\nCan't open starhub handle\n");
                bOk = false;
                }
            spcm_dwGetParam_i32 (hSync[nStarhubCount], SPC_SYNC_READ_SYNCCOUNT, &plSyncedCards[nStarhubCount]);
            }
        }



    // ------------------------------------------------------------------------
    // show cable connection info for all cards on all starhubs
    if (bOk)
        {
        for (nStarhubCount = 0; nStarhubCount < NO_OF_STARHUBS; nStarhubCount++)
            {
            printf ("Star-Hub %d information:\n", nStarhubCount);
            printf ("Star-hub is connected with %d cards\n", plSyncedCards[nStarhubCount]);
            for (nIdx = 0; nIdx < plSyncedCards[nStarhubCount]; nIdx++)
                {
                int32 lCable;
                spcm_dwGetParam_i32 (hSync[nStarhubCount], SPC_SYNC_READ_CABLECON0 + nIdx, &lCable);
                printf ("   Card Idx %d (sn %05d) is", nIdx, pstCard[nIdx]->lSerialNumber);
                if (lCable != -1) 
                    printf (" connected on cable %d\n", lCable);
                else
                    printf (" not connected with the star-hub\n");
                }
            printf ("\n\n");
            }
        }



    // ------------------------------------------------------------------------
    // setup
    if (bOk)
        {
        // all cards on all starhubs got a similar setup, trigger sources disabled as default
        for (nIdx = 0; (nIdx < nCardCount) && bOk; nIdx++)
            {
            vDoCardSetup (pstCard[nIdx]);
            bOk = bSpcMSetupTrigMask (pstCard[nIdx], 0, 0, 0, 0, 0, 0);
            }

        // 1st card is used as trigger master (un-comment the second line to have external trigger on card 0)
        bOk =bSpcMSetupTrigSoftware (pstCard[0], false);
        //bOk =bSpcMSetupTrigExternal (pstCard[0], SPC_TM_POS, false, 0);
        }

    // and run the required sequence for system starhub setups
    if (bOk)
        {
        uint32 dwError = ERR_OK;
        
        // configure all starhubs
        for (nStarhubCount = 0; nStarhubCount < NO_OF_STARHUBS; nStarhubCount++)
            {
            // sync setup, all cards on all starhubs activated
            if (!dwError) dwError = spcm_dwSetParam_i32 (hSync[nStarhubCount], SPC_SYNC_ENABLEMASK, (1 << plSyncedCards[nStarhubCount]) - 1);

            int32 lFeatures;
            spcm_dwGetParam_i32 (hSync[nStarhubCount], SPC_PCIFEATURES, &lFeatures);
            
            // last card on System-Star-Hub Master is set as clock master
            // And on each System-Star-Hub Slave also last card is set as clock master (necessary dummy)
            if (!dwError) dwError = spcm_dwSetParam_i32 (hSync[nStarhubCount], SPC_SYNC_CLKMASK, (1 << (plSyncedCards[nStarhubCount] - 1)));				
            
            // setup System-Star-Hubs to synchronize clock and trigger
            if (!dwError) dwError = spcm_dwSetParam_i32 (hSync[nStarhubCount], SPC_SYNC_MODE, SPC_SYNC_SYSTEMCLOCKTRIG);

            // define delay compensation (if required): set to default value of 4 here
            if (!dwError) dwError = spcm_dwSetParam_i32 (hSync[nStarhubCount], SPC_SYNC_SYSTEM_TRIGADJUST, 4);

            // configure a timeout of 5 seconds (5000 ms) for all cards, to avoid program deadlocks
            // in case that starting either master or slave(s) does not properly work
            if (!dwError) dwError = spcm_dwSetParam_i32 (hSync[nStarhubCount], SPC_TIMEOUT, KILO(5));
            }


        // transfer setup to sytsem master starhub card to have sampling clocks active before starting the slaves
        for (nStarhubCount = 0; nStarhubCount < NO_OF_STARHUBS; nStarhubCount++)
            {
            int32 lFeatures;
            spcm_dwGetParam_i32 (hSync[nStarhubCount], SPC_PCIFEATURES, &lFeatures);

            if (lFeatures & SPCM_FEAT_STARHUBSYSMASTER)
                {
                if (!dwError) dwError = spcm_dwSetParam_i32 (hSync[nStarhubCount], SPC_M2CMD, M2CMD_CARD_WRITESETUP);

                // leave loop once we found the master
                break;
                }
            }
        

        printf("\n... Starting all starhubs\n");


        // start all slave starhubs first
        for (nStarhubCount = 0; nStarhubCount < NO_OF_STARHUBS; nStarhubCount++)
            {
            int32 lFeatures;
            spcm_dwGetParam_i32 (hSync[nStarhubCount], SPC_PCIFEATURES, &lFeatures);
            
            // start all slaves first (but not the master itself)
            if ((lFeatures & SPCM_FEAT_STARHUBSYSSLAVE) && !(lFeatures & SPCM_FEAT_STARHUBSYSMASTER))
                {
                if (!dwError)
                    {
                    printf ("Acquisition started for slave Star-Hub%d\n", nStarhubCount);
                    dwError = spcm_dwSetParam_i32 (hSync[nStarhubCount], SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER);              
                    }

                // error message if something went wrong
                if (dwError && (dwError != ERR_TIMEOUT))
                    {
                    char szErrorMsg[ERRORTEXTLEN];
                    spcm_dwGetErrorInfo_i32 (hSync[nStarhubCount], NULL, NULL, szErrorMsg);
                    printf ("\nError:\n%s\n", szErrorMsg);
                    }
                }
            }


        // finally start master starhub
        for (nStarhubCount = 0; nStarhubCount < NO_OF_STARHUBS; nStarhubCount++)
            {
            int32 lFeatures;
            spcm_dwGetParam_i32 (hSync[nStarhubCount], SPC_PCIFEATURES, &lFeatures);
            
            // start all master last
            if (lFeatures & SPCM_FEAT_STARHUBSYSMASTER)
                {
                if (!dwError)
                    {
                    printf ("Acquisition started for master Star-Hub%d\n", nStarhubCount);
                    dwError = spcm_dwSetParam_i32 (hSync[nStarhubCount], SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY);
                    if (dwError == ERR_TIMEOUT)
                        printf ("... Master Timeout\n");
                    }

                // error message if something went wrong
                if (dwError && (dwError != ERR_TIMEOUT))
                    {
                    char szErrorMsg[ERRORTEXTLEN];
                    spcm_dwGetErrorInfo_i32 (hSync[nStarhubCount], NULL, NULL, szErrorMsg);
                    printf ("\nError:\n%s\n", szErrorMsg);
                    }
                
                // leave loop once we started the master
                break;
                }
            }

        // check if all slaves are ready
        for (nStarhubCount = 0; nStarhubCount < NO_OF_STARHUBS; nStarhubCount++)
            {
            int32 lFeatures;
            int32 lStatus;
            spcm_dwGetParam_i32 (hSync[nStarhubCount], SPC_PCIFEATURES, &lFeatures);
            
            // read status from all slaves (but not the master itself)
            if ((lFeatures & SPCM_FEAT_STARHUBSYSSLAVE) && !(lFeatures & SPCM_FEAT_STARHUBSYSMASTER))
                {
                do 
                    {
                    dwError = spcm_dwGetParam_i32 (hSync[nStarhubCount], SPC_M2STATUS, &lStatus);
                    }
                while (!(lStatus & M2STAT_CARD_READY) && !dwError);

                // error message if something went wrong
                if (dwError && (dwError != ERR_TIMEOUT))
                    {
                    char szErrorMsg[ERRORTEXTLEN];
                    spcm_dwGetErrorInfo_i32 (hSync[nStarhubCount], NULL, NULL, szErrorMsg);
                    printf ("\nError:\n%s\n", szErrorMsg);
                    }
                }
            }

        if (!dwError)
            {
            //
            // this would be the point to read data from all the cards on all starhubs
            //
            printf ("... Sucessfully finished\n");
            }
        }



    // ------------------------------------------------------------------------
    // clean up and close the driver
    for (nStarhubCount = 0; nStarhubCount < NO_OF_STARHUBS; nStarhubCount++)
        {
        if (bOk && hSync[nStarhubCount])
            spcm_vClose (hSync[nStarhubCount]);
        }

    for (nIdx = 0; nIdx < nCardCount; nIdx++)
        {
        vSpcMCloseCard (pstCard[nIdx]);
        delete (pstCard[nIdx]);
        }

    return EXIT_SUCCESS;
    }

