/*
**************************************************************************

rec_multi_binaryfile.cpp                             (c) Spectrum GmbH

**************************************************************************

this example supports all SpcmDrv based acquisition cards with the option
multiple recording installed. If timestamp is installed the timestamps are also
read.

Information about the different products and their drivers can be found
online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/platform-driver-and-series-differences

Does a continous FIFO transfer and writes data to a binary file. The file
type can be selected to simple binary file, binary file including the
timestamps at the beginning of each segment or binary file including the
SBench5 format header (*.sbs/*.sb5) file. 

Change the eFileType variable to select between the different file types 

Change the global flag g_bThread to use the threaded version or the plain
and more simple loop.

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

// ----- include of common example librarys -----
#include "../common/spcm_lib_card.h"
#include "../common/spcm_lib_data.h"
#include "../common/spcm_lib_thread.h"
#include "../sb5_file/sb5_file.h"


// ----- this is the global thread flag that defines whether we use the thread or non-thread loop -----
bool g_bThread = true;



/*
**************************************************************************
Working routine data
**************************************************************************
*/

#define FILENAME "stream_test"

typedef enum E_FILETYPE {eFT_noWrite, eFT_PlainBinary, eFT_PlainWithTimestamps, eFT_SB5_Stream} E_FILETYPE;

struct ST_WORKDATA
    {
    E_FILETYPE  eFileType;
    int64       llSamplesWritten;
    FILE*       hFile;
    char        szFileName[100];
    int32       lSegmentsize;
    };



/*
**************************************************************************
bDoCardSetup: setup matching the calculation routine
**************************************************************************
*/

bool bDoCardSetup (ST_WORKDATA* pstWorkData, ST_BUFFERDATA* pstBufferData)
    {
    int     i;

    // FIFO mode setup, we run continuously and use 128 samples of pretrigger for each segment
    pstWorkData->lSegmentsize = KILO_B(1);          // segment size
    pstWorkData->eFileType =    eFT_noWrite;        // storage mode

    // we try to set the samplerate to 1 MHz (M2i) or 20 MHz (M3i, M4i) on internal PLL, no clock output
    // increase this to test the read-out-after-overrun
	if (pstBufferData->pstCard->bM2i)
        bSpcMSetupClockPLL (pstBufferData->pstCard, MEGA(1), false);
	else if (pstBufferData->pstCard->bM2p)
        bSpcMSetupClockPLL (pstBufferData->pstCard, MEGA(10), false);
	else if (pstBufferData->pstCard->bM3i || pstBufferData->pstCard->bM4i)
        bSpcMSetupClockPLL (pstBufferData->pstCard, MEGA(20), false);
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstBufferData->pstCard->llSetSamplerate / 1000000);

    // we set external trigger for multiple recording
    bSpcMSetupTrigExternal (pstBufferData->pstCard, SPC_TM_POS, false, 0);

    // type dependent card setup
    switch (pstBufferData->pstCard->eCardFunction)
        {

        // analog acquisition card setup
        case AnalogIn:

            // we only enable 1 channel for the example
            bSpcMSetupModeRecFIFOMulti (pstBufferData->pstCard, CHANNEL0, pstWorkData->lSegmentsize, pstWorkData->lSegmentsize - 128);

            // program all input channels to +/-1 V and 50 ohm termination (if it's available)
            if (pstBufferData->pstCard->bM2i || pstBufferData->pstCard->bM2p)
                {
                for (i=0; i < pstBufferData->pstCard->lMaxChannels; i++)
                    bSpcMSetupInputChannel (pstBufferData->pstCard, i, 1000, true);
                }
            else
                {
                bool bTerm           = true;
                bool bACCoupling     = false;
                bool bBandwidthLimit = false;
                for (i=0; i < pstBufferData->pstCard->lMaxChannels; i++)
                    bSpcMSetupPathInputCh (pstBufferData->pstCard, i, 0, 1000, 0, bTerm, bACCoupling, bBandwidthLimit);
                }
            break;

        // digital acquisition card setup
        case DigitalIn:
        case DigitalIO:

            // we enable 16 channels for the example
            bSpcMSetupModeRecFIFOMulti (pstBufferData->pstCard, 0xffff, pstWorkData->lSegmentsize, pstWorkData->lSegmentsize - 128);
            
            // set all input channel groups to 110 ohm termination (if it's available)
            for (i=0; i < pstBufferData->pstCard->uCfg.stDIO.lGroups; i++)
                bSpcMSetupDigitalInput (pstBufferData->pstCard, i, true);
            break;
        }

    // set up the timestamp mode to standard if timestamp is installed
    if (pstBufferData->bStartTimestamp)
        bSpcMSetupTimestamp (pstBufferData->pstCard, SPC_TSMODE_STANDARD | SPC_TSCNT_INTERNAL, 0);

    return pstBufferData->pstCard->bSetError;
    }



/*
**************************************************************************
Setup working routine
**************************************************************************
*/

bool bWorkSetup (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA* pstWorkData = (ST_WORKDATA*) pvWorkData;

    // setup for the transfer, to avoid overrun one may use even larger data blocks
    pstBufferData->dwDataBufLen =   MEGA_B(8);
    pstBufferData->dwDataNotify =   KILO_B(512);
    pstBufferData->dwTSBufLen =     MEGA_B(1);
    pstBufferData->dwTSNotify =     KILO_B(4);

    // setup for the work
    pstWorkData->llSamplesWritten = 0;

    switch (pstWorkData->eFileType)
        {
        case eFT_PlainBinary:
        case eFT_PlainWithTimestamps:
            sprintf (pstWorkData->szFileName, "%s.bin", FILENAME);
            printf ("\nWriting to plan binary file %s\n\n", pstWorkData->szFileName);
            break;

        case eFT_SB5_Stream:
            if (pstBufferData->pstCard->lSetChannels == 1)
                sprintf (pstWorkData->szFileName, "%s.sb5", FILENAME);
            else
                sprintf (pstWorkData->szFileName, "%s.sbs", FILENAME);
            printf ("\nWriting to SBench 5 binary file %s\n\n", pstWorkData->szFileName);
            break;

        case eFT_noWrite:
            printf ("\nJust simulating, no real write\n\n");
            break;
        }

    // check some details if we're storing segments together with timestamps 
    if (pstWorkData->eFileType == eFT_PlainWithTimestamps)
        {
        uint32 dwSegments = pstBufferData->dwDataBufLen / pstWorkData->lSegmentsize / pstBufferData->pstCard->lSetChannels / pstBufferData->pstCard->lBytesPerSample;

        // check whether timestamps are installed, otherwise this mode doesn't make sense
        if (!pstBufferData->bStartTimestamp)
            {
            printf ("\nThis storing mode needs the option timestamp installed. it doesn't work with your card\n");
            return false;
            }

        // a full number of segments have to fit in our buffer otherwise the algorithm in this example won't work
        if (((uint32) (pstBufferData->dwDataBufLen / dwSegments) * dwSegments) != pstBufferData->dwDataBufLen)
            {
            printf ("\nFor this storing mode a full number of segments must fir into the buffer. Please correct setup\n");
            return false;
            }
        }

    if (pstWorkData->eFileType != eFT_noWrite)
        pstWorkData->hFile = fopen (pstWorkData->szFileName, "w+b");

    // we now have to write the SB5 header if this format has been selected
    if (pstWorkData->hFile && (pstWorkData->eFileType == eFT_SB5_Stream))
        {
        ST_SB5HEAD* pstHeader;
        bool        bReturn;

        pstHeader = pstSB5_AllocHeader (pstBufferData->pstCard->lSetChannels);
        if (!bFillSB5Header (pstBufferData->pstCard, pstHeader, "Test"))
            return false;

        bReturn = bSB5_StoreHeader (pstWorkData->hFile, pstHeader);
        vSB5_FreeHeader (pstHeader);

        return bReturn;
        }

    return ((pstWorkData->hFile != NULL) || (pstWorkData->eFileType == eFT_noWrite));
    }



/*
**************************************************************************
bWorkDo: stores data to hard disk
**************************************************************************
*/

bool bWorkDo (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA*    pstWorkData = (ST_WORKDATA*) pvWorkData;
    uint32          dwWritten, dwTimestampBytes;

    // ----- M2i+M3i timestamps are 64 bit, M4i 128 bit -----
    dwTimestampBytes = 0;
    if (pstBufferData->pstCard->bM2i || pstBufferData->pstCard->bM3i)
        dwTimestampBytes = 8;
    else if (pstBufferData->pstCard->bM4i || pstBufferData->pstCard->bM2p)
        dwTimestampBytes = 16;

    // ----- data together with timestamps -----
    if (pstWorkData->eFileType == eFT_PlainWithTimestamps)
        {

        // if we write timestamps together with data we need enough data from both here
        uint32 dwSegmentInBytes = (pstWorkData->lSegmentsize * pstBufferData->pstCard->lBytesPerSample * pstBufferData->pstCard->lSetChannels);
        uint64 qwDataAvail = pstBufferData->llDataAvailBytes;
        uint64 qwTSAvail = pstBufferData->llTSAvailBytes;

        // at least one timestamp and one segment is there
        while ((qwTSAvail >= 8) && (qwDataAvail >= dwSegmentInBytes))
            {

            // write the current timestamp and then the data 
            dwWritten = fwrite (pstBufferData->pvTSCurrentBuf, 1, dwTimestampBytes, pstWorkData->hFile);

            if (dwWritten == dwTimestampBytes)
                dwWritten = fwrite (pstBufferData->pvDataCurrentBuf, 1, dwSegmentInBytes, pstWorkData->hFile);
            if (dwWritten != dwSegmentInBytes)
                {
                printf ("\nData Write error\n");
                return false;
                }

            // count the timestamp and recalc sizes and positions
            pstWorkData->llSamplesWritten++;
            qwTSAvail -=                        dwTimestampBytes;
            pstBufferData->pvTSCurrentBuf  =    (void*) (((char*) pstBufferData->pvTSCurrentBuf) + dwTimestampBytes);
            qwDataAvail -=                      dwSegmentInBytes;
            pstBufferData->pvDataCurrentBuf =   (void*) (((char*) pstBufferData->pvDataCurrentBuf) + dwSegmentInBytes);
            }

        // subtract the bytes that we've not used so far as we still need the rest next time
        pstBufferData->llDataAvailBytes -=  qwDataAvail;
        pstBufferData->llTSAvailBytes -=    qwTSAvail;    

        // print the status
        printf ("\r%.2f MSamples written to %s, %d segments in total",
            (double) pstBufferData->qwDataTransferred / pstBufferData->pstCard->lBytesPerSample / 1024.0 / 1024.0, 
            pstWorkData->szFileName,
            (int32) pstWorkData->llSamplesWritten);

        return true;
        }



    // ----- simple data, timestamps are not recorded or ignored -----
    if (pstWorkData->eFileType != eFT_noWrite)
        {

        // we limit the data to write to chunks of notify size to avoid a jam in the fwrite function
        if (pstBufferData->llDataAvailBytes > pstBufferData->dwDataNotify)
            pstBufferData->llDataAvailBytes = pstBufferData->dwDataNotify;

        // write the data and count the samples
        dwWritten = fwrite (pstBufferData->pvDataCurrentBuf, 1, (size_t)pstBufferData->llDataAvailBytes, pstWorkData->hFile);
        pstWorkData->llSamplesWritten += dwWritten / pstBufferData->pstCard->lBytesPerSample;
        if (dwWritten != pstBufferData->llDataAvailBytes)
            {
            printf ("\nData Write error\n");
            return false;
            }

        // announce the number of data that has been written
        printf ("\r%.2f MSamples written to %s, %.2fk Timestamps received", 
            (double) pstBufferData->qwDataTransferred / pstBufferData->pstCard->lBytesPerSample / 1024.0 / 1024.0, 
            pstWorkData->szFileName,
            (double) pstBufferData->qwTSTransferred / dwTimestampBytes / 1024.0);
        }

    // ----- no data written, just announce the transferred data -----
    if (pstWorkData->eFileType == eFT_noWrite)
        printf ("\r%.2f MSamples written, %.2fk Timestamps received", 
            (double) pstBufferData->qwDataTransferred / pstBufferData->pstCard->lBytesPerSample / 1024.0 / 1024.0, 
            (double) pstBufferData->qwTSTransferred / dwTimestampBytes / 1024.0);

    return true;
    }




/*
**************************************************************************
vWorkClose: Close the work and clean up
**************************************************************************
*/

void vWorkClose (void* pvWorkData, ST_BUFFERDATA* pstBufferData)
    {
    ST_WORKDATA* pstWorkData = (ST_WORKDATA*) pvWorkData;

    if (pstWorkData->eFileType != eFT_noWrite)
        if (pstWorkData->hFile)
            {
            if (pstWorkData->eFileType == eFT_SB5_Stream)
                bSB5_UpdateSamples (pstWorkData->hFile, (int32) pstWorkData->llSamplesWritten);

            fclose (pstWorkData->hFile);
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
    ST_SPCM_CARDINFO    stCard;             // info structure of my card
    ST_WORKDATA         stWorkData;
    ST_BUFFERDATA       stBufferData;       // buffer and start information

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

    memset (&stBufferData, 0, sizeof(stBufferData));
    stBufferData.pstCard = &stCard;


    // check whether we support this card type in the example
    if ((stCard.eCardFunction != AnalogIn) && (stCard.eCardFunction != DigitalIn) && (stCard.eCardFunction != DigitalIO))
        return nSpcMErrorMessageStdOut (&stCard, "Error: Card function not supported by this example\n", false);
    if ((stCard.lFeatureMap & SPCM_FEAT_MULTI) == 0)
        return nSpcMErrorMessageStdOut (&stCard, "Error: Multiple Recording Option not installed. Examples was done especially for this option!\n", false);

    // set a flag if timestamp is installed
    stBufferData.bStartTimestamp = ((stCard.lFeatureMap & SPCM_FEAT_TIMESTAMP) != 0);



    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
    if (!stCard.bSetError)
        bDoCardSetup (&stWorkData, &stBufferData);



    // ------------------------------------------------------------------------
    // some additional information on the acquisition
    if (!stCard.bSetError)
        {
        printf ("\nData information:\n=================\n");
        printf ("Each segment is %.3f ms long\n", 1000.0 * stWorkData.lSegmentsize / stCard.llSetSamplerate);
        printf ("Maximum pulse repetition frequency to reach with this setting is %.2f Hz\n", (double) stCard.llSetSamplerate / stWorkData.lSegmentsize);
        }



    // ------------------------------------------------------------------------
    // setup the data transfer and start it
    stBufferData.bStartCard =       true;
    stBufferData.bStartData =       true;
    stBufferData.bStartTimestamp =  ((stCard.lFeatureMap & SPCM_FEAT_TIMESTAMP) != 0);

    // start the threaded loop
    stBufferData.lTimeout =         5000;
    if (!stCard.bSetError && g_bThread)
        vDoThreadMainLoop (&stBufferData, &stWorkData, bWorkSetup, bWorkDo, vWorkClose, bKeyAbortCheck);

    // this is the non threaded loop, we use a small timeout of 100 ms here as we otherwise won't return from the loop if no dat is coming
    stBufferData.lTimeout =         100;
    if (!stCard.bSetError && !g_bThread)
        {
        stBufferData.bStartExtraDMA = stBufferData.bStartTimestamp;
        vDoMainLoop (&stBufferData, &stWorkData, bWorkSetup, bWorkDo, vWorkClose, bKeyAbortCheck);
        }



    // ------------------------------------------------------------------------
    // print error information if an error occured
    if (stCard.bSetError)
        return nSpcMErrorMessageStdOut (&stCard, "An error occured while programming the card:\n", true);


    // clean up and close the driver
    vSpcMCloseCard (&stCard);

    return EXIT_SUCCESS;
    }

