/*
**************************************************************************

rep_sequence.cpp                                         (c) Spectrum GmbH

**************************************************************************

Example for all M2i, M4i, M4x, M2p analog and digital generator cards.
Shows sequence replay mode as simple sequence and with sequence
change at runtime or through external trigger.

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

// ----- operating system dependent functions for thread, event, keyboard and mutex handling -----
#include "../common/ostools/spcm_oswrap.h"
#include "../common/ostools/spcm_ostools.h"

// ----- standard c include files -----
#include <stdio.h>
#include <stdlib.h>
#include <math.h>



#define USING_EXTERNAL_TRIGGER  0   // wait for keystroke to switch to next sequence group
//#define USING_EXTERNAL_TRIGGER  1   // use external trigger to switch to next sequence group



/*
**************************************************************************
vDoCardSetup
**************************************************************************
*/

void vDoCardSetup (ST_SPCM_CARDINFO *pstCard)
    {
    int i;
    int64 llChannelMask;


    // for easy sample data handling set mask for one channel (8 ... 16 bit per sample)
    if ((pstCard->eCardFunction == DigitalOut)
     || (pstCard->eCardFunction == DigitalIO))
        llChannelMask = 0xff;
    else
        llChannelMask = CHANNEL0;


    // we try to set the samplerate to a quarter of maximum on internal PLL, no clock output
    bSpcMSetupClockPLL (pstCard, pstCard->llMaxSamplerate / 4, false);
    printf ("Sampling rate set to %.1lf MHz\n", (double) pstCard->llSetSamplerate / 1000000);


    // setup replay sequence mode with 32 segments
    bSpcMSetupModeRepSequence (pstCard, llChannelMask, 32);

    if (!USING_EXTERNAL_TRIGGER)
        // software trigger 
        bSpcMSetupTrigSoftware (pstCard, false);
    else
        // external TTL trigger (with termination off) if the "SPCSEQ_ENDLOOPONTRIG" flag is used
        // The start of the first step need then a trigger too or a "M2CMD_CARD_FORCETRIGGER" command.
        bSpcMSetupTrigExternal (pstCard, SPC_TM_POS, false);


    // type dependent card setup
    switch (pstCard->eCardFunction)
        {

        // analog generator card setup
        case AnalogOut:

            // program all output channels to +/- 1 V with no offset and hold last sample at end of output
            for (i=0; i < pstCard->lMaxChannels; i++)
                bSpcMSetupAnalogOutputChannel (pstCard, i, 1000, 0, 0, SPCM_STOPLVL_HOLDLAST);
            break;

        // digital generator card setup
        case DigitalOut:
        case DigitalIO:
            for (i=0; i < pstCard->uCfg.stDIO.lGroups; i++)
                bSpcMSetupDigitalOutput (pstCard, i, SPCM_STOPLVL_LOW, 0, 3300);
            break;
        }
    }



/*
**************************************************************************
vWriteSegmentData
**************************************************************************
*/

void vWriteSegmentData (ST_SPCM_CARDINFO *pstCard, uint32 dwSegmentIndex, uint32 dwSegmentLenSample, void* pvSegData)
    {
    uint32 dwError =      0;
    uint32 dwSegLenByte = dwSegmentLenSample * pstCard->lBytesPerSample;

    // for 8 bit sample rearange order in buffer from 16 bit to 8 bit
    if (pstCard->lBytesPerSample == 1)
        {
        int8*  pcData = (int8*)  pvSegData;
        int16* pnData = (int16*) pvSegData;

        for (uint32 i = 0; i < dwSegmentLenSample; i++)
            pcData[i] = pnData[i] & 0xff;
        }

    // setup
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCard->hDrv, SPC_SEQMODE_WRITESEGMENT, dwSegmentIndex);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCard->hDrv, SPC_SEQMODE_SEGMENTSIZE,  dwSegmentLenSample);

    // write data to board (main) sample memory
    if (!dwError) dwError = spcm_dwDefTransfer_i64 (pstCard->hDrv, SPCM_BUF_DATA, SPCM_DIR_PCTOCARD, 0, pvSegData, 0, dwSegLenByte);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCard->hDrv, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);
    }



/*
**************************************************************************
DoDataCalculation: calculates and writes the output data for all segments
**************************************************************************
*/

bool bDoDataCalculation (ST_SPCM_CARDINFO *pstCard)
    {
    uint32 i;
    uint32 dwSegmentLenSample, dwSegLenByte;

    printf ("Calculation of output data\n");


    uint32 dwFactor = 1;
    // This series has a slightly increased minimum size value.
    if (pstCard->bM4i)
        dwFactor = 6;

    // buffer for data transfer (allocate for 8 bit sample 2 byte too)
    dwSegLenByte = 2 * dwFactor * 512; // max value are from the sine calculation
    void* pvBuffer = (void*) pvAllocMemPageAligned (dwSegLenByte);
    if (!pvBuffer)
        return nSpcMErrorMessageStdOut (pstCard, "Memory allocation error\n", false) != -1;
    int16 *pnData = (int16*) pvBuffer;


    // helper values: Full Scale
    uint32 dwFS = 127;
    if (pstCard->eCardFunction == AnalogOut)
        dwFS = pstCard->uCfg.stAO.lMaxDACValue;
    uint32 dwFShalf = dwFS / 2;


    // (main) sample memory segment index:
#define SEG_RAMPUP      0 // ramp up
#define SEG_RAMPDOWN    1 // ramp down
#define SEG_SYNC        2 // negative sync puls, for example oscilloscope trigger
//                      3 // unused
#define SEG_Q1SIN       4 // first quadrant of sinus signal
#define SEG_Q2SIN       5 // second quadrant of sinus signal
#define SEG_Q3SIN       6 // third quadrant of sinus signal
#define SEG_Q4SIN       7 // fourth quadrant of sinus signal
#define SEG_STOP        8 // DC level for stop/end
//                      remainder: unused


    // --- sync puls: first half zero, second half -FS
    dwSegmentLenSample = dwFactor * 80;
    for (i = 0; i < dwSegmentLenSample / 2; i++)
        pnData[i] = 0;

    if (pstCard->eCardFunction == AnalogOut)
        for (; i < dwSegmentLenSample; i++)
            pnData[i] = -((int16) dwFS);
    else  // digital boards (no two's complement)
        for (; i < dwSegmentLenSample; i++)
            pnData[i] = 0xff;

    vWriteSegmentData (pstCard, SEG_SYNC, dwSegmentLenSample, pvBuffer);


    // --- ramp up
    dwSegmentLenSample = dwFactor * 64;
    for (i = 0; i < dwSegmentLenSample; i++)
        pnData[i] = (int16) (i * dwFShalf / dwSegmentLenSample);

    vWriteSegmentData (pstCard, SEG_RAMPUP, dwSegmentLenSample, pvBuffer);


    // --- ramp down
    dwSegmentLenSample = dwFactor * 64;
    for (i = 0; i < dwSegmentLenSample; i++)
        pnData[i] = (int16) (dwFS - (i * dwFShalf / dwSegmentLenSample));

    vWriteSegmentData (pstCard, SEG_RAMPDOWN, dwSegmentLenSample, pvBuffer);


    // --- sinus
    dwSegmentLenSample = dwFactor * 512;
    for (i = 0; i < dwSegmentLenSample; i++)
        pnData[i] = (int16) (dwFShalf + (dwFShalf * sin (2.0 * 3.14159 * i / dwSegmentLenSample) + 0.5));

    // write each quadrant in a own segment
    vWriteSegmentData (pstCard, SEG_Q1SIN, dwSegmentLenSample / 4, (void*) &pnData[0 * dwSegmentLenSample / 4]);
    vWriteSegmentData (pstCard, SEG_Q2SIN, dwSegmentLenSample / 4, (void*) &pnData[1 * dwSegmentLenSample / 4]);
    vWriteSegmentData (pstCard, SEG_Q3SIN, dwSegmentLenSample / 4, (void*) &pnData[2 * dwSegmentLenSample / 4]);
    vWriteSegmentData (pstCard, SEG_Q4SIN, dwSegmentLenSample / 4, (void*) &pnData[3 * dwSegmentLenSample / 4]);


    // --- DC level
    dwSegmentLenSample = dwFactor * 128;
    for (i = 0; i < dwSegmentLenSample; i++)
        pnData[i] = (int16) (dwFS / 2);

    vWriteSegmentData (pstCard, SEG_STOP, dwSegmentLenSample, pvBuffer);


    vFreeMemPageAligned (pvBuffer, dwSegLenByte);

    return true;
    }



/*
**************************************************************************
vWriteStepEntry
**************************************************************************
*/

void vWriteStepEntry (ST_SPCM_CARDINFO *pstCard, uint32 dwStepIndex,
                      uint32 dwStepNextIndex, uint32 dwSegmentIndex, uint32 dwLoops, uint32 dwFlags)
    {
    uint32 dwError =         0;
    uint64 qwSequenceEntry = 0;

    // setup register value
    qwSequenceEntry = (dwFlags & ~SPCSEQ_LOOPMASK) | (dwLoops & SPCSEQ_LOOPMASK);
    qwSequenceEntry <<= 32;
    qwSequenceEntry |= ((dwStepNextIndex << 16)& SPCSEQ_NEXTSTEPMASK) | (dwSegmentIndex & SPCSEQ_SEGMENTMASK);

    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCard->hDrv, SPC_SEQMODE_STEPMEM0 + dwStepIndex, qwSequenceEntry);
    }



/*
**************************************************************************
bool bConfigureSequence
**************************************************************************
*/

bool bConfigureSequence (ST_SPCM_CARDINFO *pstCard)
    {
    // sequence memory
    // four sequence loops are programmed (each with 6 steps)
    // a keystroke or ext. trigger switched to the next sequence
    // the loop value for the ramp increase in each sequence
    //  0 ...  5: sync, Q1sin, Q2sin, Q3sin, Q4sin, ramp up
    //  8 ... 13: sync, Q2sin, Q3sin, Q4sin, Q1sin, ramp down
    // 16 ... 21: sync, Q3sin, Q4sin, Q1sin, Q2sin, ramp up
    // 24 ... 29: sync, Q4sin, Q1sin, Q2sin, Q3sin, ramp down

                            // +-- StepIndex
                            // |   +-- StepNextIndex
                            // |   |   +-- SegmentIndex
                            // |   |   |              +-- Loops
                            // |   |   |              |   +-- Flags: SPCSEQ_ENDLOOPONTRIG
    // sin                  // |   |   |              |   |          For using this flag disable Software-Trigger above.
    vWriteStepEntry (pstCard,  0,  1,  SEG_SYNC,      3,  0);
    vWriteStepEntry (pstCard,  1,  2,  SEG_Q1SIN,     1,  0);
    vWriteStepEntry (pstCard,  2,  3,  SEG_Q2SIN,     1,  0);
    vWriteStepEntry (pstCard,  3,  4,  SEG_Q3SIN,     1,  0);
    vWriteStepEntry (pstCard,  4,  5,  SEG_Q4SIN,     1,  0);
    if (!USING_EXTERNAL_TRIGGER)
    vWriteStepEntry (pstCard,  5,  1,  SEG_RAMPDOWN,  1,  0);
    else
    vWriteStepEntry (pstCard,  5,  8,  SEG_RAMPDOWN,  1,  SPCSEQ_ENDLOOPONTRIG);
#define LAST_STEP_OFFSET    5

    // cos
    vWriteStepEntry (pstCard,  8,  9,  SEG_SYNC,      3,  0);
    vWriteStepEntry (pstCard,  9, 10,  SEG_Q2SIN,     1,  0);
    vWriteStepEntry (pstCard, 10, 11,  SEG_Q3SIN,     1,  0);
    vWriteStepEntry (pstCard, 11, 12,  SEG_Q4SIN,     1,  0);
    vWriteStepEntry (pstCard, 12, 13,  SEG_Q1SIN,     1,  0);
    if (!USING_EXTERNAL_TRIGGER)
    vWriteStepEntry (pstCard, 13,  9,  SEG_RAMPUP,    2,  0);
    else
    vWriteStepEntry (pstCard, 13, 16,  SEG_RAMPUP,    2,  SPCSEQ_ENDLOOPONTRIG);

    // inverted sin
    vWriteStepEntry (pstCard, 16, 17,  SEG_SYNC,      3,  0);
    vWriteStepEntry (pstCard, 17, 18,  SEG_Q3SIN,     1,  0);
    vWriteStepEntry (pstCard, 18, 19,  SEG_Q4SIN,     1,  0);
    vWriteStepEntry (pstCard, 19, 20,  SEG_Q1SIN,     1,  0);
    vWriteStepEntry (pstCard, 20, 21,  SEG_Q2SIN,     1,  0);
    if (!USING_EXTERNAL_TRIGGER)
    vWriteStepEntry (pstCard, 21, 17,  SEG_RAMPDOWN,  3,  0);
    else
    vWriteStepEntry (pstCard, 21, 24,  SEG_RAMPDOWN,  3,  SPCSEQ_ENDLOOPONTRIG);

    // inverted cos
    vWriteStepEntry (pstCard, 24, 25,  SEG_SYNC,      3,  0);
    vWriteStepEntry (pstCard, 25, 26,  SEG_Q4SIN,     1,  0);
    vWriteStepEntry (pstCard, 26, 27,  SEG_Q1SIN,     1,  0);
    vWriteStepEntry (pstCard, 27, 28,  SEG_Q2SIN,     1,  0);
    vWriteStepEntry (pstCard, 28, 29,  SEG_Q3SIN,     1,  0);
    vWriteStepEntry (pstCard, 29, 30,  SEG_RAMPUP,    4,  0);
    vWriteStepEntry (pstCard, 30, 30,  SEG_STOP,      1,  SPCSEQ_END);  // M2i: only a few sample from this segment are replayed
                                                                        // M4i: the complete segment is replayed

    // Configure the beginning (index of first seq-entry to start) of the sequence replay.
    spcm_dwSetParam_i32 (pstCard->hDrv, SPC_SEQMODE_STARTSTEP, 0);

    // dump steps if necessary
    if (0)
        {
        printf ("\n");
        for (int i = 0; i < 32; i++)
            {
            int64 llTemp;
            spcm_dwGetParam_i64 (pstCard->hDrv, SPC_SEQMODE_STEPMEM0 + i, &llTemp);
            printf ("Step %.2d: 0x%08x_%08x\n", i, (uint32) (llTemp >> 32), (uint32) llTemp);
            }
        printf ("\n\n");
        }

    return true;
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


    // ------------------------------------------------------------------------
    // init card number 0 (the first card in the system), get some information and print it
    // uncomment the second line and replace the IP address to use remote
    // cards like in a generatorNETBOX
    if (bSpcMInitCardByIdx (&stCard, 0))
    //if (bSpcMInitCardByIdx (&stCard, "192.168.1.10", 0))
        {
        printf (pszSpcMPrintDocumentationLink (&stCard, szBuffer, sizeof (szBuffer)));
        printf (pszSpcMPrintCardInfo (&stCard, szBuffer, sizeof (szBuffer)));
        }
    else
        return nSpcMErrorMessageStdOut (&stCard, "Error: Could not open card\n", true);


    // check whether we support this card type in the example
    if ((stCard.eCardFunction != AnalogOut) && (stCard.eCardFunction != DigitalOut) && (stCard.eCardFunction != DigitalIO))
        return nSpcMErrorMessageStdOut (&stCard, "Error: Card function not supported by this example\n", false);
 
    if (!(stCard.lFeatureMap & SPCM_FEAT_SEQUENCE)
      || (stCard.bM2i && (stCard.lCtrlFwVersion < 20))
      || (stCard.bM4i && (stCard.lCtrlFwVersion < 14))) // on M2i and M4i sequence mode has been released as update. on M2p it is available since first release
        return nSpcMErrorMessageStdOut (&stCard, "Error: option 'sequence replay' not installed or firmware version to old\n", false);

    printf ("\n");


    // ------------------------------------------------------------------------
    // do the card setup, error is routed in the structure so we don't care for the return values
    if (!stCard.bSetError)
        {
        vDoCardSetup (&stCard);
        }


    // ------------------------------------------------------------------------
    // calculate the amount of data we need and allocate memory buffer
    if (!stCard.bSetError)
        {
        // calculate the data
        if (!bDoDataCalculation (&stCard))
            return nSpcMErrorMessageStdOut (&stCard, "Data calculation failed\n", false);

        printf ("... data has been transferred to board memory\n");

        // setup the the sequence
        if (!bConfigureSequence (&stCard))
            return nSpcMErrorMessageStdOut (&stCard, "Sequence setup failed\n", false);

        printf ("... sequence configured\n");
        }


    // ------------------------------------------------------------------------
    // start the generation
    if (!stCard.bSetError)
        {
        // We'll start and wait until all sequences are replayed.
        spcm_dwSetParam_i32 (stCard.hDrv, SPC_TIMEOUT, 0);
        printf ("\nStarting the card\n");
        if (spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER) != ERR_OK)
            {
            spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_STOP);
            return nSpcMErrorMessageStdOut (&stCard, "... error", false);
            }


        printf ("\nsequence replay runs, switch to next sequence (3 times possible) with");
        if (!USING_EXTERNAL_TRIGGER)
            printf ("\n key: c ... change sequence\n\n");
        else
            printf ("\n a (slow) TTL signal on external trigger input connector\n\n");

        int32 lCardStatus =       0;
        uint32 dwSequenceActual = 0;    // first step in a sequence
        uint32 dwSequenceNext;
        do
            {
            if (bKbhit ())
                {
                char c = cGetch ();
                switch (c)
                    {
                    case 27: // ESC
                        spcm_dwSetParam_i32 (stCard.hDrv, SPC_M2CMD, M2CMD_CARD_STOP);
                        break;

                    case 'c':
                    case 'C':
                        if (!USING_EXTERNAL_TRIGGER)
                            {
                            dwSequenceNext = (dwSequenceActual + 8) % 32;
                            printf ("sequence %d\n", dwSequenceNext / 8);

                            // switch to next sequence
                            // (before it is possible to overwrite the segment data of the new used segments with new values)
                            uint32 dwError;
                            int64 llStep;

                            // --- change the next step value from the sequence end entry in the actual sequence
                            dwError = spcm_dwGetParam_i64 (stCard.hDrv, SPC_SEQMODE_STEPMEM0 + dwSequenceActual + LAST_STEP_OFFSET, &llStep);
                            llStep = (llStep & ~((int64)SPCSEQ_NEXTSTEPMASK)) | dwSequenceNext << 16;
                            dwError = spcm_dwSetParam_i64 (stCard.hDrv, SPC_SEQMODE_STEPMEM0 + dwSequenceActual + LAST_STEP_OFFSET, llStep);

                            dwSequenceActual = dwSequenceNext;
                            }
                        break;
                    }
                }
            else
                {
                SPCM_NAMESPACE::spcm_vSuspendThread (10); // ms

                // Demonstrate the two different sequence status values at M2i and M4i / M2p cards.
                static int32 s_lSeqStatusOld = 0;
                int32 lSeqStatus;
                spcm_dwGetParam_i32 (stCard.hDrv, SPC_SEQMODE_STATUS, &lSeqStatus);

                // Avoid a lot of outputs in none external trigger mode.
                if (USING_EXTERNAL_TRIGGER)
                    {
                    if (s_lSeqStatusOld != lSeqStatus)
                        {
                        s_lSeqStatusOld = lSeqStatus;

                        if (stCard.bM2i)
                            {
                            if (lSeqStatus & SEQSTAT_STEPCHANGE)
                                printf ("status: sequence changed\n");
                            }
                        if (stCard.bM4i
                         || stCard.bM2p)
                            {
                            // Valid values only at a startet card available.
                            if (lCardStatus & M2STAT_CARD_PRETRIGGER)
                                printf ("status: actual sequence number: %d\n", lSeqStatus);
                            }
                        }
                    }
                }

            spcm_dwGetParam_i32 (stCard.hDrv, SPC_M2STATUS, &lCardStatus);
            }
        while (!(lCardStatus & M2STAT_CARD_READY));

        printf ("\n\n programm finished (press key)");
        cGetch ();
        }


    // ------------------------------------------------------------------------
    // print error information if an error occured
    if (stCard.bSetError)
        return nSpcMErrorMessageStdOut (&stCard, "An error occured while programming the card:\n", true);


    // close the driver
    vSpcMCloseCard (&stCard);

    return EXIT_SUCCESS;
    }
