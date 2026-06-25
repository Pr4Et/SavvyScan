/*
**************************************************************************

simple_rec_segment_statistic.cpp                         (c) Spectrum GmbH

**************************************************************************

Example for all M4i analog acquisition cards. 
Shows a simple segment statistic example
using only a few necessary commands.
  
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


// define structures for more easy data access
typedef struct
    {
    SPCM_SEGSTAT_STRUCT_CHx pst_Channel[4];
    } SPCM_SEGSTAT_STRUCT_4CH;

typedef struct
    {
    SPCM_SEGSTAT_STRUCT_CHx pst_Channel[2];
    } SPCM_SEGSTAT_STRUCT_2CH;

typedef struct
    {
    SPCM_SEGSTAT_STRUCT_CHx pst_Channel[1];
    } SPCM_SEGSTAT_STRUCT_1CH;


/*
**************************************************************************
main 
**************************************************************************
*/

int main ()
    {
    drv_handle  hCard;
    int32       lCardType, lSerialNumber, lFncType, lFeatMask;
    int32       lMaxChannels, lChannelMask;
    int64       llMaxSamplerate;
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
    
    // check for necessary segment statistic feature
    spcm_dwGetParam_i32 (hCard, SPC_PCIEXTFEATURES, &lFeatMask);
    if ((lFeatMask & SPCM_FEAT_EXTFW_SEGSTAT) == 0)
        {
        printf ("This example requires segment statistic feature installed \n");
        spcm_vClose (hCard);
        return 0;
        }      


    // do a simple standard setup 
    uint32 dwSegmentsize =   KILO_B(4); // define length of each segment           
    uint32 dwPosttrigger =   KILO_B(2); // define samples as posttrigger
    uint32 dwNumOfSegments = 16;        // number of segments to record

    spcm_dwSetParam_i32 (hCard, SPC_SEGMENTSIZE,   dwSegmentsize);
    spcm_dwSetParam_i32 (hCard, SPC_POSTTRIGGER,   dwPosttrigger);

    // (all data resides in on-board memory, no streaming)
    spcm_dwSetParam_i32 (hCard, SPC_MEMSIZE,       dwSegmentsize * dwNumOfSegments);    
    
    // enable all available channels
    spcm_dwGetParam_i32 (hCard, SPC_MIINST_CHPERMODULE, &lMaxChannels);
    lChannelMask = ((int32) 1 << lMaxChannels) - 1;
    spcm_dwSetParam_i32 (hCard, SPC_CHENABLE, lChannelMask);

    // use segment statistic mode, timeout set to 5s
    spcm_dwSetParam_i32 (hCard, SPC_CARDMODE,     SPC_REC_STD_SEGSTATS);
    spcm_dwSetParam_i32 (hCard, SPC_TIMEOUT,      5000);
    
    // disable external triggers
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_ORMASK,  0); 
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_ANDMASK, 0);                     
    
    // enable channel trigger
    int32 lTriggerChannel = 0;
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_CH_ORMASK0, (SPC_TMASK0_CH0 << lTriggerChannel));        
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_CH_ANDMASK0, 0);                     
    spcm_dwSetParam_i32 (hCard, SPC_TRIG_CH0_MODE + lTriggerChannel, SPC_TM_NEG);    

    // set up maximum available sample rate, no clock output
    spcm_dwGetParam_i64 (hCard, SPC_MIINST_MAXADCLOCK, &llMaxSamplerate);
    spcm_dwSetParam_i32 (hCard, SPC_CLOCKMODE,  SPC_CM_INTPLL);
    spcm_dwSetParam_i64 (hCard, SPC_SAMPLERATE, llMaxSamplerate);
    spcm_dwSetParam_i32 (hCard, SPC_CLOCKOUT,   0);             
    
    // set timestamp to start reset 
    spcm_dwSetParam_i32 (hCard, SPC_TIMESTAMP_CMD, SPC_TSMODE_STARTRESET);



    // settings for the buffer handling
    int64 llBufferSize = sizeof (SPCM_SEGSTAT_STRUCT_CHx) * lMaxChannels * dwNumOfSegments;

    // define the data buffer
    void* pvData = (void*) pvAllocMemPageAligned ((uint64) llBufferSize);
    
    // check for memory allocation errors
    if (!pvData)
        {
        printf ("memory allocation failed\n");
        spcm_vClose (hCard);
        return 0;
        }   



    // start the acquisition
    dwError = spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_CARD_START | M2CMD_CARD_ENABLETRIGGER | M2CMD_CARD_WAITREADY);

    // check for error
    if (dwError != ERR_OK)
        {
        spcm_dwGetErrorInfo_i32 (hCard, NULL, NULL, szErrorTextBuffer);
        printf ("%s\n", szErrorTextBuffer);
        vFreeMemPageAligned (pvData, (uint64) llBufferSize);
        spcm_vClose (hCard);
        return 0;
        }
    


    // set up and start data transfer
    spcm_dwDefTransfer_i64 (hCard, SPCM_BUF_DATA, SPCM_DIR_CARDTOPC, 0, pvData , 0, llBufferSize);
    spcm_dwSetParam_i32 (hCard, SPC_M2CMD, M2CMD_DATA_STARTDMA | M2CMD_DATA_WAITDMA);
    

    // cast and printout data (depending on active channels to structure to get more easy access)
    if (lMaxChannels == 1)
        {
        SPCM_SEGSTAT_STRUCT_1CH* ppstData = (SPCM_SEGSTAT_STRUCT_1CH*) pvData;

        // read out data of every segment (one channel only to keep things simple)
        for (uint32 dwSegment = 0; dwSegment < dwNumOfSegments; dwSegment++)
            {
            printf("\nSegment %.2d: Min: %7.d  Max:%7.d, Avrg; %f",
                dwSegment,
                ppstData[dwSegment].pst_Channel[lTriggerChannel].nMin,
                ppstData[dwSegment].pst_Channel[lTriggerChannel].nMax,
                ((double) (ppstData[dwSegment].pst_Channel[lTriggerChannel].llAvrg) / (double) dwSegmentsize));
            }
        }
    else if (lMaxChannels == 2)
        {
        SPCM_SEGSTAT_STRUCT_2CH* ppstData = (SPCM_SEGSTAT_STRUCT_2CH*) pvData;

        // read out data of every segment (one channel only to keep things simple)
        for (uint32 dwSegment = 0; dwSegment < dwNumOfSegments; dwSegment++)
            {
            printf("\nSegment %.2d: Min: %7.4d  Max:%7.4d, Avrg; %f",
                dwSegment,
                ppstData[dwSegment].pst_Channel[lTriggerChannel].nMin,
                ppstData[dwSegment].pst_Channel[lTriggerChannel].nMax,
                ((double) (ppstData[dwSegment].pst_Channel[lTriggerChannel].llAvrg) / (double) dwSegmentsize));
            }
        }
    else
        {
        SPCM_SEGSTAT_STRUCT_4CH* ppstData = (SPCM_SEGSTAT_STRUCT_4CH*) pvData;

        // read out data of every segment (one channel only to keep things simple)
        for (uint32 dwSegment = 0; dwSegment < dwNumOfSegments; dwSegment++)
            {
            printf("\nSegment %.2d: Min: %7.4d  Max:%7.4d, Avrg; %f",
                dwSegment,
                ppstData[dwSegment].pst_Channel[lTriggerChannel].nMin,
                ppstData[dwSegment].pst_Channel[lTriggerChannel].nMax,
                ((double) (ppstData[dwSegment].pst_Channel[lTriggerChannel].llAvrg) / (double) dwSegmentsize));
            }
        }   
    

    // clean up
    printf ("\n\nFinished...\n");    

    vFreeMemPageAligned (pvData, (uint64) llBufferSize);
    spcm_vClose (hCard);

    return EXIT_SUCCESS;
    }

