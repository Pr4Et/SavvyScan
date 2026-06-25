/*
**************************************************************************

spcm_lib_card.cpp                                        (c) Spectrum GmbH

**************************************************************************

Supplies different common functions for C/C++ programs accessing the 
SpcM driver interface. Feel free to use this source for own projects and
modify it in any kind

Documentation for the API as well as a detailed description of the hardware
can be found in the manual for each device which can be found on our website:
https://www.spectrum-instrumentation.com/en/downloads

Further information can be found online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/knowledge-base-overview

**************************************************************************
*/



// ----- include of common example librarys -----
#include "../common/spcm_lib_card.h"

// ----- standard c include files -----
#include <cstdio>
#include <cstring>
#include <cstdlib>



/*
**************************************************************************
**************************************************************************
**************************************************************************
  Initialisation and error handling
**************************************************************************
**************************************************************************
**************************************************************************
*/



/*
**************************************************************************
bSpcMInitCardByIdx:

opens the driver with the given indes, reads out card information and
fills the CARDINFO structure
**************************************************************************
*/

bool bSpcMInitCardByName (ST_SPCM_CARDINFO *pstCardInfo, char* szDrvName);

bool bSpcMInitCardByIdx (ST_SPCM_CARDINFO *pstCardInfo, int32 lCardIdx)
    {
    if (!pstCardInfo)
        return false;

    // open the driver for card. We can use the linux notation here as the windows driver
    // only looks for the ending number. Change this line if the linux drivers are named
    // different than default
    char szDrvName[20];
    sprintf (szDrvName, "/dev/spcm%d", lCardIdx);
    return bSpcMInitCardByName (pstCardInfo, szDrvName);
    }

bool bSpcMInitCardByIdx (ST_SPCM_CARDINFO *pstCardInfo, const char* szIP, int32 lCardIdx)
    {
    if (!pstCardInfo)
        return false;

    char szVISA[50];
    sprintf (szVISA, "TCPIP::%s::inst%d::INSTR", szIP, lCardIdx);
    return bSpcMInitCardByName (pstCardInfo, szVISA);
    }

bool bSpcMInitCardByName (ST_SPCM_CARDINFO *pstCardInfo, char* szDrvName)
    {
    int32 lTmp;

    // clear the card info to have defined values
    memset ((void*) pstCardInfo, 0, sizeof(ST_SPCM_CARDINFO));
    pstCardInfo->lSetChannels =     1;
    pstCardInfo->llSetSamplerate =  1;


    pstCardInfo->hDrv = spcm_hOpen (szDrvName);
    if (!pstCardInfo->hDrv)
        {
        pstCardInfo->lErrorCode = spcm_dwGetErrorInfo_i32 (pstCardInfo->hDrv, NULL, NULL, pstCardInfo->szError);

        // card might be just "in use", and we can display some more info then
        spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_PCITYP,      &pstCardInfo->lCardType);
        spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_PCISERIALNO, &pstCardInfo->lSerialNumber);

        return false;
        }

    // ----- get index of card from name -----
    if (strncmp (szDrvName, "TCPIP", 5) == 0)
        {
        // if VISA-String contains instX, we extract the number
        // otherwise we default to zero
        char* szInst = strstr (szDrvName, "inst");
        if (szInst != NULL)
            pstCardInfo->lCardIdx = atoi (strpbrk (szInst, "0123456789"));
        else
            pstCardInfo->lCardIdx = 0;

        pstCardInfo->bRemote = true;
        }
    else
        {
        // name should be /dev/spcmX or just a number, so we locate first number in string
        // and convert it to integer
        pstCardInfo->lCardIdx = atoi (strpbrk (szDrvName, "0123456789"));
        }

    // read out card information and store it in the card info structure
    spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_PCITYP,                 &pstCardInfo->lCardType);
    spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_PCISERIALNO,            &pstCardInfo->lSerialNumber);
    spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_PCIFEATURES,            &pstCardInfo->lFeatureMap);
    spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_PCIEXTFEATURES,         &pstCardInfo->lExtFeatureMap);    
    spcm_dwGetParam_i64 (pstCardInfo->hDrv, SPC_PCIMEMSIZE,             &pstCardInfo->llInstMemBytes);
    spcm_dwGetParam_i64 (pstCardInfo->hDrv, SPC_MIINST_MINADCLOCK,      &pstCardInfo->llMinSamplerate);
    spcm_dwGetParam_i64 (pstCardInfo->hDrv, SPC_MIINST_MAXADCLOCK,      &pstCardInfo->llMaxSamplerate);
    spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_MIINST_MODULES,         &pstCardInfo->lModulesCount);
    spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_MIINST_CHPERMODULE,     &pstCardInfo->lMaxChannels);
    spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_MIINST_BYTESPERSAMPLE,  &pstCardInfo->lBytesPerSample);
    spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_GETDRVVERSION,          &pstCardInfo->lLibVersion);
    spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_GETKERNELVERSION,       &pstCardInfo->lKernelVersion);


    // fill in the examples flags
    pstCardInfo->bM2i = false;
    pstCardInfo->bM3i = false;
    pstCardInfo->bM4i = false;
    pstCardInfo->bM2p = false;
    switch (pstCardInfo->lCardType & TYP_SERIESMASK)
        {
        case TYP_M2ISERIES:     pstCardInfo->bM2i = true; break;
        case TYP_M2IEXPSERIES:  pstCardInfo->bM2i = true; break;
        case TYP_M3ISERIES:     pstCardInfo->bM3i = true; break;
        case TYP_M3IEXPSERIES:  pstCardInfo->bM3i = true; break;
        case TYP_M4IEXPSERIES:  pstCardInfo->bM4i = true; break;
        case TYP_M4XEXPSERIES:  pstCardInfo->bM4i = true; break;
        case TYP_M2PEXPSERIES:  pstCardInfo->bM2p = true; break;
        default:                break;
        }

    spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_PCIVERSION,             &lTmp);
    pstCardInfo->lBaseHwVersion = (lTmp >> 16) & 0xffff;
    pstCardInfo->lCtrlFwVersion = lTmp  & 0xffff;

    spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_PCIMODULEVERSION,       &lTmp);
    pstCardInfo->lModHwVersion = (lTmp >> 16) & 0xffff;
    pstCardInfo->lModFwVersion = lTmp  & 0xffff;

    // we need to recalculate the channels value as the driver returns channels per module
    pstCardInfo->lMaxChannels *= pstCardInfo->lModulesCount;

    // examin the type of driver
    int32 lFncType;
    spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_FNCTYPE,                &lFncType);

    switch (lFncType)
        {
        case SPCM_TYPE_AI:  pstCardInfo->eCardFunction = AnalogIn; break;
        case SPCM_TYPE_AO:  pstCardInfo->eCardFunction = AnalogOut; break;
        case SPCM_TYPE_DI:  pstCardInfo->eCardFunction = DigitalIn; break;
        case SPCM_TYPE_DO:  pstCardInfo->eCardFunction = DigitalOut; break;
        case SPCM_TYPE_DIO: pstCardInfo->eCardFunction = DigitalIO; break;
        }

    // loading the function dependant part of the CardInfo structure
    switch (pstCardInfo->eCardFunction)
        {
        case AnalogIn:
            {
            int i;
            int32 lAIFeatures;

            spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_MIINST_BITSPERSAMPLE,  &pstCardInfo->uCfg.stAI.lResolution);
            spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_READAIPATHCOUNT,  &pstCardInfo->uCfg.stAI.lPathCount);
            for (int32 lPath = 0; lPath < pstCardInfo->uCfg.stAI.lPathCount; ++lPath)
                {
                spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_READAIPATH, lPath);

                spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_READIRCOUNT, &pstCardInfo->uCfg.stAI.astPath[lPath].lRangeCount);
                for (i=0; (i<pstCardInfo->uCfg.stAI.astPath[lPath].lRangeCount) && (i<SPCM_MAX_AIRANGE); i++)
                    {
                    spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_READRANGEMIN0 + i, &pstCardInfo->uCfg.stAI.astPath[lPath].lRangeMin[i]);
                    spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_READRANGEMAX0 + i, &pstCardInfo->uCfg.stAI.astPath[lPath].lRangeMax[i]);
                    }
                spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_READAIFEATURES, &lAIFeatures);

                pstCardInfo->uCfg.stAI.astPath[lPath].bInputTermAvailable = (lAIFeatures & SPCM_AI_TERM) != 0;
                pstCardInfo->uCfg.stAI.astPath[lPath].bDiffModeAvailable =  (lAIFeatures & SPCM_AI_DIFF) != 0;
                pstCardInfo->uCfg.stAI.astPath[lPath].bACCouplingAvailable =(lAIFeatures & SPCM_AI_ACCOUPLING) != 0;
                pstCardInfo->uCfg.stAI.astPath[lPath].bBWLimitAvailable =   (lAIFeatures & SPCM_AI_LOWPASS) != 0;
                pstCardInfo->uCfg.stAI.astPath[lPath].bOffsPercentMode =    (lAIFeatures & SPCM_AI_OFFSPERCENT) != 0;
                }

            // SPC_MIINST_MAXADCVALUE added with driver version 1.34, otherwise we have to calc it from the resolution
            if (ERR_OK != spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_MIINST_MAXADCVALUE, &pstCardInfo->uCfg.stAI.lMaxADCValue))
                {
                spcm_dwGetErrorInfo_i32 (pstCardInfo->hDrv, NULL, NULL, NULL);
                switch (pstCardInfo->uCfg.stAI.lResolution)
                    {
                    case 8:  pstCardInfo->uCfg.stAI.lMaxADCValue = 128; break;
                    case 12: pstCardInfo->uCfg.stAI.lMaxADCValue = 2048; break;
                    case 14: pstCardInfo->uCfg.stAI.lMaxADCValue = 8192; break;
                    default:
                    case 16: pstCardInfo->uCfg.stAI.lMaxADCValue = 32768; break;
                    }
                }

            break;
            }

        case AnalogOut:
            {
            int32 lAOFeatures;

            spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_READAOFEATURES, &lAOFeatures);
            spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_MIINST_BITSPERSAMPLE,  &pstCardInfo->uCfg.stAO.lResolution);

            pstCardInfo->uCfg.stAO.bGainProgrammable =       (lAOFeatures & SPCM_AO_PROGGAIN) != 0;
            pstCardInfo->uCfg.stAO.bOffsetProgrammable =     (lAOFeatures & SPCM_AO_PROGOFFSET) != 0;
            pstCardInfo->uCfg.stAO.bFilterAvailable =        (lAOFeatures & SPCM_AO_PROGFILTER) != 0;
            pstCardInfo->uCfg.stAO.bStopLevelProgrammable =  (lAOFeatures & SPCM_AO_PROGSTOPLEVEL) != 0;
            pstCardInfo->uCfg.stAO.bDiffModeAvailable =      (lAOFeatures & SPCM_AO_DIFF) != 0;

            // SPC_MIINST_MAXADCVALUE added with driver version 1.34, otherwise we have to calc it from the resolution
            if (ERR_OK != spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_MIINST_MAXADCVALUE, &pstCardInfo->uCfg.stAO.lMaxDACValue))
                {
                spcm_dwGetErrorInfo_i32 (pstCardInfo->hDrv, NULL, NULL, NULL);
                switch (pstCardInfo->uCfg.stAI.lResolution)
                    {
                    case 8:  pstCardInfo->uCfg.stAO.lMaxDACValue = 127; break;
                    case 12: pstCardInfo->uCfg.stAO.lMaxDACValue = 2047; break;
                    case 14: pstCardInfo->uCfg.stAO.lMaxDACValue = 8191; break;
                    default:
                    case 16: pstCardInfo->uCfg.stAO.lMaxDACValue = 32767; break;
                    }
                }
            else
                {
                // since driver version build 3738 is the value incremented
                if ((pstCardInfo->lLibVersion & 0xffff) >= 3738)
                    pstCardInfo->uCfg.stAO.lMaxDACValue--;
                }

            break;
            }

        case DigitalIn:
        case DigitalOut:
        case DigitalIO:
            {
            if ((pstCardInfo->eCardFunction == DigitalIn) || (pstCardInfo->eCardFunction == DigitalIO))
                {
                int32 lDIFeatures;
                spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_READDIFEATURES, &lDIFeatures);

                pstCardInfo->uCfg.stDIO.bInputTermAvailable = (lDIFeatures & SPCM_DI_TERM) != 0;
                pstCardInfo->uCfg.stDIO.bDiffModeAvailable =  (lDIFeatures & SPCM_DI_DIFF) != 0;
                }

            if ((pstCardInfo->eCardFunction == DigitalOut) || (pstCardInfo->eCardFunction == DigitalIO))
                {
                int32 lDOFeatures;
                spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_READDOFEATURES, &lDOFeatures);

                pstCardInfo->uCfg.stDIO.bDiffModeAvailable =      (lDOFeatures & SPCM_DO_DIFF) != 0;
                pstCardInfo->uCfg.stDIO.bStopLevelProgrammable =  (lDOFeatures & SPCM_DO_PROGSTOPLEVEL) != 0;
                pstCardInfo->uCfg.stDIO.bOutputLevelProgrammable = (lDOFeatures & SPCM_DO_PROGOUTLEVELS) != 0;
                }

            // grouping is the number of channels in one group, we recalculate this to the number of groups
            spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_READCHGROUPING, &pstCardInfo->uCfg.stDIO.lGroups);
            pstCardInfo->uCfg.stDIO.lGroups = pstCardInfo->lMaxChannels / pstCardInfo->uCfg.stDIO.lGroups;

            break;
            }

        default:
            break;
        }
    return true;
    }



/*
**************************************************************************
vSpcMCloseCard:

closes the driver
**************************************************************************
*/

void vSpcMCloseCard (ST_SPCM_CARDINFO *pstCardInfo)
    {
    if (!pstCardInfo)
        return;

    if (pstCardInfo->hDrv)
        {
        spcm_vClose (pstCardInfo->hDrv);
        pstCardInfo->hDrv = NULL;
        }

    }



/*
**************************************************************************
nSpcMErrorMessageStdOut:

prints the error message to std out and ends the driver if it's active
program can be left with this function
**************************************************************************
*/

int nSpcMErrorMessageStdOut (ST_SPCM_CARDINFO *pstCardInfo, const char* pszMessage, bool bPrintCardErr)
    {
    if (!pstCardInfo)
        return -2;

    printf (pszMessage);

    if (bPrintCardErr)
        printf (pstCardInfo->szError);

    if (pstCardInfo->hDrv)
        vSpcMCloseCard (pstCardInfo);
    pstCardInfo->hDrv = NULL;

    return -1;
    }



/*
**************************************************************************
pszSpcMTranslateRuntimeError: translates a runtime error code and prints
                              it to a given buffer
**************************************************************************
*/

char* pszSpcMTranslateRuntimeError (uint32 dwErrorCode, char* pszBuffer)
    {
    if (!pszBuffer)
        return NULL;

    switch (dwErrorCode)
        {
        case ERR_OK:                sprintf (pszBuffer, "No Error"); break;
        case ERR_ABORT:             sprintf (pszBuffer, "Abort of Wait Function by Stop Command"); break;
        case ERR_TIMEOUT:           sprintf (pszBuffer, "Timeout"); break;
        case ERR_FIFOBUFOVERRUN:    sprintf (pszBuffer, "FIFO SW Buffer Overrun"); break;
        case ERR_FIFOHWOVERRUN:     sprintf (pszBuffer, "FIFO HW Buffer Overrun"); break;
        case ERR_FIFOFINISHED:      sprintf (pszBuffer, "FIFO Mode finished"); break;
        default:                    sprintf (pszBuffer, "Unknown Error Code %d", dwErrorCode); break;
        }

    return pszBuffer;
    }



/*
**************************************************************************
pszSpcMPrintCardInfo: prints the card information to a string for display.
**************************************************************************
*/

void vStrCatWithLen (char* pszDest, char* pszSource, int32 lStrLen)
    {
    int   nPos = 0;

    while ((nPos++ < lStrLen) && (*pszDest))
        pszDest++;

    while (nPos < lStrLen)
        {
        *pszDest++ = *pszSource++;
        if (!(*pszSource))
            {
            (*pszDest) = 0;
            return;
            }
        }
    }

// *************************************************************************

char* pszSpcMPrintCardInfo (ST_SPCM_CARDINFO *pstCardInfo, char* pszBuffer, int32 lStrLen, bool bExtended)
    {
    char szTmp[100];

    if (!pstCardInfo || !pszBuffer)
        return NULL;

    memset (pszBuffer, 0, lStrLen);

    // the card type + serial number
    switch (pstCardInfo->lCardType & TYP_SERIESMASK)
        {
        case TYP_M2ISERIES:     sprintf (szTmp, "M2i.%04x sn %05d\n",     pstCardInfo->lCardType & TYP_VERSIONMASK, pstCardInfo->lSerialNumber); break;
        case TYP_M2IEXPSERIES:  sprintf (szTmp, "M2i.%04x-Exp sn %05d\n", pstCardInfo->lCardType & TYP_VERSIONMASK, pstCardInfo->lSerialNumber); break;
        case TYP_M3ISERIES:     sprintf (szTmp, "M3i.%04x sn %05d\n",     pstCardInfo->lCardType & TYP_VERSIONMASK, pstCardInfo->lSerialNumber); break;
        case TYP_M3IEXPSERIES:  sprintf (szTmp, "M3i.%04x-Exp sn %05d\n", pstCardInfo->lCardType & TYP_VERSIONMASK, pstCardInfo->lSerialNumber); break;
        case TYP_M4IEXPSERIES:  sprintf (szTmp, "M4i.%04x-x8 sn %05d\n",  pstCardInfo->lCardType & TYP_VERSIONMASK, pstCardInfo->lSerialNumber); break;
        case TYP_M4XEXPSERIES:  sprintf (szTmp, "M4x.%04x-x4 sn %05d\n",  pstCardInfo->lCardType & TYP_VERSIONMASK, pstCardInfo->lSerialNumber); break;
        case TYP_M2PEXPSERIES:  sprintf (szTmp, "M2p.%04x-x4 sn %05d\n",  pstCardInfo->lCardType & TYP_VERSIONMASK, pstCardInfo->lSerialNumber); break;
        default:                sprintf (szTmp, "Typ: %x not supported so far\n", pstCardInfo->lCardType); break;
        }
    vStrCatWithLen (pszBuffer, szTmp, lStrLen);

    // standard details of card
    if (bExtended)
        {
        sprintf (szTmp, "  Installed memory:  %lld MByte\n", pstCardInfo->llInstMemBytes / 1024 / 1024);
        vStrCatWithLen (pszBuffer, szTmp, lStrLen);
        sprintf (szTmp, "  Max sampling rate: %.1f MS/s\n", (double) pstCardInfo->llMaxSamplerate / 1000000);
        vStrCatWithLen (pszBuffer, szTmp, lStrLen);
        sprintf (szTmp, "  Channels:          %d\n", pstCardInfo->lMaxChannels);
        vStrCatWithLen (pszBuffer, szTmp, lStrLen);
        sprintf (szTmp, "  Kernel Version:    %d.%02d build %d\n", pstCardInfo->lKernelVersion >> 24, (pstCardInfo->lKernelVersion >> 16) & 0xff, pstCardInfo->lKernelVersion & 0xffff);
        vStrCatWithLen (pszBuffer, szTmp, lStrLen);
        sprintf (szTmp, "  Library Version    %d.%02d build %d\n", pstCardInfo->lLibVersion >> 24, (pstCardInfo->lLibVersion >> 16) & 0xff, pstCardInfo->lLibVersion & 0xffff);
        vStrCatWithLen (pszBuffer, szTmp, lStrLen);
        }

    return pszBuffer;
    }

/*
**************************************************************************
pszSpcMPrintDocumentationLink: builds a link to the download area for the specific card
**************************************************************************
*/

char* pszSpcMPrintDocumentationLink (const ST_SPCM_CARDINFO* pstCardInfo, char* pszBuffer, int32 lStrLen)
    {
    if (!pstCardInfo || !pszBuffer)
        return NULL;

    memset (pszBuffer, 0, lStrLen);

    const char* szSeries = NULL;
    if (pstCardInfo->bRemote)
        szSeries = "DN2";
    else
        {
        switch (pstCardInfo->lCardType & TYP_SERIESMASK)
            {
            case TYP_M2ISERIES:    szSeries = "M2i"; break;
            case TYP_M2IEXPSERIES: szSeries = "\"M2i Express\""; break;
            case TYP_M3ISERIES:    szSeries = "M3i"; break;
            case TYP_M3IEXPSERIES: szSeries = "\"M3i Express\""; break;
            case TYP_M4IEXPSERIES: szSeries = "M4i"; break;
            case TYP_M4XEXPSERIES: szSeries = "M4x"; break;
            case TYP_M2PEXPSERIES: szSeries = "M2p"; break;
            }
        }
    int lOffset = sprintf (pszBuffer, "A detailed description of the API as well as the hardware can be found in the manual:\n");
    lOffset += sprintf (pszBuffer + lOffset, "https://www.spectrum-instrumentation.com/en/downloads/drivers?Series=%s&Families=%xxx&Tab=Documents\n\n", szSeries, ((pstCardInfo->lCardType & TYP_FAMILYMASK) >> 8));

    lOffset += sprintf (pszBuffer + lOffset, "Further information can be found online in the Knowledge Base:\n");
    lOffset += sprintf (pszBuffer + lOffset, "https://www.spectrum-instrumentation.com/en/knowledge-base-overview\n\n");

    return pszBuffer;
    }

/*
**************************************************************************
bSpcMCheckSetError: checks for error code and reads out error information
**************************************************************************
*/

bool bSpcMCheckSetError (uint32 dwError, ST_SPCM_CARDINFO *pstCardInfo)
    {
    if (dwError)
        {
        pstCardInfo->bSetError = true;
        spcm_dwGetErrorInfo_i32 (pstCardInfo->hDrv, NULL, NULL, pstCardInfo->szError);
        return false;
        }
    return true;
    }




/*
**************************************************************************
**************************************************************************
**************************************************************************
  Mode setup
**************************************************************************
**************************************************************************
**************************************************************************
*/



/*
**************************************************************************
bSpcMSetupModeRecStdSingle: record standard mode single
**************************************************************************
*/

bool bSpcMSetupModeRecStdSingle (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llMemSamples, int64 llPostSamples)
    {   
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REC_STD_SINGLE);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_MEMSIZE,     llMemSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_POSTTRIGGER, llPostSamples);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     llMemSamples;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupModeRecFIFOSingle: record FIFO mode single run
**************************************************************************
*/

bool bSpcMSetupModeRecFIFOSingle (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llPreSamples, int64 llBlockToRec, int64 llLoopToRec)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;
    
    // check for invalid block/loop combinations
    if ((llBlockToRec && !llLoopToRec) || (!llBlockToRec && llLoopToRec))
        {
        sprintf (pstCardInfo->szError, "bSpcMSetupModeRecFIFOSingle: Loop and Blocks must be either both zero or both defined to non-zero\n");
        return false;
        }

    // segment size can't be zero, we adjust it here
    if (!llBlockToRec && !llLoopToRec)
        llBlockToRec = 1024;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REC_FIFO_SINGLE);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_PRETRIGGER,  llPreSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_SEGMENTSIZE, llBlockToRec);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_LOOPS,       llLoopToRec);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     0;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }


/*
**************************************************************************
bSpcMSetupModeRecStdAverage: record standard mode Average
**************************************************************************
*/

bool bSpcMSetupModeRecStdAverage (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llMemSamples, int64 llSegmentSize, int64 llPostSamples, int32 lAverages)
    {   
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REC_STD_AVERAGE);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_MEMSIZE,     llMemSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_SEGMENTSIZE, llSegmentSize);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_POSTTRIGGER, llPostSamples);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_AVERAGES,    lAverages);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     llMemSamples;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupModeRecStdMulti: record standard mode Multiple Recording
**************************************************************************
*/

bool bSpcMSetupModeRecStdMulti (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llMemSamples, int64 llSegmentSize, int64 llPostSamples)
    {   
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REC_STD_MULTI);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_MEMSIZE,     llMemSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_SEGMENTSIZE, llSegmentSize);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_POSTTRIGGER, llPostSamples);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     llMemSamples;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupModeRecStdABA: record standard mode ABA
**************************************************************************
*/

bool bSpcMSetupModeRecStdABA (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llMemSamples, int64 llSegmentSize, int64 llPostSamples, int32 lABADivider)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REC_STD_ABA);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_MEMSIZE,     llMemSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_SEGMENTSIZE, llSegmentSize);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_POSTTRIGGER, llPostSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_ABADIVIDER,  lABADivider);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     llMemSamples;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }


/*
**************************************************************************
bSpcMSetupModeRecFIFOAverage: record FIFO mode Average
**************************************************************************
*/

bool bSpcMSetupModeRecFIFOAverage (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llSegmentSize, int64 llPostSamples, int32 lAverages, int64 llSegmentsToRec)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;    

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REC_FIFO_AVERAGE);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_SEGMENTSIZE, llSegmentSize);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_POSTTRIGGER, llPostSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_LOOPS,       llSegmentsToRec);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_AVERAGES,    lAverages);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     0;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }


/*
**************************************************************************
bSpcMSetupModeRecFIFOMulti: record FIFO mode Multi
**************************************************************************
*/

bool bSpcMSetupModeRecFIFOMulti (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llSegmentSize, int64 llPostSamples, int64 llSegmentsToRec)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REC_FIFO_MULTI);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_SEGMENTSIZE, llSegmentSize);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_POSTTRIGGER, llPostSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_LOOPS,       llSegmentsToRec);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     0;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupModeRecFIFOABA: record FIFO mode ABA
**************************************************************************
*/

bool bSpcMSetupModeRecFIFOABA (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llSegmentSize, int64 llPostSamples, int32 lABADivider, int64 llSegmentsToRec)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REC_FIFO_ABA);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_SEGMENTSIZE, llSegmentSize);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_POSTTRIGGER, llPostSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_LOOPS,       llSegmentsToRec);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_ABADIVIDER,  lABADivider);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     0;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
pvGetSegmentDataPointer: recalculates the segment start address for
                         Multiple Recording
**************************************************************************
*/

void* pvGetSegmentDataPointer (ST_SPCM_CARDINFO *pstCardInfo, void* pvDataBuffer, int32 lSegmentsize, int32 lSegmentIdx, int32 lBytesPerSample)
    {
    uint8* pcByteAdr = (uint8*) pvDataBuffer;
    return (void*) &pcByteAdr[lSegmentIdx * lSegmentsize * pstCardInfo->lSetChannels * lBytesPerSample];
    }



/*
**************************************************************************
bSpcMSetupModeRecStdGate: record standard mode gated sampling
**************************************************************************
*/

bool bSpcMSetupModeRecStdGate (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llMemSamples, int64 llPreSamples, int64 llPostSamples)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REC_STD_GATE);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_MEMSIZE,     llMemSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_PRETRIGGER,  llPreSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_POSTTRIGGER, llPostSamples);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     llMemSamples;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupModeRecFIFOGate: record FIFO mode gated sampling
**************************************************************************
*/

bool bSpcMSetupModeRecFIFOGate (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llPreSamples, int64 llPostSamples, int64 llGatesToRec)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REC_FIFO_GATE);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_PRETRIGGER,  llPreSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_POSTTRIGGER, llPostSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_LOOPS,       llGatesToRec);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     0;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupModeRepStdSingle: replay standard mode single
**************************************************************************
*/

bool bSpcMSetupModeRepStdSingle (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llMemSamples)
    {   
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REP_STD_SINGLE);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_MEMSIZE,     llMemSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_LOOPS,       1);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     llMemSamples;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupModeRepStdLoops: replay standard mode looped
**************************************************************************
*/

bool bSpcMSetupModeRepStdLoops (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llMemSamples, int64 llLoops)
    {   
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REP_STD_SINGLE); // with SPC_LOOPS == 0 this will loop continuously
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_LOOPS,       llLoops);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_MEMSIZE,     llMemSamples);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     llMemSamples;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupModeRepStdSingleRestart: replay standard mode Single Restart
**************************************************************************
*/

bool bSpcMSetupModeRepStdSingleRestart (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llMemSamples, int64 llLoops)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REP_STD_SINGLERESTART);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_MEMSIZE,     llMemSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_LOOPS,       llLoops);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     llMemSamples;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupModeRepFIFOSingle: replay FIFO mode single run
**************************************************************************
*/

bool bSpcMSetupModeRepFIFOSingle (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llBlockToRep, int64 llLoopToRep)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;
    
    // check for invalid block/loop combinations
    if ((llBlockToRep && !llLoopToRep) || (!llBlockToRep && llLoopToRep))
        {
        sprintf (pstCardInfo->szError, "bSpcMSetupModeRepFIFOSingle: Loop and Blocks must be either both zero or both defined to non-zero\n");
        return false;
        }

    // segment size can't be zero, we adjust it here
    if (!llBlockToRep && !llLoopToRep)
        llBlockToRep = 1024;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REP_FIFO_SINGLE);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_SEGMENTSIZE, llBlockToRep);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_LOOPS,       llLoopToRep);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     0;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupModeRepStdMulti: replay standard mode Multiple Replay
**************************************************************************
*/

bool bSpcMSetupModeRepStdMulti (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llMemSamples, int64 llSegmentSize, int64 llSegmentsToRep)
    {   
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REP_STD_MULTI);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_MEMSIZE,     llMemSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_SEGMENTSIZE, llSegmentSize);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_LOOPS,       llSegmentsToRep);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     llMemSamples;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupModeRepFIFOMulti: replay FIFO mode Multiple Replay
**************************************************************************
*/

bool bSpcMSetupModeRepFIFOMulti (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llSegmentSize, int64 llSegmentsToRep)
    {   
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REP_FIFO_MULTI);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_SEGMENTSIZE, llSegmentSize);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_LOOPS,       llSegmentsToRep);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     0;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupModeRepStdGate: replay standard mode Gated Replay
**************************************************************************
*/

bool bSpcMSetupModeRepStdGate (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llMemSamples, int64 llGatesToRep)
    {   
       if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REP_STD_GATE);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_MEMSIZE,     llMemSamples);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_LOOPS,       llGatesToRep);

       // store some information in the structure
    pstCardInfo->llSetMemsize =     llMemSamples;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupModeRepFIFOGate: replay FIFO mode Gated Replay
**************************************************************************
*/

bool bSpcMSetupModeRepFIFOGate (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, int64 llGatesToRep)
    {   
       if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,    SPC_REP_FIFO_GATE);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,    qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_LOOPS,       llGatesToRep);

       // store some information in the structure
    pstCardInfo->llSetMemsize =     0;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupModeRepSequence: replay sequence mode
**************************************************************************
*/

bool bSpcMSetupModeRepSequence (ST_SPCM_CARDINFO *pstCardInfo, uint64 qwChEnable, uint32 dwMaxSegments)
    {   
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CARDMODE,            SPC_REP_STD_SEQUENCE);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_CHENABLE,            qwChEnable);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_SEQMODE_MAXSEGMENTS, dwMaxSegments);

    // store some information in the structure
    pstCardInfo->llSetMemsize =     0;
    pstCardInfo->qwSetChEnableMap = qwChEnable;
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_CHCOUNT,     &pstCardInfo->lSetChannels);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
**************************************************************************
**************************************************************************
  Clock setup
**************************************************************************
**************************************************************************
**************************************************************************
*/



/*
**************************************************************************
bSpcMSetupClockPLL: internal clock using PLL
**************************************************************************
*/

bool bSpcMSetupClockPLL (ST_SPCM_CARDINFO *pstCardInfo, int64 llSamplerate, bool bClockOut)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // check for borders
    if (llSamplerate > pstCardInfo->llMaxSamplerate)
        llSamplerate = pstCardInfo->llMaxSamplerate;
    if (llSamplerate < pstCardInfo->llMinSamplerate)
        llSamplerate = pstCardInfo->llMinSamplerate;

    // setup the clock mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CLOCKMODE,   SPC_CM_INTPLL);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_SAMPLERATE,  llSamplerate);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CLOCKOUT,    bClockOut ? 1 : 0);
    if (!dwError) dwError = spcm_dwGetParam_i64 (pstCardInfo->hDrv, SPC_SAMPLERATE, &pstCardInfo->llSetSamplerate);
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_OVERSAMPLINGFACTOR, &pstCardInfo->lOversampling);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupClockQuarz: internal clock using high precision quartz
**************************************************************************
*/

bool bSpcMSetupClockQuartz (ST_SPCM_CARDINFO *pstCardInfo, int64 llSamplerate, bool bClockOut)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // check for borders
    if (llSamplerate > pstCardInfo->llMaxSamplerate)
        llSamplerate = pstCardInfo->llMaxSamplerate;
    if (llSamplerate < pstCardInfo->llMinSamplerate)
        llSamplerate = pstCardInfo->llMinSamplerate;

    // setup the clock mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CLOCKMODE,   SPC_CM_QUARTZ1);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_SAMPLERATE,  llSamplerate);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CLOCKOUT,    bClockOut ? 1 : 0);
    if (!dwError) dwError = spcm_dwGetParam_i64 (pstCardInfo->hDrv, SPC_SAMPLERATE, &pstCardInfo->llSetSamplerate);
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_OVERSAMPLINGFACTOR, &pstCardInfo->lOversampling);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupClockExternal: external clock
**************************************************************************
*/

bool bSpcMSetupClockExternal (ST_SPCM_CARDINFO *pstCardInfo, int32 lExtRange, bool bClockTerm, int32 lDivider)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    if (lDivider > 1)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CLOCKMODE, SPC_CM_EXTDIVIDER);    
    else
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CLOCKMODE, SPC_CM_EXTERNAL);

    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_EXTERNRANGE,     lExtRange);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CLOCKDIV,        lDivider);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CLOCK50OHM,      bClockTerm ? 1 : 0);

    pstCardInfo->llSetSamplerate = 1;
    pstCardInfo->lOversampling =   1;

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupClockRefClock: reference clock
**************************************************************************
*/

bool bSpcMSetupClockRefClock (ST_SPCM_CARDINFO *pstCardInfo, int32 lRefClock, int64 llSamplerate, bool bClockTerm)
    {
    if (!pstCardInfo)
        return false;

    // check for borders
    if (llSamplerate > pstCardInfo->llMaxSamplerate)
        llSamplerate = pstCardInfo->llMaxSamplerate;
    if (llSamplerate < pstCardInfo->llMinSamplerate)
        llSamplerate = pstCardInfo->llMinSamplerate;

    uint32 dwError = ERR_OK;

    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CLOCKMODE,       SPC_CM_EXTREFCLOCK);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_REFERENCECLOCK,  lRefClock);
    if (!dwError) dwError = spcm_dwSetParam_i64 (pstCardInfo->hDrv, SPC_SAMPLERATE,      llSamplerate);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CLOCK50OHM,      bClockTerm ? 1 : 0);
    if (!dwError) dwError = spcm_dwGetParam_i64 (pstCardInfo->hDrv, SPC_SAMPLERATE,     &pstCardInfo->llSetSamplerate);
    if (!dwError) dwError = spcm_dwGetParam_i32 (pstCardInfo->hDrv, SPC_OVERSAMPLINGFACTOR, &pstCardInfo->lOversampling);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
**************************************************************************
**************************************************************************
  Trigger setup
**************************************************************************
**************************************************************************
**************************************************************************
*/


/*
**************************************************************************
bSpcMSetupTrigSoftware: software trigger
**************************************************************************
*/

bool bSpcMSetupTrigSoftware (ST_SPCM_CARDINFO *pstCardInfo, bool bTrigOut)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the trigger mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ORMASK,      SPC_TMASK_SOFTWARE);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ANDMASK,     0);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ORMASK0,  0);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ORMASK1,  0);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ANDMASK0, 0);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ANDMASK1, 0);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIGGEROUT,       bTrigOut ? 1 : 0);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupTrigExternal: external trigger
                        If levels are programmable they'll be set to LVTTL
**************************************************************************
*/

bool bSpcMSetupTrigExternal (ST_SPCM_CARDINFO *pstCardInfo, int32 lExtMode, bool bTrigTerm, int32 lPulsewidth, bool bSingleSrc, int32 lExtLine)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the external trigger mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT0_MODE + lExtLine,       lExtMode);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_TERM,                       bTrigTerm ? 1 : 0);

    // we only use trigout on M2i cards as we otherwise would override the multi purpose i/o lines of M3i, M4i, M4x and M2p
    if (((pstCardInfo->lCardType & TYP_SERIESMASK) == TYP_M2ISERIES) || ((pstCardInfo->lCardType & TYP_SERIESMASK) == TYP_M2IEXPSERIES))
        {
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_OUTPUT,                 0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT0_PULSEWIDTH + lExtLine, lPulsewidth);
        }

    // on bSingleSrc flag no other trigger source is used
    if (bSingleSrc)
        {
        switch (lExtLine)
            {
            case 0 : if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ORMASK, SPC_TMASK_EXT0); break;
            case 1 : if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ORMASK, SPC_TMASK_EXT1); break;
            case 2 : if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ORMASK, SPC_TMASK_EXT2); break;
            case 3 : if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ORMASK, SPC_TMASK_EXT3); break; // X3 on M2p
            }
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ANDMASK,        0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ORMASK0,     0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ORMASK1,     0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ANDMASK0,    0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ANDMASK1,    0);
        }

    // M3i cards need trigger level to be programmed for Ext0 = analog trigger
    if (((pstCardInfo->lCardType & TYP_SERIESMASK) == TYP_M3ISERIES) || ((pstCardInfo->lCardType & TYP_SERIESMASK) == TYP_M3IEXPSERIES))
        {
        if (lExtLine == 0)
            {
            if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT0_LEVEL0, 1500); // 1500 mV
            if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT0_LEVEL1,  800); //  800 mV (rearm)
            }
        }
    // M4i/M4x cards need trigger level to be programmed for Ext0 or Ext1
    else if (((pstCardInfo->lCardType & TYP_SERIESMASK) == TYP_M4IEXPSERIES) || ((pstCardInfo->lCardType & TYP_SERIESMASK) == TYP_M4XEXPSERIES))

        {
        if (lExtLine == 0)
            {
            if ((pstCardInfo->lCardType & (TYP_FAMILYMASK | TYP_CHMASK)) == 0x7700) // single ended 77x0
                {
                if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_THRESHOLD, 1500); // 1500 mV
                }
            else
                {
                if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT0_LEVEL0, 1500);        // 1500 mV
                if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT0_LEVEL1,  800);        //  800 mV (rearm)
                if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT0_ACDC, COUPLING_DC);   // DC coupling
                }
            }
        else if (lExtLine == 1)
            {
            if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT1_LEVEL0, 1500);        // 1500 mV
            if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT1_ACDC, COUPLING_DC);   // DC coupling
            }
        }
    // M2p has a single level on Ext0
    else if ((pstCardInfo->lCardType & TYP_SERIESMASK) == TYP_M2PEXPSERIES)
        {
        if (lExtLine == 0)
            {
            if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT0_LEVEL0, 1500); // 1500 mV
            }
        }

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }

/*
**************************************************************************
bSpcMSetupTrigExternalLevel: external analog trigger with programmable levels
**************************************************************************
*/

bool bSpcMSetupTrigExternalLevel (ST_SPCM_CARDINFO *pstCardInfo, int32 lExtMode, int32 lLevel0, int32 lLevel1, bool bTrigTerm, bool bACCoupling, int32 lPulsewidth, bool bSingleSrc, int32 lExtLine)
    {
    if (!pstCardInfo)
        return false;

    // not supported by M2i and M2i Express cards as they have plain TTL trigger
    if (((pstCardInfo->lCardType & TYP_SERIESMASK) == TYP_M2ISERIES) || ((pstCardInfo->lCardType & TYP_SERIESMASK) == TYP_M2IEXPSERIES))
        return false;

    uint32 dwError = ERR_OK;

    // setup the external trigger mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT0_MODE + lExtLine,       lExtMode);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_TERM,                       bTrigTerm ? 1 : 0);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT0_ACDC + lExtLine,       bACCoupling ? COUPLING_AC : COUPLING_DC);

    // on bSingleSrc flag no other trigger source is used
    if (bSingleSrc)
        {
        switch (lExtLine)
            {
            case 0 : if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ORMASK, SPC_TMASK_EXT0); break;
            case 1 : if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ORMASK, SPC_TMASK_EXT1); break;
            case 2 : if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ORMASK, SPC_TMASK_EXT2); break;
            }
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ANDMASK,        0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ORMASK0,     0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ORMASK1,     0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ANDMASK0,    0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ANDMASK1,    0);
        }

    // M3i cards need trigger level to be programmed for Ext0 = analog trigger
    // M4i/M4x cards need trigger level to be programmed for Ext0 or Ext1
    switch (pstCardInfo->lCardType & TYP_SERIESMASK)
        {
        case TYP_M3ISERIES:
        case TYP_M3IEXPSERIES:
        case TYP_M4IEXPSERIES:
        case TYP_M4XEXPSERIES:
            {
            if (lExtLine == 0)
                {
                if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT0_LEVEL0, lLevel0);
                if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT0_LEVEL1, lLevel1);
                }
            else if (lExtLine == 1)
                {
                if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT1_LEVEL0, lLevel0);
                }
            break;
            }
        case TYP_M2PEXPSERIES:
            {
            if (lExtLine == 0)
                {
                if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_EXT0_LEVEL0, lLevel0);
                }
            break;
            }
        }

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }

/*
**************************************************************************
bSpcMSetupTrigXIO: additional BaseXIO trigger, needs option installed
**************************************************************************
*/

bool bSpcMSetupTrigXIO (ST_SPCM_CARDINFO *pstCardInfo, int32 lXIOMode, bool bSingleSrc, int32 lXIOLine)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    // setup the external trigger mode
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_XIO0_MODE + lXIOLine, lXIOMode);

    // on bSingleSrc flag no other trigger source is used
    if (bSingleSrc)
        {
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ORMASK,         lXIOLine == 0 ? SPC_TMASK_XIO0 : SPC_TMASK_XIO1);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ANDMASK,        0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ORMASK0,     0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ORMASK1,     0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ANDMASK0,    0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ANDMASK1,    0);
        }

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupTrigChannel: channel trigger is set for each channel separately
**************************************************************************
*/

bool bSpcMSetupTrigChannel (ST_SPCM_CARDINFO *pstCardInfo, int32 lChannel, int32 lTrigMode, int32 lTrigLevel0, int32 lTrigLevel1, int32 lPulsewidth, bool bTrigOut, bool bSingleSrc)
    {
    if (!pstCardInfo)
        return false;

    if ((lChannel < 0) || (lChannel >= pstCardInfo->lMaxChannels))
        {
        sprintf (pstCardInfo->szError, "bSpcMSetupTrigChannel: channel number %d not valid. Channels range from 0 to %d\n", lChannel, pstCardInfo->lMaxChannels);
        return false;
        }

    uint32 dwError = ERR_OK;

    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH0_MODE + lChannel,        lTrigMode);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH0_PULSEWIDTH + lChannel,  lPulsewidth);

    if (pstCardInfo->eCardFunction == AnalogIn)
        {
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH0_LEVEL0 + lChannel,  lTrigLevel0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH0_LEVEL1 + lChannel,  lTrigLevel1);
        }

    // we only use trigout on M2i cards as we otherwise would override the multi purpose i/o lines of M3i
    if (((pstCardInfo->lCardType & TYP_SERIESMASK) == TYP_M2ISERIES) || ((pstCardInfo->lCardType & TYP_SERIESMASK) == TYP_M2IEXPSERIES))
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_OUTPUT,                 bTrigOut ? 1 : 0);

    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_TERM,                       0);

    // on bSingleSrc flag no other trigger source is used
    if (bSingleSrc)
        {
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ORMASK,         0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ANDMASK,        0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ORMASK1,     0);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ANDMASK1,    0);

        // some cards need the and mask to use on pulsewidth mode -> to be sure we set the AND mask for all pulsewidth cards
        if ((lTrigMode & SPC_TM_PW_GREATER) || (lTrigMode & SPC_TM_PW_SMALLER))
            {
            if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ORMASK0,     0);
            if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ANDMASK0,    1 << lChannel);
            }
        else
            {
            if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ORMASK0,     1 << lChannel);
            if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ANDMASK0,    0);
            }
        }

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupTrigMask: this function sets the trigger masks (bSingleSrc 
of other commands must be false to use this)
**************************************************************************
*/

bool bSpcMSetupTrigMask (ST_SPCM_CARDINFO *pstCardInfo, uint32 dwChannelOrMask0, uint32 dwChannelOrMask1, uint32 dwChannelAndMask0, uint32 dwChannelAndMask1, uint32 dwTrigOrMask, uint32 dwTrigAndMask)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ORMASK,         dwTrigOrMask);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_ANDMASK,        dwTrigAndMask);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ORMASK0,     dwChannelOrMask0);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ORMASK1,     dwChannelOrMask1);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ANDMASK0,    dwChannelAndMask0);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TRIG_CH_ANDMASK1,    dwChannelAndMask1);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
**************************************************************************
**************************************************************************
  Input + Output channel setup
**************************************************************************
**************************************************************************
**************************************************************************
*/

/*
**************************************************************************
bSpcMSetupInputChannel: allows all input channel related settings
**************************************************************************
*/

bool bSpcMSetupInputChannel (ST_SPCM_CARDINFO *pstCardInfo, int32 lChannel, int32 lInputRange, bool bTerm, int32 lInputOffset, bool bDiffInput)
    {
    if (!pstCardInfo)
        return false;

    if ((lChannel < 0) || (lChannel >= pstCardInfo->lMaxChannels))
        {
        sprintf (pstCardInfo->szError, "SpcMSetupInputChannel: channel number %d not valid. Channels range from 0 to %d\n", lChannel, pstCardInfo->lMaxChannels);
        return false;
        }

    uint32 dwError = ERR_OK;

    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_AMP0 +   lChannel * (SPC_AMP1 - SPC_AMP0),     lInputRange);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_OFFS0 +  lChannel * (SPC_OFFS1 - SPC_OFFS0),   lInputOffset);
    if (pstCardInfo->uCfg.stAI.astPath[0].bInputTermAvailable)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_50OHM0 + lChannel * (SPC_50OHM1 - SPC_50OHM0),  bTerm ? 1 : 0);
    if (pstCardInfo->uCfg.stAI.astPath[0].bDiffModeAvailable)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_DIFF0 + lChannel * (SPC_DIFF1 - SPC_DIFF0),  bDiffInput ? 1 : 0);

    // store some information in the structure
    pstCardInfo->uCfg.stAI.lSetPath[lChannel] =   0;
    pstCardInfo->uCfg.stAI.lSetRange[lChannel] =  lInputRange;
    pstCardInfo->uCfg.stAI.lSetOffset[lChannel] = lInputOffset;

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }


/*
**************************************************************************
bSpcMSetupPathInputCh: M3i version with different settings
**************************************************************************
*/

bool bSpcMSetupPathInputCh (ST_SPCM_CARDINFO *pstCardInfo, int32 lChannel, int32 lPath, int32 lInputRange, int32 lOffset_percent, bool bTerm, bool bACCoupling, bool bBWLimit, bool bDiffInput)
    {
    if (!pstCardInfo)
        return false;

    if ((lChannel < 0) || (lChannel >= pstCardInfo->lMaxChannels))
        {
        sprintf (pstCardInfo->szError, "SpcMSetupInputChannel: channel number %d not valid. Channels range from 0 to %d\n", lChannel, pstCardInfo->lMaxChannels);
        return false;
        }

    uint32 dwError = ERR_OK;

    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_PATH0 +  lChannel * (SPC_PATH1 - SPC_PATH0),   lPath);
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_AMP0 +   lChannel * (SPC_AMP1 - SPC_AMP0),     lInputRange);
    if (pstCardInfo->uCfg.stAI.astPath[lPath].bOffsPercentMode)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_OFFS0 +  lChannel * (SPC_OFFS1 - SPC_OFFS0), lOffset_percent);
    if (pstCardInfo->uCfg.stAI.astPath[lPath].bInputTermAvailable)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_50OHM0 +  lChannel * (SPC_50OHM1  - SPC_50OHM0),   bTerm ?       1 : 0);
    if (pstCardInfo->uCfg.stAI.astPath[lPath].bDiffModeAvailable)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_DIFF0 +   lChannel * (SPC_DIFF1   - SPC_DIFF0),    bDiffInput ?  1 : 0);
    if (pstCardInfo->uCfg.stAI.astPath[lPath].bACCouplingAvailable)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_ACDC0 +   lChannel * (SPC_ACDC1   - SPC_ACDC0),    bACCoupling ? COUPLING_AC : COUPLING_DC);
    if (pstCardInfo->uCfg.stAI.astPath[lPath].bBWLimitAvailable)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_FILTER0 + lChannel * (SPC_FILTER1 - SPC_FILTER0),  bBWLimit ?    1 : 0);

    // store some information in the structure
    pstCardInfo->uCfg.stAI.lSetPath[lChannel] =   lPath;
    pstCardInfo->uCfg.stAI.lSetRange[lChannel] =  lInputRange;
    pstCardInfo->uCfg.stAI.lSetOffset[lChannel] = 0;

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }


/*
**************************************************************************
bSpcMSetupAnalogOutputChannel: allows all output channel related settings
**************************************************************************
*/

bool bSpcMSetupAnalogOutputChannel (ST_SPCM_CARDINFO *pstCardInfo, int32 lChannel, int32 lAmplitude, int32 lOutputOffset, int32 lFilter, int32 lStopMode, bool bDoubleOut, bool bDifferential)
    {
    if (!pstCardInfo)
        return false;

    if ((lChannel < 0) || (lChannel >= pstCardInfo->lMaxChannels))
        {
        sprintf (pstCardInfo->szError, "SpcMSetupAnalogOutputChannel: channel number %d not valid. Channels range from 0 to %d\n", lChannel, pstCardInfo->lMaxChannels);
        return false;
        }

    // Enable output (since M4i).
    uint32 dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_ENABLEOUT0 + lChannel * (SPC_ENABLEOUT1 - SPC_ENABLEOUT0), 1);

    // Check for programmable gain
    if (pstCardInfo->uCfg.stAO.bGainProgrammable)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_AMP0 + lChannel * (SPC_AMP1 - SPC_AMP0), lAmplitude);

    // Check for programmable offset
    if (pstCardInfo->uCfg.stAO.bOffsetProgrammable)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_OFFS0 + lChannel * (SPC_OFFS1 - SPC_OFFS0), lOutputOffset);
    
    // Check for programmable filters
    if (pstCardInfo->uCfg.stAO.bFilterAvailable)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_FILTER0 + lChannel * (SPC_FILTER1 - SPC_FILTER0), lFilter);
    
    // Check for programmable stop levels
    if (pstCardInfo->uCfg.stAO.bStopLevelProgrammable)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CH0_STOPLEVEL + lChannel * (SPC_CH1_STOPLEVEL - SPC_CH0_STOPLEVEL), lStopMode);

    // Check for programmable diffmodes
    if (pstCardInfo->uCfg.stAO.bDiffModeAvailable && !bDoubleOut)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_DIFF0 + lChannel * (SPC_DIFF1 - SPC_DIFF0),  bDifferential ? 1 : 0);

    // Check for programmable doublemodes
    if (pstCardInfo->uCfg.stAO.bDiffModeAvailable && !bDifferential)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_DOUBLEOUT0 + lChannel * (SPC_DOUBLEOUT1 - SPC_DOUBLEOUT0), bDoubleOut ? 1 : 0);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupDigitalOutput: digital output settings for one group
**************************************************************************
*/

bool bSpcMSetupDigitalOutput (ST_SPCM_CARDINFO *pstCardInfo, int32 lGroup, int32 lStopMode, int32 lLowLevel, int32 lHighLevel, bool bDiffMode)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;
    int32 lStopLvlMod = (lGroup * (pstCardInfo->lMaxChannels / pstCardInfo->uCfg.stDIO.lGroups)) >= (pstCardInfo->lMaxChannels / pstCardInfo->lModulesCount) ? 1 : 0;

    if (pstCardInfo->uCfg.stDIO.bStopLevelProgrammable)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_CH0_STOPLEVEL + lStopLvlMod * (SPC_CH1_STOPLEVEL - SPC_CH0_STOPLEVEL), lStopMode);

    if (pstCardInfo->uCfg.stDIO.bOutputLevelProgrammable)
        {
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_LOWLEVEL0 + lGroup * (SPC_LOWLEVEL1 - SPC_LOWLEVEL0), lLowLevel);
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_HIGHLEVEL0 + lGroup * (SPC_HIGHLEVEL1 - SPC_HIGHLEVEL0), lHighLevel);
        }

    if (pstCardInfo->uCfg.stDIO.bDiffModeAvailable)
        {} // to be done


    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
bSpcMSetupDigitalInput: digital input settings for one group
**************************************************************************
*/

bool bSpcMSetupDigitalInput (ST_SPCM_CARDINFO *pstCardInfo,    int32 lGroup, bool bTerm)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;

    if (pstCardInfo->uCfg.stDIO.bInputTermAvailable)
        if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_110OHM0 + lGroup * (SPC_110OHM1 - SPC_110OHM0), bTerm ? 1 : 0);

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }



/*
**************************************************************************
**************************************************************************
**************************************************************************
  Miscellaneous setup
**************************************************************************
**************************************************************************
**************************************************************************
*/



/*
**************************************************************************
bSpcMSetupTimestamp: set up the timestamp mode and performs a 
synchronisation with refernce clock if that mode is activated. Checks for
BASEXIO option if one wants to use reference clock mode
**************************************************************************
*/

bool bSpcMSetupTimestamp (ST_SPCM_CARDINFO *pstCardInfo, int32 lMode, uint32 dwRefTimeoutMS)
    {
    if (!pstCardInfo)
        return false;

    uint32 dwError = ERR_OK;
    bool bRefClockMode = ((lMode & (SPC_TSCNT_REFCLOCKPOS | SPC_TSCNT_REFCLOCKNEG)) != 0);

    // if ref clock is activated for M2i/M3i cards we check for the installation of base xio as this contains the ref clock input
    if (bRefClockMode && ((pstCardInfo->lFeatureMap & SPCM_FEAT_BASEXIO) == 0) && ((pstCardInfo->bM2i) || (pstCardInfo->bM3i)))
        {
        sprintf (pstCardInfo->szError, "Timestamp ref clock mode requires an installed BASEXIO feature!\n");
        return false;
        }

    // set the timestamp mode 
    if (!dwError) dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TIMESTAMP_CMD, lMode);

    // in ref clock mode we now try the synchronisation with external clock
    if (bRefClockMode && !dwError)
        {
        dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TIMESTAMP_TIMEOUT, dwRefTimeoutMS);
        if (!dwError)
            dwError = spcm_dwSetParam_i32 (pstCardInfo->hDrv, SPC_TIMESTAMP_CMD, SPC_TS_RESET);

        // error = synchronisation failed
        if (dwError)
            {
            sprintf (pstCardInfo->szError, "Timestamp reset: synchronisation with external ref clock failed. Check cabeling and check timeout value\n");
            return false;
            }
        }

    return bSpcMCheckSetError (dwError, pstCardInfo);
    }

