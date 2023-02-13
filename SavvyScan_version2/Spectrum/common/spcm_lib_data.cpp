/*
**************************************************************************

spcm_lib_data.cpp                                        (c) Spectrum GmbH

**************************************************************************

Offers simple data manipulation routines for the SpcMDrv data format.
Feel free to use this source for own projects and modify it in any kind

Documentation for the API as well as a detailed description of the hardware
can be found in the manual for each device which can be found on our website:
https://www.spectrum-instrumentation.com/en/downloads

Further information can be found online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/knowledge-base-overview

**************************************************************************
*/

// ----- include of common example librarys -----
#include "../common/spcm_lib_data.h"


// ----- standard c include files -----
#include <stdio.h>
#include <string.h>
#include <math.h>


/*
**************************************************************************
bSpcMDemuxDigitalInputData

demultiplexes the digital input data to seperate arrays.
The data buffers for the demultiplexed data must be allocated by the
caller. Each buffer must be of the size (LenInSamples)
**************************************************************************
*/

bool bSpcMDemuxDigitalInputDataToUInt8 (ST_SPCM_CARDINFO *pstCardInfo, void *pvMuxData, uint32 dwLenInSamples, uint8 **ppbyData)
    {
    uint32  dwSample;
    int32   lCh;
    uint8*  ppbyChPtr[SPCM_MAX_AICHANNEL];

    if (!pstCardInfo || !pvMuxData)
        return false;

    // set the sorting table for the channels
    for (lCh=0; lCh < pstCardInfo->lSetChannels; lCh++)
        ppbyChPtr[lCh] = ppbyData[lCh];

    // if two modules are active data is sorted mod0ch0, mod1ch0, mod0ch1, ...
    if (pstCardInfo->qwSetChEnableMap & ~((1 << (pstCardInfo->lMaxChannels / pstCardInfo->lModulesCount)) - 1))
        for (lCh=0; lCh < (pstCardInfo->lSetChannels >> 1); lCh++)
            {
            ppbyChPtr[2 * lCh + 0] = ppbyData[lCh];
            ppbyChPtr[2 * lCh + 1] = ppbyData[(pstCardInfo->lSetChannels >> 1) + lCh];
            }

    uint8* pbyMuxBuf = (uint8*) pvMuxData;

    for (dwSample = 0; dwSample < dwLenInSamples; dwSample++)
        for (lCh = 0; lCh < pstCardInfo->lSetChannels; lCh++)
            *ppbyChPtr[lCh]++ = *pbyMuxBuf++;

    return true;
    }


/*
**************************************************************************
bSpcMDemuxDigitalDataToInt8

demultiplexes the digital channel data to seperate arrays. At 64 bit are
the both inner words swapped to compensate the hardware-order.
The data buffers for the demultiplexed data must be allocated by the
caller.
**************************************************************************
*/

bool bSpcMDemuxDigitalDataToInt8 (ST_SPCM_CARDINFO *pstCardInfo, void *pvMuxData, uint32 dwLenInSamples, int8 **ppbyData)
    {
    int16 nGroupSample;
    int8 bySample;
    int32 lIdxOffset, lChIdx, lSampleIdx;
    uint32 dwGroupIdx;
    int16* ppnBuf;

    ppnBuf = (int16*)pvMuxData;
    lSampleIdx = -1;

    // split data
    for (dwGroupIdx=0; dwGroupIdx < (dwLenInSamples * pstCardInfo->lSetChannels) / 16; dwGroupIdx++)
        {

        nGroupSample = ppnBuf[dwGroupIdx];

        switch (pstCardInfo->lSetChannels)
        {
            // 1, 2, 4, 8 channels
            case 1:
            case 2:
            case 4:
            case 8:
                for (lChIdx=0; lChIdx < 16; lChIdx++)
                    {

                    if (!(lChIdx%pstCardInfo->lSetChannels))
                        lSampleIdx++;

                    bySample = (int8)(nGroupSample & 0x0001);
                    nGroupSample = nGroupSample >> 1;

                    ppbyData[lChIdx%pstCardInfo->lSetChannels][lSampleIdx] = bySample;
                    }

                break;

            // 16 channels
            case 16 :
                for (lChIdx=0; lChIdx < 16; lChIdx++)
                    {
                    bySample = (int8)(nGroupSample & 0x0001);
                    nGroupSample = nGroupSample >> 1;

                    ppbyData[lChIdx][dwGroupIdx] = bySample;
                    }

                break;

            // 32 channels
            case 32 :
                if (!(dwGroupIdx%2))
                    {
                    lIdxOffset = 0;
                    lSampleIdx++;
                    }
                else
                    lIdxOffset = 16;

                for (lChIdx=0; lChIdx < 16; lChIdx++)
                    {
                    bySample = (int8)(nGroupSample & 0x0001);
                    nGroupSample = nGroupSample >> 1;

                    ppbyData[lChIdx + lIdxOffset][lSampleIdx] = bySample;
                    }

                break;

            // 64 channels
            case 64 :
                switch (dwGroupIdx%4)
                    {
                    case 0 :
                        lIdxOffset = 0;
                        lSampleIdx++;
                        break;
                    case 1 :
                        lIdxOffset = 32;
                        break;
                    case 2 :
                        lIdxOffset = 16;
                        break;
                    case 3 :
                        lIdxOffset = 48;
                        break;
                    }

                for (lChIdx=0; lChIdx < 16; lChIdx++)
                    {
                    bySample = (int8)(nGroupSample & 0x0001);
                    nGroupSample = nGroupSample >> 1;

                    ppbyData[lChIdx + lIdxOffset][lSampleIdx] = bySample;
                    }

                break;

            default:
                return false;
            }
        }

    return true;
    }


/*
**************************************************************************
dSpcMIntToVoltage: recalculates an integer value to a voltage value taking
selected range and selected offset into account
**************************************************************************
*/

double dSpcMIntToVoltage (ST_SPCM_CARDINFO *pstCardInfo, int32 lChannel, double dValue)
    {
    if (!pstCardInfo)
        return 0;

    if ((lChannel < 0) || (lChannel >= pstCardInfo->lMaxChannels))
        {
        sprintf (pstCardInfo->szError, "SpcMIntToVoltage: channel number %d not valid. Channels range from 0 to %d\n", lChannel, pstCardInfo->lMaxChannels);
        return 0;
        }

    double dVoltage_mv;

    // recalculate with input range
    dVoltage_mv = pstCardInfo->uCfg.stAI.lSetRange[lChannel] * dValue / pstCardInfo->uCfg.stAI.lMaxADCValue;

    // add the signal offset
    if (pstCardInfo->uCfg.stAI.astPath[pstCardInfo->uCfg.stAI.lSetPath[lChannel]].bOffsPercentMode)
        dVoltage_mv -= pstCardInfo->uCfg.stAI.lSetRange[lChannel] * pstCardInfo->uCfg.stAI.lSetOffset[lChannel] / 100;
    else
        dVoltage_mv -= pstCardInfo->uCfg.stAI.lSetOffset[lChannel];

    return (dVoltage_mv / 1000.0);
    }



/*
**************************************************************************
bSpcMDemuxAnalogDataToVoltage

demultiplexes the analog channel data to seperate arrays.
The data buffers for the demultiplexed data must be allocated by the
caller. Each buffer must be of the size (LenInSamples * sizeof(float))

Recalculates the plain data to voltage levels taking resolution, offset
and range into account.
**************************************************************************
*/

bool bSpcMDemuxAnalogDataToVoltage (ST_SPCM_CARDINFO *pstCardInfo, void *pvMuxData, uint32 dwLenInSamples, float **ppfData)
    {
    uint32  dwSample;
    int32   lCh;
    float*  ppfChPtr[SPCM_MAX_AICHANNEL];
    double  dOffset[SPCM_MAX_AICHANNEL], dFactor[SPCM_MAX_AICHANNEL];

    if (!pstCardInfo || !pvMuxData)
        return false;

    // set the sorting table for the channels
    for (lCh=0; lCh < pstCardInfo->lSetChannels; lCh++)
        ppfChPtr[lCh] = ppfData[lCh];

    // if two modules are active data is sorted mod0ch0, mod1ch0, mod0ch1, ...
    if (pstCardInfo->qwSetChEnableMap & ~((1 << (pstCardInfo->lMaxChannels / pstCardInfo->lModulesCount)) - 1))
        for (lCh=0; lCh < (pstCardInfo->lSetChannels >> 1); lCh++)
            {
            ppfChPtr[2 * lCh + 0] = ppfData[lCh];
            ppfChPtr[2 * lCh + 1] = ppfData[(pstCardInfo->lSetChannels >> 1) + lCh];
            }

    // calculate offset and factor for re-calculation to voltage
    for (lCh=0; lCh < pstCardInfo->lSetChannels; lCh++)
        {
        dFactor[lCh] = (double) pstCardInfo->uCfg.stAI.lSetRange[lCh] / pstCardInfo->uCfg.stAI.lMaxADCValue / 1000.0;

        // add the signal offset
        if (pstCardInfo->uCfg.stAI.astPath[pstCardInfo->uCfg.stAI.lSetPath[lCh]].bOffsPercentMode)
            dOffset[lCh] = -(double) pstCardInfo->uCfg.stAI.lSetRange[lCh] * pstCardInfo->uCfg.stAI.lSetOffset[lCh] / 100.0 / 1000.0;
        else
            dOffset[lCh] = -(double) pstCardInfo->uCfg.stAI.lSetOffset[lCh] / 1000.0;
        }

    // split word data
    if (pstCardInfo->lBytesPerSample > 1)
        {
        int16* pnMuxBuf = (int16*) pvMuxData;

        for (dwSample = 0; dwSample < dwLenInSamples; dwSample++)
            for (lCh = 0; lCh < pstCardInfo->lSetChannels; lCh++)
                *ppfChPtr[lCh]++ = (float) (dOffset[lCh] + dFactor[lCh] * *pnMuxBuf++);
        }

    // split byte data
    else
        {
        int8* pbyMuxBuf = (int8*) pvMuxData;

        for (dwSample = 0; dwSample < dwLenInSamples; dwSample++)
            for (lCh = 0; lCh < pstCardInfo->lSetChannels; lCh++)
                *ppfChPtr[lCh]++ = (float) (dOffset[lCh] + dFactor[lCh] * *pbyMuxBuf++);
        }

    return true;
    }



/*
**************************************************************************
bMMuxData

multiplexes a series of channels into one buffer. The function retains
the information how much bytes one sample has from the CardInfo structure.
The source buffers must be same format and the destination buffer must be
lSetChannels * llMemsize for all the data
**************************************************************************
*/

bool bSpcMMuxData (ST_SPCM_CARDINFO *pstCardInfo, void *pvMuxData, uint32 dwLenInSamples, void **ppvData)
    {
    int32   lCh;
    void*   ppvChPtr[SPCM_MAX_AOCHANNEL];
    uint32  dwSample;

    if (!pstCardInfo || !pvMuxData)
        return false;

    // set the sorting table for the channels
    for (lCh = 0; lCh < pstCardInfo->lSetChannels; lCh++)
        ppvChPtr[lCh] = ppvData[lCh];

    // M2i does not use linear sorting if two modules are active
    if (pstCardInfo->bM2i)
        {
        // if two modules are active data is sorted mod0ch0, mod1ch0, mod0ch1, ...
        if (pstCardInfo->qwSetChEnableMap & ~((1 << (pstCardInfo->lMaxChannels / pstCardInfo->lModulesCount)) - 1))
            {
            for (lCh=0; lCh < (pstCardInfo->lSetChannels >> 1); lCh++)
                {
                ppvChPtr[2 * lCh + 0] = ppvData[lCh];
                ppvChPtr[2 * lCh + 1] = ppvData[(pstCardInfo->lSetChannels >> 1) + lCh];
                }
            }
        }

    // now start the mux loop
    uint32 dwBytesPerSumSample = (uint32) (pstCardInfo->lSetChannels * pstCardInfo->lBytesPerSample);
    uint32 dwBytesPerSample = pstCardInfo->lBytesPerSample;

    for (dwSample = 0; dwSample < dwLenInSamples; dwSample++)
        for (lCh = 0; lCh < pstCardInfo->lSetChannels; lCh++)
            memcpy (((uint8*) pvMuxData) + dwSample * dwBytesPerSumSample + lCh * dwBytesPerSample, ((uint8*) ppvChPtr[lCh]) + dwSample * dwBytesPerSample, dwBytesPerSample);

    return true;
    }



/*
**************************************************************************
bSpcMSplitAnalogAndDigitalData

If synchronous digital inputs have been used with analog data acquisition,
analog and digital data are stored in combined samples, the digital data
using the upper bits of the analog word.

This function splits analog and digital data into separate arrays. Analog
data is sign extended to int16 again to use it with any calculation
routine. Digital data

The caller is responsible to allocate buffer data for the split data. The
analog data buffer must be of the size (LenInSamples * Channels * 2), the
digital data needs (LenInSamples * Channels)
**************************************************************************
*/

bool bSpcMSplitAnalogAndDigitalData (ST_SPCM_CARDINFO *pstCardInfo, void *pvMergedData, uint32 dwLenInSamples, void *pvAnalogData, void *pvDigitalData)
    {
    int16 nSample, nAnalogVal, nDigMask;
    uint32 i;

    if (!pstCardInfo || !pvMergedData)
        return false;

    int16 *pnMergedData   = (int16*)pvMergedData;
    int16 *pnAnalogData   = (int16*)pvAnalogData;
    uint8 *pbyDigitalData = (uint8*)pvDigitalData;

    // split analog and digital part
    for (i = 0; i < dwLenInSamples; i++)
        {
        nSample = pnMergedData[i];

        if (pstCardInfo->uCfg.stAI.lResolution == 12)
            {

            // card resolution = 12 bit -> 4 digital channels
            nDigMask = nSample >> 12;
            nDigMask &= 0x000f;

            if (nSample & 0x800)
                nAnalogVal = nSample | 0xf000;
            else
                nAnalogVal = nSample & 0x0fff;
            }
        else
            {

            // card resolution = 14 bit -> 2 digital channels
            nDigMask = nSample >> 14;
            nDigMask &= 0x0003;

            if (nSample & 0x2000)
                nAnalogVal = nSample | 0xc000;
            else
                nAnalogVal = nSample & 0x3fff;
            }

        pbyDigitalData[i] = (uint8)nDigMask;
        pnAnalogData[i]   = nAnalogVal;
        }

    return true;
    }


/*
**************************************************************************
bFillSB5Header

Fills a SBench5 file header structure with the current setup of pstCard
**************************************************************************
*/

bool bFillSB5Header (ST_SPCM_CARDINFO *pstCard, ST_SB5HEAD *pstHeader, const char* pszName)
    {
    char szSource[25];
    int  i;

    if (!pstCard || !pstHeader)
        return false;

    // check the arrays that we need
    if (!pstHeader->pdSourceFS || !pstHeader->pdYOffset || !pstHeader->pdYScale || !pstHeader->plMuxIdx)
        return false;

    // signal type
    switch (pstCard->eCardFunction)
        {
        case AnalogIn:
            if (pstCard->lSetChannels == 1)
                pstHeader->lSignalType = SIGNAL_TYP_ANALOGTIME;
            else
                pstHeader->lSignalType = SIGNAL_TYP_ANALOGSTREAM;
            switch (pstCard->uCfg.stAI.lResolution)
                {
                case 8:  pstHeader->lSignalType |= SIGNAL_TYP_8BIT | SIGNAL_TYP_1BYTE; break;
                case 12: pstHeader->lSignalType |= SIGNAL_TYP_12BIT | SIGNAL_TYP_2BYTE; break;
                case 14: pstHeader->lSignalType |= SIGNAL_TYP_14BIT | SIGNAL_TYP_2BYTE; break;
                case 16: pstHeader->lSignalType |= SIGNAL_TYP_16BIT | SIGNAL_TYP_2BYTE; break;
                default: break;
                }

            break;

        case DigitalIn:
        case DigitalIO:
            if (pstCard->lSetChannels == 1)
                pstHeader->lSignalType = SIGNAL_TYP_DIGITALTIME;
            else
                pstHeader->lSignalType = SIGNAL_TYP_DIGITALSTREAM;
            break;

        default:
            return false;
        }

    // signal name
    pstHeader->pszSignalName =  new char[strlen(pszName) + 1];
    strcpy (pstHeader->pszSignalName, pszName);


    // source is card name + sn
    switch (pstCard->lCardType & TYP_SERIESMASK)
        {
        case TYP_M2ISERIES:     sprintf (szSource, "M2i.%04x sn %05d", pstCard->lCardType & TYP_VERSIONMASK, pstCard->lSerialNumber); break;
        case TYP_M2IEXPSERIES:  sprintf (szSource, "M2i.%04x-Exp sn %05d", pstCard->lCardType & TYP_VERSIONMASK, pstCard->lSerialNumber); break;
        default:                sprintf (szSource, "Unknown %04x sn %05d", pstCard->lCardType & TYP_VERSIONMASK, pstCard->lSerialNumber); break;
        }
    pstHeader->pszSource =      new char[strlen(szSource) + 1];
    strcpy (pstHeader->pszSource, szSource);


    // fill in the single information
    pstHeader->lSumSamples =    (int32) pstCard->llSetMemsize * pstCard->lSetChannels;
    pstHeader->dXScale =        1.0 / pstCard->llSetSamplerate;
    pstHeader->lChannels =      pstCard->lSetChannels;


    // channel information to be stored
    switch (pstCard->eCardFunction)
        {
        case AnalogIn:
            for (i=0; i<pstCard->lSetChannels; i++)
                {
                pstHeader->pdSourceFS[i] =  pstCard->uCfg.stAI.lSetRange[i];
                pstHeader->pdYOffset[i] =   - (double) pstCard->uCfg.stAI.lSetOffset[i] / 100 * pstCard->uCfg.stAI.lSetRange[i];
                pstHeader->pdYScale[i] =    (double) pstCard->uCfg.stAI.lSetRange[i] / pstCard->uCfg.stAI.lMaxADCValue / 1000.0;
                }

            // if two modules are active data is sorted mod0ch0, mod1ch0, mod0ch1, ...
            if (pstCard->qwSetChEnableMap & ~((1 << (pstCard->lMaxChannels / pstCard->lModulesCount)) - 1))
                for (i=0; i < (pstCard->lSetChannels >> 1); i++)
                    {
                    pstHeader->plMuxIdx[2 * i + 0] = i;
                    pstHeader->plMuxIdx[2 * i + 1] = (pstCard->lSetChannels >> 1) + i;
                    }

            // all channels on one module
            else
                for (i=0; i < pstCard->lSetChannels; i++)
                    pstHeader->plMuxIdx[i] = i;

            break;
        }

    return true;
    }



/*
**************************************************************************
bCalcSignal: calculates simple signal shapes for output card test
**************************************************************************
*/

bool bSpcMCalcSignal (ST_SPCM_CARDINFO *pstCardInfo, void *pvData, uint32 dwLenInSamples, uint32 dwByteWidth, E_SPCM_SIGSHAPE eShape, uint32 dwLoops, uint32 dwGainP)
    {
    int64   llMinFS, llMaxFS, llValue;
    uint32  i;
    int32   lResolution;
    double  dScale;
    int8*   pbyData = (int8*)  pvData;
    int16*  pnData =  (int16*) pvData;
    int32*  plData =  (int32*) pvData;
    int64*  pllData = (int64*) pvData;

    if (!pstCardInfo || !pvData || !dwLenInSamples)
        return false;


    // examine the resolution, bytewidth and min/max values
    switch (pstCardInfo->eCardFunction)
        {
        case AnalogIn:
        case AnalogOut:
            if (pstCardInfo->eCardFunction == AnalogIn)
                lResolution = pstCardInfo->uCfg.stAI.lResolution;
            else
                lResolution = pstCardInfo->uCfg.stAO.lResolution;

            switch (lResolution)
                {
                default:
                case 7:
                case 8:
                    dwByteWidth = 1;
                    break;

                case 12:
                case 14:
                case 16:
                    dwByteWidth = 2;
                    break;
                }

            llMinFS = -pstCardInfo->uCfg.stAO.lMaxDACValue - 1;
            llMaxFS =  pstCardInfo->uCfg.stAO.lMaxDACValue;
            dScale = (double) pstCardInfo->uCfg.stAO.lMaxDACValue * dwGainP / 100.0;

            break;

        case DigitalIn:
        case DigitalOut:
        case DigitalIO:
            if (dwByteWidth == 0)
                {
                sprintf (pstCardInfo->szError, "ByteWidth can't be zero for digital cards as we didn't know how much channels are activated\n");
                return false;
                }
            // two complement numbers
            llMaxFS = (((uint64) 1 << (8 * dwByteWidth - 1)) - 1);
            llMinFS = -llMaxFS - 1;
            dScale = ((double) llMaxFS) * dwGainP / 100.0;
            break;
        }


    // calculation of different signal shapes
    double dSineXScale = 2.0 * 3.14159 / dwLenInSamples * dwLoops;
    uint32 dwBlockLen = dwLenInSamples / dwLoops;
    uint32 dwBlockHalf = dwBlockLen / 2;
    uint32 dwPosInBlock;
    double dSpan = (double) ((uint64) (llMaxFS - llMinFS));
    for (i = 0; i < dwLenInSamples; i++)
        {

        dwPosInBlock = (i % dwBlockLen);

        // calculation of value
        switch (eShape)
            {

            // DC level
            case eDCZero:    llValue = 0; break;
            case eDCPlusFS:  llValue = llMaxFS; break;
            case eDCMinusFS: llValue = llMinFS; break;

            // sine
            case eSine:
            case eInvertedSine:
                llValue = (int64) (dScale * sin (dSineXScale * i));
                break;

            // rectangle
            case eRectangle:
            case eInvertedRectangle:
                if (dwPosInBlock < dwBlockHalf)
                    llValue = llMinFS;
                else
                    llValue = llMaxFS;
                break;

            // triangle
            case eTriangle:
            case eInvertedTriangle:
                if (dwPosInBlock < dwBlockHalf)
                    llValue = (int64) (llMinFS + dwPosInBlock * dSpan / dwBlockHalf);
                else
                    llValue = (int64) (llMaxFS - (dwPosInBlock - dwBlockHalf) * dSpan / dwBlockHalf);
                break;

            // sawtooth
            case eSawtooth:
            case eInvertedSawtooth:
                llValue = (int64) (llMinFS + dwPosInBlock * dSpan / dwBlockLen);
                break;

            default:
                sprintf (pstCardInfo->szError, "Unknown signal shape selected\n");
                return false;
            }

        switch (eShape)
            {
            // invert sign for inverted waveforms
            case eInvertedSine:
            case eInvertedTriangle:
            case eInvertedSawtooth:
            case eInvertedRectangle:
                llValue *= -1;
                break;
            default:
                // nothing
                break;
            }

        // write value to array
        if (llValue < llMinFS)
            llValue = llMinFS;
        else if (llValue > llMaxFS)
            llValue = llMaxFS;

        switch (dwByteWidth)
            {
            default:
            case 1: *pbyData++ = (int8)  llValue; break;
            case 2: *pnData++ =  (int16) llValue; break;
            case 4: *plData++ =  (int32) llValue; break;
            case 8: *pllData++ = (int64) llValue; break;
            }
        } // end of for-loop

    return true;
    }
