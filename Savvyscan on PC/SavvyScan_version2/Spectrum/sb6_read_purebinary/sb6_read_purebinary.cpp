/*
**************************************************************************

sb6_read_purebinary.cpp                                      (c) Spectrum GmbH

**************************************************************************

Example shows how to read a SBench6 pure binary export file
  
Feel free to use this source for own projects and modify it in any kind

Documentation for the API as well as a detailed description of the hardware
can be found in the manual for each device which can be found on our website:
https://www.spectrum-instrumentation.com/en/downloads

Further information can be found online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/knowledge-base-overview

**************************************************************************
*/

#include "sb6_read_purebinary.h"
#include <fstream>
#include <cstring>

using namespace std;

/*
**************************************************************************
SB6_bReadPureBinaryHeaderInfos
**************************************************************************
*/

bool SB6_bReadPureBinaryHeaderInfos (string sHeaderFilePath, HeaderInfo &stHeaderInfo)
    {
    string sLine, sValue; 

    ifstream oFileStream (sHeaderFilePath);
    if (!oFileStream.is_open ())
        return false;

    memset (&stHeaderInfo, 0, sizeof (stHeaderInfo));

    int32 lCurrentChIndex = -1;
    
    while (getline (oFileStream, sLine))
        {
        if (sLine.find ("[Ch") == 0)
            {
            lCurrentChIndex = -1;

            uint32 dwDigitStart = 3;
            uint32 dwDigitEnd   = sLine.find_last_of ("]");

            if (dwDigitEnd < string::npos)
                {
                sValue = sLine.substr (dwDigitStart, dwDigitEnd - dwDigitStart);
                lCurrentChIndex = stol (sValue.c_str ());
                for (int32 lIdx = stHeaderInfo.vstChInfos.size (); lIdx < lCurrentChIndex + 1; lIdx++)
                    stHeaderInfo.vstChInfos.push_back (ChannelInfo ());
                }
            }
           
        if (sLine.find ("Name") == 0 && lCurrentChIndex >= 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 2);
            stHeaderInfo.vstChInfos[lCurrentChIndex].sName = sValue;
            }

        if (sLine.find ("XUnit") == 0 && lCurrentChIndex >= 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 2);
            stHeaderInfo.vstChInfos[lCurrentChIndex].sXUnit = sValue;
            }

        if (sLine.find ("YUnit") == 0 && lCurrentChIndex >= 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 2);
            stHeaderInfo.vstChInfos[lCurrentChIndex].sYUnit = sValue;
            }

        if (sLine.find ("Description") == 0 && lCurrentChIndex >= 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 2);
            stHeaderInfo.vstChInfos[lCurrentChIndex].sDescription = sValue;
            }

        if (sLine.find ("MaxRange") == 0 && lCurrentChIndex >= 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.vstChInfos[lCurrentChIndex].dMaxRange = stod (sValue.c_str ());
            }

        if (sLine.find ("MinRange") == 0 && lCurrentChIndex >= 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.vstChInfos[lCurrentChIndex].dMinRange = stod (sValue.c_str ());
            }
            
         if (sLine.find ("OrigMaxRange") == 0 && lCurrentChIndex >= 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.vstChInfos[lCurrentChIndex].lOrigMaxRange = stol (sValue.c_str ());
            }

          if (sLine.find ("OrigMinRange") == 0 && lCurrentChIndex >= 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.vstChInfos[lCurrentChIndex].lOrigMinRange = stol (sValue.c_str ());
            }

           if (sLine.find ("UserOffset") == 0 && lCurrentChIndex >= 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.vstChInfos[lCurrentChIndex].lUserOffset = stol (sValue.c_str ());
            }

        if (sLine.find ("NumAChannels") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwNumAChannels = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("NumDChannels") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwNumDChannels = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("FileFlags") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwFileFlags = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("DataEncoding") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwDataEncoding = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("ChannelSorting") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwChannelSorting = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("TSSize") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwTSSize = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("RawDataFormat") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwRawDataFormat = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("StoreDate") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwStoreDate = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("StoreTime") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwStoreTime = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("LenH") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwLen = (static_cast <uint64> (stoul (sValue.c_str ())) << 32) | (stHeaderInfo.qwLen & 0xFFFFFFFF);
            continue;
            }

        if (sLine.find ("LenL") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwLen = (stHeaderInfo.qwLen & 0xFFFFFFFF00000000) | static_cast <uint64> (stoul (sValue.c_str ()));
            continue;
            }

        if (sLine.find ("PostH") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwPost = (static_cast <uint64> (stoul (sValue.c_str ())) << 32) | (stHeaderInfo.qwPost & 0xFFFFFFFF);
            continue;
            }

        if (sLine.find ("PostL") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwPost = (stHeaderInfo.qwPost & 0xFFFFFFFF00000000) | static_cast <uint64> (stoul (sValue.c_str ()));
            continue;
            }

        if (sLine.find ("SegmentH") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwSegment = (static_cast <uint64> (stoul (sValue.c_str ())) << 32) | (stHeaderInfo.qwSegment & 0xFFFFFFFF);
            continue;
            }

        if (sLine.find ("SegmentL") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwSegment = (stHeaderInfo.qwSegment & 0xFFFFFFFF00000000) | static_cast <uint64> (stoul (sValue.c_str ()));
            continue;
            }

        if (sLine.find ("Pretrigger") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwPretrigger = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("Resolution") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwResolution = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("Samplerate") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwSamplerate = stoull (sValue.c_str ());
            continue;
            }

        if (sLine.find ("TrigPosH") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwTrigPos = (static_cast <uint64> (stoul (sValue.c_str ())) << 32) | (stHeaderInfo.qwTrigPos & 0xFFFFFFFF);
            continue;
            }

        if (sLine.find ("TrigPosL") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwTrigPos = (stHeaderInfo.qwTrigPos & 0xFFFFFFFF00000000) | static_cast <uint64> (stoul (sValue.c_str ()));
            continue;
            }

        if (sLine.find ("TrigDelayH") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwTrigDelay = (static_cast <uint64> (stoul (sValue.c_str ())) << 32) | (stHeaderInfo.qwTrigDelay & 0xFFFFFFFF);
            continue;
            }

        if (sLine.find ("TrigDelayL") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwTrigDelay = (stHeaderInfo.qwTrigDelay & 0xFFFFFFFF00000000) | static_cast <uint64> (stoul (sValue.c_str ()));
            continue;
            }

        if (sLine.find ("OffsetH") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwOffset = (static_cast <uint64> (stoul (sValue.c_str ())) << 32) | (stHeaderInfo.qwOffset & 0xFFFFFFFF);
            continue;
            }

        if (sLine.find ("OffsetL") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwOffset = (stHeaderInfo.qwOffset & 0xFFFFFFFF00000000) | static_cast <uint64> (stoul (sValue.c_str ()));
            continue;
            }

        if (sLine.find ("Flags") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwFlags = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("ABADivider") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwABADivider = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("SlowABALenH") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwSlowABALen = (static_cast <uint64> (stoul (sValue.c_str ())) << 32) | (stHeaderInfo.qwSlowABALen & 0xFFFFFFFF);
            continue;
            }

        if (sLine.find ("SlowABALenL") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwSlowABALen = (stHeaderInfo.qwSlowABALen & 0xFFFFFFFF00000000) | static_cast <uint64> (stoul (sValue.c_str ()));
            continue;
            }

        if (sLine.find ("MaxADCValue") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwMaxADCValue = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("TSRefClock") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwTSRefClock = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("OversamplingFactor") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.dwOversamplingFactor = stoul (sValue.c_str ());
            continue;
            }

        if (sLine.find ("TSSamplerate") == 0)
            {
            sValue = sLine.substr (sLine.find_last_of ("=") + 1);
            stHeaderInfo.qwTSSamplerate = stoull (sValue.c_str ());
            continue;
            }
        }

    return true;
    }

/*
**************************************************************************
SB6_ReadPureBinaryFileAnalog
**************************************************************************
*/

vector <float*> SB6_ReadPureBinaryFileAnalog (string sFilePath, HeaderInfo stHeaderInfo, uint64 qwLengthSamples, uint64 qwOffsetSamples)
    {
    char *pbyDataBuffer = NULL;
    int16 *pnDataBuffer = NULL;
    uint64 qwLengthPerCh = 0;
    
    vector <float*> vpfChData;

    ifstream oFileStream (sFilePath, ios::binary);
    if (oFileStream.is_open ())
        {
        // get number of channels from header file
        uint32 dwNumCh = stHeaderInfo.dwNumAChannels;
           
        // check number of available channel infos
        if (dwNumCh > stHeaderInfo.vstChInfos.size ())
            return vpfChData;

        // get resolution info from header file
        uint32 dwBytesPerSample = 1;
        if (stHeaderInfo.dwResolution > 8)
            dwBytesPerSample = 2;

        // if length parameter is not set read complete data
        if (qwLengthSamples > 0)
            qwLengthPerCh = qwLengthSamples;
        else
            qwLengthPerCh = stHeaderInfo.qwLen;

        // adjust MaxADCValue if this is not set
        if (!stHeaderInfo.dwMaxADCValue)
            stHeaderInfo.dwMaxADCValue = 1;

        // calculate buffer size in bytes
        uint64 qwBufferLength = dwNumCh * dwBytesPerSample * qwLengthPerCh;

        // allocate buffer memory
        pbyDataBuffer = new char[(size_t)qwBufferLength];
    
        // allocate memory for each channel
        for (uint32 dwChIdx = 0; dwChIdx < dwNumCh; dwChIdx++)
            vpfChData.push_back (new float[(size_t)qwLengthPerCh]);

        // set data offset
        if (qwOffsetSamples > 0)
            {
            // calculate offset in bytes
            uint64 qwOffsetBytes = dwBytesPerSample * dwNumCh * qwOffsetSamples;
            oFileStream.seekg (qwOffsetBytes, ios::beg);
            }

        // read data from file
        oFileStream.read (pbyDataBuffer, qwBufferLength);

        uint64 qwDataIdx = 0;
        
        pnDataBuffer = (int16*)pbyDataBuffer;
        for (uint64 qwBufferIdx = 0; qwBufferIdx < qwBufferLength / (uint64)dwBytesPerSample; qwBufferIdx += dwNumCh)
            {
            // calculate voltage values for each channel
            for (uint32 dwChIdx = 0; dwChIdx < dwNumCh; dwChIdx++)
                {
                if (dwBytesPerSample == 2)
                    vpfChData[dwChIdx][qwDataIdx] = (float)(pnDataBuffer[qwBufferIdx + dwChIdx] * stHeaderInfo.vstChInfos[dwChIdx].lOrigMaxRange) / (float)stHeaderInfo.dwMaxADCValue;
                else
                    vpfChData[dwChIdx][qwDataIdx] = (float)(pbyDataBuffer[qwBufferIdx + dwChIdx] * stHeaderInfo.vstChInfos[dwChIdx].lOrigMaxRange) / (float)stHeaderInfo.dwMaxADCValue;
                }
                    
            qwDataIdx++;
            }
          
        // free buffer memory
        delete[] pbyDataBuffer;
        }
    
    return vpfChData;
    }

/*
**************************************************************************
SB6_ReadPureBinaryFileDigital
**************************************************************************
*/

vector <unsigned char*> SB6_ReadPureBinaryFileDigital (string sFilePath, HeaderInfo stHeaderInfo, uint64 qwLengthSamples, uint64 qwOffsetSamples)
    {
    char *pbyDataBuffer = NULL;
    uint64 qwLengthPerCh = 0;
    uint32 dwNumCh = 0;
    uint32 dwDigitalGroup = 0;
    
    vector <unsigned char*> vpbyChData;

    ifstream oFileStream (sFilePath, ios::binary);
    if (oFileStream.is_open ())
        {
        // get number of channels from header file
        dwNumCh = stHeaderInfo.dwNumDChannels;
        dwDigitalGroup = dwNumCh / 8;

        // check number of available channel infos
        if (dwNumCh > stHeaderInfo.vstChInfos.size ())
            return vpbyChData;

        // if length parameter is not set read complete data
        if (qwLengthSamples > 0)
            qwLengthPerCh = qwLengthSamples;
        else
            qwLengthPerCh = stHeaderInfo.qwLen;

        // calculate buffer size in bytes
        uint64 qwBufferLength = dwDigitalGroup * qwLengthPerCh;
        
        // allocate buffer memory
        pbyDataBuffer = new char[(size_t)qwBufferLength];
    
        // allocate memory for each digital channel group
        for (uint32 dwGroupIdx = 0; dwGroupIdx < dwDigitalGroup; dwGroupIdx++)
            vpbyChData.push_back (new unsigned char[(size_t)qwLengthPerCh]);

        // set data offset
        if (qwOffsetSamples > 0)
            {
            // calculate offset in bytes
            uint64 qwOffsetBytes = dwDigitalGroup * qwOffsetSamples;
            oFileStream.seekg (qwOffsetBytes, ios::beg);
            }

        // read data from file
        oFileStream.read (pbyDataBuffer, qwBufferLength);

        uint64 qwDataIdx = 0;

        for (uint64 qwBufferIdx = 0; qwBufferIdx < qwBufferLength; qwBufferIdx += dwDigitalGroup)
            {
            for (uint32 dwGroupIdx = 0; dwGroupIdx < dwDigitalGroup; dwGroupIdx++)
                vpbyChData[dwGroupIdx][qwDataIdx] = pbyDataBuffer[qwBufferIdx + dwGroupIdx];
            
            qwDataIdx++;
            }

        // free buffer memory
        delete[] pbyDataBuffer;
        }

    return vpbyChData;
    }

/*
**************************************************************************
SB6_vClearChData
**************************************************************************
*/

void SB6_vClearChData (vector <float*> &vpfChData)
    {
    // free memory for each channel
    for (uint32 dwChIdx = 0; dwChIdx < vpfChData.size (); dwChIdx++)
        delete[] vpfChData[dwChIdx];

    vpfChData.clear ();
    }

void SB6_vClearChData (vector <unsigned char*> &vpbyChData)
    {
    // free memory for each channel
    for (uint32 dwChIdx = 0; dwChIdx < vpbyChData.size (); dwChIdx++)
        delete[] vpbyChData[dwChIdx];

    vpbyChData.clear ();
    }
