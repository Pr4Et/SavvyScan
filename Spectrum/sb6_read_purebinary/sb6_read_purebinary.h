#ifndef SB6_READ_PUREBINARY_H
#define SB6_READ_PUREBINARY_H

#include "../c_header/dlltyp.h"
#include <string>
#include <vector>

struct ChannelInfo
    {
    std::string sName;
    std::string sXUnit;
    std::string sYUnit;
    std::string sDescription;
    double dMaxRange;
    double dMinRange;
    int32 lOrigMaxRange;
    int32 lOrigMinRange;
    int32 lUserOffset;
    };

struct HeaderInfo 
    {
    uint32 dwNumAChannels;
    uint32 dwNumDChannels;
    uint32 dwFileFlags;
    uint32 dwDataEncoding;
    uint32 dwChannelSorting;
    uint32 dwTSSize;
    uint32 dwRawDataFormat;
    uint32 dwStoreDate;
    uint32 dwStoreTime;
    uint64 qwLen;
    uint64 qwPost;
    uint64 qwSegment;
    uint32 dwPretrigger;
    uint32 dwResolution;
    int64  qwSamplerate;
    uint64 qwTrigPos;
    uint64 qwTrigDelay;
    uint64 qwOffset;
    uint32 dwFlags;
    uint32 dwABADivider;
    uint64 qwSlowABALen;
    uint32 dwMaxADCValue;
    uint32 dwTSRefClock;
    uint32 dwOversamplingFactor;
    int64  qwTSSamplerate;
    std::vector <ChannelInfo> vstChInfos;
    };

bool SB6_bReadPureBinaryHeaderInfos (std::string sHeaderFilePath, HeaderInfo &stHeaderInfo);
std::vector <float*> SB6_ReadPureBinaryFileAnalog (std::string sFilePath, HeaderInfo stHeaderInfo, uint64 qwLength = 0, uint64 qwOffset = 0);
std::vector <unsigned char*> SB6_ReadPureBinaryFileDigital (std::string sFilePath, HeaderInfo stHeaderInfo, uint64 qwLength = 0, uint64 qwOffset = 0);
void SB6_vClearChData (std::vector <float*> &vpfChData);
void SB6_vClearChData (std::vector <unsigned char*> &vpbyChData);

#endif // SB6_READ_PUREBINARY_H