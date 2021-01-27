/*
**************************************************************************

sb5_file.cpp                                              (c) Spectrum GmbH

**************************************************************************

implements the header load and write functions

Feel free to use these functions in your own programs.

**************************************************************************
*/

#include "sb5_file.h"

#include <string.h>



/*
****************************************************************************
header type information
****************************************************************************
*/

#define SH_TYP_MASK                 0xFF000000        // mask for type information
#define SH_TYP_INT32                0x01000000        // entry is 32 bit integer
#define SH_TYP_DOUBLE                0x02000000        // entry is 64 bit double 
#define SH_TYP_TEXT                    0x04000000        // entry is zero terminated string
#define SH_TYP_DATA                 0x80000000        // entry is the data area

#define SH_CHX_MASK                    0x00FF0000        // mask for channel index
#define SH_CHX_CH0                    0x00000000        // channel 0
#define SH_CHX_CH1                    0x00010000        // channel 1
#define SH_CHX_CH2                    0x00020000        // channel 2
                                                    // ...


/*
****************************************************************************
header information
****************************************************************************
*/

#define SH_ID_IDMASK                0x0000FFFF        // mask for information id
#define SH_ID_DATA                  0x80000000      // data area
#define SH_ID_XOFFSET                0x02000001        // offset x direction (trigger)
#define SH_ID_XRANGE                0x02000002        // scaling in x direction (1/sampling rate)
#define SH_ID_YRANGE                0x02000003        // scaling in y direction (voltage)
#define SH_ID_YOFFSET                0x02000004        // offset y direction
#define SH_ID_SOURCE                0x04000005      // data source (normally name of generating card)
#define SH_ID_SAMPLES               0x01000006      // number of samples in file
#define SH_ID_SIGNALNAME            0x04000007      // name of signal
#define SH_ID_SIGNALTYP             0x01000008      // type of signal as defined below
#define SH_ID_TIMESTAMP             0x04000009      // time stamp (just system time)
#define SH_ID_MULTISEGMENT          0x0100000C      // segment size on Multiple Recording files
#define SH_ID_SOURCEFS              0x0200000D      // full scale range of source signal (used on FFT display)
#define SH_ID_MUXCOUNT                0x01000015        // number of multiplexed signals in file
#define SH_ID_MUXINDEX                0x01000016        // index of multiplexed signals as sorting may be different



/*
****************************************************************************
SBench file identification
****************************************************************************
*/

#define SB5DATALEN  10
#define SB5DATA        "SB5Data___"
#define SB5STREAM    "SB5Stream_"



/*
****************************************************************************
Some simple read and write routines. All return true if sucessful
****************************************************************************
*/

void vWriteHeader (bool* pbOk, FILE* hFile, int32 lId)
    {
    if (*pbOk)
        (*pbOk) = (fwrite (&lId, sizeof(int32), 1, hFile) == 1);
    }

// *************************************************************************

void vWriteHeader (bool* pbOk, FILE* hFile, int32 lId, int32 lValue)
    {
    if (*pbOk)
        (*pbOk) = (fwrite (&lId, sizeof(int32), 1, hFile) == 1);
    if (*pbOk)
        (*pbOk) = (fwrite (&lValue, sizeof (int32), 1, hFile) == 1);
    }

// *************************************************************************

void vWriteHeader (bool* pbOk, FILE* hFile, int32 lId, double dValue)
    {
    if (*pbOk)
        (*pbOk) = (fwrite (&lId, sizeof(int32), 1, hFile) == 1);
    if (*pbOk)
        (*pbOk) = (fwrite (&dValue, sizeof (double), 1, hFile) == 1);
    }

// *************************************************************************

void vWriteHeader (bool* pbOk, FILE* hFile, int32 lId, const char* pszValue, const char* pszDefault)
    {
    int32 lLen ;

    if (pszValue)
        lLen = (int32) strlen (pszValue) + 1;
    else
        lLen = (int32) strlen (pszDefault) + 1;

    if (*pbOk)
        (*pbOk) = (fwrite (&lId, sizeof(int32), 1, hFile) == 1);
    if (*pbOk)
        (*pbOk) = (fwrite (&lLen, sizeof (int32), 1, hFile) == 1);
    if (*pbOk && pszValue)
        (*pbOk) = (fwrite (pszValue, lLen, 1, hFile) == 1);
    if (*pbOk && !pszValue)
        (*pbOk) = (fwrite (pszDefault, lLen, 1, hFile) == 1);
    }



/*
****************************************************************************
bSB5_StoreHeader: stores the given header information to an already opened 
file. The file pointer is afterwards directly at the beginning of the data
section and plain binary data can be written
****************************************************************************
*/

bool bSB5_StoreHeader (FILE* hFile, ST_SB5HEAD* pstHead)
    {
    bool    bOk = true;

    if (!hFile)
        return false;

    // some simple header checks
    if ((pstHead->lChannels < 0) || (pstHead->lChannels > SB5MAXMUX))
        return false;

    // go to beginning of file and write identification
    bOk = !fseek (hFile, 0, SEEK_SET);
    if (bOk && (pstHead->lChannels == 1))
        bOk = (fwrite (SB5DATA, SB5DATALEN, 1, hFile) == 1);
    if (bOk && (pstHead->lChannels > 1))
        bOk = (fwrite (SB5STREAM, SB5DATALEN, 1, hFile) == 1);

    // multiplex channels
    if (pstHead->lChannels > 1)
        vWriteHeader (&bOk, hFile, SH_ID_MUXCOUNT, pstHead->lChannels);

    // we write samples as second as we may need to update this again
    vWriteHeader (&bOk, hFile, SH_ID_SAMPLES, pstHead->lSumSamples);

    // write data information
    vWriteHeader (&bOk, hFile, SH_ID_SIGNALTYP, pstHead->lSignalType);
    vWriteHeader (&bOk, hFile, SH_ID_MULTISEGMENT, pstHead->lMRSegmentsize);
    vWriteHeader (&bOk, hFile, SH_ID_XOFFSET, pstHead->dXOffset);
    vWriteHeader (&bOk, hFile, SH_ID_XRANGE, pstHead->dXScale);

    // if strins aren't specified we took defaults
    vWriteHeader (&bOk, hFile, SH_ID_SIGNALNAME, pstHead->pszSignalName, "SB5File");
    vWriteHeader (&bOk, hFile, SH_ID_SOURCE, pstHead->pszSource, "<unknown>");
    vWriteHeader (&bOk, hFile, SH_ID_TIMESTAMP, pstHead->pszTimestamp, "<not specified>");

    // write channel related information
    for (int i=0; i<pstHead->lChannels; i++)
        {
        vWriteHeader (&bOk, hFile, SH_ID_MUXINDEX | ((i << 16) & SH_CHX_MASK), pstHead->plMuxIdx[i]);
        vWriteHeader (&bOk, hFile, SH_ID_YRANGE   | ((i << 16) & SH_CHX_MASK), pstHead->pdYScale[i]);
        vWriteHeader (&bOk, hFile, SH_ID_YOFFSET  | ((i << 16) & SH_CHX_MASK), pstHead->pdYOffset[i]);
        vWriteHeader (&bOk, hFile, SH_ID_SOURCEFS | ((i << 16) & SH_CHX_MASK), pstHead->pdSourceFS[i]);
        }

    // at last we write the data identifier
    vWriteHeader (&bOk, hFile, SH_ID_DATA);

    return bOk;
    }



/*
****************************************************************************
Alloc/Free Header. Allocation of header with all arrays correctly set. 
Storage for strings need to be allocated separately
****************************************************************************
*/

ST_SB5HEAD* pstSB5_AllocHeader (int32 lChannels)
    {
    ST_SB5HEAD* pstHeader;

    if ((lChannels < 0) || (lChannels > SB5MAXMUX))
        return NULL;

    pstHeader = new ST_SB5HEAD;
    memset (pstHeader, 0, sizeof (ST_SB5HEAD));

    pstHeader->pdSourceFS =     new double[lChannels];
    pstHeader->pdYOffset =      new double[lChannels];
    pstHeader->pdYScale =       new double[lChannels];
    pstHeader->plMuxIdx =       new int32[lChannels];
    pstHeader->lChannels =      lChannels;

    return pstHeader;
    }

// *************************************************************************

void vSB5_FreeHeader (ST_SB5HEAD* pstHeader)
    {
    if (!pstHeader)
        return;

    delete [] (pstHeader->pdSourceFS);
    delete [] (pstHeader->pdYOffset);
    delete [] (pstHeader->pdYScale);
    delete [] (pstHeader->plMuxIdx);
    if (pstHeader->pszSignalName)
        delete [] (pstHeader->pszSignalName);
    if (pstHeader->pszSource)
        delete [] (pstHeader->pszSource);
    if (pstHeader->pszTimestamp)
        delete [] (pstHeader->pszTimestamp);

    delete (pstHeader);
    }



/*
****************************************************************************
bSB5_UpdateSamples: tries to update the samples settings at the end of 
writing if this information wasn#t known when writing the header information
as usually on streaming files.
****************************************************************************
*/
    
bool bSB5_UpdateSamples (FILE* hFile, int32 lSamples)
    {
    bool    bOk = true;
    int32   lPos;
    int32   plEntry[2];

    // search for samples, it's either first or second
    bOk = !fseek (hFile, SB5DATALEN, SEEK_SET);
    while (bOk)
        {
        lPos = ftell (hFile);
        bOk = (fread (plEntry, sizeof(int32), 2, hFile) == 2);
        if (bOk && plEntry[0] == SH_ID_SAMPLES)
            {
            bOk = !fseek (hFile, lPos, SEEK_SET);
            vWriteHeader (&bOk, hFile, SH_ID_SAMPLES, lSamples);
            return true;
            }
        }

    return false;
    }



/*
****************************************************************************
pstSB5_LoadHeader: load the header information from an already opened file. 
The file pointer is afterwards directly at the beginning of the data
section and plain binary data can be read
****************************************************************************
*/

#define MAXTEXTLEN 1024

ST_SB5HEAD* pstSB5_LoadHeader (FILE* hFile)
    {
    ST_SB5HEAD* pstHeader;
    bool        bOk;
    int32       lChannels = 0;
    int32       lId, lValue, lLen, lCh;
    double      dValue;
    char        szBuffer[MAXTEXTLEN];
    int         i;

    if (!hFile)
        return NULL;

    // examine the identification tag
    bOk = (fread (szBuffer, SB5DATALEN, 1, hFile) == 1);
    if (!bOk)
        return NULL;

    // check for file type
    if (strncmp (szBuffer, SB5STREAM, SB5DATALEN) == 0)
        {
        bOk = (fread (&lId, sizeof(int32), 1, hFile) == 1);
        if (!bOk || (lId != SH_ID_MUXCOUNT))
            return NULL;
        bOk = (fread (&lChannels, sizeof (int32), 1, hFile) == 1);
        }
    else if (strncmp (szBuffer, SB5DATA, SB5DATALEN) == 0)
        lChannels = 1;
    else
        return NULL;

    // check number of channels
    if (!bOk || (lChannels < 0) || (lChannels > SB5MAXMUX))
        return NULL;

    // now we know that the file type is correct and we know which tpye it is
    pstHeader = pstSB5_AllocHeader (lChannels);

    // we predefine the mux signals in case that we don't find them
    for (i=0; i<lChannels; i++)
        pstHeader->plMuxIdx[i] = i;

    // now we loop through the complete header and fill our details
    while (bOk)
        {

        // ----- read the entry -----
        bOk = (fread (&lId, sizeof(int32), 1, hFile) == 1);
        if (bOk)
            switch (lId & SH_TYP_MASK)
                {

                // if we find the data tag, we're done
                case SH_TYP_DATA:
                    return pstHeader;

                // int32 entry
                case SH_TYP_INT32:
                    bOk = (fread (&lValue, sizeof(int32), 1, hFile) == 1);
                    break;

                // double entry
                case SH_TYP_DOUBLE:
                    bOk = (fread (&dValue, sizeof(double), 1, hFile) == 1);
                    break;

                // text entry
                case SH_TYP_TEXT:
                    bOk = (fread (&lLen, sizeof(int32), 1, hFile) == 1);
                    if (lLen > (MAXTEXTLEN - 1))
                        bOk = false;
                    if (bOk)
                        {
                        bOk = (fread (szBuffer, lLen, 1, hFile) == 1);
                        szBuffer[lLen] = 0;
                        }
                    break;

                // unknown entry
                default:
                    bOk = false;
                    break;
                } // switch (lId & SH_TYP_MASK)

        // ----- fill entry in header -----
        if (bOk)
            switch (lId & ~SH_CHX_MASK)
                {
                case SH_ID_XOFFSET:         
                    pstHeader->dXOffset = dValue; 
                    break;

                case SH_ID_XRANGE:          
                    pstHeader->dXScale = dValue; 
                    break;

                case SH_ID_SAMPLES:         
                    pstHeader->lSumSamples = lValue; 
                    break;

                case SH_ID_SIGNALTYP:
                    pstHeader->lSignalType = lValue; 
                    break;

                case SH_ID_MULTISEGMENT:    
                    pstHeader->lMRSegmentsize = lValue; 
                    break;

                case SH_ID_SOURCE:
                    pstHeader->pszSource = new char[strlen(szBuffer) + 1]; 
                    strcpy (pstHeader->pszSource, szBuffer);
                    break;

                case SH_ID_SIGNALNAME:
                    pstHeader->pszSignalName = new char[strlen(szBuffer) + 1]; 
                    strcpy (pstHeader->pszSignalName, szBuffer);
                    break;

                case SH_ID_TIMESTAMP:
                    pstHeader->pszTimestamp = new char[strlen(szBuffer) + 1]; 
                    strcpy (pstHeader->pszTimestamp, szBuffer);
                    break;

                case SH_ID_YRANGE:
                    lCh = (lId >> 16) & 0x0f;
                    if (lCh > (lChannels - 1))
                        bOk = false;
                    else
                        pstHeader->pdYScale[lCh] = dValue;
                    break;

                case SH_ID_YOFFSET:
                    lCh = (lId >> 16) & 0x0f;
                    if (lCh > (lChannels - 1))
                        bOk = false;
                    else
                        pstHeader->pdYOffset[lCh] = dValue;
                    break;

                case SH_ID_SOURCEFS:
                    lCh = (lId >> 16) & 0x0f;
                    if (lCh > (lChannels - 1))
                        bOk = false;
                    else
                        pstHeader->pdSourceFS[lCh] = dValue;
                    break;

                case SH_ID_MUXINDEX:
                    lCh = (lId >> 16) & 0x0f;
                    if (lCh > (lChannels - 1))
                        bOk = false;
                    else
                        pstHeader->plMuxIdx[lCh] = lValue;
                    break;

                default:
                    break;

                } // switch (lId)

        } // while (bOk)


    // something failed
    if (!bOk)
        {
        vSB5_FreeHeader (pstHeader);
        pstHeader = NULL;
        }
    return pstHeader;
    }

