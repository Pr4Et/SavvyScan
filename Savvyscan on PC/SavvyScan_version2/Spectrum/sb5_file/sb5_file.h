/*
**************************************************************************

sb5_file.h                                     (c) Spectrum GmbH , 01/2006

**************************************************************************

offers defintions and functions to handle files of SBench5 format. The
functions handle both *.sb5 data files for single channels and *.sbs data
files for streaming files with multiplexed data.

If using the storage functions please be sure to use the correct file 
names if you wish to open the files under SBench5.

Feel free to use these functions in your own programs.

**************************************************************************
*/

#ifndef SB5_FILE_H
#define SB5_FILE_H

// ----- include standard driver header from library -----
#include "../c_header/dlltyp.h"
#include "../c_header/regs.h"
#include "../c_header/spcerr.h"
#include "../c_header/spcm_drv.h"

// ----- standard c include files -----
#include <stdio.h>



/*
****************************************************************************
signal types
****************************************************************************
*/

#define SIGNAL_TYP_TYPMASK            0xFF000000      // mask for signal type
#define SIGNAL_TYP_ALL                0xFF000000
#define SIGNAL_TYP_ANALOGTIME        0x01000000        // analog signal
#define SIGNAL_TYP_DIGITALTIME        0x02000000        // digital signal
#define SIGNAL_TYP_FFT                0x04000000        // fft signal
#define SIGNAL_TYP_ANADIGTIME        0x10000000        // mixed analog/digital signal, digital bits in upper bits of analog signal
#define SIGNAL_TYP_ANALOGSTREAM        0x20000000        // analog multiplexed stream 
#define SIGNAL_TYP_DIGITALSTREAM    0x40000000        // digital multiplexed stream 
#define SIGNAL_TYP_ANADIGSTREAM        0x80000000        // mixed analog/digital multiplexed stream



/*
****************************************************************************
sample bit width (coded in signaltype)
****************************************************************************
*/

#define SIGNAL_TYP_BITMASK          0x000000FF
#define SIGNAL_TYP_2BIT             0x00000002
#define SIGNAL_TYP_4BIT             0x00000004
#define SIGNAL_TYP_6BIT             0x00000006
#define SIGNAL_TYP_8BIT             0x00000008
#define SIGNAL_TYP_10BIT            0x0000000A
#define SIGNAL_TYP_12BIT            0x0000000C
#define SIGNAL_TYP_14BIT            0x0000000E
#define SIGNAL_TYP_16BIT            0x00000010
#define SIGNAL_TYP_FLOAT            0x00000040



/*
****************************************************************************
sample byte width (coded in signal type)
****************************************************************************
*/

#define SIGNAL_TYP_BYTEMASK         0x00000F00
#define SIGNAL_TYP_1BYTE            0x00000100
#define SIGNAL_TYP_2BYTE            0x00000200
#define SIGNAL_TYP_4BYTE            0x00000400




/*
****************************************************************************
sb5 header information structure
If reading header be sure to have storage space for SB5MAXMUX for all 
channel related information.
****************************************************************************
*/

// maximum number of multiplexed channels
#define SB5MAXMUX   16

struct ST_SB5HEAD
    {
    int32   lSignalType;    // signal type as defined above, be sure to include type, bit width and byte width
    char*   pszSignalName;  // name of the signal
    char*   pszSource;      // source (normally generating card)
    char*   pszTimestamp;   // generating timestamp
    int32   lChannels;      // number of channels stored
    int32   lSumSamples;    // number of stored samples in total (SamplesPerChannel * Channels)
    int32   lMRSegmentsize; // size of one segment in multiple reoording files
    double  dXOffset;       // offset in x direction (trigger position)
    double  dXScale;        // x direction scaling (1/sampling rate)
    double* pdYScale;       // array with y scaling information for each channel
    double* pdYOffset;      // array with y offset information for each channel
    double* pdSourceFS;     // array with full scale ranges of source signal
    int32*  plMuxIdx;       // array with multiplex indexes of each channel
    };



/*
****************************************************************************
bSB5_StoreHeader: stores the given header information to an already opened 
file. The file pointer is afterwards directly at the beginning of the data
section and plain binary data can be written
****************************************************************************
*/

// allocates an empty header structure with correct arrays for channel parameters
ST_SB5HEAD* pstSB5_AllocHeader (int32 lChannels);

// frees the header structure
void vSB5_FreeHeader (ST_SB5HEAD* pstHeader);

// the store function
bool bSB5_StoreHeader (             // returns true if storage succeeded
    FILE*           hFile,          // file handle of an already opened empty file
    ST_SB5HEAD*     pstHeader);     // pointer to a filled header structure, be sure to fill all values correctly



/*
****************************************************************************
bSB5_UpdateSamples: tries to update the samples settings at the end of 
writing if this information wasn#t known when writing the header information
as usually on streaming files.
****************************************************************************
*/
    
bool bSB5_UpdateSamples (           // returns true if update succeeded
    FILE*           hFile,          // file handle of an already opened empty file
    int32           lSamples);      // samples value to update



/*
****************************************************************************
bSB5_LoadHeader: load the header information from an already opened file. 
The file pointer is afterwards directly at the beginning of the data
section and plain binary data can be read
****************************************************************************
*/

ST_SB5HEAD* pstSB5_LoadHeader (     // returns pointer to new allocated and filled header or NULL if an error occurs
    FILE*           hFile);         // file handle of an already opened empty file


#endif
