/*
**************************************************************************

spcm_lib_card.h                                (c) Spectrum GmbH , 01/2006

**************************************************************************

defines the library functions as external to use them with projects where
the lib is directly included (like dll) and where the lib has to be loaded
separately (like LabWindows or unsupported C compiler)

**************************************************************************
*/

#ifndef SPCM_LIB_CARD_H
#define SPCM_LIB_CARD_H

// ----- include standard driver header from library -----
#include "../c_header/dlltyp.h"
#include "../c_header/regs.h"
#include "../c_header/spcerr.h"
#include "../c_header/spcm_drv.h"



#define SPCM_MAX_AIRANGE    8
#define SPCM_MAX_AICHANNEL  16
#define SPCM_MAX_AIPATH     2
#define SPCM_MAX_AOCHANNEL  8
#define SPCM_MAX_DIOCHANNEL 64


/*
**************************************************************************
structure with different card information and setup information that is
used by the setup and data processing routines
**************************************************************************
*/

// different card functionalities
typedef enum E_SPCM_CARDFNC {AnalogIn, AnalogOut, DigitalOut, DigitalIn, DigitalIO} E_SPCM_CARDFNC;

// card information structure
typedef struct 
    {

    // information from the card
    drv_handle      hDrv;                   // handle to opened card driver
    int32           lCardIdx;               // index of card (from open), just for display
    int32           lCardType;              // card type as listed in the manual
    E_SPCM_CARDFNC  eCardFunction;          // function of the card
    int32           lSerialNumber;          // serial number of card
    int64           llInstMemBytes;         // installed on-board memory in bytes
    int32           lFeatureMap;            // bitmap with installed card features
    int32           lExtFeatureMap;         // bitmap with installed extended features
    int32           lMaxChannels;           // number of channels (analog or digital)
    int32           lModulesCount;          // number of installed modules for data sorting algorithm
    int32           lBytesPerSample;        // number of bytes for each sample (analog data)
    int64           llMinSamplerate;        // minimum sampling rate
    int64           llMaxSamplerate;        // maximum sampling rate

    int32           lLibVersion;            // version of the library
    int32           lKernelVersion;         // version of the kernel driver

    int32           lCtrlFwVersion;         // version of main control firmware
    int32           lBaseHwVersion;         // version of base hardware
    int32           lModHwVersion;          // version of module hardware
    int32           lModFwVersion;          // version of module firmware

    // current settings
    bool            bSetError;              // one of the functions generated an error
    int32           lErrorCode;             // error code
    char            szError[ERRORTEXTLEN];  // space for the error text
    uint64          qwSetChEnableMap;       // current channel enable map
    int64           llSetMemsize;           // programmed memory size
    int32           lSetChannels;           // number of used channels for this run
    int64           llSetSamplerate;        // current selected sampling rate (1 for external)
    int32           lOversampling;          // currently active oversampling factor

    // flags for the examples to determine card family
    bool            bM2i;                   // M2i.xxxx or M2i.xxxx-exp
    bool            bM3i;                   // M3i.xxxx or M3i.xxxx-exp
    bool            bM4i;                   // M4i.xxxx-x8 or M4x.xxxx-x4
    bool            bM2p;                   // M2p.xxxx-x4
    bool            bRemote;                // Netbox or Remote Server

    // card function dependant details
    union
        {

        // analog input cards
        struct ST_SPCM_AI
            {
            int32   lResolution;                    // resolution of analog channels
            int32   lMaxADCValue;                   // maximum range, normally 2^(Resolution-1) but can be limited
            int32   lPathCount;                     // number of input paths
            struct ST_SPCM_AI_PATH                  // the different paths may have different features
                {
                int32   lRangeCount;                    // number of analog input ranges
                int32   lRangeMin[SPCM_MAX_AIRANGE];    // analog input ranges
                int32   lRangeMax[SPCM_MAX_AIRANGE];    // ...
                bool    bInputTermAvailable;            // input termination available
                bool    bDiffModeAvailable;             // differential mode available
                bool    bACCouplingAvailable;           // AC/DC coupling softwar selectable
                bool    bBWLimitAvailable;              // bandwidth limit available
                bool    bOffsPercentMode;               // offset programmed in percent of range
                } astPath[SPCM_MAX_AIPATH];

            int32   lSetPath[SPCM_MAX_AICHANNEL];
            int32   lSetRange[SPCM_MAX_AICHANNEL];  // current used input range for each channel
            int32   lSetOffset[SPCM_MAX_AICHANNEL]; // current set input offset
            } stAI;

        // analog output cards
        struct ST_SPCM_AO
            {
            int32   lResolution;                    // resolution of analog channels
            int32   lMaxDACValue;                   // maximum range, normally 2^(Resolution-1) but can be limited
            bool    bGainProgrammable;              // programmable gain available 
            bool    bOffsetProgrammable;            // programmable offset available
            bool    bFilterAvailable;               // programmable filters available
            bool    bStopLevelProgrammable;         // programmable stop level available
            bool    bDiffModeAvailable;             // differential mode available
            } stAO;

        // digital input, outputs or i/o cards
        struct ST_SPCM_DIO
            {
            int32   lGroups;                        // number of channel groups that have individual setup
            bool    bInputTermAvailable;            // input termination available
            bool    bDiffModeAvailable;             // differential mode available
            bool    bStopLevelProgrammable;         // programmable stop level available
            bool    bOutputLevelProgrammable;       // low and high output level is programmable
            } stDIO;
        
        } uCfg;
    } ST_SPCM_CARDINFO;





/*
**************************************************************************
bSpcMInitCardByIdx:

opens the driver with the given indes, reads out card information and
fills the CARDINFO structure
**************************************************************************
*/

bool bSpcMInitCardByIdx (               // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,    // pointer to an allocated and empty card info structure
    int32               lCardIdx);      // index of card to open, index starts with zero

bool bSpcMInitCardByIdx (               // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,    // pointer to an allocated and empty card info structure
    const char*         szDrvName,      // name of card like /dev/spcm0 (local) or 123.123.123.123 (for digitizerNETBOX only)
    int32               lCardIdx);      // index of card to open, index starts with zero


/*
**************************************************************************
nErrorMessageStdOut:

prints the error message to std out and ends the driver if it's active
program can be left with this function
**************************************************************************
*/

int nSpcMErrorMessageStdOut (                   // returns -1
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    const char*         pszMessage,             // user error message,
    bool                bPrintCardErr = true);  // add card error message



/*
**************************************************************************
pszSpcMTranslateRuntimeError

translation of a runtime error code into a text message. The buffer need 
to be at least ERRORTEXTLEN long to cover any current or future messages
**************************************************************************
*/

char* pszSpcMTranslateRuntimeError (
    uint32 dwErrorCode, 
    char* pszBuffer);


/*
**************************************************************************
vSpcMCloseCard

closes the driver
**************************************************************************
*/

void vSpcMCloseCard (
    ST_SPCM_CARDINFO   *pstCardInfo);           // pointer to a filled card info structure



/*
**************************************************************************
pszSpcMPrintCardInfo

prints the card information to a string for display.
**************************************************************************
*/

char* pszSpcMPrintCardInfo (                    // returns the pointer to the printed string
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    char*               pszBuffer,              // buffer for printing
    int32               lStrLen,                // length of the buffer
    bool                bExtended = true);      // extended info, if false only name+sn

char* pszSpcMPrintDocumentationLink (
    const ST_SPCM_CARDINFO *pstCardInfo,        // pointer to a filled card info structure
    char*                   pszBuffer,          // buffer for printing
    int32                   lStrLen);           // length of the buffer


/*
**************************************************************************
bSpcMSetupModeXXX

setup one of the card modes
**************************************************************************
*/

// record standard mode single
bool bSpcMSetupModeRecStdSingle (               // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next acquisition
    int64               llMemSamples,           // recording length in samples per channel
    int64               llPostSamples);         // samples to record after trigger event

// record FIFO mode single
bool bSpcMSetupModeRecFIFOSingle (              // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next acquisition
    int64               llPreSamples,           // number of samples to be stored before the trigger event
    int64               llBlockToRec = 0,       // blocks and loops can define the maximum recording length
    int64               llLoopToRec = 0);       // in FIFO mode as Block * Loop. If zero we run continuously

// ***********************************************************************


// record standard mode average
bool bSpcMSetupModeRecStdAverage (              // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next acquisition
    int64               llMemSamples,           // recording length in samples per channel
    int64               llSegmentSize,          // size of each multiple recording segment
    int64               llPostSamples,          // samples to record after trigger event for each segment
    int32               lAverages               // number of triggered segments to average
    );         

// record standard mode multiple recording
bool bSpcMSetupModeRecStdMulti (                // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next acquisition
    int64               llMemSamples,           // recording length in samples per channel
    int64               llSegmentSize,          // size of each multiple recording segment
    int64               llPostSamples);         // samples to record after trigger event for each segment

// record standard mode ABA
bool bSpcMSetupModeRecStdABA (                  // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next acquisition
    int64               llMemSamples,           // recording length in samples per channel
    int64               llSegmentSize,          // size of each multiple recording segment
    int64               llPostSamples,          // samples to record after trigger event for each segment
    int32               lABADivider);           // divider for ABA mode slow samples

// record FIFO mode average
bool bSpcMSetupModeRecFIFOAverage (             // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next acquisition
    int64               llSegmentSize,          // size of each multiple recording segment
    int64               llPostSamples,          // samples to record after trigger event for each segment
    int32               lAverages,
    int64               llSegmentsToRec = 0);   // number of segments to record in total. If zero we reun continuously

// record FIFO mode multiple recording
bool bSpcMSetupModeRecFIFOMulti (               // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next acquisition
    int64               llSegmentSize,          // size of each multiple recording segment
    int64               llPostSamples,          // samples to record after trigger event for each segment
    int64               llSegmentsToRec = 0);   // numbe of segments to record in total. If zero we reun continuously

// record FIFO mode ABA
bool bSpcMSetupModeRecFIFOABA (                 // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next acquisition
    int64               llSegmentSize,          // size of each multiple recording segment
    int64               llPostSamples,          // samples to record after trigger event for each segment
    int32               lABADivider,            // divider for ABA mode slow samples
    int64               llSegmentsToRec = 0);   // numbe of segments to record in total. If zero we reun continuously

// recalculates the data start address of segment no. idx
void* pvGetSegmentDataPointer (                 // returns an pointer to the segment start address
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    void*               pvDataBuffer,           // pointer to the data array that holds all segments
    int32               lSegmentsize,           // size of one segment
    int32               lSegmentIdx,            // index of the segment of which we wish to get the pointer
    int32               lBytesPerSample);       // number of bytes per sample    

// ***********************************************************************

// record standard mode gated sampling
bool bSpcMSetupModeRecStdGate (                 // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next acquisition
    int64               llMemSamples,           // recording length in samples per channel
    int64               llPreSamples,           // number of samples to record before gate starts
    int64               llPostSamples);         // number of samples to record after gate ends

// record FIFO mode gated sampling
bool bSpcMSetupModeRecFIFOGate (                // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next acquisition
    int64               llPreSamples,           // number of samples to record before gate starts
    int64               llPostSamples,          // number of samples to record after gate ends
    int64               llGatesToRec = 0);      // number of gates to record


// ***********************************************************************

// replay standard mode single
bool bSpcMSetupModeRepStdSingle (               // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next generation
    int64               llMemSamples);          // samples to replay after trigger event

// replay standard mode looped
bool bSpcMSetupModeRepStdLoops (                // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next generation
    int64               llMemSamples,           // samples to replay after trigger event
    int64               llLoops = 0);           // loops to replay (0 --> infinite continuous replay)

// replay standard mode single restart
bool bSpcMSetupModeRepStdSingleRestart (        // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next generation
    int64               llMemSamples,           // samples to replay after trigger event
    int64               llLoops = 0);           // loops to replay (0 --> infinite continuous replay)

// replay FIFO mode single
bool bSpcMSetupModeRepFIFOSingle (              // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next acquisition
    int64               llBlockToRep = 0,       // blocks and loops can define the maximum replay length
    int64               llLoopToRep = 0);       // in FIFO mode as Block * Loop. If zero we run continuously

// ***********************************************************************

// standard mode multiple replay
bool bSpcMSetupModeRepStdMulti (                // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next generation
    int64               llMemSamples,           // replay length in samples per channel
    int64               llSegmentSize,          // size of each segment
    int64               llSegmentsToRep = 1);   // segments to replay (0 = infinite, 1 = memsize once, N = number of segments)

// FIFO mode multiple replay
bool bSpcMSetupModeRepFIFOMulti (               // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next generation
    int64               llSegmentSize,          // size of each segment
    int64               llSegmentsToRep = 0);   // segments to replay (0 = infinite, N = number of segments)

// ***********************************************************************

// standard mode gated replay 
bool bSpcMSetupModeRepStdGate (                 // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next generation
    int64               llMemSamples,           // recording length in samples per channel
    int64               llGatesToRep = 1);      // gates to replay (0 = infinte, 1 = memsize once, N = number of gates)

// FIFO mode gated replay 
bool bSpcMSetupModeRepFIFOGate (                // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next generation
    int64               llGatesToRep = 0);      // gates to replay (0 = infinte, N = number of gates)

// ***********************************************************************

bool bSpcMSetupModeRepSequence (                // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint64              qwChEnable,             // channel enable mask for the next generation
    uint32              dwMaxSegments);         // count of divided main sample memory segments

// ***********************************************************************



/*
**************************************************************************
bSpcMSetupClockXXX

setup the clock engine for different modes
**************************************************************************
*/

// internal clock using PLL
bool bSpcMSetupClockPLL (                       // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    int64               llSamplerate,           // desired sampling rate 
    bool                bClockOut = false);     // clock output enable

// internal clock using high precision quartz
bool bSpcMSetupClockQuartz (                    // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    int64               llSamplerate,           // sampling rate if internal clock mode
    bool                bClockOut = false);     // clock output enable

// external clock 
bool bSpcMSetupClockExternal (                  // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    int32               lExtRange,              // external clock range if external clock mode is used
    bool                bClockTerm = true,      // enable clock termination (50 ohm)
    int32               lDivider = 1);          // clock divider

// reference clock
bool bSpcMSetupClockRefClock (                  // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    int32               lRefClock,              // reference clock speed if using ref clock mode
    int64               llSamplerate,           // desired sampling rate
    bool                bClockTerm = true);     // enable clock termination (50 ohm)



/*
**************************************************************************
bSpcMSetupTriggerXXX

setup the trigger engine for different modes
**************************************************************************
*/

// software trigger
bool bSpcMSetupTrigSoftware (                   // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    bool                bTrigOut = false);      // enable trigger output

// external trigger (if input is using comparators the levels are set to TTL)
bool bSpcMSetupTrigExternal (                   // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    int32               lExtMode,               // external trigger mode
    bool                bTrigTerm = true,       // trigger termination active
    int32               lPulsewidth = 0,        // programmable pulsewidth for all external + pulsewidth modes
    bool                bSingleSrc = true,      // acts as single trigger source, all other masks cleared
    int32               lExtLine = 0);          // standard external trigger is line 0

// external analog trigger with programmable levels
bool bSpcMSetupTrigExternalLevel (              // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    int32               lExtMode,               // external trigger mode
    int32               lLevel0 = 1500,         // trigger level 0 (mV)
    int32               lLevel1 = 800,          // trigger level 1 (mV)
    bool                bTrigTerm = true,       // trigger termination active
    bool                bACCoupling = false,    // programmable AC coupling
    int32               lPulsewidth = 0,        // programmable pulsewidth for all external + pulsewidth modes
    bool                bSingleSrc = true,      // acts as single trigger source, all other masks cleared
    int32               lExtLine = 0);          // standard external trigger is line 0

// additional BaseXIO trigger
bool bSpcMSetupTrigXIO (                        // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    int32               lXIOMode,               // external trigger mode
    bool                bSingleSrc = true,      // acts as single trigger source, all other masks cleared
    int32               lXIOLine = 0);          // standard XIO trigger is line 0

// channel trigger is set for each channel separately
bool bSpcMSetupTrigChannel (                    // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    int32               lChannel,               // channel to modify
    int32               lTrigMode,              // channel trigger mode
    int32               lTrigLevel0 = 0,        // level 0
    int32               lTrigLevel1 = 0,        // level 1
    int32               lPulsewidth = 0,        // programmable pulsewidth for channel
    bool                bTrigOut = false,       // trigger output
    bool                bSingleSrc = true);     // acts as single trigger source, all other masks cleared

// this function sets the trigger masks (bSingleSrc of other commands must be false to use this)
bool bSpcMSetupTrigMask (                       // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    uint32              dwChannelOrMask0,       // or mask of all available channels
    uint32              dwChannelOrMask1 = 0,   // ...
    uint32              dwChannelAndMask0 = 0,  // and mask of all available channels
    uint32              dwChannelAndMask1 = 0,  // ...
    uint32              dwTrigOrMask = 0,       // trigger or mask (software external, basexio)
    uint32              dwTrigAndMask = 0);     // trigger and mask (software external, basexio)



/*
**************************************************************************
bSpcMSetupInputChannel

allows all input channel related settings. if one of the setup like
termination or differential inputs is not available on the card the
setting is simply ignored
**************************************************************************
*/

bool bSpcMSetupInputChannel (                   // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    int32               lChannel,               // channel to change
    int32               lInputRange,            // input range in mV = (max-min)/2, =1000 for +/-1V range
    bool                bTerm = true,           // set input termination (50 ohm) if available
    int32               lInputOffset = 0,       // programmable input offset as listed in the manual
    bool                bDiffInput = false);    // set differential input if available

// ***********************************************************************
// suitable for M3i series with enhanced inputs

bool bSpcMSetupPathInputCh (                    // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    int32               lChannel,               // channel to change
    int32               lPath,                  // input path
    int32               lInputRange,            // input range in mV = (max-min)/2, =1000 for +/-1V range
    int32               lOffset_percent = 0,    // offset in percent, if available
    bool                bTerm = true,           // set input termination (50 ohm) if available
    bool                bACCoupling = false,    // AC coupling activated
    bool                bBWLimit = false,       // bandwidth limit activated
    bool                bDiffInput = false);    // set differential input if available


/*
**************************************************************************
bSpcMSetupAnalogOutputChannel

allows all analog output channel related settings. if one of the setup like
DoubleOut is not available on the card the setting is simply ignored
**************************************************************************
*/

bool bSpcMSetupAnalogOutputChannel (                    // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO    *pstCardInfo,                   // pointer to a filled card info structure
    int32               lChannel,                       // channel to change
    int32               lAmplitude,                     // output amplitude in mV = (max-min)/2, =1000 for +/-1V range
    int32               lOutputOffset = 0,              // programmable output offset as listed in the manual
    int32               lFilter = 0,                    // programmable output filter as listed in the manual
    int32               lStopMode = SPCM_STOPLVL_ZERO,  // defines the behavior after replay or while replay is pausing         
    bool                bDoubleOut = false,             // enables identical    output on two channels of one module (if available)
    bool                bDifferential = false);         // enables differential output on two channels of one module (if available)



/*
**************************************************************************
bSpcMSetupDigitalXXXModul

allows all input and output channel related settings for one group of 
channels. If one of the setups like the programmable output levels is not
available this setup is simply ignored
**************************************************************************
*/

bool bSpcMSetupDigitalOutput (                           // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO    *pstCardInfo,                    // pointer to a filled card info structure
    int32                lGroup,                         // module/group of channels to change
    int32                lStopMode = SPCM_STOPLVL_LOW,   // defines the behavior after replay or while replay is pausing
    int32                lLowLevel = 0,                  // low level in mV if output is programmable
    int32                lHighLevel = 3300,              // high level in mV if output levels are programmable
    bool                 bDiffMode = false);             // hardware differential mode if available

bool bSpcMSetupDigitalInput (                            // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO    *pstCardInfo,                    // pointer to a filled card info structure
    int32                lGroup,                         // module/group of channels to change
    bool                 bTerm = true);                  // input termination



/*
**************************************************************************
bSpcMSetupTimestamp

set up the timestamp mode and performs a synchronisation with refernce
clock if that mode is activated
**************************************************************************
*/

bool bSpcMSetupTimestamp (                      // returns false if error occured, otherwise true
    ST_SPCM_CARDINFO   *pstCardInfo,            // pointer to a filled card info structure
    int32               lMode,                  // mode for timestamp
    uint32              dwRefTimeoutMS);        // timeout in milli seconds for synchronisation with reference clock

#endif
