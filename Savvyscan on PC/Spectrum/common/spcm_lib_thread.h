/*
**************************************************************************

spcm_lib_thread.h                              (c) Spectrum GmbH , 02/2006

**************************************************************************

Offers data transfer threads for sample data, timestamp data and ABA date.
Each thread get's a bundle of parameters including a work function that
is used to do something with the data. These thread sare used in all our
FIFO examples.

**************************************************************************
*/



#ifndef SPCM_LIB_THREAD_H
#define SPCM_LIB_THREAD_H


// ----- include standard driver header from library -----
#include "../c_header/dlltyp.h"
#include "../c_header/regs.h"
#include "../c_header/spcerr.h"
#include "../c_header/spcm_drv.h"

// ----- operating system dependent functions for thead, event and mutex handling -----
#include "../common/ostools/spcm_oswrap.h"
#include "../common/ostools/spcm_ostools.h"

// ----- include of common example librarys -----
#include "../common/spcm_lib_card.h"



// our namespace that we use inside the ostools functions
using namespace SPCM_NAMESPACE;



/*
**************************************************************************
we define three working functions that the thread calls on start, on
each data update and on close. By defining these functions one can
implement nearly every functionality using this thread function.

Be aware that all these functions are called with the mutex hold. If the
function is doing more work it is requested to release the mutex inbetween

**************************************************************************
*/

struct ST_THREADDATA;
struct ST_BUFFERDATA;

// initialisation of the working routines
typedef bool (SPCM_WORK_INIT)  (void* pvWorkData, ST_BUFFERDATA* pstBufferData);

// this function does the job. It get's the current available bytes and may modify it if it doesn't use all data
typedef bool (SPCM_WORK_DO)    (void* pvWorkData, ST_BUFFERDATA* pstBufferData);

// this function closes (if file access is done) the work and cleans up if necessary
typedef void (SPCM_WORK_CLOSE) (void* pvWorkData, ST_BUFFERDATA* pstBufferData);

// simple abort function to be implemented by user program
typedef bool (SPCM_ABORTCHECK) (void* pvWorkData, ST_BUFFERDATA* pstBufferData);



// simple abort function that waits returns true if escape is pressed
bool bKeyAbortCheck (void* pvWorkData, ST_BUFFERDATA* pstBufferData);



/*
**************************************************************************
buffer data: structure to handover current buffer status and available
bytes. The structure is modified several times:

WORK_INIT: 
  get's the general setup (pstCard object)
  defines BufLen and Notify for up to three transfer objects
  may set an error code

WORK_DO
  get's the current pointers to the data buffers and the currently available bytes
  modifies the currently available bytes if it doesn't uses all of them
  may set an error code

**************************************************************************
*/

struct ST_BUFFERDATA
    {

    // general setup and status
    ST_SPCM_CARDINFO*   pstCard;                // pointer to the shared card info structure
    uint32              dwError;                // error code 
    int32               lTimeout;               // timeout in ms for the wait commands

    // starting flags for the different transfers (all together are started)
    bool                bStartCard;             // starts the card
    bool                bStartData;             // starts the data transfer
    bool                bStartTimestamp;        // starts the timestamp transfer
    bool                bStartABA;              // starts the ABA transfer
    bool                bStartExtraDMA;         // this thread starts the extra DMA

    // these settings are to be set by WorkSetup and are used for definitions
    uint32              dwDataBufLen;           // complete buffer length 
    uint32              dwDataNotify;           // notify size 

    uint32              dwTSBufLen;             // buffer length of the timestamp
    uint32              dwTSNotify;             // notify size of timestamp

    uint32              dwABABufLen;            // ABA buffer length
    uint32              dwABANotify;            // notify size of ABA buffer

    // these settings are maintained by the main loop (thread) and are update by the working function
    void*               pvDataCurrentBuf;       // a pointer to the current data buffer starting address
    int64               llDataAvailBytes;       // current data length that is available
    uint64              qwDataTransferred;      // complete transferred amount of data

    void*               pvTSCurrentBuf;         // a pointer to the current timestamp buffer starting address
    int64               llTSAvailBytes;         // current timestamp length that is available
    uint64              qwTSTransferred;        // complete transferred amount of timestamp data

    void*               pvABACurrentBuf;        // a pointer to the current ABA buffe starting address
    int64               llABAAvailBytes;        // current ABA data length that is available
    uint64              qwABATransferred;       // complete transferred amount of ABA data
    };



/*
**************************************************************************
The thread data structure contains all exchange data between main program
and the transfer thread. It also contains the mutexes and events to 
allow synchronisation
**************************************************************************
*/
#define EVENTFLAG_INITDONE  0x0000001
#define EVENTFLAG_START     0x0000002
#define EVENTFLAG_STARTDONE 0x0000004
#define EVENTFLAG_TERMDONE  0x0000008


struct ST_THREADDATA 
    {
    // events from main to thread
    SPCM_EVENT_HANDLE   hEventStart;            // starts the acquistion

    // events from thread to main
    SPCM_EVENT_HANDLE   hEventInitDone;         // initialisation of thread is done
    SPCM_EVENT_HANDLE   hEventStartDone;        // the start has been done, tranfer is running now
    SPCM_EVENT_HANDLE   hEventUpdate;           // data is updated
    SPCM_EVENT_HANDLE   hEventTermDone;         // termination is done
    bool                bSharedUpdate;          // the update event is shared with other threads
    SPCM_EVENT_HANDLE*  phSharedUpdate;         // either the shared update event or my own one if I'm the one to share

    uint32              dwEventFlags;           // event flags to avoid Linux pthread problems

    SPCM_MUTEX_HANDLE   hMutexSharedAccess;     // controls the multiple exclusion access to the buffer data section
    ST_BUFFERDATA*      pstBufferData;          // a buffer data structure to fill and forward to working function
    void*               pvWorkData;             // the private working data

    SPCM_WORK_INIT*     pfn_bWorkInit;          // initialisation function
    SPCM_WORK_DO*       pfn_bWorkDo;            // the working function
    SPCM_WORK_CLOSE*    pfn_vWorkClose;         // closes the working routine
    };



/*
**************************************************************************
vDoThread MainLoop: 

contains the default main thread loop. This one can be used to control one 
common thread that cares for all transfers or separate threads for data, 
timestamp and ABA data.

Select the behaviour with the bSeparateThread flag. If separate threads
are used it is necessary to give init, worker and close function callbacks
for each transfer thread. 

**************************************************************************
*/

void vDoThreadMainLoop (
    ST_BUFFERDATA*      pstBufferData,          // filled card info structure
    void*               pvWorkData,             // working data structure that is passed to the calback functions

    SPCM_WORK_INIT*     pfn_bWorkInit,          // callback function for pre-run setup
    SPCM_WORK_DO*       pfn_bWorkDo,            // callback function that is called if new data is available
    SPCM_WORK_CLOSE*    pfn_vWorkClose = NULL,  // callback function at the end of run

    SPCM_ABORTCHECK*    pfn_bAbortCheck = NULL, // callback function for abort checking on each update loop
    bool                bSeparateThread = false,// use separate threads for data, timestamp and ABA

    SPCM_WORK_INIT*     pfn_bTSInit = NULL,     // callback function for timestamo pre-run setup
    SPCM_WORK_DO*       pfn_bTSDo = NULL,       // callback function that is called if new timestamps are available
    SPCM_WORK_CLOSE*    pfn_vTSClose = NULL,    // callback function at the end of timestamp run

    SPCM_WORK_INIT*     pfn_bABAInit = NULL,    // callback function for ABA pre-run setup
    SPCM_WORK_DO*       pfn_bABADo = NULL,      // callback function that is called if new ABA data is available
    SPCM_WORK_CLOSE*    pfn_vABAClose = NULL);  // callback function at the end of ABA run




/*
**************************************************************************
vDoMainLoop: 

contains the default main loop in a non-threaded version. This one can be 
used to control data transfer, timestamp transfer and ABA transfer at the 
same time without using threads.
**************************************************************************
*/

void vDoMainLoop (
    ST_BUFFERDATA*      pstBufferData,          // filled card info structure
    void*               pvWorkData,             // working data structure that is passed to the calback functions

    SPCM_WORK_INIT*     pfn_bWorkInit,          // callback function for pre-run setup
    SPCM_WORK_DO*       pfn_bWorkDo,            // callback function that is called if new data is available
    SPCM_WORK_CLOSE*    pfn_vWorkClose = NULL,  // callback function at the end of run

    SPCM_ABORTCHECK*    pfn_bAbortCheck = NULL);// callback function for abort checking on each update loop

#endif
