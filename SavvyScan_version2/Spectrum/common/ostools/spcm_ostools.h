/*
**************************************************************************

spcm_ostools.h                                 (c) Spectrum GmbH , 09/2005

**************************************************************************

tools functions that differ from OS to OS:
- min/max macros
- Keyboard
- Threads
- Events
- Mutex
- Page aligned memory allocation

**************************************************************************
*/


#ifndef SPCM_OSTOOLS_H
#define SPCM_OSTOOLS_H

#include "../../c_header/dlltyp.h"
#include "spcm_oswrap.h"


/* 
**************************************************************************
define our own namespace if it isn't defined yet. This is necessary to 
avoid name conflicts as this file is also used inside the driver 
**************************************************************************
*/

#ifndef SPCM_NAMESPACE
#   define SPCM_NAMESPACE spcmdrv
#endif

namespace SPCM_NAMESPACE {

/* 
**************************************************************************
misc helper functions
**************************************************************************
*/
uint32 dwGetTickCount();

void vSleep_ms (unsigned int dwMS);


/* 
**************************************************************************
Keyboard functions
**************************************************************************
*/

// check for key pressed
int bKbhit(void);

// get one character from keyboard
int cGetch();



/* 
**************************************************************************
thread functions
**************************************************************************
*/

// define the thread function 
typedef SPCM_THREAD_RETURN (SPCM_THREAD_CALLTYPE SPCM_THREADFUNCTION) (void* pvArguments);

// define the thread priority settings
enum SPCM_THREADPRIO {ePrioMin, ePrioNormal, ePrioMax};



// creates a thread with the given function name, returns false if creation failed
bool spcm_bCreateThread (SPCM_THREADFUNCTION* pfnThread, SPCM_THREAD_HANDLE* phThread, void* pvArguments);

// joins (waits for termination)
void spcm_vJoinThread  (SPCM_THREAD_HANDLE* phThread, uint32 dwTimeout_ms);

// closes the thread handle
void spcm_vCloseThread (SPCM_THREAD_HANDLE* phThread);

// suspend the thread for xx milli seconds
void spcm_vSuspendThread (uint32 dwMS);

// sets the priority of the thread
void spcm_vSetThreadPriority (SPCM_THREAD_HANDLE* phThread, SPCM_THREADPRIO ePriority);



/* 
**************************************************************************
event (linux: condition) functions
**************************************************************************
*/

// create an event 
bool spcm_bCreateEvent (SPCM_EVENT_HANDLE* phEvent);

// close an event
void spcm_vCloseEvent (SPCM_EVENT_HANDLE* phEvent);

// wait for one event with timeout, true if event was received, false if timeout occurs, timeout zero means that we wait forever
bool spcm_bWaitEventWithMutex (SPCM_EVENT_HANDLE* phEvent, SPCM_MUTEX_HANDLE* phMutex, uint32 dwTimeoutMS = 0);

// wait for one event without timeout and mutex
void spcm_vWaitEvent (SPCM_EVENT_HANDLE* phEvent);

// signal an event 
void spcm_vSignalEvent (SPCM_EVENT_HANDLE* phEvent);



/* 
**************************************************************************
mutex functions
**************************************************************************
*/

// create a mutex
bool spcm_bCreateMutex (SPCM_MUTEX_HANDLE* phMutex);

// close the mutex
void spcm_vCloseMutex (SPCM_MUTEX_HANDLE* phMutex);

// get the mutex
void spcm_vGetMutex (SPCM_MUTEX_HANDLE* phMutex);

// release the mutex
void spcm_vReleaseMutex (SPCM_MUTEX_HANDLE* phMutex);



/* 
**************************************************************************
memory allocation functions for page aligned allocation
**************************************************************************
*/

// allocate a memory region of the given size in bytes, data is page aligned
void* pvAllocMemPageAligned (uint64 qwBytes);

// free this data
void  vFreeMemPageAligned (void* pvMemory, uint64 qwBytes);



/* 
**************************************************************************
system information functions
It is not possible to compile this function with the default SDK of
Visual Studio 6.0.
**************************************************************************
*/

#if (defined (_MSC_VER) && (_MSC_VER >= 1300)) || defined (__GNUC__)
    uint64 qwGetTotalPhysicalMemory ();
    uint64 qwGetTotalVirtualMemory  ();
#endif

} // end of namespace SPCM_NAMESPACE


// our namespace that we use inside the ostools functions
using namespace SPCM_NAMESPACE;


#endif
