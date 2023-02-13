/*
**************************************************************************

spcm_ostools.cpp                                         (c) Spectrum GmbH

**************************************************************************

toosl functions that differ from OS to OS:
- Threads
- Events

**************************************************************************
*/

#include "../../c_header/dlltyp.h"

#include "spcm_oswrap.h"
#include "spcm_ostools.h"

#include <conio.h>


// we're using this namespace in our program
namespace SPCM_NAMESPACE {


uint32 dwGetTickCount ()
    {
    return ::GetTickCount ();
    }

void vSleep_ms (unsigned int dwMS)
    {
    ::Sleep (dwMS);
    }

/* 
**************************************************************************
Keyboard functions (wrapped for combined linux/windows use)
**************************************************************************
*/

int bKbhit(void)
    {
    return _kbhit();
    }

// ***********************************************************************

int cGetch()
    {
    return _getch();
    }



/* 
**************************************************************************
Thread functions
**************************************************************************
*/

bool spcm_bCreateThread (SPCM_THREADFUNCTION* pfnThread, SPCM_THREAD_HANDLE* phThread, void* pvArguments)
    {
    uint32 dwThreadId;

    (*phThread) = CreateThread (NULL, 0, pfnThread, pvArguments, 0, &dwThreadId);

    return ((*phThread) != NULL);
    }

// ***********************************************************************

void spcm_vCloseThread (SPCM_THREAD_HANDLE* phThread)
    {
    if (*phThread)
        CloseHandle (*phThread);
    }

// ***********************************************************************

void spcm_vJoinThread  (SPCM_THREAD_HANDLE* phThread, uint32 dwTimeout_ms)
    {
    WaitForSingleObject (*phThread, dwTimeout_ms ? dwTimeout_ms : INFINITE); 
    }

// ***********************************************************************

void spcm_vSetThreadPriority (SPCM_THREAD_HANDLE* phThread, SPCM_THREADPRIO ePriority)
    {
    switch (ePriority)
        {
        case ePrioMin:
            SetThreadPriority (*phThread, THREAD_PRIORITY_BELOW_NORMAL);
            break;
        case ePrioNormal:
            SetThreadPriority (*phThread, THREAD_PRIORITY_NORMAL);
            break;
        case ePrioMax:
            SetThreadPriority (*phThread, THREAD_PRIORITY_ABOVE_NORMAL);
            break;
        }
    }

// ***********************************************************************

void spcm_vSuspendThread (uint32 dwMS)
    {
    Sleep (dwMS);
    }


/* 
**************************************************************************
Event functions
**************************************************************************
*/

bool spcm_bCreateEvent (SPCM_EVENT_HANDLE* phEvent)
    {
    (*phEvent) = CreateEvent (NULL, false, false, NULL);
    return ((*phEvent) != NULL);
    }

// ***********************************************************************

void spcm_vCloseEvent (SPCM_EVENT_HANDLE* phEvent)
    {
    if (*phEvent)
        CloseHandle (*phEvent);
    }

// ***********************************************************************

bool spcm_bWaitEventWithMutex (SPCM_EVENT_HANDLE* phEvent, SPCM_MUTEX_HANDLE* phMutex, uint32 dwTimeoutMS)
    {
    uint32 dwReturn;

    // release the mutex, wait for the event and get the mutex again
    LeaveCriticalSection (phMutex);
    dwReturn = WaitForSingleObject ((*phEvent), dwTimeoutMS ? dwTimeoutMS : INFINITE);
    EnterCriticalSection (phMutex);

    return (dwReturn != WAIT_TIMEOUT);
    }

// ***********************************************************************
 
void spcm_vWaitEvent (SPCM_EVENT_HANDLE* phEvent)
    {
    WaitForSingleObject ((*phEvent), INFINITE);
    }

// ***********************************************************************

void spcm_vSignalEvent (SPCM_EVENT_HANDLE* phEvent)
    {
    SetEvent (*phEvent);
    }

                 

/* 
**************************************************************************
Mutex functions (we use CriticalSection here to speed it up!)
**************************************************************************
*/

bool spcm_bCreateMutex (SPCM_MUTEX_HANDLE* phMutex)
    {
    InitializeCriticalSection (phMutex);
    return true;
    }

// ***********************************************************************

void spcm_vCloseMutex (SPCM_MUTEX_HANDLE* phMutex)
    {
    if (phMutex)
        DeleteCriticalSection (phMutex);
    }

// ***********************************************************************

void spcm_vGetMutex (SPCM_MUTEX_HANDLE* phMutex)
    {
    EnterCriticalSection (phMutex);
    }

// ***********************************************************************

void spcm_vReleaseMutex (SPCM_MUTEX_HANDLE* phMutex)
    {
    LeaveCriticalSection (phMutex);
    }



/* 
**************************************************************************
Data allocation functions
**************************************************************************
*/

void* pvAllocMemPageAligned (uint64 qwBytes)
    {
    // for unknown reasons VirtualAlloc/VirtualFree leaks memory if qwBytes < 4096 (page size)
    // therefore use _aligned_malloc () to get small amounts of page aligned memory
    if (qwBytes >= 4096)
        return VirtualAlloc (NULL, (size_t) qwBytes, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
    else
        {
        void* pvMem = _aligned_malloc ((size_t) qwBytes, 4096);
        if (pvMem == NULL)
            return NULL;
        memset (pvMem, 0, (size_t) qwBytes);
        return pvMem;
        }
    }

// ***********************************************************************

void vFreeMemPageAligned (void* pvMemory, uint64 qwBytes)
    {
    // for unknown reasons VirtualAlloc/VirtualFree leaks memory if qwBytes < 4096 (page size)
    // therefore use _aligned_malloc () to get small amounts of page aligned memory    
    if (qwBytes >= 4096)
        VirtualFree (pvMemory, 0, MEM_RELEASE);
    else
        _aligned_free (pvMemory);
    }

/* 
**************************************************************************
System information functions
It is not possible to compile this function with the default SDK of
Visual Studio 6.0.
**************************************************************************
*/

#if defined (_MSC_VER) && (_MSC_VER >= 1300)
    uint64 qwGetTotalPhysicalMemory ()
        {
        MEMORYSTATUSEX stMemoryStatus;
        stMemoryStatus.dwLength = sizeof (stMemoryStatus);
        GlobalMemoryStatusEx (&stMemoryStatus);
        return stMemoryStatus.ullTotalPhys;
        }

    uint64 qwGetTotalVirtualMemory ()
        {
        MEMORYSTATUSEX stMemoryStatus;
        stMemoryStatus.dwLength = sizeof (stMemoryStatus);
        GlobalMemoryStatusEx (&stMemoryStatus);
        return stMemoryStatus.ullTotalVirtual;
        }
#endif

} // end of namespace SPCM_NAMESPACE
