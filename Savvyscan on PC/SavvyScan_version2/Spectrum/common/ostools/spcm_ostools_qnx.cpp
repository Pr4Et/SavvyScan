/*
**************************************************************************

spcm_ostools_qnx.cpp                                     (c) Spectrum GmbH

**************************************************************************

toosl functions that differ from OS to OS:
- Threads
- Events

**************************************************************************
*/

#include "../../c_header/dlltyp.h"

#include "spcm_oswrap.h"
#include "spcm_ostools.h"

#include "../kernel/spcm_kernel.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/time.h>
#include <sys/syspage.h>

#include <pthread.h>
#include <cstring>
#include <cstdio>
#include <termios.h> 
#include <ctime>
#include <unistd.h>

#include <cerrno>

// we're using this namespace in our program
namespace SPCM_NAMESPACE {

using namespace std;

/*
**************************************************************************
GetTickCount: fake windows GetTickCount function for QNX
**************************************************************************
*/

uint32 dwGetTickCount ()
    {
    struct timespec ts;
    clock_gettime (0, &ts);
    return (uint32) 1000 * ts.tv_sec + ts.tv_nsec / 1000000;
    }


/*
**************************************************************************
Sleep: suspend the thread for the number of ms
       On some systems is the usleep parameter limited to 1000000.
**************************************************************************
*/

void Sleep (uint32 dwMS)
    {
    uint32 dwRest_us, dwActual_us;

    dwRest_us = dwMS * 1000;

    // do/while construct to call usleep with zero delay too
    do
        {
        dwActual_us = dwRest_us;
        if (dwActual_us > 500000)
            dwActual_us = 500000;

        dwRest_us -= dwActual_us;

        usleep (dwActual_us);
        }
    while (dwRest_us);
    }


/* 
**************************************************************************
Keyboard workaround functions
**************************************************************************
*/


 int cGetch ()
    { 
    static int ch = -1, fd = 0;
    struct termios stTerm, stOldTerm;

    fd = fileno(stdin);
    tcgetattr(fd, &stOldTerm);
    stTerm = stOldTerm;
    stTerm.c_lflag &= ~(ICANON|ECHO);
    tcsetattr(fd, TCSANOW, &stTerm);
    ch = getchar();
    tcsetattr(fd, TCSANOW, &stOldTerm);

    return ch;
    }

// ***********************************************************************

int bKbhit ()
    {
    struct termios stTerm, stOldTerm;
    int fd = 0;
    int c = 0;

    tcgetattr(fd, &stOldTerm);
    memcpy(&stTerm, &stOldTerm, sizeof (stOldTerm));
    stTerm.c_lflag = stTerm.c_lflag & (!ICANON);
    stTerm.c_cc[VMIN] = 0;
    stTerm.c_cc[VTIME] = 1;
    tcsetattr(fd, TCSANOW, &stTerm);
    c = getchar();
    tcsetattr(fd, TCSANOW, &stOldTerm);
    if (c != -1)
        ungetc (c, stdin);

    return ((c != -1) ? 1 : 0);
    }



/*
**************************************************************************
memory allocation (page aligned)
**************************************************************************
*/

void* pvAllocMemPageAligned (uint64 qwBytes)
    {
    void* pvTmp;
    int fd = open ("/dev/zero", O_RDONLY);
    pvTmp = (void*) mmap (NULL, qwBytes, PROT_READ | PROT_WRITE, MAP_PRIVATE, fd, 0);

    // set everything to zero to get memory allocated in physical mem
    if (pvTmp != MAP_FAILED)
        memset (pvTmp, 0, qwBytes);
    else
        pvTmp = NULL;

    close (fd);
    return (pvTmp);
    }

// ***********************************************************************

void vFreeMemPageAligned (void* pvAdr, uint64 qwBytes)
    {
    munmap (pvAdr, qwBytes);
    }



/* 
**************************************************************************
Thread functions
**************************************************************************
*/

bool spcm_bCreateThread (SPCM_THREADFUNCTION* pfnThread, SPCM_THREAD_HANDLE* phThread, void* pvArguments)
    {
    return (pthread_create (phThread, NULL, pfnThread, pvArguments) == EOK);
    }

// ***********************************************************************

void spcm_vCloseThread (SPCM_THREAD_HANDLE* phThread)
    {
    pthread_detach (*phThread);
    }

// ***********************************************************************
 
void spcm_vSetThreadPriority (SPCM_THREAD_HANDLE* phThread, SPCM_THREADPRIO ePriority)
    {
    struct sched_param stSchedParams;
    int32 lPolicy;
    pthread_getschedparam (*phThread, &lPolicy, &stSchedParams);
    if (lPolicy == SCHED_OTHER)
        {
        switch (ePriority)
            {
            case ePrioMin:
                stSchedParams.sched_priority = 5; // we are nice
                break;
            case ePrioNormal:
                stSchedParams.sched_priority = 0; // we are normally nice
                break;
            case ePrioMax:
                stSchedParams.sched_priority = -5; // we are not so nice
                break;
            }
        }
    else
        {
        int lPrioMin = sched_get_priority_min (lPolicy);
        int lPrioMax = sched_get_priority_max (lPolicy);
        switch (ePriority)
            {
            case ePrioMin:
                stSchedParams.sched_priority = (lPrioMax - lPrioMin) / 4;
                break;
            case ePrioNormal:
                stSchedParams.sched_priority = (lPrioMax - lPrioMin) / 2;
                break;
            case ePrioMax:
                stSchedParams.sched_priority = 3 * (lPrioMax - lPrioMin) / 4;
                break;
            }
        }

    pthread_setschedparam (*phThread, lPolicy, &stSchedParams);
    }

// ***********************************************************************

void spcm_vJoinThread  (SPCM_THREAD_HANDLE* phThread, uint32 /*dwTimeout_ms*/)
    {
    pthread_join (*phThread, NULL);
    }

// ***********************************************************************

void spcm_vSuspendThread (uint32)
    {
    //pthread_yield ();
    sched_yield (); // TODO: guter Ersatz?
    }



/* 
**************************************************************************
Event functions
**************************************************************************
*/

bool spcm_bCreateEvent (SPCM_EVENT_HANDLE* phEvent)
    {
    // Set up the condvar attributes to use CLOCK_MONOTONIC
    pthread_condattr_t attr;
    pthread_condattr_init( &attr);
    pthread_condattr_setclock( &attr, CLOCK_MONOTONIC);

    return (pthread_cond_init (phEvent, &attr) == EOK);
    }

// ***********************************************************************

void spcm_vCloseEvent (SPCM_EVENT_HANDLE* phEvent)
    {
    pthread_cond_destroy (phEvent);
    }

// ***********************************************************************

bool spcm_bWaitEventWithMutex (SPCM_EVENT_HANDLE* phEvent, SPCM_MUTEX_HANDLE* phMutex, uint32 dwTimeoutMS)
    {
    bool bRet;

    if (dwTimeoutMS)
        {
        struct timespec stTimeSpec;
        clock_gettime (CLOCK_MONOTONIC, &stTimeSpec);

        // add my wait time
        stTimeSpec.tv_sec  +=  ((stTimeSpec.tv_sec) + dwTimeoutMS) / 1000;
        stTimeSpec.tv_nsec  = (((stTimeSpec.tv_nsec / 1000 / 1000) + dwTimeoutMS) % 1000) * 1000 * 1000;

        bRet = (pthread_cond_timedwait (phEvent, phMutex, &stTimeSpec) == 0);
        }

    // no timeout specified, we wait forever
    else
        bRet = (pthread_cond_wait (phEvent, phMutex) == 0);

    return bRet;
    }

// ***********************************************************************
 
void spcm_vWaitEvent (SPCM_EVENT_HANDLE* phEvent)
    {
    pthread_mutex_t hTmpMutex;

    pthread_mutex_init (&hTmpMutex, NULL);
    pthread_mutex_lock (&hTmpMutex);
    pthread_cond_wait (phEvent, &hTmpMutex);
    pthread_mutex_destroy (&hTmpMutex);
    }

// ***********************************************************************

void spcm_vSignalEvent (SPCM_EVENT_HANDLE* phEvent)
    {
    pthread_cond_signal (phEvent);
    }



/* 
**************************************************************************
mutex functions
**************************************************************************
*/

bool spcm_bCreateMutex (SPCM_MUTEX_HANDLE* phMutex)
    {
    return (pthread_mutex_init (phMutex, NULL) == EOK);
    }

// ***********************************************************************

void spcm_vCloseMutex (SPCM_MUTEX_HANDLE* phMutex)
    {
    pthread_mutex_destroy (phMutex);
    }

// ***********************************************************************

void spcm_vGetMutex (SPCM_MUTEX_HANDLE* phMutex)
    {
    pthread_mutex_lock (phMutex);
    }

// ***********************************************************************

void spcm_vReleaseMutex (SPCM_MUTEX_HANDLE* phMutex)
    {
    pthread_mutex_unlock (phMutex);
    }



/* 
**************************************************************************
system information functions
**************************************************************************
*/

uint64 qwGetTotalPhysicalMemory ()
    {
    struct asinfo_entry* pstEntries = SYSPAGE_ENTRY(asinfo);
    size_t count = SYSPAGE_ENTRY_SIZE(asinfo) / sizeof(struct asinfo_entry);
    char* szStrings = SYSPAGE_ENTRY(strings)->data;

    uint64 qwTotal = 0;
    for (size_t i = 0; i < count; i++)
        {
        struct asinfo_entry* pstEntry = &pstEntries[i];
        if (strcmp (szStrings + pstEntry->name, "ram") == 0)
            {
            qwTotal += pstEntry->end - pstEntry->start + 1;
            }
        }

    return qwTotal;
    }

uint64 qwGetTotalVirtualMemory ()
    {
#ifdef _LINUX64
    return (uint64)8 * 1024 * 1024 * 1024 * 1024; // 8TB, value taken from 64-Bit Windows
#else // 32 bit
    return (uint64)3 * 1024 * 1024 * 1024;        // 3GB for user-space in QNX systems
#endif
    }
} // end of SPCM_NAMESPACE
