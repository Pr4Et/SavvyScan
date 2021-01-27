#ifndef SPCM_OSWRAP_H
#define SPCM_OSWRAP_H



/*
**************************************************************************

spcm_oswrap.h                                  (c) Spectrum GmbH , 08/2005

**************************************************************************

Contains some wrapper functions, OS specific defines and OS specific 
includes to make the source compilable independant of the operating
system

**************************************************************************
*/



/*
**************************************************************************
Linux
**************************************************************************
*/

#if defined (_LINUX) || defined (_QNX)
#   include <unistd.h>
#   include <pthread.h>

// ----- Linux specific defines -----
#   define NULL_HANDLE 0
#   define _stdcall
#   define SPCM_THREAD_RETURN void*
#   define SPCM_THREAD_CALLTYPE

// ----- handles -----
#   define SPCM_THREAD_HANDLE   pthread_t
#   define SPCM_EVENT_HANDLE    pthread_cond_t
#   define SPCM_MUTEX_HANDLE    pthread_mutex_t



/*
**************************************************************************
Windows
**************************************************************************
*/

#else

// ----- Windows specific includes -----
#   include <wtypes.h>

// ----- Windows specific defines
#   define NULL_HANDLE NULL
#   define SPCM_THREAD_RETURN uint32
#   define SPCM_THREAD_CALLTYPE _stdcall

// ----- handles -----
#   define SPCM_THREAD_HANDLE   HANDLE
#   define SPCM_EVENT_HANDLE    HANDLE
#   define SPCM_MUTEX_HANDLE    CRITICAL_SECTION

#   if (_MSC_VER < 1900) // 1900 = VS2015
#       define snprintf _snprintf
#   endif
#endif

#endif //#ifndef SPCM_OSWRAP_H
