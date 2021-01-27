#ifndef DLLTYP_H
#define DLLTYP_H



/*
**************************************************************************

dlltyp.h                                           (c) Spectrum GmbH, 2006

**************************************************************************

definitions common for all Spectrum drivers and card types. This header 
tries to examine the type of compiler and then defines common data types 
that have the same length under all compilers and operating systems. 

This header should be the first header to include in all C/C++ projects

**************************************************************************

Please do not change this file as it's continuously updated with new 
driver versions. If you need other settings for your special compiler 
type please add these settings in an extra header file
  
**************************************************************************
*/

// ----- some global definitions for the drivers -----
#define MAXBRD          64
#define SPCM_MAXDEMO    100



/*
**************************************************************************
this part now tries to examine the compiler type and sets one of the
defines that is used later for the type definition
**************************************************************************
*/

// ----- Borland C+ compiler. If the version is > 0x410 it's the C++ Builder and the same types can be used as for Visual C++ -----
#ifdef __BORLANDC__
#   if (__BCPLUSPLUS__>0x410)
#       define VC_WINNT
#   elif defined(_Windows)
#       define BC_WIN31
#   else
#       define BC_DOS
#   endif
#endif

// ----- Microsoft Visual C++ compiler, either as std_call or as c_call -----
#ifdef _WINSTDCALL
#   ifdef _MSC_VER
#       ifdef _WINNT
#           define VC_STDCALLNT
#       else
#           define VC_STDCALL95
#       endif
#   endif
#elif defined(_MSC_VER)
#   ifdef _WIN32
#       ifdef _WINNT
#           define VC_WINNT
#       else
#           define VC_WIN95
#       endif
#   else
#       define VC_WIN31
#   endif
#endif
#if defined (_WIN32) || defined (_WIN64)
#   define _WIN // analog to _LINUX
#endif

// ----- GCC in various environments
#if defined (__GNUC__)
#   if defined (__MINGW32__)
#       define _GCCWIN
#   elif defined (__CYGWIN__)
#       define _GCCWIN
#   elif defined (__QNX__)
#       define _QNX
#   elif !defined(__vxworks)
#       define _LINUX
#   endif
#endif

// VxWorks
#ifdef __vxworks
#   define _VXWORKS
#endif


// ----- LabWindows/CVI
#ifdef _CVI_
#   define _LWCVI
#endif

// ----- 64 Bit Linux (in addition to _LINUX) -----
#if defined (__x86_64__) || defined (__amd64__)
#   if defined (_LINUX)
#       define _LINUX64
#   endif
#endif



/*
**************************************************************************
now we hopefully know the compiler type and define all the types matching
this compiler and the platform
**************************************************************************
*/

// ----- Borland C++ for DOS (only used with older drivers) -----
#ifdef BC_DOS
#   define int16   int
#   define uint16  unsigned int
#   define int8    char
#   define uint8   unsigned char
#   define int32   long int
#   define uint32  unsigned long int
#   define dataptr void huge*
#   define ptr8    char huge*
#   define uptr8   unsigned char huge*
#   define ptr16   int huge*
#   define uptr16  unsigned int huge*
#   define ptr32   long int huge*
#   define uptr32  unsigned long int huge*
#   define bool    int8
#   define true    1
#   define false   0
#   define EXP     extern "C" _export int16
#   define EXPC    extern  _export int16
#   define IMP     extern "C" _import int16
#   define HEAD    extern "C" int16
#endif

// ----- Borland C++ for Windows 3.1/3.11 (only used with older drivers) -----
#ifdef BC_WIN31
#   define int16   int
#   define uint16  unsigned int
#   define int8    char
#   define uint8   unsigned char
#   define int32   long int
#   define uint32  unsigned long int
#   define dataptr void huge*
#   define ptr8    char huge*
#   define uptr8   unsigned char huge*
#   define ptr16   int huge*
#   define uptr16  unsigned int huge*
#   define ptr32   long int huge*
#   define uptr32  unsigned long int huge*
#   ifdef _EasyWin
#       define EXP     extern "C" _export int16
#       define IMP     extern "C" _import int16
#       define HEAD    extern "C" int16
#   else
#       define EXP     extern "C" _export int16 FAR PASCAL
#       define EXPC    extern     _export int16 FAR PASCAL
#       define IMP     extern "C" _import int16 FAR PASCAL
#       define HEAD    extern "C" int16 FAR PASCAL
#   endif
#endif

// ----- Visual C++ for Windows 3.1/3.11 (only used with older drivers) -----
#ifdef VC_WIN31
#   define int8    char
#   define uint8   unsigned char
#   define int16   short int
#   define uint16  unsigned short int
#   define int32   long int
#   define uint32  unsigned long int
#   define dataptr void huge*
#   define ptr8    char*
#   define uptr8   unsigned char*
#   define ptr16   short int*
#   define uptr16  unsigned short int*
#   define ptr32   long int*
#   define uptr32  unsigned long int*
#   define EXP     extern "C" __declspec (dllexport) int16
#   define IMP     extern "C" __declspec (dllimport) int16
#   define HEAD    extern "C" __declspec (dllexport) int16
#endif

// ----- Visual C++ / Borland C++ Builder for Windows 32 bit starting with Windows 95 -----
#if defined(VC_WIN95) || defined(VC_WINNT)
#   ifndef NO_WTYPES_IN_DLLTYP
#       include <wtypes.h>
#   endif
#   define int8    char
#   define uint8   unsigned char
#   define int16   short int
#   define uint16  unsigned short int
#   define int32   long int
#   define uint32  unsigned long int
#   define int64   __int64
#   define uint64  unsigned int64
#   define dataptr void*
#   define ptr8    char*
#   define uptr8   unsigned char*
#   define ptr16   short int*
#   define uptr16  unsigned short int*
#   define ptr32   long int*
#   define uptr32  unsigned long int*
#   ifndef __cplusplus
#       define bool int8
#       define true 1
#       define false 0
#   endif
#   define drv_handle  void*

#   ifdef __cplusplus
#       define EXP     extern "C" __declspec (dllexport) int16
#       define EXPC    extern     __declspec (dllexport) int16
#       define IMP     extern "C" __declspec (dllimport) int16
#       define HEAD    extern "C" __declspec (dllexport) int16

#       define SPCM_EXPORT extern "C" __declspec (dllexport) 
#       define SPCM_IMPORT extern "C" __declspec (dllimport) 
#   else
#       define EXP     extern __declspec (dllexport) int16
#       define EXPC    extern     __declspec (dllexport) int16
#       define IMP     extern __declspec (dllimport) int16
#       define HEAD    extern __declspec (dllexport) int16

#       define SPCM_EXPORT extern __declspec (dllexport) 
#       define SPCM_IMPORT extern __declspec (dllimport) 
#   endif
#endif

// ----- Visual C++ using standard calls for all windows 32 bit platforms -----
#if defined(VC_STDCALL95) || defined(VC_STDCALLNT)
#   ifndef NO_WTYPES_IN_DLLTYP
#       include <wtypes.h>
#   endif
#   define int8    char
#   define uint8   unsigned char
#   define int16   short int
#   define uint16  unsigned short int
#   define int32   long int
#   define uint32  unsigned long int
#   define int64   __int64
#   define uint64  unsigned __int64
#   define dataptr void*
#   define ptr8    char*
#   define uptr8   unsigned char*
#   define ptr16   short int*
#   define uptr16  unsigned short int*
#   define ptr32   long int*
#   define uptr32  unsigned long int*
#   ifndef __cplusplus
#       define bool int8
#       define true 1
#       define false 0
#   endif
#   define drv_handle  void*
#   define EXP     extern "C" __declspec (dllexport) int16 _stdcall
#   define EXPC    extern     __declspec (dllexport) int16 _stdcall
#   define IMP     extern "C" __declspec (dllimport) int16 _stdcall
#   define HEAD    extern "C" __declspec (dllexport) int16 _stdcall

#   define SPCM_EXPORT extern "C" __declspec (dllexport) 
#   define SPCM_IMPORT extern "C" __declspec (dllimport) 
#endif

// ----- Linux -----
#if defined (_LINUX) || defined (_QNX)
#   define int8 char
#   define int16 short int
#   define int32 int
#   define int64 long long
#   define uint8 unsigned char
#   define uint16 unsigned short int
#   define uint32 unsigned int
#   define uint64 unsigned long long
#   define dataptr void *
#   define ptr8 int8*
#   define ptr16 int16*
#   define ptr32 int32*
#   define uptr8 uint8*
#   define uptr16 uint16*
#   define uptr32 uint32*
#   if !defined(bool) && !defined(__cplusplus)
#       define bool int8
#       define true 1
#       define false 0
#   endif
#   define drv_handle void*
#   define EXPC int16
#   define HEAD int16
#   define SPEC_IOC_MAGIC 's'
#   define SPEC_IOC_MAXNR 6
#   define GETPARAM  _IOR(SPEC_IOC_MAGIC,1,int32[2])
#   define SETPARAM  _IOW(SPEC_IOC_MAGIC,2,int32[2])
#   define GETCH     _IOR(SPEC_IOC_MAGIC,3,int32[1])
#   define SETCH     _IOW(SPEC_IOC_MAGIC,4,int32[1])
    typedef struct {int32 lReg; void* pvAdr;} _SETGETADR ;
#   define SETADR    _IOW(SPEC_IOC_MAGIC,5,_SETGETADR)
#   define GETADR    _IOR(SPEC_IOC_MAGIC,6,_SETGETADR)

#   ifdef __cplusplus
#       define SPCM_IMPORT extern "C"
#       if __GNUC__ >= 4
#           define SPCM_EXPORT extern "C" __attribute__ ((visibility ("default")))
#       else
#           define SPCM_EXPORT extern "C"
#       endif
#   else
#       define SPCM_IMPORT
#       define SPCM_EXPORT extern "C"
#   endif
#   define _stdcall
#endif

// ----- LabWindows/CVI -----
#if defined(_LWCVI)
#   define int8    char
#   define uint8   unsigned char
#   define int16   short int
#   define uint16  unsigned short int
#   define int32   long int
#   define uint32  unsigned long int
#   define int64   __int64
#   define uint64  unsigned int64
#   define dataptr void*
#   define ptr8    char*
#   define uptr8   unsigned char*
#   define ptr16   short int*
#   define uptr16  unsigned short int*
#   define ptr32   long int*
#   define uptr32  unsigned long int*
#   define drv_handle  void*
#   define bool    int8
#   define true    1
#   define false   0
#   define SPCM_EXPORT extern "C" __declspec (dllexport)
#   define SPCM_IMPORT 
#    define _stdcall __stdcall
#endif

// ----- Gnu C Windows -----
#if defined (_GCCWIN)
#   define int8 char
#   define int16 short int
#   define int32 int
#   define int64 long long
#   define uint8 unsigned char
#   define uint16 unsigned short int
#   define uint32 unsigned int
#   define uint64 unsigned long long
#   define dataptr void *
#   define ptr8 int8*
#   define ptr16 int16*
#   define ptr32 int32*
#   define uptr8 uint8*
#   define uptr16 uint16*
#   define uptr32 uint32*
#   if !defined(bool) && !defined(__cplusplus)
#       define bool int8
#       define true 1
#       define false 0
#   endif
#   define drv_handle void*
#   define EXPC int16
#   define HEAD int16
#   ifdef __cplusplus
#       define SPCM_IMPORT extern "C"
#   else
#       define SPCM_IMPORT
#   endif
#endif

// --- define data structure for segment statistic mode
typedef struct
    {
    int64  llAvrg:              64;
    int16  nMin:                16;
    int16  nMax:                16;
    uint32 dwMinPos:            32;
    uint32 dwMaxPos:            32;
    uint32 _Unused:             32;
    uint64 qw_Timestamp:        64;
    } SPCM_SEGSTAT_STRUCT_CHx;

#endif
