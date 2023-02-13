/*
**************************************************************************

spcm_load_dll.cpp                                       (c) Spectrum GmbH

**************************************************************************

Example for all SpcMDrv based products.

Information about the different products and their drivers can be found
online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/platform-driver-and-series-differences

This example shows how to load library functions dynamically from a
windows dll or linux shared object (.so),
e.g. if not having a matching lib-file for the compiler.

Feel free to use this source for own projects and modify it in any kind.

Documentation for the API as well as a detailed description of the hardware
can be found in the manual for each device which can be found on our website:
https://www.spectrum-instrumentation.com/en/downloads

Further information can be found online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/knowledge-base-overview

**************************************************************************
*/


// ----- include standard driver header from library -----
#include "../c_header/dlltyp.h"
#include "../c_header/regs.h"
#include "../c_header/spcerr.h"
#include "../c_header/spcm_drv.h"

// ----- standard c include files -----
#include <stdio.h>
#include <windows.h>



// ----- include the driver function definitions -----
#include "spcm_drv_def.h"


// ----- global definitions of the function pointers -----
SPCM_HOPEN*                 pfn_spcm_hOpen;
SPCM_VCLOSE*                pfn_spcm_vClose;
SPCM_DWSETPARAM_I32*        pfn_spcm_dwSetParam_i32;
SPCM_DWSETPARAM_I64M*       pfn_spcm_dwSetParam_i64m;
SPCM_DWSETPARAM_I64*        pfn_spcm_dwSetParam_i64;
SPCM_DWGETPARAM_I32*        pfn_spcm_dwGetParam_i32;
SPCM_DWGETPARAM_I64M*       pfn_spcm_dwGetParam_i64m;
SPCM_DWGETPARAM_I64*        pfn_spcm_dwGetParam_i64;
SPCM_DWDEFTRANSFER_I64M*    pfn_spcm_dwDefTransfer_i64m;
SPCM_DWDEFTRANSFER_I64*     pfn_spcm_dwDefTransfer_i64;
SPCM_DWINVALIDATEBUF*       pfn_spcm_dwInvalidateBuf;
SPCM_DWGETERRORINFO_I32*    pfn_spcm_dwGetErrorInfo_i32;
SPCM_DWGETCONTBUF_I64M*     pfn_spcm_dwGetContBuf_i64m;
SPCM_DWGETCONTBUF_I64*      pfn_spcm_dwGetContBuf_i64;


int main()
    {
    drv_handle  hDrv;
    int32       lValue;


#ifdef WIN32
    bool b32Bit = true; // set to false if you compile a 64bit program
    HMODULE     hDLL;

    // ----- 32bit -----
    if (b32Bit)
        {
        // we first load the dll and check whether it is available
        hDLL = LoadLibrary ("spcm_win32.dll");
        if (!hDLL)
            {
            printf ("Library spcm_win32.dll can't be located\n");

            return 1;
            }



        // now we load the functions from the DLL using the function names
        pfn_spcm_hOpen =              (SPCM_HOPEN*)               GetProcAddress (hDLL, "_spcm_hOpen@4");
        pfn_spcm_vClose =             (SPCM_VCLOSE*)              GetProcAddress (hDLL, "_spcm_vClose@4");
        pfn_spcm_dwSetParam_i32 =     (SPCM_DWSETPARAM_I32*)      GetProcAddress (hDLL, "_spcm_dwSetParam_i32@12");
        pfn_spcm_dwSetParam_i64m =    (SPCM_DWSETPARAM_I64M*)     GetProcAddress (hDLL, "_spcm_dwSetParam_i64m@16");
        pfn_spcm_dwSetParam_i64 =     (SPCM_DWSETPARAM_I64*)      GetProcAddress (hDLL, "_spcm_dwSetParam_i64@16");
        pfn_spcm_dwGetParam_i32 =     (SPCM_DWGETPARAM_I32*)      GetProcAddress (hDLL, "_spcm_dwGetParam_i32@12");
        pfn_spcm_dwGetParam_i64m =    (SPCM_DWGETPARAM_I64M*)     GetProcAddress (hDLL, "_spcm_dwGetParam_i64m@16");
        pfn_spcm_dwGetParam_i64 =     (SPCM_DWGETPARAM_I64*)      GetProcAddress (hDLL, "_spcm_dwGetParam_i64@12");
        pfn_spcm_dwDefTransfer_i64m = (SPCM_DWDEFTRANSFER_I64M*)  GetProcAddress (hDLL, "_spcm_dwDefTransfer_i64m@36");
        pfn_spcm_dwDefTransfer_i64 =  (SPCM_DWDEFTRANSFER_I64*)   GetProcAddress (hDLL, "_spcm_dwDefTransfer_i64@36");
        pfn_spcm_dwInvalidateBuf =    (SPCM_DWINVALIDATEBUF*)     GetProcAddress (hDLL, "_spcm_dwInvalidateBuf@8");
        pfn_spcm_dwGetErrorInfo_i32 = (SPCM_DWGETERRORINFO_I32*)  GetProcAddress (hDLL, "_spcm_dwGetErrorInfo_i32@16");
        pfn_spcm_dwGetContBuf_i64m =  (SPCM_DWGETCONTBUF_I64M*)   GetProcAddress (hDLL, "_spcm_dwGetContBuf_i64m@20");
        pfn_spcm_dwGetContBuf_i64 =   (SPCM_DWGETCONTBUF_I64*)    GetProcAddress (hDLL, "_spcm_dwGetContBuf_i64@16");
        }

    // ----- 64bit -----
    else
        {
        // we first load the dll and check whether it is available
        hDLL = LoadLibrary ("spcm_win64.dll");
        if (!hDLL)
            {
            printf ("Library spcm_win64.dll can't be located\n");

            return 1;
            }



        // now we load the functions from the DLL using the function names
        pfn_spcm_hOpen =              (SPCM_HOPEN*)               GetProcAddress (hDLL, "spcm_hOpen");
        pfn_spcm_vClose =             (SPCM_VCLOSE*)              GetProcAddress (hDLL, "spcm_vClose");
        pfn_spcm_dwSetParam_i32 =     (SPCM_DWSETPARAM_I32*)      GetProcAddress (hDLL, "spcm_dwSetParam_i32");
        pfn_spcm_dwSetParam_i64m =    (SPCM_DWSETPARAM_I64M*)     GetProcAddress (hDLL, "spcm_dwSetParam_i64m");
        pfn_spcm_dwSetParam_i64 =     (SPCM_DWSETPARAM_I64*)      GetProcAddress (hDLL, "spcm_dwSetParam_i64");
        pfn_spcm_dwGetParam_i32 =     (SPCM_DWGETPARAM_I32*)      GetProcAddress (hDLL, "spcm_dwGetParam_i32");
        pfn_spcm_dwGetParam_i64m =    (SPCM_DWGETPARAM_I64M*)     GetProcAddress (hDLL, "spcm_dwGetParam_i64m");
        pfn_spcm_dwGetParam_i64 =     (SPCM_DWGETPARAM_I64*)      GetProcAddress (hDLL, "spcm_dwGetParam_i64");
        pfn_spcm_dwDefTransfer_i64m = (SPCM_DWDEFTRANSFER_I64M*)  GetProcAddress (hDLL, "spcm_dwDefTransfer_i64m");
        pfn_spcm_dwDefTransfer_i64 =  (SPCM_DWDEFTRANSFER_I64*)   GetProcAddress (hDLL, "spcm_dwDefTransfer_i64");
        pfn_spcm_dwInvalidateBuf =    (SPCM_DWINVALIDATEBUF*)     GetProcAddress (hDLL, "spcm_dwInvalidateBuf");
        pfn_spcm_dwGetErrorInfo_i32 = (SPCM_DWGETERRORINFO_I32*)  GetProcAddress (hDLL, "spcm_dwGetErrorInfo_i32");
        pfn_spcm_dwGetContBuf_i64m =  (SPCM_DWGETCONTBUF_I64M*)   GetProcAddress (hDLL, "spcm_dwGetContBuf_i64m");
        pfn_spcm_dwGetContBuf_i64 =   (SPCM_DWGETCONTBUF_I64*)    GetProcAddress (hDLL, "spcm_dwGetContBuf_i64");
        }

#else // Linux
    void* hLib = dlopen ("libspcm_linux.so", RTLD_NOW);
    if (!hLib)
        {
        printf ("Library libspcm_linux.so can't be located\n");

        return 1;
        }

    pfn_spcm_hOpen =              (SPCM_HOPEN*)              dlsym (hLib, "spcm_hOpen");
    pfn_spcm_vClose =             (SPCM_VCLOSE*)             dlsym (hLib, "spcm_vClose");
    pfn_spcm_dwSetParam_i32 =     (SPCM_DWSETPARAM_I32*)     dlsym (hLib, "spcm_dwSetParam_i32");
    pfn_spcm_dwSetParam_i64m =    (SPCM_DWSETPARAM_I64M*)    dlsym (hLib, "spcm_dwSetParam_i64m");
    pfn_spcm_dwSetParam_i64 =     (SPCM_DWSETPARAM_I64*)     dlsym (hLib, "spcm_dwSetParam_i64");
    pfn_spcm_dwGetParam_i32 =     (SPCM_DWGETPARAM_I32*)     dlsym (hLib, "spcm_dwGetParam_i32");
    pfn_spcm_dwGetParam_i64m =    (SPCM_DWGETPARAM_I64M*)    dlsym (hLib, "spcm_dwGetParam_i64m");
    pfn_spcm_dwGetParam_i64 =     (SPCM_DWGETPARAM_I64*)     dlsym (hLib, "spcm_dwGetParam_i64");
    pfn_spcm_dwDefTransfer_i64m = (SPCM_DWDEFTRANSFER_I64M*) dlsym (hLib, "spcm_dwDefTransfer_i64m");
    pfn_spcm_dwDefTransfer_i64 =  (SPCM_DWDEFTRANSFER_I64*)  dlsym (hLib, "spcm_dwDefTransfer_i64");
    pfn_spcm_dwInvalidateBuf =    (SPCM_DWINVALIDATEBUF*)    dlsym (hLib, "spcm_dwInvalidateBuf");
    pfn_spcm_dwGetErrorInfo_i32 = (SPCM_DWGETERRORINFO_I32*) dlsym (hLib, "spcm_dwGetErrorInfo_i32");
    pfn_spcm_dwGetContBuf_i64m =  (SPCM_DWGETCONTBUF_I64M*)  dlsym (hLib, "spcm_dwGetContBuf_i64m"); 
    pfn_spcm_dwGetContBuf_i64 =   (SPCM_DWGETCONTBUF_I64*)   dlsym (hLib, "spcm_dwGetContBuf_i64"); 
#endif


    // check whether all these loads have been sucessful
    if (!pfn_spcm_hOpen || !pfn_spcm_vClose || !pfn_spcm_dwSetParam_i32 || !pfn_spcm_dwSetParam_i64m ||
        !pfn_spcm_dwSetParam_i64 || !pfn_spcm_dwGetParam_i32 || !pfn_spcm_dwGetParam_i64m || !pfn_spcm_dwGetParam_i64 ||
        !pfn_spcm_dwDefTransfer_i64m || !pfn_spcm_dwDefTransfer_i64 || !pfn_spcm_dwInvalidateBuf || !pfn_spcm_dwGetErrorInfo_i32 ||
        !pfn_spcm_dwGetContBuf_i64m || !pfn_spcm_dwGetContBuf_i64)
        {
        printf ("One of the functions wasn't found inside the library!\n");
#ifdef WIN32
        FreeLibrary (hDLL);
#else
        dlclose (hLib);
#endif

        return 1;
        }



    // as an example we open the first card, read out type and serial number and close the card again
    // uncomment the second line and replace the IP address to use remote
    // cards like in a digitizerNETBOX
    hDrv = (*pfn_spcm_hOpen) ("/dev/spcm0");
    // hDrv = (*pfn_spcm_hOpen) ("TCPIP::192.168.1.10::inst0::INSTR");
    if (hDrv)
        {
        (*pfn_spcm_dwGetParam_i32) (hDrv, SPC_PCITYP, &lValue);
        printf ("Card Type:      M2i.%04x\n", lValue & TYP_VERSIONMASK);

        (*pfn_spcm_dwGetParam_i32) (hDrv, SPC_PCISERIALNO, &lValue);
        printf ("Serial number:  %05d\n", lValue);

        (*pfn_spcm_vClose) (hDrv);
        }
    else
        printf ("no card found ...\n");


#ifdef WIN32
    FreeLibrary (hDLL);
#else
    dlclose (hLib);
#endif

    return 0;
    }

