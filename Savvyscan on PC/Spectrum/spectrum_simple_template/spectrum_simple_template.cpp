/*
**************************************************************************

spectrum_simple_template.cpp                             (c) Spectrum GmbH

**************************************************************************

Simple template as base for own programming

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

#include "../common/ostools/spcm_oswrap.h"
#include "../common/ostools/spcm_ostools.h"

// ----- standard c include files -----
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


/*
**************************************************************************
szTypeToName: doing name translation
**************************************************************************
*/

char* szTypeToName (int32 lCardType)
    {
    static char szName[50];
    switch (lCardType & TYP_SERIESMASK)
        {
        case TYP_M2ISERIES:     sprintf (szName, "M2i.%04x", lCardType & TYP_VERSIONMASK);      break;
        case TYP_M2IEXPSERIES:  sprintf (szName, "M2i.%04x-Exp", lCardType & TYP_VERSIONMASK);  break;
        case TYP_M3ISERIES:     sprintf (szName, "M3i.%04x", lCardType & TYP_VERSIONMASK);      break;
        case TYP_M3IEXPSERIES:  sprintf (szName, "M3i.%04x-Exp", lCardType & TYP_VERSIONMASK);  break;
        case TYP_M4IEXPSERIES:  sprintf (szName, "M4i.%04x-x8", lCardType & TYP_VERSIONMASK);   break;
        case TYP_M4XEXPSERIES:  sprintf (szName, "M4x.%04x-x4", lCardType & TYP_VERSIONMASK);   break;
        case TYP_M2PEXPSERIES:  sprintf (szName, "M2p.%04x-x4", lCardType & TYP_VERSIONMASK);   break;
        default:                sprintf (szName, "unknown type");                               break;
        }
    return szName;
    }



/*
**************************************************************************
main 
**************************************************************************
*/

int main ()
    {
    drv_handle  hCard;
    int32       lCardType, lSerialNumber;
    char        szName[50];

    // ------------------------------------------------------------------------
    // we try to open one card and printout some information on it

    // uncomment the second line and replace the IP address to use remote
    // cards like in a digitizerNETBOX
    sprintf (szName, "/dev/spcm0");
    // sprintf (szName, "TCPIP::192.168.1.10::inst0::INSTR");
    hCard = spcm_hOpen (szName);

    // not one card found
    if (!hCard)
        {
        printf ("no card found...\n");
        return 0;
        }

    // read out some info and print it
    spcm_dwGetParam_i32 (hCard, SPC_PCITYP,            &lCardType);
    spcm_dwGetParam_i32 (hCard, SPC_PCISERIALNO,       &lSerialNumber);
    printf ("Found: %s sn %05d\n", szTypeToName (lCardType), lSerialNumber);

    // ------------------------------------------------------------------------
    // start here with programming
    // ------------------------------------------------------------------------

    spcm_vClose (hCard);

    return EXIT_SUCCESS;
    }

