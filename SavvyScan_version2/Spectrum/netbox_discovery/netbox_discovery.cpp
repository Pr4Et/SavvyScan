// ----- include standard driver header from library -----
#include "../c_header/dlltyp.h"
#include "../c_header/regs.h"
#include "../c_header/spcerr.h"
#include "../c_header/spcm_drv.h"

#include <string>
#include <vector>
#include <iostream>
#include <cstring> // memset

#define TIMEOUT_DISCOVERY 5000 // timeout value in ms

std::vector <std::string> vsGetAvailableVisaStrings ()
    {
    const uint32 dwMaxNumRemoteCards = 50;
    const uint32 dwMaxVisaStringLen = 50;
    const uint32 dwMaxIdnStringLen = 256;

    char* pszVisa[dwMaxNumRemoteCards] = { NULL };
    char* pszIdn[dwMaxNumRemoteCards] = { NULL };

    std::vector <std::string> vsVisa;

    // allocate memory for string list
    for (uint32 i = 0; i < dwMaxNumRemoteCards; i++)
        {
        pszVisa[i] = new char [dwMaxVisaStringLen];
        pszIdn[i] = new char [dwMaxIdnStringLen];
        memset (pszVisa[i], 0, dwMaxVisaStringLen);
        memset (pszIdn[i], 0, dwMaxIdnStringLen);
        }

    // first make discovery - check if there are any LXI compatible remote devices
    uint32 dwError = spcm_dwDiscovery ((char**)pszVisa, dwMaxNumRemoteCards, dwMaxVisaStringLen, TIMEOUT_DISCOVERY);

    // second: check from which manufacturer the devices are
    spcm_dwSendIDNRequest ((char**)pszIdn, dwMaxNumRemoteCards, dwMaxIdnStringLen);

     // ----- store VISA strings for all discovered cards and open them afterwards -----
    for (uint32 i = 0; i < dwMaxNumRemoteCards; i++)
        {
        std::string sIDN (pszIdn[i]);
        
        if (sIDN.size () == 0)
            break;

        if (sIDN.find("Spectrum") != std::string::npos)
            vsVisa.push_back (pszVisa[i]);
        }

    for (uint32 i = 0; i < dwMaxNumRemoteCards; i++)
        {
        delete [] pszVisa[i];
        delete [] pszIdn[i];
        }
    
    return vsVisa;
    }

/*
**************************************************************************
main 
**************************************************************************
*/

int main ()
    {
    std::vector <std::string> vsVisaStrings;
    std::vector <int32> vlNetboxSerialnumbers;

    std::cout << "Netbox discovery running...\n\n";
    vsVisaStrings = vsGetAvailableVisaStrings ();

    for (uint32 i = 0; i < vsVisaStrings.size (); i++)
        {
        // open card
        drv_handle hCard = spcm_hOpen (vsVisaStrings[i].c_str ());
        if (!hCard)
            {
            std::cout << "no card found...\n";
            return 0;
            }

        int32 lSN;
        spcm_dwGetParam_i32 (hCard, SPC_NETBOX_SERIALNO, &lSN);
        vlNetboxSerialnumbers.push_back (lSN);

        // close card
        spcm_vClose (hCard);
        }

    if (vsVisaStrings.size () == vlNetboxSerialnumbers.size ())
        {
        if (vsVisaStrings.size ())
            std::cout << "Netboxes found:\n";
        else
            std::cout << "No Netboxes found!\n";

        for (uint32 i = 0; i < vsVisaStrings.size (); i++)
           std::cout << "SN: " << vlNetboxSerialnumbers[i] << ", Visa: " << vsVisaStrings[i] << std::endl;
        }

    return EXIT_SUCCESS;
    }

   
