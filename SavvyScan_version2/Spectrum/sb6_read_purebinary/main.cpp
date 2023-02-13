/*
**************************************************************************

sb6_read_purebinary.cpp                                      (c) Spectrum GmbH

**************************************************************************

Example shows how to read a SBench6 pure binary export file
  
Feel free to use this source for own projects and modify it in any kind

Documentation for the API as well as a detailed description of the hardware
can be found in the manual for each device which can be found on our website:
https://www.spectrum-instrumentation.com/en/downloads

Further information can be found online in the Knowledge Base:
https://www.spectrum-instrumentation.com/en/knowledge-base-overview

**************************************************************************
*/

#include "sb6_read_purebinary.h"
#include <iostream>
#include <iomanip>

#define FILE_ANLOG   0
#define FILE_DIGITAL 1 

using namespace std;

int main ()
    {
    HeaderInfo stHeaderInfo;
    vector <float*> vpfChData;
    vector <unsigned char*> vpbyChData;

    // set path to pure binary file
    string sPureBinaryFilePath = "export.bin";

    // create path to header file
    string sHeaderFilePath = sPureBinaryFilePath;
    sHeaderFilePath = sHeaderFilePath.erase (sHeaderFilePath.find_last_of (".")) + "_binheader.txt"; 

    // read infos from header file and store them in header info struct
    if (SB6_bReadPureBinaryHeaderInfos (sHeaderFilePath, stHeaderInfo))
        {
        // this example reads only pure analog or pure digital files
        if (stHeaderInfo.dwNumAChannels > 0 && stHeaderInfo.dwNumDChannels > 0)
            return 0;

        int lFileType = FILE_ANLOG;
        if (stHeaderInfo.dwNumDChannels > 0)
            lFileType = FILE_DIGITAL;

        switch (lFileType)
            {
            case FILE_ANLOG:
                // read data from binary file.
                // data is stored in vector:
                // vector[0] => data for channel 0
                // vector[1] => data for channel 1
                // ...
                // vector[n] => data for channel n
                vpfChData = SB6_ReadPureBinaryFileAnalog (sPureBinaryFilePath, stHeaderInfo);
        
                if (vpfChData.size ())
                    {
                    // plot first 32 samples for each channel
                    for (uint32 dwChIdx = 0; dwChIdx < vpfChData.size (); dwChIdx++)
                        {
                        cout << "Plot first 32 samples of " + stHeaderInfo.vstChInfos[dwChIdx].sName + ":" << endl;
                        for (uint64 qwDataIdx = 0; qwDataIdx < stHeaderInfo.qwLen && qwDataIdx < 32; qwDataIdx++)
                            cout << vpfChData[dwChIdx][qwDataIdx] << " mV\n";

                        cout << endl;
                        }

                    // free memory
                    SB6_vClearChData (vpfChData);
                    }
                break;

            case FILE_DIGITAL:
                // read data from binary file.
                // data is stored in vector:
                // vector[0] => data for digital channel group D07-D00
                // vector[1] => data for digital channel group D15-D08
                // vector[2] => data for digital channel group D23-D16
                // vector[3] => data for digital channel group D31-D24
                // ... 

                vpbyChData = SB6_ReadPureBinaryFileDigital (sPureBinaryFilePath, stHeaderInfo);

                if (vpbyChData.size ())
                    {
                    int lMaxDigChIndex = 8 * vpbyChData.size () - 1;

                    cout << "Plot first 128 samples for each channel:\n";

                    for (uint64 qwDataIdx = 0; qwDataIdx < stHeaderInfo.qwLen && qwDataIdx < 128; qwDataIdx++)
                        {
                        cout << "[D" << dec << lMaxDigChIndex << "-D0]: 0x";
                        for (int32 lChGroupIdx = vpbyChData.size () - 1; lChGroupIdx >= 0; lChGroupIdx--)
                            cout << hex << setfill ('0') << setw (2) << (uint16)vpbyChData[lChGroupIdx][qwDataIdx];
                           
                        cout << endl;
                        }

                    // free memory
                    SB6_vClearChData (vpbyChData);
                    }
                break;
            }
        }
    
    return 0;
    }
