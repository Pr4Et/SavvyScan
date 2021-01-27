/*
**************************************************************************

spc_sse_demux.cpp                                        (c) Spectrum GmbH

**************************************************************************

Implements data demultiplexing functions using SSE2 commands by intrinsic
functions. Tested with Microsoft Visual C++ (Windows) and Gnu C++ (Linux)

These functions need a CPU supporting SSE2 commands to run. Please see
the Intrinsic Guide by Intel for more information:
https://software.intel.com/sites/landingpage/IntrinsicsGuide

The demultiplexing functions speed up the process to demultiplex a single
packed data array with 2, 4 or 8 multiplexed channels into separate data
arrays, one for each channel. Compared to a simple nested for-loop these
function can speed up the demultiplexing process up to the factor of 8x

Feel free to use these functions in your own programs.

**************************************************************************
*/

#include "spc_sse_demux.h"

// ----------------------------------------------------------------------------
// ----- Demultiplexing of 2 channels.
// ----- pnData and all anChX must be 16 byte aligned.
// ----- Replace _mm_load_si128/_mm_store_si128 with _mm_loadu_si128/_mm_storeu_si128
// ----- if you want to use unaligned memory.
// ----- This function requires dwNumSamples to be a multiple of 8 because all
// ----- operations take place on 128bit = 8*16bit.
// ----- anChX needs to point to enough memory to hold dwNumSamples.
// ----------------------------------------------------------------------------
void vDemuxData2_int16 (int16* anMuxedData, uint32 dwNumSamples, int16* anCh0, int16* anCh1)
    {
    const uint32 dwNumCh = 2;

    // ----- init pointers to 128 block for easier access later on -----
    __m128i* powReadPos     = reinterpret_cast < __m128i * > (anMuxedData);
    __m128i* powWritePosCh0 = reinterpret_cast < __m128i * > (anCh0);
    __m128i* powWritePosCh1 = reinterpret_cast < __m128i * > (anCh1);

    // ----- demux data. 8 samples per channel in each run of the for loop -----
    __m128i xmmData1, xmmData2;
    __m128i xmmTmp1, xmmTmp2, xmmTmp3;
    const uint32 dwNumSamplesPer128Bit = sizeof(__m128i)/sizeof(int16);
    for (uint32 i = 0; i < dwNumSamples * dwNumCh; i += dwNumSamplesPer128Bit * dwNumCh)
        {
        // load data from memory
        xmmData1 = _mm_load_si128 (powReadPos++); // B3A3B2A2 B1A1B0A0
        xmmData2 = _mm_load_si128 (powReadPos++); // B7A7B6A6 B5A5B4A4

        // sort lower 4 samples
        xmmTmp1 = _mm_shufflelo_epi16 (xmmData1, (3 << 6) | (1 << 4) | (2 << 2) | (0 << 0)); // B3A3B2A2 B1B0A1A0
        xmmTmp1 = _mm_shufflehi_epi16 (xmmTmp1,  (3 << 6) | (1 << 4) | (2 << 2) | (0 << 0)); // B3B2A3A2 B1B0A1A0
        xmmTmp1 = _mm_shuffle_epi32   (xmmTmp1,  (3 << 6) | (1 << 4) | (2 << 2) | (0 << 0)); // B3B2B1B0 A3A2A1A0

        // sort upper 4 samples
        xmmTmp2 = _mm_shufflelo_epi16 (xmmData2, (3 << 6) | (1 << 4) | (2 << 2) | (0 << 0)); // B7A7B6A6 B5B4A5A4
        xmmTmp2 = _mm_shufflehi_epi16 (xmmTmp2,  (3 << 6) | (1 << 4) | (2 << 2) | (0 << 0)); // B7B6A7A6 B5B4A5A4
        xmmTmp2 = _mm_shuffle_epi32   (xmmTmp2,  (3 << 6) | (1 << 4) | (2 << 2) | (0 << 0)); // B7B6B5B4 A7A6A5A4

        xmmTmp3 = _mm_unpacklo_epi64 (xmmTmp1, xmmTmp2); // A7A6A5A4 A3A2A1A0
        _mm_store_si128 (powWritePosCh0, xmmTmp3);

        xmmTmp3 = _mm_unpackhi_epi64 (xmmTmp1, xmmTmp2); // B7B6B5B4 B3B2B1B0
        _mm_store_si128 (powWritePosCh1, xmmTmp3);

        // ----- increase write position in each demuxed buffer -----
        powWritePosCh0++;
        powWritePosCh1++;
        }
    }

// ----------------------------------------------------------------------------
// ----- Demultiplexing of 4 channels.
// ----- pnData and all anChX must be 16 byte aligned.
// ----- Replace _mm_load_si128/_mm_store_si128 with _mm_loadu_si128/_mm_storeu_si128
// ----- if you want to use unaligned memory.
// ----- This function requires dwNumSamples to be a multiple of 8 because all
// ----- operations take place on 128bit = 8*16bit.
// ----- anChX needs to point to enough memory to hold dwNumSamples.
// ----------------------------------------------------------------------------
void vDemuxData4_int16 (int16* anMuxedData, uint32 dwNumSamples, int16* anCh0, int16* anCh1, int16* anCh2, int16* anCh3)
    {
    const uint32 dwNumCh = 4;

    // ----- init pointers to 128 block for easier access later on -----
    __m128i* powReadPos     = reinterpret_cast < __m128i * > (anMuxedData);
    __m128i* powWritePosCh0 = reinterpret_cast < __m128i * > (anCh0);
    __m128i* powWritePosCh1 = reinterpret_cast < __m128i * > (anCh1);
    __m128i* powWritePosCh2 = reinterpret_cast < __m128i * > (anCh2);
    __m128i* powWritePosCh3 = reinterpret_cast < __m128i * > (anCh3);

    // ----- demux data. 8 samples per channel in each run of the for loop -----
    __m128i xmmData1, xmmData2, xmmData3, xmmData4;
    __m128i xmmTmp1, xmmTmp2;
    const uint32 dwNumSamplesPer128Bit = sizeof(__m128i)/sizeof(int16);
    for (uint32 i = 0; i < dwNumSamples * dwNumCh; i += dwNumSamplesPer128Bit * dwNumCh)
        {
        xmmData1 = _mm_load_si128 (powReadPos++); // D1C1B1A1 D0C0B0A0
        xmmData2 = _mm_load_si128 (powReadPos++); // D3C3B3A3 D2C2B2A2
        xmmData3 = _mm_load_si128 (powReadPos++); // D5C5B5A5 D4C4B4A4
        xmmData4 = _mm_load_si128 (powReadPos++); // D7C7B7A7 D6C6B6A6

        // lower four samples
        xmmTmp1 = _mm_unpacklo_epi16 (xmmData1, xmmData2); // D2D0 C2C0 B2B0 A2A0
        xmmTmp2 = _mm_unpackhi_epi16 (xmmData1, xmmData2); // D3D1 C3C1 B3B1 A3A1

        xmmData1 = _mm_unpacklo_epi16 (xmmTmp1, xmmTmp2); // B3B2B1B0 A3A2A1A0
        xmmData2 = _mm_unpackhi_epi16 (xmmTmp1, xmmTmp2); // D3D2D1D0 C3C2C1C0

        // upper four samples
        xmmTmp1 = _mm_unpacklo_epi16 (xmmData3, xmmData4); // D6D4 C6C4 B6B4 A6A4
        xmmTmp2 = _mm_unpackhi_epi16 (xmmData3, xmmData4); // D7D5 C7C5 B7B5 A7A5

        xmmData3 = _mm_unpacklo_epi16 (xmmTmp1, xmmTmp2); // B7B6B5B4 A7A6A5A4
        xmmData4 = _mm_unpackhi_epi16 (xmmTmp1, xmmTmp2); // D7D6D5D4 C7C6C5C4

        xmmTmp1 = _mm_unpacklo_epi64 (xmmData1, xmmData3); // A7A6A5A4 A3A2A1A0
        _mm_store_si128 (powWritePosCh0, xmmTmp1);

        xmmTmp1 = _mm_unpackhi_epi64 (xmmData1, xmmData3); // B7B6B5B4 B3B2B1B0
        _mm_store_si128 (powWritePosCh1, xmmTmp1);

        xmmTmp1 = _mm_unpacklo_epi64 (xmmData2, xmmData4); // C7C6C5C4 C3C2C1C0
        _mm_store_si128 (powWritePosCh2, xmmTmp1);

        xmmTmp1 = _mm_unpackhi_epi64 (xmmData2, xmmData4); // D7D6D5D4 D3D2D1D0
        _mm_store_si128 (powWritePosCh3, xmmTmp1); 

        // ----- increase write position in each demuxed buffer -----
        powWritePosCh0++;
        powWritePosCh1++;
        powWritePosCh2++;
        powWritePosCh3++;
        }
    }

// ----------------------------------------------------------------------------
// ----- Demultiplexing of 8 channels.
// ----- pnData and all anChX must be 16 byte aligned.
// ----- Replace _mm_load_si128/_mm_store_si128 with _mm_loadu_si128/_mm_storeu_si128
// ----- if you want to use unaligned memory.
// ----- This function requires dwNumSamples to be a multiple of 8 because all
// ----- operations take place on 128bit = 8*16bit.
// ----- anChX needs to point to enough memory to hold dwNumSamples.
// ----------------------------------------------------------------------------
void vDemuxData8_int16 (int16* anMuxedData, uint32 dwNumSamples, int16* anCh0, int16* anCh1, int16* anCh2, int16* anCh3, int16* anCh4, int16* anCh5, int16* anCh6, int16* anCh7)
    {
    const uint32 dwNumCh = 8;

    // ----- init pointers to 128 block for easier access later on -----
    __m128i* powReadPos     = reinterpret_cast < __m128i * > (anMuxedData);
    __m128i* powWritePosCh0 = reinterpret_cast < __m128i * > (anCh0);
    __m128i* powWritePosCh1 = reinterpret_cast < __m128i * > (anCh1);
    __m128i* powWritePosCh2 = reinterpret_cast < __m128i * > (anCh2);
    __m128i* powWritePosCh3 = reinterpret_cast < __m128i * > (anCh3);
    __m128i* powWritePosCh4 = reinterpret_cast < __m128i * > (anCh4);
    __m128i* powWritePosCh5 = reinterpret_cast < __m128i * > (anCh5);
    __m128i* powWritePosCh6 = reinterpret_cast < __m128i * > (anCh6);
    __m128i* powWritePosCh7 = reinterpret_cast < __m128i * > (anCh7);

    // ----- demux data. 8 samples per channel in each run of the for loop -----
    __m128i xmmData1, xmmData2, xmmData3, xmmData4;
    __m128i xmmTmp1, xmmTmp2, xmmTmp3, xmmTmp4, xmmTmp5, xmmTmp6, xmmTmp7;
    __m128i xmmResult;
    const uint32 dwNumSamplesPer128Bit = sizeof(__m128i)/sizeof(int16);
    for (uint32 i = 0; i < dwNumSamples * dwNumCh; i += dwNumSamplesPer128Bit * dwNumCh)
        {
        // ----- load and sort lower four samples -----
        xmmData1 = _mm_load_si128 (powReadPos++); // H0G0F0E0 D0C0B0A0
        xmmData2 = _mm_load_si128 (powReadPos++); // H1G1F1E1 D1C1B1A1
        xmmData3 = _mm_load_si128 (powReadPos++); // H2G2F2E2 D2C2B2A2
        xmmData4 = _mm_load_si128 (powReadPos++); // H3G3F3E3 D3C3B3A3

        xmmTmp1 = _mm_unpacklo_epi16 (xmmData1, xmmData2); // D1D0C1C0 B1B0A1A0
        xmmTmp2 = _mm_unpacklo_epi16 (xmmData3, xmmData4); // D3D2C3C2 B3B2A3A2
        xmmTmp3 = _mm_unpacklo_epi32 (xmmTmp1, xmmTmp2);   // B3B2B1B0 A3A2A1A0
        xmmTmp4 = _mm_unpackhi_epi32 (xmmTmp1, xmmTmp2);   // D3D2D1D0 C3C2C1C0

        xmmTmp1 = _mm_unpackhi_epi16 (xmmData1, xmmData2); // H1H0G1G0 F1F0E1E0
        xmmTmp2 = _mm_unpackhi_epi16 (xmmData3, xmmData4); // H3H2G3G2 F3F2E3E2
        xmmTmp5 = _mm_unpacklo_epi32 (xmmTmp1, xmmTmp2);   // F3F2F1F0 E3E2E1E0
        xmmTmp6 = _mm_unpackhi_epi32 (xmmTmp1, xmmTmp2);   // H3H2H1H0 G3G2G1G0

        // ----- load and sort upper four samples -----
        xmmData1 = _mm_load_si128 (powReadPos++); // H4G4F4E4 D4C4B4A4
        xmmData2 = _mm_load_si128 (powReadPos++); // H5G5F5E5 D5C5B5A5
        xmmData3 = _mm_load_si128 (powReadPos++); // H6G6F6E6 D6C6B6A6
        xmmData4 = _mm_load_si128 (powReadPos++); // H7G7F7E7 D7C7B7A7

        xmmTmp1 = _mm_unpacklo_epi16 (xmmData1, xmmData2); // D5D4C5C4 B5B4A5A4
        xmmTmp2 = _mm_unpacklo_epi16 (xmmData3, xmmData4); // D7D6C7C6 B7B6A7A6
        xmmTmp7 = _mm_unpacklo_epi32 (xmmTmp1, xmmTmp2);   // B7B6B5B4 A7A6A5A4

        // ----- merge lower and upper four samples and write to result buffer -----
        xmmResult = _mm_unpacklo_epi64 (xmmTmp3, xmmTmp7); // A7A6A5A4 A3A2A1A0
        _mm_store_si128 (powWritePosCh0, xmmResult);

        xmmResult = _mm_unpackhi_epi64 (xmmTmp3, xmmTmp7); // B7B6B5B4 B3B2B1B0
        _mm_store_si128 (powWritePosCh1, xmmResult);

        xmmTmp7 = _mm_unpackhi_epi32 (xmmTmp1, xmmTmp2);   // D7D6D5D4 C7C6C5C4

        xmmResult = _mm_unpacklo_epi64 (xmmTmp4, xmmTmp7); // C7C6C5C4 C3C2C1C0
        _mm_store_si128 (powWritePosCh2, xmmResult);

        xmmResult = _mm_unpackhi_epi64 (xmmTmp4, xmmTmp7); // D7D6D5D4 D3D2D1D0
        _mm_store_si128 (powWritePosCh3, xmmResult);


        xmmTmp1  = _mm_unpackhi_epi16 (xmmData1, xmmData2); // H5H4G5G4 F5F4E5E4
        xmmTmp2  = _mm_unpackhi_epi16 (xmmData3, xmmData4); // H7H6G7G6 F7F6E7E6
        xmmTmp7  = _mm_unpacklo_epi32 (xmmTmp1, xmmTmp2);   // F7F6F5F4 E7E6E5E4

        xmmResult = _mm_unpacklo_epi64 (xmmTmp5, xmmTmp7); // E7E6E5E4 E3E2E1E0
        _mm_store_si128 (powWritePosCh4, xmmResult);

        xmmResult = _mm_unpackhi_epi64 (xmmTmp5, xmmTmp7); // F7F6F5F4 F3F2F1F0
        _mm_store_si128 (powWritePosCh5, xmmResult);

        xmmTmp7 = _mm_unpackhi_epi32 (xmmTmp1, xmmTmp2);   // H7H6H5H4 G7G6G5G4

        xmmResult = _mm_unpacklo_epi64 (xmmTmp6, xmmTmp7); // G7G6G5G4 G3G2G1G0
        _mm_store_si128 (powWritePosCh6, xmmResult);

        xmmResult = _mm_unpackhi_epi64 (xmmTmp6, xmmTmp7); // H7H6H5H4 H3H2H1H0
        _mm_store_si128 (powWritePosCh7, xmmResult);

        // ----- increase write position in each demuxed buffer -----
        powWritePosCh0++;
        powWritePosCh1++;
        powWritePosCh2++;
        powWritePosCh3++;
        powWritePosCh4++;
        powWritePosCh5++;
        powWritePosCh6++;
        powWritePosCh7++;
        }
    }
