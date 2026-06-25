/*
**************************************************************************

spc_sse_demux.h                                (c) Spectrum GmbH , 03/2016

**************************************************************************

These functions demultiplex the data from one muxed DMA buffer to multiple
buffers where each contains only the samples from one channel.

Feel free to use these functions in your own programs.

**************************************************************************
*/

#ifndef SPC_SSE_DEMUX_H
#define SPC_SSE_DEMUX_H

// ----- include standard driver header from library -----
#include "../c_header/dlltyp.h"

// ----- SSE2 Intrinsics -----
#include <emmintrin.h>

// Demultiplexing of 2 channels. See source file for more details.
void vDemuxData2_int16 (int16* anMuxedData, uint32 dwNumSamples, int16* anCh0, int16* anCh1);

// Demultiplexing of 4 channels. See source file for more details.
void vDemuxData4_int16 (int16* anMuxedData, uint32 dwNumSamples, int16* anCh0, int16* anCh1, int16* anCh2, int16* anCh3);

// Demultiplexing of 8 channels. See source file for more details.
void vDemuxData8_int16 (int16* anMuxedData, uint32 dwNumSamples, int16* anCh0, int16* anCh1, int16* anCh2, int16* anCh3, int16* anCh4, int16* anCh5, int16* anCh6, int16* anCh7);

#endif // SPC_SSE_DEMUX_H
