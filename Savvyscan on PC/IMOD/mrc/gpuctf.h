/*
 *  gpuctf.h -- header file for gpuctf.cu
 *
 *  $Id$
 */
#ifndef GPUCTF_H
#define GPUCTF_H
int gpuInitializeSlice(float *sliceData, int nxFile, int nyFile, int stripXdim, int nxPad, 
                       int nyPad, bool doFullImages);
int gpuExtractAndTransform(int stripInd, int stripBegin, int stripEnd, int nxTaper,
                           int nyTaper);
int gpuCorrectCTF(int stripInd, float freq_scalex, float freq_scaley, float pointDefocus,
                  float cosAstig, float sinAstig, float focusSum, float focusDiff,
                  float cutonAngstroms, float phaseFracFactor, float phaseShift, 
                  float ampAngle, float C1, float C2, float scaleByPower, bool powerIsHalf,
                  bool generalPower, float firstZeroFreqSq, float attenStartFrac,
                  float minAttenFreqSq);
int gpuInterpolateColumns(int stripInd, int yoff, int stripStride, int stripMid,
                          int halfStrip, int curOffset, int lastOffset);
int gpuCopyColumns(int stripInd, int xoff, int yoff, int startCol, int endCol);
int gpuInterpDiagonals(int stripInd, int xoff, int yoff, int stripStride, 
                       float sinViewAxis, float cosViewAxis, float lastAxisDist,
                       float curAxisDist);
int gpuCopyDiagonals(int stripInd, int xoff, int yoff, float sinViewAxis, 
                     float cosViewAxis, float lowLim, float highLim);
int gpuReturnImage(float *finalImage);
int gpuAvailable(int nGPU, float *memory, int debug);
void gpuGetTimes(double &copy, double &prep, double &FFT, double &correct, double &interp);
#endif
