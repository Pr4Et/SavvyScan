#include "gpuframe.h"
#include "framealign.h"

int fgpuGpuAvailable(int nGPU, float *memory, int debug) 
{
  *memory = 0.;
  return 0;
}

int fgpuSetupSumming(int fullXpad, int fullYpad, int sumXpad, int sumYpad, int evenOdd) {return 1;}
int fgpuSetupAligning(int alignXpad, int alignYpad, int sumXpad, int sumYpad,
                    float *alignMask, int aliFiltSize, int groupSize, int expectStackSize,
                      int doAlignSum) {return 1;}
int fgpuAddToFullSum(float *fullArr, float shiftX, float shiftY) {return 1;}
int fgpuReturnSums(float *sumArr, float *evenArr, float *oddArr, int evenOddOnly) {return 1;}
void fgpuCleanup() {}
void fgpuRollAlignStack() {}
void fgpuRollGroupStack() {}
int fgpuSubtractAndFilterAlignSum(int stackInd, int groupRefine) {return 1;}
int fgpuNewFilterMask(float *alignMask) {return 1;}
int fgpuShiftAddToAlignSum(int stackInd, float shiftX, float shiftY, int shiftSource) {return 1;}
int fgpuCrossCorrelate(int aliInd, int refInd, float *subarea, int subXoffset,
                       int subYoffset) {return 1;}
int fgpuProcessAlignImage(float *binArr, int stackInd, int groupInd) {return 1;}
int fgpuReturnAlignFFTs(std::vector<float *> *saved, std::vector<float *> *groups,
                        float *alignSum, float *workArr) {return 1;}
void fgpuCleanSumItems() {}
void fgpuCleanAlignItems() {}
void fgpuZeroTimers() {}
void fgpuPrintTimers() {}
int fgpuClearAlignSum() {return 1;}
int fgpuSumIntoGroup(int stackInd, int groupInd) {return 1;}
void fgpuSetGroupSize(int inVal) {}
void fgpuSetPrintFunc(CharArgType func) {}
