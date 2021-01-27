/*
 *  findsection.h -- header file for findsection
 *
 *  $Id$
 */

#define MAX_FIT_COL 8
#define MAX_FIT_DATA 50

class FindSect 
{
public:
  FindSect();
  void main( int argc, char *argv[]);
  void spreadFunc(float *yvec, float *spread);
 private:

  // Methods
  void addBoxesToSample(int startX, int endX, int startY, int endY, int startZ,
                        int endZ, int binInd, int &numStat, double &meanSum);
  int findColumnMidpoint(int startX, int endX, int startY, int endY, int binInd,
                         float *buffer, int *izInside, int &izMid, float &insideMedian);
  void fitColumnBoundaries(int startX, int endX, int startY, int endY, float *boundary);
  void dumpPointModel(float *boundaries, int numPts,
                      const char *filename, MrcHeader *inHeader);
  void makeSurfaceModel(float *boundaries, const char *filename, MrcHeader *inHeader,
                        int *nxyz);
  int addToPitchModel(Imod *imod, float *boundaries, int yStart, int yEnd, int time);
  void checkBlockThicknesses();
  void setupBlocks(int numBoxes, int sclInd, int ixyz);
  void invertYifFlipped(int &start, int &end, int *dims);
  void getFittingRegion(int blockCen[], int extentNum, int maxDepth, int botTop,
                        int *blockStart, int *blockEnd, int &numBound, float &xySpacing);
  bool hasBoundaries(int axis, int rowCol, int start, int end, int botTop);
  void loadFittingMatrix(int *ixyCen, int *starts, int *ends, int botTop, float xySpacing,
                         float wgtThresh, int wgtColIn, int wgtColOut, int &numFit);
  void buildVariableList(int major, int minor, int numBound, int numGood);
  int findLowestSDEdges(int *starts, int *ends, int scl);
  int analyzeHighSD(float highSDcrit, Imod *beadModel, int *boxSize, int *startCoord,
                    int *endCoord, Imod *pitchModel);
  void makeCombinedBins(float *values, int numVals, float *bins, float *comboBins,
                        int mult, int &numBins, int &binWidth, float &maxBin, 
                        int &indMax);
  float convexAreaCovered(int numPts, int offset);
  void writeModel(const char *filename, Imod *imod);
  void computeBeadLimits(float excludePctl);

  // Member Variables
  int mNumBoxes[MAX_MBS_SCALES][3];
  int mBinning[MAX_MBS_SCALES][3];
  int mBoxSpacing[MAX_MBS_SCALES][3];
  float *mSDs, *mBuffer, *mMeans, *mColMedians, *mColSlice;
  int mStatStartInds[MAX_MBS_SCALES];
  float mEdgeMedians[MAX_MBS_SCALES], mEdgeMADNs[MAX_MBS_SCALES];
  float mCenMedians[MAX_MBS_SCALES], mCenMADNs[MAX_MBS_SCALES];
  int mNumInBlock[MAX_MBS_SCALES][3], mNumBlocks[MAX_MBS_SCALES][3];
  int mBestLowEdge[MAX_MBS_SCALES], mBestHighEdge[MAX_MBS_SCALES];
  float mFitMat[MAX_FIT_COL][MAX_FIT_DATA];
  float mFitWork[2 * MAX_FIT_DATA + MAX_FIT_COL * MAX_FIT_COL];
  int mGoodRowCol[2], mVarList[5][2], mNxyz[3];
  int mNumVars;
  int mNumThicknesses;
  float mThickMedian, mThickMADN;
  float *mBoundaries, *mBlockCenters, *mThicknesses, *mBoundRot;
  float *mProjMeans, *mProjSlice;
  int *mBestBoxSize, *mStartCoord;
  int mNumHighSD, mNumBeads;
  int mMinNumThickForCheck;
  float mCritThickMADN;
  float mBeadDiameter;
  int mSimplexIter;
  bool mAfterSimplex;
  float mBoostHighSDThickness;
  float mBoundLow, mBoundHigh;
  float mBeadLow, mBeadHigh;
  float mMinSpread, mSpreadNorms[3];
  float mBeadWeightFac, mBeadArea, mStructArea;
  float *mConvexXtmp, *mConvexYtmp;
  Icont *mConvexCont;
  float mColMaxEdgeDiffCrit;
  float mFracColMaxEdgeDiff;
  float mCritEdgeMADN;
  int mNumHighInsideCrit;
  int mThickInd;
  int mFitPitchSeparately;
  int mDebugOutput;
  int mScanBlockSize;
  int mBestScale;
  float mKfactor;
  float mMaxChange;
  float mMaxOscill;
  int mMaxIter;
  float mMaxFalloffFrac;
  float mLowFitFrac;
  float mHighFitFrac;
  float mBoundaryFrac;
  float mPitchBoundaryFrac;
  float mMinFracBoundsInCol;
  float mColumnToCenMADNCrit;
  float mTooThinCrit;
  float mFartherFromMeanCrit;
  float mMeanDiffThickFrac;
  int mNumExtraSum;
  float mExtraForPitch, mExtraPitchSum;
  int mMinForRobustPitch;
  float mBeadExclPctlSpread;
  float mBeadExclPctlThick;
  int mProjLayerStatType;
  float mProjEdgeCrit;
  int mProjUseExtrap;
  int mUseProjForSpread;
};
