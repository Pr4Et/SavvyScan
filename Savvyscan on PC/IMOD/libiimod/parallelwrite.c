/*
 * parallelwrite.c  - Functions for direct writing in parallel to one file
 *
 * Copyright (C) 2009-2016 by the Regents of the University of 
 * Colorado.  See dist/COPYRIGHT for full notice.
 *
 * $Id$
 */

/* 
 * NOTE: This file is a merger of C routines formerly in libcfshr, the C parallel writing
 * in mrcfiles.c, and C translations of Fortran routines formerly in libhvem.  It has not
 * been fully integrated and there are some redundancies.
 */

#include "imodconfig.h"
#include "b3dutil.h"
#include <string.h>
#include "iimage.h"
#include "iiunit.h"

#ifdef F77FUNCAP
#define parwrtinitialize PARWRTINITIALIZE
#define parwrtfindregion PARWRTFINDREGION
#define parwrtgetregion PARWRTGETREGION
#define parwrtproperties PARWRTPROPERTIES
#define parwrtsetcurrent PARWRTSETCURRENT
#define parwrtposn PARWRTPOSN
#define parwrtsec PARWRTSEC
#define parwrtlin PARWRTLIN
#define parwrtclose PARWRTCLOSE
#define iiuwritedummysectohdf IIUWRITEDUMMYSECTOHDF
#define iiuparwrtreclosehdf IIUPARWRTRECLOSEHDF
#define iiuparwrtflushbuffers IIUPARWRTFLUSHBUFFERS
#else
#define parwrtinitialize parwrtinitialize_
#define parwrtfindregion parwrtfindregion_
#define parwrtgetregion parwrtgetregion_
#define parwrtproperties parwrtproperties_
#define parwrtsetcurrent parwrtsetcurrent_
#define parwrtposn parwrtposn_
#define parwrtsec parwrtsec_
#define parwrtlin parwrtlin_
#define parwrtclose parwrtclose_
#define iiuwritedummysectohdf iiuwritedummysectohdf_
#define iiuparwrtreclosehdf iiuparwrtreclosehdf_
#define iiuparwrtflushbuffers iiuparwrtflushbuffers_
#endif

#define MAX_INFOS  5
#define MAXLINE 1024
#define DEFAULT_HDF_BUF_MB  20
#define HDF_BUFFER_ENV_VAR "PARALLEL_HDF_BUF_SIZE"

/* Structures for boundary region definitions */
typedef struct {
  char *file;
  int section[2];
  int startLine[2];
} BoundRegion;

typedef struct {
  BoundRegion *regions;
  int nx, ny, numBoundLines, numFiles, everySec, hdfIndex;
} BoundInfo;

/* Structure for keeping tracking of segments of data in buffers for HDF files */
typedef struct {
  char *bufStart;
  int numBytes;
  int section;
  int startLine;   /* Starting line or -1 or -2 for full section (iiu vs mrc_) */
} BufSegment;

/* Variables for all file access types */
static BoundInfo sInfos[MAX_INFOS];
static int sNumInfos = 0;
static int sCurInfo = -1;

/* HDF Buffer variables */
static char *sHdfBufs[MAX_INFOS];
static BufSegment *sSegments[MAX_INFOS];
static int sNumSegments[MAX_INFOS];
static int sMaxSegments[MAX_INFOS];
static int sBufSize = DEFAULT_HDF_BUF_MB * 1048576;
static int sBufIndex[MAX_INFOS];

/* Variables for unit-based functions */
static int sIzCur[MAX_INFOS], sIyCur[MAX_INFOS], sNxFull[MAX_INFOS], sNyFull[MAX_INFOS];
static int sNzFull[MAX_INFOS], sLinesBound[MAX_INFOS], sIunitBound[MAX_INFOS];
static int sIfOpen[MAX_INFOS], sIzBound[MAX_INFOS][2], sIyBound[MAX_INFOS][2];
static int sIfAllSec[MAX_INFOS];

/* Local functions */
static int parWrtPrepareHDF(ImodImageFile *iiFile);
static int iiuPrepareWriteHDF(int iunit);
static void pwOpenIfNeeded(int izSec, int iyLine, int nlinesWrite, int *ierr);
static int parWrtFindRegion(int secNum, int lineNum, int nlWrite, char **filename, 
                            int *sections, int *startLines);
static void advanceSection();
static void advanceLine();
static int addSegment(char *buf, int numBytes, int section, int startLine,
                      ImodImageFile *iiFile, MrcHeader *hdata, int iunit);
static int writeBufferedLines(int iunit, int section, int startLine, char *linesBuf,
                              int *numLines);
static int writeSegments(ImodImageFile *iiFile, MrcHeader *hdata, int iunit);
static void clearSegments(int infoInd);

/*!
 * Initialize parallel writing to an image file by reading in the boundary 
 * info file whose name is [filename], and storing various values in static 
 * variables. [filename] can be NULL or an empty string, in which case the count of
 * parallel writing files is incremented and this info file index is marked as having no 
 * info file.  [nxin] and [nyin] specify the dimensions of the images.  Pass the negative
 * of [nxin] to indicate that the file is an HDF file, in which case the file specified
 * by [filename] will not be read and will simply be used as a lock file.  The format of 
 * the info file is to start with a line with these entries: ^
 * Version  type  nx  numLines  numFiles  ^
 * version = 1 for now  ^
 * type = 0 for chunks in Z, 1 for chunks in Y that require boundaries on each
 * section  ^
 * nx = X dimension of image file  ^
 * numLines = number of lines written at each boundary  ^
 * numFiles = number of boundary files  ^
 * For each boundary file, there are then two lines:  ^
 * Name of boundary file  ^
 * Boundary 1 section and starting line, boundary 2 section and starting line ^
 * Sections and lines are numbered from 0.  For chunks in Z, the section number
 * should be -1 for no boundary (i.e., for boundary 1 of first chunk and 
 * boundary 2 of last chunk).  A line number of -1 indicates that the boundary
 * extends to the end of the section (thus, a setup script does not need to 
 * know the size of the images in Y).  This routine will convert a line number
 * of -1 to the appropriate starting line.  ^
 * For chunks in Y, the section number is
 * ignored (and should be -1) and the line number should be -1 for no boundary.
 * Returns 1 for failure to open file, 2 for error reading file, 3 for 
 * inappropriate values in header line of file, 4 for memory allocation errors, 
 * 5 for trying to open too many boundary info files for the array, or 6 for an error
 * opening the given file as a lock file for HDF.
 */
int parWrtInitialize(const char *filename, int nxin, int nyin)
{
  FILE *fp;
  int i, version;
  char line[MAXLINE];
  BoundInfo *bi;

  if (sNumInfos >= MAX_INFOS)
    return 5;
  bi = &sInfos[sNumInfos];
  bi->hdfIndex = -1;
  if (!filename || !strlen(filename)) {
    bi->regions = NULL;
    sCurInfo = sNumInfos++;
    return 0;
  }

  /* If nxin is negative, open the file as a lock file for HDF operations, save the 
     index and the specified size */
  if (nxin < 0) {
    bi->hdfIndex = b3dOpenLockFile(filename);
    if (bi->hdfIndex < 0)
      return 6;
    bi->nx = -nxin;
    bi->ny = nyin;
    if (getenv(HDF_BUFFER_ENV_VAR)) {
      i = atoi(getenv(HDF_BUFFER_ENV_VAR));
      if (i > 0)
        sBufSize = i * 1048576;
    }
    sHdfBufs[sNumInfos] = B3DMALLOC(char, sBufSize);
    if (!sHdfBufs[sNumInfos])
      return 7;
    clearSegments(sNumInfos);
    sCurInfo = sNumInfos++;
    return 0;
  }

  fp = fopen(filename, "r");
  if (!fp)
    return 1;
  if (fgetline(fp, line, MAXLINE) <= 0)
    return 2;

  /* Format of file is: version #, type, nx, # of bound lines, # of files */
  sscanf(line, "%d %d %d %d %d", &version, &bi->everySec, &bi->nx, &bi->numBoundLines,
         &bi->numFiles);
  if (bi->nx != nxin || bi->numBoundLines <= 0 || bi->numFiles <= 0)
    return 3;
  bi->regions = (BoundRegion *)malloc(bi->numFiles * sizeof(BoundRegion));
  if (!bi->regions)
    return 4;
  for (i = 0; i < bi->numFiles; i++) {
    if (fgetline(fp, line, MAXLINE) <= 0)
      return 2;
    bi->regions[i].file = strdup(line);
    if (!bi->regions[i].file)
      return 4;
    if (fgetline(fp, line, MAXLINE) <= 0)
      return 2;
    sscanf(line, "%d %d %d %d", &bi->regions[i].section[0], 
           &bi->regions[i].startLine[0], &bi->regions[i].section[1], 
           &bi->regions[i].startLine[1]);

    /* Replace a -1 for second start line with actual line # */
    if (bi->regions[i].section[1] >= 0 && bi->regions[i].startLine[1] < 0)
      bi->regions[i].startLine[1] = nyin - bi->numBoundLines;
  }
  fclose(fp);
  bi->ny = nyin;
  sCurInfo = sNumInfos++;
  return 0;
}

/*!
 * Returns properties read from the current boundary info file.  If the file not HDF, it
 * returns [allSec] 1 if chunks are in Y, [linesBound] with the number of lines in each 
 * boundary, [nfiles] with the number of boundary files.  If the file is HDF, the index to
 * the lock file is returned in [allSec], X size in [linesBound], and Y size in [nfiles].
 * These values are all zero if there was no boundary info file at the current index.  
 * Returns 1 for parallel writing not initialized or -1 for an HDF file.
 */
int parWrtProperties(int *allSec, int *linesBound, int *nfiles)
{
  BoundInfo *bi;
  if (!sNumInfos || sCurInfo < 0)
    return 1;
  *linesBound = 0;
  *allSec = *nfiles = 0;
  bi = &sInfos[sCurInfo];
  
  if (bi->hdfIndex >= 0) {
    *linesBound = bi->nx;
    *nfiles = bi->ny;
    *allSec = bi->hdfIndex;
    return -1;
  }
  if (!bi->regions)
    return 0;
  *linesBound = bi->numBoundLines;
  *allSec = bi->everySec;
  *nfiles = bi->numFiles;
  return 0;
}

/*!
 * Fortran wrapper for @parWrtProperties
 */
int parwrtproperties(int *allSec, int *linesBound, int *nfiles)
{
  return parWrtProperties(allSec, linesBound, nfiles);
}

/*!
 * Sets the index of the current boundary info file or file for parallel writing to 
 * [index], numbered from 0.  Nearly every function in this module applies to the current
 * file set by this index, so if there are multiple files open for parallel writing, this
 * function must be called before any operations on a different file from the previous
 * one.  Returns 1 for index out of range.
 */
int parWrtSetCurrent(int index)
{
  if (index < 0 || index >= sNumInfos)
    return 1;
  sCurInfo = index;
  return 0;
}

/*!
 * Fortran wrapper for @parWrtSetCurrent, with [index] numbered from 1.
 */
int parwrtsetcurrent(int *index)
{
  return parWrtSetCurrent(*index - 1);
}

/*!
 * Closes any lock files opened for parallel writing
 */
void parWrtClose()
{
  int ind;
  BoundInfo *bi;
  for (ind = 0; ind < sNumInfos; ind++) {
    bi = &sInfos[sCurInfo];
    if (bi->hdfIndex >= 0)
      b3dCloseLockFile(bi->hdfIndex);
  }
}

/*! Fortran wrapper for @ParWrtClose */
void parwrtclose()
{
  parWrtClose();
}

/*!
 * Writes one Z slice of data at Z = [slice] from the buffer [buf] to file
 * [fout] according to the header in [hdata].  If parallel writing has been 
 * initialized, lines will be written to a boundary file if appropriate.
 * Returns errors from writing the slice with @@mrcfiles.html#mrc_write_slice@ and 
 * also returns other non-zero values from opening the boundary file, writing its header, 
 * or writing to the file.
 */
int parallelWriteSlice(void *buf, FILE *fout, MrcHeader *hdata, int slice)
{
  static MrcHeader hbound;
  static int dsize, csize, linesBound = -1;
  static int sections[2], startLines[2];
  static FILE *fpBound;
  int err, allsec, nfiles, ib;
  char *filename;
  ImodImageFile *iiFile = (ImodImageFile *)fout;
  int parallelHDF = 0;

  if (sCurInfo >= 0 && sInfos[sCurInfo].hdfIndex >= 0) {
    parallelHDF = 1;
    mrc_getdcsize(hdata->mode, &dsize, &csize);
    err = addSegment((char *)buf, dsize * csize * hdata->nx * hdata->ny, slice, -2, 
                     iiFile, hdata, 0);
    if (err >= 0)
      return err;
    fout = iiFile->fp;
  }

  /* printf("parallelWriteSlice writing %d\n", slice);*/
  err = mrc_write_slice(buf, fout, hdata, slice, 'Z');

  if (parallelHDF && parWrtRecloseHDF(iiFile, hdata))
    return 1;
  if (err || parallelHDF)
    return err;

  if (linesBound < 0) {
    if (parWrtProperties(&allsec, &linesBound, &nfiles))
      linesBound = 0;
    if (!linesBound)
      return 0;
    mrc_head_new(&hbound, hdata->nx, linesBound, 2, hdata->mode);
    err = parWrtFindRegion(slice, 0, hdata->ny, &filename, sections, 
                           startLines);
    if (err) {
      b3dError(stdout, "ERROR: sliceWriteParallel - finding parallel writing"
               " region for slice %d (err %d)\n", slice, err);
      return err;
    }
    if (mrc_getdcsize(hdata->mode, &dsize, &csize)){
      b3dError(stdout, "ERROR: sliceWriteParallel - unknown mode.\n");
      return 1;
    }
    imodBackupFile(filename);
    fpBound = fopen(filename, "wb");
    if (!fpBound) {
      b3dError(stdout, "ERROR: sliceWriteParallel - opening boundary file %s"
               "\n", filename);
      return 1;
    }
    if (mrc_head_write(fpBound, &hbound))
      return 1;
  }

  if (!linesBound)
    return 0;
  for (ib = 0; ib < 2; ib++) {
    if (sections[ib] >= 0 && slice == sections[ib]) {
      fseek(fpBound, hbound.headerSize + ib * hbound.nx * linesBound * csize *
            dsize, SEEK_SET);
      filename = (char *)buf;
      filename += hbound.nx * startLines[ib] * csize * dsize;
      err = mrc_write_slice(filename, fpBound, &hbound, ib, 'Z');
      if (err)
        return err;
    }
  }
  return 0;
}

/*!
 * Fortran-callable function to get the parameters of parallel writing region
 * number [regionNum] (numbered from 1) and return the boundary file in
 * [filename] and the sections and starting lines in arrays [sections] and
 * [startLines].  Returns 2 if writing not initialized, or 1 if the region 
 * number is out of bounds.  This is used by Fixboundaries.
 */
int parwrtgetregion(int *regionNum, char *filename, int *sections, 
                    int *startLines, int strlen)
{
  int reg = *regionNum - 1;
  BoundRegion *regions;
  if (!sNumInfos || sCurInfo < 0 || !sInfos[sCurInfo].regions)
    return 2;
  if (reg < 0 || reg >= sInfos[sCurInfo].numFiles)
    return 1;
  regions = sInfos[sCurInfo].regions;
  sections[0] = regions[reg].section[0];
  startLines[0] = regions[reg].startLine[0];
  sections[1] = regions[reg].section[1];
  startLines[1] = regions[reg].startLine[1];
  return c2fString(regions[reg].file, filename, strlen);
}

/*! 
 * Closes the ImodImageFile [iiFile] for an HDF file opened with  @@iimage.html#iiFOpen@
 * and with header data in [hdata].  The file must be HDF and be the current parallel 
 * file.  Returns 1 for error.
 */
int parWrtRecloseHDF(ImodImageFile *iiFile, MrcHeader *hdata)
{
  static double wallSum = 0.;
  double wallStart = wallTime();
  int err;
  if (hdata && mrc_head_write(iiFile->fp, hdata)) {
    b3dError(stdout, "ERROR:parWrtRecloseHDF  - Rewriting header of HDF file\n");
    return 1;
  }
  iiClose(iiFile);
  err = b3dUnlockFile(sInfos[sCurInfo].hdfIndex);
  if (err) {
    b3dError(stdout, "ERROR: parWrtRecloseHDF - Releasing file lock for HDF file"
             " (err %d)\n", err);
    return 1;
  }
  wallSum +=  wallTime() - wallStart;
  /*printf("parWrtRecloseHDF cumulative time  %.4f\n", wallSum); */
  return 0;
}

/*!
 * Writes out data in HDF buffer for the ImodImageFile [iiFile] opened with
 * @@iimage.html#iiFOpen@ with header data in [hdata].  The file must be HDF and be the 
 * current parallel file.  Returns 1 for error.
 */
int parWrtFlushBuffers(ImodImageFile *iiFile, MrcHeader *hdata)
{
  int ierr;
  if (parWrtPrepareHDF(iiFile))
    return 1;
  ierr = writeSegments(iiFile, hdata, 0);
  if (parWrtRecloseHDF(iiFile, hdata) || ierr)
    return 1;
  return 0;
}

/*!
* Initializes the parallel writing routines for use with unit-based I/O.  The
* name of the boundary info file should be in [filename], which can be
* empty.  An available unit for opening image files should be in [iunitBound], the unit
* number for the image file to be written in [iunitOut] and the dimensions of the image 
* file in [nxin], [nyin], [nzin]. Returns the return value from the initialization routine
* @parWrtInitialize
*/
int iiuParWrtInitialize(const char *filename, int iunitBound, int nxIn, int nyIn,
                        int nzIn)
{
  int numFiles, retval;
  if (sNumInfos >= MAX_INFOS) 
    return 5;
  sIunitBound[sNumInfos] = iunitBound;
  sNxFull[sNumInfos] = B3DABS(nxIn);
  sNyFull[sNumInfos] = nyIn;
  sNzFull[sNumInfos] = nzIn;
  sIzCur[sNumInfos] = 0;
  sIyCur[sNumInfos] = 0;
  sLinesBound[sNumInfos] = 0;
  sIfOpen[sNumInfos] = 0;
  retval = parWrtInitialize(filename, nxIn, nyIn);
  if (!retval && strlen(filename) > 0) {
    parWrtProperties(&sIfAllSec[sCurInfo], &sLinesBound[sCurInfo], &numFiles);
  }
  return retval;
}

/*! Fortran wrapper for @@iiuParWrtInitialize@. Also returns 7 for failure to convert 
  Fortran string. */
int parwrtinitialize(char *filename, int *iunitBound, int *nxIn, int *nyIn, int *nzIn,
                     int namelen)
{
  int err;
  char *cstr = f2cString(filename, namelen);
  if (!cstr)
    return 7;
  err = iiuParWrtInitialize(cstr, *iunitBound, *nxIn, *nyIn, *nzIn);
  free(cstr);
  return err;
}

/*!
* Positions unit [iunit] for writing at section [iz], line [iy], numbered from 0.
*/
void parWrtPosn(int iunit, int iz, int iy)
{
  iiuSetPosition(iunit, iz, iy);
  if (sNumInfos <= 0 || sCurInfo < 0)
    return;
  sIzCur[sCurInfo] = iz;
  sIyCur[sCurInfo] = iy;
}

/*! Fortran wrapper for @parWrtPosn */
void parwrtposn(int *iunit, int *iz, int *iy)
{
  parWrtPosn(*iunit, *iz, *iy);
}

/*!
* Writes an entire section from [array] into unit [iunit] at the current
* position and writes boundary lines if appropriate
*/
int parWrtSec(int iunit, char *array)
{
  int ierr, numBytes;
  ImodImageFile *iiFile;
  float *farray = (float *)array;
  int base = sNxFull[sCurInfo] * (sNyFull[sCurInfo] - 1) + 19;

  if (sCurInfo >= 0 && sInfos[sCurInfo].hdfIndex >= 0) {
    iiFile = iiuGetIIFile(iunit);
    numBytes = iiuBufBytesPerPixel(iunit) * sNxFull[sCurInfo] * sNyFull[sCurInfo];
    ierr = addSegment(array, numBytes, sIzCur[sCurInfo], -1, iiFile, NULL, iunit);

    /* It returned 1 for error, 0 for data handled somehow and file reclosed, or -1
       to go on and write */
    if (ierr > 0)
      return 1;
    if (!ierr) {
      advanceSection();
      return 0;
    }
  }

  /*printf("parWrtSec writing %d\n", sIzCur[sCurInfo]);*/
  ierr = iiuWriteSection(iunit, array);
  if (ierr)
    return ierr;
  if (iiuParWrtRecloseHDF(iunit, 1)) {
    advanceSection();
    return 0;
  }
  if (sNumInfos <= 0 || sCurInfo < 0)
    return 0;
  if (sLinesBound[sCurInfo] == 0) 
    return 0;
  pwOpenIfNeeded(sIzCur[sCurInfo], 0, sNyFull[sCurInfo], &ierr);
  if (ierr) {
    printf("\nERROR: parWrtSec - Finding parallel write boundary region sec %d err %d\n"
           , sIzCur[sCurInfo], ierr);
    exit(1);
  }
  if (sIzCur[sCurInfo] == sIzBound[sCurInfo][0]) {
    iiuSetPosition(sIunitBound[sCurInfo], 0, 0);
    /*printf("Writing iz %d  start of sec to bound file 0,0\n", sIzCur[sCurInfo]);*/
    iiuWriteSection(sIunitBound[sCurInfo], array);
  }
  if (sIzCur[sCurInfo] == sIzBound[sCurInfo][1]) {
    /*printf("Writing iz %d  end of sec to bound file 1,0\n", sIzCur[sCurInfo]);*/
    iiuSetPosition(sIunitBound[sCurInfo], 1, 0);
    ierr = iiuWriteSection(sIunitBound[sCurInfo], 
                           &array[iiuBufBytesPerPixel(iunit) * sNxFull[sCurInfo] * 
                                  (sNyFull[sCurInfo] - sLinesBound[sCurInfo])]);
    if (ierr)
      return ierr;
  }
  advanceSection();
  return 0;
}

/*! Fortran wrapper for @parWrtSec */
int parwrtsec(int *iunit, char *array)
{
  return parWrtSec(*iunit, array);
}

/*!
* Writes one line from [array] into unit [iunit] at the current location
* and writes it to a boundary file if appropriate
*/
int parWrtLin(int iunit, char *array)
{
  int ierr, numBytes;
  ImodImageFile *iiFile;

  if (sCurInfo >= 0 && sInfos[sCurInfo].hdfIndex >= 0) {
    iiFile = iiuGetIIFile(iunit);
    numBytes = iiuBufBytesPerPixel(iunit) * sNxFull[sCurInfo];
    ierr = addSegment(array, numBytes, sIzCur[sCurInfo], sIyCur[sCurInfo], iiFile, NULL,
                      iunit);

    /* It returned 1 for error, 0 for data handled somehow and file reclosed, or -1
       to go on and write */
    if (ierr > 0)
      return 1;
    if (!ierr) {
      advanceLine();
      return 0;
    }
  }

  ierr = iiuWriteLines(iunit, array, 1);
  if (ierr)
    return ierr;
  if (iiuParWrtRecloseHDF(iunit, 1)) {
    advanceLine();
    return 0;
  }
  if (sNumInfos <= 0 || sCurInfo < 0) 
    return 0;
  if (sLinesBound[sCurInfo] == 0)
    return 0;
  if (sIfOpen[sCurInfo] == 0) {
    pwOpenIfNeeded(sIzCur[sCurInfo], sIyCur[sCurInfo], 1, &ierr);
    if (ierr) {
      printf("\nERROR: parWrtSec - Finding parallel write boundary region at %d, %d  "
             "err %d\n", sIzCur[sCurInfo], sIyCur[sCurInfo], ierr);
      exit(1);
    }
  }

  /* The find region return modified non-boundary numbers for Y chunks
   * so that these simple tests work with them too. */
  if (sIfAllSec[sCurInfo] != 0) {
    if (sIyCur[sCurInfo] < sIyBound[sCurInfo][0] + sLinesBound[sCurInfo]) {
      iiuSetPosition(sIunitBound[sCurInfo], 2 * sIzCur[sCurInfo], 
                     sIyCur[sCurInfo] - sIyBound[sCurInfo][0]);
      ierr = iiuWriteLines(sIunitBound[sCurInfo], array, 1);
    }
    if (sIyCur[sCurInfo] >= sIyBound[sCurInfo][1]) {
      iiuSetPosition(sIunitBound[sCurInfo], 2 * sIzCur[sCurInfo] + 1, 
                     sIyCur[sCurInfo] - sIyBound[sCurInfo][1]);
      ierr = iiuWriteLines(sIunitBound[sCurInfo], array, 1);
    }
  } else {
    if (sIzCur[sCurInfo] == sIzBound[sCurInfo][0] &&
        sIyCur[sCurInfo] < sIyBound[sCurInfo][0] + sLinesBound[sCurInfo]) {
      iiuSetPosition(sIunitBound[sCurInfo], 0, sIyCur[sCurInfo] - sIyBound[sCurInfo][0]);
      ierr = iiuWriteLines(sIunitBound[sCurInfo], array, 1);
    }
    if (sIzCur[sCurInfo] == sIzBound[sCurInfo][1] && 
        sIyCur[sCurInfo] >= sIyBound[sCurInfo][1]) {
      iiuSetPosition(sIunitBound[sCurInfo], 1, sIyCur[sCurInfo] - sIyBound[sCurInfo][1]);
      ierr = iiuWriteLines(sIunitBound[sCurInfo], array, 1);
    }
  }
  if (ierr)
    return ierr;
  advanceLine();
  return 0;
}

/*! Fortran wrapper for @parWrtLin */
int parwrtlin(int *iunit, char *array)
{
  return parWrtLin(*iunit, array);
}

/*!
 * Closes the ImodImageFile that is open on unit [iunit] if it is an HDF file, and writes
 * the header first if [writeHeader] is non-zero.  Returns 1 if a file was closed, 0 if 
 * not, and exits on error.
 */
int iiuParWrtRecloseHDF(int iunit, int writeHeader)
{
  ImodImageFile *iiFile;
  MrcHeader *hdata;
  if (sNumInfos <= 0 || sCurInfo < 0 || sInfos[sCurInfo].hdfIndex < 0)
    return 0;
  iiFile = iiuGetIIFile(iunit);
  hdata = iiuMrcHeader(iunit, "iiuParWrtRecloseHDF", 1, 0);
  if (parWrtRecloseHDF(iiFile, writeHeader ? hdata : NULL))
    exit(1);
  return 1;
}

/*! Fortran wrapper for @iiuParWrtRecloseHDF */
int iiuparwrtreclosehdf(int *iunit, int *writeHeader)
{
  return iiuParWrtRecloseHDF(*iunit, *writeHeader);
}

/*!
 * Adds a dataset for section 0 to the HDF file open on unit [iunit] without writing any
 * data.  Exits on error.
 */
void iiuWriteDummySecToHDF(int iunit)
{
  ImodImageFile *iiFile;
  char buf[32];
  iiuSyncWithMrcHeader(iunit);
  iiFile = iiuGetIIFile(iunit);
  if (hdfWriteDummySection(iiFile, buf, 0))
    exit(1);
}

/*! Fortran wrapper for @iiuWriteDummySecToHDF */
void iiuwritedummysectohdf(int *iunit)
{
  iiuWriteDummySecToHDF(*iunit);
}

/*! 
 * Writes output buffers for the the current file open for parallel writing, which must 
 * be an HDF file.  Returns 1 if it is not an HDF file or if there is an error writing 
 * data
 */
int iiuParWrtFlushBuffers(int iunit)
{
  int ierr;
  ImodImageFile *iiFile = iiuGetIIFile(iunit);
  if (iiuPrepareWriteHDF(iunit))
    return 1;
  ierr = writeSegments(iiFile, NULL, iunit);
  iiuParWrtRecloseHDF(iunit, 1);
  return ierr;
}

/*! Fortran wrapper for @iiuParWrtFlushBuffers */
int iiuparwrtflushbuffers(int *iunit)
{
  return iiuParWrtFlushBuffers(*iunit);
}

/**************/
/* General static functions */

/*
 * Finds the parallel writing region that contains [nlWrite] lines starting
 * with line [lineNum] on section [secNum] and returns the boundary file in
 * [filename], and the section and starting line of the first and second
 * boundaries in arrays [sections] and [startLines].  Returns 1 if parallel
 * writing not initialized, or 2 if the region is not found.  
 * Fortran wrapper no longer needed, it was in libcfshr/parallelwrite.c.
 */
static int parWrtFindRegion(int secNum, int lineNum, int nlWrite, char **filename, 
                            int *sections, int *startLines)
{
  int i, pastStart, beforeEnd;
  BoundInfo *bi;
  BoundRegion *regions;
  if (!sNumInfos || sCurInfo < 0 || !sInfos[sCurInfo].regions)
    return 1;
  bi = &sInfos[sCurInfo];
  regions = bi->regions;
  for (i = 0; i < bi->numFiles; i++) {
    pastStart = 1;
    beforeEnd = 1;
    if (bi->everySec) {
      if (regions[i].startLine[0] >= 0 &&
          regions[i].startLine[0] > lineNum + nlWrite - 1)
        pastStart = 0;
      if (regions[i].startLine[1] >= 0 &&
          regions[i].startLine[1] + bi->numBoundLines - 1 < lineNum)
        beforeEnd = 0;
    } else {
      if (regions[i].section[0] >= 0 && 
          (regions[i].section[0] > secNum || 
           (regions[i].section[0] == secNum && 
            regions[i].startLine[0] > lineNum + nlWrite - 1)))
        pastStart = 0;
      if (regions[i].section[1] >= 0 && 
          (regions[i].section[1] < secNum ||
           (regions[i].section[1] == secNum &&
            regions[i].startLine[1] + bi->numBoundLines - 1 < lineNum)))
        beforeEnd = 0;
    }
    if (pastStart && beforeEnd) {
      *filename = regions[i].file;
      sections[0] = regions[i].section[0];
      startLines[0] = regions[i].startLine[0];
      sections[1] = regions[i].section[1];
      startLines[1] = regions[i].startLine[1];

      /* Modify boundary lines for Y chunks so that simple test against the
         boundary will work even if there is no boundary */
      if (bi->everySec) {
        if (startLines[0] < 0)
          startLines[0] -= bi->numBoundLines;
        if (startLines[1] < 0)
          startLines[1] = bi->ny + 1;
      }
      return 0;
    }
  }
  return 2;
}

/*
 * Reopens the given ImodImageFile for writing, fully reanalyzing the metadata of 
 * the HDF file.
 */
static int parWrtPrepareHDF(ImodImageFile *iiFile)
{
  static double wallSum = 0.;
  double wallStart = wallTime();
  int err;
  err = b3dLockFile(sInfos[sCurInfo].hdfIndex);
  if (err) {
    b3dError(stdout, "\nERROR: parWrtPrepareHDF - Obtaining file lock for HDF file "
             "(err %d)\n", err);
    return(1);
  }

  iiFile->state = IISTATE_NOTINIT;
  err = iiReopen(iiFile);
  if (err) {
    b3dError(stdout, "\nERROR: parWrtPrepareHDF - Reopening HDF file (err %d)\n", err);
    return(1);
  }
  wallSum +=  wallTime() - wallStart;
  /*printf("parWrtPrepareHDF cumulative time  %.4f\n", wallSum);*/
  return 0;
}

/*
 * Add a section or line as a segment of data to the HDF buffer; startLine should be -2
 * for a full section to be written with mrc_write_slice, -1 for a full section to be
 * written to unit iunit, or a line number for a line of data.  Returns 1 for errors.
 */
static int addSegment(char *buf, int numBytes, int section, int startLine,
                      ImodImageFile *iiFile, MrcHeader *hdata, int iunit)
{
  BufSegment *seg;
  int err;
  if (sBufIndex[sCurInfo] + numBytes > sBufSize) {

    /* Open file now because something is bound to be written either by writeSegments or
       here because data doesn't fit */
    if (startLine < -1)
      err = parWrtPrepareHDF(iiFile);
    else 
      err = iiuPrepareWriteHDF(iunit);
    if (err)
      return 1;
    
    err = writeSegments(iiFile, hdata, iunit);

    /* Return with -1 to let the data be written if it is a section or if it doesn't fit 
       in the buffer */
    if (!err && (startLine < 0 || sBufIndex[sCurInfo] + numBytes > sBufSize))
      return -1;
    
    /* Otherwise reclose the file, return 1 if error */
    if (startLine < -1)
      parWrtRecloseHDF(iiFile, hdata);
    else
      iiuParWrtRecloseHDF(iunit, 1);
    if (err)
      return 1;
  }
  /* Now add a segment to list */
  if (sNumSegments[sCurInfo] >= sMaxSegments[sCurInfo]) {
    err = startLine >= 0 ? 128 : 8;
    B3DREALLOC(sSegments[sCurInfo], BufSegment, sMaxSegments[sCurInfo] + err);
    if (!sSegments[sCurInfo])
      return 1;
    sMaxSegments[sCurInfo] += err;
  }

  /* printf("add segment adding %d %d\n", section, startLine);*/
  seg = &(sSegments[sCurInfo][sNumSegments[sCurInfo]]);
  seg->bufStart = sHdfBufs[sCurInfo] + sBufIndex[sCurInfo];
  seg->numBytes = numBytes;
  seg->section = section;
  seg->startLine = startLine;
  memcpy(seg->bufStart, buf, numBytes);
  sNumSegments[sCurInfo]++;
  sBufIndex[sCurInfo] += numBytes;
  return 0;
}

/*
 * Write out all buffered data segments
 */
static int writeSegments(ImodImageFile *iiFile, MrcHeader *hdata, int iunit)
{
  int ind, err, numLines, startLine, lineSection;
  BufSegment *seg;
  char *linesBuf;
  if (!sNumSegments[sCurInfo])
    return 0;
  numLines = 0;
  
  /*printf("writeSegments to do %d\n", sNumSegments[sCurInfo]);*/
  for (ind = 0; ind < sNumSegments[sCurInfo]; ind++) {
    seg = &(sSegments[sCurInfo][ind]);
    if (seg->startLine < 0) {

      float *farray = (float *)seg->bufStart;
      int base = sNxFull[sCurInfo] * (sNyFull[sCurInfo] - 1) + 19;

      /* It is a section: write any lines that have been accumulated then write the 
         section */
      if (writeBufferedLines(iunit, lineSection, startLine, linesBuf, &numLines))
        return 1;
      /*printf("writeSegments writing %d\n", seg->section);*/
      if (seg->startLine == -1) {
        iiuSetPosition(iunit, seg->section, 0);
        err = iiuWriteSection(iunit, seg->bufStart);
      } else {
        err = iiWriteSection(iiFile, seg->bufStart, seg->section);
      }
      if (err)
        return 1;

    } else {
      
      /* It is a line.  If it is not sequential from the last, write all lines */
      if (numLines && (startLine + numLines != seg->startLine || 
                       lineSection != seg->section)) {
        if (writeBufferedLines(iunit, lineSection, startLine, linesBuf, &numLines))
          return 1;
      }

      /* Initialize line sequence if none yet, then increase number of lines */
      if (!numLines) {
        linesBuf = seg->bufStart;
        startLine = seg->startLine;
        lineSection = seg->section;
      }
      numLines++;
    }
  }
  if (writeBufferedLines(iunit, lineSection, startLine, linesBuf, &numLines))
    return 1;
  clearSegments(sCurInfo);
  return 0;
}

/* Clears out all the store segments of data in an HDF buffer */
static void clearSegments(int infoInd)
{
  sBufIndex[infoInd] = 0;
  sNumSegments[infoInd] = 0;
}

/**********************/
/* Unit I/O based local functions */

/* 
 * Function to determine region and open file if needed (on predefined unit), when
 * writing at section izsec, nlinesWrite lines starting at iyline.
 */
static void pwOpenIfNeeded(int izSec, int iyLine, int nlinesWrite, int *ierr)
{
  int nxyz[3];
  char *filename;
  char *titleCh = "parallel_write: boundary lines";
  float cell[6] = {0., 0., 0., 90., 90., 90.};
  *ierr = 0;
  if (sNumInfos <= 0 || sCurInfo < 0)
    return;

  if (sIfOpen[sCurInfo] != 0 || sLinesBound[sCurInfo] == 0)
    return;
  *ierr = parWrtFindRegion(izSec, iyLine, nlinesWrite, &filename, &sIzBound[sCurInfo][0],
                           &sIyBound[sCurInfo][0]);
  if (*ierr)
    return;
  nxyz[0] = sNxFull[sCurInfo];
  nxyz[1] = sLinesBound[sCurInfo];
  nxyz[2] = 2;
  if (sIfAllSec[sCurInfo] != 0)
    nxyz[2] = 2 * sNzFull[sCurInfo];
  *ierr = iiuOpen(sIunitBound[sCurInfo], filename, "NEW");
  if (*ierr)
    return;
  iiuCreateHeader(sIunitBound[sCurInfo], nxyz, nxyz, 2, nxyz, 0);
  /*printf("Opened boundary file %d %d %d\n", nxyz[0], nxyz[1], nxyz[2]);*/
  cell[0] = nxyz[0];
  cell[1] = nxyz[1];
  cell[2] = nxyz[2];
  iiuAltCell(sIunitBound[sCurInfo], cell);
  *ierr = iiuWriteHeaderStr(sIunitBound[sCurInfo], titleCh, 0, -32000., 32000., 0.);
  sIfOpen[sCurInfo] = 1;
}

/*
 * Reopens the ImodImageFile for the HDF file open on iunit.  Returns 1 for error.
 */
static int iiuPrepareWriteHDF(int iunit)
{
  ImodImageFile *iiFile = iiuGetIIFile(iunit);
  if (parWrtPrepareHDF(iiFile))
    return 1;
  iiuReassignHeaderPtr(iunit);
  iiuSyncWithMrcHeader(iunit);
  return 0;
}

/* Increment current Z and zero out current line */
static void advanceSection()
{
  sIzCur[sCurInfo]++;
  sIyCur[sCurInfo] = 0;
}

/* Increment line number and wrap to next section if at end of section */
static void advanceLine()
{
  sIyCur[sCurInfo]++;
  if (sIyCur[sCurInfo] >= sNyFull[sCurInfo]) {
    sIyCur[sCurInfo] = 0;
    sIzCur[sCurInfo]++;
  }
}

/* Writes a set of numLines lines from the buffer to the HDF file open on iunit */
static int writeBufferedLines(int iunit, int section, int startLine, char *linesBuf,
                               int *numLines)
{
  if (*numLines) {
    iiuSetPosition(iunit, section, startLine);
    if (iiuWriteLines(iunit, linesBuf, *numLines))
      return 1;
    *numLines = 0;
  }
  return 0;
}
