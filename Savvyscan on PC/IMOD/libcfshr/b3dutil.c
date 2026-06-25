/*   b3dutil.c   - utility functions for getting version and copyright, 
 *                   trimming program name, system dependent I/O, min/max, etc.
 *
 *   Copyright (C) 1995-2016 by the Regents of the University of Colorado.
 *
 *   $Id$
 */                                                                           

#include <stdarg.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "imodconfig.h"
#include "b3dutil.h"
#include "mrcslice.h"

#ifndef _WIN32
#include <sys/time.h>
#include <sys/file.h>
#include <unistd.h>
#include <errno.h>
#include <fcntl.h>
#else
#include <Windows.h>
#include <process.h>
#define getpid _getpid
typedef BOOL (WINAPI *LPFN_GLPI)(PSYSTEM_LOGICAL_PROCESSOR_INFORMATION, PDWORD);
#define putenv _putenv
#endif

#ifdef WIN32_BIGFILE
#include <io.h>
#define fileno _fileno
#define read _read
#define write _write
#define lseek _lseeki64
#define off_t __int64
#endif

#ifdef MAC103_BIGFILE
#include <sys/uio.h>
#endif

#ifdef __APPLE__
#include <sys/sysctl.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef F77FUNCAP
#define imodbackupfile IMODBACKUPFILE
#define imodgetenv IMODGETENV
#define imodgetpid IMODGETPID
#define pidtostderr PIDTOSTDERR
#define imodgetstamp IMODGETSTAMP
#define overridewritebytes OVERRIDEWRITEBYTES
#define readbytessigned READBYTESSIGNED
#define writebytessigned WRITEBYTESSIGNED
#define b3dshiftbytes B3DSHIFTBYTES
#define b3dheaderitembytes B3DHEADERITEMBYTES
#define extraisnbytesandflags EXTRAISNBYTESANDFLAGS
#define cputime CPUTIME
#define walltime WALLTIME
#define b3dmillisleep B3DMILLISLEEP
#define numompthreads NUMOMPTHREADS
#define b3dompthreadnum B3DOMPTHREADNUM
#define numberinlist NUMBERINLIST
#define balancedgrouplimits BALANCEDGROUPLIMITS
#define overrideoutputtype OVERRIDEOUTPUTTYPE
#define b3doutputfiletype B3DOUTPUTFILETYPE
#define setoutputtypefromstring SETOUTPUTTYPEFROMSTRING
#define overrideallbigtiff OVERRIDEALLBIGTIFF
#define setnextoutputsize SETNEXTOUTPUTSIZE
#define settiffcompressiontype SETTIFFCOMPRESSIONTYPE
#define set4bitoutputmode SET4BITOUTPUTMODE
#define overrideinvertmrcorigin OVERRIDEINVERTMRCORIGIN
#define b3dran B3DRAN
#define b3drand B3DRAND
#define b3dsrand B3DSRAND
#define anglewithinlimits ANGLEWITHINLIMITS
#define b3dsetlocktimeout B3DSETLOCKTIMEOUT
#define b3dopenlockfile B3DOPENLOCKFILE
#define b3dlockfile B3DLOCKFILE
#define b3dunlockfile B3DUNLOCKFILE
#define b3dcloselockfile B3DCLOSELOCKFILE
#define b3dphysicalmemory B3DPHYSICALMEMORY
#else
#define imodbackupfile imodbackupfile_
#define imodgetenv imodgetenv_
#define imodgetpid imodgetpid_
#define pidtostderr pidtostderr_
#define imodgetstamp imodgetstamp_
#define overridewritebytes overridewritebytes_
#define readbytessigned readbytessigned_
#define writebytessigned writebytessigned_
#define b3dshiftbytes b3dshiftbytes_
#define b3dheaderitembytes b3dheaderitembytes_
#define extraisnbytesandflags extraisnbytesandflags_
#define cputime cputime_
#define walltime walltime_
#define b3dmillisleep b3dmillisleep_
#define numompthreads numompthreads_
#define b3dompthreadnum b3dompthreadnum_
#define numberinlist numberinlist_
#define balancedgrouplimits balancedgrouplimits_
#define overrideoutputtype overrideoutputtype_
#define b3doutputfiletype b3doutputfiletype_
#define setoutputtypefromstring setoutputtypefromstring_
#define overrideallbigtiff overrideallbigtiff_
#define setnextoutputsize setnextoutputsize_
#define settiffcompressiontype settiffcompressiontype_
#define set4bitoutputmode set4bitoutputmode_
#define overrideinvertmrcorigin overrideinvertmrcorigin_
#define b3dran b3dran_
#define b3drand b3drand_
#define b3dsrand b3dsrand_
#define anglewithinlimits anglewithinlimits_
#define b3dsetlocktimeout b3dsetlocktimeout_
#define b3dopenlockfile b3dopenlockfile_
#define b3dlockfile b3dlockfile_
#define b3dunlockfile b3dunlockfile_
#define b3dcloselockfile b3dcloselockfile_
#define b3dphysicalmemory b3dphysicalmemory_
#endif

static int addToArgVector(const char *arg, char ***argVec, int *vecSize, int *numInVec,
                          const char *pattern, int numPrefix);


/* DNM 2/26/03: These need to be printf instead of fprintf(stderr) to not
   crash imod under Windows */
/*! Prints program name provided in [pname], IMOD version and compilation 
  date */
int imodVersion(const char *pname)
{
  if (pname)
    printf("%s Version %s %s %s\n",
           pname, VERSION_NAME, __DATE__, __TIME__);
  return(VERSION);
}

/*! Prints copyright notice */
void imodCopyright(void)
{
  char *uofc =   "Regents of the University of Colorado";
  printf("Copyright (C) %s by the %s\n", COPYRIGHT_YEARS, uofc);
  return;
}

/*! Calls imodVersion and imodCopyright */
void imodUsageHeader(const char *pname)
{
  imodVersion(pname);
  imodCopyright();
}

/*!
 * Returns the value of IMOD_DIR if it is set, or the default install location
 * for the particular OS.  If [assumed] is not NULL, it is returned with 1 if
 * the string is the default install location.
 */
char *IMOD_DIR_or_default(int *assumed)
{
#ifdef _WIN32
  static char *str = "C:\\cygwin\\usr\\local\\IMOD";
#else
#ifdef __APPLE__
  static char *str = "/Applications/IMOD";
#else
  static char *str = "/usr/local/IMOD";
#endif
#endif
  static char *envdir;
  envdir = getenv("IMOD_DIR");
  if (assumed)
    *assumed = envdir != NULL ? 0 : 1;
  if (envdir)
    return envdir;
  return str;
}

/*! Returns a program name stripped of directories and .exe,
 * given the full name (argv\[0\]) in [fullname].
 * It returns a pointer internal to argv\[0\], unless the name ends in .exe,
 * in which case it tries to return a duplicate copy.  Do not free it. */
char *imodProgName(const char *fullname)
{
  char *tail, *tailback, *exe;
  int indexe;
  tail = strrchr(fullname, '/');
  tailback = strrchr(fullname, '\\');
  if (tailback > tail)
    tail = tailback;

  if (!tail)
    return (char *)fullname;
  tail++;
  exe = strstr(tail, ".exe");
  indexe = strlen(tail) - 4;
  if (!exe || exe != tail + indexe)
    return tail;
  exe = strdup(tail);
  if (!exe)
    return tail;
  exe[indexe] = 0x00;
  return exe;
}


/*! Renames an existing file named [filename] to filename~ and
   deletes filename~ first if necessary */
int imodBackupFile(const char *filename)
{
  struct stat buf;
  int err;
  char *backname;

  /* If file does not exist, return */
  if (stat(filename, &buf))
    return 0;

  /* Get backup name */
  backname = (char *)malloc(strlen(filename) + 3);
  if (!backname)
    return -2;
  sprintf(backname, "%s~", filename);

  /* If the backup file exists, try to remove it first (Windows/Intel) */
  if (!stat(backname, &buf) && remove(backname)) {
    free(backname);
    return -1;
  }

  /* finally, rename file */
  err = rename(filename, backname);
  free(backname);
  return err;
}

/*! A fortran wrapper for imodBackupFile */

int imodbackupfile(const char *filename, int strlen)
{
  int err;
  char *cstr = f2cString(filename, strlen);
  if (!cstr)
    return -2;
  err = imodBackupFile(cstr);
  free(cstr);
  return err;
}

/*! A Fortran-callable routine to get an environment variable, [var] and return
  its value in [value].  Returns -1 for error or 1 if the variable is not 
  defined. */
int imodgetenv(char *var, char *value, int varSize, int valueSize)
{
  char *valPtr;
  char *cstr = f2cString(var, varSize);
  if (!cstr)
    return -1;
  valPtr = getenv(cstr);
  free(cstr);
  if (!valPtr)
    return 1;
  return c2fString(valPtr, value, valueSize);
}

/*! Returns the process ID */
int imodGetpid()
{
  return (int)getpid();
}

/*! Fortran-callable routine to get the process ID */
int imodgetpid()
{
  return (int)getpid();
}

/*! Prints the PID to standard error with the prefix expected by Etomo */
void pidToStderr()
{
  fprintf(stderr, "Shell PID: %d\n", (int)getpid());
  fflush(stderr);
}

/*! Fortran wrapper to @pidToStderr */ 
void pidtostderr()
{
  pidToStderr();
}

/*! Fortran-callable function to return the IMOD stamp for MRC files */
int imodgetstamp()
{
  return IMOD_MRC_STAMP;
}

static int sWriteBytesOverride = -1;

/*! 
 * Sets the value to be returned by @writeBytesSigned overriding both the default value
 * and the value of the environment variable WRITE_MODE0_SIGNED.
 */
void overrideWriteBytes(int value)
{
  sWriteBytesOverride = value;
}

/*! Fortran-callable  version of @overrideWriteBytes */
void overridewritebytes(int *value)
{
  sWriteBytesOverride = *value;
}

/*! 
 * Returns 0 if bytes are to be written unsigned, or non-zero if they are to be written,
 * based on the default value, WRITE_SBYTES_DEFAULT, the value of environment variable 
 * WRITE_MODE0_SIGNED, and whether @overrideWriteBytes has been called.
 */
int writeBytesSigned()
{
  int value = WRITE_SBYTES_DEFAULT;
  char *envPtr = getenv(WRITE_SBYTES_ENV_VAR);
  if (sWriteBytesOverride >= 0)
    return sWriteBytesOverride;
  if (envPtr)
    value = atoi(envPtr);
  return value;
}

/*! Fortran wrapper for @writeBytesSigned */
int writebytessigned()
{
  return writeBytesSigned();
}

/*!
 * Returns 1 if bytes are to be read as signed and the file mode in [mode] is 0, otherwise
 * returns 0.  [stamp] is the {imodStamp} header entry, [flags] is the {imodFlags} entry,
 * [dmin] and [dmax] are the minimum and maximum from the header.  If the [stamp] is for
 * an IMOD file, the value of the first bit of [flags] determines the result unless it
 * is overridden by environment variable READ_MODE0_SIGNED having a value < -1 or > 1.
 * Otherwise, the values of [dmin] and [dmax] determine the value unless READ_MODE0_SIGNED
 * is nonzero.
 */
int readBytesSigned(int stamp, int flags, int mode, float dmin, float dmax)
{
  int value, envVal = 0;
  char *envPtr;
  if (mode)
    return 0;
  envPtr = getenv(READ_SBYTES_ENV_VAR);
  if (envPtr)
    envVal = atoi(envPtr);
  if (stamp == IMOD_MRC_STAMP) {

    /* If the IMOD stamp exists, take the flag value but allow variable to override */
    value = flags & MRC_FLAGS_SBYTES;
    if (envVal < -1)
      value = 0;
    if (envVal > 1)
      value = 1;
  } else {
    
    /* Otherwise make decision based on min/max, and override that */
    if (dmin < 0. && dmax < 128.)
      value = 1;
    else if (dmin >=0 && dmax >= 128.)
      value = 0;
    else if (dmin < 0 && dmax >= 128)
      value = (-dmin > dmax - 128.) ? 1 : 0;
    else
      value = 0;    /* dmin >= 0 and dmax < 128.: Could be either */
    if (envVal < 0)
      value = 0;
    if (envVal > 0)
      value = 1;
  }
  return value;
}

/*! Fortran wrapper for @readBytesSigned */
int readbytessigned(int *stamp, int *flags, int *mode, float *dmin, float *dmax)
{
  return readBytesSigned(*stamp, *flags, *mode, *dmin, *dmax);
}

/*!
 * Changes bytes from unsigned to signed or back either in place or while copying to a 
 * new array, provided that [bytesSigned] is nonzero.  If [direction] is >= 0, changes 
 * unsigned bytes in [usbuf] to signed bytes in [buf] by subtracting 128; otherwise it 
 * changes signed bytes in [buf] to unsigned bytes in [usbuf] by adding 128.  The number 
 * of bytes is given by [nx] * [ny].  [usbuf] and  [buf] may be the same.
 */
void b3dShiftBytes(unsigned char *usbuf, char *sbuf, int nx, int ny, int direction,
                       int bytesSigned)
{
  size_t i, nxy;
  if (!bytesSigned)
    return;
  nxy = (size_t)nx * ny;
  if (direction >= 0)
    for (i = 0; i < nxy; i++)
      *sbuf++ = (char)((int)(*usbuf++) - 128);
  else
    for (i = 0; i < nxy; i++)
      *usbuf++ = (unsigned char)((int)(*sbuf++) + 128);
}

/*! Fortran wrapper to @b3dShiftBytes */
void b3dshiftbytes(unsigned char *usbuf, char *sbuf, int *nx, int *ny, int *direction,
                       int *bytesSigned)
{
  b3dShiftBytes(usbuf, sbuf, *nx, *ny, *direction, *bytesSigned);
}

static int sInvertMrcOriginOverride = -1;

/*!
 * Sets the value to be returned by @invertMrcOriginOnOutput overriding both the default 
 * value and the value of the environment variable INVERT_MRC_ORIGIN.
 */
void overrideInvertMrcOrigin(int value)
{
  sInvertMrcOriginOverride = value;
}

/*! Fortran wrapper for @overrideInvertMrcOrigin */
void overrideinvertmrcorigin(int *value)
{
  overrideInvertMrcOrigin(*value);
}

/*!
 * Returns 1 if the origin should be inverted on output to an MRC file or 0 if not,
 * based first on a value set by calling @@overrideInvertMrcOrigin@, if any, then on the 
 * value of environment variable INVERT_MRC_ORIGIN, if set, then on the default for this
 * version of IMOD defined by INVERT_MRC_ORIGIN_DEFAULT.
 */
int invertMrcOriginOnOutput()
{
  int value = INVERT_MRC_ORIGIN_DEFAULT;
  char *envPtr = getenv(INVERT_MRC_ORIGIN_ENV_VAR);
  if (sInvertMrcOriginOverride >= 0)
    return sInvertMrcOriginOverride;
  if (envPtr)
    value = atoi(envPtr);
  return value;
}

static int sOutputTypeOverride = -1;

/*!
 * Sets the output file type, overriding the default and the value specified by 
 * environment variable.  Valid values are OUTPUT_TYPE_MRC (2) and OUTPUT_TYPE_TIFF (1).
 */
void overrideOutputType(int type)
{
  sOutputTypeOverride = type;
}

/*! Fortran wrapper to @overrideOutputType */
void overrideoutputtype(int *type)
{
  overrideOutputType(*type);
}

/*!
 * Returns the current output file type  OUTPUT_TYPE_MRC (2), OUTPUT_TYPE_TIFF (1), or
 * OUTPUT_TYPE_HDF (5) (provided the package was build with HDF).
 * The type is either the default type (which is MRC), or the value of the environment 
 * variable IMOD_OUTPUT_FORMAT if that is set and is 'MRC', 'TIF', 'TIFF', or 'HDF', 
 * overridden by a value set with @@overrideOutputType@.
 */
int b3dOutputFileType()
{
  int type = OUTPUT_TYPE_DEFAULT;
  char *envtype = getenv(OUTPUT_TYPE_ENV_VAR);
  if (envtype) {
    if (!strcmp(envtype, "MRC"))
      type = OUTPUT_TYPE_MRC;
    else if (!strcmp(envtype, "TIFF") || !strcmp(envtype, "TIF"))
      type = OUTPUT_TYPE_TIFF;
#ifndef NO_JPEG_LIB
    else if (!strcmp(envtype, "JPEG") || !strcmp(envtype, "JPG"))
      type = OUTPUT_TYPE_JPEG;
#endif
#ifndef NO_HDF_LIB
    else if (!strcmp(envtype, "HDF"))
      type = OUTPUT_TYPE_HDF;
#endif
  }
  if (sOutputTypeOverride == OUTPUT_TYPE_MRC || sOutputTypeOverride == OUTPUT_TYPE_TIFF ||
      sOutputTypeOverride == OUTPUT_TYPE_HDF || sOutputTypeOverride == OUTPUT_TYPE_JPEG)
    type = sOutputTypeOverride;
  return type;
}

/*! Fortran wrapper to @b3dOutputFileType */
int b3doutputfiletype()
{
  return b3dOutputFileType();
}

/*!
 * Sets output filetype based on a string in [typeStr], which can contain 'MRC', 'TIF', 
 * 'TIFF', 'JPG', 'JPEG', or 'HDF' (or all lower-case equivalents).  Returns 2 for MRC, 
 * 1 for TIFF, 6 for JPEG, -6 for JPEG if there is no JPEG support, 5 for HDF, -5 for HDF 
 * if there is no HDF support, and -1 for other entries
 */
int setOutputTypeFromString(const char *typeStr)
{
  int retval = -1;
  if (!strcmp(typeStr, "MRC") || !strcmp(typeStr, "mrc"))
    retval = OUTPUT_TYPE_MRC;
  else if (!strcmp(typeStr, "TIFF") || !strcmp(typeStr, "TIF") || 
           !strcmp(typeStr, "tiff") || !strcmp(typeStr, "tif"))
    retval = OUTPUT_TYPE_TIFF;
  else if (!strcmp(typeStr, "JPEG") || !strcmp(typeStr, "JPG") || 
           !strcmp(typeStr, "jpeg") || !strcmp(typeStr, "jpg"))
#ifdef NO_JPEG_LIB
    retval = -OUTPUT_TYPE_JPEG;
#else
    retval = OUTPUT_TYPE_JPEG;
#endif
  else if (!strcmp(typeStr, "HDF") || !strcmp(typeStr, "hdf"))
#ifdef NO_HDF_LIB
    retval = -OUTPUT_TYPE_HDF;
#else
    retval = OUTPUT_TYPE_HDF;
#endif
  if (retval > 0)
    overrideOutputType(retval);
  return retval;
}

/*! Fortran wrapper for @setOutputTypeFromString */
int setoutputtypefromstring(const char *str, int strSize)
{
  char *cstr = f2cString(str, strSize);
  int ret;
  if (!cstr)
    return -2;
  ret = setOutputTypeFromString(cstr);
  free(cstr);
  return ret;
}

static int sAllBigTiffOverride = -1;

/*!
 * A [value] of 1 causes all new TIFF files to be opened in large file format with 
 * "w8", overriding both the default setting in the defined macro ALL_BIGTIFF_DEFAULT 
 * and the value of the environment value IMOD_ALL_BIG_TIFF.
 */
void overrideAllBigTiff(int value)
{
  sAllBigTiffOverride = value;
}

/*! Fortran wrapper for @overrideAllBigTiff */
void overrideallbigtiff(int *value)
{
  overrideAllBigTiff(*value);
}

/*!
 * Returns 1 if all new TIFF files should be opened in large file format, based upon the 
 * default in the defined macro ALL_BIGTIFF_DEFAULT, the value of the environment value 
 * IMOD_ALL_BIG_TIFF, and a value set with @@overrideAllBigTiff@. 
 */
int makeAllBigTiff()
{
  int value = ALL_BIGTIFF_DEFAULT;
  char *envPtr = getenv(ALL_BIGTIFF_ENV_VAR);
  if (envPtr)
    value = atoi(envPtr);
  if (sAllBigTiffOverride >= 0)
    value = sAllBigTiffOverride;
  return value;
}

/*!
 * Indicates the output size in [nx], [ny], [nz] and the data mode in [mode] of a file
 * about to be opened; the routine calls @overrideAllBigTiff to force a new TIFF file to
 * have either the large file format or the older format as appropriate.
 */
void setNextOutputSize(int nx, int ny, int nz, int mode)
{
  int bytes, channels;
  if (dataSizeForMode(mode, &bytes, &channels))
    return;
  overrideAllBigTiff(((double)nx * ny) * nz * channels * bytes > 4.0e9 ? 1 : 0);
}

/*! Fortran wrapper for @setNextOutputSize */
void setnextoutputsize(int *nx, int *ny, int *nz, int *mode)
{
  setNextOutputSize(*nx, *ny, *nz, *mode);
}

/*!
 * Sets the compression type for TIFF output to none, LZW, JPEG, or ZIP for [typeIndex]
 * values of 0, 1, 2, or 3, respectively, by setting the environment variable 
 * IMOD_TIFF_COMPRESSION.  If [overrideEnv] is 0, the variable will not be set if it is 
 * already set; if it is 1, the variable will be set unconditionally, overriding the
 * environment value.
 */
int setTiffCompressionType(int typeIndex, int overrideEnv)
{
  static char compressionStr[28];
  char typeMap[4] = {'1', '5', '7', '8'};
  if (typeIndex < 0 || typeIndex >= 4)
    return 1;
  if (!overrideEnv && getenv("IMOD_TIFF_COMPRESSION"))
    return 0;
  sprintf(compressionStr, "IMOD_TIFF_COMPRESSION=%c", typeMap[typeIndex]);
  putenv(compressionStr);
  return 0;
}

/*! Fortran wrapper to @setTiffCompressionType */
int settiffcompressiontype(int *typeIndex, int *overrideEnv)
{
  return setTiffCompressionType(*typeIndex, *overrideEnv);
}

static int sWrite4bitMode = 0;

/*!
 * If [inVal] is non-zero, byte output into an MRC file opened after this call will be 
 * packed into 4 bits and the MRC mode of the file will be 101.
 */
void set4BitOutputMode(int inVal)
{
  sWrite4bitMode = inVal;
}

/*! Fortran wrapper for @set4BitOutputMode */
void set4bitoutputmode(int *inVal)
{
  sWrite4bitMode = *inVal;
}

/*!
 * Returns a non-zero value if an MRC file being opened for writing in mode 0 should be 
 * made mode 101 and written with 4-bit values instead of bytes.
 */
int write4BitModeForBytes()
{
  return sWrite4bitMode;
}

/*! Creates a C string with a copy of a Fortran string described by [str] and 
  [strsize], using [malloc]. */
char *f2cString(const char *str, int strSize)
{
  int i;
  char *newStr;

  /* find last non-blank character */
  for (i = strSize - 1; i >= 0; i--)
    if (str[i] != ' ')
      break;

  newStr = (char *)malloc(i + 2);
  if (!newStr) {
    return NULL;
  }

  /* copy string if non-null, then put terminator at end */
  if (i >= 0)
    strncpy(newStr, str, i + 1);
  newStr[i + 1] = 0x00;
  return newStr;
}

/*! Converts a C string in [cStr] into a Fortran string [fStr] with size
  [fSize]; returns -1 for error if the string will not fit. */
int c2fString(const char *cStr, char *fStr, int fSize)
{
  while (*cStr && fSize > 0) {
    *fStr++ = *cStr++;
    fSize--;
  }

  /* Return error if there is still a non-null character */
  if (*cStr)
    return -1;

  /* Blank-pad */
  while (fSize > 0) {
    *fStr++ = ' ';
    fSize--;
  }
  return 0;
}

/* Simple error processing routines to avoid having libraries print error
   messages themselves */
static int storeError = 0;
static char errorMess[MAX_IMOD_ERROR_STRING] = "";

/*! Stores an error message and may print it as well.  The message is 
  internally printed with vsprintf.  It is printed to [fout] unless 
  b3dSetStoreError has been called with a non-zero value. */
void b3dError(FILE *fout, const char *format, ...)
{
  va_list args;
  va_start(args, format);
  
  vsprintf(errorMess, format, args);
  if (fout == stderr && storeError < 0)
    fprintf(stdout, errorMess);
  else if (fout && storeError <= 0)
    fprintf(fout, errorMess);
  va_end(args);
}

/*! Sets flag to print messages passed by [b3dError] if [ival] is 0 (the 
  default), just to store them internally if [ival] is 1, or to print messages
  destined to stderr to stdout instead if [ival] is -1. */
void b3dSetStoreError(int ival)
{
  storeError = ival;
}

/*! Returns the current error string */
char *b3dGetError()
{
  return &errorMess[0];
}

/* These routines will simply call the standard C routine under Unix, otherwise
   they will get the file descriptor and call the descriptor-based routine, or 
   its underlying equivalent in Window. */

/*! A substitute for fseek that works for large files on all systems. */
int b3dFseek(FILE *fp, int offset, int flag)
{
#if defined(WIN32_BIGFILE) || defined(MAC103_BIGFILE)
  int handle;
  off_t err;
  if (fp == stdin)
    return 0;
  handle = fileno(fp);
  err = lseek(handle, (off_t)offset, flag);
  return (err == -1 ? -1 : 0);
#else
  if (fp == stdin)
    return 0;
  return fseek(fp, offset, flag);
#endif
}

/*! A substitute for fread that works for large files on all systems.  On Windows and Mac 
 * OSX, if the file is stdin, it will call again up to 5 times to try to satisfy the full
 * request. The total bytes read can not be more than 2 GB on either system.  */
size_t b3dFread(void *buf, size_t size, size_t count, FILE *fp)
{
#if defined(WIN32_BIGFILE) || defined(MAC103_BIGFILE)
  size_t totread = 0, ntoread;
#ifdef MAC103_BIGFILE
  ssize_t nread;
#else
  int nread;
#endif
  int loop, nloop = 1;
  unsigned char *cbuf = (unsigned char *)buf;
  int handle = fileno(fp);
  ntoread = size * count;
  if (fp == stdin)
    nloop = 5;
  for (loop = 0; loop < nloop; loop++) {
    nread = read(handle, (void *)cbuf, ntoread);
    if (nread < 0) {
      nread = 0;
      break;
    }
    totread += nread;
    if (nread >= ntoread)
      break;
    ntoread -= nread;
    cbuf += nread;
  }
  return totread / size;
#else
  return fread(buf, size, count, fp);
#endif
}
 
/*! A substitute for fwrite that works for large files on all systems. */
size_t b3dFwrite(void *buf, size_t size, size_t count, FILE *fp)
{
#if defined(WIN32_BIGFILE) || defined(MAC103_BIGFILE)
  int handle = fileno(fp);
  return (size_t)(write(handle, buf, size * count) / size);
#else
  return fwrite(buf, size, count, fp);
#endif
}

/*! A substitute for rewind that works for large files on all systems. */
void b3dRewind(FILE *fp)
{
  b3dFseek(fp, 0, SEEK_SET);
}

#define SEEK_LIMIT 2000000000

/*!
 * Does a seek in a large file with pointer [fp].  The amount to seek is
 * given by [base] + [size1] * [size2].  [flag] has the standard meaning for
 * seeks, e.g., SEEK_SET, etc.  Returns the nonzero seek error if error.
 */
int mrc_big_seek(FILE *fp, int base, int size1, int size2, int flag)
{
#ifdef USE_SYSTEM_FSEEK

  /* On Mac, need to seek with lseek as long as read/write are going to
     be used for the I/O instead of fread/fwrite */
#if defined(WIN32_BIGFILE) || defined(MAC103_BIGFILE)
  int handle;
  off_t err;
  if (fp == stdin)
    return 0;
  handle = fileno(fp);
  err = lseek(handle, (off_t)size1 * (off_t)size2 + base, flag);
  return (err == -1 ? -1 : 0);
#else
#ifdef _WIN32
  return _fseeki64(fp, (__int64)size1 * (__int64)size2 + base, flag);
#else
  return fseek(fp, (long)size1 * (long)size2 + base, flag);
#endif
#endif
#else
  int smaller, bigger, ntodo, ndo, abs1, abs2;
  int steplimit, err;

  /* Do the base seek if it is non-zero, or if the rest of the seek is
     zero and we are doing a SEEK_SET */
  if (base || ((!size1 || !size2) && (flag == SEEK_SET))) {
    if ((err = b3dFseek(fp, base, flag)))
      return err;
    flag = SEEK_CUR;
  }

  if (!size1 || !size2)
    return 0;

  /* Find smaller and larger size */
  abs1 = size1 >= 0 ? size1 : -size1;
  abs2 = size2 >= 0 ? size2 : -size2;
  smaller = abs1 < abs2 ? abs1 : abs2;
  bigger = abs1 < abs2 ? abs2 : abs1;

  /* Step by multiples of the larger size, but not by more than the limit */
  steplimit = SEEK_LIMIT / bigger;
  ntodo = smaller;

  /* If one of the size entries is negative, negate the steps */
  if ((size1 < 0 && size2 >= 0) || (size1 >= 0 && size2 < 0))
    bigger = -bigger;

  while (ntodo > 0) {
    ndo = ntodo <= steplimit ? ntodo : steplimit;
    if ((err = b3dFseek(fp, ndo * bigger, flag)))
      return err;
    ntodo -= ndo;
    flag = SEEK_CUR;
  }
  return 0;
#endif
}

/*!
 * Does a seek in a large image file with pointer [fp].  The amount to seek is
 * given by [base] plus the change in position indicated by [x], [y], and
 * [z], where [nx] and [ny] are the image dimensions in X and Y and [dsize] is
 * the number of bytes per data element.  Specifically, it seeks by 
 * ^  [base] + [x] * [dsize] + [nx] * [dsize] * ([y] + [ny] * [z]).
 * ^[flag] has the standard meaning for
 * seeks, e.g., SEEK_SET, etc.  Returns the nonzero seek error if error.
 */
int mrcHugeSeek(FILE *fp, int base, int x, int y, int z, int nx, int ny, 
               int dsize, int flag)
{
  int size1, size2;
  double testSize = (double)nx * ny * dsize;

  /* Always add the X offset to the base; then test whether X * Y is OK.
   If so add the Y offset to the base and use X*Y and Z as sizes, otherwise
   use X and Y*Z as the sizes */
  base += x * dsize;
  if (fabs(testSize) < 2.e9) {
    base += nx * y * dsize;
    size1 = nx * ny * dsize;
    size2 = z;
  } else {
    size1 = nx * dsize;
    size2 = y + ny * z;
  }
  return mrc_big_seek(fp, base, size1, size2, flag);
}

/*!
 * Reads a line of characters from the file pointed to by [fp] and places it 
 * into array [s] of size [limit].  Replaces newline or return-newline with a null or 
 * terminates the string with null if reading stops because the array limit is
 * reached.  Returns the length of the string, -1 for an error, -2 for end of 
 * file after a newline, or the negative of the length 
 * of the string plus two for end of file on a line not terminated by newline.
 */
int fgetline(FILE *fp, char s[], int limit)
{
  int c, i, length;

  if (fp == NULL){
    b3dError(stderr, "fgetline: file pointer not valid\n");
    return(-1);
  }

  if (limit < 3){
    b3dError(stderr, "fgetline: limit (%d) must be > 2\n", limit);
    return(-1);
  }
     
  for (i=0; ( ((c = getc(fp)) != EOF) && (i < (limit-1)) && (c != '\n') ); i++)
    s[i]=c;

  /* 1/25/12: Take off a return too! */
  if (i > 0 && s[i-1] == '\r')
    i--;

  /* A \n or EOF on the first character leaves i at 0, so there is nothing
     special to be handled about i being 1, 9/18/09 */
               
  s[i]='\0';
  length = i;

  if (c == EOF)
    return (-1 * (length + 2));
  else
    return (length);
}

/*!
 * For the given MRC file mode or SLICE_MODE_MAX in [mode], returns the number of bytes 
 * of the basic data element in [dataSize] and the number of data channels in [channels].
 * Returns -1 for an unsupported or undefined mode.
 */
int dataSizeForMode(int mode, int *dataSize, int *channels)
{
  switch (mode) {
  case MRC_MODE_BYTE:
    *dataSize = sizeof(b3dUByte);
    *channels = 1;
    break;
  case MRC_MODE_SHORT:
  case MRC_MODE_USHORT:
    *dataSize = sizeof(b3dInt16);
    *channels = 1;
    break;
  case MRC_MODE_FLOAT:
    *dataSize = sizeof(b3dFloat);
    *channels = 1;
    break;
  case MRC_MODE_COMPLEX_SHORT:
    *dataSize = sizeof(b3dInt16);
    *channels = 2;
    break;
  case MRC_MODE_COMPLEX_FLOAT:
    *dataSize = sizeof(b3dFloat);
    *channels = 2;
    break;
  case MRC_MODE_RGB:
    *dataSize = sizeof(b3dUByte);
    *channels = 3;
    break;
  case SLICE_MODE_MAX:
    *dataSize = SLICE_MAX_DSIZE;
    *channels = SLICE_MAX_CSIZE;
    break;
  default:
    return(-1);
  }
  return(0);
}

/*! Returns the number of possible extra header items encoded as short integers
 * by SerialEM in [nflags], and the number of bytes that each occupies in 
 * [nbytes], an array that should be dimensioned to 32. */
void b3dHeaderItemBytes(int *nflags, int *nbytes)
{

  /* Keep this synced to definitions in SerialEM/EMimageExtra.cpp */
  /* 1/1/14: On Mac, icc 11.1, bytes gave 9 2's, 4 and 2 when copied to nbytes if -O2
     and put in a dylib; so make it short. */
  short extra_bytes[] = {2, 6, 4, 2, 2, 4, 2, 4, 2, 4, 2};
  int i;
  *nflags = sizeof(extra_bytes) / sizeof(short);
  for (i = 0; i < *nflags; i++)
    nbytes[i] = extra_bytes[i];
}

/*! A Fortran wrapper for @b3dHeaderItemBytes */
void b3dheaderitembytes(int *nflags, int *nbytes) 
{
  b3dHeaderItemBytes(nflags, nbytes);
}

/*!
 * Returns 1 if the {nint} and {nreal} members of an MRC header represent number of
 * bytes and flags from SerialEM, 0 otherwise.
 */
int extraIsNbytesAndFlags(int nint, int nreal)
{
  int extra_bytes[32];
  int i, extratot = 0, flag_count;
  b3dHeaderItemBytes(&flag_count, &extra_bytes[0]);

  /* DNM 12/10/01: as partial protection against mistaking other entries
     for montage information, at least make sure that the total bytes
     implied by the bits in the flag equals the nint entry. */
  /* DNM 2/3/02: make sure nreal also does not have bits beyond the flags */
  for (i = 0; i < flag_count; i++)
    if (nreal & (1 << i))
      extratot += extra_bytes[i];

  if (nint < 0 || nreal < 0 || nint + nreal == 0 || extratot != nint || 
      nreal >= (1 << flag_count))
    return 0;
  return 1;
}

/*! A Fortran wrapper for @extraIsNbytesAndFlags */
int extraisnbytesandflags(int *nint, int *nreal) 
{
  return extraIsNbytesAndFlags(*nint, *nreal);
}

/*! Set or clear bits in [flags] with the given [mask], depending on whether
 * [state] is nonzero or zero. */
void setOrClearFlags(b3dUInt32 *flags, b3dUInt32 mask, int state)
{
  if (state)
    *flags |= mask;
  else
    *flags &= ~mask;
}

/*! Return 1 if [num] is in the list of [nlist] values in [list], 0 if it is not,
  and [noListValue] if the list is empty ([nlist] 0 or [list] NULL). */
int numberInList(int num, int *list, int nlist, int noListValue)
{
  int i;
  if (!list || !nlist)
    return noListValue;
  for (i = 0; i < nlist; i++)
    if (num == list[i])
      return 1;
  return 0;
}

/*! Fortran wrapper for @numberInList */
int numberinlist(int *num, int *list, int *nlist, int *noListValue)
{
  return numberInList(*num, list, *nlist, *noListValue);
}

/*!
 * Computes inclusive limits [start] and [end] of group at index [groupInd] when dividing 
 * a total of [numTotal] items into [numGroups] groups as evenly as possible, with 
 * the remainder from the division distributed among the first groups.
 */
void balancedGroupLimits(int numTotal, int numGroups, int groupInd, int *start, int *end)
{
  int base = numTotal / numGroups;
  int rem = numTotal % numGroups;
  *start = groupInd * base + B3DMIN(groupInd, rem);
  *end = (groupInd + 1) * base + B3DMIN(groupInd + 1, rem) - 1;
}

/*! Fortran wrapper for @balancedGroupLimits */
void balancedgrouplimits(int *numTotal, int *numGroups, int *groupInd, int *start,
                         int *end)
{
  balancedGroupLimits(*numTotal, *numGroups, *groupInd, start, end);
}

/*!
 * Returns a set of pointers to the lines of [array], which has [xsize] data elements
 * per line, [ysize] lines, and data element size [dsize].  Returns NULL for memory error.
 */
unsigned char **makeLinePointers(void *array, int xsize, int ysize, int dsize)
{
  int i;
  unsigned char **linePtrs = B3DMALLOC(unsigned char *, ysize);
  if (!linePtrs)
    return NULL;
  for (i = 0; i < ysize; i++)
    linePtrs[i] = (unsigned char *)array + (size_t)xsize * i * dsize;
  return linePtrs;
}

/*! A variable argument min function for multiple integer arguments.
 * Call as:  b3dIMin(4, ival1, ival2, ival3, ival4); 
 * For only two arguments with acceptable cost of multiple evaluations,
 * use the macro, B3DMIN(val1, val2); */
/* 2/10/05: gave up on fully generic b3dMin/b3dMax as a bad idea */
int b3dIMin(int narg, ...)
{
  va_list args;
  int extreme;
  int i, val;

  va_start(args, narg);
  for (i = 0; i < narg; i++) {
    val = va_arg(args, int);
    if (!i)
      extreme = val;
    else
      extreme = val < extreme ? val : extreme;
  }
  return(extreme);
}

/*! A variable argument max function for multiple integer arguments.
 * Call as:  b3dIMax(4, ival1, ival2, ival3, ival4); 
 * For only two arguments with acceptable cost of multiple evaluations,
 * use the macro, B3DMAX(val1, val2); */
int b3dIMax(int narg, ...)
{
  va_list args;
  int extreme;
  int i, val;

  va_start(args, narg);
  for (i = 0; i < narg; i++) {
    val = va_arg(args, int);
    if (!i)
      extreme = val;
    else
      extreme = val > extreme ? val : extreme;
  }
  return(extreme);
}

/* To use the high-performance clock on linux, uncomment and link the program
   with -lrt   Otherwise the resolution is only 10 ms */
double cputime(void)
{
  /*#ifdef __linux
  struct timespec currTime;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &currTime);
  return ((double)currTime.tv_sec + currTime.tv_nsec / 1.e9);
  #else */
  return ((double)clock() / CLOCKS_PER_SEC);
  /*#endif*/
}

/*!
 * Returns a measure of time in seconds with microsecond precision on Linux and Mac
 * and the precision of the high performance counter on Windows.
 */
double wallTime(void)
{
#ifdef _WIN32
  LARGE_INTEGER freq, counts;
  QueryPerformanceFrequency(&freq);
  if (!freq.QuadPart)
    return 0.;
  QueryPerformanceCounter(&counts);
  return (((double)counts.QuadPart) / freq.QuadPart);
#else
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return ((double)tv.tv_sec + tv.tv_usec / 1000000.);
#endif
}

/*! Fortran wrapper for @wallTime */
double walltime(void)
{ 
  return wallTime();
}

/*!
 * Sleeps for [msec] milliseconds.  Returns 0 on Windows; elsewhere it returns
 * -1 for an error or the number of times the sleep was interrupted.
 */
int b3dMilliSleep(int msecs)
{
  int retval = 0;
#ifdef _WIN32
  Sleep(msecs);
#else
  struct timespec req, rem;
  req.tv_sec = msecs / 1000;
  req.tv_nsec = 1000000 * (msecs % 1000);
  while (nanosleep(&req, &rem) != 0) {
    if (errno != EINTR)
      return -1;
    req.tv_sec = rem.tv_sec;
    req.tv_nsec = rem.tv_nsec;
    retval++;
  }
#endif
  return retval;
}

/*! Fortran wrapper for @b3dMilliSleep */
int b3dmillisleep(int *msecs)
{ 
  return b3dMilliSleep(*msecs);
}

#define CPUINFO_LINE 80
#define MAX_CPU_SOCKETS 64

/*!
 * Computes the number of threads to specify in an OpenMP directive by taking
 * the minimum of the given optimal thread number [optimalThreads], the number
 * of processors, and the value of OMP_NUM_THREADS, if any.  It evaluates the
 * number of processors and value of OMP_NUM_THREADS only on the first call.  It attempts
 * to limit the number of processors to the number of physical cores by calling 
 * @@int numCoresAndLogicalProcs@.  If the environment variable IMOD_REPORT_CORES is set,
 * it reports whatever values it received from that routine.  
 * Returns 1 if there is no OpenMP available.
 */
int numOMPthreads(int optimalThreads)
{
  int numThreads = optimalThreads;
  int physicalProcs = 0;
  int logicalProcessorCount = 0;
  int processorCoreCount = 0;
#ifdef _OPENMP
  static int limThreads = -1;
  static int numProcs = -1;
  static int forceThreads = -1;
  static int ompNumProcs = -1;
  char *ompNum;

  /* One-time determination of number of physical and logical cores */
  if (numProcs < 0) {
    ompNumProcs = numProcs = omp_get_num_procs();

    /* if there are legal numbers and the logical count is the OMP
     number, set the physical processor count */
    if (!numCoresAndLogicalProcs(&processorCoreCount, &logicalProcessorCount) &&
        processorCoreCount > 0 && logicalProcessorCount == numProcs)
      physicalProcs = processorCoreCount;
    if (getenv("IMOD_REPORT_CORES"))
      printf("core count = %d  logical processors = %d  OMP num = %d => physical "
             "processors = %d\n", processorCoreCount, logicalProcessorCount, numProcs,
             physicalProcs); fflush(stdout);
      
    if (physicalProcs > 0)
      numProcs = B3DMIN(numProcs, physicalProcs);
  }

  /* Limit by number of real cores */
  numThreads = B3DMAX(1, B3DMIN(numProcs, numThreads));

  /* One-time determination of the limit set by OMP_NUM_THREADS */
  if (limThreads < 0) {
    ompNum = getenv("OMP_NUM_THREADS");
    if (ompNum)
      limThreads = atoi(ompNum);
    limThreads = B3DMAX(0, limThreads);
  }

  /* Limit to number set by OMP_NUM_THREADS and to number of real cores */
  if (limThreads > 0)
    numThreads = B3DMIN(limThreads, numThreads);

  /* One-time determination of whether user wants to force a number of threads */
  if (forceThreads < 0) {
    forceThreads = 0;
    ompNum = getenv("IMOD_FORCE_OMP_THREADS");
    if (ompNum) {
      if (!strcmp(ompNum, "ALL_CORES")) {
        if (numProcs > 0)
          forceThreads = numProcs;
      } else if (!strcmp(ompNum, "ALL_HYPER")) {
        if (ompNumProcs > 0)
          forceThreads = ompNumProcs;
      } else {
        forceThreads = atoi(ompNum);
        forceThreads = B3DMAX(0, forceThreads);
      }
    }
  }

  /* Force the number if set */
  if (forceThreads > 0)
    numThreads = forceThreads;

  if (getenv("IMOD_REPORT_CORES"))
      printf("numProcs %d  limThreads %d  numThreads %d\n", numProcs,
             limThreads, numThreads); fflush(stdout);
#else
  numThreads = 1;
#endif
  return numThreads;
}

/*! Fortran wrapper for @numOMPthreads */
int numompthreads(int *optimalThreads)
{
  return numOMPthreads(*optimalThreads);
}

/*!
 * Returns the thread number of the current thread, numbered from 0 when calling from 
 * C or numbered from 1 when calling the Fortran wrapper.
 */
int b3dOMPthreadNum()
{
#ifdef _OPENMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

int b3dompthreadnum()
{
  return b3dOMPthreadNum() + 1;
}

/*!
 * Returns the number of physical processor cores in [physical] and the number of logical
 * processors in [logical].  The return value is 1 if valid information could not be 
 * obtained for both items.  It uses sysctlbyname on Mac, GetLogicalProcessorInformation
 * on Windows, and /proc/cpuinfo on Linux.
 */
int numCoresAndLogicalProcs(int *physical, int *logical)
{
  int processorCoreCount = 0;
  int logicalProcessorCount = 0;

#ifdef __APPLE__
  int temp = 0;
  size_t lenPhys = sizeof(int);
#elif defined(_WIN32)
  LPFN_GLPI glpi;
  PSYSTEM_LOGICAL_PROCESSOR_INFORMATION buffer = NULL;
  PSYSTEM_LOGICAL_PROCESSOR_INFORMATION ptr = NULL;
  DWORD returnLength = 0;
  DWORD byteOffset = 0;
  DWORD numBits = sizeof(ULONG_PTR) * 8;
  ULONG_PTR bitTest;
  DWORD i;
#else
  FILE *fp;
  unsigned char socketFlags[MAX_CPU_SOCKETS];
  char linebuf[CPUINFO_LINE];
  int err, len, curID, curCores;
  char *colon;
#endif

#ifdef __APPLE__
  if (!sysctlbyname("hw.physicalcpu" , &temp, &lenPhys, NULL, 0))
    processorCoreCount = temp;
  lenPhys = 4;
  if (!sysctlbyname("hw.logicalcpu" , &temp, &lenPhys, NULL, 0))
    logicalProcessorCount = temp;

#elif defined(_WIN32)
  /* This is adapted from https://msdn.microsoft.com/en-us/library/ms683194 
     On systems with more than 64 logical processors, the GetLogicalProcessorInformation
     function retrieves logical processor information about processors in the processor
     group to which the calling thread is currently assigned. Use the 
     GetLogicalProcessorInformationEx function to retrieve information about processors
     in all processor groups on the system.  (Windows 7/Server 2008 or above).
  */
  glpi = (LPFN_GLPI)GetProcAddress(GetModuleHandle(TEXT("kernel32")),
                                   "GetLogicalProcessorInformation");
  if (glpi) {
    if (!glpi(buffer, &returnLength) && GetLastError() == ERROR_INSUFFICIENT_BUFFER &&
        returnLength > 0) {
      buffer = (PSYSTEM_LOGICAL_PROCESSOR_INFORMATION)malloc(returnLength);
      if (buffer) {
        if (glpi(buffer, &returnLength)) {
          ptr = buffer;
          while (byteOffset + sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION) <=
                 returnLength) {
            if (ptr->Relationship == RelationProcessorCore) {
              processorCoreCount++;
              bitTest = (ULONG_PTR)1;
              for (i = 0; i < numBits; i++) {
                if (ptr->ProcessorMask & bitTest)
                  logicalProcessorCount++;
                bitTest *= 2;
              }
            }              
            byteOffset += sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION);
            ptr++;
          }
          free(buffer);
        }
      }
    }
  }
#else

  /* Linux: look at /proc/cpuinfo */
  fp = fopen("/proc/cpuinfo", "r");
  if (fp) {
    curID = -1;
    curCores = -1;
    memset(socketFlags, 0, MAX_CPU_SOCKETS);
    while (1) {
      err = 0;
      len = fgetline(fp, linebuf, CPUINFO_LINE);
      if (!len)
        continue;
      if (len == -2)
        break;
      err = 1;
      if (len == -1)
        break;

      /* Look for a "physical id :" and a "cpu cores :" in either order */
      if (strstr(linebuf, "physical id")) {

        /* Error if already got a physical id without cpu cores */
        if (curID >= 0)
          break;
        colon = strchr(linebuf, ':');
        if (colon)
          curID = atoi(colon+1);

        /* Error if no colon or ID out of range */
        if (!colon || curID < 0 || curID >= MAX_CPU_SOCKETS)
          break;
      }
      if (strstr(linebuf, "cpu cores")) {

        /* Error if already got a cpu cores without physical id  */
        if (curCores >= 0)
          break;
        colon = strchr(linebuf, ':');
        if (colon)
          curCores = atoi(colon+1);

        /* Error if no colon or core count illegal */
        if (!colon || curCores <= 0)
          break;
      }

      /* If have both ID and core count, add one logical processor and the number of
         cores the first time this ID is seen to the core count; set ID flag and reset
         the tow numbers */
      if (curID >= 0 && curCores > 0) {
        logicalProcessorCount++;
        if (!socketFlags[curID])
          processorCoreCount += curCores;
        socketFlags[curID] = 1;
        curID = -1;
        curCores = -1;
      }
      err = 0;
      if (len < 0)
        break;
    }
    if (err)
      processorCoreCount *= -1;
    fclose(fp);
  }
#endif
  *physical = processorCoreCount;
  *logical = logicalProcessorCount;
  return (processorCoreCount <= 0 || logicalProcessorCount < 0) ? 1 : 0;
}

/*!
 * Returns the total number of CUDA cores given the compute capability values in [major] 
 * and [minor] and the number of multiprocessors in [multiprocCount].
 */
int totalCudaCores(int major, int minor, int multiprocCount)
{

  /* These values come from helper_cuda_drvapi.h in samples/common/inc for CUDA 6.0 */
  int majMin[] = {0x10, 0x20, 0x21, 0x30, 0x50, -1};
  int cores[] = {8, 32, 48, 192, 128, -1};
  int capable = (major << 4) + minor;
  int index = 0;
  while (majMin[index] > 0) {
    if (capable < majMin[index + 1] || majMin[index + 1] < 0)
      break;
    index++;
  }
  return cores[index] * multiprocCount;
}

/*!
 * Returns the total physical memory in the system in bytes, or 0 if it is not available
 */
double b3dPhysicalMemory()
{
#ifdef __APPLE__
  uint64_t temp = 0;
  size_t lenPhys = sizeof(int);
  if (!sysctlbyname("hw.memsize" , &temp, &lenPhys, NULL, 0))
    return (double)temp;
  return 0;
#elif defined(_WIN32)
  MEMORYSTATUSEX status;
  status.dwLength = sizeof(status);
  GlobalMemoryStatusEx(&status);
  return (double)status.ullTotalPhys;
#else
    long pages = sysconf(_SC_PHYS_PAGES);
    long page_size = sysconf(_SC_PAGE_SIZE);
    if (pages >= 0 && page_size >= 0)
      return (double)pages * page_size;
    return 0;
#endif
}

/*! Fortran wrapper to @b3dPhysicalMemory */
double b3dphysicalmemory()
{
  return b3dPhysicalMemory();
}

/* Function to add one string to a new argument vector, allocating as needed and 
   adding a path prefix if pattern is non-NULL and numPrefix is > 0 */
static int addToArgVector(const char *arg, char ***argVec, int *vecSize, int *numInVec,
                          const char *pattern, int numPrefix)
{
  int quantum = 8;
  if (*numInVec >= *vecSize) {
    if (*vecSize) {
      B3DREALLOC(*argVec, char *, *vecSize + quantum);
    } else {
      *argVec = B3DMALLOC(char *, quantum);
    }
    if (!*argVec)
      return 1;
    *vecSize += quantum;
  }
  if (pattern && numPrefix)
    (*argVec)[*numInVec] = B3DMALLOC(char, numPrefix + strlen(arg) + 1);
  else
    (*argVec)[*numInVec] = strdup(arg);
  if (!(*argVec)[*numInVec])
    return 1;
  if (pattern && numPrefix) {
    strncpy((*argVec)[*numInVec], pattern, numPrefix);
    strcpy((*argVec)[*numInVec] + numPrefix, arg);
  }
  (*numInVec)++;
  return 0;
}

/*!
 * Does wildcard expansion of filenames in the vector of argument strings in [argVec],
 * processing [numArg] strings.  This operation is done only in Windows.  The new number 
 * of arguments is returned in [newNum], and [ifAlloc] is set non-zero if new strings are
 * actually allocated.  In that case, the return value is the new vector; otherwise the 
 * return value is simply [argVec].  If a string fails to match any files, it is passed 
 * unchanged; [noMatchInd] is set to the index of the first such argument.  Returns NULL 
 * for memory allocation errors.
 */
char **expandArgList(const char **argVec, int numArg, int *newNum, int *ifAlloc,
                     int *noMatchInd)
{
#ifdef _WIN32
  int ind, jnd, knd, numStart, lenArg, lenMatch, valid, passNoMatch, numPrefix;
  int outSize = 0;
  char **newVec = NULL;
  WIN32_FIND_DATA FindFileData;
  char *temp;
  HANDLE hFind = NULL;
  *ifAlloc = 0;
  *newNum = 0;
  *noMatchInd = -1;
  
  for (ind = 0; ind < numArg; ind++) {
    if (strchr(argVec[ind], '*') || strchr(argVec[ind], '?')) {

      /* If there is a wild card, look for file(s) */
      passNoMatch = 0;
      hFind = FindFirstFile(argVec[ind], &FindFileData);
      if (hFind == INVALID_HANDLE_VALUE) {
        passNoMatch = 1;

      } else {

        /* Prepare to process this file and additional ones: find a path in front of
         the argument */
        numStart = *newNum;
        lenArg = strlen(argVec[ind]);
        numPrefix = 0;
        jnd = knd = -1;
        temp = strrchr(argVec[ind], '\\');
        if (temp)
          jnd = temp + 1 - argVec[ind];
        temp = strrchr(argVec[ind], '/');
        if (temp)
          knd = temp + 1 - argVec[ind];
        jnd = B3DMAX(jnd, knd);
        if (jnd > 0 && jnd < lenArg)
          numPrefix = jnd;

        /* Loop on files that are found */
        do {

          /* Find first non-matching character if any and make sure it is a ? or *,
           looking after the path prefix */
          lenMatch = strlen(FindFileData.cFileName);
          valid = 1;
          for (jnd = 0; jnd < B3DMIN(lenArg - numPrefix, lenMatch); jnd++) {
            knd = jnd + numPrefix;
            if (FindFileData.cFileName[jnd] != argVec[ind][knd]) {
              if (argVec[ind][knd] != '*' && argVec[ind][knd] != '?')
                valid = 0;
              break;
            }
          }

          /* Skip if not a valid match, otherwise add to arguments */
          if (!valid)
            continue;

          /* If a new list is not allocated yet, get any previous args allocated */
          if (!*ifAlloc) {
            for (jnd = 0; jnd < ind; jnd++)
              if (addToArgVector(argVec[jnd], &newVec, &outSize, newNum, NULL, 0)) {
                B3DFREE(newVec);                
                return NULL;
            }
            *ifAlloc = 1;
            numStart = *newNum;
          }

          if (addToArgVector(FindFileData.cFileName,  &newVec, &outSize, newNum,
                             argVec[ind], numPrefix)) {
            B3DFREE(newVec);                               
            FindClose(hFind);
            return NULL;
          }
        } while (FindNextFile(hFind, &FindFileData) != 0);
        FindClose(hFind);

        if (numStart == *newNum) {
          passNoMatch = 1;
        } else {

          /* If got any, sort the strings */
          for (jnd = numStart; jnd < *newNum - 1; jnd++) {
            for (knd  = jnd + 1; knd < *newNum; knd++) {
              if (strcmp(newVec[jnd], newVec[knd]) > 0) {
                temp = newVec[jnd];
                newVec[jnd] =  newVec[knd];
                newVec[knd] = temp;
              }
            }
          }
        }
      }

      /* If all this gave no match, add to vector if already allocating */
      if (passNoMatch) {
        if (*noMatchInd < 0)
          *noMatchInd = ind;
        if (*ifAlloc) {
          if (addToArgVector(argVec[ind], &newVec, &outSize, newNum, NULL, 0)) {
            B3DFREE(newVec);
            return NULL;
          }
        }
      }

    } else if (*ifAlloc) {
      if (addToArgVector(argVec[ind], &newVec, &outSize, newNum, NULL, 0)) {
        B3DFREE(newVec);
        return NULL;
      }
    }
  }

  if (*ifAlloc)
    return newVec;
  *newNum = numArg;
  return argVec;
  
#else
  *ifAlloc = 0;
  *noMatchInd = -1;
  *newNum = numArg;
  return (char **)argVec;
#endif
}

/*!
 * On Windows, replaces the portion of an argument vector consisting only of filenames
 * with a new list if any of the strings need wild-card expansion.  Pointers to the
 * original argument vector and number of arguments must be passed in [argv] and [argc], 
 * respectively, and [firstInd] should be a pointer to the index of the first argument 
 * after any option entries.  If no expansion is done, these entries are not modified and
 * [ifAlloc] is returned with a 0.  If any filenames are expanded, [argv] is returned with
 * a new vector consisting of just the expanded filenames, [argc] is returned with the 
 * new number of strings, and [ifAlloc] is returned with a 1.  The return value is 0
 * for success, -1 for a memory allocation error, and 1 for failure to match an entry 
 * with wildcards.  Error messages are set with @b3dError in each case.
 */
int replaceFileArgVec(const char ***argv, int *argc, int *firstInd, int *ifAlloc)
{
  int newNum, noMatchInd;
  char **newVec;
  *ifAlloc = 0;
  if (*firstInd >= *argc)
    return 0;
  newVec = expandArgList(*argv + *firstInd, *argc - *firstInd, &newNum, ifAlloc,
                         &noMatchInd);
  if (!newVec) {
    b3dError(stdout, "ERROR: %s - Allocating memory for expanded argument list\n",
             imodProgName((*argv)[0]));
    return -1;
  }
  if (noMatchInd >= 0) {
    B3DFREE(newVec);
    b3dError(stdout, "ERROR: %s - No files match entry %s\n", imodProgName((*argv)[0]),
             (*argv)[noMatchInd + *firstInd]);
    return 1;
  }
  if (*ifAlloc) {
    *argv = (const char **)newVec;
    *argc = newNum;
    *firstInd = 0;
  }
  return 0;
}

/*! 
 * Fortran-callable function to return a random number between 0 and 1 from the
 * standard C library rand() function.
 */
float b3drand()
{
  return ((float)rand()) / RAND_MAX;
}

/*!
 * Fortran-callable function to seed the @b3drand random number generator with the value
 * of [seed].
 */
void b3dsrand(int *seed)
{
  srand(*seed);
}
/*!
 * Fortran-callable function to generate a random number between 0 and 1 using the 
 * standard C library rand() function.  The random number generator is seeded by calling 
 * C library srand() function with [seed] whenever a {different} value of [seed] is
 * passed.  [seed] is never returned with a new value.  It behaves like the Intel Fortran
 * compiler ran() compatibility function in that if it is called repeatedly with [seed],
 * it will continue to generate new numbers in a sequence, unless a new value of seed is
 * passed.
 */ 
float b3dran(int *seed)
{
  static int firstTime = 1, lastSeed = 0;
  if (firstTime || *seed != lastSeed) {
    srand(*seed);
    lastSeed = *seed;
    firstTime = 0;
  }
  return ((float)rand()) / RAND_MAX;
}

/*!
 * Returns [angle] or a value adjusted by a multiple of [upperLim] - [lowerLim] so that
 * [angle] is greater than [lowerLim] and less than or equal to [upperLim].
 */
double angleWithinLimits(float angle, float lowerLim, float upperLim)
{
  float lower = B3DMIN(lowerLim, upperLim);
  float upper = B3DMAX(lowerLim, upperLim);
  float range = upper - lower;
  
  while (angle <= lower)
    angle += range;
  while (angle > upper)
    angle -= range;
  return angle;
}

/*! Fortran wrapper for @angleWithinLimits */
double anglewithinlimits(float *angle, float *lowerLim, float *upperLim)
{
  return angleWithinLimits(*angle, *lowerLim, *upperLim);
}

#define MAX_LOCK_FILES 8
#define NUM_LOCK_BYTES  1024
#ifdef _WIN32
static HANDLE sLockFiles[MAX_LOCK_FILES];
#else
static int sLockFiles[MAX_LOCK_FILES];
#endif
static int sLocksUsed[MAX_LOCK_FILES];
static float sLockTimeouts[MAX_LOCK_FILES];
static int sInitedLocks = 0;
static float sDfltLockTimeout = 30.;

/*!
 * Sets the timeout that will be assigned when a new lock file is opened, in seconds; 
 * the default is 30 sec.
 */
void b3dSetLockTimeout(float timeout)
{
  sDfltLockTimeout = timeout;
}

/*! Fortran wrapper for @b3dSetLockTimeout */
void b3dsetlocktimeout(float *timeout)
{
  sDfltLockTimeout = *timeout;
}

/*!
 * Opens the file whose name is in [filename] for use as a lock file with the following
 * routines.  Up to 8 such files can be opened.  Returns an index to be used in following
 * calls for locking, unlocking, or closing the file, numbered from 0, -1 if too many
 * lock files are open, or -2 if the file could not be opened.
 */
int b3dOpenLockFile(const char *filename)
{
  int ind;
  if (!sInitedLocks)
    for (ind = 0; ind < MAX_LOCK_FILES; ind++)
      sLocksUsed[ind] = -1;
  sInitedLocks = 1;
  for (ind = 0; ind < MAX_LOCK_FILES; ind++)
    if (sLocksUsed[ind] < 0)
      break;
  if (ind >= MAX_LOCK_FILES)
    return -1;
#ifdef _WIN32
  sLockFiles[ind] = CreateFile(filename, GENERIC_READ | GENERIC_WRITE,
                               FILE_SHARE_READ | FILE_SHARE_WRITE, NULL, OPEN_EXISTING,
                               FILE_ATTRIBUTE_NORMAL, NULL);
  if (sLockFiles[ind] == INVALID_HANDLE_VALUE)
    return -2;
#else
  sLockFiles[ind] = open(filename, O_RDWR);
  if (sLockFiles[ind] < 0)
    return -2;
#endif
  sLocksUsed[ind] = 0;
  sLockTimeouts[ind] = sDfltLockTimeout;
  return ind;
}

/*! Fortran wrapper for @@b3dOpenLockFile@.  Returns -3 for failure to convert string */
int b3dopenlockfile(char *filename, int namelen) {
  char *cstr = f2cString(filename, namelen);
  int ret;
  if (!cstr)
    return -3;
  ret = b3dOpenLockFile(cstr);
  free(cstr);
  return ret;
}

/*!
 * Obtains a lock on the file with the given [index], retrying every 50 ms until the
 * timeout defined for this file is reached.  Returns -1 for an index out of range, -2 if
 * no lock file is open with the given index, or 1 if the lock could not be obtained.  If
 * the file is already locked by this process, increments a lock count that is 
 * decremented on each call to @@b3dUnlockFile@.
 */
int b3dLockFile(int index)
{
  double wallStart = wallTime();
#ifndef _WIN32

  /* Use fcntl instead of flock out of superstition that it will be better for NFS on
     older systems */
  struct flock fcntlLock;
  fcntlLock.l_type = F_WRLCK;
  fcntlLock.l_whence = SEEK_SET;
  fcntlLock.l_start = 0;
  fcntlLock.l_len = NUM_LOCK_BYTES;
#endif
  if (index < 0 || index >= MAX_LOCK_FILES)
    return -1;
  if (sLocksUsed[index] < 0)
    return -2;
  if (sLocksUsed[index] > 0) {
    sLocksUsed[index]++;
    return 0;
  }
  while (wallTime() - wallStart < sLockTimeouts[index]) {
#ifdef _WIN32
    if (LockFile(sLockFiles[index], 0, 0, NUM_LOCK_BYTES, 0))
#else
      /*if (!flock(sLockFiles[index], LOCK_EX | LOCK_NB)) */
      if (fcntl(sLockFiles[index], F_SETLK, &fcntlLock) >= 0)
#endif
      {
        sLocksUsed[index]++;
        return 0;
      }
    b3dMilliSleep(50);
  }
  return 1;
}

/*! Fortran wrapper for @b3dLockFile */
int b3dlockfile(int *index) {
  return b3dLockFile(*index);
}

/*!
 * Releases a file lock at the given [index], or just decrements the lock count if it is 
 * greater than 1.  Returns -1 for an index out of range, -2 if no lock file is open 
 * with the given index, -3 if the file is not locked, or 1 if there is an error in the
 * call to unlock the file.
 */
int b3dUnlockFile(int index)
{
#ifndef _WIN32
  struct flock fcntlLock;
  fcntlLock.l_type = F_UNLCK;
  fcntlLock.l_whence = SEEK_SET;
  fcntlLock.l_start = 0;
  fcntlLock.l_len = NUM_LOCK_BYTES;
#endif
  if (index < 0 || index >= MAX_LOCK_FILES)
    return -1;
  if (sLocksUsed[index] < 0)
    return -2;
  if (!sLocksUsed[index])
    return -3;
  if (sLocksUsed[index] == 1) {
#ifdef _WIN32
    if (!UnlockFile(sLockFiles[index], 0, 0, NUM_LOCK_BYTES, 0))
      return 1;
#else
    /*if (flock(sLockFiles[index], LOCK_UN) < 0) */
    if (fcntl(sLockFiles[index], F_SETLK, &fcntlLock) < 0)
      return 1;
#endif
  }
  sLocksUsed[index]--;
  return 0;
}

/*! Fortran wrapper for @b3dUnlockFile */
int b3dunlockfile(int *index) {
  return b3dUnlockFile(*index);
}

/*!
 * Closes the lock file at the given [index].  Returns -1 for an index out of range, -2 
 * if no lock file is open with the given index, or 1 if there is an error in the
 * call to unlock the file.
 */
int b3dCloseLockFile(int index)
{
  if (index < 0 || index >= MAX_LOCK_FILES)
    return -1;
  if (sLocksUsed[index] < 0)
    return -2;
#ifdef _WIN32
  if (!CloseHandle(sLockFiles[index]))
    return 1;
#else
  if (close(sLockFiles[index]) < 0)
    return 1;
#endif
  sLocksUsed[index] = -1;
  return 0;
}

/*! Fortran wrapper for @b3dCloseLockFile */
int b3dcloselockfile(int *index) {
  return b3dCloseLockFile(*index);
}
