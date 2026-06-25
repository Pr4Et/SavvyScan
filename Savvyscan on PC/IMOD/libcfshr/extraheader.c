/*  extraheader.c - Functions for accessing metadata in header, mdoc or autodoc
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 2016 by the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 * $Id$
 */
#include <time.h>
#include "imodconfig.h"
#include "b3dutil.h"
#include "autodoc.h"

#ifdef F77FUNCAP
#define get_extra_header_tilts GET_EXTRA_HEADER_TILTS
#define get_extra_header_items GET_EXTRA_HEADER_ITEMS
#define get_metadata_items GET_METADATA_ITEMS
#define get_extra_header_pieces GET_EXTRA_HEADER_PIECES
#define get_metadata_pieces GET_METADATA_PIECES
#define get_metadata_by_key GET_METADATA_BY_KEY
#define getextraheadervalue GETEXTRAHEADERVALUE
#define getextraheadersecoffset GETEXTRAHEADERSECOFFSET
#define getextraheadermaxsecsize GETEXTRAHEADERMAXSECSIZE
#define getfeiextheadanglescale GETFEIEXTHEADANGLESCALE
#define copyextraheadersection COPYEXTRAHEADERSECTION
#define getmetadataweightingdoses GETMETADATAWEIGHTINGDOSES
#define priordosesfromimagedoses PRIORDOSESFROMIMAGEDOSES
#else
#ifdef G77__HACK
#define get_extra_header_tilts get_extra_header_tilts__
#define get_extra_header_items get_extra_header_items__
#define get_metadata_items get_metadata_items__
#define get_extra_header_pieces get_extra_header_pieces__
#define get_metadata_pieces get_metadata_pieces__
#define get_metadata_by_key get_metadata_by_key__
#else
#define get_extra_header_tilts get_extra_header_tilts_
#define get_extra_header_items get_extra_header_items_
#define get_metadata_items get_metadata_items_
#define get_extra_header_pieces get_extra_header_pieces_
#define get_metadata_pieces get_metadata_pieces_
#define get_metadata_by_key get_metadata_by_key_
#endif
#define getextraheadervalue getextraheadervalue_
#define getextraheadersecoffset getextraheadersecoffset_
#define getextraheadermaxsecsize getextraheadermaxsecsize_
#define getfeiextheadanglescale getfeiextheadanglescale_
#define copyextraheadersection copyextraheadersection_
#define getmetadataweightingdoses getmetadataweightingdoses_
#define priordosesfromimagedoses priordosesfromimagedoses_
#endif

static int extraHeaderSizes(void *extHead, int extSize, int numInt, int numReal,
                            int izSect, int *offset, int *size, int *maxSize);

/*!
 * Returns tilt angles from an extra header written by SerialEM or in the
 * Agard format.  It relies on the Z values of the sections in the file
 * as listed in [izpiece], which should be obtained first with
 * @getExtraHeaderPieces .
 * ^  [array] = array of extra header data
 * ^  [numExtraBytes] = number of bytes of data there
 * ^  [nbytes] = number of bytes per section
 * ^  [iflags] = flags for type of data present
 * ^  [nz] = number of sections or pieces
 * ^  [tilt] = array for tilt angles
 * ^  [numTilts] = # of tilt angles returned (or highest # if there are gaps)
 * ^  [maxTilts] = size of TILT array
 * ^  [izPiece] = Z value of each section in file  ^
 * Returns 1 for an error and sets error string with @@b3dutil.html#b3dError@.
 */
int getExtraHeaderTilts(char *array, int numExtraBytes, int nbytes, int iflags, int nz,
                         float * tilt, int *numTilts, int maxTilts, int *izPiece)

{
  return getExtraHeaderItems(array, numExtraBytes, nbytes, iflags, nz, 1, tilt,
                             tilt, numTilts, maxTilts, izPiece);
}

/*! Fortran wrapper to @@getExtraHeaderTilts@. Exits with an ERROR: string upon error */
void get_extra_header_tilts(char *array, int *numExtraBytes, int *nbytes, int *iflags,
                            int *nz, float *tilt, int *numTilts, int *maxTilts,
                            int *izPiece)
{
  b3dSetStoreError(1);
  if (getExtraHeaderTilts(array, *numExtraBytes, *nbytes, *iflags, *nz,
                          tilt, numTilts, *maxTilts, izPiece)) {
    printf("\nERROR: %s\n", b3dGetError());
    exit(1);
  }
}

/*!
 * Returns values of a defined type from an extra header written by
 * SerialEM; will also return tilt angles from a header in the
 * Agard format.  It relies on the Z values of the sections in the file
 * as listed in [izpiece], which should be obtained first with
 * @getExtraHeaderPieces .
 * ^  [array] = array of extra header data
 * ^  [numExtraBytes] = number of bytes of data there
 * ^  [nbytes] = number of bytes per section
 * ^  [iflags] = flags for type of data present
 * ^  [nz] = number of sections or pieces
 * ^  [itype] = type of data to retrieve: 1 for tilt angle, 3 for stage
 * position, 4 for magnification, 5 for intensity, 6 for exposure dose
 * ^  [val1], [val2] = arrays for one or two values to be returned
 * ^  [numVals] = # of values returned (or highest # if there are gaps)
 * ^  [maxVals] = size of [val1] and [val2] arrays
 * ^  [izPiece] = Z value of each section in file  ^
 * Returns 1 for an error and sets error string with @@b3dutil.html#b3dError@.
 */
int getExtraHeaderItems(char *array, int numExtraBytes, int nbytes, int iflags, int nz,
    int itype, float *val1, float *val2, int *numVals, int maxVals, int *izPiece)
{
  short temp, temp2, *sptr;
  float *fptr;
  double feiScale, angle;
  float ftemp;
  unsigned char btemp;
  int shorts;
  int i, ind, numByteSkip, ival, itemp, numFlags;
  int nbytes_per_item[32];

  *numVals = 0;
  if (numExtraBytes == 0 || (nbytes < 0 && nbytes != -MRC_EXT_TYPE_FEI))
    return 0;

  b3dHeaderItemBytes(&numFlags, &nbytes_per_item[0]);

  shorts = extraIsNbytesAndFlags(nbytes, iflags);
  if (shorts) {
    
    /* if data are packed as shorts, then test for the bit corresponding
     * to itype.  Skip nbyte between sections and advance starting index
     * for each entry prior to the desired one */
    if ((iflags >> (itype - 1)) % 2 == 0 || nbytes == 0)
      return 0;
    numByteSkip = nbytes;
    ind = 0;
    for (i = 1; i <= itype - 1; i++)
      if ((iflags >> (i - 1)) % 2 != 0)
        ind += nbytes_per_item[i - 1];

  } else if (nbytes == -MRC_EXT_TYPE_FEI) {
    if (itype > 1)
      return 0;

    /* Get the scaling for early bug, get the bitmask and test bit 7 for tilt angles */
    feiScale = getFeiExtHeadAngleScale(array);
    if (getExtraHeaderValue(array, 8, 3, &btemp, &temp2, &ival, &ftemp, &angle))
      return 0;
    if ((ival & (1l << 7)) == 0)
      return 0;
    ind = 0;
  } else {
    
    /* otherwise, tilt angle is the first float; need to skip over ints */
    if (iflags == 0 || itype > 1)
      return 0;
    numByteSkip = 4 * (nbytes + iflags);
    ind = 4 * nbytes;
  }
  
  for (i = 0; i < nz; i++) {
    ival = izPiece[i];
    if (ival < 0) {
      b3dError(stdout, "getExtraHeaderItems - Value array not designed for negative "
               "Z values");
      return(1);
    }
    if (ival >= maxVals) {
      b3dError(stdout, "getExtraHeaderItems - Array not big enough for data");
      return(1);
    }
    if (shorts) {
      sptr = (short *)(&array[ind]);
      temp = *sptr;
      if (itype == 1) {
        val1[ival] = temp / 100.;          /* Tilt angle * 100 */
      } else if (itype == 3) {
        val1[ival] = temp / 25.;           /* Stage X and Y * 25. */
        temp = *(sptr + 1);
        val2[ival] = temp / 25.;
      } else if (itype == 4) {
        val1[ival] = temp * 100.;            /* Magnification / 100 */
      } else if (itype == 5) {
        val1[ival] = temp / 25000.;          /* Intensity * 25000. */
      } else if (itype == 6) {
        temp2 = *(sptr + 1);
        val1[ival] = (float)SEMshortsToFloat(temp, temp2); /* Exposure dose */
      }
    } else if (nbytes == -MRC_EXT_TYPE_FEI) {

      /* Get the number of bytes in section for next skip, and get the angle */
      if (getExtraHeaderValue(array, ind, 3, &btemp, &temp2, &numByteSkip, &ftemp, 
                              &angle))
        return 0;
      if (getExtraHeaderValue(array, ind + 100, 4, &btemp, &temp2, &itemp, &ftemp,
                              &angle))
        val1[ival] = 0.;
      else
        val1[ival] = feiScale * angle;
    } else {
      fptr = (float *)(&array[ind]);
      val1[ival] = *fptr;
    }
    *numVals = B3DMAX(*numVals, ival + 1);
    ind = ind + numByteSkip;
    if (ind > numExtraBytes)
      return 0;
  }
  return 0;
}

/*! Fortran wrapper to @@getExtraHeaderItems@. Exits with an ERROR: string upon error */
void get_extra_header_items(char *array, int *numExtraBytes, int *nbytes, int *iflags,
                         int *nz, int *itype, float *val1, float *val2, int *numVals,
                         int *maxVals, int *izPiece)
{
  b3dSetStoreError(1);
  if (getExtraHeaderItems(array, *numExtraBytes, *nbytes, *iflags, *nz,
                          *itype, val1, val2, numVals, *maxVals, izPiece)) {
    printf("\nERROR: %s\n", b3dGetError());
    exit(1);
  }
}

/*!
 * Converts the two short integers [low] and [ihigh] stored by SerialEM
 * for a floating point number back into the number
 */
double SEMshortsToFloat(short low, short ihigh)
{
  int iexp, ival, ivalSign, iexpSign;
  iexpSign = 1;
  ivalSign = 1;
  if (low < 0) {
    low = -low;
    ivalSign = -1;
  }
  if (ihigh < 0) {
    ihigh = -ihigh;
    iexpSign = -1;
  }
  ival = low * 256 + (ihigh % 256);
  iexp = ihigh / 256;
  return (ivalSign * ival * pow(2., iexp * iexpSign));
}


/*!
 * Returns values of a defined type from an image metadata file written by
 * SerialEM or an HDF file in which such data have been incorporated.  It simply calls
 * @getMetadataByKey with a specific key defined by the value of [itypeData].  See
 * @getMetadataByKey for a full description; the arguments special to this call are:
 * ^  [iTypeData] = type of data to retrieve: 1 for tilt angle, 3 for
 * stage position, 4 for magnification, 5 for intensity, 6 for exposure
 * dose, 7 for pixel size, 8 for defocus, 9 for exposure time
 * ^  [val1], [val2] = arrays for one or two values to be returned  ^
 * Returns 1 for an error and sets error string with @@b3dutil.html#b3dError@.
 */
int getMetadataItems(int indAdoc, int iTypeAdoc, int nz, int iTypeData, float *val1,
                      float *val2, int *numVals, int *numFound, int maxVals, int *izPiece)
{
  float val3;
  char **valPtr = NULL;
  char keys[9][20] = {"TiltAngle", "N", "StagePosition", "Magnification", "Intensity",
                      "ExposureDose", "PixelSpacing", "Defocus", "ExposureTime"};
  int iwhich[9] = {2, 1, 3, 1, 2, 2, 2, 2, 2};

  *numVals = 0;
  *numFound = 0;
  if (iTypeData < 1 || iTypeData > 9) {
    b3dError(stdout, "getMetadataItems - type value %d is outside allowed range", 
             iTypeAdoc);
    return 1;
  }

  return getMetadataByKey(indAdoc, iTypeAdoc, nz, keys[iTypeData - 1],
                          iwhich[iTypeData - 1], val1, val2, &val3, valPtr, numVals, 
                          numFound, maxVals, izPiece);
}

/*! Fortran wrapper to @@getMetadataItems@. Exits with an ERROR: string upon error */
void get_metadata_items(int *indAdoc, int *iTypeAdoc, int *nz, int *iTypeData, 
                        float *val1, float *val2, int *numVals, int *numFound, 
                        int *maxVals, int *izPiece)
{
  b3dSetStoreError(1);
  if (getMetadataItems(*indAdoc - 1, *iTypeAdoc, *nz, *iTypeData, val1,
                       val2, numVals, numFound, *maxVals, izPiece)) {
    printf("\nERROR: %s\n", b3dGetError());
    exit(1);
  }
}

/*!
 * Returns values of almost any available type from an image metadata file written by
 * SerialEM or an HDF file in which such data have been incorporated.  The metadata file
 * should be opened first with @@autodoc.html#AdocOpenImageMetadata@, or 
 * @@autodoc.html#AdocGetImageMetaInfo@ should be called on the autodoc of an HDF file.
 * The routine relies on the Z values of the sections in the file as listed in [izPiece],
 * which should be obtained first with @@getMetadataPieces@.
 * ^  [indAdoc] = autodoc index of metadata file
 * ^  [iTypeAdoc] = type of metadata file: 1 for one file, 2 for series, 3 for other 
 * autodoc
 * ^  [nz] = number of sections or pieces
 * ^  [key] = Complete name string for data to retrieve (case-sensitive)
 * ^  [ivalType] = 1 for single integer, 2 for single float, 3 for two floats, 4 for
 * three floats, or 0 for string
 * ^  [val1], [val2], [val3] = arrays for one, two or three values to be returned.  [val2]
 * and [val3] can be single floats if not needed, but when returning strings, [val1] is
 * needed because it will be set to 0 for every index where a string is found.
 * ^  [valString] = array for value strings to be returned.  This can be NULL unless
 * [ivalType] is 0.  Specifically, this is an array of char * pointers; returned strings
 * are duplicated into there and should be freed.
 * ^  [numVals] = # of values returned (or highest # if there are gaps)
 * ^  [numFound] = # of values actually found in metadata.  If the calling program wants
 * to work with an incomplete collection of values, it should initialize the arrays
 * before calling to a value distinct from possible actual values.
 * ^  [maxVals] = size of [val1], [val2], [val3], and [valString] arrays
 * ^  [izPiece] = Z value of each section in file  ^
 * Returns 1 for an error and sets error string with @@b3dutil.html#b3dError@.
 */
int getMetadataByKey(int indAdoc, int iTypeAdoc, int nz, char *key, int ivalType,
                      float *val1, float *val2, float *val3, char **valString,
                      int *numVals, int *numFound, int maxVals, int *izPiece)
{
  int i, ind, ival, itmp;
  char sectNames[3][20] = {ADOC_ZVALUE_NAME, "Image", ADOC_ZVALUE_NAME};
  char *tempChar;
  if (AdocSetCurrent(indAdoc)) {
    b3dError(stdout, "getMetadataByKey - Failed to set autodoc index\n");
    return(1);
  }

  *numVals = 0;
  *numFound = 0;
  if (ivalType == 0) {
    if (!valString) {
      b3dError(stdout, "getMetadataByKey - Requested string items but called with NULL "
             "pointer to string array\n");
      return(1);
    }
    for (i = 0; i < maxVals; i++)
      valString[i] = NULL;
  }

  for (i = 0; i < nz; i++) {
    ival = izPiece[i];
    if (ival < 0) {
      b3dError(stdout, "getMetadataByKey - Value array not designed for negative "
             "Z values\n");
      return(1);
    }
    if (ival >= maxVals) {
      b3dError(stdout, "getMetadataByKey - Array not big enough for data\n");
      return(1);
    }

    ind = i;
    if (iTypeAdoc == 3) {
      ind  = AdocLookupByNameValue(sectNames[iTypeAdoc - 1], i);
      if (ind < 0)
        continue;
    }

    if (ivalType == 0) {
      if (!AdocGetString(sectNames[iTypeAdoc - 1], ind, key, &tempChar)) {
        B3DFREE(valString[ival]);
        valString[ival] = tempChar;
        val1[ival] = 0.;
        *numFound += 1;
      }
    } else if (ivalType == 1) {
      if (!AdocGetInteger(sectNames[iTypeAdoc - 1], ind, key, &itmp)) {
        val1[ival] = itmp;
        *numFound += 1;
      }
    } else if (ivalType == 2) {
      if (!AdocGetFloat(sectNames[iTypeAdoc - 1], ind, key, &val1[ival]))
        *numFound += 1;
    } else if (ivalType == 3) {
      if (!AdocGetTwoFloats(sectNames[iTypeAdoc - 1], ind, key, &val1[ival],&val2[ival]))
        *numFound += 1;
    } else if (ivalType == 4) {
      if (!AdocGetThreeFloats(sectNames[iTypeAdoc - 1], ind, key, &val1[ival],
                              &val2[ival], &val3[ival]))
        *numFound += 1;
    }
    *numVals = B3DMAX(*numVals, ival + 1);
  }
  return 0;
}

/*!
 * Fortran wrapper to @@getMetadataByKey@. Exits with an ERROR: string upon error.
 * [valString] can be a single character variable unless [ivalType] is 0, in which case 
 * it must be an array of character variables 
 */
void get_metadata_by_key(int *indAdoc, int *iTypeAdoc, int *nz, char *key, int *ivalType,
                         float *val1, float *val2, float *val3, char *valString,
                         int *numVals, int *numFound, int *maxVals, int *izPiece, 
                         int keySize, int valSize)
{
  char **newStrings = NULL;
  char *ckey;
  int ind;
  b3dSetStoreError(1);
  *numVals = 0;
  *numFound = 0;
  ckey = f2cString(key, keySize);
  if (ckey && *ivalType == 0) 
    newStrings = B3DMALLOC(char *, *maxVals);
  if (!ckey || (*ivalType == 0 && !newStrings)) {
    printf("\nERROR: get_metadata_by_key - Failed to allocate string arrays\n");
    exit(1);
  }
  if (getMetadataByKey(*indAdoc - 1, *iTypeAdoc, *nz, ckey, *ivalType, val1, val2, val3,
                       newStrings, numVals, numFound, *maxVals, izPiece)) {
    printf("\nERROR: %s\n", b3dGetError());
    exit(1);
  }
  
  /* Ignore strings not long enough... */
  if (*ivalType == 0) {
    for (ind = 0; ind < *numFound; ind++) {
      if (newStrings[ind])
        c2fString(newStrings[ind], &valString[ind * valSize], valSize);
      else
        c2fString(" ", &valString[ind * valSize], valSize);
      B3DFREE(newStrings[ind]);
    }
    free(newStrings);
  }
  free(ckey);
}

/*!
 * Returns piece coordinates from the extra header written by SerialEM
 * ^  [array] = array of extra header data
 * ^  [numExtraBytes] = number of bytes of data there
 * ^  [nbytes] = number of bytes per section
 * ^  [iflags] = flags for type of data present
 * ^  [nz] = number of pieces in the file
 * ^  [ixPiece], [iyPiece], [izPiece] = arrays in which coordinates are returned
 * ^  [numPieces] = number of coordinates returned (should equal [nz])
 * ^  [maxPiece] = size of [piece] arrays  ^
 * Returns 1 for an error and sets error string with @@b3dutil.html#b3dError@.
 */
int getExtraHeaderPieces(char *array, int numExtraBytes, int nbytes, int iflags, int nz,
                          int *ixPiece, int *iyPiece, int *izPiece, int *numPieces,
                          int maxPiece)
{
  unsigned short *sptr;
  int shorts;
  int i, ind;

  *numPieces = 0;
  if (numExtraBytes == 0)
    return 0;
  if (nz > maxPiece) {
    b3dError(stdout, "getExtraHeaderPieces - arrays not large enough for piece lists\n");
    return(1);
  }

  /* if data are packed as shorts, see if the montage flag is set
   * set starting index based on whether there are tilt angles too */
  shorts = extraIsNbytesAndFlags(nbytes, iflags);
  if (!nbytes || !shorts || (iflags / 2) % 2 == 0)
    return 0;
  ind = 0;
  if (iflags % 2 != 0)
    ind = 2;
  for (i = 0; i < nz; i++) {
    if (ind > numExtraBytes)
      return 0;
    sptr = (unsigned short *)(&array[ind]);
    ixPiece[i] = *sptr;
    iyPiece[i] = *(sptr + 1);
    izPiece[i] = *(sptr + 2);
    ind = ind + nbytes;
    *numPieces = i + 1;
  }
  return 0;
}

/*! Fortran wrapper to @@getExtraHeaderPieces@. Exits with an ERROR: string upon error */
void get_extra_header_pieces(char *array, int *numExtraBytes, int *nbytes, int *iflags,
                             int *nz, int *ixPiece, int *iyPiece, int *izPiece,
                             int *numPieces, int *maxPiece)
{
  b3dSetStoreError(1);
  if (getExtraHeaderPieces(array, *numExtraBytes, *nbytes, *iflags, *nz,
                           ixPiece, iyPiece, izPiece, numPieces, *maxPiece)) {
    printf("\nERROR: %s\n", b3dGetError());
    exit(1);
  }
}

/*!
 * Returns piece coordinates from a metadata autodoc file written by
 * SerialEM.  The file should be opened first with
 * @@autodoc.html#AdocOpenImageMetadata@, or information about the autodoc in an HDF file 
 * should be obtained with @@autodoc.html#AdocGetImageMetaInfo@.
 * ^  [indAdoc] = index of autodoc
 * ^  [adocType] = type of metadata file: 1 for one file, 2 for image series, 3 for 
 * autodoc
 * ^  [nz] = number of pieces in the file
 * ^  [ixPiece], [iyPiece], [izPiece] = arrays in which coordinates are
 * returned
 * ^  [maxPiece] = size of [piece] arrays
 * ^  [numFound] = number of sections piece coordinates were found for  ^
 * Returns 1 for an error and sets error string with @@b3dutil.html#b3dError@.
 */
int getMetadataPieces(int indAdoc, int adocType, int nz, int *ixPiece, int *iyPiece,
                       int *izPiece, int maxPiece, int *numFound)
{
  int i, ind;
  char sectNames[3][20] = {ADOC_ZVALUE_NAME, "Image", ADOC_ZVALUE_NAME};
  *numFound = 0;

  if (nz > maxPiece) {
    b3dError(stdout, "getMetadataPieces - Arrays not large enough for piece lists");
    return(1);
  }
  if (AdocSetCurrent(indAdoc)) {
    b3dError(stdout, "get_metadata_pieces - Failed to set autodoc index");
    return(1);
  }
  for (i = 0; i < nz; i++) {
    ind = i;
    if (adocType == 3) {
      ind = AdocLookupByNameValue(sectNames[adocType - 1], i);
      if (ind < 0)
        continue;
    }
    if (!AdocGetThreeIntegers(sectNames[adocType - 1], ind, "PieceCoordinates", 
                              &ixPiece[i], &iyPiece[i], &izPiece[i]))
      *numFound += 1;
  }
  return 0;
}

/*! Fortran wrapper to @@getMetadataPieces@. Exits with an ERROR: string upon error */
void get_metadata_pieces(int *indAdoc, int *adocType, int *nz, int *ixPiece, int *iyPiece,
                         int *izPiece, int *maxPiece, int *numFound)
{
  b3dSetStoreError(1);
  if (getMetadataPieces(*indAdoc - 1, *adocType, *nz, ixPiece, iyPiece,
                        izPiece, *maxPiece, numFound)) {
    printf("\nERROR: %s\n", b3dGetError());
    exit(1);
  }
}

/*!
 * Returns dose information from a metadata autodoc file produced by
 * SerialEM or an HDF file in which such data have been incorporated.  The metadata file
 * should be opened first with @@autodoc.html#AdocOpenImageMetadata@, or 
 * @@autodoc.html#AdocGetImageMetaInfo@ should be called on the autodoc of an HDF file.
 * The data are assumed to be from a tilt series in which dose accumulates; thus a value 
 * for the cumulative dose before each image is returned in [priorDose].  These prior 
 * doses are first sought in "PriorRecordDose" entries.  If those are not available for
 * every image, then the images are assumed to be in the standard order of a 
 * bidirectional tilt series if the value of [bidirNumInvert] is greater than 1, views 1
 * through [bidirNumInvert] are assumed to be in the file in inverted order of 
 * acquisition, while the rest are assumed to be in order.
 * are summed for the images in that order.  If that value is not greater tham 1, then
 * the routine analyzes the DateTime entries to determine the order of images.  If those 
 * are not available it just assumes they are all in order and the return value is -1.
 * The routine relies on the Z values of the sections in the file as listed in [izPiece],
 * which should probably just be filled with 0 to [nz] - 1 to obtain meaningful dose 
 * values for individual pieces of a montage.  Other arguments are:
 * ^  [indAdoc] = autodoc index of metadata file
 * ^  [iTypeAdoc] = type of metadata file: 1 for one file, 2 for series, 3 for other 
 * autodoc
 * ^  [nz] = number of sections or pieces for which data should be returned
 * ^  [izPiece] = Z value of each section in file  ^
 * Returns 1 for an error accessing the autodoc information, or 2 for insuffient 
 * ExposureDose entries, dose entries of 0, DateTime entries for only some sections, or 
 * bad month strings in the DateTime entries.  For any errors, it sets any error string 
 * with @@b3dutil.html#b3dError@.
 */
int getMetadataWeightingDoses(int indAdoc, int iTypeAdoc, int nz, int *izPiece, 
                              int bidirNumInvert, float *priorDose, float *secDose)
{
  float dummy1, dummy2;
  int numVals, numFound, ind, jnd, monInd, retVal = 0;
  char **valStrings;
  time_t *fullTimes;
  int *timeInds;
  char monthBuf[8];
  struct tm tim;
  const char *months[12] = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul",
                            "Aug", "Sep", "Oct", "Nov", "Dec"};

  /* ExposureDose HAS to be there */
  if (getMetadataByKey(indAdoc, iTypeAdoc, nz, "ExposureDose", 2, secDose, &dummy1, 
                       &dummy2, NULL, &numVals, &numFound, nz, izPiece))
    return 1;
  if (numFound < nz) {
    b3dError(stdout, "getMetadataWeightingDoses - %s entries were found in autodoc file "
             "for ExposureDose\n", numFound ? "Not enough" : "No");
    return 2;
  }

  /* And they have to be non-zero */
  for (ind = 0; ind < nz; ind++) {
    if (secDose[ind] <= 0.) {
      b3dError(stdout, "getMetadataWeightingDoses - Some sections have 0 for "
               "ExposureDose in the autodoc file\n");
      return 2;
    }
  }

  /* Look for PriorRecordDose; if they are there we are done 
     Need to ignore incomplete number of these entries thanks to a bug in SerialEM! */
  if (getMetadataByKey(indAdoc, iTypeAdoc, nz, "PriorRecordDose", 2, priorDose, &dummy1, 
                       &dummy2, NULL, &numVals, &numFound, nz, izPiece))
    return 1;
  if (numFound == nz)
    return 0;

  /* If the information about bidirectional is present, use it to compute priorDose in
     the right order and again we are done */
  if (bidirNumInvert > 1) {
    priorDosesFromImageDoses(secDose, nz, bidirNumInvert, priorDose);
    return 0;
  }

  /* Otherwise fall back to getting DateTime strings to determine section order */
  valStrings = B3DMALLOC(char *, nz);
  if (!valStrings) {
    b3dError(stdout, "getMetadataWeightingDoses - Failed to allocate array for string "
             "pointers\n");
    return 1;
  }
  memset(valStrings, 0, sizeof(char *));
  if (getMetadataByKey(indAdoc, iTypeAdoc, nz, "DateTime", 0, priorDose, &dummy1, &dummy2,
                       valStrings, &numVals, &numFound, nz, izPiece)) {
    for (ind = 0; ind < nz; ind++) 
      B3DFREE(valStrings[ind]);
    free(valStrings);
    return 1;
  }

  if (numFound > 0 && numFound < nz) {
    b3dError(stdout, "getMetadataWeightingDoses - The autodoc file has DateTime "
             "entries for only %d of %d sections\n", numFound, nz);
    return 2;
  }

  /* If no time stamps, have to rely on images being in order.  Return -1 in this case */
  if (!numFound) {
    priorDosesFromImageDoses(secDose, nz, 0, priorDose);
    free(valStrings);
    return -1;
  }

  /* Now convert the date-time stamps with sscanf, look up month, fill tm struct and 
   convert to time_t */
  fullTimes = B3DMALLOC(time_t, nz);
  timeInds = B3DMALLOC(int, nz);
  retVal = 0;
  if (!fullTimes || !timeInds) {
    b3dError(stdout, "getMetadataWeightingDoses - Failed to allocate arrays for times\n");
    retVal = 1;
  }
  if (!retVal) {
    for (ind = 0; ind < nz; ind++) {
      sscanf(valStrings[ind], "%d-%3s-%d %d:%d:%d", &tim.tm_mday, monthBuf, &tim.tm_year,
             &tim.tm_hour, &tim.tm_min, &tim.tm_sec);
      tim.tm_year += 100;
      for (monInd = 0; monInd < 12; monInd++)
        if (!strcmp(monthBuf, months[monInd]))
          break;
      if (monInd > 11) {
        b3dError(stdout, "getMetadataWeightingDoses - The DateTime entry %s in the "
                 "autodoc file has an improper month string\n", valStrings[ind]);
        retVal = 2;
        break;
      }
      tim.tm_mon = monInd;
      tim.tm_isdst = -1;
      fullTimes[ind] = mktime(&tim);
      timeInds[ind] = ind;
    }
  }

  /* Sort the times by index 
   difftime gives the time interval FROM the second time TO the first time */
  if (!retVal) {
    for (ind = 0; ind < nz - 1; ind++) {
      for (jnd = ind + 1; jnd < nz; jnd++) {
        if (difftime(fullTimes[timeInds[ind]], fullTimes[timeInds[jnd]]) > 0.) {
          monInd = timeInds[ind];
          timeInds[ind] = timeInds[jnd];
          timeInds[jnd] = monInd;
        }
      }
    }
  
    /* Now compute the prior doses */
    priorDose[timeInds[0]] = 0.;
    for (ind = 1; ind < nz; ind++) 
      priorDose[timeInds[ind]] = priorDose[timeInds[ind - 1]] + 
        secDose[timeInds[ind - 1]];
  }
  
  /* Clean up */
  for (ind = 0; ind < nz - 1; ind++)
    B3DFREE(valStrings[ind]);
  free(valStrings);
  B3DFREE(timeInds);
  B3DFREE(fullTimes);
  return retVal;
}

/*!
 * Fortran wrapper for @getMetadataWeightingDoses
 */
int getmetadataweightingdoses(int *indAdoc, int *iTypeAdoc, int *nz, int *izPiece, 
                              int *bidirNumInvert, float *priorDose, float *secDose)
{
  int err = getMetadataWeightingDoses(*indAdoc - 1, *iTypeAdoc, *nz, izPiece,
                                      *bidirNumInvert, priorDose, secDose);
  if (err > 0) {
    printf("\nERROR: %s\n", b3dGetError());
    exit(1);
  }
  return err;
}

/*!
 * Computes accumulated dose from [nz] image doses in [secDoses] and returns them in 
 * [priorDose].  If [bidirNumInvert] is greater than 1, doses 1 to [bidirNumInvert] are
 * summed in inverted order; otherwise all image doses are summed in order.
 */
void priorDosesFromImageDoses(float *secDose, int nz, int bidirNumInvert,
                              float *priorDose)
{
  int ind;
  if (bidirNumInvert < 1) {
    priorDose[0] = 0.;
    for (ind = 1; ind < nz; ind++) 
      priorDose[ind] = priorDose[ind - 1] + secDose[ind - 1];
  } else {
    priorDose[bidirNumInvert - 1] = 0.;
    for (ind = bidirNumInvert - 2; ind >= 0; ind--) 
      priorDose[ind] = priorDose[ind + 1] + secDose[ind + 1];
    priorDose[bidirNumInvert] = priorDose[0] + secDose[0];
    for (ind = bidirNumInvert + 1; ind < nz; ind++)
      priorDose[ind] = priorDose[ind - 1] + secDose[ind - 1];
  }
}

/*!
 * Fortran wrapper for @priorDosesFromImageDoses
 */
void priordosesfromimagedoses(float *secDose, int *nz, int *bidirNumInvert,
                              float *priorDose)
{
  priorDosesFromImageDoses(secDose, *nz, *bidirNumInvert, priorDose);
}

/*!
 * Gets one value from an extended header in [extHead], starting at the given [offset] and
 * of the type specified by [type].  [type] should 0 for a byte (unsigned char), 1 for a
 * short integer, 2 for a float, 3 for an integer, or 4 for a double (these are 
 * EXT_HEAD_VALUE_ defines).  Values are copied by memcpy, not by casting; the caller is 
 * responsible for casting integers to unsigned if appropriate.  Returns 1 if [type] is 
 * out of range.
 */
int getExtraHeaderValue(void *extHead, int offset, int type, unsigned char *bval,
                       short *sval, int *ival, float *fval, double *dval)
{
  unsigned char *ptr = (unsigned char *)extHead + offset;
  switch (type) {
  case EXT_HEAD_VAL_BYTE:
    *bval = *ptr;
    break;
  case EXT_HEAD_VAL_SHORT:
    memcpy(sval, ptr, 2);
    break;
  case EXT_HEAD_VAL_FLOAT:
    memcpy(fval, ptr, 4);
    break;
  case EXT_HEAD_VAL_INT:
    memcpy(ival, ptr, 4);
    break;
  case EXT_HEAD_VAL_DOUBLE:
    memcpy(dval, ptr, 8);
    break;
  default:
    return 1;
  }
  return 0;    
}

/*! Fortran wrapper for @getExtraHeaderValue */
int getextraheadervalue(void *extHead, int *offset, int *type, unsigned char *bval,
                      short *sval, int *ival, float *fval, double *dval)
{
  return getExtraHeaderValue(extHead, *offset, *type, bval, sval, ival, fval, dval);
}

/* Common routine for finding section offeset, size, and max size */;
static int extraHeaderSizes(void *extHead, int extSize, int numInt, int numReal,
                            int izSect, int *offset, int *size, int *maxSize)
{
  int iz, ival;
  unsigned char *ptr = (unsigned char *)extHead;

  /* SerialEM */
  if (extraIsNbytesAndFlags(numInt, numReal)) {
    *offset = izSect * numInt;
    *size = numInt;
    *maxSize = numInt;
    return 0;
  }

  /* Agard/old FEI */
  if (numInt >= 0 && numReal >= 0) {
    *offset = 4 * izSect  * (numInt + numReal);
    *size = 4 * (numInt + numReal);
    *maxSize = *size;
    return 0;
  }

  /* FEI1 type */
  *maxSize = 0;
  if (numInt == -MRC_EXT_TYPE_FEI) {
    *offset = 0;
    for (iz = 0; iz < izSect; iz++) {
      memcpy(&ival, ptr, 4);
      if (ival < 0 || *offset + ival > extSize)
        return 2;
      *offset += ival;
      ptr += ival;
      ACCUM_MAX(*maxSize, ival);
    }
    memcpy(size, ptr, 4);
    ACCUM_MAX(*maxSize, *size);
    return 0;
  }
  return 1;
}

/*!
 * Returns the location and size of the extended header segment for section [izSect],
 * numbered from 0.  The full extended header should be supplied in [extHead], its size in
 * [extSize], and either true {nint} and {nreal} header members in [numInt] and [numReal]
 * for a SerialEM or Agard/old FEI style header, or the negative of the MRC_EXT_TYPE_ 
 * header type in [numInt] and a possible version number in [numReal].  (These are the
 * values returned by @@unitIO.html#iiuRetExtendedType@.) The byte offset for the section
 * is returned in [offset] and the number of bytes in [size].  The return value is one 
 * for an unknown extended header type or 2 if the header is not long enough to contain
 * data for the given section.
 */
int getExtraHeaderSecOffset(void *extHead, int extSize, int numInt, int numReal,
                            int izSect, int *offset, int *size)
{
  int maxSize;
  return extraHeaderSizes(extHead, extSize, numInt, numReal, izSect, offset, size,
                          &maxSize);
}

/*! Fortran wrapper for @getExtraHeaderSecOffset */
int getextraheadersecoffset(void *extHead, int *extSize, int *numInt, int *numReal, 
                         int *izSect, int *offset, int *size)
{
  return getExtraHeaderSecOffset(extHead, *extSize, *numInt, *numReal, *izSect,
                                 offset, size);
}

/*!
 * Returns the maximum size in bytes of any extended header segment for the first 
 * [numSect] sections into [maxSize].  Other arguments and the return value are as in
 * @@getExtraHeaderSecOffset@.
 */
int getExtraHeaderMaxSecSize(void *extHead, int extSize, int numInt, int numReal, 
                          int numSect, int *maxSize)
{
  int offset, size;
  return extraHeaderSizes(extHead, extSize, numInt, numReal, numSect - 1, &offset, &size,
                          maxSize);
}

/*! Fortran wrapper for @getExtraHeaderMaxSecSize */
int getextraheadermaxsecsize(void *extHead, int *extSize, int *numInt, int *numReal, 
                         int *numSect, int *maxSize)
{
  return getExtraHeaderMaxSecSize(extHead, *extSize, *numInt, *numReal, *numSect,
                               maxSize);
}

/*!
 * Copies the extended header data for the section [izSect] from [extraIn] to [extraOut].
 * The sizes of the two arrays in bytes are in [sizeIn] and [sizeOut], and [numInt] and
 * [numReal] are either the {nint} and {nreal} members or the negative of a MRC_EXT_TYPE_
 * value and a version number.  The cumulative number of bytes copied is maintained in
 * [cumulBytesOut].  Returns 1 for unsupported type, 2 for input data not large enough to
 * contain the section, or 3 for the output array not large enough for the copy.
 */
int copyExtraHeaderSection(void *extraIn, int sizeIn, void *extraOut, int sizeOut,
                           int numInt, int numReal, int izSect, int *cumulBytesOut)
{
  int offset, size, err;
  err = getExtraHeaderSecOffset(extraIn, sizeIn, numInt, numReal, izSect, &offset, &size);
  if (err)
    return err;
  if (size + *cumulBytesOut > sizeOut)
    return 3;
  memcpy((char *)extraOut + *cumulBytesOut, (char *)extraIn + offset, size);
  *cumulBytesOut += size;
  return 0;
}

int copyextraheadersection(void *extraIn, int *sizeIn, void *extraOut, int *sizeOut,
                           int *numInt, int *numReal, int *izSect, int *cumulBytesOut)
{
  return copyExtraHeaderSection(extraIn, *sizeIn, extraOut, *sizeOut, *numInt, *numReal,
                                *izSect, cumulBytesOut);
}

/*!
 * Returns the value by which angles in the FEI1 extended header in [extHead] need to be
 * scaled; this is radians/degree for software versions before a bug was fixed.
 */
double getFeiExtHeadAngleScale(void *extHead)
{
  return strcmp((char *)extHead + 68, "4.4.0.4981") == 0 ? RADIANS_PER_DEGREE : 1.;
}

/*! Fortran wrapper for @getFeiExtHeadAngleScale */
double getfeiextheadanglescale(void *extHead)
{
  return getFeiExtHeadAngleScale(extHead);
}
