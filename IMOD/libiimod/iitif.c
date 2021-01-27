/*
 *    iitif.c    - specific routines for tiff-type ImodImageFile's
 *
 *    Authors:  James Kremer and David Mastronarde
 *
 *   Copyright (C) 1995-2016 by the Regents of the University of Colorado.
 *
 *  $Id$
 */

/************************************************************************** 
This software uses the tiff library which has the following copyright:
Copyright (c) 1988-1996 Sam Leffler
Copyright (c) 1991-1996 Silicon Graphics, Inc.

Permission to use, copy, modify, distribute, and sell this software and
its documentation for any purpose is hereby granted without fee, provided
that (i) the above copyright notices and this permission notice appear in
all copies of the software and related documentation, and (ii) the names of
Sam Leffler and Silicon Graphics may not be used in any advertising or
publicity relating to the software without the specific, prior written
permission of Sam Leffler and Silicon Graphics.

THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND,
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

IN NO EVENT SHALL SAM LEFFLER OR SILICON GRAPHICS BE LIABLE FOR
ANY SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND,
OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
WHETHER OR NOT ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF
LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
OF THIS SOFTWARE.

Additional documentation is at <ftp://ftp.sgi.com/graphics/tiff/doc>
****************************************************************************/

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <errno.h>
#include "imodconfig.h"
#ifdef NOTIFFLIBS
#include "notiffio.h"
#else
#include "tiffio.h"
#endif

#include "iimage.h"
#include "b3dutil.h"

#ifdef _WIN32
#define vsnprintf _vsnprintf
#endif

static int ReadSection(ImodImageFile *inFile, char *buf, int inSection, int convert);
static void copyLine(unsigned char *bdata, unsigned char *obuf, int xout, int convert,
                     int pixsize, int type, int format, int samples, float slope,
                     float offset, int doscale, int unpack4bits, unsigned char *map,
                     unsigned char *first4bitsMap, unsigned char *second4bitsMap,
                     unsigned char *colormap);
static int tiffWriteSectionAny(ImodImageFile *inFile, char *buf, int inSection, 
                               int ifFloat);
static TIFF *openWithoutBMode(ImodImageFile *inFile);
static int setMatchingDirectory(ImodImageFile *inFile, int dirnum);
static void closeWithError(ImodImageFile *inFile, char *message);
static int tiffSyncFromMrcHeader(ImodImageFile *inFile, MrcHeader *hdata);
static void constrainAndStoreMinMax(ImodImageFile *inFile);
typedef void (*TIFFWarningHandler)(const char *module, const char *fmt, va_list ap);

static TIFFWarningHandler oldHandler = NULL;
static void warningHandler(const char *module, const char *fmt, va_list ap);
static int sUseMapping = 2;
static int sWarningsSuppressed = 0;

int iiTIFFCheck(ImodImageFile *inFile)
{
  TIFF* tif;
  FILE *fp;
  b3dUInt16 buf;
  int dirnum = 0;
  uint16 bits, samples, photometric, sampleformat, planarConfig, compression;
  uint16 bitsIm, samplesIm, photoIm, formatIm, planarIm, resUnit;
  uint32 rowsPerStrip;
  uint32 *offsets;
  int nxim, nyim, formatDef, tileWidth, tileLength, gotMin = 0, gotMax = 0;
  int defined, i, j, hasPixelIm, hasPixel = 0, mismatch = 0, err = 0;
  float xResol, yResol, xPixelIm = 0., yPixelIm = 0., xPixel, yPixel, resScale;
  double lastMin, lastMax;
  b3dUInt16 *redp, *greenp, *bluep;
  char *resvar;
  uint16 TVIPStag = 37708;
  float pixelLimit = 3.;
  char *description, *ptr, *blank, *colon;
  int imageJ = 0, imJslices = 0, imJimages = 0;
  float imJmin = 9.e37, imJmax = -9.e37, imJpixel = -1.;
  int formatDefined, hasDateTime = 0;
  int fromSEMCCD = 0;
  RawImageInfo info;

  if (!inFile) 
    return IIERR_BAD_CALL;
  fp = inFile->fp;
  if (!fp)
    return IIERR_BAD_CALL;

  rewind(fp);
  if (fread(&buf, sizeof(b3dUInt16), 1, fp) < 1)
    err = IIERR_IO_ERROR;
  if (!err && (buf != 0x4949) && (buf != 0x4d4d))
    err = IIERR_NOT_FORMAT;
  if (!err && fread(&buf, sizeof(b3dUInt16), 1, fp) < 1)
    err = IIERR_IO_ERROR;
  if (!err && buf != 0x002a && buf != 0x2a00 && buf != 0x002b && buf != 0x2b00)
    err = IIERR_NOT_FORMAT;
  if (err) {
    if (err == IIERR_IO_ERROR)
      b3dError(stderr, "ERROR: iiTIFFCheck - Reading file %s\n", 
               inFile->filename);
    return err;
  }

  /* Close file now, but reopen it if there is a TIFF failure */
  fclose(fp);
  inFile->fp = NULL;
  tif = openWithoutBMode(inFile);
  if (!tif){
    inFile->fp = fopen(inFile->filename, inFile->fmode);
    b3dError(stderr, "ERROR: iiTIFFCheck - Calling TIFFOpen on file %s\n",
             inFile->filename);
    return(IIERR_IO_ERROR);
  }
    
  inFile->header = (char *)tif;
  inFile->nx = inFile->ny = 0;
  inFile->multipleSizes = 0;
  inFile->planesPerImage = 1;
  inFile->contigSamples = 1;
  resvar = getenv("TIFF_RES_PIXEL_LIMIT");
  if (resvar)
    pixelLimit = atof(resvar);

  /* Read each directory of the file, get properties and count usable images */
  do {
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &nxim);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &nyim);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bitsIm);
    TIFFGetField(tif, TIFFTAG_PHOTOMETRIC, &photoIm);

    /* DNM 11/18/01: field need not be defined, set a default */
    defined = TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &samplesIm);
    if (!defined)
      samplesIm = 1;

    TIFFGetField(tif, TIFFTAG_PLANARCONFIG, &planarIm);
    if (!defined)
      photoIm = PLANARCONFIG_CONTIG;

    defined = TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &formatIm);

    /* Handle pixel sizes */
    hasPixelIm = TIFFGetField(tif, TIFFTAG_XRESOLUTION, &xResol);
    resUnit = 0;
    TIFFGetField(tif, TIFFTAG_RESOLUTIONUNIT, &resUnit);
    if (resUnit > 1 && hasPixelIm) {
      if (!TIFFGetField(tif, TIFFTAG_YRESOLUTION, &yResol))
        yResol = xResol;
      resScale = 1.e8 * (resUnit == 2 ? 2.54 : 1.);
      xPixelIm = resScale / xResol;
      yPixelIm = resScale / yResol;
    } else
      hasPixelIm = 0;
    
    if (!dirnum && TIFFGetField(tif, TIFFTAG_DATETIME, &description))
      hasDateTime = 1;

    /* For the first directory, get the description and check if it ImageJ
       If so look for various values that are useful.
       Also check for SerialEMCCD title for packed 4 bit data */
    if (!dirnum && TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &description)) {
      if (strstr(description, "SerialEMCCD") && strstr(description, "4 bits packed")) {
        fromSEMCCD = 1;
        
        /* But if you find it, check for multiple titles from other programs, i.e.,
           a line that has a colon before a blank */
        ptr = description;
        while (strchr(ptr, '\n')) {
          ptr = strchr(ptr, '\n') + 1;
          blank = strchr(ptr, ' ');
          colon = strchr(ptr, ':');
          if (colon && blank && blank - ptr > colon - ptr) {
            fromSEMCCD = 0;
            break;
          }
        }
      }
      if (strstr(description, "ImageJ=") == description) {
        imageJ = 1;
        ptr = strstr(description, "slices=");
        if (ptr)
          imJslices = atoi(ptr + 7);
        ptr = strstr(description, "images=");
        if (ptr)
          imJimages = atoi(ptr + 7);
        ptr = strstr(description, "min=");
        if (ptr)
          imJmin = atoi(ptr + 4);
        ptr = strstr(description, "max=");
        if (ptr)
          imJmax = atoi(ptr + 4);
        ptr = strstr(description, "spacing=");
        if (ptr && strstr(description, "unit=micron"))
          imJpixel = 1.e4 * atof(ptr + 8);
      }
    }

    /* If this is a bigger image, it is a new standard, so set all the
       properties and reset to one directory */
    if ((float)nxim * (float)nyim > (float)inFile->nx * (float)inFile->ny) {
      inFile->nx = nxim;
      inFile->ny = nyim;

      /* Record the strip and tile size for 3dmod caching */
      inFile->tileSizeX = inFile->tileSizeY = 0;
      if (TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsPerStrip)) {
        inFile->tileSizeY = rowsPerStrip;
      } else if (TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tileWidth) &&
                 TIFFGetField(tif, TIFFTAG_TILELENGTH, &tileLength)) {
        inFile->tileSizeX = tileWidth;
        inFile->tileSizeY = tileLength;
      }
      if (TIFFGetField(tif, TIFFTAG_COMPRESSION, &compression))
        inFile->tiffCompression = compression;
      bits = bitsIm;
      photometric = photoIm;
      planarConfig = planarIm;
      samples = samplesIm;
      formatDef = defined;
      sampleformat = formatIm;
      hasPixel = hasPixelIm;
      xPixel = xPixelIm;
      yPixel = yPixelIm;
      
      if (dirnum)
        inFile->multipleSizes = 1;
      dirnum = 1;
    } else if (nxim == inFile->nx && nyim == inFile->ny) {
      dirnum++;

      /* If size matches, check that everything matches */
      if (bitsIm != bits || photoIm != photometric || planarIm != planarConfig
          || samples != samplesIm || defined != formatDef || 
          (defined && formatIm != sampleformat)) {
        mismatch = 1;
        break;
      }
    } else if (dirnum)
      inFile->multipleSizes = 1;
  } while (TIFFReadDirectory(tif));

  /* get the min and max from last directory; if we wrote it, it applies to whole file */
  if (TIFFGetField(tif, TIFFTAG_SMINSAMPLEVALUE, &lastMin))
    gotMin = 1;
  if (TIFFGetField(tif, TIFFTAG_SMAXSAMPLEVALUE, &lastMax))
    gotMax = 1;
  if (!gotMin && imageJ && imJmin < 8.e37) {
    lastMin = imJmin;
    gotMin = 1;
  }
  if (!gotMax && imageJ && imJmax > -8.e37) {
    lastMax = imJmax;
    gotMax = 1;
  }
  if (!hasPixel && imageJ && imJpixel > 0.) {
    xPixel = yPixel = imJpixel;
    hasPixel = 1;
  }

  TIFFSetDirectory(tif, 0);
  formatDefined = TIFFGetField(tif, TIFFTAG_SAMPLEFORMAT, &sampleformat);
  defined = TIFFGetField(tif, TIFFTAG_FILLORDER, &inFile->fillOrder);
  if (!defined)
    inFile->fillOrder = FILLORDER_MSB2LSB;

  /* Don't know how to get the multiple bit entries from libtiff, so can't test
     if they are all 8.  Allow 4 bit if it is not signed and one plane */
  if (mismatch || 
      !((((bits == 4 && !(formatDefined && sampleformat == SAMPLEFORMAT_INT) && 
           samples == 1) || bits == 8 || bits ==16 || bits == 32) && 
         photometric < PHOTOMETRIC_RGB) ||
        (bits == 8 && (photometric == PHOTOMETRIC_RGB || 
                       photometric == PHOTOMETRIC_PALETTE)))) {
    closeWithError(inFile, "ERROR: iiTIFFCheck - Unsupported type of TIFF file\n");
    return(IIERR_NO_SUPPORT);
  }

  /* Recognize K2 files from SerialEMCCD and mark as 4-bit: require either that it has
     the exact dimensions of the full image with half size in X, or that it has the
     description added by SerialEMCCD */
  if (bits == 8 && photometric < PHOTOMETRIC_RGB && samples == 1 && hasDateTime && 
      !(formatDefined && sampleformat == SAMPLEFORMAT_INT) &&
      (sizeCanBe4BitK2SuperRes(inFile->nx, inFile->ny) || fromSEMCCD)) {
    inFile->packed4bits = PACKED_HALF_XSIZE;
    inFile->nx *= 2;
    inFile->fillOrder = FILLORDER_LSB2MSB;
  }

  /* And switch true 4-bit files to 8 bit with the flag set */
  if (bits == 4) {
    inFile->packed4bits = PACKED_4BIT_MODE;
    bits = 8;
  }

  inFile->nz     = dirnum;
  inFile->file   = IIFILE_TIFF;
  inFile->format = IIFORMAT_LUMINANCE;

  /* Samples are assumed to be channels for RGB; record samples appropriately
     and set Z size otherwise for multiple samples */
  inFile->rgbSamples = samples;
  if (photometric < PHOTOMETRIC_RGB) {
    if (planarConfig == PLANARCONFIG_SEPARATE)
      inFile->planesPerImage = samples;
    else
      inFile->contigSamples = samples;
    inFile->nz = dirnum * samples;
  }

  /* 11/22/08: define this for all types, not just for 3-sample data */
  inFile->readSection = tiffReadSection;
  inFile->readSectionUShort = tiffReadSectionUShort;
  inFile->readSectionByte = tiffReadSectionByte;
  inFile->readSectionFloat = tiffReadSectionFloat;

  /* Set up file mode and default properties; fill in mode for raw info at same time */
  if (bits == 8) {
    inFile->type   = IITYPE_UBYTE;
    inFile->amin  = 0;
    inFile->amean  = 128;
    inFile->amax   = 255;
    inFile->mode   = MRC_MODE_BYTE;

    /* Get the max right for 4-bit data, it matters for very low count data */
    if (inFile->packed4bits) {
      inFile->amean = 7.5;
      inFile->amax = 15;
    }
    info.type = RAW_MODE_BYTE;
    if (photometric == PHOTOMETRIC_RGB) {
      inFile->format = IIFORMAT_RGB;
      inFile->mode   = MRC_MODE_RGB;
      info.type = RAW_MODE_RGB;
    } else if (photometric == PHOTOMETRIC_PALETTE) {

      /* For palette images, define as colormap, get the colormap and convert 
         it to bytes */
      inFile->format = IIFORMAT_COLORMAP;
      inFile->mode   = MRC_MODE_RGB;
      info.type = RAW_MODE_RGB;
      inFile->colormap = (unsigned char *)malloc(3 * 256 * dirnum);
      if (!inFile->colormap) {
        closeWithError(inFile, "ERROR: iiTIFFCheck - Getting memory for "
                       "colormap\n");
        return(IIERR_MEMORY_ERR);
      }
      for (j = 0; j < dirnum; j++) {
        if (setMatchingDirectory(inFile, j)) {
          closeWithError(inFile, "ERROR: iiTIFFCheck - getting directory for "
                         "colormap\n");
          return(IIERR_IO_ERROR);
        }

        TIFFGetField(tif, TIFFTAG_COLORMAP, &redp, &greenp, &bluep);
        for (i = 0; i < 256; i++) {
          inFile->colormap[j*768 + i] = (unsigned char)(redp[i] >> 8);
          inFile->colormap[j*768 + i + 256] = (unsigned char)(greenp[i] >> 8);
          inFile->colormap[j*768 + i + 512] = (unsigned char)(bluep[i] >> 8);
        }
      }
      TIFFSetDirectory(tif, 0);
    } else if (formatDefined && sampleformat == SAMPLEFORMAT_INT) {
      inFile->type   = IITYPE_BYTE;
      info.type = RAW_MODE_SBYTE;
    }
  } else {
    /* If there is a field specifying signed numbers, set up for signed;
       otherwise set up for unsigned */
    if (bits == 16) {
      if (formatDefined && sampleformat == SAMPLEFORMAT_INT) {
        inFile->type   = IITYPE_SHORT;
        inFile->amean  = 0;
        inFile->amin   = -32767;
        inFile->amax   = 32767;
        inFile->mode   = MRC_MODE_SHORT;
        info.type      = RAW_MODE_SHORT;
      } else {
        inFile->type   = IITYPE_USHORT;
        inFile->amean  = 32767;
        inFile->amin   = 0;
        inFile->amax   = 65535;
        inFile->mode   = MRC_MODE_USHORT;   /* Why was this SHORT for both? */
        info.type      = RAW_MODE_USHORT;
      }
    } else {

      /* Set up for integer data: until there is an MRC mode, set to -1 */
      if (formatDefined && sampleformat == SAMPLEFORMAT_INT) {
        inFile->type   = IITYPE_INT;
        inFile->amean  = 0;
        inFile->amin   = -65536;
        inFile->amax   = 65536;
        inFile->mode   = -1;
      } else if (formatDefined && sampleformat == SAMPLEFORMAT_UINT) {
        inFile->type   = IITYPE_UINT;
        inFile->amean  = 65536;
        inFile->amin   = 0;
        inFile->amax   = 130000;
        inFile->mode   = -1;
      } else if (formatDefined && sampleformat == SAMPLEFORMAT_IEEEFP){
        inFile->type   = IITYPE_FLOAT;
        inFile->amean  = 128.;
        inFile->amin   = 0;
        inFile->amax   = 255.;
        inFile->mode   = MRC_MODE_FLOAT;
        info.type      = RAW_MODE_FLOAT;
      } else {
        closeWithError(inFile, "ERROR: iiTIFFCheck - 32-bit TIFF "
                       "file with no data type defined\n");
        return(IIERR_NO_SUPPORT);
      }
    }
  }

  if (TIFFGetField(tif, TVIPStag, &bits, &inFile->userData) > 0) {
    inFile->userCount = bits;
    inFile->userFlags = IIFLAG_TVIPS_DATA;
    if (TIFFIsByteSwapped(tif))
      inFile->userFlags |= IIFLAG_BYTES_SWAPPED;
  }

  /* Use min and max from file if defined (better be there for float/int) */
  if (gotMin) {
    if (inFile->type == IITYPE_BYTE)
      lastMin += 128.;
    inFile->amin = lastMin;
  }
  if (gotMax) {
    if (inFile->type == IITYPE_BYTE)
      lastMax += 128.;
    inFile->amax = lastMax;
  }
  inFile->amean = (inFile->amin + inFile->amax) / 2.;

  /* Handle pixel size if it is above threshold */
  if (hasPixel && (inFile->anyTiffPixSize || xPixel / 1.e4 <= pixelLimit)) {
    inFile->xscale = xPixel;
    inFile->yscale = yPixel;
    inFile->zscale = xPixel;
  }

  inFile->smin   = inFile->amin;
  inFile->smax   = inFile->amax;

  /* Intercept ImageJ big tiff and try to convert to mrc-like */
  if (imageJ && dirnum == 1 && imJslices == imJimages && imJimages > 1 && 
      inFile->planesPerImage == 1 && photometric != PHOTOMETRIC_PALETTE && 
      inFile->mode >= 0) {
    info.nx = inFile->nx;
    info.ny = inFile->ny;
    info.nz = imJimages;
    info.swapBytes = TIFFIsByteSwapped(tif);
    info.sectionSkip = 0;
    info.yInverted = 1;
    info.amin = inFile->amin;
    info.amax = inFile->amax;
    info.pixel = hasPixel ? xPixel : 0.;
    info.zPixel = info.pixel;
    if (TIFFGetField(tif, TIFFTAG_STRIPOFFSETS, &offsets) > 0) {
      info.headerSize = *offsets;
      TIFFClose(tif);
      inFile->fp = fopen(inFile->filename, inFile->fmode);
      if (!inFile->fp) {
        b3dError(stderr, "ERROR: iiTIFFCheck - Reopening file %s for treatment as "
                 "MRC-like\n", inFile->filename);
        return(IIERR_IO_ERROR);
      }
      return iiSetupRawHeaders(inFile, &info);
    }
  }

  /* Otherwise proceed to return as a TIFF file */
  inFile->headerSize = 8;
  inFile->sectionSkip = 0;
  inFile->fp = (FILE *)tif;    
  inFile->cleanUp = tiffDelete;
  inFile->reopen = tiffReopen;
  inFile->close = tiffClose;
  inFile->fillMrcHeader = tiffFillMrcHeader;
  inFile->lastWrittenZ = inFile->nz - 1;
  return(0);
}

int tiffReopen(ImodImageFile *inFile)
{
  TIFF* tif;
  tif = openWithoutBMode(inFile);
  if (!tif)
    return 1;
  inFile->headerSize = 8;
  inFile->sectionSkip = 0;
  inFile->header = (char *)tif;    
  inFile->fp = (FILE *)tif;    
  return 0;
}

void tiffClose(ImodImageFile *inFile)
{
  TIFF* tif = (TIFF *)inFile->header;

  /* After writing a new file, if it has more than one section and there is a min and
     max, write them to last directory */
  if (inFile->newFile && inFile->format != IIFORMAT_RGB
      && tif && inFile->amax > inFile->amin) {
    constrainAndStoreMinMax(inFile);
  }
  if (tif)
    TIFFClose(tif);
  inFile->header = NULL;
  inFile->fp = NULL;
}

void tiffDelete(ImodImageFile *inFile)
{
  tiffClose(inFile);
}

int tiffFillMrcHeader(ImodImageFile *inFile, MrcHeader *hdata)
{
  char *description, *endPtr;
  TIFF* tif = (TIFF *)inFile->header;
  int startInd, endInd, descLen;
  mrc_head_new(hdata, inFile->nx, inFile->ny, inFile->nz, inFile->mode);
  hdata->bytesSigned = inFile->type == IITYPE_BYTE ? 1 : 0;
  hdata->amin = inFile->amin;
  hdata->amean = inFile->amean;
  hdata->amax = inFile->amax;
  mrc_set_scale(hdata, (double)inFile->xscale, (double)inFile->yscale,
                (double)inFile->zscale);
  hdata->fp = inFile->fp;
  hdata->packed4bits = inFile->packed4bits;
  if (hdata->packed4bits == PACKED_4BIT_MODE)
    hdata->iiuFlags |= IIUNIT_4BIT_MODE;
  if (hdata->packed4bits == PACKED_HALF_XSIZE)
    hdata->iiuFlags |= IIUNIT_HALF_XSIZE;

  /* Get description and turn into titles */
  TIFFSetDirectory(tif, 0);
  if (TIFFGetField(tif, TIFFTAG_IMAGEDESCRIPTION, &description)) {
    hdata->nlabl = 0;
    startInd = 0;
    descLen = strlen(description);
    while (hdata->nlabl < MRC_NLABELS && startInd < descLen) {
      endPtr = strchr(&description[startInd], '\n');
      if (!endPtr) {
        strncpy(hdata->labels[hdata->nlabl], &description[startInd], MRC_LABEL_SIZE);
        fixTitlePadding(hdata->labels[hdata->nlabl++]);
        break;
      } else {
        endInd = endPtr - description;
        if (endInd && description[endInd - 1] == '\r')
          endInd--;
        if (endInd - startInd > 0) {
          strncpy(hdata->labels[hdata->nlabl], &description[startInd], 
                  B3DMIN(endInd - startInd, MRC_LABEL_SIZE));
          if (endInd - startInd < MRC_LABEL_SIZE)
            hdata->labels[hdata->nlabl][endInd - startInd] = 0x00;
          fixTitlePadding(hdata->labels[hdata->nlabl++]);
        }
        if (description[endInd] == '\r')
          endInd++;
        startInd = endInd + 1;
      }
    }
  }
  
  return 0;
}

/*
 * Sync header: convert (unterminated) title strings into a single description string
 */
static int tiffSyncFromMrcHeader(ImodImageFile *inFile, MrcHeader *hdata)
{
  int ind, lab, trueLen, startInd, outInd = 0;
  char *bitStr;
  B3DFREE(inFile->description);
  if (!hdata->nlabl)
    return 0;

  /* Allocate array big enough for full lines */
  inFile->description = (char *)malloc(hdata->nlabl * (MRC_LABEL_SIZE + 1));
  if (!inFile->description)
    return 1;

  /* Find true length as last non-blank */
  for (lab = 0 ; lab < hdata->nlabl; lab++) {
    trueLen = 0;
    for (ind = 0; ind < MRC_LABEL_SIZE; ind++) {
      if (hdata->labels[lab][ind] == '\n')
        break;
      if (hdata->labels[lab][ind] != ' ')
        trueLen = ind + 1;
    }

    /* Copy string and terminate it with return or NULL at end */
    if (trueLen) {
      startInd = outInd;
      for (ind = 0; ind < trueLen; ind++) {
        inFile->description[outInd++] = hdata->labels[lab][ind];
      }

      /* Convert a string from SEMCCD so it won't be recognized by earlier versions */
      bitStr = strstr(&inFile->description[startInd], "4 bits packed");
      if (bitStr && strstr(&inFile->description[startInd], "SerialEMCCD"))
        bitStr[1] = '-';
      inFile->description[outInd++] = (lab == hdata->nlabl - 1) ? 0x00 : '\n';
    } else if (lab == hdata->nlabl - 1)
      inFile->description[outInd++] = 0x00;
  }
  tiffAddDescription(inFile->description);
  return 0;
}

/* Get the value for a field that returns a single value */
int tiffGetField(ImodImageFile *inFile, int tag, void *value)
{
  TIFF *tif;
  if (!inFile)
    return -1;
  if (!inFile->header)
    iiReopen(inFile);
  tif = (TIFF *)inFile->header;
  if (!tif)
    return -1;
  return TIFFGetField(tif, (ttag_t)tag, value);
}
/* Get the value for a field that returns the address of an array.  The count
   argument seems to be required but does not seem to return the count */
int tiffGetArray(ImodImageFile *inFile, int tag, uint16 *count, void *value)
{
  TIFF *tif;
  if (!inFile)
    return -1;
  if (!inFile->header)
    iiReopen(inFile);
  tif = (TIFF *)inFile->header;
  if (!tif)
    return -1;
  return TIFFGetField(tif, (ttag_t)tag, count, value);
}

void tiffSuppressErrors(void)
{
  TIFFSetErrorHandler(NULL);
}


void tiffSuppressWarnings(void)
{
  TIFFSetWarningHandler(NULL);
  sWarningsSuppressed = 1;
}

static void warningHandler(const char *module, const char *fmt, va_list ap)
{
  char buffer[160];
  vsnprintf(buffer, 159, fmt, ap);
  va_end(ap);

  /* It didn't work to call the old handler with some errors, so print it
     ourselves to stderr.  It was "unknown" in libtiff 3 and "Unknown" in 4 */
  if (!strstr(buffer, "nknown field with tag")) {
    if (module)
      fprintf(stderr, "%s: Warning, %s\n", module, buffer);
    else
      fprintf(stderr, "Warning, %s\n", buffer);
  }
}

/* Set up to filter out the ubiquitous unknown field warnings */
void tiffFilterWarnings(void)
{
  if (!sWarningsSuppressed)
    oldHandler = TIFFSetWarningHandler(warningHandler);
}

void tiffSetMapping(int value)
{
  sUseMapping = value;
}

/* Mode 'b' means something completely different for TIFF, so strip it */
static TIFF *openWithoutBMode(ImodImageFile *inFile)
{
  TIFF *tif;
  int len, stripped = 0;
  char *tmpmode = inFile->fmode;
  if (!tmpmode)
    return NULL;
  len = strlen(tmpmode);
  if (!len)
    return NULL;

  /* If mapping is still 2, check the variable and turn it off or set to 1 */
  if (sUseMapping == 2) {
    if (getenv("IMOD_NO_TIFF_MEM_MAP"))
      sUseMapping = 0;
    else
      sUseMapping = 1;
  }

  if (!sUseMapping || inFile->fmode[len - 1] == 'b') {
    stripped = 1;
    tmpmode = (char *)malloc(len + 4);
    if (!tmpmode)
      return NULL;
    strcpy(tmpmode, inFile->fmode);
    if (inFile->fmode[len - 1] == 'b')
      tmpmode[len - 1] = 0x00;
    if (!sUseMapping) {
      len = strlen(tmpmode);
      tmpmode[len] = 'm';
      tmpmode[len+1] = 0x00;
    }
  }
  tif = TIFFOpen(inFile->filename, tmpmode);
  if (stripped)
    free(tmpmode);
  return tif;
}

/* Find the directory at the given number that matches the proper size */
static int setMatchingDirectory(ImodImageFile *inFile, int dirnum)
{
  int nx, ny, dir = 0;
  TIFF *tif = (TIFF *)inFile->header;
  TIFFSetDirectory(tif, 0);
  do {
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &nx);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &ny);
    if (inFile->packed4bits == PACKED_HALF_XSIZE)
      nx *= 2;
    if (nx == inFile->nx && ny == inFile->ny) {
      if (dir == dirnum)
        return 0;
      else
        dir++;
    }
  } while (TIFFReadDirectory(tif));
  return 1;
}

static void closeWithError(ImodImageFile *inFile, char *message)
{
  TIFFClose((TIFF *)inFile->header);
  inFile->fp = fopen(inFile->filename, inFile->fmode);
  b3dError(stderr, message);
}

/* DNM 12/24/00: Got this working for bytes, shorts, and RGBs, for whole
   images or subsets, and using maps for scaling */
/* DNM 11/18/01: Added ability to read tiles, made tiffReadSection and 
   tiffReadSectionByte call a common routine to reduce duplicate code */

static int ReadSection(ImodImageFile *inFile, char *buf, int inSection, int convert)
{
  int nstrip, si, plane, sampOffset;
  int xout, xcopy;
  int xsize = inFile->nx;
  int ysize = inFile->ny;
  int samples = inFile->contigSamples;
  int row;
  int xstart, xend, ystart, yend, y, ofsin;
  size_t ofsout;
  int doscale;
  int byte = convert == MRSA_BYTE ? 1 : 0;
  int toShort = convert == MRSA_USHORT ? 1 : 0;
  int toFloat = convert == MRSA_FLOAT ? 1 : 0;
  int signedBytes = inFile->type == IITYPE_BYTE;
  float slope = inFile->slope;
  float offset = inFile->offset;
  int outmin = 0;
  int outmax = toShort ? 65535 : 255;
  float eps = toShort ? 0.005 / 256. : 0.005;
  int padLeft = B3DMAX(0, inFile->padLeft);
  int padRight = B3DMAX(0, inFile->padRight);
  int stripsize;
  int xmin, xmax, ymin, ymax;
  int pixsize = 1;
  int moveSize = toShort ? 2 : 1;
  unsigned char *obuf;
  unsigned char *tmp = NULL;
  unsigned char *bdata;
  unsigned char *map = NULL;
  unsigned char first4bitsMap[256];
  unsigned char second4bitsMap[256];
  unsigned char *highMap = &first4bitsMap[0];
  unsigned char *lowMap = &second4bitsMap[0];

  /* Send a color map to line converter for a colormap image if converting or
     returning RGB values, but not if the actuals bytes are wanted */
  unsigned char *colormap = B3DCHOICE(inFile->format == IIFORMAT_COLORMAP && 
                                      (convert || !inFile->rawPaletteBytes), 
                                      &inFile->colormap[768 * inSection], NULL);
  int freeMap = 0;
  uint32 rowsperstrip;
  int nread, lineBytes, lineOffset, skipHalf;
  int tilesize, tilewidth, tilelength, xtiles, ytiles, xti, yti, xDimension;

  TIFF* tif = (TIFF *)inFile->header;
  if (inFile->axis == 2) {
    b3dError(stderr, "ERROR: TIFF ReadSection - Cannot read Y planes from a TIFF file\n");
    return -1;
  }

  if (!tif)
    iiReopen(inFile);
  tif = (TIFF *)inFile->header;
  if (!tif)
    return -1;

  /* set the dimensions to read in */
  /* DNM 2/26/03: replace upper right only if negative */
  xmin   = inFile->llx;
  ymin   = inFile->lly;
  if (inFile->urx < 0)
    xmax = inFile->nx-1;
  else
    xmax = inFile->urx;
  if (inFile->ury < 0)
    ymax = inFile->ny-1;
  else
    ymax = inFile->ury;
  xout = xmax + 1 - xmin;
  xDimension = xout + padLeft + padRight;
  doscale = (convert && (offset <= -1.0 || offset >= 1.0 || 
                         slope < 1. - eps || slope > 1. + eps)) ? 1 : 0;
  lineBytes = xsize;
  lineOffset = xmin;
  skipHalf = 0;

  /* Modify bytes and offset appropriately and set flag to skip half byte at start of line
     for 4-bit, also set up byte to first/second pixel maps appropriate for fill order */
  if (inFile->packed4bits) {
    lineBytes = (xsize + 1) / 2;
    lineOffset = xmin / 2;
    skipHalf = xmin % 2;
    if (inFile->fillOrder == FILLORDER_LSB2MSB) {
      lowMap = &first4bitsMap[0];
      highMap = &second4bitsMap[0];
    }
    for (xti = 0; xti < 16; xti++) {
      for (yti = 0; yti < 16; yti++) {
        lowMap[xti + 16 * yti] = xti;
        highMap[xti + 16 * yti] = yti;
      }
    }
  }

  row = inSection / (inFile->planesPerImage * samples);
  if (setMatchingDirectory(inFile, row)) {
    b3dError(stderr, "ERROR: TIFF ReadSection - Cannot find directory %d\n", row);
    return -1;
  }
  plane = inSection % inFile->planesPerImage;
  sampOffset = inSection % samples;

  /* Set up pixsize which is the number of bytes in the input data, and 
     moveSize which is number of bytes of output.  Also get scale maps */
  if ((convert && !toFloat) || signedBytes) {
    if (inFile->type == IITYPE_SHORT) {
      pixsize = 2;
      map = get_short_map(slope, offset, outmin, outmax, MRC_RAMP_LIN, 0, 1);
      freeMap = 1;
    } else if (inFile->type == IITYPE_USHORT) {
      pixsize = 2;
      if (byte || doscale) {
        map = get_short_map(slope, offset, outmin, outmax, MRC_RAMP_LIN, 0, 0);
        freeMap = 1;
      }
    } else if (inFile->type == IITYPE_FLOAT || inFile->type == IITYPE_INT ||
               inFile->type == IITYPE_UINT) {
      pixsize = 4;
    } else if (((toShort || doscale) && !colormap) || signedBytes) {
      /* printf("slope %f  offset %f  outmin %d outmax %d\n", slope, offset, outmin, 
         outmax); */
      map = get_byte_map(slope, offset, outmin, outmax, signedBytes ? 1 : 0);
      if (inFile->format == IIFORMAT_RGB)
        pixsize = inFile->rgbSamples;
    } else if (inFile->format == IIFORMAT_RGB) {
      pixsize = inFile->rgbSamples;
    }
    if (toFloat)
      moveSize = 4;
  } else {
    if (inFile->format == IIFORMAT_RGB)
      pixsize = inFile->rgbSamples;
    else if (inFile->type == IITYPE_SHORT || inFile->type == IITYPE_USHORT)
      pixsize = 2;
    else if (inFile->type == IITYPE_FLOAT || inFile->type == IITYPE_INT ||
             inFile->type == IITYPE_UINT)
      pixsize = 4;
    if (toFloat)
      moveSize = 4;
    else
      moveSize = (inFile->format == IIFORMAT_RGB || colormap) ? 3 : pixsize;
  }

  if (freeMap && !map)
    return -1;
  /*printf("byte %d scale %d samples %d pixsize %d convert %d map %p colormap %p\n",
    byte, doscale, samples, pixsize, convert, map, colormap); */

  if (TIFFGetField(tif, TIFFTAG_ROWSPERSTRIP, &rowsperstrip)) {

    /* if data are in strips, get strip size and memory for it */
    stripsize = TIFFStripSize(tif);
    tmp = (unsigned char *)_TIFFmalloc(stripsize);
    if (!tmp) {
      if (freeMap)
        free(map);
      return -1;
    }
               
    nstrip = TIFFNumberOfStrips(tif) / inFile->planesPerImage;
    /*printf("%d %d %d %d\n", stripsize, (int)rowsperstrip, nstrip, pixsize);*/

    for (si = 0 ; si < nstrip; si++){
               
      /* Compute starting and ending Y values to use in each strip */
      ystart = ysize - 1 - (rowsperstrip * (si + 1) - 1);
      yend = ysize - 1 - (rowsperstrip * si);
      if (ymin > ystart)
        ystart = ymin;
      if (ymax < yend)
        yend = ymax;
      if (ystart > yend)
        continue;
               
      /* Read the strip if necessary */
      nread = TIFFReadEncodedStrip(tif, si + plane * nstrip, tmp, stripsize);
      /*printf("\n%d %d %d %d\n", nread, si, ystart, yend);*/
      for (y = ystart; y <= yend; y++) {

        /* for each y, compute back to row, and get offsets into
           input and output arrays */
        row = ysize - 1 - y - rowsperstrip * si;
        ofsin = samples * pixsize * (row * lineBytes + lineOffset) + sampOffset *pixsize;
        ofsout = moveSize * ((size_t)(y - ymin) * (size_t)xDimension + padLeft);
        obuf = (unsigned char *)buf + ofsout;
        bdata = tmp + ofsin;
        copyLine(bdata, obuf, xout, convert, pixsize, inFile->type, 
                 inFile->format, samples, slope, offset, doscale,
                 inFile->packed4bits ? 1 + skipHalf : 0, map, first4bitsMap, 
                 second4bitsMap, colormap);
      }    
    }
  } else {

    /* Otherwise make sure there are tiles, if not return with error */
    if (TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tilewidth)) {
      tilesize = TIFFTileSize(tif);
      tmp = (unsigned char *)_TIFFmalloc(tilesize);
    }
    if (!tmp) {
      if (freeMap)
        free(map);
      return -1;
    }
    TIFFGetField(tif, TIFFTAG_TILELENGTH, &tilelength);
    xtiles = (xsize + tilewidth - 1) / tilewidth;
    ytiles = (ysize + tilelength - 1) / tilelength;
               
    /*printf("%d %d %d %d %d %d\n", tilesize, tilewidth, tilelength,
      xtiles, ytiles, pixsize); */

    for (yti = 0; yti < ytiles; yti++) {
      for (xti = 0; xti < xtiles; xti++) {
                    
        /* Compute starting and ending Y then X values to use in 
           this tile */
        ystart = ysize - 1 - (tilelength * (yti + 1) - 1);
        yend = ysize - 1 - (tilelength * yti);
        if (ymin > ystart)
          ystart = ymin;
        if (ymax < yend)
          yend = ymax;
        if (ystart > yend)
          continue;

        xstart = xti * tilewidth;
        xend = xstart + tilewidth - 1;
        if (xmin > xstart)
          xstart = xmin;
        if (xmax < xend)
          xend = xmax;
        if (xstart > xend)
          continue;
                    
        /* Read the tile if necessary */
        si = xti + yti * xtiles + plane * xtiles * ytiles;
        nread = TIFFReadEncodedTile(tif, si, tmp, tilesize);
        xcopy = xend + 1 - xstart;

        /* Set up bytes and offset appropriately for this tile */
        lineBytes = tilewidth;
        lineOffset = xstart - xti * tilewidth;
        skipHalf = 0;
        if (inFile->packed4bits) {
          lineBytes = (tilewidth + 1) / 2;
          skipHalf = lineOffset % 2;
          lineOffset /= 2;
        }

        /* printf("%d %d %d %d\n", nread, si, ystart, yend); */
        for (y = ystart; y <= yend; y++) {
                         
          /* for each y, compute back to row, and get offsets 
             into input and output arrays */
          row = ysize - 1 - y - tilelength * yti;
          ofsin = pixsize * samples * (row * lineBytes + lineOffset) + 
            sampOffset * pixsize;
          ofsout = moveSize * 
            ((size_t)(y - ymin) * (size_t)xDimension + padLeft + xstart - xmin);
          obuf = (unsigned char *)buf + ofsout;
          bdata = tmp + ofsin;
          copyLine(bdata, obuf, xcopy, convert, pixsize, inFile->type,
                   inFile->format, samples, slope, offset, doscale, 
                   inFile->packed4bits ? 1 + skipHalf : 0, map, first4bitsMap, 
                   second4bitsMap, colormap);
        }
      }               
    }
  }
  _TIFFfree(tmp);
  if (freeMap)
    free(map);

  return 0;
}

#define RGB_TO_FLOAT()                          \
  fpixel = 0.3 * bdata[0];                      \
  fpixel += 0.59 * bdata[1];                    \
  fpixel += 0.11 * bdata[2];                    \
  bdata += pixsize;


/*
 * Copy one line of data appropriately for the data type and conversion
 */
static void copyLine(unsigned char *bdata, unsigned char *obuf, int xout, int convert,
                     int pixsize, int type, int format, int samples, float slope,
                     float offset, int doscale, int unpack4bits, unsigned char *map,
                     unsigned char *first4bitsMap, unsigned char *second4bitsMap,
                     unsigned char *colormap)
{
  b3dUInt16 *usdata;
  b3dInt16 *sdata;
  b3dUInt32 *uldata;
  b3dInt32 *ldata;
  b3dFloat *fdata;
  b3dUInt16 *usmap = (b3dUInt16 *)map;
  b3dUInt16 *usobuf = (b3dUInt16 *)obuf;
  b3dFloat *fobuf = (b3dFloat *)obuf;
  int toShort = convert == MRSA_USHORT ? 1 : 0;
  int toFloat = convert == MRSA_FLOAT ? 1 : 0;
  int outmax = toShort ? 65535 : 255;
  int signedBytes = type == IITYPE_BYTE;
  int i, j, ival;
  float fpixel;

  if (convert || signedBytes || unpack4bits) {

    /* Converted data */
    if (pixsize == 1) {
      
      /* Converted from bytes: 4 bits */
      if (unpack4bits) {
        i = 0;

        /* If there is an odd pixel at start, get it */
        if (unpack4bits > 1) {
          if (toFloat) 
            *fobuf++ = second4bitsMap[bdata[i++]];
          else if (toShort)
            *usobuf++ = usmap[second4bitsMap[bdata[i++]]];
          else if (doscale)
            *obuf++ = map[second4bitsMap[bdata[i++]]];
          else
            *obuf++ = second4bitsMap[bdata[i++]];
          xout--;
        }

        /* Process pairs of pixels */
        if (toFloat) {
          for (; i < xout / 2; i++) {
            *fobuf++ = first4bitsMap[bdata[i]];
            *fobuf++ = second4bitsMap[bdata[i]];
          }
        } else if (toShort) {
          for (; i < xout / 2; i++) {
            *usobuf++ = usmap[first4bitsMap[bdata[i]]];
            *usobuf++ = usmap[second4bitsMap[bdata[i]]];
          }
        } else if (doscale) {
          for (; i < xout / 2; i++) {
            *obuf++ = map[first4bitsMap[bdata[i]]];
            *obuf++ = map[second4bitsMap[bdata[i]]];
          }
        } else {
          for (; i < xout / 2; i++) {
            *obuf++ = first4bitsMap[bdata[i]];
            *obuf++ = second4bitsMap[bdata[i]];
          }
        }

        /* Then do odd pixel at end */
        if (xout % 2) {
          if (toFloat) 
            *fobuf++ = first4bitsMap[bdata[i]];
          else if (toShort)
            *usobuf++ = usmap[first4bitsMap[bdata[i]]];
          else if (doscale)
            *obuf++ = map[first4bitsMap[bdata[i]]];
          else
            *obuf++ = first4bitsMap[bdata[i]];
        }

      } else if (format == IIFORMAT_COLORMAP) {

        /* RGB in a colormap with conversion to floats, shorts or bytes */
        if (toFloat) {
          for (i = 0; i < xout; i++) {
            fpixel = 0.3 * colormap[*bdata];
            fpixel += 0.59 * colormap[256 + *bdata];
            fpixel += 0.11 * colormap[512 + *bdata++];
            *fobuf++ = fpixel;
          }
        } else if (toShort) {
          for (i = 0; i < xout; i++) {
            fpixel = 255. * 0.3 * colormap[*bdata];
            fpixel += 255. * 0.59 * colormap[256 + *bdata];
            fpixel += 255. * 0.11 * colormap[512 + *bdata++];
            *usobuf++ = (int)(fpixel + 0.5f);
          }
        } else {
          for (i = 0; i < xout; i++) {
            fpixel = 0.3 * colormap[*bdata];
            fpixel += 0.59 * colormap[256 + *bdata];
            fpixel += 0.11 * colormap[512 + *bdata++];
            *obuf++ = (int)(fpixel + 0.5f);
          }
        }

      } else if (samples == 1) {

        /* Single-sample conversions of bytes to float, short, mapped or unmapped bytes */
        if (toFloat && signedBytes)
          for (i = 0; i < xout; i++)
            *fobuf++ = map[*bdata++];
        else if (toFloat)
          for (i = 0; i < xout; i++)
            *fobuf++ = *bdata++;
        else if (toShort) 
          for (i = 0; i < xout; i++)
            *usobuf++ = usmap[*bdata++];
        else if (doscale || signedBytes)
          for (i = 0; i < xout; i++)
            *obuf++ = map[*bdata++];
        else
          memcpy(obuf, bdata, xout);
      } else {

        /* Interleaved sample conversions of bytes to float, short, mapped or unmapped 
           bytes bytes */
        if (toFloat && signedBytes)
          for (i = 0; i < xout; i++) {
            *fobuf++ = map[*bdata];
            bdata += samples;
          }
        else if (toFloat)
          for (i = 0; i < xout; i++) {
            *fobuf++ = *bdata;
            bdata += samples;
          }
        else if (toShort)
          for (i = 0; i < xout; i++) {
            *usobuf++ = usmap[*bdata];
            bdata += samples;
          }
        else if (doscale || signedBytes) 
          for (i = 0; i < xout; i++) {
            *obuf++ = map[*bdata];
            bdata += samples;
          }
        else
          for (i = 0; i < xout; i++) {
            *obuf++ = *bdata;
            bdata += samples;
          }
      }
      
    } else if (pixsize == 2) {
      
      /* Integers converted */
      usdata = (b3dUInt16 *)bdata;
      sdata = (b3dInt16 *)bdata;
      if (samples == 1) {

        /* single-sample integer conversions to mapped short, float, unmapped short, or
           mapped byte */
        if (toShort && map) 
          for (i = 0; i < xout; i++)
            *usobuf++ = usmap[*usdata++];
        else if (toFloat && type == IITYPE_SHORT) 
          for (i = 0; i < xout; i++)
            *fobuf++ = *sdata++;
        else if (toFloat) 
          for (i = 0; i < xout; i++)
            *fobuf++ = *usdata++;
        else if (toShort) 
          for (i = 0; i < xout; i++)
            *usobuf++ = *usdata++;
        else 
          for (i = 0; i < xout; i++)
            *obuf++ = map[*usdata++];
      } else {

        /* interleaved sample conversions to mapped short, float, unmapped short, or
           mapped byte */
        if (toShort && map) 
          for (i = 0; i < xout; i++) {
            *usobuf++ = usmap[*usdata];
            usdata += samples;
          }
        else if (toFloat && type == IITYPE_SHORT)
          for (i = 0; i < xout; i++) {
            *fobuf++ = *sdata;
            sdata += samples;
          }
        else if (toFloat)
          for (i = 0; i < xout; i++) {
            *fobuf++ = *usdata;
            usdata += samples;
          }
        else if (toShort)
          for (i = 0; i < xout; i++) {
            *usobuf++ = *usdata;
            usdata += samples;
          }
        else
          for (i = 0; i < xout; i++) {
            *obuf++ = map[*usdata];
            usdata += samples;
          }
      }
      
    } else if (type == IITYPE_INT) {

      /* Long ints */
      ldata = (b3dInt32 *)bdata;
      if (toFloat) {
        for (i = 0; i < xout; i++) {
          *fobuf++ = *ldata;
          ldata += samples;
        }
      } else {
        for (i = 0; i < xout; i++) {
          ival = slope * (*ldata) + offset;
          if (toShort)
            *usobuf++ = B3DMAX(0, B3DMIN(outmax, ival));
          else
            *obuf++ = B3DMAX(0, B3DMIN(outmax, ival));
          ldata += samples;
        }
      }

    } else if (type == IITYPE_UINT) {

      /* Unsigned Long ints */
      uldata = (b3dUInt32 *)bdata;
      if (toFloat) {
        for (i = 0; i < xout; i++) {
          *fobuf++ = *uldata;
          uldata += samples;
        }
      } else {
        for (i = 0; i < xout; i++) {
          ival = slope * (*uldata) + offset;
          if (toShort)
            *usobuf++ = B3DMAX(0, B3DMIN(outmax, ival));
          else
            *obuf++ = B3DMAX(0, B3DMIN(outmax, ival));
          uldata += samples;
        }
      }

    } else if (format == IIFORMAT_RGB) {

      /* RGB conversions to float, mapped short, mapped byte, unmapped short or byte */
      if (toFloat) {
        for (i = 0; i < xout; i++) {
          RGB_TO_FLOAT();
          *fobuf++ = fpixel;
        } 
      } else if (doscale && toShort) {
        for (i = 0; i < xout; i++) {
          RGB_TO_FLOAT();
          *usobuf++ = usmap[(int)(fpixel + 0.499f)];
        }
      } else if (doscale) {
        for (i = 0; i < xout; i++) {
          RGB_TO_FLOAT();
          *obuf++ = map[(int)(fpixel + 0.499f)];
        }
      } else if (toShort) {
        for (i = 0; i < xout; i++) {
          fpixel = 255. * 0.3 * bdata[0];
          fpixel += 255. * 0.59 * bdata[1];
          fpixel += 255. * 0.11 * bdata[2];
          bdata += pixsize;
          *usobuf++ = (int)(fpixel + 0.5f);
        }
      } else {
        for (i = 0; i < xout; i++) {
          RGB_TO_FLOAT();
          *obuf++ = (int)(fpixel + 0.5f);
        }
      }

    } else {
      
      /* Floats with interleaved samples or converted to short or byte */
      fdata = (b3dFloat *)bdata;
      if (toFloat) {
        for (i = 0; i < xout; i++) {
          *fobuf++ = *fdata;
          fdata += samples;
        }
      } else {
        if (toShort) {
          for (i = 0; i < xout; i++) {
            ival = (int)(slope * (*fdata) + offset);
            *usobuf++ = B3DMAX(0, B3DMIN(outmax, ival));
            fdata += samples;
          }
        } else {
          for (i = 0; i < xout; i++) {
            ival = (int)(slope * (*fdata) + offset);
            *obuf++ = B3DMAX(0, B3DMIN(outmax, ival));
            fdata += samples;
          }
        }
      }
    }
  } else {
    
    /* Non-converted data */
    if (samples == 1) {
      /* RGB with extra samples - skip them */
      if (format == IIFORMAT_RGB && pixsize > 3) {
        for (i = 0; i < xout; i++) {
          for (j = 0; j < 3; j++)
            *obuf++ = *bdata++;
          bdata += pixsize - 3;
        }

        /* Colormap lookup */
      } else if (colormap) {
        for (i = 0; i < xout; i++) {
          *obuf++ = colormap[*bdata];
          *obuf++ = colormap[256 + *bdata];
          *obuf++ = colormap[512 + *bdata++];
        }
          
        /* Straight copy */
      } else
        memcpy(obuf, bdata, xout * pixsize);
    } else {

      /* Multiple samples (planes) interleaved - skip the other samples */
      for (i = 0; i < xout; i++) {
        for (j = 0; j < pixsize; j++)
          *obuf++ = *bdata++;
        bdata += (samples - 1) * pixsize;
      }
    }
  }
}    

int tiffReadSectionByte(ImodImageFile *inFile, char *buf, int inSection)
{ 
  return(ReadSection(inFile, buf, inSection, MRSA_BYTE));
}

int tiffReadSectionUShort(ImodImageFile *inFile, char *buf, int inSection)
{ 
  return(ReadSection(inFile, buf, inSection, MRSA_USHORT));
}

int tiffReadSectionFloat(ImodImageFile *inFile, char *buf, int inSection)
{ 
  return(ReadSection(inFile, buf, inSection, MRSA_FLOAT));
}

int tiffReadSection(ImodImageFile *inFile, char *buf, int inSection)
{
  return(ReadSection(inFile, buf, inSection, MRSA_NOPROC));
}

/* Definitions just to make this compile for TIFF 3; they don't have to be correct */
#if TIFFLIB_VERSION < 20110401
typedef unsigned long uint64;
typedef size_t tmsize_t;
#endif

/* Parallel writing statics */
static int sFileBufSize, sSettingUpParallel = 0;
static char *sFileBuf[MAX_TIFF_THREADS];
static uint64 sCurBufInd[MAX_TIFF_THREADS];
static uint64 sMaxBufInd[MAX_TIFF_THREADS];

/*
 * Parallel writing procedures to operate on buffers instead of files
 */
static tmsize_t bufReadProc(thandle_t fd, void* buf, tmsize_t size)
{
  int fileNo = (int)fd;
  if (sCurBufInd[fileNo] + size > sFileBufSize) {
    errno = EINVAL;
    return (tmsize_t)-1;
  }
  /*printf("read %d %d %d\n", fileNo, (int)sCurBufInd[fileNo], (int)size);*/
  memcpy(buf, sFileBuf[fileNo] + sCurBufInd[fileNo], size);
  sCurBufInd[fileNo] += size;
  ACCUM_MAX(sMaxBufInd[fileNo], sCurBufInd[fileNo]);
  return size;
}

static tmsize_t bufWriteProc(thandle_t fd, void* buf, tmsize_t size)
{
  int fileNo = (int)fd;
  if (sCurBufInd[fileNo] + size >= sFileBufSize) {
    errno = EINVAL;
    return (tmsize_t)-1;
  }
  /*printf("write %d %d %d\n", fileNo, (int)sCurBufInd[fileNo], (int)size);*/
  memcpy(sFileBuf[fileNo] +  sCurBufInd[fileNo], buf, size);
  sCurBufInd[fileNo] += size;
  ACCUM_MAX(sMaxBufInd[fileNo], sCurBufInd[fileNo]);
  return size;
}

static uint64 bufSeekProc(thandle_t fd, uint64 off, int whence)
{
  int fileNo = (int)fd;
  uint64 newPos;
  off_t off_io = (off_t) off;
  int err = 1;

  switch (whence) {
  case SEEK_SET:
    if (off_io < 0)
      break;
    newPos = off_io;
    err = 0;
    break;
  case SEEK_CUR:
    if (off_io < 0 && sCurBufInd[fileNo] < -off_io)
      break;
    newPos = sCurBufInd[fileNo] + off_io;
    err = 0;
    break;
  case SEEK_END:
    newPos = sMaxBufInd[fileNo];
    err = 0;
    break;
  default:
    break;
  }
  if (err || newPos >= sFileBufSize) {
    errno = EINVAL;
    return (uint64)-1;
  }
  sCurBufInd[fileNo] = newPos;
  ACCUM_MAX(sMaxBufInd[fileNo], sCurBufInd[fileNo]);
  return newPos;
}

static int bufCloseProc(thandle_t fd)
{
  return 0;
}

static uint64 bufSizeProc(thandle_t fd)
{
  return sMaxBufInd[(int)fd];
}

/*
 * Open new file for writing
 */
int tiffOpenNew(ImodImageFile *inFile)
{
  TIFF *tif;
  int minor, psize = 1, usew8 = 0;

  if (inFile->type == IITYPE_SHORT || inFile->type == IITYPE_USHORT)
    psize = 2;
  if (inFile->type == IITYPE_INT || inFile->type == IITYPE_UINT || 
      inFile->type == IITYPE_FLOAT)
    psize = 4;
  if (inFile->format == IIFORMAT_RGB)
    psize *= 3;
  
  /* Try to open a big file in version 4 */
  if (tiffVersion(&minor) > 3 && (((double)inFile->nx * inFile->ny) * inFile->nz * psize
                                  > 4.0e9 || makeAllBigTiff()))
    usew8 = 1;
  if (sSettingUpParallel <= 0)
    tif = TIFFOpen(inFile->filename, usew8 ? "w8" : "w");
  else
    tif = TIFFClientOpen(inFile->filename, usew8 ? "w8" : "w", 
                         (thandle_t)(sSettingUpParallel - 1), bufReadProc, bufWriteProc,
                         bufSeekProc, bufCloseProc, bufSizeProc, NULL, NULL);
  if (!tif)
    return IIERR_IO_ERROR;
  inFile->header = (char *)tif;
  inFile->fp = (FILE *)tif;
  inFile->state = IISTATE_READY;
  inFile->cleanUp = tiffDelete; 
  inFile->close = tiffClose;
  inFile->fillMrcHeader = tiffFillMrcHeader;
  inFile->syncFromMrcHeader = tiffSyncFromMrcHeader;
  inFile->writeSection = iiTiffWriteSection;
  inFile->writeSectionFloat = iiTiffWriteSectionFloat;
  return 0;
}

static uint32 sRowsPerStrip, sLineBytes, sStripBytes, sXtileSize;
static int sLinesDone, sNumStrips, sAlreadyInverted, sNumXtiles, sPixSize;
static char *sTmpBuf;
static char *sDescription = NULL;

/*
 * Write next section to file with the given compression value; set inverted
 * non-zero if image is already inverted and does not need copying.  Set 
 * resolution >0 to set a resolution in DPI in the file
 */
int tiffWriteSection(ImodImageFile *inFile, void *buf, int compression, 
                     int inverted, int resolution, int quality)
{
  int strip, lines, err, tileX = 0;
  char *sbuf;
  err = tiffWriteSetup(inFile, compression, inverted, resolution, quality,
                       &strip, &lines, &tileX);
  if (err)
    return err;
  for (strip = 0; strip < sNumStrips; strip++) {
    lines = B3DMIN(sRowsPerStrip, inFile->ny - sLinesDone);
    if (inverted)
      sbuf = (char *)buf + strip * sStripBytes;
    else
      sbuf = (char *)buf + (inFile->ny - (sLinesDone + lines)) * sLineBytes;
    err = tiffWriteStrip(inFile, strip, (void *)sbuf);
    if (err)
      return err;
  }
  tiffWriteFinish(inFile);
  return 0;
}

/*
 * Set up the writing of a section, which then be written as either strips or tiles.
 * Compression is done according to the given compression and quality values.  Data will
 * be inverted before writing unless the inverted entry is nonzero.  outRows is returned
 * with the number of rows per strip or tile; outNum is returned with number of strips or
 * tiles in Y.  tileSizeX must be 0 to write in strips; in this case any incoming value 
 * of outRows is ignored.  To write in tiles, supply the desired size in X in tileSizeX 
 * and the size in Y in outRows.  These sizes will be adjusted to a multiple of 16 that
 * minimizes the padding required to make all tiles equal in size.
 */

int tiffWriteSetup(ImodImageFile *inFile, int compression, int inverted, int resolution,
                   int quality, int *outRows, int *outNum, int *tileSizeX)
{

  /* Use the recommended strip size for uncompressed data, doubling it gives some
     increase in compression at the cost of some time */
  uint32 stripTarget = compression != IICOMPRESSION_NONE ? 16384 : 8192;
  int maxStrips = 4096;
  uint16 bits, samples, photometric, sampleformat;
  time_t curtime;
  struct tm *tm;
  char datetime[40];

  TIFF *tif = (TIFF *)inFile->header;
  if (!(inFile->format == IIFORMAT_RGB || 
        (inFile->format == IIFORMAT_LUMINANCE && 
         (inFile->type == IITYPE_UBYTE || inFile->type == IITYPE_BYTE ||
          inFile->type == IITYPE_USHORT || inFile->type == IITYPE_SHORT ||
          inFile->type == IITYPE_FLOAT))))
    return(IIERR_NO_SUPPORT);
  if (inFile->state == IISTATE_BUSY)
    TIFFWriteDirectory(tif);
  inFile->state = IISTATE_READY;
  inFile->tiffCompression = compression;
  
  TIFFSetField(tif, TIFFTAG_IMAGEWIDTH, inFile->nx);
  TIFFSetField(tif, TIFFTAG_IMAGELENGTH, inFile->ny);
  TIFFSetField(tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(tif, TIFFTAG_COMPRESSION, compression);
  TIFFSetField(tif, TIFFTAG_RESOLUTIONUNIT, 
               resolution > 0 ? RESUNIT_INCH : RESUNIT_NONE);
  if (resolution > 0) {
    TIFFSetField(tif, TIFFTAG_XRESOLUTION, (float)resolution);
    TIFFSetField(tif, TIFFTAG_YRESOLUTION, (float)resolution);
  }
  if (quality >= 0 && compression == COMPRESSION_JPEG) {
    quality = B3DMIN(100, quality);
    TIFFSetField(tif, TIFFTAG_JPEGQUALITY, quality);
  }
  if (quality > 0 && compression == COMPRESSION_ADOBE_DEFLATE) {
    quality = B3DMIN(9, quality);
    TIFFSetField(tif, TIFFTAG_ZIPQUALITY, quality);
  }
  if (inFile->format == IIFORMAT_RGB) {
    samples = 3;
    photometric = PHOTOMETRIC_RGB;
    bits = 8;
    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bits, bits, bits);
  } else {
    
    samples = 1;
    photometric = PHOTOMETRIC_MINISBLACK;
    switch (inFile->type) {
    case IITYPE_BYTE:
      bits = 8;
      sampleformat = SAMPLEFORMAT_INT;
      break;
    case IITYPE_UBYTE:
      bits = 8;
      sampleformat = SAMPLEFORMAT_UINT;
      break;
    case IITYPE_SHORT:
      bits = 16;
      sampleformat = SAMPLEFORMAT_INT;
      break;
    case IITYPE_USHORT:
      bits = 16;
      sampleformat = SAMPLEFORMAT_UINT;
      break;
    case IITYPE_FLOAT:
      bits = 32;
      sampleformat = SAMPLEFORMAT_IEEEFP;
      break;
    }

    TIFFSetField(tif, TIFFTAG_BITSPERSAMPLE, bits);
    TIFFSetField(tif, TIFFTAG_SAMPLEFORMAT, sampleformat);

    /* 11/1/13: save min/max for bytes also, needed for displaying very low count data */
    if (inFile->amax > inFile->amin) {
      constrainAndStoreMinMax(inFile);
    }
  }
  TIFFSetField(tif, TIFFTAG_PHOTOMETRIC, photometric);
  TIFFSetField(tif, TIFFTAG_SAMPLESPERPIXEL, samples);
  sPixSize = samples * bits / 8;

  if (sDescription && sSettingUpParallel <= 0)
    TIFFSetField(tif, TIFFTAG_IMAGEDESCRIPTION, sDescription);
  if (!sSettingUpParallel)
    B3DFREE(sDescription);

  if (*tileSizeX) {

    /* For tiles, make sure each dimension is multiple of 16, it is what tiffcp does */
    sRowsPerStrip = *outRows;
    iiBestTileSize(inFile->nx, tileSizeX, &sNumXtiles, 16);
    iiBestTileSize(inFile->ny, (int *)(&sRowsPerStrip), &sNumStrips, 16);
    sXtileSize = *tileSizeX;
    TIFFSetField(tif, TIFFTAG_TILEWIDTH, sXtileSize);
    TIFFSetField(tif, TIFFTAG_TILELENGTH, sRowsPerStrip);
    sLineBytes = sPixSize * sXtileSize;
    
  } else {
  
    /* Increase the strip target size to avoid having too many strips, as per
       recommendations that go with libtiff4 */
    sLineBytes = sPixSize * inFile->nx;
    if (inFile->ny > maxStrips)
      stripTarget = (1 + inFile->ny / maxStrips) * sLineBytes;
    if (sSettingUpParallel <= 0)
      sRowsPerStrip = B3DMAX(1, (stripTarget + sLineBytes / 2) / sLineBytes);

    /* For JPEG compression, rows must be multiple of 8 */
    if (compression == COMPRESSION_JPEG && sRowsPerStrip % 8) {
      if (sRowsPerStrip < 5 || sRowsPerStrip % 8 > 4)
        sRowsPerStrip = 8 * ((sRowsPerStrip + 7) / 8);
      else
        sRowsPerStrip = 8 * (sRowsPerStrip / 8);
    }
    TIFFSetField(tif, TIFFTAG_ROWSPERSTRIP, sRowsPerStrip);
  }
  sNumStrips = (inFile->ny + sRowsPerStrip - 1) / sRowsPerStrip;
  sStripBytes = sRowsPerStrip * sLineBytes;
  *outNum = sNumStrips;
  *outRows = sRowsPerStrip;

  time(&curtime);
  tm = localtime(&curtime);
  sprintf(datetime, "%4d:%02d:%02d %02d:%02d:%02d", tm->tm_year + 1900, 
          tm->tm_mon, tm->tm_mday, tm->tm_hour, tm->tm_min, tm->tm_sec);
  TIFFSetField(tif, TIFFTAG_DATETIME, datetime);

  if (!inverted || sXtileSize) {
    sTmpBuf = (char *)_TIFFmalloc(sStripBytes);
    if (!sTmpBuf)
      return IIERR_MEMORY_ERR;
  }
  sLinesDone = 0;
  sAlreadyInverted = inverted;
  return 0;
}

/*
 * Write a strip of the given number.  When writing tiles, the buffer must contain the
 * whole strip of data; this routine will copy it into tiles for saving.
 */
int tiffWriteStrip(ImodImageFile *inFile, int strip, void *buf)
{
  int lines, i, xOffset, xtile, numBytes;
  size_t bufLine;
  char *inbuf;
  TIFF *tif = (TIFF *)inFile->header;
  lines = B3DMIN(sRowsPerStrip, inFile->ny - sLinesDone);

  b3dShiftBytes((unsigned char *)buf, (char *)buf, sLineBytes, lines, 1, 
                inFile->type == IITYPE_BYTE ? 1 : 0);
  if (!sXtileSize) {

    /* Ordinary strip writing */
    if (sAlreadyInverted) {
      sTmpBuf = (char *)buf;
    } else {
      for (i = 0; i < lines; i++) {
        inbuf = (char *)buf + (lines - (i + 1)) * sLineBytes;
        memcpy(sTmpBuf + i * sLineBytes, inbuf, sLineBytes);
      }
    }
    if (TIFFWriteEncodedStrip(tif, strip, sTmpBuf, sLineBytes * lines) < 0) {
      if (!sAlreadyInverted)
        _TIFFfree(sTmpBuf);
      return IIERR_IO_ERROR;
    }
  } else {
    
    /* Tile writing.  Copy just the data that exists in the strip given to us into
     the portion of the buffer corresponding to its rows and columns */
    for (xtile = 0; xtile < sNumXtiles; xtile++) {
      xOffset = xtile * sXtileSize * sPixSize;
      numBytes = B3DMIN(sXtileSize, inFile->nx - xtile * sXtileSize) * sPixSize;
      for (i = 0; i < lines; i++) {
        bufLine = sAlreadyInverted ? i : lines - (i + 1);
        inbuf = (char *)buf + bufLine * inFile->nx * sPixSize + xOffset;
        memcpy(sTmpBuf + i * sLineBytes, inbuf, numBytes);
      }

      /* But then write the entire buffer of data.  It is required to have the full 
         size */
      if (TIFFWriteEncodedTile(tif, xtile + strip * sNumXtiles, sTmpBuf, 
                               sLineBytes * sRowsPerStrip) < 0) {
        _TIFFfree(sTmpBuf);
        return IIERR_IO_ERROR;
      }
    }
  }
  b3dShiftBytes((unsigned char *)buf, (char *)buf, sLineBytes, lines, -1, 
                inFile->type == IITYPE_BYTE ? 1 : 0);
  sLinesDone += lines;
  return 0;
}

void tiffWriteFinish(ImodImageFile *inFile)
{
  inFile->state = IISTATE_BUSY;
  if (!sAlreadyInverted)
    _TIFFfree(sTmpBuf);
}

int iiTiffWriteSection(ImodImageFile *inFile, char *buf, int inSection)
{
  return (tiffWriteSectionAny(inFile, buf, inSection, 0));
}

int iiTiffWriteSectionFloat(ImodImageFile *inFile, char *buf, int inSection)
{
  return (tiffWriteSectionAny(inFile, buf, inSection, 1));
}

/* Set a description to be used on the NEXT write setup or section write */
void tiffAddDescription(const char *text)
{
  B3DFREE(sDescription);
  if (text)
    sDescription = strdup(text);
}

static int tiffWriteSectionAny(ImodImageFile *inFile, char *buf, int inSection, 
                               int ifFloat)
{
  char *useBuf = buf;
  int inverted = 0;
  int nx = inFile->nx;
  int ny = inFile->ny;
  int resolution = 0, iy;
  int compression = 1, quality = -1;

  if (inFile->padLeft || inFile->padRight) {
    b3dError(stderr, "ERROR: tiffWriteSectionAny - Cannot write from a subset of an "
             "array\n");
    return -1;
  }
  if (inSection != inFile->lastWrittenZ + 1) {
    b3dError(stderr, "ERROR: tiffWriteSectionAny - Can only write sequential sections to"
             " a TIFF file (last Z = %d, requested Z = %d)\n", inFile->lastWrittenZ, 
             inSection);
    return -1;
  }
  if (inFile->llx || inFile->lly || inFile->urx != nx -1 || inFile->ury != ny -1) {
    b3dError(stderr, "ERROR: tiffWriteSectionAny - Can only write a whole section at "
             "once\n");
    return -1;
  }
  if (inFile->format == IIFORMAT_COMPLEX) {
    b3dError(stderr, "ERROR: tiffWriteSectionAny - Cannot write complex data\n");
    return -1;
  }

  useBuf = iiMakeBufferConvertIfFloat(inFile, buf, ifFloat, &inverted,
                                      "tiffWriteSectionAny");
  if (!useBuf)
    return IIERR_MEMORY_ERR;

  if (inFile->xscale != 1.)
    resolution = 2.54e8 / inFile->xscale;
  if (getenv("IMOD_TIFF_COMPRESSION")) {
    compression = atoi(getenv("IMOD_TIFF_COMPRESSION"));
    compression = B3DMAX(1, compression);
  }
  if (getenv("IMOD_TIFF_QUALITY")) {
    quality = atoi(getenv("IMOD_TIFF_QUALITY"));
    quality = B3DMAX(1, quality);
  }
  iy = tiffWriteSection(inFile, useBuf, compression, inverted, resolution, quality);
  if (!iy)
    inFile->lastWrittenZ++;
  /* printf("Wrote TIFF %d %d\n", inFile->lastWrittenZ, inFile->file); */
  return iy;
}
  
int tiffVersion(int *minor)
{
  const char *verstrng = TIFFGetVersion();
  const char *substr;
  int update, version;
  *minor = 0;
  if (strstr(verstrng, "IMOD"))
    return 0;
  substr = strstr(verstrng, "ersion");
  if (!substr)
    return 0;
  sscanf(substr + 7, "%d.%d.%d", &version, minor, &update);
  return version;
}    

/*
 * Returns the optimal number of threads for reading the given size [nx] by [ny]
 * from a TIFF file with the compression type in [compression], and a maximum number
 * of threads given by [maxThreads].
 */
int tiffNumReadThreads(int nx, int ny, int compression, int maxThreads)
{
  int numReadThreads;
  float scale = 0.5;  /* Scale it down for unknown compression */

  /* No parallel read with no compression.  Although there is value if data are
     coming from memory, it runs slower when loading from disk/network */
  if (compression == IICOMPRESSION_NONE)
    return 1;
  if (compression == IICOMPRESSION_LZW || compression == IICOMPRESSION_ZIP)
    scale = 1.;
  else if (compression == IICOMPRESSION_JPEG)
    scale = 0.67;

  /* The logic here is 12 for K2 super-res, 8 for K2 counting, 4 for 2K, 1 for 1K */
  numReadThreads = B3DNINT(4. * scale * log(sqrt((double)nx * ny) / 944.) / log(2.));
  B3DCLAMP(numReadThreads, 1, maxThreads);
  B3DCLAMP(numReadThreads, 1, ny / 2);
  numReadThreads = numOMPthreads(numReadThreads);
  return B3DMIN(numReadThreads, MAX_TIFF_THREADS);
}

/*
 * Reads section [izRead] from a TIFF file with up to [maxThreads] threads, reading the 
 * area defined by [llx], [urx], [lly], [ury] into the buufer [readBuf].  [dataSize]
 * should have the byte size of the data elements, and [convert] should have one of
 * MRSA_NOPROC, MRSA_BYTE, MRSA_USHORT, and MRSA_FLOAT.  Returns the error from the
 * read operation.
 */
int tiffParallelRead(ImodImageFile **fileCopies, int maxThreads, int llx, int urx,
                     int lly, int ury, int dataSize, char *readBuf, int izRead, 
                     int convert)
{
  int ind, err, locErr, numThreads, ny = ury + 1 - lly, nx = urx + 1 - llx;

  /* Save the original loading limits completely and restore at the end */
  int llxSave = fileCopies[0]->llx;
  int urxSave = fileCopies[0]->urx;
  int llySave = fileCopies[0]->lly;
  int urySave = fileCopies[0]->ury;

  /* Determine right number of threads for current size and compression and assign the
     load limits to each file.  Disable multiple threads for > 2 GB buffers on Windows 
     and Mac based on results with 2 to 4 cores */
#if defined(_WIN32) || defined(__APPLE__)
  if ((float)nx * ( ny * dataSize) > 2.e9)
    maxThreads = 1;
#endif
  numThreads = tiffNumReadThreads(nx, ny, fileCopies[0]->tiffCompression, maxThreads);
  for (ind = 0; ind < numThreads; ind++) {
    fileCopies[ind]->llx = llx;
    fileCopies[ind]->urx = urx;
    fileCopies[ind]->lly = lly + ind * (ny / numThreads);
    fileCopies[ind]->ury = lly + (ind + 1) * (ny / numThreads) - 1;
    /*printf("%d %d %d\n", ind, fileCopies[ind]->lly, fileCopies[ind]->ury);*/
  }
  fileCopies[numThreads - 1]->ury = ury;
  /* Read in parallel or not */
  if (numThreads > 1) {
    err = 0;
#pragma omp parallel for num_threads(numThreads)                        \
  shared(readBuf, izRead, fileCopies, nx, lly, dataSize, convert, err)      \
  private(ind, locErr)
    for (ind = 0; ind < numThreads; ind++) {
      locErr = ReadSection(fileCopies[ind], 
                           readBuf + (fileCopies[ind]->lly - lly) * nx * dataSize,
                           izRead, convert);
      if (locErr)
        err = locErr;
    }
  } else {
    err = ReadSection(fileCopies[0], readBuf, izRead, convert);
  }
  fileCopies[0]->llx = llxSave;
  fileCopies[0]->urx = urxSave;
  fileCopies[0]->lly = llySave;
  fileCopies[0]->ury = urySave;
  return err;
}

/*
 * Parallel version of writing routine takes the same arguments as TiffWriteSection
 * will determine the number of threads and do a regular write if only one thread or
 * not appropriate compression.  It also returns a value to indicate whether it
 * actually used, or tried to use, parallel compression.
 */
int tiffParallelWrite(ImodImageFile *inFile, void *buf, int compression,
                      int inverted, int resolution, int quality, int *didParallel)
{
  ImodImageFile *tempFile[MAX_TIFF_THREADS];
  int thrNumStrips[MAX_TIFF_THREADS];
  int thrCumLines[MAX_TIFF_THREADS];
  char *thrTmpBuf[MAX_TIFF_THREADS];
  
  int linesPerFile, lastClean;
  char *sbuf, *useBuf, *inbuf;
  int numThreads, numLines, numStrips, cumLines = 0, cumStrips = 0;
  int file, i, ind = 0, psize = 1;
  char *bytep = (char *)buf;
  int strip, lines, err, numBytes, tileX = 0;
  int stripTarget = compression != IICOMPRESSION_NONE ? 16384 : 8192;
  int isBytes = inFile->type == IITYPE_BYTE ? 1 : 0;
  float overFactor = 1.2f;
  float zipScale = 1.;
  static double wallStart, compCum = 0., xferCum = 0.;

  *didParallel = 0;
  if (inFile->type == IITYPE_SHORT || inFile->type == IITYPE_USHORT)
    psize = 2;
  if (inFile->type == IITYPE_INT || inFile->type == IITYPE_UINT || 
      inFile->type == IITYPE_FLOAT)
    psize = 4;
  if (inFile->format == IIFORMAT_RGB)
    psize *= 3;
  sLineBytes = inFile->nx * psize;

  sRowsPerStrip =  B3DMAX(1, (stripTarget + sLineBytes / 2) / sLineBytes);

  /* The logic here is 16 for K2 super-res, 8 for K2 counting, 4 for 2K, 2 for 1K
     And ZIP gets a scaling ranging from 1 at 2K to 1.5 at 4K to 2 at 8K
     Note that execution time is linear in # of bytes, but thread efficiency is not */
  if (compression == IICOMPRESSION_ZIP) {
    zipScale = 1. + log(sqrt((double)inFile->nx * inFile->ny) / 2048) / log(4.);
    B3DCLAMP(zipScale, 1., 2.);
  }
  numThreads = B3DNINT(zipScale * 2. * sqrt((double)inFile->nx * inFile->ny) / 944.);
  B3DCLAMP(numThreads, 1, MAX_TIFF_THREADS);
  B3DCLAMP(numThreads, 1, (inFile->ny / sRowsPerStrip) / 4);
  numThreads = numOMPthreads(numThreads);
  numThreads = B3DMIN(numThreads, MAX_TIFF_THREADS);

  /* Do regular write if it does not qualify for parallel */
  if (numThreads < 2 || inFile->format == IIFORMAT_RGB || tiffVersion(&err) < 4 ||
      compression == IICOMPRESSION_JPEG || compression == IICOMPRESSION_NONE)
    return tiffWriteSection(inFile, buf, compression, inverted, resolution, quality);
           
  /* Set up, with flag not to delete description */
  *didParallel = 1;
  sSettingUpParallel = -1;
  err = tiffWriteSetup(inFile, compression, inverted, resolution, quality,
                       &strip, &lines, &tileX);
  if (!inverted)
    _TIFFfree(sTmpBuf);
  sSettingUpParallel = 0;
  if (err)
    return err;

  /* Set up the division of strips into files, all the extra on the last one
     Compute number of strips in last file and allow buffer size based on that */
  numStrips = sNumStrips;
  linesPerFile = inFile->ny / numThreads;
  linesPerFile = sRowsPerStrip * (linesPerFile / sRowsPerStrip);
  ind = numStrips - (numThreads - 1) * linesPerFile / sRowsPerStrip;
  sFileBufSize = 4096 + (overFactor * sStripBytes + 8) * ind;

  /* Set up everything for each file */
  for (file = 0; file < numThreads; file++) {
    numLines = B3DCHOICE(file == numThreads - 1, inFile->ny - cumLines, linesPerFile);
    lastClean = file;
    thrTmpBuf[file] = NULL;
    sFileBuf[file] = (char *)_TIFFmalloc(sFileBufSize);
    sMaxBufInd[file] = sCurBufInd[file] = 0;
    tempFile[file] = iiNew();
    if (!sFileBuf[file] || !tempFile[file]) {
      err = IIERR_MEMORY_ERR;
      break;
    }

    /* It needs a filename, although it turns out to be irrelevant */
    tempFile[file]->filename = (char *)malloc(strlen(inFile->filename) + 16);
    if (!tempFile[file]->filename) {
      err = IIERR_MEMORY_ERR;
      break;
    }
    sprintf(tempFile[file]->filename, "%s.%d.%d", inFile->filename, imodGetpid(), file);

    /* Set up other properties including the proper number of lines for this file */
    tempFile[file]->nx = inFile->nx;
    tempFile[file]->ny = numLines;
    tempFile[file]->file = IIFILE_TIFF;
    tempFile[file]->mode = inFile->mode;
    tempFile[file]->type = inFile->type;
    tempFile[file]->nz = 1;

    /* Open new file: it will call the ClientOpen with this value as client data */
    sSettingUpParallel = file + 1;
    err = tiffOpenNew(tempFile[file]);
    if (err)
      break;

    /* Do the setup of the file, and get the number of strips */
    tileX = 0;
    err = tiffWriteSetup(tempFile[file], compression, inverted, resolution, quality,
                       &strip, &lines, &tileX);
    if (err)
      break;
    thrNumStrips[file] = sNumStrips;
    thrCumLines[file] = B3DCHOICE(inverted, cumLines, inFile->ny - cumLines - numLines);
    if (!inverted)
      thrTmpBuf[file] = sTmpBuf;
    cumStrips += sNumStrips;
    cumLines += numLines;
  }
  sSettingUpParallel = 0;
  if (!err && numStrips != cumStrips) {
    b3dError(stderr, "tiffParallelWrite: number of strips for full file (%d) not equal "
             "to total for temp files (%d)\n", numStrips, cumStrips);
    err = IIERR_IO_ERROR;
  }

  wallStart = wallTime();
  if (!err) {

  /* Do the compression in parallel to the various files */
#pragma omp parallel for num_threads(numThreads)                        \
  shared(thrNumStrips, tempFile, sRowsPerStrip, inverted, bytep, thrCumLines) \
  shared(sLineBytes, isBytes)                                           \
  private(file, sLinesDone, strip, lines, sbuf, useBuf, i, inbuf)
    for (file = 0; file < numThreads; file++) {
      sLinesDone = 0;
      for (strip = 0; strip < thrNumStrips[file]; strip++) {

        /* Do all the needed operations for a strip: shift bytes, invert lines, write,
           unshift the bytes */
        lines = B3DMIN(sRowsPerStrip, tempFile[file]->ny - sLinesDone);
        if (inverted)
          sbuf = bytep + strip * sStripBytes + thrCumLines[file] * sLineBytes;
        else
          sbuf = bytep + (thrCumLines[file] + tempFile[file]->ny - (sLinesDone + lines))
            * (int)sLineBytes;
        b3dShiftBytes((unsigned char *)sbuf, sbuf, sLineBytes, lines, 1, isBytes);
        useBuf = sbuf;
        if (!inverted) {
          useBuf = thrTmpBuf[file];
          fflush(stdout);
          for (i = 0; i < lines; i++) {
            inbuf = sbuf + (lines - (i + 1)) * sLineBytes;
            memcpy(useBuf + i * sLineBytes, inbuf, sLineBytes);
          }
        }
        i = TIFFWriteEncodedStrip((TIFF *)tempFile[file]->header, strip, useBuf,
                                  sLineBytes * lines);
        b3dShiftBytes((unsigned char *)sbuf, sbuf, sLineBytes, lines, -1, isBytes);
        if (i < 0) {
          err = IIERR_IO_ERROR;
          break;
        }
        sLinesDone += lines;
      }
    }
  }
  
  compCum += wallTime() - wallStart;
  wallStart = wallTime();
  if (!err) {

    /* Transfer the raw strips to the real file; allow a generous strip buffer */
    sTmpBuf = (char *)_TIFFmalloc(2 * sStripBytes);
    cumStrips = 0;
    sLinesDone = 0;
    for (file = 0; file < numThreads; file++) {
      for (strip = 0; strip < thrNumStrips[file]; strip++) {
        lines = B3DMIN(sRowsPerStrip, inFile->ny - sLinesDone);
        numBytes = (int)TIFFReadRawStrip((TIFF *)tempFile[file]->header, strip, sTmpBuf, 
                                         2 * sStripBytes);
        if (numBytes <= 0) {
          b3dError(stderr, "tiffParallelWrite: Read error getting raw strip %d from "
                   "file %d\n", strip, file);
          err = IIERR_IO_ERROR;
          break;
        }
        err = (int)TIFFWriteRawStrip((TIFF *)inFile->header, strip + cumStrips, sTmpBuf,
                                     numBytes);
        if (err <= 0) {
          b3dError(stderr, "tiffParallelWrite: Error rewriting raw strip %d from "
                   "file %d\n", strip, file);
          err = IIERR_IO_ERROR;
          break;
        }
        sLinesDone += lines;
      }
      cumStrips += thrNumStrips[file];
      err = 0;
    }
    _TIFFfree(sTmpBuf);
  }
  
  xferCum += wallTime() - wallStart;
  /* printf("compression %.3f  copying %.3f\n", compCum, xferCum); */

  /* Clean up works for normal case and any errors */
  for (ind = 0; ind <= lastClean; ind++) {
    if (tempFile[ind]) {
      iiClose(tempFile[ind]);
      remove(tempFile[ind]->filename);
      iiDelete(tempFile[ind]);
      if (!inverted && thrTmpBuf[ind])
        _TIFFfree(thrTmpBuf[ind]);
    }
    if (sFileBuf[ind])
      _TIFFfree(sFileBuf[ind]);
  }
  inFile->state = IISTATE_BUSY;

  /* If there was an error, leave the description in case caller wants to try regular 
     save */
  if (!err)
    B3DFREE(sDescription);
  return err;
}

/*
 * Put out the min and max values: but keep them within legal range for the mode
 * since the library will simply convert them and out-of-range values will wrap around
 */
static void constrainAndStoreMinMax(ImodImageFile *inFile)
{
  double dmin, dmax;
  TIFF *tif = (TIFF *)inFile->header;
  dmin = inFile->amin;
  dmax = inFile->amax;
  if (inFile->format == IIFORMAT_RGB) {
    dmin = 0.;
    dmax = 255.;
  } else {
    switch (inFile->type) {
    case IITYPE_BYTE:
      B3DCLAMP(dmin, 0., 255.);
      B3DCLAMP(dmax, 0., 255.);
      dmin -= 128.;
      dmax -= 128.;
      break;
    case IITYPE_UBYTE:
      B3DCLAMP(dmin, 0., 255.);
      B3DCLAMP(dmax, 0., 255.);
      break;
    case IITYPE_SHORT:
      B3DCLAMP(dmin, -32768., 32767.);
      B3DCLAMP(dmax, -32768., 32767.);
      break;
    case IITYPE_USHORT:
      B3DCLAMP(dmin, 0., 65535.);
      B3DCLAMP(dmax, 0., 65535.);
      break;
    }
  }
  TIFFSetField(tif, TIFFTAG_SMINSAMPLEVALUE, dmin);
  TIFFSetField(tif, TIFFTAG_SMAXSAMPLEVALUE, dmax);
}
