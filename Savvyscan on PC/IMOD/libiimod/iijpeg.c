/*
 *    iitif.c    - specific routines for jpeg-type ImodImageFile's
 *
 *    Author:  James David Mastronarde
 *
 *   Copyright (C) 1995-2016 by the Regents of the University of Colorado.
 *
 *  $Id$
 */

#include "b3dutil.h"
#include "jpeglib.h"
#include <setjmp.h>
#include "iimage.h"

struct my_error_mgr {
  struct jpeg_error_mgr pub;  /* "public" fields */

  jmp_buf setjmp_buffer;  /* for return to caller */
};

typedef struct my_error_mgr * my_error_ptr;

static struct my_error_mgr sJerr;

static void jpegDelete(ImodImageFile *inFile);
static int jpegReadSectionByte(ImodImageFile *inFile, char *buf, int inSection);
static int jpegReadSectionUShort(ImodImageFile *inFile, char *buf, int inSection);
static int jpegReadSectionFloat(ImodImageFile *inFile, char *buf, int inSection);
static int jpegReadSection(ImodImageFile *inFile, char *buf, int inSection);
static int jpegReadSectionAny(ImodImageFile *inFile, char *buf, int inSection, int type);
static int iiJpegWriteSection(ImodImageFile *inFile, char *buf, int inSection);
static int iiJpegWriteSectionFloat(ImodImageFile *inFile, char *buf, int inSection);
static int iiJpegWriteSectionAny(ImodImageFile *inFile, char *buf, int inSection,
                                 int ifFloat);


/*
 * Here's the routine that will replace the standard error_exit method:
 */
METHODDEF(void) my_error_exit (j_common_ptr cinfo)
{
  /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
  my_error_ptr myerr = (my_error_ptr) cinfo->err;

  /* Postponed message until after returning */
  /* Return control to the setjmp point */
  longjmp(myerr->setjmp_buffer, 1);
}

/*
 * Check for and open a JPEG file
 */
int iiJPEGCheck(ImodImageFile *inFile)
{
  FILE *fp;
  unsigned char buf[4];
  j_decompress_ptr cinfoPtr;
  char messBuf[JMSG_LENGTH_MAX + 10];

  if (!inFile)
    return IIERR_BAD_CALL;
  fp = inFile->fp;
  if (!fp)
    return IIERR_BAD_CALL;

  /* Look for the magic numbers at start */
  rewind(fp);
  if (fread(&buf, 1, 4, fp) < 4) {
    b3dError(stderr, "ERROR: iiJPEGCheck - Reading file %s\n", inFile->filename);
    return IIERR_IO_ERROR;
  } 
  if (buf[0] != 0xFF || buf[1] != 0xD8 || buf[2] != 0xFF)
    return IIERR_NOT_FORMAT;

  /* Create decompression object and set up error handling */
  cinfoPtr = (j_decompress_ptr)malloc(sizeof(struct jpeg_decompress_struct));
  if (!cinfoPtr) {
    b3dError(stderr, "ERROR: iiJPEGCheck - Allocating decompression structure\n");
    return IIERR_MEMORY_ERR;
  }

  cinfoPtr->err = jpeg_std_error(&sJerr.pub);
  sJerr.pub.error_exit = my_error_exit;

  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(sJerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object
     */
    (*cinfoPtr->err->format_message) ((j_common_ptr)cinfoPtr, messBuf);
    b3dError(stderr, "iiJPEGCheck: JPEG warning/error %s\n", messBuf);
    jpeg_destroy_decompress(cinfoPtr);
    free(cinfoPtr);
    return IIERR_IO_ERROR;
  }

  /* Create object and read header to get properties */
  jpeg_create_decompress(cinfoPtr);

  rewind(fp);
  jpeg_stdio_src(cinfoPtr, fp);
  jpeg_read_header(cinfoPtr, TRUE);
  
  inFile->nx = cinfoPtr->image_width;
  inFile->ny = cinfoPtr->image_height;
  inFile->nz = 1;
  inFile->type = IITYPE_UBYTE;
  if (cinfoPtr->jpeg_color_space == JCS_GRAYSCALE) {
    inFile->format = IIFORMAT_LUMINANCE;
    inFile->mode = MRC_MODE_BYTE;
  } else if (cinfoPtr->num_components > 2) {
    inFile->format = IIFORMAT_RGB;
    inFile->mode = MRC_MODE_RGB;
  } else {
    b3dError(stderr, "ERROR: iiJPEGCheck - JPEG colorspace not GRAYSCALE and # of "
             "components is %d\n", cinfoPtr->num_components);
    jpeg_destroy_decompress(cinfoPtr);
    free(cinfoPtr);
    return IIERR_NOT_FORMAT;
  }

  /* Pixel size smaller than 400 nm cannot be encoded in the 16-bit integers of the JFIF
     so skip setting a pixel size */

  /* Set up the rest ofteh basic stuff and pointers for reading routines */
  inFile->amin = 0.;
  inFile ->amax = 255.;
  inFile->amean = 128.;
  inFile->file = IIFILE_JPEG;
  inFile->header = (char *)cinfoPtr;
  inFile->cleanUp = jpegDelete;
  inFile->fillMrcHeader = iiSimpleFillMrcHeader;

  inFile->readSection = jpegReadSection;
  inFile->readSectionUShort = jpegReadSectionUShort;
  inFile->readSectionByte = jpegReadSectionByte;
  inFile->readSectionFloat = jpegReadSectionFloat;

  /* We need abort the object and redo header later because some callers might close
   and reopen the file */
  jpeg_abort_decompress(cinfoPtr);
  return 0;
}

/*
 * When a file is being deleted, that is when to destroy and free the compression object 
 */
static void jpegDelete(ImodImageFile *inFile)
{
  if (inFile->header) {
    if (inFile->readSectionByte)
      jpeg_destroy_decompress((j_decompress_ptr)inFile->header);
    else
      jpeg_destroy_compress((j_compress_ptr)inFile->header);
  }
  B3DFREE(inFile->header);
}

/*
 * The wrappers to the reading routine
 */
static int jpegReadSectionByte(ImodImageFile *inFile, char *buf, int inSection)
{ 
  return(jpegReadSectionAny(inFile, buf, inSection, MRSA_BYTE));
}

static int jpegReadSectionUShort(ImodImageFile *inFile, char *buf, int inSection)
{ 
  return(jpegReadSectionAny(inFile, buf, inSection, MRSA_USHORT));
}

static int jpegReadSectionFloat(ImodImageFile *inFile, char *buf, int inSection)
{ 
  return(jpegReadSectionAny(inFile, buf, inSection, MRSA_FLOAT));
}

static int jpegReadSection(ImodImageFile *inFile, char *buf, int inSection)
{
  return(jpegReadSectionAny(inFile, buf, inSection, MRSA_NOPROC));
}

/*
 * The main reading routine.
 */
static int jpegReadSectionAny(ImodImageFile *inFile, char *buf, int inSection, int type)
{
  int pixSizeBuf[4] = {0, 1, 4, 2};
  char messBuf[JMSG_LENGTH_MAX + 10];
  MrcHeader hdata;
  LineProcData d;
  IloadInfo loadInfo;
  IloadInfo *li = &loadInfo;
  int padLeft, padRight;
  int  ny = inFile->ny;
  int yEnd, err = -1, jpegLine = 0, tval;
  int freeMap = 0;
  unsigned char *tmpData = NULL;
  j_decompress_ptr cinfoPtr = (j_decompress_ptr)inFile->header;
  
  if (!inFile->readSection) {
    b3dError(stderr, "ERROR: jpegReadSectionAny - Trying to read from newly created"
             " JPEG file\n");
    return IIERR_BAD_CALL;
  }

  if (inSection) {
    b3dError(stderr, "ERROR: jpegReadSectionAny - Trying to read section %d; only section"
             "0 can be read from JPEG file\n", inSection);
    return IIERR_BAD_CALL;
  }

  tmpData = B3DMALLOC(unsigned char, inFile->nx * 
                      (inFile->format == IIFORMAT_RGB ? 3 : 1));
  if (!tmpData) {
    b3dError(stderr, "ERROR: jpegReadSectionAny - Allocating line buffer for reading "
             "JPEG file\n");
    return IIERR_MEMORY_ERR;
  }

  /* Initialize any variables used in the error block to make compiler happy; set up
     the error function again */
  d.map = NULL;
  d.yStart = 0;
  cinfoPtr->err = jpeg_std_error(&sJerr.pub);
  sJerr.pub.error_exit = my_error_exit;

  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(sJerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object
     * And do everything for a normal return
     */
    if (err != IIERR_QUITTING) {
      err = (jpegLine < d.yStart && !err) ? 0 : IIERR_IO_ERROR;
      if (err) {
        (*cinfoPtr->err->format_message) ((j_common_ptr)cinfoPtr, messBuf);
        b3dError(stderr, "jpegReadSectionAny: JPEG error - %s\n", messBuf);
      }
    }
    jpeg_abort_decompress(cinfoPtr);
    B3DFREE(tmpData);
    if (freeMap)
      free(d.map);
    return err;
  }

  /* Reestablish the object, data source and header */
  rewind(inFile->fp);
  jpeg_stdio_src(cinfoPtr, inFile->fp);
  jpeg_read_header(cinfoPtr, TRUE);

  /* Set output properties if not grayscale */
  if (inFile->format == IIFORMAT_RGB)
    cinfoPtr->out_color_space = JCS_RGB;
 
  /* Translate the information to a loadInfo to call common routines and set some
     type-dependent settings as in iimrc; setup MRC header too */
  iiSimpleFillMrcHeader(inFile, &hdata);
  iiMRCsetLoadInfo(inFile, li);
  padLeft = B3DMAX(0, li->padLeft);
  padRight = B3DMAX(0, li->padRight);
  yEnd = li->ymax;
  if (type == MRSA_FLOAT || !type) {
    li->outmin = inFile->smin;
    li->outmax = inFile->smax;
  } else {
    li->outmin = 0;
    li->outmax = (type == MRSA_USHORT) ? 65535 : 255;
  }

  /* Initialize members of data structure */
  d.type = type;
  d.readY = 0;
  d.cz = 0;
  d.swapped = 0;
  err = iiInitReadSectionAny(&hdata, li, (unsigned char *)buf, &d, &freeMap, &yEnd,
                             "jpegReadSectionAny");
  if (err)
    return err;

  /* Fill in some missing pieces and adjust pointers/indexes for inversion
     No need to invert yStart/yEnd because they are used as limits for line numbers
     counting from the true bottom */
  d.xDimension = d.xsize + padLeft + padRight;
  d.needData = 1;
  pixSizeBuf[0] = d.pixSize;
  tval = (yEnd - d.yStart) * d.xDimension;
  d.pixIndex += tval;
  d.bufp += pixSizeBuf[type] * tval;
  d.usbufp += tval;
  d.fbufp += tval;
  d.deltaYsign = -1;

  /* Loop from the top of the image; it is safe to stop when the desired lines are
     obtained as the error is suppressed above */
  jpegLine = ny - 1;
  jpeg_start_decompress(cinfoPtr);
  err = 0;
  while (jpegLine >= 0 && jpegLine >= d.yStart) {
    jpeg_read_scanlines(cinfoPtr, &tmpData, 1);
    if (jpegLine <= yEnd) {
      d.bdata = tmpData + d.xStart * d.pixSize;
      err = iiProcessReadLine(&hdata, li, &d);
      if (err)
        break;
    }
    jpegLine--;
  }

  jpeg_finish_decompress(cinfoPtr);
  jpeg_abort_decompress(cinfoPtr);
  B3DFREE(tmpData);
  if (freeMap)
    free(d.map);
  return err;
}

/*!
 * Open a new JPEG file: just setup function pointers and allocate compression object
 */
int jpegOpenNew(ImodImageFile *inFile)
{
  inFile->state = IISTATE_READY;
  inFile->cleanUp = jpegDelete; 
  inFile->fillMrcHeader = iiSimpleFillMrcHeader;
  inFile->writeSection = iiJpegWriteSection;
  inFile->writeSectionFloat = iiJpegWriteSectionFloat;
  inFile->fp = fopen(inFile->filename, "wb");
  if (!inFile->fp)
    return IIERR_IO_ERROR;
  inFile->header = (char *)malloc(sizeof(struct jpeg_compress_struct));
  if (!inFile->header) {
    fclose(inFile->fp);
    return IIERR_MEMORY_ERR;
  }
  return 0;
}

/*
 * Wrappers to the real writing routine
 */
static int iiJpegWriteSection(ImodImageFile *inFile, char *buf, int inSection)
{
  return iiJpegWriteSectionAny(inFile, buf, inSection, 0);
}

static int iiJpegWriteSectionFloat(ImodImageFile *inFile, char *buf, int inSection)
{
  return iiJpegWriteSectionAny(inFile, buf, inSection, 1);
}

/*
 * Wrapper that handles the iimage calls in, checks for possible bad things in a generic
 * call, and sets quality and resolution from environment variables
 */
static int iiJpegWriteSectionAny(ImodImageFile *inFile, char *buf, int inSection,
                                 int ifFloat)
{
  int quality = -1, resolution = 0;
  int inverted = 0;
  char *useBuf;
  int err;

  if (inSection) {
    b3dError(stderr, "ERROR: iiJpegWriteSectionAny - Trying to write section %d"
             " to a JPEG file; only 0 is allowed\n", inSection);
    err = IIERR_BAD_CALL;
  }
  if (!err && (inFile->padLeft || inFile->padRight)) {
    b3dError(stderr, "ERROR: iiJpegWriteSectionAny - Cannot write from a subset of an "
             "array\n");
    err = -1;
  }
  if (!err && !((inFile->format == IIFORMAT_LUMINANCE && inFile->mode == MRC_MODE_BYTE) ||
        (inFile->format == IIFORMAT_RGB && !ifFloat))) {
    b3dError(stderr, "ERROR: iiJpegWriteSectionAny - %s\n", 
             (inFile->format == IIFORMAT_RGB && ifFloat) ? 
             "Cannot write float data to an RGB JPEG file" :
             "File mode must be byte or RGB to write to a JPEG file");
    err = -1;
  }
  if (!err && (inFile->llx || inFile->lly || 
    (inFile->urx != -1 && inFile->urx != inFile->nx - 1) || 
    (inFile->ury != -1 && inFile->ury != inFile->ny - 1))) {
    b3dError(stderr, "ERROR: iiJpegWriteSectionAny - Can only write a whole section at "
             "once\n");
    err = -1;
  }

  /* Convert floats if necessary */
  if (!err) {
    useBuf = iiMakeBufferConvertIfFloat(inFile, buf, ifFloat, &inverted,
                                        "iiJpegWriteSectionAny");
    if (!useBuf)
      err = IIERR_MEMORY_ERR;
  }

  /* Free the header, it has not been initialized yet, just created */
  if (err) {
    B3DFREE(inFile->header);
    return err;
  }

  /* Get environment variable values */
  if (getenv("IMOD_JPEG_QUALITY")) {
    quality = atoi(getenv("IMOD_JPEG_QUALITY"));
    quality = B3DCLAMP(quality, 1, 100);
  }
  if (getenv("IMOD_JPEG_RESOLUTION")) {
    err = atoi(getenv("IMOD_JPEG_RESOLUTION"));
    if (err > 65535)
      b3dError(stderr, "WARNING: iiJpegWriteSectionAny - IMOD_JPEG_RESOLUTION is set "
               "too high to be stored in 16-bit field of JPEG file\n");
    else
      resolution = err;
    resolution = B3DMAX(1, resolution);
  }

  /* Call through central routine */
  err = jpegWriteSection(inFile, useBuf, inverted, resolution, quality);
  if (!err)
    inFile->lastWrittenZ = 0;
  if (useBuf != buf)
    free(useBuf);
  return err;
}

/*!
 * Writes the section in [buf] to a JPEG file open on [inFile], which must have a mode of
 * RGB or BYTE.  Set [inverted] to 1 if line order is already inverted in Y (first line 
 * at top); [resolution] to a value up to 65535 in dots per inch or 0 for none, and
 * [quality] to a value from 1 to 100, or -1 for no setting.
 */
int jpegWriteSection(ImodImageFile *inFile, char *buf, int inverted, int resolution, 
                     int quality)
{
  char messBuf[JMSG_LENGTH_MAX + 10];
  j_compress_ptr cinfoPtr = (j_compress_ptr)inFile->header;
  int iy, useY, err = 0;
  unsigned char *rowPointer;
  
  /* Check various bad things */
  if (!inFile->writeSection) {
    b3dError(stderr, "ERROR: jpegWriteSection - Trying to write to an existing"
             " JPEG file\n");
    err = IIERR_BAD_CALL;
  }
  if (!err && !inFile->lastWrittenZ) {
    b3dError(stderr, "ERROR: jpegWriteSection - Trying to write more than one section"
             " to a JPEG file\n");
    err = IIERR_BAD_CALL;
  }
  if (!err && inFile->mode != MRC_MODE_BYTE && inFile->mode != MRC_MODE_RGB) {
    b3dError(stderr, "ERROR: jpegWriteSection - Mode for writing a JPEG file must be "
             "either byte or RGB; it is %d\n", inFile->mode);
    err = IIERR_BAD_CALL;
  }

  /* Free the header, it has not been initialized yet, just created */
  if (err) {
    B3DFREE(inFile->header);
    return err;
  }

  /* Set up error handling */
  rewind(inFile->fp);
  cinfoPtr->err = jpeg_std_error(&sJerr.pub);
  sJerr.pub.error_exit = my_error_exit;

  /* Establish the setjmp return context for my_error_exit to use. */
  if (setjmp(sJerr.setjmp_buffer)) {
    /* If we get here, the JPEG code has signaled an error.
     * We need to clean up the JPEG object
     */
    (*cinfoPtr->err->format_message) ((j_common_ptr)cinfoPtr, messBuf);
    b3dError(stderr, "jpegWriteSection: JPEG error - %s\n", messBuf);
    jpeg_abort_compress(cinfoPtr);
    return IIERR_IO_ERROR;
  }

  /* Create the compression object and set the basic properties as well as optional
     quality and resolution */
  jpeg_create_compress(cinfoPtr);
  jpeg_stdio_dest(cinfoPtr, inFile->fp);
  cinfoPtr->image_width = inFile->nx;
  cinfoPtr->image_height = inFile->ny;
  cinfoPtr->input_components = inFile->mode == MRC_MODE_RGB ? 3 : 1;
  cinfoPtr->in_color_space = inFile->mode == MRC_MODE_RGB ? JCS_RGB : JCS_GRAYSCALE;
  jpeg_set_defaults(cinfoPtr);
  if (quality > 0)
    jpeg_set_quality(cinfoPtr, B3DMIN(quality, 100), TRUE);
  if (resolution > 0) {
    cinfoPtr->write_JFIF_header = TRUE;
    cinfoPtr->density_unit = 1;
    cinfoPtr->X_density = cinfoPtr->Y_density = resolution;
  }

  /* Do the compression line by line */
  jpeg_start_compress(cinfoPtr, TRUE);
  for (iy = inFile->ny - 1; iy >= 0; iy--) {
    useY = B3DCHOICE(inverted, (inFile->ny - 1) - iy, iy);
    rowPointer = (unsigned char *)buf + cinfoPtr->input_components * useY * inFile->nx;
    jpeg_write_scanlines(cinfoPtr, &rowPointer, 1);
  }

  jpeg_finish_compress(cinfoPtr);
  jpeg_abort_compress(cinfoPtr);

  return 0;
}
