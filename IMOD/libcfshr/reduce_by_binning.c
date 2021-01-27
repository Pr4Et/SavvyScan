/* reduce_by_binning.c: has functions for extracting with or without binning, and
 * binning into a slice.
 *
 * $Id$
 *
 */

#include "imodconfig.h"
#include "mrcslice.h"

#ifdef F77FUNCAP
#define reduce_by_binning REDUCE_BY_BINNING
#define  binintoslice BININTOSLICE
#define extractwithbinning EXTRACTWITHBINNING
#define irepak IREPAK
#else
#ifdef G77__HACK
#define reduce_by_binning reduce_by_binning__
#else
#define reduce_by_binning reduce_by_binning_
#endif
#define extractwithbinning extractwithbinning_
#define irepak irepak_
#define  binintoslice binintoslice_
#endif

/*!
 * Extracts a portion of an array into a portion of another array with optional summing 
 * of pixels (binning) by the factor [nbin].  [array] has the input, with X dimension 
 * [nxDim], and the range of coordinates to extract is given in [xStart], [xEnd], 
 * [yStart], and [yEnd] * (inclusive, numbered from 0).  The slice mode is in
 * [type], and can be byte, short, unsigned short, float, or RGB.
 * [brray] receives the output, with the output size returned in [nxr] and [nyr].  
 * [nxBdim] can specify the X dimension of [brray] or be 0 to have data packed 
 * contigously (X dimension equals X size); [bxOffset] and [byOffset] can indicate X and 
 * Y offsets into [brray].  The output will have the same mode as the input and will be 
 * the average of binned values, except in two cases.  If the input is bytes or RGB and 
 * [keepByte] is 0 or -1, then the output will be signed short integers and will be 
 * the sum or the average, respectively, of the binned values, and for RGB data the output
 * will be an equal sum or average of red, green, and blue values.  [brray] can be the 
 * same as [array], provided that its X dimension is not bigger than the input size,
 * unless the binning is 1, the input is bytes and [keepByte] is not greater than 0. 
 * The output size is obtained by integer division of the input size 
 * by the binning.  If the remainder of this division is nonzero, the data are 
 * centered in the output array as nearly as possible.  Specifically, the
 * coordinates of the lower left corner of the output array are offset by 
 * ^  ((nx % nbin) / 2, (ny % nbin) / 2) ^
 * relative to the input array.  Returns 1 for an unsupported data type, 2 for an
 * input X range bigger than [nxDim], or 3 for an output range bigger than the output
 * X dimension plus offset.
 */
int extractAndBinIntoArray(void *array, int type, int nxDim, int xStart, int xEnd,
                           int yStart, int yEnd, int nbin, void *brray, int nxBdim,
                           int bxOffset, int byOffset, int keepByte, int *nxr, int *nyr)
{
  int i, j;
  int nxin = xEnd + 1 - xStart;
  int nyin = yEnd + 1 - yStart;
  int nbinsq = nbin * nbin;
  int nxout = nxin / nbin;
  int nyout = nyin / nbin;
  int ixofs = (nxin % nbin) / 2;
  int iyofs = (nyin % nbin) / 2;
  int bin234OK = (type == SLICE_MODE_RGB || (type == SLICE_MODE_BYTE && keepByte < 0)) ?
    0 : 1;
  int sum, ix, red, green, blue, pixSize, outPixSize;
  size_t iy;
  b3dFloat fsum;
  unsigned char *bdata = (unsigned char *)brray;
  b3dInt16 *sdata = (b3dInt16 *)brray;
  b3dUInt16 *usdata = (b3dUInt16 *)brray;
  b3dInt16 *sline1, *sline2, *sline3, *sline4;
  b3dUInt16 *usline1, *usline2, *usline3, *usline4;
  unsigned char *cline1, *cline2, *cline3, *cline4;
  b3dFloat *fdata = (b3dFloat *)brray;
  b3dFloat *fline1, *fline2, *fline3, *fline4;

  if (type != SLICE_MODE_BYTE && type != SLICE_MODE_SHORT && 
      type != SLICE_MODE_USHORT && type != SLICE_MODE_FLOAT && 
      type != SLICE_MODE_RGB)
    return 1;
  if (!nxBdim)
    nxBdim = nxout;
  if (nxin > nxDim)
    return 2;
  if (nxout + bxOffset > nxBdim)
    return 3;

  /* Get the bytes per pixel */
  if (type == SLICE_MODE_BYTE)
    pixSize = 1;
  else if (type == SLICE_MODE_FLOAT)
    pixSize = 4;
  else if (type == SLICE_MODE_RGB)
    pixSize = 3;
  else
    pixSize = 2;
  outPixSize = pixSize;
  if ((type == SLICE_MODE_BYTE || type == SLICE_MODE_RGB) && keepByte <= 0)
    outPixSize = 2;

  /* Advance array pointer to the beginning of the data and incorporate Y offset into X */
  array = (void *)((unsigned char *)array + (xStart + (size_t)yStart * nxDim) * pixSize);
  bxOffset += byOffset * nxBdim;

  if (nbin == 1 && pixSize == outPixSize) {

    /* Binning 1: copy the data line by line */
    for (iy = 0; iy < nyout; iy++) {
      cline1 = ((unsigned char *)array) + pixSize * iy * nxDim;
      bdata = ((unsigned char *)brray) + pixSize * (bxOffset + iy * nxBdim);
      for (ix = 0; ix < nxout * pixSize; ix++)
        *bdata++ = *cline1++;
    }

  } else if (nbin == 2 && bin234OK) {
    
    /* Binning by 2 */
    switch (type) {
    case SLICE_MODE_BYTE:
      for (iy = 0; iy < nyout; iy++) {
        cline1 = ((unsigned char *)array) + 2 * iy * nxDim;
        cline2 = cline1 + nxDim;
        if (keepByte) {
          bdata = ((unsigned char *)brray) + bxOffset + iy * nxBdim;
          for (ix = 0; ix < nxout; ix++) {
            sum = 2 + *cline1 + *(cline1 + 1) + *cline2 + *(cline2 + 1);
            *bdata++ = (unsigned char)(sum / 4);
            cline1 += 2;
            cline2 += 2;
          }
        } else {
          sdata = ((b3dInt16 *)brray) + bxOffset + iy * nxBdim;
          for (ix = 0; ix < nxout; ix++) {
            *sdata++ = (b3dInt16)(*cline1) + *(cline1 + 1) + *cline2 + 
              *(cline2 + 1);
            cline1 += 2;
            cline2 += 2;
          }
        }
      }
      break;
      
    case SLICE_MODE_SHORT:
      for (iy = 0; iy < nyout; iy++) {
        sline1 = ((b3dInt16 *)array) + 2 * iy * nxDim;
        sline2 = sline1 + nxDim;
        sdata = ((b3dInt16 *)brray) + bxOffset + iy * nxBdim;
        for (ix = 0; ix < nxout; ix++) {
          sum = 2 + *sline1 + *(sline1 + 1) + *sline2 + *(sline2 + 1);
          sline1 += 2;
          sline2 += 2;
          *sdata++ = sum / 4;
        }
      }
      break;

    case SLICE_MODE_USHORT:
      for (iy = 0; iy < nyout; iy++) {
        usline1 = ((b3dUInt16 *)array) + 2 * iy * nxDim;
        usline2 = usline1 + nxDim;
        usdata = ((b3dUInt16 *)brray) + bxOffset + iy * nxBdim;
        for (ix = 0; ix < nxout; ix++) {
          sum = 2 + *usline1 + *(usline1 + 1) + *usline2 + *(usline2 + 1);
          usline1 += 2;
          usline2 += 2;
          *usdata++ = sum / 4;
        }
      }
      break;

    case SLICE_MODE_FLOAT:
      for (iy = 0; iy < nyout; iy++) {
        fline1 = ((b3dFloat *)array) + 2 * iy * nxDim;
        fline2 = fline1 + nxDim;
        fdata = ((float *)brray) + bxOffset + iy * nxBdim;
        for (ix = 0; ix < nxout; ix++) {
          fsum = *fline1 + *(fline1 + 1) + *fline2 + *(fline2 + 1);
          fline1 += 2;
          fline2 += 2;
          *fdata++ = fsum / 4.f;
        }
      }
      break;
    }
  
  } else if (nbin == 3 && bin234OK) {

    /* Binning by 3 */
    switch (type) {
    case SLICE_MODE_BYTE:
      for (iy = 0; iy < nyout; iy++) {
        cline1 = ((unsigned char *)array) + (3 * iy + iyofs) * nxDim + ixofs;
        cline2 = cline1 + nxDim;
        cline3 = cline2 + nxDim;
        if (keepByte) {
          bdata = ((unsigned char *)brray) + bxOffset + iy * nxBdim;
          for (ix = 0; ix < nxout; ix++) {
            sum = 4 + *cline1 + *(cline1 + 1) + *(cline1 + 2) +
              *cline2 + *(cline2 + 1) + *(cline2 + 2) +
              *cline3 + *(cline3 + 1) + *(cline3 + 2);
            cline1 += 3;
            cline2 += 3;
            cline3 += 3;
            *bdata++ = sum / 9;
          }
        } else {
          sdata = ((b3dInt16 *)brray) + bxOffset + iy * nxBdim;
          for (ix = 0; ix < nxout; ix++) {
            *sdata++ = (b3dInt16)(*cline1) + *(cline1 + 1) + *(cline1 + 2) +
              *cline2 + *(cline2 + 1) + *(cline2 + 2) +
              *cline3 + *(cline3 + 1) + *(cline3 + 2);
            cline1 += 3;
            cline2 += 3;
            cline3 += 3;
          }
        }
      }
      break;
      
    case SLICE_MODE_SHORT:
      for (iy = 0; iy < nyout; iy++) {
        sline1 = ((b3dInt16 *)array) + (3 * iy + iyofs) * nxDim + ixofs;
        sline2 = sline1 + nxDim;
        sline3 = sline2 + nxDim;
        sdata = ((b3dInt16 *)brray) + bxOffset + iy * nxBdim;
        for (ix = 0; ix < nxout; ix++) {
          sum = 4 + *sline1 + *(sline1 + 1) + *(sline1 + 2) +
            *sline2 + *(sline2 + 1) + *(sline2 + 2) +
            *sline3 + *(sline3 + 1) + *(sline3 + 2);
          sline1 += 3;
          sline2 += 3;
          sline3 += 3;
          *sdata++ = sum / 9;
        }
      }
      break;
        
    case SLICE_MODE_USHORT:
      for (iy = 0; iy < nyout; iy++) {
        usline1 = ((b3dUInt16 *)array) + (3 * iy + iyofs) * nxDim + ixofs;
        usline2 = usline1 + nxDim;
        usline3 = usline2 + nxDim;
        usdata = ((b3dUInt16 *)brray) + bxOffset + iy * nxBdim;
        for (ix = 0; ix < nxout; ix++) {
          sum = 4 + *usline1 + *(usline1 + 1) + *(usline1 + 2) +
            *usline2 + *(usline2 + 1) + *(usline2 + 2) +
            *usline3 + *(usline3 + 1) + *(usline3 + 2);
          usline1 += 3;
          usline2 += 3;
          usline3 += 3;
          *usdata++ = sum / 9;
        }
      }
      break;
        
    case SLICE_MODE_FLOAT:
      for (iy = 0; iy < nyout; iy++) {
        fline1 = ((b3dFloat *)array) + (3 * iy + iyofs) * nxDim + ixofs;
        fline2 = fline1 + nxDim;
        fline3 = fline2 + nxDim;
        fdata = ((float *)brray) + bxOffset + iy * nxBdim;
        for (ix = 0; ix < nxout; ix++) {
          fsum = *fline1 + *(fline1 + 1) + *(fline1 + 2) +
            *fline2 + *(fline2 + 1) + *(fline2 + 2) +
            *fline3 + *(fline3 + 1) + *(fline3 + 2);
          fline1 += 3;
          fline2 += 3;
          fline3 += 3;
          *fdata++ = fsum / 9.f;
        }
      }
      break;
    }

  } else if (nbin == 4 && bin234OK) {

    /* Binning by 4 */
    switch (type) {
    case SLICE_MODE_BYTE:
      for (iy = 0; iy < nyout; iy++) {
        cline1 = ((unsigned char *)array) + (4 * iy + iyofs) * nxDim + ixofs;
        cline2 = cline1 + nxDim;
        cline3 = cline2 + nxDim;
        cline4 = cline3 + nxDim;
        if (keepByte) {
          bdata = ((unsigned char *)brray) + bxOffset + iy * nxBdim;
          for (ix = 0; ix < nxout; ix++) {
            sum = 8 + *cline1 + *(cline1 + 1) + *(cline1 + 2) + *(cline1 + 3) +
              *cline2 + *(cline2 + 1) + *(cline2 + 2) + *(cline2 + 3) +
              *cline3 + *(cline3 + 1) + *(cline3 + 2) + *(cline3 + 3) +
              *cline4 + *(cline4 + 1) + *(cline4 + 2) + *(cline4 + 3);
            cline1 += 4;
            cline2 += 4;
            cline3 += 4;
            cline4 += 4;
            *bdata++ = sum / 16;
          }
        } else {
          sdata = ((b3dInt16 *)brray) + bxOffset + iy * nxBdim;
          for (ix = 0; ix < nxout; ix++) {
            *sdata++ = (b3dInt16)(*cline1) + *(cline1 + 1) + *(cline1 + 2) + 
              *(cline1 + 3) +
              *cline2 + *(cline2 + 1) + *(cline2 + 2) + *(cline2 + 3) +
              *cline3 + *(cline3 + 1) + *(cline3 + 2) + *(cline3 + 3) +
              *cline4 + *(cline4 + 1) + *(cline4 + 2) + *(cline4 + 3);
            cline1 += 4;
            cline2 += 4;
            cline3 += 4;
            cline4 += 4;
          }
        }
      }
      break;
        
    case SLICE_MODE_SHORT:
      for (iy = 0; iy < nyout; iy++) {
        sline1 = ((b3dInt16 *)array) + (4 * iy + iyofs) * nxDim + ixofs;
        sline2 = sline1 + nxDim;
        sline3 = sline2 + nxDim;
        sline4 = sline3 + nxDim;
        sdata = ((b3dInt16 *)brray) + bxOffset + iy * nxBdim;
        for (ix = 0; ix < nxout; ix++) {
          sum = 8 +  *sline1 + *(sline1 + 1) + *(sline1 + 2) + *(sline1 + 3) +
            *sline2 + *(sline2 + 1) + *(sline2 + 2) + *(sline2 + 3) +
            *sline3 + *(sline3 + 1) + *(sline3 + 2) + *(sline3 + 3) +
            *sline4 + *(sline4 + 1) + *(sline4 + 2) + *(sline4 + 3);
          sline1 += 4;
          sline2 += 4;
          sline3 += 4;
          sline4 += 4;
          *sdata++ = sum / 16;
        }
      }
      break;
        
    case SLICE_MODE_USHORT:
      for (iy = 0; iy < nyout; iy++) {
        usline1 = ((b3dUInt16 *)array) + (4 * iy + iyofs) * nxDim + ixofs;
        usline2 = usline1 + nxDim;
        usline3 = usline2 + nxDim;
        usline4 = usline3 + nxDim;
        usdata = ((b3dUInt16 *)brray) + bxOffset + iy * nxBdim;
        for (ix = 0; ix < nxout; ix++) {
          sum = 8 +  *usline1 + *(usline1 + 1) + *(usline1 + 2) + 
            *(usline1 + 3) +
            *usline2 + *(usline2 + 1) + *(usline2 + 2) + *(usline2 + 3) +
            *usline3 + *(usline3 + 1) + *(usline3 + 2) + *(usline3 + 3) +
            *usline4 + *(usline4 + 1) + *(usline4 + 2) + *(usline4 + 3);
          usline1 += 4;
          usline2 += 4;
          usline3 += 4;
          usline4 += 4;
          *usdata++ = sum / 16;
        }
      }
      break;
        
    case SLICE_MODE_FLOAT:
      for (iy = 0; iy < nyout; iy++) {
        fline1 = ((b3dFloat *)array) + (4 * iy + iyofs) * nxDim + ixofs;
        fline2 = fline1 + nxDim;
        fline3 = fline2 + nxDim;
        fline4 = fline3 + nxDim;
        fdata = ((float *)brray) + bxOffset + iy * nxBdim;
        for (ix = 0; ix < nxout; ix++) {
          fsum = *fline1 + *(fline1 + 1) + *(fline1 + 2) + *(fline1 + 3) +
            *fline2 + *(fline2 + 1) + *(fline2 + 2) + *(fline2 + 3) +
            *fline3 + *(fline3 + 1) + *(fline3 + 2) + *(fline3 + 3) +
            *fline4 + *(fline4 + 1) + *(fline4 + 2) + *(fline4 + 3);
          fline1 += 4;
          fline2 += 4;
          fline3 += 4;
          fline4 += 4;
          *fdata++ = fsum / 16.f;
        }
      }
      break;
    }

  } else {

    /* Bin by arbitrary number or bin RGB */
    switch (type) {
    case SLICE_MODE_BYTE:
      for (iy = 0; iy < nyout; iy++) {
        cline1 = ((unsigned char *)array) + (nbin * iy + iyofs) * nxDim + ixofs;
        bdata = ((unsigned char *)brray) + bxOffset + iy * nxBdim;
        sdata = ((b3dInt16 *)brray) + bxOffset + iy * nxBdim;
        for (ix = 0; ix < nxout; ix++) {
          sum = 0;
          cline2 = cline1;
          for (j = 0; j < nbin; j++) {
            for (i = 0; i < nbin; i++) 
              sum += cline2[i];
            cline2 += nxDim;
          }
          if (keepByte > 0)
            *bdata++ = (unsigned char)((sum + nbinsq / 2) / nbinsq);
          else if (!keepByte)
            *sdata++ = sum;
          else
            *sdata++ = (sum + nbinsq / 2) / nbinsq;
          cline1 += nbin;
        }
      }
      break;

    case SLICE_MODE_SHORT:
      for (iy = 0; iy < nyout; iy++) {
        sline1 = ((b3dInt16 *)array) + (nbin * iy + iyofs) * nxDim+ ixofs;
        sdata = ((b3dInt16 *)brray) + bxOffset + iy * nxBdim;
        for (ix = 0; ix < nxout; ix++) {
          sum = nbinsq / 2;
          sline2 = sline1;
          for (j = 0; j < nbin; j++) {
            for (i = 0; i < nbin; i++) 
              sum += sline2[i];
            sline2 += nxDim;
          }
          *sdata++ = sum / nbinsq;
          sline1 += nbin;
        }
      }
      break;

    case SLICE_MODE_USHORT:
      for (iy = 0; iy < nyout; iy++) {
        usline1 = ((b3dUInt16 *)array) + (nbin * iy + iyofs) * nxDim+ ixofs;
        usdata = ((b3dUInt16 *)brray) + bxOffset + iy * nxBdim;
        for (ix = 0; ix < nxout; ix++) {
          sum = nbinsq / 2;
          usline2 = usline1;
          for (j = 0; j < nbin; j++) {
            for (i = 0; i < nbin; i++) 
              sum += usline2[i];
            usline2 += nxDim;
          }
          *usdata++ = sum / nbinsq;
          usline1 += nbin;
        }
      }
      break;

    case SLICE_MODE_FLOAT:
      for (iy = 0; iy < nyout; iy++) {
        fline1 = ((b3dFloat *)array) + (nbin * iy + iyofs) * nxDim + ixofs;
        fdata = ((float *)brray) + bxOffset + iy * nxBdim;
        for (ix = 0; ix < nxout; ix++) {
          fsum = 0.;
          fline2 = fline1;
          for (j = 0; j < nbin; j++) {
            for (i = 0; i < nbin; i++) 
              fsum += fline2[i];
            fline2 += nxDim;
          }
          *fdata++ = fsum / (b3dFloat)nbinsq;
          fline1 += nbin;
        }
      }
      break;

    case SLICE_MODE_RGB:
      for (iy = 0; iy < nyout; iy++) {
        cline1 = ((unsigned char *)array) + 3 * ((nbin * iy + iyofs) * nxDim + 
                                                 ixofs);
        bdata = ((unsigned char *)brray) + 3 * (bxOffset + iy * nxBdim);
        sdata = ((b3dInt16 *)brray) + bxOffset + iy * nxBdim;
        for (ix = 0; ix < nxout; ix++) {
          red = green = blue = 0;
          cline2 = cline1;
          for (j = 0; j < nbin; j++) {
            for (i = 0; i < nbin; i++) {
              red += cline2[3*i];
              green += cline2[3*i + 1];
              blue += cline2[3*i + 2];
            }
            cline2 += nxDim * 3;
          }
          if (keepByte > 0) {
            *bdata++ = (unsigned char)((red + nbinsq / 2) / nbinsq);
            *bdata++ = (unsigned char)((green + nbinsq / 2) / nbinsq);
            *bdata++ = (unsigned char)((blue + nbinsq / 2) / nbinsq);
          } else if (!keepByte)
            *sdata++ = red + green + blue;
          else
            *sdata++ = (red + green + blue + nbinsq / 2) / nbinsq;
          cline1 += nbin * 3;
        }
      }
      break;

    }
  }
  *nxr = nxout;
  *nyr = nyout;
  return 0;
}

/*!
 * Extracts a portion of an array into another array with optional summing of pixels 
 * (binning) by the factor [nbin].  It calls @extractAndBinIntoArray with output X 
 * dimension and offsets of 0 to produce a contiguously packed output in [brray].  Other 
 * arguments are as described there.  [brray] can be the  same as [array] unless the 
 * binning is 1, the input is bytes and [keepByte] is not greater than 0.  Returns 1 for 
 * an unsupported data type or 2 for an input X range bigger than [nxDim].
 */
int extractWithBinning(void *array, int type, int nxDim, int xStart, int xEnd, int yStart,
                       int yEnd, int nbin, void *brray, int keepByte, int *nxr, int *nyr)
{
  return extractAndBinIntoArray(array, type, nxDim, xStart, xEnd, yStart, yEnd, nbin,
                                brray, 0, 0, 0, keepByte, nxr, nyr);
}

/*!
 * Fortran wrapper for @extractWithBinning with floating point data
 */
int extractwithbinning(float *array, int *nxDim, int *xStart, int *xEnd,
                       int *yStart, int *yEnd, int *nbin, float *brray, int *keepByte,
                       int *nxr, int *nyr)
{
  return extractWithBinning(array, SLICE_MODE_FLOAT, *nxDim, *xStart, *xEnd, *yStart,
                            *yEnd, *nbin, brray, *keepByte, nxr, nyr);
}


/*!
 * Fortran wrapper for call to @extractWithBinning to simply repack a portion of a
 * 2D array sequentially into a 1-D array.
 */
void irepak(void *brray, void *array, int *nxin, int *nyin, int *xStart, int *xEnd,
            int *yStart, int *yEnd)
{
  int nxr, nyr;
  extractWithBinning(array, SLICE_MODE_FLOAT, *nxin, *xStart, *xEnd, *yStart,
                     *yEnd, 1, brray, 1, &nxr, &nyr);
}

/*!
 * Reduces a full array in size by summing pixels (binning) by the factor [nbin].
 * [array] has the input, with dimensions [nxin] by [nyin].  It simply calls
 * @extractWithBinning with ranges of 0 to [nxin] - 1 and 0 to [nyin] - 1;
 * other arguments have the same meaning as there.  [brray] can be the 
 * same as [array], and [nxr], [nyr] can be the same variables as [nxin],
 * [nyin].  Returns 1 for an unsupported data type.
 */
int reduceByBinning(void *array, int type, int nxin, int nyin, int nbin, 
                    void *brray, int keepByte, int *nxr, int *nyr)
{
  return extractWithBinning(array, type, nxin, 0, nxin - 1, 0, nyin - 1, nbin, brray,
                            keepByte, nxr, nyr);
}

/*!
 * Fortran wrapper for @reduceByBinning with floating point data, called as 
 * {reduce_by_binning}.  
 * Again, the output variables can safely be the same as the input variables.
 */
void reduce_by_binning(float *array, int *nx, int *ny, int *nbin,
                       float *brray, int *nxr, int *nyr)
{
  reduceByBinning(array, SLICE_MODE_FLOAT, *nx, *ny, *nbin, brray, 0, nxr, 
                  nyr);
}

/*!
 * Bins a slice of data in [array], with X dimension [nxDim], by binning factors of 
 * [binFacX] and [binFacY] in X and Y, and adds it into the slice in [brray] with a 
 * weighting of [zWeight]. The number of binned pixels to produce in X and Y is given by 
 * [nxBin] and [nyBin], and [brray] is contiguous (has X dimension [nxBin]).
 */
void binIntoSlice(float *array, int nxDim, float *brray, int nxBin, int nyBin,
                  int binFacX, int binFacY, float zWeight)
{
  int ixBin, iyBin, ix, iy, ind, ixBase, indBase, numThreads;
  float factor = zWeight / (binFacX * binFacY);

  if (binFacX * binFacY == 1) {
    
    /* Parallelizing this did more harm than good */
    for (iyBin = 0; iyBin < nyBin; iyBin ++) {
      for (ixBin = 0; ixBin < nxBin; ixBin ++) {
        ind = ixBin + nxDim * iyBin;
        brray[ind] += array[ind] * factor;
      }
    }
  } else {

    /* TODO: find good limit for this */
    numThreads = numOMPthreads(8);
#pragma omp parallel for num_threads(numThreads)                      \
  shared(nxDim, nxBin, nyBin, factor, array, brray, binFacX, binFacY) \
  private(ixBin, iyBin, ix, iy, ind, ixBase, indBase)
    for (iyBin = 0; iyBin < nyBin; iyBin ++) {
      for (iy = iyBin * binFacY; iy < (iyBin + 1) * binFacY; iy++) {
        ixBase = nxDim * iy;
        indBase = nxBin * iyBin;
        if (binFacX == 1) {
          for (ixBin = 0; ixBin < nxBin; ixBin ++)
            brray[ixBin + indBase] += array[ixBin + ixBase] * factor;

        } else if (binFacX == 2) {
          for (ixBin = 0; ixBin < nxBin; ixBin ++) {
            ind = 2 * ixBin + ixBase;
            brray[ixBin + indBase] += (array[ind] + array[ind + 1]) * factor;
          }
        } else if (binFacX == 3) {
          for (ixBin = 0; ixBin < nxBin; ixBin ++) {
            ind = 3 * ixBin + ixBase;
            brray[ixBin + indBase] += (array[ind] + array[ind + 1] + array[ind + 2]) *
              factor;
          }
        } else if (binFacX == 4) {
          for (ixBin = 0; ixBin < nxBin; ixBin ++) {
            ind = 4 * ixBin + ixBase;
            brray[ixBin + indBase] += (array[ind] + array[ind + 1] + array[ind + 2] +
                                       array[ind + 3]) * factor;
          }
        } else {
          for (ixBin = 0; ixBin < nxBin; ixBin ++) {
            ind = ixBin + indBase;
            for (ix = ixBin * binFacX; ix < (ixBin + 1) * binFacX; ix++)
              brray[ind] += array[ix + ixBase] * factor;
          }
        }
      }
    }
  }
}

/*! Fortran wrapper for @binIntoSlice */
void binintoslice(float *array, int *nxDim, float *brray, int *nxBin, int *nyBin,
                  int *binFacX, int *binFacY, float *zWeight)
{
  binIntoSlice(array, *nxDim, brray, *nxBin, *nyBin, *binFacX, *binFacY, *zWeight);
}
