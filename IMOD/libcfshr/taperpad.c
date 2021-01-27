/*
 * taperpd.c - Taper and pad a slice inside or out
 *
 * Copyright (C) 2006 by Boulder Laboratory for 3-Dimensional Electron
 * Microscopy of Cells ("BL3DEMC") and the Regents of the University of
 * Colorado.  See dist/COPYRIGHT for full notice.
 *
 * $Id$
 */

#include "mrcslice.h"
#include "imodconfig.h"
#include "b3dutil.h"

#ifdef F77FUNCAP
#define taperoutpad TAPEROUTPAD
#define smoothoutpad SMOOTHOUTPAD
#define taperinpad TAPERINPAD
#define taperinpadex TAPERINPADEX
#define slicenoisetaperpad SLICENOISETAPERPAD
#define sliceedgemean SLICEEDGEMEAN
#define splitfill SPLITFILL
#else
#define taperoutpad taperoutpad_
#define smoothoutpad smoothoutpad_
#define taperinpad taperinpad_
#define taperinpadex taperinpadex_
#define slicenoisetaperpad slicenoisetaperpad_
#define sliceedgemean sliceedgemean_
#define splitfill splitfill_
#endif

static void copyToCenter(void *array, int type, int nxbox, int nybox,
                         float *brray, int nxdim, int nx, int ny, int *ixlo,
                         int *ixhi, int *iylo, int *iyhi);
static void getRunningMeanSD(float *boxStart, int noiseRows, int noiseLength,
                             int rowStride, int lenStride, int nbox, float dmean,
                             float *samples, float *sampleSq, int numSample, float *means,
                             float *SDs);

/*!
 * Pads an image in [array] with dimensions [nxbox] by [nybox] into the center
 * of a larger array, [brray], which can be the same as [array].  The
 * SLICE_MODE is specified by [type], which must be byte, float or short
 * integer. The size of the padded image is specified by [nx] and [ny] while
 * [nxdim] specifies the X dimension of the output array.  The values of the
 * image in the padding area will be tapered from the value of a pixel at the
 * edge of the input image edge down to a common value.  That value is either
 * the mean at the edge of the input image, if [ifmean] is zero, or the value
 * supplied in [dmeanin], if [ifmean] is nonzero.  The offset to the original data in the
 * padded array is ([nx] - [nxbox]) / 2, ([ny] - [nybox]) / 2.
 */
void sliceTaperOutPad(void *array, int type, int nxbox, int nybox,
                      float *brray, int nxdim, int nx, int ny, int ifmean,
                      float dmeanin)
{
  int ixlo, ixhi, iylo, iyhi, ix, iy, nxtop, nytop, iytoplin;
  float dmean, edgel, edger, wedge, wmean, prodmean;
#ifdef PAD_TIMIMG
  static double wallCen = 0., wallSide = 0., wallTB = 0., wallEdge = 0.;
  double thisCen, now, wallStart = wallTime();
#endif

  copyToCenter(array, type, nxbox, nybox, brray, nxdim, nx, ny, &ixlo, &ixhi,
               &iylo, &iyhi);
#ifdef PAD_TIMIMG
  now = wallTime();
  thisCen = now - wallStart;
  wallCen += now - wallStart;
  wallStart = now;
#endif

  /* Do the taper if there is any padding */
  if (nxbox != nx || nybox != ny) {
    if (ifmean)
      dmean = dmeanin;
    else
      dmean = (float)sliceEdgeMean(brray, nxdim, ixlo, ixhi - 1, iylo, iyhi - 1);
    nxtop = nx - 1;
    nytop = ny - 1;

    /*  if there is a mismatch between left and right, add a column on
        right; similarly for bottom versus top, add a row on top */
    if (nx - ixhi > ixlo) {
      nxtop--;
      for (iy = 0; iy < ny; iy++)
        brray[nx - 1 + iy * nxdim] = dmean;
    }
    if (ny - iyhi > iylo) {
      nytop--;
      for (ix = 0; ix < nx; ix++)
        brray[ix + (ny - 1) * nxdim] = dmean;
    }
#ifdef PAD_TIMIMG
    now = wallTime();
    wallEdge += now - wallStart;
    wallStart = now;
#endif

    for (iy = iylo; iy < iyhi; iy++) {
      edgel = brray[ixlo + iy * nxdim];
      edger = brray[ixhi - 1 + iy * nxdim];
      for (ix = 0; ix < ixlo; ix++) {
        wedge = ((float)ix) / ixlo;
        wmean = 1.f - wedge;
        brray[ix + iy * nxdim] = wmean * dmean + wedge * edgel;
        brray[nxtop - ix + iy * nxdim] = wmean * dmean + wedge * edger;
      }
    }
#ifdef PAD_TIMIMG
    now = wallTime();
    wallSide += now - wallStart;
    wallStart = now;
#endif

    for (iy = 0; iy < iylo; iy++) {
      wedge = ((float)iy) / iylo;
      prodmean = (1.f - wedge) * dmean;
      iytoplin = nytop - iy;
      for (ix = 0; ix < nx; ix++) {
        brray[ix + iy * nxdim] = prodmean + wedge * brray[ix + iylo * nxdim];
        brray[ix + iytoplin * nxdim] = prodmean +
          wedge * brray[ix + (iyhi - 1) * nxdim];
      }
    }
#ifdef PAD_TIMIMG
    now = wallTime();
    wallTB += now - wallStart;
    printf("Cen + %.1f = %.1f edge %.1f  side %.1f  T/B %.1f\n", thisCen *1000.,
           wallCen *1000., wallEdge *1000., wallSide * 1000., wallTB *1000.);
#endif
  }
}

/* Common function for copying the box to the center of the array and returning
   limiting coordinates */
static void copyToCenter(void *array, int type, int nxbox, int nybox,
                         float *brray, int nxdim, int nx, int ny, int *ixlo,
                         int *ixhi, int *iylo, int *iyhi)
{
  b3dUByte *bytein;
  b3dInt16 *intin;
  b3dUInt16 *uintin;
  float *floatin;
  float *out;
  int  ix, iy;

  *ixlo = (nx - nxbox) / 2;
  *ixhi = *ixlo + nxbox;
  *iylo = (ny - nybox) / 2;
  *iyhi = *iylo + nybox;
  for (iy = nybox - 1; iy >= 0; iy--) {
    out = brray + *ixhi + (*iylo + iy) * nxdim - 1;
    switch (type) {
    case SLICE_MODE_BYTE:
      bytein = (b3dUByte *)array + (iy + 1) * nxbox - 1;
      for (ix = nxbox - 1; ix >= 0; ix--)
        *out-- = *bytein--;
      break;

    case SLICE_MODE_SHORT:
      intin = (b3dInt16 *)array + (iy + 1) * nxbox - 1;
      for (ix = nxbox - 1; ix >= 0; ix--)
        *out-- = *intin--;
      break;

    case SLICE_MODE_USHORT:
      uintin = (b3dUInt16 *)array + (iy + 1) * nxbox - 1;
      for (ix = nxbox - 1; ix >= 0; ix--)
        *out-- = *uintin--;
      break;

    case SLICE_MODE_FLOAT:
      floatin = (float *)array + (iy + 1) * nxbox - 1;
      for (ix = nxbox - 1; ix >= 0; ix--)
        *out-- = *floatin--;
      break;
    }
  }
}

/*!
 * Fortran wrapper to call @sliceTaperOutPad with a floating point array.
 */
void taperoutpad(void *array, int *nxbox, int *nybox, float *brray, int *nxdim,
                 int *nx, int *ny, int *ifmean, float *dmeanin)
{
  sliceTaperOutPad(array, SLICE_MODE_FLOAT, *nxbox, *nybox,
                   brray, *nxdim, *nx, *ny, *ifmean,
                   *dmeanin);
}

/*!
 * Extracts a subarea of an image, places it into the center of a potentially
 * larger array with padding, and tapers the image down to the mean value at
 * its edge, tapering pixels inside the extracted image area.  The image data
 * are in [inArray], their X dimension is [nxDimIn], and their SLICE_MODE is
 * given in [type].  The mode can be byte, signed or unsigned short, float, or
 * RGB.  The starting and ending coordinates to extract in X and Y
 * are [ixStart] to [ixEnd] and [iyStart] to [iyEnd].  The output image array is 
 * [outArray] and its X dimension is specified by [nxDimOut].  The padded image size is
 * specified by [nx] and [ny], and [nxTaper] and [nyTaper] indicate the number of
 * pixels over which to taper in X and Y.  The output array can be the same as
 * the input array if the input image fills the entire array.  The offset to the
 * original image in the padded array is [nx] / 2 - ([ixEnd] + 1 - [ixStart]) / 2,
 * [ny] / 2 - ([iyEnd] + 1 - [iyStart]) / 2.
 */
void sliceTaperInPad(void *inArray, int type, int nxDimIn, int ixStart, int ixEnd,
                     int iyStart, int iyEnd, float *outArray, int nxDimOut, int nx, int ny,
                     int nxTaper, int nyTaper)
{
  int lowBase, highBase, x1, x2, ixLow, ixHigh, iyLow, iyHigh, ix, iy, ixBase;
  int nxBox, nyBox, ixLim, y1, y2;
  float dmean, fracX, fracY, fmin;
  b3dUByte *byteIn;
  b3dInt16 *intIn;
  b3dUInt16 *uintIn;
  float *floatIn;
  float *out;
#ifdef PAD_TIMING
  double wallCen, wallSide, wallEdge, now, wallStart = wallTime();
#endif

  /* ixlo, iylo are last index below image location in output array,
     ixhi, iyhi are last index within image location */
  nxBox = ixEnd + 1 - ixStart;
  nyBox = iyEnd + 1 - iyStart;
  ixLow = nx / 2 - nxBox / 2 - 1;
  ixHigh = ixLow + nxBox;
  iyLow = ny / 2 - nyBox / 2 - 1;
  iyHigh = iyLow + nyBox;
  for (iy = iyEnd; iy >= iyStart; iy--) {
    out = outArray + ixHigh + (iyLow + 1 + iy - iyStart) * nxDimOut;
    switch (type) {
    case SLICE_MODE_BYTE:
      byteIn = (b3dUByte *)inArray + ixEnd + iy * nxDimIn;
      for (ix = ixEnd; ix >= ixStart; ix--)
        *out-- = *byteIn--;
      break;

    case SLICE_MODE_SHORT:
      intIn = (b3dInt16 *)inArray + ixEnd + iy * nxDimIn;
      for (ix = ixEnd; ix >= ixStart; ix--)
        *out-- = *intIn--;
      break;

    case SLICE_MODE_USHORT:
      uintIn = (b3dUInt16 *)inArray + ixEnd + iy * nxDimIn;
      for (ix = ixEnd; ix >= ixStart; ix--)
        *out-- = *uintIn--;
      break;

    case SLICE_MODE_FLOAT:
      floatIn = (float *)inArray + ixEnd + iy * nxDimIn;
      for (ix = ixEnd; ix >= ixStart; ix--)
        *out-- = *floatIn--;
      break;

    case SLICE_MODE_RGB:
      byteIn = (b3dUByte *)inArray + 3 * (ixEnd + iy * nxDimIn);
      for (ix = ixEnd; ix >= ixStart; ix--) {
        *out-- = byteIn[0] + byteIn[1] + byteIn[2];
        byteIn -= 3;
      }
      break;

    }
  }
#ifdef PAD_TIMING
  now = wallTime();
  wallCen = now - wallStart;
  wallStart = now;
#endif

  /* get edge mean */
  dmean = (float)sliceEdgeMean(outArray, nxDimOut, ixLow + 1, ixHigh, iyLow + 1, iyHigh);

  /* fill the rest of the array with dmean */
  if (nxBox != nx || nyBox != ny) {
    for (iy = iyLow + 1; iy <= iyHigh; iy++) {
      ixBase = iy * nxDimOut;
      for (ix = 0; ix <= ixLow; ix++)
        outArray[ix + ixBase] = dmean;
      for (ix = ixHigh + 1; ix < nx; ix++)
        outArray[ix + ixBase] = dmean;
    }
    for (iy = 0; iy <= iyLow; iy++) {
      ixBase = iy * nxDimOut;
      for (ix = 0; ix < nx; ix++)
        outArray[ix + ixBase] = dmean;
    }
    for (iy = iyHigh + 1; iy < ny; iy++) {
      ixBase = iy * nxDimOut;
      for (ix = 0; ix < nx; ix++)
        outArray[ix + ixBase] = dmean;
    }
  }
#ifdef PAD_TIMING
  now = wallTime();
  wallEdge = now - wallStart;
  wallStart = now;
#endif

  /* Taper the edges */
  for (iy = 0; iy < (nyBox + 1) / 2; iy++) {
    fracY = 1.;
    ixLim = nxTaper;
    if (iy < nyTaper) {
      fracY = (iy + 1.f) / (nyTaper + 1.f);
      ixLim = (nxBox + 1) / 2;
    }
    for (ix = 0; ix < ixLim; ix++) {
      fracX = 1.;
      if (ix < nxTaper)
        fracX = (ix + 1.f) / (nxTaper + 1.f);
      fmin = fracX < fracY ? fracX : fracY;
      if (fmin < 1.) {
        x1 = ix + 1 + ixLow;
        x2 = ixHigh - ix;
        y1 = iy + 1 + iyLow;
        y2 = iyHigh - iy;

        /*      DNM 4/28/02: for odd box sizes, deflect middle pixel to edge */
        /* to keep it from being attenuated twice */
        if (x1 == x2)
          x2 = 0;
        if (y1 == y2)
          y2 = 0;

        lowBase = y1 * nxDimOut;
        highBase = y2 * nxDimOut;
        outArray[x1 + lowBase] = fmin * (outArray[x1 + lowBase] - dmean) + dmean;
        outArray[x1 + highBase] = fmin * (outArray[x1 + highBase] - dmean) + dmean;
        outArray[x2 + lowBase] = fmin * (outArray[x2 + lowBase] - dmean) + dmean;
        outArray[x2 + highBase] = fmin * (outArray[x2 + highBase] - dmean) + dmean;
      }
    }
  }
#ifdef PAD_TIMING
  now = wallTime();
  wallSide = now - wallStart;
  printf("Cen %.1f edge %.1f  side %.1f\n", wallCen * 1000., wallEdge * 1000.,
      wallSide * 1000.);
  fflush(stdout);
#endif
}

/*!
 * Fortran wrapper to @sliceTaperInPad for padding a whole float image array
 * of size [nxbox] by [nybox].  The output array can be the same as
 * the input array.
 */
void taperinpad(void *array, int *nxbox, int *nybox, float *brray, int *nxdim,
                int *nx, int *ny, int *nxtap, int *nytap)
{
  sliceTaperInPad(array, SLICE_MODE_FLOAT, *nxbox, 0, *nxbox - 1, 0,
                  *nybox - 1, brray, *nxdim, *nx, *ny, *nxtap, *nytap);
}

/*!
 * Fortran wrapper to @sliceTaperInPad for padding a subarea from a float image array
 */
void taperinpadex(void *array, int *nxdimin, int *ix0, int *ix1,
                     int *iy0, int *iy1, float *brray, int *nxdim, int *nx, int *ny,
                     int *nxtap, int *nytap)
{
  sliceTaperInPad(array, SLICE_MODE_FLOAT, *nxdimin, *ix0, *ix1, *iy0,
                  *iy1, brray, *nxdim, *nx, *ny, *nxtap, *nytap);
}

/*!
 * Pads an image in [array] with dimensions [nxbox] by [nybox] into the center
 * of a larger array, [brray], which can be the same as [array].  The
 * SLICE_MODE is specified by [type], which must be byte, float or short
 * integer. The size of the padded image is specified by [nx] and [ny] while
 * [nxdim] specifies the X dimension of the output array.  The padding is done
 * by replicating pixels and smoothing lines progressively more for farther
 * out from the pixels in the box.  This is done in a series of progressively
 * larger rectangles; the each pixel in first rectangle contains the average
 * of the three nearest pixels on the outer edge of the actual data; each pixel
 * in the second rectangle is the average of the 5 nearest pixels in the first
 * rectangle, etc.  The offset to the original data in the
 * padded array is ([nx] - [nxbox]) / 2, ([ny] - [nybox]) / 2.
 */
void sliceSmoothOutPad(void *array, int type, int nxbox, int nybox,
                         float *brray, int nxdim, int nx, int ny)
{
  int ixlo, ixhi, iylo, iyhi, ix, iy, x, y, dirx, diry, nsum, numPad, ipad;
  int xmin, xlimlo, xmax, xlimhi, xLineLo, xLineHi, xside, xstart, xend;
  int ymin, ylimlo, ymax, ylimhi, yLineLo, yLineHi, yside, ystart, yend;
  float sum;

  copyToCenter(array, type, nxbox, nybox, brray, nxdim, nx, ny, &ixlo, &ixhi,
               &iylo, &iyhi);
  numPad = b3dIMax(4, ixlo, nx - ixhi, iylo, ny - iyhi);
  /* printf("%d %d %d %d\n", ixlo, ixhi, iylo, iyhi); */
  for (ipad = 1; ipad <= numPad; ipad++) {
    xLineLo = ixlo - ipad;
    xLineHi = ixhi + ipad - 1;
    xmin = B3DMAX(0, xLineLo);
    xlimlo = B3DMAX(0, xLineLo + 1);
    xmax = B3DMIN(nx - 1, xLineHi);
    xlimhi = B3DMIN(nx - 1, xLineHi - 1);

    yLineLo = iylo - ipad;
    yLineHi = iyhi + ipad - 1;
    ymin = B3DMAX(0, yLineLo);
    ylimlo = B3DMAX(0, yLineLo + 1);
    ymax = B3DMIN(ny - 1, yLineHi);
    ylimhi = B3DMIN(ny - 1, yLineHi - 1);
    /* printf("ipad %d X: %d %d %d %d %d %d\n", ipad, xLineLo, xLineHi, xmin,
       xmax, xlimlo, xlimhi);
       printf("ipad %d Y: %d %d %d %d %d %d\n", ipad, yLineLo, yLineHi, ymin,
       ymax, ylimlo, ylimhi); */

    diry = 1;
    for (iy = yLineLo; iy <= yLineHi; iy += yLineHi - yLineLo) {
      if (iy < 0 || iy >= ny) {
        diry = -1;
        continue;
      }
      for (ix = xmin; ix <= xmax; ix++) {
        xside = -1;
        xstart = ix - ipad;
        if (xstart < xlimlo) {
          xstart = xlimlo;
          if (xstart) {
            ystart = iy + 2 * diry;
            yend = ystart + diry * (xstart - (ix - ipad) - 1);
            yend = B3DMAX(ylimlo + 1, B3DMIN(ylimhi - 1, yend));
            if (diry * ystart <= diry * yend)
              xside = xstart;
          }
        }
        xend = ix + ipad;
        if (xend > xlimhi) {
          xend = xlimhi;
          if (xend < nx && xside < 0) {
            ystart = iy + 2 * diry;
            yend = ystart + diry * (ix + ipad - xend - 1);
            yend = B3DMAX(ylimlo + 1, B3DMIN(ylimhi - 1, yend));
            if (diry * ystart <= diry * yend)
              xside = xend;
          }
        }

        sum = 0.;
        nsum = xend + 1 - xstart;
        /* printf("ix %d iy %d xs %d xe %d xside %d diry %d\n", ix, iy,
           xstart, xend, xside, diry); */
        for (x = xstart; x <= xend; x++)
          sum += brray[x + nxdim * (iy + diry)];
        if (xside >= 0) {
          /* printf("ys %d  ye %d  nsumst %d sumst %.1f", ystart, yend, sum,
             nsum); */
          for (y = ystart; diry * y <= diry * yend; y += diry) {
            sum += brray[xside + nxdim * y];
            nsum++;
          }
          /* printf(" %.1f / %d = %.1f\n", sum, nsum, sum / nsum); */
        }
        brray[ix + nxdim * iy] = sum / nsum;
      }
      diry = -1;
    }

    dirx = 1;
    for (ix = xLineLo; ix <= xLineHi; ix += xLineHi - xLineLo) {
      if (ix < 0 || ix >= nx) {
        dirx = -1;
        continue;
      }
      for (iy = ymin; iy <= ymax; iy++) {
        yside = -1;
        ystart = iy - ipad;
        if (ystart < ylimlo) {
          ystart = ylimlo;
          if (ystart) {
            xstart = ix + 2 * dirx;
            xend = xstart + dirx * (ystart - (iy - ipad) - 1);
            xend = B3DMAX(xlimlo + 1, B3DMIN(xlimhi - 1, xend));
            if (dirx * xstart <= dirx * xend)
              yside = ystart;
          }
        }
        yend = iy + ipad;
        if (yend > ylimhi) {
          yend = ylimhi;
          if (yend < ny && yside < 0) {
            xstart = ix + 2 * dirx;
            xend = xstart + dirx * (iy + ipad - yend - 1);
            xend = B3DMAX(xlimlo + 1, B3DMIN(xlimhi - 1, xend));
            if (dirx * xstart <= dirx * xend)
              yside = yend;
          }
        }

        sum = 0.;
        nsum = yend + 1 - ystart;
        /* printf("ix %d iy %d ys %d ye %d yside %d dirx %d\n", ix, iy,
           ystart, yend, yside, dirx); */
        for (y = ystart; y <= yend; y++)
          sum += brray[ix + dirx + nxdim * y];
        if (yside >= 0) {
          /* printf("xs %d  xe %d  nsumst %d sumst %.1f", xstart, xend, sum,
             nsum); */
          for (x = xstart; dirx * x <= dirx * xend; x += dirx) {
            sum += brray[x + nxdim * yside];
            nsum++;
          }
          /* printf(" %.1f / %d = %.1f\n", sum, nsum, sum / nsum); */
        }
        brray[ix + nxdim * iy] = sum / nsum;
      }
      dirx = -1;
    }
  }
}

/*!
 * Fortran wrapper to call @sliceSmoothOutPad with a floating point array.
 */
void smoothoutpad(void *array, int *nxbox, int *nybox, float *brray,
                  int *nxdim, int *nx, int *ny)
{
  sliceSmoothOutPad(array, SLICE_MODE_FLOAT, *nxbox, *nybox,
                    brray, *nxdim, *nx, *ny);
}


/* Macro for sliceNoiseTaperPad to fill one corner of the padded area */
#define FILL_CORNER(xind, yind)                                         \
  pseudo = pseudoVals[0];                                               \
  for (iy = 0; iy < iylo; iy++) {                                       \
    for (ix = 0; ix < ixlo; ix++) {                                     \
      fracMin = B3DMIN(fracx[ix], fracy[iy]);                           \
      wgtSum = B3DMAX(0.01, fracx[ix] + fracy[iy]);                     \
      cornMean = (cornXmean * fracx[ix] + cornYmean * fracy[iy]) / wgtSum; \
      cornSD = (cornXsd * fracx[ix] + cornYsd * fracy[iy]) / wgtSum;    \
      pseudo = (197 * (pseudo + 1)) & 0xFFFFF;                          \
      brray[(xind) + (yind) * nxdim] = (cornMean + ((float)pseudo * ranFac - maxSDs) * \
                                        cornSD) * fracMin + dmean * (1. - fracMin); \
    }                                                                   \
  }                                                                     \
  pseudoVals[0] = pseudo;

#define MAX_NOISE_ROWS 5
#define MAX_NOISE_LENGTH 120
#define MAX_SAMPLES (MAX_NOISE_ROWS * MAX_NOISE_LENGTH)
#define SNTP_MAX_THREADS 16

/*!
 * Pads an image in [array] with dimensions [nxbox] by [nybox] into the center
 * of a larger array, [brray], which can be the same as [array].  The
 * SLICE_MODE is specified by [type], which must be byte, float or short
 * integer. The size of the padded image is specified by [nx] and [ny] while
 * [nxdim] specifies the X dimension of the output array.  The values of the
 * image in the padding area will be drawn randomly from an elongated box adjacent to
 * the image edge and attenuated down to the mean of the input image edge at the
 * edge of the output image.  The width of the box is supplied in [noiseRows]; the length
 * of the box is supplied in [noiseLength] and is limited to 600 divided by the width.
 * The box is centered on the pixel at the edge nearest to a pixel that is being filled.
 * An array must supplied in [temp] that is at least 2 times the maximum of [nxBox] and
 * [nyBox] plus [nx] - [nxbox] plus [ny] - [nybox].  The offset to the original data in
 * the padded array is ([nx] - [nxbox]) / 2, ([ny] - [nybox]) / 2.
 */
void sliceNoiseTaperPad(void *array, int type, int nxbox, int nybox, float *brray,
                         int nxdim, int nx, int ny, int noiseLength, int noiseRows,
                         float *temp)
{
  int ixlo, ixhi, iylo, iyhi, ix, iy, nxtop, nytop, numSamples;
  int ixBase, pseudo;
  float dmean, ranFac, fracMin, cornMean, cornSD, maxSDs = 1.73f;
  float cornXmean, cornYmean, cornXsd, cornYsd, saveCornMean, saveCornSD, wgtSum;
  float samples[MAX_SAMPLES];
  float sampleSq[MAX_SAMPLES];
  float *fracx, *fracy, *meanXpart, *meanYpart, *means, *SDs;
  static int pseudoVals[SNTP_MAX_THREADS] = {123456, 654321, 368341, 789234, 234561, 543216,
                                             683413, 892347, 345612, 432165, 834136, 923478,
                                             456123, 321654, 341368, 234789};
  int thread, maxThreads = 4, numThreads, iyStart[SNTP_MAX_THREADS], iyEnd[SNTP_MAX_THREADS];
#ifdef PAD_TIMIMG
  static double wallCen = 0., wallSide = 0., wallTB = 0., wallEdge = 0.;
  double now, wallStart = wallTime();
#endif

  /* Do not enforce those defined limits, take the given number of rows and limit the
     length of the sample by the product */
  noiseRows = b3dIMin(3, noiseRows, nxbox / 2, nybox / 2);
  noiseLength = 2 * (b3dIMin(4, noiseLength, MAX_SAMPLES / noiseRows, nxbox, nybox) / 2);
  numSamples = noiseRows * noiseLength;

  copyToCenter(array, type, nxbox, nybox, brray, nxdim, nx, ny, &ixlo, &ixhi,
               &iylo, &iyhi);
#ifdef PAD_TIMIMG
  now = wallTime();
  wallCen += now - wallStart;
  wallStart = now;
#endif

  /* Do the taper if there is any padding */
  if (nxbox == nx && nybox == ny)
    return;

  /* Need a mean from which to compute deviations */
  dmean = (float)sliceEdgeMean(brray, nxdim, ixlo, ixhi - 1, iylo, iyhi - 1);
  nxtop = nx - 1;
  nytop = ny - 1;

  /* Set up fractions and mean components */
  fracx = temp;
  meanXpart = temp + ixlo;
  fracy = meanXpart + ixlo;
  meanYpart = fracy + iylo;
  means = meanYpart + iylo;
  SDs = means + B3DMAX(nxbox, nybox);
  for (iy = 0; iy < iylo; iy++) {
    fracy[iy] = (float)iy / iylo;
    meanYpart[iy] = (1. - fracy[iy]) * dmean;
  }
  for (ix = 0; ix < ixlo; ix++) {
    fracx[ix] = (float)ix / ixlo;
    meanXpart[ix] = (1. - fracx[ix]) * dmean;
  }

  /*  if there is a mismatch between left and right, add a column on
      right; similarly for bottom versus top, add a row on top */
  if (nx - ixhi > ixlo) {
    nxtop--;
    for (iy = 0; iy < ny; iy++)
      brray[nx - 1 + iy * nxdim] = dmean;
  }
  if (ny - iyhi > iylo) {
    nytop--;
    for (ix = 0; ix < nx; ix++)
      brray[ix + (ny - 1) * nxdim] = dmean;
  }

  /* Multiply 20-bit random numbers by this factor to get a range of 2 * maxSDs */
  ranFac = maxSDs / 0x7FFFF;

  /* Set up the number of threads and divide the first Y range (iylo) into groups */
  maxThreads = B3DMIN(iylo, maxThreads);
  numThreads = B3DNINT(sqrt((ixlo + iylo) * (nx + ny)) / 170.);
  B3DCLAMP(numThreads, 1, maxThreads);
  numThreads = numOMPthreads(maxThreads);
  numThreads = B3DMIN(numThreads, SNTP_MAX_THREADS);
  numThreads = B3DMIN(iylo, numThreads);
  for (thread = 0; thread < numThreads; thread++)
    balancedGroupLimits(iylo, numThreads, thread, &iyStart[thread], &iyEnd[thread]);

  /* Get mean/SD and fill bottom  */
  getRunningMeanSD(&brray[ixlo + iylo * nxdim], noiseRows, noiseLength, nxdim, 1,
                   nxbox, dmean,  samples, sampleSq, numSamples, means, SDs);
#ifdef PAD_TIMIMG
  now = wallTime();
  wallEdge += now - wallStart;
  wallStart = now;
#endif

#pragma omp parallel for num_threads(numThreads)                                 \
  shared(brray, nxbox, means, ixlo, iylo, nytop, maxSDs, ranFac, SDs, fracy, meanYpart) \
  shared(nxdim, pseudoVals, iyStart, iyEnd, nxtop, fracx, meanXpart)       \
  private(iy, pseudo, ixBase, ix, thread)
  for (thread = 0; thread < numThreads; thread++) {
    pseudo = pseudoVals[thread];
    for (iy = iyStart[thread]; iy <= iyEnd[thread]; iy++) {
      ixBase = ixlo + iy * nxdim;
      for (ix = 0; ix < nxbox; ix++) {

        /* This is a linear (mixed?) congruential generator with a period of 2^20
           its deficiencies (small period of low order bits) are of no concern here */
        pseudo = (197 * (pseudo + 1)) & 0xFFFFF;
        brray[ix + ixBase] = (means[ix] + ((float)pseudo * ranFac - maxSDs) * SDs[ix]) *
          fracy[iy] + meanYpart[iy];
      }
    }
    pseudoVals[thread] = pseudo;
  }

#ifdef PAD_TIMIMG
  now = wallTime();
  wallTB += now - wallStart;
  wallStart = now;
#endif

  /* Save left-hand corner mean for the end, and right mean before getting new means */
  saveCornMean = means[0];
  saveCornSD = SDs[0];
  cornXmean = means[nxbox - 1];
  cornXsd = SDs[nxbox - 1];

  /* Right side mean/SD, bottom right corner then right side */
  getRunningMeanSD(&brray[ixlo + nxbox - noiseRows + iylo * nxdim], noiseRows,
                   noiseLength, 1, nxdim, nybox, dmean, samples, sampleSq, numSamples,
                   means, SDs);

  cornYmean = means[0];
  cornYsd = SDs[0];
  FILL_CORNER(nxtop - ix, iy);

#ifdef PAD_TIMIMG
  now = wallTime();
  wallEdge += now - wallStart;
  wallStart = now;
#endif

  for (thread = 0; thread < numThreads; thread++)
    balancedGroupLimits(nybox, numThreads, thread, &iyStart[thread], &iyEnd[thread]);

#pragma omp parallel for num_threads(numThreads)                                 \
  shared(brray, nxbox, means, ixlo, iylo, nytop, maxSDs, ranFac, SDs, fracy, meanYpart) \
  shared(nxdim, pseudoVals, iyStart, iyEnd, nxtop, fracx, meanXpart)       \
  private(iy, pseudo, ixBase, ix, thread)
  for (thread = 0; thread < numThreads; thread++) {
    pseudo = pseudoVals[thread];
    for (iy = iyStart[thread]; iy <= iyEnd[thread]; iy++) {
      ixBase = nxtop + (iy + iylo) * nxdim;
      for (ix = 0; ix < ixlo; ix++) {
        pseudo = (197 * (pseudo + 1)) & 0xFFFFF;
        brray[ixBase - ix] = (means[iy] + ((float)pseudo * ranFac - maxSDs) * SDs[iy]) *
          fracx[ix] + meanXpart[ix];
      }
    }
    pseudoVals[thread] = pseudo;
  }

#ifdef PAD_TIMIMG
  now = wallTime();
  wallSide += now - wallStart;
  wallStart = now;
#endif

  cornYmean = means[nybox - 1];
  cornYsd = SDs[nybox - 1];

  /* Top mean/SD, top right corner and top side */
  getRunningMeanSD(&brray[ixlo + (iylo + nybox - noiseRows) * nxdim], noiseRows,
                   noiseLength, nxdim, 1, nxbox, dmean, samples, sampleSq, numSamples,
                   means, SDs);
  cornXmean = means[nxbox - 1];
  cornXsd = SDs[nxbox - 1];
  FILL_CORNER(nxtop - ix, nytop - iy);

#ifdef PAD_TIMIMG
  now = wallTime();
  wallEdge += now - wallStart;
  wallStart = now;
#endif
  for (thread = 0; thread < numThreads; thread++)
    balancedGroupLimits(iylo, numThreads, thread, &iyStart[thread], &iyEnd[thread]);

#pragma omp parallel for num_threads(numThreads)                                 \
  shared(brray, nxbox, means, ixlo, iylo, nytop, maxSDs, ranFac, SDs, fracy, meanYpart) \
  shared(nxdim, pseudoVals, iyStart, iyEnd, nxtop, fracx, meanXpart)       \
  private(iy, pseudo, ixBase, ix, thread)
  for (thread = 0; thread < numThreads; thread++) {
    pseudo = pseudoVals[thread];
    for (iy = iyStart[thread]; iy <= iyEnd[thread]; iy++) {
      ixBase = ixlo + (nytop - iy) * nxdim;
      for (ix = 0; ix < nxbox; ix++) {
        pseudo = (197 * (pseudo + 1)) & 0xFFFFF;
        brray[ix + ixBase] = (means[ix] + ((float)pseudo * ranFac - maxSDs) * SDs[ix]) *
          fracy[iy] + meanYpart[iy];
      }
    }
    pseudoVals[thread] = pseudo;
  }

#ifdef PAD_TIMIMG
  now = wallTime();
  wallTB += now - wallStart;
  wallStart = now;
#endif
  cornXmean = means[0];
  cornXsd = SDs[0];

  /* Left side mean, then top left corner, then left side */
  getRunningMeanSD(&brray[ixlo + iylo * nxdim], noiseRows,
                   noiseLength, 1, nxdim, nybox, dmean, samples, sampleSq, numSamples,
                   means, SDs);
  cornYmean = means[nybox - 1];
  cornYsd = SDs[nybox - 1];
  FILL_CORNER(ix, nytop - iy);

#ifdef PAD_TIMIMG
  now = wallTime();
  wallEdge += now - wallStart;
  wallStart = now;
#endif
  for (thread = 0; thread < numThreads; thread++)
    balancedGroupLimits(nybox, numThreads, thread, &iyStart[thread], &iyEnd[thread]);

#pragma omp parallel for num_threads(numThreads)                        \
  shared(brray, nxbox, means, ixlo, iylo, nytop, maxSDs, ranFac, SDs, fracy, meanYpart) \
  shared(nxdim, pseudoVals, iyStart, iyEnd, nxtop, fracx, meanXpart)       \
  private(iy, pseudo, ixBase, ix, thread)
  for (thread = 0; thread < numThreads; thread++) {
    pseudo = pseudoVals[thread];
    for (iy = iyStart[thread]; iy <= iyEnd[thread]; iy++) {
      ixBase = (iy + iylo) * nxdim;
      for (ix = 0; ix < ixlo; ix++) {
        pseudo = (197 * (pseudo + 1)) & 0xFFFFF;
        brray[ixBase + ix] = (means[iy] + ((float)pseudo * ranFac - maxSDs) *
                              SDs[iy]) * fracx[ix] + meanXpart[ix];
      }
    }
    pseudoVals[thread] = pseudo;
  }

#ifdef PAD_TIMIMG
  now = wallTime();
  wallSide += now - wallStart;
  wallStart = now;
#endif

  /* Finish up with bottom left corner */
  cornYmean = means[0];
  cornYsd = SDs[0];
  cornXmean = saveCornMean;
  cornXsd = saveCornSD;

  FILL_CORNER(ix, iy);
#ifdef PAD_TIMIMG
  now = wallTime();
  wallEdge += now - wallStart;
  printf("Cen %.1f edge %.1f  side %.1f  T/B %.1f\n", wallCen * 1000., wallEdge * 1000.,
           wallSide * 1000., wallTB * 1000.);
#endif
}

/* Function to get mean and SD in a smaple box and compute it all along an edge */
static void getRunningMeanSD(float *boxStart, int noiseRows, int noiseLength,
                             int rowStride, int lenStride, int nbox, float dmean,
                             float *samples, float *sampleSq, int numSample, float *means,
                             float *SDs)
{
  double sampSum = 0., sampSqSum = 0.;
  float val, valSq, tmpMean, tmpSD;
  int sampleInd, len, row, meanInd;
  int meanBase = noiseLength / 2;

  /* Load the box and get first mean/SD */
  sampleInd = 0;
  for (len = 0; len < noiseLength; len++) {
    for (row = 0; row < noiseRows; row++) {
      val = boxStart[len * lenStride + row * rowStride] - dmean;
      valSq = val * val;
      samples[sampleInd] = val;
      sampleSq[sampleInd++] = valSq;
      sampSum += val;
      sampSqSum += valSq;
    }
  }
  sumsToAvgSDdbl(sampSum, sampSqSum, numSample, 1, &tmpMean, &tmpSD);
  means[meanBase] = dmean + tmpMean;
  SDs[meanBase] = tmpSD;

  /* Loop along the length, pulling samples out of sum then adding new ones in */
  sampleInd = 0;
  meanInd = meanBase + 1;
  for (; len < nbox - noiseLength; len++) {
    for (row = 0; row < noiseRows; row++) {
      sampSum -= samples[sampleInd];
      sampSqSum -= sampleSq[sampleInd];
      val = boxStart[len * lenStride + row * rowStride] - dmean;
      valSq = val * val;
      samples[sampleInd] = val;
      sampleSq[sampleInd++] = valSq;
      sampSum += val;
      sampSqSum += valSq;
    }
    sampleInd = sampleInd % numSample;
    sumsToAvgSDdbl(sampSum, sampSqSum, numSample, 1, &tmpMean, &tmpSD);
    means[meanInd] = dmean + tmpMean;
    SDs[meanInd] = tmpSD;
    meanInd++;
  }

  /* Copy the endpoints to complete the arrays */
  for (len = 0; len < meanBase; len++) {
    means[len] = means[meanBase];
    SDs[len] = SDs[meanBase];
  }
  meanBase = meanInd - 1;
  for (len = meanInd; len < nbox; len++) {
    means[len] = means[meanBase];
    SDs[len] = SDs[meanBase];
  }
}

/*! Fortran wrapper to call @sliceNoiseTaperPad with a floating point array */
void slicenoisetaperpad(float *array, int *nxbox, int *nybox, float *brray,
                        int *nxdim, int *nx, int *ny, int *noiseLength, int *noiseRows,
                        float *temp)
{
  sliceNoiseTaperPad(array, SLICE_MODE_FLOAT, *nxbox, *nybox, brray, *nxdim, *nx, *ny,
  *noiseLength, *noiseRows, temp);
}

/*!
 * Computes the mean along the edge for the portion of the image in [array]
 * bounded by X indices [ixlo] and [ixhi] and Y indices [iylo] and [iyhi],
 * inclusive.  [nxdim] is the X dimension of the array.  The mean is based
 * solely on the single columns at [ixlo] and [ixhi] and the single rows at
 * [iylo] and [iyhi].
 */
double sliceEdgeMean(float *array, int nxdim, int ixlo, int ixhi, int iylo,
                    int iyhi)
{
  double dmean, sum = 0.;
  int ix, iy;
  for (ix = ixlo; ix <= ixhi; ix++)
    sum += array[ix + iylo * nxdim] + array[ix + iyhi * nxdim];

  for (iy = iylo + 1; iy < iyhi; iy++)
    sum += array[ixlo + iy * nxdim] + array[ixhi + iy * nxdim];

  dmean = sum / (2 * (ixhi - ixlo + iyhi - iylo));
  return dmean;
}

/*!
 * Fortran wrapper to @sliceEdgeMean .  The indexes are numbered from 1.
 */
double sliceedgemean(float *array, int *nxdim, int *ixlo, int *ixhi, int *iylo,
                     int *iyhi)
{
  return sliceEdgeMean(array, *nxdim, *ixlo - 1, *ixhi - 1, *iylo - 1,
                       *iyhi - 1);
}

/*!
 * Splits a small image in [array] with dimensions [nxbox] by [nybox] into the
 * 4 corners of a potentially larger array.  The padded image is placed in
 * [brray] with a size of [nx] by [ny], where [nxdim] is the X dimension of
 * [brray].  The image will be padded with the mean at the edge of the image,
 * or with the value in [fillin] if [iffill] is non-zero.  In either case the
 * values will be shifted so that the mean of the whole array is zero.
 */
void sliceSplitFill(float *array, int nxbox, int nybox, float *brray,
               int nxdim, int nx, int ny, int iffill, float fillin)
{
  int i, ix, iy, ixlo, iylo, ixnew, iynew;
  float sum, fill, bias, dmean;
  sum = 0.;
  for (i = 0; i < nxbox * nybox; i++)
    sum += array[i];
  dmean = (float)(sum / (nxbox * nybox));
  fill = dmean;
  bias = dmean;
  if (nxbox != nx || nybox != ny) {
    fill = fillin;
    if (!iffill) {

      /* find mean of edge of box */
      sum = 0.;
      for (ix = 0; ix < nxbox; ix++)
        sum += array[ix] + array[ix + (nybox - 1) * nxbox];
      for (iy = 1; iy < nybox - 1; iy++)
        sum += array[iy * nxbox] + array[nxbox - 1 + iy * nxbox];
      fill = sum / (2 * nxbox + 2 * (nybox - 2));
    }

    /* Whatever the fill, subtract a bias that would produce a mean of zero */
    bias = fill + (dmean - fill) *((float)(nxbox * nybox)) / ((float)(nx * ny));

    /* fill whole brray with fill-bias  */
    for (iy = 0; iy < ny; iy++)
      for (ix = 0; ix < nx; ix++)
        brray[ix + iy * nxdim] = fill - bias;
  }

  /* move array into brray, splitting it into the 4 corners of brray */
  ixlo = -nxbox / 2;
  iylo = -nybox / 2;
  for (iy = 0; iy < nybox; iy++) {
    for (ix = 0; ix < nxbox; ix++) {
      ixnew = (ix + ixlo + nx) % nx;
      iynew = (iy + iylo + ny) % ny;
      brray[ixnew + iynew * nxdim] = array[ix + iy * nxbox] - bias;
    }
  }
}

/*!
 * Fortran wrapper to @sliceSplitFill .
 */
void splitfill(float *array, int *nxbox, int *nybox, float *brray,
               int *nxdim, int *nx, int *ny, int *iffill, float *fillin)
{
  sliceSplitFill(array, *nxbox, *nybox, brray, *nxdim, *nx, *ny, *iffill,
                 *fillin);
}
