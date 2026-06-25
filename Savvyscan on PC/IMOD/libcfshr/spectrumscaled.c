/*  spectrumscaled.c - Compute a scaled and/or reduced power spectrum
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 2016 by the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 * $Id$
 */

#include "b3dutil.h"

static double FFTMagnitude(float *array, int nx, int ny, int ix, int iy);
/*!
 * Computes a symmetrical power spectrum of the image in [image], with slice type in 
 * [type] and size in [nx] and [ny].  The image will padded to the square size in 
 * [padSize], which should be a good size for taking a FFT, and tapered down to its mean 
 * inside the edges.  If the negative of [padSize] is passed, there will be no tapering.
 * The FFT will be taken, the logarithm of the magnitudes will be
 * computed and scaled to occupy a range of 0 to 32000, and this spectrum will be 
 * mirrored through the origin.  If [finalSize] is less than [padSize], the spectrum will 
 * be reduced with antialiased reduction, using the filter type in [filtType]. If
 * [bkgdGray] is less than or equal to 0, the result will be placed in [spectrum] as
 * an array of shorts.  If [bkgdGray] is greater than 0, the spectrum
 * will be scaled to range from 0 to 255, with the mean near the edges scaled to
 * a value of [bkgdGray], and the mean of pixels in a ring with diameter [truncDiam] as a
 * fraction [finalSize] scaled to 255.  This result will be returned in [spectrum] as an 
 * array of bytes.  [twoDfft] must be supplied with the name of a 2-D FFT function with
 * that calling sequence, (todfft if using libifft in IMOD).  ^
 * The behavior is quite different if [filtType] is negative: 1) Reduction, if any, will
 * be done by Fourier cropping with @@filtxcorr.html#fourierReduceImage@; 2) No logarithm
 * is taken and amplitudes are not scaled; 3) [spectrum] is returned as floats with the
 * X dimension still padded by 2 from the FFTs, and must be allocated accordingly. ^
 * Returns -1 or -2 for the negative of an error in calling 
 * @@interpbin.html#selectZoomFilter@, 1 to 5 for an error in calling 
 * @@interpbin.html#zoomWithFilter@, -3 for a memory error, or -4 for [finalSize] greater 
 * than [padSize], [bkgdGray] more than 192, [truncDiam] negative or more than 0.75, or
 * [filtType] less than 0 and [bkgdGray] positive.
 */
int spectrumScaled(void *image, int type, int nx, int ny, void *spectrum, int padSize, 
                   int finalSize, int bkgdGray, float truncDiam, int filtType,
                   void (*twoDfft)(float *, int *, int *, int *))
{
  float fracTaper = 0.05f;
  double logScale = 5.;
  int dataType = type;
  int nxTaper, nyTaper;
  double val, cenMax;  
  double scale, sum, ringRad, rad;
  float minScale, maxScale, bkgd, shift;
  int zero = 0, one = 1;
  int i, ixbase, iyin, iyout, nsum, ix, iy, minRing, maxRing, reducing, retVal = 0;
  short int *sdata, *sdata2;
  int scaling = bkgdGray > 0 ? 1 : 0;
  int cenLim, cropping, padXdim;
  double zoom = (double)finalSize / padSize;
  float *fftArray = NULL;
  short *sTemp = (short *)spectrum;
  short *redOut = (short *)spectrum;
  float *fTemp = (float *)spectrum;
  float *fData, fData2;
  short *scaleIn;
  unsigned char *bdata = (unsigned char *)spectrum;
  unsigned char **linePtrs = NULL;

  if (padSize < 0) {
    padSize = -padSize;
    fracTaper = 0.;
  }
  padXdim = padSize + 2;
  reducing = (finalSize < padSize && filtType >= 0) ? 1 : 0;
  cropping = (finalSize < padSize && filtType < 0) ? 1 : 0;
  if (reducing) {
    retVal = -selectZoomFilter(filtType, zoom, &iyout);
    if (retVal)
      return retVal;
  }
  if (finalSize > padSize || bkgdGray > 192 || truncDiam < 0. || truncDiam > 0.75 ||
      (scaling && filtType < 0))
    return -4;

  /* Get the arrays */
  fftArray = B3DMALLOC(float, padXdim * padSize);
  if (!fftArray)
    return -3;
  if (scaling || reducing) {
    sTemp = B3DMALLOC(short, padXdim * padSize);
    if (reducing)
      linePtrs = makeLinePointers(sTemp, padSize, padSize, 2);
    if ((reducing && !linePtrs) || !sTemp) {
      B3DFREE(fftArray);
      B3DFREE(linePtrs);
      B3DFREE(sTemp);      
      return -3;
    }
    scaleIn = sTemp;
  }
  if (cropping) {
    fTemp = B3DMALLOC(float, padXdim * padSize);
    if (!fTemp) {
      B3DFREE(fftArray);
      return -3;
    }
  }
  if (scaling && reducing) {
    redOut = (short *)fftArray;
    scaleIn = redOut;
  }

  /* Get the image into the tapered, padded real array */
  nxTaper = nx * fracTaper;
  nyTaper = ny * fracTaper;
  sliceTaperInPad(image, dataType, nx, 0, nx - 1, 0, ny - 1, fftArray, padSize + 2,
                  padSize, padSize, nxTaper, nyTaper);

  /* Take the FFT */
  twoDfft(fftArray, &padSize, &padSize, &zero);

  /* Determine magnitudes around center, ignore 0,0 */
  cenMax = 0.;
  cenLim = B3DMAX(5, padSize / 10);
  if (cenLim >= padSize / 2)
    cenLim = padSize / 2 -1;
  for (i = 0; i < cenLim; i++) {
    for (iyin = 0; iyin < cenLim; iyin++) {
      if (!i && !iyin)
        continue;

      val = FFTMagnitude(fftArray, padSize, padSize, i, iyin);
      ACCUM_MAX(cenMax, val);
      val = FFTMagnitude(fftArray, padSize, padSize, i, padSize - 1 - iyin);
      ACCUM_MAX(cenMax, val);
    }
  }

  /* Process unscaled float option into padded temp and/or output arrays */
  if (filtType < 0) {
    makeAmplitudeSpectrum(fftArray, fTemp, padSize, padXdim);
    fTemp[padXdim * padSize / 2 + padSize / 2] = 0.;

    /* Fourier crop to final size */
    if (cropping) {
      twoDfft(fTemp, &padSize, &padSize, &zero);
      shift = 0.5 * ((float)padSize / (float)finalSize - 1.);
      fourierReduceImage(fTemp, padSize, padSize, (float *)redOut, finalSize, finalSize, 
                         shift, shift, fftArray);
      twoDfft((float *)redOut, &finalSize, &finalSize, &one);
      free(fTemp);
    }
    free(fftArray);
    return 0;
  }

  /* Determine inside log scaling from mean of border */
  sum = 0.;
  for (iyin = 0; iyin < padSize; iyin++)
    sum += FFTMagnitude(fftArray, padSize, padSize, padSize / 2, iyin);
  logScale /= sum / padSize;

  scale = 32000. / log(logScale * cenMax + 1.);

  /* Fill lower right of output array from top of FFT then wrap into bottom
     of FFT to fill upper right */
  iyin = padSize / 2;
  for (iyout = 0; iyout < padSize; iyout++) {
    ixbase = iyin * padXdim;
    sdata = sTemp + iyout * padSize + padSize / 2;
    for (i = ixbase; i < ixbase + padSize; i += 2) {
      val = sqrt((double)(fftArray[i] * fftArray[i] + 
                          fftArray[i + 1] * fftArray[i + 1]));
      *sdata++ = (int)(scale * log(logScale * val + 1.));
    }

    /* Mirror the next value (which could be the extra pixel) on the left edge */
    val = sqrt((double)(fftArray[i] * fftArray[i] + 
                        fftArray[i + 1] * fftArray[i + 1]));
    sTemp[((padSize - iyout) % padSize) * padSize] = 
      (int)(scale * log(logScale * val + 1.));

    /* Advance iyin and wrap when needed */
    iyin++;
    iyin %= padSize;
  }

  /* Mirror the lines in Y around the center */
  for (iyout = 0; iyout < padSize; iyout++) {
    iyin = iyout ? padSize - iyout : 0;
    sdata = sTemp + iyin * padSize + padSize / 2 + 1;
    sdata2 = sTemp + iyout * padSize + padSize / 2 - 1;
    for (i = 0; i < padSize / 2 - 1; i++)
      *sdata2-- = *sdata++;
  }
  sTemp[padSize * padSize / 2 + padSize / 2] = 32000;

  /* Reduce if called for */
  retVal = 0;
  if (reducing)
    retVal = zoomWithFilter(linePtrs, padSize, padSize, 0., 0., finalSize, finalSize, 
                            finalSize, 0, SLICE_MODE_SHORT, redOut, NULL, NULL);

  /* Scale if no error and gray value set */
  if (!retVal && scaling) {
    minScale = 0.;
    maxScale = 32000.;
    if (finalSize > 50) {

      /* Get background mean from left edge */
      sum = 0;
      for (iy = 0; iy < finalSize; iy++)
        for (ix = 2; ix < 5; ix++)
          sum += scaleIn[ix + iy * finalSize];
      bkgd = sum / (3 * finalSize);
      
      /* Average a ring at the given diameter relative to image size */
      ringRad = sqrt((double)finalSize * finalSize) * truncDiam / 2.;
      maxRing = (int)ceil(ringRad) + 1;
      minRing = (int)(0.7 * ringRad) - 1;
      sum = 0.;
      nsum = 0;
      for (iy = -maxRing; iy <= maxRing; iy++) {
        for (ix = -maxRing; ix <= maxRing; ix++) {
          if (iy < -minRing || iy > minRing || ix < -minRing || ix > minRing) {
            rad = sqrt((double)ix * ix + iy * iy);
            if (fabs(rad - ringRad) < 0.71) {
              sum +=  scaleIn[finalSize / 2 + ix + (finalSize / 2 + iy) * finalSize];
              nsum++;
            }
          }
        }
      }

      /* Set the max from the ring average and the min to place the background at the
         given gray level */
      if (nsum > 3) {
        maxScale = sum / (float)nsum;
        val = (float)bkgdGray / 256.f;
        minScale = (bkgd - maxScale * val) / (1.f - val);
      }
    }

    /* Scale the data */
    scale = 255. / (maxScale - minScale);
    for (i = 0; i < finalSize * finalSize; i++) {
      iyout = (int)(scale * (scaleIn[i] - minScale));
      B3DCLAMP(iyout, 0, 255);
      bdata[i] = iyout;
    }
  }

  if (scaling || reducing)
    free(sTemp);
  free(fftArray);
  B3DFREE(linePtrs);
  return retVal;
}

static double FFTMagnitude(float *array, int nx, int ny, int ix, int iy)
{
  int i = ix * 2 + iy * (nx + 2);
  return sqrt((double)(array[i] * array[i] + array[i + 1] * array[i + 1]));
}

/*!
 * Makes a centered amplitude spectrum from an FFT of a square image in [fftArray] and 
 * puts it in [spectrum].  The absolute value of [padSize] is the real-space size of the 
 * image and [outXdim] is the X dimension for the output array.  Set [padSize] positive if
 * an FFT is passed in [fftArray] (with an assumed X dimension of [padSize] + 2), or
 * negative if [fftArray] already contains amplitudes (with an assumed X dimension of
 *([padSize] + 2) / 2).
 */
void makeAmplitudeSpectrum(float *fftArray, float *spectrum, int padSize, int outXdim)
{
  float *fdata, *fdata2;
  int iyin, iyout, ixbase, i, inXdim = padSize + 2, amplitudes = 0;
  if (padSize < 0) {
    amplitudes = 1;
    padSize = -padSize;
    inXdim = (padSize + 2) / 2;
  }

  /* Fill lower right of output array from top of FFT then wrap into bottom
     of FFT to fill upper right */
  iyin = padSize / 2;
  for (iyout = 0; iyout < padSize; iyout++) {
    ixbase = iyin * inXdim;
    fdata = spectrum + iyout * outXdim + padSize / 2;
    if (amplitudes) {
      for (i = ixbase; i < ixbase + padSize / 2; i++)
        *fdata++ = fftArray[i];
    
      /* Mirror the next value (which could be the extra pixel) on the left edge */
      spectrum[((padSize - iyout) % padSize) * outXdim] = fftArray[i];

    } else {
      for (i = ixbase; i < ixbase + padSize; i += 2)
        *fdata++ = sqrt((double)(fftArray[i] * fftArray[i] + 
                                 fftArray[i + 1] * fftArray[i + 1]));
    
      /* Mirror the next value (which could be the extra pixel) on the left edge */
      spectrum[((padSize - iyout) % padSize) * outXdim] = 
        sqrt((double)(fftArray[i] * fftArray[i] + fftArray[i + 1] * fftArray[i + 1]));
    }
    
    /* Advance iyin and wrap when needed */
    iyin++;
    iyin %= padSize;
  }
  
  /* Mirror the lines in Y around the center */
  for (iyout = 0; iyout < padSize; iyout++) {
    iyin = iyout ? padSize - iyout : 0;
    fdata = spectrum + iyin * outXdim + padSize / 2 + 1;
    fdata2 = spectrum + iyout * outXdim + padSize / 2 - 1;
    for (i = 0; i < padSize / 2 - 1; i++)
      *fdata2-- = *fdata++;
  }
}
