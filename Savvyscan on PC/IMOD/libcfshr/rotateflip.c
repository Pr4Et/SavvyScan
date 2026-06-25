/*  rotateflip.c - Function for fast general n * 90 rotation/flip of image
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 2017 by the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 * $Id$
 */
#include "imodconfig.h"
#include "b3dutil.h"

#ifdef F77FUNCAP
#define rotateflipimage ROTATEFLIPIMAGE
#else
#define rotateflipimage rotateflipimage_
#endif

/*!
 * Performs any combination of rotation and Y flipping for the image in [array] with the
 * MRC mode given in [mode] (byte, signed/unsigned short, or float) and the size in [nx]
 * and [ny].  Set [operation] to 0-3 for counter-clockwise rotation by 90 * operation, 
 * plus 4 for flipping around the Y axis before the rotation or 8 for flipping around Y 
 * after.  Set [leftHanded] 0 for right-handed data (IMOD) or 1 for data stored inverted 
 * in Y (SerialEM).  Set [invertAfter] to include a final inversion in Y, thus making 
 * left-handed data ready to write with IMOD routines.  For short data only, [invertCon]
 * can be set non-zero to invert the contrast between the existing minimum and maximum.
 * The number of threads to use is controlled by [numThreads]; pass 0 to use the default,
 * which is 2 threads for images between 1.5 and 3K in size or for larger images when 
 * the operation does not include a rotation, or 3 threads for larger images with a 
 * rotation. The output image is returned in [brray] and its size in [nxout], [nyout].  
 * The function is relatively fast for rotations because it operates on 8 lines of input 
 * data at a time.  Returns 1 for an unsupported data mode, 2 for trying to invert 
 * contrast with non-short image, or 2 for an operation out of range.
 */
int rotateFlipImage(void *array, int mode, int nx, int ny, int operation, int leftHanded,
                    int invertAfter, int invertCon, void *brray, int *nxout, int *nyout, 
                    int numThreads)
{
  int xalong[4] = {1, 0, -1, 0};
  int yalong[4] = {0, -1, 0, 1};
  int xinter[4] = {0, 1, 0, -1};
  int yinter[4] = {1, 0, -1, 0};
  int invAfterMap[8] = {6, 5, 4, 7, 2, 1, 0, 3};
  int rightHandMap[12] = {0, 3, 2, 1, 4, 7, 6, 5, 4, 5, 6, 7};
  
  int xstart, ystart, dinter, dalong, ix, iy, strip, numStrips, rotation, flip;
  short int *bline, *blineStart, *alineStart;
  short int *aline1, *aline2, *aline3, *aline4, *aline5, *aline6, *aline7, *aline8;
  short int *bline1, *bline2, *bline3, *bline4, *bline5, *bline6, *bline7, *bline8;

  unsigned char *bubln, *bublnStart, *aublnStart;
  unsigned char *aubln1, *aubln2, *aubln3, *aubln4, *aubln5, *aubln6, *aubln7, *aubln8;
  unsigned char *bubln1, *bubln2, *bubln3, *bubln4, *bubln5, *bubln6, *bubln7, *bubln8;
  unsigned short *ualine1, *ualine2, *ualine3, *ualine4, *ualine5, *ualine6, *ualine7;
  unsigned short *ubline1, *ubline2, *ubline3, *ubline4, *ubline5, *ubline6, *ubline7;
  unsigned short *ualine8, *ubline8, *ualineStart, *ublineStart, *ubline;
  unsigned char *ubarray = (unsigned char *)array;
  unsigned char *ubbrray = (unsigned char *)brray;
  unsigned short *usarray = (unsigned short *)array;
  short *sarray = (short *)array;
  short *sbrray = (short *)brray;

  float *bfln, *bflnStart, *aflnStart;
  float *afln1, *afln2, *afln3, *afln4, *afln5, *afln6, *afln7, *afln8;
  float *bfln1, *bfln2, *bfln3, *bfln4, *bfln5, *bfln6, *bfln7, *bfln8;
  float *farray = (float *)array;
  float *fbrray = (float *)brray;
  int dmin = 1000000, dmax = -1000000;
  int nxOut;

  /* Return error for illegal mode or operation */
  if (mode != MRC_MODE_BYTE && mode != MRC_MODE_SHORT && mode != MRC_MODE_USHORT &&
      mode != MRC_MODE_FLOAT)
    return 1;
  if (invertCon && mode != MRC_MODE_SHORT && mode != MRC_MODE_USHORT)
    return 2;
  if (operation < 0 || operation > 11)
    return 3;

  /* Map a right-handed operation to the original left-handed one */
  if (!leftHanded)
    operation = rightHandMap[operation];

  /* Map the operation to produce a final flipping around X axis for output */
  if (operation < 8 && invertAfter)
    operation = invAfterMap[operation];
  rotation = operation % 4;

  /* Flip X coordinates for a flip; transpose X and Y for odd rotations */
  flip = (operation / 4) ? -1 : 1;
  *nxout = nx;
  *nyout = ny;
  if (rotation % 2) {
    *nxout = ny;
    *nyout = nx;
    
    /* For flipping before rotation, exchange the 1 and 3 rotations */
    if (flip < 0 && (operation & 4) != 0)
      rotation = 4 - rotation;
  }

  /* It is important to realize that these coordinates all start at upper left of image */
  xstart = 0;
  ystart = 0;
  if (rotation > 1)
    xstart = *nxout - 1;
  if (rotation == 1 || rotation == 2)
    ystart = *nyout - 1;
  if (flip < 0)
    xstart = *nxout - 1 - xstart;

  /* From X and Y increments along and between lines, get index increments */
  dalong = flip * xalong[rotation] + *nxout * yalong[rotation];
  dinter = flip * xinter[rotation] + *nxout * yinter[rotation];

  /* Get min/max for contrast inversion */
  if (invertCon) {
    if (mode == MRC_MODE_SHORT) {
      for (ix = 0; ix < nx * ny; ix++) {
        dmax = B3DMAX(sarray[ix], dmax);
        dmin = B3DMIN(sarray[ix], dmin);
      }
    } else {
      for (ix = 0; ix < nx * ny; ix++) {
        dmax = B3DMAX(usarray[ix], dmax);
        dmin = B3DMIN(usarray[ix], dmin);
      }
    }
  }

  /* Do the copy */
  numStrips = ny / 8;
  nxOut = *nxout;

  // Testing in ordinary compilations in Linux indicate limited efficiency especially on
  // some machines, although some did better with 4 or more threads on 8K images
  // The argument can override this, because in Windows with VC compiler, rotations are
  // slower and more threads can do much better
  if (numThreads <= 0) {
    if (nx * ny < 1500 * 1500)
      numThreads = 1;
    else if (nx * ny < 3000 * 3000 || rotation % 2 == 0)
      numThreads = 2;
    else 
      numThreads = 3;
  }
  numThreads = numOMPthreads(numThreads);
  if (mode == MRC_MODE_FLOAT) {

#pragma omp parallel for num_threads(numThreads)                           \
  shared(numStrips, nx, dinter, dalong, farray, fbrray, xstart, nxOut, ystart) \
  private(strip, bflnStart, aflnStart, ix, afln1, afln2, afln3, afln4, afln5, afln6) \
  private(afln7, afln8, bfln1, bfln2, bfln3, bfln4, bfln5, bfln6, bfln7, bfln8)
    for (strip = 0; strip < numStrips; strip++) {
      bflnStart = fbrray + xstart + nxOut * ystart + 8 * strip * dinter;
      aflnStart = farray + 8 * strip * nx;
      afln1 = aflnStart;
      afln2 = afln1 + nx;
      afln3 = afln2 + nx;
      afln4 = afln3 + nx;
      afln5 = afln4 + nx;
      afln6 = afln5 + nx;
      afln7 = afln6 + nx;
      afln8 = afln7 + nx;
      bfln1 = bflnStart;
      bfln2 = bfln1 + dinter;
      bfln3 = bfln2 + dinter;
      bfln4 = bfln3 + dinter;
      bfln5 = bfln4 + dinter;
      bfln6 = bfln5 + dinter;
      bfln7 = bfln6 + dinter;
      bfln8 = bfln7 + dinter;
      for (ix = 0; ix < nx; ix++) {
        *bfln1 = *afln1++;
        bfln1 += dalong;
        *bfln2 = *afln2++;
        bfln2 += dalong;
        *bfln3 = *afln3++;
        bfln3 += dalong;
        *bfln4 = *afln4++;
        bfln4 += dalong;
        *bfln5 = *afln5++;
        bfln5 += dalong;
        *bfln6 = *afln6++;
        bfln6 += dalong;
        *bfln7 = *afln7++;
        bfln7 += dalong;
        *bfln8 = *afln8++;
        bfln8 += dalong;
      }
    }
    bflnStart = fbrray + xstart + nxOut * ystart + 8 * numStrips * dinter;
    aflnStart = farray + 8 * numStrips *nx;
    for (iy = 8 * numStrips; iy < ny; iy++) {
      bfln = bflnStart;
      for (ix = 0; ix < nx; ix++) {
        *bfln = *aflnStart++;
        bfln += dalong;
      }
      bflnStart += dinter;
    }

  } else if (mode != MRC_MODE_BYTE) {

    if (invertCon == 0 || mode == MRC_MODE_SHORT) {

      /* SIGNED SHORTS OR NO INVERTCON */
#pragma omp parallel for num_threads(numThreads)                        \
  shared(numStrips, nx, dinter, dalong, farray, fbrray, xstart, nxOut, ystart) \
  shared(invertCon, dmin, dmax)                                         \
  private(strip, blineStart, alineStart, ix, aline1, aline2, aline3, aline4, aline5) \
  private(aline6, aline7, aline8, bline1, bline2, bline3, bline4, bline5, bline6) \
  private(bline7, bline8)
      for (strip = 0; strip < numStrips; strip++) {
        blineStart = sbrray + xstart + nxOut * ystart + 8 * strip * dinter;
        alineStart = sarray + 8 * strip * nx;
        aline1 = alineStart;
        aline2 = aline1 + nx;
        aline3 = aline2 + nx;
        aline4 = aline3 + nx;
        aline5 = aline4 + nx;
        aline6 = aline5 + nx;
        aline7 = aline6 + nx;
        aline8 = aline7 + nx;
        bline1 = blineStart;
        bline2 = bline1 + dinter;
        bline3 = bline2 + dinter;
        bline4 = bline3 + dinter;
        bline5 = bline4 + dinter;
        bline6 = bline5 + dinter;
        bline7 = bline6 + dinter;
        bline8 = bline7 + dinter;
        if (invertCon) {
          for (ix = 0; ix < nx; ix++) {
            *bline1 = (short int)(dmax + dmin - *aline1++);
            bline1 += dalong;
            *bline2 = (short int)(dmax + dmin - *aline2++);
            bline2 += dalong;
            *bline3 = (short int)(dmax + dmin - *aline3++);
            bline3 += dalong;
            *bline4 = (short int)(dmax + dmin - *aline4++);
            bline4 += dalong;
            *bline5 = (short int)(dmax + dmin - *aline5++);
            bline5 += dalong;
            *bline6 = (short int)(dmax + dmin - *aline6++);
            bline6 += dalong;
            *bline7 = (short int)(dmax + dmin - *aline7++);
            bline7 += dalong;
            *bline8 = (short int)(dmax + dmin - *aline8++);
            bline8 += dalong;
          }
        } else {
          for (ix = 0; ix < nx; ix++) {
            *bline1 = *aline1++;
            bline1 += dalong;
            *bline2 = *aline2++;
            bline2 += dalong;
            *bline3 = *aline3++;
            bline3 += dalong;
            *bline4 = *aline4++;
            bline4 += dalong;
            *bline5 = *aline5++;
            bline5 += dalong;
            *bline6 = *aline6++;
            bline6 += dalong;
            *bline7 = *aline7++;
            bline7 += dalong;
            *bline8 = *aline8++;
            bline8 += dalong;
          }
        }
      }

      /* Finish up last rows */
      blineStart = sbrray + xstart + nxOut * ystart + 8 * numStrips * dinter;
      alineStart = sarray + 8 * numStrips *nx;
      for (iy = 8 * numStrips; iy < ny; iy++) {
        bline = blineStart;
        if (invertCon) {
          for (ix = 0; ix < nx; ix++) {
            *bline = (short int)(dmax + dmin - *alineStart++);
            bline += dalong;
          }
        } else {
          for (ix = 0; ix < nx; ix++) {
            *bline = *alineStart++;
            bline += dalong;
          }
        }
        blineStart += dinter;
      }
    } else {

      /* UNSIGNED SHORTS WITH INVERTCON */
#pragma omp parallel for num_threads(numThreads)                        \
  shared(numStrips, nx, dinter, dalong, farray, fbrray, xstart, nxOut, ystart) \
  shared(invertCon, dmin, dmax)                                         \
  private(strip, ublineStart, ualineStart, ix, ualine1, ualine2, ualine3, ualine4) \
  private(ualine5, ualine6, ualine7, ualine8, ubline1, ubline2, ubline3, ubline4) \
  private(ubline5, ubline6, ubline7, ubline8)
      for (strip = 0; strip < numStrips; strip++) {
        ublineStart = (unsigned short *)sbrray + xstart + nxOut * ystart + 
          8 * strip * dinter;
        ualineStart = (unsigned short *)sarray + 8 * strip * nx;
        ualine1 = ualineStart;
        ualine2 = ualine1 + nx;
        ualine3 = ualine2 + nx;
        ualine4 = ualine3 + nx;
        ualine5 = ualine4 + nx;
        ualine6 = ualine5 + nx;
        ualine7 = ualine6 + nx;
        ualine8 = ualine7 + nx;
        ubline1 = ublineStart;
        ubline2 = ubline1 + dinter;
        ubline3 = ubline2 + dinter;
        ubline4 = ubline3 + dinter;
        ubline5 = ubline4 + dinter;
        ubline6 = ubline5 + dinter;
        ubline7 = ubline6 + dinter;
        ubline8 = ubline7 + dinter;
        for (ix = 0; ix < nx; ix++) {
          *ubline1 = (unsigned short int)(dmax + dmin - *ualine1++);
          ubline1 += dalong;
          *ubline2 = (unsigned short int)(dmax + dmin - *ualine2++);
          ubline2 += dalong;
          *ubline3 = (unsigned short int)(dmax + dmin - *ualine3++);
          ubline3 += dalong;
          *ubline4 = (unsigned short int)(dmax + dmin - *ualine4++);
          ubline4 += dalong;
          *ubline5 = (unsigned short int)(dmax + dmin - *ualine5++);
          ubline5 += dalong;
          *ubline6 = (unsigned short int)(dmax + dmin - *ualine6++);
          ubline6 += dalong;
          *ubline7 = (unsigned short int)(dmax + dmin - *ualine7++);
          ubline7 += dalong;
          *ubline8 = (unsigned short int)(dmax + dmin - *ualine8++);
          ubline8 += dalong;
        }
      }

      /* Finish up last rows */
      ublineStart = (unsigned short *)sbrray + xstart + nxOut * ystart + 
        8 * numStrips * dinter;
      ualineStart = (unsigned short *)sarray + 8 * numStrips * nx;
      for (iy = 8 * numStrips; iy < ny; iy++) {
        ubline = ublineStart;
        for (ix = 0; ix < nx; ix++) {
          *ubline = (unsigned short int)(dmax + dmin - *ualineStart++);
          ubline += dalong;
        }
        ublineStart += dinter;
      }
    }      
  } else {

    /* BYTES */
#pragma omp parallel for num_threads(numThreads)                        \
  shared(numStrips, nx, dinter, dalong, farray, fbrray, xstart, nxOut, ystart) \
  private(strip, bublnStart, aublnStart, ix, aubln1, aubln2, aubln3, aubln4, aubln5) \
  private(aubln6, aubln7, aubln8, bubln1, bubln2, bubln3, bubln4, bubln5, bubln6) \
  private(bubln7, bubln8)
    for (strip = 0; strip < numStrips; strip++) {
      bublnStart = ubbrray + xstart + nxOut * ystart + 8 * strip * dinter;
      aublnStart = ubarray + 8 * strip * nx;
      aubln1 = aublnStart;
      aubln2 = aubln1 + nx;
      aubln3 = aubln2 + nx;
      aubln4 = aubln3 + nx;
      aubln5 = aubln4 + nx;
      aubln6 = aubln5 + nx;
      aubln7 = aubln6 + nx;
      aubln8 = aubln7 + nx;
      bubln1 = bublnStart;
      bubln2 = bubln1 + dinter;
      bubln3 = bubln2 + dinter;
      bubln4 = bubln3 + dinter;
      bubln5 = bubln4 + dinter;
      bubln6 = bubln5 + dinter;
      bubln7 = bubln6 + dinter;
      bubln8 = bubln7 + dinter;
      for (ix = 0; ix < nx; ix++) {
        *bubln1 = *aubln1++;
        bubln1 += dalong;
        *bubln2 = *aubln2++;
        bubln2 += dalong;
        *bubln3 = *aubln3++;
        bubln3 += dalong;
        *bubln4 = *aubln4++;
        bubln4 += dalong;
        *bubln5 = *aubln5++;
        bubln5 += dalong;
        *bubln6 = *aubln6++;
        bubln6 += dalong;
        *bubln7 = *aubln7++;
        bubln7 += dalong;
        *bubln8 = *aubln8++;
        bubln8 += dalong;
      }
    }
    bublnStart = ubbrray + xstart + nxOut * ystart + 8 * numStrips * dinter;
    aublnStart = ubarray + 8 * numStrips *nx;
    for (iy = 8 * numStrips; iy < ny; iy++) {
      bubln = bublnStart;
      for (ix = 0; ix < nx; ix++) {
        *bubln = *aublnStart++;
        bubln += dalong;
      }
      bublnStart += dinter;
    }
  }
  return 0;
}

/*!
 * Fortran wrapper for @rotateFlipImage that assumes [array] is a floating point, 
 * right-handed image.
 */
int rotateflipimage(float *array, int *nx, int *ny, int *operation, float *brray, 
                    int *nxout, int *nyout, int *numThreads)
{
  return rotateFlipImage(array, MRC_MODE_FLOAT, *nx, *ny, *operation, 0, 0, 0, brray, 
                         nxout, nyout, *numThreads);
}
