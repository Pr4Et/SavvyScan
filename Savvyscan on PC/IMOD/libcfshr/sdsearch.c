/*  sdsearch.c - functions for measuring SD of difference between images
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 2014 by the Regents of the University of Colorado.  
 *  See dist/COPYRIGHT for full copyright notice.
 *
 * $Id$
 */

#include "imodconfig.h"
#include "b3dutil.h"

#ifdef F77FUNCAP
#define montbigsearch MONTBIGSEARCH
#define montsdcalc MONTSDCALC
#else
#define montbigsearch montbigsearch_
#define montsdcalc montsdcalc_
#endif

/*!
 * Compares images in two arrays [array] and [brray] (dimensioned
 * [nx], [ny], with image dimensions [nx] by [ny]), and finds the displacements between 
 * the two images with the minimum standard deviation of the pixel-by-pixel difference.
 * [ixBox0], [iyBox0], [ixBox1], and [iyBox1] (numbered from 0) define the array index 
 * limits of the box to be compared in ARRAY.  [dxMin] and [dyMin], upon entering,
 * should have the starting displacement between the box position in
 * [array] and the position to be compared in [brray].  Initially it will
 * search with a step size of 1, up to a distance of [limStep] from the
 * initial displacements.  It will do [numIter] iterations, cutting the step
 * size by a factor of 2 on each iteration.  [sxMin] and [dyMin] will be returned
 * with the displacement that has the minimum standard deviation [sdMin],
 * and [ddenMin] is the difference in density at that displacement (B-A). 
 */
void montBigSearch(float *array, float *brray, int nx, int ny, int ixBox0,
                   int iyBox0,int ixBox1, int iyBox1, float *dxMin, float *dyMin,
                   float *sdMin, float *ddenMin, int numIter, int limStep)
{
#define CHKLIM  100
  unsigned char checked[2 * CHKLIM + 1][2 * CHKLIM + 1];
  int xChanged, yChanged, keepon;
  int idxMin, idyMin, newIdx, newIdy, ndxy, iter, i, j, ndLim, idir;
  int absNewIdx, absNewIdy, ixDir, iyDir; 
  float dxcen, dycen, dxy, sd, dden;
  
  dxcen = *dxMin;        /* all moves relative to center */
  dycen = *dyMin;        /* at the initial position */
  montSdCalc(array, brray, nx, ny, ixBox0, iyBox0, ixBox1, iyBox1, *dxMin, *dyMin,
              sdMin, ddenMin);
  ndxy = 1;          /* initial # of steps to move by */
  for (iter = 0; iter < numIter - 1; iter++)
    ndxy *= 2;
  dxy = 1. / ndxy;         /* final true step size */
  ndLim = ndxy * limStep;      /* limiting # of steps */
  if (ndLim > CHKLIM)
    ndLim = CHKLIM;
  
  
  /*  set grid to unchecked state except at center */
  for (i = -ndLim; i <= ndLim; i++) {
    for (j = -ndLim ; j <= ndLim; j++)
      checked[j + CHKLIM][i + CHKLIM] = 0;
  }
  checked[CHKLIM][CHKLIM] = 1;
  
  idxMin = 0;
  idyMin = 0;
  for (iter = 1; iter <= numIter; iter++) {
    xChanged = 1;
    yChanged = 1;

    /* keep doing the x-y-diagonal series as long as something changes */
    while (xChanged || yChanged) {

      /*     move in x until reach minimum */
      /* first try a positive step */
      xChanged = 0;
      idir = 1;
      newIdx = idxMin + idir * ndxy;
      absNewIdx = newIdx < 0 ? -newIdx : newIdx;
      if (absNewIdx <= ndLim && !checked[newIdx + CHKLIM][idyMin + CHKLIM]) {
        montSdCalc(array, brray, nx, ny, ixBox0, iyBox0, ixBox1, iyBox1, 
                    dxcen + newIdx * dxy, dycen + idyMin * dxy, &sd, &dden);

        /* mark as checked*/
        checked[newIdx + CHKLIM][idyMin + CHKLIM] = 1; 
        if(sd < *sdMin) {
          *sdMin = sd;
          *ddenMin = dden;
          xChanged = 1;
          idxMin = newIdx;
        } else {             /* if positive step does no good */
          idir = -1;           /* switch direction */
        }
      } else {         /* or if positive has been */
        idir = -1;             /* checked or is too far, switch */
      }

      /*   now, in whichever direction is selected, keep on moving 
      until get to a higher value or edge or area */
      keepon = 1;
      while(keepon) {
        newIdx = idxMin+idir*ndxy;
        absNewIdx = newIdx < 0 ? -newIdx : newIdx;
        if(absNewIdx <= ndLim && 
          !checked[newIdx + CHKLIM][idyMin + CHKLIM]){
          montSdCalc(array, brray, nx, ny, ixBox0, iyBox0, ixBox1, iyBox1,
            dxcen + newIdx * dxy, dycen + idyMin * dxy, &sd, &dden);
          checked[newIdx + CHKLIM][idyMin + CHKLIM] = 1; 
          
          if (sd < *sdMin) {
            *sdMin = sd;
            *ddenMin = dden;
            xChanged = 1;
            idxMin = newIdx;
          } else {  /* if value is higher, stop going */
            keepon = 0;
          }
        } else {     /* or if already checked, or too */
          keepon = 0;    /* far, stop going */
        }
      }
      
      /*  follow exact same procedure in y direction */
      yChanged = 0;
      idir = 1;
      newIdy = idyMin + idir * ndxy;
      absNewIdy = newIdy < 0 ? -newIdy : newIdy;
      if(absNewIdy <= ndLim && 
        !checked[idxMin + CHKLIM][newIdy + CHKLIM]) {
        montSdCalc(array, brray, nx, ny, ixBox0, iyBox0, ixBox1, iyBox1,
                    dxcen + idxMin * dxy, dycen + newIdy * dxy, &sd, &dden);
        checked[idxMin + CHKLIM][newIdy + CHKLIM] = 1;
        if(sd < *sdMin) {
          *sdMin = sd;
          *ddenMin = dden;
          yChanged = 1;
          idyMin = newIdy;
        } else {
          idir = -1;
        }
      } else {
        idir = -1;
      }
      
      keepon = 1;
      while(keepon) {
        newIdy = idyMin+idir*ndxy;
        absNewIdy = newIdy < 0 ? -newIdy : newIdy;
        if(absNewIdy <= ndLim && 
          !checked[idxMin + CHKLIM][newIdy + CHKLIM]) {
          montSdCalc(array, brray, nx, ny, ixBox0, iyBox0, ixBox1, iyBox1,
                      dxcen + idxMin * dxy, dycen + newIdy * dxy, &sd, &dden);
          checked[idxMin + CHKLIM][newIdy + CHKLIM] = 1;
          if(sd < *sdMin) {
            *sdMin = sd;
            *ddenMin = dden;
            yChanged = 1;
            idyMin = newIdy;
          } else {
            keepon = 0;
          }
        } else {
          keepon = 0;
        }
      }
      
      /*  if no change above, check the 4 corners */
      for (ixDir = -1; ixDir <= 1; ixDir += 2) {
        for (iyDir = -1; iyDir <= 1; iyDir += 2) {

          /* keep checking corners only if nothing changed yet*/
          if(!(xChanged || yChanged)) {
            newIdx = idxMin+ixDir*ndxy;
            newIdy = idyMin+iyDir*ndxy;
            absNewIdy = newIdy < 0 ? -newIdy : newIdy;
            absNewIdx = newIdx < 0 ? -newIdx : newIdx;
            if(absNewIdy <= ndLim && absNewIdx <= ndLim &&
              !checked[newIdx + CHKLIM][newIdy + CHKLIM]) {
              montSdCalc(array, brray, nx, ny, ixBox0, iyBox0, ixBox1, iyBox1,
                          dxcen + newIdx * dxy, dycen + newIdy * dxy, &sd, &dden);
              checked[newIdx + CHKLIM][newIdy + CHKLIM] = 1;
              if(sd < *sdMin) {
                *sdMin = sd;
                *ddenMin = dden;
                xChanged = 1;
                yChanged = 1;
                idxMin = newIdx;
                idyMin = newIdy;
              }
            }
          }
        }
      }
      
      /*  if nothing ever changed in last loop, drop out of loop */
    }
    ndxy = ndxy / 2;       /*  and cut the step size */
  }
  *dxMin = dxcen + idxMin * dxy;
  *dyMin = dycen + idyMin * dxy;
}

/*!
 * Fortran wrapper to @montBigSearch, where the array index limits are numbered from 1.
 */
void montbigsearch(float *array, float *brray, int *nx, int *ny, int *ixBox0,
                   int *iyBox0,int *ixBox1, int *iyBox1, float *dxMin, float *dyMin,
                   float *sdMin, float *ddenMin, int *numIter, int *limStep)
{
  montBigSearch(array, brray, *nx, *ny, *ixBox0 - 1, *iyBox0 - 1, *ixBox1 - 1, 
                *iyBox1 - 1, dxMin, dyMin, sdMin, ddenMin, *numIter, *limStep);
}


/*!
 * Compares positions in two image arrays [array] and [brray] (each
 * dimensioned to [nx], [ny], the image dimensions as well), and within a
 * defined box, calculates the standard deviation [sd] of the pixel-by-pixel difference 
 * (B minus A).  [ixBox0], [iyBox0], [ixBox1], and [iyBox1] (numbered from 0) define the
 * array index limits of the box in [array], and [dx] and [dy] define the
 * displacement between that position and the desired position in [brray].
 * If [dx] and [dy] are integers, the computation is done rapidly between
 * corresponding pixels; otherwise, quadratic interpolation is used to
 * find the pixel value in [brray] for a given pixel in [array]. 
 */
void montSdCalc(float *array, float *brray, int nx, int ny, int ixBox0, int iyBox0,
                int ixBox1, int iyBox1, float dx, float dy, float *sd, float *dden)
{
  int nsum, ixp, iyp, idx, idy, ixb0, iyb0, ixb1, iyb1, iya, iypp1, iypm1;
  int iyb, ixb;
  float xp, yp, dxint, dyint, diff, c2, c8, c6, c4, c5, dxsq, dysq, binterp;
  double sum, sumsq;
  
  nsum = 0;
  sum = 0.;
  sumsq = 0.;
  *sd = 9999.;

  /* convert dx and dy into integer idx,idy and fractions dxint, dyint
     use high end of box to be sure to get right dxint and dyint values */
  yp = iyBox1 + dy;
  iyp = B3DNINT(yp);
  dyint = yp - iyp;
  idy = dy - dyint;
  xp = ixBox1 + dx;
  ixp = B3DNINT(xp);
  dxint = xp - ixp;
  idx = dx - dxint;
  if (dxint == 0. && dyint == 0.) {

    /* If dx, dy are integers, don't need interpolation */
    /* compute limits for box in brray */
    ixb0 = B3DMAX(0, ixBox0 + idx);
    ixb1 = B3DMIN(nx - 1, ixBox1 + idx);
    iyb0 = B3DMAX(0, iyBox0 + idy);
    iyb1 = B3DMIN(ny - 1, iyBox1 + idy);
    
    if (ixb0 <= ixb1 && ixb0 < nx && ixb1 >= 0 && 
        iyb0 <= iyb1 && iyb0 < ny && iyb1 >= 0) {
      for (iyb = iyb0; iyb <= iyb1; iyb++) {
        iya = iyb - idy;
        for (ixb = ixb0; ixb <= ixb1; ixb++) {
          diff = brray[ixb + iyb * nx] - 
            array[ixb - idx + iya * nx];
          sum = sum + diff;
          sumsq = sumsq + diff * diff;
        }
      }
      nsum = (ixb1 + 1 - ixb0) * (iyb1 + 1 - iyb0);
    }
  } else {
    
    /* compute the coefficients for quadratic interpolation */
    dysq = dyint * dyint;
    c8 = 0.5 * (dysq + dyint);
    c2 = 0.5 * (dysq - dyint);
    dxsq = dxint * dxint;
    c6 = 0.5 * (dxsq + dxint);
    c4 = 0.5 * (dxsq - dxint);
    c5 = 1. - dxsq - dysq;
    
    ixb0 = B3DMAX(1, ixBox0 + idx);
    ixb1 = B3DMIN(nx - 2, ixBox1 + idx);
    iyb0 = B3DMAX(1, iyBox0 + idy);
    iyb1 = B3DMIN(ny - 2, iyBox1 + idy);
    if(ixb0 <= ixb1 && ixb0 < nx - 1 && ixb1 > 0
      && iyb0 <= iyb1 && iyb0 < ny - 1 && iyb1 > 0) {
      for (iyp = iyb0; iyp <= iyb1; iyp++) {
        iya = iyp - idy;
        iypp1 = iyp + 1;
        iypm1 = iyp - 1;
        for (ixp = ixb0; ixp <= ixb1; ixp++) {
          binterp = c5 * brray[ixp + iyp * nx]
            + c4 * brray[ixp - 1 + iyp * nx] 
            + c6 * brray[ixp + 1 + iyp * nx]
            + c2 * brray[ixp + iypm1 * nx] 
            + c8 * brray[ixp + iypp1 * nx];
          diff = binterp - array[ixp - idx + iya * nx];
          sum = sum + diff;
          sumsq = sumsq + diff * diff;
        }
      }
      nsum = (ixb1 + 1 - ixb0) * (iyb1 + 1 - iyb0);
    }
  }
  if (nsum > 0)
    *dden = sum / nsum;
  if (nsum > 1)
    *sd = sqrt((sumsq - sum * sum / nsum) / (nsum - 1.));
}

/*!
 * Fortran wrapper to @montSdCalc, where the array index limits are numbered from 1.
 */
void montsdcalc(float *array, float *brray, int *nx, int *ny, int *ixBox0, int *iyBox0,
                int *ixBox1, int *iyBox1, float *dx, float *dy, float *sd, float *dden)
{
  montSdCalc(array, brray, *nx, *ny, *ixBox0 - 1, *iyBox0 - 1, *ixBox1 - 1, *iyBox1 - 1,
             *dx, *dy, sd, dden);
}
