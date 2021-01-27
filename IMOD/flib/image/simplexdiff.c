/*  simplexdiff.c - forms sums for the difference measure in xfsimplex
 *
 *  Copyright (C) 2015 by the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 *  $Id$
 */

#include <math.h>
#include "imodconfig.h"

#ifdef F77FUNCAP
#define simplexdiff SIMPLEXDIFF
#else
#define simplexdiff simplexdiff_
#endif

/* Components of the various cases */
#define NEAREST_IX_IY                           \
  fi = (float)i - xcen;                         \
  ix = amat[0] * fi + amat[2] * fj + xAdd;      \
  iy = amat[1] * fi + amat[3] * fj + yAdd;

#define LINEAR_IX_IY                            \
  fi = (float)i - xcen;                         \
  x = amat[0] * fi + amat[2] * fj + xAdd;       \
  y = amat[1] * fi + amat[3] * fj + yAdd;       \
  ix = (int)x;                                  \
  iy = (int)y;

#define NEAREST_ARRAY_DEN                       \
  den = crray[ix -1 +  (iy - 1) * nx];          \
  dend = drray[i - 1 +  (j - 1) * nx];

#define LINEAR_ARRAY_DEN                                                \
  ix1 = ix + 1;                                                         \
  iy1 = iy + 1;                                                         \
  fx = 1 + ix - x;                                                      \
  fx1 = 1. -fx;                                                         \
  fy = 1 + iy - y;                                                      \
  fy1 = 1. -fy;                                                         \
  ind = ix + (iy - 1) * nx;                                             \
  den = crray[ind - 1] * fx * fy + crray[ind] * fx1 * fy +              \
    crray[ind + nx - 1] * fx * fy1 + crray[ind + nx] * fx1 * fy1;       \
  dend = drray[i - 1 +  (j - 1) * nx];

#define PATCH_SUM                                                   \
  ixp = (i - nx1) / nxyPatch;                                       \
  ind = ixp + iyp * numXpatch;                                      \
  sxa[ind] += den - dend;                                           \
  sxSqa[ind] += (den - dend) * (den - dend);                        \
  numPixA[ind]++;

#define PATCH_CCC_SUM                                                   \
  ixp = (i - nx1) / nxyPatch;                                           \
  ind = ixp + iyp * numXpatch;                                          \
  sxa[ind] += den;                                                      \
  sya[ind] += dend;                                                     \
  sxya[ind] += den * dend;                                              \
  sxSqa[ind] += den * den;                                              \
  sySqa[ind] += dend * dend;                                            \
  numPixA[ind]++;


/* The function for forming simplex differences */
void simplexDiff(float *crray, float *drray, int nx, int ny, int nx1, int nx2, int ny1,
                 int ny2, float *amat, float xcen, float ycen, float xAdd, float yAdd, 
                 int ifInterp, int oldDiff, int ifCCC, double *sx, int *numPix,
                 int numXpatch, int nxyPatch, double *sxa, double *sya, double *sxya,
                 double *sxSqa, double *sySqa, int *numPixA)
{
  float fj, fi, x, y, fx, fx1, fy, fy1, den, dend;
  int j, i, ix, iy, ix1, iy1, ixp, iyp, ind;

  for (j = ny1; j <= ny2; j++) {

    /* Top of Y loop: get some constants for Y and compute source pixels for ends of 
       line */
    fj = (float)j - ycen;
    iyp = (j - ny1) / nxyPatch;
    fi = nx1 - xcen;
    ix1 = amat[0] * fi + amat[2] * fj + xAdd;
    iy1 = amat[1] * fi + amat[3] * fj + yAdd;
    fi = nx2 - xcen;
    ix = amat[0] * fi + amat[2] * fj + xAdd;
    iy = amat[1] * fi + amat[3] * fj + yAdd;
    if (ix > 1 && ix < nx && iy > 1 && iy < ny && ix1 > 1 && ix1 < nx && 
        iy1 > 1 && iy1 < ny) {

      /* NO RANGE TESTS */
      if (ifInterp == 0) {

        /* Nearest simple difference */
        if (oldDiff) {
          for (i = nx1; i <= nx2; i++) {
            NEAREST_IX_IY;
            NEAREST_ARRAY_DEN;
            *sx += fabs(den - dend);
          }
          *numPix += nx2 + 1 - nx1;
        } else if (!ifCCC) {

          /* Nearest patch difference */
          for (i = nx1; i <= nx2; i++) {
            NEAREST_IX_IY;
            NEAREST_ARRAY_DEN;
            PATCH_SUM;
          }
        } else {

          /* Nearest CCC */
          for (i = nx1; i <= nx2; i++) {
            NEAREST_IX_IY;
            NEAREST_ARRAY_DEN;
            PATCH_CCC_SUM;
          }
        }
      } else {
        if (oldDiff) {

          /* Linear simple difference */
          for (i = nx1; i <= nx2; i++) {
            LINEAR_IX_IY;
            LINEAR_ARRAY_DEN;
            *sx += fabs(den - dend);
          }
          *numPix += nx2 + 1 - nx1;
        } else if (!ifCCC) {

          /* Linear patch difference */
          for (i = nx1; i <= nx2; i++) {
            LINEAR_IX_IY;
            LINEAR_ARRAY_DEN;
            PATCH_SUM;
          }
        } else {

          /* Linear CCC */
          for (i = nx1; i <= nx2; i++) {
            LINEAR_IX_IY;
            LINEAR_ARRAY_DEN;
            PATCH_CCC_SUM;
          }
        }
        
      }
    } else {

      /* RANGE TESTS REQUIRED */
      if (ifInterp == 0) {

        /* Nearest simple difference */
        if (oldDiff) {
          for (i = nx1; i <= nx2; i++) {
            NEAREST_IX_IY;
            if (ix >= 1 && ix <= nx && iy >= 1 && iy <= ny) {
              NEAREST_ARRAY_DEN;
              *sx += fabs(den - dend);
              (*numPix)++;
            }
          }
        } else if (!ifCCC) {

          /* Nearest patch difference */
          for (i = nx1; i <= nx2; i++) {
            NEAREST_IX_IY;
            if (ix >= 1 && ix <= nx && iy >= 1 && iy <= ny) {
              NEAREST_ARRAY_DEN;
              PATCH_SUM;
            }
          }
        } else {

          /* Nearest CCC */
          for (i = nx1; i <= nx2; i++) {
            NEAREST_IX_IY;
            if (ix >= 1 && ix <= nx && iy >= 1 && iy <= ny) {
              NEAREST_ARRAY_DEN;
              PATCH_CCC_SUM;
            }
          }
        }
    
      } else {

        if (oldDiff) {

          /* Linear simple difference */
          for (i = nx1; i <= nx2; i++) {
            LINEAR_IX_IY;
            if (ix >= 1 && ix < nx && iy >= 1 && iy < ny) {
              LINEAR_ARRAY_DEN;
              *sx += fabs(den - dend);
              (*numPix)++;
            }
          }
        } else if (!ifCCC) {

          /* Linear patch difference */
          for (i = nx1; i <= nx2; i++) {
            LINEAR_IX_IY;
            if (ix >= 1 && ix < nx && iy >= 1 && iy < ny) {
              LINEAR_ARRAY_DEN;
              PATCH_SUM;
            }
          }
        } else {

          /* Linear CCC */
          for (i = nx1; i <= nx2; i++) {
            LINEAR_IX_IY;
            if (ix >= 1 && ix < nx && iy >= 1 && iy < ny) {
              LINEAR_ARRAY_DEN;
              PATCH_CCC_SUM;
            }
          }
        }
      }
      
    }
  }
}

void simplexdiff(float *crray, float *drray, int *nx, int *ny, int *nx1, int *nx2,
                 int *ny1, int *ny2, float *amat, float *xcen, float *ycen, float *xAdd,
                 float *yAdd, int *ifInterp, int *oldDiff, int *ifCCC, double *sx, 
                 int *numPix, int *numXpatch, int *nxyPatch, double *sxa, double *sya, 
                 double *sxya, double *sxSqa, double *sySqa, int *numPixA)
{
  simplexDiff(crray, drray, *nx, *ny, *nx1, *nx2, *ny1,
              *ny2, amat, *xcen, *ycen, *xAdd, *yAdd, 
              *ifInterp, *oldDiff, *ifCCC, sx, numPix, *numXpatch,
              *nxyPatch, sxa, sya, sxya, sxSqa,
              sySqa, numPixA);
}
