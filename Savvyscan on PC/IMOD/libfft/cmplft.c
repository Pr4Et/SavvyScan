/*  cmplft.c -   complex finite discrete fourier transform
 * Translated to C: David Mastronarde
 * This translation is released under the LGPL.
 *
 *  $Id$
 */

/*
  C     COMPLEX FINITE DISCRETE FOURIER TRANSFORM
  C     TRANSFORMS ONE DIMENSION OF MULTI-DIMENSIONAL DATA
  C     MODIFIED BY L. F. TEN EYCK FROM A ONE-DIMENSIONAL VERSION WRITTEN
  C     BY G. T. SANDE, 1969.
  C
  C     THIS PROGRAM CALCULATES THE TRANSFORM
  C               (X(T) + I*Y(T))*(COS(2*PI*T/N) - I*SIN(2*PI*T/N))
  C
  C     INDEXING -- THE ARRANGEMENT OF THE MULTI-DIMENSIONAL DATA IS
  C     SPECIFIED BY THE INTEGER ARRAY D, THE VALUES OF WHICH ARE USED AS
  C     CONTROL PARAMETERS IN DO LOOPS.  WHEN IT IS DESIRED TO COVER ALL
  C     ELEMENTS OF THE DATA FOR WHICH THE SUBSCRIPT BEING TRANSFORMED HAS
  C     THE VALUE I0, THE FOLLOWING IS USED.
  C
  C               I1 = (I0 - 1)*D(2) + 1
  C               DO 100 I2 = I1, D(1), D(3)
  C               I3 = I2 + D(4) - 1
  C               DO 100 I = I2, I3, D(5)
  C                  .
  C                  .
  C           100 CONTINUE
  C
  C     WITH THIS INDEXING IT IS POSSIBLE TO USE A NUMBER OF ARRANGEMENTS
  C     OF THE DATA, INCLUDING NORMAL FORTRAN COMPLEX NUMBERS (D(5) = 2)
  C     OR SEPARATE STORAGE OF REAL AND IMAGINARY PARTS.
*/

#include "cfft.h"
#include <stdlib.h>
#include <stdio.h>

void cmplft (float *x, float *y, int n, int *d)
{
  int error;
  int pmax, psym, twogrp;
  int factor[16], sym[16], unsym[16];
  double wallStart, wallSrfp = 0., wallMk = 0., wallDi = 0.;

  /*
    c     pmax is the largest prime factor that will be tolerated by this
    c     program.
    c     twogrp is the largest power of two that is treated as a special
    c     case.
  */
  pmax = 19;
  twogrp = 8;

  if (n <= 1)
    return;
  wallStart = fftStartTimer();
  srfp (n, pmax, twogrp, factor, sym, &psym, unsym, &error);
  fftAddTime(wallStart, &wallSrfp);
  if (error) {
    printf ("invalid number of points for cmplft.  n = %d\n", n);
    exit(1);
  }
  wallStart = fftStartTimer();
  mdftkd (n, factor, d, x, y);
  fftAddTime(wallStart, &wallMk);
  wallStart = fftStartTimer();
  diprp (n, sym, psym, unsym, d, x, y);
  fftAddTime(wallStart, &wallDi);
#ifdef OLD_FFT_TIMES
  printf("srfp %.1f  mdftkd %.1f   diprp %.1f\n", wallSrfp * 1000., wallMk * 1000.,
         wallDi * 1000.);
#endif
}
