/*  realft.c -   real fourier transform
 * Original author: Lynn Ten Eyck
 * Translated to C: David Mastronarde
 * This translation is released under the LGPL.
 *
 *  $Id$
 */

/*
  C     REAL FOURIER TRANSFORM
  C
  C     GIVEN A REAL SEQUENCE OF LENGTH 2N THIS SUBROUTINE CALCULATES THE
  C     UNIQUE PART OF THE FOURIER TRANSFORM.  THE FOURIER TRANSFORM HAS
  C     N + 1 UNIQUE REAL PARTS AND N - 1 UNIQUE IMAGINARY PARTS.  SINCE
  C     THE REAL PART AT X(N) IS FREQUENTLY OF INTEREST, THIS SUBROUTINE
  C     STORES IT AT X(N) RATHER THAN IN Y(0).  THEREFORE X AND Y MUST BE
  C     OF LENGTH N + 1 INSTEAD OF N.  NOTE THAT THIS STORAGE ARRANGEMENT
  C     IS DIFFERENT FROM THAT EMPLOYED BY THE HERMITIAN FOURIER TRANSFORM
  C     SUBROUTINE.
  C
  C     FOR CONVENIENCE THE DATA IS PRESENTED IN TWO PARTS, THE FIRST
  C     CONTAINING THE EVEN NUMBERED REAL TERMS AND THE SECOND CONTAINING
  C     THE ODD NUMBERED TERMS (NUMBERING STARTING AT 0).  ON RETURN THE
  C     REAL PART OF THE TRANSFORM REPLACES THE EVEN TERMS AND THE
  C     IMAGINARY PART OF THE TRANSFORM REPLACES THE ODD TERMS.
*/

#include "cfft.h"
#include <math.h>

void realft (float * even, float * odd, int n, int * dim)
{
  float a, b, c, d, e, f, co, si, twopi, two_n;
  double angle;
  int totFloats, sepAlong, limAlong, extentBetween, sepBetween;
  int i, j, k, l, nover2, i0, i1, i2;
  double wallStart, wallCum = 0.;

  twopi = 6.2831853;
  two_n = 2 * n;     /* Number of elements in X in real image */

  cmplft (even, odd, n, dim);
  totFloats = dim[1];         /* Total number of floats */
  sepAlong = dim[2];          /* Stride along one row/col */
  limAlong = dim[3];          /* Limit on index along row/col */
  extentBetween = dim[4] - 1; /* Limit on index when striding between multiple row/cols */
  sepBetween = dim[5];        /* Stride between multiple rows/cols */
  nover2 = n / 2 + 1;

  wallStart = fftStartTimer();
  if (nover2 >= 2) {
    for (i = 2; i <= nover2; i++) {
      angle = twopi * (i - 1) / two_n;
      co = cos(angle);
      si = sin(angle);
      i0 = (i - 1) * sepAlong + 1;
      j = (n + 2 - 2 * i) * sepAlong;
      for (i1 = i0; i1 <= totFloats; i1 += limAlong) {
        i2 = i1 + extentBetween;
        for (k = i1 - 1; k < i2; k += sepBetween) {
          l = k + j;
          a = (even[l] + even[k]) / 2.0;
          c = (even[l] - even[k]) / 2.0;
          b = (odd[l] + odd[k]) / 2.0;
          d = (odd[l] - odd[k]) / 2.0;
          e = c * si + b * co;
          f = c * co - b * si;
          even[k] = a + e;
          even[l] = a - e;
          odd[k] = f - d;
          odd[l] = f + d;
        }
      }
    }
  }
  if (n < 1)
    return;
  j = n * sepAlong;
  for (i1 = 1; i1 <= totFloats; i1 += limAlong) {
    i2 = i1 + extentBetween;
    for (k = i1 - 1; k < i2; k += sepBetween) {
      l = k + j;
      even[l] = even[k] - odd[k];
      odd[l] = 0.0;
      even[k] = even[k] + odd[k];
      odd[k] = 0.0;
    }
  }
  fftAddTime(wallStart, &wallCum);
#ifdef OLD_FFT_TIMES
  printf("realft loops %.1f\n", wallCum * 1000.);
#endif
}
