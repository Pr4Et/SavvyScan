/*  hermft.c -   hermitian symmetric fourier transform
 * Original author: Lynn Ten Eyck
 * Translated to C: David Mastronarde
 * This translation is released under the LGPL.
 *
 *  $Id$
 */

/*
  C     HERMITIAN SYMMETRIC FOURIER TRANSFORM
  C
  C     GIVEN THE UNIQUE TERMS OF A HERMITIAN SYMMETRIC SEQUENCE OF LENGTH
  C     2N THIS SUBROUTINE CALCULATES THE 2N REAL NUMBERS WHICH ARE ITS
  C     FOURIER TRANSFORM.  THE EVEN NUMBERED ELEMENTS OF THE TRANSFORM
  C     (0, 2, 4, . . ., 2N-2) ARE RETURNED IN X AND THE ODD NUMBERED
  C     ELEMENTS (1, 3, 5, . . ., 2N-1) IN Y.
  C
  C     A FINITE HERMITIAN SEQUENCE OF LENGTH 2N CONTAINS N + 1 UNIQUE
  C     REAL NUMBERS AND N - 1 UNIQUE IMAGINARY NUMBERS.  FOR CONVENIENCE
  C     THE REAL VALUE FOR X(N) IS STORED AT Y(0).
*/

#include "cfft.h"
#include <math.h>

void hermft(float * x, float * y, int n, int * dim)
{
  float a, b, c, d, e, f, co, si, two_n, twopi;
  double angle ;
  int totFloats, sepAlong, limAlong, extentBetween, sepBetween;
  int i, j, k, nover2, i0, i1, i2;
  int k1;
  double wallStart, wallCum = 0.;

  twopi = 6.2831853;
  two_n = 2 * n;

  totFloats = dim[1];
  sepAlong = dim[2];
  limAlong = dim[3];
  extentBetween = dim[4] - 1;
  sepBetween = dim[5];
  
  wallStart = fftStartTimer();
  for (i0 = 1; i0 <= totFloats; i0 += limAlong) {
    i1 = i0 + extentBetween;
    for (i = i0 - 1; i < i1; i += sepBetween) {
      a = x[i];
      b = y[i];
      x[i] = a + b;
      y[i] = a - b;
    }
  }

  nover2 = n / 2 + 1;
  if (nover2 < 2)
    return;
  for (i0 = 2; i0 <= nover2; i0++) {
    angle = twopi * (i0 - 1) / two_n;
    co = cos(angle);
    si = sin(angle);
    k = (n + 2 - 2 * i0) * sepAlong;
    k1 = (i0 - 1) * sepAlong + 1;
    for (i1 = k1; i1 <= totFloats; i1 += limAlong) {
      i2 = i1 + extentBetween;
      for (i = i1 - 1; i < i2; i += sepBetween) {
        j = i + k;
        a = x[i] + x[j];
        b = x[i] - x[j];
        c = y[i] + y[j];
        d = y[i] - y[j];
        e = b * co + c * si;
        f = b * si - c * co;
        x[i] = a + f;
        x[j] = a - f;
        y[i] = e + d;
        y[j] = e - d;
      }
    }
  }
  fftAddTime(wallStart, &wallCum);
#ifdef OLD_FFT_TIMES
  printf("hermft loops %.1f\n", wallCum * 1000.);
#endif
  cmplft (x, y, n, dim);
}
