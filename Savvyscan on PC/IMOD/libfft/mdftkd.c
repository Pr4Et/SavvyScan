/* mdftkd.c -   multi-dimensional complex fourier transform kernel driver
 * Original author: Lynn Ten Eyck
 * Translated to C: David Mastronarde
 * This translation is released under the LGPL.
 *
 *  $Id$
 */

#include "cfft.h"
#include <stdio.h>
#include <math.h>

void mdftkd (int n, int *factor, int *dim, float *x, float *y)
{
  int indFac, nReduced, pFac, stepAlong, sepAlong;

  sepAlong = dim[2];
  indFac = 0;
  nReduced = n;
  while (factor[indFac + 1] != 0) {
    indFac = indFac + 1;
    pFac = factor[indFac];
    nReduced = nReduced / pFac;
    stepAlong = nReduced * sepAlong;
    switch (pFac) {
    case 1:
      break;

    case 2:
      r2cftk (n, nReduced, x, y, &x[stepAlong], &y[stepAlong], dim);
      break;

    case 3:
      r3cftk (n, nReduced, x, y, &x[stepAlong], &y[stepAlong], &x[2 * stepAlong], 
              &y[2 * stepAlong], dim);
      break;

    case 4:
      r4cftk (n, nReduced, x, y, &x[stepAlong], &y[stepAlong], &x[2 * stepAlong],
              &y[2 * stepAlong], &x[3 * stepAlong], &y[3 * stepAlong], dim);
      break;

    case 5:
      r5cftk (n, nReduced, x, y, &x[stepAlong], &y[stepAlong], &x[2 * stepAlong], 
              &y[2 * stepAlong], &x[3 * stepAlong], &y[3 * stepAlong], &x[4 * stepAlong],
              &y[4 * stepAlong], dim);
      break;

    case 8:
      r8cftk (n, nReduced, x, y, &x[stepAlong], &y[stepAlong], &x[2 * stepAlong],
              &y[2 * stepAlong], &x[3 * stepAlong], &y[3 * stepAlong], &x[4 * stepAlong],
              &y[4 * stepAlong], &x[5 * stepAlong], &y[5 * stepAlong], &x[6 * stepAlong],
              &y[6 * stepAlong], &x[7 * stepAlong], &y[7 * stepAlong], dim);
      break;

    case 6:
      printf("\ntransfer error detected in mdftkd\n\n");
      return;

    case 7:
    default:
      rpcftk (n, nReduced, pFac, stepAlong, x, y, dim);
      break;
    }
  }
}

void r2cftk (int n, int nReduced, float *x0, float *y0, float *x1, float *y1,
             int *dim)
/*     radix 2 multi-dimensional complex fourier transform kernel */
{
  int fold, zero;
  int j, k, k0, nRed2, nRedOv2p1;
  int k1, sepBetween, kk, l, limAlong, sepNRed2, totFloats, extentBetween, sep;
  int nSep;
  double angle;
  float c, is, iu, rs, ru, sepAlong, twopi = 6.2831853;
  float fjm1, fnRed2;
  int itrip, ntrip;

  totFloats = dim[1];
  sep = dim[2];
  limAlong = dim[3];
  extentBetween = dim[4] - 1;
  sepBetween = dim[5];
  nSep = n * sep;
  nRed2 = nReduced * 2;
  fnRed2 = nRed2;
  nRedOv2p1 = nReduced / 2 + 1;
  sepNRed2 = sep * nRed2;

  fjm1 = -1.0;
  for (j = 1; j <= nRedOv2p1; j++) {
    fold = j > 1 && 2 * j < nReduced + 2;
    k0 = (j - 1) * sep + 1;
    fjm1 = fjm1 + 1.0;
    angle = twopi * fjm1 / fnRed2;
    zero = angle == 0.0;
    if (!zero) {
      c = cos(angle);
      sepAlong = sin(angle);
    }
    ntrip = fold ? 2 : 1;
    for (itrip = 0; itrip < ntrip; itrip++) {

      for (kk = k0; kk <= nSep; kk += sepNRed2) {
        for (l = kk; l <= totFloats; l += limAlong) {
          k1 = l + extentBetween;
          for (k = l - 1; k < k1; k += sepBetween) {
            rs = x0[k] + x1[k];
            is = y0[k] + y1[k];
            ru = x0[k] - x1[k];
            iu = y0[k] - y1[k];
            x0[k] = rs;
            y0[k] = is;
            if (!zero) {
              x1[k] = ru * c + iu * sepAlong;
              y1[k] = iu * c - ru * sepAlong;
            } else {
              x1[k] = ru;
              y1[k] = iu;
            }
          }
        }
      }
      k0 = (nReduced + 1 - j) * sep + 1;
      c = -c;
    }
  }
}

/*     radix 3 multi-dimensional complex fourier transform kernel */
void r3cftk (int n, int nReduced, float *x0, float *y0, float *x1, float *y1,
             float *x2, float *y2, int *dim)
{
  int fold, zero;
  int j, k, k0, nRed3, nRedOv2p1;
  int k1, sepBetween, kk, l, limAlong, sepNRed3, totFloats, extentBetween, sep;
  int nSep;
  double angle;
  float a = -0.5, b = 0.86602540, c1, c2, s1, s2, t, twopi = 6.2831853;
  float i0, i1, i2, ia, ib, is, r0, r1, r2, ra, rb, rs;
  float fjm1, fnRed3;
  int itrip, ntrip;
  /*      data twopi/6.2831853/, a/-0.5/, b/0.86602540/ */

  totFloats = dim[1];
  sep = dim[2];
  limAlong = dim[3];
  extentBetween = dim[4] - 1;
  sepBetween = dim[5];
  nSep = n * sep;
  nRed3 = nReduced * 3;
  fnRed3 = nRed3;
  sepNRed3 = sep * nRed3;
  nRedOv2p1 = nReduced / 2 + 1;

  fjm1 = -1.0;
  for (j = 1; j <= nRedOv2p1; j++) {
    fold = j > 1 && 2 * j < nReduced + 2;
    k0 = (j - 1) * sep + 1;
    fjm1 = fjm1 + 1.0;
    angle = twopi * fjm1 / fnRed3;
    zero = angle == 0.0;
    if (!zero) {
      c1 = cos(angle);
      s1 = sin(angle);
      c2 = c1 * c1 - s1 * s1;
      s2 = s1 * c1 + c1 * s1;
    }
    ntrip = fold ? 2 : 1;
    for (itrip = 0; itrip < ntrip; itrip++) {
      for (kk = k0; kk <= nSep; kk += sepNRed3) {
        for (l = kk; l <= totFloats; l += limAlong) {
          k1 = l + extentBetween;
          for (k = l - 1; k < k1; k += sepBetween) {
            r0 = x0[k];
            i0 = y0[k];
            rs = x1[k] + x2[k];
            is = y1[k] + y2[k];
            x0[k] = r0 + rs;
            y0[k] = i0 + is;
            ra = r0 + rs * a;
            ia = i0 + is * a;
            rb = (x1[k] - x2[k]) * b;
            ib = (y1[k] - y2[k]) * b;
            if (!zero) {
              r1 = ra + ib;
              i1 = ia - rb;
              r2 = ra - ib;
              i2 = ia + rb;
              x1[k] = r1 * c1 + i1 * s1;
              y1[k] = i1 * c1 - r1 * s1;
              x2[k] = r2 * c2 + i2 * s2;
              y2[k] = i2 * c2 - r2 * s2;
            } else {
              x1[k] = ra + ib;
              y1[k] = ia - rb;
              x2[k] = ra - ib;
              y2[k] = ia + rb;
            }
          }
        }
      }
      k0 = (nReduced + 1 - j) * sep + 1;
      t = c1 * a + s1 * b;
      s1 = c1 * b - s1 * a;
      c1 = t;
      t = c2 * a - s2 * b;
      s2 = -c2 * b - s2 * a;
      c2 = t;
    }
  }
}

/*     radix 4 multi-dimensional complex fourier transform kernel */
void r4cftk (int n, int nReduced, float *x0, float *y0, float *x1, float *y1,
             float *x2, float *y2, float *x3, float *y3, int *dim)
{
  int fold, zero;
  int j, k, k0, nRed4, nRedOv2p1;
  int k1, sepBetween, kk, l, limAlong, sepNRed4, totFloats, extentBetween, sep;
  int nSep;
  double angle;
  float c1, c2, c3, s1, s2, s3, t, twopi = 6.2831853;
  float i1, i2, i3, is0, is1, iu0, iu1, r1, r2, r3, rs0, rs1, ru0, ru1;
  float fjm1, fnRed4;
  int itrip, ntrip;

  totFloats = dim[1];
  sep = dim[2];
  limAlong = dim[3];
  extentBetween = dim[4] - 1;
  sepBetween = dim[5];
  nSep = n * sep;
  nRed4 = nReduced * 4;
  fnRed4 = nRed4;
  sepNRed4 = sep * nRed4;
  nRedOv2p1 = nReduced / 2 + 1;

  fjm1 = -1.0;
  for (j = 1; j <= nRedOv2p1; j++) {
    fold = j > 1 && 2 * j < nReduced + 2;
    k0 = (j - 1) * sep + 1;
    fjm1 = fjm1 + 1.0;
    angle = twopi * fjm1 / fnRed4;
    zero = angle == 0.0;
    if (!zero) {
      c1 = cos(angle);
      s1 = sin(angle);
      c2 = c1 * c1 - s1 * s1;
      s2 = s1 * c1 + c1 * s1;
      c3 = c2 * c1 - s2 * s1;
      s3 = s2 * c1 + c2 * s1;
    }
    ntrip = fold ? 2 : 1;
    for (itrip = 0; itrip < ntrip; itrip++) {


      for (kk = k0; kk <= nSep; kk += sepNRed4) {
        for (l = kk; l <= totFloats; l += limAlong) {

          k1 = l + extentBetween;
          for (k = l - 1; k < k1; k += sepBetween) {
            rs0 = x0[k] + x2[k];
            is0 = y0[k] + y2[k];
            ru0 = x0[k] - x2[k];
            iu0 = y0[k] - y2[k];
            rs1 = x1[k] + x3[k];
            is1 = y1[k] + y3[k];
            ru1 = x1[k] - x3[k];
            iu1 = y1[k] - y3[k];
            x0[k] = rs0 + rs1;
            y0[k] = is0 + is1;
            if (!zero) {
              r1 = ru0 + iu1;
              i1 = iu0 - ru1;
              r2 = rs0 - rs1;
              i2 = is0 - is1;
              r3 = ru0 - iu1;
              i3 = iu0 + ru1;
              x2[k] = r1 * c1 + i1 * s1;
              y2[k] = i1 * c1 - r1 * s1;
              x1[k] = r2 * c2 + i2 * s2;
              y1[k] = i2 * c2 - r2 * s2;
              x3[k] = r3 * c3 + i3 * s3;
              y3[k] = i3 * c3 - r3 * s3;
            } else {
              x2[k] = ru0 + iu1;
              y2[k] = iu0 - ru1;
              x1[k] = rs0 - rs1;
              y1[k] = is0 - is1;
              x3[k] = ru0 - iu1;
              y3[k] = iu0 + ru1;

            }
          }
        }
      }

      k0 = (nReduced + 1 - j) * sep + 1;
      t = c1;
      c1 = s1;
      s1 = t;
      c2 = -c2;
      t = c3;
      c3 = -s3;
      s3 = -t;
    }
  }
}


/*     radix 5 multi-dimensional complex fourier transform kernel */
void r5cftk (int n, int nReduced, float *x0, float *y0, float *x1, float *y1,
             float *x2, float *y2, float *x3, float *y3, float *x4, float *y4,
             int *dim)
{
  int fold, zero;
  int j, k, k0, nRed5, nRedOv2p1;
  int k1, sepBetween, kk, l, limAlong, sepNRed5, totFloats, extentBetween, sep;
  int nSep;
  double angle;
  float a1 = 0.30901699, a2 = -0.80901699, b1 = 0.95105652;
  float b2 = 0.58778525, c1, c2, c3, c4, s1, s2, s3, s4, t, twopi = 6.2831853;
  float r0, r1, r2, r3, r4, ra1, ra2, rb1, rb2, rs1, rs2, ru1, ru2;
  float i0, i1, i2, i3, i4, ia1, ia2, ib1, ib2, is1, is2, iu1, iu2;
  float fjm1, fnRed5;
  int itrip, ntrip;
  /*     data twopi/6.2831853/, a1/0.30901699/, b1/0.95105652/,
         .      a2/-0.80901699/, b2/0.58778525/ */

  totFloats = dim[1];
  sep = dim[2];
  limAlong = dim[3];
  extentBetween = dim[4] - 1;
  sepBetween = dim[5];
  nSep = n * sep;
  nRed5 = nReduced * 5;
  fnRed5 = nRed5;
  sepNRed5 = sep * nRed5;
  nRedOv2p1 = nReduced / 2 + 1;

  fjm1 = -1.0;
  for (j = 1; j <= nRedOv2p1; j++) {
    fold = j > 1 && 2 * j < nReduced + 2;
    k0 = (j - 1) * sep + 1;
    fjm1 = fjm1 + 1.0;
    angle = twopi * fjm1 / fnRed5;
    zero = angle == 0.0;
    if (!zero) {
      c1 = cos(angle);
      s1 = sin(angle);
      c2 = c1 * c1 - s1 * s1;
      s2 = s1 * c1 + c1 * s1;
      c3 = c2 * c1 - s2 * s1;
      s3 = s2 * c1 + c2 * s1;
      c4 = c2 * c2 - s2 * s2;
      s4 = s2 * c2 + c2 * s2;
    }
    ntrip = fold ? 2 : 1;
    for (itrip = 0; itrip < ntrip; itrip++) {

      for (kk = k0; kk <= nSep; kk += sepNRed5) {
        for (l = kk; l <= totFloats; l += limAlong) {
          k1 = l + extentBetween;
          for (k = l - 1; k < k1; k += sepBetween) {
            r0 = x0[k];
            i0 = y0[k];
            rs1 = x1[k] + x4[k];
            is1 = y1[k] + y4[k];
            ru1 = x1[k] - x4[k];
            iu1 = y1[k] - y4[k];
            rs2 = x2[k] + x3[k];
            is2 = y2[k] + y3[k];
            ru2 = x2[k] - x3[k];
            iu2 = y2[k] - y3[k];
            x0[k] = r0 + rs1 + rs2;
            y0[k] = i0 + is1 + is2;
            ra1 = r0 + rs1 * a1 + rs2 * a2;
            ia1 = i0 + is1 * a1 + is2 * a2;
            ra2 = r0 + rs1 * a2 + rs2 * a1;
            ia2 = i0 + is1 * a2 + is2 * a1;
            rb1 = ru1 * b1 + ru2 * b2;
            ib1 = iu1 * b1 + iu2 * b2;
            rb2 = ru1 * b2 - ru2 * b1;
            ib2 = iu1 * b2 - iu2 * b1;
            if (!zero) {
              r1 = ra1 + ib1;
              i1 = ia1 - rb1;
              r2 = ra2 + ib2;
              i2 = ia2 - rb2;
              r3 = ra2 - ib2;
              i3 = ia2 + rb2;
              r4 = ra1 - ib1;
              i4 = ia1 + rb1;
              x1[k] = r1 * c1 + i1 * s1;
              y1[k] = i1 * c1 - r1 * s1;
              x2[k] = r2 * c2 + i2 * s2;
              y2[k] = i2 * c2 - r2 * s2;
              x3[k] = r3 * c3 + i3 * s3;
              y3[k] = i3 * c3 - r3 * s3;
              x4[k] = r4 * c4 + i4 * s4;
              y4[k] = i4 * c4 - r4 * s4;
            } else {
              x1[k] = ra1 + ib1;
              y1[k] = ia1 - rb1;
              x2[k] = ra2 + ib2;
              y2[k] = ia2 - rb2;
              x3[k] = ra2 - ib2;
              y3[k] = ia2 + rb2;
              x4[k] = ra1 - ib1;
              y4[k] = ia1 + rb1;
            }
          }
        }
      }
      k0 = (nReduced + 1 - j) * sep + 1;
      t = c1 * a1 + s1 * b1;
      s1 = c1 * b1 - s1 * a1;
      c1 = t;
      t = c2 * a2 + s2 * b2;
      s2 = c2 * b2 - s2 * a2;
      c2 = t;
      t = c3 * a2 - s3 * b2;
      s3 = -c3 * b2 - s3 * a2;
      c3 = t;
      t = c4 * a1 - s4 * b1;
      s4 = -c4 * b1 - s4 * a1;
      c4 = t;
    }
  }
}


/*     radix 8 multi-dimensional complex fourier transform kernel */
void r8cftk (int n, int nReduced, float *x0, float *y0, float *x1, float *y1,
             float *x2, float *y2, float *x3, float *y3, float *x4,
             float *y4, float *x5, float *y5, float *x6, float *y6,
             float *x7, float *y7, int *dim)
{
  int fold, zero;
  int j, k, k0, nRed8, nRedOv2p1;
  int k1, sepBetween, kk, l, limAlong, sepNRed8, totFloats, extentBetween, sep;
  int nSep;
  double angle;
  float c1, c2, c3, c4, c5, c6, c7, e = 0.70710678;
  float s1, s2, s3, s4, s5, s6, s7, t, twopi = 6.2831853;
  float r1, r2, r3, r4, r5, r6, r7, rs0, rs1, rs2, rs3, ru0, ru1, ru2, ru3;
  float i1, i2, i3, i4, i5, i6, i7, is0, is1, is2, is3, iu0, iu1, iu2, iu3;
  float rss0, rss1, rsu0, rsu1, rus0, rus1, ruu0, ruu1;
  float iss0, iss1, isu0, isu1, ius0, ius1, iuu0, iuu1;
  float fjm1, fnRed8;
  int itrip, ntrip;
  /*      data twopi/6.2831853/, e/0.70710678/ */

  totFloats = dim[1];
  sep = dim[2];
  limAlong = dim[3];
  extentBetween = dim[4] - 1;
  sepBetween = dim[5];
  nSep = n * sep;
  nRed8 = nReduced * 8;
  fnRed8 = nRed8;
  sepNRed8 = sep * nRed8;
  nRedOv2p1 = nReduced / 2 + 1;

  fjm1 = -1.0;
  for (j = 1; j <= nRedOv2p1; j++) {
    fold = j > 1 && 2 * j < nReduced + 2;
    k0 = (j - 1) * sep + 1;
    fjm1 = fjm1 + 1.0;
    angle = twopi * fjm1 / fnRed8;
    zero = angle == 0.0;
    if (!zero) {
      c1 = cos(angle);
      s1 = sin(angle);
      c2 = c1 * c1 - s1 * s1;
      s2 = s1 * c1 + c1 * s1;
      c3 = c2 * c1 - s2 * s1;
      s3 = s2 * c1 + c2 * s1;
      c4 = c2 * c2 - s2 * s2;
      s4 = s2 * c2 + c2 * s2;
      c5 = c4 * c1 - s4 * s1;
      s5 = s4 * c1 + c4 * s1;
      c6 = c4 * c2 - s4 * s2;
      s6 = s4 * c2 + c4 * s2;
      c7 = c4 * c3 - s4 * s3;
      s7 = s4 * c3 + c4 * s3;
    }
    ntrip = fold ? 2 : 1;
    for (itrip = 0; itrip < ntrip; itrip++) {

      for (kk = k0; kk <= nSep; kk += sepNRed8) {
        for (l = kk; l <= totFloats; l += limAlong) {
          k1 = l + extentBetween;
          for (k = l - 1; k < k1; k += sepBetween) {
            rs0 = x0[k] + x4[k];
            is0 = y0[k] + y4[k];
            ru0 = x0[k] - x4[k];
            iu0 = y0[k] - y4[k];
            rs1 = x1[k] + x5[k];
            is1 = y1[k] + y5[k];
            ru1 = x1[k] - x5[k];
            iu1 = y1[k] - y5[k];
            rs2 = x2[k] + x6[k];
            is2 = y2[k] + y6[k];
            ru2 = x2[k] - x6[k];
            iu2 = y2[k] - y6[k];
            rs3 = x3[k] + x7[k];
            is3 = y3[k] + y7[k];
            ru3 = x3[k] - x7[k];
            iu3 = y3[k] - y7[k];
            rss0 = rs0 + rs2;
            iss0 = is0 + is2;
            rsu0 = rs0 - rs2;
            isu0 = is0 - is2;
            rss1 = rs1 + rs3;
            iss1 = is1 + is3;
            rsu1 = rs1 - rs3;
            isu1 = is1 - is3;
            rus0 = ru0 - iu2;
            ius0 = iu0 + ru2;
            ruu0 = ru0 + iu2;
            iuu0 = iu0 - ru2;
            rus1 = ru1 - iu3;
            ius1 = iu1 + ru3;
            ruu1 = ru1 + iu3;
            iuu1 = iu1 - ru3;
            t = (rus1 + ius1) * e;
            ius1 = (ius1 - rus1) * e;
            rus1 = t;
            t = (ruu1 + iuu1) * e;
            iuu1 = (iuu1 - ruu1) * e;
            ruu1 = t;
            x0[k] = rss0 + rss1;
            y0[k] = iss0 + iss1;
            if (!zero) {
              r1 = ruu0 + ruu1;
              i1 = iuu0 + iuu1;
              r2 = rsu0 + isu1;
              i2 = isu0 - rsu1;
              r3 = rus0 + ius1;
              i3 = ius0 - rus1;
              r4 = rss0 - rss1;
              i4 = iss0 - iss1;
              r5 = ruu0 - ruu1;
              i5 = iuu0 - iuu1;
              r6 = rsu0 - isu1;
              i6 = isu0 + rsu1;
              r7 = rus0 - ius1;
              i7 = ius0 + rus1;
              x4[k] = r1 * c1 + i1 * s1;
              y4[k] = i1 * c1 - r1 * s1;
              x2[k] = r2 * c2 + i2 * s2;
              y2[k] = i2 * c2 - r2 * s2;
              x6[k] = r3 * c3 + i3 * s3;
              y6[k] = i3 * c3 - r3 * s3;
              x1[k] = r4 * c4 + i4 * s4;
              y1[k] = i4 * c4 - r4 * s4;
              x5[k] = r5 * c5 + i5 * s5;
              y5[k] = i5 * c5 - r5 * s5;
              x3[k] = r6 * c6 + i6 * s6;
              y3[k] = i6 * c6 - r6 * s6;
              x7[k] = r7 * c7 + i7 * s7;
              y7[k] = i7 * c7 - r7 * s7;
            } else {
              x4[k] = ruu0 + ruu1;
              y4[k] = iuu0 + iuu1;
              x2[k] = rsu0 + isu1;
              y2[k] = isu0 - rsu1;
              x6[k] = rus0 + ius1;
              y6[k] = ius0 - rus1;
              x1[k] = rss0 - rss1;
              y1[k] = iss0 - iss1;
              x5[k] = ruu0 - ruu1;
              y5[k] = iuu0 - iuu1;
              x3[k] = rsu0 - isu1;
              y3[k] = isu0 + rsu1;
              x7[k] = rus0 - ius1;
              y7[k] = ius0 + rus1;
            }
          }
        }
      }
      k0 = (nReduced + 1 - j) * sep + 1;
      t = (c1 + s1) * e;
      s1 = (c1 - s1) * e;
      c1 = t;
      t = s2;
      s2 = c2;
      c2 = t;
      t = (-c3 + s3) * e;
      s3 = (c3 + s3) * e;
      c3 = t;
      c4 = -c4;
      t = -(c5 + s5) * e;
      s5 = (-c5 + s5) * e;
      c5 = t;
      t = -s6;
      s6 = -c6;
      c6 = t;
      t = (c7 - s7) * e;
      s7 = -(c7 + s7) * e;
      c7 = t;
    }
  }
}

/*     radix prime multi-dimensional complex fourier transform kernel */
void rpcftk (int n, int nReduced, int pFac, int stepAlong, float *x, float *y, int *dim)
{

  int fold, zero;
  double angle;
  float is, iu, rs, ru, t, twopi = 6.2831853, xt, yt;
  float fu, fp, fjm1, fnRedp;
  int j, jj, k0, k, nRedOv2p1, nRedp, pm, pp, u, v;
  int k1, sepBetween, kk, l, limAlong, sepNRedp, totFloats, extentBetween, sep;
  int nSep;

  float aa[10][10], bb[10][10];
  float a[19], b[19], c[19], sepAlong[19];
  float ia[10], ib[10], ra[10], rb[10];
  int itrip, ntrip;

  totFloats = dim[1];
  sep = dim[2];
  limAlong = dim[3];
  extentBetween = dim[4] - 1;
  sepBetween = dim[5];
  nSep = n * sep;
  nRedOv2p1 = nReduced / 2 + 1;
  nRedp = nReduced * pFac;
  fnRedp = nRedp;
  sepNRedp = sep * nRedp;
  pp = pFac / 2;
  pm = pFac - 1;
  fp = pFac;
  fu = 0.0;
  for (u = 1; u <= pp; u++) {
    fu = fu + 1.0;
    angle = twopi * fu / fp;
    jj = pFac - u;
    a[u] = cos(angle);
    b[u] = sin(angle);
    a[jj] = a[u];
    b[jj] = -b[u];
  }
  for (u = 1; u <= pp; u++) {
    for (v = 1; v <= pp; v++) {
      jj = u * v - u * v / pFac * pFac;
      aa[v][u] = a[jj];
      bb[v][u] = b[jj];

    }
  }

  fjm1 = -1.0;
  for (j = 1; j <= nRedOv2p1; j++) {
    fold = j > 1 && 2 * j < nReduced + 2;
    k0 = (j - 1) * sep + 1;
    fjm1 = fjm1 + 1.0;
    angle = twopi * fjm1 / fnRedp;
    zero = angle == 0.0;
    if (!zero) {
      c[1] = cos(angle);
      sepAlong[1] = sin(angle);
      for (u = 2; u <= pm; u++) {
        c[u] = c[u - 1] * c[1] - sepAlong[u - 1] * sepAlong[1];
        sepAlong[u] = sepAlong[u - 1] * c[1] + c[u - 1] * sepAlong[1];
      }
    }
    ntrip = fold ? 2 : 1;
    for (itrip = 0; itrip < ntrip; itrip++) {

      for (kk = k0; kk <= nSep; kk += sepNRedp) {
        for (l = kk; l <= totFloats; l += limAlong) {
          k1 = l + extentBetween;
          for (k = l - 1; k < k1; k += sepBetween) {
            xt = x[k];
            yt = y[k];
            rs = x[k + stepAlong] + x[k + stepAlong * pm];
            is = y[k + stepAlong] + y[k + stepAlong * pm];
            ru = x[k + stepAlong] - x[k + stepAlong * pm];
            iu = y[k + stepAlong] - y[k + stepAlong * pm];
            for (u = 1; u <= pp; u++) {
              ra[u] = xt + rs * aa[u][1];
              ia[u] = yt + is * aa[u][1];
              rb[u] = ru * bb[u][1];
              ib[u] = iu * bb[u][1];
            }
            xt = xt + rs;
            yt = yt + is;

            for (u = 2; u <= pp; u++) {    /* u numbers from 1 not 0 */
              jj = pFac - u;
              rs = x[k + u * stepAlong] + x[k + jj * stepAlong];
              is = y[k + u * stepAlong] + y[k + jj * stepAlong];
              ru = x[k + u * stepAlong] - x[k + jj * stepAlong];
              iu = y[k + u * stepAlong] - y[k + jj * stepAlong];
              xt = xt + rs;
              yt = yt + is;
              for (v = 1; v <= pp; v++) {
                ra[v] = ra[v] + rs * aa[v][u];
                ia[v] = ia[v] + is * aa[v][u];
                rb[v] = rb[v] + ru * bb[v][u];
                ib[v] = ib[v] + iu * bb[v][u];
              }
            }
            x[k] = xt;
            y[k] = yt;
            for (u = 1; u <= pp; u++) {
              jj = pFac - u;
              if (!zero) {
                xt = ra[u] + ib[u];
                yt = ia[u] - rb[u];
                x[k + u * stepAlong] = xt * c[u] + yt * sepAlong[u];
                y[k + u * stepAlong] = yt * c[u] - xt * sepAlong[u];
                xt = ra[u] - ib[u];
                yt = ia[u] + rb[u];
                x[k + jj * stepAlong] = xt * c[jj] + yt * sepAlong[jj];
                y[k + jj * stepAlong] = yt * c[jj] - xt * sepAlong[jj];
              } else {
                x[k + u * stepAlong] = ra[u] + ib[u];
                y[k + u * stepAlong] = ia[u] - rb[u];
                x[k + jj * stepAlong] = ra[u] - ib[u];
                y[k + jj * stepAlong] = ia[u] + rb[u];
              }
            }
          }
        }
      }
      if (!fold) break;
      k0 = (nReduced + 1 - j) * sep + 1;
      for (u = 1; u <= pm; u++) {
        t = c[u] * a[u] + sepAlong[u] * b[u];
        sepAlong[u] = -sepAlong[u] * a[u] + c[u] * b[u];
        c[u] = t;
      }
    }
  }
}
