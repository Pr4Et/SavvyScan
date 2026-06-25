/*  regression.c - replacements for old Fortran regression routines, and robust regression
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 2012 by Boulder Laboratory for 3-Dimensional Electron
 *  Microscopy of Cells ("BL3DEMC") and the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 * $Id$
 */
#include "imodconfig.h"
#include "b3dutil.h"

/* Macro to allow data matrix to be accessed with row, column indices regardless
   of its order */
#define XRC(r, c) (x[(r) * rowStride + (c) * colStride])

#ifdef F77FUNCAP
#define statmatrices STATMATRICES
#define multregress MULTREGRESS
#define robustregress ROBUSTREGRESS
#define multregressnoc MULTREGRESSNOC
#define robustregressnoc ROBUSTREGRESSNOC
#define polynomialfit POLYNOMIALFIT
#define weightedpolyfit WEIGHTEDPOLYFIT
#else
#define statmatrices statmatrices_
#define multregress multregress_
#define robustregress robustregress_
#define multregressnoc multregressnoc_
#define robustregressnoc robustregressnoc_
#define polynomialfit polynomialfit_
#define weightedpolyfit weightedpolyfit_
#endif

/*!
 * Computes basic statistical values and matrices from a data matrix representing a series
 * of measurements of multiple variables.  ^
 * Input parameters:  ^
 * [x]       - Input data matrix  ^
 * [xsize]   - Size of the fastest-progressing dimension of [x]  ^
 * [colFast] - Nonzero if the column dimension is the fastest progressing one  ^
 * [m]       - Number of columns of data, i.e., number of parameters  ^
 * [msize]   - Size of one dimension of the various square output matrices  ^
 * [ndata]   - Number of rows of data; i.e., number of separate measurements ^
 * [ifdisp]  - Value indicating either to skip computation of [d] and [r] (if 0), or to
 * treat the {m+1} column of data as weighting values (if < 0) ^
 * Outputs:
 * [sx]      - Array for the sums of the [m] data columns  ^
 * [ss]      - Array for raw sums of squares and cross-products  ^
 * [ssd]     - Array for sums of deviation squares and cross-products  ^
 * [d]       - Array for variances and covariances  (dispersion matrix) ^
 * [r]       - Array for matrix of correlation coefficients  ^
 * [xm]      - Array for the means of the [m] data columns  ^
 * [sd]      - Array for the standard deviations of the [m] data columns  ^
 * The output matrices will be treated as having leading dimension [msize], so they must 
 * be at least [msize] x [m] in size.  The data element at row i, column j would
 * be accessed as x\[i + j * xsize\] or x\[j\]\[i\] from C, or as x(i,j) from Fortran if
 * [colFast] is 0, or as x\[j + i * xsize\], x\[i\]\[j\], or x(j,i) if [colFast] is 
 * nonzero.  When  weighting is used, the returned values are weighted means, and other 
 * statistics with the weights incorporated.
 */
void statMatrices(float *x, int xsize, int colFast, int m, int msize, int ndata,
                  float *sx, float *ss, float *ssd, float *d, float *r, float *xm,
                  float *sd, int ifdisp)
{
  float fndata, den, weight, wsum;
  int i, j, k;
  int colStride = colFast ? 1 : xsize;
  int rowStride = colFast ? xsize : 1;
  
  fndata = ndata;
  
  for (i = 0; i< m; i++) {
    sx[i] = 0.;
    for (j = 0; j < m; j++) {
      ssd[msize * i + j] = 0.;
      r[msize * i + j] = 0.;
    }
  }
  
  if (ifdisp >= 0) {
    for (i = 0; i< m; i++) {
      for (k = 0; k < ndata; k++)
        sx[i] += XRC(k, i);
      xm[i] = sx[i] / fndata;
    }
  } else {
    wsum = 0.;
    for (k = 0; k < ndata; k++)
      wsum += XRC(k, m);
    for (i = 0; i< m; i++) {
      for (k = 0; k < ndata; k++)
        sx[i] += XRC(k, i) * XRC(k, m);
      xm[i] = sx[i] / wsum;
    }
  }
  
  for (k = 0; k < ndata; k++) {
    weight = 1.;
    if(ifdisp < 0)
      weight = XRC(k, m);
    for (i = 0; i< m; i++)
      for (j = 0; j < m; j++)
        ssd[i * msize + j] += (XRC(k, i) - xm[i]) * (XRC(k, j) - xm[j]) * weight;
  }
  
  for (i = 0; i< m; i++) {
    sd[i] = (float)sqrt((double)(ssd[i * msize + i] / (fndata-1.)));
    for (j = 0; j < m; j++) {
      ss[i * msize + j] = ssd[i * msize + j] + sx[i] * sx[j] / fndata;
      ss[j * msize + i] = ss[i * msize + j];
      ssd[j * msize + i] = ssd[i * msize + j];
    }
  }
  if(ifdisp == 0)
    return;
  for (i = 0; i< m; i++) {
    for (j = 0; j < m; j++) {
      d[i * msize + j] = ssd[i * msize + j] / (fndata-1.);
      d[j * msize + i] = d[i * msize + j];
      den = sd[i] * sd[j];
      r[i * msize + j] = 1.;
      if(den > 1.e-30)
        r[i * msize + j] = d[i * msize + j] / (sd[i] * sd[j]);
      r[j * msize + i] = r[i * msize + j];
    }
  }
}

/*!
 * Fortran wrapper for @statMatrices
 */
void statmatrices(float *x, int *xsize, int *colFast, int *m, int *msize, int *ndata,
                  float *sx, float *ss, float *ssd, float *d, float *r, float *xm,
                  float *sd, int *ifdisp)
{
  statMatrices(x, *xsize, *colFast, *m, *msize, *ndata, sx, ss, ssd, d, r, xm, sd,
               *ifdisp);
}


/*!
 * Computes a multiple linear regression (least-squares fit) for the relationships between
 * one or more dependent (output) variables and a set of independent (input) variables.  ^
 * Input parameters:  ^
 * [x]        - Input data matrix  ^
 * [xSize]    - Size of the fastest-progressing dimension of [x]  ^
 * [colFast]  - Nonzero if the column dimension is the fastest progressing one, i.e. if
 * successive values in the array occur in successive columns  ^
 * [numInCol]  - Number of columns of data for input variables.  This is allowed to be 0 ^
 * [numData]   - Number of rows of data; i.e., number of separate measurements ^
 * [numOutCol] - Number of columns of data for output variables; i.e., number of 
 * relationships to fit  ^
 * [wgtCol]    - Column number with weighting factors if > 0, otherwise no weighting.  
 * Columns are numbered from 0 when calling from C, or 1 calling from Fortran. ^
 * [solSize]   - Size of the fastest-progressing dimension of [sol], the array/matrix to
 * receive the solutions; must be at least [numInCol]  ^
 * [work]      - Any array for temporary use whose size must be at least 
 * (numInCol + numOutCol) * (numInCol + numOutCol)  floating point elements  ^
 * Outputs:  ^
 * [sol]       - Matrix to receive the [numOutCol] sets of [numInCol] coefficients.  Each 
 * set is placed in a column of [sol], where the row dimension is the fastest progressing 
 * one ^
 * [cons]      - Array to receive the [numOutCol] constant terms of the fits, or NULL
 * to fit an equation with no constant term  ^
 * [xMean]     - Array for the means of the [numInCol] plus [numOutCol] data columns  ^
 * [xSD]       - Array for the standard deviations of the [numInCol] plus [numOutCol] data
 * columns (these will be far from correct if [cons] is NULL)  ^
 * The input variables should be placed in the first [numInCol] columns of [x] and the
 * the output variables in the next [numOutCol] columns. The data element at row i,
 * column j would be accessed as x\[i + j * xSize\] or x\[j\]\[i\] from C, or as x(i,j) 
 * from Fortran if [colFast] is 0, or as x\[j + i * xSize\], x\[i\]\[j\], or x(j,i) if 
 * [colFast] is nonzero.  ^
 * The return value is 1 if [wgtCol] has an inappropriate value, or 3 if @gaussj returns
 * with an error.
 */
int multRegress(float *x, int xSize, int colFast, int numInCol, int numData, 
                int numOutCol, int wgtCol, float *sol, int solSize, float *cons,
                float *xMean, float *xSD, float *work)
{
  float fndata, den;
  double dsum, wsum;
  int i, j, k, mp;
  int colStride = colFast ? 1 : xSize;
  int rowStride = colFast ? xSize : 1;
  mp = numInCol + numOutCol;
  fndata = numData;
  if (wgtCol > 0 && (wgtCol < mp || (!colFast && wgtCol >= xSize)))
    return 1;

  /* Get the unweighted means */
  if (wgtCol <= 0) {
    for (i = 0; i< mp; i++) {
      dsum = 0.;
      for (k = 0; k < numData; k++)
        dsum += XRC(k, i);
      xMean[i] = dsum / fndata;
    }
  } else {

    /* Or add up the weights and get the weighted means */
    wsum = 0.;
    for (k = 0; k < numData; k++)
      wsum += XRC(k, wgtCol);
    for (i = 0; i< mp; i++) {
      dsum = 0.;
      for (k = 0; k < numData; k++)
        dsum += XRC(k, i) * XRC(k, wgtCol);
      xMean[i] = dsum / wsum;
    }
  }

  /* Get the sums of squares and cross-products of deviations */
  for (i = 0; i < mp; i++) {
    for (j = i; j < mp; j++) {
      if (i >= numInCol && i != j)
        continue;
      dsum = 0.;
      if (cons) {
        if (wgtCol <= 0) {
          for (k = 0; k < numData; k++)
            dsum += (XRC(k, i) - xMean[i]) * (XRC(k, j) - xMean[j]);
        } else {
          for (k = 0; k < numData; k++)
            dsum += (XRC(k, i) - xMean[i]) * (XRC(k, j) - xMean[j]) * XRC(k, wgtCol);
        }
      } else {
        if (wgtCol <= 0) {
          for (k = 0; k < numData; k++)
            dsum += XRC(k, i) * XRC(k, j);
        } else {
          for (k = 0; k < numData; k++)
            dsum += XRC(k, i) * XRC(k, j) * XRC(k, wgtCol);
        }
      }
      work[j * mp + i] = dsum;
    }
  }
    
  /* Get the SDs */
  for (i = 0; i< mp; i++)
    xSD[i] = (float)sqrt((double)(work[i * mp + i] / (fndata-1.)));

  /* If numInCol = 0, then just return the c values as the weighted means */
  if (!numInCol) {
    if (cons)
      for (j = 0; j < numOutCol; j++)
        cons[j] = xMean[j];
    return 0;
  }

  /* Scale it by n -1 to get covariance and by SD's to get the correlation matrix */
  for (i = 0; i < numInCol; i++) {
    for (j = i; j < mp; j++) {
      den = xSD[i] * xSD[j];
      if (den < 1.e-30) {
        /* printf("sds %d %d %f %f  den %f\n", i, j, xSD[i], xSD[j], den);
           fflush(stdout); */
        work[j * mp + i] = 1.;
      } else
        work[j * mp + i] /= den * (fndata-1.);
      if (j < numInCol)
        work[i * mp + j] = work[j * mp + i];
    }
  }
  /* for (i = 0; i < numInCol; i++) {
    for (j = 0; j < mp; j++)
      printf("  %.8f", work[j * mp + i]);
    printf("\n");
    } */
  
  /* The matrix to be solved is now completely filled and symmetric so row/column
     doesn't matter, but to call gaussj we need to transpose the final columns into sol */
  for (j = 0; j < numOutCol; j++) {
    for (i = 0; i < numInCol; i++)
      sol[j + i * numOutCol] = work[(j + numInCol) * mp + i];
    /*printf("\nsol:");
    for (i = 0; i < numInCol; i++)
    printf("    %g", sol[j + i * numOutCol]); */
  }

  if (gaussj(work, numInCol, mp, sol, numOutCol, numOutCol))
    return 3;
  
  /*for (j = 0; j < numOutCol; j++) {
    printf("\nsol sol:");
    for (i = 0; i < numInCol; i++)
      printf("    %g", sol[j + i * numOutCol]);
  }
  printf("\n");
  fflush(stdout);*/

  /* Scale the coefficients and transpose them back; get constant terms */
  memcpy(work, sol, numInCol * numOutCol * sizeof(float));
  for (j = 0; j < numOutCol; j++) {
    if (cons)
      cons[j] = xMean[numInCol + j];
    for (i = 0; i< numInCol; i++) {
      if (xSD[i] < 1.e-30)
        sol[i + solSize * j] = 0.;
      else
        sol[i + solSize * j] = work[j + i * numOutCol] * xSD[numInCol + j] / xSD[i];
      if (cons)
        cons[j] -= sol[i + solSize * j] * xMean[i];
    }
  }
  return 0;
}

/* Fortran wrapper */
int multregress(float *x, int *xSize, int *colFast, int *numInCol, int *numData, 
                int *numOutCol, int *wgtCol, float *sol, int *solSize, float *cons,
                float *xMean, float *xSD, float *work)
{
  return multRegress(x, *xSize, *colFast, *numInCol, *numData, *numOutCol, *wgtCol - 1, 
                     sol, *solSize, cons, xMean, xSD, work);
}

/*! Fortran wrapper to call @multRegress to fit with no constant term */
int multregressnoc(float *x, int *xSize, int *colFast, int *numInCol, int *numData, 
                int *numOutCol, int *wgtCol, float *sol, int *solSize,
                float *xMean, float *xSD, float *work)
{
  return multRegress(x, *xSize, *colFast, *numInCol, *numData, *numOutCol, *wgtCol - 1, 
                     sol, *solSize, NULL, xMean, xSD, work);
}

/* LAPACK NOTE:  To use ssysv, need to add
ssysv.$(O) ssytrf.$(O) ssytrs.$(O) slasyf.$(O) ssytf2.$(O)
slasyf.f  ssysv.f  ssytf2.f  ssytrf.f  ssytrs.f
to lapack, and
sscal.$(O) sger.$(O)  sgemv.$(O) sswap.$(O) scopy.$(O) isamax.$(O)  sgemm.$(O)\
        ssyr.$(O)
scopy.f  sgemm.f  sgemv.f  sger.f  sscal.f  sswap.f  ssyr.f isamax.f
to blas.  To use sspsv, need to add
ssptrs.$(O) ssptrf.$(O) sspsv.$(O)
sspsv.f  ssptrf.f  ssptrs.f  
to lapack and
sspr.$(O)
sspr.f
to blas

But hey, we could just ask for a double work array here.
The impediment is not doubles, it's having a different library.
Also, a bigger work array is needed not just for doubles but also for dsysv.
*/

/*!
 * Uses multiple linear regression to fit a polynomial of order
 * [order] to [ndata] points whose (x,y) coordinates are in the arrays [x] and [y].
 * It returns the coefficient of x to the i power in the array [slopes] and a
 * constant term in [intcpt].  The equation fit is:  ^
 * Y = intcpt + slopes\[0\] * X + slopes\[1\] * X**2 + ...  ^
 * [work] is an array whose size must be at least ([order] + 1) * ([order] + 3 + [ndata]).
 * The return value is the value returned by @@multRegress@.  Note that a Fortran 
 * function polyfit in libhvem takes care of allocating [work] to the needed size and 
 * calling the Fortran wrapper to this function.
 */
int polynomialFit(float *x, float *y, int ndata, int order, float *slopes, float *intcpt,
                  float *work)
{
  int wdim = order + 1;
  int i, j;
  float *xMean = work + wdim * ndata;
  float *xSD = xMean + wdim;
  float *mwork = xSD + wdim;
  if (!order)
    return 1;
  for (i = 0; i < ndata; i++) {
    for (j = 0; j < order; j++)
      work[i + j * ndata] = (float)pow((double)x[i], j + 1.);
    work[i + order * ndata] = y[i];
  }
  return (multRegress(work, ndata, 0, order, ndata, 1, 0, slopes, ndata, intcpt, xMean,
                      xSD, mwork));
}

/* Fortran wrapper */
int polynomialfit(float *x, float *y, int *ndata, int *order, float *slopes, 
                  float *intcpt, float *work)
{
  return polynomialFit(x, y, *ndata, *order, slopes, intcpt, work);
}

/*!
 * Uses multiple linear regression to fit a polynomial of order
 * [order] to [ndata] points whose (x,y) coordinates are in the arrays [x] and [y].
 * It returns the coefficient of x to the i power in the array [slopes] and a
 * constant term in [intcpt].  The equation fit is:  ^
 * Y = intcpt + slopes\[0\] * X + slopes\[1\] * X**2 + ...  ^
 * [work] is an array whose size must be at least ([order] + 2) * [ndata] +
 * ([order] + 1) * ([order] + 3)).  ^
 * The return value is the value returned by @@multRegress@.  This function is untested.
 */
int weightedPolyFit(float *x, float *y, float *weight, int ndata, int order,
                    float *slopes, float *intcpt, float *work)
{
  int wdim = order + 1;
  int i, j;
  float *xMean = work + (order + 2) * ndata;
  float *xSD = xMean + wdim;
  float *mwork = xSD + wdim;
  if (!order)
    return 1;
  for (i = 0; i < ndata; i++) {
    for (j = 0; j < order; j++)
      work[i + j * ndata] = (float)pow((double)x[i], j + 1.);
    work[i + order * ndata] = y[i];
    work[i + (order + 1) * ndata] = weight[i];
  }
  return (multRegress(work, ndata, 0, order, ndata, 1, order + 1, slopes, ndata, intcpt,
                      xMean, xSD, mwork));
}

/* Fortran wrapper */
int weightedpolyfit(float *x, float *y, float *weight, int *ndata, int *order,
                    float *slopes, float *intcpt, float *work)
{
  return weightedPolyFit(x, y, weight, *ndata, *order, slopes, intcpt, work);
}

/*!
 * Computes a robust least squares fit by iteratively computing a weight from the residual
 * for each data point then doing a weighted regression.  The weight is computed by
 * finding the median and normalized median absolute deviation (MADN) of the residual
 * values.  When there are multiple columns of dependent (output) variables, the square 
 * root of the sum of the squares of the residuals for the different variables is used.
 * Each residual is standardized by taking (residual - median) / MADN, and the weight is
 * computed from the Tukey bisquare equation using the standardized residual divided by 
 * the specified [kfactor].  However, with multiple output variables, all residuals 
 * are positive and ones less than the median are given a weighting of 1.
 * Input parameters:  ^
 * [x]         - Input data matrix  ^
 * [xSize]     - Size of the fastest-progressing dimension of [x]  ^
 * [colFast]   - Nonzero if the column dimension is the fastest progressing one  ^
 * [numInCol]  - Number of columns of data for input variables, which is allowed to be 0 ^
 * [numData]   - Number of rows of data; i.e., number of separate measurements ^
 * [numOutCol] - Number of columns of data for output variables; i.e., number of 
 * relationships to fit.   ^
 * [solSize]   - Size of the fastest-progressing dimension of [sol], the array/matrix to
 * receive the solutions; must be at least [numInCol]  ^
 * [work]      - Any array for temporary use whose size must be at least the maximum of
 * (numInCol + numOutCol) * (numInCol + numOutCol) and 2 * numData floating point 
 * elements ^
 * [kfactor]   - Factor by which to divide the standardized residual value in computing 
 * the Tukey bisquare weight.  4.68 is the typical value; a different value may be more 
 * appropriate with multiple dependent variables.  If the factor is negative, the 
 * absolute value is used, and a list of data row numbers to be given zero initial weights
 * can be passed in the [work] array.  The first element of [work] would have the number
 * or rows to initialize, and the following values would have row indexes, numbered from
 * 0 even when calling from Fortran.  ^
 * [maxIter]   - Maximum number of iterations, or the negative of the maximum to get 
 * trace output on each iteration.  20-50 is probably adequate.  With a negative value,
 * the program outputs the mean and maximum change in weighting and the number of points
 * with weights of 0, between 0 and 0.1, between 0.1 and 0.2, and less than 0.5.   ^
 * [maxZeroWgt] - Maximum number of points allowed to have zero weights.  When this number
 * is exceeded on an iteration, the deviations are sorted and the criterion deviation is 
 * permanently changed to midway 
 * between the deviations for the last point to be retained and first one to be 
 * eliminated.  ^
 * [maxChange]  - Maximum change in weights from one iteration to the next or across two 
 * iterations; the entry should not be smaller than 0.01.  ^
 * [maxOscill]  - Maximum change in weights from one iteration to the next, even if 
 * oscillating. The iterations are terminated when the biggest change in weights between 
 * iterations is less than [maxChange], or when it is less than [maxOscill] and the 
 * biggest change across two iterations is less than [maxChange].  ^
 * Outputs:  ^
 * [sol]       - Matrix to receive the [numOutCol] sets of [numInCol] coefficients.  Each 
 * set is placed in a column of [sol], where the row dimension is the fastest 
 * progressing one ^
 * [cons]      - Array to receive the [numOutCol] constant terms of the fits  ^
 * [xMean]     - Array for the means of the [numInCol] plus [numOutCol] data columns  ^
 * [xSD]       - Array for the standard deviations of the [numInCol] plus [numOutCol] data
 * columns  ^
 * [numIter] - Number of iterations  ^
 * The input variables should be placed in the first [numInCol] columns of [x] and
 * the output variables in the next [numOutCol] columns (see @@multRegress@).
 * Final weights are returned in the column of [x] after the output variables, and the
 * column after that is used for weights on the previous iteration, so [x]
 * must have at least [numInCol] + [numOutCol] + 2 columns.  The unweighted and
 * weighted root-mean-square residual errors are returned in the first and second 
 * elements of [work]. ^
 * The return value is 1 if there are not enough columns for the weights (detectable 
 * only when [colFast] is nonzero), 2 if there is a problem with the specification of
 * rows to be given initial weights of zero in [work], 3 if @gaussj returns with an 
 * error, or -1 if the weights do not converge.  ^
 * The procedure is described in Beaton, A.E., and Tukey, J.W., 1974, The fitting of 
 * power series, meaning polynomials, illustrated on band-spectroscopic data. 
 * Technometrics 16: 147-185, and in Gross, A. M., 1977, Confidence intervals for 
 * bisquare regression estimates. Journal of the American Statistical Association
 * 72: 341-354.  A value of 4.685 for [kfactor] gives a 95% efficiency for normally 
 * distributed errors and a single output variable (Holland, P.W., and Welsch, R.E., 1977,
 * Robust regression using iteratively reweighted least-squares. Communications in 
 * Statistics - Theory and Methods, 6: 813-827).
 */
int robustRegress(float *x, int xSize, int colFast, int numInCol, int numData,
                  int numOutCol, float *sol, int solSize, float *cons, float *xMean,
                  float *xSD, float *work, float kfactor, int *numIter, int maxIter, 
                  int maxZeroWgt, float maxChange, float maxOscill)
{
  int i, j, numOut, ierr, k, iter, num1, num2, num5, report = 0;
  int wgtCol = numInCol + numOutCol;
  int prevCol = wgtCol + 1;
  int colStride = colFast ? 1 : xSize;
  int rowStride = colFast ? xSize : 1;
  int splitLastTime = 0, keepCriterion = 0, specialWeights = 0;
  float minNonZeroWgt = 0.02;
  float diff, diffsum, diffmax, median, MADN, dev, weight, ressum, colres;
  float prev, prevmax, criterion, rmsErr, wgtRmsErr;
  if (kfactor < 1) {
    kfactor = -kfactor;
    specialWeights = 1;
  }
  if (maxIter < 0) {
    report = 1;
    maxIter = - maxIter;
  }

  if (colFast && prevCol >= xSize)
    return 1;

  /* Initialize weights to 1. */
  for (j = 0; j < numData; j++)
    XRC(j, wgtCol) = XRC(j, prevCol) = 1.;
  if (specialWeights) {
    num1 = B3DNINT(work[0]);
    if (num1 < 1 || num1 > numData / 2)
      return 2;
    for (i = 1; i <= num1; i++) {
      j = B3DNINT(work[i]);
      if (j < 0 || j >= numData)
        return 2;
      XRC(j, wgtCol) = XRC(j, prevCol) = 0.;
    }
  }

  /* Iterate */
  for (iter = 0; iter < maxIter; iter++) {

    /* Get regression solution */
    ierr = multRegress(x, xSize, colFast, numInCol, numData, numOutCol, wgtCol, sol, 
                       solSize, cons, xMean, xSD, work);
    if (ierr)
      return ierr;

    /* Compute residuals.  Note WORK has to be bigger for this! */
    rmsErr = 0.;
    wgtRmsErr = 0.;
    for (j = 0; j < numData; j++) {
      if (numOutCol == 1) {
        colres = (cons ? cons[0] : 0.) - XRC(j, numInCol);
        for (i = 0; i < numInCol; i++)
          colres += XRC(j, i) * sol[i];
        work[j] = colres;
      } else {
        ressum = 0.;
        for (k = 0; k < numOutCol; k++) {
          colres = (cons ? cons[k] : 0.) - XRC(j, (numInCol + k));
          for (i = 0; i < numInCol; i++)
            colres += XRC(j, i) * sol[i + k * solSize];
          /* printf("   %g", colres); */
          ressum += colres * colres;
        }
        work[j] = (float)sqrt(ressum);
        /* printf("   resid  %g\n", work[j]); */
      }
      rmsErr += work[j] * work[j];
      wgtRmsErr += work[j] * work[j] * XRC(j, wgtCol);
    }
    rmsErr = sqrt(rmsErr / numData);
    wgtRmsErr = sqrt(wgtRmsErr / numData);

    /* Get the median and MADN */
    rsMedian(work, numData, &work[numData], &median);
    rsMADN(work, numData, median, &work[numData], &MADN);
    if (!keepCriterion)
      criterion = kfactor * MADN;
    /* printf("Median %g   MADN %g\n", median, MADN); */

    /* Get the new weights and evaluate change from last time */
    diffsum = 0.;
    diffmax = 0.;
    numOut = 0;
    num1 = 0;
    num2 = 0;
    num5 = 0;
    prevmax = 0.;
    
    /* Convert to the deviations that will be compared to criterion and count zero's */
    for (j = 0; j < numData; j++) {
      work[j] -= median;

      /* For multiple columns, negative deviations are like no deviation */
      if (numOutCol > 1)
        work[j] = B3DMAX(0., work[j]);
      else if (work[j] < 0)
        work[j] = -work[j];
      if (work[j] > criterion)
        numOut++;
    }

    /* If there are too many out, need to adjust criterion to limit that number */
    if (numOut > maxZeroWgt) {
      for (j = 0; j < numData; j++)
        work[numData + j] = work[j];
      rsSortFloats(&work[numData], numData);
      diff = criterion;

      /* Either set it between the last one above 0 and the first one at 0, or if none
         ar allowed to be 0, set it to give the minimum weight */
      if (maxZeroWgt > 0) 
        criterion = (work[2 * numData - maxZeroWgt] + work[2 * numData - maxZeroWgt - 1])
          / 2.;
      else
        criterion = work[2 * numData - 1] / sqrt(1. - sqrt(minNonZeroWgt));
      if (report)
        printf("%d with zero weight, revising criterion from %g to %g\n", numOut,
               diff, criterion);
      keepCriterion = 1;
    }

    numOut = 0;
    for (j = 0; j < numData; j++) {

      /* Tests to handle case of MADN zero; avoid division by 0 */
      if (work[j] > criterion) {
        weight = 0.;
        numOut++;
      } else if (work[j] <= 1.e-6 * criterion) {
        weight = 1.;
      } else {
        dev = work[j] / criterion;
        weight = (1 - dev * dev) * (1 - dev * dev);
      }
      if (report) {
        if (weight > 0 && weight <= 0.1)
          num1++;
        if (weight > 0.1 && weight <= 0.2)
          num2++;
        if (weight < 0.5)
          num5++;
      }

      /* Get differences from last time and from time before that */
      diff = fabs(weight - XRC(j, wgtCol));
      diffsum += diff;
      diffmax = B3DMAX(diffmax, diff);
      prev = fabs(weight - XRC(j, prevCol));
      prevmax = B3DMAX(prevmax, prev);
      XRC(j, prevCol) = XRC(j, wgtCol);
      XRC(j, wgtCol) = weight;
    }
    if (report) {
      printf("Iter %3d del mean %.4f max %.4f prev %.4f # 0, 0.1, 0.2, <0.5: %d %d %d "
             "%d\n", iter, diffsum / numData, diffmax, prevmax, numOut, num1, num2, num5);
      fflush(stdout);
    }

    /* Quit if difference from last time is below limit, or if it is below the oscillation
       limit and difference from the previous time is below limit */
    if (!splitLastTime && (diffmax < maxChange || 
                           (diffmax < maxOscill && prevmax < maxChange)))
      break;

    /* But if we are oscillating, try half-way between and defer testing for 2 rounds */
    if (!splitLastTime && prevmax < maxChange / 2.) {
      splitLastTime = 1;
      if (report) {
        printf("        Oscillation detected: averaging previous weights\n");
        fflush(stdout);
      }
      for (j = 0; j < numData; j++)
        XRC(j, wgtCol) = (XRC(j, wgtCol) + XRC(j, prevCol)) / 2.;
    } else if (splitLastTime == 1)
      splitLastTime = 2;
    else
      splitLastTime = 0;
  }

  *numIter = iter;
  work[0] = rmsErr;
  work[1] = wgtRmsErr;
  return (*numIter < maxIter ? 0 : -1);
}

/* Fortran wrapper */
int robustregress(float *x, int *xSize, int *colFast, int *numInCol, int *numData, 
                  int *numOutCol, float *sol, int *solSize, float *cons, float *xMean, 
                  float *xSD, float *work, float *kfactor, int *numIter, int *maxIter, 
                  int *maxZeroWgt, float *maxChange, float *maxOscillate)
{
  return robustRegress(x, *xSize, *colFast, *numInCol, *numData, *numOutCol, sol,
                       *solSize, cons, xMean, xSD, work, *kfactor, numIter, 
                       *maxIter, *maxZeroWgt, *maxChange, *maxOscillate);
}

/*! Fortran wrapper to call @robustRegress to fit with no constant term */
int robustregressnoc(float *x, int *xSize, int *colFast, int *numInCol, int *numData, 
                  int *numOutCol, float *sol, int *solSize, float *xMean, 
                  float *xSD, float *work, float *kfactor, int *numIter, int *maxIter, 
                  int *maxZeroWgt, float *maxChange, float *maxOscillate)
{
  return robustRegress(x, *xSize, *colFast, *numInCol, *numData, *numOutCol, sol,
                       *solSize, NULL, xMean, xSD, work, *kfactor, numIter, 
                       *maxIter, *maxZeroWgt, *maxChange, *maxOscillate);
}

