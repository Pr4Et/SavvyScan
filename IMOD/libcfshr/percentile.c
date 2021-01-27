/*
 * percentile.c - Selecting an item at a percentile in a list
 *
 * Copyright (C) 2008 by Boulder Laboratory for 3-Dimensional Electron
 * Microscopy of Cells ("BL3DEMC") and the Regents of the University of 
 * Colorado.  See dist/COPYRIGHT for full notice.
 *
 * $Id$
 *
 */
#include "b3dutil.h"

#ifdef F77FUNCAP
#define percentilefloat PERCENTILEFLOAT
#else
#define percentilefloat percentilefloat_
#endif

/* 
 * Routines for selecting item number s (numbered from 1) out of num items
 * Based on a pascal program apparently from Handbook of Algorithms and Data Structures,
 * by Gonnet and Baeza-Yates.
 * 10/13/14: Algorithm is highly pathological for finding an extreme percentile when
 * there are many identical values at that extreme.  Inverting the order of operations
 * for low percentiles solves the problem.
 */

/*!
 * Selects item number [s] (numbered from 1) out of [num] items in the array
 * [r], where items are considered in order from low to high.  [r] is partially
 * rearranged while finding the item.  The algorithm runs in linear time and is faster
 * than sorting the array first.
 */
float percentileFloat(int s, float *r, int num)
{
  int lo = 0;
  int up = num - 1;
  int i, j;
  float temp;
  s--;
  while (up >= s && s >= lo) {
    i = lo;
    j = up;
    temp = r[s];
    if (s > num / 2) {

      /* Operations for high percentiles */
      r[s] = r[lo];
      r[lo] = temp;
      while (i < j) {
        while (r[j] > temp)
          j--;
        r[i] = r[j];
        while (i < j && r[i] <= temp)
          i++;
        r[j] = r[i];
      }
      r[i] = temp;
      if (s < i)
        up = i - 1;
      else
        lo = i + 1;
    } else {

      /* Operations for low percentiles */
      r[s] = r[up];
      r[up] = temp;
      while (i < j) {
        while (r[i] < temp)
          i++;
        r[j] = r[i];
        while (i < j && r[j] >= temp)
          j--;
        r[i] = r[j];
      }
      r[j] = temp;
      if (s > j)
        lo = j + 1;
      else
        up = j - 1;
    }
  }

  return r[s];
}

/*!
 * Fortran wrapper to @percentileFloat 
 */
double percentilefloat(int *s, float *r, int *num)
{
  return (double)percentileFloat(*s, r, *num);
}

/*!
 * Same as @percentileFloat but with an integer array [r]
 */
int percentileInt(int s, int *r, int num)
{
  int lo = 0;
  int up = num - 1;
  int i, j;
  int temp;
  s--;
  while (up >= s && s >= lo) {
    i = lo;
    j = up;
    temp = r[s];
    if (s > num / 2) {
      r[s] = r[lo];
      r[lo] = temp;
      while (i < j) {
        while (r[j] > temp)
          j--;
        r[i] = r[j];
        while (i < j && r[i] <= temp)
          i++;
        r[j] = r[i];
      }
      r[i] = temp;
      if (s < i)
        up = i - 1;
      else
        lo = i + 1;
    } else {
      r[s] = r[up];
      r[up] = temp;
      while (i < j) {
        while (r[i] < temp)
          i++;
        r[j] = r[i];
        while (i < j && r[j] >= temp)
          j--;
        r[i] = r[j];
      }
      r[j] = temp;
      if (s > j)
        lo = j + 1;
      else
        up = j - 1;
    }
  }

  return r[s];
}
