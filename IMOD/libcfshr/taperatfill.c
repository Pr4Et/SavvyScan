/*
 * taperatfill.c - Taper a slice at fill areas on edges
 *
 * Copyright (C) 2006 by Boulder Laboratory for 3-Dimensional Electron
 * Microscopy of Cells ("BL3DEMC") and the Regents of the University of 
 * Colorado.  See dist/COPYRIGHT for full notice.
 *
 * $Id$
 */

#include "b3dutil.h"
#include "mrcslice.h"
#include "imodconfig.h"

#ifdef F77FUNCAP
#define taperatfill TAPERATFILL
#define getlasttaperfillvalue GETLASTTAPERFILLVALUE
#else
#define taperatfill taperatfill_
#define getlasttaperfillvalue getlasttaperfillvalue_
#endif

struct pixdist {
  short int dx, dy;
  int edgeind;
};
struct edgelist {
  int x, y;
};
#define MAX_TAPER 256
#define MAX_AVG_OUT 16

static int sFoundFill = 0;
static float sLastFillVal = 0.;

/*!
 * Analyzes for fill areas at the edge of the slice [sl], finds the borders 
 * between actual image and fill areas, and tapers the image intensities down 
 * to the fill value over [ntaper] pixels.  If [inside] is 0 the image 
 * intensities at the edge are extended into the fill area, while if [inside] 
 * is 1 then pixels inside the image are attenuated.  Returns -1 for memory
 * allocation errors.
 */
int sliceTaperAtFill(Islice *sl, int ntaper, int inside)
{
  Ival val, val2, val3, val4;
  float fracs[MAX_TAPER * 10], fillpart[MAX_TAPER * 10];
  float fillval, lastval, weight, wsum, valsum;
  int len, longest;
  int i, j, k, newInd, mm, ix, iy, xnext, ynext, xp, yp, longix, longiy, lastix, lastiy;
  int dir, taperDir, col, row, ncol, ind, xstart, ystart, found, dirstart;
  int dist, distmin[MAX_AVG_OUT], distmax, minedge[MAX_AVG_OUT], lastout;
  int rowLim, colLim, lim1, lim2, numMins, maxAvg = 15;
  struct pixdist *plist;
  struct edgelist *elist;
  int pllimit;
  int plsize = 0;
  int elsize = 0;
  int ellimit;
  int tapersq;
  unsigned char *pixmap;
  int xsize = sl->xsize;
  int ysize = sl->ysize;
  int pmsize = xsize * ysize;
  int dxout[4] = {0, 1, 0, -1};
  int dyout[4] = {-1, 0, 1, 0};
  int dxnext[4] = {1, 0, -1, 0};
  int dynext[4] = {0, 1, 0, -1};
  int numThreads, maxThreads = 16;   /* Validated to 12 on an old dual CPU Xeon */
  int plist_chunk, elist_chunk;
  /*double wallStart = wallTime(); */

  ntaper = B3DMIN(MAX_TAPER, ntaper);
  tapersq = (ntaper + 1) * (ntaper + 1);
  inside = inside != 0 ? 1 : 0;
  elist_chunk =  (sl->xsize + sl->ysize) / 5 + 1;
  plist_chunk = ntaper * elist_chunk;
  pllimit = plist_chunk;
  ellimit = elist_chunk;
  sFoundFill = 0;

  /* find longest string of repeated values along the edge */
  longest = 0;
  for (iy = 0; iy < ysize; iy += ysize - 1) {
    len = 0;
    sliceGetVal(sl, 0, iy, val);
    lastval = val[0];
    lastix = 0;
    lastiy = iy;
    for (ix = 1; ix < xsize; ix++) {
      sliceGetVal(sl, ix, iy, val);
      if (val[0] == lastval) {
        /* same as last, add to count, see if this is a new
           best count*/
        len++;
        if (len > longest) {
          longest = len;
          fillval = lastval;
          longix = lastix;
          longiy = lastiy;
          dir = iy ? 2 : 0;
        }
      } else {
        /* different, reset count and set new lastval */
        len = 0;
        lastval = val[0];
        lastix = ix;
      }
    }
  }         
  for (ix = 0; ix < xsize; ix += xsize - 1) {
    len = 0;
    sliceGetVal(sl, ix, 0, val);
    lastval = val[0];
    lastix = ix;
    lastiy = 0;
    for (iy = 1; iy < ysize; iy++) {
      sliceGetVal(sl, ix, iy, val);
      if (val[0] == lastval) {
        /* same as last, add to count, see if this is a new
           best count */
        len++;
        if (len > longest) {
          longest = len;
          fillval = lastval;
          longix = lastix;
          longiy = lastiy;
          dir = ix ? 1 : 3;
        }
      } else {
        /* different, reset count and set new lastval */
        len = 0;
        lastval = val[0];
        lastiy = iy;
      }
    }
  }         

  /* If length below a criterion (what?) , return without error */
  if (longest < 10)
    return 0;

  /* set the slice mean so that sliceGetValue will return this outside the
     edges */
  sl->mean = fillval;

  /* Start at the middle of that long interval and walk in until non-fill
     pixel is found */
  ix = longix + dxnext[dir % 2] * longest / 2;
  iy = longiy + dynext[dir % 2] * longest / 2;
  found = 0;
  /*printf("\nlongix %d longiy %d dir %d longest %d ix %d iy %d\n", longix,
    longiy, dir, longest, ix, iy);*/
  while (!found) {
    ix -= dxout[dir];
    iy -= dyout[dir];
    if (ix < 0 || ix >= xsize || iy < 0 || iy >= ysize)
      break;
    sliceGetVal(sl, ix, iy, val);
    if (val[0] != fillval)
      found = 1;
  }
  
  if (!found)
    return 0;

  sFoundFill = 1;
  sLastFillVal = fillval;

  /* Get initial chunk of pixel list and pixmap */
  plist = (struct pixdist *)malloc(plist_chunk * sizeof(struct pixdist));
  elist = (struct edgelist *)malloc(elist_chunk * sizeof(struct pixdist));
  if (plist && elist)
    pixmap = (unsigned char *)malloc(pmsize);
  if (!plist || !pixmap || !elist) {
    if (plist)
      free(plist);
    if (elist)
      free(elist);
    return (-1);
  }

  /* clear pixmap */
  for (i = 0; i < pmsize; i++)
    pixmap[i] = 0;

  dirstart = dir;
  xstart = ix;
  ystart = iy;
  taperDir = 1 - 2 * inside;
  elist[0].x = ix;
  elist[0].y = iy;
  elsize = 1;
  lastout = 1;
  /* printf("xstart %d ystart  %d\n", xstart, ystart);*/

  do {
    ncol = 1;
    xnext = ix + dxout[dir] + dxnext[dir];
    ynext = iy + dyout[dir] + dynext[dir];
    sliceGetVal(sl, xnext, ynext, val);
    ind = 1;
    if (val[0] != fillval) {
      /* case 1: inside corner */
      ix = xnext;
      iy = ynext;
      dir = (dir + 3) % 4;
      if (inside)
        ncol = ntaper + 1;
    } else {
      xnext = ix + dxnext[dir];
      ynext = iy + dynext[dir];
      sliceGetVal(sl, xnext, ynext, val);
      if (val[0] != fillval) {
        /* case 2: straight edge to next pixel */
        ix = xnext;
        iy = ynext;
      } else {
        /* case 3: outside corner, pixel stays the same */
        dir = (dir + 1) % 4;
        if (!inside)
          ncol = ntaper + 1;
        ind = lastout;
      }
    }
    /* printf ("%d %d %d %d %d\n", xnext, ynext, ix, iy, dir); */

    /* If outside pixel is outside the data, nothing to add to lists */
    xp = ix + dxout[dir];
    yp = iy + dyout[dir];
    lastout = (xp < 0 || xp >= xsize || yp < 0 || yp >= ysize) ? 1 : 0;
    if (lastout)
      continue;

    /* Add a new point to edge list */
    if (ind) {
      if (elsize >= ellimit) {
        ellimit += elist_chunk;
        elist = (struct edgelist *) realloc
          (elist, ellimit * sizeof(struct edgelist));
        if (!elist) {
          free(plist);
          free(pixmap);
          return (-1);
        }
      }
      elist[elsize].x = ix;
      elist[elsize++].y = iy;
    }
      
    /* Precompute the limits of the loop to stay within image and save a test in loop */
    if (dxnext[dir]) {
      lim1 = ix / dxnext[dir];
      lim2 = (ix + 1 - xsize) / dxnext[dir];
    } else {
      lim1 = iy / dynext[dir];
      lim2 = (iy + 1 - ysize) / dynext[dir];
    }
    colLim = ncol - 1;
    if (lim1 >= 0 && lim1 < colLim)
      colLim = lim1;
    if (lim2 >= 0 && lim2 < colLim)
      colLim = lim2;

    if (dxout[dir]) {
      lim1 = -ix / (taperDir * dxout[dir]);
      lim2 = (xsize - 1 - ix) / (taperDir * dxout[dir]);
    } else {
      lim1 = -iy / (taperDir * dyout[dir]);
      lim2 = (ysize - 1 - iy) / (taperDir * dyout[dir]);
    }
    rowLim = ntaper - inside;
    if (lim1 >= 1 - inside && lim1 < rowLim)
      rowLim = lim1;
    if (lim2 >= 1 - inside && lim2 < rowLim)
      rowLim = lim2;

    /* Loop on all the pixels to mark */
    for (col = 0; col <= colLim; col++) {
      for (row = 1 - inside; row <= rowLim; row++) {
        xp = ix + taperDir * row * dxout[dir] - col * dxnext[dir];
        yp = iy + taperDir * row * dyout[dir] - col * dynext[dir];

        /* Skip if already on list (moved test up 9/23/16) */
        ind = xsize * yp + xp;
        if (pixmap[ind])
          continue;

        /* skip marking outside pixels for inside taper or
           inside pixels for outside taper.  Notice not good assumption that
           anything with fill value is "outside" */
        sliceGetVal(sl, xp, yp, val);
        if ((inside && val[0] == fillval) || (!inside && val[0] != fillval))
          continue;

        /* If pixel is new, mark in pixmap and make a new entry on the list */
        /* Eliminate keeping track of distance and looking pixel up 10/14/08 */

        /* But first insist that the pixel have an outside neighbor if out */
        if (!inside) {
          sliceGetVal(sl, xp + 1, yp, val);
          sliceGetVal(sl, xp, yp + 1, val2);
          sliceGetVal(sl, xp - 1, yp, val3);
          sliceGetVal(sl, xp, yp - 1, val4);
          if (val[0] != fillval && val2[0] != fillval && val3[0] != fillval
              && val4[0] != fillval)
            continue;
        }

        pixmap[ind] = 1;
        if (plsize >= pllimit) {
          pllimit += plist_chunk;
          plist = (struct pixdist *) realloc
            (plist, pllimit * sizeof(struct pixdist));
          if (!plist) {
            free (pixmap);
            free(elist);
            return (-1);
          }
        }
        plist[plsize].edgeind = elsize - 1;
        plist[plsize].dx = xp - ix;
        plist[plsize++].dy = yp - iy;
      }
    }

  } while (ix != xstart || iy != ystart || dir != dirstart);

  /* make tables of fractions and amounts to add of fillval */
  for (i = 1; i <= ntaper * 10; i++) {
    dist = inside ? i : ntaper * 10 + 1 - i;
    fracs[i] = (float)dist / (ntaper * 10 + 1);
    fillpart[i] = (1. - fracs[i]) * fillval;
  }
  /* printf("Made pixel list size %d edge list size %d at %.3f\n", plsize, elsize, 
    wallTime() - wallStart);
    wallStart = wallTime();*/

  /* Process the pixels on the list. This loop had good parallel efficiency */
  numThreads = B3DNINT(sqrt((double)plsize) / 60.);
  B3DCLAMP(numThreads, 1, maxThreads);
  numThreads = numOMPthreads(numThreads);

#pragma omp parallel for num_threads(numThreads)                        \
  shared(plsize, plist, tapersq, elsize, inside, ntaper, elist, sl, fracs, fillpart) \
  shared(maxAvg) \
  private(i, distmin, minedge, distmax, ix, iy, dir, j, k, xp, yp, dist, ind, val) \
  private(newInd, mm, numMins, weight, wsum, valsum)
  for (i = 0; i < plsize; i++) {

    /* Go in both directions from starting index looking for closest distance*/
    distmin[0] = plist[i].dx * plist[i].dx + plist[i].dy * plist[i].dy;
    minedge[0] = plist[i].edgeind;
    numMins = 1;
    distmax = 1.5 * B3DMAX(distmin[0], tapersq);
    ix = elist[minedge[0]].x + plist[i].dx;
    iy = elist[minedge[0]].y + plist[i].dy;
    /*printf("%d,%d  start dist %d to %d,%d", ix,iy,distmin[0],elist[minedge[0]].x,
      elist[minedge[0]].y);*/
    for (dir = -1; dir <= 1; dir += 2) {
      k = plist[i].edgeind;
      for (j = 0; j < elsize; j++) {
        k = (k + dir + elsize) % elsize;
        xp = ix - elist[k].x;
        yp = iy - elist[k].y;
        dist = xp * xp + yp * yp;

        /* Take it as a new point if it is less than the current biggest distance
           saved or if there is not a full set of points saved yet for outside tapering */
        if (dist < distmin[numMins - 1] || (!inside && numMins < maxAvg)) {
          if (inside) {
            distmin[0] = dist;
            minedge[0] = k;
          } else {
            
            /* Find first point that it is closer than: there might not be any */
            for (newInd = B3DMIN(numMins, maxAvg - 1); newInd > 0; newInd--)
              if (dist >= distmin[newInd - 1])
                break;

            /* Shift the existing points up */
            for (mm = B3DMIN(maxAvg - 2, numMins - 1); mm >= newInd; mm--)  {
              distmin[mm + 1] = distmin[mm];
              minedge[mm + 1] = minedge[mm];
            }
            numMins = B3DMIN(numMins + 1, maxAvg);
            distmin[newInd] = dist;
            minedge[newInd] = k;
          }
        }
        if (dist > distmax)
          break;
      }
    }

    ind = 10. * (sqrt((double)distmin[0]) + inside);
    /*printf("   end dist %d to %d,%d", distmin[0],elist[minedge[0]].x,
      elist[minedge[0]].y);*/
    
    if (ind > ntaper * 10) {
      /*printf("\n");*/
      continue;
    }
    if (inside) {
      xp = ix;
      yp = iy;
      sliceGetVal(sl, xp, yp, val);
    } else {

      /* Form weighted sum */
      wsum = 0.;
      valsum = 0.;
      for (mm = 0; mm < numMins; mm++) {
        weight =  1. / (4. + sqrt(distmin[mm]));
        wsum += weight;
        sliceGetVal(sl, elist[minedge[mm]].x, elist[minedge[mm]].y, val);
        valsum += weight * val[0];
      }
      val[0] = valsum / wsum;
    }
    /*  val[0] = fracs[plist[i].dist] * val[0] + fillpart[plist[i].dist]; */
    val[0] = fracs[ind] * val[0] + fillpart[ind];
    slicePutVal(sl, ix, iy, val);
    /*printf("  frac %f  fill %f   put val %f\n", fracs[ind], fillpart[ind], val[0]);*/
  }
  /* printf("Tapering took %.3f\n", wallTime() - wallStart); */

  free(plist);
  free(pixmap);
  return 0;
}

/*! 
 * Fortran wrapper to @sliceTaperAtFill for tapering a real*4 image in [array]
 * of size [nx] by [ny].  Returns -1 for memory error.
 */
int taperatfill(float *array, int *nx, int *ny, int *ntaper, int *inside)
{
  Islice slice;
  sliceInit(&slice, *nx, *ny, SLICE_MODE_FLOAT, array);
  return sliceTaperAtFill(&slice, *ntaper, *inside);
}

/*!
 * Returns the fill value found on last call in [value] and a return value of 1 if a value
 * was found, or 0 if not.
 */
int getLastTaperFillValue(float *value)
{
  *value = sLastFillVal;
  return sFoundFill;
}

/*! Fortran wrapper to @getLastTaperFillValue */
int getlasttaperfillvalue(float *value)
{
  return getLastTaperFillValue(value);
}
