/*
 *  $Id$
 *
 *  Author: David Mastronarde  email: mast@colorado.edu
 *
 *  Copyright (C) 1995-2014 by the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 */

#include <stdio.h>
#include <string.h>
#include "parse_params.h"
#include "imodel.h"
#include "warpfiles.h"

#define P2I_NO_FLIP          1
#define P2I_IGNORE_ZERO      2
#define P2I_COUNT_LINES      4
#define P2I_DISPLAY_VALUES   8
#define P2I_READ_WARP       16
#define P2I_TIMES_FOR_Z     32

#define MAX_VALUE_COLS 6
#define DEFAULT_SCALE 10.0

static Imod *imod_from_patches(FILE *fin, float scale, int clipSize, char *name, 
                               int flags, int valColOut);

static void usage(char *prog)
{
  imodVersion(prog);
  imodCopyright();
  printf("Usage: %s [options] patch_file output_model\n", prog);
  printf("Options:\n");
  printf("\t-s #\tScale vectors by given value (default %.1f)\n", DEFAULT_SCALE);
  printf("\t-f\tDo NOT flip the Y and Z coordinates\n");
  printf("\t-n name\tAdd given name to model object\n");
  printf("\t-c #\tSet up clipping planes enclosing area of given size\n");
  printf("\t-d #\tSet flag to display values in false color in 3dmod\n");
  printf("\t-v #[,#]\tValue column numbers (1-%d) to store as value and value2\n", 
         MAX_VALUE_COLS);
  printf("\t-z\tIgnore zero values when using SD to limit stored maximum value\n");
  printf("\t-l\tUse all lines in file; do not get line count from first line\n");
  printf("\t-w\tRead input file as a warping transformation file\n");
  printf("\t-t\tGive each contour a time equal to its Z value plus 1\n");

  exit(1);
}

int main( int argc, char *argv[])
{
  int i;
  FILE *fin = NULL, *fout;
  Imod *Model;
  float scale = 10.0;
  int clipSize = 0;
  int flags = 0;
  int nxWarp, nyWarp, nzWarp, warpFlags, version, indWarp, ibin;
  int valColOut = -1;
  float warpPixel;
  char *name = NULL;

  /* This name is hard-coded because of the script wrapper needed in Vista */
  char progname[] = "patch2imod";

  if (argc < 3)
    usage(progname);

  setExitPrefix("ERROR: patch2imod - ");
  
  for (i = 1; i < argc ; i++){
    if (argv[i][0] == '-'){
      switch (argv[i][1]){
        
      case 's':
        sscanf(argv[++i], "%f", &scale);
        break;

      case 'n':
        name = strdup(argv[++i]);
        break;

      case 'c':
        clipSize = atoi(argv[++i]);
        break;

      case 'f':
        flags |= P2I_NO_FLIP;
        break;

      case 'z':
        flags |= P2I_IGNORE_ZERO;
        break;

      case 'l':
        flags |= P2I_COUNT_LINES;
        break;

      case 'd':
        flags |= P2I_DISPLAY_VALUES;
        break;

      case 'w':
        flags |= P2I_READ_WARP;
        break;

      case 't':
        flags |= P2I_TIMES_FOR_Z;
        break;

      case 'v':
        valColOut = atoi(argv[++i]);
        if (!valColOut || valColOut < -MAX_VALUE_COLS)
          exitError("Value column # or ID must be positive and between -1 and -%d", 
                    MAX_VALUE_COLS);
        break;

      default:
        exitError("Illegal argument %s", argv[i]);
        break;
      }
    } else
      break;
  }
  if (i > (argc - 2)) {
    printf("ERROR: patch2imod - Wrong # of arguments\n");
    usage(progname);
  }

  if ((flags & P2I_COUNT_LINES) && valColOut > 0)
    exitError("You cannot enter a column ID if there is no first line with data count");

  if (flags & P2I_READ_WARP) {
    indWarp = readWarpFile(argv[i++], &nxWarp, &nyWarp, &nzWarp, &ibin, &warpPixel,
                           &version, &warpFlags);
    if (indWarp < 0)
      exitError("Reading %s as a warping file", argv[--i]);
                           
  } else {

    fin = fopen(argv[i++], "r");
    if (!fin)
      exitError("Couldn't open %s", argv[--i]);
  }

  if (imodBackupFile(argv[i])) {
    exitError("Renaming existing output file to %s~", argv[i]);
    exit(1);
  }
  fout = fopen(argv[i], "wb");
  if (!fout)
    exitError("Could not open %s", argv[i]);
  Model = (Imod *)imod_from_patches(fin, scale, clipSize, name, flags, valColOut);
     
  imodWrite(Model, fout);

  imodDelete(Model);
  exit(0);
}


#define MAXLINE 128

static Imod *imod_from_patches(FILE *fin, float scale, int clipSize, char *name,
                               int flags, int valColOut)
{
  int noflip = flags & P2I_NO_FLIP;
  int ignoreZero = flags & P2I_IGNORE_ZERO;
  int countLines = flags & P2I_COUNT_LINES;
  int displayValues = flags & P2I_DISPLAY_VALUES;
  int readWarp = flags & P2I_READ_WARP;
  int timesForZ = flags & P2I_TIMES_FOR_Z;
  int len, i, valueIDs[MAX_VALUE_COLS], valTypeMap[MAX_VALUE_COLS];
  int orderedIDs[MAX_VALUE_COLS];
  float values[MAX_VALUE_COLS];
     
  char line[MAXLINE];
  char *str, *token, *endptr;
  Imod *mod;
  Iobj *obj;
  Istore store;
  Ipoint *pts;
  int ix, iy, iz, itmp, ind, numValIDs;
  long longVal;
  int nxWarp, nyWarp, ifControl, nxGrid, nyGrid, numControl, maxProd, izWarp, pat;
  float xStart, yStart, xInterval, yInterval;
  int xmin, ymin, zmin, xmax, ymax, zmax;
  float dx, dy, dz, xx, yy, value, tmp;
  float valMin[MAX_VALUE_COLS], valMax[MAX_VALUE_COLS];
  double valSum[MAX_VALUE_COLS], valSqSum[MAX_VALUE_COLS];
  int numVals[MAX_VALUE_COLS];
  double valsd, sdmax;
  float *xVector, *yVector, *xControl, *yControl;
  float *dxGrid = NULL, *dyGrid = NULL;
  int residuals = 0;
  int dzvary = 0;
  int npatch = 0;
  int nzWarp = 1;
  int nread = 0;
  int contBase = 0;
  int maxValCols = 0;

  if (readWarp) {

    /* Read in warping file */
    getWarpFileSize(&nxWarp, &nyWarp, &nzWarp, &ifControl);
    maxProd = 0;
    for (iz = 0; iz < nzWarp; iz++) {
      if (ifControl)
        ix = getNumWarpPoints(iz, &numControl);
      else
        ix = getWarpGridSize(iz, &nxGrid, &nyGrid, &numControl);
      if (ix)
        exitError("Getting number of points in warping file at section %d", iz);
      npatch += numControl;
      maxProd = B3DMAX(maxProd, numControl);
    }
    if (!ifControl) {
      dxGrid = B3DMALLOC(float, maxProd);
      dyGrid = B3DMALLOC(float, maxProd);
      if (!dxGrid || !dyGrid)
        exitError("getting memory for warping grids");
    }
  } else if (countLines) {

    /* Count the lines in the file to get # of patches */
    while (1) {
      ix = fgetline(fin,line,MAXLINE);
      if (ix > 2)
        npatch++;
      else
        break;
    }
    if (npatch < 1)
      exitError("No usable lines in the file");
    rewind(fin);
  } else {

    /* Parse the first line for the # of patches and up to 6 value IDs */
    /* Any number of non-integer entries is skipped here */
    fgetline(fin,line,MAXLINE);
    if (strstr(line, "residuals") != NULL)
      residuals = 1;
    str = line;
    ind = 0;
    numValIDs = 0;
    for (ind = 0; ; ind++) {
      token = strtok(str, " ");
      if (!token)
        break;
      str = NULL;
      longVal = strtol(token, &endptr, 10);
      if (token[0] != 0x00 && *endptr == 0x00) {
        if (!ind)
          npatch = longVal;
        else if (numValIDs < MAX_VALUE_COLS)
          valueIDs[numValIDs++] = longVal;
      } else if (!ind)
        exitError("Converting the first item on the first line to get number of patches");
    }

    if (npatch < 1)
      exitError("Implausible number of patches = %d.", npatch);

    /* Convert the valColOut ID or -# to a column index */
    if (valColOut < 0) {
      valColOut = -valColOut - 1;
    } else {
      for (ind = 0; ind < numValIDs; ind++) {
        if (valueIDs[ind] == valColOut) {
          printf("The value with ID %d is in extra column %d\n", valColOut, ind + 1);
          valColOut = ind;
          break;
        }
      }
      if (ind == numValIDs)
        exitError("There is no value column ID # corresponding to the entered ID #");
    }
    /*printf("# vla ID %d %d %d %d\n", numValIDs, valueIDs[0],valueIDs[1],valueIDs[2]);*/

    /* Make the lists of ordered types IDs and map from column to type */
    orderedIDs[0] = valueIDs[valColOut];
    valTypeMap[valColOut] = 0;
    ind = 0;
    for (i = 1; i < MAX_VALUE_COLS; i++) {
      if (ind == valColOut)
        ind++;
      orderedIDs[i] = valueIDs[ind++];
    }
    ind = 1;
    for (i = 0; i < MAX_VALUE_COLS; i++) {
      if (i == valColOut)
        continue;
      valTypeMap[i] = ind++;
    }
  }
  /* printf("type map %d %d %d %d %d %d\n", valTypeMap[0], valTypeMap[1], valTypeMap[2],
            valTypeMap[3], valTypeMap[4], valTypeMap[5]); */

  mod = imodNew();     
  if (!mod)
    exitError("Could not get new model");

  if (imodNewObject(mod))
    exitError("Could not get new object");
  obj = mod->obj;
  obj->contsize = npatch;
  obj->cont = imodContoursNew(npatch);
  if (!obj->cont)
    exitError("Could not get contour array");
  obj->flags |= IMOD_OBJFLAG_OPEN;
  if (!residuals && !noflip && !readWarp)
    mod->flags |= IMODF_FLIPYZ;
  if (timesForZ)
    obj->flags |= IMOD_OBJFLAG_TIME;
  mod->pixsize = scale;
  xmin = ymin= zmin = 1000000;
  xmax = ymax = zmax = -1000000;
  for (i = 0; i < MAX_VALUE_COLS; i++) {
    valMin[i] = 1.e30;
    valMax[i] = -1.e30;
    numVals[i] = 0;
    valSum[i] = 0.;
    valSqSum[i] = 0.;
  }
  dz = 0.;
  store.flags = GEN_STORE_FLOAT << 2;

  for (izWarp = 0; izWarp < nzWarp; izWarp++) {
    if (readWarp) {
      iz = izWarp;
      if (ifControl) {
        if (getNumWarpPoints(iz, &npatch) || 
            getWarpPointArrays(iz, &xControl, &yControl, &xVector, &yVector))
          exitError("Getting number of control points or arrays for section %d", iz);
      } else {
        if (getWarpGrid(iz, &nxGrid, &nyGrid, &xStart, &yStart, &xInterval, &yInterval,
                        dxGrid, dyGrid, 0))
          exitError("Getting warp grid for section %d", iz);
        npatch = nxGrid * nyGrid;
      }
    }

    for (pat = 0; pat < npatch; pat++) {
      pts = (Ipoint *)malloc(2 * sizeof(Ipoint));
      if (!pts)
        exitError("Could not get new point array");

      obj->cont[pat + contBase].pts = pts;
      obj->cont[pat + contBase].psize = 2;
      if (timesForZ)
        obj->cont[pat + contBase].time = B3DMAX(1, iz + 1);
      if (readWarp) {
        if (ifControl) {
          xx = xControl[pat];
          yy = yControl[pat];
          dx = -xVector[pat];
          dy = -yVector[pat];
        } else {
          ix = pat % nxGrid;
          iy = pat / nxGrid;
          xx = xStart + ix * xInterval;
          yy = yStart + iy * yInterval;
          dx = -dxGrid[pat];
          dy = -dyGrid[pat];
        }
      } else {

        len = fgetline(fin,line, MAXLINE);
        if (len < 3)
          exitError("Error reading file at line %d.", pat + 1);
        
        /* DNM 7/26/02: read in residuals as real coordinates, without a 
           flip */
        if (residuals) {
          nread = sscanf(line, "%f %f %d %f %f", &xx, &yy, &iz, &dx, &dy);
        } else {
          /* DNM 11/15/01: have to handle either with commas or without,
             depending on whether it was produced by patchcorr3d or 
             patchcrawl3d */
          if (strchr(line, ','))
            nread = sscanf(line, "%d %d %d %f, %f, %f", &ix, &iz, &iy, &dx, &dz, &dy);
          else
            nread = sscanf(line, "%d %d %d %f %f %f %f %f %f %f %f %f", &ix, &iz, &iy, 
                           &dx, &dz, &dy, &values[0], &values[1], &values[2], &values[3],
                           &values[4], &values[5]);
          if (noflip) {
            itmp = iy;
            iy = iz;
            iz = itmp;
            tmp =dy;
            dy = dz;
            dz = tmp;
          }
          xx = ix;
          yy = iy;
        }
      }
      pts[0].x = xx;
      pts[0].y = yy;
      pts[0].z = iz;
      pts[1].x = xx + scale * dx;
      pts[1].y = yy + scale * dy;
      pts[1].z = iz + scale * dz;
      xmin = B3DMIN(xmin, xx);
      ymin = B3DMIN(ymin, yy);
      zmin = B3DMIN(zmin, iz);
      xmax = B3DMAX(xmax, xx + 1.);
      ymax = B3DMAX(ymax, yy + 1.);
      zmax = B3DMAX(zmax, iz + 1.);
      if (dz != 0.)
        dzvary = 1;

      maxValCols = B3DMAX(maxValCols, nread - 6);
      for (i = 0; i < nread - 6; i++) {
        ind = valTypeMap[i];
        value = values[i];
        valMin[ind] = B3DMIN(valMin[ind], value);
        valMax[ind] = B3DMAX(valMax[ind], value);
        store.type = GEN_STORE_VALUE1 + 2 * ind;
        if (value || !ignoreZero) {
          numVals[ind]++;
          valSum[ind] += value;
          valSqSum[ind] += value * value;
        }
        store.index.i = pat;
        store.value.f = value;
        if (istoreInsert(&obj->store, &store))
          exitError("Could not add general storage item");
      }
    }
    contBase += npatch;
  }
  /* printf("mvc %d %d %d %d %d %d %d\n", maxValCols, numVals[0], numVals[1], numVals[2], 
     numVals[3], numVals[4], numVals[5]); */
  
  if (residuals) {
    obj->symflags = IOBJ_SYMF_ARROW;
    obj->symbol = IOBJ_SYM_NONE;
    obj->symsize = 7;
  } else if (clipSize) {
    imodViewModelNew(mod);
    for (i = 0; i <  mod->viewsize; i++)
      mod->view[i].world |= WORLD_MOVE_ALL_CLIP;
    obj->clips.count = 4;
    obj->clips.flags = 0;
    obj->clips.normal[0].x = 1.;
    obj->clips.normal[1].x = -1.;
    obj->clips.normal[2].x = 0.;
    obj->clips.normal[3].x = 0.;
    obj->clips.normal[0].y = 0.;
    obj->clips.normal[1].y = 0.;
    obj->clips.normal[2].y = 1.;
    obj->clips.normal[3].y = -1.;
    for (i = 0; i < 4; i++) {
      obj->clips.normal[i].z = 0.;
      obj->clips.point[i].x = -0.5 * (xmax + xmin);
      obj->clips.point[i].y = -0.5 * (ymax + ymin);
      obj->clips.point[i].z = -0.5 * (zmax + zmin);
      if (obj->clips.normal[i].x)
        obj->clips.point[i].x += obj->clips.normal[i].x * clipSize / 2;
      if (obj->clips.normal[i].y)
        obj->clips.point[i].y += obj->clips.normal[i].y * clipSize / 2;
    }      
  }

  mod->xmax = xmax + xmin;
  mod->ymax = ymax + ymin;
  mod->zmax = zmax + zmin;
  if (dzvary)
    obj->symbol = IOBJ_SYM_CIRCLE;

  /* Set current thicken contour flag to aid deleting patches */
  obj->flags |= IMOD_OBJFLAG_THICK_CONT | IMOD_OBJFLAG_MCOLOR;
  if (displayValues)
    obj->flags |= IMOD_OBJFLAG_USE_VALUE;
  if (maxValCols) {
    maxValCols = B3DMAX(maxValCols, valColOut);
    obj->extra[IOBJ_EXSIZE - 1] = maxValCols;
  }
  for (i = 0; i < maxValCols; i++) {
    obj->extra[IOBJ_EXSIZE - 2 - i] = orderedIDs[i];
    if (numVals[i]) {
      if (numVals[i] > 10) {
        valsd = (valSqSum[i] - valSum[i] * valSum[i] / numVals[i]) / (numVals[i] - 1);
        if (valsd > 0) {
          sdmax = valSum[i] / numVals[i] + 10. * sqrt(valsd);
          valMax[i] = B3DMIN(valMax[i], sdmax);
        }
      }

      if (istoreAddMinMax(&obj->store, GEN_STORE_MINMAX1 + 2 * i, valMin[i], valMax[i]))
        exitError("Could not add general storage item");
    }
  }

  if (name)
    imodObjectSetName(obj, name);

  B3DFREE(dxGrid);
  B3DFREE(dyGrid);
  return(mod);
     
}
