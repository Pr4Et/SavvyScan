/*
 *  $Id$
 *
 *  Author: David Mastronarde  email: mast@colorado.edu
 *
 *  Copyright (C) 1995-2014 by the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 */

#include <stdio.h>
#include <string.h>

#include "imodel.h"
#include "istore.h"
#include "parse_params.h"

static int istoreFindValue(Ilist *list, int index, int type, float *value,
                           int *listInd);
#define MAX_VALUE_COLS 6

int main( int argc, char *argv[])
{
  int i, ind, col, maxValType;
  FILE *fout;
  Imod *mod;
  int npatch = 0;
  int numValues[MAX_VALUE_COLS], valueIDs[MAX_VALUE_COLS];
  int ix, iy, iz;
  int colForVal1 = 1;
  float dx, dy, dz;
  float value, maxVal[MAX_VALUE_COLS];
  int ob, co, listInd, listStart, numCols;
  Ipoint *pts;
  char format[MAX_VALUE_COLS][10];
  int colToTypeMap[MAX_VALUE_COLS];

  setExitPrefix("ERROR: imod2patch - ");


  if (argc < 3){
    if (argc != 1)
      printf("ERROR: imod2patch - wrong # of arguments\n");
    printf("imod2patch usage:\n");
    printf("imod2patch [-v col] imod_model patch_file\n");
    exit(1);
  }
  for (ind = 0; ind < MAX_VALUE_COLS; ind++) {
    numValues[ind] = 0;
    valueIDs[ind] = 0;
    colToTypeMap[ind] = -1;
    maxVal[ind] = -1.e37;
  }

  i = 1;
  if (!strcmp(argv[i], "-v")) {
    i++;
    colForVal1 = atoi(argv[i++]);
  }

  mod = imodRead(argv[i]);
  if (!mod)
    exitError("Reading model %s\n", argv[i]);
  
  if (imodBackupFile(argv[++i]))
    exitError("Renaming existing output file to %s\n", argv[i]);
  
  fout = fopen(argv[i], "w");
  if (!fout)
    exitError("Could not open %s\n", argv[i]);

  /* Count up patches and values and find max of values */
  for (ob = 0; ob < mod->objsize; ob++) {
    listInd = 0;
    for (co = 0; co < mod->obj[ob].contsize; co++) {
      listStart = listInd;
      if (mod->obj[ob].cont[co].psize >= 2) {
        npatch++;
        for (ind = 0; ind < MAX_VALUE_COLS; ind++) {
          listInd = listStart;
          if (istoreFindValue(mod->obj[ob].store, co, GEN_STORE_VALUE1 + 2 * ind, &value,
                              &listInd)) {
            maxVal[ind] = B3DMAX(maxVal[ind], value);
            numValues[ind]++;
          }
        }
      }
    }
  }
  /* printf("numval %d %d %d %f %f %f\n", numValues[0], numValues[1], numValues[2], 
     maxVal[0], maxVal[1], maxVal[2]); */

  /* Get the value ID's if any and convert an entered ID to a type #.  Also set format
     for output based on maximum value */
  maxValType = -1;
  numCols = 0;
  for (ind = 0; ind < MAX_VALUE_COLS; ind++) {
    strcpy(format[ind], "%10.4f");
    if (numValues[ind]) {
      numCols++;
      maxValType = ind;
      if (maxVal[ind] > 10.1)
        strcpy(format[ind], "%10.2f");
      else if (maxVal[ind] > 1.01)
        strcpy(format[ind], "%10.3f");
    }
  }
  for (ind = 0; ind < B3DMIN(maxValType + 1, mod->obj[0].extra[IOBJ_EXSIZE - 1]); ind++)
    valueIDs[ind] = mod->obj[0].extra[IOBJ_EXSIZE - 2 - ind];
  /*printf("# vla ID %d %d %d %d\n", maxValType, valueIDs[0],valueIDs[1],valueIDs[2]);*/

  /* Error checks on the column for value 1 */
  if (colForVal1 < 1 || colForVal1 > numCols) {
      if (numCols)
        exitError("The column for general value type must be between 1 and %d", numCols);
      exitError("There are no general values stored in the model");
  }
  colForVal1--;
  
  /* Set up map from column to value type index */
  if (numValues[0])
    colToTypeMap[colForVal1] = 0;
  col = 0;
  for (ind = 1; ind <= maxValType; ind++) {
    if (numValues[ind]) {
      if (!colToTypeMap[col])
        col++;
      colToTypeMap[col++] = ind;
    }
  }
  /*printf("%d %d %d %d\n", numCols, colToTypeMap[0], colToTypeMap[1], colToTypeMap[2]);*/

  /* Output the header line */
  fprintf(fout, "%d   edited positions", npatch);
  if (mod->obj[0].extra[IOBJ_EXSIZE - 1] > 0)
    for (col = 0; col < numCols; col++)
      fprintf(fout, "  %d", valueIDs[colToTypeMap[col]]);
  fprintf(fout, "\n");

  for (ob = 0; ob < mod->objsize; ob++) {
    listInd = 0;
    for (co = 0; co < mod->obj[ob].contsize; co++) {
      listStart = listInd;
      if (mod->obj[ob].cont[co].psize >= 2) {
        pts = mod->obj[ob].cont[co].pts;
        ix = pts[0].x + 0.5;
        iy = pts[0].y + 0.5;
        iz = pts[0].z + 0.5;
        dx = (pts[1].x - pts[0].x) / mod->pixsize;
        dy = (pts[1].y - pts[0].y) / mod->pixsize;
        dz = (pts[1].z - pts[0].z) / mod->pixsize;
        if (mod->flags & IMODF_FLIPYZ)
          fprintf(fout, "%6d %5d %5d %8.2f %8.2f %8.2f", 
                  ix, iz, iy, dx, dz, dy);
        else
          fprintf(fout, "%6d %5d %5d %8.2f %8.2f %8.2f", 
                  ix, iy, iz, dx, dy, dz);
        for (col = 0; col < numCols; col++) {
          ind = colToTypeMap[col];
          if (numValues[ind]) {
            listInd = listStart;
            value = 0.;
            istoreFindValue(mod->obj[ob].store, co, GEN_STORE_VALUE1 + 2 * ind, &value,
                            &listInd);
            fprintf(fout, format[ind], value);
          }
        }
        fprintf(fout, "\n");

      }
    }
  }
  fclose(fout);
  exit(0);
}


/* This is meant to be called sequentially for all indexes in the entity,
   not for random access */
static int istoreFindValue(Ilist *list, int index, int type, float *value,
                     int *listInd)
{
  Istore *store;
  while (*listInd < ilistSize(list)) {
    store = istoreItem(list, *listInd);
    if ((store->flags & GEN_STORE_NOINDEX) || store->index.i > index)
      break;
    (*listInd)++;

    if (store->index.i == index && store->type == type) {
      *value = store->value.f;
      return 1;
    }
  }
  return 0;
}
