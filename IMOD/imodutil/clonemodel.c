/*
 *  clonemodel.c - create a new model containing multiple copies of an input
 *                 model at specified points / orientations. 
 *
 *  Author: John Heumann   email: heumannj@colorado.edu
 *
 *  Copyright (C) 2011-2019 by the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 * $Id$
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "imodel.h"
#include "b3dutil.h"
#include "iimage.h"
#include "parse_params.h"

/* 
 * Main entry
 */
#define lineSz 1024
int main( int argc, char *argv[])
{
  Imod *inModel, *outModel, *tmpModel;
  FILE *coordFP, *outFP, *refFP;
  MrcHeader hdata;
  char line[lineSz];
  char msg[256];
  char *progname = imodProgName(argv[0]);
  char *inFile, *outFile, *coordFile, *listString, *referenceFile = NULL;
  int iOutObj = 0, numOptArgs, numNonOptArgs, numRead, temp;
  int contour, *contourList = NULL, numContours = 0;
  long posn;
  int  valuesPresent = 0, vRangeSpecified = 0;
  float xMin, xMax, yMin, yMax, zMin, zMax, value, vMin, vMax;
  float x, y, z, xAngle, yAngle, zAngle;
  int maxCtrX = 0, maxCtrY = 0, maxCtrZ = 0;
  Imat *xform;
  unsigned char cmap[3][256];

  int numOptions = 10;
  const char *options[] = {
    "input:InputFile:FN:", "output:OutputFile:FN:", "at:AtPoints:FN:",
    "reference:ReferenceImage:FN:", "x:XRange:IP:", "y:YRange:IP:", "z:ZRange:IP:",
    "v:VRange:FP:", "help:usage:B:", "contours:ContourNumbers:LI:"};
  const char *usageString = 
    "Usage: clonemodel [options] -at locationFile inputModel outputModel";

  /* Parse parameters */
  PipSetUsageString(usageString);
  PipReadOrParseOptions(argc, argv, options, numOptions, progname, 3, 1,
                        1, &numOptArgs, &numNonOptArgs, 
                        imodUsageHeader);
  if (!PipGetBoolean("usage", &temp)) {
    PipPrintHelp(progname, 0, 1, 1);
    exit(0);
  }
  if (PipGetInOutFile((char *)"InputFile", 0, &inFile))
    exitError("No input file specified");
  if (PipGetString("AtPoints", &coordFile))
    exitError("No location/orientation file specified");
  if (PipGetInOutFile("OutputFile", 1, &outFile))
    exitError("No output file specified");
  if (PipGetTwoFloats("XRange", &xMin, &xMax)) {
    xMin = 0.0F;
    xMax = FLT_MAX;
  }
  if (PipGetTwoFloats("YRange", &yMin, &yMax)) {
    yMin = 0.0F;
    yMax = FLT_MAX;
  }
  if (PipGetTwoFloats("ZRange", &zMin, &zMax)) {
    zMin = 0.0F;
    zMax = FLT_MAX;
  }
  if (PipGetTwoFloats("VRange", &vMin, &vMax)) {
    vMin = FLT_MAX;
    vMax = -FLT_MAX;
    vRangeSpecified = 0;
  }
  else 
    vRangeSpecified = 1;
  if (!PipGetString("ContourNumbers", &listString)) {
    contourList = parselist(listString, &numContours);
    free(listString);
    if (!contourList)
      exitError("Bad entry in list of contour numbers");
  }
  if (!PipGetString("ReferenceImage", &referenceFile)) {
    refFP = iiFOpen(referenceFile, "r");
    if (!refFP)
      exitError("Opening reference image file %s", referenceFile);
    if (mrc_head_read(refFP, &hdata))
      exitError("Reading header for reference image file %s", referenceFile);
    free(referenceFile);
    iiFClose(refFP);
  }
  PipDone();

  /* Retrieve 3dmods standard pseudo-color map */
  cmapConvertRamp(cmapStandardRamp(), cmap);

  /* Open the csv location/orientation file and skip the header line */
  coordFP = fopen(coordFile, "r");
  if (!coordFP) {
    sprintf(msg, "Error opening location/orientation file:\n%s", coordFile);
    exitError(msg);
  }
  free(coordFile);
  if (fgets(line, lineSz, coordFP) == NULL) {
    sprintf(msg, "Error reading location/orientation file :\n%s", coordFile);
    exitError(msg);
  }

  /* Read the input model (to be cloned) */
  inModel = imodRead(inFile);
  if (!inModel) {
    sprintf(msg, "Error reading input model:\n%s", inFile);   
    exitError(msg);
    free(inFile);
  }

  /* Create new output and temporary models */
  outModel = imodNew();
  tmpModel = imodNew();
  /* Set tmpModel max coords. Rotation center will be midpoint */
  tmpModel->xmax = inModel->xmax;
  tmpModel->ymax = inModel->ymax;
  tmpModel->zmax = inModel->zmax;

   /* DNM: Attach reference image information and copy the pixel size from input model */
  if (referenceFile && imodSetRefImage(outModel, &hdata))
      exitError("Allocating a IrefImage structure for reference information");
  outModel->pixsize = inModel->pixsize;
  outModel->units = inModel->units;
   
  /* Make a pass through the summary file checking for the optional value *
     field and remembering min/max values if range was not specified      */
  posn = ftell(coordFP);
  while (fgets(line, lineSz, coordFP)) {
    numRead = sscanf(line, "%d,%g,%g,%g,%g,%g,%g,%g", &contour, &x, &y, &z, 
		     &xAngle, &yAngle, &zAngle, &value);
    if (numRead != 7 && numRead != 8)
      exitError("Error parsing summary file");
    if (numRead == 8) {
      valuesPresent = 1;
      if (!vRangeSpecified) {
        if (value < vMin) vMin = value;
        if (value > vMax) vMax = value;
      }
    }
  }
  if (fseek(coordFP, posn, SEEK_SET) == -1)
    exitError("Error repositioning summary file!");

  /* Now really process each point in the summary file  */
  xform = imodMatNew(3);
  while (fgets(line, lineSz, coordFP)) {
    int contourOk, inRange;

    numRead = sscanf(line, "%d,%g,%g,%g,%g,%g,%g,%g", &contour, &x, &y, &z, 
		     &xAngle, &yAngle, &zAngle, &value);
    if (numRead != 7 && numRead != 8)
      exitError("Error parsing summary file");
    contourOk = (numContours == 0 || 
                 numberInList(contour, contourList, numContours, 0));
    inRange = (x >= xMin && x <= xMax && y >= yMin && y <= yMax && 
               z >= zMin && z <= zMax);
    if (contourOk && inRange) {
      Iobj *obj = imodObjectGetFirst(inModel);
      Ipoint newCenter;
      int i = 0;

      /* Remember max coords seen */
      maxCtrX = B3DMAX(maxCtrX, x);
      maxCtrY = B3DMAX(maxCtrY, y);
      maxCtrZ = B3DMAX(maxCtrZ, z);

      /* Copy all the objects from the input to the temp model */        
      while (obj != NULL) {
        Iobj *tmpobj = imodObjectDup(obj);
        /* Assign psuedo-colors if optional value field is present */
	if (valuesPresent) {
          int index = (int) (255.0 * (value - vMin) / (vMax - vMin));
          imodObjectSetColor(tmpobj, (float)cmap[0][index] / 255.0f, 
 			             (float)cmap[1][index] / 255.0f, 
			             (float)cmap[2][index] / 255.0f);
	}
        imodNewObject(tmpModel);
        imodObjectCopy(tmpobj, &(tmpModel->obj[i++]));
        free(tmpobj);
        obj = imodObjectGetNext(inModel);
      }

      newCenter.x = x;
      newCenter.y = y;
      newCenter.z = z;
      /* Construct the rotation matrix for this point */
      imodMatId(xform);
      imodMatRot(xform, zAngle, b3dZ);
      imodMatRot(xform, yAngle, b3dY);
      imodMatRot(xform, xAngle, b3dX);

      /* Transform the temp model */
      imodTransModel3D(tmpModel, xform, NULL, newCenter, 1.0, 0);

      /* Copy the transformed objects to the output model */
      obj = imodObjectGetFirst(tmpModel);
      while (obj != NULL) {
        imodNewObject(outModel);
        imodObjectCopy(obj, &outModel->obj[iOutObj++]);
        obj = imodObjectGetNext(tmpModel);
      }
      /* Prepare temp model for reuse */
      free(tmpModel->obj);
      tmpModel->obj = NULL;
      tmpModel->objsize = 0;
    }
   
  }
  fclose(coordFP);
  imodMatDelete(xform);
  if (contourList)
    free(contourList);

  /*
   * Put some max coords in the output model for IMOD's benefit.
   * These may be smaller than the actual volume limits, but will
   * get corrected automatically if the model is opened on a volume 
   */
  outModel->xmax = maxCtrX + inModel->xmax / 2;
  outModel->ymax = maxCtrY + inModel->ymax / 2;
  outModel->zmax = maxCtrZ + inModel->zmax / 2;

  /* DNM: set them bigger based on the reference volume */
  if (referenceFile) {
    ACCUM_MAX(outModel->xmax, hdata.nx);
    ACCUM_MAX(outModel->ymax, hdata.ny);
    ACCUM_MAX(outModel->zmax, hdata.nz);
  }

  /* Write the output model */
  outFP = fopen(outFile, "w");
  if (!outFP) {
    sprintf(msg, "Error writing output file:\n%s", outFile);
    exitError(msg);
  }
  imodWrite(outModel, outFP);
  fclose(outFP);

  exit(0);
}
