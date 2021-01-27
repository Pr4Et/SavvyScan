/*
 *  adocxmlconv - Interconverts autodoc and XML files
 *
 *  Author: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 1995-2006 by the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 * $Id$
 */

#include <string.h>
#include "b3dutil.h"
#include "parse_params.h"
#include "autodoc.h"

int main( int argc, char *argv[] )
{
  int adocInd, wasXML, err, sectNotElem, sectNoName, childNotElem, childAttribs;
  int valueNotText, multipleChilds;
  int argInd = 1;
  char *rootElement = NULL;
  char *inFile, *outFile;
  char *progname = imodProgName(argv[0]);
  setStandardExitPrefix(progname);

  if (argc >= 2 && !strcmp(argv[1], "-r")) {
    argInd += 2;
    rootElement = argv[2];
  }
  if (argc < argInd + 2) {
    imodVersion(progname);
    imodCopyright();
    printf("Usage:  adocxml [-r name] inputFile outputFile\n"
           "   Converts autodoc to XML or conforming parts of XML file to autodoc\n"
           "   Use -r to supply root element name when converting to XML\n");
    exit(0);
  }
  inFile = argv[argInd];
  outFile = argv[argInd + 1];

  adocInd = AdocRead(inFile);
  if (adocInd < 0)
    exitError("Reading input file %s  (return value %d)\n", inFile, adocInd);
  wasXML = AdocXmlReadStatus(&sectNotElem, &sectNoName, &childNotElem, &childAttribs,
                             &valueNotText, &multipleChilds);
  if (!wasXML) {
    AdocSetWriteAsXML(1);
    if (rootElement)
      AdocSetXmlRootElement(rootElement);
  }
  AdocGetXmlRootElement(&rootElement);
  err = AdocWrite(outFile);
  if (err > 0)
    exitError("Writing output file %s as %s\n", outFile, wasXML ? "autodoc" : "XML");
  if (err < 0)
    printf("WARNING: %s - Failed to make backup file out of existing %s\n", progname,
           outFile);
  if (wasXML < 0) {
    printf("Converted to autodoc from XML file with some unsupported features:\n");
    if (sectNotElem)
      printf("  Top-level child was not an element (%d times)\n", sectNotElem);
    if (sectNoName)
      printf("  Top-level child did not have a name attribute (%d times)\n", sectNoName);
    if (childNotElem)
      printf("  Child of top-level child was not an element (%d times)\n", childNotElem);
    if (childAttribs)
      printf("  Child of top-level child had attributes (%d times)\n", childAttribs);
    if (valueNotText)
      printf("  Child of top-level child had non-text value (%d times)\n", valueNotText);
    if (multipleChilds)
      printf("  Child of top-level child had multiple children (%d times)\n", 
             multipleChilds);

  } else if (wasXML) {
    
    printf("Converted XML file with root element \"%s\" to autodoc\n", rootElement);
  } else {
    printf("Converted autodoc to XML file\n");
  }
  AdocDone();
  free(rootElement);
  exit(0);
}
