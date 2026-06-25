/*
 * autodoc.c - A parser and manager for autodoc files
 *
 * Copyright (C) 2006-2016 by the Regents of the University of 
 * Colorado.  See dist/COPYRIGHT for full notice.
 *
 * $Id$
 */                                                                           

#include "autodoc.h"
#include "parse_params.h"
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/types.h>
#ifndef _WIN32
#include <unistd.h>
#endif
#include "b3dutil.h"
#include "mxml.h"
#include "mxmlwrap.h"

/* A section has a name and holds the key-value pairs */
typedef struct adoc_section {
  char *name;          /* value after delimiter in section header */
  char **keys;         /* Array of strings with keys */
  char **values;       /* Array of strings with values */
  int numKeys;         /* Number of key/value pairs */
  int maxKeys;         /* Current size of array */
  char **comments;     /* List of comments strings */
  int *comIndex;       /* Array of key indexes they occur before */
  int numComments;     /* Number of comments */
  unsigned char *types;  /* Array of types of keys */
} AdocSection;

/* A collection holds sections of the same type and has a name */
typedef struct adoc_collection {
  char *name;             /* section type, before delimiter in header */
  AdocSection *sections;  /* Array of sections */
  int numSections;        /* Number of sections */
  int maxSections;        /* Current size of array */
} AdocCollection;

/* An autodoc has an array of collections and a list of sections in order */
typedef struct adoc_autodoc {
  AdocCollection *collections;  
  int numCollections;
  char **finalComments;
  int numFinalCom;
  int *collList;
  int *sectList;
  int numSections;
  int maxSections;
  int inUse;
  int backedUp;
  int writeAsXML;
  char *rootElement;    /* Name for root element for XML file, in or out */
} Autodoc;

/* The static variables that can hold multiple autodocs */
static Autodoc *sAutodocs = NULL;
static int sNumAutodocs = 0;
static int sCurAdocInd = -1;
static Autodoc *sCurAdoc = NULL;
static char *sNameForOrdering = NULL;

/* Static variables for XML writing */
static int sLastWasXML = 0;
static int sNumSectNotElem, sNumSectNoName, sNumChildNotElem, sNumChildAttribs;
static int sNumValueNotText, sNumMultipleChilds;

/* Static variables for writing to file or string */
static FILE *sFile;
static char *sString;
static int sBytesLeft;
static int sBytesWritten;

#define GLOBAL_NAME ADOC_GLOBAL_NAME
#define OPEN_DELIM "["
#define CLOSE_DELIM "]"
#define XML_START "<?xml"
#define XML_COMMENT_START "!--"
#define XML_COMSTART_LEN  3

static char sDefaultDelim[] = "=";
static char *sValueDelim = sDefaultDelim;
static char *sNewDelim = NULL;

#define BIG_STR_SIZE  10240
#define ERR_STR_SIZE  1024
#define MALLOC_CHUNK 10

static void deleteAdoc(Autodoc *adoc);
static int parseKeyValue(char *line, char *end, char **key, char **value);
static int lookupKey(AdocSection *sect, const char *key);
static int lookupCollection(Autodoc *adoc, const char *name);
static int addKey(AdocSection *sect, const char *key, const char *value, int type);
static int addSection(Autodoc *adoc, int collInd, const char *name);
static int addCollection(Autodoc *adoc, const char *name);
static int addAutodoc();
static AdocSection *getSection(const char *collName, int sectInd);
static int addComments(AdocSection *sect, char **comments, int *numComments,
                       int index);
static int writeFile(FILE *afile, char *string, int stringSize, int writeAll);
static int fsPrintf(const char *format, ...);
static int findSectionInAdocList(int collInd, int sectInd);
static int setKeyValueType(const char *typeName, int sectInd, const char *key, 
                           const char *value, int type);
static int sectSetKeyValueType(AdocSection *sect, const char *key, const char *value,
                               int type, int *keyInd);
static int setArrayOfValues(const char *typeName, int sectInd, const char *key, 
                            void *vals, int numVals, int valType);
static int readXmlFile(FILE *fp);
static int writeXmlFile(const char *filename);
static int *setupSectionOrder();
static int addToCommentList(char ***commentList, int *numComments, int *maxComments,
                            char *comment);
static int testAndAddComment(mxml_node_t *node, char ***commentList, int *numComments, 
                             int *maxComments);
static void handleFinalComments(int err, char ***commentList, int numComments);
static void writeCommentToXML(mxml_node_t *parent, char *comment);
static int adocMemoryError(void *ptr, const char *routine);

/*!
 * Reads an autodoc or autodoc-compatible XML file from the file specified by [filename],
 * and returns the index of the new autodoc (numbered from zero) or -1 for an error.
 */
int AdocRead(const char *filename)
{
  int gotSection = 0;
  int err, lineLen, indst, icol, ikey, lastInd, index;
  int badline = 1234;
  AdocSection *curSect;
  AdocCollection *coll;
  char *line, *lineEnd, *key, *value;
  char bigStr[BIG_STR_SIZE];
  char errStr[ERR_STR_SIZE];
  char commentChar = '#';
  char **commentList = NULL;
  int maxComments = 0;
  int numComments = 0;
  int firstLine = 1;
  FILE *afile;

  sValueDelim = sDefaultDelim;
  afile = fopen(filename, "r");
  if (!afile) {
    b3dError(stderr, "ERROR: AdocRead - Error opening autodoc file %s", filename);
    return -1;
  }

  /* Create a new adoc, which sets up global collection/section
     and takes care of cleanup if it fails */
  if ((index = addAutodoc()) < 0) {
    fclose(afile);
    return -1;
  }

  AdocSetCurrent(index);
  curSect = &sCurAdoc->collections[0].sections[0];
  lastInd = -1;
  sLastWasXML = 0;

  while (1) {

    /* We cannot allow in-line comments so that value lines can contain 
       anything.  But do allow blank and comment lines */
    lineLen = PipReadNextLine(afile, bigStr, BIG_STR_SIZE, commentChar, 1, 0,
                              &indst);
    /*puts(bigStr);
      printf("%d %d\n", lineLen, indst); */
    if (lineLen == -3)
      break;
    if (lineLen < 0) {
      b3dError(stderr, "ERROR: AdocRead - %s %s\n", lineLen == -2 ?
               "Error reading autodoc file" : "Line too long in autodoc file", filename);
      err = -1;
      break;
    }

    /* For first line, check for XML file and read that */
    line = &bigStr[indst];
    if (firstLine && PipStartsWith(line, XML_START)) {
      err = readXmlFile(afile);
      fclose(afile);
      if (err < 0)
        return err;
      return sCurAdocInd;
    }
    firstLine = 0;

    /* First check for comment and add it to list */
    if (indst >= lineLen || line[0] == commentChar) {
      if ((err = addToCommentList(&commentList, &numComments, &maxComments, bigStr)) != 0)
        break;
      continue;
    }

    if (PipStartsWith(line, OPEN_DELIM) && strstr(line, CLOSE_DELIM)) {

      /* If this is a section start, get name - value.  Here there must be
         a value and it is an error if there is none. */
      lineEnd = strstr(line, CLOSE_DELIM);
      err = parseKeyValue(line + strlen(OPEN_DELIM), lineEnd, &key, &value);
      if (!value)
        err = 1;
      if (err) {
        err = err > 0 ? badline : err;
        break;
      }

      /* Lookup the collection under the key and create one if not found */
      icol = lookupCollection(sCurAdoc, key);
      if (icol < 0) {
        if ((err = addCollection(sCurAdoc, key)))
          break;
        icol = sCurAdoc->numCollections - 1;
      }
      coll = &sCurAdoc->collections[icol];

      /* Add a section to the collection and set it as current one */
      if ((err = addSection(sCurAdoc, icol, value)))
        break;
      curSect = &coll->sections[coll->numSections - 1];
      gotSection = 1;
      free(key);
      free(value);
      lastInd = -1;

    } else {

      /* Otherwise this is key-value inside a section.  First check for
         continuation line and append to last value. */
      if (lastInd >= 0 && !strstr(line, sValueDelim) && 
          curSect->values[lastInd]) {
        ikey = strlen(curSect->values[lastInd]);
        curSect->values[lastInd] = (char *)realloc
          (curSect->values[lastInd], ikey + lineLen - indst + 3);
        if ((err = adocMemoryError(curSect->values[lastInd], "AdocRead")))
          break;

        /* Replace null with space and new null, then append new string */
        curSect->values[lastInd][ikey] = ' ';
        curSect->values[lastInd][ikey + 1] = 0x00;
        strcat(curSect->values[lastInd], line);
        continue;
      }

      /* This should be a key-value pair now */
      lineEnd = line + lineLen - indst;
      if ((err = parseKeyValue(line, lineEnd, &key, &value))) {
        err = err > 0 ? badline : err;
        break;
      }

      /* Handle new key-value delimiter - replace previous new value if any */
      if (!gotSection && !strcmp(key, "KeyValueDelimiter") && value) {
        if (sNewDelim)
          free(sNewDelim);
        sNewDelim = strdup(value);
        if ((err = adocMemoryError(sNewDelim, "AdocRead")))
          break;
        sValueDelim = sNewDelim;
      }

      /* Handle change of comment character */
      if (!gotSection && !strcmp(key, "CommentCharacter"))
        commentChar = value[0];

      /* Look up the key first to replace an existing value */
      ikey = lookupKey(curSect, key);
      if (ikey >= 0) {
        if (curSect->values[ikey])
          free(curSect->values[ikey]);
        curSect->values[ikey] = value;
        free(key);
        lastInd = ikey;

      } else {

        /* Or just add the key-value */
        if ((err = addKey(curSect, key, value, ADOC_STRING)))
          break;
        free(key);
        if (value)
          free(value);
        lastInd = curSect->numKeys -1;
      }
    }

    /* If there are comments, attach to item just added */
    if (numComments)
      if ((err = addComments(curSect, commentList, &numComments, lastInd)))
      break;
  }

  /* END OF FILE: If error, clean out autodoc, compose message for bad line */
  if (err) {
    deleteAdoc(sCurAdoc);
    if (err == badline) {
      bigStr[ERR_STR_SIZE - 50] = 0x00;
      b3dError(stderr, "Error: AdocRead -Improperly formatted line in autodoc: %s\n",
               bigStr);
      err = -1;
    }
  }

  handleFinalComments(err, &commentList, numComments);

  fclose(afile);
  return (err ? err : index);
}

/*!
 * Returns whether the last file read in was an XML file, and information about aspects 
 * that could not be encoded as an autodoc if the last file was an XML file.  Each value 
 * is the number of occurrences of the problem: ^
 * [sectNotElem] if a top-level child is not an element, ^
 * [sectNoName] if a top-level child does not have a name attribute, ^
 * [childNotElem] if a child of a top-level child is not an element, ^
 * [childAttribs] if a child has attributes, ^
 * [valueNotText] if a child's value is not text, ^
 * [multipleChilds] if a child of a top-level child has multiple children. ^
 * If a top-level child has attributes other than "name", these will become
 * key-value pairs indistinguishable from the ones coming from child nodes, and there is
 * no indicator that this occurred.  The return value is 0 if file was not XML, 1 if it 
 * was and had no problems, or -1 if it was XML and had some of these problems.
 */
int AdocXmlReadStatus(int *sectNotElem, int *sectNoName, int *childNotElem,
                      int *childAttribs, int *valueNotText, int *multipleChilds)
{
  *sectNoName = sNumSectNoName;
  *sectNotElem = sNumSectNotElem;
  *childNotElem = sNumChildNotElem;
  *childAttribs = sNumChildAttribs;
  *valueNotText = sNumValueNotText;
  *multipleChilds = sNumMultipleChilds;
  if (!sLastWasXML)
    return 0;
  return ((sNumSectNotElem + sNumSectNoName + sNumChildNotElem + sNumChildAttribs +
           sNumValueNotText + sNumMultipleChilds) > 0 ? -1 : 1);
}

/*!
 * Try to read an image metadata file and determine its properties.  It will 
 * look for and try to read the file with name [filename] if [addMdoc] is 0,
 * otherwise it will append '.mdoc' to [filename].  Returns a non-zero in
 * [montage] if the file indicates a montage; the number of sections in
 * [numSect], and 1 for a metadata autodoc or 2 for an image series autodoc
 * in [sectType].  The return value is the index of the autodoc, or -1 for an
 * error opening or reading the file, -2 if the file does not exist, or -3 if
 * it is neither type of metadata file.
 */
int AdocOpenImageMetadata(const char *filename, int addMdoc, int *montage,
                          int *numSect, int *sectType)
{
  struct stat buf;
  char *usename = (char *)filename;
  int series, index;

  /* Attach extension to file if requested */
  if (addMdoc) {
    usename = (char *)malloc(strlen(filename) + 6);
    if (!usename)
      return -1;
    sprintf(usename, "%s.mdoc", filename);
  }

  /* Return -2 if it does not exist, -1 if error reading it */
  if (stat(usename, &buf)) {
    index = -2;
  } else {
    index = AdocRead(usename);
  }
  if (addMdoc)
    free(usename);
  if (index < 0)
    return index;

  series = AdocGetImageMetaInfo(montage, numSect, sectType);
  if (series) {
    AdocClear(index);
    return series;
  }
  return index;
}

/*!
 * Determines whether the current autodoc is a valid image metadata autodoc and
 * returns a non-zero in [montage] if the file indicates a montage; the number of 
 * sections in [numSect]; and 1 for a metadata autodoc, 2 for an image series autodoc,
 * or 3 for any other autodoc with sections named ZValue (including the autodoc holding
 * attributes for an HDF file).  The return value is -3 if it 
 * is not one of these types of autodocs.
 */
int AdocGetImageMetaInfo(int *montage, int *numSect, int *sectType) 
{
  int series;
  char *usename;

  if (!AdocGetString(GLOBAL_NAME, 0, "ImageFile", &usename)) {
    *sectType = 1;
    free(usename);
    *numSect = AdocGetNumberOfSections(ADOC_ZVALUE_NAME);
  } else if (!AdocGetInteger(GLOBAL_NAME, 0, "ImageSeries", &series) && series) {
    *sectType = 2;
    *numSect = AdocGetNumberOfSections("Image");
  } else {
    *numSect = AdocGetNumberOfSections(ADOC_ZVALUE_NAME);
    if (! *numSect)
      return -3;
    *sectType = 3;
  }
 
  *montage = 0;
  if (AdocGetInteger(GLOBAL_NAME, 0, "Montage", montage))
    AdocGetInteger(GLOBAL_NAME, 0, "IMOD.Montage", montage);
  return 0;
}

/*!
 * Finds an unused autodoc in the internal array if one exists or creates a new one if
 * not, makes it the current one, and returns its index or -1 for an error.
 */
int AdocNew()
{
  int err;
  if ((err = addAutodoc()) < 0)
    return err;
  AdocSetCurrent(err);
  return err;
}  

/*!
 * Returns the current autodoc index, or -1 if there is none
 */
int AdocGetCurrentIndex()
{
  return sCurAdocInd;
}

/*!
 * Makes the autodoc at [index] be the current autodoc.  The index is numbered
 * from zero.  Returns -1 for an index out of bounds.
 */
int AdocSetCurrent(int index)
{
  if (index < 0 || index >= sNumAutodocs)
    return -1;
  sCurAdocInd = index;
  sCurAdoc = &sAutodocs[sCurAdocInd];
  return 0;
}

/*!
 * Deletes all data from the autodoc at [index] and marks it as unused.
 */
void AdocClear(int index)
{
  if (index >= 0 && index < sNumAutodocs)
    deleteAdoc(&sAutodocs[index]);
}

/*!
 * Deletes all autodocs and returns the module to its initial state. 
 */
void AdocDone()
{
  int i;
  for (i = 0; i < sNumAutodocs; i++)
    deleteAdoc(&sAutodocs[i]);
  if (sAutodocs)
    free(sAutodocs);
  B3DFREE(sNameForOrdering);
  sAutodocs = NULL;
  sNumAutodocs = 0;
  sCurAdocInd = -1;
  sCurAdoc = NULL;
}

/*!
 * Writes the current autodoc to the file specified by [filename].  Returns 1
 * for failure to back up a previous file, and -1 for other errors 
 */
int AdocWrite(const char *filename)
{
  int i,backerr, retval = 0;
  FILE *afile;

  if (!sCurAdoc)
    return -1;
  if (!sCurAdoc->backedUp)
    backerr = imodBackupFile(filename);
  sCurAdoc->backedUp = 1;
  if (sCurAdoc->writeAsXML)
    return writeXmlFile(filename);
  afile = fopen(filename, "w");
  if (!afile)
    return -1;
  if (writeFile(afile, NULL, 0, 1))
    retval = -1;
  else
    for (i = 0; i < sCurAdoc->numFinalCom; i++)
      fprintf(afile, "%s\n", sCurAdoc->finalComments[i]);

  fclose(afile);
  return retval ? retval : backerr;
}

/*!
 * Sets the current autodoc to be written as XML if [asXML] is non-zero.
 */
void AdocSetWriteAsXML(int asXML)
{
  if (sCurAdoc)
    sCurAdoc->writeAsXML = asXML ? 1 : 0;
}

/*!
 * Returns 1 if current autodoc is to be written as XML, 0 if not, -1 for invalid 
 * current autodoc
 */
int AdocGetWriteAsXML()
{
  if (!sCurAdoc)
    return -1;
  return sCurAdoc->writeAsXML;
}

/*!
 * Returns a copy of the root element name into [string] if a file was loaded from XML;
 * it should be freed with {free} if it is non-NULL.  Returns -1 for
 * no current autodoc or 1 for memory error.
 */
int AdocGetXmlRootElement(char **string)
{
  if (!sCurAdoc)
    return -1;
  *string = NULL;
  if (sCurAdoc->rootElement) {
    *string = strdup(sCurAdoc->rootElement);
    if (adocMemoryError(*string, "AdocGetXmlRootElement"))
      return 1;
  }
  return 0;
}

/*!
 * Sets the root element name to [element] when writing as an XML file; returns -1 for
 * no current autodoc or 1 for memory error.
 */
int AdocSetXmlRootElement(const char *element)
{
  if (!sCurAdoc || !element)
    return -1;
  B3DFREE(sCurAdoc->rootElement);
  sCurAdoc->rootElement = strdup(element);
  if (adocMemoryError(sCurAdoc->rootElement, "AdocSetXmlRootElement"))
    return 1;
  return 0;
}

/*!
 * Appends the last section in the current autodoc to the file specified by 
 * [filename]; or writes the autodoc completely if it is being written as XML.
 * Returns -1 for errors.
 */
int AdocAppendSection(const char *filename)
{
  FILE *afile;
  int retval;
  if (!sCurAdoc)
    return -1;
  if (sCurAdoc->writeAsXML)
    return writeXmlFile(filename);
  afile = fopen(filename, "a");
  if (!afile)
    return -1;
  retval = writeFile(afile, NULL, 0, 0);
  fclose(afile);
  return retval;
}

/*!
 * Writes either the entire current autodoc, if [writeAll] is non-zero, or just the last 
 * section otherwise, into [string], whose size must be specified by [stringSize].
 * Returns -1 for all errors, including not fitting into the string.
 */
int AdocPrintToString(char *string, int stringSize, int writeAll)
{
  return writeFile(NULL, string, stringSize, writeAll);
}

/*!
 * Sets the next write so that all sections in the collection with name [typeName]
 * are output after the global section and in order by the value of their names.
 * Returns 1 for memory error.
 */
int AdocOrderWriteByValue(const char *typeName)
{
  B3DFREE(sNameForOrdering);
  if (!typeName)
    return 0;
  sNameForOrdering = strdup(typeName);
  if (adocMemoryError(sNameForOrdering, "AdocOrderWriteByValue"))
    return 1;
  return 0;
}

/* Function to actually write the file or last section only */
static int writeFile(FILE *afile, char *string, int stringSize, int writeAll)
{
  int i,j,k, ind, comInd, write, lastBlank, useInd, retval = 0;
  AdocCollection *coll;
  AdocSection *sect;
  int *ordSectInds;
  int orderedWrite = writeAll && sNameForOrdering && sCurAdoc->numSections > 1;
  sFile = afile;
  sString = string;
  sBytesLeft = stringSize;
  sBytesWritten = 0;
  if (!afile && !string)
    return -1;

  /* For ordered writing, get arrays for indexes and float values, and set up values
     for all the sections of the given type */
  if (orderedWrite) {
    ordSectInds = setupSectionOrder();
    if (!ordSectInds)
      return -1;
  }

  /* Initialize delimiter, loop on indexes in the autodoc */
  sValueDelim = sDefaultDelim;
  for (ind = 0; ind < sCurAdoc->numSections; ind++) {
    write = (writeAll || ind == sCurAdoc->numSections - 1) ? 1 : 0;
    useInd = orderedWrite ? ordSectInds[ind] : ind;
    i = sCurAdoc->collList[useInd];
    j = sCurAdoc->sectList[useInd];
    coll = &sCurAdoc->collections[i];
    sect = &coll->sections[j];

    /* dump comments before section */
    comInd = 0;
    lastBlank = 0;
    while (write && comInd < sect->numComments && sect->comIndex[comInd] == -1) {
      if (write) {
        lastBlank = sect->comments[comInd][0] == 0x00 ? 1 : 0;
        if (fsPrintf("%s\n", sect->comments[comInd++])) {
          retval = -1;
          break;
        }
      }
    }
    if (retval)
      break;

    /* Write section name unless we're in global */
    if ((i || j || strcmp(sect->name, GLOBAL_NAME)) && write) {
      if (fsPrintf("%s[%s %s %s]\n", lastBlank ? "" : "\n", coll->name,
                   sValueDelim, sect->name)) {
        retval = -1;
        break;
      }
    }

    /* Loop on key-values */
    for (k = 0; k < sect->numKeys; k++) {

      /* dump comments associated with this index */
      while (write && comInd < sect->numComments && sect->comIndex[comInd] == k) {
        if (write) {
          if (fsPrintf("%s\n", sect->comments[comInd++])) {
            retval = -1;
            break;
          }
        }
      }
      if (retval)
        break;

      /* Print key-value pairs with non-null values */
      if (sect->keys[k] && sect->values[k]) {
        if (write) {
          if (fsPrintf("%s %s %s\n", sect->keys[k], sValueDelim, 
                       sect->values[k])) {
            retval = -1;
            break;
          }
        }
        
        /* After a new delimiter is written, need to set delimiter */
        if (!i && !j && !strcmp("KeyValueDelimiter", sect->keys[k])) {
          if (sNewDelim)
            free(sNewDelim);
          sNewDelim = strdup(sect->values[k]);
          if (adocMemoryError(sNewDelim, "AdocWrite")) {
            retval = -1;
            break;
          }
          sValueDelim = sNewDelim;
        }

        /* Print keys without values too */
      } else if (sect->keys[k] && write) {
        if (fsPrintf("%s %s \n", sect->keys[k], sValueDelim)) {
          retval = -1;
          break;
        }
      }
    }
    if (retval)
      break;
  }
  if (orderedWrite)
    free(ordSectInds);
  return retval;
}

/* Function to print to file or string using the static variables set up by writeFile */
static int fsPrintf(const char *format, ...)
{
  int retval = 0, numWritten;
  va_list args;
  va_start(args, format);
  if (sFile) {
    if (vfprintf(sFile, format, args) < 0) {
      b3dError(stderr, "ERROR: AdocWrite - writing element to file\n");
      retval = -1;
    }
  } else {
    numWritten = vsnprintf(&sString[sBytesWritten], sBytesLeft, format, args);
    if (numWritten >= sBytesLeft) {
      b3dError(stderr, "ERROR: AdocWrite - writing element to string\n");
      retval = -1;
    } else {
      sBytesLeft -= numWritten;
      sBytesWritten += numWritten;
    }
  }
  va_end(args);
  return retval;
}

/* 
 * Returns a set of indexes to the sections to be used when writing by ordered values.
 * Indexes should be freed when done.
 */
static int *setupSectionOrder()
{
  int i, j, ind;
  AdocCollection *coll;
  float maxValue = -1.e37;
  float *ordSectValues;
  int *ordSectInds;

  ordSectValues = B3DMALLOC(float, sCurAdoc->numSections);
  ordSectInds = B3DMALLOC(int, sCurAdoc->numSections);
  if (!ordSectInds || !ordSectValues) {
    adocMemoryError(NULL, "setupOrderedWrite");
    return NULL;
  }
  ordSectValues[0] = maxValue;
  ordSectInds[0] = 0;
  for (ind = 1; ind < sCurAdoc->numSections; ind++) {
    i = sCurAdoc->collList[ind];
    j = sCurAdoc->sectList[ind];
    ordSectInds[ind] = ind;
    coll = &sCurAdoc->collections[i];
    if (!strcmp(coll->name, sNameForOrdering) && coll->sections[j].name) {
      ordSectValues[ind] = atof(coll->sections[j].name);
      maxValue = B3DMAX(maxValue, ordSectValues[ind]);
    }
  }

  /* Then set up values for the rest of the sections, above the real values and in
     the order they occur */
  maxValue = B3DMAX(0., maxValue);
  for (ind = 1; ind < sCurAdoc->numSections; ind++) {
    i = sCurAdoc->collList[ind];
    coll = &sCurAdoc->collections[i];
    if (strcmp(coll->name, sNameForOrdering) || !coll->sections[j].name) {
      maxValue += 1.;
      ordSectValues[ind] = maxValue;
    }
  }

  /* Sort, use the sorted indexes below */
  rsSortIndexedFloats(ordSectValues, ordSectInds, sCurAdoc->numSections);
  free(ordSectValues);
  return ordSectInds;
}

/*!
 * Adds a section of type specified by [typeName] and name given by [name].  
 * Returns the index of the new section in the collection of sections of that
 * type, or -1 for error. 
 */
int AdocAddSection(const char *typeName, const char *name)
{
  AdocCollection *coll;
  int collInd;

  if (!sCurAdoc || !typeName || !name)
    return -1;
  collInd = lookupCollection(sCurAdoc, typeName);
  if (collInd < 0) {
    if (addCollection(sCurAdoc, typeName))
      return -1;
    collInd = sCurAdoc->numCollections - 1;
  }
  coll = &sCurAdoc->collections[collInd];
  if (addSection(sCurAdoc, collInd, name))
    return -1;
  return coll->numSections - 1;
}

/*!
 * Inserts a section with name given by [name] into the collection of sections of type 
 * [typeName], at the index [sectInd].  The collection must exist unless [sectInd] is 0,
 * and [sectInd] must be less than or equal to the number of sections in that collection.
 * Returns -1 for error.
 */
int AdocInsertSection(const char *typeName, int sectInd, const char *name)
{
  AdocCollection *coll;
  int i, collInd, masterInd, numSect = 0;
  AdocSection newSect;
  if (!sCurAdoc || !typeName || !name)
    return -1;
  collInd = lookupCollection(sCurAdoc, typeName);
  if (collInd >= 0)
    numSect = sCurAdoc->collections[collInd].numSections;
  if (sectInd < 0 || sectInd > numSect)
    return -1;
  
  /* Find the index of this section in the master list if it needs to be shuffled */
  if (sectInd < numSect) {
    masterInd = findSectionInAdocList(collInd, sectInd);
    if (masterInd < 0)
      return -1;
  }

  /* Add section to end regardless, then return if that is all that is needed */
  if (AdocAddSection(typeName, name) < 0)
    return -1;
  if (sectInd == numSect)
    return 0;

  /* Fix collection index if a new collection had to be added */
  if (collInd < 0)
    collInd = sCurAdoc->numCollections - 1;
  coll = &sCurAdoc->collections[collInd];

  /* Save the new section then move existing sections up and copy new one into place */
  memcpy(&newSect, &coll->sections[coll->numSections - 1], sizeof(AdocSection));
  for (i = coll->numSections - 1; i > sectInd; i--)
    memcpy(&coll->sections[i], &coll->sections[i - 1], sizeof(AdocSection));
  memcpy(&coll->sections[sectInd], &newSect, sizeof(AdocSection));
  
  /* Move the master lists up and decrement any other indices in this collection */
  for (i = sCurAdoc->numSections - 1; i > masterInd; i--) {
    sCurAdoc->collList[i] = sCurAdoc->collList[i - 1];
    sCurAdoc->sectList[i] = sCurAdoc->sectList[i - 1];
    if (sCurAdoc->collList[i] == collInd && sCurAdoc->sectList[i] >= sectInd)
      sCurAdoc->sectList[i]++;
  }

  return 0;
}

/*!
 * Deletes the section at index [sectInd] from the collection of sections of type 
 * [typeName].  Returns -1 for error.
 */
int AdocDeleteSection(const char *typeName, int sectInd)
{
  AdocCollection *coll;
  int collInd, i, masterInd;
  if (!sCurAdoc || !typeName)
    return -1;
  collInd = lookupCollection(sCurAdoc, typeName);
  if (collInd < 0)
    return -1;
  coll = &sCurAdoc->collections[collInd];
  if (sectInd < 0 || sectInd >= coll->numSections)
    return -1;

  /* Find the index of this section in the master list */
  masterInd = findSectionInAdocList(collInd, sectInd);
  if (masterInd < 0)
    return -1;

  /* Repack the sections */
  for (i = sectInd + 1; i < coll->numSections; i++)
    memcpy(&coll->sections[i - 1], &coll->sections[i], sizeof(AdocSection));
  coll->numSections--;

  /* Repack the master list and decrement any other indices in this collection */
  for (i = masterInd + 1; i < sCurAdoc->numSections; i++) {
    if (sCurAdoc->collList[i] == collInd && sCurAdoc->sectList[i] > sectInd)
      sCurAdoc->sectList[i]--;
    sCurAdoc->collList[i - 1] = sCurAdoc->collList[i];
    sCurAdoc->sectList[i - 1] = sCurAdoc->sectList[i];
  }
  sCurAdoc->numSections--;
  return 0;
}

/*!
 * Changes the name of the section at index [sectInd] from the collection of sections of 
 * type [typeName] to [newName].  Returns -1 for error.
 */
int AdocChangeSectionName(const char *typeName, int sectInd, const char *newName)
{
  AdocCollection *coll;
  int collInd;
  char *newCopy;
  if (!sCurAdoc || !typeName || !newName)
    return -1;
  collInd = lookupCollection(sCurAdoc, typeName);
  if (collInd < 0)
    return -1;
  coll = &sCurAdoc->collections[collInd];
  if (sectInd < 0 || sectInd >= coll->numSections)
    return -1;
  newCopy = strdup(newName);
  if (!newCopy)
    return -1;
  B3DFREE(coll->sections[sectInd].name);
  coll->sections[sectInd].name = newCopy;
  return 0;
}

/*!
 * Looks up a section of type specified by [typeName] and name given by [name].  
 * Returns the index of that section in the collection of sections of that
 * type, -1 if there is none, or -2 for error.
 */
int AdocLookupSection(const char *typeName, const char *name)
{
  AdocCollection *coll;
  int collInd, sectInd;

  if (!sCurAdoc || !typeName || !name)
    return -2;
  collInd = lookupCollection(sCurAdoc, typeName);
  if (collInd < 0)
    return -2;
  coll = &sCurAdoc->collections[collInd];
  for (sectInd = 0; sectInd < coll->numSections; sectInd++)
    if (!strcmp(coll->sections[sectInd].name, name))
      return sectInd;
  return -1;
}

/*!
 * Looks up a section of type specified by [typeName] whose name converts to the integer
 * value [nameValue]. Returns the index of that section in the collection of sections of
 * that type, -1 if there is none, or -2 for error.  When calling from Fortran, be sure
 * to pass the actual value in the file, as this routine will not adjust the passed
 * value down by 1.
 */
int AdocLookupByNameValue(const char *typeName, int nameValue)
{
  char buf[15];
  sprintf(buf, "%d", nameValue);
  return AdocLookupSection(typeName, buf);
}

/*!
 * Looks in the collection of sections of type [typeName], converts their name strings
 * to integers, and returns the index of the first section whose name is greater than
 * [nameValue], the number of sections if there is no such section, or -1 for error
 * (including if a section exists whose name converts to [nameValue]).  This returned
 * index can thus be used in a call to @AdocInsertSection to maintain the sections in
 * numeric order by name value.
 */
int AdocFindInsertIndex(const char *typeName, int nameValue)
{
  AdocCollection *coll;
  int collInd, sectInd, sectValue;

  if (!sCurAdoc || !typeName)
    return -1;
  collInd = lookupCollection(sCurAdoc, typeName);
  if (collInd < 0)
    return 0;
  coll = &sCurAdoc->collections[collInd];
  for (sectInd = 0; sectInd < coll->numSections; sectInd++) {
    sectValue = atoi(coll->sections[sectInd].name);
    if (nameValue == sectValue)
      return -1;
    if (nameValue < sectValue)
      return sectInd;
  }
  return coll->numSections;
}

/*!
 * Transfers all key/value pairs from the section of index [sectInd] in the collection of
 * type [typeName] of the current autodoc to a section of the same type in the autodoc
 * with index [toAdocInd], with a name [newName].  If [byValue] is non-zero, the new
 * name is converted to an integer and the section is inserted so as to maintain sections
 * in order.  The collection in the receiving autodoc is created if necessary; the new
 * section may already exist or will be created if necessary; existing values with the
 * same key are overwritten.  To transfer the global data, set [typeName] to 
 * ADOC_GLOBAL_NAME, [sectInd] to 0, and [byValue] to 0; [newName] is ignored.  Returns
 * -1 if the originating section does not exist, -2 for incorrect index for new autodoc
 * or [newName] NULL (unless [typeName] is ADOC_GLOBAL_NAME), -3 for failure to add a new 
 * section, or -4 for failure to add a key-value.
 */
int AdocTransferSection(const char *typeName, int sectInd, int toAdocInd, 
                        const char *newName, int byValue)
{
  int err = 0, ind, newSectInd, collInd, nameVal;
  AdocSection *sect;
  int curIndSave = sCurAdocInd;

  /* Get the section then switch adocs */
  if (!(sect = getSection(typeName, sectInd)))
    return -1;
  if (toAdocInd < 0 || toAdocInd == sCurAdocInd || toAdocInd >= sNumAutodocs)
    return -2;
  AdocSetCurrent(toAdocInd);

  /* Set index to 0 for global section or go on to look up and/or add section */
  if (!strcmp(typeName, GLOBAL_NAME)) {
    newSectInd = 0;
  } else {
    if (!newName)
      return -2;

    /* If the section does not exist, add it, using insert if the collection does exist */
    newSectInd = AdocLookupSection(typeName, newName);
    if (newSectInd < 0) {
      collInd = lookupCollection(sCurAdoc, typeName);
      if (collInd < 0) {
        newSectInd = 0;
        err = AdocAddSection(typeName, newName);
      } else {
        newSectInd = sCurAdoc->collections[collInd].numSections;
        if (byValue) {
          nameVal = atoi(newName);
          newSectInd = AdocFindInsertIndex(typeName, nameVal);
        }
        err = AdocInsertSection(typeName, newSectInd, newName);
      }
      if (err < 0) {
        AdocSetCurrent(curIndSave);
        return -3;
      }
    }
  }

  /* Copy the key/values and their types.  Skip NULL ones, which happen with HDF adoc */
  for (ind = 0; ind < sect->numKeys && !err; ind++) {
    if (sect->keys[ind] && sect->values[ind] && setKeyValueType
        (typeName, newSectInd, sect->keys[ind], sect->values[ind], sect->types[ind]) < 0)
      err = -4;
  }
  AdocSetCurrent(curIndSave);
  return err;
}

/*!
 * Sets a key-value pair to [key] and [value] in the section with index 
 * [sectInd] in the collection of sections of type [typeName].  The section 
 * must already exist.  Replaces an existing value if any.  [value] may be 
 * NULL.  Returns -1 for error.
 */
int AdocSetKeyValue(const char *typeName, int sectInd, const char *key, const char *value)
{
  return setKeyValueType(typeName, sectInd, key, value, ADOC_STRING);
}

/*
 * The function that actually sets a value, given the type as well.
 */
static int setKeyValueType(const char *typeName, int sectInd, const char *key, 
                           const char *value, int type)
{
  AdocSection *sect;
  int keyInd;

  if (!(sect = getSection(typeName, sectInd)))
    return -1;
  if (!key || !value)
    return -1;
  return sectSetKeyValueType(sect, key, value, type, &keyInd);
}

static int sectSetKeyValueType(AdocSection *sect, const char *key, const char *value,
                               int type, int *keyInd)
{
  *keyInd = lookupKey(sect, key);

  /* If key already exists, clear out value and set it again */
  if (*keyInd >= 0) {
    if (sect->values[*keyInd])
      free(sect->values[*keyInd]);
    if (value) {
      sect->values[*keyInd] = strdup(value);
      if (adocMemoryError(sect->values[*keyInd], "AdocSetKeyValue"))
        return -1;
      sect->types[*keyInd] = type;
    } else {
      sect->values[*keyInd] = NULL;
      sect->types[*keyInd] = ADOC_NO_VALUE;
    }
  } else {
    *keyInd = sect->numKeys;
    return addKey(sect, key, value, type);
  }
  return 0;
}


/*!
 * Sets the value of [key] to the integer [ival] in the section with index 
 * [sectInd] in the collection of sections of type [typeName].  The section 
 * must already exist.  Replaces an existing value if any.
 * Returns -1 for error.
 */
int AdocSetInteger(const char *typeName, int sectInd, const char *key, int ival)
{
  char str[30];
  sprintf(str, "%d", ival);
  return(setKeyValueType(typeName, sectInd, key, str, ADOC_ONE_INT));
}

/*!
 * Like @AdocSetInteger, except that the value is set to the two integers
 * [ival1] [ival2].
 */
int AdocSetTwoIntegers(const char *typeName, int sectInd, const char *key, int ival1, 
                       int ival2)
{
  char str[60];
  sprintf(str, "%d %d", ival1, ival2);
  return(setKeyValueType(typeName, sectInd, key, str, ADOC_TWO_INTS));
}

/*!
 * Like @AdocSetInteger, except that the value is set to the three integers
 * [ival1] [ival2] [ival3].
 */
int AdocSetThreeIntegers(const char *typeName, int sectInd, const char *key, int ival1,
                         int ival2, int ival3)
{
  char str[90];
  sprintf(str, "%d %d %d", ival1, ival2, ival3);
  return(setKeyValueType(typeName, sectInd, key, str, ADOC_THREE_INTS));
}

/*!
 * Like @AdocSetInteger, except that the value is set to the array of [numVals] integers 
 * in [ivals].
 */
int AdocSetIntegerArray(const char *typeName, int sectInd, const char *key, int *ivals,
                        int numVals)
{
  return(setArrayOfValues(typeName, sectInd, key, ivals, numVals, ADOC_INT_ARRAY));
}

/*!
 * Sets the value of [key] to the float [val] in the section with index 
 * [sectInd] in the collection of sections of type [typeName].  The section 
 * must already exist.  Replaces an existing value if any.
 * Returns -1 for error.
 */
int AdocSetFloat(const char *typeName, int sectInd, const char *key, float val)
{
  char str[30];
  sprintf(str, "%g", val);
  return(setKeyValueType(typeName, sectInd, key, str, ADOC_ONE_FLOAT));
}

/*!
 * Like @AdocSetFloat, except that the value is set to the two floats
 * [val1] [val2].
 */
int AdocSetTwoFloats(const char *typeName, int sectInd, const char *key, float val1, 
                     float val2)
{
  char str[60];
  sprintf(str, "%g %g", val1, val2);
  return(setKeyValueType(typeName, sectInd, key, str, ADOC_TWO_FLOATS));
}

/*!
 * Like @AdocSetFloat, except that the value is set to the three floats
 * [val1] [val2] [val3].
 */
int AdocSetThreeFloats(const char *typeName, int sectInd, const char *key, float val1,
                       float val2, float val3)
{
  char str[90];
  sprintf(str, "%g %g %g", val1, val2, val3);
  return(setKeyValueType(typeName, sectInd, key, str, ADOC_THREE_FLOATS));
}

/*!
 * Like @AdocSetFloat, except that the value is set to the array of [numVals] floats in 
 * [vals].
 */
int AdocSetFloatArray(const char *typeName, int sectInd, const char *key, float *vals,
                      int numVals)
{
  return(setArrayOfValues(typeName, sectInd, key, vals, numVals, ADOC_FLOAT_ARRAY));
}
    
/*! Like @AdocSetFloat, except that [val] is a double */
int AdocSetDouble(const char *typeName, int sectInd, const char *key, double val)
{
  char str[30];
  sprintf(str, "%g", val);
  return(setKeyValueType(typeName, sectInd, key, str, ADOC_ONE_DOUBLE));
}

/*
 * Function to set a line of floats or integers as the string value; valType should be
 * ADOC_INT_ARRAY or ADOC_FLOAT_ARRAY.
 */
static int setArrayOfValues(const char *typeName, int sectInd, const char *key, 
                            void *vals, int numVals, int valType)
{
  char tmp[40];
  char *fullStr;
  int *ivals = (int *)vals;
  float *fvals = (float *)vals;
  int ind, totLen = 0;

  /* Add up the characters needed for the eahc value */
  for (ind = 0; ind < numVals; ind++) {
    if (valType == ADOC_INT_ARRAY)
      sprintf(tmp, "%d ", ivals[ind]);
    else
      sprintf(tmp, "%g ", fvals[ind]);
    totLen += strlen(tmp) + 1;
  }

  /* Get the string and build it up by writing again */
  fullStr = B3DMALLOC(char, totLen);
  if (!fullStr)
    return -1;
  fullStr[0] = 0x00;
  for (ind = 0; ind < numVals; ind++) {
    if (valType == ADOC_INT_ARRAY)
      sprintf(tmp, "%s%d", ind ? " " : "", ivals[ind]);
    else
      sprintf(tmp, "%s%g", ind ? " " : "", fvals[ind]);
    strcat(fullStr, tmp);
  }
  ind = setKeyValueType(typeName, sectInd, key, fullStr, valType);
  free(fullStr);
  return ind;
}

/*!
 * Deletes the key-value pair matching [key] in the section with index
 * [sectInd] in the collection of sections of type [typeName].  Clears out
 * both the key and the value.  Returns -1 for error. 
 */  
int AdocDeleteKeyValue(const char *typeName, int sectInd, const char *key)
{
  AdocSection *sect;
  int keyInd;

  if (!(sect = getSection(typeName, sectInd)))
    return -1;
  keyInd = lookupKey(sect, key);
  if (keyInd < 0)
    return -1;
  B3DFREE(sect->values[keyInd]);
  B3DFREE(sect->keys[keyInd]);
  sect->types[keyInd] = ADOC_NO_VALUE;
  return 0;
}

/*
 * Routines for getting data from or modifying the current autodoc
 */

/*!
 * Returns the number of collections of sections, excluding the global data collection,
 * or -1 if there is no current autodoc.
 */
int AdocGetNumCollections()
{
  if (!sCurAdoc)
    return -1;
  return sCurAdoc->numCollections - 1;
}

/*!
 * Returns the name of the collection with index [collInd] into [string], which is 
 * allocated and should be freed with {free}.  Returns -1 for errors.
 */
int AdocGetCollectionName(int collInd, char **string)
{
  if (!sCurAdoc || collInd < 0 || collInd >= sCurAdoc->numCollections - 1)
    return -1;
  *string = strdup(sCurAdoc->collections[collInd].name);
  return (adocMemoryError(*string, "AdocGetCollectionName"));
}

/*!
 * Gets the name of the section with index [sectInd] in the collection of 
 * sections of type [typeName].  Returns the name in [string], which is allocated and 
 * should be freed with {free}.  Returns -1 for errors.
 */
int AdocGetSectionName(const char *typeName, int sectInd, char **string)
{
  AdocSection *sect;

  if (!(sect = getSection(typeName, sectInd)))
    return -1;
  *string = strdup(sect->name);
  return (adocMemoryError(*string, "AdocGetSectionName"));
}

/*!
 * Returns the number of sections of type [typeName].  Returns -1 for errors,
 * and 0 if there are no sections of the given type.
 */
int AdocGetNumberOfSections(const char *typeName)
{
  int collInd;
  if (!sCurAdoc || !typeName)
    return -1;
  collInd = lookupCollection(sCurAdoc, typeName);
  if (collInd < 0)
    return 0;
  return sCurAdoc->collections[collInd].numSections;
}

/*!
 * Returns the number of key-value pairs in the section with index [sectInd]
 * in the collection of sections of type [typeName].  Returns -1 for errors.
 */
int AdocGetNumberOfKeys(const char *typeName, int sectInd)
{
  AdocSection *sect;
  if (!(sect = getSection(typeName, sectInd)))
    return -1;
  return sect->numKeys;
}

/*!
 * Gets the key at index [keyInd] in the section with index [sectInd] in the collection 
 * of sections of type [typeName].  Returns a copy of the key in
 * [key], which should be freed with {free}, or returns NULL if this key has been deleted.
 *  Returns -1 for errors.
 */
int AdocGetKeyByIndex(const char *typeName, int sectInd, int keyInd, char **key)
{
  AdocSection *sect;
  if (!(sect = getSection(typeName, sectInd)))
    return -1;
  if (keyInd < 0 || keyInd >= sect->numKeys)
    return -1;
  *key = NULL;
  if (sect->keys[keyInd])
    *key = strdup(sect->keys[keyInd]);
  return 0;
}

/*!
 * Returns information about the value string matching [key] in the section with index 
 * [sectInd] in the collection of sections of type [typeName]: the value type 
 * (ADOC_ONE_INT, etc) in [valType] and the 
 * number of space-separated tokens in the string in [numTokens].  Returns -1 if the key 
 * is null, the section does not exist, or for a memory error; returns 1 if the
 * key does not occur in the given section or if the value is null.
 */
int AdocGetValTypeAndSize(const char *typeName, int sectInd, const char *key, 
                          int *valType, int *numTokens)
{
  AdocSection *sect;
  int keyInd;
  char *valstr, *parsed;
  *valType = ADOC_NO_VALUE;
  *numTokens = 0;
  if (!key)
    return -1;
  if (!(sect = getSection(typeName, sectInd)))
    return -1;
  keyInd = lookupKey(sect, key);
  if (keyInd < 0 || !sect->values[keyInd])
    return 1;
  *valType = sect->types[keyInd];

  /* Get a copy of the string and use the dreadful strtok */
  valstr = strdup(sect->values[keyInd]);
  if (!valstr) {
    adocMemoryError(NULL, "AdocGetValTypeAndSize");
    return -1;
  }
  parsed = valstr;
  while (strtok(parsed, " ") != NULL) {
    parsed = NULL;
    (*numTokens)++;
  }
  free(valstr);
  return 0;
}

/*!
 * Gets the value string matching [key] in the section with index [sectInd]
 * in the collection of sections of type [typeName].  Returns a copy of the
 * value in [string]; it should be freed with {free}.  Returns -1 if the key is
 * null, the section does not exist, or for a memory error; returns 1 if the 
 * key does not occur in the given section or if the value is null.  
 */
int AdocGetString(const char *typeName, int sectInd, const char *key, char **string)
{
  AdocSection *sect;
  int keyInd;

  if (!key)
    return -1;
  if (!(sect = getSection(typeName, sectInd)))
    return -1;
  keyInd = lookupKey(sect, key);
  if (keyInd < 0 || !sect->values[keyInd])
    return 1;
  *string = strdup(sect->values[keyInd]);
  return (adocMemoryError(*string, "AdocGetString"));
}

/*!
 * Like @AdocGetString, except that it extracts one integer from the value
 * string and returns its value in [val1].
 */
int AdocGetInteger(const char *typeName, int sectInd, const char *key, int *val1)
{
  int err;
  int num = 1;
  int tmp[1];
  if ((err = AdocGetIntegerArray(typeName, sectInd, key, tmp, &num, 1)) != 0)
    return err;
  *val1 = tmp[0];
  return 0;
}

/*! Like @AdocGetInteger except that it returns a float */
int AdocGetFloat(const char *typeName, int sectInd, const char *key, float *val1)
{
  int err;
  int num = 1;
  float tmp[1];
  if ((err = AdocGetFloatArray(typeName, sectInd, key, tmp, &num, 1)) != 0)
    return err;
  *val1 = tmp[0];
  return 0;
}

/*! Like @AdocGetInteger except that it returns a double */
int AdocGetDouble(const char *typeName, int sectInd, const char *key, double *val1)
{
  int err;
  int numToGet = 1;
  char *string;
  if ((err = AdocGetString(typeName, sectInd, key, &string)))
    return err;
  err = PipGetLineOfValues(string, string, (void *)val1, PIP_DOUBLE, &numToGet, 1);
  free(string);
  return err;
}

/*!
 * Like @AdocGetString, except that it extracts two integers from the value
 * string and returns their values in [val1] and [val2].
 */
int AdocGetTwoIntegers(const char *typeName, int sectInd, const char *key, int *val1,
                         int *val2)
{
  int err;
  int num = 2;
  int tmp[2];
  if ((err = AdocGetIntegerArray(typeName, sectInd, key, tmp, &num, 2)) != 0)
    return err;
  *val1 = tmp[0];
  *val2 = tmp[1];
  return 0;
}

/*! Like @AdocGetTwoIntegers except that it returns floats */
int AdocGetTwoFloats(const char *typeName, int sectInd, const char *key, float *val1,
                         float *val2)
{
  int err;
  int num = 2;
  float tmp[2];
  if ((err = AdocGetFloatArray(typeName, sectInd, key, tmp, &num, 2)) != 0)
    return err;
  *val1 = tmp[0];
  *val2 = tmp[1];
  return 0;
}

/*!
 * Like @AdocGetString, except that it extracts three integers from the value
 * string and returns their values in [val1], [val2], and [val3].
 */
int AdocGetThreeIntegers(const char *typeName, int sectInd, const char *key, int *val1,
                         int *val2, int *val3)
{
  int err;
  int num = 3;
  int tmp[3];
  if ((err = AdocGetIntegerArray(typeName, sectInd, key, tmp, &num, 3)) != 0)
    return err;
  *val1 = tmp[0];
  *val2 = tmp[1];
  *val3 = tmp[2];
  return 0;
}

/*! Like @AdocGetThreeIntegers except that it returns floats */
int AdocGetThreeFloats(const char *typeName, int sectInd, const char *key, float *val1,
                         float *val2, float *val3)
{
  int err;
  int num = 3;
  float tmp[3];
  if ((err = AdocGetFloatArray(typeName, sectInd, key, tmp, &num, 3)) != 0)
    return err;
  *val1 = tmp[0];
  *val2 = tmp[1];
  *val3 = tmp[2];
  return 0;
}

/*!
 * Like @AdocGetString, except that it extracts a set of integers from the
 * value string and returns their values in [array], whose size is given in 
 * [arraySize].  Set [numToGet] to the number of values to get, or 0 to get
 * all the values; in the latter case [numToGet] is returned with the number
 * of values retrieved.  Returns -1 for errors in parsing, too few values on 
 * the line, or not enough space in the array, as well as for failures in 
 * getting the value string.
 */
int AdocGetIntegerArray(const char *typeName, int sectInd, const char *key, int *array,
                        int *numToGet, int arraySize)
{
  char *string;
  int err;
  if ((err = AdocGetString(typeName, sectInd, key, &string)))
    return err;
  err = PipGetLineOfValues(string, string, (void *)array, PIP_INTEGER,
                           numToGet, arraySize);
  free(string);
  return err;
}

/*! Like @AdocGetIntegerArray except that it returns floats */
int AdocGetFloatArray(const char *typeName, int sectInd, const char *key, float *array,
                        int *numToGet, int arraySize)
{
  char *string;
  int err;
  if ((err = AdocGetString(typeName, sectInd, key, &string)))
    return err;
  err = PipGetLineOfValues(string, string, (void *)array, PIP_FLOAT,
                           numToGet, arraySize);
  free(string);
  return err;
}

/*!
 * DIRECT WRITING FUNCTIONS
 */
/*!
 * This is the first of a series of function that write directly to an open file instead
 * of building or modifying an autodoc structure; they have no effect on the current
 * autodoc.  This function writes one key-value pair with [key] and the integer [ival]
 * directly to file [fp]. Returns 1 for a write error.
 */
int AdocWriteInteger(FILE *fp, const char *key, int ival)
{
  if (fprintf(fp, "%s = %d\n", key, ival) < 0)
    return 1;
  return 0;
}

/*!
 * Like @AdocWriteInteger, except that the value is set to the two integers
 * [ival1] [ival2].
 */
int AdocWriteTwoIntegers(FILE *fp, const char *key, int ival1, int ival2)
{
  if (fprintf(fp, "%s = %d %d\n", key, ival1, ival2) < 0)
    return 1;
  return 0;
}

/*!
 * Like @AdocWriteInteger, except that the value is set to the three integers
 * [ival1] [ival2] [ival3].
 */
int AdocWriteThreeIntegers(FILE *fp, const char *key, int ival1, int ival2, int ival3)
{
  if (fprintf(fp, "%s = %d %d %d\n", key, ival1, ival2, ival3) < 0)
    return 1;
  return 0;
}

/*!
 * Like @AdocWriteInteger, except that the value is set to the array of [numVals]
 * integers in [ivals].
 */
int AdocWriteIntegerArray(FILE *fp, const char *key, int *ivals, int numVals)
{
  int ind;
  if (fprintf(fp, "%s =", key) < 0)
    return 1;
  for (ind = 0; ind < numVals; ind++)
    if (fprintf(fp, " %d", ivals[ind]) < 0)
      return 1;
  if (fprintf(fp, "\n") < 0)
    return 1;
  return 0;
}

/*! Like @AdocWriteInteger, except that the value is set to the float [val]. */
int AdocWriteFloat(FILE *fp, const char *key, float val)
{
  if (fprintf(fp, "%s = %g\n", key, val) < 0)
    return 1;
  return 0;
}

/*!
 * Like @AdocWriteInteger, except that the value is set to the two floats [val1] [val2].
 */
int AdocWriteTwoFloats(FILE *fp, const char *key, float val1, float val2)
{
  if (fprintf(fp, "%s = %g %g\n", key, val1, val2) < 0)
    return 1;
  return 0;
}

/*!
 * Like @AdocWriteInteger, except that the value is set to the thre floats
 * [val1] [val2] [val3].
 */
int AdocWriteThreeFloats(FILE *fp, const char *key, float val1, float val2, float val3)
{
  if (fprintf(fp, "%s = %g %g %g\n", key, val1, val2, val3) < 0)
    return 1;
  return 0;
}

/*!
 * Like @AdocWriteInteger, except that the value is set to the array of [numVals]
 * floats in [vals].
 */
int AdocWriteFloatArray(FILE *fp, const char *key, float *vals, int numVals)
{
  int ind;
  if (fprintf(fp, "%s =", key) < 0)
    return 1;
  for (ind = 0; ind < numVals; ind++)
    if (fprintf(fp, " %g", vals[ind]) < 0)
      return 1;
  if (fprintf(fp, "\n") < 0)
    return 1;
  return 0;
}

/* Like @AdocWriteInteger, except that the value is set to the double [val]. */
int AdocWriteDouble(FILE *fp, const char *key, double val)
{
  if (fprintf(fp, "%s = %g\n", key, val) < 0)
    return 1;
  return 0;
}

/*!
 * Writes one key-value pair directly to the file in [fp] from the strings in [key] and
 * [value].  Returns 1 for write error.
 */
int AdocWriteKeyValue(FILE *fp, const char *key, const char *value)
{
  if (fprintf(fp, "%s = %s\n", key, value) < 0)
    return 1;
  return 0;
}

/*!
 * Writes a section header directly to the file in [fp] from the strings in [key] and
 * [value].  Returns 1 for write error.
 */
int AdocWriteSectionStart(FILE *fp, const char *key, const char *value)
{
  if (fprintf(fp, "[%s = %s]\n", key, value ? value : "") < 0)
    return 1;
  return 0;
}

/*
 * XML HANDLING ROUTINES AND SOME COMMENT READING UTILITIES
 */
/*
 * Read the open file as an XML file
 */
static int readXmlFile(FILE *fp)
{
  mxml_node_t *xml, *node, *top, *sectNode, *child;
  const char *key, *value;
  int icol, ind, global, lastInd, err = 0;
  AdocSection *curSect;
  AdocCollection *coll;
  char **commentList = NULL;
  int maxComments = 0;
  int numComments = 0;
  /*FILE *ofp; */

  sLastWasXML = 1;
  rewind(fp);

  sNumSectNotElem = sNumSectNoName = sNumChildNotElem = sNumChildAttribs = 0;
  sNumValueNotText = sNumMultipleChilds = 0;
  
  xml = mxmlLoadFile(NULL, fp, MXML_OPAQUE_CALLBACK);
  if (!xml) {
    b3dError(stderr, "ERROR: AdocRead - Loading file as XML\n");
    return -2;
  }

  /* To walk through nodes and list them */
  /*
  node = mxmlWalkNext(xml, xml, MXML_DESCEND);
  while (node != NULL) {
    printf("node %p type %d ", node, node->type);
    if (node->type == MXML_ELEMENT) 
      printf("value :%s: numatt %d next %p prev %p parent %p child %p lastc %p\n",
             node->value.element.name, 
             node->value.element.num_attrs,
             node->next, node->prev, node->parent, node->child, node->last_child);
    else if (node->type == MXML_TEXT)
      printf("value :%s: next %p prev %p parent %p child %p lastc %p\n",
             node->value.text.string,
             node->next, node->prev, node->parent, node->child, node->last_child);
    else if (node->type == MXML_OPAQUE) {
      printf("value :");
      if (node->value.opaque) {
        unsigned char *chtemp = (unsigned char *)node->value.opaque;
        for (err = 0; err < strlen(node->value.opaque); err++)
          if (chtemp[err] > 31 && chtemp[err] < 128)
            printf("%c", chtemp[err]);
      }
      printf(": next %p prev %p parent %p child %p lastc %p\n",
             node->next, node->prev, node->parent, node->child, node->last_child);
    } else
      printf("\n");
    node = mxmlWalkNext(node, xml, MXML_DESCEND);
  }
  */

  /* Get the top node and then get the next if it is opaque node */
  top = mxmlWalkNext(xml, xml, MXML_DESCEND);
  if (top->type == MXML_OPAQUE)
    top = mxmlWalkNext(top, xml, MXML_DESCEND);
  key = mxmlGetElement(top);
  if (key)
    sCurAdoc->rootElement = strdup(key);

  /* Walk through the children of the top node, (autodoc) */
  sectNode = mxmlGetFirstChild(top);
  while (sectNode != NULL) {

    /* If a child of top is opaque, it is presumed whitespace and skip it */
    if (sectNode->type == MXML_OPAQUE) {
      sectNode = mxmlGetNextSibling(sectNode);
      continue;
    }
    
    if (testAndAddComment(sectNode, &commentList, &numComments, &maxComments)) {
      sectNode = mxmlGetNextSibling(sectNode);
      continue;
    }

    /* This is the section type if it an element node, which it really will be, and
     it must have the "name" attribute for the value */
    key = mxmlGetElement(sectNode);
    global = (key && !strcmp(key, GLOBAL_NAME)) ? 1 : 0;
    value = mxmlElementGetAttr(sectNode, "name");
    if (!key)
      sNumSectNotElem++;
    if (!global && !value)
      sNumSectNoName++;
    if (!key || (!global && !value)) {
      sectNode = mxmlGetNextSibling(sectNode);
      continue;
    }

    /* Lookup the collection under the key and create one if not found */
    icol = lookupCollection(sCurAdoc, key);
    if (icol < 0) {
      if ((err = addCollection(sCurAdoc, key)))
        break;
      icol = sCurAdoc->numCollections - 1;
    }
    coll = &sCurAdoc->collections[icol];
    
    /* Add a section to the collection and set it as current one */
    if (!global && (err = addSection(sCurAdoc, icol, value)))
      break;
    curSect = &coll->sections[coll->numSections - 1];
    lastInd = -1;
    if (numComments)
      if ((err = addComments(curSect, commentList, &numComments, lastInd)))
        break;
    
    /* Assign any other attributes as key-values in the section */
    for (ind = 0; ind < sectNode->value.element.num_attrs; ind++) {
      if (strcmp(sectNode->value.element.attrs[ind].name, "name")) {
        if ((err = sectSetKeyValueType(curSect, sectNode->value.element.attrs[ind].name,
                                       sectNode->value.element.attrs[ind].value, 
                                       ADOC_STRING, &lastInd)))
          break;
      }
    }

    /* Now walk through the children of section node */
    node = mxmlGetFirstChild(sectNode);
    while (node != NULL) {

      /* Again, skip children that are opaque and definitely white space */
      if (node->type == MXML_OPAQUE && node->value.opaque && 
          node->value.opaque[0] == '\n') {
        node = mxmlGetNextSibling(node);
        continue;
      }

      if (testAndAddComment(node, &commentList, &numComments, &maxComments)) {
        node = mxmlGetNextSibling(node);
        continue;
      }
      
      /* It must be an element */
      key = mxmlGetElement(node);
      if (!key) 
        sNumChildNotElem++;
      else {
        if (node->value.element.num_attrs)
          sNumChildAttribs++;
        child = mxmlGetFirstChild(node);

        /* The first child should be opaque and there should be only one */
        if (child && mxmlGetType(child) != MXML_OPAQUE)
          sNumValueNotText++;
        else {
          if (child && mxmlGetLastChild(node) != child)
            sNumMultipleChilds++;
          value = child ? child->value.opaque : NULL;
          if ((err = sectSetKeyValueType(curSect, key, value, ADOC_STRING, &lastInd)))
            break;
          if (numComments)
            if ((err = addComments(curSect, commentList, &numComments, lastInd)))
              break;
        }
      }

      /* Step to next key-value in section, if any */
      node = mxmlGetNextSibling(node);
    }
    if (err)
      break;

    /* Add any comments that have accumulated */
    if (numComments)
      if ((err = addComments(curSect, commentList, &numComments, lastInd)))
        break;

    /* Step to next section if any */
    sectNode = mxmlGetNextSibling(sectNode);
  }

  /* This is what needs to be done after reading in an xml to save as an xml with the
     white space callback  and opaque reading */
  /*
  ofp = fopen("read-in.xml", "w");
  if (ofp) {
  mxml_node_t *lastNode;
    node = mxmlWalkNext(xml, xml, MXML_DESCEND);
    lastNode = xml;
    while (node != NULL) {
      if (node->type == MXML_OPAQUE && node->value.opaque &&
          node->value.opaque[0] == '\n') {
        mxmlRemove(node);
        mxmlDelete(node);
        //node->value.opaque[0] = 0x00;
          node = lastNode;
      }
      lastNode = node;
      node = mxmlWalkNext(node, xml, MXML_DESCEND);
    }
    ixmlResetLastLevel();
    mxmlSaveFile(xml, ofp, (mxml_save_cb_t)ixmlWhitespace_cb);
    fclose(ofp);
  }
  */

  if (err)
    deleteAdoc(sCurAdoc);
  handleFinalComments(err, &commentList, numComments);

  mxmlDelete(xml);
  return err;
}

/*
 * Add a new comment to the list of comments for an item
 */
static int addToCommentList(char ***commentList, int *numComments, int *maxComments,
                            char *comment)
{
  if (*numComments >= *maxComments) {
    if (*maxComments)
      *commentList = (char **)realloc(*commentList,
                                      (*maxComments + 1) * sizeof(char *));
    else
      *commentList = (char **)malloc(sizeof(char *));
    if (adocMemoryError(*commentList, "AdocRead"))
      return 1;
    (*maxComments)++;
  }
  (*commentList)[*numComments] = strdup(comment);
  if (adocMemoryError((*commentList)[(*numComments)++], "AdocRead"))
    return 1;
  return 0;
}

/*
 * Test for an XML comment node and if do, strip the comment caharcts, add a # on front
 * and add to list of comments
 */
static int testAndAddComment(mxml_node_t *node, char ***commentList, int *numComments, 
                               int *maxComments)
{
  const char *key;
  char *tmpStr;
  int len;
  key = mxmlGetElement(node);
  if (node->type == MXML_ELEMENT && key && PipStartsWith(key, XML_COMMENT_START)) {
    tmpStr = strdup(&key[XML_COMSTART_LEN - 1]);
    if (tmpStr) {
      tmpStr[0] = '#';
      len = strlen(tmpStr);
      if (len >= 2 && tmpStr[len - 2] == '-' && tmpStr[len -1] == '-') {
        len -= 2;
        tmpStr[len] = 0x00;
      }
      if (len)
        addToCommentList(commentList, numComments, maxComments, tmpStr);
      free(tmpStr);
    }
    return 1;
  }
  return 0;
}

/*
 * At end of reading, clean up comment list if there was an error, or assign it to the 
 * finalComments element of the autodoc
 */
static void handleFinalComments(int err, char ***commentList, int numComments)
{
  int i;
  if (err) {

    /* Clean out comment list */ 
    for (i = 0; i < numComments; i++)
      if (*commentList[i])
        free(*commentList[i]);
    if (*commentList)
      free(*commentList);
  } else if (*commentList) {

    /* If good, transfer any comments to the autodoc */
    if (numComments) {
      sCurAdoc->finalComments = *commentList;
      sCurAdoc->numFinalCom = numComments;
    } else
      free(*commentList);
  }
}

/*
 * Writes the current autodoc as an XML file.
 */
static int writeXmlFile(const char *filename)
{
  FILE *afile;
  mxml_node_t *xml, *node, *elem, *top;
  int i, j, k, ind, useInd, comInd;
  AdocCollection *coll;
  AdocSection *sect;
  int *ordSectInds;
  int orderedWrite = sNameForOrdering && sCurAdoc->numSections > 1;
  /*double wallStart = wallTime(); */

  if (!sCurAdoc)
    return -1;
  afile = fopen(filename, "w");
  if (!afile)
    return -1;
  if (orderedWrite) {
    ordSectInds = setupSectionOrder();
    if (!ordSectInds)
      return -1;
  }

  xml = mxmlNewXML("1.0");
  top = mxmlNewElement(xml, B3DCHOICE(sCurAdoc->rootElement != NULL, 
                                      sCurAdoc->rootElement, "autodoc"));
  for (ind = 0; ind < sCurAdoc->numSections; ind++) {
    useInd = orderedWrite ? ordSectInds[ind] : ind;
    comInd = 0;
    i = sCurAdoc->collList[useInd];
    j = sCurAdoc->sectList[useInd];
    coll = &sCurAdoc->collections[i];
    sect = &coll->sections[j];
    for (; comInd < sect->numComments && sect->comIndex[comInd] == -1; comInd++)
      writeCommentToXML(top, sect->comments[comInd]);
    node = mxmlNewElement(top, coll->name);
    if (i || j || strcmp(sect->name, GLOBAL_NAME)) {
      mxmlElementSetAttr(node, "name", sect->name);
    }

    /* Loop on key-values */
    for (k = 0; k < sect->numKeys; k++) {
      for (; comInd < sect->numComments && sect->comIndex[comInd] == k; comInd++)
        writeCommentToXML(node, sect->comments[comInd]);
      elem = mxmlNewElement(node, sect->keys[k]);
      if (sect->values[k])
        mxmlNewText(elem, 0, sect->values[k]);
    }
  }
  for (i = 0; i < sCurAdoc->numFinalCom; i++)
    writeCommentToXML(top, sCurAdoc->finalComments[i]);

  /*printf("xml build time %.1f\n", (wallTime() - wallStart) * 1000.);
    wallStart = wallTime();*/
  /*
  node = mxmlWalkNext(xml, xml, MXML_DESCEND);
  while (node != NULL) {
    printf("type %d text %s value %s numatt %d next %s prev %s parent %s child %s lastc %s\n",
           node->type, mxmlGetText(node, NULL), node->type? node->value.text.string : node->value.element.name, 
           node->type ? 0 : node->value.element.num_attrs,
           node->next ? node->next->value.element.name : "NULL",
           node->prev ? node->prev->value.element.name : "NULL",
           node->parent ? node->parent->value.element.name : "NULL",
           node->child ? node->child->value.element.name : "NULL",
           node->last_child ? node->last_child->value.element.name : "NULL");
      node = mxmlWalkNext(node, xml, MXML_DESCEND);
      }
  */
  mxmlSetWrapMargin(0);
  ixmlResetLastLevel();
  ind = mxmlSaveFile(xml, afile, (mxml_save_cb_t)ixmlWhitespace_cb);
  fclose(afile);
  /*printf("xml write time %.1f\n", (wallTime() - wallStart) * 1000.);*/
  mxmlDelete(xml);
  return ind;
}

/*
 * Write one comment to the XML file
 */
static void writeCommentToXML(mxml_node_t *parent, char *comment)
{
  int len;
  char *tmpStr;
  mxml_node_t *node;

  len = strlen(comment);
  if (!len)
    return;
  tmpStr = (char *)malloc(len + 10);
  if (tmpStr) {
    sprintf(tmpStr, XML_COMMENT_START"%s%s--", &comment[1], 
            comment[len - 1] == ' ' ? "" : " ");
    node = mxmlNewElement(parent, tmpStr);
    free(tmpStr);
  }
}

/*
 * Routines for adding keys, sections, collections, autodocs 
 */

/* Adds a key-value pair to the given section, without checking for 
   duplication */
static int addKey(AdocSection *sect, const char *key, const char *value, int type)
{

  /* First allocate enough memory if needed */
  if (!sect->maxKeys) {
    sect->keys = (char **)malloc(MALLOC_CHUNK * sizeof(char *));
    sect->values = (char **)malloc(MALLOC_CHUNK * sizeof(char *));
    sect->types = B3DMALLOC(unsigned char, MALLOC_CHUNK);
    sect->maxKeys = MALLOC_CHUNK;
  } else if (sect->numKeys >= sect->maxKeys) {
    sect->keys = (char **)realloc(sect->keys, (sect->maxKeys + MALLOC_CHUNK) * 
                                  sizeof(char *));
    sect->values = (char **)realloc(sect->values, 
                                    (sect->maxKeys + MALLOC_CHUNK) *
                                    sizeof(char *));
    B3DREALLOC(sect->types, unsigned char, sect->maxKeys + MALLOC_CHUNK);
    sect->maxKeys += MALLOC_CHUNK;
  }
  if (!sect->keys || !sect->values || !sect->types) {
    adocMemoryError(NULL, "addKey");
    return -1;
  }

  /* Copy key and value and increment count */
  sect->keys[sect->numKeys] = strdup(key);
  if (value)
    sect->values[sect->numKeys] = strdup(value);
  else
    sect->values[sect->numKeys] = NULL;
  sect->types[sect->numKeys] = value ? type : ADOC_NO_VALUE;
  if (! sect->keys[sect->numKeys] || (value && !sect->values[sect->numKeys])) {
    adocMemoryError(NULL, "addKey");
    return -1;
  }
  sect->numKeys++;
  return 0;
}

/* Adds a section of the given name to the collection */
static int addSection(Autodoc *adoc, int collInd, const char *name)
{
  AdocCollection *coll = &adoc->collections[collInd];
  AdocSection *sect;

  /* First allocate enough memory if needed for the sections in the collection
     and for the master lists in the autodoc */
  if (!coll->maxSections) {
    coll->sections = (AdocSection *)malloc(MALLOC_CHUNK * sizeof(AdocSection));
    coll->maxSections = MALLOC_CHUNK;
  } else if (coll->numSections >= coll->maxSections) {
    coll->sections = (AdocSection *)
      realloc(coll->sections, (coll->maxSections + MALLOC_CHUNK)
              * sizeof(AdocSection));
    coll->maxSections += MALLOC_CHUNK;
  }
  if (adocMemoryError(coll->sections, "addSection"))
    return -1;

  if (!adoc->maxSections) {
    adoc->collList = (int *)malloc(MALLOC_CHUNK * sizeof(int));
    adoc->sectList = (int *)malloc(MALLOC_CHUNK * sizeof(int));
    adoc->maxSections = MALLOC_CHUNK;
  } else if (adoc->numSections >= adoc->maxSections) {
    adoc->collList = (int *)realloc
      (adoc->collList, (adoc->maxSections + MALLOC_CHUNK) * sizeof(int));
    adoc->sectList = (int *)realloc
      (adoc->sectList, (adoc->maxSections + MALLOC_CHUNK) * sizeof(int));
    adoc->maxSections += MALLOC_CHUNK;
  }
  if (!adoc->collList || !adoc->sectList) {
    adocMemoryError(NULL, "addSection");
    return -1;
  }

  /* Copy the name and initialize to empty keys */
  sect = &coll->sections[coll->numSections];
  sect->name = strdup(name);
  if (adocMemoryError(sect->name, "addSection"))
    return -1;
  sect->keys = NULL;
  sect->values = NULL;
  sect->types = NULL;
  sect->numKeys = 0;
  sect->maxKeys = 0;
  sect->comments = NULL;
  sect->comIndex = NULL;
  sect->numComments = 0;

  /* Add the collection and section # to master list */
  adoc->collList[adoc->numSections] = collInd;
  adoc->sectList[adoc->numSections++] = coll->numSections++;
  return 0;
}

/* Adds a collection of the given name to the autodoc */
static int addCollection(Autodoc *adoc, const char *name)
{
  AdocCollection *coll;

  /* Allocate just one at a time when needed */
  if (!adoc->numCollections)
    adoc->collections = (AdocCollection *)malloc(sizeof(AdocCollection));
  else
    adoc->collections = (AdocCollection *)realloc(adoc->collections,
                                                  (adoc->numCollections + 1) *
                                                  sizeof(AdocCollection));
  if (adocMemoryError(adoc->collections, "addCollection"))
    return -1;
  coll = &adoc->collections[adoc->numCollections];
  coll->name = strdup(name);
  if (adocMemoryError(coll->name, "addCollection"))
    return -1;
  coll->numSections = 0;
  coll->maxSections = 0;
  coll->sections = NULL;
  adoc->numCollections++;
  return 0;
}

/* Adds an autodoc to the array.  Initializes it with a global collection and
   section, and clears it out if this fails */
static int addAutodoc()
{
  Autodoc *adoc;
  int index = -1, i;
  
  /* Search for a free autodoc in array */
  for (i = 0; i < sNumAutodocs; i++) {
    if (!sAutodocs[i].inUse) {
      index = i;
      break;
    }
  }

  if (index < 0) {

    /* Allocate just one at a time when needed */
    if (!sNumAutodocs)
      sAutodocs = (Autodoc *)malloc(sizeof(Autodoc));
    else
      sAutodocs = (Autodoc *)realloc(sAutodocs, (sNumAutodocs + 1) *
                                    sizeof(Autodoc));
    if (adocMemoryError(sAutodocs, "addAutodoc"))
      return -1;
    index = sNumAutodocs++;
  }

  /* Initialize collections */
  adoc = &sAutodocs[index];
  adoc->collections = NULL;
  adoc->numCollections = 0;
  adoc->finalComments = NULL;
  adoc->numFinalCom = 0;
  adoc->collList = NULL;
  adoc->sectList = NULL;
  adoc->numSections = 0;
  adoc->maxSections = 0;
  adoc->inUse = 1;
  adoc->backedUp = 0;
  adoc->writeAsXML = 0;
  adoc->rootElement = NULL;

  /* Add a collection and section for global data */
  if (addCollection(adoc, GLOBAL_NAME))
    return -1;
  if (addSection(adoc, 0, GLOBAL_NAME)) {
    deleteAdoc(adoc);
    return -1;
  }
  return index;
}

/* 
 * Utility routines for freeing, parsing and lookup
 */
/* Frees all data in an autodoc and zero out the top-level items so this can
   be called more than once */
static void deleteAdoc(Autodoc *adoc)
{
  AdocSection *sect;
  AdocCollection *coll;
  int i, j, k;
  for (i = 0; i < adoc->numCollections; i++) {
    coll = &adoc->collections[i];
    for (j = 0; j < coll->numSections; j++) {
      sect = &coll->sections[j];

      /* Clean key/values out of sections */
      for (k = 0; k < sect->numKeys; k++) {
        B3DFREE(sect->keys[k]);
        B3DFREE(sect->values[k]);
      }
      B3DFREE(sect->keys);
      B3DFREE(sect->values);
      B3DFREE(sect->types);
      B3DFREE(sect->name);

      /* Clean comments out of sections */
      for (k = 0; k < sect->numComments; k++)
        B3DFREE(sect->comments[k]);
      B3DFREE(sect->comments);
      B3DFREE(sect->comIndex);
    }


    /* Free sections */
    B3DFREE(coll->sections);
    B3DFREE(coll->name);
  }

  /* Free collections */
  B3DFREE(adoc->collections);
  adoc->numCollections = 0;

  /* Free lists of sections */
  B3DFREE(adoc->collList);
  B3DFREE(adoc->sectList);
  adoc->numSections = 0;
  adoc->maxSections = 0;
  
  /* Free final comments */
  for (i = 0; i < adoc->numFinalCom; i++)
    B3DFREE(adoc->finalComments[i]);
  B3DFREE(adoc->finalComments);
  adoc->numFinalCom = 0;
  B3DFREE(adoc->rootElement);
  adoc->inUse = 0;
}

/* Parses the characters in line up to (not including) end for the construct
   key = value.  Returns an allocated copy of key in key and value in value.
   Returns NULL in value if there is no value.  Returns 1 for other malformed
   lines, and -1 for other errors. */
static int parseKeyValue(char *line, char *end, char **key, char **value)
{
  char *valStart, *keyEnd;
  int keyLen, valLen;

  /* Eat spaces at start and end */
  while (line < end && (*line == ' ' || *line == '\t'))
    line++;
  while (line < end && (*(end - 1) == ' ' || *(end - 1) == '\t'))
    end--;
  if (line == end)
    return 1;
         
  /* Find delimiter.  If it is not there or no text before it, error */
  valStart = strstr(line, sValueDelim);
  if (!valStart || valStart == line)
    return 1;
  
  /* Eat spaces after key */
  keyEnd = valStart;
  while (keyEnd > line && (*(keyEnd - 1) == ' ' || *(keyEnd - 1) == '\t'))
    keyEnd--;
  
  /* Eat spaces after the delimiter.  Allow an empty value */
  valStart += strlen(sValueDelim);
  while (valStart < end && (*valStart == ' ' || *valStart == '\t'))
    valStart++;

  /* Allocate for strings and copy them */
  keyLen = keyEnd - line;
  *key = (char *)malloc(keyLen + 1);
  if (adocMemoryError(*key, "parseKeyValue"))
    return -1;
  memcpy(*key, line, keyLen);
  (*key)[keyLen] = 0x00;

  valLen = end - valStart;
  *value = NULL;
  if (valLen) {
    *value = (char *)malloc(valLen + 1);
    if (adocMemoryError(*value, "parseKeyValue"))
      return -1;
    memcpy(*value, valStart, valLen);
    (*value)[valLen] = 0x00;
  }
  return 0;
}

/* Looks up a key in a section and returns its index, or -1 if not present */
static int lookupKey(AdocSection *sect, const char *key)
{
  int i;
  if (!key)
    return -1;
  for (i = 0; i < sect->numKeys; i++)
    if (sect->keys[i] && !strcmp(key, sect->keys[i]))
      return i;
  return -1;
}

/* Looks up a collection in the autodoc by name; returns index or -1 if not
   there */
static int lookupCollection(Autodoc *adoc, const char *name)
{
  int i;
  if (!name || !adoc)
    return -1;
  for (i = 0; i < adoc->numCollections; i++)
    if (!strcmp(name, adoc->collections[i].name))
      return i;
  return -1;
}

/* Returns the section in the given collection and with given index, or NULL
   for error */
static AdocSection *getSection(const char *typeName, int sectInd)
{
  AdocCollection *coll;
  int collInd;
  if (!sCurAdoc || !typeName || sectInd < 0)
    return NULL;
  collInd = lookupCollection(sCurAdoc, typeName);
  if (collInd < 0)
    return NULL;
  coll = &sCurAdoc->collections[collInd];
  if (sectInd >= coll->numSections)
    return NULL;
  return (&coll->sections[sectInd]);
}

/* Adds the comments in the list of comments to the section, taking the
   strings without copying and setting numComments to zero */
static int addComments(AdocSection *sect, char **comments, int *numComments,
                       int index)
{
  int i, newNum = sect->numComments + *numComments;
  if (sect->numComments) {
    sect->comments = (char **)realloc(sect->comments, newNum * sizeof(char *));
    sect->comIndex = (int *)realloc(sect->comIndex, newNum * sizeof(int));
  } else {
    sect->comments = (char **)malloc(newNum * sizeof(char *));
    sect->comIndex = (int *)malloc(newNum * sizeof(int));
  }
  if (!sect->comments || !sect->comIndex) {
    adocMemoryError(NULL, "addComments");
    return -1;
  }
  for (i = 0; i < *numComments; i++) {
    sect->comments[sect->numComments] = comments[i];
    sect->comIndex[sect->numComments++] = index;
  }
  *numComments = 0;
  return 0;
}

/* Find the index of a section in the master list */
static int findSectionInAdocList(int collInd, int sectInd)
{
  int i;
  for (i = 0; i < sCurAdoc->numSections; i++)
    if (sCurAdoc->collList[i] == collInd && sCurAdoc->sectList[i] == sectInd)
      return i;
  return -1;
}

static int adocMemoryError(void *ptr, const char *routine)
{
  if (ptr)
    return 0;
  b3dError(stderr, "ERROR: %s - Allocating memory for string or autodoc component\n",
           routine);
  return -1;
}
