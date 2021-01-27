/* writelist - For converting an integer list to ranges and printing it
 * $Id$
 */

#include "b3dutil.h"
#include "imodconfig.h"

#ifdef F77FUNCAP
#define writelist WRITELIST
#define wrlist WRLIST
#else
#define writelist writelist_
#define wrlist wrlist_
#endif

static char *addRangeToLine(int numStart, int numEnd, char *line);

/*!
 * Prints a list of comma-separated ranges for the [numList] values in [list], breaking
 * the list into multiple lines of maximum length [lineLen] if necessary.  Only 
 * ranges in increasing order are recognized.  Returns 1 for memory error.
 */
int writeList(int *list, int numList, int lineLen)
{
  int ind, len;
  char *lineStr = listToString(list, numList);
  char *saveStr = lineStr;
  if (!lineStr)
    return 1;
  while (strlen(lineStr) > lineLen) {
    for (ind = lineLen - 1; ind > 0; ind--)
      if (lineStr[ind] == ',')
        break;
    lineStr[ind] = '\n';
    lineStr += ind + 1;
  }
  printf(saveStr);
  printf("\n");
  fflush(stdout);
  free(saveStr);
  return 0;
}

/*! Fortran wrapper for writeList */
int writelist(int *list, int *numList, int *lineLen)
{
  return writeList(list, *numList, *lineLen);
}

/*! Fortran wrapper for old wrlist subroutine that used a line length of 80 */
void wrlist(int *list, int *numList)
{
  writeList(list, *numList, 80);
}

/*!
 * Converts the [numList] values in [list] into a string of comma-separated ranges, where
 * only ranges in increasing order are recognized.  The returned string is allocated with
 * malloc() and should be freed with free().  Returns NULL for memory error.
 */
char *listToString(int *list, int numList)
{
  char *retString = NULL;
  int ind;
  int numStart = list[0];
  for (ind = 1; ind < numList; ind++) {
    if (list[ind] != list[ind - 1] + 1) {
      retString = addRangeToLine(numStart, list[ind - 1], retString);
      if (!retString)
        return NULL;
      numStart = list[ind];
    }
  }
  retString = addRangeToLine(numStart, list[numList - 1], retString);
  return retString;
}

static char *addRangeToLine(int numStart, int numEnd, char *line)
{
  char rangeString[32];
  char endString[16];
char *retLine;
sprintf(rangeString, "%d", numStart);
  if (numEnd > numStart) {
    strcat(rangeString, "-");
    sprintf(endString, "%d", numEnd);
    strcat(rangeString, endString);
  }

  if (line) {
    retLine = realloc(line, strlen(line) + strlen(rangeString) + 4);
    if (!retLine)
      return NULL;
    strcat(retLine, ",");
    strcat(retLine, rangeString);
  } else {
    retLine = (char *)malloc(strlen(rangeString) + 4);
    if (!retLine)
      return NULL;
    strcpy(retLine, rangeString);
  }
  return retLine;  
}
