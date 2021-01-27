#ifndef HVEMTYPES_H
#define HVEMTYPES_H
/*
  $Id$
*/


#include <sys/types.h>

/* Read the definitions of the b3d... data types of defined bit size */
#include "imodconfig.h"

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#define b3dX 0
#define b3dY 1
#define b3dZ 2

#include <limits.h>

#ifndef _WIN32
#ifndef FLT_MAX 
#define FLT_MAX         1.E+37f
#endif
#endif

#ifndef SEEK_SET
#define SEEK_SET 0
#endif
#ifndef SEEK_CUR
#define SEEK_CUR 1
#endif
#ifndef SEEK_END
#define SEEK_END 2
#endif

#endif /* hvemtypes.h */
