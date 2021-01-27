/*
 *  imodmesh.c -- Add mesh data to a model.
 *
 *  Original author: James Kremer
 *  Revised by: David Mastronarde   email: mast@colorado.edu
 *
 *  Copyright (C) 1995-2006 by Boulder Laboratory for 3-Dimensional Electron
 *  Microscopy of Cells ("BL3DEMC") and the Regents of the University of 
 *  Colorado.  See dist/COPYRIGHT for full copyright notice.
 *
 *  $Id$
 */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <math.h>

#include "imodel.h"
#include "b3dutil.h"
#include "mkmesh.h"

static void imodMeshesRescaleNormal(Imesh *inMesh, int *meshsize,
                                    Ipoint *spnt);
static void imodMeshesDeleteZRange(Imesh **meshp, int *meshsize, int minz,
                                   int maxz, int resol);
static void imodMeshDeleteZRange(Imesh *mesh, int minz, int maxz);
static int ObjOnList(int ob, int list[], int nlist);
static int deletez(Ipoint point, int minz, int maxz);

#define DEFAULT_TOL  0.25
#define DEFAULT_FLAT 1.5

static int imodmesh_usage(char *prog, int retcode)
{
  fprintf(stderr, "%s usage: %s [options] [input files...]\n",prog,prog);
  fprintf(stderr, "options:\n");
  fprintf(stderr, "\t-u  \tUse stored meshing parameters.\n");
  fprintf(stderr, "\t-c  \tCap ends of object.\n");
  fprintf(stderr, "\t-C  \tCap all unconnected ends of object.\n");
  fprintf(stderr, "\t-D list\tDo not cap ends unconnected to Z values in the "
          "list.\n");
  fprintf(stderr, "\t-p #\tOnly connect contours that overlap by the "
          "given percentage.\n");
  fprintf(stderr, "\t-s  \tConnect contours across sections with no contours."
          "\n");
  fprintf(stderr, "\t-P #\tDo multiple passes to connect contours farther "
          "apart in Z.\n");
  fprintf(stderr, "\t-S  \tUse surface numbers for connections.\n");
  fprintf(stderr, "\t-F #\tCriterion Z difference for analyzing for tilted "
          "contours.\n");
  fprintf(stderr, "\t-I  \tIgnore time values and connect across times.\n");
  fprintf(stderr, "\t-f  \tForce more connections to non-overlapping "
          "contours.\n");
  fprintf(stderr, "\t-t list\tRender open contour objects in list as tubes."
          "\n");
  fprintf(stderr, "\t-d #\tSet diameter for tubes (default is 3D line width)."
          "\n");
  fprintf(stderr, "\t-E  \tCap ends of tubes.\n");
  fprintf(stderr, "\t-H  \tCap ends of tubes with hemispheres.\n");
  fprintf(stderr, "\t-T  \tDo time consuming calculations, "
          "may help reduce artifacts.\n");
    
  fprintf(stderr, "\t-o list\tDo operations only on objects in list"
          " (ranges allowed).\n");
  fprintf(stderr, "\t-R #\tTolerance (maximum error) for point reduction"
          " (default %.2f).\n", DEFAULT_TOL);
    
  fprintf(stderr, "\t-i #\tOnly mesh sections at the given z increment.\n");
  fprintf(stderr, "\t-z #,#,#  Only mesh sections within the given z range "
          "at a z increment.\n");
  fprintf(stderr, "\t-l \tMark new meshes as low-resolution and keep "
          "high res meshes.\n");
  fprintf(stderr, "\t-x #,#\tClip mesh outside given lower and upper limits "
          "in X.\n");
  fprintf(stderr, "\t-y #,#\tClip mesh outside given lower and upper limits "
          "in Y.\n");
  fprintf(stderr, "\t-no* \tOverride a stored parameter; * is one of "
          "cCDsSIftEHTxyz.\n");

  fprintf(stderr, "\t-a \tAppend mesh data to object, replacing mesh"
          " in same z range.\n");
  fprintf(stderr, "\t-e  \tErase meshes, overrides other options.\n");
  fprintf(stderr, "\t-N \tRecompute normals for existing meshes.\n");
  fprintf(stderr, "\t-n \tRescale normals using the -Z value "
          "for existing meshes.\n");
  fprintf(stderr, "\t-Z #\tNormal scaling z multiplier.\n");
  fprintf(stderr, "\t-B  \tMake mesh backward-compatible to IMOD before 3.6.14.\n");
  if (retcode)
    exit(retcode);
  return(retcode);
}

#define DELETE_INDEX  -123456
#define LOWRES_INCZ  4
#define LOWRES_TOL   2.

int main(int argc, char **argv)
{
  Imod *imod;
  Iobj *obj;
  Ipoint spnt;
  Ipoint max;
  FILE *fout;
  Iobjview *obv;
  int ob, iarg, iv;
  int erase_mesh = FALSE;
  int renorm     = FALSE;
  int remeshNorm = FALSE;
  int cap = IMESH_CAP_OFF;
  int times = TRUE;
  int useOldParam = FALSE;
  int notimes = FALSE;
  int nocap = FALSE;
  int noSkipList = FALSE;
  int noxmin = FALSE;
  int noymin = FALSE;
  int nozmin = FALSE;
  int notube = FALSE;
  int dowarn1 = TRUE;
  int dowarn2 = TRUE;
  int dowarn3 = TRUE;
  unsigned int setFlags = IMESH_MK_NORM | IMESH_MK_TIME;
  unsigned int clearFlags = 0;
  unsigned int flags;
  float overlap = -1.;
  float tubeDiameter = -99.;
  float flatCrit = -1.;
  float tol = -1.;
  int resol = 0;
  int passes = -1;

  /* Added 2.00 mesh only given sections. */
  int minz = DEFAULT_VALUE;
  int maxz = DEFAULT_VALUE;
  int incz = -1;
  float zscale = 1.0;
  int append = 0;
  Imesh *tmsh;
  int    tmshsize, m;
  Ipoint triMin = {-DEFAULT_FLOAT, -DEFAULT_FLOAT, -DEFAULT_FLOAT};
  Ipoint triMax = {DEFAULT_FLOAT, DEFAULT_FLOAT, DEFAULT_FLOAT};
  
  int newPolyNorm = 1;
  int nlist = 0;
  int *list;
  int cap_skip_nz = 0;
  int *cap_skip_zlist = NULL;
  int ilist, useParam;
  int ntube_list = 0;
  int *tube_list;
  MeshParams *param = NULL;
  char *progname = imodProgName(argv[0]);
  char oldvar[] = "IMODMESH_OLDMESH";

  if (argc < 1){
    imodVersion(progname);
    imodCopyright();
    imodmesh_usage(progname, -1);
  }
    
  for (iarg = 1; iarg < argc ; iarg++){
    if (argv[iarg][0] == '-'){
      if (!(argv[iarg][1] == 'n' && argv[iarg][2] == 'o')) {

        switch (argv[iarg][1]){
        case 'd':
          tubeDiameter = atof(argv[++iarg]);
          break;
        case 'o': /* select object list */
          list = parselist(argv[++iarg], &nlist);
          if (!list) {
            fprintf(stderr, "%s: Error parsing object list\n", progname);
            exit(3);
          }
          break;
          /* 9/7/06: Removed -r option */
        case 'R':
          tol = atof(argv[++iarg]);
          if (tol < 0.)
            tol = 0.;
          break;
        case 'F':
          flatCrit = atof(argv[++iarg]);
          break;
        case 'c':
          cap = IMESH_CAP_END;
          break;
        case 'C':
          cap = IMESH_CAP_ALL;
          break;
        case 'D': /* Do not cap at Z in list. */
          cap_skip_zlist = parselist(argv[++iarg], &cap_skip_nz);
          if (!cap_skip_zlist) {
            fprintf(stderr, "%s: Error parsing Z list\n", progname);
            exit(3);
          }
          break;
        case 'h':
          imodmesh_usage(progname, -1);
          break;
        case 'e': /* erase meshes */
          erase_mesh = TRUE;
          break;
        case 'f':
          setFlags |= IMESH_MK_STRAY;
          break;
        case 'E':
          setFlags |= IMESH_MK_CAP_TUBE;
          break;
        case 'H':
          setFlags |= IMESH_MK_CAP_DOME;
          break;
        case 'T':
          setFlags |= IMESH_MK_FAST;
          break;
        case 's':
          setFlags |= IMESH_MK_SKIP;
          break;
        case 'S':
          setFlags |= IMESH_MK_SURF;
          break;
        case 'I':
          times = FALSE;
          break;
        
        case 'n':
          renorm = TRUE;
          break;
        case 'N':
          remeshNorm = TRUE;
          break;
        
        case 't':
          tube_list = parselist(argv[++iarg], &ntube_list);
          if (!tube_list) {
            fprintf(stderr, "%s: Error parsing list of objects to "
                    "mesh as tubes\n", progname);
            exit(3);
          }
          break;
            
        
        case 'z':
          sscanf(argv[++iarg], "%d%*c%d%*c%d", &minz, &maxz, &incz);
          break;
        case 'i':
          incz = atoi(argv[++iarg]);
          if (incz < 1)
            incz = 1;
          break;
        case 'x':
          sscanf(argv[++iarg], "%f%*c%f", &triMin.x, &triMax.x);
          break;
        case 'y':
          sscanf(argv[++iarg], "%f%*c%f", &triMin.y, &triMax.y);
          break;
        case 'Z':
          zscale = atof(argv[++iarg]);
          break;
        case 'a':
          append = TRUE;
          break;
        case 'p':
          overlap = atof(argv[++iarg]) / 100.;
          break;
        case 'l':
          resol = 1;
          break;
        case 'P':
          passes = atoi(argv[++iarg]);
          if (passes < 1)
            passes =1;
          break;
        case 'B':
          newPolyNorm = 0;
          break;
        case 'u':
          useOldParam = TRUE;
          break;
        
        case '?':
          imodmesh_usage(progname, -1);
          break;
        default:
          fprintf(stderr, "%s: unknown option %s\n", progname, argv[iarg]);
          imodmesh_usage(progname, -1);
          break;
               
        }
      } else {
        switch (argv[iarg][3]) {
        case 'c':
        case 'C':
          nocap = TRUE;
          break;
        case 'D':
          noSkipList = TRUE;
          break;
        case 'f':
          clearFlags |= IMESH_MK_STRAY;
          break;
        case 'E':
          clearFlags |= IMESH_MK_CAP_TUBE;
          break;
        case 'H':
          clearFlags |= IMESH_MK_CAP_DOME;
          break;
        case 'T':
          clearFlags |= IMESH_MK_FAST;
          break;
        case 's':
          clearFlags |= IMESH_MK_SKIP;
          break;
        case 'S':
          clearFlags |= IMESH_MK_SURF;
          break;
        case 'I':
          notimes = TRUE;
          break;
        case 'x':
          noxmin = TRUE;
          break;
        case 'y':
          noymin = TRUE;
          break;
        case 'z':
          nozmin = TRUE;
          break;
        case 't':
          notube = TRUE;
          break;
        default:
          fprintf(stderr, "%s: unknown option %s\n", progname, argv[iarg]);
          imodmesh_usage(progname, -1);
          break;
        }               
      }  

    } else
      break;
  }

  if (iarg >= argc)
    imodmesh_usage(progname, -1);
  if (replaceFileArgVec((const char ***)&argv, &argc, &iarg, &ob))
    exit(1);

  if ((setFlags & clearFlags) || (!times && notimes) || (ntube_list && notube)
      || (nocap && cap != IMESH_CAP_OFF) || (noSkipList && cap_skip_nz) ||
      (nozmin && minz != DEFAULT_VALUE) || 
      (noxmin && triMin.x > -DEFAULT_FLOAT) || 
      (noymin && triMin.y > -DEFAULT_FLOAT)) {
    fprintf(stderr, "%s: You cannot enter an option and its 'no' variation\n",
            progname);
    exit(1);
  }

  if (times) {
    setFlags |= IMESH_MK_TIME;
    clearFlags &= ~IMESH_MK_TIME;
  } else {
    clearFlags |= IMESH_MK_TIME;
    setFlags &= ~IMESH_MK_TIME;
  }

  if (newPolyNorm && (getenv(oldvar) != NULL)) {
    newPolyNorm = 0;
    printf("Making backward-compatible meshes because "
           "%s is set in environment\n", oldvar);
  }
  imeshSetNewPolyNorm(newPolyNorm);

  /* Skin multiple models in serial. */
  while (iarg < argc) {

    printf("Model: %s\n", argv[iarg]);
        
    imod = imodRead(argv[iarg]);
    if (!imod){
      fprintf(stderr, "%s: Error reading model %s\n", progname, argv[iarg]);
      exit(3);
    }

    spnt.x = imod->xscale;
    spnt.y = imod->yscale;
    spnt.z = imod->zscale * zscale;
    imodel_maxpt(imod, &max);

    if (erase_mesh)
      printf("Erasing mesh from objects.\n");

    for (ob = 0; ob < imod->objsize; ob++) {
      obj = &(imod->obj[ob]);
      if (!obj->contsize && !erase_mesh)
        continue;
      if (ObjOnList(ob, list, nlist)) {

        if (erase_mesh) {
          if (obj->meshsize) {
            imodMeshesDeleteZRange(&obj->mesh, &obj->meshsize,
                                   minz, maxz, resol);

            /* Turn off standard mesh view flags if they are all on and the
               regular resolution is being erased.  Have to modify views too */
            flags = IMOD_OBJFLAG_MESH | IMOD_OBJFLAG_NOLINE |
              IMOD_OBJFLAG_FILL;
            if (!resol && (obj->flags & flags) == flags)
              setOrClearFlags(&obj->flags, flags, 0);
            for (iv = 1; iv < imod->viewsize; iv++) {
              if (ob < imod->view[iv].objvsize) {
                obv = &imod->view[iv].objview[ob];
                if (!resol && (obv->flags & flags) == flags)
                  setOrClearFlags(&obv->flags, flags, 0);
              }
            }
          }
        } else if (renorm) {
          spnt.x = spnt.y = 1.0f;
          spnt.z = zscale;
          imodMeshesRescaleNormal(obj->mesh, &obj->meshsize, &spnt);
        } else if (remeshNorm) {
          if (obj->meshsize)
            obj->mesh = imeshReMeshNormal(obj->mesh, &obj->meshsize, &spnt,
                                          resol);
        } else {
          if (!iobjScat(obj->flags)) {

            useParam = (useOldParam && obj->meshParam) ? 1 : 0;
            if (!obj->meshParam) {
              obj->meshParam = imeshParamsNew();
              if (!obj->meshParam) {
                fprintf(stderr, "%s: Error creating new parameter structure\n",
                        progname);
                exit(3);
              }
            }
            param = obj->meshParam;

            if (useParam)
              param->flags = (param->flags & ~clearFlags) | setFlags;
            else
              param->flags = setFlags;

            if (!useParam || ntube_list > 0 || notube) {
              if (ntube_list > 0 && ObjOnList(ob, tube_list, ntube_list))
                param->flags |= IMESH_MK_TUBE;
              else
                param->flags &= ~IMESH_MK_TUBE;
            }

            if (!useParam || minz != DEFAULT_VALUE || nozmin) {
              param->minz = minz;
              param->maxz = maxz;
            }

            if (!useParam || incz >= 0) {
              if (resol) {
                param->inczLowRes = incz;
                if (param->inczLowRes < 0) {
                  param->inczLowRes = LOWRES_INCZ;
                  if (dowarn1)
                    printf("Setting Z increment to %d for low resolution "
                           "mesh\n" , param->inczLowRes);
                  dowarn1 = FALSE;
                }
              } else {
                param->inczHighRes = B3DMAX(1, incz);
              }
            }

            if (!useParam || tol >= 0) {
              if (resol) {
                param->tolLowRes = tol;
                if (param->tolLowRes < 0.) {
                  param->tolLowRes = LOWRES_TOL;
                  if (dowarn2)
                    printf("Setting tolerance to %.2f for low resolution "
                           "mesh\n" , param->tolLowRes);
                  dowarn2 = FALSE;
                }
              } else {
                param->tolHighRes = tol;
                if (param->tolHighRes < 0.) {
                  param->tolHighRes = DEFAULT_TOL;
                  if (dowarn3)
                    printf("Setting tolerance to %.2f for point "
                           "reduction\n", param->tolHighRes);
                  dowarn3 = FALSE;
                }
              }
            }

            if ((!useParam || tubeDiameter >= -10.) && 
                (param->flags & IMESH_MK_TUBE))
              param->tubeDiameter = B3DMAX(-2., tubeDiameter);

            if (!useParam || flatCrit >= 0)
              param->flatCrit = flatCrit;
            if (param->flatCrit < 0.)
              param->flatCrit = DEFAULT_FLAT;

            if (!useParam || overlap >= 0)
              param->overlap = B3DMAX(0., overlap);

            if (!useParam || passes >= 0)
              param->passes = B3DMAX(1, passes);

            if (!useParam || cap != IMESH_CAP_OFF || nocap)
              param->cap = cap;

            if (!useParam || triMin.x > -DEFAULT_FLOAT || noxmin) {
              param->xmin = triMin.x;
              param->xmax = triMax.x;
            }
            if (!useParam || triMin.y > -DEFAULT_FLOAT || noymin) {
              param->ymin = triMin.y;
              param->ymax = triMax.y;
            }
        
            if (!useParam || cap_skip_nz || noSkipList) {
              if (imeshCopySkipList(cap_skip_zlist, cap_skip_nz,
                                    &param->capSkipZlist, &param->capSkipNz))
                exit(3);

              /* If there is a list of capping exclusion Z values, check if 0
                 and model max are on list, and extend range by 1 */
              for (ilist = 0; ilist < param->capSkipNz; ilist++) {
                if (!param->capSkipZlist[ilist]) {
                  param->capSkipZlist[param->capSkipNz++] = -1;
                  break;
                }
              }
              for (ilist = 0; ilist < param->capSkipNz; ilist++) {
                if (param->capSkipZlist[ilist] == floor (max.z + 0.5) ) {
                  param->capSkipZlist[param->capSkipNz++] = floor(max.z + 1.5);
                  break;
                }
              }
            }

            tmshsize = 0;
            if (obj->meshsize){
              /* DNM: if appending, delete defined range 
                 and save mesh to add to later */
              if (append) {
                imodMeshesDeleteZRange(&obj->mesh, 
                                       &obj->meshsize,
                                       minz, maxz, resol);
              } else 
                imodMeshesDeleteRes(&obj->mesh,
                                    &obj->meshsize, resol);
              tmsh = obj->mesh;
              tmshsize = obj->meshsize;
              obj->meshsize = 0;
              obj->mesh = NULL;
            }
            printf("Meshing object  %d\n",ob+1);
            fflush(stdout);
            if (analyzePrepSkinObj(obj, resol, &spnt, NULL)) {
              fprintf(stderr, "%s: meshing error in model %s\n", progname,
                      argv[iarg]);
              exit(3);

            }

            /* set resolution flag now */
            for(m = 0; m < obj->meshsize; m++) {
              obj->mesh[m].flag |= (resol << IMESH_FLAG_RES_SHIFT);
              /*printf("mesh %d  vsize %d  lsize %d storesize %d\n", m, 
                     obj->mesh[m].vsize, obj->mesh[m].lsize, 
                     ilistSize(obj->mesh[m].store)); */
              /*istoreDump(obj->mesh[m].store);*/
            }

            if (tmshsize) {
              /* If appending, add new mesh to saved mesh 
                 then put it back into this object's mesh */
              for (m = 0; m < obj->meshsize; m++) {
                tmsh = imodel_mesh_add(&obj->mesh[m], tmsh, &tmshsize);
                if (!tmsh) {
                  fprintf(stderr, "%s: Error allocating new mesh array\n", progname);
                  exit(3);
                }
              }
              if (obj->meshsize)
                free(obj->mesh);
              if (append) 
                obj->mesh = imeshReMeshNormal(tmsh, &tmshsize, &spnt, resol);
              else
                obj->mesh = tmsh;
              obj->meshsize = tmshsize;
              /*for(m = 0; m < obj->meshsize; m++)
                 printf("mesh %d  vsize %d  lsize %d storesize %d\n", m, 
                     obj->mesh[m].vsize, obj->mesh[m].lsize, 
                     ilistSize(obj->mesh[m].store)); */
              /*if (m)istoreDump(obj->mesh[m].store);*/
            }

            /* Turn on standard mesh view flags if none of them are on, 
               otherwise turn on mesh flag */
            flags = IMOD_OBJFLAG_MESH | IMOD_OBJFLAG_NOLINE |
              IMOD_OBJFLAG_FILL;
            if (obj->flags & flags)
              setOrClearFlags(&obj->flags, IMOD_OBJFLAG_MESH, 1);
            else
              setOrClearFlags(&obj->flags, flags, 1);
            for (iv = 1; iv < imod->viewsize; iv++) {
              if (ob < imod->view[iv].objvsize) {
                obv = &imod->view[iv].objview[ob];
                if (obv->flags & flags)
                  setOrClearFlags(&obv->flags, IMOD_OBJFLAG_MESH, 1);
                else
                  setOrClearFlags(&obv->flags, flags, 1);
              }
            }
          }else
            printf("Skipping object %d\n", ob+1);
        }
      }

    }
        
    /* Save backup of Model to Model~ */
    if (imodBackupFile(argv[iarg])) {
      fprintf(stderr, "%s: Error, couldn't create backup", progname);
      exit(3);
    }

    fout = fopen(argv[iarg], "wb");
    if (!fout){
      fprintf(stderr, "%s: Error, couldn't open output.", progname);
      exit(3);
    }

    imodWrite(imod, fout);
    fclose(fout);
    imodDelete(imod);
    iarg++;
  }
  exit(0);
}


     
static void imodMeshesDeleteZRange(Imesh **meshp, 
                                   int *meshsize, 
                                   int minz, int maxz, int resol)
{
  int m;
  Imesh *mesh = *meshp;
  if ((minz == DEFAULT_VALUE) && (maxz == DEFAULT_VALUE)){
    imodMeshesDeleteRes(meshp, meshsize, resol);
    return;
  }

  for(m = 0; m < *meshsize; m++){
    if (imeshResol(mesh[m].flag) == resol)
      imodMeshDeleteZRange(&mesh[m], minz, maxz);
  }
  return;
}

static int deletez(Ipoint point, int minz, int maxz)
{
  int rmpoint = FALSE;
  int rmlow   = FALSE;
  int rmup    = FALSE;
  int z = point.z + 0.5f;

  if ((minz != DEFAULT_VALUE) && (z >= minz))
    rmlow = TRUE;
    
  if ((maxz != DEFAULT_VALUE) && (z <= maxz))
    rmup = TRUE;
  if ( (rmlow) && (rmup))
    rmpoint = TRUE;

  return(rmpoint);
}


/* Delete geometries in the range of minz and maxz.
 * The vertex data is not changed, only the index array
 * is changed.
 */
static void imodMeshDeleteZRange(Imesh *mesh, int minz, int maxz)
{
  int i;
  Ipoint pt[6];
  int newsize = 0;

  for(i = 0; i < mesh->lsize; i++){
        
    switch(mesh->list[i]){

    case IMOD_MESH_NORMAL:
      i++;
      break;


    case IMOD_MESH_BGNPOLYNORM:
      while(mesh->list[++i] != IMOD_MESH_ENDPOLY){
        if ((i + 7) > mesh->lsize) 
          break;
        pt[0] = mesh->vert[mesh->list[i++]];
        pt[1] = mesh->vert[mesh->list[i++]];
        pt[2] = mesh->vert[mesh->list[i++]];
        pt[3] = mesh->vert[mesh->list[i++]];
        pt[4] = mesh->vert[mesh->list[i++]];
        pt[5] = mesh->vert[mesh->list[i]];

        if ((deletez(pt[1], minz, maxz)) &&
            (deletez(pt[3], minz, maxz)) &&
            (deletez(pt[5], minz, maxz))
            ){
          i -= 5;
          mesh->list[i++] = DELETE_INDEX;
          mesh->list[i++] = DELETE_INDEX;
          mesh->list[i++] = DELETE_INDEX;
          mesh->list[i++] = DELETE_INDEX;
          mesh->list[i++] = DELETE_INDEX;
          mesh->list[i] = DELETE_INDEX;
        }
      }
      break;

    case IMOD_MESH_BGNPOLYNORM2:
      while(mesh->list[++i] != IMOD_MESH_ENDPOLY){
        if ((i + 4) > mesh->lsize) 
          break;
        pt[1] = mesh->vert[mesh->list[i++]];
        pt[3] = mesh->vert[mesh->list[i++]];
        pt[5] = mesh->vert[mesh->list[i]];

        if ((deletez(pt[1], minz, maxz)) &&
            (deletez(pt[3], minz, maxz)) &&
            (deletez(pt[5], minz, maxz))
            ){
          i -= 2;
          mesh->list[i++] = DELETE_INDEX;
          mesh->list[i++] = DELETE_INDEX;
          mesh->list[i] = DELETE_INDEX;
        }
      }
      break;


    case IMOD_MESH_END:
    case IMOD_MESH_SWAP:
    case IMOD_MESH_ENDTRI:
    case IMOD_MESH_BGNTRI:
      break;
            
    default:
      /*            if (deletez(mesh->vert[mesh->list[i]], minz, maxz))
                    mesh->list[i] = DELETE_INDEX;  seems like a bad idea */
      /*                imodMeshDeleteIndex(mesh, i); */

      break;

    }

  }

  for(i = 0; i < mesh->lsize; i++){
    if (mesh->list[i] != DELETE_INDEX)
      mesh->list[newsize++] = mesh->list[i];
  }
  /*    printf("original and new indices %d %d\n", mesh->lsize, newsize); */
  mesh->lsize = newsize;
}




static void imodMeshesRescaleNormal
(Imesh *inMesh, int *meshsize, Ipoint *spnt)
{
  int m;
  int i, j, ind, lsize;
  int listInc, vertBase, normAdd;
  Imesh *mesh;
    
  if ((!inMesh) || (*meshsize < 1)) return;


  for(m = 0; m < *meshsize; m++){
    mesh = &inMesh[m];
    lsize = mesh->lsize;
        
    for(i  = 0; i < lsize; i++){
      switch(mesh->list[i]){
      case IMOD_MESH_BGNPOLY:
      case IMOD_MESH_BGNBIGPOLY:
        while(mesh->list[++i] != IMOD_MESH_ENDPOLY);
        break;
                
      case IMOD_MESH_NORMAL:
        i++;
        if ((mesh->list[i] < mesh->vsize) && (mesh->list[i] > -1)){
          mesh->vert[mesh->list[i]].x *= spnt->x;
          mesh->vert[mesh->list[i]].y *= spnt->y;
          mesh->vert[mesh->list[i]].z *= spnt->z;
          imodPointNormalize(&mesh->vert[mesh->list[i]]);
        }
        break;
                
      case IMOD_MESH_BGNPOLYNORM:
      case IMOD_MESH_BGNPOLYNORM2:
        imodMeshPolyNormFactors(mesh->list[i++], &listInc, &vertBase,
                                &normAdd);
        while (mesh->list[i] != IMOD_MESH_ENDPOLY) {
          for (j = 0; j < 3; j++) {
            ind = mesh->list[i] + normAdd;
            mesh->vert[ind].x *= spnt->x;
            mesh->vert[ind].y *= spnt->y;
            mesh->vert[ind].z *= spnt->z;
            imodPointNormalize(&mesh->vert[ind]);
            i += listInc;
          }
        }
        break;
                
      case IMOD_MESH_BGNTRI:
      case IMOD_MESH_ENDTRI:
      case IMOD_MESH_SWAP:
        break;
                
      default:
        break;
      }
    }
  }
  return;
}


static int ObjOnList(int ob, int list[], int nlist)
{
  int i;
  if (!nlist)
    return 1;
  else
    for (i = 0; i < nlist; i++) {
      if (ob + 1 == list[i])
        return 1;
    }
  return 0;
}
