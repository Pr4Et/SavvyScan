! $Id$
module mtkvars
  integer LIMVERTS, LIMTRIANG, LIMPOLY, LIM_MESH_OBJ, LIMSURF, LIMCONT, LIMSHIFT
  integer  LIM_MESH_IND, LIMBINSAVE, LIMXYZ, LIMSIZES
  parameter (LIMVERTS = 1800000, LIMTRIANG = LIMVERTS * 2, LIMPOLY = 150000)
  parameter (LIM_MESH_OBJ = 255, LIMSURF = 30000, LIMCONT = LIMPOLY, LIMXYZ = 500000)
  parameter (LIMSHIFT = 20000, LIMSIZES = LIMXYZ, LIMGRAPHS = 50, LIM_WIMP_OBJ = 90000)
  parameter (LIM_MESH_IND = LIMVERTS * 6, LIMBINSAVE = LIM_MESH_IND * 2)
  !       Old notes on memory usage
  !       limverts * 12 + limtri * 80 = limverts * 172
  !       limpoly * 40 + limconts * 24 = limpoly * 64
  !       limsurf * 44 = 0.4
  !       limshift * 16 = 0.3
  real*4, allocatable :: verts(:,:)                   !vertex array
  integer*4, allocatable :: indVert(:,:)              !index to 3 vertices of triangl
  real*4, allocatable :: triXYrot(:,:,:)              !rotated X and Y coords
  real*4, allocatable :: triZrot(:)                   !rotated Z coord
  real*4, allocatable :: cosBeta(:), sinBeta(:)       !cosine and sine of beta
  real*4, allocatable :: cosGamma(:), sinGamma(:)     !cosine and sine of gamma
  real*4, allocatable :: triangXmin(:), triangXmax(:)
  real*4, allocatable :: triangYmin(:), triangYmax(:)
  real*4, allocatable :: triangZmin(:), triangZmax(:)
  real*4 polyXmin(LIMPOLY), polyXmax(LIMPOLY)
  real*4 polyYmin(LIMPOLY), polyYmax(LIMPOLY)
  real*4 polyZmin(LIMPOLY), polyZmax(LIMPOLY)
  real*4 surfXmin(LIMSURF), surfXmax(LIMSURF)
  real*4 surfYmin(LIMSURF), surfYmax(LIMSURF)
  real*4 surfZmin(LIMSURF), surfZmax(LIMSURF)
  real*4 contXmin(LIMCONT), contXmax(LIMCONT)
  real*4 contYmin(LIMCONT), contYmax(LIMCONT)
  real*4 contZval(LIMCONT)
  real*4 polyArea(LIMPOLY), surfArea(LIMSURF)
  real*4 shifts(3,LIMSHIFT)
  logical*1 shifted(LIMSHIFT)
  integer*4 numInPoly(LIMPOLY)                !# of triangles in polygon
  integer*4 indStartPoly(LIMPOLY)             !starting triangle on Z-plane
  integer*4 iobjMesh(LIM_MESH_OBJ)            !objects with meshes loaded
  integer*4 iobjPoly(LIM_MESH_OBJ)            !starting polygon of object
  integer*4 numPolyObj(LIM_MESH_OBJ)          !# of polygons in object
  integer*4 listSurf(LIMPOLY)                 !list of polygons by surface
  integer*4 numInSurf(LIMSURF)                !# of polygons in surface
  integer*4 indStartSurf(LIMSURF)             !starting polygon of surface
  integer*4 iobjSurf(LIM_MESH_OBJ)            !starting surface of object
  integer*4 numSurfObj(LIM_MESH_OBJ)          !# of surfaces in object
  integer*4 listCont(LIMCONT)                 !list of contours by surface
  integer*4 indStartCont(LIMSURF)             !starting contour of surface
  integer*4 numContInSurf(LIMSURF)            !# of contours in surface
  integer*4 iobjShift(LIM_MESH_OBJ)           !list of shifted objects
  integer*4 numItemShifted(LIM_MESH_OBJ)      !# of items shifted in object
  integer*4 indStartShift(LIM_MESH_OBJ)       !index to first shift of object
  integer*4 numVerts                          !# of vertices loaded
  integer*4 numTriang                         !# of triangles
  integer*4 numPoly                           !# of polygons
  integer*4 numSurf                           !# of surfaces
  integer*4 numConts                          !# of contours
  integer*4 numMeshLoaded                     !# of meshes loaded
  integer*4 numObjShifted                     !# of objects shifted
  integer*2, allocatable :: ibinSave(:)
  real*4, allocatable :: aaLine(:), bbLine(:), ccLine(:), ddLine(:)
  real*4, allocatable :: xMin(:), xMax(:), yMin(:), yMax(:)
  real*4, allocatable :: zMin(:), zMax(:)
  real*4, allocatable :: sizes(:)
  real*4 globalXmin(LIM_WIMP_OBJ), globalYmin(LIM_WIMP_OBJ), globalZmin(LIM_WIMP_OBJ)
  real*4 globalXmax(LIM_WIMP_OBJ), globalYmax(LIM_WIMP_OBJ), globalZmax(LIM_WIMP_OBJ)
  logical*1 isNeighPt(LIMGRAPHS,LIM_WIMP_OBJ)
  logical*1 needShift(LIMXYZ), needCheck(LIMXYZ)
  integer*2 numTrials(LIMXYZ)

end module mtkvars
