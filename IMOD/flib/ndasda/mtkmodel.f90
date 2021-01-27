! $Id$
!
! READ_MODEL reads a WIMP model.
! if MODELFILE is blank, it requests a model file name; otherwise it
! attempts to open the file specified by that name.
! If there is an error, it loops back on the request for a file name.
! If a new file name is needed, it then asks for a file of tilt info,
! then requests scaling parameters
! for this file, and returns them in XYSCAL and ZSCAL.
! The file names are returned in MODELFILE and TILTFILE.
!
!
subroutine read_model(modelFile, tiltFile, xyScale, zTotScale, xOffset, yOffset, &
    zOffset, ifFlip, iobjFlag, limFlag, zGapStart, zGapEnd, numGaps)
  use mtkvars
  implicit none
  include 'model.inc90'
  character*(*) modelFile, tiltFile
  integer*4 iobjFlag(*)
  logical exist, readw_or_imod, newFile
  real*4 tiltZstart(1000), tilt(1000), cosTilt(1000), remapZ(1000)
  integer*4 listGap(1000)
  real*4 zGapStart(*), zGapEnd(*)
  integer*4 i, icur, ierr, ierr2, ifFlip, itilt, limFlag, numGaps, numListGap, numTilts
  real*4 dfltScale, fTemp, secThick, umPerPixel, xImScale, xMag, xOffset, xTemp, xyScale
  real*4 yImScale, yOffset, yTemp, zImScale, zNew, zOffset, zScale, zTemp, zTotScale, zz
  integer*4 getImodHead, getImodFlags, getImodScales
  real*4 cosd
  integer*4 in5
  common /nmsInput/ in5
  !
  newFile = .false.
  if (modelFile .ne. ' ') go to 92
91 print *,'Enter name of input model file, or Return to skip', &
      ' to entering options'
  read(in5, '(a)') modelFile
  if (modelFile == ' ') return
  newFile = .true.
  !
92 exist = readw_or_imod(modelFile)
  if (.not.exist) go to 91
  !
  if (newFile) then
    write(*,'(1x,a,$)') &
        'Name of file of tilt info (Return if none): '
    read(in5, '(a)') tiltFile
  endif
  !
  if (tiltFile .ne. ' ') then
    call dopen(3, tiltFile, 'ro', 'f')
    numTilts = 0
3   i = numTilts + 1
    read(3,*,end = 5) tiltZstart(i), tilt(i)
    numTilts = i
    cosTilt(numTilts) = cosd(tilt(numTilts))
    go to 3
5   remapZ(1) = tiltZstart(1)
    do i = 1, numTilts - 1
      remapZ(i + 1) = remapZ(i) + (tiltZstart(i + 1) - tiltZstart(i)) &
          / cosTilt(i)
    enddo
    zMax = -1.e10
    do i = 1, n_point
      zz = p_coord(3, i)
      if (zz >= tiltZstart(1)) then
        itilt = numTilts
        do while(zz < tiltZstart(itilt))
          itilt = itilt - 1
        enddo
        zNew = remapZ(itilt) + (zz - tiltZstart(itilt)) / &
            cosTilt(itilt)
        p_coord(3, i) = zNew
        zMax = max(zMax, zNew)
      endif
    enddo
    print *,'new maximum z is', zMax
  endif
  !
  ! DNM 10/7/04: Better get scale and offset regardless of new or not
  !
  ierr = getImodHead(xyScale, zScale, xOffset, yOffset, zOffset, ifFlip)
  ierr2 = getImodScales(xImScale, yImScale, zImScale)
  if (newFile) then
    dfltScale = 1.e6
    if (ierr == 0 .and. abs(xyScale - dfltScale) / dfltScale > 1.e-5) then
      write(*,'(a,f10.6,a)') ' Scale set from model header at', &
          xyScale, ' microns/pixel'
      zTotScale = xyScale * zScale
    else
      ifFlip = 0
      write(*,'(1x,a,$)') 'Magnification of negatives: '
      read(in5,*) xMag
      write(*,'(1x,a,$)') 'Scale at which negatives were digitized' &
          //' (microns per pixel): '
      read(in5,*) umPerPixel
      write(*,'(1x,a,$)') 'Nominal section thickness (nm): '
      read(in5,*) secThick
      !
      xyScale = umPerPixel / xMag
      zTotScale = secThick / 1000.
    endif
  endif
  !
  ! flip, shift back to imod coordinates, and scale points
  !
  do i = 1, n_point
    xTemp = (p_coord(1, i) - xOffset) / xImScale
    yTemp = (p_coord(2, i) - yOffset) / yImScale
    zTemp = (p_coord(3, i) - zOffset) / zImScale
    if (ifFlip .ne. 0) then
      fTemp = yTemp
      yTemp = zTemp
      zTemp = fTemp
    endif
    p_coord(1, i) = xTemp * xyScale
    p_coord(2, i) = yTemp * xyScale
    p_coord(3, i) = zTemp * zTotScale
  enddo
  !
  ! get the object types
  !
  do i = 1, limFlag
    iobjFlag(i) = 1
  enddo
  ierr = getImodFlags(iobjFlag, limFlag)
  if (ierr .ne. 0) then
    print *,'Error getting object types; assuming all objects '// &
        'are open contours'
  else
    !
    ! strip mesh flag off of 1's and 2's
    !
    do i = 1, limFlag
      if (iobjFlag(i) > 4) iobjFlag(i) = iobjFlag(i) - 4
    enddo
  endif
  !
  ! get list of z gaps to connect surfaces across
  !
  if (newFile) then
    print *,'Enter list of Z values (numbered from 0) across ', &
        'which surfaces should be connected (Return for none)'
    call rdlist(in5, listGap, numListGap)
    numGaps = 0
    if (numListGap > 0) then
      !
      ! really dumb, but convert the list back to start and end values
      !
      icur = 1
      numGaps = 1
      zGapStart(1) = listGap(1) * zTotScale
      do while(icur <= numListGap)
        if (icur == numListGap .or. listGap(icur) + 1 .ne. listGap(icur + 1)) &
            then
          zGapEnd(numGaps) = listGap(icur) * zTotScale
          icur = icur + 1
          if (icur <= numListGap) then
            numGaps = numGaps + 1
            zGapStart(numGaps) = listGap(icur) * zTotScale
          endif
        else
          icur = icur + 1
        endif
      enddo
    endif
  endif
  !
  ! initialize # of loaded meshes and number of shifted objects
  !
  numMeshLoaded = 0
  numObjShifted = 0
  return
end subroutine read_model

!
! return open contours and scattered points in the "mt" arrays
!
subroutine get_objects(zStart, zEnd, xModPt, yModPt, zModPt, indStart, numPtInObj, &
    itype, numWobj, iobjImod, iobjFlag, limxyz, limwobj)
  implicit none
  include 'model.inc90'
  real*4 xModPt(*), yModPt(*), zModPt(*)
  integer*4 indStart(*), numPtInObj(*), itype(*), iobjImod(*), iobjFlag(*)
  integer*4 ibase, iflag, imodObj, indFree, iobj, ip, ipEnd, ipStart, ipt, limwobj
  integer*4 limxyz, numWobj
  real*4 zEnd, zStart
  !
  indFree = 1
  numWobj = 0
  !
  do iobj = 1, max_mod_obj
    imodObj = 256 - obj_color(2, iobj)
    if (imodObj <= 0) call errorExit('TOO MANY MODEL OBJECTS FOR ARRAYS')
    iflag = iobjFlag(imodObj)
    if (npt_in_obj(iobj) > 0) then
      ibase = ibase_obj(iobj)
      !
      ! keep "coplanar" contours now (huh?) ; if z limits are set, trim
      ! open contours until start and end within limits
      !
      ipStart = 1
      ipEnd = npt_in_obj(iobj)
      if (iflag == 1 .and. (zStart .ne. 0. .or. zEnd .ne. 0.)) then
        do while (((zEnd .ne. 0 .and. &
            p_coord(3, abs(object(ipEnd + ibase))) > zEnd) &
            .or. (zStart .ne. 0 .and. &
            p_coord(3, abs(object(ipEnd + ibase))) < zStart)))
          ipEnd = ipEnd-1
          if (ipEnd <= 0) exit
        enddo
        do while (ipStart < ipEnd .and. &
            ((zEnd .ne. 0 .and. p_coord(3, abs(object(ipStart + ibase))) > zEnd) &
            .or. (zStart .ne. 0 .and. &
            p_coord(3, abs(object(ipStart + ibase))) < zStart)))
          ipStart = ipStart + 1
        enddo
      endif
      if (iflag == 2 .or. (iflag == 1 .and. ipEnd + 1 - ipStart >= 2)) then
        numWobj = numWobj + 1
        if (numWobj > limwobj) call errorExit( &
            'TOO MANY LINE & SCATTERED POINT CONTOURS FOR ARRAYS')
        indStart(numWobj) = indFree
        numPtInObj(numWobj) = ipEnd + 1 - ipStart
        iobjImod(numWobj) = iobj
        itype(numWobj) = imodObj
        if (obj_color(1, iobj) == 0) itype(numWobj) = -itype(numWobj)
        if (indFree + ipEnd - ipStart >= limxyz) call errorExit( &
            'TOO MANY POINTS FOR LINE/SCATTERED POINT ARRAYS')
        do ipt = ipStart, ipEnd
          ip = abs(object(ipt + ibase))
          xModPt(indFree) = p_coord(1, ip)
          yModPt(indFree) = p_coord(2, ip)
          zModPt(indFree) = p_coord(3, ip)
          indFree = indFree + 1
        enddo
      endif
    endif
  enddo
  return
end subroutine get_objects

!
!
! SAVE_MODEL saves the current model after unscaling it.  If NINWIN
! is nonzero, it can also save new objects marking the points of
! closest approach
!
subroutine save_model(xyScale, zTotScale, xOffset, yOffset, zOffset, ifFlip, &
    iobjFlag, xModPt, yModPt, zModPt, indStart, &
    numPtInObj, itype, numWobj, numInWin, iobjInWin, numObjInWin, iobjImod, endSeparation)
  use mtkvars
  implicit none
  include 'model.inc90'
  integer LIMCHANGE, LIM_CONT_SIZES
  parameter (LIMCHANGE = 2000, LIM_CONT_SIZES = 10000)
  real*4 xModPt(*), yModPt(*), zModPt(*), endSeparation(*), changeLow(LIMCHANGE)
  real*4 changeHigh(LIMCHANGE)
  integer*4 indStart(*), numPtInObj(*), itype(*), iobjInWin(*), iobjImod(*)
  real*4 xyScale, zTotScale, xOffset, yOffset, zOffset, fTemp, xTemp, yTemp, zTemp
  real*4 xImScale, yImScale, zImScale, contSizes(LIM_CONT_SIZES)
  integer*4 ifFlip, numWobj, numInWin, numObjInWin, ibaseScat, locForObj, indCont, isurf
  integer*4 itypeOld(LIMCHANGE), itypeNew(LIMCHANGE), iobjFlag(*), itypeUsed(LIMCHANGE)
  integer*4 i, lastObject, numTypeChange, itypeConLine, itypeScat, ierr2, iflag, iobj
  integer*4 itypeTemp, maybeNew, ityp, itypeInd, iow, ibase, ipt, ii, jj, ibaseConLine
  integer*4 iobjNew, numSizes, imodObj, imodCont, indMesh, iorig
  character*120 lastModel
  integer*4 getImodObjSize, getImodScales
  integer*4 getContPointSizes, putContPointSizes
  integer*4 in5
  common /nmsInput/ in5
  !
95 write(*,'(1x,a,$)') 'Name of output model file: '
  read(in5, '(a)') lastModel
  !
  lastObject = getImodObjSize()
  numTypeChange = 0
  if (numObjInWin .ne. 0) then
    print *,'Enter list of objects for which to put', &
        'contours in window into new objects (Return for none)'
    call rdlist(in5, itypeOld, numTypeChange)
    do i = 1, numTypeChange
      itypeNew(i) = lastObject + i
      itypeUsed(i) = 0
    enddo
    if (numTypeChange .ne. 0 .and. numObjInWin < 0) then
      print *,'Enter lower and upper limits of end separation', &
          ' for each object change'
      read(in5,*) (changeLow(i), changeHigh(i), i = 1, numTypeChange)
    endif
  endif
  itypeConLine = lastObject + numTypeChange + 1
  itypeScat = lastObject + numTypeChange + 2
  !
  ! zero out the model first; all non-mesh objects
  !
  ierr2 = getImodScales(xImScale, yImScale, zImScale)
  do i = 1, max_mod_obj
    if (npt_in_obj(i) > 0) then
      iflag = iobjFlag(256 - obj_color(2, i))
      if (iflag == 1 .or. iflag == 2) npt_in_obj(i) = 0
    endif
  enddo
  !
  ! save each separate object.  Keep track of which new objects are
  ! being accessed
  !
  do iobj = 1, numWobj
    iorig = iobjImod(iobj)
    ! ibase_obj(iobj) =indstrt(iobj) -1
    itypeTemp = itype(iobj)
    if (numTypeChange .ne. 0) then
      maybeNew = -999
      if (numObjInWin >= 0) then
        !
        ! First see if this contour is in an object being changed
        !
        do ityp = 1, numTypeChange
          if (itypeTemp == itypeOld(ityp)) then
            maybeNew = itypeNew(ityp)
            itypeInd = ityp
          endif
        enddo
        !
        ! Then see if this contour is in the window (positive iobjwin)
        !
        if (maybeNew .ne. -999) then
          do iow = 1, numObjInWin
            if (iobj == iobjInWin(iow)) then
              itypeTemp = maybeNew
              itypeUsed(itypeInd) = 1
            endif
          enddo
        endif
      else
        !
        ! For dealing with end separations
        !
        do iow = 1, -numObjInWin
          if (iobjImod(iobj) == iobjInWin(iow)) then
            do ityp = 1, numTypeChange
              if (itypeTemp == itypeOld(ityp) .and. endSeparation(iow) >= &
                  changeLow(ityp) .and. endSeparation(iow) <= changeHigh(ityp)) then
                maybeNew = itypeNew(ityp)
                itypeUsed(ityp) = 1
              endif
            enddo
          endif
        enddo
        if (maybeNew .ne. -999) itypeTemp = maybeNew
      endif
    endif
    !
    ! reassign contour to new object by making a new contour with the
    ! new color value, leaving old one empty, and copy the
    ! contour's points in regardless of whether it is new or old
    ! Also transfer point sizes to a new contour
    !
    iobjNew = iorig
    if (itypeTemp .ne. itype(iobj)) then
      max_mod_obj = max_mod_obj + 1
      iobjNew = max_mod_obj
      obj_color(2, iobjNew) = 256 - abs(itypeTemp)
      call objToCont(iorig, obj_color, imodObj, imodCont)
      ierr2 = getContPointSizes(imodObj, imodCont, contSizes, LIM_CONT_SIZES, &
          numSizes)
      if (ierr2 == 0 .and. numSizes > 0) then
        call objToCont(iobjNew, obj_color, imodObj, imodCont)
        ierr2 = putContPointSizes(imodObj, imodCont, contSizes, numSizes) ;
      endif
    endif
    obj_color(2, iobjNew) = 256 - abs(itypeTemp)
    npt_in_obj(iobjNew) = numPtInObj(iobj)
    if (itypeTemp >= 0) then
      obj_color(1, iobjNew) = 1
    else
      obj_color(1, iobjNew) = 0
    endif
    ibase = ibase_obj(iorig)
    ibase_obj(iobjNew) = ibase

    do ipt = 1, numPtInObj(iobj)
      ii = indStart(iobj) + ipt - 1
      jj = object(ibase + ipt)
      p_coord(1, jj) = xModPt(ii)
      p_coord(2, jj) = yModPt(ii)
      p_coord(3, jj) = zModPt(ii)
    enddo
  enddo
  !
  ! Next look at surfaces in window (negative iobjwin) and find the objects
  ! they originate from
  !
  if (numTypeChange > 0) then
    do iow = 1, numObjInWin
      isurf = -iobjInWin(iow)
      if (isurf > 0) then
        itypeTemp = 0
        do indMesh = 1, numMeshLoaded
          if (isurf >= iobjSurf(indMesh) .and. isurf < &
              iobjSurf(indMesh) + numSurfObj(indMesh)) then
            !
            ! See if this object is on change list
            !
            itypeInd = -1
            do ityp = 1, numTypeChange
              if (iobjMesh(indMesh) == itypeOld(ityp)) itypeInd = ityp
            enddo
            if (itypeInd > 0) then
              !
              ! change all the contours in the surface: again, make a new
              ! contour, assign points to it and zero out old one, and
              ! transfer point sizes
              !
              itypeUsed(itypeInd) = 1
              do indCont = indStartCont(isurf), indStartCont(isurf) +  &
                  numContInSurf(isurf) - 1
                iobj = listCont(indCont)
                max_mod_obj = max_mod_obj + 1
                iobjNew = max_mod_obj
                npt_in_obj(iobjNew) = npt_in_obj(iobj)
                npt_in_obj(iobj) = 0
                ibase_obj(iobjNew) = ibase_obj(iobj)
                obj_color(2, iobjNew) = 256 - itypeNew(itypeInd)
                obj_color(1, iobjNew) = 1
                call objToCont(iobj, obj_color, imodObj, imodCont)
                ierr2 = getContPointSizes(imodObj, imodCont, contSizes, &
                    LIM_CONT_SIZES, numSizes)
                if (ierr2 == 0 .and. numSizes > 0) then
                  call objToCont(iobjNew, obj_color, imodObj, imodCont)
                  ierr2 = putContPointSizes(imodObj, imodCont, contSizes, &
                      numSizes) ;
                endif
              enddo
            endif
          endif
        enddo
      endif
    enddo
  endif
  !
  ! create connectors
  !
  n_point = ibase_free + 3 * numInWin
  do i = 1, 3 * numInWin
    jj = ibase_free + i
    ii = i
    if (numWobj > 0) ii = indStart(numWobj) + numPtInObj(numWobj) + i - 1
    object(jj) = jj
    p_coord(1, jj) = xModPt(ii)
    p_coord(2, jj) = yModPt(ii)
    p_coord(3, jj) = zModPt(ii)
  enddo
  !
  ! transfer object flag to new objects if any
  !
  do ityp = 1, numTypeChange
    if (itypeUsed(ityp) .ne. 0) then
      lastObject = lastObject + 1
      iflag = iobjFlag(itypeOld(ityp))
      if (iflag == 4) iflag = 0
      call putImodFlag(lastObject, iflag)
    endif
  enddo
  !
  ! save connecting lines if any
  !
  ibaseConLine = ibase_free
  ibaseScat = ibase_free + 2 * numInWin
  do i = 1, numInWin
    iobj = max_mod_obj + i
    npt_in_obj(iobj) = 2
    obj_color(1, iobj) = 1
    obj_color(2, iobj) = 256 - itypeConLine
    locForObj = ibase_free + 1 + 3 * (i - 1)
    ibase_obj(iobj) = ibaseConLine + 2 * (i - 1)
    object(ibase_obj(iobj) + 1) = locForObj
    object(ibase_obj(iobj) + 2) = locForObj + 2
    object(ibaseScat + i) = locForObj + 1
  enddo
  if (numInWin > 0) then
    n_object = n_object + numInWin + 1
    max_mod_obj = max_mod_obj + numInWin + 1
    npt_in_obj(max_mod_obj) = numInWin
    ibase_obj(max_mod_obj) = ibaseScat
    obj_color(1, max_mod_obj) = 1
    obj_color(2, max_mod_obj) = 256 - itypeScat

    call putImodFlag(lastObject + 1, 1)
    lastObject = lastObject + 2
    call putImodFlag(lastObject, 2)
    call putScatSize(lastObject, 5)
  endif
  print *,'There are', lastObject, ' objects in the output model'
  print *,'Connecting lines and points are in objects', lastObject - 1, &
      ' and', lastObject
  !
  ! unscale, unflip
  !
  do i = 1, n_point
    xTemp = p_coord(1, i) / xyScale
    yTemp = p_coord(2, i) / xyScale
    zTemp = p_coord(3, i) / zTotScale
    if (ifFlip .ne. 0) then
      fTemp = yTemp
      yTemp = zTemp
      zTemp = fTemp
    endif
    p_coord(1, i) = xImScale * xTemp + xOffset
    p_coord(2, i) = yImScale * yTemp + yOffset
    p_coord(3, i) = zImScale * zTemp + zOffset
  enddo
  !
  call write_wmod(lastModel)
  close(20)
  return
end subroutine save_model


! PROCESS_MESH reads in the mesh for object IMODOBJ, scales coordinates
! by XYSCAL in X/Y and ZSCAL in Z, and fills arrays with indices to
! triangles and polygons.  It determines the sines and cosines needed
! to rotate each triangle to a plane, and stores coordinates for
! rotated triangles.  It then combines adjacent, overlapping polygons
! into surfaces, continuing a surface across gaps if any are specified.
! MODIND is for temporary storage of indices
! ZGAPST and ZGAPND have starting and ending Z values of NGAP gaps.
! IFERR returns 1 if there is an error
! Diagnostic reports and error messages are suppressed if IFHUSH = 1

subroutine process_mesh(imodObj, xyScale, zTotScale, modIndTemp, zGapStart, &
    zGapEnd, numGaps, ifErr, ifHush)
  use mtkvars
  implicit none
  include 'model.inc90'
  integer LIMCXY, LIMIPMIN
  parameter (LIMCXY = 2000, LIMIPMIN = 100)
  integer*4 modIndTemp(LIM_MESH_IND), ipolyMin(LIMIPMIN)
  real*4 aVector(3), bVector(3), cVector(3), zGapStart(*), zGapEnd(*)
  real*4 cont1X(LIMCXY), cont2X(LIMCXY), cont1Y(LIMCXY), cont2Y(LIMCXY)
  real*4 tempXmin(LIMCONT), tempXmax(LIMCONT)
  real*4 tempYmin(LIMCONT), tempYmax(LIMCONT)
  real*4 tempZval(LIMCONT)

  integer*4 i, ibase, ic1, ic2, icontSurf, idoneCode, iendCode, ierr, ifErr
  integer*4 ifFound, ifHush, ifMostIn, ifOk, ifZeroedShift, ii, ilist, imodObj, ind, ind1
  integer*4 ind2, ind3, indCont, indGap, indMesh, iobj, ip, ipCur, ipoint, ipoly
  integer*4 ipt, ishift, istartCode, isum, isurf, itest, itr, itriang, itypeObj
  integer*4 iv, ivert, jlist, jnd1, jnd2, jnd3, jpoly, jtr, list, listCur, listEnd
  integer*4 newInd, newVerts, numAddCont, numDupTriang, numGaps, numInCont1
  integer*4 numInCont2, numInObj, numNonTriang, numOrphan, numPolyMin, numTemp
  real*4 cLenInXY, cLength, cosBet, cosGam, delZ, delZmin, sinBet, sinGam, toler, txmax
  real*4 txmin, tymax, tymin, tzmax, tzmin, tzval, xyScale, zExtract1, zExtract2, zIncr
  real*4 zRotSum, zTotScale
  logical anyinside, inside
  integer*4 getImodVerts

  istartCode = -23
  iendCode = -22
  idoneCode = -1
  ifErr = 1
  numDupTriang = 0
  numNonTriang = 0
  ierr = getImodVerts(imodObj, verts(1, numVerts + 1), &
      modIndTemp, LIMVERTS - numVerts, LIM_MESH_IND, newVerts, newInd)
  if (ierr .ne. 0) then
    if (ifHush == 0) print *,'Error or insufficient space '// &
        'loading mesh for object', imodObj
    return
  endif
  !
  ! scale the vertices
  !
  do iv = numVerts + 1, numVerts + newVerts
    verts(1, iv) = verts(1, iv) * xyScale
    verts(2, iv) = verts(2, iv) * xyScale
    verts(3, iv) = verts(3, iv) * zTotScale
  enddo
  !
  ! process polygons, get indices to triangles
  !
  numMeshLoaded = numMeshLoaded + 1
  iobjMesh(numMeshLoaded) = imodObj
  iobjPoly(numMeshLoaded) = numPoly + 1
  ilist = 1
  do while(modIndTemp(ilist) .ne. idoneCode)
    !
    ! advance to start of a polygon
    !
    do while(modIndTemp(ilist) .ne. istartCode .and. &
        modIndTemp(ilist) .ne. idoneCode)
      ilist = ilist + 1
    enddo
    if (modIndTemp(ilist) == istartCode) then
      ilist = ilist + 1
      if (numPoly == LIMPOLY) then
        if (ifHush == 0) print *,'Too many polygons for arrays'
        return
      endif
      numPoly = numPoly + 1
      numInPoly(numPoly) = 0
      polyXmin(numPoly) = 1.e10
      polyYmin(numPoly) = 1.e10
      polyZmin(numPoly) = 1.e10
      polyXmax(numPoly) = -1.e10
      polyYmax(numPoly) = -1.e10
      polyZmax(numPoly) = -1.e10
      polyArea(numPoly) = 0.
      do while (modIndTemp(ilist) .ne. iendCode)
        !
        ! start a new triangle: first check that its a triangle
        ! shift the indices up by amount already in verts array
        ! plus one because they are C indices
        !
        ind1 = modIndTemp(ilist) + numVerts + 1
        ind2 = modIndTemp(ilist + 1) + numVerts + 1
        ind3 = modIndTemp(ilist + 2) + numVerts + 1
        ifOk = 1
        if (ind1 == ind2 .or. ind2 == ind3 .or. ind1 == ind3) then
          ! print *,'NON-TRIANGLE ignored', npoly, ntriang
          numNonTriang = numNonTriang + 1
          ifOk = 0
        endif
        !
        ! check for duplicate triangles in this polygon
        !
        if (ifOk == 1 .and. numInPoly(numPoly) > 0) then
          isum = ind1 + ind2 + ind3
          do itriang = indStartPoly(numPoly), numTriang
            jnd1 = indVert(1, itriang)
            jnd2 = indVert(2, itriang)
            jnd3 = indVert(3, itriang)
            if (jnd1 + jnd2 + jnd3 == isum) then
              if ((ind1 == jnd1 .and. (ind2 == jnd2 .or. ind2 == jnd3)) .or. &
                  (ind1 == jnd2 .and. (ind2 == jnd1 .or. ind2 == jnd3)) .or. &
                  (ind1 == jnd3 .and. (ind2 == jnd1 .or. ind2 == jnd2))) &
                  then
                ifOk = 0
                numDupTriang = numDupTriang + 1
                ! print *,'DUPLICATE TRIANGLE ignored', &
                ! npoly, ntriang, itri
              endif
            endif
          enddo
        endif
        !
        if (ifOk .ne. 0) then
          !
          ! save the triangle indices
          !
          if (numTriang == LIMTRIANG) then
            if (ifHush == 0) print *,'Too many triangles in mesh' &
                //' to fit in arrays'
            return
          endif
          numTriang = numTriang + 1
          if (numInPoly(numPoly) == 0) indStartPoly(numPoly) = numTriang
          numInPoly(numPoly) = numInPoly(numPoly) + 1
          indVert(1, numTriang) = ind1
          indVert(2, numTriang) = ind2
          indVert(3, numTriang) = ind3
          !
          ! get the bounding boxes
          !
          triangXmin(numTriang) = min(verts(1, ind1), verts(1, ind2), &
              verts(1, ind3))
          triangYmin(numTriang) = min(verts(2, ind1), verts(2, ind2), &
              verts(2, ind3))
          triangZmin(numTriang) = min(verts(3, ind1), verts(3, ind2), &
              verts(3, ind3))
          triangXmax(numTriang) = max(verts(1, ind1), verts(1, ind2), &
              verts(1, ind3))
          triangYmax(numTriang) = max(verts(2, ind1), verts(2, ind2), &
              verts(2, ind3))
          triangZmax(numTriang) = max(verts(3, ind1), verts(3, ind2), &
              verts(3, ind3))
          polyXmin(numPoly) = min(polyXmin(numPoly), triangXmin(numTriang))
          polyYmin(numPoly) = min(polyYmin(numPoly), triangYmin(numTriang))
          polyZmin(numPoly) = min(polyZmin(numPoly), triangZmin(numTriang))
          polyXmax(numPoly) = max(polyXmax(numPoly), triangXmax(numTriang))
          polyYmax(numPoly) = max(polyYmax(numPoly), triangYmax(numTriang))
          polyZmax(numPoly) = max(polyZmax(numPoly), triangZmax(numTriang))
          !
          ! compute rotation to a plane and rotated coordinates
          !
          do i = 1, 3
            aVector(i) = verts(i, ind2) - verts(i, ind1)
            bVector(i) = verts(i, ind3) - verts(i, ind2)
          enddo
          call crossproduct(aVector, bVector, cVector)
          cLenInXY = sqrt(cVector(1)**2 + cVector(2)**2)
          cLength = sqrt(cVector(1)**2 + cVector(2)**2 + cVector(3)**2)
          polyArea(numPoly) = polyArea(numPoly) + cLength / 2
          !
          ! get rotation angles, save sines and cosines
          !
          if (cLenInXY < 1.e-10) then
            sinGam = 0.
            cosGam = 1.
          else
            sinGam = -cVector(2) / cLenInXY
            cosGam = cVector(1) / cLenInXY
          endif
          sinBet = cLenInXY / cLength
          cosBet = cVector(3) / cLength
          cosGamma(numTriang) = cosGam
          sinGamma(numTriang) = sinGam
          cosBeta(numTriang) = cosBet
          sinBeta(numTriang) = sinBet
          !
          ! calculate and save rotated coordinates
          !
          zRotSum = 0.
          do i = 1, 3
            iv = indVert(i, numTriang)
            triXYrot(i, 1, numTriang) = (verts(1, iv) * cosGam - &
                verts(2, iv) * sinGam) * cosBet - verts(3, iv) * sinBet
            triXYrot(i, 2, numTriang) = verts(1, iv) * sinGam + &
                verts(2, iv) * cosGam
            zRotSum = zRotSum + (verts(1, iv) * cosGam - &
                verts(2, iv) * sinGam) * sinBet + verts(3, iv) * cosBet
          enddo
          triZrot(numTriang) = zRotSum / 3.
        endif
        ilist = ilist + 3
      enddo
    endif
  enddo
  numVerts = numVerts + newVerts
  numPolyObj(numMeshLoaded) = numPoly + 1 - iobjPoly(numMeshLoaded)
  !
  ! sort the polygons into surfaces: use modind as temporary surface #
  !
  ! print *,'Sorting surfaces'
  listEnd = iobjPoly(numMeshLoaded) - 1
  do i = listEnd + 1, numPoly
    modIndTemp(i) = 0
  enddo
  iobjSurf(numMeshLoaded) = numSurf + 1
  !
  ! do until all polygons are on the list
  !
  do while(listEnd < numPoly)
    !
    ! find next unassigned polygon, assign it next surface #
    !
    ifFound = 0
    ip = iobjPoly(numMeshLoaded)
    do while(ip <= numPoly .and. ifFound == 0)
      if (modIndTemp(ip) == 0) then
        ifFound = 1
      endif
      ip = ip + 1
    enddo
    !
    ! start a new surface, set current point to the new one on list
    !
    ip = ip - 1
    if (numSurf >= LIMSURF) then
      if (ifHush == 0) print *,'Too many surfaces for arrays'
      return
    endif
    numSurf = numSurf + 1
    modIndTemp(ip) = numSurf
    listEnd = listEnd + 1
    listSurf(listEnd) = ip
    indStartSurf(numSurf) = listEnd
    numContInSurf(numSurf) = 0
    surfXmin(numSurf) = polyXmin(ip)
    surfYmin(numSurf) = polyYmin(ip)
    surfZmin(numSurf) = polyZmin(ip)
    surfXmax(numSurf) = polyXmax(ip)
    surfYmax(numSurf) = polyYmax(ip)
    surfZmax(numSurf) = polyZmax(ip)
    surfArea(numSurf) = polyArea(ip)
    ! print *,'starting surface ', nsurf, ' with polygon', ip
    listCur = listEnd
    !
    ! search for more to add to list until current is caught up to end
    !
    do while(listCur <= listEnd)
      ipCur = listSurf(listCur)
      do ip = iobjPoly(numMeshLoaded), numPoly
        if (modIndTemp(ip) == 0) then
          !
          ! check for bounding box overlap, don't worry about whether
          ! on same plane in Z.  THIS WILL FAIL FOR MESHES WITH FORCED
          ! CONNECTIONS BETWEEN CONTOURS THAT DON'T OVERLAP
          !
          if (polyXmin(ip) <= polyXmax(ipCur) .and. &
              polyXmin(ipCur) <= polyXmax(ip) .and. &
              polyYmin(ip) <= polyYmax(ipCur) .and. &
              polyYmin(ipCur) <= polyYmax(ip) .and. &
              polyZmin(ip) <= polyZmax(ipCur) .and. &
              polyZmin(ipCur) <= polyZmax(ip)) then
            !
            ! confirm connection by searching for common vertex
            !
            ifFound = 0
            itr = indStartPoly(ipCur)
            do while(itr < indStartPoly(ipCur) + numInPoly(ipCur) .and. &
                ifFound == 0)
              ind1 = indVert(1, itr)
              ind2 = indVert(2, itr)
              ind3 = indVert(3, itr)
              jtr = indStartPoly(ip)
              do while(jtr < indStartPoly(ip) + numInPoly(ip) .and. ifFound == 0)
                do iv = 1, 3
                  itest = indVert(iv, jtr)
                  if (ind1 == itest .or. ind2 == itest .or. ind3 == itest) ifFound = 1
                enddo
                jtr = jtr + 1
              enddo
              itr = itr + 1
            enddo
            !
            ! if there is a common vertex, then add the polygon being
            ! tested to the surface and the list.
            ! THIS COULD BE STRENGTHENED BY REQUIRING MULTIPLE COMMON
            ! POINTS, BUT THAT MIGHT NOT ALWAYS WORK EITHER
            !
            if (ifFound > 0) then
              modIndTemp(ip) = numSurf
              listEnd = listEnd + 1
              listSurf(listEnd) = ip
              surfXmin(numSurf) = min(surfXmin(numSurf), polyXmin(ip))
              surfYmin(numSurf) = min(surfYmin(numSurf), polyYmin(ip))
              surfZmin(numSurf) = min(surfZmin(numSurf), polyZmin(ip))
              surfXmax(numSurf) = max(surfXmax(numSurf), polyXmax(ip))
              surfYmax(numSurf) = max(surfYmax(numSurf), polyYmax(ip))
              surfZmax(numSurf) = max(surfZmax(numSurf), polyZmax(ip))
              surfArea(numSurf) = surfArea(numSurf) + polyArea(ip)
            endif
          endif
        endif
      enddo
      !
      ! end of scan through polygons; advance the current list position
      !
      listCur = listCur + 1
      !
      ! if there are gaps to connect across and the surface appears done
      ! then see if this surface has an edge in the region.
      !
      if (numGaps > 0 .and. listCur > listEnd) then
        zIncr = 1.01 * zTotScale
        do indGap = 1, numGaps
          if (surfZmax(numSurf) + zIncr >= zGapStart(indGap) .and. &
              surfZmin(numSurf) - zIncr <= zGapEnd(indGap)) then
            do ilist = indStartSurf(numSurf), listEnd
              ipoly = listSurf(ilist)
              ifMostIn = 0
              if (polyZmax(ipoly) + zIncr >= zGapStart(indGap) .and. &
                  polyZmin(ipoly) - zIncr <= zGapEnd(indGap)) then
                !
                ! for a polygon on the gap, first search for the closest
                ! polygon in free list whose bounding box overlaps and
                ! that is in gap
                !
                delZmin = 1.e10
                numPolyMin = 0
                do ip = iobjPoly(numMeshLoaded), numPoly
                  if (modIndTemp(ip) == 0) then
                    if (polyXmin(ip) <= polyXmax(ipoly) .and. &
                        polyXmin(ipoly) <= polyXmax(ip) .and. &
                        polyYmin(ip) <= polyYmax(ipoly) .and. &
                        polyYmin(ipoly) <= polyYmax(ip) .and. &
                        polyZmin(ip) - zIncr <= zGapEnd(indGap) .and. &
                        polyZmax(ip) + zIncr >= zGapStart(indGap)) then
                      delZ = abs(polyZmin(ip) - polyZmin(ipoly))
                      if (delZ < delZmin) then
                        delZmin = delZ
                        numPolyMin = 1
                        ipolyMin(1) = ip
                      elseif (delZ == delZmin .and. numPolyMin < LIMIPMIN) then
                        numPolyMin = numPolyMin + 1
                        ipolyMin(numPolyMin) = ip
                      endif
                    endif
                  endif
                enddo
                !
                ! now that direction is known, make sure that the
                ! polygon in the surface is actually the farthest one in
                !
                if (numPolyMin > 0) then
                  ifMostIn = 1
                  do jlist = indStartSurf(numSurf), listEnd
                    jpoly = listSurf(jlist)
                    if (jpoly .ne. ipoly .and. &
                        ((polyZmax(ipolyMin(1)) > polyZmax(ipoly) .and. &
                        polyZmax(jpoly) > polyZmax(ipoly) .and. &
                        polyZmax(jpoly) + zIncr <= zGapEnd(indGap)) .or. &
                        (polyZmax(ipolyMin(1)) < polyZmax(ipoly) .and. &
                        polyZmin(jpoly) < polyZmin(ipoly) .and. &
                        polyZmin(jpoly) - zIncr >= zGapStart(indGap))) .and. &
                        polyXmin(ipoly) <= polyXmax(jpoly) .and. &
                        polyXmin(jpoly) <= polyXmax(ipoly) .and. &
                        polyYmin(ipoly) <= polyYmax(jpoly) .and. &
                        polyYmin(jpoly) <= polyYmax(ipoly)) ifMostIn = 0
                  enddo
                  !
                  ! if this is the farther one into the gap, extract the
                  ! contours
                  !
                  if (ifMostIn .ne. 0) then
                    do ipt = 1, numPolyMin
                      ip = ipolyMin(ipt)
                      zExtract1 = polyZmax(ipoly)
                      zExtract2 = polyZmin(ip)
                      if (polyZmax(ip) < polyZmax(ipoly)) then
                        zExtract1 = polyZmin(ipoly)
                        zExtract2 = polyZmax(ip)
                      endif
                      call contour_from_poly(ipoly, zExtract1, cont1X, cont1Y, &
                          numInCont1, LIMCXY)
                      call contour_from_poly(ip, zExtract2, cont2X, cont2Y, &
                          numInCont2, LIMCXY)
                      !
                      ! check each contour for points inside the other
                      !
                      anyinside = .false.
                      if (numInCont2 > 2 .and. numInCont1 > 2) then
                        ic1 = 1
                        do while(ic1 <= numInCont1 .and. .not.anyinside)
                          anyinside = inside(cont2X, cont2Y, numInCont2, cont1X(ic1), &
                              cont1Y(ic1))
                          ic1 = ic1 + 1
                        enddo
                        ic2 = 1
                        do while(ic2 <= numInCont2 .and. .not.anyinside)
                          anyinside = inside(cont1X, cont1Y, numInCont1, cont2X(ic2), &
                              cont2Y(ic2))
                          ic2 = ic2 + 1
                        enddo
                      endif
                      !
                      ! add to surface
                      !
                      if (anyinside) then
                        modIndTemp(ip) = numSurf
                        listEnd = listEnd + 1
                        listSurf(listEnd) = ip
                        surfXmin(numSurf) = min(surfXmin(numSurf), &
                            polyXmin(ip))
                        surfYmin(numSurf) = min(surfYmin(numSurf), &
                            polyYmin(ip))
                        surfZmin(numSurf) = min(surfZmin(numSurf), &
                            polyZmin(ip))
                        surfXmax(numSurf) = max(surfXmax(numSurf), &
                            polyXmax(ip))
                        surfYmax(numSurf) = max(surfYmax(numSurf), &
                            polyYmax(ip))
                        surfZmax(numSurf) = max(surfZmax(numSurf), &
                            polyZmax(ip))
                        surfArea(numSurf) = surfArea(numSurf) + &
                            polyArea(ip)
                      endif
                    enddo
                  endif
                endif
              endif
            enddo
          endif
        enddo
      endif
    enddo
    !
    ! surface all found, finish assignments for it
    !
    numInSurf(numSurf) = listEnd + 1 - indStartSurf(numSurf)
  enddo
  numSurfObj(numMeshLoaded) = numSurf + 1 - iobjSurf(numMeshLoaded)
  !
  ! GIVEN SURFACES, NEED TO APPLY SHIFTS NOW IF ANY EXIST
  !
  ifZeroedShift = 0
  do iobj = 1, numObjShifted
    if (iobjShift(iobj) == imodObj) then
      if (numItemShifted(iobj) == numSurfObj(numMeshLoaded)) then
        do ii = 1, numItemShifted(iobj)
          ishift = ii + indStartShift(iobj) - 1
          isurf = iobjSurf(numMeshLoaded) + ii - 1
          call shiftSurfMesh(isurf, shifts(1, ishift), 1, modIndTemp, &
              ifZeroedShift, zTotScale)
        enddo
      else
        if (ifHush == 0) print *, &
            'Number of shifts and surfaces do not match'
      endif
    endif
  enddo
  !
  ! NOW CROSS-REFERENCE CONTOURS TO SURFACES
  !
  ! print *,'Cross-referencing contours'
  numTemp = 0
  numOrphan = 0
  toler = 0.1 * xyScale
  do iobj = 1, max_mod_obj
    if (obj_color(2, iobj) == 256 - imodObj .and. npt_in_obj(iobj) > 2) &
        then
      !
      ! look for contours in object with at least 3 points
      !
      numInObj = npt_in_obj(iobj)
      ibase = ibase_obj(iobj)
      txmin = 1.e10
      txmax = -txmin
      tymin = txmin
      tymax = txmax
      !
      ! get min and max, copy contour
      !
      do ipt = 1, min(numInObj, LIMCXY)
        ipoint = abs(object(ibase + ipt))
        cont1X(ipt) = p_coord(1, ipoint)
        cont1Y(ipt) = p_coord(2, ipoint)
        txmin = min(txmin, cont1X(ipt))
        tymin = min(tymin, cont1Y(ipt))
        txmax = max(txmax, cont1X(ipt))
        tymax = max(tymax, cont1Y(ipt))
      enddo
      tzval = p_coord(3, ipoint)
      tzmax = tzval + toler
      tzmin = tzval - toler
      !
      ! scan surfaces of this object
      !
      isurf = iobjSurf(numMeshLoaded)
      icontSurf = 0
      do while(isurf <= iobjSurf(numMeshLoaded) + numSurfObj(numMeshLoaded) - 1 &
          .and. icontSurf == 0)
        if (txmin <= surfXmax(isurf) .and. surfXmin(isurf) <= txmax .and. &
            tymin <= surfYmax(isurf) .and. surfYmin(isurf) <= tymax &
            .and. tzmax >= surfZmin(isurf) .and. tzmin <= surfZmax(isurf)) &
            then
          !
          ! if intersect, check through polygons
          !
          list = indStartSurf(isurf)
          do while(list <= indStartSurf(isurf) + numInSurf(isurf) - 1 &
              .and. icontSurf == 0)
            ipoly = listSurf(list)
            if (txmin <= polyXmax(ipoly) .and. polyXmin(ipoly) <= txmax .and. &
                tymin <= polyYmax(ipoly) .and. polyYmin(ipoly) <= tymax &
                .and. tzmax >= polyZmin(ipoly) .and. &
                tzmin <= polyZmax(ipoly)) then
              !
              ! then check triangles
              !
              itriang = indStartPoly(ipoly)
              do while(itriang <= indStartPoly(ipoly) + numInPoly(ipoly) - 1 &
                  .and. icontSurf == 0)
                !
                ! check each vertex of triangle against each of contour
                !
                do ivert = 1, 3
                  iv = indVert(ivert, itriang)
                  do ipt = 1, min(numInObj, LIMCXY)
                    if (abs(verts(1, iv) - cont1X(ipt)) < toler .and. &
                        abs(verts(2, iv) - cont1Y(ipt)) < toler .and. &
                        abs(verts(3, iv) - tzval) < toler) icontSurf = isurf
                  enddo
                enddo
                itriang = itriang + 1
              enddo
            endif
            list = list + 1
          enddo
        endif
        isurf = isurf + 1
      enddo
      !
      ! got out, found or not
      !
      if (icontSurf .ne. 0) then
        if (numConts + numTemp >= LIMCONT) then
          if (ifHush == 0) print *,'Too many contours for arrays'
          return
        endif
        numTemp = numTemp + 1
        modIndTemp(numTemp) = iobj
        modIndTemp(numTemp + max_mod_obj) = icontSurf
        numContInSurf(icontSurf) = numContInSurf(icontSurf) + 1
        tempXmin(numTemp) = txmin
        tempYmin(numTemp) = tymin
        tempXmax(numTemp) = txmax
        tempYmax(numTemp) = tymax
        tempZval(numTemp) = tzval
      else
        numOrphan = numOrphan + 1
        ! print *,iobj, ' is orphan'
      endif
    endif
  enddo
  !
  ! put the contours in lists, and everything else
  !
  indStartCont(iobjSurf(numMeshLoaded)) = numConts + 1
  do isurf = iobjSurf(numMeshLoaded), &
      iobjSurf(numMeshLoaded) + numSurfObj(numMeshLoaded) - 1
    indStartCont(isurf + 1) = indStartCont(isurf) + numContInSurf(isurf)
    numContInSurf(isurf) = 0
    numConts = indStartCont(isurf + 1) - 1
  enddo
  do i = 1, numTemp
    isurf = modIndTemp(i + max_mod_obj)
    ind = indStartCont(isurf) + numContInSurf(isurf)
    numContInSurf(isurf) = numContInSurf(isurf) + 1
    listCont(ind) = modIndTemp(i)
    contXmin(ind) = tempXmin(i)
    contYmin(ind) = tempYmin(i)
    contXmax(ind) = tempXmax(i)
    contYmax(ind) = tempYmax(i)
    contZval(ind) = tempZval(i)
  enddo
  numAddCont = numConts + 1 - indStartCont(iobjSurf(numMeshLoaded))

  ! check integrity of all loaded mesh and contour data
  !
  do indMesh = 1, numMeshLoaded
    imodObj = iobjMesh(indMesh)
    itypeObj = 256 - imodObj
    do isurf = iobjSurf(indMesh), iobjSurf(indMesh) + numSurfObj(indMesh) - 1
      do indCont = indStartCont(isurf), indStartCont(isurf) + numContInSurf(isurf) - 1
        iobj = listCont(indCont)
        if (itypeObj .ne. obj_color(2, iobj)) then
          print *,'object mismatch', indMesh, imodObj, isurf, indCont, &
              iobj, itypeObj, obj_color(2, iobj)
        endif
      enddo
    enddo
  enddo

  ifErr = 0
  if (ifHush .ne. 0) return
  write(*,109) imodObj, newVerts, newInd, &
      numTriang + 1 - indStartPoly(iobjPoly(numMeshLoaded)), &
      numPolyObj(numMeshLoaded), numSurfObj(numMeshLoaded), numAddCont
109 format('Object',i4,':',i9,' vertices,',i9,' indices,' &
      ,i9,' triangles,',/,i7,' polygons,',i7,' surfaces,',i7, &
      ' contours')
  if (numDupTriang > 0) print *,numDupTriang, &
      ' duplicate triangles ignored'
  if (numNonTriang > 0) print *,numNonTriang, &
      ' non-triangles ignored'
  if (numOrphan .ne. 0) print *,numOrphan, &
      ' orphan contours not referenced to mesh'
  return
end subroutine process_mesh



! CONTOUR_FROM_POLY extracts a contour at the z level given by ZEXTR
! from polygon IP; NC points are placed in CX, CY; LIMC specifies the
! limiting size of the CX and CY arrays

subroutine contour_from_poly(ip, zExtract, contX, contY, numC, limc)
  use mtkvars
  implicit none
  real*4 contX(*), contY(*), zExtract
  integer*4 ip,  numC, limc
  integer*4 last, itr, iv, ind
  last = -1
  numC = 0
  do itr = indStartPoly(ip), indStartPoly(ip) + numInPoly(ip) - 1
    do iv = 1, 3
      ind = indVert(iv, itr)
      if (ind .ne. last .and. verts(3, ind) == zExtract .and. numC < limc) then
        numC = numC + 1
        contX(numC) = verts(1, ind)
        contY(numC) = verts(2, ind)
        last = ind
      endif
    enddo
  enddo
  return
end subroutine contour_from_poly

subroutine errorexit(message)
  character*(*) message
  print *
  print *,'ERROR: MTK - ', message
  call exit(1)
end subroutine errorexit

