! RESAMPLEMOD
!
! A program to turn a fiber model so that it is perpendicular to the
! Z axis, then resample at section intervals so that it can be
! analyzed with NDA.
! See man page for details
!
! $Id$
!
program resamplemod
  implicit none
  include 'model.inc90'
  integer LIMSEC, LIMOBJ
  parameter (LIMSEC = 100000, LIMOBJ = 10000)
  character*320 modelFile, newModel
  logical exist, readw_or_imod
  real*4 rmat(3,3), vecMean(3), vecNormal(3)
  integer*4 iobjUse(LIMOBJ), iobjExclude(LIMOBJ)
  integer getImodHead
  real*4 xmt(LIMSEC), ymt(LIMSEC), zmt(LIMSEC)
  logical failed
  integer*4 ierr, numObjUse, numObjExclude, ifFlipped, indZ, indY, mainAxis
  real*4 xyScale, zScale, xOffset, yOffset, zOffset, xMin, yMin, zMin, xTemp, yTemp, zTemp
  integer*4 iobj, numInObj, ibase, ipt, itmp, numVectors, i, ifUse, ifExclude, j
  real*4 xShift, yShift, zShift, temp, frac, xMax, yMax, zMax
  integer*4 izStart, izEnd, newNumInObj, iz, maxX, maxY, maxZ, imodObj, imodCont
  integer*4 skipInvert
  integer*4 numberInList

  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetInteger, PipGetBoolean
  integer*4 PipGetString
  integer*4 PipGetInOutFile
  !
  ! fallbacks from ../../manpages/autodoc2man -2 2  resamplemod
  !
  integer numOptions
  parameter (numOptions = 8)
  character*(40 * numOptions) options(1)
  options(1) = &
      'input:InputFile:FN:@output:OutputFile:FN:@'// &
      'exclude:ExcludeObjects:LI:@'// &
      'direction:DirectionObjects:LI:@main:MainAxis:I:@'// &
      'skip:SkipInversion:I:@param:ParameterFile:PF:@'// &
      'help:usage:B:'

  skipInvert = 0
  mainAxis = 3
  !
  ! Pip startup: set error, parse options, check help, set flag if used
  !
  call PipReadOrParseOptions(options, numOptions, 'resamplemod', &
      'ERROR: RESAMPLEMOD - ', .false., 2, 1, 1, numOptArg, numNonOptArg)
  !
  if (PipGetInOutFile('InputFile', 1, ' ', modelFile) .ne. 0) &
      call exitError('No input file specified')
  !
  ! read in the model
  !
  exist = readw_or_imod(modelFile)
  if (.not.exist) call exitError('Reading model')

  if (PipGetInOutFile('OutputFile', 2, ' ', newModel) .ne. 0) &
      call exitError('No output file specified')

  numObjUse = 0
  if (PipGetString('DirectionObjects', modelFile) == 0) &
      call parselist(modelFile, iobjUse, numObjUse)

  numObjExclude = 0
  if (PipGetString('ExcludeObjects', modelFile) == 0) &
      call parselist(modelFile, iobjExclude, numObjExclude)

  ifFlipped = 0
  ierr = getImodHead(xyScale, zScale, xOffset, yOffset, zOffset, ifFlipped)
  if (ierr .ne. 0) call exitError('Getting scaling from model header')
  call scale_model(0)

  if (ifFlipped == 0) then
    indY = 2
    indZ = 3
  else
    indZ = 2
    indY = 3
  endif

  ierr = PipGetInteger('SkipInversion', skipInvert)
  ierr = PipGetInteger('MainAxis', mainAxis)
  mainAxis = max(1, min(3, mainAxis))
  if (mainAxis == 2) then
    mainAxis = indY
  elseif (mainAxis == 3) then
    mainAxis = indZ
  endif
  call PipDone()
  !
  ! invert objects if starting z > ending z.
  !
  if (skipInvert == 0) then

    do iobj = 1, max_mod_obj
      numInObj = npt_in_obj(iobj)
      if (numInObj > 0) then
        ibase = ibase_obj(iobj)
        if (p_coord(mainAxis, abs(object(1 + ibase))) > &
            p_coord(mainAxis, abs(object(numInObj + ibase)))) then
          do ipt = 1, numInObj / 2
            itmp = object(ibase + ipt)
            object(ibase + ipt) = object(ibase + numInObj + 1 - ipt)
            object(ibase + numInObj + 1 - ipt) = itmp
          enddo
        endif
      endif
    enddo
  endif

  do i = 1, 3
    vecMean(i) = 0.
    vecNormal(i) = 0.
  enddo
  vecNormal(indZ) = 1.
  numVectors = 0
  !
  ! get mean vector to turn to vertical, including Z scale
  !
  do iobj = 1, max_mod_obj
    call objToCont(iobj, obj_color, imodObj, imodCont)
    ifUse = numberInList(imodObj, iobjUse, numObjUse, 1)
    numInObj = npt_in_obj(iobj)
    if (numInObj > 1 .and. ifUse == 1) then
      ibase = ibase_obj(iobj)
      i = abs(object(ibase + 1))
      ipt = abs(object(ibase + numInObj))
      !
      ! If inversion was skipped, check if this object needs to be
      ! inverted for adding into the sum
      !
      if (skipInvert > 0 .and. &
          p_coord(mainAxis, i) > p_coord(mainAxis, ipt)) then
        i = abs(object(ibase + numInObj))
        ipt = abs(object(ibase + 1))
      endif

      vecMean(1) = vecMean(1) + p_coord(1, ipt) - p_coord(1, i)
      vecMean(indY) = vecMean(indY) + p_coord(indY, ipt) - p_coord(indY, i)
      vecMean(indZ) = vecMean(indZ) + zScale * ( &
          p_coord(indZ, ipt) - p_coord(indZ, i))
      numVectors = numVectors + 1
    endif
  enddo
  if (numVectors == 0) call exitError('No contours found')
  do i = 1, 3
    vecMean(i) = vecMean(i) / numVectors
  enddo

  call vectorsToRotMat(vecMean, vecNormal, rmat)

  write(*,'(a,3f10.1)') 'Mean vector:', (vecMean(i), i = 1, 3)
  print *,'Rotation matrix:'
  write(*,'(3f8.4)') ((rmat(i, ipt), i = 1, 3), ipt = 1, 3)
  xTemp = vecMean(1)
  yTemp = vecMean(2)
  zTemp = vecMean(3)
  xmt(1) = xTemp * rmat(1, 1) + yTemp * rmat(1, 2) + zTemp * rmat(1, 3)
  ymt(1) = xTemp * rmat(2, 1) + yTemp * rmat(2, 2) + zTemp * rmat(2, 3)
  zmt(1) = (xTemp * rmat(3, 1) + yTemp * rmat(3, 2) + zTemp * rmat(3, 3)) / zScale
  write(*,'(a,3f10.1)') 'Rotated vector:', xmt(1), ymt(1), zmt(1)

  xMin = 1.e10
  yMin = 1.e10
  zMin = 1.e10
  do iobj = 1, max_mod_obj
    call objToCont(iobj, obj_color, imodObj, imodCont)
    ifExclude = numberInList(imodObj, iobjExclude, numObjExclude, 0)
    numInObj = npt_in_obj(iobj)
    if (numInObj > 0) then
      ibase = ibase_obj(iobj)
      !
      ! rotate the points in 3-D, scaling then unscaling by Z
      !
      do i = 1, numInObj
        ipt = abs(object(i + ibase))
        xTemp = p_coord(1, ipt)
        yTemp = p_coord(indY, ipt)
        zTemp = zScale * p_coord(indZ, ipt)
        xmt(i) = xTemp * rmat(1, 1) + yTemp * rmat(1, 2) + zTemp * rmat(1, 3)
        ymt(i) = xTemp * rmat(2, 1) + yTemp * rmat(2, 2) + zTemp * rmat(2, 3)
        zmt(i) = (xTemp * rmat(3, 1) + yTemp * rmat(3, 2) + zTemp * rmat(3, 3)) / zScale
        xMin = min(xmt(i), xMin)
        yMin = min(ymt(i), yMin)
        zMin = min(zmt(i), zMin)
      enddo
      xShift = 10 - xMin
      yShift = 10 - yMin
      zShift = 1 - nint(zMin)
      !
      ! if not resampling, just place points back into contour
      !
      if (ifExclude == 1 .or. numInObj == 1) then
        do i = 1, numInObj
          ipt = abs(object(i + ibase))
          p_coord(1, ipt) = xmt(i)
          p_coord(indY, ipt) = ymt(i)
          p_coord(indZ, ipt) = zmt(i)
        enddo
      else
        !
        ! invert contour if this has gotten it backwards
        !
        if (zmt(1) > zmt(numInObj)) then
          do i = 1, numInObj / 2
            j = numInObj + 1 - i
            temp = xmt(i)
            xmt(i) = xmt(j)
            xmt(j) = temp
            temp = ymt(i)
            ymt(i) = ymt(j)
            ymt(j) = temp
            temp = zmt(i)
            zmt(i) = zmt(j)
            zmt(j) = temp
          enddo
        endif
        !
        ! get z limits for resampling, set first point into contour
        !
        call object_mover(iobj, failed)
        if (failed) call exitError('Insufficient object space')
        ibase = ibase_obj(iobj)
        izStart = nint(zmt(1))
        if (izStart <= zmt(1) + 0.01) izStart = izStart + 1
        izEnd = nint(zmt(numInObj))
        if (izEnd >= zmt(numInObj) - 0.01) izEnd = izEnd-1
        ipt = abs(object(1 + ibase))
        p_coord(1, ipt) = xmt(1)
        p_coord(indY, ipt) = ymt(1)
        p_coord(indZ, ipt) = zmt(1)
        newNumInObj = 1
        !
        ! interpolate a point at each new Z level
        !
        do iz = izStart, izEnd
          j = 1
          do while (j < numInObj .and. &
              .not.(zmt(j) <= iz .and. zmt(j + 1) > iz))
            j = j + 1
          enddo
          newNumInObj = newNumInObj + 1
          if (newNumInObj > numInObj) then
            n_point = n_point + 1
            ntot_in_obj = ntot_in_obj + 1
            object(ibase + newNumInObj) = n_point
            pt_label(n_point) = 0
          endif
          ipt = object(ibase + newNumInObj)
          frac = 0.
          if (zmt(j + 1) - zmt(j) > 0.01) &
              frac = (iz - zmt(j)) / (zmt(j + 1) - zmt(j))
          p_coord(1, ipt) = xmt(j) + frac * (xmt(j + 1) - xmt(j))
          p_coord(indY, ipt) = ymt(j) + frac * (ymt(j + 1) - ymt(j))
          p_coord(indZ, ipt) = iz
        enddo
        !
        ! finish with last point
        !
        newNumInObj = newNumInObj + 1
        if (newNumInObj > numInObj) then
          n_point = n_point + 1
          ntot_in_obj = ntot_in_obj + 1
          object(ibase + newNumInObj) = n_point
          pt_label(n_point) = 0
        endif
        ipt = object(ibase + newNumInObj)
        p_coord(1, ipt) = xmt(numInObj)
        p_coord(indY, ipt) = ymt(numInObj)
        p_coord(indZ, ipt) = zmt(numInObj)
        npt_in_obj(iobj) = newNumInObj
        ibase_free = max(ibase_free, ibase + newNumInObj)
      endif
    endif
  enddo
  !
  ! shift the points to be positive
  !
  xMax = -1000.
  yMax = -1000.
  zMax = -1000.
  do iobj = 1, max_mod_obj
    numInObj = npt_in_obj(iobj)
    if (numInObj > 0) then
      ibase = ibase_obj(iobj)
      do ipt = 1, numInObj
        i = abs(object(ipt + ibase))
        p_coord(1, i) = p_coord(1, i) + xShift
        p_coord(indY, i) = p_coord(indY, i) + yShift
        p_coord(indZ, i) = p_coord(indZ, i) + zShift
        xMax = max(xMax, p_coord(1, i))
        yMax = max(yMax, p_coord(indY, i))
        zMax = max(zMax, p_coord(indZ, i))
      enddo
    endif
  enddo
  maxX = 10 * (int(xMax + 19) / 10)
  maxY = 10 * (int(yMax + 19) / 10)
  maxZ = nint(zMax) + 2
  call scale_model(1)
  call putImodMaxes(maxX, maxY, maxZ)
  call write_wmod(newModel)
  call exit
end program resamplemod


subroutine vectorsToRotMat(avec, bvec, rotm)
  implicit none
  real*4 rotm(3,3), avec(*), bvec(*), cvec(3)
  real*4 x, y, z, sinAngle, cosAngle, angle, sind, acosd, oneMinusCosA, vecMagnitude
  equivalence (x, cvec(1)), (y, cvec(2)), (z, cvec(3))
  integer*4 i, j

  call crossProduct(avec, bvec, cvec)
  do i = 1, 3
    do j = 1, 3
      rotm(i, j) = 0.
      if (i == j) rotm(i, j) = 1.
    enddo
  enddo
  if (vecMagnitude(cvec) == 0.) return
  cosAngle = (avec(1) * bvec(1) + avec(2) * bvec(2) + avec(3) * bvec(3)) / &
      (vecMagnitude(avec) * vecMagnitude(bvec))
  angle = acosd(cosAngle)
  sinAngle = -sind(angle)
  oneMinusCosA = 1.0 - cosAngle
  call normalize(cvec)
  rotm(1, 1) = x * x * oneMinusCosA + cosAngle
  rotm(1, 2) = y * x * oneMinusCosA + (sinAngle * z)
  rotm(1, 3) = z * x * oneMinusCosA - (sinAngle * y)

  rotm(2, 1) = x * y * oneMinusCosA - (sinAngle * z)
  rotm(2, 2) = y * y * oneMinusCosA + cosAngle
  rotm(2, 3) = z * y * oneMinusCosA + (sinAngle * x)

  rotm(3, 1) = x * z * oneMinusCosA + (sinAngle * y)
  rotm(3, 2) = y * z * oneMinusCosA - (sinAngle * x)
  rotm(3, 3) = z * z * oneMinusCosA + cosAngle
  return
end subroutine vectorsToRotMat



real*4 function vecMagnitude(vec)
  implicit none
  real*4 vec(3)
  vecMagnitude = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
  return
end function vecMagnitude


subroutine crossProduct(avec, bvec, cvec)
  implicit none
  real*4 avec(*), bvec(*), cvec(*)
  cvec(1) = avec(2) * bvec(3) - avec(3) * bvec(2)
  cvec(2) = avec(3) * bvec(1) - avec(1) * bvec(3)
  cvec(3) = avec(1) * bvec(2) - avec(2) * bvec(1)
  return
end subroutine crossProduct



subroutine normalize(vect)
  implicit none
  real*4 vect(3), sum, fac
  integer*4 i
  sum = 0.
  do i = 1, 3
    sum = sum + vect(i)**2
  enddo
  fac = sqrt(sum)
  do i = 1, 3
    vect(i) = vect(i) / fac
  enddo
  return
end subroutine normalize
