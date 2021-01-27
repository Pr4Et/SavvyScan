! * * * * * SOLVEMATCH * * * * * *
!
! SOLVEMATCH will solve for a 3-dimensional linear transformation
! relating the tomogram volumes resulting from tilt series around two
! different axes.  It uses information from the 3-D coordinates of
! fiducials found by TILTALIGN.  See man page for further information.
!
! $Id$
!
program solvematch
  implicit none
  include 'model.inc90'
  integer IDIM, MAT_SIZE
  parameter (IDIM = 10000, MAT_SIZE = 20)
  real*4 xMat(MAT_SIZE,IDIM)
  real*4 pointsA(3,IDIM), pointsB(3,IDIM), a(3,4), delXYZ(3)
  real*4 devXyzMax(3), pointRot(3), cenMeanLoc(3), amatLocal(3,4), dxyzLocal(3)
  integer*4 idrop(IDIM), indOrig(IDIM), mapped(IDIM)
  integer*4 listCorrA(IDIM), listCorrB(IDIM), icontToPointB(IDIM)
  integer*4 mapAtoB(IDIM), nxyz(3,2), jxyz(3) /1, 3, 2/
  integer*4 icontA(IDIM), icontB(IDIM), icontToPointA(IDIM)
  integer*4 icontAB(IDIM,2), icontToPointAB(IDIM,2), listCorrAB(IDIM,2)
  equivalence (icontA, icontAB), (icontB, icontAB(1, 2))
  equivalence (icontToPointA, icontToPointAB), (icontToPointB, icontToPointAB(1, 2))
  equivalence (listCorrA, listCorrAB), (listCorrB, listCorrAB(1, 2))
  integer*4 numABpoints(2), numApoints, numBpoints, listUse(IDIM)
  equivalence (numApoints, numABpoints(1)), (numBpoints, numABpoints(2))
  real*4 transferX(IDIM,2), transferY(IDIM,2), fidModX(IDIM,2)
  real*4 fidModY(IDIM,2)
  integer*4 modObj(IDIM,2), modCont(IDIM,2), izBest(2), numFid(2)
  integer*4 modObjFid(IDIM,2), modContFid(IDIM,2)
  character*320 filename
  character*1 abText(2) /'A', 'B'/
  character*1 badAxis1, badAxis2
  integer*4 minNumToStart, numList, numListA, numListB, numData, ifZshifts
  integer*4 ia, numSurf, model, numModPts, ipt, ip, iobj, i, ndata, ifAngleOfs
  integer*4 j, ifAdded, ipntMax, iptA, iptB, iptAatMin, iptBatMin
  integer*4 maxDrop, numDrop, iofs, ixyz, ierrFactor
  real*4 addRatio, addCrit, distMin, devAvg, devSd
  real*4 devMax, dx, dist, critProb, elimMin, absProbCrit, stopLimit, xtilta, xtiltB
  real*4 dy, dz, aScale, bScale
  integer*4 ierr, ifFlip, numColFit, maxContA, maxContB, icolFixed
  real*4 xyScale, zScale, xOffset, yOffset, zOffset, xImScale, yImScale, zImScale
  real*4 loMaxAvgRatio, hiMaxAvgRatio, loMaxLimRatio, hiMaxLimRatio
  real*4 aDelta, bDelta, aPixelSize, bPixelSize, xcen, ycen, transTol
  real*4 angleOffsetA, angleOffsetB, zShiftA, zShiftB
  integer*4 nxFidA, nyFidA, nxFidB, nyFidB, iTransA, iTransB, ifTransBtoA
  integer*4 numTransCoord, izInd, indA, indB, nListUse, localNum, invertedInDepth
  integer*4 numLocalX, numLocalY, numBig, ixl, iyl, indOrigAtMax, limRaised
  real*4 xmin, xmax, ymin, ymax, size, targetSize, dxLocal, dyLocal
  real*4 sumMean, sumMax, shiftLimit, cenShiftX, cenShiftY, cenShiftZ, devAvgLoc
  real*4 devMaxLoc, devAllMax, globLocAvgRatio, axisCrit, reportDiff
  real*4 yzScaleDiff, sumSq, xyScaleDiff, xzScaleDiff, axisScale(3)
  real*4 transPixel, freeInput(20), determPositive, determInverted
  logical*4 relativeFids, matchAtoB, invertEntered

  logical readw_or_imod
  integer*4 getImodHead, getImodScales, imodGetenv
  real*4 determ3
  !
  logical pipInput
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetInteger
  integer*4 PipGetString, PipGetTwoFloats, PipGetFloat
  integer*4 PipGetInOutFile, PipGetLogical

  character*6 modelOption(2) /'AMatch', 'BMatch'/
  character*9 tomoOption(2) /'ATomogram', 'BTomogram'/
  character*10240 listString
  !
  ! fallbacks from ../../manpages/autodoc2man -3 2  solvematch
  !
  integer numOptions
  parameter (numOptions = 26)
  character*(40 * numOptions) options(1)
  options(1) = &
      'output:OutputFile:FN:@afiducials:AFiducialFile:FN:@'// &
      'bfiducials:BFiducialFile:FN:@alist:ACorrespondenceList:LI:@'// &
      'blist:BCorrespondenceList:LI:@transfer:TransferCoordinateFile:FN:@'// &
      'amodel:AFiducialModel:FN:@bmodel:BFiducialModel:FN:@use:UsePoints:LI:@'// &
      'atob:MatchingAtoB:B:@xtilts:XAxisTilts:FP:@angles:AngleOffsetsToTilt:FP:@'// &
      'zshifts:ZShiftsToTilt:FP:@surfaces:SurfacesOrUseModels:I:@'// &
      'inverted:InvertedInDepth:I:@maxresid:MaximumResidual:F:@local:LocalFitting:I:@'// &
      'center:CenterShiftLimit:F:@amatch:AMatchingModel:FN:@'// &
      'bmatch:BMatchingModel:FN:@atomogram:ATomogramOrSizeXYZ:CH:@'// &
      'btomogram:BTomogramOrSizeXYZ:CH:@scales:ScaleFactors:FP:@'// &
      'aniso:AnisotropicLimit:F:@param:ParameterFile:PF:@help:usage:B:'
  !
  ! don't add a point if it's this much higher than the limit for
  ! quitting
  !
  addRatio = 1.
  !
  ! require at least this many points to start
  !
  minNumToStart = 4
  !
  ! initialize these items here in case no fiducials are entered
  !
  numApoints = 0
  numBpoints = 0
  numListA = 0
  numListB = 0
  numData = 0
  numSurf = 0
  numColFit = 3
  icolFixed = 0
  filename = ' '
  xtilta = 0.
  xtiltB = 0.
  angleOffsetA = 0.
  angleOffsetB = 0.
  localNum = 0
  shiftLimit = 10.
  zShiftA = 0.
  zShiftB = 0.
  stopLimit = 8.
  aScale = 1.
  bScale = 1.
  aDelta = 0.
  bDelta = 0.
  aPixelSize = 0.
  bPixelSize = 0.
  relativeFids = .true.
  indA = 1
  indB = 2
  transTol = 3.
  numTransCoord = 0
  matchAtoB = .false.
  axisCrit = 10.
  transPixel = 0.
  invertEntered = .false.
  !
  ! Pip startup: set error, parse options, check help, set flag if used
  !
  call PipReadOrParseOptions(options, numOptions, 'solvematch', &
      'ERROR: SOLVEMATCH - ', .true., 3, 0, 1, numOptArg, numNonOptArg)
  pipInput = numOptArg + numNonOptArg > 0
  !
  ! Get first fid filename if any
  ! Also get the surfaces option; if it's -2, then force models only
  ! even if there are fids
  !
  if (pipInput) then
    ierr = PipGetLogical('MatchingAtoB', matchAtoB)
    if (matchAtoB) indA = 2
    indB = 3 - indA
    ierr = PipGetString('AFiducialFile', filename)
    numSurf = 2
    ierr = PipGetInteger('SurfacesOrUseModels', numSurf)
    if (numSurf == -2) then
      numSurf = 0
      filename = ' '
    endif
    invertEntered = PipGetInteger('InvertedInDepth', invertedInDepth) == 0
  else
    write(*,'(1x,a,/,a$)') 'Name of file with 3-D fiducial'// &
        ' coordinates for first tilt series,', &
        ' or Return to align with matching model files only: '
    read(*,'(a)') filename
  endif
  !
  ! If no fiducials are to be used, skip down to maximum residual entry
  ! after initializing some variables
  maxContA = 0
  maxContB = 0
  if (filename == ' ') go to 40
  !
  numColFit = 4
  call getFiducials(filename, icontA, pointsA, numApoints, modObj(1, 1), &
      modCont(1, 1), IDIM, 'first tomogram', aPixelSize,  nxFidA, nyFidA)
  !
  if (pipInput) then
    if (PipGetString('BFiducialFile', filename) > 0) call exitError &
        ('No file specified for fiducials in second tilt series')
  else
    write(*,'(1x,a,$)') 'Name of file with 3-D fiducial'// &
        ' coordinates for second tilt series: '
    read(*,'(a)') filename
  endif
  call getFiducials(filename, icontB, pointsB, numBpoints, modObj(1, 2), &
      modCont(1, 2), IDIM, 'second tomogram', bPixelSize, nxFidB, nyFidB)
  !
  ! fill listcorr with actual contour numbers, and find maximum contours
  !
  do i = 1, numApoints
    listCorrA(i) = icontA(i)
    mapAtoB(i) = 0
    maxContA = max(maxContA, icontA(i))
  enddo
  do i = 1, numBpoints
    listCorrB(i) = icontB(i)
    mapped(i) = 0
    maxContB = max(maxContB, icontB(i))
  enddo
  !
  ! build index from contour numbers to points in array
  !
  if (maxContA > IDIM .or. maxContB > IDIM) call exitError( &
      'Contour numbers too high for arrays')
  do j = 1, 2
    do i = 1, max(maxContA, maxContB)
      icontToPointAB(i, j) = 0
    enddo
    do i = 1, numABpoints(j)
      icontToPointAB(icontAB(i, j), j) = i
    enddo
  enddo
  !
  if (pipInput .and. &
      PipGetString('TransferCoordinateFile', filename) == 0) then
    !
    ! Read transfer file, get best section and coordinates
    !
    call dopen(1, filename, 'ro', 'f')
    read(1, '(a)') listString
    call frefor(listString, freeInput, j)
    izBest(1) = nint(freeInput(1))
    izBest(2) = nint(freeInput(2))
    ifTransBtoA = nint(freeInput(3))
    !
    ! If there is a fourth value
    ! then get the scaling to apply to model coords to match transfer coord
    if (j > 3) transPixel = freeInput(4)
    !
    i = 0
11  read(1,*,err = 12, end = 12) (transferX(i + 1, j), transferY(i + 1, j), j = 1, 2)
    i = i + 1
    if (i >= IDIM) call exitError('Too many transfer coordinates for arrays')
    goto 11
12  numTransCoord = i
    ! print *,nTransCoord, ' transfer coords'
    !
    ! Set up indexes to best z for the first and second model, and for
    ! accessing the absolute A or B axis transfer coords and fid coords
    !
    izInd = 1
    if (matchAtoB .neqv. ifTransBtoA > 0) izInd = 2
    iTransA = 1
    if (ifTransBtoA > 0) iTransA = 2
    iTransB = 3 - iTransA

    !
    ! Read the two fiducial models and get coordinates at best Z and
    ! their objects and contours
    !
    call readFidModelFile('AFiducialModel', izBest(izInd), abText(indA), &
        fidModX(1, indA), fidModY(1, indA), modObjFid(1, indA), &
        modContFid(1, indA), numFid(indA), IDIM, transPixel)
    call readFidModelFile('BFiducialModel', izBest(3 - izInd), &
        abText(indB), fidModX(1, indB), fidModY(1, indB), &
        modObjFid(1, indB), modContFid(1, indB), numFid(indB), IDIM, transPixel)
    !
    ! Get the list of points to use from the true A series
    !
    nListUse = numABpoints(indA)
    do i = 1, numABpoints(indA)
      listUse(i) = icontAB(i, indA)
    enddo
    if (PipGetString('UsePoints', listString) == 0) &
        call parselist2(listString, listUse, nListUse, IDIM)
    ! print *,nListUse, ' points:', (listUse(i), i = 1, nListUse)
    !
    ! Build list of corresponding points
    !
    numList = 0
    do i = 1, nListUse
      if (listUse(i) > 0 .and. listUse(i) <= max(maxContA, maxContB) &
          .and. icontToPointAB(listUse(i), indA) > 0) then
        !
        ! IF the point in A is a legal point from fid file, find the
        ! same object/contour in the fiducial model
        !
        iptA = icontToPointAB(listUse(i), indA)
        ! print *,'legal point', i, listUse(i), ipta
        ip = 0
        do j = 1, numFid(1)
          if (modObj(iptA, indA) == modObjFid(j, 1) .and. &
              modCont(iptA, indA) == modContFid(j, 1)) ip = j
        enddo
        if (ip > 0) then
          !
          ! Now match the model coords to the transfer coords
          !
          ! print *,'matches in a fid model:', ip
          ia = 0
          do j = 1, numTransCoord
            if (sqrt((fidModX(ip, 1) - transferX(j, iTransA))**2 + &
                (fidModY(ip, 1) - transferY(j, iTransA))**2) < transTol) ia = j
          enddo
          if (ia > 0) then
            ! print *,'matches transfer:', ia
            !
            ! Found one, then look for a match to the B transfer coords
            ! in the B fiducial model
            !
            ipt = 0
            do j = 1, numTransCoord
              if (sqrt((fidModX(j, 2) - transferX(ia, iTransB))**2 + &
                  (fidModY(j, 2) - transferY(ia, iTransB))**2) < transTol) ipt = j
            enddo
            if (ipt > 0) then
              ! print *,'matches in B fid model:', ipt
              !
              ! Find the obj/cont in the points of the fid file
              !
              iptB = 0
              do j = 1, numFid(2)
                if (modObj(j, indB) == modObjFid(ipt, 2) .and. &
                    modCont(j, indB) == modContFid(ipt, 2)) iptB = j
              enddo
              if (iptB > 0) then
                !
                ! check for uniqueness
                !
                numListB = 0
                do j = 1, numList
                  if (listCorrAB(j, indB) == icontAB(iptB, indB)) numListB = 1
                enddo
                if (numListB == 0) then
                  !
                  ! Bingo: we have corresponding points
                  !
                  ! print *,'Matching:', listUse(i), icontAB(iptb, indB)
                  numList = numList + 1
                  listCorrAB(numList, indA) = listUse(i)
                  listCorrAB(numList, indB) = icontAB(iptB, indB)
                endif
              endif
            endif
          endif
        endif
      endif
    enddo
    if (numList == 0) call exitError('No corresponding points found using transfer '// &
        'coords')
    numListA = numList
    numListB = numList
  else
    !
    ! No transfer coordinates, get corresponding lists the old way
    !
    numList = min(numApoints, numBpoints)
    numListA = numList
    if (pipInput) then
      if (PipGetString('ACorrespondenceList', listString) == 0) &
          call parselist2(listString, listCorrA, numListA, IDIM)
    else
      write (*,113) numList
113   format('Enter a list of points in the first series for which', &
          ' you are sure of the',/,' corresponding point in the', &
          ' second series (Ranges are OK;',/' enter / if the first', &
          i3, ' points are in one-to-one correspondence between', &
          ' the series')
      call rdlist2(5, listCorrA, numListA, IDIM)
    endif
    if (numListA > IDIM) call exitError('Too many points for arrays')
    if (numListA > numList) call exitError( &
        'You have entered more numbers than the minimum number of points in A and B')
    !
    numListB = numListA
    if (pipInput) then
      if (PipGetString('BCorrespondenceList', listString) == 0) &
          call parselist2(listString, listCorrB, numListB, IDIM)
    else
      write (*,114)
114   format('Enter a list of the corresponding points in the', &
          ' second series ',/,' - enter / for ', $)
      call wrlist(listCorrB, numList)
      call rdlist2(5, listCorrB, numListB, IDIM)
    endif
  endif
  !
  if (numListA < minNumToStart) then
    print *
    print *,'ERROR: SOLVEMATCH - Need at least', minNumToStart, ' points to get started'
    call exit(1)
  endif
  if (numListB .ne. numListA) then
    print *
    print *,'ERROR: SOLVEMATCH - You must have the same number ' &
        //' of entries in each list', 'you made', numListA, ' and', &
        numListB, ' entries for lists '//abText(indA) //' and '//abText(indB)
    call exit(1)
  endif
  !
  ! check legality and build map lists
  !
  ! call wrlist(listcorra, nlista)
  ! call wrlist(listcorrb, nlistb)
  do i = 1, numListA
    if (listCorrA(i) <= 0 .or. listCorrB(i) <= 0) call exitError( &
        'You entered a point number less than or equal to zero')
    if (listCorrA(i) > maxContA) call exitError(' You entered a point number higher'// &
        ' than the number of points in '//abText(indA))
    iptA = icontToPointA(listCorrA(i))
    if (iptA == 0) call exitError(' You entered a point number that is not included'// &
        ' in the points from '//abText(indA))
    if (listCorrB(i) > maxContB) call exitError(' You entered a point number higher'// &
        ' than the number of points in '//abText(indB))
    iptB = icontToPointB(listCorrB(i))
    if (iptB == 0) call exitError(' You entered a point number that is not included'// &
        ' in the points from '//abText(indB))
    if (mapAtoB(iptA) .ne. 0) then
      print *
      print *,'ERROR: SOLVEMATCH - Point #', listCorrA(i), &
          ' IN '//abText(indA) //' referred to twice'
      call exit(1)
    elseif (mapped(iptB) .ne. 0) then
      print *
      print *,'ERROR: SOLVEMATCH - Point #', listCorrB(i), &
          ' IN '//abText(indB) //' referred to twice'
      call exit(1)
    endif
    mapAtoB(iptA) = iptB
    mapped(iptB) = iptA
  enddo
  ! print *,(mapAtoB(i), i=1, numApoints)
  ! print *,(mapped(i), i=1, numBpoints)
  !
  ! Get tilt axis angle, plus try to get tomogram pixel size (delta)
  ! and compute scaling factors, overriden by an entry
  !
  if (pipInput) then
    ierr = PipGetTwoFloats('XAxisTilts', xtilta, xtiltB)
    ifZshifts = 1 - PipGetTwoFloats('ZShiftsToTilt', zShiftA, zShiftB)
    ifAngleOfs = 1 - PipGetTwoFloats('AngleOffsetsToTilt', angleOffsetA, &
        angleOffsetB)
    ierr = PipGetInteger('LocalFitting', localNum)
    ierr = PipGetFloat('CenterShiftLimit', shiftLimit)
    ierr = PipGetFloat('AnisotropicLimit', axisCrit)
    call getDelta('ATomogramOrSizeXYZ', aDelta, nxyz(1, 1))
    call getDelta('BTomogramOrSizeXYZ', bDelta, nxyz(1, 2))
    if (aDelta * aPixelSize > 0.) aScale = aPixelSize / aDelta
    if (bDelta * bPixelSize > 0.) bScale = bPixelSize / bDelta
    ierrFactor = PipGetTwoFloats('ScaleFactors', aScale, bScale)
    ! print *,aDelta, bDelta, aPixelSize, bPixelSize, ierrFactor
    !
    ! Conditions for absolute fiducials are that the pixel sizes be
    ! available (meaning new absolute coordinates) and that either
    ! the tomogram pixel sizes were available too or scale factors
    ! were entered and tomogram sizes were provided by getDelta
    !
    if (numSurf .ne. 0 .and. aPixelSize > 0. .and. bPixelSize > 0. &
        .and. (aDelta > 0. .or. (ierrFactor == 0 .and. &
        aDelta == 0.)) .and. (bDelta > 0 .or. &
        (ierrFactor == 0 .and. bDelta == 0.))) then
      !
      ! Shift the original X and Y to the center of the volume
      ! have to divide size by scale because points aren't scaled up yet
      ! Use the actual fiducial size if it exists
      !
      xcen = 0.5 * nxyz(1, 1) / aScale
      ycen = 0.5 * nxyz(jxyz(2), 1) / aScale
      if (nxFidA > 0 .and. nyFidA > 0) then
        xcen = 0.5 * nxFidA
        ycen = 0.5 * nyFidA
      endif
      do i = 1, numApoints
        pointsA(1, i) = pointsA(1, i) - xcen
        pointsA(2, i) = pointsA(2, i) - ycen
      enddo
      !
      xcen = 0.5 * nxyz(1, 2) / bScale
      ycen = 0.5 * nxyz(jxyz(2), 2) / bScale
      if (nxFidB > 0 .and. nyFidB > 0) then
        xcen = 0.5 * nxFidB
        ycen = 0.5 * nyFidB
      endif
      do i = 1, numBpoints
        pointsB(1, i) = pointsB(1, i) - xcen
        pointsB(2, i) = pointsB(2, i) - ycen
      enddo
      relativeFids = .false.
    endif
  else
    write(*,'(1x,a,$)') 'Tilts around the X-axis applied in generating tomograms A and B: '
    read(5,*) xtilta, xtiltB
  endif
  !
  ! Adjust the positions of the points for x axis tilt, scaling, angle
  ! offset and z shift when building the tomogram
  !
  call rotateFids(pointsA, numApoints, xtilta, aScale, angleOffsetA, zShiftA)
  call rotateFids(pointsB, numBpoints, xtiltB, bScale, angleOffsetB, zShiftB)
  !
40 if (pipInput) then
    ierr = PipGetFloat('MaximumResidual', stopLimit)
  else
    write(*,'(1x,a,/,a,$)') 'Maximum residual value above which '// &
        'this program should', ' exit with an error: '
    read(5,*) stopLimit
  endif
  !
  ! if no fiducials, now skip to model entry section
  !
  if (numApoints == 0) go to 50
  !
  ! fill array for regression
  !
  do ia = 1, numApoints
    if (mapAtoB(ia) .ne. 0) then
      numData = numData + 1
      do j = 1, 3
        xMat(j, numData) = pointsB(jxyz(j), mapAtoB(ia))
        xMat(j + 5, numData) = pointsA(jxyz(j), ia)
      enddo
      xMat(4, numData) = 1.
      indOrig(numData) = icontA(ia)
      ! write(*,'(2i4,6f8.1)')ia,mapAtoB(ia),(xMat(j, numData), j=1, 3),  &
      !    (xMat(j, numData), j=6, 8)
    endif
  enddo
  !
  print *,numData, ' pairs of fiducial points originally specified'
  if (.not. pipInput) then
    numSurf = 2
    write(*,'(1x,a,/,a,/,a,/,a,$)') 'Enter 0 to solve for '// &
        'displacements using matching model files, or -1, 1 or 2' &
        , '  to solve only for 3x3 matrix (2 if fiducials are on 2' &
        //' surfaces, 1 if they are', &
        '  on one surface and tomograms are NOT inverted,'// &
        ' or -1 if fiducials are on one', '  surface' &
        //' and tomograms ARE inverted relative to each other): '
    read(5,*) numSurf
    if (abs(numSurf) == 1) then
      invertEntered = .true.
      invertedInDepth = (1 - numSurf) / 2
    endif
  endif
  !
  ! Enter matching models if nsurf is 0
  !
50 if (numSurf == 0) then
    iofs = numColFit + 1
    do model = 1, 2
      !
      ! DNM 7/20/02: Changed to just get the nx, ny, nz and not the
      ! origin from the image file; to use model header information
      ! to scale back to image index coordinates; and to not use or mess
      ! up the y-z transposition variable
      !
      if (.not.pipInput) write(*,'(1x,a,i2,a)') 'Enter NX, NY, NZ of '// &
          'tomogram', model, ', or name of tomogram file'
      call get_nxyz(pipInput, tomoOption(model), 'SOLVEMATCH', 1, nxyz(1, model))
      ! print *,(nxyz(i, model), i=1, 3)
      if (pipInput) then
        if (PipGetString(modelOption(model), filename) > 0) then
          print *
          print *,'ERROR: SOLVEMATCH - No matching model for tomogram', model
          call exit(1)
        endif
      else
        write(*,'(1x,a,i2,a,$)') 'Name of model file from tomogram', model, ': '
        read(*,'(a)') filename
      endif
      if (.not.readw_or_imod(filename)) call exitError('Reading model file')

      ierr = getImodHead(xyScale, zScale, xOffset, yOffset, zOffset, ifFlip)
      ierr = getImodScales(xImScale, yImScale, zImScale)

      numModPts = 0
      do iobj = 1, max_mod_obj
        do ip = 1, npt_in_obj(iobj)
          ipt = abs(object(ibase_obj(iobj) + ip))
          numModPts = numModPts + 1
          if (numData + numModPts > IDIM) call exitError('Too many points for arrays')
          xMat(1 + iofs, numData + numModPts) = (p_coord(1, ipt) - xOffset) / xImScale &
              - 0.5 * nxyz(1, model)
          xMat(2 + iofs, numData + numModPts) = (p_coord(2, ipt) - yOffset) / yImScale &
              - 0.5 * nxyz(2, model)
          xMat(3 + iofs, numData + numModPts) = (p_coord(3, ipt) - zOffset) / zImScale &
              - 0.5 * nxyz(3, model)
          xMat(4, numData + numModPts) = 0.
          indOrig(numData + numModPts) = -numModPts
          if (numApoints == 0) indOrig(numData + numModPts) = numModPts
        enddo
      enddo
      if (model == 1) ndata = numModPts
      iofs = 0
    enddo
    if (numModPts .ne. ndata) call exitError( &
        '# of points does not match between matching models')
    print *,numModPts, ' point pairs from models'
    iofs = numColFit + 1
  else
    !
    ! If no model points
    ! get rid of dummy column: fit 3 columns, pack dependent vars to
    ! the left, and set the offset for adding more dependent var data
    !
    numModPts = 0
    numColFit = 3
    do i = 1, numData
      do j = 5, 7
        xMat(j, i) = xMat(j + 1, i)
      enddo
    enddo
    iofs = 4
    !
    ! if only one surface, "fix" column 2, encode sign in icolfix
    !
    if (abs(numSurf) == 1) then
      icolFixed = 2
      if (invertEntered .and. invertedInDepth > 0) icolFixed = -2
    endif
  endif
  numData = numData + numModPts
  !
  ! loop to add points that weren't initially indicated - as long as
  ! there are any left to add and the last minimum distance was still
  ! low enough
  !
  ifAdded = 0
  addCrit = addRatio * stopLimit
  distMin = 0.
  !
  do while(numData - numModPts < min(numApoints, numBpoints) .and. distMin < addCrit)
    call do3multr(xMat, MAT_SIZE, numData, numColFit, numData, icolFixed, a, delXYZ, &
        cenMeanLoc, devAvg, devSd, devMax, ipntMax, devXyzMax)
    distMin = 1.e10
    !
    ! apply to each point in B that is not mapped to
    !
    do iptB = 1, numBpoints
      if (mapped(iptB) == 0) then
        do ixyz = 1, 3
          pointRot(ixyz) = delXYZ(ixyz)
          if (numColFit == 4) pointRot(ixyz) = delXYZ(ixyz) + a(ixyz, 4)
          do j = 1, 3
            pointRot(ixyz) = pointRot(ixyz) + &
                a(ixyz, j) * pointsB(jxyz(j), iptB)
          enddo
        enddo
        !
        ! search through points in A that don't have map
        !
        do iptA = 1, numApoints
          if (mapAtoB(iptA) == 0) then
            dx = pointsA(jxyz(1), iptA) - pointRot(1)
            if (abs(dx) < distMin) then
              dy = pointsA(jxyz(2), iptA) - pointRot(2)
              if (abs(dy) < distMin) then
                dz = pointsA(jxyz(3), iptA) - pointRot(3)
                if (abs(dz) < distMin) then
                  dist = sqrt(dx**2 + dy**2 + dz**2)
                  if (dist < distMin) then
                    iptAatMin = iptA
                    iptBatMin = iptB
                    distMin = dist
                  endif
                endif
              endif
            endif
          endif
        enddo
      endif
    enddo
    !
    ! add the closest fitting point-pair
    !
    ! print *,iptBatMin, iptAatMin, distmin, devavg, devmax
    if (distMin < addCrit) then
      ifAdded = 1
      numData = numData + 1
      mapAtoB(iptAatMin) = iptBatMin
      mapped(iptBatMin) = iptAatMin
      do j = 1, 3
        xMat(j, numData) = pointsB(jxyz(j), iptBatMin)
        xMat(j + iofs, numData) = pointsA(jxyz(j), iptAatMin)
      enddo
      xMat(4, numData) = 1.
      indOrig(numData) = icontA(iptAatMin)
    endif
  enddo
  !
  print *,numData, ' pairs of points are available for fitting'
  if (ifAdded .ne. 0 .or. numTransCoord > 0) then
    !
    ! rebuild lists of actual contour numbers
    !
    numListA = 0
    do i = 1, numApoints
      if (mapAtoB(i) .ne. 0) then
        numListA = numListA + 1
        listCorrA(numListA) = icontA(i)
        listCorrB(numListA) = icontB(mapAtoB(i))
      endif
    enddo
    print *,'In the final list of correspondences used for', &
        ' fits, points from '//abText(indA) //' are:'
    call wrlist(listCorrA, numListA)
    print *,'Points from '//abText(indB) //' are:'
    call wrlist(listCorrB, numListA)
  endif

  !write(*,105) ((xMat(i, j), i=1, 4), (xMat(i, j), i=1+iofs, 3+iofs), j=1, numData)
  !105 format(7f9.2)
  !
  ! Now figure out whether Z inversion is needed or not for numSurf = +/-1
  if (abs(numSurf) == 1 .and. .not. invertEntered) then
    call do3multr(xMat, MAT_SIZE, numData, numColFit, numData, icolFixed, a, delXYZ, &
        cenMeanLoc, devAvg, devSd, devMax, ipntMax, devXyzMax)
    determPositive = determ3(a)
    call do3multr(xMat, MAT_SIZE, numData, numColFit, numData, -icolFixed, a, delXYZ, &
        cenMeanLoc, devAvg, devSd, devMax, ipntMax, devXyzMax)
    determInverted = determ3(a)
    if (determPositive > 0. .eqv. determInverted > 0.)  &
        call exitError('Cannot determine if there is an inversion in the depth '// &
        'dimension; add InvertedInDepth option')
    if (determInverted > 0.) icolFixed = -icolFixed
  endif
  !
  maxDrop = nint(0.1 * (numData - 1))
  critProb = 0.01
  elimMin = 3.
  absProbCrit = 0.002
  call solve_wo_outliers(xMat, MAT_SIZE, numData, numColFit, icolFixed, maxDrop, critProb, &
      absProbCrit, elimMin, idrop, numDrop, a, delXYZ, cenMeanLoc, devAvg, devSd, &
      devMax, ipntMax, devXyzMax)
  !
  if (numDrop .ne. 0) then
    write(*,104) numDrop, devAvg, devSd, abText(indA), (indOrig(idrop(i)), i = 1, numDrop)
    write(*,115) (xMat(numColFit + 1, i), i = numData + 1 - numDrop, numData)
104 format(/,i3,' points dropped by outlier elimination; ', &
        'residual mean =',f7.2,', SD =',f7.2,/, &
        ' point # in ',a1,':',(9i7))
115 format(' deviations  :',(9f7.1))
  endif
  !
  indOrigAtMax = indOrig(ipntMax)
  if (indOrigAtMax <= maxContA .and. modObj(1, 1) > 0) then
    iptA = icontToPointA(indOrigAtMax)
    !
    ! BRT is using 'Mean residual' as tag
    write(*,1011) devAvg, devMax, indOrigAtMax, modObj(iptA, 1), &
        modCont(iptA, 1), abText(indA)
1011 format(//,' Mean residual',f8.3,',  maximum',f9.3, &
        ' at point #',i4,' (Obj',i3,' cont',i4,' in ',a1,')')
  else
    write(*,1012) devAvg, devMax, indOrigAtMax, abText(indA)
1012 format(//,' Mean residual',f8.3,',  maximum',f9.3, &
        ' at point #',i4,' (in ',a1,')')
  endif
  write(*,1013) (devXyzMax(i), i = 1, 3)
1013 format('  Deviations:',3f9.3)
  !
  if (numColFit > 3) then
    !
    ! fit to both: report the dummy variable offset, make sure that
    ! the dxyz are non-zero for matchshifts scanning
    !
    write(*,103) (a(i, 4), i = 1, 3)
103 format(/,' X, Y, Z offsets for fiducial dummy variable:',3f10.3)
    if (abs(delXYZ(1)) < 0.001 .and. abs(delXYZ(2)) < 0.001 .and. &
        abs(delXYZ(3)) < 0.001) delXYZ(1) = 0.0013
    !
  else if (numModPts == 0 .and. relativeFids) then
    !
    ! No matching models and relative fiducials: set dxyz to zero
    ! and report the offset from the fit
    !
    write(*,107) (delXYZ(i), i = 1, 3)
107 format(/,' X, Y, Z offsets from fiducial fit:',3f10.3)
    delXYZ(1) = 0.
    delXYZ(2) = 0.
    delXYZ(3) = 0.
    if (imodGetenv('SOLVEMATCH_TEST', filename) .ne. 0) call exitError( &
        'Relative fiducial coordinates are no longer allowed; rerun '// &
        'Tiltalign for both axes to get absolute coordinates')
  else
    !
    ! Absolute fiducials: just make sure they are not exactly zero
    !
    if (abs(delXYZ(1)) < 0.001 .and. abs(delXYZ(2)) < 0.001 .and. &
        abs(delXYZ(3)) < 0.001) delXYZ(1) = 0.0013
  endif
  !
  print *
  print *,'Transformation matrix for matchvol:'
  write(*,102) ((a(i, j), j = 1, 3), delXYZ(i), i = 1, 3)
102 format(3f10.6,f10.3)
  !
  ! Compute scaling of vectors and report it
  ! Then give warnings of unequal scalings
  do j = 1, 3
    sumSq = 0.
    do i = 1, 3
      sumSq = sumSq + a(i, j)**2
    enddo
    axisScale(j) = sqrt(sumSq)
  enddo
  !
  ! BRT is using 'Scaling along' as tag
  write(*,118) (axisScale(i), i = 1, 3)
118 format(/, 'Scaling along the three axes - X:',f7.3,'  Y:',f7.3,'  Z:',f7.3)
  xyScaleDiff = 100. * abs((axisScale(1) - axisScale(2)) / axisScale(1))
  xzScaleDiff = 100. * abs((axisScale(1) - axisScale(3)) / axisScale(1))
  yzScaleDiff = 100. * abs((axisScale(2) - axisScale(3)) / axisScale(3))
  badAxis1 = ' '
  badAxis2 = ' '
  if (axisCrit > 0) then
    if (xyScaleDiff > axisCrit .and. xzScaleDiff > axisCrit .and. &
        yzScaleDiff > axisCrit) then
      reportDiff = min(xyScaleDiff, xzScaleDiff, yzScaleDiff)
      write(*,'(/,a,f3.0,a)') 'WARNING: The scalings along all three axes'// &
          ' differ from each other by more than', reportDiff, '%'
      badAxis1 = 'A'
    elseif (xyScaleDiff > axisCrit .and. xzScaleDiff > axisCrit) then
      badAxis1 = 'X'
      reportDiff = 0.5 * (xyScaleDiff + xzScaleDiff)
    elseif (xyScaleDiff > axisCrit .and. yzScaleDiff > axisCrit) then
      badAxis1 = 'Y'
      reportDiff = 0.5 * (xyScaleDiff + yzScaleDiff)
    elseif (xzScaleDiff > axisCrit .and. yzScaleDiff > axisCrit) then
      badAxis1 = 'Z'
      reportDiff = 0.5 * (xzScaleDiff + yzScaleDiff)
    elseif (xyScaleDiff > axisCrit) then
      badAxis1 = 'X'
      badAxis2 = 'Y'
      reportDiff = xyScaleDiff
    elseif (xzScaleDiff > axisCrit) then
      badAxis1 = 'X'
      badAxis2 = 'Z'
      reportDiff = xzScaleDiff
    elseif (yzScaleDiff > axisCrit) then
      badAxis1 = 'Y'
      badAxis2 = 'Z'
      reportDiff = yzScaleDiff
    endif
  endif
  if (badAxis2 .ne. ' ') then
    write(*,'(/,a,a,a,a,a,f4.0,a)') 'WARNING: The scaling along the ', &
        badAxis1, ' and ', badAxis2, ' axes differ by', reportDiff, '%'
  elseif (badAxis1 .ne. 'A' .and. badAxis1 .ne. ' ') then
    write(*,'(/,a,a,a,f4.0,a)') 'WARNING: The scaling along the ', badAxis1, &
        ' axis differs from the other two axes by', reportDiff, '%'
  endif
  !
  ! Issue specific dual-axis warning with advice
  ! BRT is looking for 'Try specifying' and 'on one surface'
  if (badAxis1 == 'Y' .and. numSurf == 2 .and. (ifAngleOfs .ne. 0 .or. &
      ifZshifts .ne. 0)) write(*,119)
119 format('WARNING: Y scaling is probably wrong because you specified that', &
      ' points are on',/,'WARNING:    two surfaces but there are too few ', &
      'points on one surface.',/,'WARNING:    Try specifying that points ', &
      'are on one surface')

  filename = ' '
  if (pipInput) then
    if (PipGetInOutFile('OutputFile', 1, ' ', filename) .ne. 0)  &
        call exitError('No output file specified')
  else
    print *,'Enter name of file to place transformation in, or ', 'Return for none'
    read(5, '(a)') filename
  endif
  if (filename .ne. ' ') then
    call dopen(1, filename, 'new', 'f')
    write(1, 102) ((a(i, j), j = 1, 3), delXYZ(i), i = 1, 3)
    close(1)
  endif
  !
  ierr = 0
  limRaised = devMax + 1.2
  devAllMax = 0.
  if (devMax > stopLimit .and. numModPts == 0 .and. localNum > 0) then
    if (localNum < 6 .or. numData <= localNum) then
      if (localNum < 6) write(*,'(/,a)') 'ERROR: SOLVEMATCH - Local'// &
          ' fits must have a minimum of 6 points'
      if (numData <= localNum) write(*,'(/,a,/,a)') 'Local fitting is not'// &
          ' available because the number of matched points', &
          ' is no bigger than the minimum for local fitting'
    else
      !
      ! For local fits, get extent of the data
      !
      xmin = 1.e10
      xmax = -xmin
      ymin = 1.e10
      ymax = -ymin
      do ia = 1, numApoints
        if (mapAtoB(ia) > 0) then
          xmin = min(xmin, pointsA(1, ia))
          ymin = min(ymin, pointsA(2, ia))
          xmax = max(xmax, pointsA(1, ia))
          ymax = max(ymax, pointsA(2, ia))
        endif
      enddo
      ! print *,'minmax', xmin, xmax, ymin, ymax
      !
      ! Set up number and interval between local areas
      !
      targetSize = sqrt(localNum * (xmax - xmin) * (ymax - ymin) / numData)
      numLocalX = max(1., 2. * (xmax - xmin) / targetSize + 0.1)
      numLocalY = max(1., 2. * (ymax - ymin) / targetSize + 0.1)
      dxLocal = 0
      dyLocal = 0
      if (numLocalX > 1) &
          dxLocal = (xmax - xmin - targetSize) / (numLocalX - 1)
      if (numLocalY > 1) &
          dyLocal = (ymax - ymin - targetSize) / (numLocalY - 1)
      !
      ! Loop on local areas, getting mean residual and number with max
      ! above limit
      !
      sumMean = 0.
      sumMax = 0.
      numBig = 0
      do ixl = 1, numLocalX
        xcen = xmin + targetSize / 2. + (ixl - 1) * dxLocal
        do iyl = 1, numLocalY
          ycen = ymin + targetSize / 2. + (iyl - 1) * dyLocal
          size = targetSize
          call fillLocalData(xcen, ycen, size, localNum, pointsA, pointsB, &
              numApoints, mapAtoB, jxyz, xMat, MAT_SIZE, numData)
          ! do j = 1, ndat
          ! write(*,'(6f8.1)') (xMat(i, j), i=1, 3), (xMat(i, j), i=5, 7)
          ! enddo
          maxDrop = nint(0.1 * numData)
          if (numData <= 6) maxDrop = 0
          call solve_wo_outliers(xMat, MAT_SIZE, numData, numColFit, icolFixed, maxDrop, &
              critProb, absProbCrit, elimMin, idrop, numDrop, amatLocal, dxyzLocal, &
              cenMeanLoc, devAvgLoc, devSd, devMaxLoc, ipntMax, devXyzMax)
          ! print *,xcen, ycen, size, ndat, devavgloc, devmaxloc
          sumMean = sumMean + devAvgLoc
          sumMax = sumMax + devMaxLoc
          devAllMax = max(devAllMax, devMaxLoc)
          if (devMaxLoc > stopLimit) numBig = numBig + 1
        enddo
      enddo
      !
      ! BRT is using 'Local fits' and 'Average mean' as tags
      write(*,1015) localNum, sumMean / (numLocalX * numLocalY), &
          sumMax / (numLocalX * numLocalY), devAllMax, numBig, &
          numLocalX * numLocalY, stopLimit
1015  format(/,'Local fits to a minimum of',i4,' points give:',/, &
          '   Average mean residual',f8.2,/,'   Average max residual',f8.2, &
          /,'   Biggest max residual',f9.2,/,i5,' of',i5, &
          ' local fits with max residual above',f8.1)

      if (devAllMax < 1.5 * stopLimit .and. &
          numBig <= 0.05 * numLocalX * numLocalY) then
        write(*,'(/,a,/)') 'The local fits indicate that THIS SOLUTION IS GOOD ENOUGH'
        devMax = stopLimit - 1
      endif
      limRaised = devAllMax + 1.2

      if (shiftLimit > 0.) then
        !
        ! Get the data for the central area and transform the points by
        ! the global transformation and measure shift in the points.
        ! 3/19/07: switched to this from doing a fit because the fit can
        ! screw up the shifts if there are not points on both sides
        xcen = 0.
        ycen = 0.
        size = targetSize
        call fillLocalData(xcen, ycen, size, localNum, pointsA, pointsB, numApoints, &
            mapAtoB, jxyz, xMat, MAT_SIZE, numData)
        dxyzLocal(1) = 0.
        dxyzLocal(2) = 0.
        dxyzLocal(3) = 0.
        do ip = 1, numData
          do j = 1, 3
            xMat(10 + j, ip) = delXYZ(j)
            do i = 1, 3
              xMat(10 + j, ip) = xMat(10 + j, ip) + a(j, i) * xMat(i, ip)
            enddo
            dxyzLocal(j) = dxyzLocal(j) +  &
                (xMat(10 + j, ip) - xMat(numColFit + 1 + j, ip)) / numData
          enddo
          ! write(*,132) (xMat(j, ip), j=1, 3), (xMat(j, ip), j=5, 7), (xMat(j,ip),j=11,13)
          !132           format(9f8.1)
        enddo
        cenShiftX = dxyzLocal(1)
        cenShiftY = dxyzLocal(2)
        cenShiftZ = dxyzLocal(3)
        dist = sqrt(cenShiftX**2 + cenShiftY**2 + cenShiftZ**2)
        ! print *,'distance', dist, csdx, csdy, csdz
        if (dist >= shiftLimit) then
          !
          ! BRT is looking for 'InitialShiftXYZ' and 'needs'
          write(*, 1016) dist, nint(cenShiftX), nint(cenShiftY), nint(cenShiftZ), &
              nint(cenShiftX), nint(cenShiftZ), nint(cenShiftY)
1016      format(/,'Center shift indicated by local fit is',f6.0, &
              ', bigger than the specified limit',/, &
              '   The InitialShiftXYZ for corrsearch3d needs to be',3i5,/, &
              '   In Etomo, set Patchcorr Initial shifts in X, Y, Z to',3i5)
          if (nxyz(jxyz(3), 2) > nxyz(jxyz(3), 1)) &
              write(*,1017) nxyz(2, jxyz(3))
1017      format('   You should also set thickness of initial ', &
              'matching file to at least', i5,/, &
              '     (In Etomo, Initial match size for Matchvol1)')
          !
          ! BRT is looking for 'CenterShiftLimit' and 'avoid stopping'
          write(*,1018) nint(dist) + 1
1018      format('   To avoid stopping with this error, set CenterShift', 'Limit to',i4, &
              /, '     (In Etomo, Limit on center shift for Solvematch)')
          if (devMax < stopLimit) then
            !
            ! BRT is looking for INITIAL SHIFT' and 'SOLUTION IS OK'
            write(*,'(/,a)') 'ERROR: SOLVEMATCH - Initial shift needs to'// &
                ' be set for patch correlation (but solution is OK)'
          else
            write(*,'(/,a)') 'ERROR: SOLVEMATCH - Initial shift needs to'// &
                ' be set for patch correlation'
          endif
          ierr = 1
        endif
      endif
    endif
  endif


  if (devMax > stopLimit) then
    !
    ! Give some guidance based upon the ratios between max and mean
    ! deviation and max deviation and stopping limit
    !
    write(*,106) devMax
106 format(/, 'The maximum residual is',f8.2,', too high to proceed')
    loMaxAvgRatio = 4.
    hiMaxAvgRatio = 12.
    loMaxLimRatio = 2.
    hiMaxLimRatio = 3.
    globLocAvgRatio = 3.
    if ((devMax < loMaxAvgRatio * devAvg .and. &
        devMax < loMaxLimRatio * stopLimit) .or. (devAllMax > 0 .and. &
        devAvgLoc * globLocAvgRatio < devAvg .and. &
        devAllMax < loMaxLimRatio * stopLimit)) then
      if (numTransCoord > 0) then
        write(*,1115) limRaised
1115    format('Since corresponding points were picked using coordinates ', &
            'from transferfid,',/,'this is almost certainly due to ', &
            'distortion between the volumes,',/, &
            ' and you should just raise the residual limit to',i4)
      else
        if (devAllMax <= 0.) then
          write(*,111) loMaxAvgRatio, devAvg, loMaxLimRatio, limRaised
        else
          write(*,1116) globLocAvgRatio, devAllMax, loMaxLimRatio, limRaised
        endif
111     format('Since the maximum residual is less than',f6.1, &
            ' times the mean residual (', f8.2,')',/,' and less than', &
            f6.1, ' times the specified residual limit,',/, &
            ' this is probably due to distortion between the volumes,',/, &
            'and you should probably just raise the residual limit to',i4)
1116    format('Since the local fits improved the mean residual by more', &
            ' than a factor of',f6.1,/,' and the local maximum residual (', &
            f8.2,') is less than', f6.1, ' times the',/, &
            ' specified residual limit, this is probably due to ', &
            'distortion between the',/, ' volumes, and you should ', &
            'probably just raise the residual limit to',i4)
      endif
    elseif (numTransCoord > 0) then
      write(*,112) limRaised, indOrigAtMax
112   format('Since corresponding points were picked using coordinates ', &
          'from transferfid,',/,'this is probably due to distortion ', &
          'between the volumes.',/,'You could raise the residual limit to', &
          i4,' or start with a subset of points.',/,'Bad correspondence ', &
          'is unlikely but you could check points (especially',i4,').')
    elseif (devMax > hiMaxAvgRatio * devAvg .or. &
        devMax > hiMaxLimRatio * stopLimit) then
      if (devMax > hiMaxAvgRatio * devAvg) &
          write(*,108) hiMaxAvgRatio, devAvg
108   format('The maximum residual is more than', f6.1, &
          ' times the mean residual (', f8.2,')')
      if (devMax > hiMaxLimRatio * stopLimit) &
          write(*,109) hiMaxLimRatio
109   format('The maximum residual is more than', f6.1, &
          ' times the specified residual limit')
      print *,'This is probably due to a bad correspondence list.'
      write (*,110) indOrigAtMax
110   format('Check the points (especially',i4, ') or start with a subset of the list')
    else
      print *,'The situation is ambiguous but could be due to a', &
          ' bad correspondence list.'
      write (*,110) indOrigAtMax
    endif
    call exitError('Maximum residual is too high to proceed')
  endif
  call exit(ierr)
end program solvematch


! GETDELTA returns the "delta" or pixel size from a tomogram as well
! as the nxyz
! OPTION is the option that specifies the tomogram or nx, ny, nz
! DELTA is returned with -1 if the option was not entered at all,
! 0 if sizes were entered, or the actual delta value from the file
!
subroutine getDelta(option, delta, nxyz)
  implicit none
  integer*4 nxyz(3), mxyz(3), mode, i, PipGetString
  character*(*) option
  character*80 line
  real*4 dmin, dmax, dmean, delTmp(3), delta
  logical line_is_filename
  !
  delta = -1.
  if (PipGetString(option, line) > 0)  return
  delta = 0.
  if (line_is_filename(line)) go to 10
  read(line,*,err = 10, end = 10) (nxyz(i), i = 1, 3)
  return
  !
10 call ialprt(.false.)
  call imOpen(5, line, 'ro')
  call irdhdr(5, nxyz, mxyz, mode, dmin, dmax, dmean)
  call iiuRetDelta(5, delTmp)
  call iiuClose(5)
  call ialprt(.true.)
  delta = delTmp(1)
  return
end subroutine getDelta


! getFiducials reads the fiducials from FILENAME, storing the
! sequential point number in ICONTA, the X, Y, Z coordinates in PNTA,
! and the number of points in NPNTA
! IDIM specifies the dimensions of the arrays
! AXIS specifies the axis for a message
! If a pixel size is found on the first line, it is returned in
! PIXELSIZE, otherwise this is set to 0
!
subroutine getFiducials(filename, icontA, pointsA, numApoints, modObj, &
    modCont, IDIM, axis, pixelSize, nxFidDim, nyFid)
  implicit none
  character*(*) filename
  character*(*) axis
  character*100 line
  integer*4 icontA(*), IDIM, numApoints, len, i, j, nxFidDim, nyFid, unmFields
  integer*4 modObj(*), modCont(*), numeric(10)
  real*4 pointsA(3, IDIM), pixelSize, xnumIn(10)
  !
  numApoints = 0
  pixelSize = 0.
  nxFidDim = 0
  nyFid = 0
  call dopen(1, filename, 'old', 'f')
  !
  ! get the first line and search for PixelSize: (first version)
  ! or Pix: and Dim: (second version)
  !
  read(1, '(a)', err = 10, end = 15) line
  len = len_trim(line)
  i = 1
  do while (i < len - 11)
    if (line(i:i + 10) == 'Pixel size:') &
        read(line(i + 11:len),*, err = 8) pixelSize
    if (line(i:i + 3) == 'Pix:') &
        read(line(i + 4:len),*, err = 8) pixelSize
    if (line(i:i + 3) == 'Dim:') &
        read(line(i + 4:len),*, err = 8) nxFidDim, nyFid
    i = i + 1
  enddo
  !
  ! process first line then loop until end
  !
  do while(.true.)
    numApoints = numApoints + 1
    if (numApoints > IDIM) call exitError('Too many points for arrays')
    call frefor2(line, xnumIn, numeric, unmFields, 6)
    if (unmFields == 6 .and. numeric(5) == 1) then
      read(line, *) icontA(numApoints), (pointsA(j, numApoints), j = 1, 3), &
          modObj(numApoints), modCont(numApoints)
    else
      read(line, *) icontA(numApoints), (pointsA(j, numApoints), j = 1, 3)
      modObj(numApoints) = 0
      modCont(numApoints) = 0
    endif
    !
    ! invert the Z coordinates because the tomogram built by TILT is
    ! inverted relative to the solution produced by TILTALIGN
    ! 12/1/09: THIS IS NOT TRUE AFTER ACCOUNTING FOR FLIPPING IN 3dMOD
    pointsA(3, numApoints) = -pointsA(3, numApoints)
    read(1, '(a)', end = 15, err = 10) line
  enddo
  !
15 print *,numApoints, ' points from ', axis
  close(1)
  return
10 call exitError('Reading fiducial point file')
8 call exitError('Reading pixel size or dimensions from point file')
end subroutine getFiducials


! readFidModel gets the filename with the given OPTION, reads fiducial
! model, finds points with Z = IZBEST, and returns their X/Y coords
! in FIDMODX, FIDMODY and their object and contour numbers in MODOBJFID
! and MODCONTFID.  NUMFID is the number of fiducials returned; iDIM
! specifies the dimensions of the arrays; ABTEXT is A or B.  transPixel
! is a pixel size for transfer coords if non-zero.
!
subroutine readFidModelFile(option, izBest, abText, fidModX, &
    fidModY, modObjFid, modContFid, numFid, IDIM, transPixel)
  implicit none
  include 'model.inc90'
  character*(*) option
  character*1 abText
  integer*4 izBest, IDIM, numFid, modObjFid(*), modContFid(*)
  real*4 fidModX(*), fidModY(*), fidScale, transPixel
  real*4 xImScale, yImScale, zImScale
  character*160 filename
  integer*4 iobj, ip, ipt
  logical*4 looking
  logical readw_or_imod
  integer*4 PipGetString, getImodScales

  if (PipGetString(option, filename) .ne. 0) call exitError('Fiducial '// &
      'models for both axes must be entered to use transfer coords')
  if (.not.readw_or_imod(filename)) call exitError( &
      'Reading fiducial model for axis '//abText)
  !
  ! If a pixel size was defined for transfer, scale coords to match that
  fidScale = 1.
  ip = getImodScales(xImScale, yImScale, zImScale)
  if (ip == 0 .and. transPixel > 0.) fidScale = xImScale / transPixel
  call scale_model(0)
  numFid = 0
  do iobj = 1, max_mod_obj
    looking = .true.
    ip = 1
    !
    ! Find point at best z value, record it and its obj/cont
    !
    do while (looking .and. ip <= npt_in_obj(iobj))
      ipt = abs(object(ibase_obj(iobj) + ip))
      if (nint(p_coord(3, ipt)) == izBest) then
        numFid = numFid + 1
        if (numFid > IDIM) call exitError('Too many fiducial '// &
            'contours for arrays in model for axis '//abText)
        fidModX(numFid) = p_coord(1, ipt) * fidScale
        fidModY(numFid) = p_coord(2, ipt) * fidScale
        call objToCont(iobj, obj_color, modObjFid(numFid), &
            modContFid(numFid))
        ! print *, modObjFid(numFid), modContFid(numFid), fidModX(numFid), &
        ! fidModY(numFid)
        looking = .false.
      endif
      ip = ip + 1
    enddo
  enddo
  return
end subroutine readFidModelFile


! ROTATEFIDS rotates, scales and shifts the NPNTA fiducials in PNTA
! XTILT is the x axis tilt, aScale is the scaling, angleOffset
!
subroutine rotateFids(points, numPoints, xtilt, scale, angleOffset, zShift)
  implicit none
  real*4 points(3,*), xtilt, scale, angleOffset, zShift
  real*4 cosa, cosb, sinAlpha, sinBeta, ytmp
  integer*4 numPoints, i
  real*4 sind, cosd
  !
  ! use the negative of the angle to account for the inversion of the
  ! tomogram; rotate the 3-d points about the Y axis first, then the Z
  ! axis; add Z shift
  ! 12/1/09: THE NEGATIVE COMPENSATES FOR THE FLIPPING OF Z ABOVE, BUT
  ! ARE THERE OTHER CONSEQUENCES?  IS ALL THIS RIGHT?
  cosa = cosd(-xtilt) * scale
  sinAlpha = sind(-xtilt) * scale
  sinBeta = sind(-angleOffset)
  cosb = cosd(-angleOffset)
  do i = 1, numPoints
    ytmp = cosb * points(1, i) - sinBeta * points(3, i)
    points(3, i) = sinBeta * points(1, i) + cosb * points(3, i)
    points(1, i) = ytmp
    ytmp = cosa * points(2, i) - sinAlpha * points(3, i)
    points(3, i) = sinAlpha * points(2, i) + cosa * points(3, i) + zShift
    points(2, i) = ytmp
    points(1, i) = points(1, i) * scale
  enddo
  return
end subroutine rotateFids


! fillLocalData fills the data array XMAT with at least MINDAT points
! from an area centered at XCEN, YCEN.  SIZE is called with an initial
! trial size, and returned with the final size needed to include the
! required points.  PNTA and PNTB have the points, NPNTA is the total
! number in PNTA, MAPAB is the mapping to indices in PNTB, JXYZ is the
! dimension mapping.  NDAT is returned with number of points.
!
subroutine fillLocalData(xcen, ycen, size, minData, pointsA, pointsB, numApoints, &
    mapAtoB, jxyz, xMat, MAT_SIZE, numData)
  implicit none
  integer*4 numApoints, mapAtoB(*), jxyz(3), MAT_SIZE, numData, ia, j, minData
  real*4 xcen, ycen, size, pointsA(3,*), pointsB(3,*), xMat(MAT_SIZE,*)
  real*4 xmin, ymin, xmax, ymax, sizeIn
  !
  sizeIn = size
  numData = 0
  do while (numData < minData)
    numData = 0
    xmin = xcen - size / 2.
    xmax = xcen + size / 2.
    ymin = ycen - size / 2.
    ymax = ycen + size / 2.
    do ia = 1, numApoints
      if (mapAtoB(ia) > 0 .and. pointsA(1, ia) >= xmin .and. &
          pointsA(1, ia) <= xmax .and. pointsA(2, ia) >= ymin &
          .and. pointsA(2, ia) <= ymax) then
        numData = numData + 1
        do j = 1, 3
          xMat(j, numData) = pointsB(jxyz(j), mapAtoB(ia))
          xMat(j + 4, numData) = pointsA(jxyz(j), ia)
        enddo
      endif
    enddo
    if (numData < minData) size = size + 0.02 * sizeIn
  enddo
  return
end subroutine fillLocalData

! Determinant of a 3x3 matrix
real*4 function determ3(a)
  implicit none
  real*4 a(3,3)
  determ3 = a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) + a(1,2)*a(2,3)*a(3,1) -  &
      a(1,2)*a(2,1)*a(3,3) + a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)
  return
end function determ3
