! * * * * * REFINEMATCH * * * * * *
!
! REFINEMATCH will solve for a general 3-dimensional linear
! transformation to align two volumes to each other.  It performs
! multiple linear regression on the  displacements between the volumes
! determined at a matrix of positions.  The displacements must be
! contained in a file with the following form:
!
! Number of displacements
! One line for each displacement consisting of the X, Y, and Z
! .  coordinates in the first volume, then the displacements in X, Y
! .  and Z involved in moving from the first to the second volume
!
! See man page for details
!
! David Mastronarde, 1995
!
! $Id$
!
program refinematch
  implicit none
  integer LIMEXTRA, LIMCRIT
  parameter (LIMEXTRA = 20, LIMCRIT = 20)
  real*4 aMat(3,3), delXYZ(3), devXYZ(3), aMatInit(3,3), delInit(3), prodMat(3,3)
  real*4 devXYZmax(3), fitCenter(3), cxlast(3), freinp(LIMEXTRA), prodDel(3)
  integer*4 nxyz(3), numeric(LIMEXTRA), IDextra(LIMEXTRA)
  real*4 selectCrit(LIMCRIT, LIMEXTRA)
  integer*4 numSelectCrit, indSelect, icolSelect(LIMEXTRA), isignSelect(LIMEXTRA)  
  integer*4, allocatable :: idrop(:), indVert(:), numVerts(:), indToPatch(:)
  real*4, allocatable :: fitMat(:,:), xVerts(:), yVerts(:), contourZ(:), extraVals(:,:)
  real*4, allocatable :: cenXYZ(:,:), vecXYZ(:,:), dropRes(:)
  logical inside, oneLayer(3), haveCCC
  character*320 filename, line, initialFile, productFile
  !
  integer*4 numConts, ierr, ifFlip, i, indY, indZ, indcur, iobj, numData, numFit
  real*4 fracDrop, cccRes, prodScaleFac, cenDist, cenDistMin
  integer*4 ifUse, icont, icontMin, j, ind, maxDrop, ndrop, ipntMax, ip, ipt, numFields
  real*4 dzMin, critProb, elimMin, absProbCrit, devMean, devSD, devMax, stopLim, dz
  integer*4 icolFixed, limCont, limVert, maxExtra, numExtraID, limPatch, matCols
  integer*4 numColSelect, ipat, indClosest
  integer*4 idResidCol/2/
  integer*4 readNumPatches

  logical pipInput
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetInteger
  integer*4 PipGetString, PipGetFloat, PipGetTwoFloats
  integer*4 PipGetInOutFile
  !
  ! fallbacks from ../../manpages/autodoc2man -3 2  refinematch
  !
  integer numOptions
  parameter (numOptions = 17)
  character*(40 * numOptions) options(1)
  options(1) = &
      'patch:PatchFile:FN:@output:OutputFile:FN:@region:RegionModel:FN:@'// &
      'volume:VolumeOrSizeXYZ:FN:@residual:ResidualPatchOutput:FN:@'// &
      'reduced:ReducedVectorOutput:FN:@limit:MeanResidualLimit:F:@'// &
      'extra:ExtraValueSelection:IPM:@select:SelectionCriteria:FAM:@'// &
      'maxfrac:MaxFractionToDrop:F:@minresid:MinResidualToDrop:F:@'// &
      'prob:CriterionProbabilities:FP:@initial:InitialTransformFile:FN:@'// &
      'product:ProductTransformFile:FN:@scale:ScaleShiftByFactor:F:@'// &
      'param:ParameterFile:PF:@help:usage:B:'
  !
  fracDrop = 0.1
  critProb = 0.01
  elimMin = 0.5
  absProbCrit = 0.002
  matCols = 20
  numColSelect = 0
  numSelectCrit = 0
  initialFile = ' '
  productFile = ' '
  !
  call PipReadOrParseOptions(options, numOptions, 'refinematch', &
      'ERROR: REFINEMATCH - ', .true., 3, 1, 1, numOptArg, numNonOptArg)
  pipInput = numOptArg + numNonOptArg > 0
  !
  ! Open patch file
  !
  if (PipGetInOutFile('PatchFile', 1, &
      'Name of file with correlation positions and results', filename) .ne. 0) &
      call exitError('No input patch file specified')
  call dopen(1, filename, 'old', 'f')

  if (.not. pipInput) print *,'Enter file name or NX, NY, NZ of tomogram being matched to'
  call get_nxyz(pipInput, 'VolumeOrSizeXYZ', 'REFINEMATCH', 1, nxyz)

  if (pipInput) then
    filename = ' '
    ierr = PipGetString('RegionModel', filename)
    ierr = PipGetFloat('MaxFractionToDrop', fracDrop)
    ierr = PipGetFloat('MinResidualToDrop', elimMin)
    ierr = PipGetTwoFloats('CriterionProbabilities', critProb, absProbCrit)
    if (PipGetFloat('MeanResidualLimit', stopLim) > 0) call exitError( &
        'You must enter a mean residual limit')
    ierr = PipGetString('InitialTransformFile', initialFile)
    ind = PipGetString('ProductTransformFile', productFile)
    if (ierr + ind == 1) call exitError( &
        'You must enter both -initial and -product files, not just one')
    prodScaleFac = 1.
    ierr = PipGetFloat('ScaleShiftByFactor', prodScaleFac)
    if (initialFile .ne. ' ') then
      call dopen(3, initialFile, 'ro', 'f')
      read(3,*, err = 99) ((aMatInit(i, j), j = 1, 3), delInit(i), i = 1, 3)
      close(3)
    endif
  else
    write(*,'(1x,a,/,a,$)') &
        'Enter name of model file with contour enclosing area to ' &
        //'use,', ' or Return to use all patches: '
    read(*,'(a)') filename
  endif

  numConts = 0
  if (filename .ne. ' ') then
    call getContourArraySizes(filename, 0, limCont, limVert)
    allocate(xVerts(limVert), yVerts(limVert), contourZ(limCont), indVert(limCont),  &
        numVerts(limCont), stat = ierr)
    call memoryError(ierr, 'arrays for boundary contours')
    call get_region_contours(filename, 'REFINEMATCH', xVerts, yVerts, numVerts, &
        indVert, contourZ, numConts, ifFlip, limCont, limVert, 0)
  else
    ifFlip = 0
    if (nxyz(2) < nxyz(3)) ifFlip = 1
  endif
  indY = 2
  if (ifFlip .ne. 0) indY = 3
  indZ = 5 - indY
  !
  IDextra(:) = 0
  if (readNumPatches(1, numData, numExtraID, IDextra, LIMEXTRA) .ne. 0) call exitError &
      ('Reading first line of patch file')
  maxExtra = 0
  do i = 1, numData
    read(1,'(a)') line
    call frefor2(line, freinp, numeric, numFields, LIMEXTRA)
    maxExtra = max(maxExtra, numFields - 6)
  enddo
  rewind(1)
  read(1,'(a)') line

  limPatch = numData + 10
  allocate(fitMat(matCols, limPatch),  cenXYZ(limPatch,3), vecXYZ(limPatch,3), &
      indToPatch(limPatch), idrop(limPatch), dropRes(nint(fracDrop * limPatch) + 10), &
      stat = ierr)
  call memoryError(ierr, 'arrays for data matrix and patches')
  if (maxExtra > 0) then
    allocate(extraVals(maxExtra, limPatch), stat=ierr)
    call memoryError(ierr, 'array for extra values')
    extraVals(:,:) = 0
  endif

  if (pipInput) call getExtraSelections(icolSelect, isignSelect, numColSelect, &
      numSelectCrit, IDextra, maxExtra, selectCrit, LIMEXTRA, LIMCRIT)

  ! Read in the patch data and extra data
  do j = 1, 3
    oneLayer(j) = .true.
    cxlast(j) = 0.
  enddo
  do i = 1, numData
    !
    ! these are center coordinates and location of the second volume
    ! relative to the first volume
    !
    read(1, '(a)') line
    call frefor2(line, freinp, numeric, numFields, LIMEXTRA)
    do j = 1, 3
      cenXYZ(i, j) = freinp(j)
      vecXYZ(i, j) = freinp(j + 3)
      if (i > 1 .and. cenXYZ(i, j) .ne. cxlast(j)) oneLayer(j) = .false.
      cxlast(j) = cenXYZ(i, j)
    enddo
    do j = 7, numFields
      extraVals(j - 6, i) = freinp(j)
    enddo
  enddo
  close(1)
  !
  if (.not.pipInput) then
    write(*,'(1x,a,$)') 'Mean residual above which to STOP and exit with an error: '
    read(5,*) stopLim
  endif
  !
  icolFixed = 0
  do i = 1, 3
    if (icolFixed .ne. 0 .and. oneLayer(i)) call exitError( &
        'Cannot fit to patches that extend in only one dimension')
    if (oneLayer(i)) icolFixed = i
  enddo
  if (icolFixed > 0) print *,'There is only one layer of patches in the ',  &
      char(ichar('W') + icolFixed), ' dimension'
  !
  ! Loop on the selection criteria
  if (numColSelect > 0) write(*,'(i8,a,/)') numData, ' total patches in patch file'
  do indSelect = 1, max(1, numSelectCrit)
    if (numColSelect > 0) write(*,'(a,f9.3)') &
        'Applying extra column selection criteria', selectCrit(indSelect, 1)
    numFit = 0
    do i = 1, numData
      ifUse = 1
      if (numConts > 0) call checkBoundaryConts(cenXYZ(i, 1), cenXYZ(i, indY),  &
          cenXYZ(i, indZ), ifUse, numConts, numVerts, xVerts, yVerts, contourZ, indVert)
      !
      ! Eliminate points that do not pass selection criteria
      if (numColSelect > 0 .and. ifUse > 0) then
        do ind = 1, numColSelect
          if ((extraVals(icolSelect(ind), i) - selectCrit(indSelect, ind)) * &
              isignSelect(ind) < 0) ifUse = 0
        enddo
      endif
      !
      if (ifUse > 0) then
        numFit = numFit + 1
        indToPatch(numFit) = i
        ! write(*,'(3f6.0)') cenXYZ(i, 1), cenXYZ(i, 3), cenXYZ(i, 2)
        do j = 1, 3
          !
          ! the regression requires coordinates of second volume as
          ! independent variables (columns 1-3), those in first volume
          ! as dependent variables (stored in 5-7), to obtain
          ! transformation to get from second to first volume
          ! cx+dx in second volume matches cx in first volume
          !
          fitMat(j + 4, numFit) = cenXYZ(i, j) - 0.5 * nxyz(j)
          fitMat(j, numFit) = fitMat(j + 4, numFit) + vecXYZ(i, j)
        enddo
      endif
    enddo

    if (numFit < 4) call exitError('Too few data points for fitting')
    print *,numFit, ' data points will be used for fit'
    !
    ! Get the solution
    maxDrop = nint(fracDrop * numFit)
    call solve_wo_outliers(fitMat, matCols, numFit, 3, icolFixed, maxDrop, critProb, &
        absProbCrit, elimMin, idrop, ndrop, aMat, delXYZ, fitCenter, devMean, devSD, &
        devMax, ipntMax, devXYZmax)
    !
    ! Leave the loop if pass the limit
    if (devMean <= stopLim) then
      exit
    endif
    if (indSelect < numSelectCrit) write(*,101) devMean, devMax
  enddo

  if (ndrop .ne. 0) then
    write(*,104) ndrop
104 format(/,i3,' patches dropped by outlier elimination:')
    if (ndrop <= 10) then
      print *,'    patch position     residual'
      do i = 1, ndrop
        write(*,106) (fitMat(j, numFit + i - ndrop), j = 1, 4)
106     format(3f7.0,f9.2)
      enddo
    else
      do i = 1, ndrop
        dropRes(i) = fitMat(4, numFit + i - ndrop)
      enddo
      call summarizeDrops(dropRes, ndrop, ' ')
    endif
  endif
  !
  write(*,101) devMean, devMax
101 format(/,' Mean residual',f8.3,',  maximum',f8.3,/)
  !
  print *,'Refining transformation:'
  write(*,102) ((aMat(i, j), j = 1, 3), delXYZ(i), i = 1, 3)
102 format(3f10.6,f10.3)
  !
  if (pipInput) then
    !
    ! For patch residual, first make a proper index from original to ordered rows
    do i = 1, numFit
      fitMat(6, nint(fitMat(5, i))) = i
    enddo
    if (PipGetString('ResidualPatchOutput', filename) == 0) then
      call dopen(1, filename, 'new', 'f')
      if (maxExtra > 0) then
        write(1, 120) numFit, ' positions', idResidCol, (IDextra(i), i = 1, maxExtra)
      else
        write(1, 120) numFit, ' positions', idResidCol
120     format(i7,a,20i8)
      endif
      do i = 1, numFit
        ipat = indToPatch(i)
        if (maxExtra > 0) then
          write(1, 110) (nint(cenXYZ(ipat, j)), j = 1, 3), (vecXYZ(ipat, j), j= 1, 3), &
              fitMat(4, nint(fitMat(6, i))), (extraVals(j, ipat), j = 1, maxExtra)
        else
          write(1, 110) (nint(cenXYZ(ipat, j)), j = 1, 3), (vecXYZ(ipat, j), j= 1, 3), &
              fitMat(4, nint(fitMat(6, i)))
110       format(3i6,3f9.2,f10.2,20f12.4)
        endif
      enddo
      close(1)
    endif
    !
    ! For reduced vectors, put out the residual vector stored in cols 15-17
    ! plus either the residual or original extra values
    if (PipGetString('ReducedVectorOutput', filename) == 0) then
      call dopen(1, filename, 'new', 'f')
      if (maxExtra > 0) then
        write(1, 120) numFit, ' positions', (IDextra(i), i = 1, maxExtra)
      else
        write(1, 120) numFit, ' positions', idResidCol
      endif
      do i = 1, numFit
        ind = nint(fitMat(6, i))
        ipat = indToPatch(i)
        if (maxExtra > 0) then
          write(1, 111) (nint(cenXYZ(ipat, j)), j = 1, 3), (fitMat(j, ind), j = 15, 17), &
              (extraVals(j, ipat), j = 1, maxExtra)
111       format(3i6,3f9.2,20f12.4)
        else
          write(1, 110) (nint(cenXYZ(ipat, j)), j = 1, 3), (fitMat(j, ind), j = 15, 17), &
              fitMat(4, ind)
        endif
      enddo
      close(1)
    endif
    filename = ' '
    ierr = PipGetString('OutputFile', filename)
  else
    print *,'Enter name of file to place transformation in, or Return for none'
    read(5, '(a)') filename
  endif
  if (filename .ne. ' ') then
    call dopen(1, filename, 'new', 'f')
    write(1, 102) ((aMat(i, j), j = 1, 3), delXYZ(i), i = 1, 3)
    close(1)
  endif

  ! Take product with initial file and output that if requested
  if (initialFile .ne. ' ') then
    call xfmult3d(aMatInit, delInit, aMat, delXYZ, prodMat, prodDel)
    call dopen(1, productFile, 'new', 'f')
    write(1, 102) ((prodMat(i, j), j = 1, 3), prodScaleFac * prodDel(i), i = 1, 3)
    close(1)
    !
    ! This triggers reporting of center shift as Y residual of used patch  nearest center
    cenDistMin = 1.e37
    do i = 1, numFit - ndrop
      ipat = nint(fitMat(5, i))
      cenDist = sqrt((cenXYZ(ipat, 1) - nxyz(1) / 2.)**2 + &
          (cenXYZ(ipat, 3) - nxyz(3) / 2.)**2)
      if (cenDist < cenDistMin) then
        cenDistMin = cenDist
        indClosest = i
      endif
    enddo
    write(*,'(a,f8.1)')'Implied center shift is ', fitMat(16, indClosest) * prodScaleFac
  endif
  
  if (devMean > stopLim) then
    write(*,'(/,a)') 'REFINEMATCH - Mean residual too high;'// &
        ' either raise the limit or use warping'
    call exit(2)
  endif
  call exit(0)
99 call exitError('Reading initial transform file')
end program refinematch
