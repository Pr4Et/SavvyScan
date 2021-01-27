! * * * * * FINDWARP * * * * * *
!
! FINDWARP will solve for a series of general 3-dimensional linear
! transformations that can then be used by WARPVOL to align two
! volumes to each other.  It performs a series of multiple linear
! regression on subsets of the displacements between the volumes
! determined at a matrix of positions (patches) .
!
! See man page for details
!
! $Id$
!
program findwarp
  implicit none
  integer LIMTARG, LIMEXTRA
  parameter (LIMTARG = 30, LIMEXTRA = 20)
  real*4 firstAmat(3,3), firstDelta(3), a(3,3), delXYZ(3), cenSaveMin(3), cenSaveMax(3)
  real*4 devXYZmax(3), cenLocal(3), amatTmp(3,3), delTmp(3), cenXYZsum(3)
  integer*4 nxyzVol(3)
  real*4 debugXYZ(3), dxLocal, dyLocal, dzLocal
  integer*4, allocatable :: indDropped(:), numTimesDropped(:), idrop(:)
  real*4, allocatable :: xVerts(:), yVerts(:), contourZ(:), dropSum(:), fitMat(:,:)
  integer*4, allocatable :: indVertStart(:), numVerts(:), listPositions(:,:)
  real*4, allocatable :: amatSave(:,:,:), cenToSave(:,:), cenXYZ(:,:), vecXYZ(:,:)
  real*4, allocatable :: delXYZsave(:,:), residSum(:), devMeanAuto(:,:), devMaxAuto(:,:)
  logical, allocatable :: solved(:), exists(:)
  integer*4, allocatable :: numXautoFit(:), numYautoFit(:), numZautoFit(:), inDiag2(:)
  integer*4, allocatable :: inRowX(:), inRowY(:), inRowZ(:), numResid(:), inDiag1(:)
  real*4, allocatable :: extraVals(:,:), autoMeanPatches(:)
  real*4 freNum(LIMEXTRA)
  character*320 filename, residFile, line
  real*4 cenXYZin(3), vecXYZin(3), targetResid(LIMTARG), selectCrit(LIMTARG, LIMEXTRA)
  integer*4 numSelectCrit, indSelect, icolSelect(LIMEXTRA), isignSelect(LIMEXTRA)
  integer*4 numXYZpatch(3), indXYZ(3), numZpatchUse, limPatch, limAxis, limDiag, limFit
  integer*4 numXfit, numYfit, numZfit, numXYXfit(3), numXpatchUse, numYpatchUse
  integer*4 numXfitIn, numYfitIn, numZfitIn, numXYZfitIn(3), numXYZpatchUse(3)
  integer*4 numXYZfit(3), numColSelect, numeric(LIMEXTRA), IDextra(LIMEXTRA)
  equivalence (numXpatchTot, numXYZpatch(1)), (numYpatchTot, numXYZpatch(2)), &
      (numZpatchTot, numXYZpatch(3))
  equivalence (numXfit, numXYXfit(1)), (numYfit, numXYXfit(2)), (numZfit, numXYXfit(3))
  equivalence (numXpatchUse, numXYZpatchUse(1)), (numYpatchUse, numXYZpatchUse(2)), &
      (numZpatchUse, numXYZpatchUse(3))
  equivalence (numXfitIn, numXYZfitIn(1)), (numYfitIn, numXYZfitIn(2)), &
      (numZfitIn, numXYZfitIn(3))
  equivalence (numXfit, numXYZfit(1)), (numYfit, numXYZfit(2)), (numZfit, numXYZfit(3))
  integer*4 numData, numPosInFile, numXpatchTot, numYpatchTot, numZpatchTot, i, j, ind
  integer*4 numConts, ierr, indY, indZ, numTarget, indTarget, intCenPos, k, itmp
  integer*4 numXoffset, numYoffset, numZoffset, ifSubset, ifLocalSlabs, ifdebug
  real*4 ratioMin, ratioMax, fracDrop, probCrit, absProbCrit, elimMinResid
  integer*4 ifAuto, numAuto, ix, iy, iz, indAuto, limVert, limCont, matColDim, maxExtra
  integer*4 numLocalDone, numXexclHigh, numYexclhigh, numZexclHigh, numFields, numExtraID
  integer*4 ifDiddle, numXlocal, numYlocal, numZlocal, numZero, numDevSum, minExtent
  integer*4 indUse, nlistDropped, numDropTot, locX, locY, locZ, lx, ly, lz, ifUse
  real*4 devMeanMin, devMeanSum, devMaxSum, devMaxMax, patchMeanNum
  real*4 dzMin, dz, devMax, devMean, devSD, discount, desiredMaxMax
  real*4 devMeanAvg, devMaxAvg, devMeanMax, determMean
  integer*4 icontMin, icont, indv, indLcl, ipntMax, maxDrop, numLowDeterm, maxAutoPatch
  integer*4 ifInDrop, numDrop, ifFlip, indPatch, indLocal, icolFixed, nyDiag
  integer*4 numElimSelect
  character*5 rowSlabText(2) /'rows ', 'slabs'/
  character*5 rowSlabCapText(2) /'ROWS ', 'SLABS'/
  character*1 yzText(2) /'Y', 'Z'/
  logical*4 debugHere, legacyRatios
  integer*4 numberInList, readNumPatches, imodGetEnv
  !
  logical pipInput
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetInteger, PipGetTwoIntegers
  integer*4 PipGetString, PipGetFloat, PipGetFloatArray, PipGetTwoFloats
  integer*4 PipGetInOutFile, PipGetThreeFloats, PipNumberOfEntries
  !
  ! fallbacks from ../../manpages/autodoc2man -3 2  findwarp
  !
  integer numOptions
  parameter (numOptions = 22)
  character*(40 * numOptions) options(1)

  indPatch(ix, iy, iz) = ix + (iy - 1) * numXpatchTot + (iz - 1) * numXpatchTot * &
      numYpatchTot
  indLocal(ix, iy, iz) = ix + (iy - 1) * numXlocal + (iz - 1) * numXlocal * numYlocal

  options(1) = &
      'patch:PatchFile:FN:@output:OutputFile:FN:@region:RegionModel:FN:@'// &
      'volume:VolumeOrSizeXYZ:FN:@initial:InitialTransformFile:FN:@'// &
      'residual:ResidualPatchOutput:FN:@target:TargetMeanResidual:FA:@'// &
      'measured:MeasuredRatioMinAndMax:FP:@xskip:XSkipLeftAndRight:IP:@'// &
      'yskip:YSkipLowerAndUpper:IP:@zskip:ZSkipLowerAndUpper:IP:@'// &
      'extra:ExtraValueSelection:IPM:@select:SelectionCriteria:FAM:@'// &
      'rowcol:LocalRowsAndColumns:IP:@slabs:LocalSlabs:I:@'// &
      'maxfrac:MaxFractionToDrop:F:@minresid:MinResidualToDrop:F:@'// &
      'prob:CriterionProbabilities:FP:@discount:DiscountIfZeroVectors:F:@'// &
      'debug:DebugAtXYZ:FT:@param:ParameterFile:PF:@help:usage:B:'
  !
  matColDim = 20
  numXoffset = 0
  numYoffset = 0
  numZoffset = 0
  numXexclHigh = 0
  numYexclhigh = 0
  numZexclHigh = 0
  ifSubset = 0
  ifDebug = 0
  ratioMin = 4.0
  ratioMax = 20.0
  fracDrop = 0.1
  probCrit = 0.01
  absProbCrit = 0.002
  elimMinResid = 0.5
  desiredMaxMax = 0.
  ifLocalSlabs = 0
  residFile = ' '
  discount = 0.
  minExtent = 2
  limFit = 40000
  numAuto = 0
  numSelectCrit = 0
  numColSelect = 0
  indSelect = 1
  !
  call PipReadOrParseOptions(options, numOptions, 'findwarp', &
      'ERROR: FINDWARP - ', .true., 3, 1, 1, numOptArg, numNonOptArg)
  pipInput = numOptArg + numNonOptArg > 0
  !
  ! Open patch file and figure out sampled positions
  !
  if (PipGetInOutFile('PatchFile', 1, &
      'Name of file with correlation positions and results', filename) .ne. 0) &
      call exitError('No input patch file specified')

  call dopen(1, filename, 'old', 'f')
  IDextra(:) = 0
  if (readNumPatches(1, numData, numExtraID, IDextra, LIMEXTRA) .ne. 0) call exitError &
      ('Reading first line of patch file')
  allocate(listPositions(numData + 4,3), stat = ierr)
  call memoryError(ierr, 'array for position list')
  legacyRatios = numExtraID == 0

  numXpatchTot = 0
  numYpatchTot = 0
  numZpatchTot = 0
  maxExtra = 0
  do i = 1, numData
    read(1,'(a)') line
    call frefor2(line, freNum, numeric, numFields, LIMEXTRA)
    maxExtra = max(maxExtra, numFields - 6)
    do j = 1, 3
      intCenPos = nint(freNum(j))
      if (numberInList(intCenPos, listPositions(1, j), numXYZpatch(j), 0) == 0) then
        numXYZpatch(j) = numXYZpatch(j) + 1
        listPositions(numXYZpatch(j), j) = intCenPos
      endif
    enddo
  enddo
  !
  ! sort the position lists
  !
  do i = 1, 3
    do j = 1, numXYZpatch(i) - 1
      do k = j, numXYZpatch(i)
        if (listPositions(j, i) > listPositions(k, i)) then
          itmp = listPositions(j, i)
          listPositions(j, i) = listPositions(k, i)
          listPositions(k, i) = itmp
        endif
      enddo
    enddo
  enddo
  !
  ! Allocate arrays based on full set of patch positions
  limPatch = numXpatchTot * numYpatchTot * numZpatchTot + 10
  limAxis = maxval(numXYZpatch) + 10
  limDiag = 2 * limAxis
  allocate(amatSave(3,3,limPatch), cenToSave(3,limPatch), delXYZsave(3,limPatch), &
      residSum(limPatch), solved(limPatch), exists(limPatch), cenXYZ(limPatch,3), &
      vecXYZ(limPatch,3), inRowX(limAxis), inRowY(limAxis), inRowZ(limAxis), &
      numResid(limPatch), inDiag1(0:limDiag), inDiag2(0:limDiag), indDropped(limPatch), &
      numTimesDropped(limPatch), dropSum(limPatch), stat = ierr)
  call memoryError(ierr, 'arrays for patches')
  if (maxExtra > 0) then
    allocate(extraVals(maxExtra, limPatch), stat=ierr)
    call memoryError(ierr, 'array for extra values')
    extraVals(:,:) = 0
  endif
  !
  numXpatchUse = numXpatchTot
  numYpatchUse = numYpatchTot
  numZpatchUse = numZpatchTot
  print *,'Number of patches in X, Y and Z is:', numXpatchTot, numYpatchTot, numZpatchTot
  rewind(1)
  read(1,*) numData
  !
  if (.not. pipInput) print *,'Enter NX, NY, NZ or name of file ', &
      'for tomogram being matched to'
  call get_nxyz(pipInput, 'VolumeOrSizeXYZ', 'FINDWARP', 1, nxyzVol)
  !
  ! mark positions as nonexistent and fill cx from list positions
  !
  do k = 1, numZpatchTot
    do j = 1, numYpatchTot
      do i = 1, numXpatchTot
        ind = indPatch(i, j, k)
        exists(ind) = .false.
        cenXYZ(ind, 1) = listPositions(i, 1)
        cenXYZ(ind, 2) = listPositions(j, 2)
        cenXYZ(ind, 3) = listPositions(k, 3)
      enddo
    enddo
  enddo
  !
  ! read each line, look up positions in list and store in right place
  !
  do i = 1, numData
    !
    ! these are center coordinates and location of the second volume
    ! relative to the first volume
    !
    read(1,'(a)') line
    call frefor2(line, freNum, numeric, numFields, LIMEXTRA)
    do j = 1, 3
      intCenPos = nint(freNum(j))
      k = 1
      do while(k <= numXYZpatch(j) .and. intCenPos > listPositions(k, j))
        k = k + 1
      enddo
      indXYZ(j) = k
    enddo
    ind = indPatch(indXYZ(1), indXYZ(2), indXYZ(3))
    exists(ind) = .true.
    do j = 1, 3
      cenXYZ(ind, j) = freNum(j)
      vecXYZ(ind, j) = freNum(j + 3)
    enddo
    do j = 7, numFields
      extraVals(j - 6, ind) = freNum(j)
    enddo
  enddo
  !
  numPosInFile = numData
  close(1)
  deallocate(listPositions, stat=ierr)
  !
  ! get patch region model, which should give flip value; otherwise set flip from patches
  !
  if (pipInput) then
    filename = ' '
    ierr = PipGetString('RegionModel', filename)
  else
    write(*,'(1x,a,/,a,$)') 'Enter name of model file with contour enclosing area to ' &
        //'use,', ' or Return to use all patches: '
    read(*,'(a)') filename
  endif
  numConts = 0
  if (filename .ne. ' ') then
    call getContourArraySizes(filename, 0, limCont, limVert)
    allocate(xVerts(limVert), yVerts(limVert), contourZ(limCont),  &
        indVertStart(limCont), numVerts(limCont), stat = ierr)
    call memoryError(ierr, 'arrays for boundary contours')
    call get_region_contours(filename, 'FINDWARP', xVerts, yVerts, numVerts, &
        indVertStart, contourZ, numConts, ifFlip, limCont, limVert, 0)
  else
    ifFlip = 0
    if (numYpatchTot < numZpatchTot) ifFlip = 1
  endif
  !
  ! Set indexes to Y-extent (row) and thickness (slab) variables
  indY = 2
  if (ifFlip .ne. 0) indY = 3
  indZ = 5 - indY
  numXYZfitIn(indZ) = numXYZpatchUse(indZ)
  numXYZfit(indZ) = numXYZpatchUse(indZ)
  !
  ! Get most other PIP entries
  if (pipInput) then
    ierr = PipGetFloat('MaxFractionToDrop', fracDrop)
    ierr = PipGetFloat('MinResidualToDrop', elimMinResid)
    ierr = PipGetTwoFloats('CriterionProbabilities', probCrit, absProbCrit)
    ierr = PipGetString('ResidualPatchOutput', residFile)
    ierr = PipGetFloat('DiscountIfZeroVectors', discount)
    ierr = PipGetInteger('MinExtentToFit', minExtent)
    minExtent = max(2, minExtent)
    if (imodGetEnv('FINDWARP_LEGACY_RATIOS', line) == 0) then
      read(line, *, iostat = ierr) ix
      if (ierr == 0 .and. ix .ne. 0) legacyRatios = ix > 0
    endif
    if (PipGetInteger('LegacyRatioEvaluation', ix) == 0) then
      if (ix .ne. 0) legacyRatios = ix > 0
    endif
    ifDebug = 1 - PipGetThreeFloats('DebugAtXYZ', debugXYZ(1), debugXYZ(2), &
        debugXYZ(3))
    !
    ! Get, translate and validate the extra selection entries
    call getExtraSelections(icolSelect, isignSelect, numColSelect, numSelectCrit, &
        IDextra, maxExtra, selectCrit, LIMEXTRA, LIMTARG)
  endif
  !
  ! aspectmax=3.
  !
  ! THIS POINT IS LOOPED BACK TO INTERACTIVELY TO CHANGE PARAMETERS OR GO TO AUTOFIT
  !
8 if (pipInput) then
    ifAuto = PipGetTwoIntegers('LocalRowsAndColumns', numXfitIn, numXYZfitIn(indY))
  else
    write(*,'(1x,a,$)') '1 to find best warping automatically, 0 '// &
        'to proceed interactively: '
    read(5,*) ifAuto
  endif
  !
  if (ifAuto .ne. 0) then
    !
    ! initialize auto 1: get target and ratio parameters
    !
    if (pipInput) then
      numTarget = 0
      if (PipGetFloatArray('TargetMeanResidual', targetResid, numTarget, LIMTARG) > 0) &
          call exitError('Target mean residual must be entered for automatic fits')
      desiredMaxMax = 8. * targetResid(numTarget)
      ierr = PipGetFloat('DesiredMaxResidual', desiredMaxMax)
      ierr = PipGetTwoFloats('MeasuredRatioMinAndMax', ratioMin, ratioMax)
    else
      write(*,'(1x,a,$)') 'One or more mean residuals to achieve: '
      read(5, '(a)') filename
      call frefor(filename, targetResid, numTarget)
      write(*,'(1x,a,/,a,2f5.1,a,$)') 'Minimum and maximum ratio of '// &
          'measurements to unknowns to test', '  (/ for' &
          , ratioMin, ratioMax, '): '
      read(5,*) ratioMin, ratioMax
    endif
    ! write(*,'(1x,a,f5.0,a,$)') &
    ! 'Maximum aspect ratio to allow in an area to be fit (/ for' &
    ! , aspectmax, '): '
    ! read(5,*) aspectmax
  endif
  !
  ! Get rows, columns, slabs to exclude and check entries
  !
  if (pipInput) then
    ierr = PipGetTwoIntegers('XSkipLeftAndRight', numXoffset, numXexclHigh) + &
        PipGetTwoIntegers('YSkipLowerAndUpper', numYoffset, numYexclhigh) + &
        PipGetTwoIntegers('ZSkipLowerAndUpper', numZoffset, numZexclHigh)
    ifSubset = 3 - ierr
  else
    write(*,'(1x,a,$)') '0 to include all positions, or 1 to '// &
        'exclude rows or columns of patches: '
    read(5,*) ifSubset
  endif
  !
  if (ifSubset .ne. 0) then
10  if (.not. pipInput) then
      numXexclHigh = numXpatchTot - numXoffset - numXpatchUse
      write(*,'(1x,a,2i3,a,$)') '# of columns to exclude' &
          //' on the left and right in X (/ for', numXoffset, numXexclHigh, '): '
      read(5,*) numXoffset, numXexclHigh
    endif
    numXpatchUse = numXpatchTot - numXoffset - numXexclHigh
    if (numXpatchUse + numXoffset > numXpatchTot .or. numXpatchUse < 2) then
      if (ifAuto .ne. 0) call exitError( &
          'Illegal entry for number of columns to exclude in X')
      print *,'Illegal entry'
      numXoffset = 0
      numXpatchUse = numXpatchTot
      go to 10
    endif

12  if (.not. pipInput) then
      numYexclHigh = numYpatchTot - numYoffset - numYpatchUse
      write(*,'(1x,a,a,a,2i3,a,$)') '# of ', rowSlabText(1 + ifFlip), ' to exclude' &
          //' on the bottom and top in Y (/ for', numYoffset, numYexclHigh, '): '
      read(5,*) numYoffset, numYexclHigh
    endif
    numYpatchUse = numYpatchTot - numYoffset - numYexclHigh
    if (numYpatchUse + numYoffset > numYpatchTot .or. numYpatchUse < 1 ) then
      if (ifAuto .ne. 0) call exitError('Illegal entry for number of '// &
          rowSlabCapText(1 + ifFlip) //' to exclude in Y')
      print *,'Illegal entry'
      numYoffset = 0
      numYpatchUse = numYpatchTot
      go to 12
    endif

14  if (.not. pipInput) then
      numZexclHigh = numZpatchTot - numZoffset - numZpatchUse
      write(*,'(1x,a,a,a,2i3,a,$)') '# of ', rowSlabText(2 - ifFlip), ' to exclude' &
          //' on the bottom and top in Z (/ for', numZoffset, numZexclHigh, '): '
      read(5,*) numZoffset, numZexclHigh
    endif
    numZpatchUse = numZpatchTot - numZoffset - numZexclHigh
    if (numZpatchUse + numZoffset > numZpatchTot .or. numZpatchUse < 2) then
      if (ifAuto .ne. 0) call exitError( 'Illegal entry for number of '// &
          rowSlabCapText(2 - ifFlip) //' to exclude in Z')
      print *,'Illegal entry'
      numZoffset = 0
      numZpatchUse = numZpatchTot
      go to 14
    endif
    print *,'Remaining # of patches in X, Y and Z is:', &
        numXpatchUse, numYpatchUse, numZpatchUse
  else
    numXpatchUse = numXpatchTot
    numYpatchUse = numYpatchTot
    numZpatchUse = numZpatchTot
    numXoffset = 0
    numYoffset = 0
    numZoffset = 0
  endif
  !
  ! Initialize nfit for slabs before checking on slabs
  if (pipInput .and. ifAuto .ne. 0) then
    numXYXfit(indZ) = numXYZpatchUse(indZ)
    numXYZfitIn(indZ) = numXYXfit(indZ)
  endif
  !
  ! Figure out if doing subsets of slabs but not for interactive when
  ! auto was already selected
  if (pipInput .or. ifAuto == 0) then
    if (numXYZpatchUse(indZ) > 2) then
      if (pipInput) then
        ifLocalSlabs = 1 - PipGetInteger('LocalSlabs', numXYXfit(indZ))
        if (ifLocalSlabs > 0) then
          if (numXYXfit(indZ) < 2 .or. numXYXfit(indZ) > numXYZpatch(indZ)) &
              call exitError('Number of local slabs out of allowed range')
          numXYZfitIn(indZ) = numXYXfit(indZ)
        endif
      else if (ifAuto == 0) then
        write(*,'(1x,a,a,a,a,a,$)') '0 to fit to all patches in ', yzText(2 - ifFlip), &
            ', or 1 to fit to subsets in ', yzText(2 - ifFlip), ': '
        read(5,*) ifLocalSlabs
      endif
    endif
  endif

  if (ifAuto .ne. 0) then
    !
    ! set up parameters for automatic finding: make list of possible
    ! nxfit, nyfit, nzfit values
    !
    indAuto = -1
    numXYZfit(indZ) = min(numXYZfit(indZ), numXYZpatchUse(indZ))
    call setAutoFits(indAuto)
    if (numAuto == 0) call exitError('No fitting parameters give the required ratio '// &
        'of measurements to unknowns - there are probably too few patches')
    allocate(numXautoFit(indAuto), numYautoFit(indAuto), numZautoFit(indAuto), &
        autoMeanPatches(indAuto), devMeanAuto(indAuto, max(1,numSelectCrit)), &
        devMaxAuto(indAuto, max(1,numSelectCrit)), stat = ierr)
    call memoryError(ierr, 'arrays for autofits')
    call setAutoFits(indAuto)
    ! write(*,'(4i5, i6, f7.1)') (i, numXautoFit(i), numYautoFit(i), numZautoFit(i), &
    ! numXautoFit(i) * numYautoFit(i) * numZautoFit(i), autoMeanPatches(i), i=1, numAuto)
    !
    ! set up for first round and skip to set up this round's patches
    !
    indAuto = 1
    numXfit = -100
    numYfit = -100
    numZfit = -100
    numLocalDone = 0
    devMeanMin = 10000.
    indTarget = 1
    write(*,112) targetResid(1)
112 format(/,'Seeking a warping with mean residual below',f9.3)
    if (numSelectCrit > 0) write(*,132)(selectCrit(1, i), i = 1, numColSelect)
132 format(2x,'with extra column selection criteria', 10f9.3)
    call countExtraEliminations()
    go to 20
  endif
  !
  ! Get parameter to control outlier elimination
  !
  if (.not.pipInput) then
    write(*,'(1x,a,$)') '1 to enter parameters to control outlier elimination, 0 not to: '
    read(5,*) ifDiddle
    if (ifDiddle .ne. 0) then
      write(*,'(1x,a,f5.2,a,$)') 'Maximum fraction of patches to eliminate (/ for', &
          fracDrop, '): '
      read(5,*) fracDrop
      write(*,'(1x,a,f5.2,a,$)') 'Minimum residual needed to do any elimination (/ for', &
          elimMinResid, '): '
      read(5,*) elimMinResid
      write(*,'(1x,a,f6.3,a,$)') 'Criterion probability for ' // &
          'candidates for elimination (/ for', probCrit, '): '
      read(5,*) probCrit
      write(*,'(1x,a,f6.3,a,$)') 'Criterion probability for enforced' &
          //' elimination (/ for', absProbCrit, '): '
      read(5,*) absProbCrit
    endif
  endif
  !
  ! Initialize for first time in or back after parameter setting
  !
  numXfit = -100
  numYfit = -100
  numZfit = -100
  numLocalDone = 0
  if (.not.pipInput) write(*,111)
111 format(/,' Enter 0 for the number of patches in X to loop', &
      ' back and find best fit ',/, &
      '  automatically, include a different subset of patches ', &
      /,'  or specify new outlier elimination parameters',/)
  !
  ! THIS POINT IS LOOPED BACK TO AFTER EVERY INTERACTIVE OR AUTO FIT
  !
  ! Get the input number of local patches unless doing auto
  !
20 if (ifAuto == 0 .and. ifLocalSlabs .ne. 0) then
    if (.not.pipInput) then
      if (numLocalDone == 0) then
        write(*,'(1x,a,$)') 'Number of local patches for fit in X, Y and Z: '
      else
        write(*,'(1x,a,/,a,$)') 'Number of local patches for fit in X, Y and Z,', &
            '    or / to redo and save last result: '
      endif
      read(5,*) numXfitIn, numYfitIn, numZfitIn
    endif
  elseif (ifAuto == 0) then
    if (.not.pipInput) then
      if (numLocalDone == 0) then
        write(*,'(1x,a,a,a,$)') &
            'Number of local patches for fit in X and ', yzText(1 + ifFlip), ': '
      else
        write(*,'(1x,a,a,a,/,a,$)') &
            'Number of local patches for fit in X and ', yzText(1 + ifFlip), ',', &
            '    or / to redo and save last result: '
      endif
      read(5,*) numXfitIn, numXYZfitIn(indY)
    endif
    numXYZfitIn(indZ) = numXYZpatchUse(indZ)
  else
    numXfitIn = numXautoFit(indAuto)
    numYfitIn = numYautoFit(indAuto)
    numZfitIn = numZautoFit(indAuto)
    patchMeanNum = autoMeanPatches(indAuto)
  endif
  !
  ! Loop back to earlier parameter entries on a zero entry
  if (numXfitIn == 0) go to 8
  !
  ! Save data and terminate on a duplicate entry
  if (numXfitIn == numXfit .and. numYfitIn == numYfit .and. numZfitIn == numZfit .and. &
      numLocalDone > 0) &
      call saveAndTerminate()
  !
  ! Otherwise set up number of locations to fit and check entries
  numXlocal = numXpatchUse + 1 - numXfitIn
  numYlocal = numYpatchUse + 1 - numYfitIn
  numZlocal = numZpatchUse + 1 - numZfitIn
  if (numXfitIn < 2 .or. numXYZfitIn(indY) < 2 .or. numXlocal < 1 .or. numZlocal < 1 &
      .or. numXYZfitIn(indZ) < 1 .or. numYlocal < 1) then
    if (ifAuto .ne. 0) call exitError('Improper number to include in fit')
    print *,'Illegal entry, try again'
    go to 20
  endif
  !
  ! If arrays not allocated yet, do them to big size for interactive or current size
  ! for one-shot run from PIP or biggest size from auto setup
  if (.not. allocated(fitMat)) then
    if (pipInput) then
      limFit = numXfitIn * numYfitIn * numZfitIn + 10
      if (ifAuto .ne. 0) limFit = maxAutoPatch + 10
    endif
    allocate(fitMat(matColDim, limFit), idrop(limFit), stat = ierr)
    call memoryError(ierr, 'arrays for fitting')
  endif
  if (numXfitIn * numYfitIn * numZfitIn > limFit) then
    if (ifAuto .ne. 0) call exitError('Too many patches for array sizes')
    print *,'Too many patches for array sizes, try again'
    go to 20
  endif
  numXfit = numXfitIn
  numYfit = numYfitIn
  numZfit = numZfitIn
  if (ifAuto == 0) patchMeanNum = numXfit * numYfit * numZfit
  !
  ! Do the fits finally
  !
  call fitLocalPatches()
  !
  ! check for auto control
  !
  if (ifAuto .ne. 0) then
    devMeanAvg = 10000.
    if (numDevSum > 0) devMeanAvg = devMeanSum / numDevSum
    devMeanMin = min(devMeanMin, devMeanAvg)
    devMeanAuto(indAuto, indSelect) = devMeanAvg
    devMaxAuto(indAuto, indSelect) = devMaxMax
    if (devMeanAvg <= targetResid(indTarget) .and. (indSelect >= numSelectCrit .or.  &
        desiredMaxMax <= 0. .or. devMaxMax <= desiredMaxMax)) then
      if (numSelectCrit > 0 .and. indTarget == 1) &
          write(*,'(4x,a,i6,a,i6,a)')'- eliminates', numElimSelect,' of', numPosInFile, &
          ' patches'
      !
      ! done: set nauto to zero to allow printing of results
      !
      numAuto = 0
      write(*,107) numXfit, numYfit, numZfit
107   format(/,'Desired residual achieved with fits to',i3,',',i3, &
          ', and',i3,' patches in X, Y, Z',/)
    else
      indAuto = indAuto + 1
      ! print *,iauto, ' Did fit to', nfitx, nfity, nfitz
      if (indAuto > numAuto) then
        if (numSelectCrit > 0 .and. indTarget == 1) &
            write(*,'(4x,a,i6,a,i6,a,f8.3)')'- eliminates', numElimSelect,' of',  &
            numPosInFile, ' patches, gives best mean residual', devMeanMin
        !
        ! If there are multiple selection criteria, first redo all fits with the next
        if (indSelect < numSelectCrit) then
          indAuto = 1
          numXfit = -100
          numLocalDone = 0
          indSelect = indSelect + 1
          write(*,132)(selectCrit(indSelect, i), i = 1, numColSelect)
          call countExtraEliminations()
          go to 20
        endif
        !
        ! See if any other criteria have been met and loop back if so
        !
        do j = indTarget + 1, numTarget
          write(*,112) targetResid(j)
          lx = 1
          ly = max(1, numSelectCrit)
          do ix = lx, ly
            if (numSelectCrit > 0) write(*,132)(selectCrit(ix, i), i = 1, numColSelect)
            do i = 1, numAuto
              if (devMeanAuto(i, ix) <= targetResid(j) .and. (ix == ly .or. &
                  desiredMaxMax <= 0. .or. devMaxAuto(i, ix) <= desiredMaxMax)) then
                indTarget = j
                indSelect = ix
                indAuto = i
                numXfit = -100
                numLocalDone = 0
                go to 20
              endif
            enddo
          enddo
        enddo
        !
        ! write patch file if desired on last fit before error message
        !
        call outputPatchRes(residFile, numPosInFile, numXpatchTot * numYpatchTot * &
            numZpatchTot, exists, residSum, numResid, indDropped, numTimesDropped, &
            nlistDropped, cenXYZ, vecXYZ, limPatch, maxExtra, extraVals, IDextra)
        if (discount > 0. .and. numLocalDone > 0 .and. numDevSum == 0) &
            call exitError('All fits had too many zero vectors: raise '// &
            '-discount fraction or set it to zero')
        write(*,108) devMeanMin
108     format(/,'ERROR: FINDWARP - Failed to find a warping with a ', &
            'mean residual below',f9.3)
        call exit(2)
      endif
    endif
  endif

  if (nlistDropped > 0 .and. numAuto == 0) then
    write(*,105) nlistDropped, numDropTot
105 format(i5,' separate patches eliminated as outliers a total' &
        ,' of',i6,' times:')
    if (nlistDropped <= 10) then
      print *,'    patch position   # of times   mean residual'
      do i = 1, nlistDropped
        write(*,106) (cenXYZ(indDropped(i), j), j = 1, 3), numTimesDropped(i), &
            dropSum(i) / numTimesDropped(i)
106     format(3f7.0,i8,f12.2)
      enddo
    else
      do i = 1, nlistDropped
        dropSum(i) = dropSum(i) / numTimesDropped(i)
      enddo
      call summarizeDrops(dropSum, nlistDropped, 'mean ')
    endif
    write(*,*)
  endif
  !
  if (numLocalDone > 0 .and. numAuto == 0) then
    devMeanAvg = devMeanSum / max(1, numDevSum)
    devMaxAvg = devMaxSum / max(1, numDevSum)
    write(*,101) devMeanAvg, devMeanMax, devMaxAvg, devMaxMax
101 format('Mean residual has an average of',f8.3, &
        ' and a maximum of',f8.3,/ , &
        'Max  residual has an average of',f8.3,' and a maximum of' &
        ,f8.3)
    if (discount > 0.) write(*,'(/,a,i7,a,i7,a)') 'These averages are '// &
        'based on', numDevSum, ' of', numLocalDone, ' fits'
    !
    ! finish up after auto fits: this call is unneeded and is just here to make it clear
    if (ifAuto .ne. 0) call saveAndTerminate()
  elseif (ifAuto == 0) then
    print *,'No locations could be solved for'
  endif
  go to 20


CONTAINS

  ! fitLocalPatches does the fits to all the local patches 
  !
  subroutine fitLocalPatches()
    logical inside
    real*4 determ
    devMeanSum = 0.
    devMaxSum = 0.
    devMaxMax = 0.
    devMeanMax = 0.
    nlistDropped = 0
    numDropTot = 0
    numLocalDone = 0
    numDevSum = 0.
    determMean = 0.
    cenSaveMin(:) = 1.e30
    cenSaveMax(:) = -1.e30
    numResid(:) = 0
    residSum(:) = 0.

    do locZ = 1, numZlocal
      do locY = 1, numYlocal
        do locX = 1, numXlocal
          numData = 0
          numZero = 0
          do i = 1, 3
            cenXYZsum(i) = 0.
          enddo
          !
          ! count up number in each row in each dimension and on each
          ! diagonal in the major dimensions
          !
          do i = 1, max(numXfit, numYfit, numZfit)
            inRowX(i) = 0
            inRowY(i) = 0
            inRowZ(i) = 0
          enddo
          nyDiag = numYfit
          if (ifFlip == 1) nyDiag = numZfit
          !
          ! Zero the parts of the arrays corresponding to corners, which
          ! won't be tested
          do i = 0, numXfit + nyDiag - 2
            inDiag1(i) = 0
            inDiag2(i) = 0
          enddo
          do lz = locZ + numZoffset, locZ + numZoffset + numZfit - 1
            do ly = locY + numYoffset, locY + numYoffset + numYfit - 1
              do lx = locX + numXoffset, locX + numXoffset + numXfit - 1
                ifUse = 0
                ind = indPatch(lx, ly, lz)
                if (exists(ind)) ifUse = 1
                do i = 1, 3
                  cenXYZsum(i) = cenXYZsum(i) + cenXYZ(ind, i) - 0.5 * nxyzVol(i)
                enddo
                if (numConts > 0 .and. ifUse > 0) call checkBoundaryConts(cenXYZ(ind, 1),&
                    cenXYZ(ind, indY), cenXYZ(ind, indZ), ifUse, numConts, numVerts, &
                    xVerts, yVerts, contourZ, indVertStart)
                !
                ! Eliminate points that do not pass selection criteria
                if (numColSelect > 0 .and. ifUse > 0) then
                  do i = 1, numColSelect
                    if ((extraVals(icolSelect(i), ind) - selectCrit(indSelect, i)) * &
                        isignSelect(i) < 0) ifUse = 0
                  enddo
                endif
                !
                if (ifUse > 0) then
                  numData = numData + 1
                  if (vecXYZ(ind, 1) == 0. .and. vecXYZ(ind, 2) == 0 .and. &
                      vecXYZ(ind, 3) == 0) numZero = numZero + 1
                  do j = 1, 3
                    !
                    ! the regression requires coordinates of second volume as
                    ! independent variables (columns 1-3), those in first
                    ! volume as dependent variables (stored in 5-7), to
                    ! obtain transformation to get from second to first
                    ! volume cx+dx in second volume matches cx in first
                    ! volume
                    !
                    fitMat(j + 4, numData) = cenXYZ(ind, j) - 0.5 * nxyzVol(j)
                    fitMat(j, numData) = fitMat(j + 4, numData) + vecXYZ(ind, j)
                  enddo
                  !
                  ! Solve_wo_outliers uses columns 8-17; save indexi in 18
                  ! Add to row counts and to diagonal counts
                  !
                  fitMat(18, numData) = ind
                  ix = lx + 1 - locX - numXoffset
                  iy = ly + 1 - locY - numYoffset
                  iz = lz + 1 - locZ - numZoffset
                  inRowX(ix) = inRowX(ix) + 1
                  inRowY(iy) = inRowY(iy) + 1
                  inRowZ(iz) = inRowZ(iz) + 1
                  if (ifFlip == 1) iy = iz
                  iz = (ix - iy) + nyDiag - 1
                  inDiag1(iz) = inDiag1(iz) + 1
                  iz = ix + iy - 2
                  inDiag2(iz) = inDiag2(iz) + 1
                endif
              enddo
            enddo
          enddo
          indLcl = indLocal(locX, locY, locZ)
          !
          ! Need regular array of positions, so use the xyzsum to get
          ! censave, not the cenloc values from the regression
          !
          debugHere = ifDebug .ne. 0
          do i = 1, 3
            cenToSave(i, indLcl) = cenXYZsum(i) / (numXfit * numYfit * numZfit)
            cenSaveMin(i) = min(cenSaveMin(i), cenToSave(i, indLcl))
            cenSaveMax(i) = max(cenSaveMax(i), cenToSave(i, indLcl))
            if (debugHere) &
                debugHere = abs(cenToSave(i, indLcl) - debugXYZ(i)) < 1.
          enddo
          !
          ! solve for this location if there are at least half of the
          ! normal number of patches present and if there are guaranteed
          ! to be at least 3 patches in a different row from the dominant
          ! one, even if the max are dropped from other rows
          ! But treat thickness differently: if there are not enough data
          ! on another layer, or if there is only one layer being fit,
          ! then set the appropriate column as fixed in the fits
          !
          solved(indLcl) = numData >= patchMeanNum / 2
          maxDrop = nint(fracDrop * numData)
          icolFixed = 0
          do i = 1, max(numXfit, numYfit, numZfit)
            if (debugHere) print *,'in row', i, ':', inRowX(i), inRowY(i), inRowZ(i)
            if (ifFlip == 1) then
              if (inRowX(i) > numData - 3 - maxDrop .or. &
                  inRowZ(i) > numData - 3 - maxDrop) solved(indLcl) = .false.
              if (inRowY(i) > numData - 3 - maxDrop .or. numYfit == 1) icolFixed = 2
            else
              if (inRowX(i) > numData - 3 - maxDrop .or. &
                  inRowY(i) > numData - 3 - maxDrop) solved(indLcl) = .false.
              if (inRowZ(i) > numData - 3 - maxDrop .or. numZfit == 1) icolFixed = 3
            endif
          enddo
          do i = 1, numXfit + nyDiag - 3
            if (inDiag1(i) > numData - 3 - maxDrop .or. &
                inDiag2(i) > numData - 3 - maxDrop) solved(indLcl) = .false.
          enddo
          if (solved(indLcl)) then
            call solve_wo_outliers(fitMat, matColDim, numData, 3, icolFixed, maxDrop, &
                probCrit, absProbCrit, elimMinResid, idrop, numDrop, a, delXYZ, &
                cenLocal, devMean, devSD, devMax, ipntMax, devXYZmax)
            !
            if (debugHere) then
              do i = 1, numData
                write(*,'(8f9.2)') (fitMat(j, i), j = 1, 7)
              enddo
              print *,'cenloc', (cenLocal(i), i = 1, 3)
              print *,'censave', (cenToSave(i, indLcl), i = 1, 3)
              write(*,'(9f8.3)') ((a(i, j), i = 1, 3), j = 1, 3)
            endif
            !
            ! Accumulate information about dropped points
            !
            do i = 1, numDrop
              ifInDrop = 0
              do j = 1, nlistDropped
                if (nint(fitMat(18, idrop(i))) == indDropped(j)) then
                  ifInDrop = 1
                  numTimesDropped(j) = numTimesDropped(j) + 1
                  dropSum(j) = dropSum(j) + fitMat(4, numData + i - numDrop)
                endif
              enddo
              if (ifInDrop == 0 .and. ifInDrop < limFit) then
                nlistDropped = nlistDropped + 1
                indDropped(nlistDropped) = nint(fitMat(18, idrop(i)))
                numTimesDropped(nlistDropped) = 1
                dropSum(nlistDropped) = fitMat(4, numData + i - numDrop)
              endif
            enddo
            numDropTot = numDropTot + numDrop
            !
            ! if residual output asked for, accumulate info about all resids
            !
            if (residFile .ne. ' ') then
              do i = 1, numData
                ind = nint(fitMat(18, nint(fitMat(5, i))))
                numResid(ind) = numResid(ind) + 1
                residSum(ind) = residSum(ind) + fitMat(4, i)
              enddo
            endif
            !
            if (discount == 0. .or. float(numZero) / numData <= discount) &
                then
              devMeanSum = devMeanSum + devMean
              devMaxSum = devMaxSum + devMax
              numDevSum = numDevSum + 1
            endif
            devMeanMax = max(devMeanMax, devMean)
            devMaxMax = max(devMaxMax, devMax)
            !
            ! mark this location as solved and save the solution.
            !
            if (debugHere) write(*,'(6i4,3f8.1)') indLcl, locX, locY, locZ, numData, &
                numDrop, (delXYZ(i), i = 1, 3)
            do i = 1, 3
              delXYZsave(i, indLcl) = delXYZ(i)
              do j = 1, 3
                amatSave(i, j, indLcl) = a(i, j)
              enddo
              if (debugHere) write(*,122) (a(i, j), j = 1, 3), delXYZ(i)
122             format(3f10.6,f10.3)
            enddo
            numLocalDone = numLocalDone + 1
            determMean = determMean + abs(determ(a))
          elseif (debugHere) then
            print *,'Not solved', numData
          endif
        enddo
      enddo
    enddo
    determMean = determMean / max(1, numLocalDone)
    return
  end subroutine fitLocalPatches

  !
  ! Count up eliminations by the selection criteria
  !
  subroutine countExtraEliminations()
    if (numColSelect > 0) then
      numElimSelect = 0
      do ind = 1, numXpatchTot * numYpatchTot * numZpatchTot
        if (exists(ind)) then
          ifUse = 1
          do i = 1, numColSelect
            if ((extraVals(icolSelect(i), ind) - selectCrit(indSelect, i)) * &
                isignSelect(i) < 0) ifUse = 0
          enddo
          if (ifUse == 0) numElimSelect = numElimSelect + 1
        endif
      enddo
    endif
    return
  end subroutine countExtraEliminations


  ! Sets up a list of number of local patches for autofits, where each one has a ratio of 
  ! measured to unknown within the min and max range.  Call with limAuto <= 0 to just
  ! count up the number and return a good value for limAuto
  !
  subroutine setAutoFits(limAuto)
    implicit none
    integer*4 minInXYZ(3), nxFit, nyFit, nzFit, nXYZusable(3), minLocXYZ(3), maxLocXYZ(3)
    integer*4 limAuto, numExists, maxInXYZ(3)
    real*4 ratio, ratioFac, ratioAvg, patchMean, rtmp, avgRelax
    avgRelax = 0.8
    numAuto = 0
    maxAutoPatch = 0
    maxInXYZ(:) = numXYZpatchUse(:)
    minInXYZ(:) = min(numXYZpatchUse(:), minExtent)
    minInXYZ(indZ) = numXYZpatchUse(indZ)
    if (ifLocalSlabs .ne. 0) minInXYZ(indZ) = numXYZfit(indZ)
    ratioFac = 4.0
    if (numXYZpatchUse(indZ) == 1) ratioFac = 3.0
    maxInXYZ(1) = nint(5. * ratioMax * ratioFac / (minInXYZ(2) * minInXYZ(3)))
    maxInXYZ(2) = nint(5. * ratioMax * ratioFac / (minInXYZ(1) * minInXYZ(3)))
    maxInXYZ(3) = nint(5. * ratioMax * ratioFac / (minInXYZ(1) * minInXYZ(2)))
    maxInXYZ(:) = min(numXYZpatchUse(:), max(maxInXYZ(:), minInXYZ(:)))
    do nxFit = maxInXYZ(1), minInXYZ(1), -1
      do nyFit = maxInXYZ(2), minInXYZ(2), -1
        do nzFit = maxInXYZ(3), minInXYZ(3), -1
          numXlocal = numXpatchUse + 1 - nxFit
          numYlocal = numYpatchUse + 1 - nyFit
          numZlocal = numZpatchUse + 1 - nzFit
          !
          ! For legacy ratios, simply take the nominal # to fit as usable number
          if (legacyRatios) then
            nXYZusable(1) = nxFit
            nXYZusable(2) = nyFit
            nXYZusable(3) = nzFit
            numExists = nxFit * nyFit * nzFit * numXlocal * numYlocal * numZlocal
            avgRelax = 1.
          else
            !
            ! Otherwise, loop on all local patches for this size and determine maximum
            ! existing size in each dimension
            nXYZusable(:) = 0
            numExists = 0
            do locZ = 1, numZlocal
              do locY = 1, numYlocal
                do locX = 1, numXlocal
                  maxLocXYZ(:) = 0
                  minLocXYZ(:) = numXYZpatch(:)
                  do lz = locZ + numZoffset, locZ + numZoffset + nzFit - 1
                    do ly = locY + numYoffset, locY + numYoffset + nyFit - 1
                      do lx = locX + numXoffset, locX + numXoffset + nxFit - 1
                        if (exists(indPatch(lx, ly, lz))) then
                          numExists = numExists + 1
                          minLocXYZ(1) = min(minLocXYZ(1), lx)
                          maxLocXYZ(1) = max(maxLocXYZ(1), lx)
                          minLocXYZ(2) = min(minLocXYZ(2), ly)
                          maxLocXYZ(2) = max(maxLocXYZ(2), ly)
                          minLocXYZ(3) = min(minLocXYZ(3), lz)
                          maxLocXYZ(3) = max(maxLocXYZ(3), lz)
                        endif
                      enddo
                    enddo
                  enddo
                  nXYZusable(:) = max(nXYZusable(:), maxLocXYZ(:) + 1 - minLocXYZ(:))
                enddo
              enddo
            enddo
          endif
          !
          ! Get the mean number of patches in each area and a ratio based on that,
          ! plus a ratio based on the nominal patch number.  The fit is acceptable if
          ! the nominal ratio is within limits and the average ratio is almost within
          ! the limits; also allow larger fits where the average ratio is well within
          ! the upper limit
          patchMean = float(numExists) / (numXlocal * numYlocal * numZlocal) 
          ratioAvg = patchMean / ratioFac
          ratio = nXYZusable(1) * nXYZusable(2) * nXYZusable(3) / ratioFac
          if ((nxFit < numXpatchUse .or. (ifFlip == 0 .and. nyFit < numYpatchUse) .or. &
              (ifFlip .ne. 0 .and. nzFit < numZpatchUse)) .and.  &
              ratio >= ratioMin .and. ratioAvg >= avgRelax * ratioMin .and.  &
              (ratio <= ratioMax .or. ratioAvg <= avgRelax * ratioMax .or.  &
              (nxFit == minInXYZ(1) .and. nyFit == minInXYZ(2) .and.  &
              nzFit == minInXYZ(3)))) then
            numAuto = numAuto + 1
            if (limAuto > 0) then
              numXautoFit(numAuto) = nxFit
              numYautoFit(numAuto) = nyFit
              numZautoFit(numAuto) = nzFit
              autoMeanPatches(numAuto) = patchMean
              maxAutoPatch = max(maxAutoPatch, nxFit * nyFit * nzFit)
            endif
          endif
        enddo
      enddo
    enddo
    if (limAuto <= 0) then
      limAuto = numAuto + 10
      return
    endif
    !
    ! sort the list by size of area in inverted order
    do i = 1, numAuto - 1
      do j = i + 1, numAuto
        if (autoMeanPatches(i) < autoMeanPatches(j)) then
        !if (numXautoFit(i) * numYautoFit(i) * numZautoFit(i) < &
        !    numXautoFit(j) * numYautoFit(j) * numZautoFit(j)) then
          itmp = numXautoFit(i)
          numXautoFit(i) = numXautoFit(j)
          numXautoFit(j) = itmp
          itmp = numYautoFit(i)
          numYautoFit(i) = numYautoFit(j)
          numYautoFit(j) = itmp
          itmp = numZautoFit(i)
          numZautoFit(i) = numZautoFit(j)
          numZautoFit(j) = itmp
          rtmp = autoMeanPatches(i)
          autoMeanPatches(i) = autoMeanPatches(j)
          autoMeanPatches(j) = rtmp
        endif
      enddo
    enddo
    return
  end subroutine setAutoFits


  ! saveAndTerminate saves the results (getting initial file and output file) and exits
  !
  subroutine saveAndTerminate()
    real*4 determ
    if (numXlocal > 1 .or. numYlocal > 1 .or. numZlocal > 1) then
      !
      ! Eliminate locations with low determinants.  This is very
      ! conservative measure before observed failures were diagonal
      ! degeneracies in 2x2 fits
      numLowDeterm = 0
      do locZ = 1, numZlocal
        do locY = 1, numYlocal
          do locX = 1, numXlocal
            ind = indLocal(locX, locY, locZ)
            if (solved(ind)) then
              if (abs(determ(amatSave(1, 1, ind))) < 0.01 * determMean) then
                solved(ind) = .false.
                numLowDeterm = numLowDeterm + 1
              endif
            endif
          enddo
        enddo
      enddo
      !
      if (numLowDeterm > 0) write(*,'(/,i4,a)') numLowDeterm, &
          ' fits were eliminated due to low matrix determinant'
      !
      if (ifAuto .ne. 0) write(*,*)
      if (pipInput) then
        ierr = PipGetString('InitialTransformFile', filename)
      else
        print *,'Enter name of file with initial transformation, typically solve.xf', &
            '   (Return if none)'
        read(*,'(a)') filename
      endif
      if (filename .ne. ' ') then
        call dopen(1, filename, 'old', 'f')
        read(1,*) ((firstAmat(i, j), j = 1, 3), firstDelta(i), i = 1, 3)
        close(1)
      else
        do i = 1, 3
          do j = 1, 3
            firstAmat(i, j) = 0.
          enddo
          firstAmat(i, i) = 1.
          firstDelta(i) = 0.
        enddo
      endif
      !
      if (PipGetInOutFile('OutputFile', 2, &
          'Name of file to place warping transformations in', filename) == 0) then
        call dopen(1, filename, 'new', 'f')
        !
        ! Output new style header to allow missing data
        dxLocal = 1.
        dyLocal = 1.
        dzLocal = 1.
        if (numXlocal > 1) dxLocal = (cenSaveMax(1) - cenSaveMin(1)) / (numXlocal - 1)
        if (numYlocal > 1) dyLocal = (cenSaveMax(2) - cenSaveMin(2)) / (numYlocal - 1)
        if (numZlocal > 1) dzLocal = (cenSaveMax(3) - cenSaveMin(3)) / (numZlocal - 1)
        write(1, 104) numXlocal, numYlocal, numZlocal, (cenSaveMin(i), i = 1, 3), &
            dxLocal, dyLocal, dzLocal
104     format(i5,2i6,3f11.2,3f10.4)
        !
        do locZ = 1, numZlocal
          do locY = 1, numYlocal
            do locX = 1, numXlocal
              ind = indLocal(locX, locY, locZ)
              indUse = ind
              !
              ! if this location was solved, combine and invert
              if (solved(ind)) then
                call xfmult3d(firstAmat, firstDelta, amatSave(1, 1, indUse), &
                    delXYZsave(1, indUse), amatTmp, delTmp)
                call xfinv3d(amatTmp, delTmp, a, delXYZ)
                write(1, 103) (cenToSave(i, ind), i = 1, 3)
103             format(3f9.1)
                write(1, 102) ((a(i, j), j = 1, 3), delXYZ(i), i = 1, 3)
102             format(3f10.6,f10.3)
              endif
            enddo
          enddo
        enddo
        close(1)
      endif
    else
      !
      ! Save a single transform if fit to whole area; won't happen if auto
      !
      print *,'Enter name of file in which to place single refining transformation'
      read(5, '(a)') filename
      call dopen(1, filename, 'new', 'f')
      write(1, 102) ((a(i, j), j = 1, 3), delXYZ(i), i = 1, 3)
      close(1)
    endif

    call outputPatchRes(residFile, numPosInFile, numXpatchTot * numYpatchTot * &
        numZpatchTot, exists, residSum, numResid, indDropped, numTimesDropped, &
        nlistDropped, cenXYZ, vecXYZ, limPatch, maxExtra, extraVals, IDextra)

    call exit(0)
  end subroutine saveAndTerminate

end program findwarp


!
! Output new patch file with mean residuals if requested
!
subroutine outputPatchRes(residFile, numPosInFile, numPatchTot, exists, residSum, &
    numResid, indDropped, numTimesDropped, nlistDropped, cenXYZ, vecXYZ, limPatch, &
    maxExtra, extraVals, IDextra)
  implicit none
  character*(*) residFile
  integer*4 numPosInFile, numPatchTot, limPatch, numResid(*), ind, i, indDropped(*)
  integer*4 numTimesDropped(*), nlistDropped, maxExtra, IDextra(*)
  real*4 cenXYZ(limPatch,3), vecXYZ(limPatch,3), residSum(*), dist, dropFrac
  real*4 extraVals(maxExtra, limPatch)
  logical exists(*)
  integer*4 idResidCol/3/, idOutlierCol/4/
  if (residFile == ' ') return
  call dopen(1, residFile, 'new', 'f')
  if (maxExtra > 0) then
    write(1, '(i7,a,20i8)') numPosInFile, ' positions', idResidCol, idOutlierCol,  &
        (IDextra(i), i = 1, maxExtra)
  else
    write(1, '(i7,a,20i8)') numPosInFile, ' positions', idResidCol, idOutlierCol
  endif
  do ind = 1, numPatchTot
    if (exists(ind)) then
      dropFrac = 0.
      do i = 1, nlistDropped
        if (ind == indDropped(i)) then
          dropFrac = float(numTimesDropped(i)) / max(numTimesDropped(i), 1, numResid(ind))
          exit
        endif
      enddo
      dist = residSum(ind) / max(1, numResid(ind))
      if (maxExtra > 0) then
        write(1, 110) (nint(cenXYZ(ind, i)), i = 1, 3), (vecXYZ(ind, i), i = 1, 3), &
            dist, dropFrac, (extraVals(i, ind), i = 1, maxExtra)
      else
        write(1, 110) (nint(cenXYZ(ind, i)), i = 1, 3), (vecXYZ(ind, i), i = 1, 3), &
            dist, dropFrac
      endif
110   format(3i6,3f9.2,f10.2,f7.3,20f12.4)
    endif
  enddo
  close(1)
  return
end subroutine outputPatchRes

real*4 function determ(a)
  implicit none
  real * 4 a(3, 3)
  determ = a(1, 1) * a(2, 2) * a(3, 3) + a(2, 1) * a(3, 2) * a(1, 3) +  &
      a(1, 2) * a(2, 3) * a(3, 1) - a(1, 3) * a(2, 2) * a(3, 1) -  &
      a(2, 1) * a(1, 2) * a(3, 3) - a(1, 1) * a(3, 2) * a(2, 3)
  return
end function determ
