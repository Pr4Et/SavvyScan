! * * * * * CORRSEARCH3D * * * * *
!
! CORRSEARCH3D will determine the 3D displacement between two image
! volumes at a regular array of positions.  At each position, it
! extracts a patch of a specified size from each volume, then
! searches for a local maximum in the cross-correlation between the
! two patches.  The starting point of the search is based upon the
! displacements of adjacent patches that have already been determined.
! If there are no such adjacent patches, or if no maximum is found for
! displacements within a specified range, the program uses FFT-based
! cross-correlation instead.
!
! See man page for details.
!
! David Mastronarde, 7/16/01
!
! $Id$
!
program corrsearch3d
  implicit none
  integer LIM_LOCAL_BIN
  parameter (LIM_LOCAL_BIN = 6)
  integer*4 nx, ny, nz, nx2, ny2, nz2, nxyz(3), nxyz2(3)
  real*4 dxInitial, dyInitial, dzInitial, dxyzInit(3)
  real*4 dxVolume, dyVolume, dzVolume, dxyzVol(3)
  equivalence (dxInitial, dxyzInit(1)), (dyInitial, dxyzInit(2)), (dzInitial, dxyzInit(3))
  equivalence (dxVolume, dxyzVol(1)), (dyVolume, dxyzVol(2)), (dzVolume, dxyzVol(3))
  equivalence (nx, nxyz(1)), (ny, nxyz(2)), (nz, nxyz(3))
  equivalence (nx2, nxyz2(1)), (ny2, nxyz2(2)), (nz2, nxyz2(3))
  integer*4 nxPatch, nyPatch, nzPatch, nxyzPatch(3)
  equivalence (nxPatch, nxyzPatch(1)), (nyPatch, nxyzPatch(2)), (nzPatch, nxyzPatch(3))
  !
  integer*4 mxyz(3), nxyzBsource(3), limVert, limCont
  real*4 ctf(8193)
  !
  logical inside, exist, readSmallMod, found, flipMessages, loadFullWidth
  integer getImodHead, getImodScales
  real*4, allocatable :: xVertex(:), yVertex(:), zcont(:)
  integer*4, allocatable :: indVert(:), numVertex(:), indVertB(:), numVertB(:)
  real*4, allocatable :: xVertexB(:), yVertexB(:), zcontB(:)
  character*320 fileA, fileB, outputFile, modelFile, tempFile, xfFile, bModel
  !
  ! limPatch
  real*4, allocatable :: dxPatch(:), dyPatch(:), dzPatch(:), directErr(:)
  integer*4, allocatable :: ifDone(:), ixSequence(:), iySequence(:)
  integer*4, allocatable :: izSequence(:), idirSequence(:)
  real*4, allocatable :: dxDirect(:), dyDirect(:), dzDirect(:), ccDirect(:), corrCoef(:)
  real*4, allocatable :: work(:), buffer(:)
  real*4 xvertBsource(4), yvertBsource(4), xvertSort(4), yvertSort(4)
  real*4 xvertTilt(4), yvertTilt(4), xvertTiltB(4), yvertTiltB(4)
  !
  integer*4 indPatch, ix, iy, iz, maxShift, limPatch, limWork
  integer*4 ifDebug, numFourierPatch, numCorrs, mode, numXpatch, numYpatch, numZpatch
  integer*4 nBordXlow, nbordXhigh, nbordYlow, nbordYhigh, nbordZlow, nbordZhigh, nxTaper
  integer*4 ixStart, iyStart, izStart, ixDelta, iyDelta, izDelta, nyTaper, nzTaper
  integer*4 nbordSourceXlow, nbordSourceXhigh, nbordSourceZlow, nbordSourceZhigh, indV
  real*4 dmin, dmax, dmean, dum, a11, a12, a21, a22, dxSource, dzSource, xx, zz
  real*4 tmp, radius2, sigma1, sigma2, sigmaKernel, delta, aScale, bScale
  integer*4 i, j, indYb, numSequence, nxCTF, nyCTF, nzCTF, nxPatchTmp, nyPatchTmp
  integer*4 nxPad, nyPad, nzPad, numPatchPixels, indPatchA, nzPatchTmp, idim
  integer*4 indPatchB, indLoadA, indLoadB, maxXload, numCont, idimOptimal, idimHigher
  integer*4 ierr, ifFlip, numPosTotal, numBcont, ifFlipB
  integer*4 nxSeries, nySeries, nxSeriesB, nySeriesB, ifAxis
  real*4 axisRotation
  real*4 xcen, ycen, zcen, xpatchLow, xpatchHigh, zpatchLow, zpatchHigh, dz, dzMin
  integer*4 indP, ifUse, icont, icontMin, ixSpan, izSpan, kernelSize, ixPatchEnd
  integer*4 izPatchStart, izPatchEnd, izDir, izPatch, iz0, iz1, izCen, izAdjacent, indA
  integer*4 iyPatchStart, iyPatchEnd, iyDir, iyPatch, iy0, iy1, iyCen, ixPatchStart
  integer*4 ixDir, ixPatch, ix0, ix1, ixCen, numAdjacent, ixAdjacent, iyAdjacent
  integer*4 loadY0, loadY1, loadZ0, loadZ1, loadX0, loadX1, numMore, nyLoadEx, nzLoadEx
  integer*4 nxLoad, nyLoad, nzLoad, indSequence, numNear, loadExtra, numSkip
  integer*4 loadYb0, loadYb1, loadZb0, loadZb1, loadXb0, loadXb1, indScratch
  integer*4 nxLoadB, nyLoadb, nzLoadB, ixB0, ixB1, iyB0, iyB1, izB0, izB1
  integer*4 ixMin, ixMax, iyMin, iyMax, izMin, izMax, ifShiftIn, numAdjacentLook
  integer*4 nxXCpad, nyXCpad, nzXCpad, nxXCbord, nyXCbord, nzXCbord, niceLim
  integer*4 num3CorrThreads, numSmoothThreads, numCCCthreads, ifHistVerbose
  real*4 aSource(3,3), dxyzSource(3), dxNew, dyNew, dzNew, peak, wsumAdjacent, wsumNear
  real*4 dxSum, dySum, dzSum, err, perPos, dxAdjacent, dyAdjacent, dzAdjacent, zmodCen
  real*4 dxSumNear, dySumNear, dzSumNear, dxNear, dyNear, dzNear, distSq, distNear
  real*4 ccRatio, wcc, sizeSwitch, ymodCen, distAdjacent, badHistPeakCrit
  real*4, allocatable :: statMeans(:), statSDs(:), statBuffer(:), patchFracHighSD(:)
  real*4, allocatable :: patchMeanStruct(:)
  integer*4 lsdBoxSize(3, LIM_LOCAL_BIN), numLsdBinnings, ibinBest, numInSample
  integer*4 lsdBinning(3, LIM_LOCAL_BIN) /1,1,1, 2,2,2, 3,3,3, 4,4,4, 6,6,6, 8,8,8/
  integer*4 lsdStartCoord(3), lsdEndCoord(3), lsdBoxStart(3, LIM_LOCAL_BIN)
  integer*4 numLsdBoxes(3, LIM_LOCAL_BIN), lsdBufferStarts(LIM_LOCAL_BIN + 1)
  integer*4 lsdStatStarts(LIM_LOCAL_BIN + 1), lsdSpacing(3, LIM_LOCAL_BIN)
  integer*4 ibStart(3), ibEnd(3), numHistBins, numElimBySD, loopElim
  real*4 edgeSampleFrac, fallbackStructFrac, betterScaleCrit, sdElimMaxFrac
  real*4 fallbackPctFrac, fracForHist, elimBySdCrit, elimBySdType, highSDlevel
  real*4 bottomHistFrac
  integer*4 idCCCcol/1/, idFracCol/5/, idStructCol/6/
  integer*4 niceFrame, usingFFTW, niceFFTlimit, numOMPthreads
  real*8 wallTime, wallStart, wallbest, wallCCC, wallCum(5)

  logical pipInput
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetInteger, PipGetTwoIntegers, PipGetThreeIntegers
  integer*4 PipGetString, PipGetFloat, PipGetThreeFloats, PipgetTwoFloats
  integer*4 PipGetInOutFile, PipGetLogical, PipGetBoolean
  !
  ! fallbacks from ../../manpages/autodoc2man -3 2  corrsearch3d
  !
  integer numOptions
  parameter (numOptions = 33)
  character*(40 * numOptions) options(1)

  indPatch(ix, iy, iz) = ix + (iy - 1) * numXpatch + (iz - 1) * numXpatch * numYpatch

  options(1) = &
      'ref:ReferenceFile:FN:@align:FileToAlign:FN:@output:OutputFile:FN:@'// &
      'region:RegionModel:FN:@size:PatchSizeXYZ:IT:@number:NumberOfPatchesXYZ:IT:@'// &
      'xminmax:XMinAndMax:IP:@yminmax:YMinAndMax:IP:@zminmax:ZMinAndMax:IP:@'// &
      'tilt:TiltSeriesSizeXY:IP:@btilt:BTiltSeriesSizeXY:IP:@'// &
      'axis:AxisRotationAngle:F:@taper:TapersInXYZ:IT:@pad:PadsInXYZ:IT:@'// &
      'maxshift:MaximumShift:I:@volume:VolumeShiftXYZ:FT:@'// &
      'initial:InitialShiftXYZ:FT:@bsource:BSourceOrSizeXYZ:CH:@'// &
      'bxform:BSourceTransform:FN:@bxborder:BSourceBorderXLoHi:IP:@'// &
      'byzborder:BSourceBorderYZLoHi:IP:@bregion:BRegionModel:FN:@'// &
      'binnings:LocalSDNumBinnings:I:@box:BoxSizeForLocalSD:IT:@'// &
      'elim:EliminateByLocalSD:FP:@kernel:KernelSigma:F:@ksize:KernelSize:F:@'// &
      'lowpass:LowPassRadiusSigma:FP:@sigma1:HighPassSigma:F:@'// &
      'messages:FlipYZMessages:B:@debug:DebugMode:I:@param:ParameterFile:PF:@'// &
      'help:usage:B:'
  !
  ! IFDEBUG 1 for debugging output, 2 for dummy patch output at all
  ! positions, 3 to do both search and fourier cross-correlation everywhere
  !
  ifDebug = 0
  numFourierPatch = 0
  numCorrs = 0
  nxyzBsource(1) = 0
  nbordSourceXlow = 36
  nbordSourceXhigh = 36
  nbordSourceZlow = 36
  nbordSourceZhigh = 36
  maxShift = 10
  do i = 1, 3
    dxyzVol(i) = 0.
    dxyzInit(i) = 0.
  enddo
  modelFile = ' '
  xfFile = ' '
  bModel = ' '
  radius2 = 0.
  sigma1 = 0.
  sigma2 = 0.
  sigmaKernel = 0.
  kernelSize = 3
  sizeSwitch = 1.49
  delta = 0.
  ifShiftIn = 0
  distAdjacent = 1.5
  distNear = 4.
  ccRatio = 0.33
  flipMessages = .false.
  numLsdBinnings = 0
  wallBest = 0.
  wallCum(:) = 0.
  wallCCC = 0.
  idimOptimal = 250000000
  idimHigher = 490000000
  numHistBins = 1000
  edgeSampleFrac = 0.5     ! Fraction of edge pixels to use for median/MADN
  fracForHist = 0.98       ! Top fraction of boxes to use for histogram (trim above)
  bottomHistFrac = 0.005   ! Bottom fraction of boxes for start of histogram
  fallbackPctFrac = 0.75   ! Percentile for fallback "typical high SD value"
  fallbackStructFrac = 0.5 ! Fraction of "typical high" to use like histo dip
  betterScaleCrit = 0.33   ! Difference of ratios must rise by this much to use
  elimBySdType = 0.
  elimBySdCrit = 0.5
  nxSeries = 0
  nxSeriesB = 0
  badHistPeakCrit = 5.
  sdElimMaxFrac = 0.95
  !
  ! Pip startup: set error, parse options, check help, set flag if used
  !
  call PipReadOrParseOptions(options, numOptions, 'corrsearch3d', &
      'ERROR: CORRSEARCH3D - ', .true., 3, 2, 1, numOptArg, &
      numNonOptArg)
  pipInput = numOptArg + numNonOptArg > 0
  !
  ! Open image files.
  !
  if (PipGetInOutFile('ReferenceFile', 1, 'Image file to align to', fileA) &
      .ne. 0) call exitError('No reference file specified')
  if (PipGetInOutFile('FileToAlign', 2, 'Image file being aligned', fileB) &
      .ne. 0) call exitError('No file to align specified')
  if (PipGetInOutFile('OutputFile', 3, 'Output file for displacements', &
      outputFile) .ne. 0) &
      call exitError('No output file for displacements specified')

  if (pipInput) then
    ierr = PipGetString('RegionModel', modelFile)
    ierr = PipGetString('BRegionModel', bModel)
  else
    write(*,'(1x,a,$)') 'Enter blank line : '
    read(5, 50) tempFile
    print *,'Enter model file with contours enclosing areas to ', &
        'analyze, or Return for none'
    read(5, 50) modelFile
50  format(a)
  endif
  !
  call imopen(1, fileA, 'RO')
  call irdhdr(1, nxyz, mxyz, mode, dmin, dmax, dmean)
  aScale = 1.
  if (dmax - dmin > 0.) aScale = 1.e10
  if (dmax - dmin > 1.e-10) aScale = 100. / (dmax - dmin)
  !
  call imopen(2, fileB, 'RO')
  call irdhdr(2, nxyz2, mxyz, mode, dmin, dmax, dmean)
  bScale = 1.
  if (dmax - dmin > 0.) bScale = 1.e10
  if (dmax - dmin > 1.e-10) bScale = 100. / (dmax - dmin)
  ! print *,'scaling factors', ascale, bscale
  !
  if (pipInput) then
    if (PipGetThreeIntegers('PatchSizeXYZ', nxPatch, nyPatch, nzPatch) .ne. &
        0) call exitError('No patch size specified')
    if (PipGetThreeIntegers('NumberOfPatchesXYZ', numXpatch, numYpatch, numZpatch) &
        .ne. 0) call exitError('No number of patches specified')
    if (PipGetTwoIntegers('XMinAndMax', ixMin, ixMax) + &
        PipGetTwoIntegers('YMinAndMax', iyMin, iyMax) + &
        PipGetTwoIntegers('ZMinAndMax', izMin, izMax) .ne. 0) call &
        exitError('Min and max coordinates in X, Y, and Z must be entered')
    j = 0
    ierr = PipGetBoolean('InvertYLimits', j)
    if (j > 0) then
      i = ny + 1 - iyMin
      iyMin = ny + 1 - iyMax
      iyMax = i
    endif
    nBordXlow = ixMin - 1
    nbordXhigh = nx - ixMax
    nbordYlow = iyMin - 1
    nbordYhigh = ny - iyMax
    nbordZlow = izMin - 1
    nbordZhigh = nz - izMax
    nxTaper = (nxPatch + 9) / 10
    nyTaper = (nyPatch + 9) / 10
    nzTaper = (nzPatch + 9) / 10
    !
    ! If maximum shift is not entered, scale it by square root of largest size over 1000.
    ! so it gets bigger for bigger sets
    ierr = PipGetThreeIntegers('TapersInXYZ', nxTaper, nyTaper, nzTaper)
    if (PipGetInteger('MaximumShift', maxShift) .ne. 0)  &
        maxShift = maxShift * max(1., sqrt(max(nx2, ny2, nz2) / 1000.))
    call get_nxyz(.true., 'BSourceOrSizeXYZ', ' ', 5, nxyzBsource)
    ierr = PipGetTwoIntegers('BSourceBorderXLoHi', nbordSourceXlow, nbordSourceXhigh)
    ierr = PipGetTwoIntegers('BSourceBorderYZLoHi', nbordSourceZlow, nbordSourceZhigh)
    ierr = PipGetString('BSourceTransform', xfFile)
    ifShiftIn = 1 - PipGetThreeFloats('VolumeShiftXYZ', dxVolume, &
        dyVolume, dzVolume)
    ierr = PipGetThreeFloats('InitialShiftXYZ', dxInitial, &
        dyInitial, dzInitial)
    ierr = PipGetTwoFloats('LowPassRadiusSigma', radius2, sigma2)
    ierr = PipGetFloat('HighPassSigma', sigma1)
    ierr = PipGetFloat('KernelSigma', sigmaKernel)
    if (sigmaKernel > sizeSwitch) kernelSize = 5
    ierr = PipGetInteger('KernelSize', kernelSize)
    if (kernelSize .ne. 3 .and. kernelSize .ne. 5) call exitError( &
        'Kernel size must be 3 or 5')
    ierr = PipGetLogical('FlipYZMessages', flipMessages)
    ierr = PipGetInteger('LocalSDNumBinnings', numLsdBinnings)
    if (numLsdBinnings > LIM_LOCAL_BIN) &
        call exitError('Number of local binnings is too high')
    if (numLsdBinnings > 0) then
      if (PipGetThreeIntegers('BoxSizeForLocalSD', lsdBoxSize(1, 1), lsdBoxSize(2, 1), &
          lsdBoxSize(3, 1)) .ne. 0) &
          call exitError('Box size must be entered for local SD analysis')
      ierr = PipGetTwoFloats('EliminateByLocalSD', elimBySdType, elimBySdCrit)
      if (2 * lsdBoxSize(1, 1) > nxPatch .or. 2 * lsdBoxSize(2, 1) > nyPatch .or.  &
          2 * lsdBoxSize(3, 1) > nzPatch) then
        do i = 1, 3
          lsdBoxSize(i, 1) = min(lsdBoxSize(i, 1), nxyzPatch(i) / 2)
        enddo
        write(*,'(a,a,3i4)')'WARNING: Box size for SD analysis too large for patches,', &
            ' changed to',(lsdBoxSize(i, 1), i = 1, 3)
      endif
    endif

    ifAxis = 1 - PipGetFloat('AxisRotationAngle', axisRotation)
    ierr = 1 - PipGetTwoIntegers('TiltSeriesSizeXY', nxSeries, nySeries)
    ix = 1 - PipGetTwoIntegers('BTiltSeriesSizeXY', nxSeriesB, nySeriesB)
    if (ifAxis == 0 .and. ix + ierr > 0) &
        call exitError('Axis rotation angle must be entered for tilt series size'// &
        ' to be useful')
    if (nxyzBsource(1) == 0) nxSeriesB = 0
    ierr = PipGetInteger('DebugMode', ifDebug)
  else
    write(*,'(1x,a,$)') 'X, Y, and Z size of patches: '
    read(5,*) nxPatch, nyPatch, nzPatch
    write(*,'(1x,a,$)') 'Number of patches in X, Y and Z: '
    read(5,*) numXpatch, numYpatch, numZpatch
    write(*,'(1x,a,$)') 'Border sizes on lower and upper sides in X,' &
        //' Y, and Z: '
    read(5,*) nBordXlow, nbordXhigh, nbordYlow, nbordYhigh, nbordZlow, nbordZhigh
    write(*,'(1x,a,$)') 'Number of pixels to taper in X, Y, Z: '
    read(5,*) nxTaper, nyTaper, nzTaper
    write(*,'(1x,a,$)') 'Maximum shift to analyze for by searching: '
    read(5,*) maxShift
    !
    ! get inputs for analyzing existence of data in B file
    !
    print *,'Enter name or NX, NY, NZ of the untransformed source', &
        ' for the image file being aligned, or Return for none'
    call get_nxyz(.false., ' ', ' ', 5, nxyzBsource)
    print *,'Enter name of file used to transform the image file', &
        ' being aligned, or Return for none'
    read(5, 50) xfFile
    write(*,'(1x,a,/,a,$)') 'Border sizes on lower and upper sides' &
        //' in X and in Y or Z in the untransformed', &
        ' source for the image file being aligned: '
    read(5,*) nbordSourceXlow, nbordSourceXhigh, nbordSourceZlow, nbordSourceZhigh
  endif
  ifHistVerbose = 0;
  do while (ifDebug >= 10)
    ifHistVerbose = ifHistVerbose + 1
    ifDebug = ifDebug - 10
  enddo
  !
  ! If there is no initial shift entered, enforce centered transform
  ! by setting initial shift to half the difference in size
  !
  if (ifShiftIn == 0) then
    dxVolume = (nxyz2(1) - nxyz(1)) / 2.
    dyVolume = (nxyz2(2) - nxyz(2)) / 2.
    dzVolume = (nxyz2(3) - nxyz(3)) / 2.
  endif
  !
  ! get padding for cross correlations
  !
  nxXCbord = (nxPatch + 4) / 5
  nyXCbord = (nyPatch + 4) / 5
  nzXCbord = (nzPatch + 4) / 5
  if (pipInput) &
      ierr = PipGetThreeIntegers('PadsInXYZ', nxXCbord, nyXCbord, nzXCbord)
  niceLim = niceFFTlimit()
  nxXCpad = niceFrame(nxPatch + 2 * nxXCbord, 2, niceLim)
  nyXCpad = niceFrame(nyPatch + 2 * nyXCbord, 2, niceLim)
  nzXCpad = niceFrame(nzPatch + 2 * nzXCbord, 2, niceLim)
  limWork = 10
  if (usingFFTW() == 0) limWork = (nxXCpad + 2) * nzXCpad + 10
  nxCTF = nxXCpad
  !
  ! Initialize transform: since this is a center-based transform,
  ! neutralize component of volume coordinate shift due to size difference
  do i = 1, 3
    do j = 1, 3
      aSource(i, j) = 0.
    enddo
    dxyzSource(i) = (nxyz2(i) - nxyz(i)) / 2. - dxyzVol(i)
    aSource(i, i) = 1.
  enddo
  !
  ! require b source dimensions if one is entered or
  ! fall back to current B volume size; read in transform and add shifts
  if (nxyzBsource(1) == 0) then
    if (xfFile .ne. ' ') call exitError( &
        'B source file dimensions must be entered to use 3D transforms')
    do i = 1, 3
      nxyzBsource(i) = nxyz2(i)
    enddo
  endif
  !
  ! Read transform, add volume-shift component of center-based transform
  if (xfFile .ne. ' ') then
    call dopen(1, xfFile, 'ro', 'f')
    do i = 1, 3
      read(1,*,err = 95, end = 95) (aSource(i, j), j = 1, 3), tmp
      dxyzSource(i) = dxyzSource(i) + tmp
    enddo
    close(1)
  endif
  !
  ! get model contours as primary source of flip information and fallback
  ! to Y/Z dimensions
  !
  numCont = 0
  if (modelFile .ne. ' ') then
    call getContourArraySizes(modelFile, 1, limCont, limVert)
    allocate(xVertex(limVert), yVertex(limVert), zcont(limCont), indVert(limCont), &
        numVertex(limCont), stat = ierr)
    call memoryError(ierr, 'arrays for contours')
    call get_region_contours(modelFile, 'CORRSEARCH3D', xVertex, yVertex, &
        numVertex, indVert, zcont, numCont, ifFlip, limCont, limVert, 1)


    if (ifDebug > 0) write(*,'(5(f7.0,f8.0))') &
        (xVertex(i), yVertex(i), i = indVert(1), indVert(1) + numVertex(1) - 1)
  else
    ifFlip = 0
    if (nz > ny) ifFlip = 1
  endif
  !
  ! now that flip is known, get patch starting positions and delta between
  ! patches
  !
  indYb = 2
  if (ifFlip .ne. 0 .and. flipMessages) indYb = 3
  call checkAndSetPatches(nx, nBordXlow, nbordXhigh, nxPatch, numXpatch, &
      ixStart, ixDelta, 1)
  call checkAndSetPatches(ny, nbordYlow, nbordYhigh, nyPatch, numYpatch, &
      iyStart, iyDelta, indYb)
  call checkAndSetPatches(nz, nbordZlow, nbordZhigh, nzPatch, numZpatch, &
      izStart, izDelta, 5 - indYb)


  numBcont = 0
  if (bModel .ne. ' ') then
    !
    ! For b model, get its coordinates in the source native section plane
    ! then transform to coordinates in A volume native plane
    !
    print *,'Processing model on source for second volume'
    call getContourArraySizes(bModel, 2, limCont, limVert)
    allocate(xVertexB(limVert), yVertexB(limVert), zcontB(limCont), indVertB(limCont), &
        numVertB(limCont), stat = ierr)
    call memoryError(ierr, 'arrays for contours')
    call get_region_contours(bModel, 'CORRSEARCH3D', xVertexB, yVertexB, &
        numVertB, indVertB, zcontB, numBcont, ifFlipB, limCont, limVert, 2)
    do j = 1, numBcont
      if (ifDebug > 0) write(*,'(5(f7.0,f8.0))') &
          (xVertexB(i), yVertexB(i), i = indVertB(j), indVertB(j) + numVertB(j) - 1)
      do i = indVertB(j), indVertB(j) + numVertB(j) - 1
        call xformBsourceToA(xVertexB(i), yVertexB(i), &
            nxyzBsource, nxyz, ifFlipB, ifFlip, aSource, dxyzSource, &
            xVertexB(i), yVertexB(i))
      enddo
      if (ifDebug > 0) write(*,'(5(f7.0,f8.0))') &
          (xVertexB(i), yVertexB(i), i = indVertB(j), indVertB(j) + numVertB(j) - 1)
    enddo
  else
    ifFlipB = 0
    if (nxyzBsource(3) > nxyzBsource(2)) ifFlipB = 1
  endif
  indYb = 2
  if (ifFlipB .ne. 0) indYb = 3
  if (ifDebug > 0) print *,'flips', ifFlip, ifFlipB
  !
  ! compute transformed locations of corners of b source volume
  ! 10/14/09: This used to be done inside the if below, but it is needed
  ! for all evaluations of the patches.
  call xformBsourceToA(float(nbordSourceXlow), float(nbordSourceZlow), &
      nxyzBsource, nxyz, ifFlipB, ifFlip, aSource, dxyzSource, &
      xvertBsource(1), yvertBsource(1))
  call xformBsourceToA(float(nxyzBsource(1) - nbordSourceXhigh), float(nbordSourceZlow), &
      nxyzBsource, nxyz, ifFlipB, ifFlip, aSource, dxyzSource, &
      xvertBsource(2), yvertBsource(2))
  call xformBsourceToA(float(nxyzBsource(1) - nbordSourceXhigh), &
      float(nxyzBsource(indYb) - nbordSourceZhigh), &
      nxyzBsource, nxyz, ifFlipB, ifFlip, aSource, dxyzSource, &
      xvertBsource(3), yvertBsource(3))
  call xformBsourceToA(float(nbordSourceXlow), float(nxyzBsource(indYb) -  &
      nbordSourceZhigh), nxyzBsource, nxyz, ifFlipB, ifFlip, aSource, dxyzSource, &
      xvertBsource(4), yvertBsource(4))
  if (ifDebug > 0) &
      print *,'bverts', (xvertBsource(i), yvertBsource(i), i = 1, 4)
  if (ifShiftIn .ne. 0 .or. xfFile .ne. ' ') then
    !
    ! now order the coordinates to find middle values that can be
    ! used to adjust the entered lower and upper limits, so that the
    ! basic grid can be set up to span between those limits
    !
    do i = 1, 4
      xvertSort(i) = xvertBsource(i)
      yvertSort(i) = yvertBsource(i)
    enddo
    do i = 1, 3
      do j = i + 1, 4
        if (xvertSort(i) > xvertSort(j)) then
          tmp = xvertSort(i)
          xvertSort(i) = xvertSort(j)
          xvertSort(j) = tmp
        endif
        if (yvertSort(i) > yvertSort(j)) then
          tmp = yvertSort(i)
          yvertSort(i) = yvertSort(j)
          yvertSort(j) = tmp
        endif
      enddo
    enddo
    !
    ! compute revised lower and upper limits for X and Y/Z
    !
    call revisePatchRange(nx, nBordXlow, nbordXhigh, xvertSort(2), xvertSort(3), &
        nxPatch, numXpatch, ixStart, ixDelta)
    if (ifFlip .ne. 0) then
      call revisePatchRange(nz, nbordZlow, nbordZhigh, yvertSort(2), yvertSort(3), &
          nzPatch, numZpatch, izStart, izDelta)
    else
      call revisePatchRange(ny, nbordYlow, nbordYhigh, yvertSort(2), yvertSort(3), &
          nyPatch, numYpatch, iyStart, iyDelta)
    endif
  endif
  limPatch = numXpatch * numYpatch * numZpatch + 10
  allocate(dxPatch(limPatch), dyPatch(limPatch), dzPatch(limPatch), directErr(limPatch), &
      ifDone(limPatch), ixSequence(limPatch), iySequence(limPatch), work(limWork), &
      izSequence(limPatch), idirSequence(limPatch), dxDirect(limPatch), &
      dyDirect(limPatch), dzDirect(limPatch), ccDirect(limPatch), corrCoef(limPatch), &
      stat = ierr)
  call memoryError(ierr, 'arrays for patch variables')
  !
  if (ifDebug .ne. 0) then
    print *,'Scan limits in X:', ixStart, ixStart + (numXpatch - 1) * ixDelta &
        + nxPatch
    print *,'Scan limits in Y:', iyStart, iyStart + (numYpatch - 1) * iyDelta &
        + nyPatch
    print *,'Scan limits in Z:', izStart, izStart + (numZpatch - 1) * izDelta &
        + nzPatch
  endif
  !
  ! Set up another set of boundaries if the tilt series size was entered
  ierr = 2
  if (ifFlip > 0) ierr = 3
  call tiltVertexSetup(nxSeries, nySeries, nxyz(1), nxyz(ierr), xvertTilt, yvertTilt)
  call tiltVertexSetup(nxSeriesB, nySeriesB, nxyzBsource(1), nxyzBsource(indYb),  &
      xvertTiltB, yvertTiltB)
  if (nxSeriesb > 0) then
    do i = 1, 4
      call xformBsourceToA(xvertTiltB(i), yvertTiltB(i), nxyzBsource, nxyz, ifFlipB, &
          ifFlip, aSource, dxyzSource, xvertTiltB(i), yvertTiltB(i))
    enddo
    if (ifDebug > 1) write(*,'(a, 8f9.1)')' B to A', (xvertTiltB(ix), yvertTiltB(ix), &
        ix = 1, 4)
  endif
  !
  ! Do analysis of local SDs
  if (numLsdBinnings > 0) then
    call analyzeLocalSDs()
  endif
  !
  ! prescan for patches inside boundaries, to get total count and flags
  ! for whether to do. Loop twice if doing SD analysis so that it can redo it without if
  ! too many get eliminated
  !
  do loopElim = 1, max(1, min(2, numLsdBinnings + 1))
    numPosTotal = 0
    numElimBySD = 0
    do iz = 1, numZpatch
      zcen = izStart + (iz - 1) * izDelta + nzPatch / 2
      do iy = 1, numYpatch
        ycen = iyStart + (iy - 1) * iyDelta + nyPatch / 2
        ymodCen = ycen
        zmodCen = zcen
        if (ifFlip .ne. 0) then
          zmodCen = ycen
          ymodCen = zcen
        endif
        !
        do ix = 1, numXpatch
          indP = indPatch(ix, iy, iz)
          xcen = ixStart + (ix - 1) * ixDelta + nxPatch / 2
          ifUse = 1
          !
          ! If SD analysis was done, eliminate first based on that
          if (numLsdBinnings > 0 .and. nint(elimBySdType) > 0) then
            if ((nint(elimBySdType) > 1 .and. &
                patchMeanStruct(indP) < elimBySdCrit * highSDlevel) .or. &
                (nint(elimBySdType) == 1 .and. patchFracHighSD(indP) < elimBySdCrit)) then
              ifUse = 0
              numElimBySD = numElimBySD + 1
              if (ifDebug == 2) write(*,'(3f8.1,a,f7.4)')xcen, ycen, zcen, &
                  ' eliminated by SD criterion, frac', patchFracHighSD(indP)
            endif
          endif
          !
          ! If A model was entered, find nearest contour in Z and see if
          ! patch is inside it
          if (numCont > 0 .and. ifUse == 1) then
            call checkBoundaryConts(xcen, ymodCen, zmodCen, ifUse, numCont, &
                numVertex, xVertex, yVertex, zcont, indVert)
            if (ifUse == 0 .and. ifDebug == 2) &
                print *,xcen, ymodCen, ' eliminated by A model contour', icontMin
          endif
          !
          ! Do the same if still ok and B model was entered
          if (numBcont > 0 .and. ifUse == 1) then
            call checkBoundaryConts(xcen, ymodCen, zmodCen, ifUse, numBcont, &
                numVertB, xVertexB, yVertexB, zcontB, indVertB)
            if (ifUse == 0 .and. ifDebug == 2) &
                print *,xcen, ymodCen, ' eliminated by B model'
          endif
          !
          if (ifUse == 1) then
            ifUse = 0
            !
            ! now make sure all corners of the patch are inside transformed
            ! area from B borders
            !
            xpatchLow = xcen - (nxPatch - nxTaper) / 2
            xpatchHigh = xcen + (nxPatch - nxTaper) / 2
            if (ifFlip .ne. 0) then
              zpatchLow = zcen - (nzPatch - nzTaper) / 2
              zpatchHigh = zcen + (nzPatch - nzTaper) / 2
            else
              zpatchLow = ycen - (nyPatch - nyTaper) / 2
              zpatchHigh = ycen + (nyPatch - nyTaper) / 2
            endif
            if (inside(xvertBsource, yvertBsource, 4, xpatchLow, zpatchLow) .and. &
                inside(xvertBsource, yvertBsource, 4, xpatchLow, zpatchHigh) .and. &
                inside(xvertBsource, yvertBsource, 4, xpatchHigh, zpatchLow) .and. &
                inside(xvertBsource, yvertBsource, 4, xpatchHigh, zpatchHigh)) ifUse = 1
            if (ifUse == 0 .and. ifDebug == 2) &
                print *,xcen, ymodCen, ' eliminated by B boundaries'
          endif
          if (ifUse == 1 .and. nxSeries > 0) then
            if (.not. inside(xvertTilt, yvertTilt, 4, xcen, ymodCen)) ifuse = 0
            if (ifUse == 0 .and. ifDebug == 2) &
                print *,xcen, ymodCen, ' eliminated by A tilt series borders'
          endif
          if (ifUse == 1 .and. nxSeriesB > 0) then
            if (.not. inside(xvertTiltB, yvertTiltB, 4, xcen, ymodCen)) ifuse = 0
            if (ifUse == 0 .and. ifDebug == 2) &
                print *,xcen, ymodCen, ' eliminated by B tilt series borders'
          endif

          if (ifUse == 0) then
            ifDone(indP) = -1
          else
            ifDone(indP) = 0
            numPosTotal = numPosTotal + 1
          endif
          directErr(indP) = -1.
        enddo
      enddo
    enddo
    !
    if (numLsdBinnings > 0) then
      i = numXpatch * numYpatch * numZpatch
      write(*,'(i7,a,i7,a)') numElimBySD, ' of', i, &
          ' total possible patches eliminated by SD criterion'
      if ((numPosTotal < 1 .and. numElimBySD > 0) .or. &
          (numPosTotal < 4 .and. numElimBySD > i / 2) .or. &
          numElimBySD > sdElimMaxFrac * i) then
        elimBySdType = 0.
        print *,'That is too much elimination; trying again without SD criterion [CSD1]'
      else
        exit
      endif
    endif
  enddo
  if (numPosTotal < 1) call exitError('No patches fit within all of the constraints')
  !
  ! set indexes at which to load data and compose patches
  ! Here is the padding for the B patch direct correlation
  !
  nxPad = nxPatch + 2 * (maxShift + 1)
  nyPad = nyPatch + 2 * (maxShift + 1)
  nzPad = nzPatch + 2 * (maxShift + 1)
  numPatchPixels = nxPatch * nyPatch * nzPatch
  !
  ! Patch A space is set by Fourier correlation padding, B space by
  ! max of Fourier and direct padded volumes
  !
  indPatchA = 1
  indPatchB = indPatchA + (nxXCpad + 2) * nyXCpad * nzXCpad
  indScratch = indPatchB
  idim = max(nxPad * nyPad * nzPad, (nxXCpad + 2) * nyXCpad * nzXCpad)
  
  ! Separate scratch for filtering only needed if filtering B also
  ! if (sigmaKernel > 0.) indScratch = indPatchB + idim
  indLoadA = indScratch + idim
  !
  ! get maximum load size: exload could be a separate parameter
  ! First get it for the optimal memory usage; if that is not big enough,
  ! then get it for a higher usage; but in any case make it big enough
  loadExtra = 2 * (maxShift + 1)
  nyLoadEx = min(ny2, nyPatch + loadExtra)
  nzLoadEx = min(nz2, nzPatch + loadExtra)
  maxXload = (idimOptimal - indLoadA - loadExtra * nyLoadEx * nzLoadEx) / &
      (nyPatch * nzPatch + nyLoadEx * nzLoadEx)
  if (maxXload < nxPatch) then
    maxXload = (idimHigher - indLoadA - loadExtra * nyLoadEx * nzLoadEx) / &
        (nyPatch * nzPatch + nyLoadEx * nzLoadEx)
    maxXload = max(nxPatch + 4, maxXload)
  endif
  maxXload = min(nx + 2, maxXload)
  idim = indLoadA + maxXload * nyPatch * nzPatch +  &
      min(nx2 + 2, maxXload + loadExtra) * nyLoadEx * nzLoadEx + 10
  indLoadB = indLoadA + maxXload * nyPatch * nzPatch
  if (ifDebug > 0)  &
      print *,numPatchPixels, indpatcha, indpatchb, indloada, indloadb, maxxload
  allocate(buffer(idim), stat=ierr)
  call memoryError(ierr, 'image buffer')
  !
  call dopen(1, outputFile, 'new', 'f')
  write(1, '(i7,a,i3)') numPosTotal, ' positions', idCCCcol

  ! Set up number of threads based on analysis relative to 40x20x40 patches
  if (kernelSize == 3) then
    numSmoothThreads = max(1, min(6, nint(2 * (numPatchPixels / 32000.)**0.4)))
  else
    numSmoothThreads =  max(1, min(8, nint(4 * (numPatchPixels / 32000.)**0.6)))
  endif
  if (ifDebug > 0) print *,'smooth threads',numSmoothThreads
  numSmoothThreads = numOMPthreads(numSmoothThreads)
  num3CorrThreads = max(1, min(8, nint(2 * (numPatchPixels / 32000.)**0.45)))
  if (ifDebug > 0) print *,'search threads',num3CorrThreads
  num3CorrThreads = numOMPthreads(num3CorrThreads)
  numCCCThreads = max(1, min(8, nint(4 * (numPatchPixels / 32000.)**0.6)))
  if (ifDebug > 0) print *,'CCC threads',numCCCThreads
  numCCCThreads = numOMPthreads(numCCCThreads)
  if (ifDebug > 0) print *,'Actual threads:', numSmoothThreads, num3CorrThreads, &
      numCCCThreads
  !
  ! loop from center out in all directions, X inner, then short dimension,
  ! then long dimension of Y and Z
  !
  if (ifFlip .ne. 0) then
    call sequencePatches(numXpatch, numYpatch, numZpatch, ixSequence, iySequence, &
        izSequence, idirSequence, numSequence)
  else
    call sequencePatches(numXpatch, numZpatch, numYpatch, ixSequence, izSequence, &
        iySequence, idirSequence, numSequence)
  endif

  loadX1 = -1
  loadXb1 = -1
  numSkip = 0
  do indSequence = 1, numSequence
    ixPatch = ixSequence(indSequence)
    iyPatch = iySequence(indSequence)
    izPatch = izSequence(indSequence)
    ixDir = idirSequence(indSequence)
    loadFullWidth = ixPatch == ixSequence(1)
    iz0 = izStart + (izPatch - 1) * izDelta
    iz1 = iz0 + nzPatch - 1
    izCen = iz0 + nzPatch / 2
    iy0 = iyStart + (iyPatch - 1) * iyDelta
    iy1 = iy0 + nyPatch - 1
    iyCen = iy0 + nyPatch / 2
    ix0 = ixStart + (ixPatch - 1) * ixDelta
    ix1 = ix0 + nxPatch - 1
    ixCen = ix0 + nxPatch / 2
    ! print *,'doing', ixcen, iycen, izcen
    indP = indPatch(ixPatch, iyPatch, izPatch)
    if (ifDebug == 2 .and. ifDone(indP) == 0) then
      write(1, 105) ixCen, iyCen, izCen, 2., 2., 2.
      ifDone(indP) = 1
    endif

    if (ifDone(indP) == 0) then
      !
      ! find and average adjacent and nearby patches, weighting by the
      ! correlation coefficient
      !
      numAdjacent = 0
      numNear = 0
      numAdjacentLook = distNear + 1.
      dxSum = 0.
      dySum = 0.
      dzSum = 0.
      dxSumNear = 0.
      dySumNear = 0.
      dzSumNear = 0.
      wsumAdjacent = 0.
      wsumNear = 0.
      do ix = -numAdjacentLook, numAdjacentLook
        do iy = -numAdjacentLook, numAdjacentLook
          do iz = -numAdjacentLook, numAdjacentLook
            distSq = ix**2 + iy**2 + iz**2
            if (distSq > 0 .and. distSq <= distNear**2) then
              ixAdjacent = ixPatch + ix
              iyAdjacent = iyPatch + iy
              izAdjacent = izPatch + iz
              if (ixAdjacent > 0 .and. ixAdjacent <= numXpatch .and. &
                  iyAdjacent > 0 .and. iyAdjacent <= numYpatch .and. &
                  izAdjacent > 0 .and. izAdjacent <= numZpatch) then
                indA = indPatch(ixAdjacent, iyAdjacent, izAdjacent)
                if (ifDone(indA) > 0) then
                  numNear = numNear + 1
                  wcc = max(0.01, corrCoef(indA))
                  dxSumNear = dxSumNear + dxPatch(indA) * wcc
                  dySumNear = dySumNear + dyPatch(indA) * wcc
                  dzSumNear = dzSumNear + dzPatch(indA) * wcc
                  wsumNear = wsumNear + wcc
                  if (distSq <= distAdjacent**2) then
                    numAdjacent = numAdjacent + 1
                    dxSum = dxSum + dxPatch(indA) * wcc
                    dySum = dySum + dyPatch(indA) * wcc
                    dzSum = dzSum + dzPatch(indA) * wcc
                    wsumAdjacent = wsumAdjacent + wcc
                  endif
                endif
              endif
            endif
          enddo
        enddo
      enddo
      ! nadj=max(1, nadj)
      dxAdjacent = dxSum / max(.01, wsumAdjacent)
      dyAdjacent = dySum / max(.01, wsumAdjacent)
      dzAdjacent = dzSum / max(.01, wsumAdjacent)
      dxNear = dxSumNear / max(.01, wsumNear)
      dyNear = dySumNear / max(.01, wsumNear)
      dzNear = dzSumNear / max(.01, wsumNear)
      !
      ! Do direct search if there is something adjacent and its correlation
      ! coefficients are good enough and it is not very deviant
      !
      if (numAdjacent > 0 .and.  &
          wsumAdjacent / numAdjacent > ccRatio * wsumNear / numNear &
          .and. abs(dxNear - dxAdjacent) < maxShift &
          .and. abs(dyNear - dyAdjacent) < maxShift &
          .and. abs(dzNear - dzAdjacent) < maxShift) then

        call setBload(ix0, ix1, nx2, dxAdjacent, dxVolume, ixB0, ixB1)
        call setBload(iy0, iy1, ny2, dyAdjacent, dyVolume, iyB0, iyB1)
        call setBload(iz0, iz1, nz2, dzAdjacent, dzVolume, izB0, izB1)
        if (ix1 > ix0 .and. iy1 > iy0 .and. iz1 > iz0 .and. &
            (ix1 + 1 - ix0) * (iy1 + 1 - iy0) * (iz1 + 1 - iz0) >= numPatchPixels / 2) &
            then
          !
          ! revise parameters based on load
          !
          nxPatchTmp = ix1 + 1 - ix0
          nyPatchTmp = iy1 + 1 - iy0
          nzPatchTmp = iz1 + 1 - iz0
          nxPad = nxPatchTmp + 2 * (maxShift + 1)
          nyPad = nyPatchTmp + 2 * (maxShift + 1)
          nzPad = nzPatchTmp + 2 * (maxShift + 1)
          ! print *,'direct', ixcen, iycen, izcen, nxptmp, nyptmp, nzptmp
          call flush(6)
          !
          ! get the a patch from the loaded data into an exact fit,
          ! taper, pad, kernel filter optionally and set to zero mean
          !
          call loadExtractProcess(1, aScale, indLoadA, indPatchA, ix0, ix1, iy0, iy1, &
              iz0, iz1, 0, nxyz, loadX0, loadX1, nxLoad, loadY0, loadY1, nyLoad, &
              loadZ0, loadZ1, nzLoad, nxPatchTmp, nxPatchTmp, nyPatchTmp, nzPatchTmp,  &
              sigmaKernel)
          !
          ! get the b patch from the loaded data padded to maximum shift
          ! We only filter one patch because that is all that is needed to smooth the
          ! CCF in which we are extracting a peak; there, this smoothing is equivalent
          ! to smoothing both by half as much.  The CCC would be better if both were
          ! smoothed equivalently, but this is not worth it because the B patch takes
          ! much longer to smooth
          call loadExtractProcess(2, bScale, indLoadB, indPatchB, ixB0, ixB1, iyB0, &
              iyB1, izB0, izB1, loadExtra, nxyz2, loadXb0, loadXb1, nxLoadB, loadYb0, &
              loadYb1, nyLoadb, loadZb0, loadZb1, nzLoadB, nxPad, nxPad, nyPad, nzPad, &
              0.)
          !
          dxNew = 0.
          dyNew = 0.
          dzNew = 0.
          ! call dumpVolume(buf(indpatcha), nxptmp, nxptmp, nyptmp, nzptmp, &
          ! 'dumpa.')
          ! call dumpVolume(buf(indpatchb), nxpad, nxpad, nypad, nzpad, &
          ! 'dumpb.')
          wallStart = wallTime()
          call findBestCorr(buffer(indPatchA), nxPatchTmp, nyPatchTmp, nzPatchTmp, ix0, &
              iy0, iz0, buffer(indPatchB), nxPad, nyPad, nzPad, ix0 - maxShift - 1, &
              iy0 - maxShift - 1, iz0 - maxShift - 1, ix0, ix1, iy0, iy1, &
              iz0, iz1, dxNew, dyNew, dzNew, maxShift, found, numCorrs, num3CorrThreads)
          wallBest = wallBest + wallTime() - wallStart
          if (found) then
            ifDone(indP) = 1
            dxPatch(indP) = dxNew + nint(dxAdjacent)
            dyPatch(indP) = dyNew + nint(dyAdjacent)
            dzPatch(indP) = dzNew + nint(dzAdjacent)
            dxDirect(indP) = dxPatch(indP)
            dyDirect(indP) = dyPatch(indP)
            dzDirect(indP) = dzPatch(indP)
            !
            ! compute correlation coefficient
            !
            wallStart = wallTime()
            call oneCorrCoeff(buffer(indPatchA), nxPatchTmp, nyPatchTmp, nzPatchTmp, &
                buffer(indPatchB), nxPad, nyPad, nzPad, nxPatchTmp, nyPatchTmp, &
                nzPatchTmp, dxNew, dyNew, dzNew, corrCoef(indP), numCCCThreads)
            ccDirect(indP) = corrCoef(indP)
            wallCCC = wallCCC + wallTime() - wallStart
          endif
        else
          ifDone(indP) = -1
          numSkip = numSkip + 1
        endif
      endif
      !
      ! If there are no adjacent patches, or something was fishy,
      ! do a full cross corr
      !
      if (ifDone(indP) == 0 .or. ifDebug == 3) then

        if (loadX1 < 0 .and. loadXb1 < 0) then
          dxNear = dxInitial
          dyNear = dyInitial
          dzNear = dzInitial
        endif
        call setBload(ix0, ix1, nx2, dxNear, dxVolume, ixB0, ixB1)
        call setBload(iy0, iy1, ny2, dyNear, dyVolume, iyB0, iyB1)
        call setBload(iz0, iz1, nz2, dzNear, dzVolume, izB0, izB1)
        if (ix1 > ix0 .and. iy1 > iy0 .and. iz1 > iz0 .and. &
            (ix1 + 1 - ix0) * (iy1 + 1 - iy0) * (iz1 + 1 - iz0) >= numPatchPixels / 2) &
            then
          !
          ! Adjust load sizes and make new ctf if needed
          !
          nxPatchTmp = ix1 + 1 - ix0
          nyPatchTmp = iy1 + 1 - iy0
          nzPatchTmp = iz1 + 1 - iz0
          nxXCpad = niceFrame(nxPatchTmp + 2 * nxXCbord, 2, niceLim)
          nyXCpad = niceFrame(nyPatchTmp + 2 * nyXCbord, 2, niceLim)
          nzXCpad = niceFrame(nzPatchTmp + 2 * nzXCbord, 2, niceLim)
          ! print *,'XC', nxptmp, nyptmp, nzptmp, nxXCpad, nyXCpad, nzXCpad, &
          ! ix0, ix1, iy0, iy1, iz0, iz1, ixb0, ixb1, iyb0, iyb1, izb0, izb1
          call flush(6)
          if ((radius2 > 0. .or. sigma1 .ne. 0. .or. sigma2 .ne. 0) .and. &
              (nxCTF .ne. nxXCpad .or. nyCTF .ne. nyXCpad .or. &
              nzCTF .ne. nzXCpad)) then
            call setCTFwSR(sigma1, sigma2, 0., radius2, ctf, nxXCpad, &
                max(nyXCpad, nzXCpad), delta)
            nxCTF = nxXCpad
            nyCTF = nyXCpad
            nzCTF = nzXCpad
          endif

          !
          ! get the both patches from the loaded data padded for XCorr
          !
          call loadExtractProcess(1, aScale, indLoadA, indPatchA, ix0, ix1, iy0, iy1, &
              iz0, iz1, 0, nxyz, loadX0, loadX1, nxLoad, loadY0, loadY1, nyLoad,  &
              loadZ0, loadZ1, nzLoad, nxXCpad + 2, nxXCpad, nyXCpad, nzXCpad, &
              sigmaKernel)
          ! call dumpVolume(buf(indpatcha), nxXCpad + 2, nxXCpad, &
          ! nyXCpad, nzXCpad, 'dumpa.')
          call loadExtractProcess(2, bScale, indLoadB, indPatchB, ixB0, ixB1, iyB0, &
              iyB1, izB0, izB1, loadExtra, nxyz2, loadXb0, loadXb1, nxLoadB, loadYb0, &
              loadYb1, nyLoadb, loadZb0, loadZb1, nzLoadB, nxXCpad + 2, nxXCpad, &
              nyXCpad, nzXCpad, 0.)

          ! call dumpVolume(buf(indpatchb), nxXCpad + 2, nxXCpad, &
          ! nyXCpad, nzXCpad, 'dumpb.')
          call fourierCorr(buffer(indPatchA), buffer(indPatchB), (nxXCpad + 2) / 2, &
              nyXCpad, nzXCpad, work, ctf, delta)
          ! call dumpVolume(buf(indpatcha), nxXCpad + 2, nxXCpad, &
          ! nyXCpad, nzXCpad, 'dumpcorr.')

          call findXcorrPeak(buffer(indPatchA), nxXCpad + 2, nyXCpad, nzXCpad, &
              dxNew, dyNew, dzNew, peak)
          dxPatch(indP) = dxNew + nint(dxNear)
          dyPatch(indP) = dyNew + nint(dyNear)
          dzPatch(indP) = dzNew + nint(dzNear)
          !
          ! For coef, reload the patches without the 2-pixel X padding
          !
          call loadExtractProcess(1, aScale, indLoadA, indPatchA, ix0, ix1, iy0, iy1, &
              iz0, iz1, 0, nxyz, loadX0, loadX1, nxLoad, loadY0, loadY1, nyLoad, &
              loadZ0, loadZ1, nzLoad, nxXCpad, nxXCpad, nyXCpad, nzXCpad, sigmaKernel)
          call loadExtractProcess(2, bScale, indLoadB, indPatchB, ixB0, ixB1, iyB0, &
              iyB1, izB0, izB1, loadExtra, nxyz2, loadXb0, loadXb1, nxLoadB, loadYb0, &
              loadYb1, nyLoadb, loadZb0, loadZb1, nzLoadB, nxXCpad, nxXCpad, nyXCpad, &
              nzXCpad, 0.)
          call oneCorrCoeff(buffer(indPatchA), nxXCpad, nyXCpad, nzXCpad, &
              buffer(indPatchB), nxXCpad, nyXCpad, nzXCpad, nxPatchTmp, &
              nyPatchTmp, nzPatchTmp, dxNew, dyNew, dzNew, corrCoef(indP), numCCCThreads)

          if (ifDone(indP) == 0) numFourierPatch = numFourierPatch + 1
          if (ifDebug == 3 .and. ifDone(indP) > 0) &
              directErr(indP) = sqrt((dxDirect(indP) - dxPatch(indP))**2 + &
              (dyDirect(indP) - dyPatch(indP))**2 + &
              (dzDirect(indP) - dzPatch(indP))**2)
          ifDone(indP) = 1
        else if (ifDone(indP) == 0) then
          ifDone(indP) = -1
          numSkip = numSkip + 1
        endif
      endif
      if (ifDone(indP) > 0) then
        write(1, 105) ixCen, iyCen, izCen, dxPatch(indP), dyPatch(indP), &
            dzPatch(indP), corrCoef(indP)
105     format(3i6,3f9.2,f10.4,3f9.2,f12.4,f8.3)
        call flush(1)
      endif
    endif
  enddo
  !
  ! If any ones were skipped, need to rewrite the file
  !
  if ((numSkip > 0 .or. numLsdBinnings > 0 .or. ifDebug == 3) .and. ifDebug .ne. 2) then
    rewind(1)
    if (numLsdBinnings > 0) then
      write(1, '(i7,a,3i3)') numPosTotal - numSkip, ' positions', idCCCcol, idFracCol, &
          idStructCol
    else
      write(1, '(i7,a,3i3)') numPosTotal - numSkip, ' positions', idCCCcol
    endif
    do indSequence = 1, numSequence
      ixPatch = ixSequence(indSequence)
      iyPatch = iySequence(indSequence)
      izPatch = izSequence(indSequence)
      indP = indPatch(ixPatch, iyPatch, izPatch)
      if (ifDone(indP) > 0) then
        izCen = izStart + (izPatch - 1) * izDelta + nzPatch / 2
        iyCen = iyStart + (iyPatch - 1) * iyDelta + nyPatch / 2
        ixCen = ixStart + (ixPatch - 1) * ixDelta + nxPatch / 2
        if (ifDebug == 3) then
          if (directErr(indP) < 0.) then
            write(1, 105) ixCen, iyCen, izCen, dxPatch(indP), dyPatch(indP), &
                dzPatch(indP), corrCoef(indP)
          else
            write(1, 105) ixCen, iyCen, izCen, dxPatch(indP), dyPatch(indP), &
                dzPatch(indP), corrCoef(indP), dxDirect(indP), dyDirect(indP), &
                dzDirect(indP), ccDirect(indP), directErr(indP)
          endif
        else
          if (numLsdBinnings > 0) then
            write(1, 115) ixCen, iyCen, izCen, dxPatch(indP), dyPatch(indP), &
                dzPatch(indP), corrCoef(indP), patchFracHighSD(indP),  &
                patchMeanStruct(indP) / highSDlevel
115         format(3i6,3f9.2,f10.4,f8.4,f10.5)
          else
            write(1, 105) ixCen, iyCen, izCen, dxPatch(indP), dyPatch(indP), &
                dzPatch(indP), corrCoef(indP)
          endif
        endif
      endif
    enddo
  endif
  close(1)
  perPos = (3. *numCorrs) / max(1, numPosTotal - numFourierPatch)
  write(*,'(f8.2,a,i5,a)') perPos, ' correlations per position,'// &
      ' Fourier correlations computed', numFourierPatch, ' times'
  if (ifDebug > 0) write(*,'(a,8f8.4)')'LETKZ, search, oneCCC:',  &
      (wallCum(i),i=1,5), wallBest, wallCCC
  call exit(0)
95 call exitError('Reading transform file')

CONTAINS


  ! tiltVertexSetup constructs a bounding box from the original data and rotates it
  ! by the negative of the tilt axis angle.  In the data set that this was first tried
  ! with, the good area appeared to go about 60 pixels BEYOND the borders of this area.
  ! Even if the good data stops nearer to this boundary, it should be OK to let the patch
  ! centers go out to these borders.  Thus this just has xbord and ybord zero.
  !
  subroutine tiltVertexSetup(nxTilt, nyTilt, nxRec, nyRec, xvert, yvert)
    integer*4 nxTilt, nyTilt, nxRec, nyRec
    real*4 xvert(4), yvert(4), xbord, ybord
    real*4 rotMat(2,3)
    real*4 cosd, sind
    if (nxTilt == 0 .or. nxRec == 0) return
    xbord = 0.
    ybord = 0.
    xvert(1) = nxRec / 2. - (nxTilt - xbord) / 2.
    yvert(1) = nyRec / 2. - (nyTilt - ybord) / 2.
    xvert(2) = nxRec / 2. + (nxTilt - xbord) / 2.
    yvert(2) = yvert(1)
    xvert(3) = xvert(2)
    yvert(3) = nyRec / 2. + (nyTilt - ybord) / 2.
    xvert(4) = xvert(1)
    yvert(4) = yvert(3)
    rotMat(1, 1) = cosd(axisRotation)
    rotMat(1, 2) = sind(axisRotation)
    rotMat(2, 1) = -rotMat(1, 2)
    rotMat(2, 2) = rotMat(1, 1)
    rotMat(1:2, 3) = 0.;
    if (ifDebug > 0) write(*,'(a, 8f9.1)')'raw verts', (xvert(ix), yvert(ix), ix = 1, 4)
    do ix = 1, 4
      call xfApply(rotMat, nxRec / 2., nyRec / 2., xvert(ix), yvert(ix), xvert(ix),  &
          yvert(ix))
    enddo
    if (ifDebug > 0) write(*,'(a, 8f9.1)')'tilt verts', (xvert(ix), yvert(ix), ix = 1, 4)
    return
  end subroutine tiltVertexSetup

  ! loadExtractProcess takes care of loading data as needed from the
  ! given unit IUNIT into the right area of BUFFER, extracting the desired
  ! patch from there, tapering it, and smoothing via a scratch patch
  ! area if specified
  !
  subroutine loadExtractProcess(iunit, scale, indLoadA, indPatchA, ix0, ix1, iy0, iy1, &
      iz0, iz1, loadExtra, nxyz, loadX0, loadX1, nxLoad, loadY0, loadY1, nyLoad, &
      loadZ0, loadZ1, nzLoad, nxPadDim, nxPad, nyPad, nzPad, sigmaKernel)
    implicit none
    real*4 sigmaKernel, scale
    integer*4 indLoadA, indPatchA, ix0, ix1, iy0, iy1, iz0, iz1
    integer*4 loadExtra, nxyz(3), loadX0, loadX1, nxLoad
    integer*4 loadY0, loadY1, nyLoad, loadZ0, loadZ1, nzLoad
    integer*4 nxPadDim, nxPad, nyPad, nzPad
    integer*4 indTaper, iunit
    real*8 wallTime, wallNow

    if (ifDebug > 0) wallStart = wallTime()
    call manageLoad(iunit, buffer(indLoadA), ix0, ix1, iy0, iy1, iz0, iz1, &
        loadExtra / 2, ixDir, ixDelta, maxXload, loadFullWidth, nxyz, loadX0, loadX1, &
        nxLoad, loadY0, loadY1, nyLoad, loadZ0, loadZ1, nzLoad)
    if (ifDebug > 0) then
      wallNow = wallTime()
      wallCum(1) = wallCum(1) + wallNow - wallStart
      wallStart = wallNow
    endif
    !
    ! get the patch from the loaded data and scale it at the same time; taper inside
    ! and shift mean to zero
    !
    call extractPatch(buffer(indLoadA), nxLoad, nyLoad, nzLoad, loadX0, loadY0, loadZ0, &
        ix0, iy0, iz0, buffer(indPatchA), nxPatchTmp, nyPatchTmp, nzPatchTmp, scale)
    if (ifDebug > 0) then
      wallNow = wallTime()
      wallCum(2) = wallCum(2) + wallNow - wallStart
      wallStart = wallNow
    endif
    indTaper = indPatchA
    if (sigmaKernel > 0.) indTaper = indScratch
    call taperInVol(buffer(indPatchA), nxPatchTmp, nyPatchTmp, nzPatchTmp, &
        buffer(indTaper), nxPadDim, nxPad, nyPad, nzPad, nxTaper, nyTaper, nzTaper)
    if (ifDebug > 0) then
      wallNow = wallTime()
      wallCum(3) = wallCum(3) + wallNow - wallStart
      wallStart = wallNow
    endif

    if (sigmaKernel > 0.) then
      !call dumpVolume(buffer(indtaper), nxpadDim, nxpad, nypad, nzpad, 'taper.')
      call kernelSmooth(buffer(indTaper), buffer(indPatchA), nxPadDim, nxPad, nyPad, &
          nzPad, kernelSize, sigmaKernel, numSmoothThreads)
      !call dumpVolume(buffer(indpatcha), nxpadDim, nxpad, nypad, nzpad, 'smooth.')
      if (ifDebug > 0) then
        wallNow = wallTime()
        wallCum(4) = wallCum(4) + wallNow - wallStart
        wallStart = wallNow
      endif
    endif
    call volMeanZero(buffer(indPatchA), nxPadDim, nxPad, nyPad, nzPad)
    if (ifDebug > 0) then
      wallNow = wallTime()
      wallCum(5) = wallCum(5) + wallNow - wallStart
    endif
    return
  end subroutine loadExtractProcess


  ! analyzeLocalSDs sets up parameters for the multi-binning local SD analysis,
  ! calls the routines to obtain the statistics, and analyzes the edge values and the
  ! overall distribution to set threshold values
  !
  subroutine analyzeLocalSDs()
    integer*4 numPix, lsdZind, numBoxes, ibin, numForMed, lsdFuncData(4), idelSamp
    integer*4 numSample, numPatchStructLoop, loopPatchStruct
    real*4 edgeMedian(LIM_LOCAL_BIN), edgeMADN(LIM_LOCAL_BIN), fallbackSD(LIM_LOCAL_BIN)
    real*4 histDip(LIM_LOCAL_BIN), peakBelow(LIM_LOCAL_BIN), peakAbove(LIM_LOCAL_BIN)
    real*4 histStart, histEnd, ratio, ratioLast, ratioDiff, boxStructCrit
    real*4 diffLast
    logical useFallback
    integer*4 multiBinSetup, multiBinStats, lsdLoadFunc, findHistogramDip
    real*8 percentileFloat
    external lsdLoadFunc
    !
    ! Set up sizes and spacings
    do j = 1, numLsdBinnings
      numPix = 1
      do i = 1, 3
        lsdBoxSize(i, j) = max(1, nint((float(lsdBoxSize(i, 1)) * lsdBinning(i, 1)) / &
            lsdBinning(i, j)))
        numPix = numPix * lsdBoxSize(i, j)
        lsdSpacing(i, j) = max(1, nint(min(lsdBoxSize(i, j) / 2., nxyzPatch(i) / 10.)))
      enddo
      if (numPix < 15) call exitError('Box size for local SDs is too small')
    enddo

    lsdStartCoord(1) = ixStart
    lsdStartCoord(2) = iyStart
    lsdStartCoord(3) = izStart
    lsdEndCoord(1) = ixStart + (numXpatch - 1) * ixDelta + nxPatch
    lsdEndCoord(2) = iyStart + (numYpatch - 1) * iyDelta + nyPatch
    lsdEndCoord(3) = izStart + (numZpatch - 1) * izDelta + nzPatch
    lsdZind = 3
    if (ifFlip > 0) lsdZind = 2
    lsdStartCoord(lsdZind) = 0
    lsdEndCoord(lsdZind) = nxyz(lsdZind) - 1
    ierr = multiBinSetup(lsdBinning, lsdBoxSize, lsdSpacing, numLsdBinnings,  &
        lsdStartCoord, lsdEndCoord, lsdBoxStart, numLsdBoxes, lsdBufferStarts,  &
        lsdStatStarts)
    if (ierr .ne. 0) call exitError('Setting up multi-bin analysis')
    numBoxes = 0
    do ibin = 1, numLsdBinnings
      numBoxes = max(numBoxes, lsdStatStarts(ibin + 1) - lsdStatStarts(ibin))
    enddo
    numSample = min(numBoxes, 2000 * numHistBins)
    ibin = max(lsdBufferStarts(numLsdBinnings + 1),  &
        numSample + max(numSample, numHistBins))
    allocate(statBuffer(ibin), statMeans(lsdStatStarts(numLsdBinnings + 1)),  &
        statSDs(lsdStatStarts(numLsdBinnings + 1)), patchFracHighSD(limPatch),  &
        patchMeanStruct(limPatch), stat = ierr)
    call memoryError(ierr, 'arrays for multi-bin analysis')
    lsdFuncData(1) = lsdStartCoord(1)
    lsdFuncData(2) = lsdEndCoord(1)
    lsdFuncData(3) = lsdStartCoord(2)
    lsdFuncData(4) = lsdEndCoord(2)
    wallStart = wallTime()
    ierr = multiBinStats(lsdBinning, lsdBoxSize, lsdSpacing, numLsdBinnings,  &
        lsdStartCoord, lsdEndCoord, lsdBoxStart, numLsdBoxes, lsdBufferStarts,  &
        lsdStatStarts, statBuffer, statMeans, statSDs, lsdFuncData, lsdLoadFunc)
    if (ierr .ne. 0) call exitError('Doing multi-bin analysis of SDs')
    if (ifDebug > 0) write(*,'(a,f9.3)')'multiBinStats time', wallTime() - wallStart
    
    !
    ! Collect values on the top and bottom planes
    useFallback = .true.
    do ibin = 1, numLsdBinnings
      !
      ! Collect values on the top and bottom planes and analyze lowest fraction for 
      ! median/MADN
      do i = 1, 3
        ibStart(i) = 1
        ibEnd(i) = numLsdBoxes(i, ibin)
      enddo
      ibStart(lsdZind) = ibEnd(lsdZind)
      numInSample = 0
      call addBoxesToSample(ibin)
      ibStart(lsdZind) = 1
      ibEnd(lsdZind) = 1
      call addBoxesToSample(ibin)
      !
      ! Sort the edge values and take median/MADN of a fraction of them
      call rsSortFloats(statBuffer, numInSample)
      numForMed = max(2., edgeSampleFrac * numInSample)
      call rsMedianOfSorted(statBuffer, numForMed, edgeMedian(ibin))
      call rsFastMADN(statBuffer, numForMed, edgeMedian(ibin), statBuffer(numForMed + 1),&
          edgeMADN(ibin))
      write(*,'(/,a,i2,a,f9.2,a,f8.2,/,a)')'For scaling', ibin, ': edge SD median = ', &
          edgeMedian(ibin), ', MADN = ', edgeMADN(ibin), &
          'Analyzing histogram of square root of SDs:'

      ! Get a sample of data for large volumes to do percentile/histogram analysis
      numBoxes = lsdStatStarts(ibin + 1) - lsdStatStarts(ibin)
      numSample = min(numBoxes, 2000 * numHistBins)
      idelSamp = numBoxes / numSample
      ix = lsdStatStarts(ibin) + 1
      do i = 1, numSample
        statBuffer(i) = statSDs(ix)
        ix = ix + idelSamp
      enddo

      ! Find a moderate high percentile as a measure of typical strong structure
      fallbackSD(ibin) = percentileFloat(max(1, nint(fallbackPctFrac * numSample)), &
          statBuffer, numSample)

      ! Take the square root to spread out the lower end of the histogram
      ! Find a high-percentile cutoff value for limiting the histogram range
      statBuffer(1:numSample) = sqrt(statBuffer(1:numSample))
      !histStart = minval(statBuffer(1:numSample))
      histStart = percentileFloat(max(1, nint(bottomHistFrac * numSample)), statBuffer, &
          numSample)
      histEnd = percentileFloat(max(1, nint(fracForHist * numSample)), statBuffer, &
          numSample)
      
      wallStart = wallTime()
      if (findHistogramDip(statBuffer, numSample, 0, statBuffer(numSample + 1),  &
          numHistBins, histStart, histEnd, histDip(ibin),  &
          peakBelow(ibin), peakAbove(ibin), ifHistVerbose) .ne. 0) then
        histDip(ibin) = -1
      else
        histDip(ibin) = histDip(ibin)**2
        peakBelow(ibin) = peakBelow(ibin)**2
        peakAbove(ibin) = peakAbove(ibin)**2
        write(*,'(a,f9.2,a,2f9.2)')'SD value of dip is ', histDip(ibin), ', peaks at ', &
            peakBelow(ibin), peakAbove(ibin)
        tmp = (peakBelow(ibin) - edgeMedian(ibin)) / edgeMADN(ibin)
        if (tmp > badHistPeakCrit) then
          histDip(ibin) = -1
          write(*,'(a,f9.2,a)')'Lower peak is', tmp, &
              ' MADNs above edge median so this dip is assumed to be wrong'
        else
          useFallback = .false.
        endif
      endif
      if (ifDebug > 0) write(*,'(i8,a,i12, a,f9.3)')numSample, ' of',numBoxes,  &
          ' boxes; histogram time', wallTime() - wallStart
    enddo
    print *

    ! If not using fallback and eliminating by SD, loop twice, first time see if too
    ! many are elimated and if so use fallbacks
    numPatchStructLoop = 2
    if (useFallback .or. nint(elimBySdType) < 1) numPatchStructLoop = 1
    do loopPatchStruct = 1, numPatchStructLoop

      ! Find most distinguishable binning based on ones with dips, or based on all
      ! fallbacks
      ratioLast = 0.
      ratioDiff = 0.
      ibinBest = -1
      do ibin = 1, numLsdBinnings
        if (useFallback) then
          ratio = (fallbackSD(ibin) - edgeMedian(ibin)) / edgeMADN(ibin)
        else
          if (histDip(ibin) < 0) cycle
          ratio = (peakAbove(ibin) - edgeMedian(ibin)) / edgeMADN(ibin)
        endif
        if (ibinBest < 0) then
          ibinBest = ibin
        else
          ratioDiff = (ratio - ratioLast) / (ibin - ibinBest)
          if (ratioDiff < betterScaleCrit * diffLast) cycle
          ibinBest = ibin
        endif
        ratioLast = ratio;
        diffLast = ratioDiff;
      enddo
      if (numLsdBinnings > 1) &
          write(*,'(a,i2,a,/)')'Selected scaling #', ibinBest, ' for further analysis'
      !
      ! Set up criteria now
      if (useFallback) then
        boxStructCrit = fallbackSD(ibinBest) * fallbackStructFrac
        highSDlevel = fallbackSD(ibinBest)
        write(*,'(a,f9.2,a,f9.2)')'No histogram dips found; using fallback SD of ', &
            fallbackSD(ibinBest), ' for strong structure and criterion of ',boxStructCrit
      else
        boxStructCrit = histDip(ibinBest)
        highSDlevel = peakAbove(ibinBest)
      endif
    !
      ! Now evaluate each patch, find mean SD of blocks in patch and fraction of blocks
      ! above the structure criterion.  Yes, this was actually good up to 10 threads
      wallStart = wallTime()
      ierr = numOMPthreads(10)
      
      !$OMP PARALLEL DO NUM_THREADS(ierr) &
      !$OMP& SHARED(numZpatch, izDelta, nzPatch, numYpatch, iyDelta, nyPatch, numXpatch)&
      !$OMP& SHARED(ixDelta, nxPatch, patchMeanStruct, patchFracHighSD, lsdStatStarts)&
      !$OMP& SHARED(ibinBest, numLsdBoxes, statSDs, boxStructCrit, izStart, iyStart)&
      !$OMP& SHARED(ixStart)&
      !$OMP& PRIVATE(izPatch, iz0, iyPatch, iy0, ixPatch, ix0, ibStart, ibEnd, numBoxes)&
      !$OMP& PRIVATE(ix, iy, iz, i, indP)
      do izPatch = 1, numZpatch
        iz0 = izStart + (izPatch - 1) * izDelta
        call findBoxesInsidePatch(3, iz0, iz0 + nzPatch - 1, ibStart(3), ibEnd(3))
        do iyPatch = 1, numYpatch
          iy0 = iyStart + (iyPatch - 1) * iyDelta
          call findBoxesInsidePatch(2, iy0, iy0 + nyPatch - 1, ibStart(2), ibEnd(2))
          do ixPatch = 1, numXpatch
            ix0 = ixStart + (ixPatch - 1) * ixDelta
            call findBoxesInsidePatch(1, ix0, ix0 + nxPatch - 1, ibStart(1), ibEnd(1))
            indP = indPatch(ixPatch, iyPatch, izPatch)
            numBoxes = 0
            patchMeanStruct(indP) = 0.;
            patchFracHighSD(indP) = 0.;
            do iz = ibStart(3), ibEnd(3)
              do iy = ibStart(2), ibEnd(2)
                do ix = ibStart(1), ibEnd(1)
                  i = lsdStatStarts(ibinBest) + ((iz - 1) * numLsdBoxes(2, ibinBest) + &
                      (iy - 1)) * numLsdBoxes(1, ibinBest) + ix
                  patchMeanStruct(indP) = patchMeanStruct(indP) + statSDs(i)
                  if (statSDs(i) >= boxStructCrit) &
                      patchFracHighSD(indP) = patchFracHighSD(indP) + 1
                  numBoxes = numBoxes + 1
                enddo
              enddo
            enddo
            if (numBoxes > 0) then
              patchFracHighSD(indP) = patchFracHighSD(indP) / numBoxes
              patchMeanStruct(indP) = patchMeanStruct(indP) / numBoxes
            else
              patchMeanStruct(indP) = -1.;
              patchFracHighSD(indP) = -1.;
            endif
          enddo
        enddo
      enddo
      !$OMP END PARALLEL DO
      if (ifDebug > 0) write(*,'(a,f9.3)') ' patch structure time', wallTime() - wallStart

      ! Count up the eliminations by the criterion
      if (loopPatchStruct < numPatchStructLoop) then
        numElimBySD = 0
        do iz = 1, numZpatch
          do iy = 1, numYpatch
            do ix = 1, numXpatch
              indP = indPatch(ix, iy, iz)
              !
              if ((nint(elimBySdType) > 1 .and. &
                  patchMeanStruct(indP) < elimBySdCrit * highSDlevel) .or. &
                  (nint(elimBySdType) == 1 .and. patchFracHighSD(indP) < elimBySdCrit)) &
                  then
                numElimBySD = numElimBySD + 1
              endif
            enddo
          enddo
        enddo

        ! Use fallback if too many
        i = numXpatch * numYpatch * numZpatch
        if (numElimBySD > sdElimMaxFrac * i) then
          useFallback = .true.
          write(*,'(i7,a,i7,a,/,a)') numElimBySD, ' of', i, &
              ' total possible patches would be eliminated by these SD criterion', &
              'That is too much elimination; recomputing with fallback criteria'
        else
          exit
        endif
      endif
    enddo

  end subroutine analyzeLocalSDs


  ! addBoxesToSample adds boxes in the range ibStart, ibEnd into a sample
  !
  subroutine addBoxesToSample(iscl)
    integer*4 iscl, ixBox, iyBox, izBox, ibase
    do izBox = ibStart(3), ibEnd(3)
      do iyBox = ibStart(2), ibEnd(2)
        ibase = lsdStatStarts(iscl) + ((izBox - 1) * numLsdBoxes(2, iscl) + (iyBox - 1)) &
            * numLsdBoxes(1, iscl)
        do ixBox = ibStart(1), ibEnd(1)
          numInSample = numInSample + 1
          statBuffer(numInSample) = statSDs(ibase + ixBox)
        enddo
      enddo
    enddo
  end subroutine addBoxesToSample


  ! findBoxesInsidePatch finds the box numbers in a patch whose starting and ending
  ! coordinates are in ipStart, ipEnd, and puts the range of boxes in ibStart, ibEnd
  !
  subroutine findBoxesInsidePatch(ixyz, ipStart, ipEnd, ibStart, ibEnd)
    integer*4 iscl, ixyz, ipStart, ipEnd, ibStart, ibEnd
    iscl = ibinBest
    ibStart = ((ipStart - lsdStartCoord(ixyz)) / float(lsdBinning(ixyz, iscl)) -  &
        lsdBoxStart(ixyz, iscl)) / lsdSpacing(ixyz, iscl) + 1
    do while (lsdStartCoord(ixyz) + lsdBinning(ixyz, iscl) * (lsdBoxStart(ixyz, iscl) + &
        (ibStart - 1) * lsdSpacing(ixyz, iscl)) < ipStart)
      ibStart = ibStart + 1
    enddo
    ibEnd = ceiling(((ipEnd - lsdStartCoord(ixyz)) / float(lsdBinning(ixyz, iscl)) &
        + 1 - lsdBoxStart(ixyz, iscl) - lsdBoxSize(ixyz, iscl)) /  &
        lsdSpacing(ixyz, iscl)) + 1
    ibEnd = min(numLsdBoxes(ixyz, iscl), ibEnd)
    do while (lsdStartCoord(ixyz) + lsdBinning(ixyz, iscl) * (lsdBoxStart(ixyz, iscl) + &
        (ibEnd - 1) * lsdSpacing(ixyz, iscl) + lsdBoxSize(ixyz, iscl)) < ipEnd)
      ibEnd = ibEnd - 1
    enddo
    return
  end subroutine findBoxesInsidePatch

end program corrsearch3d


! FIND_BEST_CORR will search for the best 3-D displacement between
! two volumes in ARRAY and BRRAY.  ARRAY is dimensioned to NXA by NYA
! by NZA, and its starting index coordinates are LOADAX0, LOADAY0,
! LOADAZ0.  BRRAY is dimensioned to NXB by NYB by NZB, and its
! starting index coordinates are LOADBX0, LOADBY0, LOADBZ0.  The volume
! to be correlated is specified by starting and ending index
! coordinates IX0, IX1, IY0, IY1, IZ0, IZ1.  DXADJ, DYADJ, DZADJ
! should contain a starting shift upon entry, and will return the
! shift with the best correlation.  MAXSHIFT specifies the maximum
! shift that is allowed.  FOUND is returned as TRUE if a correlation
! peak is found within the maximum shift.  NCORRS is a variable that
! allows the calling program to maintain a count of the total number
! of correlations computed.
!
! This was originally written to work with patches embedded in larger
! loaded volumes, with the two volumes potentially loaded separately,
! hence the unnecessary complexity for the program at hand.
!
subroutine findBestCorr(array, nxA, nyA, nza, loadAx0, loadAy0, loadAz0, brray, nxB, &
    nyB, nzb, loadBx0, loadBy0, loadBz0, ix0, ix1, iy0, iy1, iz0, iz1, dxAdjacent, &
    dyAdjacent, dzAdjacent, maxShift, found, numCorrs, numThreads)
  implicit none

  integer*4 nxA, nxB, nyA, nyB, nza, nzb, numThreads
  real*4 array(nxA,nyA,nza), brray(nxB,nyB,nzb)
  real*8 corrs(-1:1,-1:1,-1:1), corrTmp(-2:2,-2:2,-2:2)
  real*8 corrMax
  logical done(-1:1, -1:1, -1:1), doneTmp(-2:2, -2:2, -2:2)
  integer*4 idySequence(9) /0, -1, 1, 0, 0, -1, 1, -1, 1/
  integer*4 idzSequence(9) /0, 0, 0, 1, -1, -1, -1, 1, 1/
  logical found
  integer*4 ix0, ix1, iy0, iy1, iz0, iz1, maxShift, numCorrs
  integer*4 loadAx0, loadBx0, loadAy0, loadBy0, loadAz0, loadBz0
  real*4 dxAdjacent, dyAdjacent, dzAdjacent
  !
  integer*4 idxGlobal, idyGlobal, idzGlobal, ix, iy, iz, indSequence, idy, idz
  integer*4 idyCorr, idzCorr, ix0corr, ix1corr, iy0corr, iy1Corr, iz0corr, iz1cor
  integer*4 indMax
  real*4 cx, cy, cz, y1, y2, y3
  real*8 parabolicFitPosition
  !
  ! get global displacement of b, including the load offset
  !
  ! print *,'findbest', nxa, nya, nza, &
  ! loadax0, loaday0, loadaz0, nxb, nyb, nzb, &
  ! loadbx0, loadby0, loadbz0, ix0, ix1, iy0, iy1, iz0, iz1, dxadj, &
  ! dyadj, dzadj, maxshift
  idxGlobal = nint(dxAdjacent) + loadAx0 - loadBx0
  idyGlobal = nint(dyAdjacent) + loadAy0 - loadBy0
  idzGlobal = nint(dzAdjacent) + loadAz0 - loadBz0
  !
  ! clear flags for existence of corr
  !
  do iz = -1, 1
    do iy = -1, 1
      do ix = -1, 1
        done(ix, iy, iz) = .false.
      enddo
    enddo
  enddo

  corrMax = -1.e30
  indSequence = 1
  do while(indSequence <= 9)
    idy = idySequence(indSequence)
    idz = idzSequence(indSequence)
    ! print *,iseq, idy, idz
    if (.not.(done(-1, idy, idz) .and. done(0, idy, idz) .and. &
        done(1, idy, idz))) then
      !
      ! if the whole row does not exist, do the correlations
      ! limit the extent if b is displaced and near an edge
      !
      idyCorr = idyGlobal + idy
      idzCorr = idzGlobal + idz
      ix0corr = max(ix0 - loadAx0, -(idxGlobal - 1))
      ix1corr = min(ix1 - loadAx0, nxB - 1 - (idxGlobal + 1))
      iy0corr = max(iy0 - loadAy0, -idyGlobal)
      iy1Corr = min(iy1 - loadAy0, nyB - 1 - idyGlobal)
      iz0corr = max(iz0 - loadAz0, -idzGlobal)
      iz1cor = min(iz1 - loadAz0, nzb - 1 - idzGlobal)
      numCorrs = numCorrs + 1
      call threeCorrs(array, nxA, nyA, brray, nxB, nyB, ix0corr, ix1corr, &
          iy0corr, iy1Corr, iz0corr, iz1cor, idxGlobal, idyCorr, idzCorr, &
          corrs(-1, idy, idz), corrs(0, idy, idz), corrs(1, idy, idz), numThreads)
      done(-1, idy, idz) = .true.
      done(0, idy, idz) = .true.
      done(1, idy, idz) = .true.
    endif
    if (corrs(0, idy, idz) > corrs(-1, idy, idz) .and. &
        corrs(0, idy, idz) > corrs(1, idy, idz)) then
      indMax = 0
    elseif (corrs(-1, idy, idz) > corrs(0, idy, idz) .and. &
        corrs(-1, idy, idz) > corrs(1, idy, idz)) then
      indMax = -1
    else
      indMax = 1
    endif
    if (corrs(indMax, idy, idz) > corrMax) then
      ! print *,'moving by', indmax, idy, idz
      corrMax = corrs(indMax, idy, idz)
      !
      ! if there is a new maximum, shift the done flags and the existing
      ! correlations, and reset the sequence
      !
      idxGlobal = idxGlobal + indMax
      idyGlobal = idyGlobal + idy
      idzGlobal = idzGlobal + idz
      !
      ! but if beyond the limit, return failure
      !
      if (max(abs(idxGlobal + loadBx0 - loadAx0), abs(idyGlobal + loadBy0 - loadAy0), &
          abs(idzGlobal + loadBz0 - loadAz0)) > maxShift) then
        found = .false.
        return
      endif
      do iz = -1, 1
        do iy = -1, 1
          do ix = -1, 1
            doneTmp(ix, iy, iz) = .false.
          enddo
        enddo
      enddo
      do iz = -1, 1
        do iy = -1, 1
          do ix = -1, 1
            doneTmp(ix - indMax, iy - idy, iz - idz) = done(ix, iy, iz)
            corrTmp(ix - indMax, iy - idy, iz - idz) = corrs(ix, iy, iz)
          enddo
        enddo
      enddo
      do iz = -1, 1
        do iy = -1, 1
          do ix = -1, 1
            done(ix, iy, iz) = doneTmp(ix, iy, iz)
            corrs(ix, iy, iz) = corrTmp(ix, iy, iz)
          enddo
        enddo
      enddo
      if (indMax .ne. 0 .or. idy .ne. 0 .or. idz .ne. 0) indSequence = 0
    endif
    indSequence = indSequence + 1
  enddo
  !
  ! do independent parabolic fits in 3 dimensions
  !
  y1 = corrs(-1, 0, 0)
  y2 = corrs(0, 0, 0)
  y3 = corrs(1, 0, 0)
  cx = parabolicFitPosition(y1, y2, y3)
  y1 = corrs(0, -1, 0)
  y3 = corrs(0, 1, 0)
  cy = parabolicFitPosition(y1, y2, y3)
  y1 = corrs(0, 0, -1)
  y3 = corrs(0, 0, 1)
  cz = parabolicFitPosition(y1, y2, y3)
  !
  dxAdjacent = idxGlobal + cx + loadBx0 - loadAx0
  dyAdjacent = idyGlobal + cy + loadBy0 - loadAy0
  dzAdjacent = idzGlobal + cz + loadBz0 - loadAz0
  found = .true.
  ! print *,'returning a peak'
  return
end subroutine findBestCorr


! THREECORRS computes three correlations between volumes in ARRAY
! and BRRAY.  The volume in ARRAY is dimensioned to NXA by NYA, that
! in BRRAY is dimensioned to NXB by NYB.  The starting and ending
! index coordinates (numbered from 0) in ARRAY over which the
! correlations are to be computed are IX0, IX1, IY0, IY1, IZ0, IZ1.
! The shift between coordinates in ARRAY and coordinates in B is given
! by IDX, IDY, IDZ.  The three correlations are returned in CORR1 (for
! IDX-1), CORR2 (for IDX), and CORR3 (for IDX+1) .
!
! On both the SGI and the PC, three correlations can be computed at
! once almost as fast as one can.  This technique thus speeds up the
! overall search by almost a factor of 3.
!
subroutine threeCorrs(array, nxA, nyA, brray, nxB, nyB, ix0, ix1, iy0, iy1, iz0, iz1, &
    idx, idy, idz, corr1, corr2, corr3, numThreads)
  implicit none
  real*4 array(*), brray(*)
  real*8 sum1, sum2, sum3, corr1, corr2, corr3
  integer*4 ix0, ix1, iy0, iy1, iz0, iz1, idx, idy, idz, nxA, nxB, nyA, nyB, numThreads
  integer*4 iz, izB, iy, iyb, indBaseA, indDelB, ix, ixB, nsum
  sum1 = 0.
  sum2 = 0.
  sum3 = 0.

  !$OMP PARALLEL DO NUM_THREADS(numThreads) REDUCTION (+ : sum1, sum2, sum3) &
  !$OMP& SHARED(iz0, iz1, idz, iy0, iy1, idy, nxA, nyA, nxB, nyB, ix0, ix1, array, brray)&
  !$OMP& PRIVATE(iz, izB, iy, iyb, indBaseA, indDelB, ix, ixB)
  do iz = iz0, iz1
    izB = iz + idz
    do iy = iy0, iy1
      iyb = iy + idy
      indBaseA = 1 + iy * nxA + iz * nxA * nyA
      indDelB = 1 + iyb * nxB + izB * nxB * nyB + idx - indBaseA

      do ix = indBaseA + ix0, indBaseA + ix1
        ixB = ix + indDelB
        sum1 = sum1 + array(ix) * brray(ixB - 1)
        sum2 = sum2 + array(ix) * brray(ixB)
        sum3 = sum3 + array(ix) * brray(ixB + 1)
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO
  nsum = (iz1 + 1 - iz0) * (iy1 + 1 - iy0) * (ix1 + 1 - ix0)
  corr1 = sum1 / nsum
  corr2 = sum2 / nsum
  corr3 = sum3 / nsum
  ! print *,idx, idy, idz, corr1, corr2, corr3
  return
end subroutine threeCorrs

! This was a more straightforward version for testing THREECORRS
!
!!$    subroutine safecorrs(array, nxa, nya, brray, nxb, nyb, ix0, ix1, &
!!$                iy0, iy1, iz0, iz1, idx, idy, idz, corr1, corr2, corr3)
!!$    real*4 array(nxa,nya,*), brray(nxb,nyb,*)
!!$    real*8 sum1, sum2, sum3, corr1, corr2, corr3
!!$    c
!!$    c       print *,nxa, nya, nxb, nyb, ix0, ix1, iy0, iy1, iz0, iz1
!!$    sum1=0.
!!$    sum2=0.
!!$    sum3=0.
!!$    do iz=iz0+1, iz1+1
!!$    izb=iz+idz
!!$    do iy=iy0+1, iy1+1
!!$    iyb=iy+idy
!!$    c             print *,array(ix0+1, iy, iz), brray(ix0+idx+1, iyb, izb)
!!$    do ix=ix0+1, ix1+1
!!$    ixb=ix+idx
!!$    c             print *,array(ix, iy, iz), brray(ixb, iyb, izb)
!!$    sum1=sum1+array(ix, iy, iz)*brray(ixb-1, iyb, izb)
!!$    sum2=sum2+array(ix, iy, iz)*brray(ixb, iyb, izb)
!!$    sum3=sum3+array(ix, iy, iz)*brray(ixb+1, iyb, izb)
!!$    enddo
!!$    enddo
!!$    enddo
!!$    nsum=(iz1+1-iz0)*(iy1+1-iy0)*(ix1+1-ix0)
!!$    corr1=sum1/nsum
!!$    corr2=sum2/nsum
!!$    corr3=sum3/nsum
!!$    c       print *,idx, idy, idz, corr1, corr2, corr3
!!$    return
!!$    end


! oneCorrCoeff computes one correlation coefficient between volumes in
! ARRAY and BRRAY.  The volume in ARRAY is dimensioned to NXA by NYA by
! NZA, that in BRRAY is dimensioned to NXB by NYB by NZB.  The size of
! volume to be correlated is NXPATCH by NYPATCH by NZPATCH and it is
! assumed to be centered in the arrays.  The shift between coordinates
! in ARRAY and coordinates in B is given by DX, DY, DZ.  The correlation
! coefficient is returned in CORR2.
!
subroutine oneCorrCoeff(array, nxA, nyA, nza, brray, nxB, nyB, nzb, nxPatch, nyPatch, &
    nzPatch, dx, dy, dz, corr2, numThreads)
  implicit none
  real*4 array(*), brray(*), corr2, dx, dy, dz, denom
  real*8 aSum, sum2, bsum2, asumSq, bsumSq2
  integer*4 ix0, ix1, iy0, iy1, iz0, iz1, idx, idy, idz, nxA, nxB, nyA, nyB, nza, nzb
  integer*4 iz, izB, iy, iyb, indBaseA, indDelB, ix, ixB, nsum
  integer*4 nxPatch, nyPatch, nzPatch, numThreads
  sum2 = 0.
  aSum = 0.
  bsum2 = 0.
  asumSq = 0.
  bsumSq2 = 0.
  idx = nint(dx) + (nxB - nxA) / 2
  ix0 = (nxA - nxPatch) / 2
  ix1 = min(nxPatch + ix0 - 1, nxB - 1 - idx)
  ix0 = max(ix0, -idx)
  idy = nint(dy) + (nyB - nyA) / 2
  iy0 = (nyA - nyPatch) / 2
  iy1 = min(nyPatch + iy0 - 1, nyB - 1 - idy)
  iy0 = max(iy0, -idy)
  idz = nint(dz) + (nzb - nza) / 2
  iz0 = (nza - nzPatch) / 2
  iz1 = min(nzPatch + iz0 - 1, nzb - 1 - idz)
  iz0 = max(iz0, -idz)

  !$OMP PARALLEL DO NUM_THREADS(numThreads) &
  !$OMP& REDUCTION(+ : sum2, aSum, asumSq, bsum2, bsumSq2) &
  !$OMP& SHARED(ix0, ix1, iy0, iy1, iz0, iz1, idx, idy, idz, nxA, nxB, nyA, nyB) &
  !$OMP& SHARED(array, brray) &
  !$OMP& PRIVATE(iz, izB, iy, iyb, indBaseA, indDelB, ix, ixB)
  do iz = iz0, iz1
    izB = iz + idz
    do iy = iy0, iy1
      iyb = iy + idy
      indBaseA = 1 + iy * nxA + iz * nxA * nyA
      indDelB = 1 + iyb * nxB + izB * nxB * nyB + idx - indBaseA

      do ix = indBaseA + ix0, indBaseA + ix1
        ixB = ix + indDelB
        sum2 = sum2 + array(ix) * brray(ixB)
        aSum = aSum + array(ix)
        asumSq = asumSq + array(ix)**2
        bsum2 = bsum2 + brray(ixB)
        bsumSq2 = bsumSq2 + brray(ixB)**2
      enddo
    enddo
  enddo
  !$OMP END PARALLEL DO
  nsum = (iz1 + 1 - iz0) * (iy1 + 1 - iy0) * (ix1 + 1 - ix0)
  denom = (nsum * asumSq - aSum**2) * (nsum * bsumSq2 - bsum2**2)
  !
  ! Set the ccc to 0 if the denominator is illegal, otherwise limit it
  ! to +/-1
  if (denom <= 0.) then
    corr2 = 0.
  else
    denom = sqrt(denom)
    corr2 = (nsum * sum2 - aSum * bsum2)
    if (denom < corr2) then
      corr2 = sign(1., corr2)
    else
      corr2 = corr2 / denom
    endif
  endif
  ! print *,idx, idy, idz, corr2, denom
  return
end subroutine oneCorrCoeff


! kernelSmooth applies a 3D gaussian kernel to the data in ARRAY and
! places the result in BRRAY.  The image size is NX x NY x NZ and the
! dimensions of the arrays are NXDIM x NY x NZ.  ISIZE specifies the
! kernel size (3 or 5) and sigma is the standard deviation of the
! Gaussian
!
subroutine kernelSmooth(array, brray, nxDim, nx, ny, nz, isize, sigma, numThreads)
  implicit none
  integer*4 nx, nxDim, ny, nz, ix, iy, iz, i, j, k, isize, mid, less, numThreads
  real*4 array(nxDim, ny, nz), brray(nxDim, ny, nz), sigma
  real*4 w(5,5,5), wsum
  !
  ! Make up the gaussian kernel
  !
  mid = (isize + 1) / 2
  wsum = 0.
  do i = 1, isize
    do j = 1, isize
      do k = 1, isize
        w(i, j, k) = exp(-((i - mid)**2 + (j - mid)**2 + (k - mid)**2) / &
            sigma**2)
        wsum = wsum + w(i, j, k)
      enddo
    enddo
  enddo
  do i = 1, isize
    do j = 1, isize
      do k = 1, isize
        w(i, j, k) = w(i, j, k) / wsum
      enddo
    enddo
  enddo
  !
  ! Form the weighted sums: moving this to a subroutine was needed to keep the
  ! single-thread performance from being 2x slower than non-OMP
  !
  !$OMP PARALLEL DO NUM_THREADS(numThreads) &
  !$OMP& SHARED(nz, isize, ny, nx, nxDim, mid, array, brray, w)
  do iz = 0, nz - isize
    call smoothOnePlane(array, brray, nxDim, nx, ny, iz, isize, w)
  enddo
  !$OMP END PARALLEL DO
  !
  ! Copy the walls of the volume
  !
  less = (mid - 1) / 2
  do iz = 1, nz - less, nz - less - 1
    do iy = 1, ny
      brray(1:nx, iy, iz) = array(1:nx, iy, iz)
      if (isize > 3) then
        brray(1:nx, iy, iz + 1) = array(1:nx, iy, iz + 1)
      endif
    enddo
  enddo
  do iy = 1, ny - less, ny - less - 1
    do iz = 1, nz
      brray(1:nx, iy, iz) = array(1:nx, iy, iz)
      if (isize > 3) then
        brray(1:nx, iy + 1, iz) = array(1:nx, iy + 1, iz)
      endif
    enddo
  enddo
  do ix = 1, nx - less, nx - less - 1
    do iz = 1, nz
      brray(ix, 1:ny, iz) = array(ix, 1:ny, iz)
      if (isize > 3) then
        brray(ix + 1, 1:ny, iz) = array(ix + 1, 1:ny, iz)
      endif
    enddo
  enddo
  return
end subroutine kernelSmooth

! Form the weighted sums for one Z plane: this formulation is 3 times faster
! than adding everything into one output pixel at a time
subroutine smoothOnePlane(array, brray, nxDim, nx, ny, iz, isize, w)
  implicit none
  integer*4 nx, nxDim, ny, iy, iz, j, k, isize, mid
  real*4 array(nxDim, ny, *), brray(nxDim, ny, *)
  real*4 w(5,5,5)
  mid = (isize + 1) / 2
  do iy = 0, ny - isize
    brray(mid:mid+nx-isize, iy + mid, iz + mid) = 0.
    do k = 1, isize
      do j = 1, isize
        if (isize == 3) then
          brray(2:nx-1, iy + 2, iz + 2) =  brray(2:nx-1, iy + 2, iz + 2) + &
              array(1:nx-2, iy + j, iz + k) * w(1, j, k) + &
              array(2:nx-1, iy + j, iz + k) * w(2, j, k) + &
              array(3:nx, iy + j, iz + k) * w(3, j, k)
        else
          brray(3:nx-2, iy + 3, iz + 3) =  brray(3:nx-2, iy + 3, iz + 3) + &
              array(1:nx-4, iy + j, iz + k) * w(1, j, k) + &
              array(2:nx-3, iy + j, iz + k) * w(2, j, k) + &
              array(3:nx-2, iy + j, iz + k) * w(3, j, k) + &
              array(4:nx-1, iy + j, iz + k) * w(4, j, k) + &
              array(5:nx, iy + j, iz + k) * w(5, j, k)
        endif
      enddo
    enddo
  enddo
end subroutine smoothOnePlane


! Computes a cross-correlation between volumes in ARRAY and BRRAY
! via fourier transforms and filters by the values in CTF if DELTA,
! the frequency spacing of the points in CTF, is nonzero.  NXDIM is
! the X dimension of the FFT, NY and NZ are the Y and Z image dimensions.
! WORK must be dimensioned at least NXDIM x NZ.
!
subroutine fourierCorr(array, brray, nxDim, ny, nz, work, ctf, delta)
  implicit none
  integer*4 nxDim, ny, nz
  real*4 ctf(*), delta, work
  complex array(nxDim, ny, nz), brray(nxDim, ny, nz)
  integer*4 jx, jy, jz, nxReal, indF
  real*4 delX, delY, delZ, xa, ya, za, s
  !
  ! Get nx in real space and take 3D FFT's
  !
  nxReal = 2 * (nxDim - 1)
  call thrdfft(array, work, nxReal, ny, nz, 0)
  call thrdfft(brray, work, nxReal, ny, nz, 0)
  !
  ! multiply complex conjugate of array by brray, put back in array
  ! This is different from usual so that we will get the amount b is
  ! displaced from A, not the amount to shift B to align to A
  !
  do jz = 1, nz
    do jy = 1, ny
      do jx = 1, nxDim
        array(jx, jy, jz) = conjg(array(jx, jy, jz)) * brray(jx, jy, jz)
      enddo
    enddo
  enddo
  !
  ! Filter if delta set
  !
  if (delta > 0.) then
    delX = 0.5 / (nxDim - 1.)
    delY = 1. / ny
    delZ = 1. / nz
    do jz = 1, nz
      za = (jz - 1) * delZ
      if (za > 0.5) za = 1. - za
      do jy = 1, ny
        ya = (jy - 1) * delY
        if (ya > 0.5) ya = 1. - ya
        do jx = 1, nxDim
          xa = (jx - 1) * delX
          s = sqrt(xa**2 + ya**2 + za**2)
          indF = s / delta + 1.5
          array(jx, jy, jz) = array(jx, jy, jz) * ctf(indF)
        enddo
      enddo
    enddo
  endif
  call thrdfft(array, work, nxReal, ny, nz, -1)
  return
end subroutine fourierCorr


! findXcorrPeak finds the peak in a cross-correlation in ARRAY,
! dimensioned to NXDIM x NY x NZ and image size NXDIM-2, NY, NZ.
! It fits a parabola in each dimension to get interpolated peak
! positions in XPEAK, YPEAK, ZPEAK, and returns peak magnitude in PEAK.
!
subroutine findXcorrPeak(array, nxDim, ny, nz, xpeak, ypeak, zpeak, peak)
  implicit none
  integer*4 nxDim, ny, nz
  real*4 array(nxDim, ny, nz), xpeak, ypeak, zpeak, peak
  integer*4 ix, iy, iz, ixPeak, iyPeak, izPeak, nx
  real*4 cx, cy, cz, y1, y2, y3
  integer*4 indmap
  real*8 parabolicFitPosition
  !
  nx = nxDim - 2
  peak = -1.e30
  do iz = 1, nz
    do iy = 1, ny
      do ix = 1, nx
        if (array(ix, iy, iz) > peak) then
          peak = array(ix, iy, iz)
          ixPeak = ix
          iyPeak = iy
          izPeak = iz
        endif
      enddo
    enddo
  enddo
  !
  ! simply fit a parabola to the two adjacent points in X or Y or Z
  !
  y1 = array(indmap(ixPeak - 1, nx), iyPeak, izPeak)
  y2 = peak
  y3 = array(indmap(ixPeak + 1, nx), iyPeak, izPeak)
  cx = parabolicFitPosition(y1, y2, y3)

  y1 = array(ixPeak, indmap(iyPeak - 1, ny), izPeak)
  y3 = array(ixPeak, indmap(iyPeak + 1, ny), izPeak)
  cy = parabolicFitPosition(y1, y2, y3)

  y1 = array(ixPeak, iyPeak, indmap(izPeak - 1, nz))
  y3 = array(ixPeak, iyPeak, indmap(izPeak + 1, nz))
  cz = parabolicFitPosition(y1, y2, y3)
  !
  ! return adjusted pixel coordinate minus 1
  !
  xpeak = ixPeak + cx - 1.
  ypeak = iyPeak + cy - 1.
  zpeak = izPeak + cz - 1.
  if (xpeak > nx / 2) xpeak = xpeak - nx
  if (ypeak > ny / 2) ypeak = ypeak - ny
  if (zpeak > nz / 2) zpeak = zpeak - nz
  return
end subroutine findXcorrPeak

! setBload takes the desired coordinates on an axis, IX0 and IX1, the
! size in that dimension, NX2, the incremental and initial offsets
! DXADJ and DXINITIAL, and computes the limits for data that need
! to be loaded from B in IXB0, IXB1.  It adjusts IX0 and IX1 as
! necessary to keep everything within limits
!
subroutine setBload(ix0, ix1, nx2, dxAdjacent, dxVolume, ixB0, ixB1)
  implicit none
  integer*4 ix0, ix1, nx2, ixB0, ixB1, idxAdjacent, idxVolume
  real*4 dxAdjacent, dxVolume
  idxAdjacent = nint(dxAdjacent)
  idxVolume = nint(dxVolume)
  ixB0 = max(0, ix0 + idxAdjacent + idxVolume)
  ix0 = ixB0 - idxAdjacent - idxVolume
  ixB1 = min(nx2 - 1, ix1 + idxAdjacent + idxVolume)
  ix1 = ixB1 - idxAdjacent - idxVolume
  return
end subroutine setBload

! MANAGELOAD tests whether the desired volume specified by IX0, IX1,
! IY0, IY1, IZ0, IZ1 is already loaded, given the loaded limits in
! LOADX0, etc.  If not, it loads the data, with extra amounts specified
! by LOADEXH and a maximum load in X specified by MAXXLOAD
!
subroutine manageLoad(iunit, buffer, ix0, ix1, iy0, iy1, iz0, iz1, loadExtraHalf, &
    ixDir, ixDelta, maxXload, loadFullWidth, nxyz, loadX0, loadX1, nxLoad, loadY0, &
    loadY1, nyLoad, loadZ0, loadZ1, nzLoad)
  implicit none
  real*4 buffer(*)
  integer*4 ix0, ix1, iy0, iy1, iz0, iz1, loadExtraHalf, ixDir, ixDelta, maxXload
  integer*4 loadX0, loadX1, nxLoad, loadY0, loadY1, nyLoad, loadZ0, loadZ1
  integer*4 nzLoad, numMore, nxyz(3), iunit
  logical loadFullWidth
  !
  if (ix0 >= loadX0 .and. ix1 <= loadX1 .and. &
      iy0 >= loadY0 .and. iy1 <= loadY1 .and. &
      iz0 >= loadZ0 .and. iz1 <= loadZ1) return
  !
  ! need to load new data
  !
  loadY0 = max(0, iy0 - loadExtraHalf)
  loadY1 = min(nxyz(2) - 1, iy1 + loadExtraHalf)
  loadZ0 = max(0, iz0 - loadExtraHalf)
  loadZ1 = min(nxyz(3) - 1, iz1 + loadExtraHalf)
  !
  ! compute limits in X, loading as much as possible
  ! but limiting to edge of data and then truncating
  ! to the end of a patch
  !
  if (maxXload >= nxyz(1) .and. loadFullWidth) then
    loadX0 = 0
    loadX1 = nxyz(1) - 1
  else if (ixDir > 0) then
    loadX0 = max(0, ix0 - loadExtraHalf)
    loadX1 = min(nxyz(1) - 1, ix0 + maxXload-1 + loadExtraHalf)
    numMore = (loadX1 - loadExtraHalf - ix1) / ixDelta
    loadX1 = min(nxyz(1) - 1, ix1 + ixDelta * numMore + loadExtraHalf)
  else
    loadX1 = min(nxyz(1) - 1, ix1 + loadExtraHalf)
    loadX0 = max(0, ix1 + 1 - maxXload - loadExtraHalf)
    numMore = (ix0 - loadX0 + loadExtraHalf) / ixDelta
    loadX0 = max(0, ix0 - ixDelta * numMore - loadExtraHalf)
  endif
  ! write(*,'(a,i2,12i5)')'loading data', iunit, ix0, ix1, iy0, iy1, iz0, iz1, loadx0, &
  !  loadx1,  loady0, loady1, loadz0, loadz1
  nxLoad = loadX1 + 1 - loadX0
  nyLoad = loadY1 + 1 - loadY0
  nzLoad = loadZ1 + 1 - loadZ0
  call loadVol(iunit, buffer, nxLoad, nyLoad, loadX0, loadX1, loadY0, loadY1, &
      loadZ0, loadZ1)
  return
end subroutine manageLoad


! LOADVOL loads a subset of the volume from unit IUNIT, into ARRAY
! assuming dimensions of NXDIM by NYDIM, from index coordinates
! IX0, IX1, IY0, IY1, IZ0, IZ1.
!
subroutine loadVol(iunit, array, nxDim, nyDim, ix0, ix1, iy0, iy1, iz0, iz1)
  implicit none
  integer*4 nxDim, nyDim, ix0, ix1, iy0, iy1, iz0, iz1, iunit, indZ, iz
  real*4 array(nxDim,nyDim,*)
  !
  ! print *,iunit, nxdim, nydim, ix0, ix1, iy0, iy1, iz0, iz1
  indZ = 0
  do iz = iz0, iz1
    indZ = indZ + 1
    call imposn(iunit, iz, 0)
    call irdpas(iunit, array(1, 1, indZ), nxDim, nyDim, ix0, ix1, iy0, iy1,*99)
  enddo
  return
99 call exitError('Error reading file')
end subroutine loadVol


! EXTRACT_PATCH extracts a patch of dimensions NXPATCH by NYPATCH by
! NZPATCH into ARRAY from the loaded volume in BUF, whose dimensions
! are NXLOAD by NYLOAD by NZLOAD.  BUF is loaded from starting index
! coordinates LOADX0, LOADY0, LOADZ0, and the starting index
! coordinates of the patch are IX0, IY0, IZ0.  Values are multipled by SCALE.
!
subroutine extractPatch(buffer, nxLoad, nyLoad, nzLoad, loadX0, loadY0, loadZ0, ix0, &
    iy0, iz0, array, nxPatch, nyPatch, nzPatch, scale)
  implicit none
  integer*4 nxLoad, nyLoad, nzLoad, loadX0, loadY0, loadZ0, ix0, iy0
  integer*4 iz0, nxPatch, nyPatch, nzPatch
  real*4 buffer(nxLoad,nyLoad,nzLoad), array(nxPatch,nyPatch,nzPatch), scale
  integer*4 iz, iy, ix
  !
  ix = ix0 - loadX0
  iy = iy0 - loadY0
  iz = iz0 - loadZ0
  array(1:nxPatch, 1:nyPatch, 1:nzPatch) =  &
      buffer(1+ix:nxPatch+ix, 1+iy:nyPatch+iy, 1+iz:nzPatch+iz) * scale
  return
end subroutine extractPatch


! VOLMEANZERO shifts the mean to zero of the volume in ARRAY
! dimensioned NXDIM by NY by NZ, image size NX by NY by NZ
!
subroutine volMeanZero(array, nxDim, nx, ny, nz)
  implicit none
  integer*4 nxDim, nx, ny, nz
  real*4 array(nxDim,ny,nz)
  real*8 arsum
  real*4 dmean
  arsum = sum(array(1:nx, 1:ny, 1:nz))
  dmean = arsum / (nx * ny * nz)
  array(1:nx, 1:ny, 1:nz) = array(1:nx, 1:ny, 1:nz) - dmean
  return
end subroutine volMeanZero


! checkAndSetPatches does error checks and sets the basic start
! and delta for the patches in one dimension.
!
subroutine checkAndSetPatches(nx, nBordXlow, nbordXhigh, nxPatch, numXpatch, &
    ixStart, ixDelta, iaxis)
  implicit none
  integer*4 nx, nBordXlow, nbordXhigh, nxPatch, numXpatch, ixStart, ixDelta, iaxis
  character*6 axis
  !
  ! check basic input properties
  !
  axis = char(ichar('W') + iaxis) //' axis'
  if (nBordXlow < 0 .or. nbordXhigh < 0) call exitError( &
      'A negative border was entered for the '//axis)
  if (nxPatch <= 4) call exitError( &
      'Patch size negative or too small for the '//axis)
  if (numXpatch <= 0) call exitError( &
      'Number of patches must be positive for the '//axis)
  if (nxPatch > nx - (nBordXlow + nbordXhigh)) then
    write(*, '(/,a,i4,a,i4,a,a)') 'ERROR: CORRSEARCH3D -  Patch size (', &
        nxPatch, ') is bigger than specified range (', nx - (nBordXlow + nbordXhigh), &
        ') for the ', axis
    call exit(1)
  endif
  !
  ! If multiple patches, compute the delta and then adjust the number
  ! of patches down to require a delta of at least 2
  !
  if (numXpatch > 1) then
    ixStart = nBordXlow
    ixDelta = (nx - (nBordXlow + nbordXhigh + nxPatch)) / (numXpatch - 1)
    do while (numXpatch > 1 .and. ixDelta < 2)
      numXpatch = numXpatch - 1
      if (numXpatch > 1) &
          ixDelta = (nx - (nBordXlow + nbordXhigh + nxPatch)) / (numXpatch - 1)
    enddo
  endif
  !
  ! If only one patch originally or now, center it in range
  !
  if (numXpatch == 1) then
    ixStart = (nBordXlow + nx - nbordXhigh) / 2 - nxPatch / 2
    ixDelta = 1
  endif
  return
end subroutine checkAndSetPatches

! revisePatchRange adjusts the starting position and number of
! patches based on the vertex constraints in xvert2, xvert3
!
subroutine revisePatchRange(nx, nBordXlow, nbordXhigh, xvert2, xvert3, nxPatch, &
    numXpatch, ixStart, ixDelta)
  implicit none
  integer*4 nx, nBordXlow, nbordXhigh, numXpatch, ixStart, ixDelta, nxPatch
  real*4 xvert2, xvert3
  integer*4 ixlo2, ixhi2, ixSpan, numx2, newXdelta
  ixlo2 = max(nBordXlow, nint(xvert2))
  ixhi2 = min(nx - nbordXhigh, nint(xvert3))
  if (numXpatch > 1) then
    !
    ! get new number of intervals inside the limits, and a new delta
    ! to span the limits
    !
    ixSpan = max(0, ixhi2 - ixlo2 - nxPatch)
    numx2 = ixSpan / ixDelta + 1
    newXdelta = ixSpan / numx2
    if (newXdelta < 0.6 * ixDelta) then
      !
      ! but if delta is too small, drop the number of intervals
      !
      numx2 = numx2 - 1
      if (numx2 > 0) then
        newXdelta = ixSpan / numx2
      else
        newXdelta = ixDelta
      endif
    endif
    !
    ! now adjust ixlo2 up by half of remainder to center the patches
    ! in the span and find true start and end that fits inside the
    ! original low-hi limits
    !
    ixlo2 = ixlo2 + mod(ixSpan, newXdelta) / 2
    ixDelta = newXdelta
    ixStart = ixlo2 - ixDelta * ((ixlo2 - nBordXlow) / ixDelta)
    numXpatch = (nx - nbordXhigh - nxPatch - ixStart) / ixDelta + 1
  else
    !
    ! if only one patch, put it in new middle
    !
    ixStart = (ixlo2 + ixhi2) / 2 - nxPatch / 2
  endif
  return
end subroutine revisePatchRange

! Transforms a position XB, YB in B source to XA, YA in A.  NXYZBSRC and
! NXYZ are the dimensions of B and A, IFFLIPB and IFFLIP are 1 if the
! long dimension is Z in B or A, ASRC(3, 3) and DXYZSRC(3) have the
! 3D transformation
!
subroutine xformBsourceToA(xb, yb, nxyzBsource, nxyz, ifFlipB, ifFlip, aSource, &
    dxyzSource, xa, ya)
  implicit none
  real*4 xa, ya, xb, yb, aSource(3,3), dxyzSource(3), tmpB(3), tmpA(3)
  integer*4 nxyzBsource(3), nxyz(3), ifFlipB, ifFlip
  integer*4 i, j, indYb

  indYb = 2
  if (ifFlipB .ne. 0) indYb = 3
  tmpB(1) = xb
  tmpB(indYb) = yb
  tmpB(5 - indYb) = nxyzBsource(5 - indYb) / 2.
  do i = 1, 3
    tmpA(i) = dxyzSource(i) + nxyz(i) / 2.
    do j = 1, 3
      tmpA(i) = tmpA(i) + aSource(i, j) * (tmpB(j) - nxyzBsource(j) / 2.)
    enddo
    ! print *,(asrc(i, j), j = 1, 3), dxyzsrc(i), tmpb(i), tmpa(i)
  enddo
  xa = tmpA(1)
  ya = tmpA(2)
  if (ifFlip .ne. 0) ya = tmpA(3)
  return
end subroutine xformBsourceToA

! Fills arrays IXSEQ, IYSEQ, IZSEQ with a sequence of patch numbers
! starting from the center outward, progressing in X, then Y, then Z
! NUMXPAT, NUMYPAT, NUMZPAT is number of patches in each direction;
! NUMSEQ is returned with total number to loop on
!
subroutine sequencePatches(numXpatch, numYpatch, numZpatch, ixSequence, iySequence, &
    izSequence, idirSequence, numSequence)
  implicit none
  integer*4 numXpatch, numYpatch, numZpatch, ixSequence(*), iySequence(*), izSequence(*)
  integer*4 izPatchStart, izPatchEnd, izDir, idirSequence(*), numSequence
  integer*4 iyPatchStart, iyPatchEnd, iyDir, ixPatchStart, ixPatchEnd, ixDir, ixPatch
  integer*4 iyPatch, izPatch
  !
  izPatchStart = numZpatch / 2 + 1
  izPatchEnd = numZpatch
  numSequence = 0
  do izDir = 1, -1, -2
    do izPatch = izPatchStart, izPatchEnd, izDir
      iyPatchStart = numYpatch / 2 + 1
      iyPatchEnd = numYpatch
      do iyDir = 1, -1, -2
        do iyPatch = iyPatchStart, iyPatchEnd, iyDir
          ixPatchStart = numXpatch / 2 + 1
          ixPatchEnd = numXpatch
          do ixDir = 1, -1, -2
            do ixPatch = ixPatchStart, ixPatchEnd, ixDir
              numSequence = numSequence + 1
              ixSequence(numSequence) = ixPatch
              iySequence(numSequence) = iyPatch
              izSequence(numSequence) = izPatch
              idirSequence(numSequence) = ixDir
            enddo
            ixPatchStart = ixPatchStart - 1
            ixPatchEnd = 1
          enddo
        enddo
        iyPatchStart = iyPatchStart - 1
        iyPatchEnd = 1
      enddo
    enddo
    izPatchStart = izPatchStart - 1
    izPatchEnd = 1
  enddo
  return
end subroutine sequencePatches


! Dumps a volume from the contents of CRRAY with the given rootname
!
subroutine dumpVolume(crray, nxDim, nxPad, nyPad, nzPad, rootname)
  implicit none
  integer*4 nxDim, nxPad, nyPad, nzPad
  real*4 crray(nxDim,nyPad,nzPad)
  real*4 title(20), scale, dmin, dmax, tmin, tmax, dmt, brray(1000), dsub
  integer*4 kxyz(3), ix, iy, iz, mode
  character*(*) rootname
  integer*4 numFiles /1/
  save numFiles
  character*4 buffer
  character*160 filename

  mode = 0
  call int_iwrite(buffer, numFiles, iz)
  numFiles = numFiles + 1
  filename = trim(rootname)//trim(adjustl(buffer(1:iz)))
  call imopen(4, filename, 'new')
  !
  kxyz(1) = nxPad
  kxyz(2) = nyPad
  kxyz(3) = nzPad
  dmin = 1.e30
  dmax = -1.e30
  call icrhdr(4, kxyz, kxyz, mode, title, 0)
  call ialsiz_sam_cel(4, nxPad, nyPad, nzPad)
  do iz = 1, nzPad
    call iclden(crray(1, 1, iz), nxDim, nyPad, 1, nxPad, 1, nyPad, tmin, tmax, dmt)
    dmin = min(dmin, tmin)
    dmax = max(dmax, tmax)
  enddo
  if (mode == 0) then
    scale = 255. / (dmax - dmin)
    dsub = dmin
    dmin = 0
    dmax = 255
  else
    scale = 1.
    dsub = 0.
  endif
  do iz = 1, nzPad
    do iy = 1, nyPad
      do ix = 1, nxPad
        brray(ix) = scale * (crray(ix, iy, iz) - dsub)
      enddo
      call iwrlin(4, brray)
    enddo
  enddo
  !
  call iwrhdr(4, title, -1, dmin, dmax, 128.)
  call imclose(4)
  return
end subroutine dumpVolume


! The function for loading data for multibinstats
!
integer*4 function lsdLoadFunc(iz, idata, buffer)
  implicit none
  integer*4 iz, idata(*), iiuReadSecPart
  real*4 buffer
  call iiuSetPosition(1, iz, 0)
  lsdLoadFunc = iiuReadSecPart(1, buffer, idata(2) + 1 - idata(1), idata(1), idata(2), &
      idata(3), idata(4))
  return
end function lsdLoadFunc
