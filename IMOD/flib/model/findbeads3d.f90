!****   FINDBEADS3D    *********************************************
!
! Will find gold beads in a tomogram, given the size of the beads
!
! For details, see man page
!
! $Id$
!
program findbeads3d
  implicit none
  integer maxPiece, limHisto
  parameter (maxPiece = 400)
  parameter (limHisto = 10000)
  integer*4 lenPiece(maxPiece,3), ind0(maxPiece,3), ind1(maxPiece,3)
  integer*4 nx, ny, nz, nxyz(3), nxCorr, nyCorr, nzCorr, limPeak, maxArray
  equivalence (nx, nxyz(1)), (ny, nxyz(2)), (nz, nxyz(3))
  integer*4, allocatable :: indPeak(:), indCorr(:)
  real*4, allocatable :: peakVal(:), peakPos(:,:), corrVal(:), corrPos(:, :), array(:)
  real*4 origin(3), delta(3), histo(limHisto)
  integer*4 mxyz(3), mode, maxShift, i, index, loopCorr, numLook
  integer*4 nsum, nxOverlap, nyOverlap, nzOverlap, ierr, maxVol
  integer*4 minYsize, minZsize, maxZ, maxYZ, numXpieces, numYpieces
  integer*4 numZpieces, numPeaks, maxPeaks, numCorrs, indPlanes
  integer*4 ixStart, ixEnd, iyStart, iyEnd, izStart, izEnd
  real*4 elongation, radius, distMin, xpeak, ypeak, zpeak
  real*4 dmin, dmax, dmean, polarity, peakCorr, fracAvg, avgFallback, storeFallback
  integer*4 ixMin, ixMax, iyMin, iyMax, izMin, izMax
  integer*4 ix, ixPiece, iyPiece, izPiece, minInside, maxXsize
  integer*4 ixPeak, iyPeak, izPeak, numSave, ix0, ix1, iy0, iy1, iz0, iz1
  integer*4 numPass, minGuess, iVerbose, ibinning, num3CorrThreads
  real*4 sumXoffset, sumYoffset, sumZoffset, storeThresh, beadSize, degToRad
  real*4 histDip, peakBelow, peakAbove, dxAdjacent, dyAdjacent, dzAdjacent, avgThresh
  real*4 peakRelMin, sepMin, blackThresh, findCorr, aboveAvg, aboveSD

  logical loaded, yElongated, lightBeads, found, cleanBoth
  character*320 imageFile, modelFile, firstFile, tiltFile
  integer*4 findHistogramDip, numOMPthreads
  !
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetInteger, PipGetLogical
  integer*4 PipGetString, PipGetFloat, PipGetTwoFloats, PipGetTwoIntegers
  integer*4 PipGetInOutFile
  !
  ! fallbacks from ../../manpages/autodoc2man -3 2  findbeads3d
  !
  integer numOptions
  parameter (numOptions = 23)
  character*(40 * numOptions) options(1)
  options(1) = &
      'input:InputFile:FN:@output:OutputFile:FN:@candidate:CandidateModel:FN:@'// &
      'size:BeadSize:F:@binning:BinningOfVolume:I:@xminmax:XMinAndMax:IP:@'// &
      'yminmax:YMinAndMax:IP:@zminmax:ZMinAndMax:IP:@light:LightBeads:B:@'// &
      'angle:AngleRange:FP:@tilt:TiltFile:FN:@ylong:YAxisElongated:B:@'// &
      'peakmin:MinRelativeStrength:F:@threshold:ThresholdForAveraging:F:@'// &
      'store:StorageThreshold:F:@fallback:FallbackThresholds:FP:@'// &
      'spacing:MinSpacing:F:@both:EliminateBoth:B:@guess:GuessNumBeads:I:@'// &
      'max:MaxNumBeads:I:@verbose:VerboseOutput:I:@param:ParameterFile:PF:@help:usage:B:'
  !
  degToRad = 0.0174532
  modelFile = ' '
  firstFile = ' '
  elongation = 1.
  yElongated = .false.
  polarity = -1.
  lightBeads = .false.
  cleanBoth = .false.
  minInside = 64
  fracAvg = 0.5
  peakRelMin = 0.05
  numPass = 1
  maxPeaks = 50000
  sepMin = 0.9
  avgThresh = -2.
  storeThresh = 0.
  avgFallback = 0.
  storeFallback = 0.
  iVerbose = 0
  minGuess = 0
  ibinning = 1
  maxArray = 200000000
  !
  ! Pip startup: set error, parse options, do help output
  !
  call PipReadOrParseOptions(options, numOptions, 'findbeads3d', &
      'ERROR: FINDBEADS3D - ', .false., 2, 1, 1, numOptArg, numNonOptArg)
  !
  ! Open image file
  !
  if (PipGetInOutFile('InputFile', 1, ' ', imageFile) &
      .ne. 0) call exitError('No input file specified')

  ierr = PipGetInOutFile('OutputFile', 2, ' ', modelFile)
  call imopen(1, imageFile, 'RO')
  call irdhdr(1, nxyz, mxyz, mode, dmin, dmax, dmean)
  call iiuRetDelta(1, delta)
  call iiuRetOrigin(1, origin(1), origin(2), origin(3))

  if (PipGetFloat('BeadSize', beadSize) .ne. 0) call exitError( &
      'Bead diameter must be entered')
  ierr = PipGetInteger('BinningOfVolume', ibinning)
  if (ibinning < 1) call exitError('Binning must be positive')
  beadSize = beadSize / ibinning
  radius = beadSize / 2.

  ierr = PipGetLogical('YAxisElongated', yElongated)
  !
  ierr = PipGetFloat('MinRelativeStrength', peakRelMin)
  peakRelMin = sqrt(max(0.,peakRelMin))
  ierr = PipGetFloat('StorageThreshold', storeThresh)
  ierr = PipGetFloat('ThresholdForAveraging', avgThresh)
  !
  ierr = PipGetLogical('LightBeads', lightBeads)
  if (lightBeads) polarity = 1.
  ierr = PipGetString('CandidateModel', firstFile)
  !
  ierr = PipGetFloat('MinSpacing', sepMin)
  distMin = sepMin * beadSize
  ierr = PipGetLogical('EliminateBoth', cleanBoth)
  ierr = PipGetInteger('GuessNumBeads', minGuess)
  ierr = PipGetInteger('MaxNumBeads', maxPeaks)
  ierr = PipGetInteger('VerboseOutput', iVerbose)
  ierr = PipGetTwoFloats('FallbackThresholds', avgFallback, storeFallback)
  ierr = PipGetTwoFloats('AngleRange', dxAdjacent, dyAdjacent)
  ix = PipGetString('TiltFile', tiltFile)
  if (ierr == 0 .or. ix == 0) then
    if (ierr + ix == 0) call exitError( &
        'You cannot enter both an angle range and a tilt file')
    if (ix == 0) then
      call dopen(1, tiltFile, 'ro', 'f')
      dxAdjacent = 2000.
      dyAdjacent = -2000.
10    read(1,*,err = 97, end = 20) dzAdjacent
      dxAdjacent = min(dzAdjacent, dxAdjacent)
      dyAdjacent = max(dzAdjacent, dyAdjacent)
      go to 10
20    continue
    endif
    !
    ! Elongation factor from Radermacher 1988 paper
    ! cryoposition looks for 'Elongation factor is'
    dzAdjacent = 0.5 * (abs(dxAdjacent) + abs(dyAdjacent)) * degToRad
    elongation = sqrt((dzAdjacent + cos(dzAdjacent) * sin(dzAdjacent)) / &
        (dzAdjacent - cos(dzAdjacent) * sin(dzAdjacent)))
    print *,'Elongation factor is', elongation
  endif
  
  ixMin = 1
  iyMin = 1
  izMin = 1
  ixMax = nx
  iyMax = ny
  izMax = nz
  ierr = PipGetTwoIntegers('XMinAndMax', ixMin, ixMax)
  ierr = PipGetTwoIntegers('YMinAndMax', iyMin, iyMax)
  ierr = PipGetTwoIntegers('ZMinAndMax', izMin, izMax)
  if (ixMin < 1 .or. ixMax > nx .or. ixMax - ixMin < 2 * beadSize .or. &
      iyMin < 1 .or. iyMax > ny .or. iyMax - iyMin < 2 * beadSize .or. &
      izMin < 1 .or. izMax > nz .or. izMax - izMin < 2 * beadSize) call exitError( &
      'Coordinate min and max values are out of range or define too small a volume')
  if (maxPeaks < 3) call exitError('The -MaxNumBeads entry is negative or too small')
  limPeak = maxPeaks + 10
  allocate(indPeak(limPeak), indCorr(limPeak), peakVal(limPeak), peakPos(3, limPeak), &
      corrVal(limPeak), corrPos(3, limPeak), stat = ierr)
  call memoryError(ierr, 'arrays for peaks')
  !
  call PipDone()
  !
  ! Determine the correlation box size
  !
  nxCorr = 2 * nint(1.5 * radius + 0.5)
  nyCorr = nxCorr
  nzCorr = 2 * nint((0.5 + elongation) * radius + 0.5)
  if (yElongated) then
    nyCorr = nzCorr
    nzCorr = nxCorr
  endif
  num3CorrThreads = max(1, min(8, nint(2 * (nxCorr * nyCorr * nzCorr / 32000.)**0.45)))
  num3CorrThreads = numOMPthreads(num3CorrThreads);

  maxArray = min(float(maxArray), (nx + 4.) * (ny + 4.) * (nz + 4.) +  &
      nxCorr * nyCorr * nzCorr)
  allocate(array(maxArray), stat = ierr)
  call memoryError(ierr, 'array for image data')
  maxVol = maxArray - nxCorr * nyCorr * nzCorr
  !
  ! Given correlation dimensions, set up the overlaps and minimum sizes
  !
  nxOverlap = 2 * ((nxCorr + 1) / 2 + 1)
  nyOverlap = 2 * ((nyCorr + 1) / 2 + 1)
  nzOverlap = 2 * ((nzCorr + 1) / 2 + 1)
  minYsize = minInside + nyOverlap
  minZsize = minInside + nzOverlap
  maxZ = maxVol / (nx * ny) - 3
  maxYZ = (sqrt(9. + 4. * maxVol / nx) - 3.) / 2.
  if (iVerbose > 0) print *,maxVol, nxOverlap, minYsize, minZsize, maxZ, maxYZ
  !
  if (maxZ >= nz) then
    !
    ! the entire load will fit at once, set min's to ny and nz
    !
    minZsize = nz
    minYsize = ny
  else if (maxZ / 2 >= minZsize) then
    !
    ! X/Y planes will fit; set miny to ny and increase minZsize
    !
    minYsize = ny
    minZsize = maxZ / 2
  else if (maxYZ / 2 >= minZsize .and. maxYZ / 2 >= minYsize) then
    !
    ! X rows will fit in their entirety; increase the min sizes
    !
    minYsize = maxYZ / 2
    minZsize = maxYZ / 2
  endif
  !
  ! Get the definition of the pieces
  !
  call definePieces(izMin, izMax, nzOverlap, minZsize, 0, maxPiece, numZpieces, &
      lenPiece(1, 3), ind0(1, 3), ind1(1, 3), ierr)
  if (ierr .ne. 0) call exitError('Too many pieces for array in Z dimension')
  call definePieces(iyMin, iyMax, nyOverlap, minYsize, 0, maxPiece, numYpieces, &
      lenPiece(1, 2), ind0(1, 2), ind1(1, 2), ierr)
  if (ierr .ne. 0) call exitError('Too many pieces for array in Y dimension')
  maxXsize = maxVol / (lenPiece(1, 2) * (lenPiece(1, 3) + 3))
  call definePieces(ixMin, ixMax, nxOverlap, 0, maxXsize, maxPiece, numXpieces, &
      lenPiece(1, 1), ind0(1, 1), ind1(1, 1), ierr)
  if (ierr .ne. 0) call exitError('Too many pieces for array in X dimension')
  !
  ! Redefine the longest dimension of Y or Z
  if (ny > nz) then
    maxXsize = maxVol / (lenPiece(1, 1) * (lenPiece(1, 3) + 3))
    if (iVerbose > 0) print *,'New max for y', maxXsize
    call definePieces(iyMin, iyMax, nyOverlap, 0, maxXsize, maxPiece, &
        numYpieces, lenPiece(1, 2), ind0(1, 2), ind1(1, 2), ierr)
  else
    maxXsize = maxVol / (lenPiece(1, 1) * (lenPiece(1, 2) + 3))
    if (iVerbose > 0) print *,'New max for z', maxXsize
    call definePieces(izMin, izMax, nzOverlap, 0, maxXsize, maxPiece, &
        numZpieces, lenPiece(1, 3), ind0(1, 3), ind1(1, 3), ierr)
  endif

  if (ierr .ne. 0 .or. lenPiece(1, 1) * lenPiece(1, 2) * (lenPiece(1, 3) + 3) &
      > maxVol) call exitError('Bug in dividing volume into pieces')

  nsum = max(1, nint(0.75 * radius))
  if (iVerbose > 0) then
    print *,numXpieces, numYpieces, numZpieces, nxOverlap, nyOverlap, nzOverlap
    print *,(ind0(i, 1), ind1(i, 1), lenPiece(i, 1), i = 1, numXpieces)
    print *,(ind0(i, 2), ind1(i, 2), lenPiece(i, 2), i = 1, numYpieces)
    print *,(ind0(i, 3), ind1(i, 3), lenPiece(i, 3), i = 1, numZpieces)
    print *,'nsum = ', nsum
  endif
  numPeaks = 0
  indPlanes = lenPiece(1, 1) * lenPiece(1, 2) * lenPiece(1, 3) + 1
  !
  ! Scan for peaks in simple pixel sums
  !
  do izPiece = 1, numZpieces
    do iyPiece = 1, numYpieces
      do ixPiece = 1, numXpieces
        call loadvol(1, array, lenPiece(ixPiece, 1), lenPiece(iyPiece, 2), &
            ind0(ixPiece, 1), ind1(ixPiece, 1), ind0(iyPiece, 2), &
            ind1(iyPiece, 2), ind0(izPiece, 3), ind1(izPiece, 3))
        call getAnalysisLimits(ind0(ixPiece, 1), ind1(ixPiece, 1), nx, &
            nxOverlap, nsum, ixStart, ixEnd)
        call getAnalysisLimits(ind0(iyPiece, 2), ind1(iyPiece, 2), ny, &
            nyOverlap, nsum, iyStart, iyEnd)
        call getAnalysisLimits(ind0(izPiece, 3), ind1(izPiece, 3), nz, &
            nzOverlap, nsum, izStart, izEnd)
        call findPixelSumPeaks(array, array(indPlanes), &
            lenPiece(ixPiece, 1), lenPiece(iyPiece, 2), lenPiece(izPiece, 3), &
            ixStart, ixEnd, iyStart, iyEnd, izStart, izEnd, &
            ind0(ixPiece, 1), ind0(iyPiece, 2), ind0(izPiece, 3), nsum, &
            polarity, nxCorr, nyCorr, nzCorr, indPeak, peakVal, peakPos, &
            maxPeaks, numPeaks, peakRelMin)
      enddo
    enddo
  enddo
  !
  ! remove points that are too close to each other
  !
  print *,numPeaks, ' candidate peaks found'
  peakCorr = peakVal(indPeak(1))
  do i = 1, numPeaks
    peakVal(indPeak(i)) = peakVal(indPeak(i)) / peakCorr
    corrVal(i) = peakVal(indPeak(i))
    if (iVerbose > 1) write(*,'(f8.4,3f8.1)') peakVal(indPeak(i)), &
        (peakPos(ix, indPeak(i)), ix = 1, 3)
  enddo

  call cleanSortedList(indPeak, peakVal, peakPos, numPeaks, distMin, &
      cleanBoth, peakRelMin)
  print *,numPeaks, ' candidate peaks left after eliminating close points'
  !
  call flush(6)
  ierr = findHistogramDip(corrVal, numPeaks, minGuess, histo, limHisto, &
      0., 1., histDip, peakBelow, peakAbove, max(0, iVerbose - 2))
  blackThresh = histDip
  !
  if (ierr .ne. 0 .and. avgThresh < 0 .and. avgFallback > 0.) then
    print *,'No histogram dip found for initial peaks, using fallback averaging '// &
        'threshold'
    avgThresh = avgFallback
  endif
  if (avgThresh < 0) then
    if (ierr .ne. 0) call exitError('No histogram dip found for initial'// &
        ' peaks; enter positive -thresh to proceed')
    findCorr = histDip
    if (avgThresh < -1) then
      call findValueInList(corrVal, numPeaks, &
          histDip + (peakAbove - histDip) / 4., numLook)
    else
      call findValueInList(corrVal, numPeaks, histDip, numLook)
      numLook = max(1, nint(-avgThresh * numLook))
    endif
  else
    !
    ! Find number to average for number between 0 and 1
    if (avgThresh <= 1.) then
      call findValueInList(corrVal, numPeaks, avgThresh, numLook)
    else
      numLook = nint(avgThresh)
    endif
    if (ierr .ne. 0) blackThresh = 0
  endif
  call writePeakModel(firstFile, indPeak, peakVal, peakPos, &
      numPeaks, blackThresh, nxyz, delta, origin, radius)
  !
  ! Now loop through the pieces again looking for points and getting
  ! correlation positions
  !
  !
  maxShift = max(8, nint(radius))
  numCorrs = 0
  numSave = numPeaks
  print *,numLook, ' peaks being averaged to make reference for correlation'
  do loopCorr = 1, 2 * numPass
    !
    ! Zero the array on odd loops
    if (mod(loopCorr, 2) == 1) then
      do i = 1, nxCorr * nyCorr * nzCorr
        array(maxVol + i) = 0.
      enddo
    endif
    !
    ! Copy correlation positions to peak positions on loop 3
    if (loopCorr == 3) then
      do i = 1, numCorrs
        ix  = indCorr(i)
        indPeak(i) = ix
        peakVal(ix) = corrVal(ix)
        peakPos(1, ix) = corrPos(1, ix)
        peakPos(2, ix) = corrPos(2, ix)
        peakPos(3, ix) = corrPos(3, ix)
      enddo
      numPeaks = numCorrs
      numCorrs = 0
    endif
    if (loopCorr == 2 * numPass) numLook = numSave
    if (iVerbose > 0) print *,loopCorr, numLook
    do izPiece = 1, numZpieces
      do iyPiece = 1, numYpieces
        do ixPiece = 1, numXpieces
          loaded = numXpieces * numYpieces * numZpieces == 1
          call getAnalysisLimits(ind0(ixPiece, 1), ind1(ixPiece, 1), nx, &
              nxOverlap, nxCorr, ixStart, ixEnd)
          call getAnalysisLimits(ind0(iyPiece, 2), ind1(iyPiece, 2), ny, &
              nyOverlap, nyCorr, iyStart, iyEnd)
          call getAnalysisLimits(ind0(izPiece, 3), ind1(izPiece, 3), nz, &
              nzOverlap, nzCorr, izStart, izEnd)
          !
          ! Loop on points, for ones inside the box, get correlation
          !
          if (iVerbose > 0) print *,ixStart, ixEnd, iyStart, iyEnd, &
              izStart, izEnd
          do i = 1, numLook
            ix = indPeak(i)
            xpeak = peakPos(1, ix) - ind0(ixPiece, 1)
            ypeak = peakPos(2, ix) - ind0(iyPiece, 2)
            zpeak = peakPos(3, ix) - ind0(izPiece, 3)
            ixPeak = nint(xpeak)
            iyPeak = nint(ypeak)
            izPeak = nint(zpeak)
            if (ixPeak >= ixStart .and. ixPeak <= ixEnd .and. &
                iyPeak >= iyStart .and. iyPeak <= iyEnd .and. &
                izPeak >= izStart .and. izPeak <= izEnd) then

              ! print *,i, xpeak, ypeak, zpeak
              if (.not.loaded) then
                call loadvol(1, array, lenPiece(ixPiece, 1), &
                    lenPiece(iyPiece, 2), &
                    ind0(ixPiece, 1), ind1(ixPiece, 1), ind0(iyPiece, 2), &
                    ind1(iyPiece, 2), ind0(izPiece, 3), ind1(izPiece, 3))
                loaded = .true.
              endif

              if (mod(loopCorr, 2) == 1) then
                !
                ! On the first round from pixel sums, get centroid before
                ! adding peak in
                if (loopCorr == 1) then
                  ix0 = max(1, ixPeak - nxCorr / 2)
                  ix1 = min(lenPiece(ixPiece, 1), ix0 + nxCorr - 1)
                  iy0 = max(1, iyPeak - nyCorr / 2)
                  iy1 = min(lenPiece(iyPiece, 2), iy0 + nyCorr - 1)
                  iz0 = max(1, izPeak - nzCorr / 2)
                  iz1 = min(lenPiece(izPiece, 3), iz0 + nzCorr - 1)
                  call findPeakCenter(array, lenPiece(ixPiece, 1), &
                      lenPiece(iyPiece, 2), lenPiece(izPiece, 3), ix0, ix1, &
                      iy0, iy1, iz0, iz1, polarity, dxAdjacent, dyAdjacent, dzAdjacent)
                  xpeak = (ix1 + ix0 - 1) / 2. + dxAdjacent
                  ypeak = (iy1 + iy0 - 1) / 2. + dyAdjacent
                  zpeak = (iz1 + iz0 - 1) / 2. + dzAdjacent
                endif
                call addPeakToSum(array(maxVol + 1), nxCorr, &
                    nyCorr, nzCorr, array, lenPiece(ixPiece, 1), &
                    lenPiece(iyPiece, 2), lenPiece(izPiece, 3), xpeak, &
                    ypeak, zpeak)
              else
                call find_best_corr(array(maxVol + 1), nxCorr, nyCorr, nzCorr, array, &
                    lenPiece(ixPiece, 1), lenPiece(iyPiece, 2), lenPiece(izPiece, 3), &
                    ixStart, ixEnd, iyStart, iyEnd, izStart, izEnd, xpeak, ypeak, zpeak, &
                    maxShift, found, peakCorr, num3CorrThreads)
                if (found) then
                  call addToSortedList(indCorr, corrVal, numCorrs, &
                      maxPeaks, peakRelMin, sqrt(peakCorr), index)
                  if (index > 0) then
                    corrPos(1, index) = xpeak + ind0(ixPiece, 1) + sumXoffset
                    corrPos(2, index) = ypeak + ind0(iyPiece, 2) + sumYoffset
                    corrPos(3, index) = zpeak + ind0(izPiece, 3) + sumZoffset
                    ! write(*,'(i6,f8.4,i3,f10.5,6f7.1)') i, peakVal(ix), &
                    ! peakCorr*1.e-9, (peakPos(ix1, ix), ix1=1, 3), &
                    ! (corrPos(ix1, index), ix1=1, 3)
                  endif
                endif
              endif
            endif
          enddo
        enddo
      enddo
    enddo
    if (mod(loopCorr, 2) == 1) then
      call findPeakCenter(array(maxVol + 1), nxCorr, nyCorr, nzCorr, &
          1, nxCorr, 1, nyCorr, 1, nzCorr, &
          polarity, sumXoffset, sumYoffset, sumZoffset)
      if (iVerbose > 0) print *,'Offsets:', sumXoffset, sumYoffset, &
          sumZoffset
      ! do iz = 1, nzCorr
      ! do iy = 1, nyCorr
      ! write(*,'(2i3,(9i8))') iy, iz, (nint(array(maxVol + ix + (iy - 1) * &
      ! nxCorr + (iz - 1) * nxCorr * nyCorr)), ix =1, nxCorr)
      ! enddo
      ! enddo
      call meanzero(array(maxVol + 1), nxCorr, nxCorr, nyCorr * nzCorr)
      if (iVerbose > 0) call volWrite(17, 'volsum.st', array(maxVol + 1), &
          nxCorr, nyCorr, nzCorr)
    endif
  enddo
  !
  print *,numCorrs, ' peaks found by correlation'
  call cleanSortedList(indCorr, corrVal, corrPos, numCorrs, distMin, &
      cleanBoth, peakRelMin)
  print *,numCorrs, ' peaks left after eliminating close points'
  peakCorr = corrVal(indCorr(1))
  do i = 1, numCorrs
    corrVal(indCorr(i)) = corrVal(indCorr(i)) / peakCorr
    peakVal(i) = corrVal(indCorr(i))
    if (iVerbose > 1) write(*,'(f8.4,3f8.1)') corrVal(indCorr(i)), &
        (corrPos(ix, indCorr(i)), ix = 1, 3)
  enddo
  !
  call flush(6)
  ierr = findHistogramDip(peakVal, numCorrs, minGuess, histo, limHisto, &
      0., 1., histDip, peakBelow, peakAbove, max(0, iVerbose - 2))
  blackThresh = histDip
  !
  ! cryoposition looks for 'using fallback storage threshold'
  ! and for 'Storing' and 'peaks in model' in one line
  if (ierr .ne. 0 .and. storeThresh <= 0. .and. storeFallback > 0) then
    print *,'No dip found in histogram of correlation peaks, using fallback storage '// &
        'threshold'
    storeThresh = storeFallback
  endif
  if (storeThresh <= 0.) then
    if (ierr .ne. 0) call exitError('No dip found in histogram of '// &
        'correlation peaks; enter -store with positive value to proceed')
    call findValueInList(peakVal, numCorrs, histDip, numLook)
    print *,numLook, ' peaks are above the histogram dip'
    if (storeThresh < 0) then
      numPeaks = max(1, min(numCorrs, nint(-storeThresh * numLook)))
      print *,'Storing', numPeaks, ' peaks in model'
    else
      call avgsd(peakVal, numLook, aboveAvg, aboveSD, peakBelow)
      peakBelow = aboveAvg - 5. * aboveSD
      call findValueInList(peakVal, numCorrs, peakBelow, numPeaks)
      numPeaks = max(numLook, min(2 * numLook, numPeaks))
      write(*,'(a,i7,a,f7.4)') 'Storing an additional', numPeaks - numLook, &
          ' peaks in model down to a value of', peakBelow
    endif

  else
    call findValueInList(peakVal, numCorrs, min(1., storeThresh), numPeaks)
    if (ierr .ne. 0) blackThresh = 0
    write(*,'(a,i7,a,f7.4)') 'Storing', numPeaks, &
        ' peaks in model above threshold of', storeThresh
  endif

  call writePeakModel(modelFile, indCorr, corrVal, corrPos, &
      numPeaks, blackThresh, nxyz, delta, origin, radius)
  call exit(0)
97 call exitError('Reading tilt file')
end program findbeads3d



! DEFINEPIECES will determine how to divide an extent into overlapping
! pieces.  Inputs:
! indMin = minimum coordinate, numbered from 1
! indMax = maximum coordinate, numbered from 1
! nOverlap = overlap between pieces
! minSize = minimum size of pieces, or 0 to apply maximum size instead
! maxSize = maximum size of pieces ( if minSize is 0)
! maxPieces = maximum number of pieces allowed
! Outputs:
! nPieces = number of pieces
! lenPiece = array with lengths of pieces
! ind0, ind1 = arrays of starting, ending indices numbered from zero
! ierr = 1 if too many pieces for arrays
!
subroutine definePieces(indMin, indMax, nOverlap, minSize, maxSize, maxPieces, &
    numPieces, lenPiece, ind0, ind1, ierr)
  implicit none
  integer*4 nTotal, nOverlap, minSize, maxSize, maxPieces, numPieces, indMin, indMax
  integer*4 lenPiece(*), ind0(*), ind1(*), isize, i, ierr, irem
  !
  ! If a minimum size is defined, find the biggest number of pieces that
  ! divide into sizes bigger than the minimum
  !
  nTotal = indMax + 1 - indMin
  ierr = 1
  if (minSize > 0) then
    isize = minSize + 1
    numPieces = 0
    do while (numPieces < maxPieces .and. isize >= minSize)
      numPieces = numPieces + 1
      isize = (nTotal + (numPieces - 1) * nOverlap) / numPieces
    enddo
    if (isize >= minSize) return
    if (numPieces > 1) numPieces = numPieces - 1
  else
    !
    ! Otherwise just compute the number from the maximum size
    !
    numPieces = (nTotal - maxSize) / (maxSize - nOverlap) + 1
    if (numPieces <= 0) numPieces = 1
    isize = (nTotal + (numPieces - 1) * nOverlap + numPieces - 1) / numPieces
    if (isize > maxSize) numPieces = numPieces + 1
    if (numPieces > maxPieces) return
  endif
  !
  ! get basic size and remainder to distribute, loop on pieces
  !
  isize = (nTotal + (numPieces - 1) * nOverlap) / numPieces
  irem = mod(nTotal + (numPieces - 1) * nOverlap, numPieces)
  ind0(1) = indMin - 1
  do i = 1, numPieces
    lenPiece(i) = isize
    if (irem >= i) lenPiece(i) = lenPiece(i) + 1
    if (i > 1) ind0(i) = ind0(i - 1) + lenPiece(i - 1) - nOverlap
    ind1(i) = ind0(i) + lenPiece(i) - 1
  enddo
  ierr = 0
  return
end subroutine definePieces

! LOADVOL loads a subset of the volume from unit IUNIT, into ARRAY
! assuming dimensions of NXDIM by NYDIM, from index coordinates
! IX0, IX1, IY0, IY1, IZ0, IZ1.
!
subroutine loadvol(iunit, array, nxDim, nyDim, ix0, ix1, iy0, iy1, iz0, iz1)
  implicit none
  integer*4 nxDim, nyDim, ix0, ix1, iy0, iy1, iz0, iz1, iunit, indz, iz
  real*4 array(nxDim,nyDim,*)
  !
  ! print *,iunit, nxdim, nydim, ix0, ix1, iy0, iy1, iz0, iz1
  indz = 0
  do iz = iz0, iz1
    indz = indz + 1
    call iiuSetPosition(iunit, iz, 0)
    call irdpas(iunit, array(1, 1, indz), nxDim, nyDim, ix0, ix1, iy0, iy1,*99)
  enddo
  return
99 call exitError('Reading image file')
end subroutine loadvol


! AddToSortedList maintains a list of VALUES with an ordering in INDEX
! numVals is the number currently in the list
! maxVals is the maximum number allowed
! peakRelMin is a threshold for adding a value relative to current
! maximum
! valNew is the new value to add
! newIndex is returned with the index in VALUES at which it was placed,
! or 0 if the value is too small to fit on the list
!
subroutine addToSortedList(index, values, numVals, maxVals, peakRelMin, &
    valNew, newIndex)
  implicit none
  integer*4 index(*), numVals, maxVals, newIndex, i
  real*4 values(*), valNew, peakRelMin
  integer*4 less, more, itest, newOrder
  !
  ! quick test for full list and less than the last item on it
  ! And for value less than threshold for storing
  newIndex = 0
  if (numVals == maxVals) then
    if (values(index(maxVals)) >= valNew) return
  endif
  if (numVals > 0) then
    if (valNew < peakRelMin * values(index(1))) return
  endif
  !
  ! Handle simple cases of inserting at the front or end of the list
  !
  if (numVals == 0) then
    newOrder = 1
  else if (valNew >= values(index(1))) then
    newOrder = 1
  else if (values(index(numVals)) >= valNew) then
    newOrder = numVals + 1
  else
    !
    ! Otherwise search for position.  Set up index of one more and one
    ! less than the value.  Divide the interval in two and revise less
    ! or more until there is no longer an interval between them
    !
    more = 1
    less = numVals
    do while (less - more > 1)
      itest = (less + more) / 2
      if (values(index(itest)) == valNew) then
        more = itest
        less = itest
      else if (values(index(itest)) < valNew) then
        less = itest
      else
        more = itest
      endif
    enddo
    newOrder = less
  endif
  !
  ! If space exists, position is at end of list; otherwise it is the
  ! position occupied by the one being bumped off
  !
  if (numVals < maxVals) then
    newIndex = numVals + 1
    numVals = numVals + 1
  else
    newIndex = index(maxVals)
  endif
  values(newIndex) = valNew
  !
  ! shift index array up if necessary
  !
  do i = numVals, newOrder + 1, -1
    index(i) = index(i - 1)
  enddo
  index(newOrder) = newIndex
  ! print *,'Inserted ', valNew, ' at', newOrder, ', index', newIndex
  ! write(*,'(4(2i5,f9.0))') (i, index(i), values(index(i)), i =1, numVals)
  return
end subroutine addToSortedList

! Cleans the sorted list by eliminating ones that are too close
!
subroutine cleanSortedList(index, peakVal, peakPos, numVals, distMin, &
    cleanBoth, peakRelMin)
  implicit none
  integer*4 index(*), peakVal(*), numVals
  real*4 peakPos(3,*), distMin, peakRelMin, thresh
  integer*4 i, j, indi, indj
  real*4 xx, yy, zz, dx, dy, dz, sqrMin
  logical cleanBoth, knockout
  !
  sqrMin = distMin**2
  thresh = peakRelMin * peakVal(index(1))
  !
  ! Scan from strongest down, knocking out peaks if point is too close
  !
  do j = 1, numVals - 1
    indj = index(j)
    knockout = .false.
    if (indj > 0) then
      if (peakVal(indj) < thresh) then
        indj = 0
        knockout = .true.
      endif
    endif
    if (indj > 0) then
      xx = peakPos(1, indj)
      yy = peakPos(2, indj)
      zz = peakPos(3, indj)
      do i = j + 1, numVals
        indi = index(i)
        if (indi > 0) then
          dx = xx - peakPos(1, indi)
          if (abs(dx) < distMin) then
            dy = yy - peakPos(2, indi)
            if (abs(dy) < distMin) then
              dz = zz - peakPos(3, indi)
              if (abs(dz) < distMin .and. &
                  dx**2 + dy**2 + dz**2 < sqrMin) then
                index(i) = 0
                knockout = .true.
              endif
            endif
          endif
        endif
      enddo
    endif
    if (knockout .and. cleanBoth) index(j) = 0
  enddo
  !
  ! repack index
  !
  j = 0
  do i = 1, numVals
    if (index(i) > 0) then
      j = j + 1
      index(j) = index(i)
    endif
  enddo
  numVals = j
  return
end subroutine cleanSortedList


! Finds the limits for analyzing one chunk of the tomogram
!
subroutine getAnalysisLimits(ind0, ind1, nTotal, nOverlap, ncorr, &
    iStart, iEnd)
  implicit none
  integer*4 ind0, ind1, nTotal, nOverlap, ncorr, iStart, iEnd
  if (ind0 > 0) then
    iStart = nOverlap / 2
  else
    iStart = 2 + ncorr / 2
  endif
  if (ind1 < nTotal - 1) then
    iEnd = ind1 + 1 - ind0 - nOverlap / 2
  else
    iEnd = ind1 - ind0 - ncorr / 2
  endif
  return
end subroutine getAnalysisLimits


! findPixelSumPeaks will find peaks within a subvolume consisting of
! rolling sums over cubes NSUM on a side
! ARRAY contains the volume, dimensioned NX x nY x NZ
! PLANES is an array NX x NY x 3 for storing sums of pixels
! ixStart, ixEnd are the starting indixes in X at which to look for
! peaks, similarly for iyStart, iyEnd and izStart, izEnd.
! ixOffset, iyOffset, izOffset are the offsets of the subvolume in the
! whole volume
! POLARITY is -1 to look for dark peaks
! indPeak in an array with indexes to sorted peak values in peakVal
! peakPos will receive the X, y, Z coordinates of the peaks
! maxPeaks is the maximum number of peaks to store
! numPeaks is the number currently stored
! peakRelMin is a threshold for storing peaks relative to current max
!
subroutine findPixelSumPeaks(array, planes, nx, ny, nz, ixStart, ixEnd, &
    iyStart, iyEnd, izStart, izEnd, ixOffset, iyOffset, izOffset, &
    nsum, polarity, nxCorr, nyCorr, nzCorr, indPeak, peakVal, peakPos, &
    maxPeaks, numPeaks, peakRelMin)
  implicit none
  integer*4 nx, ny, nz, ixStart, ixEnd, iyStart, iyEnd, izStart, izEnd
  real*4 array(nx, ny, nz), planes(nx, ny, 3), peakVal(*), peakPos(3, *)
  integer*4 indPeak(*), nxCorr, nyCorr, nzCorr
  integer*4 ixOffset, iyOffset, izOffset, nsum, maxPeaks, numPeaks
  real*4 polarity, cen, posShift, sum, sumScale, peakRelMin
  integer*4 ix, iy, iz, jx, jy, jz, nextPlane, numBefore, numAfter, ip1, ip2
  integer*4 ip3, index, ix0, ix1, iy0, iy1, iz0, iz1

  !
  nextPlane = 1
  numBefore = (nsum - 1) / 2
  numAfter = nsum / 2
  posShift = 0.5 * (numAfter - numBefore - 1)
  sumScale = polarity / nsum**3
  !
  ! Loop on the Z planes with extended limits
  !
  do iz = izStart - 1, izEnd + 1
    !
    ! first compute the sums on the next plane
    !
    do iy = iyStart - 1, iyEnd + 1
      do ix = ixStart - 1, ixEnd + 1
        sum = 0.
        do jz = iz - numBefore, iz + numAfter
          do jy = iy - numBefore, iy + numAfter
            do jx = ix - numBefore, ix + numAfter
              sum = sum + array(jx, jy, jz)
            enddo
          enddo
        enddo
        planes(ix, iy, nextPlane) = sum * sumScale
      enddo
    enddo
    !
    ! If at least 3 planes have been computed, now search for peaks
    !
    if (iz > izStart) then
      ip1 = mod(nextPlane, 3) + 1
      ip2 = mod(nextPlane + 1, 3) + 1
      ip3 = nextPlane
      do iy = iyStart, iyEnd
        do ix = ixStart, ixEnd
          cen = planes(ix, iy, ip2)
          !
          ! Test for whether greater than all diagonal and ones ahead
          ! in X, Y, or Z, and >= ones behind in X, Y, or Z
          !
          if (cen >= planes(ix - 1, iy, ip2) .and. cen > planes(ix + 1, iy, ip2) &
              .and. cen >= planes(ix, iy - 1, ip2) .and. &
              cen > planes(ix, iy + 1, ip2) .and. cen >= planes(ix, iy, ip1) &
              .and. cen > planes(ix, iy, ip3) .and. &
              cen > planes(ix, iy + 1, ip1) .and. cen >= planes(ix, iy - 1, ip1) &
              .and. cen >= planes(ix - 1, iy, ip1) .and. &
              cen >= planes(ix + 1, iy, ip1) .and. cen >= planes(ix + 1, iy, ip3) &
              .and. cen >= planes(ix - 1, iy, ip3) .and. &
              cen > planes(ix, iy + 1, ip3) .and. cen >= planes(ix, iy - 1, ip3) &
              .and. cen >= planes(ix - 1, iy - 1, ip2) .and. &
              cen >= planes(ix + 1, iy - 1, ip2) .and. &
              cen >= planes(ix + 1, iy + 1, ip2) .and. &
              cen >= planes(ix - 1, iy + 1, ip2) .and. &
              cen >= planes(ix - 1, iy - 1, ip1) .and. &
              cen >= planes(ix + 1, iy - 1, ip1) .and. &
              cen >= planes(ix + 1, iy + 1, ip1) .and. &
              cen >= planes(ix - 1, iy + 1, ip1) .and. &
              cen >= planes(ix - 1, iy - 1, ip3) .and. &
              cen >= planes(ix + 1, iy - 1, ip3) .and. &
              cen >= planes(ix + 1, iy + 1, ip3) .and. &
              cen >= planes(ix - 1, iy + 1, ip3)) then
            !
            ! Got a peak; try to add it to the list
            !
            ix0 = ix - nxCorr / 2
            ix1 = ix0 + nxCorr - 1
            iy0 = iy - nyCorr / 2
            iy1 = iy0 + nyCorr - 1
            iz0 = iz - nzCorr / 2 - 1
            iz1 = iz0 + nzCorr - 1
            ! print *,cen, ix0, ix1, iy0, iy1, iz0, iz1, nx, ny, nz
            if (ix0 >= 1 .and. ix1 <= nx .and. iy0 >= 1 .and. &
                iy1 <= ny .and. iz0 >= 1 .and. iz1 <= nz) then
              call integratePeak(array, nx, ny, nz, ix0, ix1, iy0, iy1, &
                  iz0, iz1, polarity, cen)
              ! print *,cen
              if (cen > 0) then
                call addToSortedList(indPeak, peakVal, numPeaks, maxPeaks, &
                    peakRelMin, sqrt(cen), index)
                if (index > 0) then
                  peakPos(1, index) = ix + ixOffset + posShift
                  peakPos(2, index) = iy + iyOffset + posShift
                  peakPos(3, index) = iz + izOffset + posShift - 1
                endif
              endif
            endif
          endif
        enddo
      enddo

    endif
    nextPlane = mod(nextPlane, 3) + 1
  enddo
  return
end subroutine findPixelSumPeaks

! addPeakToSum adds one peak located at dxadj, dyadj, dzadj in the
! volume brray (size nxb, nyb, nzb) to the sum in array (size nxa, nya,
! nza)
subroutine addPeakToSum(array, nxa, nya, nza, brray, nxb, nyb, nzb, &
    dxAdjacent, dyAdjacent, dzAdjacent)
  implicit none
  integer*4 nxa, nxb, nya, nyb, nza, nzb
  real*4 array(nxa,nya,nza), brray(nxb,nyb,nzb)
  real*4 dxAdjacent, dyAdjacent, dzAdjacent, dx, dy, dz, d11, d12, d21, d22
  integer*4 ix, iy, iz, idx, idy, idz, ixp, iyp, izp, ixpP1, iypP1, izpP1
  idx = int(dxAdjacent) - nxa / 2
  dx = dxAdjacent - int(dxAdjacent)
  idy = int(dyAdjacent) - nya / 2
  dy = dyAdjacent - int(dyAdjacent)
  idz = int(dzAdjacent) - nza / 2
  dz = dzAdjacent - int(dzAdjacent)
  d11 = (1. - dx) * (1. - dy)
  d12 = (1. - dx) * dy
  d21 = dx * (1. - dy)
  d22 = dx * dy
  do iz = 1, nza
    izp = max(1, iz + idz)
    izpP1 = min(izp + 1, nzb)
    do iy = 1, nya
      iyp = max(1, iy + idy)
      iypP1 = min(iyp + 1, nyb)
      do ix = 1, nxa
        ixp = max(1, ix + idx)
        ixpP1 = min(ixp + 1, nxb)
        array(ix, iy, iz) = array(ix, iy, iz) + &
            (1. - dz) * (d11 * brray(ixp, iyp, izp) &
            + d12 * brray(ixp, iypP1, izp) &
            + d21 * brray(ixpP1, iyp, izp) &
            + d22 * brray(ixpP1, iypP1, izp)) + &
            dz * (d11 * brray(ixp, iyp, izpP1) &
            + d12 * brray(ixp, iypP1, izpP1) &
            + d21 * brray(ixpP1, iyp, izpP1) &
            + d22 * brray(ixpP1, iypP1, izpP1))
      enddo
    enddo
  enddo
  return
end subroutine addPeakToSum


! Finds the centroid of the pixels in box above the background at the
! edges of the box
!
subroutine findPeakCenter(array, nxa, nya, nza, ix0, ix1, iy0, iy1, iz0, iz1, &
    polarity, dxAdjacent, dyAdjacent, dzAdjacent)
  implicit none
  integer*4 nxa, nya, nza, ix0, ix1, iy0, iy1, iz0, iz1
  real*4 array(nxa,nya,nza)
  real*4 polarity, dxAdjacent, dyAdjacent, dzAdjacent, edge, xsum, ysum, zsum, wsum, diff
  integer*4 ix, iy, iz
  !
  call peakEdgeMean(array, nxa, nya, nza, ix0, ix1, iy0, iy1, iz0, iz1, edge)
  !
  ! Get weighted sum of pixel indexes
  xsum = 0.
  ysum = 0.
  zsum = 0.
  wsum = 0.
  do iz = iz0 + 1, iz1 - 1
    do iy = iy0 + 1, iy1 - 1
      do ix = ix0 + 1, ix1 - 1
        diff = polarity * (array(ix, iy, iz) - edge)
        if (diff > 0.) then
          wsum = wsum + diff
          xsum = xsum + diff * (ix - ix0)
          ysum = ysum + diff * (iy - iy0)
          zsum = zsum + diff * (iz - iz0)
        endif
      enddo
    enddo
  enddo
  !
  ! Get offset from center
  dxAdjacent = xsum / wsum - (ix1 - ix0) / 2.
  dyAdjacent = ysum / wsum - (iy1 - iy0) / 2.
  dzAdjacent = zsum / wsum - (iz1 - iz0) / 2.
  return
end subroutine findPeakCenter

! Finds the integral of a peak relative to background at the edges
!
subroutine integratePeak(array, nx, ny, nz, ix0, ix1, iy0, iy1, iz0, iz1, &
    polarity, peak)
  implicit none
  integer*4 nx, ny, nz, ix0, ix1, iy0, iy1, iz0, iz1, ix, iy, iz
  real*4 array(nx,ny,nz), polarity, peak, edge
  !
  call peakEdgeMean(array, nx, ny, nz, ix0, ix1, iy0, iy1, iz0, iz1, edge)
  peak = 0.
  do iz = iz0 + 1, iz1 - 1
    do iy = iy0 + 1, iy1 - 1
      do ix = ix0 + 1, ix1 - 1
        peak = peak + array(ix, iy, iz) - edge
      enddo
    enddo
  enddo
  peak = peak * polarity
  return
end subroutine integratePeak


! Finds the mean along the walls of a box
!
subroutine peakEdgeMean(array, nx, ny, nz, ix0, ix1, iy0, iy1, iz0, iz1, &
    edge)
  implicit none
  integer*4 nx, ny, nz, ix0, ix1, iy0, iy1, iz0, iz1, ix, iy, iz
  real*4 array(nx,ny,nz), edgeSum, edge
  edgeSum = 0
  do iy = iy0, iy1
    do ix = ix0, ix1
      edgeSum = edgeSum + array(ix, iy, iz0) + array(ix, iy, iz1)
    enddo
  enddo
  do iz = iz0 + 1, iz1 - 1
    do ix = ix0, ix1
      edgeSum = edgeSum + array(ix, iy0, iz) + array(ix, iy1, iz)
    enddo
  enddo
  do iz = iz0 + 1, iz1 - 1
    do iy = iy0 + 1, iy1 - 1
      edgeSum = edgeSum + array(ix0, iy, iz) + array(ix1, iy, iz)
    enddo
  enddo
  edge = edgeSum / (2 * ((ix1 + 1 - ix0) * (iy1 + 1 - iy0) + &
      (ix1 + 1 - ix0) * (iz1 - 1 - iz0) + (iy1 - 1 - iy0) * (iz1 - 1 - iz0)))
  return
end subroutine peakEdgeMean

! FIND_BEST_CORR will search for the location in the volume BRRAY
! with the highest correlation to the
! volume in ARRAY.  ARRAY is dimensioned to NXA by NYA
! by NZA; BRRAY is dimensioned to NXB by NYB by NZB.
! IXSTART, IXEND, IYSTART, IYEND, IZSTART, IZEND are the limits for
! valid center locations (array coordinates) in BRRAY.
! DXADJ, DYADJ, DZADJ
! should contain starting location upon entry, and will return the
! location with the best correlation.  MAXSHIFT specifies the maximum
! shift that is allowed.  FOUND is returned as TRUE if a correlation
! peak is found within the maximum shift.  The maximum cross-product is
! returned in PEAKCORR
!
subroutine find_best_corr(array, nxa, nya, nza, brray, nxb, nyb, nzb, &
    ixStart, ixEnd, iyStart, iyEnd, izStart, izEnd, dxAdjacent, dyAdjacent, &
    dzAdjacent, maxShift, found, peakCorr, numThreads)
  implicit none

  integer*4 nxa, nxb, nya, nyb, nza, nzb, numThreads
  real*4 array(nxa,nya,nza), brray(nxb,nyb,nzb), peakCorr
  real*8 corrs(-1:1,-1:1,-1:1), corrTmp(-2:2,-2:2,-2:2)
  real*8 corrMax
  logical done(-1:1, -1:1, -1:1), doneTmp(-2:2, -2:2, -2:2)
  integer*4 idySequence(9) /0, -1, 1, 0, 0, -1, 1, -1, 1/
  integer*4 idzSequence(9) /0, 0, 0, 1, -1, -1, -1, 1, 1/
  logical found
  integer*4 maxShift, ixStart, ixEnd, iyStart, iyEnd, izStart, izEnd

  real*4 dxAdjacent, dyAdjacent, dzAdjacent
  !
  integer*4 idxGlobal, idyGlobal, idzGlobal, ix, iy, iz, indSequence, idy, idz
  integer*4 idyCorr, idzCorr, idxCorr, minSeq
  integer*4 indMax, idxMax, idyMax, idzMax
  real*4 cx, y1, y2, y3, cy, cz
  real*8 parabolicFitPosition
  !
  ! Minimum # of rows to do in sequence before shifting center
  minSeq = 5
  !
  ! get global displacement of b
  !
  ! print *,'findbest', nxa, nya, nza, &
  ! nxb, nyb, nzb, ixStart, ixEnd, iyStart, iyEnd, izStart, izEnd, &
  ! dxadj, dyadj, dzadj, maxshift
  idxGlobal = nint(dxAdjacent)
  idyGlobal = nint(dyAdjacent)
  idzGlobal = nint(dzAdjacent)
  ! print *,'starting', idxglb, idyglb, idzglb
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
      idxCorr = idxGlobal - nxa / 2
      idyCorr = idyGlobal + idy - nya / 2
      idzCorr = idzGlobal + idz - nza / 2
      call threecorrs(array, nxa, nya, brray, nxb, nyb, 0, nxa - 1, &
          0, nya - 1, 0, nza - 1, idxCorr, idyCorr, idzCorr, &
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
      ! print *,'New max', corrmax
      corrMax = corrs(indMax, idy, idz)
      idxMax = indMax
      idyMax = idy
      idzMax = idz
    endif

    if (indSequence >= minSeq .and. &
        (idxMax .ne. 0 .or. idyMax .ne. 0 .or. idzMax .ne. 0)) then
      !
      ! if there is a new maximum, after a minimum number of rows has
      ! been done, shift the done flags and the existing
      ! correlations, and reset the sequence
      !
      idxGlobal = idxGlobal + idxMax
      idyGlobal = idyGlobal + idyMax
      idzGlobal = idzGlobal + idzMax
      ! print *,'moving by', idxmax, idymax, idzmax
      !
      ! but if beyond the limit, return failure
      !
      if (max(abs(idxGlobal - dxAdjacent), abs(idyGlobal - dyAdjacent),  &
          abs(idzGlobal - dzAdjacent)) > maxShift .or.  &
          idxGlobal < ixStart .or. idxGlobal > ixEnd .or. &
          idyGlobal < iyStart .or. idyGlobal > iyEnd .or. &
          idzGlobal < izStart .or. idzGlobal > izEnd) then
        found = .false.
        ! print *,'outside'
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
      indSequence = 0
      idxMax = 0
      idyMax = 0
      idzMax = 0
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
  ! print *,'X:', y1, y2, y3, cx
  y1 = corrs(0, -1, 0)
  y3 = corrs(0, 1, 0)
  cy = parabolicFitPosition(y1, y2, y3)
  ! print *,'Y:', y1, y2, y3, cy
  y1 = corrs(0, 0, -1)
  y3 = corrs(0, 0, 1)
  cz = parabolicFitPosition(y1, y2, y3)
  ! print *,'Z:', y1, y2, y3, cz
  !
  dxAdjacent = idxGlobal + cx
  dyAdjacent = idyGlobal + cy
  dzAdjacent = idzGlobal + cz
  peakCorr = y2
  found = .true.
  ! print *,'returning a peak', peakCorr, dxadj, dyadj, dzadj
  return
end subroutine find_best_corr


! THREECORRS computes three correlations between volumes in ARRAY
! and BRRAY.  The volume in ARRAY is dimensioned to NXA by NYA, that
! in BRRAY is dimensioned to NXB by NYB.  The starting and ending
! index coordinates (numbered from 0) in ARRAY over which the
! correlations are to be computed are IX0, IX1, IY0, IY1, IZ0, IZ1.
! The shift between coordinates in ARRAY and coordinates in B is given
! by IDX, IDY, IDZ.  The three correlations are returned in CORR1 (for
! IDX-1), CORR2 (for IDX), and CORR3 (for IDX+1) .
!
subroutine threecorrs(array, nxa, nya, brray, nxb, nyb, ix0, ix1, &
    iy0, iy1, iz0, iz1, idx, idy, idz, corr1, corr2, corr3, numThreads)
  implicit none
  real*4 array(*), brray(*)
  ! real*4 first, prev
  real*8 sum1, sum2, sum3, corr1, corr2, corr3
  integer*4 ix0, ix1, iy0, iy1, iz0, iz1, idx, idy, idz, nxa, nxb, nya, nyb, numThreads
  integer*4 iz, izb, iy, iyb, indBaseA, indDelB, ix, ixb, nsum
  sum1 = 0.
  sum2 = 0.
  sum3 = 0.

  !$OMP PARALLEL DO NUM_THREADS(numThreads) REDUCTION (+ : sum1, sum2, sum3) &
  !$OMP& SHARED(iz0, iz1, idz, iy0, iy1, idy, nxA, nyA, nxB, nyB, ix0, ix1, array, brray)&
  !$OMP& PRIVATE(iz, izB, iy, iyb, indBaseA, indDelB, ix, ixB)
  do iz = iz0, iz1
    izb = iz + idz
    do iy = iy0, iy1
      iyb = iy + idy
      indBaseA = 1 + iy * nxa + iz * nxa * nya
      indDelB = 1 + iyb * nxb + izb * nxb * nyb + idx - indBaseA

      ! 11/28/14: deleted all the stuff for correlation coefficients, which are not
      ! appropriate for featureless particles
      do ix = indBaseA + ix0, indBaseA + ix1
        ixb = ix + indDelB
        sum1 = sum1 + array(ix) * brray(ixb - 1)
        sum2 = sum2 + array(ix) * brray(ixb)
        sum3 = sum3 + array(ix) * brray(ixb + 1)
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
end subroutine threecorrs

subroutine findValueInList(peakVal, numPeaks, findVal, index)
  implicit none
  real*4 peakVal(*), findVal
  integer*4 numPeaks, index, i
  i = 1
  do while (i <= numPeaks .and. peakVal(i) >= findVal)
    i = i + 1
  enddo
  index = i - 1
  return
end subroutine findValueInList

! Write a models with the given number of peaks to the file and set the
! threshold if any
!
subroutine writePeakModel(firstFile, indPeak, peakVal, peakPos, &
    numPeaks, blackThresh, nxyz, delta, origin, radius)
  implicit none
  include 'model.inc90'
  character*(*) firstFile
  integer*4 indPeak(*), numPeaks, iobj, i, j, ierr, nxyz(3)
  real*4 peakVal(*), peakPos(3,*), delta(*), origin(*), radius, blackThresh
  real*4  peak, peakMin, peakMax
  integer*4 ind(3) /1, 2, 3/
  integer*4 putContValue, putImodFlag, putValBlackWhite, putImageRef
  integer*4 putImodMaxes

  if (firstFile == ' ') return

  n_point = numPeaks
  iobj = 0
  max_mod_obj = numPeaks
  call newImod()
  peakMin = peakVal(indPeak(numPeaks))
  peakMax = peakVal(indPeak(1))
  ierr = putImodFlag(1, 2)
  ierr = putImodFlag(1, 7)
  call putScatSize(1, max(1, 1 + ceiling(radius)))
  ierr = putImodMaxes(nxyz(1), nxyz(2), nxyz(3))
  do i = 1, n_point
    peak = peakVal(indPeak(i))
    !
    ! set up object if it is time for next one
    !
    iobj = i
    obj_color(2, iobj) = 255
    npt_in_obj(iobj) = 1
    ibase_obj(iobj) = i - 1

    do j = 1, 3
      p_coord(j, i) = peakPos(ind(j), indPeak(i))
    enddo
    p_coord(3, i) = p_coord(3, i) - 0.5
    ierr = putContValue(1, iobj, peak)
    object(i) = i
  enddo
  iobj = 0
  if (blackThresh > 0) iobj = max(0, min(255, &
      nint(255. * (blackThresh - peakMin) / peakMax - peakMin)))
  ierr = putValBlackWhite(1, iobj, 255)
  ierr = putImageRef(delta, origin)
  call putImodObjName(1, '3D bead positions')
  call scale_model(1)
  call write_wmod(firstFile)
  return
end subroutine writePeakModel


! Routine to write a volume in a contiguous array
!
subroutine volWrite(iunit, fileOut, array, nx, ny, nz)
  implicit none
  character*(*) fileOut
  integer*4 iunit, nx, ny, nz, iz
  real*4 array(nx,ny,nz), tmin, tmax, tmean, dmin, dmax, dmean, title(20)
  integer*4 nxyz(3)
  real*4 cell(6) /0., 0., 1., 90., 90., 90./
  character*80 titlech
  !
  call imopen(iunit, fileOut, 'NEW')
  write(titlech, 3000)
3000 format('DEBUG VOLUME')
  nxyz(1) = nx
  nxyz(2) = ny
  nxyz(3) = nz
  cell(1) = nx
  cell(2) = ny
  cell(3) = nz
  call iiuCreateHeader(iunit, nxyz, nxyz, 2, title, 0)
  call iiuAltCell(iunit, cell)
  dmean = 0.
  tmin = 1.e30
  tmax = -1.e30
  do iz = 1, nz
    call iclden(array(1, 1, iz), nx, ny, 1, nx, 1, ny, tmin, tmax, tmean)
    call iiuWriteSection(iunit, array(1, 1, iz))
    dmin = min(tmin, dmin)
    dmax = max(tmax, dmax)
    dmean = dmean + tmean
  enddo
  call iiuWriteHeaderStr(iunit, titlech, 1, dmin, dmax, dmean / nz)
  call iiuClose(iunit)
  return
end subroutine volWrite
