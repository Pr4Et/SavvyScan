! * * * * * * * TAPEROUTVOL * * * * * * * *
!
! TAPEROUTVOL will cut a subset out of an image volume, pad it into a
! larger volume, and taper the intensity down to the mean value of the
! volume over the extent of the padding region, i.e., from the edge of
! the actual excised pixels to the edge of the new volume.  None of
! the original excised pixels are attenuated by this method.
!
! See man page for details.
!
! David Mastronarde, 3/1/01
!
! $Id$
!
program taperoutvol
  implicit none
  !
  integer*4 nxyz(3), mxyz(3), nxyzOut(3), mxyzOut(3), nxIn, nyIn, nzIn
  real*4 title(20), cellOut(6), tilt(3), origTilt(3)
  real*4, allocatable :: array(:), padWork(:)
  !
  character*320 inputFile, outputFile
  character*9 dateStr
  character*8 timeStr
  character*80 titleCh
  !
  equivalence (nxIn, nxyz(1)), (nyIn, nxyz(2)), (nzIn, nxyz(3))
  !
  integer*4 ixLow, iyLow, izLow, ixHigh, izHigh, iyHigh, nxBox, nyBox, nzBox
  integer*4 nxOut, nyOut, nzOut, izStart, izEnd, mode, kti
  integer*4 iz, izRead, i, ierr
  real*4 dmin2, dmax2, dmean2, dmin, dmax, dmean, tmpMin, tmpMax, sumMean
  real*4 atten, base, tmpMean, originX, originY, originZ
  logical*4 pipInput, noFFT, noisePad
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetTwoIntegers, PipGetThreeIntegers
  integer*4 PipGetInOutFile, PipGetLogical
  !
  ! fallbacks from ../../manpages/autodoc2man -3 2  taperoutvol
  !
  integer numOptions
  parameter (numOptions = 10)
  character*(40 * numOptions) options(1)
  options(1) = &
      'input:InputFile:FN:@output:OutputFile:FN:@xminmax:XMinAndMax:IP:@'// &
      'yminmax:YMinAndMax:IP:@zminmax:ZMinAndMax:IP:@taper:TaperPadsInXYZ:IT:@'// &
      'nofft:NoFFTSizes:B:@noise:NoisePadding:B:@param:ParameterFile:PF:@help:usage:B:'
  !
  noFFT = .false.
  noisePad = .false.
  !
  ! Pip startup: set error, parse options, check help, set flag if used
  !
  call PipReadOrParseOptions(options, numOptions, 'taperoutvol', &
      'ERROR: TAPEROUTVOL - ', .true., 2, 1, 1, numOptArg, &
      numNonOptArg)
  pipInput = numOptArg + numNonOptArg > 0

  if (PipGetInOutFile('InputFile', 1, 'Name of image input file', inputFile) &
      .ne. 0) call exitError('No input file specified')
  !
  call imopen(1, inputFile, 'RO')
  call irdhdr(1, nxyz, mxyz, mode, dmin2, dmax2, dmean2)
  !
  if (PipGetInOutFile('OutputFile', 2, 'Name of output file', outputFile) &
      .ne. 0) call exitError('No output file specified')
  !
  if (pipInput) ierr = PipGetLogical('NoFFTSizes', noFFT)
  if (pipInput) ierr = PipGetLogical('NoisePadding', noisePad)
  call taperPrep(pipInput, noFFT, nxyz, ixLow, ixHigh, iyLow, iyHigh, izLow, &
      izHigh, nxBox, nyBox, nzBox, nxOut, nyOut, nzOut, nxyzOut, mxyzOut, cellOut, &
      originX, originY, originZ)
  !
  iz = 1
  if (noisePad) iz = 2 * max(nxBox, nyBox) + (nxOut - nxBox) + (nyOut - nyBox)
  allocate(array(nxOut * nyOut), padWork(iz), stat = ierr)
  if (ierr .ne. 0) call exitError('Allocating array for image plane')
  !
  call imopen(2, outputFile, 'NEW')
  call time(timeStr)
  call b3ddate(dateStr)
  write(titleCh, 301) dateStr, timeStr
  read(titleCh, '(20a4)') (title(kti), kti = 1, 20)
301 format('TAPEROUTVOL: Taper outside of excised volume',t57,a9,2x,a8)
  call iiuCreateHeader(2, nxyzOut, mxyzOut, mode, title, 0)
  call iiuTransLabels(2, 1)
  call iiuAltCell(2, cellOut)
  call iiuAltOrigin(2, originX, originY, originZ)
  call iiuRetTilt(1, tilt)
  call iiuRetTiltOrig(1, origTilt)
  call iiuAltTilt(2, tilt)
  call iiuAltTiltOrig(2, origTilt)
  dmin = 1.e30
  dmax = -1.e30
  sumMean = 0.
  !
  izStart = izLow - (nzOut - nzBox) / 2
  izEnd = izStart + nzOut - 1
  do iz = izStart, izEnd
    izRead = max(izLow, min(izHigh, iz))
    call iiuSetPosition(1, izRead, 0)
    call irdpas(1, array, nxBox, nyBox, ixLow, ixHigh, iyLow, iyHigh, *99)
    if (noisePad) then
      call sliceNoiseTaperPad(array, nxBox, nyBox, array, nxOut, nxOut, nyOut, &
          max(20, min(120, max(nxBox, nyBox) / 50)), 4, padWork)
    else
      call taperoutpad(array, nxBox, nyBox, array, nxOut, nxOut, nyOut, 1, dmean2)
    endif
    if (iz < izLow .or. iz > izHigh) then
      if (iz < izLow) then
        atten = float(iz - izStart) / (izLow - izStart)
      else
        atten = float(izEnd - iz) / (izEnd - izHigh)
      endif
      base = (1. -atten) * dmean2
      do i = 1, nxOut * nyOut
        array(i) = base + atten * array(i)
      enddo
    endif
    call iclden(array, nxOut, nyOut, 1, nxOut, 1, nyOut, tmpMin, tmpMax, tmpMean)
    dmin = min(dmin, tmpMin)
    dmax = max(dmax, tmpMax)
    sumMean = sumMean + tmpMean
    call iiuWriteSection(2, array)
  enddo
  dmean = sumMean / nzOut
  call iiuWriteHeader(2, title, 1, dmin, dmax, dmean)
  call iiuClose(2)
  call exit(0)
99 call exitError('Reading file')
end program taperoutvol
