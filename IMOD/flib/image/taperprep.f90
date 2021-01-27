! $Id$
!
! TAPERPREP does common work of getting subvolume parameters and
! computing the real image box size for taperoutvol and combinefft
!
! PIPINPUT indicates whether PIP input is being used
! noFFT is a flag that sizes are not to be increased for FFTs
! NXYZ(3) has the input file sizes
! IXLOW, IXHIGH, IYLOW, IYHIGH, IZLOW, IZHIGH are returned with the index
! coordinates of the subvolume in the file
! NXBOX, NYBOX, NZBOX are returned with the size of the box being read in
! NXOUT, NYOUT, NZOUT are returned with the size of the padded volume
! NXYZOUT, MXYZOUT are filled with the values of nxOUT, nyOUT, nzOUT
! CELLOUT is returned with the correct cell values
! ORIGINX, ORIGINY, ORIGINZ are returned with the subvolume origin values
!
subroutine taperPrep(pipInput, noFFT, nxyz, ixLow, ixHigh, iyLow, iyHigh, izLow, &
    izHigh, nxBox, nyBox, nzBox, nxOut, nyOut, nzOut, nxyzOut, mxyzOut, cellOut, &
    originX, originY, originZ)
  implicit none
  logical pipInput, noFFT
  integer*4 nxyz(3), ixLow, ixHigh, iyLow, iyHigh, izLow, izHigh, nxBox, nyBox, nzBox
  integer*4 nxOut, nyOut, nzOut, nxyzOut(3), mxyzOut(3)
  real*4 cellOut(6), originX, originY, originZ
  integer*4 numPadX, mumPadY, numPadZ, ierr
  real*4 delta(3)
  integer*4 niceframe, PipGetTwoIntegers, PipGetThreeIntegers, niceFFTlimit

  ixLow = 0
  iyLow = 0
  izLow = 0
  ixHigh = nxyz(1) - 1
  iyHigh = nxyz(2) - 1
  izHigh = nxyz(3) - 1
  numPadX = 0
  mumPadY = 0
  numPadZ = 0
  if (pipInput) then
    ierr = PipGetTwoIntegers('XMinAndMax', ixLow, ixHigh)
    ierr = PipGetTwoIntegers('YMinAndMax', iyLow, iyHigh)
    ierr = PipGetTwoIntegers('ZMinAndMax', izLow, izHigh)
    ierr = PipGetThreeIntegers('TaperPadsInXYZ', numPadX, mumPadY, numPadZ)
  else
    write(*,'(1x,a,/,a,$)') 'Starting and ending X, then Y, then '// &
        'Z index coordinates to extract', ' (/ for whole volume): '
    read(5,*) ixLow, ixHigh, iyLow, iyHigh, izLow, izHigh
    write(*,'(1x,a,$)') 'Width of pad/taper borders in X, Y, and Z: '
    read(5,*) numPadX, mumPadY, numPadZ
  endif
  !
  if (ixLow < 0 .or. ixHigh >= nxyz(1) .or. iyLow < 0 .or. iyHigh >= nxyz(2) &
      .or. izLow < 0 .or. izHigh >= nxyz(3)) call exitError( &
      'Block not all inside volume')
  nxBox = ixHigh + 1 - ixLow
  nyBox = iyHigh + 1 - iyLow
  nzBox = izHigh + 1 - izLow
  !
  if (noFFT) then
    nxOut = nxBox + 2 * numPadX
    nyOut = nyBox + 2 * mumPadY
    nzOut = nzBox + 2 * numPadZ
  else
    nxOut = niceframe(2 * ((nxBox + 1) / 2 + numPadX), 2, niceFFTlimit())
    nyOut = niceframe(2 * ((nyBox + 1) / 2 + mumPadY), 2, niceFFTlimit())
    nzOut = nzBox
    if (nzOut > 1 .or. numPadZ > 0) &
        nzOut = niceframe(2 * ((nzBox + 1) / 2 + numPadZ), 2, niceFFTlimit())
  endif

  call iiuRetDelta(1, delta)
  call iiuRetOrigin(1, originX, originY, originZ)
  mxyzOut(1) = nxOut
  mxyzOut(2) = nyOut
  mxyzOut(3) = nzOut
  nxyzOut = mxyzOut
  cellOut(1:3) = mxyzOut * delta
  cellOut(4:6) = 90.
  originX = originX - delta(1) * (ixLow - (nxOut - nxBox) / 2)
  originY = originY - delta(2) * (iyLow - (nyOut - nyBox) / 2)
  originZ = originZ - delta(3) * (izLow - (nzOut - nzBox) / 2)
  return
end subroutine taperPrep
