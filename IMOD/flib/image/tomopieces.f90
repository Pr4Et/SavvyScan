! * * * * * TOMOPIECES * * * * * *
!
! TOMOPIECES figures out how to chop up a tomogram into pieces so that
! the Fourier transforms of each piece can be done in memory.
!
! See manual page for more details.
!
! $Id$
!
program tomopieces
  implicit none
  integer LIMPIECES
  parameter (LIMPIECES = 5000)
  character*80 outputFile
  integer*4 nxyz(3)
  real*4 megaVox/80/                        !MAXIMUM MEGAVOXELS
  integer*4 izBackStart(LIMPIECES), izBackEnd(LIMPIECES), izOutStart(LIMPIECES), &
      izOutEnd(LIMPIECES)
  integer*4 iyBackStart(LIMPIECES), iyBackEnd(LIMPIECES), itOutStart(LIMPIECES), &
      iyOutEnd(LIMPIECES)
  integer*4 ixBackStart(LIMPIECES), ixBackEnd(LIMPIECES), ixOutStart(LIMPIECES), &
      ixOutEnd(LIMPIECES)
  logical*4 tooBig, noFFT
  integer*4 minOverlap, nxPad, nyPad, nzPad, maxPiecesX, maxPiecesY, numZpieces
  integer*4 maxPiecesZ, maxLayerPieces, ierr, nx, ny, nz, numXpieces, numYpieces
  integer*4 nxExtra, nyExtra, nzExtra, nxOut, nyOut, nzOut, ix, iy, iz, i
  real*4 perim, perimMin, piecePerim, piecePerimMin
  integer*4 numXpieceMin, numYpieceMin, numZpieceMin, numOut, limitYpieces
  integer*4 padNiceIfFFT

  logical pipinput
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetFloat, PipGetInteger, PipGetLogical
  integer*4 PipGetInOutFile
  !
  ! fallbacks from ../../manpages/autodoc2man -3 2  tomopieces
  !
  integer numOptions
  parameter (numOptions = 11)
  character*(40 * numOptions) options(1)
  options(1) = &
      'tomogram:TomogramOrSizeXYZ:CH:@megavox:MegaVoxels:F:@'// &
      'xpad:XPadding:I:@ypad:YPadding:I:@zpad:ZPadding:I:@'// &
      'xmaxpiece:XMaximumPieces:I:@ymaxpiece:YMaximumPieces:I:@'// &
      'zmaxpiece:ZMaximumPieces:I:@'// &
      'minoverlap:MinimumOverlap:I:@param:ParameterFile:PF:@'// &
      'help:usage:B:'
  !
  minOverlap = 4                              !MINIMUM OVERLAP TO OUTPUT
  nxPad = 8                                   !PADDING / TAPER EXTENT IN X
  nzPad = 8                                   !PADDING / TAPER EXTENT IN Z
  nyPad = 4                                   !PADDING IN Y
  maxPiecesX = 0                              !MAX PIECES IN X
  maxPiecesY = 0                              !MAX PIECES IN X
  maxPiecesZ = -1                             !MAX PIECES IN X
  maxLayerPieces = 100
  noFFT = .false.
  !
  ! Pip startup: set error, parse options, do help output, get the
  ! one obligatory argument to get size
  ! But turn off the entry printing first!
  call PipEnableEntryOutput(0)
  call PipReadOrParseOptions(options, numOptions, 'tomopieces', &
      'ERROR: TOMOPIECES - ', .false., 1, 1, 0, numOptArg, &
      numNonOptArg)
  call get_nxyz(.true., 'TomogramOrSizeXYZ', 'TOMOPIECES', 1, nxyz)

  ierr = PipGetFloat('MegaVoxels', megaVox)
  ierr = PipGetInteger('MinimumOverlap', minOverlap)
  ierr = PipGetInteger('XPadding', nxPad)
  ierr = PipGetInteger('YPadding', nyPad)
  ierr = PipGetInteger('ZPadding', nzPad)
  ierr = PipGetInteger('XMaximumPieces', maxPiecesX)
  ierr = PipGetInteger('YMaximumPieces', maxPiecesY)
  ierr = PipGetInteger('ZMaximumPieces', maxPiecesZ)
  ierr = PipGetLogical('NoFFTSizes', noFFT)
  call PipDone()
  !
  nx = nxyz(1)
  ny = nxyz(2)
  nz = nxyz(3)
  !
  ! Set defaults very big if no maxima entered - but if X and Y maxes
  ! are still 0, set up for 19 pieces on a layer
  ! Also, if one is zero and one is 1, set up for the 19-piece limit
  !
  if (maxPiecesX == 0 .and. (maxPiecesY == 0 .or. maxPiecesY == 1)) then
    maxPiecesX = maxLayerPieces
  else if (maxPiecesX == 1 .and. maxPiecesY == 0) then
    maxPiecesY = maxLayerPieces
  else
    if (maxPiecesX <= 0) maxPiecesX = nx / 2 - 1
    if (maxPiecesY <= 0) maxPiecesY = ny / 2 - 1
  endif
  if (maxPiecesZ < 0) maxPiecesZ = nz / 2 - 1
  !
  ! loop on the possible pieces in X; for each one, get the padded size
  !
  perimMin = (10. * nx) * ny * nz
  piecePerimMin = perimMin
  numXpieceMin = 0
  do numXpieces = 1, maxPiecesX
    nxExtra = (nx + (numXpieces - 1) * minOverlap + numXpieces - 1) / numXpieces
    nxOut = padNiceIfFFT(nxExtra, nxPad, noFFT)
    !
    ! loop on the possible pieces in Y
    !
    limitYpieces = maxPiecesY
    if (limitYpieces == 0) limitYpieces = maxLayerPieces / numXpieces

    do numYpieces = 1, limitYpieces
      nyExtra = (ny + (numYpieces - 1) * minOverlap + numYpieces - 1) / numYpieces
      nyOut = padNiceIfFFT(nyExtra, nyPad, noFFT)
      numZpieces = 1
      tooBig = .true.
      !
      ! then loop on pieces in Z, get padded size and compute padded
      ! volume size - until it is no longer too big for maxmem
      !
      do while(numZpieces <= maxPiecesZ .and. tooBig)
        nzExtra = (nz + (numZpieces - 1) * minOverlap + numZpieces - 1) / numZpieces
        nzOut = padNiceIfFFT(nzExtra, nzPad, noFFT)
        if (float(nxOut * nyOut) * nzOut > megaVox * 1.e6) then
          numZpieces = numZpieces + 1
        else
          tooBig = .false.
        endif
      enddo
      !
      ! compute perimeter of pieces and keep track of minimum
      ! If total perimeter is equal, favor one with minimum individual
      ! perimeter (at the expense of more pieces)
      ! Of the equivalent sets when nx = nz, this favors the ones with
      ! fewer X pieces
      !
      perim = float(nx * ny) * numZpieces + float(nx * nz) * numYpieces + float(ny * nz) &
          * numXpieces
      piecePerim = nxOut * nyOut + nxOut * nzOut + nyOut * nzOut
      ! print *,nxp, nyp, nzp, perim, piecePerim, nxout, nyout, nzout
      if (.not.tooBig .and. (perim < perimMin .or. &
          (perim == perimMin .and. piecePerim < piecePerimMin))) then
        perimMin = perim
        piecePerimMin = piecePerim
        numXpieceMin = numXpieces
        numYpieceMin = numYpieces
        numZpieceMin = numZpieces
      endif
    enddo
  enddo

  if (numXpieceMin == 0) call exitError( &
      'Pieces are all too large with given maximum numbers')
  numXpieces = numXpieceMin
  numZpieces = numZpieceMin
  numYpieces = numYpieceMin
  !
  ! get the starting and ending limits to extract and coordinates
  ! for getting back from the padded volume
  call getRanges(nx, numXpieces, minOverlap, nxPad, ixOutStart, ixOutEnd, &
      ixBackStart, ixBackEnd, noFFT)
  call getRanges(ny, numYpieces, minOverlap, nyPad, itOutStart, iyOutEnd, &
      iyBackStart, iyBackEnd, noFFT)
  call getRanges(nz, numZpieces, minOverlap, nzPad, izOutStart, izOutEnd, &
      izBackStart, izBackEnd, noFFT)
  !
  write(*,101) numXpieces, numYpieces, numZpieces
101 format(3i4)
  do iz = 1, numZpieces
    do iy = 1, numYpieces
      do ix = 1, numXpieces
        call rangeOut(ixOutStart(ix), ixOutEnd(ix), ',', outputFile, numOut)
        call rangeAdd(itOutStart(iy), outputFile, numOut)
        call rangeAdd(iyOutEnd(iy), outputFile, numOut)
        call rangeAdd(izOutStart(iz), outputFile, numOut)
        call rangeAdd(izOutEnd(iz), outputFile, numOut)
        write(*,102) outputFile(1:numOut)
102     format(a)
      enddo
    enddo
  enddo
  do i = 1, numXpieces
    call rangeOut(ixBackStart(i), ixBackEnd(i), ',', outputFile, numOut)
    write(*,102) outputFile(1:numOut)
  enddo
  do i = 1, numYpieces
    call rangeOut(iyBackStart(i), iyBackEnd(i), ',', outputFile, numOut)
    write(*,102) outputFile(1:numOut)
  enddo
  do i = 1, numZpieces
    call rangeOut(izBackStart(i), izBackEnd(i), ',', outputFile, numOut)
    write(*,102) outputFile(1:numOut)
  enddo
  !
  call exit(0)
end program tomopieces

! get size to extract, and size of padded output, and offset for
! getting back from the padded volume
!
subroutine getRanges(nx, numXpieces, minOverlap, nxPad, ixOutStart, ixOutEnd, &
    ixBackStart, ixBackEnd, noFFT)
  implicit none
  integer*4 nx, numXpieces, minOverlap, nxPad
  integer*4 ixOutStart(*), ixOutEnd(*), ixBackStart(*), ixBackEnd(*)
  integer*4 nxExtra, nxOut, nxBackOffset, lpaTotal, lapBase, lapExtra, ip, lap
  integer*4 lapBottom, lapTop, padNiceIfFFT
  logical*4 noFFT
  !
  nxExtra = (nx + (numXpieces - 1) * minOverlap + numXpieces - 1) / numXpieces
  nxOut = padNiceIfFFT(nxExtra, nxPad, noFFT)
  nxBackOffset = (nxOut - nxExtra) / 2
  !
  ! divide the total overlap into nearly equal parts
  !
  lpaTotal = nxExtra * numXpieces - nx
  lapBase = lpaTotal / max(1, numXpieces - 1)
  lapExtra = mod(lpaTotal, max(1, numXpieces - 1))
  ixOutStart(1) = 0
  ixBackStart(1) = nxBackOffset
  !
  ! get coordinates to extract, and coordinates for reassembly
  !
  do ip = 1, numXpieces
    ixOutEnd(ip) = ixOutStart(ip) + nxExtra - 1
    if (ip < numXpieces) then
      lap = lapBase
      if (ip <= lapExtra) lap = lap + 1
      lapTop = lap / 2
      lapBottom = lap - lapTop
      ixOutStart(ip + 1) = ixOutStart(ip) + nxExtra - lap
      ixBackEnd(ip) = nxBackOffset + nxExtra - 1 - lapTop
      ixBackStart(ip + 1) = nxBackOffset + lapBottom
    else
      ixBackEnd(ip) = nxBackOffset + nxExtra - 1
    endif
  enddo
  return
end subroutine getRanges

integer*4 function padNiceIfFFT(nxExtra, nxPad, noFFT)
  implicit none
  integer*4 nxExtra, nxPad
  logical*4 noFFT
  integer*4 niceFrame
  if (noFFT) then
    padNiceIfFFT = nxExtra + 2 * nxPad
  else
    padNiceIfFFT = niceFrame(2 * ((nxExtra + 1) / 2 + nxPad), 2, 19)
  endif
  return
end function padNiceIfFFT

subroutine rangeOut(izStart, izEnd, link, buf, numOut)
  character*(*) buf
  character*1 link
  call int_iwrite(buf, izStart, numOut)
  buf(numOut + 1:numOut + 1) = link
  call int_iwrite(buf(numOut + 2:), izEnd, numAdd)
  numOut = numOut + numAdd + 1
  return
end subroutine rangeOut

subroutine rangeAdd(izEnd, buf, numOut)
  character*(*) buf
  buf(numOut + 1:numOut + 1) = ','
  call int_iwrite(buf(numOut + 2:), izEnd, numAdd)
  numOut = numOut + numAdd + 1
  return
end subroutine rangeAdd

