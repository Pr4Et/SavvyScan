! BSUBS.FOR has subroutines needed by BLEND only:
!
! *READ_LIST
! *FUNCTION ONEINTRP
! *FASTINTERP
! JOINT_TO_ROTRANS
! EDGE_TO_ROTRANS
! SOLVE_ROTRANS
! FUNCTION STANDEV
! EDGESWAP
! function standev
! readEdgeFunc
! DOEDGE
! LINCOM_ROTRANS
! RECEN_ROTRANS
! COUNTEDGES
! mostInteriorFrames
! constrainFrames
! DXYDGRINTERP
! positionInPiece
! initNearList
! findNearestPiece
! DXYDGRINTERP
! CROSSVALUE
! *XCORREDGE
! FIND_BEST_SHIFTS
! includeEdge
! checkGroup
! sortVarsIntoGroups
! findBestGradient
! gradfunc
! findEdgeToUse
! getDataLimits
! IWRBINNED
! GETEXTRAINDENTS
! readExclusionModel
! DUMPEDGE
!
! $Id$

! Ones marked with * HAVE had crufty substitution.  Be sure to exclude bvvars next time

! READ_LIST get the name of a piece list file, reads the file
! and returns the # of pieces NPCLIST, the min and max z values MINZPC
! and MAXZPC, the piece coordinates I[XYZ]PCLIST, the optional negative
! number within each section, NEGLIST, and a logical variable MULTINEG
! that is true if the section
!
subroutine read_list(ixPcList, iyPcList, izPcList, negList, &
    multiNeg, numPieceList, minZpiece, maxZpiece, anyNeg, pipInput)
  !
  implicit none
  integer*4 ixPcList(*), iyPcList(*), izPcList(*), negList(*)
  logical multiNeg(*), anyNeg, pipInput
  integer*4 numPieceList, minZpiece, maxZpiece
  logical gotFirst, anyZero
  real*4 freeInput(10)
  character*320 fileName
  character*120 dummy
  character*32 errMess(4) / &
      'error opening file', &
      'error reading file, line', &
      'bad number of values on line', &
      'bad multinegative specification'/
  integer*4 ierr, lenAct, numInput, ineg, ipc, iz, listFirst
  integer*4 PipGetString
  !
  ! get file name and open file
  !
  if (pipInput) then
    if (PipGetString('PieceListInput', fileName) .ne. 0) call exitError &
        ('NO INPUT PIECE LIST FILE SPECIFIED')
  else
    write(*,'(1x,a,$)') 'name of input piece list file: '
    read(5, '(a)') fileName
  endif
  ierr = 1
  numPieceList = 0
  anyNeg = .false.
  ! 7/14/00 CER: remove carriagecontrol for LINUX
  open(3, file = fileName, form = 'formatted', status = 'old' &
      , err = 20)
  minZpiece = 100000
  maxZpiece = -100000
  !
  ! read each line in turn, get numbers with free format input
  !
12 ierr = 2
  read(3, '(a)', err = 20, end = 14) dummy
  lenAct = len(dummy)
  do while(dummy(lenAct:lenAct) == ' ' .or. dummy(lenAct:lenAct) == char(0))
    lenAct = lenAct - 1
    if (lenAct <= 0) go to 12
  enddo
  call frefor(dummy, freeInput, numInput)
  ierr = 3
  if (numInput < 3 .or. numInput > 4) go to 20        !error if fewer than 3 numbers
  numPieceList = numPieceList + 1
  ixPcList(numPieceList) = freeInput(1)
  iyPcList(numPieceList) = freeInput(2)
  izPcList(numPieceList) = freeInput(3)
  negList(numPieceList) = 0                        !if 4th number, its a neg #
  if (numInput == 4) negList(numPieceList) = freeInput(4)
  minZpiece = min(minZpiece, izPcList(numPieceList))
  maxZpiece = max(maxZpiece, izPcList(numPieceList))
  go to 12
  !
  ! now check for multinegative montaging: all pieces in a section must
  ! be labeled by negative if any are
  !
14 ierr = 4
  do iz = minZpiece, maxZpiece
    ineg = iz + 1 - minZpiece
    multiNeg(ineg) = .false.
    anyZero = .false.
    gotFirst = .false.
    do ipc = 1, numPieceList
      if (izPcList(ipc) == iz) then
        anyZero = anyZero .or. (negList(ipc) == 0)
        if (gotFirst) then
          multiNeg(ineg) = &
              multiNeg(ineg) .or. (negList(ipc) .ne. listFirst)
        else
          listFirst = negList(ipc)
          gotFirst = .true.
        endif
      endif
    enddo
    if (multiNeg(ineg) .and. anyZero) go to 20
    anyNeg = anyNeg .or. multiNeg(ineg)
  enddo
  close(3)
  return
20 write(*,'(1x,a,a,a,a,i6)') 'ERROR: BLENDMONT - ', &
      trim(fileName), ' - ', errMess(ierr), numPieceList + 1
  close(3)
  call exit(1)
end subroutine read_list



! function ONEINTRP used interpolation in a real array CRRAY
! (dimensions NX, NY) to find the pixel value of the point at
! coordinates (X, Y) .  Coordinates run from 0 to NX-1.  A real value
! is returned.  If the pixel is not within the array, DFILL is returned
! (from value in common) .  IZPC is the piece number (numbered from 1)
! The interpolation order is set from the value in common.
!
real*4 function oneintrp(arrIn, nx, ny, x, y, izPiece)
  !
  use blendvars
  implicit none
  integer*4 nx, ny, izPiece
  real*4 x, y
  real*4 arrIn(nx,ny)
  integer*4 ixp, iyp, ixpP1, ixpM1, iypP1, iypM1, ixpP2, iypP2
  real*4 dx, dy, xp, yp, v2, v4, v6, v8, v5, a, b, c, d, vmin, vmax
  real*4 dxM1, dxDxM1, fx1, fx2, fx3, fx4, dyM1, dyDyM1, v1, v3
  !
  oneintrp = dfill
  xp = x + 1.
  yp = y + 1.
  if (doFields) then
    call interpolateGrid(xp, yp, fieldDx(1, 1, memIndex(izPiece)), &
        fieldDy(1, 1, memIndex(izPiece)), lmField, nxField, nyField, &
        xFieldStrt, yFieldStrt, xFieldIntrv, yFieldIntrv, dx, dy)
    xp = xp + dx
    yp = yp + dy
  endif

  if (interpOrder <= 1) then
    !
    ! Linear interpolation
    !
    ixp = xp
    iyp = yp
    if (ixp >= 1 .and. ixp <= nx .and. iyp >= 1 .and. &
        iyp <= ny) then
      dx = xp - ixp
      dy = yp - iyp
      ixpP1 = min(nx, ixp + 1)
      iypP1 = min(ny, iyp + 1)
      oneintrp = (1. - dy) * ((1. - dx) * arrIn(ixp, iyp) + &
          dx * arrIn(ixpP1, iyp)) + &
          dy * ((1. - dx) * arrIn(ixp, iypP1) + &
          dx * arrIn(ixpP1, iypP1))
    endif
  elseif (interpOrder == 2) then
    !
    ! Old quadratic interpolation
    !
    ixp = nint(xp)
    iyp = nint(yp)
    if (ixp < 1 .or. ixp > nx .or. iyp < 1 .or. iyp > ny) return
    dx = xp - ixp
    dy = yp - iyp
    ! but if on an integer boundary already, done
    if (dx == 0 .and. dy == 0.) then
      oneintrp = arrIn(ixp, iyp)
      return
    endif
    !
    ixpP1 = min(nx, ixp + 1)
    ixpM1 = max(1, ixp - 1)
    iypP1 = min(ny, iyp + 1)
    iypM1 = max(1, iyp - 1)
    !
    ! Set up terms for quadratic interpolation
    !
    v2 = arrIn(ixp, iypM1)
    v4 = arrIn(ixpM1, iyp)
    v5 = arrIn(ixp, iyp)
    v6 = arrIn(ixpP1, iyp)
    v8 = arrIn(ixp, iypP1)
    ! find min and max of all 5 points
    vmax = max(v2, v4, v5, v6, v8)
    vmin = min(v2, v4, v5, v6, v8)
    !
    a = (v6 + v4) * .5 - v5
    b = (v8 + v2) * .5 - v5
    c = (v6 - v4) * .5
    d = (v8 - v2) * .5
    !
    ! limit the new density to between the min and max of original points
    oneintrp = max(vmin, min(vmax, a * dx * dx + b * dy * dy + c * dx + d * dy + v5))
  else
    !
    ! cubic interpolation
    !
    ixp = xp
    iyp = yp
    if (ixp >= 1 .and. ixp <= nx .and. iyp >= 1 .and. &
        iyp <= ny) then

      dx = xp - ixp
      dy = yp - iyp
      ixpP1 = min(nx, ixp + 1)
      ixpM1 = max(1, ixp - 1)
      iypP1 = min(ny, iyp + 1)
      iypM1 = max(1, iyp - 1)
      ixpP2 = min(nx, ixp + 2)
      iypP2 = min(ny, iyp + 2)

      dxM1 = dx - 1.
      dxDxM1 = dx * dxM1
      fx1 = -dxM1 * dxDxM1
      fx4 = dx * dxDxM1
      fx2 = 1 + dx**2 * (dx - 2.)
      fx3 = dx * (1. -dxDxM1)

      dyM1 = dy - 1.
      dyDyM1 = dy * dyM1

      v1 = fx1 * arrIn(ixpM1, iypM1) + fx2 * arrIn(ixp, iypM1) + &
          fx3 * arrIn(ixpP1, iypM1) + fx4 * arrIn(ixpP2, iypM1)
      v2 = fx1 * arrIn(ixpM1, iyp) + fx2 * arrIn(ixp, iyp) + &
          fx3 * arrIn(ixpP1, iyp) + fx4 * arrIn(ixpP2, iyp)
      v3 = fx1 * arrIn(ixpM1, iypP1) + fx2 * arrIn(ixp, iypP1) + &
          fx3 * arrIn(ixpP1, iypP1) + fx4 * arrIn(ixpP2, iypP1)
      v4 = fx1 * arrIn(ixpM1, iypP2) + fx2 * arrIn(ixp, iypP2) + &
          fx3 * arrIn(ixpP1, iypP2) + fx4 * arrIn(ixpP2, iypP2)
      oneintrp = -dyM1 * dyDyM1 * v1 + (1. +dy**2 * (dy - 2.)) * v2 + &
          dy * (1. -dyDyM1) * v3 + dy * dyDyM1 * v4

    endif
  endif
  return
end function oneintrp



! FASTINTERP uses quadratic interpolation to fill a portion of an
! output array DRRAY (dimensions NXD x NYD) from pixels in
! the input real array CRRAY (dimensions NXC x NYC) .  It will fill an
! area within the coordinates from INDXLO to INDXHI and from INDYLO to
! INDYHI, where the x coordinate of the lower left corner of the output
! array is NEWPCXLL.  The transformation is specified by the
! 2x2 matrix AMAT and FDX and FDY, where these must be prepared so that
! the coordinates in the input array (range 0 to NXC-1 etc) are
! obtained simply by X = A11 * INDX + A12 * INDY + DX et (i.e. no
! center offsets are needed) .  Pixels outside the input array are
! filled with DFILL from common.  IZPC is the piece number.
!
subroutine fastInterp(arrOut, nxOutput, nyOutput, arrIn, nxInput, nyInput, indXlow, &
    inXhigh, indYlow, indYhigh, newPcXlowLeft, amat, fdx, fdy, izPiece)
  !
  use blendvars
  implicit none
  integer*4 nxOutput, nyOutput, nxInput, nyInput, indXlow, inXhigh, indYlow, indYhigh
  integer*4 newPcXlowLeft, izPiece
  real*4 arrOut(nxOutput,nyOutput)
  real*4 arrIn(nxInput,nyInput)
  real*4 amat(2,2), dx, dy, fdx, fdy, xbase, ybase
  integer*4 indY, iyOut, indX, ixOut, ixp, iyp, ixpP1, ixpM1, iypP1, iypM1
  integer*4 ixpP2, iypP2
  real*4 pixVal, xp, yp, v2, v4, v6, v8, v5, a, b, c, d, vmin, vmax, rx, ry
  real*4 dxM1, dxDxM1, fx1, fx2, fx3, fx4, dyM1, dyDyM1, v1, v3
  !
  do indY = indYlow, indYhigh
    iyOut = indY + 1 - indYlow
    ry = indY
    xbase = amat(1, 2) * ry + fdx
    ybase = amat(2, 2) * ry + fdy
    if (interpOrder <= 1) then
      !
      ! Linear interpolation
      !
      do indX = indXlow, inXhigh
        ixOut = indX + 1 - newPcXlowLeft
        pixVal = dfill
        rx = indX
        if (secHasWarp) then
          call interpolateGrid(rx - 0.5, indY - 0.5, warpDx, warpDy, lmWarpX, nxWarp, &
              nyWarp, xWarpStrt, yWarpStrt, xWarpIntrv, yWarpIntrv, dx, dy)
          ! if (indy == 3845) print *,'fi', indx, rx, ry, dx, dy
          rx = rx + dx
          ry = indY + dy
          xbase = amat(1, 2) * ry + fdx
          ybase = amat(2, 2) * ry + fdy
        endif
        xp = amat(1, 1) * rx + xbase
        yp = amat(2, 1) * rx + ybase
        if (doFields) then
          call interpolateGrid(xp, yp, fieldDx(1, 1, memIndex(izPiece)), &
              fieldDy(1, 1, memIndex(izPiece)), lmField, nxField, nyField, &
              xFieldStrt, yFieldStrt, xFieldIntrv, yFieldIntrv, dx, dy)
          xp = xp + dx
          yp = yp + dy
        endif
        ixp = xp
        iyp = yp
        if (ixp >= 1 .and. ixp < nxInput .and. iyp >= 1 .and. &
            iyp < nyInput) then
          dx = xp - ixp
          dy = yp - iyp
          ixpP1 = ixp + 1
          iypP1 = iyp + 1
          pixVal = (1. - dy) * ((1. - dx) * arrIn(ixp, iyp) + &
              dx * arrIn(ixpP1, iyp)) + &
              dy * ((1. - dx) * arrIn(ixp, iypP1) + &
              dx * arrIn(ixpP1, iypP1))
        endif
        arrOut(ixOut, iyOut) = pixVal
      enddo
      !
    elseif (interpOrder == 2) then
      !
      ! Old quadratic interpolation
      !
      do indX = indXlow, inXhigh
        ixOut = indX + 1 - newPcXlowLeft
        pixVal = dfill
        rx = indX
        if (secHasWarp) then
          call interpolateGrid(rx - 0.5, indY - 0.5, warpDx, warpDy, lmWarpX, nxWarp, &
              nyWarp, xWarpStrt, yWarpStrt, xWarpIntrv, yWarpIntrv, dx, dy)
          rx = rx + dx
          ry = indY + dy
          xbase = amat(1, 2) * ry + fdx
          ybase = amat(2, 2) * ry + fdy
        endif
        xp = amat(1, 1) * rx + xbase
        yp = amat(2, 1) * rx + ybase
        if (doFields) then
          call interpolateGrid(xp, yp, fieldDx(1, 1, memIndex(izPiece)), &
              fieldDy(1, 1, memIndex(izPiece)), lmField, nxField, nyField, &
              xFieldStrt, yFieldStrt, xFieldIntrv, yFieldIntrv, dx, dy)
          xp = xp + dx
          yp = yp + dy
        endif
        ixp = nint(xp)
        iyp = nint(yp)
        if (ixp < 1 .or. ixp > nxInput .or. iyp < 1 .or. iyp > nyInput) &
            go to 80
        !
        ! Do quadratic interpolation
        !
        dx = xp - ixp
        dy = yp - iyp
        v5 = arrIn(ixp, iyp)
        !
        ! but if on an integer boundary already, done
        !
        if (dx == 0 .and. dy == 0.) then
          pixVal = v5
          go to 80
        endif
        !
        ixpP1 = min(nxInput, ixp + 1)
        ixpM1 = max(1, ixp - 1)
        iypP1 = min(nyInput, iyp + 1)
        iypM1 = max(1, iyp - 1)
        !
        ! Set up terms for quadratic interpolation
        !
        v2 = arrIn(ixp, iypM1)
        v4 = arrIn(ixpM1, iyp)
        v6 = arrIn(ixpP1, iyp)
        v8 = arrIn(ixp, iypP1)
        !
        ! find min and max of all 5 points
        !
        vmax = max(v2, v4, v5, v6, v8)
        vmin = min(v2, v4, v5, v6, v8)
        !
        a = (v6 + v4) * .5 - v5
        b = (v8 + v2) * .5 - v5
        c = (v6 - v4) * .5
        d = (v8 - v2) * .5
        !
        ! limit the new density to between min and max of original points
        !
        pixVal = max(vmin, min(vmax, &
            a * dx * dx + b * dy * dy + c * dx + d * dy + v5))
80      arrOut(ixOut, iyOut) = pixVal
      enddo
    else
      !
      ! cubic interpolation
      !
      do indX = indXlow, inXhigh
        ixOut = indX + 1 - newPcXlowLeft
        pixVal = dfill
        rx = indX
        if (secHasWarp) then
          call interpolateGrid(rx - 0.5, indY - 0.5, warpDx, warpDy, lmWarpX, nxWarp, &
              nyWarp, xWarpStrt, yWarpStrt, xWarpIntrv, yWarpIntrv, dx, dy)
          rx = rx + dx
          ry = indY + dy
          xbase = amat(1, 2) * ry + fdx
          ybase = amat(2, 2) * ry + fdy
        endif
        xp = amat(1, 1) * rx + xbase
        yp = amat(2, 1) * rx + ybase
        if (doFields) then
          call interpolateGrid(xp, yp, fieldDx(1, 1, memIndex(izPiece)), &
              fieldDy(1, 1, memIndex(izPiece)), lmField, nxField, nyField, &
              xFieldStrt, yFieldStrt, xFieldIntrv, yFieldIntrv, dx, dy)
          xp = xp + dx
          yp = yp + dy
        endif
        ixp = xp
        iyp = yp
        if (ixp >= 2 .and. ixp < nxInput - 1 .and. iyp >= 2 .and. &
            iyp < nyInput - 1) then

          dx = xp - ixp
          dy = yp - iyp
          ixpP1 = ixp + 1
          ixpM1 = ixp - 1
          iypP1 = iyp + 1
          iypM1 = iyp - 1
          ixpP2 = ixp + 2
          iypP2 = iyp + 2

          dxM1 = dx - 1.
          dxDxM1 = dx * dxM1
          fx1 = -dxM1 * dxDxM1
          fx4 = dx * dxDxM1
          fx2 = 1 + dx**2 * (dx - 2.)
          fx3 = dx * (1. -dxDxM1)

          dyM1 = dy - 1.
          dyDyM1 = dy * dyM1

          v1 = fx1 * arrIn(ixpM1, iypM1) + fx2 * arrIn(ixp, iypM1) + &
              fx3 * arrIn(ixpP1, iypM1) + fx4 * arrIn(ixpP2, iypM1)
          v2 = fx1 * arrIn(ixpM1, iyp) + fx2 * arrIn(ixp, iyp) + &
              fx3 * arrIn(ixpP1, iyp) + fx4 * arrIn(ixpP2, iyp)
          v3 = fx1 * arrIn(ixpM1, iypP1) + fx2 * arrIn(ixp, iypP1) + &
              fx3 * arrIn(ixpP1, iypP1) + fx4 * arrIn(ixpP2, iypP1)
          v4 = fx1 * arrIn(ixpM1, iypP2) + fx2 * arrIn(ixp, iypP2) + &
              fx3 * arrIn(ixpP1, iypP2) + fx4 * arrIn(ixpP2, iypP2)
          pixVal = -dyM1 * dyDyM1 * v1 + (1. +dy**2 * (dy - 2.)) * v2 + &
              dy * (1. -dyDyM1) * v3 + dy * dyDyM1 * v4
          !
        endif
        arrOut(ixOut, iyOut) = pixVal
      enddo
    endif
  enddo
  return
end subroutine fastInterp





! JOINT_TO_ROTRANS uses the values in the list of NEDGE edge functions
! in DXGRID and DYGRID (dimensions IXDIM x IYDIM, # of values in X and
! Y directions NXGRID x NYGRID, pixel interval between values INTXGRID
! and INTYGRID in X and Y directions, x and y coordinates of the lower-
! left corner of the edge in the lower and upper piece in IXPCLO,
! IYPCLO and IXPCHI, IYPCHI) to derive a rotation/translation
! around the center of the joint; it returns the angle in THETAMIN and
! the translation in DXMIN and DYMIN.
!
subroutine joint_to_rotrans(dxgrid, dygrid, ixdim, iydim, nxgrid, &
    nygrid, intxgrid, intygrid, ixpclo, iypclo, ixpchi, iypchi, nedge &
    , r)
  !
  real*4 dxgrid(ixdim,iydim,*), dygrid(ixdim,iydim,*)
  integer*4 nxgrid(*), nygrid(*), ixpclo(*), ixpchi(*), iypclo(*) &
      , iypchi(*)
  !
  ! structure /rotrans/
  ! real*4 theta, dx, dy, xcen, ycen
  ! end structure
  real*4 r(6)
  !
  parameter (npnts = 5000)
  real*4 x1(npnts), y1(npnts), x2(npnts), y2(npnts)
  !
  ! make set of coordinate pairs in each grid, add up sums
  x1sum = 0.
  x2sum = 0.
  y1sum = 0.
  y2sum = 0.
  nn = 0
  do ied = 1, nedge
    do ix = 1, nxgrid(ied)
      do iy = 1, nygrid(ied)
        nn = nn + 1
        x1(nn) = (ix - 1) * intxgrid + ixpclo(ied)
        y1(nn) = (iy - 1) * intygrid + iypclo(ied)
        x2(nn) = (ix - 1) * intxgrid + ixpchi(ied) + dxgrid(ix, iy, ied)
        y2(nn) = (iy - 1) * intygrid + iypchi(ied) + dygrid(ix, iy, ied)
        x1sum = x1sum + x1(nn)
        y1sum = y1sum + y1(nn)
      enddo
    enddo
  enddo
  !
  ! get center in lower image and shift points to center around there
  !
  r(4) = x1sum / nn
  r(5) = y1sum / nn
  do i = 1, nn
    x1(i) = x1(i) - r(4)
    x2(i) = x2(i) - r(4)
    y1(i) = y1(i) - r(5)
    y2(i) = y2(i) - r(5)
  enddo
  call solve_rotrans(x1, y1, x2, y2, nn, r(1), r(2), r(3))
  return
end subroutine joint_to_rotrans



! EDGE_TO_ROTRANS uses the values in an edge function specified by
! DXGRID and DYGRID (dimensions IXDIM x IYDIM, # of values in X and Y
! directions NXGRID x NYGRID, pixel interval between values INTXGRID
! and INTYGRID in X and Y directions) to derive a rotation/translation
! around the center of the edge; it returns the angle in THETAMIN and
! the translation in DXMIN and DYMIN.
!
subroutine edge_to_rotrans(dxgrid, dygrid, ixdim, iydim, nxgrid, &
    nygrid, intxgrid, intygrid, thetamin, dxmin, dymin)
  !
  real*4 dxgrid(ixdim,iydim), dygrid(ixdim,iydim)
  parameter (npnts = 1000)
  real*4 x1(npnts), y1(npnts), x2(npnts), y2(npnts)
  !
  ! make a set of coordinate pairs centered around center of grid
  !
  xcen = (nxgrid-1) / 2. + 1.
  ycen = (nygrid-1) / 2. + 1.
  nn = 0
  do ix = 1, nxgrid
    do iy = 1, nygrid
      nn = nn + 1
      x1(nn) = (ix - xcen) * intxgrid
      y1(nn) = (iy - ycen) * intygrid
      x2(nn) = x1(nn) + dxgrid(ix, iy)
      y2(nn) = y1(nn) + dygrid(ix, iy)
    enddo
  enddo
  call solve_rotrans(x1, y1, x2, y2, nn, thetamin, dxmin, dymin)
  return
end subroutine edge_to_rotrans



! SOLVE_ROTRANS takes a list of NN coordinate pairs in X1, Y1 and X2, Y2
! and finds the rotation THETAMIN about the origin and translation
! DXMIN and DYMIN that best superimposes the points
!
subroutine solve_rotrans(x1, y1, x2, y2, nn, thetamin, dxmin, dymin)
  real*4 x1(*), y1(*), x2(*), y2(*)
  !
  ! the real way to do it
  !
  x1s = 0.
  x2s = 0.
  y1s = 0.
  y2s = 0.
  do i = 1, nn
    x1s = x1s + x1(i)
    x2s = x2s + x2(i)
    y1s = y1s + y1(i)
    y2s = y2s + y2(i)
  enddo
  x1s = x1s / nn
  x2s = x2s / nn
  y1s = y1s / nn
  y2s = y2s / nn
  ssd23 = 0.
  ssd14 = 0.
  ssd13 = 0.
  ssd24 = 0.
  do i = 1, nn
    ssd23 = ssd23 + (y1(i) - y1s) * (x2(i) - x2s)
    ssd14 = ssd14 + (x1(i) - x1s) * (y2(i) - y2s)
    ssd13 = ssd13 + (x1(i) - x1s) * (x2(i) - x2s)
    ssd24 = ssd24 + (y1(i) - y1s) * (y2(i) - y2s)
  enddo
  thetamin = 0.
  if (abs(ssd13 + ssd24) > 1.e-20) &
      thetamin = atand(-(ssd23 - ssd14) / (ssd13 + ssd24))
  costh = cosd(thetamin)
  sinth = sind(thetamin)
  dxmin = x2s - x1s * costh + y1s * sinth
  dymin = y2s - x1s * sinth - y1s * costh

!!$    c
!!$    c         search for angle that minimizes the sd of differences between pairs
!!$    c         set up for big range of angles
!!$    c
!!$    thetalo = -10.
!!$    thetahi = 10.
!!$    dtheta = 1.
!!$    cuttheta = 10.
!!$    niter = 4
!!$    sdmin = 1.e10
!!$    do iter = 1, niter
!!$    c
!!$    c           get number of steps between lo and hi
!!$    c
!!$    nsteps = (thetahi - thetalo) / dtheta + 1
!!$    dtheta = (thetahi - thetalo) / (nsteps + 1)
!!$    do istep = 1, nsteps
!!$    theta = thetalo + dtheta * (istep - 1)
!!$    costh = cosd(theta)
!!$    sinth = sind(theta)
!!$    dxsum = 0.
!!$    dxsumsq = 0.
!!$    dysum = 0.
!!$    dysumsq = 0.
!!$    c
!!$    c             for each theta value, back - rotate the points in X2, y2 and compute
!!$    c             sd of their difference from points in x1, y1
!!$    c
!!$    do i = 1, nn
!!$    dx = x2(i) * costh + y2(i) * sinth - x1(i)
!!$    dy = -x2(i) * sinth + y2(i) * costh - y1(i)
!!$    dxsum = dxsum + dx
!!$    dxsumsq = dxsumsq + dx**2
!!$    dysum = dysum + dy
!!$    dysumsq = dysumsq + dy**2
!!$    enddo
!!$    xsd = standev(dxsum, dxsumsq, nn)
!!$    ysd = standev(dysum, dysumsq, nn)
!!$    c
!!$    c             just minimize the sum of the two sd's
!!$    c
!!$    if (xsd + ysd<sdmin) then
!!$    sdmin = xsd + ysd
!!$    dxmin = dxsum / nn
!!$    dymin = dysum / nn
!!$    thetamin = theta
!!$    endif
!!$    enddo
!!$    c
!!$    c           for next round, go between two points around minimum, cut step size
!!$    c
!!$    thetalo = thetamin - dtheta
!!$    thetahi = thetamin + dtheta
!!$    dtheta = dtheta / cuttheta
!!$    enddo
  return
end subroutine solve_rotrans
!
function standev(sum, sumsq, nsum)
  standev = 0
  if (nsum <= 1) return
  diff = sumsq - sum**2 / nsum
  if (diff > 0.) standev = sqrt(diff / (nsum - 1))
  return
end function standev



! EDGESWAP manages the edge function buffers.  Given the edge # and
! type with  IEDGE and IXY, it looks to see if that function is already
! in the buffers, reads it in if not (replacing the edge function that
! has the longest time since last use), and returns the index of the
! edge in the buffers in INDBUF.
!
subroutine edgeswap(iedge, ixy, indbuf)
  !
  use blendvars
  implicit none
  integer*4 iedge, ixy, indbuf
  integer*4 minused, ioldest, i
  !
  indbuf = ibufedge(iedge, ixy)
  if (indbuf == 0) then
    !
    ! find oldest
    !
    minused = jusedgct + 1
    do i = 1, limedgbf
      if (minused > lasedguse(i)) then
        minused = lasedguse(i)
        ioldest = i
      endif
    enddo
    !
    ! mark oldest as no longer present, mark this as present
    !
    if (iedgbflist(ioldest) > 0) &
        ibufedge(iedgbflist(ioldest), ixybflist(ioldest)) = 0
    iedgbflist(ioldest) = iedge
    ixybflist(ioldest) = ixy
    indbuf = ioldest
    ibufedge(iedge, ixy) = indbuf
    !
    ! read edge into  buffer
    !
    call readEdgeFunc(iedge, ixy, ioldest)
  endif
  !
  ! mark this one as used most recently
  !
  jusedgct = jusedgct + 1
  lasedguse(indbuf) = jusedgct
  return
end subroutine edgeswap


! READEDGEFUNC actually reads in the edge function IEDGE, direction
! IXY, to the given buffer INDBUF, and swaps bytes as necessary
!
subroutine readEdgeFunc(iedge, ixy, indbuf)
  use blendvars
  implicit none
  integer*4 iedge, ixy, indbuf
  integer*4 ix, iy, nxgr, nygr, idum1, idum2
  !
  if (needbyteswap == 0) then
    read(iunedge(ixy), rec = 1 + iedge) &
        nxgr, nygr, ixgrdstbf(indbuf), ixofsbf(indbuf) &
        , iygrdstbf(indbuf), iyofsbf(indbuf) &
        , ((dxgrbf(ix, iy, indbuf), dygrbf(ix, iy, indbuf) &
        , ddengrbf(ix, iy, indbuf), ix = 1, nxgr), iy = 1, nygr)
  else
    read(iunedge(ixy), rec = 1 + iedge) nxgr, nygr
    call convert_longs(nxgr, 1)
    call convert_longs(nygr, 1)
    read(iunedge(ixy), rec = 1 + iedge) &
        idum1, idum2, ixgrdstbf(indbuf), ixofsbf(indbuf) &
        , iygrdstbf(indbuf), iyofsbf(indbuf) &
        , ((dxgrbf(ix, iy, indbuf), dygrbf(ix, iy, indbuf) &
        , ddengrbf(ix, iy, indbuf), ix = 1, nxgr), iy = 1, nygr)
    call convert_longs(ixgrdstbf(indbuf), 1)
    call convert_longs(ixofsbf(indbuf), 1)
    call convert_longs(iygrdstbf(indbuf), 1)
    call convert_longs(iyofsbf(indbuf), 1)
    do iy = 1, nygr
      call convert_floats(dxgrbf(1, iy, indbuf), nxgr)
      call convert_floats(dygrbf(1, iy, indbuf), nxgr)
      call convert_floats(ddengrbf(1, iy, indbuf), nxgr)
    enddo
  endif
  if (ifskipEdge(iedge, ixy) > 0) then
    dxgrid(1:nxgr, 1:nygr) = 0.
    dygrid(1:nxgr, 1:nygr) = 0.
    ddengrid(1:nxgr, 1:nygr) = 0.
  endif
  ! print *,'read all of edge', ixy, iedge, ixpclist(ipiecelower(iedge, ixy)), &
  ! iypclist(ipiecelower(iedge, ixy)), nxgr, nygr
  nxgrbf(indbuf) = nxgr
  nygrbf(indbuf) = nygr
  intxgrbf(indbuf) = intgrcopy(ixy)
  intygrbf(indbuf) = intgrcopy(3 - ixy)
  return
end subroutine readEdgeFunc


! DOEDGE "does" edge # IEDGE, whose direction is given by IXY.  It
! looks to see if this edge is on a joint between negatives; if so it
! composes the whole list of edges along that joint and arranges to
! find the edge functions for all of those edges, from the center
! outward.  Edge functions are found, then smoothed using the other
! arguments as parameters, then written to the appropriate edge file
!
subroutine doedge(iedge, ixy, edgedone, sdcrit, devcrit, nfit, &
    norder, nskip, docross, xcreadin, xclegacy, edgedispx, edgedispy, idimedge)
  !
  use blendvars
  implicit none
  ! real*4 array(*)
  integer*4 nfit(2), nskip(2)
  logical docross, xcreadin, xclegacy
  integer*4 iedge, ixy, norder, idimedge
  real*4 sdcrit, devcrit
  !
  logical edgedone(idimedge, 2)
  real*4 edgedispx(idimedge,2), edgedispy(idimedge,2)
  !
  integer limpneg
  parameter (limpneg = 20)
  integer*4 multcoord(limpneg), multedge(limpneg), multmp(limpneg) &
      , mcotmp(limpneg), igrstr(2), igrofs(2)
  !
  integer*4 intxgrid, intygrid, nmult, intscan, ipclo, ipcup, ipc, mltco
  integer*4 i, j, itmp, midcoord, mindiff, imult, imid, middone, indlow
  integer*4 indup, ixdisp, iydisp, ixdispmid, iydispmid, lastedge
  integer*4 lastxdisp, lastydisp, idiff, jedge, nxgr, nygr, ix, iy, indentXcorr
  real*4 xdisp, ydisp, theta, edgedx, edgedy, dxmid, dymid, xdispl, ydispl
  real*4 costh, sinth, xrel, yrel, thetamid, delIndent(2)
  integer*4 indentUse(2), limitLo, limitHi, limitLo2, limitHi2
  real*4 cosd, sind
  real*8 wallstart, walltime
  !
  ! make list of edges to be done
  !
  intxgrid = intgrid(ixy)
  intygrid = intgrid(3 - ixy)
  nmult = 1
  multedge(1) = iedge
  intscan = 6
  ipclo = ipiecelower(iedge, ixy)
  ipcup = ipieceupper(iedge, ixy)
  ipcBelowEdge = ipclo
  if (neglist(ipclo) .ne. neglist(ipcup)) then
    !
    ! if edge is across a negative boundary, need to look for all such
    ! edges and add them to list
    !
    nmult = 0
    intscan = 9
    do i = 1, nedge(ixy)
      ipc = ipiecelower(i, ixy)
      if (izpclist(ipclo) == izpclist(ipc) .and. &
          neglist(ipclo) == neglist(ipc) .and. &
          neglist(ipcup) == neglist(ipieceupper(i, ixy))) then
        nmult = nmult + 1
        multedge(nmult) = i
        ! get coordinate of edge in ortho direction
        if (ixy == 1) then
          mltco = iypclist(ipc)
        else
          mltco = ixpclist(ipc)
        endif
        multcoord(nmult) = mltco
      endif
    enddo
    !
    ! order list to go out from center of edge.  GROSS.. first order it
    !
    do i = 1, nmult - 1
      do j = i, nmult
        if (multcoord(i) > multcoord(j)) then
          itmp = multcoord(i)
          multcoord(i) = multcoord(j)
          multcoord(j) = itmp
          itmp = multedge(i)
          multedge(i) = multedge(j)
          multedge(j) = itmp
        endif
      enddo
    enddo
    midcoord = (multcoord(nmult) + multcoord(1)) / 2
    !
    ! find element closest to center
    !
    mindiff = 100000
    do i = 1, nmult
      idiff = abs(multcoord(i) - midcoord)
      if (idiff < mindiff) then
        mindiff = idiff
        imid = i
      endif
    enddo
    !
    ! set up order from there to top then back from middle to bottom
    !
    imult = 0
    do i = imid, nmult
      imult = imult + 1
      multmp(imult) = multedge(i)
      mcotmp(imult) = multcoord(i)
    enddo
    middone = imult
    do i = imid-1, 1, -1
      imult = imult + 1
      multmp(imult) = multedge(i)
      mcotmp(imult) = multcoord(i)
    enddo
    do i = 1, nmult
      multedge(i) = multmp(i)
      multcoord(i) = mcotmp(i)
    enddo
  endif
  !
  ! finally ready to set up to get edge
  !
  do imult = 1, nmult
    jedge = multedge(imult)
    !
    call shuffler(ipiecelower(jedge, ixy), indlow)
    call shuffler(ipieceupper(jedge, ixy), indup)
    !
    if (imult == 1) then
      !
      ! for first time, set these parameters to
      ! their default values with 0 offsets
      !
      ixdisp = 0
      iydisp = 0
      if (docross) then
        if (xcreadin) then
          xdisp = edgedispx(jedge, ixy)
          ydisp = edgedispy(jedge, ixy)
        else
          call getExtraIndents(ipiecelower(jedge, ixy), &
              ipieceupper(jedge, ixy), ixy, delIndent)
          indentXcorr = 2
          if (delIndent(ixy) > 0. .and. ifillTreatment == 1) &
              indentXcorr = int(delIndent(ixy)) + 2
          if (ifskipEdge(jedge, ixy) > 0) then
            xdisp = 0.
            ydisp = 0.
          else
            call xcorrEdge(array(indlow), array(indup), &
                ixy, xdisp, ydisp, xclegacy, indentXcorr)
          endif
          edgedispx(jedge, ixy) = xdisp
          edgedispy(jedge, ixy) = ydisp
        endif
        ixdisp = nint(xdisp)
        iydisp = nint(ydisp)
        !write(*,'(1x,a,2i4,a,2i5)') char(ixy + ichar('W')) //' edge, pieces' &
        !    , ipiecelower(jedge, ixy), ipieceupper(jedge, ixy), &
        !    '  ixydisp:', ixdisp, iydisp
      endif
      ixdispmid = ixdisp
      iydispmid = iydisp
      !
    else
      if (imult == middone + 1) then
        !
        theta = thetamid                      !at midway point, restore
        edgedx = dxmid                        !values from first (middle)
        edgedy = dymid                        !edge
        lastedge = multedge(1)
        lastxdisp = ixdispmid
        lastydisp = iydispmid
      else
        call edge_to_rotrans(dxgrid, dygrid, ixgdim, iygdim, nxgr, &
            nygr, intxgrid, intygrid, theta, edgedx, edgedy)
        if (imult == 2) then
          thetamid = theta                    !if that was first edge, save
          dxmid = edgedx                      !the value
          dymid = edgedy
        endif
        lastedge = multedge(imult - 1)
      endif
      !
      ! find displacement of center of next edge relative to center of
      ! last edge.  First get x/y displacements between edges
      !
      xdispl = ixpclist(ipieceupper(jedge, ixy)) - &
          ixpclist(ipieceupper(lastedge, ixy))
      ydispl = iypclist(ipieceupper(jedge, ixy)) - &
          iypclist(ipieceupper(lastedge, ixy))
      costh = cosd(theta)
      sinth = sind(theta)
      ! rotate vector by theta and displace by dx, dy; the movement in
      ! the tip of the displacement vector is the expected relative
      ! displacement between this frame and the last
      xrel = xdispl * costh - ydispl * sinth + edgedx - xdispl
      yrel = xdispl * sinth + ydispl * costh + edgedy - ydispl
      ! add pixel displacement of last frame to get total expected pixel
      ! displacement of this frame
      ixdisp = nint(xrel) + lastxdisp
      iydisp = nint(yrel) + lastydisp
    endif
    !
    ! Determine extra indentation if distortion corrections
    !
    call getExtraIndents(ipiecelower(jedge, ixy), ipieceupper(jedge, ixy), &
        ixy, delIndent)
    ! if (doFields) write(*,'(1x,a,2i4,a,2f5.1)') &
    ! char(ixy+ichar('W')) //' edge, pieces' &
    ! , ipiecelower(jedge, ixy), ipieceupper(jedge, ixy), &
    ! '  extra indents:', delIndent(1), delIndent(2)
    !
    indentUse(1) = indent(1) + nint(delIndent(1))
    indentUse(2) = indent(2) + nint(delIndent(2))
    !
    ! Determine data limits for edge in long dimension if flag set
    limitLo = 0
    limitHi = 0
    if (limitData) then
      call getDataLimits(ipiecelower(jedge, ixy), 3 - ixy, 2, limitLo, &
          limitHi)
      call getDataLimits(ipieceupper(jedge, ixy), 3 - ixy, 1, limitLo2, &
          limitHi2)
      limitLo = max(limitLo, limitLo2)
      limitHi = min(limitHi, limitHi2)
    endif
    !
    call setgridchars(nxyzin, noverlap, iboxsiz, indentUse, intgrid, &
        ixy, ixdisp, iydisp, limitLo, limitHi, nxgr, nygr, igrstr, igrofs)
    if (nxgr > nxgrid(ixy) .or. nygr > nygrid(ixy)) call exitError( &
        'ONE GRID HAS MORE POINTS THAN ORIGINALLY EXPECTED')
    lastxdisp = ixdisp
    lastydisp = iydisp
    !
    ! write(*,'(1x,a,2i4,a,2i4,a,2i5,a,2i4)') &
    ! char(ixy+ichar('W')) //' edge, pieces' &
    ! , ipiecelower(jedge, ixy), ipieceupper(jedge, ixy), &
    ! '  ngrid:', nxgr, nygr, '  start lower:', igrstr, &
    ! '  upper:', igrofs
    if (ifskipEdge(jedge, ixy) == 0) then
      wallstart = walltime()
      call findedgefunc(array(indlow), array(indup), nxin, nyin, &
          igrstr(1), igrstr(2), igrofs(1), igrofs(2), nxgr, nygr, &
          intxgrid, intygrid, iboxsiz(ixy), iboxsiz(3 - ixy), intscan, &
          dxgrid, dygrid, sdgrid, ddengrid, ixgdim, iygdim)
      !
      if (izUnsmoothedPatch >= 0) then
        do iy = 1, nygr
          do ix = 1, nxgr
            write(10, '(3i6,3f9.2,f12.4)') igrstr(1) + (ix - 1) * intxgrid, &
                igrstr(2) + (iy - 1) * intygrid, izUnsmoothedPatch, dxgrid(ix, iy), &
                dygrid(ix, iy), 0., sdgrid(ix, iy)
          enddo
        enddo
        izUnsmoothedPatch = izUnsmoothedPatch + 1
      endif
      call smoothgrid(dxgrid, dygrid, sdgrid, ddengrid, ixgdim, &
          iygdim, nxgr, nygr, sdcrit, devcrit, nfit(ixy), nfit(3 - ixy), &
          norder, nskip(ixy), nskip(3 - ixy))
      ! write(*,'(a,f10.6)') 'Edge function time', walltime() -wallstart
      if (izSmoothedPatch >= 0) then
        do iy = 1, nygr
          do ix = 1, nxgr
            write(11, '(3i6,3f9.2,f12.4)') igrstr(1) + (ix - 1) * intxgrid, &
                igrstr(2) + (iy - 1) * intygrid, izSmoothedPatch, dxgrid(ix, iy), &
                dygrid(ix, iy), 0., sdgrid(ix, iy)
          enddo
        enddo
        izSmoothedPatch = izSmoothedPatch + 1
      endif
    else
      dxgrid(1:nxgr, 1:nygr) = 0.
      dygrid(1:nxgr, 1:nygr) = 0.
      ddengrid(1:nxgr, 1:nygr) = 0.
    endif
    !
!!$        xrel = 0.
!!$        yrel = 0.
!!$        xdispl = 0.
!!$        do ix = 1, nxgr
!!$          do iy = 1, nygr
!!$            costh = sqrt(dxgrid(ix, iy)**2 + dygrid(ix, iy)**2)
!!$            xrel = xrel + costh
!!$            yrel = max(yrel, costh)
!!$            xdispl = xdispl + ddengrid(ix, iy)
!!$          enddo
!!$        enddo
!!$        xdispl = xdispl / (nxgr * nygr)
!!$        write(*,'(1x,a,2i4,a,2f6.2,a,f8.1)') &
!!$            char(ixy + ichar('W')) //' edge, pieces' &
!!$            , ipiecelower(jedge, ixy), ipieceupper(jedge, ixy), &
!!$            '  mean, max vector:', xrel / (nxgr * nygr), yrel, ' den diff:', xdispl

    if (needByteSwap == 0) then
      write(iunedge(ixy), rec = jedge + 1) nxgr, nygr, (igrstr(i), igrofs(i) &
          , i = 1, 2) , ((dxgrid(ix, iy), dygrid(ix, iy), &
          ddengrid(ix, iy), ix = 1, nxgr), iy = 1, nygr)
    else
      !
      ! In case there were incomplete edges, be able to write swapped
      ixdisp = nxgr
      iydisp = nygr
      call convert_longs(ixdisp, 1)
      call convert_longs(iydisp, 1)
      call convert_longs(igrstr, 2)
      call convert_longs(igrofs, 2)
      do iy = 1, nygr
        call convert_floats(dxgrid(1, iy), nxgr)
        call convert_floats(dygrid(1, iy), nxgr)
        call convert_floats(ddengrid(1, iy), nxgr)
      enddo
      write(iunedge(ixy), rec = jedge + 1) ixdisp, iydisp, (igrstr(i), igrofs(i) &
          , i = 1, 2) , ((dxgrid(ix, iy), dygrid(ix, iy), &
          ddengrid(ix, iy), ix = 1, nxgr), iy = 1, nygr)
    endif
    !
    edgedone(jedge, ixy) = .true.
    !
    ! Write records for any edges not done yet
    ixdisp = -1
    if (needByteSwap .ne. 0) call convert_longs(ixdisp, 1)
    do ix = lastWritten(ixy) + 1, jedge - 1
      if (.not.edgedone(ix, ixy)) &
          write(iunedge(ixy), rec = ix + 1) ixdisp, ixdisp
    enddo
    lastWritten(ixy) = jedge
  enddo
  return
end subroutine doedge



! LINCOM_ROTRANS forms a linear combination of two rotation/
! translations R1 and R2, with weights W1 and W2, shifts the result to
! the center of R2, and puts the result in S
!
subroutine lincom_rotrans(r1, w1, r2, w2, s)
  !
  ! structure /rotrans/
  ! real*4 theta, dx, dy, xcen, ycen
  ! end structure
  real*4 r1(6), r2(6), rt1(6), rt2(6), s(6)
  !
  call recen_rotrans(r1, r2(4), r2(5), rt1)
  call xfcopy(r2, rt2)
  rt2(1) = w1 * rt1(1) + w2 * r2(1)
  rt2(2) = w1 * rt1(2) + w2 * r2(2)
  rt2(3) = w1 * rt1(3) + w2 * r2(3)
  call xfcopy(rt2, s)
  return
end subroutine lincom_rotrans



! RECEN_ROTRANS shift the rotation/translation R to the new center
! XCENEW, YCENEW and returns the result in S
!
subroutine recen_rotrans(r, xcenew, ycenew, s)
  !
  ! structure /rotrans/
  ! real*4 theta, dx, dy, xcen, ycen
  ! end structure
  real*4 r(6), s(6)
  !
  sinth = sind(r(1))
  cosm1 = cosd(r(1)) - 1.
  s(1) = r(1)
  s(2) = r(2) + cosm1 * (xcenew - r(4)) - sinth * (ycenew - r(5))
  s(3) = r(3) + cosm1 * (ycenew - r(5)) + sinth * (xcenew - r(4))
  s(4) = xcenew
  s(5) = ycenew
  return
end subroutine recen_rotrans

!!$    subroutine recen_xform(f, dxcen, dycen, g)
!!$    structure / xform/
!!$    real*4 a(2,2), dx, dy
!!$    end structure
!!$    record / xform / f, g
!!$    g = f
!!$    g.dx = f.dx + (f.a(1, 1) - 1.) * dxcen + f.a(1, 2) * dycen
!!$    g.dy = f.dy + (f.a(2, 2) - 1.) * dycen + f.a(2, 1) * dxcen
!!$    return
!!$    end



! COUNTEDGES takes coordinate INDX, INDY in the output image, converts
! it to XG, YG with the inverse of optional g transform, and analyses
! the edges and pieces that the point is in, or at least near
!
subroutine countedges(indx, indy, xg, yg, useEdges)
  !
  use blendvars
  implicit none
  integer*4 indx, indy
  real*4 xg, yg
  !
  logical edgeonlist, needcheck(maxInPc, 2), ngframe, useEdges, inLimit(2)
  real*4 xycur(2)
  integer*4 movedPiece(4), k, axisin, maxInside, ipcCross
  integer*4 ixframe, iyframe, ipc, ixfrm, iyfrm, limitLo, j, lookxfr, lookyfr
  integer*4 indinp, newedge, newpiece, iflo, listno, ixy, i, idSearch, limitHi
  real*4 xtmp, xframe, yframe, ytmp, xbak, ybak, distmin, dist
  real*4 xpcCross, ypcCross
  logical b3dxor
  !
  numpieces = 0
  ipcCross = 0
  numedges(1) = 0
  numedges(2) = 0
  idSearch = 1
  !
  ! get frame # that it is nominally in: actual one or nearby valid frame
  !
  xg = indx
  yg = indy
  if (secHasWarp) then
    call interpolateGrid(xg - 0.5, yg - 0.5, warpDx, warpDy, lmWarpX, nxWarp, nyWarp, &
        xWarpStrt, yWarpStrt, xWarpIntrv, yWarpIntrv, xtmp, ytmp)
    ! if (indy == 3845) print *,'ce', indx, xg, yg, xtmp, ytmp
    xg = xg + xtmp
    yg = yg + ytmp
  endif
  if (dogxforms) then
    xtmp = xg
    xg = ginv(1, 1) * xtmp + ginv(1, 2) * yg + ginv(1, 3)
    yg = ginv(2, 1) * xtmp + ginv(2, 2) * yg + ginv(2, 3)
  endif
  xframe = (xg - minxpiece - nxoverlap / 2) / (nxin - nxoverlap)
  yframe = (yg - minypiece - nyoverlap / 2) / (nyin - nyoverlap)
  ixframe = xframe + 1.                         !truncation gives proper frame
  iyframe = yframe + 1.
  ngframe = ixframe < 1 .or. ixframe > nxpieces .or. iyframe < 1 .or. &
      iyframe > nypieces .or. mappiece(min(nxpieces, max(1, ixframe)), &
      min(nypieces, max(1, iyframe))) == 0
  if (multng) then
    !
    ! if there are multineg h's (including piece shifting!), need to make
    ! sure point is actually in frame, but if frame no good, switch to
    ! nearest frame first
    if (ngframe) call findNearestPiece(ixframe, iyframe)
    ipc = mappiece(ixframe, iyframe)
    call positionInPiece(xg, yg, ipc, xbak, ybak)
    ngframe = xbak < 0 .or. xbak > nxin - 1 .or. ybak < 0 .or. ybak > nyin - 1

    if (ngframe) then
      !
      ! Use the error to switch to a nearby frame as center of search
      !
      xtmp = xbak + ixpclist(ipc)
      ytmp = ybak + iypclist(ipc)
      ixframe = (xtmp - minxpiece - nxoverlap / 2) / (nxin - nxoverlap) + 1.
      iyframe = (ytmp - minypiece - nyoverlap / 2) / (nyin - nyoverlap) + 1.
      !
      ! If still not in a frame, start in the nearest and expand the search
      !
      if (ixframe < 1 .or. ixframe > nxpieces .or. iyframe < 1 .or. &
          iyframe > nypieces .or. mappiece(min(nxpieces, max(1, ixframe)), &
          min(nypieces, max(1, iyframe))) == 0) then
        call findNearestPiece(ixframe, iyframe)
        idSearch = 2
      endif
    else
      !
      ! but even if frame is good, if in a corner, switch to the 9-piece
      ! search to start with most interior point
      !
      ngframe = (xbak < edgelonear(1) .and. ybak < edgelonear(2)) .or. &
          (xbak < edgelonear(1) .and. ybak > edgehinear(2)) .or. &
          (xbak > edgehinear(1) .and. ybak < edgelonear(2)) .or. &
          (xbak > edgehinear(1) .and. ybak > edgehinear(2))
    endif
  endif
  !
  ! if not a good frame, look in square of nine potential pieces,
  ! switch to the one that the point is a minimum distance from
  !
  if (ngframe) then
    distmin = 1.e10
    !
    ! continue loop to find most piece where point is most interior,
    ! not just to find the first one.  This should be run rarely so
    ! it is not a big drain
    !

    ixfrm = max(1, min(nxpieces, ixframe - idSearch))
    do while( ixfrm <= min(nxpieces, max(1, ixframe + idSearch)))
      iyfrm = max(1, min(nypieces, iyframe - idSearch))
      do while( iyfrm <= min(nypieces, max(1, iyframe + idSearch)))
        ipc = mappiece(ixfrm, iyfrm)
        if (ipc .ne. 0) then
          !
          ! get real coordinate in piece, adjusting for h if present
          !
          call positionInPiece(xg, yg, ipc, xtmp, ytmp)
          !
          ! distance is negative for a piece that point is actually in;
          ! it is negative of distance from nearest edge
          !
          dist = max(-xtmp, xtmp - (nxin - 1), -ytmp, ytmp - (nyin - 1))
          if (dist < distmin) then
            distmin = dist
            minxframe = ixfrm
            minyframe = iyfrm
          endif
        endif
        iyfrm = iyfrm + 1
      enddo
      ixfrm = ixfrm + 1
    enddo
    !
    ! return if no pieces in this loop.  This seems odd but there are
    ! weird edge effects if it returns just because it is not inside any
    if (distmin == 1.e10) return
    ixframe = minxframe
    iyframe = minyframe
  endif
  !
  ! initialize list of pieces with this piece # on it, then start looping
  ! over the pieces present in the list
  ! Keep track of min and max frame numbers on list
  numpieces = 1
  inpiece(1) = mappiece(ixframe, iyframe)
  indinp = 1
  needcheck(1, 1) = .true.
  needcheck(1, 2) = .true.
  minxframe = ixframe
  maxxframe = ixframe
  minyframe = iyframe
  maxyframe = iyframe
  inpxframe(1) = ixframe
  inpyframe(1) = iyframe
  call positionInPiece(xg, yg, inpiece(1), xinpiece(1), yinpiece(1))
  !
  do while (indinp <= numpieces)
    !
    ! come into this loop looking at a piece onlist; use true
    ! coordinates in piece to see if point is near/in an edge to another
    !
    ipc = inpiece(indinp)
    xycur(1) = xinpiece(indinp)
    xycur(2) = yinpiece(indinp)
    ! print *,'looking at', ipc, xycur(1), xycur(2)
    !
    ! check the x and y directions to see if point is near an edge
    !
    do ixy = 1, 2
      do iflo = 0, 1
        if (needcheck(indinp, ixy)) then
          newedge = 0
          if (iflo == 1 .and. xycur(ixy) < edgelonear(ixy) .and. &
              iedgelower(ipc, ixy) > 0) then
            newedge = iedgelower(ipc, ixy)
            newpiece = ipiecelower(newedge, ixy)
          endif
          if (iflo == 0 .and. xycur(ixy) > edgehinear(ixy) .and. &
              iedgeupper(ipc, ixy) > 0) then
            newedge = iedgeupper(ipc, ixy)
            newpiece = ipieceupper(newedge, ixy)
          endif
          !
          ! But if there are 4 pieces, make sure this picks up an edge
          ! between this piece and an existing piece
          if (newedge == 0 .and. numpieces >= 4) then
            lookxfr = inpxframe(indinp) + (2 - ixy) * (1 - 2 * iflo)
            lookyfr = inpyframe(indinp) + (ixy - 1) * (1 - 2 * iflo)
            do i = 1, numpieces
              if (inpxframe(i) == lookxfr .and. &
                  inpyframe(i) == lookyfr) then
                newpiece = inpiece(i)
                if (iflo == 0) newedge = iedgeupper(ipc, ixy)
                if (iflo == 1) newedge = iedgelower(ipc, ixy)
                exit
              endif
            enddo
          endif
          !
          ! if check picked up a new edge, see if edge is on list already
          !
          if (newedge .ne. 0) then
            edgeonlist = .false.
            do i = 1, numedges(ixy)
              edgeonlist = edgeonlist .or. (inedge(i, ixy) == newedge)
            enddo
            !
            ! if not, add it, and the implied piece, to list
            !
            if (.not.edgeonlist) then
              listno = 0
              do i = 1, numpieces
                if (newpiece == inpiece(i)) listno = i
              enddo
              if (listno == 0) then
                !
                ! but if adding a new piece, check for point actually in
                ! piece first
                !
                call positionInPiece(xg, yg, newpiece, xbak, ybak)
                if (xbak >= 0 .and. xbak <= nxin - 1 .and. &
                    ybak >= 0 .and. ybak <= nyin - 1) then
                  numpieces = numpieces + 1
                  inpiece(numpieces) = newpiece
                  xinpiece(numpieces) = xbak
                  yinpiece(numpieces) = ybak
                  !
                  ! Get new frame numbers and maintain min/max
                  if (ixy == 1) then
                    inpxframe(numpieces) = inpxframe(indinp) + 1 - 2 * iflo
                    inpyframe(numpieces) = inpyframe(indinp)
                  else
                    inpxframe(numpieces) = inpxframe(indinp)
                    inpyframe(numpieces) = inpyframe(indinp) + 1 - 2 * iflo
                  endif
                  minxframe = min(minxframe, inpxframe(numpieces))
                  minyframe = min(minyframe, inpyframe(numpieces))
                  maxxframe = max(maxxframe, inpxframe(numpieces))
                  maxyframe = max(maxyframe, inpyframe(numpieces))
                  !
                  ! If there are crossed limits to edges, then still need
                  ! to check this axis for this new piece
                  needcheck(numpieces, ixy) = &
                      edgelonear(ixy) >= edgehinear(ixy)
                  needcheck(numpieces, 3 - ixy) = .true.
                  listno = numpieces
                else if (numpieces == 1 .and. anyDisjoint( &
                    inpxframe(indinp), inpyframe(indinp))) then
                  !
                  ! If the point was NOT in this overlapping piece, and
                  ! any corners are disjoint, find cross-corner piece
                  ! Look across upper and lower edges on the other axis
                  ! from the rejected piece to see if point in other piece
                  newedge = iedgelower(newpiece, 3 - ixy)
                  if (newedge .ne. 0) then
                    call positionInPiece(xg, yg, &
                        ipiecelower(newedge, 3 - ixy), xbak, ybak)
                    if (xbak >= 0 .and. xbak <= nxin - 1 .and. &
                        ybak >= 0 .and. ybak <= nyin - 1) then
                      ipcCross = ipiecelower(newedge, 3 - ixy)
                      xpcCross = xbak
                      ypcCross = ybak
                    endif
                  endif
                  !
                  newedge = iedgeupper(newpiece, 3 - ixy)
                  if (newedge .ne. 0) then
                    call positionInPiece(xg, yg, &
                        ipieceupper(newedge, 3 - ixy), xbak, ybak)
                    if (xbak >= 0 .and. xbak <= nxin - 1 .and. &
                        ybak >= 0 .and. ybak <= nyin - 1) then
                      ipcCross = ipieceupper(newedge, 3 - ixy)
                      xpcCross = xbak
                      ypcCross = ybak
                    endif
                  endif
                endif
              endif
              !
              ! add edge to list only if legal piece found
              !
              if (listno > 0) then
                numedges(ixy) = numedges(ixy) + 1
                inedge(numedges(ixy), ixy) = newedge
                if (iflo == 0) then
                  inedupper(numedges(ixy), ixy) = listno
                  inedlower(numedges(ixy), ixy) = indinp
                else
                  inedupper(numedges(ixy), ixy) = indinp
                  inedlower(numedges(ixy), ixy) = listno
                endif
              endif
            endif
          endif
        endif
      enddo
    enddo
    indinp = indinp + 1
  enddo
  ! if (indy==998) write(*,'(2i5,i3,a,4i3,10i5)') indx, indy, numPieces, &
  ! ' pieces', minxframe, maxxframe, minyframe, maxyframe, (inPiece(i), i=1, numPieces)
  !
  ! If the pieces extend too far in either direction, we need to reduce the
  ! list to the ones where the point is most interior
  if (maxxframe > minxframe + 1 .or. maxyframe > minyframe + 1) then
    !
    ! first find axis where point is most interior
    maxInside = -10000000
    do i = 1, numPieces
      if (min(xinpiece(i), nxin - 1 - xinpiece(i)) > maxInside) then
        maxInside = min(xinpiece(i), nxin - 1 - xinpiece(i))
        axisin = 1
      endif
      if (min(yinpiece(i), nyin - 1 - yinpiece(i)) > maxInside) then
        maxInside = min(yinpiece(i), nyin - 1 - xinpiece(i))
        axisin = 2
      endif
    enddo
    ! if (indy==685) print *,'maxInside, axisin', maxInside, axisin
    !
    ! do the most interior axis first and find best set of frames for it
    ! Then do the other axis, with new constraint on first axis
    do ixy = 1, 2
      if (axisin == ixy) then
        call mostInteriorFrames(1, xinpiece, nxin, ixframe, iyframe)
        minxframe = ixframe
        maxxframe = min(ixframe + 1, maxxframe)
        if (ixy == 1) call constrainFrames(minxframe, maxxframe, &
            minyframe, maxyframe, inpxframe, inpyframe)
      else
        call mostInteriorFrames(2, yinpiece, nyin, ixframe, iyframe)
        minyframe = iyframe
        maxyframe = min(iyframe + 1, maxyframe)
        if (ixy == 1) call constrainFrames(minyframe, maxyframe, &
            minxframe, maxxframe, inpyframe, inpxframe)
      endif
      ! if (indx==37) print *,ixy, minxframe, maxxframe, minyframe, maxyframe
    enddo
    !
    ! Repack the frames, retaining only the ones within range and keeping
    ! track of former numbers
    j = 0
    do i = 1, numPieces
      if (inpxframe(i) >= minxframe .and. inpxframe(i) <= maxxframe .and. &
          inpyframe(i) >= minyframe .and. inpyframe(i) <= maxyframe) then
        j = j + 1
        inPiece(j) = inPiece(i)
        xinPiece(j) = xinPiece(i)
        yinPiece(j) = yinPiece(i)
        inpxFrame(j) = inpxFrame(i)
        inpyFrame(j) = inpyFrame(i)
        movedPiece(j) = i
      endif
    enddo
    numPieces = j
    !
    ! Repack edges too
    do ixy = 1, 2
      j = 0
      do i = 1, numedges(ixy)
        ixfrm = 0
        iyfrm = 0
        !
        ! Find the pieces that the lower and upper pieces became; if both
        ! exist, retain the edge and reassign the numbers
        do k = 1, numPieces
          if (inedlower(i, ixy) == movedPiece(k)) ixfrm = k
          if (inedupper(i, ixy) == movedPiece(k)) iyfrm = k
        enddo
        if (ixfrm > 0 .and. iyfrm > 0) then
          j = j + 1
          inedlower(j, ixy) = ixfrm
          inedupper(j, ixy) = iyfrm
          inedge(j, ixy) = inedge(i, ixy)
        endif
      enddo
      numedges(ixy) = j
    enddo
  endif
  !
  ! If there is one piece and a potential cross-piece, consider switching
  if (numPieces == 1 .and. ipcCross .ne. 0) then
    ixframe = (ixpclist(ipcCross) - minxpiece) / (nxin - nxoverlap) + 1
    iyframe = (iypclist(ipcCross) - minypiece) / (nyin - nyoverlap) + 1
    ixy = mapDisjoint(min(ixframe, inpxFrame(1)), min(iyframe, inpyFrame(1)))
    !
    ! Use cross piece if point is more interior along the other axis from
    ! the disjoint edges (ixy is 1, 2 for X, 3, 4 for Y)
    if (ixy > 0) then
      if ((ixy <= 2 .and. min(ypcCross, nyin - 1 - ypcCross) > &
          min(yinpiece(1), nyin - 1 - yinpiece(1))) .or. &
          (ixy > 2 .and. min(xpcCross, nxin - 1 - xpcCross) > &
          min(xinpiece(1), nxin - 1 - xinpiece(1)))) then
        inpiece(1) = ipcCross
        xinpiece(1) = xpcCross
        yinpiece(1) = ypcCross
        inpxFrame(1) = ixframe
        inpyFrame(1) = iyframe
      endif
    endif
  endif
  !
  ! If there are two pieces and the limit flag is set, find out if one
  ! should be thrown away
  if (limitData .and. numPieces == 2) then
    ixy = 1
    if (numedges(2) > 0) ixy = 2
    do i = 1, 2
      iflo = 1
      !
      ! if the piece is lower, get the limits on upper side
      if (i == inedlower(1, ixy)) iflo = 2
      call getDataLimits(inpiece(i), 3 - ixy, iflo, limitLo, limitHi)
      inLimit(i) = (ixy == 1 .and. yinPiece(i) >= limitLo .and. &
          yinPiece(i) <= limitHi) .or. (ixy == 2 .and. &
          xinPiece(i) >= limitLo .and. xinPiece(i) <= limitHi)
    enddo
    if (b3dxor(inLimit(1), inLimit(2))) then
      numPieces = 1
      numEdges(ixy) = 0
      if (inLimit(2)) then
        xinPiece(1) = xinPiece(2)
        yinPiece(1) = yinPiece(2)
        inPiece(1) = inPiece(2)
      endif
    endif
  endif
  !
  ! Replace edges with ones to use if called for
  if (useEdges) then
    do ixy = 1, 2
      do i = 1, numedges(ixy)
        call findEdgeToUse(inedge(i, ixy), ixy, newedge)
        if (newedge .ne. 0) inedge(i, ixy) = newedge
      enddo
    enddo
  endif
  return

CONTAINS
  !
  ! Find which set of pieces has the point most interior in a given
  ! direction
  subroutine mostInteriorFrames(ixy, xyinpiece, nxyin, ixbest, iybest)
    real*4 xyinpiece(*), closest, distin(2,2), dist1, dist2
    integer*4 nxyin, ixbest, iybest, ix, iy, ixy, maxOneCol, ixOneCol, iyOneCol
    maxInside = -100000
    maxOneCol = -100000
    ixbest = minxframe
    iybest = minyframe
    ixOneCol = minxframe
    iyOneCol = minyframe
    do iyfrm = minyframe, max(minyframe, maxyframe-1)
      do ixfrm = minxframe, max(minxframe, maxxframe-1)
        !
        ! Find the minimum distance to edge for frames that fit in this range
        distin = 1000000.
        do i = 1, numPieces
          ix = inpxframe(i) + 1 - ixfrm
          iy = inpyframe(i) + 1 - iyfrm
          if ((ix + 1) / 2 == 1 .and. (iy + 1) / 2 == 1)  distin(ix, iy) = &
              min(xyinpiece(i), nxyin - xyinpiece(i))
        enddo
        if (ixy == 1) then
          dist1 = min(distin(1, 1), distin(1, 2))
          dist2 = min(distin(2, 1), distin(2, 2))
        else
          dist1 = min(distin(1, 1), distin(2, 1))
          dist2 = min(distin(1, 2), distin(2, 2))
        endif
        !
        ! Keep track of which range maximizes this distance separately for
        ! ones with one column and two
        if (dist1 < 999999. .and. dist2 < 999999.) then
          closest = (dist1 + dist2) / 2.
          if (closest <= nxyin .and. closest > maxInside) then
            maxInside = closest
            ixbest = ixfrm
            iybest = iyfrm
          endif
        else
          closest = min(dist1, dist2)
          if (closest <= nxyin .and. closest > maxOneCol) then
            maxOneCol = closest
            ixOneCol = ixfrm
            iyOneCol = iyfrm
          endif
        endif
        ! if (indy==685.and.(indx+1) /2==798/2) print *,ixfrm, iyfrm, dist1, dist2, closest
      enddo
    enddo
    if (maxInside < 0 .and. maxOneCol > 0) then
      ixbest = ixOneCol
      iybest = iyOneCol
    endif
    return
  end subroutine mostInteriorFrames

  subroutine constrainFrames(minAxis, maxAxis, minOther, maxOther, &
      inpfAxis, inpfOther)
    integer*4 minAxis, maxAxis, minOther, maxOther, inpfAxis(*), inpfOther(*)
    i = minOther
    minOther = maxOther
    maxOther = i
    do i = 1, numPieces
      if (inpfAxis(i) >= minAxis .and. inpfAxis(i) <= maxAxis) then
        minOther = min(minOther, inpfOther(i))
        maxOther = max(maxOther, inpfOther(i))
      endif
    enddo
    return
  end subroutine constrainFrames

end subroutine countedges


! Computes position of global point xg, yg (after g transforms) in piece
! ipc, returns result in xinpc, yinpc
!
subroutine positionInPiece(xg, yg, ipc, xinpc, yinpc)
  use blendvars
  implicit none
  real*4 xg, yg, xinpc, yinpc, xtmp
  integer*4 ipc
  xinpc = xg - ixpclist(ipc)
  yinpc = yg - iypclist(ipc)
  if (multng) then
    xtmp = xinpc
    xinpc = hinv(1, 1, ipc) * xtmp + hinv(1, 2, ipc) * yinpc + hinv(1, 3, ipc)
    yinpc = hinv(2, 1, ipc) * xtmp + hinv(2, 2, ipc) * yinpc + hinv(2, 3, ipc)
  endif
  return
end subroutine positionInPiece


! initNearList sets up an ordered list of dx, dy values to nearby pieces
! up to the maximum distance in maxDistNear
!
subroutine initNearList()
  use blendvars
  implicit none
  integer*4 limx, limy, idx, idy, i, j, itmp
  real*4 temp

  limx = min(nxpieces - 1, maxDistNear)
  limy = min(nypieces - 1, maxDistNear)
  numPcNear = 0
  do idx = -limx, limx
    do idy = -limy, limy
      if (idx .ne. 0 .or. idy .ne. 0) then
        numPcNear = numPcNear + 1
        idxPcNear(numPcNear) = idx
        idyPcNear(numPcNear) = idy
        array(numPcNear) = idx**2 + idy**2
      endif
    enddo
  enddo
  do i = 1, numPcNear - 1
    do j = i + 1, numPcNear
      if (array(i) > array(j)) then
        temp = array(i)
        array(i) = array(j)
        array(j) = temp
        itmp = idxPcNear(i)
        idxPcNear(i) = idxPcNear(j)
        idxPcNear(j) = itmp
        itmp = idyPcNear(i)
        idyPcNear(i) = idyPcNear(j)
        idyPcNear(j) = itmp
      endif
    enddo
  enddo
  return
end subroutine initNearList


! findNearestPiece first modifies the frame numbers to be within the
! range for the montage, then checks whether this piece exists.  If not
! it switches to the nearest piece that does exist
!
subroutine findNearestPiece(ixFrame, iyFrame)
  use blendvars
  implicit none
  integer*4 ixFrame, iyFrame, i, ixnew, iynew
  ixFrame = max(1, min(nxpieces, ixFrame))
  iyFrame = max(1, min(nypieces, iyFrame))
  if (mappiece(ixFrame, iyFrame) .ne. 0) return
  do i = 1, numPcNear
    ixnew = ixFrame + idxPcNear(i)
    iynew = iyFrame + idyPcNear(i)
    if (ixnew >= 1 .and. ixnew <= nxpieces .and. iynew >= 1 .and. &
        iynew <= nypieces) then
      if (mappiece(ixnew, iynew) .ne. 0) then
        ixFrame = ixnew
        iyFrame = iynew
        return
      endif
    endif
  enddo
  return
end subroutine findNearestPiece

! DXYDGRINTERP takes a coordinate X1, Y1 in the lower piece of the edge
! at index INDEDG in the edge buffer, finds the coordinate within the
! edge function grid, uses bilinear interpolation to find the values of
! DX, DY and DDEN, and returns the coordinate in the upper piece, X2, Y2
! and the average density difference DDEN
!
subroutine dxydgrinterp(x1, y1, indedg, x2, y2, dden)
  !
  use blendvars
  implicit none
  real*4 x1, x2, y1, y2, dden
  integer*4 indedg
  !
  real*4 xingrid, yingrid, xgrid, ygrid, fx1, fx, c00, c10, c01, c11
  real*4 fy1, fy, dxinterp, dyinterp
  integer*4 ixg, iyg, ixg1, iyg1
  !
  ! find fractional coordinate within edge grid
  !
  xingrid = x1 - ixgrdstbf(indedg)
  yingrid = y1 - iygrdstbf(indedg)
  xgrid = 1. +(xingrid) / intxgrbf(indedg)
  ygrid = 1. +(yingrid) / intygrbf(indedg)
  !
  ! find all fractions and indices needed for bilinear interpolation
  !
  ixg = xgrid
  ixg = max(1, min(nxgrbf(indedg) - 1, ixg))
  iyg = ygrid
  iyg = max(1, min(nygrbf(indedg) - 1, iyg))
  fx1 = max(0., min(1., xgrid - ixg))             !NO EXTRAPOLATIONS ALLOWED
  fx = 1. -fx1
  ixg1 = ixg + 1
  fy1 = max(0., min(1., ygrid - iyg))
  fy = 1. -fy1
  iyg1 = iyg + 1
  c00 = fx * fy
  c10 = fx1 * fy
  c01 = fx * fy1
  c11 = fx1 * fy1
  !
  ! interpolate
  !
  dxinterp = c00 * dxgrbf(ixg, iyg, indedg) + c10 * dxgrbf(ixg1, iyg, indedg) &
      + c01 * dxgrbf(ixg, iyg1, indedg) + c11 * dxgrbf(ixg1, iyg1, indedg)
  dyinterp = c00 * dygrbf(ixg, iyg, indedg) + c10 * dygrbf(ixg1, iyg, indedg) &
      + c01 * dygrbf(ixg, iyg1, indedg) + c11 * dygrbf(ixg1, iyg1, indedg)
  dden = c00 * ddengrbf(ixg, iyg, indedg) + c10 * ddengrbf(ixg1, iyg, indedg) &
      + c01 * ddengrbf(ixg, iyg1, indedg) + c11 * ddengrbf(ixg1, iyg1, indedg)
  !
  x2 = xingrid + dxinterp + ixofsbf(indedg)
  y2 = yingrid + dyinterp + iyofsbf(indedg)
  return
end subroutine dxydgrinterp

subroutine crossvalue(xinlong, nxpieces, nypieces, nshort, nlong)
  logical xinlong
  if (xinlong) then
    nshort = nypieces
    nlong = nxpieces
  else
    nshort = nxpieces
    nlong = nypieces
  endif
  return
end subroutine crossvalue


! Does cross-correlation on an edge
!
subroutine xcorrEdge(arrLower, arrUpper, ixy, xDisplace, yDisplace, legacy, indentXC)
  use blendvars
  implicit none
  real*4 arrLower(*), arrUpper(*), xDisplace, yDisplace
  integer*4 indentXC, ixy
  logical legacy
  integer*4 nxyBox(2), ind0(2), ind1(2), iDisplace(2)
  real*4 ctf(8193), rDisplace(2)
  real*4 overFrac, delta, sdMin, delDenMin, ccc
  integer*4 indentSD, numIter, limStep, iyx, nxPad, nyPad, indentUse
  integer*4 ixDisplace, iyDisplace, i, numExtra(2), nbin, ierr, idum
  integer*4 numSmooth, nxSmooth, nySmooth, indPeak, ifEvalCCC, ifLegacy
  integer*4 numPixel, maxLongShift, niceFFTlimit, ind0Upper, ind1Upper
  integer*4 taperAtFill, extractWithBinning
  real*8 wallTime, wallStart
  external todfft, dumpEdge

  indentSD = 5                                !indent for sdsearch
  overFrac = 0.9                              !fraction of overlap to use
  numIter = 4                                   !iterations for sdsearch
  limStep = 10                                !limiting distance
  nbin = nbinXcorr
  numSmooth = 6
  wallStart = wallTime()
  !
  ! find size and limits of box in overlap zone to cut out
  !
  iyx = 3 - ixy
  ifLegacy = 0
  if (legacy) ifLegacy = 1
  ifEvalCCC = 0
  if (numXcorrPeaks > 1 .and. .not.legacy) ifEvalCCC = 1

  call montXCBasicSizes(ixy, nbin, indentXC, nxyzIn, noverlap, aspectMax, extraWidth, &
      padFrac, niceFFTlimit(), indentUse, nxyBox, numExtra, nxPad, nyPad, maxLongShift)
  call montXCIndsAndCTF(ixy, nxyzIn, noverlap, nxyBox, nbin, indentUse, numExtra, nxPad, &
      nyPad, numSmooth, sigma1, sigma2, radius1, radius2, ifEvalCCC, ind0, ind1,  &
      ind0Upper, ind1Upper, nxSmooth, nySmooth, ctf, delta)
  !print *,ixy, indentUse, nxybox(ixy), nxybox(iyx), ind0(iyx), ind1(iyx), &
  !    ind0(ixy), ind1(ixy)
  !print *,nxpad, nypad, nxSmooth, nySmooth

  if (nxyBox(1) * nxyBox(2) * nbin**2 > maxbsiz / 2 .or. nxPad * nyPad > idimc) &
      call exitError('CORRELATION ARRAYS WERE NOT MADE LARGE ENOUGH')
  !
  ! get the first image, lower piece
  ierr = extractWithBinning(arrLower, nxIn, ind0(1), ind1(1), ind0(2), ind1(2), nbin, &
      brray, 1, i, idum)
  if (ifillTreatment == 2) &
      ierr = taperAtFill(brray, nxyBox(1), nxyBox(2), 64, 0)
  !
  ! get the second image, upper piece
  ind0(ixy) = ind0Upper
  ind1(ixy) = ind1Upper
  !print *,'upper', ind0(ixy), ind1(ixy)
  ierr = extractWithBinning(arrUpper, nxIn, ind0(1), ind1(1), ind0(2), ind1(2), nbin, &
      brray(maxbsiz / 2 + 1), 1, i, idum)
  if (ifillTreatment == 2) &
      ierr = taperAtFill(brray, nxyBox(1), nxyBox(2), 64, 0)
  !
  ! Do the correlation
  call montXCorrEdge(brray, brray(maxbsiz / 2 + 1), nxyBox, nxyzIn, noverlap, nxSmooth, &
      nySmooth,nxPad, nyPad, xcray, xdray, xeray, numXcorrPeaks, ifLegacy, ctf, delta, &
      numExtra, nbin, ixy, maxLongShift, 0, xDisplace, yDisplace, ccc, todfft, dumpEdge, &
      0)

  if (legacy) return
  !
  ! the following is adopted largely from setgridchars
  !
  ixDisplace = nint(xDisplace)
  iyDisplace = nint(yDisplace)
  if (ixy == 1) then
    iDisplace(1) = nxIn - nxOverlap + ixDisplace
    iDisplace(2) = iyDisplace
  else
    iDisplace(1) = ixDisplace
    iDisplace(2) = nyIn - nyOverlap + iyDisplace
  endif
  iyx = 3 - ixy
  !
  ! get size of box, limit to size of overlap zone
  !
  nxyBox(ixy) = min(noverlap(ixy), nxyzIn(ixy) - iDisplace(ixy))
  nxyBox(ixy) = min(nxyBox(ixy) - indentSD * 2, nint(overFrac * nxyBox(ixy)))
  nxyBox(iyx) = nxyzIn(iyx) - abs(iDisplace(iyx))
  nxyBox(iyx) = min(nxyBox(iyx) - indentSD * 2, nint(aspectMax * nxyBox(ixy)))
  do i = 1, 2
    ind0(i) = (nxyzIn(i) + iDisplace(i) - nxyBox(i)) / 2
    ind1(i) = ind0(i) + nxyBox(i)
    rDisplace(i) = -iDisplace(i)
  enddo
  !
  ! integer scan is not needed, but to use it uncomment this
  !
  ! intscan=6
  ! call sdintscan(crray, drray, nxin, nyin, ind0(1), ind0(2), ind1(1), &
  ! ind1(2), -idispl(1) -intscan, -idispl(2) -intscan, &
  ! -idispl(1) +intscan, -idispl(2) +intscan, sdmin, ddenmin, &
  ! idxmin, idymin)
  ! rdispl(1) =idxmin
  ! rdispl(2) =idymin
  !
  
  call montBigSearch(arrLower, arrUpper, nxIn, nyIn, ind0(1), ind0(2), ind1(1), &
      ind1(2), rDisplace(1), rDisplace(2), sdMin, delDenMin, numIter, limStep)
  if (ixy == 1) then
    xDisplace = -rDisplace(1) - (nxIn - nxOverlap)
    yDisplace = -rDisplace(2)
  else
    xDisplace = -rDisplace(1)
    yDisplace = -rDisplace(2) - (nyIn - nyOverlap)
  endif
  ! write(*,'(a,f10.6)') 'time after big search', walltime() -wallstart
  ! write(*,'(2f8.2,2f8.2)') xpeak(indPeak), ypeak(indPeak), xdisp, yDisplace
  return
end subroutine xcorrEdge


! Solve for the shifts of all the pieces based on the displacements
! across their edges
!
subroutine find_best_shifts(dxgridmean, dygridmean, idimedge, idir, &
    izsect, h, nsum, bavg, bmax, aavg, amax)

  use blendvars
  implicit none
  integer*4 idir, izsect, nsum, idimedge
  real*4 bavg, bmax, aavg, amax
  !
  real*4 h(2,3,*)
  real*4 dxgridmean(idimedge,2), dygridmean(idimedge,2)
  integer*4 numInRow, newGroup, lowup
  !
  ! Set maxvar higher to get comparisons
  integer maxvar, maxGaussj
  parameter (maxGaussj = 10, maxvar = maxGaussj)
  real*4 a(maxvar,maxvar)
  !
  real*4 critmaxmove, critMoveDiff, wErrMean, wErrMax
  integer*4 intervalForTest, numAvgForTest, findPieceShifts, maxiter
  integer*4 nvar, ipc, ivar, m, ixy, iedge, neighpc, neighvar, ipclo, i, j, numGroups
  integer*4 numPrev, ndxy, nallvar, nextvar, nextCheck, numToCheck, igroup
  real*4 asum, bsum, xsum, ysum, bdist, adist, dxgroup, dygroup
  real*8 wallstart, walltime, wallAdj, wallGaussj
  integer*4 gaussj
!!$c
!!$c       variables for SVD
!!$      integer maxgels
!!$      parameter (maxgels = 12000)
!!$      real*8 daa(maxgels, maxgels), dbb(maxgels, 2), dwork(50*maxgels)
!!$      real*8 singval(maxgels), acond
!!$      common / sngval / daa, dbb, dwork, singval
  !
  ! The data coming in are the displacements of upper piece from
  ! being in alignment with the lower piece if idir = 1, or the
  ! shift needed to align upper piece with lower if idir = -1
  !
  ! build list of variables: ALL means all pieces that have an edge
  !
  nallvar = 0
  do ipc = 1, npclist
    if (izpclist(ipc) == izsect) then
      call xfunit(h(1, 1, ipc), 1.)
      call xfunit(hinv(1, 1, ipc), 1.)
      indvar(ipc) = 0
      if (iedgelower(ipc, 1) > 0 .or. &
          iedgelower(ipc, 2) > 0 .or. &
          iedgeupper(ipc, 1) > 0 .or. &
          iedgeupper(ipc, 2) > 0) then
        nallvar = nallvar + 1
        if (nallvar > limvar) then
          print *,'nallvar, limvar', nallvar, limvar
          call exitError( &
              'ARRAYS WERE NOT MADE LARGE ENOUGH FOR FIND_BEST_SHIFTS')
        endif
        iallVarpc(nallvar) = ipc
        ivarGroup(nallvar) = 0
        indvar(ipc) = nallvar
      endif
    endif
  enddo
  !
  ! Classify pieces into separate groups if any by following connections
  ! between them
  call sortVarsIntoGroups()
  !
  nsum = 0
  bsum = 0.
  bmax = 0.
  asum = 0.
  amax = 0.
  bavg = 0.
  aavg = 0.
  numPrev = 0
  if (nallVar == 0) call xfunit(h(1, 1, 1), 1.0)
  if (nallvar < 2) return
  !
  ! Loop on groups, set up to do fit for each group
  do igroup = 1, numGroups
    nvar = 0
    do ivar = 1, nallVar
      if (ivarGroup(ivar) == igroup) then
        nvar = nvar + 1
        ivarpc(nvar) = iallVarpc(ivar)
        indvar(ivarpc(nvar)) = nvar
      endif
    enddo
    ! print *,nvar
    if (nvar > 1) then
      !
      ! build matrix of simultaneous equations for minimization solution
      ! by matrix inversion or SVD
      do ivar = 1, nvar - 1
        ipc = ivarpc(ivar)
        do m = 1, nvar - 1
          rowTmp(m) = 0.
          bb(1, ivar) = 0.
          bb(2, ivar) = 0.
        enddo
        !
        do ixy = 1, 2
          if (includeEdge(1, ipc, ixy, iedge)) then
            rowTmp(ivar) = rowTmp(ivar) + 1
            neighpc = ipiecelower(iedge, ixy)
            neighvar = indvar(neighpc)
            !
            ! for a regular neighbor, enter a -1 in its term; but for the
            ! last variable being eliminated, enter a +1 for ALL other
            ! variables instead
            !
            if (neighvar .ne. nvar) then
              rowTmp(neighvar) = rowTmp(neighvar) - 1
            else
              do m = 1, nvar - 1
                rowTmp(m) = rowTmp(m) + 1
              enddo
            endif
            !
            ! when this piece is an upper piece, subtract displacements
            ! from constant term
            !
            bb(1, ivar) = bb(1, ivar) - idir * dxgridmean(iedge, ixy)
            bb(2, ivar) = bb(2, ivar) - idir * dygridmean(iedge, ixy)
          endif
          !
          if (includeEdge(2, ipc, ixy, iedge)) then
            rowTmp(ivar) = rowTmp(ivar) + 1
            neighpc = ipieceupper(iedge, ixy)
            neighvar = indvar(neighpc)
            if (neighvar .ne. nvar) then
              rowTmp(neighvar) = rowTmp(neighvar) - 1
            else
              do m = 1, nvar - 1
                rowTmp(m) = rowTmp(m) + 1
              enddo
            endif
            !
            ! when a lower piece, add displacements to constant terms
            !
            bb(1, ivar) = bb(1, ivar) + idir * dxgridmean(iedge, ixy)
            bb(2, ivar) = bb(2, ivar) + idir * dygridmean(iedge, ixy)
          endif
        enddo
!!$c
!!$c         LOAD FOR SVD
!!$            do m = 1, nvar - 1
!!$              daa(ivar, m) = rowTmp(m)
!!$            enddo
!!$            dbb(ivar, 1) = bb(1, ivar)
!!$            dbb(ivar, 2) = bb(2, ivar)
        !
        ! Load the row data in a if below maxvar
        if (nvar <= maxvar) then
          do m = 1, nvar - 1
            a(m, ivar) = rowTmp(m)
          enddo
        endif
      enddo
!!$c
!!$c       Solve SVD
!!$          wallstart = walltime()
!!$          acond = -1.
!!$          call dgelss(nvar - 1, nvar - 1, 2, daa, maxgels, dbb, maxgels, singval, &
!!$              acond, m, dwork, 50 * maxgels, ixy)
!!$          write(*,'(a,i5,a,i5,a,f10.4)') 'dgelss  rank', m, &
!!$              '  info', ixy, '  time', walltime() - wallstart

      !
      ! Solve by iteration first
      if (nvar > maxGaussj) then
        critMaxMove = 1.e-4
        numAvgForTest = 10
        intervalForTest = 100
        critMoveDiff = 1.e-6
        maxiter = 100 + nvar * 10
        wallstart = walltime()
        wallstart = walltime()
        if (findPieceShifts(ivarpc, nvar, indvar, ixpclist, iypclist, &
            dxgridmean, dygridmean, idir, ipiecelower, ipieceupper, &
            ifskipEdge, limedge, dxyvar, limvar, iedgelower, iedgeupper, &
            limnpc, fpsWork, 1, 0, 2, robustCrit, critMaxMove, &
            critMoveDiff, maxiter, numAvgForTest, intervalForTest, i, &
            wErrMean, wErrMax) &
            .ne. 0) call exitError('CALLING findPieceShifts')
        wallAdj = walltime() - wallstart
        ! write(*,'(i6,a,f8.4,a,f15.6)') i, ' iterations, time', wallAdj

!!$            ysum = 0.
!!$            xsum = 0.
!!$            do ivar = 1, nvar - 1
!!$              ipc = ivarpc(ivar)
!!$              do ixy = 1, 2
!!$                xsum = max(xsum, abs(dxyvar(ivar, ixy) - dbb(ivar, ixy)))
!!$                do j = 1, 2
!!$                  if (includeEdge(2, ipc, j, iedge)) then
!!$                    numInRow = 0
!!$                    neighpc = ipieceupper(iedge, j)
!!$                    neighvar = indvar(neighpc)
!!$                    if (neighvar .ne. nvar) ysum = max(ysum, abs( &
!!$                        (dbb(ivar, ixy) - dbb(neighvar, ixy)) &
!!$                        - (dxyvar(ivar, ixy) - dxyvar(neighvar, ixy))))
!!$                  endif
!!$                enddo
!!$              enddo
!!$            enddo
!!$            write(*,'(a,f12.7,a,f12.7)') 'Shift adj - SVD max position diff' &
!!$                , xsum, '  edge diff', ysum
      endif
      !
      ! solve the equations with gaussj if within range

      !
      ! write(*,'(9i5)') (ivarpc(i), i=1, nvar)
      ! write(*,'(8f7.1)') ((a(j, i), i=1, nvar-1), j=1, nvar-1)
      ! write(*,'(8f9.2)') ((bb(j, i), i=1, nvar-1), j=1, 2)
      if (nvar <= maxvar) then
        wallstart = walltime()
        i = gaussj(a, nvar - 1, maxvar, bb, 2, 2)
        if (i > 0) call exitError( &
            'SINGULAR MATRIX WHEN SOLVING LINEAR EQUATIONS IN GAUSSJ')
        if (i < 0) call exitError( &
            'TOO MANY VARIABLES TO SOLVE LINEAR EQUATIONS WITH GAUSSJ')
        wallGaussj = walltime() - wallstart
        ! write(*,'(8f9.2)') ((bb(j, i), i=1, nvar-1), j=1, 2)
      endif
      !
      ! Use the iteration solution, do comparisons if both were done
      if (nvar > maxGaussj) then
        do j = 1, 2
          xsum = 0.
          do i = 1, nvar - 1
            if (nvar <= maxvar) xsum = &
                max(xsum, abs(bb(j, i) - dxyvar(i, j)))
            bb(j, i) = dxyvar(i, j)
          enddo
          if (nvar <= maxvar) write(*,'(a,i2,a,f15.7)'), 'axis', j, &
              ' max shift adj - gaussj difference', xsum
          ! write(*,'(8f9.2)') (bb(j, i), i=1, nvar-1)
        enddo
        if (nvar <= maxvar) write(*,'(a,f12.6, a, f12.6)') 'gaussj time', &
            wallGaussj, '   shift adj time', wallAdj
      endif
    endif
    !
    ! take the b values as dx and dy; compute the
    ! sum to get the shift for the final piece
    xsum = 0.
    ysum = 0.
    do i = 1, nvar - 1
      h(1, 3, ivarpc(i)) = bb(1, i)
      h(2, 3, ivarpc(i)) = bb(2, i)
      xsum = xsum + bb(1, i)
      ysum = ysum + bb(2, i)
    enddo
    h(1, 3, ivarpc(nvar)) = -xsum
    h(2, 3, ivarpc(nvar)) = -ysum
    !
    ! For multiple groups, find the mean displacement between this group
    ! and all previous ones and adjust shifts to make that mean be zero
    ! but to keep the overall mean zero
    if (igroup > 1) then
      dxgroup = 0.
      dygroup = 0.
      ndxy = 0
      !
      ! rebuild index to variables
      do ivar = 1, nallvar
        indvar(iallVarpc(ivar)) = ivar
      enddo
      !
      ! Loop on pieces in this group, and for each edge to a piece in a
      ! lower group, add up the displacement across the edge
      do ivar = 1, nvar
        ipc = ivarpc(ivar)
        do lowup = 1, 2
          do ixy = 1, 2
            if (.not. includeEdge(lowup, ipc, ixy, iedge)) then
              if (iedge > 0) then
                if (lowup == 1) then
                  ipclo = ipiecelower(iedge, ixy)
                else
                  ipclo = ipieceupper(iedge, ixy)
                endif
                if (ivarGroup(indvar(ipclo)) < igroup) then
                  dxgroup = dxgroup + h(1, 3, ipc) - h(1, 3, ipclo)
                  dygroup = dygroup + h(2, 3, ipc) - h(2, 3, ipclo)
                  ndxy = ndxy + 1
                endif
              endif
            endif
          enddo
        enddo
      enddo
      !
      ! Adjust positions by a weighted fraction of the mean displacement
      if (ndxy > 0) then
        dxgroup = dxgroup / ndxy
        dygroup = dygroup / ndxy
        ! print *,'Shift relative to previous group', dxgroup, dygroup
        ! print *,'adjusting this', -(dygroup * numPrev) / (nvar +numPrev), &
        ! '   this', (dygroup * nvar) / (nvar + numPrev)
        do ivar = 1, nallvar
          ipc = iallVarpc(ivar)
          if (ivarGroup(ivar) < igroup) then
            h(1, 3, ipc) = h(1, 3, ipc) + (dxgroup * nvar) / (nvar + numPrev)
            h(2, 3, ipc) = h(2, 3, ipc) + (dygroup * nvar) / (nvar + numPrev)
          else if (ivarGroup(ivar) == igroup) then
            h(1, 3, ipc) = h(1, 3, ipc) - (dxgroup * numPrev) / (nvar + numPrev)
            h(2, 3, ipc) = h(2, 3, ipc) - (dygroup * numPrev) / (nvar + numPrev)
          endif
        enddo
      endif
    endif
    numPrev = numPrev + nvar
  enddo
  !
  ! compute and return the results
  !
  do ivar = 1, nallvar
    ipc = iallVarpc(ivar)
    call xfinvert(h(1, 1, ipc), hinv(1, 1, ipc))
    ! write(*,'(i7,2i4,2f8.1)') ipc, 1+(ixpclist(ipc) -minxpiece) /(nxin-nxoverlap), &
    ! 1+(iypclist(ipc) -minypiece) /(nyin-nyoverlap), h(1, 3, ipc), h(2, 3, ipc)
    do ixy = 1, 2
      if (includeEdge(1, ipc, ixy, iedge)) then
        ipclo = ipiecelower(iedge, ixy)
        bdist = sqrt(dxgridmean(iedge, ixy)**2 + &
            dygridmean(iedge, ixy)**2)
        bsum = bsum + bdist
        adist = sqrt((idir * dxgridmean(iedge, ixy) + h(1, 3, ipc) - h(1, 3, ipclo)) &
            **2 + (idir * dygridmean(iedge, ixy) + h(2, 3, ipc) - h(2, 3, ipclo))**2)
        asum = asum + adist
        amax = max(amax, adist)
        bmax = max(bmax, bdist)
        nsum = nsum + 1
        ! write(*,'(i7,2i4,f8.2)') ipc, ixy-1+(ixpclist(ipc) -minxpiece) /(nxin-nxoverlap), &
        ! 2-ixy+(iypclist(ipc) -minypiece) /(nyin-nyoverlap), adist
      endif
    enddo
  enddo
  bavg = bsum / nsum
  aavg = asum / nsum
  ! write(*,'(i3,a,2f6.2,a,2f6.2)') nsum, ' edges, mean, max '// &
  ! 'displacement before:', bsum/nsum, bmax, ', after:', asum/nsum, amax
  !
  return

CONTAINS

  ! Test for whether an edge should be included, lowup = 1/2 for lower/
  ! upper edge, kpc is piece index, kxy is x/y direction, the edge number
  ! if any is returned in ked
  !
  logical function includeEdge(lowup, kpc, kxy, ked)
    integer*4 lowup, kpc, kxy, ked
    includeEdge = .false.
    if (lowup == 1) then
      ked = iedgelower(kpc, kxy)
    else
      ked = iedgeupper(kpc, kxy)
    endif
    if (ked == 0) return
    if (ifskipEdge(ked, kxy) > 1) return
    includeEdge = .true.
    return
  end function includeEdge


  ! Check for whether to add a piece to the current group
  !
  subroutine checkGroup(kpc)
    integer*4 kpc
    ! print *,'checking group of ', kpc, ivarGroup(indvar(kpc))
    if (ivarGroup(indvar(kpc)) > 0) return
    ivarGroup(indvar(kpc)) = numGroups
    numToCheck = numToCheck + 1
    listCheck(numToCheck) = kpc
    return
  end subroutine checkGroup


  ! Classify pieces into separate groups if any by following connections
  ! between them
  !
  subroutine sortVarsIntoGroups()
    numGroups = 0
    nextvar = 1
    do while (nextvar <= nallvar)
      !
      ! Look for next piece that is unclassified
      do while (nextvar <= nallvar)
        if (ivarGroup(nextvar) == 0) exit
        nextvar = nextvar + 1
      enddo
      if (nextvar > nallvar) exit
      !
      ! Start a new group and initialize check list
      numGroups = numGroups + 1
      ivarGroup(nextvar) = numGroups
      nextCheck = 1
      numToCheck = 1
      listCheck(1) = iallVarpc(nextvar)
      !
      ! Loop on the check list, for next piece, check its 4 edges and
      ! assign and add to check list other pieces not in group yet
      do while (nextCheck <= numToCheck)
        ipc = listCheck(nextCheck)
        do ixy = 1, 2
          if (includeEdge(1, ipc, ixy, iedge)) &
              call checkGroup(ipiecelower(iedge, ixy))
          if (includeEdge(2, ipc, ixy, iedge)) &
              call checkGroup(ipieceupper(iedge, ixy))
        enddo
        nextCheck = nextCheck + 1
      enddo
      nextvar = nextvar + 1
    enddo
    ! print *,'groups:', numGroups
    !
    ! Have to make sure the groups are going to be done in a connected order
    ! if there are more than 2
    do igroup = 1, numGroups - 2
      newGroup = 0
      do ivar = 1, nallVar
        !
        ! Look at pieces in all the groups already done
        if (ivarGroup(ivar) <= igroup) then
          ipc = iallVarpc(ivar)
          !
          ! Look for an excluded edge to a piece in a higher group
          do ixy = 1, 2
            if (.not. includeEdge(1, ipc, ixy, iedge)) then
              if (iedge > 0) then
                ipclo = ipiecelower(iedge, ixy)
                m = ivarGroup(indvar(ipclo))
                if (m > igroup) newGroup = m
              endif
            endif
            if (.not. includeEdge(2, ipc, ixy, iedge)) then
              if (iedge > 0) then
                ipclo = ipieceupper(iedge, ixy)
                m = ivarGroup(indvar(ipclo))
                if (m > igroup) newGroup = m
              endif
            endif
          enddo
          !
          ! Swap the group numbers of that group and the next one
          if (newGroup > 0) then
            ! print *,'swapping groups', igroup + 1, newGroup
            do i = 1, nallVar
              if (ivarGroup(i) == igroup + 1) then
                ivarGroup(i) = newGroup
              else if (ivarGroup(i) == newGroup) then
                ivarGroup(i) = igroup + 1
              endif
            enddo
            exit
          endif
        endif
      enddo
    enddo
    return
  end subroutine sortVarsIntoGroups

end subroutine find_best_shifts


subroutine findBestGradient(dxgridmean, dygridmean, idimEdge, idir, izsect, &
    gradnew, rotnew)
  use blendvars
  implicit none
  integer*4 idir, izsect, idimEdge
  real*4 gradnew, rotnew
  real*4 dxgridmean(idimedge,2), dygridmean(idimedge,2)
  !
  ! Stuff for amoeba: ftol2 and ptol2 are used the FIRST time
  !
  integer nvar
  parameter (nvar = 2)
  real*4 pp(nvar+1,nvar+1), yy(nvar+1), ptol(nvar), da(nvar)
  real*4 ptol1, ftol1, ptol2, ftol2, delfac, var(nvar)
  data da/0.5, 0.2/
  integer*4 jmin, iter, i
  external gradfunc
  !
  integer*4 nedg, ifTrace, nTrial, izedge
  real*4 errMin
  common / funccom / nedg, ifTrace, nTrial, izedge, errMin
  !
  integer*4 ipc, ipclo, ixy, iedge

  ifTrace = 0
  ptol1 = 1.e-5
  ftol1 = 1.e-5
  ptol2 = 1.e-3
  ftol2 = 1.e-3
  delfac = 2.

  izedge = izsect
  nedg = 0
  do ipc = 1, npclist
    if (izpclist(ipc) == izsect .and. (iedgelower(ipc, 1) > 0 .or. &
        iedgelower(ipc, 2) > 0)) then
      do ixy = 1, 2
        iedge = iedgelower(ipc, ixy)
        if (iedge > 0) then
          nedg = nedg + 1
          ipclo = ipiecelower(iedge, ixy)
          !
          ! Get the residual at the edge: the amount the upper piece
          ! is still displaced away from alignment with lower
          !
          dxedge(iedge, ixy) = idir * dxgridmean(iedge, ixy)
          dyedge(iedge, ixy) = idir * dygridmean(iedge, ixy)
          !
          ! Get center for gradient in each piece
          !
          if (focusAdjusted) then
            gradXcenLo(nedg) = nxin / 2.
            gradYcenLo(nedg) = nyin / 2.
            gradXcenHi(nedg) = nxin / 2.
            gradYcenHi(nedg) = nyin / 2.
          else
            gradXcenLo(nedg) = (minxpiece + nxpieces * (nxin - nxoverlap) &
                + nxoverlap) / 2. - ixpclist(ipclo)
            gradYcenLo(nedg) = (minypiece + nypieces * (nyin - nyoverlap) &
                + nyoverlap) / 2. - iypclist(ipclo)
            gradXcenHi(nedg) = (minxpiece + nxpieces * (nxin - nxoverlap) &
                + nxoverlap) / 2. - ixpclist(ipc)
            gradYcenHi(nedg) = (minypiece + nypieces * (nyin - nyoverlap) &
                + nyoverlap) / 2. - iypclist(ipc)
          endif
          !
          ! Get center point of overlap zone in each piece
          !
          if (ixy == 1) then
            overXcenLo(nedg) = nxin - nxoverlap / 2
            overXcenHi(nedg) = nxoverlap / 2
            overYcenLo(nedg) = nyin / 2
            overYcenHi(nedg) = nyin / 2
          else
            overYcenLo(nedg) = nyin - nyoverlap / 2
            overYcenHi(nedg) = nyoverlap / 2
            overXcenLo(nedg) = nxin / 2
            overXcenHi(nedg) = nxin / 2
          endif

          ! write(*, '(i2,2i4,2f6.1, 8f7.0)') ixy, ipclo, ipc, dxedge(iedge, ixy), &
          ! dyedge(iedge, ixy), overXcenLo(nedg), overYcenLo(nedg), &
          ! overXcenHi(nedg), overYcenHi(nedg), gradXcenLo(nedg), &
          ! gradYcenLo(nedg), gradXcenHi(nedg), gradYcenHi(nedg)
        endif
      enddo
    endif
  enddo

  !
  ! set up for minimization
  !
  errMin = 1.e30
  nTrial = 0
  var(1) = 0.
  var(2) = 0.
  call amoebaInit(pp, yy, nvar + 1, nvar, delfac, ptol2, var, da, &
      gradfunc, ptol)
  call amoeba(pp, yy, nvar + 1, nvar, ftol2, gradfunc, iter, ptol, jmin)
  !
  ! per Press et al. recommendation, just restart at current location
  !
  do i = 1, nvar
    var(i) = pp(jmin, i)
  enddo
  call amoebaInit(pp, yy, nvar + 1, nvar, delfac / 4., ptol1, var, da, &
      gradfunc, ptol)
  call amoeba(pp, yy, nvar + 1, nvar, ftol1, gradfunc, iter, ptol, jmin)
  !
  ! recover result
  !
  do i = 1, nvar
    var(i) = pp(jmin, i)
  enddo
  call gradfunc(var, errMin)
  write(*,73) var(1), var(2), errMin
73 format(' Implied incremental gradient:',2f9.4,'  mean error:',f10.4)
  call flush(6)
  gradnew = var(1)
  rotnew = var(2)
  return
end subroutine findBestGradient


subroutine gradfunc(p, funcErr)
  use blendvars
  implicit none
  real*4 p(*), funcErr

  integer*4 nedg, ifTrace, nTrial, izedge
  real*4 errMin
  common / funccom / nedg, ifTrace, nTrial, izedge, errMin
  !
  real*4 tiltang, errSum, dxlo, dxhi, dylo, dyhi
  real*4 bmean, bmax, aftmean, aftmax
  integer*4 ied, ixy, iedge, ipc
  character*1 starout

  tiltang = tiltAngles(min(ilistz, numAngles))
  nTrial = nTrial + 1
  errSum = 0.

  ied = 0
  do ipc = 1, npclist
    if (izpclist(ipc) == izedge .and. (iedgelower(ipc, 1) > 0 .or. &
        iedgelower(ipc, 2) > 0)) then
      do ixy = 1, 2
        iedge = iedgelower(ipc, ixy)
        if (iedge > 0) then
          ied = ied + 1
          call magGradientShift(overXcenLo(ied), overYcenLo(ied), &
              nxin, nyin, gradXcenLo(ied), gradYcenLo(ied), &
              pixelMagGrad, axisRot, tiltang, p(1), p(2), dxlo, dylo)
          call magGradientShift(overXcenHi(ied), overYcenHi(ied), &
              nxin, nyin, gradXcenHi(ied), gradYcenHi(ied), &
              pixelMagGrad, axisRot, tiltang, p(1), p(2), dxhi, dyhi)
          !
          ! Point moves by negative of shift to get to undistorted image
          ! Displacement changes by negative of upper shift and positive
          ! of  lower shift
          !
          dxadj(iedge, ixy) = dxedge(iedge, ixy) + dxlo - dxhi
          dyadj(iedge, ixy) = dyedge(iedge, ixy) + dylo - dyhi
        endif
      enddo
    endif
  enddo

  call find_best_shifts(dxadj, dyadj, limedge, 1, izedge, htmp, &
      iedge, bmean, bmax, aftmean, aftmax)

  funcErr = aftmean
  if (ifTrace > 0) then
    starout = ' '
    if (funcErr < errMin) then
      starout = '*'
      errMin = funcErr
    endif
    if (iftrace > 1 .or. starout == '*') &
        write(*,72) starout, nTrial, funcErr, (p(ied), ied = 1, 2)
72  format(1x,a1,i4,f15.5, 2f9.4)
    call flush(6)
  endif
  !
  ! Amoeba flails on hard zeros so add something (no longer, 6/17/06)
  !
  funcErr = funcErr
  return
end subroutine gradfunc


! findEdgeToUse looks up an edge to use in place of the given edge
! if Z limits are defined on the edges to use.  It returns the original
! edge number if no limits have been set, or if no substitte edge should
! be used.  It returns 0 if a different Z level should be used but no
! edge exists at it.
!
subroutine findEdgeToUse(iedge, ixy, iuse)
  use blendvars
  implicit none
  integer*4 iedge, iuse, ixy, izlow, izhigh, ipclo, ixfrm, iyfrm, izuse, i
  iuse = iedge
  if (numUseEdge == 0 .and. izUseDefLow < 0) return
  !
  ! Find frame number of piece then look up frame in the use list
  ipclo = ipiecelower(iedge, ixy)
  ixfrm = 1 + (ixpclist(ipclo) - minxpiece) / (nxin - nxoverlap)
  iyfrm = 1 + (iypclist(ipclo) - minypiece) / (nyin - nyoverlap)
  izlow = -1
  do i = 1, numUseEdge
    if (ixFrmUseEdge(i) == ixfrm .and. iyFrmUseEdge(i) == iyfrm .and. &
        ixy == ixyUseEdge(i)) then
      izlow = izLowUse(i)
      izhigh = izHighUse(i)
    endif
  enddo
  !
  ! If frame not found, use the default or skip out if none
  ! Then find out if there is another Z to use or not
  if (izlow < 0) then
    if (izUseDefLow < 0) return
    izlow = izUseDefLow
    izhigh = izUseDefHigh
  endif
  if (izpclist(ipclo) < izlow) then
    izuse = izlow
  else if (izpclist(ipclo) > izhigh) then
    izuse = izhigh
  else
    return
  endif
  !
  ! Seek these coordinates on that Z with an edge above
  do i = 1, npclist
    if (izpclist(i) == izuse .and. ixpclist(i) == ixpclist(ipclo) .and. &
        iypclist(i) == iypclist(ipclo) .and. iedgeUpper(i, ixy) > 0) &
        then
      iuse = iedgeUpper(i, ixy)
      return
    endif
  enddo
  iuse = 0
  return
end subroutine findEdgeToUse


! getDataLimits returns the limits of actual data (as opposed to gray
! fill data) in the dimension given by IXY (1 for X, 2 for Y) and for
! an edge on the side given by LOHI (1 for lower edge, 2 for upper edge) .
! It returns coordinates numbered from 0 in limitLo and limitHi, looking
! them up in the limDataLo and limDataHi arrays or finding the values and
! storing them in those arrays when done.
!
subroutine getDataLimits(ipc, ixy, lohi, limitLo, limitHi)
  use blendvars
  implicit none
  integer*4 ipc, ixy, limitLo, limitHi, ind, incPix, incLine, numPix, line
  integer*4 lineEnd, lineStart, idir, i, lineBase, lineGood, maxSame, numSame
  integer*4 limInd, lohi, iPixStr, iPixEnd
  real*4 value
  !
  limInd = limDataInd(ipc)
  if (limInd <= 0) then
    limitLo = 0
    limitHi = nxin - 1
    if (ixy == 2) limitHi = nyin - 1
    return
  endif
  if (limDataLo(limInd, ixy, lohi) >= 0 .and. &
      limDataHi(limInd, ixy, lohi) >= 0) then
    limitLo = limDataLo(limInd, ixy, lohi)
    limitHi = limDataHi(limInd, ixy, lohi)
    return
  endif
  call shuffler(ipc, ind)

  lineEnd = 0
  if (ixy == 1) then
    incPix = nxin
    incLine = 1
    numPix = min(nyin, (3 * nyoverlap) / 2)
    iPixEnd = nyin - 1
    lineStart = nxin - 1
  else
    incPix = 1
    incLine = nxin
    numPix = min(nxin, (3 * nxoverlap) / 2)
    iPixEnd = nxin - 1
    lineStart = nyin - 1
  endif
  if (lohi == 1) then
    iPixStr = 0
    iPixEnd = numPix - 1
  else
    iPixStr = iPixEnd + 1 - numPix
  endif
  maxSame = numPix / 4
  do idir = -1, 1, 2
    !
    ! First find out if first line has a common value
    i = iPixStr + 1
    lineBase = ind + lineStart * incLine
    value = array(lineBase + iPixStr * incPix)
    do while (i <= iPixEnd)
      if (array(lineBase + i * incPix) .ne. value) exit
      i = i + 1
    enddo
    if (i <= iPixEnd) then
      lineGood = lineStart
    else
      !
      ! If it made it to the end, next search for line with less than
      ! maximum number of pixels at this value
      lineGood = -1
      line = lineStart + idir
      do while (lineGood < 0 .and. idir * (line - lineEnd) < 0)
        i = iPixStr
        numSame = 0
        lineBase = ind + line * incLine
        do while (i <= iPixEnd .and. numSame < maxSame)
          if (array(lineBase + i * incPix) == value) numSame = numSame + 1
          i = i + 1
        enddo
        if (numSame < maxSame) lineGood = line
        line = line + idir
      enddo
    endif
    !
    ! If go to end, set it to next to last line, then save limit
    if (lineGood < 0) lineGood = lineEnd - idir

    if (idir == -1) then
      limitHi = lineGood
    else
      limitLo = lineGood
    endif
    lineEnd = lineStart
    lineStart = 0
  enddo
  limDataLo(limInd, ixy, lohi) = limitLo
  limDataHi(limInd, ixy, lohi) = limitHi
  return
end subroutine getDataLimits


! iwrBinned writes the data in ARRAY to unit IUNIT, with binning
! given by IBINNING, using BRRAY as a scratch line.  The data in ARRAY
! are NY lines of length NX.  The output will consist of NYOUT lines
! of length NXOUT.  IXST and IYST are possible starting pixels, where
! negative numbers are used to indicate non-existent pixels.
! The routine also maintains the min and max density in DMIN and DMAX
! and adds to the sum of densities in DSUM8.
!
subroutine iwrBinned(iunit, array, brray, nx, nxout, ixst, ny, nyout, &
    iyst, iBinning, dmin, dmax, dsum8)
  implicit none
  integer*4 nx, nxout, ixst, ny, nyout, iyst, iBinning, iunit
  real*4 array(nx,ny), brray(*), sum, dmin, dmax
  integer*4 ix, iy, jx, jy, jxbase, jybase, nsum
  real*8 dsum8
  !
  if (iBinning == 1) then
    do iy = 1, nyout
      do ix = 1, nxout
        dmin = min(dmin, array(ix, iy))
        dmax = max(dmax, array(ix, iy))
        dsum8 = dsum8 + array(ix, iy)
      enddo
      call parWrtLin(iunit, array(1, iy))
    enddo
    return
  endif

  do iy = 1, nyout
    jybase = (iy - 1) * iBinning + iyst
    do ix = 1, nxout
      jxbase = (ix - 1) * iBinning + ixst
      sum = 0.
      if (ix == 1 .or. ix == nxout .or. iy == 1 .or. iy == nyout) &
          then
        nsum = 0
        do jy = max(1, jybase + 1), min(ny , jybase + iBinning)
          do jx = max(1, jxbase + 1), min(nx, jxbase + iBinning)
            sum = sum + array(jx, jy)
            nsum = nsum + 1
          enddo
        enddo
        brray(ix) = sum / nsum
      else
        do jy = jybase + 1, jybase + iBinning
          do jx = jxbase + 1, jxbase + iBinning
            sum = sum + array(jx, jy)
          enddo
        enddo
        brray(ix) = sum / iBinning**2
      endif
      dmin = min(dmin, brray(ix))
      dmax = max(dmax, brray(ix))
      dsum8 = dsum8 + brray(ix)
    enddo
    call parWrtLin(iunit, brray)
  enddo
  return
end subroutine iwrBinned


! getExtraIndents computes a border or indent for correlations and
! edge functions when there are distortion corrections, based upon the
! maximum  distortion vector along any edge of the two pieces, whose
! numbers are given by IPCLOW and IPCUP.  IXY specifies the direction
! of the overlap zone.  The indents in X and Y are returned in
! delIndent.
!
subroutine getExtraIndents(ipclow, ipcup, ixy, delIndent)
  use blendvars
  implicit none
  integer*4 ipclow, ipcup, ixy, ix, memlow, memup, iy
  real*4 delIndent(2)
  !
  delIndent(1) = 0.
  delIndent(2) = 0.
  if (.not. doFields) return
  memlow = memIndex(ipclow)
  memup = memIndex(ipcup)
  !
  ! The undistorted image moves in the direction opposite to the
  ! field vector, so positive vectors at the right edge of the lower
  ! piece move the border in to left and require more indent in
  ! short direction.
  !
  if (ixy == 1) then
    do iy = 1, nyField
      delIndent(1) = max(delIndent(1), fieldDx(nxField, iy, memlow), &
          - fieldDx(1, iy, memup))
    enddo
    delIndent(2) = max(0., -fieldDy(nxField, 1, memlow), &
        - fieldDy(1, 1, memup), fieldDy(nxField, nyField, memlow), &
        fieldDy(1, nyField, memup))
  else
    do ix = 1, nxField
      delIndent(1) = max(delIndent(1), fieldDy(ix, nyField, memlow), &
          - fieldDy(ix, 1, memup))
    enddo
    delIndent(2) = max(0., -fieldDx(1, nyField, memlow), &
        - fieldDx(1, 1, memup), fieldDx(nxField, nyField, memlow), &
        fieldDx(nxField, 1, memup))
  endif
  return
end subroutine getExtraIndents


! Reads a model of edges to exclude, analyzes it and marks edges
! for exclusion.  FILNAM has the name of the model, EDGEDISPX and
! EDGEDISPY are the edge displacements (zero unless read from file) with
! direst dimension IDIMEDGE, findEFforAdjusted is a flag that edge
! functions should be found for excluded edges that have nonzero
! displacements, and numSkip is returned with the number of edges skipped
!
subroutine readExclusionModel(filnam, edgedispx, edgedispy, idimedge, &
    ifUseAdjusted, mapAllPc, nxmap, nymap, minzpc, numSkip)
  use blendvars
  implicit none
  include 'smallmodel.inc90'
  character*(*) filnam
  integer*4 nxmap, nymap, minzpc, mapAllPc(nxmap,nymap,*)
  integer*4 ixy, ied, ipc, iobj, ipt, ipnt, numEFonly, idimedge, numSkip
  integer*4 ixpc, iypc, lenx, leny, numNear, ifUseAdjusted, ixframe, iyframe
  integer*4 ixright, iytop, izpc
  real*4 edgedispx(idimedge,2), edgedispy(idimedge,2), vertex(4,2)
  logical exist, readSmallMod, inside
  !
  exist = readSmallMod(filnam)
  if (.not.exist) call exitError('READING EDGE EXCLUSION MODEL FILE')
  call scaleModelToImage(1, 0)
  numEFonly = 0
  numSkip = 0
  numnear = 0
  !
  ! Loop on the edges, make a quadrangle for area nearest to each
  do ixy = 1, 2
    do ied = 1, nedge(ixy)
      ipc = ipieceupper(ied, ixy)
      ixpc = ixpclist(ipc)
      iypc = iypclist(ipc)
      izpc = izpclist(ipc) + 1 - minzpc
      lenx = nxin - nxoverlap
      leny = nyin - nyoverlap
      ixframe = 1 + (ixpc - minxpiece) / lenx
      iyframe = 1 + (iypc - minypiece) / leny
      ixright = ixpc + lenx
      iytop = iypc + leny
      !
      ! WARNING: This code tracks how 3dmod lays out data when displaying
      ! a montage: it copies each entire piece into the buffer in Z order
      ! A piece extends farther in one direction if it is either the last
      ! piece or it is laid down after the next in that direction
      if (ixframe == nxpieces) then
        ixright = ixpc + nxin
      else if (mapAllPc(ixframe + 1, iyframe, izpc) < ipc) then
        ixright = ixpc + nxin
      endif
      if (iyframe == nypieces) then
        iytop = iypc + nyin
      else if (mapAllPc(ixframe, iyframe + 1, izpc) < ipc) then
        iytop = iypc + nyin
      endif
      !
      ! A piece starts farther in one direction if the piece before it
      ! is laid down after this piece
      if (ixframe > 1) then
        if (mapAllPc(ixframe-1, iyframe, izpc) > ipc) &
            ixpc = ixpc + nxoverlap
      endif
      if (iyframe > 1) then
        if (mapAllPc(ixframe, iyframe-1, izpc) > ipc) &
            iypc = iypc + nyoverlap
      endif
      lenx = ixright - ixpc
      leny = iytop - iypc
      !
      ! This puts the vertices on the edge at its visible ends and one
      ! vertex in the middle of this visible piece; the other vertex is
      ! just a fixed distance into the piece below
      vertex(1, 1) = ixpc + lenx / 2.
      vertex(1, 2) = iypc + leny / 2.
      vertex(2, 1) = ixpc
      vertex(2, 2) = iypc
      if (ixy == 1) then
        vertex(3, 1) = ixpc - (nxin - nxoverlap) / 2.
        vertex(3, 2) = vertex(1, 2)
        vertex(4, 1) = vertex(2, 1)
        vertex(4, 2) = vertex(2, 2) + leny
      else
        vertex(3, 1) = vertex(1, 1)
        vertex(3, 2) = iypc - (nyin - nyoverlap) / 2.
        vertex(4, 1) = vertex(2, 1) + lenx
        vertex(4, 2) = vertex(2, 2)
      endif
      ! print *,'edge', ixy, ixpc, iypc, int(vertex(4, 1)), int(vertex(4, 2))
      !
      ! test each model point for being inside this quadrangle
      do iobj = 1, max_mod_obj
        do ipt = 1, npt_in_obj(iobj)
          ipnt = abs(object(ipt + ibase_obj(iobj)))
          if (nint(p_coord(3, ipnt)) == izpclist(ipc)) then
            if (inside(vertex(1, 1), vertex(1, 2), 4, p_coord(1, ipnt), &
                p_coord(2, ipnt)) .and. ifskipEdge(ied, ixy) == 0) then
              ! write(*,'(3i3,10f6.0)') ied, ixy, ipc, p_coord(1, ipnt), &
              ! p_coord(2, ipnt), ((vertex(i, j), j=1, 2), i=1, 4)
              numNear = numNear + 1
              ifskipEdge(ied, ixy) = 2
              if ((edgedispx(ied, ixy) .ne. 0. .or. &
                  edgedispy(ied, ixy) .ne. 0.) .and. ifUseAdjusted > 0) then
                ifskipEdge(ied, ixy) = 1
                if (ifUseAdjusted > 1) ifskipEdge(ied, ixy) = 0
                numEFonly = numEFonly + 1
              endif
              if (ifskipEdge(ied, ixy) > 0) numSkip = numSkip + 1
            endif
          endif
        enddo
      enddo
    enddo
  enddo
  if (numNear > 0) then
    if (numEFonly == 0) then
      write(*,'(/,i7,a,/,a)') numSkip, ' edges will be given zero edge '// &
          'functions and', '   excluded when solving for shifts'
    else if (ifUseAdjusted > 1) then
      write(*,'(/,i7,a,i7,a,/,a)') numSkip, &
          ' edges will be given zero edge functions but', numEFonly, &
          ' others have', '  non-zero displacements and edge functions '// &
          'will be found for them'
    else
      write(*,'(/,i7,a,i7,a,/,a)') numSkip, &
          ' edges will be given zero edge functions but', numEFonly, &
          ' of them have', '  non-zero displacements and will be '// &
          'included when solving for shifts'
    endif
  endif
  if (n_point > numNear) &
      write(*,'(/,a,i7,a,i7,a)') 'WARNING: only', numSkip, ' of the ', n_point, &
      ' model points were near an edge'
  return
end subroutine readExclusionModel


! dumpEdge writes a padded image or correlation to a file.  The
! image in CRRAY, NXDIM is the X array dimension and NXPAD and NYPAD
! are image sizes.  IXY is 1 for X or 2 for Y edge.
!
subroutine dumpedge(crray, nxdim, nxpad, nypad, ixy, ifcorr)
  use blendvars
  implicit none
  integer maxline
  parameter (maxline = 4096)
  integer*4 nxdim, nxpad, nypad, ixy, ifcorr
  real*4 crray(nxdim,nypad)
  real*4 title(20), scale, dmin, dmax, dmt, bline(maxline)
  integer*4 kxyz(3), ix, iy
  !
  if (ifDumpXY(ixy) < 0 .or. nxpad > maxline) return
  if (ifDumpXY(ixy) == 0) then
    nzOutXY(ixy) = 0
    nxOutXY(ixy) = nxpad
    nyOutXY(ixy) = nypad
    kxyz(1) = nxpad
    kxyz(2) = nypad
    kxyz(3) = 0
    call icrhdr(2 + ixy, kxyz, kxyz, 2, title, 0)
    ifDumpXY(ixy) = 1
  endif
  call imposn(2 + ixy, nzOutXY(ixy), 0)
  nzOutXY(ixy) = nzOutXY(ixy) + 1
  call ialsiz_sam_cel(2 + ixy, nxOutXY(ixy), nyOutXY(ixy), nzOutXY(ixy))

  call iclden(crray, nxdim, nypad, 1, nxpad, 1, nypad, dmin, dmax, dmt)
  scale = 255. / (dmax - dmin)
  do iy = 1, nyOutXY(ixy)
    if (iy <= nypad) then
      if (ifcorr == 0) then
        do ix = 1, nxOutXY(ixy)
          bline(ix) = scale * (crray(min(ix, nxpad), iy) - dmin)
        enddo
      else
        do ix = 1, nxOutXY(ixy)
          bline(ix) = scale * (crray(min(mod(ix + nxpad / 2 - 1, nxpad) + 1, nxpad), &
              min(mod(iy + nypad / 2 - 1, nypad) + 1, nypad)) - dmin)
        enddo
      endif
    endif
    call iwrlin(2 + ixy, bline)
  enddo
  if (ifcorr .ne. 0) then
    ix = (ixpclist(ipcBelowEdge) - minxpiece) / (nxin - nxoverlap)
    iy = (iypclist(ipcBelowEdge) - minypiece) / (nyin - nyoverlap)
    if (ixy == 1) ix = iy * (nxpieces - 1) + ix + 1
    if (ixy == 2) ix = ix * (nypieces - 1) + iy + 1
    write (*,'(1x,a,i4,a,i5)') char(ixy + ichar('W')) //' edge', ix, &
        ' corr at Z', nzOutXY(ixy)
    call flush(6)
  endif
  !
  call iwrhdr(2 + ixy, title, -1, 0., 255., 128.)
  return
end subroutine dumpedge

