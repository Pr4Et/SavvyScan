! * * * * * * * ASSEMBLEVOL * * * * * * * *
!
! Assemblevol will assemble a single MRC file from separate files
! that form an array of subvolumes in X, Y, and Z.  In effect, it can
! take montaged images and compose them into single images, just as
! Reducemont can, but its advantage is that it can take the images
! from multiple files.  Its primary use is for reassembling a tomogram
! after it has been chopped into pieces, using the coordinates output
! by Tomopieces.
!
! See man page for details.
!
! David Mastronarde, 3/1/01
!
! $Id$
!
program assemblevol
  implicit none
  integer*4 nx, ny, nz
  integer*4 nxyz(3), mxyz(3), nxyzst(3), nxyz2(3), mxyz2(3)
  equivalence (nx, nxyz(1)), (ny, nxyz(2)), (nz, nxyz(3))
  real*4 title(20), cell2(6), delta(3), tilt(3), originalTilt(3)
  real*4, allocatable :: array(:), brray(:)
  !
  integer*4, allocatable :: ixLow(:), ixHigh(:), iyLow(:), iyHigh(:)
  integer*4, allocatable :: izLow(:), izHigh(:)
  logical*4, allocatable :: xDefined(:), yDefined(:), zDefined(:)
  character*320 outFile
  character*320, allocatable :: files(:)
  character*9 curDate
  character*8 curTime
  character*80 titlech
  !
  integer*4 numXfiles, numYfiles, numZfiles, nxOut, nyOut, nzOut, maxX, maxY, i, ifile
  integer*4 numFiles, limRange, maxFilesToOpen, ix, iy, iz
  integer(kind = 8) idimIn, idimOut
  real*4 dmin2, dmax2, dmean2, dmin, dmax, dmean, tmin, tmax, tmean, tempMin
  real*4 xOrigin, yOrigin, zOrigin
  integer*4 mode, modeFirst, kti, indZfile, iunit, iyOffset, ixOffset, nyBox, nxBox
  integer*4 indXfile, indYfile, layerFile, ierr, numOptFiles
  integer*4 numXrangeIn, numYrangeIn, numZrangeIn, maxZ
  logical openLayer
  data nxyzst/0, 0, 0/
  !
  logical pipInput
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetTwoIntegers, PipGetString, PipNumberOfEntries, PipGetNonOptionArg
  integer*4 PipGetInteger
  !
  ! fallbacks from ../../manpages/autodoc2man -3 2  assemblevol
  !
  integer numOptions
  parameter (numOptions = 10)
  character*(40 * numOptions) options(1)
  options(1) = &
      'input:InputFile:FNM:@output:OutputFile:FN:@'// &
      'xextract:StartEndToExtractInX:IPM:@yextract:StartEndToExtractInY:IPM:@'// &
      'zextract:StartEndToExtractInZ:IPM:@nxfiles:NumberOfFilesInX:I:@'// &
      'nyfiles:NumberOfFilesInY:I:@nzfiles:NumberOfFilesInZ:I:@'// &
      'param:ParameterFile:PF:@help:usage:B:'
  !
  maxFilesToOpen = 256
  !
  ! Pip startup: set error, parse options, do help output
  !
  call PipReadOrParseOptions(options, numOptions, 'assemblevol', &
      'ERROR: ASSEMBLEVOL - ', .true., 2, 2, 1, numOptArg, numNonOptArg)
  pipInput = numOptArg + numNonOptArg > 0

  if (pipInput) then
    !
    ! Get output file from either place, adjust non-option count down if it is that
    if (PipGetString('OutputFile', outFile) .ne. 0) then
      if (numNonOptArg == 0) call exitError('Output filename must be entered '// &
          'either with option or as last non-option argument')
      ierr = PipGetNonOptionArg(numNonOptArg, outFile)
      numNonOptArg = numNonOptArg - 1
    endif
    !
    ! Input files are option entries and remaining non-opt entries
    ierr = PipNumberOfEntries('InputFile', numOptFiles)
    numFiles = numOptFiles + numNonOptArg
    call getPipFileNumbers(numXfiles, numXrangeIn, 'X')
    call getPipFileNumbers(numYfiles, numYrangeIn, 'Y')
    call getPipFileNumbers(numZfiles, numZrangeIn, 'Z')
    if (numXfiles * numYfiles * numZfiles .ne. numFiles) call exitError &
        ('The number of input files entered does not equal the product of the number '// &
        'in X, Y and Z implied by the other entries')
  else
    
    write(*,'(1x,a,$)') 'Output file for assembled volume: '
    read(5, 101) outFile
101 format(a)
    !
    write(*,'(1x,a,$)') 'Numbers of input files in X, Y, and Z: '
    read(5,*) numXfiles, numYfiles, numZfiles
    if (numXfiles < 0 .or. numYfiles < 0 .or. numZfiles < 0) &
        call exitError('Number of files must be positive')
  endif

  limRange = max(numXfiles, numYfiles, numZfiles)
  numFiles = numXfiles * numYfiles * numZfiles
  allocate(ixLow(limRange), ixHigh(limRange), iyLow(limRange), iyHigh(limRange),  &
      izLow(limRange), izHigh(limRange), files(numFiles), xDefined(limRange), &
      yDefined(limRange), zDefined(limRange), stat = indXfile)
  call memoryError(indXfile, 'arrays for ranges or filenames')

  openLayer = numXfiles * numYfiles <= maxFilesToOpen
  !
  if (.not. pipInput) print *,'Enter the starting and ending index coordinates '// &
      'for the pixels to extract', ' from the files at successive positions '// &
      'in each dimension', ' (0,0 for full extent or -1,-1 for first pixel only)'
  nxOut = 0
  nyOut = 0
  nzOut = 0
  maxX = 0
  maxY = 0
  maxZ = 0
  ixLow(:) = 0
  iyLow(:) = 0
  izLow(:) = 0
  ixHigh(:) = 0
  iyHigh(:) = 0
  izHigh(:) = 0
  xDefined(:) = .false.
  yDefined(:) = .false.
  zDefined(:) = .false.

  ! Get the lower and upper limits for each position on each axis
  call getCheckLowHighLimits(numXfiles, numXrangeIn, pipInput, nxOut, ixLow, ixHigh, &
      xDefined, maxX, 'X')
  call getCheckLowHighLimits(numYfiles, numYrangeIn, pipInput, nyOut, iyLow, iyHigh, &
      yDefined, maxY, 'Y')
  call getCheckLowHighLimits(numZfiles, numZrangeIn, pipInput, nzOut, izLow, izHigh, &
      zDefined, maxZ, 'Z')
  !
  if (.not. pipInput) print *,'Enter the input file names at successive positions '// &
      'in X, then Y, then Z'
  ifile = 1
  call iiuAltPrint(0)
  do iz = 1, numZfiles
    do iy = 1, numYfiles
      do ix = 1, numXfiles
        if (pipInput) then
          if (ifile <= numOptFiles) then
            ierr = PipGetString('InputFile', files(ifile))
          else
            ierr = PipGetNonOptionArg(ifile - numOptFiles, files(ifile))
          endif
        else
          write(*,'(1x,a,3i4,a,$)') 'Name of file at', ix, iy, iz, ': '
          read(5, 101) files(ifile)
        endif
        call imopen(2, files(ifile), 'ro')
        call irdhdr(2, nxyz, mxyz, mode, dmin2, dmax2, dmean2)
        if (ix == 1 .and. iy == 1 .and. iz == 1) then
          call iiuRetDelta(2, delta)
          call iiuRetOrigin(2, xOrigin, yOrigin, zOrigin)
          call iiuRetTilt(2, tilt)
          call iiuRetTiltOrig(2, originalTilt)
        endif
        if (ifile == 1) modeFirst = mode
        if (mode .ne. modeFirst) call exitError( 'Mode mismatch for this file')
        !
        ! Collect coordinates if they are not defined yet
        call setUndefinedLimits(ix, nx, nxOut, ixLow, ixHigh, xDefined, maxX)
        call setUndefinedLimits(iy, ny, nyOut, iyLow, iyHigh, yDefined, maxY)
        call setUndefinedLimits(iz, nz, nzOut, izLow, izHigh, zDefined, maxZ)
        if (ixHigh(ix) >= nx .or. iyHigh(iy) >= ny .or. izHigh(iz) >= nz) &
            call exitError('Upper coordinate too high for this file')
        call iiuClose(2)
        ifile = ifile + 1
      enddo
    enddo
  enddo
  call PipDone()
  !
  idimOut = int(nxOut, kind = 8) * nyOut + 10
  idimIn = int(maxX, kind = 8) * maxY + 10
  allocate(array(idimIn), brray(idimOut), stat = indXfile)
  call memoryError(indXfile, 'arrays for image data')
  !
  call imopen(1, outFile, 'NEW')
  nxyz2(1) = nxOut
  nxyz2(2) = nyOut
  nxyz2(3) = nzOut
  mxyz2(1) = nxOut
  mxyz2(2) = nyOut
  mxyz2(3) = nzOut
  cell2(1) = nxOut * delta(1)
  cell2(2) = nyOut * delta(2)
  cell2(3) = nzOut * delta(3)
  cell2(4) = 90.
  cell2(5) = 90.
  cell2(6) = 90.
  xOrigin = xOrigin - ixLow(1) * delta(1)
  yOrigin = yOrigin - iyLow(1) * delta(2)
  zOrigin = zOrigin - izLow(1) * delta(3)
  !
  call time(curTime)
  call b3ddate(curDate)
  write(titlech, 301) curDate, curTime
  read(titlech, '(20a4)') (title(kti), kti = 1, 20)
301 format('ASSEMBLEVOL: Reassemble a volume from pieces',t57,a9,2x,a8)
  call iiuCreateHeader(1, nxyz2, mxyz2, mode, title, 0)
  call iiuAltCell(1, cell2)
  call iiuAltOrigin(1, xOrigin, yOrigin, zOrigin)
  call iiuAltTilt(1, tilt)
  call iiuAltTiltOrig(1, originalTilt)
  dmin = 1.e30
  dmax = -1.e30
  tmean = 0.
  !
  ifile = 1
  do indZfile = 1, numZfiles
    !
    ! open the files on this layer if possible
    !
    if (openLayer) then
      do i = 1, numYfiles * numXfiles
        call imopen(i + 1, files(ifile), 'ro')
        call irdhdr(i + 1, nxyz, mxyz, mode, dmin2, dmax2, dmean2)
        !
        ! DNM 7/31/02: transfer labels from first file
        !
        if (ifile == 1) call iiuTransLabels(1, i + 1)
        ifile = ifile + 1
      enddo
    endif
    !
    ! loop on the sections to be composed
    !
    do iz = izLow(indZfile), izHigh(indZfile)
      iunit = 2
      iyOffset = 0
      !
      ! loop on the files in X and Y
      !
      layerFile = ifile
      do indYfile = 1, numYfiles
        ixOffset = 0
        nyBox = iyHigh(indYfile) + 1 - iyLow(indYfile)
        do indXfile = 1, numXfiles
          !
          ! open files one at a time if necessary
          !
          if (.not.openLayer) then
            call imopen(2, files(layerFile), 'ro')
            call irdhdr(2, nxyz, mxyz, mode, dmin2, dmax2, dmean2)
            if (indXfile == 1 .and. indYfile == 1 .and. indZfile == 1 .and. &
                iz == izLow(indZfile)) call iiuTransLabels(1, 2)
            layerFile = layerFile + 1
          endif
          !
          ! read the section, and insert into big array
          !
          call iiuSetPosition(iunit, iz, 0)
          nxBox = ixHigh(indXfile) + 1 - ixLow(indXfile)
          call irdpas(iunit, array, nxBox, nyBox, ixLow(indXfile), ixHigh(indXfile), &
              iyLow(indYfile), iyHigh(indYfile), *99)
          call insert_array(array, nxBox, nyBox, brray, nxOut, nyOut, ixOffset, &
              iyOffset)
          ixOffset = ixOffset + nxBox
          if (openLayer) then
            iunit = iunit + 1
          else
            call iiuClose(2)
          endif
        enddo
        iyOffset = iyOffset + nyBox
      enddo
      !
      ! section done, get density and write it
      !
      call iclden(brray, nxOut, nyOut, 1, nxOut, 1, nyOut, tmin, tmax, tempMin)
      dmin = min(dmin, tmin)
      dmax = max(dmax, tmax)
      tmean = tmean + tempMin
      call iiuWriteSection(1, brray)
    enddo
    !
    ! close layer files if opened; otherwise set file number for next
    ! layer
    !
    if (openLayer) then
      do i = 1, numYfiles * numXfiles
        call iiuClose(i + 1)
      enddo
    else
      ifile = layerFile
    endif
  enddo
  dmean = tmean / nzOut
  call iiuWriteHeader(1, title, 1, dmin, dmax, dmean)
  call iiuClose(1)
  write(*,'(/,i6,a)') numXfiles * numYfiles * numZfiles, ' files reassembled'
  call exit(0)
99 call exitError('Reading file')
end program assemblevol


! Copy an array with an offset into a larger array
!
subroutine insert_array(array, nxBox, nyBox, brray, nxOut, nyOut, ixOffset, &
    iyOffset)
  implicit none
  integer*4 nxBox, nyBox, nxOut, nyOut, ixOffset, iyOffset, ix, iy
  real*4 array(nxBox,nyBox), brray(nxOut,nyOut)
  do iy = 1, nyBox
    do ix = 1, nxBox
      brray(ix + ixOffset, iy + iyOffset) = array(ix, iy)
    enddo
  enddo
  return
end subroutine insert_array


! Get the number of files for an axis one way or another
!
subroutine getPipFileNumbers(numXfiles, numXrangeIn, axis)
  implicit none
  integer*4 numXfiles, numXrangeIn, ierr
  character*1 axis
  integer*4 PipNumberOfEntries, PipGetInteger
  ierr = PipNumberOfEntries('StartEndToExtractIn'//axis, numXrangeIn)
  if (PipGetInteger('NumberOfFilesIn'//axis, numXfiles) == 0) then
    if (numXrangeIn > 0) call exitError('You cannot enter both StartEndToExtractIn'// &
        axis//' and NumberOfFilesIn'//axis)
  else
    numXfiles = max(1, numXrangeIn)
  endif
  return
end subroutine getPipFileNumbers


! Get the low and high limits for each position on an axis
!
subroutine getCheckLowHighLimits(numXfiles, numXrangeIn, pipInput, nxOut, ixLow, &
    ixHigh, xDefined, maxX, axis)
  implicit none 
  character*1 axis
  logical*4 pipInput, xDefined(*)
  integer*4 nxOut, ixLow(*), ixHigh(*), maxX, ind, numXfiles, numXrangeIn, ierr
  integer*4 PipGetTwoIntegers
  do ind = 1, numXfiles
    if (pipInput) then
      if (numXrangeIn > 0) &
          ierr = PipGetTwoIntegers('StartEndToExtractIn'//axis, ixLow(ind), ixHigh(ind))
    else
      write(*,'(1x,a,a,i3,a,a,a,$)') axis, ' coordinates for files at position #', ind,  &
          ' in ', axis, ': '
      read(5,*) ixLow(ind), ixHigh(ind)
    endif
    !
    ! Special case of -1, -1: mark axis as defined and set to 0,0
    ! Otherwise axis is defined if either is nonzero
    if (ixLow(ind) == -1 .and. ixHigh(ind) == -1) then
      xDefined(ind) = .true.
      ixLow(ind) = 0
      ixHigh(ind) = 0
    else
      xDefined(ind) = ixLow(ind) > 0 .or. ixHigh(ind) > 0
    endif
    if (ixLow(ind) < 0 .or. ixHigh(ind) < ixLow(ind)) call exitError( &
        'Illegal '//axis//' coordinate less than zero or out of order')
    if (xDefined(ind)) then
      nxOut = nxOut + ixHigh(ind) + 1 - ixLow(ind)
      maxX = max(maxX, ixHigh(ind) + 1 - ixLow(ind))
    endif
  enddo
  return
end subroutine getCheckLowHighLimits


! If the limits are still undefined for a position, set them to full range
!
subroutine setUndefinedLimits(ix, nx, nxOut, ixLow, ixHigh, xDefined, maxX)
  implicit none 
  logical*4 pipInput, xDefined(*)
  integer*4 nxOut, ixLow(*), ixHigh(*), maxX, ix, nx
  if (xDefined(ix)) return
  ixLow(ix) = 0
  ixHigh(ix) = nx - 1
  nxOut = nxOut + nx
  maxX = max(maxX, nx)
  xDefined(ix) = .true.
  return
end subroutine setUndefinedLimits
