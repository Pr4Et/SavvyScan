! FIXBOUNDARIES
!
! Rewrites data near chunk boundaries after direct writing in parallel
! to an output file. Takes two arguments, the name of the main image
! file and the name of the boundary info file.  Uses the entries in the
! info file to rewrite all of the data in the boundary files into the
! main file.
!
! $Id$
!
program fixboundaries
  implicit none
  integer IDIM
  parameter (IDIM = 1000000)
  real*4 array(IDIM)
  character*320 mainFile, infoFile
  integer*4 nxyz(3), mxyz(3), nxyz2(3), mode, i, j, ifile, iz, ierr
  integer*4 izSecs(2), lineStart(2), numFiles, ifAllSec, linesGuard
  real*4 dmin, dmax, dmean
  integer*4 parWrtInitialize, parWrtGetRegion, parWrtProperties, iiuFileType
  !
  call getinout(2, mainFile, infoFile)
  call setExitPrefix('ERROR: fixboundaries - ')
  call imopen(1, mainFile, 'old')
  call irdhdr(1, nxyz, mxyz, mode, dmin, dmax, dmean)
  iz = nxyz(1)
  if (iiuFileType(1) == 5) iz = -iz
  ierr = parWrtInitialize(infoFile, 2, iz, nxyz(2), nxyz(3))
  if (ierr .ne. 0) call exitError('Initializing from the parallel write information file')
  ierr = parWrtProperties(ifAllSec, linesGuard, numFiles)
  if (ierr > 0) call exitError( &
      'Getting properties from the parallel write information file')
  if (ierr < 0) then
    print *, 'Fixboundaries: Nothing to do for an HDF file'
    call exit(0)
  endif
  !
  ! Loop on all the files
  do ifile = 1, numFiles
    if (parWrtGetRegion(ifile, infoFile, izSecs, lineStart) .ne. 0) &
        call exitError('Getting information for one region')
    call imopen(2, infoFile, 'ro')
    call irdhdr(2, nxyz2, mxyz, mode, dmin, dmax, dmean)
    if (ifAllSec .ne. 0) then
      !
      ! If writing all sections, loop on all sections
      do iz = 1, nxyz(3)
        do j = 1, 2
          if (lineStart(j) >= 0) then
            call iiuSetPosition(2, iz * 2 + j - 3, 0)
            call iiuSetPosition(1, iz - 1, lineStart(j))
            do i = 1, linesGuard
              call irdlin(2, array, *99)
              call iiuWriteLines(1, array, 1)
            enddo
          endif
        enddo
      enddo
    else
      !
      ! Otherwise just write the two sections
      do j = 1, 2
        if (izSecs(j) >= 0) then
          call iiuSetPosition(2, j - 1, 0)
          call iiuSetPosition(1, izSecs(j), lineStart(j))
          do i = 1, linesGuard
            call irdlin(2, array, *99)
            call iiuWriteLines(1, array, 1)
          enddo
        endif
      enddo
    endif
    call iiuClose(2)
  enddo
  call iiuClose(1)
  call exit(0)
99 call exitError('Reading from guard file')
end program
