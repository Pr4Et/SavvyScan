!*************EXTRACTPIECES**********************************************
!
! EXTRACTPIECES will extract piece coordinates from the header of
! an image file, if they are present, and produce a file with those
! coordinates (a piece list file) .
!
! David Mastronarde, 1/2/00
!
program extractpieces
  implicit none
  integer*4 nxyz(3), mxyz(3), nz, maxPiece, maxExtra, numExtraBytes, nbytes, iflags
  real*4, allocatable :: array(:)
  integer*4, allocatable :: ixPiece(:), iyPiece(:), izPiece(:)
  integer*4 mode, numPieces, i, ierr, ifAddMdoc, indAdoc, iTypeAdoc, numSect
  integer*4 montage
  real*4 dmin, dmax, dmean
  logical useMdoc, hdfFile
  integer*4 iiuFileType, iiuRetAdocIndex, AdocSetCurrent, AdocGetImageMetaInfo
  !
  character*320 inFile, outFile, metaFile
  character*10 metaOrHDF /'metadata'/
  !
  equivalence (nz, nxyz(3))

  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetString, AdocOpenImageMetadata
  integer*4 PipGetInOutFile, PipParseInput, PipGetBoolean, PipGetLogical
  integer numOptions
  parameter (numOptions = 5)
  character*(120 * numOptions) options(1)
  options(1) = &
      'input:InputFile:FN:Name of input image file@'// &
      'output:OutputFile:FN:Name of output piece list file@'// &
      'mdoc:MdocMetadataFile:B:Take coordinates from metadata file named'// &
      ' inputfile.mdoc (omit name)@'// &
      'other:OtherMetadataFile:FN:Other metadata file to take '// &
      'coordinates from (follow with name)@'// &
      'help:usage:B:Print help output'

  metaFile = ' '
  inFile = ' '
  useMdoc = .false.
  hdfFile = .false.
  nz = 0
  maxPiece = 0
  !
  call PipExitOnError(0, 'ERROR: EXTRACTPIECES - ')
  call PipAllowCommaDefaults(1)
  ierr = PipParseInput(options, numOptions, '@', numOptArg, numNonOptArg)
  if (PipGetBoolean('help', ierr) == 0) then
    call PipPrintHelp('extractpieces', 0, 1, 1)
    call exit(0)
  endif
  ifAddMdoc = PipGetString('OtherMetadataFile', metaFile)
  ierr = PipGetLogical('MdocMetadataFile', useMdoc)
  if (useMdoc .and. ifAddMdoc == 0) call exitError( &
      'You cannot enter both -mdoc and -other')
  if (PipGetInOutFile('InputFile', 1, 'Name of input image file', inFile) &
      .ne. 0 .and. ifAddMdoc .ne. 0) call exitError( &
      'You must enter either an input image file or a metadata file with -other')
  if (PipGetInOutFile('OutputFile', 2, 'Name of output piece list file', &
      outFile) .ne. 0) call exitError('No output file specified')
  call PipDone()
  !
  if (inFile .ne. ' ') then
    call imopen(1, inFile, 'RO')
    call irdhdr(1, nxyz, mxyz, mode, dmin, dmax, dmean)
    hdfFile = iiuFileType(1) == 5
    if (hdfFile .and. ifAddMdoc .ne. 0) then
      indAdoc = iiuRetAdocIndex(1, 0, 0)
      if (indAdoc <= 0) call exitError('Getting autodoc index for hdf file')
      if (AdocSetCurrent(indAdoc) .ne. 0) call exitError( &
          'Setting autodoc structure of HDF file as current autodoc')
      if (AdocGetImageMetaInfo(montage, numSect, iTypeAdoc) < 0) then
        print *,'This HDF file does not have metadata about each section'
        call exit(0)
      endif
      metaOrHDF = 'HDF'
    endif
  endif
  !
  ! try to open a metadata file regardless, unless got one from HDF file
  if (metaFile == ' ') metaFile = inFile
  if (.not.hdfFile .or. ifAddMdoc == 0) &
      indAdoc = AdocOpenImageMetadata(metaFile, ifAddMdoc, montage, numSect, iTypeAdoc)
  if (indAdoc > 0 .and. inFile == ' ') nz = numSect
  numPieces = 0

  if (.not.useMdoc .and. ifAddMdoc .ne. 0 .and. .not.hdfFile) then
    !
    ! Get data from the image header
    call iiuRetNumExtended(1, numExtraBytes)
    call iiuRetExtendedType(1, nbytes, iflags)
    maxExtra = numExtraBytes + 1024
    maxPiece = nz + 1024
    allocate(array(maxExtra / 4), ixPiece(maxPiece), iyPiece(maxPiece), &
        izPiece(maxPiece), stat = ierr)
    call memoryError(ierr, 'arrays for extra header or piece data')
    call iiuRetExtendedData(1, numExtraBytes, array)
    call get_extra_header_pieces(array, numExtraBytes, nbytes, iflags, nz, &
        ixPiece, iyPiece, izPiece, numPieces, maxPiece)
    if (numPieces == 0) then
      print *,'There are no piece coordinates in this image file'
      if (indAdoc .ne. -2)  &
          print *,'Looking for piece coordinates in the associated metadata file'
    endif
  endif
  
  if (numPieces == 0) then
    !
    ! Or get data from the autodoc file
    !
    ! Thanks to a bug in SerialEM, Montage flag may be missing.  So just plow
    ! ahead regardless
    if (indAdoc > 0 .and. numSect == nz) then
      if (maxPiece == 0) then
        maxPiece = nz + 1024
        allocate(ixPiece(maxPiece), iyPiece(maxPiece), izPiece(maxPiece), stat = ierr)
        call memoryError(ierr, 'arrays for piece data')
      endif
      call get_metadata_pieces(indAdoc, iTypeAdoc, nz, ixPiece, iyPiece, &
          izPiece, maxPiece, numPieces)
    endif
    !
    ! Give lots of different error messages
    if (indAdoc < 0 .or. numPieces < nz) then
      if (indAdoc == -2) then
        if (useMdoc .or. ifAddMdoc == 0) print *,'The metadata file does not exist'
      else if (indAdoc == -3) then
        print *,'The autodoc file is not a recognized type of image metadata file'
      else if (indAdoc == -1) then
        print *,'There was an error opening or reading the metadata file'
      else if (numPieces > 0 .or. numSect .ne. nz) &
          then
        write(*,'(a,a,a)')'The ', trim(metaOrHDF),' file does not have piece '// &
            'coordinates for every image in the file'
      else
        write(*,'(a,a,a)')'There are no piece coordinates in the ', trim(metaOrHDF), &
            ' file'
      endif
    endif
  endif
  !
  ! It used to output whatever is there, even if it is short, so restore that 6/24/16
  if (numPieces > 0) then
    call dopen(1, outFile, 'new', 'f')
    write(1, '(2i9,i7)') (ixPiece(i), iyPiece(i), izPiece(i), i = 1, numPieces)
    close(1)
    print *,numPieces, ' piece coordinates output to file'
  endif

  call iiuClose(1)
  !
  call exit(0)
end program extractpieces

