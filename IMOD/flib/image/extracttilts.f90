!*************EXTRACTTILTS.f**********************************************
!
! EXTRACTTILTS will extract tilt angles or other per-section
! information from the header of an image file, if they are present,
! and produce a file with a list of the values.
!
! See man page for more details
!
! David Mastronarde, 1/2/00
!
! $Id$
!
program extracttilts
  implicit none
  integer*4 maxExtra, maxTilts, maxPiece
  integer*4 nxyz(3), mxyz(3)
  real*4, allocatable :: tilt(:), val2(:), array(:)
  integer*4, allocatable :: ixPiece(:), iyPiece(:), izPiece(:)
  character*80, allocatable :: valString(:)
  !
  character*320 inFile, outFile, metafile, keyName
  character*10 zvalueName, globalName
  character*40 typeText(9) /'tilt angle', ' ', 'stage position', &
      'magnification', 'intensity value', 'exposure dose', 'pixel spacing', &
      'defocus', 'exposure time'/
  !
  integer*4 nz, ierr, ifTilt, ifMag, ifStage, numPieces, i, nbytes, iflags
  integer*4 maxz, iunitOut, lenText, mode, itype, ifC2, numTilts, numExtraBytes
  integer*4 numTiltOut, ifDose, ifWarn, ifAll, ifCamExp, ifDefocus, ifPixel, numFound
  integer*4 ifAddMdoc, indAdoc, iTypeAdoc, numSect, montage, ifMdoc, ifNoImage, ifKey
  integer*4 keyValType, ifAttrib
  logical*4 mdocPrimary, hdfFile
  integer*4 AdocOpenImageMetadata, iiuFileType, iiuRetAdocIndex, AdocSetCurrent
  integer*4 AdocGetImageMetaInfo, AdocWrite, AdocOrderWriteByValue
  real*4 dmin, dmax, dmean
  equivalence (nz, nxyz(3))
  !
  logical pipinput
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetBoolean, PipGetInOutFile, PipGetString
  !
  ! fallbacks from ../../manpages/autodoc2man -3 2  extracttilts
  !
  integer numOptions
  parameter (numOptions = 17)
  character*(40 * numOptions) options(1)
  options(1) = &
      'input:InputFile:FN:@output:OutputFile:FN:@mdoc:MdocMetadataFile:B:@'// &
      'other:OtherMetadataFile:FN:@tilts:TiltAngles:B:@stage:StagePositions:B:@'// &
      'mag:Magnifications:B:@intensities:Intensities:B:@exp:ExposureDose:B:@'// &
      'camera:CameraExposure:B:@pixel:PixelSpacing:B:@defocus:Defocus:B:@'// &
      'key:KeyName:CH:@warn:WarnIfTiltsSuspicious:B:@'// &
      'all:AllPieces:B:@help:usage:B:'
  !
  outFile = ' '
  ifMag = 0
  ifStage = 0
  ifC2 = 0
  ifTilt = 0
  ifDose = 0
  ifCamExp = 0
  ifPixel = 0
  ifDefocus = 0
  ifWarn = 0
  numPieces = 0
  ifAll = 0
  ifMdoc = 0
  ifAddMdoc = 1
  metafile = ' '
  ifKey = 0
  hdfFile = .false.
  ifAttrib = 0
  !
  ! Pip startup: set error, parse options, check help, set flag if used
  !
  call PipReadOrParseOptions(options, numOptions, 'extracttilts', &
      'ERROR: EXTRACTTILTS - ', .true., 1, 1, 1, numOptArg, numNonOptArg)
  pipinput = numOptArg + numNonOptArg > 0

  ifNoImage = PipGetInOutFile('InputFile', 1, 'Image input file', inFile)
  ierr = PipGetInOutFile('OutputFile', 2, &
      'Name of output file, or return to print out values', outFile)
  !
  if (pipinput) then
    ierr = PipGetBoolean('AllPieces', ifAll)
    ierr = PipGetBoolean('TiltAngles', ifTilt)
    ierr = PipGetBoolean('Magnifications', ifMag)
    ierr = PipGetBoolean('StagePositions', ifStage)
    ierr = PipGetBoolean('Intensities', ifC2)
    ierr = PipGetBoolean('ExposureDose', ifDose)
    ierr = PipGetBoolean('CameraExposure', ifCamExp)
    ierr = PipGetBoolean('Defocus', ifDefocus)
    ierr = PipGetBoolean('PixelSpacing', ifPixel)
    ierr = PipGetBoolean('MdocMetadataFile', ifMdoc)
    ierr = PipGetBoolean('WarnIfTiltsSuspicious', ifWarn)
    ierr = PipGetBoolean('AttributesInHDFfile', ifAttrib)
    ifAddMdoc = PipGetString('OtherMetadataFile', metafile)
    if (ifMdoc .ne. 0 .and. metafile .ne. ' ')  &
        call exitError('You cannot enter both -mdoc and -other')
    if (PipGetString('KeyName', keyName) == 0) then
      ifKey = 1
      itype = 2
      typeText(2) = keyName
    endif
    ierr = ifTilt + ifMag + ifStage + ifC2 + ifDose + ifDefocus + ifCamExp + ifPixel + &
        ifKey
    if (ifAttrib > 0 .and. ierr > 0) call exitError( &
        'You cannot specify anything to extract when writing attributes to a file')
    if (ierr == 0) ifTilt = 1
    if (ierr > 1) call exitError('You must enter only one option for data to extract')
    if (ifTilt .ne. 0) itype = 1
    if (ifStage .ne. 0) itype = 3
    if (ifMag .ne. 0) itype = 4
    if (ifC2 .ne. 0) itype = 5
    if (ifDose .ne. 0) itype = 6
    if (ifPixel .ne. 0) itype = 7
    if (ifDefocus .ne. 0) itype = 8
    if (ifCamExp .ne. 0) itype = 9
  else
    ifTilt = 1
    itype = 1
  endif

  if (ifAttrib > 0 .and. ifNoImage .ne. 0) call exitError( &
      'You must enter an input image file with -attrib')
  if (ifAttrib > 0 .and. (ifMdoc .ne. 0 .or. metafile .ne. ' ')) call exitError( &
      'You cannot enter -mdoc or -other with -attrib')
  if (ifAttrib > 0 .and. outFile == ' ') call exitError( &
      'You must enter an output file to extract attributes')

  ! Require an image unless -other is used with an mdoc file
  if (metafile == ' ' .and. ifNoImage .ne. 0) call exitError( &
      'You must enter either an input image file or a metadata file with -other')
  call PipDone()

  maxExtra = 40
  if (ifNoImage == 0) then
    call imopen(1, inFile, 'RO')
    call irdhdr(1, nxyz, mxyz, mode, dmin, dmax, dmean)
    !
    hdfFile = iiuFileType(1) == 5
    if (ifAttrib > 0 .and. .not.hdfFile) call exitError( &
        'The input image file is not an HDF file; cannot extract attributes')
    if (hdfFile .and. ifAddMdoc .ne. 0) then
      indAdoc = iiuRetAdocIndex(1, 0, 0)
      if (indAdoc <= 0) call exitError('Getting autodoc index for HDF file')
      if (AdocSetCurrent(indAdoc) .ne. 0) call exitError( &
          'Setting autodoc structure of HDF file as current autodoc')
      numSect = nz
      !
      ! Write out autodoc and exit if requested
      if (ifAttrib > 0) then
        call AdocGetStandardNames(globalName, zvalueName)
        if (AdocOrderWriteByValue(zvalueName) .ne. 0 .or. AdocWrite(outFile) .ne. 0) &
            call exitError('Writing autodoc structure to file')
        print *,'Attributes written to ', trim(outFile)
        call exit(0)
      endif
      if (AdocGetImageMetaInfo(montage, numSect, iTypeAdoc) < 0 .or. numSect < nz)  &
          call exitError('This HDF file does not have metadata about each section')
    else
      call iiuRetNumExtended(1, numExtraBytes)
      call iiuRetExtendedType(1, nbytes, iflags)
      maxExtra = numExtraBytes + 1024
    endif
  endif

  !
  ! Always try to open a metadata file, and check for errors if it was required
  if (metafile == ' ') metafile = inFile
  if (.not.hdfFile .or. ifAddMdoc == 0) &
      indAdoc = AdocOpenImageMetadata(metafile, ifAddMdoc, montage, numSect, iTypeAdoc)
  mdocPrimary = ifMdoc .ne. 0 .or. ifAddMdoc == 0 .or. hdfFile
  if (mdocPrimary) then
    if (indAdoc == -1) call exitError('Opening or reading the metadata file')
    if (indAdoc == -2) call exitError('The metadata file does not exist')
    if (indAdoc == -3) call exitError('The autodoc file is not '// &
        'a recognized form of image metadata file')
    if (ifNoImage > 0) nz = numSect
    if (numSect .ne. nz) call exitError('The image and metadata files'// &
        ' do not have the same number of sections')
  endif

  ! If the adoc is just a fallback, just disallow it if it doesn't match nz
  if (numSect .ne. nz) indAdoc = -4
  if ((ifKey > 0 .or. itype > 6) .and. indAdoc < 0) call exitError( &
      'This kind of information can be obtained only from a metadata file')

  maxPiece = nz + 1024
  maxTilts = maxPiece
  allocate(tilt(maxTilts), val2(maxTilts), array(maxExtra / 4), &
      ixPiece(maxPiece), iyPiece(maxPiece), izPiece(maxPiece), stat = ierr)
  call memoryError(ierr, 'arrays for extra header data')
  if (ifKey > 0) then
    allocate(valString(maxTilts), stat = ierr)
    call memoryError(ierr, 'array for string metadata')
  endif

  if (ifNoImage == 0) call iiuRetExtendedData(1, numExtraBytes, array)

  !
  ! Get piece coordinates unless doing "all".  First get from image file unless mdoc is
  ! specified; then get from the mdoc if it is primary or if image file gave nothing
  if (ifAll == 0) then
    if (.not. mdocPrimary) then
      call get_extra_header_pieces (array, numExtraBytes, nbytes, iflags, &
          nz, ixPiece, iyPiece, izPiece, numPieces, maxPiece)
    endif
    if ((montage .ne. 0 .or. hdfFile) .and. (mdocPrimary .or.  &
        (indAdoc >= 0 .and. numPieces == 0))) then
      !
      ! Get from any HDF file in case montage flag is missing; but if there are no
      ! pieces then move on without error
      call get_metadata_pieces(indAdoc, iTypeAdoc, nz, ixPiece, iyPiece, izPiece, &
          maxPiece, numPieces)
      if (numPieces > 0) montage = 1
      if (montage > 0 .and. numPieces < nz) then
        if (hdfFile .and. ifAddMdoc .ne. 0) call exitError('The HDF file is marked as'//&
            ' a montage but does not have piece coordinates for every section')
        if (mdocPrimary) call exitError('The metadata file does not '// &
            'have piece coordinates for every section')
        call exitError('There are no piece coordinates in the image file; '// & 
            'the metadata file indicates a montage but does not have piece '// &
            'coordinates for every section')
      endif
    endif
  endif
  if (numPieces == 0) then
    do i = 1, nz
      izPiece(i) = i - 1
    enddo
    maxz = nz
  else
    maxz = 0
    do i = 1, numPieces
      maxz = max(maxz, izPiece(i) + 1)
    enddo
  endif
  !
  ! set up a marker value for empty slots
  !
  do i = 1, maxz
    tilt(i) = -999.
  enddo
  !
  numTilts = 0
  if (.not. mdocPrimary .and. itype <= 6 .and. ifKey == 0) then
    call get_extra_header_items(array, numExtraBytes, nbytes, iflags, nz, itype, &
        tilt, val2, numTilts, maxTilts, izPiece)
    if (numTilts == 0 .and. indAdoc < 0) then
      write(*,'(/,a,a,a)') 'ERROR: EXTRACTTILTS - No ', trim(typeText(itype)), &
          ' information in this image file'
      call exit(1)
    endif
  endif
  if (mdocPrimary .or. itype > 6 .or. numTilts == 0 .or. ifKey > 0) then
    if (.not. mdocPrimary) print *,'Taking information from associated metadata file'
    if (ifKey == 0) then
      call get_metadata_items(indAdoc, iTypeAdoc, nz, itype, tilt, val2, &
          numTilts, numFound, maxTilts, izPiece)
    else
      call get_metadata_by_key(indAdoc, iTypeAdoc, nz, keyName, 0, tilt, val2, val2, &
          valString, numTilts, numFound, maxTilts, izPiece)
    endif
    if (numFound .ne. nz) then
      if (.not. mdocPrimary .and. itype <= 6) then
        write(*,'(/,a,a,a)') 'ERROR: EXTRACTTILTS - ', trim(typeText(itype)), &
            ' information is not present in image file and is missing for all or '// &
            'some sections in metadata file'
      else if (hdfFile .and. ifAddMdoc .ne. 0) then
        write(*,'(/,a,a,a)') 'ERROR: EXTRACTTILTS - ', trim(typeText(itype)), &
            ' is missing for all or some sections in HDF file'
      else
        write(*,'(/,a,a,a)') 'ERROR: EXTRACTTILTS - ', trim(typeText(itype)), &
            ' is missing for all or some sections in metadata file'
      endif
      call exit(1)
    endif
  endif
  !
  ! pack the tilts down
  !
  numTiltOut = 0
  do i = 1, numTilts
    if (tilt(i) .ne. -999.) then
      numTiltOut = numTiltOut + 1
      tilt(numTiltOut) = tilt(i)
      val2(numTiltOut) = val2(i)
      if (ifKey > 0) valString(numTiltOut) = valString(i)
    endif
  enddo
  !
  iunitOut = 6
  if (outFile .ne. ' ') then
    call dopen(1, outFile, 'new', 'f')
    iunitOut = 1
  endif

  if (ifTilt .ne. 0) write(iunitOut, '(f7.2)') (tilt(i), i = 1, numTiltOut)
  if (ifMag .ne. 0) write(iunitOut, '(i7)') (nint(tilt(i)), i = 1, numTiltOut)
  if (ifC2 .ne. 0) write(iunitOut, '(f8.5)') (tilt(i), i = 1, numTiltOut)
  if (ifStage .ne. 0)  write(iunitOut, '(2f9.2)') (tilt(i), val2(i), i = 1, numTiltOut)
  if (ifDose .ne. 0) write(iunitOut, '(f13.5)') (tilt(i), i = 1, numTiltOut)
  if (ifDefocus .ne. 0 .or. ifPixel .ne. 0 .or. ifCamExp .ne. 0) &
      write(iunitOut, '(f11.3)') (tilt(i), i = 1, numTiltOut)
  if (ifKey > 0) write(iunitOut, '(a)'), (valString(i), i = 1, numTiltOut)

  if (iunitOut == 1) then
    close(1)
    print *,numTiltOut, ' ', trim(typeText(itype)), 's output to file'
  endif

  if (ifTilt > 0 .and. ifWarn > 0) then
    numTilts = 0
    do i = 1, numTiltOut
      if (abs(tilt(i)) < 0.1) numTilts = numTilts + 1
      if (abs(tilt(i)) > 95.) ifMag = ifMag + 1
    enddo
    if (numTiltOut > 2 .and. numTilts > numTiltOut / 2) &
        write(*,103) numTilts, 'near zero'
103 format('WARNING: extracttilts - ',i4, ' of the extracted tilt angles are ',a)
    if (ifMag > 0) write(*,103) ifMag, 'greater than 95 degrees'
  endif

  call iiuClose(1)
  !
  call exit(0)
end program extracttilts
