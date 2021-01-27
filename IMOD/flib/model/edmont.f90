!*************EDMONT.F**********************************************
!
! A GENERAL MONTAGE EDITOR TO MOVE IMAGES INTO, OUT OF, OR BETWEEN
! MONTAGES It can float the images to a common range or mean of density,
! scaling all of the pieces in a section by the same amount. It can
! output only the pieces that intersect a defined subset of the image
! area.
!
! $Id$
program edmont
  !
  implicit none
  include 'smallmodel.inc90'
  !
  integer*4 nxyz(3), mxyz(3), nxyzst(3), nxyz2(3), mxyz2(3), nx, ny, nz
  real*4 title(20), cell2(6), cell(6), delta(3), xOrigin, yOrigin, zOrigin
  !
  character*320, allocatable :: inFile(:), outFile(:), pieceFileIn(:), pieceFileOut(:)
  !
  equivalence (nx, nxyz(1)), (ny, nxyz(2)), (nz, nxyz(3))
  character*320 modelFile
  !
  data nxyzst/0, 0, 0/
  character*20 floatText/' '/, truncText/' '/
  integer*4, allocatable :: inZlist(:), inlistTmp(:), nlist(:), listInd(:) &
      , numSecOut(:), listz(:) , ixPcList(:), iyPcList(:) , izPcList(:) &
      , ixPieceOut(:), iyPieceOut(:), izPieceOut(:), ixko(:), iyko(:), izko(:)
  real*4 optimalMax(16)
  real*4, allocatable :: dminSec(:), dmaxSec(:), avgsec(:), sdsec(:), array(:)
  real*4, allocatable :: tempLine(:)
  integer*1, allocatable :: extraIn(:), extraOut(:)
  data optimalMax/255., 32767., 4*255., 65535., 2*255., 511., 1023., 2047., 4095., &
      8191., 16383., 32767./
  logical*4 rescale, allHaveHeaderCoords, useMdoc, outDocChanged
  integer(kind = 8) i8, numPixels
  character dat*9, tim*8
  integer*2 temp
  character*80 titlech
  integer*4 maxExtra, maxPieces, limExtra, limPiece, idim, limList
  integer*4 ierr, numFilesIn, numImageInFiles, numPieceFiles, numSecLists
  integer*4 listTotal, ifile, i, mode, nbyteExtraIn, numInZlist, minXpiece , numXpieces, j
  integer*4 mxOverlap, minYpiece, numYpieces, myOverlap, newMode, nxOverlap
  integer*4 nyOverlap, minAllXpiece, minAllYpiece, maxXpiece, maxYpiece
  integer*4 numOutValues, numToGet, numOutFiles, numOutTotal, numOutEntries
  integer*4 numImageOutFiles, noutXpiece, noutYpiece, minXinc, maxXinc, minYinc
  integer*4 maxYinc, minFrame, maxFrame, knocks, ifKnock, iobj, ipt, ifMean
  integer*4 ifFloat, ifRenumber, ifShiftXy, ibinning, maxBinning, nxBin, numberOffset
  integer*4 nxbOverlap, ixStart, nyBin, nybOverlap, iyStart, lenTemp, ilis, indSec
  real*4 fracZero, zmin, zmax, dmin2, dmax2, avg, sd, tmpMin, tmpMax, dmaxOut
  real*8 sum8, sumsq8, tsum8, tsumsq8, dmean2
  real*4 optimalIn, optimalOut, bottomIn, bottomOut, sum, zminSec, zmaxSec, dminOut
  real*4 denOutMin, scaleFac, const, den, dmax, dmean, dmin, dminIn
  integer*4 numSecRead, ifAnyOut, ipc, isec, isecOut, ifileOut, ipieceOut, kti, numPcList
  integer*4 nbytePerSecOut, nbyteExtraOut, indXout, minSubXpiece, minSubYpiece
  integer*4 nbyteClear, nbytePerSecIn, iflagExtraIn, numPLfileOut, maxSecOut, ind
  integer*4 maxExtraOut, maxNumXpieces, maxNumYpieces, maxByteExtraIn, numByteCopy
  integer*4 indAdocIn, indAdocOut, indSectOut, indSectIn
  character*100000 listString
  character*10 zvalueName, globalName

  logical readSmallMod, notKnock
  integer*4 AdocSetThreeIntegers, iiuFileType, b3dOutputFileType
  integer*4 iiuRetAdocIndex, AdocLookupByNameValue, AdocTransferSection, AdocWrite
  !
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetInteger, PipGetBoolean, PipGetLogical
  integer*4 PipGetString, PipGetTwoIntegers
  integer*4 PipGetIntegerArray, PipGetNonOptionArg
  !
  ! fallbacks from ../../manpages/autodoc2man -3 2 edmont
  integer numOptions
  parameter (numOptions = 18)
  character*(40 * numOptions) options(1)
  options(1) = &
      'imin:ImageInputFile:FNM:@plin:PieceListInput:FNM:@'// &
      'imout:ImageOutputFile:FNM:@plout:PieceListOutput:FNM:@'// &
      'secs:SectionsToRead:LIM:@numout:NumberToOutput:IAM:@'// &
      'mode:ModeToOutput:I:@xminmax:XMinAndMax:IP:@'// &
      'Yminmax:YMinAndMax:IP:@xframes:XFrameMinAndMax:IP:@'// &
      'yframes:YFrameMinAndMax:IP:@float:FloatDensities:I:@'// &
      'bin:BinByFactor:I:@exclude:ExclusionModel:FN:@'// &
      'renumber:RenumberZFromZero:B:@shift:ShiftXYToZero:B:@'// &
      'param:ParameterFile:PF:@help:usage:B:'
  !
  !
  limExtra = 20000000
  limPiece = 2500000
  limList = 10000000
  maxBinning = 32
  fracZero = 0.
  ifMean = 0
  ifFloat = 0
  ifRenumber = 0
  ifShiftXy = 0
  ibinning = 1
  minXinc = 0
  maxXinc = 0
  minYinc = 0
  maxYinc = 0
  numberOffset = 0
  allHaveHeaderCoords = .true.
  useMdoc = .false.
  call AdocGetStandardNames(globalName, zvalueName)

  call PipReadOrParseOptions(options, numOptions, 'edmont', &
      'ERROR: EDMONT - ', .false., 2, 2, 1, numOptArg, &
      numNonOptArg)
  !
  ! Allocate oversized temporary arrays
  allocate(inlistTmp(limList), listz(limPiece), ixPcList(limPiece), &
      iyPcList(limPiece) , izPcList(limPiece), extraIn(limExtra), &
      stat = ierr)
  call memoryError(ierr, 'initial oversized arrays')

  !
  ! Get the input files
  call PipNumberOfEntries('ImageInputFile', numImageInFiles)
  numFilesIn = numImageInFiles + max(0, numNonOptArg - 1)
  if (numFilesIn == 0) call exitError('No input image file name entered')
  call PipNumberOfEntries('PieceListInput', numPieceFiles)
  if (numPieceFiles > 0 .and. numPieceFiles .ne. numFilesIn) call exitError &
      ('There must be one piece list entry for each image file if '// &
      'there are any piece list files')
  call PipNumberOfEntries('SectionsToRead', numSecLists)
  if (numSecLists > numFilesIn) call exitError( &
      'There are more section lists than input files')
  ierr = PipGetBoolean('NumberedFromOne', numberOffset)
  ierr = PipGetLogical('UseMdocFiles', useMdoc)

  allocate(inFile(numFilesIn), pieceFileIn(numFilesIn), nlist(numFilesIn), &
      listInd(numFilesIn), stat = ierr)
  call memoryError(ierr, 'arrays for files')

  listTotal = 0
  maxExtra = 1
  maxPieces = 0
  maxByteExtraIn = 0
  maxNumXpieces = 0
  maxNumYpieces = 0
  nxOverlap = -10000
  nyOverlap = -10000
  do ifile = 1, numFilesIn
    if (ifile <= numImageInFiles) then
      ierr = PipGetString('ImageInputFile', inFile(ifile))
    else
      ierr = PipGetNonOptionArg(ifile - numImageInFiles, inFile(ifile))
    endif
    pieceFileIn(ifile) = ' '
    if (numPieceFiles > 0) then
      ierr = PipGetString('PieceListInput', pieceFileIn(ifile))
      if (pieceFileIn(ifile) == 'none') pieceFileIn(ifile) = ' '
    endif
    !
    ! open the files to get properties and check them
    call openAndAnalyzeFiles(inFile(ifile), pieceFileIn(ifile), mxyz, dminIn, &
        .false., mode, extraIn, limExtra, nbyteExtraIn, nbytePerSecIn, ixPcList, &
        iyPcList, izPcList, numPcList, limPiece, listz, numInZlist, minXpiece, &
        numXpieces, mxOverlap, minYpiece, numYpieces, myOverlap, useMdoc)
    call iiuClose(1)
    maxExtra = max(maxExtra, nbyteExtraIn)
    maxPieces = max(maxPieces, numPcList)
    maxByteExtraIn = max(maxByteExtraIn, nbytePerSecIn)
    maxNumXpieces = max(maxNumXpieces, numXpieces)
    maxNumYpieces = max(maxNumYpieces, numYpieces)
    if (nbyteExtraIn == 0) allHaveHeaderCoords = .false.
    !
    ! The first file with multiple pieces defines the overlap on an axis
    if (nxOverlap == -10000 .and. numXpieces > 1) nxOverlap = mxOverlap
    if (nyOverlap == -10000 .and. numYpieces > 1) nyOverlap = myOverlap
    !
    ! The first file defines output mode and size; initialize for mins/maxs
    if (ifile == 1) then
      newMode = mode
      nxyz = mxyz
      !
      ! Keep track of overall minimum and maximum piece coordinate
      minAllXpiece = minXpiece
      minAllYpiece = minYpiece
      maxXpiece = minXpiece + (numXpieces - 1) * (nx - nxOverlap)
      maxYpiece = minYpiece + (numYpieces - 1) * (ny - nyOverlap)
    else
      !
      ! if overlap still not defined on an axis and there is now a
      ! disparity, try to set up an overlap that makes sense
      if (nxOverlap == -10000 .and. minXpiece .ne. minAllXpiece) then
        i = 1
        do while (abs(minXpiece - minAllXpiece) / i > nx)
          i = i + 1
        enddo
        nxOverlap = nx - abs(minXpiece - minAllXpiece) / i
      endif
      if (nyOverlap == -10000 .and. minYpiece .ne. minAllYpiece) then
        i = 1
        do while (abs(minYpiece - minAllYpiece) / i > ny)
          i = i + 1
        enddo
        nyOverlap = ny - abs(minYpiece - minAllYpiece) / i
      endif
      if (nx .ne. mxyz(1) .or. ny .ne. mxyz(2)) &
          call exitError('All image files must have the same size pieces')
      if ((numXpieces > 1 .and. nxOverlap .ne. mxOverlap) .or. &
          (numYpieces > 1 .and. nyOverlap .ne. myOverlap)) &
          call exitError('All montages must have the same overlaps')
      if (mod(abs(minXpiece - minAllXpiece), nx - nxOverlap) > 0 .or. &
          mod(abs(minYpiece - minAllYpiece), ny - nyOverlap) > 0) &
          call exitError('All montages must have piece coordinates on'// &
          ' the same regular grid')
      minAllXpiece = min(minAllXpiece, minXpiece)
      minAllYpiece = min(minAllYpiece, minYpiece)
      maxXpiece = max(maxXpiece, minXpiece + (numXpieces - 1) * (nx - nxOverlap))
      maxYpiece = max(maxYpiece, minYpiece + (numYpieces - 1) * (ny - nyOverlap))
    endif
    !
    ! Initialize the section list then get it if there is one
    if (listTotal + numInZlist > limList) &
        call exitError('Too many sections for initial arrays')
    inlistTmp(listTotal + 1:listTotal + numInZlist) = listz(1:numInZlist)
    nlist(ifile) = numInZlist
    if (ifile <= numSecLists) then
      ierr = PipGetString('SectionsToRead', listString)
      call parselist2(listString, inlistTmp(listTotal + 1), nlist(ifile), &
          limList - listTotal)
      !
      ! Check for legality
      do i = 1, nlist(ifile)
        inlistTmp(listTotal + i) = inlistTmp(listTotal + i) - numberOffset
        ierr = 1
        do j = 1, numInZlist
          if (listz(j) == inlistTmp(listTotal + i)) then
            ierr = 0
            exit
          endif
        enddo
        if (ierr == 1) then
          write(*,'(/,a,i5,a,i6,a)') 'ERROR: EDMONT - SECTION LIST #', &
              ifile, ' CONTAINS ', &
              inlistTmp(listTotal + i), ', WHICH IS NOT A SECTION IN THAT FILE'
          call exit(1)
        endif
      enddo
    endif
    listInd(ifile) = listTotal + 1
    listTotal = listTotal + nlist(ifile)
  enddo
  !
  ! Get output files and numbers to output
  call PipNumberOfEntries('ImageOutputFile', numImageOutFiles)
  numOutFiles = numImageOutFiles + min(1, numNonOptArg)
  if (numOutFiles == 0) call exitError('No output image file name entered')
  allocate(outFile(numOutFiles), pieceFileOut(numOutFiles), numSecOut(numOutFiles), &
      stat = ierr)
  call memoryError(ierr, 'arrays for output files')
  call PipNumberOfEntries('PieceListOutput', numPLfileOut)
  if ((numPLfileOut == 0 .and. numPieceFiles > 0) .or. &
      (numPLfileOut .ne. 0 .and. numPLfileOut .ne. numOutFiles)) &
      call exitError('There must be an output piece list file for each '// &
      'output image file if there are input piece list files')
  !
  ! Take care of output numbers first so they can be totaled
  if (numOutFiles == 1) then
    numSecOut(1) = listTotal
  elseif (numOutFiles == listTotal) then
    do i = 1, numOutFiles
      numSecOut(i) = 1
    enddo
  else
    call PipNumberOfEntries('NumberToOutput', numOutEntries)
    if (numOutEntries == 0) &
        call exitError('You must specify number of sections '// &
        'to write to each output file')

    numOutValues = 0
    do i = 1, numOutEntries
      numToGet = 0
      ierr = PipGetIntegerArray('NumberToOutput', &
          numSecOut(numOutValues + 1), numToGet, numOutFiles - numOutValues)
      numOutValues = numOutValues + numToGet
    enddo
    if (numOutValues .ne. numOutFiles) call exitError( &
        'The number of values for sections to output does' &
        //' not equal the number of output files')
  endif

  numOutTotal = 0
  maxSecOut = 0
  do ifile = 1, numOutFiles
    if (ifile <= numImageOutFiles) then
      ierr = PipGetString('ImageOutputFile', outFile(ifile))
    else
      ierr = PipGetNonOptionArg(numNonOptArg, outFile(ifile))
    endif
    pieceFileOut(ifile) = ' '
    if (numPLfileOut > 0) &
        ierr = PipGetString('PieceListOutput', pieceFileOut(ifile))
    numOutTotal = numOutTotal + numSecOut(ifile)
    maxSecOut = max(maxSecOut, numSecOut(ifile))
  enddo
  if (numOutTotal .ne. listTotal) call exitError( &
      'Number of input and output sections does not match')
  noutXpiece = 1 + (maxXpiece - minAllXpiece) / (nx - nxOverlap)
  noutYpiece = 1 + (maxYpiece - minAllYpiece) / (ny - nyOverlap)
  !
  ! get limits on output size
  !
  ierr = PipGetTwoIntegers('XMinAndMax', minXinc, maxXinc)
  if (PipGetTwoIntegers('XFrameMinAndMax', minFrame, maxFrame) == 0) then
    if (ierr == 0) call exitError('You cannot enter both -xminmax and -xframes')
    if (minFrame < 1 .or. maxFrame > noutXpiece .or. &
        minFrame > maxFrame) call exitError &
        ('Minimum and maximum frames in X out of range or out of order')
    minXinc = minAllXpiece + (minFrame - 2) * (nx - nxOverlap) + nx
    maxXinc = minAllXpiece + maxFrame * (nx - nxOverlap) - 1
    maxNumXpieces = maxFrame + 1 - minFrame
  endif
  ierr = PipGetTwoIntegers('YMinAndMax', minYinc, maxYinc)
  if (PipGetTwoIntegers('YFrameMinAndMax', minFrame, maxFrame) == 0) then
    if (ierr == 0) call exitError('You cannot enter both -yminmax and -yframes')
    if (minFrame < 1 .or. maxFrame > noutYpiece .or. &
        minFrame > maxFrame) call exitError &
        ('Minimum and maximum frames in Y out of range or out of order')
    minYinc = minAllYpiece + (minFrame - 2) * (ny - nyOverlap) + ny
    maxYinc = minAllYpiece + maxFrame * (ny - nyOverlap) - 1
    maxNumYpieces = maxFrame + 1 - minFrame
  endif
  !
  ! Get model to knock out pieces
  knocks = 0
  ifKnock = 1 - PipGetString('ExclusionModel', modelFile)
  if (ifKnock .ne. 0) then
    if (.not.readSmallMod(modelFile)) call exitError( &
        'Reading piece exclusion model')
    call scale_model(0)
    do iobj = 1, max_mod_obj
      knocks = knocks + npt_in_obj(iobj)
    enddo
  endif
  allocate(ixko(knocks + 1), iyko(knocks + 1), izko(knocks + 1), stat = ierr)
  call memoryError(ierr, 'arrays for excluding pieces')
  if (ifKnock .ne. 0) then
    knocks = 0
    do iobj = 1, max_mod_obj
      do ipt = 1, npt_in_obj(iobj)
        knocks = knocks + 1
        ixko(knocks) = nint(p_coord(1, ibase_obj(iobj) + ipt))
        iyko(knocks) = nint(p_coord(2, ibase_obj(iobj) + ipt))
        izko(knocks) = nint(p_coord(3, ibase_obj(iobj) + ipt))
      enddo
    enddo
  endif
  !
  ierr = PipGetInteger('ModeToOutput', newMode)
  !
  ierr = PipGetInteger('FloatDensities', ifFloat)
  if (ifFloat > 1) ifMean = 1
  !
  ! renumber sections starting at zero: and if not, check for duplicate Z
  ierr = PipGetBoolean('RenumberZFromZero', ifRenumber)
  ierr = PipGetBoolean('ShiftXYToZero', ifShiftXy)
  ierr = PipGetInteger('BinByFactor', ibinning)
  if (ibinning < 1 .or. ibinning > maxBinning) &
      call exitError('Binning is outside of allowed range')
  minSubXpiece = minAllXpiece / ibinning
  minSubYpiece = minAllYpiece / ibinning
  call PipDone()
  !
  ! Manage memory
  allocate(inZlist(listTotal), stat = ierr)
  call memoryError(ierr, 'array for sections to do')
  inZlist(1:listTotal) = inlistTmp(1:listTotal)
  deallocate(inlistTmp, listz, ixPcList, iyPcList, izPcList, extraIn)
  maxExtraOut = maxSecOut * maxNumXpieces * maxNumYpieces * maxByteExtraIn + 1
  i = maxSecOut * maxNumXpieces * maxNumYpieces
  allocate(listz(maxPieces), ixPcList(maxPieces), &
      iyPcList(maxPieces) , izPcList(maxPieces), extraIn(maxExtra), &
      extraOut(maxExtraOut), dminSec(listTotal), dmaxSec(listTotal), avgsec(listTotal), &
      sdsec(listTotal), ixPieceOut(i), iyPieceOut(i), izPieceOut(i), stat = ierr)
  call memoryError(ierr, 'arrays for pieces and extra header')
  nxBin = nx / ibinning
  nxbOverlap = nint((nxOverlap - mod(nx, ibinning)) / float(ibinning))
  ixStart = mod(nx, ibinning) / 2
  nyBin = ny / ibinning
  nybOverlap = nint((nyOverlap - mod(ny, ibinning)) / float(ibinning))
  iyStart = mod(ny, ibinning) / 2
  lenTemp = 2 * ibinning * nx
  idim = nxBin * nyBin
  allocate(tempLine(lenTemp), array(idim), stat = ierr)
  call memoryError(ierr, 'arrays for image data')
  if (ifRenumber == 0 .or. ifFloat .ne. 0 .or. knocks > 0 .or. &
      minXinc .ne. 0 .or. &
      maxXinc .ne. 0 .or. minYinc .ne. 0 .or. maxYinc .ne. 0) then
    !
    ! need to go through all the files and check for unique output pieces
    ! if not renumbering in Z; get the actual range of pieces and check
    ! for existence of pieces on each section if any pieces are being
    ! removed, and get means if floating
    isecOut = 1
    ifileOut = 1
    ipieceOut = 0
    if (ifFloat .ne. 0) then
      floatText = ', floated to range'
      if (ifFloat < 0) floatText = ', scaled to range'
      if (fracZero .ne. 0.) &
          write(truncText, '(a,f6.3)') ', truncated by', fracZero
      if (ifMean .ne. 0) floatText = ', floated to means'
    endif
    zmin = 1.e30
    zmax = -1.e30
    dminOut = 1.e30
    dmaxOut = -1.e30
    minSubXpiece = maxXpiece + nx
    minSubYpiece = maxYpiece + ny
    do ifile = 1, numFilesIn
      call openAndAnalyzeFiles(inFile(ifile), pieceFileIn(ifile), mxyz, dminIn, &
          .false., mode, extraIn, maxExtra, nbyteExtraIn, nbytePerSecIn, ixPcList, &
          iyPcList, izPcList, numPcList, maxPieces, listz, numInZlist, &
          minXpiece , numXpieces, mxOverlap, minYpiece, numYpieces, myOverlap, useMdoc)

      do ilis = 1, nlist(ifile)
        indSec = ilis + listInd(ifile) - 1
        numSecRead = inZlist(indSec)
        dmin2 = 1.e30
        dmax2 = -1.e30
        sum8 = 0.
        sumsq8 = 0.
        ifAnyOut = 0
        do ipc = 1, numPcList
          if (izPcList(ipc) == numSecRead .and. &
              ((minXinc == 0 .and. maxXinc == 0) .or. &
              (max(ixPcList(ipc), minXinc) <= &
              min(ixPcList(ipc) + nx - 1, maxXinc))) .and. &
              ((minYinc == 0 .and. maxYinc == 0) .or. &
              (max(iyPcList(ipc), minYinc) <= &
              min(iyPcList(ipc) + ny - 1, maxYinc))) .and. &
              notKnock(ixPcList(ipc), iyPcList(ipc), numSecRead, nx, ny, nxOverlap, &
              nyOverlap, ixko, iyko, izko, knocks)) then
            minSubXpiece = min(minSubXpiece, ixPcList(ipc))
            minSubYpiece = min(minSubYpiece, iyPcList(ipc))
            if (ifFloat .ne. 0) then
              !
              ! if floating, need to read all sections to get stats
              ! find the minimum of the ratio (dmean-dmin) /(dmax-dmin)
              call irdbinned(1, ipc - 1, array, nxBin, nyBin, ixStart, &
                  iyStart, ibinning, nxBin, nyBin, tempLine, lenTemp, ierr)
              if (ierr .ne. 0) call exitError('Reading image file')
              if (ifMean == 0) then
                call iclden(array, nxBin, nyBin, 1, nxBin, 1, nyBin, &
                    tmpMin, tmpMax, sum)
              else
                call iclavgsd(array, nxBin, nyBin, 1, nxBin, 1, nyBin, &
                    tmpMin, tmpMax, tsum8, tsumsq8, avg, sd)

                sum8 = sum8 + tsum8
                sumsq8 = sumsq8 + tsumsq8
              endif
              dmin2 = min(dmin2, tmpMin)
              dmax2 = max(dmax2, tmpMax)
            endif
            ifAnyOut = ifAnyOut + 1
            if (ifRenumber == 0) then
              do i = 1, ipieceOut
                if (ixPieceOut(i) == ixPcList(ipc) .and. iyPieceOut(i) == &
                    iyPcList(ipc) .and. izPieceOut(i) == izPcList(ipc)) then
                  write(*,'(/,a,i5,a,3i8)') 'ERROR: EDMONT - YOU MUST '// &
                      'RENUMBER Z; OUTPUT FILE #', ifile, &
                      ' WOULD CONTAIN TWO SECTIONS WITH X,Y,Z OF', &
                      ixPieceOut(i), iyPieceOut(i), izPieceOut(i)
                  call exit(1)
                endif
              enddo
            endif
            ipieceOut = ipieceOut + 1
            ixPieceOut(ipieceOut) = ixPcList(ipc)
            iyPieceOut(ipieceOut) = iyPcList(ipc)
            izPieceOut(ipieceOut) = izPcList(ipc)
          endif
        enddo
        !
        ! Insist that some pieces  be found on this section!
        if (ifAnyOut == 0) then
          write(*,'(/,a,i6,a,i5)') 'ERROR: EDMONT - NO PIECES ARE '// &
              'INCLUDED FROM SECTION', numSecRead, &
              ' WHICH IS IN THE LIST OF SECTIONS TO USE FROM FILE #', ifile
          call exit(1)
        endif
        if (ifFloat .ne. 0) then
          dminSec(indSec) = dmin2
          dmaxSec(indSec) = dmax2
          dminOut = min(dminOut, dmin2)
          dmaxOut = max(dmaxOut, dmax2)
          if (ifMean .ne. 0) then
            call sums_to_avgsd8(sum8, sumsq8, nxBin, nyBin * ifAnyOut, avg, sd)
            avgsec(indSec) = avg
            sdsec(indSec) = sd
            zmin = min(zmin, (dmin2 - avg) / sd)
            zmax = max(zmax, (dmax2 - avg) / sd)
          endif
        endif
        isecOut = isecOut + 1
        if (isecOut > numSecOut(ifileOut)) then
          isecOut = 1
          ifileOut = ifileOut + 1
          ipieceOut = 0
        endif
      enddo
      call iiuClose(1)
    enddo
    !
    ! Convert minimum coordinates to binned coordinates
    i = (minSubXpiece - minAllXpiece) / (nx - nxOverlap)
    minSubXpiece = i * (nxBin - nxbOverlap) + minAllXpiece / ibinning
    i = (minSubYpiece - minAllYpiece) / (ny - nyOverlap)
    minSubYpiece = i * (nyBin - nybOverlap) + minAllYpiece / ibinning
  endif
  !
  ! start looping over input images
  !
  call time(tim)
  call b3ddate(dat)
  isec = 1
  isecOut = 1
  ifileOut = 1
  ipieceOut = 0
  do ifile = 1, numFilesIn
    call openAndAnalyzeFiles(inFile(ifile), pieceFileIn(ifile), mxyz, dminIn, &
        .true., mode, extraIn, maxExtra, nbyteExtraIn, nbytePerSecIn, ixPcList, &
        iyPcList, izPcList, numPcList, maxPieces, listz, numInZlist, minXpiece , &
        numXpieces, mxOverlap, minYpiece, numYpieces, myOverlap, useMdoc)
    call iiuRetSize(1, nxyz, mxyz, nxyzst)
    call iiuRetCell(1, cell)
    call iiuRetDelta(1, delta)
    call iiuRetOrigin(1, xOrigin, yOrigin, zOrigin)
    !
    ! get extra header information if any
    !
    if (nbyteExtraIn > 0) then
      call iiuRetExtendedData(1, nbyteExtraIn, extraIn)
      call iiuRetExtendedType(1, nbytePerSecIn, iflagExtraIn)
    endif
    indAdocIn = -1;
    if (iiuFileType(1) == 5 .or. useMdoc) then
      indAdocIn = iiuRetAdocIndex(1, 0, 1)
      if (indAdocIn < 0) then
        write(*,'(/,a,a)') 'ERROR: EDMONT - COULD NOT OPEN AUTODOC INFORMATION FOR '// &
            'INPUT FILE ', trim(inFile(ifile))
        call exit(1)
      endif
    endif
    !
    ! get each section in input file
    do ilis = 1, nlist(ifile)
      indSec = ilis + listInd(ifile) - 1
      numSecRead = inZlist(indSec)
      ifAnyOut = 0
      do ipc = 1, numPcList
        if (izPcList(ipc) == numSecRead .and. &
            ((minXinc == 0 .and. maxXinc == 0) .or. &
            (max(ixPcList(ipc), minXinc) <= &
            min(ixPcList(ipc) + nx - 1, maxXinc))) .and. &
            ((minYinc == 0 .and. maxYinc == 0) .or. &
            (max(iyPcList(ipc), minYinc) <= &
            min(iyPcList(ipc) + ny - 1, maxYinc))) .and. &
            notKnock(ixPcList(ipc), iyPcList(ipc), numSecRead, nx, ny, nxOverlap, &
            nyOverlap, ixko, iyko, izko, knocks)) then
          call irdbinned(1, ipc - 1, array, nxBin, nyBin, ixStart, &
              iyStart, ibinning, nxBin, nyBin, tempLine, lenTemp, ierr)
          if (ierr .ne. 0) call exitError('Reading image file')
          numPixels = int(nxBin, kind = 8) * nyBin
          !
          ! calculate new min and max after rescaling under various
          ! possibilities
          !
          optimalIn = optimalMax(mode + 1)
          optimalOut = optimalMax(newMode + 1)
          !
          ! set bottom of input range to 0 unless mode 1 or 2; set bottom
          ! of output range to 0 unless not changing modes
          !
          bottomIn = 0.
          if (dminIn < 0. .and. (mode == 1 .or. mode == 2)) bottomIn = -optimalIn
          bottomOut = 0.
          if (mode == newMode) bottomOut = bottomIn
          rescale = .false.
          if (ifFloat <= 0) then
            !
            ! get min and max
            call iclden(array, nxBin, nyBin, 1, nxBin, 1, nyBin, &
                tmpMin, tmpMax, sum)
            dmin2 = tmpMin
            dmax2 = tmpMax
          endif
          !
          ! for no float: if mode = 2, no rescale
          !
          if (ifFloat <= 0 .and. newMode .ne. 2 .and. &
              (mode .ne. 2 .or. ifFloat < 0)) then
            !
            ! if within proper input range, rescale from input range to
            ! output range only if mode is changing
            !
            rescale = mode .ne. newMode .and. (mode .ne. 2 .or. ifFloat < 0)
            if (ifFloat < 0 .and. rescale) then
              bottomIn = dminOut
              optimalIn = dmaxOut
            endif
            if (tmpMin >= bottomIn .and. tmpMax <= optimalIn) then
              dmin2 = (tmpMin - bottomIn) * (optimalOut - bottomOut) / &
                  (optimalIn - bottomIn) + bottomOut
              dmax2 = (tmpMax - bottomIn) * (optimalOut - bottomOut) / &
                  (optimalIn - bottomIn) + bottomOut
            elseif (rescale) then
              ! :if outside proper range, tell user to start over
              call exitError('Input data outside expected range. '// &
                  'Start over, specifying float to range')
            endif
          elseif (ifFloat > 0) then
            ! If floating: scale to a dmin2 that will knock out fraczero
            ! of the range after truncation to zero
            dmin2 = -optimalOut * fracZero / (1. -fracZero)
            rescale = .true.
            tmpMin = dminSec(indSec)
            tmpMax = dmaxSec(indSec)
            if (ifMean == 0) then
              ! :float to range, new dmax2 is the max of the range
              dmax2 = optimalOut
            else
              ! :float to mean, it's very hairy, first need mean again
              zminSec = (tmpMin - avgsec(indSec)) / sdsec(indSec)
              zmaxSec = (tmpMax - avgsec(indSec)) / sdsec(indSec)
              dmin2 = (zminSec - zmin) * optimalOut / (zmax - zmin)
              dmax2 = (zmaxSec - zmin) * optimalOut / (zmax - zmin)
              dmin2 = max(0., dmin2)
              ! but in case of problems, just limit dmax2 to the range
              dmax2 = min(dmax2, optimalOut)
            endif
          endif
          dmean2 = 0.
          ! set up minimum value to output based on mode
          if (newMode == 1) then
            denOutMin = -32767
          elseif (newMode == 2) then
            denOutMin = -1.e30
            optimalOut = 1.e30
          else
            denOutMin = 0.
          endif
          !
          if (rescale) then
            ! if scaling, set up equation, scale and compute new mean
            scaleFac = (dmax2 - dmin2) / (tmpMax - tmpMin)
            const = dmin2 - scaleFac * tmpMin
            dmin2 = 1.e20
            dmax2 = -1.e20
            do i8 = 1, numPixels
              den = scaleFac * array(i8) + const
              if (den < denOutMin) then
                ! ntrunclo=ntrunclo+1
                den = denOutMin
              elseif (den > optimalOut) then
                ! ntrunchi=ntrunchi+1
                den = optimalOut
              endif
              array(i8) = den
              dmean2 = dmean2 + den
              dmin2 = min(dmin2, den)
              dmax2 = max(dmax2, den)
            enddo
          else
            ! if not scaling, just need new mean
            do i8 = 1, numPixels
              dmean2 = dmean2 + array(i8)
            enddo
          endif
          dmean2 = dmean2 / numPixels
          print *,'frame', isec - 1, ': min&max before and after, mean:'
          write(*,'(5f10.2)') tmpMin, tmpMax, dmin2, dmax2, dmean2
          ! see if need to open an output file
          if (ipieceOut == 0) then
            !
            ! Create output file, transfer header from currently open
            ! file, fix it enough to get going
            call imopen(2, outFile(ifileOut), 'NEW')
            call iiuTransHeader(2, 1)
            call iiuAltMode(2, newMode)
            nxyz2(1) = nxBin
            nxyz2(2) = nyBin
            nxyz2(3) = 1
            call iiuAltSize(2, nxyz2, nxyzst)
            !
            ! adjust extra header information if current file has it
            !
            nbyteExtraOut = 0
            if (nbyteExtraIn > 0 .and. allHaveHeaderCoords) then
              nbytePerSecOut = nbytePerSecIn
              nbyteExtraOut = numSecOut(ifileOut) * maxByteExtraIn * maxNumXpieces * &
                  maxNumYpieces
              call iiuAltNumExtended(2, nbyteExtraOut)
              call imposn(2, 0, 0)
              indXout = 0
            endif
            !
            ! Open autodoc for output file if appropriate
            ! Transfer global data.  
            outDocChanged = .false.
            indAdocOut = -1
            if ((b3dOutputFileType() == 5 .or. useMdoc) .and. indAdocIn >= 0) then
              indAdocOut = iiuRetAdocIndex(2, 0, -1)
              if (indAdocOut <= 0) call exitError('Cannot get autodoc index for output')
              call setCurrentAdocOrExit(indAdocIn, 'input')
              if (AdocTransferSection(globalName, 1, indAdocOut, globalName, 0) .ne. 0) &
                  call exitError('Transferring global data between autodocs') 
            endif
            write(titlech, 301) floatText, dat, tim
301         format('EDMONT: Images transferred',a18,t57,a9,2x,a8)
            read(titlech, '(20a4)') (title(kti), kti = 1, 20)
            dmax = -100000.
            dmin = 100000.
            dmean = 0.
          endif
          !
          dmin = min(dmin, dmin2)
          dmax = max(dmax, dmax2)
          dmean = dmean + dmean2
          !
          call iiuWriteSection(2, array)
          isec = isec + 1
          ipieceOut = ipieceOut + 1
          i = (ixPcList(ipc) - minAllXpiece) / (nx - nxOverlap)
          ixPieceOut(ipieceOut) = i * (nxBin - nxbOverlap) + minAllXpiece / ibinning
          i = (iyPcList(ipc) - minAllYpiece) / (ny - nyOverlap)
          iyPieceOut(ipieceOut) = i * (nyBin - nybOverlap) + minAllYpiece / ibinning
          izPieceOut(ipieceOut) = izPcList(ipc)
          if (ifRenumber .ne. 0) izPieceOut(ipieceOut) = isecOut - 1
          if (ifShiftXy .ne. 0) then
            ixPieceOut(ipieceOut) = ixPieceOut(ipieceOut) - minSubXpiece
            iyPieceOut(ipieceOut) = iyPieceOut(ipieceOut) - minSubYpiece
          endif
          ifAnyOut = ifAnyOut + 1
          !
          ! transfer extra header bytes if present
          !
          if (nbyteExtraOut .ne. 0 .and. indXout < nbyteExtraOut) then
            numByteCopy = min(nbytePerSecOut, nbytePerSecIn, nbyteExtraIn)
            nbyteClear = nbytePerSecOut - numByteCopy
            do i = 1, numByteCopy
              indXout = indXout + 1
              extraOut(indXout) = extraIn((ipc - 1) * nbytePerSecIn + i)
            enddo
            do i = 1, nbyteClear
              indXout = indXout + 1
              extraOut(indXout) = 0
            enddo
            !
            ! Need to replace values if there is renumbering, shifting,
            ! binning, or the piece coordinates actually came from pl file
            if ((ifRenumber .ne. 0 .or. ifShiftXy .ne. 0 .or. ibinning > 1 &
                .or. pieceFileIn(ifile) .ne. ' ') &
                .and. numByteCopy >= 6) then
              temp = izPieceOut(ipieceOut)
              ind = indXout + 5 - nbytePerSecOut
              if (mod(iflagExtraIn, 2) .ne. 0) ind = ind + 2
              if (ifRenumber .ne. 0 .or. pieceFileIn(ifile) .ne. ' ') &
                  call move(extraOut(ind), temp, 2)
              if (ifShiftXy .ne. 0 .or. ibinning > 1 .or. &
                  pieceFileIn(ifile) .ne. ' ') then
                temp = ixPieceOut(ipieceOut)
                call move(extraOut(ind-4), temp, 2)
                temp = iyPieceOut(ipieceOut)
                call move(extraOut(ind-2), temp, 2)
              endif
            endif
          endif
          !
          ! Transfer an adoc section and set the piece coordinates
          if (indAdocIn > 0 .and. indAdocOut > 0) then
            call setCurrentAdocOrExit(indAdocIn, 'input')
            indSectIn = AdocLookupByNameValue(zvalueName, ipc - 1)
            if (indSectIn > 0) then
              call int_iwrite(listString, ipieceOut - 1, ierr)
              if (AdocTransferSection(zvalueName, indSectIn, indAdocOut, listString, 1) &
                  .ne. 0) call exitError('Transferring section data between autodocs')
              outDocChanged = .true.
              call setCurrentAdocOrExit(indAdocOut, 'output')
              indSectOut = AdocLookupByNameValue(zvalueName, ipieceOut - 1)
              if (indSectOut < 0 .or. AdocSetThreeIntegers(zvalueName, indSectOut,  &
                  'PieceCoordinates', ixPieceOut(ipieceOut), iyPieceOut(ipieceOut),  &
                  izPieceOut(ipieceOut)) .ne. 0) call exitError( &
                  'Setting piece coordinates in output autodoc')
            endif
          endif
        endif
      enddo
      if (ifAnyOut > 0) isecOut = isecOut + 1
      ! see if need to close stack file
      if (isecOut > numSecOut(ifileOut) .or. &
          (ifile == numFilesIn .and. ilis == nlist(ifile))) then
        ! set new size, keep old nxyzst
        if (pieceFileOut(ifileOut) .ne. ' ') then
          call dopen(3, pieceFileOut(ifileOut), 'new', 'f')
          write(3, '(2i9,i7)')(ixPieceOut(i), iyPieceOut(i), izPieceOut(i), i = 1, &
              ipieceOut)
          close(3)
        endif
        nxyz2(3) = ipieceOut
        call iiuAltSize(2, nxyz2, nxyzst)
        ! if mxyz=nxyz, keep this relationship
        if (mxyz(1) == nx .and. mxyz(2) == ny .and. mxyz(3) == nz) &
            then
          mxyz2(1) = nxBin
          mxyz2(2) = nyBin
          mxyz2(3) = ipieceOut
          call iiuAltSample(2, mxyz2)
        endif
        ! keep delta the same by scaling cell size from change in mxyz
        do i = 1, 3
          cell2(i) = mxyz2(i) * (cell(i) / mxyz(i))
          if (i < 3) cell2(i) = cell2(i) * ibinning
          cell2(i + 3) = 90.
        enddo
        call iiuAltCell(2, cell2)
        !
        ! adjust origin if shifting piece coords to 0
        if (ifShiftXy .ne. 0) &
            call iiuAltOrigin(2, xOrigin - ibinning * minSubXpiece * delta(1), &
            yOrigin - ibinning * minSubYpiece * delta(2), zOrigin)
        if (nbyteExtraOut > 0) call iiuAltExtendedData(2, nbyteExtraOut, extraOut)
        dmean = dmean / ipieceOut

        if (outDocChanged .and. iiuFileType(2) .ne. 5) then
          call setCurrentAdocOrExit(indAdocOut, 'output')
          if (AdocWrite(trim(outFile(ifileOut))//'.mdoc') .ne.0) call exitError( &
              'Writing mdoc file for output file')
          call AdocClear(indAdocOut)
        endif

        !
        call iwrhdr(2, title, 1, dmin, dmax, dmean)
        call iiuClose(2)
        ipieceOut = 0
        isecOut = 1
        ifileOut = ifileOut + 1
      endif
    enddo
    call iiuClose(1)
  enddo
  !
  call exit(0)
end program edmont


! Tests for whether a piece should be included given its coordinates
! and the set of points to knock out
!
logical function notKnock(ixPiece, iyPiece, izPiece, nx, ny, nxOverlap, nyOverlap, ixko, &
    iyko, izko, knocks)
  implicit none
  integer*4 ixko(*), iyko(*), izko(*), ixPiece, iyPiece, izPiece, nx, ny, knocks, i
  integer*4 nxOverlap, nyOverlap, nxBorder, nyBorder
  notKnock = .false.
  nxBorder = min(nxOverlap, nx / 4)
  nyBorder = min(nyOverlap, ny / 4)
  do i = 1, knocks
    if (izko(i) == izPiece .and. ixPiece + nxBorder <= ixko(i) .and. &
        ixPiece + nx - nxBorder > ixko(i) .and. iyPiece + nyBorder <= iyko(i) .and. &
        iyPiece + ny - nyBorder > iyko(i)) return
  enddo
  notKnock = .true.
  return
end function notKnock


! Given a piece list file, read the piece list; given no piece file,
! attempt to read piece coordinates from image header
!
subroutine read_pl_or_header(pieceFileIn, inFile, extraIn, maxExtra, nbyteExtraIn, &
    nbytePerSecIn, ixPcList, iyPcList, izPcList, numPcList, limPcList, useMdoc)
  implicit none
  character*(*) pieceFileIn, inFile
  integer*1 extraIn(*)
  logical*4 useMdoc
  integer*4 ixPcList(*), iyPcList(*), izPcList(*), maxExtra, nbyteExtraIn, numPcList
  integer*4 nxyz(3), mxyz(3), limPcList, mode, nbytePerSecIn, iflagExtraIn, indAdoc
  integer*4 montage, numSect, iTypeAdoc
  real*4 dminIn, dmaxin, dmeanIn
  logical*4 useFile, useAutodoc
  integer*4 iiuFileType, iiuRetAdocIndex, AdocGetImageMetaInfo
  !
  nbyteExtraIn = 0
  nbytePerSecIn = 0
  useFile = pieceFileIn .ne. ' ' .and. pieceFileIn .ne. 'none'
  !
  ! Use piece list as primary source if defined
  if (useFile) then
    call read_piece_list2(pieceFileIn, ixPcList, iyPcList, izPcList, numPcList, limPcList)
  endif
  !
  ! Open file and autodoc if indicated or for HDF file
  call ialprt(.false.)
  call imopen(4, inFile, 'RO')
  call irdhdr(4, nxyz, mxyz, mode, dminIn, dmaxin, dmeanIn)
  useAutodoc = useMdoc .or. iiuFileType(4) == 5
  if (useAutodoc) then
    indAdoc = iiuRetAdocIndex(4, 0, 1)
    if (indAdoc < 0) then
      write(*,'(/,a,a)') 'ERROR: EDMONT - Could not open autodoc information for '// &
          'input file ', trim(inFile)
      call exit(1)
    endif
    call setCurrentAdocOrExit(indAdoc, 'input')
    if (AdocGetImageMetaInfo(montage, numSect, iTypeAdoc) == 0 .and. .not.useFile) then
      call get_metadata_pieces(indAdoc, iTypeAdoc, nxyz(3), ixPcList, iyPcList, &
          izPcList, limPcList, numPcList)
    endif
    if (iiuFileType(4) .ne. 5) call AdocClear(indAdoc)
  endif
  !
  ! get extra header information if any
  !
  call iiuRetNumExtended(4, nbyteExtraIn)
  if (nbyteExtraIn > 0) then
    if (nbyteExtraIn > maxExtra .and. .not.useFile) then
      write(*,'(/,a,a,a)') 'ERROR: EDMONT - No piece list file was given for input '// &
          'file ', trim(inFile), ' and the extra header data are too large for the array'
      call exit(1)
    else
      call iiuRetExtendedType(4, nbytePerSecIn, iflagExtraIn)
      if (useFile) then
        if (mod(iflagExtraIn / 2, 2) == 0) nbyteExtraIn = 0
      else
        call iiuRetExtendedData(4, nbyteExtraIn, extraIn)
        call get_extra_header_pieces(extraIn, nbyteExtraIn, nbytePerSecIn, iflagExtraIn, &
            nxyz(3), ixPcList, iyPcList, izPcList, numPcList, limPcList)
      endif
    endif
  endif
  if (.not.(useFile .or. useAutodoc) .and. (nbyteExtraIn == 0 .or. numPcList == 0)) then
    write(*,'(/,a,a,a)')'ERROR: EDMONT - No piece list file was given for input file ', &
        trim(inFile), ' and the header does not contain piece coordinates'
    call exit(1)
  endif
  call iiuClose(4)
  return
end subroutine read_pl_or_header


! Open an image file, get piece coordinates one way or another, analyze
! the list of Z values in the file, and determine the montage
! characteristics in each dimension
!
subroutine openAndAnalyzeFiles(imageFile, pieceFile, nxyz, dmin2, printing, &
    mode, extraIn, maxExtra, nbyteExtraIn, nbytePerSecIn, ixPcList, iyPcList, &
    izPcList, numPcList, limPcList, listz, numInZlist, minXpiece , numXpieces, &
    nxOverlap, minYpiece, numYpieces, nyOverlap, useMdoc)
  implicit none
  character*(*) imageFile, pieceFile
  logical*4 printing, useMdoc
  integer*4 nxyz(3), mode, maxExtra, nbyteExtraIn, ixPcList(*), iyPcList(*)
  integer*4 izPcList(*), numPcList, limPcList, listz(*), numInZlist, minXpiece
  integer*4 numXpieces, nxOverlap, minYpiece, numYpieces, nyOverlap
  integer*1 extraIn(*)
  real*4 dmin2, dmax2, dmean2
  integer*4 mxyz(3), nbytePerSecIn
  !
  ! Read the piece data first in case there is a problem opening the
  ! same file on two channels
  call read_pl_or_header(pieceFile, imageFile, extraIn, maxExtra, nbyteExtraIn, &
      nbytePerSecIn, ixPcList, iyPcList, izPcList, numPcList, limPcList, useMdoc)
  call ialprt(printing)
  call imopen(1, imageFile, 'RO')
  call irdhdr(1, nxyz, mxyz, mode, dmin2, dmax2, dmean2)

  call fill_listz(izPcList, numPcList, listz, numInZlist)
  call checklist(ixPcList, numPcList, 1, nxyz(1), minXpiece , numXpieces, nxOverlap)
  call checklist(iyPcList, numPcList, 1, nxyz(2), minYpiece , numYpieces, nyOverlap)
  return
end subroutine openAndAnalyzeFiles

!
! Early history
! David Mastronarde for VAX       5/9/89
! 4/19/90 Added line-by-line prompts for files and section entry, added
! time/date stamp, fixed bugs on rescaling in place and missing secs
! 1999: added ability to knock out pieces.
! 1/3/00: made it handle extra header data, made scaling logic more
! like NEWSTACK and made sure it could handle negative integers.
! 10/24/00: made it actually use coordinates in header and renumber
! sections sequentially.
!
