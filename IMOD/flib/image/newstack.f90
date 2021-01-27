! ************* NEWSTACK **********************************************
!
! NEWSTACK is a general stack editor to move images into, out of, or
! between stacks.  It can float the images to a common range or mean of
! density. It can apply a general linear transform specified as a line
! in a file. It can put the output into a smaller or larger array and
! independently recenter each image separately from the transform.
! Images can be taken from multiple input files and placed into multiple
! output files.
!
! for all details see the man page.
!
! $Id$
!
program newstack
  implicit none
  integer MAXTEMP, maxChunks, LIMGRADSEC
  parameter (maxChunks = 250, LIMGRADSEC = 10000)
  parameter (MAXTEMP = 5000000)
  integer*4 nx, ny, nz
  real*4, allocatable :: array(:)
  !
  integer*4 nxyz(3), mxyz(3), nxyzst(3), nxyz2(3), mxyz2(3), maxExtraIn, maxExtraOut
  real*4 cell2(6), cell(6), title(20), delta(3), xOrigin, yOrigin, zOrigin, deltafirst(3)
  !
  character*320 xfFile, inFileList, outFileList
  character*320, allocatable :: inFile(:), outFile(:)
  character*320 idfFile, magGradFile, newAngleFile
  character*320 tempName, temp_filename, seriesExt
  character*100000 listString
  character*6 convFormat
  character*10 convNum, zvalueName, globalName
  equivalence (nx, nxyz(1)), (ny, nxyz(2)), (nz, nxyz(3))
  !
  data nxyzst/0, 0, 0/
  character*20 floatText/' '/, xfText/' '/, trunctText/' '/
  real*4 frot(2,3), fexp(2,3), fprod(2,3)
  integer*4, allocatable :: inList(:)
  integer*4, allocatable :: nlist(:), listInd(:), numSecOut(:), lineTmp(:), isecExclude(:)
  real*4 optimalMax(17)
  integer*4 lineOutSt(maxChunks+1), numLinesOut(maxChunks)
  integer*4 lineInSt(maxChunks+1), numLinesIn(maxChunks)
  real*4, allocatable :: scaleFacs(:), scaleConsts(:), secZmins(:), secZmaxes(:), ztemp(:)
  real*4, allocatable :: secMins(:), secMaxes(:)
  integer*1, allocatable :: extraIn(:), extraOut(:)
  integer*1 btiltTemp(4)
  integer*2 itiltTemp(2)
  real*4 rtiltTemp
  equivalence (rtiltTemp, btiltTemp), (itiltTemp, btiltTemp)
  data optimalMax/255., 32767., 255., 32767., 255., 255., 65535., 255., 255., &
      511., 1023., 2047., 4095., 8191., 16383., 32767., 255./
  !
  integer(kind = 8) idimInOut, limToAlloc, i8, numPix, iChunkBase, iBufOutBase, istart
  integer(kind = 8) numMove, moveOffset, limIfFail, loadBaseInd
  integer*4 ifDistort, idfBinning, iBinning, idfNx, idfNy, iWarpFlags
  integer*4 nxGrid, nyGrid, numFields, numIdfUse
  real*4 xGridIntrv, yGridIntrv, pixelSize, xGridStart, yGridStart, warpScale
  real*4 xnBig, ynbig
  !
  integer*4 ifMagGrad, numMagGrad, magUse
  real*4 pixelMagGrad, axisRot
  integer*4, allocatable :: lineUse(:), listReplace(:), idfUse(:), nControl(:)
  real*4, allocatable :: xcen(:), ycen(:), secMean(:), f(:,:,:), extraTilts(:)
  real*4, allocatable :: allFileTilts(:)
  real*4, allocatable :: tmpDx(:,:), tmpDy(:,:), fieldDx(:,:), fieldDy(:,:)
  real*4, allocatable :: warpXoffsets(:), warpYoffsets(:)
  integer*4, allocatable :: listVolumes(:), intForBigString(:), izNotPiece(:)
  real*4 sdChunk(maxChunks+1)
  real*8 pixChunk(maxChunks+1), dsumChunk(maxChunks+1)

  real*4 tiltAngles(LIMGRADSEC), dmagPerMicron(LIMGRADSEC), rotPerMicron(LIMGRADSEC)
  !
  logical rescale, blankOutput, adjustOrigin, hasWarp, fillTmp, fillNeeded, stripExtra
  logical readShrunk, numberedFromOne, twoDirections, useMdocFiles, outDocChanged, quiet
  logical saveTilts, serialEMtype, packed4bitInput, pack4bitOutput, noisePad, sizesMatch
  logical phaseShift, sameSecList, FEI1type, needHeaderAngles, needRangeFix, rangeFixOK
  logical scaleRangeFixOK, scaleFixAloneOK, needScaleFix, doRangeScale, doRangeOnly
  logical doScaleOnly, scaleFixShiftedOK, processInPlace, preSetScaling, inPlace
  character dat*9, timeStr*8, tempExt*9
  logical nbytes_and_flags
  character*80 titlech
  integer*4 inUnit, numInFiles, listTotal, numOutTot, numOutFiles, nxOut, nyOut, lmGrid
  integer*4 newMode, ifOffset, ifXform, nXforms, nLineUse, ifMean, ifFloat, ifWarping
  integer*4 nsum, ilist, iFile, iSecRead, loadYstart, loadYend, isec, isecOut, itype
  real*4 xOffsAll, yOffsAll, fracZero, dminSpecified, dmaxSpecified, contrastLo
  real*4 zmin, zmax, diffMinMean, diffMaxMean, grandSum, sdSec, contrastHi
  real*4 grandMean, shiftMin, shiftMax, shiftMean, dminIn, dmaxIn, dmeanIn
  integer*4 iOutFile, numTruncLow, numTruncHigh, ifHeaderOut, ifTempOpen, nByteSymIn
  integer*4 nByteExtraIn, iFlagExtraIn, mode, nByteExtraOut, nByteSymOut, indExtraOut
  real*4 dmin, dmax, dmean, dmin2, dmax2, dmean2, optimalIn, optimalOut, bottomIn
  real*4 bottomOut, xcenIn, ycenIn, dx, dy, fieldMaxX, ystart, readReduction
  integer*4 linesLeft, numChunks, nextLine, iChunk, ifOutChunk, iscan, iyTest, iVerbose
  integer*4 iyBase, iy1, iy2, lnu, maxin, numScaleFacs, maxFieldX, needYfirst, needYlast
  real*4 dmeanSec, tmpMin, tmpMax, tmpMin2, tmpMax2, val, scaleFactor, ftReduceFac
  integer*4 needYstart, needYend, numLinesLoad, numYload, numYchunk, iseriesBase, nyNeeded
  integer*4 ix1, ix2, nByteCopy, nByteClear, ifLinear, limEntered, insideTaper, indFile
  real*4 constAdd, densOutMin, dens, tmin2, tmax2, tmean2, avgSec, enteredSD, enteredMean
  integer*4 numInputFiles, numSecLists, numOutputFiles, numToGet, maxNxGrid, maxNyGrid
  integer*4 numOutValues, numOutEntries, ierr, ierr2, i, kti, iy, ind, numSecTrunc
  integer*4 maxFieldY, nxFirst, nyFirst, nxBin, nyBin, indGlobalAdoc
  integer*4 lenTemp, ierr3, applyFirst, numTaper, numberOffset, numExclude, ifChunkIn
  integer*4 ifOnePerFile, ifUseFill, listIncrement, indOut, ifMeanSdEntered, nzChunkIn
  integer*4 numReplace, isecReplace, modeOld, loadYoffset, listAlloc
  integer*4 indFilter, linesShrink, numAllSec, maxNumXF, nxMax, nyMax, ifControl, nzChunk
  integer*4 indFiltTemp, ifFiltSet, ifShrink, numVolRead, if3dVolumes, nxTile, nyTile
  integer*4 indAdocIn, indAdocOut, indSectIn, needClose1, needClose2, nxTileIn, nyTileIn
  integer*4 nxFSpad, nyFSpad, maxFSpad, minFSpad, nxDimNeed, nyDimNeed, numReverse
  integer*4 indArg, ifirstExtraType, ianyExtraType, numForTilt, iiuFlags, intStringSize
  integer*4 nxFCropPad, nyFCropPad, nxRepakDim, nyRepakDim, ioutBase, modeOfFirst
  integer*4 ifContrast, ifScaleMM, ifFloatDen
  real*4 rxOffset, ryOffset, fsPadFrac, actualFac, xstart, binningOfInput, rangeBase
  real*4 fieldMaxY, rotateAngle, expandFactor, fillVal, shrinkFactor, sdCritForScaleFix
  real*4 boostForMin, boostForMax, fixRangeSDs, rangeExpandFac, rangeAdd, fixScaleFac
  real*4 physicalMem
  integer*4 limSec, numIntOrBytesIn, maxInFileAngles, ifReorderByTilt, numRangeFixScans
  real*8 dsum, dsumSq, tsum, tsum2, tsumSq, wallStart, wallTime, loadTime, saveTime
  real*8 rotTime, taperTime
  !
  ! Function names
  real*4 cosd, sind
  integer*4 taperAtFill, selectZoomFilter, zoomFiltInterp, numberInList, niceFFTlimit
  integer*4 readCheckWarpFile, AdocLookupByNameValue, AdocTransferSection, AdocSetCurrent
  integer*4 AdocFindInsertIndex, AdocSetFloat, AdocInsertSection
  integer*4 iiuRetAdocIndex, iiuVolumeOpen, iiuAltChunkSizes, AdocWrite, niceFrame
  integer*4 getLinearTransform, findMaxGridSize, getSizeAdjustedGrid, iiuFileType
  integer*4 setOutputTypeFromString, b3dOutputFileType, AdocSetInteger, iiuWriteGlobalAdoc
  integer*4 fourierCropSizes, getExtraHeaderMaxSecSize, getExtraHeaderSecOffset
  integer*4 copyExtraHeaderSection
  real*8 b3dPhysicalMemory
  !
  logical pipinput
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetInteger, PipGetBoolean, PipGetLogical, PipGetThreeIntegers
  integer*4 PipGetString, PipGetTwoIntegers, PipGetFloatArray, PipGetFloat
  integer*4 PipGetIntegerArray, PipGetNonOptionArg, PipGetTwoFloats
  !
  ! fallbacks from ../../manpages/autodoc2man -3 2  newstack
  !
  integer numOptions
  parameter (numOptions = 66)
  character*(40 * numOptions) options(1)
  options(1) = &
      'input:InputFile:FNM:@output:OutputFile:FNM:@fileinlist:FileOfInputs:FN:@'// &
      'fileoutlist:FileOfOutputs:FN:@reverse:ReverseInputFileOrder:I:@'// &
      'split:SplitStartingNumber:I:@append:AppendExtension:CH:@'// &
      'format:FormatOfOutputFile:CH:@volumes:VolumesToRead:LI:@3d:Store3DVolumes:I:@'// &
      'chunk:ChunkSizesInXYZ:IT:@mdoc:UseMdocFiles:B:@tilt:TiltAngleFile:FN:@'// &
      'reorder:ReorderByTiltAngle:I:@angle:AngleFileToReorder:FN:@'// &
      'newangle:NewAngleOutputFile:I:@secs:SectionsToRead:LIM:@'// &
      'samesec:SameSectionsToRead:B:@fromone:NumberedFromOne:B:@'// &
      'exclude:ExcludeSections:LI:@twodir:TwoDirectionTiltSeries:B:@'// &
      'skip:SkipSectionIncrement:I:@numout:NumberToOutput:IAM:@'// &
      'replace:ReplaceSections:LI:@blank:BlankOutput:B:@offset:OffsetsInXandY:FAM:@'// &
      'applyfirst:ApplyOffsetsFirst:B:@xform:TransformFile:FN:@'// &
      'uselines:UseTransformLines:LIM:@onexform:OneTransformPerFile:B:@'// &
      'phase:PhaseShiftFFT:B:@rotate:RotateByAngle:F:@expand:ExpandByFactor:F:@'// &
      'shrink:ShrinkByFactor:F:@antialias:AntialiasFilter:I:@bin:BinByFactor:I:@'// &
      'ftreduce:FourierReduceByFactor:F:@noise:NoisePadForFFT:B:@'// &
      'distort:DistortionField:FN:@imagebinned:ImagesAreBinned:F:@'// &
      'fields:UseFields:LIM:@subarea:SubareaOffsetsXandY:FAM:@'// &
      'gradient:GradientFile:FN:@origin:AdjustOrigin:B:@'// &
      'linear:LinearInterpolation:B:@nearest:NearestNeighbor:B:@'// &
      'size:SizeToOutputInXandY:IP:@mode:ModeToOutput:I:@'// &
      'bytes:BytesSignedInOutput:I:@strip:StripExtraHeader:B:@'// &
      'float:FloatDensities:I:@meansd:MeanAndStandardDeviation:FP:@'// &
      'contrast:ContrastBlackWhite:IP:@scale:ScaleMinAndMax:FP:@'// &
      'multadd:MultiplyAndAdd:FPM:@fixrange:FixRangeIfNeeded:FP:@'// &
      'rfparam:RangeFixingParams:FP:@fill:FillValue:F:@taper:TaperAtFill:IP:@'// &
      'memory:MemoryLimit:I:@test:TestLimits:IP:@megasec:MaxMegaSections:I:@'// &
      'quiet:QuietOutput:B:@verbose:VerboseOutput:I:@param:ParameterFile:PF:@'// &
      'help:usage:B:'
  !
  ! Pip startup: set error, parse options, check help, set flag if used
  !
  call PipReadOrParseOptions(options, numOptions, 'newstack', &
      'ERROR: NEWSTACK - ', .true., 2, 2, 1, numOptArg, &
      numNonOptArg)
  pipinput = numOptArg + numNonOptArg > 0
  !
  ! defaults
  !
  inFileList = ' '
  outFileList = ' '
  numInputFiles = 0
  numOutputFiles = 0
  numSecLists = 0
  ifOnePerFile = 0
  ifDistort = 0
  ifMagGrad = 0
  ifMeanSdEntered = 0
  maxFieldY = 0
  idfFile = ' '
  magGradFile = ' '
  rotateAngle = 0.
  expandFactor = 0.
  iBinning = 1
  binningOfInput = 1
  lenTemp = MAXTEMP
  idimInOut = limToAlloc - 1
  limEntered = 0
  applyFirst = 0
  ifLinear = 0
  numScaleFacs = 0
  ifUseFill = 0
  blankOutput = .false.
  stripExtra = .false.
  adjustOrigin = .false.
  listIncrement = 1
  numReplace = 0
  iVerbose = 0
  limIfFail = 1950000000 / 4
  numTaper = 0
  insideTaper = 0
  loadTime = 0.
  saveTime = 0.
  rotTime = 0.
  taperTime = 0.
  maxExtraIn = 4
  maxExtraOut = 0
  indFilter = 6
  shrinkFactor = 1.
  linesShrink = 0
  iseriesBase = -1
  seriesExt = ' '
  maxNumXF = 40000
  ifWarping = 0
  ifControl = 0
  lmGrid = 200
  readReduction = 1.
  readShrunk = .false.
  numberedFromOne = .false.
  twoDirections = .false.
  numberOffset = 0
  numExclude = 0
  numSecTrunc = 0
  numVolRead = 0
  if3dVolumes = 0
  nxTile = 0
  nyTile = 0
  nzChunk = 1
  useMdocFiles = .false.
  quiet = .false.
  phaseShift = .false.
  maxFSpad = 100
  minFSpad = 10
  fsPadFrac = 0.1
  nxFSpad = 0
  nyFSpad = 0
  ftReduceFac = 0.
  noisePad = .false.
  numReverse = 0
  saveTilts = .false.
  ianyExtraType = 0
  needClose1 = 0
  needClose2 = 0
  sameSecList = .false.
  newAngleFile = ' '
  maxInFileAngles = 0
  needHeaderAngles = .false.
  processInPlace = .false.
  preSetScaling = .false.
  ifReorderByTilt = 0
  fixRangeSDs = 0.
  numRangeFixScans = 9
  fixScaleFac = 1.
  sdCritForScaleFix = 10.
  rangeExpandFac = 1.2
  limSec = 1000000
  limToAlloc = 4 * limSec
  physicalMem = b3dPhysicalMemory() / 4.

  if (pipInput) then
    if (PipGetInteger('MaxMegaSections', ierr) == 0) limSec = 1000000 * max(1, ierr)
  endif
  intStringSize = limSec * 8
  !
  ! Preliminary allocation of array
  allocate(array(limToAlloc), inList(limSec), stat = ierr)
  call memoryError(ierr, 'Allocating main array')
  call AdocGetStandardNames(globalName, zvalueName)
  !
  ! read in list of input files
  !
  call iiuAltPrint(0)
  inUnit = 5
  !
  ! get number of input files and other preliminary items
  !
  if (pipinput) then
    ierr = PipGetInteger('VerboseOutput', iVerbose)
    ierr = PipGetString('FileOfInputs', inFileList)
    call PipNumberOfEntries('InputFile', numInputFiles)
    numInFiles = numInputFiles + max(0, numNonOptArg - 1)
    if (numInFiles > 0 .and. inFileList .ne. ' ') call exitError( &
        'You cannot enter both input files and an input list file')
    if (inFileList .ne. ' ') numInFiles = -1
    call PipNumberOfEntries('SectionsToRead', numSecLists)
    ierr = PipGetInteger('SkipSectionIncrement', listIncrement)
    ierr = PipGetLogical('BlankOutput', blankOutput)
    ierr = PipGetLogical('StripExtraHeader', stripExtra)
    ierr = PipGetLogical('NumberedFromOne', numberedFromOne)
    ierr = PipGetLogical('TwoDirectionTiltSeries', twoDirections)
    if (numberedFromOne) numberOffset = 1
    ierr = PipGetLogical('SameSectionsToRead', sameSecList)
    if (sameSecList .and. numSecLists > 1) call exitError( &
        'You cannot enter -samesec with multiple section list entries')
    if (sameSecList .and. twoDirections) call exitError( &
        'You cannot enter -samesec with -twodir')
    if (PipGetInteger('ReverseInputFileOrder', numReverse) == 0) then
      if (numInFiles < 0)  &
          call exitError('You cannot enter -reverse with an input file list')
      if (abs(numReverse) > numInFiles) &
          call exitError('The entry to -reverse is bigger than the number of input files')
      if (numReverse == 0) numReverse = numInFiles
    endif
    if (PipGetInteger('BytesSignedInOutput', i) == 0) call overrideWriteBytes(i)
    if (PipGetString('ExcludeSections', listString) == 0) then
      call parseList2(listString, inList, numExclude, limSec)
    endif
  else
    write(*,'(1x,a,$)') '# of input files (or -1 to read list'// &
        ' of input files from file): '
    read(5,*) numInFiles
  endif
  
  allocate(isecExclude(max(1, numExclude)), stat = ierr)
  call memoryError(ierr, 'array for excluded sections')
  if (numExclude > 0) then
    isecExclude(1:numExclude) = inList(1:numExclude) - numberOffset
  endif
  
  !
  ! if it is negative, open a list file, set up input from 7
  !
  if (numInFiles == 0) call exitError('No input file specified')
  if (numInFiles < 0) then
    inUnit = 7
    if (.not.pipinput) then
      write(*,'(1x,a,$)') 'Name of input list file: '
      read(5, 101) inFileList
    endif
    call dopen(7, inFileList, 'ro', 'f')
    read(inUnit,*) numInFiles
  endif
  if (.not.pipinput .or. inUnit == 7) then
    allocate(intForBigString(intStringSize / 4 + 10), stat = ierr)
    call memoryError(ierr, 'Array for giant string')
  endif
  listTotal = 0
  numAllSec = 0
  allocate(inFile(numInFiles), nlist(numInFiles), listInd(numInFiles),  &
      lineTmp(numInFiles), scaleFacs(numInFiles), scaleConsts(numInFiles),  &
      listVolumes(numInFiles), extraIn(maxExtraIn), stat = ierr)
  call memoryError(ierr, 'Arrays for input files')
  !
  ! Get HDF and volume related options
  if (pipInput) then
    if (PipGetString('FormatOfOutputFile', listString) == 0) then
      ierr = setOutputTypeFromString(listString)
      if (ierr == -5) &
          call exitError('HDF files are not supported by this IMOD package')
      if (ierr == -6) &
          call exitError('JPEG files are not supported by this IMOD package')
      if (ierr < 0) call exitError('Unrecognized entry for output file format')
    endif
    if (PipGetString('VolumesToRead', listString) == 0) then
      call parseList2(listString, inList, numVolRead, limSec)
      numVolRead = min(numVolRead, numInFiles)
      listVolumes(1:numVolRead) = inList(1:numVolRead)
    endif
    ifChunkIn = 1 - PipGetThreeIntegers('ChunkSizesInXYZ', nxTile, nyTile, nzChunk)
    if (ifChunkIn > 0) if3dVolumes = 1
    ierr = PipGetInteger('Store3DVolumes', if3dVolumes)
    if (ifChunkIn > 0 .and. if3dVolumes < 0) call exitError( &
        'You cannot enter chunk sizes and forbid volume output with -3d -1')
    if (if3dVolumes > 0) call overrideOutputType(5)
    ierr = PipGetLogical('UseMdocFiles', useMdocFiles)
  endif
  !
  nxMax = 0
  nyMax = 0
  if (twoDirections .and. numInFiles .ne. 2)  &
      call exitError('There must be exactly two input files to use -twodir')
  !
  ! For pip input, get all the filenames now in case they need to be reversed
  if (pipinput .and. inUnit .ne. 7) then
    do indArg = 1, numInFiles
      indFile = indArg
      if (numReverse > 0 .and. indArg <= numReverse)  &
          indFile = numReverse + 1 - indArg
      if (numReverse < 0 .and. indArg > numInFiles + numReverse)  &
          indFile = numInFiles + numReverse + (numInFiles + 1 - indArg)
      if (indArg <= numInputFiles) then
        ierr = PipGetString('InputFile', inFile(indFile))
      else
        ierr = PipGetNonOptionArg(indArg - numInputFiles, inFile(indFile))
      endif
    enddo
  endif

  sizesMatch = .true.
  do indFile = 1, numInFiles
    !
    ! get the next filename if it wasn't gotten in previous loop
    if (.not. (pipinput .and. inUnit .ne. 7)) then
      if (inUnit .ne. 7) then
        if (numInFiles == 1) then
          write(*,'(1x,a,$)') 'Name of input file: '
        else
          write(*,'(1x,a,i3,a,$)') 'Name of input file #', indFile, ': '
        endif
      endif
      read(inUnit, 101) inFile(indFile)
    endif
    !
    ! open file to make sure it exists and get default section list
    !
    call openInputFile(indFile)
    call irdhdr(1, nxyz, mxyz, mode, dmin2, dmax2, dmean2)
    if (indFile == 1) then
      nxFirst = nx
      nyFirst = ny
      modeOfFirst = mode
      call iiuRetDelta(1, deltafirst)
      !
      ! Retain first input file volume structure in output if already doing HDF unless
      ! user said not  to; adopt its chunk sizes unless user entered them
      call iiuRetChunkSizes(1, nxTileIn, nyTileIn, nzChunkIn)
      if (nzChunkIn > 0 .and. b3dOutputFileType() == 5) then
        if (if3dVolumes == 0) if3dVolumes = 1
        if (ifChunkIn == 0 .and. if3dVolumes > 0) then
          nxTile = nxTileIn
          nyTile = nyTileIn
          nzChunk = nzChunkIn
        endif
      endif
    endif
    call iiuClose(1)
    if (needClose1 > 0) call iiuClose(needClose1)
    if (nx .ne. nxFirst .or. ny .ne. nyFirst) sizesMatch = .false.
    nxMax = max(nx, nxMax)
    nyMax = max(ny, nyMax)
    nlist(indFile) = nz
    do isec = 1, nz
      iy = isec - 1
      if (twoDirections .and. indFile == 1) iy = nz - isec
      inList(min(listTotal + isec, limSec)) = iy + numberOffset
    enddo
    numAllSec = numAllSec + nz
    !
    ! get section list
    !
    if (.not.pipinput .or. inUnit == 7) then
      if (inUnit .ne. 7) print *,'Enter list of sections to read from' &
          //' file (/ for all, 1st sec is 0; ranges OK)'
      ierr = listTotal - limSec
      call readBigList(inUnit, inList(listTotal + 1), nlist(indFile), ierr,  &
          intForBigString, intStringSize)
      if (ierr > 0) then
        if (ierr == 4) call exitError('There must be a readable section list after'// &
            ' each filename in list of input files')
        call exitError('Processing section list in list of input files')
      endif
    elseif (indFile <= numSecLists) then
      ierr = PipGetString('SectionsToRead', listString)
      if (ierr == 0 .and. twoDirections)  &
          call exitError('You cannot enter section lists with -twodir')
      call parseList2(listString, inList(listTotal + 1), nlist(indFile),  &
          limSec - listTotal)
      !
      ! If there is only one input file and more section lists, accumulate the lists now
      if (numInFiles == 1 .and. numSecLists > 1) then
        do i = 2, numSecLists
          ierr = PipGetString('SectionsToRead', listString)
          call parseList2(listString, inList(listTotal + 1 + nlist(indFile)), ierr, &
              limSec - listTotal - nlist(indFile))
          nlist(indFile) = nlist(indFile) + ierr
        enddo
      endif
    elseif (sameSecList .and. indFile > 1) then
      nlist(indFile) = nlist(1)
      inList(listTotal + 1 : listTotal +  nlist(1)) = inList(1 : nlist(1))
    endif
    !
    ! check list legality and whether excluded; copy over if not excluded
    !
    listInd(indFile) = listTotal + 1
    indOut = listInd(indFile)
    ix1 = max(1, listIncrement)
    ix2 = numberOffset
    if (sameSecList .and. indFile > 1) then
      ix1 = 1
      ix2 = 0
    endif
    do isec = listTotal + 1, listTotal + nlist(indFile), ix1
      inList(isec) = inList(isec) - ix2
      if (.not.blankOutput .and. &
          (inList(isec) < 0 .or. inList(isec) >= nz)) then
        write(*,'(/,a,i9,a,a)') 'ERROR: NEWSTACK -', inList(isec) + numberOffset, &
            ' is an illegal section number for ', trim(inFile(indFile))
        call exit(1)
      endif
      if (numberInList(inList(isec), isecExclude, numExclude, 0) == 0) then
        inList(indOut) = inList(isec)
        indOut = indOut + 1
      endif
    enddo
    nlist(indFile) = indOut - listInd(indFile)
    listTotal = listTotal + nlist(indFile)
  enddo
  close(7)
101 format(a)
  !
  maxNumXF = max(maxNumXF, numAllSec)
  listAlloc = listTotal + 10
  if (.not. pipinput .or. inUnit == 7) deallocate(intForBigString)
  allocate(lineUse(listAlloc), listReplace(listAlloc), idfUse(listAlloc),  &
      xcen(listAlloc), ycen(listAlloc), secMean(listAlloc), f(2, 3, maxNumXF), &
      stat = ierr)
  call memoryError(ierr, 'Arrays for input files')
  !
  ! read in list of output files
  !
  inUnit = 5
  numOutTot = 0
  !
  ! get number of output files
  !
  if (pipinput) then
    ierr = PipGetString('FileOfOutputs', outFileList)
    call PipNumberOfEntries('OutputFile', numOutputFiles)
    numOutFiles = numOutputFiles + min(1, numNonOptArg)
    if (numOutFiles > 0 .and. outFileList .ne. ' ') call exitError( &
        'You cannot enter both output files and an output'// &
        ' list file')
    if (outFileList .ne. ' ') numOutFiles = -1
    ierr = PipGetInteger('SplitStartingNumber', iseriesBase)
    if (iseriesBase >= 0 .and. numOutFiles .ne. 1) call exitError('There'// &
        ' must be only one output file name for series of numbered files')
    if (iseriesBase >= 0) numOutFiles = listTotal
    ierr = PipGetString('AppendExtension', seriesExt)
  else
    write(*,'(1x,a,$)') '# of output files (or -1 to read list'// &
        ' of output files from file): '
    read(5,*) numOutFiles
  endif
  if (numOutFiles == 0) call exitError('No output file specified')
  !
  if (numOutFiles > 0) then
    allocate(outFile(numOutFiles), numSecOut(numOutFiles), stat = ierr)
    call memoryError(ierr, 'Arrays for output files')
  endif
  !
  ! get list input
  !
  if (numOutFiles < 0) then
    inUnit = 7
    if (.not.pipinput) then
      write(*,'(1x,a,$)') 'Name of output list file: '
      read(5, 101) outFileList
    endif
    call dopen(7, outFileList, 'ro', 'f')
    read(inUnit,*) numOutFiles
    if (numOutFiles <= 0) call exitError('The output list file must start'// &
        ' with a positive number of files to output')
    allocate(outFile(numOutFiles), numSecOut(numOutFiles), stat = ierr)
    call memoryError(ierr, 'Arrays for output files')
  elseif (numOutFiles == 1 .and. .not.pipinput) then
    !
    ! get single output file
    !
    write(*,'(1x,a,$)') 'Name of output file: '
    read(5, 101) outFile(1)
    numSecOut(1) = listTotal
    numOutTot = listTotal
  endif
  !
  ! or get all the output files and the number of sections
  !
  if (numOutTot == 0) then
    if (pipinput .and. inUnit .ne. 7) then
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
      do i = 1, numOutFiles
        if (i <= numOutputFiles) then
          ierr = PipGetString('OutputFile', outFile(i))
        else
          ierr = PipGetNonOptionArg(numNonOptArg, outFile(i))
        endif
        numOutTot = numOutTot + numSecOut(i)
      enddo
    else
      do i = 1, numOutFiles
        if (inUnit .ne. 7) write(*,'(1x,a,i3,a,$)') &
            'Name of output file #', i, ': '
        read(inUnit, 101) outFile(i)
        if (inUnit .ne. 7) write(*,'(1x,a,$)') &
            'Number of sections to store in that file: '
        read(inUnit,*) numSecOut(i)
        numOutTot = numOutTot + numSecOut(i)
      enddo
    endif
  endif
  !
  ! if series, now take the one filename as root name and make filenames
  if (iseriesBase >= 0) then
    ierr = alog10(10. * (iseriesBase + listTotal - 1))
    tempName = outFile(1)
    write(convFormat, 132) ierr, ierr
132 format('(i',i1,'.',i1,')')
    do i = 1, listTotal
      numSecOut(i) = 1
      write(convNum, convFormat) i + iseriesBase - 1
      if (seriesExt == ' ') then
        outFile(i) = trim(tempName) //'.'//trim(adjustl(convNum))
      else
        outFile(i) = trim(tempName) //trim(adjustl(convNum)) //'.'//trim(seriesExt)
      endif
    enddo
    numOutTot = listTotal
  endif

  if (numOutTot .ne. listTotal) call exitError( &
      'Number of input and output sections does not match')
  !
  ! get new size and mode and offsets
  !
  nxOut = -1
  nyOut = -1
  newMode = modeOfFirst
  xOffsAll = 0.
  yOffsAll = 0.
  ifOffset = 0
  if (pipinput) then
    ierr = PipGetTwoIntegers('SizeToOutputInXandY', nxOut, nyOut)
    ierr = PipGetInteger('ModeToOutput', newMode)
    call getOffsetEntries('OffsetsInXandY', 'offset', xcen, ycen, ifOffset)
  else
    write(*,'(1x,a,$)') 'Output file X and Y dimensions'// &
        ' (/ for same as first input file): '
    read(5,*) nxOut, nyOut
    write(*,'(1x,a,$)') 'Output file data mode (/ for same as first input file): '
    read(5,*) newMode
    !
    ! get list of x, y coordinate offsets
    !
    write(*,'(1x,a,/,a,$)') '1 to offset centers of individual '// &
        'images,', '  -1 to apply same offset to all sections,'// &
        ' or 0 for no offsets: '
    read(5,*) ifOffset
    if (ifOffset > 0) then
      print *,'Enter X and Y center offsets for each section'
      read(5,*) (xcen(i), ycen(i), i = 1, listTotal)
    elseif (ifOffset < 0) then
      write(*,'(1x,a,$)') 'X and Y center offsets for all sections: '
      read(5,*) xOffsAll, yOffsAll
    endif
    !
    ! fill offset list if only one or none
    if (ifOffset <= 0) then
      do i = 1, listTotal
        xcen(i) = xOffsAll
        ycen(i) = yOffsAll
      enddo
    endif
  endif

  if (newMode == 16) call exitError('Cannot output'// &
      ' color data (mode 16); use colornewst instead or set another output mode')

  ! set up for making mode 101 4-bit output
  pack4bitOutput = newMode == 101
  if (pack4bitOutput) then
    newMode = 0
    call set4BitOutputMode(1)
  endif
  !
  ! Get list of transforms
  ! ifXform initially indicates if transforms entered, including warping ones, but is
  !   set to 1 for rotating expanding, distorting, or applying mag gradients
  !   It is set back to 0 if Fourier shift or reduce is used (phaseShift, ftReduceFac)
  ! ifWarping is set if the transforms are warps
  ! If ShrinkByFactor, the entry is translated to readShrunk if it happens on read-in
  !   or its inverse to the same expandFactor used for expand, with linesShrunk as the
  !   flag for post-read shrinking
  ifXform = 0
  if (pipinput) then
    ierr = PipGetBoolean('LinearInterpolation', ifLinear)
    if (PipGetString('TransformFile', xfFile) == 0) ifXform = 1
    ierr = PipGetBoolean('OneTransformPerFile', ifOnePerFile)
    ix1 = 0
    ierr = PipGetBoolean('NearestNeighbor', ix1)
    if (ix1 .ne. 0 .and. ifLinear .ne. 0) call exitError( &
        'You cannot enter both -linear and -nearest')
    if (ix1 .ne. 0) ifLinear = -1
    ierr = PipGetLogical('PhaseShiftFFT', phaseShift)
    ierr = PipGetFloat('FourierReduceByFactor', ftReduceFac)
    ierr = PipGetLogical('NoisePadForFFT', noisePad)
    if ((phaseShift .or. ftReduceFac > 0 ) .and. ifLinear .ne. 0) call exitError( &
        'You cannot enter -phase or -ftreduce with -linear or -nearest')
    if (ftReduceFac > 0 .and. ftReduceFac <= 1.) call exitError( &
        'Factor for reducing with FFT must be > 1')
    if (ftReduceFac > 0) phaseShift = .false.
  else
    write(*,'(1x,a,$)') '1 or 2 to transform images with cubic or' &
        //' linear interpolation, 0 not to: '
    read(5,*) ifXform
    if (ifXform .ne. 0) then
      write(*,'(1x,a,$)') 'Name of transform file: '
      read(5, 101) xfFile
    endif
    if (ifXform > 1) ifLinear = 1
  endif
  if (ifXform .ne. 0) then
    !
    ! read transforms, set default to section list unless there is just
    ! one line in file
    !
    ierr = readCheckWarpFile(xfFile, 0, 1, idfNx, idfNy, nXforms, &
        idfBinning, pixelSize,  iWarpFlags, listString)
    if (ierr < -1) call exitError(listString)
    if (ierr >= 0) then
      write(*,'(a,a)') 'Warping file opened: ', trim(xfFile)
      ifWarping = 1
      if (mod(iWarpFlags / 2, 2) .ne. 0) ifControl = 1
      if (nXforms > maxNumXF) call exitError( &
          'Too many sections in warping file for transform array')
      warpScale = pixelSize / deltafirst(1)
      do i = 1, nXforms
        if (getLinearTransform(i, f(1, 1, i)) .ne. 0) &
            call exitError('Getting linear transform from warp file')
        f(1, 3, i) = f(1, 3, i) * warpScale
        f(2, 3, i) = f(2, 3, i) * warpScale
      enddo
      numFields = nXforms
    else
      call dopen(3, xfFile, 'ro', 'f')
      call xfrdall2(3, f, nXforms, maxNumXF, ierr)
      if (ierr == 2) call exitError('Reading transform file')
      if (ierr == 1) call exitError( &
          'Too many transforms in file for transform array')
      close(3)
    endif
    if (nXforms == 0) call exitError('The transform file contains no transforms')

    call getItemsToUse(nXforms, listTotal, inList, 'UseTransformLines', &
        listString, pipinput, 'TRANSFORM LINE', ifOnePerFile, numInFiles, &
        lineUse, nLineUse, numberOffset, listAlloc)

    if (ifOnePerFile > 0) then
      if (nLineUse < numInFiles) call exitError( &
          'Not enough transforms specified for the input files')
      !
      ! Copy list to temp array and build list with line for each sec
      !
      do iy = 1, numInFiles
        lineTmp(iy) = lineUse(iy)
      enddo
      nLineUse = 0
      do iy = 1, numInFiles
        do i = 1, nlist(iy)
          nLineUse = nLineUse + 1
          lineUse(nLineUse) = lineTmp(iy)
        enddo
      enddo
    endif
    !
    ! use single number for all sections
    !
    xfText = ', transformed'
    if (ifWarping .ne. 0) xfText = ', warped'
    if (nLineUse == 1) then
      do i = 2, listTotal
        lineUse(i) = lineUse(1)
      enddo
      nLineUse = listTotal
    endif
    if (nLineUse .ne. listTotal) call exitError( &
        'Specified # of transform lines does not match # of sections')
  endif
  !
  ! find out if float or other density modification
  ! ifFloat ends up with the entered value (0 to 4), or 2 if mean/sd entered,
  !   -1 for scaleminmax, multadd or range fixing, -2 for contrast
  ! ifMean is 1 for ifFloat > 1
  ! dminSpecified, dmaxSpecified are the new min and max for -scale
  !   the contrast entry is immediately converted to this: namely the min and max values
  !   relative to a byte range that the input file min and max should map to
  ! ifMeanSdEntered is set if that is how the target mean/sd is determined
  fracZero = 0.
  ifMean = 0
  ifFloat = 0
  dminSpecified = 0
  dmaxSpecified = 0
  contrastLo = 0
  contrastHi = 255
  if (pipinput) then
    ifContrast = 1 - PipGetTwoFloats('ContrastBlackWhite', contrastLo, contrastHi)
    ifScaleMM = 1 - PipGetTwoFloats('ScaleMinAndMax', dminSpecified, dmaxSpecified)
    ifFloatDen = 1 - PipGetInteger('FloatDensities', ifFloat)
    ifMeanSdEntered = 1 - PipGetTwoFloats('MeanAndStandardDeviation', enteredMean, &
        enteredSD)
    call PipNumberOfEntries('MultiplyAndAdd', numScaleFacs)
    if (ifMeanSdEntered .ne. 0) then
      if ((ifFloatDen .ne. 0 .and. ifFloat .ne. 2) .or.  &
          ifContrast + ifScaleMM + numScaleFacs > 0)  &
          call exitError('You cannot use -meansd with any scaling option except -float 2')
      ifFloatDen = 1
      ifFloat = 2
    endif
    ifUseFill = 1 - PipGetFloat('FillValue', fillVal)
    if (ifFloat >= 4) then
      if (ifScaleMM == 0) call exitError('You must enter -scale with -float 4')
    else
      if (ifContrast + ifScaleMM + ifFloatDen + min(numScaleFacs, 1) > 1) &
          call exitError('The -scale, -contrast, -multadd, and -float '// &
          'options are mutually exclusive except with -float 4')
      if (ifFloat < 0) call exitError('You must use -contrast or ' &
          //'-scale instead of a negative -float entry')
      if (ifContrast .ne. 0) ifFloat = -2
      if (ifScaleMM .ne. 0 .or. numScaleFacs .ne. 0) ifFloat = -1
      !
      ! get scale factors, make sure there are right number
      !
      if (numScaleFacs > 0) then
        if (numScaleFacs .ne. 1 .and. numScaleFacs .ne. numInFiles) &
            call exitError('You must enter -multadd either once '// &
            'or once per input file')
        do i = 1, numScaleFacs
          ierr = PipGetTwoFloats('MultiplyAndAdd', scaleFacs(i), &
              scaleConsts(i))
        enddo
      endif
    endif
    !
    ! Get entry for fixing ranges by shifting or changing mode
    if (PipGetTwoFloats('FixRangeIfNeeded', fixRangeSDs, fixScaleFac) == 0) then
      if (ifFloat .ne. 0 .or. numScaleFacs > 0)  &
          call exitError('You cannot enter -fixrange with any scaling options')
      if (fixRangeSDs .ne. 0 .and. fixRangeSDs < 2)  &
          call exitError('The entry for -fixrange must be at least 2')
      if (numInFiles > 1) &
          call exitError('You cannot use -fixrange with more than one input file')
      if (newMode .ne. modeOfFirst)  &
          call exitError('You cannot use -fixrange if you enter an output mode')
      rangeFixOK =  &
          (modeOfFirst == 1 .and. dmean2 < 0) .or. (modeOfFirst == 6 .and. dmean2 <16000)
      scaleFixAloneOK = fixScaleFac > 1. .and. dmean2 >= 0.
      scaleFixShiftedOK = fixScaleFac > 1. .and. dmean2 < 0.
      scaleRangeFixOK = fixScaleFac > 1.
      if (.not. (rangeFixOK .or. scaleFixAloneOK .or. scaleRangeFixOK)) fixRangeSDs = 0.
      needRangeFix = .false.
      needScaleFix = .false.
      ierr = PipGetTwoFloats('RangeFixingParams', sdCritForScaleFix, rangeExpandFac)
    endif
  else
    write(*,102)
102 format(' Enter 0 for no floating',/,8x, &
        '-2 to scale to bytes based on black and white contrast ', &
        'levels',/,8x, &
        '-1 to specify a single rescaling of all sections' &
        ,/,8x,' 1 to float all sections to same range',/,8x,' ', &
        '2 to float all sections to same mean & standard deviation' &
        ,/,8x,' 3 to shift sections to same mean without scaling',/ &
        ,6x,'or 4 to shift to same mean and specify a single', &
        ' rescaling: ',$)
    read(5,*) ifFloat
  endif
  if (ifFloat > 1) ifMean = 1
  if (ifFloat < 0) then
    floatText = ', densities scaled'
    if (ifFloat == -1 .and. .not.pipinput) then
      write(*,'(1x,a,/,a,$)') 'Values to scale input file''s'// &
          ' min and max to,', '   or / to scale to maximum range,' &
          //' or 1,1 to override mode scaling: '
      read(5,*) dminSpecified, dmaxSpecified
    elseif (ifFloat < -1) then
      if (.not.pipinput) then
        write(*,'(1x,a,$)') 'Contrast ramp black and white settings ' &
            //'(values between 0 and 255): '
        read(5,*) contrastLo, contrastHi
      endif
      contrastHi = max(contrastHi, contrastLo + 1.)
      dminSpecified = -contrastLo * 255 / (contrastHi - contrastLo)
      dmaxSpecified = dminSpecified + 65025 / (contrastHi - contrastLo)
      ! print *,contrastLo, contrastHi, dminSpecified, dmaxSpecified
    endif
  endif
  !
  ! get new options
  !
  if (pipinput) then
    ierr = PipGetBoolean('ApplyOffsetsFirst', applyFirst)
    ierr = PipGetFloat('RotateByAngle', rotateAngle)
    ierr = PipGetFloat('ExpandByFactor', expandFactor)
    ierr = PipGetInteger('BinByFactor', iBinning)
    ierr = PipGetString('DistortionField', idfFile)
    ierr = PipGetString('GradientFile', magGradFile)
    ierr = PipGetLogical('AdjustOrigin', adjustOrigin)
    ierr = PipGetTwoIntegers('TaperAtFill', numTaper, insideTaper)
    ierr = PipGetLogical('QuietOutput', quiet)
    ierr = PipGetInteger('ReorderByTiltAngle', ifReorderByTilt)
    !
    ! Tilt angles.  Allocate arrays for tilt angles to insert and to reorder if doing
    ! either operation
    ix1 = PipGetString('TiltAngleFile', tempName)
    if (ix1 == 0 .and. ifReorderByTilt .ne. 0) call exitError('You cannot enter '// &
        '-tilt with angles to insert and -reorder to reorder by angle')
    if (ifReorderByTilt .ne. 0) then
      if (numOutFiles > 1) call exitError('You cannot use -reorder with more than '// &
          'one output file')
      if (blankOutput .or. twoDirections) call exitError('You cannot use -reorder '// &
          'with the -blank or -twodir option')
      ierr = PipGetString('AngleFileToReorder', tempName)
      needHeaderAngles = ierr .ne. 0
      ierr = PipGetString('NewAngleOutputFile', newAngleFile)
    endif
    if (ix1 == 0 .or. ifReorderByTilt .ne. 0) then
      allocate(extraTilts(listTotal), stat = ierr)
      call memoryError(ierr, 'Arrays for tilt angles')
    endif
    !
    ! Now read the angles from either kind of file
    if (ix1 == 0 .or. (ifReorderByTilt .ne. 0 .and. .not. needHeaderAngles)) then
      if (ix1 == 0) then
        if (stripExtra) call exitError('You cannot enter both -tilt and -strip')
        if (index(tempName, '.', .true.) == 1) then
          ind = index(inFile(1), '.', .true.) - 1
          if (ind <= 0) ind = len_trim(inFile(1))
          tempName = inFile(1)(1:ind) // tempName
        endif
        saveTilts = .true.
      else
        saveTilts = abs(ifReorderByTilt) > 1
      endif
      call dopen(3, tempName, 'ro', 'f')
      do ind = 1, listTotal
        read(3, *, iostat=ierr) extraTilts(ind)
        if (ierr .ne. 0) call exitError('Reading tilt angle file: it must have as '// &
            'many lines as sections being written')
      enddo
      close(3)
    endif
    !
    ! Memory limits
    limEntered = 1 - PipGetTwoIntegers('TestLimits', ierr, lenTemp)
    if (limEntered > 0) limToAlloc = ierr
    if (PipGetInteger('MemoryLimit', ierr) == 0) then
      limEntered = 2
      limToAlloc = int(ierr, kind = 8) * 1024 * 256
    endif
    if (limEntered > 0) then
      if (limToAlloc < 1000 .or. lenTemp < 1 .or. lenTemp > &
          limToAlloc / 2) call exitError('Inappropriate memory limits entered')
      idimInOut = limToAlloc - lenTemp
      call reallocateArray()
    endif
    !
    if (ifWarping .ne. 0 .and. (idfFile .ne. ' ' .or. magGradFile .ne. ' ')) call &
        exitError('You cannot use distortion corrections with warping transforms')
    !
    if (ifWarping .ne. 0 .and. (rotateAngle .ne. 0 .or. expandFactor .ne. 0.)) call &
        exitError('You cannot use -expand or -rotate with warping transforms')
    if (iBinning <= 0) call exitError('Binning factor must be a positive number')
    readReduction = iBinning
    !
    ! Get filter entry and if there is binning and no separate shrink entry, convert
    ! the binning to a shrinkage.  Also allow a negative filter entry to set the default
    indFiltTemp = indFilter
    ifFiltSet = 1 - PipGetInteger('AntialiasFilter', indFiltTemp)
    if (indFiltTemp < 0) indFiltTemp = indFilter
    ifShrink = 1 - PipGetFloat('ShrinkByFactor', shrinkFactor)
    if (ifFiltSet > 0 .and. ifShrink == 0 .and. iBinning > 1 .and. indFiltTemp > 0) then
      shrinkFactor = iBinning
      iBinning = 1
      print *,'Doing antialias-filtered image reduction instead of ordinary binning'
    endif
    indFilter = max(0, indFiltTemp - 1)
    !
    ! Handle shrinkage
    if (ifShrink > 0 .or. shrinkFactor > 1.) then

      !
      ! Do shrinkage on input unless there is binning specified, since this will be
      ! more memory-efficient by default and it will produce a correct origin by default
      ! with no size change
      readShrunk = iBinning == 1
      if (iBinning > 1 .and. (ifXform .ne. 0 .or. rotateAngle .ne. 0. .or. &
          expandFactor .ne. 0. .or.  idfFile .ne. ' ' .or. magGradFile .ne. ' ' .or. &
          ifWarping .ne. 0)) call exitError('You cannot use both -shrink '// &
          'and -bin with -xform, -rotate, -expand, -distort, or -gradient')
      if (shrinkFactor <= 1.) call exitError('Factor for -shrink must be greater than 1')
      if (ifWarping .ne. 0 .and. abs(nint(shrinkFactor) - shrinkFactor) > 1.e-4) call  &
          exitError('You cannot use -shrink with warping unless the factor is an integer')
      ierr = 1
      i = indFilter
      do while (ierr == 1)
        ierr = selectZoomFilter(indFilter, 1. / shrinkFactor, linesShrink)
        if (ierr == 1) indFilter = indFilter - 1
        if (ierr > 1) call exitError( 'Selecting antialiasing filter')
      enddo
      if (indFilter < i) print *,'Using the last antialiasing filter, #', indFilter + 1
      if (readShrunk) then
        readReduction = shrinkFactor
        linesShrink = 0
      else
        !
        ! Post-read shrinkage: provide extra buffer of what needs reading in for a chunk
        ! and set an expansion factor
        linesShrink = linesShrink / 2 + 2
        expandFactor = 1. / shrinkFactor
      endif
      if (iVerbose > 0) print *,'Shrinking; readShrunk ', readShrunk
    endif
    !
    ! Check validity of phase shifting now that all requested actions are processed
    if (phaseShift .or. ftReduceFac > 0) then
      if (.not. sizesMatch) call exitError( &
          'All input files must have the same size in x and y with -phase or -ftreduce')
      if (ftReduceFac > 0 .and. (rotateAngle .ne. 0. .or. expandFactor .ne. 0. .or. &
          idfFile .ne. ' ' .or. magGradFile .ne. ' ' .or. ifWarping .ne. 0 .or. &
          readReduction > 1.))  &
          call exitError('You cannot use -ftreduce with -rotate, -expand, -distort, '// &
          '-gradient, -shrink, -bin OR warping')

      if (phaseShift .and. (rotateAngle .ne. 0. .or. expandFactor .ne. 0. .or.  &
          idfFile .ne. ' ' .or. magGradFile .ne. ' ' .or. ifWarping .ne. 0))  &
          call exitError('You cannot use -phase with -rotate, -expand, -distort, '// &
          '-gradient, or warping, or with -shrink unless it is the only other operation')
      if (applyFirst .ne. 0)  &
          call exitError('You cannot use -applyfirst with -phase or -ftreduce')
      !
      ! Check transforms and convert to xcen/ycen shifts
      if (ifXform .ne. 0) then
        tmpMax = .01 / max(nx, ny)
        do i = 1, listTotal
          lnu = lineUse(i) + 1
          if (abs(f(1, 1, lnu) - 1.) > tmpMax .or. abs(f(2, 2, lnu) - 1.) > tmpMax .or. &
              abs(f(2, 1, lnu)) > tmpMax .or. abs(f(1, 2, lnu)) > tmpMax) &
              call exitError('With -phase or -reduce, transforms must contain only '// &
              'shifts (first four terms must be 1 0 0 1)')
          xcen(i) = xcen(i) - f(1, 3, lnu)
          ycen(i) = ycen(i) - f(2, 3, lnu)
        enddo
        ifXform = 0
      endif
      !
      ! Set up for fourier cropping
      if (ftReduceFac > 0) then
        ierr = fourierCropSizes(nxFirst, ftReduceFac, 0.01, 16, niceFFTlimit(), nxFSpad, &
            nxFCropPad, actualFac)
        if (ierr > 0) call exitError('Reduction factor must be an integer or integer '// &
            'divided by 2, 3, 4, 5, 6, 8, or 10 to 3 decimal places')
        ierr = fourierCropSizes(nyFirst, ftReduceFac, 0.01, 16, niceFFTlimit(), nyFSpad, &
            nyFCropPad, actualFac)
        if (iVerbose > 0) print *,'Actual factor',actualFac
        !
        ! Reduce the shifts here for use when extracting image area; the fractional part
        ! will be boosted back up for the shift in the crop
        do i = 1, listTotal
          xcen(i) = xcen(i) / actualFac
          ycen(i) = ycen(i) / actualFac
        enddo
      endif
    endif
    !
    ! Section replacement
    if (PipGetString('ReplaceSections', listString) == 0) then
      call parseList2(listString, listReplace, numReplace, listTotal)
      if (numReplace > 0) then
        ! print *,'replacing', (listReplace(i), i = 1, numReplace)
        if (numOutFiles > 1) call exitError( &
            'There must be only one output file to use -replace')
        if (if3dVolumes > 0) call exitError('You cannot use -3d or -chunk with -replace')
        if (saveTilts .or. ifReorderByTilt .ne. 0) call exitError('You cannot use '// &
            '-tilt or -reorder with -replace')
        if (.not. quiet) call iiuAltPrint(1)
        call iiAllowMultiVolume(0)
        call imopen(2, outFile(1), 'OLD')
        call irdhdr(2, nxyz2, mxyz2, modeOld, dmin, dmax, dmean)
        call iiuAltPrint(0)
        call iiuFileInfo(2, ix1, ix2, iiuFlags)
        pack4bitOutput = btest(iiuFlags, 5) .or. btest(iiuFlags, 6)
        do i = 1, numReplace
          listReplace(i) = listReplace(i) - numberOffset
          if (listReplace(i) < 0 .or. listReplace(i) >= nxyz2(3)) &
              call exitError('Replacement section number out of range')
        enddo
      endif
    endif
    !
    ! Distortion field
    if (idfFile .ne. ' ') then
      ifDistort = 1
      xfText = ', undistorted'

      ierr = readCheckWarpFile(idfFile, 1, 1, idfNx, idfNy, numFields, &
          idfBinning, pixelSize, iWarpFlags, listString)
      if (ierr < 0) call exitError(listString)
      !
      if (PipGetFloat('ImagesAreBinned', binningOfInput) .ne. 0) then
        if (nxFirst <= idfNx * idfBinning / 2 .and. &
            nyFirst <= idfNy * idfBinning / 2) call exitError &
            ('you must specify current binning of images with -imagebinned because '// &
            'they are not larger than half the camera size')
      endif
      if (binningOfInput <= 0) call exitError('Image binning must be a positive number')
      warpScale = float(idfBinning) / binningOfInput
      !
      ! Set up default field numbers to use then process use list if any
      !
      call getItemsToUse(numFields, listTotal, inList, 'UseFields', listString, &
          pipinput, 'FIELD', 0, 0, idfUse, numIdfUse, numberOffset, listAlloc)
      if (numIdfUse == 1) then
        do i = 2, listTotal
          idfUse(i) = idfUse(1)
        enddo
        numIdfUse = listTotal
      endif
      if (numIdfUse .ne. listTotal) call exitError( &
          'Specified # of fields does not match # of sections')

    endif
    !
    ! get mag gradient information; multiply pixel size by binning
    !
    if (magGradFile .ne. ' ') then
      ifMagGrad = 1
      xfText = ', undistorted'
      call readMagGradients(magGradFile, LIMGRADSEC, pixelMagGrad, axisRot, &
          tiltAngles, dmagPerMicron, rotPerMicron, numMagGrad)
      pixelMagGrad = pixelMagGrad * readReduction
    endif
    !
    ! Get offsets from center of image relative to warping or distortion field    
    if (ifDistort > 0) then
      allocate(warpXoffsets(listAlloc), warpYoffsets(listAlloc), stat = ierr)
      call memoryError(ierr, 'Arrays for subarea offsets')
      warpXoffsets(:) = 0.
      warpYoffsets(:) = 0.
      call getOffsetEntries('SubareaOffsetsXandY', 'subarea offset', warpXoffsets, &
          warpYoffsets, i)
    endif
  endif
  call PipDone()
  !
  ! if not transforming and distorting, rotating, or expanding, set up
  ! a unit transform
  !
  if (ifXform == 0 .and. (ifDistort .ne. 0 .or. ifMagGrad .ne. 0 .or. &
      rotateAngle .ne. 0. .or. expandFactor .ne. 0.)) then
    ifXform = 1
    call xfunit(f(1, 1, 1), 1.)
    do i = 1, listTotal
      lineUse(i) = 0
    enddo
    nLineUse = listTotal
    nXforms = 1
    do while (rotateAngle > 180.01 .or. rotateAngle < -180.01)
      rotateAngle = rotateAngle - sign(360., rotateAngle)
    enddo
  endif
  !
  ! set up rotation and expansion transforms and multiply by transforms
  !
  if (rotateAngle .ne. 0. .or. expandFactor .ne. 0.) then
    call xfunit(frot, 1.)
    if (rotateAngle .ne. 0.) then
      frot(1, 1) = cosd(rotateAngle)
      frot(1, 2) = -sind(rotateAngle)
      frot(2, 2) = frot(1, 1)
      frot(2, 1) = -frot(1, 2)
      ! This was needed to correct for cubInterp rotating off center,
      ! changed 10/12/07
      ! frot(1, 3) = 0.5 * (frot(1, 1) + frot(1, 2)) - 0.5
      ! frot(2, 3) = 0.5 * (frot(2, 1) + frot(2, 2)) - 0.5
    endif
    if (expandFactor == 0.) expandFactor = 1.
    call xfunit(fexp, expandFactor)
    call xfmult(frot, fexp, fprod)
    ! print *,'transform', ((fprod(i, iy), i=1, 2), iy=1, 3)
    do i = 1, nXforms
      call xfmult(f(1, 1, i), fprod, frot)
      call xfcopy(frot, f(1, 1, i))
    enddo
  endif
  if (ftReduceFac > 1.) expandFactor = 1. / actualFac
  if (expandFactor == 0.) expandFactor = 1.
  !
  ! Cancel range-fixing if no interpolations are being done
  if (ifXform == 0 .and. ifDistort == 0 .and. ifMagGrad == 0 .and. .not. phaseShift)  &
      fixRangeSDs = 0.
  !
  ! adjust xcen, ycen and transforms if binning and allocate temp space
  ! NOTE: THE MULTIPLE SETTINGS OF idimInOut ARE ATROCIOUS AND NEED TO BE RATIONALIZED
  if (readReduction > 1.) then
    do i = 1, listTotal
      xcen(i) = xcen(i) / readReduction
      ycen(i) = ycen(i) / readReduction
    enddo
    if (ifXform .ne. 0) then
      do i = 1, nXforms
        f(1, 3, i) = f(1, 3, i) / readReduction
        f(2, 3, i) = f(2, 3, i) / readReduction
      enddo
    endif
    idimInOut = limToAlloc - lenTemp
  else if (ftReduceFac > 0. .and. nxFSpad > 0) then
    idimInOut = limToAlloc - (nxFCropPad + 2) * nyFCropPad
  else
    idimInOut = limToAlloc - 1
  endif
  
  ! For any kind of floating or fixing range, set some parameters and title entry
  if (ifFloat > 0 .or. fixRangeSDs > 0) then
    if (ifFloat > 0) then
      floatText = ', floated to range'
      if (fracZero .ne. 0.) &
          write(trunctText, '(a,f6.3)') ', truncated by', fracZero
      if (ifMean .ne. 0) then
        if (ifFloat == 2) then
          floatText = ', floated to means'
          zmin = 1.e10
          zmax = -1.e10
          allocate(secZmins(numOutTot), secZmaxes(numOutTot), ztemp(numOutTot),  &
              stat = ierr)
          call memoryError(ierr, 'Arrays for Z min/max')
        else
          diffMinMean = 0.
          diffMaxMean = 0.
          grandSum = 0.
          nsum = 0.
          if (ifFloat == 3) then
            floatText = ',  shifted to means'
          else
            floatText = ', mean shift&scaled'
            if (.not.pipinput) then
              write(*,'(1x,a,/,a,$)') 'Values to scale the shifted'// &
                  ' min and max to,', '   or / to scale to maximum'// &
                  ' range: '
              read(5,*) dminSpecified, dmaxSpecified
            endif
          endif
        endif
      allocate(secMins(numOutTot), secMaxes(numOutTot), stat = ierr)
      call memoryError(ierr, 'ARRAYS FOR MINS/MAXES')
      endif
    endif
      !
      ! if means and/or SDs needed, now read all sections to get them
      !
    if ((ifMean .ne. 0 .and. ifMeanSdEntered == 0) .or. fixRangeSDs > 0) then
      do iFile = 1, numInFiles
        call openInputFile(iFile)
        call irdhdr(1, nxyz, mxyz, mode, dmin2, dmax2, dmean2)
        !
        ! get the binned size to read
        !
        call getReducedSize(nx, readReduction, readShrunk, nxBin, rxOffset)
        call getReducedSize(ny, readReduction, readShrunk, nyBin, ryOffset)
        !
        nyNeeded = nyBin
        call reallocateIfNeeded()
        !
        do ilist = 1, nlist(iFile)
          ind = ilist + listInd(iFile) - 1
          iSecRead = inList(ind)
          if (iSecRead >= 0 .and. iSecRead < nz .and. (fixRangeSDs <= 0. .or. &
              (mod(ind, nlist(iFile) / min(nlist(iFile), numRangeFixScans)) == 0 .and. &
              (rangeFixOK .or. scaleRangeFixOK .or. scaleFixAloneOK .or. &
              scaleFixShiftedOK)))) then
            !
            if (iVerbose > 0) print *,'scanning for mean/sd', iSecRead
            call scanSection(array, idimInOut, nxBin, nyBin, 0, readReduction, rxOffset, &
                ryOffset, iSecRead, ifFloat, dmin2, dmax2, dmean2, sdSec, loadYstart, &
                loadYend, indFilter, readShrunk, fixRangeSDs, array(idimInOut + 1), &
                lenTemp)
            !
            ! Evaluate range fixing and skip the rest; first get the min and max values 
            ! to be considered, expanded range limited by a number of SDs from mean
            if (fixRangeSDs > 0) then
              rangeAdd = (rangeExpandFac - 1.) * (dmax2 - dmin2)
              tmpMin = max(dmean2 - fixRangeSDs * sdSec, dmin2 - rangeAdd)
              tmpMax = min(dmean2 + fixRangeSDs * sdSec, dmax2 + rangeAdd)
              if (sdSec < sdCritForScaleFix) needScaleFix = .true.
              rangeBase = 0.
              if (modeOfFirst == 1) rangeBase = 32768.
              !
              ! Range fixing is needed if min goes below base; only allowed if max is OK
              if (tmpMin < -rangeBase) needRangeFix = .true.
              if (tmpMax + rangeBase > 32767.) rangeFixOK = .false.
              if ((tmpMax + rangeBase) * fixScaleFac > 32767.) scaleRangeFixOK = .false.
              !
              ! scaling alone is done and thus evaluated differently for mode 1 and 6
              if (dmean2 < 0.) then
                if ((tmpMax + rangeBase) * fixScaleFac > 32767.)  &
                    scaleFixShiftedOK = .false.
              else
                if (tmpMax * fixScaleFac > 65535. - rangeBase) scaleFixAloneOK = .false.
              endif
              if (iVerbose > 0) print *,dmin2, dmax2, dmean2, sdSec, tmpMin, tmpMax, &
                  rangeBase, needScaleFix, needRangeFix, rangeFixOK, scaleRangeFixOK, &
                  scaleFixAloneOK, scaleFixShiftedOK
              cycle
            endif
            !
            ! Otherwise store the values
            secMean(ind) = dmean2
            secMins(ind) = dmin2
            secMaxes(ind) = dmax2
            !
            if (ifFloat == 2) then
              !
              ! find the min and max Z values ((density-mean) /sd)
              !
              secZmins(ind) = 0.
              secZmaxes(ind) = 0.
              if (dmax2 > dmin2 .and. sdSec > 0.) then
                secZmins(ind) = (dmin2 - dmean2) / sdSec
                secZmaxes(ind) = (dmax2 - dmean2) / sdSec
              endif
            else
              !
              ! or, if shifting, get maximum range from mean
              !
              diffMinMean = min(diffMinMean, dmin2 - dmean2)
              diffMaxMean = max(diffMaxMean, dmax2 - dmean2)
              grandSum = grandSum + dmean2
              nsum = nsum + 1
            endif
          endif
        enddo     ! end of loop on files
        call iiuClose(1)
        if (needClose1 > 0) call iiuClose(needClose1)
      enddo     ! End of loop on sections
      !
      ! for shift to mean, figure out new mean, min and max and whether
      ! scaling will be needed to fit range
      !
      if (ifFloat > 2) then
        grandMean = grandSum / nsum
        shiftMin = max(0., grandMean + diffMinMean)
        shiftMean = shiftMin - diffMinMean
        shiftMax = shiftMean + diffMaxMean
        if (ifFloat == 3 .and. mode .ne. 2 .and. newMode .ne. 2 .and. &
            optimalMax(mode + 1) < shiftMax) then
          print *,'Densities will be compressed by', &
              optimalMax(mode + 1) / shiftMax, ' to fit in range'
          floatText = ', mean shift&scaled'
        endif
      endif
      !
      ! For float to mean, find outliers and get zmin and zmax without them
      if (ifFloat == 2) then
        !call rsMedian(secZmaxes, numOutTot, ztemp, zmaxMed)
        !call rsMADN(secZmaxes, numOutTot, zmaxMed, ztemp, zmaxMADN)
        !print *,'median', zmaxMed, '   MADN', zmaxMADN
        !write(*,'(8f9.2)') (ztemp(i) / zmaxMADN, i = 1, numOutTot)
        call rsMadMedianOutliers(secZmins, numOutTot, 8., ztemp)
        do i = 1, numOutTot
          if (ztemp(i) >= 0.) then
            zmin = min(zmin, secZmins(i))
          else
            numSecTrunc = numSecTrunc + 1
          endif
        enddo
        call rsMadMedianOutliers(secZmaxes, numOutTot, 8., ztemp)
        do i = 1, numOutTot
          if (ztemp(i) <= 0.) then
            zmax = max(zmax, secZmaxes(i))
          else
            numSecTrunc = numSecTrunc + 1
          endif
        enddo
        deallocate(secZmins, secZmaxes, ztemp)
      endif
      !
      ! Finalize the range scaling/fixing actions, set ifFloat -1 for any kind
      ! Set numScaleFacs, scaleFacs(1), and scaleConsts(1) for a single scaling, set
      ! newMode if switching to 1 from 6
      if (fixRangeSDs > 0.) then
        doRangeScale = needRangeFix .and. needScaleFix .and. rangeFixOK .and. &
            scaleRangeFixOK
        doRangeOnly = .not. doRangeScale .and. needRangeFix .and. rangeFixOK
        doScaleOnly = .not. doRangeScale .and. needScaleFix .and.  &
            (scaleFixAloneOK .or. scaleFixShiftedOK)
        if (doRangeScale .or. doRangeOnly .or. doScaleOnly) then
          numScaleFacs = 1
          scaleFacs(1) = 1.
          if (doRangeScale .or. doScaleOnly) then
            scaleFacs(1) = fixScaleFac
            write(*,'(a,f6.1,a,/,a,f6.1,/)')'INFO: Newstack scaling values by ', &
                fixScaleFac, ' to preserve intensity resolution', &
                '  because SD of values is below', sdCritForScaleFix
          endif
          ifFloat = -1
          if (doRangeScale .or. doRangeOnly) then
            if (newMode == 6) then
              newMode = 1
              scaleConsts(1) = 0.
              write(*,'(a,/)')'INFO: Newstack changing mode to 1 (signed integer) to '// &
                  'avoid truncation'
            else
              scaleConsts(1) = 32768. * scaleFacs(1)
              write(*,'(a,a,/,a,/)')'INFO: Newstack shifting values up by 32768 to ', &
                  'avoid truncation',  &
                  '  If taking logarithm in Tilt program, make the offset 0 instead 32768.'
            endif
          else if (newMode == 1 .and. scaleFixShiftedOK) then
            scaleConsts(1) = 32768. * (scaleFacs(1) - 1.)
          endif
        endif
      endif
    endif
  endif
  !
  ! start looping over input images for processing and output
  !
  isec = 1
  isecOut = 1
  isecReplace = 1
  iOutFile = 1
  if (.not. quiet) call iiuAltPrint(1)
  call time(timeStr)
  call b3dDate(dat)
  numTruncLow = 0
  numTruncHigh = 0
  ifHeaderOut = 0
  ifTempOpen = 0
  do iFile = 1, numInFiles
    call openInputFile(iFile)
    call irdhdr(1, nxyz, mxyz, mode, dminIn, dmaxIn, dmeanIn)
    call iiuRetSize(1, nxyz, mxyz, nxyzst)
    call iiuRetCell(1, cell)
    call iiuFileInfo(1, ix1, ix2, iiuFlags)
    packed4bitInput = btest(iiuFlags, 5) .or. btest(iiuFlags, 6)
    !
    ! get the binned size to read
    !
    call getReducedSize(nx, readReduction, readShrunk, nxBin, rxOffset)
    call getReducedSize(ny, readReduction, readShrunk, nyBin, ryOffset)
    if (iVerbose > 0) &
        print *,'Size and offsets X:', nxBin, rxOffset, ', Y:', nyBin, ryOffset
    !
    ! get extra header information if any
    !
    call iiuRetNumExtended(1, nByteSymIn)
    itype = 0
    numIntOrBytesIn = 0
    iFlagExtraIn = 0
    FEI1type = .false.
    if (needHeaderAngles .and. nByteSymIn == 0) call exitError('There is no extended'// &
        ' header; tilt angles for -reorder must be entered with the -angles option')
    if (nByteSymIn > 0) then
      !
      ! Deallocate array if it was allocated and is not big enough
      if (maxExtraIn > 0 .and. nByteSymIn > maxExtraIn) then
        deallocate(extraIn, stat = ierr)
        maxExtraIn = 0
      endif
      !
      ! Allocate array if needed
      if (maxExtraIn == 0) then
        maxExtraIn = nByteSymIn + 1024
        allocate(extraIn(maxExtraIn), stat = ierr)
        if (ierr .ne. 0) call exitError('Allocating memory for extra header arrays')
      endif

      call iiuRetExtendedData(1, nByteSymIn, extraIn)
      call iiuRetExtendedType(1, numIntOrBytesIn, iFlagExtraIn)
      !
      ! DNM 4/18/02: if these numbers do not represent bytes and
      ! flags, then number of bytes is 4 times nint + nreal
      !
      serialEMtype = nbytes_and_flags(numIntOrBytesIn, iFlagExtraIn)
      itype = 1
      if (serialEMtype) itype = -1
      !
      ! Mark new FEI type as 2 and no type for unknown
      if (numIntOrBytesIn < 0) then
        if (numIntOrBytesIn == -3) then
          itype = 2
          FEI1type = .true.
        else
          itype = 0
          nByteSymIn = 0
        endif
      endif
      if (serialEMtype .and. saveTilts .and. mod(iFlagExtraIn, 2) == 0)  &
          call exitError('You cannot save tilt angles into a SerialEM extended '// &
          'header that was not saved with tilt angles')
      !
      ! Get tilt angles for reordering
      if (needHeaderAngles) then
        ! 
        ! Take care of tilt angle array 
        if (nxyz(3) > maxInFileAngles) then
          if (maxInFileAngles > 0) deallocate(allFileTilts, izNotPiece, stat = ierr)
          maxInFileAngles = nxyz(3) + 8
          allocate(allFileTilts(maxInFileAngles), izNotPiece(maxInFileAngles),  &
              stat = ierr)
          call memoryError(ierr, 'arrays for all tilt angles in file')
        endif
        !
        ! Set up dummy piece lists and get the tilt angles
        do ind = 1, nxyz(3)
          izNotPiece(ind) = ind - 1
        enddo
        call get_extra_header_tilts(extraIn, nByteSymIn, numIntOrBytesIn, iFlagExtraIn, &
            nxyz(3), allFileTilts, ind, maxInFileAngles, izNotPiece)
        if (ind < nxyz(3)) call exitError('There are either not enough or no tilt '// &
            'angles in extended header; tilt angles for -reorder must be entered '// &
            'with the -angles option')
        !
        ! Fill the angle arrays from the section list
        do ind = listInd(iFile), listInd(iFile) + nlist(iFile) - 1
          iSecRead = inList(ind)
          extraTilts(ind) = allFileTilts(iSecRead + 1)
        enddo
      endif
    endif

    if (ianyExtraType == 0) ianyExtraType = itype
    if (ifile == 1) ifirstExtraType = itype
    if (ianyExtraType .ne. 0 .and. itype .ne. 0 .and. ianyExtraType .ne. itype .and.  &
        .not. stripExtra) call exitError('You cannot include files with different '// &
        'types of extended headers in the same run; add the -strip option')
    if (saveTilts .and. itype .ne. ifirstExtraType) call exitError( &
        'To store tilt angles, all input files must have the same type of extended'// &
        ' header or no extended header')
    if (saveTilts .and. itype .eq. 2) call exitError( &
        'You cannot store tilt angles back into a new FEI1-style extended header')
    ierr = 0
    if (useMdocFiles) ierr = 1
    indAdocIn = iiuRetAdocIndex(1, 0, ierr)
    !
    ! Reorder the section list and output angles for this file
    if (ifReorderByTilt .ne. 0) then
       do ind = listInd(iFile), listInd(iFile) + nlist(iFile) - 2
         do ilist = ind + 1,  listInd(iFile) + nlist(iFile) - 1
           if (sign(1, ifReorderByTilt) * (extraTilts(ind) - extraTilts(ilist)) &
               > 0.01) then
             ix1 = inList(ind)
             inList(ind) = inList(ilist)
             inList(ilist) = ix1
             tmpMax = extraTilts(ind)
             extraTilts(ind) = extraTilts(ilist)
             extraTilts(ilist) = tmpMax
           endif
         enddo
       enddo
    endif
    !
    ! get each section in input file
    !
    do ilist = 1, nlist(iFile)
      iSecRead = inList(ilist + listInd(iFile) - 1)
      !
      ! set output characteristics from first section, transposing size
      ! for 90 degree rotation and resizing as necessary
      !
      if (isec == 1) then
        if (abs(abs(rotateAngle) - 90.) < 0.1 .and. nxOut <= 0 .and. &
            nyOut <= 0) then
          nxOut = nint(nyBin * expandFactor)
          nyOut = nint(nxBin * expandFactor)
        endif
        if (nxOut <= 0) nxOut = nint(nxBin * expandFactor)
        if (nyOut <= 0) nyOut = nint(nyBin * expandFactor)
        !
        ! if warping or distortions, figure out how big to allocate the arrays
        if (ifDistort + ifWarping .ne. 0) then
          allocate(nControl(numFields), stat = ierr)
          call memoryError(ierr, 'Array for number of control points')
          dx = 0.
          dy = 0.
          !
          ! Expanded grid size is based on the input size for distortion and the
          ! output size for warping
          if (ifDistort .ne. 0) then
            xnBig = nxMax / warpScale
            ynbig = nyMax / warpScale
          else
            xnBig = readReduction * nxOut / warpScale
            ynbig = readReduction * nyOut / warpScale
            if (applyFirst == 0) then
              dx = 1.e20
              dy = 1.e20
              xnBig = -dx
              xnBig = -dx
              do i = 1, listTotal
                xnBig = max(xnBig, (nxOut + xcen(i)) * readReduction / warpScale)
                dx = min(dx, xcen(i) * readReduction / warpScale)
                ynbig = max(ynbig, (nyOut + ycen(i)) * readReduction / warpScale)
                dy = min(dy, ycen(i) * readReduction / warpScale)
              enddo
            endif
          endif
          if (findMaxGridSize(dx, xnBig, dy, ynbig, nControl, maxNxGrid, maxNyGrid, &
              listString) .ne. 0) call exitError(listString)
          ! print *,maxNxGrid, maxNyGrid
          lmGrid = max(lmGrid, maxNxGrid, maxNyGrid)
        endif
        if (ifDistort + ifMagGrad + ifWarping .ne. 0) then
          allocate(fieldDx(lmGrid, lmGrid), fieldDy(lmGrid, lmGrid), &
              tmpDx(lmGrid, lmGrid), tmpDy(lmGrid, lmGrid), stat = ierr)
          call memoryError(ierr, 'Arrays for warping fields')
        endif
      endif
      !
      if (numTaper == 1) then
        numTaper = min(127, max(16, nint((nxOut + nyOut) / 200.)))
        write(*,'(/,a,i4,a)') 'Tapering will be done over', numTaper, ' pixels'
      endif
      !
      ! First see if this is the first section to replace
      if (numReplace > 0 .and. isecReplace == 1) then
        if (nxOut .ne. nxyz2(1) .or. nyOut .ne. nxyz2(2)) call exitError( &
            'Existing output file does not have right size in X or Y')
        if (newMode .ne. modeOld) call exitError( &
            'Output mode does not match existing output file')
        isecOut = listReplace(isecReplace) + 1
      endif
      !
      ! then see if need to open an output file
      if (numReplace == 0 .and. isecOut == 1) then
        !
        ! Create output file, transfer header from currently open file, fix it
        !
        if (if3dVolumes > 1) then
          call iiAllowMultiVolume(1)
          call imopen(12, outFile(iOutFile), 'OLD')
          if (iiuVolumeOpen(2, 12, -1) .ne. 0) call exitError( &
              'Opening new volume in existing file')
          needClose2 = 12
          call iiAllowMultiVolume(0)
          if (if3dVolumes == 3) then
            indGlobalAdoc = iiuRetAdocIndex(12, 1, 0)
            if (indGlobalAdoc >= 0) then
              if (AdocSetCurrent(indGlobalAdoc) .ne. 0) then
                indGlobalAdoc = -1
              else if (AdocSetInteger(globalName, 1, 'image_pyramid', 1) .ne. 0) then
                indGlobalAdoc = -1
              else if (iiuWriteGlobalAdoc(12) .ne. 0) then
                indGlobalAdoc = -1
              endif
            endif
            if (indGlobalAdoc < 0) call exitError( &
                'Setting image_pyramid attribute in global section of HDF file')
          endif
        else
          call setNextOutputSize(nxOut, nyOut, numSecOut(iOutFile), newMode);
          call imopen(2, outFile(iOutFile), 'NEW')
          needClose2 = 0
        endif
        if (if3dVolumes > 0) then
          nxTileIn = nxTile
          nyTileIn = nyTile
          nzChunkIn = nzChunk
          call iiBestTileSize(nxOut, nxTileIn, ierr, 1)
          call iiBestTileSize(nyOut, nyTileIn, ierr, 1)
          call iiBestTileSize(numSecOut(iOutFile), nzChunkIn, ierr, 1)
          if (iiuAltChunkSizes(2, nxTileIn, nyTileIn, nzChunkIn) .ne. 0) call exitError( &
              'Setting chunk sizes in new volume')
          if (nxTileIn .ne. nxOut .or. nyTileIn .ne. nyOut .or. nzChunkIn .ne.  &
              numSecOut(iOutFile)) write(*,'(a,i7,a,i7,a,i4)') 'Actual chunk size: ', &
              nxTileIn, ' by', nyTileIn, ' by', nzChunkIn
        endif

        call iiuTransHeader(2, 1)
        call iiuAltMode(2, newMode)
        !
        ! set new size, keep old nxyzst
        !
        nxyz2(1) = nxOut
        nxyz2(2) = nyOut
        nxyz2(3) = numSecOut(iOutFile)
        call iiuAltSize(2, nxyz2, nxyzst)
        !
        ! if mxyz=nxyz, keep this relationship
        !
        if (mxyz(1) == nx .and. mxyz(2) == ny .and. mxyz(3) == nz) then
          mxyz2(1) = nxOut
          mxyz2(2) = nyOut
          mxyz2(3) = numSecOut(iOutFile)
          call iiuAltSample(2, mxyz2)
        else
          mxyz2(1) = mxyz(1)
          mxyz2(2) = mxyz(2)
          mxyz2(3) = mxyz(3)
        endif
        !
        ! keep delta the same by scaling cell size from change in mxyz
        !
        do i = 1, 3
          cell2(i) = mxyz2(i) * (cell(i) / mxyz(i)) * readReduction / expandFactor
          cell2(i + 3) = 90.
        enddo
        cell2(3) = mxyz2(3) * (cell(3) / mxyz(3))
        call iiuAltCell(2, cell2)
        !
        ! shift origin by the fraction pixel offset when binning or reducing with read
        ! a positive change is needed to indicate origin inside image
        ! When reduction was added, removed a convoluted subpixel adjustment to this
        ! that was wrong as basis for adjusting origin for size changes
        call iiuRetDelta(1, delta)
        call iiuRetOrigin(1, xOrigin, yOrigin, zOrigin)
        if (readReduction > 1) then
          xOrigin = xOrigin - delta(1) * rxOffset
          yOrigin = yOrigin - delta(2) * ryOffset
        endif
        !
        ! Adjust origin if requested: it is different depending on whether
        ! there are transforms and whether offset was applied before or
        ! after.  delta can be modified, it will be reread
        if (adjustOrigin) then
          zOrigin = zOrigin - iSecRead * delta(3)
          delta(1) = delta(1) * readReduction / expandFactor
          delta(2) = delta(2) * readReduction / expandFactor
          if (ifXform == 0 .and. ftReduceFac <= 0. .and. .not. phaseShift) then
            xOrigin = xOrigin - (nxBin / 2 + nint(xcen(isec)) - nxOut / 2) * delta(1)
            yOrigin = yOrigin - (nyBin / 2 + nint(ycen(isec)) - nyOut / 2) * delta(2)
          elseif (applyFirst .ne. 0) then
            xOrigin = xOrigin - (expandFactor * (nxBin / 2. + xcen(isec))  - nxOut / 2.) &
                * delta(1)
            yOrigin = yOrigin - (expandFactor * (nyBin / 2. + ycen(isec)) - nyOut / 2.)  &
                * delta(2)
          else
            xOrigin = xOrigin - (expandFactor * nxBin / 2. + xcen(isec) - nxOut / 2.) *  &
                delta(1)
            yOrigin = yOrigin - (expandFactor * nyBin / 2. + ycen(isec) - nyOut / 2.) *  &
                delta(2)
          endif
        endif
        if (adjustOrigin .or. readReduction > 1) &
            call iiuAltOrigin(2, xOrigin, yOrigin, zOrigin)
        !
        if (trunctText == ' ') then
          write(titlech, 302) xfText, floatText, dat, timeStr
        else
          write(titlech, 301) xfText, floatText, trunctText
        endif
        read(titlech, '(20a4)') (title(kti), kti = 1, 20)
301     format('NEWSTACK: Images copied',a13,a18,a20)
302     format('NEWSTACK: Images copied',a13,a18,t57,a9,2x,a8)
        dmax = -1.e30
        dmin = 1.e30
        dmean = 0.
        !
        ! adjust extra header information if currently open file has it
        !
        nByteSymOut = 0
        outDocChanged = .false.
        if ((nByteSymIn > 0 .or. saveTilts) .and. .not.stripExtra .and. &
            b3dOutputFileType() == 2) then
          if (nByteSymIn == 0) then
            nByteExtraOut = 4
            serialEMtype = .false.
            call iiuAltExtendedType(2, 0, 1)
          else
            ierr = getExtraHeaderMaxSecSize(extraIn, nByteSymIn, numIntOrBytesIn, &
                iFlagExtraIn, nxyz(3), nByteExtraOut)
          endif
          nByteSymOut = numSecOut(iOutFile) * nByteExtraOut
          if (maxExtraOut > 0 .and. nByteSymOut > maxExtraOut) then
            deallocate(extraOut, stat = ierr)
            maxExtraOut = 0
          endif
          !
          ! Allocate array if needed
          if (maxExtraOut == 0) then
            maxExtraOut = nByteSymOut + 1024
            allocate(extraOut(maxExtraOut), stat = ierr)
            call memoryError(ierr, 'Arrays for extra header data')
          endif
          call iiuAltNumExtended(2, nByteSymOut)
          call iiuSetPosition(2, 0, 0)
          indExtraOut = 0
        else
          call iiuAltNumExtended(2, 0)
        endif

        ierr = 0
        if (useMdocFiles) ierr = -1
        indAdocOut = iiuRetAdocIndex(2, 0, ierr)
        !
        ! Transfer global section if either file is not HDF; itrhdr takes care of HDF->HDF
        if ((iiuFileType(2) .ne. 5 .or. iiuFileType(1) .ne. 5) .and. indAdocOut > 0  &
            .and. indAdocIn > 0) then
          call setCurrentAdocOrExit(indAdocIn, 'input')
          if (AdocTransferSection(globalName, 1, indAdocOut, globalName, 0) .ne. 0) &
              call exitError('Transferring global data between autodocs')
        endif
      endif
      !
      ! handle complex images here and skip out
      !
      if ((newMode + 1) / 2 == 2 .or. (mode + 1) / 2 == 2) then
        if ((mode + 1) / 2 .ne. 2 .or. (newMode + 1) / 2 .ne. 2) call exitError( &
            'All input files must be complex if any are')

        if (limEntered == 0 .and. nx * ny * 2 > idimInOut) then
          idimInOut = nx * ny * 2
          lenTemp = 1
          limToAlloc = idimInOut + 1
          call reallocateArray()
        endif
        if (nx * ny * 2 > idimInOut) call exitError('Input image too large for array.')

        call iiuSetPosition(1, iSecRead, 0)
        call irdsec(1, array,*99)
        call iclcdn(array, nx, ny, 1, nx, 1, ny, dmin2, dmax2, dmean2)
        call iiuSetPosition(2, isecOut - 1, 0)
        call iiuWriteSection(2, array)
        go to 80
      endif
      !
      ! determine whether rescaling will be needed
      !
      ! for no float: if either mode = 2, no rescale;
      ! otherwise rescale from input range to output range only if
      ! mode is changing
      !
      rescale = .false.
      if (ifFloat == 0 .and. newMode .ne. 2 .and. mode .ne. 2) then
        rescale = mode .ne. newMode .or. (pack4bitOutput .and. .not.packed4bitInput)
      elseif (ifFloat .ne. 0) then
        rescale = .true.
      endif
      preSetScaling = ifFloat <= 0
      !
      optimalIn = optimalMax(mode + 1)
      optimalOut = optimalMax(newMode + 1)
      if (packed4bitInput) optimalIn = 15.
      if (pack4bitOutput) optimalOut = 15.
      !
      ! set bottom of input range to 0 unless mode 1 or 2 and already negative or not
      ! rescaling; set bottom of output range to 0 unless not changing modes
      !
      bottomIn = 0.
      if (dminIn < 0. .or. .not. rescale) then
        if (mode == 1) bottomIn = -optimalIn - 1
        if (mode == 2) bottomIn = -optimalIn
      endif
      bottomOut = 0.
      if (mode == newMode) bottomOut = bottomIn
      !
      ! Handle blank images here and skip out
      !
      if (iSecRead < 0 .or. iSecRead >= nz) then
        tmpMin = dmeanIn
        if (ifUseFill .ne. 0) tmpMin = fillVal
        tmpMax = tmpMin
        dsumSq = 0.
        dsum = tmpMin * (float(nxOut) * nyOut)
        nyNeeded = 1
        call reallocateIfNeeded()
        call findScaleFactors(tmpMin, tmpMax)
        do i = 1, nxOut
          array(i) = tmpMin * scaleFactor + constAdd
        enddo
        call iiuSetPosition(2, isecOut - 1, 0)
        do i = 1, nyOut
          call iwrlin(2, array)
        enddo
        dmin2 = array(1)
        dmax2 = dmin2
        dmean2 = dmin2
        go to 80
      endif
      if (iVerbose > 0) print *,'rescale', rescale
      !
      ! Get the index of the transform
      if (ifXform .ne. 0) then
        lnu = lineUse(isec) + 1
        call xfcopy(f(1, 1, lnu), fprod)
      endif
      !
      ! if doing distortions or warping, get the grid
      !
      hasWarp = .false.
      dx = 0.
      dy = 0.
      if (ifDistort > 0) then
        iy = idfUse(isec) + 1
        hasWarp = .true.
        xnBig = nx / warpScale
        ynbig = ny / warpScale
      else if (ifWarping .ne. 0) then
        iy = lnu
        hasWarp = nControl(iy) > 2
        xnBig = readReduction * nxOut / warpScale
        ynbig = readReduction * nyOut / warpScale
        !
        ! for warping with center offset applied after, it will subtract the
        ! offset from the grid start and add it to the grid displacements
        if (applyFirst == 0) then
          dx = readReduction * xcen(isec) / warpScale
          dy = readReduction * ycen(isec) / warpScale
        endif
      endif
      if (hasWarp) then
        if (getSizeAdjustedGrid(iy, xnBig, ynbig, dx, dy, 1, warpScale,  &
            nint(readReduction), nxGrid, nyGrid, xGridStart, yGridStart, xGridIntrv, &
            yGridIntrv, fieldDx, fieldDy, lmGrid, lmGrid, listString) .ne. 0) &
            call exitError(listString)
        !
        !print *,nxGrid, nyGrid, xGridIntrv, yGridIntrv, xGridStart, xGridStart + (nxGrid &
        ! - 1) * xGridIntrv, yGridStart, yGridStart + (nyGrid - 1)*yGridIntrv
        ! do ierr = 1, nyGrid
        ! write(*,'(f7.2,9f8.2)') (fieldDx(i, ierr), fieldDy(i, ierr), i=1, nxGrid)
        ! enddo
        !
        ! copy field to tmpDx, y in case there are mag grads
        tmpDx(1:nxGrid, 1:nyGrid) = fieldDx(1:nxGrid, 1:nyGrid)
        tmpDy(1:nxGrid, 1:nyGrid) = fieldDy(1:nxGrid, 1:nyGrid)
      endif
      !
      ! if doing mag gradients, set up or add to distortion field
      !
      if (ifMagGrad .ne. 0) then
        magUse = min(iSecRead + 1, numMagGrad)
        if (ifDistort .ne. 0) then
          call addMagGradField(tmpDx, tmpDy, fieldDx, fieldDy, lmGrid, nxBin, &
              nyBin, nxGrid, nyGrid, xGridStart, yGridStart, xGridIntrv, &
              yGridIntrv, nxBin / 2., nyBin / 2., pixelMagGrad, axisRot, &
              tiltAngles(magUse), dmagPerMicron(magUse), rotPerMicron(magUse))
        else
          call makeMagGradField(tmpDx, tmpDy, fieldDx, fieldDy, lmGrid, &
              nxBin, nyBin, nxGrid, nyGrid, xGridStart, yGridStart, xGridIntrv, &
              yGridIntrv, nxBin / 2., nyBin / 2., &
              pixelMagGrad, axisRot, tiltAngles(magUse), &
              dmagPerMicron(magUse), rotPerMicron(magUse))
        endif
      endif
      !
      ! if transforming, and apply first is selected, get the shifts by
      ! applying the offset first then multiplying that by the transform
      ! Otherwise do it the old way, subtract offset from final xform, but only
      ! if not warping (now that warping is known for this section)
      if (ifXform .ne. 0) then
        if (applyFirst .ne. 0) then
          call xfunit(frot, 1.)
          frot(1, 3) = -xcen(isec)
          frot(2, 3) = -ycen(isec)
          call xfmult(frot, f(1, 1, lnu), fprod)
        elseif (.not. (ifWarping .ne. 0 .and. hasWarp)) then
          fprod(1, 3) = fprod(1, 3) - xcen(isec)
          fprod(2, 3) = fprod(2, 3) - ycen(isec)
        endif
      endif
      !
      ! get maximum Y deviation with current field to adjust chunk
      ! limits with
      !
      if (ifMagGrad .ne. 0 .or. hasWarp) then
        fieldMaxX = 0.
        fieldMaxY = 0.
        do iy = 1, nyGrid
          do i = 1, nxGrid
            fieldMaxX = max(fieldMaxX, abs(fieldDx(i, iy)))
            fieldMaxY = max(fieldMaxY, abs(fieldDy(i, iy)))
          enddo
        enddo
        maxFieldX = int(fieldMaxX + 1.5)
        maxFieldY = int(fieldMaxY + 1.5)
      endif
      !
      ! Determine starting and ending lines needed from the input, and whether any
      ! fill is needed
      call linesNeededForOutput(0, nyOut - 1, needYfirst, needYlast, fillNeeded, inPlace)
      nyNeeded = needYlast + 1 - needYfirst
      !
      ! Get padded input size for phase shifting
      if (phaseShift) then
        nxFSpad = niceFrame(min(maxFSpad, max(minFSpad, nint(fsPadFrac * nxBin))) + &
            nxBin, 2, niceFFTlimit())
        nyFSpad = niceFrame(min(maxFSpad, max(minFSpad, nint(fsPadFrac * nyNeeded))) + &
            nyNeeded, 2, niceFFTlimit())
        nxFCropPad = nxFSpad
        nyFCropPad = nyFSpad
      endif
      nxDimNeed = max(nxBin, nxFSpad + 2)
      nyDimNeed = max(nyNeeded, nyFSpad + 1)
      !
      ! It is possible to process in place instead of having a separate output buffer
      ! if there are no transformations, the output fits into the input space, the
      ! linesNeeded function indicates that the shifts in xcen/ycen and the output size
      ! do not result in padding to the left or below the image.  It can work with any
      ! scaling if the whole input image fits, so set it this way first
      ! Fourier shifting alone IS done in place with whatever portion of input data is
      ! needed to provide for a fractional shift of actual data, so an output buffer
      ! is needed only if padding below/to left
      ! Fourier reduction is done from the whole input image into the temporary space
      ! and that can be processed back into the input buffer as long as the output isn't
      ! really padded bigger than the input
      ! 
      processInPlace = ifXform == 0 .and. nxOut <= nxDimNeed .and. nyOut <= nyDimNeed &
          .and. inPlace
      !
      ! Now that needed input and output size is finally known, make sure memory is
      ! enough.  This routine will turn off processInPlace before trying to allocate
      ! if the input image doesn't fit within the memory limit, but if it doesn't get
      ! what it asks for and has to fall back, we need to test here if the processInPlace
      ! should be turned off because of no preSetScaling
      call reallocateIfNeeded()
      if (idimInOut / nxDimNeed <= nyDimNeed .and. .not. preSetScaling) &
          processInPlace = .false.
      if (iVerbose > 0) print *,'preSetScaling ', preSetScaling, '   processInPlace ',  &
          processInPlace
      !
      ! figure out how the data will be loaded and saved; first see if both input and
      ! output images will fit in one chunk, or if entire input image will fit
      ! Leave numChunks as 0 if input image will not fit
      ! If it does fit, linesLeft is the number of output lines that leave space for
      ! and numChunks is the number of output chunks needed
      ! If processing in place, linesLeft is the amount of input it can hold
      numChunks = 0
      if (idimInOut / nxDimNeed > nyDimNeed) then
        linesLeft = (idimInOut - nxDimNeed * int(nyDimNeed, kind = 8)) / nxOut
        if (processInPlace) linesLeft = idimInOut / nxDimNeed
        numChunks = (nyOut + linesLeft - 1) / linesLeft
        if (iVerbose > 0) print *,'linesleft', linesLeft, '  nchunk', numChunks
      endif
      if (numChunks > 1 .and. numTaper > 0) &
          call exitError('Cannot taper output image - it does not '// &
          'fit completely in memory')
      if (numChunks > 1 .and. (phaseShift .or. ftReduceFac > 0.)) &
          call exitError('Cannot apply Fourier shift - input and output images do not '//&
          'fit completely in memory')
      if (numChunks > maxChunks .and. processInPlace) call exitError( &
          'The images are too large to process with the current memory limit')

      if (numChunks == 1 .or. (numChunks > 0 .and. numChunks <= maxChunks .and. &
          preSetScaling)) then
        !
        ! Use entire input and multi-chunk output only if not rescaling in way that 
        ! requires min/max/mean of output.  This encompasses case of processing in place,
        ! where only the identical chunk of input is needed for output
        ! set up chunks for output, and set up that
        ! whole image is needed for every output chunk.
        ! Note that lines are numbered from 0 here
        !
        lineOutSt(1) = 0
        do iChunk = 1, numChunks
          nextLine = (nyOut / numChunks) * iChunk + min(iChunk, mod(nyOut, numChunks))
          numLinesOut(iChunk) = nextLine - lineOutSt(iChunk)
          lineOutSt(iChunk + 1) = nextLine
          if (processInPlace .and. ftReduceFac == 0 .and. .not. phaseShift) then
            lineInSt(iChunk) = needYfirst + lineOutSt(iChunk)
            numLinesIn(iChunk) = numLinesOut(iChunk)
          else
            lineInSt(iChunk) = needYfirst
            numLinesIn(iChunk) = nyNeeded
          endif
        enddo
        maxin = nyDimNeed
        ifOutChunk = 1
      else
        !
        ! otherwise, break output into successively more chunks, and
        ! find out how many lines are needed of input for each chunk,
        ! Scan twice: first trying to fit all the output to avoid a
        ! read-write to temp file; then breaking equally
        !
        ifOutChunk = -1
        iscan = 1
        iyTest = nyOut
        do while(iscan <= 2 .and. ifOutChunk < 0)
          numChunks = 1
          do while(numChunks <= maxChunks .and. ifOutChunk < 0)
            lineOutSt(1) = 0
            maxin = 0
            do iChunk = 1, numChunks
              nextLine = (nyOut / numChunks) * iChunk + min(iChunk, mod(nyOut, numChunks))
              numLinesOut(iChunk) = nextLine - lineOutSt(iChunk)
              lineOutSt(iChunk + 1) = nextLine
              !
              call linesNeededForOutput(lineOutSt(iChunk), nextLine - 1, iy1, iy2, &
                  fillTmp, inPlace)
              lineInSt(iChunk) = iy1
              numLinesIn(iChunk) = iy2 + 1 - iy1
              maxin = max(maxin, numLinesIn(iChunk))
            enddo
            !
            ! Will the input and output data now fit?  then terminate.
            !
            if (iscan == 2) iyTest = numLinesOut(1)
            if (idimInOut / maxin > nxBin .and. idimInOut / iyTest > nxOut .and. &
                 maxin * int(nxBin, kind = 8) + iyTest * int(nxOut, kind = 8) <= &
                 idimInOut) then
              ifOutChunk = iscan - 1
            else
              numChunks = numChunks + 1
            endif
          enddo
          iscan = iscan + 1
        enddo
        if (ifOutChunk < 0) call exitError(' Input image too large for array.')
      endif
      if (iVerbose > 0) then
        print *,'number of chunks:', numChunks, ifOutChunk
        do i = 1, numChunks
          print *,i, lineInSt(i), numLinesIn(i), lineOutSt(i), numLinesOut(i)
        enddo
      endif
      !
      ! open temp file if one is needed
      !
      if (rescale .and. .not. preSetScaling .and. ifOutChunk > 0 .and. numChunks > 1  &
          .and. ifTempOpen == 0) then
        tempExt = 'nws      '
        tempExt(4:5) = timeStr(1:2)
        tempExt(6:7) = timeStr(4:5)
        tempExt(8:9) = timeStr(7:8)
        tempName = temp_filename(outFile(iOutFile), ' ', tempExt)
        !
        call imopen(3, tempName, 'scratch')
        call iiuCreateHeader(3, nxyz2, nxyz2, 2, title, 0)
        ifTempOpen = 1
      endif
      !
      iBufOutBase = int(maxin, kind = 8) * nxDimNeed + 1
      !
      ! get the mean of section from previous scan, or a new scan
      !
      loadYstart = -1
      loadYend = -1
      if (ifUseFill .ne. 0) then
        dmeanSec = fillVal
      elseif (ifMean .ne. 0 .and. ifMeanSdEntered == 0) then
        dmeanSec = secMean(ilist + listInd(iFile) - 1)
      elseif (.not. fillNeeded) then
        dmeanSec = dmeanIn
      else
        if (iVerbose > 0) print *,'scanning for mean for fill', iSecRead, nyNeeded, &
            needYfirst,lenTemp
        wallStart = wallTime()
        call scanSection(array, idimInOut, nxBin, nyNeeded, needYfirst, readReduction, &
            rxOffset, ryOffset, iSecRead, 0, dmin2, dmax2, dmeanSec, sdSec, loadYstart, &
            loadYend, indFilter, readShrunk, fixRangeSDs, array(idimInOut + 1), lenTemp)
        loadYend = min(loadYend, loadYstart + maxin - 1)
        loadTime = loadTime + wallTime() - wallStart
      endif
      !
      ! loop on chunks
      ! tmpMin and max accumulate values for processed chunks before scaling,
      ! dmin2/dmax2/dmean2 accumulate values for final output and chunks are written
      dsum = 0.
      dsumSq = 0.
      tmpMin = 1.e30
      tmpMax = -1.e30
      dmin2 = 1.e30
      dmax2 = -1.e30
      dmean2 = 0.
      do iChunk = 1, numChunks
        needYstart = lineInSt(iChunk)
        needYend = needYstart + numLinesIn(iChunk) - 1
        !
        ! first load data that is needed if not already loaded
        !
        wallStart = wallTime()
        if (needYstart < loadYstart .or. needYend > loadYend) then
          loadYoffset = needYstart
          loadBaseInd = 0
          if (loadYstart <= needYstart .and. loadYend >= needYstart) then
            !
            ! move data down if it will fill a bottom region
            !
            numMove = (loadYend + 1 - needYstart) * int(nxBin, kind = 8)
            moveOffset = (needYstart - loadYstart) * int(nxBin, kind = 8)
            if (iVerbose > 0) print *,'moving data down', numMove, moveOffset
            do i8 = 1, numMove
              array(i8) = array(i8 + moveOffset)
            enddo
            numLinesLoad = needYend - loadYend
            loadYoffset = loadYend + 1
            loadBaseInd = numMove
          elseif (needYstart <= loadYstart .and. needYend >= loadYstart) then
            !
            ! move data up if it will fill top
            !
            numMove = (needYend + 1 - loadYstart) * int(nxBin, kind = 8)
            moveOffset = (loadYstart - needYstart) * int(nxBin, kind = 8)
            if (iVerbose > 0) print *,'moving data up', numMove, moveOffset
            do i8 = numMove, 1, -1
              array(i8 + moveOffset) = array(i8)
            enddo
            numLinesLoad = loadYstart - needYstart
          else
            !
            ! otherwise just get whole needed region
            !
            numLinesLoad = needYend + 1 - needYstart
            if (iVerbose > 0) print *,'loading whole region', needYstart, needYend, &
                numLinesLoad
          endif
          call readBinnedOrReduced(1, iSecRead, array(1 + loadBaseInd), nxBin, &
              numLinesLoad, rxOffset, ryOffset + readReduction * loadYoffset, &
              readReduction, nxBin, numLinesLoad, indFilter, readShrunk, &
              array(idimInOut + 1), lenTemp)
          loadYstart = needYstart
          loadYend = needYend
        endif
        numYload = loadYend + 1 - loadYstart
        numYchunk = numLinesOut(iChunk)
        numPix = int(nxOut, kind = 8) * numYchunk
        iChunkBase = iBufOutBase
        if (ifOutChunk == 0)  &
            iChunkBase = iBufOutBase + lineOutSt(iChunk) * int(nxOut, kind = 8)
        if (processInPlace) iChunkBase = 1
        loadTime = loadTime + wallTime() - wallStart

        if (ifXform .ne. 0) then
          !
          ! do transform if called for
          !
          wallStart = wallTime()
          xcenIn = nxBin / 2.
          ycenIn = nyBin / 2. -loadYstart
          dx = fprod(1, 3)
          dy = (nyOut - numYchunk) / 2. +fprod(2, 3) - lineOutSt(iChunk)
          ! dx=f(1, 3, lnu) -xcen(isec)
          ! dy=(nyOut-numYchunk) /2.+f(2, 3, lnu) - ycen(isec) - lineOutSt(iChunk)
          if (linesShrink > 0) then
            ierr = zoomFiltInterp(array, array(iChunkBase), nxBin, numYload, nxOut, &
                numYchunk, xcenIn , ycenIn, dx, dy, dmeanSec)
            if (ierr .ne. 0) then
              write(listString, '(a,i3)') &
                  'CALLING zoomFiltInterp FOR IMAGE REDUCTION, ERROR', ierr
              call exitError(listString)
            endif
          elseif (.not. hasWarp .and. ifMagGrad == 0) then
            call cubInterp(array, array(iChunkBase), nxBin, numYload, nxOut, numYchunk, &
                fprod, xcenIn , ycenIn, dx, dy, 1., dmeanSec, ifLinear)
          else
            !
            ! if undistorting, adjust the grid start down by first loaded input line
            ! and by offset of subarea
            ! if warping, adjust it down by first output line
            ystart = yGridStart - loadYstart
            xstart = xGridStart
            if (ifWarping .ne. 0) ystart = yGridStart - lineOutSt(iChunk)
            if (ifDistort > 0) then
              ystart = ystart - warpYoffsets(isec) / readReduction
              xstart = xstart - warpXoffsets(isec) / readReduction
            endif
            call warpInterp(array, array(iChunkBase), nxBin, numYload, nxOut, numYchunk, &
                fprod, xcenIn , ycenIn, dx, dy, 1., dmeanSec, ifLinear, ifWarping, &
                fieldDx, fieldDy, lmGrid, nxGrid, nyGrid, xstart, ystart, &
                xGridIntrv, yGridIntrv)
          endif
          rotTime = rotTime + wallTime() - wallStart
        else
          !
          ! otherwise repack array into output space nxOut by nyOut, with
          ! offset as specified, using the special repack routine
          !
          ! But first apply phase shift or reduction in FFT
          ioutBase = 1
          if (phaseShift .or. ftReduceFac > 0.) then
            wallStart = wallTime()
            if (noisePad) then
              call sliceNoiseTaperPad(array, nxBin, numYload, array, nxFSpad + 2, &
                  nxFSpad, nyFSpad, max(20, min(120, max(nxBin, numYload) / 50)), 4,  &
                  array(1 + (nxFSpad + 2) * nyFSpad))
            else
              call taperOutPad(array, nxBin, numYload, array, nxFSpad + 2, nxFSpad, &
                  nyFSpad, 1, dmeanSec)
            endif
            taperTime = taperTime + wallTime() - wallStart
            call todfft(array, nxFSpad, nyFSpad, 0)
            dx = nint(xcen(isec)) - xcen(isec)
            dy = nint(ycen(isec)) - ycen(isec)
            if (phaseShift) then
              call fourierShiftImage(array, nxFSpad, nyFSpad, dx, dy,  &
                  array(1 + (nxFSpad + 2) * nyFSpad))
            else
              ioutBase = 1 + idimInOut
              call fourierReduceImage(array, nxFSpad, nyFSpad, array(ioutBase), &
                  nxFCropPad, nyFCropPad, actualFac * dx, actualFac * dy, &
                  array(1 + (nxFSpad + 2) * nyFSpad))
            endif
            call todfft(array(ioutBase), nxFCropPad, nyFCropPad, 1)
            !
            ! Replicate last real column of image into the extra two elements 
            do i = 1, nyFCropPad
              ix1 = ioutBase + i * (nxFCropPad + 2) - 3
              array(ix1 + 1) = array(ix1)
              array(ix1 + 2) = array(ix1)
            enddo
            rotTime = rotTime + wallTime() - wallStart
          endif
          !
          ! Then repack, adjusting starting coordinates for padded array
          ix1 = nxBin / 2 - nxOut / 2 + nint(xcen(isec))
          iyBase = nyBin / 2 - nyOut / 2 + nint(ycen(isec))
          iy1 = iyBase + lineOutSt(iChunk) - loadYstart
          nxRepakDim = nxDimNeed
          nyRepakDim = max(numYload, nyFSpad)
          if (phaseShift) then
            ix1 = ix1 + (nxFSpad - nxBin) / 2
            iy1 = iy1 + (nyFSpad - numYload) / 2
          elseif (ftReduceFac > 0.) then
            ix1 = nxFCropPad / 2 - nxOut / 2 + nint(xcen(isec))
            iy1 = nyFCropPad / 2 - nyOut / 2 + nint(ycen(isec))
            nxRepakDim = nxFCropPad + 2
            nyRepakDim = nyFCropPad
          endif
          ix2 = ix1 + nxOut - 1
          iy2 = iy1 + numYchunk - 1
          !
          call irepak2(array(iChunkBase), array(ioutBase), nxRepakDim, nyRepakDim, &
              ix1, ix2, iy1, iy2, dmeanSec)
          if (iVerbose > 0) print *,'did repack', iChunkBase, ioutBase, nxRepakDim, &
              nyRepakDim, ix1, ix2, iy1, iy2
        endif
        if (numTaper > 0) then
          wallStart = wallTime()
          if (taperAtFill(array(iChunkBase), nxOut, nyOut, numTaper, insideTaper) .ne. &
              0) call exiterror('Memory allocation error tapering image')
          taperTime = taperTime + wallTime() - wallStart
        endif
        !
        ! if no rescaling, or if mean is needed now, accumulate sums for
        ! mean
        !
        if (.not.rescale .or. ifMean .ne. 0) then
          if (ifFloat == 2) then
            call iclAvgSd(array(iChunkBase), nxOut, numYchunk, 1, nxOut, 1, numYchunk, &
                tmin2, tmax2, tsum, tsumSq, avgSec, sdChunk(ichunk))
            pixChunk(ichunk) = int(nxOut , kind = 8) * numYchunk
            dsumChunk(ichunk) = tsum
            if (iVerbose > 0) print *,'chunk mean&sd min/max', iChunk, avgSec,  &
                sdChunk(ichunk), tmin2, tmax2
            dsumSq = dsumSq + tsumSq
          else
            call iclden(array(iChunkBase), nxOut, numYchunk, 1, nxOut, 1, numYchunk, &
                tmin2, tmax2, tmean2)
            tsum = tmean2 * numPix
          endif
          tmpMin = min(tmpMin, tmin2)
          tmpMax = max(tmpMax, tmax2)
          dsum = dsum + tsum
          if (iVerbose > 0) print *,'did iclden ', tmin2, tmax2, tmpMin, tmpMax
        else
          !
          ! otherwise get new min and max quickly
          !
          tsum = 0.
          do i8 = 1, numPix
            val = array(i8 + iChunkBase - 1)
            if (val < tmpMin) tmpMin = val
            if (val > tmpMax) tmpMax = val
            tsum = tsum + val
          enddo
          dsum = dsum + tsum
        endif
        if (ifFloat == 0 .and. newMode .ne. 2 .and. mode .ne. 2) then
          !
          ! 6/27/01: really want to truncate rather than rescale; so if
          ! the min or max is now out of range for the input mode,
          ! truncate the data and adjust the min and max
          !
          if (tmpMin < bottomIn .or. tmpMax > optimalIn) then
            tsum2 = 0.
            do i8 = 1, numPix
              val = array(i8 + iChunkBase - 1)
              if (val < bottomIn) numTruncLow = numTruncLow + 1
              if (val > optimalIn) numTruncHigh = numTruncHigh + 1
              val = max(bottomIn, min(optimalIn, val))
              tsum2 = tsum2 + val
              array(i8 + iChunkBase - 1) = val
            enddo
            tmpMin = max(tmpMin, bottomIn)
            tmpMax = min(tmpMax, optimalIn)
            dsum = dsum + tsum2 - tsum
          endif
        endif
        !        
        ! If the scaling is pre-set, the scaling should be the same regardless of these
        ! min and max values passed in, and the chunk can be scaled and written
        ! The passed in values are not modified, but tmpMin/Max are, so save and restore
        wallStart = wallTime()
        if (preSetScaling) then
          tmpMin2 = tmpMin
          tmpMax2 = tmpMax
          call findScaleFactors(dminIn, dmaxIn)
          tmpMin = tmpMin2
          tmpMax = tmpMax2
          call scaleAndWriteChunk()
        !
        ! write all but last chunk
        !
        elseif (.not.rescale .or. preSetScaling) then
          if (iVerbose > 0) print *,'writing to real file', iChunk
          call iiuSetPosition(2, isecOut - 1, lineOutSt(iChunk))
          call iiuWriteLines(2, array(iChunkBase), numLinesOut(iChunk))
        elseif (iChunk .ne. numChunks .and. ifOutChunk > 0) then
          if (iVerbose > 0) print *,'writing to temp file', iChunk
          call iiuSetPosition(3, 0, lineOutSt(iChunk))
          call iiuWriteLines(3, array(iBufOutBase), numLinesOut(iChunk))
        endif
        saveTime = saveTime + wallTime() - wallStart
      enddo

      call findScaleFactors(tmpMin, tmpMax)

      if (rescale .and. .not. preSetScaling) then
        dmean2 = 0.
        !
        ! loop backwards on chunks, reloading all but last, scaling
        ! and rewriting
        !
        do iChunk = numChunks, 1, -1
          iChunkBase = iBufOutBase
          if (ifOutChunk == 0)  &
              iChunkBase = iBufOutBase + lineOutSt(iChunk) * int(nxOut, kind = 8)
          if (processInPlace) iChunkBase = 1
          if (iChunk .ne. numChunks .and. ifOutChunk > 0) then
            if (iVerbose > 0) print *,'reading', iChunk
            call iiuSetPosition(3, 0, lineOutSt(iChunk))
            call irdsecl(3, array(iBufOutBase), numLinesOut(iChunk),*99)
          endif
          call scaleAndWriteChunk()
        enddo
      else if (.not. preSetScaling) then
        !
        ! if not scaling
        !
        dmin2 = tmpMin
        dmax2 = tmpMax
        dmean2 = dsum
      endif
      !
      dmean2 = dmean2 / (float(nxOut) * nyOut)
      if (.not. quiet) then
        if (ifHeaderOut == 0) print *, &
            'section   input min&max       output min&max  &  mean'
        ifHeaderOut = 1
        write(*,'(i8,5f10.2)') isec - 1, tmpMin, tmpMax, dmin2, dmax2, dmean2
      endif
      !
80    isecOut = isecOut + 1
      dmin = min(dmin, dmin2)
      dmax = max(dmax, dmax2)
      if (numReplace == 0) then
        dmean = dmean + dmean2
        !
        ! transfer extra header bytes if present
        !
        if (nByteSymOut .ne. 0 .and. indExtraOut < nByteSymOut) then
          !
          ! get this section's size and offset, make sure its OK; only FEI gives error
          if (getExtraHeaderSecOffset(extraIn, nByteSymIn, numIntOrBytesIn, &
              iFlagExtraIn, isecRead, moveOffset, nByteExtraIn) .ne. 0) call exitError( &
              'FEI1 extended header does not contain data for all sections being read')
          if (FEI1type) then
            !
            ! FEI type, call the copy function
            if (copyExtraHeaderSection(extraIn, nByteSymIn, extraOut, nByteSymOut, &
                numIntOrBytesIn, iFlagExtraIn, isecRead, indExtraOut) .ne. 0)  &
                call exitError('Space allowed for copying FEI1 extra header not big'// &
                ' enough; this is due to either program error or mixing files with '// &
                'different sizes per section')
          else
            nByteCopy = min(nByteExtraOut, nByteExtraIn, nByteSymIn)
            numForTilt = 0
            if (saveTilts) then
              !
              ! To save tilt angles, put angle in the integer or real then copy the right
              ! number of bytes; adjust the number to copy and number to clear
              if (serialEMtype) then
                itiltTemp(1) = nint(100. * extraTilts(isec))
                numForTilt = 2
              else
                rtiltTemp = extraTilts(isec)
                numForTilt = 4
              endif
              do i = 1, numForTilt
                indExtraOut = indExtraOut + 1
                extraOut(indExtraOut) = btiltTemp(i)
              enddo
              nByteCopy = max(nByteCopy, numForTilt) - numForTilt
            endif
            nByteClear = nByteExtraOut - numForTilt - nByteCopy
            !
            ! Copy bytes, then clear out the rest if any
            do i = 1, nByteCopy
              indExtraOut = indExtraOut + 1
              extraOut(indExtraOut) = extraIn(iSecRead * nByteExtraIn + i + numForTilt)
            enddo
            do i = 1, nByteClear
              indExtraOut = indExtraOut + 1
              extraOut(indExtraOut) = 0
            enddo
          endif
        endif
        !
        ! Transfer an adoc section
        call int_iwrite(listString, isecOut - 2, ierr)
        if (indAdocIn > 0 .and. indAdocOut > 0) then
          call setCurrentAdocOrExit(indAdocIn, 'input')
          indSectIn = AdocLookupByNameValue(zvalueName, isecRead)
          if (indSectIn > 0) then
            if (AdocTransferSection(zvalueName, indSectIn, indAdocOut, listString, 1) &
                .ne. 0) call exitError('Transferring section data between autodocs')
            outDocChanged = .true.
          endif
        endif
          !
          ! Save tilts in the adoc section: if section does not exist, create it at the
          ! right index; add value to section
        if (indAdocOut > 0 .and. saveTilts) then
          call setCurrentAdocOrExit(indAdocOut, 'output')
          indSectIn = AdocLookupByNameValue(zvalueName, isecOut - 2)
          if (indSectIn <= 0) then
            indSectIn = AdocFindInsertIndex(zvalueName, isecOut - 2)
            if (indSectIn > 0) then
              if (AdocInsertSection(zvalueName, indSectIn, listString) < 0) &
                  indSectIn = -1
            endif
            if (indSectIn < 0)  &
                call exitError('Adding an autodoc section for saving tilt angle')
          endif
          if (AdocSetFloat(zvalueName, indSectIn, 'TiltAngle', extraTilts(isec)) &
              .ne. 0) call exitError('Adding tilt angle to autodoc')
          outDocChanged = .true.
        endif
      else if (isecReplace < listTotal) then
        isecReplace = isecReplace + 1
        isecOut = listReplace(isecReplace) + 1
      endif
      !
      ! see if need to close stack file
      !
      if (numReplace == 0 .and. isecOut > numSecOut(iOutFile)) then
        if (outDocChanged .and. iiuFileType(2) .ne. 5) then
          call setCurrentAdocOrExit(indAdocOut, 'output')
          if (AdocWrite(trim(outFile(iOutFile)) //'.mdoc') .ne. 0) call exitError( &
              'Writing mdoc file for output file')
          call AdocClear(indAdocOut)
        endif
        if (nByteSymOut > 0) call iiuAltExtendedData(2, nByteSymOut, extraOut)
        dmean = dmean / numSecOut(iOutFile)
        call iiuWriteHeader(2, title, 1, dmin, dmax, dmean)
        call iiuClose(2)
        if (needClose2 > 0) call iiuClose(needClose2)
        isecOut = 1
        iOutFile = iOutFile + 1
      endif
      isec = isec + 1
    enddo
    call iiuClose(1)
    if (needClose1 > 0) call iiuClose(needClose1)
  enddo
  if (numReplace > 0) then
    call iiuWriteHeader(2, title, -1, dmin, dmax, dmean)
    call iiuClose(2)
  endif
  !
  if (ifTempOpen .ne. 0) call iiuClose(3)
  if (newAngleFile .ne. ' ') then
    call dopen(3, newAngleFile, 'new', 'f')
    write(3, '(f9.2)')(extraTilts(ind), ind = 1, listTotal)
    close(3)
  endif

  if (numTruncLow + numTruncHigh > 0) write(*,103) numTruncLow, numTruncHigh
103 format(' TRUNCATIONS OCCURRED:',i11,' at low end,',i11, &
      ' at high end of range')
  if (numSecTrunc > 0 .and. &
      numTruncLow + numTruncHigh > numSecTrunc * 4. * (1. + nxOut * (nyOut / 1.e6))) then
    write(*,'(/,a,i4,a,i11,a)') 'WARNING: NEWSTACK - ', numSecTrunc,  &
        ' sections had extreme ranges and were truncated to preserve dynamic range '// &
        '(overall, ', numTruncLow + numTruncHigh, ' pixels were truncated)'
  else if (numSecTrunc > 0) then
    write(*,'(/,a,i4,a)') 'NOTE: ', numSecTrunc,  &
        ' sections had extreme ranges and were truncated to preserve dynamic range '
  endif
  if (iVerbose > 0) write(*,'(a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,f8.4)') 'loadtime', &
      loadTime, '  savetime', saveTime, '  sum', loadTime + saveTime, &
      '  rottime', rotTime, '  taper', taperTime
  call exit(0)
99 call exitError(' End of image while reading')

CONTAINS

  ! Evaluates how much memory is needed for temporary and input/output arrays if no limits
  ! are entered, and reallocates array if necessary
  subroutine reallocateIfNeeded()
    integer(kind = 8) needDim
    integer*4 needTemp, minChunkLines
    real*4 useLimit, inPlaceFac/1./, defLimit/3.75e9/
    !
    ! If system memory available, divide by 4 to get float size
    ! set limit to at least 400 MB, and the minimum of 3/4 of memory and memory - 1GB
    ! Once that gets over the default limit, let it be the max of the default and half of
    ! memory.  But if processing in place and spooling is possible (i.e., not with
    ! FFT's), use 1/4 as much of physical memory
    physicalMem = b3dPhysicalMemory() / 4.
    useLimit = defLimit
    if (processInPlace .and. ftReduceFac == 0 .and. .not. phaseShift) inPlaceFac = 0.25
    if (physicalMem > 0) then
      useLimit = max(0.1e9, min(0.75 * inPlaceFac * physicalMem, physicalMem - 0.25e9))
      if (useLimit > defLimit) useLimit = max(defLimit, 0.5 * inPlaceFac * physicalMem)
    endif
    if (iVerbose > 0) print *,'MB of physical memory', physicalMem / 250000., '  limit', &
        useLimit / 250000.
    !
    ! Set a temp size if no limits entered or if memory limit only entered
    if (limEntered .ne. 1) then
      needTemp = 1
      !
      ! Anticipate irdReduced's needs; 6 is the biggest support width needed for any filter
      if (readShrunk) then
        minChunkLines = 10
        if (readReduction > 32) minChunkLines = 3
        needTemp = max(nx * (ceiling((minChunkLines + 6) * readReduction) + 20),  &
            int(min(int(MAXTEMP, kind = 8), int(nx, kind = 8) * ny)))
      endif
      if (iBinning > 1) needTemp = nx * iBinning
      if (ftReduceFac > 0. .and. nxFSpad > 0) &
          needTemp = max(needTemp, (int(nxFCropPad, kind = 8) + 2) * nyFCropPad)
      if ((phaseShift .or. ftReduceFac > 0.) .and. nxFSpad > 0 .and. noisePad) &
          needTemp = max(needTemp, 2 * max(nxBin, nyBin) + (nxFSpad - nxBin) + &
          (nyFSpad - nyBin))
      lenTemp = needTemp
      !
      ! But if temp size was entered, make sure it is big enough for FT cropping
    else if (ftReduceFac > 0. .and. nxFSpad > 0 .and.  &
        lenTemp < (nxFCropPad + 2) * nyFCropPad) then
      call exitError('Too small a temporary array size entered for Fourier reduction')
    endif
    !
    ! Computed needed size for input and output if no limits entered
    ! Start with the basic input dimension, or padded size for FFT work
    ! Add the output size, or if processing in place, add 2 input lines to keep the test
    ! for input fitting into the array happy
    if (limEntered == 0) then
      needDim = int(nxBin, kind = 8) * nyNeeded
      if ((phaseShift .or. ftReduceFac > 0.) .and. nxFSpad > 0)  &
          needDim = int(nxFSpad + 2, kind = 8) * (nyFSpad + 1)
      if (nxOut > 0 .and. nyOut > 0) then
        if (processInPlace .and. .not. preSetScaling .and.  &
            needDim + 2 * nxBin > useLimit) processInPlace = .false.
        if (processInPlace) then
          needDim = needDim + 2 * nxBin
        else
          needDim = needDim + int(nxOut, kind = 8) * nyOut
        endif
      endif
      if (needDim > useLimit) needDim = useLimit
      if (iVerbose > 0) print *,'reallocate sizes:', nxOut, nyOut, nxBin, nyBin,  &
          needDim, needTemp
      if (needDim + needTemp > limToAlloc) then
        limToAlloc = needDim + needTemp
        call reallocateArray()
        idimInOut = limToAlloc - lenTemp
      else
        idimInOut = needDim
      endif
    else if (limEntered == 2) then
      !
      ! Or at least adjust the in/out size available based on the current limits if only
      ! the memory limit was entered
      idimInOut = limToAlloc - lenTemp
    endif
  end subroutine reallocateIfNeeded

  ! Reallocates the main array to the current required or specified size, and if that
  ! fails, it allocates the fallback amount and makes sure that leaves some input/output
  ! space
  !
  subroutine reallocateArray()
    if (iVerbose > 0) print *, 'reallocating array to', limToAlloc / (1024 * 256.), &
        ' MB'
    deallocate(array)
    allocate(array(limToAlloc), stat = ierr)
    if (ierr .ne. 0 .and. limToAlloc > limIfFail) then
      limToAlloc = limIfFail
      idimInOut = limToAlloc - lenTemp
      if (iVerbose > 0) print *, 'failed, dropping reallocation to', &
          limToAlloc / (1024 * 256), ' MB'
      allocate(array(limToAlloc), stat = ierr)
    endif
    if (ierr .ne. 0) call exitError('Reallocating memory for main array')
    if (limToAlloc - lenTemp < 100) call exitError('With achievable memory'//  &
        ' allocation, the temporary array does not leave enough space for input/output')
    return
  end subroutine reallocateArray


  ! Finds what lines of input are needed to produce lines from LINEOUTFIRST through
  ! LINEOUTLAST of output (numbered from 0), and also determines if x or Y input goes
  ! out of range so that fill is needed
  !
  subroutine linesNeededForOutput(lineOutFirst, lineOutLast, iyIn1, iyIn2, needFill, &
      inPlace)
    integer*4 lineOutFirst, lineOutLast, iyIn1, iyIn2, ixIn1, ixIn2, ixBase
    logical needFill, inPlace
    real*4 xp1, yp1, xp2, yp2, xp3, yp3, xp4, yp4
    if (ftReduceFac > 0.) then
      !
      ! case of fourier cropping: need whole input image
      iyIn1 = 0
      iyIn2 = nyBin - 1
      ixIn1 = 0
      ixIn2 = nxBin - 1
      needFill = abs(xcen(isec)) > 1. .or. abs(ycen(isec)) > 1.
      inPlace = .true.
    elseif (ifXform == 0) then
      !
      ! simple case of no transform
      !
      iyBase = nyBin / 2 + ycen(isec) - (nyOut / 2)
      ixBase = nxBin / 2 + xcen(isec) - (nxOut / 2)
      iyIn1 = max(0, iyBase + lineOutFirst)
      iyIn2 = min(nyBin - 1, iyBase + lineOutLast)
      needFill = iyBase + lineOutFirst < 0 .or. iyBase + lineOutLast >= nyBin .or. &
          ixBase < 0 .or. ixBase + nxOut > nxBin
      inPlace = ixBase >= 0 .and. iyBase + lineOutFirst >= 0
    else
      !
      ! transform: get input needs of 4 corners
      ! pass and get back coordinates numbered from 1, subtract
      ! an extra 1 to get to lines numbered from 0
      ! Allow extra for distortion field Y component
      !
      xcenIn = nxBin / 2.
      ycenIn = nyBin / 2.
      dx = fprod(1, 3)
      dy = fprod(2, 3)
      ! dx=f(1, 3, lnu) -xcen(isec)
      ! dy=f(2, 3, lnu) -ycen(isec)
      call backXform(nxOut, nyOut, fprod, xcenIn , ycenIn, dx, dy, 1, lineOutFirst + 1, &
          xp1, yp1)
      call backXform(nxOut, nyOut, fprod, xcenIn , ycenIn, dx, dy, nxOut,  &
          lineOutFirst + 1, xp2, yp2)
      call backXform(nxOut, nyOut, fprod, xcenIn , ycenIn, dx, dy, 1, lineOutLast + 1, &
          xp3, yp3)
      call backXform(nxOut, nyOut, fprod, xcenIn , ycenIn, dx, dy, nxOut,  &
          lineOutLast + 1, xp4, yp4)

      iyIn1 = int(min(yp1, yp2, yp3, yp4)) - 2 - maxFieldY - linesShrink
      iyIn2 = int(max(yp1, yp2, yp3, yp4)) + 1 + maxFieldY + linesShrink
      ixIn1 = int(min(xp1, xp2, xp3, xp4)) - 2 - maxFieldX - linesShrink
      ixIn2 = int(max(xp1, xp2, xp3, xp4)) + 1 + maxFieldX + linesShrink
      needFill = ixIn1 < 0 .or. ixIn1 >= nxBin .or. ixIn2 < 0 .or. ixIn2 >= nxBin .or. &
          iyIn1 < 0 .or. iyIn1 >= nyBin .or. iyIn2 < 0 .or. iyIn2 >= nyBin
      iyIn1 = min(nyBin - 1, max(0, iyIn1))
      iyIn2 = min(nyBin - 1, max(0, iyIn2))
      inPlace = .false.
    endif

  end subroutine linesNeededForOutput


  ! Determine the scale factors scaleFactor and constAdd from a host of option
  ! settings, controlling values, and values determined for the particular
  ! section.  This code is really dreadful since first it calculates new
  ! min and max and uses that to determine scaling.  It also modifies tmpMin/tmpMax
  !
  subroutine findScaleFactors(tmpMinIn, tmpMaxIn)
    !
    real*4 avgSec, sdSec, dminNew, dmaxNew, dminOut, dmaxOut, zminSec, zmaxSec
    real*4 tmpMean, tmpMinShift, tmpMaxShift, tmpMinIn, tmpMaxIn
    !
    ! calculate new min and max after rescaling under various possibilities
    !
    scaleFactor = 1.
    constAdd = 0.
    tmpMin = tmpMinIn
    tmpMax = tmpMaxIn
    dminOut = tmpMin
    dmaxOut = tmpMax
    !
    if (ifFloat == 0 .and. rescale) then
      !
      ! no float but mode change (not to mode 2) :
      ! rescale from input range to output range
      !
      dminOut = (tmpMin - bottomIn) * (optimalOut - bottomOut) / &
          (optimalIn - bottomIn) + bottomOut
      dmaxOut = (tmpMax - bottomIn) * (optimalOut - bottomOut) / &
          (optimalIn - bottomIn) + bottomOut
    elseif (ifFloat < 0 .and. numScaleFacs == 0) then
      !
      ! if specified global rescale, set values that dminIn and dmaxIn
      ! map to, either the maximum range or the values specified
      !
      if (dminSpecified == 0 .and. dmaxSpecified == 0) then
        dminNew = 0.
        dmaxNew = optimalOut
      else if (dminSpecified == dmaxSpecified) then
        dminNew = dminIn
        dmaxNew = dmaxIn
      else
        dminNew = dminSpecified
        dmaxNew = dmaxSpecified
      endif
      !
      ! then compute what this section's tmpMin and tmpMax map to
      !
      dminOut = (tmpMin - dminIn) * (dmaxNew - dminNew) / (dmaxIn - dminIn) + dminNew
      dmaxOut = (tmpMax - dminIn) * (dmaxNew - dminNew) / (dmaxIn - dminIn) + dminNew
    elseif (ifFloat > 0 .and. ifMeanSdEntered == 0) then
      !
      ! if floating: scale to a dminOut that will knock out fracZero of
      ! the range after truncation to zero
      !
      dminOut = -optimalOut * fracZero / (1. -fracZero)
      if (ifMean == 0) then
        !
        ! float to range, new dmaxOut is the max of the range
        !
        dmaxOut = optimalOut
      elseif (ifFloat == 2) then
        ! :float to mean, it's very hairy
        ! For blank image, use same approach as before to get a 0 sd, otherwise get the
        ! accurate SD from the chunk sums
        if (iSecRead < 0 .or. iSecRead >= nz) then
          call sums_to_avgsd8(dsum, dsumSq, nxOut, nyOut, avgSec, sdSec)
        else
          call chunkSumsToAvgsd(dsumChunk, sdChunk, pixChunk, numChunks, nxOut, nyOut, &
              avgSec, sdSec)
        endif
        if (tmpMin == tmpMax .or. sdSec == 0.) sdSec = 1.
        !
        ! If either the min or the max has become MORE extreme due to interpolation,
        ! then the scaling needs to be reduced by increasing the SD value used to
        ! compute it.  The increase in SD is that which would make the new min or max be
        ! the same number of SDs away from the mean as the old one was
        ind = ilist + listInd(iFile) - 1
        if (.not.(iSecRead < 0 .or. iSecRead >= nz) .and.  &
            (tmpMin < secMins(ind) .or. tmpMax > secMaxes(ind))) then
          boostForMin = 1.
          boostForMax = 1.
          if (avgSec - secMins(ind) > 0.25 * (avgSec - tmpMin)  .and.  &
              avgSec - secMins(ind) < avgSec - tmpMin)  &
              boostForMin = (avgSec - tmpMin) / (avgSec - secMins(ind))
          if (secMaxes(ind) - avgSec > 0.25 * (tmpMax - avgSec) .and.  &
              secMaxes(ind) - avgSec < tmpMax - avgSec)  &
              boostForMax = (tmpMax - avgSec) / (secMaxes(ind) - avgSec)
          sdSec = sdSec * max(boostForMin, boostForMax)
        endif

        ! Truncate the min and max to what will fit in the common range determined after
        ! outlier elimination, then compute the min and max that those map to.
        tmpMin = max(tmpMin, zmin * sdsec + avgSec)
        tmpMax = min(tmpMax, zmax * sdsec + avgSec)
        zminSec = (tmpMin - avgSec) / sdSec
        zmaxSec = (tmpMax - avgSec) / sdSec
        dminOut = (zminSec - zmin) * optimalOut / (zmax - zmin)
        dmaxOut = (zmaxSec - zmin) * optimalOut / (zmax - zmin)
        dminOut = max(0., dminOut)
        dmaxOut = min(dmaxOut, optimalOut)
      else
        !
        ! shift to mean
        !
        tmpMean = dsum / (float(nxOut) * nyOut)
        !
        ! values that min and max shift to
        !
        tmpMinShift = tmpMin + shiftMean - tmpMean
        tmpMaxShift = tmpMax + shiftMean - tmpMean
        !
        if (ifFloat == 3) then
          !
          ! for no specified scaling, set new min and max to
          ! shifted values
          !
          dminOut = tmpMinShift
          dmaxOut = tmpMaxShift
          if (newMode .ne. 2) then
            !
            ! then, if mode is not 2, set up for scaling if range is
            ! too large and/or if there is a modal shift
            !
            optimalIn = max(optimalIn, shiftMax)
            dminOut = tmpMinShift * optimalOut / optimalIn
            dmaxOut = tmpMaxShift * optimalOut / optimalIn
          endif
        else
          !
          ! for specified scaling of shifted means
          !
          if (dminSpecified == dmaxSpecified) then
            dminNew = 0.5
            dmaxNew = optimalOut - 0.5
          else
            dminNew = dminSpecified
            dmaxNew = dmaxSpecified
          endif
          !
          ! for specified scaling, compute what this section's tmpMin
          ! and tmpMax map to
          !
          dminOut = (tmpMinShift - shiftMin) * (dmaxNew - dminNew) / &
              (shiftMax - shiftMin) + dminNew
          dmaxOut = (tmpMaxShift - shiftMin) * (dmaxNew - dminNew) / &
              (shiftMax - shiftMin) + dminNew
        endif
      endif
    endif
    !
    if (rescale) then
      !
      ! if scaling, set up equation, scale and compute new mean
      ! or use scaling factors directly
      !
      if (numScaleFacs > 0) then
        scaleFactor = scaleFacs(min(numScaleFacs, iFile))
        constAdd = scaleConsts(min(numScaleFacs, iFile))
      else if (ifMeanSdEntered .ne. 0) then
        call sums_to_avgsd8(dsum, dsumSq, nxOut, nyOut, avgSec, sdSec)
        if (sdSec == 0) sdSec = 1.
        scaleFactor = enteredSD / sdSec
        constAdd = enteredMean - scaleFactor * avgSec
      else
        !
        ! 2/9/05: keep scale factor 1 if image has no range
        !
        scaleFactor = 1.
        if (dmaxOut .ne. dminOut .and. tmpMax .ne. tmpMin) &
            scaleFactor = (dmaxOut - dminOut) / (tmpMax - tmpMin)
        constAdd = dminOut - scaleFactor * tmpMin
      endif
    endif
    return
  end subroutine findScaleFactors

  subroutine openInputFile(indInFile)
    integer*4 indInFile
    if (indInFile <= numVolRead) then
      call iiAllowMultiVolume(1)
      if (listVolumes(indInFile) > 1) then
        call imopen(11, inFile(indInFile), 'RO')
        if (iiuVolumeOpen(1, 11, listVolumes(indInFile) - 1) .ne. 0)  &
            call exitError('Opening volume in multi-volume file')
        needClose1 = 11
      else
        call imopen(1, inFile(indInFile), 'RO')
        needClose1 = 0
      endif
    else
      call iiAllowMultiVolume(0)
      call imopen(1, inFile(indInFile), 'RO')
    endif
  end subroutine openInputFile


  ! getOffsetEntries gets one or many entries for an offset option
  !
  subroutine getOffsetEntries(option, nameText, xOffset, yOffset, ifOneManyOffsets)
    character*(*) option, nameText
    real*4 xOffset(*), yOffset(*)
    integer*4 ifOneManyOffsets, PipGetFloatArray
    ifOneManyOffsets = 0
    xOffsAll = 0.
    yOffsAll = 0.
    call PipNumberOfEntries(option, numOutEntries)
    if (numOutEntries > 0) then
      ifOneManyOffsets = 1
      numOutValues = 0
      do i = 1, numOutEntries
        numToGet = 0
        ierr = PipGetFloatArray(option, array(numOutValues + 1), numToGet, &
            limSec * 2 - numOutValues)
        numOutValues = numOutValues + numToGet
      enddo
      if (numOutValues .ne. 2 .and. numOutValues .ne. 2 * listTotal) &
          call exitError('There must be either one '//trim(nameText)//' or an '// &
          trim(nameText)//' for each section')
      do i = 1, numOutValues / 2
        xOffset(i) = array(2 * i - 1)
        yOffset(i) = array(2 * i)
      enddo
      if (numOutValues == 2) ifOneManyOffsets = -1
      xOffsAll = xOffset(1)
      yOffsAll = yOffset(1)
    endif
    !
    ! Zero out the lists if one or no entry
    if (ifOneManyOffsets <= 0) then
      do i = 1, listTotal
        xOffset(i) = xOffsAll
        yOffset(i) = yOffsAll
      enddo
    endif
  end subroutine getOffsetEntries


  ! scaleAndWriteChunk applies the defined scaling and truncates data to range, and
  ! accumulates min/max/sum.  Then it writes the chunk
  !
  subroutine scaleAndWriteChunk()
    ! set up minimum value to output based on mode
    if (newMode == 1) then
      densOutMin = -32768
    elseif (newMode == 2) then
      densOutMin = -1.e30
      optimalOut = 1.e30
    else
      densOutMin = 0.
    endif
    do iy = 1, numLinesOut(iChunk)
      istart = iChunkBase + (iy - 1) * int(nxOut, kind = 8)
      tsum = 0.
      do i8 = istart, istart + nxOut - 1
        dens = scaleFactor * array(i8) + constAdd
        if (dens < densOutMin) then
          numTruncLow = numTruncLow + 1
          dens = densOutMin
        elseif (dens > optimalOut) then
          numTruncHigh = numTruncHigh + 1
          dens = optimalOut
        endif
        array(i8) = dens
        tsum = tsum + dens
        dmin2 = min(dmin2, dens)
        dmax2 = max(dmax2, dens)
      enddo
      dmean2 = dmean2 + tsum
    enddo
    wallStart = wallTime()
    if (iVerbose > 0) print *,'writing', iChunk
    call iiuSetPosition(2, isecOut - 1, lineOutSt(iChunk))
    call iiuWriteLines(2, array(iChunkBase), numLinesOut(iChunk))
    saveTime = saveTime + wallTime() - wallStart
  end subroutine scaleAndWriteChunk

end program newstack


! getItemsToUse gets the list of transform line numbers or distortion fields to apply,
! given the number available in NXFORMS, the list of LISTTOTAL section numbers in
! INLIST, the PIP option in OPTION, the PIPINPUT flag if doing pip input, and a scratch
! string in LISTSTRING. ERROR should have a base error string, IFONEPERFILE > 0 to do
! one transform per file, and NUMINFILES should have number of files. NUMBEROFFSET
! should be 0 or 1 depending on whether items are numbered from 1.  The list
! of items is returned in LINEUSE, and the number in the list in NLINEUSE
!
subroutine getItemsToUse(nXforms, listTotal, inList, option, listString, pipinput, &
    error, ifOnePerFile, numInFiles, lineUse, nLineUse, numberOffset, limSec)
  implicit none
  integer*4 nXforms, listTotal, inList(*), numXfLines, lineUse(*), nLineUse, numberOffset
  integer*4 ifOnePerFile, numInFiles
  integer*4 limSec, numLinesTemp, ierr, i, PipGetString
  character*(*) error, option, listString
  character*80 errString
  logical*4 pipinput
  !
  ! Set up default list, add one back if numbered from one
  write(errString, '(a,a,a)') 'Too many ', error, ' numbers for arrays'
  if (nXforms == 1) then
    !
    ! for one transform, set up single line for now
    !
    nLineUse = 1
    lineUse(1) = numberOffset
  elseif (ifOnePerFile > 0) then
    !
    ! for one transform per file, default is 0 to nfile - 1
    !
    nLineUse = numInFiles
    do i = 1, numInFiles
      lineUse(i) = i + numberOffset - 1
    enddo
  else
    !
    ! Otherwise default comes from section list
    !
    nLineUse = listTotal
    do i = 1, listTotal
      lineUse(i) = inList(i) + numberOffset
    enddo
  endif
  !
  numXfLines = 0
  if (pipinput) then
    call PipNumberOfEntries(option, numXfLines)
    if (numXfLines > 0) then
      numLinesTemp = nLineUse
      nLineUse = 0
      do i = 1, numXfLines
        ierr = PipGetString(option, listString)
        call parseList2(listString, lineUse(nLineUse + 1), numLinesTemp, &
            limSec - nLineUse)
        nLineUse = nLineUse + numLinesTemp
      enddo
    endif
  else
    !
    print *,'Enter list of lines to use in file, or a single line number to apply that'
    print *,' transform to all sections (1st line is 0; ranges OK; / for section list)'
    call rdlist2(5, lineUse, nLineUse, limSec)
  endif

  do i = 1, nLineUse
    lineUse(i) = lineUse(i) - numberOffset
    if (lineUse(i) < 0 .or. lineUse(i) >= nXforms) then
      write(errString, '(a, a,i5)') error, ' number out of bounds:', &
          lineUse(i) + numberOffset
      call exitError(trim(errString))
    endif
  enddo
  return
end subroutine getItemsToUse


! irepak2 repacks an image from a portion of a 2-d array sequentially
! into a 1-d array (which should not be the same array) .  Pixels
! outside the range of the original array will be filled with the
! supplied value of dmean.  brray is the repacked array,
! everything else follows definition of iwrpas; i.e. array is
! dimensioned mx by my, and the starting and ending index coordinates
! (numbered from 0) are given by nx1, nx2, ny1, ny2
!
subroutine irepak2(brray, array, mx, my, nx1, nx2, ny1, ny2, dmean)
  implicit none
  integer*4 mx, my, nx1, nx2, ny1, ny2
  real*4 brray(*), array(mx,my), dmean
  integer*4 iy, ix
  integer(kind = 8) ind
  ind = 1
  do iy = ny1 + 1, ny2 + 1
    if (iy >= 1 .and. iy <= my) then
      do ix = nx1 + 1, nx2 + 1
        if (ix >= 1 .and. ix <= mx) then
          brray(ind) = array(ix, iy)
        else
          brray(ind) = dmean
        endif
        ind = ind + 1
      enddo
    else
      do ix = nx1 + 1, nx2 + 1
        brray(ind) = dmean
        ind = ind + 1
      enddo
    endif
  enddo
  return
end subroutine irepak2


! scanSection will determine the min DMIN2, max DMAX2, and mean DMEAN2 of section
! ISECREAD.  It will also determine the standard deviation SDSEC if IFFLOAT = 2.
! It uses ARRAY for storage.  IDIMINOUT specifies the size of ARRAY, while NX is the
! binned image size in x, NYNEEDED is the number of lines to scan, and NEEDYFIRST
! is the first line.  The image will be loaded in chunks if necessary.  LOADYSTART
! and LOADYEND are the starting and ending lines (numbered from 0) that are left
! in ARRAY.
!
subroutine scanSection(array, idimInOut, nx, nyNeeded, needYfirst, reduction, rxOffset, &
    ryOffset, iSecRead, ifFloat, dmin2, dmax2, dmean2, sdSec, loadYstart, loadYend, &
    indFilter, readShrunk, fixRangeSDs, temp, lenTemp)
  implicit none
  integer(kind = 8) idimInOut
  integer*4 nx, iSecRead, ifFloat, loadYstart, loadYend, lenTemp, needYfirst
  integer*4 indFilter
  real*4 array(idimInOut), temp(lenTemp), dmin2, dmax2, dmean2, sdSec
  real*4 reduction, rxOffset, ryOffset, fixRangeSDs
  real*8 sumLoad(1 + nyNeeded / (idimInOut / nx))
  real*4 sdLoad(1 + nyNeeded / (idimInOut / nx))
  real*8 pixLoad(1 + nyNeeded / (idimInOut / nx))
  logical readShrunk
  integer*4 maxLines, numLoads, iline, iload, numLines, nyNeeded
  real*4 tmin2, tmax2, tmean2, avgSec
  real*8 dsum, tsumSq
  !
  ! load in chunks if necessary, based on the maximum number
  ! of lines that will fit in the array
  !
  maxLines = idimInOut / nx
  numLoads = (nyNeeded + maxLines - 1) / maxLines
  iline = needYfirst
  dmin2 = 1.e30
  dmax2 = -dmin2
  dsum = 0.
  do iload = 1, numLoads
    numLines = nyNeeded / numLoads
    if (iload <= mod(nyNeeded, numLoads)) numLines = numLines + 1
    call readBinnedOrReduced(1, iSecRead, array, nx, numLines, rxOffset, ryOffset + &
        reduction * iline, reduction, nx, numLines, indFilter, readShrunk, temp, lenTemp)
    !
    ! accumulate sums for mean and sd if float 2, otherwise
    ! just the mean.  Store the floating SD's and # of pixels, compute more accurate mean
    !
    if (ifFloat == 2 .or. fixRangeSDs > 0.) then
      call iclAvgSd(array, nx, numLines, 1, nx, 1, numLines, tmin2, tmax2, &
          sumLoad(iload), tsumSq, dmean2, sdLoad(iload))
      pixLoad(iload) = int(nx, kind = 8) * numLines
    else
      call iclden(array, nx, numLines, 1, nx, 1, numLines, tmin2, &
          tmax2, tmean2)
      sumLoad(iload) = (tmean2 * nx) * numLines
    endif
    dmin2 = min(dmin2, tmin2)
    dmax2 = max(dmax2, tmax2)
    dsum = dsum + sumLoad(iload)
    iline = iline + numLines
  enddo
  !
  ! Compute accurate overall mean, and then combine all the mean and SD data to get the SD
  dmean2 = (dsum / nx) / nyNeeded
  if (ifFloat == 2 .or. fixRangeSDs > 0.) then
    call chunkSumsToAvgsd(sumLoad, sdLoad, pixLoad, numLoads, nx, nyNeeded, dmean2, sdSec)
  endif
  loadYend = iline-1
  loadYstart = iline - numLines
  return
end subroutine scanSection


! backXform will determine the array coordinates XP and YP in an input
! image that transform to the given array indexes IX, IY in an output
! image.  Other arguments match those of cubInterp: the input image is
! NXA by NYA; the output image is NXB by NYB; XFCEN, YFCEN is the center
! coordinate of the input image; AMAT is the 2x2 transformation matrix
! and XTRANS, YTRANS are the translations.
!
subroutine backXform(nxb, nyb, amat, xfcen, yfcen, xtrans, ytrans, ix, iy, xp, yp)
  implicit none
  integer*4 nxb, nyb, ix, iy
  real*4 amat(2,2), xfcen, yfcen, xtrans, ytrans, xp, yp
  real*4 xcen, ycen, xcenOut, ycenOut, denom, a11, a12, a21, a22, dyo, dxo
  !
  ! Calc inverse transformation
  !
  xcen = nxb / 2. + xtrans + 0.5
  ycen = nyb / 2. + ytrans + 0.5
  xcenOut = xfcen + 0.5
  ycenOut = yfcen + 0.5
  denom = amat(1, 1) * amat(2, 2) - amat(1, 2) * amat(2, 1)
  a11 =  amat(2, 2) / denom
  a12 = -amat(1, 2) / denom
  a21 = -amat(2, 1) / denom
  a22 =  amat(1, 1) / denom
  !
  ! get coordinate transforming to ix, iy
  !
  dyo = iy - ycen
  dxo = ix - xcen
  xp = a11 * dxo + a12 * dyo + xcenOut
  yp = a21 * dxo + a22 * dyo + ycenOut
  return
end subroutine backXform


! getReducedSize returns the input size and offset in one dimension when there is binning
! or shrinkage reduction
!
subroutine getReducedSize(nx, reduction, doShrink, nxBin, xOffset)
  implicit none
  integer*4 nx, nxBin, ixOffset
  real*4 reduction, xOffset
  logical doShrink
  !
  ! For non-integer shrinkage, just cut the size, otherwise match the binned
  ! size for consistency in tilt series processing
  if (doShrink .and. abs(nint(reduction) - reduction) > 1.e-4) then
    nxBin = nx / reduction
    xOffset = (nx - nxBin * reduction) / 2.
  else
    call getBinnedSize(nx, nint(reduction), nxBin, ixOffset)
    xOffset = ixOffset
  endif
  return
end subroutine getReducedSize


! readBinnedOrReduced will read an area from the file with either binning or reduction
!
subroutine readBinnedOrReduced(imUnit, iz, array, nxDim, nyDim, xUBstart, yUBstart,  &
    redFac, nxRed, nyRed, ifiltType, doShrink, temp, lenTemp)
  implicit none
  integer*4 imUnit, iz, nxDim, nyDim, nxRed, nyRed, ifiltType, lenTemp, ierr
  real*4 array(*), xUBstart, yUBstart, temp(*), redFac
  logical doShrink
  character*100 listString

  if (doShrink) then
    call irdReduced(imUnit, iz, array, nxDim, xUBstart, yUBstart, redFac, nxRed, &
        nyRed, ifiltType, temp, lenTemp, ierr)
    if (ierr > 0) then
      write(listString, '(a,i2,a)') 'Calling irdReduced to read image (error code', &
          ierr, ')'
      call exitError(listString)
    endif
  else
    call irdBinned(imUnit, iz, array, nxDim, nyDim, nint(xUBstart), nint(yUBstart), &
        nint(redFac), nxRed, nyRed, temp, lenTemp, ierr)
  endif
  if (ierr .ne. 0) call exitError('Reading image file')
end subroutine readBinnedOrReduced


! chunkSumsToAvgsd takes a set of sums (as doubles) and SD values (as floats) and the
! number of pixels in each chunk (as doubles) and two numbers whose product is the total
! number of pixels, and computes the mean and an accurate SD
!
subroutine chunkSumsToAvgsd(dsumChunk, sdChunk, pixChunk, numChunks, nx, ny, avgSec, &
    sdSec)
  implicit none
  real*8 dsumChunk(*), pixChunk(*)
  real*4 sdChunk(*), avgSec, sdSec
  integer*4 numChunks, nx, ny, ichunk
  real*8 pixTot, dmean8, dsum, dsumSq, avgChunk(numChunks)
  dsum = 0.
  pixTot = int(nx, kind = 8) * ny
  do ichunk = 1, numChunks
    dsum = dsum + dsumChunk(ichunk)
    avgChunk(ichunk) = dsumChunk(ichunk) / pixChunk(ichunk)
  enddo
  dmean8 = dsum / pixTot
  avgSec = dmean8
  dsumSq = 0.
  do ichunk = 1, numChunks
    dsumSq = dsumSq + pixChunk(ichunk) * (avgChunk(ichunk)**2 - dmean8**2) + &
        (pixChunk(ichunk) - 1.) * sdChunk(ichunk)**2
  enddo
  sdSec = sqrt(max(0., dsumSq / max(1., pixTot - 1.)))
  return
end subroutine chunkSumsToAvgsd
