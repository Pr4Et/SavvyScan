! * * * * * * * * BLENDMONT.FOR * * * * * * *
!
! BLENDMONT will take montaged images, blend their overlapping edges
! together, and output the blended images with essentially no overlap.
! Transformations can be applied to align serial sections, and
! densities can be floated over whole sections so that density
! occupies the maximum range for each section.  The program can also
! correct for minor or substantial displacements between pieces, using
! cross-correlation to find substantial shifts.  Multinegative montages
! can also be handled, with each negative transformed so as to produce
! the best fit between the negatives.
!
! See the man page for further details.
!
! David Mastronarde, February 1989
! 12/21/98 added capability to do initial cross-correlation of overlaps
! and correct for sloppy montages
! 12/99: revamped the treatment of overlap zones to avoid artifacts in
! sloppy montages
! 8/22/00: made it rename existing edge function files if 'new' ones
! are requested and they already exist
! 12/11/00: implemented byte-swapping when reading edge function files
! built on other machines
! 6/1/01: implemented ability to write and read in edge correlation
! displacements; added a search step to improve on cross-correlations
!
! $Id$
!
program blendmont
  use blendvars
  implicit none
  integer nbytes_recl_item, LIMNEG
  include 'recl_bytes.inc'
  character*320 fileName, edgeName, plOutFile, outFile
  character*320 rootName, edgeName2, aliCoordFile
  character*4 edgeExtension(2) /'.xef', '.yef'/
  character*5 xcorrExtension(0:2) /'.ecd', '.xecd', '.yecd'/
  character*18 actionStr
  integer*4 mxyzIn(3), nxyzst(3) /0, 0, 0/
  real*4 cell(6) /1., 1., 1., 0., 0., 0./
  real*4 delta(3), xOrigin, yOrigin, zOrigin
  real*4, allocatable :: binLine(:)
  character*80 titleStr
  !
  ! structure /rotrans/
  ! real*4 theta, dx, dy, xcen, ycen
  ! end structure
  !
  integer*4 nskipRegress(2) /1, 2/                   !regression position skip s&l
  integer*4 numFit(2) /5, 7/                    !# to include in regression s&l
  integer*4 igridStart(2), iOffset(2)           !just temp usage here
  !
  parameter (LIMNEG = 30)
  integer*4 minNegX(LIMNEG), maxNegX(LIMNEG) !min and max coords of
  integer*4 minNegY(LIMNEG), maxNegY(LIMNEG) !each neg
  integer*4 jointLower(LIMNEG,2), jointUpper(LIMNEG,2) !joints around neg
  integer*4 negLower(LIMNEG,2), negUpper(LIMNEG,2) !negs around joint
  integer*4 negIndex(LIMNEG)                  !index to arbitrary neg #'s
  integer*4 numInHAdj(LIMNEG), numJoint(2)
  integer*4 numEdgeOnJoint(15), listOnJoint(15,LIMNEG) !#, list edges on joint
  integer*4 ixPcLower(15), ixPcUpper(15), iyPcLower(15), iyPcUpper(15)
  !
  integer*4 numEdgeTmp(5,2)
  integer*4, allocatable :: listZ(:), izWant(:), izAllWant(:)
  logical skipXforms, undistortOnly, xcWriteOut, sameEdgeShifts
  character dat * 9, tim * 8
  real*4 edgeFrac(2), title(20)
  real*4 edgeFracX, edgeFracY
  equivalence (edgeFrac(1), edgeFracX), (edgeFrac(2), edgeFracY)
  integer*4, allocatable :: mapAllPiece(:,:,:), numControl(:)
  integer*4, allocatable :: ixPcTemp(:), iyPcTemp(:), izPcTemp(:), negTmp(:)
  logical, allocatable :: multiNeg(:), multiTemp(:), edgeDone(:,:)
  logical active4(3, 2)
  logical anyPixels, inFrame, doFast, anyNeg, anyLinesOut, xIsLongDim, testMode
  logical shiftEach, doCross, fromEdge, exist, xcReadIn, xcLegacy, outputpl, parallelHDF
  logical xcorrDebug, verySloppy, useEdges, edgesIncomplete, adjustOrigin, useFill
  logical edgesSeparated, fromCorrOnly, samePieces, noFFTsizes, yChunks, doWarp
  real*4, allocatable :: dxGridMean(:,:), dyGridMean(:,:)
  real*4, allocatable :: edgeDisplaceX(:,:), edgeDisplaceY(:,:)
  real*4 afterMean(2), afterMax(2)
  character*6 edgeXcorrText(2) /'xcorr:', 'edges:'/
  real*4 fastf(2,3), gxfTemp(2,3), fastTemp(2,3)
  real*4, allocatable :: hxf(:,:,:), gxf(:,:,:)
  ! the former rotrans structures
  real*4 hCum(6,LIMNEG), hAdj(6,LIMNEG), r(6,LIMNEG,2), hFrame(6), rotxNet(6)
  integer*4 modePower(0:15) /8, 15, 8, 0, 0, 0, 16, 0, 0, 9, 10, 11, 12, 13, 14, 15/
  integer*4  minXwant, maxXwant, modeParallel, nzAllWant, numOut, numChunks
  integer*4  minYwant, maxYwant, numExtra, nxyBox(2), nExtra(2)
  !
  character*320 boundFile
  integer*4 idfNx, idfNy, idfBinning
  real*4 pixelIdf, binRatio, binningOfInput
  !
  integer*4 modeIn, numListZ, minZpiece, maxPieceZ, numSect, nxTotalPix, nyTotalPix
  integer*4 modeOut, ifFloat, i, minXoverlap
  integer*4 minYoverlap, numTrials, ifSloppy, ioptAbs, numXmissing, numYmissing
  integer*4 nxTotalWant, nytotwant, newXpieces, newYpieces, ifOldEdge
  integer*4 newXtotalPix, newYtotalPix, newXframe, newYframe, newMinXpiece
  integer*4 newMinYpiece, ifWant, numGxforms, ierr, numXframesPerNeg, numYframesPerNeg
  real*4 dmin, dmax, outMin, outMax, dfltInMin, pixelTot, fillVal
  real*4 dfltInMax, curInMin, curInMax, pixelScale, pixelAdd
  real*8 tsum, curSum, grandSum, realNsum
  integer*4 ixFrame, iyFrame, ipc, numShort, numLong, ishort, ilong, indArray
  integer*4 newUseCount, newPcXlowLeft, newPcYlowLeft, nxFast, nyFast, iyFast, ixFast
  integer*4 indYlow, indYhigh, numLinesOut, indXlow, indXhigh, inOnePiece, ioptNeg
  integer*4 indx, indy, ixNeg, iyNeg, numXneg, iz, ineg, listFirst, ipcHigh
  integer*4 ix, iy, ilineOut, ifill, numAlong, numAcross, numEd, ixy, lenRecord, j
  real*4 dminOut, dmaxOut, tmean, sdCrit, devCrit, gridScale
  integer*4 numZwant, newXoverlap, newYoverlap, izSect, iwant, ifDiddle
  integer*4 ixOut, iyOut, ifRevise, ipolyOrder, iedgeDir, iyx, imem, jedge, numNegatives
  integer*4 iwhich, negLow, negUp, indLow, indUpper, joint, ied, iedge, ipcLower, ipcUpper
  integer*4 maxSwing, ishift, numShiftNeg, inde, numBestEdge, indBest
  real*4 errAdjust, errLim, dxSum, dySum, sumX, sumY, xDisplace, yDisplace
  real*4 beforeMean, beforeMax, warpScale, warpXoffset, warpYoffset
  integer*4 indGxf, ilis, numIter, lineBase, numEdgesIn, indEdge, indLower
  real*4 x1, y1, w1, w2, c11, c12, c21
  real*4 xSrcConst, ySrcConst, c22, denom, fb11, fb12, fb21, fb22, x1last, y1last
  integer*4 indPcOne, indPcLo, indPcUp, ipiece1, ipiece2, iter, jndP1, jndP2, jndP3
  real*4 dx2, dy2, x2, y2, pixVal, w3, wMax, f2b11, f2b12, f2b21, f2b22, val
  integer*4 ind12edg, ind13Edge, iCornType, ind23edg, ind14edg, ind43edg
  real*4 f3b11, f3b12, f3b21, f3b22, dden, dden12, dden13, dden23, x3, y3, x3t, y3t
  real*4 dy3, dx3, bx, by, emin, w4, f4b11, f4b12, f4b21, f4b22, dden14, dden43
  integer*4 ipiece3, ipiece4, nxGridIn, nyGridIn, indArray1, indArray2, jndP4
  real*4 x4, y4, dx4, dy4, xg, ysrc, delDmagPerUm, delRotPerUm, delIndent(2), xtmp
  real*4 dmagNew, drotNew, tiltOffset, edgeStart, edgeEnd, bwOffset, ecdBinning, warpPixel
  integer*4 iBinning, nyWrite, indentXC, numZero, mostNeeded
  integer*4 lineOffset, ixOffset, iyOffset, linesBuffered, iBufferBase
  integer*4 maxXcorrBinning, nxyXcorrTarget, numSkip, lineStart, lineEnd
  integer*4 nxyPadded, nxyBoxed, ifUseAdjusted, iyOutOffset, ifEdgeFuncOnly
  integer*4 ipFirst(4), numFirst, ixyFuncStart, ixyFuncEnd, numLinesWrite
  integer*4 iedgeDelX, iedgeDelY, ixUnaliStart, iyUnaliStart
  integer*4 iwarpNx, iwarpNy, indWarpFile
  integer*4 imodBackupFile, parWrtInitialize, getWarpGrid, readCheckWarpFile
  integer*4 getGridParameters, findMaxGridSize, getSizeAdjustedGrid, getLinearTransform
  integer*4 setCurrentWarpFile, montXCFindBinning, niceFFTlimit
  integer*4 b3dLockFile, b3dOutputFileType, iiuFileType, iiTestIfHDF,iiuParWrtFlushBuffers
  real*8 wallTime, wallStart, fastCum, slowCum
  !
  logical pipInput
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetInteger, PipGetBoolean, PipGetThreeFloats
  integer*4 PipGetString, PipGetFloat, PipGetIntegerArray
  integer*4 PipGetTwoIntegers, PipGetTwoFloats, PipGetLogical
  integer*4 PipGetInOutFile, PipNumberOfEntries
  !
  ! fallbacks from ../../manpages/autodoc2man -3 2 blendmont
  integer numOptions
  parameter (numOptions = 73)
  character*(40 * numOptions) options(1)
  options(1) = &
      'imin:ImageInputFile:FN:@plin:PieceListInput:FN:@imout:ImageOutputFile:FN:@'// &
      'plout:PieceListOutput:FN:@aligned:AlignedPieceCoordFile:FN:@'// &
      'rootname:RootNameForEdges:CH:@oldedge:OldEdgeFunctions:B:@'// &
      'perneg:FramesPerNegativeXandY:IP:@missing:MissingFromFirstNegativeXandY:IP:@'// &
      'mode:ModeToOutput:I:@float:FloatToRange:B:@fill:FillValue:F:@'// &
      'xform:TransformFile:FN:@center:TransformCenterXandY:FP:@'// &
      'unaligned:UnalignedStartingXandY:IP:@order:InterpolationOrder:I:@'// &
      'sections:SectionsToDo:LI:@xminmax:StartingAndEndingX:IP:@'// &
      'yminmax:StartingAndEndingY:IP:@nofft:NoResizeForFFT:B:@origin:AdjustOrigin:B:@'// &
      'bin:BinByFactor:I:@maxsize:MaximumNewSizeXandY:IP:@'// &
      'minoverlap:MinimumOverlapXandY:IP:@distort:DistortionField:FN:@'// &
      'imagebinned:ImagesAreBinned:F:@gradient:GradientFile:FN:@'// &
      'adjusted:AdjustedFocus:B:@addgrad:AddToGradient:FP:@tiltfile:TiltFile:FN:@'// &
      'offset:OffsetTilts:F:@geometry:TiltGeometry:FT:@'// &
      'justUndistort:JustUndistort:B:@test:TestMode:B:@sloppy:SloppyMontage:B:@'// &
      'very:VerySloppyMontage:B:@shift:ShiftPieces:B:@edge:ShiftFromEdges:B:@'// &
      'xcorr:ShiftFromXcorrs:B:@readxcorr:ReadInXcorrs:B:@'// &
      'ecdbin:BinningForEdgeShifts:F:@overlap:OverlapForEdgeShifts:IP:@'// &
      'skip:SkipEdgeModelFile:FN:@nonzero:NonzeroSkippedEdgeUse:I:@'// &
      'robust:RobustFitCriterion:F:@width:BlendingWidthXandY:IP:@'// &
      'boxsize:BoxSizeShortAndLong:IP:@grid:GridSpacingShortAndLong:IP:@'// &
      'indents:IndentShortAndLong:IP:@goodedge:GoodEdgeLowAndHighZ:IP:@'// &
      'onegood:OneGoodEdgeLimits:IAM:@same:SameEdgeShifts:B:@'// &
      'exclude:ExcludeFillFromEdges:B:@unsmooth:UnsmoothedPatchFile:FN:@'// &
      'smooth:SmoothedPatchFile:FN:@parallel:ParallelMode:IP:@subset:SubsetToDo:LI:@'// &
      'lines:LineSubsetToDo:IP:@boundary:BoundaryInfoFile:FN:@'// &
      'functions:EdgeFunctionsOnly:I:@aspect:AspectRatioForXcorr:F:@'// &
      'pad:PadFraction:F:@extra:ExtraXcorrWidth:F:@numpeaks:NumberOfXcorrPeaks:I:@'// &
      'radius1:FilterRadius1:F:@radius2:FilterRadius2:F:@sigma1:FilterSigma1:F:@'// &
      'sigma2:FilterSigma2:F:@treat:TreatFillForXcorr:I:@xcdbg:XcorrDebug:B:@'// &
      'taper:TaperFraction:F:@param:ParameterFile:PF:@help:usage:B:'
  !
  ! initialization of many things
  !
  ixDebug = -12000
  iyDebug = -12000
  inPiece(0) = 0
  iunEdge(1) = 7
  iunEdge(2) = 8
  indent(1) = 3                               !minimum indent short & long
  indent(2) = 3
  intGrid(1) = 6                              !grid interval short & long
  intGrid(2) = 6
  iboxSiz(1) = 10                             !box size short & long
  iboxSiz(2) = 15
  ifFloat = 0
  ifSloppy = 0
  ioptAbs = 0
  numXmissing = 0
  numYmissing = 0
  ifOldEdge = 0
  shiftEach = .false.
  xcLegacy = .false.
  fromEdge = .false.
  xcReadIn = .false.
  doMagGrad = .false.
  undistort = .false.
  focusAdjusted = .false.
  binningOfInput = 1
  interpOrder = 2
  testMode = .false.
  undistortOnly = .false.
  limitData = .false.
  adjustOrigin = .false.
  noFFTsizes = .false.
  yChunks = .false.
  sameEdgeShifts = .false.
  parallelHDF = .false.
  iBinning = 1
  numAngles = 0
  numUseEdge = 0
  izUseDefLow = -1
  izUseDefHigh = -1
  modeParallel = 0
  outFile = ' '
  numXcorrPeaks = 1
  numSkip = 0
  ifUseAdjusted = 0
  iyOutOffset = 0
  ifEdgeFuncOnly = 0
  edgeName2 = ' '
  boundFile = ' '
  izUnsmoothedPatch = -1
  izSmoothedPatch = -1
  robustCrit = 0.
  ecdBinning = 1.
  !
  ! Xcorr parameters
  ! 11/5/05: increased taper fraction 0.05->0.1 to protect against
  ! edge effects with default filter
  !
  xcorrDebug = .false.
  ifDumpXY(1) = -1
  ifDumpXY(2) = -1
  ifillTreatment = 1
  aspectMax = 2.0                             !maximum aspect ratio of block
  padFrac = 0.45
  extraWidth = 0.
  radius1 = 0.
  radius2 = 0.35
  sigma1 = 0.05
  sigma2 = 0.05
  maxXcorrBinning = 3
  nxyXcorrTarget = 1024
  verySloppy = .false.
  lmField = 1
  maxFields = 1
  useFill = .false.
  memLim = 256
  !
  ! Pip startup: set error, parse options, check help, set flag if used
  !
  call PipReadOrParseOptions(options, numOptions, 'blendmont', &
      'ERROR: BLENDMONT - ', .true., 0, 0, 0, numOptArg, &
      numNonOptArg)
  pipInput = numOptArg + numNonOptArg > 0
  !
  if (PipGetInOutFile('ImageInputFile', numNonOptArg + 1, &
      'Input image file', fileName) .ne. 0) call exitError( &
      'No input image file specified')
  call imopen(1, fileName, 'ro')
  call irdhdr(1, nxyzIn, mxyzIn, modeIn, dmin, dmax, dmean)
  call iiuRetDelta(1, delta)
  if (pipInput) ierr = PipGetInteger('EdgeFunctionsOnly', ifEdgeFuncOnly)
  ifEdgeFuncOnly = max(0, min(3, ifEdgeFuncOnly))
  if (ifEdgeFuncOnly == 1 .or. ifEdgeFuncOnly == 2) then
    ixyFuncStart = ifEdgeFuncOnly
    ixyFuncEnd = ifEdgeFuncOnly
  else
    ixyFuncStart = 1
    ixyFuncEnd = 2
  endif
  !
  if (PipGetInOutFile('ImageOutputFile', numNonOptArg + 1, &
      'Output image file', outFile) .ne. 0 .and. ifEdgeFuncOnly == 0) &
      call exitError('No output image file specified')
  !
  modeOut = modeIn
  if (pipInput) then
    ierr = PipGetInteger('ModeToOutput', modeOut)
  else
    write(*,'(1x,a,i2,a,$)') 'Mode for output file [/ for', modeOut, ']: '
    read(5,*) modeOut
  endif
  if (modeOut < 0 .or. modeOut > 15 .or. &
      (modeOut >= 3 .and. modeOut <= 8 .and. modeOut .ne. 6)) call exitError( &
      'Bad mode value')
  !
  ! set up output range, and default input range and minimum.  If real
  ! input, use actual input range and min; otherwise use theoretical
  ! range for the input mode.  This will preserve data values if no float
  !
  outMin = 0.
  outMax = 2**modePower(modeOut) - 1.
  if (modeOut == 1) outMin = -outMax
  if (modeIn == 2) then
    dfltInMax = dmax
    dfltInMin = dmin
    if (modeOut == 2) then
      outMin = dmin
      outMax = dmax
    endif
  else
    dfltInMax = 2**modePower(modeIn) - 1.
    dfltInMin = 0.
    if (modeIn == 1) dfltInMin = -dfltInMax
  endif
  !
  fileName = ' '
  if (pipInput) then
    ierr = PipGetBoolean('FloatToRange', ifFloat)
    ierr = PipGetString('TransformFile', fileName)
    ierr = PipGetLogical('JustUndistort', undistortOnly)
  else
    write(*,'(1x,a,$)') &
        '1 to float each section to maximum range, 0 not to: '
    read(5,*) ifFloat
    !
    write(*,'(1x,a,$)') &
        'File of g transforms to apply (Return if none): '
    read(5, '(a)') fileName
  endif
  !
  ! Preserve legacy behavior of floating to positive range for mode 1
  if (ifFloat .ne. 0 .and. modeIn == 1 .and. modeOut == 1) then
    outMin = 0.
    dfltInMin = 0.
  endif

  limSect = 100000
  allocate(ixPcTemp(limInit), iyPcTemp(limInit), izPcTemp(limInit), &
      negTmp(limInit), multiTemp(limSect), stat = ierr)
  if (ierr .ne. 0) call exitError('Allocating initial arrays')
  !
  call read_list(ixPcTemp, iyPcTemp, izPcTemp, negTmp, &
      multiTemp, npcList, minZpiece, maxPieceZ, anyNeg, pipInput)
  limNpc = npcList + 10
  allocate(ixPcList(limNpc), iyPcList(limNpc), izPcList(limNpc), &
      negList(limNpc), limDataInd(limNpc), iedgeLower(limNpc, 2), &
      iedgeUpper(limNpc, 2), hinv(2, 3, limNpc), memIndex(limNpc), &
      htmp(2, 3, limNpc), hxf(2, 3, limNpc), indVar(limNpc), stat = ierr)
  if (ierr .ne. 0) call exitError('Allocating arrays per piece')
  ixPcList(1:npcList) = ixPcTemp(1:npcList)
  iyPcList(1:npcList) = iyPcTemp(1:npcList)
  izPcList(1:npcList) = izPcTemp(1:npcList)
  negList(1:npcList) = negTmp(1:npcList)
  !
  numSect = maxPieceZ + 1 - minZpiece
  limSect = numSect + 1
  allocate(tiltAngles(limSect), dmagPerUm(limSect), rotPerUm(limSect), &
      listZ(limSect), izWant(limSect), izAllWant(limSect), gxf(2, 3, limSect), &
      multiNeg(limSect), izMemList(memLim), lastUsed(memLim), stat = ierr)
  if (ierr .ne. 0) call exitError('Allocating arrays per piece')
  multiNeg(1:numSect) = multiTemp(1:numSect)
  !
  call fill_listz(izPcList, npcList, listZ, numListZ)
  !
  doGxforms = .false.
  doWarp = .false.
  if (fileName .ne. ' ' .and. .not.undistortOnly) then
    doGxforms = .true.
    !
    ! Open as warping file if possible
    indWarpFile = readCheckWarpFile(fileName, 0, 1, iwarpNx, iwarpNy, numGxforms, ix, &
        warpPixel, iy, edgeName)
    if (indWarpFile < -1) call exitError(edgeName)
    doWarp = indWarpFile >= 0
    if (doWarp) then
      if (numGxforms > limSect)  call exitError( &
          'Too many sections in warping file for transform array')
      warpScale = warpPixel / delta(1)
      do i = 1, numGxforms
        if (getLinearTransform(i, gxf(1, 1, i)) .ne. 0) &
            call exitError('Getting linear transform from warp file')
        gxf(1, 3, i) = gxf(1, 3, i) * warpScale
        gxf(2, 3, i) = gxf(2, 3, i) * warpScale
      enddo
    else
      !
      ! Get regular transforms
      call dopen(3, fileName, 'ro', 'f')
      call xfrdall2(3, gxf, numGxforms, limSect, ierr)
      close(3)
      !
      ! It is OK if ierr=-1 : more transforms than sections
      if (ierr > 0) call exitError ('Reading transforms')
    endif
  endif


  if (doGxforms) then
    if (numListZ > numGxforms) call exitError('More sections than G transforms')

    skipXforms = .false.
    if (numListZ < numSect) then
      if (numGxforms == numListZ) then
        print *,'Looks like there are transforms only for '// &
            'the sections that exist in file'
      elseif (numGxforms == numSect) then
        print *,'There seem to be transforms for each Z value,' &
            //' including ones missing from file'
        skipXforms = .true.
      else
        call exitError('Cannot tell how transforms match up to '// &
            'sections, because of missing sections')
      endif
    endif
  endif
  !
  ! now check lists and get basic properties of overlap etc
  !
  call checklist(ixPcList, npcList, 1, nxin, minXpiece, nxPieces, &
      nXoverlap)
  call checklist(iyPcList, npcList, 1, nyin, minYpiece, nyPieces, &
      nyOverlap)
  if (nxPieces <= 0 .or. nyPieces <= 0) call exitError &
      ('Checklist reported a problem with the piece list in one direction')
  !
  nxTotalPix = nxPieces * (nxin - nXoverlap) + nXoverlap
  nyTotalPix = nyPieces * (nyin - nyOverlap) + nyOverlap
  if (pipInput) print *,'Input file:'
  write(*,115) nxTotalPix, 'X', nxPieces, nxin, nXoverlap
  write(*,115) nyTotalPix, 'Y', nyPieces, nyin, nyOverlap
115 format(i7,' total ',a1,' pixels in',i4,' pieces of', &
      i6, ' pixels, with overlap of',i5)
  !
  limEdge = limSect * max(1, nxPieces * nyPieces)
  allocate(dxGridMean(limEdge, 2), dyGridMean(limEdge, 2), &
      edgeDisplaceX(limEdge, 2), edgeDisplaceY(limEdge, 2), edgeDone(limEdge, 2), &
      ipieceLower(limEdge, 2), ipieceUpper(limEdge, 2), ibufEdge(limEdge, 2), &
      ifSkipEdge(limEdge, 2), stat = ierr)
  if (ierr .ne. 0) call exitError('Allocating arrays per edge')
  !
  ! find out if global multi-neg specifications are needed or desired
  ! But first deal with correlation control parameters
  ! Here are the defaults for VerySloppy
  !
  if (pipInput) then
    ierr = PipGetLogical('VerySloppyMontage', verySloppy)
    if (verySloppy) then
      ifSloppy = 1
      aspectMax = 5.
      radius1 = -0.01
      extraWidth = 0.25
      numXcorrPeaks = 16
    else
      ierr = PipGetBoolean('SloppyMontage', ifSloppy)
    endif
    ierr = PipGetFloat('AspectRatio', aspectMax)
    ierr = PipGetInteger('NumberOfXcorrPeaks', numXcorrPeaks)
    numXcorrPeaks = max(1, min(limXcorrPeaks, numXcorrPeaks))
    ierr = PipGetFloat('PadFraction', padFrac)
    ierr = PipGetFloat('ExtraXcorrWidth', extraWidth)
    ierr = PipGetFloat('FilterSigma1', sigma1)
    ierr = PipGetFloat('FilterSigma2', sigma2)
    ierr = PipGetFloat('FilterRadius1', radius1)
    ierr = PipGetFloat('FilterRadius2', radius2)
    shiftEach = ifSloppy .ne. 0
    if (PipGetTwoIntegers('FramesPerNegativeXandY', numXframesPerNeg, &
        numYframesPerNeg) == 0) ioptAbs = 1
    if (.not.shiftEach) ierr = PipGetLogical('ShiftPieces', shiftEach)
    ierr = PipGetLogical('ShiftFromEdges', fromEdge)
    ierr = PipGetLogical('ShiftFromXcorrs', xcLegacy)
    ierr = PipGetLogical('ReadInXcorrs', xcReadIn)
    ierr = PipGetInteger('NonzeroSkippedEdgeUse', ifUseAdjusted)
    ierr = PipGetFloat('RobustFitCriterion', robustCrit)
    ierr = PipGetLogical('TestMode', testMode)
    useFill = PipGetFloat('FillValue', fillVal) == 0
    ierr = PipGetTwoIntegers('GoodEdgeLowAndHighZ', izUseDefLow, &
        izUseDefHigh)
    ierr = PipNumberOfEntries('OneGoodEdgeLimits', numUseEdge)
    do i = 1, numUseEdge
      ierr = PipGetIntegerArray('OneGoodEdgeLimits', ixPcLower, 5, 15)
      ixFrmUseEdge(i) = ixPcLower(1)
      iyFrmUseEdge(i) = ixPcLower(2)
      ixyUseEdge(i) = ixPcLower(3)
      izLowUse(i) = ixPcLower(4)
      izHighUse(i) = ixPcLower(5)
    enddo
    ierr = PipGetLogical('SameEdgeShifts', sameEdgeShifts)
    if ((anyNeg .or. ioptAbs .ne. 0) .and. (shiftEach .or. xcReadIn) &
        .and. .not.undistortOnly) &
        call exitError('you cannot use ShiftPieces or '// &
        'ReadInXcorrs with multiple negatives')
    if (fromEdge .and. xcLegacy .and. .not.undistortOnly) call exitError &
        ('You cannot use both ShiftFromEdges and ShiftFromXcorrs')
    if ((izUseDefLow >= 0 .or. numUseEdge > 0) .and. shiftEach) then
      if (.not. sameEdgeShifts) call &
          exitError('You cannot use good edge limits when shifting pieces')
      if (fromEdge) &
          call exitError('You cannot use ShiftFromEdges with good edge limits')
      xcLegacy = .true.
    endif
  else
    if (anyNeg) then
      print *,'There are multi-negative specifications in list file'
      write(*,'(1x,a,$)') '1 to do initial cross-correlations '// &
          'in overlap zones, 0 not to: '
      read(5,*) ifSloppy
    else
      print *,'Enter the negative of one of the following '// &
          'options to do initial', 'cross-correlations in overlap'// &
          ' zones in combination with the particular option'
      write(*,'(1x,a,/,a,/,a,/,a,/,a,/,a,$)') &
          'Enter 1 to specify division'// &
          ' into negatives to apply to all sections,', &
          '      2 to use edge functions to find a shift for each ' &
          //'frame to align frames', &
          '      3 to use cross-correlation only to find a shift '// &
          'for each frame', &
          '      4 to use only cross-correlation displacements '// &
          'read from a file', &
          '      5 to use best shifts from edge functions and '// &
          'correlations', &
          '      6 to use best shifts from edge functions and '// &
          'displacements read from a file: '
      read(5,*) ioptNeg
      ifSloppy = 0
      if (ioptNeg < 0) ifSloppy = 1
      ioptAbs = abs(ioptNeg)
      shiftEach = ioptAbs >= 2
      fromEdge = ioptAbs == 2
      xcLegacy = ioptAbs == 3 .or. ioptAbs == 4
      xcReadIn = ioptAbs == 4 .or. ioptAbs == 6
    endif
  endif
  !
  if (ioptAbs == 1) then
    if (pipInput) then
      ierr = PipGetTwoIntegers('FramesPerNegativeXandY', &
          numXmissing, numYmissing)
    else
      write(*,'(1x,a,$)') '# of frames per negative in X;'// &
          ' # missing from left-most negative: '
      read(5,*) numXframesPerNeg, numXmissing
      write(*,'(1x,a,$)') '# of frames per negative in Y;'// &
          ' # missing from bottom-most negative: '
      read(5,*) numYframesPerNeg, numYmissing
    endif
    numXneg = (nxPieces + (numXframesPerNeg - 1) + numXmissing) / numXframesPerNeg
    !
    ! derive frame number of each piece and assign negative #
    !
    do ipc = 1, npcList
      ixFrame = (ixPcList(ipc) - minXpiece) / (nxin - nXoverlap)
      iyFrame = (iyPcList(ipc) - minYpiece) / (nyin - nyOverlap)
      ixNeg = (ixFrame + numXmissing) / numXframesPerNeg
      iyNeg = (iyFrame + numYmissing) / numYframesPerNeg
      negList(ipc) = 1 + ixNeg + iyNeg * numXneg
    enddo
    !
    ! now deduce true multi-neg character of each section
    !
    do iz = minZpiece, maxPieceZ
      ineg = iz + 1 - minZpiece
      multiNeg(ineg) = .false.
      listFirst = -100000
      do ipc = 1, npcList
        if (izPcList(ipc) == iz) then
          if (listFirst == -100000) listFirst = negList(ipc)
          multiNeg(ineg) = &
              multiNeg(ineg) .or. (negList(ipc) .ne. listFirst)
        endif
      enddo
      anyNeg = anyNeg .or. multiNeg(ineg)
    enddo
    !
  endif
  !
  plOutFile = ' '
  if (pipInput) then
    ierr = PipGetString('PieceListOutput', plOutFile)
  else
    write(*,'(1x,a,$)') &
        'Name of new piece list file (Return for none): '
    read(5, '(a)') plOutFile
  endif
  outputpl = plOutFile .ne. ' ' .and. ifEdgeFuncOnly == 0
  aliCoordFile = ' '
  ierr = PipGetString('AlignedPieceCoordFile', aliCoordFile)
  !
  ! find out center of transforms
  !
  gxCen = minXpiece + nxTotalPix / 2
  gyCen = minYpiece + nyTotalPix / 2
  if (pipInput) then
    ierr = PipGetTwoFloats('TransformCenterXandY', gxCen, gyCen)
  else
    if (doGxforms) then
      write(*,'(1x,a,/,a,f6.1,a,f6.1,a,$)')  &
          'Enter true center coordinates of the transforms,', &
          '    or / for the default', gxCen, ',', gyCen, &
          ', the center of the image area: '
      read(5,*) gxCen, gyCen
    endif
  endif
  !
  ! Add 0.5 to get the center of rotation to be around the center of image
  gxCen = gxCen + 0.5
  gyCen = gyCen + 0.5
  !
  ! get list of sections desired, set up default as all sections
  !
  do i = 1, numListZ
    izPcTemp(i) = listZ(i)
  enddo
  numZwant = numListZ
  if (pipInput) then
    if (PipGetString('SectionsToDo', fileName) == 0) &
        call parselist(fileName, izPcTemp, numZwant)
  else
    print *,'Enter list of sections to be included in output '// &
        'file (ranges ok)', '   or / to include all sections'
    call rdlist(5, izPcTemp, numZwant)
  endif
  !
  ! copy/pack list, eliminating non-existent sections and duplicates
  ixOut = 0
  do ix = 1, numZwant
    ierr = 0
    do i = 1, numListZ
      if (listZ(i) == izPcTemp(ix)) ierr = 1
    enddo
    do i = 1, ixOut
      if (izWant(i) == izPcTemp(ix)) ierr = 0
    enddo
    if (ierr == 1) then
      ixOut = ixOut + 1
      izWant(ixOut) = izPcTemp(ix)
    endif
  enddo
  numZwant  = ixOut
  if (numZwant == 0) call exitError('Section output list does not '// &
      'include any actual sections')
  nzAllWant = numZwant
  !
  ! Set flag for whether to write out edge displacements
  ! write edge correlations if they are not being read in and if they are
  ! computed in their entirety: but this needs to be modified later
  xcWriteOut = .not.testMode .and. .not.xcReadIn .and. &
      ((ifSloppy == 1 .and. shiftEach) .or. &
      (ifSloppy == 0 .and. shiftEach .and. .not.fromEdge))
  !
  deallocate(ixPcTemp, iyPcTemp, izPcTemp, negTmp, multiTemp, stat = ierr)
  !
  if (undistortOnly) then
    !
    ! If undistorting, copy sizes etc.
    !
    newXframe = nxin
    newXoverlap = nXoverlap
    newXpieces = nxPieces
    newMinXpiece = minXpiece
    newYframe = nyin
    newYoverlap = nyOverlap
    newYpieces = nyPieces
    newMinYpiece = minYpiece
    actionStr = 'undistorted only'
  else
    !
    ! Set up output size
    !
    minXoverlap = 2
    minYoverlap = 2
    newXframe = 100000000
    newYframe = 100000000
    numTrials = 0
32  minXwant = minXpiece
    minYwant = minYpiece
    maxXwant = minXwant + nxTotalPix - 1
    maxYwant = minYwant + nyTotalPix - 1
    if (pipInput) then
      ierr = PipGetLogical('XcorrDebug', xcorrDebug)
      ierr = PipGetLogical('AdjustOrigin', adjustOrigin)
      ierr = PipGetLogical('ExcludeFillFromEdges', limitData)
      ierr = PipGetLogical('NoResizeForFFT', noFFTsizes)
      ierr = PipGetTwoIntegers('StartingAndEndingX', minXwant, maxXwant)
      ierr = PipGetTwoIntegers('StartingAndEndingY', minYwant, maxYwant)
      ierr = PipGetTwoIntegers('MaximumNewSizeXandY', newXframe, newYframe)
      ierr = PipGetTwoIntegers('MinimumOverlapXandY', minXoverlap, minYoverlap)
      ierr = PipGetInteger('BinByFactor', iBinning)
      if (iBinning < 1) call exitError('Binning must be positive')
      if (iBinning > maxBin) call exitError('Binning is too large')
      if (PipGetString('UnsmoothedPatchFile', edgeName) == 0) then
        call dopen(10, edgeName, 'new', 'f')
        izUnsmoothedPatch = 0
      endif
      if (PipGetString('SmoothedPatchFile', edgeName) == 0) then
        call dopen(11, edgeName, 'new', 'f')
        izSmoothedPatch = 0
      endif
    else
      write(*,'(1x,a,/,a,4i6,a,$)') 'Enter Min X, Max X, Min Y, and'// &
          ' Max Y coordinates of desired output section,' &
          , '    or / for whole input section [=' &
          , minXwant, maxXwant, minYwant, maxYwant, ']: '
      read(5,*) minXwant, maxXwant, minYwant, maxYwant
      write(*,'(1x,a,$)') 'Maximum new X and Y frame size, minimum overlap: '
      read(5,*) newXframe, newYframe, minXoverlap, minYoverlap
    endif
    if (numTrials <= 1) then                     !on first 2 trials, enforce min
      minXoverlap = max(2, minXoverlap)        !overlap of 2 so things look
      minYoverlap = max(2, minYoverlap)        !nice in wimp.  After that, let
    endif                                   !the user have it.
    numTrials = numTrials + 1
    !
    ! If no resizing desired, take exactly what is requested in one frame
    if (noFFTsizes .and. maxXwant + 1 - minXwant <= newXframe .and. &
        maxYwant + 1 - minYwant <= newYframe) then
      nxTotalWant = maxXwant + 1 - minXwant
      newXframe = nxTotalWant
      newXtotalPix = nxTotalWant
      newXpieces = 1
      newXoverlap = 2
      nytotwant = maxYwant + 1 - minYwant
      newYframe = nytotwant
      newYtotalPix = nytotwant
      newYpieces = 1
      newYoverlap = 2
    else
      nxTotalWant = 2 * ((maxXwant + 2 - minXwant) / 2)
      nytotwant = 2 * ((maxYwant + 2 - minYwant) / 2)
      call setOverlap(nxTotalWant, minXoverlap, noFFTsizes, newXframe, 2, newXpieces, &
          newXoverlap, newXtotalPix)
      call setOverlap(nytotwant, minYoverlap, noFFTsizes, newYframe, 2, newYpieces, &
          newYoverlap, newYtotalPix)
    endif
    !
    if (.not.outputpl .and. ifEdgeFuncOnly == 0 .and. &
        (newXpieces > 1 .or. newYpieces > 1)) &
        call exitError('You must specify an output piece list file' &
        //' to have more than one output frame')
    if (iBinning > 1 .and. (newXpieces > 1 .or. newYpieces > 1)) &
        call exitError('With binning, output must be into a single frame')
    if (pipInput) print *,'Output file:'
    write(*,115) newXtotalPix, 'X', newXpieces, newXframe, newXoverlap
    write(*,115) newYtotalPix, 'Y', newYpieces, newYframe, newYoverlap
    !
    if (.not.pipInput) then
      write(*,'(1x,a,$)') '1 to revise frame size/overlap: '
      read(5,*) ifRevise
      if (ifRevise .ne. 0) go to 32
    endif
    !
    newMinXpiece = minXwant - (newXtotalPix - nxTotalWant) / 2
    newMinYpiece = minYwant - (newYtotalPix - nytotwant) / 2
    actionStr = 'blended and recut'
  endif
  write(*,'(a,2i7)') 'Starting coordinates of output in X and Y =', newMinXpiece, &
      newMinYpiece
  !
  nxOut = newXframe
  nyOut = newYframe
  dminOut = 1.e10
  dmaxOut = -1.e10
  grandSum = 0.
  call getBinnedSize(nxOut, iBinning, nxBin, ixOffset)
  call getBinnedSize(nyOut, iBinning, nyBin, iyOffset)
  nzBin = 0
  numLinesWrite = nyBin
  hxCen = nxin / 2.
  hyCen = nyin / 2.
  !
  ! get edge indexes for pieces and piece indexes for edges
  ! first build a map in the array of all pieces present
  !
  allocate(mapAllPiece(nxPieces, nyPieces, numSect), stat = ierr)
  if (ierr .ne. 0) call exitError('Allocating array to map pieces')
  do iz = 1, numSect
    do iy = 1, nyPieces
      do ix = 1, nxPieces
        mapAllPiece(ix, iy, iz) = 0
      enddo
    enddo
  enddo
  !
  do ipc = 1, npcList
    mapAllPiece(1 + (ixPcList(ipc) - minXpiece) / (nxin - nXoverlap), &
        1 + (iyPcList(ipc) - minYpiece) / (nyin - nyOverlap), &
        izPcList(ipc) + 1 - minZpiece) = ipc
    do i = 1, 2
      iedgeLower(ipc, i) = 0
      iedgeUpper(ipc, i) = 0
      limDataInd(ipc) = -1
    enddo
  enddo
  !
  ! look at all the edges in turn, add to list if pieces on both sides
  !
  numAlong = nyPieces
  numAcross = nxPieces
  do ixy = 1, 2
    numEd = 0
    do iz = 1, numSect
      do iy = 1, numAlong
        do ix = 2, numAcross
          if (ixy == 1) then
            ipcLower = mapAllPiece(ix - 1, iy, iz)
            ipcHigh = mapAllPiece(ix, iy, iz)
          else
            ipcLower = mapAllPiece(iy, ix - 1, iz)
            ipcHigh = mapAllPiece(iy, ix, iz)
          endif
          if (ipcLower .ne. 0 .and. ipcHigh .ne. 0) then
            numEd = numEd + 1
            ipieceLower(numEd, ixy) = ipcLower
            ipieceUpper(numEd, ixy) = ipcHigh
            iedgeLower(ipcHigh, ixy) = numEd
            iedgeUpper(ipcLower, ixy) = numEd
            edgeDone(numEd, ixy) = .false.
            ifSkipEdge(numEd, ixy) = 0
            edgeDisplaceX(numEd, ixy) = 0.
            edgeDisplaceY(numEd, ixy) = 0.
          endif
        enddo
      enddo
    enddo
    nedge(ixy) = numEd
    numAlong = nxPieces
    numAcross = nyPieces
  enddo
  !
  ! Allocate data depending on number of pieces (limvar)
  limVar = nxPieces * nyPieces
  if (.not. undistortOnly) then
    allocate(bb(2, limVar), ivarPc(limVar), iallVarPc(limVar), &
        ivarGroup(limVar), listCheck(limVar), fpsWork(20 * limVar + limVar / 4 + 4), &
        dxyVar(limVar, 2), rowTmp(limVar * 2), stat = ierr)
    if (ierr .ne. 0) call exitError('Allocating arrays for finding shifts')
    !
    if (testMode) then
      i = 2 * limVar
      allocate(gradXcenLo(i), gradXcenHi(i), gradYcenLo(i), &
          gradYcenHi(i), overXcenLo(i), overXcenHi(i), &
          overYcenLo(i), overYcenHi(i), dxEdge(limEdge, 2), &
          dyEdge(limEdge, 2), dxAdj(limEdge, 2), dyAdj(limEdge, 2), stat = ierr)
      if (ierr .ne. 0) call exitError('Allocating arrays for test mode')
    endif
  endif
  !
  ! Determine size needed for output and correlation arrays
  maxLineLength = newXframe + 32
  maxBsiz = (ifastSiz + maxBin) * maxLineLength
  !
  ! Find binning up to limit that will get padded size down to target
  nbinXcorr = montXCFindBinning(maxXcorrBinning, nxyXcorrTarget, 0, nxyzIn, &
      nOverlap, aspectMax, extraWidth, padFrac, niceFFTlimit(), nxyPadded, nxyBoxed)
  idimc = nxyPadded
  maxBsiz = max(maxBsiz, 2 * nxyBoxed * nbinXcorr**2)
  ! print *,'nbinxcorr, dims', nbinXcorr, idimc, maxbsiz
  allocate(binLine(maxLineLength), brray(maxBsiz), stat = ierr)
  if (ierr .ne. 0) call exitError('Allocating output line arrays')
  !
  ! get edge file name root, parameters if doing a new one
  !
  if (pipInput .and. .not.undistortOnly) then
    ierr = PipGetBoolean('OldEdgeFunctions', ifOldEdge)
    if (PipGetString('RootNameForEdges', rootName) .ne. 0) call &
        exitError('No root name for edge functions specified')
  else if (.not.undistortOnly) then
    write(*,'(1x,a,$)') '0 for new edge files, 1 to read old ones: '
    read(5,*) ifOldEdge
    write(*,'(1x,a,$)') 'Root file name for edge function files: '
    read(5, '(a)') rootName
  endif
  if (ifOldEdge .ne. 0 .and. ifEdgeFuncOnly .ne. 0) call exitError( &
      'You cannot use old edge functions when just computing edge functions')
  !
  ! Find out about parallel mode
  if (pipInput) then
    ierr = PipGetTwoIntegers('ParallelMode', modeParallel, ix)
    if (modeParallel .ne. 0) then
      yChunks = ix .ne. 0
      if (outputpl .or. newXpieces > 1 .or. newYpieces > 1) &
          call exitError('Parallel mode requires output in one piece with no piece list')
      if (undistortOnly .or. testMode .or. xcorrDebug .or. &
          ifEdgeFuncOnly .ne. 0) call exitError('No parallel mode allowed' &
          //'with UndistortOnly, EdgeFunctionsOnly, TestMode, or XcorrDebug')
      if (modeParallel < 0 .and. ifOldEdge == 0) call exitError( &
          'Parallel mode allowed only when using old edge functions')
      if (modeParallel < 0 .and. xcWriteOut) call exitError( &
          'Parallel mode not allowed if writing edge correlation displacements')
      if (modeParallel > 0) then
        !
        ! Figure out lists to output
        if (modeParallel <= 1) call exitError('Target number of '// &
            'chunks should be at least 2')
        !
        ! Set up for breaking into Z chunks and modify for chunks in Y
        ix = numZwant
        ixOut = 1
        if (yChunks) then
          ix = newYtotalPix / iBinning
          ixOut = newMinYpiece
        endif
        numChunks = min(modeParallel, ix)
        numExtra = mod(ix, numChunks)
        do i = 1, numChunks
          numOut = ix / numChunks
          if (i <= numExtra) numOut = numOut + 1
          if (yChunks) then
            iy = ixOut + iBinning * numOut - 1
            if (i == numChunks) iy = newMinYpiece + newYtotalPix - 1
            write(*,'(a,2i9)') 'LineSubsetToDo ', ixOut, iy
            indXlow = (ixOut - newMinYpiece) / iBinning
            indXhigh = (iy - newMinYpiece) / iBinning
            ixOut = ixOut + iBinning * numOut
          else
            write(*,'(a,$)') 'SubsetToDo '
            call wrlist(izWant(ixOut), numOut)
            indXlow = ixOut - 1
            indXhigh = ixOut + numOut - 2
            ixOut = ixOut + numOut
          endif
          if (i == 1) indXlow = -1
          if (i == numChunks) indXhigh = -1
          write(*,'(a,2i9)') 'ChunkBoundary  ', indXlow, indXhigh
        enddo
        !
        ! Output the image size for convenient parsing
        write(*,'(a,2i9)') 'Output image size:', nxBin, nyBin
        !
        ! Or get subset list
      elseif (modeParallel < -1) then
        if (yChunks) then
          if (PipGetTwoIntegers('LineSubsetToDo', lineStart, lineEnd) .ne. &
              0) call exitError('You must enter LineSubsetToDo '// &
              'with this parallel mode')
          if (lineStart < newMinYpiece .or. lineEnd > &
              newMinYpiece + newYtotalPix - 1) call exitError( &
              'The starting or ending line of the subset is out of range')
          if (modeParallel < -2) call exitError('Only direct writing'// &
              'to a single file is allowed with subsets of lines')
          if (mod(lineStart - newMinYpiece, iBinning) .ne. 0 .or. &
              (mod(lineEnd + 1 - lineStart, iBinning) .ne. 0 .and. &
              lineEnd .ne. newMinYpiece + newYtotalPix - 1)) call exitError &
              ('Starting line and number of lines must be a multiple '// &
              'of the binning')
          iyOutOffset = (lineStart - newMinYpiece) / iBinning
          numLinesWrite = (lineEnd + 1 - lineStart) / iBinning
        else
          if (PipGetString('SubsetToDo', fileName) .ne. 0) call exitError( &
              'You must enter SubsetToDo with this parallel mode')
          !
          ! For direct writing, save full want list as all want list
          if (modeParallel == -2) izAllWant(1:numZwant) = izWant(1:numZwant)
          !
          ! In either case, replace the actual want list
          call parselist(fileName, izWant, numZwant)
        endif
        parallelHDF = modeParallel == -2 .and. iiTestIfHDF(outFile) > 0
        ierr = PipGetString('BoundaryInfoFile', boundFile)
        if (parallelHDF .and. boundFile == ' ') call exitError('A boundary info file '// &
            'must be entered for parallel mode -2 if output file type is HDF')
      else
        !
        ! Open output file: set Z size for the file (nzbin is kept as a
        ! running count of sections output otherwise)
        nzBin = numZwant
        parallelHDF = b3dOutputFileType() == 5
      endif
    endif
  endif
  !
  ! Initialize parallel writing if there is a boundary file
  ixy = nxBin
  if (parallelHDF) ixy = -nxBin
  ierr = parWrtInitialize(boundFile, 5, ixy, nyBin, nzAllWant)
  if (ierr .ne. 0) then
    write(*,'(/,a,i3)') 'ERROR: BLENDMONT - Initializing parallel write '// &
        'boundary file, error', ierr
    call exit(1)
  endif
  !
  edgesIncomplete = .false.
  if (ifOldEdge .ne. 0 .and. .not.undistortOnly) then
    !
    ! for old files, open, get edge count and # of grids in X and Y
    !
    do ixy = 1, 2
      lenRecord = 24 / nbytes_recl_item
      edgeName = trim(rootName) //edgeExtension(ixy)
      open(iunEdge(ixy), file = edgeName, status = 'old', &
          form = 'unformatted', access = 'direct', recl = lenRecord, err = 53)
      read(iunEdge(ixy), rec = 1) (numEdgeTmp(i, ixy), i = 1, 5)
      close(iunEdge(ixy))
    enddo
    !
    ! make sure edge counts match and intgrid was consistent
    !
    if (numEdgeTmp(1, 1) .ne. nedge(1) .or. numEdgeTmp(1, 2) .ne. nedge(2)) then
      call convert_longs(numEdgeTmp, 10)
      if (numEdgeTmp(1, 1) .ne. nedge(1) .or. numEdgeTmp(1, 2) .ne. nedge(2)) &
          call exitError('Wrong # of edges in edge function file')
      needByteSwap = 1
    endif
    if (numEdgeTmp(4, 1) .ne. numEdgeTmp(5, 2) .or. &
        numEdgeTmp(5, 1) .ne. numEdgeTmp(4, 2)) call exitError( &
        'Inconsistent grid spacings between edge function files')
    !
    ! set up record size and reopen with right record size
    !
    do ixy = 1, 2
      nxGrid(ixy) = numEdgeTmp(2, ixy)
      nyGrid(ixy) = numEdgeTmp(3, ixy)
      lenRecord = 4 * max(6, 3 * (nxGrid(ixy) * nyGrid(ixy) + 2)) / nbytes_recl_item
      edgeName = trim(rootName) //edgeExtension(ixy)
      open(iunEdge(ixy), file = edgeName, status = 'old', &
          form = 'unformatted', access = 'direct', recl = lenRecord, err = 53)
    enddo
    !
    ! Read edge headers to set up edgedone array
    numZero = 0
    do ixy = 1, 2
      do iedge = 1, nedge(ixy)
        read(iunEdge(ixy), rec = 1 + iedge, err = 52) ixPcLower(1), ixPcLower(2)
        if (needByteSwap .ne. 0) call convert_longs(ixPcLower, 2)
        ! print *,'edge', ixy, iedge, ixpclo(1), ixpclo(2)
        if (ixPcLower(1) >= 0 .and. ixPcLower(2) >= 0) &
            edgeDone(iedge, ixy) = .true.
      enddo
      numZero = numZero + 1
52    continue
    enddo
    !
    ! If either set of edge functions is incomplete, set flag; only set
    ! the interval now if edges are complete
    if (numZero < 2) then
      edgesIncomplete = .true.
    else
      intGrid(1) = numEdgeTmp(4, 1)
      intGrid(2) = numEdgeTmp(5, 1)
    endif
    go to 54
    !
    ! 8/8/03: if there is an error opening old files, just build new ones
    ! This is to allow command files to say use old functions even if
    ! a previous command file might not get run
    !
53  ifOldEdge = 0
    if (ixy == 2) close(iunEdge(1))
    write(*,'(/,a)') 'WARNING: BLENDMONT - Error opening old edge '// &
        'function file; new edge functions will be computed'
    do ixy = 1, 2
      do iedge = 1, nedge(ixy)
        edgeDone(iedge, ixy) = .false.
      enddo
    enddo
  endif

54 if (modeParallel < 0 .and. (ifOldEdge == 0 .or. edgesIncomplete)) &
      call exitError('Parallel mode not allowed if edge functions '// &
      'are incomplete or nonexistent')
  if ((ifOldEdge == 0 .or. edgesIncomplete) .and. .not.undistortOnly) then
    ifDiddle = 0
    ! write(*,'(1x,a,$)') &
    ! '1 to diddle with edge function parameters: '
    ! read(5,*) ifdiddle
    sdCrit = 2.
    devCrit = 2.
    ipolyOrder = 2
    !
    ! make smart defaults for grid parameters
    !
    gridScale = min(8., max(1., max(nxin, nyin) / 512.))
    do ixy = 1, 2
      iboxSiz(ixy) = nint(iboxSiz(ixy) * gridScale)
      indent(ixy) = nint(indent(ixy) * gridScale)
      intGrid(ixy) = nint(intGrid(ixy) * gridScale)
      lastWritten(ixy) = 0
    enddo
    !
    if (pipInput) then
      ierr = PipGetIntegerArray('BoxSizeShortAndLong', iboxSiz, 2, 2)
      ierr = PipGetIntegerArray('IndentShortAndLong', indent, 2, 2)
      ierr = PipGetIntegerArray('GridSpacingShortAndLong', intGrid, 2, 2)
    endif
    ! print *,'box size', (iboxsiz(i), i=1, 2), '  grid', (intgrid(i), i=1, 2)

    if (ifDiddle .ne. 0) then
      write(*,'(1x,a,$)') 'criterion # of sds away from mean'// &
          ' of sd and deviation: '
      read(5,*) sdCrit, devCrit
      write(*,'(1x,a,$)') '# grid positions in short and long'// &
          ' directions to include in regression: '
      read(5,*) (numFit(i), i = 1, 2)
      write(*,'(1x,a,$)') 'order of polynomial: '
      read(5,*) ipolyOrder
      write(*,'(1x,a,$)') 'intervals at which to do regression' &
          //' in short and long directions: '
      read(5,*) (nskipRegress(i), i = 1, 2)
    endif
    !
    ! get edge function characteristics for x and y edges
    !
    ! print *,(iboxsiz(I), indent(i), intgrid(i), i=1, 2)
    do ixy = 1, 2
      call setgridchars(nxyzIn, nOverlap, iboxSiz, indent, intGrid, &
          ixy, 0, 0, 0, 0, nxGrid(ixy), nyGrid(ixy), igridStart, iOffset)
      ! print *,ixy, nxgrid(ixy), nygrid(ixy), nedgetmp(2, ixy), nedgetmp(3, ixy)
      if (edgesIncomplete .and. (nxGrid(ixy) .ne. numEdgeTmp(2, ixy) .or. &
          nyGrid(ixy) .ne. numEdgeTmp(3, ixy))) call exitError('Cannot use'// &
          ' incomplete old edge function file with current parameters')
    enddo
  endif
  if (ifOldEdge == 0  .and. .not.undistortOnly .and. modeParallel <= 0) then
    needByteSwap = 0
    ! open file, write header record
    ! set record length to total bytes divided by system-dependent
    ! number of bytes per item
    !
    do ixy = ixyFuncStart, ixyFuncEnd
      lenRecord = 4 * max(6, 3 * (nxGrid(ixy) * nyGrid(ixy) + 2)) / nbytes_recl_item
      edgeName = trim(rootName) //edgeExtension(ixy)
      ierr = imodBackupFile(edgeName)
      if (ierr .ne. 0) write(6, '(/,a)') ' WARNING: BLENDMONT - error renaming'// &
          ' existing edge function file'
      open(iunEdge(ixy), file = edgeName, status = 'new' &
          , form = 'unformatted', access = 'direct', recl = lenRecord)
      ! write header record
      write(iunEdge(ixy), rec = 1) nedge(ixy), nxGrid(ixy), nyGrid(ixy) &
          , intGrid(ixy), intGrid(3 - ixy)
    enddo
  endif
  !
  ! Do not write ecd file if not all sections are being computed
  if (ifOldEdge == 1 .and. numZwant .ne. numListZ) xcWriteOut = .false.
  !
  ! Read old .ecd file(s)
  if (xcReadIn .and. .not.undistortOnly) then
    ierr = PipGetFloat('BinningForEdgeShifts', ecdBinning)
    iedgeDelX = 0
    iedgeDelY = 0
    if (PipGetTwoIntegers('OverlapForEdgeShifts', iedgeDelX, iedgeDelY) &
        == 0) then
      iedgeDelX = nXoverlap - nint(ecdBinning * iedgeDelX)
      iedgeDelY = nyOverlap - nint(ecdBinning * iedgeDelY)
      print *,'Adjusting by ', iedgeDelX, iedgeDelY
    endif
    edgeName = trim(rootName) //trim(xcorrExtension(0))
    inquire(file = edgeName, exist = exist)
    iy = 4
    if (.not.exist) then
      !
      ! If the file does not exist, look for the two separate files
      ! and put second name in a different variable.  Set unit number to
      ! read one file then the other.
      edgeName = trim(rootName) //trim(xcorrExtension(1))
      inquire(file = edgeName, exist = exist)
      if (exist) then
        edgeName2 = trim(rootName) //trim(xcorrExtension(2))
        inquire(file = edgeName, exist = exist)
      endif
      if (.not.exist) then
        edgeName = trim(rootName) //trim(xcorrExtension(0))
        write(*,'(/,a,a)') 'ERROR: BLENDMONT - Edge correlation file does'// &
            ' not exist: ', trim(edgeName)
        call exit(1)
      endif
      call dopen(5, edgeName2, 'ro', 'f')
      iy = 5
    endif
    call dopen(4, edgeName, 'ro', 'f')
    read(4,*) numEdgeTmp(1, 1), numEdgeTmp(1, 2)
    if (numEdgeTmp(1, 1) .ne. nedge(1) .or. numEdgeTmp(1, 2) .ne. nedge(2)) &
        call exitError('Wrong # of edges in edge correlation file')
    ix = 4

    do ixy = 1, 2
      do i = 1, nedge(ixy)
        read(ix, '(a)') titleStr
        call frefor(titleStr, title, ixOut)
        edgeDisplaceX(i, ixy) = title(1) * ecdBinning
        edgeDisplaceY(i, ixy) = title(2) * ecdBinning
        if (ixy == 1) edgeDisplaceX(i, ixy) = edgeDisplaceX(i, ixy) + iedgeDelX
        if (ixy == 2) edgeDisplaceY(i, ixy) = edgeDisplaceY(i, ixy) + iedgeDelY
        !
        ! Read in skip edge flag and treat it same as when using an
        ! exclusion model
        if (ixOut > 2) then
          if (title(3) .ne. 0) then
            ifSkipEdge(i, ixy) = 2
            if ((edgeDisplaceX(i, ixy) .ne. 0. .or. &
                edgeDisplaceY(i, ixy) .ne. 0.) .and. ifUseAdjusted > 0) then
              ifSkipEdge(i, ixy) = 1
              if (ifUseAdjusted > 1) ifSkipEdge(i, ixy) = 0
            endif
          endif
        endif
      enddo
      ix = iy
    enddo
    close(4)
    close(5)
  endif
  ixgDim = 0
  iygDim = 0
  do ixy = 1, 2
    ixgDim = max(ixgDim, nxGrid(ixy))
    iygDim = max(iygDim, nyGrid(ixy))
  enddo
  !
  ! make default blending width be 80% of overlap up to 50, then
  ! half of overlap above 50
  iblend(1) = max(max(nXoverlap / 2, 50), min(4 * nXoverlap / 5, 50))
  iblend(2) = max(max(nyOverlap / 2, 50), min(4 * nyOverlap / 5, 50))
  if (pipInput) then
    ierr = PipGetIntegerArray('BlendingWidthXandY', iblend, 2, 2)
  else
    write(*,'(1x,a,2i5,a,$)') 'Blending width in X & Y (/ for', &
        iblend(1), iblend(2), '): '
    read(5,*) (iblend(i), i = 1, 2)
  endif
  !
  ! Do other Pip-only options,
  ! Read in mag gradient and distortion field files if specified
  !
  if (pipInput) then
    interpOrder = 3
    ierr = PipGetInteger('InterpolationOrder', interpOrder)
    ierr = PipGetLogical('AdjustedFocus', focusAdjusted)
    if (PipGetString('GradientFile', fileName) == 0) then
      doMagGrad = .true.
      call readMagGradients(fileName, limSect, pixelMagGrad, axisRot, tiltAngles, &
          dmagPerUm, rotPerUm, numMagGrad)
      if (numMagGrad .ne. numListZ) &
          print *,'WARNING: BLENDMONT - # of mag gradients (', &
          numMagGrad, ') does not match # of sections (', numListZ, ')'
      numAngles = numMagGrad
    endif
    !
    ! Look for tilt angles if no mag gradients, then adjust if any
    !
    if (numAngles == 0 .and. PipGetString('TiltFile', fileName) == 0) then
      call read_tilt_file(numAngles, 14, fileName, tiltAngles, limSect)
      if (numAngles .ne. numListZ) &
          print *,'WARNING: BLENDMONT - # of tilt angles (', &
          numAngles, ') does not match # of sections (', numListZ, ')'
    endif
    if (numAngles > 0 .and. PipGetFloat('OffsetTilts', tiltOffset) == 0) then
      do i = 1, numAngles
        tiltAngles(i) = tiltAngles(i) + tiltOffset
      enddo
    endif
    !
    ! If doing added gradients, use a tiltgeometry entry only if no
    ! gradient file
    !
    if (PipGetTwoFloats('AddToGradient', delDmagPerUm, delRotPerUm) == 0) then
      if (doMagGrad) then
        do i = 1, numMagGrad
          dmagPerUm(i) = dmagPerUm(i) + delDmagPerUm
          rotPerUm(i) = rotPerUm(i) + delRotPerUm
        enddo
      else
        if (PipGetThreeFloats('TiltGeometry', pixelMagGrad, axisRot, tiltOffset) &
            .ne. 0) call exitError('-tilt or -gradient must be entered with -add')
        pixelMagGrad = pixelMagGrad * 10.
        numMagGrad = 1
        dmagPerUm(1) = delDmagPerUm
        rotPerUm(1) = delRotPerUm
        doMagGrad = .true.
        if (numAngles == 0) then
          numAngles = 1
          tiltAngles(1) = tiltOffset
        endif
      endif
    endif
    if (doMagGrad) then
      lmField = 200
      maxFields = 16
    endif
    !
    if (PipGetString('DistortionField', fileName) == 0) then
      undistort = .true.
      ierr = readCheckWarpFile(fileName, 1, 1, idfNx, idfNy, ix, idfBinning, &
          pixelIdf, iy, edgeName)
      if (ierr < 0) call exitError(edgeName)

      ierr = getGridParameters(1, nxField, nyField, xFieldStrt, yFieldStrt, &
          xFieldIntrv, yFieldIntrv)
      lmField = max(lmField, nxField, nyField)

      if (PipGetFloat('ImagesAreBinned', binningOfInput) .ne. 0) then
        !
        ! If input binning was not specified object if it is ambiguous
        !
        if (nxin <= idfNx * idfBinning / 2 .and. nyin <= idfNy * idfBinning / 2) &
            call exitError('You must specify binning of images because they '// &
            'are not larger than half the camera size')
      endif
      if (binningOfInput <= 0) call exitError &
          ('Image binning must be a positive number')
    endif
    !
    ! Allocate field arrays and load distortion field
    if (doMagGrad .or. undistort) then
      allocate(distDx(lmField, lmField), distDy(lmField, lmField), &
          fieldDx(lmField, lmField, maxFields), fieldDy(lmField, lmField, maxFields), &
          stat = ierr)
      call memoryError(ierr, 'arrays for distortion fields')
    endif
    if (undistort) then
      if (getWarpGrid(1, nxField, nyField, xFieldStrt, yFieldStrt, xFieldIntrv, &
          yFieldIntrv, distDx, distDy, lmField) .ne. 0) call exitError( &
          'Getting distortion field from warp file')
      !
      ! Adjust grid start and interval and field itself for the
      ! overall binning
      !
      binRatio = idfBinning / binningOfInput
      xFieldStrt = xFieldStrt * binRatio
      yFieldStrt = yFieldStrt * binRatio
      xFieldIntrv = xFieldIntrv * binRatio
      yFieldIntrv = yFieldIntrv * binRatio
      !
      ! if images are not full field, adjust grid start by half the
      ! difference between field and image size
      !
      xFieldStrt = xFieldStrt - (idfNx * binRatio -  nxin) / 2.
      yFieldStrt = yFieldStrt - (idfNy * binRatio -  nyin) / 2.
      !
      ! scale field
      do iy = 1, nyField
        do i = 1, nxField
          distDx(i, iy) = distDx(i, iy) * binRatio
          distDy(i, iy) = distDy(i, iy) * binRatio
        enddo
      enddo
    endif
    print *,nxField, nyField, lmField
    ! print *,xFieldStrt, yfieldStrt, xFieldIntrv, yFieldIntrv
    ! write(*,'(10f7.2)') (distDx(i, 5), distDy(i, 5), i=1, min(nxField, 10))
    !
    ! Handle debug output - open files and set flags
    !
    if (xcorrDebug) then
      if (nxPieces > 1) then
        edgeName = trim(rootName) //'.xdbg'
        call imopen(3, edgeName, 'new')
        ifDumpXY(1) = 0
      endif
      if (nyPieces > 1) then
        edgeName = trim(rootName) //'.ydbg'
        call imopen(4, edgeName, 'new')
        ifDumpXY(2) = 0
      endif
    endif
    !
    ! Get fill treatment; if not entered and very sloppy and distortion,
    ! set for taper
    if (PipGetInteger('TreatFillForXcorr', ifillTreatment) .ne. 0 .and. &
        verySloppy .and. (undistort .or. doMagGrad)) ifillTreatment = 2
    !
    ! Check for model of edges to exclude
    if (PipGetString('SkipEdgeModelFile', fileName) == 0 .and. &
        .not. undistortOnly)  call readExclusionModel(fileName, edgeDisplaceX, &
        edgeDisplaceY, limEdge, ifUseAdjusted, mapAllPiece, nxPieces, nyPieces, &
        minZpiece, numSkip)
  endif
  !
  ! Now that skip flags have been set, write ecd file if
  ! it was read in from two halves
  if (edgeName2 .ne. ' ') call writeEdgeCorrelations()
  !
  doFields = undistort .or. doMagGrad
  if (undistortOnly .and. .not. doFields) call exitError('You must'// &
      'enter -gradient and/or -distort with -justUndistort')
  if (undistortOnly .and. testMode) call exitError( &
      'You cannot enter both -test and -justUndistort')

  !
  ! Set up for warping
  if (doWarp) then
    if (PipGetTwoIntegers('UnalignedStartingXandY', ixUnaliStart, iyUnaliStart) &
        .ne. 0) then
      !
      ! If there is no unaligned start entered and the output area matches the
      ! warp file area, then assume the same start as the current output
      ixUnaliStart = newMinXpiece
      iyUnaliStart = newMinYpiece
      !
      ! Otherwise issue warning if sizes don't match and assume it was centered on
      ! input / full output
      if (nint(warpScale * iwarpNx) .ne. newXtotalPix .or. &
          nint(warpScale * iwarpNy) .ne. newYtotalPix) then
        ixUnaliStart = nint(minXpiece + nxTotalPix / 2. - warpScale * iwarpNx / 2.)
        iyUnaliStart = nint(minYpiece + nyTotalPix / 2. - warpScale * iwarpNy / 2.)
        write(edgeName, '(a,a,a,2i7)') 'WARNING: BLENDMONT - Area being output is ', &
            'different size from unaligned area; you may need to enter ', &
            'UnalignedStartingXandY for warping to work right; assuming starts of', &
            ixUnaliStart, iyUnaliStart
        write(*,'(/,a)') trim(edgeName)
      endif
    endif
    !
    ! Given starting coordinate for warping coordinates, get an offset from that
    ! Just divide the offset by warpScale here, it's never needed otherwise
    warpXoffset = (newMinXpiece + newXtotalPix / 2. - &
        (ixUnaliStart + warpScale * iwarpNx / 2.)) / warpScale
    warpYoffset = (newMinYpiece + newYtotalPix / 2. - &
        (iyUnaliStart + warpScale * iwarpNy / 2.)) / warpScale
    if (setCurrentWarpFile(indWarpFile) .ne. 0) call exitError( &
        'Setting current warp file')
    allocate(numControl(numGxforms), stat = ierr)
    call memoryError(ierr, 'array for number of control points')
    if (findMaxGridSize(warpXoffset, warpXoffset + newXtotalPix / warpScale, &
        warpYoffset, warpYoffset + newYtotalPix / warpScale, numControl, lmWarpX, &
        lmWarpY, edgeName) .ne. 0) call exitError(edgeName)
    allocate(warpDx(lmWarpX, lmWarpY), warpDy(lmWarpX, lmWarpY), stat = ierr)
    call memoryError(ierr, 'arrays for warping grids')
  endif

  call PipDone()
  !
  ! Allocate more arrays now
  deallocate(mapAllPiece)
  allocate(xcray(idimc / 2), xdray(idimc / 2), mapPiece(nxPieces, nyPieces), &
      mapDisjoint(nxPieces, nyPieces), anyDisjoint(nxPieces, nyPieces), &
      limDataLo(limVar, 2, 2), limDataHi(limVar, 2, 2), &
      dxGrBf(ixgDim, iygDim, limEdgBf), dyGrBf(ixgDim, iygDim, limEdgBf), &
      ddenGrBf(ixgDim, iygDim, limEdgBf), dxGrid(ixgDim, iygDim), &
      dyGrid(ixgDim, iygDim), ddenGrid(ixgDim, iygDim), &
      sdGrid(ixgDim, iygDim), stat = ierr)
  call memoryError(ierr, 'correlation or edge buffer arrays')
  ix = 10
  if (numXcorrPeaks > 1 .and. .not. xcLegacy) ix = idimc / 2
  allocate(xeray(ix), stat = ierr)
  call memoryError(ierr, 'correlation array')
  !
  ! Set maximum load; allow one extra slot if doing fields
  ! Limit by the maximum number of pieces that could be needed, which is
  ! the full array, or just one if undistorting
  npixIn = nxin * int(nyin, kind = 8)
  mostNeeded = max(2, nxPieces * nyPieces)
  if (undistortOnly) mostNeeded = 1
  !
  ! First compute a limit based on the minimum memory we would like to use
  if (doFields) then
    maxLoad = min(memMinimum / npixIn - 1, memLim, maxFields, mostNeeded)
    ix = 1
  else
    maxLoad = min(memMinimum / npixIn, memLim, mostNeeded)
    ix = 0
  endif
  !
  ! Increase that to 4 or whatever is needed if it is below that
  ! but if this doesn't fit into the preferred memory limit, then
  ! just go for 2 pieces and make sure that fits too!
  maxLoad = max(maxLoad, min(4, mostNeeded))
  if ((maxLoad + ix) * npixIn > memPreferred) maxLoad = min(2, mostNeeded)
  if ((maxLoad + ix) * npixIn > memMaximum) call exitError( &
      'Images too large for 32-bit array indexes')
  maxSiz = (maxLoad + ix) * npixIn
  allocate(array(maxSiz), stat = ierr)
  if (ierr .ne. 0 .and. maxLoad > min(2, mostNeeded)) then
    maxLoad = min(2, mostNeeded)
    maxSiz = (maxLoad + ix) * npixIn
    allocate(array(maxSiz), stat = ierr)
  endif
  if (ierr .ne. 0) call exitError('Allocating main image array')
  write(*,'(a,i5,a)') 'Allocated', maxSiz / (1024 * 256), &
      ' MB of memory for main image array'
  !
  ! initialize memory allocator and dx, dy lists for looking for near piece
  !
  call clearShuffle()
  call initNearList()
  !
  intGrCopy(1) = intGrid(1)
  intGrCopy(2) = intGrid(2)
  !
  ! All errors checked, exit if setting up parallel
  if (modeParallel > 0) call exit(0)
  !
  ! Now open output files after errors have been checked, unless direct
  ! writing in parallel mode
  if (modeParallel .ne. -2 .and. ifEdgeFuncOnly == 0) then
    !
    ! Do as much of header as possible, shift origin same as in newstack
    call imopen(2, outFile, 'new')
    call iiuTransHeader(2, 1)
    call iiuAltNumExtended(2, 0)
    call iiuAltExtendedType(2, 0, 0)
    call iiuAltMode(2, modeOut)
    call iiuAltSize(2, nxyzBin, nxyzst)
    call iiuRetOrigin(1, xOrigin, yOrigin, zOrigin)
    xOrigin = xOrigin - delta(1) * ixOffset
    yOrigin = yOrigin - delta(2) * iyOffset
    if (adjustOrigin) then
      xOrigin = xOrigin - delta(1) * (newMinXpiece - minXpiece)
      yOrigin = yOrigin - delta(2) * (newMinYpiece - minYpiece)
      zOrigin = zOrigin - delta(3) * (izWant(1) - listZ(1))
    endif
    call iiuAltOrigin(2, xOrigin, yOrigin, zOrigin)
    cell(1) = nxBin * delta(1) * iBinning
    cell(2) = nyBin * delta(2) * iBinning
    !
    ! Finish it and exit if setting up for direct writing
    if (modeParallel == -1) then
      call iiuAltSample(2, nxyzBin)
      cell(3) = nzBin * delta(3)
      call iiuAltCell(2, cell)
    endif
    !
    call b3ddate(dat)
    call time(tim)
    write(titleStr, 90) actionStr, dat, tim
90  format( 'BLENDMONT: Montage pieces ',a, t57, a9, 2x, a8 )
    if (parallelHDF) call iiuWriteDummySecToHDF(2)
    call iiuWriteHeaderStr(2, titleStr, 1, dmin, dmax, dmean)
    if (modeParallel == -1) then
      call iiuClose(2)
      call exit(0)
    endif
  else if (ifEdgeFuncOnly == 0) then
    if (parallelHDF) then
      call parWrtProperties(ixy, ipc, ix)
      if (b3dLockFile(ixy) .ne. 0) &
          call exitError('Could not get lock for opening HDF file')
    endif
    call imopen(2, outFile, 'old')
    call irdhdr(2, ixPcLower, ixPcUpper, ixOut, wll, wlr, wul)
    if (ixPcLower(1) .ne. nxBin .or. ixPcLower(2) .ne. nyBin .or. ixOut .ne. modeOut) &
        call exitError('Existing output file does not have right size or mode')
    if (parallelHDF) call iiuParWrtRecloseHDF(2, 1)
  endif
  !
  if (outputpl) call dopen(3, plOutFile, 'new', 'f')
  fastCum = 0.
  slowCum = 0.
  !
  ! loop on z: do everything within each section for maximum efficiency
  !
  do ilistz = 1, numListZ
    doingEdgeFunc = .true.
    izSect = listZ(ilistz)
    !
    ! test if this section is wanted in output: if not, skip in test mode
    !
    ifWant = 0
    do iwant = 1, numZwant
      if (izWant(iwant) == izSect) ifWant = 1
    enddo
    if (ifWant == 0 .and. (testMode .or. undistortOnly)) cycle
    !
    ! Look up actual section to write if doing direct parallel writes
    if (ifWant .ne. 0 .and. modeParallel == -2 .and. .not. yChunks) then
      do i = 1, nzAllWant
        if (izAllWant(i) == izSect) numOut = i - 1
      enddo
    endif
    !
    write(*,'(a,i5)') ' working on section #', izSect
    call flush(6)
    multng = multiNeg(izSect + 1 - minZpiece)
    xIsLongDim = nxPieces > nyPieces
    !
    ! make a map of pieces in this section and set up index to data limits
    !
    do iyFrame = 1, nyPieces
      do ixFrame = 1, nxPieces
        mapPiece(ixFrame, iyFrame) = 0
      enddo
    enddo
    ipcLower = 0
    do ipc = 1, npcList
      if (izPcList(ipc) == izSect) then
        ixFrame = 1 + (ixPcList(ipc) - minXpiece) / (nxin - nXoverlap)
        iyFrame = 1 + (iyPcList(ipc) - minYpiece) / (nyin - nyOverlap)
        mapPiece(ixFrame, iyFrame) = ipc
        ipcLower = ipcLower + 1
        limDataInd(ipc) = ipcLower
        do i = 1, 2
          limDataLo(ipcLower, i, 1) = -1
          limDataLo(ipcLower, i, 2) = -1
          limDataHi(ipcLower, i, 1) = -1
          limDataHi(ipcLower, i, 2) = -1
        enddo
      endif
    enddo
    !
    if (.not. undistortOnly) then
      !
      ! First get edge functions for the section if they haven't been done
      call findSectionEdgeFunctions()
    endif
    !
    ! now if doing multinegatives, need to solve for h transforms
    !
    if (multng .and. .not. undistortOnly) then
      call findMultinegTransforms()
    endif                                   !end of multineg stuff

    if (.not. undistortOnly) then
      !
      ! initialize edge buffer allocation
      !
      jusEdgCt = 0
      do ixy = 1, 2
        do i = 1, nedge(ixy)
          ibufEdge(i, ixy) = 0
        enddo
      enddo
      do i = 1, limEdgBf
        iedgBfList(i) = -1
        lasEdgUse(i) = 0
      enddo
      !
      ! if this section is not wanted in output, skip out
      !
      if (ifWant == 0) cycle
      !
      ! scan through all edges in this section to determine limits for
      ! when a point is near an edge
      !
      do ixy = ixyFuncStart, ixyFuncEnd
        edgeLoNear(ixy) = 0.
        edgeHiNear(ixy) = nxyzIn(ixy) - 1.
        do iedge = 1, nedge(ixy)
          if (izPcList(ipieceLower(iedge, ixy)) == izSect) then
            jedge = iedge
            if (useEdges) call findEdgeToUse(iedge, ixy, jedge)
            ! print *,'reading header of edge', ixy, jedge, ' for edge', ixy, iedge
            if (needByteSwap == 0) then
              read(iunEdge(ixy), rec = 1 + jedge) nxGridIn, nyGridIn, (igridStart(i), &
                  iOffset(i), i = 1, 2)
            else
              read(iunEdge(ixy), rec = 1 + jedge) nxGridIn, nyGridIn, &
                  (igridStart(i), iOffset(i), i = 1, 2)
              call convert_longs(nxGridIn, 1)
              call convert_longs(nyGridIn, 1)
              call convert_longs(igridStart, 2)
              call convert_longs(iOffset, 2)
            endif
            edgeHiNear(ixy) = min(edgeHiNear(ixy), &
                float(igridStart(ixy)))
            edgeLoNear(ixy) = max(edgeLoNear(ixy), float &
                ((min(nxGridIn, nyGridIn) - 1) * intGrid(1) + iOffset(ixy) + 10))
            ! print *,nxgr, nygr, (igridstr(i), iofset(i), i=1, 2)
          endif
        enddo
        ! print *,ixy, edgelonear(ixy), edgehinear(ixy)
      enddo
      edgesSeparated = edgeHiNear(1) - edgeLoNear(1) > 0.05 * nxin .and. &
          edgeHiNear(2) - edgeLoNear(2) > 0.05 * nyin
    endif
    !
    ! Now analyze for h transforms if need to shift each piece - this sets
    ! multng because it is a flag that hinv exists and needs to be used
    !
    mapDisjoint = 0
    anyDisjoint = .false.
    if (shiftEach .and. .not. undistortOnly .and. ifEdgeFuncOnly .ne. 1 &
        .and. ifEdgeFuncOnly .ne. 2) then
      call getBestPieceShifts()
      !
      ! Analyze for disjoint edges on this section
      ! Look for non-overlap between cross-corner pieces and code by type
      do iyFrame = 1, nyPieces - 1
        do ixFrame = 1, nxPieces - 1
          if (mapPiece(ixFrame, iyFrame) > 0 .and. mapPiece(ixFrame + 1, iyFrame) &
              > 0 .and. mapPiece(ixFrame, iyFrame + 1) > 0 .and. &
              mapPiece(ixFrame + 1, iyFrame + 1) > 0) then
            ix = mapPiece(ixFrame, iyFrame)
            iy = mapPiece(ixFrame + 1, iyFrame + 1)
            if (ixPcList(ix) + hxf(1, 3, ix) + nxin <= &
                ixPcList(iy) + hxf(1, 3, iy)) mapDisjoint(ixFrame, iyFrame) = 1
            if (iyPcList(ix) + hxf(2, 3, ix) + nyin <= &
                iyPcList(iy) + hxf(2, 3, iy)) mapDisjoint(ixFrame, iyFrame) = 3
            ix = mapPiece(ixFrame, iyFrame + 1)
            iy = mapPiece(ixFrame + 1, iyFrame)
            if (ixPcList(ix) + hxf(1, 3, ix) + nxin <= &
                ixPcList(iy) + hxf(1, 3, iy)) mapDisjoint(ixFrame, iyFrame) = 2
            if (iyPcList(iy) + hxf(2, 3, iy) + nyin <= &
                iyPcList(ix) + hxf(2, 3, ix)) mapDisjoint(ixFrame, iyFrame) = 4
            if (mapDisjoint(ixFrame, iyFrame) .ne. 0) &
                anyDisjoint(ixFrame:ixFrame + 1, iyFrame:iyFrame + 1) = .true.
            ! if (mapDisjoint(ixfrm, iyfrm) .ne. 0) print *,'disjoint', ixfrm, iyfrm
          endif
        enddo
      enddo
    endif
    !
    ! if doing g transforms, get inverse and recenter it from center
    ! of output image to corner of image
    !
    if (doGxforms) then
      if (skipXforms) then
        indGxf = izSect + 1 - minZpiece
      else
        do ilis = 1, numListZ
          if (listZ(ilis) == izSect) indGxf = ilis
        enddo
      endif
      call xfcopy(gxf(1, 1, indGxf), gxfTemp)
      gxfTemp(1, 3) = gxfTemp(1, 3) + (1. -gxfTemp(1, 1)) * gxCen - gxfTemp(1, 2) * gyCen
      gxfTemp(2, 3) = gxfTemp(2, 3) + (1. -gxfTemp(2, 2)) * gyCen - gxfTemp(2, 1) * gxCen
      call xfinvert(gxfTemp, ginv)
    endif
    !
    ! Adjust flag for shuffler to know whether to undistort, and clear
    ! out the undistorted pieces in memory
    !
    doingEdgeFunc = .false.
    if (doFields) call clearShuffle()
    if (testMode .or. ifEdgeFuncOnly .ne. 0) cycle ! To end of section loop
    !
    ! If warping, find out if there is warping on this section and get grid
    secHasWarp = .false.
    if (doWarp) then
      secHasWarp = numControl(indGxf) > 2
    endif
    if (secHasWarp) then
      !
      ! Get the grid without adjustment of coordinates for any offset, then adjust
      ! the grid start for the old starting X coordinate
      if (getSizeAdjustedGrid(indGxf, newXtotalPix / warpScale, newYtotalPix / warpScale, &
          warpXoffset, warpYoffset, 0, warpScale, 1, nxWarp, nyWarp, xWarpStrt, &
          yWarpStrt, xWarpIntrv, yWarpIntrv, warpDx, warpDy, lmWarpX, lmWarpY, &
          edgeName) .ne. 0) call exitError(edgeName)
      xWarpStrt = xWarpStrt + ixUnaliStart
      yWarpStrt = yWarpStrt + iyUnaliStart
      ! print *,nxWarp, nyWarp, xWarpStrt, yWarpStrt, xWarpIntrv, yWarpIntrv
      ! write(*,'(10f7.2)') ((warpDx(ix, iy), warpDy(ix, iy), ix=1, 10), iy=3, 4)
    endif
    !
    ! if floating, need to get current input min and max
    ! To pad edges properly, need current mean: put it into dmean
    !
    call  crossValue(xIsLongDim, nxPieces, nyPieces, numShort, numLong)
    curInMax = -1.e10
    curInMin = 1.e10
    curSum = 0.
    realNsum = 0.
    do ilong = numLong, 1, -1
      do ishort = numShort, 1, -1
        call crossValue(xIsLongDim, ishort, ilong, ixFrame, iyFrame)
        if (mapPiece(ixFrame, iyFrame) > 0) then
          call shuffler(mapPiece(ixFrame, iyFrame), indArray)
          do iy = 1, nyin
            tsum = 0.
            do i = indArray + (iy - 1) * nxin, indArray + iy * nxin - 1
              curInMin = min(curInMin, array(i))
              curInMax = max(curInMax, array(i))
              tsum = tsum + array(i)
            enddo
            curSum = curSum + tsum
            realNsum = realNsum + nxin
          enddo
        endif
      enddo
    enddo
    dmean = curSum / realNsum
    dfill = dmean
    if (useFill) dfill = fillVal
    if (ifFloat == 0) then
      curInMin = dfltInMin
      curInMax = dfltInMax
    endif
    !
    ! now get output scaling factor and additive factor
    !
    pixelScale = (outMax - outMin) / max(1., curInMax - curInMin)
    pixelAdd = outMin - pixelScale * curInMin
    !
    ! look through memory list and renumber them with priorities
    ! backwards from the first needed piece
    !
    newUseCount = juseCount - 1
    do ilong = 1, numLong
      do ishort = 1, numShort
        call crossValue(xIsLongDim, ishort, ilong, ixFrame, iyFrame)
        if (mapPiece(ixFrame, iyFrame) > 0) then
          do i = 1, maxLoad
            if (izMemList(i) == mapPiece(ixFrame, iyFrame)) then
              lastUsed(i) = newUseCount
              newUseCount = newUseCount-1
            endif
          enddo
        endif
      enddo
    enddo
    !
    ! UNDISTORTING ONLY: loop on all frames in section, in order as
    ! they were in input file
    !
    if (undistortOnly) then
      call clearShuffle()
      doingEdgeFunc = .true.
      do ipc = 1, npcList
        if (izPcList(ipc) == izSect) then
          call shuffler(ipc, indArray)
          tsum = 0.
          do i = 1, npixIn
            val = pixelScale * array(i + indArray - 1) + pixelAdd
            array(i + indArray - 1) = val
            dminOut = min(dminOut, val)
            dmaxOut = max(dmaxOut, val)
            tsum = tsum + val
          enddo
          grandSum = grandSum + tsum
          nzBin = nzBin + 1
          call iiuWriteSection(2, array(indArray))
          !
          if (outputpl) write(3, '(2i9,i7)') newPcXlowLeft, newPcYlowLeft, izSect
        endif
      enddo
      cycle                                 ! To end of section loop
    endif
    !
    ! GET THE PIXEL OUT
    ! -  loop on output frames; within each frame loop on little boxes
    !
    call  crossValue(xIsLongDim, newXpieces, newYpieces, numShort, numLong)
    !
    do ilong = 1, numLong
      do ishort = 1, numShort
        call crossValue(xIsLongDim, ishort, ilong, ixOut, iyOut)
        ! write(*,'(a,2i4)') ' composing frame at', ixout, iyout
        newPcYlowLeft = newMinYpiece + (iyOut - 1) * (newYframe - newYoverlap)
        newPcXlowLeft = newMinXpiece + (ixOut - 1) * (newXframe - newXoverlap)
        anyPixels = iBinning > 1 .or. newXpieces * newYpieces == 1
        anyLinesOut = .false.
        tsum = 0.
        lineOffset = iyOffset
        linesBuffered = 0
        iBufferBase = 0
        !
        ! do fast little boxes
        !
        nxFast = (newXframe + (ifastSiz - 1)) / ifastSiz
        nyFast = (newYframe + (ifastSiz - 1)) / ifastSiz
        if (yChunks) nyFast = (lineEnd - lineStart + ifastSiz) / ifastSiz
        !
        ! loop on boxes, get lower & upper limits in each box
        !
        do iyFast = 1, nyFast
          indYlow = newPcYlowLeft + (iyFast - 1) * ifastSiz
          indYhigh = min(indYlow + ifastSiz, newPcYlowLeft + newYframe) - 1
          if (yChunks) then
            indYlow = lineStart + (iyFast - 1) * ifastSiz
            indYhigh = min(indYlow + ifastSiz - 1, lineEnd)
          endif
          numLinesOut = indYhigh + 1 - indYlow
          !
          ! fill array with dfill
          !
          do i = 1, nxOut * numLinesOut
            brray(i + iBufferBase) = dfill
          enddo
          !
          do ixFast = 1, nxFast
            indXlow = newPcXlowLeft + (ixFast - 1) * ifastSiz
            indXhigh = min(indXlow + ifastSiz, newPcXlowLeft + newXframe) - 1
            !
            ! check # of edges, and prime piece number, for each corner
            !
            doFast = .true.
            inFrame = .false.
            inOnePiece = 0
            samePieces = .true.
            debug = (ixDebug >= indXlow .and. ixDebug <= indXhigh) .or. &
                (iyDebug >= indYlow .and. iyDebug <= indYhigh)
            do indy = indYlow, indYhigh, max(1, indYhigh - indYlow)
              do indx = indXlow, indXhigh, max(1, indXhigh - indXlow)
                call countEdges(indx, indy, xg, ysrc, useEdges)
                if (numPieces > 0) then
                  if (inOnePiece == 0) then
                    inOnePiece = inPiece(1)
                    numFirst = numPieces
                    ipFirst(1:4) = inPiece(1:4)
                  endif
                  doFast = doFast .and. numPieces == 1 .and. &
                      inPiece(1) == inOnePiece
                  samePieces = samePieces .and. numPieces == numFirst
                  do i = 1, numPieces
                    inFrame = inFrame .or. (xInPiece(i) >= 0. .and. &
                        xInPiece(i) <= nxin - 1. .and. yInPiece(i) &
                        >= 0. .and. yInPiece(i) <= nyin - 1.)
                    samePieces = samePieces .and. ipFirst(i) == inPiece(i)
                  enddo
                else
                  samePieces = .false.
                endif
              enddo
            enddo
            if (debug) print *,indYlow, indYhigh, numPieces, xInPiece(1), yInPiece(1)
            !
            ! ALL ON ONE PIECE: do a fast transform of whole box
            !
            wallStart = wallTime()
            if (doFast .and. inFrame) then
              if (debug) print *,  'fast box', indXlow, indXhigh, indYlow, indYhigh
              call xfunit(fastf, 1.)         !start with unit xform
              !
              ! if doing g xforms, put operations into xform that will
              ! perform inverse of xform
              !
              if (doGxforms) then
                call xfcopy(ginv, fastf)
              endif
              !
              ! shift coordinates down to be within piece
              !
              fastf(1, 3) = fastf(1, 3) - ixPcList(inOnePiece)
              fastf(2, 3) = fastf(2, 3) - iyPcList(inOnePiece)
              !
              ! if doing h's, implement the h inverse
              !
              if (multng) then
                call xfmult(fastf, hinv(1, 1, inOnePiece), fastTemp)
                call xfcopy(fastTemp, fastf)
              endif
              !
              ! now add 1 to get array index
              !
              fastf(1, 3) = fastf(1, 3) + 1.
              fastf(2, 3) = fastf(2, 3) + 1.
              !
              call shuffler(inOnePiece, indArray)
              call fastInterp(brray(iBufferBase + 1), nxOut, numLinesOut, &
                  array(indArray), nxin, nyin, indXlow, indXhigh, indYlow, &
                  indYhigh, newPcXlowLeft, fastf, fastf(1, 3), &
                  fastf(2, 3), inOnePiece)
              anyPixels = .true.
              fastCum = fastCum + wallTime() - wallStart
              !
            elseif (inFrame) then
              !
              ! in or near an edge: loop on each pixel
              !
              numIter = 10
              lineBase = iBufferBase + 1 - newPcXlowLeft
              do indy = indYlow, indYhigh
                do indx = indXlow, indXhigh
                  ! Set debug .or. or .and. here
                  debug = indx == ixDebug .or. indy == iyDebug
                  if (samePieces) then
                    !
                    ! If it is all the same pieces in this box, then update
                    ! the positions in the pieces
                    xg = indx
                    ysrc = indy
                    if (secHasWarp) then
                      call interpolateGrid(xg - 0.5, ysrc - 0.5, warpDx, warpDy, &
                          lmWarpX, nxWarp, nyWarp, xWarpStrt, yWarpStrt, xWarpIntrv, &
                          yWarpIntrv, xg, ysrc)
                      xg = xg + indx
                      ysrc = ysrc + indy
                    endif
                    if (doGxforms) then
                      xtmp = ginv(1, 1) * xg + ginv(1, 2) * ysrc + ginv(1, 3)
                      ysrc = ginv(2, 1) * xg + ginv(2, 2) * ysrc + ginv(2, 3)
                      xg = xtmp
                    endif
                    do i = 1, numPieces
                      call positionInPiece(xg, ysrc, inPiece(i), xInPiece(i), &
                          yInPiece(i))
                    enddo
                  else
                    call countEdges(indx, indy, xg, ysrc, useEdges)
                  endif
                  if (debug) write(*,'(2i6,i7,i2,a,i3,a,5i5)') indx, indy, &
                      numEdges(1), numEdges(2), ' edges', numPieces, &
                      ' pieces', (inPiece(i), i = 1, numPieces)
                  !
                  ! load the edges and compute edge fractions
                  !
                  call computeEdgeFractions()
                  !
                  ! get indices of pieces and edges and the weighting
                  ! of each piece: for now, numbers
                  call getPieceIndicesAndWeighting(useEdges)
                  !
                  ! NOW SORT OUT THE CASES OF 1, 2, 3 or 4 PIECES

                  call getPixelFromPieces()
                  !
                  ! stick limited pixval into array and mark as output
                  !
                  brray(lineBase + indx) = &
                      max(curInMin, min(curInMax, pixVal))
                  anyPixels = .true.
                enddo
                lineBase = lineBase + nxOut
              enddo
              slowCum = slowCum + wallTime() - wallStart
            endif
          enddo
          !
          ! if any pixels have been present in this frame, write line out
          !
          if (modeParallel .ne. -2 .or. yChunks) numOut = nzBin
          if (anyPixels) then
            do i = 1, nxOut * numLinesOut
              brray(i + iBufferBase) = pixelScale * brray(i + iBufferBase) + pixelAdd
            enddo
            !
            ! Set line position based on binned pixels output
            ! Set lines to write based on what is in buffer already
            ! Get new number of lines left in buffer; if at top of
            ! frame, increment lines to write if there are any lines left
            !
            ilineOut = ((iyFast - 1) * ifastSiz - iyOffset) / iBinning
            call parWrtPosn(2, numOut, ilineOut + iyOutOffset)
            ilineOut = linesBuffered + numLinesOut
            nyWrite = (ilineOut - lineOffset) / iBinning
            linesBuffered = mod(ilineOut - lineOffset, iBinning)
            if (indYhigh == newPcYlowLeft + newYframe - 1 .and. &
                linesBuffered > 0) nyWrite = nyWrite + 1
            !
            ! write data
            !
            call iwrBinned(2, brray, binLine, nxOut, nxBin, ixOffset, &
                ilineOut, nyWrite, lineOffset, iBinning, dminOut, &
                dmaxOut, tsum)
            !
            ! Set Y offset to zero after first time, set base and
            ! move remaining lines to bottom
            !
            lineOffset = 0
            ifill = (ilineOut - linesBuffered) * nxOut
            iBufferBase = nxOut * linesBuffered
            do i = 1, iBufferBase
              brray(i) = brray(i + ifill)
            enddo
            !
            ! if this is the first time anything is written, and it
            ! wasn't the first set of lines, then need to go back and
            ! fill the lower part of frame with mean values
            !
            if (.not.anyLinesOut .and. iyFast > 1) then
              val = dfill * pixelScale + pixelAdd
              brray(1:nxOut) = val
              call parWrtPosn(2, numOut, iyOutOffset)
              do ifill = 1, (iyFast - 1) * ifastSiz
                call parWrtLin(2, brray)
              enddo
              tsum = tsum + val * nxOut * (iyFast - 1) * ifastSiz
            endif
            anyLinesOut = .true.
          endif
        enddo
        !
        ! if any pixels present, write piece coordinates
        !
        if (anyPixels) then
          grandSum = grandSum + tsum
          ! write(*,'(a,i5)') ' wrote new frame #', nzbin
          nzBin = nzBin + 1
          !
          if (outputpl) write(3, '(2i9,i7)') newPcXlowLeft, newPcYlowLeft, izSect
        endif
        !
      enddo                                 ! Loop on frames - short dim
    enddo                                   ! Loop on frames - long dim
  enddo                                     ! Loop on sections
  !
  close(3)
  ! write(*,'(a,2f12.6)') 'fast box and single pixel times:', fastcum, slowcum
  if (ifEdgeFuncOnly == 0 .and. .not. testMode) then
    !
    ! If direct parallel, output stats
    pixelTot = (float(nxBin) * numLinesWrite) * nzBin
    tmean = grandSum / pixelTot
    if (modeParallel == -2) then
      write(*,'(a,3g15.7,f15.0)') 'Min, max, mean, # pixels=', dminOut, &
          dmaxOut, tmean, pixelTot
    else
      !
      ! otherwise finalize the header
      !
      call iiuAltSize(2, nxyzBin, nxyzst)
      call iiuAltSample(2, nxyzBin)
      cell(3) = nzBin * delta(3)
      call iiuAltCell(2, cell)
      call iiuWriteHeader(2, title, -1, dminOut, dmaxOut, tmean)
    endif
    if (parallelHDF) then
      if (iiuParWrtFlushBuffers(2) .ne. 0) &
          call exitError("Finishing writing to output HDF file")
      call parWrtClose()
    endif
    call iiuClose(2)
  endif
  if (undistortOnly) call exit(0)
  !
  ! write edge correlations
  !
  if (xcWriteOut) call writeEdgeCorrelations()
  !
  ! Write aligned piece coordinates
  if (aliCoordFile .ne. ' ') then
    call dopen(14, aliCoordFile, 'new', 'f')
    do ipc = 1, npcList
      write(14, '(2i9,i6)') ixPcList(ipc) + nint(hxf(1, 3, ipc)),  iyPcList(ipc) +  &
          nint(hxf(2, 3, ipc)), izPcList(ipc)
    enddo
    close(14)
  endif
      
  !
  ! rewrite header for new edge functions so that they have later date
  ! than the edge correlations; close files
  !
  do ixy = ixyFuncStart, ixyFuncEnd
    if (ifOldEdge == 0) write(iunEdge(ixy), rec = 1) nedge(ixy), &
        nxGrid(ixy), nyGrid(ixy) , intGrid(ixy), intGrid(3 - ixy)
    close(iunEdge(ixy))
  enddo
  do ixy = 1, 2
    if (ifDumpXY(ixy) > 0) then
      call iiuWriteHeader(2 + ixy, title, -1, 0., 255., 128.)
      call imclose(2 + ixy)
    endif
  enddo
  if (izUnsmoothedPatch >= 0) close(10)
  if (izSmoothedPatch >= 0) close(11)
  !
  call exit(0)

CONTAINS

  ! Write the full edge correlation file, or the X or Y component only
  !
  subroutine writeEdgeCorrelations()
    edgeName = trim(rootName) //trim(xcorrExtension(mod(ifEdgeFuncOnly, 3)))
    call dopen(4, edgeName, 'new', 'f')
    if (ixyFuncStart == 1) write(4, '(2i7)') nedge(1), nedge(2)
    do ixy = ixyFuncStart, ixyFuncEnd
      if (numSkip == 0 .and. nedge(ixy) > 0) then
        write(4, '(f9.3,f10.3)') (edgeDisplaceX(i, ixy), edgeDisplaceY(i, ixy), &
            i = 1, nedge(ixy))
      elseif (nedge(ixy) > 0) then
        write(4, '(f9.3,f10.3,i4)') (edgeDisplaceX(i, ixy), edgeDisplaceY(i, ixy), &
            ifSkipEdge(i, ixy), i = 1, nedge(ixy))
      endif
    enddo
    close(4)
    return
  end subroutine writeEdgeCorrelations


  ! GETS EDGE FUNCTIONS FOR A SECTION IF THEY HAVEN'T BEEN DONE YET
  !
  subroutine findSectionEdgeFunctions()
    !
    ! Check if any edges are to be substituted on this section
    ! Do so if any use values are different from original and none are 0
    useEdges = .false.
    if (numUseEdge > 0 .or. izUseDefLow >= 0) then
      numZero = 0
      do ixy = 1, 2
        do iedge = 1, nedge(ixy)
          if (izSect == izPcList(ipieceLower(iedge, ixy))) then
            call findEdgeToUse(iedge, ixy, jedge)
            if (jedge == 0) numZero = numZero + 1
            if (jedge .ne. iedge) useEdges = .true.
          endif
        enddo
      enddo
      if (numZero > 0) useEdges = .false.
      ! print *,'numzero, useedges', numZero, useEdges
    endif
    !
    ! loop on short then long direction
    !
    do iedgeDir = 1, 2
      call  crossValue(xIsLongDim, iedgeDir, 3 - iedgeDir, ixy, iyx)
      if (iyx == ifEdgeFuncOnly) cycle
      !
      ! loop on all edges of that type with pieces in section that are
      ! not done yet
      !
      do iedge = 1, nedge(ixy)
        if (izSect == izPcList(ipieceLower(iedge, ixy))) then
          jedge = iedge
          if (useEdges) call findEdgeToUse(iedge, ixy, jedge)
          ! print *,'checking edge', ixy, jedge, ' for edge', ixy, iedge
          if (.not.edgeDone(jedge, ixy)) then
            !
            ! do cross-correlation if the sloppy flag is set and either
            ! we are shifting each piece or the pieces are't on same neg
            !
            ! print *,'Doing edge ', ixy, jedge
            doCross = ifSloppy .ne. 0 .and. (shiftEach .or. ( &
                anyNeg .and. negList(ipieceLower(jedge, ixy)) &
                .ne. negList(ipieceUpper(jedge, ixy))))
            call doEdge(jedge, ixy, edgeDone, sdCrit, devCrit, &
                numFit, ipolyOrder, nskipRegress, doCross, xcReadIn, xcLegacy, &
                edgeDisplaceX, edgeDisplaceY, limEdge)
            !
            ! after each one, check memory list to see if there's any
            ! pieces with undone lower edge in orthogonal direction
            !
            if (.not.useEdges .and. ixy .ne. ifEdgeFuncOnly) then
              do imem = 1, maxLoad
                ipc = izMemList(imem)
                if (ipc > 0) then
                  jedge = iedgeLower(ipc, iyx)
                  if (jedge > 0) then
                    if (.not.edgeDone(jedge, iyx)) then
                      ! print *,'Doing edge ', iyx, jedge
                      doCross = ifSloppy .ne. 0 .and. (shiftEach .or. ( &
                          anyNeg .and. negList(ipieceLower(jedge, iyx)) &
                          .ne. negList(ipieceUpper(jedge, iyx))))
                      call doEdge(jedge, iyx, edgeDone, sdCrit, devCrit, &
                          numFit, ipolyOrder, nskipRegress, doCross, xcReadIn, &
                          xcLegacy, edgeDisplaceX, edgeDisplaceY, limEdge)
                    endif
                  endif
                endif
              enddo
            endif
          endif
        endif
      enddo
    enddo
    return
  end subroutine findSectionEdgeFunctions


  ! FINDS H TRANSFORMS IF MULTIPLE NEGATIVES (UNTESTED!)
  !
  subroutine findMultinegTransforms()
    real*4 sind, cosd
    !
    ! first make index of negatives and find coordinates of each
    !
    numNegatives = 0
    do i = 1, LIMNEG
      maxNegX(i) = -1000000
      minNegX(i) = 1000000
      maxNegY(i) = -1000000
      minNegY(i) = 1000000
    enddo
    do ipc = 1, npcList
      if (izPcList(ipc) == izSect) then
        iwhich = 0
        do i = 1, numNegatives
          if (negList(ipc) == negIndex(i)) iwhich = i
        enddo
        if (iwhich == 0) then
          numNegatives = numNegatives + 1
          negIndex(numNegatives) = negList(ipc)
          iwhich = numNegatives
        endif
        ! point max coordinate 1 past end: it'll be easier
        minNegX(iwhich) = min(minNegX(iwhich), ixPcList(ipc))
        maxNegX(iwhich) = max(maxNegX(iwhich), ixPcList(ipc) + nxin)
        minNegY(iwhich) = min(minNegY(iwhich), iyPcList(ipc))
        maxNegY(iwhich) = max(maxNegY(iwhich), iyPcList(ipc) + nyin)
      endif
    enddo
    !
    ! look at all edges, make tables of joints between negs
    !
    do ixy = 1, 2
      numJoint(ixy) = 0
      do i = 1, numNegatives
        jointUpper(i, ixy) = 0
        jointLower(i, ixy) = 0
      enddo
      !
      do iedge = 1, nedge(ixy)
        negLow = negList(ipieceLower(iedge, ixy))
        negUp = negList(ipieceUpper(iedge, ixy))
        if (izPcList(ipieceLower(iedge, ixy)) == izSect .and. &
            negLow .ne. negUp) then
          ! convert neg #'s to neg indexes
          do i = 1, numNegatives
            if (negIndex(i) == negLow) indLow = i
            if (negIndex(i) == negUp) indUpper = i
          enddo
          ! see if joint already on list
          joint = 0
          do j = 1, numJoint(ixy)
            if (negLower(j, ixy) == indLow .and. &
                negUpper(j, ixy) == indUpper) joint = j
          enddo
          ! if not, add to list
          if (joint == 0) then
            numJoint(ixy) = numJoint(ixy) + 1
            joint = numJoint(ixy)
            negLower(joint, ixy) = indLow
            negUpper(joint, ixy) = indUpper
            ! point the negatives to the joint
            jointLower(indLow, ixy) = joint
            jointUpper(indUpper, ixy) = joint
            numEdgeOnJoint(joint) = 0
          endif
          ! add edge to list of ones on that joint
          numEdgeOnJoint(joint) = numEdgeOnJoint(joint) + 1
          listOnJoint(numEdgeOnJoint(joint), joint) = iedge
        endif
      enddo
      !
      ! now process all edges on each joint to get a rotrans
      !
      do joint = 1, numJoint(ixy)
        do ied = 1, numEdgeOnJoint(joint)
          iedge = listOnJoint(ied, joint)
          call readEdgeFunc(iedge, ixy, ied)

          ipcLower = ipieceLower(iedge, ixy)
          ixPcLower(ied) = ixGrdStBf(ied) + ixPcList(ipcLower)
          iyPcLower(ied) = iyGrdStBf(ied) + iyPcList(ipcLower)
          ipcUpper = ipieceUpper(iedge, ixy)
          ixPcUpper(ied) = ixOfsBf(ied) + ixPcList(ipcUpper)
          iyPcUpper(ied) = iyOfsBf(ied) + iyPcList(ipcUpper)
        enddo
        call joint_to_rotrans(dxGrBf, dyGrBf, ixgDim, iygDim, &
            nxGrBf, nyGrBf, intGrid(ixy), intGrid(3 - ixy) &
            , ixPcLower, iyPcLower, ixPcUpper, iyPcUpper, numEdgeOnJoint(joint) &
            , r(1, joint, ixy))
        write(*,'(1x,a,2i4,a,f6.2,a,2f6.1)') &
            char(ixy + ichar('W')) //' joint, negatives' &
            , negList(ipcLower), negList(ipcUpper), '   theta =' &
            , r(1, joint, ixy), '   dx, dy =' &
            , r(2, joint, ixy), r(3, joint, ixy)
      enddo
    enddo
    !
    ! Resolve the edges into rotrans centered on each negative
    ! Initially set the rotrans for each negative to null
    !
    maxSwing = 0
    do i = 1, numNegatives
      hCum(1, i) = 0.
      hCum(2, i) = 0.
      hCum(3, i) = 0.
      hCum(4, i) = 0.5 * (maxNegX(i) + minNegX(i))
      hCum(5, i) = 0.5 * (maxNegY(i) + minNegY(i))
      maxSwing = max(maxSwing, (maxNegX(i) - minNegX(i)) / 2, &
          (maxNegY(i) - minNegY(i)) / 2)
    enddo
    !
    ! start loop - who knows how long this will take
    !
    ishift = 1
    errAdjust = 1.e10
    numShiftNeg = 100
    errLim = 0.1 * numNegatives
    do while(ishift <= numShiftNeg .and. errAdjust > errLim)
      !
      ! set sum of adjustment h's to null
      !
      errAdjust = 0.
      do i = 1, numNegatives
        call lincom_rotrans(hCum(1, i), 0., hCum(1, i), 0., hAdj(1, i))
        numInHAdj(i) = 0
      enddo
      !
      ! loop on joints, adding up net adjustments still needed
      !
      do ixy = 1, 2
        do joint = 1, numJoint(ixy)
          !
          ! get net rotrans still needed at this joint: subtract upper
          ! h and add lower h to joint rotrans
          !
          negUp = negUpper(joint, ixy)
          negLow = negLower(joint, ixy)
          call lincom_rotrans(hCum(1, negUp), -1., r(1, joint, ixy), -1., &
              rotxNet)
          call lincom_rotrans(hCum(1, negLow), 1., rotxNet, 1., rotxNet)
          !
          ! now add half of the net needed to the upper adjustment
          ! subtract half from the lower adjustment
          !
          call lincom_rotrans(rotxNet, 0.5, hAdj(1, negUp), 1. &
              , hAdj(1, negUp))
          numInHAdj(negUp) = numInHAdj(negUp) + 1
          call lincom_rotrans(rotxNet, -0.5, hAdj(1, negLow), 1. &
              , hAdj(1, negLow))
          numInHAdj(negLow) = numInHAdj(negLow) + 1
        enddo
      enddo
      !
      ! get the average adjustment, adjust hcum with it
      !
      dxSum = 0.
      dySum = 0.
      do i = 1, numNegatives
        errAdjust = errAdjust + (abs(hAdj(2, i)) + abs(hAdj(3, i)) + &
            abs(hAdj(1, i)) * maxSwing) / numInHAdj(i)
        call lincom_rotrans(hAdj(1, i), 1. / numInHAdj(i), hCum(1, i), 1. &
            , hCum(1, i))
        dxSum = dxSum + hCum(2, i)
        dySum = dySum + hCum(3, i)
      enddo
      !
    enddo                                     !end of cycle
    !
    ! shift all dx and dy to have a mean of zero
    !
    do i = 1, numNegatives
      hCum(1, i) = hCum(1, i) - dxSum / numNegatives
      hCum(2, i) = hCum(2, i) - dySum / numNegatives
    enddo
    !
    ! compute the h function (and hinv) centered on corner of frame
    !
    do ipc = 1, npcList
      if (izPcList(ipc) == izSect) then
        do i = 1, numNegatives
          if (negIndex(i) == negList(ipc)) then
            call recen_rotrans(hCum(1, i), float(ixPcList(ipc)) &
                , float(iyPcList(ipc)), hFrame)
            hxf(1, 1, ipc) = cosd(hFrame(1))
            hxf(2, 1, ipc) = sind(hFrame(1))
            hxf(1, 2, ipc) = -hxf(2, 1, ipc)
            hxf(2, 2, ipc) = hxf(1, 1, ipc)
            hxf(1, 3, ipc) = hFrame(2)
            hxf(2, 3, ipc) = hFrame(3)
            call xfinvert(hxf(1, 1, ipc), hinv(1, 1, ipc))
          endif
        enddo
      endif
    enddo
    return
  end subroutine findMultinegTransforms


  ! FIND THE SHIFTS OF EACH PIECE THAT BEST ALIGN THE PIECES
  !
  subroutine getBestPieceShifts()
    real*4 edgeDispMean
    !
    ! Set the multineg flag to indicate that there are h transforms
    multng = .true.
    do ixy = 1, 2
      edgeDispMean = 0.
      do iedge = 1, nedge(ixy)
        if (izPcList(ipieceLower(iedge, ixy)) == izSect) then
          !
          ! need displacements implied by edges, unless this is to be
          ! done by old cross-correlation only
          !
          if (.not.xcLegacy) then
            call edgeSwap(iedge, ixy, inde)
            !
            ! compute the mean current displacement of upper relative
            ! to lower implied by the d[xy]mean
            !
            sumX = 0.
            sumY = 0.
            do ix = 1, nxGrBf(inde)
              do iy = 1, nyGrBf(inde)
                sumX = sumX + dxGrBf(ix, iy, inde)
                sumY = sumY + dyGrBf(ix, iy, inde)
              enddo
            enddo
            dxGridMean(iedge, ixy) = sumX / (nxGrBf(inde) * nyGrBf(inde))
            dyGridMean(iedge, ixy) = sumY / (nxGrBf(inde) * nyGrBf(inde))
            !
            ! adjust these by the current displacements implied by
            ! starting and offset coordinates of the edge areas
            !
            if (ixy == 1) then
              dxGridMean(iedge, ixy) = dxGridMean(iedge, ixy) + &
                  ixOfsBf(inde) - ixGrdStBf(inde) + nxin - nXoverlap
              dyGridMean(iedge, ixy) = dyGridMean(iedge, ixy) + &
                  iyOfsBf(inde) - iyGrdStBf(inde)
            else
              dxGridMean(iedge, ixy) = dxGridMean(iedge, ixy) + &
                  ixOfsBf(inde) - ixGrdStBf(inde)
              dyGridMean(iedge, ixy) = dyGridMean(iedge, ixy) + &
                  iyOfsBf(inde) - iyGrdStBf(inde) + nyin - nyOverlap
            endif
            !
            ! But if this edge is skipped for edge functions and included
            ! for finding shifts, substitute the (read-in) correlation shift
            if (ifSkipEdge(iedge, ixy) == 1) then
              dxGridMean(iedge, ixy) = -edgeDisplaceX(iedge, ixy)
              dyGridMean(iedge, ixy) = -edgeDisplaceY(iedge, ixy)
            endif
          endif
          !
          if (.not.fromEdge .and. .not.xcReadIn .and. &
              .not.(ifSloppy == 1 .and. ifOldEdge == 0)) then
            !
            ! If the edges of the image are lousy, it's better to use
            ! correlation, so here is this option.  Compute the
            ! correlations unless doing this by edges only,
            ! if they aren't already available
            !
            if (ifSkipEdge(iedge, ixy) > 0) then
              xDisplace = 0.
              yDisplace = 0.
            else
              call shuffler(ipieceLower(iedge, ixy), indLow)
              call shuffler(ipieceUpper(iedge, ixy), indUpper)

              call getExtraIndents(ipieceLower(iedge, ixy), &
                  ipieceUpper(iedge, ixy), ixy, delIndent)
              indentXC = 0
              if (delIndent(ixy) > 0. .and. ifillTreatment == 1) &
                  indentXC = int(delIndent(ixy)) + 1
              call xcorrEdge(array(indLow), array(indUpper), &
                  ixy, xDisplace, yDisplace, xcLegacy, indentXC)
            endif
            edgeDisplaceX(iedge, ixy) = xDisplace
            edgeDisplaceY(iedge, ixy) = yDisplace
          endif
          ! write(*,'(1x,a,2i4,a,2f8.2,a,2f8.2)') &
          ! char(ixy+ichar('W')) //' edge, pieces' &
          ! , ipiecelower(iedge, ixy), ipieceupper(iedge, ixy), &
          ! '  dxygridmean:', dxgridmean(iedge, ixy), &
          ! dygridmean(iedge, ixy), '  xcorr:', -xdisp, -ydisp
          ! dxgridmean(iedge, ixy) =-xdisp
          ! dygridmean(iedge, ixy) =-ydisp
          ! endif
        endif
        if (ixy == 1) &
            edgeDispMean = edgeDispMean + edgeDisplaceX(iedge, ixy) / nedge(ixy)
        if (ixy == 2) &
            edgeDispMean = edgeDispMean + edgeDisplaceY(iedge, ixy) / nedge(ixy)
      enddo
      if ((ixy == 1 .and. abs(edgeDispMean * nxPieces) > nxin * 3) .or. &
          (ixy == 2 .and. abs(edgeDispMean * nyPieces) > nyin * 3)) &
          write(*,'(/,a,f8.0,a)') 'WARNING: mean edge shift of', edgeDispMean, &
          ' in '// char(ixy + ichar('W')) //' may give artifacts; consider'// &
          ' adjusting overlaps with edpiecepoint and -overlap option'
    enddo
    !
    ! If there is only one piece in one direction, then do it from
    ! correlation only unless directed to use the edge, because the error
    ! is zero in either case and it is impossible to tell which is better
    fromCorrOnly = xcLegacy .or. &
        (.not.fromEdge .and. (nxPieces == 1 .or. nyPieces == 1))
    !
    if (.not. fromEdge) then
      call find_best_shifts(edgeDisplaceX, edgeDisplaceY, limEdge, -1, izSect, hxf, &
          numBestEdge, beforeMean, beforeMax, afterMean(1), afterMax(1))
      indBest = 1
    endif
    if (.not. fromCorrOnly) then
      call find_best_shifts(dxGridMean, dyGridMean, limEdge, 1, izSect, hxf, &
          numBestEdge, beforeMean, beforeMax, afterMean(2), afterMax(2))
      indBest = 2
    endif
    !
    ! if first one was better based upon mean, redo it and reset the
    ! index to 1
    !
    if (.not.(fromCorrOnly .or. fromEdge) .and. &
        (afterMean(1) < afterMean(2))) then
      call find_best_shifts(edgeDisplaceX, edgeDisplaceY, limEdge, -1, izSect, hxf, &
          numBestEdge, beforeMean, beforeMax, afterMean(1), afterMax(1))
      indBest = 1
    endif
    if (numBestEdge > 0) write(*,'(i5,a,2f7.1,a,a,2f7.2)') numBestEdge, &
        ' edges, mean&max error before:', beforeMean, beforeMax, &
        ', after by ', edgeXcorrText(indBest), afterMean(indBest), afterMax(indBest)
    if (testMode) then
      write(*,'(a,i4,a,2f8.3,a,2f9.4)') ' section:', &
          izSect, '  gradient:', dmagPerUm(min(ilistz, numMagGrad)), &
          rotPerUm(min(ilistz, numMagGrad)), &
          '  mean, max error:', afterMean(indBest), afterMax(indBest)
      if (indBest == 1) then
        call findBestGradient(edgeDisplaceX, edgeDisplaceY, limEdge, -1, izSect, &
            dmagNew, drotNew)
      else
        call findBestGradient(dxGridMean, dyGridMean, limEdge, 1, izSect, &
            dmagNew, drotNew)
      endif
      write(*,'(a,2f9.4)') ' Total gradient implied by displacements:', &
          dmagPerUm(min(ilistz, numMagGrad)) + dmagNew, &
          rotPerUm(min(ilistz, numMagGrad)) + drotNew
    endif
    return
  end subroutine getBestPieceShifts


  ! LOADS THE EDGES AND COMPUTES EDGE FRACTIONS FOR A PIXEL
  !
  subroutine computeEdgeFractions()
    numEdgesIn = 0
    do ixy = 1, 2
      do ied = 1, numEdges(ixy)
        iedge = inEdge(ied, ixy)
        call edgeSwap(iedge, ixy, indEdge)
        indLower = inEdLower(ied, ixy)
        indEdge4(ied, ixy) = indEdge
        if (edgesSeparated) then
          if (ixy == 1) then
            edgeFrac4(ied, ixy) = 0.5 + (xInPiece(indLower) - &
                (ixGrdStBf(indEdge) + (nxGrBf(indEdge) - 1) * intGrid(1) / 2.)) / &
                min(nxGrBf(indEdge) * intGrid(1), iblend(1))
          else
            edgeFrac4(ied, ixy) = 0.5 + (yInPiece(indLower) &
                - (iyGrdStBf(indEdge) + (nyGrBf(indEdge) - 1) * intGrid(1) / 2.)) / &
                min(nyGrBf(indEdge) * intGrid(1), iblend(2))
          endif
        else
          if (ixy == 1) then
            bwOffset = max(0, nxGrBf(indEdge) * intGrid(1) - iblend(1)) / 2
            edgeStart = max(0.55 * nxin, ixGrdStBf(indEdge) - intGrid(1) / 2. +bwOffset)
            edgeEnd = min(0.45 * nxin, ixOfsBf(indEdge) + (nxGrBf(indEdge) - 0.5)* &
                intGrid(1) - bwOffset) + ixGrdStBf(indEdge) - ixOfsBf(indEdge)
            bwOffset = (xInPiece(indLower) - edgeStart) / (edgeEnd - edgeStart)
          else
            bwOffset = max(0, nyGrBf(indEdge) * intGrid(1) - iblend(2)) / 2
            edgeStart = max(0.55 * nyin, iyGrdStBf(indEdge) - intGrid(1) / 2. +bwOffset)
            edgeEnd = min(0.45 * nyin, iyOfsBf(indEdge) + (nyGrBf(indEdge) - 0.5)* &
                intGrid(1) - bwOffset) + iyGrdStBf(indEdge) - iyOfsBf(indEdge)
            bwOffset = (yInPiece(indLower) - edgeStart) / (edgeEnd - edgeStart)
          endif
          if (debug) print *,ied, ixy, edgeFrac4(ied, ixy), bwOffset
          edgeFrac4(ied, ixy) = bwOffset
        endif
        active4(ied, ixy) = edgeFrac4(ied, ixy) < .999 &
            .and. edgeFrac4(ied, ixy) > .001
        if (active4(ied, ixy)) numEdgesIn = numEdgesIn + 1
        if (edgeFrac4(ied, ixy) < 0.) edgeFrac4(ied, ixy) = 0.
        if (edgeFrac4(ied, ixy) > 1.) edgeFrac4(ied, ixy) = 1.
        if (debug) print *,ied, ixy, edgeFrac4(ied, ixy), active4(ied, ixy)
      enddo
    enddo
    return
  end subroutine computeEdgeFractions



  ! SORT OUT THE CASES OF 1, 2, 3 or 4 PIECES AND COMPUTE PIXEL
  !
  subroutine getPixelFromPieces()
    real*4 oneIntrp
    !
    if (nActiveP <= 1) then                     !ONE PIECE
      if (wll > 0.) then
        indPcOne = indP1
      elseif (wlr > 0.) then
        indPcOne = indP2
      elseif (wul > 0.) then
        indPcOne = indP3
      else
        indPcOne = indP4
      endif
      if (nActiveP == 0) indPcOne = 1
      !
      call shuffler(inPiece(indPcOne), indArray)
      pixVal = oneIntrp(array(indArray), nxin, nyin, xInPiece(indPcOne), &
          yInPiece(indPcOne), inPiece(indPcOne))
      !
      ! ONE EDGE, TWO PIECES
      !
    elseif (nActiveP == 2) then
      !
      ! find the pieces around the edge, set edge index
      !
      if (wll > 0. .and. wlr > 0.) then
        indPcLo = indP1
        indPcUp = indP2
        w1 = wll
        indEdge = indEdge4(inde12, 1)
      elseif (wul > 0. .and. wur > 0.) then
        indPcLo = indP3
        indPcUp = indP4
        w1 = wul
        indEdge = indEdge4(inde34, 1)
      elseif (wll > 0. .and. wul > 0.) then
        indPcLo = indP1
        indPcUp = indP3
        w1 = wll
        indEdge = indEdge4(inde13, 2)
      else
        indPcLo = indP2
        indPcUp = indP4
        w1 = wlr
        indEdge = indEdge4(inde24, 2)
      endif
      !
      ! set up pieces #'s, and starting coords
      !
      ipiece1 = inPiece(indPcLo)
      ipiece2 = inPiece(indPcUp)
      x1 = xInPiece(indPcLo)
      y1 = yInPiece(indPcLo)
      w2 = 1. -w1
      !
      ! set up to solve equation for (x1, y1) and (x2, y2)
      ! given their difference and the desired xg, yg
      !
      xSrcConst = xg - w1 * ixPcList(ipiece1) - w2 * ixPcList(ipiece2)
      ySrcConst = ysrc - w1 * iyPcList(ipiece1) - w2 * iyPcList(ipiece2)
      if (multng) then
        xSrcConst = xSrcConst - w1 * hxf(1, 3, ipiece1) - w2 * hxf(1, 3, ipiece2)
        ySrcConst = ySrcConst - w1 * hxf(2, 3, ipiece1) - w2 * hxf(2, 3, ipiece2)
        c11 = w1 * hxf(1, 1, ipiece1) + w2 * hxf(1, 1, ipiece2)
        c12 = w1 * hxf(1, 2, ipiece1) + w2 * hxf(1, 2, ipiece2)
        c21 = w1 * hxf(2, 1, ipiece1) + w2 * hxf(2, 1, ipiece2)
        c22 = w1 * hxf(2, 2, ipiece1) + w2 * hxf(2, 2, ipiece2)
        denom = c11 * c22 - c12 * c21
        fb11 = w2 * hxf(1, 1, ipiece2)
        fb12 = w2 * hxf(1, 2, ipiece2)
        fb21 = w2 * hxf(2, 1, ipiece2)
        fb22 = w2 * hxf(2, 2, ipiece2)
      endif
      !
      ! in this loop, use the difference between (x1, y1)
      ! and (x2, y2) at the one place to solve for values of
      ! those coords; iterate until stabilize
      !
      x1last = -100.
      y1last = -100.
      iter = 1
      do while (iter <= numIter .and. (abs(x1 - x1last) > 0.01 &
          .or. abs(y1 - y1last) > 0.01))
        call dxydGrInterp(x1, y1, indEdge, x2, y2, dden)
        x1last = x1
        y1last = y1
        dx2 = x2 - x1
        dy2 = y2 - y1
        if (multng) then
          bx = xSrcConst - fb11 * dx2 - fb12 * dy2
          by = ySrcConst - fb21 * dx2 - fb22 * dy2
          x1 = (bx * c22 - by * c12) / denom
          y1 = (by * c11 - bx * c21) / denom
        else
          x1 = xSrcConst - w2 * dx2
          y1 = ySrcConst - w2 * dy2
        endif
        x2 = x1 + dx2
        y2 = y1 + dy2
        iter = iter + 1
      enddo
      !
      ! get the pixel from the piece with the bigger weight
      ! and adjust density by the difference across edge
      ! Except average them within 4 pixels of center!
      !
      if (abs(w1 - w2) * max(iblend(1), iblend(2)) < 2.5) then
        call shuffler(ipiece1, indArray1)
        call shuffler(ipiece2, indArray2)
        pixVal = w1 * oneIntrp(array(indArray1), nxin, nyin, x1, y1, ipiece1) + &
            w2 * oneIntrp(array(indArray2), nxin, nyin, x2, y2, ipiece2)
      elseif (w1 > w2) then
        call shuffler(ipiece1, indArray)
        pixVal = oneIntrp(array(indArray), nxin, nyin, x1, y1, ipiece1) + w2 * dden
      else
        call shuffler(ipiece2, indArray)
        pixVal = oneIntrp(array(indArray), nxin, nyin, x2, y2, ipiece2) - w1 * dden
      endif
      !
      ! THREE PIECES AND TWO EDGES
      !
    elseif (nActiveP == 3) then
      !
      ! now decide between the three cases of divergent, serial, or
      ! convergent functions, and assign pieces and weights accordingly
      !
      ! DIVERGENT: if lower piece is same for x and y edge
      !
      if (wur <= 0.) then
        w1 = wll
        w2 = wlr
        w3 = wul
        jndP1 = indP1
        jndP2 = indP2
        jndP3 = indP3
        ind12edg = indEdge4(inde12, 1)
        ind13Edge = indEdge4(inde13, 2)
        iCornType = 1
        !
        ! CONVERGENT: if upper piece same for x and y edge
        !
      elseif (wll <= 0.) then
        w3 = wur
        w1 = wul
        w2 = wlr
        jndP1 = indP3
        jndP2 = indP2
        jndP3 = indP4
        ind13Edge = indEdge4(inde34, 1)
        ind23edg = indEdge4(inde24, 2)
        iCornType = 2
        !
        ! SERIAL: lower of one is the upper of the other
        !
      else
        w1 = wll
        w3 = wur
        jndP1 = indP1
        jndP3 = indP4
        if (wlr <= 0.) then
          w2 = wul
          jndP2 = indP3
          ind12edg = indEdge4(inde13, 2)
          ind23edg = indEdge4(inde34, 1)
        else
          w2 = wlr
          jndP2 = indP2
          ind12edg = indEdge4(inde12, 1)
          ind23edg = indEdge4(inde24, 2)
        endif
        iCornType = 3
      endif

      ipiece1 = inPiece(jndP1)
      ipiece2 = inPiece(jndP2)
      ipiece3 = inPiece(jndP3)
      x1 = xInPiece(jndP1)
      y1 = yInPiece(jndP1)
      wMax = max(w1, w2, w3)
      !
      ! set up to solve equations for new (x1, y1), (x2, y2)
      ! and (x3, y3) given the differences between them and
      ! the desired weighted coordinate (xg, yg)
      !
      xSrcConst = xg - w1 * ixPcList(ipiece1) - &
          w2 * ixPcList(ipiece2) - w3 * ixPcList(ipiece3)
      ySrcConst = ysrc - w1 * iyPcList(ipiece1) - &
          w2 * iyPcList(ipiece2) - w3 * iyPcList(ipiece3)
      if (multng) then
        xSrcConst = xSrcConst - w1 * hxf(1, 3, ipiece1) - w2 * hxf(1, 3, ipiece2) -  &
            w3 * hxf(1, 3, ipiece3)
        ySrcConst = ySrcConst - w1 * hxf(2, 3, ipiece1) - w2 * hxf(2, 3, ipiece2) -  &
            w3 * hxf(2, 3, ipiece3)
        c11 = w1 * hxf(1, 1, ipiece1) + w2 * hxf(1, 1, ipiece2) + w3 * hxf(1, 1, ipiece3)
        c12 = w1 * hxf(1, 2, ipiece1) + w2 * hxf(1, 2, ipiece2) + w3 * hxf(1, 2, ipiece3)
        c21 = w1 * hxf(2, 1, ipiece1) + w2 * hxf(2, 1, ipiece2) + w3 * hxf(2, 1, ipiece3)
        c22 = w1 * hxf(2, 2, ipiece1) + w2 * hxf(2, 2, ipiece2) + w3 * hxf(2, 2, ipiece3)
        denom = c11 * c22 - c12 * c21
        f2b11 = w2 * hxf(1, 1, ipiece2)
        f2b12 = w2 * hxf(1, 2, ipiece2)
        f2b21 = w2 * hxf(2, 1, ipiece2)
        f2b22 = w2 * hxf(2, 2, ipiece2)
        f3b11 = w3 * hxf(1, 1, ipiece3)
        f3b12 = w3 * hxf(1, 2, ipiece3)
        f3b21 = w3 * hxf(2, 1, ipiece3)
        f3b22 = w3 * hxf(2, 2, ipiece3)
      endif
      !
      ! do iteration, starting with coordinates and solving
      ! for new coordinates until convergence
      !
      x1last = -100.
      y1last = -100.
      iter = 1
      do while(iter <= numIter .and. (abs(x1 - x1last) > 0.01 &
          .or. abs(y1 - y1last) > 0.01))
        if (iCornType == 1) then
          !
          ! divergent case
          !
          call dxydGrInterp(x1, y1, ind12edg, x2, y2, dden12)
          call dxydGrInterp(x1, y1, ind13Edge, x3, y3, dden13)
        elseif (iCornType == 2) then
          !
          ! convergent case
          !
          call dxydGrInterp(x1, y1, ind13Edge, x3, y3, dden13)
          !
          if (iter == 1) then
            x2 = x3 + ixGrdStBf(ind23edg) - ixOfsBf(ind23edg)
            y2 = y3 + iyGrdStBf(ind23edg) - iyOfsBf(ind23edg)
          endif
          call dxydGrInterp(x2, y2, ind23edg, x3t, y3t, dden23)
          x2 = x2 + x3 - x3t
          y2 = y2 + y3 - y3t
        else
          !
          ! serial case
          !
          call dxydGrInterp(x1, y1, ind12edg, x2, y2 , dden12)
          call dxydGrInterp(x2, y2, ind23edg, x3, y3 , dden23)
        endif
        !
        ! solve equations for new coordinates
        !
        x1last = x1
        y1last = y1
        dx2 = x2 - x1
        dy2 = y2 - y1
        dx3 = x3 - x1
        dy3 = y3 - y1
        if (multng) then
          bx = xSrcConst - f2b11 * dx2 - f2b12 * dy2 - f3b11 * dx3 - f3b12 * dy3
          by = ySrcConst - f2b21 * dx2 - f2b22 * dy2 - f3b21 * dx3 - f3b22 * dy3
          x1 = (bx * c22 - by * c12) / denom
          y1 = (by * c11 - bx * c21) / denom
        else
          x1 = xSrcConst - w2 * dx2 - w3 * dx3
          y1 = ySrcConst - w2 * dy2 - w3 * dy3
        endif
        x2 = x1 + dx2
        y2 = y1 + dy2
        x3 = x1 + dx3
        y3 = y1 + dy3
        iter = iter + 1
      enddo
      !
      ! take pixel from the piece with the highest weight
      !
      if (w1 == wMax) then
        call shuffler(ipiece1, indArray)
        pixVal = oneIntrp(array(indArray), nxin, nyin, x1, y1, ipiece1)
      elseif (w2 == wMax) then
        call shuffler(ipiece2, indArray)
        pixVal = oneIntrp(array(indArray), nxin, nyin, x2, y2, ipiece2)
      else
        call shuffler(ipiece3, indArray)
        pixVal = oneIntrp(array(indArray), nxin, nyin, x3, y3, ipiece3)
      endif
      !
      ! Adjust for differences in mean density: divergent
      !
      if (iCornType == 1) then
        if (w1 == wMax) then
          pixVal = pixVal + w2 * dden12 + w3 * dden13
        elseif (w2 == wMax) then
          pixVal = pixVal + (w2 - 1.) * dden12 + w3 * dden13
        else
          pixVal = pixVal + w2 * dden12 + (w3 - 1.) * dden13
        endif
        !
        ! convergent
        !
      elseif (iCornType == 2) then
        if (w1 == wMax) then
          pixVal = pixVal - (w1 - 1.) * dden13 - w2 * dden23
        elseif (w2 == wMax) then
          pixVal = pixVal - w1 * dden13 - (w2 - 1.) * dden23
        else
          pixVal = pixVal - w1 * dden13 - w2 * dden23
        endif
        !
        ! serial
        !
      else
        !
        if (w1 == wMax) then
          pixVal = pixVal - (w1 - 1.) * dden12 + w3 * dden23
        elseif (w2 == wMax) then
          pixVal = pixVal - w1 * dden12 + w3 * dden23
        else
          pixVal = pixVal - w1 * dden12 + (w3 - 1.) * dden23
        endif
      endif
    else
      !
      ! FOUR PIECES, THREE EDGES USED FOR SOLUTION
      !
      ! First, need to have only 3 active edges, so if
      ! there are four, knock one out based on ex and ey
      !
      if (numEdgesIn == 4) then
        emin = min(ex, 1. -ex, ey, 1. -ey)
        if (ex == emin) then
          active4(inde34, 1) = .false.
        elseif (1. -ex == emin) then
          active4(inde12, 1) = .false.
        elseif (ey == emin) then
          active4(inde24, 2) = .false.
        else
          active4(inde13, 2) = .false.
        endif
      endif
      !
      ! here there is always a serial chain from ll to ur, through either ul
      ! or lr, and in each case the fourth piece is either divergent from
      ! the first or converging on the fourth
      !
      w1 = wll
      w3 = wur
      if (.not.active4(inde12, 1) .or. &
          .not.active4(inde24, 2)) then
        w2 = wul
        w4 = wlr
        jndP2 = indP3
        jndP4 = indP2
        ind12edg = indEdge4(inde13, 2)
        ind23edg = indEdge4(inde34, 1)
        if (.not.active4(inde24, 2)) then
          ind14edg = indEdge4(inde12, 1)
          iCornType = 1                           !divergent
        else
          ind43edg = indEdge4(inde24, 2)
          iCornType = 2                           !convergent
        endif
      else
        w2 = wlr
        w4 = wul
        jndP2 = indP2
        jndP4 = indP3
        ind12edg = indEdge4(inde12, 1)
        ind23edg = indEdge4(inde24, 2)
        if (.not.active4(inde34, 1)) then
          ind14edg = indEdge4(inde13, 2)
          iCornType = 1
        else
          ind43edg = indEdge4(inde34, 1)
          iCornType = 2
        endif
      endif
      !
      ipiece1 = inPiece(indP1)
      ipiece2 = inPiece(jndP2)
      ipiece3 = inPiece(indP4)
      ipiece4 = inPiece(jndP4)
      x1 = xInPiece(indP1)
      y1 = yInPiece(indP1)
      wMax = max(w1, w2, w3, w4)
      !
      ! set up to solve equations for new (x1, y1), (x2, y2)
      ! (x3, y3), and (x4, y4) given the differences between
      ! them and the desired weighted coordinate (xg, yg)
      !
      xSrcConst = xg - w1 * ixPcList(ipiece1) - w2 * ixPcList(ipiece2) &
          - w3 * ixPcList(ipiece3) - w4 * ixPcList(ipiece4)
      ySrcConst = ysrc - w1 * iyPcList(ipiece1) - w2 * iyPcList(ipiece2) &
          - w3 * iyPcList(ipiece3) - w4 * iyPcList(ipiece4)
      if (multng) then
        xSrcConst = xSrcConst - w1 * hxf(1, 3, ipiece1) - w2 * hxf(1, 3, ipiece2) &
            - w3 * hxf(1, 3, ipiece3) - w4 * hxf(1, 3, ipiece4)
        ySrcConst = ySrcConst - w1 * hxf(2, 3, ipiece1) - w2 * hxf(2, 3, ipiece2) &
            - w3 * hxf(2, 3, ipiece3) - w4 * hxf(2, 3, ipiece4)
        c11 = w1 * hxf(1, 1, ipiece1) + w2 * hxf(1, 1, ipiece2) &
            + w3 * hxf(1, 1, ipiece3) + w4 * hxf(1, 1, ipiece4)
        c12 = w1 * hxf(1, 2, ipiece1) + w2 * hxf(1, 2, ipiece2) &
            + w3 * hxf(1, 2, ipiece3) + w4 * hxf(1, 2, ipiece4)
        c21 = w1 * hxf(2, 1, ipiece1) + w2 * hxf(2, 1, ipiece2) &
            + w3 * hxf(2, 1, ipiece3) + w4 * hxf(2, 1, ipiece4)
        c22 = w1 * hxf(2, 2, ipiece1) + w2 * hxf(2, 2, ipiece2) &
            + w3 * hxf(2, 2, ipiece3) + w4 * hxf(2, 2, ipiece4)
        denom = c11 * c22 - c12 * c21
        f2b11 = w2 * hxf(1, 1, ipiece2)
        f2b12 = w2 * hxf(1, 2, ipiece2)
        f2b21 = w2 * hxf(2, 1, ipiece2)
        f2b22 = w2 * hxf(2, 2, ipiece2)
        f3b11 = w3 * hxf(1, 1, ipiece3)
        f3b12 = w3 * hxf(1, 2, ipiece3)
        f3b21 = w3 * hxf(2, 1, ipiece3)
        f3b22 = w3 * hxf(2, 2, ipiece3)
        f4b11 = w4 * hxf(1, 1, ipiece4)
        f4b12 = w4 * hxf(1, 2, ipiece4)
        f4b21 = w4 * hxf(2, 1, ipiece4)
        f4b22 = w4 * hxf(2, 2, ipiece4)
      endif
      !
      ! do iteration, starting with coordinates and solving
      ! for new coordinates until convergence
      !
      x1last = -100.
      y1last = -100.
      iter = 1
      do while(iter <= numIter .and. (abs(x1 - x1last) > 0.01 &
          .or. abs(y1 - y1last) > 0.01))
        call dxydGrInterp(x1, y1, ind12edg, x2, y2 , dden12)
        call dxydGrInterp(x2, y2, ind23edg, x3, y3 , dden23)
        if (iCornType == 1) then
          !
          ! divergent case
          !
          call dxydGrInterp(x1, y1, ind14edg, x4, y4 , dden14)
        else
          !
          ! convergent case
          !
          if (iter == 1) then
            x4 = x3 + ixGrdStBf(ind43edg) - ixOfsBf(ind43edg)
            y4 = y3 + iyGrdStBf(ind43edg) - iyOfsBf(ind43edg)
          endif
          call dxydGrInterp(x4, y4, ind43edg, x3t, y3t , dden43)
          x4 = x4 + x3 - x3t
          y4 = y4 + y3 - y3t
        endif
        !
        ! solve equations for new coordinates
        !
        x1last = x1
        y1last = y1
        dx2 = x2 - x1
        dy2 = y2 - y1
        dx3 = x3 - x1
        dy3 = y3 - y1
        dx4 = x4 - x1
        dy4 = y4 - y1
        if (multng) then
          bx = xSrcConst - f2b11 * dx2 - f2b12 * dy2 &
              - f3b11 * dx3 - f3b12 * dy3 - f4b11 * dx4 - f4b12 * dy4
          by = ySrcConst - f2b21 * dx2 - f2b22 * dy2 &
              - f3b21 * dx3 - f3b22 * dy3 - f4b21 * dx4 - f4b22 * dy4
          x1 = (bx * c22 - by * c12) / denom
          y1 = (by * c11 - bx * c21) / denom
        else
          x1 = xSrcConst - w2 * dx2 - w3 * dx3 - w4 * dx4
          y1 = ySrcConst - w2 * dy2 - w3 * dy3 - w4 * dy4
        endif
        x2 = x1 + dx2
        y2 = y1 + dy2
        x3 = x1 + dx3
        y3 = y1 + dy3
        x4 = x1 + dx4
        y4 = y1 + dy4
        iter = iter + 1
      enddo
      !
      ! take pixel from the piece with the highest weight
      ! and adjust density appropriately for case
      !
      if (w1 == wMax) then
        call shuffler(ipiece1, indArray)
        pixVal = oneIntrp(array(indArray), nxin, nyin, x1, y1, ipiece1)
        if (iCornType == 1) then
          pixVal = pixVal + (w2 + w3) * dden12 + w3 * dden23 + w4 * dden14
        else
          pixVal = pixVal + (1. -w1) * dden12 + (w3 + w4) * dden23 - w4 * dden43
        endif
      elseif (w2 == wMax) then
        call shuffler(ipiece2, indArray)
        pixVal = oneIntrp(array(indArray), nxin, nyin, x2, y2, ipiece2)
        if (iCornType == 1) then
          pixVal = pixVal + (w2 + w3 - 1.) * dden12 + w3 * dden23 + w4 * dden14
        else
          pixVal = pixVal - w1 * dden12 + (w3 + w4) * dden23 - w4 * dden43
        endif
      elseif (w3 == wMax) then
        call shuffler(ipiece3, indArray)
        pixVal = oneIntrp(array(indArray), nxin, nyin, x3, y3, ipiece3)
        if (iCornType == 1) then
          pixVal = pixVal + (w2 + w3 - 1.) * dden12 + (w3 - 1.) * dden23 + w4 * dden14
        else
          pixVal = pixVal - w1 * dden12 + (w3 + w4 - 1.) * dden23 - w4 * dden43
        endif
      else
        call shuffler(ipiece4, indArray)
        pixVal = oneIntrp(array(indArray), nxin, nyin, x4, y4, ipiece4)
        if (iCornType == 1) then
          pixVal = pixVal + (w2 + w3) * dden12 + w3 * dden23 + (w4 - 1.) * dden14
        else
          pixVal = pixVal - w1 * dden12 + (w3 + w4 - 1.) * dden23 - (w4 - 1.) * dden43
        endif
      endif
    endif
  end subroutine getPixelFromPieces

end program blendmont

! GET INDICES OF PIECES AND EDGES AND THE WEIGHTING OF EACH PIECE
!
subroutine getPieceIndicesAndWeighting(useEdges)
  use blendvars
  implicit none
  logical*4 useEdges
  real*4 er3, eb3, el4, eb4, fx, fy, dr1, dt1, dl2, dt2, dr3, db3, dl4, db4, dla, dra
  real*4 dba, dta, ax, ay, f12, f13, f34, f24, er1, et1, el2, et2
  integer*4 ipFrom, ipXto, ipYto, ixyDisjoint, ied, ipc, i, indEdge, ixy, jedge
  real*4 bwOffset, edgeStart, edgeEnd, wSum

  ! for now, numbers 1 to 4 represent lower left, lower right, upper left,
  ! and upper right
  if (numPieces > 1) then
    indP1 = 0                                 !index of piece 1, 2, 3, or 4
    indP2 = 0                                 !if the point is in them
    indP3 = 0
    indP4 = 0
    inde12 = 0                                !index of edges between 1 and 2
    inde13 = 0                                !1 and 3, etc
    inde34 = 0
    inde24 = 0
    f12 = 0.                                  !edge fractions between 1 and 2
    f13 = 0.                                  !1 and 3, etc
    f34 = 0.
    f24 = 0.
    er1 = 0.                                  !end fractions to right and top
    et1 = 0.                                  !of piece 1, left and bottom of
    el2 = 0.                                  !piece 3, etc
    et2 = 0.
    er3 = 0.
    eb3 = 0.
    el4 = 0.
    eb4 = 0.
    fx = 0.                                   !the composite edge fractions
    fy = 0.                                   !in x and y directions
    if (numPieces > 2) then
      do ipc = 1, numPieces
        i = 2 * (inpYframe(ipc) - minYframe) + 1 + &
            inpXframe(ipc) - minXframe
        indP1234(i) = ipc
      enddo
      do ied = 1, numEdges(1)
        if (inEdLower(ied, 1) == indP1) inde12 = ied
        if (inEdLower(ied, 1) == indP3) inde34 = ied
      enddo
      do ied = 1, numEdges(2)
        if (inEdLower(ied, 2) == indP1) inde13 = ied
        if (inEdLower(ied, 2) == indP2) inde24 = ied
      enddo
      if (debug) print *,inde12, inde34, inde13, inde24
      if (inde12 .ne. 0) f12 = edgeFrac4(inde12, 1)
      if (inde34 .ne. 0) f34 = edgeFrac4(inde34, 1)
      if (inde13 .ne. 0) f13 = edgeFrac4(inde13, 2)
      if (inde24 .ne. 0) f24 = edgeFrac4(inde24, 2)
      if (debug) write(*,'(a,4i6,a,4f8.4)') 'piece 1234', &
          inPiece(indP1), inPiece(indP2), inPiece(indP3), &
          inPiece(indP4), '  edge fractions', f12, f34, f13, f24
    else
      !
      ! two piece case - identify as upper or lower to
      ! simplify computation of end fractions
      !
      if (numEdges(1) > 0) then
        fx = edgeFrac4(1, 1)
        if (0.5 * (yInPiece(1) + yInPiece(2)) > nyin / 2) then
          inde12 = 1
          if (xInPiece(1) > xInPiece(2)) then
            indP1 = 1
            indP2 = 2
          else
            indP1 = 2
            indP2 = 1
          endif
        else
          inde34 = 1
          fy = 1.
          if (xInPiece(1) > xInPiece(2)) then
            indP3 = 1
            indP4 = 2
          else
            indP3 = 2
            indP4 = 1
          endif
        endif
      else
        fy = edgeFrac4(1, 2)
        if (0.5 * (xInPiece(1) + xInPiece(2)) > nxin / 2) then
          inde13 = 1
          if (yInPiece(1) > yInPiece(2)) then
            indP1 = 1
            indP3 = 2
          else
            indP1 = 2
            indP3 = 1
          endif
        else
          inde24 = 1
          fx = 1.
          if (yInPiece(1) > yInPiece(2)) then
            indP2 = 1
            indP4 = 2
          else
            indP2 = 2
            indP4 = 1
          endif
        endif
      endif
    endif
    !
    ! get distance to top, right, bottom, or left edges
    ! as needed for each piece, and compute end fractions
    !
    if (indP1 > 0) then
      dr1 = nxin - 1. -xInPiece(indP1)
      dt1 = nyin - 1. -yInPiece(indP1)
      er1 = (dr1 + 1.) / iblend(1)
      et1 = dt1 / iblend(2)
    endif
    if (indP2 > 0) then
      dl2 = xInPiece(indP2)
      dt2 = nyin - 1. -yInPiece(indP2)
      el2 = (dl2 + 1.) / iblend(1)
      et2 = (dt2 + 1.) / iblend(2)
    endif
    if (indP3 > 0) then
      dr3 = nxin - 1. -xInPiece(indP3)
      db3 = yInPiece(indP3)
      er3 = (dr3 + 1.) / iblend(1)
      eb3 = (db3 + 1.) / iblend(2)
    endif
    if (indP4 > 0) then
      dl4 = xInPiece(indP4)
      db4 = yInPiece(indP4)
      el4 = (dl4 + 1.) / iblend(1)
      eb4 = (db4 + 1.) / iblend(2)
    endif
    !
    ! If 3 pieces, check for a disjoint edge: use previous state if pieces
    ! match the last time
    if (numPieces == 3 .and. anyDisjoint(minXframe, minYframe)) then
      if (inPiece(indP1) == lastp1 .and. inPiece(indP2) == lastp2 .and. &
          inPiece(indP3) == lastp3 .and. inPiece(indP4) == lastp4) then
        ixyDisjoint = lastxyDisjoint
      else
        !
        ! if pieces have changed, analyze edge starts and ends.
        ! First find the missing edges and set some indexes
        ixyDisjoint = 0
        if (indP1 == 0) then
          inEdge(2, 1) = iedgeLower(inPiece(indP2), 1)
          inEdge(2, 2) = iedgeLower(inPiece(indP3), 2)
          ipFrom = indP4
          ipXto = indP2
          ipYto = indP3
        else if (indP2 == 0) then
          inEdge(2, 1) = iedgeUpper(inPiece(indP1), 1)
          inEdge(2, 2) = iedgeLower(inPiece(indP4), 2)
          ipFrom = indP3
          ipXto = indP1
          ipYto = indP4
        else if (indP3 == 0) then
          inEdge(2, 1) = iedgeLower(inPiece(indP4), 1)
          inEdge(2, 2) = iedgeUpper(inPiece(indP1), 2)
          ipFrom = indP2
          ipXto = indP4
          ipYto = indP1
        else if (indP4 == 0) then
          inEdge(2, 1) = iedgeUpper(inPiece(indP3), 1)
          inEdge(2, 2) = iedgeUpper(inPiece(indP2), 2)
          ipFrom = indP1
          ipXto = indP3
          ipYto = indP2
        endif
        !
        ! Replace with used edge numbers
        if (useEdges) then
          do ixy = 1, 2
            call findEdgeToUse(inEdge(2, ixy), ixy, jedge)
            if (jedge .ne. 0) inEdge(2, ixy) = jedge
          enddo
        endif
        !
        ! Analyze X edge first then Y edge - don't worry if both are bad
        if (inEdge(2, 1) > 0 .and. inEdge(2, 2) > 0) then
          !
          ! Get start of already included X edge and translate it into
          ! piece on missing X edge, and get end of the missing X edge
          ! in that piece.  If end is before start, it is disjoint
          if (indP1 == 0 .or. indP3 == 0) then
            indEdge = indEdge4(1, 1)
            bwOffset = max(0, nxGrBf(indEdge) * intGrid(1) - iblend(1)) / 2
            edgeStart = ixOfsBf(indEdge) - intGrid(1) / 2. + bwOffset + &
                xInPiece(ipXto) - xInPiece(ipFrom)
            call edgeSwap(inEdge(2, 1), 1, indEdge)
            bwOffset = max(0, nxGrBf(indEdge) * intGrid(1) - iblend(1)) / 2
            edgeEnd = ixOfsBf(indEdge) + (nxGrBf(indEdge) - 0.5) * intGrid(1) - bwOffset
            if (edgeEnd <= edgeStart) ixyDisjoint = 1
          else
            !
            ! Other two cases, compare end of existing to start of missing
            indEdge = indEdge4(1, 1)
            bwOffset = max(0, nxGrBf(indEdge) * intGrid(1) - iblend(1)) / 2
            edgeEnd = ixGrdStBf(indEdge) + (nxGrBf(indEdge) - 0.5) * intGrid(1) - &
                bwOffset +  xInPiece(ipXto) - xInPiece(ipFrom)
            call edgeSwap(inEdge(2, 1), 1, indEdge)
            bwOffset = max(0, nxGrBf(indEdge) * intGrid(1) - iblend(1)) / 2
            edgeStart = ixGrdStBf(indEdge) - intGrid(1) / 2. + bwOffset
            if (edgeEnd <= edgeStart) ixyDisjoint = 1
          endif
          !
          ! Similarly two cases for Y
          if (indP1 == 0 .or. indP2 == 0) then
            indEdge = indEdge4(1, 2)
            bwOffset = max(0, nyGrBf(indEdge) * intGrid(1) - iblend(2)) / 2
            edgeStart = iyOfsBf(indEdge) - intGrid(1) / 2. + bwOffset + &
                yInPiece(ipYto) - yInPiece(ipFrom)
            call edgeSwap(inEdge(2, 2), 2, indEdge)
            bwOffset = max(0, nyGrBf(indEdge) * intGrid(1) - iblend(2)) / 2
            edgeEnd = iyOfsBf(indEdge) + (nyGrBf(indEdge) - 0.5) * intGrid(1) - bwOffset
            if (edgeEnd <= edgeStart) ixyDisjoint = 2
          else
            indEdge = indEdge4(1, 2)
            bwOffset = max(0, nyGrBf(indEdge) * intGrid(1) - iblend(2)) / 2
            edgeEnd = iyGrdStBf(indEdge) + (nyGrBf(indEdge) - 0.5) * intGrid(1) - &
                bwOffset +  yInPiece(ipYto) - yInPiece(ipFrom)
            call edgeSwap(inEdge(2, 2), 2, indEdge)
            bwOffset = max(0, nyGrBf(indEdge) * intGrid(1) - iblend(2)) / 2
            edgeStart = iyGrdStBf(indEdge) - intGrid(1) / 2. + bwOffset
            if (edgeEnd <= edgeStart) ixyDisjoint = 2
          endif
          !
          ! If one was found, compute the start and end for the half of
          ! the edge to be used on the other axis
          if (ixyDisjoint == 1) then
            indEdge = indEdge4(1, 2)
            bwOffset = max(0, nyGrBf(indEdge) * intGrid(1) - iblend(2)) / 2
            if (indP1 == 0 .or. indP2 == 0) then
              startSkew = iyOfsBf(indEdge) + (nyGrBf(indEdge) - 1) * intGrid(1) / 2.
              endSkew = iyOfsBf(indEdge) + (nyGrBf(indEdge) - 0.5) * intGrid(1) - &
                  bwOffset
            else
              startSkew = iyGrdStBf(indEdge) - intGrid(1) / 2. + bwOffset
              endSkew = iyGrdStBf(indEdge) + (nyGrBf(indEdge) - 1) * intGrid(1) / 2.
            endif
          else if (ixyDisjoint == 2) then
            indEdge = indEdge4(1, 1)
            bwOffset = max(0, nxGrBf(indEdge) * intGrid(1) - iblend(1)) / 2
            if (indP1 == 0 .or. indP3 == 0) then
              startSkew = ixOfsBf(indEdge) + (nxGrBf(indEdge) - 1) * intGrid(1) / 2.
              endSkew = ixOfsBf(indEdge) + (nxGrBf(indEdge) - 0.5) * intGrid(1) &
                  - bwOffset
            else
              startSkew = ixGrdStBf(indEdge) - intGrid(1) / 2. + bwOffset
              endSkew = ixGrdStBf(indEdge) + (nxGrBf(indEdge) - 1) * intGrid(1) / 2.
            endif
          endif
        endif
        !
        ! Save for future use
        lastp1 = inPiece(indP1)
        lastp2 = inPiece(indP2)
        lastp3 = inPiece(indP3)
        lastp4 = inPiece(indP4)
        lastxyDisjoint = ixyDisjoint
      endif
      !
      ! Now if there is disjoint edge, modify edge and end fractions to
      ! start at the midpoint of the piece
      if (ixyDisjoint == 1) then
        if (indP1 == 0) then
          edgeFrac4(1, 2) = max(0., min(1., (yInPiece(indP4) - &
              startSkew) / (endSkew - startSkew)))
          db4 = max(0., yInPiece(indP4) - startSkew)
          eb4 = (db4 + 1.) / iblend(2)
          if (debug) write(*,'(a,3f8.1,a,f8.4,f8.1,f8.4)') 's&e skew, y', &
              startSkew, endSkew, yInPiece(indP4), &
              ' mod ef, db4, eb4', edgeFrac4(1, 2), db4, eb4
        else if (indP2 == 0) then
          edgeFrac4(1, 2) = max(0., min(1., (yInPiece(indP3) - &
              startSkew) / (endSkew - startSkew)))
          db3 = max(0., yInPiece(indP3) - startSkew)
          eb3 = (db3 + 1.) / iblend(2)
        else if (indP3 == 0) then
          edgeFrac4(1, 2) = max(0., min(1., (yInPiece(indP2) - &
              startSkew) / (endSkew - startSkew)))
          dt2 = max(0., endSkew - yInPiece(indP2))
          et2 = (dt2 + 1.) / iblend(2)
        else
          edgeFrac4(1, 2) = max(0., min(1., (yInPiece(indP1) - &
              startSkew) / (endSkew - startSkew)))
          dt1 = max(0., endSkew - yInPiece(indP1))
          et1 = (dt1 + 1.) / iblend(2)
        endif
      else if (ixyDisjoint == 2) then
        if (indP1 == 0) then
          edgeFrac4(1, 1) = max(0., min(1., (xInPiece(indP4) - &
              startSkew) / (endSkew - startSkew)))
          dl4 = max(0., xInPiece(indP4) - startSkew)
          el4 = (dl4 + 1.) / iblend(1)
        else if (indP2 == 0) then
          edgeFrac4(1, 1) = max(0., min(1., (xInPiece(indP3) - &
              startSkew) / (endSkew - startSkew)))
          dr3 = max(0., endSkew - xInPiece(indP3))
          er3 = (dr3 + 1.) / iblend(1)
        else if (indP3 == 0) then
          edgeFrac4(1, 1) = max(0., min(1., (xInPiece(indP2) - &
              startSkew) / (endSkew - startSkew)))
          dl2 = max(0., xInPiece(indP2) - startSkew)
          el2 = (dl2 + 1.) / iblend(1)
        else
          edgeFrac4(1, 1) = max(0., min(1., (xInPiece(indP1) - &
              startSkew) / (endSkew - startSkew)))
          dr1 = max(0., endSkew - xInPiece(indP1))
          er1 = (dr3 + 1.) / iblend(1)
        endif
      endif
    endif
    !
    ! If there are 4 pieces, fx and fy are weighted sum of the f's in the
    ! two overlapping edges.  The weights ex and ey are modified
    ! fractional distances across the overlap zone.  First get the
    ! distances to the borders of the overlap zone (dla, etc) and
    ! use them to get absolute fractional distances across zone, ax in the
    ! y direction and ay in the x direction.  They are used to make ex
    ! slide from being a distance across the lower overlap to being
    ! a distance across the upper overlap.  This gives continuity with the
    ! edge in 2 or 3-piece cases
    !
    if (numPieces == 4) then
      dla = min(dl2, dl4)
      dra = min(dr1, dr3)
      dba = min(db3, db4)
      dta = min(dt1, dt2)
      ax = dba / (dta + dba)
      ay = dla / (dra + dla)
      ex = ((1 - ay) * db3 + ay * db4) / ((1 - ay) * (dt1 + db3) + ay * (dt2 + db4))
      ey = ((1 - ax) * dl2 + ax * dl4) / ((1 - ax) * (dr1 + dl2) + ax * (dr3 + dl4))
      fx = (1 - ex) * f12 + ex * f34
      fy = (1 - ey) * f13 + ey * f24
    elseif (numPieces == 3) then
      !
      ! Three-piece case is simple, only two edges
      !
      fx = edgeFrac4(1, 1)
      fy = edgeFrac4(1, 2)
    endif
    !
    ! weighting factors are a product of the two f's,
    ! attenuated if necessary by fractional distance to
    ! end of piece, then normalized to sum to 1.
    !
    wll = min(1 - fx, et1) * min(1 - fy, er1)
    wlr = min(fx, et2) * min(1 - fy, el2)
    wul = min(1 - fx, eb3) * min(fy, er3)
    wur = min(fx, eb4) * min(fy, el4)
    wSum = wll + wlr + wul + wur
    if (wSum > 0.) then
      wll = wll / wSum
      wlr = wlr / wSum
      wul = wul / wSum
      wur = wur / wSum
    endif
    !
    ! count up active pieces implied by the w's
    !
    nActiveP = 0
    if (wll > 0.) nActiveP = nActiveP + 1
    if (wlr > 0.) nActiveP = nActiveP + 1
    if (wul > 0.) nActiveP = nActiveP + 1
    if (wur > 0.) nActiveP = nActiveP + 1
    !
    ! filter out cross-corner 2-piece cases, pick the
    ! piece where the point is most interior
    !
    if (nActiveP == 2 .and. wll * wur > 0.) then
      if (min(dt1, dr1) < min(db4, dl4)) then
        wll = 0.
      else
        wur = 0.
      endif
      nActiveP = 1
    elseif (nActiveP == 2 .and. wul * wlr > 0.) then
      if (min(dt2, dl2) < min(db3, dr3)) then
        wlr = 0.
      else
        wul = 0.
      endif
      nActiveP = 1
    endif
  else
    !
    ! the one-piece case, avoid all that computation
    !
    nActiveP = 1
    indP1 = 1
    wll = 1.
  endif
  if (debug) write(*,'(a,i3,a,4f8.4)') 'Active', nActiveP, '  weights', wll, wlr, wul, wur
  return
end subroutine getPieceIndicesAndWeighting
