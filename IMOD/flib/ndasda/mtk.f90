! MTK analyzes the distribution of distances of closest approach
! between objects in 3 dimensions.
!
! See the man page for details
!
! David Mastronarde  November 1991
! General object version March 2000
!
! $Id$
  !
program mtk
  call plax_initialize('mtk')
  call exit(0)
end program mtk

! 10/29/15: notes after fixing crashes on Mac.
! Stack is very limited in threads on OS X, probably limited on Linux too.
! The key is to minimize arrays on the stack.
! Variables in modules do not count so must be on the heap as for one in commons
! The dynamic allocation is probably unnecessary but would be a better way to fail if
! there isn't enough memory, and is a first step to accommodating any model size

module mtkredovars
  implicit none
  integer LIMGRAPHS, LIMBINS, LIMWOBJ, LIMXYZ, LIMREGION, LIMTYPE, LIMFLAG
  integer LIMPROBS, LIMPROBSETS
  parameter (LIMGRAPHS = 50, LIMBINS = 1001, LIMWOBJ = 300000, LIMXYZ = 20000000, &
      LIMREGION = 200)
  ! Old notes on memory
  ! limwobj * 36 + limxyz * 12
  ! limgraphs * 4 * (12 + 3 * limtyp + 2 * limbins)
  parameter (LIMTYPE = 250, LIMFLAG = 512)
  parameter (LIMPROBS = 50, LIMPROBSETS = 50)
  real*4 graphs(LIMBINS,LIMGRAPHS), areas(LIMBINS,LIMGRAPHS)
  real*4, allocatable :: xModPt(:), yModPt(:), zModPt(:)
  integer*4 iobjFlag(LIMFLAG)
  integer*4 indStart(LIMWOBJ), numPointsObj(LIMWOBJ)
  integer*4 itype(LIMWOBJ)                 !types of sample points
  integer*4 numRefTypes(LIMGRAPHS), numNeighTypes(LIMGRAPHS) !# of types
  integer*4 itypeRef(LIMTYPE,LIMGRAPHS), itypeNeigh(LIMTYPE,LIMGRAPHS)
  integer*4 iwhichEnd(LIMTYPE,LIMGRAPHS), numInClass(LIMTYPE,LIMREGION)
  integer*4 ifFlipRegion(LIMREGION)
  integer*4 numGapsRegion(LIMREGION), indGapRegion(LIMREGION)
  real*4 xyScaleRegion(LIMREGION), zScaleRegion(LIMREGION)
  real*4 xOffsetRegion(LIMREGION), yOffsetRegion(LIMREGION)
  real*4 zStartRegion(LIMREGION), zEndRegion(LIMREGION)
  real*4 xyMaskRegion(4,LIMREGION), zOffsetRegion(LIMREGION)
  real*4 zGapStart(10*LIMREGION), zGapEnd(10*LIMREGION)
  character*320 modelFile, modelRegion(LIMREGION)
  character*320 tiltFile, tiltRegion(LIMREGION)
  logical forceLoad, graphEach, shuffled, onlyShifted, nearestOnly
  logical sampled, converted

  real*4 xmaxDisplay(4), ymaxDisplay(4), changeFrac(LIMTYPE)
  integer*4 igraphDisplay(100), itypeFrom(LIMTYPE), itypeTo(LIMTYPE)
  integer*4 itypeShift(LIMTYPE), itypeCheck(LIMTYPE)
  integer*4 numBinsByGrf(LIMGRAPHS), iobjWin(LIMWOBJ), iobjMod(LIMWOBJ)
  real*4 deltaRadByGrf(LIMGRAPHS), powerGraph(LIMGRAPHS)
  real*4 endSeparation(LIMWOBJ), xyzEnd(3,LIMWOBJ)
  real*4 probNear(LIMPROBS,LIMPROBSETS), deltaNear(LIMPROBSETS)
  integer*4 numNearBins(LIMPROBSETS), iuseProb(LIMPROBSETS)
  integer*4 numRegions, ifFlip, numGaps, indGap, ifShuffle, ifSample, numWobj
  integer*4 ishiftFlag, numCheckType, iobjBound, ifCheckUnshifted, maxTrials
  integer*4 numTrialCycle, ifScatSurf, manyRandom, ifExcludEout, ifConvert, numChange
  integer*4 ifBundEnd, numBins, numGraphs, limPtsFit, numInWin, numObjInWin, ifCloseg
  integer*4 irefFlag, neighFlag, maxGraph, numShiftTypes
  real*4 zScale, xyScale, xOffset, yOffset, zOffset, zStart, zEnd, ranMin, ranMax, boxTol
  real*4 ranZrel, cycleMaxShiftFac, deltaRad, power, winMin, winMax, sampleLen, padBound
  real*4 fracOmit
end module mtkredovars

subroutine realGraphicsMain(comFileArg)
  use mtkredovars
  implicit none
  character*(*) comFileArg
  integer itypeAll, LIMRAND, nOptNeedModel
  parameter (itypeAll = 999)
  parameter (LIMRAND = 1000)
  parameter (nOptNeedModel = 10)
  !
  integer*4 numInTmp(LIMREGION)
  integer*4 iobjRegion(LIMREGION), numPointsRegion(LIMREGION)
  integer*4 itypeCrossInd(-LIMFLAG:LIMFLAG)
  real*4 xyMask(4)
  !
  real*4 probSave(LIMBINS), dum(LIMGRAPHS)
  integer*4 listExtraGraph(LIMGRAPHS), idum(LIMGRAPHS)
  !
  integer*4 ibaseStart(LIMGRAPHS), ibaseEnd(LIMGRAPHS)
  integer*4 integStart(LIMGRAPHS), integEnd(LIMGRAPHS), numRandAbove(LIMRAND)
  real*4 baseline(LIMGRAPHS), realInteg(LIMGRAPHS), sumInteg(LIMRAND), sumSqInteg(LIMRAND)
  integer*4 ioptNeedModel(nOptNeedModel) / 15, 16, 19, 26, 20, 21, 22, 27, 23, 40/
  integer*4 i, ierr, ifAnyPlot, numExtraGraph, numNearBinsSave, ibinAvgStart, ibinAvgEnd
  integer*4 ibinStart, ibinEnd, iwin, iout, iopt, ifNoPlax, numTypes, ifMaskScaled
  integer*4 ity, iclass, ii, ibinTypeStart, ibinTypeEnd, ifRaw, jgrf, lineStart, lineEnd
  integer*4 line, ib1, ib2, ib, nbinAvg, ibin, jj, numDisp, ifxy, iwn, iplot, ifPage
  integer*4 iregionStart, jgraphBase, numInTot, ir, itcr, iow, ireg, ifDoRand, ifNewConv
  integer*4 ityp, notOnList, j, numProbSets, ifUseSave, jgraphAdd, ifAllSame, ifDoMeanSd
  integer*4 maxTmp, numTotControl, numControlDo, ifReally,  icont, jgrMean, jgrSd, iunit
  integer*4 igraphExtra, numBinsHist, numOut, numAtMax, numIn, indSt, ind, iobjShift
  real*4 baseVal, avgDensity, sdDen, semDen, areaSum, countSum, trueAvg, sum, centroid
  real*4 radMax, yMax, yAllMax, dyTmp, yLowTmp, c1, c2, c3, powerNew, rr2, facOld, facNew
  real*4 deltaNearSave, base, randInteg, binSum, binSumSq, avgIntegral, sdIntegral
  real*4 pctAbove, deltaHist, sumsq, dist, connectAvg, connectSd

  logical OptionNeedsModel
  logical checkGrf, checkExtra
  character*40 objName
  integer*4 getImodObjName
  integer*4 in5
  common /nmsInput/ in5
  !
  call setExitPrefix('ERROR: mtk - ')
  allocate(xModPt(LIMXYZ), yModPt(LIMXYZ), zModPt(LIMXYZ), stat = ierr)
  call memoryErrorUC(ierr, 'ARRAYS FOR POINTS')
  call allocateMtkvars()
  in5 = 5
  ifAnyPlot = 0
  numExtraGraph = 0
  winMin = 0.
  winMax = 0.
  manyRandom = 0
  ifExcludEout = 1
  onlyShifted = .false.
  nearestOnly = .false.
  numNearBinsSave = 0
  ibinAvgStart = 1
  ibinAvgEnd = 1
  ibinStart = 1
  ibinEnd = 1
  baseVal = 0
  padBound = 0.
  iwin = 1

  call opencomfile(comFileArg)
  !
  print *,'Enter name of output file to store density values in ', &
      '(Return for none)'
  read(in5, '(a)') modelFile
  if (modelFile == ' ') then
    iout = 6
  else
    iout = 7
    call dopen(7, modelFile, 'new', 'f')
  endif
  iopt = -1
  !
  write(*,'(1x,a,$)') &
      '0 for graphs in Plax window, 1 to suppress graphs: '
  read(in5,*) ifNoPlax
  call scrnOpen(ifNoPlax)
  !
  write(*,'(1x,a,/,a,$)') '0 for 3D density/closest approach '// &
      'analysis', '  or 1 for ends versus bundle analysis: '
  read(in5,*) ifBundEnd
  !
  modelFile = ' '
  call read_model(modelFile, tiltFile, xyScale, zScale, xOffset, yOffset, zOffset, &
      ifFlip, iobjFlag, LIMFLAG, zGapStart, zGapEnd, numGaps)
  if (modelFile == ' ') go to 38
  indGap = 1
  !
  ! initialize for first region
  !
10 numRegions = 1
  numTypes = 0
  do i = -LIMFLAG, LIMFLAG
    itypeCrossInd(i) = 0
  enddo
  irefFlag = -1
  neighFlag = -1
  !
12 if (modelFile .ne. ' ') then
    write(*,'(1x,a,$)') 'Starting and ending Z to include (0,0 '// &
        'for all; 0,-1 for new model): '
    read(in5,*) zStart, zEnd
  else
    zStart = 0
    zEnd = -1
  endif
  if (zStart == 0. .and. zEnd == -1.) then
    if (numRegions == 1) then
      indGap = 1
    else
      indGap = indGap + numGaps
    endif
    modelFile = ' '
    call read_model(modelFile, tiltFile, xyScale, zScale, xOffset, yOffset, &
        zOffset, ifFlip, iobjFlag, LIMFLAG, zGapStart(indGap), &
        zGapEnd(indGap), numGaps)
    if (modelFile == ' ') go to 40
    if (numRegions == 1) go to 12
    call checkFlags(numGraphs, itypeRef, numRefTypes, itypeNeigh, &
        numNeighTypes, iobjFlag, LIMFLAG, irefFlag, neighFlag)
    if (irefFlag >= 0) go to 12
    print *,'Previous regions are being discarded'
    go to 10
  endif
  !
  if (ifBundEnd .ne. 0) then
    xyMask(1) = -1.e10
    xyMask(2) = 1.e10
    xyMask(3) = -1.e10
    xyMask(4) = 1.e10
    ifMaskScaled = 1
    write(*,'(1x,a,/,a,/,a,$)') 'Enter lower & upper limits of X, ' &
        //'and lower & upper limits of Y,', '   within which to'// &
        ' include points for analysis;', '   and 0 if in pixels'// &
        ' or 1 if in microns (/ for no limits): '
    read(in5,*) (xyMask(i), i = 1, 4), ifMaskScaled
    do i = 1, 4
      if (ifMaskScaled == 0) xyMask(i) = xyScale * xyMask(i)
      xyMaskRegion(i, numRegions) = xyMask(i)
    enddo
  endif
  !
  call get_objects(zStart, zEnd, xModPt, yModPt, zModPt, indStart, numPointsObj, &
      itype, numWobj, iobjMod, iobjFlag, LIMXYZ, LIMWOBJ)
  modelRegion(numRegions) = modelFile
  tiltRegion(numRegions) = tiltFile
  zScaleRegion(numRegions) = zScale
  xyScaleRegion(numRegions) = xyScale
  xOffsetRegion(numRegions) = xOffset
  yOffsetRegion(numRegions) = yOffset
  zOffsetRegion(numRegions) = zOffset
  ifFlipRegion(numRegions) = ifFlip
  numPointsRegion(numRegions) = numWobj
  zStartRegion(numRegions) = zStart
  zEndRegion(numRegions) = zEnd
  numGapsRegion(numRegions) = numGaps
  indGapRegion(numRegions) = indGap
  !
  ! count the types
  !
  do i = 1, LIMTYPE
    numInClass(i, numRegions) = 0
  enddo
  do i = 1, numWobj
    ity = itype(i)
    iclass = itypeCrossInd(ity)
    if (iclass == 0) then
      numTypes = numTypes + 1
      iclass = numTypes
      itypeCrossInd(ity) = numTypes
    endif
    if (iobjFlag(abs(ity)) == 1) then
      numInClass(iclass, numRegions) = numInClass(iclass, numRegions) + 1
    else
      numInClass(iclass, numRegions) = numInClass(iclass, numRegions) + numPointsObj(i)
    endif
  enddo
  !
  write(*,*) 'Object   kind    number     name'
  do ity = -LIMFLAG, LIMFLAG
    i = itypeCrossInd(ity)
    objName = ' '
    if (i > 0) then
      ierr = getImodObjName(abs(ity), objName)
      if (iobjFlag(abs(ity)) == 1) then
        write(*,'(i5,2x,a,i8,4x,a)') ity, '  lines', numInClass(i, numRegions), &
            objName
      else
        write(*,'(i5,2x,a,i8,4x,a)') ity, ' points', numInClass(i, numRegions), &
            objName
      endif
    endif
    if (ity > 0 .and. iobjFlag(max(1, ity)) == 4) then
      ierr = getImodObjName(abs(ity), objName)
      write(*,'(i5,2x,a,12x,a)') ity, ' meshes', objName
    endif
  enddo
  write(*,*)
  !
  if (numRegions > 1) go to 35
  !
  call getBinSpec(deltaRad, numBins, power, limPtsFit, padBound, fracOmit, &
      ifBundEnd, sampleLen, ifCloseg, ifScatSurf)
  if (padBound < 0.) padBound = -xyScale * padBound
  !
  call getGraphSpec(numGraphs, itypeRef, numRefTypes, itypeNeigh, numNeighTypes, &
      iwhichEnd, ifBundEnd, iobjFlag, LIMFLAG, irefFlag, neighFlag)
  maxGraph = numGraphs
  !
  !
35 if (ifBundEnd == 0) then
    call closedist(xModPt, yModPt, zModPt, indStart, numPointsObj, itype, numWobj, &
        deltaRad, numBins, numGraphs, numRefTypes, numNeighTypes, itypeRef, &
        itypeNeigh, power, limPtsFit, winMin, winMax, numInWin, graphs, areas, &
        iobjWin, numObjInWin, iobjMod, xyzEnd, endSeparation, sampleLen, &
        ifCloseg, ifScatSurf, irefFlag, neighFlag, xyScale, &
        zScale, powerGraph, zGapStart(indGap), &
        zGapEnd(indGap), numGaps, manyRandom, onlyShifted, nearestOnly)
  else
    call bundledist(xModPt, yModPt, zModPt, indStart, numPointsObj, itype, numWobj, &
        deltaRad, numBins, numGraphs, numRefTypes, numNeighTypes, itypeRef, itypeNeigh, &
        iwhichEnd, limPtsFit, padBound, fracOmit, graphs, areas, zScale, &
        zStart, zEnd, xyMask, winMin, winMax, iobjWin, numObjInWin, xyzEnd)
  endif
  do ii = 1, maxGraph
    ! powergrf(ii) =power
    numBinsByGrf(ii) = numBins
    deltaRadByGrf(ii) = deltaRad
  enddo
  sampled = .false.
  shuffled = .false.
  converted = .false.
  !
  ! if doing multiple regions, average these graphs into the ones
  ! for previous regions.
  !
  if (numRegions > 1) &
      call addToAvg(graphs, areas, LIMBINS, numGraphs, numBins, numRegions)
  !
  ! display up to the first 4 new graphs
  !
  call fourDsp(graphs, LIMBINS, 1, min(4, numGraphs), numBinsByGrf, deltaRadByGrf, &
      xmaxDisplay, ymaxDisplay, igraphDisplay)
  !
  if (iopt .ne. -1) go to 40
38 write(*,104)
104 format(' 1/2: Type/Average selected bins of the graph in a', &
      ' specified window',/,' 3: Compute integrated number', &
      ' of (excess/missing) items in selected bins',/, &
      ' 4/5: Display one graph in a window/Enter list of graphs', &
      ' to display',/, &
      ' 6/7: Rescale X or Y axis of one window/Y axis of all', &
      ' windows',/, ' 8/9: Plot one window/all windows to', &
      ' PostScript graphics file',/, &
      ' 10/11: Output PostScript file to screen/printer',/, &
      ' 12: Output single or average graph to file      45: Set PostScript filename',/, &
      ' 13: Loop back to specify new range of Z to analyze', &
      ' (or new model)',/, ' 14: Change radial weighting of', &
      ' a graph',/, &
      ' 15: Analyze new region and average with ', &
      'previous region(s)',/, ' 16: Redo current region(s) with', &
      ' new bin size, # of bins, or types for graphs',/, &
      ' 17: Set min & max distances at which to compute angles', &
      ' and add lines to model ',/, ' 18: Save bins of a graph' &
      ,' to specify rejection probabilities for random points',/, &
      ' 19/26/20: Do current region(s) with shuffled/converted', &
      ' types or random shifts',/,' 21: Save current set of ', &
      'objects and their types as an IMOD model',/, &
      ' 22/27/23: Do many sets & integrals with ', &
      'shuffled/converted types/random shifts', &
      /,' 24: Take command input from file        25: Exit', &
      /,' 28/29/30 Save a graph/Average/Combine 2 graphs into', &
      ' an extra graph location',/,' 31/32: Save graph in', &
      ' file/Read from file into an extra graph location',/, &
      ' 33: Replace some sets of bins by their averages' &
      ,/,' 37/38/39 Add list of graphs/Read', &
      ' list of graphs from file/Read&Add from file',/, &
      ' 40: Unshift an object',/, &
      ' 41: Toggle between including and excluding items that ', &
      'failed to shift',/, &
      ' 42: Export graph values or points for drawing to file',/, &
      ' 43: List distances of close approach between min/max limits',/, &
      ' 44: Toggle between recording distances to all and nearest neighbors')
  !
40 write(*,'(1x,a,$)') 'Option, or -1 for list of choices: '
  read(in5,*,err = 40, end = 225) iopt
  if (iopt == -1) go to 38
  if (OptionNeedsModel(modelFile, ioptNeedModel, &
      nOptNeedModel, iopt)) go to 40
  go to(201, 202, 203, 204, 205, 206, 207, 208, 209, 210, 210, 212, 213, 214, &
      215, 216, 217, 218, 219, 220, 221, 222, 222, 224, 225, 226, 222, 228, &
      228, 228, 228, 228, 228, 213, 235, 235, 228, 228, 228, 240, 241, 228, 243, 244,  &
      245), iopt
  go to 40
  !
  ! type bins
  !
201 ibinTypeStart = 1
  ibinTypeEnd = -1
  write(*,'(1x,a,$)') 'Window # (- for raw counts), starting,' &
      //' ending bins to type (/ for all): '
  read(in5,*) iwin, ibinTypeStart, ibinTypeEnd
  ifRaw = iwin
  iwin = abs(iwin)
  if (iwin <= 0 .or. iwin > 4) go to 40
  if (igraphDisplay(iwin) == 0) go to 40
  jgrf = igraphDisplay(iwin)
  if (ibinTypeEnd == -1) ibinTypeEnd = numBinsByGrf(jgrf)
  ibinTypeEnd = min(ibinTypeEnd, numBinsByGrf(jgrf))
  ibinTypeStart = max(1, ibinTypeStart)
  lineStart = 1 + (ibinTypeStart - 1) / 5
  lineEnd = 1 + (ibinTypeEnd-1) / 5
  write(*,107)
107 format(12x,'1',14x,'2',14x,'3',14x,'4',14x,'5')
  do line = lineStart, lineEnd
    ib1 = 5 * (line - 1)
    ib2 = min(ibinTypeEnd, ib1 + 5)
    if (ifRaw > 0) then
      write(*,108) ib1, (graphs(ib, jgrf), ib = ib1 + 1, ib2)
108   format(i4,'/',f14.8,4f15.8)
    else
      write(*,1108) ib1, (nint(graphs(ib, jgrf) * areas(ib, jgrf)) &
          , ib = ib1 + 1, ib2)
1108  format(i4,'/',i7,4i15)
    endif
  enddo
  go to 40
  !
  ! average bins
  !
202 write(*,'(1x,a,i2,2i4,a$)') 'Window #, starting, ending bins'// &
      ' to average [', iwin, ibinAvgStart, ibinAvgEnd, ']: '
  read(in5,*) iwin, ibinAvgStart, ibinAvgEnd
  if (iwin <= 0 .or. iwin > 4) go to 40
  if (igraphDisplay(iwin) == 0) go to 40
  jgrf = igraphDisplay(iwin)
  ibinAvgEnd = min(ibinAvgEnd, numBinsByGrf(jgrf))
  ibinAvgStart = max(1, ibinAvgStart)
  nbinAvg = ibinAvgEnd + 1 - ibinAvgStart
  call avgsd(graphs(ibinAvgStart, jgrf), nbinAvg, avgDensity, sdDen, semDen)
  write(*,109) nbinAvg, avgDensity, sdDen, semDen
109 format(i4,' bins, mean =',f15.8,', SD =',f15.9,', SEM =',f15.9)
  baseVal = avgDensity
  areaSum = 0.
  countSum = 0.
  do ibin = ibinAvgStart, ibinAvgEnd
    areaSum = areaSum + areas(ibin, jgrf)
    countSum = countSum + areas(ibin, jgrf) * graphs(ibin, jgrf)
  enddo
  trueAvg = 0.
  if (areaSum > 0.) trueAvg = countSum / areaSum
  write(*,1109) trueAvg
1109 format(' True area-weighted average =',f15.8)
  go to 40
  !
  ! integrate bins
  !
203 write(*,'(1x,a,/,a,i2,2i4,f15.8,a,$)') 'Enter window #, starting' &
      //' & ending bins to integrate,', '     and base value '// &
      'to subtract [', iwin, ibinStart, ibinEnd, baseVal, ']: '
  read(in5,*) iwin, ibinStart, ibinEnd, baseVal
  if (iwin <= 0 .or. iwin > 4) go to 40
  if (igraphDisplay(iwin) == 0) go to 40
  jgrf = igraphDisplay(iwin)
  ibinEnd = min(ibinEnd, numBinsByGrf(jgrf))
  ibinStart = max(1, ibinStart)
  call integrate(graphs(1, jgrf), areas(1, jgrf), numBinsByGrf(jgrf), &
      deltaRadByGrf(jgrf), powerGraph(jgrf), ibinStart, ibinEnd, 0, 0, baseVal, sum, &
      centroid)
  write(*,110) ibinEnd + 1 - ibinStart, (ibinStart - 1) * deltaRadByGrf(jgrf), &
      ibinEnd * deltaRadByGrf(jgrf), sum, centroid
110 format(' For the',i4,' bins from', &
      f10.3,' to',f10.3,/,5x, &
      'the integrated number of (excess/missing) items is',f10.5/,5x, &
      'the centroid of their distances is', f11.4)
  go to 40
  !
  ! Display graph in window
  !
204 write(*,'(1x,a,$)') &
      'Display graph in window; Enter graph # and window #: '
  read(in5,*) jgrf, iwin
  if (iwin <= 0 .or. iwin > 4 .or. &
      .not.checkGrf(jgrf, maxGraph, numExtraGraph, listExtraGraph)) go to 40
  xmaxDisplay(iwin) = -1.
  ymaxDisplay(iwin) = -1.
  igraphDisplay(iwin) = jgrf
  call graphDsp(graphs(1, jgrf), numBinsByGrf(jgrf), deltaRadByGrf(jgrf), iwin, &
      jgrf, xmaxDisplay(iwin), ymaxDisplay(iwin))
  go to 40
  !
  ! display list of graphs in windows
  !
205 write(*,'(1x,a,$)') 'List of graphs to display (ranges OK): '
  do jj = 1, 4
    igraphDisplay(jj) = 0
  enddo
  call rdlist(in5, igraphDisplay, numDisp)
  if (numDisp == 0) go to 40
  do jj = 1, min(4, numDisp)
    if (.not.checkGrf(igraphDisplay(jj), maxGraph, numExtraGraph, listExtraGraph)) &
        go to 205
  enddo
  call someDsp(graphs, LIMBINS, numBinsByGrf, deltaRadByGrf, xmaxDisplay, ymaxDisplay, &
      igraphDisplay)
  go to 40
  !
  ! Rescale X or Y axis of graph in window
  !
206 write(*,'(1x,a,$)') &
      'Number of window to rescale; 0 or 1 to rescale X or Y: '
  read(in5,*) iwin, ifxy
  if (iwin <= 0 .or. iwin > 4) go to 40
  jgrf = igraphDisplay(iwin)
  if (jgrf == 0) go to 40
  if (ifxy == 0) then
    radMax = numBinsByGrf(jgrf) * deltaRadByGrf(jgrf)
    write(*,106) ' X', radMax, xmaxDisplay(iwin)
106 format(a,' ranges up to',f13.6,' and current full scale', &
        ' value is',f13.6,/,'   Enter new full scale value: ',$)
    read(in5,*) xmaxDisplay(iwin)
  else
    yMax = 0.
    do i = 1, numBinsByGrf(jgrf)
      yMax = max(yMax, graphs(i, igraphDisplay(iwin)))
    enddo
    write(*,106) ' Y', yMax, ymaxDisplay(iwin)
    read(in5,*) ymaxDisplay(iwin)
  endif
  call graphDsp(graphs(1, jgrf), numBinsByGrf(jgrf), deltaRadByGrf(jgrf), iwin, &
      jgrf, xmaxDisplay(iwin), ymaxDisplay(iwin))
  go to 40
  !
  ! rescale all graphs to same maximum
  !
207 yAllMax = 0
  do iwn = 1, 4
    jgrf = igraphDisplay(iwn)
    if (jgrf > 0) then
      do ibin = 1, numBinsByGrf(jgrf)
        yAllMax = max(yAllMax, graphs(ibin, jgrf))
      enddo
    endif
    call scale(0, 0.9999 * yAllMax, dyTmp, yLowTmp)
    yAllMax = 10. *dyTmp
  enddo
  do iwn = 1, 4
    jgrf = igraphDisplay(iwn)
    if (jgrf > 0) then
      ymaxDisplay(iwn) = yAllMax
      call graphDsp(graphs(1, jgrf), numBinsByGrf(jgrf), deltaRadByGrf(jgrf), &
          iwn, jgrf, xmaxDisplay(iwn), ymaxDisplay(iwn))
    endif
  enddo
  go to 40
  !
  ! Plot one window to metacode file
  !
208 write(*,'(1x,a,$)') &
      'Window number to plot, plot number or 0 to specify plot: '
  read(in5,*) iwin, iplot
  if (iwin <= 0 .or. iwin > 4) go to 40
  jgrf = igraphDisplay(iwin)
  if (jgrf == 0) go to 40
  if (ifAnyPlot .ne. 0) then
    write(*,2081)
2081 format(' 0 for plot on same page as previous plot(s),', &
        ' 1 for new page: ',$)
    read(in5,*) ifPage
    call psSetup(1, c1, c2, c3, 0)
    if (ifPage .ne. 0) call psFrame
  endif
  call graphPlt(graphs(1, jgrf), numBinsByGrf(jgrf), deltaRadByGrf(jgrf), iplot, &
      jgrf, xmaxDisplay(iwin), ymaxDisplay(iwin))
  ifAnyPlot = 1
  go to 40
  !
  ! Quick plot all 4 windows
  !
209 if (ifAnyPlot .ne. 0) then
    write(*,2081)
    read(in5,*) ifPage
    call psSetup(1, c1, c2, c3, 0)
    if (ifPage .ne. 0) call psFrame
  endif
  do iwn = 1, 4
    jgrf = igraphDisplay(iwn)
    if (jgrf > 0) then
      call graphPlt(graphs(1, jgrf), numBinsByGrf(jgrf), deltaRadByGrf(jgrf), &
          iwn, jgrf, xmaxDisplay(iwn), ymaxDisplay(iwn))
    endif
  enddo
  ifAnyPlot = 1
  go to 40
  !
  ! Set filename for graphics output
245 write(*,'(1x,a)') 'Enter name of file for PostScript graphics output'
  read(in5, '(a)') modelFile
  call psSetFilename(modelFile)
  go to 40
  !
  ! Output plot file to workstation or printer
  !
210 call pltout(11 - iopt)
  ifAnyPlot = 0
  go to 40
  !
  ! write graph to file with lots of info
  !
212 write(*,'(1x,a,$)') 'Graph #: '
  read(in5,*) jgrf
  if (.not.checkGrf(jgrf, maxGraph, numExtraGraph, listExtraGraph)) go to 40
  if (jgrf <= numGraphs) then
    iregionStart = numRegions
    jgraphBase = jgrf
  else
    iregionStart = 1
    jgraphBase = jgrf - numGraphs
  endif
  write(iout, 111) modelFile
111 format(' Model file:  ',a)
  numInTot = 0
  do ir = iregionStart, numRegions
    numInTmp(ir) = 0
    do ity = 1, numRefTypes(jgraphBase)
      if (itypeRef(ity, jgraphBase) == itypeAll) then
        numInTmp(ir) = numPointsRegion(ir)
        go to 362
      endif
      itcr = itypeCrossInd(itypeRef(ity, jgraphBase))
      if (itcr .ne. 0) numInTmp(ir) = numInTmp(ir) + numInClass(itcr, ir)
    enddo
362 numInTot = numInTot + numInTmp(ir)
  enddo
  write(iout, 112) ' Reference items:', (numInTmp(i), i = iregionStart, numRegions)
112 format(a,(t20,6i10))
  if (iregionStart .ne. numRegions) write(iout, 112) '         Total:', numInTot
  !
  numInTot = 0
  do ir = iregionStart, numRegions
    numInTmp(ir) = 0
    do ity = 1, numNeighTypes(jgraphBase)
      if (itypeNeigh(ity, jgraphBase) == itypeAll) then
        numInTmp(ir) = numPointsRegion(ir)
        go to 364
      endif
      itcr = itypeCrossInd(itypeNeigh(ity, jgraphBase))
      if (itcr .ne. 0) numInTmp(ir) = numInTmp(ir) + numInClass(itcr, ir)
    enddo
364 numInTot = numInTot + numInTmp(ir)
  enddo
  write(iout, 112) ' Neighboring items:' &
      , (numInTmp(i), i = iregionStart, numRegions)
  if (iregionStart .ne. numRegions) write(iout, 112) '         Total:', numInTot
  !
  write(iout, 114) ' Reference objects:', &
      (itypeRef(ity, jgraphBase), ity = 1, numRefTypes(jgraphBase))
114 format(a,(12i5))
  write(iout, 114) ' Neighboring objects:', &
      (itypeNeigh(ity, jgraphBase), ity = 1, numNeighTypes(jgraphBase))
  write(iout, '(a,f12.5)') ' Delta r =', deltaRadByGrf(jgrf)
  write(iout, '(5g15.8)') (graphs(i, jgrf), i = 1, numBinsByGrf(jgrf))
  go to 40
  !
  ! Specify new boundary object: but check first if averaging
  ! 8/5/03: remove confirmation to make it easier to do command files
  !
!!$    213     if (nregion>1) then
!!$    write(*,'(1x,a,$)') 'This will destroy the stored average.'// &
!!$                  '  Enter 1 to do so: '
!!$    read(in5,*) ifreally
!!$    if (ifreally .ne. 1) go to 40
!!$    endif
213 if (iopt == 13) go to 10
  if (numObjInWin > 0 .and. winMin < winMax) then
    do iow = 1, numObjInWin
      if (iobjWin(iow) > 0) iobjWin(iow) = iobjMod(iobjWin(iow))
    enddo
    numObjInWin = -numObjInWin
  endif
  if (winMin < winMax) then
    write(*,'(1x,a,$)') 'New min and max distances: '
    read(in5,*) winMin, winMax
  endif
  ifBundEnd = max(1, ifBundEnd) - ifBundEnd
  go to 10
  !
  ! change power of a graph in a window
  !
214 write(*,'(1x,a,$)') 'Change power of graph; Enter window #'// &
      ' and new power: '
  read(in5,*) iwin, powerNew
  if (iwin <= 0 .or. iwin > 4) go to 40
  jgrf = igraphDisplay(iwin)
  if (jgrf == 0) go to 40
  xmaxDisplay(iwin) = -1.
  ymaxDisplay(iwin) = -1.
  do ibin = 1, numBinsByGrf(jgrf)
    rr2 = 2. *deltaRadByGrf(jgrf) * (ibin - 0.5)
    facOld = 1.
    if (powerGraph(jgrf) .ne. 0.) facOld = rr2**powerGraph(jgrf)
    facNew = 1.
    if (powerNew .ne. 0.) facNew = rr2**powerNew
    areas(ibin, jgrf) = areas(ibin, jgrf) * facNew / facOld
    graphs(ibin, jgrf) = graphs(ibin, jgrf) * facOld / facNew
  enddo
  powerGraph(jgrf) = powerNew
  call graphDsp(graphs(1, jgrf), numBinsByGrf(jgrf), deltaRadByGrf(jgrf), iwin, &
      jgrf, xmaxDisplay(iwin), ymaxDisplay(iwin))
  go to 40
  !
  ! Specify new region to analyze and average with previous region(s)
  !
215 if (2 * numGraphs > LIMGRAPHS) then
    print *,'Too many graphs to accumulate averages'
    go to 40
  endif
  if (sampled .or. shuffled .or. converted) then
    print *,'Data have been altered; need to redo previous'// &
        ' regions before proceeding...'
    forceLoad = .true.
    ifShuffle = 0
    ifSample = 0
    ifConvert = 0
    graphEach = .true.
    call redoRegions()
  endif
  if (numRegions == 1) then
    write(*,'(a,i2,a,i2)') ' The averages will be'// &
        ' stored in graphs ', numGraphs + 1, ' to ', 2 * numGraphs
    call addToAvg(graphs, areas, LIMBINS, numGraphs, numBins, numRegions)
  endif
  maxGraph = 2 * numGraphs
  numRegions = numRegions + 1
  go to 12
  !
  ! do other kind of histogram on region(s)
  !
217 write(*,'(1x,a,$)') 'Min & max distance for adding connecting' &
      //' lines and computing angles: '
  read(in5,*) winMin, winMax
  !
  ! Redo region(s) with new bins and/or graphs
  !
216 call getBinSpec(deltaRad, numBins, power, limPtsFit, padBound, fracOmit, &
      ifBundEnd, sampleLen, ifCloseg, ifScatSurf)
  if (padBound < 0.) padBound = -xyScale * padBound
  !
2161 call getGraphSpec(numGraphs, itypeRef, numRefTypes, itypeNeigh, numNeighTypes, &
      iwhichEnd, ifBundEnd, iobjFlag, LIMFLAG, irefFlag, neighFlag)
  !
  ! check legality with any other model files
  !
  do ireg = 1, numRegions
    if (modelFile .ne. modelRegion(ireg)) then
      modelFile = modelRegion(ireg)
      tiltFile = tiltRegion(ireg)
      zScale = zScaleRegion(ireg)
      xyScale = xyScaleRegion(ireg)
      xOffset = xOffsetRegion(ireg)
      yOffset = yOffsetRegion(ireg)
      zOffset = zOffsetRegion(ireg)
      ifFlip = ifFlipRegion(ireg)
      numGaps = numGapsRegion(ireg)
      indGap = indGapRegion(ireg)
      call read_model(modelFile, tiltFile, xyScale, zScale, xOffset, yOffset, &
          zOffset, ifFlip, iobjFlag, LIMFLAG, zGapStart(indGap), &
          zGapEnd(indGap), numGaps)
      call checkFlags(numGraphs, itypeRef, numRefTypes, itypeNeigh, &
          numNeighTypes, iobjFlag, LIMFLAG, irefFlag, neighFlag)
      if (irefFlag < 0) then
        print *,'You need to enter graph specifications that ', &
            'will work for all models'
        go to 2161
      endif
    endif
  enddo
  maxGraph = numGraphs * min(2, numRegions)
  forceLoad = numRegions == 1 .and. (shuffled .or. sampled .or. converted)
  ifShuffle = 0
  ifSample = 0
  ifConvert = 0
  !
  ! general place to redo regions,
  !
315 if (forceLoad) then
    write(*,'(1x,a,/,a,$)') 'The data to be analyzed have been'// &
        ' modified/randomized.', ' Enter 1 to work with the'// &
        ' altered data or 0 to reload original data: '
    read(in5,*) ifDoRand
    forceLoad = ifDoRand == 0
  endif
  graphEach = .true.
  call redoRegions()
  if (numRegions > 1) call fourDsp(graphs, LIMBINS, numGraphs + 1, numGraphs &
      + min(4, numGraphs), numBinsByGrf, deltaRadByGrf, xmaxDisplay, ymaxDisplay, &
      igraphDisplay)
  go to 40
  !
  ! save bins of some graph to specify restriction on distances
  !
218 write(*,'(1x,a,$)') &
      'Graph #, baseline level corresponding to probability 1.0: '
  read(in5,*) jgrf, base
  if (.not.checkGrf(jgrf, maxGraph, numExtraGraph, listExtraGraph)) go to 40
  numNearBinsSave = 0
  deltaNearSave = deltaRadByGrf(jgrf)
  do while(numNearBinsSave < numBinsByGrf(jgrf) .and. &
      graphs(numNearBinsSave + 1, jgrf) < base)
    numNearBinsSave = numNearBinsSave + 1
    probSave(numNearBinsSave) = graphs(numNearBinsSave, jgrf) / base
  enddo
  if (numNearBinsSave > LIMPROBS) then
    print *,numNearBinsSave, &
        ' bins would be saved - too many for the arrays'
    numNearBinsSave = 0
  else
    print *,numNearBinsSave, ' probability bins saved'
  endif
  go to 40
  !
  ! shuffle types
  !
219 ifShuffle = 1
  forceLoad = numRegions == 1 .and. (sampled .or. converted)
  ifSample = 0
  ifConvert = 0
  if (iopt == 19) go to 315
  go to 319
  !
  ! convert types
  !
226 ifShuffle = 0
  forceLoad = numRegions == 1 .and. (sampled .or. converted .or. shuffled)
  ifSample = 0
  ifConvert = 1
  if (numChange > 0) then
    write(*,'(1x,a,$)') &
        '0 to use last conversions, 1 to specify new ones: '
    read(in5,*) ifNewConv
    if (ifNewConv .ne. 0) numChange = 0
  endif
  if (numChange == 0) then
    write(*,'(1x,a,$)') 'Number of objects to convert: '
    read(in5,*) numChange
    print *,'For each conversion, enter the object to change', &
        ' from, the object to change to,', &
        ' and the fraction of contours of that object to convert.'
    do i = 1, numChange
      write(*,'(1x,a,i3,a$)') 'Conversion', i, ': '
      read(in5,*) itypeFrom(i), itypeTo(i), changeFrac(i)
    enddo
  endif
  if (iopt == 26) go to 315
  go to 319
  !
  ! do series of control sets: first see if need to recompute graphs
  !
222 if (shuffled .or. sampled .or. converted) then
    print *,'Data were modified/randomized; must rebuild '// &
        'actual graphs to get integrals...'
    ifShuffle = 0
    ifSample = 0
    ifConvert = 0
    forceLoad = .true.
    graphEach = .false.
    call redoRegions()
  endif
  if (iopt == 22) go to 219
  if (iopt == 27) go to 226
  !
  ! shift objects randomly in X/Y
  !
220 write(*,'(1x,a,$)') 'Minimum and maximum distance to shift'// &
      'in X/Y plane: '
  read(in5,*) ranMin, ranMax
  !
  write(*,'(1x,a,$)') 'Maximum amount to shift in Z relative ' &
      //'to X/Y shift (0 for no Z shift): '
  read(in5,*) ranZrel
  !
  print *,'Enter list of objects to shift', &
      ' (Return for all)'
  call rdlist(in5, itypeShift, numShiftTypes)
  if (numShiftTypes == 0) then
    itypeShift(1) = itypeAll
    do i = 1, LIMFLAG
      if (iobjFlag(i) >= 0) then
        numShiftTypes = numShiftTypes + 1
        itypeShift(numShiftTypes) = i
      endif
    enddo
  endif
  ishiftFlag = -1
  do i = 1, numShiftTypes
    ityp = abs(itypeShift(i))
    if (iobjFlag(ityp) >= 0) then
      if (ishiftFlag < 0) then
        ishiftFlag = iobjFlag(ityp)
      elseif (iobjFlag(ityp) .ne. ishiftFlag) then
        print *,'Object types do not all match'
        go to 40
      endif
    endif
  enddo
  !
  print *,'Enter list of other objects to check'// &
      ' distances from' , ' (Return for all other)'
  call rdlist(in5, itypeCheck, numCheckType)
  if (numCheckType == 0) then
    itypeCheck(1) = itypeAll
    do i = 1, LIMFLAG
      if (iobjFlag(i) >= 0) then
        notOnList = 1
        do j = 1, numShiftTypes
          if (itypeShift(j) == i) notOnList = 0
        enddo
        if (notOnList == 1) then
          numCheckType = numCheckType + 1
          itypeCheck(numCheckType) = i
        endif
      endif
    enddo
  endif
  do i = 1, numShiftTypes
    do j = 1, numCheckType
      if (itypeShift(i) == itypeCheck(j)) then
        print *,itypeShift(i), ' is in both lists'
        go to 40
      endif
    enddo
  enddo
  !
  write(*,'(1x,a,$)') '# of probability curves to use '// &
      'for rejection of close spacings: '
  read(in5,*) numProbSets
  if (numProbSets > LIMPROBSETS) then
    print *,'Too many curves for arrays'
    go to 40
  endif
  if (numProbSets == 0) then
    numNearBins(1) = 0
    numProbSets = 1
  else
    if (numProbSets > 1) print *,'The first set will be ', &
        'applied to the distances between shifted items'
    do j = 1, numProbSets

      ifUseSave = 0
      if (numNearBinsSave > 0 .and. j == 1) then
        write(*,'(1x,a,$)') '1 to use saved probability '// &
            'bins for the first curve, 0 not to: '
        read(in5,*) ifUseSave
      endif
      !
      if (ifUseSave .ne. 0) then
        numNearBins(j) = numNearBinsSave
        deltaNear(j) = deltaNearSave
        do ii = 1, numNearBinsSave
          probNear(ii, j) = probSave(ii)
        enddo
      else
        !
        write(*,'(1x,a,i3,a,$)') 'Curve #', j, &
            ': Enter # of bins for rejection and  bin size: '
        read(in5,*) numNearBins(j), deltaNear(j)
        if (numNearBins(j) > LIMPROBS) then
          print *,'Too many bins for arrays'
          go to 40
        endif
        if (numNearBins(j) > 0) then
          write(*,'(1x,a,$)') 'Probability values: '
          read(in5,*) (probNear(i, j), i = 1, numNearBins(j))
        endif
      endif
    enddo
  endif
  if (numProbSets > 1 .and. numCheckType > 0) then
    write(*,'(1x,a,i3,/,a,$)') 'Enter the # of the curve to '// &
        'use for each of the', numCheckType, &
        'objects being checked against: '
    read(in5,*) (iuseProb(i), i = 1, numCheckType)
    do i = 1, numCheckType
      if (iuseProb(i) < 1 .or. iuseProb(i) > numProbSets) then
        print *,'Illegal curve number'
        go to 40
      endif
    enddo
  else
    do i = 1, numCheckType
      iuseProb(i) = 1
    enddo
  endif
  !
  write(*,'(1x,a,$)') 'Maximum distance to shift outside '// &
      'bounding box of original data: '
  read(in5,*) boxTol
  write(*,'(1x,a,$)') &
      'Object # of object with bounding contours, or 0 if none: '
  read(in5,*) iobjBound
  write(*,'(1x,a,/,a,$)') '1 to check shifted items against'// &
      ' ones yet to be shifted, or 0 to check', &
      ' only against ones that have been shifted already: '
  read(in5,*) ifCheckUnshifted
  write(*,'(1x,a,$)') 'Maximum number of trials: '
  read(in5,*) maxTrials
  write(*,'(1x,a,$)') '# of trials per cycle, factor to change ' &
      //'maximum shift by per cycle: '
  read(in5,*) numTrialCycle, cycleMaxShiftFac
  ifShuffle = 0
  ifSample = 1
  ifConvert = 0
  forceLoad = numRegions == 1 .and. (sampled .or. converted .or. shuffled)
  if (iopt == 20) go to 315
  !
  ! general setup to do series of control runs and gather statistics
  !
319 jgraphAdd = 0
  if (numRegions > 1) jgraphAdd = numGraphs
  write(*,'(1x,a,$)') '0 to specify integral for each graph '// &
      'separately, 1 to use same bins for all: '
  read(in5,*) ifAllSame
  !
  do jj = 1, numGraphs
    if (ifAllSame == 0) write(*,'(10x,a,i3)') &
        'Specify bins and baseline for graph #', jj + jgraphAdd
    !
    if (ifAllSame == 0 .or. jj == 1) then
321   write(*,'(1x,a,$)') 'Starting and ending bins to integrate: '
      read(in5,*) integStart(jj), integEnd(jj)
      if (integEnd(jj) < integStart(jj) .or. integStart(jj) < 1 &
          .or. integEnd(jj) > numBins) go to 321
      !
322   write(*,'(1x,a,$)') 'Start & end bins to compute baseline'// &
          ' from, or 0,0 to use fixed value: '
      read(in5,*) ibaseStart(jj), ibaseEnd(jj)
      if ((ibaseStart(jj) .ne. 0 .or. ibaseEnd(jj) .ne. 0) .and. &
          (ibaseEnd(jj) < ibaseStart(jj) .or. ibaseEnd(jj) > numBins &
          .or. ibaseStart(jj) < 1)) go to 322
      !
      if (ibaseStart(jj) == 0 .and. ibaseEnd(jj) == 0) then
        write(*,'(1x,a,$)') 'Fixed baseline value: '
        read(in5,*) baseline(jj)
      endif
      !
    else
      integStart(jj) = integStart(1)
      integEnd(jj) = integEnd(1)
      ibaseStart(jj) = ibaseStart(1)
      ibaseEnd(jj) = ibaseEnd(1)
      baseline(jj) = baseline(1)
    endif
    !
    numRandAbove(jj) = 0
    sumInteg(jj) = 0
    sumSqInteg(jj) = 0
    call integrate(graphs(1, jj + jgraphAdd), areas(1, jj + jgraphAdd), numBins, &
        deltaRad, powerGraph(jj + jgraphAdd), integStart(jj), integEnd(jj), &
        ibaseStart(jj), ibaseEnd(jj), baseline(jj), realInteg(jj), centroid)
    write(*,116) jj + jgraphAdd, realInteg(jj)
116 format(' Graph #',i3,', real integral =',f10.5)
  enddo
  !
  write(*,'(1x,a,$)') 'Enter 1 to accumulate mean and standard ' &
      //'deviation graphs or 0 not to: '
  read(in5,*) ifDoMeanSd

  maxTmp = jgraphAdd + 3 * numGraphs
  if (ifDoMeanSd .ne. 0) then
    if (maxTmp <= LIMGRAPHS) then
      write(*,'(1x,a,i3,a,i3)') 'Mean and standard deviation' &
          //' graphs will be stored in graphs', &
          jgraphAdd + numGraphs + 1, ' to', maxTmp
    else
      ifDoMeanSd = 0
      print *,'Too many graphs to accumulate mean and SD graphs'
    endif
  endif
  if (ifDoMeanSd .ne. 0) maxGraph = maxTmp
  numTotControl = 0
  !
324 write(*,'(1x,a,$)') &
      'Number of control sets to run, or 0 to enter new option: '
  read(in5,*) numControlDo
  if (numControlDo <= 0) go to 40
  if (numTotControl .ne. 0 .and. numControlDo .ne. 0) then
    write(*,'(1x,a,$)') 'Enter 1 if you really want to do more'// &
        ' sets, or 0 if that was a mistake: '
    read(in5,*) ifReally
    if (ifReally .ne. 1) go to 40
  endif
  print *,' '
  !
  icont = 1
  !
  ! Get graphs from all regions as needed
  !
424 numTotControl = numTotControl + 1
  write(*,'(a,i4,$)') char(13) //'Working on set', numTotControl
  call flush(6)
  ! if (sampled) forceload=nregion==1
  manyRandom = 1
  graphEach = .false.
  call redoRegions()
  !
  ! Get integrals from graphs and accumulate statistics
  !
  manyRandom = 0
  do jj = 1, numGraphs
    call integrate(graphs(1, jj + jgraphAdd), areas(1, jj + jgraphAdd), numBins, &
        deltaRad, powerGraph(jj + jgraphAdd), integStart(jj), integEnd(jj), &
        ibaseStart(jj), ibaseEnd(jj), baseline(jj), randInteg, centroid)
    if (randInteg > realInteg(jj)) &
        numRandAbove(jj) = numRandAbove(jj) + 1
    sumInteg(jj) = sumInteg(jj) + randInteg
    sumSqInteg(jj) = sumSqInteg(jj) + randInteg**2
  enddo
  !
  ! accumulate mean&SD if desired
  !
  if (ifDoMeanSd .ne. 0) then
    do jj = 1 + jgraphAdd, numGraphs + jgraphAdd
      jgrMean = jj + numGraphs
      jgrSd = jgrMean + numGraphs
      do ibin = 1, numBins
        if (numTotControl == 1) then
          graphs(ibin, jgrMean) = graphs(ibin, jj)
          graphs(ibin, jgrSd) = 0.
          areas(ibin, jgrMean) = areas(ibin, jj)
          areas(ibin, jgrSd) = areas(ibin, jj)
          powerGraph(jgrMean) = powerGraph(jj)
          powerGraph(jgrSd) = powerGraph(jj)
        else
          binSum = graphs(ibin, jgrMean) * (numTotControl - 1)
          binSumSq = graphs(ibin, jgrSd) * sqrt(numTotControl - 2.) &
              + binSum**2 / (numTotControl - 1)
          binSum = binSum + graphs(ibin, jj)
          binSumSq = binSumSq + graphs(ibin, jj)**2
          call sums_to_avgsd(binSum, binSumSq, numTotControl, &
              graphs(ibin, jgrMean), graphs(ibin, jgrSd))
          areas(ibin, jgrMean) = areas(ibin, jgrMean) + areas(ibin, jj)
          areas(ibin, jgrSd) = areas(ibin, jgrMean)
        endif
      enddo
    enddo
  else
  endif
  icont = icont + 1
  if (icont <= numControlDo) go to 424
  !
  do iunit = 6, iout
    write(iunit, 117) numTotControl
117 format(/,' Total # of controls:',i6,/,' Graph    real   ', &
        '    random_integrals     higher_than_real',/, &
        '   #    integral      mean     S.D.          #      %')
    do jj = 1, numGraphs
      call sums_to_avgsd(sumInteg(jj), sumSqInteg(jj), numTotControl, &
          avgIntegral, sdIntegral)
      pctAbove = 100. *numRandAbove(jj) / numTotControl
      write(iunit, 118) jj + jgraphAdd, realInteg(jj), avgIntegral, sdIntegral, &
          numRandAbove(jj), pctAbove
118   format(i5,f10.5,f12.5,f10.6,i10,f8.2)
    enddo
  enddo
  go to 324
  !
  ! save current objects as model file
  !
221 call save_model(xyScale, zScale, xOffset, yOffset, zOffset, ifFlip, iobjFlag, &
        xModPt, yModPt, zModPt, indStart, numPointsObj, itype, numWobj, numInWin, &
        iobjWin, numObjInWin, iobjMod, endSeparation)
  call read_model(modelFile, tiltFile, xyScale, zScale, xOffset, yOffset, zOffset, &
      ifFlip, iobjFlag, LIMFLAG, zGapStart(indGap), &
      zGapEnd(indGap), numGaps)
  go to 40
  !
  ! Switch to reading in a command file
  !
224 call opencomfile(' ')
  go to 40
  !
  ! exit
  !
225 call scrnClose
  call psExit
  !
  ! call to manipulate graphs
  !
228 do i = 1, LIMGRAPHS
    dum(i) = 0.
    idum(i) = 0
  enddo
  call manipGraphs(iopt, 'mtk', graphs, areas, numBinsByGrf, deltaRadByGrf, &
      idum, dum, powerGraph, maxGraph, numExtraGraph, listExtraGraph, &
      igraphDisplay, xmaxDisplay, ymaxDisplay)
  go to 40
  !
  ! process end separations
  !
235 if (numObjInWin >= 0) go to 40
  write(*,'(1x,a,$)') &
      'Graph # to place histogram in, bin width, # of bins: '
  read(in5,*) igraphExtra, deltaHist, numBinsHist
  if (numBinsHist <= 0 .or. deltaHist <= 0.) go to 40
  if (.not.checkExtra(igraphExtra, LIMGRAPHS, listExtraGraph, &
      numExtraGraph)) go to 40
  do ibin = 1, numBinsHist
    areas(ibin, igraphExtra) = 1.
    graphs(ibin, igraphExtra) = 0.
  enddo
  numBinsByGrf(igraphExtra) = numBinsHist
  deltaRadByGrf(igraphExtra) = deltaHist
  powerGraph(igraphExtra) = 1.
  if (iopt == 35) then
    numOut = 0
    numAtMax = 0
    numIn = 0
    do iow = 1, -numObjInWin
      if (endSeparation(iow) < 0.) then
        numOut = numOut + 1
      elseif (endSeparation(iow) == 1.e10) then
        numAtMax = numAtMax + 1
      else
        numIn = numIn + 1
        ibin = endSeparation(iow) / deltaHist + 1.
        if (ibin <= numBinsHist) &
            graphs(ibin, igraphExtra) = graphs(ibin, igraphExtra) + 1
      endif
    enddo
    write(*,1235) - numObjInWin, numIn, numAtMax, numOut
1235 format(i5,' total ends,',i4,' with close approaches,',i4, &
        ' without,',i4,' not in model')
  else
    indSt = indStart(numWobj) + numPointsObj(numWobj)
    do ind = indSt, indSt + numInWin - 1
      ibin = xModPt(ind) / deltaHist + 1.
      if (ibin <= numBinsHist) &
          graphs(ibin, igraphExtra) = graphs(ibin, igraphExtra) + 1
    enddo
    print *,numInWin, ' total distances between close approaches', &
        ' and ends'
  endif
  go to 40
  !
  ! unshift selected object
  !
240 write(*,'(1x,a,$)') &
      'Object # to unshift: '
  read(in5,*) iobjShift
  call unshift_object(iobjShift, iobjFlag(iobjShift), &
      xModPt, yModPt, zModPt, itype, numWobj, indStart, numPointsObj, xyScale, zScale, &
      zGapStart(indGap), zGapEnd(indGap), numGaps)
  go to 40
241 onlyShifted = .not.onlyShifted
  if (onlyShifted) then
    print *,'Items that failed to shift will be EXCLUDED'
  else
    print *,'Items that failed to shift will be INCLUDED'
  endif
  go to 40

243 if (numInWin <= 0) then
    print *,'No distances have been found in window'
    go to 40
  endif
  sum = 0.
  sumsq = 0.
  print *,'Length of connector lines:'
  do i = 1, numInWin
    ii = indStart(numWobj) + numPointsObj(numWobj) + 3 * (i - 1)
    dist = sqrt((xModPt(ii) - xModPt(ii + 2))**2 + (yModPt(ii) - yModPt(ii + 2))**2 + &
        (zModPt(ii) - zModPt(ii + 2))**2)
    sum = sum + dist
    sumsq = sumsq + dist**2
    write(*,'(f12.5)') dist
  enddo
  call sums_to_avgsd(sum, sumsq, numInWin, connectAvg, connectSd)
  write(*,1243) numInWin, connectAvg, connectSd
1243 format('N = ',i7,'   mean = ',f12.5,'   SD = ',f12.5)
  go to 40

244 nearestOnly = .not. nearestOnly
  if (nearestOnly) then
    print *,'Only distances to nearest neighbor will be recorded; redoing regions'
  else
    print *,'Distances to all neighbors will be recorded; redoing regions'
  endif
  forceLoad = numRegions == 1 .and. (sampled .or. converted .or. shuffled)
  go to 315

end subroutine realGraphicsMain

!
!
!
! general code to redo regions, optionally shuffling points or randomly
! sampling them.  It manages the logicals SHUFFLED and SAMPLED to
! reflect the actual state of points in memory.
! It is expecting:
! IFSAMPLE = 1 to randomly sample points
! IFSHUFFLE =1 to shuffle types
! IFCONVERT = 1 to convert fraction of some types to other types
! FORCELOAD true to force model to get reloaded even if only one region
! GRAPHEACH true to display graphs after each region
!
subroutine redoRegions()
  use mtkredovars
  implicit none

  integer*4 ireg, ii

  do ireg = 1, numRegions
    if (numRegions > 1 .or. forceLoad) then
      if (modelFile .ne. modelRegion(ireg) .or. forceLoad) then
        modelFile = modelRegion(ireg)
        tiltFile = tiltRegion(ireg)
        zScale = zScaleRegion(ireg)
        xyScale = xyScaleRegion(ireg)
        xOffset = xOffsetRegion(ireg)
        yOffset = yOffsetRegion(ireg)
        zOffset = zOffsetRegion(ireg)
        ifFlip = ifFlipRegion(ireg)
        numGaps = numGapsRegion(ireg)
        indGap = indGapRegion(ireg)
        call read_model(modelFile, tiltFile, xyScale, zScale, xOffset, yOffset, &
            zOffset, ifFlip, iobjFlag, LIMFLAG, zGapStart(indGap), &
            zGapEnd(indGap), numGaps)
        onlyShifted = .false.
      endif
      zStart = zStartRegion(ireg)
      zEnd = zEndRegion(ireg)
      call get_objects(zStart, zEnd, xModPt, yModPt, zModPt, indStart, numPointsObj, &
          itype, numWobj, iobjMod, iobjFlag, LIMXYZ, LIMWOBJ)
      shuffled = .false.
      sampled = .false.
      converted = .false.
    endif
    if (ifShuffle .ne. 0) then
      call shuffle(itype, numWobj)
      shuffled = .true.
    endif
    if (ifSample .ne. 0) then
      call random_shifts(xModPt, yModPt, zModPt, indStart, numPointsObj, itype, numWobj, &
          iobjFlag, ranMin, ranMax, probNear, LIMPROBS, deltaNear, &
          numNearBins, numShiftTypes, itypeShift, ishiftFlag, &
          numCheckType, itypeCheck, iuseProb, zGapStart(indGap), &
          zGapEnd(indGap), numGaps, boxTol, iobjBound, &
          ifCheckUnshifted, ranZrel, &
          maxTrials, numTrialCycle, cycleMaxShiftFac, ifScatSurf, xyScale, zScale, &
          manyRandom, ifExcludEout)
      sampled = .true.
    endif
    if (ifConvert .ne. 0) then
      call change_type(itype, numWobj, itypeFrom, itypeTo, changeFrac, numChange)
      converted = .true.
    endif
    !
    if (ifBundEnd == 0) then
      call closedist(xModPt, yModPt, zModPt, indStart, numPointsObj, itype, numWobj, &
          deltaRad, numBins, numGraphs, numRefTypes, numNeighTypes, itypeRef, &
          itypeNeigh, power, limPtsFit, winMin, winMax, numInWin, graphs, &
          areas, iobjWin, numObjInWin, iobjMod, xyzEnd, endSeparation, &
          sampleLen, ifCloseg, ifScatSurf, irefFlag, neighFlag, &
          xyScale, zScale, powerGraph, zGapStart(indGap), &
          zGapEnd(indGap), numGaps, manyRandom, onlyShifted, nearestOnly)
    else
      call bundledist(xModPt, yModPt, zModPt, indStart, numPointsObj, itype, numWobj, &
          deltaRad, numBins, numGraphs, numRefTypes, numNeighTypes, itypeRef, &
          itypeNeigh, iwhichEnd, limPtsFit, padBound, fracOmit, graphs, areas, zScale, &
          zStart, zEnd, xyMaskRegion(1, ireg), winMin, winMax, iobjWin, &
          numObjInWin, xyzEnd)
    endif
    do ii = 1, maxGraph
      ! powergrf(ii) =power
      numBinsByGrf(ii) = numBins
      deltaRadByGrf(ii) = deltaRad
    enddo
    !
    if (numRegions > 1) &
        call addToAvg(graphs, areas, LIMBINS, numGraphs, numBins, ireg)
    if (graphEach) call fourDsp(graphs, LIMBINS, 1, min(4, numGraphs), &
        numBinsByGrf, deltaRadByGrf, xmaxDisplay, ymaxDisplay, igraphDisplay)
  enddo
end subroutine redoRegions

subroutine allocateMtkvars()
  use mtkvars
  integer*4 ierr
  allocate(verts(3, LIMVERTS), indVert(3, LIMTRIANG), triXYrot(3, 2, LIMTRIANG),  &
      triZrot(LIMTRIANG), cosBeta(LIMTRIANG), sinBeta(LIMTRIANG), cosGamma(LIMTRIANG), &
      sinGamma(LIMTRIANG), triangXmin(LIMTRIANG), triangXmax(LIMTRIANG), &
      triangYmin(LIMTRIANG), triangYmax(LIMTRIANG), triangZmin(LIMTRIANG), &
      triangZmax(LIMTRIANG), ibinSave(LIMBINSAVE), aaLine(LIMXYZ), bbLine(LIMXYZ), &
      ccLine(LIMXYZ), ddLine(LIMXYZ), xMin(LIMXYZ), xMax(LIMXYZ), yMin(LIMXYZ), &
      yMax(LIMXYZ), zMin(LIMXYZ), zMax(LIMXYZ), sizes(LIMSIZES), stat = ierr)
  call memoryErrorUC(ierr, 'ARRAYS FOR MTKVARS')
  return
end subroutine allocateMtkvars
