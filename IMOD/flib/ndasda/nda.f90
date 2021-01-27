! NDA (Neighbor density analysis) is a program for analysis of spatial
! point patterns where points may be of different types.  It produces
! graphs of point density versus distance from a point and graphs of
! point density as a function of angular separation between neighbors
! to a reference point.  It can shuffle the types of points and create
! random patterns of points within the boundary containing the real
! points so as to evaluate the statistical significance of apparent
! pattersn in the real points.  See the file NDA.DOC for instructions.
!
! David Mastronarde  7/31/90
!
! $Id$
!
program nda
  call plax_initialize('nda')
  call exit(0)
end program nda

module ndaredovars
  implicit none
  integer LIMBINS, LIMTYPE, LIMVERT
  integer LIMGRAPHS, LIMREGION, LIMPNTS
  parameter (LIMTYPE = 50,  LIMBINS = 1001,  LIMVERT = 5000)
  parameter (LIMGRAPHS = 50, LIMPNTS = 50000, LIMREGION = 200)
  real*4 bx(LIMVERT), by(LIMVERT)            !boundary vertices
  real*4 sx(LIMPNTS), sy(LIMPNTS)            !sample points
  real*4 graphs(LIMBINS,LIMGRAPHS), areas(LIMBINS,LIMGRAPHS)
  integer*4 itype(LIMPNTS)                  !types of sample points
  integer*4 numRefTypes(LIMGRAPHS), numNeighTypes(LIMGRAPHS) !# of types
  integer*4 itypeRef(LIMTYPE,LIMGRAPHS), itypeNeigh(LIMTYPE,LIMGRAPHS)
  integer*4 numAngTypes(LIMGRAPHS), itypeAng(LIMTYPE,LIMGRAPHS)
  integer*4 itypeShift(LIMTYPE)
  integer*4 numInClass(LIMTYPE,LIMREGION)
  integer*4 iobjregion(LIMREGION)
  integer*4 ifScaleRegion(LIMREGION), numTypeBoundRegion(LIMREGION)
  integer*4 itypeBoundRegion(LIMTYPE,LIMREGION)
  integer*4 itypeCrossInd(-256:256), ifConvRegion(LIMREGION)
  integer*4 igraphDisplay(100), itypeFrom(LIMTYPE), itypeTo(LIMTYPE)
  integer*4 numBinsByGrf(LIMGRAPHS), ifAngDiffByGrf(LIMGRAPHS)
  real*4 deltaRadByGrf(LIMGRAPHS), radMinByGrf(LIMGRAPHS), radMaxByGrf(LIMGRAPHS)
  real*4 zzRegion(LIMREGION)
  real*4 xyScaleRegion(LIMREGION)
  real*4 padBoundRegion(LIMREGION), fracomregion(LIMREGION)
  real*4 xmaxDisplay(4), ymaxDisplay(4), changeFrac(LIMTYPE)
  real*4 probNear(LIMBINS)
  integer*4 numRegions, numVert, numTypes
  integer*4 numPoints, numChange
  integer*4 numNearBins, numTypeShift
  integer*4 numBins, numGraphs
  integer*4 maxGraph, iunitAllDist
  integer*4 ifSample, ifShuffle, ifConvert, ifAngDiff
  real*4 deltaNear, boundDistMin, deltaRad
  real*4 radMin, radMax
  logical sampled, converted, forceLoad, graphEach, shuffled
  character*320 modelFile, modelregion(LIMREGION), outFile
  
end module ndaredovars

subroutine realGraphicsMain(comFileArg)
  use ndaredovars
  implicit none
  character*(*) comFileArg
  integer NOPTNEEDMODEL, LIMRAND, ITYPEALL
  parameter (ITYPEALL = 999, LIMRAND = 1000, NOPTNEEDMODEL = 12)
  integer*4 itypeBound(LIMTYPE)
  !
  integer*4 numInTmp(LIMREGION), numPointsRegion(LIMREGION)
  real*4 wholeAreaRegion(LIMREGION)
  !
  logical inside, checkGrf
  !
  real*4 probNearSave(LIMBINS), probNearIn(LIMBINS)
  integer*4 listExtraGraph(LIMGRAPHS)
  character*9 distangtxt
  !
  integer*4 ibaseStart(LIMGRAPHS), ibaseEnd(LIMGRAPHS), integStart(LIMGRAPHS)
  integer*4 integEnd(LIMGRAPHS), numRandAbove(LIMRAND)
  real*4 baseline(LIMGRAPHS), realInteg(LIMGRAPHS), sumInteg(LIMRAND), sumSqInteg(LIMRAND)
  integer*4 nrow, ncol, irow, icol, ifByCol, nxTick, nyTick
  real*4 xgutter, ygutter, ymaxFix
  common /grfarr/ nrow, ncol, irow, icol, ifByCol, nxTick, nyTick, &
      xgutter, ygutter, ymaxFix
  !
  integer*4 ioptNeedModel(NOPTNEEDMODEL) /14, 15, 16, 17, 19, 26, 20, 21, 22, 27, 23, 40/
  logical optionNeedsModel
  character*40 objName
  integer*4 getImodObjName, iobjFromCont
  integer*4 in5
  common /nmsinput/ in5
  !
  real*4 wholeArea, areaSum, avgden, avgint, base, baseVal, c1, c2, c3
  real*4 countSum, deltaNearSave, dens, density, dyTmp, err, fracOmit
  real*4 padBound, padBoundIn, pctAbove, power, randInteg, range, centroid
  real*4 sdDen, sdIntegral, semDen, sum, binSum, binSumSq, to, trueavg, x, xmax
  real*4 xmin, xyScale, yAllMax, yLowTmp, ymax, ymin, zz

  integer*4 i, ib, ib1, ib2, ibinAvgEnd, ibinAvgStart, ibin, ibinEnd, ibinStart
  integer*4 ibinTypeEnd, ibinTypeStart
  integer*4 icont, icontBound, ierr, ifAllSame, ifAnyPlot
  integer*4 ifConvex, ifDoMeanSd, ifDoRand, ifNoPlax, ifNewConv, ifPage, ifRaw, ifReally
  integer*4 ifScale, ifSpecial, ifUseSave, ifxy, ii
  integer*4 line, iobjBound, iobjBoundIn, iopt, iout, iplot, ir, ireg, iregionStart
  integer*4 itcr, ity, iunit, iv, iwin, iwn, jgrbas, jgraph, jgraphAdd, jj, jgrMean
  integer*4 jgrSd, lastAngDiff, numInTot
  integer*4 lineEnd, lineStart, maxTmp, nbinAvg
  integer*4 numNearBinsSave, numDisp, numControlDo, numExtraGraph
  integer*4 numTotControl, numTypeBound

  in5 = 5
  ifAnyPlot = 0
  numExtraGraph = 0
  numTypeBound = 0
  ibinAvgStart = 1
  ibinAvgEnd = 1
  ibinStart = 1
  ibinEnd = 1
  baseVal = 0
  numNearBinsSave = 0
  iwin = 1
  zz = -1.
  ifSpecial = 0
  iunitAllDist = 0

  call openComFile(comFileArg)
  !
  print *,'Enter name of output file to store density values in (Return for none)'
  read(in5, '(a)') modelFile
  if (modelFile == ' ') then
    iout = 6
  else
    iout = 7
    call dopen(7, modelFile, 'new', 'f')
  endif
  iopt = -1
  !
  write(*,'(1x,a,$)') '0 for graphs in a Plax window, 1 to suppress graphs: '
  read(in5,*) ifNoPlax
  call scrnOpen(ifNoPlax)
  !
  modelFile = ' '
  call read_model(modelFile, ifScale, xyScale)
  if (modelFile == ' ') go to 38
  !
  lastAngDiff = -9999
8 write(*,'(1x,a,/,a,$)') 'Enter 0 for graphs of density versus'// &
      ' radial distance', '    or 1 for graphs of density versus'// &
      ' angular difference within an annulus: '
  read(in5,*) ifAngDiff
  !
  ! initialize for first region
  !
10 numRegions = 1
  numTypes = 0
  do i = -256, 256
    itypeCrossInd(i) = 0
  enddo
  !
  ! get boundary of region
  !
12 if (modelFile == ' ') then
    iobjBoundIn = -1
  else
    write(*,'(1x,a,/,a,/,a,$)') 'To specify region boundary, enter' &
        //' IMOD object and contour number,', ' or WIMP object ' &
        //'number AND 0, or 0,0 to take all points on one section,' &
        , '   ( or -1,0 for new model): '
    read(in5,*) iobjBoundIn, icontBound
  endif
  if (iobjBoundIn < 0) then
    modelFile = ' '
    call read_model(modelFile, ifScale, xyScale)
    if (modelFile == ' ') go to 40
    go to 12
  endif
  if (iobjBoundIn == 0) then
    zz = zz + 1.
    write(*,'(1x,a,/,a,/,a,f7.1,i3,a,$)') 'Enter section # (from 0); and 0 ' &
        //'for unpadded rectangle around all points,', '     1 to' &
        //' select special options (subsets of points, padding,' &
        //' convex polygon),', ' or -1 to use previously selected' &
        //' special options (/ for', zz, ifSpecial, '): '
    read(in5,*) zz, ifSpecial
    if (ifSpecial > 0) then
      print *,'Enter list of types to consider in finding boundary (Return for all types)'
      call rdlist(in5, itypeBound, numTypeBound)
      write(*,'(1x,a,$)') 'Distance to pad boundary beyond '// &
          'points (- if in pixels to be scaled): '
      read(in5,*) padBoundIn
      padBound = padBoundIn
      write(*,'(1x,a,$)') '0 for rectangle or 1 to find smallest convex polygon: '
      read(in5,*) ifConvex
      if (ifConvex .ne. 0) then
        write(*,'(1x,a,$)') 'Fraction of far out points to omit from consideration: '
        read(in5,*) fracOmit
      endif
    elseif (ifSpecial < 0) then
      padBound = padBoundIn
    else
      ifConvex = 0
      padBound = 0.
    endif
    if (padBound < 0.) then
      padBound = -padBound
      if (ifScale .ne. 0) padBound = xyScale * padBound
    endif
    if (ifSpecial == 0 .or. numTypeBound == 0) then
      numTypeBound = 1
      itypeBound(1) = ITYPEALL
    endif
    iobjBound = 0
  else
    iobjBound = iobjFromCont(iobjBoundIn, icontBound)
    if (iobjBound == 0) go to 12
  endif
  call get_boundary_obj(iobjBound, bx, by, numVert, zz, itypeBound, &
      numTypeBound, padBound, ifConvex, fracOmit, sx, sy, LIMVERT)
  if (numVert == 0) go to 12
  !
  iobjregion(numRegions) = iobjBound
  zzRegion(numRegions) = zz
  modelregion(numRegions) = modelFile
  ifScaleRegion(numRegions) = ifScale
  xyScaleRegion(numRegions) = xyScale
  numTypeBoundRegion(numRegions) = numTypeBound
  do i = 1, numTypeBound
    itypeBoundRegion(i, numRegions) = itypeBound(i)
  enddo
  padBoundRegion(numRegions) = padBound
  ifconvregion(numRegions) = ifConvex
  fracomregion(numRegions) = fracOmit
  !
  ! get points in region
  !
  call get_points(bx, by, numVert, zz, itypeCrossInd, numTypes, sx, sy, itype, &
      numPoints, numInClass(1, numRegions))
  !
  ! compute area and extent
  !
  wholeArea = 0.
  xmin = bx(1)
  ymin = by(1)
  xmax = bx(1)
  ymax = by(1)
  do iv = 1, numVert
    wholeArea = wholeArea + 0.5 * (by(iv + 1) + by(iv)) * (bx(iv + 1) - bx(iv))
    xmin = min(xmin, bx(iv))
    ymin = min(ymin, by(iv))
    xmax = max(xmax, bx(iv))
    ymax = max(ymax, by(iv))
  enddo
  wholeArea = abs(wholeArea)
  wholeAreaRegion(numRegions) = wholeArea
  numPointsRegion(numRegions) = numPoints
  density = numPoints / wholeArea
  range = max(xmax - xmin, ymax - ymin)
  write(*,101) numPoints, wholeArea, density, range
101 format(i5,' points, area =',f13.6,',  density =',f13.6,/, &
      ' Maximum extent of region is',f12.5)
  !
  write(*,*) 'Type   number    density'
  do ity = -256, 256
    i = itypeCrossInd(ity)
    if (i > 0) then
      objName = ' '
      ierr = getImodObjName(abs(ity), objName)
      dens = numInClass(i, numRegions) / wholeArea
      write(*,'(i5,i8,f13.6,4x,a)') ity, numInClass(i, numRegions), dens, objName
    endif
  enddo
  write(*,*)
  !
  if (numRegions > 1) go to 35
  !
  call getBinSpec(ifAngDiff, deltaRad, numBins, radMin, radMax)
  !
  call getGraphSpec(ifAngDiff, lastAngDiff, numGraphs, itypeRef, numRefTypes, &
      itypeNeigh, numNeighTypes, itypeAng, numAngTypes)
  maxGraph = numGraphs
  !
  !
35 if (ifAngDiff == 0) then
    call dengraph(bx, by, numVert, sx, sy, itype, numPoints, deltaRad, numBins, &
        numGraphs, numRefTypes, numNeighTypes, itypeRef, itypeNeigh, graphs, areas, &
        iunitAllDist)
  else
    call angDist(bx, by, numVert, sx, sy, itype, numPoints, radMin, radMax, numBins, &
        numGraphs, numRefTypes, numNeighTypes, numAngTypes, itypeRef, itypeNeigh, &
        itypeAng, graphs, areas, iunitAllDist)
  endif
  do ii = 1, maxGraph
    numBinsByGrf(ii) = numBins
    deltaRadByGrf(ii) = deltaRad
    ifAngDiffByGrf(ii) = ifAngDiff
    radMinByGrf(ii) = radMin
    radMaxByGrf(ii) = radMax
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
      ' PostScript graphics file',/,' 10/11: Output PostScript', &
      ' graphics file to screen window/printer',/, &
      ' 12: Output single or average graph to file      45: Set PostScript filename',/, &
      ' 13: Loop back to specify model contour defining new area', &
      ' to analyze',/, ' 14: Loop back to specify radial or ', &
      'angular graph and new boundary contour ',/, &
      ' 15: Analyze new region and average with ', &
      'previous region(s)',/, ' 16: Redo current region(s) with', &
      ' new bin size, # of bins, or types for graphs',/, &
      ' 17: Redo current region(s) with angular instead of ', &
      'radial graph or vice versa',/, ' 18: Save bins of a graph' &
      ,' to specify rejection probabilities for random points',/, &
      ' 19/26/20: Do current region(s) with shuffled/converted', &
      ' types or random points',/,' 21: Save current set of ', &
      'points and their types as an IMOD model',/, &
      ' 22/27/23: Do many sets with shuffled/converted', &
      ' types/random pnts for integrals', &
      /,' 24: Take command input from file        25: Exit', &
      /,' 28/29/30 Save a graph/Average/Combine 2 graphs into', &
      ' an extra graph location',/,' 31/32: Save graph in', &
      ' file/Read from file into an extra graph location',/, &
      ' 33: Replace some sets of bins by their averages',/, &
      ' 34/35: Set up special big array for plots/Plot all ', &
      'windows in array',/,' 37/38/39 Add list of graphs/Read', &
      ' list of graphs from file/Read&Add from file',/, &
      ' 42: Export graph values or points for drawing to file',/, &
      ' 43: Start or stop outputting all distances or angles to a file')
  !
40 write(*,'(1x,a,$)') 'Option, or -1 for list of choices: '
  read(in5,*,err = 40) iopt
  if (iopt == -1) go to 38
  if (optionNeedsModel(modelFile, ioptNeedModel, nOptNeedModel, iopt)) go to 40
  go to(201, 202, 203, 204, 205, 206, 207, 208, 209, 210, &
      210, 212, 213, 213, 215, 216, 217, 218, 219, 220, &
      221, 222, 222, 224, 225, 226, 222, 228, 228, 228, &
      228, 228, 228, 234, 235, 40, 228, 228, 228, 40, 40, 228, 243, 40, 245) &
      , iopt
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
  jgraph = igraphDisplay(iwin)
  if (ibinTypeEnd == -1) ibinTypeEnd = numBinsByGrf(jgraph)
  ibinTypeEnd = min(ibinTypeEnd, numBinsByGrf(jgraph))
  ibinTypeStart = max(1, ibinTypeStart)
  lineStart = 1 + (ibinTypeStart - 1) / 5
  lineEnd = 1 + (ibinTypeEnd-1) / 5
  write(*,107)
107 format(12x,'1',14x,'2',14x,'3',14x,'4',14x,'5')
  do line = lineStart, lineEnd
    ib1 = 5 * (line - 1)
    ib2 = min(ibinTypeEnd, ib1 + 5)
    if (ifRaw > 0) then
      write(*,108) ib1, (graphs(ib, jgraph), ib = ib1 + 1, ib2)
108   format(i4,'/',f14.8,4f15.8)
    else
      write(*,1108) ib1, (nint(graphs(ib, jgraph) * areas(ib, jgraph)) &
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
  jgraph = igraphDisplay(iwin)
  ibinAvgEnd = min(ibinAvgEnd, numBinsByGrf(jgraph))
  ibinAvgStart = max(1, ibinAvgStart)
  nbinAvg = ibinAvgEnd + 1 - ibinAvgStart
  call avgsd(graphs(ibinAvgStart, jgraph), nbinAvg, avgden, sdDen, semDen)
  write(*,109) nbinAvg, avgden, sdDen, semDen
109 format(i4,' bins, mean =',f15.8,', SD =',f15.9,', SEM =',f15.9)
  baseVal = avgden
  areaSum = 0.
  countSum = 0.
  do ibin = ibinAvgStart, ibinAvgEnd
    areaSum = areaSum + areas(ibin, jgraph)
    countSum = countSum + areas(ibin, jgraph) * graphs(ibin, jgraph)
  enddo
  trueavg = 0.
  if (areaSum > 0.) trueavg = countSum / areaSum
  write(*,1109) trueavg
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
  jgraph = igraphDisplay(iwin)
  ibinEnd = min(ibinEnd, numBinsByGrf(jgraph))
  ibinStart = max(1, ibinStart)
  call integrate(graphs(1, jgraph), ifAngDiffByGrf(jgraph), deltaRadByGrf(jgraph), &
      radMinByGrf(jgraph), radMaxByGrf(jgraph), numBinsByGrf(jgraph), &
      ibinStart, ibinEnd, 0, 0, baseVal, sum, centroid)
  distangtxt = 'distances'
  if (ifAngDiff .ne. 0) distangtxt = 'angles'
  write(*,110) ibinEnd + 1 - ibinStart, distangtxt, (ibinStart - 1) *  &
      deltaRadByGrf(jgraph), ibinEnd * deltaRadByGrf(jgraph), sum, distangtxt, centroid
110 format(' For the',i4,' bins covering ',a,' from', f10.3,' to',f10.3,/,5x, &
      'the integrated number of (excess/missing) items is',f10.5,/,5x, &
      'the centroid of their ',a,' is', f11.4)
  go to 40
  !
  ! Display graph in window
  !
204 write(*,'(1x,a,$)') 'Display graph in window; Enter graph # and window #: '
  read(in5,*) jgraph, iwin
  if (iwin <= 0 .or. iwin > 4 .or. &
      .not.checkGrf(jgraph, maxGraph, numExtraGraph, listExtraGraph)) go to 40
  xmaxDisplay(iwin) = -1.
  ymaxDisplay(iwin) = -1.
  igraphDisplay(iwin) = jgraph
  call graphDsp(graphs(1, jgraph), numBinsByGrf(jgraph), deltaRadByGrf(jgraph), iwin, &
      jgraph, xmaxDisplay(iwin), ymaxDisplay(iwin))
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
206 write(*,'(1x,a,$)') 'Number of window to rescale; 0 or 1 to rescale X or Y: '
  read(in5,*) iwin, ifxy
  if (iwin <= 0 .or. iwin > 4) go to 40
  jgraph = igraphDisplay(iwin)
  if (jgraph == 0) go to 40
  if (ifxy == 0) then
    radMax = numBinsByGrf(jgraph) * deltaRadByGrf(jgraph)
    write(*,106) ' X', radMax, xmaxDisplay(iwin)
106 format(a,' ranges up to',f13.6,' and current full scale', &
        ' value is',f13.6,/,'   Enter new full scale value: ',$)
    read(in5,*) xmaxDisplay(iwin)
  else
    ymax = 0.
    do i = 1, numBinsByGrf(jgraph)
      ymax = max(ymax, graphs(i, igraphDisplay(iwin)))
    enddo
    write(*,106) ' Y', ymax, ymaxDisplay(iwin)
    read(in5,*) ymaxDisplay(iwin)
  endif
  call graphDsp(graphs(1, jgraph), numBinsByGrf(jgraph), deltaRadByGrf(jgraph), iwin, &
      jgraph, xmaxDisplay(iwin), ymaxDisplay(iwin))
  go to 40
  !
  ! rescale all graphs to same maximum
  !
207 yAllMax = 0
  do iwn = 1, 4
    jgraph = igraphDisplay(iwn)
    if (jgraph > 0) then
      do ibin = 1, numBinsByGrf(jgraph)
        yAllMax = max(yAllMax, graphs(ibin, jgraph))
      enddo
    endif
    call scale(0, 0.9999 * yAllMax, dyTmp, yLowTmp)
    yAllMax = 10. *dyTmp
  enddo
  do iwn = 1, 4
    jgraph = igraphDisplay(iwn)
    if (jgraph > 0) then
      ymaxDisplay(iwn) = yAllMax
      call graphDsp(graphs(1, jgraph), numBinsByGrf(jgraph), deltaRadByGrf(jgraph), &
          iwn, jgraph, xmaxDisplay(iwn), ymaxDisplay(iwn))
    endif
  enddo
  go to 40
  !
  ! Plot one window to metacode file
  !
208 write(*,'(1x,a,$)') 'Window number to plot, plot number or 0 to specify plot: '
  read(in5,*) iwin, iplot
  if (iwin <= 0 .or. iwin > 4) go to 40
  jgraph = igraphDisplay(iwin)
  if (jgraph == 0) go to 40
  if (ifAnyPlot .ne. 0) then
    write(*,2081)
2081 format(' 0 for plot on same page as previous plot(s), 1 for new page: ',$)
    read(in5,*) ifPage
    call psSetup(1, c1, c2, c3, 0)
    if (ifPage .ne. 0) call psFrame
  endif
  call graphPlt(graphs(1, jgraph), numBinsByGrf(jgraph), deltaRadByGrf(jgraph), iplot, &
      jgraph, xmaxDisplay(iwin), ymaxDisplay(iwin))
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
    jgraph = igraphDisplay(iwn)
    if (jgraph > 0) then
      call graphPlt(graphs(1, jgraph), numBinsByGrf(jgraph), deltaRadByGrf(jgraph), &
          iwn, jgraph, xmaxDisplay(iwn), ymaxDisplay(iwn))
    endif
  enddo
  ifAnyPlot = 1
  go to 40
  !
  ! Set filename for graphics output
245 write(*,'(1x,a)') 'Enter name of file for PostScript graphics output'
  read(in5,'(a)')outFile
  call psSetFilename(outFile)
  go to 40
  !
  ! Set up big array of graphs
  !
234 call setGraphArray
  go to 40
  !
  ! Plot all windows to big graph array
  !
235 if (irow <= 0 .or. icol <= 0) go to 40
  do iwn = 1, 4
    jgraph = igraphDisplay(iwn)
    if (jgraph > 0) then
      call psSetup(1, c1, c2, c3, 0)
      if (ifAnyPlot .ne. 0 .and. irow == 1 .and. icol == 1) call psFrame
      if (ymaxFix > 0 .and. ymaxFix .ne. ymaxDisplay(iwn)) then
        ymaxDisplay(iwn) = ymaxFix
        call graphDsp(graphs(1, jgraph), numBinsByGrf(jgraph), deltaRadByGrf(jgraph), &
            iwn, jgraph, xmaxDisplay(iwn), ymaxDisplay(iwn))
      endif
      write(*,2351) jgraph, irow, icol
2351  format(' Plotting graph #',i3,' to row',i3,', column',i3)
      call graphPlt(graphs(1, jgraph), numBinsByGrf(jgraph), deltaRadByGrf(jgraph), &
          - 1, jgraph, xmaxDisplay(iwn), ymaxDisplay(iwn))
      ifAnyPlot = 1
    endif
  enddo
  go to 40
  !
  ! Output plot file to workstation or printer
  !
210 call pltout(11 - iopt)
  ifAnyPlot = 0
  if (irow > 0 .and. icol > 0) then
    irow = 1
    icol = 1
  endif
  go to 40
  !
  ! write graph to file with lots of info
  !
212 write(*,'(1x,a,$)') 'Graph #: '
  read(in5,*) jgraph
  if (.not.checkGrf(jgraph, maxGraph, numExtraGraph, listExtraGraph)) go to 40
  if (jgraph <= numGraphs) then
    iregionStart = numRegions
    jgrbas = jgraph
  else
    iregionStart = 1
    jgrbas = jgraph - numGraphs
  endif
  write(iout, 111) modelFile
111 format(' Model file:  ',a)
  write(iout, 112) ' Boundary object #:', &
      (iobjregion(i), i = iregionStart, numRegions)
112 format(a,(t20,6i10))
  write(iout, 113) ' Z values', (zzRegion(i), i = iregionStart, numRegions)
113 format(a,(t20,6f10.3))
  write(iout, 113) ' Area:', (wholeAreaRegion(i), i = iregionStart, numRegions)
  numInTot = 0
  do ir = iregionStart, numRegions
    numInTmp(ir) = 0
    do ity = 1, numRefTypes(jgrbas)
      if (itypeRef(ity, jgrbas) == ITYPEALL) then
        numInTmp(ir) = numPointsRegion(ir)
        go to 362
      endif
      itcr = itypeCrossInd(itypeRef(ity, jgrbas))
      if (itcr .ne. 0) numInTmp(ir) = numInTmp(ir) + numInClass(itcr, ir)
    enddo
362 numInTot = numInTot + numInTmp(ir)
  enddo
  write(iout, 112) ' Reference items:', (numInTmp(i), i = iregionStart, numRegions)
  if (iregionStart .ne. numRegions) write(iout, 112) '         Total:', numInTot
  !
  numInTot = 0
  do ir = iregionStart, numRegions
    numInTmp(ir) = 0
    do ity = 1, numNeighTypes(jgrbas)
      if (itypeNeigh(ity, jgrbas) == ITYPEALL) then
        numInTmp(ir) = numPointsRegion(ir)
        go to 364
      endif
      itcr = itypeCrossInd(itypeNeigh(ity, jgrbas))
      if (itcr .ne. 0) numInTmp(ir) = numInTmp(ir) + numInClass(itcr, ir)
    enddo
364 numInTot = numInTot + numInTmp(ir)
  enddo
  write(iout, 112) ' Neighboring items:' &
      , (numInTmp(i), i = iregionStart, numRegions)
  if (iregionStart .ne. numRegions) write(iout, 112) '         Total:', numInTot
  !
  if (ifAngDiff > 0) then
    numInTot = 0
    do ir = iregionStart, numRegions
      numInTmp(ir) = 0
      do ity = 1, numAngTypes(jgrbas)
        if (itypeAng(ity, jgrbas) == ITYPEALL) then
          numInTmp(ir) = numPointsRegion(ir)
          go to 366
        endif
        itcr = itypeCrossInd(itypeAng(ity, jgrbas))
        if (itcr .ne. 0) numInTmp(ir) = numInTmp(ir) + numInClass(itcr, ir)
      enddo
366   numInTot = numInTot + numInTmp(ir)
    enddo
    write(iout, 112) ' Annular nbr items:' &
        , (numInTmp(i), i = iregionStart, numRegions)
    if (iregionStart .ne. numRegions) write(iout, 112) '         Total:', numInTot
  endif
  write(iout, 114) ' Reference types:', &
      (itypeRef(ity, jgrbas), ity = 1, numRefTypes(jgrbas))
114 format(a,(12i5))
  write(iout, 114) ' Neighboring types:', &
      (itypeNeigh(ity, jgrbas), ity = 1, numNeighTypes(jgrbas))
  if (ifAngDiff > 0) then
    write(iout, 114) ' Annular nbr types:', &
        (itypeAng(ity, jgrbas), ity = 1, numAngTypes(jgrbas))
    write(iout, 115) radMinByGrf(jgraph), radMaxByGrf(jgraph), deltaRadByGrf(jgraph)
115 format(' Rmin =',f10.3,', rmax =',f10.3, &
        ', delta theta = ',f6.1)
  else
    write(iout, '(a,f12.5)') ' Delta r =', deltaRadByGrf(jgraph)
  endif
  write(iout, '(5g15.8)') (graphs(i, jgraph), i = 1, numBinsByGrf(jgraph))
  go to 40
  !
  ! Toggle output of distances to file
  !
243 write(*,'(1x,a)')  &
        'Enter name of file to save all distances to, or just Enter to stop saving'
  read(in5,'(a)')outFile
  if (iunitAllDist > 0) close(iunitAllDist)
  iunitAllDist = 0
  if (outFile == ' ') go to 40
  call dopen(11, outFile, 'new', 'f')
  iunitAllDist = 11
  go to 40
  !
  ! Specify new boundary object: but check first if averaging
  !
213 if (numRegions > 1) then
    write(*,'(1x,a,$)') 'This will destroy the stored average.'// &
        '  Enter 1 to do so: '
    read(in5,*) ifReally
    if (ifReally .ne. 1) go to 40
  endif
  if (iopt == 13) go to 10
  go to 8
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
217 if (ifAngDiff == 0) then
    ifAngDiff = 1
  else
    ifAngDiff = 0
  endif
  !
  ! Redo region(s) with new bins and/or graphs
  !
216 call getBinSpec(ifAngDiff, deltaRad, numBins, radMin, radMax)
  !
  call getGraphSpec(ifAngDiff, lastAngDiff, numGraphs, itypeRef, numRefTypes, &
      itypeNeigh, numNeighTypes, itypeAng, numAngTypes)
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
  read(in5,*) jgraph, base
  if (.not.checkGrf(jgraph, maxGraph, numExtraGraph, listExtraGraph)) go to 40
  numNearBinsSave = 0
  deltaNearSave = deltaRadByGrf(jgraph)
  do while(numNearBinsSave < numBinsByGrf(jgraph) .and. &
      graphs(numNearBinsSave + 1, jgraph) < base)
    numNearBinsSave = numNearBinsSave + 1
    probNearSave(numNearBinsSave) = graphs(numNearBinsSave, jgraph) / base
  enddo
  print *,numNearBinsSave, ' probability bins saved'
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
    write(*,'(1x,a,$)') 'Number of types to convert: '
    read(in5,*) numChange
    print *,'For each conversion, enter the type to change', &
        ' from, the type to change to,', &
        ' and the fraction of points of that type to convert.'
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
  ! make up control points
  !
220 write(*,'(1x,a,$)') 'Minimum distance from boundary (- if in'// &
      ' pixels to be scaled): '
  read(in5,*) boundDistMin
  if (boundDistMin < 0.) then
    boundDistMin = -boundDistMin
    if (ifScale .ne. 0) boundDistMin = xyScale * boundDistMin
  endif
  !
  print *,'Enter list of types of points to shift ', &
      '(Return for all types)'
  call rdlist(in5, itypeShift, numTypeShift)
  if (numTypeShift == 0) then
    numTypeShift = 1
    itypeShift(1) = ITYPEALL
  endif
  !
  ifUseSave = 0
  if (numNearBinsSave > 0) then
    write(*,'(1x,a,$)') &
        '1 to use saved probability bins, 0 not to: '
    read(in5,*) ifUseSave
  endif
  !
  if (ifUseSave .ne. 0) then
    numNearBins = numNearBinsSave
    deltaNear = deltaNearSave
    do ii = 1, numNearBinsSave
      probNearIn(ii) = probNearSave(ii)
    enddo
  else
    !
    write(*,'(1x,a,$)') '# of bins for rejection, bin size: '
    read(in5,*) numNearBins, deltaNear
    if (numNearBins > 0) then
      write(*,'(1x,a,$)') 'Probability values: '
      read(in5,*) (probNearIn(i), i = 1, numNearBins)
    endif
  endif
  !
  if (numNearBins > 0) then
    ! write(*,'(1x,a,$)') 'Power to raise the probability values to' &
    ! //' for use (try about 3): '
    ! read(in5,*) power
    ! if (power==0.) power=1.
    power = 1.
    do ii = 1, numNearBins
      probNear(ii) = probNearIn(ii)**power
    enddo
  endif
  !
  ifShuffle = 0
  ifSample = 1
  ifConvert = 0
  forceLoad = .false.
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
          ' from, or 0,0 to used fixed value: '
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
    call integrate(graphs(1, jj + jgraphAdd), ifAngDiff, deltaRad, radMin, radMax, &
        numBins, integStart(jj), integEnd(jj), ibaseStart(jj), &
        ibaseEnd(jj), baseline(jj), realInteg(jj), centroid)
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
  graphEach = .false.
  call redoRegions()
  !
  ! Get integrals from graphs and accumulate statistics
  !
  do jj = 1, numGraphs
    call integrate(graphs(1, jj + jgraphAdd), ifAngDiff, deltaRad, radMin, &
        radMax, numBins, integStart(jj), integEnd(jj), ibaseStart(jj), &
        ibaseEnd(jj), baseline(jj), randInteg, centroid)
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
          avgint, sdIntegral)
      pctAbove = 100. *numRandAbove(jj) / numTotControl
      write(iunit, 118) jj + jgraphAdd, realInteg(jj), avgint, sdIntegral, &
          numRandAbove(jj), pctAbove
118   format(i5,f10.5,f12.5,f10.6,i10,f8.2)
    enddo
  enddo
  go to 324
  !
  ! save current points and boundary as model file
  !
221 call save_model(bx, by, numVert, sx, sy, itype, numPoints, zz, ifScale, xyScale)
  call read_model(modelFile, ifScale, xyScale)
  go to 40
  !
  ! Switch to reading in a command file
  !
224 call openComFile(' ')
  go to 40
  !
  ! exit
  !
225 call scrnClose
  call psExit
  !
  ! call to manipulate graphs
  !
228 call manipGraphs(iopt, 'nda', graphs, areas, numBinsByGrf, deltaRadByGrf, &
      ifAngDiffByGrf, radMinByGrf, radMaxByGrf, maxGraph, numExtraGraph, listExtraGraph, &
      igraphDisplay, xmaxDisplay, ymaxDisplay)
  go to 40
end subroutine realGraphicsMain

!CONTAINS

!
! subroutine to redo regions, optionally shuffling points or randomly
! sampling them.  It manages the logicals SHUFFLED and SAMPLED to
! reflect the actual state of points in memory. It is expecting:
! IFSAMPLE = 1 to randomly sample points
! IFSHUFFLE =1 to shuffle types
! IFCONVERT = 1 to convert fraction of some types to other types
! FORCELOAD true to force model to get reloaded even if only one region
! GRAPHEACH true to display graphs after each region
!
subroutine redoRegions()
  use ndaredovars
  implicit none
  integer*4 ireg, ifScale, ii
  real*4 zz, xyScale
  do ireg = 1, numRegions
    if (numRegions > 1 .or. forceLoad) then
      if (modelFile .ne. modelregion(ireg)) then
        call read_model(modelregion(ireg), ifScaleRegion(ireg), &
            xyScaleRegion(ireg))
        modelFile = modelregion(ireg)
        ifScale = ifScaleRegion(ireg)
        xyScale = xyScaleRegion(ireg)
      endif
      zz = zzRegion(ireg)
      call get_boundary_obj(iobjregion(ireg), bx, by, numVert, zz, &
          itypeBoundRegion(1, ireg), numTypeBoundRegion(ireg), &
          padBoundRegion(ireg), ifconvregion(ireg), &
          fracomregion(ireg), sx, sy, LIMVERT)
      call get_points(bx, by, numVert, zz, itypeCrossInd, numTypes, sx, sy, &
          itype, numPoints, numInClass(1, ireg))
      shuffled = .false.
      sampled = .false.
      converted = .false.
    endif
    if (ifShuffle .ne. 0) then
      call shuffle(itype, numPoints)
      shuffled = .true.
    endif
    if (ifSample .ne. 0) then
      call random_points(bx, by, numVert, probNear, deltaNear, numNearBins, &
          boundDistMin, numPoints, itype, itypeShift, numTypeShift, sx, sy)
      sampled = .true.
    endif
    if (ifConvert .ne. 0) then
      call change_type(itype, numPoints, itypeFrom, itypeTo, changeFrac, numChange)
      converted = .true.
    endif
    !
    if (ifAngDiff == 0) then
      call dengraph(bx, by, numVert, sx, sy, itype, numPoints, deltaRad, numBins, &
          numGraphs, numRefTypes, numNeighTypes, itypeRef, itypeNeigh, graphs, areas,  &
          iunitAllDist)
    else
      call angDist(bx, by, numVert, sx, sy, itype, numPoints, radMin, radMax, numBins, &
          numGraphs, numRefTypes, numNeighTypes, numAngTypes, itypeRef, itypeNeigh, &
          itypeAng , graphs, areas, iunitAllDist)
    endif
    do ii = 1, maxGraph
      numBinsByGrf(ii) = numBins
      deltaRadByGrf(ii) = deltaRad
      ifAngDiffByGrf(ii) = ifAngDiff
      radMinByGrf(ii) = radMin
      radMaxByGrf(ii) = radMax
    enddo
    !
    if (numRegions > 1) &
        call addToAvg(graphs, areas, LIMBINS, numGraphs, numBins, ireg)
    if (graphEach) call fourDsp(graphs, LIMBINS, 1, min(4, numGraphs), &
        numBinsByGrf, deltaRadByGrf, xmaxDisplay, ymaxDisplay, igraphDisplay)
  enddo
  return
end subroutine redoRegions
  
