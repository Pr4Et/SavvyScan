! CLOSEDIST produces a series of "graphs" of spatial density
! of neighboring items as a function of radial distance
! from an average reference item.  The reference points and the
! neighboring points may be of various types or sets of types, as
! defined separately for each graph.
! DELTARAD is the bin width (in distance units) desired for the graphs
! NUMBINS is the number of bins for each graph
! NUMGRAPHS is the number of graphs
! NUMREFTYPE is the number of types for reference points for each graph
! NUMNEIGHTYPE is the number of types for neighbor points for each graph
! ITYPEREF(I, J) is the Ith reference type for the Jth graph
! ITYPENEIGH(I, J) is the Ith neighbor type for the Jth graph
! GRAPHS(I, J) is returned with the density at the Ith bin of Jth graph
! FRACSUM(I, J) is returned with total area contributing to that bin
! XMODPT, YMODPT, ZMODPT are arrays with point or line coordinates
! INDSTART is the starting index for each object (contour)
! NUMPNTINOBJ is the number of point in each object
! NUMWOBJ is the number of those objects
! ITYPEOBJ has the type for eavh object
! POWER is the exponent for distance in the volume component
! POWERGRAPH is the array of power values for each graph
! LIMFIT is number of points to fit to lines over (Z based lines)
! WINDOWMIN, WINDOWMAX are window limits for saving connectors
! NUMINWIN  returned with # of distances in window saved in arrays
! IOBJINWIN has the object or negative surface? number of items in window
! NUMOBJINWIN is number of objects in that list (pass in - to find end sep)
! IREFFLAG and NEIGHFLAG have values for the object types (closed, etc)
! XYZEND  has coordinates of endpoints
! ENDSEP is returned with end separation
! SAMPLELEN sampling length for lines
! IFDOCLOSESEG to measure to closest point on segment
! IFUSESCATSURF  flag to measure scattered points from surface
! ZGAPSTART, ZAPGEND are starting and ending Z of gaps
! NUMGAPS is the number of gaps
! XYSCALE, ZSCALE are the pixel size and Z scaling applied to the model
! MANYRANDOM is a flag to suppress output if many random sets being done
! ONLYSHIFTED is logical to compute only agains sucessfully shifted items
! NEARESTONLY is a flag to count nearest neighbors only
!
! $Id$

subroutine closedist(xModPt, yModPt, zModPt, indStart, numPtInObj, itypeObj, numWobj, &
    deltaRad, numBins, numGraphs, numRefType, numNeighType, itypeRef, itypeNeigh, &
    power, limFit, windowMin, windowMax, numInWin, graphs, fracSum, iobjInWin, &
    numObjInWin, iobjImod, xyzEnd, endSep, sampleLen, ifDoCloseSeg, &
    ifUseScatSurf, irefFlag, neighFlag, xyScale, zScale, powerGraph, zGapStart, &
    zGapEnd, numGaps, manyRandom, onlyShifted, nearestOnly)
  use mtkvars
  implicit none
  integer LIMBINS, LIMTEMP, LIMTYPE, itypeAll
  parameter (LIMBINS = 1001, LIMTYPE = 50, itypeAll = 999)
  ! limwobj * (78 + limgraphs) + limxyz * 44 = 3.8 + 11
  ! limverts * 24 =  22  (note the equivalence, and this is shared with mtkrandom)
  parameter (LIMTEMP = 1000)
  real*4 xModPt(*), yModPt(*), zModPt(*)
  integer*4 indStart(*), numPtInObj(*), iobjImod(*)
  real*4 graphs(LIMBINS,*), endSep(*), xyzEnd(3,*)
  integer*4 itypeObj(*)                       !types of sample points
  integer*4 numRefType(*), numNeighType(*)         !# of types for ref and neigh
  integer*4 itypeRef(LIMTYPE,*), itypeNeigh(LIMTYPE,*), iobjInWin(*)
  logical onlyShifted, nearestOnly
  integer*4 igraphRef(LIMGRAPHS), numRefObj(LIMGRAPHS), numClose(2,2), minBin(LIMGRAPHS)
  real*4 xTemp(LIMTEMP), yTemp(LIMTEMP), zTemp(LIMTEMP)
  real*4 xTempMin(LIMTEMP), yTempMin(LIMTEMP), zTempMin(LIMTEMP)
  real*4 xTempMax(LIMTEMP), yTempMax(LIMTEMP), zTempMax(LIMTEMP)
  !
  real*4 fracSum(LIMBINS,*), powerGraph(*)
  real*4 zGapStart(*), zGapEnd(*)
  logical farUp, farDown, findSep, zBasedLines, useBinSave, refDone
  integer*4 numSizeLoaded/0/
  integer*4 indSize(LIM_WIMP_OBJ), nextSize(LIM_MESH_OBJ)
  integer*4 itypeSize(LIM_MESH_OBJ)
  integer*4 iNeigh, ib, ibin, ibinTriang, idir, ierr, ifAny, ifDoCloseSeg, ifLoaded
  integer*4 ifNeighInWin, ifRefInWin, ifSave, ifSkipFail, ifUseScatSurf, iferr
  integer*4 ifreesize, ii, imodObj, inay, ind, indFree, indMesh, indNeed, indNeigh
  integer*4 indRef, indRefInd, indSurfInList, indTri, indTriRot, indexShift
  integer*4 indOfRefEnd, iobjMax, iobjNeigh, iobjNeighEnd, iobjRef, iout, iow, ip
  integer*4 ipoly, iref, irefFlag, is, isNeigh, ishift, isurf, isurfMax, itmp, j1, j3
  integer*4 jj, jnd, jpoly, jshift, jsurf, jtri, kk, limFit, limNeigh, limRef, list
  integer*4 manyRandom, minNeigh, minStart, needFlag, needRef, needed, neighFlag
  integer*4 neighSkip, neighbor, numBinSave, numBins, numWobj, numDoing
  integer*4 numEndObjInWin, numFitNeigh, numFitRef, numGaps, numGraphs, numInWin
  integer*4 numInWinTot, numLoaded, numNeighInWin, numObjInWin, numRefInWin, numTemp

  real*4 angle, angleAvg, angleSD, angleSum, angleSumSq, cosAngle, deltaRad, dist
  real*4 distLim, distMin, distMinSq, distSqr, dx1, dx2, dy1, dy2, dz1, dz2, frac
  real*4 fracEnd, oneMinusT, power, powerUse, radPower, refXmax, refXmin, refYmax
  real*4 refYmin, refZmax, refZmin, sampleFrac, sampleLen, segmentLen
  real*4 sepNeigh, sepRef, shellVol, sizeNeigh, sizeRef, sizeTemp, ttNeigh, ttRef
  real*4 ttn, windowMax, windowMin, x1, x1Min, x2, x2min, xEnd1, xMaxNeigh, xMinNeigh
  real*4 xPrime, xRot, xRot1, xRot2, xStart1, xyScale, y1, y1Min, y2, y2Min, yMaxNeigh
  real*4 yMinNeigh, yRot, yRot1, yRot2, yStart1, yEnd1, z1, z1Min, z2, z2min
  real*4 zEnd1, zEnd2, zMaxNeigh, zMinNeigh, zRot1, zScale, zStart1, zStart2
  integer*4 getImodVerts, getImodSizes
  real*4 acosd

  save numSizeLoaded, indSize, nextSize, itypeSize
  !
  findSep = windowMin < windowMax .and. numObjInWin < 0
  if (findSep) then
    do ii = 1, -numObjInWin
      endSep(ii) = -1.
    enddo
  endif
  numSizeLoaded = 0
  ifreesize = 1
  zBasedLines = irefFlag == 1 .and. neighFlag == 1 .and. &
      sampleLen <= 0. .and. limFit >= 2
  numBinSave = numWobj + (numWobj - 1) * (numWobj - 2) / 2
  useBinSave = irefFlag == 1 .and. neighFlag == 1 .and. sampleLen <= 0. &
      .and. numBinSave < LIMBINSAVE
  !
  ! Analyze which points are valid neighbor points for each graph
  !
  do ii = 1, numWobj
    needed = 0
    indRef = indStart(ii)
    globalXmin(ii) = 1.e10
    globalXmax(ii) = -1.e10
    globalYmin(ii) = 1.e10
    globalYmax(ii) = -1.e10
    globalZmin(ii) = 1.e10
    globalZmax(ii) = -1.e10
    do jj = 1, numGraphs
      isNeighPt(jj, ii) = .false.
      do kk = 1, numNeighType(jj)
        if (itypeNeigh(kk, jj) == itypeAll .or. &
            itypeNeigh(kk, jj) == itypeObj(ii)) isNeighPt(jj, ii) = .true.
      enddo
      if (isNeighPt(jj, ii)) then
        needed = 1
        needFlag = neighFlag
      endif
      do kk = 1, numRefType(jj)
        if (itypeRef(kk, jj) == itypeAll .or. &
            itypeRef(kk, jj) == itypeObj(ii)) then
          needed = 1
          needFlag = irefFlag
        endif
      enddo
    enddo
    if (needed .ne. 0) then
      if (needFlag == 1) then
        !
        ! the old code for line to line fits, assuming progression in Z
        ! precompute all the needed line fits for this MT
        !
        numFitRef = 2
        if (zBasedLines) numFitRef = min(limFit, numPtInObj(ii))
        limRef = indRef + numPtInObj(ii) - numFitRef
        do iref = indRef, limRef
          if (zBasedLines) then
            zStart1 = zModPt(iref)
            zEnd2 = zModPt(iref + numFitRef - 1)
            call linfit(zModPt(iref), xModPt(iref), numFitRef, aaLine(iref), bbLine(iref))
            call linfit(zModPt(iref), yModPt(iref), numFitRef, ccLine(iref), ddLine(iref))
            xStart1 = aaLine(iref) * zStart1 + bbLine(iref)
            yStart1 = ccLine(iref) * zStart1 + ddLine(iref)
            xEnd1 = aaLine(iref) * zEnd2 + bbLine(iref)
            yEnd1 = ccLine(iref) * zEnd2 + ddLine(iref)
          else
            xStart1 = xModPt(iref)
            xEnd1 = xModPt(iref + 1)
            yStart1 = yModPt(iref)
            yEnd1 = yModPt(iref + 1)
            zStart1 = zModPt(iref)
            zEnd2 = zModPt(iref + 1)
          endif
          xMax(iref) = max(xStart1, xEnd1)
          yMax(iref) = max(yStart1, yEnd1)
          zMax(iref) = max(zStart1, zEnd2)
          xMin(iref) = min(xStart1, xEnd1)
          yMin(iref) = min(yStart1, yEnd1)
          zMin(iref) = min(zStart1, zEnd2)
          globalXmin(ii) = min(globalXmin(ii), xMin(iref))
          globalXmax(ii) = max(globalXmax(ii), xMax(iref))
          globalYmin(ii) = min(globalYmin(ii), yMin(iref))
          globalYmax(ii) = max(globalYmax(ii), yMax(iref))
          globalZmin(ii) = min(globalZmin(ii), zMin(iref))
          globalZmax(ii) = max(globalZmax(ii), zMax(iref))
        enddo
      elseif (needFlag == 2) then
        !
        ! scattered points; first load sizes if needed
        !
        limRef = indRef + numPtInObj(ii) - 1
        if (ifUseScatSurf .ne. 0) then
          !
          ! see if sizes are loaded yet for this obj; if not load them
          !
          ifLoaded = 0
          ind = 1
          do while(ind <= numSizeLoaded .and. ifLoaded == 0)
            if (abs(itypeObj(ii)) == itypeSize(ind)) then
              ifLoaded = 1
              indSize(ii) = nextSize(ind)
              nextSize(ind) = nextSize(ind) + numPtInObj(ii)
            endif
            ind = ind + 1
          enddo
          if (ifLoaded == 0) then
            numSizeLoaded = numSizeLoaded + 1
            itypeSize(numSizeLoaded) = abs(itypeObj(ii))
            ierr = getImodSizes(abs(itypeObj(ii)), sizes(ifreesize), &
                LIMSIZES + 1 - ifreesize, numLoaded)
            if (ierr .ne. 0) then
              print *,'Too many point sizes for arrays'
              return
            endif
            indSize(ii) = ifreesize
            nextSize(numSizeLoaded) = ifreesize + numPtInObj(ii)
            ifreesize = ifreesize + numLoaded
          endif
        endif
        !
        ! scale the sizes, set bounding boxes
        !
        do iref = indRef, limRef
          sizeTemp = 0.
          if (ifUseScatSurf .ne. 0) then
            is = indSize(ii) + iref - indRef
            sizes(is) = sizes(is) * xyScale
            sizeTemp = sizes(is)
          endif
          xMin(iref) = xModPt(iref) - sizeTemp
          yMin(iref) = yModPt(iref) - sizeTemp
          zMin(iref) = zModPt(iref) - sizeTemp
          xMax(iref) = xModPt(iref) + sizeTemp
          yMax(iref) = yModPt(iref) + sizeTemp
          zMax(iref) = zModPt(iref) + sizeTemp
          globalXmin(ii) = min(globalXmin(ii), xMin(iref))
          globalXmax(ii) = max(globalXmax(ii), xMax(iref))
          globalYmin(ii) = min(globalYmin(ii), yMin(iref))
          globalYmax(ii) = max(globalYmax(ii), yMax(iref))
          globalZmin(ii) = min(globalZmin(ii), zMin(iref))
          globalZmax(ii) = max(globalZmax(ii), zMax(iref))
        enddo
        !
      else
        print *,'Improper object type; cannot proceed'
        return
      endif
      if (findSep) then
        do jj = 1, -numObjInWin
          if (iobjInWin(jj) == iobjImod(ii)) endSep(jj) = 1.e10
        enddo
      endif
    endif
  enddo
  !
  ! Now see if all objects already have meshes loaded
  !
  if (irefFlag == 4 .or. neighFlag == 4) then
    ifAny = 0
    if (irefFlag == 4) then
      do jj = 1, numGraphs
        do kk = 1, numRefType(jj)
          ifLoaded = 0
          do ind = 1, numMeshLoaded
            if (iobjMesh(ind) == abs(itypeRef(kk, jj))) ifLoaded = 1
          enddo
          if (ifLoaded == 0) ifAny = 1
        enddo
      enddo
    endif
    if (neighFlag == 4) then
      do jj = 1, numGraphs
        do kk = 1, numNeighType(jj)
          ifLoaded = 0
          do ind = 1, numMeshLoaded
            if (iobjMesh(ind) == abs(itypeNeigh(kk, jj))) ifLoaded = 1
          enddo
          if (ifLoaded == 0) ifAny = 1
        enddo
      enddo
    endif
    !
    ! if not, need to load them all over again
    !
    if (ifAny .ne. 0) then
      if (manyRandom == 0) print *,'Processing meshes...'
      numMeshLoaded = 0
      numVerts = 0
      numTriang = 0
      numPoly = 0
      numSurf = 0
      numConts = 0
      do jj = 1, numGraphs
        if (irefFlag == 4) then
          do kk = 1, numRefType(jj)
            ifLoaded = 0
            imodObj = abs(itypeRef(kk, jj))
            do ind = 1, numMeshLoaded
              if (iobjMesh(ind) == imodObj) ifLoaded = 1
            enddo
            if (ifLoaded == 0) then
              call process_mesh(imodObj, xyScale, zScale, ibinSave, &
                  zGapStart, zGapEnd, numGaps, iferr, manyRandom)
              if (iferr .ne. 0) go to 99
            endif
          enddo
        endif
        if (neighFlag == 4) then
          do kk = 1, numNeighType(jj)
            ifLoaded = 0
            imodObj = abs(itypeNeigh(kk, jj))
            do ind = 1, numMeshLoaded
              if (iobjMesh(ind) == imodObj) ifLoaded = 1
            enddo
            if (ifLoaded == 0) then
              call process_mesh(imodObj, xyScale, zScale, ibinSave, &
                  zGapStart, zGapEnd, numGaps, iferr, manyRandom)
              if (iferr .ne. 0) go to 99
            endif
          enddo
        endif
      enddo
    endif
    !
    ! set neighpt for graphs a mesh is needed in
    !
    if (neighFlag == 4) then
      do ii = 1, numMeshLoaded
        do jj = 1, numGraphs
          isNeighPt(jj, ii) = .false.
          do kk = 1, numNeighType(jj)
            if (abs(itypeNeigh(kk, jj)) == iobjMesh(ii)) &
                isNeighPt(jj, ii) = .true.
          enddo
        enddo
      enddo
    endif
    numBinSave = numSurf + (numSurf - 1) * (numSurf - 2) / 2
    useBinSave = irefFlag == 4 .and. neighFlag == 4 &
        .and. numBinSave <= LIMBINSAVE
  endif
  if (manyRandom == 0) print *,'Analyzing closest distances...'
  !
  ! call setfrac(fracsum(1, 1), power, delr, nbins)
  !
  ! zero out the graphs
  !
  do ii = 1, numGraphs
    do jj = 1, numBins
      graphs(jj, ii) = 0.
      fracSum(jj, ii) = 0.
    enddo
    numRefObj(ii) = 0
    minBin(ii) = 0
  enddo
  if (useBinSave) then
    do ind = 1, numBinSave
      ibinSave(ind) = -1
    enddo
  endif
  distLim = deltaRad * numBins
  numInWin = 0
  numInWinTot = 0
  if (.not.findSep) numObjInWin = 0
  indFree = 1
  if (numWobj > 0) indFree = indStart(numWobj) + numPtInObj(numWobj)
  angleSum = 0.
  angleSumSq = 0.
  numRefInWin = 0
  numNeighInWin = 0
  numEndObjInWin = 1 - numObjInWin
  numDoing = 0
  numClose(1, 1) = 0
  numClose(1, 2) = 0
  numClose(2, 1) = 0
  numClose(2, 2) = 0
  indOfRefEnd = numWobj
  if (irefFlag == 4) indOfRefEnd = numMeshLoaded
  !
  ! loop through each sample point considered as reference
  !
  do iobjRef = 1, indOfRefEnd
    !
    ! first make list of graphs the reference point is needed in
    !
    needRef = 0
    do jj = 1, numGraphs
      needed = 0
      if (irefFlag < 4) then
        do kk = 1, numRefType(jj)
          if (itypeRef(kk, jj) == itypeAll .or. &
              itypeRef(kk, jj) == itypeObj(iobjRef)) needed = 1
        enddo
      else
        do kk = 1, numRefType(jj)
          if (abs(itypeRef(kk, jj)) == iobjMesh(iobjRef)) needed = 1
        enddo
      endif
      if (needed > 0) then
        needRef = needRef + 1
        igraphRef(needRef) = jj
        if (.not. nearestOnly) numRefObj(jj) = numRefObj(jj) + 1
      endif
    enddo
    !
    ! get shift index if relevant, cancel needref for failed lines
    !
    if (onlyShifted) then
      ishift = indexShift(iobjRef, irefFlag, itypeObj, numPtInObj, numWobj)
      if (irefFlag == 1 .and. ishift > 0) then
        if (.not.shifted(ishift)) needRef = 0
      endif
    endif
    !
    if (needRef > 0 .and. irefFlag < 4) then
      !
      refDone = .false.
      indRef = indStart(iobjRef)
      limRef = indRef + numPtInObj(iobjRef) - 1
      if (zBasedLines) then
        numFitRef = min(limFit, numPtInObj(iobjRef))
        limRef = indRef + numPtInObj(iobjRef) - numFitRef
      endif
      !
      ! set up if doing samples along the line
      !
      iref = indRef
      sampleFrac = 0.
      if (irefFlag == 1 .and. sampleLen > 0.) call get_next_sample &
          (xModPt, yModPt, zModPt, iref, sampleFrac, limRef, sampleLen, ifDoCloseSeg, &
          xTemp, yTemp, zTemp, numTemp, indRefInd, fracEnd, segmentLen, &
          xTempMin, xTempMax, yTempMin, yTempMax, zTempMin, zTempMax, &
          refXmin, refXmax, refYmin, refYmax, refZmin, refZmax)

      if (irefFlag == 1 .and. manyRandom == 0) then
        numDoing = numDoing + 1
        write(*,'(a,i5,$)') char(13) //'Doing line #', numDoing
        call flush(6)
      endif

      do while(.not.refDone)
        !
        ! if doing lines
        !
        ifSkipFail = 0
        if (irefFlag == 1) then
          !
          ! if doing whole lines: use global min/maxes
          !
          if (sampleLen <= 0.) then
            refXmin = globalXmin(iobjRef)
            refYmin = globalYmin(iobjRef)
            refZmin = globalZmin(iobjRef)
            refXmax = globalXmax(iobjRef)
            refYmax = globalYmax(iobjRef)
            refZmax = globalZmax(iobjRef)
            !
            ! if not z based, copy line and xmin, etc, to tmp
            !
            if (.not.zBasedLines) then
              do iref = indRef, limRef
                iout = iref + 1 - indRef
                xTemp(iout) = xModPt(iref)
                yTemp(iout) = yModPt(iref)
                zTemp(iout) = zModPt(iref)
                xTempMin(iout) = xMin(iref)
                yTempMin(iout) = yMin(iref)
                zTempMin(iout) = zMin(iref)
                xTempMax(iout) = xMax(iref)
                yTempMax(iout) = yMax(iref)
                zTempMax(iout) = zMax(iref)
              enddo
              numTemp = iout
            endif
          else
            !
            ! doing sample segments: should be all set up
            !
          endif
          sizeRef = 0.
        else
          !
          ! doing points, set min's and max's
          !
          refXmin = xMin(iref)
          refYmin = yMin(iref)
          refZmin = zMin(iref)
          refXmax = xMax(iref)
          refYmax = yMax(iref)
          refZmax = zMax(iref)
          sizeRef = 0.
          xTemp(1) = xModPt(iref)
          yTemp(1) = yModPt(iref)
          zTemp(1) = zModPt(iref)
          numTemp = 1
          if (ifUseScatSurf .ne. 0) then
            is = indSize(iobjRef) + iref - indRef
            sizeRef = sizes(is)
          endif
          if (neighFlag .ne. 2 .and. manyRandom == 0) then
            numDoing = numDoing + 1
            write(*,'(a,i7,$)') char(13) //'Doing point #', numDoing
            call flush(6)
          endif
          if (onlyShifted .and. ishift > 0) then
            if (.not.shifted(ishift)) ifSkipFail = 1
          endif
        endif
        !
        if (neighFlag < 4 .and. ifSkipFail == 0) then

          do iobjNeigh = 1, numWobj
            isNeigh = 0
            do indNeed = 1, needRef
              jj = igraphRef(indNeed)
              if (isNeighPt(jj, iobjNeigh)) isNeigh = 1
            enddo
            !
            ! cancel failed line neighbors if appropriate
            !
            if (onlyShifted) then
              jshift = indexShift(iobjNeigh, neighFlag, itypeObj, numPtInObj, numWobj)
              if (neighFlag == 1 .and. jshift > 0) then
                if (.not.shifted(jshift)) isNeigh = 0
              endif
            endif

            if ((iobjNeigh .ne. iobjRef .or. irefFlag == 2) &
                .and. isNeigh .ne. 0) then
              if (useBinSave) then
                iobjMax = max(iobjRef, iobjNeigh)
                ibinTriang = (iobjMax - 1) * (iobjMax - 2) / 2 + min(iobjRef, iobjNeigh)
                ibin = ibinSave(ibinTriang)
              endif

              if (iobjRef < iobjNeigh .or. ibin == -1 .or. &
                  .not.useBinSave) then
                ibin = 0
                distMin = 1.01 * distLim
                distMinSq = distMin**2
                indNeigh = indStart(iobjNeigh)
                limNeigh = indNeigh + numPtInObj(iobjNeigh) - 1
                if (neighFlag == 1) limNeigh = limNeigh - 1
                !
                ! check against global limits
                !
                if (refXmin - globalXmax(iobjNeigh) < distMin .and. &
                    globalXmin(iobjNeigh) - refXmax < distMin .and. &
                    refYmin - globalYmax(iobjNeigh) < distMin .and. &
                    globalYmin(iobjNeigh) - refYmax < distMin .and. &
                    refZmin - globalZmax(iobjNeigh) < distMin .and. &
                    globalZmin(iobjNeigh) - refZmax < distMin) then
                  if (zBasedLines) then
                    !
                    ! OLD CODE FOR Z-BASED LINES
                    !
                    numFitNeigh = min(limFit, numPtInObj(iobjNeigh))
                    limNeigh = indNeigh + numPtInObj(iobjNeigh) - numFitNeigh
                    minNeigh = indNeigh
                    do iref = indRef, limRef
                      zStart1 = zModPt(iref)
                      zEnd2 = zModPt(iref + numFitRef - 1)
                      minStart = minNeigh
                      iNeigh = minStart
                      idir = -1
                      do while(iNeigh <= limNeigh)
                        zStart2 = zModPt(iNeigh)
                        zEnd1 = zModPt(iNeigh + numFitNeigh - 1)
                        farDown = zStart1 - zEnd1 > distMin
                        farUp = zStart2 - zEnd2 > distMin
                        if (.not.(farDown .or. farUp)) then
                          if (xMin(iNeigh) - xMax(iref) < distMin .and. &
                              xMin(iref) - xMax(iNeigh) < distMin) then
                            if (yMin(iNeigh) - yMax(iref) < distMin .and. &
                                yMin(iref) - yMax(iNeigh) < distMin) then
                              call segment_dist_z(aaLine(iref), bbLine(iref), &
                                  ccLine(iref), ddLine(iref), aaLine(iNeigh), &
                                  bbLine(iNeigh), ccLine(iNeigh), ddLine(iNeigh), &
                                  zStart1, zEnd2, zStart2, zEnd1, z1, z2, distSqr)
                              ! print *,dsqr
                              if (distSqr < distMinSq) then
                                distMinSq = distSqr
                                distMin = sqrt(distSqr)
                                minNeigh = iNeigh
                                if (distMin < windowMax .and. &
                                    distMin >= windowMin) then
                                  z1Min = z1
                                  z2min = z2
                                  x1Min = z1 * aaLine(iref) + bbLine(iref)
                                  y1Min = z1 * ccLine(iref) + ddLine(iref)
                                  x2min = z2 * aaLine(iNeigh) + bbLine(iNeigh)
                                  y2Min = z2 * ccLine(iNeigh) + ddLine(iNeigh)
                                  cosAngle = (1. +aaLine(iref) * aaLine(iNeigh) + &
                                      ccLine(iref) * ccLine(iNeigh)) / sqrt( &
                                      (1. +aaLine(iref)**2 + ccLine(iref)**2)* &
                                      (1. +aaLine(iNeigh)**2 + ccLine(iNeigh)**2))
                                endif
                              endif
                            endif
                          endif
                          iNeigh = iNeigh + idir
                        else
                          if (idir == -1 .and. farDown) then
                            iNeigh = indNeigh - 1
                          elseif (farUp) then
                            iNeigh = limNeigh + 1
                          else
                            iNeigh = iNeigh + idir
                          endif
                          if (farDown) minNeigh = max(minNeigh, iNeigh + 1)
                        endif
                        if (iNeigh < indNeigh) then
                          idir = 1
                          iNeigh = minStart + 1
                        endif
                        ! write(*,'(5f10.5)') zs1, zn1, zs2, zn2, distmin
                      enddo
                    enddo
                    !
                    ! NEW LINES TO LINES: seek global minimum
                    !
                  elseif (neighFlag == 1 .and. numTemp > 1) then
                    do itmp = 1, numTemp - 1
                      do iNeigh = indNeigh, limNeigh
                        if (xMin(iNeigh) - xTempMax(itmp) < distMin .and. &
                            xTempMin(itmp) - xMax(iNeigh) < distMin .and. &
                            yMin(iNeigh) - yTempMax(itmp) < distMin .and. &
                            yTempMin(itmp) - yMax(iNeigh) < distMin .and. &
                            zMin(iNeigh) - zTempMax(itmp) < distMin .and. &
                            zTempMin(itmp) - zMax(iNeigh) < distMin) then
                          call segment_dist(xTemp(itmp), yTemp(itmp), &
                              zTemp(itmp), xTemp(itmp + 1), yTemp(itmp + 1), &
                              zTemp(itmp + 1), xModPt(iNeigh), yModPt(iNeigh), &
                              zModPt(iNeigh), xModPt(iNeigh + 1), &
                              yModPt(iNeigh + 1), zModPt(iNeigh + 1), ttRef, &
                              ttNeigh, distSqr)
                          if (distSqr < distMinSq) then
                            distMinSq = distSqr
                            distMin = sqrt(distSqr)
                            if (distMin < windowMax .and. &
                                distMin >= windowMin) then
                              dx1 = xTemp(itmp + 1) - xTemp(itmp)
                              dy1 = yTemp(itmp + 1) - yTemp(itmp)
                              dz1 = zTemp(itmp + 1) - zTemp(itmp)
                              dx2 = xModPt(iNeigh + 1) - xModPt(iNeigh)
                              dy2 = yModPt(iNeigh + 1) - yModPt(iNeigh)
                              dz2 = zModPt(iNeigh + 1) - zModPt(iNeigh)
                              cosAngle = (dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / &
                                  sqrt((dx1**2 + dy1**2 + dz1**2)* &
                                  (dx2**2 + dy2**2 + dz2**2))
                              x1Min = xTemp(itmp) + ttRef * dx1
                              y1Min = yTemp(itmp) + ttRef * dy1
                              z1Min = zTemp(itmp) + ttRef * dz1
                              x2min = xModPt(iNeigh) + ttNeigh * dx2
                              y2Min = yModPt(iNeigh) + ttNeigh * dy2
                              z2min = zModPt(iNeigh) + ttNeigh * dz2
                            endif
                          endif
                        endif
                      enddo
                    enddo
                    !
                    ! SINGLE POINTS TO LINES: again, global minimum
                    !
                  elseif (neighFlag == 1) then
                    do iNeigh = indNeigh, limNeigh
                      if (xMin(iNeigh) - refXmax < distMin .and. &
                          refXmin - xMax(iNeigh) < distMin .and. &
                          yMin(iNeigh) - refYmax < distMin .and. &
                          refYmin - yMax(iNeigh) < distMin .and. &
                          zMin(iNeigh) - refZmax < distMin .and. &
                          refZmin - zMax(iNeigh) < distMin) then
                        inay = iNeigh
                        ip = inay + 1
                        call point_line_dist(xModPt(inay), yModPt(inay), &
                            zModPt(inay), xModPt(ip), yModPt(ip), zModPt(ip), &
                            xTemp(1), yTemp(1), zTemp(1), ttNeigh, distSqr)
                        dist = sqrt(distSqr) - sizeRef
                        if (dist < distMin) then
                          distMin = dist
                          if (distMin < windowMax .and. &
                              distMin >= windowMin) then
                            x1Min = xTemp(1)
                            y1Min = yTemp(1)
                            z1Min = zTemp(1)
                            oneMinusT = 1. -ttNeigh
                            x2min = xModPt(inay) * oneMinusT + xModPt(ip) * ttNeigh
                            y2Min = yModPt(inay) * oneMinusT + yModPt(ip) * ttNeigh
                            z2min = zModPt(inay) * oneMinusT + zModPt(ip) * ttNeigh
                            call trim_win_seg(x1Min, y1Min, z1Min, &
                                x2min, y2Min, z2min, sizeRef, 0.)
                          endif
                        endif
                      endif
                    enddo
                    !
                    ! LINE SEGMENTS TO POINTS: consider each neighbor
                    !
                  elseif (neighFlag == 2 .and. numTemp > 1) then
                    do iNeigh = indNeigh, limNeigh
                      neighSkip = 0
                      if (onlyShifted .and. jshift > 0) then
                        if (.not.shifted(jshift)) neighSkip = 1
                        jshift = jshift + 1
                      endif
                      if (neighSkip == 0) then
                        inay = iNeigh
                        sizeNeigh = 0.
                        if (ifUseScatSurf .ne. 0) then
                          is = indSize(iobjNeigh) + iNeigh - indNeigh
                          sizeNeigh = sizes(is)
                        endif
                        distMin = 1.01 * distLim
                        do itmp = 1, numTemp - 1
                          if (xMin(iNeigh) - xTempMax(itmp) < distMin .and. &
                              xTempMin(itmp) - xMax(iNeigh) < distMin .and. &
                              yMin(iNeigh) - yTempMax(itmp) < distMin .and. &
                              yTempMin(itmp) - yMax(iNeigh) < distMin .and. &
                              zMin(iNeigh) - zTempMax(itmp) < distMin .and. &
                              zTempMin(itmp) - zMax(iNeigh) < distMin) then
                            call point_line_dist(xTemp(itmp), yTemp(itmp), &
                                zTemp(itmp), xTemp(itmp + 1), yTemp(itmp + 1), &
                                zTemp(itmp + 1), xModPt(inay), yModPt(inay), &
                                zModPt(inay), ttRef, distSqr)
                            dist = sqrt(distSqr) - sizeNeigh
                            if (dist < distMin) then
                              distMin = dist
                              oneMinusT = 1. -ttRef
                              x1Min = xTemp(itmp) * oneMinusT + xTemp(itmp + 1) * ttRef
                              y1Min = yTemp(itmp) * oneMinusT + yTemp(itmp + 1) * ttRef
                              z1Min = zTemp(itmp) * oneMinusT + zTemp(itmp + 1) * ttRef
                            endif
                          endif
                        enddo
                        !
                        ! bin the distance to this point
                        !
                        if (distMin < distLim) then
                          ibin = max(1., distMin / deltaRad + 1.)
                          call addDistanceToGraphs(iobjNeigh)
                        endif
                        !
                        ! add to list if in window
                        !
                        if (distMin < windowMax .and. distMin >= windowMin) then
                          x2 = xModPt(inay)
                          y2 = yModPt(inay)
                          z2 = zModPt(inay)
                          call trim_win_seg(x1Min, y1Min, z1Min, x2, &
                              y2, z2, 0., sizeNeigh)
                          call save_connector(x1Min, y1Min, z1Min, x2, &
                              y2, z2, xModPt, yModPt, zModPt, LIMXYZ, indFree, &
                              iobjRef, iobjNeigh, iobjInWin, numObjInWin, &
                              numInWin, numInWinTot, numRefInWin, numNeighInWin)
                        endif
                      endif
                    enddo
                    ! set bin back to zero to avoid getting counted
                    ! below
                    ibin = 0
                    distMin = 1.01 * distLim
                  else
                    !
                    ! POINTS TO POINTS: consider each neighbor
                    !
                    do iNeigh = indNeigh, limNeigh
                      neighSkip = 0
                      if (onlyShifted .and. jshift > 0) then
                        if (.not.shifted(jshift)) neighSkip = 1
                        jshift = jshift + 1
                      endif
                      sizeNeigh = 0.
                      if (ifUseScatSurf .ne. 0) then
                        is = indSize(iobjNeigh) + iNeigh - indNeigh
                        sizeNeigh = sizes(is)
                      endif
                      if (iref .ne. iNeigh .and. neighSkip == 0 .and. &
                          xMin(iNeigh) - refXmax < distLim .and. &
                          refXmin - xMax(iNeigh) < distLim .and. &
                          yMin(iNeigh) - refYmax < distLim .and. &
                          refYmin - yMax(iNeigh) < distLim .and. &
                          zMin(iNeigh) - refZmax < distLim .and. &
                          refZmin - zMax(iNeigh) < distLim) then
                        inay = iNeigh
                        dist = sqrt((xTemp(1) - xModPt(inay))**2 + &
                            (yTemp(1) - yModPt(inay))**2 + &
                            (zTemp(1) - zModPt(inay))**2) - sizeRef - sizeNeigh
                        !
                        ! bin the distance to this point
                        !
                        if (dist < distLim) then
                          ibin = max(1., dist / deltaRad + 1.)
                          call addDistanceToGraphs(iobjNeigh)
                        endif
                        !
                        ! add to list if in window
                        !
                        if (dist < windowMax .and. dist >= windowMin) then
                          x1 = xTemp(1)
                          y1 = yTemp(1)
                          z1 = zTemp(1)
                          x2 = xModPt(inay)
                          y2 = yModPt(inay)
                          z2 = zModPt(inay)
                          call trim_win_seg(x1, y1, z1, x2, y2, z2, &
                              sizeRef, sizeNeigh)
                          ifSave = 1
                          if (irefFlag == 2) then
                            ind = 1
                            do while(ifSave == 1 .and. ind <= numInWin)
                              j1 = indStart(numWobj) + numPtInObj(numWobj) + (ind - 1) * 3
                              j3 = j1 + 2
                              if (x1 == xModPt(j3) .and. x2 == xModPt(j1) .and. &
                                  y1 == yModPt(j3) .and. y2 == yModPt(j1) .and. &
                                  z1 == zModPt(j3) .and. z2 == zModPt(j1)) &
                                  ifSave = 0
                              ind = ind + 1
                            enddo
                          endif
                          if (ifSave .ne. 0) call save_connector(x1, y1, z1, x2, &
                              y2, z2, xModPt, yModPt, zModPt, LIMXYZ, indFree, iobjRef, &
                              iobjNeigh, iobjInWin, numObjInWin, numInWin, numInWinTot, &
                              numRefInWin, numNeighInWin)
                        endif
                      endif
                    enddo
                    ibin = 0
                    distMin = 1.01 * distLim
                  endif
                  !

                  if (distMin < distLim) ibin = max(1., distMin / deltaRad + 1.)
                  if (distMin < windowMax .and. distMin >= windowMin) then
                    angle = acosd(cosAngle)
                    angleSum = angleSum + angle
                    angleSumSq = angleSumSq + angle**2
                    if (findSep) then
                      ifRefInWin = 0
                      ifNeighInWin = 0
                      numInWin = numInWin + 1
                      do iow = 1, -numObjInWin
                        if (iobjImod(iobjRef) == iobjInWin(iow)) then
                          ifRefInWin = 1
                          sepRef = sqrt((x1Min - xyzEnd(1, iow))**2 + &
                              (y1Min - xyzEnd(2, iow))**2 + &
                              (z1Min - xyzEnd(3, iow))**2)
                          endSep(iow) = min(endSep(iow), sepRef)
                          xModPt(indFree) = sepRef
                          indFree = indFree + 1
                        endif
                        if (iobjImod(iobjNeigh) == iobjInWin(iow)) then
                          ifNeighInWin = 1
                          sepNeigh = sqrt((x2min - xyzEnd(1, iow))**2 + &
                              (y2Min - xyzEnd(2, iow))**2 + &
                              (z2min - xyzEnd(3, iow))**2)
                          endSep(iow) = min(endSep(iow), sepNeigh)
                          xModPt(indFree) = sepNeigh
                          indFree = indFree + 1
                        endif
                      enddo
                      numClose(ifRefInWin + 1, ifNeighInWin + 1) = &
                          numClose(ifRefInWin + 1, ifNeighInWin + 1) + 1
                      if (ifRefInWin == 0) then
                        do iow = 1 - numObjInWin, numEndObjInWin
                          if (iobjRef == iobjInWin(iow)) ifRefInWin = 1
                        enddo
                        if (ifRefInWin == 0) then
                          numEndObjInWin = numEndObjInWin + 1
                          iobjInWin(numEndObjInWin) = iobjRef
                          numRefInWin = numRefInWin + 1
                        endif
                      endif
                      if (ifNeighInWin == 0) then
                        do iow = 1 - numObjInWin, numEndObjInWin
                          if (iobjNeigh == iobjInWin(iow)) ifNeighInWin = 1
                        enddo
                        if (ifNeighInWin == 0) then
                          numEndObjInWin = numEndObjInWin + 1
                          iobjInWin(numEndObjInWin) = iobjNeigh
                          numNeighInWin = numNeighInWin + 1
                        endif
                      endif
                    else
                      call save_connector(x1Min, y1Min, z1Min, x2min, &
                          y2Min, z2min, xModPt, yModPt, zModPt, LIMXYZ, indFree, &
                          iobjRef, iobjNeigh, iobjInWin, numObjInWin, numInWin, &
                          numInWinTot, numRefInWin, numNeighInWin)
                    endif
                  endif
                else
                  ! print *,'Eliminated based on global min/max'
                endif
              endif
              if (ibin > 0) then
                call addDistanceToGraphs(iobjNeigh)
              endif
              if (useBinSave) ibinSave(ibinTriang) = ibin
            endif
          enddo
        elseif (ifSkipFail == 0) then
          !
          ! NEIGHBOR IS A MESH
          !
          do indMesh = 1, numMeshLoaded
            isNeigh = 0
            do jj = 1, needRef
              if (isNeighPt(igraphRef(jj), indMesh)) isNeigh = 1
            enddo
            if (isNeigh == 1) then
              !
              ! find a separate distance for each surface
              !
              if (onlyShifted) jshift = indexShift(indMesh, &
                  4, itypeObj, numPtInObj, numWobj)
              do isurf = iobjSurf(indMesh), iobjSurf(indMesh) + numSurfObj(indMesh) - 1
                !
                ! check against surface global limits, eliminate failed
                !
                distMin = 1.01 * distLim
                distMinSq = distMin**2
                neighSkip = 0
                if (onlyShifted .and. jshift > 0) then
                  if (.not.shifted(jshift)) neighSkip = 1
                  jshift = jshift + 1
                endif
                if (neighSkip == 0 .and. &
                    refXmin - surfXmax(isurf) < distMin .and. &
                    surfXmin(isurf) - refXmax < distMin .and. &
                    refYmin - surfYmax(isurf) < distMin .and. &
                    surfYmin(isurf) - refYmax < distMin .and. &
                    refZmin - surfZmax(isurf) < distMin .and. &
                    surfZmin(isurf) - refZmax < distMin) then
                  do list = indStartSurf(isurf), &
                      indStartSurf(isurf) + numInSurf(isurf) - 1
                    ipoly = listSurf(list)
                    !
                    ! check against polygon global limits
                    !
                    if (refXmin - polyXmax(ipoly) < distMin .and. &
                        polyXmin(ipoly) - refXmax < distMin .and. &
                        refYmin - polyYmax(ipoly) < distMin .and. &
                        polyYmin(ipoly) - refYmax < distMin .and. &
                        refZmin - polyZmax(ipoly) < distMin .and. &
                        polyZmin(ipoly) - refZmax < distMin) then
                      !
                      ! scan through triangles
                      !
                      do indTri = indStartPoly(ipoly), indStartPoly(ipoly) + &
                          numInPoly(ipoly) - 1
                        if (refXmin - triangXmax(indTri) < distMin .and. &
                            triangXmin(indTri) - refXmax < distMin .and. &
                            refYmin - triangYmax(indTri) < distMin .and. &
                            triangYmin(indTri) - refYmax < distMin .and. &
                            refZmin - triangZmax(indTri) < distMin .and. &
                            triangZmin(indTri) - refZmax < distMin) then
                          !
                          ! SINGLE POINT TO TRIANGLE
                          !
                          if (numTemp == 1) then
                            call point_to_triangle(xTemp(1), &
                                yTemp(1), zTemp(1), indTri, xRot, yRot, dist)
                            dist = dist - sizeRef
                            if (dist < distMin) then
                              distMin = dist
                              if (distMin < windowMax .and. &
                                  distMin >= windowMin) then
                                x1Min = xTemp(1)
                                y1Min = yTemp(1)
                                z1Min = zTemp(1)
                                xPrime = xRot * cosBeta(indTri) + triZrot(indTri)* &
                                    sinBeta(indTri)
                                x2min = xPrime * cosGamma(indTri) + yRot *  &
                                    sinGamma(indTri)
                                y2Min = -xPrime * sinGamma(indTri) + yRot *  &
                                    cosGamma(indTri)
                                z2min = -xRot * sinBeta(indTri) + triZrot(indTri) * &
                                    cosBeta(indTri)
                                call trim_win_seg(x1Min, y1Min, z1Min, &
                                    x2min, y2Min, z2min, sizeRef, 0.)
                              endif

                            endif
                          else
                            !
                            ! SERIES OF LINE SEGMENTS TO TRIANGLE
                            !
                            do itmp = 1, numTemp - 1
                              if (xTempMin(itmp) - triangXmax(indTri) < &
                                  distMin .and. &
                                  triangXmin(indTri) - xTempMax(itmp) < &
                                  distMin .and. &
                                  yTempMin(itmp) - triangYmax(indTri) < &
                                  distMin .and. &
                                  triangYmin(indTri) - yTempMax(itmp) < &
                                  distMin .and. &
                                  zTempMin(itmp) - triangZmax(indTri) < &
                                  distMin .and. &
                                  triangZmin(indTri) - zTempMax(itmp) < &
                                  distMin) then
                                call segment_to_triangle(xTemp(itmp), &
                                    yTemp(itmp), zTemp(itmp), xTemp(itmp + 1), &
                                    yTemp(itmp + 1), zTemp(itmp + 1), &
                                    indTri, ttRef, xRot, yRot, dist)
                                if (dist < distMin) then
                                  distMin = dist
                                  if (distMin < windowMax .and. &
                                      distMin >= windowMin) then
                                    oneMinusT = 1. -ttRef
                                    x1Min = xTemp(itmp) * oneMinusT + xTemp(itmp + 1) *  &
                                        ttRef
                                    y1Min = yTemp(itmp) * oneMinusT + yTemp(itmp + 1) *  &
                                        ttRef
                                    z1Min = zTemp(itmp) * oneMinusT + zTemp(itmp + 1) *  &
                                        ttRef
                                    xPrime = xRot * cosBeta(indTri) + &
                                        triZrot(indTri) * sinBeta(indTri)
                                    x2min = xPrime * cosGamma(indTri) + &
                                        yRot * sinGamma(indTri)
                                    y2Min = -xPrime * sinGamma(indTri) + &
                                        yRot * cosGamma(indTri)
                                    z2min = -xRot * sinBeta(indTri) + &
                                        triZrot(indTri) * cosBeta(indTri)
                                  endif
                                endif
                              endif
                            enddo
                          endif
                          !
                        endif
                      enddo

                    endif
                  enddo
                endif
                !
                ! add to bins
                !
                if (distMin < distLim) then
                  ibin = max(1., distMin / deltaRad + 1.)
                  call addDistanceToGraphs(indMesh)
!!$                    if (distmin<0.005) then
!!$                    distmin2=0.03
!!$                    distabs=0.015
!!$                    call check_line_surface(isurf, indref, limref, &
!!$                                                  xmt, ymt, zmt, &
!!$                                                  distmin2, distabs, zscal)
!!$                    print *,glbxmin(iobjref), glbxmax(iobjref), &
!!$                                                  glbymin(iobjref), glbymax(iobjref), &
!!$                                                  glbzmin(iobjref), glbzmax(iobjref)
!!$                    print *,surfxmin(isurf), surfxmax(isurf), &
!!$                                                  surfymin(isurf), surfymax(isurf), &
!!$                                                  surfzmin(isurf), surfzmax(isurf)
!!$                    print *,isurf, iobjref, iref, distmin, distmin2
!!$                    endif

                  !
                endif
                !
                ! take care of window stuff here
                !
                if (distMin < windowMax .and. distMin >= windowMin) &
                    call save_connector(x1Min, y1Min, z1Min, x2min, y2Min, &
                    z2min, xModPt, yModPt, zModPt, LIMXYZ, indFree, iobjRef, -isurf, &
                    iobjInWin, numObjInWin, numInWin, numInWinTot, numRefInWin, numNeighInWin)
              enddo
            endif
          enddo
        endif
        !
        ! advance at end of loop and test for termination
        !
        if (irefFlag == 1) then
          !
          ! whole line is done
          !
          refDone = sampleLen <= 0.
          if (.not.refDone) then
            !
            ! segmented line: add cylindrical/spherical shells to
            ! fracsum volume; advance to next sample, done if less than
            ! half of sample length traversed
            !
            do ib = 1, numBins
              shellVol = 3.14159 * deltaRad**2 * (4. *deltaRad * (ib**2 - ib + 0.333333) &
                  + 2 * (ib - 0.5))
              do indNeed = 1, needRef
                jj = igraphRef(indNeed)
                fracSum(ib, jj) = fracSum(ib, jj) + shellVol
              enddo
            enddo
            iref = indRefInd
            sampleFrac = fracEnd
            call get_next_sample &
                (xModPt, yModPt, zModPt, iref, sampleFrac, limRef, sampleLen, &
                ifDoCloseSeg, xTemp, yTemp, zTemp, numTemp, indRefInd, fracEnd, &
                segmentLen, xTempMin, xTempMax, yTempMin, yTempMax, zTempMin, zTempMax, &
                refXmin, refXmax, refYmin, refYmax, refZmin, refZmax)
            refDone = segmentLen < 0.5 * sampleLen
          endif
        else
          !
          ! points: advance iref, add to fracsum
          !
          iref = iref + 1
          refDone = iref > limRef
          if (onlyShifted .and. ishift > 0) ishift = ishift + 1
          do ib = 1, numBins
            shellVol = 1.33333 * 3.14159* &
                ((sizeRef + ib * deltaRad)**3 - (sizeRef + (ib - 1) * deltaRad)**3)
            do indNeed = 1, needRef
              jj = igraphRef(indNeed)
              fracSum(ib, jj) = fracSum(ib, jj) + shellVol
            enddo
          enddo
        endif
        call addNearestToGraphs()
      enddo
    elseif (needRef > 0) then
      !
      ! REFERENCE IS A MESH
      ! loop on surfaces
      !
      do isurf = iobjSurf(iobjRef), iobjSurf(iobjRef) + numSurfObj(iobjRef) - 1
        !
        ! if failed shift, set loop limit to 0
        !
        ifSkipFail = 0
        if (onlyShifted .and. ishift > 0) then
          if (.not.shifted(ishift)) ifSkipFail = 1
          ishift = ishift + 1
        endif
        iobjNeighEnd = 0
        !
        if (ifSkipFail == 0) then
          if (manyRandom == 0) then
            write(*,'(a,i5,$)') char(13) //'Doing surface #', isurf
            call flush(6)
          endif
          refXmin = surfXmin(isurf)
          refYmin = surfYmin(isurf)
          refZmin = surfZmin(isurf)
          refXmax = surfXmax(isurf)
          refYmax = surfYmax(isurf)
          refZmax = surfZmax(isurf)
          iobjNeighEnd = numWobj
          !
          ! add spherical shells to fracsum, with the equivalent radius
          ! taken from the area of the surface
          !
          sizeRef = sqrt(0.25 * surfArea(isurf) / 3.14159)
          do ib = 1, numBins
            shellVol = 1.33333 * 3.14159* &
                ((sizeRef + ib * deltaRad)**3 - (sizeRef + (ib - 1) * deltaRad)**3)
            do indNeed = 1, needRef
              jj = igraphRef(indNeed)
              fracSum(ib, jj) = fracSum(ib, jj) + shellVol
            enddo
          enddo
          !
          if (neighFlag == 4) iobjNeighEnd = numMeshLoaded
        endif
        !
        ! scan through objects for neighbors, unless failed shift
        !
        do iobjNeigh = 1, iobjNeighEnd
          isNeigh = 0
          do indNeed = 1, needRef
            jj = igraphRef(indNeed)
            if (isNeighPt(jj, iobjNeigh)) isNeigh = 1
          enddo
          distMin = 1.01 * distLim
          if (neighFlag < 4 .and. isNeigh .ne. 0) then
            !
            ! POINTS OR LINES
            !
            indNeigh = indStart(iobjNeigh)
            limNeigh = indNeigh + numPtInObj(iobjNeigh) - 1
            if (neighFlag == 1) limNeigh = limNeigh - 1
            neighSkip = 0
            if (onlyShifted) then
              jshift = indexShift(iobjNeigh, neighFlag, itypeObj, &
                  numPtInObj, numWobj)
              if (neighFlag == 1 .and. jshift > 0) then
                if (.not.shifted(jshift)) neighSkip = 1
              endif
            endif
            !
            ! check against global limits of neighbor
            !
            if (neighSkip == 0 .and. &
                refXmin - globalXmax(iobjNeigh) < distMin .and. &
                globalXmin(iobjNeigh) - refXmax < distMin .and. &
                refYmin - globalYmax(iobjNeigh) < distMin .and. &
                globalYmin(iobjNeigh) - refYmax < distMin .and. &
                refZmin - globalZmax(iobjNeigh) < distMin .and. &
                globalZmin(iobjNeigh) - refZmax < distMin) then
              do iNeigh = indNeigh, limNeigh
                if (neighFlag == 2) then
                  distMin = 1.01 * distLim
                  sizeNeigh = 0.
                  if (ifUseScatSurf .ne. 0) then
                    is = indSize(iobjNeigh) + iNeigh - indNeigh
                    sizeNeigh = sizes(is)
                  endif
                  neighSkip = 0
                  if (onlyShifted .and. jshift > 0) then
                    if (.not.shifted(jshift)) neighSkip = 1
                    jshift = jshift + 1
                  endif
                endif
                !
                ! loop on points/segments, check against surface limits
                !
                xMinNeigh = xMin(iNeigh)
                yMinNeigh = yMin(iNeigh)
                zMinNeigh = zMin(iNeigh)
                xMaxNeigh = xMax(iNeigh)
                yMaxNeigh = yMax(iNeigh)
                zMaxNeigh = zMax(iNeigh)
                if (neighSkip == 0 .and. &
                    xMinNeigh - refXmax < distMin .and. &
                    refXmin - xMaxNeigh < distMin .and. &
                    yMinNeigh - refYmax < distMin .and. &
                    refYmin - yMaxNeigh < distMin .and. &
                    zMinNeigh - refZmax < distMin .and. &
                    refZmin - zMaxNeigh < distMin) then
                  !
                  ! loop on polygons in surface
                  !
                  do list = indStartSurf(isurf), &
                      indStartSurf(isurf) + numInSurf(isurf) - 1
                    ipoly = listSurf(list)
                    if (xMinNeigh - polyXmax(ipoly) < distMin .and. &
                        polyXmin(ipoly) - xMaxNeigh < distMin .and. &
                        yMinNeigh - polyYmax(ipoly) < distMin .and. &
                        polyYmin(ipoly) - yMaxNeigh < distMin .and. &
                        zMinNeigh - polyZmax(ipoly) < distMin .and. &
                        polyZmin(ipoly) - zMaxNeigh < distMin) then
                      !
                      ! loop on triangles in polygon
                      !
                      do indTri = indStartPoly(ipoly), &
                          indStartPoly(ipoly) + numInPoly(ipoly) - 1
                        if (xMinNeigh - triangXmax(indTri) < distMin .and. &
                            triangXmin(indTri) - xMaxNeigh < distMin .and. &
                            yMinNeigh - triangYmax(indTri) < distMin .and. &
                            triangYmin(indTri) - yMaxNeigh < distMin .and. &
                            zMinNeigh - triangZmax(indTri) < distMin .and. &
                            triangZmin(indTri) - zMaxNeigh < distMin) then
                          if (neighFlag == 1) then
                            !
                            ! MESH TO LINES
                            !
                            call segment_to_triangle(xModPt(iNeigh), &
                                yModPt(iNeigh), zModPt(iNeigh), xModPt(iNeigh + 1), &
                                yModPt(iNeigh + 1), zModPt(iNeigh + 1), &
                                indTri, ttn, xRot, yRot, dist)
                            if (dist < distMin) then
                              distMin = dist
                              if (distMin < windowMax .and. &
                                  distMin >= windowMin) then
                                oneMinusT = 1. -ttn
                                x1Min = xModPt(iNeigh) * oneMinusT + xModPt(iNeigh + 1) &
                                    * ttn
                                y1Min = yModPt(iNeigh) * oneMinusT + yModPt(iNeigh + 1) &
                                    * ttn
                                z1Min = zModPt(iNeigh) * oneMinusT + zModPt(iNeigh + 1) &
                                    * ttn
                                xPrime = xRot * cosBeta(indTri) + &
                                    triZrot(indTri) * sinBeta(indTri)
                                x2min = xPrime * cosGamma(indTri) + &
                                    yRot * sinGamma(indTri)
                                y2Min = -xPrime * sinGamma(indTri) + &
                                    yRot * cosGamma(indTri)
                                z2min = -xRot * sinBeta(indTri) + &
                                    triZrot(indTri) * cosBeta(indTri)
                              endif
                            endif
                          else
                            !
                            ! MESH TO POINTS
                            !
                            call point_to_triangle(xModPt(iNeigh), &
                                yModPt(iNeigh), zModPt(iNeigh), indTri, &
                                xRot, yRot, dist)
                            dist = dist - sizeNeigh
                            if (dist < distMin) then
                              distMin = dist
                              if (distMin < windowMax .and. &
                                  distMin >= windowMin) then
                                x1Min = xModPt(iNeigh)
                                y1Min = yModPt(iNeigh)
                                z1Min = zModPt(iNeigh)
                                xPrime = xRot * cosBeta(indTri) + triZrot(indTri)* &
                                    sinBeta(indTri)
                                x2min = xPrime * cosGamma(indTri) + yRot * &
                                    sinGamma(indTri)
                                y2Min = -xPrime * sinGamma(indTri) + yRot * &
                                    cosGamma(indTri)
                                z2min = -xRot * sinBeta(indTri) + triZrot(indTri) * &
                                    cosBeta(indTri)
                                call trim_win_seg(x1Min, y1Min, z1Min, &
                                    x2min, y2Min, z2min, sizeNeigh, 0.)
                              endif
                            endif
                          endif
                        endif
                      enddo
                    endif
                  enddo
                endif
                !
                ! doing points: bin distance for each point
                !
                if (neighFlag == 2 .and. distMin < distLim) then
                  ibin = max(1., distMin / deltaRad + 1.)
                  call addDistanceToGraphs(iobjNeigh)
                  !
                  if (distMin < windowMax .and. distMin >= windowMin) &
                      call save_connector(x1Min, y1Min, z1Min, x2min, &
                      y2Min, z2min, xModPt, yModPt, zModPt, LIMXYZ, indFree, -isurf, &
                      iobjNeigh, iobjInWin, numObjInWin, numInWin, numInWinTot, &
                      numRefInWin, numNeighInWin)
                endif

              enddo
              !
              ! doing lines: bin distance after getting global minimum
              !
              if (neighFlag == 1 .and. distMin < distLim) then
                ibin = max(1., distMin / deltaRad + 1.)
                call addDistanceToGraphs(iobjNeigh)
                !
                if (distMin < windowMax .and. distMin >= windowMin) &
                    call save_connector(x1Min, y1Min, z1Min, x2min, &
                    y2Min, z2min, xModPt, yModPt, zModPt, LIMXYZ, indFree, -isurf, &
                    iobjNeigh, iobjInWin, numObjInWin, numInWin, numInWinTot, &
                    numRefInWin, numNeighInWin)
              endif
            endif
          elseif (isNeigh .ne. 0) then
            !
            ! MESH TO MESH
            ! loop through surfaces of neighbor, skipping failed shifts
            !
            if (onlyShifted) jshift = indexShift(iobjNeigh, &
                4, itypeObj, numPtInObj, numWobj)

            do jsurf = iobjSurf(iobjNeigh), iobjSurf(iobjNeigh) + &
                numSurfObj(iobjNeigh) - 1
              neighSkip = 0
              if (onlyShifted .and. jshift > 0) then
                if (.not.shifted(jshift)) neighSkip = 1
                jshift = jshift + 1
              endif
              if (isurf .ne. jsurf .and. neighSkip == 0) then
                if (useBinSave) then
                  isurfMax = max(isurf, jsurf)
                  ibinTriang = (isurfMax - 1) * (isurfMax - 2) / 2 + min(isurf, jsurf)
                  ibin = ibinSave(ibinTriang)
                endif

                if (isurf < jsurf .or. ibin == -1 .or. &
                    .not.useBinSave) then

                  distMin = 1.01 * distLim
                  ibin = 0
                  if (refXmin - surfXmax(jsurf) < distMin .and. &
                      surfXmin(jsurf) - refXmax < distMin .and. &
                      refYmin - surfYmax(jsurf) < distMin .and. &
                      surfYmin(jsurf) - refYmax < distMin .and. &
                      refZmin - surfZmax(jsurf) < distMin .and. &
                      surfZmin(jsurf) - refZmax < distMin) then
                    !
                    ! loop on polygons of neighbor surface
                    !
                    do indSurfInList = indStartSurf(jsurf), &
                        indStartSurf(jsurf) + numInSurf(jsurf) - 1
                      jpoly = listSurf(indSurfInList)
                      if (refXmin - polyXmax(jpoly) < distMin .and. &
                          polyXmin(jpoly) - refXmax < distMin .and. &
                          refYmin - polyYmax(jpoly) < distMin .and. &
                          polyYmin(jpoly) - refYmax < distMin .and. &
                          refZmin - polyZmax(jpoly) < distMin .and. &
                          polyZmin(jpoly) - refZmax < distMin) then
                        !
                        ! loop on triangles of neighbor
                        !
                        do jtri = indStartPoly(jpoly), &
                            indStartPoly(jpoly) + numInPoly(jpoly) - 1
                          xMinNeigh = triangXmin(jtri)
                          yMinNeigh = triangYmin(jtri)
                          zMinNeigh = triangZmin(jtri)
                          xMaxNeigh = triangXmax(jtri)
                          yMaxNeigh = triangYmax(jtri)
                          zMaxNeigh = triangZmax(jtri)
                          if (refXmin - xMaxNeigh < distMin .and. &
                              xMinNeigh - refXmax < distMin .and. &
                              refYmin - yMaxNeigh < distMin .and. &
                              yMinNeigh - refYmax < distMin .and. &
                              refZmin - zMaxNeigh < distMin .and. &
                              zMinNeigh - refZmax < distMin) then
                            !
                            ! now loop on polygons of reference surface
                            !
                            do list = indStartSurf(isurf), &
                                indStartSurf(isurf) + numInSurf(isurf) - 1
                              ipoly = listSurf(list)
                              if (xMinNeigh - polyXmax(ipoly) < distMin .and. &
                                  polyXmin(ipoly) - xMaxNeigh < distMin .and. &
                                  yMinNeigh - polyYmax(ipoly) < distMin .and. &
                                  polyYmin(ipoly) - yMaxNeigh < distMin .and. &
                                  zMinNeigh - polyZmax(ipoly) < distMin .and. &
                                  polyZmin(ipoly) - zMaxNeigh < distMin) then
                                !
                                ! loop on triangles in reference polygon
                                !
                                do indTri = indStartPoly(ipoly), &
                                    indStartPoly(ipoly) + numInPoly(ipoly) - 1
                                  if (xMinNeigh - triangXmax(indTri) < distMin .and. &
                                      triangXmin(indTri) - xMaxNeigh < distMin .and. &
                                      yMinNeigh - triangYmax(indTri) < distMin .and. &
                                      triangYmin(indTri) - yMaxNeigh < distMin .and. &
                                      zMinNeigh - triangZmax(indTri) < distMin .and. &
                                      triangZmin(indTri) - zMaxNeigh < distMin) then
                                    call triangle_to_triangle(indTri, jtri, xRot1, &
                                        yRot1, zRot1, xRot2, yRot2, indTriRot, dist)
                                    if (dist < distMin) then
                                      distMin = dist
                                      if (distMin < windowMax .and. &
                                          distMin >= windowMin) then
                                        xPrime = xRot1 * cosBeta(indTriRot) + &
                                            zRot1 * sinBeta(indTriRot)
                                        x1Min = xPrime * cosGamma(indTriRot) + &
                                            yRot1 * sinGamma(indTriRot)
                                        y1Min = -xPrime * sinGamma(indTriRot) + yRot1* &
                                            cosGamma(indTriRot)
                                        z1Min = -xRot1 * sinBeta(indTriRot) + zRot1* &
                                            cosBeta(indTriRot)
                                        xPrime = xRot2 * cosBeta(indTriRot) + &
                                            triZrot(indTriRot) * sinBeta(indTriRot)
                                        x2min = xPrime * cosGamma(indTriRot) + &
                                            yRot2 * sinGamma(indTriRot)
                                        y2Min = -xPrime * sinGamma(indTriRot) + yRot2* &
                                            cosGamma(indTriRot)
                                        z2min = -xRot2 * sinBeta(indTriRot) + &
                                            triZrot(indTriRot) * cosBeta(indTriRot)
                                      endif
                                    endif
                                  endif
                                enddo
                              endif
                            enddo
                          endif
                        enddo
                      endif
                    enddo
                    !
                    ! doing meshes: bin distance after getting global
                    ! minimum
                    !
                    if (distMin < distLim) then
                      ibin = max(1., distMin / deltaRad + 1.)
                      !
                      if (distMin < windowMax .and. distMin >= windowMin) &
                          call save_connector(x1Min, y1Min, z1Min, x2min, &
                          y2Min, z2min, xModPt, yModPt, zModPt, LIMXYZ, indFree, -isurf, &
                          - jsurf, iobjInWin, numObjInWin, numInWin, numInWinTot, &
                          numRefInWin, numNeighInWin)
                    endif

                  endif
!!$                    if (distmin<0.005) then
!!$                    distmin2=0.03
!!$                    distabs=0.015
!!$                    call check_two_meshes(isurf, jsurf, distmin2, &
!!$                                                  distabs, zscal)
!!$                    print *,'close', isurf, jsurf, distmin, distmin2
!!$                    endif

                endif
                if (ibin > 0) then
                  call addDistanceToGraphs(iobjNeigh)
                endif
                if (useBinSave) ibinSave(ibinTriang) = ibin
              endif
            enddo
          endif
        enddo
        call addNearestToGraphs()
      enddo
    endif
  enddo
  !
  ! scale counts: use power parameter for whole lines to compute fracsum,
  ! otherwise set power to 2 (fracsum should be complete)
  !
  if ((irefFlag == 1 .and. sampleLen <= 0.) .or. nearestOnly) then
    do ibin = 1, numBins
      if (power == 0.) then
        radPower = 1.
      else
        radPower = (2. *(ibin - 0.5) * deltaRad)**power
      endif
      frac = radPower * deltaRad * 3.14159
      if (nearestOnly) frac = 1.
      do jj = 1, numGraphs
        fracSum(ibin, jj) = frac * numRefObj(jj)
        graphs(ibin, jj) = graphs(ibin, jj) / fracSum(ibin, jj)
      enddo
    enddo
    powerUse = power
    if (nearestOnly) powerUse = 0
  else
    do ibin = 1, numBins
      do jj = 1, numGraphs
        graphs(ibin, jj) = graphs(ibin, jj) / fracSum(ibin, jj)
      enddo
    enddo
    powerUse = 2
  endif
  do jj = 1, numGraphs
    powerGraph(jj) = powerUse
  enddo
  !
  if (numInWinTot > 0.) then
    call sums_to_avgsd(angleSum, angleSumSq, numInWinTot, angleAvg, angleSD)
    write(*,106) numInWinTot, angleAvg, angleSD, numRefInWin, numNeighInWin
106 format(i4,' distances in window; angle mean=',f6.2,',  SD=', &
        f6.2,/,i6,' reference and',i5,' neighbor objects')
    if (numInWin < numInWinTot) print *,numInWinTot - numInWin, &
        ' connectors could not be saved - arrays full'
  endif
  if (findSep) then
    write(*,108) ((numClose(ind, jnd), ind = 1, 2), jnd = 1, 2)
108 format(i5,' &',i4,' distances in window for neighbor OFF end', &
        ' list, reference OFF & ON',/,i5,' &',i4, &
        ' distances in window for neighbor ON end', &
        ' list, reference OFF & ON')
    numInWin = indFree - (indStart(numWobj) + numPtInObj(numWobj))
  endif
  return
99 print *,'Data not loaded; try fewer objects'
  numMeshLoaded = 0
  return

CONTAINS

  ! Add the bin to each of the graphs that includes it, unless doing nearest neighbor
  ! then just maintain the minimum for the graph
  subroutine addDistanceToGraphs(neighbor)
    integer*4 neighbor
    do indNeed = 1, needRef
      jj = igraphRef(indNeed)
      if (isNeighPt(jj, neighbor)) then
        if (nearestOnly) then
          if (minBin(jj) == 0 .or. ibin < minBin(jj)) minBin(jj) = ibin
        else
          graphs(ibin, jj) = graphs(ibin, jj) + 1.
        endif
      endif
    enddo
    return
  end subroutine addDistanceToGraphs

  ! After getting a minimu distance to a reference, add it to the graphs
  subroutine addNearestToGraphs()
    if (.not. nearestOnly) return
    do indNeed = 1, needRef
      jj = igraphRef(indNeed)
      if (minBin(jj) > 0) graphs(minBin(jj), jj) = graphs(minBin(jj), jj) + 1.
      numRefObj(jj) = numRefObj(jj) + 1
      minBin(jj) = 0
    enddo
    return
  end subroutine addNearestToGraphs

end subroutine closedist
