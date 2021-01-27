! $Id$
!
subroutine random_shifts(xModPt, yModPt, zModPt, indStart, numPtInObj, itypeObj, &
    numWobj, iobjFlag, ranMin, ranMax, probNear, limProbs, deltaNear, numNearBins, &
    numShiftTypes, itypeShift, ishiftFlag, numCheckTypes, itypeCheck, indProbCurve, &
    zGapStart, zGapEnd, numGaps, boxToler, iobjBound, ifCheckUnshifted, ranZrel, &
    maxTrials, numTrialCycle, cycleFactor, ifUseScatSurf, xyScale, zScale, manyRandom, &
    ifExcludeOut)
  use mtkvars
  implicit none
  integer LIMBOUND, LIMTYPE, itypeAll
  parameter (LIMTYPE = 50, itypeAll = 999, LIMBOUND = 1000)
  ! limwobj * 28 = 0.8
  ! limxyz * 40 = 10.0
  ! limverts * 24 = 22
  real*4 xModPt(*), yModPt(*), zModPt(*), probNear(limProbs,*), deltaNear(*)
  integer*4 indStart(*), numPtInObj(*), iobjFlag(*)
  integer*4 itypeObj(*)                       !types of sample points
  integer*4 itypeShift(*), itypeCheck(*), numNearBins(*), indProbCurve(*)
  real*4 zGapStart(*), zGapEnd(*)
  real * 4 delXYZ(3)
  real*4 probUse(250)
  integer*4 listBound(LIMBOUND)
  real*4 zBound(LIMBOUND)
  integer*4 numSizeLoaded/0/
  integer*4 indSize(LIM_WIMP_OBJ), nextSize(LIM_MESH_OBJ)
  integer*4 itypeSize(LIM_MESH_OBJ)
  integer*4 ierr, ifCheckFixed, ifCheckUnshifted, ifErr, ifExcludeOut, ifFirst
  integer*4 ifLoaded, ifOnList, ifUseScatSurf, ifZeroedShift, ifreeSize, ii, imo
  integer*4 imodObj, inListToShift, inay, ind, indCheckShift, indMesh, indNeigh
  integer*4 indRef, ineigh, iobjBound, iobjNeigh, iobjRef, ip, iref, iroundNum, is
  integer*4 isBase, iseed, ishift, ishiftFlag, isurf, itypeCheckTemp, jmesh, jprob
  integer*4 jsurf, kk, limNeigh, limProbs, limRef, locShiftInd, manyRandom, maxTrials
  integer*4 meshCheck, needFlag, needTotal, needed, numBound, numCheckTemp
  integer*4 numCheckTypes, numContTrim, numExclude, numFail, numFinal, numGaps
  integer*4 numItemSum, numItemsTemp, numItemsTot, numLoaded, numMeshTrim
  integer*4 numPolyTrim, numShiftTypes, numSurfTrim, numTrialCycle, numTrialTot
  integer*4 numTriangTrim, numVertsTrim, numWobj, ibaseShift, icall, icheck

  real*4 ranMax, ranMin, ranZrel, refXmax, refXmin, refYmax, refYmin, refZmax
  real*4 refZmin, shiftXmax
  real*4 shiftXmin, shiftYmax, shiftYmin, shiftZmax, shiftZmin, sizeNeigh, sizeRef
  real*4 sizeTemp, ttNeigh, ttRef, xEnd1, xStart1, xyScale, yStart1, zEnd1, zScale
  real*4 avgNumTrial, boxToler, cycleFactor, delX, delY, delZ, deltaUse, dist
  real*4 distAbs, distLimit, distMin, dminSq, dsqr, zStart1, yEnd1
  save numSizeLoaded, indSize, nextSize, itypeSize
  equivalence (delX, delXYZ(1)), (delY, delXYZ(2)), (delZ, delXYZ(3))
  logical doCheck, badShift, passing, outside_boundary
  real*4 b3dran
  integer*4 getImodSizes
  data ifFirst/1/
  save ifFirst, iseed
  data icall/0/
  save icall
  !
  character*8 seedTime
  !
  numSizeLoaded = 0
  ifreeSize = 1
  if (ifFirst .ne. 0) then
    call time(seedTime)
    iseed = 2 * ((ichar(seedTime(8:8)) + 128 * (ichar(seedTime(7:7)) + &
        128 * (ichar(seedTime(5:5)) + 128 * ichar(seedTime(4:4))))) / 2) + 1
    ifFirst = 0
  endif
  ! print *,' '
  ! print *,'seed =', iseed
  ! iseed = 112007479
  ! if (icall==1) iseed=112023601
  ! icall=icall+1
  shiftXmin = 1.e10
  shiftXmax = -1.e10
  shiftYmin = 1.e10
  shiftYmax = -1.e10
  shiftZmin = 1.e10
  shiftZmax = -1.e10
  call get_random_boundary(iobjBound, listBound, zBound, numBound, LIMBOUND)
  !
  ! scan through the objects to be shifted, setting up the table of
  ! shifts and undoing shifts that already exist, for lines and points
  !
  ibaseShift = 1
  numItemsTot = 0
  if (numObjShifted > 0) ibaseShift = indStartShift(numObjShifted) + &
      numItemShifted(numObjShifted)
  do kk = 1, numShiftTypes
    if (iobjFlag(abs(itypeShift(kk))) == 1) then
      ifOnList = 0
      do imo = 1, numObjShifted
        if (iobjShift(imo) == itypeShift(kk)) ifOnList = imo
      enddo
      if (ifOnList .ne. 0) then
        !
        ! if this line object is already on the list, shift it back and
        ! zero shifts
        !
        ishift = indStartShift(ifOnList)
        do ii = 1, numWobj
          if (itypeShift(kk) == itypeObj(ii)) then
            do ip = indStart(ii), indStart(ii) + numPtInObj(ii) - 1
              xModPt(ip) = xModPt(ip) - shifts(1, ishift)
              yModPt(ip) = yModPt(ip) - shifts(2, ishift)
              zModPt(ip) = zModPt(ip) - shifts(3, ishift)
            enddo
            shifts(1, ishift) = 0.
            shifts(2, ishift) = 0.
            shifts(3, ishift) = 0.
            shifted(ishift) = .true.
            ishift = ishift + 1
          endif
        enddo
      else
        !
        ! otherwise, set up the right number of shifts for the lines
        !
        numObjShifted = numObjShifted + 1
        indStartShift(numObjShifted) = ibaseShift
        iobjShift(numObjShifted) = itypeShift(kk)
        numItemsTemp = 0
        do ii = 1, numWobj
          if (itypeShift(kk) == itypeObj(ii)) then
            numItemsTemp = numItemsTemp + 1
            shifts(1, ibaseShift) = 0.
            shifts(2, ibaseShift) = 0.
            shifts(3, ibaseShift) = 0.
            shifted(ibaseShift) = .true.
            ibaseShift = ibaseShift + 1
            if (ibaseShift > LIMSHIFT) then
              print *,'Too many shifts for arrays'
              return
            endif
          endif
        enddo
        numItemShifted(numObjShifted) = numItemsTemp
        ifOnList = numObjShifted
      endif
      numItemsTot = numItemsTot + numItemShifted(ifOnList)
    elseif (iobjFlag(abs(itypeShift(kk))) == 2) then
      ifOnList = 0
      do imo = 1, numObjShifted
        if (iobjShift(imo) == itypeShift(kk)) ifOnList = imo
      enddo
      if (ifOnList > 0) then
        !
        ! if points are already on the list, shift them back and zero
        ! shifts
        !
        ishift = indStartShift(ifOnList)
        do ii = 1, numWobj
          if (itypeShift(kk) == itypeObj(ii)) then
            do ip = indStart(ii), indStart(ii) + numPtInObj(ii) - 1
              xModPt(ip) = xModPt(ip) - shifts(1, ishift)
              yModPt(ip) = yModPt(ip) - shifts(2, ishift)
              zModPt(ip) = zModPt(ip) - shifts(3, ishift)
              shifts(1, ishift) = 0.
              shifts(2, ishift) = 0.
              shifts(3, ishift) = 0.
              shifted(ishift) = .true.
              ishift = ishift + 1
            enddo
          endif
        enddo
      else
        !
        ! otherwise, set up the right number of shifts for the points
        !
        numObjShifted = numObjShifted + 1
        indStartShift(numObjShifted) = ibaseShift
        iobjShift(numObjShifted) = itypeShift(kk)
        numItemsTemp = 0
        do ii = 1, numWobj
          if (itypeShift(kk) == itypeObj(ii)) then
            do ip = 1, numPtInObj(ii)
              numItemsTemp = numItemsTemp + 1
              shifts(1, ibaseShift) = 0.
              shifts(2, ibaseShift) = 0.
              shifts(3, ibaseShift) = 0.
              shifted(ibaseShift) = .true.
              ibaseShift = ibaseShift + 1
              if (ibaseShift > LIMSHIFT) then
                print *,'Too many shifts for arrays'
                return
              endif
            enddo
          endif
        enddo
        numItemShifted(numObjShifted) = numItemsTemp
        ifOnList = numObjShifted
      endif
      numItemsTot = numItemsTot + numItemShifted(ifOnList)
    endif
  enddo
  !
  ! scan through loaded "mt" objects, load all sizes if needed
  !
  do ii = 1, numWobj
    needed = 0
    indRef = indStart(ii)
    limRef = indRef + numPtInObj(ii) - 1
    globalXmin(ii) = 1.e10
    globalXmax(ii) = -1.e10
    globalYmin(ii) = 1.e10
    globalYmax(ii) = -1.e10
    globalZmin(ii) = 1.e10
    globalZmax(ii) = -1.e10
    do kk = 1, numShiftTypes
      if (itypeShift(kk) == itypeObj(ii)) needed = 1
    enddo
    inListToShift = needed
    do kk = 1, numCheckTypes
      if (itypeCheck(kk) == itypeObj(ii)) needed = 1
    enddo
    if (needed .ne. 0) then
      needFlag = iobjFlag(abs(itypeObj(ii)))
      if (needFlag == 1) then
        do iref = indRef, limRef - 1
          xStart1 = xModPt(iref)
          xEnd1 = xModPt(iref + 1)
          yStart1 = yModPt(iref)
          yEnd1 = yModPt(iref + 1)
          zStart1 = zModPt(iref)
          zEnd1 = zModPt(iref + 1)
          xMax(iref) = max(xStart1, xEnd1)
          yMax(iref) = max(yStart1, yEnd1)
          zMax(iref) = max(zStart1, zEnd1)
          xMin(iref) = min(xStart1, xEnd1)
          yMin(iref) = min(yStart1, yEnd1)
          zMin(iref) = min(zStart1, zEnd1)
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
            ierr = getImodSizes(abs(itypeObj(ii)), sizes(ifreeSize), &
                LIMSIZES + 1 - ifreeSize, numLoaded)
            if (ierr .ne. 0) then
              print *,'Too many point sizes for arrays'
              return
            endif
            indSize(ii) = ifreeSize
            nextSize(numSizeLoaded) = ifreeSize + numPtInObj(ii)
            ifreeSize = ifreeSize + numLoaded
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
        print *,'Inappropriate object type in MT lists'
        return
      endif
      if (inListToShift == 1) then
        shiftXmin = min(shiftXmin, globalXmin(ii))
        shiftYmin = min(shiftYmin, globalYmin(ii))
        shiftZmin = min(shiftZmin, globalZmin(ii))
        shiftXmax = max(shiftXmax, globalXmax(ii))
        shiftYmax = max(shiftYmax, globalYmax(ii))
        shiftZmax = max(shiftZmax, globalZmax(ii))
      endif
    endif
  enddo
  !
  ! now initialize mesh loading if any meshes are used, then load
  ! shifted ones, if any, and set up to trim back to that
  !
  meshCheck = 0
  do kk = 1, numCheckTypes
    if (iobjFlag(abs(itypeCheck(kk))) == 4) meshCheck = 1
  enddo
  if (meshCheck == 1 .or. ishiftFlag == 4) then
    numMeshLoaded = 0
    numVerts = 0
    numTriang = 0
    numPoly = 0
    numSurf = 0
    numConts = 0
  endif
  if (ishiftFlag == 4) then
    ifZeroedShift = 0
    do kk = 1, numShiftTypes
      imodObj = abs(itypeShift(kk))
      call process_mesh(imodObj, xyScale, zScale, ibinSave, &
          zGapStart, zGapEnd, numGaps, ifErr, 1)
      if (ifErr .ne. 0) go to 99
      ifOnList = 0
      do imo = 1, numObjShifted
        if (iobjShift(imo) == itypeShift(kk)) ifOnList = imo
      enddo
      if (ifOnList > 0) then
        !
        ! if mesh is already on the list, shift surfaces back to zero
        ! shift and zero the shifts
        !
        ishift = indStartShift(ifOnList)
        do isurf = iobjSurf(numMeshLoaded), &
            iobjSurf(numMeshLoaded) + numSurfObj(numMeshLoaded) - 1
          call shiftSurface(isurf, shifts(1, ishift), -1, shifted, &
              ifZeroedShift, zScale)
          shifts(1, ishift) = 0.
          shifts(2, ishift) = 0.
          shifts(3, ishift) = 0.
          shifted(ishift) = .true.
          ishift = ishift + 1
        enddo
      else
        !
        ! otherwise, set up the right number of shifts for the points
        !
        numObjShifted = numObjShifted + 1
        indStartShift(numObjShifted) = ibaseShift
        iobjShift(numObjShifted) = itypeShift(kk)
        numItemsTemp = 0
        do isurf = iobjSurf(numMeshLoaded), &
            iobjSurf(numMeshLoaded) + numSurfObj(numMeshLoaded) - 1
          numItemsTemp = numItemsTemp + 1
          shifts(1, ibaseShift) = 0.
          shifts(2, ibaseShift) = 0.
          shifts(3, ibaseShift) = 0.
          shifted(ibaseShift) = .true.
          ibaseShift = ibaseShift + 1
          if (ibaseShift > LIMSHIFT) then
            print *,'Too many shifts for arrays'
            return
          endif
        enddo
        numItemShifted(numObjShifted) = numItemsTemp
        ifOnList = numObjShifted
      endif
      numItemsTot = numItemsTot + numItemShifted(ifOnList)
      !
      ! get mins/maxes
      !
      do isurf = iobjSurf(numMeshLoaded), &
          iobjSurf(numMeshLoaded) + numSurfObj(numMeshLoaded) - 1
        shiftXmin = min(shiftXmin, surfXmin(isurf))
        shiftYmin = min(shiftYmin, surfYmin(isurf))
        shiftZmin = min(shiftZmin, surfZmin(isurf))
        shiftXmax = max(shiftXmax, surfXmax(isurf))
        shiftYmax = max(shiftYmax, surfYmax(isurf))
        shiftZmax = max(shiftZmax, surfZmax(isurf))
      enddo
    enddo
  endif
  shiftXmin = shiftXmin - boxToler
  shiftXmax = shiftXmax + boxToler
  shiftYmin = shiftYmin - boxToler
  shiftYmax = shiftYmax + boxToler
  shiftZmin = shiftZmin - boxToler
  shiftZmax = shiftZmax + boxToler
  numMeshTrim = numMeshLoaded
  numVertsTrim = numVerts
  numTriangTrim = numTriang
  numPolyTrim = numPoly
  numSurfTrim = numSurf
  numContTrim = numConts
  !
  ! initialize the shift flag and count arrays
  !
  do ind = 1, numItemsTot
    needShift(ind) = .true.
    numTrials(ind) = 0
  enddo
  numFinal = 0
  !
  numFail = 0
  numFinal = 0
  iroundNum = 0
  numExclude = 0
  needTotal = 0

  do while(numFinal + numFail + numExclude < numItemsTot)
    if (iroundNum > 0 .and. manyRandom == 0) then
      write(*,'(a,i4,a,i5,a,i5,a,i5,a,$)') char(13) //'Round', iroundNum, &
          ':', numFinal, ' moved,', numFail, ' failed,', numExclude, &
          ' excluded'
      call flush(6)
    endif
    iroundNum = iroundNum + 1
    locShiftInd = 1
    do kk = 1, numShiftTypes
      ifOnList = 0
      do imo = 1, numObjShifted
        if (iobjShift(imo) == itypeShift(kk)) ifOnList = imo
      enddo
      ishift = indStartShift(ifOnList)
      if (ishiftFlag == 1) then
        !
        ! loop through the objects looking for ones to shift
        !
        do iobjRef = 1, numWobj
          !
          if (itypeObj(iobjRef) == itypeShift(kk)) then
            indRef = indStart(iobjRef)
            limRef = indRef + numPtInObj(iobjRef) - 2
            !
            ! on first trial, check if outside boundaries and fail
            ! DNM 7/3/02: this was a single if test, but SGI was calling
            ! outside_boudary even when nbound was zero
            !
            if (iroundNum == 1 .and. ifExcludeOut .ne. 0 .and. &
                numBound > 0) then
              if (outside_boundary(numBound, listBound, zBound, xModPt, yModPt, &
                  zModPt, indRef, limRef + 1, 0., 0., 0.)) then
                shifted(ishift) = .false.
                needShift(locShiftInd) = .false.
                needCheck(locShiftInd) = .false.
                numExclude = numExclude + 1
              endif
            endif
            !
            do while(needShift(locShiftInd))
              !
              call getRanShift(numTrials(locShiftInd), numTrialCycle, cycleFactor, &
                  ranMin, ranMax, ranZrel, shifts(1, ishift), delXYZ, &
                  globalXmin(iobjRef), globalXmax(iobjRef), globalYmin(iobjRef), &
                  globalYmax(iobjRef), globalZmin(iobjRef), globalZmax(iobjRef), &
                  shiftXmin, shiftXmax, shiftYmin, shiftYmax, shiftZmin, &
                  shiftZmax, zScale, iseed, numBound, listBound, zBound, xModPt, yModPt, &
                  zModPt, indRef, limRef + 1)
              !
              ! now shift the line and all its components
              !
              if (delX .ne. 0. .or. delY .ne. 0. .or. delZ .ne. 0.) then
                call shiftLine(iobjRef, indRef, limRef, delX, delY, delZ, &
                    xModPt, yModPt, zModPt, xMin, xMax, yMin, yMax, zMin, zMax, &
                    globalXmin, globalXmax, globalYmin, globalYmax, globalZmin, &
                    globalZmax, zScale)
                refXmin = globalXmin(iobjRef)
                refYmin = globalYmin(iobjRef)
                refZmin = globalZmin(iobjRef)
                refXmax = globalXmax(iobjRef)
                refYmax = globalYmax(iobjRef)
                refZmax = globalZmax(iobjRef)
                !
                ! set up for loops through shifted and fixed objects
                ! skip through all and restore if at end of trials
                !
                passing = numTrials(locShiftInd) <= maxTrials
              else
                passing = .false.
              endif
              !
              ifCheckFixed = 0
              numCheckTemp = numShiftTypes
              indCheckShift = 1
              do while(passing .and. ifCheckFixed <= 1)
                icheck = 1
                do while(passing .and. icheck <= numCheckTemp)
                  if (ifCheckFixed == 0) then
                    itypeCheckTemp = itypeShift(icheck)
                    jprob = 1
                  else
                    itypeCheckTemp = itypeCheck(icheck)
                    jprob = indProbCurve(icheck)
                  endif
                  call get_distlims(numNearBins, deltaNear, probNear, &
                      limProbs, jprob, distLimit, distAbs, probUse, deltaUse)
                  if (iobjFlag(itypeCheckTemp) == 1) then
                    !
                    ! THE LINE AGAINST LINES
                    !
                    iobjNeigh = 1
                    do while(passing .and. iobjNeigh <= numWobj)
                      if (itypeObj(iobjNeigh) == itypeCheckTemp) then
                        distMin = 1.01 * distLimit
                        doCheck = &
                            globalXmin(iobjNeigh) - refXmax < distMin .and. &
                            refXmin - globalXmax(iobjNeigh) < distMin .and. &
                            globalYmin(iobjNeigh) - refYmax < distMin .and. &
                            refYmin - globalYmax(iobjNeigh) < distMin .and. &
                            globalZmin(iobjNeigh) - refZmax < distMin .and. &
                            refZmin - globalZmax(iobjNeigh) < distMin
                        if (ifCheckFixed == 0) then
                          doCheck = doCheck .and. iobjNeigh .ne. iobjRef
                          if (ifCheckUnshifted == 0) doCheck = doCheck .and. &
                              .not.needShift(indCheckShift)
                          indCheckShift = indCheckShift + 1
                        endif
                        if (doCheck) then
                          dminSq = distMin**2

                          indNeigh = indStart(iobjNeigh)
                          limNeigh = indNeigh + numPtInObj(iobjNeigh) - 2
                          ineigh = indNeigh
                          do while(passing .and. ineigh <= limNeigh)
                            if (xMin(ineigh) - refXmax < distMin .and. &
                                refXmin - xMax(ineigh) < distMin .and. &
                                yMin(ineigh) - refYmax < distMin .and. &
                                refYmin - yMax(ineigh) < distMin .and. &
                                zMin(ineigh) - refZmax < distMin .and. &
                                refZmin - zMax(ineigh) < distMin) then
                              iref = indRef
                              inay = ineigh
                              do while(passing .and. iref <= limRef)
                                if (xMin(inay) - xMax(iref) < distMin .and. &
                                    xMin(iref) - xMax(inay) < distMin .and. &
                                    yMin(inay) - yMax(iref) < distMin .and. &
                                    yMin(iref) - yMax(inay) < distMin .and. &
                                    zMin(inay) - zMax(iref) < distMin .and. &
                                    zMin(iref) - zMax(inay) < distMin) then
                                  call segment_dist(xModPt(iref), yModPt(iref), &
                                      zModPt(iref), xModPt(iref + 1), yModPt(iref + 1), &
                                      zModPt(iref + 1), xModPt(inay), yModPt(inay), &
                                      zModPt(inay), xModPt(inay + 1), &
                                      yModPt(inay + 1), zModPt(inay + 1), ttRef, &
                                      ttNeigh, dsqr)
                                  if (dsqr < dminSq) then
                                    dminSq = dsqr
                                    distMin = sqrt(dsqr)
                                    passing = distMin >= distAbs
                                  endif
                                endif
                                iref = iref + 1
                              enddo
                            endif
                            ineigh = ineigh + 1
                          enddo
                          !
                          ! got a distance now - see if need to reject
                          !
                          if (passing .and. distMin < distLimit) &
                              passing = b3dran(iseed) < &
                              probUse(int(distMin / deltaUse + 1.))
                        endif
                      endif
                      iobjNeigh = iobjNeigh + 1
                    enddo
                  elseif (iobjFlag(itypeCheckTemp) == 2) then
                    !
                    ! THE LINE AGAINST POINTS
                    !
                    iobjNeigh = 1
                    do while(passing .and. iobjNeigh <= numWobj)
                      if (itypeObj(iobjNeigh) == itypeCheckTemp) then
                        distMin = 1.01 * distLimit
                        doCheck = &
                            globalXmin(iobjNeigh) - refXmax < distMin .and. &
                            refXmin - globalXmax(iobjNeigh) < distMin .and. &
                            globalYmin(iobjNeigh) - refYmax < distMin .and. &
                            refYmin - globalYmax(iobjNeigh) < distMin .and. &
                            globalZmin(iobjNeigh) - refZmax < distMin .and. &
                            refZmin - globalZmax(iobjNeigh) < distMin
                        if (doCheck) then

                          indNeigh = indStart(iobjNeigh)
                          limNeigh = indNeigh + numPtInObj(iobjNeigh) - 2
                          inay = indNeigh
                          do while(passing .and. inay <= limNeigh)
                            distMin = 1.01 * distLimit
                            if (xMin(inay) - refXmax < distMin .and. &
                                refXmin - xMax(inay) < distMin .and. &
                                yMin(inay) - refYmax < distMin .and. &
                                refYmin - yMax(inay) < distMin .and. &
                                zMin(inay) - refZmax < distMin .and. &
                                refZmin - zMax(inay) < distMin) then

                              sizeNeigh = 0.
                              if (ifUseScatSurf .ne. 0) then
                                is = indSize(iobjNeigh) + inay - indNeigh
                                sizeNeigh = sizes(is)
                              endif

                              iref = indRef
                              do while(passing .and. iref <= limRef)
                                if (xMin(inay) - xMax(iref) < distMin .and. &
                                    xMin(iref) - xMax(inay) < distMin .and. &
                                    yMin(inay) - yMax(iref) < distMin .and. &
                                    yMin(iref) - yMax(inay) < distMin .and. &
                                    zMin(inay) - zMax(iref) < distMin .and. &
                                    zMin(iref) - zMax(inay) < distMin) then
                                  call point_line_dist(xModPt(iref), &
                                      yModPt(iref), zModPt(iref), &
                                      xModPt(iref + 1), yModPt(iref + 1), &
                                      zModPt(iref + 1), xModPt(inay), yModPt(inay), &
                                      zModPt(inay), ttRef, dsqr)
                                  dist = sqrt(dsqr) - sizeNeigh
                                  if (dist < distMin) then
                                    distMin = dist
                                    passing = distMin >= distAbs
                                  endif
                                endif
                                iref = iref + 1
                              enddo
                              !
                              ! got a distance now - see if reject
                              !
                              if (passing .and. distMin < distLimit) &
                                  passing = b3dran(iseed) < &
                                  probUse(int(distMin / deltaUse + 1.))
                            endif
                            inay = inay + 1
                          enddo
                        endif
                      endif
                      iobjNeigh = iobjNeigh + 1
                    enddo
                  endif
                  icheck = icheck + 1
                enddo
                ifCheckFixed = ifCheckFixed + 1
                numCheckTemp = numCheckTypes
              enddo
              !
              ! got through all stage 1 tests: if passed, set flag; if
              ! not, test for trial limit and restore to zero shift
              !
              if (passing) then
                needShift(locShiftInd) = .false.
                needCheck(locShiftInd) = .true.
                needTotal = needTotal + 1
              elseif (numTrials(locShiftInd) >= maxTrials) then
                call shiftLine(iobjRef, indRef, limRef, -shifts(1, ishift), &
                    - shifts(2, ishift), -shifts(3, ishift), &
                    xModPt, yModPt, zModPt, xMin, xMax, yMin, yMax, zMin, zMax, &
                    globalXmin, globalXmax, globalYmin, globalYmax, globalZmin, &
                    globalZmax, zScale)
                shifts(1, ishift) = 0.
                shifts(2, ishift) = 0.
                shifts(3, ishift) = 0.
                shifted(ishift) = .false.
                needShift(locShiftInd) = .false.
                needCheck(locShiftInd) = .false.
                numFail = numFail + 1
              endif
            enddo
            !
            ! advance to next item to shift
            !
            ishift = ishift + 1
            locShiftInd = locShiftInd + 1
          endif
        enddo
      elseif (ishiftFlag == 2) then
        !
        ! POINTS TO SHIFT: LOOP THROUGH OBJECTS
        !
        do iobjRef = 1, numWobj
          if (itypeObj(iobjRef) == itypeShift(kk)) then
            !
            indRef = indStart(iobjRef)
            limRef = indRef + numPtInObj(iobjRef) - 1
            do iref = indRef, limRef
              !
              ! on first trial, check if outside boundaries and fail
              !
              if (iroundNum == 1 .and. ifExcludeOut .ne. 0 .and. numBound > 0) then
                if (outside_boundary(numBound, listBound, zBound, xModPt, yModPt, &
                    zModPt, iref, iref, 0., 0., 0.)) then
                  shifted(ishift) = .false.
                  needShift(locShiftInd) = .false.
                  needCheck(locShiftInd) = .false.
                  numExclude = numExclude + 1
                endif
              endif
              !
              do while(needShift(locShiftInd))
                !
                call getRanShift(numTrials(locShiftInd), numTrialCycle, cycleFactor, &
                    ranMin, ranMax, ranZrel, shifts(1, ishift), delXYZ, &
                    xMin(iref), xMax(iref), yMin(iref), &
                    yMax(iref), zMin(iref), zMax(iref), &
                    shiftXmin, shiftXmax, shiftYmin, shiftYmax, shiftZmin, &
                    shiftZmax, zScale, iseed, numBound, listBound, zBound, &
                    xModPt, yModPt, zModPt, iref, iref)
                !
                ! now shift the point and all its components
                !
                if (delX .ne. 0. .or. delY .ne. 0. .or. delZ .ne. 0.) then
                  call shiftPoint(iobjRef, iref, delX, delY, delZ, &
                      xModPt, yModPt, zModPt, xMin, xMax, yMin, yMax, zMin, zMax, &
                      globalXmin, globalXmax, globalYmin, globalYmax, globalZmin, &
                      globalZmax, zScale)
                  refXmin = xMin(iref)
                  refYmin = yMin(iref)
                  refZmin = zMin(iref)
                  refXmax = xMax(iref)
                  refYmax = yMax(iref)
                  refZmax = zMax(iref)
                  sizeRef = 0.
                  if (ifUseScatSurf .ne. 0) then
                    is = indSize(iobjRef) + iref - indRef
                    sizeRef = sizes(is)
                  endif
                  !
                  ! set up for loops through shifted and fixed objects
                  ! skip through all and restore if at end of trials
                  !
                  passing = numTrials(locShiftInd) <= maxTrials
                else
                  passing = .false.
                endif
                !
                ifCheckFixed = 0
                numCheckTemp = numShiftTypes
                indCheckShift = 1
                do while(passing .and. ifCheckFixed <= 1)
                  icheck = 1
                  do while(passing .and. icheck <= numCheckTemp)
                    if (ifCheckFixed == 0) then
                      itypeCheckTemp = itypeShift(icheck)
                      jprob = 1
                    else
                      itypeCheckTemp = itypeCheck(icheck)
                      jprob = indProbCurve(icheck)
                    endif
                    call get_distlims(numNearBins, deltaNear, probNear, &
                        limProbs, jprob, distLimit, distAbs, probUse, deltaUse)
                    if (iobjFlag(itypeCheckTemp) == 1) then
                      !
                      ! THE POINT AGAINST LINES
                      !
                      iobjNeigh = 1
                      do while(passing .and. iobjNeigh <= numWobj)
                        if (itypeObj(iobjNeigh) == itypeCheckTemp) then
                          distMin = 1.01 * distLimit
                          doCheck = &
                              globalXmin(iobjNeigh) - refXmax < distMin .and. &
                              refXmin - globalXmax(iobjNeigh) < distMin .and. &
                              globalYmin(iobjNeigh) - refYmax < distMin .and. &
                              refYmin - globalYmax(iobjNeigh) < distMin .and. &
                              globalZmin(iobjNeigh) - refZmax < distMin .and. &
                              refZmin - globalZmax(iobjNeigh) < distMin
                          if (doCheck) then
                            dminSq = distMin**2

                            indNeigh = indStart(iobjNeigh)
                            limNeigh = indNeigh + numPtInObj(iobjNeigh) - 2
                            inay = indNeigh
                            do while(passing .and. inay <= limNeigh)
                              if (xMin(inay) - refXmax < distMin .and. &
                                  refXmin - xMax(inay) < distMin .and. &
                                  yMin(inay) - refYmax < distMin .and. &
                                  refYmin - yMax(inay) < distMin .and. &
                                  zMin(inay) - refZmax < distMin .and. &
                                  refZmin - zMax(inay) < distMin) then
                                call point_line_dist(xModPt(inay), yModPt(inay), &
                                    zModPt(inay), xModPt(inay + 1), yModPt(inay + 1), &
                                    zModPt(inay + 1), xModPt(iref), yModPt(iref), &
                                    zModPt(iref), ttRef, dsqr)
                                dist = sqrt(dsqr) - sizeRef
                                if (dist < distMin) then
                                  distMin = dist
                                  passing = distMin >= distAbs
                                endif
                              endif
                              inay = inay + 1
                            enddo
                            !
                            ! got a distance now - see if need to reject
                            !
                            if (passing .and. distMin < distLimit) &
                                passing = b3dran(iseed) < &
                                probUse(int(distMin / deltaUse + 1.))
                          endif
                        endif
                        iobjNeigh = iobjNeigh + 1
                      enddo

                    elseif (iobjFlag(itypeCheckTemp) == 2) then
                      !
                      ! THE POINT AGAINST POINTS
                      !
                      iobjNeigh = 1
                      do while(passing .and. iobjNeigh <= numWobj)
                        distMin = 1.01 * distLimit
                        if (itypeObj(iobjNeigh) == itypeCheckTemp .and. &
                            globalXmin(iobjNeigh) - refXmax < distMin .and. &
                            refXmin - globalXmax(iobjNeigh) < distMin .and. &
                            globalYmin(iobjNeigh) - refYmax < distMin .and. &
                            refYmin - globalYmax(iobjNeigh) < distMin .and. &
                            globalZmin(iobjNeigh) - refZmax < distMin .and. &
                            refZmin - globalZmax(iobjNeigh) < distMin) then
                          indNeigh = indStart(iobjNeigh)
                          limNeigh = indNeigh + numPtInObj(iobjNeigh) - 1
                          inay = indNeigh
                          if (ifUseScatSurf .ne. 0) &
                              isBase = indSize(iobjNeigh) - indNeigh
                          do while(passing .and. inay <= limNeigh)
                            doCheck = xMin(inay) - refXmax < distMin .and. &
                                refXmin - xMax(inay) < distMin .and. &
                                yMin(inay) - refYmax < distMin .and. &
                                refYmin - yMax(inay) < distMin .and. &
                                zMin(inay) - refZmax < distMin .and. &
                                refZmin - zMax(inay) < distMin
                            if (ifCheckFixed == 0) then
                              doCheck = doCheck .and. inay .ne. iref
                              if (ifCheckUnshifted == 0) doCheck = &
                                  doCheck .and. .not.needShift(indCheckShift)
                              indCheckShift = indCheckShift + 1
                            endif
                            if (doCheck) then
                              sizeNeigh = 0.
                              if (ifUseScatSurf .ne. 0) sizeNeigh = sizes(inay + isBase)
                              dist = sqrt((xModPt(inay) - xModPt(iref))**2 + &
                                  (yModPt(inay) - yModPt(iref))**2 + &
                                  (zModPt(inay) - zModPt(iref))**2) - sizeNeigh - sizeRef
                              !
                              ! got a distance now - see if reject
                              !
                              passing = dist >= distAbs
                              if (passing .and. dist < distLimit) &
                                  passing = b3dran(iseed) < &
                                  probUse(int(dist / deltaUse + 1.))
                            endif
                            inay = inay + 1
                          enddo
                        endif
                        iobjNeigh = iobjNeigh + 1
                      enddo
                    endif
                    icheck = icheck + 1
                  enddo
                  ifCheckFixed = ifCheckFixed + 1
                  numCheckTemp = numCheckTypes
                enddo
                !
                ! got through all stage 1 tests: if passed, set flag; if
                ! not, test for trial limit and restore to zero shift
                !
                if (passing) then
                  needShift(locShiftInd) = .false.
                  needCheck(locShiftInd) = .true.
                  needTotal = needTotal + 1
                elseif (numTrials(locShiftInd) >= maxTrials) then
                  call shiftPoint(iobjRef, iref, -shifts(1, ishift), &
                      - shifts(2, ishift), -shifts(3, ishift), &
                      xModPt, yModPt, zModPt, xMin, xMax, yMin, yMax, zMin, zMax, &
                      globalXmin, globalXmax, globalYmin, globalYmax, globalZmin, &
                      globalZmax, zScale)
                  shifts(1, ishift) = 0.
                  shifts(2, ishift) = 0.
                  shifts(3, ishift) = 0.
                  shifted(ishift) = .false.
                  needShift(locShiftInd) = .false.
                  needCheck(locShiftInd) = .false.
                  numFail = numFail + 1
                endif
              enddo
              !
              ! advance to next item to shift
              !
              ishift = ishift + 1
              locShiftInd = locShiftInd + 1
            enddo
          endif
        enddo
      else
        !
        ! MESHES TO SHIFT: find the mesh in the loaded ones
        !
        do ind = 1, numMeshLoaded
          if (iobjMesh(ind) == abs(itypeShift(kk))) indMesh = ind
        enddo
        !
        ! loop through surfaces
        !
        do isurf = iobjSurf(indMesh), iobjSurf(indMesh) + numSurfObj(indMesh) - 1
          !
          ! on first trial, check if outside boundaries and fail
          !
          if (iroundNum == 1 .and. ifExcludeOut .ne. 0 .and. numBound > 0) then
            if (outside_boundary(numBound, listBound, zBound, xModPt, yModPt, &
                zModPt, isurf, 0, 0., 0., 0.)) then
              shifted(ishift) = .false.
              needShift(locShiftInd) = .false.
              needCheck(locShiftInd) = .false.
              numExclude = numExclude + 1
            endif
          endif
          !
          do while(needShift(locShiftInd))
            !
            call getRanShift(numTrials(locShiftInd), numTrialCycle, cycleFactor, &
                ranMin, ranMax, ranZrel, shifts(1, ishift), delXYZ, &
                surfXmin(isurf), surfXmax(isurf), surfYmin(isurf), &
                surfYmax(isurf), surfZmin(isurf), surfZmax(isurf), &
                shiftXmin, shiftXmax, shiftYmin, shiftYmax, shiftZmin, &
                shiftZmax, zScale, iseed, numBound, listBound, zBound, &
                xModPt, yModPt, zModPt, isurf, 0)
            !
            ! now shift the surface and all its components
            !
            if (delX .ne. 0. .or. delY .ne. 0. .or. delZ .ne. 0.) then
              ifZeroedShift = 0
              call shiftSurface(isurf, delXYZ, 1, ibinSave, ifZeroedShift, zScale)
              refXmin = surfXmin(isurf)
              refYmin = surfYmin(isurf)
              refZmin = surfZmin(isurf)
              refXmax = surfXmax(isurf)
              refYmax = surfYmax(isurf)
              refZmax = surfZmax(isurf)
              !
              ! set up for loops through shifted and fixed objects
              ! skip through all and restore if at end of trials
              !
              passing = numTrials(locShiftInd) <= maxTrials
            else
              passing = .false.
            endif
            !
            icheck = 1
            do while(passing .and. icheck <= numCheckTypes)
              jprob = indProbCurve(icheck)
              call get_distlims(numNearBins, deltaNear, probNear, &
                  limProbs, jprob, distLimit, distAbs, probUse, deltaUse)
              itypeCheckTemp = itypeCheck(icheck)
              iobjNeigh = 1
              if (iobjFlag(itypeCheckTemp) == 4) iobjNeigh = numWobj + 1
              do while(passing .and. iobjNeigh <= numWobj)
                distMin = 1.01 * distLimit
                indNeigh = indStart(iobjNeigh)
                if (itypeObj(iobjNeigh) == itypeCheckTemp .and. &
                    globalXmin(iobjNeigh) - refXmax < distMin .and. &
                    refXmin - globalXmax(iobjNeigh) < distMin .and. &
                    globalYmin(iobjNeigh) - refYmax < distMin .and. &
                    refYmin - globalYmax(iobjNeigh) < distMin .and. &
                    globalZmin(iobjNeigh) - refZmax < distMin .and. &
                    refZmin - globalZmax(iobjNeigh) < distMin) then
                  if (iobjFlag(itypeCheckTemp) == 1) then
                    !
                    ! THE MESH AGAINST A LINE
                    !
                    limNeigh = indNeigh + numPtInObj(iobjNeigh) - 2
                    call check_line_surface(isurf, indNeigh, limNeigh, &
                        xModPt, yModPt, zModPt, &
                        distMin, distAbs, zScale)
                    passing = distMin >= distAbs
                    if (passing .and. distMin < distLimit) &
                        passing = b3dran(iseed) < &
                        probUse(int(distMin / deltaUse + 1.))
                  else
                    !
                    ! THE MESH AGAINST A SET OF POINTS
                    !
                    limNeigh = indNeigh + numPtInObj(iobjNeigh) - 1
                    inay = indNeigh
                    if (ifUseScatSurf .ne. 0) isBase = indSize(iobjNeigh) - indNeigh
                    do while(passing .and. inay <= limNeigh)
                      distMin = 1.01 * distLimit
                      if (xMin(inay) - refXmax < distMin .and. &
                          refXmin - xMax(inay) < distMin .and. &
                          yMin(inay) - refYmax < distMin .and. &
                          refYmin - yMax(inay) < distMin .and. &
                          zMin(inay) - refZmax < distMin .and. &
                          refZmin - zMax(inay) < distMin) then
                        sizeNeigh = 0.
                        if (ifUseScatSurf .ne. 0) sizeNeigh = sizes(inay + isBase)
                        call check_point_surface(isurf, xModPt(inay), yModPt(inay), &
                            zModPt(inay), xMin(inay), xMax(inay), yMin(inay),  &
                            yMax(inay), zMin(inay), zMax(inay), sizeNeigh, distMin, &
                            distAbs, zScale)
                        passing = distMin >= distAbs
                        if (passing .and. distMin < distLimit) &
                            passing = b3dran(iseed) < &
                            probUse(int(distMin / deltaUse + 1.))
                      endif
                      inay = inay + 1
                    enddo
                  endif
                endif
                iobjNeigh = iobjNeigh + 1
              enddo
              icheck = icheck + 1
            enddo
            !
            ! NOW CHECK THIS MESH AGAINST OTHER SHIFTED MESHES
            !
            icheck = 1
            indCheckShift = 1
            jprob = 1
            call get_distlims(numNearBins, deltaNear, probNear, &
                limProbs, jprob, distLimit, distAbs, probUse, deltaUse)
            do while(passing .and. icheck <= numShiftTypes)
              do ind = 1, numMeshLoaded
                if (iobjMesh(ind) == abs(itypeShift(icheck))) jmesh = ind
              enddo
              jsurf = iobjSurf(jmesh)
              do while(passing .and. &
                  jsurf <= iobjSurf(jmesh) + numSurfObj(jmesh) - 1)
                distMin = 1.01 * distLimit
                doCheck = surfXmin(jsurf) - refXmax < distMin .and. &
                    refXmin - surfXmax(jsurf) < distMin .and. &
                    surfYmin(jsurf) - refYmax < distMin .and. &
                    refYmin - surfYmax(jsurf) < distMin .and. &
                    surfZmin(jsurf) - refZmax < distMin .and. &
                    refZmin - surfZmax(jsurf) < distMin
                doCheck = doCheck .and. jsurf .ne. isurf
                if (ifCheckUnshifted == 0) doCheck = doCheck .and. &
                    .not.needShift(indCheckShift)
                indCheckShift = indCheckShift + 1

                if (doCheck) then
                  call check_two_meshes(isurf, jsurf, distMin, &
                      distAbs, zScale)
                  passing = distMin >= distAbs
                  if (passing .and. distMin < distLimit) &
                      passing = b3dran(iseed) < &
                      probUse(int(distMin / deltaUse + 1.))
                endif
                jsurf = jsurf + 1
              enddo

              icheck = icheck + 1
            enddo
            ! got through all stage 1 tests: if passed, set flag; if
            ! not, test for trial limit and restore to zero shift
            !
            if (passing) then
              needShift(locShiftInd) = .false.
              needCheck(locShiftInd) = .true.
              needTotal = needTotal + 1
            elseif (numTrials(locShiftInd) >= maxTrials) then
              ifZeroedShift = 0
              call shiftSurface(isurf, shifts(1, ishift), -1, ibinSave, &
                  ifZeroedShift, zScale)
              shifts(1, ishift) = 0.
              shifts(2, ishift) = 0.
              shifts(3, ishift) = 0.
              shifted(ishift) = .false.
              needShift(locShiftInd) = .false.
              needCheck(locShiftInd) = .false.
              numFail = numFail + 1
            endif
          enddo
          !
          ! advance to next item to shift
          !
          ishift = ishift + 1
          locShiftInd = locShiftInd + 1
        enddo
      endif
    enddo
    !
    ! FINISHED STAGE 1, NOW LOOP ON THE MESHES TO BE CHECKED
    !
    do icheck = 1, numCheckTypes
      jprob = indProbCurve(icheck)
      call get_distlims(numNearBins, deltaNear, probNear, &
          limProbs, jprob, distLimit, distAbs, probUse, deltaUse)
      itypeCheckTemp = itypeCheck(icheck)
      if (iobjFlag(itypeCheckTemp) == 4 .and. needTotal > 0) then
        indMesh = 0
        imodObj = abs(itypeCheckTemp)
        do ind = 1, numMeshLoaded
          if (iobjMesh(ind) == imodObj) indMesh = ind
        enddo
        if (indMesh == 0) then
          call process_mesh(imodObj, xyScale, zScale, ibinSave, &
              zGapStart, zGapEnd, numGaps, ifErr, 1)
          if (ifErr .ne. 0) then
            !
            ! if not enough space, trim back, try again
            !
            numMeshLoaded = numMeshTrim
            numVerts = numVertsTrim
            numTriang = numTriangTrim
            numPoly = numPolyTrim
            numSurf = numSurfTrim
            numConts = numContTrim
            call process_mesh(imodObj, xyScale, zScale, ibinSave, &
                zGapStart, zGapEnd, numGaps, ifErr, 1)
            if (ifErr .ne. 0) go to 99
          endif
          indMesh = numMeshLoaded
        endif
        !
        ! go through surfaces checking for distances
        !
        do isurf = iobjSurf(indMesh), iobjSurf(indMesh) + numSurfObj(indMesh) - 1
          refXmin = surfXmin(isurf)
          refYmin = surfYmin(isurf)
          refZmin = surfZmin(isurf)
          refXmax = surfXmax(isurf)
          refYmax = surfYmax(isurf)
          refZmax = surfZmax(isurf)
          locShiftInd = 1
          do kk = 1, numShiftTypes
            ifOnList = 0
            do imo = 1, numObjShifted
              if (iobjShift(imo) == itypeShift(kk)) ifOnList = imo
            enddo
            ishift = indStartShift(ifOnList)
            if (ishiftFlag == 1) then
              !
              ! CHECK LINES AGAINST THE SURFACE
              !
              do iobjRef = 1, numWobj
                if (itypeObj(iobjRef) == itypeShift(kk)) then
                  if (needCheck(locShiftInd)) then
                    indRef = indStart(iobjRef)
                    limRef = indRef + numPtInObj(iobjRef) - 2
                    distMin = 1.01 * distLimit
                    if (globalXmin(iobjRef) - refXmax < distMin .and. &
                        refXmin - globalXmax(iobjRef) < distMin .and. &
                        globalYmin(iobjRef) - refYmax < distMin .and. &
                        refYmin - globalYmax(iobjRef) < distMin .and. &
                        globalZmin(iobjRef) - refZmax < distMin .and. &
                        refZmin - globalZmax(iobjRef) < distMin) then
                      call check_line_surface(isurf, indRef, limRef, &
                          xModPt, yModPt, zModPt, &
                          distMin, distAbs, zScale)
                      passing = distMin >= distAbs
                      if (passing .and. distMin < distLimit) &
                          passing = b3dran(iseed) < &
                          probUse(int(distMin / deltaUse + 1.))
                      if (.not.passing) then
                        !
                        ! If it doesn't pass, toggle the two flags and
                        ! shift the line back to zero shift
                        !
                        needCheck(locShiftInd) = .false.
                        needShift(locShiftInd) = .true.
                        needTotal = needTotal - 1
                        call shiftLine(iobjRef, indRef, limRef, &
                            - shifts(1, ishift), -shifts(2, ishift), &
                            - shifts(3, ishift), xModPt, yModPt, zModPt, &
                            xMin, xMax, yMin, yMax, zMin, zMax, globalXmin, &
                            globalXmax, globalYmin, globalYmax, globalZmin, &
                            globalZmax, zScale)
                        shifts(1, ishift) = 0.
                        shifts(2, ishift) = 0.
                        shifts(3, ishift) = 0.
                      endif
                    endif
                  endif
                  ishift = ishift + 1
                  locShiftInd = locShiftInd + 1
                endif
              enddo
            elseif (ishiftFlag == 2) then
              !
              ! CHECK POINTS AGAINST THE SURFACE
              !
              do iobjRef = 1, numWobj
                if (itypeObj(iobjRef) == itypeShift(kk)) then
                  indRef = indStart(iobjRef)
                  limRef = indRef + numPtInObj(iobjRef) - 1
                  distMin = 1.01 * distLimit
                  if (globalXmin(iobjRef) - refXmax < distMin .and. &
                      refXmin - globalXmax(iobjRef) < distMin .and. &
                      globalYmin(iobjRef) - refYmax < distMin .and. &
                      refYmin - globalYmax(iobjRef) < distMin .and. &
                      globalZmin(iobjRef) - refZmax < distMin .and. &
                      refZmin - globalZmax(iobjRef) < distMin) then
                    if (ifUseScatSurf .ne. 0) isBase = indSize(iobjRef) - indRef

                    do iref = indRef, limRef
                      if (needCheck(locShiftInd)) then
                        distMin = 1.01 * distLimit
                        if (xMin(iref) - refXmax < distMin .and. &
                            refXmin - xMax(iref) < distMin .and. &
                            yMin(iref) - refYmax < distMin .and. &
                            refYmin - yMax(iref) < distMin .and. &
                            zMin(iref) - refZmax < distMin .and. &
                            refZmin - zMax(iref) < distMin) then
                          sizeRef = 0.
                          if (ifUseScatSurf .ne. 0) sizeRef = sizes(iref + isBase)

                          call check_point_surface(isurf, xModPt(iref), &
                              yModPt(iref), zModPt(iref), &
                              xMin(iref), xMax(iref), yMin(iref), &
                              yMax(iref), zMin(iref), zMax(iref), &
                              sizeRef, distMin, distAbs, zScale)
                          passing = distMin >= distAbs
                          if (passing .and. distMin < distLimit) &
                              passing = b3dran(iseed) < &
                              probUse(int(distMin / deltaUse + 1.))
                          if (.not.passing) then
                            !
                            ! if it doesn't pass, shift point back to 0
                            !
                            needCheck(locShiftInd) = .false.
                            needShift(locShiftInd) = .true.
                            needTotal = needTotal - 1
                            call shiftPoint(iobjRef, iref, &
                                - shifts(1, ishift), -shifts(2, ishift), &
                                - shifts(3, ishift), xModPt, yModPt, zModPt, &
                                xMin, xMax, yMin, yMax, zMin, zMax, globalXmin, &
                                globalXmax, globalYmin, globalYmax, globalZmin, &
                                globalZmax, zScale)
                            shifts(1, ishift) = 0.
                            shifts(2, ishift) = 0.
                            shifts(3, ishift) = 0.
                          endif
                        endif
                      endif
                      locShiftInd = locShiftInd + 1
                      ishift = ishift + 1
                    enddo
                  else
                    !
                    ! If skip over a whole contour, advance both indexes
                    ! to account for all the points skipped
                    !
                    locShiftInd = locShiftInd + numPtInObj(iobjRef)
                    ishift = ishift + numPtInObj(iobjRef)
                  endif
                endif
              enddo
            else
              !
              ! CHECK OTHER SURFACES AGAINST THIS SURFACE
              !
              ifZeroedShift = 0
              do ind = 1, numMeshLoaded
                if (iobjMesh(ind) == itypeShift(kk)) jmesh = ind
              enddo
              do jsurf = iobjSurf(jmesh), iobjSurf(jmesh) + numSurfObj(jmesh) - 1
                distMin = 1.01 * distLimit
                doCheck = needCheck(locShiftInd) .and. &
                    surfXmin(jsurf) - refXmax < distMin .and. &
                    refXmin - surfXmax(jsurf) < distMin .and. &
                    surfYmin(jsurf) - refYmax < distMin .and. &
                    refYmin - surfYmax(jsurf) < distMin .and. &
                    surfZmin(jsurf) - refZmax < distMin .and. &
                    refZmin - surfZmax(jsurf) < distMin

                if (doCheck) then
                  call check_two_meshes(isurf, jsurf, distMin, &
                      distAbs, zScale)
                  passing = distMin >= distAbs
                  if (passing .and. distMin < distLimit) &
                      passing = b3dran(iseed) < &
                      probUse(int(distMin / deltaUse + 1.))
                  if (.not.passing) then
                    needCheck(locShiftInd) = .false.
                    needShift(locShiftInd) = .true.
                    needTotal = needTotal - 1
                    call shiftSurface(jsurf, shifts(1, ishift), -1, ibinSave, &
                        ifZeroedShift, zScale)
                    shifts(1, ishift) = 0.
                    shifts(2, ishift) = 0.
                    shifts(3, ishift) = 0.
                  endif
                endif
                locShiftInd = locShiftInd + 1
                ishift = ishift + 1
              enddo
            endif
          enddo
        enddo
      endif
    enddo
    !
    ! finalize anyone still listed as "needcheck"
    !
    do ind = 1, numItemsTot
      if (needCheck(ind)) then
        needCheck(ind) = .false.
        numFinal = numFinal + 1
        ! write(*,'(2i5,3f8.4)') i, ntrials(i), (shifts(j, i), j=1, 3)
      endif
    enddo
  enddo
  !
  numTrialTot = 0
  numItemSum = 0
  do ind = 1, numItemsTot
    if (numTrials(ind) > 0) then
      numItemSum = numItemSum + 1
      numTrialTot = numTrialTot + numTrials(ind)
    endif
  enddo
  avgNumTrial = float(numTrialTot) / numItemSum
  !
  write(*,104) numFinal, numFail, numExclude, avgNumTrial
104 format(i5,' items shifted;',i5,' could not be;',i5, &
      ' excluded;',f7.1, ' trials/item')
  return
99 print *,'Unable to load meshes'
  return
end subroutine random_shifts

subroutine get_distlims(numNearBins, deltaNear, probNear, &
    limProbs, jprob, distLimit, distAbs, probUse, deltaUse)
  real*4 probNear(limProbs,*), probUse(*), deltaNear(*)
  integer*4 numNearBins(*)
  distLimit = deltaNear(jprob) * numNearBins(jprob)
  deltaUse = deltaNear(jprob)
  distAbs = 0.
  ind = 1
  do while (ind <= numNearBins(jprob) .and. probNear(ind, jprob) == 0)
    distAbs = ind * deltaUse
    ind = ind + 1
  enddo
  do ind = 1, numNearBins(jprob)
    probUse(ind) = probNear(ind, jprob)
  enddo
  return
end subroutine get_distlims
