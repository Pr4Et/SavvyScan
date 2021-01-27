! SUBROUTINES FOR MTK ONLY
!
! $Id$
!
! GETBINSPEC gets a bin specification appropriate for the type of
! graphs being done
!
subroutine getBinSpec(deltaRad, numBins, power, limFit, padBound, fracOmit, &
    ifBundEnd, sampleLen, ifCloseSeg, ifUseScatSurf)
  implicit none
  integer LIMBINS
  parameter (LIMBINS = 1001)
  real*4 deltaRad, power, padBound, fracOmit, sampleLen
  integer*4 numBins, ifChange, ifBundEnd, limFit, ifCloseSeg, ifUseScatSurf
  logical didOnce / .false./
  save didOnce
  integer*4 in5
  common /nmsInput/ in5
  !
  write(*,'(1x,a,$)') &
      'Bin width (radial distance), number of bins: '
  read(in5,*) deltaRad, numBins
  if (numBins >= LIMBINS) print *,'# of bins truncated to', LIMBINS - 1
  numBins = min(numBins, LIMBINS - 1)
  !
  ifChange = 1
  if (didOnce) then
    write(*,'(1x,a,$)') '0 to keep other parameters the same,'// &
        ' 1 to change them: '
    read(in5,*) ifChange
  endif
  didOnce = .true.
  if (ifBundEnd == 0 .and. ifChange .ne. 0) then
    write(*,'(1x,a,$)') 'Sampling length for lines, or 0 to find ' &
        //'closest approach to whole line: '
    read(in5,*) sampleLen
    write(*,'(1x,a,/,a,$)') '(For whole lines) Power for radial ' &
        //'weighting,', '  number of points to fit over: '
    read(in5,*) power, limFit
    limFit = max(1, limFit)
    write(*,'(1x,a,/,a,$)') '(For line segments) 0 to find '// &
        'distance from start of sample segment,', &
        '  or 1 to find closest approach to segment: '
    read(in5,*) ifCloseSeg
    write(*,'(1x,a,$)') '0 to measure from center of scattered ' &
        //'points or 1 to measure from surface: '
    read(in5,*) ifUseScatSurf
  elseif (ifChange .ne. 0) then
    write(*,'(1x,a,/,a,$)') 'Enter distance to pad boundaries (- ' &
        //'if in pixels to be scaled),', '   and fraction '// &
        'of farthest points to omit from bundle: '
    read(in5,*) padBound, fracOmit
    write(*,'(1x,a,$)') '0 or 1 to analyze distances to all '// &
        'bundles or only to nearest bundle: '
    read(in5,*) limFit
  endif
  return
end subroutine getBinSpec



! GETGRAPHSPEC gets specifications of the types involved in each
! desired graph.
! NGRAPH is returned with the number of graphs
! NxxxTYP is returned with the number of types for reference items
! (xxx=REF), or neighboring items (xxx=NEIGH) ;
! ITYPxxx is returned with the list of types for each case
!
subroutine getGraphSpec(numGraphs, itypeRef, &
    numRefType, itypeNeigh, numNeighType, iwhichEnd, ifBundEnd, &
    iobjFlag, limFlag, irefFlag, neighFlag)
  parameter (limType = 50, itypeAll = 999)
  integer*4 numRefType(*), numNeighType(*)         !# of types for ref and neigh
  integer*4 itypeRef(limType,*), itypeNeigh(limType,*)
  integer*4 iwhichEnd(limType,*), iobjFlag(*)
  logical didOnce / .false./
  save didOnce
  integer*4 in5
  common /nmsInput/ in5
  !
  ifChange = 1
  if (didOnce) then
    write(*,'(1x,a,$)') '0 to keep same graph specifications or'// &
        ' 1 to specify new graphs: '
    read(in5,*) ifChange
  endif
  didOnce = .true.
10 if (ifChange .ne. 0) then
    irefFlag = -1
    neighFlag = -1
    write(*,'(1x,a,$)') 'Number of different graphs to compute: '
    read(in5,*) numGraphs
    !
    do ii = 1, numGraphs
      write(*,102) ii, 'from (reference'
102   format(' For graph #',i3,', enter list of objects to ', &
          'measure distances',/,5x,a,' objects)', &
          ' (Return for all, ranges OK)')
      call rdlist(in5, itypeRef(1, ii), numRefType(ii))
      if (numRefType(ii) == 0) then
        numRefType(ii) = 1
        itypeRef(1, ii) = itypeAll
      endif
      !
      write(*,102) ii, 'to (neighboring'
      call rdlist(in5, itypeNeigh(1, ii), numNeighType(ii))
      if (numNeighType(ii) == 0) then
        numNeighType(ii) = 1
        itypeNeigh(1, ii) = itypeAll
      endif
      !
      if (ifBundEnd .ne. 0) then
        print *,'For each neighbor object, enter 0 to count ', &
            'neither end, 1 to count the end', '    at low Z, 2 to' &
            , ' count the end at high Z, or 3 to count both ends'
        read(in5,*) (iwhichEnd(j, ii), j = 1, numNeighType(ii))
      endif
    enddo
  endif
  call checkFlags(numGraphs, itypeRef, numRefType, itypeNeigh, &
      numNeighType, iobjFlag, limFlag, irefFlag, neighFlag)
  if (irefFlag >= 0 .and. neighFlag >= 0) return
  if (ifChange .ne. 0) then
    print *,'Start over with graph specifications'
  else
    print *,'You need to enter new graph specifications'
  endif
  go to 10
end subroutine getGraphSpec


subroutine checkFlags(numGraphs, itypeRef, numRefType, itypeNeigh, &
    numNeighType, iobjFlag, limFlag, irefFlag, neighFlag)
  parameter (limType = 50, itypeAll = 999)
  integer*4 numRefType(*), numNeighType(*)         !# of types for ref and neigh
  integer*4 itypeRef(limType,*), itypeNeigh(limType,*), iobjFlag(*)
  logical feasible
  integer*4 icanDoPairs(2,10)
  data icanDoPairs/1, 1, 1, 2, 1, 4, 2, 1, 2, 2, 2, 4, 4, 1, 4, 2, 4, 4, -1, -1/
  do ii = 1, numGraphs
    if (numRefType(ii) == 1 .and. itypeRef(1, ii) == itypeAll) then
      do i = 1, limFlag
        if (iobjFlag(i) >= 0) then
          if (irefFlag < 0) then
            irefFlag = iobjFlag(i)
          elseif (iobjFlag(i) .ne. irefFlag) then
            print *,'Reference object types do not all match'
            go to 10
          endif
        endif
      enddo
    else
      do i = 1, numRefType(ii)
        itype = abs(itypeRef(i, ii))
        if (iobjFlag(itype) >= 0) then
          if (irefFlag < 0) then
            irefFlag = iobjFlag(itype)
          elseif (iobjFlag(itype) .ne. irefFlag) then
            print *,'Reference object types do not all match'
            go to 10
          endif
        else
          print *,'There are no contours in reference object', itype
          go to 10
        endif
      enddo
    endif
    if (numNeighType(ii) == 1 .and. itypeNeigh(1, ii) == itypeAll) then
      do i = 1, limFlag
        if (iobjFlag(i) >= 0) then
          if (neighFlag < 0) then
            neighFlag = iobjFlag(i)
          elseif (iobjFlag(i) .ne. neighFlag) then
            print *,'Neighboring object types do not all match'
            go to 10
          endif
        endif
      enddo
    else
      do i = 1, numNeighType(ii)
        itype = abs(itypeNeigh(i, ii))
        if (iobjFlag(itype) >= 0) then
          if (neighFlag < 0) then
            neighFlag = iobjFlag(itype)
          elseif (iobjFlag(itype) .ne. neighFlag) then
            print *,'Neighboring object types do not all match'
            go to 10
          endif
        else
          print *,'There are no contours in neighbor object', itype
          go to 10
        endif
      enddo
    endif
  enddo
  !
  ! check for feasibility
  !
  i = 1
  feasible = .false.
  do while(icanDoPairs(1, i) >= 0 .and. .not.feasible)
    feasible = icanDoPairs(1, i) == irefFlag .and. &
        icanDoPairs(2, i) == neighFlag
    i = i + 1
  enddo
  if (feasible) return
  print *,'Distances cannot be measured from the reference to the neighbor types'
10 irefFlag = -1
  neighFlag = -1
  return
end subroutine checkFlags



! INTEGRATE computes the integral of the items in some bins of a graph,
! subtracting a baseline value as well.
! GRAPHS has the values of the graph
! numBins is number of bins
! integrateStart, integrateEnd are starting and ending bins to integrate
! ibaseStart, ibaseEnd are starting and ending bins to compute the baseline
! value from, or 0, 0 to use the value supplied in baseVal
! SUM is returned with the integral
! CENTROID is returned with centroid of distance
subroutine integrate(graphs, areas, numBins, deltaRad, power, integrateStart, &
    integrateEnd, ibaseStart, ibaseEnd, baseVal, sum, centroid)
  implicit none
  real*4 graphs(*), areas(*), deltaRad, power, baseVal, sum, centroid
  integer*4 numBins, integrateStart, integrateEnd, ibaseStart, ibaseEnd
  real*4 baseline, distSum, counts, frac
  integer*4 ibin
  baseline = baseVal
  if (ibaseStart <= ibaseEnd .and. (ibaseStart .ne. 0 .or. ibaseEnd .ne. 0)) then
    baseline = 0.
    do ibin = ibaseStart, ibaseEnd
      baseline = baseline + graphs(ibin)
    enddo
    baseline = baseline / (ibaseEnd + 1 - ibaseStart)
  endif
  sum = 0.
  distSum = 0.
  centroid = 0.
  do ibin = integrateStart, integrateEnd
    ! if (power==0.) then
    ! radpow=1.
    ! else
    ! radpow=(2.*(ibin-0.5)*delr)**power
    ! endif
    ! frac=radpow*delr*3.14159
    frac = areas(ibin)
    counts = (graphs(ibin) - baseline) * frac
    sum = sum + counts
    distSum = distSum + deltaRad * (ibin - 0.5) * counts
  enddo
  if (sum > 0.) centroid = distSum / sum
  return
end subroutine integrate

! LINFIT fits a straight line to the N points in arrays X and Y by the
! method of least squares, returning SLOPE, and intercept BINT
!
subroutine linfit(x, y, n, slope, bIntercept)
  dimension x(*), y(*)
  if (n == 2) then
    slope = (y(2) - y(1)) / (x(2) - x(1))
    bIntercept = y(1) - slope * x(1)
  else
    sumX = 0.
    sumXsq = 0.
    sumY = 0.
    sumXY = 0.
    sumYsq = 0.
    do i = 1, n
      sumX = sumX + x(i)
      sumY = sumY + y(i)
      sumXY = sumXY + x(i) * y(i)
      sumXsq = sumXsq + x(i)**2
      sumYsq = sumYsq + y(i)**2
    enddo
    d = n * sumXsq - sumX**2
    slope = (n * sumXY - sumX * sumY) / d
    bIntercept = (sumXsq * sumY - sumX * sumXY) / d
  endif
  return
end subroutine linfit


! SEGMENT_DIST_Z computes the distance between segments of two
! "z-based lines", characterized by x=A1*z+B1 and y=C1*z+D1 for z
! between ZS1 and ZN1, and x=A2*z+B2 and y=C2*z+D2 for z between ZS2
! and ZN2.  It returns the sqaure of the distance in DIST, and the
! z values of the points of closest approach on the two segments
! in Z1 and Z2.
!
subroutine segment_dist_z(a1, b1, c1, d1, a2, b2, c2, d2, zStart1, zEnd1, zStart2, zEnd2 &
    , z1, z2, dist)
  logical out1, out2
  integer*4 ifTrace/0/, numTrace/0/
  point_to_line(a, b, c, d, x, y, z) = (a * x + c * y + z - a * b - c * d) /  &
      (1. +a**2 + c**2)
  !
  ! first find z1, z2 position of global minimum
  !
  sqr1 = 1 + a1**2 + c1**2
  sqr2 = 1 + a2**2 + c2**2
  crossTerm = -(1 + a1 * a2 + c1 * c2)
  const1 = -(a1 * (b1 - b2) + c1 * (d1 - d2))
  const2 = a2 * (b1 - b2) + c2 * (d1 - d2)
  den = sqr1 * sqr2 - crossTerm**2
  z1Numer = const1 * sqr2 - const2 * crossTerm
  z2num = sqr1 * const2 - const1 * crossTerm
  if (abs(den) < 1.e-20 .or. &
      abs(den) < 1.e-6 * max(abs(z1Numer), abs(z2num))) then
    !
    ! parallel lines: "just" check the 4 endpoints
    !
    xStart1 = a1 * zStart1 + b1
    xEnd1 = a1 * zEnd1 + b1
    yStart1 = c1 * zStart1 + d1
    yEnd1 = c1 * zEnd1 + d1
    xStart2 = a2 * zStart2 + b2
    xEnd2 = a2 * zEnd2 + b2
    yStart2 = c2 * zStart2 + d2
    yEnd2 = c2 * zEnd2 + d2
    z2Trunc = max(zStart2, min(zEnd2, point_to_line(a2, b2, c2, d2, xStart1, yStart1, &
        zStart1)))
    dist = (a2 * z2Trunc + b2 - xStart1)**2 + (c2 * z2Trunc + d2 - yStart1)**2 + &
        (z2Trunc - zStart1)**2
    z1 = zStart1
    z2 = z2Trunc
    z2Trunc = max(zStart2, min(zEnd2, point_to_line(a2, b2, c2, d2, xEnd1, yEnd1, zEnd1)))
    tdist = (a2 * z2Trunc + b2 - xEnd1)**2 + (c2 * z2Trunc + d2 - yEnd1)**2 + &
        (z2Trunc - zEnd1)**2
    if (tdist < dist) then
      dist = tdist
      z1 = zEnd1
      z2 = z2Trunc
    endif
    z1Trunc = max(zStart1, min(zEnd1, point_to_line(a1, b1, c1, d1, xStart2, yStart2, &
        zStart2)))
    tdist = (a1 * z1Trunc + b1 - xStart2)**2 + (c1 * z1Trunc + d1 - yStart2)**2 + &
        (z1Trunc - zStart2)**2
    if (tdist < dist) then
      dist = tdist
      z1 = z1Trunc
      z2 = zStart2
    endif
    z1Trunc = max(zStart1, min(zEnd1, point_to_line(a1, b1, c1, d1, xEnd2, yEnd2, zEnd2)))
    tdist = (a1 * z1Trunc + b1 - xEnd2)**2 + (c1 * z1Trunc + d1 - yEnd2)**2 + &
        (z1Trunc - zEnd2)**2
    if (tdist < dist) then
      dist = tdist
      z1 = z1Trunc
      z2 = zEnd2
    endif
    go to 20
  endif
  z1 = z1Numer / den
  z2 = z2num / den
  out1 = z1 < zStart1 .or. z1 > zEnd1
  out2 = z2 < zStart2 .or. z2 > zEnd2
  if (out1 .and. out2) then
    z1 = max(zStart1, min(zEnd1, z1))
    z2Trunc = max(zStart2, min(zEnd2, point_to_line(a2, b2, c2, d2, &
        a1 * z1 + b1, c1 * z1 + d1, z1)))
    z2 = max(zStart2, min(zEnd2, z2))
    z1Trunc = max(zStart1, min(zEnd1, point_to_line(a1, b1, c1, d1, &
        a2 * z2 + b2, c2 * z2 + d2, z2)))
    if (z1 .ne. z1Trunc .or. z2 .ne. z2Trunc) then
      dist = (a1 * z1 + b1 - a2 * z2Trunc - b2)**2 +  &
          (c1 * z1 + d1 - c2 * z2Trunc - d2)**2 + (z1 - z2Trunc)**2
      tdist = (a1 * z1Trunc + b1 - a2 * z2 - b2)**2 +  &
          (c1 * z1Trunc + d1 - c2 * z2 - d2)**2 + (z1Trunc - z2)**2
      if (dist < tdist) then
        z2 = z2Trunc
      else
        dist = tdist
        z1 = z1Trunc
      endif
      go to 20
    endif
  elseif (out1) then
    z1 = max(zStart1, min(zEnd1, z1))
    z2 = max(zStart2, min(zEnd2, point_to_line(a2, b2, c2, d2, &
        a1 * z1 + b1, c1 * z1 + d1, z1)))
  elseif (out2) then
    z2 = max(zStart2, min(zEnd2, z2))
    z1 = max(zStart1, min(zEnd1, point_to_line(a1, b1, c1, d1, &
        a2 * z2 + b2, c2 * z2 + d2, z2)))
  endif
  dist = (a1 * z1 + b1 - a2 * z2 - b2)**2 + (c1 * z1 + d1 - c2 * z2 - d2)**2 + &
      (z1 - z2)**2
20 if (ifTrace == 0) return
  numTrace = numTrace + 1
  if (numTrace < ifTrace) return
  numTrace = 0
  dmin = 1.e20
  do i = 0, 100
    z1Trunc = zStart1 + i * (zEnd1 - zStart1) / 100.
    x1trunc = a1 * z1Trunc + b1
    y1Trunc = c1 * z1Trunc + d1
    do j = 0, 100
      z2Trunc = zStart2 + j * (zEnd2 - zStart2) / 100.
      d = (a2 * z2Trunc + b2 - x1trunc)**2 + (c2 * z2Trunc + d2 - y1Trunc)**2 +  &
          (z2Trunc - z1Trunc)**2
      if (d < dmin) then
        dmin = d
        z1m = z1Trunc
        z2m = z2Trunc
      endif
    enddo
  enddo
  write(*,101) a1, b1, c1, d1, zStart1, zEnd1, a2, b2, c2, d2, zStart2, zEnd2, z1, z2, &
      dist, z1m, z2m, dmin
101 format(6f10.5,/,6f10.5,/,3f8.5,/,3f8.5)
  return
end subroutine segment_dist_z



subroutine get_next_sample &
    (xModPt, yModPt, zModPt, iref, sampleFrac, limRef, sampleLen, ifCloseSeg, &
    xTemp, yTemp, zTemp, numTemp, irefEnd, fracEnd, segmentLen, &
    xTempMin, xTempMax, yTempMin, yTempMax, zTempMin, zTempMax, &
    refXmin, refXmax, refYmin, refYmax, refZmin, refZmax)
  real*4 xModPt(*), yModPt(*), zModPt(*), xTemp(*), yTemp(*), zTemp(*)
  real*4 xTempMin(*), yTempMin(*), zTempMin(*)
  real*4 xTempMax(*), yTempMax(*), zTempMax(*)
  !
  segmentLen = 0.
  numTemp = 0
  if (iref >= limRef) return
  !
  ! start by putting out the interpolated first point
  !
  numTemp = 1
  next = iref + 1
  xTemp(1) = (1. -sampleFrac) * xModPt(iref) + sampleFrac * xModPt(iref + 1)
  yTemp(1) = (1. -sampleFrac) * yModPt(iref) + sampleFrac * yModPt(iref + 1)
  zTemp(1) = (1. -sampleFrac) * zModPt(iref) + sampleFrac * zModPt(iref + 1)
  !
  ! loop until end of line or sample length reached
  !
  do while(next <= limRef .and. segmentLen < sampleLen)
    !
    ! length of current segment: if it is too short, just skip it
    !
    curLen = sqrt((xModPt(next) - xTemp(numTemp))**2 + &
        (yModPt(next) - yTemp(numTemp))**2 + &
        (zModPt(next) - zTemp(numTemp))**2)
    if (curLen < 1.e-6 * sampleLen) then
      next = next + 1
    else if (segmentLen + curLen < sampleLen) then
      !
      ! if that still falls short, add point to output list and advance
      !
      numTemp = numTemp + 1
      xTemp(numTemp) = xModPt(next)
      yTemp(numTemp) = yModPt(next)
      zTemp(numTemp) = zModPt(next)
      segmentLen = segmentLen + curLen
      next = next + 1
    else
      !
      ! otherwise, find interpolated point, add that to list
      !
      irefEnd = next - 1
      fracEnd = (sampleLen - segmentLen) / curLen
      if (irefEnd == iref) fracEnd = sampleFrac + (1. -sampleFrac) * fracEnd
      segmentLen = sampleLen
      numTemp = numTemp + 1
      xTemp(numTemp) = (1. -fracEnd) * xModPt(irefEnd) + fracEnd * xModPt(irefEnd + 1)
      yTemp(numTemp) = (1. -fracEnd) * yModPt(irefEnd) + fracEnd * yModPt(irefEnd + 1)
      zTemp(numTemp) = (1. -fracEnd) * zModPt(irefEnd) + fracEnd * zModPt(irefEnd + 1)
    endif
  enddo
  !
  ! if termination was because at end of segment, set ending values
  !
  if (segmentLen < sampleLen) then
    irefEnd = limRef
    fracEnd = 0.
  endif
  !
  ! do min's and max's as necessary
  !
  if (ifCloseSeg == 0) then
    refXmin = xTemp(1)
    refYmin = yTemp(1)
    refZmin = zTemp(1)
    refXmax = xTemp(1)
    refYmax = yTemp(1)
    refZmax = zTemp(1)
    numTemp = 1
  else
    refXmin = 1.e10
    refYmin = 1.e10
    refZmin = 1.e10
    refXmax = -1.e10
    refYmax = -1.e10
    refZmax = -1.e10
    do i = 1, numTemp - 1
      xTempMin(i) = min(xTemp(i), xTemp(i + 1))
      yTempMin(i) = min(yTemp(i), yTemp(i + 1))
      zTempMin(i) = min(zTemp(i), zTemp(i + 1))
      xTempMax(i) = max(xTemp(i), xTemp(i + 1))
      yTempMax(i) = max(yTemp(i), yTemp(i + 1))
      zTempMax(i) = max(zTemp(i), zTemp(i + 1))
      refXmin = min(refXmin, xTempMin(i))
      refYmin = min(refYmin, yTempMin(i))
      refZmin = min(refZmin, zTempMin(i))
      refXmax = max(refXmax, xTempMax(i))
      refYmax = max(refYmax, yTempMax(i))
      refZmax = max(refZmax, zTempMax(i))
    enddo
  endif
  return
end subroutine get_next_sample


subroutine trim_win_seg(x1, y1, z1, x2, y2, z2, s1, s2)
  if (s1 == 0 .and. s2 == 0.) return
  untrim = sqrt((x2 - x1)**2 + (y2 - y1)**2 + (z2 - z1)**2)
  f1 = min(1., s1 / untrim)
  f2 = min(1., s2 / untrim)
  if (f1 + f2 >= 1.) then
    fs = f1 + f2
    f1 = f1 / fs
    f2 = f2 / fs
  endif
  x1trunc = x1 + f1 * (x2 - x1)
  x2 = x2 + f2 * (x1 - x2)
  x1 = x1trunc
  y1Trunc = y1 + f1 * (y2 - y1)
  y2 = y2 + f2 * (y1 - y2)
  y1 = y1Trunc
  z1Trunc = z1 + f1 * (z2 - z1)
  z2 = z2 + f2 * (z1 - z2)
  z1 = z1Trunc
  return
end subroutine trim_win_seg

subroutine save_connector(x1min, y1min, z1min, x2min, y2min, z2min, &
    xModPt, yModPt, zModPt, LIMXYZ, indFree, iobjRef, iobjNeigh, iobjWin, &
    numObjInWin, numInWin, numInWinTot, numRefInWin, numNeighInWin)
  real*4 xModPt(*), yModPt(*), zModPt(*)
  integer*4 iobjWin(*)
  !
  numInWinTot = numInWinTot + 1
  ifRefInWin = 0
  ifNeighInWin = 0
  do iow = 1, numObjInWin
    if (iobjRef == iobjWin(iow)) ifRefInWin = 1
    if (iobjNeigh == iobjWin(iow)) ifNeighInWin = 1
  enddo
  if (ifRefInWin == 0) then
    numObjInWin = numObjInWin + 1
    iobjWin(numObjInWin) = iobjRef
    numRefInWin = numRefInWin + 1
  endif
  if (ifNeighInWin == 0) then
    numObjInWin = numObjInWin + 1
    iobjWin(numObjInWin) = iobjNeigh
    numNeighInWin = numNeighInWin + 1
  endif
  if (indFree + 2 > LIMXYZ) return
  numInWin = numInWin + 1
  xModPt(indFree) = x1min
  xModPt(indFree + 1) = 0.5 * (x1min + x2min)
  xModPt(indFree + 2) = x2min
  yModPt(indFree) = y1min
  yModPt(indFree + 1) = 0.5 * (y1min + y2min)
  yModPt(indFree + 2) = y2min
  zModPt(indFree) = z1min
  zModPt(indFree + 1) = 0.5 * (z1min + z2min)
  zModPt(indFree + 2) = z2min
  indFree = indFree + 3
  return
end subroutine save_connector



! SEGMENT_DIST measures the distance between two line segments, one
! from XS1, YS1, ZS1 to XN1, YN1, ZN1 and one from XS2, YS2, ZS2 to
! XN2, YN2, ZN2.  It returns the square of the distance in DIST and two
! parameters, T1 and T2, specifying the fraction of the distance along
! the segments where the points of closest approach are.
!
subroutine segment_dist(xStart1, yStart1, zStart1, xEnd1, yEnd1, zEnd1 &
    , xStart2, yStart2, zStart2, xEnd2, yEnd2, zEnd2, t1, t2, dist)
  logical out1, out2
  integer*4 ifTrace/0/, numTrace/0/
  save numTrace
  !
  ! function returns the parameter t for the point along the line
  ! closest to x, y, z
  !
  point_to_line(xStart, yStart, zStart, a, b, c, x, y, z) = &
      (a * (x - xStart) + b * (y - yStart) + c * (z - zStart)) / (a**2 + b**2 + c**2)
  !
  ! first find z1, z2 position of global minimum
  !
  a1 = xEnd1 - xStart1
  b1 = yEnd1 - yStart1
  c1 = zEnd1 - zStart1
  a2 = xEnd2 - xStart2
  b2 = yEnd2 - yStart2
  c2 = zEnd2 - zStart2
  sqr1 = a1**2 + b1**2 + c1**2
  sqr2 = a2**2 + b2**2 + c2**2
  crossTerm = -(a1 * a2 + b1 * b2 + c1 * c2)
  const1 = a1 * (xStart2 - xStart1) + b1 * (yStart2 - yStart1) + c1 * (zStart2 - zStart1)
  const2 = -(a2 * (xStart2 - xStart1) + b2 * (yStart2 - yStart1) +  &
      c2 * (zStart2 - zStart1))
  den = sqr1 * sqr2 - crossTerm**2
  t1Numer = const1 * sqr2 - const2 * crossTerm
  t2Numer = sqr1 * const2 - const1 * crossTerm
  if (abs(den) < 1.e-20 .or. &
      abs(den) < 1.e-6 * max(abs(t1Numer), abs(t2Numer))) then
    !
    ! parallel lines: "just" check the 4 endpoints
    ! start of line 1 versus line 2
    !
    t2Trunc = max(0., min(1., point_to_line &
        (xStart2, yStart2, zStart2, a2, b2, c2, xStart1, yStart1, zStart1)))
    dist = (a2 * t2Trunc + xStart2 - xStart1)**2 + (b2 * t2Trunc + yStart2 - yStart1)**2 &
        + (c2 * t2Trunc + zStart2 - zStart1)**2
    t1 = 0.
    t2 = t2Trunc
    !
    ! end of line 1 versus line 2
    !
    t2Trunc = max(0., min(1., point_to_line &
        (xStart2, yStart2, zStart2, a2, b2, c2, xEnd1, yEnd1, zEnd1)))
    tdist = (a2 * t2Trunc + xStart2 - xEnd1)**2 + (b2 * t2Trunc + yStart2 - yEnd1)**2 + &
        (c2 * t2Trunc + zStart2 - zEnd1)**2
    if (tdist < dist) then
      dist = tdist
      t1 = 1.
      t2 = t2Trunc
    endif
    !
    ! start of line 2 versus line 1
    !
    t1Trun = max(0., min(1., point_to_line &
        (xStart1, yStart1, zStart1, a1, b1, c1, xStart2, yStart2, zStart2)))
    tdist = (a1 * t1Trun + xStart1 - xStart2)**2 + (b1 * t1Trun + yStart1 - yStart2)**2 &
        + (c1 * t1Trun + zStart1 - zStart2)**2
    if (tdist < dist) then
      dist = tdist
      t1 = t1Trun
      t2 = 0.
    endif
    !
    ! end of line 2 versus line 1
    !
    t1Trun = max(0., min(1., point_to_line &
        (xStart1, yStart1, zStart1, a1, b1, c1, xEnd2, yEnd2, zEnd2)))
    tdist = (a1 * t1Trun + xStart1 - xEnd2)**2 + (b1 * t1Trun + yStart1 - yEnd2)**2 + &
        (c1 * t1Trun + zStart1 - zEnd2)**2
    if (tdist < dist) then
      dist = tdist
      t1 = t1Trun
      t2 = 1.
    endif
    go to 20
  endif
  t1 = t1Numer / den
  t2 = t2Numer / den
  out1 = t1 < 0. .or. t1 > 1.
  out2 = t2 < 0. .or. t2 > 1.
  if (out1 .and. out2) then
    !
    ! if both closest points are out of bounds, truncate each one to
    ! its segment, then find closest point on other segment to that
    ! truncated point.  If this gives different answers, pick the
    ! pair with the closest approach
    !
    t1 = max(0., min(1., t1))
    t2Trunc = max(0., min(1., point_to_line(xStart2, yStart2, zStart2, a2, b2, c2, &
        a1 * t1 + xStart1, b1 * t1 + yStart1, c1 * t1 + zStart1)))
    t2 = max(0., min(1., t2))
    t1Trun = max(0., min(1., point_to_line(xStart1, yStart1, zStart1, a1, b1, c1, &
        a2 * t2 + xStart2, b2 * t2 + yStart2, c2 * t2 + zStart2)))
    if (t1 .ne. t1Trun .or. t2 .ne. t2Trunc) then
      dist = (a1 * t1 + xStart1 - a2 * t2Trunc - xStart2)**2 +  &
          (b1 * t1 + yStart1 - b2 * t2Trunc - yStart2)**2 + &
          (c1 * t1 + zStart1 - c2 * t2Trunc - zStart2)**2
      tdist = (a1 * t1Trun + xStart1 - a2 * t2 - xStart2)**2 +  &
          (b1 * t1Trun + yStart1 - b2 * t2 - yStart2)**2 + &
          (c1 * t1Trun + zStart1 - c2 * t2 - zStart2)**2
      if (dist < tdist) then
        t2 = t2Trunc
      else
        dist = tdist
        t1 = t1Trun
      endif
      go to 20
    endif
  elseif (out1) then
    !
    ! If outside one segment but in the other, truncate the one it is
    ! outside, then find closest point to other segment
    !
    t1 = max(0., min(1., t1))
    t2 = max(0., min(1., point_to_line(xStart2, yStart2, zStart2, a2, b2, c2, &
        a1 * t1 + xStart1, b1 * t1 + yStart1, c1 * t1 + zStart1)))
  elseif (out2) then
    t2 = max(0., min(1., t2))
    t1 = max(0., min(1., point_to_line(xStart1, yStart1, zStart1, a1, b1, c1, &
        a2 * t2 + xStart2, b2 * t2 + yStart2, c2 * t2 + zStart2)))
  endif
  dist = (a1 * t1 + xStart1 - a2 * t2 - xStart2)**2 + &
      (b1 * t1 + yStart1 - b2 * t2 - yStart2)**2 + &
      (c1 * t1 + zStart1 - c2 * t2 - zStart2)**2
20 if (ifTrace == 0) return
  numTrace = numTrace + 1
  if (numTrace < ifTrace) return
  numTrace = 0
  dmin = 1.e20
  do i = 0, 100
    t1Trun = i / 100.
    x1trunc = a1 * t1Trun + xStart1
    y1Trunc = b1 * t1Trun + yStart1
    z1Trunc = c1 * t1Trun + zStart1
    do j = 0, 100
      t2Trunc = j / 100.
      d = (a2 * t2Trunc + xStart2 - x1trunc)**2 + (b2 * t2Trunc + yStart2 - y1Trunc)**2 &
          + (c2 * t2Trunc + zStart2 - z1Trunc)**2
      if (d < dmin) then
        dmin = d
        t1m = t1Trun
        t2m = t2Trunc
      endif
    enddo
  enddo
  write(*,101) xStart1, yStart1, zStart1, xEnd1, yEnd1, zEnd1, &
      xStart2, yStart2, zStart2, xEnd2, yEnd2, zEnd2, t1, t2, &
      dist, t1m, t2m, dmin
101 format(6f10.5,/,6f10.5,/,3f8.5,/,3f8.5)
  return
end subroutine segment_dist


! INSIDE_TRIANGLE returns true if the point XT, YT is inside or on the
! boundary of the triangle whose corners are in arrays BX and BY.
! It is somewhat faster than INSIDE.  It is based on "Computational
! Geometry in C" by Joseph O'Rourke, 1998.
!
logical function inside_triangle(bx, by, xTest, yTest)
  real*4 bx(*), by(*)
  logical inside
  !
  ! compute signed area of each triangle between the point and an edge
  !
  area0 = (bx(1) - xTest) * (by(2) - yTest) - (bx(2) - xTest) * (by(1) - yTest)
  area1 = (bx(2) - xTest) * (by(3) - yTest) - (bx(3) - xTest) * (by(2) - yTest)
  area2 = (bx(3) - xTest) * (by(1) - yTest) - (bx(1) - xTest) * (by(3) - yTest)

  if (area0 .ne. 0. .and. area1 .ne. 0. .and. area2 .ne. 0.) then
    inside_triangle = area0 > 0. .and. area1 > 0. .and. area2 > 0. &
        .or. area0 < 0. .and. area1 < 0. .and. area2 < 0.
    return
  endif
  if (area0 == 0. .and. area1 == 0. .and. area2 == 0.) then
    inside_triangle = inside(bx, by, 3, xTest, yTest)
    return
  endif

  inside_triangle = area0 == 0. .and. area1 > 0. .and. area2 > 0. &
      .or. area0 == 0. .and. area1 < 0. .and. area2 < 0. &
      .or. area1 == 0. .and. area0 > 0. .and. area2 > 0. &
      .or. area1 == 0. .and. area0 < 0. .and. area2 < 0. &
      .or. area2 == 0. .and. area0 > 0. .and. area1 > 0. &
      .or. area2 == 0. .and. area0 < 0. .and. area1 < 0. &
      .or. area0 == 0. .and. area1 == 0. &
      .or. area0 == 0. .and. area2 == 0. &
      .or. area0 == 1. .and. area2 == 0.
  return
end function inside_triangle




! POINT_LINE_DIST measures the distance from the point X, Y, Z to the
! line segment from XS, YS, ZS ato XN, YN, ZN.  It returns the square of
! the distance in DSQR and the  parameter T specifying the position
! along the segment of the point of closest approach (between 0 at
! XS, YS, ZS and 1 at XN, YN, ZN) .
!
subroutine point_line_dist(xStart, yStart, zStart, xEnd, yEnd, zEnd, x, y, z, t, dsqr)
  a = xEnd - xStart
  b = yEnd - yStart
  c = zEnd - zStart
  t = 0
  if (a**2 + b**2 + c**2 > 1.e-24) &
      t = max(0., min(1., (a * (x - xStart) + b * (y - yStart) +  &
      c * (z - zStart)) / (a**2 + b**2 + c**2)))
  dsqr = (a * t + xStart - x)**2 + (b * t + yStart - y)**2 + (c * t + zStart - z)**2
  return
end subroutine point_line_dist



! POINT_TO_TRIANGLE measures the distance of closest approach between
! the point XP, YP, ZP and the triangle ITRI in common MTKCOM.  It
! returns the separation in DIST, and the coordinates in the rotated
! triangle, XROT and YROT, of the point of closest approach to the
! triangle.
!
subroutine point_to_triangle(xp, yp, zp, itriang, xRot, yRot, dist)
  use mtkvars
  logical inside_triangle
  !
  ! rotate the point
  !
  temp = xp * cosGamma(itriang) - yp * sinGamma(itriang)
  xr = temp * cosBeta(itriang) - zp * sinBeta(itriang)
  yr = xp * sinGamma(itriang) + yp * cosGamma(itriang)
  zr = temp * sinBeta(itriang) + zp * cosBeta(itriang)
  if (inside_triangle(triXYrot(1, 1, itriang), triXYrot(1, 2, itriang), xr, yr)) then
    !
    ! if x, y is inside planar triangle, that's the answer
    !
    xRot = xr
    yRot = yr
    dist = abs(zr - triZrot(itriang))
  else
    !
    ! otherwise, get distance to each edge of triangle
    !
    dmin = 1.e10
    do iv = 1, 3
      ivNext = iv + 1
      if (ivNext == 4) ivNext = 1
      call point_line_dist(triXYrot(iv, 1, itriang), triXYrot(iv, 2, itriang), &
          triZrot(itriang), triXYrot(ivNext, 1, itriang), triXYrot(ivNext, 2, itriang), &
          triZrot(itriang), xr, yr, zr, t, dsqr)
      if (dsqr < dmin) then
        ivMin = iv
        ivNextMin = ivNext
        dmin = dsqr
        tmin = t
      endif
    enddo
    dist = sqrt(dmin)
    xRot = (1. -tmin) * triXYrot(ivMin, 1, itriang) + tmin * triXYrot(ivNextMin, 1, &
        itriang)
    yRot = (1. -tmin) * triXYrot(ivMin, 2, itriang) + tmin * triXYrot(ivNextMin, 2, &
        itriang)
  endif
  return
end subroutine point_to_triangle



! SEGMENT_TO_TRIANGLE measures the distance of closest approach between
! the line segment connecting XS, YS, ZS to XN, YN, ZN and the triangle
! ITRI in common MTKCOM.  It returns the separation in DIST, the
! parameter T specifying the position along the segment of the point of
! closest approach (between 0 at XS, YS, ZS and 1 at XN, YN, ZN), and the
! coordinates in the rotated triangle, XROT and YROT of the point of
! closest approach to the triangle.
!
subroutine segment_to_triangle(xStart, yStart, zStart, xEnd, yEnd, zEnd, itriang, t, &
    xRot, yRot, dist)
  use mtkvars
  logical inside_triangle, startIn, endIn, b3dxor
  !
  ! rotate the endpoints
  !
  temp = xStart * cosGamma(itriang) - yStart * sinGamma(itriang)
  xStartRot = temp * cosBeta(itriang) - zStart * sinBeta(itriang)
  yStartRot = xStart * sinGamma(itriang) + yStart * cosGamma(itriang)
  zStartRot = temp * sinBeta(itriang) + zStart * cosBeta(itriang)
  temp = xEnd * cosGamma(itriang) - yEnd * sinGamma(itriang)
  xEndRot = temp * cosBeta(itriang) - zEnd * sinBeta(itriang)
  yEndRot = xEnd * sinGamma(itriang) + yEnd * cosGamma(itriang)
  zEndRot = temp * sinBeta(itriang) + zEnd * cosBeta(itriang)
  startIn = inside_triangle(triXYrot(1, 1, itriang), triXYrot(1, 2, itriang), xStartRot, &
      yStartRot)
  endIn = inside_triangle(triXYrot(1, 1, itriang), triXYrot(1, 2, itriang), xEndRot, &
      yEndRot)
  !
  if (startIn .and. endIn) then
    !
    ! if both endpoints are over triangle, then one must be closest
    ! unless line passes through triangle
    !
    if (b3dxor(zStartRot > triZrot(itriang), zEndRot > triZrot(itriang))) then
      dist = 0.
      t = (triZrot(itriang) - zStartRot) / (zEndRot - zStartRot)
      xRot = (1. -t) * xStartRot + t * xEndRot
      yRot = (1. -t) * yStartRot + t * yEndRot
    else
      distEnd = abs(zStartRot - triZrot(itriang))
      distStart = abs(zEndRot - triZrot(itriang))
      if (distEnd < distStart) then
        dist = distEnd
        xRot = xStartRot
        yRot = yStartRot
        t = 0.
      else
        dist = distStart
        xRot = xEndRot
        yRot = yEndRot
        t = 1.
      endif
    endif
    return
  endif
  !
  ! if one endpoint is over, it is a candidate
  !
  dist = 1.e10
  if (startIn) then
    xRot = xStartRot
    yRot = yStartRot
    dist = abs(zStartRot - triZrot(itriang))
    t = 0.
  endif
  if (endIn) then
    xRot = xEndRot
    yRot = yEndRot
    dist = abs(zEndRot - triZrot(itriang))
    t = 1.
  endif
  !
  ! but still need to check each line segment
  !
  dmin = dist**2
  ivMin = 0
  do iv = 1, 3
    ivNext = iv + 1
    if (ivNext > 3) ivNext = 1
    call segment_dist(triXYrot(iv, 1, itriang), triXYrot(iv, 2, itriang), &
        triZrot(itriang), triXYrot(ivNext, 1, itriang), triXYrot(ivNext, 2, itriang), &
        triZrot(itriang), xStartRot, yStartRot, zStartRot, xEndRot, yEndRot, zEndRot,  &
        t1, t2, dsqr)
    if (dsqr < dmin) then
      ivMin = iv
      ivNextMin = ivNext
      dmin = dsqr
      t1AtMin = t1
      t2AtMin = t2
    endif
  enddo
  !
  ! if a segment was it, need square root and rotated position
  !
  if (ivMin > 0) then
    dist = sqrt(dmin)
    xRot = (1. -t1AtMin) * triXYrot(ivMin, 1, itriang) + t1AtMin *  &
        triXYrot(ivNextMin, 1, itriang)
    yRot = (1. -t1AtMin) * triXYrot(ivMin, 2, itriang) + t1AtMin *  &
        triXYrot(ivNextMin, 2, itriang)
    t = t2AtMin
  endif
  return
end subroutine segment_to_triangle



! TRIANGLE_TO_TRIANGLE computes the minimum distance between two
! triangles.  It uses VERTS, INDVERT, XYROT, ZROT, SBET, CBET, SGAM,
! and CGAM from the common MTKCOM.  ITRI1 and ITRI2 are the indices of
! the two triangles.  The routine returns DIST, the distance between
! triangles; ITRIR, the index of the triangle in whose rotated
! coordinate system the closest approaching points are located;
! XR1, YR1, ZR1, the coordinates of the point on the other triangle;
! and XR2 and YR2, the coordinates of the closest point on the
! triangle ITRIR.
!
subroutine triangle_to_triangle(itri1, itri2, xRot1, yRot1, zRot1, xRot2, yRot2, &
    indTriRot, dist)
  use mtkvars
  logical inside_triangle, vertinTri(3, 2), b3dxor
  real*4 xr(3,2), yr(3,2), zr(3,2)
  integer*4 jtri(2)
  jtri(1) = itri1
  jtri(2) = itri2
  distMin = 1.e10
  do it = 1, 2
    !
    ! rotate the endpoints
    !
    do iv = 1, 3
      ind = indVert(iv, jtri(it))
      itriang = jtri(3 - it)
      temp = verts(1, ind) * cosGamma(itriang) - verts(2, ind) * sinGamma(itriang)
      xr(iv, it) = temp * cosBeta(itriang) - verts(3, ind) * sinBeta(itriang)
      yr(iv, it) = verts(1, ind) * sinGamma(itriang) + verts(2, ind) * cosGamma(itriang)
      zr(iv, it) = temp * sinBeta(itriang) + verts(3, ind) * cosBeta(itriang)
      vertinTri(iv, it) = inside_triangle(triXYrot(1, 1, itriang),  &
          triXYrot(1, 2, itriang), xr(iv, it), yr(iv, it))
    enddo
    !
    ! test for paired vertices on opposite sides of the other triangle
    !
    do iv = 1, 3
      ivNext = mod(iv, 3) + 1
      if (vertinTri(iv, it) .and. vertinTri(ivNext, it) .and. &
          b3dxor(zr(iv, it) > triZrot(itriang), zr(ivNext, it) > triZrot(itriang))) then
        dist = 0.
        t = (triZrot(itriang) - zr(iv, it)) / (zr(iv, it) - zr(ivNext, it))
        xRot1 = xr(iv, it) * (1. -t) + xr(ivNext, it) * t
        yRot1 = yr(iv, it) * (1. -t) + yr(ivNext, it) * t
        zRot1 = triZrot(itriang)
        xRot2 = xRot1
        yRot2 = yRot1
        indTriRot = itriang
        return
      endif
      !
      ! if any one endpoint is over, it is a candidate for minimum point
      ! being on the face of the triangle
      !
      dist = abs(zr(iv, it) - triZrot(itriang))
      if (vertinTri(iv, it) .and. dist < distMin) then
        distMin = dist
        xRot1 = xr(iv, it)
        yRot1 = yr(iv, it)
        zRot1 = zr(iv, it)
        xRot2 = xRot1
        yRot2 = yRot1
        indTriRot = itriang
      endif

    enddo
  enddo
  !
  ! but still need to check each line segment pair, in coordinates of
  ! rotated second triangle
  !
  dmin = distMin**2
  dist = distMin
  ivMin = 0
  itriang = itri2
  do iv = 1, 3
    ivNext = mod(iv, 3) + 1
    do jv = 1, 3
      jvNext = mod(jv, 3) + 1
      call segment_dist(triXYrot(iv, 1, itriang), triXYrot(iv, 2, itriang), &
          triZrot(itriang), triXYrot(ivNext, 1, itriang), triXYrot(ivNext, 2, itriang), &
          triZrot(itriang), xr(jv, 1), yr(jv, 1), zr(jv, 1), xr(jvNext, 1), &
          yr(jvNext, 1), zr(jvNext, 1), t1, t2, dsqr)
      if (dsqr < dmin) then
        ivMin = iv
        ivNextMin = ivNext
        jvMin = jv
        jvNextMin = jvNext
        dmin = dsqr
        t1AtMin = t1
        t2AtMin = t2
      endif
    enddo
  enddo
  !
  ! if a segment was it, need square root and interpolated positions
  ! relative to rotated second triangle
  !
  if (ivMin > 0) then
    dist = sqrt(dmin)
    indTriRot = itriang
    xRot1 = (1. -t2AtMin) * xr(jvMin, 1) + t2AtMin * xr(jvNextMin, 1)
    yRot1 = (1. -t2AtMin) * yr(jvMin, 1) + t2AtMin * yr(jvNextMin, 1)
    zRot1 = (1. -t2AtMin) * zr(jvMin, 1) + t2AtMin * zr(jvNextMin, 1)
    xRot2 = (1. -t1AtMin) * triXYrot(ivMin, 1, itriang) + &
        t1AtMin * triXYrot(ivNextMin, 1, itriang)
    yRot2 = (1. -t1AtMin) * triXYrot(ivMin, 2, itriang) + &
        t1AtMin * triXYrot(ivNextMin, 2, itriang)
  endif
  return
end subroutine triangle_to_triangle


subroutine crossproduct(avec, bvec, cvec)
  real*4 avec(*), bvec(*), cvec(*)
  cvec(1) = avec(2) * bvec(3) - avec(3) * bvec(2)
  cvec(2) = avec(3) * bvec(1) - avec(1) * bvec(3)
  cvec(3) = avec(1) * bvec(2) - avec(2) * bvec(1)
  return
end subroutine crossproduct


function indexshift(iobj, iflag, itypeObj, numPtInObj, numWobj)
  use mtkvars
  integer*4 itypeObj(*), numPtInObj(*)
  if (iflag < 4) then
    itype = itypeObj(iobj)
  else
    itype = iobjMesh(iobj)
  endif
  inList = 0
  do i = 1, numObjShifted
    if (iobjShift(i) == itype) inList = i
  enddo
  if (inList == 0) then
    indexshift = -1
    return
  endif
  indexshift = indStartShift(inList)
  if (iflag == 4) return
  do i = 1, iobj - 1
    if (itypeObj(i) == itype) then
      if (iflag == 1) then
        indexshift = indexshift + 1
      else
        indexshift = indexshift + numPtInObj(i)
      endif
    endif
  enddo
  return
end function indexshift
