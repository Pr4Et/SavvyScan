! GRAPHDEN produces a series of "graphs" of average density (per
! unit area) of neighboring items as a function of radial distance
! from an average reference point.  The reference points and the
! neighboring points may be of various types or sets of types, as
! defined separately for each graph.  The subroutine accurately
! corrects for edge effects, so that the points may be contained
! within an arbitrary, irregular boundary and the routine will compute
! density based only on the area within the boundary.
! NUMVERTS is the number of boundary points
! XBOUND, YBOUND are the coordinates of the boundary points
! NUMPOINTS is the number of items
! XPOINT, YPOINT, ITYPE are the coordinates and types of the items
! DELRAD is the bin width (in distance units) desired for the graphs
! NUMBINS is the number of bins for each graph
! NUMGRAPH is the number of graphs
! NUMREFTYPES is the number of types for reference points for each graph
! NUMNEIGHTYPES is the number of types for neighbor points for each graph
! ITYPEREF(I, J) is the Ith reference type for the Jth graph
! ITYPENEIGH(I, J) is the Ith neighbor type for the Jth graph
! GRAPHS(I, J) is returned with the density at the Ith bin of Jth graph
! FRACSUM(I, J) is returned with total area contributing to that bin
!
! $Id$

subroutine denGraph(xbound, ybound, numVerts, xpoint, ypoint, itype, numPoints, delRad, &
    numBins, numGraphs, numRefTypes, numNeighTypes, itypeRef, itypeNeigh, graphs,  &
    fracSum, iunitAllDist)
  integer ITYPEALL, LIMBINS, LIMCROSS, LIMGRAPHS, LIMPNTS, LIMPTS_X_GRFS, LIMTYPE, LIMVERT
  parameter (LIMGRAPHS = 50, LIMBINS = 1001, LIMPNTS = 50000, LIMCROSS = 20, &
      LIMVERT = 500, LIMTYPE = 50, ITYPEALL = 999, LIMPTS_X_GRFS = 100000)
  real*4 xbound(*), ybound(*)                        !boundary vertices
  real*4 xpoint(*), ypoint(*)                        !sample points
  real*4 graphs(LIMBINS,LIMGRAPHS)
  integer*4 itype(*)                        !types of sample points
  integer*4 numRefTypes(*), numNeighTypes(*)         !# of types for ref and neigh
  integer*4 itypeRef(LIMTYPE,*), itypeNeigh(LIMTYPE,*)
  !
  real*4 fracSum(LIMBINS,LIMGRAPHS), crossAngle(LIMCROSS,LIMBINS)
  real*4 edgeSq(LIMVERT)
  integer*4 numCross(LIMBINS), igraphRef(LIMGRAPHS), iunitAllDist
  logical*1 isNeighInGrf(LIMPTS_X_GRFS)
  real*4 annulusArea, ctmp, delRad, delX, delY, diffSum, distMaxSq
  real*4 distMin, fracadd, rr, radMax, rmaxSq, radMidMaxSq
  real*4 t1Sum, t2Sum, totalArea, x0, x1, x2, xsolve, y0, y1, y2, ysolve
  integer*4 i, ibinEnd, ibin, ibinMin, ibinStart, idir, ii, ineed, indRef
  integer*4 ivert, jj, kk, indNay, numBins, needed, neededAsRef, numGraphs
  integer*4 numPoints, numSolve, numVerts, ifOut
  real * 8 fsol, rad, rsq, distSq, d1Sq, d2Sq, dbSq, bOver2a, radPart, distMaxSqr, root
  logical inside
  real*4 atan2d
  !
  radMax = numBins * delRad
  rmaxSq = radMax**2
  radMidMaxSq = ((numBins - 0.5) * delRad)**2
  if (iunitAllDist > 0) write(iunitAllDist, '(a, g16.6)')'All graphed distances up to ', &
      radMax
  !
  ! Analyze which points are valid neighbor points for each graph
  !
  do ii = 1, numPoints
    do jj = 1, numGraphs
      isNeighInGrf((ii - 1) * numGraphs + jj) = .false.
      do kk = 1, numNeighTypes(jj)
        if (itypeNeigh(kk, jj) == ITYPEALL .or. &
            itypeNeigh(kk, jj) == itype(ii)) &
            isNeighInGrf((ii - 1) * numGraphs + jj) = .true.
      enddo
    enddo
  enddo
  !
  ! compute squares of edge lengths
  !
  do i = 1, numVerts
    edgeSq(i) = (xbound(i + 1) - xbound(i))**2 + (ybound(i + 1) - ybound(i))**2
  enddo

  !
  ! zero out the graphs and the fractional area tables
  !
  do ii = 1, numGraphs
    do jj = 1, numBins
      graphs(jj, ii) = 0.
      fracSum(jj, ii) = 0.
    enddo
  enddo
  !
  ! loop through each sample point considered as reference
  !
  t1Sum = 0.
  t2Sum = 0.
  do indRef = 1, numPoints
    !
    ! first make list of graphs the reference point is needed in
    !
    neededAsRef = 0
    do jj = 1, numGraphs
      needed = 0
      do kk = 1, numRefTypes(jj)
        if (itypeRef(kk, jj) == ITYPEALL .or. &
            itypeRef(kk, jj) == itype(indRef)) needed = 1
      enddo
      if (needed > 0) then
        neededAsRef = neededAsRef + 1
        igraphRef(neededAsRef) = jj
      endif
    enddo
    !
    if (neededAsRef > 0) then
      !
      ! t1=secnds(0.)
      x0 = xpoint(indRef)
      y0 = ypoint(indRef)
      !
      ! otherwise need to analyze number of crossings of annulus at
      ! each bin;  start by zeroing list of crossings
      !
      do ii = 1, numBins
        numCross(ii) = 0
      enddo
      !
      ! consider each edge of boundary in turn, find range of annuli
      ! that cross it
      !
      x1 = xbound(1)
      y1 = ybound(1)
      d1Sq = (x0 - x1)**2 + (y0 - y1)**2
      do ivert = 1, numVerts
        !
        ! get distance to nearest point within the edge segment
        !
        dbSq = edgeSq(ivert)
        x2 = xbound(ivert + 1)
        y2 = ybound(ivert + 1)
        d2Sq = (x0 - x2)**2 + (y0 - y2)**2
        if (d2Sq > d1Sq + dbSq) then
          distSq = d1Sq                     !if (x1, y1) is closest
        elseif (d1Sq > d2Sq + dbSq) then
          distSq = d2Sq                     !if (x2, y2) is closest
        else                              !if somewhere between is closer
          distSq = d2Sq - (d2Sq + dbSq - d1Sq)**2 / (4. *dbSq)
        endif
        !
        ! if that distance is smaller than radius of highest bin,
        ! get # of first bin that might cross edge segment
        !
        if (distSq < radMidMaxSq) then
          distMin = sqrt(max(0.d0, distSq))
          ibinMin = distMin / delRad + 1.5
          bOver2a = -(d2Sq - d1Sq - dbSq) / (2. *dbSq)
          radPart = bOver2a**2 - d1Sq / dbSq
          distMaxSq = max(d1Sq, d2Sq)
          !
          ! solve for crossings of each bin and store angles
          !
          do ibin = ibinMin, numBins
            rr = (ibin - 0.5) * delRad
            rsq = rr**2
            if (rsq > distMaxSq) go to 20
            rad = radPart + rsq / dbSq
            if (rad >= 0) then
              root = sqrt(max(0.d0, rad))
              numSolve = 0
              do idir = -1, 1, 2
                fsol = bOver2a + idir * root
                !
                ! make test open on one end to prevent double-counting
                ! of segment endpoint
                !
                if (fsol > 0 .and. fsol <= 1) then
                  numSolve = numSolve + 1
                  xsolve = fsol * (x2 - x1) + x1 - x0
                  ysolve = fsol * (y2 - y1) + y1 - y0
                  numCross(ibin) = numCross(ibin) + 1
                  crossAngle(numCross(ibin), ibin) = atan2d(ysolve, xsolve)
                endif
              enddo
              if (numSolve == 0) print *,'neither solution valid'
            else
              print *,'imaginary solution'
            endif
          enddo
20        continue
        endif
        x1 = x2
        y1 = y2
        d1Sq = d2Sq
      enddo
      !
      ! find first and last bins with non-zero #'s of crossings
      !
      ibinStart = numBins + 1
      ibinEnd = numBins
      do i = 1, numBins
        if (numCross(i) .ne. 0) then
          ibinEnd = i
          if (ibinStart == numBins + 1) ibinStart = i
        endif
      enddo
      !
      ! add 1 to fracsum for all bins up to the first, then analyze
      ! crossings for the rest
      !
      fracadd = 1.
      do ibin = 1, ibinEnd
        if (ibin >= ibinStart) then
          if (numCross(ibin) == 0) then
            print *,'no boundary crossings for bin between bins' &
                , ' with some crossings'
          elseif (mod(numCross(ibin), 2) .ne. 0) then
            print *,'odd number of crossings'
          else
            !
            ! order the crossing angles then add up difference between
            ! non-overlapping pairs
            !
            do ii = 1, numCross(ibin) - 1
              do jj = ii + 1, numCross(ibin)
                if (crossAngle(ii, ibin) > crossAngle(jj, ibin)) then
                  ctmp = crossAngle(ii, ibin)
                  crossAngle(ii, ibin) = crossAngle(jj, ibin)
                  crossAngle(jj, ibin) = ctmp
                endif
              enddo
            enddo
            !
            diffSum = 0.
            do ii = 2, numCross(ibin), 2
              diffSum = diffSum + crossAngle(ii, ibin) - &
                  crossAngle(ii - 1, ibin)
            enddo
            !
            ! if point at +/-180 degrees is inside, need the complement
            !
            if (inside(xbound, ybound, numVerts, x0 - (ibin - 0.5) * delRad, y0)) &
                diffSum = 360. -diffSum
            fracadd = diffSum / 360.
          endif
        endif
        !
        do ineed = 1, neededAsRef
          jj = igraphRef(ineed)
          fracSum(ibin, jj) = fracSum(ibin, jj) + fracadd
        enddo
      enddo
      ! t1sum=t1sum+secnds(t1)
      ! t2=secnds(0.)
      !
      ! now add neighboring points to bins as needed
      !
      do indNay = 1, numPoints
        if (indNay .ne. indRef) then
          delX = abs(x0 - xpoint(indNay))
          if (delX < radMax) then
            delY = abs(y0 - ypoint(indNay))
            if (delY < radMax) then
              rsq = delX**2 + delY**2
              if (rsq < rmaxSq) then
                ibin = sqrt(rsq) / delRad + 1
                ifOut = iunitAllDist
                do ineed = 1, neededAsRef
                  jj = igraphRef(ineed)
                  if (isNeighInGrf((indNay - 1) * numGraphs + jj)) then
                    graphs(ibin, jj) = graphs(ibin, jj) + 1.
                    if (ifOut > 0) then
                      write(iunitAllDist, '(2i10, g16.6)') indRef, indNay, sqrt(rsq)
                      ifOut = 0
                    endif
                  endif
                enddo
              endif
            endif
          endif
        endif
      enddo
      ! t2sum=t2sum+secnds(t2)
      !
    endif
  enddo
  ! print *,t1sum, t2sum
  !
  ! convert counts to densities
  !
!!$    do jj=1, ngraph
!!$    print *,'Fracsum #', jj
!!$    write(*,'(7f11.4)') (fracsum(ii, jj), ii=1, nbins)
!!$    print *,'Graph #', jj
!!$    write(*,'(7f11.4)') (graphs(ii, jj), ii=1, nbins)
!!$    if (jj.ne.ngraph) call mypause (' ')
!!$    enddo
  do ibin = 1, numBins
    annulusArea = 3.14159 * delRad**2 * (2 * ibin - 1)
    do jj = 1, numGraphs
      totalArea = annulusArea * fracSum(ibin, jj)
      if (totalArea == 0.) then
        graphs(ibin, jj) = 0.
      else
        graphs(ibin, jj) = graphs(ibin, jj) / totalArea
      endif
      fracSum(ibin, jj) = totalArea
    enddo
  enddo
  return
end subroutine denGraph
