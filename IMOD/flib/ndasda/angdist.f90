! ANGDIST computes a set of "graphs" of the density of items as a
! function of angular separation, for items within a specified annulus
! centered on a reference item.  The reference points, the
! neighboring points and their angular neighbors may be of various
! types or sets of types, as defined separately for each graph.  The
! subroutine  corrects approximately for edge effects, so that the
! points may be contained within an arbitrary, irregular boundary and
! the routine will compute density based only on the area within the
! boundary.
! NVERT is the number of boundary points
! XBOUND, YBOUND are the coordinates of the boundary points
! NPOINT is the number of items
! XPOINT, YPOINT, ITYPE are the coordinates and types of the items
! RADMIN, RADMAX are the inner and outer radii of the annulus
! NBINS is the number of bins for each graph
! NGRAPH is the number of graphs
! NREFTYP is the number of types for reference points for each graph
! NNEIGHTYP is the number of types for neighbor points for each graph
! NANGTYP is the number of types for angular neighbors to neighbor
! points for each graph
! ITYPREF(I, J) is the Ith reference type for the Jth graph
! ITYPNEIGH(I, J) is the Ith neighbor type for the Jth graph
! ITYPANG(I, J) is the Ith angluar neighbor type for the Jth graph
! GRAPHS(I, J) is returned with the density at the Ith bin of Jth graph
! FRACSUM(I, J) is returned with total area contributing to that bin
!
! $Id$
!
subroutine angDist(xbound, ybound, numVerts, xpoint, ypoint, itype, numPoints, radMin, &
    radMax, numBins, numGraphs, numRefTypes, numNeighTypes, numAngTypes, itypeRef, &
    itypeNeigh, itypeAngle, graphs, fracSum, iunitAllDist)
  integer ITYPEALL, LIMBINS, LIMCROSS, LIMGRAPHS, LIMNAY, LIMPNTS, LIMPTXGRF, LIMTYP
  integer LIMVERT
  parameter (LIMGRAPHS = 50, LIMBINS = 1001, LIMPNTS = 50000, LIMCROSS = 20, &
      LIMVERT = 500, LIMTYP = 50, ITYPEALL = 999, LIMNAY = 100, &
      LIMPTXGRF = 100000)
  real*4 xbound(*), ybound(*)                        !boundary vertices
  real*4 xpoint(*), ypoint(*)                        !sample points
  real*4 graphs(LIMBINS,LIMGRAPHS)
  integer*4 itype(*)                        !types of sample points
  integer*4 numRefTypes(*), numNeighTypes(*)         !# of types for ref and neigh
  integer*4 itypeRef(LIMTYP,*), itypeNeigh(LIMTYP,*)
  integer*4 numAngTypes(*), itypeAngle(LIMTYP,*)
  !
  real*4 fracSum(LIMBINS,LIMGRAPHS), crossAngle(LIMCROSS)
  real*4 edgeSq(LIMVERT), cwCcwMin(2), thetaNeigh(LIMNAY)
  integer*4 igraphRef(LIMGRAPHS), neighbor(LIMNAY), igraphNeigh(LIMNAY), iunitAllDist
  LOGICAL*1 isNeighInGrf(LIMPTXGRF), isAngPtInGraph(LIMPTXGRF)
  real*4 abs, angleDiff, binArea, bOver2a, ccwAngle, cwAngle, d1Sq
  real*4 d2Sq, dbSq, delX, delY, distSq, fsol, rad, radMax, radMaxSq, radMidMaxSq
  real*4 radMin, radMinSq, root, rsq, t1Sum
  real*4 t2Sum, totalArea, x0, x1, x2, xsolve, y0, y1, y2, ysolve
  integer*4 i, ibin, ic, icw, idir, ifix, ii, iptNum, indRef, ivert, ifOut
  integer*4 jgrf, jj, kk, max, indAng, indNay, numBinExtra, numBinAdd, numBins, numCross
  integer*4 needed, neededAsNeigh, neededAsRef, numGraphs, numInAnnulus
  integer*4 numPoints, numSolve, numVerts
  real*4 atan2d
  !
  radMaxSq = radMax**2
  radMinSq = radMin**2
  radMidMaxSq = (0.5 * (radMin + radMax))**2
  if (iunitAllDist > 0) write(iunitAllDist, '(a, g16.6, a, g16.6)') &
      'All graphed angles between distances of', radMin,'  and ', radMax
  !
  ! Analyze which points are valid neighbor and angle points for each
  ! graph
  !
  do ii = 1, numPoints
    do jj = 1, numGraphs
      isNeighInGrf((ii - 1) * numGraphs + jj) = .false.
      do kk = 1, numNeighTypes(jj)
        if (itypeNeigh(kk, jj) == ITYPEALL .or. &
            itypeNeigh(kk, jj) == itype(ii)) &
            isNeighInGrf((ii - 1) * numGraphs + jj) = .true.
      enddo
      isAngPtInGraph((ii - 1) * numGraphs + jj) = .false.
      do kk = 1, numAngTypes(jj)
        if (itypeAngle(kk, jj) == ITYPEALL .or. &
            itypeAngle(kk, jj) == itype(ii)) &
            isAngPtInGraph((ii - 1) * numGraphs + jj) = .true.
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
      numCross = 0
      !
      ! need to analyze number of boundary crossings by annulus
      ! consider each edge of boundary in turn, find crossings
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
        ! if that distance is smaller than middle radius of annulus,
        ! and max distance is higher, solve for crossing angles
        !
        if (distSq < radMidMaxSq .and. max(d1Sq, d2Sq) >= radMidMaxSq) &
            then
          bOver2a = -(d2Sq - d1Sq - dbSq) / (2. *dbSq)
          rad = bOver2a**2 + (radMidMaxSq - d1Sq) / dbSq
          if (rad >= 0) then
            root = sqrt(rad)
            numSolve = 0
            do idir = -1, 1, 2
              fsol = bOver2a + idir * root
              if (fsol >= 0 .and. fsol <= 1) then
                numSolve = numSolve + 1
                xsolve = fsol * (x2 - x1) + x1 - x0
                ysolve = fsol * (y2 - y1) + y1 - y0
                numCross = numCross + 1
                crossAngle(numCross) = atan2d(ysolve, xsolve)
              endif
            enddo
            if (numSolve == 0) print *,'neither solution valid'
          else
            print *,'imaginary solution'
          endif
        endif
        x1 = x2
        y1 = y2
        d1Sq = d2Sq
      enddo
      !
      if (mod(numCross, 2) .ne. 0) print *,'odd number of crossings'
      ! t1sum=t1sum+secnds(t1)
      ! t2=secnds(0.)
      !
      ! now make list of all neighboring points in annulus
      !
      numInAnnulus = 0
      do indNay = 1, numPoints
        if (indNay .ne. indRef) then
          delX = x0 - xpoint(indNay)
          if (abs(delX) < radMax) then
            delY = y0 - ypoint(indNay)
            if (abs(delY) < radMax) then
              rsq = delX**2 + delY**2
              if (rsq == 0.) then
                print *,'Warning: duplicate point - run CHECKMTMOD'
              elseif (rsq < radMaxSq .and. rsq >= radMinSq) then
                numInAnnulus = numInAnnulus + 1
                neighbor(numInAnnulus) = indNay
                thetaNeigh(numInAnnulus) = atan2d(delY, delX)
              endif
            endif
          endif
        endif
      enddo
      !
      ! loop through each neighbor point and consider it as reference
      !
      do indNay = 1, numInAnnulus
        !
        ! first make list of graphs this point is a valid neighbor for
        !
        neededAsNeigh = 0
        do jgrf = 1, neededAsRef
          if (isNeighInGrf(igraphRef(jgrf) + (neighbor(indNay) - 1) * numGraphs)) &
              then
            neededAsNeigh = neededAsNeigh + 1
            igraphNeigh(neededAsNeigh) = igraphRef(jgrf)
          endif
        enddo
        !
        if (neededAsNeigh > 0) then
          !
          ! if it's needed in any graphs, and there are no border
          ! crossings, just increment every bin of fracsum twice
          !
          if (numCross == 0) then
            do jgrf = 1, neededAsNeigh
              jj = igraphNeigh(jgrf)
              do ibin = 1, numBins
                fracSum(ibin, jj) = fracSum(ibin, jj) + 2.
              enddo
            enddo
          else
            !
            ! otherwise, need to determine which bins are inside border
            ! find minimum clockwise and countercw distance to a crossing
            !
            cwCcwMin(1) = 400.
            cwCcwMin(2) = 400.
            do ic = 1, numCross
              cwAngle = crossAngle(ic) - thetaNeigh(indNay)
              ccwAngle = -cwAngle
              if (cwAngle < 0.) cwAngle = cwAngle + 360.
              if (ccwAngle < 0.) ccwAngle = ccwAngle + 360.
              cwCcwMin(1) = min(cwCcwMin(1), cwAngle)
              cwCcwMin(2) = min(cwCcwMin(2), ccwAngle)
            enddo
            !
            ! add 1 to bins within boundary on each side
            !
            do icw = 1, 2
              numBinAdd = numBins * cwCcwMin(icw) / 180.
              numBinExtra = 2 * numBins + 1 - numBinAdd
              do jgrf = 1, neededAsNeigh
                jj = igraphNeigh(jgrf)
                do ibin = 1, min(numBins, numBinAdd)
                  fracSum(ibin, jj) = fracSum(ibin, jj) + 1.
                enddo
                do ibin = numBinExtra, numBins
                  fracSum(ibin, jj) = fracSum(ibin, jj) + 1.
                enddo
              enddo
            enddo
          endif
          !
          ! now process other neighbors in annulus
          !

          do indAng = 1, numInAnnulus
            if (indAng .ne. indNay) then
              iptNum = neighbor(indAng)
              angleDiff = abs(thetaNeigh(indAng) - thetaNeigh(indNay))
              if (angleDiff > 180.) angleDiff = 360. -angleDiff
              ibin = min(numBins, ifix(numBins * angleDiff / 180. +1.))
              ifOut = iunitAllDist
              do jgrf = 1, neededAsNeigh
                jj = igraphNeigh(jgrf)
                if (isAngPtInGraph((iptNum - 1) * numGraphs + jj)) then
                  graphs(ibin, jj) = graphs(ibin, jj) + 1.
                  if (ifOut > 0) then
                    write(iunitAllDist, '(3i10, f11.3)') indRef, neighbor(indNay), &
                        iptNum, angleDiff
                    ifOut = 0
                  endif
                endif
              enddo
            endif
          enddo
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
  binArea = 3.14159 * (radMaxSq - radMinSq) / (2. *numBins)
  do ibin = 1, numBins
    do jj = 1, numGraphs
      totalArea = binArea * fracSum(ibin, jj)
      if (totalArea == 0.) then
        graphs(ibin, jj) = 0.
      else
        graphs(ibin, jj) = graphs(ibin, jj) / totalArea
      endif
      fracSum(ibin, jj) = totalArea
    enddo
  enddo
  return
end subroutine angDist

