! $Id$

subroutine shiftSurface(isurf, shift, idir, vertShifted, ifZeroedShift, zScale)
  use mtkvars
  real*4 shift(3)
  logical*1 vertShifted(*)
  include 'model.inc90'
  xShift = idir * shift(1)
  yShift = idir * shift(2)
  zShift = idir * shift(3)
  !
  ! loop through contours, adjusting min/max's and all points
  !
  do icont = indStartCont(isurf), indStartCont(isurf) + numContInSurf(isurf) - 1
    iobj = listCont(icont)
    ibase = ibase_obj(iobj)
    contXmin(icont) = contXmin(icont) + xShift
    contXmax(icont) = contXmax(icont) + xShift
    contYmin(icont) = contYmin(icont) + yShift
    contYmax(icont) = contYmax(icont) + yShift
    contZval(icont) = contZval(icont) + zShift
    do ipt = 1, npt_in_obj(iobj)
      ipoint = object(ibase + ipt)
      p_coord(1, ipoint) = p_coord(1, ipoint) + xShift
      p_coord(2, ipoint) = p_coord(2, ipoint) + yShift
      p_coord(3, ipoint) = zScale * nint((p_coord(3, ipoint) + zShift) / zScale)
    enddo
  enddo
  call shiftSurfMesh(isurf, shift, idir, vertShifted, ifZeroedShift, zScale)
  return
end subroutine shiftSurface


subroutine shiftSurfMesh(isurf, shift, idir, vertShifted, ifZeroedShift, zScale)
  use mtkvars
  real*4 shift(3)
  logical*1 vertShifted(*)
  !
  ! set up array to keep track of shifted vertices
  !
  if (ifZeroedShift == 0) then
    do i = 1, numVerts
      vertShifted(i) = .false.
    enddo
  endif
  ifZeroedShift = 1
  xShift = idir * shift(1)
  yShift = idir * shift(2)
  zShift = idir * shift(3)
  surfXmin(isurf) = surfXmin(isurf) + xShift
  surfXmax(isurf) = surfXmax(isurf) + xShift
  surfYmin(isurf) = surfYmin(isurf) + yShift
  surfYmax(isurf) = surfYmax(isurf) + yShift
  surfZmin(isurf) = surfZmin(isurf) + zShift
  surfZmax(isurf) = surfZmax(isurf) + zShift
  !
  ! loop on polygons in surface
  !
  do list = indStartSurf(isurf), &
      indStartSurf(isurf) + numInSurf(isurf) - 1
    ipoly = listSurf(list)
    polyXmin(ipoly) = polyXmin(ipoly) + xShift
    polyXmax(ipoly) = polyXmax(ipoly) + xShift
    polyYmin(ipoly) = polyYmin(ipoly) + yShift
    polyYmax(ipoly) = polyYmax(ipoly) + yShift
    polyZmin(ipoly) = polyZmin(ipoly) + zShift
    polyZmax(ipoly) = polyZmax(ipoly) + zShift
    !
    ! loop on triangles in polygon
    !
    do indTri = indStartPoly(ipoly), &
        indStartPoly(ipoly) + numInPoly(ipoly) - 1
      triangXmin(indTri) = triangXmin(indTri) + xShift
      triangXmax(indTri) = triangXmax(indTri) + xShift
      triangYmin(indTri) = triangYmin(indTri) + yShift
      triangYmax(indTri) = triangYmax(indTri) + yShift
      triangZmin(indTri) = triangZmin(indTri) + zShift
      triangZmax(indTri) = triangZmax(indTri) + zShift
      !
      ! get rotated shift to add to rotated coordinates
      !
      temp = xShift * cosGamma(indTri) - yShift * sinGamma(indTri)
      xShiftRot = temp * cosBeta(indTri) - zShift * sinBeta(indTri)
      yShiftRot = xShift * sinGamma(indTri) + yShift * cosGamma(indTri)
      zShiftRot = temp * sinBeta(indTri) + zShift * cosBeta(indTri)
      triZrot(indTri) = triZrot(indTri) + zShiftRot
      do i = 1, 3
        triXYrot(i, 1, indTri) = triXYrot(i, 1, indTri) + xShiftRot
        triXYrot(i, 2, indTri) = triXYrot(i, 2, indTri) + yShiftRot
        !
        ! shift vertices if not already done
        !
        ind = indVert(i, indTri)
        if (.not.vertShifted(ind)) then
          verts(1, ind) = verts(1, ind) + xShift
          verts(2, ind) = verts(2, ind) + yShift
          verts(3, ind) = zScale * nint((verts(3, ind) + zShift) / zScale)
          vertShifted(ind) = .true.
        endif
      enddo
    enddo
  enddo
  return
end subroutine shiftSurfMesh

subroutine shiftLine(iobjRef, indRef, limRef, delatX, deltaY, deltaZ, &
    xModPt, yModPt, zModPt, xMin, xMax, yMin, yMax, zMin, zMax, &
    globalXmin, globalXmax, globalYmin, globalYmax, globalZmin, globalZmax, zScale)
  real*4 xModPt(*), yModPt(*), zModPt(*), xMin(*), xMax(*), yMin(*), yMax(*)
  real*4 zMin(*), zMax(*), globalXmin(*), globalXmax(*), globalYmin(*)
  real*4  globalYmax(*), globalZmin(*), globalZmax(*)
  globalXmin(iobjRef) = globalXmin(iobjRef) + delatX
  globalYmin(iobjRef) = globalYmin(iobjRef) + deltaY
  globalZmin(iobjRef) = globalZmin(iobjRef) + deltaZ
  globalXmax(iobjRef) = globalXmax(iobjRef) + delatX
  globalYmax(iobjRef) = globalYmax(iobjRef) + deltaY
  globalZmax(iobjRef) = globalZmax(iobjRef) + deltaZ
  do iref = indRef, limRef
    xMin(iref) = xMin(iref) + delatX
    yMin(iref) = yMin(iref) + deltaY
    zMin(iref) = zMin(iref) + deltaZ
    xMax(iref) = xMax(iref) + delatX
    yMax(iref) = yMax(iref) + deltaY
    zMax(iref) = zMax(iref) + deltaZ
  enddo
  do iref = indRef, limRef + 1
    xModPt(iref) = xModPt(iref) + delatX
    yModPt(iref) = yModPt(iref) + deltaY
    zModPt(iref) = zScale * nint((zModPt(iref) + deltaZ) / zScale)
  enddo
  return
end subroutine shiftLine

subroutine shiftPoint(iobjRef, iref, delatX, deltaY, deltaZ, &
    xModPt, yModPt, zModPt, xMin, xMax, yMin, yMax, zMin, zMax, &
    globalXmin, globalXmax, globalYmin, globalYmax, globalZmin, globalZmax, zScale)
  real*4 xModPt(*), yModPt(*), zModPt(*), xMin(*), xMax(*), yMin(*), yMax(*)
  real*4 zMin(*), zMax(*), globalXmin(*), globalXmax(*), globalYmin(*)
  real*4  globalYmax(*), globalZmin(*), globalZmax(*)
  xMin(iref) = xMin(iref) + delatX
  yMin(iref) = yMin(iref) + deltaY
  zMin(iref) = zMin(iref) + deltaZ
  xMax(iref) = xMax(iref) + delatX
  yMax(iref) = yMax(iref) + deltaY
  zMax(iref) = zMax(iref) + deltaZ
  xModPt(iref) = xModPt(iref) + delatX
  yModPt(iref) = yModPt(iref) + deltaY
  zModPt(iref) = zScale * nint((zModPt(iref) + deltaZ) / zScale)
  globalXmin(iobjRef) = min(globalXmin(iobjRef), xMin(iref))
  globalYmin(iobjRef) = min(globalYmin(iobjRef), yMin(iref))
  globalZmin(iobjRef) = min(globalZmin(iobjRef), zMin(iref))
  globalXmax(iobjRef) = max(globalXmax(iobjRef), xMax(iref))
  globalYmax(iobjRef) = max(globalYmax(iobjRef), yMax(iref))
  globalZmax(iobjRef) = max(globalZmax(iobjRef), zMax(iref))
  return
end subroutine shiftPoint

subroutine unshift_object(ityepShift, iobjFlag, xModPt, yModPt, &
    zModPt, itype, numWobj, indStart, numPtInObj, xyScale, zScale, &
    zGapStart, zGapEnd, numGaps)
  use mtkvars
  real*4 xModPt(*), yModPt(*), zModPt(*), zGapStart, zGapEnd
  integer*4 itype(*), indStart(*), numPtInObj(*)
  include 'model.inc90'
  integer*4 modInd(LIM_MESH_IND)
  !
  ifInList = 0
  do imo = 1, numObjShifted
    if (iobjShift(imo) == ityepShift) ifInList = imo
  enddo
  if (ifInList == 0) then
    print *,'This object is not currently shifted'
    return
  endif
  ishift = indStartShift(ifInList)
  if (iobjFlag == 1) then
    !
    ! unshift lines
    !
    do iobj = 1, numWobj
      if (itype(iobj) == ityepShift) then
        if (shifted(ishift)) then
          delatX = -shifts(1, ishift)
          deltaY = -shifts(2, ishift)
          deltaZ = -shifts(3, ishift)
          do ind = indStart(iobj), indStart(iobj) + numPtInObj(iobj) - 1
            xModPt(ind) = xModPt(ind) + delatX
            yModPt(ind) = yModPt(ind) + deltaY
            zModPt(ind) = zScale * nint((zModPt(ind) + deltaZ) / zScale)
          enddo
          shifts(1, ishift) = 0.
          shifts(2, ishift) = 0.
          shifts(3, ishift) = 0.
        endif
        ishift = ishift + 1
      endif
    enddo
  elseif (iobjFlag == 2) then
    !
    ! unshift points
    !
    do iobj = 1, numWobj
      if (itype(iobj) == ityepShift) then
        do ind = indStart(iobj), indStart(iobj) + numPtInObj(iobj) - 1
          if (shifted(ishift)) then
            delatX = -shifts(1, ishift)
            deltaY = -shifts(2, ishift)
            deltaZ = -shifts(3, ishift)
            xModPt(ind) = xModPt(ind) + delatX
            yModPt(ind) = yModPt(ind) + deltaY
            zModPt(ind) = zScale * nint((zModPt(ind) + deltaZ) / zScale)
            shifts(1, ishift) = 0.
            shifts(2, ishift) = 0.
            shifts(3, ishift) = 0.
          endif
          ishift = ishift + 1
        enddo
      endif
    enddo
  else
    !
    ! unshift meshes
    !
    indMesh = 0
    do i = 1, numMeshLoaded
      if (iobjMesh(i) == ityepShift) indMesh = i
    enddo
    if (indMesh == 0) then
      numMeshLoaded = 0
      numVerts = 0
      numTriang = 0
      numPoly = 0
      numSurf = 0
      numConts = 0
      call process_mesh(ityepShift, xyScale, zScale, modInd, &
          zGapStart, zGapEnd, numGaps, ifErr, 0)
      if (ifErr .ne. 0) return
      indMesh = 1
    endif
    ifZeroedShift = 0
    do isurf = iobjSurf(indMesh), iobjSurf(indMesh) + numSurfObj(indMesh) - 1
      if (shifted(ishift)) then
        call shiftSurface(isurf, shifts(1, ishift), -1, modInd, &
            ifZeroedShift, zScale)
        shifts(1, ishift) = 0.
        shifts(2, ishift) = 0.
        shifts(3, ishift) = 0.
      endif
      ishift = ishift + 1
    enddo

  endif
  return
end subroutine unshift_object

subroutine check_line_surface(isurf, indRef, limRef, &
    xModPt, yModPt, zModPt, distMin, distAbs, zScale)
  use mtkvars
  real*4 xModPt(*), yModPt(*), zModPt(*)
  logical inside_wimpobj, crosses_wimpobj
  !
  ! loop on line segments, check against bounding box of surface
  !
  sLimit = .01 * zScale
  do iref = indRef, limRef
    refXmin = xMin(iref)
    refYmin = yMin(iref)
    refZmin = zMin(iref)
    refXmax = xMax(iref)
    refYmax = yMax(iref)
    refZmax = zMax(iref)
    if (refXmin - surfXmax(isurf) < distMin .and. &
        surfXmin(isurf) - refXmax < distMin .and. &
        refYmin - surfYmax(isurf) < distMin .and. &
        surfYmin(isurf) - refYmax < distMin .and. &
        refZmin - surfZmax(isurf) < distMin .and. &
        surfZmin(isurf) - refZmax < distMin) then
      !
      ! better search now for contours that overlap
      !
      do icont = indStartCont(isurf), indStartCont(isurf) + numContInSurf(isurf) - 1
        if (refXmin <= contXmax(icont) .and. &
            contXmin(icont) <= refXmax .and. &
            refYmin <= contYmax(icont) .and. &
            contYmin(icont) <= refYmax .and. &
            refZmin <= contZval(icont) + sLimit .and. &
            contZval(icont) - sLimit <= refZmax) then
          iobj = listCont(icont)
          if (refZmax - refZmin > sLimit) then
            !
            ! if line passes through the zplane, get point on the
            ! contour's plane and check if inside
            !
            frac = (contZval(icont) - zModPt(iref)) / (zModPt(iref + 1) - zModPt(iref))
            xInterp = xModPt(iref) + frac * (xModPt(iref + 1) - xModPt(iref))
            yInterp = yModPt(iref) + frac * (yModPt(iref + 1) - yModPt(iref))
            if (inside_wimpobj(iobj, xInterp, yInterp)) then
              distMin = -1.
              return
            endif
          else
            !
            ! otherwise, need to check for 2-d intersection between the
            ! line and each of the segments of the contour
            !
            if (crosses_wimpobj(iobj, xModPt(iref), yModPt(iref), &
                xModPt(iref + 1), yModPt(iref + 1))) then
              distMin = -1.
              return
            endif
          endif
        endif
      enddo

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
              call segment_to_triangle(xModPt(iref), &
                  yModPt(iref), zModPt(iref), xModPt(iref + 1), &
                  yModPt(iref + 1), zModPt(iref + 1), &
                  indTri, ttRef, xrot, yrot, dist)
              if (dist < distMin) then
                distMin = dist
                if (dist < distAbs) return
              endif
            endif
          enddo
        endif
      enddo
    endif
  enddo
  return
end subroutine check_line_surface

!
subroutine check_point_surface(isurf, xModPt, yModPt, zModPt, xMinIn, xMaxIn, &
    yMinIn, yMaxIn, zMinIn, zMaxIn, pointSize, distMin, distAbs, zScale)
  use mtkvars
  logical inside_wimpobj
  !
  ! loop on contours, see if point is inside
  !
  sLimit = .01 * zScale
  do icont = indStartCont(isurf), indStartCont(isurf) + numContInSurf(isurf) - 1
    if (contXmin(icont) <= xModPt .and. &
        contXmax(icont) >= xModPt .and. &
        contYmin(icont) <= yModPt .and. &
        contYmax(icont) >= yModPt .and. &
        abs(contZval(icont) - zModPt) <= sLimit) then
      iobj = listCont(icont)
      if (inside_wimpobj(iobj, xModPt, yModPt)) then
        distMin = -1.
        return
      endif
    endif
  enddo
  !
  ! find minimum distance of point from surface now
  !
  do list = indStartSurf(isurf), &
      indStartSurf(isurf) + numInSurf(isurf) - 1
    ipoly = listSurf(list)
    !
    ! check against polygon global limits
    !
    if (xMinIn - polyXmax(ipoly) < distMin .and. &
        polyXmin(ipoly) - xMaxIn < distMin .and. &
        yMinIn - polyYmax(ipoly) < distMin .and. &
        polyYmin(ipoly) - yMaxIn < distMin .and. &
        zMinIn - polyZmax(ipoly) < distMin .and. &
        polyZmin(ipoly) - zMaxIn < distMin) then
      !
      ! scan through triangles
      !
      do indTri = indStartPoly(ipoly), indStartPoly(ipoly) + &
          numInPoly(ipoly) - 1
        if (xMinIn - triangXmax(indTri) < distMin .and. &
            triangXmin(indTri) - xMaxIn < distMin .and. &
            yMinIn - triangYmax(indTri) < distMin .and. &
            triangYmin(indTri) - yMaxIn < distMin .and. &
            zMinIn - triangZmax(indTri) < distMin .and. &
            triangZmin(indTri) - zMaxIn < distMin) then
          !
          !
          call point_to_triangle(xModPt, yModPt, zModPt, indTri, xrot, yrot, dist)
          dist = dist - pointSize
          if (dist < distMin) then
            distMin = dist
            if (dist < distAbs) return
          endif
        endif
      enddo
    endif
  enddo
  return
end subroutine check_point_surface


subroutine check_two_meshes(isurf, jsurf, distMin, distAbs, zScale)
  use mtkvars
  include 'model.inc90'
  logical crosses_wimpobj, inside_wimpobj
  sLimit = .1 * zScale
  !
  ! first check contours against each other to find overlap
  !
  refXmin = surfXmin(jsurf)
  refYmin = surfYmin(jsurf)
  refZmin = surfZmin(jsurf)
  refXmax = surfXmax(jsurf)
  refYmax = surfYmax(jsurf)
  refZmax = surfZmax(jsurf)
  do icont = indStartCont(isurf), indStartCont(isurf) + numContInSurf(isurf) - 1
    if (refXmin <= contXmax(icont) .and. &
        contXmin(icont) <= refXmax .and. &
        refYmin <= contYmax(icont) .and. &
        contYmin(icont) <= refYmax .and. &
        refZmin <= contZval(icont) + sLimit .and. &
        contZval(icont) - sLimit <= refZmax) then
      coXmin = contXmin(icont)
      coYmin = contYmin(icont)
      coXmax = contXmax(icont)
      coYmax = contYmax(icont)
      do jcont = indStartCont(jsurf), indStartCont(jsurf) + numContInSurf(jsurf) - 1
        if (abs(contZval(icont) - contZval(jcont)) < sLimit .and. &
            contXmin(jcont) <= coXmax .and. &
            coXmin <= contXmax(jcont) .and. &
            contYmin(jcont) <= coYmax .and. &
            coYmin <= contYmax(jcont)) then
          iobj = listCont(icont)
          jobj = listCont(jcont)
          !
          ! check each segment of j contour to see if crosses i
          !
          jbase = ibase_obj(jobj)
          jObjNumPt = npt_in_obj(jobj)
          ipt = object(jbase + jObjNumPt)
          x1 = p_coord(1, ipt)
          y1 = p_coord(2, ipt)
          jpt = object(ibase_obj(iobj) + 1)
          x2 = p_coord(1, jpt)
          y2 = p_coord(2, jpt)
          !
          ! but first check if one point of either is inside the other
          ! this handles case of completely enclosed contours
          !
          if (inside_wimpobj(iobj, x1, y1) .or. &
              inside_wimpobj(jobj, x2, y2)) then
            distMin = -1.
            return
          endif
          !
          do j = 1, jObjNumPt
            jpt = object(jbase + j)
            x2 = p_coord(1, jpt)
            y2 = p_coord(2, jpt)
            if (min(x1, x2) <= coXmax .and. coXmin <= max(x1, x2) .and. &
                min(y1, y2) <= coYmax .and. coYmin <= max(y1, y2)) then
              if (crosses_wimpobj(iobj, x1, y1, x2, y2)) then
                distMin = -1.
                return
              endif
            endif
            x1 = x2
            y1 = y2
          enddo
        endif
      enddo
    endif
  enddo
  !
  ! now find actual distance between surfaces
  !
  !
  ! loop on polygons of neighbor surface
  !
  refXmin = surfXmin(isurf)
  refYmin = surfYmin(isurf)
  refZmin = surfZmin(isurf)
  refXmax = surfXmax(isurf)
  refYmax = surfYmax(isurf)
  refZmax = surfZmax(isurf)
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
        xmaxNeigh = triangXmax(jtri)
        yMaxNeigh = triangYmax(jtri)
        zMaxNeigh = triangZmax(jtri)
        if (refXmin - xmaxNeigh < distMin .and. &
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
                polyXmin(ipoly) - xmaxNeigh < distMin .and. &
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
                    triangXmin(indTri) - xmaxNeigh < distMin .and. &
                    yMinNeigh - triangYmax(indTri) < distMin .and. &
                    triangYmin(indTri) - yMaxNeigh < distMin .and. &
                    zMinNeigh - triangZmax(indTri) < distMin .and. &
                    triangZmin(indTri) - zMaxNeigh < distMin) then
                  call triangle_to_triangle(indTri, jtri, &
                      xRot1, yRot1, zRot1, xRot2, yROt2, indTriRot, dist)
                  if (dist < distMin) then
                    distMin = dist
                    if (dist < distAbs) return
                  endif
                endif
              enddo
            endif
          enddo
        endif
      enddo
    endif
  enddo
  return
end subroutine check_two_meshes



logical function inside_wimpobj (iobj, x, y)
  include 'model.inc90'
  !
  inside_wimpobj = .FALSE.
  numPoints = npt_in_obj(iobj)
  IF (numPoints <= 2) RETURN
  ibase = ibase_obj(iobj)
  ipt = object(ibase + numPoints)
  x0 = x - p_coord(1, ipt)
  y0 = y - p_coord(2, ipt)
  IF (x0 == 0. .AND. y0 == 0.) then
    inside_wimpobj = .TRUE.
    RETURN
  endif
  sumAng = 0.
  DO j = 1, numPoints
    ipt = object(ibase + j)
    xa = x - p_coord(1, ipt)
    ya = y - p_coord(2, ipt)
    IF (xa == 0. .AND. ya == 0.) then
      inside_wimpobj = .TRUE.
      RETURN
    endif
    sinAng = (ya * x0 - xa * y0)
    cosAng = (xa * x0 + ya * y0)
    angle = atan2(sinAng, cosAng)
    sumAng = sumAng + angle
    x0 = xa
    y0 = ya
  enddo
  IF (ABS(sumAng) >= 4.0) inside_wimpobj = .TRUE.
  RETURN
end function inside_wimpobj


logical function crosses_wimpobj (iobj, x1Start, y1Start, x1End, y1End)
  include 'model.inc90'
  !
  crosses_wimpobj = .FALSE.
  numPoints = npt_in_obj(iobj)
  IF (numPoints <= 1) RETURN
  ibase = ibase_obj(iobj)
  ipt = object(ibase + numPoints)
  x2Start = p_coord(1, ipt)
  y2Start = p_coord(2, ipt)
  dx1 = x1End - x1Start
  dy1 = y1End - y1Start
  crosses_wimpobj = .true.
  DO j = 1, numPoints
    ipt = object(ibase + j)
    x2End = p_coord(1, ipt)
    y2End = p_coord(2, ipt)
    dx2 = x2Start - x2End
    dy2 = y2Start - y2End
    dxs = x2Start - x1Start
    dys = y2Start - y1Start
    denom = dx1 * dy2 - dy1 * dx2
    tNumer = dxs * dy2 - dys * dx2
    uNumer = dx1 * dys - dy1 * dxs
    if (denom < 0) then
      if (tNumer <= 0. .and. uNumer <= 0. .and. tNumer >= denom .and. uNumer >= denom) &
          return
    else
      if (tNumer >= 0. .and. uNumer >= 0. .and. tNumer <= denom .and. uNumer <= denom) &
          return
    endif
    x2Start = x2End
    y2Start = y2End
  enddo
  crosses_wimpobj = .false.
  RETURN
end function crosses_wimpobj


subroutine getRanShift(numTrialsIn, numTrialCycle, cycleFactor, ranMin, &
    ranMax, ranZrel, curShift, deltaXYZ, testXmin, testXmax, testYmin, &
    testYmax, testZmin, testZmax, shiftXmin, shiftXmax, shiftYmin, &
    shiftYmax, shiftZmin, shiftZmax, zScale, iseed, numBound, listBound, &
    zBound, xModPt, yModPt, zModPt, indStart, indEnd)
  use mtkvars
  real*4 curShift(3), deltaXYZ(3), xModPt(*), yModPt(*), zModPt(*), zBound(*)
  integer*4 listBound(*)
  integer*2 numTrialsIn
  include 'model.inc90'
  logical badShift, outside_boundary
  real*4 b3dran
  !
  numAttempt = 0
  numBoundChecks = 0
  numTrialsIn = numTrialsIn + 1
  numChange = numTrialsIn / numTrialCycle
  changeFac = cycleFactor**numChange
  ranMaxTemp = ranMax * changeFac
  ranMinTemp = ranMin * changeFac
  ranZtemp = ranMaxTemp * ranZrel
  ranZmin = ranMinTemp * ranZrel
  badShift = .true.
  do while(badShift)
    x0 = b3dran(iseed) + b3dran(iseed) + b3dran(iseed) + b3dran(iseed)
    shiftX = ranMaxTemp * (2. *b3dran(iseed) - 1.)
    x0 = b3dran(iseed) + b3dran(iseed) + b3dran(iseed) + b3dran(iseed)
    shiftY = ranMaxTemp * (2. *b3dran(iseed) - 1.)
    x0 = b3dran(iseed) + b3dran(iseed) + b3dran(iseed) + b3dran(iseed)
    shiftZ = zScale * nint(ranZtemp * (2. *b3dran(iseed) - 1.) / zScale)
    distShift = shiftX**2 + shiftY**2
    delatX = shiftX - curShift(1)
    deltaY = shiftY - curShift(2)
    deltaZ = shiftZ - curShift(3)
    !
    ! if shift is out of range or moves the item outside the
    ! bounding box, repeat without counting a trial
    !
    badShift = distShift < ranMinTemp**2 .or. &
        distShift > ranMaxTemp**2 .or. &
        abs(shiftZ) < ranZmin .or. abs(shiftZ) > ranZtemp .or. &
        testXmin + delatX < shiftXmin .or. &
        testXmax + delatX > shiftXmax .or. &
        testYmin + deltaY < shiftYmin .or. &
        testYmax + deltaY > shiftYmax .or. &
        testZmin + deltaZ < shiftZmin .or. &
        testZmax + deltaZ > shiftZmax
    numAttempt = numAttempt + 1
    if (numAttempt > 10000) then
      print *,'too many trials with bad shift'
      deltaXYZ(1) = 0.
      deltaXYZ(2) = 0.
      deltaXYZ(3) = 0.
      return
    endif
    if (.not.badShift .and. numBound > 0) then
      numBoundChecks = numBoundChecks + 1
      badShift = outside_boundary(numBound, listBound, &
          zBound, xModPt, yModPt, zModPt, indStart, indEnd, delatX, deltaY, deltaZ)

      if (badShift .and. numBoundChecks > 500) then
        print *,'too many shifts outside boundaries'
        deltaXYZ(1) = 0.
        deltaXYZ(2) = 0.
        deltaXYZ(3) = 0.
        return
      endif
    endif
  enddo
  !
  curShift(1) = shiftX
  curShift(2) = shiftY
  curShift(3) = shiftZ
  deltaXYZ(1) = delatX
  deltaXYZ(2) = deltaY
  deltaXYZ(3) = deltaZ
  return
end subroutine getRanShift


! OUTSIDE_BOUNDARY tests whether a point, line, or surface, shifted
! by delx, dely, delz, crosses or is outside the boundary contours

logical function outside_boundary(numBound, listBound, &
    zBound, xModPt, yModPt, zModPt, indStart, indEnd, delatX, deltaY, deltaZ)
  use mtkvars
  real*4 xModPt(*), yModPt(*), zModPt(*), zBound(*)
  integer*4 listBound(*)
  include 'model.inc90'
  logical badShift, crosses_wimpobj, inside_wimpobj
  badShift = .false.

  if (indEnd > 0) then
    !
    ! check single point or line of points for whether they are
    ! all inside the nearest boundary in z
    !
    ind = indStart
    do while(.not.badShift .and. ind <= indEnd)
      xPoint = xModPt(ind) + delatX
      yPoint = yModPt(ind) + deltaY
      zPoint = zModPt(ind) + deltaZ
      icont = nearest_boundary(zPoint, zBound, numBound)
      badShift = .not.inside_wimpobj(listBound(icont), xPoint, yPoint)
      ind = ind + 1
    enddo
  else
    !
    ! for mesh, loop on contours of surface, checking against
    ! nearest boundary contour in z
    !
    ind = indStartCont(indStart)
    do while(.not.badShift .and. ind < &
        indStartCont(indStart) + numContInSurf(indStart))
      jobj = listCont(ind)
      zPoint = contZval(ind) + deltaZ
      icont = nearest_boundary(zPoint, zBound, numBound)
      iobj = listBound(icont)
      jbase = ibase_obj(jobj)
      jObjNumPt = npt_in_obj(jobj)
      ipt = object(jbase + jObjNumPt)
      x1 = p_coord(1, ipt) + delatX
      y1 = p_coord(2, ipt) + deltaY
      !
      ! first check that at least one point is inside the boundary
      ! then check each segment for crossing the boundary
      !
      badShift = .not.inside_wimpobj(iobj, x1, y1)
      !
      j = 1
      do while(.not.badShift .and. j <= jObjNumPt)
        jpt = object(jbase + j)
        x2 = p_coord(1, jpt) + delatX
        y2 = p_coord(2, jpt) + deltaY
        badShift = crosses_wimpobj(iobj, x1, y1, x2, y2)
        x1 = x2
        y1 = y2
        j = j + 1
      enddo
      !
      ind = ind + 1
    enddo
  endif
  outside_boundary = badShift
  return
end function outside_boundary


subroutine get_random_boundary(iobjBound, listBound, zBound, &
    numBound, limBound)
  real*4 zBound(*)
  integer*4 listBound(*)
  include 'model.inc90'
  numBound = 0
  if (iobjBound <= 0) return
  do iobj = 1, max_mod_obj
    if (256 - obj_color(2, iobj) == iobjBound .and. npt_in_obj(iobj) > 2) then
      if (numBound == limBound) then
        print *,'WARNING: NOT ALL BOUNDARY CONTOURS WERE ', &
            'LOADED - TOO MANY FOR ARRAYS'
        return
      endif
      numBound = numBound + 1
      listBound(numBound) = iobj
      zBound(numBound) = p_coord(3, object(ibase_obj(iobj) + 1))
    endif
  enddo
  return
end subroutine get_random_boundary

function nearest_boundary(z, zBound, numBound)
  real*4 zBound(*)
  delZmin = 1.e20
  do i = 1, numBound
    if (abs(z - zBound(i)) < delZmin) then
      delZmin = abs(z - zBound(i))
      nearest_boundary = i
    endif
  enddo
  return
end function nearest_boundary
