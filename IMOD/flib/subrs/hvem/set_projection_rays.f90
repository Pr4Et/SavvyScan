! $Id$
!
! !
! Sets up information about projection rays at an angle whose sine and
! cosine are [sinang] and [cosang], through
! an image slice whose size is [nxslice] by [nyslice].  The number of
! rays, i.e. the size of centered output in X, is given by [nxout].
! The arrays [xraystr] and [yraystr] are returned with the coordinates
! at which to start each ray, the array [nrayinc] is returned with the
! number of points in each ray, and [nraymax] is returned with the
! maximum number of points.  [sinang] and [cosang] may be modified for
! projections very near vertical or horizontal.
! !
subroutine set_projection_rays(sinAng, cosAng, nxSlice, nySlice, nxOut, xrayStart, &
    yrayStart, numRayIncrem, nrayMax)
  implicit none
  real*4 xrayStart(*), yrayStart(*)
  integer*4 nxSlice, nySlice, nxOut, numRayIncrem(*), nrayMax
  real*4 sinAng, cosAng, xrayTmp, yrayTmp, xgood(4), ygood(4)
  integer*4 nrayTmp, numGoodInter, indGood, ixOut
  real*4 tanAng, rayIntercept, xleft, xright, ybot, ytop, yleft, yright, xtop, xbot
  logical b3dxor
  !
  nrayMax = 0
  do ixOut = 1, nxOut
    nrayTmp = 0
    xrayTmp = 0.
    yrayTmp = 0.
    !
    ! if a near-vertical projection, set up to be exactly vertical
    ! (limit was 0.01, try it at 0.001)
    if (abs(sinAng) < 0.001) then
      sinAng = 0.
      cosAng = sign(1., cosAng)
      if (cosAng > 0.) then
        yrayTmp = 2.
        xrayTmp = nxSlice / 2 + ixOut - nxOut / 2
      else
        yrayTmp = nySlice-1.
        xrayTmp = nxSlice / 2 + 1 + nxOut / 2 - ixOut
      endif
      if (xrayTmp >= 1 .and. xrayTmp <= nxSlice) &
          nrayTmp = nySlice-2
      !
      ! if a near-horizontal projection, set up to be exactly horizontal
      !
    elseif (abs(cosAng) < 0.001) then
      sinAng = sign(1., sinAng)
      cosAng = 0.
      if (sinAng > 0.) then
        xrayTmp = 2.
        yrayTmp = nySlice / 2 + 1 + nxOut / 2 - ixOut
      else
        xrayTmp = nxSlice-1.
        yrayTmp = nySlice / 2 + ixOut - nxOut / 2
      endif
      if (yrayTmp >= 1 .and. yrayTmp <= nySlice) &
          nrayTmp = nxSlice-2
      !
      ! otherwise need to look at intersections with slice box
      !
    else
      numGoodInter = 0
      tanAng = sinAng / cosAng
      rayIntercept = -(ixOut - 0.5 - nxOut / 2) / sinAng
      !
      ! coordinates of edges of slice box
      !
      xleft = 1.01 - (nxSlice / 2)
      xright = nxSlice - 1.01 - (nxSlice / 2)
      ybot = 1.01 - (nySlice / 2)
      ytop = nySlice - 1.01 - (nySlice / 2)
      !
      ! corresponding intersections of the ray with extended edges
      !
      yleft = xleft / tanAng + rayIntercept
      yright = xright / tanAng + rayIntercept
      xbot = (ybot - rayIntercept) * tanAng
      xtop = (ytop - rayIntercept) * tanAng
      !
      ! make list of intersections that are actually within slice box
      !
      if (yleft >= ybot .and. yleft <= ytop) then
        numGoodInter = numGoodInter + 1
        xgood(numGoodInter) = xleft
        ygood(numGoodInter) = yleft
      endif
      if (yright >= ybot .and. yright <= ytop) then
        numGoodInter = numGoodInter + 1
        xgood(numGoodInter) = xright
        ygood(numGoodInter) = yright
      endif
      if (xbot >= xleft .and. xbot <= xright) then
        numGoodInter = numGoodInter + 1
        xgood(numGoodInter) = xbot
        ygood(numGoodInter) = ybot
      endif
      if (xtop >= xleft .and. xtop <= xright) then
        numGoodInter = numGoodInter + 1
        xgood(numGoodInter) = xtop
        ygood(numGoodInter) = ytop
      endif
      !
      ! if there are real intersections, use them to set up ray start
      !
      if (numGoodInter > 0) then
        indGood = 1
        if ((ygood(2) < ygood(1)) .neqv. (cosAng < 0)) indGood = 2
        xrayTmp = xgood(indGood) + nxSlice / 2 + 0.5
        yrayTmp = ygood(indGood) + nySlice / 2 + 0.5
        nrayTmp = 1 + &
            sqrt((xgood(1) - xgood(2))**2 + (ygood(1) - ygood(2))**2)
        if (nrayTmp < 3) nrayTmp = 0
      endif
    endif
    !
    ! store final ray information
    !
    xrayStart(ixOut) = xrayTmp
    yrayStart(ixOut) = yrayTmp
    numRayIncrem(ixOut) = nrayTmp
    nrayMax = max(nrayMax, nrayTmp)
  enddo
  return
end subroutine set_projection_rays
