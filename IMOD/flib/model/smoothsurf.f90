! * * * * * * SMOOTHSURF * * * * * * *
!
! This program will smooth a surface defined by model contours,
! adjusting the positions of points in each contour as needed.  At each
! point on the surface, it fits a 3-D polynomial to all points within
! a defined range of Z-levels and within a specified distance of the
! central point.  That point's position is then replaced by the fitted
! position from the polynomial.  After this surface smoothing
! operation, each contour is independently smoothed by local fitting to
! 2-D (ordinary) polynomials.
!
! See man page for more details
!
! David Mastronarde, 9/9/97
!
! $Id$
!
program smoothsurf
  implicit none
  integer idim, LIMPOINTS, iDzLimit, LIMFLAGS
  include 'model.inc90'
  parameter (idim = 20000, LIMPOINTS = 100000, iDzLimit = 5000, LIMFLAGS = 20000)
  !
  real*4 xtmp(LIMPOINTS), ytmp(LIMPOINTS), ztmp(LIMPOINTS), pt(LIMPOINTS,3)
  real*4 p_new(2,max_pt)
  equivalence (xtmp, pt), (ytmp, pt(1, 2)), (ztmp, pt(1, 3))
  !
  character*320 inputFile, outputFile
  !
  logical readw_or_imod, failed, getModelObjectRange, retainMeshes
  include 'statsize.inc'
  real*4 xmat(msiz,idim), xm(msiz), sd(msiz), ssd(msiz,msiz), b1(msiz)
  integer*4 izObj(max_obj_num), isurf(max_obj_num)
  integer*4 iobjBest(-iDzLimit:iDzLimit), iptBest(-iDzLimit:iDzLimit)
  integer*4 iobjDo(LIMFLAGS), iflags(LIMFLAGS)
  logical*1 objectOK(max_obj_num)

  integer*4 minPoints, numObjDo, iorder, numToFitInZ, iorder2, ierr, ifFlipped
  real*4 tolCross, closeThresh, xImScale, yImScale, zImScale, xOffset, yOffset, zOffset
  real*4 sepMin, distLimit, xyScale, zScale, dzSq, distMin, dx, dy, distSq
  integer*4 iobj, numInObj, ibase, ipnt1, i, ifOnList, imodObj, ichk, ipt, ifSurfs
  integer*4 ifSpan, iDzLow, iDzHigh, iDzMin, iDzMax, loop, iz, jobj, jbase, numFitsDone
  integer*4 jpnt, jpt, ipLast, ipNext, numZfit, numFit, ifTooFar, ifSave, j, numFitTotal
  real*4 xx, yy, segLen, sinTheta, cosTheta, xLast, yLast, xrot, yRot, dist, frac
  real*4 xcen, ycen, sep, distLast, distCen, yNew, c1, xMiddle
  real*4 yMiddle, tmp, bIntercept
  integer*4 numInsert, norder, numIndepVars, iCheck, ipnt, idz, numInJ, idir, iptCen
  integer*4 ipc, iflag, numObjSmooth, numTotalPoints, numObjTotal, modelObj, modelCont
  integer*4 getImodHead, indMap, getImodFlags, getImodObjSize, getObjSurfaces
  integer*4 numberInList
  common /bigarr/ p_new
  !
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetInteger
  integer*4 PipGetString, PipGetFloat, PipGetLogical
  integer*4 PipGetInOutFile
  !
  ! fallbacks from ../../manpages/autodoc2man -3 2  smoothsurf
  !
  integer numOptions
  parameter (numOptions = 10)
  character*(40 * numOptions) options(1)
  options(1) = &
      'input:InputFile:FN:@output:OutputFile:FN:@objects:ObjectsToSmooth:LI:@'// &
      'nz:NumberOfSections:I:@distance:MaximumDistance:F:@contorder:ContourOrder:I:@'// &
      'surforder:SurfaceOrder:I:@sort:SortSurfaces:I:@retain:RetainMeshes:B:@'// &
      'help:usage:B:'
  !
  ! defaults
  !
  minPoints = 4
  tolCross = 4.
  closeThresh = 1.
  numObjDo = 0
  iorder = 3
  iorder2 = 2
  numToFitInZ = 7
  distLimit = 15.
  ifSurfs = 1
  retainMeshes = .false.
  !
  ! Pip startup: set error, parse options, check help, set flag if used
  !
  call PipReadOrParseOptions(options, numOptions, 'smoothsurf', &
      'ERROR: SMOOTHSURF - ', .false., 2, 1, 1, numOptArg, &
      numNonOptArg)
  !
  if (PipGetInOutFile('InputFile', 1, ' ', inputFile) .ne. 0) &
      call exitError('No input file specified')
  !
  if (PipGetInOutFile('OutputFile', 2, ' ', outputFile) .ne. 0) &
      call exitError('No output file specified')
  !
  call imodPartialMode(1)
  if (.not.readw_or_imod(inputFile)) call exitError('Reading model')
  !

  if (PipGetString('ObjectsToSmooth', inputFile) == 0) &
      call parseList(inputFile, iobjDo, numObjDo)
  if (numObjDo > LIMFLAGS) call exitError('Too many objects in list for arrays')

  ierr = PipGetInteger('NumberOfSections', numToFitInZ)
  ierr = PipGetInteger('ContourOrder', iorder2)
  ierr = PipGetInteger('SurfaceOrder', iorder)
  ierr = PipGetFloat('MaximumDistance', distLimit)
  ierr = PipGetInteger('SortSurfaces', ifSurfs)
  ierr = PipGetLogical('RetainMeshes', retainMeshes)
  call PipDone()

  numObjTotal = getImodObjSize()
  if (numToFitInZ < 1) call exitError('Number of sections is too small')
  if (numToFitInZ > iDzLimit) call exitError('Number of sections is too large for arrays')
  if (iorder2 < 0 .or. iorder2 > 4) call exitError( &
      'Contour smoothing order is outside of allowed range')
  if (iorder < 1 .or. iorder > 5) call exitError( &
      'Surface smoothing order is outside of allowed range')
  if (distLimit <= 1.) call exitError('Maximum distance is too small')
  !
  ! unflip data if it was flipped
  !
  if (getImodHead(xyScale, zScale, xOffset, yOffset, zOffset, ifFlipped) .ne. 0) &
      call exitError('Reading model header information for scaling')
  !
  ! set minimum separation between points to get at least minpts
  ! each side of center
  !
  sepMin = distLimit / minPoints
  !
  ! Get object flags
  !
  do i = 1, LIMFLAGS
    iflags(i) = 0
  enddo
  ierr = getImodFlags(iflags, LIMFLAGS)
  if (ierr .ne. 0) print *,'Error getting object types, assuming', &
      ' all are closed contours'
  !
  ! loop on model objects and decide whether to load
  !
  do imodObj = 1, numObjTotal
    !
    ! is IMOD object on list and object not scattered
    !
    ifOnList = numberInList(imodObj, iobjDo, numObjDo, 1)
    iflag = 0
    if (imodObj <= LIMFLAGS) iflag = mod(iflags(imodObj), 4)
    if (ifOnList == 1 .and. iflag < 2) then
      !
      ! Mark all objects to include by requiring at leat 2 points and
      ! closed contour or open coplanar contour
      !
      if (.not.getModelObjectRange(imodObj, imodObj)) then
        print *
        print *, 'ERROR: SMOOTHSURF - Loading data for object #', imodObj
        call exit(1)
      endif
      call scale_model(0)
      if (ifFlipped .ne. 0) then
        do i = 1, n_point
          tmp = p_coord(2, i)
          p_coord(2, i) = p_coord(3, i)
          p_coord(3, i) = tmp
        enddo
      endif

      if (ifSurfs .ne. 0) then
        if (getObjSurfaces(imodObj, ifSurfs, isurf) .ne. 0) call exitError( &
            'Getting surface values of contours')
      else
        do iobj = 1, max_mod_obj
          isurf(iobj) = 0
        enddo
      endif
      numObjSmooth = 0
      numTotalPoints = 0
      numFitTotal = 0
      numFitsDone = 0
      do iobj = 1, max_mod_obj
        numInObj = npt_in_obj(iobj)
        ibase = ibase_obj(iobj)
        ipnt1 = object(ibase + 1)
        !
        ! require at least 2 points
        !
        objectOK(iobj) = numInObj > 2
        if (numInObj > 2) izObj(iobj) = nint(p_coord(3, ipnt1))
        !
        ! for open objects, make sure contour is coplanar
        !
        if (objectOK(iobj) .and. iflag == 1) then
          do ipt = 2, numInObj
            ipnt = object(ipt + ibase)
            if (nint(p_coord(3, ipnt)) .ne. izObj(iobj)) &
                objectOK(iobj) = .false.
          enddo
        endif
        !
        ! uncross/straighten out start and end points after removing
        ! duplicates, set Z of object negative if not an OK one
        !
        if (objectOK(iobj)) then
          call elimClose(p_coord, 3, object, ibase, numInObj, closeThresh / 4., 3)
          call uncrossCont(p_coord, 3, object, ibase, numInObj, numInObj, 20., &
              tolCross * 2.)
          ! if (ninobj < npt_in_obj(iobj)) &
          ! print *,'removed', npt_in_obj(iobj) -ninobj, ' from', iobj
          npt_in_obj(iobj) = numInObj
          numObjSmooth = numObjSmooth + 1
          numTotalPoints = numTotalPoints + numInObj
        else
          izObj(iobj) = -100000.
        endif
      enddo
      write(*,102) imodObj, numTotalPoints, numObjSmooth
102   format('Doing object',i5,', ',i8,' points in',i6, &
          ' contours being smoothed')
      !
      ! make copy into new coordinates
      !
      do i = 1, n_point
        p_new(1, i) = p_coord(1, i)
        p_new(2, i) = p_coord(2, i)
      enddo
      !
      do iobj = 1, max_mod_obj
        numInObj = npt_in_obj(iobj)
        ibase = ibase_obj(iobj)
        if (objectOK(iobj)) then
          call objToCont(iobj, obj_color, modelObj, modelCont)
          ! print *,'contour', modcont
          do ipt = 1, numInObj
            ! if (modcont == 1) print *,'point', ipt
            ipnt = object(ipt + ibase)
            !
            ! consider this point in this object; look for the closest
            ! point in objects of the same color (i.e. contours of the
            ! same object) within the z range and the distance limit
            !
            xx = p_coord(1, ipnt)
            yy = p_coord(2, ipnt)
            do i = -numToFitInZ, numToFitInZ
              iobjBest(i) = 0
            enddo
            iobjBest(0) = iobj
            iptBest(0) = ipt
            ifSpan = 0
            iDzLow = -numToFitInZ / 2
            iDzHigh = iDzLow + numToFitInZ - 1
            iDzMin = 0
            iDzMax = 0
            loop = 1
            do while(ifSpan == 0 .and. loop <= 2)
              do idz = iDzLow, iDzHigh
                if (iobjBest(idz) == 0) then
                  dzSq = idz**2
                  iz = izObj(iobj) + idz
                  distMin = distLimit**2
                  do jobj = 1, max_mod_obj
                    if (objectOK(jobj) .and. izObj(jobj) == iz .and. &
                        obj_color(2, jobj) == obj_color(2, iobj) .and. &
                        isurf(iobj) == isurf(jobj)) then
                      jbase = ibase_obj(jobj)
                      do jpt = 1, npt_in_obj(jobj)
                        jpnt = object(jpt + jbase)
                        dx = xx - p_coord(1, jpnt)
                        if (abs(dx) <= distLimit) then
                          dy = yy - p_coord(2, jpnt)
                          if (abs(dy) <= distLimit) then
                            distSq = dx**2 + dy**2 + dzSq
                            if (distSq < distMin) then
                              distMin = distSq
                              iobjBest(idz) = jobj
                              iptBest(idz) = jpt
                            endif
                          endif
                        endif
                      enddo
                    endif
                  enddo
                endif
                !
                if (iobjBest(idz) .ne. 0) then
                  iDzMin = min(iDzMin, idz)
                  iDzMax = max(iDzMax, idz)
                endif
              enddo
              !
              ! if didn't get the full range in both directions in Z,
              ! redo the search at farther Z values to find a full span
              ! of Z values that includes the point in question near its
              ! middle
              !
              if (loop == 1) then
                if (iDzMin == iDzLow .and. iDzMax == iDzHigh) then
                  ifSpan = 1
                else
                  if (iDzMax < iDzHigh) iDzLow = iDzMax - (numToFitInZ - 1)
                  if (iDzMin > -numToFitInZ / 2) iDzHigh = iDzMin + numToFitInZ - 1
                endif
              endif
              loop = loop + 1
            enddo
            !
            ! use next and previous point to define an angle to rotate to
            ! the horizontal
            !
            segLen = 0.
            do idir = 1, numInObj / 4
              ipLast = object(ibase + indMap(ipt - idir, numInObj))
              ipNext = object(ibase + indMap(ipt + idir, numInObj))
              dx = p_coord(1, ipNext) - p_coord(1, ipLast)
              dy = p_coord(2, ipNext) - p_coord(2, ipLast)
              segLen = sqrt(dx**2 + dy**2)
              if (segLen > 2) then
                exit
              endif
            enddo
            if (segLen > 2. .and. ((iDzMin < 0 .and. iDzMax > 0) .or. &
                iDzMax - iDzMin >= (3 * numToFitInZ) / 4)) then
              sinTheta = -dy / segLen
              cosTheta = dx / segLen
              numZfit = 0
              numFit = 0
              !
              ! for each Z level, start at the closest point and add
              ! points within the distance limit or until X starts to
              ! fold back toward the central point
              !
              do idz = iDzMin, iDzMax
                if (iobjBest(idz) .ne. 0) then
                  numZfit = numZfit + 1
                  jobj = iobjBest(idz)
                  jbase = ibase_obj(jobj)
                  xLast = -1.e10
                  jpt = iptBest(idz)
                  numInJ = npt_in_obj(jobj)
                  do idir = 1, -1, -2
                    ifTooFar = 0
                    do while(ifTooFar == 0 .and. numFit < idim)
                      jpnt = object(jbase + jpt)
                      dx = p_coord(1, jpnt) - xx
                      dy = p_coord(2, jpnt) - yy
                      xrot = cosTheta * dx - sinTheta * dy
                      yRot = sinTheta * dx + cosTheta * dy
                      if (idir * (xrot - xLast) < 0.) then
                        ifTooFar = 1
                      else
                        dist = sqrt(dx**2 + dy**2 + idz**2)
                        ifSave = 1
                        if (dist > distLimit) then
                          !
                          ! if go past the distance limit, add a point
                          ! within the limit if is far from the last
                          ! point
                          !
                          if (dist - distLast > 0.1 * distLimit) then
                            frac = (distLimit - distLast) / (dist - distLast)
                            xrot = xLast + frac * (xrot - xLast)
                            yRot = yLast + frac * (yRot - yLast)
                          else
                            ifSave = 0
                          endif
                          ifTooFar = 1
                        endif
                        if (ifSave == 1) then
                          numFit = numFit + 1
                          xtmp(numFit) = xrot
                          ytmp(numFit) = yRot
                          ztmp(numFit) = idz
                          if (jpt == iptBest(idz)) then
                            xcen = xrot
                            ycen = yRot
                            distCen = dist
                          else
                            !
                            ! also add points to maintain a maximum
                            ! separation between point = sepmin
                            !
                            sep = sqrt((xrot - xLast)**2 + (yRot - yLast)**2)
                            numInsert = sep / sepMin
                            do j = 1, numInsert
                              frac = j / (numInsert + 1.)
                              if (numFit < idim) then
                                numFit = numFit + 1
                                xtmp(numFit) = xLast + frac * (xrot - xLast)
                                ytmp(numFit) = yLast + frac * (yRot - yLast)
                                ztmp(numFit) = idz
                              endif
                            enddo
                          endif
                        endif
                        xLast = xrot
                        yLast = yRot
                        distLast = dist
                        jpt = indMap(jpt + idir, numInJ)
                      endif
                    enddo
                    jpt = indMap(iptBest(idz) - 1, numInJ)
                    xLast = xcen
                    yLast = ycen
                    distLast = distCen
                  enddo
                endif
              enddo
              !
              ! figure out a valid order for the fit and do the fit
              !
              norder = 0
              do while(norder < iorder .and. norder < numZfit - 1 .and. &
                  (norder + 1) * (norder + 4) <= numFit)
                norder = norder + 1
              enddo
              if (norder >= 0) then
                numFitTotal = numFitTotal + numFit
                numFitsDone = numFitsDone + 1
                numIndepVars = norder * (norder + 3) / 2  !# of independent variables
                do i = 1, numFit
                  call polytermReal(xtmp(i), ztmp(i), norder, xmat(1, i))
                  xmat(numIndepVars + 1, i) = ytmp(i)
                  ! if (modcont==1 .and. ipt == 128) print *,xt(i), yt(i), zt(i)
                enddo
                call multRegress(xmat, msiz, 1, numIndepVars, numFit, 1, 0, b1, msiz, &
                    c1, xm, sd, ssd)

                yNew = c1
                !
                ! back rotate the fitted point to get the new value
                !
                p_new(1, ipnt) = sinTheta * yNew + xx
                p_new(2, ipnt) = cosTheta * yNew + yy
              endif
              ! else
              ! print *,'skipping', iobj, ipt, dlen, idzmin, idzmax
            endif
          enddo
        endif
      enddo
      ! print *,'Fits at ', nDidFit, ' points, fit to', nFitTot
      !
      ! now treat each individual contour
      !
      do iobj = 1, max_mod_obj
        numInObj = npt_in_obj(iobj)
        if (numInObj > 4 .and. objectOK(iobj)) then
          call objToCont(iobj, obj_color, modelObj, modelCont)
          ! print *,'contour', modcont
          ibase = ibase_obj(iobj)
          !
          ! smoothing the same way as above: take each point as a center
          !
          do iptCen = 1, numInObj
            ipc = object(ibase + iptCen)
            xMiddle = p_new(1, ipc)
            yMiddle = p_new(2, ipc)
            if (iorder2 > 0) then
              segLen = 0.
              do idir = 1, numInObj / 4
                ipLast = object(ibase + indMap(iptCen - idir, numInObj))
                ipNext = object(ibase + indMap(iptCen + idir, numInObj))
                dx = p_new(1, ipNext) - p_new(1, ipLast)
                dy = p_new(2, ipNext) - p_new(2, ipLast)
                segLen = sqrt(dx**2 + dy**2)
                if (segLen > 2) then
                  exit
                endif
              enddo
              if (segLen > 2.) then
                sinTheta = -dy / segLen
                cosTheta = dx / segLen
                numFit = 0
                xLast = -1.e10
                jpt = iptCen
                !
                ! work from the center outward until get past limit or X
                ! folds back
                !
                do idir = 1, -1, -2
                  ifTooFar = 0
                  do while(ifTooFar == 0 .and. numFit < idim)
                    !
                    ! rotate so tangent is horizontal
                    !
                    jpnt = object(ibase + jpt)
                    dx = p_new(1, jpnt) - xMiddle
                    dy = p_new(2, jpnt) - yMiddle
                    xrot = cosTheta * dx - sinTheta * dy
                    yRot = sinTheta * dx + cosTheta * dy
                    if (idir * (xrot - xLast) < 0.) then
                      ifTooFar = 1
                    else
                      dist = sqrt(dx**2 + dy**2)
                      ifSave = 1
                      if (dist > distLimit) then
                        if (dist - distLast > 0.1 * distLimit) then
                          frac = (distLimit - distLast) / (dist - distLast)
                          xrot = xLast + frac * (xrot - xLast)
                          yRot = yLast + frac * (yRot - yLast)
                        else
                          ifSave = 0
                        endif
                        ifTooFar = 1
                      endif
                      if (ifSave == 1) then
                        numFit = numFit + 1
                        xtmp(numFit) = xrot
                        ytmp(numFit) = yRot
                        if (jpt == iptCen) then
                          xcen = xrot
                          ycen = yRot
                          distCen = dist
                        else
                          !
                          ! add points to maintain maximum separation:
                          ! this is superior to a fit to a fixed number
                          ! of points
                          sep = sqrt((xrot - xLast)**2 + (yRot - yLast)**2)
                          numInsert = sep / sepMin
                          do j = 1, numInsert
                            frac = j / (numInsert + 1.)
                            if (numFit < idim) then
                              numFit = numFit + 1
                              xtmp(numFit) = xLast + frac * (xrot - xLast)
                              ytmp(numFit) = yLast + frac * (yRot - yLast)
                            endif
                          enddo
                        endif
                      endif
                      xLast = xrot
                      yLast = yRot
                      distLast = dist
                      jpt = indMap(jpt + idir, numInObj)
                    endif
                  enddo
                  jpt = indMap(iptCen - 1, numInObj)
                  xLast = xcen
                  yLast = ycen
                  distLast = distCen
                enddo
                !
                ! get order, set up and do fit, substitute fitted point
                !
                norder = min(iorder2, numFit - 2)
                ! if (modcont==1) print *,iptcen, nfit, norder
                if (norder > 0) then
                  do i = 1, numFit
                    do j = 1, norder
                      xmat(j, i) = xtmp(i)**j
                    enddo
                    xmat(norder + 1, i) = ytmp(i)
                    ! if (modcont==1 .and. iptcen == 74) print *,(xr(j, i), j=1, norder), yt(i)
                  enddo
                  call multRegress(xmat, msiz, 1, norder, numFit, 1, 0, b1, msiz, &
                      bIntercept, xm, sd, ssd)
                  ! call polyfit(xt, yt, nfit, norder, slop, bint)
                  xMiddle = sinTheta * bIntercept + xMiddle
                  yMiddle = cosTheta * bIntercept + yMiddle
                endif
              endif
            endif
            p_coord(1, ipc) = xMiddle
            p_coord(2, ipc) = yMiddle
          enddo
          !
          ! eliminate close points
          !
          call elimClose(p_coord, 3, object, ibase, numInObj, closeThresh, 3)

          npt_in_obj(iobj) = numInObj
          !
          ! fix up points that are crossed
          !
          do ipt = 1, numInObj
            call uncrossCont(p_coord, 3, object, ibase, numInObj, ipt, 30., &
                tolCross)
          enddo
        endif
      enddo
      !
      ! reflip data if it was unflipped
      !
      if (ifFlipped .ne. 0) then
        do i = 1, n_point
          tmp = p_coord(2, i)
          p_coord(2, i) = p_coord(3, i)
          p_coord(3, i) = tmp
        enddo
      endif
      if (.not. retainMeshes) call deleteImodMeshes(imodObj)
      call scale_model(1)
      call putModelObjects()
    endif
  enddo
  n_point = -1
  call write_wmod(outputFile)
  if (retainMeshes) then
    print *,'DONE - Be sure to remesh the smoothed objects, especially before iterating'
  else
    print *,'DONE - Meshes have been deleted and the smoothed objects need remeshing'
  endif
  call exit(0)
end program smoothsurf

!
real*4 function goodAngle(thin)
  implicit none
  real*4 thin, theta
  theta = thin
  if (theta < -180) theta = theta + 360.
  if (theta >= 180.) theta = theta - 360.
  goodAngle = theta
  return
end function goodAngle



subroutine uncrossCont(p_coord, nxyz, object, ibase, numInObj, ip4, angle, &
    tol)
  implicit none
  integer*4 object(*), nxyz, ibase, numInObj, ip4
  real*4 p_coord(nxyz,*), angle, tol
  integer*4 ip3, ip1, ip2, ipnt1, ipnt2, ipnt3, ipnt4, itry
  real*4 x1, y1, x2, y2, x3, y3, x4, y4, tmin1, tmin2, distSq1, distSq2
  real*4 tStart, tEnd, tcon, halfDiff, tmid, halfCrit, fracLimStart
  real*4 fracLimEnd, fracStart, fracEnd, x1t, x4t, y1t, y4t
  integer*4 indMap
  real*4 atan2d, goodAngle
  !
  if (numInObj <= 3) return
  ip3 = indMap(ip4 - 1, numInObj)
  ip1 = indMap(ip4 + 1, numInObj)
  ip2 = indMap(ip4 + 2, numInObj)
  ipnt1 = object(ibase + ip1)
  ipnt2 = object(ibase + ip2)
  ipnt3 = object(ibase + ip3)
  ipnt4 = object(ibase + ip4)
  x1 = p_coord(1, ipnt1)
  y1 = p_coord(2, ipnt1)
  x2 = p_coord(1, ipnt2)
  y2 = p_coord(2, ipnt2)
  x3 = p_coord(1, ipnt3)
  y3 = p_coord(2, ipnt3)
  x4 = p_coord(1, ipnt4)
  y4 = p_coord(2, ipnt4)
  ! print *,sqrt((x1-x2)**2 +(y1-y2)**2), sqrt((x3-x4)**2 +(y3-y4)**2)
  call point_to_line(x4, y4, x1, y1, x2, y2, tmin1, distSq1)
  call point_to_line(x1, y1, x3, y3, x4, y4, tmin2, distSq2)
  if (min(distSq1, distSq2) <= tol**2) then
    tStart = atan2d(y2 - y1, x2 - x1)
    tEnd = atan2d(y4 - y3, x4 - x3)
    tcon = atan2d(y1 - y4, x1 - x4)
    halfDiff = 0.5 * goodAngle(tStart - tEnd)
    tmid = goodAngle(tEnd + halfDiff)
    halfCrit = abs(halfDiff) + angle
    if (abs(goodAngle(tcon - tmid)) > halfCrit) then
      fracLimStart = 1. -1. / sqrt((x2 - x1)**2 + (y2 - y1)**2)
      fracLimEnd = 1. -1. / sqrt((x4 - x3)**2 + (y4 - y3)**2)
      itry = 1
      ! print *,'shifting point', ipnt4, ' of', ninobj
      do while(itry <= 10)
        fracStart = max(0., itry * fracLimStart / 10.)
        fracEnd = max(0., itry * fracLimEnd / 10.)
        x1t = x1 + fracStart * (x2 - x1)
        x4t = x4 + fracEnd * (x3 - x4)
        y1t = y1 + fracStart * (y2 - y1)
        y4t = y4 + fracEnd * (y3 - y4)
        tcon = atan2d(y1t - y4t, x1t - x4t)
        if (abs(goodAngle(tcon - tmid)) <= halfCrit) itry = 10
        itry = itry + 1
      enddo
      p_coord(1, ipnt1) = x1t
      p_coord(2, ipnt1) = y1t
      p_coord(1, ipnt4) = x4t
      p_coord(2, ipnt4) = y4t
    endif
  endif
  return
end subroutine uncrossCont




subroutine point_to_line(x0, y0, x1, y1, x2, y2, tmin, distSq)
  implicit none
  real*4 x0, y0, x1, y1, x2, y2, tmin, distSq
  tmin = ((x0 - x1) * (x2 - x1) + (y0 - y1) * (y2 - y1)) / ((x2 - x1)**2 + (y2 - y1)**2)
  tmin = max(0., min(1., tmin))
  distSq = (x1 + tmin * (x2 - x1) - x0)**2 + (y1 + tmin * (y2 - y1) - y0)**2
  return
end subroutine point_to_line



subroutine elimClose(p_new, nxyz, object, ibase, numInObj, closeThresh, &
    minPoints)
  implicit none
  integer*4 object(*), nxyz, ibase, numInObj, minPoints
  real*4 p_new(nxyz,*), closeThresh
  real*4 closesq
  integer*4 iseg, idel, ip3, ipt3, ipSeg, ip1, ip4, ipt1, ipt4, i
  integer*4 indMap, numIn
  !
  closesq = closeThresh**2
  iseg = 1
  idel = 0
  numIn = numInObj
  do while(iseg <= numInObj .and. numInObj > minPoints)
    ip3 = indMap(iseg + 1, numInObj)
    ipt3 = object(ibase + ip3)
    ipSeg = object(ibase + iseg)
    if ((p_new(1, ipt3) - p_new(1, ipSeg))**2 + &
        (p_new(2, ipt3) - p_new(2, ipSeg))**2 < closesq) then
      ! print *,'eliminating segment', iseg, ' of', ninobj
      ip1 = indMap(iseg - 1, numInObj)
      ip4 = indMap(iseg + 2, numInObj)
      ipt1 = object(ibase + ip1)
      ipt4 = object(ibase + ip4)
      idel = iseg
      if ((p_new(1, ipt1) - p_new(1, ipSeg))**2 + &
          (p_new(2, ipt1) - p_new(2, ipSeg))**2 > &
          (p_new(1, ipt3) - p_new(1, ipt4))**2 + &
          (p_new(2, ipt3) - p_new(2, ipt4))**2) idel = ip3
      do i = ibase + idel + 1, ibase + numInObj
        object(i - 1) = object(i)
      enddo
      numInObj = numInObj - 1
    else
      iseg = iseg + 1
    endif
  enddo
  return
end subroutine elimClose
