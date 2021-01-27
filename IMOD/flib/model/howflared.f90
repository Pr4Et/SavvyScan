! HOWFLARED
!
! Computes various measures of microtubule flaring and curvature from a
! model with tracings of the microtubule walls.
! See man page for details.
!
! Written by David Mastronarde, 1996
!
! $Id$
!
program howflared
  implicit none
  include 'smallmodel.inc90'
  integer LIMDATA, LIMCOL
  parameter (LIMDATA = 1000, LIMCOL = 30)
  logical readSmallMod, exist, started
  character*120 modelFile, fileOut, colString, curveOut
  integer*4 indObjAtZ(99), numWalls(2), indPackStart(99)
  real*4 xObj(99), xWall(LIMDATA,2), yWall(LIMDATA,2), xxTemp(LIMDATA), yyTemp(LIMDATA)
  real*4 xxLR(LIMDATA), areaSum(3), arSumSqrt(3), fitYmin(2), angleSum(3)
  real*4 width(LIMDATA), col(LIMDATA,LIMCOL), colTemp(LIMCOL), finalAng(2), totalLen(2)
  integer*4 izSave(2,LIMDATA), icolOut(LIMDATA), iobjSave(LIMDATA), isurfSave(LIMDATA)
  integer*4 itimes(max_obj_num), isurfs(max_obj_num), npTimes(LIMDATA)
  integer*4 icontSave(LIMDATA)
  logical objIsUsed(max_obj_num)
  real*4 rmat(3,3), packXYZ(3,LIMDATA), packedRot(3,LIMDATA), straightEnd(3,2,LIMDATA)
  real*4 xCurve(LIMDATA), yCurve(LIMDATA), vec(3), endMean(3)
  integer*4 numColOut, ident, ierr, iobj, jobj, imodObj, imodCont, numWidths, iModel
  integer*4 jpt, izObj, numAtZ, markerObj, numPacked, ipt, markerPackInd, itemp, lrInd
  integer*4 numData, indOpposite, ibase, ifFlip, numModels, numStartEnds, numIDs, idStart
  integer*4 iobjLim, numCurve, i, j, izTime, numVec, numSDs, maxNPtimes
  real*4 fitTop, fitBot, fitYmax, yMax, yMin, ytmp, xIntrp, polarity, xx, yy, xLast
  real*4 baseLast, yLast, xFit, deltaX, deltaXwidth, deltaY, avgWidth, sdWidth, sem
  real*4 avgCol, sdCol, semCol, base, slope, xLeft, deltaAng, segmentLen, pixSize, dxLast
  real*4 pixSizeDef, xyScale, zScale, xOffset, yOffset, zOffset, widthDef, defScale
  real*4 fitStartDef, fitEndDef, cosAngle, sinAngle, centroid(3), dyLast
  real*4 dotLength, dotAvg, dotSD, dotSum, dotSumSq, sdSum, sdSumSq
  logical*4 useSurf, noPairs
  integer*4 getImodTimes, getImodHead, getImodSurfaces
  real*4 xinterp, acosd, atan2d, cosd, sind, atand
  !
  logical pipInput
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetInteger, PipGetBoolean, PipGetTwoFloats
  integer*4 PipGetString, PipGetFloat, PipGetLogical
  integer*4 PipGetInOutFile, PipNumberOfEntries
  !
  ! fallbacks from ../../manpages/autodoc2man -2 2  howflared
  !
  integer numOptions
  parameter (numOptions = 11)
  character*(40 * numOptions) options(1)
  options(1) = &
      'output:OutputFile:FN:@columns:ColumnsToOutput:LI:@'// &
      'point:PointOutputFile:FN:@pixel:PixelSizeDefault:F:@'// &
      'width:WidthDefault:F:@surface:UseSurfaceNumbers:B:@'// &
      'nopairs:NoPairsOrMarkers:B:@model:ModelFile:FNM:@'// &
      'fit:FitTopAndBottom:FPM:@id:Identifier:IM:@help:usage:B:'
  !
  pixSizeDef = 1.
  widthDef = 20.
  idStart = 1
  fitStartDef = 0.
  fitEndDef = 0.
  useSurf = .false.
  noPairs = .false.
  curveOut = ' '

  print *,'Columns 1-6 are for LEFT, 7-12 for RIGHT, ', &
      '13-18 for SUM, 19-22 LEFT, 23-26 RIGHT'
  print *,'The pattern for 1-18 is: 1 area, 2 square root area, ', &
      '3 angular sum (radians),', ' 4 width-normalized', &
      ' area, 5 normalized root, 6 normalized angular sum'
  print *,'19/23 is total length, 20/24 is final angle (degrees),', &
      '21/25 is angle change per unit length, ', &
      '22/26 is average radius of curvature'
  !
  call PipReadOrParseOptions(options, numOptions, 'howflared', &
      'ERROR: HOWFLARED - ', .false., 3, 0, 1, numOptArg, &
      numNonOptArg)
  if (PipGetInOutFile('OutputFile', 1, ' ', fileOut) .ne. 0) &
      call exitError('NO OUTPUT FILE SPECIFIED')

  if (PipGetString('ColumnsToOutput', colString) .ne. 0) &
      call exitError('NO LIST OF COLUMNS TO OUTPUT ENTERED')

  call parseList(colString, icolOut, numColOut)
  do i = 1, numColOut
    if (icolOut(i) < 0 .or. icolOut(i) > 26) call exitError( &
        'ILLEGAL COLUMN NUMBER ENTERED')
  enddo

  ierr = PipGetFloat('PixelSizeDefault', pixSizeDef)
  ierr = PipGetFloat('WidthDefault', widthDef)
  ierr = PipGetLogical('UseSurfaceNumbers', useSurf)
  ierr = PipGetLogical('NoPairsOrMarkers', noPairs)
  if (useSurf .and. noPairs) call exitError( &
      'YOU CANNOT ENTER BOTH -surface AND -nopairs')

  ierr = PipNumberOfEntries('ModelFile', numModels)
  if (ierr .ne. 0 .or. numModels == 0) &
      call exitError('NO INPUT MODELS ENTERED')

  ierr = PipNumberOfEntries('FitTopAndBottom', numStartEnds)
  if (numStartEnds > 1 .and. numStartEnds .ne. numModels) call exitError &
      ('NUMBER OF ENTRIES FOR START AND END VALUES MUST EQUAL 0, 1, OR'// &
      ' THE NUMBER OF MODELS')

  ierr = PipNumberOfEntries('Identifier', numIDs)
  if (numIDs > 1 .and. numIDs .ne. numModels) call exitError &
      ('NUMBER OF IDENTIFIER VALUES MUST EQUAL 0, 1, OR'// &
      ' THE NUMBER OF MODELS')

  !
  ! open the output file
  !
  call dopen(1, fileOut, 'new', 'f')
  i = numColOut + 2
  if (useSurf) i = i + 1
  write(1, '(i4)') i
  !
  ! Open file for curve output if present
  !
  ierr = PipGetString('PointOutputFile', curveOut)
  if (curveOut .ne. ' ') call dopen(2, curveOut, 'new', 'f')
  if (curveOut .ne. '' .and. noPairs) call dopen(3, 'starts.'//curveout, 'new', 'f')
  !
  ! loop on the models and get the model name and optional fit start and
  ! end and identifier for each model
  !
  sdSum = 0.
  sdSumSq = 0.
  numSDs = 0
  do iModel = 1, numModels
    ierr = PipGetString('ModelFile', modelFile)
    print *,'Processing model: ', trim(modelFile)

    if (.not.readSmallMod(modelFile)) call exitError('READING MODEL FILE')
    !
    if (iModel <= numStartEnds) ierr = PipGetTwoFloats('FitTopAndBottom', &
        fitStartDef, fitEndDef)
    fitTop = fitStartDef
    fitBot = fitEndDef

    if (iModel <= numIDs) ierr = PipGetInteger('Identifier', idStart)
    ident = idStart
    idStart = idStart + 1
    !
    pixSize = pixSizeDef
    defScale = 1.e6
    ierr = getImodHead(xyScale, zScale, xOffset, yOffset, zOffset, ifFlip)
    if (ierr == 0 .and. abs(xyScale - defScale) / defScale > 1.e-5) then
      pixSize = 1000. * xyScale
      print *,'Pixel size set from model header to', pixSize
    else
      print *,'Using default pixel size of', pixSize
    endif
    !
    ! look at each object and try to find match, and find the left and
    ! right walls and the connector
    !
    ierr = getImodTimes(itimes)
    ierr = getImodSurfaces(isurfs)
    do iobj = 1, max_mod_obj
      objIsUsed(iobj) = .false.
    enddo
    numWidths = 0
    maxNPtimes = -1
    call scale_model(0)

    do jobj = 1, max_mod_obj
      if (npt_in_obj(jobj) > 1 .and. .not.objIsUsed(jobj)) then
        jpt = object(ibase_obj(jobj) + 1)
        izObj = nint(p_coord(3, jpt))
        !
        ! have a candidate for matching, now scan through objects
        ! starting with this one
        !
        numAtZ = 0
        markerObj = 0
        numPacked = 0
        izTime = izObj + 1
        if (itimes(jobj) > 0) izTime = itimes(jobj)
        iobjLim = max_mod_obj
        if (noPairs) iobjLim = jobj
        do iobj = jobj, iobjLim
          if (npt_in_obj(iobj) > 1 .and. .not.objIsUsed(iobj)) then
            ipt = object(ibase_obj(iobj) + 1)
            !
            ! object matches if it is same imod object, same time, and
            ! either at an actual time or at same Z
            !
            if (obj_color(2, iobj) == obj_color(2, jobj) .and. &
                itimes(iobj) == itimes(jobj) .and. &
                (.not.useSurf .or. isurfs(iobj) == isurfs(jobj)) .and. &
                (itimes(iobj) > 0. .or. &
                nint(p_coord(3, ipt)) == izObj)) then
              !
              ! store data about connector or about wall line
              !
              if (npt_in_obj(iobj) == 2) then
                markerObj = iobj
                markerPackInd = numPacked + 1
              else
                numAtZ = min(99, numAtZ + 1)
                indObjAtZ(numAtZ) = iobj
                xObj(numAtZ) = p_coord(1, ipt)
                indPackStart(numAtZ) = numPacked + 1
              endif
              objIsUsed(iobj) = .true.
              if (numPacked + npt_in_obj(iobj) > LIMDATA) &
                  call exitError('TOO MANY POINTS FOR ARRAYS')
              do i = 1, npt_in_obj(iobj)
                ipt = object(ibase_obj(iobj) + i)
                do j = 1, 2
                  packXYZ(j, i + numPacked) = p_coord(j, ipt)
                enddo
                packXYZ(3, i + numPacked) = zScale * p_coord(3, ipt)
              enddo
              numPacked = numPacked + npt_in_obj(iobj)
            endif
          endif
        enddo

        call objToCont(jobj, obj_color, imodObj, imodCont)
        if (numAtZ > 2) then
          write(*,'(/,a,i4,a,i4,a,i4,a,i4)') 'ERROR: HOWFLARED -', &
              numAtZ - 1, ' contours match up with object', imodObj, &
              ', contour', imodCont, ', Z/time', izTime
          call exit(1)
        endif
        !
        ! got one or two: fit to a plane and get rotated points
        call planeFit(packXYZ, numPacked, centroid, rmat, packedRot)
        ! do i=1, numpacked
        ! do j=1, 3
        ! packrot(j, i) =packxyz(j, i)
        ! enddo
        ! print *,(packxyz(j, i), j=1, 3), (packedrot(j, i), j=1, 3)
        ! enddo

        ! switch so left is first, copy into wall arrays
        !
        if (numAtZ == 2 .and. xObj(2) < xObj(1)) then
          itemp = indObjAtZ(1)
          indObjAtZ(1) = indObjAtZ(2)
          indObjAtZ(2) = itemp
          itemp = indPackStart(1)
          indPackStart(1) = indPackStart(2)
          indPackStart(2) = itemp
        endif
        yMin = 1.e10
        yMax = -1.e10
        do lrInd = 1, numAtZ
          iobj = indObjAtZ(lrInd)
          numWalls(lrInd) = npt_in_obj(iobj)
          do i = 1, numWalls(lrInd)
            ipt = indPackStart(lrInd) + i - 1
            xWall(i, lrInd) = packedRot(1, ipt)
            yWall(i, lrInd) = packedRot(2, ipt)
            yMin = min(yMin, yWall(i, lrInd))
            yMax = max(yMax, yWall(i, lrInd))
          enddo
        enddo
        !
        ! compute top and bottom limits to fit
        !
        if (fitTop == 0.) then
          fitYmax = yMax
        elseif (fitTop < 1.) then
          fitYmax = yMax - fitTop * (yMax - yMin)
        else
          fitYmax = yMax - fitTop
        endif
        do lrInd = 1, numAtZ
          !
          ! bottom = 0 means use the second point or the marker,
          ! otherwise use a fractional distance or an absolute distance
          !
          if (fitBot == 0.) then
            fitYmin(lrInd) = yWall(2, lrInd)
          elseif (fitBot < 1.) then
            fitYmin(lrInd) = yMax - fitBot * (yMax - yMin)
          else
            fitYmin(lrInd) = yMax - fitBot
          endif
          if (markerObj .ne. 0) then
            i = markerPackInd + 2 - lrInd
            if (packedRot(1, markerPackInd) < packedRot(1, markerPackInd + 1)) &
                i = markerPackInd + lrInd - 1
            if (fitBot > 0) then
              fitYmin(lrInd) = max(fitYmin(lrInd), packedRot(2, i))
            else
              fitYmin(lrInd) = packedRot(2, i)
            endif
          endif
        enddo
        !
        numData = 0
        ! print *,'fitmax & min', fitmax, (fitmin(lr), lr =1 , natz)
        do lrInd = 1, numAtZ
          !
          ! do fit for one wall: use all points from limit of this wall,
          ! interpolate at lower limit and maybe at upper limit
          !
          do i = 1, numWalls(lrInd)
            ytmp = yWall(i, lrInd)
            if (ytmp >= fitYmin(lrInd) .and. ytmp <= fitYmax) then
              numData = numData + 1
              xxTemp(numData) = xWall(i, lrInd)
              yyTemp(numData) = ytmp
              xxLR(numData) = lrInd - 1
            endif
          enddo
          ! print *,lr, ' native', ndat
          ytmp = fitYmin(lrInd)
          xIntrp = xinterp(xWall(1, lrInd), yWall(1, lrInd), numWalls(lrInd), ytmp)
          if (xIntrp .ne. -99999. .and. abs(ytmp - yyTemp(numData)) > 1.e-3) then
            numData = numData + 1
            xxTemp(numData) = xIntrp
            yyTemp(numData) = ytmp
            xxLR(numData) = lrInd - 1
          endif
          ! print *,lr, ' bottom', ndat
          if (fitTop .ne. 0.) then
            ytmp = fitYmax
            xIntrp = xinterp(xWall(1, lrInd), yWall(1, lrInd), numWalls(lrInd), ytmp)
            if (xIntrp .ne. -99999.) then
              numData = numData + 1
              xxTemp(numData) = xIntrp
              yyTemp(numData) = ytmp
              xxLR(numData) = lrInd - 1
            endif
          endif
          ! print *,lr, ' top', ndat
          !
          ! Then add points interpolated from Y positions of other wall
          !
          if (numAtZ == 2) then
            indOpposite = 3 - lrInd
            do j = 1, numWalls(indOpposite)
              ytmp = yWall(j, indOpposite)
              if (ytmp >= fitYmin(indOpposite) .and. ytmp <= fitYmax) then
                xIntrp = xinterp(xWall(1, lrInd), yWall(1, lrInd), numWalls(lrInd), ytmp)
                if (xIntrp .ne. -99999.) then
                  numData = numData + 1
                  xxTemp(numData) = xIntrp
                  yyTemp(numData) = ytmp
                  xxLR(numData) = lrInd - 1
                endif
              endif
            enddo
          endif
          ! print *,lr, ' other', ndat
        enddo
        !
        ! fit points from both walls to a single 3-parameter fit
        !
        if (numData == 0) then
          write(*,'(/,a,i4,a,i4,a,i4)') &
              'ERROR: HOWFLARED - no points to fit to for obj', &
              imodObj, ' surf', isurfs(jobj), ' cont', imodCont
          call exit(1)
        endif
        if (numAtZ == 2) then
          call lsfit2(yyTemp, xxLR, xxTemp, numData, slope, deltaXwidth, xLeft)
        else
          call lsfit(yyTemp, xxTemp, numData, slope, xLeft, deltaXwidth)
          deltaXwidth = widthDef / pixSize
        endif
        ! print *,'fit', ndat, slope, delxw, xleft, ywall(1, 1)
        numWidths = numWidths + 1
        width(numWidths) = deltaXwidth * pixSize
        do i = 1, 26
          col(numWidths, i) = 0.
        enddo
        areaSum(2) = 0.
        angleSum(2) = 0.
        do lrInd = 1, numAtZ
          polarity = 2 * lrInd - 3
          areaSum(lrInd) = 0.
          angleSum(lrInd) = 0.
          totalLen(lrInd) = 0.
          finalAng(lrInd) = 0.
          dxLast = -slope
          dyLast = -1.
          !
          ! integrate either from bottom of fit or from the zero-marker
          !
          if (markerObj == 0) then
            yLast = fitYmin(lrInd)
            xFit = (lrInd - 1) * deltaXwidth + xLeft
            baseLast = slope * yLast + xFit
          else
            if (packedRot(1, markerPackInd) < packedRot(1, markerPackInd + 1)) then
              yLast = packedRot(2, markerPackInd + lrInd - 1)
            else
              yLast = packedRot(2, markerPackInd + 2 - lrInd)
            endif
            baseLast = xinterp(xWall(1, lrInd), yWall(1, lrInd), numWalls(lrInd), yLast)
          endif
          xLast = xinterp(xWall(1, lrInd), yWall(1, lrInd), numWalls(lrInd), yLast)
          numCurve = 1
          xCurve(1) = baseLast
          yCurve(1) = yLast
          started = .false.
          if (xLast .ne. -99999.) then
            do i = 1, numWalls(lrInd)
              yy = yWall(i, lrInd)
              xx = xWall(i, lrInd)
              if (yy < yLast .or. started) then
                !
                ! first time, if not a paired item, figure out polarity
                !
                if (.not.started .and. numAtZ == 1) &
                    polarity = sign(1., xWall(numWalls(1), 1) - xx)
                !
                ! sum area between segment and baseline if still going
                ! down
                !
                base = baseLast + slope * (yy - yLast)
                if (yy < yLast) areaSum(lrInd) = areaSum(lrInd) + &
                    polarity * (yLast - yy) * (xx - base + xLast - baseLast) / 2.
                baseLast = base
                !
                ! sum product of angular deviation and length
                !
                deltaX = xx - xLast
                deltaY = yy - yLast
                segmentLen = sqrt(deltaX**2 + deltaY**2)
                totalLen(lrInd) = totalLen(lrInd) + segmentLen
                if (segmentLen > 1.e-6) then
                  deltaAng = atan2d(deltaY, deltaX) - atan2d(dyLast, dxLast)
                  if (deltaAng < -180.) deltaAng = deltaAng + 360.
                  if (deltaAng > 180.) deltaAng = deltaAng - 360.
                  finalAng(lrInd) = finalAng(lrInd) + polarity * deltaAng
                  angleSum(lrInd) = angleSum(lrInd) + segmentLen * finalAng(lrInd) * &
                      3.14159 / 180.
                endif
                ! print *,itimes(jobj), lr, pol, i, delx, dely, seglen, delang, finalang(lr)
                xLast = xx
                yLast = yy
                dxLast = deltaX
                dyLast = deltaY
                started = .true.
                numCurve = numCurve + 1
                xCurve(numCurve) = xx
                yCurve(numCurve) = yy
              endif
            enddo
            !
            ! Output points for this wall if more than one
            !
            if (curveOut .ne. ' ' .and. numCurve > 1) then
              deltaAng = atand(slope)
              cosAngle = pixSize * cosd(deltaAng)
              sinAngle = pixSize * sind(deltaAng)
              if (noPairs) then
                write(2, 111) numCurve, ident, imodObj, imodCont, izTime, lrInd
              else if (useSurf) then
                write(2, 111) numCurve, ident, imodObj, isurfs(jobj), izTime, lrInd
              else
                write(2, 111) numCurve, ident, imodObj, izTime, lrInd
              endif
111           format(6i5)
              do i = 1, numCurve
                deltaX = xCurve(i) - xCurve(1)
                deltaY = yCurve(i) - yCurve(1)
                xx = cosAngle * deltaX - sinAngle * deltaY
                yy = sinAngle * deltaX + cosAngle * deltaY
                write(2, '(2f9.2)') xx, yy
              enddo
            endif
          endif
        enddo
        areaSum(3) = areaSum(1) + areaSum(2)
        angleSum(3) = angleSum(1) + angleSum(2)
        do i = 1, 3
          arSumSqrt(i) = sign(sqrt(abs(areaSum(i))), areaSum(i))
          ibase = 6 * i - 5
          col(numWidths, ibase) = areaSum(i) * pixSize**2
          col(numWidths, ibase + 1) = arSumSqrt(i) * pixSize
          col(numWidths, ibase + 2) = angleSum(i) * pixSize
        enddo
        do i = 1, numAtZ
          ibase = 19 + (i - 1) * 4
          col(numWidths, ibase) = totalLen(i) * pixSize
          col(numWidths, ibase + 1) = finalAng(i)
          if (totalLen(i) > 0.01) col(numWidths, ibase + 2) = &
              finalAng(i) / (totalLen(i) * pixSize)
          if (abs(finalAng(i)) > 0.1) col(numWidths, ibase + 3) = &
              (totalLen(i) * pixSize * 360.) / (2. *3.14159 * finalAng(i))
        enddo
        izSave(1, numWidths) = izTime
        !
        ! Save "no pair" times separately, in case there is only one file we don't want iz
        npTimes(numWidths) = itimes(jobj)
        maxNPtimes = max(maxNPtimes, itimes(jobj))

        iobjSave(numWidths) = imodObj
        isurfSave(numWidths) = isurfs(jobj)
        icontSave(numWidths) = imodCont
        if (noPairs) then
          ipt = indPackStart(1)
          do i = 1, 3
            straightEnd(i, 1, numWidths) = packXYZ(i, ipt)
            straightEnd(i, 2, numWidths) = packXYZ(i, ipt + 1)
          enddo
        endif
        ! write(1, 101) ident, izsav(1, nwidth) &
        ! , (arsum(i), arsqrt(i), angsum(i), i=1, 3)
        ! 101           format(2i5,(6f10.2))
      endif
    enddo
    !
    call avgsd(width, numWidths, avgWidth, sdWidth, sem)
    write(*,107) avgWidth, sdWidth, numWidths
107 format('Width mean =',f8.2,', s.d. =',f8.2,', n =',i4)
    do j = 1, numWidths
      do i = 1, 3
        ibase = 6 * i - 5
        col(j, ibase + 3) = col(j, ibase) / avgWidth**2
        col(j, ibase + 4) = col(j, ibase + 1) / avgWidth
        col(j, ibase + 5) = col(j, ibase + 2) / avgWidth
      enddo
      do i = 1, numColOut
        colTemp(i) = col(j, icolOut(i))
      enddo
      if (useSurf .or. noPairs) then
        write(1, 102) ident, iobjSave(j), isurfSave(j), izSave(1, j), &
            (colTemp(i), i = 1, numColOut)
102     format(4i5,(14f10.3))
      else
        write(1, 103) ident, iobjSave(j), izSave(1, j), &
            (colTemp(i), i = 1, numColOut)
103     format(3i5,(14f10.3))
      endif
    enddo
    do i = 1, 26
      call avgsd(col(1, i), numWidths, avgCol, sdCol, semCol)
      write(*,108) i, avgCol, sdCol, numWidths
108   format(i3,2f10.3,i4)
    enddo
    !
    ! If there is curve output and no pairs, automatically put out special "starts" file
    if (curveOut .ne. ' ' .and. noPairs) then
      do izTime = 0, maxNPtimes
        vec(1:3) = 0.
        endMean(1:3) = 0.
        numVec = 0
        !
        ! For each time, find all the protofilaments and get the average vector of the
        ! the first two points and the average position of the endpoint
        do  j = 1, numWidths
          if (npTimes(j) == izTime) then
            vec(1:3) = vec(1:3) + (straightEnd(1:3, 2, j) - straightEnd(1:3, 1, j))
            endMean(1:3) = endMean(1:3) + straightEnd(1:3, 1, j)
            numVec = numVec + 1
          endif
        enddo
        if (numVec > 1) then
          vec(1:3) = vec(1:3) / sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
          endMean(1:3) = endMean(1:3) / numVec
          dotSum = 0.
          dotSumSq = 0.
          !
          ! The projection length of the vector from the average start point to the 
          ! second point for each PF along the average vector is the dot product of that
          ! vector and the average vector
          do  j = 1, numWidths
            if (npTimes(j) == izTime) then
              dotLength = vec(1) * (straightEnd(1, 2, j) - endMean(1)) + &
                  vec(2) * (straightEnd(2, 2, j) - endMean(2)) + &
                  vec(3) * (straightEnd(3, 2, j) - endMean(3))
              write(3, '(i4,i5,i6,f10.3)') imodel, izTime, icontSave(j), dotLength
              dotSum = dotSum + dotLength
              dotSumSq = dotSumSq + dotLength**2
            endif
          enddo
          call sums_to_avgsd(dotSum, dotSumSq, numVec, dotAvg, dotSD, ytmp)
          write(3, '(i4,i5,a,f10.3)') imodel, izTime, '  SD  ', dotSD
          numSDs = numSDs + numVec
          sdSum = sdSum + numVec * dotSD
          sdSumSq = sdSumSq + numVec * dotSD**2
        endif
      enddo
    endif
  enddo
  close(1)
  if (curveOut .ne. ' ') close(2)
  if (curveOut .ne. ' ' .and. noPairs) then
    if (numSDs > 0 .and. numSDs > numVec) then
      call sums_to_avgsd(sdSum, sdSumSq, numSDs, dotAvg, dotSD, ytmp)
      write(3, '(a,i0,a,2f10.3)') 'Weighted mean and SD of SDs from', numSDs, &
          ' protofilaments', dotAvg, dotSD
    endif
    close(3)
  endif
  call exit(0)
end program howflared


real*4 function xinterp(x,y,nn,yIntcp)
  implicit none
  real*4 x(*), y(*), yIntcp
  integer*4 nn, i

  xinterp = -99999.
  do i = 1, nn - 1
    if (yIntcp >= min(y(i + 1), y(i)) .and. yIntcp <= max(y(i + 1), y(i))) then
      xinterp = x(i) + (x(i + 1) - x(i)) * (yIntcp - y(i)) / (y(i + 1) - y(i))
      return
    endif
  enddo
  return
end function xinterp
