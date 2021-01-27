! -----------------------------------------------------------------
!
! TILT......A program for reconstructing a three-dimensional object
! from a series of two-dimensional projections.
! The projections are assumed to arise from rotation about
! a fixed tilt axis.
!
! See man page for details
!
! $Id$
!
program tilt
  use tiltvars
  implicit none
  integer*4 nxyzTmp(3), nxyzst(3)
  data nxyzst/0., 0., 0./
  character*20 radtxt1/'Radial weighting'/
  character*18 radtxt2/'   function'/
  !
  integer*4 intervalHeadSave, numSliceOut, numSlices, nxprj2
  integer*4 inLoadStart, inLoadEnd, lastReady, lastCalc, nextFreeVertSlice
  integer*4 lvertSliceStart, lvertSliceEnd, numVertSliceInRing, ni, loadLimit, lsliceOut
  integer*4 lsliceStart, lsliceEnd, lslice, needStart, needEnd, itryEnd, ibaseSIRT
  integer*4 itry, ifEnough, lastStart, lastEnd, i, lsliceMin, lsliceMax, nextReadFree
  integer*4 iv, numAlready, lsliceProjEnd, lReadStart, lReadEnd, numReadInRing
  real*4 dmin, dmax, ycenFix, absSinAlf, tanAlpha, dmin4, dmax4, dmin5, dmax5, dmin6
  real*4 valMin, vertSliceCen, vertYcenFix, riBot, riTop, tmax, tmean, tmin, dmax6
  integer*4 lriMin, lriMax, loadStart, loadEnd, lri, isirtIter, ierr, j, k
  integer*4 ibase, lstart, nv, istart, nl, iyload, ix, ioffset
  integer*4 iringStart, mode, needGpuStart, needGpuEnd, keepOnGpu, numLoadGpu
  real*4 unscaledMin, unscaledMax, recScale, recAdd, dmean, pixelTot
  real*4 reprojFill, curMean, firstMean, vertSum, numVertSum, edgeFillOrig
  real*4 composeFill
  integer*4 iset, mapEnd, nz5, lfillStart, lfillEnd
  real*8 dtot8
  logical*4 shiftedGpuLoad, composedOne, truncations, extremes
  integer*4 gpuLoadProj, gpuReprojOneSlice, parWrtSetCurrent, iiuParWrtFlushBuffers
  real*8 wallTime, timeStart, dsum, dpix
  !
  intervalHeadSave = 20
  numSliceOut = 0
  dtot8 = 0.
  dmin = 1.e30
  dmax = -1.e30
  dmin4 = 1.e30
  dmax4 = -1.e30
  dmin5 = 1.e30
  dmax5 = -1.e30
  dmin6 = 1.e30
  dmax6 = -1.e30
  nz5 = 0
  debug = .false.
  !
  ! Open files and read control data
  call inputParameters
  numSlices = (isliceEnd - isliceStart) / idelSlice + 1
  valMin = 1.e-3 * (dmaxIn - dminIn)
  edgeFill = dmeanIn
  if (ifLog .ne. 0) edgeFill = alog10(max(valMin, dmeanIn + baseForLog))
  edgeFill = edgeFill * zeroWeight
  edgeFillOrig = edgeFill
  mapEnd = indOutSlice + isliceSizeBP - 1

  if (debug)  print *,'iflog=', ifLog, ' scale=', outScale, '  edgefill=', edgeFill
  composeFill = edgeFill * numViews
  !
  ! recompute items not in common
  nxprj2 = nxProj + 2 + numPad
  !
  ! initialize variables for loaded slices
  inLoadStart = 0
  inLoadEnd = 0
  lastReady = 0
  lastCalc = 0
  !
  ! initialize variables for ring buffer of vertical slices: # in the ring, next free
  ! position, starting and ending slice number of slices in ring
  numVertSliceInRing = 0
  nextFreeVertSlice = 1
  lvertSliceStart = -1
  lvertSliceEnd = -1
  !
  ! initialize similar variables for ring buffer of read-in slices
  numReadInRing = 0
  nextReadFree = 1
  lReadStart = -1
  lReadEnd = -1
  lfillStart = -1
  ycenFix = ycenOut
  absSinAlf = abs(sinAlpha(1))
  composedOne = ifAlpha >= 0
  numVertSum = 0
  vertSum = 0.
  !
  ! Calculate and report stack loading
  ! The stack is constructed as follows:
  ! (Radial weighting function, 1 or two planes) --->
  ! (current output plane) --->
  ! (vertical reconstructed planes for new X-axis tilting) --->
  ! (plane(s) read in from rec for SIRT) --->
  ! (working reprojection plane for SIRT)
  ! (first slice, first view) ....(first slice, last view) --->
  ! ..................................................--->
  ! (last slice, first view ) ....(last slice, last view)
  ! (Extra space for another set of slices, if cosine stretching)
  ! After each view two additional elements are inserted for
  ! the transform.
  !
  ni = numPlanes * inPlaneSize
  if (.not. recReproj) then
    write(6, 900) maxStack, radtxt1, indOutSlice - 1, radtxt2, isliceSizeBP
    if (ifAlpha < 0) write(6, 902) numVertNeeded, numVertNeeded * ithickBP * iwidth
    if (numSIRTiter > 0) then
      if (indWorkPlane - ireadBase > 0) &
          write(6, 904) max(1, numReadNeed), indWorkPlane - ireadBase
      write(6, 907) inPlaneSize
    endif
    write(6, 903) numPlanes, ni
    if (ipExtraSize .ne. 0) write(6, 901) ipExtraSize
    !
    ! Allocate maskEdges array and prepare fixed maskEdges
    allocate(ixUnmaskedSE(2, ithickBP), stat = ierr)
    call memoryErrorUC(ierr, 'MASK ARRAY')
    if (ifAlpha == 0 .or. .not.maskEdges) call maskPrep(isliceStart)
  endif
  if (debug) print *,'slicen', centerSlice, ', imap', indOutSlice, ', nbase', indLoadBase
  if (ifAlpha >= 0) then
    loadLimit = isliceEnd
  else
    loadLimit = centerSlice + (isliceEnd - centerSlice) * cosAlpha(1) +  &
        yOffset * sinAlpha(1) + 0.5 * ithickOut * absSinAlf + 2.
  endif
  !
  ! Main loop over slices perpendicular to tilt axis
  ! ------------------------------------------------
  lsliceOut = isliceStart
  do while (lsliceOut <= isliceEnd)
    !
    ! get limits for slices that are needed: the slice itself for regular
    ! work, or required vertical slices for new-style X tilting
    !
    if (debug) print *,'working on', lsliceOut
    if (ifAlpha >= 0) then
      lsliceStart = lsliceOut
      lsliceEnd = lsliceOut
    else
      tanAlpha = sinAlpha(1) / cosAlpha(1)
      lsliceMin = centerSlice + (lsliceOut - centerSlice) * cosAlpha(1) +  &
          yOffset * sinAlpha(1) - 0.5 * ithickOut * absSinAlf - 1.
      lsliceMax = centerSlice + (lsliceOut - centerSlice) * cosAlpha(1) +  &
          yOffset * sinAlpha(1) + 0.5 * ithickOut * absSinAlf + 2.
      if (debug) print *,'need slices', lsliceMin, lsliceMax
      lsliceMin = max(1, lsliceMin)
      lsliceMax = min(nyProj, lsliceMax)
      lsliceStart = lsliceMin
      if (lsliceMin >= lvertSliceStart .and. lsliceMin <= lvertSliceEnd)  &
          lsliceStart = lvertSliceEnd + 1
      lsliceEnd = lsliceMax
    endif
    !
    ! loop on needed vertical slices
    !
    if (debug) print *,'looping to get', lsliceStart, lsliceEnd

    do lslice = lsliceStart, lsliceEnd
      shiftedGpuLoad = .false.
      !
      ! Load stack with as many lines from projections as will
      ! fit into the remaining space.
      !
      ! Enter loading procedures unless the load is already set for the
      ! current slice
      !
      if (inLoadStart == 0 .or. lslice > lastReady) then
        needStart = 0
        needEnd = 0
        itryEnd = 0
        lastReady = 0
        itry = lslice
        ifEnough = 0
        !
        ! loop on successive output slices to find what input slices are
        ! needed for them; until all slices would be loaded or there
        ! would be no more room
        do while (itry > 0 .and. itry <= nyProj .and. ifEnough == 0 .and. &
            itry <= loadLimit .and. itryEnd - needStart + 1 <= numPlanes)
          lastStart = neededStarts(itry - indNeededBase)
          lastEnd = neededEnds(itry - indNeededBase)
          !
          ! values here must work in single-slice case too
          ! if this is the first time, set needstart
          ! if this load still fits, set needend and lastready
          !
          itryEnd = lastEnd
          if (needStart == 0) needStart = lastStart
          if (itryEnd - needStart + 1 <= numPlanes) then
            needEnd = itryEnd
            lastReady = itry
          else
            ifEnough = 1
          endif
          itry = itry + idelSlice
        enddo
        if (debug) print *,'itryend,needstart,needend,lastready', itryEnd, &
            needStart, needEnd, lastReady
        if (needEnd == 0) call exitError &
            ('INSUFFICIENT STACK SPACE TO DO A SINGLE SLICE')
        !
        ! if some are already loaded, need to shift them down
        !
        numAlready = max(0, inLoadEnd - needStart + 1)
        if (inLoadEnd == 0) numAlready = 0
        ioffset = (needStart - inLoadStart) * inPlaneSize
        do i = indLoadBase, indLoadBase + numAlready * inPlaneSize-1
          array(i) = array(i + ioffset)
        enddo
        if (numAlready .ne. 0 .and. debug) print *,'shifting', needStart, inLoadEnd
        !
        ! If it is also time to load something on GPU, shift existing
        ! data if appropriate and enable copy of filtered data
        if (numGpuPlanes > 0) then
          if (loadGpuStart <= 0 .or. neededEnds(lsliceOut - indNeededBase) > &
              loadGpuEnd) call shiftGpuSetupCopy()
        endif
        !
        ! load the planes in one plane high if stretching
        !
        ibase = indLoadBase + numAlready * inPlaneSize + ipExtraSize
        lstart = needStart + numAlready
        !
        if (debug) print *,'loading', lstart, needEnd
        if (.not. recReproj) then
          do nv = 1, numViews
            istart = ibase
            do nl = lstart, needEnd, idelSlice
              ! Position to read from projection NV at record NL
              iyload = max(0, min(nyProj - 1, nl - 1))
              call iiuSetPosition(1, mapUsedView(nv) - 1, iyload)
              call irdlin(1, array(istart),*999)
              ! Take log if requested
              ! 3/31/04: limit values to .001 time dynamic range
              if (ifLog .ne. 0) then
                do ix = istart, istart + nxProj - 1
                  array(ix) = alog10(max(valMin, array(ix) + baseForLog))
                enddo
              else
                do ix = istart, istart + nxProj - 1
                  array(ix) = array(ix) * exposeWeight(nv)
                enddo
              endif
              !
              ! pad with taper between start and end of line
              call taperEndToStart(istart)
              !
              istart = istart + inPlaneSize
            enddo
            ibase = ibase + nxprj2
          enddo
          !
        else
          !
          ! Load reconstruction for projections
          do nl = lstart, needEnd, idelSlice
            call iiuSetPosition(3, nl - 1, 0)
            call irdpas(3, array(ibase), maxXload + 1 - minXload, &
                ithickReproj, minXload - 1, maxXload - 1, minYreproj - 1, &
                maxYreproj - 1, *999)
            !
            ! Undo the scaling that was used to write, and apply a
            ! scaling that will make the data close for projecting
            array(ibase:ibase + isliceSizeBP - 1) = (array(ibase:ibase + isliceSizeBP -  &
                1) / outScale - outAdd) / filterScale
            ibase = ibase + inPlaneSize
          enddo
        endif
        if (.not.recReproj .and. numSIRTiter <= 0) then
          ibase = indLoadBase + numAlready * inPlaneSize
          do  nl = lstart, needEnd, idelSlice
            call transform(ibase, nl, 1)
            ibase = ibase + inPlaneSize
          enddo
        endif
        inLoadStart = needStart
        inLoadEnd = needEnd
        indLoadEnd = indLoadBase + inPlaneSize *  &
            ((inLoadEnd - inLoadStart) / idelSlice + 1) - 1
      end if
      !
      ! Stack is full.  Now check if GPU needs to be loaded
      if (numGpuPlanes > 0) then
        if (loadGpuStart <= 0 .or. neededEnds(lsliceOut - indNeededBase) > &
            loadGpuEnd) then
          !
          ! Load as much as possible.  If the shift fails the first time
          ! it will set loadGpuStart to 0, call again to recompute for
          ! full load
          if (.not. shiftedGpuLoad) call shiftGpuSetupCopy()
          if (.not. shiftedGpuLoad) call shiftGpuSetupCopy()
          ibase = indLoadBase + (needGpuStart + keepOnGpu - inLoadStart) * inPlaneSize
          if (debug) write(*,'(a,i4,a,i4,a,i4,i10)') 'Loading GPU, #', &
              numLoadGpu, '  lstart', needGpuStart + keepOnGpu, &
              ' start pos base', keepOnGpu + 1, ibase - indLoadBase
          timeStart = wallTime()
          if (gpuLoadProj(array(ibase), numLoadGpu, needGpuStart + keepOnGpu, &
              keepOnGpu + 1) == &
              0) then
            loadGpuStart = needGpuStart
            loadGpuEnd = needGpuEnd
          else
            loadGpuStart = 0
          endif
          if (debug) write(*,'(a,f8.4)') 'Loading time', wallTime() - timeStart
        endif
      endif
      !
      istart = indLoadBase + inPlaneSize * (lslice - inLoadStart) / idelSlice
      if (.not.recReproj) then
        !
        ! Backprojection: Process all views for current slice
        ! If new-style X tilt, set  the Y center based on the slice
        ! number, and adjust the y offset slightly
        !
        if (ifAlpha < 0) ycenOut = ycenFix + (1. / cosAlpha(1) - 1.) * yOffset &
            - nint(tanAlpha * (lslice - centerSlice))
        if (numSIRTiter == 0) then
          !
          ! print *,indLoadBase, lslice, inloadstr
          ! print *,'projecting', lslice, ' at', istart, ', ycen =', ycenOut
          call project(istart, lslice)
          if (maskEdges) call maskSlice(indOutSlice, ithickBP)
          !
          ! move vertical slice into ring buffer, adjust ring variables
          !
          if (ifAlpha < 0) then
            ! print *,'moving slice to ring position', nextfreevs
            ioffset = isliceSizeBP + (nextFreeVertSlice - 1) * ithickBP * iwidth
            do i = indOutSlice, indOutSlice + ithickBP * iwidth - 1
              array(i + ioffset) = array(i)
            enddo
            call manageRing(numVertNeeded, numVertSliceInRing, nextFreeVertSlice, &
                lvertSliceStart, lvertSliceEnd, lslice)
          endif
        else
          !
          ! for SIRT, either do the zero iteration here, load single slice
          ! into read-in slice spot, or decompose vertical slice from read-
          ! in slices into spot in  vertical slice ring buffer.
          ! Descale them by output scaling only
          ! Set up for starting from read-in data
          recScale = filterScale * numViews
          iset = 1
          reprojFill = dmeanIn
          if (sirtFromZero) then
            !
            ! Doing zero iteration: set up location to place slice
            ibaseSIRT = ireadBase
            if (ifAlpha < 0) then
              ibaseSIRT = indOutSlice + isliceSizeBP + (nextFreeVertSlice - 1) * &
                  ithickBP * iwidth
              call manageRing(numVertNeeded, numVertSliceInRing, nextFreeVertSlice, &
                  lvertSliceStart, lvertSliceEnd, lslice)
            endif
            !
            ! Copy and filter the lines needed by filter set 1
            do iv = 1, numViews
              j = indWorkPlane + (iv - 1) * nxprj2
              k = istart + (iv - 1) * nxprj2
              do i = 0, nxprj2 - 3
                array(j + i) = array(k + i)
              enddo
            enddo
            call transform(indWorkPlane, lslice, 1)
            edgeFill = edgeFillOrig
            !
            ! Backproject and maskEdges/taper edges
            call project(indWorkPlane, lslice)
            if (maskEdges) call maskSlice(indOutSlice, ithickBP)
            array(ibaseSIRT:ibaseSIRT + isliceSizeBP - 1) = array(indOutSlice:mapEnd) / &
                recScale
            iset = 2
            if (ifOutSirtRec == 3) then
              print *,'writing zero iteration slice', lslice
              call writeInternalSlice(5, ibaseSIRT)
            endif

          else if (ifAlpha == 0) then
            !
            ! Read in single slice
            call iiuSetPosition(3, lslice - 1, 0)
            call irdsec(3, array(ireadBase), *999)
            dsum = 0.
            dpix = 0.
            do i = 0, ithickOut * iwidth - 1
              array(i + ireadBase) = (array(i + ireadBase) / outScale - outAdd)
              dsum = dsum + array(i + ireadBase)
              dpix = dpix + 1.
            enddo
            if (debug) print *,'Loaded mean= ', dsum / dpix, '   pmean= ', dmeanIn
            ibaseSIRT = ireadBase
          else
            !
            ! Vertical slices: read file directly if it is vertical slices
            ibaseSIRT = indOutSlice + isliceSizeBP + (nextFreeVertSlice - 1) * ithickBP &
                * iwidth
            if (vertSirtInput) then
              if (debug) print *,'loading vertical slice', lslice
              call iiuSetPosition(3, lslice - 1, 0)
              call irdsec(3, array(ibaseSIRT), *999)
            else
              !
              ! Decompose vertical slice from read-in slices
              ! First see if necessary slices are loaded in ring
              vertSliceCen = lslice - centerSlice
              vertYcenFix = (ithickBP / 2 + 0.5) - nint(tanAlpha * (lslice - centerSlice)) &
                  + yOffset / cosAlpha(1)
              riBot = centerSlice + vertSliceCen * cosAlpha(1) + (1 - vertYcenFix) * &
                  sinAlpha(1)
              riTop = riBot + (ithickBP - 1) * sinAlpha(1)
              lriMin = max(1, floor(min(riBot, riTop)))
              lriMax = min(nyProj, ceiling(max(riBot, riTop)))

              loadStart = lriMin
              if (lriMin >= lReadStart .and. lriMin <= lReadEnd) loadStart = lReadEnd + 1
              loadEnd = lriMax
              !
              ! Read into ring and manage the ring pointers
              if (debug .and. loadEnd >= loadStart) &
                  print *,'reading into ring', loadStart, loadEnd
              do lri = loadStart, loadEnd
                call iiuSetPosition(3, lri - 1, 0)
                ioffset = ireadBase + (nextReadFree - 1) * ithickOut * iwidth
                call irdsec(3, array(ioffset), *999)
                do i = 0, ithickOut * iwidth - 1
                  array(i + ioffset) = (array(i + ioffset) / outScale - outAdd)
                enddo
                call manageRing(numReadNeed, numReadInRing, nextReadFree, lReadStart, &
                    lReadEnd, lri)
              enddo
              iringStart = 1
              if (numReadInRing == numReadNeed) iringStart = nextReadFree
              call decompose(lslice, lReadStart, lReadEnd, iringStart, ibaseSIRT)
            endif
            if (ifOutSirtRec == 4) then
              print *,'writing decomposed slice', lslice
              call writeInternalSlice(5, ibaseSIRT)
            endif
            call manageRing(numVertNeeded, numVertSliceInRing, nextFreeVertSlice, &
                lvertSliceStart, lvertSliceEnd, lslice)
          endif
          !
          ! Ready to iterate.  First get the starting slice mean
          call iclden(array(ibaseSIRT), iwidth, ithickBP, 1, iwidth, 1, &
              ithickBP, unscaledMin, unscaledMax, firstMean)
          ! It seems to give less edge artifacts using the mean all the time
          ! if (sirtFromZero) rpfill = firstmean
          reprojFill = firstMean
          !
          do isirtIter = 1, numSIRTiter
            ierr = 1
            timeStart = wallTime()

            if (useGPU) then
              ierr = gpuReprojOneSlice(array(ibaseSIRT), array(indWorkPlane), &
                  sinReproj, cosReproj, ycenOut, numViews, reprojFill)
            endif
            if (ierr .ne. 0) then
              do iv = 1, numViews
                i = (iv - 1) * iwidth + 1
                j = indWorkPlane + (iv - 1) * nxprj2
                ! call reproject(array(ibaseSIRT), iwidth, ithickBP, iwidth, &
                ! sinReproj(iv), cosReproj(iv), xRayStart(i), yRayStart(i), &
                ! numPixInRay(i), maxRayPixels(iv), dmeanIn, array(j), 1, 1)
                call reprojOneAngle(array(ibaseSIRT), array(j), 1, 1, 1, &
                    cosReproj(iv), -sinReproj(iv), 1., 0., cosReproj(iv), &
                    iwidth, ithickBP, iwidth * ithickBP, iwidth, 1, 1, 0., 0., &
                    iwidth / 2 + 0.5, ycenOut, iwidth / 2 + 0.5, &
                    centerSlice, 0, 0., 0., reprojFill)
              enddo
            endif
            if (debug) write(*,'(a,f8.4)') 'Reproj time = ', wallTime() - timeStart
            if (isirtIter == numSIRTiter .and. ifOutSirtProj == 1) &
                call writeInternalSlice(4, indWorkPlane)
            !
            ! Subtract input projection lines from result.  No scaling
            ! needed here; these intensities should be in register
            ! Taper ends
            ! dsum = 0.
            do iv = 1, numViews
              j = indWorkPlane + (iv - 1) * nxprj2
              k = istart + (iv - 1) * nxprj2
              do i = 0, iwidth - 1
                array(j + i) = array(j + i) - array(k + i)
              enddo
              call taperEndToStart(j)
            enddo
            if (isirtIter == numSIRTiter .and. ifOutSirtProj == 2) &
                call writeInternalSlice(4, indWorkPlane)
            !
            ! filter the working difference lines and get a rough value for
            ! the edgeFill.  Trying to get this better did not prevent edge
            ! artifacts - the masking was needed
            call transform(indWorkPlane, lslice, iset)
            dsum = 0.
            do iv = 1, numViews
              j = indWorkPlane + (iv - 1) * nxprj2
              dsum = dsum + array(j + iwidth + numPad / 2 - 1)
            enddo
            edgeFill = dsum / numViews
            ! print *,'   (for diff bp) edgefill =', edgeFill
            !
            ! Backproject the difference
            call project(indWorkPlane, lslice)
            if (isirtIter == numSIRTiter .and. ifOutSirtRec == 1) then
              array(indOutSlice:mapEnd) = (array(indOutSlice:mapEnd)  + outAdd * &
                  recScale) * (outScale / recScale)
              print *,'writing bp difference', lslice
              call writeInternalSlice(5, indOutSlice)
              array(indOutSlice:mapEnd) = array(indOutSlice:mapEnd) /  &
                  (outScale / recScale) - outAdd * recScale
            endif
            !
            ! Accumulate report values
            if (iterForReport > 0) call sampleForReport(array(indOutSlice), &
                lslice, ithickBP, isirtIter, outScale / recScale, outAdd * recScale)
            !
            ! Subtract from input slice
            ! outScale to account for difference between ordinary scaling
            ! for output by 1 / numViews and the fact that input slice was
            ! descaled by 1/filterScale
            i = ibaseSIRT + isliceSizeBP - 1
            if (isignConstraint == 0) then
              array(ibaseSIRT:i) = array(ibaseSIRT:i) - array(indOutSlice:mapEnd) / &
                  recScale
            else if (isignConstraint < 0) then
              array(ibaseSIRT:i) = min(0., array(ibaseSIRT:i) - &
                  array(indOutSlice:mapEnd) / recScale)
            else
              array(ibaseSIRT:i) = max(0., array(ibaseSIRT:i) - &
                  array(indOutSlice:mapEnd) / recScale)
            endif
            if (maskEdges) call maskSlice(ibaseSIRT, ithickBP)
            if (isirtIter == numSIRTiter - 2 .and. ifOutSirtRec == 2) then
              print *,'writing internal slice', lslice
              call writeInternalSlice(5, ibaseSIRT)
            endif
            if (isirtIter == numSIRTiter .and. saveVertSlices .and. &
                lslice >= isliceStart .and. lslice <= isliceEnd) then
              if (debug) print *,'writing vertical slice', lslice
              ierr = parWrtSetCurrent(2)
              call iclden(array(ibaseSIRT), iwidth, ithickBP, 1, iwidth, 1, ithickBP, &
                  tmin, tmax, tmean)
              dmin6 = min(dmin6, tmin)
              dmax6 = max(dmax6, tmax)
              call parWrtSec(6, array(ibaseSIRT))
              ierr = parWrtSetCurrent(1)
            endif
            !
            ! Adjust fill value by change in mean
            ! Accumulate mean of starting vertical slices until one output
            ! slice has been composed, and set composeFill from mean
            if (isirtIter < numSIRTiter .or. ifAlpha < 0) then
              call iclden(array(ibaseSIRT), iwidth, ithickBP, 1, iwidth, 1, &
                  ithickBP, unscaledMin, unscaledMax, curMean)
              !
              ! Doing it the same way with restarts as with going from 0 seems
              ! to prevent low frequency artifacts
              ! if (sirtFromZero) then
              reprojFill = curMean
              ! else
              ! rpfill = dmeanIn + curmean - firstmean
              ! endif
              if (isirtIter == numSIRTiter .and. .not. composedOne) then
                vertSum = vertSum + curMean
                numVertSum = numVertSum + 1
              endif
            endif
            if (numVertSum > 0) composeFill = vertSum / numVertSum
          enddo
        endif
      endif
    enddo
    !
    if (.not.recReproj) then
      if (ifAlpha < 0) then
        !
        ! interpolate output slice from vertical slices
        !
        iringStart = 1
        if (numVertSliceInRing == numVertNeeded) iringStart = nextFreeVertSlice
        if (numVertSliceInRing > 0) then
          !
          ! If there is anything in ring to use, then we can compose an
          ! output slice, first dumping any deferred fill slices
          call dumpFillSlices()
          ! print *,'composing', lsliceout, ' from', lvsstart, lvsend, iringstart
          call compose(lsliceOut, lvertSliceStart, lvertSliceEnd, 1, iringStart,  &
              composeFill)
          composedOne = .true.
        else
          !
          ! But if there is nothing there, keep track of a range of slices
          ! that need to be filled - after the composeFill value is set
          if (lfillStart < 0) lfillStart = lsliceOut
          lfillEnd = lsliceOut
        endif
      endif
      !
      ! Dump slice unless there are fill slices being held
      if (lfillStart < 0) call dumpUpdateHeader(lsliceOut)
      lsliceOut = lsliceOut + idelSlice
    else
      !
      ! REPROJECT all ready slices to minimize file mangling
      lsliceProjEnd = lastReady
      call reprojectRec(lsliceOut, lsliceProjEnd, inLoadStart, inLoadEnd, dmin, &
          dmax, dtot8)
      lsliceOut = lsliceProjEnd + 1
    endif
  enddo
  !
  ! End of main loop
  !-----------------
  !
  call dumpFillSlices()
  ! Close files
  call iiuClose(1)
  pixelTot = float(numSlices) * iwidth * ithickOut
  if (reprojBP .or. recReproj) pixelTot = float(numSlices) * iwidth * numReproj
  dmean = dtot8 / pixelTot
  !
  ! get numbers used to compute scaling report, and then fix min/max if
  ! necessary
  unscaledMin = dmin / outScale - outAdd
  unscaledMax = dmax / outScale - outAdd
  truncations = .false.
  extremes = .false.
  if (newMode == 0) then
    truncations = dmin < 0. .or. dmax > 256.
    dmin = max(0., dmin)
    dmax = min(255., dmax)
  else if (newMode == 1) then
    extremes = useGPU .and. max(abs(dmin), abs(dmax)) > effectiveScale * 3.e5
    truncations = dmin < -32768. .or. dmax > 32768.
    dmin = max(-32768., dmin)
    dmax = min(32767., dmax)
  endif
  if (minTotSlice <= 0) then
    if (perpendicular .and. intervalHeadSave > 0 .and. .not.(reprojBP .or. recReproj)) &
        then
      nxyzTmp(3) = numSlices
      call iiuAltSize(2, nxyzTmp, nxyzst)
    endif
    call iiuWriteHeader(2, title, 1, dmin, dmax, dmean)
    if (.not.(reprojBP .or. recReproj)) then
      write(6, 930) 'reconstruction'
    else
      write(6, 930) 'reprojection'
    endif
    call irdhdr(2, nxyzTmp, nxyzst, mode, dmin, dmax, dmean)
    if (saveVertSlices) call iiuWriteHeader(6, title, 1, dmin6, dmax6,  &
        (dmin6 + dmax6) / 2.)
  else
    write(*,'(a,3g15.7,f15.0)') 'Min, max, mean, # pixels=', dmin, dmax, dmean, pixelTot
  endif
  if (parallelHDF) then
    if (iiuParWrtFlushBuffers(2) .ne. 0) &
        call exitError("Finishing writing to output HDF file")
    if (saveVertSlices) then
      ierr = parWrtSetCurrent(2)
      if (iiuParWrtFlushBuffers(6) .ne. 0) &
          call exitError("Finishing writing to vertical slice HDF file")
      ierr = parWrtSetCurrent(1)
    endif
    call parWrtClose()
  endif
  call iiuClose(2)
  if (saveVertSlices) call iiuClose(6)
  if (.not.(reprojBP .or. recReproj .or. numSIRTiter > 0)) then
    recScale = numViews * 235. / (unscaledMax - unscaledMin)
    recAdd = (10. *(unscaledMax - unscaledMin) / 235. -unscaledMin) / numViews
    write(6, 905) 'bytes (10-245)', recAdd, recScale
    recScale = numViews * 30000. / (unscaledMax - unscaledMin)
    recAdd = (-15000. *(unscaledMax - unscaledMin) / 30000. -unscaledMin) / numViews
    write(6, 905) '-15000 to 15000', recAdd, recScale
    write(6, 910) numSlices
  endif
  if (extremes) write(6, 911)
  if (truncations .and. .not. extremes) write(6, 912)
  if (useGPU) call gpuDone()
  if (iterForReport > 0) then
    do i = 1, max(1, numSIRTiter)
      if (reportVals(3, i) > 0) write(6, 908) i + iterForReport - 1, &
          isliceStart, isliceEnd, reportVals(1, i) / reportVals(3, i), &
          reportVals(2, i) / reportVals(3, i)
    enddo
  endif
  if (ifOutSirtProj > 0) then
    call iiuWriteHeader(4, title, 1, dmin4, dmax4, (dmin4 + dmax4) / 2.)
    call iiuClose(4)
  endif
  if (ifOutSirtRec > 0) then
    call iiuWriteHeader(5, title, 1, dmin5, dmax5, (dmin5 + dmax5) / 2.)
    call iiuClose(5)
  endif
  call exit(0)
999 write(6, 920) mapUsedView(nv), nl
  call exit(1)
  !
  !
900 format(//,' STACK LOADING' &
      /,' -------------' &
      //,' Total stack size           ',i11,/ &
      /,1x,a20,'         ',i9/,a,//, &
      ' Output slice                 ',i9,/)
901 format(' Stretching buffer            ',i9,/)
902 format(1x,i4,  ' Untilted slices         ',i9,/)
903 format(1x,i4,' Transposed projections  ',i9,/)
904 format(1x,i4,' Slice(s) for SIRT       ',i9,/)
907 format(' Reprojection lines for SIRT  ',i9,/)
905 format(/,' To scale output to ',a,', use SCALE to add',f12.3, &
      ' and scale by',f12.5)
908 format(/,'Iter ', i4.3, ', slices', 2i6, ', diff rec mean&sd:', 2f15.3)
910 format(//' Reconstruction of',i5,' slices complete.', &
      //,1x,78('-'))
911 format(/,'WARNING: TILT - EXTREMELY LARGE VALUES OCCURRED AND VALUES ', &
      'WERE TRUNCATED WHEN OUTPUT TO FILE; THERE COULD BE ERRORS IN GPU ', &
      'COMPUTATION.  RUN gputilttest',/)
912 format(/,'WARNING: TILT - SOME VALUES WERE TRUNCATED WHEN OUTPUT TO THE', &
      ' FILE; CHECK THE OUTPUT SCALING FACTOR',/)
920 format(//' ERROR: TILT -  reading in view',i3,' for slice' &
      ,i5,/)
930 format(//' Header on ',a,' file' / &
      ' --------------------------------'//)

CONTAINS
  !
  ! Determine parameters for loading GPU and make call to shift existing
  ! data if any is to be retained
  subroutine shiftGpuSetupCopy()
    integer*4 gpuShiftProj
    needGpuStart = neededStarts(lsliceOut - indNeededBase)
    needGpuEnd = min(needGpuStart + numGpuPlanes - 1,  needEnd)
    keepOnGpu = 0
    if (loadGpuStart > 0) keepOnGpu = loadGpuEnd + 1 - needGpuStart
    numLoadGpu = needGpuEnd + 1 - needGpuStart - keepOnGpu
    if (debug) write(*,'(a,i4,a,i4,a,i4)') 'Shifting GPU, #', numLoadGpu, &
        '  lstart', needGpuStart + keepOnGpu, ' start pos', keepOnGpu + 1
    timeStart = wallTime()
    if (gpuShiftProj(numLoadGpu, needGpuStart + keepOnGpu, keepOnGpu + 1) &
        == 0) then
      shiftedGpuLoad = .true.
    else
      loadGpuStart = 0
    endif
    if (debug) write(*,'(a,f8.4)') 'Shifting time', wallTime() - timeStart
    return
  end subroutine shiftGpuSetupCopy

  ! For writing two kinds of test output
  !
  subroutine writeInternalSlice(isUnit, iwriteStart)
    integer*4 isUnit, iwriteStart
    real*4 istmin, istmax, istmean
    if (isUnit == 4) then
      do iv = 1, numViews
        j = iwriteStart + (iv - 1) * nxprj2
        call iiuSetPosition(4, iv - 1, lslice - isliceStart)
        call iwrlin(4, array(j))
      enddo
      call iclden(array(iwriteStart), nxprj2, iv, 1, iwidth, 1, numViews, &
          istmin, istmax, istmean)
      dmin4 = min(dmin4, istmin)
      dmax4 = max(dmax4, istmax)
    else
      call iiuWriteSection(5, array(iwriteStart))
      call iclden(array(iwriteStart), iwidth, ithickBP, 1, iwidth, 1, ithickBP, &
          istmin, istmax, istmean)
      dmin5 = min(dmin5, istmin)
      dmax5 = max(dmax5, istmax)
      nz5 = nz5 + 1
      call ialsiz_sam_cel(5, iwidth, ithickBP, nz5)
    endif
    return
  end subroutine writeInternalSlice

  ! Dump a slice and update the header periodically if appropriate
  !
  subroutine dumpUpdateHeader(lsliceDump)
    integer*4 lsliceDump
    !
    ! Write out current slice
    call dumpSlice(lsliceDump, dmin, dmax, dtot8)
    ! DNM 10/22/03:  Can't use flush in Windows/Intel because of sample.com
    ! call flush(6)
    !
    ! write out header periodically, restore writing position
    if (perpendicular .and. intervalHeadSave > 0 .and. .not.reprojBP .and. &
        minTotSlice <= 0) then
      numSliceOut = numSliceOut + 1
      nxyzTmp(1) = iwidth
      nxyzTmp(2) = ithickOut
      nxyzTmp(3) = numSliceOut
      if (mod(numSliceOut, intervalHeadSave) == 1) then
        call iiuAltSize(2, nxyzTmp, nxyzst)
        dmean = dtot8 / (float(numSliceOut) * iwidth * ithickOut)
        call iiuWriteHeader(2, title, -1, dmin, dmax, dmean)
        call parWrtPosn(2, numSliceOut, 0)
      endif
    endif
  end subroutine dumpUpdateHeader

  ! If there are fill slices that haven't been output yet, dump them now
  ! and turn off signal that they exist
  subroutine dumpFillSlices()
    if (lfillStart >= 0) then
      do i = lfillStart, lfillEnd, idelSlice
        array(indOutSlice:indOutSlice + iwidth * ithickOut - 1) = composeFill
        call dumpUpdateHeader(i)
      enddo
      lfillStart = -1
    endif
  end subroutine dumpFillSlices

  ! END OF MAIN PROGRAM UNIT
end program tilt


! Taper intensities across pad region from end of line to start
! This is not a perfect solution, since it turns a single-pixel deviation at the edge
! of the image into a step, which gets accentuated by the radial filter and leaves a
! little artifact at the edge.  It used to average over 2 pixels, but this is actually
! slightly worse than just using the single pixels at each end.  Mirroring a few pixels
! on each edge before the taper was also worse. In any case, it is much better than
! padding with a constant.
!
subroutine taperEndToStart(indStart)
  use tiltvars
  implicit none
  real*4 xsum, startMean, endMean, f
  integer*4 ipad, nsum, ix, indStart, numAverage
  numAverage = 1
  if (numAverage > 1) then
    nsum = 0
    xsum = 0.
    do ix = indStart, indStart + min(numAverage, nxProj - 1)
      nsum = nsum + 1
      xsum = xsum + array(ix)
    enddo
    startMean = xsum / nsum
    nsum = 0
    xsum = 0.
    do ix = indStart + max(0, nxProj - numAverage - 1), indStart + nxProj - 1
      nsum = nsum + 1
      xsum = xsum + array(ix)
    enddo
    endMean = xsum / nsum
  else
    startMean = array(indStart)
    endMean = array(indStart + nxProj - 1)
  endif
  do ipad = 1, numPad
    f = ipad / (numPad + 1.)
    array(indStart + nxProj + ipad-1) = f * startMean + (1. -f) * endMean
  enddo
  return
end subroutine taperEndToStart

subroutine manageRing(numVertNeeded, numVertSliceInRing, nextFreeVertSlice, &
    lvertSliceStart, lvertSliceEnd, lslice)
  implicit none
  integer*4 numVertNeeded, numVertSliceInRing, nextFreeVertSlice, lvertSliceStart
  integer*4 lvertSliceEnd, lslice
  if (numVertSliceInRing < numVertNeeded) then
    if (numVertSliceInRing == 0) lvertSliceStart = lslice
    numVertSliceInRing = numVertSliceInRing + 1
  else
    lvertSliceStart = lvertSliceStart + 1
  endif
  lvertSliceEnd = lslice
  nextFreeVertSlice = nextFreeVertSlice + 1
  if (nextFreeVertSlice > numVertNeeded) nextFreeVertSlice = 1
  return
end subroutine manageRing



! --------------------------------------------------------------------
subroutine radialWeights(iradMaxIn, radFallIn, ifilterSet, ifMultByGaussian)
  ! -----------------------------
  !
  ! Set Radial Transform weighting
  ! Linear ramp plus Gaussian fall off
  use tiltvars
  implicit none
  integer*4 iradMaxIn, ifilterSet, ifMultByGaussian
  real*4 radFallIn
  integer*4 nxProjPad, iradEnd, iradMax, iv, iw, indBase, i, iradLimit, jv, imax
  integer*4 irampEnd
  real*4 stretch, avgInterval, atten, sumInterval, wsum, z, arg, sirtFrac, zmax, freq
  real*4 diffMin, diff, attenSum, tabFac, sinDiff, fact, exactReachTopMean, fakeAlpha
  real*4 fakeIterUse, fakeMatchAdd, radFall
  real*4, allocatable :: wgtAtten(:)
  real*4 degToRad /0.0174533/

  !
  allocate(wgtAtten(limView), stat = i)
  call memoryErrorUC(i, 'ARRAY FOR WEIGHTS PER VIEW')
  sirtFrac = 0.99
  !
  ! Empirically determined parameters for fake SIRT.  The alpha value is good for
  ! low iteration numbers but then iteration needs to be scaled to match SIRT
  fakeMatchAdd = 0.3
  fakeAlpha = 0.00195
  fakeIterUse = numFakeSIRTiter
  if (numFakeSIRTiter > 15) fakeIterUse = 15 + 0.8 * (numFakeSIRTiter - 15)
  if (numFakeSIRTiter > 30) fakeIterUse = 27 + 0.6 * (numFakeSIRTiter - 30)
  nxProjPad = nxProj + 2 + numPad
  iradEnd = nxProjPad / 2
  stretch = float(nxProj + numPad) / nxProj
  iradMax = nint(iradMaxIn * stretch)
  radFall = radFallIn * stretch
  !
  ! Compute the ramp to the start of the falloff by default, or all the way out if 
  ! multiplying by Gaussian
  irampEnd = min(iradMax, iradEnd)
  if (ifMultByGaussian > 0) irampEnd = iradEnd
  avgInterval = 1.
  attenSum = 0.
  zeroWeight = 0.
  exactReachTopMean = 0.
  if (numTiltIncWgt > 0 .and. numWgtAngles > 1) then
    avgInterval = (wgtAngles(numWgtAngles) - wgtAngles(1)) / (numWgtAngles - 1)
    if (debug) write(6, 401)
401 format(/' View  Angle Weighting')
  endif
  !
  ! Set up the attenuations for the weighting angles
  do iv = 1, numWgtAngles
    atten = 1.
    if (numTiltIncWgt > 0. .and. numWgtAngles > 1) then
      sumInterval = 0
      wsum = 0.
      do iw = 1, numTiltIncWgt
        if (iv - iw > 0) then
          wsum = wsum + tiltIncWgts(iw)
          sumInterval = sumInterval + tiltIncWgts(iw) * (wgtAngles(iv + 1 - iw) - &
              wgtAngles(iv - iw))
        endif
        if (iv + iw <= numViews) then
          wsum = wsum + tiltIncWgts(iw)
          sumInterval = sumInterval + tiltIncWgts(iw) * (wgtAngles(iv + iw) - &
              wgtAngles(iv + iw - 1))
        endif
      enddo
      atten = atten * (sumInterval / wsum) / avgInterval
    endif
    wgtAtten(iv) = atten
  enddo
  !
  ! Set up linear ramp
  if (numExactCycles == 0 .or. flatFrac > 0) then
    do iv = 1, numViews
      !
      ! Get weighting from nearest weighting angle
      atten = 1.
      if (numTiltIncWgt > 0 .and. numWgtAngles > 1) then
        diffMin = 1.e10
        do iw = 1, numWgtAngles
          diff = abs(angles(iv) - wgtAngles(iw))
          if (diff < diffMin) then
            diffMin = diff
            atten = wgtAtten(iw)
          endif
        enddo
        if (debug) write(6, 402) iv, angles(iv) / degToRad, atten
402     format(i4,f8.2,f10.5)
      endif
      !
      ! Take negative if subtracting
      if (numViewSubtract > 0) then
        do i = 1, numViewSubtract
          if (iviewSubtract(i) == 0 .or. mapUsedView(iv) == iviewSubtract(i)) &
              atten = -atten
        enddo
      endif
      !
      attenSum = attenSum + atten
      indBase = (iv - 1 + (ifilterSet - 1) * numViews) * nxProjPad
      do  i = 1, irampEnd
        ! This was the basic filter
        ! ARRAY(ibase+2*I-1) =atten*(I-1)
        ! This is the mixture of the basic filter and a flat filter with
        ! a scaling that would give approximately the same output magnitude
        freq = (i - 1.) / nxProjPad
        if (numFakeSIRTiter <= 0 .or. i == 1 .or. freq < fakeAlpha) then
          array(indBase + 2 * i - 1) = atten * ((1. - flatFrac) * (i - 1.) + &
              flatFrac * filterScale)
        else
          array(indBase + 2 * i - 1) = atten * (i - 1.) *  &
              (1. - (1. - fakeAlpha / freq)**(fakeIterUse + fakeMatchadd))
        endif
        !
        ! This 0.2 is what Kak and Slaney's weighting function gives at 0
        if (i == 1) array(indBase + 2 * i - 1) = atten * 0.2
        !
        ! This is the SIRT filter, which divides the error equally among
        ! the pixels on a ray.
        if (flatFrac > 1) array(indBase + 2 * i - 1) = sirtFrac * atten * &
            filterScale / (ithickBP / cosBeta(iv))
        !
        ! And just the value and compute the mean zero weighting
        array(indBase + 2 * i) = array(indBase + 2 * i - 1)
        if (i == 1) zeroWeight = zeroWeight + array(indBase + 2 * i - 1) / numViews
      enddo
    enddo
    if (debug) print *,'Mean weighting factor', attenSum / numViews
  else
    !
    ! Exact filters
    ! For each actual view, loop on weighting angles and compute influence factors
    array(1: numViews * nxProjPad) = 0.
    do iv = 1, numViews
      indBase = (iv - 1) * nxProjPad
      do jv = 1, numWgtAngles
        sinDiff = sin(abs(wgtAngles(jv) - angles(iv)))
        fact = sinDiff * exactObjSize / (nxProj + numPad)
        iradLimit = numExactCycles / max(1.e-6, fact)
        tabFac = exactSamples * fact
        do i = 1, min(iradLimit, irampEnd)
          array(indBase + 2 * i) = array(indBase + 2 * i) +  &
              exactTable(nint((i - 1) * tabFac))
        enddo
      enddo
    enddo
    !
    ! Invert the sum of influence factors
    do iv = 1, numViews
      indBase = (iv - 1) * nxProjPad
      zmax = 0.
      do i = 1, irampEnd
        if (array(indBase + 2 * i) == 0.) array(indBase + 2 * i) = 1.
        z = iradEnd / array(indBase + 2 * i)
        array(indBase + 2 * i - 1) = z
        array(indBase + 2 * i) = z
        if (z > zmax + 0.001) then
          zmax = z
          imax = i
        endif
        !
        ! Limit the zero frequency weighting as above to avoid different mean and max
        if (i == 1) then
          array(indBase + 2 * i - 1) = min(0.2, z)
          array(indBase + 2 * i) = min(0.2, z)
          zeroWeight = zeroWeight + array(indBase + 2 * i - 1) / numViews
        endif
      enddo
      exactReachTopMean = exactReachTopMean + 0.5 * (imax - 1.) / (iradEnd * numViews)
    enddo
    write(*,'(/,a,f6.3,a)')'Exact filters reach highest point at mean frequency of',  &
        exactReachTopMean, ' cycles/pixel'
  endif
  !
  ! Set up Gaussian
  do i = iradMax + 1, iradEnd
    atten = 0.
    arg = float(i - iradMax) / max(0.01, radFall)
    if (arg < 8) atten = exp(-arg * arg / 2.)
    !
    ! Scale current value if multiplying, or last value if just falling off
    if (ifMultByGaussian > 0) then
      imax = i
    else
      imax = iradMax
    endif
    indBase = 0
    do iv = 1, numViews
      z = atten * array(indBase + 2 * imax)
      array(indBase + 2 * i - 1) = z
      array(indBase + 2 * i) = z
      indBase = indBase + nxProjPad
    enddo
  enddo
  if (debug) then
    ! do iv = 1, numViews
    iv = numViews / 2
    indBase = (iv - 1) * nxProjPad
    do i = 1, 19
      write(*,'(a,2i5,f12.5)')'RF:',iv,i,array(indBase + 2 * i)
    enddo
    do i = 20, iradEnd, 10
      write(*,'(a,2i5,f12.5)')'RF:',iv,i,array(indBase + 2 * i)
    enddo
    ! enddo
  endif
  return
end subroutine radialWeights
!
!
! ---------------------------------------------------------------------
subroutine maskPrep(lslice)
  ! ----------------
  !
  ! This subroutine prepares the limits of the slice width to be computed
  ! if masking is used
  use tiltvars
  implicit none
  real*4 radiusLeft, radiusRight, y, yy, ycenUse
  integer*4 i, ixLeft, ixRight, lslice
  !
  ! Compute left and right edges of unmasked area
  if (maskEdges) then
    !
    ! Adjust the Y center for alpha tilt (already adjusted for ifAlpha < 0)
    ycenUse = ycenOut
    if (ifAlpha > 0) ycenUse = ycenOut + (1. / cosAlpha(1) - 1.) * yOffset &
        - nint((lslice - centerSlice) * sinAlpha(1) / cosAlpha(1))
    !
    ! Get square of radius of arcs of edge of input data from input center
    radiusLeft = (xcenIn + axisXoffset - 1)**2
    radiusRight = (nxProj - xcenIn - axisXoffset)**2
    do i = 1, ithickBP
      y = i - ycenUse
      yy = min(y * y, radiusLeft)
      !
      ! get distance of X coordinate from center, subtract from or add to
      ! center and round up on left, down on right, plus added maskEdges pixels
      ixLeft = xcenOut + 1. -sqrt(radiusLeft - yy)
      ixUnmaskedSE(1, i) = max(1, ixLeft + numExtraMaskPix)
      yy = min(y * y, radiusRight)
      ixRight = xcenOut + sqrt(radiusRight - yy)
      ixUnmaskedSE(2, i) = min(iwidth, ixRight - numExtraMaskPix)
    enddo
    !-------------------------------------------------
    ! If no maskEdges
  else
    do i = 1, ithickBP
      ixUnmaskedSE(1, i) = 1
      ixUnmaskedSE(2, i) = iwidth
    enddo
  END if
  return
end subroutine maskprep
!
!
! ---------------------------------------------------------------------
subroutine transform(ibaseInd, lslice, ifilterSet)
  ! ----------------------------
  !
  ! This subroutine applies a one-dimensional Fourier transform to
  ! all views corresponding to a given slice, applies the radial
  ! weighting function and then applies an inverse Fourier transform.
  !
  use tiltvars
  implicit none
  integer*4 nxProjPad, indStart, index, indRad, iv, i, ibFrom, ibaseTo, ixp
  integer*4 ixpP1, ixpM1, ixpP2, ibaseInd, lslice, ifilterSet
  real*4 x, xp, dx, dxM1, v4, v5, v6, a, c, denNew, dxDxM1, diffMax
  real*4 fx1, fx2, fx3, fx4
  real*8 wallTime, tstart
  integer*4 gpuFilterLines

  nxProjPad = nxProj + 2 + numPad
  indStart = ibaseInd + ipExtraSize
  tstart = wallTime()
  index = 1
  if (useGPU) then
    index = gpuFilterLines(array(indStart), lslice, ifilterSet)
  endif
  if (index .ne. 0) then
    !
    ! Apply forward Fourier transform
    call odfft(array(indStart), nxProj + numPad, numViews, 0)
    !
    ! Apply Radial weighting
    index = indStart
    indRad = 1 + (ifilterSet - 1) * nxProjPad * numViews
    v4 = 0.
    v5 = 0.
    do  iv = 1, numViews
      do i = 1, nxProjPad
        v4 = v4 + array(index)
        v5 = v5 + array(indRad)
        array(index) = array(index) * array(indRad)
        indRad = indRad + 1
        index = index + 1
      enddo
    enddo
    !
    ! Apply inverse transform
    call odfft(array(indStart), nxProj + numPad, numViews, 1)
  endif
  if (debug) write(*,'(a,f8.4)') 'Filter time', wallTime() - tstart
  if (ipExtraSize == 0) return
  !
  ! do cosine stretch and move down one plane
  ! Use cubic interpolation a la cubinterp
  !
  ! print *,'istart, ibase', istart, ibase
  do iv = 1, numViews
    ibFrom = indStart + (iv - 1) * nxProjPad - 1
    ibaseTo = ibaseInd + indStretchLine(iv) - 1
    ! print *,nv, ibfrom, ibto, nxStretched(nv)
    diffMax = 0.
    if (interpOrdStretch == 1) then
      !
      ! linear interpolation
      !
      do i = 1, nxStretched(iv)
        x = i / float(interpFacStretch) + stretchOffset(iv)
        xp = min(max(1., x * cosBeta(iv)), float(nxProj))
        ixp = xp
        dx = xp - ixp
        ixp = ixp + ibFrom
        ixpP1 = min(ixp + 1, nxProj + ibFrom)
        dxM1 = dx - 1.
        array(ibaseTo + i) = -dxM1 * array(ixp) + dx * array(ixpP1)
      enddo
    else if (interpOrdStretch == 2) then
      !
      ! quadratic
      !
      do i = 1, nxStretched(iv)
        x = i / float(interpFacStretch) + stretchOffset(iv)
        xp = min(max(1., x * cosBeta(iv)), float(nxProj))
        ixp = nint(xp)
        dx = xp - ixp
        ixp = ixp + ibFrom
        ixpP1 = min(ixp + 1, nxProj + ibFrom)
        ixpM1 = max(ixp - 1, 1 + ibFrom)
        v4 = array(ixpM1)
        v5 = array(ixp)
        v6 = array(ixpP1)
        !
        a = (v6 + v4) * .5 - v5
        c = (v6 - v4) * .5
        !
        denNew = a * dx * dx + c * dx + v5
        ! dennew=min(dennew, max(v4, v5, v6))
        ! dennew=max(dennew, min(v4, v5, v6))
        array(ibaseTo + i) = denNew
      enddo
    else
      !
      ! cubic
      !
      do i = 1, nxStretched(iv)
        x = i / float(interpFacStretch) + stretchOffset(iv)
        xp = min(max(1., x * cosBeta(iv)), float(nxProj))
        ixp = xp
        dx = xp - ixp
        ixp = ixp + ibFrom
        ixpP1 = min(ixp + 1, nxProj + ibFrom)
        ixpM1 = max(ixp - 1, 1 + ibFrom)
        ixpP2 = min(ixp + 2, nxProj + ibFrom)

        dxM1 = dx - 1.
        dxDxM1 = dx * dxM1
        fx1 = -dxM1 * dxDxM1
        fx4 = dx * dxDxM1
        fx2 = 1 + dx**2 * (dx - 2.)
        fx3 = dx * (1. -dxDxM1)
        denNew = fx1 * array(ixpM1) + fx2 * array(ixp) + &
            fx3 * array(ixpP1) + fx4 * array(ixpP2)
        ! dennew=min(dennew, max(array(ixpm1), array(ixp), array(ixpp1), &
        ! array(ixpp2)))
        ! dennew=max(dennew, min(array(ixpm1), array(ixp), array(ixpp1), &
        ! array(ixpp2)))
        array(ibaseTo + i) = denNew
      enddo
    endif
  enddo

  return
end subroutine transform


! ---------------------------------------------------------------------
subroutine project(indStart, lslice)
  ! --------------------------
  !
  ! This subroutine assembles one reconstructed slice perpendicular
  ! to the tilt axis, using a back projection method.
  !
  use tiltvars
  implicit none
  integer*4 jStart(3), jEnd(3)
  real*8 xproj8, tstart
  integer*4 nxProjPad, iprojDelta, ipoint, iv, index, i, j
  real*4 cbeta, sbeta, zz, zPart, yy, yproj, yfrac, oneMyFrac
  integer*4 jProj, jLeft, jRight, iproj, ip1, ip2, ind, iprojBase, ifYtest
  integer*4 jTestLeft, jTestRight, indStart, lslice, jregion
  real*4 xLeft, xRight, x, xfrac, oneMxFrac, zBottom, zTop, xproj, yEndTol
  integer*4 gpubpNox, gpubpXtilt, gpubpLocal
  real*8 wallTime
  !
  ! A note on the ubiquitous ytol: It is needed to keep artifacts from
  ! building up at ends of data set through SIRT, from reprojection of the
  ! line between real data and fill, or backprojection of an edge in the
  ! projection difference.  2.05 was sufficient for X-axis tilt cases but
  ! 3.05 was needed for local alignments, thus it is set to 3.05 everywhere
  ! (Here, reprojection routines, and in GPU routines)
  yEndTol = 3.05
  !
  ! Determine maskEdges extent if it is variable
  if (ifAlpha .ne. 0 .and. maskEdges) call maskPrep(lslice)
  nxProjPad = nxProj + 2 + numPad
  tstart = wallTime()
  !
  ! GPU backprojection
  if (useGPU) then
    ind = 1
    if (ifAlpha <= 0 .and. nxWarp == 0) then
      ind = gpubpNox(array(indOutSlice), array(indStart), sinBeta, cosBeta, nxProj, &
          xcenIn + axisXoffset, xcenOut, ycenOut, edgeFill)
    else if (nxWarp == 0 .and. loadGpuStart > 0) then
      ind = gpubpXtilt(array(indOutSlice), sinBeta, cosBeta, sinAlpha, cosAlpha, xzfac, &
          yzfac, nxProj, nyProj, xcenIn + axisXoffset, xcenOut, ycenOut, lslice, &
          centerSlice, edgeFill)
    else if (loadGpuStart > 0) then
      ind = gpubpLocal(array(indOutSlice), lslice, nxWarp, nyWarp, ixStartWarp, &
          iyStartWarp, idelXwarp, idelYwarp, nxProj, xcenOut, xcenIn, axisXoffset, &
          ycenOut, centerSlice, edgeFill)
    endif
    if (ind == 0) then
      if (debug) write(*, '(a,f9.5)') 'GPU backprojection time', &
          wallTime() - tstart
      return
    endif
  endif
  !
  ! CPU backprojection: clear out the slice
  do i = 0, isliceSizeBP - 1
    array(indOutSlice + i) = 0.
  enddo
  iprojDelta = idelSlice * inPlaneSize
  !
  if (nxWarp == 0) then
    !
    ! Loop over all views
    ipoint = indStart - 1
    do iv = 1, numViews
      !
      ! Set view angle
      cbeta = cosBeta(iv)
      sbeta = sinBeta(iv)
      !
      ! Loop over all points in output slice
      index = indOutSlice
      !
      do i = 1, ithickBP
        zz = (i - ycenOut) * compress(iv)
        if (ifAlpha <= 0) then
          zPart = zz * sbeta + xcenIn + axisXoffset
        else
          !
          ! If x-axis tilting, find interpolation factor between the
          ! slices
          !
          yy = lslice - centerSlice
          zPart = yy * sinAlpha(iv) * sbeta + zz * (cosAlpha(iv) * sbeta + xzfac(iv)) + &
              xcenIn + axisXoffset
          yproj = yy * cosAlpha(iv) - zz * (sinAlpha(iv) - yzfac(iv)) + centerSlice
          !
          ! if inside the tolerance, clamp it to the endpoints
          if (yproj >= 1. - yEndTol .and. yproj <= nyProj + yEndTol) &
              yproj = max(1., min(float(nyProj), yproj))
          jProj = yproj
          jProj = min(nyProj - 1, jProj)
          yfrac = yproj - jProj
          oneMyFrac = 1. -yfrac
        endif
        !
        ! compute left and right limits that come from legal data
        !
        x = cbeta
        if (abs(cbeta) < 0.001) x = sign(0.001, cbeta)
        xLeft = (1. -zPart) / x + xcenOut
        xRight = (nxProj - zPart) / x + xcenOut
        if (xRight < xLeft) then
          x = xLeft
          xLeft = xRight
          xRight = x
        endif
        jLeft = xLeft
        if (jLeft < xLeft) jLeft = jLeft + 1
        jLeft = max(jLeft, ixUnmaskedSE(1, i))
        jRight = xRight
        if (jRight == xRight) jRight = jRight - 1
        jRight = min(jRight, ixUnmaskedSE(2, i))
        !
        ! If the limits are now crossed, just skip to full fill at end
        if (jLeft <= jRight) then
          !
          ! set up starting index and projection position
          !
          do ind = index + ixUnmaskedSE(1, i) - 1, index + jLeft - 2
            array(ind) = array(ind) + edgeFill
          enddo
          index = index + (jLeft - 1)
          x = jLeft - xcenOut
          if (interpFacStretch .ne. 0) then
            !
            ! Computation with prestretched data
            !
            xproj8 = interpFacStretch * (zPart / cbeta + x - stretchOffset(iv))
            iproj = xproj8
            xfrac = xproj8 - iproj
            iproj = iproj + ipoint + indStretchLine(iv)
            oneMxFrac = 1. -xfrac
            if (ifAlpha <= 0) then
              !
              ! interpolation in simple case of no x-axis tilt
              !
              do ind = index, index + jRight - jLeft
                array(ind) = array(ind) + &
                    oneMxFrac * array(iproj) + xfrac * array(iproj + 1)
                iproj = iproj + interpFacStretch
              enddo
              index = index + jRight + 1 - jLeft
            else
              !
              ! If x-axis tilting, interpolate from two lines
              !
              ip1 = iproj + (jProj - lslice) * iprojDelta
              ip2 = ip1 + iprojDelta
              if (yproj >= 1. .and. yproj <= nyProj .and. &
                  ip1 >= indLoadBase .and. ip2 >= indLoadBase .and. ip1 < indLoadEnd &
                  .and. ip2 < indLoadEnd) then

                do ind = index, index + jRight - jLeft
                  array(ind) = array(ind) + &
                      oneMxFrac * (oneMyFrac * array(ip1) + yfrac * array(ip2)) + &
                      xfrac * (oneMyFrac * array(ip1 + 1) + yfrac * array(ip2 + 1))
                  ip1 = ip1 + interpFacStretch
                  ip2 = ip2 + interpFacStretch
                enddo
              else
                do ind = index, index + jRight + 1 - jLeft - 1
                  array(ind) = array(ind) + edgeFill
                enddo
              endif
              index = index + jRight + 1 - jLeft
            endif
          else
            !
            ! Computation direct from projection data
            !
            xproj8 = zPart + x * cbeta
            if (ifAlpha <= 0) then
              !
              ! interpolation in simple case of no x-axis tilt
              !
              call bpsumnox(array, index, ipoint, jRight + 1 - jLeft, xproj8, cbeta)
            else
              !
              ! If x-axis tilting
              !
              iproj = xproj8
              iprojBase = ipoint + (jProj - lslice) * iprojDelta
              ip1 = iprojBase + iproj
              ip2 = ip1 + iprojDelta
              if (yproj >= 1. .and. yproj <= nyProj .and. &
                  ip1 >= indLoadBase .and. ip2 >= indLoadBase .and. ip1 < indLoadEnd &
                  .and. ip2 < indLoadEnd) then

                call bpsumxtilt(array, index, iprojBase, iprojDelta, jRight + 1 - jLeft, &
                    xproj8, cbeta, yfrac, oneMyFrac)
              else
                do ind = index, index + jRight + 1 - jLeft - 1
                  array(ind) = array(ind) + edgeFill
                enddo
                index = index + jRight + 1 - jLeft
              endif
            endif
          endif
        else
          jRight = 0
        endif
        do ind = index, index + iwidth - jRight - 1
          array(ind) = array(ind) + edgeFill
        enddo
        index = index + iwidth - jRight
      enddo
      !
      !-------------------------------------------
      !
      ! End of projection loop
      if (interpFacStretch == 0) ipoint = ipoint + nxProjPad
    enddo
  else
    !
    ! LOCAL ALIGNMENTS
    !
    ! Loop over all views
    ipoint = indStart - 1
    do iv = 1, numViews
      !
      ! precompute the factors for getting xproj and yproj all the
      ! way across the slice
      !
      ifYtest = 0
      zBottom = (1 - ycenOut) * compress(iv)
      zTop = (ithickBP - ycenOut) * compress(iv)
      do j = 1, iwidth
        !
        ! get the fixed and z-dependent component of the
        ! projection coordinates

        call localProjFactors(j, lslice, iv, xprojfs(j),  xprojzs(j), &
            yprojfs(j), yprojzs(j))
        !
        ! see if any y testing is needed in the inner loop by checking
        ! yproj at top and bottom in Z
        !
        yproj = yprojfs(j) + yprojzs(j) * zBottom
        jProj = yproj
        ip1 = ipoint + (jProj - lslice) * iprojDelta + 1
        ip2 = ip1 + iprojDelta
        if (ip1 <= indLoadBase .or. ip2 <= indLoadBase .or. ip1 >= indLoadEnd &
            .or. ip2 >= indLoadEnd .or. jProj < 1 .or. jProj >= nyProj) &
            ifYtest = 1
        yproj = yprojfs(j) + yprojzs(j) * zTop
        jProj = yproj
        ip1 = ipoint + (jProj - lslice) * iprojDelta + 1
        ip2 = ip1 + iprojDelta
        if (ip1 <= indLoadBase .or. ip2 <= indLoadBase .or. ip1 >= indLoadEnd &
            .or. ip2 >= indLoadEnd .or. jProj < 1 .or. jProj >= nyProj) &
            ifYtest = 1
      enddo
      !
      ! walk in from each end until xproj is safely within bounds
      ! to define region where no x checking is needed
      !
      jTestLeft = 0
      j = 1
      do while(jTestLeft == 0 .and. j < iwidth)
        if (min(xprojfs(j) + zBottom * xprojzs(j), &
            xprojfs(j) + zTop * xprojzs(j)) >= 1) jTestLeft = j
        j = j + 1
      enddo
      if (jTestLeft == 0) jTestLeft = iwidth
      !
      jTestRight = 0
      j = iwidth
      do while(jTestRight == 0 .and. j > 1)
        if (max(xprojfs(j) + zBottom * xprojzs(j), &
            xprojfs(j) + zTop * xprojzs(j)) < nxProj) jTestRight = j
        j = j - 1
      enddo
      if (jTestRight == 0) jTestRight = 1
      if (jTestRight < jTestLeft) then
        jTestRight = iwidth / 2
        jTestLeft = jTestRight + 1
      endif
      !
      index = indOutSlice
      !
      ! loop over the slice, outer loop on z levels
      !
      do i = 1, ithickBP
        zz = (i - ycenOut) * compress(iv)
        jLeft = max(jTestLeft, ixUnmaskedSE(1, i))
        jRight = min(jTestRight, ixUnmaskedSE(2, i))
        index = index + ixUnmaskedSE(1, i) - 1
        !
        ! set up to do inner loop in three regions of X
        !
        jStart(1) = ixUnmaskedSE(1, i)
        jEnd(1) = jLeft - 1
        jStart(2) = jLeft
        jEnd(2) = jRight
        jStart(3) = jRight + 1
        jEnd(3) = ixUnmaskedSE(2, i)
        do jregion = 1, 3
          if (jregion .ne. 2 .or. ifYtest == 1) then
            !
            ! loop involving full testing - either left or right
            ! sides needing x testing, or anywhere if y testing
            ! needed
            !
            do j = jStart(jregion), jEnd(jregion)
              xproj = xprojfs(j) + zz * xprojzs(j)
              yproj = yprojfs(j) + zz * yprojzs(j)
              if (yproj >= 1. - yEndTol .and. yproj <= nyProj + yEndTol) &
                  yproj = max(1., min(float(nyProj), yproj))
              if (xproj >= 1 .and. xproj <= nxProj .and. &
                  yproj >= 1. .and. yproj <= nyProj) then
                !
                iproj = xproj
                iproj = min(nxProj - 1, iproj)
                xfrac = xproj - iproj
                jProj = yproj
                jProj = min(nyProj - 1, jProj)
                yfrac = yproj - jProj
                !
                ip1 = ipoint + (jProj - lslice) * iprojDelta + iproj
                ip2 = ip1 + iprojDelta
                if (ip1 >= indLoadBase .and. ip2 >= indLoadBase .and. &
                    ip1 < indLoadEnd .and. ip2 < indLoadEnd) then
                  array(index) = array(index) + &
                      (1. -yfrac) * ((1. -xfrac) * array(ip1) &
                      + xfrac * array(ip1 + 1)) + &
                      yfrac * ((1. -xfrac) * array(ip2) &
                      + xfrac * array(ip2 + 1))
                else
                  array(index) = array(index) + edgeFill
                endif
              else
                array(index) = array(index) + edgeFill
              endif
              index = index + 1
            enddo
            !
            ! loop for no x-testing and no y testing
            !
          else
            call bpsumLocal(array, index, zz, xprojfs, xprojzs, yprojfs, &
                yprojzs, ipoint, iprojDelta, lslice, jStart(jregion), &
                jEnd(jregion))
          endif
        enddo
        index = index + iwidth - ixUnmaskedSE(2, i)
      enddo
      !-------------------------------------------
      !
      ! End of projection loop
      ipoint = ipoint + nxProjPad
    enddo
  endif
  if (debug) write(*, '(a,f9.5)') 'CPU backprojection time', &
      wallTime() - tstart
  return
end subroutine project
!
!-------------------------------------------------------------------------
!
! COMPOSE will interpolate the output slice LSLICEOUT from vertical
! slices in the ring buffer, where lvertSliceStart and lVertSliceEnd are the starting
! and ending slices in the ring buffer, IDIR is the direction of
! reconstruction, and IRINGSTART is the position of lvertSliceStart in the
! ring buffer.
!
subroutine compose(lsliceOut, lvertSliceStart, lVertSliceEnd, idir, iringStart, &
    composeFill)
  use tiltvars
  implicit none
  integer*4 lsliceOut, lvertSliceStart, lVertSliceEnd, idir, iringStart
  real*4 composeFill
  integer*4 ind1(4), ind2(4), ind3(4), ind4(4)
  real*4 tanAlpha, vertCen, centeredJ, centeredL, vSlice, vyCentered, fx, vy, fy, f22, f23
  real*4 f32, f33
  integer*4 ivSlice, ifMiss, i, lvSlice, iring, ibaseInd, ivy, indCen, jnd5, jnd2, j, k
  real*4 fx1, fx2, fx3, fx4, fy1, fy2, fy3, fy4, v1, v2, v3, v4, f5, f2, f8, f4, f6
  integer*4 jnd8, jnd4, jnd6, numFill
  !
  ! 12/12/09: stopped reading base here, read on output; eliminate zeroing
  !
  tanAlpha = sinAlpha(1) / cosAlpha(1)
  vertCen = ithickBP / 2 + 0.5
  fx = 0.
  fy = 0.
  numFill = 0
  !
  ! loop on lines of data
  !
  do j = 1, ithickOut
    centeredJ = j - (ithickOut / 2 + 0.5) - yOffset
    centeredL = lsliceOut - centerSlice
    !
    ! calculate slice number and y position in vertical slices
    !
    vSlice = centeredL * cosAlpha(1) - centeredJ * sinAlpha(1) + centerSlice
    vyCentered = centeredL * sinAlpha(1) + centeredJ * cosAlpha(1)
    ivSlice = vSlice
    fx = vSlice - ivSlice
    ifMiss = 0
    !
    ! for each of 4 slices needed for cubic interpolation, initialize
    ! data indexes at zero then see if slice exists in ring
    !
    do i = 1, 4
      ind1(i) = 0
      ind2(i) = 0
      ind3(i) = 0
      ind4(i) = 0
      lvSlice = ivSlice + i - 2
      if (idir * (lvSlice - lvertSliceStart) >= 0 .and.  &
          idir * (lVertSliceEnd - lvSlice) >= 0) then
        !
        ! if slice exists, get base index for the slice, compute the
        ! y index in the slice, then set the 4 data indexes if they
        ! are within the slice
        !
        iring = idir * (lvSlice - lvertSliceStart) + iringStart
        if (iring > numVertNeeded) iring = iring - numVertNeeded
        ibaseInd = indOutSlice + isliceSizeBP + (iring - 1) * ithickBP * iwidth
        vy = vyCentered + vertCen - nint(tanAlpha * (lvSlice - centerSlice)) +  &
            yOffset / cosAlpha(1)
        ivy = vy
        fy = vy - ivy
        if (ivy - 1 >= 1 .and. ivy - 1 <= ithickBP) &
            ind1(i) = ibaseInd + iwidth * (ivy - 2)
        if (ivy >= 1 .and. ivy <= ithickBP) ind2(i) = ibaseInd + iwidth * (ivy - 1)
        if (ivy + 1 >= 1 .and. ivy + 1 <= ithickBP) ind3(i) = ibaseInd + iwidth * ivy
        if (ivy + 2 >= 1 .and. ivy + 2 <= ithickBP) &
            ind4(i) = ibaseInd + iwidth * (ivy + 1)
      endif
      if (ind1(i) == 0 .or. ind2(i) == 0 .or. ind3(i) == 0 .or. &
          ind4(i) == 0) ifMiss = 1
    enddo
    ibaseInd = indOutSlice + (j - 1) * iwidth - 1
    if (interpOrdXtilt > 2 .and. ifMiss == 0) then
      !
      ! cubic interpolation if selected, and no data missing
      !
      fx1 = 2. *fx**2 - fx**3 - fx
      fx2 = fx**3 - 2. *fx**2 + 1
      fx3 = fx**2 + fx - fx**3
      fx4 = fx**3 - fx**2
      fy1 = 2. *fy**2 - fy**3 - fy
      fy2 = fy**3 - 2. *fy**2 + 1
      fy3 = fy**2 + fy - fy**3
      fy4 = fy**3 - fy**2
      do i = 1, iwidth
        v1 = fx1 * array(ind1(1)) + fx2 * array(ind1(2)) + &
            fx3 * array(ind1(3)) + fx4 * array(ind1(4))
        v2 = fx1 * array(ind2(1)) + fx2 * array(ind2(2)) + &
            fx3 * array(ind2(3)) + fx4 * array(ind2(4))
        v3 = fx1 * array(ind3(1)) + fx2 * array(ind3(2)) + &
            fx3 * array(ind3(3)) + fx4 * array(ind3(4))
        v4 = fx1 * array(ind4(1)) + fx2 * array(ind4(2)) + &
            fx3 * array(ind4(3)) + fx4 * array(ind4(4))
        array(ibaseInd + i) = fy1 * v1 + fy2 * v2 + fy3 * v3 + fy4 * v4
        do k = 1, 4
          ind1(k) = ind1(k) + 1
          ind2(k) = ind2(k) + 1
          ind3(k) = ind3(k) + 1
          ind4(k) = ind4(k) + 1
        enddo
      enddo
    elseif (interpOrdXtilt == 2 .and. ifMiss == 0) then
      !
      ! quadratic interpolation if selected, and no data missing
      ! shift to next column or row if fractions > 0.5
      !
      indCen = 2
      if (fx > 0.5) then
        indCen = 3
        fx = fx - 1.
      endif
      if (fy <= 0.5) then
        jnd5 = ind2(indCen)
        jnd2 = ind1(indCen)
        jnd8 = ind3(indCen)
        jnd4 = ind2(indCen - 1)
        jnd6 = ind2(indCen + 1)
      else
        fy = fy - 1.
        jnd5 = ind3(indCen)
        jnd2 = ind2(indCen)
        jnd8 = ind4(indCen)
        jnd4 = ind3(indCen - 1)
        jnd6 = ind3(indCen + 1)
      endif
      !
      ! get coefficients and do the interpolation
      !
      f5 = 1. -fx**2 - fy**2
      f2 = (fy**2 - fy) / 2.
      f8 = f2 + fy
      f4 = (fx**2 - fx) / 2.
      f6 = f4 + fx
      do i = 1, iwidth
        array(ibaseInd + i) = f5 * array(jnd5) + f2 * array(jnd2) + &
            f4 * array(jnd4) + f6 * array(jnd6) + f8 * array(jnd8)
        jnd5 = jnd5 + 1
        jnd2 = jnd2 + 1
        jnd4 = jnd4 + 1
        jnd6 = jnd6 + 1
        jnd8 = jnd8 + 1
      enddo
    else
      !
      ! linear interpolation
      !
      ! print *,j, ind2(2), ind2(3), ind3(2), ind3(3)
      if (ind2(2) == 0 .or. ind2(3) == 0 .or. ind3(2) == 0 .or. &
          ind3(3) == 0) then
        !
        ! if there is a problem, see if it can be rescued by shifting
        ! center back to left or below
        !
        if (fx < 0.02 .and. ind2(1) .ne. 0 .and. ind3(1) .ne. 0 .and. &
            ind2(2) .ne. 0 .and. ind3(2) .ne. 0) then
          fx = fx + 1
          ind2(3) = ind2(2)
          ind2(2) = ind2(1)
          ind3(3) = ind3(2)
          ind3(2) = ind3(1)
        elseif (fy < 0.02 .and. ind1(2) .ne. 0 .and. ind1(3) .ne. 0 .and. &
            ind2(2) .ne. 0 .and. ind3(2) .ne. 0) then
          fy = fy + 1
          ind3(2) = ind2(2)
          ind2(2) = ind1(2)
          ind3(3) = ind2(3)
          ind2(3) = ind1(3)
        endif
      endif
      !
      ! do linear interpolation if conditions are right, otherwise fill
      !
      if (ind2(2) .ne. 0 .and. ind2(3) .ne. 0 .and. ind3(2) .ne. 0 .and. &
          ind3(3) .ne. 0) then
        f22 = (1. -fy) * (1. -fx)
        f23 = (1. -fy) * fx
        f32 = fy * (1. -fx)
        f33 = fy * fx
        do i = 1, iwidth
          array(ibaseInd + i) = f22 * array(ind2(2)) + &
              f23 * array(ind2(3)) + f32 * array(ind3(2)) + f33 * array(ind3(3))
          ind2(2) = ind2(2) + 1
          ind2(3) = ind2(3) + 1
          ind3(2) = ind3(2) + 1
          ind3(3) = ind3(3) + 1
        enddo
      else
        ! print *,'filling', j
        do i = 1, iwidth
          array(i + ibaseInd) = composeFill
        enddo
        numFill = numFill + 1
      endif
    endif
  enddo
  ! if (nfill .ne. 0) print *,nfill, ' lines filled, edgefill =', composeFill
  return
end subroutine compose

!
! DECOMPOSE will interpolate a vertical slice LSLICE from input slices
! in the read-in ring buffer, where lReadStart and lReadEnd are the
! starting and ending slices in the ring buffer, IRINGSTART is the
! position of lReadStart in the ring buffer, and ibaseSIRT is the index
! in stack at which to place the slice.
!
subroutine decompose(lslice, lReadStart, lreadEnd, iringStart, ibaseSIRT)
  use tiltvars
  implicit none
  integer*4 lslice, lReadStart, lreadEnd, iringStart, ibaseSIRT
  real*4 tanAlpha, outCen, centeredL, vSliceCentered, vyCentered, outSlice, outYpos
  real*4 f11, f12, f21, f22, fx, fy
  integer*4 ibaseVert, ibase1, ibase2, j, ioutSlice, jOut, i, iring

  tanAlpha = sinAlpha(1) / cosAlpha(1)
  outCen = ithickOut / 2 + 0.5
  centeredL = lslice - centerSlice
  vSliceCentered = centeredL
  !
  ! loop on lines of data
  do j = 1, ithickBP
    ibaseVert = ibaseSIRT + (j - 1) * iwidth
    !
    ! calculate slice number and y position in input slices
    !
    vyCentered = j - (ithickBP / 2 + 0.5 - nint(tanAlpha * centeredL) +  &
        yOffset / cosAlpha(1))
    outSlice = centerSlice + vSliceCentered * cosAlpha(1) + vyCentered * sinAlpha(1)
    outYpos = outCen + yOffset - vSliceCentered * sinAlpha(1) + vyCentered * cosAlpha(1)
    ! print *,j, vycen, outsl, outj
    ! if (outsl >= lreadStart - 0.5 .and. outsl <= lreadEnd + 0.5 &
    ! .and. outj >= 0.5 .and.outj <= ithickOut + 0.5) then
    !
    ! For a legal position, get interpolation integers and fractions,
    ! adjust if within half pixel of end
    ioutSlice = outSlice
    fx = outSlice - ioutSlice
    if (ioutSlice < lReadStart) then
      ioutSlice = lReadStart
      fx = 0.
    else if (ioutSlice > lreadEnd - 1) then
      ioutSlice = lreadEnd - 1
      fx = 1.
    endif
    jOut = outYpos
    fy = outYpos - jOut
    if (jOut < 1) then
      jOut = 1
      fy = 0.
    else if (jOut > ithickOut - 1) then
      jOut = ithickOut - 1
      fy = 1.
    endif
    !
    ! Get slice indexes in ring
    iring = ioutSlice - lReadStart + iringStart
    if (iring > numReadNeed) iring = iring - numReadNeed
    ibase1 = ireadBase + (iring - 1) * ithickOut * iwidth + &
        (jOut - 1) * iwidth
    iring = ioutSlice + 1 - lReadStart + iringStart
    if (iring > numReadNeed) iring = iring - numReadNeed
    ibase2 = ireadBase + (iring - 1) * ithickOut * iwidth + &
        (jOut - 1) * iwidth
    !
    ! Interpolate line
    f11 = (1. -fy) * (1. -fx)
    f12 = (1. -fy) * fx
    f21 = fy * (1. -fx)
    f22 = fy * fx
    do i = 0, iwidth - 1
      array(ibaseVert + i) = f11 * array(ibase1 + i) + f12 * array(ibase2 + i) + &
          f21 * array(ibase1 + i + iwidth) + f22 * array(ibase2 + i + iwidth)
    enddo
    ! else
    !
    ! Otherwise fill line
    ! do i = 0, iwidth - 1
    ! array(ibasev+i) = dmeanIn
    ! enddo
    ! endif
  enddo
  return
end subroutine decompose

!
!-------------------------------------------------------------------------
subroutine dumpSlice(lslice, dmin, dmax, dtot8)
  ! --------------------------------------
  !
  use tiltvars
  implicit none
  integer*4 lslice, numParExtraLines, iend, index, i, j, imapOut, indDel
  real*4 dmin, dmax, fill
  real*8 dtot8, dtmp8
  !
  ! If adding to a base rec, read in each line and add scaled values
  imapOut = indOutSlice
  if (numSIRTiter > 0 .and. ifAlpha >= 0) imapOut = ireadBase
  if (readBaseRec) then
    index = imapOut
    call iiuSetPosition(3, lslice - 1, 0)
    if (iterForReport > 0) call sampleForReport(array(imapOut), &
        lslice, ithickOut, 1, outScale, outAdd)
    do j = 1, ithickOut
      call irdlin(3, projLine)
      if (recSubtraction) then
        !
        ! SIRT subtraction from base with possible sign constraints
        if (isignConstraint == 0) then
          array(index:index + iwidth - 1) = projLine(1:iwidth) / outScale - outAdd - &
              array(index:index + iwidth - 1)
        else if (isignConstraint < 0) then
          array(index:index + iwidth - 1) = min(0., projLine(1:iwidth) / outScale - &
              outAdd - array(index:index + iwidth - 1))
        else
          array(index:index + iwidth - 1) = max(0., projLine(1:iwidth) / outScale - &
              outAdd - array(index:index + iwidth - 1))
        endif
      else
        !
        ! Generic addition of the scaled base data
        array(index:index + iwidth - 1) = array(index:index + iwidth - 1) + &
            projLine(1:iwidth) / baseOutScale - baseOutAdd
      endif
      index = index + iwidth
    enddo
  endif
  !
  numParExtraLines = 100
  iend = imapOut + ithickOut * iwidth - 1
  !
  ! outScale
  ! DNM simplified and fixed bug in getting min/max/mean
  dtmp8 = 0.
  !
  ! DNM 9/23/04: incorporate reprojBP option
  !
  if (reprojBP) then
    !--------------outScale
    do i = imapOut, iend
      array(i) = (array(i) + outAdd) * outScale
    enddo
    !
    ! Fill value assumes edge fill value
    !
    fill = (edgeFill + outAdd) * outScale
    do j = 1, numReproj
      i = (j - 1) * iwidth + 1
      call reproject(array(imapOut), iwidth, ithickBP, iwidth, sinReproj(j), &
          cosReproj(j), xRayStart(i), yRayStart(i), &
          numPixInRay(i), maxRayPixels(j), fill, projLine, 0, 0)
      do i = 1, iwidth
        dmin = amin1(projLine(i), dmin)
        dmax = amax1(projLine(i), dmax)
        dtmp8 = dtmp8 + projLine(i)
      enddo
      i = (lslice - isliceStart) / idelSlice
      if (minTotSlice > 0) i = lslice - minTotSlice
      call parWrtPosn(2, j - 1, i)
      call parWrtLin(2, projLine)
    enddo
    dtot8 = dtot8 + dtmp8
    return
  endif
  !
  !--------------outScale and get min / max / sum
  do i = imapOut, iend
    array(i) = (array(i) + outAdd) * outScale
    dmin = amin1(array(i), dmin)
    dmax = amax1(array(i), dmax)
    dtmp8 = dtmp8 + array(i)
    !
  enddo
  dtot8 = dtot8 + dtmp8
  !
  ! Dump slice
  if (perpendicular) then
    ! ....slices correspond to sections of map
    call parWrtSec(2, array(imapOut))
  else
    ! ....slices must be properly stored
    ! Take each line of array and place it in the correct section of the map.
    index = imapOut
    indDel = 1
    if (rotateBy90) then
      index = imapOut + (ithickOut - 1) * iwidth
      indDel = -1
    endif
    do j = 1, ithickOut
      call iiuSetPosition(2, j - 1, (lslice - isliceStart) / idelSlice)
      call iwrlin(2, array(index))
      !
      ! DNM 2/29/01: partially demangle the parallel output by writing
      ! up to 100 lines at a time in this plane
      !
      if (mod((lslice - isliceStart) / idelSlice, numParExtraLines) == 0) then
        do i = 1, min(numParExtraLines - 1, (isliceEnd - lslice) / idelSlice)
          call iwrlin(2, array(index))
        enddo
      endif
      index = index + indDel * iwidth
    enddo
  endif
  !
  return
end subroutine dumpSlice

! maskEdges out (blur) the edges of the slice at the given index
!
subroutine maskSlice(ibaseMask, ithickMask)
  use tiltvars
  implicit none
  integer*4 ibaseMask, ithickMask, i, j, idir, limit, lr, nsum, numSmooth
  integer*4 numTaper, index, ibase
  real*4 sum, edgeMean, frac
  !
  numSmooth = 10
  numTaper = 10
  idir = -1
  limit = 1
  do lr = 1, 2
    !
    ! find mean along edge
    sum = 0.
    do j = 1, ithickMask
      sum = sum + array(ibaseMask + (j - 1) * iwidth + ixUnmaskedSE(lr, j) - 1)
    enddo
    edgeMean = sum / ithickMask
    !
    ! For each line, sum progressively more pixels along edge out to a
    ! limit
    do j = 1, ithickMask
      ibase = ibaseMask + (j - 1) * iwidth - 1
      index = ixUnmaskedSE(lr, j) + idir
      sum = array(ibase + ixUnmaskedSE(lr, j))
      nsum = 1
      do i = 1, numSmooth
        if (idir * (limit - index) < 0) exit
        if (j + i <= ithickMask) then
          nsum = nsum + 1
          sum = sum + array(ibase + i * iwidth + ixUnmaskedSE(lr, j + i))
        endif
        if (j - i >= 1) then
          nsum = nsum + 1
          sum = sum + array(ibase - i * iwidth + ixUnmaskedSE(lr, j - i))
        endif
        !
        ! Taper partway to mean over this smoothing distance
        frac = i / (numTaper + numSmooth + 1.)
        array(ibase + index) = (1. - frac) * sum / nsum + frac * edgeMean
        index = index + idir
      enddo
      !
      ! Then taper rest of way down to mean over more pixels, then fill
      ! with mean
      do i = 1, numTaper
        if (idir * (limit - index) < 0) exit
        frac = (i + numSmooth) / (numTaper + numSmooth + 1.)
        array(ibase + index) = (1. - frac) * sum / nsum + frac * edgeMean
        index = index + idir
      enddo
      do i = index, limit, idir
        array(ibase + i) = edgeMean
      enddo
    enddo
    !
    ! Set up for other direction
    limit = iwidth
    idir = 1
  enddo
  return
end subroutine maskSlice

! The input routine, gets input and sets up almost all parameters
! -----------------------------------------------------------------
subroutine inputParameters()
  ! ----------------

  use tiltvars
  implicit none
  integer LIMNUM
  parameter (LIMNUM = 100)
  !
  integer*4 mpxyz(3), noutXyz(3), nrecXyz(3), nvsXyz(3), maxNeeds(LIMNUM)
  integer*4 maxTex2D(2), maxTexLayer(3), maxTex3D(3)
  real*4 outHdrTilt(3), cell(6), degToRad
  data outHdrTilt/90., 0., 0./
  data cell/0., 0., 0., 90., 90., 90./
  data degToRad/0.0174533/
  character dat*9, tim*8
  real*4 delta(3)
  character*80 titlech
  !
  character*1024 card
  character*320 inputFile, outputFile, recFile, baseFile, boundFile, vertBoundFile
  character*320 vertOutFile, angleOutput, transformFile, defocusFile, gpuErrorStr
  integer*4 nfields, inum(LIMNUM)
  real*4 xnum(LIMNUM)
  !
  integer*4, allocatable :: ivExclude(:), ivReproj(:)
  real*4, allocatable :: packLocal(:,:), angReproj(:)
  integer*4 mode, newAngles, ifTiltFile, numViewUse, numViewExclude, numNeedEval
  real*4 delAngle, compFactor, globalAlpha, xOffset, scaleLocal, radMax, radFall, xoffAdj
  integer*4 iradFall, iradMax, numCompress, nxfull, nyFull, ixSubset, iySubset, ixSubsetIn
  integer*4 kti, indBase, ipos, idType, lens, nxFullIn
  integer*4 nd1, nd2, nv, numSlices, indi, i, iex, numViewOrig, iv
  real*4 vd1, vd2, delTheta, theta, thetaView
  integer*4 nxProjPad, numInput, numExcludeList, j, ind, numEnt
  integer*4 numPadTmp, nprojPad, ifZfactors, localZfacs
  integer*4 ifThickIn, ifSliceIn, ifWidthIn, imageBinned, ifSubsetIn, ierr
  real*4 pixelLocal, dminTmp, dmaxTmp, dmeanTmp, frac, originX, originY, originZ
  real*4 gpuMemoryFrac, gpuMemory, pixForDefocus, focusInvert
  integer*4 nViewsReproj, iwideReproj, k, ind1, ind2, ifExpWeight
  integer*4 minMemory, indGPU, iactGpuFailOption, iactGpuFailEnviron
  integer*4 ifGpuByEnviron, indDelta, ifExit, ifMultByGaussian, ifHammingLike, if3Dtexture
  logical*4 adjustOrigin, projModel, readw_or_imod
  integer*4 niceFrame, parWrtInitialize, gpuAvailable, imodGetEnv, parWrtSetCurrent
  integer*4 gpuAllocArrays, allocateArray, gpuLoadLocals, gpuLoadFilter, niceFFTlimit
  integer*4 b3dLockFile, b3dOutputFileType, iiuFileType, iiTestIfHDF
  real*8 wallTime, wallStart
  !
  integer*4 numOptArg, numNonOptArg
  integer*4 PipGetInteger, PipGetBoolean, PipGetLogical, PipGetTwoFloats
  integer*4 PipGetString, PipGetFloat, PipGetTwoIntegers, PipGetFloatArray
  integer*4 PipGetInOutFile, PipGetIntegerArray, PipNumberOfEntries
  integer*4 PipGetThreeFloats
  !
  ! fallbacks from ../../manpages/autodoc2man -3 2  tilt
  !
  integer numOptions
  parameter (numOptions = 76)
  character*(40 * numOptions) options(1)
  options(1) = &
      'input:InputProjections:FN:@output:OutputFile:FN:@:TILTFILE:FN:@'// &
      ':XTILTFILE:FN:@:ZFACTORFILE:FN:@:LOCALFILE:FN:@:BoundaryInfoFile:FN:@'// &
      ':WeightAngleFile:FN:@:WeightFile:FN:@:WIDTH:I:@:SLICE:IA:@:TOTALSLICES:IP:@'// &
      ':THICKNESS:I:@:OFFSET:FA:@:SHIFT:FA:@:ANGLES:FAM:@:XAXISTILT:F:@'// &
      ':COMPFRACTION:F:@:COMPRESS:FAM:@:FULLIMAGE:IP:@:SUBSETSTART:IP:@'// &
      ':IMAGEBINNED:I:@:LOCALSCALE:F:@:LOG:F:@:RADIAL:FP:@:FalloffIsTrueSigma:B:@'// &
      ':MultiplyByGaussian:B:@:HammingLikeFilter:F:@:DENSWEIGHT:FA:@'// &
      ':ExactFilterSize:I:@:FakeSIRTiterations:I:@:INCLUDE:LIM:@:EXCLUDELIST2:LIM:@'// &
      ':COSINTERP:IA:@:XTILTINTERP:I:@:UseGPU:I:@:ActionIfGPUFails:IP:@:MODE:I:@'// &
      ':SCALE:FP:@:MASK:I:@:PERPENDICULAR:B:@:PARALLEL:B:@:RotateBy90:B:@'// &
      ':AdjustOrigin:B:@:TITLE:CH:@:BaseRecFile:FN:@:BaseNumViews:I:@'// &
      ':SubtractFromBase:LI:@:MinMaxMean:IT:@:REPROJECT:FAM:@:ViewsToReproject:LI:@'// &
      'recfile:RecFileToReproject:FN:@xminmax:XMinAndMaxReproj:IP:@'// &
      'yminmax:YMinAndMaxReproj:IP:@zminmax:ZMinAndMaxReproj:IP:@'// &
      'threshold:ThresholdedReproj:FT:@:FlatFilterFraction:F:@:SIRTIterations:I:@'// &
      ':SIRTSubtraction:B:@:StartingIteration:I:@:VertBoundaryFile:FN:@'// &
      ':VertSliceOutputFile:FN:@:VertForSIRTInput:B:@:ConstrainSign:I:@'// &
      ':ProjectModel:FN:@:AngleOutputFile:FN:@:AlignTransformFile:FN:@'// &
      ':DefocusFile:FN:@:PixelForDefocus:FP:@:DONE:B:@:FBPINTERP:I:@:REPLICATE:FPM:@'// &
      'debug:DebugOutput:B:@internal:InternalSIRTSlices:IP:@param:ParameterFile:PF:@'// &
      'help:usage:B:'
  !
  recReproj = .false.
  nViewsReproj = 0
  useGPU = .false.
  indGPU = -1;
  numGpuPlanes = 0
  numSIRTiter = 0
  sirtFromZero = .false.
  recSubtraction = .false.
  projSubtraction = .false.
  vertSirtInput = .false.
  saveVertSlices = .false.
  angleOutput = ' '
  flatFrac = 0.
  iterForReport = 0
  threshPolarity = 0.
  focusInvert = 0.
  defocusFile = ' '
  gpuErrorStr = ' '
  !
  ! Minimum array size to allocate, desired number of slices to allocate
  ! for if it exceeds that minimum size
  minMemory = 20000000
  numNeedEval = 10
  gpuMemoryFrac = 0.8
  !
  ! Pip startup: set error, parse options, check help, set flag if used
  !
  call PipSetSpecialFlags(1, 1, 1, 2, 0)
  call PipReadOrParseOptions(options, numOptions, 'tilt', &
      'ERROR: TILT - ', .false., 0, 1, 1, numOptArg, numNonOptArg)

  if (PipGetInOutFile('InputProjections', 1, ' ', inputFile) .ne. 0) &
      call exitError('NO INPUT FILE WITH PROJECTIONS SPECIFIED')
  if (PipGetInOutFile('OutputFile', 2, ' ', outputFile) .ne. 0) &
      call exitError('NO OUTPUT FILE SPECIFIED')
  !
  ! Allocate array a little bit for temp use
  limReproj = 100000
  allocate(array(limReproj), angReproj(limReproj), stat = ierr)
  call memoryErrorUC(ierr, 'SMALL TEMPORARY ARRAYS')
  ierr = PipGetLogical('DebugOutput', debug)
  !
  ! Open input projection file
  call imOpen(1, inputFile, 'RO')
  call irdhdr(1, nprojXyz, mpxyz, mode, dminIn, dmaxIn, dmeanIn)
  call iiuRetDataType(1, idType, lens, nd1, nd2, vd1, vd2)
  numViews = nprojXyz(3)
  limView = numViews + 10
  allocate(sinBeta(limView), cosBeta(limView), sinAlpha(limView), cosAlpha(limView),  &
      alpha(limView), angles(limView), xzfac(limView), yzfac(limView),  &
      compress(limView), nxStretched(limView), indStretchLine(limView),  &
      stretchOffset(limView), exposeWeight(limView), mapUsedView(limView), &
      iviewSubtract(limView), wgtAngles(limView), ivExclude(limView), &
      ivReproj(limView), stat = ierr)
  call memoryErrorUC(ierr, 'ARRAYS FOR VARIABLES PER VIEW')
  !
  ! The approximate implicit scaling caused by the default radial filter
  filterScale = nprojXyz(1) / 2.2
  !
  newAngles = 0
  ifTiltFile = 0
  !
  ! Get model file to project and other optional files
  projModel = PipGetString('ProjectModel', recFile) == 0
  if (projModel) then
    if (.not.readw_or_imod(recFile)) call exitError('READING MODEL FILE TO REPROJECT')
    ierr = PipGetString('AngleOutputFile', angleOutput)
    ierr = PipGetString('AlignTransformFile', transformFile)
    ierr = PipGetString('DefocusFile', defocusFile)
    ierr = PipGetTwoFloats('PixelForDefocus', pixForDefocus, focusInvert)
  endif
  !
  ! Get entries for reprojection from rec file
  ierr = PipGetInteger('SIRTIterations', numSIRTiter)
  ierr = PipGetThreeFloats('ThresholdedReproj', threshForReproj, threshPolarity, &
      threshSumFac)
  if (PipGetString('RecFileToReproject', recFile) == 0) then
    if (projModel) call exitError( &
        'YOU CANNOT USE -RecFileToReproject with -ProjectModel')
    projMean = dmeanIn
    call imOpen(3, recFile, 'RO')
    call irdhdr(3, nrecXyz, mpxyz, mode, dminIn, dmaxIn, dmeanIn)
    if (numSIRTiter <= 0) then
      recReproj = .true.
      minXreproj = 0
      minYreproj = 0
      minZreproj = 0
      maxXreproj = nrecXyz(1) - 1
      maxYreproj = nrecXyz(2) - 1
      maxZreproj = nrecXyz(3) - 1
      ierr = PipGetTwoIntegers('XMinAndMaxReproj', minXreproj, maxXreproj)
      ierr = PipGetTwoIntegers('YMinAndMaxReproj', minYreproj, maxYreproj)
      ierr = PipGetTwoIntegers('ZMinAndMaxReproj', minZreproj, maxZreproj)
      if (minXreproj < 0 .or. minYreproj < 0 .or. maxXreproj >= &
          nrecXyz(1) .or. maxYreproj >= nrecXyz(2)) call exitError( &
          'Min or Max X, or Y coordinate to project is out of range')
      if (PipGetString('ViewsToReproj', card) == 0) then
        call parseList(card, ivReproj, nViewsReproj)
        if (nViewsReproj > limView) call exitError( &
            'TOO MANY VIEWS IN LIST TO REPROJECT FOR ARRAYS')
      endif
      ierr = PipGetLogical('SIRTSubtraction', projSubtraction)
    endif
    minXreproj = minXreproj + 1
    maxXreproj = maxXreproj + 1
    minYreproj = minYreproj + 1
    maxYreproj = maxYreproj + 1
    minZreproj = minZreproj + 1
    maxZreproj = maxZreproj + 1
    !
    ! If not reading from a rec and doing sirt, then must be doing from 0
  else if (numSIRTiter > 0) then
    sirtFromZero = .true.
    flatFrac = 1.
  endif
  if (threshPolarity .ne. 0. .and. .not. recReproj) call exitError( &
      'THRESHOLDED REPROJECTION CAN BE USED ONLY WITH THE -RecFileToReproject OPTION')
  if (threshPolarity .ne. 0. .and. (numSIRTiter > 0 .or. projSubtraction)) &
      call exitError('THRESHOLDED REPROJECTION CANNOT BE USED WITH -SIRTSubtraction '// &
      'OR -SIRTIterations')

  if (.not. recReproj .and. .not.projModel .and. numSIRTiter <= 0 .and. &
      PipGetString('BaseRecFile', baseFile) == 0) then
    readBaseRec = .true.
    if (PipGetString('SubtractFromBase', card) == 0) then
      call parseList(card, iviewSubtract, numViewSubtract)
      if (numViewSubtract > limView) call exitError( &
          'TOO MANY VIEWS IN LIST TO SUBTRACT FOR ARRAYS')
      if (numViewSubtract == 1 .and. iviewSubtract(1) < 0) then
        recSubtraction = .true.
        numViewSubtract = 0
      endif
    endif
    if (.not. recSubtraction .and. &
        PipGetInteger('BaseNumViews', numViewBase) .ne. 0) call exitError( &
        'YOU MUST ENTER -BaseNumViews with -BaseRecFile')
    call imOpen(3, baseFile, 'RO')
    call irdhdr(3, nrecXyz, mpxyz, mode, dminTmp, dmaxTmp, dmeanTmp)
  endif
  !
  !-------------------------------------------------------------
  ! Set up defaults:
  !
  !...... Default slice is all rows in a projection plane
  isliceStart = 1
  isliceEnd = nprojXyz(2)
  idelSlice = 1
  !...... and has the same number of columns.
  iwidth = nprojXyz(1)
  !...... Default is no maskEdges and extra pixels to maskEdges is 0
  maskEdges = .false.
  numExtraMaskPix = 0
  !...... Default is no scaling of output map
  outAdd = 0.
  outScale = 1.
  !...... Default is no offset or rotation
  delAngle = 0.
  !...... Default is output mode 2
  newMode = 2
  !...... Start with no list of views to use or exclude
  numViewUse = 0
  numViewExclude = 0
  !...... Default is no logarithms
  ifLog = 0
  !...... Default radial weighting parameters - no filtering
  iradMax = nprojXyz(1) / 2 + 1
  radFall = 0.
  ifMultByGaussian = 0
  !...... Default overall and individual compression of 1; no alpha tilt
  numCompress = 0
  compFactor = 1.
  do nv = 1, numViews
    compress(nv) = 1.
    alpha(nv) = 0.
    xzfac(nv) = 0.
    yzfac(nv) = 0.
    exposeWeight(nv) = 1.
  enddo
  ifAlpha = 0
  globalAlpha = 0.
  ifZfactors = 0
  !
  !...... Default weighting by density of adjacent views
  numTiltIncWgt = 2
  do i = 1, numTiltIncWgt
    tiltIncWgts(i) = 1. / (i - 0.5)
  enddo
  numWgtAngles = 0
  !
  xOffset = 0
  yOffset = 0
  axisXoffset = 0.
  nxWarp = 0
  nyWarp = 0
  scaleLocal = 0.
  nxFullIn = 0
  nyFull = 0
  ixSubsetIn = 0
  iySubset = 0
  ithickBP = 10
  imageBinned = 1
  ifThickIn = 0
  ifWidthIn = 0
  ifSliceIn = 0
  ifSubsetIn = 0
  ifExpWeight = 0
  minTotSlice = -1
  maxTotSlice = -1
  exactSamples = 100.
  numExactCycles = 0
  numFakeSIRTiter = 0
  !
  !...... Default double - width linear interpolation in cosine stretching
  interpFacStretch = 2
  interpOrdStretch = 1
  interpOrdXtilt = 1
  perpendicular = .true.
  rotateBy90 = .false.
  reprojBP = .false.
  numReproj = 0
  adjustOrigin = .false.
  numViewSubtract = 0
  iactGpuFailOption = 0
  iactGpuFailEnviron = 0
  ifOutSirtProj = 0
  ifOutSirtRec = 0
  isignConstraint = 0
  vertOutFile = ' '
  vertBoundFile = ' '
  !
  !...... Default title
  call b3dDate(dat)
  call time(tim)
  !
  write(titlech, 49) 'Tomographic reconstruction', dat, tim
  if (recReproj) write(titlech, 49) 'Reprojection from tomogram', dat, tim
  read(titlech, '(20a4)') (title(kti), kti = 1, 20)
  !
  !

  if (PipGetString('TITLE', card) == 0) then
    write(titlech, 49) card(1:50), dat, tim
    read(titlech, '(20a4)') (title(kti), kti = 1, 20)
    write(6, 101) title
  endif
  !
  nfields = 0
  if (PipGetIntegerArray('SLICE', inum, nfields, LIMNUM) == 0) then
    if (nfields / 2 .ne. 1) &
        call exitError('Wrong number of fields on SLICE line')
    isliceStart = inum(1) + 1
    isliceEnd = inum(2) + 1
    if (nfields > 2) idelSlice = inum(3)
    if (idelSlice <= 0) call exitError( &
        'Negative slice increments are not allowed')
    ifSliceIn = 1
  endif
  !
  if (PipGetInteger('THICKNESS', ithickBP) == 0) ifThickIn = 1
  !
  if (PipGetInteger('MASK', numExtraMaskPix) == 0) then
    maskEdges = .true.
    numExtraMaskPix = max(-nprojXyz(1) / 50, min(nprojXyz(1) / 50, numExtraMaskPix))
    write(6, 401) numExtraMaskPix
  endif
  if3Dtexture = -999
  ierr = PipGetInteger('TextureType', if3Dtexture)
  !
  ! Hamming and RADIAL entry are mutually exclusive
  ifHammingLike = 1 - PipGetFloat('HammingLikeFilter', radMax)
  if (ifHammingLike > 0) then
    iradMax = nxProj * radMax
    radFall = 0.438 * (0.5 * nxProj - iradMax)
    ifMultByGaussian = 1
    write(6, 503) iradMax / float(nxProj), radFall / nxProj
  endif

  if (PipGetTwoFloats('RADIAL', radMax, radFall) == 0) then
    if (ifHammingLike > 0)call exitError('YOU CANNOT ENTER HammingLikeFilter WITH RADIAL')
    iradMax = radMax
    iradFall = radFall
    if (iradMax == 0) iradMax = nxProj * radMax
    if (iradFall == 0) radFall = nxProj * radFall
    !
    ! Adjust the falloff if it is NOT true sigma, the function now uses value correctly
    iv = 0
    ierr = PipGetBoolean('FalloffIsTrueSigma', iv)
    if (iv == 0) radFall = sqrt(0.5) * radFall
    write(6, 501) iradMax / float(nxProj), radFall / nxProj
  endif
  !
  if (PipGetBoolean('MultiplyByGaussian', ifMultByGaussian) == 0) then
    if (ifHammingLike > 0 .and. ifMultByGaussian == 0) call exitError( &
        'YOU CANNOT ENTER HammingLikeFilter AND MultiplyByGaussian 0')
    if (ifHammingLike == 0) write(6, 502)
  endif
  !
  nfields = 0
  if (PipGetFloatArray('OFFSET', xnum, nfields, LIMNUM) == 0) then
    if (nfields == 0 .or. nfields >= 3) &
        call exitError('Wrong number of fields on OFFSET line')
    if (nfields == 2) axisXoffset = xnum(2)
    delAngle = xnum(1)
  endif
  !
  if (PipGetTwoFloats('SCALE', outAdd, outScale) == 0) &
      write(6, 701) outAdd, outScale
  !
  iv = PipGetLogical('PERPENDICULAR', perpendicular)
  j = PipGetLogical('RotateBy90', rotateBy90)
  i = 0
  k = PipGetBoolean('PARALLEL', i)
  if ((iv == 0 .and. perpendicular .and. (i > 0 .or. rotateBy90)) .or.  &
      (i > 0 .and. rotateBy90)) &
      call exitError('YOU CAN SELECT ONLY ONE OF PERPENDICULAR, PARALLEL, AND RotateBy90')
  if (i > 0 .or. rotateBy90) perpendicular = .false.
  if (perpendicular) then
    write(6, 801)
  elseif (rotateBy90) then
    write(6, 902)
  else
    write(6, 901)
  endif
  !
  if (PipGetInteger('MODE', newMode) == 0) then
    if (newMode < 0 .or. newMode > 15 .or. (newMode > 2 .and. newMode < 9)) &
        call exitError('Illegal output mode')
    write(6, 1001) newMode
  endif
  !
  ierr = PipNumberOfEntries('INCLUDE', numEnt)
  do j = 1, numEnt
    ierr = PipGetString('INCLUDE', card)
    call parseList(card, mapUsedView(numViewUse + 1), numExcludeList)
    if (numViewUse + numExcludeList > limView) call exitError( &
        'TOO MANY INCLUDED VIEWS FOR ARRAYS')
    do i = numViewUse + 1, numViewUse + numExcludeList
      if (mapUsedView(i) < 1 .or. mapUsedView(i) > numViews) call exitError( &
          'Illegal view number in INCLUDE list')
    enddo
    numViewUse = numViewUse + numExcludeList
  enddo
  !
  ierr = PipNumberOfEntries('EXCLUDELIST2', numEnt)
  do j = 1, numEnt
    ierr = PipGetString('EXCLUDELIST2', card)
    call parseList(card, ivExclude(numViewExclude + 1), numExcludeList)
    if (numViewExclude + numExcludeList > limView) call exitError( &
        'TOO MANY EXCLUDED VIEWS FOR ARRAYS')
    do i = numViewExclude + 1, numViewExclude + numExcludeList
      if (ivExclude(i) < 1 .or. ivExclude(i) > numViews) call exitError( &
          'Illegal view number in EXCLUDE list')
    enddo
    numViewExclude = numViewExclude + numExcludeList
  enddo
  !
  if (numViewUse > 0 .and. numViewExclude > 0) call exitError( &
      'Illegal to have both INCLUDE and EXCLUDE entries')
  !
  if (PipGetFloat('LOG', baseForLog) == 0) then
    ifLog = 1
    write(6, 1301) baseForLog
  endif
  !
  ! Removed replications
  !
  ierr = PipNumberOfEntries('ANGLES', numEnt)
  do j = 1, numEnt
    nfields = 0
    ierr = PipGetFloatArray('ANGLES',  angles(newAngles + 1), nfields, &
        limView - newAngles)
    newAngles = newAngles + nfields
  enddo
  !
  ierr = PipNumberOfEntries('COMPRESS', numEnt)
  do j = 1, numEnt
    nfields = 0
    ierr = PipGetFloatArray('COMPRESS',  angles(numCompress + 1), nfields, &
        limView - numCompress)
    numCompress = numCompress + nfields
  enddo
  !
  if (PipGetFloat('COMPFRACTION', compFactor) == 0) &
      write(6, 1701) compFactor
  !
  nfields = 0
  if (PipGetFloatArray('DENSWEIGHT', xnum, nfields, LIMNUM) == 0) then
    numTiltIncWgt = nint(xnum(1))
    if (numTiltIncWgt > 0) then
      do i = 1, numTiltIncWgt
        tiltIncWgts(i) = 1. / (i - 0.5)
      enddo
      if (nfields == numTiltIncWgt + 1) then
        do i = 1, numTiltIncWgt
          tiltIncWgts(i) = xnum(i + 1)
        enddo
      elseif (nfields .ne. 1) then
        call exitError('Wrong number of fields on DENSWEIGHT line')
      endif
      write(6, 1801) numTiltIncWgt, (tiltIncWgts(i), i = 1, numTiltIncWgt)
    else
      write(6, 1802)
    endif
  endif
  !
  ! Keep this with the above entry so that nfields indicates if it was entered
  if (PipGetInteger('ExactFilterSize', numEnt) == 0) then
    if (nfields > 0) call exitError('YOU CANNOT ENTER DENSWEIGHT WITH ExactFilterSize')
    exactObjSize = numEnt
    allocate(exactTable(0 : nint(exactSamples)), stat = ierr)
    call memoryErrorUC(ierr, 'ARRAY FOR EXACT FILTER SAMPLES')
    numExactCycles = 1
    do i = 0, nint(exactSamples)
      exactTable(i) = 1. - i / exactSamples
    enddo
  endif
  !
  if (PipGetInteger('FakeSIRTiterations', numFakeSIRTiter) == 0 .and. numExactCycles > 0)&
      call exitError('YOU CANNOT ENTER FakeSIRTiterations WITH ExactFilterSize')
  !
  if (PipGetString('TILTFILE', card) == 0) then
    call dopen(3, card, 'ro', 'f')
    read(3,*,err = 2411, end = 2411) (angles(i), i = 1, numViews)
    close(3)
    ifTiltFile = 1
  endif
  !
  if (PipGetInteger('WIDTH', iwidth) == 0)  ifWidthIn = 1
  !
  nfields = 0
  if (PipGetFloatArray('SHIFT', xnum, nfields, LIMNUM) == 0) then
    if ((nfields + 1) / 2 .ne. 1) &
        call exitError('Wrong number of fields on SHIFT line')
    xOffset = xnum(1)
    if (nfields == 2) yOffset = xnum(2)
  endif
  !
  if (PipGetString('XTILTFILE', card) == 0) then
    call dopen(3, card, 'ro', 'f')
    read(3,*,err = 2412, end = 2412) (alpha(i), i = 1, numViews)
    close(3)
    ifAlpha = 1
    do i = 1, numViews
      if (abs(alpha(i) - alpha(1)) > 1.e-5) ifAlpha = 2
    enddo
    if (ifAlpha == 2) write(6, 2201)
    if (ifAlpha == 1) write(6, 2202) - alpha(1)
  endif
  !
  if (PipGetString('WeightFile', card) == 0) then
    ifExpWeight = 1
    call dopen(3, card, 'ro', 'f')
    read(3,*,err = 2414, end = 2414) (exposeWeight(i), i = 1, numViews)
    close(3)
  endif
  !
  boundFile = ' '
  ierr = PipGetString('BoundaryInfoFile', boundFile)
  !
  if (PipGetFloat('XAXISTILT', globalAlpha) == 0) then
    write(6, 2301) globalAlpha
    if (abs(globalAlpha) > 1.e-5 .and. ifAlpha == 0) ifAlpha = 1
  endif
  !
  ! REPROJECT entry must be read in before local alignments
  ! violates original unless a blank entry is allowed
  ierr = PipNumberOfEntries('REPROJECT', numEnt)
  do j = 1, numEnt
    nfields = 0
    ierr = PipGetFloatArray('REPROJECT',  xnum, nfields, LIMNUM)
    if (nfields == 0) then
      nfields = 1
      xnum(1) = 0.
    endif
    if (nfields + numReproj > limReproj) call exitError( &
        'TOO MANY REPROJECTION ANGLES FOR ARRAYS')
    do i = 1, nfields
      numReproj = numReproj + 1
      angReproj(numReproj) = xnum(i)
    enddo
    if (j == 1) write(6, 3101)
    reprojBP = .not.recReproj
  enddo
  !
  if (PipGetString('LOCALFILE', card) == 0) then
    call dopen(3, card, 'ro', 'f')
    read(3, '(a)', err = 2410, end = 2410) titlech
    ! read(3,*) nxWarp, nyWarp, ixStartWarp, iyStartWarp, idelXwarp, idelYwarp
    call frefor(titlech, xnum, numInput)
    ifDelAlpha = 0
    if (numInput > 6) ifDelAlpha = nint(xnum(7))
    pixelLocal = 0.
    if (numInput > 7) pixelLocal = xnum(8)
    localZfacs = 0
    if (numInput > 8) localZfacs = xnum(9)
    nxWarp = nint(xnum(1))
    nyWarp = nint(xnum(2))
    ixStartWarp = nint(xnum(3))
    iyStartWarp = nint(xnum(4))
    idelXwarp = nint(xnum(5))
    idelYwarp = nint(xnum(6))
    numWarpPos = nxWarp * nyWarp
    limWarp = nxWarp * nyWarp * numViews
    ipos = limWarp
    indDelta = numViews
    !
    ! If reprojecting rec, make sure arrays are big enough for all the
    ! reprojections and  allocate extra space at top of some arrays for
    ! temporary use
    if (recReproj) then
      indDelta = max(numViews, numReproj)
      limWarp = nxWarp * nyWarp * indDelta
      ipos = limWarp + indDelta
    endif
    if (nxWarp < 2 .or. nyWarp < 2) call exitError( &
        'THERE MUST BE AT LEAST TWO LOCAL ALIGNMENT AREAS IN X AND IN Y')
    allocate(indWarp(numWarpPos), delAlpha(ipos), cWarpBeta(limWarp), &
        sWarpBeta(limWarp), cWarpAlpha(limWarp), sWarpAlpha(limWarp), fwarp(2, 3, ipos), &
        delBeta(ipos), warpXZfac(ipos), warpYZfac(ipos), stat = ierr)
    call memoryErrorUC(ierr, 'ARRAYS FOR LOCAL ALIGNMENT DATA')
    indBase = 0
    do ipos = 1, nxWarp * nyWarp
      indWarp(ipos) = indBase
      read(3,*,err = 2410, end = 2410) (delBeta(i), i = indBase + 1, indBase + numViews)
      if (ifDelAlpha > 0) then
        read(3,*,err = 2410, end = 2410) (delAlpha(i), i = indBase + 1, indBase + &
            numViews)
      else
        do i = indBase + 1, indBase + numViews
          delAlpha(i) = 0.
        enddo
      endif
      !
      ! Set z factors to zero, read in if supplied, then negate them
      !
      do i = indBase + 1, indBase + numViews
        warpXZfac(i) = 0.
        warpYZfac(i) = 0.
      enddo
      if (localZfacs > 0) read(3,*,err = 2410, end = 2410) (warpXZfac(i), &
          warpYZfac(i), i = indBase + 1, indBase + numViews)
      do i = indBase + 1, indBase + numViews
        warpXZfac(i) = -warpXZfac(i)
        warpYZfac(i) = -warpYZfac(i)
      enddo
      do i = 1, numViews
        call xfread(3, fwarp(1, 1, i + indBase), 2410, 2410)
      enddo
      indBase = indBase + indDelta
    enddo
    close(3)
    write(6, 2401)
  endif
  !
  if (PipGetFloat('LOCALSCALE', scaleLocal) == 0) write(6, 2501) scaleLocal
  !
  ierr = PipGetTwoIntegers('FULLIMAGE', nxFullIn, nyFull)
  if (PipGetTwoIntegers('SUBSETSTART', ixSubsetIn, iySubset) == 0) &
      ifSubsetIn = 1
  !
  nfields = 0
  if (PipGetIntegerArray('COSINTERP', inum, nfields, LIMNUM) == 0) then
    interpOrdStretch = inum(1)
    if (nfields > 1) interpFacStretch = inum(2)
    interpOrdStretch = max(0, min(3, interpOrdStretch))
    if (interpOrdStretch == 0) interpFacStretch = 0
    if (interpFacStretch == 0) then
      print *,'Cosine stretching is disabled'
    else
      write(6, 2801) interpOrdStretch, interpFacStretch
    endif
  endif
  !
  if (PipGetInteger('XTILTINTERP', interpOrdXtilt) == 0) then
    if (interpOrdXtilt > 2) interpOrdXtilt = 3
    if (interpOrdXtilt <= 0) then
      print *,'New-style X-tilting with vertical slices is disabled'
    else
      write(6, 3001) interpOrdXtilt
    endif
  endif
  !
  ! Read environment variable first, then override by entry
  if (imodGetEnv('IMOD_USE_GPU', card) == 0) read(card,*) indGPU
  ifGpuByEnviron = PipGetInteger('UseGPU', indGPU)
  if (imodGetEnv('IMOD_USE_GPU2', card) == 0) then
    read(card,*) indGPU
    ifGpuByEnviron = 1
  endif
  useGPU = indGPU >= 0 .and. threshPolarity == 0
  ierr = PipGetTwoIntegers('ActionIfGPUFails', iactGpuFailOption, &
      iactGpuFailEnviron)
  !
  if (PipGetString('ZFACTORFILE', card) == 0) then
    call dopen(3, card, 'ro', 'f')
    read(3,*,err = 2413, end = 2413) (xzfac(i), yzfac(i), i = 1, numViews)
    close(3)
    ifZfactors = 1
    write(6, 3201)
  endif
  !
  if (PipGetInteger('IMAGEBINNED', imageBinned) == 0) then
    imageBinned = max(1, imageBinned)
    if (imageBinned > 1) write(6, 3301) imageBinned
  endif
  !
  if (PipGetTwoIntegers('TOTALSLICES', inum(1), inum(2)) == 0) then
    minTotSlice = inum(1) + 1
    maxTotSlice = inum(2) + 1
  endif
  !
  ierr = PipGetFloat('FlatFilterFraction', flatFrac)
  flatFrac = max(0.,  flatFrac)
  if (flatFrac > 0 .and. numFakeSIRTiter > 0) call exitError( &
      'YOU CANNOT ENTER FlatFilterFraction WITH FakeSIRTiterations')
  !
  ierr = PipGetLogical('AdjustOrigin', adjustOrigin)
  !
  if (.not. recReproj) &
      ierr = PipGetThreeFloats('MinMaxMean', dminIn, dmaxIn, dmeanIn)
  !
  if (PipGetString('WeightAngleFile', card) == 0) then
    call dopen(3, card, 'ro', 'f')
313 read(3,*,err = 2415, end = 314) wgtAngles(numWgtAngles + 1)
    numWgtAngles = numWgtAngles + 1
    go to 313
    !
    ! Sort the angles
314 do i = 1, numWgtAngles - 1
      do j = i + 1, numWgtAngles
        if (wgtAngles(i) > wgtAngles(j)) then
          dminTmp = wgtAngles(i)
          wgtAngles(i) = wgtAngles(j)
          wgtAngles(j) = dminTmp
        endif
      enddo
    enddo
  endif
  !
  ! SIRT-related options
  ierr = PipGetInteger('ConstrainSign', isignConstraint)
  !
  ierr = PipGetTwoIntegers('InternalSIRTSlices', ifOutSirtProj, ifOutSirtRec)
  if (numSIRTiter > 0 .or. recSubtraction) &
      ierr = PipGetInteger('StartingIteration', iterForReport)
  ierr = PipGetLogical('VertForSIRTInput', vertSirtInput)
  ierr = PipGetString('VertSliceOutputFile', vertOutFile)
  ierr = PipGetString('VertBoundaryFile', vertBoundFile)
  !
  call PipDone()
  !
  ! END OF OPTION READING
  !
  write(6, 48)
  if (ifAlpha .ne. 0 .and. idelSlice .ne. 1) call exitError( &
      'Cannot do X axis tilt with non-consecutive slices')
  if (nxWarp .ne. 0 .and. idelSlice .ne. 1) call exitError( &
      'Cannot do local alignments with non-consecutive slices')
  if (minTotSlice > 0 .and. idelSlice .ne. 1) call exitError( &
      'Cannot do chunk writing with non-consecutive slices')
  if (nxFullIn == 0 .and. nyFull == 0 .and. &
      (ixSubsetIn .ne. 0 .or. iySubset .ne. 0)) call exitError( &
      'YOU MUST ENTER THE FULL IMAGE SIZE IF YOU HAVE A SUBSET')
  if (.not.perpendicular .and. minTotSlice > 0) call exitError( &
      'Cannot do chunk writing with parallel slices')
  if (numReproj > 0 .and. nViewsReproj > 0) call exitError( &
      'You cannot enter both views and angles to reproject')
  if (projModel .and. (numReproj > 0 .or. &
      nViewsReproj > 0)) call exitError('You cannot do projection '// &
      'from a model with image reprojection')
  if (numSIRTiter > 0 .and. (numReproj > 0 .or. &
      nViewsReproj > 0)) call exitError('You cannot do SIRT '// &
      'with entries for angles/views to reproject')
  !
  ! Scale dimensions down by binning then report them
  !
  ixSubset = ixSubsetIn
  nxFull = nxFullIn
  if (imageBinned > 1) then
    if (ifSliceIn .ne. 0 .and. (minTotSlice <= 0 .or. isliceStart > 0)) then
      isliceStart = max(1, min(nprojXyz(2), &
          (isliceStart + imageBinned - 1) / imageBinned))
      isliceEnd = max(1, min(nprojXyz(2), &
          (isliceEnd + imageBinned - 1) / imageBinned))
    endif
    if (ifThickIn .ne. 0) ithickBP = ithickBP / imageBinned
    axisXoffset = axisXoffset / imageBinned
    nxFull = (nxFullIn + imageBinned - 1) / imageBinned
    nyFull = (nyFull + imageBinned - 1) / imageBinned
    ixSubset = ixSubsetIn / imageBinned
    iySubset = iySubset / imageBinned
    xOffset = xOffset / imageBinned
    yOffset = yOffset / imageBinned
    if (ifWidthIn .ne. 0) iwidth = iwidth / imageBinned
    if (minTotSlice >= 0 .and. .not. recReproj) then
      minTotSlice = max(1, min(nprojXyz(2), &
          (minTotSlice + imageBinned - 1) / imageBinned))
      maxTotSlice = max(1, min(nprojXyz(2), &
          (maxTotSlice + imageBinned - 1) / imageBinned))
    endif
    if (numExactCycles > 0) exactObjSize = exactObjSize / imageBinned
  endif
  !
  ! Set up an effective scaling factor for non-log scaling to test output
  ! The equation in ( ) is based on fitting to the amplification of the range
  ! with linear scaling for different input sizes
  effectiveScale = 1.
  if (ifLog == 0) effectiveScale = outScale * (0.011 * nprojXyz(1) + 6) / 2.
  !
  if (recReproj) then
    if (debug) print *,minTotSlice, maxTotSlice, minZreproj, maxZreproj, nrecXyz(3)
    if ((minTotSlice <= 0 .and. (minZreproj <= 0 .or. maxZreproj > &
        nrecXyz(3))) .or. (minTotSlice >= 0 .and. &
        maxTotSlice > nrecXyz(3))) call exitError( &
        'Min or Max Z coordinate to project is out of range')

    if (.not.perpendicular) call exitError( &
        'Cannot reproject from reconstruction output with PARALLEL')
    if (idelSlice .ne. 1) call exitError( &
        'Cannot reproject from reconstruction with a slice increment')
    if (iwidth .ne. nrecXyz(1) .or. isliceEnd + 1 - isliceStart .ne. nrecXyz(3) .or. &
        ithickBP .ne. nrecXyz(2)) call exitError( &
        'Dimensions of rec file do not match expected values')
  else
    !
    ! Check conditions of SIRT (would slice increment work?)
    if (numSIRTiter > 0) then
      if (.not.perpendicular .or. idelSlice .ne. 1 .or. &
          xOffset .ne. 0.) call exitError('Cannot do SIRT with PARALLEL'// &
          ' output, slice increment, or X shifts')
      if (nxWarp .ne. 0 .or. ifZfactors .ne. 0 .or. ifAlpha > 1 .or. &
          (ifAlpha == 1 .and. interpOrdXtilt == 0)) &
          call exitError('Cannot do SIRT with  local alignments, Z '// &
          'factors, or variable or old-style X tilt')
      if (iwidth .ne. nxProj .or. (.not. sirtFromZero .and. (iwidth .ne. nrecXyz(1) .or. &
          (ithickBP .ne. nrecXyz(2) .and. .not.vertSirtInput) .or.  &
          nrecXyz(3) .ne. nyProj))) call exitError( 'For SIRT, sizes of '// &
          'input projections, rec file, and width/thickness entries must match')
      write(6, 3501) numSIRTiter
      interpFacStretch = 0
    endif
    if (ifSliceIn .ne. 0) write(6, 201) isliceStart, isliceEnd, idelSlice
    if (ifThickIn .ne. 0) write(6, 301) ithickBP
    if (delAngle .ne. 0. .or. axisXoffset .ne. 0.) write(6, 601) delAngle, axisXoffset
    if (nxFull .ne. 0 .or. nyFull .ne. 0) write(6, 2601) nxFull, nyFull
    if (ifSubsetIn .ne. 0) write(6, 2701) ixSubset, iySubset
    if (ifWidthIn .ne. 0) write(6, 2001) iwidth
    if (xOffset .ne. 0 .or. yOffset .ne. 0) write(6, 2101) yOffset, xOffset
    if (minTotSlice > 0) write(6, 3401) minTotSlice, maxTotSlice
  endif
  !
  ! If NEWANGLES is 0, get angles from file header.  Otherwise check if angles OK
  !
  if (newAngles == 0 .and. ifTiltFile == 0) then
    !
    ! Tilt information is stored in stack header. Read into angles
    ! array. All sections are assumed to be equally spaced. If not,
    ! you need to set things up differently. In such a case, great
    ! care should be taken, since missing views may have severe
    ! effects on the quality of the reconstruction.
    !
    !
    ! call irtdat(1, idtype, lens, nd1, nd2, vd1, vd2)
    !
    if (idType .ne. 1) call exitError( 'There are no tilt angles from ANGLES or '// &
        'TILTFILE entries or from the image file header')
    !
    if (nd1 .ne. 2) call exitError(' Tilt axis not along Y.')
    !
    delTheta = vd1
    theta = vd2
    !
    do nv = 1, numViews
      angles(nv) = theta
      theta = theta + delTheta
    enddo
    !
  else
    if (ifTiltFile == 1 .and. newAngles .ne. 0) then
      call exitError('Tried to enter angles with both ANGLES and TILTFILE')
    elseif (ifTiltFile == 1) then
      write(6,*) ' Tilt angles were entered from a tilt file'
    elseif (newAngles == numViews) then
      write(6,*) ' Tilt angles were entered with ANGLES card(s)'
    else
      call exitError('If using ANGLES, a value must be entered for each view')
    endif
  endif
  !
  if (numCompress > 0) then
    if (numCompress == numViews) then
      write(6,*) ' Compression values were entered with COMPRESS card(s)'
    else
      call exitError('If using COMPRESS, a value must be entered for each view')
    endif
    do nv = 1, numViews
      compress(nv) = 1. +(compress(nv) - 1.) / compFactor
    enddo
  endif
  !
  if (globalAlpha .ne. 0.) then
    do iv = 1, numViews
      alpha(iv) = alpha(iv) - globalAlpha
    enddo
  endif
  !
  if (ifExpWeight .ne. 0) then
    if (ifLog == 0) write(6,*) ' Weighting factors were entered from a file'
    if (ifLog .ne. 0) write(6,*) ' Weighting factors were entered '// &
        'but will be ignored because log is being taken'
  endif
  !
  ! if no INCLUDE cards, set up map to views, excluding any specified by
  ! EXCLUDE cards
  if (numViewUse == 0) then
    do i = 1, numViews
      ierr = 0
      do iex = 1, numViewExclude
        if (i == ivExclude(iex)) ierr = 1
      enddo
      if (ierr == 0) then
        numViewUse = numViewUse + 1
        mapUsedView(numViewUse) = i
      endif
    enddo
  endif
  !
  ! Replace angles at +/-90 with 89.95 etc
  do i = 1, numViewUse
    j = mapUsedView(i)
    if (abs(abs(angles(j)) - 90.) < 0.05) &
        angles(j) = sign(90. - sign(0.05, 90 - abs(angles(j))), angles(j))
  enddo
  !
  ! If reprojecting from rec and no angles entered, copy angles in original
  ! order
  if (recReproj .and. numReproj == 0) then
    if (nViewsReproj == 1 .and. ivReproj(1) == 0) then
      numReproj = numViews
      do i = 1, numViews
        angReproj(i) = angles(i)
      enddo
    else if (nViewsReproj > 0) then
      numReproj = nViewsReproj
      do i = 1, numReproj
        if (ivReproj(i) < 1 .or. ivReproj(i) > numViews) call exitError( &
            'View number to reproject is out of range')
        angReproj(i) = angles(ivReproj(i))
      enddo
    else
      !
      ! For default set of included views, order them by view number by
      ! first ordering the mapUsedView array by view number.
      do i = 1, numViewUse-1
        do j = i + 1, numViewUse
          if (mapUsedView(i) > mapUsedView(j)) then
            indi = mapUsedView(i)
            mapUsedView(i) = mapUsedView(j)
            mapUsedView(j) = indi
          endif
        enddo
      enddo
      numReproj = numViewUse
      do i = 1, numViewUse
        angReproj(i) = angles(mapUsedView(i))
      enddo
    endif
  endif
  !
  ! order the MAPUSE array by angle
  do i = 1, numViewUse-1
    do j = i + 1, numViewUse
      indi = mapUsedView(i)
      if (angles(indi) > angles(mapUsedView(j))) then
        mapUsedView(i) = mapUsedView(j)
        mapUsedView(j) = indi
        indi = mapUsedView(i)
      endif
    enddo
  enddo
  !
  ! For SIRT, now copy the angles in the ordered list; this is the order
  ! in which the reprojections are needed internally
  ! Also adjust the recon mean for a fill value
  if (numSIRTiter > 0) then
    numReproj = numViewUse
    do i = 1, numViewUse
      angReproj(i) = angles(mapUsedView(i))
    enddo
    outAdd = outAdd / filterScale
    outScale = outScale * filterScale
    dmeanIn = dmeanIn / outScale - outAdd
  endif
  if (numReproj > 0) then
    allocate(cosReproj(numReproj), sinReproj(numReproj), stat = ierr)
    call memoryErrorUC(ierr, 'ARRAYS FOR REPROJECTION SINES/COSINES')
    do i = 1, numReproj
      cosReproj(i) = cos(degToRad * angReproj(i))
      sinReproj(i) = sin(degToRad * angReproj(i))
    enddo
  endif

  !
  ! Open output map file
  call iiuRetDelta(1, delta)
  call iiuRetOrigin(1, originX, originY, originZ)
  if (.not. recReproj) then
    if ((minTotSlice <= 0 .and. (isliceStart < 1 .or. isliceEnd < 1)) &
        .or. isliceStart > nprojXyz(2) .or. isliceEnd > nprojXyz(2)) call exitError( &
        'SLICE NUMBERS OUT OF RANGE')
    numSlices = (isliceEnd - isliceStart) / idelSlice + 1
    if (minTotSlice > 0 .and. isliceStart < 1) &
        numSlices = maxTotSlice + 1 - minTotSlice
    ! print *,'NSLICE', minTotSlice, maxTotSlice, isliceStart, numSlices
    if (numSlices <= 0) call exitError( 'SLICE NUMBERS REVERSED')
    if (reprojBP .or. readBaseRec) then
      allocate(projLine(iwidth), stat = ierr)
      call memoryErrorUC(ierr, 'ARRAY FOR PROJECTION LINE')
    endif
    if (defocusFile .ne. ' ' .and. pixForDefocus == 0.) pixForDefocus = delta(1) / 10.
    !
    ! DNM 7/27/02: transfer pixel sizes depending on orientation of output
    !
    noutXyz(1) = iwidth
    cell(1) = iwidth * delta(1)
    if (perpendicular) then
      noutXyz(2) = ithickBP
      noutXyz(3) = numSlices
      cell(2) = ithickBP * delta(1)
      cell(3) = numSlices * idelSlice * delta(2)
    else
      noutXyz(2) = numSlices
      noutXyz(3) = ithickBP
      cell(3) = ithickBP * delta(1)
      cell(2) = numSlices * idelSlice * delta(2)
    END if
    if (reprojBP) then
      noutXyz(2) = numSlices
      noutXyz(3) = numReproj
      cell(2) = numSlices * idelSlice * delta(2)
      cell(3) = delta(1) * numReproj
      j = iwidth * numReproj
      allocate(xRayStart(j), yRayStart(j), numPixInRay(j), maxRayPixels(numReproj),  &
          stat = ierr)
      call memoryErrorUC(ierr, 'ARRAYS FOR PROJECTION RAY DATA')
      do i = 1, numReproj
        j = (i - 1) * iwidth + 1
        !
        ! Note that this will set cosReproj to 0 after you carefully kept it from being 0
        call set_projection_rays(sinReproj(i), cosReproj(i), iwidth, ithickBP, &
            iwidth, xRayStart(j), yRayStart(j), numPixInRay(j), maxRayPixels(i))
      enddo
    endif
  else
    !
    ! recReproj stuff
    noutXyz(1) = maxXreproj + 1 - minXreproj
    noutXyz(2) = maxZreproj + 1 - minZreproj
    ithickReproj = maxYreproj + 1 - minYreproj
    noutXyz(3) = numReproj
    if (minTotSlice > 0) noutXyz(2) = maxTotSlice + 1 - minTotSlice
    if (noutXyz(1) < 1 .or. noutXyz(2) < 1 .or. ithickReproj < 1) &
        call exitError('Min and max limits for output are reversed for X, Y, or Z')
    cell(1) = noutXyz(1) * delta(1)
    cell(2) = noutXyz(2) * delta(1)
    cell(3) = delta(1) * numReproj
    if (projSubtraction .and. (noutXyz(1) .ne. nxProj .or. maxZreproj > &
        nyProj .or. noutXyz(3) .ne. nprojXyz(3))) call exitError('OUTPUT SIZE'// &
        ' MUST MATCH ORIGINAL PROJECTION FILE SIZE FOR SIRT SUBTRACTION')
  endif
  !
  ! Check compatibility of base rec file
  if ((readBaseRec .and. .not. recReproj) .and. (nrecXyz(1) .ne. iwidth .or. &
      nrecXyz(2) .ne. ithickBP .or. nrecXyz(3) < isliceEnd)) call exitError( &
      'BASE REC FILE IS NOT THE SAME SIZE AS OUTPUT FILE')
  !
  ! Initialize parallel writing routines if bound file entered
  parallelHDF = .false.
  if (.not.projModel) then
    if (minTotSlice > 0 .and. ((.not.recReproj .and. isliceStart > 0) .or. &
        (recReproj .and. minZreproj > 0))) then
      parallelHDF = boundFile .ne. ' ' .and. iiTestIfHDF(outputFile) > 0
    else if (minTotSlice > 0) then
      parallelHDF = b3dOutputFileType() == 5
    endif
  endif
  ind1 = noutXyz(1)
  if (parallelHDF) ind1 = -ind1
  ierr = parWrtInitialize(boundFile, 7, ind1, noutXyz(2), noutXyz(3))
  if (ierr .ne. 0) then
    write(*,'(a,i3)') 'ERROR: TILT - INITIALIZING PARALLEL WRITE BOUNDARY FILE, ERROR', &
        ierr
    call exit(1)
  endif
  !
  ! open old file if in chunk mode and there is real starting slice
  ! otherwise open new file
  !
  if (.not.projModel) then
    if (minTotSlice > 0 .and. ((.not.recReproj .and. isliceStart > 0) .or. &
        (recReproj .and. minZreproj > 0))) then
      if (parallelHDF) then
        call parWrtProperties(ind1, ind2, k)
        if (b3dLockFile(ind1) .ne. 0) &
            call exitError("COULD NOT GET LOCK FOR OPENING HDF FILE")
        ind2 = noutXyz(3)
      endif
      call imOpen(2, outputFile, 'OLD')
      call irdhdr(2, noutXyz, mpxyz, newMode, dminTmp, dmaxTmp, dmeanTmp)
      if (parallelHDF) noutXyz(3) = ind2
    else
      call imOpen(2, outputFile, 'NEW')
      call iiuCreateHeader(2, noutXyz, noutXyz, newMode, title, 0)
      !print *,'created', noutXyz
    endif
    if (parallelHDF .and. iiuFileType(2) .ne. 5)  &
        call exitError('EXPECTED OUTPUT FILE TO BE AN HDF FILE BUT IT IS NOT')
    call iiuTransLabels(2, 1)
    call iiuAltCell(2, cell)
  endif
  !
  ! if doing perpendicular slices, set up header info to make coordinates
  ! congruent with those of tilt series
  !
  if (recReproj) then
    call iiuAltOrigin(2, 0., 0., 0.)
    outHdrTilt(1) = 0.
    call iiuAltTilt(2, outHdrTilt)
  else if (perpendicular) then
    outHdrTilt(1) = 90
    if (adjustOrigin) then
      !
      ! Full adjustment if requested
      originX = originX  - delta(1) * (nprojXyz(1) / 2 - iwidth / 2 - xOffset)
      originZ = originY - delta(1) * float(max(0, isliceStart - 1))
      if (minTotSlice > 0 .and. isliceStart <= 0) &
          originZ = originY - delta(1) * (minTotSlice-1)
      originY = delta(1) * (ithickBP / 2 + yOffset)
    else
      !
      ! Legacy origin.  All kinds of wrong.
      originX = cell(1) / 2. +axisXoffset
      originY = cell(2) / 2.
      originZ = -float(max(0, isliceStart - 1))
    endif

    if (.not.projModel) then
      call iiuAltOrigin(2, originX, originY, originZ)
      call iiuAltTilt(2, outHdrTilt)
      call iiuAltSpaceGroup(2, 1)
    endif
  endif
  !
  ! chunk mode starter run: write header and exit
  !
  ifExit = 0
  if (.not.projModel) then
    if (minTotSlice > 0 .and. ((.not.recReproj .and. isliceStart <= 0) .or. &
        (recReproj .and. minZreproj <= 0))) then
      if (parallelHDF) call iiuWriteDummySecToHDF(2)
      call iiuWriteHeader(2, title, 1, dminIn, dmaxIn, dmeanIn)
      call iiuClose(2)
      ifExit = 1
    elseif (minTotSlice > 0) then
      if (.not. recReproj) call parWrtPosn(2, isliceStart - minTotSlice, 0)
      if (readBaseRec .and. .not. recReproj) &
          call iiuSetPosition(3, isliceStart - minTotSlice, 0)
      if (parallelHDF) call iiuParWrtRecloseHDF(2, 1)
    endif
  endif
  !
  ! If reprojecting, need to look up each angle in full list of angles and
  ! find ones to interpolate from, then pack data into arrays that are
  ! otherwise used for packing these factors down
  if (recReproj) then
    do i = 1, numReproj
      sinBeta(i) = angReproj(i)
      call lookupAngle(angReproj(i), angles, numViews, ind1, ind2, frac)
      cosBeta(i) = (1. - frac) * compress(ind1) + frac * compress(ind2)
      sinAlpha(i) = (1. - frac) * alpha(ind1) + frac * alpha(ind2)
      cosAlpha(i) = (1. - frac) * xzfac(ind1) + frac * xzfac(ind2)
      array(i) = (1. - frac) * yzfac(ind1) + frac * yzfac(ind2)
      array(i + numViews) = (1. - frac) * exposeWeight(ind1) + frac * exposeWeight(ind2)
    enddo
    !
    ! Do the same thing with all the local data: pack it into the spot at
    ! the top of the local data then copy it back into the local area
    if (nxWarp > 0) then
      do i = 1, nxWarp * nyWarp
        do iv = 1, numReproj
          call lookupAngle(angReproj(iv), angles, numViews, ind1, ind2, frac)
          ind1 = ind1 + indWarp(i)
          ind2 = ind2 + indWarp(i)
          delBeta(indBase + iv) = (1. -frac) * delBeta(ind1) + frac * delBeta(ind2)
          delAlpha(indBase + iv) = (1. -frac) * delAlpha(ind1) + frac * delAlpha(ind2)
          warpXZfac(indBase + iv) = (1. -frac) * warpXZfac(ind1) + frac * warpXZfac(ind2)
          warpYZfac(indBase + iv) = (1. -frac) * warpYZfac(ind1) + frac * warpYZfac(ind2)
          do j = 1, 2
            do k = 1, 3
              fwarp(j, k, indBase + iv) = (1. -frac) * fwarp(j, k, ind1) + &
                  frac * fwarp(j, k, ind2)
            enddo
          enddo
        enddo
        do iv = 1, numReproj
          ind1 = indBase + iv
          ind2 = indWarp(i) + iv
          delBeta(ind2) = delBeta(ind1)
          delAlpha(ind2) = delAlpha(ind1)
          warpXZfac(ind2) = warpXZfac(ind1)
          warpYZfac(ind2) = warpYZfac(ind1)
          do j = 1, 2
            do k = 1, 3
              fwarp(j, k, ind2) = fwarp(j, k, ind1)
            enddo
          enddo
        enddo
      enddo
    endif
    !
    ! Replace the mapUsedView array
    numViewUse = numReproj
    do i = 1, numViewUse
      mapUsedView(i) = i
    enddo
  else
    !
    ! pack angles and other data down as specified by MAPUSE
    ! Negate the z factors since things are upside down here
    ! Note that local data is not packed but always referenced by mapUsedView
    !
    do i = 1, numViewUse
      sinBeta(i) = angles(mapUsedView(i))
      cosBeta(i) = compress(mapUsedView(i))
      sinAlpha(i) = alpha(mapUsedView(i))
      cosAlpha(i) = xzfac(mapUsedView(i))
      array(i) = yzfac(mapUsedView(i))
      array(i + numViews) = exposeWeight(mapUsedView(i))
    enddo
  endif
  do i = 1, numViewUse
    angles(i) = sinBeta(i)
    compress(i) = cosBeta(i)
    alpha(i) = sinAlpha(i)
    xzfac(i) = -cosAlpha(i)
    yzfac(i) = -array(i)
    exposeWeight(i) = array(i + numViews)
  enddo
  numViewOrig = numViews
  numViews = numViewUse
  !
  write(6, 51) (angles(nv), nv = 1, numViews)
  write(6, 52)
  !
  ! Turn off cosine stretch for high angles
  if (angles(1) < -80. .or. angles(numViews) > 80.) then
    if (interpFacStretch > 0) write(*,662)
662 format(/,'Tilt angles are too high to use cosine stretching')
    interpFacStretch = 0
  endif
  !
  ! Set up trig tables -  Then convert angles to radians
  !
  do  iv = 1, numViews
    thetaView = angles(iv) + delAngle
    if (thetaView > 180.) thetaView = thetaView - 360.
    if (thetaView <= -180.) thetaView = thetaView + 360.
    cosBeta(iv) = cos(thetaView * degToRad)
    !
    ! Keep cosine from going to zero so it can be divided by
    if (abs(cosBeta(iv)) < 1.e-6) cosBeta(iv) = sign(1.e-6, cosBeta(iv))
    !
    ! Take the negative of the sine of tilt angle to account for the fact
    ! that all equations are written for rotations in the X/Z plane,
    ! viewed from  the negative Y axis
    ! Take the negative of alpha because the entered value is the amount
    ! that the specimen is tilted and we need to rotate by negative of that
    sinBeta(iv) = -sin(thetaView * degToRad)
    cosAlpha(iv) = cos(alpha(iv) * degToRad)
    sinAlpha(iv) = -sin(alpha(iv) * degToRad)
    angles(iv) = -degToRad * (angles(iv) + delAngle)
  enddo
  !
  ! If there are weighting angles, convert those the same way, otherwise
  ! copy the main angles to weighting angles
  if (numWgtAngles > 0) then
    do iv = 1, numWgtAngles
      wgtAngles(iv) = -degToRad * (wgtAngles(iv) + delAngle)
    enddo
  else
    numWgtAngles = numViews
    do iv = 1, numViews
      wgtAngles(iv) = angles(iv)
    enddo
  endif
  !
  ! if fixed x axis tilt, set up to try to compute vertical planes
  ! and interpolate output planes: adjust thickness that needs to
  ! be computed, and find number of vertical planes that are needed
  !
  if (ifZfactors > 0 .and. ifAlpha == 0) ifAlpha = 1
  ithickOut = ithickBP
  ycenModProj = ithickBP / 2 + 0.5 + yOffset
  if (ifAlpha == 1 .and. nxWarp == 0 .and. interpOrdXtilt > 0 .and. &
      ifZfactors == 0 .and. .not.recReproj) then
    ifAlpha = -1
    ithickOut = ithickBP
    ithickBP = ithickBP / cosAlpha(1) + 4.5
    numVertNeeded = ithickOut * abs(sinAlpha(1)) + 5.
    numReadNeed = 0
    if (numSIRTiter > 0) then
      saveVertSlices = vertOutFile .ne. ' '
      if (.not.sirtFromZero .and. .not.vertSirtInput) &
          numReadNeed = ithickBP * abs(sinAlpha(1)) + 4.
      if (.not.sirtFromZero .and. vertSirtInput .and. ithickBP .ne. nrecXyz(2)) &
          call exitError( &
          'THICKNESS OF VERTICAL SLICE INPUT FILE DOES NOT MATCH NEEDED THICKNESS')
    endif
  endif
  if (ifAlpha .ne. -1 .and. numSIRTiter > 0 .and. (vertSirtInput .or. vertOutFile &
      .ne. ' ')) call exitError('VertForSIRTInput OR VertSliceOutputFile '// &
      'CANNOT BE ENTERED UNLESS VERTICAL SLICES ARE BEING COMPUTED')
  !
  ! Now that we know vertical slices, open or set up output file for them under SIRT
  if (numSIRTiter > 0 .and. saveVertSlices) then
    !
    ! Initialize parallel writing if vertical bound file
    nvsXyz = noutXyz
    nvsXyz(2) = ithickBP
    if (parallelHDF .and. minTotSlice > 0 .and. isliceStart > 0 .and.  &
        vertBoundFile == ' ')  &
        call exitError('A BOUNDARY FILE FOR VERTICAL SLICES MUST BE ENTERED IF THERE '// &
        'IS ONE FOR PRIMARY OUTPUT AND OUTPUT FILE TYPE IS HDF')
    ind1 = nvsXyz(1)
    if (parallelHDF) ind1 = -ind1
    ierr = parWrtInitialize(vertBoundFile, 8, ind1, nvsXyz(2), nvsXyz(3))
    if (ierr .ne. 0) then
      write(*,'(a,i3)') 'ERROR: TILT - INITIALIZING PARALLEL WRITE '// &
          'BOUNDARY FILE FOR VERTICAL SLICES, ERROR', ierr
      call exit(1)
    endif
    if (minTotSlice > 0 .and. isliceStart > 0) then
      if (parallelHDF) then
        call parWrtProperties(ind1, ind2, k)
        if (b3dLockFile(ind1) .ne. 0)  &
            call exitError('COULD NOT GET LOCK FOR OPENING HDF FILE')
        ind2 = nvsXyz(3)
      endif
      call imOpen(6, vertOutFile, 'OLD')
      call irdhdr(6, nvsXyz, mpxyz, newMode, dminTmp, dmaxTmp, dmeanTmp)
      if (parallelHDF) nvsXyz(3) = ind2
    else
      call imOpen(6, vertOutFile, 'NEW')
      call iiuCreateHeader(6, nvsXyz, nvsXyz, 2, title, 0)
    endif
    if (parallelHDF .and. iiuFileType(6) .ne. 5)  &
        call exitError('EXPECTED VERTICAL SLICE FILE TO BE AN HDF FILE BUT IT IS NOT')
    !
    ! chunk mode: either write header, or set up to write correct location
    if (ifExit .ne. 0) then
      if (parallelHDF) call iiuWriteDummySecToHDF(6)
      call iiuWriteHeader(6, title, 0, dminIn, dmaxIn, dmeanIn)
      call iiuClose(6)
    elseif (minTotSlice > 0) then
      call parWrtPosn(6, isliceStart - minTotSlice, 0)
      if (parallelHDF) call iiuParWrtRecloseHDF(6, 1)
    endif
    ierr = parWrtSetCurrent(1)
  endif


  if (ifExit .ne. 0) then
    print *,'Exiting after setting up output file for chunk writing'
    call exit(0)
  endif
  !
  ! Set center of output plane and center of input for transformations
  ! Allow the full size to be less than the aligned stack, with a negative
  ! subset start
  !
  if (nxFull == 0) then
    nxFull = nprojXyz(1)
    nxFullIn = nprojXyz(1) * imageBinned
  endif
  if (nyFull == 0) nyFull = nprojXyz(2)
  xcenIn = nxFull / 2. + 0.5 - ixSubset
  centerSlice = nyFull / 2. + 0.5 - iySubset
  xoffAdj = xoffset - ((nprojXyz(1) * imageBinned - nxFullIn) / 2 + ixSubsetIn)
  xcenOut = iwidth / 2 + 0.5 + axisXoffset + xoffAdj
  ycenOut = ithickBP / 2 + 0.5 + yOffset
  if (numSIRTiter > 0 .and. abs(xoffAdj + axisXoffset) > 0.1) then
    if (mod(nprojXyz(1), 2) > 0) call exitError( &
        'CANNOT DO INTERNAL SIRT WITH AN ODD INPUT SIZE IN X')
    call exitError( &
        'CANNOT DO INTERNAL SIRT WITH A TILT AXIS OFFSET FROM CENTER OF INPUT IMAGES')
  endif
  !
  ! if doing warping, convert the angles to radians and set sign
  ! Also cancel the z factors if global entry was not made
  !
  if (nxWarp > 0) then
    do i = 1, numViewOrig * nxWarp * nyWarp
      delBeta(i) = -degToRad * delBeta(i)
    enddo
    do iv = 1, numViews
      do i = 1, nxWarp * nyWarp
        ind = indWarp(i) + mapUsedView(iv)
        cWarpBeta(ind) = cos(angles(iv) + delBeta(ind))
        sWarpBeta(ind) = sin(angles(iv) + delBeta(ind))
        cWarpAlpha(ind) = cos(degToRad * (alpha(iv) + delAlpha(ind)))
        sWarpAlpha(ind) = -sin(degToRad * (alpha(iv) + delAlpha(ind)))
        if (ifZfactors == 0) then
          warpXZfac(ind) = 0.
          warpYZfac(ind) = 0.
        endif
      enddo
    enddo
    !
    ! See if local outScale was entered; if not see if it can be set from
    ! pixel size and local align pixel size
    if (scaleLocal <= 0.) then
      scaleLocal = 1.
      if (pixelLocal > 0) then
        scaleLocal = pixelLocal / delta(1)
        if (abs(scaleLocal - 1.) > 0.001) write(6, 53) scaleLocal
      endif
    endif
    !
    ! outScale the x and y dimensions and shifts if aligned data were
    ! shrunk relative to the local alignment solution
    ! 10/16/04: fixed to use mapUsedView to outScale used views properly
    !
    if (scaleLocal .ne. 1.) then
      ixStartWarp = nint(ixStartWarp * scaleLocal)
      iyStartWarp = nint(iyStartWarp * scaleLocal)
      idelXwarp = nint(idelXwarp * scaleLocal)
      idelYwarp = nint(idelYwarp * scaleLocal)
      do iv = 1, numViews
        do i = 1, nxWarp * nyWarp
          ind = indWarp(i) + mapUsedView(iv)
          fwarp(1, 3, ind) = fwarp(1, 3, ind) * scaleLocal
          fwarp(2, 3, ind) = fwarp(2, 3, ind) * scaleLocal
        enddo
      enddo
    endif
    !
    ! if the input data is a subset in X or Y, subtract starting
    ! coordinates from ixStartWarp and iyStartWarp
    !
    ixStartWarp = ixStartWarp - ixSubset
    iyStartWarp = iyStartWarp - iySubset
  endif
  !
  ! Done with array in its small form and with angReproj
  deallocate(array, angReproj, ivExclude, ivReproj, stat = ierr)
  !
  ! Here is the place to project model points and exit
  if (projModel) call projectModel(outputFile, angleOutput, transformFile, delta, &
      numViewOrig, defocusFile, pixForDefocus, focusInvert)
  !
  ! If reprojecting, set the pointers and return
  if (recReproj) then
    indOutSlice = 1
    minXload = minXreproj
    maxXload = maxXreproj
    if (nxWarp .ne. 0) then
      minXload = max(1, minXload - 100)
      maxXload = min(nrecXyz(1), maxXload + 100)
    endif
    iwideReproj = maxXload + 1 - minXload
    isliceSizeBP = iwideReproj * ithickReproj
    inPlaneSize = isliceSizeBP
    if (nxWarp .ne. 0) then
      dxWarpDelz = idelXwarp / 2.
      numWarpDelz = max(2., (iwideReproj - 1) / dxWarpDelz) + 1
      dxWarpDelz = (iwideReproj - 1.) / (numWarpDelz - 1.)
    endif
    !
    ! Get projection offsets for getting from coordinate in reprojection
    ! to coordinate in original projections.  The X coordinate must account
    ! for the original offset in building the reconstruction plus any
    ! additional offset.  But the line number here is the line # in the
    ! reconstruction so we only need to adjust Y by the original starting
    ! line.  Also replace the slice limits.
    xprojOffset = minXreproj - 1 + nprojXyz(1) / 2 - iwidth / 2 - xOffset
    yprojOffset = isliceStart - 1
    isliceStart = minZreproj
    isliceEnd = maxZreproj
    iwidth = noutXyz(1)
    nyProj = nrecXyz(3)
    dmeanIn = (dmeanIn / outScale - outAdd) / filterScale
    if (threshPolarity .ne. 0.) then
      threshForReproj = (threshForReproj / outScale - outAdd) / filterScale
      threshMarkVal = threshForReproj
      threshFillVal = dmeanIn
      if (threshSumFac < 1.) threshMarkVal = threshMarkVal * ithickReproj
    endif
    indLoadBase = 1
    ipExtraSize = 0
    numPad = 0
    if (debug) print *,'scale: ', outScale, outAdd

    numNeedEval = min(numNeedEval, isliceEnd + 1 - isliceStart)
    call setNeededSlices(maxNeeds, numNeedEval)
    if (allocateArray(maxNeeds, numNeedEval, 1, minMemory) == 0) &
        call exitError('THE MAIN ARRAY CANNOT BE ALLOCATED LARGE ENOUGH'// &
        ' TO REPROJECT A SINGLE Y VALUE')
    allocate(reprojLines(iwidth * numPlanes), stat = ierr)
    call memoryErrorUC(ierr, 'ARRAY FOR REPROJECTED LINES')
    !
    if (projSubtraction) then
      allocate(origLines(iwidth * numPlanes), stat = ierr)
      call memoryErrorUC(ierr, 'ARRAY FOR ORIGINAL LINES')
    endif
    if (useGPU) then
      ind = 0
      if (debug) ind = 1
      useGPU = gpuAvailable(indGPU, gpuMemory, maxTex2D, maxTexLayer, maxTex3D, ind) &
          .ne. 0
      if (.not. useGPU) gpuErrorStr = 'No GPU is available'
      if (.not. useGPU .and. .not. debug) gpuErrorStr = gpuErrorStr// &
          ', run gputilttest for more details'
      ind = maxNeeds(1) * inPlaneSize + iwidth * numPlanes
      iex = iwidth * numPlanes
      kti = 0
      if (useGPU .and. nxWarp > 0) then
        ind = ind + maxNeeds(1) * (8 * iwidth + numWarpDelz) + &
            12 * numWarpPos * numViews
        iex = iex + 12 * numWarpPos * numViews
        kti = numWarpDelz
        call packLocalData()
      endif
      if (useGPU) then
        useGPU = 4 * ind <= gpuMemoryFrac * gpuMemory
        if (.not. useGPU) gpuErrorStr = 'GPU is available but it has insufficient '// &
            'memory to reproject with current parameters'
      endif
      if (useGPU)  call allocateGpuPlanes(iex, nxWarp * nyWarp, kti, 0, &
          numPlanes, iwideReproj, ithickReproj)
      if (useGPU .and. nxWarp .ne. 0) then
        useGPU = gpuLoadLocals(packLocal, nxWarp * nyWarp) == 0
        if (.not. useGPU) gpuErrorStr = 'Failed to load local alignment data into GPU'
      endif
      if (useGPU) print *,'Using the GPU for reprojection'
      if (allocated(packLocal)) deallocate(packLocal)
      call warnOrExitIfNoGPU()
      if (.not. useGPU) print *,'The GPU cannot be used, using the CPU for reprojection'
    endif
    !
    ! Finally allocate the warpDelz now that number of lines is known,
    ! and projecton factors now that number of planes is known
    if (nxWarp .ne. 0) then
      kti = iwideReproj * numPlanes
      allocate(warpDelz(numWarpDelz * max(1, numGpuPlanes)), xprojfs(kti), &
          xprojzs(kti), yprojfs(kti), yprojzs(kti), stat = ierr)
      call memoryErrorUC(ierr, 'LOCAL PROJECTION FACTOR OR warpDelz ARRAYS')
    endif
    return
  endif
  !
  ! BACKPROJECTION ONLY.  First allocate the projection factor array
  if (nxWarp .ne. 0) then
    allocate(xprojfs(iwidth), xprojzs(iwidth), yprojfs(iwidth), &
        yprojzs(iwidth), stat = ierr)
    call memoryErrorUC(ierr, 'ARRAYS FOR LOCAL PROJECTION FACTORS')
  endif
  !
  ! If reading base, figure out total views being added and adjust scales
  if (readBaseRec .and. .not. recSubtraction) then
    iv = numViews
    do j = 1, numViews
      k = 0
      do i = 1, numViewSubtract
        if (iviewSubtract(i) == 0 .or. mapUsedView(j) == iviewSubtract(i)) k = 1
      enddo
      iv = iv - 2 * k
    enddo
    baseOutScale = outScale / numViewBase
    baseOutAdd = outAdd * numViewBase
    outScale = outScale / (iv + numViewBase)
    outAdd = outAdd * (iv + numViewBase)
    if (debug)  print *,'base: ', baseOutScale, baseOutAdd
  else if (numSIRTiter <= 0) then
    outScale = outScale / numViews
    outAdd = outAdd * numViews
  endif
  if (debug) print *,'scale: ', outScale, outAdd
  numNeedEval = min(numNeedEval, isliceEnd + 1 - isliceStart)
  call setNeededSlices(maxNeeds, numNeedEval)
  if (debug) print *,(maxNeeds(i), i = 1, numNeedEval)
  if (iterForReport > 0) then
    allocate(reportVals(3, max(1, numSIRTiter)), stat = ierr)
    call memoryErrorUC(ierr, 'REPORT VALUE ARRAY')
    reportVals(1:3, 1:max(1, numSIRTiter)) = 0.
  endif
  !
  ! 12/13/09: removed fast backprojection code
  !
  ! Set up padding: 10% of X size or minimum of 16, max of 50
  numPadTmp = min(50, 2 * max(8, nprojXyz(1) / 20))
  ! npadtmp = 2 * nprojXyz(1)
  nprojPad = niceFrame(2 * ((nprojXyz(1) + numPadTmp) / 2), 2, niceFFTlimit())
  numPad = nprojPad - nprojXyz(1)
  !
  ! Set up defaults for plane size and start of planes of input data
  isliceSizeBP = iwidth * ithickBP
  ipExtraSize = 0
  nxProjPad = nxProj + 2 + numPad
  inPlaneSize = nxProjPad * numViews
  indOutSlice = inPlaneSize + 1
  indLoadBase = indOutSlice + isliceSizeBP
  if (numSIRTiter > 0) then
    if (sirtFromZero) indOutSlice = 2 * inPlaneSize + 1
    indLoadBase = indOutSlice + isliceSizeBP
    ireadBase = indLoadBase
    indWorkPlane = indLoadBase + isliceSizeBP
    indLoadBase = indWorkPlane + inPlaneSize
    nvsXyz(1) = iwidth
    nvsXyz(2) = isliceEnd + 1 - isliceStart
    nvsXyz(3) = numViews
    if (ifOutSirtProj > 0) then
      call imOpen(4, 'sirttst.prj', 'NEW')
      call iiuCreateHeader(4, nvsXyz, nvsXyz, 2, title, 0)
      call iiuWriteHeader(4, title, 0, -1.e6, 1.e6, 0.)
    endif
    nvsXyz(2) = ithickBP
    nvsXyz(3) = isliceEnd + 1 - isliceStart
    if (ifOutSirtRec > 0) then
      call imOpen(5, 'sirttst.drec', 'NEW')
      call iiuCreateHeader(5, nvsXyz, nvsXyz, 2, title, 0)
      call iiuWriteHeader(5, title, 0, -1.e6, 1.e6, 0.)
    endif
  endif
  maxStack = 0
  !
  ! Determine if GPU can be used, but don't try to allocate yet
  if (useGPU) then
    ind = 0
    if (debug) ind = 1
    wallStart = wallTime()
    useGPU = gpuAvailable(indGPU, gpuMemory, maxTex2D, maxTexLayer, maxTex3D, ind) .ne. 0
    if (debug) write(*,'(a,f8.4)') 'Time to test if GPU available: ', &
        wallTime() - wallStart
    if (.not. useGPU) gpuErrorStr =  &
        'No GPU is available, run gputilttest for more details'
    if (useGPU) then
      !
      ! Basic need is input planes for reconstructing one slice plus 2
      ! slices for radial filter and planes being filtered, plus output
      ! slice. Local alignment adds 4 arrays for local proj factors
      iv = (maxNeeds(1) + 2) * inPlaneSize + isliceSizeBP
      if (nxWarp .ne. 0) &
          iv = iv + 4 * (iwidth * numViews + 12 * numWarpPos * numViews)
      if (numSIRTiter > 0) iv = iv + isliceSizeBP + inPlaneSize
      if (sirtFromZero) iv = iv + inPlaneSize
      useGPU = 4 * iv <= gpuMemoryFrac * gpuMemory
      if (useGPU) then
        interpFacStretch = 0
      else
        gpuErrorStr = 'GPU is available but it has insufficient '// &
            'memory to backproject with current parameters'
      endif
    endif
    call warnOrExitIfNoGPU()
    if (.not. useGPU) print *, 'The GPU cannot be used; using CPU for backprojection'
  endif
  !
  ! next evaluate cosine stretch and new-style tilt for memory
  !
  if (ifAlpha < 0) then
    !
    ! new-style X-axis tilt
    !
    indLoadBase = indOutSlice + isliceSizeBP * (numVertNeeded + 1)
    if (numSIRTiter > 0) then
      ireadBase = indLoadBase
      indWorkPlane = indLoadBase + numReadNeed * ithickOut * iwidth
      indLoadBase = indWorkPlane + inPlaneSize
    endif
    !
    ! find out what cosine stretch adds if called for
    !
    if (interpFacStretch > 0) then
      call set_cos_stretch()
      inPlaneSize = max(inPlaneSize, indStretchLine(numViews + 1))
      ipExtraSize = inPlaneSize
    endif
    !
    ! Does everything fit?  If not, drop back to old style tilting
    !
    if (allocateArray(maxNeeds, numNeedEval, 4, minMemory) == 0) then
      if (numSIRTiter > 0) call exitError( &
          'ALLOCATING ARRAYS NEEDED FOR IN-MEMORY SIRT ITERATIONS')
      ifAlpha = 1
      ithickBP = ithickOut
      isliceSizeBP = iwidth * ithickBP
      ycenOut = ithickBP / 2 + 0.5 + yOffset
      indLoadBase = indOutSlice + isliceSizeBP
      ipExtraSize = 0
      inPlaneSize = nxProjPad * numViews
      write(*,63)
63    format(/,'Failed to allocate an array big enough ', &
          'to use new-style X-axis tilting')
      call setNeededSlices(maxNeeds, numNeedEval)
    endif
  endif
  !
  ! If not allocated yet and not warping, try cosine stretch here
  !
  if (maxStack == 0 .and. nxWarp == 0 .and. interpFacStretch > 0) then
    call set_cos_stretch()
    !
    ! set size of plane as max of loading size and stretched size
    ! also set that an extra plane is needed
    ! if there is not enough space for the planes needed, then
    ! disable stretching and drop back to regular code
    !
    inPlaneSize = max(inPlaneSize, indStretchLine(numViews + 1))
    ipExtraSize = inPlaneSize
    if (allocateArray(maxNeeds, numNeedEval, 1, minMemory) == 0) then
      ipExtraSize = 0
      inPlaneSize = nxProjPad * numViews
      interpFacStretch = 0
      write(*,62)
62    format(/,'Failed to allocate an array big enough ', &
          'to use cosine stretching')
    endif
  endif
  !
  ! If array still not allocated (failure, or local alignments), do it now
  if (maxStack == 0) then
    if (allocateArray(maxNeeds, numNeedEval, 1, minMemory) == 0) &
        call exitError('COULD NOT ALLOCATE MAIN ARRAY LARGE ENOUGH TO '// &
        'RECONSTRUCT A SINGLE SLICE')
  endif
  !
  ! Set up radial weighting
  if (numSIRTiter > 0 .and. .not.sirtFromZero) flatFrac = 2.
  call radialWeights(iradMax, radFall, 1, ifMultByGaussian)
  if (sirtFromZero) then
    frac = zeroWeight
    flatFrac = 2.
    call radialWeights(iradMax, radFall, 2, ifMultByGaussian)
    zeroWeight = frac
  endif
  !
  ! If Using GPU, make sure memory is OK now and allocate and load things
  ind = isliceSizeBP + 2 * inPlaneSize
  if (useGPU .and. nxWarp .ne. 0) then
    ind = ind + 4 * numViews * iwidth + 12 * numWarpPos * numViews
    call packLocalData()
  endif
  wallStart = wallTime()
  if (useGPU) then
    if (ifAlpha <= 0 .and. nxWarp == 0) then
      iv = 0
      j = 1
      if (numSIRTiter > 0) iv = numViews
      if (sirtFromZero) j = 2
      useGPU = gpuAllocArrays(iwidth, ithickBP, nxProjPad, numViews, 1, numViews, 0, &
          0, j, iv, 1, 1, 0) == 0
    else
      call allocateGpuPlanes(ind, nxWarp * nyWarp, 0, 1, ithickBP, nxProjPad, numViews)
    endif
    if (debug .and. useGPU) write(*,'(a,f8.4)') 'Time to allocate on GPU: ', &
        wallTime() - wallStart

    ! print *,useGPU
    wallStart = wallTime()
    if (useGPU) then
      useGPU = gpuLoadFilter(array) == 0
      if (.not. useGPU) gpuErrorStr = 'Failed to load filter array into GPU'
    endif
    ! print *,useGPU
    if (useGPU .and. nxWarp .ne. 0) then
      useGPU = gpuLoadLocals(packLocal, nxWarp * nyWarp) == 0
      if (.not. useGPU) gpuErrorStr = 'Failed to load local alignment data into GPU'
    endif
    ! print *,useGPU
    if (debug .and. useGPU) write(*,'(a,f8.4)') 'Time to load filter/locals on GPU: ', &
        wallTime() - wallStart
    if (allocated(packLocal)) deallocate(packLocal)
    call warnOrExitIfNoGPU()
    if (useGPU) then
      print *,'Using GPU for backprojection'
    else
      print *, 'The GPU cannot be used - using CPU for backprojection'
    endif
  endif

  ! print *,interpFacStretch, ipExtraSize, ifAlpha, numVertNeeded
  !
  return
  !
2410 call exitError('READING LOCAL TILT ALIGNMENT DATA FROM FILE')
2411 call exitError('READING TILT ANGLES FROM FILE')
2412 call exitError('READING X-AXIS TILT ANGLES FROM FILE')
2413 call exitError('READING Z FACTORS FROM FILE')
2414 call exitError('READING WEIGHTING FACTORS FROM FILE')
2415 call exitError('READING ANGLES FOR WEIGHTING FROM FILE')
  !
  !
48 format(/,/,1xx,78('-'))
49 format('TILT: ',a,t57,a9,2xx,a8)
51 format(/,' Projection angles:',//,(8f9.2))
52 format(/,/,1xx,78('-'))
53 format(/,'Scaling of local alignments by ',f8.3, ' determined from pixel sizes')
101 format(/,' Title:    ',20a4)
201 format(/,' Rows',i5,' to',i5,', (at intervals of',i4,') of the' &
      ,' projection planes will be reconstructed.')
301 format(/,' Thickness of reconstructed slice is',i5, ' pixels.')
401 format(/,' Mask applied to edges of output slices with',i5, ' extra pixels masked')
501 format(/,' Radial weighting function Gaussian starts at cutoff =',f7.4, '  sigma =', &
        f7.4)
502 format(/,' Multiplying radial function by Gaussian')
503 format(/,' Doing Hamming-like filter by multiplying by Gaussian with cutoff =',f7.4, &
        '  sigma =', f7.4)
601 format(/,' Output map rotated by',f6.1,' degrees about tilt axis', &
      ' with respect to tilt origin', /, ' Tilt axis displaced by',f9.2, &
      ' pixels from centre of projection')
701 format(/,' Output map densities incremented by',f8.2, ' and then multiplied by',f8.2)
801 format(/,' Output map is sectioned perpendicular to the ' ,'tilt axis')
901 format(/,' Output map is sectioned parallel to the' &
      ,' zero tilt projection, inverting handedness')
902 format(/,' Output map is sectioned parallel to the' &
      ,' zero tilt projection, retaining handedness')
1001 format(/,' Data mode of output file is',i3)
1301 format(/,' Taking logarithm of input data plus',f10.3)
1701 format(/,' Compression was confined to',f6.3, &
      ' of the distance over which it was measured')
1801 format(/,' Weighting by tilt density computed to distance of' &
      ,i3,' views',/,'  weighting factors:',(10f6.3))
1802 format(/,' No weighting by tilt density')
2001 format(/,' Width of reconstruction is',i6,' pixels')
2101 format(/,' Output slice shifted up',f7.1,' and to right',f7.1, ' pixels')
2201 format(/,' Alpha tilting to be applied with angles from file')
2202 format(/,' Constant alpha tilt of',f6.1,' to be applied based on ', &
         'angles from file')
2301 format(/,' Global alpha tilt of',f6.1,' will be applied')
2401 format(/,' Local tilt alignment information read from file')
2501 format(/,' Local alignment positions and shifts reduced by',f7.4)
2601 format(/,' Full aligned stack will be assumed to be',i6,' by', i6,' pixels')
2701 format(/,' Aligned stack will be assumed to be a subset starting at',2i6)
2801 format(/,' Cosine stretching, if any, will have interpolation order', i2, &
         ', sampling factor',i2)
3001 format(/,' X-tilting with vertical slices, if any, will have interpolation order', &
         i2)
3101 format(/,' Output will be one or more reprojections')
3201 format(/,' Z-dependent shifts to be applied with factors from file')
3301 format(/,' Dimensions and coordinates will be scaled down by a factor of ',i2)
3401 format(/,' Computed slices are part of a total volume from slice', i6,' to',i6)
3501 format(/,i4,' iterations of SIRT algorithm will be done')

CONTAINS

  !
  ! Allocate as many planes as possible on the GPU up to the number
  ! allowed in ARRAY; nonPlane gives the number of floats needed for
  ! other data and numWarps is number of local positions, numDelz is
  ! size of warpDelz array, numfilts is number of filter
  ! filter arrays needed, nygout is size of output in y, nxGPlane and
  ! nxyplane are size of input data
  !
  subroutine allocateGpuPlanes(nonPlane, numWarps, numDelz, numFilts, &
      nygout, nxGPlane, nyGPlane)
    integer*4 nonPlane, numWarps, maxGpuPlane, nygout, nxGPlane, nyGPlane
    integer*4 numFilts, numDelz, iperPlane, indEval3D

    ! Evaluate texture types in this order: 3D, layered, 2D.
    ! But for no X tilt/Z-factors or locals, only 2D textures are allowed
    ! And for no locals or no reprojection, only 2D layers or 2D textures are allowed
    ! (Unimplemented options had no advantage)
    if (ifAlpha <= 0 .and. ifZfactors == 0 .and. nxWarp == 0) then
      indEval3D = 2
    else if (nxWarp == 0 .and. .not. recReproj) then
      indEval3D = 1
    else
      indEval3D = 0
    endif
    if (debug) then
      print *, 'In allocateGpuPlanes:', ifAlpha,ifZfactors, nxWarp, recReproj, indEval3D
      print *, nxGPlane, maxTex3D(1), nyGPlane, maxTex3D(2), maxNeeds(1), maxTex3D(3)
    endif
    !
    ! Check that the test entry was not for an unimplemented type 
    if ((indEval3D > 0 .and. if3Dtexture > 0) .or. (indEval3D > 1 .and.  &
        if3Dtexture > -999 .and. if3Dtexture .ne. 0)) call exitError( &
        'Entry for TextureType not legal with current type of computation')
    !
    ! If 3D is allowed and fits, take it; then if Layered is allowed and fits, take it
    ! then fall back to 2D
    if (if3Dtexture == -999 .and. indEval3d == 0 .and. nxGPlane < maxTex3D(1) .and. &
        nyGPlane < maxTex3D(2) .and. maxNeeds(1) < maxTex3D(3)) if3Dtexture = 1
    if (if3Dtexture == -999 .and. indEval3d <= 1 .and. nxGPlane < maxTexLayer(1) .and. &
        nyGPlane < maxTexLayer(2) .and. maxNeeds(1) < maxTexLayer(3)) if3Dtexture = -1
    if (if3Dtexture == -999) if3Dtexture = 0
    if (debug) print *,'Proceeding to allocate on GPU with texture type ', if3Dtexture
    
    !
    ! Start with as many planes as possible but no more than in array
    ! and if doing 2D textures, no more lines for array on GPU than the Y limit for
    ! 2D textures (32767 before Fermi, typically 65536)
    iperPlane = inPlaneSize
    if (numDelz > 0) iperPlane = inPlaneSize + 8 * iwidth + numWarpDelz
    maxGpuPlane = (gpuMemoryFrac * gpuMemory / 4. - nonPlane) / iperPlane
    if (maxGpuPlane < maxNeeds(1)) then
      write(gpuErrorStr, '(a,i4,a,i4,a)') 'The GPU only has enough memory to load', &
          maxGpuPlane, ' planes of data and, with current parameters, up to',  &
          maxNeeds(1), ' input planes are required to reconstruct a single plane'
    else if (if3Dtexture == 0 .and. maxTex2D(2) / nyGPlane < maxNeeds(1)) then
      write(gpuErrorStr, '(a,i4,a,i4,a)') 'With current parameters, up to', maxNeeds(1), &
          ' input planes are required to reconstruct a single plane, while '// &
          'the allowed number of lines in GPU texture memory will fit only', &
          maxTex2D(2) / nyGPlane, ' input planes'
    endif
    !
    ! Limit maximum planes by the type allowed for the texture
    if (if3Dtexture > 0) then
      maxGpuPlane = min(maxGpuPlane, numPlanes, maxTex3D(3))
    else if (if3Dtexture < 0) then
      maxGpuPlane = min(maxGpuPlane, numPlanes, maxTexLayer(3))
    else
      maxGpuPlane = min(maxGpuPlane, numPlanes, maxTex2D(2) / nyGPlane)
    endif

    ! FOR TESTING MEMORY SHIFTING ETC
    ! ind = max(maxNeeds(1), min(ind, numPlanes / 3))
    numGpuPlanes = 0
    ! print *, ind, numPlanes, maxNeeds(1), maxGpuPlane
    do i = maxGpuPlane, maxNeeds(1), -1
      if (gpuAllocArrays(iwidth, nygout, nxGPlane, nyGPlane, i, numViews, &
          numWarps,  numDelz, numFilts, 0, maxGpuPlane, maxNeeds(1), if3Dtexture) == 0) &
          then
        numGpuPlanes = i
        loadGpuStart = 0
        loadGpuEnd = 0
        exit
      endif
    enddo
    useGPU = numGpuPlanes > 0
    if (.not. useGPU .and. maxGpuPlane >= maxNeeds(1)) gpuErrorStr = 'Failed to '// &
        'allocate enough memory on the GPU to recontruct even a single plane'
    return
  end subroutine allocateGpuPlanes

  ! Issue desired warning or exit on error if GPU not available
  subroutine warnOrExitIfNoGPU()
    if (useGPU) return
    if ((ifGpuByEnviron .ne. 0 .and. iactGpuFailEnviron == 2) .or.  &
        (ifGpuByEnviron == 0 .and. iactGpuFailOption == 2)) then
      write(*, '(/,a,a)') 'ERROR: tilt - ', trim(gpuErrorStr)
      call exitError('Use of the GPU was requested but a GPU cannot be used')
    endif
    write(*, '(a)') trim(gpuErrorStr)
    if ((ifGpuByEnviron .ne. 0 .and. iactGpuFailEnviron == 1) .or.  &
        (ifGpuByEnviron == 0 .and. iactGpuFailOption == 1)) write(*,'(a)') &
        'MESSAGE: Use of the GPU was requested but a GPU will not be used'
    call flush(6)
    return
  end subroutine warnOrExitIfNoGPU


  subroutine packLocalData()
    allocate(packLocal(numWarpPos * numViews, 12), stat = ierr)
    if (ierr .ne. 0) then
      useGPU = .false.
      gpuErrorStr = 'Failed to allocate array needed for using GPU with local alignments'
      call warnOrExitIfNoGPU()
    else
      !
      ! Pack data into one array
      do ipos = 1, nxWarp * nyWarp
        do iv = 1, numViews
          i = indWarp(ipos) + mapUsedView(iv)
          j = (ipos - 1) * numViews + iv
          packLocal(j, 1) = fwarp(1, 1, i)
          packLocal(j, 2) = fwarp(2, 1, i)
          packLocal(j, 3) = fwarp(1, 2, i)
          packLocal(j, 4) = fwarp(2, 2, i)
          packLocal(j, 5) = fwarp(1, 3, i)
          packLocal(j, 6) = fwarp(2, 3, i)
          packLocal(j, 7) = cWarpAlpha(i)
          packLocal(j, 8) = sWarpAlpha(i)
          packLocal(j, 9) = cWarpBeta(i)
          packLocal(j, 10) = sWarpBeta(i)
          packLocal(j, 11) = warpXZfac(i)
          packLocal(j, 12) = warpYZfac(i)
        enddo
      enddo
    endif
    return
  end subroutine packLocalData

  ! END OF INPUT ROUTINE
end subroutine inputParameters


! Finds two nearest angles to PROJANGLE and returns their indices IND1 and
! IND2 and an interpolation fraction FRAC
!
subroutine lookupAngle(projAngle, angles, numViews, ind1, ind2, frac)
  implicit none
  real*4 projAngle, angles(*), frac
  integer*4 numViews, ind1, ind2, i
  ind1 = 0
  ind2 = 0
  do i = 1, numViews
    if (projAngle >= angles(i)) then
      if (ind1 == 0) ind1 = i
      if (angles(i) > angles(ind1)) ind1 = i
    else
      if (ind2 == 0) ind2 = i
      if (angles(i) < angles(ind2)) ind2 = i
    endif
  enddo
  frac = 0.
  if (ind1 == 0) then
    ind1 = ind2
  elseif (ind2 == 0) then
    ind2 = ind1
  else
    frac = (projAngle - angles(ind1)) / (angles(ind2) - angles(ind1))
  endif
  return
end subroutine lookupAngle


! Compute mean and SD of interior of a slice for SIRT report
!
subroutine sampleForReport(slice, lslice, kthick, iteration, sampScale, sampAdd)
  use tiltvars
  implicit none
  real*4 slice(*), sampScale, sampAdd, avg, SD
  integer*4 lslice, kthick, iteration, iskip, ixLow, iyLow, ierr
  integer*4 sampleMeanSD
  !
  iskip = (isliceEnd - isliceStart) / 10
  if (lslice < isliceStart + iskip .or. lslice > isliceEnd - iskip) return
  ixLow = iwidth / 10
  iyLow = kthick / 4
  ierr = sampleMeanSD(slice, iwidth, kthick, 0.02, ixLow, iyLow, iwidth - 2 * &
      ixLow, kthick - 2 * iyLow, avg, SD)
  if (ierr .ne. 0) return
  reportVals(1, iteration) = reportVals(1, iteration) + &
      (avg + sampAdd) * sampScale
  reportVals(2, iteration) = reportVals(2, iteration) + SD * sampScale
  reportVals(3, iteration) = reportVals(3, iteration) + 1
  return
end subroutine sampleForReport


! return indices to four local areas, and fractions to apply for each
! at location ix, iy in view iv, where ix and iy are indexes in
! the reconstruction adjusted to match coordinates of projections
!
subroutine local_factors(ix, iy, iv, ind1, ind2, ind3, ind4, f1, f2, f3, f4)
  !
  use tiltvars
  implicit none
  integer*4 ix, iy, iv, ind1, ind2, ind3, ind4, ixt, ixPos, iyt, iyPos
  real*4 f1, f2, f3, f4, fx, fy
  !
  ixt = min(max(ix - ixStartWarp, 0), (nxWarp - 1) * idelXwarp)
  ixPos = min(ixt / idelXwarp + 1, nxWarp - 1)
  fx = float(ixt - (ixPos - 1) * idelXwarp) / idelXwarp
  iyt = min(max(iy - iyStartWarp, 0), (nyWarp - 1) * idelYwarp)
  iyPos = min(iyt / idelYwarp + 1, nyWarp - 1)
  fy = float(iyt - (iyPos - 1) * idelYwarp) / idelYwarp

  ind1 = indWarp(nxWarp * (iyPos - 1) + ixPos) + iv
  ind2 = indWarp(nxWarp * (iyPos - 1) + ixPos + 1) + iv
  ind3 = indWarp(nxWarp * iyPos + ixPos) + iv
  ind4 = indWarp(nxWarp * iyPos + ixPos + 1) + iv
  f1 = (1. -fy) * (1. -fx)
  f2 = (1. -fy) * fx
  f3 = fy * (1. -fx)
  f4 = fy * fx
  return
end subroutine local_factors


! Compute local projection factors at a position in a column for view iv:
! j is X index in the reconstruction, lslice is slice # in aligned stack
!
subroutine localProjFactors(j, lslice, iv, xprojFix, xprojZ, yprojFix, yprojZ)
  use tiltvars
  implicit none
  integer*4 j, lslice, iv
  real*4 xprojFix, xprojZ, yprojFix, yprojZ
  integer*4 ind1, ind2, ind3, ind4, ixc
  real*4 f1, f2, f3, f4, xx, yy
  real*4 cosAlph, sinAlph, a11, a12, a21, a22, xAdd, yAdd, xAllAdd, yAllAdd
  real*4 cosAlph2, sinAlph2, a112, a122, a212, a222, xAdd2, yAdd2
  real*4 cosAlph3, sinAlph3, a113, a123, a213, a223, xAdd3, yAdd3
  real*4 cosAlph4, sinAlph4, a114, a124, a214, a224, xAdd4, yAdd4
  real*4 f1x, f2x, f3x, f4x, f1xy, f2xy, f3xy, f4xy
  real*4 f1y, f2y, f3y, f4y, f1yy, f2yy, f3yy, f4yy
  real*4 xp1f, xp1z, yp1f, xp2f, xp2z, yp2f, xp3f, xp3z, yp3f, xp4f, xp4z, yp4f
  real*4 cosBet, sinBet, cosBet2, sinBet2, cosBet3, sinBet3, cosBet4, sinBet4
  !
  ! get transform and angle adjustment
  !
  ixc = nint(j - xcenOut + xcenIn + axisXoffset)
  call local_factors(ixc, lslice, mapUsedView(iv), ind1, ind2, ind3, ind4, f1, f2, f3, f4)
  !
  ! get all the factors needed to compute a projection position
  ! from the four local transforms
  !
  cosBet = cWarpBeta(ind1)
  sinBet = sWarpBeta(ind1)
  cosAlph = cWarpAlpha(ind1)
  sinAlph = sWarpAlpha(ind1)
  a11 = fwarp(1, 1, ind1)
  a12 = fwarp(1, 2, ind1)
  a21 = fwarp(2, 1, ind1)
  a22 = fwarp(2, 2, ind1)
  xAdd = fwarp(1, 3, ind1) + xcenIn - xcenIn * a11 - centerSlice * a12
  yAdd = fwarp(2, 3, ind1) + centerSlice - xcenIn * a21 - centerSlice * a22
  !
  cosBet2 = cWarpBeta(ind2)
  sinBet2 = sWarpBeta(ind2)
  cosAlph2 = cWarpAlpha(ind2)
  sinAlph2 = sWarpAlpha(ind2)
  a112 = fwarp(1, 1, ind2)
  a122 = fwarp(1, 2, ind2)
  a212 = fwarp(2, 1, ind2)
  a222 = fwarp(2, 2, ind2)
  xAdd2 = fwarp(1, 3, ind2) + xcenIn - xcenIn * a112 - centerSlice * a122
  yAdd2 = fwarp(2, 3, ind2) + centerSlice - xcenIn * a212 - centerSlice * a222
  !
  cosBet3 = cWarpBeta(ind3)
  sinBet3 = sWarpBeta(ind3)
  cosAlph3 = cWarpAlpha(ind3)
  sinAlph3 = sWarpAlpha(ind3)
  a113 = fwarp(1, 1, ind3)
  a123 = fwarp(1, 2, ind3)
  a213 = fwarp(2, 1, ind3)
  a223 = fwarp(2, 2, ind3)
  xAdd3 = fwarp(1, 3, ind3) + xcenIn - xcenIn * a113 - centerSlice * a123
  yAdd3 = fwarp(2, 3, ind3) + centerSlice - xcenIn * a213 - centerSlice * a223
  !
  cosBet4 = cWarpBeta(ind4)
  sinBet4 = sWarpBeta(ind4)
  cosAlph4 = cWarpAlpha(ind4)
  sinAlph4 = sWarpAlpha(ind4)
  a114 = fwarp(1, 1, ind4)
  a124 = fwarp(1, 2, ind4)
  a214 = fwarp(2, 1, ind4)
  a224 = fwarp(2, 2, ind4)
  xAdd4 = fwarp(1, 3, ind4) + xcenIn - xcenIn * a114 - centerSlice * a124
  yAdd4 = fwarp(2, 3, ind4) + centerSlice - xcenIn * a214 - centerSlice * a224
  !
  f1x = f1 * a11
  f2x = f2 * a112
  f3x = f3 * a113
  f4x = f4 * a114
  f1xy = f1 * a12
  f2xy = f2 * a122
  f3xy = f3 * a123
  f4xy = f4 * a124
  ! fxfromy=f1*a12+f2*a122+f3*a123+f4*a124
  f1y = f1 * a21
  f2y = f2 * a212
  f3y = f3 * a213
  f4y = f4 * a214
  f1yy = f1 * a22
  f2yy = f2 * a222
  f3yy = f3 * a223
  f4yy = f4 * a224
  ! fyfromy=f1*a22+f2*a222+f3*a223+f4*a224
  xAllAdd = f1 * xAdd + f2 * xAdd2 + f3 * xAdd3 + f4 * xAdd4
  yAllAdd = f1 * yAdd + f2 * yAdd2 + f3 * yAdd3 + f4 * yAdd4
  !
  ! Each projection position is a sum of a fixed factor ("..f")
  ! and a factor that multiplies z ("..z")
  !
  xx = j - xcenOut
  yy = lslice - centerSlice
  xp1f = xx * cosBet + yy * sinAlph * sinBet + xcenIn + axisXoffset
  xp1z = cosAlph * sinBet + warpXZfac(ind1)
  xp2f = xx * cosBet2 + yy * sinAlph2 * sinBet2 + xcenIn + axisXoffset
  xp2z = cosAlph2 * sinBet2 + warpXZfac(ind2)
  xp3f = xx * cosBet3 + yy * sinAlph3 * sinBet3 + xcenIn + axisXoffset
  xp3z = cosAlph3 * sinBet3 + warpXZfac(ind3)
  xp4f = xx * cosBet4 + yy * sinAlph4 * sinBet4 + xcenIn + axisXoffset
  xp4z = cosAlph4 * sinBet4 + warpXZfac(ind4)

  yp1f = yy * cosAlph + centerSlice
  yp2f = yy * cosAlph2 + centerSlice
  yp3f = yy * cosAlph3 + centerSlice
  yp4f = yy * cosAlph4 + centerSlice
  !
  ! store the fixed and z-dependent component of the
  ! projection coordinates
  !
  xprojFix = f1x * xp1f + f2x * xp2f + f3x * xp3f + f4x * xp4f + &
      f1xy * yp1f + f2xy * yp2f + f3xy * yp3f + f4xy * yp4f + xAllAdd
  xprojZ = f1x * xp1z + f2x * xp2z + f3x * xp3z + f4x * xp4z - &
      (f1xy * (sinAlph - warpYZfac(ind1)) + f2xy * (sinAlph2 - warpYZfac(ind2)) + &
      f3xy * (sinAlph3 - warpYZfac(ind3)) + f4xy * (sinAlph4 - warpYZfac(ind4)))
  yprojFix = f1y * xp1f + f2y * xp2f + f3y * xp3f + f4y * xp4f + &
      f1yy * yp1f + f2yy * yp2f + f3yy * yp3f + f4yy * yp4f + yAllAdd
  yprojZ = f1y * xp1z + f2y * xp2z + f3y * xp3z + f4y * xp4z - &
      (f1yy * (sinAlph - warpYZfac(ind1)) + f2yy * (sinAlph2 - warpYZfac(ind2)) + &
      f3yy * (sinAlph3 - warpYZfac(ind3)) + f4yy * (sinAlph4 - warpYZfac(ind4)))
  return
end subroutine localProjFactors


! This is the former code for assessing backprojection positions,
! new method using localProjFactors verified to give the same result
! But this shows how the BP position can be computed directly
!!$  call local_factors (ixsam, itry, mapUsedView(iv), ind1, ind2, &
!!$      ind3, ind4, f1, f2, f3, f4)
!!$c
!!$c   for each position, find back - projection location
!!$c   transform if necessary, and use to get min and
!!$c   max slices needed to get this position
!!$c
!!$  xx = ixsam - xcenOut
!!$  yy = itry - centerSlice
!!$  zz = iy - ycenOut
!!$c   Global bp position
!!$  xp = xx * cosBeta(iv) + yy * sinAlpha(iv) * sinBeta(iv) + zz * (cosAlpha(iv) * sinBeta(iv) + xzfac(iv)) + &
!!$      xcenIn + axisXoffset
!!$  yp = yy * cosAlpha(iv) - zz * (sinAlpha(iv) - yzfac(iv)) + centerSlice
!!$c   Local position:
!!$  xp = xx * cWarpBeta(ind1) + yy * sWarpAlpha(ind1) * sWarpBeta(ind1) + &
!!$      zz * (cWarpAlpha(ind1) * sWarpBeta(ind1) + warpXZfac(ind1)) + xcenIn + axisXoffset
!!$  yp = yy * cWarpAlpha(ind1) - zz * (sWarpAlpha(ind1) - warpYZfac(ind1)) + centerSlice
!!$  call xfapply(fwarp(1, 1, ind1), xcenIn, centerSlice, xp, yp, xp, yp)
!!$  xp2 = xx * cWarpBeta(ind2) + yy * sWarpAlpha(ind2) * sWarpBeta(ind2) + &
!!$      zz * (cWarpAlpha(ind2) * sWarpBeta(ind2) + warpXZfac(ind2)) + xcenIn + axisXoffset
!!$  yp2 = yy * cWarpAlpha(ind2) - zz * (sWarpAlpha(ind2) - warpYZfac(ind2)) + centerSlice
!!$  call xfapply(fwarp(1, 1, ind2), xcenIn, centerSlice, xp2, yp2, xp2, yp2)
!!$  xp3 = xx * cWarpBeta(ind3) + yy * sWarpAlpha(ind3) * sWarpBeta(ind3) + &
!!$      zz * (cWarpAlpha(ind3) * sWarpBeta(ind3) + warpXZfac(ind3)) + xcenIn + axisXoffset
!!$  yp3 = yy * cWarpAlpha(ind3) - zz * (sWarpAlpha(ind3) - warpYZfac(ind3)) + centerSlice
!!$  call xfapply(fwarp(1, 1, ind3), xcenIn, centerSlice, xp3, yp3, xp3, yp3)
!!$  xp4 = xx * cWarpBeta(ind4) + yy * sWarpAlpha(ind4) * sWarpBeta(ind4) + &
!!$      zz * (cWarpAlpha(ind4) * sWarpBeta(ind4) + warpXZfac(ind4)) + xcenIn + axisXoffset
!!$  yp4 = yy * cWarpAlpha(ind4) - zz * (sWarpAlpha(ind4) - warpYZfac(ind4)) + centerSlice
!!$  call xfapply(fwarp(1, 1, ind4), xcenIn, centerSlice, xp4, yp4, xp4, yp4)
!!$  xp = f1 * xp + f2 * xp2 + f3 * xp3 + f4 * xp4
!!$  yp = f1 * yp + f2 * yp2 + f3 * yp3 + f4 * yp4


! Finds the point at centered Z coordinate zz projecting to
! xproj, yproj in view iv of original projections.  xx is X index in
! reconstruction, yy is slice number in original projections
!
subroutine findProjectingPoint(xproj, yproj, zz, iv, xx, yy)
  use tiltvars
  implicit none
  real*4 xproj, yproj, zz, xx, yy
  integer*4 iv, iter, ifDone, ixAssay, iyAssay
  real*4 xprojFix11, xprojZ11, yprojFix11, yprojZ11, xprojFix21, xprojZ21, &
      yprojFix21, yprojZ21, xprojFix12, xprojZ12, yprojFix12, yprojZ12
  real*4 xp11, yp11, xp12, yp12, xp21, yp21, xerr, yerr, dxpx, dxpy, dypx
  real*4 dypy, fx, fy, den
  !
  iter = 1
  ifDone = 0
  do while (ifDone == 0 .and. iter <= 5)
    ixAssay = floor(xx)
    iyAssay = floor(yy)
    call localProjFactors(ixAssay, iyAssay, iv, xprojFix11, xprojZ11, yprojFix11, &
        yprojZ11)
    call localProjFactors(ixAssay + 1, iyAssay, iv, xprojFix21, xprojZ21, yprojFix21, &
        yprojZ21)
    call localProjFactors(ixAssay, iyAssay + 1, iv, xprojFix12, xprojZ12, yprojFix12, &
        yprojZ12)
    xp11 = xprojFix11 + xprojZ11 * zz
    yp11 = yprojFix11 + yprojZ11 * zz
    xp21 = xprojFix21 + xprojZ21 * zz
    yp21 = yprojFix21 + yprojZ21 * zz
    xp12 = xprojFix12 + xprojZ12 * zz
    yp12 = yprojFix12 + yprojZ12 * zz
    xerr = xproj - xp11
    yerr = yproj - yp11
    dxpx = xp21 - xp11
    dxpy = xp12 - xp11
    dypx = yp21 - yp11
    dypy = yp12 - yp11
    den = dxpx * dypy - dxpy * dypx
    fx = (xerr * dypy - yerr * dxpy) / den
    fy = (dxpx * yerr - dypx * xerr) / den
    xx = ixAssay + fx
    yy = iyAssay + fy
    if (fx > -0.1 .and. fx < 1.1 .and. fy > -0.1 .and. fy < 1.1) &
        ifDone = 1
    iter = iter + 1
  enddo
  return
end subroutine findProjectingPoint

!
! Compute space needed for cosine stretched data
!
subroutine set_cos_stretch()
  use tiltvars
  implicit none
  integer*4 lsliceMin, lsliceMax, iv, ix, iy, lslice
  real*4 tanAlph, xpMax, xpMin, zz, zPart, yy, xproj
  ! make the indexes be bases, numbered from 0
  !
  indStretchLine(1) = 0
  lsliceMin = min(isliceEnd, isliceStart)
  lsliceMax = max(isliceEnd, isliceStart)
  if (ifAlpha < 0) then
    !
    ! New-style X tilting: SET MINIMUM NUMBER OF INPUT SLICES HERE
    !
    lsliceMin = centerSlice + (lsliceMin - centerSlice) * cosAlpha(1) +  &
        yOffset * sinAlpha(1) - 0.5 * ithickOut * abs(sinAlpha(1)) - 1.
    lsliceMax = centerSlice + (lsliceMax - centerSlice) * cosAlpha(1) +  &
        yOffset * sinAlpha(1) + 0.5 * ithickOut * abs(sinAlpha(1)) + 2.
    tanAlph = sinAlpha(1) / cosAlpha(1)
    lsliceMin = max(1, lsliceMin)
    lsliceMax = min(lsliceMax, nyProj)
  endif
  do iv = 1, numViews
    xpMax = 1
    xpMin = nxProj
    !
    ! find min and max position of 8 corners of reconstruction
    !
    do ix = 1, iwidth, iwidth - 1
      do iy = 1, ithickBP, ithickBP - 1
        do lslice = lsliceMin, lsliceMax, max(1, lsliceMax - lsliceMin)
          zz = (iy - ycenOut) * compress(iv)
          if (ifAlpha < 0) zz = compress(iv) * &
              (iy - (ycenOut - nint(tanAlph * (lslice - centerSlice))))
          if (ifAlpha <= 0) then
            zPart = zz * sinBeta(iv) + xcenIn + axisXoffset
          else
            yy = lslice - centerSlice
            zPart = yy * sinAlpha(iv) * sinBeta(iv) + zz * (cosAlpha(iv) * sinBeta(iv) + &
                xzfac(iv)) + xcenIn + axisXoffset
          endif
          xproj = zPart + (ix - xcenOut) * cosBeta(iv)
          xpMin = max(1., min(xpMin, xproj))
          xpMax = min(float(nxProj), max(xpMax, xproj))
        enddo
      enddo
    enddo
    ! print *,iv, xpmin, xpmax
    !
    ! set up extent and offset of stretches
    !
    stretchOffset(iv) = xpMin / cosBeta(iv) - 1. / interpFacStretch
    nxStretched(iv) = interpFacStretch * (xpMax - xpMin) / cosBeta(iv) + 2.
    indStretchLine(iv + 1) = indStretchLine(iv) + nxStretched(iv)
    ! print *,iv, xpmin, xpmax, stretchOffset(iv), nxStretched(iv), indStretchLine(iv)
  enddo
  return
end subroutine set_cos_stretch


! Determine starting and ending input slice needed to reconstruct
! each output slice, as well as the maximum needed over all slices
! for a series of numbers of output slices up to numEval
!
subroutine setNeededSlices(maxNeeds, numEval)
  use tiltvars
  implicit none
  integer*4 numEval, maxNeeds(*)
  integer*4 lsliceMin, lsliceMax, ierr, itry, nxAssay, minSlice, ixAssay
  integer*4 maxSlice, iassay, ixSample, iv, iy, iyp
  real*4 dxAssay, dxTemp, xx, yy, zz, xp, yp
  real*4 xprojFix, xprojZ, yprojFix, yprojZ, xproj, yproj
  lsliceMin = isliceStart
  lsliceMax = isliceEnd
  if (ifAlpha < 0) then
    lsliceMin = centerSlice + (isliceStart - centerSlice) * cosAlpha(1) +  &
        yOffset * sinAlpha(1) - 0.5 * ithickOut * abs(sinAlpha(1)) - 1.
    lsliceMax = centerSlice + (isliceEnd - centerSlice) * cosAlpha(1) +  &
        yOffset * sinAlpha(1) + 0.5 * ithickOut * abs(sinAlpha(1)) + 2.
    lsliceMin = max(1, lsliceMin)
    lsliceMax = min(lsliceMax, nyProj)
  endif
  indNeededBase = lsliceMin - 1
  numNeedSE = lsliceMax - indNeededBase
  allocate(neededStarts(numNeedSE), neededEnds(numNeedSE), stat = ierr)
  call memoryErrorUC(ierr, 'ARRAYS needStarts/needEnds')

  do itry = lsliceMin, lsliceMax

    if (ifAlpha <= 0 .and. nxWarp == 0) then
      !
      ! regular case is simple: just need the current slice
      !
      neededStarts(itry - indNeededBase) = itry
      neededEnds(itry - indNeededBase) = itry
    else
      !
      ! for old-style X-tilt or local alignment, determine what
      ! slices are needed by sampling
      ! set up sample points: left and right if no warp,
      ! or half the warp spacing
      !
      if (nxWarp == 0) then
        nxAssay = 2
        dxAssay = iwidth - 1
      else
        dxTemp = idelXwarp / 2
        nxAssay = max(2., iwidth / dxTemp + 1.)
        dxAssay = (iwidth - 1.) / (nxAssay - 1.)
      endif
      !
      ! sample top and bottom at each position
      !
      minSlice = nyProj + 1
      maxSlice = 0
      do iassay = 1, nxAssay
        ixAssay = nint(1 + (iassay - 1) * dxAssay)
        do iv = 1, numViews
          if (.not. recReproj) then
            ixSample = nint(ixAssay - xcenOut + xcenIn + axisXoffset)
            if (nxWarp .ne. 0) then
              call localProjFactors(ixAssay, itry, iv, xprojFix, &
                  xprojZ, yprojFix, yprojZ)
            endif
            do iy = 1, ithickBP, ithickBP - 1
              !
              ! for each position, find back-projection location
              ! transform if necessary, and use to get min and
              ! max slices needed to get this position
              !
              xx = ixSample - xcenOut
              yy = itry - centerSlice
              zz = iy - ycenOut
              xp = xx * cosBeta(iv) + yy * sinAlpha(iv) * sinBeta(iv) + &
                  zz * (cosAlpha(iv) * sinBeta(iv) + xzfac(iv)) + xcenIn + axisXoffset
              yp = yy * cosAlpha(iv) - zz * (sinAlpha(iv) - yzfac(iv)) + centerSlice
              if (nxWarp .ne. 0) then
                xp = xprojFix + xprojZ * zz
                yp = yprojFix + yprojZ * zz
              endif
              iyp = max(1., yp)
              minSlice = min(minSlice, iyp)
              maxSlice = max(maxSlice, min(nyProj, iyp + 1))
              ! if (debug) print *,xx, yy, zz, iyp, minslice, maxslice
            enddo
          else
            !
            ! Projections: get Y coordinate in original projection
            ! if local, get the X coordinate in reconstruction too
            ! then get the refinement
            xproj = ixAssay + xprojOffset
            yproj = itry + yprojOffset
            do iy = 1, ithickReproj, ithickReproj - 1
              zz = iy + minYreproj - 1 - ycenOut
              yy = (yproj + zz * (sinAlpha(iv) - yzfac(iv)) - centerSlice) /  &
                  cosAlpha(iv) + centerSlice
              if (nxWarp .ne. 0) then
                xx = (xproj - yy * sinAlpha(iv) * sinBeta(iv) - zz * (cosAlpha(iv) * &
                    sinBeta(iv) + xzfac(iv)) - xcenIn - axisXoffset) / cosBeta(iv) + &
                    xcenOut
                call findProjectingPoint(xproj, yproj, zz, iv, xx, yy)
              endif
              iyp = max(1., yy - yprojOffset)
              minSlice = min(minSlice, iyp)
              maxSlice = max(maxSlice, min(nyProj, iyp + 1))
            enddo
          endif
        enddo
      enddo
      !
      ! set up starts and ends
      !
      neededStarts(itry - indNeededBase) = max(1, minSlice)
      neededEnds(itry - indNeededBase) = min(nyProj, maxSlice)
    endif
  enddo
  !
  ! Count maximum # of slices needed for number of slices to be computed
  do iv = 1, numEval
    maxNeeds(iv) = 0
    do iy = lsliceMin, lsliceMax + 1 - iv
      maxSlice = neededEnds(iy + iv - 1 - indNeededBase) + 1 -  &
          neededStarts(iy - indNeededBase)
      maxNeeds(iv) = max(maxNeeds(iv), maxSlice)
    enddo
  enddo
  return
end subroutine setNeededSlices


! Allocate array, trying to get enough to do numEval slices without
! reloading any data, based on the numbers in maxNeeds, and trying fewer
! slices down to minLoad if that fails.  MinMemory is the minimum amount
! it will allocate.
!
integer*4 function allocateArray(maxNeeds, numEval, minLoad, minMemory)
  use tiltvars
  implicit none
  integer*4 maxNeeds(*), numEval, minLoad, minMemory, ierr, i
  integer(kind = 8) memNeed, minNeed
  minNeed = inPlaneSize
  minNeed = indLoadBase + ipExtraSize + minNeed * &
      (neededEnds(numNeedSE) + 1 - neededStarts(1)) + 32
  if (minNeed > minMemory) minNeed = minMemory
  do i = numEval, minLoad, -1
    memNeed = inPlaneSize
    memNeed = indLoadBase + ipExtraSize + memNeed * maxNeeds(i) + 32
    memNeed = max(memNeed, minNeed)
    if (memNeed < 2147000000) then
      allocate(array(memNeed), stat = ierr)
      if (ierr == 0) then
        maxStack = memNeed
        numPlanes = (maxStack - indLoadBase - ipExtraSize + 1) / inPlaneSize
        allocateArray = i
        write(*,'(/,a,i5,a)') 'Allocated', nint(maxStack / (1024 * 256.)), &
            ' MB for stack array'
        return
      endif
    endif
  enddo
  allocateArray = 0
  return
end function allocateArray


subroutine reproject(array, nxs, nys, nxOut, sinAngle, cosAngle, xRayStart, &
    yRayStart, numPixInRay, maxRayPixels, fill, projLine, linear, noScale)
  implicit none
  integer*4 nxs, nys, nxOut, numPixInRay(*), maxRayPixels, linear, noScale
  real*4 array(nxs, nys), xRayStart(*), yRayStart(*), fill, projLine(*)
  integer*4 ixOut, iray, ixr, iyr, numRayPts, idir
  real*4 sinAngle, cosAngle, rayFac, rayAdd, xRay, yRay, pixTemp, fullFill
  real * 4 dx, dy, v2, v4, v6, v8, v5, a, b, c, d
  !
  rayFac = 1. / maxRayPixels
  fullFill = fill
  if (noScale .ne. 0) then
    rayFac = 1.
    fullFill = fill * maxRayPixels
  endif
  do ixOut = 1, nxOut
    projLine(ixOut) = fullFill
    numRayPts = numPixInRay(ixOut)
    if (numRayPts > 0) then
      pixTemp = 0.
      if (sinAngle .ne. 0.) then
        if (linear == 0) then
          do iray = 0, numRayPts - 1
            xRay = xRayStart(ixOut) + iray * sinAngle
            yRay = yRayStart(ixOut) + iray * cosAngle
            ixr = nint(xRay)
            iyr = nint(yRay)
            dx = xRay - ixr
            dy = yRay - iyr
            v2 = array(ixr, iyr - 1)
            v4 = array(ixr - 1, iyr)
            v5 = array(ixr, iyr)
            v6 = array(ixr + 1, iyr)
            v8 = array(ixr, iyr + 1)
            !
            a = (v6 + v4) * .5 - v5
            b = (v8 + v2) * .5 - v5
            c = (v6 - v4) * .5
            d = (v8 - v2) * .5
            pixTemp = pixTemp + a * dx * dx + b * dy * dy + c * dx + d * dy + v5
          enddo
        else
          do iray = 0, numRayPts - 1
            xRay = xRayStart(ixOut) + iray * sinAngle
            yRay = yRayStart(ixOut) + iray * cosAngle
            ixr = xRay
            iyr = yRay
            dx = xRay - ixr
            dy = yRay - iyr
            pixTemp = pixTemp + (1 - dy) * &
                ((1. -dx) * array(ixr, iyr) + dx * array(ixr + 1, iyr)) + dy * &
                ((1. -dx) * array(ixr, iyr + 1) + dx * array(ixr + 1, iyr + 1))
          enddo
        endif
      else
        !
        ! vertical projection
        !
        ixr = nint(xRayStart(ixOut))
        iyr = nint(yRayStart(ixOut))
        idir = sign(1., cosAngle)
        do iray = 0, numRayPts - 1
          pixTemp = pixTemp + array(ixr, iyr + idir * iray)
        enddo
      endif

      rayAdd = rayFac * (maxRayPixels - numRayPts) * fill
      projLine(ixOut) = rayFac * pixTemp + rayAdd
    endif
  enddo
  return
end subroutine reproject

! Reprojects slices from lsliceStart to lsliceEnd and writes the reprojections
! Fewer slices may be done if on the GPU, and the ending slice done is
! returned in lsliceEnd.  INLOADSTART and INLOADEND are slice numbers of
! started and ending slices loaded; min/max/sum densities are maintained
! in DMIN, DMAX, and DTOT8
!
subroutine reprojectRec(lsliceStart, lsliceEnd, inLoadStart, inLoadEnd, dmin, dmax, dtot8)
  use tiltvars
  implicit none
  integer*4 lsliceStart, lsliceEnd, inLoadStart, inLoadEnd
  real*4 dmin, dmax
  integer*4 iv, ix, iy, iz, ixp, line, i, iys, ind
  integer*4 ind1, ind2, ind3, ind4, load
  real*4 cosAlph, sinAlph, cosBet, sinBet, delZ, delX, fz, oneMfz, zz, xx, fx
  real*4 oneMfx, yy, fy, oneMfy, xproj, yproj, d11, d12, d21, d22
  real*4 f1, f2, f3, f4, xxGood, yyGood, zzGood
  real*4 yEndTol, xprojMin, xprojMax, xJump, zJump, delY, diffXmax, diffYmax
  integer*4 indBase, nxLoad, ixc, lastZdone
  integer*4 indJump, numJump, lgpuEnd, lineBase
  real*4 ycenAdj
  real*8 sum, dtot8, wallTime, wallStart, wallCumul
  logical*4 tryJump
  real*4 reprojDelZ
  integer*4 gpuReproject, gpuReprojLocal

  yEndTol = 3.05
  xJump = 5.0
  nxLoad = maxXload + 1 - minXload
  wallStart = wallTime()
  wallCumul = 0.
  ycenAdj = ycenOut - (minYreproj - 1)
  !
  if (useGPU .and. loadGpuStart > 0) then
    indJump = 1
    !
    ! GPU REPROJECTION: Find last slice that can be done
    do lgpuEnd = lsliceEnd, lsliceStart, -1
      if (neededEnds(lgpuEnd - indNeededBase) <= loadGpuEnd) then
        indJump = 0
        exit
      endif
    enddo
    if (indJump == 0) then
      !
      ! Loop on views; do non-local case first
      do iv = 1, numViews
        if (nxWarp == 0) then
          delZ = reprojDelZ(sinBeta(iv), cosBeta(iv), sinAlpha(iv), cosAlpha(iv), &
              xzfac(iv), yzfac(iv))
          indJump = gpuReproject(reprojLines, sinBeta(iv), cosBeta(iv), sinAlpha(iv), &
              cosAlpha(iv), xzfac(iv), yzfac(iv), delZ, lsliceStart, lgpuEnd, &
              ithickReproj, xcenOut, xcenIn + axisXoffset, minXreproj, xprojOffset, &
              ycenOut, minYreproj, yprojOffset, centerSlice, ifAlpha, dmeanIn)
        else
          !
          ! GPU with local alignments: fill warpDelz array for all lines
          do line = lsliceStart, lgpuEnd
            call fillWarpDelz(warpDelz(1 + (line - lsliceStart) * numWarpDelz))
          enddo
          !
          ! Get the xprojmin and max adjusted by 5
          xprojMin = 10000000.
          xprojMax = 0.
          do load = inLoadStart, inLoadEnd
            iys = nint(load + yprojOffset)
            do ix = 1, nxLoad, nxLoad - 1
              call localProjFactors(ix + minXload - 1, iys, iv, sinBet, &
                  cosBet, sinAlph, cosAlph)
              sinAlph = sinBet + (1 - ycenAdj) * cosBet
              cosAlph = sinBet + (ithickReproj - ycenAdj) * cosBet
              xprojMin = min(xprojMin, sinAlph - 5., cosAlph - 5.)
              xprojMax = max(xprojMax, sinAlph + 5., cosAlph + 5.)
            enddo
          enddo
          ! print *,'xprojmin, max', xprojMin, xprojMax
          !
          ! Do it
          indJump = gpuReprojLocal(reprojLines, sinBeta(iv), cosBeta(iv), sinAlpha(iv), &
              cosAlpha(iv), xzfac(iv), yzfac(iv), nxWarp, nyWarp, ixStartWarp, &
              iyStartWarp, idelXwarp, idelYwarp, warpDelz, numWarpDelz, &
              dxWarpDelz, xprojMin, xprojMax, lsliceStart, lgpuEnd, &
              ithickReproj, iv, xcenOut, xcenIn, axisXoffset, minXload, xprojOffset, &
              ycenAdj, yprojOffset, centerSlice, dmeanIn)
        endif
        if (indJump .ne. 0) exit
        wallCumul = wallCumul + wallTime() - wallStart
        call writeReprojLines(iv, lsliceStart, lgpuEnd, dmin, dmax, dtot8)
        wallStart = wallTime()
      enddo
      if (indJump == 0) then
        if (debug) write(*, '(a,f9.5)') 'GPU reprojection time', wallCumul
        lsliceEnd = lgpuEnd
        return
      endif
    endif
  endif
  !
  ! CPU REPROJECTION: loop on views; first handle non-local alignments
  do iv = 1, numViews
    if (threshPolarity .ne. 0.) then
      call thresholdedReproj()
    elseif (nxWarp == 0) then
      !
      ! Get the delta z for this view
      cosAlph = cosAlpha(iv)
      sinAlph = sinAlpha(iv)
      cosBet = cosBeta(iv)
      sinBet = sinBeta(iv)
      delZ = reprojDelZ(sinBet, cosBet, sinAlph, cosAlph, xzfac(iv), yzfac(iv))
      ! print *,sbeta, cbeta, salf, calf, xzfac(iv), yzfac(iv)
      ! print *,delx, delz
      !
      ! Loop on the output lines to be done
      do line = lsliceStart, lsliceEnd
        lineBase = (line - lsliceStart) * iwidth + 1
        call reprojOneAngle(array(indLoadBase), reprojLines(lineBase), &
            inLoadStart, inLoadEnd, line, cosBet, sinBet, cosAlph, sinAlph, &
            delZ, iwidth, ithickReproj, inPlaneSize, nxLoad, minXreproj, &
            minYreproj, xprojOffset, yprojOffset, xcenOut, ycenOut, &
            xcenIn + axisXoffset, centerSlice, ifAlpha, xzfac(iv), yzfac(iv), dmeanIn)
      enddo
    else
      !
      ! LOCAL ALIGNMENTS
      !
      ! first step: precompute all the x/yprojf/z  for all slices
      ! general BUG ycenAdj replaces ycenOut - minYreproj (off by 1)
      xprojMin = 10000000.
      xprojMax = 0.
      do load = inLoadStart, inLoadEnd
        indBase = nxLoad * (load - inLoadStart)
        iys = nint(load + yprojOffset)
        do ix = 1, nxLoad
          ind = indBase + ix
          call localProjFactors(ix + minXload - 1, iys, iv, xprojfs(ind), &
              xprojzs(ind), yprojfs(ind), yprojzs(ind))
          if (ix == 1) xprojMin = min(xprojMin, &
              xprojfs(ind) + (1 - ycenAdj) * xprojzs(ind),  xprojfs(ind) &
              + (ithickReproj - ycenAdj) * xprojzs(ind))
          if (ix == nxLoad) xprojMax = max(xprojMax, &
              xprojfs(ind) + (1 - ycenAdj) * xprojzs(ind),  xprojfs(ind) &
              + (ithickReproj - ycenAdj) * xprojzs(ind))
        enddo
      enddo
      ! print *,'xprojmin, max', xprojMin, xprojMax
      !
      ! loop on lines to be done
      do line = lsliceStart, lsliceEnd
        lineBase = (line - lsliceStart) * iwidth
        call fillWarpDelz(warpDelz)
        ! print *,iv, line, inloadstr, inloadend
        !
        ! loop on pixels across line
        yproj = line + yprojOffset
        do ixp = 1, iwidth
          !
          ! Get x projection coord, starting centered Z coordinate, and
          ! approximate x and y coordinates
          ! xproj, yproj are coordinates in original projections
          ! Equations relate them to coordinates in reconstruction
          ! and then X coordinate is adjusted to be a loaded X index
          ! and Y coordinate is adjusted to be a slice of reconstruction
          xproj = ixp + xprojOffset
          zz = 1. -ycenAdj
          sum = 0.
          ! print *,ixp, xproj, yproj, xx, yy
          ! BUG these lines needed to be swapped and yprojOffset deferred
          yy = (yproj + zz * (sinAlpha(iv) - yzfac(iv)) - centerSlice) / cosAlpha(iv) + &
              centerSlice
          xx = (xproj - yy * sinAlpha(iv) * sinBeta(iv) -  &
              zz * (cosAlpha(iv) * sinBeta(iv) + &
              xzfac(iv)) - xcenIn - axisXoffset) / cosBeta(iv) + xcenOut - (minXload - 1)
          yy = yy - yprojOffset
          !
          ! Move on ray up in Z
          lastZdone = 0
          tryJump = .true.
          diffXmax = 0
          diffYmax = 0
          ! BUG ?? Surely this should be abs(sinBeta(iv))
          zJump = xJump * cosBeta(iv) / max(0.2, abs(sinBeta(iv)))
          do while (zz < ithickReproj + 1 - ycenAdj .and. &
              lastZdone == 0)
            if (xproj < xprojMin - 5. .or. xproj > xprojMax + 5.) then
              sum = sum + dmeanIn
            else
              call loadedProjectingPoint(xproj, yproj, zz, nxLoad, &
                  inLoadStart, inLoadEnd, xx, yy)
              !
              ! If X or Y is out of bounds, fill with mean
              if (yy < inLoadStart - yEndTol .or. yy > inLoadEnd + yEndTol &
                  .or. xx < 1. .or. xx >= nxLoad) then
                sum = sum + dmeanIn
              else
                !
                ! otherwise, get x, y, z indexes, clamp y to limits, allow
                ! a fractional Z pixel at top of volume
                ix = xx
                fx = xx - ix
                oneMfx = 1. - fx
                yy = max(float(inLoadStart), min(inLoadEnd - 0.01, yy))
                iy = yy
                fy = yy - iy
                oneMfy = 1. - fy
                ! BUG ????  Shouldn't this be + ycenOut - minYreproj?
                iz = max(1., zz + ycenAdj)
                fz = zz + ycenAdj - iz
                oneMfz = 1. - fz
                if (iz == ithickReproj) then
                  iz = iz - 1
                  fz = oneMfz
                  oneMfz = 0.
                  lastZdone = 1
                endif
                !
                ! Do the interpolation
                d11 = oneMfx * oneMfy
                d12 = oneMfx * fy
                d21 = fx * oneMfy
                d22 = fx * fy
                ind = indLoadBase + inPlaneSize * (iy - inLoadStart) + (iz - 1) * nxLoad &
                    + ix - 1
                sum = sum + oneMfz * (d11 * array(ind) &
                    + d12 * array(ind + inPlaneSize) + d21 * array(ind + 1) &
                    + d22 * array(ind + inPlaneSize + 1)) &
                    + fz * (d11 * array(ind + nxLoad) &
                    + d12 * array(ind + inPlaneSize + nxLoad) &
                    + d21 * array(ind + 1 + nxLoad) &
                    + d22 * array(ind + inPlaneSize + 1 + nxLoad))
                !
                do while(tryJump)
                  !
                  ! If jumping is OK, save the current position and compute
                  ! how many steps can be jumped, stopping below the top
                  xxGood = xx
                  yyGood = yy
                  zzGood = zz
                  ind = max(1., min(float(numWarpDelz), xx / dxWarpDelz))
                  delZ = warpDelz(ind)
                  numJump = zJump / delZ
                  if (zz + zJump > ithickReproj - ycenAdj - 1) then
                    numJump = (ithickReproj - ycenAdj - 1 - zz) / delZ
                    tryJump = .false.
                  endif
                  if (numJump > 0) then
                    !
                    ! Make the jump, find the projecting point;
                    ! if it's out of bounds restore last point
                    zz = zz + numJump * delZ
                    xx = xx + numJump * sinBeta(iv)
                    call loadedProjectingPoint(xproj, yproj, zz, &
                        nxLoad, inLoadStart, inLoadEnd, xx, yy)
                    if (yy < inLoadStart .or. yy > inLoadEnd .or. &
                        xx < 1. .or. xx >= nxLoad) then
                      numJump = 0
                      xx = xxGood
                      yy = yyGood
                      zz = zzGood
                      tryJump = .false.
                    else
                      delX = (xx - xxGood) / numJump
                      delY = (yy - yyGood) / numJump
                    endif
                  endif
                  !
                  ! Loop on points from last one to final one
                  do indJump = 1, numJump
                    xx = xxGood + indJump * delX
                    yy = yyGood + indJump * delY
                    zz = zzGood + indJump * delZ
                    ix = xx
                    fx = xx - ix
                    oneMfx = 1. - fx
                    iy = min(int(yy), inLoadEnd - 1)
                    fy = yy - iy
                    oneMfy = 1. - fy
                    ! BUG again, need ycenAdj
                    iz = zz + ycenAdj
                    fz = zz + ycenAdj - iz
                    oneMfz = 1. - fz
                    d11 = oneMfx * oneMfy
                    d12 = oneMfx * fy
                    d21 = fx * oneMfy
                    d22 = fx * fy
                    ind = indLoadBase + inPlaneSize * (iy - inLoadStart) + (iz - 1) * &
                        nxLoad + ix - 1
                    sum = sum + oneMfz * (d11 * array(ind) &
                        + d12 * array(ind + inPlaneSize) + d21 * array(ind + 1) &
                        + d22 * array(ind + inPlaneSize + 1)) &
                        + fz * (d11 * array(ind + nxLoad) &
                        + d12 * array(ind + inPlaneSize + nxLoad) &
                        + d21 * array(ind + 1 + nxLoad) &
                        + d22 * array(ind + inPlaneSize + 1 + nxLoad))
                    ! fx = xx
                    ! fy = yy
                    ! call loadedProjectingPoint(xproj, yproj, zz, &
                    ! nxload, inloadstr, inloadend, fx, fy)
                    ! diffxmax = max(diffxmax , abs(fx - xx))
                    ! diffymax = max(diffymax , abs(fy - yy))
                  enddo
                enddo
              endif
            endif
            !
            ! Adjust Z by local factor, move X approximately for next pixel
            ind = max(1., min(float(numWarpDelz), xx / dxWarpDelz))
            zz = zz + warpDelz(ind)
            xx = xx + sinBeta(iv)
          enddo
          reprojLines(lineBase + ixp) = sum
          ! write (*,'(i5,2f10.4)') ixp, diffxmax, diffymax
        enddo
      enddo
    endif
    wallCumul = wallCumul + wallTime() - wallStart
    call writeReprojLines(iv, lsliceStart, lsliceEnd, dmin, dmax, dtot8)
    wallStart = wallTime()
  enddo
  if (debug) write(*, '(a,f9.5)') 'CPU reprojection time', wallCumul

CONTAINS
  !
  ! compute delta z as function of X across the loaded slice
  ! which is not ideal since the data will not be coming from slice
  subroutine fillWarpDelz(warpDelZ)
    real*4 warpDelZ(*)
    iys = nint(line + yprojOffset)
    do i = 1, numWarpDelz
      xx = 1 + dxWarpDelz * (i - 1)
      ixc = nint(xx + minXload - 1 - xcenOut + xcenIn + axisXoffset)
      call local_factors(ixc, iys, iv, ind1, ind2, ind3, ind4, f1, f2, &
          f3, f4)
      warpDelZ(i) = f1 * reprojDelZ(sWarpBeta(ind1), cWarpBeta(ind1), &
          sWarpAlpha(ind1), cWarpAlpha(ind1), warpXZfac(ind1), warpYZfac(ind1)) &
          + f2 * reprojDelZ(sWarpBeta(ind2), cWarpBeta(ind2), &
          sWarpAlpha(ind2), cWarpAlpha(ind2), warpXZfac(ind2), warpYZfac(ind2)) &
          + f3 * reprojDelZ(sWarpBeta(ind3), cWarpBeta(ind3), &
          sWarpAlpha(ind3), cWarpAlpha(ind3), warpXZfac(ind3), warpYZfac(ind3)) &
          + f4 * reprojDelZ(sWarpBeta(ind4), cWarpBeta(ind4), &
          sWarpAlpha(ind4), cWarpAlpha(ind4), warpXZfac(ind4), warpYZfac(ind4))
    enddo
    ! print *,'got delz eg:', wrpdlz(1), wrpdlz(numWarpDelz/2), &
    ! wrpdlz(numWarpDelz)
  end subroutine fillWarpDelz


  ! Does a reprojection only of discrete points beyond a threshold
  !
  subroutine thresholdedReproj()
    real*4 polarity, f11, f12, f21, f22, rlX, rlSlice, rlZ
    integer*4 numLines, iyp, loadedSlice
    numLines = lsliceEnd + 1 - lsliceStart
    reprojLines(1 : numLines * iwidth) = threshFillVal * ithickReproj
    polarity = sign(1., threshPolarity)
    do loadedSlice = inLoadStart, inLoadEnd
      ind = indLoadBase + (loadedSlice - inLoadStart) * inPlaneSize
      do iy = 1, ithickReproj
        do ix = 1, nxLoad
          if (polarity * (array(ind) - threshForReproj) >= 0) then
            !
            ! Get real coordinate position within full reconstruction and projection
            ! position in full aligned stack
            rlX = ix + minXreproj - 1
            rlSlice = loadedSlice
            rlZ = iy + minYreproj - 1
            call projectionPosition(iv, rlX, rlZ, rlSlice, xproj, yproj, ind1, ind2, &
                f11, f12, f21, f22)
            !
            ! Adjust for position in reprojection being produced then adjust Y to 
            ! be an index in the lines being produced
            xproj = xproj - xprojOffset
            yproj = yproj - (minZreproj - 1 + yprojOffset)
            yproj = yproj + isliceStart - lsliceStart
            ixp = xproj
            iyp = yproj
            if (ixp >= 0 .and. ixp <= iwidth .and. iyp >= 0 .and. iyp <= numLines) then
              !
              ! If any of the 4 actual pixels is in range, get the interpolation factors
              ! for those 4 surrounding pixels
              fx = xproj - ixp
              fy = yproj - iyp
              f11 = (1. - fx) * (1. - fy)
              f12 = (1. - fx) * fy
              f21 = fx * (1. - fy)
              f22 = fx * fy
              ind1 = ixp + iwidth * (iyp - 1)
              !
              ! If NOT summing, just mark any pixel with fraction above threshold
              ! Otherwise add the fraction times the value
              if (threshSumFac < 1.) then
                if (f11 >= threshSumFac .and. ixp > 0 .and. iyp > 0) &
                    reprojLines(ind1) = threshMarkVal
                if (f12 >= threshSumFac .and. ixp > 0 .and. iyp < numLines) &
                    reprojLines(ind1 + iwidth) = threshMarkVal
                if (f21 >= threshSumFac .and. ixp < iwidth .and. iyp > 0) &
                    reprojLines(ind1 + 1) = threshMarkVal
                if (f22 >= threshSumFac .and. ixp < iwidth .and. iyp < numLines) &
                    reprojLines(ind1 + 1 + iwidth) = threshMarkVal
              else
                if (ixp > 0 .and. iyp > 0) reprojLines(ind1) = &
                    reprojLines(ind1) + f11 * threshMarkVal - threshFillVal
                if (ixp > 0 .and. iyp < numLines) reprojLines(ind1 + iwidth) = &
                    reprojLines(ind1 + iwidth) + f12 * threshMarkVal - threshFillVal
                if (ixp < iwidth .and. iyp > 0) reprojLines(ind1 + 1) = &
                    reprojLines(ind1 + 1) + f21 * threshMarkVal - threshFillVal
                if (ixp < iwidth .and. iyp < numLines) reprojLines(ind1 + 1 + iwidth) = &
                    reprojLines(ind1 + 1 + iwidth) + f22 * threshMarkVal - threshFillVal
              endif
            endif
          endif
          ind = ind + 1
        enddo
      enddo
    enddo
  end subroutine thresholdedReproj

end subroutine reprojectRec

! Finds loaded point that projects to xproj, yproj at centered Z value
! zz, using stored values for [xy]zfac[fv].
! Takes starting value in xx, yy and returns found value.
! X coordinate needs to be a loaded X index
! Y coordinate yy is in slices of reconstruction, yproj in original proj
!
subroutine loadedProjectingPoint(xproj, yproj, zz, nxLoad, inLoadStart, inLoadEnd, xx, yy)
  use tiltvars
  implicit none
  real*4 xproj, yproj, zz, xx, yy
  integer*4 nxLoad, inLoadStart, inLoadEnd
  integer*4 iter, ifDone, ind, ix, iy, ifOut, i
  real*4 xp11, yp11, xp12, yp12, xp21, yp21, xErr, yErr, dypx, dxpy, dxpx
  real*4 dypy, den, fx, fy
  ! logical*4 dbout
  ! dbout = abs(xproj - 700.) < 3 .and. abs(zz) < 3
  ! dbout = .false.
  ! print *,'Finding proj pt to', xproj, yproj, zz
  iter = 0
  ifDone = 0
  do while (ifDone == 0 .and. iter < 5)
    ix = floor(xx)
    iy = floor(yy)
    ifOut = 0
    if (ix < 1 .or. ix >= nxLoad .or. iy < inLoadStart .or. &
        iy >= inLoadEnd) then
      ifOut = 1
      ix = min(nxLoad - 1, max(1, ix))
      iy = min(inLoadEnd - 1, max(inLoadStart, iy))
    endif
    ind = nxLoad * (iy - inLoadStart) + ix
    xp11 = xprojfs(ind) + xprojzs(ind) * zz
    yp11 = yprojfs(ind) + yprojzs(ind) * zz
    xp21 = xprojfs(ind + 1) + xprojzs(ind + 1) * zz
    yp21 = yprojfs(ind + 1) + yprojzs(ind + 1) * zz
    xp12 = xprojfs(ind + nxLoad) + xprojzs(ind + nxLoad) * zz
    yp12 = yprojfs(ind + nxLoad) + yprojzs(ind + nxLoad) * zz
    ! write(*,101) 'facs', (xprojfs(i), xprojzs(i), yprojfs(i), &
    ! yprojzs(i), i=ind, ind+1)
    ! write(*,101) 'xps', xx, yy, zz, xp11, yp11, xp21, yp21, xp12, yp12
    !101     format(a,9f8.2)
    xErr = xproj - xp11
    yErr = yproj - yp11
    dxpx = xp21 - xp11
    dxpy = xp12 - xp11
    dypx = yp21 - yp11
    dypy = yp12 - yp11
    den = dxpx * dypy - dxpy * dypx
    fx = (xErr * dypy - yErr * dxpy) / den
    fy = (dxpx * yErr - dypx * xErr) / den
    ! write(*,101) 'dx,err,f', dxpx, dxpy, dypx, dypy, den, xerr, yerr, fx, fy
    xx = ix + fx
    yy = iy + fy
    if (fx > -0.1 .and. fx < 1.1 .and. fy > -0.1 .and. &
        fy < 1.1) ifDone = 1
    if (ifOut .ne. 0 .and. (iter > 0 .or.  xx < 0. .or. &
        xx > nxLoad + 1 .or. yy < inLoadStart - 1. .or. &
        yy > inLoadEnd + 1.)) ifDone = 1
    iter = iter + 1
  enddo
  return
end subroutine loadedProjectingPoint

! reprojOneAngle reprojects a line at one angle from projection data
! in ARRAY into reprojLines.  LINE is the Y value in projections, Z value
! in reconstructed slices.  ARRAY is loaded with slices from inLoadStart to
! inLoadEnd, with a slice size of inPlaneSize and NXLOAD values on each line
! in X.  COSBET, SINBET, COSALPH, SINALPH are cosines and sines of tilt angle
! and alpha tilt.  IFALPHA is non-zero for tilt araound the X axis.
! iwidth is the width and ithickReproj is the thickness to reprojection;
! the coordinates of the region in the slice to reproject start at
! minXreproj, minYreproj and are offset from center by xprojOffset,
! yprojOffset.  DELZIN is the spacing in Y at which to sample points along
! a projection ray. XCENout, YCENout are center coordinates of output;
! xcenPdelxx is xcenIn + axisXoffset; centerSlice is the middle coordinate in Z
! XZFACVIEW, YZFACVIEW are Z factors for this view, and PMEAN is the mean of
! the slices for filling.
! Many of these names are the same as in the rest of the program and
! have the same meaning, but they are passed in, not taken from the
! tiltvars module, so that this can be called in cases other than the
! reprojectRec where they are defined.
!
subroutine reprojOneAngle(array, reprojLines, inLoadStart, inLoadEnd, line, &
    cosBet, sinBet, cosAlph, sinAlph, delzIn, iwidth, ithickReproj, inPlaneSize, &
    nxLoad, minXreproj, minYreproj, xprojOffset, yprojOffset, xcenOut, ycenOut, &
    xcenPdelxx, centerSlice, ifAlpha, xzFacView, yzFacView, dmeanIn)
  implicit none
  real*4 array(*), reprojLines(*), cosBet, sinBet, cosAlph, sinAlph, delX, delzIn
  integer*4 inLoadStart, inLoadEnd, line, iwidth, ithickReproj, inPlaneSize, nxLoad
  integer*4 minXreproj, minYreproj, ifAlpha
  real*4 xprojOffset, yprojOffset, centerSlice, xzFacView, yzFacView, dmeanIn
  real*4 xcenOut, ycenOut, xcenPdelxx
  !
  integer*4 ix, iz, i, numZ, kz, iys, ixEnd, ixStart, ind, indBase
  real*4 zNum, fz, oneMfz, zz, xx, fx, yEndTol, pfill, salfSbetOverCalf, xcenAdj, yslc
  real*4 oneMfx, yy, fy, oneMfy, xproj, yproj, ySlice, d11, d12, d21, d22, delY
  real*8 xx8, yy8, zz8
  integer*4 numX, kx, ixyOKstart, ixyOKend, iyFix, ifixStart, ifixEnd
  real*4 delZ, xNum, eps

  yEndTol = 3.05
  eps = 0.01
  reprojLines(1:iwidth) = 0.
  if (abs(sinBet * ithickReproj) <= abs(cosBet * iwidth)) then
    !
    delZ = delzIn
    delX = 1. / cosBet
    zNum = 1. + (ithickReproj - 1) / delZ
    numZ = zNum
    if (zNum - numZ >= 0.1) numZ = numZ + 1
    !
    ! Loop up in Z through slices, adding in lines of data to the
    ! output line
    do kz = 1, numZ
      zz = 1 + (kz - 1) * delZ
      iz = zz
      fz = zz - iz
      oneMfz = 1. - fz
      pfill = dmeanIn
      !
      ! If Z is past the top, drop back one line and set up fractions
      ! to take just a fraction of the top line
      if (zz >= ithickReproj) then
        zz = ithickReproj
        iz = ithickReproj - 1
        fz = oneMfz
        oneMfz = 0.
        pfill = dmeanIn * fz
      endif
      zz = zz + minYreproj - 1 - ycenOut
      !
      ! Get y slice for this z value
      yproj = line + yprojOffset
      yy = (yproj + zz * (sinAlph - yzFacView) - centerSlice) / cosAlph
      ySlice = yy + centerSlice - yprojOffset
      if (ifAlpha == 0) ySlice = line
      ! if (line==591) print *,kz, zz, iz, fz, omfz, yproj, yy, yslice
      if (ySlice < inLoadStart - yEndTol .or. &
          ySlice > inLoadEnd + yEndTol) then
        !
        ! Really out of bounds, do fill
        ! if (line==591) print *,'Out of bounds, view, line, zz', line, zz
        reprojLines(1:iwidth) = reprojLines(1:iwidth) + pfill
      else
        !
        ! otherwise set up iy and interpolation factors
        iys = floor(ySlice)
        if (ifAlpha .ne. 0) then
          if (iys < inLoadStart) then
            iys = inLoadStart
            fy = 0.
          else if (iys >= inLoadEnd) then
            iys = inLoadEnd - 1
            fy = 1.
          else
            fy = ySlice - iys
          endif
          oneMfy = 1. - fy
        endif
        !
        ! Now get starting X coordinate, fill to left
        xproj = 1 + xprojOffset
        xx = (xproj - (yy * sinAlph * sinBet + zz * (cosAlph * sinBet + &
            xzFacView) + xcenPdelxx)) / cosBet + xcenOut - (minXreproj - 1)
        ixStart = 1
        if (xx < 1) then
          ixStart = ceiling((1. - xx) / delX + 1.)
        elseif (xx >= iwidth) then
          ixStart = ceiling((iwidth - xx) / delX + 1.)
        endif
        xx = xx + (ixStart - 1) * delX
        if (xx < 1 .or. xx >= iwidth) then
          ixStart = ixStart + 1
          xx = xx + delX
        endif
        if (ixStart > 1) reprojLines(1:ixStart - 1) = reprojLines(1:ixStart - 1) + pfill
        !
        ! get ending X coordinate, fill to right
        ixEnd = iwidth
        if (xx + (ixEnd - ixStart) * delX >= iwidth - eps) then
          ixEnd = (iwidth - xx) / delX + ixStart
          if (xx + (ixEnd - ixStart) * delX >= iwidth - eps) ixEnd = ixEnd - 1
        elseif (xx + (ixEnd - ixStart) * delX < 1 + eps) then
          ixEnd = (1. - xx) / delX + ixStart
          if (xx + (ixEnd - ixStart) * delX < 1 + eps) ixEnd = ixEnd - 1
        endif
        if (ixEnd < iwidth) &
            reprojLines(ixEnd + 1:iwidth) = reprojLines(ixEnd + 1:iwidth) + pfill

        ! if (line == lsStart) write(*,'(3i6,3f11.3)') iv, ixst, ixnd, xx, &
        ! (ixst+xprojOffset - &
        ! (yy * salf * sbeta + zz * (calf * sbeta + &
        ! xzfacv) + xcenPdelxx)) / cbeta + xcenOut - (minXreproj-1) &
        ! , (ixnd+xprojOffset - &
        ! (yy * salf * sbeta + zz * (calf * sbeta + &
        ! xzfacv) + xcenPdelxx)) / cbeta + xcenOut - (minXreproj-1)
        !
        ! Add the line in: do simple 2x2 interpolation if no alpha
        indBase = 1 + inPlaneSize * (iys - inLoadStart) + (iz - 1) * nxLoad
        ! if (line==591) print *,ixst, ixnd
        xx8 = xx
        if (ifAlpha == 0) then
          do i = ixStart, ixEnd
            ix = xx8
            fx = xx8 - ix
            oneMfx = 1. - fx
            ind = indBase + ix - 1
            reprojLines(i) = reprojLines(i) + &
                oneMfz * oneMfx * array(ind) + &
                oneMfz * fx * array(ind + 1) + &
                fz * oneMfx * array(ind + nxLoad) + &
                fz * fx * array(ind + nxLoad + 1)
            ! if (line==591.and.i==164) print *,reprojLines(i), array(ind), &
            ! array(ind + 1), array(ind + nxload), array(ind + nxload + 1)
            xx8 = xx8 + delX
          enddo
        else
          !
          ! Or do the full 3D interpolation if any variation in Y
          do i = ixStart, ixEnd
            ix = xx8
            fx = xx8 - ix
            oneMfx = 1. - fx
            d11 = oneMfx * oneMfy
            d12 = oneMfx * fy
            d21 = fx * oneMfy
            d22 = fx * fy
            ind = indBase + ix - 1
            reprojLines(i) = reprojLines(i) + &
                oneMfz * (d11 * array(ind) &
                + d12 * array(ind + inPlaneSize) + d21 * array(ind + 1) &
                + d22 * array(ind + inPlaneSize + 1)) &
                + fz * (d11 * array(ind + nxLoad) &
                + d12 * array(ind + inPlaneSize + nxLoad) &
                + d21 * array(ind + 1 + nxLoad) &
                + d22 * array(ind + inPlaneSize + 1 + nxLoad))
            xx8 = xx8 + delX
          enddo
        endif
      endif
    enddo

  else
    !
    ! angles higher than the corner angle need to be done in vertical lines,
    ! outer loop on X instead of z
    ! Spacing between vertical lines is now sine beta
    ! The step between pixels along a line is 1/sin beta with no alpha tilt,
    ! The alpha tilt compresses it by the delta Z factor divided by cosine
    ! beta, the amount that delta Z factor is compressed from cosine beta.
    delX = abs(sinBet)
    xNum = 1. + (iwidth - 1) / delX
    numX = xNum
    if (xNum - numX >= 0.1) numX = numX + 1
    delZ = delzIn / (sinBet * abs(cosBet))
    delY = delZ * (sinAlph - yzFacView) / cosAlph
    if (abs(delY) < 1.e-10) delY = 0.
    ! print *,'delx, dely, delz', delx, dely, delz
    !
    ! Loop in X across slices, adding in vertical lines of data to the output line
    do kx = 1, numX
      xx = 1 + (kx - 1) * delX
      ix = xx
      fx = xx - ix
      oneMfx = 1. - fx
      pfill = dmeanIn
      !
      ! If X is past the end, drop back one line and set up fractions
      ! to take just a fraction of the right column
      if (xx >= iwidth) then
        xx = iwidth
        ix = iwidth - 1
        fx = oneMfx
        oneMfx = 0.
        pfill = dmeanIn * fx
      endif

      ! get starting Z coordinate
      salfSbetOverCalf = sinAlph * sinBet / cosAlph
      xcenAdj = xcenOut - (minXreproj - 1)
      xproj = 1 + xprojOffset
      yproj = line + yprojOffset

      zz = (xproj - (yproj - centerSlice) * salfSbetOverCalf - xcenPdelxx -  &
          (xx - xcenAdj) * cosBet) / ((sinAlph - yzFacView) * salfSbetOverCalf +  &
          cosAlph * sinBet + xzFacView)
      !
      ! Get y slice for this z value, then convert Z to be index coordinates in slice
      yy = (yproj + zz * (sinAlph - yzFacView) - centerSlice) / cosAlph
      ySlice = yy + centerSlice - yprojOffset
      if (ifAlpha == 0) ySlice = line
      zz = zz - (minYreproj - 1 - ycenOut)
      !
      ! Get starting X proj limit based on Z
      ixStart = 1
      if (zz < 1.) then
        ixStart = ceiling((1. - zz) / delZ + 1.)
      elseif (zz >= ithickReproj) then
        ixStart = ceiling((ithickReproj - zz) / delZ + 1.)
      endif
      !
      ! Revise starting limit for Y
      if (ifAlpha .ne. 0 .and. delY .ne. 0) then
        yslc = ySlice + (ixStart - 1.) * delY
        if (yslc < inLoadStart - yEndTol) then
          ixStart = ceiling((inLoadStart - yEndTol - ySlice) / delY + 1.)
        elseif (yslc > inLoadEnd + yEndTol) then
          ixStart = ceiling((inLoadEnd + yEndTol - ySlice) / delY + 1.)
        endif
      endif
      !
      ! Adjust Z start for final start and make sure it works, adjust Y also
      zz = zz + (ixStart - 1.) * delZ
      if (zz < 1. .or. zz >= ithickReproj) then
        zz = zz + delZ
        ixStart = ixStart + 1
      endif
      ySlice = ySlice + (ixStart - 1.) * delY
      if (ifAlpha == 0) ySlice = line
      !
      ! get ending coordinate based on limits in Z and Y
      ixEnd = iwidth
      if (zz + (ixEnd - ixStart) * delZ >= ithickReproj - eps) then
        ixEnd = (ithickReproj - zz) / delZ + ixStart
        if (zz + (ixEnd - ixStart) * delZ >= ithickReproj - eps) ixEnd = ixEnd - 1
      elseif ( zz + (ixEnd - ixStart) * delZ < 1. + eps) then
        ixEnd = (1. - zz) / delZ + ixStart
        if (zz + (ixEnd - ixStart) * delZ < 1. + eps) ixEnd = ixEnd - 1
      endif
      if (ifAlpha .ne. 0 .and. delY .ne. 0) then
        yslc = ySlice + (ixEnd - ixStart) * delY
        if (yslc < inLoadStart - yEndTol) then
          ixEnd = (inLoadStart - yEndTol - ySlice) / delY + ixStart
        elseif (yslc > inLoadEnd + yEndTol) then
          ixEnd = (inLoadEnd + yEndTol - ySlice) / delY + ixStart
        endif
      endif
      !
      ! Now get X indexes within which Y can safely be varied
      ixyOKstart = ixStart
      ixyOKend = ixEnd
      if (ifAlpha .ne. 0 .and. delY .ne. 0.) then
        if (ySlice < inLoadStart) then
          ixyOKstart = ceiling((inLoadStart - ySlice) / delY + ixStart)
          if (ySlice + (ixyOKstart - ixStart) * delY < inLoadStart + eps) &
              ixyOKstart = ixyOKstart + 1
        elseif (ySlice >= inLoadEnd) then
          ixyOKstart = ceiling((inLoadEnd - ySlice) / delY + ixStart)
          if (ySlice + (ixyOKstart - ixStart) * delY >= inLoadEnd - eps) &
              ixyOKstart = ixyOKstart + 1
        endif
        ySlice = ySlice + (ixyOKstart - ixStart) * delY
        !
        yslc = ySlice + (ixEnd - ixyOKstart) * delY
        if (yslc < inLoadStart) then
          ixyOKend = (inLoadStart - ySlice) / delY + ixyOKstart
          if (ySlice + (ixyOKend - ixyOKstart) * delY < inLoadStart + eps) &
              ixyOKend = ixyOKend-1
        elseif (yslc >= inLoadEnd) then
          ixyOKend = (inLoadEnd - ySlice) / delY + ixyOKstart
          if (ySlice + (ixyOKend - ixyOKstart) * delY >= inLoadEnd - eps) &
              ixyOKend = ixyOKend-1
        endif
      endif
      ! write( *,'(i5,f7.1,4i5,2f7.1)') kx, xx, ixst, ixyOKst, ixyOKnd, ixnd, zz, &
      ! zz+(ixnd-ixst)*delz
      !
      ! Do the fills
      if (ixStart > 1) reprojLines(1:ixStart - 1) = reprojLines(1:ixStart - 1) + pfill
      if (ixEnd < iwidth) &
          reprojLines(ixEnd + 1:iwidth) = reprojLines(ixEnd + 1:iwidth) + pfill
      !
      ! Add the line in: do simple 2x2 interpolation if no alpha
      ! if (line==591) print *,ixst, ixnd
      if (ifAlpha == 0) then
        zz8 = zz
        indBase = 1 + inPlaneSize * (line - inLoadStart) + ix - 1
        do i = ixStart, ixEnd
          iz = zz8
          fz = zz8 - iz
          oneMfz = 1. - fz
          ind = indBase + (iz - 1) * nxLoad
          reprojLines(i) = reprojLines(i) + &
              oneMfz * oneMfx * array(ind) + &
              oneMfz * fx * array(ind + 1) + &
              fz * oneMfx * array(ind + nxLoad) + &
              fz * fx * array(ind + nxLoad + 1)
          ! if (i==70) print *,reprojLines(i)
          ! if (line==591.and.i==164) print *,reprojLines(i), array(ind), &
          ! array(ind + 1), array(ind + nxload), array(ind + nxload + 1)
          zz8 = zz8 + delZ
        enddo
      else
        !
        ! Or do the full 3D interpolation if any variation in Y, starting with
        ! the loop where Y varies
        yy8 = ySlice
        zz8 = zz + (ixyOKstart - ixStart) * delZ
        indBase = 1 - inPlaneSize * inLoadStart + ix - 1
        do i = ixyOKstart, ixyOKend
          iz = zz8
          fz = zz8 - iz
          oneMfz = 1. - fz
          iys = yy8
          fy = yy8 - iys
          oneMfy = 1. - fy
          d11 = oneMfx * oneMfy
          d12 = oneMfx * fy
          d21 = fx * oneMfy
          d22 = fx * fy
          ind = indBase + inPlaneSize * iys + (iz - 1) * nxLoad
          reprojLines(i) = reprojLines(i) + &
              oneMfz * (d11 * array(ind) &
              + d12 * array(ind + inPlaneSize) + d21 * array(ind + 1) &
              + d22 * array(ind + inPlaneSize + 1)) &
              + fz * (d11 * array(ind + nxLoad) &
              + d12 * array(ind + inPlaneSize + nxLoad) &
              + d21 * array(ind + 1 + nxLoad) &
              + d22 * array(ind + inPlaneSize + 1 + nxLoad))
          zz8 = zz8 + delZ
          yy8 = yy8 + delY
        enddo
        !
        ! Now do special loops with Y fixed - do the one at the end first
        ! since Y and Z are all set for that
        ifixStart = ixyOKend + 1
        ifixEnd = ixEnd
        do iyFix = 1, 2
          do i = ifixStart, ifixEnd
            iz = zz8
            fz = zz8 - iz
            oneMfz = 1. - fz
            d11 = oneMfx * oneMfy
            d12 = oneMfx * fy
            d21 = fx * oneMfy
            d22 = fx * fy
            ind = indBase + inPlaneSize * iys + (iz - 1) * nxLoad
            reprojLines(i) = reprojLines(i) + &
                oneMfz * (d11 * array(ind) &
                + d12 * array(ind + inPlaneSize) + d21 * array(ind + 1) &
                + d22 * array(ind + inPlaneSize + 1)) &
                + fz * (d11 * array(ind + nxLoad) &
                + d12 * array(ind + inPlaneSize + nxLoad) &
                + d21 * array(ind + 1 + nxLoad) &
                + d22 * array(ind + inPlaneSize + 1 + nxLoad))
            zz8 = zz8 + delZ
          enddo
          !
          ! Set up for loop with Y fixed at start, reset y and z
          yy8 = ySlice
          zz8 = zz
          iys = yy8
          fy = yy8 - iys
          oneMfy = 1. - fy
          ifixStart = ixStart
          ifixEnd = ixyOKstart - 1
        enddo
      endif
    enddo
  endif
  return
end subroutine reprojOneAngle

! Computes the change in Z that moves by 1 pixel along a projection
! ray given the sines and cosines of alpha and beta and the z factors
!
real*4 function reprojDelz(sinBet, cosBet, sinAlph, cosAlph, xzFac, yzFac)
  implicit none
  real*4 sinBet, cosBet, sinAlph, cosAlph, xzFac, yzFac,  dyFac
  dyFac = (sinAlph - yzFac) / cosAlph
  reprojDelz = 1. / sqrt(1. + dyFac**2 + &
      ((dyFac * sinAlph * sinBet + cosAlph * sinBet + xzFac) / cosBet)**2)
  return
end function reprojDelz


! Writes line LINE for view IV of a reprojection
!
subroutine writeReprojLines(iv, lineStart, lineEnd, dmin, dmax, dtot8)
  use tiltvars
  implicit none
  integer*4 line, i, iyOut, iv, lineStart, lineEnd, numVals
  real*4 dmin, dmax, val
  real*8 dtot8
  !
  ! Write the line after scaling.  Scale log data to give approximately
  ! constant mean levels.  Descale non-log data by exposure weights
  numVals = iwidth * (lineEnd + 1 - lineStart)
  if (ifLog .ne. 0) then
    ! Hopefully this works for local as well
    if (threshPolarity .ne. 0.) then
      val = alog10(projMean + baseForLog) - ithickReproj * dmeanIn
    else if (abs(sinBeta(iv) * ithickReproj) <= abs(cosBeta(iv) * iwidth)) then
      val = alog10(projMean + baseForLog) - ithickReproj * dmeanIn / abs(cosBeta(iv))
    else
      val = alog10(projMean + baseForLog) - iwidth * dmeanIn / abs(sinBeta(iv))
    endif
    if (debug) print *,iv, lineStart, lineEnd, val
    do i = 1, numVals
      reprojLines(i) = 10**(reprojLines(i) + val) - baseForLog
    enddo
  else
    val = exposeWeight(iv)
    if (threshPolarity .ne. 0.) val = exposeWeight((numViews + 1) / 2)
    do i = 1, numVals
      reprojLines(i) = reprojLines(i) / val
    enddo
  endif
  if (projSubtraction) then
    call iiuSetPosition(1, iv - 1, lineStart - 1)
    call irdsecl(1, origLines, lineEnd + 1 - lineStart, *99)
    reprojLines(1:numVals) = reprojLines(1:numVals) - origLines(1:numVals)
  endif
  do i = 1, numVals
    val = reprojLines(i)
    ! if (debug .and. val < dmin) print *,'min:', i, val
    ! if (debug .and. val > dmax) print *,'max:', i, val
    dmin = min(dmin, val)
    dmax = max(dmax, val)
    dtot8 = dtot8 + val
  enddo
  do line = lineStart, lineEnd
    iyOut = line - isliceStart
    if (minTotSlice > 0) iyOut = line - minTotSlice
    call parWrtPosn(2, iv - 1, iyOut)
    call parWrtLin(2, reprojLines(1 + (line - lineStart) * iwidth))
  enddo
  return
99 call exitError('READING FROM ORIGINAL PROJECTION FILE')
end subroutine writeReprojLines


! Projects model points onto the included views
!
subroutine projectModel(outModel, outAngles, transformFile, delta, numViewsOrig,  &
    defocusFile, pixForDefocus, focusInvert)
  use tiltvars
  implicit none
  include 'model.inc90'
  character*(*) outModel, outAngles, transformFile, defocusFile
  real*4 delta(3), origin(3), pixForDefocus, focusInvert
  character*120 objName
  integer*4 numViewsOrig, ibase, numPoints, iobj, ipt, ip1, iv, nfv, numXfs
  real*4 oneVal, rlJ, rlI, rlSlice, yy, xproj, yproj, zOffset, f11, f12, f21, f22
  integer*4 j, lslice, imodObj, imodCont, ierr, size
  real*4 yz22, xprojf, xprojz, yprojf, yprojz, degToRad, ptDefocus, betaInv
  real*4 fjlMat(2,2), f1234(4), gamma, stretch, strPhi, smag, gammaSum, betaSum, alphaSum
  integer*4 ind1234(4), jInc, lsInc, ic
  equivalence (fjlMat(1, 1), f11), (fjlMat(1, 2), f12), (fjlMat(2, 1), f21),  &
      (fjlMat(2, 2), f22)
  data degToRad/0.0174533/
  real*4, allocatable :: values(:), coords(:,:), alixf(:,:,:), aliRot(:), defocus(:)
  integer*4, allocatable :: mapFileToView(:)
  integer*4 getContValue, putImageRef, putContValue, putImodFlag
  integer*4 getScatSize, putScatSize, putImodMaxes, getImodObjName, getZFromMinusPt5
  !
  call iiuRetOrigin(1, origin(1), origin(2), origin(3))
  call scale_model(0)
  if (getScatSize(1, size) .ne. 0) size = 5
  !
  ! Models are defined as having Z coordinates ranging from -0.5 to NZ - 0.5 so Z
  ! will need to be shifted up by 1 to get to pixel index coordinates
  ! But in old beadtrack models, Z started at 0 so the offset needs to be only 0.5
  ! Recognize old beadtrack model by lack of flag AND original object name so that
  ! other old models will work
  zOffset = 1.0
  if (getZFromMinusPt5() == 0 .and. getImodObjName(1, objName) == 0) then
    if (objName(1:8) == 'Wimp no.') zOffset = 0.5
  endif
  j = numViewsOrig + 10
  allocate(values(n_point), coords(3, n_point), mapFileToView(limView),  &
      alixf(2, 3, j), aliRot(j), defocus(j), stat = ierr)
  call memoryErrorUC(ierr, 'ARRAYS FOR REPROJECTING MODEL')
  if (outAngles .ne. ' ') then
    call dopen(12, outAngles, 'new', 'f')
    if (transformFile .ne. ' ') then
      call dopen(13, transformFile, 'old', 'f')
      call xfrdall2(13, alixf, numXfs, numViewsOrig + 10, j)
      close(13)
      if (numXfs > numViewsOrig) write(*,'(/,a)') &
          'WARNING: TILT - MORE ALIGNMENT TRANSFORMS THAN VIEWS IN INPUT STACK'
      if (j == 2) call exitError('READING FILE OF ALIGNMENT TRANSFORMS')
      if (numXfs < numViewsOrig) call exitError( &
          'FEWER ALIGNMENT TRANSFORMS THAN VIEWS IN INPUT STACK')
      do j = 1, numViewsOrig
        call amat_to_rotmagstr(alixf(1, 1, j), aliRot(j), smag, stretch, strPhi)
      enddo
    endif
    !
    if (defocusFile .ne. ' ') then
      call dopen(13, defocusFile, 'old', 'f')
      read(13,*,err = 98, end = 98) (defocus(j), j = 1, numViewsOrig)
      close(13)
      if (focusInvert == 0) then
        focusInvert = 1.
      else
        focusInvert = -1.
      endif
    endif
  endif
  !
  ! get each point and its contour value into the arrays
  numPoints = 0
  do iobj = 1, max_mod_obj
    call objToCont(iobj, obj_color, imodObj, imodCont)
    if (getContValue(imodObj, imodCont, oneVal) .ne. 0) oneVal = -1.
    ibase = ibase_obj(iobj)
    do ipt = 1, npt_in_obj(iobj)
      numPoints = numPoints + 1
      values(numPoints) = oneVal
      ip1 = abs(object(ipt + ibase))
      coords(1, numPoints) = p_coord(1, ip1)
      coords(2, numPoints) = p_coord(2, ip1)
      coords(3, numPoints) = p_coord(3, ip1)
    enddo
  enddo
  !
  ! Start a new model
  call newImod()
  n_point = 0
  iobj = 0
  if (putImageRef(delta, origin) .ne. 0 .or. putImodMaxes(nprojXyz(1), nprojXyz(2),  &
      nprojXyz(3)) .ne. 0) call exitError( &
      'Putting image reference or maximum size information in output model')
  !
  ! Build a map from views in file to ordered views in program
  do nfv = 1, numViewsOrig
    mapFileToView(nfv) = 0
  enddo
  do iv = 1, numViews
    mapFileToView(mapUsedView(iv)) = iv
  enddo
  !
  ! Loop on the points, start new contour for each
  do ipt = 1, numPoints
    iobj = iobj + 1
    obj_color(1, iobj) = 1
    obj_color(2, iobj) = 255
    ierr = putContValue(1, iobj, values(ipt))
    ibase_obj(iobj) = n_point
    npt_in_obj(iobj) = 0
    !
    ! Get real pixel coordinates in tomogram file
    rlJ = coords(1, ipt) + 0.5
    rlI = coords(2, ipt) + 0.5
    rlSlice = coords(3, ipt) + zOffset
    !
    ! This may never be tested but seems simple enough
    if (.not.perpendicular) then
      rlI = coords(3, ipt) + zOffset
      rlSlice = coords(2, ipt) + 0.5
    endif
    !
    ! Loop on the views in the file
    do nfv = 1, numViewsOrig
      iv = mapFileToView(nfv)
      if (iv > 0) then
        call projectionPosition(iv, rlJ, rlI, rlSlice, xproj, yproj, j, lslice,  &
            f11, f12, f21, f22)
        if (outAngles .ne. ' ') then
          betaSum = 0.
          alphaSum = 0.
          gammaSum = 0.
          if (nxWarp == 0) then
            alphaSum = alpha(iv)
            betaSum = -angles(iv) / degToRad
          else
            do jInc = 1, 2
              do lsInc = 1, 2
                call local_factors(nint(j + jInc - 1 - xcenOut + xcenIn + axisXoffset), &
                    lslice + lsInc - 1,  nfv, ind1234(1), ind1234(2), &
                    ind1234(3), ind1234(4), f1234(1), f1234(2), f1234(3), f1234(4))
                f1234(1:4) = f1234(1:4) * fjlMat(jInc, lsInc)
                do ic = 1, 4
                  call amat_to_rotmagstr(fwarp(1, 1, ind1234(ic)), gamma, smag, &
                      stretch, strPhi)
                  !
                  ! alpha was left as degrees and with its native sign, with the negative
                  ! taken when taking the sin; but tilt angle and delBeta were converted
                  ! to radians and made negative
                  ! The transform has the amount image needs to be rotated; need
                  ! negative for rotation of specimen
                  alphaSum = alphaSum + f1234(ic) * (alpha(iv) + delAlpha(ind1234(ic)))
                  betaSum = betaSum - f1234(ic) *  (angles(iv) + delBeta(ind1234(ic))) / &
                      degToRad
                  gammaSum = gammaSum - f1234(ic) * gamma
                enddo
              enddo
            enddo
          endif
          if (transformFile .ne. ' ') gammaSum = gammaSum - aliRot(nfv)
          ptDefocus = 0.
          if (defocusFile .ne. ' ') then
            yy = (yproj - centerSlice) * tan(alphaSum * degToRad)
            betaInv = focusInvert * betaSum * degToRad
            ptDefocus = defocus(nfv) - ((xproj - xcenIn) - yy * sin(betaInv)) * &
                tan(betaInv) + yy * cos(betaInv)
          endif
          write(12, '(3i6,3f9.3,f10.0)') iobj, iv, nfv, alphaSum, betaSum, gammaSum,  &
              ptDefocus
        endif
        !
        ! Store model coordinates
        n_point = n_point + 1
        if (n_point > max_pt) call exitError( &
            'Too many projection points for small model arrays')
        npt_in_obj(iobj) = npt_in_obj(iobj) + 1
        object(n_point) = n_point
        p_coord(1, n_point) = xproj - 0.5
        p_coord(2, n_point) = yproj - 0.5
        p_coord(3, n_point) = nfv - 1.
      endif
    enddo
  enddo
  !
  ! Save model
  max_mod_obj = iobj
  !
  ! Set to open contour, show values etc., and show sphere on section only
  ierr = putImodFlag(1, 1)
  ierr = putImodFlag(1, 7)
  ierr = putImodFlag(1, 9)
  ierr = putScatSize(1, size)
  call scale_model(1)
  call write_wmod(outModel)
  print *,n_point, ' points written to output model'
  if (outAngles .ne. ' ') close(12)
  call exit(0)
98 call exitError('READING FILE OF DEFOCUS VALUES')
end subroutine projectModel

! Computes the projection position xproj,yproj on view iv of position rlX, rlZ,
! rlSlice in the reconstruction, where these are real pixel coordinates, equal to 1
! in the middle of the first pixel.  Other returned values are the rounded down X and 
! slice numbers in J and lslice, and the four interpolation factors for the real pixel
! position
!
subroutine projectionPosition(iv, rlX, rlZ, rlSlice, xproj, yproj, j, lslice, &
    f11, f12, f21, f22)
  use tiltvars
  implicit none
  real*4 rlX, rlZ, rlSlice, zz, yy, zpart, xproj, yproj
  integer*4 j, lslice, iv
  real*4 fj, fls, f11, f12, f21, f22, xf11, xz11, yf11, yz11
  real*4 xf21, xz21, yf21, yz21, xf12, xz12, yf12, yz12, xf22, xz22, yf22
  real*4 yz22, xprojf, xprojz, yprojf, yprojz
  
  zz = (rlZ - ycenModProj) * compress(iv)
  yy = rlSlice - centerSlice
  if (nxWarp == 0) then
    zpart = yy * sinAlpha(iv) * sinBeta(iv) + zz * (cosAlpha(iv) * sinBeta(iv) +  &
        xzfac(iv)) + xcenIn + axisXoffset
    yproj = yy * cosAlpha(iv) - zz * (sinAlpha(iv) - yzfac(iv)) + centerSlice
    xproj = zpart + (rlX - xcenOut) * cosBeta(iv)
  else
    !
    ! local alignments
    j = rlX
    fj = rlX - j
    lslice = rlSlice
    fls = rlSlice - lslice
    f11 = (1. -fj) * (1. -fls)
    f12 = (1. -fj) * fls
    f21 = fj * (1. -fls)
    f22 = fj * fls
    call localProjFactors(j, lslice, iv, xf11, xz11, yf11, yz11)
    call localProjFactors(j + 1, lslice, iv, xf21, xz21, yf21, yz21)
    call localProjFactors(j, lslice + 1, iv, xf12, xz12, yf12, yz12)
    call localProjFactors(j + 1, lslice + 1, iv, xf22, xz22, yf22, yz22)
    xprojf = f11 * xf11 + f12 * xf12 + f21 * xf21 + f22 * xf22
    xprojz = f11 * xz11 + f12 * xz12 + f21 * xz21 + f22 * xz22
    yprojf = f11 * yf11 + f12 * yf12 + f21 * yf21 + f22 * yf22
    yprojz = f11 * yz11 + f12 * yz12 + f21 * yz21 + f22 * yz22
    xproj = xprojf + zz * xprojz
    yproj = yprojf + zz * yprojz
    !
  endif
end subroutine projectionPosition
