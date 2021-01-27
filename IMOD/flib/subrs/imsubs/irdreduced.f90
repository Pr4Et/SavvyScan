! IRDREDUCED: read in image reduced in size with anitalias filter
!
! $Id$
!
!!
! Reads a specified portion of a section while shrinking the input by an arbitrary
! amount with antialias reduction.  Arguments are:  ^
! [imUnit] - image file unit number  ^
! [iz]     - section #, numbered from 0  ^
! [array]  - array in which to return image  ^
! [nxDim]  - X dimension of array (set to [nxRed] for a packed, one-D array)  ^
! [xUBstart] - starting unbinned X coordinate, 0 is at left edge of first pixel  ^
! [yUBstart] - starting unbinned Y coordinate, 0 is at bottom edge of first pixel  ^
! [redFac]   - factor to reduce by, must by > 1  ^
! [nxRed]    - reduced size of image in X  ^
! [nyRin]    - reduced size of image in Y  ^
! [ifiltType] - Type of filter (0-5), as defined by selectZoomFilter  ^
! [temp]      - temporary real array for reading unbinned data into  ^
! [lenTemp]   - Size of [temp] array  ^
! [ierr]      - return value, 0 if succeeded, -1 on read error, 1 if filter type is out
! of range, 2 for zoom out of range, 3 for insufficient temporary array space, 4 for
! starting or ending unbinned coordinates out of range, or 5 for memory allocation error ^
!
! As with irdBinned, [xUBstart] and [yUBstart] can be negative as long as they are 
! greater than or equal to 1 - [redFac], and the output extent can also be up to
! [redFac] - 1 beyond what can be produced by the input array.  The extra pixels on any
! edge will be produced by binning of available pixels (matching the output of irdBinned)
! if redFac is integral; otherwise they will be copied from adjacent pixels.
! The temporary array needs to be big enough to fit enough input to compose 10 lines of
! output, but data copies will be reduced by larger arrays.
!!
subroutine irdReduced(imUnit, iz, array, nxDim, xUBstart, yUBstart, redFac, nxRed, &
    nyRed, ifiltType, temp, lenTemp, ierr)
  implicit none
  integer*4 imUnit, nxDim, ix0, ix1, iy0, iy1, nx, ny, ifiltType
  integer*4 lenTemp, nxRed, nyRed, ierr, iz, nxyz(3), mxyz(3), nxyzst(3)
  real*4 xUBstart, yUBstart, redFac, chunkYstart, zoomFac, xUseStart, yUseStart
  real*4 array(nxDim, *), temp(lenTemp)
  integer*4 nxLoad, nyLoad, maxLineLoad, loadYstart, ifiltWidth, ihalfWidth, ix
  integer*4 iyStart, maxChunkLines, lastY1, lastY0, indStart, numCopy, iyEnd, ierr2
  integer*4 ibXoffset, ibYoffset, nxRedUse, nyRedUse, nbin, loadXoffset, loadXextra
  integer*4 loadYoffset, loadYextra, iyBinStart, iyBinEnd, iyEdgeStart, iyEdgeOffset
  integer*4 ixEdgeOffset
  logical fillXend, fillYend
  integer*4 selectZoomFilter, zoomWithFilter, iiuReadSecPart
  !
  call iiuRetSize(imUnit, nxyz, mxyz, nxyzst)
  nx = nxyz(1)
  ny = nxyz(2)
  !
  ! Check the validity of the limits here and get the adjusted limits to use
  zoomFac = 1. / redFac
  call irdRedSizesForLoad(xUBstart, redFac, nxRed, nx, xUseStart, nxRedUse,  &
      ibXoffset, fillXend, nbin, loadXoffset, loadXextra, ixEdgeOffset, ierr)
  if (ierr .ne. 0) return
  call irdRedSizesForLoad(yUBstart, redFac, nyRed, ny, yUseStart, nyRedUse,  &
      ibYoffset, fillYend, nbin, loadYoffset, loadYextra, iyEdgeOffset, ierr)
  if (ierr .ne. 0) return
  iyEdgeStart = 1
  !
  ! Set up the filter and get its width
  ierr = selectZoomFilter(ifiltType, zoomFac, ifiltWidth)
  if (ierr .ne. 0) return

  ! get the limits of the load in X and the maximum # of lines that can be done in one
  ! chunk; insist on being able to do 10 lines per load
  ihalfWidth = (ifiltWidth + 3) / 2
  ix0 = max(0, floor(xUseStart - ihalfWidth))
  ix1 = min(nx - 1, ceiling(xUseStart + redFac * nxRedUse + ihalfWidth))
  nxLoad = ix1 + 1 - ix0
  maxLineLoad = lenTemp / nxLoad
  maxChunkLines = zoomFac * (maxLineLoad - 2 * ihalfWidth)
  ierr = 3
  !print *,lenTemp, maxLineLoad, maxChunkLines, nyRedUse/maxChunkLines
  if ((redFac <= 32 .and. maxChunkLines < 10) .or. maxChunkLines < 2) return
  !
  ! Loop on chunks
  iyStart = 0
  lastY1 = -1
  do while (iyStart < nyRedUse)
    
    ! get limits of loaded data needed for next chunk, if the arithmetic is a bit off,
    ! reduce the end of the chunk to get within limits
    iyEnd = min(nyRedUse, iyStart + maxChunkLines)
    iy0 = max(0, floor(yUseStart + redFac * iyStart - ihalfWidth))
    iy1 = min(ny - 1, ceiling(yUseStart + redFac * iyEnd + ihalfWidth))
    do while (iy1 >= iy0 + maxLineLoad)
      iyEnd = iyEnd - 1
      iy1 = min(ny - 1, ceiling(yUseStart + redFac * iyEnd + ihalfWidth))
    enddo
    indStart = 1
    loadYstart = iy0

    ! if there is overlap, shift the data down and adjust the starting line and index
    ! But if we already have the last line, forget it and just read that again
    if (iy0 <= lastY1 .and. lastY1 < ny - 1) then
      indStart = (iy0 - lastY0) * nxLoad
      numCopy = (lastY1 + 1 - iy0) * nxLoad
      do ix = 1, numCopy
        temp(ix) = temp(ix + indStart)
      enddo
      loadYstart = lastY1 + 1
      indStart = numCopy + 1
    endif
    lastY0 = iy0
    lastY1 = iy1
    call iiuSetPosition(imUnit, iz, 0)
    ierr = -1
    ierr2 = iiuReadSecPart(imUnit, temp(indStart), nxLoad, ix0, ix1, loadYstart, iy1)
    if (ierr2 .ne. 0) return

    chunkYstart = yUseStart + iyStart * redFac - iy0
    ierr = zoomWithFilter(temp, nxLoad, iy1 + 1 - iy0, xUseStart - ix0, chunkYstart, &
        nxRedUse, iyEnd - iyStart, nxDim, ibXoffset, array(1, iyStart + ibYoffset + 1))
    if (ierr .ne. 0) return

    if (nbin > 0) then

      ! Set parameters and adjust for start or end if binning extra regions on edge
      iyEdgeStart = nint(chunkYstart) + 1
      if (iyStart == 0 .and. loadYoffset > 0) iyEdgeStart = 1
      iyBinStart = iyStart + ibYoffset + 1
      if (iyStart == 0) iyBinStart = 1
      iyBinEnd = iyEnd + ibYoffset
      if (iyEnd >= nyRedUse) iyBinEnd = nyRed

      ! Bottom row first time only
      if (iyStart == 0 .and. loadYoffset > 0) &
          call irdRedBinEdge(temp, nxLoad, iy1 + 1 - iy0, 1, ixEdgeOffset, 1, 0, &
          nbin, loadYoffset, array, nxDim, 1, nxRed, 1, 1)
      
      ! Left side always
      if (loadXoffset > 0) &
          call irdRedBinEdge(temp, nxLoad, iy1 + 1 - iy0, 1, 0, iyEdgeStart, &
          iyEdgeOffset, loadXoffset, nbin, array, nxDim, 1, 1, iyBinStart, iyBinEnd)
      !
      ! Right side always 
      if (loadXextra > 0) &
          call irdRedBinEdge(temp, nxLoad, iy1 + 1 - iy0, nxLoad + 1 - loadXextra, 0,  &
          iyEdgeStart, iyEdgeOffset, loadXoffset, nbin, array, nxDim, nxRed, nxRed, &
          iyBinStart, iyBinEnd)
      !
      ! Top edge last time
      if (iyEnd >= nyRedUse .and. loadYextra > 0) &
          call irdRedBinEdge(temp, nxLoad, iy1 + 1 - iy0, 1, ixEdgeOffset,  &
          iy1 + 2 - iy0 - loadYextra, 0, nbin, loadYextra, array, nxDim, 1, nxRed,  &
          nyRed, nyRed)
      iyEdgeOffset = 0
    endif
    iyStart = iyEnd
  enddo
  !
  ! If not binning edges, fill all the edges of the array by copying if necessary
  if (nbin == 0) then
    if (yUBstart < 0) array(1:nxRed, 1) = array(1:nxRed, 2)
    if (fillYend) array(1:nxRed, nyRed) = array(1:nxRed, nyRed - 1) 
    if (xUBstart < 0) array(1, 1:nyRed) = array(2, 1:nyRed)
    if (fillXend) array(nxRed, 1:nyRed) = array(nxRed - 1, 1:nyRed) 
  endif
  ierr = 0
  return
end subroutine irdReduced


! irdRedSizesForLoad checks whether the load parameters are legal for one dimension,
! allowing up to the reduction minus 1 of extra output on the either end
! If there is extra output, it adjusts the output side to use and the starting offset to
! use so that they will be legal for calling zoomWithFilter.  It sets ibXoffset with an
! offset for the array coordinates in which the output is placed by zoomWithFilter, and
! sets fillEnd if there is filling on the end of the array too; if reduction is integral,
! it sets nbin to this value and returns the number of pixels to be binned at the start
! and end in loadOffset and loadExtra, and a complementary offset needed for the
! edge in the opposite dimension
!
subroutine irdRedSizesForLoad(xUBstart, redFac, nxRed, nx, xUseStart, nxRedUse,  &
    ibXoffset, fillEnd, nbin, loadOffset, loadExtra, ixEdgeOffset, ierr)
  implicit none
  integer*4  nxRed, nx, nxRedUse, ierr, ibXoffset, nbin, loadOffset, loadExtra
  integer*4 ixEdgeOffset
  real*4 xUBstart, redFac, xUseStart
  logical fillEnd
  xUseStart = xUBstart
  nxRedUse = nxRed
  ibXoffset = 0
  nbin = 0
  loadOffset = 0
  loadExtra = 0
  ixEdgeOffset = 0
  fillEnd = .false.
  ierr = 4
  !
  ! Test for legality
  if (xUBstart < -(redFac - 0.99) .or.  &
      int(xUBstart + redFac * nxRed) > nx + redFac - 0.99) return
  ierr = 0
  !
  ! test for integer binning
  if (abs(nint(redFac) - redFac) < 1.e-4) nbin = nint(redFac)
  !
  ! Adjust for negative offset
  if (xUBstart < 0) then
    xUseStart = xUBstart + redFac
    nxRedUse = nxRed -1
    ibXoffset = 1
    if (nbin > 0) then
      loadOffset = nint(xUseStart)
      ixEdgeOffset = nbin - loadOffset
    endif
  endif
  !
  ! Adjust for extra pixels at end
  if (int(xUseStart + redFac * nxRedUse) > nx) then
    nxRedUse = nxRedUse - 1
    fillEnd = .true.
    if (nbin > 0) loadExtra = nx - int(xUseStart + redFac * nxRedUse)
  endif
  ! write(*,'(i5,f8.1,i5,f8.1,i5,4i3)')nx, xUBstart, nxRed, xUseStart, nxRedUse, &
  !    ibXoffset, nbin, loadOffset, loadExtra
  return
end subroutine irdRedSizesForLoad


! irdRedBinEdge bins data along one edge of the image with the factors nbinX and nbinY, 
! for the extent, starting at the actual index in the input array ixInStart, iyInStart,
! but with potential offsets for the first pixel of ixInOffset and iyInOffset.
! ixOutStart, ixOutEnd, iyOutStart, iyOutEnd define the limits of the output array to
! be filled.
! 
subroutine irdRedBinEdge(temp, nxIn, nyIn, ixInStart, ixInOffset, iyInStart, iyInOffset, &
    nbinX, nbinY, array, nxDimOut, ixOutStart, ixOutEnd, iyOutStart, iyOutEnd)
  implicit none 
  integer*4 nxIn, nyIn, ixInStart, ixInOffset, iyInStart, iyInOffset, nbinX, nbinY
  integer*4 nxDimOut, ixOutStart, ixOutEnd, iyOutStart, iyOutEnd
  real*4 array(nxDimOut, *), temp(nxIn, nyIn)
  integer*4 ix, iy, ixb, iyb, ixStart, ixEnd, iyStart, iyEnd, ixUse, iyUse
  
  !write(*,'(13i6)')nxIn, nyIn, ixInStart, ixInOffset, iyInStart, iyInOffset, nbinX,  &
  !    nbinY, nxDimOut, ixOutStart, ixOutEnd, iyOutStart, iyOutEnd
  iyStart = iyInStart - iyInOffset
  do iyb = iyOutStart, iyOutEnd
    ixStart = ixInStart - ixInOffset
    iyEnd = min(nyIn, iyStart + nbinY - 1)
    iyUse = max(1, iyStart)
    do ixb = ixOutStart, ixOutEnd
      ixEnd = min(nxIn, ixStart + nbinX - 1)
      ixUse = max(1, ixStart)
      array(ixb, iyb) = sum(temp(ixUse:ixEnd, iyUse:iyEnd)) / ((ixEnd + 1 - ixUse) * &
          (iyEnd + 1 - iyUse))
      ! write(*,'(6i5, f10.2)')ixb,iyb,ixUse, ixEnd, iyUse, iyEnd, array(ixb, iyb)
      ixStart = ixStart + nbinX
    enddo
    iyStart = iyStart + nbinY
  enddo
  return
end subroutine irdRedBinEdge
