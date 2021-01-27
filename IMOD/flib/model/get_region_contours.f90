! This file has some routines shared by corrsearch3d, findwarp, and refinematch
!
! GET_REGION_CONTOURS will read the MODELFILE as a small model, and
! extract contours.  It first determines whether the contours lie in
! X/Y or X/Z planes, setting IFFLIP to 1 if they are in Y/Z planes.
! It extracts the planar points into XVERT, YVERT, where INDVERT is an
! index to the start of each contour, NVERT contains the number of
! points in each contour, ZCONT is the contour Y or Z value, NCONT is
! the number of contours.  LIMCONT and LIMVERT specify the limiting
! dimensions for contours and vertices.  Pass LIMCONT as <= 0 on an initial call
! to have them returned with the needed sizes.  If IMUNIT is greater than 0
! the model is scaled to the index coordinates of the image file open
! on that unit
!
! If contours are not planar in either direction it issues a warning
! using PROGNAME as the program name and and adopts the flip value
! from the model header, which is not necessarily correct.
!
! $Id$
!
subroutine get_region_contours(modelFile, progname, xVerts, yVerts, numVerts, &
    indVert, contourZ,  numConts, ifFlip, limCont, limVert, imUnit)
  implicit none
  include 'smallmodel.inc90'
  character*(*) modelFile, progname
  real*4 xVerts(*), yVerts(*), contourZ(*)
  integer*4 numVerts(*), indVert(*), numConts, ifFlip, limCont, limVert, imUnit
  integer*4 ip, ipt, iobj, iy, iz, indY, indZ, indCur, ierr, getImodHead
  logical*4 exist, readSmallMod, planarInZ, planarInY
  real*4 xyScale, zScale, xOffset, yOffset, zOffset
  !
  exist = readSmallMod(modelFile)
  if (.not.exist) call exitError('Error reading model')
  if (imUnit > 0) then
    call scaleModelToImage(imUnit, 0)
  else
    call scale_model(0)
  endif
  numConts = 0
  !
  ! figure out if all contours are coplanar Z or Y
  ! Loop through contours until both planar flags go false or all done
  !
  iobj = 1
  planarInZ = .true.
  planarInY = .true.
  do while (iobj <= max_mod_obj .and. (planarInZ .or. planarInY))
    if (npt_in_obj(iobj) >= 3) then
      ip = 2
      ipt = abs(object(ibase_obj(iobj) + 1))
      iy = nint(p_coord(2, ipt))
      iz = nint(p_coord(3, ipt))
      do while (ip <= npt_in_obj(iobj) .and. (planarInZ .or. planarInY))
        ipt = abs(object(ibase_obj(iobj) + ip))
        if (nint(p_coord(2, ipt)) .ne. iy) planarInY = .false.
        if (nint(p_coord(3, ipt)) .ne. iz) planarInZ = .false.
        ip = ip + 1
      enddo
    endif
    iobj = iobj + 1
  enddo
  !
  ! Set flip flag or fallback to what is in model header
  !
  if (planarInZ .and. .not. planarInY) then
    ifFlip = 0
  else if (planarInY .and. .not. planarInZ) then
    ifFlip = 1
  else
    ierr = getImodHead(xyScale, zScale, xOffset, yOffset, zOffset, ifFlip)
    if (limCont > 0) write(*,'(/,a,a,a)') 'WARNING: ', progname,  &
        ' - CONTOURS NOT ALL COPLANAR, USING HEADER FLIP FLAG'
  endif
  indY = 2
  indZ = 3
  if (ifFlip .ne. 0) then
    indY = 3
    indZ = 2
  endif
  !
  ! If limCont <= 0, add up number of contours and points and return
  if (limCont <= 0) then
    limCont = 0
    limVert = 0
    do iobj = 1, max_mod_obj
      if (npt_in_obj(iobj) >= 3) then
        limCont = limCont + 1
        limVert = limVert + npt_in_obj(iobj)
      endif
    enddo
    return
  endif
  !
  ! Load planar data into x/y arrays
  !
  indCur = 0
  do iobj = 1, max_mod_obj
    if (npt_in_obj(iobj) >= 3) then
      numConts = numConts + 1
      if (numConts > limCont) call exitError('Too many contours in model')
      numVerts(numConts) = npt_in_obj(iobj)
      if (indCur + numVerts(numConts) > limVert) call exitError( &
          'Too many points in contours')
      do ip = 1, numVerts(numConts)
        ipt = abs(object(ibase_obj(iobj) + ip))
        xVerts(ip + indCur) = p_coord(1, ipt)
        yVerts(ip + indCur) = p_coord(indY, ipt)
      enddo
      contourZ(numConts) = p_coord(indZ, ipt)
      indVert(numConts) = indCur + 1
      indCur = indCur + numVerts(numConts)
    endif
  enddo
  if (progname .ne. 'FILLTOMO') &
      print *,numConts, ' contours available for deciding which patches to analyze'
  return
end subroutine get_region_contours


! getContourArraySizes is a convenience function for getting the array sizes
!
subroutine getContourArraySizes(modelFile, imUnit, limCont, limVert)
  implicit none
  character*(*) modelFile
  integer*4 numConts, ifFlip, limCont, limVert, imUnit, numVerts, indVert
  real*4 dummy
  limCont = -1
  limVert = -1
  call get_region_contours(modelFile, 'DUMMY', dummy, dummy, numVerts, &
      indVert, dummy,  numConts, ifFlip, limCont, limVert, imUnit)
  limVert = limVert + 10
  limCont = limCont + 10
  return
end subroutine getContourArraySizes


! checkBoundaryConts checks a patch center against the boundary contours by finding the
! contour at the nearest Z level and testing whether the center is inside the contour
!
subroutine checkBoundaryConts(cenX, cenY, cenZ, ifUse, numConts, numVerts, &
    xVerts, yVerts, contourZ, indVertStart)
  implicit none
  integer*4 ifUse, numConts, numVerts(*), indVertStart(*)
  real*4 xVerts(*), yVerts(*), contourZ(*), cenX, cenY, cenZ
  integer*4 icont, icontMin, indv
  real*4 dz, dzMin
  logical inside

  ifUse = 0
  !
  ! find nearest contour in Z and see if patch is inside it
  !
  dzMin = 100000.
  do icont = 1, numConts
    dz = abs(cenZ - contourZ(icont))
    if (dz < dzMin) then
      dzMin = dz
      icontMin = icont
    endif
  enddo
  indv = indVertStart(icontMin)
  if (inside(xVerts(indv), yVerts(indv), numVerts(icontMin), cenX, cenY)) ifUse = 1
  return
end subroutine checkBoundaryConts


! summarizeDrops outputs a summary of residuals dropped as outliers,
! breaking them into 10 bins from the minimum to either the maximum or
! 5 SDs above the mean, whichever is less.  DROPSUM has the residuals,
! NLISTD is the number of values, and meanText is text to indicate if
! they are mean residuals or not.
!
subroutine summarizeDrops(dropSum, numListDrop, meanText)
  implicit none
  character*(*) meanText
  real*4 dropSum(*)
  integer*4 numListDrop, i, j, inBin
  real*4 dropMin, dropMax, dropAvg, dropSD, dropSEM, binLow, binHigh, dropBin
  !
  dropMin = 1.e10
  dropMax = 0.
  do i = 1, numListDrop
    dropMin = min(dropMin, dropSum(i))
    dropMax = max(dropMax, dropSum(i))
  enddo
  call avgsd(dropSum, numListDrop, dropAvg, dropSD, dropSEM)
  dropBin = (min(dropMax, dropAvg + 5 * dropSD) - dropMin) / 10.
  !
  do j = 1, 10
    binLow = dropMin + (j - 1) * dropBin
    binHigh = binLow + dropBin
    if (j == 10) binHigh = dropMax + 0.001
    inBin = 0
    do i = 1, numListDrop
      if (dropSum(i) >= binLow .and. dropSum(i) < binHigh) &
          inBin = inBin + 1
    enddo
    if (inBin > 0) write(*,109) inBin, meanText, binLow, binHigh
109 format(i8,' with ',a,'residuals in',f10.2,' -',f10.2)
  enddo
  return
end subroutine summarizeDrops


! getExtraSelections gets specifications for selecting patches based on the values in
! extra columns and checks them for validity
!
subroutine getExtraSelections(icolSelect, isignSelect, numColSelect, numSelectCrit, &
    IDextra, maxExtra, selectCrit, limSelect, limCrit)
  implicit none
  integer*4 icolSelect(*), isignSelect(*), numColSelect, numSelectCrit, IDextra(*)
  integer*4 maxExtra, limSelect, limCrit
  real*4 selectCrit(limCrit, limSelect)
  integer*4 ierr, ix, i
  integer*4 PipNumberOfEntries, PipGetTwoIntegers, PipGetFloatArray

  ierr = PipNumberOfEntries('ExtraValueSelection', numColSelect)
  ierr = PipNumberOfEntries('SelectionCriteria', ix)
  if (ix .ne. numColSelect) call exitError( &
      'There must be one -select entry for each -extra entry')
  if (numColSelect > 0 .and. maxExtra == 0) call exitError( &
      'There are no extra value columns in the patch file')
  if (numColSelect > limSelect) call exitError('Too many selections for arrays')
  do i = 1, numColSelect
    ierr = PipGetTwoIntegers('ExtraValueSelection', icolSelect(i), isignSelect(i))
    ix = 0
    ierr = PipGetFloatArray('SelectionCriteria', selectCrit(1, i), ix, limCrit)
    if (i > 1 .and. ix .ne. numSelectCrit) call exitError( &
        'Every entry of -select must have the same number of criteria')
    numSelectCrit = ix
    if (abs(isignSelect(i)) .ne. 1) call exitError( &
        'The second value on the -select entry must be 1 or -1')
    if (icolSelect(i) < 0) then
      if (icolSelect(i) > maxExtra) call exitError('There is no extra value column'// &
          ' corresponding to the column number entered with -extra')
      icolSelect(i) = -icolSelect(i)
    else
      if (icolSelect(i) == 0) call exitError('0 is not a value column number or id #')
      ierr = 1
      do ix = 1, maxExtra
        if (icolSelect(i) == IDextra(ix)) then
          ierr = 0
          icolSelect(i) = ix
          exit
        endif
      enddo
      if (ierr > 0) call exitError( &
          'There is no extra value column with the column id # entered with -extra')
    endif
  enddo
  return
end subroutine getExtraSelections
