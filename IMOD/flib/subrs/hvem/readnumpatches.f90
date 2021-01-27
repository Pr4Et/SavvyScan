! !
! Reads in the first line of a patch file opened on the unit [iunit], and returns the
! first value on the line in [numPatch].  It looks for further numeric values on the
! line, up to [limVals] of them, and returns them in [values] and the number found in
! [numVals].  The return value is 1 for error reading the file, or 2 if the first item
! on the line is not an integer.
! !
integer*4 function readNumPatches(iunit, numPatch, numVals, idValues, limVals)
  implicit none
  integer*4 iunit, numPatch, numVals, limVals, idValues(*)
  real*4 xnum(2 * limVals + 9)
  integer*4 numeric(2 * limVals + 9), ierr, ind, numFields
  character*1024 line
  readNumPatches = 1
  read(iunit, '(a)', iostat = ierr) line
  if (ierr .ne. 0) return
  call frefor2(line, xnum, numeric, numFields, 2 * limVals + 9)
  readNumPatches = 2
  if (line(1:1) == '/' .or. numFields < 1 .or. numeric(1) == 0) return
  numPatch = nint(xnum(1))
  if (abs(numPatch - xnum(1)) > 1.e3) return
  numVals = 0
  do ind = 2, numFields
    if (numeric(ind) > 0 .and. numVals < limVals) then
      numVals = numVals + 1
      idValues(numVals) = nint(xnum(ind))
    endif
  enddo
  readNumPatches = 0
  return
end function readNumPatches
