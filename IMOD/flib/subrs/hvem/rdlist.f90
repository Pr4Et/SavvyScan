! RDLIST reads a line from unit IUNit specifying ranges of values and
! returns a list of all of the values (NLIST values in array LIST) .
! E.g. 1-5, 7, 9, 11, 15-20.  Numbers separated by dashed are replaced by
! all of the numbers in the range.  NUmbers need not be in any order,
! and backward ranges (10-5) are handled.  Any characters besides
! digits are valid separators. A / at the beginning of the line will
! return an unmodified list. Negative numbers can be entered provided
! that the minus sign immediately precedes the number.  E.g.: -3 - -1
! or -3--1 will give -3, -2, -1; -3, -1, 1 or -3, -1, 1 will give -3, -1, 1.
!
! $Id$
!
!
! Basic (unsafe) function provides unit number, returns list and number in list
!
subroutine rdlist(iunit, list, numInList)
  implicit none
  integer*4 list(*), iunit, numInList
  call rdlist2(iunit, list, numInList, 0)
  return
end subroutine rdlist


! Safer function provides limlist to specify array size, which can be negative to
! return with an error in limlist: 1 for too many, 3 for allocation error,
! 4 for bad character in entry
!
subroutine rdlist2(iunit, list, numInList, limList)
  implicit none
  integer*4 list(*), iunit, numInList, limList
  character*10240 line
  read(iunit, '(a)') line
  call parselist2(line, list, numInList, limList)
  return
end subroutine rdlist2


! readBigList takes an allocated array of whatever type and treats it as character
! variable for reading in a potentially very big list.  It is a workaround because
! gfortran did not support allocated character variables until 4.8 or so.
!
subroutine readBigList(inUnit, inList, nList, limList, bigString, length)
  implicit none
  integer*4 inUnit, inList, nList, limList,length
  character(len=length) bigString
  read(inUnit, '(a)') bigString
  call parselist2(bigString, inList, nList, limList)
end subroutine readBigList


! Unsafe version: parse the list in line, return in list and # in nlist
!
subroutine parselist(line, list, numInList)
  implicit none
  character*(*) line
  integer*4 list(*), numInList
  call parselist2(line, list, numInList, 0)
  return
end subroutine parselist

! safer function provides limList to specify array size, which can be negative to
! return with an error in limlIst
subroutine parselist2(line, list, numInList, limList)
  implicit none
  character*(*) line
  integer*4 list(*), numInList, limList, parselistfw, ierr
  ierr = parselistfw(line, list, numInList, abs(limList))
  if (ierr == 0) then
    if (limList < 0) limList = 0
    return
  endif
  if (ierr < 0) then
    write(*,'(/,a)') 'ERROR: PARSELIST - TOO MANY LIST VALUES FOR ARRAY'
  elseif (ierr == 1) then
    write(*,'(/,a)') 'ERROR: PARSELIST - FAILED TO ALLOCATE MEMORY FOR LIST'
  else
    write(*,'(/,a)') 'ERROR: PARSELIST - BAD CHARACTER IN ENTRY'
  endif
  if (limList > 0) call exit(1)
  limList = ierr + 2
  return
end subroutine parselist2
