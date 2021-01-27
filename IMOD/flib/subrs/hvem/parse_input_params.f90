! PipParseInput provides a simple front end to the "pip" package
! for parsing input parameters
!
! "options" is an array of character strings, one for each option,
! listing the short name, long name, type code, and help string,
! separated by commas.  In this case "separator" should be a space.
! Alternatively, options can be an array consisting of one long
! character string, with each option separated from the next by
! the character in "separator"
! "numOptions" is the number of options in the array or string
! The routine returns numOptArg, the number of option arguments
! on the command line; and numNonOptArg, the number of non-option
! arguments.
!
! $Id$
!

integer*4 function PipParseInput(options, numOptions, separator, numOptArg, numNonOptArg)
  implicit none
  character*(*) options(*)
  character separator
  integer*4 numOptArg, numOptions, numNonOptArg
  integer*4 PipInitialize, PipAddOption, PipParseEntries
  integer*4 i, j, indStr, indEnd, lenAll
  !
  ! initialize then pass the options one by one
  !
  PipParseInput = PipInitialize(numOptions)
  if (PipParseInput .ne. 0) return
  if (separator == ' ') then
    do i = 1, numOptions
      PipParseInput = PipAddOption(options(i))
      if (PipParseInput .ne. 0) return
    enddo
  else
    !
    ! if options are all in one string with a separator
    !
    indStr = 1
    lenAll = len_trim(options(1))
    do i = 1, numOptions
      j = indStr
      indEnd = 0
      do while (j <= lenAll .and. indEnd == 0)
        if (options(1) (j:j) == separator) indEnd = j - 1
        j = j + 1
      enddo
      if (j > lenAll) indEnd = lenAll
      if (indStr > indEnd) then
        call PipSetError('Too few options in string')
        PipParseInput = -1
        return
      endif
      PipParseInput = PipAddOption(options(1) (indStr:indEnd))
      if (PipParseInput .ne. 0) return
      indStr = j
    enddo
  endif

  PipParseInput = PipParseEntries(numOptArg, numNonOptArg)
  return
end function PipParseInput


! Function to parse the entries to the program after options have
! been defined
!
integer*4 function PipParseEntries(numOptArg, numNonOptArg)
  implicit none
  integer bufferSize
  parameter (bufferSize = 1024)
  integer*4 numOptArg, numNonOptArg
  integer*4 iargc, i
  integer*4 PipNextArg, PipReadStdinIfSet
  character*(bufferSize) string
  !
  ! pass the arguments in one by one
  !
  do i = 1, iargc()
    call getarg(i, string)
    if (len_trim(string) == bufferSize) then
      call PipSetError('Input argument too long for buffer in PipParseEntries')
      PipParseEntries = -1
      return
    endif
    PipParseEntries = PipNextArg(string)
    if (PipParseEntries < 0) return
    if (PipParseEntries > 0 .and. i == iargc()) then
      call PipSetError('A value was expected but not found for'// &
          ' the last option on the command line')
      PipParseEntries = -1
      return
    endif
  enddo
  !
  ! Or read stdin if the flag is set to do it when no arguments
  !
  if (iargc() == 0) then
    PipParseEntries = PipReadStdinIfSet()
    if (PipParseEntries .ne. 0) return
  endif
  !
  ! get numbers to return
  !
  call PipNumberOfArgs(numOptArg, numNonOptArg)
  call PipPrintEntries()
  PipParseEntries = 0
  return
end function PipParseEntries


! PipGetLogical: Function to get a Fortran logical directly
!
integer*4 function PipGetLogical(option, value)
  implicit none
  character*(*) option
  logical value
  integer*4 intval, ierr
  integer*4 PipGetBoolean

  intval = 0
  PipGetLogical = 0
  ierr = PipGetBoolean(option, intval)
  if (ierr .ne. 0) then
    PipGetLogical = ierr
  else
    value = intval .ne. 0
  endif
  return
end function PipGetLogical


! PipReadOrParseOptions will first try to read an options file,
! then fallback to a list of options supplied to it
!
subroutine PipReadOrParseOptions(options, numOptions, progName, exitString, interactive, &
    minArgs, numInFiles, numOutFiles, numOptArg, numNonOptArg)
  implicit none
  character*(*) options(*)
  character*(*) progName
  character*(*) exitString
  logical interactive
  integer*4 minArgs, numInFiles, numOutFiles
  integer*4 numNonOptArg, numOptArg, numOptions, ierr
  character*240 errString
  integer*4 PipGetError, PipParseInput, PipReadOptionFile, PipParseEntries
  integer*4 PipGetBoolean, PipPrintHelp
  !
  ! First try to read autodoc file
  !
  call PipAllowCommaDefaults(1)
  ierr = PipReadOptionFile(progName, 1, 0)
  call PipExitOnError(0, exitString)
  !
  ! If that is OK, go parse the entries;
  ! otherwise print error message and use fallback option list
  !
  if (ierr == 0) then
    ierr = PipParseEntries(numOptArg, numNonOptArg)
  else
    ierr = PipGetError(errString)
    write(*, '(a, a)') 'PIP WARNING: ', errString
    print *,'Using fallback options in main program'
    ierr = PipParseInput(options, numOptions, '@', numOptArg, numNonOptArg)
    if (ierr == 0) call PipReadProgDefaults(progName)
  endif
  !
  ! Process help input
  !
  if (interactive .and. numOptArg + numNonOptArg == 0) return
  if (numOptArg + numNonOptArg < minArgs .or. PipGetBoolean('help', ierr) == 0) then
    ierr = PipPrintHelp(progName, 0, numInFiles, numOutFiles)
    call exit(0)
  endif
  return
end subroutine PipReadOrParseOptions


! PipGetInOutFile gets one of the input or output files specified by
! "option", or found at the "nonOptArgNo' non-option argument position
! The name is returned in "filename"
! If there no such entry, it returns with an error
! If PIP input is not being used, it used "prompt" to ask for the name
!
integer*4 function PipGetInOutFile(option, nonOptArgNo, prompt, filename)
  implicit none
  character*(*) option
  character*(*) prompt
  character*(*) filename
  integer*4 nonOptArgNo, numOptArg, numNonOptArg
  integer*4 PipGetString, PipGetNonOptionArg
  !
  PipGetInOutFile = 0
  call PipNumberOfArgs(numOptArg, numNonOptArg)
  !
  ! if there is PIP input, first look for explicit option by the name
  ! then get the given non-option argument if there are enough
  !
  if (numOptArg + numNonOptArg > 0) then
    PipGetInOutFile = PipGetString(option, filename)
    if (PipGetInOutFile .ne. 0) then
      if (numNonOptArg < nonOptArgNo) return
      PipGetInOutFile = PipGetNonOptionArg(nonOptArgNo, filename)
    endif
  else
    !
    ! Otherwise get interactive input with the prompt
    !
    write(*,'(1x,a,a,$)') prompt, ': '
    read(5, '(a)') filename
  endif
  return
end function PipGetInOutFile


! Controls whether to exit on error
!
subroutine PipExitOnError(ifUseStderr, message)
  implicit none
  character*(*) message
  integer*4 ifUseStderr
  call PipExitOnErrorFW(ifUseStderr, message)
  call setExitPrefix(message)
  return
end subroutine PipExitOnError


! Exits with error status after issuing the given message, with the
! prefix set by calling setExitPrefix
!
subroutine exitError(message)
  implicit none
  character*(*) message
  character*32 prefix
  common / exitprefix / prefix
  write(*,'(/,a,a,a)') trim(prefix), ' ', trim(message)
  call pipexit(1)
end subroutine exitError


! Sets prefix for exiting with error
!
subroutine setExitPrefix(message)
  character*(*) message
  character*32 prefix
  common / exitprefix / prefix
  prefix = message
  return
end subroutine setExitPrefix


! Exits with allocation error message if ierr not zero
!
subroutine memoryError(ierr, message)
  implicit none
  integer*4 ierr
  character*(*) message
  if (ierr .ne. 0) call exitError('Failure to allocate '//message)
  return
end subroutine memoryError

! Exits with upper case allocation error message if ierr not zero
!
subroutine memoryErrorUC(ierr, message)
  implicit none
  integer*4 ierr
  character*(*) message
  if (ierr .ne. 0) call exitError('FAILURE TO ALLOCATE '//message)
  return
end subroutine memoryErrorUC


! Exits with message if current Adoc cannot be set
!
subroutine setCurrentAdocOrExit(indAdoc, message)
  implicit none
  integer*4 indAdoc, AdocSetCurrent
  character*(*) message
  if (AdocSetCurrent(indAdoc) .ne. 0) call exitError('Selecting '//trim(message) // &
      ' autodoc as current one')
  return
end subroutine setCurrentAdocOrExit
