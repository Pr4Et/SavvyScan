c       CONVERTMOD  will read in a model file (WIMP ASCII or IMOD) and
c       output a WIMP-style ASCII model file

      include 'model.inc'

      character*320 oldfile,newfile
      character*6 binasc1,binasc2
      logical readw_or_imod
      integer*4 ierr
      integer*4 imodWriteAsWimp

c       get specifications
      call getinout(2, oldfile, newfile)
      if(.not.readw_or_imod(oldfile))then
        print *,'Error reading mode file'
        call exit(1)
      endif
c       
      close(20)
      ierr = imodWriteAsWimp(newfile)
      if (ierr .eq. -5) then
        open(20,file=newfile,status='new', form='formatted')
        call store_mod(newfile)
        close(20)
      elseif (ierr .ne. 0) then
        print *,'Error', ierr,' writing to WIMP model file'
        call exit(1)
      endif
      call exit(0)
      end

