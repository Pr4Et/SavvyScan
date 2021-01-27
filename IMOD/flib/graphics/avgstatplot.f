*       * * * * AVGSTATPLOT * * * * *
c       
c       AVGSTATPLOT is an interactive program for displaying and plotting the
c       output of the program IMAVGSTAT.  See man page for details.
c       

      call plax_initialize('avgstatplot')
      call exit(0)
      end

      subroutine realgraphicsmain(comFileArg)
      character*(*) comFileArg
      parameter (limarea=5000,limset=400,limdat=100000)
      real*4 avg(limarea,limset),sd(limarea,limset)
     &    ,sem(limarea,limset),xx(limdat),yy(limdat),semadd(limdat),
     &    sclfac(limset),scladd(limset),xxft(limarea),yyft(limarea)
      integer*4 nsumarea(limarea),indregion(limarea),ngx(limdat),
     &    nsymb(500),kset(500),ktype(500),kreg(500),nregress(limset),
     &    nreguse(limset),ireguse(limarea,limset),isrescl(limset),
     &    npixarea(limarea),normreg(limset),normarea(limset),
     &    nsampl(limset)
      character*80 statname
      character*23 defsetxt,fmt
      character*8 outregarea(37)
c       real*8 g01caf
c       
5     write(*,'(1x,a,$)')'Name of statistics file: '
      read(*,'(a)')statname
      close(1)
      call dopen(1,statname,'ro','f')
c       
      write(*,'(1x,a,$)')
     &    '0 for plots on parallax, 1 for terminal only: '
      read(*,*)iffil
      call scrnOpen(iffil)
c       
      read(1,*)nregion
      read(1,*)(nsumarea(i),i=1,nregion)
c       
      ntotarea=0
      do i=1,nregion
        indregion(i)=ntotarea+1
        ntotarea=ntotarea+nsumarea(i)
      enddo
c       
      read(1,*)ntotin
      if(ntotin.ne.ntotarea)print *,
     &    'Warning - number of areas does not add up correctly'
      read(1,*)(npixarea(i),i=1,ntotin)
      read(1,*)nsets
c       
      write(*,101)nsets,nregion,ntotarea,(nsumarea(i),i=1,nregion)
101   format(i5,' data sets',/,i5,' summing regions with a total of',
     &    i5,' summing areas',/,' Number of areas in these regions:',
     &    /,(20i4))
c       
      do iset=1,nsets
        read(1,*)nsampl(iset)
        do iarea=1,ntotarea
          read(1,*)idum,avg(iarea,iset),sd(iarea,iset)
     &        ,sem(iarea,iset)
        enddo
      enddo
c       
10    print *,'Enter list of numbers of sets to plot (ranges ok)'
      call rdlist(5,kset,nsetplot)
      ntytot=0
      do while(ntytot.lt.nsetplot)
        write(*,'(a,i3,a)')' Enter list of symbol types for',
     &      nsetplot-ntytot, ' sets (ranges ok)'
        call rdlist(5,ktype(ntytot+1),ntyadd)
        ntytot=ntytot+ntyadd
      enddo
      nrescale=nsetplot
      defsetxt='all sets)'
      do i=1,nsetplot
        isrescl(i)=kset(i)
      enddo
c       
15    print *,'Enter list of numbers of regions to plot (ranges ok)'
      call rdlist(5,kreg,nregplot)
c       
      ifsubset=0
      if(nregplot.eq.1.and.nsumarea(kreg(1)).gt.20)then
        ifsubset=1
        istrsub=1
        iendsub=nsumarea(kreg(1))
        write(*,'(1x,a,$)')
     &      'Beginning and ending areas to plot, or / for all: '
        read(5,*)istrsub,iendsub
        istrsub=max(1,istrsub)
        iendsub=min(nsumarea(kreg(1)),iendsub)
        avgspace=0
        write(*,'(1x,a,/,a,$)')'Enter interval to average over'//
     &      ' (fractions OK), # of areas to roll average right,',
     &      '  and # of areas to replicate symmetrically'//
     &      ' (or / for no averaging): '
        read(5,*)avgspace,nrollavg,nreplic
      endif
c       
20    write(*,'(1x,a,/,a,$)')'For error bars, enter a small positive'
     &    //' # for that # of SEM, a negative # for',
     &    ' that # of SD, or a large positive # for '//
     &    'confidence limits at that % level: '
      read(*,*)facsem
c       
      write(*,'(1x,a,$)')'0 to use means, 1 to use integrals: '
      read(*,*)integral
      print *,'Enter list of numbers of sets to rescale'
      print *,'     (ranges OK, Return for none, / for ',defsetxt
      call rdlist(5,isrescl,nrescale)
      defsetxt='same sets as last time)'
      if(nrescale.gt.0.and.(nrescale.gt.1.or.isrescl(1).ne.0))then
        write(*,'(1x,a,$)')'0 to specify each independently, 1 to '
     &      //'use same specification for all: '
        read(*,*)ifsamescl
        nspecify=nrescale
        if(ifsamescl.ne.0)nspecify=1
        do ires=1,nspecify
          write(*,'(1x,a,i3,a,/,a,/,a,$)')'For set #',isrescl(ires),
     &        ', enter 0 to specify scaling, 999 to divide by one'//
     &        ' area,','  or # of other set to compare with',
     &        '     (+ set # for regression, - set # '//
     &        'to shift to same mean): '
          read(*,*)nregress(ires)
          if(nregress(ires).eq.0)then
            write(*,'(1x,a,$)')'Factors to multiply by, then add: '
            read(*,*)sclfac(ires),scladd(ires)
          elseif(nregress(ires).eq.999)then
            write(*,'(1x,a,$)')'Region #, and # of area within'//
     &          ' region to divide by: '
            read(*,*)normreg(ires),normarea(ires)
          else
            print *,'Enter list of regions to use for comparing sets'//
     &          ' (ranges OK)'
            call rdlist(5,ireguse(1,ires),nreguse(ires))
          endif
        enddo
      else
        nrescale=0
      endif
c       
      fracofs=0
      if(nsetplot.gt.1)then
        write(*,'(1x,a,$)')
     &      'Fraction to offset sets from each other in X: '
        read(*,*)fracofs
      endif
c       
      nx=0
      ngrps=0
      ymax=0.
      do iset=1,nsetplot
        xar=(iset-nsetplot/2)*fracofs
        jset=kset(iset)
        ires=0
        sclf=1.
        scla=0.
        do i=1,nrescale
          if(jset.eq.isrescl(i))ires=i
        enddo
        if(ires.gt.0)then
          if(ifsamescl.ne.0)ires=1
          if(nregress(ires).eq.0)then
            sclf=sclfac(ires)
            scla=scladd(ires)
          elseif(nregress(ires).eq.999)then
            scla=0.
            indar=indregion(normreg(ires))+normarea(ires)-1
            sclf=1./avg(indar,jset)
            if(integral.ne.0)sclf=sclf/npixarea(indar)
          else
            npnts=0
            do ireg=1,nreguse(ires)
              jreg=ireguse(ireg,ires)
              indstr=indregion(jreg)
              if(ifsubset.eq.0)then
                indend=indstr+nsumarea(jreg)-1
              else
                indend=indstr+iendsub-1
                indstr=indstr+istrsub-1
              endif
              do indar=indstr,indend
                npnts=npnts+1
                yyft(npnts)=avg(indar,abs(nregress(ires)))
                xxft(npnts)=avg(indar,jset)
              enddo
            enddo
            if(npnts.gt.0)then
              if(nregress(ires).gt.0)then
                call lsfit(xxft,yyft,npnts,sclf,scla,rho)
              else
                sclf=1.
                call avgsd(xxft,npnts,xxavg,xxsd,xxsem)
                call avgsd(yyft,npnts,yyavg,yysd,yysem)
                scla=yyavg-xxavg
              endif
            endif
            write(*,103)jset,npnts,rho,sclf,scla
103         format(' Set #',i3,', n=',i4,', r=',f6.3,', multiply by',
     &          f9.5,' and add',f7.2)
          endif
        endif
c         
        do ireg=1,nregplot
          ngrps=ngrps+1
          nsymb(ngrps)=ktype(iset)
          jreg=kreg(ireg)
          indstr=indregion(jreg)
          if(ifsubset.eq.0)then
            indend=indstr+nsumarea(jreg)-1
          else
            indend=indstr+iendsub-1
            indstr=indstr+istrsub-1
          endif
          do indar=indstr,indend
            xar=xar+1.
            nx=nx+1
            if(facsem.gt.30.)then
              ifail=0
c               tcrit=g01caf(dble((1.+0.01*facsem)/2.),nsampl(jset)-1,ifail)
              tcrit=tvalue((1.+0.01*facsem)/2.,nsampl(jset)-1)
              semadd(nx)=sclf*tcrit*sem(indar,jset)
            elseif(facsem.ge.0)then
              semadd(nx)=sclf*facsem*sem(indar,jset)
            else
              semadd(nx)=-sclf*facsem*sd(indar,jset)
            endif
            xx(nx)=xar
            yy(nx)=sclf*avg(indar,jset)+scla
            if(integral.ne.0)then
              yy(nx)=yy(nx)*npixarea(indar)
              semadd(nx)=semadd(nx)*npixarea(indar)
            endif
            ngx(nx)=ngrps
            ymax=max(ymax,abs(yy(nx)))
          enddo
          if(avgspace.gt.0.)then
            navgspace=nint(avgspace)
            ntotval=indend+1-indstr
            ipt=nx+1-ntotval
            do iavg=1,navgspace
              navlim=1+ntotval/avgspace
              navg=0
              do itmp=1,navlim
                indyy=ipt+nint((itmp-1)*avgspace)
                if(indyy.le.nx)then
                  yyft(itmp)=yy(indyy)
                  navg=navg+1
                endif
              enddo
              call avgsd(yyft,navg,spacavg,spacsd,spacsem)
              yy(ipt)=spacavg
              if(facsem.ge.0)semadd(ipt)=facsem*spacsem
              if(facsem.lt.0)semadd(ipt)=-facsem*spacsd
              ipt=ipt+1
            enddo
            do iavg=1,navgspace
              ipt=nx+iavg-ntotval
              xxft(iavg)=semadd(ipt)
              yyft(iavg)=yy(ipt)
            enddo
            ifrom=navgspace-nreplic/2-nrollavg
            if(ifrom.lt.1)ifrom=ifrom+navgspace
            if(ifrom.lt.1)ifrom=ifrom+navgspace
            if(ifrom.gt.navgspace)ifrom=ifrom-navgspace
            do iavg=1,navgspace+nreplic
              ipt=nx+iavg-ntotval
              yy(ipt)=yyft(ifrom)
              semadd(ipt)=xxft(ifrom)
              ifrom=mod(ifrom,navgspace)+1
            enddo
            nx=nx+navgspace+nreplic-ntotval
          endif
        enddo
      enddo
      call errplt(xx,yy,ngx,nx,nsymb,ngrps,semadd,0,0)
30    write(*,'(1x,a,/,a,$)')'Respecify error bars (1), regions (2)'//
     &    ', data sets (3), or input file (4),',
     &    ' Plot metacode on screen (5) or printer (6),'//
     &    ' type values (7), or exit (8): '
      read(*,*)iopt
      if(iopt.lt.1)go to 30
      go to (20,15,10,5,35,35,40,60),iopt
      go to 30
c       
35    call pltout(6-iopt)
      go to 30
40    print *,'Enter file name to store values in, or Return',
     &    ' to type values on screen'
      read(*,'(a)')statname
      if(statname.eq.' ')then
        iout=6
      else
        iout=9
        ifappend=0
        close(9)
        open(9,file=statname,err=40,status='unknown')
42      read(9,'(a4)',end=44)statname
        if(ifappend.eq.0)print *,'Appending to existing file...'
        ifappend=1
        go to 42
      endif
44    write(iout,'(a,20x,a)')' Set','Region-area'
      nareaout=0
      do ireg=1,nregplot
        do iarea=1,nsumarea(jreg)
          nareaout=nareaout+1
          write(outregarea(nareaout),'(i5,''-'',i2)')kreg(ireg),iarea
        enddo
      enddo
      write(iout,'(5x,9a8)')(outregarea(i),i=1,nareaout)
      write(iout,*)
      ipow=max(0.,min(4.,alog10(10.*ymax)))
      fmt='(i4,(t7,9f8.'//char(52-ipow)//'))'
      nperset=nx/nsetplot
      ibase=0
      do iset=1,nsetplot
        write(iout,fmt)kset(iset),(yy(i),i=ibase+1,ibase+nperset)
        ibase=ibase+nperset
      enddo
      write(iout,'(/)')
      go to 30
60    call scrnClose()
      call psClose()
      end
