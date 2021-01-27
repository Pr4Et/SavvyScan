! !
! Find_surfaces will analyze a set of [numRealPt] points, with coordinates
! in [xyz] (dimensioned to (3,*)), correlate the Z with the X and Y
! coordinates of those points, and determine the angles that the points
! would have to be tilted first around the X axis then around the Y axis
! in order for them to lie parallel to the X - Y plane.   If [numSurface] is
! 1, it will simply fit a plane to all of the points and base its
! estimates on that fit.  If [numSurface] is 2, it will attempt to divide
! the points into two groups occupying two surfaces, then provide
! separate estimates of the new tilt angle based on the slope of the
! plane through either set of points, or on the average slope.  The
! array [igroup] is returned with a value for each point of 1 if on lower
! surface, 2 if on upper, but only if 2 - surface analysis is done.  All
! outputs are printed for units 6 through [iunit], so set [iunit] to 6
! for output to standard out only, to 7 for output to a file on unit 7,
! or < 6 for no output.  If [ifComp] is non - zero, it assumes data were
! obtained at a single tilt angle [tiltMax] and at zero tilt and will
! estimate the true tilt angle of the section, returned in [tiltNew].
! [tiltAdd] should be set to an existing change in tilt angles so that
! the total tilt angle change can be output.  The values in [znew] and [znewInput], the
! actual and input amounts to shift the tilt axis in Z, and in [imageBinned] allow it to
! report on the unbinned thickness between fiducials and shift needed to
! center them.
! !
! $Id$
!
subroutine find_surfaces(xyz, numRealPt, numSurface, tiltMax, &
    iunit2, tiltNew, igroup, ifComp, tiltAdd, znew, znewInput, imageBinned)
  implicit none
  real*4 xyz(3,*), tiltMax, tiltNew, tiltAdd, bintcpMinus, bintcp, znew, znewInput
  integer*4 numRealPt, numSurface, iunit2, ifComp, imageBinned
  real*4 xmat(numRealPt, 4)
  integer*4 igroup(*)
  integer*4 i, iun, numPntMinus, numPntPlus
  real*4 aSlope, bSlope, alpha, slope, resid, truePlus
  real*4 residMinus, residPlus, thick
  real*4 shiftInc, shiftTot, botExtreme, topExtreme
  integer*4 surfaceSort
  !
  ! first fit a line to all of the points to get starting angle
  !
  do i = 1, numRealPt
    xmat(i, 1) = xyz(1, i)
    xmat(i, 2) = xyz(2, i)
    xmat(i, 3) = xyz(3, i)
    igroup(i) = 0
  enddo
  call lsfit2Resid(xmat(1, 1), xmat(1, 2), xmat(1, 3), numRealPt, aSlope, bSlope, &
      bintcp, alpha, slope, resid, botExtreme, topExtreme)
  !
  ! Batchruntomo is looking for '# of points' and the number after '=' and expects
  ! either one number for all, or 3 numbers for all, bottom, top
  do iun = 6, iunit2
    write(iun, '(/,a,//,a,//,a,/,a,i16,/,a,f11.2)')' SURFACE ANALYSIS:', &
        ' The following parameters are appropriate if fiducials' &
        //' are NOT on two surfaces:', &
        ' Fit of one plane to all fiducials:', &
        ' # of points = ', numRealPt, &
        ' Mean residual =    ', resid
    if (numRealPt < 3) then
      write(iun, 99) ' Too few points to estimate tilt angle change or X-axis tilt'
99    format(a)
      truePlus = 0.
    else
      write(iun, '(a, f8.4)') ' Adjusted slope =      ', slope
      if (numRealPt < 4) then
        write(iun, 99) ' Too few points to estimate X-axis tilt'
      else
        write(iun, '(a, f10.2)') ' X axis tilt needed =', -alpha
      endif
      call calcTiltNew(slope, tiltMax, iunit2, truePlus, ifComp, tiltAdd)
    endif
  enddo
  !
  ! if there are supposed to be 2 surfaces, call the surfaceSort routine to get them
  ! sorted out.  Then fit a pair of parallel planes to both surfaces
  !
  if (numSurface > 1 .and. numRealPt > 3) then
    if (surfaceSort(xyz, numRealPt, 0, igroup) .ne. 0) &
      call exitError('Allocating memory in surfaceSort')

    call twoSurfaceFits(xyz, igroup, numRealPt, xmat, numPntMinus, bintcpMinus, &
        residMinus, numPntPlus, aSlope, bSlope, bintcp, residPlus, alpha, slope, resid, &
        botExtreme, topExtreme)

    do iun = 6, iunit2
      write(iun, '(a,/,a)') ' The following '// &
          'parameters are appropriate if fiducials ARE on two surfaces,', &
          ' based on fit of two parallel planes to fiducials sorted onto two surfaces:'
      write(iun, 101) 'bottom', numPntMinus, residMinus, bintcpMinus
101   format(/, ' On ',a,' surface:',/, &
          ' # of points = ',i12,/, &
          ' Mean residual =',f11.2,/, &
          ' Z axis intercept =',f8.1)
      write(iun, 101) 'top', numPntPlus, residPlus, bintcp
      thick = bintcp - bintcpMinus
      write(iun, 102)resid, thick, slope
102   format(/,' Overall mean residual =',f15.2,/, &
          ' Thickness at Z intercepts =',f11.1,/, &
          ' Adjusted slope =',f22.4)
      if (numRealPt > 4) then
        write(iun, '(a, f18.2)') ' X axis tilt needed =', -alpha
      else
        write(iun, 99) ' Too few points to estimate X-axis tilt'
      endif
    enddo
    call calcTiltNew(slope, tiltMax, iunit2, tiltNew, ifComp, tiltAdd)
  else
    numPntPlus = numRealPt
    tiltNew = truePlus
    if (numSurface > 1) write(*,'(/,a,/)') &
        ' There are too few points to analyze for distribution on two surfaces'
  endif
  !
  ! Get the unbinned thickness and shifts needed to center the gold
  ! The direction is opposite to expected because the positive Z points
  ! come out on the bottom of the tomogram, presumably due to rotation
  thick = imageBinned * (topExtreme - botExtreme)
  shiftTot = imageBinned * (topExtreme + botExtreme) / 2.
  shiftInc = shiftTot - imageBinned * znew
  shiftTot = shiftInc + imageBinned * znewInput
  do iun = 6, iunit2
    write(iun, 103) thick, shiftInc, shiftTot
103 format(' Unbinned thickness needed to contain centers of all ', &
        'fiducials =', f13.0,/, &
        ' Incremental unbinned shift needed to center range of fi', &
        'ducials in Z =',f8.1,/, &
        ' Total unbinned shift needed to center range of fiducial', &
        's in Z =',f14.1)
  enddo
  return
end subroutine find_surfaces


subroutine twoSurfaceFits(xyz, igroup, numRealPt, xmat, numPntMinus, &
    bintcpMinus, residMinus, numPntPlus, aSlope, &
    bSlope, bintcp, residPlus, alpha, slope, resid, botExtreme, topExtreme)
  implicit none
  real*4 xyz(3,*), bintcpMinus
  real*4 xmat(numRealPt, 4), work(16), coeff(4), xmean(4), xsd(4)
  real*4 residMinus, residPlus, aSlope, bSlope, bintcp, alpha, slope, resid
  real*4 botExtreme, topExtreme, dev
  integer*4 igroup(*), numPntMinus, numRealPt, numPntPlus, numCol
  integer*4 ipt,i
  integer*4 multRegress
  real*4 atand, cosd, sind
  !
  ! first fit a plane to points in the first group
  !
  numCol = 3
  if (numRealPt < 4) numCol = 2
  numPntMinus = 0
  numPntPlus = 0

  ! Load the data matrix, putting the Y data in column 2 and overlaying with group
  ! number if there are only two columns
  do ipt = 1, numRealPt
    if (igroup(ipt) == 1) then
      numPntMinus = numPntMinus + 1
    else
      numPntPlus = numPntPlus + 1
    endif
    xmat(ipt, 1) = xyz(1, ipt)
    xmat(ipt, 2) = xyz(2, ipt)
    xmat(ipt, numCol) = igroup(ipt) - 1
    xmat(ipt, numCol + 1) = xyz(3, ipt)
  enddo

  if (multRegress(xmat, numRealPt, 0, numCol, numRealPt, 1, 0, coeff, 4, bintcpMinus,  &
      xmean, xsd, work) .ne. 0) call exitError( &
      'IN MATRIX INVERSION FOR FITTING PLANES TO 3-D POINTS')

  ! Compute a single slope and angles
  bintcp = bintcpMinus + coeff(numCol)
  coeff(numCol) = 0.
  aSlope = coeff(1)
  bSlope = 0.
  if (numCol == 3) bSlope = coeff(2)
  alpha = atand(bSlope)
  slope = coeff(1) / (cosd(alpha) - bSlope * sind(alpha))

  ! Get averall mean residual and mean for each surface, and the extreme residual above
  ! top and below bottom assuming pitch gets corrected
  resid = 0.
  residPlus = 0.
  residMinus = 0.
  botExtreme = 1.e30
  topExtreme = -1.e30
  do ipt = 1, numRealPt
    if (igroup(ipt) == 1) then
      dev = xyz(3, ipt) - (xyz(1, ipt) * aSlope + xyz(2, ipt) * bSlope + bintcpMinus)
      residMinus = residMinus + abs(dev)
      botExtreme = min(botExtreme, bintcpMinus + dev)
    else
      dev = xyz(3, ipt) - (xyz(1, ipt) * aSlope + xyz(2, ipt) * bSlope + bintcp)
      residPlus = residPlus + abs(dev)
      topExtreme = max(topExtreme, bintcp + dev)
    endif
    resid = resid + abs(dev)
  enddo
  resid = resid / numRealPt
  residPlus = residPlus / numPntPlus
  residMinus = residMinus / numPntMinus

  return
end subroutine twoSurfaceFits

! DNM 5/3/02: give the new maximum tilt angle only for compression
! solutions
!
subroutine calcTiltNew(slopeMinus, tiltMax, iunit2, trueMinus, ifComp, tiltAdd)
  implicit none
  real*4 slopeMinus, tiltMax, trueMinus, tiltAdd
  integer*4 iunit2, ifComp
  real*4 cosNew, globalDelta
  integer*4 iun
  real*4 atand, acosd, cosd, sind

  cosNew = slopeMinus * sind(tiltMax) + cosd(tiltMax)
  globalDelta = atand(-slopeMinus)
  do iun = 6, iunit2
    write(iun, 102) globalDelta, globalDelta + tiltAdd
102 format(' Incremental tilt angle change =',f7.2,/, &
        ' Total tilt angle change =',f13.2)
    if (ifComp .ne. 0) then
      if (abs(cosNew) <= 1.) then
        trueMinus = sign(acosd(cosNew), tiltMax)
        write(iun, 103) tiltMax, trueMinus
103     format('     or, change a fixed maximum tilt angle from' &
            ,f8.2,' to',f8.2,/)
      else
        write(iun, '(a,/)') '      or. . . But cannot derive '// &
            'implied tilt angle: |cos| > 1'
      endif
    else
      write(iun,*)
    endif
  enddo
  return
end subroutine calcTiltNew


subroutine lsfit2Resid(x, y, z, n, a, b, c, alpha, slope, resid, devMin, devMax)
  implicit none
  real*4 x(*), y(*), z(*), a, b, c, alpha, slope, resid, ro, devMin, devMax, dev
  integer*4 n, i
  real*4 sind, cosd, atand

  if (n > 3) then
    call lsfit2(x, y, z, n, a, b, c)
  elseif (n > 2) then
    b = 0.
    call lsfit(x, z, n, a, c, ro)
  else
    a = 0.
    b = 0.
    c = (z(1) + z(n)) / 2.
  endif
  resid = 0.
  devMin = 1.e10
  devMax = -1.e10
  do i = 1, n
    dev = z(i) - (x(i) * a + y(i) * b + c)
    resid = resid + abs(dev)
    devMin = min(devMin, dev)
    devMax = max(devMax, dev)
  enddo
  resid = resid / n
  alpha = atand(b)
  slope = a / (cosd(alpha) - b * sind(alpha))
  return
end subroutine lsfit2Resid
