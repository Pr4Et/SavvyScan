! $Id$

! A module to replace common for communicating between func and do3multr, but also to 
! allow solve_wo_outliers to initialize do3multr on first call for reproducible searches
module d3multvars
  implicit none
  real*8 sumSq(4,6), aa(3,3), errMin
  integer*4 nullAxis, magSign, ifTrace, numTrials
  logical firstTime/.true./
  real*8 a11, a12, a13, a21, a22, a23, a31, a32, a33
  equivalence (aa(1, 1), a11), (aa(1, 2), a12), (aa(1, 3), a13)
  equivalence (aa(2, 1), a21), (aa(2, 2), a22), (aa(2, 3), a23)
  equivalence (aa(3, 1), a31), (aa(3, 2), a32), (aa(3, 3), a33)
end module d3multvars
  
!
! SOLVE_WO_OUTLIERS solves for a fit between sets of positions in 3D
! and eliminates outlying position-pairs from the solution
!
! XMAT is the data matrix, MATCOLS is its first dimension and must be at least 17.
! NUMDATA is the full amount of data
! NUMCOL is the number of columns of independent variables
! ICOLFIXED specifies which column has insufficient variance (is fixed)
! and should be negative to indicate inversion in that dimension
! MAXDROP is the maximum number of points to drop as outliers
! CRITPROB and ABSPROBCRIT are the criterion probabilities for considering
! a point an outlier
! If the maximum residual is below ELIMMIN nothing will be eliminated
! NUMDROP is number of points dropped, point numbers returned in IDROP
! AMAT3D is the 3x3 matrix computed
! DELXYZ has the displacements
! CENMEANLOC is the mean location in the dependent variables
! DEVMEAN, DEVSD, DEVMAX are mean, SD, and maximum deviations
! IPNTMAX and DEVXYZMAX given the point number and deviation in X, y, Z
! at which the maximum occurred
!
! The routine copies original data to columns 8 to 14 (or numCol + 5 to
! 2*(numCol + 4)), puts a cross index from the ordered data to the
! original row in column 5 (or numCol + 2), returns a mean residual in
! column 4 (numCol + 1) for the ordered data, and puts residual vector
! components for the ordered data in columns 15 to 17 (2*(numCol + 4) + 1, 2, 3) .
!
subroutine solve_wo_outliers(xMat, matCols, numData, numCol, icolFixed, maxDrop, &
    critProb, absProbCrit, elimMin, idrop, numDrop, aMat3d, delXYZ, cenMeanLoc, devMean, &
    devSD, devMax, ipntMax, devXYZmax)
  use d3multvars
  implicit none
  integer maxErr
  parameter (maxErr = 100000)
  integer*4 idrop(*), matCols
  real*4 xMat(matCols,*)
  real*4 aMat3d(3,*), delXYZ(3), devXYZmax(3), cenMeanLoc(3), cenTmp(3)
  integer*4 index(maxErr)
  integer*4 numData, numCol, maxDrop, numDrop, ipntMax, icolFixed
  real*4 critProb, absProbCrit, elimMin, devMean, devSD, devMax
  integer*4 i, j, lastDrop, itmp, numKeep, jdrop
  real*4 probPerPoint, absPerPoint, sigmaFromMean, sigmaFromSD, sigma, z, prob
  real*4 erfcc, gprob
  gprob(z) = 1. -0.5 * erfcc(z / 1.414214)
  firstTime = .true.
  !
  ! get probability per single point from the overall criterion prob
  !
  probPerPoint = (1. -critProb)**(1. / numData)
  absPerPoint = (1. -absProbCrit)**(1. / numData)
  !
  ! copy the data into columns 8-14 (for ncol = 3)
  !
  do i = 1, numData
    do j = 1, numCol + 4
      xMat(j + numCol + 4, i) = xMat(j, i)
    enddo
  enddo

  call do3multr(xMat, matCols, numData, numCol, numData, icolFixed, aMat3d, delXYZ, &
      cenMeanLoc, devMean, devSD, devMax, ipntMax, devXYZmax)
  numDrop = 0
  !
  ! If returning right away, load cross-indexes into col 5
  if (maxDrop == 0 .or. devMax < elimMin) then
    do i = 1, numData
      xMat(numCol + 2, i) = i
    enddo
    return
  endif
  !
  ! order the residuals
  !
  lastDrop = 0
  if (numData > maxErr) then
    print *, 'CANNOT FIND OUTLIERS: ARRAYS NOT LARGE ENOUGH'
    return
  endif
  !
  ! Sort the residuals and keep index back to initial values
  !
  do i = 1, numData
    index(i) = i
  enddo
  do i = 1, numData - 1
    do j = i + 1, numData
      if (xMat(numCol + 1, index(i)) > xMat(numCol + 1, index(j))) then
        itmp = index(i)
        index(i) = index(j)
        index(j) = itmp
      endif
    enddo
  enddo
  !
  ! load the data in this order
  !
  do i = 1, numData
    do j = 1, numCol + 4
      xMat(j, i) = xMat(j + numCol + 4, index(i))
    enddo
  enddo
  !
  ! Drop successively more points: get mean and S.D. of the remaining
  ! points and check how many of the points pass the criterion
  ! for outliers.
  !
  do jdrop = 1, maxDrop + 1
    call do3multr(xMat, matCols, numData, numCol, numData - jdrop, icolFixed, aMat3d, &
        delXYZ, cenTmp, devMean, devSD, devMax, ipntMax, devXYZmax)
    !
    ! estimate the sigma for the error distribution as the maximum of
    ! the values implied by the mean and the SD of the deviations
    !
    sigmaFromMean = devMean / sqrt(8. / 3.14159)
    sigmaFromSD = devSD / sqrt(3. -8. / 3.14159)
    sigma = max(sigmaFromMean, sigmaFromSD)

    numKeep = 0
    do j = numData - jdrop + 1, numData
      if (sigma < 0.1 * abs(xMat(numCol + 1, j)) .or. sigma < 1.e-5) then
        z = sign(10., xMat(numCol + 1, j))
      else
        z = xMat(numCol + 1, j) / sigma
      endif
      prob = 2 * (gprob(z) - 0.5) - sqrt(2. / 3.14159) * z * exp(-z**2 / 2.)
      if (prob < probPerPoint) numKeep = numKeep + 1
      if (prob >= absPerPoint) numDrop = min(maxDrop, max(numDrop, numData + 1 - j))
    enddo
    !
    ! If all points are outliers, this is a candidate for a set to drop
    ! When only the first point is kept, and all the rest of the points
    ! were outliers on the previous round, then this is a safe place to
    ! draw the line between good data and outliers.  In this case, set
    ! ndrop; and at end take the biggest ndrop that fits these criteria
    !
    if (numKeep == 0) lastDrop = jdrop
    if (numKeep == 1 .and. lastDrop == jdrop - 1 .and. lastDrop > 0) &
        numDrop = lastDrop
    ! print *,'drop', jdrop, ', keep', nkeep, ', lastdrop', lastdrop, &
    ! ',  ndrop =', ndrop
  enddo
  !
  ! when finish loop, need to redo with right amount of data and save
  ! indices in column past the residuals
  !
  do i = 1, numDrop
    idrop(i) = index(numData + i - numDrop)
  enddo
  call do3multr(xMat, matCols, numData, numCol, numData - numDrop, icolFixed, aMat3d, &
      delXYZ, cenTmp, devMean, devSD, devMax, ipntMax, devXYZmax)
  ipntMax = index(ipntMax)
  do i = 1, numData
    xMat(numCol + 2, i) = index(i)
  enddo
  return
end subroutine solve_wo_outliers


! DO3MULTR does the three regressions to determine a transformation
! matrix, or it does a search for a reduced set of parameters if
! ICOLFIXIN is non-zero.
! XMAT is the data matrix, MATCOLS is its first dimension
! NUMDATA is the full amount of data
! NUMCOLIN is the number of columns of independent variables
! NUMFIT is the number of data points to use
! ICOLFIXIN specifies which column has insufficient variance (is fixed)
! and should be negative to indicate inversion in that dimension
! AMAT3D is the 3x3 matrix computed
! DELXYZ has the displacements
! CENMEANLOC is the mean location in the dependent variables
! DEVMEAN, DEVSD, DEVMAX are mean, SD, and maximum deviations
! IPNTMAX and DEVXYZMAX given the point number and deviation in X, y, Z
! at which the maximum occurred
!
subroutine do3multr(xMat, matCols, numData, numColIn, numFit, icolFixIn, aMat3d, delXYZ, &
    cenMeanLoc, devMean, devSD, devMax, ipntMax, devXYZmax)
  use d3multvars
  implicit none
  real*4 xMat(matCols,*), xMeans(matCols), sd(matCols)
  real*4 work(100), b1(matCols), b3(matCols, 3)
  real*4 aMat3d(3,*), delXYZ(3), devXYZ(3), devXYZmax(3), cenMeanLoc(3)
  integer*4 numData, numColIn, numFit, ipntMax, icolFixed, icolFixIn, matCols
  real*4 devMean, devSD, devMax, xMeanSave(6)
  integer*4 ixyz, i, j, ipnt, numColDo, k, km
  real*4 const, devSum, devSq, devPnt, amat(2,2), funcErr
  real*8 sum
  real*4 pp(7,7), yy(7), ptol(6), da(6), ptol1, ftol1, delfac, var(6)
  real*4 ftol2, ptol2
  data da/2., 2., .02, .02, 2., 2./
  integer*4 jmin, iter, icol
  save var

  external func
  !
  ! values for simplex fit
  !
  ptol2 = 5.e-4
  ftol2 = 5.e-4
  ptol1 = 1.e-5
  ftol1 = 1.e-5
  delfac = 2.
  errMin = 1.e30
  numTrials = 0
  ! set to 1 for results of each fit, 2 for trace of new minima, 3 for
  ! full trace of simplex search
  ifTrace = 0

  numColDo = numColIn
  icolFixed = abs(icolFixIn)
  if (icolFixed == 0) then
    !
    ! Simple fit, shove the data down to fill the empty spot
    do i = 1, numFit
      do j = numColDo + 1, numColDo + 3
        xMat(j, i) = xMat(j + 1, i)
      enddo
    enddo
    !
    ! Do the fit and fill the matrix
    call multRegress(xMat, matCols, 1, numColDo, numFit, 3, 0, b3, matCols, delXYZ, &
        xMeans, sd, work)
    do ixyz = 1, 3
      do j = 1, numColDo
        aMat3d(ixyz, j) = b3(j, ixyz)
      enddo
    enddo
    !
    ! restore the data
    do i = 1, numFit
      do j = numColDo + 3, numColDo + 1, -1
        xMat(j + 1, i) = xMat(j, i)
      enddo
    enddo
  else
    nullAxis = icolFixed
    magSign = sign(1, icolFixIn)
    !
    ! if one column is fixed first get sums of squares and cross products
    !
    do icol = 1, 6
      j = icol
      if (icol > 3) j = icol + 1
      sum = 0.
      do i = 1, numFit
        sum = sum + xMat(j, i)
      enddo
      xMeanSave(icol) = sum / numFit
    enddo
    do icol = 1, 6
      j = icol
      if (icol > 3) j = icol + 1
      do ixyz = 1, 4
        k = ixyz
        km = ixyz
        if (ixyz == 4) then
          k = j
          km = icol
        endif
        sum = 0.
        do i = 1, numFit
          sum = sum + (xMat(j, i) - xMeanSave(icol)) * (xMat(k, i) - xMeanSave(km))
        enddo
        sumSq(ixyz, icol) = sum
      enddo
    enddo
    ! print *,'xmsav', (xmsav(i), i=1, 6)
    ! print *,'smsq'
    ! write(*,'(6f12.1)') ((smsq(i, j), j=1, 6), i=1, 4)
    !
    ! now decrement the number of columns to do
    ! subtract that column's independent var from its dependent var (FOR NO CURRENT REASON)
    ! and pack the independent vars into the smaller number of columns
    !
    numColDo = numColIn - 1
    do i = 1, numFit
      xMat(numColIn + 1 + icolFixed, i) = xMat(numColIn + 1 + icolFixed, i) - &
          xMat(icolFixed, i)
      xMat(numColIn + 1, i) = xMat(icolFixed, i)
      do j = icolFixed, numColDo
        xMat(j, i) = xMat(j + 1, i)
      enddo
    enddo
    !
    ! do two multr's, moving the appropriate column of dependent
    ! var data into the one past the packed independent vars
    !
    icol = 1
    do ixyz = 1, 3
      ! print *,'FIT #', ixyz
      if (ixyz .ne. icolFixed) then
        do i = 1, numFit
          xMat(numColDo + 1, i) = xMat(numColIn + 1 + ixyz, i)
          ! if (mod(i, 10) ==1) write(*,'(i4,5f9.2)') i, (xr(j, i), j=1, &
          ! ncoldo+1)
        enddo
        call multRegress(xMat, matCols, 1, numColDo, numFit, 1, 0, b1, matCols, const, &
            xMeans, sd, work)
        do j = 1, numColDo
          aMat3d(icol, j) = b1(j)
        enddo
        ! print *,'solution:', (b1(j), j=1, ncoldo), const
        delXYZ(ixyz) = const
        cenMeanLoc(ixyz) = xMeans(numColDo + 1)
        icol = icol + 1
      endif
    enddo
    !
    ! restore the data
    !
    do i = 1, numFit
      do j = numColDo, icolFixed, -1
        xMat(j + 1, i) = xMat(j, i)
      enddo
      xMat(icolFixed, i) = xMat(numColIn + 1, i)
      xMat(numColIn + 1 + icolFixed, i) = xMat(numColIn + 1 + icolFixed, i) + &
          xMat(icolFixed, i)
    enddo
    !
    ! get starting variables for minimization from matrix the first time,
    ! or just restart from previous run
    !
    amat(1, 1) = aMat3d(1, 1)
    amat(1, 2) = aMat3d(1, 2)
    amat(2, 1) = aMat3d(2, 1)
    amat(2, 2) = aMat3d(2, 2)
    if (icolFixed == 2) then
      amat(1, 2) = -aMat3d(1, 2)
      amat(2, 1) = -aMat3d(2, 1)
    endif

    if (firstTime) then
      call amat_to_rotmag(amat, var(1), var(2), var(3), var(4))
      var(5) = 0.
      var(6) = 0.
    endif
    firstTime = .false.
    call amoebaInit(pp, yy, 7, 6, delfac, ptol2, var, da, func, ptol)
    call amoeba(pp, yy, 7, 6, ftol2, func, iter, ptol, jmin)
    !
    ! per Press et al. recommendation, just restart at current location
    !
    do i = 1, 6
      var(i) = pp(jmin, i)
    enddo
    call amoebaInit(pp, yy, 7, 6, delfac, ptol1, var, da, func, ptol)
    call amoeba(pp, yy, 7, 6, ftol1, func, iter, ptol, jmin)
    !
    ! recover result, get aa matrix from best result, compute dxyz
    !
    do i = 1, 6
      var(i) = pp(jmin, i)
    enddo
    funcErr = sqrt(max(0., yy(jmin) / numFit))
    if (ifTrace .ne. 0) write(*,'(i5,f10.4, 2f9.2,2f9.4,2f9.2)') &
        iter, funcErr, (var(i), i = 1, 6)
    do i = 1, 3
      cenMeanLoc(i) = xMeanSave(i + 3)
      delXYZ(i) = xMeanSave(i + 3)
      do j = 1, 3
        aMat3d(i, j) = aa(i, j)
        delXYZ(i) = delXYZ(i) - aMat3d(i, j) * xMeanSave(j)
      enddo
      if (ifTrace .ne. 0) write(*,102) (aMat3d(i, j), j = 1, 3), delXYZ(i)
102   format(3f10.6,f10.3)
    enddo
  endif
  !
  ! compute deviations for all the data (numData), keep track of max
  ! for the ones in the fit (ndo)
  ! return residual vector components in columns 15 to 17
  !
  devSum = 0.
  devMax = -1.
  devSq = 0.
  do ipnt = 1, numData
    do ixyz = 1, 3
      devXYZ(ixyz) = delXYZ(ixyz) - xMat(numColIn + 1 + ixyz, ipnt)
      do j = 1, numColIn
        devXYZ(ixyz) = devXYZ(ixyz) + aMat3d(ixyz, j) * xMat(j, ipnt)
      enddo
      xMat(2 * (numColIn + 4) + ixyz, ipnt) = devXYZ(ixyz)
    enddo
    devPnt = sqrt(devXYZ(1)**2 + devXYZ(2)**2 + devXYZ(3)**2)
    xMat(numColIn + 1, ipnt) = devPnt
    if (ipnt <= numFit) then
      devSum = devSum + devPnt
      devSq = devSq + devPnt**2
      if (devPnt > devMax) then
        devMax = devPnt
        ipntMax = ipnt
        do i = 1, 3
          devXYZmax(i) = devXYZ(i)
        enddo
      endif
    endif
  enddo
  call sums_to_avgsd(devSum, devSq, numFit, devMean, devSD)
  return
end subroutine do3multr


!
! FUNC is called by the simplex routine to compute the sum squared
! error for the solution vector in solVec
!
subroutine func(solVec, funcErr)
  use d3multvars
  implicit none
  real*4 solVec(*), funcErr
  real*8 b(2,2), c(2,2), d(2,2), smag(3), all(2,2,3)
  equivalence (all, b), (all(1, 1, 2), c), (all(1, 1, 3), d)
  real*8 cosAng1, cosAng2, sinAng1, sinAng2, dtor
  real*4 r4mat(2,2)
  integer*4 ind5, ind6, i
  character*1 starOut
  !
  dtor = 3.1415926536 / 180.
  numTrials = numTrials + 1
  cosAng1 = cos(solVec(5) * dtor)
  sinAng1 = sin(solVec(5) * dtor)
  cosAng2 = cos(solVec(6) * dtor)
  sinAng2 = sin(solVec(6) * dtor)
  !
  ! set up the X, Y, and Z rotation matrices
  !
  if (nullAxis == 1) then
    ind5 = 2
    ind6 = 3
  else if (nullAxis == 2) then
    ind5 = 1
    ind6 = 3
  else
    ind5 = 1
    ind6 = 2
  endif
  !
  call rotmag_to_amat(solVec(1), solVec(2), solVec(3), solVec(4), r4mat)
  all(1, 1, nullAxis) = r4mat(1, 1)
  all(1, 2, nullAxis) = r4mat(1, 2)
  all(2, 1, nullAxis) = r4mat(2, 1)
  all(2, 2, nullAxis) = r4mat(2, 2)
  smag(nullAxis) = magSign * solVec(3)
  smag(ind5) = 1.
  smag(ind6) = 1.
  all(1, 1, ind5) = cosAng1
  all(1, 2, ind5) = -sinAng1
  all(2, 1, ind5) = sinAng1
  all(2, 2, ind5) = cosAng1
  all(1, 1, ind6) = cosAng2
  all(1, 2, ind6) = -sinAng2
  all(2, 1, ind6) = sinAng2
  all(2, 2, ind6) = cosAng2
  c(1, 2) = -c(1, 2)
  c(2, 1) = -c(2, 1)
  !
  ! get products for full a matrix.  If Y is the null axis, put it
  ! first to avoid Z and X rotations playing off against each other
  ! with a 90 degree rotation around Y
  !
  if (nullAxis == 2) then
    !
    ! Product C * B * D:
    !
    a11 = b(2, 1) * c(1, 2) * d(2, 1) + smag(1) * c(1, 1) * d(1, 1)
    a12 = b(2, 1) * c(1, 2) * d(2, 2) + smag(1) * c(1, 1) * d(1, 2)
    a13 = smag(3) * b(2, 2) * c(1, 2)
    a21 = smag(2) * b(1, 1) * d(2, 1)
    a22 = smag(2) * b(1, 1) * d(2, 2)
    a23 = smag(2) * b(1, 2) * smag(3)
    a31 = b(2, 1) * c(2, 2) * d(2, 1) + smag(1) * c(2, 1) * d(1, 1)
    a32 = b(2, 1) * c(2, 2) * d(2, 2) + smag(1) * c(2, 1) * d(1, 2)
    a33 = smag(3) * b(2, 2) * c(2, 2)
  else
    !
    ! Product B * C * D:
    !
    a11 = smag(1) * c(1, 1) * d(1, 1)
    a12 = smag(1) * c(1, 1) * d(1, 2)
    a13 = smag(1) * c(1, 2) * smag(3)
    a21 = b(1, 2) * c(2, 1) * d(1, 1) + smag(2) * b(1, 1) * d(2, 1)
    a22 = b(1, 2) * c(2, 1) * d(1, 2) + smag(2) * b(1, 1) * d(2, 2)
    a23 = smag(3) * b(1, 2) * c(2, 2)
    a31 = b(2, 2) * c(2, 1) * d(1, 1) + smag(2) * b(2, 1) * d(2, 1)
    a32 = b(2, 2) * c(2, 1) * d(1, 2) + smag(2) * b(2, 1) * d(2, 2)
    a33 = smag(3) * b(2, 2) * c(2, 2)
  endif
  !
  ! get error sum
  !
  funcErr = (a11**2 + a21**2 + a31**2) * sumSq(1, 1) + &
      2. * (a11 * a12 + a21 * a22 + a31 * a32) * sumSq(1, 2) + &
      2. * (a11 * a13 + a21 * a23 + a31 * a33) * sumSq(1, 3) + &
      (a12**2 + a22**2 + a32**2) * sumSq(2, 2) + &
      2. * (a12 * a13 + a22 * a23 + a32 * a33) * sumSq(2, 3) + &
      (a13**2 + a23**2 + a33**2) * sumSq(3, 3) - &
      2. * (a11 * sumSq(1, 4) + a21 * sumSq(1, 5) + a31 * sumSq(1, 6)) - &
      2. * (a12 * sumSq(2, 4) + a22 * sumSq(2, 5) + a32 * sumSq(2, 6)) - &
      2. * (a13 * sumSq(3, 4) + a23 * sumSq(3, 5) + a33 * sumSq(3, 6)) + &
      sumSq(4, 4) + sumSq(4, 5) + sumSq(4, 6)

  if (ifTrace > 1) then
    starOut = ' '
    if (funcErr < errMin) then
      starOut = '*'
      errMin = funcErr
    endif
    if (ifTrace > 2 .or. funcErr < errMin) &
        write(*,72) starOut, numTrials, funcErr, (solVec(i), i = 1, 6)
72  format(1x,a1,i4,f15.5, 2f9.3,2f9.5,2f9.3)
  endif

  return
end subroutine func
