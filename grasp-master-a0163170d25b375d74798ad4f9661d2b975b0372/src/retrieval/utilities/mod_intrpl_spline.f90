! **
! **  Copyright 2016 CNRS & Universite Lille 1. All rights reserved.  
! **  Licensed under the GRASP Open Source License V1.0 (see LICENSE file)
! **

        !> @file mod_intrpl_spline.f90
        !> Module contains spline interpolation routines.
        !>
        !>

module mod_intrpl_spline
	 contains
! **********************************************************************
      subroutine intrpl_spline(MM, XX, YY, X, Y, ibeg, icurr, b, c, d)

      integer, intent(in) :: MM
      integer, intent(in) :: ibeg, icurr ! routine entry numbers
      double precision, intent(in) :: XX(MM), YY(MM), X
      double precision, intent(out) :: Y
      double precision, intent(inout) :: b(MM), c(MM), d(MM)
! ibeg, icurr - routine entry numbers
! ibeg = icurr  - compute the coefficients b, c, and d for a cubic interpolating spline
!                 at fist routine call
! icurr > ibeg - use the coefficients computed at first routine call
!
      if(ibeg .gt. icurr) then
        write(*,'(2(a,i0,2x))') 'ibeg = ',ibeg,'icurr = ',icurr
        write(*,'(a)') 'ibeg must be less than icurr'
        stop 'stop in intrpl_spline'
      elseif(ibeg .eq. icurr) then
        call spline(MM, XX, YY, b, c, d)
      endif

      Y = seval(MM, X, XX, YY, b, c, d)

      return
      end subroutine intrpl_spline

! **********************************************************************
! http://www.netlib.org/fmm/
      subroutine spline (n, x, y, b, c, d)
      integer n
      double precision x(n), y(n), b(n), c(n), d(n)
!
!  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  for a cubic interpolating spline
!
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!
!    for  x(i) .le. x .le. x(i+1)
!
!  input..
!
!    n = the number of data points or knots (n.ge.2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!
!  output..
!
!    b, c, d  = arrays of spline coefficients as defined above.
!
!  using  p  to denote differentiation,
!
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!
!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.
!
!
      integer nm1, ib, i
      double precision t
!
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side.
!
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
!
!  forward elimination
!
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
!
!  back substitution
!
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
   40 continue
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return
!
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
      return
      end

! **********************************************************************
! http://www.netlib.org/fmm/
      double precision function seval(n, u, x, y, b, c, d)
      integer n
      double precision  u, x(n), y(n), b(n), c(n), d(n)
!
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
      integer i, j, k
      double precision dx
      data i/1/
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
!
!  binary search
!
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
!
!  evaluate spline
!
   30 dx = u - x(i)
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
      end

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine allocates and initializes arrays for spline interpolation
      subroutine alloc_arrays_spline_intrpl(N, XS, YS, b, c, d)

      implicit none
! ----------------------------------------------------------------------
      integer, intent(in) :: N
      double precision, dimension(:), intent(inout) :: XS, YS, b, c, d
      allocatable :: XS, YS, b, c, d
      integer :: ierr
! ----------------------------------------------------------------------
      allocate(XS(N), stat=ierr)
      if(ierr/=0) stop 'Can not allocate XS array in alloc_arrays_spline_intrpl.'
      XS(:)= 0.0
      allocate(YS(N), stat=ierr)
      if(ierr/=0) stop 'Can not allocate YS array in alloc_arrays_spline_intrpl.'
      YS(:)= 0.0
      allocate(b(N), stat=ierr)
      if(ierr/=0) stop 'Can not allocate b array in alloc_arrays_spline_intrpl.'
      b(:)= 0.0
      allocate(c(N), stat=ierr)
      if(ierr/=0) stop 'Can not allocate c array in alloc_arrays_spline_intrpl.'
      c(:)= 0.0
      allocate(d(N), stat=ierr)
      if(ierr/=0) stop 'Can not allocate d array in alloc_arrays_spline_intrpl.'
      d(:)= 0.0

      end subroutine alloc_arrays_spline_intrpl

! ssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssssss
        !> @brief Routine deallocates arrays for spline interpolation
      subroutine dealloc_arrays_spline_intrpl(XS, YS, b, c, d)

      implicit none
! ----------------------------------------------------------------------
      double precision, dimension(:), intent(inout) :: XS, YS, b, c, d
      allocatable :: XS, YS, b, c, d
      integer :: ierr
! ----------------------------------------------------------------------
      if(allocated(XS)) then
        deallocate(XS, stat=ierr)
        if(ierr/=0) stop 'Can not deallocate XS array in dealloc_arrays_spline_intrpl.'
      endif
      if(allocated(YS)) then
        deallocate(YS, stat=ierr)
        if(ierr/=0) stop 'Can not deallocate YS array in dealloc_arrays_spline_intrpl.'
      endif
      if(allocated(b)) then
        deallocate(b, stat=ierr)
        if(ierr/=0) stop 'Can not dallocate b array in dealloc_arrays_spline_intrpl.'
      endif
      if(allocated(c)) then
        deallocate(c, stat=ierr)
        if(ierr/=0) stop 'Can not deallocate c array in dealloc_arrays_spline_intrpl.'
      endif
      if(allocated(d)) then
        deallocate(d, stat=ierr)
        if(ierr/=0) stop 'Can not deallocate d array in dealloc_arrays_spline_intrpl.'
      endif

      end subroutine dealloc_arrays_spline_intrpl

end module mod_intrpl_spline
