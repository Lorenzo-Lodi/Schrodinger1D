program convergence
implicit none
! 1 November 2021
! Some tests to produce plots of the convergence of some approximations for the equation
! y''(x) = f(x) y(x)
! give two starting points
integer, parameter :: fp = selected_real_kind( 6) !single precision          32-bit
!integer, parameter :: fp = selected_real_kind(15) !double precision          64-bit
!integer, parameter :: fp = selected_real_kind(18) !extended precision        80-bit
!integer, parameter :: fp = selected_real_kind(33) !double-double precision  128-bit
real(kind=fp), parameter :: Pi = acos(-1._fp)
integer :: i, j, npoints
real(kind=fp) :: h ! step size
real(kind=fp) :: xmin, xmax, x
real(kind=fp) :: y0, y1, y2
real(kind=fp) :: eps

!do j = 3, 40
do j = 7, 7
    npoints = nint( 1.6_fp**real(j, fp) )
    xmin = -0.5_fp
    xmax = 2.0_fp
    h = (xmax-xmin)/real(npoints -1, fp)

    y0 = y_exact(xmin)
    x = xmin+h
    y1 = y_exact(x)

    eps = h**2 / 12._fp
    do i= 2, npoints-1
        ! NB x is pointing to x_n, not x_{n+1}
        !y2 = y1*(2._fp + h**2 * f(x)) - y0 ! 3-point formula
        
        ! Numerov (v1)
!        y2 = (y0*(12._fp - h**2 * f(x-h)) -2._fp*y1*(12._fp + 5._fp * h**2  *f(x) )) / (h**2 * f(x+h) -12._fp)

        ! Numerov (v2) - a rewriting
!        y2 = (2._fp*y1*(1._fp + 5._fp*eps*f(x)) - y0*(1._fp - eps*f(x-h))) / (1._fp - eps*f(x+h))
        
        ! Numerov (v3) - another rewriting
        ! This should be the best for minimizing round-off error
        y2 = ( 2._fp*y1 - y0 + eps*(10._fp*y1*f(x) + y0*f(x-h)) ) / (1._fp - eps*f(x+h))
write(*,*) y2
        x = xmin + real(i, fp)*h
        y0 = y1
        y1 = y2
    enddo

!    write(*,*) npoints, abs(y_exact(x)-y2)
enddo

contains

!------------------------------
pure function f(x)
real(kind=fp), intent(in) :: x
real(kind=fp) :: f 
  f = - pi**2
end function
!------------------------------
pure function y_exact(x)
real(kind=fp), intent(in) :: x
real(kind=fp) :: y_exact 
   y_exact  = cos(x*Pi)
end function
!------------------------------

end program convergence
