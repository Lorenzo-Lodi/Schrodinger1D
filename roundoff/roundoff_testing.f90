program ThreePoint
implicit none
! Results for gfortran 9.3.0 under Linux
!        (decimal_digits) => kind (bytes)
! selected_real_kind( 6)  =>  4   single precision         32-bit
! selected_real_kind(15)  =>  5   double precision         64-bit
! selected_real_kind(18)  => 10   extended precision       80-bit
! selected_real_kind(33)  => 16   double-double precision 128-bit
!
integer, parameter :: sp = selected_real_kind(6)
integer, parameter :: dp = selected_real_kind(15)
integer, parameter :: xp = selected_real_kind(18)
integer, parameter :: qp = selected_real_kind(33)

integer, parameter :: fp = dp
real(kind=fp) :: h
real(kind=fp) :: E
real(kind=fp) :: mass = 2.0_fp
real(kind=fp) :: xmin, xmax, x
! integer :: i
real(kind=fp) :: f0, f1, f2
character(len=100) :: fmt

xmin = -5.5_fp
xmax = -xmin
h=0.0001_fp    ! step; don't go lower than 0.01 for single precision
!write(*,'(A, I10)') 'Number of steps = ', int( (xmax-xmin)/h )

E = 0.854_fp

fmt = '(F40.30, ES40.30)'

! first (leftmost) point
x = xmin
f0 = 0._fp
write(*, fmt) x, f0

! second point
x = x + h
f1 = 0.0000001_fp ! arbitrary initial value
write(*, fmt) x, f1

do
   x = x + h
   f2 = 2.0_fp*f1*(1.0_fp-h*h*mass*(E-V(x))) - f0
   write(*, fmt) x, f2
   f0 = f1
   f1 = f2
   if( x >= xmax) exit
enddo


contains

!--------------------------------------------------------
! This is the potential function
pure function V(x)
real(kind=fp), intent(in) :: x
real(kind=fp) :: V

V = x*x

end function V

end program ThreePoint
