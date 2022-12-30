program ThreePoint
implicit none
! Results for gfortran 9.3.0 under Linux
!        (decimal_digits) => kind (bytes)
! selected_real_kind( 6)  =>  4   single precision         32-bit
! selected_real_kind(15)  =>  5   double precision         64-bit
! selected_real_kind(18)  => 10   extended precision       80-bit
! selected_real_kind(33)  => 16   double-double precision 128-bit
!
integer, parameter :: fp = selected_real_kind(18)
real(kind=fp) :: h
real(kind=fp) :: E_trial, E_k, E_exact, E_low, E_high
real(kind=fp) :: mass = 2.0_fp
real(kind=fp) :: xmin, xmax, x
integer :: n_of_desired_nodes, i

xmin = -5.5_fp
xmax = -xmin
h=0.00001_fp    ! step; don't go lower than 0.01 for single precision
write(*,'(A, I10)') 'Number of steps = ', int( (xmax-xmin)/h )

! here we set the minimum and maximum
!  find minumum and maximum of the potential function
E_high = -huge(0._fp)
E_low  =  huge(0._fp)
x = xmin
do
  if( x >= xmax) exit
  if( V(x) > E_high) E_high = V(x)
  if( V(x) < E_low) E_low = V(x)
  x = x + h
enddo


! just to be safe, enlarge the range by approx 10%
E_low  = E_low  - 0.05_fp*(E_high-E_low)
E_high = E_high + 0.05_fp*(E_high-E_low)
write(*,*) 'E_low, E_high = ', E_low, E_high

E_trial = (E_high+E_low)/2._fp   ! trial energy (same for all states)

do n_of_desired_nodes = 0, 0
  E_exact = 0.5_fp + n_of_desired_nodes
  E_k = find_eigenvalue(n_of_desired_nodes, E_low, E_high, E_trial)
!  write(*,'(A, ES12.4,I6, 2F22.14, ES14.2)') 'h, n_of_desired_nodes, approx, exact, rel. err. = ', &
!     h, n_of_desired_nodes, E_k, E_exact, (E_k - E_exact)/E_exact
    write(*,'(ES12.4, 2F40.30)') h, E_k
enddo

E_k = V(-1.e10_fp)

contains

!-------------------------------------
real(kind=fp) function find_eigenvalue(n_of_desired_nodes, E_low_input, E_high_input, E_trial)
integer, intent(in) :: n_of_desired_nodes
real(kind=fp), intent(in) :: E_trial, E_low_input, E_high_input
real(kind=fp) :: E, E_high, E_low
integer :: i, imax
! There only a weak (logarithmic) dependence on target_relative_error: halving target_relative_error needs one more iteration,
! ie one more call to is_energy_too_high or count_nodes
! Decreasing target_relative_error by 1000 times increases the iteration by 10
! Decreasing target_relative_error by 1e6 times increases the iteration by 20
real(kind=fp), parameter :: target_relative_error = 2._fp*epsilon(0._fp)

E = E_trial
E_low = E_low_input
E_high = E_high_input

! now we can bisect the energy
imax = 100 ! maximum of 100 bisection, reduces by 2**100 
do i=1, imax
    E = (E_low + E_high)*0.5_fp
    write(*,*) 'i, E_low, E_high, E, count_nodes(xmin, xmax, h, E)', i, E_low, E_high, E, count_nodes(xmin, xmax, h, E)
    if( (E_high - E_low)/abs(E)  < target_relative_error) exit
    if( count_nodes(xmin, xmax, h, E) > n_of_desired_nodes) then
    !if( is_energy_too_high(xmin, xmax, h, E, n_of_desired_nodes) ) then
      E_high = E
    else
      E_low = E
    endif
    
 !  write(*,*) 'Energy has been bracketed by ', E_low, E_high
enddo

find_eigenvalue = E

end function find_eigenvalue
!--------------------------------------------------------
integer function count_nodes(xmin, xmax, h, E)
real(kind=fp), intent(in) :: xmin, xmax, h, E
real(kind=fp) :: f0, f1, f2
real(kind=fp) :: x

    ! first (leftmost) point
    x = xmin
    f0 = 0._fp

    ! second point
    x = x + h
    f1 = 0.0000001_fp ! arbitrary initial value

    count_nodes = 0
    do
      x = x + h
      f2 = 2.0_fp*f1*(1.0_fp-h*h*mass*(E-V(x))) - f0
      if( sign(1.0_fp,f2) /= sign(1.0_fp,f1) )  count_nodes = count_nodes +1
      f0 = f1
      f1 = f2
    !  write(*,*) 'f0, f1, f2', f0, f1, f2
      if( x >= xmax) exit
    enddo

end function count_nodes
!--------------------------------------------------------
logical function is_energy_too_high(xmin, xmax, h, E, n_of_desired_nodes)
real(kind=fp), intent(in) :: xmin, xmax, h, E
integer, intent(in) :: n_of_desired_nodes
integer :: count_nodes
real(kind=fp) :: f0, f1, f2
real(kind=fp) :: x

    ! first (leftmost) point
    x = xmin
    f0 = 0._fp

    ! second point
    x = x + h
    f1 = 0.0000001_fp ! arbitrary initial value

    count_nodes = 0
    do
      x = x + h
      f2 = 2.0_fp*f1*(1.0_fp-h*h*mass*(E-V(x))) - f0
      if( sign(1.0_fp,f2) /= sign(1.0_fp,f1) ) then
         count_nodes = count_nodes +1
         if(count_nodes > n_of_desired_nodes) exit
      endif
      f0 = f1
      f1 = f2
    !  write(*,*) 'f0, f1, f2', f0, f1, f2
      if( x >= xmax) exit
    enddo
    
    is_energy_too_high = count_nodes > n_of_desired_nodes

end function is_energy_too_high

!--------------------------------------------------------
! This is the potential function
function V(x)
real(kind=fp), intent(in) :: x
real(kind=fp) :: V
integer :: nOfCalls =0
save nOfCalls

nOfCalls = nOfCalls + 1

V = x*x

! ugly hack
if(x == -1.e10_fp) write(*,*) 'Total number of calls to V(x) = ', nOfCalls

end function V

end program ThreePoint
