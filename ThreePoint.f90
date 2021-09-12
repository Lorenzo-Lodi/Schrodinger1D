program ThreePoint
implicit none
double precision :: h
double precision :: E
double precision :: mass = 2.0d0
double precision :: xmin, xmax
integer :: n_of_nodes

E = 0.7d0 ! trial energy
xmin = -5.0d0
xmax = -xmin
h=0.0000001d0 ! step

write(*,*) 'n_of_nodes = ', count_nodes(xmin, xmax, h, E)

contains

pure integer function count_nodes(xmin, xmax, h, E)
double precision, intent(in) :: xmin, xmax, h, E
double precision :: f0, f1, f2
double precision :: x

    ! first (leftmost) point
    x = xmin
    f0 = 0.d0

    ! second point
    x = x + h
    f1 = 0.00001d0

    count_nodes = 0
    do
      x = x + h
      f2 = 2.0d0*f1*(1.0d0-h*h*mass*(E-V(x))) - f0
      if( sign(1.0d0,f0) /= sign(1.0d0,f1) )  count_nodes = count_nodes +1
      f0 = f1
      f1 = f2
      if( x >= xmax) exit
    enddo

end function count_nodes

! This is the potential function
pure function V(x)
double precision, intent(in) :: x
double precision :: V

V = x*x

end function V

end program ThreePoint
