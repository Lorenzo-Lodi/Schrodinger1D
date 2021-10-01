program compare
implicit none
integer, parameter :: u1=101, u2=102
integer, parameter :: fp=selected_real_kind(15)
real(kind=fp) :: x1, x2
real(kind=fp) :: f1, f2
integer :: ierr

open(unit=u1, file='qp_0.1.txt', status='old', action='read')
open(unit=u2, file='dp_0.1.txt', status='old', action='read')

do 
 read(u1,*, iostat=ierr) x1, f1
 if(ierr /=0) exit
 read(u2,*, iostat=ierr) x2, f2
 if(ierr /=0) exit
 
! write(*,*) (x1-x2)/x1, (f1-f2)/f1
  if( f1 /=0 ) write(*,*) real(x1), abs(real( (f1-f2)/f1))
enddo

close(u1)
close(u2)
end program compare
