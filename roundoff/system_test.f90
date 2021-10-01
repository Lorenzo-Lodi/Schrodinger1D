program test
implicit none
integer :: i
integer, parameter :: f32  = selected_real_kind( 6) ! single precision    32-bit
integer, parameter :: f64  = selected_real_kind(15) ! double precision    64-bit
integer, parameter :: f80  = selected_real_kind(18) ! extended precision  80-bit
integer, parameter :: f128 = selected_real_kind(33) ! quad precision     128-bit

write(*,*) 'single   precision eps = ', epsilon(0._f32)
write(*,*) 'double   precision eps = ', epsilon(0._f64)
write(*,*) 'extended precision eps = ', epsilon(0._f80)
write(*,*) 'quad     precision eps = ', epsilon(0._f128)

end program test
