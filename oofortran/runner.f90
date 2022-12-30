program runner
use UserModule
implicit none
type(User) :: myUser
    
call myUser%constructor('Lorenzo', 'Lodi')
call myUser%printDetails()

end program runner