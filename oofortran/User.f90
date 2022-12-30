module UserModule
   private

   type, public :: User
      character(len=20), private :: name
      character(len=20), private :: surname
   contains
      procedure :: constructor
      procedure :: printDetails
   end type

contains

   subroutine constructor(this, name, surname)
    class(User), intent(out) :: this
    character(len=20), intent(in) :: name
    character(len=20), intent(in) :: surname
    this%name = name
    this%surname = surname
   end subroutine

   subroutine printDetails(this)
      class(User), intent(in) :: this
      write (*, *) 'Name is '//trim(this%name)//' and surname is '//trim(this%surname)
   end subroutine

end module
