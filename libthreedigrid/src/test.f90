module m_test

    implicit none

    type :: ClassTest
        character(LEN=1024) :: description
        double precision, pointer :: vals(:)
    contains
        procedure :: init_test
        procedure :: get_description
        procedure :: get_values
        procedure :: print_values
    end type ClassTest

    contains

    subroutine init_test(self, description, values)

        class(ClassTest) :: self
        character(LEN=1024), intent(in) :: description
        double precision, intent(in) :: values(:)

        self%description = trim(description)
        allocate(self%vals, source=values)

    end subroutine init_test

    function get_description(self) result(descr)

        class(ClassTest) :: self
        character(LEN=1024) :: descr

        descr = trim(self%description)

    end function get_description

    subroutine get_values(self, m0, m1, vals)

        class(ClassTest) :: self
        integer, intent(in) :: m0
        integer, intent(in) :: m1
        double precision, intent(out) :: vals(m1-m0+1)

        write(*,*) 'Digging deeper ', m0, m1
        vals = self%vals(m0:m1)

    end subroutine get_values

    subroutine get_values_by_ref(self, m0, m1, vals_ptr)

        use iso_c_binding

        class(ClassTest), target :: self
        integer, intent(in) :: m0
        integer, intent(in) :: m1
        TYPE(c_ptr), intent(out) :: vals_ptr
        double precision, pointer :: vals(:)

        write(*,*) 'Digging deeper 2', m0, m1
        vals => self%vals(m0:m1)
        write(*,*) 'VALUES: ', vals
        vals_ptr = C_LOC(self%vals(m0:m1))
        write(*,*) 'PTR: ', vals_ptr
        write(*,*) 'LOC DATA: ', C_LOC(self%vals(m0:m1))

    end subroutine get_values_by_ref


    subroutine print_values(self)

        class(ClassTest) :: self

        write(*,*) 'Values are: ', self%vals

    end subroutine print_values

end module m_test