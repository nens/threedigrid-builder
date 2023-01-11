module m_array_utils

    implicit none 


    interface check_bounds
        module procedure :: check_bounds_d_1d
        module procedure :: check_bounds_d_2d
        module procedure :: check_bounds_i_1d
        module procedure :: check_bounds_i_2d
    end interface
    
    contains

    function check_bounds_d_1d(array, i0, i1) result(inbound)

        double precision, intent(in) :: array(:)
        integer, intent(in) :: i0
        integer, intent(in) :: i1
        logical :: inbound

        inbound = .FALSE.
        if (i0 >= lbound(array, 1) .and. i1 <= ubound(array, 1)) then
            if(i0 <= i1) then
                inbound = .TRUE.
            else
                inbound = .FALSE.
            endif
        else
            write(*,*) '** WARNING: Array out of bounds'
        endif

    end function check_bounds_d_1d

    function check_bounds_d_2d(array, i0, i1, axis_in) result(inbound)

        double precision, intent(in) :: array(:,:)
        integer, intent(in) :: i0
        integer, intent(in) :: i1
        integer, intent(in), optional :: axis_in
        integer :: axis
        logical :: inbound

        if (present(axis_in)) then
            axis = axis_in
        else
            axis = 1
        endif

        inbound = .FALSE.
        if (i0 >= lbound(array, axis) .and. i1 <= ubound(array, axis)) then
            if(i0 <= i1) then
                inbound = .TRUE.
            else
                write(*,*) '** WARNING: Indices in incorrect order.'
                inbound = .FALSE.
            endif
        else
            write(*,*) '** WARNING: Array out of bounds'
        endif

    end function check_bounds_d_2d

    function check_bounds_i_1d(array, i0, i1) result(inbound)

        integer, intent(in) :: array(:)
        integer, intent(in) :: i0
        integer, intent(in) :: i1
        logical :: inbound

        inbound = .FALSE.
        if (i0>=lbound(array, 1) .and. i1<=ubound(array, 1)) then
            if(i0 <= i1) then
                inbound = .TRUE.
            else
                inbound = .FALSE.
            endif
        else
            write(*,*) '** WARNING: Array out of bounds'
        endif

    end function check_bounds_i_1d

    function check_bounds_i_2d(array, i0, i1, axis_in) result(inbound)

        integer, intent(in) :: array(:,:)
        integer, intent(in) :: i0
        integer, intent(in) :: i1
        integer, intent(in), optional :: axis_in
        integer :: axis
        logical :: inbound

        if (present(axis_in)) then
            axis = axis_in
        else
            axis = 1
        endif

        inbound = .FALSE.
        if (i0>=lbound(array, axis) .and. i1<=ubound(array, axis)) then
            if(i0 <= i1) then
                inbound = .TRUE.
            else
                inbound = .FALSE.
            endif
        else
            write(*,*) '** WARNING: Array out of bounds'
        endif

    end function check_bounds_i_2d



end module m_array_utils