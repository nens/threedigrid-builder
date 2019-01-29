module simulation

    implicit none
    
    contains

    
    subroutine update(dtc, time)
        
        use timeloop, only : run_sim

        implicit none

        double precision, intent(in) :: dtc !< Custom timestep size, use -1 to use model default.
        double precision, intent(out) :: time 
        integer :: python_iets

        time = run_sim(dtc, func_for)

    end subroutine update

    subroutine func_for(python_iets)
        
        !f2py intent(callback, hide) py_func
        external :: py_func    
        integer, intent(out) :: python_iets

        call py_func(python_iets)

    end subroutine func_for
    

    subroutine stop_sim(time)
    
        use sim_state, only : do_run, t1
        
        double precision, intent(out) :: time
    
        do_run = 0
        time = t1
    
    end subroutine stop_sim

    
    subroutine set_state(name, duration, values, interpolation, status)

        use sim_state, only : assign_state_var, forcing_assigned
        use m_forcing

        character(LEN=*), intent(in) :: name
        integer, intent(in) :: duration
        double precision, intent(in) :: values(:)
        integer, intent(in), optional :: interpolation
        integer, intent(out) :: status
        class(forcing), pointer :: rain_instance
        logical :: assigned

        status = 1
        assigned = forcing_assigned(trim(name))
        if (assigned) then
            write(*,*) '**INFO State already assigned: ', trim(name)
            return
        else
            select case(trim(name))
            case('rain')
                allocate(rain_instance)
                if (present(interpolation)) then
                    call rain_instance%init(name, duration, interpolation)
                else
                    call rain_instance%init(name, duration)
                endif
                call rain_instance%set_data(values)
                call assign_state_var(trim(name), rain_instance)
            case default
                write(*,*) '**ERROR Unknown State variable: ', trim(name)
            end select
        endif
        status = 0

    end subroutine set_state


    subroutine update_state_duration(name, duration, status)

        !use sim_state, only : rain_instance

        character(LEN=*), intent(in) :: name
        integer, intent(in) :: duration
        integer, intent(out) :: status
    
        status = 1
        select case(trim(name))
        case('rain')
            !call rain_instance%update_duration(duration)
        case default
            write(*,*) '**ERROR Unknown State variable: ', trim(name)
        end select
        status = 0

    end subroutine update_state_duration


    subroutine num_active_forcings(n)

        use sim_state, only : forcings

        integer, intent(out) :: n

        if (associated(forcings)) then
            n = forcings%len()
        else
            n = 0
        endif

    end subroutine num_active_forcings

    subroutine num_active_furcings(n)

        integer, intent(out) :: n

        call num_active_forcings(n)

    end subroutine num_active_furcings

    
    subroutine get_active_forcings(size, str_size, active_forcings)

        use sim_state, only : forcings
        use m_forcing

        class(*), pointer :: p
        class(forcing), pointer :: temp
        integer, intent(in) :: str_size
        integer, intent(in) :: size
        character, intent(out) :: active_forcings(size, str_size)
        integer :: n


        if (associated(forcings)) then
            call forcings%iter_restart()
            do n=1,forcings%len()
                p => forcings%iter_next()
                if (associated(p)) then
                    select type (p)
                        type is (forcing)
                        temp => p
                    end select
                endif
                active_forcings(n,:) = string_to_array(trim(temp%type()), str_size)
            enddo
        else
            active_forcings = ''
        endif

    end subroutine get_active_forcings

    pure function string_to_array(string, string_size) result(array)

    character(LEN=*), intent(in) :: string
    integer, intent(in) :: string_size
    character(LEN=1) :: array(string_size)
    integer :: i 

    array = ''
    do i=1, len(trim(string))
        array(i) = string(i:i)
    enddo

end function string_to_array
    

end module simulation
    