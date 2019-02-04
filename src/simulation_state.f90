module simulation

    implicit none
    
    contains

    
    subroutine update(dtc, has_interaction, time)
        
        use timeloop, only : run_sim

        implicit none

        double precision, intent(in) :: dtc !< Custom timestep size, use -1 to use model default.
        double precision, intent(out) :: time 
        integer :: has_interaction

        if (has_interaction == 1) then
            time = run_sim(dtc, check_interaction)
        else
            time = run_sim(dtc)
        endif

    end subroutine update

    subroutine check_interaction(status)
        
        external :: py_check_interaction
        !f2py intent(callback, hide) py_check_interaction
        integer, intent(out) :: status

        call py_check_interaction(status)

    end subroutine check_interaction
    

    subroutine stop_sim(time)
    
        use sim_state, only : do_run, t1
        
        double precision, intent(out) :: time
    
        do_run = 0
        time = t1
    
    end subroutine stop_sim

    
    subroutine set_state_forcing(name, t_start, t_end, values, interpolation, status)

        use sim_state, only : assign_forcing_var, select_forcing, get_forcings
        use m_lists
        use m_forcing

        character(LEN=*), intent(in) :: name
        double precision, intent(in) :: t_start
        double precision, intent(in) :: t_end
        double precision, intent(in) :: values(:)
        integer, intent(in), optional :: interpolation
        integer, intent(out) :: status
        class(forcing), pointer :: forcing_instance
        class(list), pointer :: forcings
        integer :: list_idx

        status = 1
        forcings => get_forcings()
        if (associated(forcings)) then
            call select_forcing(forcings, trim(name), forcing_instance, list_idx)
            if (associated(forcing_instance)) then
                write(*,*) '**INFO : Forcing already assigned: ', trim(forcing_instance%type())
                return
            endif
        endif
        
        allocate(forcing_instance)
        if (present(interpolation)) then
            call forcing_instance%init(name, t_start, t_end, interpolation)
        else
            call forcing_instance%init(name, t_start, t_end)
        endif
        call forcing_instance%set_data(values)
        call assign_forcing_var(trim(name), forcing_instance)
        status = 0

    end subroutine set_state_forcing


    subroutine remove_forcing(name, status) 

        use sim_state, only : select_forcing, get_forcings
        use m_lists
        use m_forcing

        character(LEN=*), intent(in) :: name
        integer, intent(out) :: status
        class(forcing), pointer :: forcing_instance
        class(list), pointer :: forcings
        integer :: list_idx
        
        status = 1
        forcings => get_forcings()
        if (associated(forcings)) then
            call select_forcing(forcings, trim(name), forcing_instance, list_idx)
            if (associated(forcing_instance)) then
                call forcings%pop(list_idx)
                status = 0
            endif
        endif

        
    end subroutine

    subroutine update_forcing_end(name, t_end, status)

        use sim_state, only : select_forcing, get_forcings
        use m_lists
        use m_forcing

        character(LEN=*), intent(in) :: name
        double precision, intent(in) :: t_end
        integer, intent(out) :: status
        class(forcing), pointer :: forcing_instance
        class(list), pointer :: forcings
        integer :: list_idx
    
        status = 1
        forcings => get_forcings()
        if (associated(forcings)) then
            call select_forcing(forcings, trim(name), forcing_instance, list_idx)
            call forcing_instance%update_end(t_end)
            status = 0
        endif

    end subroutine update_forcing_end


    subroutine num_active_forcings(n)

        use sim_state, only : get_forcings
        use m_lists

        integer, intent(out) :: n
        class(list), pointer :: forcings

        forcings => get_forcings()
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

    
    subroutine get_active_forcings(forcings_size, str_size, active_forcings)

        use sim_state, only : get_forcings
        use m_forcing
        use m_lists

        class(*), pointer :: p
        class(forcing), pointer :: temp
        class(list), pointer :: forcings
        integer, intent(in) :: str_size
        integer, intent(in) :: forcings_size
        character, intent(out) :: active_forcings(forcings_size, str_size)
        integer :: n

        forcings => get_forcings()
        if (associated(forcings)) then
            call forcings%iter_restart()
            do n=1,forcings%len()
                p => forcings%next()
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

    pure function string_to_array(string_in, string_size) result(array)

    character(LEN=*), intent(in) :: string_in
    integer, intent(in) :: string_size
    character(LEN=1) :: array(string_size)
    integer :: i 

    array = ''
    do i=1, len(trim(string_in))
        array(i) = string_in(i:i)
    enddo

end function string_to_array
    

end module simulation
    