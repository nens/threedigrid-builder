module simulation

    implicit none
    
    contains
    
    subroutine set_state(name, duration, values, interpolation, status)

        use sim_state, only : rain_instance, assign_state_var

        character(LEN=*), intent(in) :: name
        integer, intent(in) :: duration
        double precision, intent(in) :: values(:)
        integer, intent(in), optional :: interpolation
        integer, intent(out) :: status
    
        status = 1
        
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
        status = 0

    end subroutine set_state


    subroutine update_state_duration(name, duration, status)

        use sim_state, only : rain_instance

        character(LEN=*), intent(in) :: name
        integer, intent(in) :: duration
        integer, intent(out) :: status
    
        status = 1
        
        select case(trim(name))
        case('rain')
            call rain_instance%update_duration(duration)
        case default
            write(*,*) '**ERROR Unknown State variable: ', trim(name)
        end select
        status = 0

    end subroutine update_state_duration


    subroutine update(dtc, time)
        
        use timeloop, only : run_sim
        
        implicit none
        
        double precision, intent(in) :: dtc !< Custom timestep size, use -1 to use model default.
        double precision, intent(out) :: time !< Custom timestep size, use -1 to use model default.

        call run_sim(dtc, time)
        
    end subroutine update
    
    subroutine stop_sim(time)
    
        use sim_state, only : do_run, t1
        
        double precision, intent(out) :: time
    
        do_run = 0
        time = t1
    
    end subroutine stop_sim
    
end module simulation
    