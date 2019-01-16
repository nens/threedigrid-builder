module simulation

    implicit none
    
    contains
    
    subroutine set_state(name, duration, values, interpolation, status)

        use sim_state

        character(LEN=*), intent(in) :: name
        integer, intent(in) :: duration
        double precision, intent(in) :: values(:)
        integer, intent(in), optional :: interpolation
        integer, intent(out) :: status
    
        status = 1
        call init_state_var(name, duration, values)
        status = 0

    end subroutine set_state

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
    