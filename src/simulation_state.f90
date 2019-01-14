module pysimstate

    use sim_state
    !use data_definitions, only : nodall
    
    implicit none
    
    double precision, pointer :: input_rain(:)
    
    contains
    
    subroutine init_state_vars(ierr)
        
        integer, intent(out) :: ierr
    
        ierr = 0
        input_rain => rain_in
        if(.not.associated(input_rain)) then
            ierr = 1
        endif
        write(*,*) input_rain
            
    
    
    end subroutine init_state_vars
    
    subroutine set_state_var_duration(var_name, duration)
    
    character(LEN=*), intent(in) :: var_name
    integer, intent(in) :: duration
    
    select case(trim(var_name))
    case('rain')
        call apply_state_var(trim(var_name), duration)
    end select
    
    end subroutine set_state_var_duration
    
    
end module pysimstate
    