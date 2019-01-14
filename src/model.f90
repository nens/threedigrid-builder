module model
    
    use iso_c_binding
    
    implicit none

    contains

        subroutine initialize(config_file, status)

            use m_config, only : read_config

            character(LEN=*), intent(in) :: config_file
            logical :: default_paths
            integer, intent(out) :: status

            !MAYBE DO SOME CHECKS ON CONFIG_FILE. LET'S SEE WHAT IS WISE TO DO.

            status = 1
            default_paths = .TRUE. ! FOR NOW THIS TRUE UNTIL WE CAN SET THIS MORE EASILY FROM THE OUTSIDE.
            call read_config(trim(config_file), default_paths)
            status = 0

        end subroutine initialize

        subroutine set_paths(base_dir, project_dir, result_dir, log_dir, input_generated_dir, preprocessed_dir, boundary_dir, &
               rain_dir, external_input_dir, state_dir, status)
			   
            use m_config, only : set_folder_structure
        
            implicit none
            
            character(LEN=*), intent(in) :: base_dir
            character(LEN=*), intent(in) :: project_dir
            character(LEN=*), intent(in) :: result_dir
            character(LEN=*), intent(in) :: log_dir
            character(LEN=*), intent(in) :: input_generated_dir
            character(LEN=*), intent(in) :: preprocessed_dir
            character(LEN=*), intent(in) :: boundary_dir
            character(LEN=*), intent(in) :: rain_dir
            character(LEN=*), intent(in) :: external_input_dir
            character(LEN=*), intent(in) :: state_dir
            integer, intent(out) :: status
            
            status = 1
            call set_folder_structure(base_dir, project_dir, input_generated_dir, log_dir, result_dir, preprocessed_dir, rain_dir, boundary_dir, external_input_dir, state_dir)
            status = 0
            
        end subroutine set_paths

        subroutine py_load_model(gridadmin_file, griddata_file, status)
            
            use timeloop, only : init_sim

            implicit none

            character(LEN=*), intent(in) :: gridadmin_file
            character(LEN=*), intent(in) :: griddata_file
            integer, intent(out) :: status
            
            status = 1
            call load_model(gridadmin_file, griddata_file)
            call init_sim() !TODO THIS ROUTINE NEEDS TO BE REWRITTEN WITH SIM_STAT STUFF
            status = 0

        end subroutine py_load_model

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
        
        
end module model

