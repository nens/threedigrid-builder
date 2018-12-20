module api
    
    use iso_c_binding
    
    
    implicit none
    integer :: MAXDIMS = 6
    character(len=1024) :: config_file

    contains


        function set_paths(c_base_dir, c_result_dir, c_log_dir, c_input_generated_dir, c_preprocessed_dir, c_boundary_dir, &
               c_rain_dir, c_external_input_dir, c_state_dir) result(ierr)
        
            use data_definitions
        
            implicit none
            
            integer :: ierr
            character(LEN=1024), intent(in) :: c_base_dir
            character(LEN=1024), intent(in) :: c_result_dir
            character(LEN=1024), intent(in) :: c_log_dir
            character(LEN=1024), intent(in) :: c_input_generated_dir
            character(LEN=1024), intent(in) :: c_preprocessed_dir
            character(LEN=1024), intent(in) :: c_boundary_dir
            character(LEN=1024), intent(in) :: c_rain_dir
            character(LEN=1024), intent(in) :: c_external_input_dir
            character(LEN=1024), intent(in) :: c_state_dir

            !call mess(LEVEL_DEBUG, 'Now we are going to set paths to model files')
            root = c_base_dir
            !call mess(LEVEL_DEBUG, 'basedir is:', root)
            result_path = c_result_dir
            !call mess(LEVEL_DEBUG, 'resultdir is:', result_path)
            log_path = c_log_dir
            !call mess(LEVEL_DEBUG, 'logdir is:', log_path)
            input_path = c_input_generated_dir
            !call mess(LEVEL_DEBUG, 'inputfilesdir is:', input_path)
            admin_path = c_preprocessed_dir
            !call mess(LEVEL_DEBUG, 'preprocesseddir is:', admin_path)
            boundary_path = c_boundary_dir
            !call mess(LEVEL_DEBUG, 'boundarydir is:', boundary_path)
            rain_path = c_rain_dir
            !call mess(LEVEL_DEBUG, 'raindir is:', rain_path)
            external_input_path = c_external_input_dir
            !call mess(LEVEL_DEBUG, 'external input dir is:', external_input_path)
            state_path = c_state_dir
            !call mess(LEVEL_DEBUG, 'external input dir is:', external_input_path)
            ierr = 0
            
        end function set_paths

        function initialize(c_config_file) result(ierr)
            
            use data_definitions
            use timeloop
            implicit none
            ! Variables
            character(LEN=1024), intent(in) :: c_config_file
            integer :: ierr

            !call mess(LEVEL_INFO, 'initialize')
            
            ierr = 0
            nogrid = 1
            interactive = 1

            config_file = c_config_file
            call read_config(config_file)
            ! NOTE: filepaths() is not called here anymore, paths MUST be set externally via set_paths()
            call load_model(0)
            call init_model()
            call init_sim()

        end function initialize

        function update(dtc) result(time)
        
            use timeloop
            
            implicit none
            
            double precision, intent(in) :: dtc !< Custom timestep size, use -1 to use model default.
            double precision :: time !< Custom timestep size, use -1 to use model default.

            call run_sim(dtc, time)
            
        end function update
        
        function stop_sim() result(time)
        
            use sim_state
            use data_definitions
            
            double precision :: time
        
            do_run = 0
            time = t1
        
        end function stop_sim
        
        
end module api

