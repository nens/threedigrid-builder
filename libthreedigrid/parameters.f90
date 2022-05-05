module parameters
    
    integer, parameter :: STRINGLEN = 1024
    double precision, parameter   :: NODATA = -9999.0d0
    real, parameter   :: NODATA_R = -9999.0
    
    integer, parameter :: CHUNK_SIZE = 1024
    integer, parameter :: COMPRESSION_LEVEL = 9    
    integer, parameter :: NODE_2D = 1,&
                          NODE_2D_GW = 2,&
                          NODE_1D = 3,&
                          NODE_1D_STOR = 4,&
                          NODE_2D_BOUND = 5,&
                          NODE_2D_GW_BOUND = 6,&
                          NODE_1D_BOUND = 7
    
end module parameters
