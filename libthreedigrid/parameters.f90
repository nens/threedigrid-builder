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

    integer, parameter :: LINE_1D_EMBEDDED = 0,&
                          LINE_1D_ISOLATED = 1,&
                          LINE_1D_CONNECTED = 2,&
                          LINE_1D_LONG_CRESTED = 3,&
                          LINE_1D_SHORT_CRESTED = 4,&
                          LINE_1D_DOUBLE_CONNECTED = 5,&
                          LINE_1D2D_SINGLE_CONNECTED_CLOSED = 51,&
                          LINE_1D2D_SINGLE_CONNECTED_OPEN_WATER = 52,&
                          LINE_1D2D_DOUBLE_CONNECTED_CLOSED = 53,&
                          LINE_1D2D_DOUBLE_CONNECTED_OPEN_WATER = 54,&
                          LINE_1D2D_POSSIBLE_BREACH = 55,&
                          LINE_1D2D_ACTIVE_BREACH = 56,&
                          LINE_1D2D_GROUNDWATER = 57,&
                          LINE_2D = 100,&
                          LINE_2D_OBSTACLE = 101,&
                          LINE_2D_VERTICAL = 150,&
                          LINE_2D_GROUNDWATER = -150,&
                          LINE_2D_BOUNDARY_WEST = 200,&
                          LINE_2D_BOUNDARY_EAST = 300,&
                          LINE_2D_BOUNDARY_SOUTH = 400,&
                          LINE_2D_BOUNDARY_NORTH = 500,&
                          LINE_2D_GROUNDWATER_BOUNDARY_WEST = 600,&
                          LINE_2D_GROUNDWATER_BOUNDARY_EAST = 700,&
                          LINE_2D_GROUNDWATER_BOUNDARY_SOUTH = 800,&
                          LINE_2D_GROUNDWATER_BOUNDARY_NORTH = 900
    
end module parameters
