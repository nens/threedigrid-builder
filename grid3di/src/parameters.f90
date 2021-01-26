module parameters
    
    integer, parameter :: STRINGLEN = 1024
    integer, parameter :: MANNING = 2,&
                          CHEZY = 1
    double precision, parameter   :: NODATA = -9999.0d0
    real, parameter   :: NODATA_R = -9999.0
    double precision, parameter   :: ddnop = -9999.0d0
    
    integer, parameter :: MAXIMUM = 0,&
                          MINIMUM = 1,&
                          AVERAGE = 2,&
                          CUMULATIVE = 3,&
                          COUNTING = 4
    
    integer, parameter :: CHUNK_SIZE = 1024
    integer, parameter :: COMPRESSION_LEVEL = 9    
    integer, parameter :: NODE_2D = 1,&
                          NODE_2D_GW = 2,&
                          NODE_1D = 3,&
                          NODE_1D_STOR = 4,&
                          NODE_2D_BOUND = 5,&
                          NODE_2D_GW_BOUND = 6,&
                          NODE_1D_BOUND = 7

    double precision, parameter :: G = 9.81d0
    double precision, parameter :: PI = 4.0d0*datan(1.0d0)
    double precision, parameter :: D1BY6 = 1.0d0 / 6.0d0
    double precision, parameter :: RHOREL = 1.293d-3 ! relative density of air versus water
#include <version.inc>

!DEC$ ATTRIBUTES DLLEXPORT, ALIAS: "STRINGLEN"::STRINGLEN
    
end module parameters
