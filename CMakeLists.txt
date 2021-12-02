cmake_minimum_required (VERSION 3.5)

PROJECT (3Di)

###############################################################
##CMAKE OUTPUT PATHS                                         ##
###############################################################
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})


enable_language(Fortran)
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
message (STATUS "COMPILER NAME IS: " ${Fortran_COMPILER_NAME})

###############################################################
##PLATFORM SPECIFICS                                         ##
###############################################################
if (WIN32)
  #get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
  #message (STATUS "COMPILER NAME IS: " ${Fortran_COMPILER_NAME})
  add_definitions(-DWIN32)
elseif (UNIX)
  #get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
  #message (STATUS "COMPILER NAME IS: " ${Fortran_COMPILER_NAME})
  add_definitions(-DUNIX)
elseif (MSVC)
  add_definitions(-DWIN32)
endif()
###############################################################
##PLATFORM SPECIFICS                                         ##
###############################################################
include(FortranCInterface)
FortranCInterface_HEADER(FCMangle.h
                         MACRO_NAMESPACE "FC_"
                         SYMBOL_NAMESPACE "FC_"
                         SYMBOLS mysub mymod:my_sub)

###############################################################
##BUILD TYPE                                                 ##
###############################################################
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)
message (STATUS "BUILD TYPE IS: " ${CMAKE_BUILD_TYPE})

# default installation
get_filename_component (default_prefix "./libthreedigrid" ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory"
      FORCE)
SET(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/include)

###############################################################
##BUILD TYPE                                                 ##
###############################################################



###############################################################
##COMPILER OPTIONS                                           ##
###############################################################
IF (Fortran_COMPILER_NAME MATCHES "gfortran*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -cpp -ffree-line-length-none -ffree-form -fimplicit-none -fall-intrinsics -fbackslash --param max-unroll-times=4 -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -cpp -ffree-line-length-none -ffree-form -fimplicit-none -ffpe-trap=zero,overflow,underflow,invalid -fall-intrinsics -fbackslash --param max-unroll-times=4 -g -Og -Wconversion -Wmaybe-uninitialized")
ELSEIF (Fortran_COMPILER_NAME MATCHES "ifort*")
  # ifort (UNFINISHED)
  # TODO
ELSE (Fortran_COMPILER_NAME MATCHES "gfortran*")
  message (STATUS "Fortran compiler: " ${Fortran_COMPILER_NAME})
  message (STATUS "No optimized Fortran compiler flags are known, we just try something...")
  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -fPIC -ffree-line-length-none -ffree-form -fimplicit-none -fall-intrinsics -fbackslash -funroll-loops --param max-unroll-times=4 -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -fPIC -ffree-line-length-none -ffree-form -fimplicit-none -fall-intrinsics -fbackslash -funroll-loops --param max-unroll-times=4 -g -Og")
ENDIF (Fortran_COMPILER_NAME MATCHES "gfortran*")

IF(CMAKE_BUILD_TYPE MATCHES RELEASE)
  message (STATUS "COMPILER FLAGS ARE: " ${CMAKE_Fortran_FLAGS_RELEASE})
ENDIF()
IF(CMAKE_BUILD_TYPE MATCHES DEBUG)
  message (STATUS "COMPILER FLAGS ARE: " ${CMAKE_Fortran_FLAGS_DEBUG})
ENDIF()
###############################################################
##COMPILER OPTIONS                                           ##
###############################################################


###############################################################
##CMAKE PATHS                                                ##
###############################################################
add_subdirectory(${CMAKE_SOURCE_DIR}/libthreedigrid) # here we use core_makegrid/CMakeLists.txt
