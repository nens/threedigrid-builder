cmake_minimum_required (VERSION 3.5)

PROJECT (3Di)

###############################################################
##CMAKE OUTPUT PATHS                                         ##
###############################################################
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})
SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
configure_file (
  "${CMAKE_SOURCE_DIR}/version.inc.in"
  "${CMAKE_Fortran_MODULE_DIRECTORY}/version.inc"
  )


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
##VERSION INFO                                               ##
###############################################################
file (STRINGS "VERSION.rst" VERSION_NUMBER)
message (STATUS 3Di version: ${VERSION_NUMBER}) 
string(REPLACE "." ";" VERSION_LIST ${VERSION_NUMBER})
message (STATUS "list is " ${VERSION_LIST}) 
list(LENGTH VERSION_LIST LEN_VERSION)
list(GET VERSION_LIST 0 3di_VERSION_MAJOR)
list(GET VERSION_LIST 1 3di_VERSION_MINOR)
list(GET VERSION_LIST 2 3di_VERSION_PATCH)
message (STATUS "Major " ${3di_VERSION_MAJOR}) 
message (STATUS "Minor " ${3di_VERSION_MINOR}) 
message (STATUS "Patch " ${3di_VERSION_PATCH}) 
if(LEN_VERSION EQUAL 4)
  list(GET VERSION_LIST 3 3di_VERSION_TWEAK)
  message (STATUS "Tweak " ${3di_VERSION_TWEAK}) 
  set (3di_VERSION ${3di_VERSION_MAJOR}.${3di_VERSION_MINOR}.${3di_VERSION_PATCH}.${3di_VERSION_TWEAK})
else(LEN_VERSION EQUAL 4)
  set (3di_VERSION ${3di_VERSION_MAJOR}.${3di_VERSION_MINOR}.${3di_VERSION_PATCH})
endif()

#if (3di_VERSION_TWEAK 
set (VERSION ${3di_VERSION})
message (STATUS "3Di version: " ${3di_VERSION}) 
###############################################################
##VERSION INFO                                               ##
###############################################################





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
get_filename_component (default_prefix "/usr/local/" ABSOLUTE)
set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
      "Choose the installation directory"
      FORCE)

###############################################################
##BUILD TYPE                                                 ##
###############################################################



###############################################################
##COMPILER OPTIONS                                           ##
###############################################################
IF (Fortran_COMPILER_NAME MATCHES "gfortran*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -fPIC -cpp -ffree-line-length-none -ffree-form -fimplicit-none -fall-intrinsics -fbackslash --param max-unroll-times=4 -O3 -ffast-math -fopenmp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -fPIC -cpp -ffree-line-length-none -ffree-form -fimplicit-none -ffpe-trap=zero,overflow,underflow,invalid -fall-intrinsics -fbackslash --param max-unroll-times=4 -g -Og -ffast-math -fopenmp -Wconversion -Wmaybe-uninitialized")
ELSEIF (Fortran_COMPILER_NAME MATCHES "ifort*")
  # ifort (UNFINISHED)
  # TODO
ELSE (Fortran_COMPILER_NAME MATCHES "gfortran*")
  message (STATUS "Fortran compiler: " ${Fortran_COMPILER_NAME})
  message (STATUS "No optimized Fortran compiler flags are known, we just try something...")
  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -fPIC -ffree-line-length-none -ffree-form -fimplicit-none -fall-intrinsics -fbackslash -funroll-loops --param max-unroll-times=4 -O3 -ffast-math -fopenmp")
  set (CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -fPIC -ffree-line-length-none -ffree-form -fimplicit-none -fall-intrinsics -fbackslash -funroll-loops --param max-unroll-times=4 -g -Og -ffast-math")
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

#add_subdirectory(${CMAKE_SOURCE_DIR}/lgpl) # here we use lgpl/CMakeLists.txt

add_subdirectory(${CMAKE_SOURCE_DIR}/lgpl) # here we use lgpl/CMakeLists.txt  ##FIX THIS TO BE SEPARATE LIBRARY FOR LINKING

add_subdirectory(${CMAKE_SOURCE_DIR}/libgrid) # here we use core_makegrid/CMakeLists.txt


###############################################################
##PACKAGING                                                  ##
###############################################################
if (WIN32)
    if (USE_WIX_TOOLSET)
        set(CPACK_GENERATOR "WIX") # this need WiX Tooset installed and a path to candle.exe
    else ()
        set(CPACK_GENERATOR "NSIS") # this needs NSIS installed, and available
    endif ()
    set(CPACK_SOURCE_GENERATOR "ZIP")
elseif(UNIX)
    set(CPACK_GENERATOR "DEB")
    set(CPACK_SOURCE_GENERATOR "TGZ")
else(WIN32)
	MESSAGE(SEND_ERROR "No Platform found for package generation")
endif(WIN32)

set(CPACK_PACKAGE_NAME "3digrid")
set(CPACK_SYSTEM_NAME "${CMAKE_SYSTEM_PROCESSOR}")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Martijn Siemerink") #required
set(CPACK_PACKAGING_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")
set(CPACK_PACKAGE_VERSION "${3di_VERSION}")
set(CPACK_PACKAGE_VERSION_MAJOR "${3di_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${3di_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${3di_VERSION_PATCH}")

set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_SOURCE_DIR}/CHANGES.rst")
set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.rst")
set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_SOURCE_DIR}/build/;${CMAKE_SOURCE_DIR}/.git/")

include(CPack)