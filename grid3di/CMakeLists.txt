cmake_minimum_required(VERSION 3.5)

file(GLOB grid_SOURCES src/*.f90)

add_library(gridgen3di SHARED ${grid_SOURCES} ) #${shared_SOURCES}
add_dependencies(gridgen3di lgpl)
target_link_libraries(gridgen3di lgpl ${EXTRA_LIBS})

SET_TARGET_PROPERTIES(gridgen3di PROPERTIES SOVERSION ${3di_VERSION_MAJOR})

INSTALL(TARGETS gridgen3di
  LIBRARY DESTINATION lib)

INSTALL (
    DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    FILES_MATCHING PATTERN "*.mod*")