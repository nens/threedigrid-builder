cmake_minimum_required(VERSION 3.5)

file(GLOB grid_SOURCES src/*.f90)

add_library(threedigrid SHARED ${grid_SOURCES} ) #${shared_SOURCES}
add_dependencies(threedigrid lgpl)
target_link_libraries(threedigrid lgpl ${EXTRA_LIBS})

SET_TARGET_PROPERTIES(threedigrid PROPERTIES SOVERSION ${3di_VERSION_MAJOR})

INSTALL(TARGETS threedigrid
  LIBRARY DESTINATION lib)

INSTALL (
    DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    FILES_MATCHING PATTERN "*.mod*")