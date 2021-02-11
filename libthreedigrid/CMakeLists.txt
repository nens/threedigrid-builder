cmake_minimum_required(VERSION 3.5)

file(GLOB grid_SOURCES src/*.f90)

add_library(libthreedigrid SHARED ${grid_SOURCES} ) #${shared_SOURCES}
add_dependencies(libthreedigrid lgpl)
target_link_libraries(libthreedigrid lgpl ${EXTRA_LIBS})

SET_TARGET_PROPERTIES(libthreedigrid PROPERTIES SOVERSION ${3di_VERSION_MAJOR})

INSTALL(TARGETS libthreedigrid
  LIBRARY DESTINATION lib)

INSTALL (
    DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}
    DESTINATION ${CMAKE_INSTALL_PREFIX}
    FILES_MATCHING PATTERN "*.mod*")