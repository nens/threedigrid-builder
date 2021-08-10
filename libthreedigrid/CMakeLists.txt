cmake_minimum_required(VERSION 3.5)

file(GLOB grid_SOURCES src/*.f90)

add_library(threedigrid SHARED ${grid_SOURCES} ) #${shared_SOURCES}
target_link_libraries(threedigrid)

SET_TARGET_PROPERTIES(threedigrid PROPERTIES SOVERSION ${3di_VERSION_MAJOR})

INSTALL(TARGETS threedigrid
  PUBLIC_HEADER ARCHIVE LIBRARY DESTINATION lib)
