cmake_minimum_required(VERSION 3.5)

file(GLOB lgpl_SOURCES src/*.f90)

add_library(lgpl SHARED ${lgpl_SOURCES})

INSTALL(TARGETS lgpl
  LIBRARY DESTINATION lib
)
