cmake_minimum_required(VERSION 3.17.2...3.26)
project(fgrid LANGUAGES C Fortran)

enable_language(Fortran)

if (WIN32)
  add_definitions(-DWIN32)
elseif (UNIX)
  add_definitions(-DUNIX)
elseif (MSVC)
  add_definitions(-DWIN32)
endif()

include(FortranCInterface)
FortranCInterface_HEADER(FCMangle.h
                         MACRO_NAMESPACE "FC_"
                         SYMBOL_NAMESPACE "FC_"
                         SYMBOLS mysub mymod:my_sub)

if (WIN32)
  find_package(Python COMPONENTS Interpreter Development REQUIRED)
  if(Python_FOUND)
      message(STATUS "Python Found: ${Python_EXECUTABLE}")
      message(STATUS "Python Found: ${Python_INCLUDE_DIRS}")
      message(STATUS "Python Found: ${Python_LIBRARIES}")
      message(STATUS "Python Found: ${Python_LIBRARY_DIRS}")
  endif()
endif()

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)

# NumPy headers
execute_process(
  COMMAND "${Python_EXECUTABLE}"
  -c "import numpy; print(numpy.get_include())"
  OUTPUT_VARIABLE NumPy_INCLUDE_DIRS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
# F2PY headers
execute_process(
  COMMAND "${Python_EXECUTABLE}"
  -c "import os, numpy.f2py; print(os.path.join(os.path.dirname(numpy.f2py.__file__), 'src'))"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Make them cmake_paths for windows compatibility.
cmake_path(SET Python_INCLUDE_DIRS ${Python_INCLUDE_DIRS})
cmake_path(SET NumPy_INCLUDE_DIRS ${NumPy_INCLUDE_DIRS})
cmake_path(SET F2PY_INCLUDE_DIR ${F2PY_INCLUDE_DIR})

# Compiler flags
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -cpp -s -ffree-line-length-none -ffree-form -fimplicit-none -fbackslash --param max-unroll-times=4 -O3")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -cpp -s -fPIC -ffree-line-length-none -ffree-form -fimplicit-none -ffpe-trap=zero,overflow,underflow -fall-intrinsics -fbackslash --param max-unroll-times=4 -g -O0 -Wconversion -Wmaybe-uninitialized -fcheck=all -fimplicit-none -Wextra")

# Use a private module name for the extension to avoid collisions and match package import
# This must stay in sync with threedigrid_builder/grid/fgrid/__init__.py
set(f2py_module_name "_fgrid")
set(fortran_src_dir "${CMAKE_SOURCE_DIR}/libthreedigrid")
set(fortran_src_file
  "${fortran_src_dir}/parameters.f90"
  "${fortran_src_dir}/log.f90"
  "${fortran_src_dir}/array_utils.f90"
  "${fortran_src_dir}/geo_utils.f90"
  "${fortran_src_dir}/cells.f90"
  "${fortran_src_dir}/quadtree.f90"
)
set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
set(f2py_module_c "${f2py_module_name}module.c")
set(f2py_wrapper "${f2py_module_name}-f2pywrappers2.f90")
set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

include_directories(${NumPy_INCLUDE_DIRS} ${Python_INCLUDE_DIRS} ${F2PY_INCLUDE_DIR})

# Put generated sources in a dedicated folder to avoid name clashes
set(generated_dir "${CMAKE_CURRENT_BINARY_DIR}/generated")
file(MAKE_DIRECTORY "${generated_dir}")

# Generate f2py sources (C API and Fortran wrappers)
add_custom_command(
  OUTPUT "${generated_dir}/${f2py_module_c}" "${generated_dir}/${f2py_wrapper}"
  COMMAND ${Python_EXECUTABLE} -m numpy.f2py
          -m ${f2py_module_name}
          --lower
          "${fortran_src_dir}/cells.f90"
          "${fortran_src_dir}/quadtree.f90"
  WORKING_DIRECTORY "${generated_dir}"
  DEPENDS ${fortran_src_file}
)

# Mark the generated source
set_source_files_properties("${generated_dir}/${f2py_module_c}" PROPERTIES GENERATED TRUE)
set_source_files_properties("${generated_dir}/${f2py_wrapper}" PROPERTIES GENERATED TRUE)

# Build the Python extension from the generated C and Fortran sources
python_add_library(${f2py_module_name} MODULE
  "${generated_dir}/${f2py_module_c}"
  "${generated_dir}/${f2py_wrapper}"
  "${F2PY_INCLUDE_DIR}/fortranobject.c"
  ${fortran_src_file}
  WITH_SOABI
)

target_include_directories(${f2py_module_name} PUBLIC ${F2PY_INCLUDE_DIR})
target_link_libraries(${f2py_module_name} PUBLIC Python::NumPy)

install(TARGETS ${f2py_module_name} DESTINATION "threedigrid_builder/grid/fgrid/")

# For developer convenience: copy the built extension back into the source tree
option(FGRID_LOCAL_COPY "Copy built extension into source tree for local dev" ON)
if(FGRID_LOCAL_COPY)
  add_custom_command(TARGET ${f2py_module_name} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_SOURCE_DIR}/threedigrid_builder/grid/fgrid"
    COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:${f2py_module_name}>" "${CMAKE_SOURCE_DIR}/threedigrid_builder/grid/fgrid/"
  )
endif()