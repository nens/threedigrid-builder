### setup project ###
cmake_minimum_required(VERSION 3.9)

project(fgrid)
enable_language(Fortran)


# Compiler flags
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -ffree-form -fimplicit-none -fbackslash --param max-unroll-times=4 -O3")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -ffree-form -fimplicit-none -fbounds-check -ffpe-trap=zero,overflow,underflow,invalid -fall-intrinsics -fbackslash --param max-unroll-times=4 -g -O0 -fopenmp -Wconversion -Wmaybe-uninitialized")

# Safety net
if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(
    FATAL_ERROR
      "In-source builds not allowed. Please make a new directory (called a build directory) and run CMake from there.\n"
  )
endif()

# Ensure scikit-build modules
if (NOT SKBUILD)
  find_package(PythonInterp 3.8 REQUIRED)
  # Kanged --> https://github.com/Kitware/torch_liberator/blob/master/CMakeLists.txt
  # If skbuild is not the driver; include its utilities in CMAKE_MODULE_PATH
  execute_process(
    COMMAND "${PYTHON_EXECUTABLE}"
    -c "import os, skbuild; print(os.path.dirname(skbuild.__file__))"
    OUTPUT_VARIABLE SKBLD_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  list(APPEND CMAKE_MODULE_PATH "${SKBLD_DIR}/resources/cmake")
  message(STATUS "Looking in ${SKBLD_DIR}/resources/cmake for CMake modules")
endif()

# scikit-build style includes
find_package(PythonExtensions REQUIRED) # for ${PYTHON_EXTENSION_MODULE_SUFFIX}

# Grab the variables from a local Python installation
# NumPy headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}"
  -c "import numpy; print(numpy.get_include())"
  OUTPUT_VARIABLE NumPy_INCLUDE_DIRS
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
# F2PY headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}"
  -c "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Prepping the module
set(f2py_module_name "_fgrid")
set(fortran_src_dir "${CMAKE_SOURCE_DIR}/libthreedigrid")
set(fortran_src_file "${fortran_src_dir}/parameters.f90"
                      "${fortran_src_dir}/array_utils.f90"
                      "${fortran_src_dir}/geo_utils.f90"
                      "${fortran_src_dir}/cells.f90"
                      "${fortran_src_dir}/quadtree.f90"
                      "${fortran_src_dir}/_fgrid.pyf"
                    )
set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})
set(f2py_module_c "${f2py_module_name}module.c")
set(f2py_wrapper "${f2py_module_name}-f2pywrappers2.f90")

set(generated_module_file ${f2py_module_name}${PYTHON_EXTENSION_MODULE_SUFFIX})

add_custom_target(${f2py_module_name} ALL
  DEPENDS "${fortran_src_file}"
  )

add_custom_command(
  OUTPUT "${f2py_module_name}module.c" "${f2py_module_name}-f2pywrappers2.f90"
  COMMAND ${PYTHON_EXECUTABLE} 
      -m "numpy.f2py"
      -m ${f2py_module_name}
      --lower
      ${fortran_src_file}
      #only: make_quadtree set_2d_computational_nodes_lines
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  DEPENDS ${fortran_src_file}
  )

add_library(${generated_module_file} MODULE
            "${f2py_module_name}module.c"
            "${F2PY_INCLUDE_DIR}/fortranobject.c"
            "${f2py_module_name}-f2pywrappers2.f90"
             ${fortran_src_file})

target_include_directories(${generated_module_file} PUBLIC
                            ${NumPy_INCLUDE_DIRS}
                            ${PYTHON_INCLUDE_DIR}                            
                            ${F2PY_INCLUDE_DIR}
                           )
set_target_properties(${generated_module_file} PROPERTIES SUFFIX "")
set_target_properties(${generated_module_file} PROPERTIES PREFIX "")
target_compile_definitions(${generated_module_file} PRIVATE "NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION" 
)


if (UNIX)
  if (APPLE)
    set_target_properties(${generated_module_file} PROPERTIES
    LINK_FLAGS  '-Wl,-rpath,-dylib,-undefined,dynamic_lookup')
  else()
    set_target_properties(${generated_module_file} PROPERTIES
  LINK_FLAGS  '-Wl,--allow-shlib-undefined')
  endif()
endif()

if (SKBUILD)
  install(TARGETS ${generated_module_file} DESTINATION threedigrid_builder/grid/fgrid/)
else()
  install(TARGETS ${generated_module_file} DESTINATION "${CMAKE_SOURCE_DIR}/threedigrid_builder/grid/fgrid/")
endif()