from numpy.distutils.core import setup, Extension
import platform

#import ipdb;ipdb.set_trace()
#comp = distutils.customized_fcompiler()

#ignore_flags = ["/O1", "/assume:underscore", "/names:lowercase"]

#comp.compiler_f90 = [flag for flag in comp.compiler_f90 if flag not in ignore_flags]

if platform.system() == 'Linux':
    f_args = ["-O3", "-ffree-line-length-none", "-ffree-form", "-fimplicit-none"]
    f_macros = [('F2PY_REPORT_ON_ARRAY_COPY', '1'), 
				('NPY_NO_DEPRECATED_API', '1'),
                ('NPY_1_7_API_VERSION', '1')]
    f_include = ['/usr/include', '/usr/include/hdf5/serial', '/opt/threedicore/include']
    f_lib = ['/opt/threedicore/lib']
    f_libs = ['flow']
else:
    f_args = ["/Od", "/debug:full", "/assume:nounderscore"]
    f_macros = [('NO_APPEND_FORTRAN','1'),
                ('F2PY_REPORT_ON_ARRAY_COPY', '1'), 
				('NPY_NO_DEPRECATED_API', '1')]
    f_include = ['../threedi-calculationcore/core_flow1d2d/x64/Debug', '../threedi-calculationcore/lgpl/x64/Debug']
    f_lib = ['../threedi-calculationcore/core_flow1d2d/x64/Debug', '../threedi-calculationcore/lgpl/x64/Debug']
    f_libs = ['libflow']

api_module = Extension(
    name= 'py3di',
    sources = ["src/model.f90", 
               "src/simulation_state.f90"],
    define_macros = f_macros,
    include_dirs=f_include,
    library_dirs=f_lib,
    libraries=f_libs,
    extra_f90_compile_args = f_args,
    )

setup(
    name = 'py3di',
    version = '0.1',
    description = 'Python package with fortran API code',
    ext_modules = [api_module],
)
