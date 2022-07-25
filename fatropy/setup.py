from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

fatrop_extension = Extension(
    name="fatropy",
    sources=["fatropy.pyx"],
    libraries=["fatrop"],
    library_dirs=["../build/fatrop"],
    include_dirs=["../fatrop/ocp","../fatrop/aux","../fatrop/solver","../fatrop/blasfeo_wrapper","../fatrop"],
    language="c++",
    define_macros=[("LEVEL1_DCACHE_LINE_SIZE","64")]
)
setup(
    name="fatropy",
    ext_modules=cythonize([fatrop_extension],compiler_directives={'language_level' : "3"})
)