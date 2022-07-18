from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

fatrop_extension = Extension(
    name="fatropy",
    sources=["fatropy.pyx"],
    libraries=["fatrop"],
    library_dirs=["../build/fatrop"],
    include_dirs=["../fatrop/ocp","../fatrop"],
    language="c++",
    define_macros=[("LEVEL1_DCACHE_LINE_SIZE","8")]
)
setup(
    name="fatropy",
    ext_modules=cythonize([fatrop_extension])
)