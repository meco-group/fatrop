#
# Fatrop - A fast trajectory optimization solver
# Copyright (C) 2022, 2023 Lander Vanroye, KU Leuven. All rights reserved.
#
# This file is part of Fatrop.
#
# Fatrop is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Fatrop is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Fatrop.  If not, see <http://www.gnu.org/licenses/>.#
import setuptools
from Cython.Build import cythonize

with open("README.md", 'r') as f:
    long_description = f.read()

fatrop_extension = setuptools.Extension(
    name="fatrop.fatropy",
    sources=["src/fatrop/fatropy/fatropy.pyx"],
    libraries=["fatrop"],
    library_dirs=["../build/fatrop"],
    # runtime_library_dirs=["INSTALLATION FOLDER"],
    include_dirs=["../fatrop/ocp","../fatrop/auxiliary","../fatrop/solver","../fatrop/blasfeo_wrapper","../fatrop/templates","../fatrop", "/opt/blasfeo/include", "../external/blasfeo/include","src/fatrop/fatropy"],
    language="c++",
    define_macros=[("LEVEL1_DCACHE_LINE_SIZE","64")]
)
setuptools.setup(
    package_dir={"": "src"},
    long_description=long_description,
    packages=setuptools.find_namespace_packages(where='src'),
    ext_modules=cythonize([fatrop_extension],compiler_directives={'language_level' : "3"})
)
