from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import numpy

extensions = [
    Extension("_simulateQDSM_core",
              ["_simulateQDSM_core.py",
               "_simulateQDSM_core.pxd"],
        include_dirs = [numpy.get_include()],
        libraries = [],
        library_dirs = []),
#    Extension("*", ["*.pyx"],
#        include_dirs = [numpy.get_include()],
#        libraries = [],
#        library_dirs = []),
]
setup(
    name = "My hello app",
    ext_modules = extensions,
    #ext_modules = cythonize(extensions),
    cmdclass = {'build_ext':build_ext}
)
