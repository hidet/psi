from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    name = 'exttrigloop',
    ext_modules = cythonize("exttrigloop.pyx"),
    include_dirs = [numpy.get_include()]
)
