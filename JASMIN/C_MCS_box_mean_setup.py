from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=cythonize("JASMIN/C_MCS_box_compute_mean_v2.pyx"),
    include_dirs=[numpy.get_include()]
)
