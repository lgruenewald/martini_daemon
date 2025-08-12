#!/usr/bin/env python3

from setuptools import Extension, setup, find_packages, Command
from Cython.Build import cythonize
import numpy as np

NAME = "martini_daemon"
INCLUDES = [f"{NAME}/", np.get_include()]
COMPILE_ARGS = ["-O3", "-flto", "-fno-math-errno", "-fno-trapping-math"]
LIBRARIES = ["m"]


pyx_files = ["utils", "detection"]

ext_modules = [
    Extension(
       f"{NAME}.{x}",
       sources=[f"{NAME}/{x}.pyx"],
       include_dirs=INCLUDES,
       extra_compile_args=COMPILE_ARGS,
       libraries=LIBRARIES
    ) for x in pyx_files
]

setup(
    name=NAME,
    version="0.0.1",
    ext_modules=cythonize(ext_modules),
    packages=find_packages(),
    zip_safe=False,
)
