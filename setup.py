#!/usr/bin/env python3

from setuptools import Extension, setup, find_packages
from Cython.Build import cythonize

ext_modules = [
    Extension(
        "martini_daemon.utils",
        sources=["martini_daemon/utils.pyx"],
        include_dirs=["."],
        extra_compile_args=["-O3", "-Wall"],
        libraries=["m"]
    )
]

setup(
    name="martini_daemon",
    version="0.0.1",
    ext_modules=cythonize(ext_modules),
    packages=find_packages()
)
