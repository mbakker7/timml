from __future__ import division, absolute_import, print_function
import sys
import os

try:
    from numpy.distutils.core import Extension, setup
except ImportError:
    sys.exit("install requires: 'numpy'.")

f_compile_args = []

ext = Extension(name = 'besselaes',
                sources = ['src/besselaes.f90'],
                extra_f90_compile_args=f_compile_args,
                )

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name = 'timml',
          version="4.0.4",
          packages = ['timml'],
          install_requires = ['numpy>=1.9', 'scipy', 'matplotlib>=1.4'],
          description = "TimML, A Multi-Layer, Analytical Element Model",
          author = 'Mark Bakker',
          author_email = 'mark.bakker@tudelft.nl',
          url = 'https://github.com/mbakker7/timml',
          license = 'MIT',
          ext_modules = [ext]
          )
