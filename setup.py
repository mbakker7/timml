from __future__ import division, absolute_import, print_function
import sys
import os
import platform

try:
    from numpy.distutils.core import Extension, setup
except ImportError:
    sys.exit("install requires: 'numpy'.")

cputune = ['-march=native',]

if os.name == 'nt':
    compile_args = ['-static']
else: 
    compile_args = ['-static-libgcc', '-Wl', '-lgfortran', '-lquadmath']
    cputune = []

def setup_package():
    
    metadata = dict(
        name = 'timml',
        version="5.0.0",
        install_requires = ['numpy>=1.12', 'scipy>=0.19', 'matplotlib>=2.0'],
        description = "TimML, A Multi-Layer, Analytical Element Model",
        author = 'Mark Bakker',
        author_email = 'mark.bakker@tudelft.nl',
        url = 'https://github.com/mbakker7/timml',
        license = 'MIT',
        packages = ['timml'],
        ext_modules = [
          Extension(
              'timml.besselaesnew', ['timml/src/besselaesnew.f95'],
              extra_compile_args=compile_args + cputune
              )
          ]
        )

    setup(**metadata)

if __name__ == '__main__':
    setup_package()
