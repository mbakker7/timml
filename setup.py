from __future__ import division, absolute_import, print_function
import sys
import os

try:
    from numpy.distutils.core import Extension, setup
except ImportError:
    sys.exit("install requires: 'numpy'.")

ext1 = Extension(name = 'timml.besselaesnew', 
                 sources = ['timml/src/besselaesnew.f95'])

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
        ext_modules = [ext1]
        )

    setup(**metadata)

if __name__ == '__main__':
    setup_package()
