from __future__ import division, absolute_import, print_function
from timml import __version__
from numpy.distutils.core import Extension

ext = Extension(name = 'besselaes',
                sources = ['src/besselaes.f90'])

if __name__ == "__main__":
    from numpy.distutils.core import setup
    setup(name = 'timml',
          version=__version__,
          packages = ['timml'],
          install_requires = ['numpy>=1.9', 'scipy', 'matplotlib>=1.4'],
          description = "TimML, A Multi-Layer, Analytical Element Model",
          author = 'Mark Bakker',
          author_email = 'mark.bakker@tudelft.nl',
          url = 'https://github.com/mbakker7/timml',
          license = 'MIT',
          ext_modules = [ext]
          )
