from __future__ import division, print_function, absolute_import

import os
from os.path import join

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('besselaesnew', parent_package, top_path)

    sources = ['src/besselaesnew.f95']

    config.add_extension('besselaesnew', sources=sources)

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
