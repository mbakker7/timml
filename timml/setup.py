from __future__ import division, print_function, absolute_import

import sys

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('timml', parent_package, top_path)
    config.add_subpackage('besselaesnew')

    return config

if __name__ == '__main__':
    from numpy.disutils.core import setup
    setup(**configuration(top_path='').todict())
