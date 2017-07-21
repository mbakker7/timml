import sys
import os

try:
    from numpy.distutils.core import Extension, setup
except ImportError:
    sys.exit("install requires: 'numpy'.")

def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration('timml', parent_package, top_path)
    config.add_extension('besselaesnew',
                         sources = ['timml/src/besselaesnew.f95'])

    config.get_version(os.path.join('timml', 'version.py'))

    config.make_config_py()
    return config

def setup_package():

    metadata = dict(
        name = 'timml',
        version="5.0.0",
        install_requires = ['numpy>=1.9', 'scipy', 'matplotlib>=1.4'],
        description = "TimML, A Multi-Layer, Analytical Element Model",
        author = 'Mark Bakker',
        author_email = 'mark.bakker@tudelft.nl',
        url = 'https://github.com/mbakker7/timml',
        license = 'MIT',
        packages = ['timml'],
        )

    metadata['configuration'] = configuration

    setup(**metadata)

if __name__ == '__main__':
    setup_package()
