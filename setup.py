from os import path

from setuptools import setup

version = {}
with open("timml/version.py") as fp:
    exec(fp.read(), version)

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    l_d = f.read()

setup(
    name="timml",
    version=version["__version__"],
    description="Steady multi-layer AEM Model",
    long_description=l_d,
    long_description_content_type='text/markdown',
    author="Mark Bakker",
    author_email="markbak@gmail.com",
    url="https://github.com/mbakker7/timml",
    license="MIT",
    packages=["timml", "timml/besselaesnumba"],
    python_requires='>3.5',
    install_requires=["numpy>=1.12", "scipy>=0.19",
                      "numba>=0.39", "matplotlib>=3.0"],
    classifiers=['Topic :: Scientific/Engineering :: Hydrology'],
)
