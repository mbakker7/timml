from __future__ import division, absolute_import, print_function
import sys
import os
import platform

version = {}
with open("timml/version.py") as fp:
    exec(fp.read(), version)

l_d = ""
try:
    import pypandoc

    l_d = pypandoc.convert("README.md", "rst")
except:
    pass

try:
    from numpy.distutils.core import Extension, setup
except ImportError:
    sys.exit("install requires: 'numpy'.")

cputune = ["-march=native"]

if os.name == "nt":
    compile_args = ["-static-libgcc", "-Wall", "-shared"]
else:
    compile_args = ["-static-libgcc", "-Wall", "-lgfortran", "-lquadmath"]
    cputune = []


def setup_package():

    metadata = dict(
        name="timml",
        version=version["__version__"],
        description="Steady multi-layer AEM Model",
        long_description=l_d,
        author="Mark Bakker",
        author_email="markbak@gmail.com",
        url="https://github.com/mbakker7/timml",
        license="MIT",
        packages=["timml"],
        ext_modules=[
            Extension(
                "timml.besselaesnew",
                ["timml/src/besselaesnew.f95"],
                extra_compile_args=compile_args + cputune,
            )
        ],
    )

    setup(**metadata)


if __name__ == "__main__":
    setup_package()
