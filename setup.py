from setuptools import setup
version = {}
with open("timml/version.py") as fp:
    exec(fp.read(), version)

l_d = ""
try:
    import pypandoc

    l_d = pypandoc.convert("README.md", "rst")
except:
    pass


setup(
    name="timml",
    version=version["__version__"],
    description="Steady multi-layer AEM Model",
    long_description=l_d,
    author="Mark Bakker",
    author_email="markbak@gmail.com",
    url="https://github.com/mbakker7/timml",
    license="MIT",
    packages=["timml", "timml/besselaesnumba"],
    python_requires='>3.5',
    install_requires=["numpy>=1.12", "scipy>=0.19", "numba>=0.40", "matplotlib>=2.0"],
)