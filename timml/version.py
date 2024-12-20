from importlib import metadata
from platform import python_version

__version__ = "6.6.0"


def show_versions() -> None:
    """Print the version of dependencies."""
    msg = (
        f"TimML version      : {__version__}\n\n"
        f"Python version     : {python_version()}\n"
        f"Numpy version      : {metadata.version('numpy')}\n"
        f"Numba version      : {metadata.version('numba')}\n"
        f"Scipy version      : {metadata.version('scipy')}\n"
        f"Pandas version     : {metadata.version('pandas')}\n"
        f"Matplotlib version : {metadata.version('matplotlib')}"
    )
    print(msg)
