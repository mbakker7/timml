"""Configure Bessel function backends (numba/fortran).

Switches the internal Bessel implementation used by TimML.
"""
from importlib import import_module
from warnings import warn


def set_bessel_method(method="numba"):
    global bessel
    if method == "fortran":
        try:
            besselaesnew = import_module("timml.src.besselaesnew")
            bessel = besselaesnew.besselaesnew
            bessel.initialize()
        except ImportError:
            warn(
                "Cannot import compiled fortran bessel module! Defaulting to numba!",
                category=ImportWarning,
                stacklevel=1,
            )
            bessel = import_module("timml.besselaesnumba.besselaesnumba")
    elif method == "numba":
        bessel = import_module("timml.besselaesnumba.besselaesnumba")
    else:
        raise ValueError("method must be one of ['fortran', 'numba']")


bessel = None  # is set in timml.__init__ or modified by set_bessel_method()
