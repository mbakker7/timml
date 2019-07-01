import numpy as np
import pytest
from numpy.testing import assert_allclose
try:
    from timml.besselaesnew import besselaesnew
    besselaesnew.initialize()
except ModuleNotFoundError:
    pass

"""
x=2.0
y=1.0
z1=-3.0 - 1.0j
z2=2.0 + 2.0j
lambda=[0.0, 2.0, 11.0]
order = (0,1)
ilap=1
naq=3

lambda is a keyword in Python, so use order of function arguments

"""


@pytest.mark.skip(reason="no fortran extension by default")
def potbesldho(x):
    pot = besselaesnew.potbesldho(2.0, 1.0, np.complex(-3.0, -1.0), np.complex(2.0, 2.0),
                                  [0.0, 2.0, 11.0], x, 1, 3)
    return pot


@pytest.mark.skip(reason="no fortran extension by default")
def test_potbesldho():
    assert_allclose(potbesldho(0), np.array(
        [-0.31055947, -0.23498503, -0.30327438]))
    assert_allclose(potbesldho(1), np.array(
        [-0.17694283, -0.15257055, -0.17583515]))


@pytest.mark.skip(reason="no fortran extension by default")
def test_potbesldv():
    potv = besselaesnew.potbesldv(2.0, 1.0, np.complex(-3.0, -1.0), np.complex(2.0, 2.0),
                                  [0.0, 2.0, 11.0], 1, 1, 3)
    assert_allclose(potv[0], np.array([-0.31055947, -0.23498503, -0.30327438]))
    assert_allclose(potv[1], np.array([-0.17694283, -0.15257055, -0.17583515]))


@pytest.mark.skip(reason="no fortran extension by default")
def test_disbesldho():
    qxqy_zero = besselaesnew.disbesldho(2.0, 1.0, np.complex(-3.0, -1.0), np.complex(2.0, 2.0),
                                        [0.0, 2.0, 11.0], 0, 1, 3)
    assert_allclose(qxqy_zero[0], np.array(
        [-0.170131146, -0.18423853, -0.173157849]))
    assert_allclose(qxqy_zero[1], np.array(
        [0.0274405074, 0.0888068675, 0.0342656083]))

    qxqy_one = besselaesnew.disbesldho(2.0, 1.0, np.complex(-3.0, -1.0), np.complex(2.0, 2.0),
                                       [0.0, 2.0, 11.0], 1, 1, 3)
    assert_allclose(qxqy_one[0], np.array(
        [-0.10412493, -0.1084466406, -0.104477618]))
    assert_allclose(qxqy_one[1], np.array(
        [0.106176131, 0.1162738781, 0.1067421121]))


@pytest.mark.skip(reason="no fortran extension by default")
def test_disbesldv():
    qxqyv = besselaesnew.disbesldv(2.0, 1.0, np.complex(-3.0, -1.0), np.complex(2.0, 2.0),
                                   [0.0, 2.0, 11.0], 1, 1, 3)
    assert_allclose(qxqyv[0], np.array(
        [-0.17013114606375021, -0.18423853257632447, -0.17315784943727297]))
    assert_allclose(qxqyv[2], np.array(
        [2.7440507429637e-002, 8.880686745447e-002, 3.426560831291e-002]))
    assert_allclose(qxqyv[1], np.array(
        [-0.10412493484448178, -0.10844664064434061, -0.10447761803194042]))
    assert_allclose(qxqyv[3], np.array(
        [0.10617613097471285, 0.11627387807684744, 0.10674211206906066]))
