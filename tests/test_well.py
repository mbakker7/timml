import numpy as np

import timml as tml


def test_wellstring_layers_int():
    # int
    ml = tml.Model3D(z=[0, -1, -2, -3])
    ws = tml.WellStringBase(ml, [(0, 0), (10, 0), (0, 10)], layers=2)
    assert (ws.layers == 2).all() and ws.layers.shape[0] == 3
    assert ws.nlayers == 1


def test_wellstring_layers_list():
    # list
    ml = tml.Model3D(z=[0, -1, -2, -3])
    ws = tml.WellStringBase(ml, [(0, 0), (10, 0), (0, 10)], layers=[1, 2])
    assert (ws.layers == np.array([[1, 2], [1, 2], [1, 2]])).all()
    assert ws.nlayers == 2


def test_wellstring_layers_list_of_tuples():
    # list of differently sized tuples
    ml = tml.Model3D(z=[0, -1, -2, -3])
    ws = tml.WellStringBase(
        ml, [(0, 0), (10, 0), (0, 10)], layers=[(0, 1, 2), (2, 3), (3,)]
    )
    assert ws.layers == [(0, 1, 2), (2, 3), (3,)]
    assert ws.nlayers == 3


def test_wellstring_layers_1d_array():
    # 1d array
    ml = tml.Model3D(z=[0, -1, -2, -3])
    ws = tml.WellStringBase(ml, [(0, 0), (10, 0), (0, 10)], layers=np.arange(1, 4))
    assert (ws.layers == np.array([[1, 2, 3], [1, 2, 3], [1, 2, 3]])).all()
    assert ws.nlayers == 3


def test_wellstring_layers_2d_array():
    # 2d array
    ml = tml.Model3D(z=[0, -1, -2, -3])
    ws = tml.WellStringBase(
        ml, [(0, 0), (10, 0), (0, 10)], layers=np.arange(1, 7).reshape((3, 2))
    )
    assert (ws.layers == np.arange(1, 7).reshape((3, 2))).all()
    assert ws.nlayers == 2
