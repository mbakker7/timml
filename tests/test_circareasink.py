from timml import *
import numpy as np
from numpy.testing import assert_allclose
from numpy.testing import assert_equal

class TestCircAreaSink():
    
    def test_circ_large(self):
        
        ml = ModelMaq(kaq=[1, 2], z=[10, 5, 4, 0], c=20)
        c = CircAreaSink(ml, 0, 0, 100, 0.01)
        ml.solve()

        assert_equal(c.islarge, [False])

        d = 1e-5
        x = 20
        y = 10

        qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
        qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)
        qx, qy = ml.disvec(x, y)
        
        assert_allclose(qx, qxnum)
        assert_allclose(qy, qynum)

        x = 200
        y = -50

        qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
        qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)
        qx, qy = ml.disvec(x, y)
        
        assert_allclose(qx, qxnum)
        assert_allclose(qy, qynum)


    def test_small_resistance(self):

        ml = ModelMaq(kaq=[1, 2], z=[10, 5, 4, 0], c=1e-3)
        c = CircAreaSink(ml, 0, 0, 100, 0.01)
        ml.solve()

        assert_equal(c.islarge, [True])

        d = 1e-5
        x = 99.5 * np.cos(np.pi / 3)
        y = 99.5 * np.sin(np.pi / 3)
        
        qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
        qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)
        qx, qy = ml.disvec(x, y)

        assert_allclose(qx, qxnum)
        assert_allclose(qy, qynum)

        x = 105.5 * np.cos(np.pi / 3)
        y = 105.5 * np.sin(np.pi / 3)
        
        qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
        qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)
        qx, qy = ml.disvec(x, y)
        
        assert_allclose(qx, qxnum)
        assert_allclose(qy, qynum)


    def test_mult_resistance(self):

        ml = ModelMaq(kaq=[1, 2, 4, 2], z=[20, 16, 14, 12, 10, 5, 4, 0], c=[1, 100, 4])
        c = CircAreaSink(ml, 0, 0, 1000, 0.01)
        ml.solve()

        assert_equal(c.islarge, [False, False, True])

        d = 1e-5
        x = 99.5 * np.cos(np.pi / 3)
        y = 99.5 * np.sin(np.pi / 3)
        
        qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
        qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)
        qx, qy = ml.disvec(x, y)

        assert_allclose(qx, qxnum)
        assert_allclose(qy, qynum)

        x = 105.5 * np.cos(np.pi / 3)
        y = 105.5 * np.sin(np.pi / 3)
    
        qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
        qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)
        qx, qy = ml.disvec(x, y)

        assert_allclose(qx, qxnum)
        assert_allclose(qy, qynum)

        x = 995 * np.cos(np.pi / 3)
        y = -995 * np.sin(np.pi / 3)

        qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
        qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)
        qx, qy = ml.disvec(x, y)

        assert_allclose(qx, qxnum)
        assert_allclose(qy, qynum)

        x = 1005 * np.cos(np.pi / 3)
        y = -1005 * np.sin(np.pi / 3)

        qxnum = -ml.aq.T * (ml.head(x + d, y) - ml.head(x - d, y)) / (2 * d)
        qynum = -ml.aq.T * (ml.head(x, y + d) - ml.head(x, y - d)) / (2 * d)
        qx, qy = ml.disvec(x, y)

        assert_allclose(qx, qxnum)
        assert_allclose(qy, qynum)

    def test_storeinput(self):
        
        ml = ModelMaq(kaq=[1, 2, 4, 2], z=[20, 16, 14, 12, 10, 5, 4, 0], c=[1, 100, 4])
        c = CircAreaSink(ml, 0, 0, 1000, 0.01)
        ml.solve()
        
        assert c.inputvalues['xc'] == 0
        assert c.inputvalues['yc'] == 0
        assert c.inputvalues['R'] == 1000
        assert c.inputvalues['N'] == 0.01
        assert c.inputvalues['layer'] == 0
        assert c.inputvalues['name'] == 'CircAreaSink'
        assert c.inputvalues['label'] == None
        assert c.inputvalues['model'].modelname == 'ml'