from timml import *
import numpy as np

class TestUflow():

    def test_storeinput(self):
    
        ml = ModelMaq()
        uf = Uflow(ml, 0.01, 45)
    
        assert uf.angle == 45
        assert uf.slope == 0.01
        assert uf.inputvalues['model'].modelname == 'ml'
        assert uf.inputvalues['angle'] == 45
        assert uf.inputvalues['slope'] == 0.01