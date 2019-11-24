from timml import *

class TestWell():

    def test_storeinput(self):

        ml = ModelMaq()
        well = Well(ml)

        assert well.inputvalues['xw'] == 0
        assert well.inputvalues['yw'] == 0
        assert well.inputvalues['Qw'] == 100.0
        assert well.inputvalues['rw'] == 0.1
        assert well.inputvalues['res'] == 0.0
        assert well.inputvalues['model'].modelname == 'ml'


class TestHeadWell():

    def test_storeinput(self):

        ml = ModelMaq()
        well = HeadWell(ml)

        assert well.inputvalues['xw'] == 0
        assert well.inputvalues['yw'] == 0
        assert well.inputvalues['hw'] == 10
        assert well.inputvalues['res'] == 0
        assert well.inputvalues['model'].modelname == 'ml'
