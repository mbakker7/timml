from timml import *

class TestConstant():

    def test_storeinput(self):

        ml = ModelMaq()
        c = Constant(ml)

        assert c.inputvalues['model'].modelname == 'ml'


class TestConstantStar():

    def test_storeinput(self):

        ml = ModelMaq()
        c = ConstantStar(ml)

        assert c.inputvalues['model'].modelname == 'ml'
