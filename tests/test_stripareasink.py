from timml import *

class TestStripAreaSink():

    def test_storeinput(self):

        ml = ModelMaq()
        s = StripAreaSink(ml)

        assert s.inputvalues['model'].modelname == 'ml'


class TestStripAreaSinkInhom():

    def test_storeinput(self):

        ml = ModelMaq()
        s = StripAreaSinkInhom(ml)

        assert s.inputvalues['model'].modelname == 'ml'
