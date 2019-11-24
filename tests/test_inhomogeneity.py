from timml import *
from numpy.testing import assert_allclose

class TestPolygonInhom():

    def test_xy_vertices(self):
        # model with vertices entered counter-clockwise
        ccxy = [(-1, -1), (1, -1), (1, 1), (-1, 1)]

        ml1 = ModelMaq(z=[1, 0, -1], c=[100], topboundary="semi", hstar=1.0)
        inhom = PolygonInhomMaq(ml1, ccxy)
        w = Well(ml1)
        ml1.solve()
        h1 = ml1.headgrid2(-2, 2, 10, -2, 2, 10)

        # model with vertices entered clockwise
        cwxy = [(-1, -1), (-1, 1), (1, 1), (1, -1)]

        ml2 = ModelMaq(z=[1, 0, -1], c=[100], topboundary="semi", hstar=1.0)
        inhom = PolygonInhomMaq(ml2, xy=cwxy)
        w = Well(ml2)
        ml2.solve()
        h2 = ml2.headgrid2(-2, 2, 10, -2, 2, 10)

        assert_allclose(h1, h2)

    def test_storeinput(self):

        xy = [(-1, -1), (-1, 1), (1, 1), (1, -1)]
        ml = ModelMaq()
        inhom = PolygonInhomMaq(ml, xy)

        assert inhom.inputvalues['kaq'] == 1
        assert inhom.inputvalues['npor'] == 0.3
        assert inhom.inputvalues['order'] == 3
        assert inhom.inputvalues['ndeg'] == 3
        assert inhom.inputvalues['hstar'] == None
        assert inhom.inputvalues['topboundary'] == 'conf'
        assert inhom.inputvalues['model'].modelname == 'ml'

class TestBuildingPit():

    def test_storeinput(self):
        
        xy = [(-1, -1), (-1, 1), (1, 1), (1, -1)]
        ml = ModelMaq()
        pit = BuildingPit(ml, xy)

        assert pit.inputvalues['kaq'] == 1
