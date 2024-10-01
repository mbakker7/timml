import numpy as np
import timml as tml


def modelmaq():
    # model parameters
    kh = 10  # m/day
    ctop = 1000.0  # resistance top leaky layer in days
    ztop = 0.0  # surface elevation
    zbot = -20.0  # bottom elevation of the model
    z = np.array([ztop + 1, ztop, -10, -10, zbot])
    ml = tml.ModelMaq(kaq=kh, z=z, c=[ctop, 1], topboundary="semi", hstar=0.0)
    return ml


def model3d():
    # model parameters
    kh = 10  # m/day
    kzoverkh = 0.25
    ctop = 1000.0  # resistance top leaky layer in days
    ztop = 0.0  # surface elevation
    zbot = -20.0  # bottom elevation of the model
    z = np.array([ztop, -10, zbot])
    ml = tml.Model3D(
        kaq=kh,
        kzoverkh=kzoverkh,
        z=z,
        topres=ctop,
        topthick=1.0,
        topboundary="semi",
        hstar=0.0,
    )
    return ml


def test_refine_n_segments_line():
    x1, x2, y1, y2 = -5, 5, 0, 0
    xy = np.array([(x1, y1), (x2, y2)])
    xyr, reindexer = tml.util.refine_n_segments(xy, "line", 3)
    assert np.allclose(xyr[:, 0], np.array([-5.0, -2.5, 2.5, 5.0]))
    assert (reindexer == 0).all() and len(reindexer) == 3


def test_refine_n_segments_polygon():
    xy = [
        (-10, -5),
        (10, -5),
        (10, 5),
        (-10, 5),
        (-10, -5),
    ]
    xyr, reindexer = tml.util.refine_n_segments(xy, "polygon", 2)
    assert np.all(reindexer == np.array([0, 0, 1, 1, 2, 2, 3, 3]))
    assert len(xyr) == 9


def test_refine_linesink():
    ml = modelmaq()
    tml.LineSinkBase(ml, refine_level=3)
    ml.solve(silent=True)
    assert np.allclose(ml.elementlist[-3].Qls, [25.0])
    assert np.allclose(ml.elementlist[-2].Qls, [50.0])
    assert np.allclose(ml.elementlist[-1].Qls, [25.0])
    assert len(ml.elementlist) == 4


def test_refine_headlinesink():
    ml = modelmaq()
    tml.HeadLineSink(ml, refine_level=2)
    ml.solve(silent=True)
    assert len(ml.elementlist) == 3


def test_refine_headlinesinkstring():
    ml = modelmaq()
    hls = tml.HeadLineSinkString(ml, refine_level=2)
    ml.solve(silent=True)
    assert len(hls.lslist) == 2
    assert ml.head(10, 10, layers=[0]) == 0.0
    assert np.sum(hls.discharge()) == 0.0


def test_refine_leakylinedoublet():
    ml = modelmaq()
    tml.LeakyLineDoublet(ml, res=100, refine_level=2)
    ml.solve(silent=True)
    assert len(ml.elementlist) == 3
    assert ml.head(10, 10, layers=[0]) == 0.0


def test_refine_leakylinedoubletstring():
    ml = modelmaq()
    llds = tml.LeakyLineDoubletString(ml, res=100, refine_level=2)
    ml.solve(silent=True)
    assert len(llds.ldlist) == 2
    assert ml.head(10, 10, layers=[0]) == 0.0


def test_refine_implinedoublet():
    ml = modelmaq()
    tml.ImpLineDoublet(ml, refine_level=2)
    ml.solve(silent=True)
    assert len(ml.elementlist) == 3
    assert ml.head(10, 10, layers=[0]) == 0.0


def test_refine_implinedoubletstring():
    ml = modelmaq()
    llds = tml.ImpLineDoubletString(ml, refine_level=2)
    ml.solve(silent=True)
    assert len(llds.ldlist) == 2
    assert ml.head(10, 10, layers=[0]) == 0.0


def test_refine_polygonimhommaq():
    ml = modelmaq()
    xy = [
        (-10, -5),
        (10, -5),
        (10, 5),
        (-10, 5),
        (-10, -5),
    ]
    inhom = tml.PolygonInhomMaq(
        ml,
        xy,
        kaq=ml.aq.kaq,
        z=ml.aq.z[1:],
        c=ml.aq.c[1:],
        topboundary="conf",
        refine_level=2,
    )
    tml.Well(ml, 0, 0)
    ml.solve(silent=True)
    eps = 1e-6
    xyin = [
        (-10 + eps, -5 + eps),
        (10 - eps, -5 + eps),
        (10 - eps, 5 - eps),
        (-10 + eps, 5 - eps),
        (-10 + eps, -5 + eps),
    ]
    assert len(ml.elementlist) == 19
    assert np.allclose(np.sum(ml.intnormflux(xyin, ndeg=99)), [100.0], rtol=1e-3)


def test_refine_polygonimhom3d():
    ml = model3d()
    xy = [
        (-10, -5),
        (10, -5),
        (10, 5),
        (-10, 5),
        (-10, -5),
    ]
    inhom = tml.PolygonInhom3D(
        ml,
        xy,
        kaq=ml.aq.kaq,
        kzoverkh=0.25,
        z=ml.aq.z[1:],
        topboundary="conf",
        refine_level=2,
    )
    tml.Well(ml, 0, 0)
    ml.solve(silent=True)
    eps = 1e-6
    xyin = [
        (-10 + eps, -5 + eps),
        (10 - eps, -5 + eps),
        (10 - eps, 5 - eps),
        (-10 + eps, 5 - eps),
        (-10 + eps, -5 + eps),
    ]
    assert len(ml.elementlist) == 19
    assert np.allclose(np.sum(ml.intnormflux(xyin, ndeg=99)), [100.0], rtol=1e-3)


def test_refine_buildingpitmaq():
    ml = modelmaq()
    xy = [
        (-10, -5),
        (10, -5),
        (10, 5),
        (-10, 5),
        (-10, -5),
    ]
    tml.BuildingPitMaq(
        ml,
        xy,
        kaq=ml.aq.kaq,
        z=ml.aq.z[1:],
        c=ml.aq.c[1:],
        topboundary="conf",
        refine_level=2,
    )
    tml.Well(ml, 0, 0)
    ml.solve(silent=True)
    eps = 1e-6
    xyin = [
        (-10 + eps, -5 + eps),
        (10 - eps, -5 + eps),
        (10 - eps, 5 - eps),
        (-10 + eps, 5 - eps),
        (-10 + eps, -5 + eps),
    ]
    assert len(ml.elementlist) == 35
    # accuracy of intnormflux around inner boundary is reasonable but not perfect
    assert np.allclose(
        np.sum(ml.intnormflux(xyin, ndeg=99), axis=1),
        [0.0, 100.0],
        atol=1e-1,
        rtol=1e-3,
    )


def test_refine_buildingpit3d():
    ml = model3d()
    xy = [
        (-10, -5),
        (10, -5),
        (10, 5),
        (-10, 5),
        (-10, -5),
    ]
    tml.BuildingPit3D(
        ml,
        xy,
        kaq=ml.aq.kaq,
        kzoverkh=0.25,
        z=ml.aq.z[1:],
        topboundary="conf",
        refine_level=2,
    )
    tml.Well(ml, 0, 0)
    ml.solve(silent=True)
    eps = 1e-6
    xyin = [
        (-10 + eps, -5 + eps),
        (10 - eps, -5 + eps),
        (10 - eps, 5 - eps),
        (-10 + eps, 5 - eps),
        (-10 + eps, -5 + eps),
    ]
    assert len(ml.elementlist) == 35
    # NOTE: accuracy of intnormflux around inner boundary isn't great...
    assert np.allclose(
        np.sum(ml.intnormflux(xyin, ndeg=99), axis=1),
        [0.0, 100.0],
        atol=1.0,
        rtol=1e-3,
    )


def test_refine_leakybuildingpitmaq():
    ml = modelmaq()
    xy = [
        (-10, -5),
        (10, -5),
        (10, 5),
        (-10, 5),
        (-10, -5),
    ]
    tml.LeakyBuildingPitMaq(
        ml,
        xy,
        kaq=ml.aq.kaq,
        z=ml.aq.z[1:],
        c=ml.aq.c[1:],
        res=[100, 100, 1, 100],
        topboundary="conf",
        refine_level=2,
    )
    tml.Well(ml, 0, 0)
    ml.solve(silent=True)
    eps = 1e-6
    xyin = [
        (-10 + eps, -5 + eps),
        (10 - eps, -5 + eps),
        (10 - eps, 5 - eps),
        (-10 + eps, 5 - eps),
        (-10 + eps, -5 + eps),
    ]
    assert len(ml.elementlist) == 35
    # accuracy of intnormflux around inner boundary is reasonable but not perfect
    assert np.allclose(np.sum(ml.intnormflux(xyin, ndeg=99)), [100.0], rtol=1e-3)


def test_refine_leakybuildingpit3d():
    ml = model3d()
    xy = [
        (-10, -5),
        (10, -5),
        (10, 5),
        (-10, 5),
        (-10, -5),
    ]
    tml.LeakyBuildingPit3D(
        ml,
        xy,
        kaq=ml.aq.kaq,
        kzoverkh=0.25,
        z=ml.aq.z[1:],
        res=[100, 100, 1, 100],
        topboundary="conf",
        refine_level=2,
    )
    tml.Well(ml, 0, 0)
    ml.solve(silent=True)
    eps = 1e-6
    xyin = [
        (-10 + eps, -5 + eps),
        (10 - eps, -5 + eps),
        (10 - eps, 5 - eps),
        (-10 + eps, 5 - eps),
        (-10 + eps, -5 + eps),
    ]
    assert len(ml.elementlist) == 35
    # accuracy of intnormflux around inner boundary is reasonable but not perfect
    assert np.allclose(np.sum(ml.intnormflux(xyin, ndeg=99)), [100.0], rtol=1e-3)


def test_global_refine_option():
    ml = modelmaq()
    tml.HeadLineSink(ml, refine_level=1)
    ml.solve(refine_level=3, silent=True)
    assert len(ml.elementlist) == 4


def test_multiple_solves():
    ml = modelmaq()
    tml.HeadLineSink(ml, refine_level=3)
    ml.solve(silent=True)
    assert len(ml.elementlist) == 4
    ml.solve(silent=True, refine_level=1)
    assert len(ml.elementlist) == 2


def test_reset_headlinesinkstring():
    ml = modelmaq()
    hls = tml.HeadLineSinkString(ml, refine_level=2)
    ml.initialize()
    assert len(hls.lslist) == 2
    ml.initialize(refine_level=1)
    assert len(hls.lslist) == 1


def test_reset_leakybuildingpitmaq():
    ml = modelmaq()
    xy = [
        (-10, -5),
        (10, -5),
        (10, 5),
        (-10, 5),
        (-10, -5),
    ]
    tml.LeakyBuildingPitMaq(
        ml,
        xy,
        kaq=ml.aq.kaq,
        z=ml.aq.z[1:],
        c=ml.aq.c[1:],
        res=[100, 100, 1, 100],
        topboundary="conf",
        refine_level=2,
    )
    tml.Well(ml, 0, 0)
    ml.initialize()
    assert len(ml.elementlist) == 35
    ml.initialize(refine_level=1)
    assert len(ml.elementlist) == 19
