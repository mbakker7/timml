import numpy as np


def param_maq(kaq, z, c, npor, top):
    # Computes the parameters for a ModelBase from input for a maq model
    kaq = np.atleast_1d(kaq).astype("d")
    z = np.atleast_1d(z).astype("d")
    c = np.atleast_1d(c).astype("d")
    npor = np.atleast_1d(npor).astype("d")
    if top == "conf":
        Naq = int(len(z) / 2)
        ltype = np.array(list((Naq - 1) * "al" + "a"))
    else:  # leaky layer on top
        Naq = int((len(z) - 1) / 2)
        ltype = np.array(list(Naq * "la"))
    if len(kaq) == 1:
        kaq = kaq * np.ones(Naq)
    assert len(kaq) == Naq, "Error: length of kaq needs to be 1 or" + str(Naq)
    H = z[:-1] - z[1:]
    assert np.all(H >= 0), "Error: Not all layers thicknesses are non-negative" + str(H)
    if top == "conf":
        if len(c) == 1:
            c = c * np.ones(Naq - 1)
        if len(npor) == 1:
            npor = npor * np.ones(2 * Naq - 1)
        assert len(c) == Naq - 1, "Error: Length of c needs to be 1 or" + str(Naq - 1)
        assert len(npor) == 2 * Naq - 1, "Error: Length of npor needs to be 1 or" + str(
            2 * Naq - 1
        )
        c = np.hstack((1e100, c))
    else:  # leaky layer on top
        if len(c) == 1:
            c = c * np.ones(Naq)
        if len(npor) == 1:
            npor = npor * np.ones(2 * Naq)
        assert len(c) == Naq, "Error: Length of c needs to be 1 or" + str(Naq)
        assert len(npor) == 2 * Naq, "Error: Length of npor needs to be 1 or" + str(
            2 * Naq
        )
    return kaq, c, npor, ltype


def param_3d(kaq, z, kzoverkh, npor, top="conf", topres=0):
    # Computes the parameters for a ModelBase from input for a model3d model
    kaq = np.atleast_1d(kaq).astype("d")
    z = np.atleast_1d(z).astype("d")
    kzoverkh = np.atleast_1d(kzoverkh).astype("d")
    npor = np.atleast_1d(npor).astype("d")
    if top == "conf":
        Naq = len(z) - 1
        ltype = np.array(Naq * ["a"])
    elif top == "semi":
        Naq = len(z) - 1
        ltype = np.hstack(("l", Naq * ["a"]))
    if len(kaq) == 1:
        kaq = kaq * np.ones(Naq)
    assert len(kaq) == Naq, "Error: length of kaq needs to be 1 or" + str(Naq)
    if len(kzoverkh) == 1:
        kzoverkh = kzoverkh * np.ones(Naq)
    assert len(kzoverkh) == Naq, "Error: length of kzoverkh needs to be 1 or" + str(Naq)
    if len(npor) == 1:
        if top == "conf":
            npor = npor * np.ones(Naq)
        elif top == "semi":
            npor = npor * np.ones(Naq + 1)
    if top == "conf":
        assert len(npor) == Naq, "Error: length of npor needs to be 1 or" + str(Naq)
    elif top == "semi":
        assert len(npor) == Naq + 1, "Error: length of npor needs to be 1 or" + str(
            Naq + 1
        )
    H = z[:-1] - z[1:]
    assert np.all(H >= 0), "Error: Not all layers thicknesses are non-negative" + str(H)
    c = 0.5 * H[:-1] / (kzoverkh[:-1] * kaq[:-1]) + 0.5 * H[1:] / (
        kzoverkh[1:] * kaq[1:]
    )
    if top == "conf":
        c = np.hstack((1e100, c))
    elif top == "semi":
        c = np.hstack((topres + 0.5 * H[0] / (kzoverkh[0] * kaq[0]), c))
    return kaq, c, npor, ltype
