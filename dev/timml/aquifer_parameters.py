import numpy as np

def param_maq(kaq, z, c, npor, top):
    # Computes the parameters for a ModelBase from input for a maq model
    kaq = np.atleast_1d(kaq).astype('d')
    z = np.atleast_1d(z).astype('d')
    c = np.atleast_1d(c).astype('d')
    npor = np.atleast_1d(npor).astype('d')
    if top == 'conf':
        Naq = len(z) / 2
        ltype = np.array( list( (Naq-1) * 'al' + 'a' ) )
    else: # leaky layer on top
        Naq = (len(z) - 1) / 2
        ltype = np.array( list( Naq * 'la' ) )
    if len(kaq) == 1: kaq = kaq * np.ones(Naq)
    assert len(kaq) == Naq, 'Error: length of kaq needs to be 1 or' + str(Naq)
    H = z[:-1] - z[1:]
    assert np.all(H >= 0), 'Error: Not all layers thicknesses are non-negative' + str(H) 
    if top == 'conf':
        if len(c) == 1: c = c * np.ones(Naq - 1)
        if len(npor) == 1: npor = npor * np.ones(2 * Naq - 1)
        assert len(c) == Naq-1, 'Error: Length of c needs to be 1 or' + str(Naq-1)
        assert len(npor) == 2 * Naq - 1, 'Error: Length of npor needs to be 1 or' + str(2*Naq-1)
        Haq = H[::2]
        c = np.hstack((1e100,c))
    else: # leaky layer on top
        if len(c) == 1: c = c * np.ones(Naq)
        if len(npor) == 1: npor = npor * np.ones(2 * Naq)
        assert len(c) == Naq, 'Error: Length of c needs to be 1 or' + str(Naq)
        assert len(npor) == 2 * Naq, 'Error: Length of npor needs to be 1 or' + str(2*Naq)
        Haq = H[1::2]
    return kaq, Haq, c, npor, ltype