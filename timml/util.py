import numpy as np

from timml.controlpoints import controlpoints


def compute_z1z2(xy):
    # Returns z1 and z2 of polygon, in clockwise order
    x, y = list(zip(*xy, strict=False))
    if x[0] == x[-1] and y[0] == y[-1]:  # In case last point is repeated
        x = x[:-1]
        y = y[:-1]
    z1 = np.array(x) + np.array(y) * 1j
    index = list(range(1, len(z1))) + [0]
    z2 = z1[index]
    Z = 1e-6j
    z = Z * (z2[0] - z1[0]) / 2.0 + 0.5 * (z1[0] + z2[0])
    bigZ = (2.0 * z - (z1 + z2)) / (z2 - z1)
    bigZmin1 = bigZ - 1.0
    bigZplus1 = bigZ + 1.0
    angle = np.sum(np.log(bigZmin1 / bigZplus1).imag)
    if angle < np.pi:  # reverse order
        z1 = z1[::-1]
        z2 = z1[index]
    return z1, z2


def refine_n_segments(xy, shape_type, n_segments):
    """Refine line segments into n_segments each.

    Use cosine-rule to determine new segment lengths.

    Parameters
    ----------
    xy : list of tuple or np.array
        list of coordinates or 2d-array containing x-coordinates in the first column
        and y-coordinates in the second column
    shape_type : str
        shape type, either "line" or "polygon".
    n_segments : int
        number of segments to split each line segment into.

    Returns
    -------
    xynew : np.array
        array containing refined coordinates
    reindexer : np.array
        array containing index to original line segment, useful for obtaining element
        parameters from original input.
    """
    if shape_type == "polygon":
        z1, z2 = compute_z1z2(xy)
    elif shape_type == "line":
        z = xy[:, 0] + 1j * xy[:, 1]
        z1 = z[:-1]
        z2 = z[1:]
    else:
        raise ValueError("shptype must be one of 'polygon' or 'line'.")
    xpts = []
    ypts = []
    reindexer = []
    for i in range(len(z1)):
        xcpi, ycpi = controlpoints(n_segments - 1, z1[i], z2[i], include_ends=True)
        if i < len(z1) - 1:
            xpts.append(xcpi[:-1])
            ypts.append(ycpi[:-1])
        else:
            xpts.append(xcpi)
            ypts.append(ycpi)
        reindexer.append(i * np.ones(len(xcpi[:-1]), dtype=int))
    return (
        np.vstack([np.concatenate(xpts), np.concatenate(ypts)]).T,
        np.concatenate(reindexer),
    )
