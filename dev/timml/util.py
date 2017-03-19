import numpy as np
import matplotlib.pyplot as plt

def timcontour(ml, x1, x2, nx, y1, y2, ny, layers, levels, newfig=True):
    layers = np.atleast_1d(layers)
    xg = np.linspace(x1, x2, nx)
    yg = np.linspace(y1, y2, ny)
    h = ml.headgrid2(x1, x2, nx, y1, y2, ny, layers)
    if newfig: plt.figure()
    for ilayer in range(len(layers)):
        plt.contour(xg, yg, h[ilayer], levels)
    plt.show()