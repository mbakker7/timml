import numpy as np

def timtraceline(ml, xstart, ystart, zstart, hstepmax, vstepfrac=0.2, tmax=1e12, nstepmax=100, win=[-1e30, 1e30, -1e30, 1e30], silent=False):
    # treating aquifer layers and leaky layers the same way
    xw1, xw2, yw1, yw2 = win
    terminate = False
    message = 'no message'
    eps = 1e-10  # used to place point just above or below aquifer top or bottom
    direction = np.sign(hstepmax)  # negative means backwards
    hstepmax = np.abs(hstepmax)
    aq = ml.aq.find_aquifer_data(xstart, ystart)
    if zstart > aq.z[0] or zstart < aq.z[-1]:
        terminate = True
        message = 'starting z value not inside aquifer'
    layer, ltype, modellayer = aq.findlayer(zstart)
    # slightly alter starting location not to get stuck in surpring points
    xyzt = [np.array([xstart * (1 + eps), ystart * (1 + eps), zstart, 0])]  # starting at time 0
    for i in range(nstepmax):
        if terminate: break
        x0, y0, z0, t0 = xyzt[-1]
        aq = ml.aq.find_aquifer_data(x0, y0)  # find new aquifer 
        layer, ltype, modellayer = aq.findlayer(z0)
        v0 = ml.velocomp(x0, y0, z0, aq, [layer, ltype]) * direction  # wordt nog gebruikt
        vx, vy, vz = v0
        if ltype == 'l':  # in leaky layer
            if vz > 0:  # upward through leaky layer
                if modellayer == 0:  # steps out of the top
                    z1 = aq.z[modellayer]
                    terminate = True
                else:
                    modellayer -= 1
                    z1 = aq.z[modellayer + 1] + eps * aq.Hlayer[modellayer]  # just above new bottom
            elif vz < 0:
                if modellayer == aq.Nlayers - 1:  # steps out of bottom
                    z1 = aq.z[modellayer + 1]
                    terminate = True
                else:
                    modellayer += 1
                    z1 = aq.z[modellayer] - eps * aq.Hlayer[modellayer]  # just below new top
            else:
                message = 'at point of zero leakage in leaky layer'
                terminate = True
                break
            t1 = t0 + abs((z1 - z0) / vz)
            xyztnew = [np.array([x0, y0, z1, t1])]
        else:  # in aquifer layer         
            vh = np.sqrt(vx ** 2 + vy ** 2)
            if vz > 0:  # flows upward
                if aq.z[modellayer] - z0 < vstepfrac * aq.Haq[layer]:
                    z1 = aq.z[modellayer] - eps * aq.Hlayer[modellayer]  # just below top
                else:
                    z1 = z0 + vstepfrac * aq.Haq[layer]
                tvstep = (z1 - z0) / vz
            elif vz < 0:
                if z0 - aq.z[modellayer + 1] < vstepfrac * aq.Haq[layer]:
                    z1 = aq.z[modellayer + 1] + eps * aq.Hlayer[modellayer]  # just above bot
                else:
                    z1 = z0 - vstepfrac * aq.Haq[layer]
                tvstep = (z0 - z1) / abs(vz)
            else:  # vz=0
                tvstep = np.inf
                z1 = z0
            if tvstep == np.inf and vh == 0:  # this should never happen anymore
                message = 'at point of zero velocity'
                terminate = True
                break
            if vh * tvstep > hstepmax:  # max horizonal step smaller than max vertical step
                thstep = hstepmax / vh
                z1 = z0 + thstep * vz
            else:
                thstep = tvstep
                # z1 is already computed
            x1 = x0 + thstep * vx
            y1 = y0 + thstep * vy
            t1 = t0 + thstep
            xyzt1 = np.array([x1, y1, z1, t1])
            # check if point needs to be changed
            correction = True
            for e in aq.elementlist:
                changed, terminate, xyztnew = e.changetrace(xyzt[-1], xyzt1, aq, layer, ltype, modellayer, direction)
                if changed or terminate:
                    correction = False
                    break
            if correction:  # correction step
                vx, vy, vz = 0.5 * (v0 + ml.velocomp(x1, y1, z1, aq, [layer, ltype]) * direction)
                vh = np.sqrt(vx ** 2 + vy ** 2)
                if vz > 0:  # flows upward
                    tvstep = min(aq.z[modellayer] - z0, vstepfrac * aq.Haq[layer]) / vz
                elif vz < 0:
                    tvstep = min(z0 - aq.z[modellayer + 1], vstepfrac * aq.Haq[layer]) / abs(vz)
                else:  # vz=0
                    tvstep = np.inf
                if vh * tvstep > hstepmax:  # max horizonal step smaller than vertical step
                    thstep = hstepmax / vh
                    x1 = x0 + thstep * vx
                    y1 = y0 + thstep * vy
                    z1 = z0 + thstep * vz        
                else:
                    thstep = tvstep
                    x1 = x0 + thstep * vx
                    y1 = y0 + thstep * vy
                    if vz > 0:  # flows upward
                        if aq.z[modellayer] - z0 < vstepfrac * aq.Haq[layer]:
                            if modellayer == 0:  # steps out of the top
                                z1 = aq.z[modellayer]
                                terminate = True
                                message = 'flowed out of top'
                            else:
                                modellayer -= 1
                                z1 = aq.z[modellayer + 1] + eps * aq.Hlayer[modellayer]  # just above new bottom
                        else:
                            z1 = z0 + thstep * vz
                    else:
                        if z0 - aq.z[modellayer + 1] < vstepfrac * aq.Haq[layer]:
                            if modellayer == aq.Nlayers - 1:  # steps out of bottom
                                z1 = aq.z[modellayer + 1]
                                terminate = True
                                message = 'flowed out of bottom'
                            else:
                                modellayer += 1
                                z1 = aq.z[modellayer] - eps * aq.Hlayer[modellayer]  # just below new top
                        else:
                            z1 = z0 + thstep * vz
                    if not terminate:
                        layer = aq.layernumber[modellayer]
                        ltype = aq.ltype[modellayer]
                t1 = t0 + thstep
                xyztnew = [np.array([x1, y1, z1, t1])]
                # check again if point needs to be changed
                for e in aq.elementlist:
                    changed, terminate, xyztchanged = e.changetrace(xyzt[-1], xyztnew[0], aq, layer, ltype, modellayer, direction)
                    if changed or terminate:
                        xyztnew = xyztchanged
                        break
        # check if outside window
        x1, y1, z1, t1 = xyztnew[0]
        frac = -1  # used to check later whether something changed
        if x1 < xw1:
            frac = abs((x0 - xw1) / (x1 - x0))
            x1, y1, z1, t1 = xyzt[-1] + frac * (xyztnew[0] - xyzt[-1])
            message = 'reached window boundary'
        if x1 > xw2:
            frac = abs((x0 - xw2) / (x1 - x0))
            x1, y1, z1, t1 = xyzt[-1] + frac * (xyztnew[0] - xyzt[-1])
            message = 'reached window boundary'
        if y1 < yw1:
            frac = abs((y0 - yw1) / (y1 - y0))
            x1, y1, z1, t1 = xyzt[-1] + frac * (xyztnew[0] - xyzt[-1])
            message = 'reached window boundary'
        if y1 > yw2:
            frac = abs((y0 - yw2) / (y1 - y0))
            x1, y1, z1, t1 = xyzt[-1] + frac * (xyztnew[0] - xyzt[-1])
            message = 'reached window boundary'
        if t1 > tmax:
            frac = abs((tmax - t0) / (t1 - t0))
            x1, y1, z1, t1 = xyzt[-1] + frac * (xyztnew[0] - xyzt[-1])
            message = 'reached tmax'
        if frac > 0:  # at least one of the above 5 ifs was true 
            terminate = True
            xyztnew = [np.array([x1, y1, z1, t1])]
        xyzt.extend(xyztnew)
        # check if element is reached
        for e in aq.elementlist:
            if e.stoptrace(xyzt[-1], layer, ltype, hstepmax, direction)[0]:
                terminate, xyztnew, message = e.stoptrace(xyzt[-1], layer, ltype, hstepmax, direction)
                xyzt.append(xyztnew)
    if not silent:
        print(message)
    return np.array(xyzt)

def timtracelines(ml, xstart, ystart, zstart, hstepmax, vstepfrac=0.2, tmax=1e12, nstepmax=100, silent='.', win=[-1e30, 1e30, -1e30, 1e30]):
    xyztlist = []
    for i in range(len(xstart)):
        xyztlist.append(timtraceline(ml, xstart[i], ystart[i], zstart[i], hstepmax=hstepmax, vstepfrac=vstepfrac, tmax=tmax, nstepmax=nstepmax, silent=silent, win=win))
        if silent == '.':
            print('.', end='', flush=True)
    if silent == '.':
        print('')
    return xyztlist

def crossline(xa, ya, xb, yb, z1, z2):
    eps = 1e-8
    za = xa + ya * 1j
    zb = xb + yb * 1j
    Za = (2 * za - (z1 + z2)) / (z2 - z1)
    Zb = (2 * zb - (z1 + z2)) / (z2 - z1)
    if Za.imag * Zb.imag < 0:
        Xa, Ya = Za.real, Za.imag
        Xb, Yb = Zb.real, Zb.imag
        X = Xa - Ya * (Xb - Xa) / (Yb - Ya)
        if abs(X) <= 1:
            Z = X + eps * np.sign(Yb) * 1j
            z = 0.5 * ((z2 - z1) * Z + z1 + z2)
            return True, z.real, z.imag
    return False
            
    

                