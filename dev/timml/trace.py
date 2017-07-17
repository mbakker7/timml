import numpy as np

def timtraceline(ml, xstart, ystart, zstart, hstepmax, nstepmax=100):
    direction = np.sign(hstepmax)
    aq = ml.aq.find_aquifer_data(xstart, ystart)
    assert zstart <= aq.z[0] and z >= aq.z[-1], "z value not inside aquifer"
    layer, ltype, modellayer = aq.findlayer(z)
    xyzt = np.array([[xstart, ystart, zstart, 0.0]])
    for i in range(nstepmax):
        x0, y0, z0 = xyzt[-1, :3]
        v0, vzbot, vztop = ml.velocitytrace(x0, y0, z0, aq, layer, ltype)
        if ltype == 'l':  # in leaky layer
            xnew, ynew = x0, y0
            if v0[2] * direction > 0:  # going up
                znew = aq.zlltop[layer]
                delt = (znew - ztrace[-1]) / abs(v0[2])
                modellayer -= 1  # layer above
            else:  # going down
                znew = aq.zllbot[layer]
                delt = (ztrace[-1] - znew) / abs(v0[2])
                modellayer += 1  # layer below
            tnew = ttrace[-1] + direction * delt
            layer = aq.layernumber[modellayer] 
            ltype = aq.ltype[modellayer]
        else:  # in aquifer
            tstep = hstepmax / np.sqrt(v0[0] ** 2 + v[1] ** 2)
            xtry, ytry, ztry = xyzt[-1, :3] + direction * tstep * v0
            layertry, ltypetry, modellayertry = aq.findlayer(ztry)
            if modellayertry < modellayer:  # stepping to layer above
                frac = (ml.aq.z[modellayer] - z0) / (ztry - z0)
                xtry, ytry, ztry = xyzt[-1, :3] + frac * direction * tstep * v0
            elif modellayertry > modellayer:  # stepping to layer below
                frac = (z0 - ml.aq.z[modellayer + 1]) / (ztry - z0)
                xtry, ytry, ztry = xyzt[-1, :3] + frac * direction * tstep * v0
            vtry, vzbot, vztop = ml.velocitytrace(xtry, ytry, ztry, aq, layer, ltype)
            vavg = 0.5 * (v0 + vtry)
            tstep = hstepmax / np.sqrt(vavg[0] ** 2 + vavg[1] ** 2)
            xnew, ynew, znew = xyzt[-1, :3] + direction * tstep * vavg
            layernew, ltypenew, modellayernew = aq.findlayer(znew)
            if modellayernew < modellayer:  # stepping to layer above
                frac = (ml.aq.z[modellayer] - z0) / (znew - z0)
                xnew, ynew, znew = xyzt[-1, :3] + frac * direction * tstep * vavg
                modellayernew = modellayer - 1
                layernew = aq.layernumber[modellayer]
                ltypenew = aq.ltype[modellayer] 
            elif modellayertry > modellayer:  # stepping to layer below
                frac = (z0 - ml.aq.z[modellayer + 1]) / (znew - z0)
                xnew, ynew, znew = xyzt[-1, :3] + frac * direction * tstep * vavg
            

                