from ml import *
from numpy import *
from mlaquifer import *
def traceline(ml,xstart,ystart,zstart,stepin,tmax=1e30,maxsteps=10,tstart=0.0,window=[-1e30,-1e30,1e30,1e30],\
              labfrac = 2.0, Hfrac = 2.0, verbose = True):
    '''Tracing routine
    Input:
    - xyzstart: XYZ object
    - step: 3D space step size
    - tmax: maximum time (will stop at tmax if no other element has been reached)
    - maxsteps: maximum number of steps (default = 10)
    Returns:
    - list[xyztrace,tt,[stopreason,stopelementlabel,totaltime],lt]'''

    # Set variable labfrac (should be input variable). Requirement is step < lab / labfrac
    labfrac = float(labfrac)
    # Set variable Hfrac (should be input variable). Requirement is vertical step < H / Hfrac
    Hfrac = float(Hfrac)   

    # Initialize variables
    stopreason = None; stopelement = None
    xyztrace = [(xstart,ystart,zstart)]; tt = [tstart]
    aq = ml.aq.findAquiferData(xstart,ystart)
    pyLayer = aq.inWhichPyLayer(zstart)
    lt = [ pyLayer ]
    tiny = 1e-12  # Used to figure out if we hit a stagnation point
    [xmin,ymin,xmax,ymax] = window

    # Check if starting inside aquifer system
    if pyLayer == 9999 or pyLayer == -9999:
        stopreason = 'Starting point outside aquifer'
        if verbose: print stopreason
        return [array(xyztrace),tt,[stopreason,stopelement,tt[-1]],lt]

    # Check if starting inside window
    if xstart < xmin or xstart > xmax or ystart < ymin or ystart > ymax :
        stopreason = 'Starting outside window'
        if verbose: print stopreason
        return [array(xyztrace),tt,[stopreason,stopelement,tt[-1]],lt]

    # Determine direction of trace
    if stepin > 0:  # Forward trace
        idir = 1.0
        stepsize = stepin
    else:
        idir = -1.0  # Backward trace
        stepsize = -stepin  # So stepsize is always a positive number !
    step = stepsize # Used first time through loop

    # Check if starting inside leaky layer
    if pyLayer < 0:  # Starting in leaky layer
        vxyz = ml.velocity2(xstart,ystart,zstart,aq,pyLayer)
        if ( vxyz[2] > 0 and idir > 0 ) or ( vxyz[2] < 0 and idir < 0 ) :  # Flowing up
            znew = aq.zb[-pyLayer-1]
            tnew = (znew-zstart)/vxyz[2]
            pyLayer = -pyLayer-1
        elif ( vxyz[2] < 0 and idir > 0 ) or ( vxyz[2] > 0 and idir < 0 ) : # Flowing down
            znew = aq.zt[-pyLayer]
            tnew = (zstart-znew)/(-vxyz[2])
            pyLayer = -pyLayer
        else:
            stopreason = 'Starting in leaky layer with zero vertical flow'
            return [array(xyztrace),tt,[stopreason,stopelement,tt[-1]],lt]
        xyztrace.append( (xstart, ystart, znew) )
        tt.append(tstart+idir*abs(tnew)); lt.append(pyLayer)

    # Store current variables as old
    xyzold = array([ xyztrace[-1][0], xyztrace[-1][1], xyztrace[-1][2] ])    
    told = tt[-1]
    aqOld = aq
    nstep = 0

    # Set some parameters
    Nelements = len(ml.elementList)
    NelementsRange = range(Nelements)

    # Compute list of distances to elements
    distanceArray = zeros(Nelements,'d')
    for i in NelementsRange:
        distanceArray[i] = ml.elementList[i].distanceSquaredToElement(xyzold[0],xyzold[1])
    distanceArray = sqrt(distanceArray)

### Start of tracing loop
    while nstep < maxsteps:

        changed = 0; stop = 0  # Variables to keep track if point location is changed due to element, or trace should be stopped

        # Compute new distance if less than zero
        for i in NelementsRange:
            if distanceArray[i] < step:
                distanceArray[i] = sqrt( ml.elementList[i].distanceSquaredToElement(xyzold[0],xyzold[1]) )

        # Check if stepsize is larger than lab_min / labfrac
        step = stepsize
        if aqOld.Naquifers > 1:
            if stepsize > aqOld.labsorted[-1] / labfrac:  # Then we have to check whether the stepsize should be reduced
                distancemin = min( distanceArray )
                if distancemin < aqOld.threelabsorted[0]:
                    for i in range(len(aqOld.lab)):
                        if distancemin < aqOld.threelabsorted[i]:
                            ilab = i
                        else:
                            break
                    step = aqOld.labsorted[ilab] / labfrac
                    if step > stepsize: step = stepsize  # For large lab, the above condition gives bigger step

        # Compute velocity
        vxyz = ml.velocity2( xyzold[0], xyzold[1], xyzold[2], aqOld, pyLayer )
        velo = sqrt( vxyz[0] * vxyz[0] + vxyz[1] * vxyz[1] + vxyz[2] * vxyz[2] )  # it may be quicker to return vxyz as list
        # Check if at stagnation point. For now, stop the trace
        if velo < tiny:
            if nstep == 0:
                stop = 1
                stopreason = 'Error: Starting at stagnation point. Move starting point'
                break
            else:
                stop = 1
                stopreason = 'reached stagnation point'
                break

        # Compute time step
        tstep = step / velo             

        # Check vertical step distance
        verticalstep = abs( tstep * vxyz[2] )
        if verticalstep > aqOld.H[pyLayer] / Hfrac:
            step = aqOld.H[pyLayer] / Hfrac / verticalstep * step
            tstep = step / velo

        # Compute candidate new point and aquifer
        xyznew = xyzold + vxyz * tstep * idir
        aqNew = ml.aq.findAquiferData( xyznew[0], xyznew[1] )

        # If new aquifer, step nicely into new aquifer
        backToOld = 0 # Else it doesn't get set when there is no aquifer change. This can be nicer
        if aqNew != aqOld:
            if aqNew != ml.aq:  # Thus inhomogeneity
                [changed, stop, xyznew, backToOld] = aqNew.crossBoundary(pyLayer,xyzold,xyznew,aqOld,step,idir)
            else:
                [changed, stop, xyznew, backToOld] = aqOld.crossBoundary(pyLayer,xyzold,xyznew,aqNew,step,idir)
            step = sqrt( sum( (xyznew-xyzold)**2 ) )
            tstep = step / velo
            if not backToOld:
                xyztrace.append(copy.copy(xyznew))
                tt.append( tt[-1] + tstep * idir )
                lt.append( pyLayer )
                changed = 1
            if backToOld:
                aqNew = aqOld
            else:
                if aqNew.fakesemi and not aqOld.fakesemi: pyLayer = pyLayer + 1
                if not aqNew.fakesemi and aqOld.fakesemi: pyLayer = pyLayer - 1
                if aqNew.Naquifers != aqOld.Naquifers: pyLayer = aqNew.inWhichPyLayer(xyznew[2])  # Layer number may have changed

        # If still in same aquifer, check if moving to new layer.
        # Also check if same number of aquifers on each side; assumed one to one connection
        # This was if, then it was elif, now it is back to elif. But we really do have to check whether we go to a new layer
        if ( aqNew == aqOld ) or backToOld: # Not sure why this used to be ( aqNew.Naquifers == aqOld.Naquifers ):

            if xyznew[2] < aqNew.zb[pyLayer]:  # Stepping down
                [xyztrace,tt,lt,xyznew,pyLayer,stop] = stepDownLayer(ml,aqOld,pyLayer,xyzold,xyznew,tstep,idir,xyztrace,tt,lt)
                changed = 1
                if stop:
                    stopreason = 'flowed out of bottom of aquifer system'
                    break                
            elif xyznew[2] > aqNew.zt[pyLayer]:  # Stepping up
                [xyztrace,tt,lt,xyznew,pyLayer,stop] = stepUpLayer(ml,aqOld,pyLayer,xyzold,xyznew,tstep,idir,xyztrace,tt,lt)
                changed = 1
                if stop:
                    stopreason = 'flowed out of top of aquifer system'
                    break
            elif backToOld:  # Need to add point, and continue
                xyztrace.append(copy.copy(xyznew))
                tt.append( tt[-1] + tstep * idir )
                lt.append( pyLayer )
                changed = 1

        # Need to fix when jumps to new layer OVER inhomogeneity boundary

        if changed == 0:  # Staying within same layer (it is assumed that the check of elements will take care of new aquifer                
            # Check if it reaches or crosses an element
            for i in NelementsRange:
                echanged = 0
                if distanceArray[i] < step:  # Then element is within reach
                    [echanged, stop, xyzchanged] = ml.elementList[i].nearElement(pyLayer,xyzold,xyznew,step,idir)
                if echanged:  
                    xyznew = xyzchanged
                    step = sqrt( sum( (xyznew-xyzold)**2 ) )
                    tstep = step / velo
                    xyztrace.append(copy.copy(xyznew))
                    tt.append( tt[-1] + tstep * idir )
                    lt.append( pyLayer )
                    changed = 1
                    # Used to allow for only one change, but not anymore. Not sure if this always works
                if stop:
                    stopreason = 'reached element of type ' + ml.elementList[i].type
                    stopelement = ml.elementList[i].label
                    break
            if changed == 1:  # Compute new aquifer and layer
                aqNew = ml.aq.findAquiferData(xyznew[0],xyznew[1])
                pyLayer = aqNew.inWhichPyLayer(xyznew[2])
        if stop: break

    ### Correction step
        if changed == 0:  # No change was made, do correction step

            # Compute new velocity and stepsize
            vxyznew = ml.velocity2(xyznew[0],xyznew[1],xyznew[2],aqNew,pyLayer)
            vxyz = 0.5 * (vxyz + vxyznew)  # New velocity is average of the two
            velo = sqrt( vxyz[0] * vxyz[0] + vxyz[1] * vxyz[1] + vxyz[2] * vxyz[2] )

            # Check if at stagnation point. For now, stop the trace
            if velo < tiny:
                stop = 1
                stopreason = 'reached stagnation point'
                break

            # Compute time step
            tstep = step / velo
            
            # Compute new candidate new point
            xyznew = xyzold + vxyz * tstep * idir

            # Check again if still in same aquifer
            aqNew = ml.aq.findAquiferData(xyznew[0],xyznew[1])

            # If new aquifer, step nicely into new aquifer
            if aqNew != aqOld:
                if aqNew != ml.aq:  # Thus inhomogeneity
                    [changed, stop, xyznew, backToOld] = aqNew.crossBoundary(pyLayer,xyzold,xyznew,aqOld,step,idir)
                else:
                    [changed, stop, xyznew, backToOld] = aqOld.crossBoundary(pyLayer,xyzold,xyznew,aqNew,step,idir)
                step = sqrt( sum( (xyznew-xyzold)**2 ) )
                tstep = step / velo
                xyztrace.append(copy.copy(xyznew))
                tt.append( tt[-1] + tstep * idir )
                lt.append( pyLayer )
                changed = 1
                if backToOld:
                    aqNew = aqOld
                else:
                    if aqNew.fakesemi and not aqOld.fakesemi: pyLayer = pyLayer + 1
                    if not aqNew.fakesemi and aqOld.fakesemi: pyLayer = pyLayer - 1
                    if aqNew.Naquifers != aqOld.Naquifers: pyLayer = aqNew.inWhichPyLayer(xyznew[2])  # Layer number may have changed

            # If still in same aquifer, check if moving to new layer.
            # Also check if same number of aquifers on each side; assumed one to one connection
            # Made into elif
            elif ( aqNew == aqOld ) or ( aqNew.Naquifers == aqOld.Naquifers ):
                
                if xyznew[2] < aqNew.zb[pyLayer]:  # Stepping down
                    [xyztrace,tt,lt,xyznew,pyLayer,stop] = stepDownLayer(ml,aqOld,pyLayer,xyzold,xyznew,tstep,idir,xyztrace,tt,lt)
                    changed = 1
                    if stop:
                        stopreason = 'flowed out of bottom of aquifer system'
                        break                
                elif xyznew[2] > aqNew.zt[pyLayer]:  # Stepping up
                    [xyztrace,tt,lt,xyznew,pyLayer,stop] = stepUpLayer(ml,aqOld,pyLayer,xyzold,xyznew,tstep,idir,xyztrace,tt,lt)
                    changed = 1
                    if stop:
                        stopreason = 'flowed out of top of aquifer system'
                        break

            if changed == 0:  # Staying within layer                 
                # Check if it reaches or crosses an element
                for i in NelementsRange:
                    echanged = 0
                    if distanceArray[i] < step:
                        [echanged, stop, xyzchanged] = ml.elementList[i].nearElement(pyLayer,xyzold,xyznew,step,idir)
                    if echanged:  
                        xyznew = xyzchanged
                        step = sqrt( sum( (xyznew-xyzold)**2 ) )
                        tstep = step / velo
                        xyztrace.append(copy.copy(xyznew))
                        tt.append( tt[-1] + tstep * idir )
                        lt.append( pyLayer )
                        changed = 1
                        # Used to allow for only one change, but not anymore. Not sure if this always works
                    if stop:
                        stopreason = 'reached element of type ' + ml.elementList[i].type
                        stopelement = ml.elementList[i].label
                        break
                if changed == 1:  # Compute new aquifer and layer
                    aqNew = ml.aq.findAquiferData(xyznew[0],xyznew[1])
                    pyLayer = aqNew.inWhichPyLayer(xyznew[2])
            if stop: break
            
    ### End of correction step

        if changed == 0:
            xyztrace.append(copy.copy(xyznew))
            tt.append( tt[-1] + tstep * idir )
            lt.append( pyLayer )

        # Check if time expired
        if abs(tt[-1]) > abs(tmax):
            if changed != 0:  # Could have stepped to other layer
                if abs(tt[-2]) > abs(tmax):  # Time at beginning of leaky layer was already too large
                    xyztrace.pop(-1)
                    tt.pop(-1)
                    lt.pop(-1)
            tstepold = abs(tt[-1]) - abs(tt[-2])
            tstepnew = abs(tmax) - abs(tt[-2])
            xyztrace[-1] = xyztrace[-2] + tstepnew/tstepold * ( xyztrace[-1] - xyztrace[-2] )
            tt[-1] = idir * tmax
            stop = 1
            stopreason = 'reached tmax'
            break  # Nicelye stops at tmax; assumes smooth transition between last two points
        
        # Check if new point inside window (Should be modified to nicely stop at window boundary!)
        if xyznew[0] < xmin or xyznew[0] > xmax or xyznew[1] < ymin or xyznew[1] > ymax :
            if xyznew[1] > ymax:
                frac = ( ymax - xyztrace[-2][1] ) / ( xyztrace[-1][1] - xyztrace[-2][1] )
                xyztrace[-1] = xyztrace[-2] + frac * ( xyztrace[-1] - xyztrace[-2] )
                tt[-1] = tt[-2] + frac * ( tt[-1] - tt[-2] ) 
            stop = 1
            stopreason = 'reached window boundary'
            break

        nstep = nstep + 1
        xyzold = xyznew
        aqOld = aqNew
        distanceArray = distanceArray - step

    if not stop: stopreason = 'reached maximum number of steps'
    if verbose: print stopreason

    return [array(xyztrace),tt,[stopreason,stopelement,tt[-1]],lt]

def stepDownLayer(ml,aq,pyLayer,xyzold,xyznew,tstep,idir,xyztrace,tt,lt):
    stop = 0
    if pyLayer == aq.Naquifers-1:  
        # Occurs when step is too large and velocity very vertical; take step to halfway bottom
        frac = 0.5 * ( xyzold[2] - aq.zb[pyLayer] ) / ( xyzold[2] - xyznew[2] )
        xyznew = xyzold + frac * ( xyznew - xyzold )
        xyztrace.append(copy.copy(xyznew)); tt.append( tt[-1] + tstep * idir * frac ); lt.append( pyLayer )
    else:
        frac = ( xyzold[2] - aq.zb[pyLayer] ) / ( xyzold[2] - xyznew[2] )
        xyznew = xyzold + frac * ( xyznew - xyzold )
        xyztrace.append(copy.copy(xyznew)); tt.append( tt[-1] + tstep * idir * frac ); lt.append( pyLayer )
        # Put next point at top of layer below
        xyznew[2] = aq.zt[pyLayer+1]
        xyztrace.append(copy.copy(xyznew))
        if aq.HLeakyLayer[pyLayer+1] > 0:  # If leaky layer thickness > 0, add travel time through leaky layer
            vxyzLeakyLayer = ml.velocity2( xyznew[0], xyznew[1], 0.0, aq, -(pyLayer+1) )
            extratime = aq.HLeakyLayer[pyLayer+1] / abs( vxyzLeakyLayer[2] )
            tt.append( tt[-1] + idir * extratime )  # Gotta check for tmax here
        else:
            tt.append( tt[-1] )
        pyLayer = pyLayer + 1  # Force new layer to be the new layer
        lt.append( pyLayer )
    return [xyztrace,tt,lt,xyznew,pyLayer,stop]

def stepUpLayer(ml,aq,pyLayer,xyzold,xyznew,tstep,idir,xyztrace,tt,lt):
    stop = 0
    # Put next point at bottom of layer above
    if ( not aq.fakesemi and pyLayer > 0 ) or ( aq.fakesemi and pyLayer > 1 ):  # Then there is an aquifer above
        frac = ( aq.zt[pyLayer] - xyzold[2] ) / ( xyznew[2] - xyzold[2] )
        xyznew = xyzold + frac * ( xyznew - xyzold )
        xyztrace.append(copy.copy(xyznew)); tt.append( tt[-1] + tstep * idir * frac ); lt.append(pyLayer)
        xyznew[2] = aq.zb[pyLayer-1]
        xyztrace.append(copy.copy(xyznew))
        if aq.HLeakyLayer[pyLayer] > 0:  # If leaky layer thickness > 0, add travel time through leaky layer
            vxyzLeakyLayer = ml.velocity2( xyznew[0], xyznew[1], 0.0, aq, -pyLayer )
            extratime = aq.HLeakyLayer[pyLayer] / abs( vxyzLeakyLayer[2] )
            tt.append( tt[-1] + idir * extratime )
        else:
            tt.append( tt[-1] )
        pyLayer = pyLayer - 1  # Force new layer to be the new layer
        lt.append(pyLayer)
    else:  # Flowing out of top. Should check if this is legit
        if aq.fakesemi:
            itop = 1
        else:
            itop = 0
        vxyz1 = ml.velocity2(xyzold[0],xyzold[1],xyzold[2],aq,pyLayer)
        vxyz2 = ml.velocity2(xyzold[0],xyzold[1],aq.zt[itop],aq,pyLayer)
        if ( vxyz2[2] <= 0 and idir > 0 ) or ( vxyz2[2] >= 0 and idir < 0 ):  # Take smaller step else it goes out of top
            frac = 0.5 * ( aq.zt[itop] - xyzold[2] ) / ( xyznew[2] - xyzold[2] )
            xyznew = xyzold + frac * ( xyznew - xyzold )
            xyztrace.append(copy.copy(xyznew)); tt.append( tt[-1] + tstep * idir * frac ); lt.append( pyLayer )
        elif abs( vxyz2[2] - vxyz1[2] ) > 0.001:  # Then too big a jump in velocity; half distance to top; still experimental
            frac = 0.5 * ( aq.zt[itop] - xyzold[2] ) / ( xyznew[2] - xyzold[2] )
            xyznew = xyzold + frac * ( xyznew - xyzold )
            xyztrace.append(copy.copy(xyznew)); tt.append( tt[-1] + tstep * idir * frac ); lt.append( pyLayer )
        else:  # I guess it really goes out
            frac = ( aq.zt[itop] - xyzold[2] ) / ( xyznew[2] - xyzold[2] )
            xyznew = xyzold + frac * ( xyznew - xyzold )
            xyztrace.append(copy.copy(xyznew)); tt.append( tt[-1] + tstep * idir * frac ); lt.append( pyLayer )
            if aq.fakesemi:
                xyznew[2] = aq.zb[0]
                xyztrace.append(copy.copy(xyznew))
                vxyzLeakyLayer = ml.velocity2( xyznew[0], xyznew[1], 0.0, aq, -pyLayer )
                extratime = aq.HLeakyLayer[pyLayer] / abs( vxyzLeakyLayer[2] )
                tt.append( tt[-1] + idir * extratime ); lt.append( pyLayer )
            stop = 1
    return [xyztrace,tt,lt,xyznew,pyLayer,stop]

def tracewrite(self,trace):
    '''Routine to write trace (TraceLine instance) to file'''
    out = open('/temp/trace.bln','w')
    out.write(str(len(trace.xyzt))+' 1\n')
    for i in range(len(trace.xyzt)):
        out.write(str(trace.xyzt[i].x)+' '+str(trace.xyzt[i].y)+'\n')
    out.close()

def matlabtracelines(ml,xrange,yrange,zrange,step,filename,twoD=1,tmax=1e30,Nmax=20,labfrac=2.0,Hfrac=5.0,window=[-1e30,-1e30,1e30,1e30]):
    '''Routine for calculating multiple tracelines and writing them to file.
    xyz is list of XYZ objects, filename is character string, such as 'mark.dat'
    and must be between quotes
    '''
    out = open(filename,'w')
    tend = []
    for i in range(len(xrange)):
        x = xrange[i]; y = yrange[i]; z = zrange[i]
        [xyz,t,stop,lt] = traceline(ml,x,y,z,step,tmax,Nmax,labfrac=labfrac,Hfrac=Hfrac,window=window)
        tend = tend + [ t[-1] ]
        out.write('xyz'+str(i)+'=[')
        out.write( str(xyz[0,0])+','+str(xyz[0,1])+','+str(xyz[0,2]) )
        for j in range(1,xyz.shape[0]):
            out.write( ';'+str(xyz[j,0])+','+str(xyz[j,1])+','+str(xyz[j,2]) )
        out.write('];\n')
        if twoD==1:
            out.write( 'plot(xyz'+str(i)+'(:,1),xyz'+str(i)+'(:,2))\n' )
        elif twoD==2:
            out.write( 'plot(xyz'+str(i)+'(:,1),xyz'+str(i)+'(:,3))\n' )
        else:
            out.write( 'plot3(xyz'+str(i)+'(:,1),xyz'+str(i)+'(:,2),xyz'+str(i)+'(:,3))\n' )
    out.close()
    return tend

def surfertrace(ml,xstart,ystart,zstart,stepin,tmax=1e30,maxsteps=10,tstart=0.0,filename='/temp/dump.bln'):
    '''filename between quotes with .bln extension. Layout of all elements, independent of layer'''
    [xyz,t,stop,lt]=traceline(ml,xstart,ystart,zstart,stepin,tmax,maxsteps,tstart)
    out = open(filename,'w')
    out.write(str(len(xyz))+' 1\n')
    for p in xyz:
        out.write(str(p[0])+' '+str(p[1])+' '+str(p[2])+'\n')
    out.close

def matlabtracelines2(ml,xrange,yrange,zrange,step,filename,twoD=1,tmax=1e30,Nmax=20):
    '''Routine for calculating multiple tracelines and writing them to file.
    xyz is list of XYZ objects, filename is character string, such as 'mark.dat'
    and must be between quotes
    '''
    out = open(filename,'w')
    for i in range(len(xrange)):
        print (str(i)+' of '+str(len(xrange)))
        x = xrange[i]; y = yrange[i]; z = zrange[i]
#        try:        
        [xyz,t,stop,lt] = traceline(ml,x,y,z,step,tmax,Nmax)
#        except:
#           print 'Problem with point',x,y,z
#           out.write('xyz'+str(i)+'=['+str(x)+','+str(y)+','+str(z)+','+str(0))
#           out.write('];\n')
#           out.write('stopreason'+str(i)+'=')
#           out.write("'")
#           out.write("error in python")
#           out.write("';")        
#           out.write('\n')
#           continue
        out.write('xyz'+str(i)+'=[')
        out.write( str(xyz[0,0])+','+str(xyz[0,1])+','+str(xyz[0,2])+','+str(t[0]) )
        for j in range(1,xyz.shape[0]):
            out.write('\n')
            out.write( str(xyz[j,0])+','+str(xyz[j,1])+','+str(xyz[j,2])+','+str(t[j]) )
        out.write('];\n')
        
        out.write('stopreason'+str(i)+'=')
        out.write("'")
        out.write(stop[0])
        out.write("';")        
        out.write('\n')
        
        
        if twoD==1:
            out.write( 'plot(xyz'+str(i)+'(:,1),xyz'+str(i)+'(:,2))\n' )
        elif twoD==2:
            out.write( 'plot(xyz'+str(i)+'(:,1),xyz'+str(i)+'(:,3))\n' )
        else:
            out.write( '\n' )
    out.close()

