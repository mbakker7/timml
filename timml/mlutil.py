from matplotlib.pyplot import *
from matplotlib.colors import colorConverter
from matplotlib.collections import LineCollection
from numpy import *
from mltrace import *

def timcontour( ml, xmin, xmax, nx, ymin, ymax, ny, layers = 1, levels = 10, color = None, \
               width = 0.5, style = '-', separate = 0, layout = 1, newfig = 1, labels = 0, labelfmt = '%1.2f', xsec = 0,
               returnheads = 0, returncontours = 0, fillcontour = 0, size=(8,8), mayavi=False, verbose = True):
    '''Contours head with pylab'''
    rcParams['contour.negative_linestyle']='solid'

    rows = []   # Store every matrix in one long row

    # Determine aquifers to plot
    if layers == 'all':
        aquiferRange = range(ml.aq.Naquifers)
    if type(layers) == list:
        aquiferRange = layers
    else:
        if layers <= ml.aq.Naquifers:
            aquiferRange = range(layers)
        else:
            aquiferRange = range(ml.aq.Naquifers)
    Naquifers = len(aquiferRange)

    if mayavi:
        xg,yg = mgrid[xmin:xmax:nx*1j, ymin:ymax:ny*1j]
    else:
        xg,yg = meshgrid( linspace( xmin, xmax, nx ), linspace( ymin, ymax, ny ) )
        
    if verbose: print 'grid of '+str((nx,ny))+'. gridding in progress. hit ctrl-c to abort'
    
    Nrow,Ncol = shape(xg)

    head = zeros( ( Naquifers, ny, nx ), 'd' )
    for irow in range(Nrow):
        for jcol in range(Ncol):
            # Figure out aquifer first, so we can process fake aquifers; no repeated effort as aq is passed to headVector
            
            aq = ml.aq.findAquiferData(xg[irow,jcol], yg[irow,jcol])
            h = ml.headVector( xg[irow,jcol], yg[irow,jcol], aq )
            for k in range(Naquifers):
                if not aq.fakesemi:
                    if len(h) < aquiferRange[k]:
                        head[k,irow,jcol] = h[0]
                    else:
                        head[k,irow,jcol] = h[aquiferRange[k]]
                else:
                    head[k,irow,jcol] = h[aquiferRange[i]+1]
    # Contour
    # Manage colors
    if type( color ) is str:
        color = Naquifers * [color]
    elif type( color ) is list:
        Ncolor = len(color)      
        if Ncolor < Naquifers:
            color = color + Naquifers * [ color[0] ]
    elif type( color ) is type(None):
        color = ['b','r','g','m','c']
        if Naquifers > 5:
            color = int(ceil(Naquifers/5.)) * color
    # Manage line styles
    if type( style ) is str:
        style = Naquifers * [style]
    elif type( style ) is list:
        Nstyle = len(style)      
        if Nstyle < Naquifers:
            style = style + Naquifers * [ style[0] ]
    # Manage levels to draw
    if type(levels) is list:
        levdum = arange( levels[0],levels[1],levels[2] )
        levels = len(aquiferRange)*[levdum]
    elif levels == 'ask':
        levels = []
        for k in range(len(aquiferRange)):
            hmin = amin(head[k,:,:].flat)
            hmax = amax(head[k,:,:].flat)
            print 'Layer ',aquiferRange[k],' min,max: ',hmin,', ',hmax,'. Enter: hmin,hmax,step '
            h1,h2,delh = eval( raw_input() )
            levels = levels + [ arange(h1,h2+1e-8,delh) ]
    elif type(levels) is int:
        levels = len(aquiferRange)*[levels]
    elif isinstance(levels,ndarray):
        levels = len(aquiferRange)*[levels]
        
    # Drawing separate figures for each head
    xsize,ysize = size
    if separate:
        for k in range(len(aquiferRange)):
            figure( figsize=size )
            axis('scaled')
            axis( [xmin,xmax,ymin,ymax] )
            if layout: timlayout(ml,overlay=True,autoscale=False)
            if fillcontour:
                contourf( xg, yg, head[k,:,:], levels[k], colors = color[k] )
            else:
                contour( xg, yg, head[k,:,:], levels[k], colors = color[k], linewidths = width )
    else:
    # Drawing all heads on one figure
        if newfig:
            fig = figure( figsize=size )
            if xsec:
                #xc,yc = .1 + 0.5*0.8, .38 + 0.5*0.55
                #aspect = float((ymax-ymin)) / (xmax-xmin)
                #print 'aspect ',aspect
                #if aspect < .55 / .8:  # Wider than high
                #    dx,dy = .8, aspect * 0.8
                #    ax1 = axes([0.1,yc-0.5*dy,.8,dy])
                #else:
                #    dx,dy = 0.55 / aspect, 0.55  # Higher than wide
                #    print 'hello ',xc-0.5*dx,0.38,dx,.55
                #    ax1 = axes([xc-0.5*dx,0.38,dx,.55])
                ax1 = axes([.1,.38,.8,.55])
                #ax1.set_aspect(aspect='equal',adjustable='box-forced')
                ax1.set_aspect(aspect='equal')
                setp( ax1.get_xticklabels(), visible=False)
            else:
                ax1 = subplot(111)
                ax1.set_aspect(aspect='equal',adjustable='box')
        else:
            ax1 = gca()
        if layout: timlayout(ml,overlay=True,autoscale=False)
        for k in range(len(aquiferRange)):
            if fillcontour:
                contourset = ax1.contourf( xg, yg, head[k,:,:], levels[k], fmt=labelfmt )
            else:
                contourset = ax1.contour( xg, yg, head[k,:,:], levels[k], colors = color[k], linewidths = width, fmt=labelfmt )
                if verbose: print 'done with contouring'
                if style[k] != '-':
                    print 'mpl bug in setting line styles of collections'
##                for l in contourset.collections:
##                    if style[k] == '-':
##                        l.set_linestyle( (0, (1.0, 0.0)) )
##                    elif style[k] == '--':
##                        l.set_linestyle( (0, (6.0, 6.0)) )
##                    elif style[k] == ':':
##                        l.set_linestyle( (0, (1.0, 3.0)) )
##                    elif style[k] == '-.':
##                        l.set_linestyle( (0, (3.0, 5.0, 1.0, 5.0)) )
            if labels:
                clabel( contourset, inline = 1, fmt = labelfmt )
        ax1.set_xlim(xmin,xmax)
        ax1.set_ylim(ymin,ymax)
        if xsec:
            ax2 = axes([.1,.1,.8,.25],sharex=ax1)
            for i in range(ml.aq.Naquifers-1):
                fill( [xmin,xmax,xmax,xmin], [ml.aq.zt[i+1],ml.aq.zt[i+1],ml.aq.zb[i],ml.aq.zb[i]], fc=[.8,.8,.8],ec=[.8,.8,.8])
            ax2.set_xlim(xmin,xmax)
            ax2.set_ylim(ml.aq.zb[-1],ml.aq.zt[0])
        else:
            ax1.set_xlim(xmin,xmax)
            ax1.set_ylim(ymin,ymax)
    draw()
    #if not noshow: show()
    if not hasattr(ml,'trace'): ml.trace = TraceSettings(ml)
    ml.trace.xsec = xsec     
    if returnheads and returncontours:
        return head,contourset
    elif returnheads:
        return xg,yg,head
    elif returncontours:
        return contourset
    
def timcontourlocal(ml, nx=50, ny=50, Naquifers=1, levels = 10, color = None, \
               width = 0.5, style = '-', newfig = 1, labels = 0, labelfmt = '%1.2f', fillcontour = 0, xsec = False):
    ax = gcf().axes[0]
    x1,x2 = ax.get_xlim()
    y1,y2 = ax.get_ylim()
    if xsec & (not newfig):
        print 'Adding cross-section can only be done in new figure'
        newfig = 1
    layout = 0
    if newfig: layout = 1
    timcontour( ml, x1, x2, nx, y1, y2, ny, Naquifers=Naquifers, levels=levels, color=color, \
               width=width, style=style, separate=0, layout=layout, newfig=newfig, labels=labels, labelfmt=labelfmt, xsec=xsec,
               returnheads=0, returncontours=0, fillcontour=fillcontour)
    
class TraceSettings:
    def __init__(self,ml):
        self.forward = True
        self.tmax = 1e20
        self.zbegin = [0.5*(ml.aq.zb[0]+ml.aq.zt[0])]
        self.xsec = False
    
def setTrace(ml,forward=None,tmax=None,zbegin=None):
    if not hasattr(ml,'trace'): ml.trace = TraceSettings()
    if not forward is None:
        ml.trace.forward = forward
    if not tmax is None:
        ml.trace.tmax = tmax
    if not zbegin is None:
        if not iterable(zbegin): zbegin = [zbegin]
        ml.trace.zbegin = zbegin
        
def traceOn(ml,forward=None,tmax=None,zbegin=None):
    ml.ip = InteractivePathline(ml) # Needed so that ip remains after function is called; just ip.press isn't enough
    gcf().canvas.mpl_connect('button_press_event',ml.ip.press)
    setTrace(ml,forward,tmax,zbegin)

class InteractivePathline:
    def __init__(self,ml):
        self.ml = ml
    def press(self,event):
        #print 'Hello'
        if event.inaxes is None: return
        #print 'event.button ',event.button
        if event.button != 3: return
        ax = gcf().axes[0]
        x1,x2,y1,y2 = ax.axis()
        step = (x2 - x1) / 100.0
        if not self.ml.trace.forward: step = -step
        tmax = self.ml.trace.tmax
        zbegin = self.ml.trace.zbegin
        xsec = self.ml.trace.xsec
        Npoints = len(zbegin)
        timtracelines( self.ml, Npoints*[event.xdata], Npoints*[event.ydata], zbegin, step, tmax=tmax, Nmax=200, window=(x1,y1,x2,y2), xsec=xsec )
   

def timlayout( ml, color = 'k', overlay = 0, width = 0.5, style = 1, autoscale = True ):
    if not overlay:
        figure(figsize=(8,8))
        subplot(111)
    ax = gca()
    winx = ax.get_xlim()
    winy = ax.get_ylim()
    for e in ml.elementList:
        a = e.layout()
        nterms = len(a)
        for i in range(0,nterms,3):
            if a[i] > 1:
                if e.aquiferParent.fakesemi and e.pylayers[0] == 0: # In top layer (the fake layer) of fake semi
                    ax.plot( a[i+1], a[i+2], color = [.8,.8,.8] )
                else:
                    if style == 1:
                        ax.plot( a[i+1], a[i+2], color, linewidth = width )
                    elif style == 2:
                        ax.fill( a[i+1], a[i+2], facecolor = color, edgecolor = color )
            elif a[i] == 1:
                ax.plot( a[i+1], a[i+2], color+'o', markersize=3 ) 
    intensity = linspace(0.7,0.9,len(ml.aq.inhomList))
    for (col,inhom) in zip(intensity,ml.aq.inhomList):
        corners = inhom.layout()
        ax.fill( corners[0], corners[1], facecolor = cm.bone(col), edgecolor = [.8,.8,.8])
    if not autoscale:
        ax.set_xlim(winx)
        ax.set_ylim(winy)
        draw()
    if not overlay:
        axis('scaled')


def timtracelines(ml,xlist,ylist,zlist,step,twoD=1,tmax=1e30,Nmax=200,labfrac=2.0,
                    Hfrac=5.0,window=[-1e30,-1e30,1e30,1e30],overlay=1,color=None,
                    width=0.5,style='-',xsec=0,layout=True, verbose = True):
    '''Routine for plotting multiple tracelines using pylab'''
    # Set colors
    if type( color ) is str:
        color = ml.aq.Naquifers * [color]
    elif type( color ) is list:
        Ncolor = len(color)      
        if Ncolor < ml.aq.Naquifers:
            color = color + ml.aq.Naquifers * [ color[0] ]
    elif color is None:
        color = ['b','r','g','m','c']
        if ml.aq.Naquifers > 5:
            color = int(ceil(ml.aq.Naquifers/5.)) * color
    # Set figure
    if not overlay:
        fig = figure()
        ax1 = subplot(111)
    if overlay:
        fig = gcf()
        ax1 = fig.axes[0]
        if xsec:
            ax2 = fig.axes[1]
    xmin,xmax = ax1.get_xlim()
    ymin,ymax = ax1.get_ylim()
    for i in range(len(xlist)):
        x = xlist[i]; y = ylist[i]; z = zlist[i]
        [xyz,t,stop,pylayers] = traceline(ml,x,y,z,step,tmax,Nmax,labfrac=labfrac,Hfrac=Hfrac,window=window,verbose=verbose)
        pylayers = array(pylayers)
        if xsec:
            ax2.plot(xyz[:,0],xyz[:,2],color=[.7,.7,.7])
        for j in range(pylayers.min(),pylayers.max()+1):
            ax1.plot( where(pylayers==j,xyz[:,0],nan), where(pylayers==j,xyz[:,1],nan), color[j])
            if xsec:
                ax2.plot( where(pylayers==j,xyz[:,0],nan), where(pylayers==j,xyz[:,2],nan), color[j])
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(ymin,ymax)
    if xsec:
        ax2.set_ylim(ml.aq.zb[-1],ml.aq.zt[0])
        ax2.set_xlim(xmin,xmax)
    draw()
    return

def capturezone( ml, w, N, z, tmax, xsec=False ):
    xstart = w.xw + 1.01*w.rw * cos( arange(0.01,2*pi,2*pi/N) )
    ystart = w.yw + 1.01*w.rw * sin( arange(0.01,2*pi,2*pi/N) )
    zstart = z * ones(len(xstart))
    ax = gcf().axes[0]
    x1,x2,y1,y2 = ax.axis()
    step = (x2 - x1) / 100.0
    timtracelines(ml,xstart,ystart,zstart,-step,tmax=tmax,xsec=xsec)
    

def timvertcontour( ml, x1, y1, x2, y2, nx, zmin, zmax, nz, levels = 10, color = None, \
               width = 0.5, style = '-', newfig = 1, labels = 0, labelfmt = '%1.3f',
               returnheads = 0, returncontours = 0, fill = 0, size=None):
    #global CurrentModel
    #CurrentModel = ml
    # Compute grid
    xg = linspace(x1,x2,nx)
    yg = linspace(y1,y2,nx)
    zg = linspace(zmin,zmax,nz)
    horlength = sqrt( (x2-x1)**2 + (y2-y1)**2 )
    hor, vert = meshgrid( linspace( 0, horlength, nx ), zg )
    print 'grid of '+str((nx,nz))+'. gridding in progress. hit ctrl-c to abort'
    head = zeros( ( nz, nx ), 'd' )
    for irow in range(nz):
        for jcol in range(nx):
            # Figure out aquifer first, so we can process fake aquifers; no repeated effort as aq is passed to headVector
            head[irow,jcol] = ml.head3D( xg[jcol], yg[jcol], zg[irow] )
    # Contour
    # Manage levels to draw
    if type(levels) is list:
        levels = arange( levels[0],levels[2]+1e-8,levels[1] )
    elif levels == 'ask':
        hmin = min(head)
        hmax = max(head)
        print 'min,max: ',hmin,', ',hmax,'. Enter: hmin,step,hmax '
        input = raw_input(); input = string.split(input,',')
        levels = arange(float(eval(input[0])),float(eval(input[2]))+1e-8,float(eval(input[1])))
    # Drawing all heads on one figure
    if newfig:
        if size == None: size = (8,8)
        fig = figure( figsize=size )
    if fill:
        contourset = contourf( hor, vert, head[:,:], levels, linewidths = width, fmt=labelfmt )
    else:
        contourset = contour( hor, vert, head[:,:], levels, colors = color, linewidths = width, fmt=labelfmt )
        if style != '-':
            print 'mpl bug in setting line styles of collections'
##        for l in contourset.collections:
##            if style == '-':
##                l.set_linestyle( (0, (1.0, 0.0)) )
##            elif style == '--':
##                l.set_linestyle( (0, (6.0, 6.0)) )
##            elif style == ':':
##                l.set_linestyle( (0, (1.0, 3.0)) )
##            elif style == '-.':
##                l.set_linestyle( (0, (3.0, 5.0, 1.0, 5.0)) )
        if labels:
            clabel(contourset,fmt=labelfmt)
    draw()
    if returnheads and returncontours:
        return head,L
    elif returnheads:
        return xg,yg,head
    elif returncontours:
        return L
