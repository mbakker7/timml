from mllinesink import *
from mlaquifer import *
from mlinhom import *
from mlpolyareasink import *
from shapely.geometry import Polygon as Polygonsh
from matplotlib.pyplot import *

class Lake:
    
    def __init__(self,ml,xylist,hlake,cbot,Hbot,nbot=0.3,Lmax=None,areatol=0.2):
        #assert order <= 2, 'TimML Input Error: Maximum order of a Lake is 2'
        self.ml = ml
        self.ml.lakeList.append(self)
        if Lmax == None:
            self.xylist = xylist
        else:
            self.xylist = [xylist[0]]
            xylist = xylist + [xylist[0]] # Add first point at end
            for i in range(1,len(xylist)):
                x0,y0 = xylist[i-1]; x1,y1 = xylist[i]
                L = sqrt( (x0-x1)**2 + (y0-y1)**2 )
                if L > Lmax:
                    Np = int(L/Lmax) + 2
                    xp = linspace( x0, x1, Np )
                    yp = linspace( y0, y1, Np )
                else:
                    xp = [x0,x1]; yp = [y0,y1]
                for xy in zip(xp[1:],yp[1:]):
                    self.xylist.append(xy)
            self.xylist = self.xylist[:-1] # Remove point at end
        self.hlake = float(hlake)
        self.cbot = float(cbot)
        self.Hbot = float(Hbot)
        self.nbot = float(nbot)
        self.order = 2  # Hardcoded until short higher-order line elements are operational
        self.areatol = areatol
        self.polygon = Polygonsh(self.xylist)
        # xmin,xmax used for gridding when checking for percolation
        x,y = zip(*self.xylist)
        self.xmin = min(x); self.xmax = max(x)
        self.ymin = min(y); self.ymax = max(y)
        self.inhomList = []
        self.elementList = []
        self.areasinkList = []; self.areasinkarea = 0.0
        self.setup(self.xylist)

    def setup(self,xylist):
        k = [100*self.ml.aq.Tcomp] + self.ml.aq.k.tolist()
        zb = [self.ml.aq.zt[0]+self.Hbot] + self.ml.aq.zb.tolist()
        zt = [zb[0]+1.0] + self.ml.aq.zt.tolist()
        c = self.ml.aq.c.tolist()
        if len(c) == 2:  # Only one aquifer
            c = [self.cbot]
        else:
            c.insert(1,self.cbot)
            c = c[1:-1]  # Strip off first and last value
        n = [self.nbot] + self.ml.aq.n.tolist()
        nll = [self.nbot] + self.ml.aq.nll.tolist()
        self.inhomList.append( PolygonInhom(self.ml,k,zb,zt,c,xylist,n,nll,True) )
        self.elementList.extend( MakeInhomPolySide(self.ml, xylist, self.order, closed = True) )
        # Put HeadLineSink along first straight part of boundary (even if broken-up in multiple segments)
        xylist.append( xylist[0] ) # add first point again at the end
        istart = 0
        x0,y0 = xylist[0]; x1,y1 = xylist[1]; an0 = arctan2((y1-y0),(x1-x0))
        for i in range(istart,len(xylist)-1):
            x0,y0 = xylist[i]; x1,y1 = xylist[i+1]; an = arctan2((y1-y0),(x1-x0))
            if an != an0:
                x0,y0 = xylist[istart]; x1,y1 = xylist[i]
                self.elementList.append( HeadLineSink(self.ml,x0,y0,x1,y1,self.hlake,[1],self.inhomList[-1]) )
                istart = i
                an0 = an
        # Add HeadLineSink along last segment
        x0,y0 = xylist[istart]; x1,y1 = xylist[i+1]
        self.elementList.append( HeadLineSink(self.ml,x0,y0,x1,y1,self.hlake,[1],self.inhomList[-1]) )
        ## Place head line-sink exactly inside. Nicer would be to put the control point just inside, but that is not an option in the current code
        #z0 = x0 + 1j * y0; z1 = x1 + 1j * y1
        #Z0 = complex( -1+1e-6, 1e-6 ); Z1 = complex( 1-1e-6, 1e-6 )
        #z0new = 0.5 * Z0 * (z1-z0) + 0.5*(z0+z1); z1new = 0.5 * Z1 * (z1-z0) + 0.5*(z0+z1)
        #x0 = z0new.real; y0 = z0new.imag; x1 = z1new.real; y1 = z1new.imag        
        assert self.elementList[0].aqout == self.ml.aq, 'TimML Input Error: Outside of Lake must be background aquifer of model'

    def check_percolate(self,hmin):
        # Number of points is the bigger of 5 or the number of lambda's in x or y direction
        if len(self.inhomList) > 0:
            lab = self.inhomList[0].lab[0]
        else:
            lab = 1e6 # Only used to determine nx,ny
        nx = max( [5, int( (self.xmax-self.xmin)/lab )] )
        ny = max( [5, int( (self.ymax-self.ymin)/lab )] )
        # Make grid in top layer of aquifer
        x,y = meshgrid(linspace(self.xmin,self.xmax,nx),linspace(self.ymin,self.ymax,ny))
        f = vectorize(self.ml.head)
        h = f(0,x,y)  # fixed to base zero
        # Set interactive false, draw contours, set interactive back to what it was
        interactive = rcParams['interactive']
        rcParams['interactive'] = False
        cobj = contourf(x,y,h,[hmin-1e6,hmin])
        clf() # Clear current figure
        rcParams['interactive'] = interactive
        # We should only get one collection of contour objects, as we only contour one level
        if len(cobj.collections) > 1: print 'cobj.collections has more than 1 item'
        a = cobj.collections[0]
        b = a.get_paths()
        polygonlist = []
        for bb in b:
            xy = bb.vertices
            polygonlist.append( Polygonsh(xy.tolist()) )
        # areasink is intersection of lake with contour level
        # newlake is difference between lake and new areasink
        areasinklist = []
        for p in polygonlist:
            areasink = self.polygon.intersection(p)
            if ~areasink.is_empty:
                areasinklist.append( areasink )
        if len(areasinklist) == 0:
            areasink = None
            newlake = None
        else:
            newlake = self.polygon
            for p in polygonlist:
                newlake = newlake.difference(p)
            if newlake.is_empty: newlake = None
        return areasinklist,newlake

    def change_lake(self):
        areasinklist,newlake = self.check_percolate(self.ml.aq.zt[0])
        if len(areasinklist) == 0:  # No percolating part
            if len(self.areasinkList) > 0:  # Need to remove areasink of previous iteration (won't happen much I guess)
                print "Lake doesn't percolate"
                self.remove_elements()
                self.setup( self.xylist )
                return True
            else:  # No change from previous iteration
                return False 
        elif newlake == None:  # Entire lake percolates
            if len(self.elementList) > 0:  # Remove lake and replace by areasink
                print 'Entire lake percolates'
                self.remove_elements()
                self.create_areasink(areasinklist[0])
                return True
            else:  # No change from previous iteration
                return False
        else:
            for a in areasinklist: # Check if any areasink entirely inside lake (cannot do embedded inhoms yet)
                if not a.touches(newlake):
                    print 'Error: percolating area-sink entirely inside lake'
                    return False  # No change is made
        if len(self.areasinkList) > 0:
            area = 0.0
            for a in areasinklist:
                area += a.area
            change = abs( area - self.areasinkarea ) / self.areasinkarea           
            if change < self.areatol:
                print 'Relative change in percolation area: ',change,' within tolerance of ',self.areatol
                return False  # No change required
            else:
                print 'Relative change in percolation area: ',change
        else:
            print 'Part of lake percolates'
        # Remove elements from model and inhom from aquifer
        self.remove_elements()
        # Create new elements
        if isinstance(newlake,Polygonsh):  
            newlakelist = [newlake]
        else: # It is a multi-polygon
            newlakelist = []
            for n in newlake.geoms:
                newlakelist.append(n)
        for n in newlakelist:
            xylake = array(n.exterior); xylake = xylake.tolist()
            xylake = check_direction(xylake)
            self.setup( xylake )
        for a in areasinklist:
            self.create_areasink(a)
        return True      

    def remove_elements(self):
        for e in self.elementList:
            self.ml.elementList.remove(e)
        for a in self.areasinkList:
            self.ml.elementList.remove(a)
        self.elementList = []
        self.areasinkList = []; self.areasinkarea = 0.0
        for e in self.inhomList:
            self.ml.aq.inhomList.remove( e )
        self.inhomList = []
        
    def create_areasink(self,newareasink):
        xyareasink = array(newareasink.exterior); xyareasink.tolist()
        xyareasink = check_direction(xyareasink)
        infil = ( self.hlake - self.ml.aq.zt[0] ) / self.cbot
        self.areasinkList.append( PolyAreaSink(self.ml,xyareasink,infil) )
        self.areasinkarea += newareasink.area


    
