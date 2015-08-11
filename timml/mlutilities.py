'''
mlutilities.py contain scripts for the creation of Surfer and Matlab
grids and and layouts.
This file is part of the TimML library and is distributed under
the GNU LPGL. See the TimML.py file for more details.
(c) Mark Bakker, 2002-2007
'''
from ml import *
from numpy import *
import pickle

def surfgrid(ml,xmin,xmax,nx,ymin,ymax,ny,filename='/temp/dump',Naquifers=1):
    '''Give filename without extension'''
    xstep = float(xmax-xmin)/nx
    ystep = float(ymax-ymin)/ny
    rows = []   # Store every matrix in one long row
    if Naquifers == 'all':
        aquiferRange = range(ml.aq.Naquifers)
    if type(Naquifers) == list:
        aquiferRange = Naquifers
    else:
        aquiferRange = range(Naquifers)
    for i in range(max(aquiferRange)+1):
        rows = rows + [[]]        
    for j in range(ny+1):
        y=ymin + j*ystep
        for i in range(nx+1):
            x=xmin + i*xstep;
            hvec = ml.headVector(x,y)
            for k in aquiferRange:
                if k < len(hvec):
                    rows[k] = rows[k] + [hvec[k]]
                else:
                    rows[k] = rows[k] + [hvec[0]]
    for k in aquiferRange:
        zmin = min(rows[k]); zmax = max(rows[k])
        out = open(filename+'.'+str(k+1)+'.grd','w')
        out.write('DSAA\n')
        out.write(str(nx+1)+' '+str(ny+1)+'\n')
        out.write(str(xmin)+' '+str(xmin+nx*xstep)+'\n')
        out.write(str(ymin)+' '+str(ymin+ny*ystep)+'\n')
        out.write(str(zmin)+' '+str(zmax)+'\n')
        for j in range(ny+1):
            row = rows[k][j*(nx+1) : (j+1)*(nx+1)]
            for i in range(nx):
                out.write(str(row[i])+' ')
            out.write(str(row[nx])+'\n')
        out.close

def surfergrid(ml,xmin,xstep,xmax,ymin,ystep,ymax,filename='/temp/dump',Naquifers=1):
    '''Give filename without extension'''
    nx=int(round((xmax-xmin)/xstep))  # nx is number of intervals
    ny=int(round((ymax-ymin)/ystep))  # ny is number of intervals
    surfgrid(ml,xmin,xmax,nx,ymin,ymax,ny,filename,Naquifers)
        
def matgrid(ml,xmin,xmax,nx,ymin,ymax,ny,filename='/temp/dump.m',Naquifers=1):
    '''Spits out matlab m file; x and y are stored in xg,yg; heads in h1, h2, ...'''
    xstep = float(xmax-xmin)/nx
    ystep = float(ymax-ymin)/ny
    rows = []   # Store every matrix in one long row
    if Naquifers == 'all':
        aquiferRange = range(ml.aq.Naquifers)
    if type(Naquifers) == list:
        aquiferRange = Naquifers
    else:
        aquiferRange = range(Naquifers)
    out = open(filename,'w')
#    out.write('[xg,yg]=meshgrid('+str(xmin)+'+'+str(xstep)+'*[0:'+str(nx)+']'\
#              ','+str(ymin)+'+'+str(ystep)+'*[0:'+str(ny)+']);\n')
    out.write('[xg,yg]=meshgrid(linspace('+str(xmin)+','+str(xmax)+','+str(nx+1)+\
              '),linspace('+str(ymin)+','+str(ymax)+','+str(ny+1)+'));\n')
    print 'grid of '+str((nx,ny))+'. gridding in progress. hit ctrl-c to abort'
    rows = zeros((nx+1,len(aquiferRange)),'d')
    for j in range(ny+1):
        y=ymin + j*ystep
        for i in range(nx+1):
            x=xmin + i*xstep;
            h = ml.headVector(x,y)
            for k in range(len(aquiferRange)):
                if aquiferRange[k] >= len(h):
                    rows[i,k] = h[0]
                else:
                    rows[i,k] = h[aquiferRange[k]]
        for k in range(len(aquiferRange)):
            # Need to convert rows to list, else the str function inserts \n-s
            out.write( 'h' + str(aquiferRange[k]+1) + '(' + str(j+1) + ',:)=' + str(list(rows[:,k])) + ';\n' )
    out.close

def matlabgrid(ml,xmin,xstep,xmax,ymin,ystep,ymax,filename='/temp/dump.m',Naquifers=1):
    '''Spits out matlab m file; x and y are stored in xg,yg; heads in h1, h2, ...'''
    nx=int(round((xmax-xmin)/xstep))  # nx is number of intervals
    ny=int(round((ymax-ymin)/ystep))  # ny is number of intervals
    matgrid(ml,xmin,xmax,nx,ymin,ymax,ny,filename,Naquifers)

def matlabqxgrid(ml,xmin,xstep,xmax,ymin,ystep,ymax,filename='/temp/dump.m',Naquifers=1):
    '''Spits out matlab m file; x and y are stored in xg,yg; heads in h1, h2, ...'''
    nx=int((xmax-xmin)/xstep)  # nx is number of intervals
    ny=int((ymax-ymin)/ystep)  # ny is number of intervals
    rows = []   # Store every matrix in one long row
    for i in range(Naquifers):
        rows = rows + [[]]        
    for j in range(ny+1):
        y=ymin + j*ystep
        for i in range(nx+1):
            x=xmin + i*xstep; 
            [qxvec,qyvec] = ml.dischargeVector(x,y)
            for k in range(Naquifers):
                if k < len(qxvec):
                    rows[k] = rows[k] + [qxvec[k]]
                else:
                    rows[k] = rows[k] + [qxvec[0]]
    out = open(filename,'w')
    out.write('[xg,yg]=meshgrid('+str(xmin)+':'+str(xstep)+':'+str(xmax)+\
              ','+str(ymin)+':'+str(ystep)+':'+str(ymax)+');\n')
    for k in range(Naquifers):
        out.write('h'+str(k+1)+'=[\n')
        for j in range(ny+1):
            row = rows[k][j*(nx+1) : (j+1)*(nx+1)]
            for i in range(nx):
                out.write(str(row[i])+' ')
            out.write(str(row[nx])+';\n')
        out.write('];\n')
        out.close

def matlabqygrid(ml,xmin,xstep,xmax,ymin,ystep,ymax,filename='/temp/dump.m',Naquifers=1):
    '''Spits out matlab m file; x and y are stored in xg,yg; heads in h1, h2, ...'''
    nx=int((xmax-xmin)/xstep)  # nx is number of intervals
    ny=int((ymax-ymin)/ystep)  # ny is number of intervals
    rows = []   # Store every matrix in one long row
    for i in range(Naquifers):
        rows = rows + [[]]        
    for j in range(ny+1):
        y=ymin + j*ystep
        for i in range(nx+1):
            x=xmin + i*xstep; 
            [qxvec,qyvec] = ml.dischargeVector(x,y)
            for k in range(Naquifers):
                if k < len(qyvec):
                    rows[k] = rows[k] + [qyvec[k]]
                else:
                    rows[k] = rows[k] + [qyvec[0]]
    out = open(filename,'w')
    out.write('[xg,yg]=meshgrid('+str(xmin)+':'+str(xstep)+':'+str(xmax)+\
              ','+str(ymin)+':'+str(ystep)+':'+str(ymax)+');\n')
    for k in range(Naquifers):
        out.write('h'+str(k+1)+'=[\n')
        for j in range(ny+1):
            row = rows[k][j*(nx+1) : (j+1)*(nx+1)]
            for i in range(nx):
                out.write(str(row[i])+' ')
            out.write(str(row[nx])+';\n')
        out.write('];\n')
        out.close

def matlabvertgrid(ml,xmin,xstep,xmax,zmin,zstep,zmax,ycrosssection,filename='/temp/dump.m'):
    '''Spits out matlab m file; x and y are stored in xg,yg; heads in h1, h2, ...'''
    nx=int((xmax-xmin)/xstep)  # nx is number of intervals
    nz=int((zmax-zmin)/zstep)  # nz is number of intervals
    rows = []   # Store matrix in one long row     
    for j in range(nz+1):
        z=zmin + j*zstep
        for i in range(nx+1):
            x=xmin + i*xstep; 
            head = ml.head3D(x,ycrosssection,z)
            rows = rows + [head]
    out = open(filename,'w')
    out.write('[xg,yg]=meshgrid('+str(xmin)+':'+str(xstep)+':'+str(xmax)+\
              ','+str(zmin)+':'+str(zstep)+':'+str(zmax)+');\n')
    out.write('h=[\n')
    for j in range(nz+1):
        row = rows[j*(nx+1) : (j+1)*(nx+1)]
        for i in range(nx):
            out.write(str(row[i])+' ')
        out.write(str(row[nx])+';\n')
    out.write('];\n')
    out.close

def surfvertgrid(ml,xmin,ymin,xmax,ymax,nh,zmin,zmax,nz,filename='/temp/dump',interp=1):
    '''Give filename without extension'''
    L = sqrt( (xmax-xmin)**2 + (ymax-ymin)**2 )
    xstep = float(xmax-xmin) / nh
    ystep = float(ymax-ymin) / nh
    zstep = float(zmax-zmin) / nz
    grid = zeros((nz+1,nh+1),'d')
    for i in range(nz+1):
        z = zmin + i*zstep
        for j in range(nh+1):
            x = xmin + j*xstep; y = ymin + j*ystep
            if interp:
                grid[i,j] = ml.head3Dinterp(x,y,z)
            else:
                grid[i,j] = ml.head3D(x,y,z)
    hmin = min(min(grid)); hmax = max(max(grid))
    out = open(filename+'.grd','w')
    out.write('DSAA\n')
    out.write(str(nh+1)+' '+str(nz+1)+'\n')
    out.write('0 '+str(L)+'\n')
    out.write(str(zmin)+' '+str(zmin+nz*zstep)+'\n')
    out.write(str(hmin)+' '+str(hmax)+'\n')
    for i in range(nz+1):
        for j in range(nh):
            out.write(str(grid[i,j])+' ')
        out.write(str(grid[i,-1])+'\n')
    out.close

def matlablayout(ml,filename='/temp/dump.m',color='k'):
    '''filename between quotes with .m extension. Layout of all elements, independent of layer'''
    out = open(filename,'w')
    for e in ml.elementList:
        a = e.layout()
        nterms = len(a)
        for i in range(0,nterms,3):
            if a[i] > 0:
                out.write('plot('+str(a[i+1])+','+str(a[i+2])+",'"+color+"')\n")
    out.close

def surferlayout(ml,filename='/temp/dump.bln'):
    '''filename between quotes with .bln extension. Layout of all elements, independent of layer'''
    out = open(filename,'w')
    for e in ml.elementList:
        a = e.layout()
        nterms = len(a)
        for k in range(0,nterms,3):
            if a[k] > 0:
                out.write(str(a[k])+' 1\n')
                for i in range(a[k]):
                    out.write(str(a[k+1][i])+' '+str(a[k+2][i])+'\n')
    out.close

def saveModel(ml,filename):
    file = open(filename,'w')
    pickle.dump(ml,file)
    file.close()

def loadModel(filename):
    file = open(filename,'r')
    ml = pickle.load(file)
    file.close()
    return ml

def NumLap(f,x,y,xstep):
    f0 = f(x,y)
    f1 = f(x+xstep,y)
    f2 = f(x,y+xstep)
    f3 = f(x-xstep,y)
    f4 = f(x,y-xstep)
    numlap = (f1 + f2 + f3 + f4 - 4.0 * f0) / (xstep**2)
    return numlap

def integratedLeftDischarge(ml,xy,filename):
    Q = ml.integratedLeftDischarge(xy)
    savetxt(filename,Q)

