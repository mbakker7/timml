from mllinesink import *
from mllinedoublet import *
from mllinesinkgeneral import *
from mlaquifer import *

def MakeInhomPolySide(ml,xylist,order,closed=False):
    ''' Creates analytic elements along boundary segment of inhomogeneity.
    Will close loop is closed is set to True'''

    elementList = []
    aqleft,aqright = FindAquiferLeftRight(ml,xylist[0],xylist[1])
    if aqleft.Naquifers == 1 and aqright.Naquifers == 1:
        for i in range(len(xylist)-1):
            ld = LineDoublet(ml,xylist[i][0],xylist[i][1],xylist[i+1][0],xylist[i+1][1], aqleft, aqright, order)
            elementList.append(ld)
        if closed:
            ld = LineDoublet(ml,xylist[-1][0],xylist[-1][1],xylist[0][0],xylist[0][1], aqleft, aqright, order)
            elementList.append(ld)
    else:        
        for i in range(len(xylist)-1):
            ls = DoubleLineSinkGen(ml,xylist[i][0],xylist[i][1],xylist[i+1][0],xylist[i+1][1],order, aqleft, aqright )
            ld = LineDoublet(ml,xylist[i][0],xylist[i][1],xylist[i+1][0],xylist[i+1][1], aqleft, aqright, order)
            elementList.extend((ls,ld))
        if closed:
            ls = DoubleLineSinkGen(ml,xylist[-1][0],xylist[-1][1],xylist[0][0],xylist[0][1],order, aqleft, aqright )
            ld = LineDoublet(ml,xylist[-1][0],xylist[-1][1],xylist[0][0],xylist[0][1], aqleft, aqright, order)
            elementList.extend((ls,ld))
    return elementList
            
def FindAquiferLeftRight(ml,xy1,xy2):
    z1 = complex( xy1[0], xy1[1] )
    z2 = complex( xy2[0], xy2[1] )
    Zcp = 1e-6j  # control point just on inside
    zcp = Zcp * (z2 - z1) / 2.0 + 0.5 * (z1 + z2)
    xcpin = zcp.real; ycpin = zcp.imag
    Zcp = -1e-6j  # control point just on outside
    zcp = Zcp * (z2 - z1) / 2.0 + 0.5 * (z1 + z2)
    xcpout = zcp.real; ycpout = zcp.imag
    aqleft = ml.aq.findAquiferData(xcpin,ycpin)
    aqright = ml.aq.findAquiferData(xcpout,ycpout)
    return aqleft,aqright

def MakeInhomSide(ml,xylist,aqleft,aqright,order,closed=False):
    ''' Creates analytic elements along boundary segment of inhomogeneity.
    Will close loop is closed is set to True'''
    print 'MakeInhomSide is deprecated and replaced by MakeInhomPolySide'
    
    if aqleft.Naquifers == 1 and aqright.Naquifers == 1:
        for i in range(len(xylist)-1):
            LineDoublet(ml,xylist[i][0],xylist[i][1],xylist[i+1][0],xylist[i+1][1], aqleft, aqright, order)
        if closed:
            LineDoublet(ml,xylist[-1][0],xylist[-1][1],xylist[0][0],xylist[0][1], aqleft, aqright, order)
    else:        
        for i in range(len(xylist)-1):
            DoubleLineSinkGen(ml,xylist[i][0],xylist[i][1],xylist[i+1][0],xylist[i+1][1],order, aqleft, aqright )
            LineDoublet(ml,xylist[i][0],xylist[i][1],xylist[i+1][0],xylist[i+1][1], aqleft, aqright, order)
        if closed:
            DoubleLineSinkGen(ml,xylist[-1][0],xylist[-1][1],xylist[0][0],xylist[0][1],order, aqleft, aqright )
            LineDoublet(ml,xylist[-1][0],xylist[-1][1],xylist[0][0],xylist[0][1], aqleft, aqright, order)



def MakeCompInhomogeneity(modelParent,xylist,aqin,aqout,Ndegree=2,overspec=1):
    # Create a comprehensive inhomogeneity (all layers jump same, preferably have equal k)
    # Create line-doublets
    ldList = []
    for i in range(len(xylist)-1):
        ld = LineDoublet(modelParent,xylist[i][0],xylist[i][1],xylist[i+1][0],xylist[i+1][1],aqin,aqout,Ndegree,overspec)
        ldList = ldList + [ld]
    ld = LineDoublet(modelParent,xylist[len(xylist)-1][0],xylist[len(xylist)-1][1],xylist[0][0],xylist[0][1],aqin,aqout,Ndegree,overspec)
    ldList = ldList + [ld]
    ldList[0].ldLeft = ldList[len(ldList)-1]
    for i in range(1,len(ldList)):
        ldList[i].ldLeft = ldList[i-1]
    for i in range(len(ldList)-1):
        ldList[i].ldRight = ldList[i+1]
    ldList[len(ldList)-1].ldRight = ldList[0]

def MakeInhomogeneity(modelParent,xylist,aqin,aqout,refineLS=1):
    # Create an inhomogeneity
    # Create line-doublets

    assert (aqin.Naquifers==1 or aqout.Naquifers==1), 'TimML error: either inside or outside aquifer must consist of only one layer'

    LineDoubletString(modelParent,xylist,aqin,aqout)

    # Create line-sinks

    if aqin.Naquifers > 1 or aqout.Naquifers > 1:
        xlist = []; ylist = []
        for xy in xylist:
            xlist = xlist + [ xy[0] ]
            ylist = ylist + [ xy[1] ]
        xlist = array(xlist,'d')
        ylist = array(ylist,'d')
        
        if aqin.Naquifers > 1: aquifer = aqin
        if aqout.Naquifers > 1: aquifer = aqout
        aqlist = range(1,aquifer.Naquifers+1)    
        for i in range(len(xlist)-1):
            LineSinkGen(modelParent,xlist[i],ylist[i],xlist[i+1],ylist[i+1],[0,0,0],2,aqlist,aquifer)
        LineSinkGen(modelParent,xlist[-1],ylist[-1],xlist[0],ylist[0],[0,0,0],2,aqlist,aquifer)

def AquiferSystemInhomogeneity(modelParent,xylist,aqin,aqout):
    '''Create aquifer system inhomogeneity'''

    assert (aqin.Naquifers==1 or aqout.Naquifers==1), 'TimML error: either inside or outside aquifer must consist of only one layer'
        
    xlist = []; ylist = []
    for xy in xylist:
        xlist = xlist + [ xy[0] ]
        ylist = ylist + [ xy[1] ]
    xlist = array(xlist,'d')
    ylist = array(ylist,'d')

    for i in range(len(xlist)-1):
        LineDoublet(modelParent,xlist[i],ylist[i],xlist[i+1],ylist[i+1],aqin,aqout)
    LineDoublet(modelParent,xlist[-1],ylist[-1],xlist[0],ylist[0],aqin,aqout)

    if aqin.Naquifers > 1 or aqout.Naquifers > 1: # create line-sinks
    
        if aqin.Naquifers > 1: aquifer = aqin
        if aqout.Naquifers > 1: aquifer = aqout
        aqlist = range(1,aquifer.Naquifers+1)    
        for i in range(len(xlist)-1):
            LineSinkGen(modelParent,xlist[i],ylist[i],xlist[i+1],ylist[i+1],[0,0,0],2,aqlist,aquifer)
        LineSinkGen(modelParent,xlist[-1],ylist[-1],xlist[0],ylist[0],[0,0,0],2,aqlist,aquifer)

##def MakeInhomogeneityOld(modelParent,xylist,aqin,aqout,refineLS=1):
##    # Create an inhomogeneity
##    # Create line-doublets
##
##    assert (aqin.Naquifers==1 or aqout.Naquifers==1), 'TimML error: either inside or outside aquifer must consist of only one layer'
##    ldList = []
##    for i in range(len(xylist)-1):
##        ld = LineDoublet(modelParent,xylist[i][0],xylist[i][1],xylist[i+1][0],xylist[i+1][1],aqin,aqout)
##        ldList = ldList + [ld]
##    ld = LineDoublet(modelParent,xylist[len(xylist)-1][0],xylist[len(xylist)-1][1],xylist[0][0],xylist[0][1],aqin,aqout)
##    ldList = ldList + [ld]
##    ldList[0].ldLeft = ldList[len(ldList)-1]
##    for i in range(1,len(ldList)):
##        ldList[i].ldLeft = ldList[i-1]
##    for i in range(len(ldList)-1):
##        ldList[i].ldRight = ldList[i+1]
##    ldList[len(ldList)-1].ldRight = ldList[0]
##
##    # Create line-sinks
##
##    if aqin.Naquifers > 1 or aqout.Naquifers > 1:
##        xlist = []; ylist = []
##        for xy in xylist:
##            xlist = xlist + [ xy[0] ]
##            ylist = ylist + [ xy[1] ]
##        xlist = array(xlist,'d')
##        ylist = array(ylist,'d')
##        
##        if refineLS > 1:
##            xlistnew = []; ylistnew = []
##            for k in range(len(xlist)-1):
##                for i in range(int(refineLS)):
##                    xlistnew = xlistnew + [ xlist[k] + i * (xlist[k+1]-xlist[k]) / int(refineLS) ]
##                    ylistnew = ylistnew + [ ylist[k] + i * (ylist[k+1]-ylist[k]) / int(refineLS) ]
##            for i in range(int(refineLS)):
##                    xlistnew = xlistnew + [ xlist[-1] + i * (xlist[0]-xlist[-1]) / int(refineLS) ]
##                    ylistnew = ylistnew + [ ylist[-1] + i * (ylist[0]-ylist[-1]) / int(refineLS) ] 
##            xlist = xlistnew; ylist = ylistnew
##
##        if aqin.Naquifers > 1: aquifer = aqin
##        if aqout.Naquifers > 1: aquifer = aqout
##        aqlist = range(1,aquifer.Naquifers+1)    
##        for i in range(len(xlist)-1):
##            LineSink(modelParent,xlist[i],ylist[i],xlist[i+1],ylist[i+1],0,aqlist,aquifer)
##        LineSink(modelParent,xlist[-1],ylist[-1],xlist[0],ylist[0],0,aqlist,aquifer)

