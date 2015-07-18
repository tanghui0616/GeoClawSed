
"""
Module to create topo and qinit data files for this example.
"""

from clawpack.geoclaw import topotools
from numpy import *

def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints = 2001
    nypoints = 2001
    xlower = 0
    xupper = 200000.e0
    yupper = 200000.e0
    ylower = 0
    outfile= "bowl.topotype2"     
    topotools.topo2writer(outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def makeqinit():
    """
    Create qinit data file
    """
    nxpoints = 2001
    nypoints = 2001
    xlower = 0.e0
    xupper = 200000.e0
    yupper = 200000.e0
    ylower = 0.e0
    outfile= "hump.xyz"     
    topotools.topo1writer(outfile,qinit,xlower,xupper,ylower,yupper,nxpoints,nypoints)

def topo(x,y):
    """
    Parabolic bowl
    """
    # value of z at origin:  Try zmin = 80 for shoreline or 250 for no shore
    from pylab import plot,show,ylim
    zmin = 500.0 #
    nxpoints = 2001
    nypoints = 2001
    slope = 0.0075 #slope
    endpoint = 80000 #endpoint
    z = zeros(shape=shape(y))
    for i in range(nxpoints):
        for j in range(nypoints):
            if (y[i,j]<endpoint):
                z[i,j] = slope*abs(y[i,j]-endpoint)-zmin
            else:
                z[i,j] = -zmin
    plot(y,z,'-')
    ylim(-600,200)
    show()
#z = 1.e-2* (x-x) - zmin
    return z


def qinit(x,y):
    """
    Gaussian hump:
    """
    from numpy import where
    from pylab import plot,show
    nxpoints = 2001
    nypoints = 2001
    #import sympy
        #def sech(x):
        #return sympy.cosh(x)**(-1)
    zmin = 1000.0
    zmin1 = 1000
    g = 9.81
    A = 10.0 # wave height
    k = sqrt(3*A/(4*zmin**3))
    z = zeros(shape=shape(y))
    u = zeros(shape=shape(y))
    hu = zeros(shape=shape(y))
    y0 = 180000
    c = sqrt(g*(A+zmin))
    rho = 1.0e3
    for i in range(nxpoints):
        for j in range(nypoints):
            z[i,j] = A*cosh(k*(y[i,j]-y0))**(-2)
            u[i,j] = -sqrt(g/zmin)*z[i,j]
            hu[i,j] =(z[i,j]+zmin)*u[i,j]
    plot(y,u,'-')
    show()
    #ze = -((y+0e0)**2)/10.
    #z = where(ze>-400000., 400.e0*exp(ze), 0.)
    return hu

if __name__=='__main__':
    maketopo()
    makeqinit()
