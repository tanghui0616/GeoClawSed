
    """
    Module to create topo and qinit data files for this example.
    """
    from  clawpack.geoclaw import sedtools
    from numpy import *
    def makesed():
        """
        Output topography file for the entire domain
        """
        nxpoints = 2001
        nypoints = 2001
        xlower = 0
        xupper = 200000.e0
        yupper = 200000.e0
        ylower = 0
        outfile= "sed.type2"
        topotools.sed2writer(outfile,sed,xlower,xupper,ylower,yupper,nxpoints,nypoints)
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
        topotools.sed1writer(outfile,qinit,xlower,xupper,ylower,yupper,nxpoints,nypoints)

    def sed(x,y):
        """
        Test case
        """
        # value of z at origin:  Try zmin = 80 for shoreline or 250 for no shore
        from pylab import plot,show,ylim
        nxpoints = 2001
        nypoints = 2001
        slope = 0.0275 #slope
        endpoint = 80000 #endpoint
        th = zeros(shape=shape(nxpoints,nypoints))
        for i in range(nxpoints):
            for j in range(nypoints):
                th(i,j) = 0.5
        return th

    def qinit(x,y):
        """
        Gaussian hump:
        """
        from numpy import where
        from pylab import plot,show
        nxpoints = 2001
        nypoints = 2001
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
        return hu

    if __name__=='__main__':
        makesed()
        makeqinit()
