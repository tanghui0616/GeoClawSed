
    """
    Module to create topo and qinit data files for this example.
    """
    from  clawpack.geoclaw import sedtools
    from numpy import *
    def makesed():
        """
        Output sediment file for the entire domain including sediment thickness and grainsize distribution
        """
        nxpoints = 2001
        nypoints = 2001
        xlower = 0
        xupper = 200000.e0
        yupper = 200000.e0
        ylower = 0
        outfile1= "thick.type2"
        outfile2= "perc.type2"
        sedtools.sed2writer(outfile1,thickness,xlower,xupper,ylower,yupper,nxpoints,nypoints)
        sedtools.sed2writer(outfile2,perc,xlower,xupper,ylower,yupper,nxpoints,nypoints,nclass)
    


    def thickness(x,y):
        """
        Test case for sediment transport
        """
        # value of z at origin:  Try zmin = 80 for shoreline or 250 for no shore
        nxpoints = 2001
        nypoints = 2001
        endpoint = 80000 #endpoint
        th = zeros(shape=shape(nxpoints,nypoints))
        for i in range(nxpoints):
            for j in range(nypoints):
                th(i,j) = 0.5
        return th
        
    def perc(x,y):
        """
        Test case for sediment transport
        """
        # value of z at origin:  Try zmin = 80 for shoreline or 250 for no shore
        nxpoints = 2001
        nypoints = 2001
        nclass = 2
        endpoint = 80000 #endpoint
        perc = zeros(shape=shape(nxpoints,nypoints,nclass))
        for i in range(nxpoints):
            for j in range(nypoints):
                for k in range (nclass):
                    if (i > 1000):
                        perc(i,j,0) = 1.0
                        perc(i,j,1) = 0.0
                    else:
                        perc(i,j,0) = 0.0
                        perc(i,j,1) = 1.0
        return perc


    if __name__=='__main__':
        makesed()

