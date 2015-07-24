
"""    Module to create topo and qinit data files for this example.
"""
from clawpack.geoclaw.sediment.sedpertools import per1writer
from clawpack.geoclaw.sediment.sedthtools import sed1writer
from numpy import *
def makesed():
    """
    Output sediment file for the entire domain including sediment thickness and grainsize distribution
    """
    nxpoints = 11
    nypoints = 11
    nclass = 2
    xlower = 0
    xupper = 200000.e0
    yupper = 200000.e0
    ylower = 0
    outfile1= "thick.type2"
    outfile2= "perc.type2"
    sed1writer(outfile1,thickness,xlower,xupper,ylower,yupper,nxpoints,nypoints)
    per1writer(outfile2,perc,xlower,xupper,ylower,yupper,nxpoints,nypoints,nclass)
    


def thickness(x,y):
    """
    Test case for sediment transport
    """
    # value of z at origin:  Try zmin = 80 for shoreline or 250 for no shore
    nxpoints = 11
    nypoints = 11
    th = zeros((nxpoints,nypoints))
    for i in range(nxpoints):
        for j in range(nypoints):
            th[i,j] = 0.5
    #th(i,j) = 0.5
    return th
        
def perc(x,y):
    """
    Test case for sediment transport
    """
    # value of z at origin:  Try zmin = 80 for shoreline or 250 for no shore
    nxpoints = 11
    nypoints = 11
    nclass = 2
    endpoint = 80000 #endpoint
    percentage = zeros((nxpoints,nypoints,nclass))
    for i in range(nxpoints):
        for j in range(nypoints):
            for k in range (nclass):
                if (i > 5):
                    percentage[i,j,0] = 1.0
                    percentage[i,j,1] = 0.0
                else:
                    percentage[i,j,0] = 0.0
                    percentage[i,j,1] = 1.0
    return percentage


if __name__=='__main__':
    makesed()

