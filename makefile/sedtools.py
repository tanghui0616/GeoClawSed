#!/usr/bin/env python
# encoding: utf-8

    """GeoClaw Sediments Tools Module
    
    Module provides several functions for reading, writing and manipulating
    sediments (deposit) files.
    
    :Functions:
    - dms2decimal - Convert (degrees, minutes, seconds) to decimal degrees
    - dist_meters2latlong - Convert dx, dy distance in meters to degrees
    - dist_latlong2meters - Convert dx, dy distance in degrees to meters
    - haversine - Calculate the haversine based great circle distance
    - inv_haversine - Inverts the haversine distance
    
    :Topography Class:
    
    :TODO:
    - Tests are implemented but not passing, should we expect the arrays to be
    identical?
    - Add sub and super sampling capababilities
    - Add functions for creating sediment based off a sediment function, incorporate
    the create_sed_func into sediment class, maybe allow more broad
    initialization ability to the class to handle this?
    - Fix `in_poly` function
    - Add remove/fill no data value
    - Probably should better handle remote files (fetching from http)
    - Add more robust plotting capabilities
    """
    import os
    import urllib
    import types

    import numpy

    def topo1writer (outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints):
        """
        Function topo1writer will write out the topofiles by evaluating the
        function topo on the grid specified by the other parameters.
        
        Assumes topo can be called on arrays X,Y produced by numpy.meshgrid.
        
        Output file is of "topotype1," which we use to refer to a file with
        (x,y,z) values on each line, progressing from upper left corner across
        rows, then down.
        """
        topography = Topography(topo_func=topo)
    
        topography.x = numpy.linspace(xlower,xupper,nxpoints)
        topography.y = numpy.linspace(ylower,yupper,nypoints)
    
        topography.write(outfile, topo_type=1)


    def topo2writer (outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints, \
                 nodata_value=-99999):
        """Write out a topo type 2 file by evaluating the function *topo*.
        
        This routine is here for backwards compatibility and simply creates a new
        topography object and writes it out.
        
        """
    
        topography = Topography(topo_func=topo)
    
        topography.x = numpy.linspace(xlower,xupper,nxpoints)
        topography.y = numpy.linspace(ylower,yupper,nypoints)
    
        topography.write(outfile, topo_type=2)


    def topo3writer (outfile,topo,xlower,xupper,ylower,yupper,nxpoints,nypoints, \
                 nodata_value=-99999):
        """Write out a topo type 3 file by evaluating the function *topo*.
        
        This routine is here for backwards compatibility and simply creates a new
        topography object and writes it out.
        
        """
    
        topography = Topography(topo_func=topo)
    
        topography.x = numpy.linspace(xlower,xupper,nxpoints)
        topography.y = numpy.linspace(ylower,yupper,nypoints)
    
        topography.write(outfile, topo_type=3)

