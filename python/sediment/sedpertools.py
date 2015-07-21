#!/usr/bin/env python
# encoding: utf-8
    
    """GeoClaw Sediments Grainsize distribution Tools Module
        
        Module provides several functions for reading, writing and ploting
        sediments grainsize distribution files.
        
        :Functions:
        - per1writer
        - per2writer

        
        
        
        :Percentage Class:
        
        
    """
import os
import urllib
import types
import numpy

def per1writer (outfile,perc,xlower,xupper,ylower,yupper,nxpoints,nypoints,nclasses):
    """
        Function sed1writer will write out the sediment thickness profile by evaluating the
        function thickness on the grid specified by the other parameters.
        
        Assumes thickness can be called on arrays X,Y produced by numpy.meshgrid.
        
        Output file is of "sedtype1," which we use to refer to a file with
        (x,y,z) values on each line, progressing from upper left corner across
        rows, then down.
        """
    percentage = Percentage(Per_func=perc)
    percentage.x = numpy.linspace(xlower,xupper,nxpoints)
    percentage.y = numpy.linspace(ylower,yupper,nypoints)
    percentage.write(outfile, Per_type=1)

def per2writer (outfile,perc,xlower,xupper,ylower,yupper,nxpoints,nypoints,nclasses):
    """
        Function sed1writer will write out the sediment thickness profile by evaluating the
        function thickness on the grid specified by the other parameters.
        
        Assumes thickness can be called on arrays X,Y produced by numpy.meshgrid.
        
        Output file is of "sedtype1," which we use to refer to a file with
        (x,y,z) values on each line, progressing from upper left corner across
        rows, then down.
        """
    percentage = Percentage(Per_func=perc)
    percentage.x = numpy.linspace(xlower,xupper,nxpoints)
    percentage.y = numpy.linspace(ylower,yupper,nypoints)
    percentage.write(outfile, Per_type=2)

# ==============================================================================
#  Percentage class
# ==============================================================================
class Percentage(object):

    """Base sediment class.
    
        A class representing a single sediment thickness file
    
        :Properties:
    
        :Initialization:
        -
    
        :Examples:
    
        >>> import clawpack.geoclaw.tools as sed
        >>> sed_file = sed.Sediment('./sed.tt3')
        >>> sed_file.plot()
    
    """
    @property
    def per(self):
        r"""A representation of the data as an 1d array."""
        if (self._per is None) and self.unstructured:
            self.read(mask=True)
        return self._per
    @per.setter
    def per(self, value):
        self._per = value
    @per.deleter
    def per(self):
        del self._per
    
    @property
    def Per(self):
        r"""A representation of the data as a 2d array."""
        if self._Per is None:
            self.generate_2d_percentage(mask=True)
        return self._Per
    @Per.setter
    def Per(self, value):
        self._Per = value
    @Per.deleter
    def Per(self):
        del self._Per
    
    @property
    def x(self):
        """One dimensional coorindate array in x direction."""
        if self._x is None:
            self.read(mask=True)
        return self._x
    @x.setter
    def x(self, value):
        self._extent = None
        self._x = value
    @x.deleter
    def x(self):
        del self._x
    
    @property
    def X(self):
        """Two dimensional coordinate array in x direction."""
        if self._X is None:
            self.generate_2d_coordinates(mask=True)
        return self._X
    @X.deleter
    def X(self):
        del self._X
    
    @property
    def y(self):
        r"""One dimensional coordinate array in y direction."""
        if self._y is None:
            self.read(mask=True)
        return self._y
    @y.setter
    def y(self, value):
        self._extent = None
        self._y = value
    @y.deleter
    def y(self):
        del self._y
    
    @property
    def Y(self):
        r"""Two dimensional coordinate array in y direction."""
        if self._Y is None:
            self.generate_2d_coordinates(mask=True)
        return self._Y
    @Y.deleter
    def Y(self):
        del self._Y
    
    @property
    def extent(self):
        r"""Extent of the Sediment."""
        if self._extent is None:
            self._extent = ( numpy.min(self.x), numpy.max(self.x),
                            numpy.min(self.y), numpy.max(self.y) )
        return self._extent
    @extent.setter
    def extent(self, value):
        self._extent = value
    
    @property
    def delta(self):
        r"""Spacing of data points."""
        if self._delta is None:
            if self.unstructured:
                
                # Calculate the smallest spacing between grid points
                dx = numpy.infty
                dy = numpy.infty
                num_comparisons = self.x.shape[0] - 1
                for i in xrange(self.x.shape[0]):
                    for j in xrange(num_comparisons):
                        dx = min(dx, numpy.abs(self.x[i + j + 1] - self.x[i]))
                        dy = min(dy, numpy.abs(self.y[i + j + 1] - self.y[i]))
                    
                    num_comparisons -= 1
                self._delta = [dx, dy]
            else:
                # All other Sediment types should have equally spaced grid
                # points in each direction
                begin_delta = numpy.array([abs(self.x[1] - self.x[0]),
                                           abs(self.y[1] - self.y[0])])
                end_delta =   numpy.array([abs(self.x[-2] - self.x[-1]),
                                            abs(self.y[-2] - self.y[-1])])
                assert numpy.allclose(begin_delta, end_delta, 1e-8),   \
                            "Grid spacing delta not constant, %s != %s." %  \
                            (begin_delta, end_delta)
                self._delta = numpy.round(begin_delta[0], 15)
        return self._delta

    def __init__(self, path=None, per_func=None, per_type=None,unstructured=False):
        """
            
            Sediment Grainsize Distribution initialization routine.
        
        """
    
        super(Sediment, self).__init__()
    
        self.path = path
        self.per_func = per_func
        self.per_type = per_type
        self.unstructured = unstructured
        self.no_data_value = -9999
    
        # Data storage for only calculating array shapes when needed
        self._per = None
        self._Per = None
        self._x = None
        self._X = None
        self._y = None
        self._Y = None
        self._extent = None
        self._delta = None


    def generate_2d_percentage(self, mask=True):
        """Generate a 2d array of the sediment thickness."""
    
        # Check to see if we need to generate these
        if self._Per is None:
        
            if self.unstructured:
            
                raise ValueError("Unstructured data does not allow for use of" \
                                 + " 2d arrays")
        
        if self.path is not None:
            if self._per is None:
                # Try to read the data, may not have done this yet
                self.read(path=self.path, mask=mask)
                if self._Per is not None:
                    # We are done, the read function did our work
                    return
                
                # See if self._X and self._Y are already computed and use them if
                # available, otherwise just use self._x and self._y
                if self._X is not None and self._Y is not None:
                    new_shape = self._X.shape
                else:
                    new_shape = (self._x.shape[0], self._y.shape[0])
                # Reshape, note that the mask follows along with the new array
                self._Per = numpy.reshape(self._per, new_shape)
        
        elif self.per_func is not None:
            # Generate sediment thickness profile via sed_func
            self._Per = self.per_func(self.X, self.Y)

    def generate_2d_coordinates(self, mask=True):
        """Generate 2d coordinate arrays."""
    
        # Check to see if we need to generate these
        if self._X is None and self._Y is None:
        
            if (self._x is not None) and (self._y is not None):
                self._X,self._Y = numpy.meshgrid(self._x, self._y)
    
        if self._X is None and self._Y is None:
            if self.unstructured:
                # Really no way to do this here with performing interpolation via
                # extract.  Note that if the interpolation is performed these
                # arrays are already stored in self._X and self._Y
                raise ValueError("Unstructured data does not allow for use of" \
                            + " 2d coordinates")
        
            if self.path is not None:
                if abs(self.per_type) == 1:
                # Reading this sed_type should produce the X and Y arrays
                    self.read(mask=mask)
                elif abs(self.per_type) in [2,3]:
                    if self._x is None or self._y is None:
                    # Try to read the data to get these, may not have been done yet
                        self.read(mask=mask)
                # Generate arrays
                        self._X, self._Y = numpy.meshgrid(self._x, self._y)
                    else:
                        raise ValueError("Unrecognized sed_type: %s" % self.sed_type)
        
            elif self.sed_func is not None:
                if self._x is None or self._y is None:
                    raise ValueError("The x and y arrays must be set to ",
                                 "create 2d coordinate arrays.")
                self._X, self._Y = numpy.meshgrid(self._x, self._y)
        
        
            # If masking has been requested try to get the mask first from
            # self._Sed and then self._sed
            if mask:
                if self._Per is None:
                # Check to see if we really need to do anything here
                if isinstance(self._per, numpy.ma.MaskedArray):
                    # Try to create self._Z
                    self.generate_2d_percentage(mask=mask)
            
            if isinstance(self._Per, numpy.ma.MaskedArray):
                # Use Th's mask for the X and Y coordinates
                self._x = numpy.ma.MaskedArray(self._X, mask=self._Per.mask,
                                               copy=False)
                self._Y = numpy.ma.MaskedArray(self._Y, mask=self._Per.mask,
                                               copy=False)




    def read(self, path=None, sed_type=None, unstructured=False,
             mask=True, filter_region=None, force=False):
            """Read in the data from the object's *path* attribute.
        
            Stores the resulting data in one of the sets of *x*, *y*, and *th* or
            *X*, *Y*, and *Th*.
        
            :Input:
            - *path* (str)  file to read
            - *Per_type* (int)
            - *unstructured* (bool)
            - *mask* (bool)
            - *filter_region* (tuple)
            The first three might have already been set when instatiating object.
        
            """
    
        if (path is None) and (self.path is None):
                raise ValueError("*** Need to set path for file to read")
    
        if path:
            self.path = path   # set or perhaps reset
            self.Per_type = None  # force resetting below
    
        if unstructured:
                self.unstructured = unstructured
    
        if self.Per_type is None:
            if Per_type is not None:
                self.per_type = per_type
        else:
            # Try to look at suffix for type
            extension = os.path.splitext(self.path)[1][1:]
            if extension[:2] == "tt":
                self.per_type = int(extension[2])
            elif extension == 'per':
                self.per_type = 1
        else:
            # Default to 1
            self.sed_type = 1
    
        if self.unstructured:
            # Read in the data as series of tuples
            data = numpy.loadtxt(self.path)
            points = []
            values = []
        
            # Filter region if requested
            if filter_region is not None:
                for coordinate in data:
                    if filter_region[0] <= coordinate[0] <= filter_region[1]:
                        if filter_region[2] <= coordinate[1] <= filter_region[3]:
                            points.append(coordinate[0:2])
                            values.append(coordinate[2])
            
                if len(points) == 0:
                    raise Exception("No points were found inside requested " \
                                    + "filter region.")
            
                # Cast lists as ndarrays
                self._x = numpy.array(points[:,0])
                self._y = numpy.array(points[:,1])
                self._per = numpy.array(values)
        
            else:
                self._x = data[:,0]
                self._y = data[:,1]
                self._per = data[:,2]
    
        else:
            # Data is in one of the GeoClaw supported formats
            if abs(self.per_type) == 1:
                data = numpy.loadtxt(self.path)
                N = [0,0]
                y0 = data[0,1]
                for (n, y) in enumerate(data[1:,1]):
                    if y != y0:
                        N[1] = n + 1
                        break
                N[0] = data.shape[0] / N[1]
            
            self._x = data[:N[1],0]
            self._y = data[::N[1],1]
            self._Th = numpy.flipud(data[:,2].reshape(N))
            self._delta = self.X[0,1] - self.X[0,0]
        
        elif abs(self.sed_type) in [2,3]:
            # Get header information
            N = self.read_header()
            self._x = numpy.linspace(self.extent[0], self.extent[1], N[0])
            self._y = numpy.linspace(self.extent[2], self.extent[3], N[1])
            
            if abs(self.sed_type) == 2:
                # Data is read in as a single column, reshape it
                self._Th = numpy.loadtxt(self.path, skiprows=6).reshape(N[1],N[0])
                self._Th = numpy.flipud(self._Th)
            elif abs(self.sed_type) == 3:
                # Data is read in starting at the top right corner
                self._Th = numpy.flipud(numpy.loadtxt(self.path, skiprows=6))
            
            if mask:
                self._Th = numpy.ma.masked_values(self._Th, self.no_data_value, copy=False)
        
        else:
            raise IOError("Unrecognized sed_type: %s" % self.sed_type)
        
        # Perform region filtering
        if filter_region is not None:
            # Find indices of region
            region_index = [None, None, None, None]
            region_index[0] = (self.x >= filter_region[0]).nonzero()[0][0]
            region_index[1] = (self.x <= filter_region[1]).nonzero()[0][-1]
            region_index[2] = (self.y >= filter_region[2]).nonzero()[0][0]
            region_index[3] = (self.y <= filter_region[3]).nonzero()[0][-1]
            
            self._x = self._x[region_index[0]:region_index[1]]
            self._y = self._y[region_index[2]:region_index[3]]
            
            # Force regeneration of 2d coordinate arrays and extent
            if self._X is not None or self._Y is not None:
                del self._X, self._Y
                self._X = None
                self._Y = None
            self._extent = None
            
            # Modify Z array as well
            self._Th = self._Th[region_index[2]:region_index[3],
                                region_index[0]:region_index[1]]







