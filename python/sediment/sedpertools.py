#!/usr/bin/env python
# encoding: utf-8
    
    """GeoClaw Sediments Grainsize distribution Tools Module
        
        Module provides several functions for reading, writing and ploting
        sediments grainsize distribution files.
        
        :Functions:
        - per1writer
        - per2writer
        - per3writer

        
        
        
        :Percentage Class:
        
        
    """
import os
import urllib
import types
import numpy

# ==============================================================================
#  Functions
# ==============================================================================
def per1writer (outfile,perc,xlower,xupper,ylower,yupper,nxpoints,nypoints,nclasses):
    """
        Function per1writer will write out the sediment thickness profile by evaluating the
        function thickness on the grid specified by the other parameters.
        
        Assumes thickness can be called on arrays X,Y produced by numpy.meshgrid.
        
        Output file is of "pertype1," which we use to refer to a file with
        (x,y,z) values on each line, progressing from upper left corner across
        rows, then down.
        """
    percentage = Percentage(Per_func=perc)
    percentage.x = numpy.linspace(xlower,xupper,nxpoints)
    percentage.y = numpy.linspace(ylower,yupper,nypoints)
    percentage.NC = nclasses
    percentage.write(outfile, Per_type=1)

def per2writer (outfile,perc,xlower,xupper,ylower,yupper,nxpoints,nypoints,nclasses):

    percentage = Percentage(Per_func=perc)
    percentage.x = numpy.linspace(xlower,xupper,nxpoints)
    percentage.y = numpy.linspace(ylower,yupper,nypoints)
    percentage.NC = nclasses
    percentage.write(outfile, Per_type=2)

def per3writer (outfile,perc,xlower,xupper,ylower,yupper,nxpoints,nypoints,nclasses):
    
    percentage = Percentage(Per_func=perc)
    percentage.x = numpy.linspace(xlower,xupper,nxpoints)
    percentage.y = numpy.linspace(ylower,yupper,nypoints)
    percentage.NC = nclasses
    percentage.write(outfile, Per_type=3)

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
        >>> per_file = per.Sediment('./sed.tt3')
        >>> per_file.plot()
    
    """
    # PROPERTIES
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
    def NC(self):
        """Classes of sediment"""
        if self._NC is None:
            self.NC = read(mask=True)
        return self._NC
    @NC.setter
    def NC(self,value):
        self._NC = value
    @NC.deleter
    def NC(self)
        del self._NC

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
                
    # INITIALIZE

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
        self._NC = None
        self._extent = None
        self._delta = None

       self.coordinate_transform = lambda x,y: (x,y)

    def generate_2d_percentage(self, mask=True):
        """Generate a 2d array of the sediment grainsize distribution."""
    
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
                    new_shape = (self._X.shape[0],self._Y.shape[0],self._NC)
                else:
                    new_shape = (self._x.shape[0], self._y.shape[0],self._NC)
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
                # Reading this per_type should produce the X and Y arrays
                    self.read(mask=mask)
                elif abs(self.per_type) in [2,3]:
                    if self._x is None or self._y is None:
                        # Try to read the data to get these, may not have been done yet
                        self.read(mask=mask)
                        # Generate arrays
                        self._X, self._Y = numpy.meshgrid(self._x, self._y)
                else:
                    raise ValueError("Unrecognized per_type: %s" % self.per_type)
        
            elif self.per_func is not None:
                if self._x is None or self._y is None:
                    raise ValueError("The x and y arrays must be set to ",
                                 "create 2d coordinate arrays.")
                self._X, self._Y = numpy.meshgrid(self._x, self._y)
        
        
            # If masking has been requested try to get the mask first from
            # self._per and then self._per
            if mask:
                if self._Per is None:
                    # Check to see if we really need to do anything here
                    if isinstance(self._per, numpy.ma.MaskedArray):
                        # Try to create self._Per
                        self.generate_2d_percentage(mask=mask)
            
                if isinstance(self._Per, numpy.ma.MaskedArray):
                        # Use Th's mask for the X and Y coordinates
                    self._X = numpy.ma.MaskedArray(self._X, mask=self._Per.mask,
                                               copy=False)
                    self._Y = numpy.ma.MaskedArray(self._Y, mask=self._Per.mask,
                                               copy=False)




    def read(self, path=None, per_type=None, unstructured=False,
             mask=True, filter_region=None, force=False):
            """Read in the data from the object's *path* attribute.
        
            Stores the resulting data in one of the sets of *x*, *y*, and *per* or
            *X*, *Y*, and *Per*.
        
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
    
        if self.per_type is None:
            if per_type is not None:
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
                    self.per_type = 3
    
        if self.unstructured:
            # Read in the data as series of tuples
            data = numpy.loadtxt(self.path)
            points = []
            values = []
            self.NC = NC
        

            self._x = data[:,0]
            self._y = data[:,1]
            self._per = data[:,2:self.NC+1]
    
        else:
            # Data is in one of the GeoClaw supported formats
            if abs(self.per_type) == 1:
                data = numpy.loadtxt(self.path)
                N = [0,0,0]
                y0 = data[0,1]
                for (n, y) in enumerate(data[1:,1]):
                    if y != y0:
                        N[1] = n + 1
                        break
                N[2] = self.NC
                N[0] = data.shape[0] / N[1]
            
                self._x = data[:N[1],0]
                self._y = data[::N[1],1]
                self._Per = numpy.flipud(data[:,2:N[3]+1].reshape(N))
                self._delta = self.X[0,1] - self.X[0,0]
        
            elif abs(self.per_type) in [2,3]:
                # Get header information
                N = self.read_header()
                self._x = numpy.linspace(self.extent[0], self.extent[1], N[0])
                self._y = numpy.linspace(self.extent[2], self.extent[3], N[1])
            
                if abs(self.per_type) == 2:
                    # Data is read in as a single column, reshape it
                    self._Per = numpy.loadtxt(self.path, skiprows=6).reshape(N[1],N[0],N[3])
                    self._Per = numpy.flipud(self._Per)
                elif abs(self.per_type) == 3:
                    # Data is read in starting at the top right corner
                    self._Per = numpy.flipud(numpy.loadtxt(self.path, skiprows=6))
            
                if mask:
                    self._Per = numpy.ma.masked_values(self._Per, self.no_data_value, copy=False)
        
            else:
                raise IOError("Unrecognized per_type: %s" % self.per_type)


    def read_header(self):
    r"""Read in header of sediment file at path.
        
        If a value returns numpy.nan then the value was not retrievable.  Note
        that this routine can read in headers whose values and labels are
        swapped.
        
        """
    
        if abs(self.per_type) in [2,3]:
        
            # Default values to track errors
            num_cells = [numpy.nan,numpy.nan,numpy.nan]
            self._extent = [numpy.nan,numpy.nan,numpy.nan,numpy.nan]
            self._delta = numpy.nan
        
            with open(self.path, 'r') as per_file:
                # Check to see if we need to flip the header values
                first_line = per_file.readline()
            try:
                num_cells[0] = int(first_line.split()[0])
            except ValueError:
                # Assume the header is flipped from what we expect
                num_cells[0] = int(first_line.split()[-1])
                value_index = -1
            else:
                value_index = 0
            
            num_cells[1] = int(per_file.readline().split()[value_index])
            num_cells[2] = int(per_file.readline().split()[value_index+1])
            self._extent[0] = float(per_file.readline().split()[value_index])
            self._extent[2] = float(per_file.readline().split()[value_index])
            self._delta = float(per_file.readline().split()[value_index])
            self.no_data_value = float(per_file.readline().split()[value_index])
            
            self._extent[1] = self._extent[0] + num_cells[0] * self.delta
            self._extent[3] = self._extent[2] + num_cells[1] * self.delta
    
        else:
                raise IOError("Cannot read header for per_type %s" % self.per_type)
    
        return num_cells

    def write(self, path, no_data_value=None, per_type=None, masked=True):
        """Write out a Sediment file to path of type *per_type*.
        
            Writes out a Sediment file of per type specified with *per_type* or
            inferred from the output file's extension, defaulting to 3, to path
            from data in Th.  The rest of the arguments are used to write the header
            data.
        
        """
    
        # Determine per type if not specified
        if per_type is None:
            if self.per_type is not None:
                per_type = self.per_type
            else:
                # Try to look at suffix for type
                extension = os.path.splitext(path)[1][1:]
                if extension[:2] == "tt" or extension[:2] == 'pertype':
                    per_type = int(extension[2])
                elif extension == 'per':
                    per_type = 1
                else:
                    # Default to 3
                    per_type = 3
    
        if no_data_value is None:
            no_data_value = self.no_data_value
    
        # Check to see if masks have been applied to sediment, if so use them
        # if masked is True
        if isinstance(self.Z, numpy.ma.MaskedArray) and masked:
            pass
        else:
            pass
    
        with open(path, 'w') as outfile:
            if self.unstructured:
                for (i, per) in enumerate(self.per):
                    outfile.write("%s %s %s\n" % (self.x[i], self.y[i], per))
        
            elif per_type == 1:
                for j in range(len(self.y)-1, -1, -1):
                    latitude = self.y[j]
                    for (i, longitude) in enumerate(self.x):
                        outfile.write("%s %s %s\n" % (longitude, latitude, self.per[j,i,:]))
        
            elif per_type == 2 or per_type == 3:
                # Write out header
                outfile.write('%6i                              ncols\n' % self.Per.shape[1])
                outfile.write('%6i                              nrows\n' % self.Per.shape[0])
                outfile.write('%6i                              nclasses\n' % self.Per.shape[3])
                outfile.write('%22.15e              xlower\n' % self.extent[0])
                outfile.write('%22.15e              ylower\n' % self.extent[2])
                outfile.write('%22.15e              cellsize\n' % self.delta)
                outfile.write('%10i                          nodata_value\n' % no_data_value)
            
                masked_Per = isinstance(self.Per, numpy.ma.MaskedArray)
            
            # Write out Sediment data
                if per_type == 2:
                    if masked_Per:
                        Per_filled = numpy.flipud(self.Per.filled())
                    else:
                        Per_filled = numpy.flipud(self.Per)
                    for i in xrange(self.Per.shape[0]):
                        for j in xrange(self.Per.shape[1]):
                            outfile.write("%22.15e\n" % Per_filled[i,j,:])
                    if masked_Per:
                        del Per_filled
                elif per_type == 3:
                    if masked_Per:
                        Per_flipped = numpy.flipud(self.Per.filled())
                    else:
                        Per_flipped = numpy.flipud(self.Per)
                    for i in xrange(self.Per.shape[0]):
                        for j in xrange(self.Per.shape[1]):
                            outfile.write("%22.15e   " % (Per_flipped[i,j]))
                            outfile.write("\n")
                    if masked_Per:
                        del Per_flipped
        
            else:
                raise NotImplemented("Output type %s not implemented." % per_type)










