#!/usr/bin/env python
# encoding: utf-8

    """GeoClaw Sediments thickness Tools Module
    
    Module provides several functions for reading, writing and ploting
    sediments thickness (deposit) files.
    
    :Functions:
    - sed1writer
    - sed2writer
    - sed3writer
    

    
    :Sediments Class:
    

    """
import os
import urllib
import types

import numpy

# ==============================================================================
#  Functions
# ==============================================================================

def sed1writer (outfile,thick,xlower,xupper,ylower,yupper,nxpoints,nypoints):
    """
    Function sed1writer will write out the sediment thickness profile by evaluating the
    function thickness on the grid specified by the other parameters.
    
    Assumes thickness can be called on arrays X,Y produced by numpy.meshgrid.
    
    Output file is of "sedtype1," which we use to refer to a file with
    (x,y,z) values on each line, progressing from upper left corner across
    rows, then down.
    """
    sediment = Sediment(Sed_func=thick)
    sediment.x = numpy.linspace(xlower,xupper,nxpoints)
    sediment.y = numpy.linspace(ylower,yupper,nypoints)
    sediment.write(outfile, sed_type=1)


def sed2writer (outfile,thick,xlower,xupper,ylower,yupper,nxpoints,nypoints):
    """
        Function sed2writer will write out the sediment thickness profile by evaluating the
        function thickness on the grid specified by the other parameters.
        This routine is here for simply creates a new sediment object and writes it out.
    """
    
    sediment = Sediment(sed_func=thick)
    sediment.x = numpy.linspace(xlower,xupper,nxpoints)
    sediment.y = numpy.linspace(ylower,yupper,nypoints)
    sediment.write(outfile, sed_type=2)

def sed3writer (outfile,thick,xlower,xupper,ylower,yupper,nxpoints,nypoints):
    """        
        Function sed2writer will write out the sediment thickness profile by evaluating the
        function thickness on the grid specified by the other parameters.
        This routine is here for simply creating a new sediment object and writes it out.
        
    """
    
    sediment = Sediment(sed_func=thick)
    sediment.x = numpy.linspace(xlower,xupper,nxpoints)
    sediment.y = numpy.linspace(ylower,yupper,nypoints)
    sediment.write(outfile, sed_type=3)

# Get sediment data directory maybe not needed in here, not sure now.

# ==============================================================================
#  Sediment class
# ==============================================================================
class Sediment(object):
    """Base sediment class.
        
        A class representing a single sediment thickness file
        
        :Properties:
        
        :Initialization:
        -
        
        :Examples:
        
        >>> import clawpack.geoclaw.sedthtools as sed
        >>> sed_file = sed.Sediment('./sed.tt3')
        >>> sed_file.plot()
        
    """
    @property
    def th(self):
        """A representation of the data as an 1d array."""
        if (self._th is None) and self.unstructured:
            self.read(mask=True)
        return self._th
    @th.setter
    def th(self, value):
        self._th = value
    @th.deleter
    def th(self):
        del self._th
    
    @property
    def Th(self):
        """A representation of the data as a 2d array."""
        if self._Th is None:
            self.generate_2d_thickness(mask=True)
        return self._Th
    @Th.setter
    def Th(self, value):
        self._Th = value
    @Th.deleter
    def Th(self):
        del self._Th
    
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
        """One dimensional coordinate array in y direction."""
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

    def __init__(self, path=None, sed_func=None, sed_type=None,
             unstructured=False):
        r"""Sediment initialization routine.
        
        See :class:`Sediment` for more info.
        
        """
    
        super(Sediment, self).__init__()
    
        self.path = path
        self.sed_func = sed_func
        self.sed_type = sed_type
    
        self.unstructured = unstructured
        self.no_data_value = -9999
    
        # Data storage for only calculating array shapes when needed
        self._th = None
        self._Th = None
        self._x = None
        self._X = None
        self._y = None
        self._Y = None
        self._extent = None
        self._delta = None
    
        self.coordinate_transform = lambda x,y: (x,y)

    def generate_2d_thickness(self, mask=True):
        r"""Generate a 2d array of the sediment thickness."""
    
        # Check to see if we need to generate these
        if self._Th is None:
        
            if self.unstructured:

                raise ValueError("Unstructured data does not allow for use of" \
                                 + " 2d arrays")
        
            if self.path is not None:
                if self._th is None:
                    # Try to read the data, may not have done this yet
                    self.read(path=self.path, mask=mask)
                    if self._Th is not None:
                        # We are done, the read function did our work
                        return
            
                    # See if self._X and self._Y are already computed and use them if
                    # available, otherwise just use self._x and self._y
                    if self._X is not None and self._Y is not None:
                        new_shape = self._X.shape
                    else:
                        new_shape = (self._x.shape[0], self._y.shape[0])
                    # Reshape, note that the mask follows along with the new array
                    self._Th = numpy.reshape(self._th, new_shape)
        
            elif self.sed_func is not None:
                # Generate sediment thickness profile via sed_func
                self._Th = self.sed_func(self.X, self.Y)

    def generate_2d_coordinates(self, mask=True):
        r"""Generate 2d coordinate arrays."""
    
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
                if abs(self.sed_type) == 1:
                    # Reading this sed_type should produce the X and Y arrays
                    self.read(mask=mask)
                elif abs(self.sed_type) in [2,3]:
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
            # self._TH and then self._th
            if mask:
                if self._Th is None:
                    # Check to see if we really need to do anything here
                    if isinstance(self._th, numpy.ma.MaskedArray):
                    # Try to create self._Z
                        self.generate_2d_thickness(mask=mask)
            
                if isinstance(self._Th, numpy.ma.MaskedArray):
                    # Use Th's mask for the X and Y coordinates
                    self._X = numpy.ma.MaskedArray(self._X, mask=self._Th.mask,
                                               copy=False)
                    self._Y = numpy.ma.MaskedArray(self._Y, mask=self._Th.mask,
                                                copy=False)



    def read(self, path=None, sed_type=None, unstructured=False,
             mask=True, filter_region=None, force=False):
        """Read in the data from the object's *path* attribute.
        
        Stores the resulting data in one of the sets of *x*, *y*, and *th* or
        *X*, *Y*, and *Th*.
        
        :Input:
        - *path* (str)  file to read
        - *sed_type* (int)
        - *unstructured* (bool)
        - *mask* (bool)
        The first three might have already been set when instatiating object.
        
        """
    
        if (path is None) and (self.path is None):
            raise ValueError("*** Need to set path for file to read")
    
        if path:
            self.path = path   # set or perhaps reset
            self.sed_type = None  # force resetting below
    
        if unstructured:
            self.unstructured = unstructured
        
        if self.sed_type is None:
            if sed_type is not None:
            self.sed_type = sed_type
            else:
                # Try to look at suffix for type
                extension = os.path.splitext(self.path)[1][1:]
                if extension[:2] == "tt":
                    self.sed_type = int(extension[2])
                elif extension == 'sed':
                    self.sed_type = 1
                elif extension == 'asc':
                    self.sed_type = 3
                else:
                    # Default to 3
                    self.sed_type = 3
    
        if self.unstructured:
            # Read in the data as series of tuples
            data = numpy.loadtxt(self.path)

            self._x = data[:,0]
            self._y = data[:,1]
            self._th = data[:,2]
    
        else:
            # Data is in one of the GeoClaw supported formats
            if abs(self.sed_type) == 1:
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
        



    def read_header(self):
        r"""Read in header of sediment file at path.
        
        If a value returns numpy.nan then the value was not retrievable.  Note
        that this routine can read in headers whose values and labels are
        swapped.
        
        """
    
        if abs(self.sed_type) in [2,3]:
        
            # Default values to track errors
            num_cells = [numpy.nan,numpy.nan]
            self._extent = [numpy.nan,numpy.nan,numpy.nan,numpy.nan]
            self._delta = numpy.nan
        
            with open(self.path, 'r') as sed_file:
                # Check to see if we need to flip the header values
                first_line = sed_file.readline()
                try:
                    num_cells[0] = int(first_line.split()[0])
                except ValueError:
                    # Assume the header is flipped from what we expect
                    num_cells[0] = int(first_line.split()[-1])
                    value_index = -1
                else:
                    value_index = 0
            
                num_cells[1] = int(sed_file.readline().split()[value_index])
                self._extent[0] = float(sed_file.readline().split()[value_index])
                self._extent[2] = float(sed_file.readline().split()[value_index])
                self._delta = float(sed_file.readline().split()[value_index])
                self.no_data_value = float(sed_file.readline().split()[value_index])
            
                self._extent[1] = self._extent[0] + num_cells[0] * self.delta
                self._extent[3] = self._extent[2] + num_cells[1] * self.delta
    
        else:
            raise IOError("Cannot read header for sed_type %s" % self.sed_type)

        return num_cells


    def write(self, path, no_data_value=None, sed_type=None, masked=True):
        """Write out a Sediment file to path of type *sed_type*.
        
        Writes out a Sediment file of sed type specified with *sed_type* or
        inferred from the output file's extension, defaulting to 3, to path
        from data in Th.  The rest of the arguments are used to write the header
        data.
        
        """
    
        # Determine sed type if not specified
        if sed_type is None:
            if self.sed_type is not None:
                sed_type = self.sed_type
            else:
                # Try to look at suffix for type
                extension = os.path.splitext(path)[1][1:]
                if extension[:2] == "tt" or extension[:2] == 'sedtype':
                    sed_type = int(extension[2])
                elif extension == 'sed':
                    sed_type = 1
                else:
                    # Default to 3
                    sed_type = 3
    
        # Default to this object's no_data_value if the passed is None,
        # otherwise the argument will override the object's value or it will
        # default to -9999 (default for the class)
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
                for (i, sed) in enumerate(self.th):
                    outfile.write("%s %s %s\n" % (self.x[i], self.y[i], th))
        
            elif sed_type == 1:
                # longitudes = numpy.linspace(lower[0], lower[0] + delta * Z.shape[0], Z.shape[0])
                # latitudes = numpy.linspace(lower[1], lower[1] + delta * Z.shape[1], Z.shape[1])
                for j in range(len(self.y)-1, -1, -1):
                    latitude = self.y[j]
                    for (i, longitude) in enumerate(self.x):
                        outfile.write("%s %s %s\n" % (longitude, latitude, self.Th[j,i]))
        
            elif sed_type == 2 or sed_type == 3:
                # Write out header
                outfile.write('%6i                              ncols\n' % self.Th.shape[1])
                outfile.write('%6i                              nrows\n' % self.Th.shape[0])
                outfile.write('%22.15e              xlower\n' % self.extent[0])
                outfile.write('%22.15e              ylower\n' % self.extent[2])
                outfile.write('%22.15e              cellsize\n' % self.delta)
                outfile.write('%10i                          nodata_value\n' % no_data_value)
            
                masked_Th = isinstance(self.Th, numpy.ma.MaskedArray)
            
                # Write out Sediment data
                if sed_type == 2:
                    if masked_Th:
                        Th_filled = numpy.flipud(self.Th.filled())
                    else:
                        Th_filled = numpy.flipud(self.Th)
                    for i in xrange(self.Th.shape[0]):
                        for j in xrange(self.Th.shape[1]):
                            outfile.write("%22.15e\n" % Th_filled[i,j])
                    if masked_Th:
                        del Th_filled
                elif sed_type == 3:
                    if masked_Th:
                        Th_flipped = numpy.flipud(self.Th.filled())
                    else:
                        Th_flipped = numpy.flipud(self.Th)
                    for i in xrange(self.Th.shape[0]):
                        for j in xrange(self.Th.shape[1]):
                            outfile.write("%22.15e   " % (Th_flipped[i,j]))
                        outfile.write("\n")
                    if masked_Th:
                        del Th_flipped
        
            else:
                raise NotImplemented("Output type %s not implemented." % sed_type)



































