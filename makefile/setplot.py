
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 
from pylab import *

try:
    from setplotfg import setplotfg
except:
    print "Did not find setplotfg.py"
    setplotfg = None


#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps, geoplot
    

    plotdata.clearfigures()  # clear any old figures,axes,items data

    def set_drytol(current_data):
        # The drytol parameter is used in masking land and water and
        # affects what color map is used for cells with small water depth h.
        # The cell will be plotted as dry if h < drytol.
        # The best value to use often depends on the application and can
        # be set here (measured in meters):
        current_data.user['drytol'] = 1.e-2

    plotdata.beforeframe = set_drytol

    # To plot gauge locations on pcolor or contour plot, use this as
    # an afteraxis function:

    def addgauges(current_data):
        from clawpack.visclaw import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)
    

    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='pcolor', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Surface'
    plotaxes.scaled = True

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -20
    plotitem.pcolor_cmax = 20
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [1,1,0]
    plotitem.amr_patchedges_show = [1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [1,1,0]
    plotaxes.xlimits = [0,200000]
    plotaxes.ylimits = [0,200000]



    #-----------------------------------------
    # Figure for zoom
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Zoom', figno=10)
    #plotfigure.show = False
    plotfigure.kwargs = {'figsize':[12,7]}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('diag zoom')
    plotaxes.axescmd = 'axes([0.0,0.1,0.6,0.6])'
    plotaxes.title = 'On diagonal'
    plotaxes.scaled = True
    plotaxes.xlimits = [55,66]
    plotaxes.ylimits = [55,66]
    plotaxes.afteraxes = addgauges

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -2.9
    plotitem.pcolor_cmax = 2.9
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [1,1,0]
    plotitem.amr_patchedges_show = [1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [1,1,0]

    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = arange(-10., 0., 1.)
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,0,1]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True

    # Add contour lines of topography:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = arange(0., 11., 1.)
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,0,1]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True

    # Add dashed contour line for shoreline
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = [0.]
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'dashed'}
    plotitem.amr_contour_show = [0,0,1]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True



    #-----------------------------------------
    # Figure for zoom near axis
    #-----------------------------------------
    #plotfigure = plotdata.new_plotfigure(name='Zoom2', figno=11)
    # now included in same figure as zoom on diagonal

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('x zoom')
    plotaxes.show = True
    plotaxes.axescmd = 'axes([0.5,0.1,0.6,0.6])'
    plotaxes.title = 'On x-axis'
    plotaxes.scaled = True
    plotaxes.xlimits = [82,93]
    plotaxes.ylimits = [-5,6]
    plotaxes.afteraxes = addgauges

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    #plotitem.plot_var = geoplot.surface
    plotitem.plot_var = geoplot.surface_or_depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = -0.9
    plotitem.pcolor_cmax = 0.9
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [1,1,0]
    plotitem.amr_patchedges_show = [1]

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [1,1,0]


    # Add contour lines of bathymetry:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = arange(-10., 0., 1.)
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,0,1]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True

    # Add contour lines of topography:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    from numpy import arange, linspace
    plotitem.contour_levels = arange(0., 11., 1.)
    plotitem.amr_contour_colors = ['g']  # color on each level
    plotitem.kwargs = {'linestyles':'solid'}
    plotitem.amr_contour_show = [0,0,1]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True

    # Add dashed contour line for shoreline
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = geoplot.topo
    plotitem.contour_levels = [0.]
    plotitem.amr_contour_colors = ['k']  # color on each level
    plotitem.kwargs = {'linestyles':'dashed'}
    plotitem.amr_contour_show = [0,0,1]  # show contours only on finest level
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 0
    plotitem.show = True



    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Surface & topo', figno=300, \
                    type='each_gauge')

    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-2.0, 2.0]
    plotaxes.title = 'Surface'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3
    plotitem.plotstyle = 'b-'

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo
        
    plotitem.plot_var = gaugetopo
    plotitem.plotstyle = 'g-'
    def add_zeroline(current_data):
        from pylab import plot, legend
        t = current_data.t
        legend(('surface','topography'),loc='lower left')
        plot(t, 0*t, 'k')

    plotaxes.afteraxes = add_zeroline


    #-----------------------------------------
    # Figure for patches alone
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='patches', figno=2)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,1]
    plotaxes.ylimits = [0,1]
    plotaxes.title = 'patches'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.amr_celledges_show = [1,1,0]   
    plotitem.amr_patchedges_show = [1]     

    #-----------------------------------------
    # Scatter plot of surface for radially symmetric
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Scatter', figno=200)
    plotfigure.show = False

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0., 100.]
    plotaxes.ylimits = [-.5, 1.]
    plotaxes.title = 'Scatter plot of surface'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    plotitem.plot_var = geoplot.surface
    def q_vs_radius(current_data):
        from numpy import sqrt
        x = current_data.x
        y = current_data.y
        r = sqrt(x**2 + y**2)
        q = current_data.var
        return r,q
    plotitem.map_2d_to_1d = q_vs_radius
    plotitem.plotstyle = 'o'
    plotitem.amr_color=['b','r','g']
    plotaxes.afteraxes = "pylab.legend(['Level 1','Level 2'])"
    

    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = [4,5,104,105]  # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.format = 'ascii'                # Format of output
    # plotdata.format = 'netcdf'
    
    #-----------------------------------------
    # Figure for cross section
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='cross-section', figno=1)
    
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0,100000.e0]
    plotaxes.ylimits = [-10,100]
    plotaxes.title = 'Cross section at x=0'
    def plot_topo_xsec(current_data):
        from pylab import plot, hold, cos,sin,where,legend,nan
        t = current_data.t
        hold(True)
        #eta1 = sigma*h0/a**2 * (2.*x*cos(omega*t) + 2.*y*sin(omega*t) -sigma)
        #etatrue = where(eta1>B, eta1, nan)
        #plot(x, etatrue, 'r', label="true solution", linewidth=2)
        #plot(x, B, 'g', label="bathymetry")
        ## plot([0],[-1],'kx',label="Level 1")  # shouldn't show up in plots,
        ## plot([0],[-1],'bo',label="Level 2")  # but will produced desired legend
        #plot([0],[-1],'bo',label="Computed")  ## need to fix plotstyle
        legend()
        hold(False)
    plotaxes.afteraxes = plot_topo_xsec

    
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    
    def ysec(current_data):
        # Return x value and surface eta at this point, along y=0
        from pylab import find,ravel
        x = current_data.x
        y = current_data.y
        dx = current_data.dx
        q = current_data.q
        t = current_data.t
        xlower = current_data.xlower
        xupper = current_data.xupper
        
        ij = find((x <= dx/2.) & (x > -dx/2.))
        y_slice = ravel(y)[ij]
        eta_slice = ravel(q[3,:,:])[ij]
        #print 'data'
        #print y_slice
        #fname = 'surface%d'%int(t)+'.txt'
        fname = r'/Users/huitang/clawpack-5.2.0/geoclaw/examples/tsunami/bowl-radial/_output/'+'surface_t=%d'%int(t)+'_ylower=%d'%int(y[1,1])+'.txt'
        fname2 = r'/Users/huitang/clawpack-5.2.0/geoclaw/examples/tsunami/bowl-radial/_output/location'+'_ylower=%d'%int(y[1,1])+'.txt'
        fname3 = r'/Users/huitang/clawpack-5.2.0/geoclaw/examples/tsunami/bowl-radial/_output/'+'momentums_t=%d'%int(t)+'_ylower=%d'%int(y[1,1])+'.txt'
        fname4 = r'/Users/huitang/clawpack-5.2.0/geoclaw/examples/tsunami/bowl-radial/_output/'+'depth_t=%d'%int(t)+'_ylower=%d'%int(y[1,1])+'.txt'
        savetxt(fname3,ravel(q[2,:,:])[ij])
        savetxt(fname,ravel(q[3,:,:])[ij])
        savetxt(fname4,ravel(q[0,:,:])[ij])
        savetxt(fname2,ravel(y)[ij])
        print '***********************************************'
        #print y[1]
        print ravel(q[0,:,:])[ij]
        #print y_slice
        #print y[1,1]
        #print xlower
        #print eta_slice
        return y_slice, eta_slice
  
    plotitem.map_2d_to_1d = ysec
    #data = ysec(current_data)
    #print plotdata
    #savetxt('ysec.txt',ysec)
    plotitem.plotstyle = 'bo'     ## need to be able to set amr_plotstyle
    plotitem.kwargs = {'markersize':5}
    plotitem.amr_show = [1]  # plot on all levels
    

    return plotdata

    
