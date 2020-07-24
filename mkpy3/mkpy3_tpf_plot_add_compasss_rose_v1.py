#!/usr/bin/env python3

# Kenneth John Mighell
# Kepler/K2 Support Scientist
# Kepler/K2 Science Office
# NASA Ames Research Center / SETI Institute


def mkpy3_tpf_plot_add_compass_rose_v1(\
  ax = None,
  tpf = None,                                       
  cx = None,
  cy = None,                               
  north_arrow_length_arcsec = None,
  edge_color = None,
  inside_color = None,
  edge_lw = None,
  inside_lw = None,                                      
  verbose=None
):
    """
Function : mkpy3_tpf_plot_add_compass_rose_v1()

Purpose: Add a compass rose to a tpf.plot() graph.  
    NOTE: Long arm points North (increasing declination)
    NOTE: Short arm points East (increasing right ascension)

Parameters
----------
ax : matplotlib Axes object
tpf : lightkurve Kepler/K2/TESS TPF object
cx : float (optional)
    pixel column number (X position) of the center of the compass rose 
cy : float (optional)
    pixel row number (Y position) of the center of the compass rose 
north_arrow_length_arcsec : float (optional) [default: 6 arcsec]
    length of the North arrow arm of the compass rose in arcsec 
    [N.B. A Kepler photometer CCD pixel is 3.98 arcsec/pixel].
    [N.B. A TESS photometer CCD pixel is 21.00 arcsec/pixel].
edge_color : matplotlib color name (optional) [default: 'blue']
    color of the edge of the compass rose 
inside_color : matplotlib color name (optional) [default: 'yellow']
    color of the inside of the compass rose 
edge_lw : float (optional) [default: 7]
    line width of the edge of the compass rose
inside_lw : float (optional) [default: 4]
    line width of the inside of the compass rose
verbose : bool (optional)
    if True, print extra information

Returns: nothing 

Example:

#==========
import matplotlib.pyplot as plt 
import lightkurve as lk
#
tpf = lk.search_targetpixelfile(target='Kepler-138b',mission='Kepler',quarter=10).download()
#  ---> Exoplanet Kelper-138b is "KIC 7603200"
#
# Plot the 2nd frame of the TPF
ax = tpf.plot(frame=1)
#
# add a compass rose 
mkpy3_tpf_plot_add_compass_rose_v1(ax=ax, tpf=tpf)
#
plt.savefig('mkpy3_spud.png', bbox_inches="tight")
#==========

# Kenneth John Mighell
# Kepler/K2 Support Scientist
# Kepler/K2 Science Office
# NASA Ames Research Center / SETI Institute
    """
    import matplotlib.pyplot as plt 
    import numpy as np

    assert(tpf is not None)
    assert(ax is not None)
    if (cx is None):
        xlim = ax.get_xlim()
        cx = (xlim[0]+xlim[1])/2.0  # center of the X axis
    if (cy is None):
        ylim = ax.get_ylim()
        cy = (ylim[0]+ylim[1])/2.0  # center of the Y axis   
    if (north_arrow_length_arcsec is None): 
        north_arrow_length_arcsec = 6  # good default for Kepler/K2 observations
    if (edge_color is None): edge_color = 'blue'
    if (inside_color is None): inside_color = 'yellow'
    if (edge_lw is None): edge_lw = 7
    if (inside_lw is None): inside_lw = 4
    if (verbose is None): verbose = False
    if (verbose):
        print(ax,'=ax')
        print(tpf,'=tpf')
        print(cx,'=cx')
        print(cy,'=cy')
        print(north_arrow_length_arcsec,'=north_arrow_length_arcsec')
        print(edge_color,'=edge_color')
        print(inside_color,'=inside_color')
        print(edge_lw,'=edge_lw')
        print(inside_lw,'=inside_lw')
        print(verbose,'=verbose')
    pass#if

    wcs = tpf.wcs
    cx0 = cx
    cy0 = cy
    north_arrow_length_deg = north_arrow_length_arcsec/3600.0 
    east_arrow_length_deg  = north_arrow_length_deg/2.0  # east arm is half as long as north arm

    # north arm of compass rose
    pixcrd0 = np.array( [[cx0,cy0]], dtype=np.float_)  # KJM: N.B. pixcrd0 must be a numpy 2-d array
    world = wcs.wcs_pix2world(pixcrd0, 0)  # pixels --> right ascension and declination
    world[0][1] += north_arrow_length_deg
    pixcrd1 = wcs.wcs_world2pix(world, 0)  # right ascension and declination --> pixels
    n_x0 = pixcrd0[0][0]
    n_y0 = pixcrd0[0][1]
    n_x1 = pixcrd1[0][0]
    n_y1 = pixcrd1[0][1]
    
    # east arm of compass rose
    pixcrd0 = np.array( [[cx0,cy0]], dtype=np.float_)  # KJM: N.B. pixcrd0 must be a numpy 2-d array
    world = wcs.wcs_pix2world(pixcrd0, 0)  # pixels --> right ascension and declination
    declination = world[0][1]
    world[0][0] += east_arrow_length_deg/np.cos(np.deg2rad(declination))
    pixcrd1 = wcs.wcs_world2pix(world, 0)  # right ascension and declination --> pixels
    e_x0 = pixcrd0[0][0]
    e_y0 = pixcrd0[0][1]
    e_x1 = pixcrd1[0][0]
    e_y1 = pixcrd1[0][1]
    
    # draw edge of compass rose with thick lines
    line = plt.Line2D( (n_x0,n_x1), (n_y0,n_y1), lw=edge_lw, color=edge_color) 
    ax.add_line(line)
    line = plt.Line2D( (e_x0,e_x1), (e_y0,e_y1), lw=edge_lw, color=edge_color) 
    ax.add_line(line)
    
    # draw middle of compass rose with thin lines
    line = plt.Line2D((n_x0,n_x1), (n_y0,n_y1), lw=inside_lw, color=inside_color) 
    ax.add_line(line)
    line = plt.Line2D((e_x0,e_x1), (e_y0,e_y1), lw=inside_lw, color=inside_color)
    ax.add_line(line)

    return    
pass#def


if (__name__ == '__main__'): 
    import matplotlib.pyplot as plt 
    import lightkurve as lk
    #
    tpf = lk.search_targetpixelfile(target='Kepler-138b',mission='Kepler',quarter=10).download()
    #  ---> Exoplanet Kelper-138b is "KIC 7603200"
    #
    # Plot the 2nd frame of the TPF
    ax = tpf.plot(frame=1)
    #
    # add a compass rose 
    mkpy3_tpf_plot_add_compass_rose_v1(ax=ax, tpf=tpf)
    #
    plot_file = 'mkpy3_plot.png'
    plt.savefig(plot_file, bbox_inches="tight")
    plt.close()
    print(plot_file,' <--- new PNG file written')
pass#if
#EOF
