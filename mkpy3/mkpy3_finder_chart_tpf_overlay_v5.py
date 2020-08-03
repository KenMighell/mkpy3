#!/usr/bin/env python3

# Kenneth John Mighell
# Kepler Support Scientist
# Kepler / K2 Science Office
# NASA Ames Research Center / SETI Institute


def mkpy3_finder_chart_tpf_overlay_v5(
  ax=None,
  survey_wcs=None,
  tpf=None,
  frame=None,
  colors=[None,'cornflowerblue','red'],
  lws=[0,3,4],
  zorders=[0,1,2],
  verbose=None
):
    """
Function: mkpy3_finder_chart_tpf_overlay_v5()

Purpose: 
    
Plot the aperture overlay on top of the sky survey image.

Parameters
----------
ax : matplotlib axes.Axes object
    axis object
survey_wcs : 
    World Coordinate System from FITS survey image HDU
tpf : lightkurve tpf object (optional)
    a lightkurve object of a Kepler/K2/TESS Target Pixel File
frame : int (optional)
    frame number of the Target Pixel File [starting at zero]
colors : 3-item list of color names [Matplotlib] (optional)
    Default: [None, 'cornflowerblue', 'red']
lws : 3-item list of line widths [Matplotlib] (optional)
    Default: [0,3,4]
zorders : 3-item list of zorder values (ioptional)
    Default: [0,1,2]
verbose : bool (optional)
    if True, print extra information

Returns: nothing

# Kenneth John Mighell
# Kepler/K2 Support Scientist
# Kepler/K2 Science Office
# NASA Ames Research Center / SETI Institute
    """
    func_ = 'mkpy3_finder_chart_tpf_overlay_v5'
    date_ = '2020AUG03'
    version_ = 'xe'
    #
    import numpy as np
    import astropy.units as u
    import ntpath
    import os
    #
    import mkpy3_tpf_get_coordinates_v1 as km1
    #
    assert(ax is not None)
    assert(tpf is not None)
    if (frame is None): frame = 0
    if (verbose is None): verbose = False
    if (verbose):
        print(frame,'=frame')
        print(colors,'=colors')
        print(lws,'=lws')
        print(zorders,'=zorders')
        print(verbose,'=verbose')
        print(ntpath.basename(tpf.path),'<--- TPF filename')
        print(os.path.dirname(tpf.path),'<--- TPF dirname')
    pass#if
    #
    # ===== add overlay to plot =====
    tpf_data = tpf.flux[frame]  # get frame data
    #
    # determine which pixels have data (not nans) or are in aperture mask
    # valid values: 0 = no data (nans), 1 = data, 2 = mask
    d = np.zeros(tpf_data.size,dtype=int)
    d[np.isfinite(tpf_data).flatten()] += 1
    d[tpf.pipeline_mask.flatten()] += 1
    #
    #**********
    #-----
    """
    # original method for lightkurve.__version__ <= 2.0.dev
    # get the RA,DEC values for the pixel centers using the tpf.get_coordinates() method
    pxrav  = tpf.get_coordinates()[0][frame].flatten()  # pixel RA vector
    pxdecv = tpf.get_coordinates()[1][frame].flatten()  # pixel DEC vector
    # ^--- both will fail if the version of lightkurve used has the tpf.get_coordinates bug
    # failure mode: the pixel RA,DEC values are not correct: off by one row and one column
    """
    #-----
    """
    # get the RA,DEC values for the pixel centers using the km1.mkpy3_tpf_get_coordinates_v1() function
    pxrav  = km1.mkpy3_tpf_get_coordinates_v1(tpf=tpf,recreate_bug=True)[0][frame].flatten()  # pixel RA vector
    pxdecv = km1.mkpy3_tpf_get_coordinates_v1(tpf=tpf,recreate_bug=True)[1][frame].flatten()  # pixel DEC vector
    # ^--- both will fail as the tpf.get_coordinates bug is recreated
    # failure mode: the pixel RA,DEC values are not correct: off by one row and one column
    """
    #-----
    # this will work even if using lightkurve.__version__ <= 2.0.dev
    # get the RA,DEC values for the pixel centers using the km1.mkpy3_tpf_get_coordinates_v1() function
    pxrav  = km1.mkpy3_tpf_get_coordinates_v1(tpf=tpf)[0][frame].flatten()  # pixel RA vector
    pxdecv = km1.mkpy3_tpf_get_coordinates_v1(tpf=tpf)[1][frame].flatten()  # pixel DEC vector
    # ^--- both should succeed
    #-----
    #**********
    #
    #==========================================================================
    #
    # See comments by Keaton Bell: 
    # https://github.com/KeplerGO/lightkurve/issues/14
    #
    # convert RA,DEC to pixel coordinates of the *survey* image
    origin0 = 0
    pixels = survey_wcs.wcs_world2pix(pxrav*u.degree, pxdecv*u.degree, origin0)
    # wcs_world2pix documentation: origin=0 (ZERO) when using Numpy ---------^ 
    #
    xpx = pixels[0]  # useful alias
    ypx = pixels[1]  # useful alias
    #
    # reshape for plotting
    xy = np.reshape(pixels,(2,pixels[0].size)).T
    npixels = len(xy)
    #
    # compute median offsets [in *survey* pixels] between TPF pixels
    dx = np.nanmedian(np.diff(xpx))
    dy = np.nanmedian(np.diff(ypx))        
    #
    # define locations of corners relative to pixel centers
    corners = np.array([[1.,1.],[1.,-1.],[-1.,-1.],[-1.,1],[1.,1.]])        
    #
    # offsetmatrix is a rotation/scaling matrix:
    offsetmatrix = np.array(((dx,-dy), (dy, dx)))/2.
    #      dx=cosine(theta) ---^         ^--- dy=sine(theta)  
    # where theta is the rotation angle of offsetmatrix
    for i in range(len(corners)):
        corners[i] = np.cross(offsetmatrix,corners[i])
    #
    # plot boundaries of each pixel
    for i in range(npixels):
        ccoords = xy[i]+corners
        k = d[i]
        c = colors[k]
        lw = lws[k]
        zorder = zorders[k]
        xxx = ccoords[:,0]
        yyy = ccoords[:,1]
        ax.plot(xxx, yyy, c=c, lw=lw, zorder=zorder)
    pass#for
    #==========================================================================
    #
    if (verbose): 
        cadenceno = tpf.cadenceno[frame]
        print()
        print('%d =cadenceno  <---  %d=frame' % (cadenceno,frame))
        print()
        print('%s %s %s' % (func_,date_,version_))
    pass#if
pass#def


if (__name__ == '__main__'):
    import matplotlib.pyplot as plt
    import astropy.units as u
    import os
    import ntpath

    import lightkurve as lk
    lk.log.setLevel('INFO')
    
    import mkpy3_finder_chart_survey_fits_image_get_v1 as km1
    import mkpy3_finder_chart_image_show_v1 as km2
    
    # Exoplanet Kelper-138b is "KIC 7603200":
    tpf = lk.search_targetpixelfile(target='kepler-138b',mission='kepler',\
      quarter=10).download(quality_bitmask=0)  
    print('TPF filename:', ntpath.basename(tpf.path))
    print('TPF dirname: ', os.path.dirname(tpf.path))
    
    target = 'Kepler-138b'
    title_ = tpf.hdu[0].header['object']
    title_ += ' : Exoplanet '+target

    ra_deg = tpf.ra
    dec_deg = tpf.dec
    
    # get survey image data
    width_height_arcmin = 2.00
    survey = '2MASS-J'
    survey_hdu, survey_hdr, survey_data, survey_wcs, survey_cframe = \
      km1.mkpy3_finder_chart_survey_fits_image_get_v1(ra_deg, dec_deg, 
      radius_arcmin=width_height_arcmin, survey=survey, verbose=True)
    
    # create a matplotlib figure object
    fig = plt.figure(figsize=(12,12));
    
    # create a matplotlib axis object with right ascension and declination axes
    ax = plt.subplot(projection=survey_wcs)
            
    # show the survey image
    percentile = 99.0
    km2.mkpy3_finder_chart_image_show_v1(ax=ax, image_data=survey_data, \
      percentile=percentile, verbose=True)
    
    # replace title with a custom title
    plt.suptitle(title_,size=24)
    
    # show the TPF overlay                                                                                                                            
    frame = 0
    mkpy3_finder_chart_tpf_overlay_v5(ax=ax, survey_wcs=survey_wcs, tpf=tpf,\
      frame=frame, verbose=True)
    
    # put a yellow circle at the target position
    ax.scatter(ra_deg*u.deg, dec_deg*u.deg, \
      transform=ax.get_transform(survey_cframe), \
      s=600, edgecolor='yellow', facecolor='None', lw=3, zorder=10);
    
    # adjust the plot margins
    plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)

    pname = 'mkpy3_plot.png'
    if (pname != ''): 
        plt.savefig(pname, bbox_inches = "tight")
        print()
        print(pname,' <--- plot filename has been written.\n')
    pass#if
pass#if
#EOF
