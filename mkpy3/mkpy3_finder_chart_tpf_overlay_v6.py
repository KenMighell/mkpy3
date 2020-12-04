#!/usr/bin/env python3

# file://mkpy3_finder_chart_tpf_overlay_v6.py

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute

"""
import sys
#
pyver = (sys.version_info.major*10) + (sys.version_info.minor)
if (pyver < 27):
    print('*** ERROR *** This application needs Python 2.7 or higher.')
    sys.exit(1)
# pass:if
#
try:
    import lightkurve as lk
except Exception:
    print('\n***** ERROR *****\n')
    print('The Python package lightkurve needs to be installed.\n')
    print('This is the installation command for lightkurve using pip:\n')
    print('pip install lightkurve --upgrade\n')
    print('For further installation details see the lightkurve homepage:\n')
    print('http://lightkurve.keplerscience.org/install.html\n')
    sys.exit(1)
# pass:try
"""


def mkpy3_finder_chart_tpf_overlay_v6(
  ax=None,
  survey_wcs=None,
  tpf=None,
  frame=None,
  colors=[None, 'cornflowerblue', 'red'],
  lws=[0, 3, 4],
  zorders=[0, 1, 2],
  verbose=None
):
    """
Function: mkpy3_finder_chart_tpf_overlay_v6()

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

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute
    """
    func_ = 'mkpy3_finder_chart_tpf_overlay_v6'
    date_ = '2020NOV23'
    version_ = 'xh'
    #
    import numpy as np
    import astropy.units as u
    import ntpath
    import os
    #
    import mkpy3
    #
    assert(ax is not None)
    assert(tpf is not None)
    if (frame is None):
        frame = 0
    # pass:if
    if (verbose is None):
        verbose = False
    # pass:if
    if (verbose):
        print(frame, '=frame')
        print(colors, '=colors')
        print(lws, '=lws')
        print(zorders, '=zorders')
        print(verbose, '=verbose')
        print(ntpath.basename(tpf.path), '<--- TPF filename')
        print(os.path.dirname(tpf.path), '<--- TPF dirname')
        print()
    # pass:if
    #
    # ===== add overlay to plot =====
    tpf_data = tpf.flux[frame]  # get frame data
    #
    # determine which pixels have data (not nans) or are in aperture mask
    # valid values: 0 = no data (nans), 1 = data, 2 = mask
    d = np.zeros(tpf_data.size, dtype=int)
    d[np.isfinite(tpf_data).flatten()] += 1
    d[tpf.pipeline_mask.flatten()] += 1
    #
    # =========================================================================
    # -----
    """
    # original method for lightkurve.__version__ <= 2.0.dev
    # get the RA,DEC values for the pixel centers using the
    # tpf.get_coordinates() method
    pxrav  = tpf.get_coordinates()[0][frame].flatten()  # pixel RA vector
    pxdecv = tpf.get_coordinates()[1][frame].flatten()  # pixel DEC vector
    # ^--- both will fail if the version of lightkurve used has the
    # tpf.get_coordinates bug
    # failure mode: the pixel RA,DEC values are not correct: off by one row
    # and one column
    """
    # -----
    """
    # get the RA,DEC values for the pixel centers using the
    # km1.mkpy3_tpf_get_coordinates_v1() function
    pxrav  = km1.mkpy3_tpf_get_coordinates_v1(
            tpf=tpf,recreate_bug=True)[0][frame].flatten()  # pixel RA vector
    pxdecv = km1.mkpy3_tpf_get_coordinates_v1(
            tpf=tpf,recreate_bug=True)[1][frame].flatten()  # pixel DEC vector
    # ^--- both will fail as the tpf.get_coordinates bug is recreated
    # failure mode: the pixel RA,DEC values are not correct: off by one row
    # and one column
    """
    # -----
    # this will work even if using lightkurve.__version__ <= 2.0.dev
    # get the RA,DEC values for the pixel centers:
    pxrav = mkpy3.mkpy3_tpf_get_coordinates_v1(tpf=tpf)[0][frame].flatten()
    # ^--- pixel RA vector
    pxdecv = mkpy3.mkpy3_tpf_get_coordinates_v1(tpf=tpf)[1][frame].flatten()
    # ^--- pixel DEC vector
    # ^--- both function calls should succeed
    # -----
    # =========================================================================
    #
    # See comments by Keaton Bell:
    # https://github.com/KeplerGO/lightkurve/issues/14
    #
    # convert RA,DEC to pixel coordinates of the *survey* image
    origin0 = 0
    pixels = survey_wcs.wcs_world2pix(pxrav * u.degree, pxdecv * u.degree, origin0)
    # wcs_world2pix documentation: origin=0 (ZERO) when using Numpy ---------^
    #
    xpx = pixels[0]  # useful alias
    ypx = pixels[1]  # useful alias
    #
    # reshape for plotting
    xy = np.reshape(pixels, (2, pixels[0].size)).T
    npixels = len(xy)
    #
    # compute median offsets [in *survey* pixels] between TPF pixels
    dx = np.nanmedian(np.diff(xpx))
    dy = np.nanmedian(np.diff(ypx))
    #
    # define locations of corners relative to pixel centers
    corners = np.array([[1., 1.], [1., -1.], [-1., -1.], [-1., 1], [1., 1.]])
    #
    # offsetmatrix is a rotation/scaling matrix:
    offsetmatrix = np.array(((dx, -dy), (dy, dx))) / 2.
    #      dx=cosine(theta) ---^         ^--- dy=sine(theta)
    # where theta is the rotation angle of offsetmatrix
    for i in range(len(corners)):
        corners[i] = np.cross(offsetmatrix, corners[i])
    #
    # plot boundaries of each pixel
    for i in range(npixels):
        ccoords = xy[i] + corners
        k = d[i]
        c = colors[k]
        lw = lws[k]
        zorder = zorders[k]
        xxx = ccoords[:, 0]
        yyy = ccoords[:, 1]
        ax.plot(xxx, yyy, c=c, lw=lw, zorder=zorder)
    # pass:for
    # =========================================================================
    #
    if (verbose):
        cadenceno = tpf.cadenceno[frame]
        print('%d =cadenceno  <---  %d=frame' % (cadenceno, frame))
        print('\n%s %s %s  :-)' % (func_, date_, version_))
    # pass:if
# pass:def


def xmkpy3_finder_chart_tpf_overlay_v6():
    import matplotlib.pyplot as plt
    import astropy.units as u
    import sys
    import os
    import ntpath
    import argparse
    import ast

    import lightkurve as lk  # ignore PEP8 warning of redefinition
    lk.log.setLevel('INFO')

    import mkpy3

    #
    # argparse: BEGIN =========================================================
    #
    parser = argparse.ArgumentParser()
    #
    parser.add_argument(
        '--tpf_filename', action="store", type=str, default=None,
        help="Filename of the Target Pixel File (TPF) [default: None]")
    parser.add_argument(
        '--frame', action="store", type=int, default=0,
        help='Frame number (integer) [default: 0]')
    parser.add_argument(
        '--survey', action="store", type=str, default='2MASS-J',
        help="Survey name (str) [default: '2MASS-J']")
    parser.add_argument(
        '--width_height_arcmin', action="store", type=float, default=2.0,
        help='Width and height size in arcmin (float) [default: 2.0]')
    parser.add_argument(
        '--show_plot', type=mkpy3.mkpy3_util_str2bool, default=True,
        help='If True, show the plot [default=True]')
    parser.add_argument(
        '--plotfile', action="store", type=str, default='mkpy3_plot.png',
        help="Filename of the output plotfile [default: 'mkpy3_plot.png']")
    parser.add_argument(
        '--overwrite', type=mkpy3.mkpy3_util_str2bool, default=False,
        help='If True, overwrite ("clobber") an existing output file '
        '[default: False.')
    parser.add_argument(
        '--figsize', action="store", type=ast.literal_eval, default="[9,9]",
        help="2-item list of figure width and height [Matplotlib] (str) "
        "[default: \"[9,9]\"")
    parser.add_argument(
        '--title', action="store", type=str, default=None,
        help='Title of the finder chart (str) [default: None]')
    parser.add_argument(
        '--percentile', action="store", type=float, default=99.0,
        help='Percentile [percentage of pixels to keep: 0.0 to 100.0] '
        '(float) [default: 99.0]')
    parser.add_argument(
        '--cmap', action="store", type=str, default=None,
        help="Colormap name [Matplotlib] (str) [default: 'gray']")
    parser.add_argument(
        '--colors', action="store", type=ast.literal_eval,
        default="[None,'cornflowerblue','red']",
        help="3-item list of overlay color names [Matplotlib] (str) "
        "[default: \"['None','cornflowerblue','red']\"")
    parser.add_argument(
        '--lws', action="store", type=ast.literal_eval, default="[0,3,4]",
        help="3-item list of overlay line widths [Matplotlib] (str) "
        "[default: \"[0,3,4]\"")
    parser.add_argument(
        '--zorders', action="store", type=ast.literal_eval, default="[0,1,2]",
        help="3-item list of overlay zorder values [Matplotlib] (str) "
        "[default: \"[0,1,2]\"")
    parser.add_argument(
        '--marker_dict', action="store", type=ast.literal_eval,
        default="{'edgecolor':'yellow', 's':600, 'facecolor':'None', 'lw':3, "
        "'zorder':10}",
        help="marker kwargs (dictonary string) for ax.scatter() [Matplotlib] "
        "(str) [default: \"{'edgecolor':'yellow', 's':600, 'facecolor':'None',"
        " 'lw':3, 'zorder':10}\"")
    parser.add_argument(
        '--verbose', type=mkpy3.mkpy3_util_str2bool, default=False,
        help='Print extra information if True (bool) [default=False]')
    #
    args = parser.parse_args()

    tpf_filename = args.tpf_filename
    frame = args.frame
    survey = args.survey
    width_height_arcmin = args.width_height_arcmin
    show_plot = args.show_plot
    plotfile = args.plotfile
    overwrite = args.overwrite
    figsize = args.figsize
    title_ = args.title
    percentile = args.percentile
    cmap = args.cmap
    colors = args.colors
    lws = args.lws
    zorders = args.zorders
    marker_dict = args.marker_dict
    verbose = args.verbose
    #
    if (verbose):
        print('%s =args.tpf_filename' % (args.tpf_filename))
        print('%s =args.frame' % (args.frame))
        print("'%s' =args.survey" % (args.survey))
        print('%s =args.width_height_arcmin' % (args.width_height_arcmin))
        print('%s =args.show_plot' % (args.show_plot))
        print("'%s' =args.plotfile" % (args.plotfile))
        print('%s =args.overwrite' % (args.overwrite))
        print('%s =args.figsize' % (args.figsize))
        print('%s =args.title' % (args.title))
        print('%s =args.percentile' % (args.percentile))
        print('%s =args.cmap' % (args.cmap))
        print('%s =args.colors' % (args.colors))
        print('%s =args.lws' % (args.lws))
        print('%s =args.zorders' % (args.zorders))
        print('%s =args.marker_dict' % (args.marker_dict))
        print('%s =args.verbose' % (args.verbose))
        print()
    # pass:if
    #
    # argparse: END ===========================================================
    #

    ok = (type(marker_dict) is dict) or (marker_dict is None)
    if (not ok):
        print()
        print('**** ERROR ***** BAD ARGUMENT VALUE *****')
        print()
        print('marker_dict must be a dictionary or None:')
        print(marker_dict, '=marker_dict')
        print()
        sys.exit(1)
    # pass:if

    if (tpf_filename is not None):
        mkpy3.mkpy3_util_check_file_exists(tpf_filename, True)
        tpf = lk.read(tpf_filename)
    else:
        tpf = lk.search_targetpixelfile(
            target='kepler-138b', mission='kepler', cadence='long',
            quarter=10).download(quality_bitmask=0)
        #   6--- exoplanet Kelper-138b is "KIC 7603200"
        print()
        print('No TargetPixelFile (TPF) filename given.')
        print()
        print('Using default TPF [Kepler Q10 observations of exoplanet Kepler-'
              '138b (KIC 760320)]:')
    # pass#if
    print()
    print('TPF filename:', ntpath.basename(tpf.path))
    print('TPF dirname: ', os.path.dirname(tpf.path))
    print()

    ra_deg = tpf.ra
    dec_deg = tpf.dec
    if (verbose):
        print()
        print(ra_deg, '=ra_deg')
        print(dec_deg, '=dec_deg')
    # pass#if

    print()

    # get survey image data
    # survey = '2MASS-J'  # hard-wired option
    survey_hdu, survey_hdr, survey_data, survey_wcs, survey_cframe = \
        mkpy3.mkpy3_finder_chart_survey_fits_image_get_v1(
            ra_deg, dec_deg, radius_arcmin=width_height_arcmin,
            survey=survey, verbose=verbose)

    # create a matplotlib figure object
    fig = plt.figure(figsize=figsize)

    # create a matplotlib axis object with right ascension and declination axes
    ax = plt.subplot(projection=survey_wcs)

    # show the survey image
    mkpy3.mkpy3_finder_chart_image_show_v1(
        ax=ax, image_data=survey_data,
        percentile=percentile, cmap=cmap, verbose=verbose)

    # show the TPF overlay
    mkpy3_finder_chart_tpf_overlay_v6(
        ax=ax, survey_wcs=survey_wcs, tpf=tpf,
        frame=frame, colors=colors, lws=lws, zorders=zorders, verbose=verbose)

    # add title
    if (title_ is None):
        title_ = tpf.hdu[0].header['object']
    # pass:if
    plt.suptitle(title_, size=25)

    # put a yellow circle at the target position
    if (type(marker_dict) is dict):
        ax.scatter(
            ra_deg * u.deg, dec_deg * u.deg,
            transform=ax.get_transform(survey_cframe),
            **marker_dict)

    # adjust the plot margins
    plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)

    if (plotfile != ''):
        if (plotfile != 'mkpy3_plot.png'):
            mkpy3.mkpy3_util_check_file_exists(plotfile, overwrite)
        plt.savefig(plotfile, dpi=300)  # , bbox_inches = "tight")
        print('\n%s <--- plotfile written  :-)\n' % (plotfile))
    # pass:if

    if (show_plot):
        plt.ioff()
        plt.show()
    # pass:if

    plt.close()

# pass#def

# =============================================================================


if (__name__ == '__main__'):
    xmkpy3_finder_chart_tpf_overlay_v6()
# pass:if

# EOF
