#!/usr/bin/env python3

# file://mkpy3_tpf_overlay_v6.py

__version__ = '2020SEP29T1138  v0.39'

# Kenneth John Mighell
# Kepler Support Scientist
# Kepler / K2 Science Office
# NASA Ames Research Center / SETI Institute

###############################################################################


def check_file_exists(filename, overwrite):
    """
Utility function.
    """
    import sys
    import os
    assert(isinstance(filename, str))
    assert(isinstance(overwrite, bool))
    msg = 'Requested output file already exists (overwrite=False):\n'
    if (not overwrite):
        if (os.path.isfile(filename)):
            print('\n***** ERROR *****\n\n%s' % (msg))
            print("new_filename='%s'\n" % filename)
            sys.exit(1)
        # pass:if
    # pass:if
# pass:def

###############################################################################


def mkpy3_tpf_overlay_v6(
  tpf=None,
  frame=0,
  survey='2MASS-J',  # '2MASS-J' or 'DSS2 Red'
  width_height_arcmin=2.0,
  rotate_deg_str='None',
  shrink=1.0,
  show_plot=True,
  plot_file='mkpy3_plot.png',
  overwrite=False,
  figsize_str='[9,9]',
  title=None,
  percentile=99.0,
  cmap='gray_r',
  colors_str="[None,'cornflowerblue','red']",
  lws_str='[0,3,4]',
  zorders_str='[0,2,4]',
  marker_kwargs_str="{'edgecolor':'yellow', 's':600, 'facecolor':'None',\
    'lw':3, 'zorder':10}",  # or 'None'
  print_gaia_dr2=True,
  gaia_dr2_kwargs_str="{'edgecolor':'cyan', 's':300, 'facecolor':'None',\
    'lw':3, 'zorder':20}",  # or 'None'
  print_vsx=True,
  vsx_kwargs_str="{'s':900, 'color':'lawngreen', 'marker':'x', 'lw':5,\
    'zorder':30}",  # or 'None'
  sexagesimal=False,
  verbose=False
):
    """
Function: mkpy3_tpf_overlay_v6()

Purpose:

Plot a Kepler/K2/TESS TargetPixelFile (TPF) overlay on a sky survey image.

Parameters
----------
tpf : (lightkurve TargetPixelFile object) (optional)
    A lightkurve TargetPixelFile (TPF) object.
    [default: None]
frame : (int) (optional)
    Frame number to use.
    [range: 0 to number of cadences in the TPF minus 1]
    [default: 0]
survey : (str) (optional)
    A sky survey name.
    [default: '2MASS-J'] [verified: '2MASS-J', 'DSS2 Red']
rotate_deg_str : (str) (optional)
    Angle in degrees to rotate the sky survey image.
    [default: 'None']
    [examples: 'None' or '12.345' (a float) or "'tpf'" (a str)]
width_height_arcmin : (float) (optional)
    Width and height of the survey image [arcmin].
    [default: 2.0]
shrink : (float) (optional)
    Survey search radius shrink factor.
    [range: 0.0 to 1.0]
    [default: 1.0]
show_plot : (bool) (optional)
    If True, show the plot.
    [default=True]
plot_file : (str) (optional)
    Filename of the output plot.
    [default: 'mkpy3_plot.png']
overwrite : (bool) (optional)
    If True, overwrite ("clobber") an existing output file.
    If False, do *not* create output file when plot_file != 'mkpy3_plot.png'.
    [default: False]
figsize_str : (str) (optional)
    A string of a 2-time list of figure widht and height [Matplotlib].
    [default: '[9,9]']
title : (str) (optional)
    Title of the plot.
    If None, a title will be created.
    An empty string ('') will produce a blank title.
    [default: None]
percentile : (float) (optional)
    Percentile [percentage of pixels to keep] used to set the colorbar.
    [range: 0.0 to 100.0]
    [default: 99.0]
cmap : (str) (optional)
    Colormap name [Matplotlib].
    [default: 'gray_r']
colors_str : (str) (optional)
    A string of a 3-item list of overlay color names [Matplotlib].
    [default: "['None','cornflowerblue','red']"]
lws_str : (str) (optional)
    A string of a 3-item list of overlay line widths [Matplotlib].
    [default: '[0,3,4]']
zorders_str : (str) (optional)
    A string of a 3-item list of overlay zorder values [Matplotlib].
    [default: '[0,1,2]']
marker_kwargs_str : (str) (optional)
    A string of a dictionary of arguments for ax.scatter() [Matplotlib].
    The target is marked according to the kwarg values.
    If set to None, the target is *not* marked.
    [default: "{'edgecolor':'yellow', 's':600, 'facecolor':'None', 'lw':3,
    'zorder':10}"]
print_gaia_dr2 : (bool) (optional)
    If True, print the GAIA DR2 catalog results.
    [default=True]
gaia_dr2_kwargs_str : (str) (optional)
    A string of a dictionary of arguments for ax.scatter() [Matplotlib].
    GAIA DR2 stars are marked accordinbg to the kwarg values.
    If set to None, no GAIA DR2 data are shown and plotted.
    [default: "{'edgecolor':'cyan', 's':300, 'facecolor':'None', 'lw':3,
    'zorder':20}"]
print_vsx : (bool) (optional)
    If True, print the VSX catalog results.
vsx_kwargs_str : (str) (optional)
    A string of a dictionary of arguments for ax.scatter() [Matplotlib].
    VSX varaible stars are marked accordinbg to the kwarg values.
    If set to None, no VSX data are shown and plotted.
    [default: "{'s':900, 'color':'lawngreen', 'marker':'x', 'lw':5,
    'zorder':30}"]
sexagesimal : (bool) (optional)
    If True, print catalog positions as sexagesimal [hms dms].
    [default=False]
verbose : (bool) (optional)
    If True, print extra information.
    [default: False]

Returns
-------
ax : (matplotlib axes object) or (None)
    A matplotlib axes object *if* show_plot is False *else* None .

# Kenneth John Mighell
# Kepler Support Scientist
# Kepler / K2 Science Office
# NASA Ames Research Center / SETI Institute
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import copy as cp
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import ast
    import sys
    #
    import reproject as rp  # pip install reproject
    #
    import lightkurve as lk
    lk.log.setLevel('INFO')
    #
    import mkpy3_finder_chart_survey_fits_image_get_v1 as km1
    import mkpy3_finder_chart_image_show_v1 as km2
    import mkpy3_finder_chart_tpf_overlay_v6 as km3
    import mkpy3_vizier_gaia_dr2_cone_get_v2 as km4
    import mkpy3_vizier_vsx_cone_get_v2 as km5

    # assert(tpf is not None)
    if (tpf is None):
        tpf = lk.search_targetpixelfile(
            target='kepler-138b', mission='kepler', quarter=10).download(
                quality_bitmask=0)
        # ^--- exoplanet Kelper-138b is "KIC 7603200"
    # pass:if

    title_ = title
    rotate_deg = ast.literal_eval(rotate_deg_str)
    figsize = ast.literal_eval(figsize_str)
    colors = ast.literal_eval(colors_str)
    lws = ast.literal_eval(lws_str)
    zorders = ast.literal_eval(zorders_str)
    marker_kwargs = ast.literal_eval(marker_kwargs_str)
    gaia_dr2_kwargs = ast.literal_eval(gaia_dr2_kwargs_str)
    vsx_kwargs = ast.literal_eval(vsx_kwargs_str)

    assert(shrink >= 0.0)
    assert((percentile > 0.0) and (percentile <= 100.0))
    assert(isinstance(figsize, list))
    assert(len(figsize) == 2)
    assert(isinstance(colors, list))
    assert(len(colors) == 3)
    assert(isinstance(lws, list))
    assert(len(lws) == 3)
    assert(isinstance(zorders, list))
    assert(len(zorders) == 3)
    assert(isinstance(marker_kwargs, dict) or (marker_kwargs is None))
    assert(isinstance(gaia_dr2_kwargs, dict) or (gaia_dr2_kwargs is None))
    assert(isinstance(vsx_kwargs, dict) or (vsx_kwargs is None))

    if (verbose):
        print(tpf, 'tpf')
        print(frame, '=frame')
        print(survey, '=survey')
        print(rotate_deg, '=rotate_deg')
        print('^--- ', type(rotate_deg), '=type(rotate_deg)')
        print(width_height_arcmin, '=width_height_arcmin')
        print(shrink, '=shrink')
        print(show_plot, '=show_plot')
        print(plot_file, '=plot_file')
        print(overwrite, '=overwrite')
        print(figsize, '=figsize')
        print(title_, '=title')
        print(percentile, '=percentile')
        print(cmap, '=cmap')
        print(colors, '=colors')
        print(lws, '=lws')
        print(zorders, '=zorders')
        print(marker_kwargs, '=marker_kwargs')
        print(gaia_dr2_kwargs, '=gaia_dr2_kwargs')
        print(vsx_kwargs, '=vsx_kwargs')
        print(sexagesimal, '=sexagesimal')
        print(verbose, '=verbose')
    # pass:if

    ra_deg = tpf.ra
    dec_deg = tpf.dec
    if (verbose):
        print()
        print(ra_deg, '=ra_deg')
        print(dec_deg, '=dec_deg')
    # pass:if

    # get survey image data
    # survey = '2MASS-J'  # 'DSS2 Red' # hard-wired options
    survey_hdu, survey_hdr, survey_data, survey_wcs, survey_cframe = \
        km1.mkpy3_finder_chart_survey_fits_image_get_v1(
          ra_deg, dec_deg,
          radius_arcmin=width_height_arcmin, survey=survey, verbose=verbose)

    skip = (rotate_deg is None) or (rotate_deg == 0.0)
    if (not skip):
        print('\n***** ROTATE IMAGE *****:')
        if (verbose):
            print(rotate_deg, '=rotate_deg')
            print('^---', type(rotate_deg))
        # pass:if
        is_number = isinstance(rotate_deg, (float, int))
        if (is_number):
            if (verbose):
                print('rotate_deg is a number!')
            # pass:if
        else:
            is_str = isinstance(rotate_deg, str)
            if (is_str):
                if (verbose):
                    print('rotate_deg is a string')
                # pass:if
            else:
                print(type(rotate_deg), '=type(rotate_deg')
                print('^--- can not process this type of object')
                sys.exit(1)
            # pass:if
            if (rotate_deg == 'tpf'):
                if (verbose):
                    print(rotate_deg, ' is an acceptable string')
                # pass:if
            else:
                print('^--- ***ERROR*** NOT AN ACCEPTABLE STRING!')
                sys.exit(1)
            # pass:if
        # pass:if
        if (rotate_deg == 'tpf'):
            #
            # determine if X axis needs to be flipped =========================
            ny = tpf.flux.shape[1]
            nx = tpf.flux.shape[2]
            cx = nx / 2.0  # center X
            cy = ny / 2.0  # center Y
            #
            # HACK: BEGIN =====================================================
            #
            # determine if X axis needs to be flipped:
            # does RA (EAST arm) increase to the right?
            pixcrd0 = np.array([[cx, cy]], dtype=np.float_)  # center of compass rose
            # ^--- pixcrd0 must be a numpy 2-d array
            world = tpf.wcs.wcs_pix2world(pixcrd0, 0)  # pixels --> right ascension and declination
            world[0][0] += 1.0  # increase RA by 1.0 deg [*not* scaled by cos(DEC)]
            pixcrd1 = tpf.wcs.wcs_world2pix(world, 0)  # right ascension and declination --> pixels:
            e_x0 = pixcrd0[0][0]
            e_x1 = pixcrd1[0][0]
            flip_x = False
            if (e_x1 > e_x0):
                flip_x = True
            # pass:if
            #
            # determine if Y axis needs to be flipped =========================
            # does DEC (NORTH arm) increase to the top?
            pixcrd0 = np.array([[cx, cy]], dtype=np.float_)  # center of compass rose
            # ^--- pixcrd0 must be a numpy 2-d array
            world = tpf.wcs.wcs_pix2world(pixcrd0, 0)  # pixels --> right ascension and declination
            world[0][1] += 1.0  # increase DEC by 1.0 deg
            pixcrd1 = tpf.wcs.wcs_world2pix(world, 0)  # right ascension and declination --> pixels:
            n_y0 = pixcrd0[0][1]
            n_y1 = pixcrd1[0][1]
            flip_y = False
            if (n_y1 < n_y0):
                flip_y = True
            # pass:if
            #
            flip = flip_x or flip_y  # HACK
            print(flip_x, '=flip_x')
            print(flip_y, '=flip_y')
            print(flip, '=flip')
            #
            # HACK: END =======================================================
            #
            # compute rotation angle
            delta_y = tpf.wcs.wcs.pc[1, 0]
            delta_x = tpf.wcs.wcs.pc[1, 1]
            angle_rad = np.arctan2(delta_y, delta_x)
            # HACK: BEGIN =====================================================
            if (flip):
                angle_rad *= -1.0
            # pass:if
            # HACK: END =======================================================
            angle_deg = np.rad2deg(angle_rad)
            rotate_deg = angle_deg
            is_number = True
        # pass:if
        print('%.3f =rotate_deg [deg] : rotation angle' % rotate_deg)
        if (is_number):
            # modify survey World Coordinate System (WCS)
            wcs = cp.deepcopy(survey_wcs)
            theta = np.deg2rad(rotate_deg)
            # create rotation matrix
            rotation_matrix = np.matrix(
              [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
            # adjust the survey WCS pc coefficients
            rotated_pc = np.dot(rotation_matrix, wcs.wcs.pc)
            # update the survey WCS pc coefficients
            wcs.wcs.pc = rotated_pc
            # rotate survey image data using updated WCS
            reprojected_survey_data, reprojected_footprint = \
              rp.reproject_interp(
                survey_hdu, wcs, shape_out=survey_data.shape)
            # aliases
            survey_data = reprojected_survey_data
            survey_wcs = wcs
        # pass:if
        print()
    # pass:if

    # create a matplotlib figure object
    plt.figure(figsize=figsize)

    # create a matplotlib axis object with right ascension and declination axes
    ax = plt.subplot(projection=survey_wcs)

    # show the survey image
    km2.mkpy3_finder_chart_image_show_v1(
      ax=ax, image_data=survey_data,
      percentile=percentile, cmap=cmap, verbose=verbose)

    # show the TPF overlay
    km3.mkpy3_finder_chart_tpf_overlay_v6(
      ax=ax, survey_wcs=survey_wcs,
      tpf=tpf, frame=frame, colors=colors, lws=lws, zorders=zorders,
      verbose=verbose)

    # add title
    if (title_ is None):
        hdr = tpf.hdu[0].header  # alias
        try:  # Kepler?
            quarter = hdr['quarter']
            tag3_ = ('Quarter %02d' % quarter)
            tag2_ = hdr['object']
            tag1_ = 'Kepler'
            title_ = tag2_ + ' : ' + tag1_ + ' : ' + tag3_
        except Exception:
            try:  # K2?
                campaign = hdr['campaign']
                tag3_ = ('Campaign %02d' % campaign)
                tag2_ = hdr['object']
                tag1_ = 'K2'
                title_ = tag2_ + ' : ' + tag1_ + ' : ' + tag3_
            except Exception:  # TESS!
                sector = hdr['sector']
                tag2_ = ('Sector %02d' % sector)
                tag1_ = 'TESS'
                title_ = tag1_ + ' : ' + tag2_
            # pass:try
        # pass:try
        assert(title_ is not None)
    # pass:if
    plt.suptitle(title_, size=25)

    # option: mark the target
    if (isinstance(marker_kwargs, dict)):
        ax.scatter(
          ra_deg * u.deg, dec_deg * u.deg,
          transform=ax.get_transform(survey_cframe),
          **marker_kwargs)

    # CATALOGS:BEGIN  =========================================================

    ra_deg = tpf.ra
    dec_deg = tpf.dec

    fudge = np.sqrt(2) / 2.0
    radius_arcsec = width_height_arcmin * 60.0 * fudge * shrink
    print()
    print('%.6f =radius_arcsec  (%.6f =shrink)' % (radius_arcsec, shrink))

    # ===== GAIA DR2 CATALOG ==================================================

    proceed = isinstance(gaia_dr2_kwargs, dict) and (shrink > 0.0)
    while (proceed):

        raj2000, dej2000, sep_arcsec, gaia_dr2_result = \
            km4.mkpy3_vizier_gaia_dr2_cone_get_v2(
              ra_deg=ra_deg, dec_deg=dec_deg,
              radius_arcsec=radius_arcsec, verbose=verbose)
        if (gaia_dr2_result is None):  # nothing found
            break

        if (type(gaia_dr2_kwargs) is dict):
            ax.scatter(
                raj2000, dej2000,
                transform=ax.get_transform(survey_cframe), **gaia_dr2_kwargs)

        # numpy vectors of useful columns
        xra = raj2000  # alias
        yde = dej2000  # alias
        gmag = np.array(gaia_dr2_result['Gmag'])
        src = np.array(gaia_dr2_result['Source'])
        pmra = np.array(gaia_dr2_result['pmRA'])
        pmde = np.array(gaia_dr2_result['pmDE'])
        plx = np.array(gaia_dr2_result['Plx'])
        pmra = np.array(gaia_dr2_result['pmRA'])
        pmde = np.array(gaia_dr2_result['pmDE'])

        print()
        print(print_gaia_dr2, '=print_gaia_dr2')
        if (print_gaia_dr2):
            print('^--- set this keyword argument to False to *not* print', end='')
        else:
            print('^--- set this keyword argument to True to print', end='')
        # pass:if
        print(' the GAIA DR2 catalog results.')
        if (print_gaia_dr2):
            print()
            print()
            print('# GAIA DR2 : Global Astrometric Interferometer for Astrophy'
                  'sics-- Data Release 2')
            print('# n GAIA2_Source             sep    RA_ICRS      DE_ICRS   '
                  '    pmRA     pmDE      Plx     Gmag')
            print('#                       [arcsec]    [deg]        [deg]     '
                  ' [mas/yr] [mas/yr]    [mas]    [mag]')
            for k in range(len(xra)):
                j = k  # idx[k]
                raj = xra[j]
                dej = yde[j]
                gmagj = gmag[j]
                sepj = sep_arcsec[j]
                kk = k + 1
                srcj = src[j]
                pmraj = pmra[j]
                pmdej = pmde[j]
                plxj = plx[j]
                print('%3d %d %8.3f %12.7f %12.7f %8.3f %8.3f %8.3f %8.3f' % (
                  kk, srcj, sepj, raj, dej, pmraj, pmdej, plxj, gmagj))
            #  pass:for

            if (sexagesimal):
                print()
                print('# GAIA DR2 : Global Astrometric Interferometer for Astr'
                      'ophysics -- Data Release 2')
                print('# n GAIA2_Source          RA_ICRS       DE_ICRS      RA'
                      '_ICRS         DE_ICRS')
                print('#                         [deg]         [deg]        [h'
                      'ms]           [dms]')
                for k in range(len(xra)):
                    j = k  # idx[k]
                    xraj = xra[j]
                    ydej = yde[j]
                    gmagj = gmag[j]
                    sc1 = SkyCoord(
                        ra=xraj, dec=ydej, frame='icrs', unit='degree')
                    ra_ = sc1.ra.to_string(u.hour)
                    dec_ = sc1.dec
                    sepj = sep_arcsec[j]
                    kk = k + 1
                    srcj = src[j]
                    print('%3d %d %12.7f %12.7f  %15s %15s' % (
                      kk, srcj, xraj, ydej, ra_, dec_))
                # pass:for
            # pass:if
        # pass:if
        break
    # pass:while

    # ===== VSX CATALOG =======================================================

    proceed = isinstance(vsx_kwargs, dict) and (shrink > 0)
    while (proceed):

        raj2000, dej2000, sep_arcsec, vsx_result = \
          km5.mkpy3_vizier_vsx_cone_get_v2(
            ra_deg=ra_deg, dec_deg=dec_deg,
            radius_arcsec=radius_arcsec, verbose=verbose)
        if (vsx_result is None):  # nothing found
            break

        if (type(vsx_kwargs) is dict):
            ax.scatter(
              raj2000, dej2000,
              transform=ax.get_transform(survey_cframe), **vsx_kwargs)

        # numpy vectors of useful columns
        name = np.array(vsx_result['Name'], dtype=np.str)
        mag_max = np.array(vsx_result['max'])
        mag_min = np.array(vsx_result['min'])
        period = np.array(vsx_result['Period'])
        vsx_type = np.array(vsx_result['Type'], dtype=np.str)

        print()
        print(print_vsx, '=print_vsx')
        if (print_vsx):
            print('^--- set this keyword argument to False to *not* print', end='')
        else:
            print('^--- set this keyword argument to True to print', end='')
        # pass:if
        print(' the VSX catalog results.')
        if (print_vsx):
            print()
            print()
            print('# VSX : AAVSO International Variable Star indeX')
            print('# n      sep    RAJ2000      DEJ2000       Period     '
                  'VSX_max   VSX_min  VSX_Name      VSX_Type')
            print('#   [arcsec]    [deg]        [deg]         [days]     '
                  '[mag]     [mag]')
            for j in range(raj2000.size):
                k = j + 1
                raj = raj2000[j]
                dej = dej2000[j]
                sepj = sep_arcsec[j]
                pj = period[j]
                mxj = mag_max[j]
                mnj = mag_min[j]
                namej = name[j]
                vsxtypej = vsx_type[j]
                print(
                    "%3d %8.3f %12.7f %12.7f %12.6f %9.3f %9.3f '%s' '%s'" %
                    (k, sepj, raj, dej, pj, mxj, mnj, namej, vsxtypej))
            # pass:for

            if (sexagesimal):
                print()
                print('# VSX : AAVSO International Variable Star indeX')
                print('# n   RAJ2000       DEJ2000        RAJ2000        DEJ20'
                      '00')
                print('#     [deg]         [deg]          [hms]          '
                      '[dms]')
                for j in range(raj2000.size):
                    k = j + 1
                    xraj = raj2000[j]
                    ydej = dej2000[j]
                    sc1 = SkyCoord(
                        ra=xraj, dec=ydej, frame='icrs', unit='degree')
                    ra_ = sc1.ra.to_string(u.hour)
                    dec_ = sc1.dec
                    print('%3d %12.7f %12.7f  %15s %15s' % (
                      k, xraj, ydej, ra_, dec_))
                # pass:for
            # pass:if
        # pass:if
        break
    # pass:while

    # =========================================================================
    # CATALOGS: END ===========================================================
    # =========================================================================

    print()
    print(ra_deg, '=ra_deg')
    print(dec_deg, '=dec_deg')
    print('%d =cadenceno' % (tpf.cadenceno[frame]))
    print(frame, '=frame')

    # =========================================================================

    # reset plotting area to show only the survey image range in X and Y
    nx = survey_data.shape[1]
    ny = survey_data.shape[0]
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)

    # adjust the plot margins
    plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)

    if (plot_file == ''):
        plot_file = None
    if (plot_file is not None):
        if (plot_file != 'mkpy3_plot.png'):
            check_file_exists(plot_file, overwrite)
        plt.savefig(plot_file, dpi=300)  # , bbox_inches = "tight")
        print('\n%s <--- plot_file written\n' % (plot_file))
    # pass:if

    if (show_plot):
        plt.ioff()
        plt.show()
        ax = None
    # pass:if

    return ax
# pass:def

###############################################################################


if (__name__ == '__main__'):
    import matplotlib.pyplot as plt
    import mkpy3_plot_add_compass_rose_v4 as km
    #
    rotate = False
    show_compass = True
    show_plot = False
    interactive = True
    #
    if (not rotate):
        ax = mkpy3_tpf_overlay_v6(shrink=0.6, rotate_deg_str='0.0',
          show_plot=show_plot, plot_file='',
          print_gaia_dr2=False, sexagesimal=True)
    else:
        ax = mkpy3_tpf_overlay_v6(shrink=0.6, rotate_deg_str="'tpf'",
          show_plot=show_plot, plot_file='',
          print_gaia_dr2=False, sexagesimal=True)
    # pass:if
    ax.grid(True, color='palegreen', lw=2, zorder=1)
    if (show_compass):
        km.mkpy3_plot_add_compass_rose_v4(
          ax=ax, north_arm_arcsec=12)
    # pass:if
    plot_file = 'mkpy3_plot.png'
    plt.savefig(plot_file, bbox_inches="tight")
    print(plot_file, ' <--- new PNG file written')
    if (interactive):  # NOTE: does not work in a Jupyter notebook
        plt.ioff()
        plt.show()
    # pass:if
    plt.close()
# pass:if

# EOF
