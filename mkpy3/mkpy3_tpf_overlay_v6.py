#!/usr/bin/env python3

# file://mkpy3_tpf_overlay_v6.py

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute

# =============================================================================


def mkpy3__wcs_compass_check_v2(wcs, verbose=False):
    """
Utility function.
    """
    import numpy as np
    import copy as cp
    import inspect

    assert (wcs is not None)

    func_ = inspect.stack()[0][3]  # function name

    pc11 = wcs.wcs.pc[0][0]
    pc12 = wcs.wcs.pc[0][1]
    pc21 = wcs.wcs.pc[1][0]
    pc22 = wcs.wcs.pc[1][1]
    cdelt1 = wcs.wcs.cdelt[0]
    cdelt2 = wcs.wcs.cdelt[1]
    cd11 = cdelt1 * pc11
    cd12 = cdelt1 * pc12
    cd21 = cdelt2 * pc21
    cd22 = cdelt2 * pc22
    mirrored = ((cd11 * cd22) - (cd12 * cd21)) < 0.0
    positionAngle_deg = np.rad2deg(np.arctan2(cd12, cd11))
    if (mirrored):
        if (positionAngle_deg >= 0):
            positionAngle_deg += (-180.0)
        else:
            positionAngle_deg += (+180.0)
    # pass:if
    if (verbose):
        print()
        print('%s(): BEGIN ==================================================' %
          (func_))
        print()
        print(wcs)
        print('^--- wcs')
        print()
        print(positionAngle_deg, '=positionAngle_deg')
        print(mirrored, '=mirrored')
    # pass:if
    #
    cx = wcs.wcs.crpix[0] - 1  # CRPIX1 is one-offset
    cy = wcs.wcs.crpix[1] - 1  # CRPIX2 is one-offset
    #
    # does DEC (NORTH arm) increase to the TOP (increasing Y) ?
    pixcrd0 = np.array([[cx, cy]], dtype=np.float_)  # center of compass rose
    # ^--- pixcrd0 must be a numpy 2-d array
    world = wcs.wcs_pix2world(pixcrd0, 0)  # pixels --> RA and DEC
    world0 = cp.deepcopy(world)
    north_arm_deg = 0.1  # 6 arcmin
    world[0][1] += north_arm_deg  # increase DEC
    world1 = cp.deepcopy(world)
    pixcrd1 = wcs.wcs_world2pix(world, 0)  # RA and DEC --> pixels
    n_x0 = pixcrd0[0][0]
    n_y0 = pixcrd0[0][1]
    n_x1 = pixcrd1[0][0]
    n_y1 = pixcrd1[0][1]
    negate = -1.0
    n_dt = n_x1 - n_x0  # top
    n_db = n_y1 - n_y0  # bottom
    n_pa_deg = negate * np.rad2deg(np.arctan2(n_dt, n_db))
    #
    north_top_half = (n_y1 > n_y0)  # NORTH ARM top is ABOVE of center?
    if (verbose):
        print()
        print('***NORTH***:')
        print(cx, cy, '=cx,cy')
        print(wcs.wcs.crpix, '=wcs.wcs.crpix')
        print(pixcrd0, '=pixcrd0')
        print(pixcrd1, '=pixcrd1')
        print(n_dt, '=n_dt =(n_x1-n_x0)  [top]')
        print(n_db, '=n_db =(n_y1-n_y0)  [bottom]')
        print(world0, '=world0 : center (RA,DEC) [deg]')
        print(world1, '=world1 : North arm tip (RA,DEC) [deg]')
        # print(n_dv, n_dh, '=n_dv, n_dh')  # new
        print(n_pa_deg, '=n_pa_deg')  # new
        print(north_top_half, '=north_top_half')
    # pass:if
    #
    # does RA (EAST arm) increase to the LEFT (*decreasing* X) ?
    pixcrd0 = np.array([[cx, cy]], dtype=np.float_)  # center of compass rose
    # ^--- pixcrd0 must be a numpy 2-d array
    world = wcs.wcs_pix2world(pixcrd0, 0)  # pixels --> RA and DEC
    world0 = cp.deepcopy(world)
    declination = world[0][1]
    east_arm_deg = north_arm_deg
    world[0][0] += east_arm_deg / np.cos(np.deg2rad(declination))  # increase RA
    world1 = cp.deepcopy(world)
    pixcrd1 = wcs.wcs_world2pix(world, 0)  # RA and DEC --> pixels
    e_x0 = pixcrd0[0][0]
    e_x1 = pixcrd1[0][0]
    east_left_half = (e_x1 < e_x0)  # EAST ARM tip LEFT is LEFT of center?
    if (verbose):
        print('\n***EAST***:')
        print(cx, cy, '=cx,cy')
        print(wcs.wcs.crpix, '=wcs.wcs.crpix')
        print(pixcrd0, '=pixcrd0')
        print(pixcrd1, '=pixcrd1')
        print(world0, '=world0 : center (RA,DEC) [deg]')
        print(world1, '=world1 : North arm tip (RA,DEC) [deg]')
        print(east_left_half, '=east_left_half')
        print()
    # pass:if

    # sanity check
    delta_deg = (positionAngle_deg - n_pa_deg)
    if (np.fabs(delta_deg) >= 0.1):
        print('\n\n')
        print('***WARNING*** BEGIN ==========================================')
        print('***INFO***  %s()' % (func_))
        print('***INFO***:\n', wcs, '\n^-- wcs\n')
        print('***INFO***', positionAngle_deg, '=positionAngle_deg')
        print('***INFO***', n_pa_deg, '=n_pa_deg  <--------------------------')
        print('***INFO***', delta_deg,
          '= delta_deg =(positionAngle_deg - n_pa_deg)')
        print('***INFO***', (np.fabs(delta_deg) < 0.1),
          '=(np.fabs(delta_deg) < 0.1)  [*** WARNING *** : should be True  8=X]')
        print('***INFO***', n_dt, '=n_dt  [px]  (delta X axis)')
        print('***INFO***', n_db, '=n_db  [px]  (delta Y axis)')
        print('***WARNING*** END ============================================')
        print()
    # pass:if

    # that's all folks!
    if (verbose):
        print('%s(): END ====================================================' %
          (func_))
        print()
    # pass:if

    chatty = False
    if (chatty):
        print()
        print(func_, positionAngle_deg, '=positionAngle_deg')
        print(func_, n_pa_deg, '=n_pa_deg')
        print(func_, mirrored, '=mirrored')
        print(func_, north_top_half, '=north_top_half')
        print(func_, east_left_half, '=east_left_half')
    # pass:if

    return positionAngle_deg, n_pa_deg, mirrored, north_top_half, east_left_half
# pass:def

# =============================================================================


def mkpy3_tpf_overlay_v6(
  tpf=None,
  frame=0,
  survey='2MASS-J',  # '2MASS-J' or 'DSS2 Red'
  width_height_arcmin=2.0,
  rotationAngle_deg=None,
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
    '''
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
rotationAngle_deg : (None) or (float) or ('tpf') (optional)
    Angle in degrees to rotate the sky survey image
    [default for Kepler & K2 missions: None]
    [default for TESS mission: 'tpf']
    [example values: None, 12.345, 'tpf']
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
    A string of a 2-element list of figure width and height [Matplotlib].
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
    Returns a matplotlib axes object *if* show_plot is False *else* None .

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute
    '''
    import matplotlib.pyplot as plt
    import numpy as np
    import copy as cp
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import ast
    import sys
    import inspect
    #
    import reproject as rp  # pip install reproject
    #
    import lightkurve as lk
    lk.log.setLevel('INFO')
    #
    import mkpy3

    func_ = inspect.stack()[0][3]  # function name

    # assert(tpf is not None)
    if (tpf is None):
        tpf = lk.search_targetpixelfile(
            target='kepler-138b', mission='kepler', quarter=10).download(
                quality_bitmask=0)
        # ^--- exoplanet Kelper-138b is "KIC 7603200"
    # pass:if

    if ((rotationAngle_deg is None) and (tpf.mission.lower() == 'tess')):
        rotationAngle_deg = 'tpf'
    # pass:if

    title_ = title
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
        print(func_, '=func_')
        print(tpf, '=tpf')
        print(frame, '=frame')
        print(survey, '=survey')
        print(rotationAngle_deg, '=rotationAngle_deg')
        print('^---', type(rotationAngle_deg), '=type(rotationAngle_deg)  ***INFO***')
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
        print(print_gaia_dr2, '=print_gaia_dr2')
        print(gaia_dr2_kwargs, '=gaia_dr2_kwargs')
        print(print_vsx, '=print_vsx')
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

    tpf_positionAngle_deg,\
      tpf_n_pa_deg,\
      tpf_mirrored,\
      tpf_north_top_half,\
      tpf_east_left_half \
      = mkpy3__wcs_compass_check_v2(wcs=tpf.wcs, verbose=verbose)
    if (verbose):
        print(tpf_positionAngle_deg, '=tpf_positionAngle_deg')
        print(tpf_n_pa_deg, '=tpf_n_pa_deg')
        print(tpf_mirrored, '=tpf_mirrored')
        print(tpf_north_top_half, '=tpf_north_top_half')
        print(tpf_east_left_half, '=tpf_east_left_half')
    # pass:if

    # get survey image data
    # survey = '2MASS-J'  # 'DSS2 Red' # hard-wired options
    survey_hdu, survey_hdr, survey_data, survey_wcs, survey_cframe = \
        mkpy3.mkpy3_finder_chart_survey_fits_image_get_v1(
          ra_deg, dec_deg,
          radius_arcmin=width_height_arcmin, survey=survey, verbose=verbose)

    survey_positionAngle_deg,\
      survey_n_pa_deg,\
      survey_mirrored,\
      survey_north_top_half,\
      survey_east_left_half \
      = mkpy3__wcs_compass_check_v2(wcs=survey_wcs, verbose=verbose)
    if (verbose):
        print(survey_positionAngle_deg, '=survey_positionAngle_deg')
        print(survey_n_pa_deg, '=survey_n_pa_deg')
        print(survey_mirrored, '=survey_mirrored')
        print(survey_north_top_half, '=survey_north_top_half')
        print(survey_east_left_half, '=survey_east_left_half')
    # pass:if

    survey_rotate_deg = None
    rotate_deg = cp.deepcopy(rotationAngle_deg)
    skip = (rotate_deg is None) or (rotate_deg == 0.0)
    if (not skip):
        print('\n***** ROTATE IMAGE *****:')
        if (verbose):
            print(rotationAngle_deg, '=rotationAngle_deg')
            print('^---', type(rotationAngle_deg), '=type(rotationAngle_deg)')
        # pass:if
        is_number = isinstance(rotationAngle_deg, (float, int))
        if (is_number):
            if (verbose):
                print('rotationAngle_deg is a number!')
            # pass:if
        else:
            is_str = isinstance(rotationAngle_deg, str)
            if (is_str):
                if (verbose):
                    print('rotationAngle_deg is a string')
                # pass:if
            else:
                print(type(rotationAngle_deg), '=type(rotationAngle_deg')
                print('^--- can not process this type of object')
                sys.exit(1)
            # pass:if
            if (rotationAngle_deg != 'tpf'):  # a 3-char str with a value of tpf
                print(rotationAngle_deg, '=rotationAngle_deg')
                print("^--- ***ERROR*** BAD STRING VALUE!  ONLY 'tpf' IS ALLOWED!")
                sys.exit(1)
            # pass:if
        # pass:if
        if (rotationAngle_deg == 'tpf'):  # a 3-char str with a value of tpf
            tpf_pa_deg = tpf_n_pa_deg  # less reliable: tpf_positionAngle_deg
            if (not tpf_east_left_half):
                tpf_pa_deg *= (-1.0)
            # pass:if
            tpf_pa_fabs_deg = np.fabs(tpf_pa_deg)
            tpf_pa_sign = np.sign(tpf_pa_deg)
            if (verbose):
                print('[!]', tpf_pa_deg, '=tpf_pa_deg')
                print('[!]', tpf_pa_fabs_deg, '=tpf_pa_fabs_deg')
                print('[!]', tpf_pa_sign, '=tpf_pa_sign')
            # pass:if
            reverse = (-1.0)
            if (tpf_pa_fabs_deg <= 90.0):
                rotate_deg = reverse * tpf_pa_deg
                if (verbose):
                    print('[!]', rotate_deg, '= rotate_deg = reverse * tpf_pa_deg')
                # pass:if
            else:
                rotate_deg = reverse * tpf_pa_sign * (180.0 - tpf_pa_fabs_deg)
                if (verbose):
                    print('[!]', rotate_deg,
                      '= reverse * tpf_pa_sign * (180.0 - tpf_pa_fabs_deg)')
                # pass:if
            is_number = True
        # pass:if
        survey_rotate_deg = rotate_deg
        if (verbose):
            print('%.3f =survey_rotate_deg [deg] : rotation angle '
              '<------------' % survey_rotate_deg)
        # pass:if
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
            #
            survey_positionAngle_deg,\
              survey_n_pa_deg,\
              survey_mirrored,\
              survey_north_top_half,\
              survey_east_left_half \
              = mkpy3__wcs_compass_check_v2(wcs=survey_wcs, verbose=verbose)
            print(survey_rotate_deg, '=survey_rotate_deg')
            if (verbose):
                print(survey_positionAngle_deg, '=survey_positionAngle_deg [revised]')
                print(survey_n_pa_deg, '=survey_n_pa_deg [revised]')
                print(survey_mirrored, '=survey_mirrored [revised]')
                print(survey_north_top_half, '=survey_north_top_half [revised]')
                print(survey_east_left_half, '=survey_east_left_half [revised]')
            # pass:if
        # pass:if
    # pass:if

    # create a matplotlib figure object
    plt.figure(figsize=figsize)

    # create a matplotlib axis object with right ascension and declination axes
    ax = plt.subplot(projection=survey_wcs)

    # HACK: BEGIN : *NEW* attributes of the axis
    ax.tpf_positionAngle_deg = tpf_positionAngle_deg
    ax.tpf_n_pa_deg = tpf_n_pa_deg
    ax.tpf_mirrored = tpf_mirrored
    ax.tpf_north_top_half = tpf_north_top_half
    ax.tpf_east_left_half = tpf_east_left_half
    ax.survey_rotationAngle_deg = rotationAngle_deg
    ax.survey_rotate_deg = survey_rotate_deg
    ax.survey_positionAngle_deg = survey_positionAngle_deg
    ax.survey_n_pa_deg = survey_n_pa_deg
    ax.survey_mirrored = survey_mirrored
    ax.survey_north_top_half = survey_north_top_half
    ax.survey_east_left_half = survey_east_left_half
    # HACK: END

    # show the survey image
    mkpy3.mkpy3_finder_chart_image_show_v1(
      ax=ax, image_data=survey_data,
      percentile=percentile, cmap=cmap, verbose=verbose)

    # show the TPF overlay
    mkpy3.mkpy3_finder_chart_tpf_overlay_v6(
      ax=ax, survey_wcs=survey_wcs,
      tpf=tpf, frame=frame, colors=colors, lws=lws, zorders=zorders,
      verbose=verbose)

    # add suptitle
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
    print(proceed, '=proceed [GAIA2 catalog]')
    while (proceed):

        raj2000, dej2000, sep_arcsec, gaia_dr2_result = \
            mkpy3.mkpy3_vizier_gaia_dr2_cone_get_v2(
              ra_deg=ra_deg, dec_deg=dec_deg,
              radius_arcsec=radius_arcsec, verbose=verbose)
        if (gaia_dr2_result is None):  # nothing found
            print(gaia_dr2_result,
              '=gaia_dr2_result [***WARNING*** NOTHING FOUND]')
            break
        # pass:if
        if (type(gaia_dr2_kwargs) is dict):
            ax.scatter(
                raj2000, dej2000,
                transform=ax.get_transform(survey_cframe), **gaia_dr2_kwargs)
        # pass:if

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
          mkpy3.mkpy3_vizier_vsx_cone_get_v2(
            ra_deg=ra_deg, dec_deg=dec_deg,
            radius_arcsec=radius_arcsec, verbose=verbose)
        if (vsx_result is None):  # nothing found
            print(vsx_result, '=vsx_result [***WARNING*** NOTHING FOUND]')
            break
        # pass:if
        if (type(vsx_kwargs) is dict):
            ax.scatter(
              raj2000, dej2000,
              transform=ax.get_transform(survey_cframe), **vsx_kwargs)
        # pass:if

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
    print()
    print('%d =cadenceno' % (tpf.cadenceno[frame]))
    print(frame, '=frame')

    # reset plotting area to show only the survey image range in X and Y
    nx = survey_data.shape[1]
    ny = survey_data.shape[0]
    ax.set_xlim(0, nx)
    ax.set_ylim(0, ny)

    # adjust the plot margins
    plt.subplots_adjust(left=0.2, right=0.9, top=0.9, bottom=0.2)

    # match orienttion of tpf.plot() graph ====================================
    ax.xaxis_inverted = False
    ax.yaxis_inverted = False
    if (rotationAngle_deg == 'tpf'):  # a 3-char str with a value of tpf
        if (ax.survey_east_left_half != ax.tpf_east_left_half):
            ax.invert_xaxis()
            ax.xaxis_inverted = True
        # pass:if
        if (ax.survey_north_top_half != ax.tpf_north_top_half):
            ax.invert_yaxis()
            ax.yaxis_inverted = True
        # pass:if
    # pass:if
    # =========================================================================

    if (plot_file == ''):
        plot_file = None
    if (plot_file is not None):
        if (plot_file != 'mkpy3_plot.png'):
            mkpy3.mkpy3_util_check_file_exists(plot_file, overwrite)
        plt.savefig(plot_file, dpi=300)  # , bbox_inches = "tight")
        print('\n%s <--- plot_file written\n' % (plot_file))
    # pass:if

    if (show_plot):
        plt.ioff()
        plt.show()
        ax = None
    # pass:if

    print()
    print('DONE:', func_)
    print()

    return ax
# pass:def

# =============================================================================


def xmkpy3_tpf_overlay_v6():
    """
Unit test
    """
    import matplotlib.pyplot as plt
    import os
    import ntpath
    import datetime
    import lightkurve as lk
    from astropy.visualization import ImageNormalize, PercentileInterval,\
      SqrtStretch
    #
    import mkpy3

    cmap = 'gray_r'
    verbose = False

    mission = None
    target = None
    quarter = None
    campaign = None
    sector = None

    print('\n\nDATA ==========================================================')

    obj = 6  # <--- USER CUSTOMIZE

    if (obj == 1):
        target = 'CD Ind'
        mission = 'TESS'
        sector = 1
        width_height_arcmin = 6  # use with TESS  <--- USER CUSTOMIZE
        north_arm_arcsec = 42  # USER CUSTOMIZE
        frame = 1  # USER CUSTOMIZE
        percentile = 99.9  # USER CUSTOMIZE
        title_ = target + ' : TESS : Sector ' + str(sector)  # USER CUSTOMIZE
        print('DOWNLOAD START:', str(datetime.datetime.now()))
        search_results = lk.search_tesscut(target=target, sector=sector)[0]
        tpf = search_results[0].download(cutout_size=(11, 11), quality_bitmask=0)
        print('DOWNLOAD  STOP:', str(datetime.datetime.now()))
    # pass:if

    if (obj == 2):
        target = 'XZ Cyg'
        mission = 'TESS'
        sector = 14
        width_height_arcmin = 6  # use with TESS  <--- USER CUSTOMIZE
        north_arm_arcsec = 42  # USER CUSTOMIZE
        frame = 1  # USER CUSTOMIZE
        percentile = 99.9  # USER CUSTOMIZE
        title_ = target + ' : TESS : Sector ' + str(sector)  # USER CUSTOMIZE
        print('DOWNLOAD START:', str(datetime.datetime.now()))
        search_results = lk.search_tesscut(target=target, sector=sector)[0]
        tpf = search_results[0].download(cutout_size=(11, 11), quality_bitmask=0)
        print('DOWNLOAD  STOP:', str(datetime.datetime.now()))
    # pass:if

    if (obj == 3):
        target = 'V1460 Her'
        mission = 'TESS'
        sector = 24
        width_height_arcmin = 6  # use with TESS  <--- USER CUSTOMIZE
        north_arm_arcsec = 42  # USER CUSTOMIZE
        frame = 1  # USER CUSTOMIZE
        percentile = 99.9  # USER CUSTOMIZE
        title_ = target + ' : TESS : Sector ' + str(sector)  # USER CUSTOMIZE
        print('DOWNLOAD START:', str(datetime.datetime.now()))
        search_results = lk.search_tesscut(target=target, sector=sector)[0]
        tpf = search_results[0].download(cutout_size=(11, 11), quality_bitmask=0)
        print('DOWNLOAD  STOP:', str(datetime.datetime.now()))
    # pass:if

    if (obj == 4):
        target = 'Kepler-138b'
        mission = 'Kepler'
        quarter = 10
        width_height_arcmin = 1.8  # use with Kepler/K2  <--- USER CUSTOMIZE
        north_arm_arcsec = 6  # USER CUSTOMIZE
        frame = 1  # USER CUSTOMIZE
        percentile = 99.0  # USER CUSTOMIZE
        title_ = target + ' : Kepler : Quarter ' + str(quarter)  # USER CUSTOMIZE
        print('DOWNLOAD START:', str(datetime.datetime.now()))
        tpf = lk.search_targetpixelfile(
          target=target, mission=mission, quarter=quarter).download(
          quality_bitmask=0)
        # ^--- exoplanet Kelper-138b is "KIC 7603200"
        print('DOWNLOAD  STOP:', str(datetime.datetime.now()))
    # pass:if

    if (obj == 5):
        target = 'K2-34b'
        mission = 'k2'
        campaign = 18
        width_height_arcmin = 1.8  # use with Kepler/K2  <--- USER CUSTOMIZE
        north_arm_arcsec = 6  # USER CUSTOMIZE
        frame = 1  # USER CUSTOMIZE
        percentile = 99.0  # USER CUSTOMIZE
        title_ = target + ' : K2 : Campaign ' + str(campaign)  # USER CUSTOMIZE
        print('DOWNLOAD START:', str(datetime.datetime.now()))
        tpf = lk.search_targetpixelfile(
          target=target, mission=mission, campaign=campaign).download(
          quality_bitmask=0)
        # ^--- exoplanet K2-34b is "EPIC 212110888"
        print('DOWNLOAD  STOP:', str(datetime.datetime.now()))
    # pass:if

    if (obj == 6):
        target = 'CD Ind'
        sector = 1
        mission = 'TESS'
        frame = 1
        width_height_arcmin = 6  # USER CUSTOMIZE
        north_arm_arcsec = 42    # USER CUSTOMIZE
        percentile = 99.5        # USER CUSTOMIZE
        title_ = target + ' : TESS : Sector ' + str(sector)  # USER CUSTOMIZE
        #
        radius = 120  # arcsec
        print('DOWNLOAD START:', str(datetime.datetime.now()))
        search_results = lk.search_targetpixelfile(target,
          radius=radius, mission=mission, sector=sector)
        tpf = search_results[0].download(quality_bitmask=0)
        print('DOWNLOAD  STOP:', str(datetime.datetime.now()))
        print(search_results)
        print(tpf)
        print(tpf.wcs)
        print(':-)')
    # pass:if

    print()
    print(mission, '=mission')
    print(target, '=target')
    if (mission == 'kepler'):
        print(quarter, '=quarter')
    # pass:if
    if (mission == 'k2'):
        print(campaign, '=campaign')
    # pass:if
    if (mission == 'tess'):
        print(sector, '=sector')
    # pass:if

    tpf_dir = os.path.dirname(tpf.path)
    tpf_file = ntpath.basename(tpf.path)
    cwd = os.getcwd()
    print('\n         CWD:', cwd)
    print(' TPF dirname:', tpf_dir)
    print('TPF filename:', tpf_file)

    print()
    print(tpf)
    print('^--- tpf')
    if (verbose):
        print()
        print(tpf.wcs)
        print('^--- tpf.wcs')
    # pass:if

    print('\n\nPLOT#1 =======================================================')
    ax = tpf.plot(frame=frame, cmap=cmap)
    ax.set_title(title_)
    mkpy3.mkpy3_plot_add_compass_rose_v5(ax=ax, north_arm_arcsec=north_arm_arcsec,
      wcs=tpf.wcs, verbose=verbose)
    # mark the target with a yellow circle:
    target_wx0, target_wy0 = tpf.wcs.wcs_world2pix(tpf.ra, tpf.dec, 0)
    marker_kwargs = \
      {'edgecolor': 'yellow', 's': 600, 'facecolor': 'None', 'lw': 3, 'zorder': 10}
    ax.scatter(target_wx0 + tpf.column, target_wy0 + tpf.row, **marker_kwargs)
    oplot1 = 'mkpy3_plot1.png'
    plt.savefig(oplot1, bbox_inches="tight")
    print('\n', oplot1, ' <--- new PNG file written')
    plt.close()

    print('\n\nPLOT#2 =======================================================')
    fig = plt.figure(figsize=(7, 7))
    ax = plt.subplot(projection=tpf.wcs)
    image_data = tpf.hdu[1].data['flux'][frame]
    norm = ImageNormalize(image_data, interval=PercentileInterval(percentile),
      stretch=SqrtStretch())
    ax.imshow(image_data, norm=norm, cmap=cmap)
    ax.tick_params(axis='x', labelsize=16, length=5, width=2, labeltop=True,
      labelbottom=True)
    ax.tick_params(axis='y', labelsize=16, length=5, width=2, labelright=True,
      labelleft=True)
    ax.coords[0].set_major_formatter('d.dd')
    ax.coords[1].set_major_formatter('d.dd')
    ax.set_xlabel('Right Ascension (J2000)', size=24)
    ax.set_ylabel('Declination (J2000)', size=24)
    fig.suptitle(title_, size=24)
    ax.grid(True, color='palegreen', lw=2, zorder=1)
    mkpy3.mkpy3_plot_add_compass_rose_v5(ax=ax, north_arm_arcsec=north_arm_arcsec,
      wcs=tpf.wcs, verbose=verbose)
    marker_kwargs =\
      {'edgecolor': 'yellow', 's': 600, 'facecolor': 'None', 'lw': 3, 'zorder': 10}
    ax.scatter(tpf.ra, tpf.dec, transform=ax.get_transform('icrs'), **marker_kwargs)
    oplot2 = 'mkpy3_plot2.png'
    plt.savefig(oplot2, dpi=150, bbox_inches="tight")
    print('\n', oplot2, ' <--- plot file written')
    plt.close()

    # =========================================================================

    shrink = 0.75
    print_gaia_dr2 = False
    print_vsx = False

    print('\n\nPLOT#4 =======================================================')
    rotationAngle_deg = 0.0  # no rotation
    ax = mkpy3_tpf_overlay_v6(tpf=tpf, rotationAngle_deg=rotationAngle_deg,
      width_height_arcmin=width_height_arcmin, percentile=percentile,
      shrink=shrink, show_plot=False, plot_file='', title=title_,
      print_gaia_dr2=print_gaia_dr2, print_vsx=print_vsx, verbose=verbose)
    ax.coords[0].set_major_formatter('d.dd')
    ax.coords[1].set_major_formatter('d.dd')
    ax.tick_params(axis='x', labelsize=16, length=5, width=2, labeltop=True,
      labelbottom=True)
    ax.tick_params(axis='y', labelsize=16, length=5, width=2, labelright=True,
      labelleft=True)
    ax.grid(True, color='palegreen', lw=2, zorder=1)
    mkpy3.mkpy3_plot_add_compass_rose_v5(ax=ax, north_arm_arcsec=2 * north_arm_arcsec)
    oplot4 = 'mkpy3_plot4.png'
    plt.savefig(oplot4, bbox_inches="tight")
    print('\n', oplot4, ' <--- new PNG file written')
    plt.close()

    print('\n\nPLOT#3 =======================================================')
    #
    # compute rotation based on the WCS of the TPF:
    rotationAngle_deg = 'tpf'
    #
    ax = mkpy3_tpf_overlay_v6(tpf=tpf, rotationAngle_deg=rotationAngle_deg,
      width_height_arcmin=width_height_arcmin,
      shrink=shrink, show_plot=False, plot_file='', title=title_,
      percentile=percentile,
      print_gaia_dr2=print_gaia_dr2, print_vsx=print_vsx, verbose=verbose)
    ax.coords[0].set_major_formatter('d.dd')
    ax.coords[1].set_major_formatter('d.dd')
    ax.tick_params(axis='x', labelsize=16, length=5, width=2, labeltop=True,
      labelbottom=True)
    ax.tick_params(axis='y', labelsize=16, length=5, width=2, labelright=True,
      labelleft=True)
    ax.grid(True, color='palegreen', lw=2, zorder=1)
    mkpy3.mkpy3_plot_add_compass_rose_v5(
      ax=ax, north_arm_arcsec=2 * north_arm_arcsec)
    oplot3 = 'mkpy3_plot3.png'
    plt.savefig(oplot3, bbox_inches="tight")
    print('\n', oplot3, ' <--- new PNG file written')
    plt.close()

    print()
    print('==================================================================')
    print()
    print('[*]', ax.tpf_positionAngle_deg, '=ax.tpf_positionAngle_deg')
    print('[*]', ax.tpf_n_pa_deg, '=ax.tpf_n_pa_deg')
    print('[*]', ax.tpf_mirrored, '=ax.tpf_mirrored')
    print('[*]', ax.tpf_north_top_half, '=ax.tpf_north_top_half')
    print('[*]', ax.tpf_east_left_half, '=ax.tpf_east_left_half')
    print('[*]')
    print('[*] <---', ax.survey_rotationAngle_deg, '=ax.survey_rotationAngle_deg')
    print('[*] --->', ax.survey_rotate_deg, '=ax.survey_rotate_deg')
    print('[*]')
    print('[*]', ax.survey_positionAngle_deg, '=ax.survey_positionAngle_deg')
    print('[*]', ax.survey_n_pa_deg, '=ax.survey_n_pa_deg')
    print('[*]', ax.survey_mirrored, '=ax.survey_mirrored')
    print('[*]', ax.survey_north_top_half, '=ax.survey_north_top_half')
    print('[*]', ax.survey_east_left_half, '=ax.survey_east_left_half')
    print('[*]')
    print('[*]', ax.xaxis_inverted, '=ax.xaxis_inverted')
    print('[*]', ax.yaxis_inverted, '=ax.yaxis_inverted')

    print()
    print('[-]', target, '=target')
    print('[-]', mission, '=mission')
    if (mission == 'TESS'):
        print('[-]', sector, '=sector')
    elif (mission == 'kepler'):
        print('[-]', quarter, '=quarter')
    else:
        print('[-]', campaign, '=campaign')
    # pass:if
    print('[-]', tpf.ra, '=tpf.ra [deg]')
    print('[-]', tpf.dec, '=tpf.dec [deg]')
    print('[-]', obj, '=obj')

    print()
    print()
    print('PLOTS DONE =======================================================')
    print()
    print('plot files written:\n')
    print(' ', oplot1, oplot2, oplot3, oplot4)

    print()
    print('FILE:', __file__)
    print('TIME:', str(datetime.datetime.now()))
    print('DONE:', __name__)
# pass:def

# =============================================================================


if (__name__ == '__main__'):
    xmkpy3_tpf_overlay_v6()
# pass:if

# EOF
