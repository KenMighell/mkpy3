#!/usr/bin/env python3

# file://mkpy3_plot_add_compass_rose_v5.py

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute

# =============================================================================


def mkpy3_plot_add_compass_rose_v5(
  ax=None,
  wcs=None,
  draw=None,
  cx=None,
  cy=None,
  north_arm_arcsec=None,
  edge_color=None,
  inside_color=None,
  edge_lw=None,
  inside_lw=None,
  verbose=None
):
    """
Function : mkpy3_plot_add_compass_rose_v5()

Purpose: Add a compass rose to a matplotlib axis.
    NOTE: Long arm points North (increasing declination)
    NOTE: Short arm points East (increasing right ascension)

Parameters
----------
ax : matplotlib Axes object
wcs : astropy FITS World Coordinate System (WCS) object (optional)
    [default: None --> ax.wcs]
    [NOTE: lightkurve TPF plots: use wcs=tpf.wcs]
cx : (float) (optional)
    pixel column number (X position) of the center of the compass rose
cy : (float) (optional)
    pixel row number (Y position) of the center of the compass rose
draw : (bool) (optional)
    draw the compass if True [default: True]
north_arm_arcsec : (float) (optional) [default: 6 arcsec]
    length of the North arrow arm of the compass rose in arcsec
    [N.B. A Kepler photometer pixel is 3.98 arcsec/pixel]
    [N.B. A TESS photometer pixel is 21.00 arcsec/pixel]
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
tpf = lk.search_targetpixelfile(
  target='Kepler-138b', mission='Kepler', quarter=10).download()
#         ^--- Exoplanet Kelper-138b is "KIC 7603200"
#
# Plot the 2nd frame of the TPF
ax = tpf.plot(frame=1)
#
# add a compass rose using the WCS from the TargetPixelFile
mkpy3_plot_add_compass_rose_v5(ax=ax, wcs=tpf.wcs, verbose=True)
#
plot_file = 'mkpy3_plot.png'
plt.savefig(plot_file, bbox_inches="tight")
plt.close()
print(plot_file, ' <--- new PNG file written')
#==========

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.coordinates as ac
    import astropy.units as u
    import inspect

    func_ = inspect.stack()[0][3]  # function name

    assert(ax is not None), '***ERROR***: Requires a matplotlib axes object'
    if (wcs is None):
        assert(hasattr(ax, 'wcs')),\
          '***ERROR***: Requires that ax.wcs exists (is an attribute)'
        wcs = ax.wcs  # alias
    # pass:if
    # parse the last line of the WCS summary to (maybe) get NAXIS1 and NAXIS2:
    naxis_list = repr(wcs).splitlines()[-1].split()
    nxw = int(naxis_list[2])  # =? NAXIS1 [***WARNING***: sometimes False]
    nyw = int(naxis_list[3])  # =? NAXIS2 [***WARNING***: sometimes False]
    #
    xmin, xmax = ax.get_xlim()  # X-axis min, X-axis max
    ymin, ymax = ax.get_ylim()  # Y-axis min, Y-axis max
    llx = np.int(np.rint(xmin + 0.5))  # should be equal to tpf.column
    lly = np.int(np.rint(ymin + 0.5))  # should be equal to tpf.row
    if (draw is None):
        draw = True
    # pass:if
    if (cx is None):
        xlim = ax.get_xlim()
        cx = ((xlim[0] + xlim[1]) / 2.0)  # center of the X axis
    # pass:if
    cxw = cx - llx  # window X coordinate
    if (cy is None):
        ylim = ax.get_ylim()
        cy = ((ylim[0] + ylim[1]) / 2.0)  # center of the Y axis
    # pass:if
    cyw = cy - lly  # window Y coordinate
    if (north_arm_arcsec is None):
        north_arm_arcsec = 6  # default for Kepler/K2 observations
        # north_arm_arcsec = 63  # default for TESS observations
    # pass:if
    if (edge_color is None):
        edge_color = 'blue'
    # pass:if
    if (inside_color is None):
        inside_color = 'yellow'
    # pass:if
    if (edge_lw is None):
        edge_lw = 7
    # pass:if
    if (inside_lw is None):
        inside_lw = 4
    # pass:if
    if (verbose is None):
        verbose = False
    # pass:if
    if (verbose):
        print()
        print(ax, '=ax')
        print(wcs)
        print('^--- =wcs')
        print(draw, '=draw')
        print(cx, '=cx')
        print(cy, '=cy')
        print(xmin, '=xmin ***INFO***')
        print(xmax, '=xmax ***INFO***')
        print(ymin, '=ymin ***INFO***')
        print(ymax, '=ymax ***INFO***')
        print(llx, '=llx ***INFO***')
        print(lly, '=lly ***INFO***')
        print(cxw, '=cxw ***INFO***')
        print(cyw, '=cyw ***INFO***')
        print(nxw, '=nxw ***INFO***')
        print(nyw, '=nyw ***INFO***')
        print(xmin, xmax, '=xmin,xmax ***INFO***')
        print(north_arm_arcsec, '=north_arm_arcsec')
        print(edge_color, '=edge_color')
        print(inside_color, '=inside_color')
        print(edge_lw, '=edge_lw')
        print(inside_lw, '=inside_lw')
        print(verbose, '=verbose')
        print()
        print('==============================================================')
    # pass:if

    cx0 = cxw  # alias: zero-offset window coordinates
    cy0 = cyw  # alias: zero-offset window coordinates

    north_arm_deg = north_arm_arcsec / 3600.0  # long compass arm point North
    east_arm_deg = north_arm_deg / 2.0  # short compass arm points East
    east_arm_arcsec = east_arm_deg * 3600.0
    if (verbose):
        print()
        print(north_arm_arcsec, '=north_arm_arcsec')
        print(east_arm_arcsec, '=east_arm_arcsec')
    # pass:if

    # NORTH arm of compass rose
    pixcrd0 = np.array([[cx0, cy0]], dtype=np.float_)
    # ^--- pixcrd0 must be a numpy 2-d array
    # pixels --> right ascension and declination:
    world = wcs.wcs_pix2world(pixcrd0, 0)
    world[0][1] += north_arm_deg
    # right ascension and declination --> pixels:
    pixcrd1 = wcs.wcs_world2pix(world, 0)
    n_x0 = pixcrd0[0][0] + llx  # window --> CCD
    n_y0 = pixcrd0[0][1] + lly  # window --> CCD
    n_x1 = pixcrd1[0][0] + llx  # window --> CCD
    n_y1 = pixcrd1[0][1] + lly  # window --> CCD
    negate = -1.0
    n_dt = n_x1 - n_x0  # top
    n_db = n_y1 - n_y0  # bottom
    n_pa_deg = negate * np.rad2deg(np.arctan2(n_dt, n_db))

    # sanity check
    if (verbose):
        world0 = wcs.wcs_pix2world(pixcrd0, 0)
        world1 = wcs.wcs_pix2world(pixcrd1, 0)
        sc0 = ac.SkyCoord(world0, unit=u.deg)
        sc1 = ac.SkyCoord(world1, unit=u.deg)
        sep_arcsec = sc0.separation(sc1).arcsec[0]
        diff_mas = np.abs(sep_arcsec - north_arm_arcsec) * 1000
        ok = (diff_mas < 1)  # difference less than one milliarcsec?
        assert(ok), '***ERROR*** diff_mas >= 1 mas [0]'
        print()
        print('*****NORTH ARM*****')
        print(pixcrd0, '=pixcrd0')
        print(pixcrd1, '=pixcrd1')
        print(n_dt, '=n_dt =(n_x1-n_x0)  [top]')
        print(n_db, '=n_db =(n_y1-n_y0)  [bottom]')
        print(n_pa_deg, '=n_pa_deg')
        print(world0, '=world0 : center (RA,DEC) [deg]')
        print(world1, '=world1 : North arm tip (RA,DEC) [deg]')
        print(sc0, '=sc0')
        print(sc1, '=sc1')
        print(sep_arcsec, '=sep_arcsec')
        print(north_arm_arcsec, '=north_arm_arcsec')
        print(diff_mas, '=diff_mas')
        print(ok, '=ok')
    # pass:if

    # EAST arm of compass rose
    pixcrd0 = np.array([[cx0, cy0]], dtype=np.float_)
    # ^--- pixcrd0 must be a numpy 2-d array
    # pixels --> right ascension and declination:
    world = wcs.wcs_pix2world(pixcrd0, 0)
    declination = world[0][1]
    world[0][0] += east_arm_deg / np.cos(np.deg2rad(declination))
    # right ascension and declination --> pixels:
    pixcrd1 = wcs.wcs_world2pix(world, 0)
    e_x0 = pixcrd0[0][0] + llx  # window --> CCD
    e_y0 = pixcrd0[0][1] + lly  # window --> CCD
    e_x1 = pixcrd1[0][0] + llx  # window --> CCD
    e_y1 = pixcrd1[0][1] + lly  # window --> CCD

    # sanity check
    if (verbose):
        world0 = wcs.wcs_pix2world(pixcrd0, 0)
        world1 = wcs.wcs_pix2world(pixcrd1, 0)
        sc0 = ac.SkyCoord(world0, unit=u.deg)
        sc1 = ac.SkyCoord(world1, unit=u.deg)
        sep_arcsec = sc0.separation(sc1).arcsec[0]
        diff_mas = np.abs(sep_arcsec - east_arm_arcsec) * 1000
        ok = (diff_mas < 1)  # difference less than one milliarcsec?
        assert(ok), '***ERROR*** diff_mas >= 1 mas [1]'
        print()
        print('*****EAST ARM*****')
        print(pixcrd0, '=pixcrd0')
        print(pixcrd1, '=pixcrd1')
        print(world0, '=world0 : center (RA,DEC) [deg]')
        print(world1, '=world1 : East arm tip (RA,DEC) [deg]')
        print(sc0, '=sc0')
        print(sc1, '=sc1')
        print(sep_arcsec, '=sep_arcsec')
        print(east_arm_arcsec, '=east_arm_arcsec')
        print(diff_mas, '=diff_mas')
        print(ok, '=ok')
    # pass:if

    if (draw):
        # draw edge of compass rose with thick lines
        line = plt.Line2D((n_x0, n_x1), (n_y0, n_y1), lw=edge_lw, color=edge_color)
        ax.add_line(line)
        line = plt.Line2D((e_x0, e_x1), (e_y0, e_y1), lw=edge_lw, color=edge_color)
        ax.add_line(line)
        # draw middle of compass rose with thin lines
        line = plt.Line2D(
          (n_x0, n_x1), (n_y0, n_y1), lw=inside_lw, color=inside_color)
        ax.add_line(line)
        line = plt.Line2D(
          (e_x0, e_x1), (e_y0, e_y1), lw=inside_lw, color=inside_color)
        ax.add_line(line)
    # pass:if

    # sanity checks:
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
    north_top_half = (n_y1 > n_y0)  # NORTH ARM top is ABOVE of center?
    east_left_half = (e_x1 < e_x0)  # EAST ARM tip LEFT is LEFT of center?
    delta_deg = (positionAngle_deg - n_pa_deg)

    # HACK: BEGIN : create *new* attributes for the axis ======================
    ax.compass_positionAngle_deg = positionAngle_deg  # HACK: from CD matrix
    ax.compass_n_pa_deg = n_pa_deg  # HACK: from north compass arm
    ax.compass_mirrored = mirrored  # HACK" from CD matrix
    ax.compass_north_top_half = north_top_half  # HACK: from north compass arm
    ax.compass_east_left_half = east_left_half  # HACK: from east compass arm
    # HACK: END ===============================================================
    if (verbose):
        print()
        print('--->', ax.compass_positionAngle_deg, '=ax.compass_positionAngle_deg')
        print('--->', ax.compass_n_pa_deg, '=ax.compass_n_pa_deg')
        print('--->', ax.compass_mirrored, '=ax.compass_mirrored')
        print('--->', ax.compass_north_top_half, '=ax.compass_north_top_half')
        print('--->', ax.compass_east_left_half, '=ax.compass_east_left_half')
    # pass:if

    # sanity check
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
    # pass:if

    # that's all folks!
    if (verbose):
        print('\n==============================================================')
        print('END: %s()\n' % (func_))
    # pass:if
    return
# pass:def


# =============================================================================


def xmkpy3_plot_add_compass_rose_v5():
    import matplotlib.pyplot as plt
    import lightkurve as lk
    #
    obj = 3  # USER CUSTOMIZE
    #
    tpf = None
    if (obj == 1):
        north_arm_arcsec = 8
        tpf = lk.search_targetpixelfile(
          target='Kepler-138b', mission='Kepler', quarter=10)\
          .download(quality_bitmask=0)
        #         ^--- Exoplanet Kelper-138b is "KIC 7603200"
    elif (obj == 2):
        north_arm_arcsec = 42
        search_results = lk.search_tesscut(target='CD Ind', sector=1)
        tpf = search_results[0].download(cutout_size=(11, 11), quality_bitmask=0)
    elif (obj == 3):
        north_arm_arcsec = 42
        tpf = lk.search_targetpixelfile(target='CD Ind', mission='TESS',
          sector=1).download(quality_bitmask=0)
    else:
        print('***ERROR*** BAD OBJECT NUMBER:', obj)
    # pass:if
    assert(tpf is not None)
    #
    # Plot the 2nd frame of the TPF
    ax = tpf.plot(frame=1, cmap='gray_r')
    #
    # add a compass rose using the WCS from the TargetPixelFile
    mkpy3_plot_add_compass_rose_v5(ax=ax, wcs=tpf.wcs,
        north_arm_arcsec=north_arm_arcsec, verbose=True)
    #
    print(tpf.path, '\n^--- tpf.path\n')
    print(tpf.ra, '=tpf.ra')
    print(tpf.dec, '=tpf.dec')
    print(obj, '=obj')
    #
    plot_file = 'mkpy3_plot1.png'
    plt.savefig(plot_file, bbox_inches="tight")
    plt.close()
    print('\n', plot_file, ' <--- new PNG file written')
# pass:def


# =============================================================================


if (__name__ == '__main__'):
    xmkpy3_plot_add_compass_rose_v5()
# pass:if

# EOF
