#!/usr/bin/env python3

# Kenneth John Mighell
# Kepler Support Scientist
# Kepler / K2 Science Office
# NASA Ames Research Center / SETI Institute

__version__ = '2020SEP27T1821 v0.23'


def mkpy3_plot_add_compass_rose_v4(
  ax=None,
  wcs=None,
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
Function : mkpy3_plot_add_compass_rose_v4()

Purpose: Add a compass rose to a matplotlib axis.
    NOTE: Long arm points North (increasing declination)
    NOTE: Short arm points East (increasing right ascension)

Parameters
----------
ax : matplotlib Axes object
wcs : astropy FITS World Coordinate System (WCS) object (optional)
    [default: None --> ax.wcs]
    [N.B. with lightkurve TPF plots use wcs=tpf.wcs]
cx : float (optional)
    pixel column number (X position) of the center of the compass rose
cy : float (optional)
    pixel row number (Y position) of the center of the compass rose
north_arm_arcsec : float (optional) [default: 6 arcsec]
    length of the North arrow arm of the compass rose in arcsec
    [N.B. A Kepler photometer CCD pixel is 3.98 arcsec/pixel]
    [N.B. A TESS photometer CCD pixel is 21.00 arcsec/pixel]
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
mkpy3_plot_add_compass_rose_v4(ax=ax, wcs=tpf.wcs, verbose=True)
#
plot_file = 'mkpy3_plot.png'
plt.savefig(plot_file, bbox_inches="tight")
plt.close()
print(plot_file, ' <--- new PNG file written')
#==========

# Kenneth John Mighell
# Kepler Support Scientist
# Kepler / K2 Science Office
# NASA Ames Research Center / SETI Institute
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import astropy.coordinates as ac
    import astropy.units as u

    assert(ax is not None), '***ERROR***: Requires a matplotlib axes object'
    if (wcs is None):
        assert(hasattr(ax, 'wcs')),\
          '***ERROR***: Requires that ax.wcs exists (is an attribute)'
        wcs = ax.wcs
    # pass:if
    if (cx is None):
        xlim = ax.get_xlim()
        cx = (xlim[0] + xlim[1]) / 2.0  # center of the X axis
    # pass:if
    if (cy is None):
        ylim = ax.get_ylim()
        cy = (ylim[0] + ylim[1]) / 2.0  # center of the Y axis
    # pass:if
    if (north_arm_arcsec is None):
        north_arm_arcsec = 6  # default for Kepler/K2 observations
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
        print('^--- =wcs  (NOTE: ignore NAXIS values)')
        print(cx, '=cx')
        print(cy, '=cy')
        print(north_arm_arcsec, '=north_arm_arcsec')
        print(edge_color, '=edge_color')
        print(inside_color, '=inside_color')
        print(edge_lw, '=edge_lw')
        print(inside_lw, '=inside_lw')
        print(verbose, '=verbose')
        print()
        print('==============================================================')
    # pass:if

    cx0 = cx  # zero-offset coordinates
    cy0 = cy  # zero-offset coordinates
    north_arm_deg = north_arm_arcsec / 3600.0
    # east arm is half as long as north arm:
    east_arm_deg = north_arm_deg / 2.0
    east_arm_arcsec = east_arm_deg * 3600.0

    if (verbose):
        print()
        print(north_arm_arcsec, '=north_arm_arcsec')
        print(east_arm_arcsec, '=east_arm_arcsec')
    # pass:if

    # north arm of compass rose
    pixcrd0 = np.array([[cx0, cy0]], dtype=np.float_)
    # ^--- pixcrd0 must be a numpy 2-d array
    # pixels --> right ascension and declination:
    world = wcs.wcs_pix2world(pixcrd0, 0)
    world[0][1] += north_arm_deg
    # right ascension and declination --> pixels:
    pixcrd1 = wcs.wcs_world2pix(world, 0)
    n_x0 = pixcrd0[0][0]
    n_y0 = pixcrd0[0][1]
    n_x1 = pixcrd1[0][0]
    n_y1 = pixcrd1[0][1]

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
        print(pixcrd0, '=pixcrd0')
        print(pixcrd1, '=pixcrd1')
        print(world0, '=world0 : center (RA,DEC) [deg]')
        print(world1, '=world1 : North arm tip (RA,DEC) [deg]')
        print(sc0, '=sc0')
        print(sc1, '=sc1')
        print(sep_arcsec, '=sep_arcsec')
        print(north_arm_arcsec, '=north_arm_arcsec')
        print(diff_mas, '=diff_mas')
        print(ok, '=ok')
    # pass:if

    # east arm of compass rose
    pixcrd0 = np.array([[cx0, cy0]], dtype=np.float_)
    # ^--- pixcrd0 must be a numpy 2-d array
    # pixels --> right ascension and declination:
    world = wcs.wcs_pix2world(pixcrd0, 0)
    declination = world[0][1]
    world[0][0] += east_arm_deg / np.cos(np.deg2rad(declination))
    # right ascension and declination --> pixels:
    pixcrd1 = wcs.wcs_world2pix(world, 0)
    e_x0 = pixcrd0[0][0]
    e_y0 = pixcrd0[0][1]
    e_x1 = pixcrd1[0][0]
    e_y1 = pixcrd1[0][1]

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

    if (verbose):
        print()
        print('==============================================================')
        print()
    # pass:if

    return
# pass:def


if (__name__ == '__main__'):
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
    mkpy3_plot_add_compass_rose_v4(ax=ax, wcs=tpf.wcs, verbose=True)
    #
    plot_file = 'mkpy3_plot.png'
    plt.savefig(plot_file, bbox_inches="tight")
    plt.close()
    print(plot_file, ' <--- new PNG file written')
# pass:if
# EOF
