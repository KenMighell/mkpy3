#!/usr/bin/env python3

# Kenneth John Mighell
# Kepler Support Scientist
# Kepler / K2 Science Office
# NASA Ames Research Center / SETI Institute

__version__ = '2020SEP19T0851 0.17'


def mkpy3_tpf_plot_add_compass_rose_v2(
  ax=None,
  tpf=None,
  cx=None,
  cy=None,
  north_arrow_length_arcsec=None,
  edge_color=None,
  inside_color=None,
  edge_lw=None,
  inside_lw=None,
  verbose=None
):
    """
Function : mkpy3_tpf_plot_add_compass_rose_v2()

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
  target='Kepler-138b',mission='Kepler',quarter=10).download()
#         ^---> Exoplanet Kelper-138b is "KIC 7603200"
#
# Plot the 2nd frame of the TPF
ax = tpf.plot(frame=1)
#
# add a compass rose
mkpy3_tpf_plot_add_compass_rose_v2(ax=ax, tpf=tpf)
#
plt.savefig('mkpy3_spud.png', bbox_inches="tight")
#==========

# Kenneth John Mighell
# Kepler Support Scientist
# Kepler / K2 Science Office
# NASA Ames Research Center / SETI Institute
    """
    import mkpy3_wcs_plot_add_compasss_rose_v1 as km1

    assert(ax is not None)
    assert(tpf is not None)
    if (cx is None):
        xlim = ax.get_xlim()
        cx = (xlim[0] + xlim[1]) / 2.0  # center of the X axis
    # pass:if
    if (cy is None):
        ylim = ax.get_ylim()
        cy = (ylim[0] + ylim[1]) / 2.0  # center of the Y axis
    # pass:if
    if (north_arrow_length_arcsec is None):
        north_arrow_length_arcsec = 6  # default for Kepler/K2 observations
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
        print(ax, '=ax')
        print(tpf, '=tpf')
        print(cx, '=cx')
        print(cy, '=cy')
        print(north_arrow_length_arcsec, '=north_arrow_length_arcsec')
        print(edge_color, '=edge_color')
        print(inside_color, '=inside_color')
        print(edge_lw, '=edge_lw')
        print(inside_lw, '=inside_lw')
        print(verbose, '=verbose')
    # pass:if

    wcs = tpf.wcs
    km1.mkpy3_wcs_plot_add_compass_rose_v1(
      ax=ax,
      wcs=wcs,
      cx=cx,
      cy=cy,
      north_arrow_length_arcsec=north_arrow_length_arcsec,
      edge_color=edge_color,
      inside_color=inside_color,
      edge_lw=edge_lw,
      inside_lw=inside_lw,
      verbose=verbose
    )

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
    # add a compass rose
    mkpy3_tpf_plot_add_compass_rose_v2(ax=ax, tpf=tpf)
    #
    plot_file = 'mkpy3_plot.png'
    plt.savefig(plot_file, bbox_inches="tight")
    plt.close()
    print(plot_file, ' <--- new PNG file written')
# pass:if
# EOF
