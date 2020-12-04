#!/usr/bin/env python3

# file://mkpy3_tess_tpf_overlay_v6.py

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute

# =============================================================================

def mkpy3_tess_tpf_overlay_v6(
    tpf=None,
    frame=0,
    survey='2MASS-J',
    rotationAngle_deg=None,  # None or '123.456' (float) or 'tpf'
    width_height_arcmin=6.0,
    shrink=1.0,
    show_plot=True,
    plot_file='mkpy3_plot.png',
    overwrite=False,
    figsize_str='[9,9]',
    title=None,
    percentile=99.5,
    cmap='gray_r',
    colors_str="[None,'dodgerblue','red']",
    lws_str='[0,3,4]',
    zorders_str='[0,2,4]',
    marker_kwargs_str="{'edgecolor':'yellow', 's':600, 'facecolor':'None', "
    "'lw':3, 'zorder':10}",  # or 'None'
    print_gaia_dr2=True,
    gaia_dr2_kwargs_str="{'edgecolor':'cyan', 's':150, 'facecolor':'None', "
    "'lw':3, 'zorder':20}",  # or 'None'
    print_vsx=True,
    vsx_kwargs_str="{'s':900, 'color':'lawngreen', 'marker':'x', 'lw':5, "
    "'zorder':30}",  # or 'None'
    sexagesimal=False,
    verbose=False
):
    """
Function: mkpy3_tess_tpf_overlay_v6()

Purpose: Plot a TESS TargetPixelFile (TPF) overlay on a sky survey image.

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
rotationAngle_deg : None, (float), or 'tpf' (3-char str) (optional)
    Angle in degrees to rotate the sky survey image.
    [default: None ---> 'tpf']
    [example values: 'None' or 12.345 (float) or 'tpf']
width_height_arcmin : (float) (optional)
    Width and height of the survey image [arcmin].
    [default: 6.0]
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
    [default: 99.5]
cmap : (str) (optional)
    Colormap name [Matplotlib].
    [default: 'gray_r']
colors_str : (str) (optional)
    A string of a 3-item list of overlay color names [Matplotlib].
    [default: "['None','dodgerblue','red']"]
lws_str : (str) (optional)
    A string of a 3-item list of overlay line widths [Matplotlib].
    [default: '[0,3,4]']
zorders_str : (str) (optional)
    A string of a 3-item list of overlay zorder values [Matplotlib].
    [default: '[0,2,4]']
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
    [default: "{'edgecolor':'cyan', 's':150, 'facecolor':'None', 'lw':3,
    'zorder':20}"]
print_vsx : (bool) (optional)
    If True, print the VSX catalog results.
    [default=True]
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

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute
    """
    import ntpath
    import os
    import sys
    import lightkurve as lk

    import mkpy3

    if (tpf is None):
        target = 'V1460 Her'
        sector = 24
        # alternative:
        # target = 'XZ Cyg'
        # sector = 14
        title = target + ' : TESS : Sector ' + str(sector)
        search_result = lk.search_tesscut(target=target, sector=sector)[0]
        tpf = search_result.download(cutout_size=(11, 13), quality_bitmask=0)
    # pass:if

    assert(tpf is not None)
    tpf_filename = ntpath.basename(tpf.path)
    tpf_dirname = os.path.dirname(tpf.path)
    try:
        print()
        print('TPF filename:', tpf_filename)
        print('TPF dirname: ', tpf_dirname)
        assert(tpf.mission == 'TESS')
        print()
    except Exception:
        print(tpf_filename, '=tpf_filename')
        print('^--- *** ERROR *** This file does not appear to be a TESS '
          'TargetPixelFile')
        print()
        print('Bye...\n', flush=True)
        sys.exit(1)
    # pass:try

    ax = mkpy3.mkpy3_tpf_overlay_v6(
        tpf=tpf,
        frame=frame,
        survey=survey,
        rotationAngle_deg=rotationAngle_deg,
        width_height_arcmin=width_height_arcmin,
        shrink=shrink,
        show_plot=show_plot,
        plot_file=plot_file,
        overwrite=overwrite,
        figsize_str=figsize_str,
        title=title,
        percentile=percentile,
        cmap=cmap,
        colors_str=colors_str,
        lws_str=lws_str,
        zorders_str=zorders_str,
        marker_kwargs_str=marker_kwargs_str,
        print_gaia_dr2=print_gaia_dr2,
        gaia_dr2_kwargs_str=gaia_dr2_kwargs_str,
        print_vsx=print_vsx,
        vsx_kwargs_str=vsx_kwargs_str,
        sexagesimal=sexagesimal,
        verbose=verbose)

    return ax
# pass:def


# =============================================================================


def xmkpy3_tess_tpf_overlay_v6():
    """
Unit test
    """
    import argparse
    import ast
    import lightkurve as lk
    import mkpy3
    #
    #
    # ===== argparse:BEGIN ====================================================
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
        '--rotationAngle_deg', action="store",
        type=ast.literal_eval, default=None,
        help="Rotation angle in degrees (string) [default: None] "
        "[examples: None or 12.345 (float) or 'tpf'")
    parser.add_argument(
        '--width_height_arcmin', action="store", type=float, default=6.0,
        help='Width and height size in arcmin (float) [default: 6.0]')
    parser.add_argument(
        '--shrink', type=float, default=1.0,
        help='Survey search radius shrink factor (float) [default: 1.0]')
    parser.add_argument(
        '--show_plot', type=mkpy3.mkpy3_util_str2bool, default=True,
        help='If True, show the plot [default=True]')
    parser.add_argument(
        '--plot_file', action="store", type=str, default='mkpy3_plot.png',
        help='Filename of the output plot [default: "mkpy3_plot.png"]')
    parser.add_argument(
        '--overwrite', type=mkpy3.mkpy3_util_str2bool, default=False,
        help='If True, overwrite ("clobber") an existing output file '
        '[default: False]')
    parser.add_argument(
        '--figsize_str', action="store",
        type=ast.literal_eval, default="[9,9]",
        help="string of a 2-item list of figure width and height [Matplotlib] "
        "(str) [default: '[9,9]'")
    parser.add_argument(
        '--title', action="store", type=str, default=None,
        help='Title of the finder chart (str) [default: None]')
    parser.add_argument(
        '--percentile', action="store", type=float, default=99.5,
        help='Percentile [percentage of pixels to keep: 0.0 to 100.0] '
        '(float) [default: 99.5]')
    parser.add_argument(
        '--cmap', action="store", type=str, default=None,
        help="Colormap name [Matplotlib] (str) [default: 'gray_r']")
    parser.add_argument(
        '--colors_str', action="store",
        type=ast.literal_eval, default="[None,'dodgerblue','red']",
        help="string of a 3-item list of overlay color names [Matplotlib] "
        "(str) [default: \"['None','dodgerblue','red']\"")
    parser.add_argument(
        '--lws_str', action="store",
        type=ast.literal_eval, default="[0,3,4]",
        help="string of a 3-item list of overlay line widths [Matplotlib] "
        "(str) [default: \"[0,3,4]\"")
    parser.add_argument(
        '--zorders_str', action="store",
        type=ast.literal_eval, default="[0,2,4]",
        help="string of a 3-item list of overlay zorder values [Matplotlib] "
        "(str) [default: \"[0,2,4]\"")
    kwargs_ = "{'edgecolor':'yellow', 's':600, 'facecolor':'None', 'lw':3, "\
        "'zorder':10}"
    parser.add_argument(
        '--marker_kwargs_str', action="store",
        type=ast.literal_eval, default=kwargs_,
        help="marker kwargs (string of a dictonary) for ax.scatter() "
        "[Matplotlib] " + '(str) [default: "' + kwargs_ + '"')
    kwargs_ = "{'edgecolor':'cyan', 's':150, 'facecolor':'None', 'lw':3, "\
        "'zorder':20}"
    parser.add_argument(
        '--print_gaia_dr2', type=mkpy3.mkpy3_util_str2bool, default=True,
        help='If True, print the GAIA DR2 catalog results [default=True]')
    parser.add_argument(
        '--gaia_dr2_kwargs_str', action="store",
        type=ast.literal_eval, default=kwargs_,
        help="GAIA DR2 marker kwargs (string of a dictonary) for ax.scatter() "
        "[Matplotlib] "'(str) [default: "' + kwargs_ + '"')
    kwargs_ = "{'s':900, 'color':'lawngreen', 'marker':'x', 'lw':5, "\
        "'zorder':30}"
    parser.add_argument(
        '--print_vsx', type=mkpy3.mkpy3_util_str2bool, default=True,
        help='If True, print the VSX catalog results [default=True]')
    parser.add_argument(
        '--vsx_kwargs_str', action="store",
        type=ast.literal_eval, default=kwargs_,
        help="VSX marker kwargs (string of a dictonary) for ax.scatter() "
        "[Matplotlib] (str) [default: '" + kwargs_ + "'")
    parser.add_argument(
        '--sexagesimal', type=mkpy3.mkpy3_util_str2bool, default=False,
        help='Print catalog positions as sexagesimal [hms dms] if True (bool) '
        '[default=False]')
    parser.add_argument(
        '--verbose', type=mkpy3.mkpy3_util_str2bool, default=False,
        help='Print extra information if True (bool) [default=False]')
    #
    args = parser.parse_args()
    #
    # ===== argparse:END ======================================================
    #

    tpf_filename = args.tpf_filename
    frame = args.frame
    survey = args.survey
    rotationAngle_deg = args.rotationAngle_deg
    width_height_arcmin = args.width_height_arcmin
    shrink = args.shrink
    show_plot = args.show_plot
    plot_file = args.plot_file
    overwrite = args.overwrite
    figsize_str = str(args.figsize_str)
    title = args.title
    percentile = args.percentile
    cmap = args.cmap
    colors_str = str(args.colors_str)
    lws_str = str(args.lws_str)
    zorders_str = str(args.zorders_str)
    marker_kwargs_str = str(args.marker_kwargs_str)
    print_gaia_dr2 = args.print_gaia_dr2
    gaia_dr2_kwargs_str = str(args.gaia_dr2_kwargs_str)
    print_vsx = args.print_vsx
    vsx_kwargs_str = str(args.vsx_kwargs_str)
    sexagesimal = args.sexagesimal
    verbose = args.verbose

    if (tpf_filename is not None):
        mkpy3.mkpy3_util_check_file_exists(tpf_filename, True)
        tpf = lk.open(tpf_filename, quality_bitmask=0)
    # pass:if

    if (tpf_filename is None):
        tpf = None
    # pass:if

    shrink = 0.4
    ax = mkpy3_tess_tpf_overlay_v6(
      tpf=tpf,
      frame=frame,
      survey=survey,
      rotationAngle_deg=rotationAngle_deg,
      width_height_arcmin=width_height_arcmin,
      shrink=shrink,
      show_plot=show_plot,
      plot_file=plot_file,
      overwrite=overwrite,
      figsize_str=figsize_str,
      title=title,
      percentile=percentile,
      cmap=cmap,
      colors_str=colors_str,
      lws_str=lws_str,
      zorders_str=zorders_str,
      marker_kwargs_str=marker_kwargs_str,
      print_gaia_dr2=print_gaia_dr2,
      gaia_dr2_kwargs_str=gaia_dr2_kwargs_str,
      print_vsx=print_vsx,
      vsx_kwargs_str=vsx_kwargs_str,
      sexagesimal=sexagesimal,
      verbose=verbose
    )
# pass:def


if (__name__ == '__main__'):
    xmkpy3_tess_tpf_overlay_v6()
# pass:if

# EOF
