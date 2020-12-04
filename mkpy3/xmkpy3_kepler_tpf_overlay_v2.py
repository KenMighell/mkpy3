#!/usr/bin/env python3

# file://xmkpy3_kepler_tpf_overlay_v2.py

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute

# =============================================================================


def xmkpy3_kepler_tpf_overlay_v2():
    """
Unit test
    """
    import os
    import sys
    import ntpath
    import argparse
    import ast
    import lightkurve as lk
    #
    import mkpy3
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
        '--width_height_arcmin', action="store", type=float, default=2.0,
        help='Width and height size in arcmin (float) [default: 2.0]')
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
        '--percentile', action="store", type=float, default=99.0,
        help='Percentile [percentage of pixels to keep: 0.0 to 100.0] '
        '(float) [default: 99.0]')
    parser.add_argument(
        '--cmap', action="store", type=str, default=None,
        help="Colormap name [Matplotlib] (str) [default: 'gray_r']")
    parser.add_argument(
        '--colors_str', action="store",
        type=ast.literal_eval, default="[None,'cornflowerblue','red']",
        help="string of a 3-item list of overlay color names [Matplotlib] "
        "(str) [default: \"['None','cornflowerblue','red']\"")
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
    kwargs_ = "{'edgecolor':'cyan', 's':300, 'facecolor':'None', 'lw':3, "\
        "'zorder':20}"
    parser.add_argument(
        '--gaia_dr2_kwargs_str', action="store",
        type=ast.literal_eval, default=kwargs_,
        help="GAIA DR2 marker kwargs (string of a dictonary) for ax.scatter() "
        "[Matplotlib] "'(str) [default: "' + kwargs_ + '"')
    kwargs_ = "{'s':900, 'color':'lawngreen', 'marker':'x', 'lw':5, "\
        "'zorder':30}"
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
    gaia_dr2_kwargs_str = str(args.gaia_dr2_kwargs_str)
    vsx_kwargs_str = str(args.vsx_kwargs_str)
    sexagesimal = args.sexagesimal
    verbose = args.verbose

    print()
    if (tpf_filename is not None):
        mkpy3.mkpy3_util_check_file_exists(tpf_filename, True)
        tpf = lk.open(tpf_filename)
    else:
        print('No TargetPixelFile (TPF) filename given.\n')
        tpf = lk.search_targetpixelfile(
            target='kepler-138b', mission='kepler', cadence='long',
            quarter=10).download(quality_bitmask=0)
        # ^--- exoplanet Kelper-138b is "KIC 7603200"
        print()
        print(
            'Using default TPF [Kepler Q10 observations of exoplanet Kepler'
            '-138b (KIC 760320)].')
        print()
        shrink *= 0.8
    # pass:if
    try:
        print('TPF filename:', ntpath.basename(tpf.path))
        print('TPF dirname: ', os.path.dirname(tpf.path))
        assert(tpf.mission == 'Kepler')
        print()
    except Exception:
        print(tpf_filename, '=tpf_filename')
        print('^--- *** ERROR *** This file does not appear to be a Kepler '
              'TargetPixelFile')
        print()
        print('Bye...\n', flush=True)
        sys.exit(1)
    # pass:try

    ax = mkpy3.mkpy3_tpf_overlay_v6(
      tpf=tpf,
      frame=frame,
      survey=survey,
      width_height_arcmin=width_height_arcmin,
      shrink=shrink,
      show_plot=show_plot,
      plot_file=plot_file,
      overwrite=overwrite,
      figsize_str=figsize_str,
      title=title,
      percentile=percentile,
      cmap=args.cmap,
      colors_str=colors_str,
      lws_str=lws_str,
      zorders_str=zorders_str,
      marker_kwargs_str=marker_kwargs_str,
      gaia_dr2_kwargs_str=gaia_dr2_kwargs_str,
      vsx_kwargs_str=vsx_kwargs_str,
      sexagesimal=sexagesimal,
      verbose=verbose
    )
# pass:def


if (__name__ == '__main__'):
    xmkpy3_kepler_tpf_overlay_v2()
# pass:if

# EOF
