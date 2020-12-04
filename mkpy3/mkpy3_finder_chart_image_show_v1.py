#!/usr/bin/env python3

# file://mkpy3_finder_chart_image_show_v1.py

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute

# =============================================================================


def mkpy3_finder_chart_image_show_v1(
  ax=None,
  image_data=None,
  percentile=None,
  cmap=None,
  title=None,
  verbose=None
):
    """
Function: mkpy3_finder_chart_image_show()

Purpose:

Plot the sky survey image data using matplotlib's ax.imshow.

Parameters
----------
ax : matplotlib axes.Axes object
    axis object
image_data : 2-d numpy.ndarray
    2-d survey image data array
percentile : float [0.0 to 100.0] (optional)
    percentile (percent pixels to keep)
cmap : string [matplotlib colormap name] (optional)
    cmap name [e.g., 'gray_r']
title : string (optional)
    title string of the plot
verbose : bool (optional)
    if True, print extra information

Returns: nothing

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute
    """
    from astropy.visualization import ImageNormalize, PercentileInterval,\
        SqrtStretch
    #
    assert(ax is not None)
    assert(image_data is not None)
    if (percentile is None):
        percentile = 99
    # pass:if
    if (cmap is None):
        cmap = 'gray_r'
    # pass:if
    if (title is None):
        title = ''
    # pass:if
    if (verbose is None):
        verbose = False
    # pass:if
    if (verbose):
        print(percentile, '=percentile')
        print("'%s' =cmap" % (cmap))
        print(verbose, ' =verbose')
        print()
    # pass#if
    #
    # =========================================================================
    # show the survey image with right ascension and declination coordinates
    norm = ImageNormalize(
        image_data, interval=PercentileInterval(percentile),
        stretch=SqrtStretch())
    ax.imshow(image_data, origin='lower', norm=norm, cmap=cmap)
    #
    # set plot parameters
    ax.tick_params(axis='x', labelsize=16, length=6, width=2)
    ax.tick_params(axis='y', labelsize=16, length=6, width=2)
    ax.set_xlabel('Right Ascension (J2000)', size=24)
    ax.set_ylabel('Declination (J2000)', size=24)
    if (title != ''):
        ax.set_title(title, size=20)
    # pass:if
    ax.set_aspect(1)
    #
    return
# pass:def


# =============================================================================


def xmkpy3_finder_chart_image_show_v1():
    import matplotlib.pyplot as plt
    import astropy.units as u
    import os
    import ntpath
    import lightkurve as lk
    lk.log.setLevel('INFO')

    import mkpy3

    # Exoplanet Kelper-138b is "KIC 7603200":
    tpf = lk.search_targetpixelfile(
        target='kepler-138b', mission='kepler', cadence='long',
        quarter=10).download(quality_bitmask=0)
    print('TPF filename:', ntpath.basename(tpf.path))
    print('TPF dirname: ', os.path.dirname(tpf.path))

    target = 'Kepler-138b'
    title_ = tpf.hdu[0].header['object']
    title_ += ' : Exoplanet ' + target

    ra_deg = tpf.ra
    dec_deg = tpf.dec

    # get survey image data
    width_height_arcmin = 3.00
    survey = '2MASS-J'
    survey_hdu, survey_hdr, survey_data, survey_wcs, survey_cframe = \
        mkpy3.mkpy3_finder_chart_survey_fits_image_get_v1(
            ra_deg, dec_deg, radius_arcmin=width_height_arcmin, survey=survey,
            verbose=True)

    # create a matplotlib figure object
    fig = plt.figure(figsize=(12, 12))

    # create a matplotlib axis object with right ascension and declination axes
    ax = plt.subplot(projection=survey_wcs)

    # show the survey image
    percentile = 99.0
    mkpy3_finder_chart_image_show_v1(
        ax=ax, image_data=survey_data, percentile=percentile, verbose=True)

    # set suptitle
    ax.set_title('', size=24)  # NUM, ME VEXO?
    plt.suptitle(title_, size=24)

    # put a yellow circle at the target position
    ax.scatter(
        ra_deg * u.deg, dec_deg * u.deg, transform=ax.get_transform(survey_cframe),
        s=600, edgecolor='yellow', facecolor='None', lw=3, zorder=10)

    pname = 'mkpy3_plot.png'
    if (pname != ''):
        plt.savefig(pname, bbox_inches="tight")
        print(pname, ' <--- plot filename has been written!  :-)\n')
    # pass:if
# pass:def


# =============================================================================


if (__name__ == '__main__'):
    xmkpy3_finder_chart_image_show_v1()
# pass:if

# EOF
