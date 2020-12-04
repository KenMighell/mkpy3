#!/usr/bin/env python3

# file://mkpy3_finder_chart_survey_fits_image_get_v1.py

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute

# =============================================================================


def mkpy3_finder_chart_survey_fits_image_get_v1(
  ra_deg=None,
  dec_deg=None,
  radius_arcmin=None,
  survey=None,
  cframe=None,
  verbose=None
):
    """
Function: mkpy3_finder_chart_survey_fits_image_get_v1()

Purpose:

Gets sky survey image data around a position on the sky.

Parameters
----------
ra_deg : float (optional)
    right ascencsion [deg]
dec_deg : float (optional)
    declination [deg]
radius_arcmin : float (optional)
    radius (halfwidth and halfheight of image) [arcmin]
survey : string (optional) [e.g., '2MASS-J', 'DSS2 Red', etc.]
    survey string name
cframe : str (optional)
    coordinate frame name [e.g., 'fk5', 'icrs', etc.]
verbose : bool (optional)
    if True, print extra information

Returns
-------
hdu :
    Header/Data Unit (HDU) of the survey FITS file
hdr :
    header associated with hdu
data :
    data associated with hdu
wcs :
    World Coordinate System from hdu
cframe :
    coordinate frame of the survey data

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute
    """
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.skyview import SkyView
    from astropy.wcs import WCS
    #
    if (ra_deg is None):
        ra_deg = 291.41829  # Kepler-93b
    if (dec_deg is None):
        dec_deg = 38.67236  # Kepler-93b
    if (radius_arcmin is None):
        radius_arcmin = 1.99
    if (survey is None):
        survey = '2MASS-J'  # alternate: 'DSS2 Red'
        # ^--- to see all surveys: astroquery.skyview.SkyView.list_surveys()
    # pass:if
    if (cframe is None):
        cframe = 'fk5'  # N.B.: '2MASS-J' uses 'fk5'
    if (verbose is None):
        verbose = False
    #
    if (verbose):
        print(ra_deg, '=ra_deg')
        print(dec_deg, '=dec_deg')
        print(radius_arcmin, '=radius_arcmin')
        print("'%s' =survey" % (survey))
        print("'%s' =cframe" % (cframe))
        print(verbose, '=verbose')
        print()
    # pass#if
    #
    # sc <--- astropy sky coordinates
    sc = SkyCoord(ra=ra_deg * u.degree, dec=dec_deg * u.degree, frame=cframe)
    # image list  # assume that the list contains a single image
    imgl = SkyView.get_images(position=sc, survey=survey, radius=radius_arcmin * u.arcmin)
    #
    # outputs:
    hdu = imgl[0]        # Header/Data Unit of the FITS image
    hdr = hdu[0].header  # header associated with the HDU
    data = hdu[0].data   # data associated with the HDU
    wcs = WCS(hdr)       # World Coordinate System from the FITS header of the survey image
    #
    return hdu, hdr, data, wcs, cframe
# pass#def


def xmkpy3_finder_chart_survey_fits_image_get_v1():
    import lightkurve as lk
    lk.log.setLevel('INFO')
    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.visualization import ImageNormalize, PercentileInterval, SqrtStretch
    import os
    import ntpath

    # Exoplanet Kelper-138b is "KIC 7603200":
    tpf = lk.search_targetpixelfile(target='kepler-138b', mission='kepler',
      cadence='long', quarter=10).download(quality_bitmask=0)
    print('TPF filename:', ntpath.basename(tpf.path))
    print('TPF dirname: ', os.path.dirname(tpf.path))

    target = 'Kepler-138b'
    ra_deg = tpf.ra
    dec_deg = tpf.dec

    # get survey image data
    width_height_arcmin = 3.00
    survey = '2MASS-J'
    survey_hdu, survey_hdr, survey_data, survey_wcs, survey_cframe = \
      mkpy3_finder_chart_survey_fits_image_get_v1(ra_deg, dec_deg,
      radius_arcmin=width_height_arcmin, survey=survey, verbose=True)

    # create a matplotlib figure object
    fig = plt.figure(figsize=(12, 12))

    # create a matplotlib axis object with right ascension and declination axes
    ax = plt.subplot(projection=survey_wcs)

    norm = ImageNormalize(survey_data, interval=PercentileInterval(99.0),
      stretch=SqrtStretch())
    ax.imshow(survey_data, origin='lower', norm=norm, cmap='gray_r')

    ax.set_xlabel('Right Ascension (J2000)')
    ax.set_ylabel('Declination (J2000)')
    ax.set_title('')
    plt.suptitle(target)

    # put a yellow circle at the target position
    ax.scatter(ra_deg * u.deg, dec_deg * u.deg,
      transform=ax.get_transform(survey_cframe),
      s=600, edgecolor='yellow', facecolor='None', lw=3, zorder=100)

    pname = 'mkpy3_plot.png'
    if (pname != ''):
        plt.savefig(pname, bbox_inches="tight")
        print(pname, ' <--- plot filename has been written!  :-)\n')
    # pass:if
# pass:def


if (__name__ == '__main__'):
    xmkpy3_finder_chart_survey_fits_image_get_v1()
# pass:if

# EOF
