#!/usr/bin/env python3

# file://mkpy3_vizier_catalog_cone_get_v4.py

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute

# =============================================================================


def mkpy3_vizier_catalog_cone_get_v4(
  ra_deg=None,
  dec_deg=None,
  radius_arcsec=None,
  vizier_catalog=None,
  verbose=None
):
    """
Function : mkpy3_vizier_catalog_cone_get_v4()

Purpose: Perform a cone search of a Vizier catalog using astroquery.

Parameters
----------
ra_deg : float (optional)
    right ascencsion of target [deg]
dec_deg : float (optional)
    declination of target [deg]
radius_arcsec : float (optional)
    search radius [arcsec]
vizier_catalog : str (optional)
    Vizier catalog name. Examples:
    'B/vsx/vsx'    <--- International Variable Star Index (VSX)
    'I/345/gaia2'  <--- GAIA DR2
verbose : bool (optional)
    verbose output [default: False]

Returns
-------
raj2000 : float array
    right ascension (J2000) [deg]
dej2000 : float array
    declination (J2000) [deg]
sep_arcsec : float array
    separation from target [arcsec]
vizier_catalog_result :
    catalog table returned by Vizier

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute
    """
    import numpy as np
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from astroquery.vizier import Vizier
    import copy
    #
    #
    assert(ra_deg is not None)
    assert(dec_deg is not None)
    assert(radius_arcsec is not None)
    assert(vizier_catalog is not None)
    if (verbose is None):
        verbose = False
    assert(ra_deg >= 0.0)
    assert(ra_deg < 360.0)
    assert(dec_deg >= -90.0)
    assert(dec_deg <= +90.0)
    assert(radius_arcsec > 0.0)
    assert(vizier_catalog is not None)
    if (verbose is None):
        verbose = False
    #
    # HACK: BEGIN =============================================================
    # fixes a nasty bug found in astroquery.vizier.Vizier:
    chunk = 10000.0
    radius_arcsec = np.rint(radius_arcsec * chunk) / chunk  # HACK A FURBALL!
    # search radius now rounded to 4 decimal places (~100 microdegrees)
    # HACK: END ===============================================================
    #
    if (verbose):
        print()
        print(ra_deg, '=ra_deg')
        print(dec_deg, '=dec_deg')
        print(radius_arcsec, '=radius_arcsec')
        print(vizier_catalog, '=vizier_catalog')
    # pass:if
    #
    # target
    sc = SkyCoord(ra_deg, dec_deg, frame='icrs', unit='deg')
    #
    radius = radius_arcsec * u.arcsec
    #
    v = Vizier(columns=["**", "+_r"], catalog=vizier_catalog)
    #     all columns ---^     ^---> add column for (increasing) angular
    #                                separation
    # N.B.: "*" <--- a *single* asterisk in the first argument gets only the
    # *default* columns
    v.ROW_LIMIT = -1  # no row limit
    result_query = v.query_region(sc, radius=radius)
    no_targets_found_message = ValueError(
        'Either no sources were found in the query region or Vizier is unavail'
        'able')
    if result_query is None:
        raise no_targets_found_message
    elif len(result_query) == 0:
        raise no_targets_found_message
    result = result_query[0]
    #
    # create output objects:
    RAJ2000 = np.array(result['RAJ2000'])
    DEJ2000 = np.array(result['DEJ2000'])
    sep_arcsec = np.array(result['_r'])
    vizier_result = copy.deepcopy(result)
    #
    if (verbose):
        print()
        print(vizier_result)
        print('^--- =vizier_result [dumps partial table (first ... last column'
              's)]')
        print()
        print(vizier_result.info)
        print('^--- vizier_result.info [show all columns]')
        print()
        print(RAJ2000)
        print('^--- =RAJ2000')
        print()
        print(DEJ2000)
        print('^--- =DEJ2000')
        print()
        print(sep_arcsec)
        print('^--- =sep_arcsec')
        print('\n#%s:' % (vizier_catalog))
        print('#index RAJ2000 DEJ000 sep_arcsec')
        for j in range(RAJ2000.size):
            print(j, RAJ2000[j], DEJ2000[j], sep_arcsec[j])
        # pass:for
    # pass:if
    return (RAJ2000, DEJ2000, sep_arcsec, vizier_result)
# pass:def


# =============================================================================


def xmkpy3_vizier_catalog_cone_get_v4():
    print()
    print('#')
    print('# RR Lyrae test ==================================================')
    print('#')
    ra_deg = 291.3663013467642   # RR Lyr
    dec_deg = +42.7843585094725  # RR Lyr
    radius_arcsec = 300
    vizier_catalog = 'B/vsx/vsx'
    verbose = True
    #
    proceed = True
    try:
        raj2000, dej2000, sep_arcsec, vizier_vsx_result = \
            mkpy3_vizier_catalog_cone_get_v4(
              ra_deg=ra_deg, dec_deg=dec_deg, radius_arcsec=radius_arcsec,
              vizier_catalog=vizier_catalog, verbose=verbose)
    except Exception:
        proceed = False
    # pass:try
    #
    assert(proceed is True)
    print('\n#VSX:')
    print('#index raj2000 dej2000 sep_arcsec')
    for j in range(raj2000.size):
        print(j, raj2000[j], dej2000[j], sep_arcsec[j])
    # pass:for
    #
    # =========================================================================
    print()
    print('#')
    print('# empty catalog test  ============================================')
    print('#')
    ra_deg = 0.0
    dec_deg = 0.0
    radius_arcsec = 0.1
    vizier_catalog = 'B/vsx/vsx'
    verbose = True
    #
    proceed = True
    try:
        raj2000, dej2000, sep_arcsec, vizier_vsx_result = \
            mkpy3_vizier_catalog_cone_get_v4(
              ra_deg=ra_deg, dec_deg=dec_deg, radius_arcsec=radius_arcsec,
              vizier_catalog=vizier_catalog, verbose=verbose)
    except Exception:
        proceed = False
    # pass:try
    #
    if (proceed):
        print('\n#VSX:')
        print('#index raj2000 dej2000 sep_arcsec')
        for j in range(raj2000.size):
            print(j, raj2000[j], dej2000[j], sep_arcsec[j])
        # pass:for
    else:
        print('\nSearch catalog is empty!  (as expected)\n')
    # pass:if
# pass:def


# =============================================================================


if (__name__ == '__main__'):
    xmkpy3_vizier_catalog_cone_get_v4()
# pass:if

# EOF
