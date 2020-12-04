#!/usr/bin/env python3

# file://mkpy3_vizier_gaia_dr2_cone_get_v2.py

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute

# =============================================================================


def mkpy3_vizier_gaia_dr2_cone_get_v2(
  ra_deg=None,
  dec_deg=None,
  radius_arcsec=None,
  verbose=None
):
    """
Function : mkpy3_vizier_gaia_dr2_cone_get_v2()

Purpose: Perform a cone search of the GAIA DR2 catalog using Vizier.

Parameters
----------
ra_deg : float (optional)
    right ascencsion of target[deg]
dec_deg : float (optional)
    declination of target [deg]
radius_arcsec : float (optional)
    search radius [arcsec]
verbose : bool (optional)
    verbose output

Returns
-------
raj2000 : float array
    right ascension (J2000) [deg]
dej2000 : float array
    declination (J2000) [deg]
sep_arcsec : float array
    separation from target [arcsec]
vizier_gaia_dr2_result :
    VSX table returned by Vizier

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute
    """
    import mkpy3
    if (verbose is None):
        verbose = False
    if ((ra_deg is None) or (dec_deg is None)):
        ra_deg = 291.3663013467642  # RR Lyr
        dec_deg = +42.7843585094725  # RR Lyr
        verbose = True
    # pass:if
    if (radius_arcsec is None):
        radius_arcsec = 300
    #
    vizier_catalog = 'I/345/gaia2'
    try:
        raj2000, dej2000, sep_arcsec, vizier_gaia_dr2_result = \
            mkpy3.mkpy3_vizier_catalog_cone_get_v4(
              ra_deg=ra_deg,
              dec_deg=dec_deg,
              radius_arcsec=radius_arcsec,
              vizier_catalog=vizier_catalog,
              verbose=verbose)
    except Exception:
        raj2000 = None
        dej2000 = None
        sep_arcsec = None
        vizier_gaia_dr2_result = None
    # pass:try
    #
    return (raj2000, dej2000, sep_arcsec, vizier_gaia_dr2_result)
# pass:def


# =============================================================================


def xmkpy3_vizier_gaia_dr2_cone_get_v2():
    import numpy as np
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    #
    ra_deg = 291.3663013467642   # RR Ly
    dec_deg = +42.7843585094725  # RR Lyr
    radius_arcsec = 20
    raj2000, dej2000, sep_arcsec, gaia_dr2_result = \
        mkpy3_vizier_gaia_dr2_cone_get_v2(
          ra_deg=ra_deg, dec_deg=dec_deg, radius_arcsec=radius_arcsec,
          verbose=True)

    assert(gaia_dr2_result is not None)
    name = np.array(gaia_dr2_result['DR2Name'], dtype=np.str)
    print('\n#GAIA2_DR2:')
    print('#index raj2000 dej2000 sep_arcsec')
    for j in range(raj2000.size):
        print(j, raj2000[j], dej2000[j], sep_arcsec[j], "'" + name[j] + "'")
    # pass:for

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
    print('# n GAIA2_Source             sep  RA_ICRS       DE_ICRS        pmRA'
          '     pmDE      Plx')
    print('#                       [arcsec]  [deg]         [deg]      [mas/yr]'
          ' [mas/yr]    [mas]')
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
        print('%3d %d %8.3f %12.7f %12.7f %8.3f %8.3f %8.3f' %
              (kk, srcj, sepj, raj, dej, pmraj, pmdej, plxj))

    print()
    print('# n GAIA2_Source          Gmag    RA_ICRS       DE_ICRS      RA_ICR'
          'S         DE_ICRS')
    print('#                         [mag]   [deg]         [deg]        [hms] '
          '          [dms]')
    for k in range(len(xra)):
        j = k  # idx[k]
        xraj = xra[j]
        ydej = yde[j]
        gmagj = gmag[j]
        sc1 = SkyCoord(ra=xraj, dec=ydej, frame='icrs', unit='degree')
        ra_ = sc1.ra.to_string(u.hour)
        dec_ = sc1.dec
        sepj = sep_arcsec[j]
        kk = k + 1
        srcj = src[j]
        print('%3d %d %8.3f %12.7f %12.7f  %15s %15s' %
              (kk, srcj, gmagj, xraj, ydej, ra_, dec_))

    print()
    print(ra_deg, dec_deg, '=ra_deg, dec_deg')

    # empty catalog test
    ra_deg = 0
    dec_deg = 0
    radius_arcsec = 0.1
    raj2000, dej2000, sep_arcsec, gaia_dr2_result = \
        mkpy3_vizier_gaia_dr2_cone_get_v2(
          ra_deg=ra_deg, dec_deg=dec_deg, radius_arcsec=radius_arcsec,
          verbose=True)
    assert(gaia_dr2_result is None)
    print('\n:-)')
# pass:if


# =============================================================================


if (__name__ == '__main__'):
    xmkpy3_vizier_gaia_dr2_cone_get_v2()
# pass:if

# EOF
