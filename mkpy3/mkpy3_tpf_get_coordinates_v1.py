#!/usr/bin/env python3

# file://mkpy3_tpf_get_coordinates_v1.py

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute

# =============================================================================


def mkpy3_tpf_get_coordinates_v1(
  tpf=None,
  cadence='all',
  recreate_bug=False):
    """
Function: mkpy3_tpf_get_coordinates_v1()

Purpose:

    Returns two 3D arrays of RA and Dec values in decimal degrees.

    If cadence number is given, returns 2D arrays for that cadence (frame). If
    cadence is 'all' returns one RA, Dec value for each pixel in every cadence.

    Uses the FITS WorldCoordinateSsystem solution and the POS_CORR data from
    the TPF header.

Parameters
----------
tpf : lightkurve TargetPixelFile (TPF)
    If None, use kplr007603200-2011271113734_lpd-targ.fits
cadence : 'all' or int
    Which cadences (frames) to return the RA , Dec coordinates for.

Returns
-------
ra : numpy array, same shape as tpf.flux[cadence]
    Array containing RA values at every pixel location
dec : numpy array, same shape as tpf.flux[cadence]
    Array containing Dec values at pixel location.'

# Kenneth Mighell
# Kepler Support Scientist
# NASA Ames Research Center / SETI Institute

ORIGINAL CODE: BEGIN
Call example:
    ra, dec = tpf.get_coordinates()
Source:
    package: lightkurve
    version: 2.0.dev
    file:    lightkurve/lightkurve/targetpixelfile.py
    class:   TargetPixelFile(object)
    method:  get_coordinates(self, cadence='all'):
ORIGINAL CODE: END
    """
    import numpy as np
    import warnings
    import ntpath
    import os
    import sys
    import lightkurve as lk
    #
    if (tpf is None):
        # Exoplanet Kelper-138b is "KIC 7603200":
        tpf = lk.search_targetpixelfile(
            target='kepler-138b', mission='kepler',
            quarter=10).download(quality_bitmask=0)
        print('TPF filename:', ntpath.basename(tpf.path))
        print('TPF dirname: ', os.path.dirname(tpf.path))
    # pass:if
    if (cadence is None):
        cadence = 0  # frame=0
    # pass:if
    #
    self = tpf  # alias
    w = self.wcs
    X, Y = np.meshgrid(np.arange(self.shape[2]), np.arange(self.shape[1]))
    pos_corr1_pix = np.copy(self.hdu[1].data['POS_CORR1'])
    pos_corr2_pix = np.copy(self.hdu[1].data['POS_CORR2'])

    # We zero POS_CORR* when the values are NaN or make no sense (>50px)
    with warnings.catch_warnings():  # Comparing NaNs to numbers is OK here
        warnings.simplefilter("ignore", RuntimeWarning)
        bad = np.any(
            [~np.isfinite(pos_corr1_pix),
             ~np.isfinite(pos_corr2_pix),
             np.abs(pos_corr1_pix - np.nanmedian(pos_corr1_pix)) > 50,
             np.abs(pos_corr2_pix - np.nanmedian(pos_corr2_pix)) > 50], axis=0)
    pos_corr1_pix[bad], pos_corr2_pix[bad] = 0, 0

    # Add in POSCORRs
    X = (np.atleast_3d(X).transpose([2, 0, 1]) +
         np.atleast_3d(pos_corr1_pix).transpose([1, 2, 0]))
    Y = (np.atleast_3d(Y).transpose([2, 0, 1]) +
         np.atleast_3d(pos_corr2_pix).transpose([1, 2, 0]))

    # Pass through WCS
    origin = 0
    if (recreate_bug):
        origin = 1
        sys.stdout.flush()
        print()
        print('***********************************')
        print('* mkpy3_tpf_get_coordinates_v1(): *')
        print('* WARNING: recreate_bug=True      *')
        print('* BAD RA & DEC VALUES COMPUTED!   *')
        print('***********************************')
        print()
    # pass:if
    ra, dec = w.wcs_pix2world(X.ravel(), Y.ravel(), origin)
    ra = ra.reshape((pos_corr1_pix.shape[0], self.shape[1], self.shape[2]))
    dec = dec.reshape((pos_corr2_pix.shape[0], self.shape[1], self.shape[2]))
    ra, dec = ra[self.quality_mask], dec[self.quality_mask]
    if cadence != 'all':
        return ra[cadence], dec[cadence]
    # pass:if
    return ra, dec
# pass:if


# =============================================================================


def xmkpy3_tpf_get_coordinates_v1():
    """
Unit test
    """
    import numpy as np
    import lightkurve as lk
    #
    print(lk.__version__, '=lk.__version__')

    def msg(ok, tag_):
        print('***' + tag_ + ': ', end="")
        if (ok):
            print('PASS***')
        else:
            print("FAIL***")
    # pass:def

    tpf = lk.search_targetpixelfile(
        target='kepler-138b', mission='kepler', cadence='long',
        quarter=10).download(quality_bitmask=0)

    w = tpf.wcs  # alias

    ll_x0 = 0
    ll_y0 = 0
    print(ll_x0, '=ll_x0')
    print(ll_y0, '=ll_y0')

    origin0 = 0
    ra_ll_x0, dec_ll_y0 = w.wcs_pix2world(ll_x0, ll_y0, origin0)
    print(ra_ll_x0, dec_ll_y0, '=ra_ll_x0, dec_ll_y0')

    print()
    x0_ra_ll_x0, y0_dec_ll_y0 = w.wcs_world2pix(ra_ll_x0, dec_ll_y0, origin0)
    print(x0_ra_ll_x0, y0_dec_ll_y0,
          '=x0_ra_ll_x0, y0_dec_ll_y0 [should be about (0,0)]')

    ra_x0_ra_ll_x0, dec_y0_dec_ll_y0 = w.wcs_pix2world(
        x0_ra_ll_x0, y0_dec_ll_y0, origin0)
    print(ra_x0_ra_ll_x0, dec_y0_dec_ll_y0,
          '=ra_x0_ra_ll_x0, dec_y0_dec_ll_y0')

    print('\nra_x0_ra_ll_x0 is_close_to ra_ll_x0 ?')
    ok = np.abs(ra_x0_ra_ll_x0 - ra_ll_x0) < 0.000001
    msg(ok, 'TEST1')
    print('^--- THIS BETTER PASS!')

    print('\ndec_y0_dec_ll_y0 is_close_to dec_ll_y0 ?')
    ok = np.abs(dec_y0_dec_ll_y0 - dec_ll_y0) < 0.000001
    msg(ok, 'TEST2')
    print('^--- THIS BETTER PASS!')

    print()
    frame0 = 0
    # Set one of the next 3 if statements to TRUE depending on the function to
    # be tested
    if (False):
        print('---> check tpf.get_coordinates()')
        rax_ll_x0 = tpf.get_coordinates()[0][frame0][0][0]
        decx_ll_y0 = tpf.get_coordinates()[1][frame0][0][0]
        print('NOTE: next two tests will PASS --- if the tpf.get_coordinates b'
              'ug has been fixed')
    # pass:if
    if (True):
        print('---> check mkpy3_tpf_get_coordinates_v1()')
        rax_ll_x0 = mkpy3_tpf_get_coordinates_v1(tpf=tpf)[0][frame0][0][0]
        decx_ll_y0 = mkpy3_tpf_get_coordinates_v1(tpf=tpf)[1][frame0][0][0]
        print('NOTE: next two tests should PASS')
    # pass:if
    if (False):
        print('---> check mkpy3_tpf_get_coordinates_v1(...,recreate_bug=True)')
        rax_ll_x0 = mkpy3_tpf_get_coordinates_v1(
            tpf=tpf, recreate_bug=True)[0][frame0][0][0]
        decx_ll_y0 = mkpy3_tpf_get_coordinates_v1(
            tpf=tpf, recreate_bug=True)[1][frame0][0][0]
        print('NOTE: next two tests should FAIL')
    # pass:if

    print(rax_ll_x0, decx_ll_y0, '=rax_ll_x0, decx_ll_y0')

    print()
    x0_rax_ll_x0, y0_decx_ll_y0 = \
        w.wcs_world2pix(rax_ll_x0, decx_ll_y0, origin0)
    print(x0_rax_ll_x0, y0_decx_ll_y0, '=x0_rax_ll_x0, y_decx_ll_y0')

    tpf_pos_corr1_frame0 = tpf.pos_corr1[frame0]
    tpf_pos_corr2_frame0 = tpf.pos_corr2[frame0]
    print(tpf_pos_corr1_frame0, tpf_pos_corr2_frame0,
          '=tpf_pos_corr1_frame0, tpf_pos_corr2_frame0')

    xx0_rax_ll_x0 = x0_rax_ll_x0 - tpf_pos_corr1_frame0
    yy0_decx_ll_y0 = y0_decx_ll_y0 - tpf_pos_corr2_frame0
    print(xx0_rax_ll_x0, yy0_decx_ll_y0,
          '=xx0_rax_ll_x0, yy0_decx_ll_y0 [should be about (0,0)]')

    ra_xx0_rax_ll_x0, dec_yy0_decx_ll_y0 = \
        w.wcs_pix2world(xx0_rax_ll_x0, yy0_decx_ll_y0, origin0)
    print(ra_xx0_rax_ll_x0, dec_yy0_decx_ll_y0,
          '=ra_xx0_rax_ll_x0, dec_yy0_decx_ll_y0')

    print('\nra_xx0_rax_ll_x0 is_close_to ra_ll_x0 ?')
    ok = np.abs(ra_xx0_rax_ll_x0 - ra_ll_x0) < 0.000001
    msg(ok, 'TEST3')

    print('\ndec_yy0_decx_ll_y0 is_close_to dec_ll_y0 ?')
    ok = np.abs(dec_yy0_decx_ll_y0 - dec_ll_y0) < 0.000001
    msg(ok, 'TEST4')
# pass:def


# =============================================================================


if (__name__ == '__main__'):
    xmkpy3_tpf_get_coordinates_v1()
# pass:if

# EOF
