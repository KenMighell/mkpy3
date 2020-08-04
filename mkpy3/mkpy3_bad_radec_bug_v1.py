#!/usr/bin/env python3

# Kenneth John Mighell
# Kepler/K2 Support Scientist
# Kepler/K2 Science Office
# NASA Ames Research Center / SETI Institute


def mkpy3_bad_radec_bug_v1(tpf=None, verbose=None):
    """
Function : mkpy3_bad_radec_bug_v1()

Purpose: Check to see if the tpf.get_coordinates radec bug is present

Parameters
----------
tpf : lightkurve TargetPixelFile (TPF) object (optional)
    if None, use kplr007603200-2011271113734_lpd-targ.fits
verbose : bool (optional)
    if True, print extra information

Returns 
-------
bad_radec : bool
    if True, the tpf.get_coordinates radec bug exists

# Kenneth John Mighell
# Kepler/K2 Support Scientist
# Kepler/K2 Science Office
# NASA Ames Research Center / SETI Institute
    """
    import numpy as np
    import os
    import ntpath
    import lightkurve as lk
    if (verbose is None): verbose = False
    if (tpf is None):
        tpf = lk.search_targetpixelfile(target='kepler-138b',mission='kepler',\
          quarter=10).download(quality_bitmask=0)  # Exoplanet Kelper-138b is "KIC 7603200" 
    pass#if
    if (verbose): 
        print(lk.__version__,'=lk.__version__') 
        print('TPF filename:', ntpath.basename(tpf.path))
        print('TPF dirname: ', os.path.dirname(tpf.path))
    pass#if
    ra, dec = tpf.wcs.wcs_pix2world(0.0, 0.0, 0)
    if (verbose): print(ra, dec, '=ra, dec')  # use ds9 to verify ra and dev values
    rax1, decx1 = tpf.wcs.wcs_pix2world(tpf.pos_corr1[0], tpf.pos_corr2[0], 0) # first (zeroth) frame
    if (verbose): print(rax1, decx1, '=rax1, decx1')
    x1, y1 = tpf.wcs.wcs_world2pix(rax1, decx1,0)
    if (verbose): print(x1, y1, '=x1, y1')
    rax2 = tpf.get_coordinates()[0][0][0][0]  # first (zeroth) frame
    decx2 = tpf.get_coordinates()[1][0][0][0] # first (zeroth) frame
    if (verbose): print(rax2, decx2, '=rax2, decx2')
    x2, y2 = tpf.wcs.wcs_world2pix(rax2, decx2,0)
    if (verbose): print(x2, y2, '=x2, y2')
    bad_ra = (np.abs(rax2-rax1) > 0.00001)
    if (verbose): print(bad_ra, '=bad_ra')
    bad_dec = (np.abs(decx2-decx1) > 0.00001)
    if (verbose): print(bad_dec, '=bad_dec')
    bad_radec = bad_ra or bad_dec
    if (verbose): print(bad_radec,'=bad_radec')
    return bad_radec # if true, radec pairs returned by tpf.get_cooridnates() are not correct
pass#def

if (__name__ == '__main__'):
    bug_found  = mkpy3_bad_radec_bug_v1(verbose=True)
    print(flush=True)
    if (bug_found):
        print('*** BUG FOUND ***  8=X')
    else:
        print('OK')
    pass#if
pass#if
#EOF