#!/usr/bin/env python3

# file://mkpy3_vizier_vsx_cone_get_v2.py

# Kenneth Mighell
# SETI Institute

# =============================================================================


def mkpy3_vizier_vsx_cone_get_v2(
    ra_deg=None, dec_deg=None, radius_arcsec=None, verbose=None
):
    """
Function : mkpy3_vizier_vsx_cone_get_v1()

Purpose: Perform a cone search of the VSX catalog using Vizier.

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
vizier_vsx_result :
    VSX table returned by Vizier

Kenneth Mighell
SETI Institute
    """
    import mkpy3

    if verbose is None:
        verbose = False
    if (ra_deg is None) or (dec_deg is None):
        ra_deg = 291.3663013467642  # RR Lyr
        dec_deg = +42.7843585094725  # RR Lyr
        verbose = True
        # fi
    if radius_arcsec is None:
        radius_arcsec = 300
    vizier_catalog = "B/vsx/vsx"
    try:
        (
            raj2000,
            dej2000,
            sep_arcsec,
            vizier_vsx_result,
        ) = mkpy3.mkpy3_vizier_catalog_cone_get_v4(
            ra_deg=ra_deg,
            dec_deg=dec_deg,
            radius_arcsec=radius_arcsec,
            vizier_catalog=vizier_catalog,
            verbose=verbose,
        )
    except Exception:
        raj2000 = None
        dej2000 = None
        sep_arcsec = None
        vizier_vsx_result = None
        # yrt
    #
    return (raj2000, dej2000, sep_arcsec, vizier_vsx_result)
    # fed


# =============================================================================


def xmkpy3_vizier_vsx_cone_get_v2():
    import numpy as np

    ra_deg = 291.3663013467642  # RR Lyr
    dec_deg = +42.7843585094725  # RR Lyr
    radius_arcsec = 300
    raj2000, dej2000, sep_arcsec, vsx_result = mkpy3_vizier_vsx_cone_get_v2(
        ra_deg=ra_deg, dec_deg=dec_deg, radius_arcsec=radius_arcsec, verbose=True
    )
    #
    assert vsx_result is not None
    name = np.array(vsx_result["Name"], dtype=str)
    print("\n#VSX:")
    print("#index raj2000 dej2000 sep_arcsec name")
    for j in range(raj2000.size):
        print(j, raj2000[j], dej2000[j], sep_arcsec[j], "'" + name[j] + "'")
    #
    # empty catalog test
    ra_deg = 0.0
    dec_deg = 0.0
    radius_arcsec = 0.1
    raj2000, dej2000, sep_arcsec, vsx_result = mkpy3_vizier_vsx_cone_get_v2(
        ra_deg=ra_deg, dec_deg=dec_deg, radius_arcsec=radius_arcsec, verbose=True
    )
    assert vsx_result is None
    print(":-)")
    return None
    # fed


# =================================================


if __name__ == "__main__":
    xmkpy3_vizier_vsx_cone_get_v2()
    # fi

# EOF
