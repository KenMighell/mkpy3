#!/usr/bin/env python3

# Kenneth John Mighell
# Kepler Support Scientist
# Kepler / K2 Science Office
# NASA Ames Research Center / SETI Institute


def mkpy3_tess_tpf_overlay_v2():
    import matplotlib.pyplot as plt
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import os
    import ntpath
    import numpy as np

    import lightkurve as lk
    lk.log.setLevel('INFO')
    print(lk.__version__,'=lk.__version__')
    
    import mkpy3_finder_chart_survey_fits_image_get_v1 as km1
    import mkpy3_finder_chart_image_show_v1 as km2
    import mkpy3_finder_chart_tpf_overlay_v5 as km3
    import mkpy3_vizier_vsx_cone_get_v1 as km4
    import mkpy3_vizier_gaia_dr2_cone_get_v1 as km5

    # XZ Cyg (RR Lyrae type variable) [GAIA DR2 2142052889490819328]
    search_result = lk.search_tesscut('XZ Cyg',sector=14)
    print(search_result,'\n^--- search_result')
    tpf = search_result.download(cutout_size=10,quality_bitmask=0)
    print('TPF filename:', ntpath.basename(tpf.path))
    print('TPF dirname: ', os.path.dirname(tpf.path))
    
    target = 'XZ Cyg'
    
    title_ = tpf.hdu[0].header['object']
    title_ += ' : Exoplanet '+target

    ra_deg = tpf.ra
    dec_deg = tpf.dec
    
    title_ = ('%s  [TESS SECTOR=%02d]' % (target,tpf.sector))

    # get survey image data
    width_height_arcmin = 6.00
    survey = '2MASS-J'
    survey_hdu, survey_hdr, survey_data, survey_wcs, survey_cframe = \
      km1.mkpy3_finder_chart_survey_fits_image_get_v1(ra_deg, dec_deg, 
      radius_arcmin=width_height_arcmin, survey=survey)
    
    # create a matplotlib figure object
    fig = plt.figure(figsize=(12,12));
    
    # create a matplotlib axis object with right ascension and declination axes
    ax = plt.subplot(projection=survey_wcs)
            
    # show the survey image
    percentile = 99.5
    km2.mkpy3_finder_chart_image_show_v1(ax=ax, image_data=survey_data, \
      percentile=percentile, verbose=True)
    
    # replace title with a custom title
    plt.suptitle(title_,size=24)
    
    # show the TPF overlay                                                                                                                            
    frame = 0
    km3.mkpy3_finder_chart_tpf_overlay_v5(ax=ax, survey_wcs=survey_wcs, \
      tpf=tpf, frame=frame, verbose=True)
    
    # put a yellow circle at the target position
    ax.scatter(ra_deg*u.deg, dec_deg*u.deg, \
      transform=ax.get_transform(survey_cframe), \
      s=600, edgecolor='yellow', facecolor='None', lw=3, zorder=10);
    
    #import mkpy3_vizier_vsx_cone_get_v1 as km4
    raj2000, dej2000, sep_arcsec, vsx_result = \
      km4.mkpy3_vizier_vsx_cone_get_v1(\
      ra_deg=ra_deg, dec_deg=dec_deg, radius_arcsec=300)
            
    name = np.array(vsx_result['Name'],dtype=np.str)
    print()
    print('\n#VSX:')
    print('#index raj2000 dej2000 sep_arcsec name')
    for j in range(raj2000.size):
        print(j,raj2000[j],dej2000[j],sep_arcsec[j],"'"+name[j]+"'")
    pass#for

    mark_vsx_variables = True
    if (mark_vsx_variables):
        # plot VSX stars with a green X
        ax.scatter(raj2000, dej2000, transform=ax.get_transform(survey_cframe),\
          s=600, color='red', marker='x', lw=3, zorder=20)
    pass#if

    raj2000, dej2000, sep_arcsec, gaia_dr2_result = \
      km5.mkpy3_vizier_gaia_dr2_cone_get_v1(\
      ra_deg=ra_deg, dec_deg=dec_deg, radius_arcsec=180)
    
    mark_gaia_dr2 = True
    if (mark_gaia_dr2):
        name = np.array(gaia_dr2_result['DR2Name'],dtype=np.str)
        print()
        print('\n#GAIA DR2:')
        print('#index raj2000 dej2000 sep_arcsec, DR2Name')
        for j in range(raj2000.size):
            print('%3d %14.9f %14.9f %8.3f %s' % \
              (j,raj2000[j],dej2000[j],sep_arcsec[j],"'"+name[j]+"'"))
            extra = False
            if (extra):
                sc1 = SkyCoord(ra=raj2000[j], dec=dej2000[j], frame=survey_cframe, unit='degree')
                ra_ = sc1.ra.to_string(u.hour)
                dec_ = sc1.dec
                print('%3d %s %s %8.3f' % (j,ra_,dec_,sep_arcsec[j]))
            pass#if
        pass#for
    pass#if

    #ax.scatter(raj2000, dej2000, transform=ax.get_transform(survey_cframe),\
    #      s=600, color='cyan', marker='+', lw=3, zorder=20)
    ax.scatter(raj2000, dej2000, transform=ax.get_transform(survey_cframe),\
      s=300, edgecolor='cyan', facecolor='None', lw=3, zorder=10);

    # reset plotting area to show only the survey image range in X and Y
    nx = survey_data.shape[1]
    ny = survey_data.shape[0]
    ax.set_xlim(0,nx)
    ax.set_ylim(0,ny)
    
    pname = 'mkpy3_plot.png'
    if (pname != ''): 
        plt.savefig(pname, bbox_inches = "tight")
        print(pname,' <--- plot filename has been written!  :-)\n')
    pass#if
pass#def


if (__name__ == '__main__'):
    mkpy3_tess_tpf_overlay_v2()
pass#if

