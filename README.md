# mkpy3
mkpy3 : Python3 tools for working with data from NASA's TESS/K2/Kepler missions.

Kenneth Mighell  
Kepler Support Scientist   
NASA Ames Research Center / SETI Institute

---

#### Making a finder chart for 2-min TESS observations of CD Ind using mkpy3

<!-- original: works with GitHub
<p float="left">
  <img src="./mkpy3_plot_figa.png" width="320" />
  <img src="./mkpy3_plot_figb.png" width="240" /> 
</p>
-->

<!-- revised: works with GitHub *and* test.pypi.org -->
<p float="left">
  <img src="https://github.com/KenMighell/mkpy3/raw/master/mkpy3_plot_figa.png" width="320" />
  <img src="https://github.com/KenMighell/mkpy3/raw/master/mkpy3_plot_figb.png" width="240" /> 
</p>


* *blue squares* show the TESS target pixel file overlay.
* *red squares* show the TESS target piexl file aperture overlay.
* *green lines* show the Right Ascension / Declination grid.
* *yellow circle* marks the target (CD Ind).
* *cyan circles* show some of the GAIA DR2 catalog stars in the field.
* *green X* shows the only VSX catalog star in the field (the target).
* *Compass rose*:
    * *long arm* of the compass rose points NORTH.
    * *short arm* of the compass rose points EAST.

**Click the following link to see a Jupyter notebook that creates the above two plots:**

[https://github.com/KenMighell/mkpy3/blob/master/mkpy3/docs/source/tutorials/nb_tess_finder_chart_cd_ind.ipynb](https://github.com/KenMighell/mkpy3/blob/master/mkpy3/docs/source/tutorials/nb_tess_finder_chart_cd_ind.ipynb)

Alternatively, you can create these plots using the following Python 3 code snippets.

**The following code snippet creates the finder chart (above left) using mkpy3:**

```
# download TESS 2-min observation of CD Ind from STScI
import lightkurve as lk
search_results = lk.search_targetpixelfile('CD Ind', radius=60, 
  mission='TESS', sector=1)    
tpf = search_results[0].download(quality_bitmask=0)
#
# show the TPF overlay on rotated "DSS2 Red" survey image
import mkpy3_tpf_overlay_v6 as km1
import mkpy3_plot_add_compass_rose_v5 as km2
ax = km1.mkpy3_tpf_overlay_v6(tpf=tpf, survey='DSS2 Red', 
  rotationAngle_deg='tpf', width_height_arcmin=6, percentile=99.5,
  shrink=0.4, show_plot=False, plot_file='',
  title='CD Ind : TESS : Sector 1', print_gaia_dr2=False)
ax.coords[0].set_major_formatter('d.dd')
ax.coords[1].set_major_formatter('d.dd')
ax.tick_params(axis='x', labelsize=16, length=5, width=2,
  labeltop=True, labelbottom=True)
ax.tick_params(axis='y', labelsize=16, length=5, width=2,
  labelright=True, labelleft=True)
ax.grid(True, color='palegreen', lw=2, zorder=1)  # show RA/DEC grid
km2.mkpy3_plot_add_compass_rose_v5(ax=ax, north_arm_arcsec=50)
#
# save, show, and close the plot
import matplotlib.pyplot as plt
plt.savefig('mkpy3_plot1.png', bbox_inches="tight");
plt.show(); plt.close()
```

**The following code snippet creates the annotated TPF image plot (above right) using mkpy3:**

```
# download TESS 2-min observation of CD Ind from STScI
import lightkurve as lk
search_results = lk.search_targetpixelfile('CD Ind', radius=60,  
  mission='TESS', sector=1)    
tpf = search_results[0].download(quality_bitmask=0)
#
# show the first frame of the TPF with RA/DEC axis and grid
import matplotlib.pyplot as plt
import astropy.visualization as av
import mkpy3_plot_add_compass_rose_v5 as km1
fig = plt.figure(figsize=(7, 7))
ax = plt.subplot(projection=tpf.wcs)
interval = av.PercentileInterval(99.9)
stretch = av.SqrtStretch()
frame = 0
image_data = tpf.flux[frame]
norm = av.ImageNormalize(image_data, interval=interval,  
  stretch=stretch)
ax.imshow(image_data, norm=norm, cmap='gray_r')
ax.set_xlabel('Right Ascension (J2000)', size=24)
ax.set_ylabel('Declination (J2000)', size=24)
ax.tick_params(axis='x', labelsize=16, length=5, width=2,
  labeltop=True, labelbottom=True)
ax.tick_params(axis='y', labelsize=16, length=5, width=2,
  labelright=True, labelleft=True)
ax.coords[0].set_major_formatter('d.dd')
ax.coords[1].set_major_formatter('d.dd')
ax.grid(True, color='palegreen', lw=2, zorder=1)
km1.mkpy3_plot_add_compass_rose_v5(ax=ax,  # add a compass rose
  north_arm_arcsec=50, wcs=tpf.wcs)
marker_kwargs = {'edgecolor': 'yellow',  
  's': 600, 'facecolor': 'None', 'lw': 3, 'zorder': 10}
ax.scatter(tpf.ra, tpf.dec, transform=ax.get_transform('icrs'),  
  **marker_kwargs)
fig.suptitle('CD Ind : TESS : Sector 1 : Frame 0', size=24)
#
# Save, show, and close the plot
plt.savefig('mkpy3_plot2.png', bbox_inches="tight");  
plt.show(); plt.close()
```

[//]: # (EOF)
