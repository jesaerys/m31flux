m31flux
=======

Synthetic broadband flux maps of M31 from resolved optical photometry.

(Jacob E. Simones, Julianne J. Dalcanton, Andrew E. Dolphin, Alexia R.
Lewis, Evan D. Skillman, Daniel R. Weisz, Benjamin F. Williams)

This package contains all of the code used for the M31 flux maps project.
Details regarding the data processing and analysis steps are documented in
the individual modules.

`m31flux` requires these packages:

- `astrogrid <http://github.com/jesaerys/astrogrid>`_
- `astropy <http://www.astropy.org>`_
- `match-wrapper <http://github.com/jesaerys/match-wrapper>`_
- `montage-wrapper <http://www.astropy.org/montage-wrapper>`_
- `numpy <http://www.numpy.org>`_


SFH data
========
Alexia has divided all of the PHAT bricks (except 1 and 3, which cover the
bulge) into 15x30 grids and measured the SFH in each grid cell (a.k.a.,
"pixel" or "region"). The SFHs were measured from the PHAT F475W,F814W
".gst" photometry using MATCH assuming the following:

- Kroupa IMF
- Padova isochrones (the default models in MATCH 2.5)
- A distance modulus of 24.47 (McConnachie et al. 2005)
- Metallicity ranging from -2.3 to 0.1, constrained to increase over time
  (the 'zinc' flag in the calcsfh routine in MATCH)
- 34 log-spaced age bins from log10(age) = 6.60 to 10.15.
- A CMD exclusion area to mask out evolved stars from 1.25 to 5.0 in
  F457W-F814W and from 21.0 to 27.2 in F475W.
- The Av,dAv extinction model, where the parameters were optimized for each
  cell.

The details as reported by Alexia:

  "

  I originally ran a grid of SFHs in Av, dAv space with steps of 0.3,
  between A[0,1] and Av+dAv<=2.5. I then ran a routine to interpolate that
  grid so that I could find rough 1-, 2-, and 3-sigma contours of the fit
  values. I then ran all fits that fell within the 2-sigma contours in
  steps of 0.1. I continued this process until there were no new runs. I
  also ran all Av,dAv pairs in the 8-point box immediately surrounding the
  best fit if it wasn't already complete. These were all done without use
  of the -mcdata flag. ::

    $ calcsfh master.sfhparam_0.5 M31-B05_15x30-126.gst.match \
        M31-B05_15x30-126_gst.matchfake M31-B05_15x30-126_0.5-1.9.sfh \
        -zinc -dAvy=0 -dAv=1.9 > M31-B05_15x30-126_0.5-1.9.scrn

  When I completed this process, I determined the best fit and reran that
  with the -mcdata flag. ::

    $ calcsfh master.sfhparam_0.5 M31-B09_15x30-001.gst.match \
        M31-B09_15x30-001_gst.matchfake M31-B09_15x30-001_0.2-2.1_best.sfh \
        -zinc -mcdata -dAvy=0 -dAv=2.1 > M31-B09_15x30-001_0.2-2.1_best.scrn

  I then ran hybridMC::

    $ hybridMC M31-B09_15x30-001_0.2-2.1_best.sfh.dat \
        M31-B09_15x30-001_0.2-2.1_best.mcmc -tint=2.0 -nmc=10000 -dt=0.015

  I originally zcombined just the best fit files and the hmc results and
  then merged them::

    $ zcombine -bestonly M31-B09_15x30-001_0.2-2.1_best.sfh \
        > M31-B09_15x30-001_0.2-2.1_best.zc.fit_unbinned.sfh

    $ zcombine -unweighted -medbest -jeffreys \
        -best=M31-B09_15x30-001_0.2-2.1_best.zc.fit_unbinned.sfh \
        M31-B09_15x30-001_0.2-2.1_best.mcmc \
        > M31-B09_15x30-001_0.2-2.1_best.zc.fit_unbinned.mcmc.sfh

    $ zcmerge M31-B09_15x30-001_0.2-2.1_best.zc.fit_unbinned.sfh \
        M31-B09_15x30-001_0.2-2.1_best.zc.fit_unbinned.mcmc.sfh -absolute \
        > M31-B09_15x30-001_0.2-2.1_best.zc.fit_unbinned.total.sfh

  More recently, after consulting with Andy, I also added in the
  uncertainties due to the search in Av and dAv. To do this, I combined all
  SFHs computed for a given region into one file::

    $ cat SFH1 SFH2 ... SFHn > allfits.sfh

  Then I used zcombine and merged with the hmc results::

    $ zcombine M31-B09_15x30-001_allfits.sfh \
        > M31-B09_15x30-001_0.2-2.1_best.zc.fit_unbinned.best_dust.sfh

    $ zcmerge M31-B09_15x30-001_0.2-2.1_best.zc.fit_unbinned.best_dust.sfh \
        M31-B09_15x30-001_0.2-2.1_best.zc.fit_unbinned.mcmc.sfh -absolute \
        > M31-B09_15x30-001.zc.fit_unbinned.final.sfh

  "

I received the following files from Alexia for each brick:

- List of RA,dec coordinates of the corners of all cells in the brick (file
  kind [1]_ 'corners')
- List of best-fit Av,dAv extinction parameters of all cells in the brick
  (file kind 'extpar')
- For each cell,

  - F475W,F814W ".gst" photometry file (kind 'phot'), and the SFH and CMD
    calcsfh output files (kinds 'sfh' and 'cmd') from the best-fit SFH
    solution

.. [1] "kind" refers to the different types of files tracked by
   `m31flux.config.path`.


Process raw calcsfh data
========================
**process_sfh.py:** The 'sfh' files were processed with zcombine to create
tables of SFR versus age ('bestzcb').


Create synthetic flux maps
==========================
**make_maps.py:** The SFH of each cell was used to calculate its flux in a
given filter. The flux values were assembled into modeled flux maps for
each brick, and the brick maps were added together using Montage to produce
a modeled flux map of M31 over the PHAT survey area.

The following bands were considered:

- GALEX FUV ('mod_fuv_int' for intrinsic flux, 'mod_fuv_red' for reddened
  flux based on the best-fit extinction parameters) and NUV ('mod_nuv_int'
  for intrinsic, 'mod_nuv_red' for reddened)

Images from the GALEX Deep Imaging Survey (DIS) were used to produce maps
of observed FUV ('galex_fuv') and NUV ('galex_nuv') flux.


Issues
======
- An exclude gate was used to mask out the RGB while fitting the CMDs. Has
  anyone assessed the maximum reliable age limit this places on the SFHs?
  A quick method would be to find the oldest isochrone that is visible
  outside the exclude gate.

- Create Spitzer 24um images, too?
