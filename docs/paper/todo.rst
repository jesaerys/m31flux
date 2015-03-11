Bugs in the code
================
* fluxmap.py: make_mod_fuv_int() and make_mod_nuv_int() both have the kwargs
  ``agelimdmod=config.DIST.distmod``.

  - obviously ``agelimdmod`` is not supposed to be an argument.
  - is ``agelim`` supposed to be set or not? It only shows up in
    make_mod_xuv_int(), not in make_mod_xuv_red().
  - Note that dmod is set within util.calc_flux, so providing dmod as a kwarg
    doesn't current do anything.

  * This is just a cosmetic bug, the results are not erroneous. The modeled
    fluxes have been derived using the full SFHs and the correct dmod.


Metallicity evolution
=====================
Is it feasible to tackle this? Probably won't make much difference.


Observational uncertainties
===========================
For statistical/Poisson uncertainties, divide the \fuv{} and \nuv{} intensity
(count) maps by the CCD gain to convert the units into photons, and take the
square root. Flux is proportional to photon rate, so the flux uncertainties are
just a constant times the raw Poisson uncertainties. The uncertainties won't be
useful for the map visualizations, but they will be important for scatter plots
comparing the synthetic and observed fluxes.


Background subtraction
======================
ACTION:
    Don't worry about this. What I've done is pretty standard. Just highlight
    that one of the shortcomings is its influence on faint pixels.

What I've done:
    The background subtraction was based on a measurement of the mean flux level
    in a nearby, off-galaxy area of the sky. The measured FUV background value
    was insufficient, causing an obvious downturn in the flux ratio distribution
    (think f/(f+bg) vs. f for f <~ bg). The opposite effect was found for the
    measured NUV background. In each filter, I compensated for this by adjusting
    the background value until the flux ratio distribution was "flat", i.e.,
    until the median flux ratio was no longer a function of observed flux.

What I would like to do:
    Actually fit the background, e.g., using a chi-square minimization of f_sfh
    - f_obs. Although this method is certainly more justifiable than what I've
    done, I don't expect that the background subtraction will ultimately make
    much difference for the overall results of the paper. However, I will still
    try to implement the fitting method as time allows.


NUV timescale
=============
ACTION:
    Respond to everyone, clarifying that both the FUV and the NUV flux
    calibrations under consideration assume a constant SFR over the last 100
    Myr. This is mentioned in the second to last paragraph in Section 4. Using
    200 Myr for NUV would ruin the comparison.

Andy:
    Figure 11 that SFR(NUV) should be compared with the mean SFR over the past 200
    Myr, not that over the past 100 Myr.

Evan:
    The comment about averaging over 200 Myr for the NUV seems appropriate.
