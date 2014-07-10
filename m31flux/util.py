"""

==============
`m31flux.util`
==============

General utilities.


Functions
---------

================== =========================================================
`safe_mkdir`       Create a directory only if it does not already exist.
`calc_flux`        Calculate flux in a given filter for the SFH of the given
                   cell.
`make_brick_image` Apply a function to all cells in given brick, then
                   assign each cell to a pixel to assemble an image.
`galex_pre`        Mask border pixels and convert into flux units.
`galex_post`       Measure the background flux level and subtract it from
                   the image.
================== =========================================================

"""
import astrogrid
import astropy.io.fits
import astropy.stats.funcs
import errno
import match_wrapper as match
import numpy as np
import os

from . import config



# Shell stuff
# -----------

def safe_mkdir(path):
    """Create a directory only if it does not already exist."""
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise



# Maps
# ----

def calc_flux(brick, cell, band, **kwargs):
    """Calculate flux in a given filter for the SFH of the given cell.

    Parameters
    ----------
    brick, cell : int
        Brick and cell numbers.
    band : str
        Filter in which flux is calculated. See `fsps.find_filter` or
        `fsps.fsps.FILTERS` for valid filter names.
    agelim : float, optional
        Maximum age (yr) considered in the SFH; data for older ages are
        ignored in the flux calculation. If None (default), the full SFH is
        used.
    dmod : float, optional
        Distance modulus for calculating apparent magnitude. Default is
        None. See `astrogrid.flux.calc_mag`.
    av, dav : float, optional
        Av and dAv extinction parameters. Default is None. See
        `astrogrid.flux.calc_sed`.
    dust_curve : string, optional
        Dust curve. Should be the name of a key in
        `astrogrid.flux.DUST_TYPE`. Default is 'cardelli'.
    fsps_kwargs : dict, optional
        Dictionary of keyword arguments for `fsps.StellarPopulation`.
        Default is an empty dictionary. See `astrogrid.flux.calc_sed`.


    Returns
    -------
    float
        Flux (in erg s-1 cm-2 A-1) in the given filter modeled from the SFH
        of the given brick and cell.

    """
    agelim = kwargs.get('agelim')
    dmod = kwargs.get('dmod')
    av, dav = kwargs.get('av'), kwargs.get('dav')
    dust_curve = kwargs.get('dust_curve', 'cardelli')

    fsps_kwargs = kwargs.get('fsps_kwargs')
    if fsps_kwargs is None:
        fsps_kwargs = {}

    zcbfile = config.path('bestzcb', field=brick, subfield=cell)

    if zcbfile:
        sfh = match.io.open_zcbfile(zcbfile)
        sfh['age_i'] = 10**sfh['log(age_i)']  # Linearize ages
        sfh['age_f'] = 10**sfh['log(age_f)']
        if agelim is not None:
            i = sfh['age_f'] <= agelim
            sfh = sfh[i]

        # Use the mean metallicity over the last 100 Myr, or the most
        # recent metallicity available before then.
        i = sfh['SFR'] > 0
        j = sfh['log(age_i)'] < 8
        if not np.sum(i & j):
            if not np.sum(i):
                logZ = 0  # If metallcities aren't available for some reason, guess solar
            else:
                logZ = sfh['[M/H]'][i][0]  # Most recent available
        else:
            logZ = np.log10(np.mean(10**sfh['[M/H]'][i & j]))  # 100 Myr mean value
        fsps_kwargs['zmet'] = astrogrid.flux.get_zmet(logZ)

        sfh['SFR'][0] *= 1 - sfh['age_i'][0]/sfh['age_f'][0]  # Rescale 1st age bin
        sfh['age_i'][0] = 0

        age, sfr = (sfh['age_i'], sfh['age_f']), sfh['SFR']
        wave, spec = astrogrid.flux.calc_sed(
            sfr, age, av=av, dav=dav, dust_curve=dust_curve,
            fsps_kwargs=fsps_kwargs)
        mag = astrogrid.flux.calc_mag(wave, spec, band, dmod=dmod)
        flux = astrogrid.flux.mag2flux(mag, band)

    else:
        flux = np.nan

    return flux


def calc_mean_sfr(brick, cell):
    """100 Myr mean SFR."""
    zcbfile = config.path('bestzcb', field=brick, subfield=cell)

    if zcbfile:
        sfh = match.io.open_zcbfile(zcbfile)

        # Get SFH over last 100 Myr
        i = sfh['log(age_f)'] <= 8.00
        sfr = sfh['SFR'][i]
        agei, agef = 10**table['log(age_i)'][i], 10**table['log(age_f)'][i]

        # Calculate mean SFR
        dt = agef - agei
        mean = np.average(sfr, weights=dt)  # Equivalent to np.sum(mass)/1e8

    else:
        mean = np.nan

    return mean


def make_brick_image(brick, func):
    """Apply a function to all cells in given brick, then assign each cell
    to a pixel to assemble an image.

    Parameters
    ----------
    brick : int
        The brick number for which to create an image.
    func : function
        Any function that takes a brick number and a cell number as the
        first and second arguments and returns the value for the cell.

    Returns
    -------
    astropy.io.fits.PrimaryHDU
        The image data and header.

    """
    # Make the grid
    cell_list = config.CELL_LIST[config.SORT]  # Arrange cells in grid order
    grid = astrogrid.Grid((config.NROW, config.NCOL), func,
                          [(brick, cell) for cell in cell_list],
                          update=True)

    # Create a header
    coordfile = config.path('corners', field=brick)
    lon, lat = config.cornergrid(config.open_cornersfile(coordfile))
    x, y = grid.edges
    hdr = astrogrid.wcs.make_header(x, y, lon, lat)

    hdu = astropy.io.fits.PrimaryHDU(grid.data_grid, header=hdr)

    return hdu


def galex_pre(data, hdr, band):
    """Mask border pixels and convert into flux units.

    The chip is not perfectly circular, so a small number of pixels near
    the border that are legitimately zero may be masked.

    Parameters
    ----------
    data : array
        Image data.
    hdr : astropy.io.fits.Header
        FITS header for the image.
    band : {'galex_fuv', 'galex_nuv'}
        The GALEX filter of the data. This is used for converting count
        rates into fluxes.

    Returns
    -------
    tuple
        The processed data array and header.

    """
    # Set border pixels to np.nan
    y, x = np.indices(data.shape)  # Array coordinates
    x, y = x + 1, y + 1  # Pixel coordinates of the pixel centers
    r = np.sqrt((x - config.GALEX_CHIP_X0)**2 + (y - config.GALEX_CHIP_Y0)**2)
    i = (r > config.GALEX_CHIP_RAD) & (data == 0)
    data = np.where(i, np.nan, data)  # A new array (instead of data[i])

    # Convert to flux
    data = astrogrid.flux._galex_cps2flux(data, band)

    return data, hdr


def galex_post(data, hdr, band):
    """Measure the background flux level and subtract it from the image.

    The background (foreground) flux is measured from the mean value of
    pixels within a predefined rectangle, using sigma clipping to reject
    non-background pixels.

    .. note:: The rectangle coordinates and clipping limit are hardcoded.
       The chosen rectangle is a small area apparently devoid of stars in
       GALEX NUV.

    Parameters
    ----------
    data : array
        Image data.
    hdr : astropy.io.fits.Header
        FITS header for the image.
    band : {'galex_fuv', 'galex_nuv'}
        The GALEX filter of the data. This is used for converting count
        rates into fluxes.

    Returns
    -------
    tuple
        The processed data array and header.

    """
    sig = 1  # Sigma clipping limit

    # Rectangle
    x1, x2 = 0, 13
    y1, y2 = 165, 211

    ### maybe this only works for small areas where the background doesn't
    ### have a slope?

    ### could also try min, or bottom x% of the pixels

    # Estimate the foreground level and subtract
    sample = data[y1:y2,x1:x2].ravel()
    #sample = sample[sample < np.median(sample)]  # bottom 50%
    mask = astropy.stats.funcs.sigma_clip(
        sample, sig=sig, iters=None, cenfunc=np.mean).mask
    val = np.mean(sample[mask])
    data = data - val

    return data, hdr



# Deprecated?
# ```````````
def calc_flux_scombine(brick, cell, band, imf_type=2, dmod=None, av=None, dav=None):
    """Deprecated?"""
    # Lookup logZ by zmet parameter
    logZ_dict = {
        16: -0.39,
        17: -0.30,
        18: -0.20,
        19: -0.10,
        20: 0.00,
        }

    zcbfile = config.path('bestzcb', field=brick, subfield=cell)

    if zcbfile:
        sfh = match.io.open_zcbfile(zcbfile)

        # Use the mean metallicity over the last 100 Myr, or the most
        # recent metallicity available before then.
        i = sfh['SFR'] > 0
        j = sfh['log(age_i)'] < 8
        if not np.sum(i & j):
            if not np.sum(i):
                logZ = 0  # If metallcities aren't available for some reason, guess solar
            else:
                logZ = sfh['[M/H]'][i][0]  # Most recent available
        else:
            logZ = np.mean(sfh['[M/H]'][i & j])  # 100 Myr mean value
        zmet = astrogrid.flux.get_zmet(logZ)
        logZ = logZ_dict[zmet]

        astrogrid.flux.calc_mag_scombine(zcbfile, band, SPEC_DIR, imf_type,
                                         logZ, dmod=dmod, av=av, dav=dav)
        flux = astrogrid.flux.mag2flux(mag, band)

    else:
        flux = np.nan

    return flux
