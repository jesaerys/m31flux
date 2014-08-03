"""

==============
`m31flux.util`
==============

General utilities.


Functions
---------

================== =========================================================
`calc_flux`        Calculate flux in a given filter for the SFH of the
                   given cell.
`calc_mean_sfr`    Calculate the 100 Myr mean SFR for the SFH of the given
                   cell.
`make_brick_image` Apply a function to all cells in given brick, then
                   assign each cell to a pixel and assemble an image.
`make_counter`     Create a generator that wraps an arbitrary iterator and
                   prints a progress counter.
`galex_pre`        Mask border pixels in a GALEX image and convert into
                   flux units.
`galex_post`       Measure and subtract the background flux level for a
                   GALEX mosaic.
================== =========================================================

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import astrogrid
import astropy.io.fits
import match_wrapper as match
import numpy as np
import os
import sys

from . import config


def make_counter(prefix=None, end=None, newline=True):
    """Create a generator that wraps an arbitrary iterator and prints a
    progress counter.

    The returned generator takes an iterator as its only argument and
    updates the progress string every time it yields the next value in the
    iterator. The progress string has the form::

      {prefix}{i}/{total}{end}

    where ``i`` is the iteration number. All output is printed to stdout,
    and the stdout buffer is flushed after each iteration.

    Parameters
    ----------
    prefix : str, optional
        String to print before the counter. Default is None.
    end : str, optional
        String to print after the counter once all iterations have ocurred.
        Default is None.
    newline : bool, optional
        Whether to add a newline to the end of `end`. Default is True.

    Returns
    -------
    generator

    Examples
    --------
    >>> import time
    >>> counter = make_counter(prefix='progress: ', end='  finished')
    >>> for i in counter(range(10)):
    ...     time.sleep(0.5)

    """
    prefix = prefix if prefix else ''
    end = end if end else ''
    end = '{}\n'.format(end) if newline else end
    def counter(iterable):
        total = len(iterable)
        width = len(str(total))
        kwargs = dict(prefix=prefix, width=width, total=total)
        line = '\r{prefix}{0:{width}d}/{total}'
        sys.stdout.write(line.format(0, **kwargs))
        sys.stdout.flush()
        for i, val in enumerate(iterable):
            yield val
            sys.stdout.write(line.format(i+1, **kwargs))
            sys.stdout.flush()
        sys.stdout.write(end)
    return counter


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
    av, dav : float, optional
        Av and dAv extinction parameters. Default is None. See
        `astrogrid.flux.calc_sed`.
    dust_curve : str, optional
        Dust curve. Should be the name of a key in
        `astrogrid.flux.DUST_TYPE`. Default is 'cardelli'.
    fsps_kwargs : dict, optional
        Dictionary of keyword arguments for `fsps.StellarPopulation`.
        Default is an empty dictionary. See `astrogrid.flux.calc_sed`.

    Returns
    -------
    float
        Flux (in erg s-1 cm-2 A-1) or magnitude in the given filter modeled
        from the SFH of the given brick and cell.

    """
    agelim = kwargs.get('agelim')
    av, dav = kwargs.get('av'), kwargs.get('dav')
    dust_curve = kwargs.get('dust_curve', 'cardelli')
    fsps_kwargs = kwargs.get('fsps_kwargs', {})

    zcbfile = config.path('bestzcb', field=brick, subfield=cell)

    if zcbfile:
        sfh = match.io.open_zcbfile(zcbfile)
        sfh['age_i'] = 10**sfh['log(age_i)']  # Linearize ages
        sfh['age_f'] = 10**sfh['log(age_f)']
        if agelim is not None:  # Don't test truthiness, agelim could be 0
            i = sfh['age_f'] <= agelim
            sfh = sfh[i]

        # Use the mean metallicity over the last 100 Myr, or the most
        # recent metallicity available before then.
        i = sfh['SFR'] > 0
        j = sfh['log(age_i)'] < 8
        if not np.any(i & j):
            # Use most recent metallicity available, or else assume solar
            logZ = sfh['[M/H]'][i][0] if np.any(i) else 0
        else:
            # 100 Myr mean value
            logZ = np.log10(np.mean(10**sfh['[M/H]'][i & j]))
        fsps_kwargs['zmet'] = astrogrid.flux.get_zmet(logZ)

        # Rescale 1st age bin
        sfh['SFR'][0] *= 1 - sfh['age_i'][0]/sfh['age_f'][0]
        sfh['age_i'][0] = 0

        age, sfr = (sfh['age_i'], sfh['age_f']), sfh['SFR']
        wave, spec, lum_ir = astrogrid.flux.calc_sed(
            sfr, age, av=av, dav=dav, dust_curve=dust_curve,
            fsps_kwargs=fsps_kwargs)
        mag = astrogrid.flux.calc_mag(
            wave, spec, band, dmod=config.DIST.distmod)
        flux = astrogrid.flux.mag2flux(mag, band)

    else:
        flux = np.nan

    return flux


def calc_mean_sfr(brick, cell):
    """Calculate the 100 Myr mean SFR for the SFH of the given cell.

    Parameters
    ----------
    brick, cell : int
        Brick and cell numbers.

    Returns
    -------
    float
        100 Myr mean SFR.

    """
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
    to a pixel and assemble an image.

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
    grid = astrogrid.Grid(
        (config.NROW, config.NCOL), func,
        [(brick, cell) for cell in cell_list], update=True)

    # Create a header
    coordfile = config.path('corners', field=brick)
    lon, lat = config.cornergrid(config.open_cornersfile(coordfile))
    x, y = grid.edges
    hdr = astrogrid.wcs.make_header(x, y, lon, lat)

    hdu = astropy.io.fits.PrimaryHDU(grid.data_grid, header=hdr)

    return hdu


def galex_pre(data, hdr, band):
    """Mask border pixels in a GALEX image and convert into flux units.

    The GALEX chip is not perfectly circular, so a small number of pixels
    near the border that are legitimately zero may be masked (set to NaN).

    Parameters
    ----------
    data : ndarray
        GALEX image data.
    hdr : astropy.io.fits.Header
        FITS header for the image.
    band : {'galex_fuv', 'galex_nuv'}
        The GALEX filter of the data. This is used for converting count
        rates into fluxes.

    Returns
    -------
    ndarray
    astropy.io.fits.Header

    """
    # Set border pixels to np.nan
    x = np.arange(data.shape[1]).reshape(1, -1) + 1  # Pixel coords of centers
    y = np.arange(data.shape[0]).reshape(-1, 1) + 1
    r = np.sqrt((x - config.GALEX_CHIP_X0)**2 + (y - config.GALEX_CHIP_Y0)**2)
    i = (r > config.GALEX_CHIP_RAD) & (data == 0)
    data = np.where(i, np.nan, data)  # A new array (instead of data[i])

    # Convert to flux
    data = astrogrid.flux.galex_cps2flux(data, band)

    return data, hdr


def galex_post(data, hdr, band, image_file, rectangle):
    """Measure and subtract the background flux level for a GALEX mosaic.

    The background is measured in the given image, assumed to be in units
    of flux per arcsec2, as the median pixel value within a rectangular
    aperture. This is converted into a total background flux per pixel in
    `data`, using `hdr` to determine the pixel area.

    Parameters
    ----------
    data : ndarray
        Mosaic image data.
    hdr : astropy.io.fits.Header
        FITS header for the image.
    band : {'galex_fuv', 'galex_nuv'}
        The GALEX filter of the data.
    image_file : str
        Path to the image in which to measure the background. The units are
        assumed to be in flux per arcsec2.
    rectangle : sequence
        The limits of the rectangular aperture as xmin, xmax, ymin, ymax.

    Returns
    -------
    ndarray
    astropy.io.fits.Header

    """
    # Pixel area
    dx, dy = astrogrid.wcs.calc_pixscale(hdr).arcsec
    pixarea = dx * dy  # arcsec2

    # Estimate the foreground level and subtract
    data2 = astropy.io.fits.getdata(image_file)  # flux per arcsec2
    x1, x2, y1, y2 = rectangle
    sample = data2[y1:y2,x1:x2].ravel()
    median = np.median(sample)
    val = median * pixarea  # flux
    data = data - val

    print('Median background level in {0:s}: {1:.3e} flux/arcsec2\n'
          'Total flux subtracted per pixel: {2:.3e}'
          .format(os.path.basename(image_file), median, val))

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
