"""

==============
`m31flux.util`
==============

General utilities.


Functions
---------

================== =========================================================
`safe_mkdir`       Create a directory only if it does not already exist.
`safe_symlink`     Create a symlink only if it does not already exist.
`calc_flux`        Calculate flux in a given filter for the SFH of the given
                   cell.
`make_brick_image` Apply a function to all cells in given brick, then
                   assign each cell to a pixel to assemble an image.
================== =========================================================

"""
import astrogrid
import astropy.io.fits
import astropy.stats.funcs
import errno
import match_wrapper as match
import montage_wrapper as montage
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


def safe_symlink(src, dst):
    """Create a symlink only if it does not already exist."""
    try:
        os.symlink(src, dst)
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
        wave, spec = astrogrid.flux.calc_sed(sfr, age, av=av, dav=dav, fsps_kwargs=fsps_kwargs)
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
    grid = astrogrid.Grid((config.NROW, config.NCOL), func,
                          [(brick, cell) for cell in config.CELL_LIST],
                          update=True)

    # Create a header
    coordfile = config.path('corners', field=brick)
    table = config.open_cornersfile(coordfile)
    lon, lat = config.list2grid(table)
    x, y = grid.edges
    hdr = astrogrid.wcs.make_header(x, y, lon, lat)

    hdu = astropy.io.fits.PrimaryHDU(data, header=hdr)

    return hdu


def galex_pre(data, hdr, band):
    # Mask border pixels (chip is not perfectly circular, so some pixels
    # near the border that are legitimately zero may be masked)
    y, x = np.indices(data.shape)  # array coordinates
    x, y = x + 1, y + 1  # pixel coordinates of the pixel centers
    r = np.sqrt((x - config.GALEX_CHIP_X0)**2 + (y - config.GALEX_CHIP_Y0)**2)
    i = (r > config.GALEX_CHIP_RAD) & (data == 0)
    data[i] = np.nan

    # Convert to flux
    data = astrogrid.flux._galex_cps2flux(data, band)

    return data, hdr


def galex_post(data, hdr, band):
    sig = 1  # Sigma clipping limit

    # This somewhat small area is apparently devoid of stars in the NUV image
    x1, x2 = 1550, 1715
    y1, y2 = 2290, 2475

    # Estimate the foreground level and subtract
    data_box = data[y1:y2,x1:x2].ravel()
    mask = sigma_clip(data_box, sig=sig, iters=None, cenfunc=np.mean).mask
    val = np.mean(data_box[mask])
    data = data - val

    return data, hdr


# REPLACE MOST OF THIS WITH montage_wrapper.mosaic!
def make_mosaic(kind, make_header=False, cdelt=None, background_match=False,
                level_only=True, preprocess=None, postprocess=None):
    """Make a mosiac.

    Parameters
    ----------
    kind : str
        File kind (see `m31flux.config`).
    make_header : bool, optional
        If True, a template header for the mosaic will be created based on
        the input images. If False (default), then the file pointed to by
        `config.path` is used instead.
    cdelt : float, optional
        Set the pixel scale (deg/pix) of the mosaic. If None, the mosaic of
        Alexia's pixel regions will have a pixel scale of 23.75 arcsec (in
        both x and y).
    background_match, level_only : bool, optional
        Background matching parameters. See `montage_wrapper.mosaic` for
        details.
    preprocess, postprocess : function, optional
        Functions for processing the raw input images before the input
        density images are created and after the final mosaic is created.
        The function arguments should be the image data array and the image
        header (`astropy.io.fits.Header`), and the return values should be
        the same. Default is None.

    Notes
    -----
    - Input and reprojected maps were compared for 21 PHAT bricks. For a
      given brick, the percent difference between the pixel sum of the
      native map and the pixel sum of the reprojected map was less than
      ~0.01%. The average percent difference was ~0.001%.

    - mProjExec and mProject produce nearly identical results (within
      ~0.01%).

    - Settings that affect reprojection: mProject has the 'z' flag, which
      controls the drizzle factor. 1 seems to be the default, and appears
      to give the best results anyway so there is no reason to change it.
      mProjExec doesn't really have any settings that change the
      reprojection results.

    - Pixel area varies slightly accross brick 15 in both the native and
      the reprojected images. The effect causes only ~0.001% difference
      between the minimum and maximum areas, however, so the pixels can be
      safely treated as having constant area.

    - Three pixel area calculation methods were tested: 1) calculate the
      areas of all pixels from the coordinates of their corners assuming
      spherical rectangles, 2) calculate the areas of all pixels from their
      x and y scales assuming planar rectangles, and 3) same as 2, but only
      calculate the area for the reference pixel and assume that value for
      all of the pixels. All three area calculation methods agree to within
      ~0.01% to ~0.001%. Might as well just use the simplest method (3).

    """
    input_files = config.path(kind)
    density_files = config.path('{:s}.density'.format(kind))
    for input_file, density_file in zip(input_files, density_files):
        # Load input image, perform initial processing
        data, hdr = astropy.io.fits.getdata(input_file, header=True)
        if preprocess is not None:
            data, hdr = preprocess(data, hdr)

        # Divide by pixel area to create a density image
        dx, dy = astrogrid.wcs.calc_pixscale(hdr, ref='crpix')
        area = dx * dy  # arcsec2
        data = data / area

        # Write
        hdu = astropy.io.fits.PrimaryHDU(data, header=hdr)
        dirname = os.path.dirname(density_file)
        safe_mkdir(dirname)
        if os.path.exists(density_file):
            os.remove(density_file)
        hdu.writeto(density_file)

    # Template header
    header_file = config.path('{:s}.hdr'.format(kind))
    if make_header:
        dirname = os.path.dirname(header_file)
        safe_mkdir(dirname)
        montage.mMakeHdr(meta_file1, header_file, cdelt=cdelt)

    # Make mosaic
    input_dir = os.path.dirname(density_files[0])
    montage.mosaic(input_dir, output_dir, header=header_file,
                   background_match=background_match, level_only=level_only,
                   exact_size=True)

    ## Metadata for input images
    #meta_file1 = os.path.join(input_dir, 'input.tbl')
    #montage.mImgtbl(input_dir, meta_file1, corners=True)


    # Reproject to the template
    reproject_files = config.path('{:s}.reproject'.format(kind))
    reproject_dir = os.path.dirname(reproject_files[0])
    safe_mkdir(reproject_dir)
    stat_file = os.path.join(reproject_dir, 'stats.tbl')
    montage.mProjExec(meta_file1, header_file, reproject_dir,
                      stat_file, raw_dir=input_dir)

    # Metadata for reprojected images
    meta_file2 = os.path.join(reproject_dir, 'reproject.tbl')
    montage.mImgtbl(reproject_dir, meta_file2, corners=True)

    # Build density mosaic
    mosaic_file1 = config.path('{:s}.reproject.add'.format(kind))
    dirname = os.path.dirname(mosaic_file1)
    safe_mkdir(dirname)
    montage.mAdd(meta_file2, header_file, mosaic_file1,
                 img_dir=reproject_dir, exact=True)

    # Final mosaic
    data, hdr = astropy.io.fits.getdata(mosaic_file1, header=True)
    area_file = config.path('{:s}.area.add'.format(kind))
    area = astropy.io.fits.getdata(area_file) * (180/np.pi*3600)**2 # arcsec2
    data = data * area
    mosaic_file2 = config.path('{:s}.add'.format(kind))
    dirname = os.path.dirname(mosaic_file2)
    safe_mkdir(dirname)
    hdu = astropy.io.fits.PrimaryHDU(data, header=hdr)
    if os.path.exists(mosaic_file2):
        os.remove(mosaic_file2)
    hdu.writeto(mosaic_file2)

    return None


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
