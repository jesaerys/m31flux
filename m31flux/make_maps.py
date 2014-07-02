"""

===================
`m31flux.make_maps`
===================

Create modeled flux maps from SFH data.


Constants
---------

============= ======================================================
`FSPS_KWARGS` The IMF to use for all FSPS stellar population models.
============= ======================================================


Functions
---------

======================= ===================================================
`calc_flux_mod_fuv_red` Calculate reddened FUV flux for the given brick and
                        cell.
`calc_flux_mod_fuv_int` Calculate intrinsic FUV flux for the given brick
                        and cell.
======================= ===================================================

"""
import astrogrid
import os

from . import config, util



FSPS_KWARGS = {'imf_type': astrogrid.flux.IMF_TYPE[config.IMF]}
"""The IMF to use for all FSPS stellar population models."""



# FUV
# ---

def calc_flux_mod_fuv_red(brick, cell):
    """Calculate reddened FUV flux for the given brick and cell."""
    av, dav = config.EXTPAR_DICT[(brick, cell)]
    return util.calc_flux(brick, cell, 'galex_fuv', dmod=config.DMOD,
                          av=av, dav=dav, fsps_kwargs=FSPS_KWARGS)


def calc_flux_mod_fuv_int(brick, cell):
    """Calculate intrinsic FUV flux for the given brick and cell."""
    return util.calc_flux(brick, cell, 'galex_fuv', agelimdmod=config.DMOD,
                          fsps_kwargs=FSPS_KWARGS)


def prep_galex_fuv(data, hdr):
    return util.prep_galex(data, hdr, 'galex_fuv')



# NUV
# ---
def calc_flux_mod_nuv_red(brick, cell):
    av, dav = config.EXTPAR_DICT[(brick, cell)]
    return util.calc_flux(brick, cell, 'galex_nuv', dmod=config.DMOD,
                          av=av, dav=dav, fsps_kwargs=FSPS_KWARGS)


def calc_flux_mod_nuv_int(brick, cell):
    return util.calc_flux(brick, cell, 'galex_nuv', dmod=config.DMOD,
                          fsps_kwargs=FSPS_KWARGS)


def prep_galex_nuv(data, hdr):
    return util.prep_galex(data, hdr, 'galex_nuv')



# Main
# ----
def make_brick_images(kind, func):
    """Make one type of image for all bricks.

    Parameters
    ----------
    kind : str
        File kind (see `m31flux.config`).
    func : function
        Function that calculates the value of a cell (see
        `m31flux.make_brick_image`).

    """
    for brick in config.BRICK_LIST:
        hdu = util.make_brick_image(brick, func)
        filename = config.path(kind, field=brick)
        dirname = os.path.dirname(filename)
        util.safe_mkdir(dirname)
        if os.path.exists(filename):
            os.remove(filename)
        hdu.writeto(filename)
    return None


def main():
    make_brick_images('mod_fuv_red', calc_flux_mod_fuv_red)
    #util.make_mosaic('mod_fuv_red', make_header=True)
    #make_brick_images('mod_fuv_int', calc_flux_mod_fuv_int)
    #util.make_mosaic('mod_fuv_int')
    #util.make_mosaic('galex_fuv', procinput=prep_galex_fuv)

    #make_brick_images('mod_nuv_red', calc_flux_mod_nuv_red)
    #util.make_mosaic('mod_nuv_red')
    #make_brick_images('mod_nuv_int', calc_flux_mod_nuv_int)
    #util.make_mosaic('mod_nuv_int')
    #util.make_mosaic('galex_nuv', procinput=prep_galex_nuv)

    return None


if __name__ == '__main__':
    main()




def make_mean_sfr_100myr_brickmaps():
    make_brickmaps_base('mean_sfr_100myr', calc_mean_sfr)
    return None


def make_mean_sfr_100myr_mosaic(cdelt=None):
    work_dir = os.path.join(config.ANALYSIS_DIR, '_mean_sfr_100myr')
    make_mosaic_base('mean_sfr_100myr', work_dir, cdelt=cdelt)
    return None


def __old__make_mod_fuv_brickmaps(reddened=False):
    # wrap function based on `reddened`
    def func(brick, pixel):
        return calc_flux(brick, pixel, 'galex_FUV', reddened=reddened)

    if reddened:
        kind = 'mod_fuv_red'
    else:
        kind = 'mod_fuv'

    make_brickmaps_base(kind, func)
    return None


def __old__make_mod_fuv_mosaic(attenuated=False, cdelt=None):
    """
    reddened : bool
        If True, create a mosaic for reddened FUV flux. Else make a
        mosaic for intrinsic FUV flux (default).

    """
    if attenuated:
        kind = 'mod_fuv_red'
    else:
        kind = 'mod_fuv'

    work_dir = os.path.join(config.ANALYSIS_DIR, '_{:s}'.format(kind))
    work_dir = config.path('{0:s}_montage_dir'.append(kind))
    make_mosaic_base(kind, work_dir, cdelt=cdelt)
    return None


def make_mod_fuv_brickmaps(attenuated=False):
    sps = fsps.StellarPopulation()
    sps.params['sfh'] = 0
    sps.params['sfh'] = 0
    sps.params['zmet'] = 4

    # wrap function based on `sps` and `attenuated`
    def func(brick, pixel):
        return calc_flux(brick, pixel, sps, 'galex_FUV', attenuated=attenuated)

    if attenuated:
        kind = 'mod_fuv_attenuated'
    else:
        kind = 'mod_fuv_intrinsic'

    make_brickmaps_base(kind, func)
    return None


def make_galex_uv_mosaic(band):
    """`band` is either 'fuv' or 'nuv'."""
    if band == 'fuv':
        a = 1.40e-15
    elif band == 'nuv':
        a = 2.06e-16

    # Set up a working directory and subdirectories
    work_dir = os.path.join(config.MAP_DIR, '_galex_{:s}'.format(band))
    input_dir = os.path.join(work_dir, 'input')
    proj_dir = os.path.join(work_dir, 'reproject')
    for path in (work_dir, input_dir, proj_dir):
        util.safe_mkdir(path)

    # Prepare input images
    for field in config.GALEX_FIELD_LIST:
        filename = 'PS_M31_MOS{0:s}-{1:s}d-int.fits'.format(field, band[0])
        mapfile = os.path.join(config.GALEX_DIR, filename)
        hdr = fits.getheader(mapfile)
        data = fits.getdata(mapfile)

        # Set border pixels to NaN
        x0, y0 = 1920, 1920
        r = 1400
        y, x = np.indices(data.shape)
        incircle = (x+1)**2 + (y+1)**2 <= r**2
        iszero = data == 0
        data[-incircle & iszero] = np.nan

        # Convert input image units from cps to flux
        data = data * a

        # Write image
        mapfile = os.path.join(input_dir, filename)
        hdu = fits.PrimaryHDU(data, header=hdr)
        hdu.writeto(mapfile)

    # Metadata for input images
    metafile1 = os.path.join(work_dir, 'input.tbl')
    montage.mImgtbl(input_dir, metafile1, corners=True)

    # Copy template header
    hdrfile = '/Users/Jake/Research/PHAT/sfhmaps/data/map/_mod_fuv_attenuated/template.hdr'

    # Reproject to the template
    statfile = os.path.join(work_dir, 'stats.tbl')
    montage.mProjExec(metafile1, hdrfile, proj_dir, statfile, raw_dir=input_dir)

    # Metadata for reprojected images
    metafile2 = os.path.join(work_dir, 'reproject.tbl')
    montage.mImgtbl(proj_dir, metafile2, corners=True)

    # Build final mosaic
    mosaicfile = '/Users/Jake/Research/PHAT/sfhmaps/data/map/galex_{:s}.fits'.format(band)
    montage.mAdd(metafile2, hdrfile, mosaicfile, img_dir=proj_dir, exact=True)

    return None


def make_galex_uv_brickmaps(band, clean=False):
    """`band` is either 'fuv' or 'nuv'."""
    for brick in config.BRICK_LIST:

        # Write template header
        brickmap = config.get_file(brick=brick, kind='mod_fuv_attenuated')
        work_dir1 = os.path.dirname(brickmap)
        hdrfile = os.path.join(work_dir1, 'template.hdr')
        if clean:
            subprocess.call('rm {:s}'.format(hdrfile), shell=True)
        else:
            montage.mGetHdr(brickmap, hdrfile)

        # Metadata for input images
        work_dir2 = os.path.join(config.MAP_DIR, '_galex_{:s}'.format(band))
        input_dir = os.path.join(work_dir2, 'input')
        metafile1 = os.path.join(work_dir1, 'input.tbl')
        if clean:
            subprocess.call('rm {:s}'.format(metafile1), shell=True)
        else:
            montage.mImgtbl(input_dir, metafile1, corners=True)

        # Reproject to the template
        proj_dir = os.path.join(work_dir1, 'reproject')
        statfile = os.path.join(work_dir1, 'stats.tbl')
        if clean:
            subprocess.call('rm -rf {:s}'.format(proj_dir), shell=True)
            subprocess.call('rm {:s}'.format(statfile), shell=True)
        else:
            util.safe_mkdir(proj_dir)
            montage.mProjExec(metafile1, hdrfile, proj_dir, statfile, raw_dir=input_dir)

        # Metadata for reprojected images
        metafile2 = os.path.join(work_dir1, 'reproject.tbl')
        if clean:
            subprocess.call('rm {:s}'.format(metafile2), shell=True)
        else:
            montage.mImgtbl(proj_dir, metafile2, corners=True)

        # Build final mosaic
        mosaicfile = os.path.join(work_dir1, 'b{0:02d}_galex_{1:s}.fits'.format(brick, band))
        if clean:
            mosaicfile = os.path.splitext(mosaicfile)[0]
            subprocess.call('rm {:s}_area.fits'.format(mosaicfile), shell=True)
        else:
            montage.mAdd(metafile2, hdrfile, mosaicfile, img_dir=proj_dir, exact=True)

    return None


def main():
    #make_mean_sfr_100myr_brickmaps()
    #make_mean_sfr_100myr_mosaic()

    #make_mod_fuv_brickmaps(attenuated=True)
    #make_mod_fuv_mosaic(attenuated=True)

    #make_galex_uv_mosaic('nuv')
    #make_galex_uv_brickmaps('fuv', clean=True)

    # Tests
    # -----
    #for brick in config.BRICK_LIST:
    #    run_imwcs(brick)
    #compare_wcs()

    return None



# Tests
# -----

def flux_conservation_test():
    """The following test shows that SFR is *not* conserved in reprojection."""
    filename0 = '/Users/Jake/Research/PHAT/sfh_maps/data/map/b15/b15_mean_sfr_100myr.fits'
    filename = '/Users/Jake/Research/PHAT/sfh_maps/data/map/_mean_sfr_100myr/reproject/hdu0_b15_mean_sfr_100myr.fits'

    data0 = fits.getdata(filename0)  # input image
    data = fits.getdata(filename)  # reprojected image

    np.nansum(data0)  # 0.03619
    np.nansum(data)  # 0.043404; 20% larger!!!


    """According to this blog post
    (http://montageblog.wordpress.com/2011/06/24/does-montage-conserve-flux-when-changing-the-image-resolution/),
    it may be necessary to distinguish between total/absolute quantities
    and densities (the post refers to "flux" vs. "flux density", though it
    is not clear what these terms actually mean... perhaps luminosity vs.
    flux, where flux is luminosity per unit area?). If I understand this
    correctly, then I should be able to normalize the input SFR maps by
    pixel area, and then recover SFR in the reprojected maps by multiplying
    by the area of a reprojected pixel. I test this below."""

    # set up directories for reprojection
    work_dir = '/Users/Jake/Research/PHAT/sfh_maps/data/map/test'
    input_dir = os.path.join(work_dir, 'input')
    proj_dir = os.path.join(work_dir, 'reproject')
    for path in (work_dir, input_dir, proj_dir):
        util.safe_mkdir(path)

    # make SFR density map
    testname0 = os.path.join(work_dir, 'input/test.fits')
    hdr0 = fits.getheader(filename0)
    from astropy import wcs
    w0 = wcs.WCS(hdr0)
    xy0 = np.mgrid[0:hdr0['naxis2']+1,0:hdr0['naxis1']+1][::-1] + 0.5
    ad0 = w0.wcs_pix2world(xy0.T.reshape(-1, 2), 1).reshape(xy0.T.shape).T
    area0 = []
    for i in range(data0.shape[0]):
        for j in range(data0.shape[1]):
            lon = np.array([ad0[0,i,j], ad0[0,i+1,j], ad0[0,i+1,j+1], ad0[0,i,j+1]])
            lat = np.array([ad0[1,i,j], ad0[1,i+1,j], ad0[1,i+1,j+1], ad0[1,i,j+1]])
            area0.append(util.sparea(lon, lat))
    area0 = np.array(area0).reshape(data0.shape)
    datat0 = data0 / area0  # SFR deg-2
    hdu = fits.PrimaryHDU(datat0, header=hdr0)
    hdu.writeto(testname0)

    # input image metadata
    metafile1 = os.path.join(work_dir, 'input.tbl')
    montage.mImgtbl(input_dir, metafile1, corners=True)

    # use the same header
    hdrfile = os.path.join(work_dir, 'template.hdr')
    import subprocess
    subprocess.call('cp /Users/Jake/Research/PHAT/sfh_maps/data/map/_mean_sfr_100myr/template.hdr {:s}'.format(hdrfile), shell=True)

    # reproject
    statfile = os.path.join(work_dir, 'stats.tbl')
    montage.mProjExec(metafile1, hdrfile, proj_dir, statfile, raw_dir=input_dir)

    # convert SFR density back to SFR
    testname = os.path.join(work_dir, 'reproject/hdu0_test.fits')
    datat = fits.getdata(testname)
    hdr = fits.getheader(testname)
    w = wcs.WCS(hdr)
    xy = np.mgrid[0:hdr['naxis2']+1,0:hdr['naxis1']+1][::-1] + 0.5
    ad = w.wcs_pix2world(xy.T.reshape(-1, 2), 1).reshape(xy.T.shape).T
    area= []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            lon = np.array([ad[0,i,j], ad[0,i+1,j], ad[0,i+1,j+1], ad[0,i,j+1]])
            lat = np.array([ad[1,i,j], ad[1,i+1,j], ad[1,i+1,j+1], ad[1,i,j+1]])
            area.append(util.sparea(lon, lat))
    area = np.array(area).reshape(data.shape)
    datat *= area

    # conservation test
    np.nansum(data0)  # 0.03619
    np.nansum(datat)  # 0.03840; 6% larger

    """
    6% is much better, but it's still large enough that I'm concerned that
    I'm doing something wrong.

    ***UPDATE***
    I think the key is to multiply the reprojected image by the
    accompanying 'area' file. This accounts for the fact that edge pixels
    in the reprojected image don't have 100% coverage in the input image.

    """
    return None


def plot_img(mode):
    from astropy.io import fits
    from matplotlib import pyplot as plt

    if mode == 'mean_sfr_100myr':
        imgfile = 'mean_sfr_100myr.fits'
        outfile = 'mean_sfr_100myr.pdf'
        fmin, fmax = 0, 1
        a = 1e3
    elif mode == 'mod_fuv_attenuated':
        imgfile = 'mod_fuv_attenuated.fits'
        outfile = 'mod_fuv_attenuated.pdf'
        fmin, fmax = 0, 1
        a = 1e15
    elif mode == 'galex_fuv':
        imgfile = 'galex_fuv.fits'
        outfile = 'galex_fuv.pdf'
        fmin, fmax = 0, 1
        a = 5e17
    elif mode == 'galex_nuv':
        imgfile = 'galex_nuv.fits'
        outfile = 'galex_nuv.pdf'
        fmin, fmax = 0, 1
        a = 5e17
    data = fits.getdata(imgfile)

    data = data[::-1].T  # landscape instead of portrait
    data = np.log10(a*data + 1)

    fig_dx = 10.0
    aspect = float(data.shape[0])/data.shape[1]
    fig_dy = fig_dx * aspect
    fig = plt.figure(figsize=(fig_dx, fig_dy))
    ax = fig.add_axes([0, 0, 1, 1])

    cmap = plt.cm.gist_heat_r
    vmin = (np.nanmax(data)-np.nanmin(data))*fmin + np.nanmin(data)
    vmax = (np.nanmax(data)-np.nanmin(data))*fmax + np.nanmin(data)
    img = ax.imshow(data, interpolation='nearest', origin='lower',
                    cmap=cmap, vmin=vmin, vmax=vmax)

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.savefig(outfile)

    return None
#plot_img('mean_sfr_100myr')
#plot_img('galex_fuv')
#plot_img('galex_nuv')
#plot_img('mod_fuv_attenuated')


def plot_img2():
    from astropy.io import fits
    from matplotlib import pyplot as plt

    imgfile1 = 'mod_fuv_attenuated.fits'
    imgfile2 = 'galex_nuv.fits'
    outfile = 'nuv_ratio'
    data1 = fits.getdata(imgfile1)
    data2 = fits.getdata(imgfile2)

    data = data1 / data2
    hdr = fits.getheader(imgfile1)
    hdu = fits.PrimaryHDU(data, header=hdr)
    hdu.writeto('{:s}.fits'.format(outfile))

    data = data[::-1].T  # landscape instead of portrait

    fig_dx = 10.0
    aspect = float(data.shape[0])/data.shape[1]
    fig_dy = fig_dx * aspect
    fig = plt.figure(figsize=(fig_dx, fig_dy))
    ax = fig.add_axes([0, 0, 1, 1])

    fmin, fmax = 0, 1
    cmap = plt.cm.gist_heat_r
    vmin = (np.nanmax(data)-np.nanmin(data))*fmin + np.nanmin(data)
    vmax = (np.nanmax(data)-np.nanmin(data))*fmax + np.nanmin(data)
    img = ax.imshow(data, interpolation='nearest', origin='lower',
                    cmap=cmap, vmin=vmin, vmax=vmax)

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.savefig('{:s}.pdf'.format(outfile))

    return None
#plot_img2()


def plot_grid(**kwargs):
    plt.clf()
    gs = gridspec.GridSpec(3, 1, **kwargs)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax3 = plt.subplot(gs[2])
