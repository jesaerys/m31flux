"""

=================
`m31flux.fluxmap`
=================

Create modeled flux maps from SFH data.


Constants
---------

============= ======================================================
`FSPS_KWARGS` The IMF to use for all FSPS stellar population models.
============= ======================================================


Functions
---------

=================== ======================================================
`make_brick_images` Make one type of image for all bricks.
`make_mod_fuv_red`  Reddened FUV flux.
`make_mod_fuv_int`  Intrinsic FUV flux.
`galex_pre_fuv`     Mask border pixels and convert into flux units.
`galex_post_fuv`    Measure the background flux level and subtract it from
                    the image.
`make_galex_fuv`    GALEX FUV flux.
`make_mod_nuv_red`  Reddened NUV flux.
`make_mod_nuv_int`  Intrinsic NUV flux.
`galex_pre_nuv`     Mask border pixels and convert into flux units.
`galex_post_nuv`    Measure the background flux level and subtract it from
                    the image.
`make_galex_nuv`    GALEX NUV flux.
=================== ======================================================

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import astrogrid
import os

from . import config, util



FSPS_KWARGS = {'imf_type': astrogrid.flux.IMF_TYPE[config.IMF]}
"""The IMF to use for all FSPS stellar population models."""


def make_brick_images(kind, func):
    """Make one type of image for all bricks.

    Parameters
    ----------
    kind : str
        File kind (see `m31flux.config`).
    func : function
        Function that calculates the value of a cell. See
        `m31flux.util.make_brick_image`.

    """
    prefix = 'Calculating images (all bricks): '
    c = util.make_counter(prefix=prefix, end='  done')
    for brick in c(config.BRICK_LIST):
        hdu = util.make_brick_image(brick, func)
        filename = config.path(kind, field=brick)
        dirname = os.path.dirname(filename)
        try:
            os.makedirs(dirname)
        except OSError:
            pass
        try:
            hdu.writeto(filename)
        except IOError:
            os.remove(filename)
            hdu.writeto(filename)
    return



# FUV
# ---

def make_mod_fuv_red():
    """Reddened FUV flux."""
    def func(brick, cell):
        av, dav = config.EXTPAR_DICT[(brick, cell)]
        return util.calc_flux(
            brick, cell, 'galex_fuv', dmod=config.DIST.distmod, av=av, dav=dav,
            dust_curve=config.DUST_CURVE, fsps_kwargs=FSPS_KWARGS)

    print(
        '\n'
        'm31flux.fluxmap.make_mod_fuv_red\n'
        '--------------------------------'
        )
    make_brick_images('mod_fuv_red', func)

    print('Mosaicking images', end='')
    util.sys.stdout.flush()
    input_files = config.path('mod_fuv_red')
    mosaic_file = config.path('mod_fuv_red.mosaic')
    work_dir = config.path('mod_fuv_red.montage')
    header = None  # Have Montage create a header; this will be the master
    weights_file = config.path('weights')  # Save weights for this header
    astrogrid.mwe.mosaic(input_files, mosaic_file, work_dir, header=header,
                         weights_file=weights_file)
    print('  done')

    return


def make_mod_fuv_int():
    """Intrinsic FUV flux."""
    def func(brick, cell):
        return util.calc_flux(
            brick, cell, 'galex_fuv', agelimdmod=config.DIST.distmod,
            dust_curve=config.DUST_CURVE, fsps_kwargs=FSPS_KWARGS)

    print(
        '\n'
        'm31flux.fluxmap.make_mod_fuv_int\n'
        '--------------------------------'
        )
    make_brick_images('mod_fuv_int', func)

    print('Mosaicking images', end='')
    util.sys.stdout.flush()
    input_files = config.path('mod_fuv_int')
    mosaic_file = config.path('mod_fuv_int.mosaic')
    work_dir = config.path('mod_fuv_int.montage')
    header = config.path('mod_fuv_int.hdr')
    astrogrid.mwe.mosaic(input_files, mosaic_file, work_dir, header=header)
    print('  done')

    return


def galex_pre_fuv(data, hdr):
    """Mask border pixels and convert into flux units."""
    return util.galex_pre(data, hdr, 'galex_fuv')


def galex_post_fuv(data, hdr):
    """Measure the background flux level and subtract it from the image."""
    return util.galex_post(
        data, hdr, 'galex_fuv', config.path('galex_fuv.bg'),
        config.GALEX_BG_RECTANGLE)


def make_galex_fuv():
    """GALEX FUV flux."""
    print(
        '\n'
        'm31flux.fluxmap.make_galex_fuv\n'
        '------------------------------'
        )
    print('Mosaicking images')
    util.sys.stdout.flush()
    input_files = config.path('galex_fuv')
    mosaic_file = config.path('galex_fuv.mosaic')
    work_dir = config.path('galex_fuv.montage')
    header = config.path('galex_fuv.hdr')
    astrogrid.mwe.mosaic(input_files, mosaic_file, work_dir,
                         background_match=True, header=header,
                         postprocess=galex_post_fuv, preprocess=galex_pre_fuv)
    print('done')
    return



# NUV
# ---

def make_mod_nuv_red():
    """Reddened NUV flux."""
    def func(brick, cell):
        av, dav = config.EXTPAR_DICT[(brick, cell)]
        return util.calc_flux(
            brick, cell, 'galex_nuv', dmod=config.DIST.distmod, av=av, dav=dav,
            dust_curve=config.DUST_CURVE, fsps_kwargs=FSPS_KWARGS)

    print(
        '\n'
        'm31flux.fluxmap.make_mod_nuv_red\n'
        '--------------------------------'
        )
    make_brick_images('mod_nuv_red', func)

    print('Mosaicking images', end='')
    util.sys.stdout.flush()
    input_files = config.path('mod_nuv_red')
    mosaic_file = config.path('mod_nuv_red.mosaic')
    work_dir = config.path('mod_nuv_red.montage')
    header = config.path('mod_nuv_red.hdr')
    astrogrid.mwe.mosaic(input_files, mosaic_file, work_dir, header=header)
    print('  done')

    return


def make_mod_nuv_int():
    """Intrinsic NUV flux."""
    def func(brick, cell):
        return util.calc_flux(
            brick, cell, 'galex_nuv', agelimdmod=config.DIST.distmod,
            dust_curve=config.DUST_CURVE, fsps_kwargs=FSPS_KWARGS)

    print(
        '\n'
        'm31flux.fluxmap.make_mod_nuv_int\n'
        '--------------------------------'
        )
    make_brick_images('mod_nuv_int', func)

    print('Mosaicking images', end='')
    util.sys.stdout.flush()
    input_files = config.path('mod_nuv_int')
    mosaic_file = config.path('mod_nuv_int.mosaic')
    work_dir = config.path('mod_nuv_int.montage')
    header = config.path('mod_nuv_int.hdr')
    astrogrid.mwe.mosaic(input_files, mosaic_file, work_dir, header=header)
    print('  done')

    return


def galex_pre_nuv(data, hdr):
    """Mask border pixels and convert into flux units."""
    return util.galex_pre(data, hdr, 'galex_nuv')


def galex_post_nuv(data, hdr):
    """Measure the background flux level and subtract it from the image."""
    return util.galex_post(
        data, hdr, 'galex_nuv', config.path('galex_nuv.bg'),
        config.GALEX_BG_RECTANGLE)


def make_galex_nuv():
    """GALEX NUV flux."""
    print(
        '\n'
        'm31flux.fluxmap.make_galex_nuv\n'
        '------------------------------'
        )
    print('Mosaicking images')
    util.sys.stdout.flush()
    input_files = config.path('galex_nuv')
    mosaic_file = config.path('galex_nuv.mosaic')
    work_dir = config.path('galex_nuv.montage')
    header = config.path('galex_nuv.hdr')
    astrogrid.mwe.mosaic(input_files, mosaic_file, work_dir,
                         background_match=True, header=header,
                         postprocess=galex_post_nuv, preprocess=galex_pre_nuv)
    print('done')
    return
