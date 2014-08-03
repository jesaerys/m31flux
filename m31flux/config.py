"""

================
`m31flux.config`
================

Configuration for `m31flux`.

This module contains constants and functions that locate data (both
existing files and files created during analysis) and define how it is
structured. The main feature is a `path` function to keep track of all
project files. Every tracked file has an associated "kind" [1]_, a string
that uniquely identifies the file or the group/category of files it belongs
to. This system completely eliminates the need to hardcode filenames in any
module or analysis script; only this config file (namely, the `_KIND_DICT`
dictionary) has to be edited to use the `m31flux` package.

Also defined are various project-specific things like brick grid
parameters, galaxy and modeling parameters, and parsers for nonstandard
files.

.. [1] Use of the word "kind" will generally refer to the `path` argument.
   Also, a string in single quotes used as an adjective for a file
   implicitly refers to the file's kind. For example, the phrase "'sfh'
   files" means files of kind 'sfh', i.e., the files returned by
   ``path('sfh')``.


Constants
---------

==================== ========================================================
`BRICK_LIST`         Names of all PHAT bricks available in the project in
                     numerical order.
`NROW`               Number of rows in a brick grid.
`NCOL`               Number of columns in a brick grid.
`CELL_LIST`          Names of all cells in a brick grid in numerical order.
`SORT`               Indices that sort a list of cells in numerical order
                     into grid order, and vice versa.
`GALEX_TILE_LIST`    Names of all GALEX DIS tiles/fields covering the PHAT
                     survey area in numerical order.
`GALEX_CHIP_X0`      Approximate x pixel coordinate of chip center.
`GALEX_CHIP_Y0`      Approximate y pixel coordinate of chip center.
`GALEX_CHIP_RAD`     Approximate chip radius in pixels, slightly
                     undersized.
`GALEX_BG_RECTANGLE` xmin, xmax, ymin, and ymax (array) coordinates of the
                     rectangular aperture for measuring the FUV and NUV
                     background/foreground level.
`DIST`               Distance to M31.
`INC`                M31 disk inclination angle.
`IMF`                Initial mass function.
`DUST_CURVE`         Dust curve.
`_KIND_DICT`         Paths to various file kinds.
`EXTPAR_DICT`        Dictionary of the best-fit Av and dAv extinction
                     parameters for every cell.
`MISSING_CELLS`      List of (brick, cell) tuples for where there is no data.
==================== ========================================================


Functions
---------

================== =======================================================
`open_extparfile`  Create a table from an 'extpar' file.
`open_cornersfile` Create a table from a 'corners' file.
`cornergrid`       Reshape a table of cell corner coordinates into a grid.
================== =======================================================

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import astrogrid
import astropy.coordinates
import astropy.table
import numpy as np
import os



# PHAT bricks and the brick grid system
# =====================================

BRICK_LIST = np.array([2]+range(4,24))
"""Names of all PHAT bricks available in the project in numerical order."""

NROW = 15
"""Number of rows in a brick grid."""

NCOL = 30
"""Number of columns in a brick grid."""

CELL_LIST = np.arange(NROW*NCOL) + 1
"""Names of all cells in a brick grid in numerical order."""

SORT = np.arange(NROW*NCOL).reshape((NROW, NCOL))[::-1].ravel()
"""Indices that sort a list of cell names in numerical order into grid
order, and vice versa.

"Grid order" refers to the order of the cells in a flattened array (e.g.,
`numpy.ravel`).

"""



# GALEX
# =====

GALEX_TILE_LIST = [0, 7, 8, 9, 10]
"""Names of all GALEX DIS tiles/fields covering the PHAT survey area in
numerical order.

"""

GALEX_CHIP_X0 = 1920
"""Approximate x pixel coordinate of chip center."""

GALEX_CHIP_Y0 = 1920
"""Approximate y pixel coordinate of chip center."""

GALEX_CHIP_RAD = 1400
"""Approximate chip radius in pixels, slightly undersized."""

GALEX_BG_RECTANGLE = (75, 145, 65, 135)
"""xmin, xmax, ymin, and ymax (array) coordinates of the rectangular
aperture for measuring the FUV and NUV background/foreground levels.

The coordinates refer specifically to the 'galex_fuv.bg' and 'galex_nuv.bg'
files.

"""



# M31 and other parameters
# ========================

DIST = astropy.coordinates.Distance(distmod=24.47)
"""Distance to M31.

McConnachie, A. W., Irwin, M. J., Ferguson, A. M. N., et al. 2005, MNRAS,
356, 979.

"""


INC = astropy.coordinates.Angle(78.0, unit='deg')
"""M31 disk inclination angle.

Tully, R. B. 1994, yCat, 7145, 0.

"""

IMF = 'Kroupa'
"""Initial mass function."""

DUST_CURVE = 'cardelli'
"""Dust curve.

Cardelli, J.~A., Clayton, G.~C., \& Mathis, J.~S.\ 1989, \apj, 345, 245.

"""



# Paths
# =====

_PROJECT_DIR = '/Users/Jake/Research/PHAT/m31flux'
_SFH_DIR = os.path.join(_PROJECT_DIR, 'sfh')
_ANALYSIS_DIR = os.path.join(_PROJECT_DIR, 'analysis')
_GALEX_DIR = '/Users/Jake/Research/Storage/M31/GALEX/DIS'


_KIND_DICT = {
    'extpar':  (_SFH_DIR, 'B{0:02d}', 'b{0:02d}_region_AvdAv.dat'),
    'corners': (_SFH_DIR, 'B{0:02d}',
                'M31-B{0:02d}_15x30_subregion-exact-vertices.dat'),
    'phot':    (_SFH_DIR, 'B{0:02d}', 'phot',
                'M31-B{0:02d}_15x30-{1:03d}.gst.match'),
    'sfh':     (_SFH_DIR, 'B{0:02d}', 'sfh',
                'M31-B{0:02d}_15x30-{1:03d}_{2:.1f}-{3:.1f}_best.sfh'),
    'cmd':     (_SFH_DIR, 'B{0:02d}', 'cmd',
                'M31-B{0:02d}_15x30-{1:03d}_{2:.1f}-{3:.1f}_best.sfh.cmd'),
    'bestzcb': (_ANALYSIS_DIR, 'b{0:02d}',
                'bestzcb', 'b{0:02d}-{1:03d}_best.zcb'),

    'mod_fuv_int': (_ANALYSIS_DIR, 'b{0:02d}', 'b{0:02d}_mod_fuv_int.fits'),
    'mod_fuv_int.mosaic': (_ANALYSIS_DIR, 'mod_fuv_int.fits'),
    'mod_fuv_int.montage': (_ANALYSIS_DIR, '_mod_fuv_int'),
    'mod_fuv_int.hdr': 'mod_fuv_red.hdr',

    'mod_fuv_red': (_ANALYSIS_DIR, 'b{0:02d}', 'b{0:02d}_mod_fuv_red.fits'),
    'mod_fuv_red.mosaic': (_ANALYSIS_DIR, 'mod_fuv_red.fits'),
    'mod_fuv_red.montage': (_ANALYSIS_DIR, '_mod_fuv_red'),
    'mod_fuv_red.hdr': (_ANALYSIS_DIR, '_mod_fuv_red', 'template.hdr'),
    'weights': (_ANALYSIS_DIR, 'weights.fits'),

    'galex_fuv': (_GALEX_DIR, 'PS_M31_MOS{0:02d}-fd-int.fits'),
    'galex_fuv.mosaic': (_ANALYSIS_DIR, 'galex_fuv.fits'),
    'galex_fuv.montage': (_ANALYSIS_DIR, '_galex_fuv'),
    'galex_fuv.hdr': 'mod_fuv_red.hdr',
    'galex_fuv.bg': (_ANALYSIS_DIR, '_galex_fuv', 'corrected',
                     'hdu0_PS_M31_MOS07-fd-int_density.fits'),

    'mod_nuv_int': (_ANALYSIS_DIR, 'b{0:02d}', 'b{0:02d}_mod_nuv_int.fits'),
    'mod_nuv_int.mosaic': (_ANALYSIS_DIR, 'mod_nuv_int.fits'),
    'mod_nuv_int.montage': (_ANALYSIS_DIR, '_mod_nuv_int'),
    'mod_nuv_int.hdr': 'mod_fuv_red.hdr',

    'mod_nuv_red': (_ANALYSIS_DIR, 'b{0:02d}', 'b{0:02d}_mod_nuv_red.fits'),
    'mod_nuv_red.mosaic': (_ANALYSIS_DIR, 'mod_nuv_red.fits'),
    'mod_nuv_red.montage': (_ANALYSIS_DIR, '_mod_nuv_red'),
    'mod_nuv_red.hdr': 'mod_fuv_red.hdr',

    'galex_nuv': (_GALEX_DIR, 'PS_M31_MOS{0:02d}-nd-int.fits'),
    'galex_nuv.mosaic': (_ANALYSIS_DIR, 'galex_nuv.fits'),
    'galex_nuv.montage': (_ANALYSIS_DIR, '_galex_nuv'),
    'galex_nuv.hdr': 'mod_fuv_red.hdr',
    'galex_nuv.bg': (_ANALYSIS_DIR, '_galex_nuv', 'corrected',
                     'hdu0_PS_M31_MOS07-nd-int_density.fits'),
    }
"""Paths to various file kinds.

This dictionary is used by `_path` and defines the locations of all file
kinds in the project. Each kind is specified by a tuple of path elements
(strings) which are later combined into a single string using
`os.path.join`.

The elements may contain format strings as placeholders for inserting one
or more of the values in a standard tuple, ``(field, subfield, av, dav)``,
when `_path` is actually called. The first two values are documented in
`_path`, and the last two refer to the Av and dAv extinction parameters.
The values are referenced in a format string according to their tuple
index. For example, putting ``'{1:02d}'`` in a path element will insert the
value of ``subfield`` formatted as a padded two-column integer.

To link one file kind to another, set the kind's value to the name of the
other kind. This is useful if two kinds actually share the same file.

"""


def _path(kind, extpar_dict, **kwargs):
    """File paths for the `m31flux` project.

    Parameters
    ----------
    kind : str
        File kind. Every tracked file has an associated kind, a string that
        uniquely identifies the file or the group/category of files it
        belongs to. The exact set of files returned is controlled using
        keyword arguments. Valid kinds, listed below, are defined in
        scripts/main.py in the main m31flux repository (optional suffixes
        are given in brackets).

        - 'corners'
        - 'extpar'
        - 'phot'
        - 'sfh'
        - 'cmd'
        - 'bestzcb'
        - 'mod_fuv_int[.mosaic|.montage|.hdr]'
        - 'mod_fuv_red[.mosaic|.montage|.hdr]'
        - 'galex_fuv[.mosaic|.montage|.hdr|.bg]'
        - 'mod_nuv_int[.mosaic|.montage|.hdr]'
        - 'mod_nuv_red[.mosaic|.montage|.hdr]'
        - 'galex_nuv[.mosaic|.montage|.hdr|.bg]'
        - 'weights'

    extpar_dict : dict
        Dictionary of brick grid extinction parameter values as (Av, dAv)
        tuples keyed by (brick, cell) tuples. This dictionary is required
        because some paths contain Av and dAv parameter values.
    field, subfield : int or sequence of ints, optional
        Return the file(s) of type `kind` for a combination of fields and
        subfields. None (default) is equivalent to a list of all fields or
        all subfields. "field" and "subfield" are generic terms that
        represent more specific organizational terms for different kinds of
        files:

        ================ ===== =============================
        type of data     field subfield
        ================ ===== =============================
        PHAT             brick field
        PHAT brick grids brick cell (a.k.a. pixel or region)
        GALEX            tile
        ================ ===== =============================

        These keywords are ignored for kinds that do not have fields or
        subfields.

    missing : sequence, optional
        (brick, cell) tuples for which no SFH data are available. This is
        used to avoid returning paths that do not exist. Default is an
        empty list.
    fillsubfield : bool, optional
        If True (default), any subfields listed in `missing` are added to
        the returned path list as None. If False, all missing subfields are
        omitted from the path list.

    Returns
    -------
    str or list of strs
        Paths matching the given `kind`, `field`, and `subfield`. A list is
        returned when there are multiple matches. A single string is
        returned when there is only one match, unless `field` and/or
        `subfield` are specifically given as a sequence.

    """
    # Follow any links in _KIND_DICT
    while isinstance(_KIND_DICT[kind], basestring):
        kind = _KIND_DICT[kind]

    # List of fields
    single_field = False
    # kinds and suffixes of kinds without fields
    kinds = ['weights']
    suffixes = ['.mosaic', '.montage', '.hdr', '.bg']
    if kind in kinds or any(kind.endswith(suffix) for suffix in suffixes):
        field = [None]
        single_field = True
    else:
        if 'galex' in kind:
            allfields = GALEX_TILE_LIST
        else:
            allfields = BRICK_LIST
        field =  kwargs.get('field')
        if field is None:
            field = allfields
        else:
            # Only fields in allfields are available
            try:  # field is a list of field numbers
                field = [fld for fld in field if fld in allfields]
            except TypeError:  # field is a single field number
                field = [field] if field in allfields else []
                single_field = True

    # List of subfields
    single_subfield = False
    subfield_kinds = ['phot', 'sfh', 'cmd', 'bestzcb']  # kinds that have subfields
    if kind in subfield_kinds:
        allsubfields = CELL_LIST
        subfield = kwargs.get('subfield')
        if subfield is None:
            subfield = allsubfields
        else:
            # Only subfields in allsubfields are available
            try:  # subfield is a list of subfields
                subfield = [sfld for sfld in subfield if sfld in allsubfields]
            except TypeError:  # subfield is a single subfield number
                subfield = [subfield] if subfield in allsubfields else []
                single_subfield = True
    else:
        subfield = [None]
        single_subfield = True

    # List of paths
    missing = kwargs.get('missing', [])
    path_list = []
    pth = os.path.join(*_KIND_DICT[kind])
    for fld in field:
        for sfld in subfield:
            if (fld, sfld) in missing:
                if kwargs.get('fillsubfield', True):
                    path_list.append(None)
                else:
                    continue
            else:
                av, dav = extpar_dict.get((fld, sfld), (None, None))
                path_list.append(pth.format(fld, sfld, av, dav))

    # Return None if no paths; don't return a list of length 1 unless a
    # list was specifically given for field or subfield
    if not path_list:
        path_list = None
    elif len(path_list) == 1 and single_field and single_subfield:
        path_list = path_list[0]

    return path_list



# Parsers
# =======

def open_extparfile(filename):
    """Create a table from an 'extpar' file.

    An 'extpar' file contains the best-fit Av and dAv extinction parameters
    for all cells in a single brick. There are three columns: cell
    name/number, Av, and dAv. There is one row per cell, and the rows are
    sorted by cell name in numerical order. Cells without SFH data have
    their number, Av, and dAv values set to 0, 99, and 99, respectively.

    Parameters
    ----------
    filename : str
        Path to the file.

    Returns
    -------
    astropy.table.Table
        See Notes for the columns.

    Notes
    -----
    The columns in the output table are,

    ====== ====================================================
    column description
    ====== ====================================================
    cell   Cell number
    av     Av (foreground V-band) extinction parameter value
    dav    dAv (differential V-band) extinction parameter value
    ====== ====================================================

    """
    table = astropy.table.Table.read(filename, format='ascii')
    table['col1'].name = 'cell'
    table['col2'].name = 'av'
    table['col3'].name = 'dav'
    return table


def open_cornersfile(filename):
    """Create a table from a 'corners' file.

    A 'corners' file contains eight columns for the RA and dec coordinates
    of each corner of each cell in a single brick. The corners are numbered
    clockwise with 1 at the xmin,ymax corner and 4 at the xmin,ymin corner.
    There is one row per cell, and the rows are sorted by cell name/number
    in numerical order (although cell name is not actually given as a
    column in the file).

    Parameters
    ----------
    filename : str
        Path to the file.

    Returns
    -------
    astropy.table.Table
        See Notes for the columns.

    Notes
    -----
    The columns in the output table are,

    ====== ===============
    column description
    ====== ===============
    RA1    RA of corner 1
    dec1   dec of corner 1
    RA2    RA of corner 2
    dec2   dec of corner 2
    RA3    RA of corner 3
    dec3   dec of corner 3
    RA4    RA of corner 4
    dec4   dec of corner 4
    ====== ===============

    """
    table = astropy.table.Table.read(filename, format='ascii')
    table['col1'].name = 'RA1'
    table['col2'].name = 'dec1'
    table['col3'].name = 'RA2'
    table['col4'].name = 'dec2'
    table['col5'].name = 'RA3'
    table['col6'].name = 'dec3'
    table['col7'].name = 'RA4'
    table['col8'].name = 'dec4'
    return table


def cornergrid(table):
    """Reshape a table of cell corner coordinates into a grid.

    Parameters
    ----------
    table : ndarray
        Table from `open_cornersfile` containing the RA,dec coordinates of
        the corners of cells in a given brick. It is assumed that the rows
        are sorted by cell name/number in numerical order.

    Returns
    -------
    ((NROW+1),(NCOL+1)) ndarray
        RA coordinates of the cell corners in the brick grid, i.e., all
        combinations of i and j integers.
    ((NROW+1),(NCOL+1)) ndarray
        dec coordinates of the cell corners.

    """
    table = table[SORT]
    shape = NROW, NCOL

    agrid = np.zeros((NROW+1, NCOL+1))
    agrid[1:,:-1] = table['RA1'].reshape(shape)
    agrid[1:,1:] = table['RA2'].reshape(shape)
    agrid[:-1,1:] = table['RA3'].reshape(shape)
    agrid[:-1,:-1] = table['RA4'].reshape(shape)

    dgrid = np.zeros((NROW+1, NCOL+1))
    dgrid[1:,:-1] = table['dec1'].reshape(shape)
    dgrid[1:,1:] = table['dec2'].reshape(shape)
    dgrid[:-1,1:] = table['dec3'].reshape(shape)
    dgrid[:-1,:-1] = table['dec4'].reshape(shape)

    return agrid, dgrid



# Run at import
# =============


# Extinction parameter dictionary
# -------------------------------
# The bottleneck limiting the import speed of this module is in
# `_load_extpar_dict`, probably in creating the `Table` instances.

def _load_extpar_dict():
    """Load Av and dAv extinction parameter values for all cells.

    Loop through all bricks, open their 'extpar' files, and save the Av and
    dAv values for all cells as (Av, dAv) tuples in a dictionary keyed by
    (brick, cell). Cells without SFH data have their cell number, Av, and
    dAv values set to 0, 99, and 99, respectively.

    Returns
    -------
    dict

    """
    pth = os.path.join(*_KIND_DICT['extpar'])
    extpar_dict = {}
    for brick in BRICK_LIST:
        extparfile = pth.format(brick)
        table = open_extparfile(extparfile)
        # Have to iterate over CELL_LIST -- table['cell'] is set to 0 where
        # there is no data instead of giving the actual cell number.
        extpar_dict.update({(brick, cell): (row['av'], row['dav'])
                            for cell, row in zip(CELL_LIST, table)})
    return extpar_dict

EXTPAR_DICT = _load_extpar_dict()
"""Dictionary of the best-fit Av and dAv extinction parameters for every cell.

The values are (Av, dAv) tuples and are keyed by (brick, cell). Both Av and
dAv are equal to 99 if the cell does not actually have any SFH data.

"""


# List of cells without data
# --------------------------

def _find_missing_cells(extpar_dict, missing=None):
    """Find specific cells for which no data are available.

    If a cell did not have SFH data available when its brick's 'extpar'
    file was written, the Av and dAv values were recorded as 99. Cells
    without data can therefore be identified by searching `extpar_dict` for
    Av and dAv values equal to 99. Any such cells are added to a list as a
    (brick, cell) tuple.

    Parameters
    ----------
    extpar_dict : dict
        Extinction parameters as (Av, dAv) tuples keyed by (brick, cell).
    missing : list, optional
        Existing list of missing cells. Missing items are appended directly
        (i.e., the original list is modified). If None (default), a new
        list is created.

    Returns
    -------
    list
        (brick, cell) tuples for where there is no SFH data.

    """
    badval = 99
    missing = [] if missing is None else missing
    missing += [brickcell for brickcell, avdav in extpar_dict.items()
                if avdav == (badval, badval) and brickcell not in missing]
    return missing

MISSING_CELLS = _find_missing_cells(EXTPAR_DICT)
"""List of (brick, cell) tuples where there is no data."""


# Wrap `_path`
# ------------
# Create a wrapper around `_path` that specifically uses `EXTPAR_DICT` and
# `MISSING_CELLS`. Also inherit the original docstring with some
# modifications.

def __experimental__inheritdocstr(func):
    """Automatic parameter filtering.

    .. note:: Still need a way to retain the newline before the 'Returns'
       section when filtering out the last parameter.

    """
    indent = '    '
    param_list = ['extpar_dict']  # skip these parameters

    docstr1 = _path.__doc__
    lines1 = docstr1.split('\n{:s}'.format(indent))
    lines2 = []
    for line in lines1:
        if ' : ' in line:
            param = line.split(' : ')[0]
        elif line and not line[0].isspace():
            param = ''
        if param in param_list:
            continue
        else:
            lines2.append(line)
    docstr2 = '\n{:s}'.format(indent).join(lines2)

    func.__doc__ = docstr2
    return func

def _inheritdocstr(func):
    """Manually select pieces of docstring to inherit.

    Make sure to update this when the `_path` docstring is updated!

    """
    docstr1 = _path.__doc__
    lines1 = docstr1.split('\n')
    lines2 = lines1[:26] + lines1[30:48] + lines1[52:]
    docstr2 = '\n'.join(lines2)
    func.__doc__ = docstr2
    return func

@_inheritdocstr
def path(kind, **kwargs):
    """Wrap `_path` to explicitly use `EXTPAR_DICT` and `MISSING_CELLS`."""
    return _path(kind, EXTPAR_DICT, missing=MISSING_CELLS, **kwargs)
