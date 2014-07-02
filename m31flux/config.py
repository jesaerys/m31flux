"""

=====================
`m31flux.config`
=====================

Configuration for `m31flux`.

This module contains constants and functions that locate data (both
existing files and files created during analysis) and define how it is
structured. The main features are,

1. A `path` function, which keeps track of all project files. Paths to
   particular files can be changed and new files can be added by modifying the
   `_kind_dict`.

2. Other project-specific things like galaxy and modeling parameters and
   parsers for nonstandard files.

3. A description of how the PHAT bricks are gridded and a model for
   mapping the data to the grids. See `Brick grid model`_ below.


Constants
=========


Functions
=========


Brick grid model
================

There are 23 bricks in the PHAT survey, labeled by integer starting with 1.
Each brick is divided into an n*m grid, and the grid cells are labeled by
integer from 1 to n*m. Viewing a brick with north and west approximately up
and to the right, respectively, the cells are arranged according to the
diagram::

    n+0.5     +-----------+-----------+-----------+-----------+
    n      n-1|         1 |         2 |    ...    |      m    |
    n-0.5     +-----------+-----------+-----------+-----------+
    n-1    n-2|     1*m+1 |     1*m+2 |    ...    |    2*m    |      ^
    n-1.5     +-----------+-----------+-----------+-----------+      N
    n-2    n-3|     2*m+1 |     2*m+2 |    ...    |    3*m    |  < E   W >
    n-2.5     +-----------+-----------+-----------+-----------+      S
    ...    ...|    ...    |    ...    |    ...    |    ...    |      v
    1.5       +-----------+-----------+-----------+-----------+
    1      i=0| (n-1)*m+1 | (n-1)*m+2 |    ...    |    n*m    |
  y=0.5       +-----------+-----------+-----------+-----------+
                  j=0           1          ...         m-1
           x=0.5    1    1.5    2    2.5   ...   m-0.5   m   m+0.5

The convention defined in this module is to set the grid origin at cell
number ``(n-1)*m+1`` with the rows (i) increasing upward and the columns
(j) increasing to the right in the diagram. In addition, a pixel coordinate
system is assigned such that the center of the origin cell has pixel
coordinates x,y = 1,1, with x and y increasing with j and i, respectively.

*The cell numbers do not reflect their actual order in the grid array under
this convention!* This is can be confusing because cell 1 is not the first
cell. Rather, in a flattened n*m array (e.g., `numpy.ravel`), cell
``(n-1)*m+1`` is first and cell m is last. The cell numbers are merely
labels that happen to be integers.



###
SFHs were derived for the cells of all bricks listed in `BRICK_LIST`,
except for any cells listed in `MISSING_CELLS` (see `_find_missing_cells`).
The SFHs were calculated using MATCH, assuming the Av,dAv extinction model.
The `EXTPAR_DICT` dictionary contains the extinction parameter values for
each cell (see also `_load_extpar_dict`).

The `path` function locates all data in the project so that no path names
have to be hardcoded in any modules or analyis scripts.
###

"""
import astrogrid
import astropy.table
import numpy as np
import os



# PHAT bricks and the brick grid system
# =====================================
BRICK_LIST = np.array([2]+range(4,24))
"""List of all PHAT bricks available in the project in numerical order."""

NROW = 15
"""Number of rows in a brick grid."""

NCOL = 30
"""Number of columns in a brick grid."""

CELL_LIST = np.arange(NROW*NCOL) + 1
"""List of all cells in a brick grid in numerical order."""

SORT = np.arange(NROW*NCOL).reshape((NROW, NCOL))[::-1].ravel()
"""Indices that sort a list of cells in numerical order to grid order, and
vice versa.

"Grid order" refers to the order of the cells in a flattened array (e.g.,
`numpy.ravel`). See the main documentation of `m31flux.config` for a
description of brick grid coordinate conventions.

"""



# GALEX
# =====
GALEX_FIELD_LIST = [0, 7, 8, 9, 10]
"""List of all GALEX DIS fields (or tiles) covering the PHAT survey area,
in numerical order.

"""

#GALEX_FUV_FOREGROUND_CPS = 8.3573867e-4  # m31flux.util._calc_galex_foreground
#GALEX_NUV_FOREGROUND_CPS = 5.393513e-3  # m31flux.util._calc_galex_foreground

GALEX_CHIP_X0 = 1920
"""Approximate x pixel coordinate of chip center."""

GALEX_CHIP_Y0 = 1920
"""Approximate y pixel coordinate of chip center."""

GALEX_CHIP_RAD = 1400
"""Approximate chip radius in pixels, slightly undersized."""



# M31 and other parameters
# ========================
DMOD = 24.47
"""Distance modulus.

McConnachie, A. W., Irwin, M. J., Ferguson, A. M. N., et al. 2005, MNRAS,
356, 979.

"""

DIST_PC = 10**(DMOD/5. + 1)
"""Distance in pc."""

DIST_CM = DIST_PC * 3.08567758e18
"""Distance in cm."""

IMF = 'Kroupa'
"""Initial mass function."""



# Paths
# =====
_PROJECT_DIR = '/Users/Jake/Research/PHAT/sfhmaps-phat'
_SFH_DIR = os.path.join(_PROJECT_DIR, 'sfh')
_ANALYSIS_DIR = os.path.join(_PROJECT_DIR, 'analysis')
_GALEX_DIR = '/Users/Jake/Research/Storage/M31/GALEX/DIS'


_kind_dict = {
    'extpar':  (_SFH_DIR, 'B{0:02d}', 'b{0:02d}_region_AvdAv.dat'),
    'corners': (_SFH_DIR, 'B{0:02d}', 'M31-B{0:02d}_15x30_subregion-exact-vertices.dat'),
    'phot':    (_SFH_DIR, 'B{0:02d}', 'phot', 'M31-B{0:02d}_15x30-{1:03d}.gst.match'),
    'sfh':     (_SFH_DIR, 'B{0:02d}', 'sfh',  'M31-B{0:02d}_15x30-{1:03d}_{2:.1f}-{3:.1f}_best.sfh'),
    'cmd':     (_SFH_DIR, 'B{0:02d}', 'cmd',  'M31-B{0:02d}_15x30-{1:03d}_{2:.1f}-{3:.1f}_best.sfh.cmd'),
    'bestzcb': (_ANALYSIS_DIR, 'b{0:02d}', 'bestzcb', 'b{0:02d}-{1:03d}_best.zcb'),

    'mod_fuv_int': (_ANALYSIS_DIR, 'b{0:02d}', 'b{0:02d}_mod_fuv_int.fits'),
    'mod_fuv_int.add': (_ANALYSIS_DIR, 'mod_fuv_int.fits'),
    'mod_fuv_int.density': (_ANALYSIS_DIR, '_mod_fuv_int', 'input', 'b{0:02d}_mod_fuv_int_density.fits'),
    'mod_fuv_int.reproject': (_ANALYSIS_DIR, '_mod_fuv_int', 'reproject', 'hdu0_b{0:02d}_mod_fuv_int_density.fits'),
    'mod_fuv_int.reproject.add': (_ANALYSIS_DIR, '_mod_fuv_int', 'add', 'mod_fuv_int_density.fits'),
    'mod_fuv_int.area': (_ANALYSIS_DIR, '_mod_fuv_int', 'reproject', 'hdu0_b{0:02d}_mod_fuv_int_density_area.fits'),
    'mod_fuv_int.area.add': (_ANALYSIS_DIR, '_mod_fuv_int', 'add', 'mod_fuv_int_density_area.fits'),
    'mod_fuv_int.hdr': 'mod_fuv_red.hdr',

    'mod_fuv_red': (_ANALYSIS_DIR, 'b{0:02d}', 'b{0:02d}_mod_fuv_red.fits'),
    'mod_fuv_red.add': (_ANALYSIS_DIR, 'mod_fuv_red.fits'),
    'mod_fuv_red.density': (_ANALYSIS_DIR, '_mod_fuv_red', 'input', 'b{0:02d}_mod_fuv_red_density.fits'),
    'mod_fuv_red.reproject': (_ANALYSIS_DIR, '_mod_fuv_red', 'reproject', 'hdu0_b{0:02d}_mod_fuv_red_density.fits'),
    'mod_fuv_red.reproject.add': (_ANALYSIS_DIR, '_mod_fuv_red', 'add', 'mod_fuv_red_density.fits'),
    'mod_fuv_red.area': (_ANALYSIS_DIR, '_mod_fuv_red', 'reproject', 'hdu0_b{0:02d}_mod_fuv_red_density_area.fits'),
    'mod_fuv_red.area.add': (_ANALYSIS_DIR, '_mod_fuv_red', 'add', 'mod_fuv_red_density_area.fits'),
    'mod_fuv_red.hdr': (_ANALYSIS_DIR, '_mod_fuv_red', 'template.hdr'),

    'galex_fuv': (_GALEX_DIR, 'PS_M31_MOS{0:02d}-fd-int.fits'),
    'galex_fuv.add': (_ANALYSIS_DIR, 'galex_fuv.fits'),
    'galex_fuv.density': (_ANALYSIS_DIR, '_galex_fuv', 'input', 'MOS{0:02d}_galex_fuv_density.fits'),
    'galex_fuv.reproject': (_ANALYSIS_DIR, '_galex_fuv', 'reproject', 'hdu0_MOS{0:02d}_galex_fuv_density.fits'),
    'galex_fuv.reproject.add': (_ANALYSIS_DIR, '_galex_fuv', 'add', 'galex_fuv_density.fits'),
    'galex_fuv.area': (_ANALYSIS_DIR, '_galex_fuv', 'reproject', 'hdu0_MOS{0:02d}_galex_fuv_density_area.fits'),
    'galex_fuv.area.add': (_ANALYSIS_DIR, '_galex_fuv', 'add', 'galex_fuv_density_area.fits'),
    'galex_fuv.hdr': 'mod_fuv_red.hdr',

    'mod_nuv_int': (_ANALYSIS_DIR, 'b{0:02d}', 'b{0:02d}_mod_nuv_int.fits'),
    'mod_nuv_int.add': (_ANALYSIS_DIR, 'mod_nuv_int.fits'),
    'mod_nuv_int.density': (_ANALYSIS_DIR, '_mod_nuv_int', 'input', 'b{0:02d}_mod_nuv_int_density.fits'),
    'mod_nuv_int.reproject': (_ANALYSIS_DIR, '_mod_nuv_int', 'reproject', 'hdu0_b{0:02d}_mod_nuv_int_density.fits'),
    'mod_nuv_int.reproject.add': (_ANALYSIS_DIR, '_mod_nuv_int', 'add', 'mod_nuv_int_density.fits'),
    'mod_nuv_int.area': (_ANALYSIS_DIR, '_mod_nuv_int', 'reproject', 'hdu0_b{0:02d}_mod_nuv_int_density_area.fits'),
    'mod_nuv_int.area.add': (_ANALYSIS_DIR, '_mod_nuv_int', 'add', 'mod_nuv_int_density_area.fits'),
    'mod_nuv_int.hdr': 'mod_fuv_red.hdr',

    'mod_nuv_red': (_ANALYSIS_DIR, 'b{0:02d}', 'b{0:02d}_mod_nuv_red.fits'),
    'mod_nuv_red.add': (_ANALYSIS_DIR, 'mod_nuv_red.fits'),
    'mod_nuv_red.density': (_ANALYSIS_DIR, '_mod_nuv_red', 'input', 'b{0:02d}_mod_nuv_red_density.fits'),
    'mod_nuv_red.reproject': (_ANALYSIS_DIR, '_mod_nuv_red', 'reproject', 'hdu0_b{0:02d}_mod_nuv_red_density.fits'),
    'mod_nuv_red.reproject.add': (_ANALYSIS_DIR, '_mod_nuv_red', 'add', 'mod_nuv_red_density.fits'),
    'mod_nuv_red.area': (_ANALYSIS_DIR, '_mod_nuv_red', 'reproject', 'hdu0_b{0:02d}_mod_nuv_red_density_area.fits'),
    'mod_nuv_red.area.add': (_ANALYSIS_DIR, '_mod_nuv_red', 'add', 'mod_nuv_red_density_area.fits'),
    'mod_nuv_red.hdr': 'mod_fuv_red.hdr',

    'galex_nuv': (_GALEX_DIR, 'PS_M31_MOS{0:02d}-nd-int.fits'),
    'galex_nuv.add': (_ANALYSIS_DIR, 'galex_nuv.fits'),
    'galex_nuv.density': (_ANALYSIS_DIR, '_galex_nuv', 'input', 'MOS{0:02d}_galex_nuv_density.fits'),
    'galex_nuv.reproject': (_ANALYSIS_DIR, '_galex_nuv', 'reproject', 'hdu0_MOS{0:02d}_galex_nuv_density.fits'),
    'galex_nuv.reproject.add': (_ANALYSIS_DIR, '_galex_nuv', 'add', 'galex_nuv_density.fits'),
    'galex_nuv.area': (_ANALYSIS_DIR, '_galex_nuv', 'reproject', 'hdu0_MOS{0:02d}_galex_nuv_density_area.fits'),
    'galex_nuv.area.add': (_ANALYSIS_DIR, '_galex_nuv', 'add', 'galex_nuv_density_area.fits'),
    'galex_nuv.hdr': 'mod_fuv_red.hdr'
    }
"""Path elements used to build paths to the different file kinds.

This dictionary is used by `_path` and defines the locations of all file
kinds in the project. Each file kind is specified by a tuple of path
elements (strings) which are later combined into a single string using
`os.path.join`.

The elements may contain format strings as placeholders for inserting one
or more of the values in a standard tuple, ``(field, subfield, av, dav)``
when `_path` is actually called. The first two values are documented in
`_path`, and the last two refer to Av and dAv extinction parameters. The
values are referenced in a format string according to their tuple index.
For example, putting ``'{1:02d}'`` in a path element will insert the value
of ``subfield`` formatted as a padded two-column integer.

To link one file kind to another, set the kind's value to the name of the
other kind. This is useful if two kinds actually share the same file.

"""


def _path(kind, extpar_dict, **kwargs):
    """File paths for the `m31flux` project.

    Parameters
    ----------
    kind : str
        The type of file. The exact set of files returned is controlled
        using keyword arguments. See the Notes section below.
    extpar_dict : dictionary
        Dictionary of brick grid extinction parameter values as (av, dav)
        tuples, keyed by (brick, cell) tuples. This dictionary is mandatory
        as some file names contain av and dav parameter values.
    field, subfield : int or list of ints, optional
        Return the file(s) of type `kind` for a combination of fields and
        subfields. None (default) is equivalent to a list of all fields or
        all subfields. "field" and "subfield" are generic terms that
        represent more specific organizational terms for different kinds of
        files:

        ================ ===== ========
        kind category    field subfield
        ================ ===== ========
        PHAT             brick field
        PHAT brick grids brick cell
        GALEX            tile
        ================ ===== ========

    missing : list, optional
        List of (brick, cell) tuples for which no SFH data are available
        (useful when automatically generating paths so that only valid
        paths are returned).
    fillsubfield : bool, optional
        If True (default), any subfields listed in `missing` are added to
        the returned path list as None. If False, all missing subfields are
        omitted from the path list.

    Returns
    -------
    string or list of strings
        List of file paths.

    Notes
    -----
    Possible values for the `kind` parameter are:

    corners
        RA,dec coordinates of the corners of each cell in a brick. See
        `open_cornersfile` for details.
    extpar
        Best-fit Av and dAv extinction parameters for each cell in a brick.
        See `open_extparfile` for details.
    phot, sfh, cmd, bestzcb
        MATCH-related files for a given cell, produced by the following
        procedure (kind names are in brackets; refer to the MATCH README
        for further details about these file types)::

          calcsfh par [phot] fake [sfh] -zinc -mcdata -dAvy=0 -dAv=dAv  # other output: [cmd], hmcdat
          zcombine [sfh] -bestonly > [bestzcb]

        .. note:: In case additional kinds are supported in the future,
           should probably stick to a naming scheme similar to this::

             hybridMC hmcdat [hmcsfh] -tint=2.0 -nmc=10000 -dt=0.015
             zcombine -unweighted -medbest -jeffreys -best=[bestzcb] [hmcsfh] > [hmczcb]
             zcmerge [bestzcb] [hmczcb] -absolute > [besthmczcb]

             zcombine [extsyssfh]* > [extsyszcb]  # extinction systematics
             zcmerge [bestzcb] [extsyszcb] -absolute > [bestextzcb]

             zcombine [isosyssfh]* > [isosyszcb]  # isochrone systematics
             zcmerge [bestzcb] [isosyszcb] -absolute > [bestisozcb]

    <band>
        Field image.
    <band>.add
        Final mosaic for the header <band>.hdr. Product of
        <band>.reproject.add and <band>.area.add.
    <band>.density
        Field density image (signal per unit pixel area).
    <band>.reproject, <band>.area
        Field density image reprojected to the header <band>.hdr, and the
        corresponding area file.
    <band>.reproject.add, <band>.area.add
        Density image mosaic for the header <band>.hdr, and the
        corresponding area file.
    <band>.hdr
        The WCS header template for the mosaic.

    where <band> may be one of the following:

    mod_fuv_int
        Intrinsic FUV flux modeled from the SFHs.
    mod_fuv_red
        Reddened FUV flux modeled from the SFHs and the Av,dAv extinction
        model.
    galex_fuv
        Observed FUV flux from GALEX.
    mod_nuv_int
        Intrinsic NUV flux modeled from the SFHs.
    mod_nuv_red
        Reddened NUV flux modeled from the SFHs and the Av,dAv extinction
        model.
    galex_nuv
        Observed NUV flux from GALEX.

    """
    # Follow any links in _kind_dict
    while not astrogrid.util.islistlike(_kind_dict[kind]):
        kind = _kind_dict[kind]

    # List of fields
    if kind.endswith('.add') or kind.endswith('.hdr'):
        field = [None]
    else:
        if 'galex' in kind:
            allfields = GALEX_FIELD_LIST
        else:
            allfields = BRICK_LIST
        field =  kwargs.get('field')
        if field is None:
            field = allfields
        elif not astrogrid.util.islistlike(field):
            field = [field]  # convert to list
        # Only fields in allfields are available
        field = [fld for fld in field if fld in allfields]

    # List of subfields
    subfield_kinds = ['phot', 'sfh', 'cmd', 'bestzcb']  # kinds that have subfields
    if kind in subfield_kinds:
        allsubfields = CELL_LIST
        subfield = kwargs.get('subfield')
        if subfield is None:
            subfield = allsubfields
        elif not astrogrid.util.islistlike(subfield):
            subfield = [subfield]  # convert to list
        # Only subfields in allsubfields are available
        subfield = [sfld for sfld in subfield if sfld in allsubfields]
    else:
        subfield = [None]

    # List of paths
    missing = kwargs.get('missing', [])
    path_list = []
    pth = os.path.join(*_kind_dict[kind])
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
    cond = len(path_list) == 1 and not astrogrid.util.islistlike(kwargs.get('field'))
    if not path_list:
        path_list = None
    elif (cond and not astrogrid.util.islistlike(kwargs.get('subfield'))
            if kind in subfield_kinds else cond):
        path_list = path_list[0]

    return path_list



# Parsers
# =======
def open_extparfile(filename):
    """Create a table from a file of kind 'extpar'.

    An 'extpar' file contains three columns: the cell number and the Av and
    dAv parameter values for the cell. There is one row per cell, and the
    cells are listed in increasing order consistent with output from
    `grid2list`. Cells without SFH data have Av and dAv both set to 99.

    Returns an `astropy.table.Table` instance with columns 'cell', 'av',
    and 'dav'.

    """
    table = astropy.table.Table.read(filename, format='ascii')
    table['col1'].name = 'cell'
    table['col2'].name = 'av'
    table['col3'].name = 'dav'
    return table


def open_cornersfile(filename):
    """Create a table from a file of kind 'corners'.

    A 'corners' file contains eight columns for the RA and dec coordinates
    of each corner of each cell in a brick. The corners are numbered
    clockwise with 1 at the xmin,ymax corner and 4 at the xmin,ymin corner.
    There is one row per cell, and the cells are listed in increasing order
    consistent with output from `grid2list`.

    Returns an `astropy.table.Table` instance with columns 'RA1', 'dec1',
    'RA2', 'dec2', 'RA3', 'dec3', 'RA4', and 'dec4'.

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


def cornergrid(arr, shape=(NROW, NCOL)):
    """Reshape a list of cell corners into a grid.

    Parameters
    ----------
    arr : array
        1d array of corner values from `open_cornersfile`.
    shape : tuple, optional
        The shape of the output grid. Default is (`NROWS`, `NCOLS`).

    Returns
    -------
    tuple
        A tuple of two 2d arrays of the given shape: one for RA, the other
        for dec.

    """
    ny, nx = shape

    agrid = np.zeros((ny+1, nx+1))
    agrid[1:,:-1] = arr['RA1'].reshape(shape)[::-1]
    agrid[1:,1:] = arr['RA2'].reshape(shape)[::-1]
    agrid[:-1,1:] = arr['RA3'].reshape(shape)[::-1]
    agrid[:-1,:-1] = arr['RA4'].reshape(shape)[::-1]

    dgrid = np.zeros((ny+1, nx+1))
    dgrid[1:,:-1] = arr['dec1'].reshape(shape)[::-1]
    dgrid[1:,1:] = arr['dec2'].reshape(shape)[::-1]
    dgrid[:-1,1:] = arr['dec3'].reshape(shape)[::-1]
    dgrid[:-1,:-1] = arr['dec4'].reshape(shape)[::-1]

    return agrid, dgrid



# Run at import
# =============


# Extinction parameter dictionary
# -------------------------------
# The bottleneck limiting import speed is in _load_extpar_dict, probably in
# creating the `Table` instances.
def _load_extpar_dict():
    """Load av and dav extinction parameter values for all cells.

    Loop through all bricks, open their 'extpar' files, and save the av and
    dav values for all cells as (av, dav) tuples in dictionary keyed by
    (brick, cell). Cells without SFH data have av and dav both set to 99.

    """
    pth = os.path.join(*_kind_dict['extpar'])
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


# List of cells without data
# --------------------------
def _find_missing_cells(extpar_dict, missing=None, badval=99):
    """Find specific cells for which no data are available.

    If a cell did not have SFH data available when its brick's 'extpar'
    file was written, the av and dav values were recorded as `badval`.
    Cells without data can therefore be identified by searching
    `extpar_dict` for av and dav values equal to `badval`. Any such cells
    are added to a list as a (brick, cell) tuple.

    Parameters
    ----------
    extpar_dict : dictionary
        Dictionary of extinction parameters as (av, dav) tuples, keyed by
        (brick, cell).
    missing : list, optional
        Existing list of missing cells. Missing items are appended
        directly (i.e., the list is modified). If None (default), a new
        list is created.
    badval : int or float, optional
        Value of av or dav which signifies that a given cell has no data.
        Default is 99.

    Returns
    -------
    list
        List of (brick, cell) tuples that have no data.

    """
    if missing is None:
        missing = []

    missing += [brickcell for brickcell, avdav in extpar_dict.items()
                if avdav == (badval, badval) and brickcell not in missing]

    return missing

MISSING_CELLS = _find_missing_cells(EXTPAR_DICT)


# Wrap `_path`
# ------------
# Want a wrapper around `_path` that specifically uses `EXTPAR_DICT` and
# `MISSING_CELLS`. Also want to inherit the original docstring with some
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
    lines2 = lines1[:7] + lines1[11:26] + lines1[30:]
    docstr2 = '\n'.join(lines2)
    func.__doc__ = docstr2
    return func

@_inheritdocstr
def path(kind, **kwargs):
    """Wrap `_path` to explicitly use `EXTPAR_DICT` and `MISSING_CELLS`."""
    return _path(kind, EXTPAR_DICT, missing=MISSING_CELLS, **kwargs)
