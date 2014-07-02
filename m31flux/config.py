"""

=====================
`sfhmaps_phat.config`
=====================

Constants and functions that locate data (both where it is and where it
will be) and define how it is structured.


Brick grid model
----------------
There are 23 bricks in the PHAT survey, numbered by integer starting with
1. Each brick is divided into a grid with ny (`NROW`) rows and nx (`NCOL`)
columns. Viewing a brick with north and west approximately up and to the
right, respectively, the cells in the grid are numbered from 1 to ny*nx
(`CELL_LIST`) according to the diagram::

   +---------------------------------+
   |          1           2 ...    nx|      ^
  ^|     1*nx+1      1*nx+2 ...  2*nx|      N
   |     2*nx+1      2*nx+2 ...  3*nx|  < E   W >
  y|        ...         ... ...   ...|      S
   |(ny-1)*nx+1 (ny-1)*nx+2 ... ny*nx|      v
   +---------------------------------+
  0    x   >

The grid is attached to a pixel coordinate system with the x,y origin at
cell number ``(ny-1)*nx+1`` so that rows increase upward and columns
increase to the right in the diagram. The center of the origin cell has
pixel coordinates (x, y) = (1, 1) and the outer corner is at (0.5, 0.5).

SFHs were derived for the cells of all bricks listed in `BRICK_LIST`,
except for any cells listed in `MISSING_CELLS` (see `_find_missing_cells`).
The SFHs were calculated using MATCH, assuming the Av,dAv extinction model.
The `EXTPAR_DICT` dictionary contains the extinction parameter values for
each cell (see also `_load_extpar_dict`).

The `path` function locates all data in the project so that no path names
have to be hardcoded in any modules or analyis scripts.

"""
import astrogrid
import astropy.table
import numpy as np
import os



# PHAT bricks and the brick grid system
# =====================================
BRICK_LIST = np.array([2]+range(4,24))
NROW, NCOL = 15, 30
# This list of indices will 
REORDER = np.arange(NROW*NCOL).reshape((NROW, NCOL))[::-1].ravel()
# List of cells, ordered as a flattened array for astrogrid.Grid
CELL_LIST = np.arange(NROW*NCOL)[REORDER] + 1



# GALEX
# =====
GALEX_FIELD_LIST = [0, 7, 8, 9, 10]
GALEX_FUV_FOREGROUND_CPS = 8.3573867e-4  # sfhmaps_phat.util._calc_galex_foreground
GALEX_NUV_FOREGROUND_CPS = 5.393513e-3  # sfhmaps_phat.util._calc_galex_foreground
GALEX_CHIP_X0 = 1920  # Approximate x pixel coordinate of chip center
GALEX_CHIP_Y0 = 1920  # Approximate y pixel coordinate of chip center
GALEX_CHIP_RAD = 1400  # Approximate chip radius in pixels, slightly undersized



# M31 and other parameters
# ========================
DMOD = 24.47  # Distance modulus; McConnachie, A. W., Irwin, M. J., Ferguson, A. M. N., et al. 2005, MNRAS, 356, 979
DIST_PC = 10**(DMOD/5. + 1)  # Distance in pc
DIST_CM = DIST_PC * 3.08567758e18  # Distance in cm
IMF = 'Kroupa'
IMF_TYPE = 2  # FSPS IMF code



# Paths
# =====
PROJECT_DIR = '/Users/Jake/Research/PHAT/sfhmaps-phat'
SFH_DIR = os.path.join(PROJECT_DIR, 'sfh')
ANALYSIS_DIR = os.path.join(PROJECT_DIR, 'analysis')
GALEX_DIR = '/Users/Jake/Research/Storage/M31/GALEX/DIS'


"""Each file kind is specified by a tuple of path elements which are later
combined using `os.path.join`. All paths are formatted using a tuple of
values, (field, subfield, av, dav). Values are selected and placed in a
path string according to their indices (a value is omitted simply by not
using its index in any format strings in a path).

"""
_kind_dict = {
    'extpar':  (SFH_DIR, 'B{0:02d}', 'b{0:02d}_region_AvdAv.dat'),
    'corners': (SFH_DIR, 'B{0:02d}', 'M31-B{0:02d}_15x30_subregion-exact-vertices.dat'),
    'phot':    (SFH_DIR, 'B{0:02d}', 'phot', 'M31-B{0:02d}_15x30-{1:03d}.gst.match'),
    'sfh':     (SFH_DIR, 'B{0:02d}', 'sfh',  'M31-B{0:02d}_15x30-{1:03d}_{2:.1f}-{3:.1f}_best.sfh'),
    'cmd':     (SFH_DIR, 'B{0:02d}', 'cmd',  'M31-B{0:02d}_15x30-{1:03d}_{2:.1f}-{3:.1f}_best.sfh.cmd'),
    'bestzcb': (ANALYSIS_DIR, 'b{0:02d}', 'bestzcb', 'b{0:02d}-{1:03d}_best.zcb'),

    'mod_fuv_int': (ANALYSIS_DIR, 'b{0:02d}', 'b{0:02d}_mod_fuv_int.fits'),
    'mod_fuv_int.add': (ANALYSIS_DIR, 'mod_fuv_int.fits'),
    'mod_fuv_int.density': (ANALYSIS_DIR, '_mod_fuv_int', 'input', 'b{0:02d}_mod_fuv_int_density.fits'),
    'mod_fuv_int.reproject': (ANALYSIS_DIR, '_mod_fuv_int', 'reproject', 'hdu0_b{0:02d}_mod_fuv_int_density.fits'),
    'mod_fuv_int.reproject.add': (ANALYSIS_DIR, '_mod_fuv_int', 'add', 'mod_fuv_int_density.fits'),
    'mod_fuv_int.area': (ANALYSIS_DIR, '_mod_fuv_int', 'reproject', 'hdu0_b{0:02d}_mod_fuv_int_density_area.fits'),
    'mod_fuv_int.area.add': (ANALYSIS_DIR, '_mod_fuv_int', 'add', 'mod_fuv_int_density_area.fits'),
    'mod_fuv_int.hdr': 'mod_fuv_red.hdr',

    'mod_fuv_red': (ANALYSIS_DIR, 'b{0:02d}', 'b{0:02d}_mod_fuv_red.fits'),
    'mod_fuv_red.add': (ANALYSIS_DIR, 'mod_fuv_red.fits'),
    'mod_fuv_red.density': (ANALYSIS_DIR, '_mod_fuv_red', 'input', 'b{0:02d}_mod_fuv_red_density.fits'),
    'mod_fuv_red.reproject': (ANALYSIS_DIR, '_mod_fuv_red', 'reproject', 'hdu0_b{0:02d}_mod_fuv_red_density.fits'),
    'mod_fuv_red.reproject.add': (ANALYSIS_DIR, '_mod_fuv_red', 'add', 'mod_fuv_red_density.fits'),
    'mod_fuv_red.area': (ANALYSIS_DIR, '_mod_fuv_red', 'reproject', 'hdu0_b{0:02d}_mod_fuv_red_density_area.fits'),
    'mod_fuv_red.area.add': (ANALYSIS_DIR, '_mod_fuv_red', 'add', 'mod_fuv_red_density_area.fits'),
    'mod_fuv_red.hdr': (ANALYSIS_DIR, '_mod_fuv_red', 'template.hdr'),

    'galex_fuv': (GALEX_DIR, 'PS_M31_MOS{0:02d}-fd-int.fits'),
    'galex_fuv.add': (ANALYSIS_DIR, 'galex_fuv.fits'),
    'galex_fuv.density': (ANALYSIS_DIR, '_galex_fuv', 'input', 'MOS{0:02d}_galex_fuv_density.fits'),
    'galex_fuv.reproject': (ANALYSIS_DIR, '_galex_fuv', 'reproject', 'hdu0_MOS{0:02d}_galex_fuv_density.fits'),
    'galex_fuv.reproject.add': (ANALYSIS_DIR, '_galex_fuv', 'add', 'galex_fuv_density.fits'),
    'galex_fuv.area': (ANALYSIS_DIR, '_galex_fuv', 'reproject', 'hdu0_MOS{0:02d}_galex_fuv_density_area.fits'),
    'galex_fuv.area.add': (ANALYSIS_DIR, '_galex_fuv', 'add', 'galex_fuv_density_area.fits'),
    'galex_fuv.hdr': 'mod_fuv_red.hdr',

    'mod_nuv_int': (ANALYSIS_DIR, 'b{0:02d}', 'b{0:02d}_mod_nuv_int.fits'),
    'mod_nuv_int.add': (ANALYSIS_DIR, 'mod_nuv_int.fits'),
    'mod_nuv_int.density': (ANALYSIS_DIR, '_mod_nuv_int', 'input', 'b{0:02d}_mod_nuv_int_density.fits'),
    'mod_nuv_int.reproject': (ANALYSIS_DIR, '_mod_nuv_int', 'reproject', 'hdu0_b{0:02d}_mod_nuv_int_density.fits'),
    'mod_nuv_int.reproject.add': (ANALYSIS_DIR, '_mod_nuv_int', 'add', 'mod_nuv_int_density.fits'),
    'mod_nuv_int.area': (ANALYSIS_DIR, '_mod_nuv_int', 'reproject', 'hdu0_b{0:02d}_mod_nuv_int_density_area.fits'),
    'mod_nuv_int.area.add': (ANALYSIS_DIR, '_mod_nuv_int', 'add', 'mod_nuv_int_density_area.fits'),
    'mod_nuv_int.hdr': 'mod_fuv_red.hdr',

    'mod_nuv_red': (ANALYSIS_DIR, 'b{0:02d}', 'b{0:02d}_mod_nuv_red.fits'),
    'mod_nuv_red.add': (ANALYSIS_DIR, 'mod_nuv_red.fits'),
    'mod_nuv_red.density': (ANALYSIS_DIR, '_mod_nuv_red', 'input', 'b{0:02d}_mod_nuv_red_density.fits'),
    'mod_nuv_red.reproject': (ANALYSIS_DIR, '_mod_nuv_red', 'reproject', 'hdu0_b{0:02d}_mod_nuv_red_density.fits'),
    'mod_nuv_red.reproject.add': (ANALYSIS_DIR, '_mod_nuv_red', 'add', 'mod_nuv_red_density.fits'),
    'mod_nuv_red.area': (ANALYSIS_DIR, '_mod_nuv_red', 'reproject', 'hdu0_b{0:02d}_mod_nuv_red_density_area.fits'),
    'mod_nuv_red.area.add': (ANALYSIS_DIR, '_mod_nuv_red', 'add', 'mod_nuv_red_density_area.fits'),
    'mod_nuv_red.hdr': 'mod_fuv_red.hdr',

    'galex_nuv': (GALEX_DIR, 'PS_M31_MOS{0:02d}-nd-int.fits'),
    'galex_nuv.add': (ANALYSIS_DIR, 'galex_nuv.fits'),
    'galex_nuv.density': (ANALYSIS_DIR, '_galex_nuv', 'input', 'MOS{0:02d}_galex_nuv_density.fits'),
    'galex_nuv.reproject': (ANALYSIS_DIR, '_galex_nuv', 'reproject', 'hdu0_MOS{0:02d}_galex_nuv_density.fits'),
    'galex_nuv.reproject.add': (ANALYSIS_DIR, '_galex_nuv', 'add', 'galex_nuv_density.fits'),
    'galex_nuv.area': (ANALYSIS_DIR, '_galex_nuv', 'reproject', 'hdu0_MOS{0:02d}_galex_nuv_density_area.fits'),
    'galex_nuv.area.add': (ANALYSIS_DIR, '_galex_nuv', 'add', 'galex_nuv_density_area.fits'),
    'galex_nuv.hdr': 'mod_fuv_red.hdr'
    }


def _path(kind, extpar_dict, **kwargs):
    """File paths for the `sfhmaps_phat` project.

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
# .. note:: The bottleneck limiting import speed is in _load_extpar_dict,
#    probably in creating the `Table` instances.
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


EXTPAR_DICT = _load_extpar_dict()
MISSING_CELLS = _find_missing_cells(EXTPAR_DICT)


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
