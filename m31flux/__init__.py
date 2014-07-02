"""

=========
`m31flux`
=========

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


Modules
-------

============= ============================================
`config`      Configuration for `m31flux`.
`util`        General utilities.
`process_sfh` Process calcsfh output files using zcombine.
`make_maps`   Create modeled flux maps from SFH data.
============= ============================================

"""
from . import config
from . import util
