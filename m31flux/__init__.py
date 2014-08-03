"""

=========
`m31flux`
=========

Data processing and analysis package for the M31 flux maps project.

See scripts/main.py in the main m31flux repository for more complete
documentation.

`m31flux` requires these packages:

- `astrogrid <http://github.com/jesaerys/astrogrid>`_
- `astropy <http://www.astropy.org>`_
- `match-wrapper <http://github.com/jesaerys/match-wrapper>`_
- `numpy <http://www.numpy.org>`_


Modules
-------

============= ============================================
`config`      Configuration for `m31flux`.
`fluxmap`     Create modeled flux maps from SFH data.
`process_sfh` Process calcsfh output files using zcombine.
`util`        General utilities.
============= ============================================

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

from . import config
from . import fluxmap
from . import process_sfh
from . import util
