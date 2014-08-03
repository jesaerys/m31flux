"""

=====================
`m31flux.process_sfh`
=====================

Process calcsfh output files using zcombine.


Functions
---------

============== =========================================================
`make_bestzcb` Create zcombine files for the best SFH solutions only, no
               uncertainties.
============== =========================================================

"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import os
import match_wrapper as match

from . import config, util


def make_bestzcb():
    """Create zcombine files for the best SFH solutions only, no
    uncertainties.

    """
    print(
        '\n'
        'm31flux.process_sfh.make_bestzcb\n'
        '--------------------------------'
        )
    prefix = 'Processing .sfh files with zcombine (all bricks, all cells): '
    c = util.make_counter(prefix=prefix, end='  done')
    sfhfile_list = config.path('sfh', fillsubfield=False)
    zcbfile_list = config.path('bestzcb', fillsubfield=False)
    for sfhfile, zcbfile in c(zip(sfhfile_list, zcbfile_list)):
        dirname = os.path.dirname(zcbfile)
        try:
            os.makedirs(dirname)
        except OSError:
            pass
        match.zcombine(sfhfile, zcbfile, bestonly=True)
    return
