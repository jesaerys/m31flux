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
import os
import match_wrapper as match

from . import config, util


def make_bestzcb(norun=False):
    """Create zcombine files for the best SFH solutions only, no
    uncertainties.

    Parameters
    ----------
    norun : bool, optional
        If True, zcombine won't actually be executed. Default is False.

    Returns
    -------
    None

    """
    sfhfile_list = config.path('sfh', fillsubfield=False)
    zcbfile_list = config.path('bestzcb', fillsubfield=False)
    for sfhfile, zcbfile in zip(sfhfile_list, zcbfile_list):
        dirname = os.path.dirname(zcbfile)
        util.safe_mkdir(dirname)
        match.zcombine(sfhfile, zcbfile, bestonly=True, norun=norun)


def main():
    make_bestzcb()


if __name__ == "__main__":
    main()
