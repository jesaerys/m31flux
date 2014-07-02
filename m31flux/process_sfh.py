"""
Process calcsfh output files using zcombine.

"""
import os
import match_wrapper as match

from sfhmaps_phat import config, util


def make_bestzcb(norun=False):
    """Create zcombine files for the best SFH solutions only, no
    uncertainties.

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
