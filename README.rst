m31flux
=======

Synthetic broadband flux maps of M31 from resolved optical photometry.

Jacob E. Simones, Julianne J. Dalcanton, Andrew E. Dolphin, Benjamin D.
Johnson, Alexia R. Lewis, Evan D. Skillman, Daniel R. Weisz, Benjamin F.
Williams

Star formation histories derived from resolved optical photometry from PHAT
were used to create maps of broadband flux in the ultraviolet and mid
infrared. The maps cover roughly one third of the disk of M31 and have a
resolution of ~430 pc. The maps were compared with observations from GALEX
and Spitzer.

This repository contains the following:

- `m31flux`: the data processing and analysis python package for the M31
  flux maps project. These packages are required:

  - `astrogrid <http://github.com/jesaerys/astrogrid>`_
  - `astropy <http://www.astropy.org>`_
  - `match-wrapper <http://github.com/jesaerys/match-wrapper>`_
  - `numpy <http://www.numpy.org>`_

- scripts: python scripts, including main, which serves as the frontend
  to `m31flux` and is the complete research and analysis record.

- maps: modeled and observed flux maps in FITS format produced by
  scripts/main.

- figs: figures based on the flux maps and their analyses.

- docs: sphinx documentation and paper drafts. Coming soon!
