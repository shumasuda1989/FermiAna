#!/usr/bin/env python

import sys
try:
    import pyfits
    ow=dict(clobber=True)
except ImportError:
    import astropy.io.fits as pyfits
    ow=dict(overwrite=True)

with pyfits.open(sys.argv[1]) as f1, pyfits.open(sys.argv[2]) as f2:

    for hdu in f2[3:]:
        f1.append(hdu)

    if len(sys.argv) > 3:
        f1.writeto(sys.argv[3],**ow)
    else:
        f1.writeto(sys.argv[1],**ow)

