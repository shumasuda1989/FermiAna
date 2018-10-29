#!/usr/bin/env python

import sys
import pyfits

f1=pyfits.open(sys.argv[1])
f2=pyfits.open(sys.argv[2])

for hdu in f2[3:]:
    f1.append(hdu)

if len(sys.argv) > 3:
    f1.writeto(sys.argv[3],clobber=True)
else:
    f1.writeto(sys.argv[1],clobber=True)

