# usege: python ShiftPhase.py INPUTFILE PHASESHIFT

import sys
import pyfits
import numpy as np

hdulist = pyfits.open(sys.argv[1],mode='update')

# hdulist.info()
# prihdr = hdulist[0].header
# print repr(prihdr)

Nev = hdulist[1].header['NAXIS2']
print Nev
tbdata = hdulist[1].data
col=tbdata.field('PULSE_PHASE')
print col.shape[0]
print col

col[:]+=float(sys.argv[2])
col[:]%=1.0

hdulist.flush()

print col

hdulist.close()
