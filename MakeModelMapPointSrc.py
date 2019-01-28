#!/usr/bin/env python

import os
import pyfits
from pywcs import WCS
from BinnedAnalysis import BinnedObs, BinnedAnalysis

env=os.environ

content_xml='''<?xml version="1.0" standalone="no"?>
<source_library title="source library">
  <source name="testsource" type="PointSource">
    <spectrum apply_edisp="false" type="PowerLaw">
      <parameter free="0" max="10000" min="0.0001" name="Prefactor" scale="1e-10" value="1.0" />
      <parameter free="0" max="10" min="0" name="Index" scale="-1" value="%.1f" />
      <parameter free="0" max="500000" min="30" name="Scale" scale="1" value="1000" />
    </spectrum>
    <spatialModel type="SkyDirFunction">
      <parameter free="0" max="360" min="-360" name="RA" scale="1" value="%f" />
      <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="%f" />
    </spatialModel>
  </source>
</source_library>
'''


def main(srcMaps=None,expCube=None,binnedExpMap=None,outname='modelcube.fits',
         IRFs='CALDB',index=2):
    if srcMaps is None:
        srcMaps=env['ccube']
        expCube=env['lvtime']
        binnedExpMap=env['bexpcube']
        IRFs=env.get('irfs','CALDB')

    obs = BinnedObs(srcMaps=srcMaps,expCube=expCube,
                    binnedExpMap=binnedExpMap,irfs=IRFs)

    with pyfits.open(srcMaps) as f:
        w=WCS(f[0].header)
        cpix=int(f[0].header['NAXIS1'])/2
    RA,DEC,e = w.wcs_pix2sky([[cpix,cpix,0]],0)[0]

    xmlfile=outname.replace('.fits','')+'.xml'

    with open(xmlfile,'w') as f:
        f.write(content_xml%(index,RA,DEC))

    like1=BinnedAnalysis(obs,xmlfile,optimizer='')
    os.remove(xmlfile)

    print 'creating Model Map...'
    like1.writeModelMap(outname)

    print 'Done!'

if __name__ == '__main__':
    import sys
    if len(sys.argv)>1:
        ### count cube, lvtime cube, exposure cube, output model
        main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
    else:
        out=env.get('mdlmap','modelcube.fits').replace('.fits','_ps.fits')
        main(outname=out)
