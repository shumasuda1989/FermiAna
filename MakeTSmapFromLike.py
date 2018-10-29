#!/usr/bin/env python

import os
from time import sleep
import subprocess
import numpy as np
import pyfits
from pywcs import WCS
from ReadResult import ReadResult

env=os.environ

# ra=float(env['xref'])
# dec=float(env['yref'])

# nxpix=int(env['nxpix'])
# nypix=int(env['nypix'])
# delta=float(env['binsz'])

ra=float(env['Ra'])
dec=float(env['Dec'])

nxpix=int(env['npix3'])
nypix=int(env['npix3'])
delta=float(env['binszts'])
proj=env['proj']

Npix=1

if nxpix%Npix != 0 or nypix%Npix != 0:
    print 'invalid pixel number!!'
    exit(1)

option=env.get('option',',24h').split(',')
if len(option) != 2:
    exit(1)

w = WCS(naxis=2)

w.wcs.crpix = [(nxpix+1)/2.,(nypix+1)/2.]
w.wcs.cdelt = np.array([-delta,delta])
w.wcs.crval = [ra,dec]

w.wcs.ctype = ['RA---%s'%proj, 'DEC--%s'%proj]
w.wcs.equinox = 2000.0

head=w.to_header()



x=np.array([ Npix*(i+0.5)-0.5 for i in range(nxpix/Npix)])
y=np.array([ Npix*(i+0.5)-0.5 for i in range(nypix/Npix)])
XX,YY=np.meshgrid(x,y)

NN=nxpix*nypix/Npix**2
RA,DEC = w.wcs_pix2sky(XX.reshape(NN),YY.reshape(NN),0)
RA=RA.reshape((nxpix/Npix,nypix/Npix))
DEC=DEC.reshape((nxpix/Npix,nypix/Npix))

lockdir='tslockLike'
filedir='tsfileLike'

if not os.path.isdir(lockdir):
    os.mkdir(lockdir)
if not os.path.isdir(filedir):
    os.mkdir(filedir)

with pyfits.open(env['tsmap']) as fi:
    wh=np.where(fi[0].data==0)

if env.get('SKIP','') == '':

    with open(env['srcmdlin']) as f:
        lines=f.read()
        lines=lines.replace('</source_library>\n','')
        lines=lines.replace('</source_library>','')

    # for i in range(nxpix/Npix):
    #     for j in range(nypix/Npix):

    for i,j in zip(wh[0],wh[1]):
        print i,j,'RA:',RA[i,j], 'DEC:',DEC[i,j]

        testsrc='  <source name="putative source" type="PointSource">\n    <spectrum apply_edisp="false" type="PowerLaw">\n      <parameter free="1" max="10000" min="0.0001" name="Prefactor" scale="1e-17" value="1.00" />\n      <parameter free="1" max="10" min="0" name="Index" scale="-1" value="2.0" />\n      <parameter free="0" max="500000" min="30" name="Scale" scale="1" value="10000" />\n    </spectrum>\n    <spatialModel type="SkyDirFunction">\n      <parameter free="0" max="360" min="-360" name="RA" scale="1" value="%f" />\n      <parameter free="0" max="90" min="-90" name="DEC" scale="1" value="%f" />\n    </spatialModel>\n  </source>\n</source_library>\n' % (RA[i,j],DEC[i,j])

        lock='%s/lock_%d_%d'%(lockdir,i,j)
        log='%s/log_%d_%d.txt'%(filedir,i,j)
        xmlin='%s/%s'%(filedir,env['srcmdlin'].replace('.xml','_%d_%d.xml'%(i,j)))
        xmlout='%s/%s'%(filedir,env['srcmdlout'].replace('.xml','_%d_%d.xml'%(i,j)))
        results='%s/%s'%(filedir,env['results'].replace('.dat','_%d_%d.dat'%(i,j)))
        with open(xmlin, mode='w') as f:
            f.write(lines+testsrc)

        command='subcmd.sh %s "source %s/fermi-init.sh; touch %s; srcmap=%s lvtime=%s bexpcube=%s srcmdlin=%s srcmdlout=%s optimizer=%s METHOD=%s specfile=None results=%s USE_BL_EDISP=%s refit=%s plot=%s slist=\'putative source\' ../Likelihood.py; rm %s" %s %s' %(option[0],env['FERMI_DIR'], lock, env['srcmap'], env['lvtime'], env['bexpcube'], xmlin, xmlout, env['optimizer'], env['METHOD'], results, env['USE_BL_EDISP'], env['refit'], env['plot'], lock, log, option[1])

            # if j%2==0:
            #     cmd1='touch %s; srcmap=%s lvtime=%s bexpcube=%s srcmdlin=%s srcmdlout=%s optimizer=%s METHOD=%s specfile=None results=%s USE_BL_EDISP=%s refit=%s plot=%s slist=\'putative source\' ../Likelihood.py; rm %s' % (lock, env['srcmap'], env['lvtime'], env['bexpcube'], xmlin, xmlout, env['optimizer'], env['METHOD'], results, env['USE_BL_EDISP'], env['refit'], env['plot'], lock)

            # else:
            #     command='subcmd.sh %s "source %s/fermi-init.sh; %s; touch %s; srcmap=%s lvtime=%s bexpcube=%s srcmdlin=%s srcmdlout=%s optimizer=%s METHOD=%s specfile=None results=%s USE_BL_EDISP=%s refit=%s plot=%s slist=\'putative source\' ../Likelihood.py; rm %s" %s %s' %(option[0],env['FERMI_DIR'], cmd1, lock, env['srcmap'], env['lvtime'], env['bexpcube'], xmlin, xmlout, env['optimizer'], env['METHOD'], results, env['USE_BL_EDISP'], env['refit'], env['plot'], lock, log, option[1])
                # print command
                # exit(0)
        subprocess.call(command, shell=True)


print 'waiting for TS map calculation in each bin...'

while len(os.listdir(lockdir)) > 0:
    sleep(10)

print 'Done!!'

# data=np.zeros((nxpix,nypix))
data=np.full((nxpix,nypix),-1)

for i in range(nxpix/Npix):
    for j in range(nypix/Npix):
        results='%s/%s'%(filedir,env['results'].replace('.dat','_%d_%d.dat'%(i,j)))
        xmlout='%s/%s'%(filedir,env['srcmdlout'].replace('.xml','_%d_%d.xml'%(i,j)))
        if os.path.isfile(results) and os.path.isfile(xmlout):
            print '%d_%d'%(i,j)
            res,_=ReadResult(results)
            data[Npix*i:Npix*(i+1),Npix*j:Npix*(j+1)]=res['putative source']['TS value']
            if data[Npix*i:Npix*(i+1),Npix*j:Npix*(j+1)][0]==0:
                print '=0'

with pyfits.open(env['tsmap']) as fi:
    head=fi[0].header

# head['NAXIS1']=nxpix
# head['NAXIS2']=nypix
# head['CRPIX1']=(nxpix+1)/2.
# head['CRPIX2']=(nypix+1)/2.
# head['CRVAL1']=ra
# head['CRVAL2']=dec

newhdu=pyfits.PrimaryHDU(data=data,header=head)
newhdu.writeto(env['tsmap'].replace('.fits','Like.fits'),clobber=True,checksum=True)




