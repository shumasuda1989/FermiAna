#!/usr/bin/env python

import os
from time import sleep
import subprocess
import numpy as np
import pyfits
from pywcs import WCS

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

Npix=int(env['Npix'])

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

lockdir='tslock'
filedir='tsfile'

if not os.path.isdir(lockdir):
    os.mkdir(lockdir)
if not os.path.isdir(filedir):
    os.mkdir(filedir)


if env.get('SKIP') is None or env.get('SKIP') == '':

    for i in range(nxpix/Npix):
        for j in range(nypix/Npix):

            print i,j,'RA:',RA[i,j], 'DEC:',DEC[i,j]

            lock='%s/lock_%d_%d'%(lockdir,i,j)
            log='%s/log_%d_%d.txt'%(filedir,i,j)
            tsmap='%s/%s'%(filedir,env['tsmap'].replace('.fits','_%d_%d.fits'%(i,j)))
            command=('subcmd.sh %s "source %s/fermi-init.sh; touch %s; gttsmap statistic=%s evfile=%s scfile=%s bexpmap=%s expcube=%s srcmdl=%s cmap=%s outfile=%s irfs=CALDB optimizer=%s nxpix=%s nypix=%s binsz=%s xref=%s yref=%s coordsys=%s proj=%s; rm %s" %s %s' %
                     (option[0],env['FERMI_DIR'], lock, env['METHOD'], env['evfile'], env['scfile'], env['bexpcube'], env['lvtime'], env['srcmdlfix'], env['ccube'], tsmap, env['optimizer'], Npix, Npix, env['binszts'], RA[i,j], DEC[i,j], env['coordsys'], env['proj'], lock, log, option[1]))

            subprocess.call(command, shell=True)

print 'waiting for TS map calculation in each bin...'

while len(os.listdir(lockdir)) > 0:
    sleep(10)

print 'Done!!'

data=np.zeros((nxpix,nypix))

for i in range(nxpix/Npix):
    for j in range(nypix/Npix):
        tsmap='%s/%s'%(filedir,env['tsmap'].replace('.fits','_%d_%d.fits'%(i,j)))
        print '%d_%d'%(i,j)
        with pyfits.open(tsmap) as hdul:
            data[Npix*i:Npix*(i+1),Npix*j:Npix*(j+1)]=hdul[0].data


with pyfits.open('%s/%s'%(filedir,env['tsmap'].replace('.fits','_0_0.fits'))) as fi:
    head=fi[0].header

head['NAXIS1']=nxpix
head['NAXIS2']=nypix
head['CRPIX1']=(nxpix+1)/2.
head['CRPIX2']=(nypix+1)/2.
head['CRVAL1']=ra
head['CRVAL2']=dec

newhdu=pyfits.PrimaryHDU(data=data,header=head)
newhdu.writeto(env['tsmap'],clobber=True,checksum=True)
