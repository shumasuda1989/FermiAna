#!/usr/bin/env python

import sys
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.visualization import simple_norm
import regions
from scipy.optimize import curve_fit

plt.rcParams['savefig.directory']=''


def GetSkyMap(files):
    def get(filename):
        ihdu=0
        if filename.endswith(']'):
            tmp=filename.split('[')
            try:
                ihdu=int(tmp[-1][:-1])
            except ValueError:
                ihdu=tmp[-1][:-1]
            filename=tmp[0]

        hdu = fits.open(filename)
        data=hdu[ihdu].data.squeeze()
        head=hdu[ihdu].header
        try:
            head['CRPIX2']
        except KeyError:
            head=hdu[0].header
        w = WCS(head).sub(2)
        return data,w

    if not isinstance(files,list):
        files=files.split(',')

    print files
    d,ww=get(files[0])
    # print d.shape
    for fstack in files[1:]:
        dstack,_=get(fstack)
        # print dstack.shape
        d=np.vstack([d,dstack])

    # print d.shape
    return d,ww


def poisson_lnl(nc, mu):
    lnl = np.zeros(nc.shape)
    msk = nc > 0

    lnl[msk] = nc[msk] * np.log(mu[msk]) - mu[msk]
    lnl[~msk] = -mu[~msk]
    return lnl


def MakeResidmap(fcm,fmm,kernel=0.1,regfs=[],savefig=False,pointmodel=None):

    datacm,wc = GetSkyMap(fcm)
    datamm,w  = GetSkyMap(fmm)
    dataem = np.ones(datacm.shape)

    if datacm.shape != datamm.shape:
        print 'Error: Invalid inputs. Shapes of Count map and Model map must be same.'
        exit(1)
    if wc.wcs != w.wcs:
        print w.wcs
        print 'Error: Invalid inputs. Wcs of Count map and Model map must be same.'
        exit(1)
    xmax=datacm.shape[-1]; ymax=datacm.shape[-2]

    datacm=datacm.astype(datamm.dtype)
    if pointmodel is None:
        psfker = Gaussian2DKernel(kernel/w.wcs.cdelt[1])
        norm=np.sum(psfker.array**2)

        cm = convolve(datacm, psfker, normalize_kernel=False)
        mm = convolve(datamm, psfker, normalize_kernel=False)
        em = convolve(dataem, psfker, normalize_kernel=False)

        cmnorm=cm/norm
        mmnorm=mm/norm

    else:
        psfker,_ = GetSkyMap(pointmodel)
        if datacm.shape[0] != psfker.shape[0]:
            print 'Error: Invalid inputs. # of Ebin of Count map and PS model map must be same.'
            exit(1)

        psfker_cnt = psfker / np.sum(psfker, axis=(1,2),keepdims=True)
        psfker /= np.sum(psfker)
        psfker /= np.sum(psfker**2)

        cpix=xmax/2
        thrs=0.001
        plt.figure(figsize=(11.5,9))

        cm = np.zeros(datacm.shape)
        mm = np.zeros(datamm.shape)
        em = np.zeros(dataem.shape)
        cmnorm = np.zeros(datacm.shape)
        mmnorm = np.zeros(datamm.shape)

        for i in range(datacm.shape[0]):
            ax=plt.subplot(9,9,i+1)
            ker_c = psfker_cnt[i]
            ker = psfker[i]

            mpix=ker[cpix,:]>ker[cpix,cpix]*thrs
            npix = int(max(3, np.round(np.sum(mpix) / 2.)))
            spix = slice(cpix - npix, cpix + npix + 1)

            ker_c=ker_c[spix,spix]
            ker=ker[spix,spix]
            # print ker.shape

            cm[i] = convolve(datacm[i], ker_c, normalize_kernel=False)
            mm[i] = convolve(datamm[i], ker_c, normalize_kernel=False)
            em[i] = convolve(dataem[i], ker_c, normalize_kernel=False)
            cmnorm[i] = convolve(datacm[i], ker, normalize_kernel=False)
            mmnorm[i] = convolve(datamm[i], ker, normalize_kernel=False)

            ax.imshow(cm[i], origin='lower',norm=simple_norm(cm[i],'sqrt'))


        plt.tight_layout()
        cm=np.sum(cm,axis=0)
        mm=np.sum(mm,axis=0)
        em=np.sum(em,axis=0)
        cmnorm=np.sum(cmnorm,axis=0)
        mmnorm=np.sum(mmnorm,axis=0)

        if False:
            cm=cmnorm
            mm=mmnorm
            

    em /= np.max(em)

    excess = cm-mm
    excnorm = cmnorm-mmnorm

    ts = 2.0 * (poisson_lnl(cmnorm, cmnorm) - poisson_lnl(cmnorm, mmnorm))
    ts[ts<0] = 0
    sigma = np.sqrt(ts)
    sigma[excnorm < 0] *= -1
    # sigma = excnorm/np.sqrt(mmnorm)

    figlist=[]
    axlist=[]

    for skymap,scale,cmap in zip([cm, mm, excess],['sqrt','sqrt','linear'],['magma','magma','RdBu_r']):
        fig=plt.figure(figsize=(11.5,9))
        figlist.append(fig)
        ax = plt.subplot(projection=w)
        axlist.append(ax)
        kwargs=dict(origin='lower', interpolation='nearest')
        arr=skymap/em
        aa=np.min(arr); bb=np.max(arr)
        if aa < 0 :#and not 2./5 < -aa/(bb-aa) < 3./5:
            # tmp=(-aa+bb)/2
            tmp=max(-aa,bb)
            minmax=[-tmp,tmp]
        else:
            minmax=arr
        im=ax.imshow(arr, origin='lower', interpolation='spline16',norm=simple_norm(minmax,scale),cmap=cmap)

        ax.grid(color='white', ls='solid')
        if w.wcs.ctype[0].find('GLON') >-1:
            xlabel='Galactic Longitude'
            ylabel='Galactic Latitude'
        elif w.wcs.ctype[0].find('RA') >-1:
            xlabel=r'Right Ascension$_\mathrm{J%d}$' % int(w.wcs.equinox)
            ylabel=r'Declination$_\mathrm{J%d}$' % int(w.wcs.equinox)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        cbar=fig.colorbar(im)
        cbar.set_label('Counts')


    fig=plt.figure(figsize=(11.5,9))
    figlist.append(fig)
    ax = plt.subplot(projection=w)
    axlist.append(ax)
    im=ax.imshow(sigma, origin='lower', interpolation='spline16',
                 cmap='RdBu_r',vmax=5,vmin=-5)
    ax.contour(sigma, levels=[-3,3], colors='k', alpha=0.5)
    ax.grid(color='white', ls='solid')
    if w.wcs.ctype[0].find('GLON') >-1:
        xlabel='Galactic Longitude'
        ylabel='Galactic Latitude'
    elif w.wcs.ctype[0].find('RA') >-1:
        xlabel=r'Right Ascension$_\mathrm{J%d}$' % int(w.wcs.equinox)
        ylabel=r'Declination$_\mathrm{J%d}$' % int(w.wcs.equinox)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    cbar=fig.colorbar(im)
    cbar.set_label('Significance [$\sigma$]')

    for regf in regfs:
        regs = regions.read_ds9(regf)
        for reg in regs:
            meta=reg.meta; vis =reg.visual
            # print meta; print vis

            if hasattr(reg,'to_pixel'):
                reg=reg.to_pixel(w)

            x=reg.center.x; y=reg.center.y

            if -0.5 <= x < xmax-0.5  and  -0.5 <= y < ymax-0.5: pass
            else:   continue

            fontfam=vis.get('font','normal')
            if fontfam in ('roman','normal','helvetica'):
                fontfam='sans-serif'
            fontfam=fontfam.replace('times','serif')
            vis['color']=vis.get('color','lime').replace('green','lime')
            for i_ax,ax in enumerate(axlist):

                if isinstance(reg,regions.PointPixelRegion):
                    ax.plot(x,y,vis.get('symbol','o'), markerfacecolor='none',
                            markersize=float(vis.get('symsize','11'))/2.,
                            color=vis.get('color','lime'))
                else:
                    reg.plot(ax=ax,zorder=2)
                    if i_ax>0: pass
                    elif hasattr(reg,'radius'):  y+=reg.radius
                    elif hasattr(reg,'height'):
                        # print reg.width;print reg.height;print reg.angle
                        alpha=reg.angle.to('rad').value
                        theta=np.linspace(0,np.pi,100)
                        yy=max(np.abs( reg.width/2* np.sin(alpha)*np.cos(theta)
                                       +reg.height/2*np.cos(alpha)*np.sin(theta)))
                        y+=yy

                ax.text(x,y+ymax/100.,meta.get('label',''), 
                        color=vis.get('color','lime'),
                        family=fontfam,
                        size=float(vis.get('fontsize','10')),
                        style=vis.get('fontweight','roman').replace('roman','normal'),
                        weight=vis.get('fontstyle','normal'),
                        horizontalalignment='center',
                        verticalalignment='bottom')
    # plt.tight_layout()

    plt.rcParams['font.size']=12

    d=sigma.flatten()
    nbin=100
    xmax=max(6.4, np.max(d), -np.min(d))+0.1; xmin=-xmax

    figlist.append(plt.figure(figsize=(8,6)))
    ax=plt.subplot(111)

    n,bins,_=ax.hist(d,nbin,(xmin,xmax),density=True)

    mean  = np.mean(d)
    sigma = np.std(d)
    print mean,sigma

    y=lambda x,mu,sigma: np.exp(-((x-mu)/sigma)**2/2)/np.sqrt(2*np.pi)/sigma
    x=np.linspace(xmin,xmax,200)
    ax.plot(x,y(x,0,1),'-',color='k')

    x0=(bins[:-1]+bins[1:])/2.

    wh=np.where((n>0))
    popt,pcov=curve_fit(y,x0[wh],n[wh],[mean,sigma],np.sqrt(n[wh]))

    print popt
    plt.plot(x,y(x,popt[0],popt[1]),'--',color='r')


    ax.set_xlabel('Significance ($\sigma$)')
    ax.set_ylabel('Probability')

    # ax.text(0.05, 0.95, 'Gaussian fit:\n$\mu=%.2f$\n$\sigma=%.2f$'%(para[1],para[2]), transform=ax.transAxes,verticalalignment='top')
    ax.text(0.05, 0.95, 'Gaussian fit:\n$\mu=%.2f$\n$\sigma=%.2f$'%(popt[0],popt[1]), transform=ax.transAxes,verticalalignment='top',color='r')
    ax.text(0.77, 0.95, 'statistics:\nmean: $%.2f$\nrms:    %.2f'%(mean,sigma), transform=ax.transAxes,verticalalignment='top')

    figlist.append(plt.figure(figsize=(8,6)))
    ax=plt.subplot(111)
    ax.hist(d,nbin,(xmin,xmax),density=True,log=True)
    ax.plot(x,y(x,0,1),'-',color='k')
    plt.plot(x,y(x,popt[0],popt[1]),'--',color='r')
    ax.set_xlabel('Significance ($\sigma$)')
    ax.set_ylabel('Probability')

    if savefig is not None:
        fignamelist=[savefig+'_count',savefig+'_model',savefig+'_resid',
                     savefig+'_resid_sig',savefig+'_resid_sighist',
                     savefig+'_resid_sighistlog']
        if len(figlist) != len(fignamelist):
            print 'Error: len of figlist and fignamelist must be same.'
        for fig,figname in zip(figlist,fignamelist):
            fig.savefig(figname+'.png')
            fig.savefig(figname+'.pdf')

    plt.show()


if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('cmfits',help='count map fits file')
    parser.add_argument('mmfits',help='model map fits file')
    parser.add_argument('-r','--region', help='DS9 region file', action='append',default=[])
    parser.add_argument('-k','--kernel', type=float, help='map is smoothed with 2D Gaussian with this sigma value in unit of deg.',default=0.1)
    parser.add_argument('-s','--savefig', help='save images')#,action='store_true')
    parser.add_argument('-m','--modelpntsrc', help='point source model map')

    args = parser.parse_args()

    # # plt.rc('text', usetex=True)
    plt.rc('font', size=14)#, family='serif')


    MakeResidmap(args.cmfits,args.mmfits,args.kernel,args.region,args.savefig,args.modelpntsrc)
