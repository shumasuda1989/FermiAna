import matplotlib
from matplotlib import pyplot
import numpy as np
from numpy import array, nan

from SpectralModel import SpectralModelXML, SpectralModel, SpectralModelSEDdict

_ver=matplotlib.__version__.split('.')
_ver=float('.'.join(_ver[:2]))
_lims={'uplims':True} if _ver>=1.4 else {'lolims':True}

def setzorder(Art,zorder):
    if hasattr(Art, '__iter__'):
        for it in Art:
            setzorder(it,zorder)
    else:
        if hasattr(Art,'set_zorder'):
            Art.set_zorder(zorder)

def Spectrum_XML(ax,E,xml,src,dat=None,col='b',SED=True,kwspec=None,kwerr=None):

    Ep=np.logspace( np.log10(E[0]), np.log10(E[-1]), num=100 )

    if dat is None:
        f,_=SpectralModelXML(xml,src,lamb=True); ferr=lambda ee,ff: 0
    else:
        f,ferr,_=SpectralModelXML(xml,src,lamb=True,butt=True,resultsdat=dat)
    F=f(Ep)
    Ferr=ferr(Ep,F)

    if SED:
        F*=Ep**2
        Ferr*=Ep**2

    if dat is None:
        kwargs=dict(color=col)
        if kwspec is not None: kwargs.update(kwspec)                      
        return ax.plot(Ep,F,**kwargs)
    else:
        kwargs=dict(facecolor=col,alpha=0.2)
        if kwerr is not None: kwargs.update(kwerr)
        a=ax.fill_between(Ep,F-Ferr,F+Ferr,**kwargs)
        kwargs=dict(color=col)
        if kwspec is not None: kwargs.update(kwspec)                      
        b=ax.plot(Ep,F,**kwargs)
        return a,b

def Spectrum(ax,E,dat,src=None,col='b',xfac=1.,yfac=1.,SED=True,kwspec=None):
    Ep=np.logspace( np.log10(E[0]), np.log10(E[-1]), num=100 )

    if src is None:
        f,_=SpectralModelSEDdict(dat,lamb=True)
    else:
        f,_=SpectralModel(dat,src,lamb=True)
    F=f(Ep)
    if SED:
        F*=Ep**2
    kwargs=dict(color=col)
    if kwspec is not None: kwargs.update(kwspec)                      
    return ax.plot(Ep*xfac,F*yfac,**kwargs)


def Butterfly(ax,E,dat,src,col='b',xfac=1.,yfac=1.,SED=True,kwspec=None,kwerr=None):
    Ep=np.logspace( np.log10(E[0]), np.log10(E[-1]), num=100 )

    f,ferr,_=SpectralModel(dat,src,lamb=True,butt=True)
    F=f(Ep)
    Ferr=ferr(Ep,F)
    if SED:
        F*=Ep**2
        Ferr*=Ep**2

    kwargs=dict(facecolor=col,alpha=0.2)
    if kwerr is not None: kwargs.update(kwerr)
    a=ax.fill_between(Ep*xfac,(F-Ferr)*yfac,(F+Ferr)*yfac,**kwargs)

    kwargs=dict(color=col)
    if kwspec is not None: kwargs.update(kwspec)                      
    b=ax.plot(Ep*xfac,F*yfac,**kwargs)
    return a,b

def SEDGraph(ax,fdatapoint,col='b',xfac=1.,yfac=1.,SED=True,kwdata=None,kwul=None):
    pf=np.loadtxt(fdatapoint).T
    ul =np.where(pf[3]==-1)
    nul=np.where(pf[3]!=-1)

    if SED:
        pf[1]*=pf[0]**2
        pf[3]*=pf[0]**2

    kwargs=dict(fmt='o',linestyle='None',capsize=4,color=col)
    if kwdata is not None: kwargs.update(kwdata)
    a=ax.errorbar(pf[0][nul]*xfac,pf[1][nul]*yfac,yerr=pf[3][nul]*yfac,**kwargs)

    kwargs=dict(fmt='',linestyle='None',elinewidth=0.75,capsize=4,color=col)
    kwargs.update(_lims)
    if kwul is not None:   kwargs.update(kwul)
    b=ax.errorbar(pf[0][ul]*xfac,pf[1][ul]*yfac,yerr=pf[1][ul]*0.3*yfac,**kwargs)

    return a,b

def SEDGraphFromSEDdict(ax,dat,col='b',xfac=1.,yfac=1.,tsmin=4.,SED=True,kwdata=None,kwul=None):
    with open(dat) as f:
        dic=eval(f.read())
    E   = np.array(dic['Energy']['Value'])
    F   = np.array(dic['dNdE']['Value'])
    ErrF= np.array(dic['dNdE']['Average_Error'])
    UL  = np.array(dic['dNdE']['Upper_Limit'])
    TS  = np.array(dic['Test_Statistic'])

    ul =np.where(TS<=tsmin)
    nul=np.where(TS>tsmin)

    if SED:
        F   *=E**2
        ErrF*=E**2
        UL  *=E**2

    kwargs=dict(fmt='o',linestyle='None',capsize=4,color=col)
    if kwdata is not None: kwargs.update(kwdata)
    a=ax.errorbar(E[nul]*xfac,F[nul]*yfac,yerr=ErrF[nul]*yfac,**kwargs)

    kwargs=dict(fmt='',linestyle='None',elinewidth=0.75,capsize=4,color=col)
    kwargs.update(_lims)
    if kwul is not None:   kwargs.update(kwul)
    b=ax.errorbar(E[ul]*xfac,UL[ul]*yfac,yerr=UL[ul]*0.3*yfac,**kwargs)

    return a,b

def DrawSED(ax,E,resdat=None,src=None,seddat=None,fdatapoint=None,col='b',xfac=1.,yfac=1.,SED=True,kwspec=None,kwerr=None,kwdata=None,kwul=None):
    if   seddat is not None  and  resdat is src is fdatapoint is None:
        a=Spectrum(ax,E,seddat,col=col,xfac=xfac,yfac=yfac,SED=SED,kwspec=kwspec)
        b=SEDGraphFromSEDdict(ax,seddat,col=col,xfac=xfac,yfac=yfac,SED=SED,kwdata=kwdata,kwul=kwul)
        return a,b
    elif not None in (resdat,src,seddat)  and  fdatapoint is None:
        a=Butterfly(ax,E,resdat,src,col=col,xfac=xfac,yfac=yfac,SED=SED,kwspec=kwspec,kwerr=kwerr)
        b=SEDGraphFromSEDdict(ax,seddat,col=col,xfac=xfac,yfac=yfac,SED=SED,kwdata=kwdata,kwul=kwul)
        return a,b
    elif not None in (resdat,src,fdatapoint)  and  seddat is None:
        a=Butterfly(ax,E,resdat,src,col=col,xfac=xfac,yfac=yfac,SED=SED,kwspec=kwspec,kwerr=kwerr)
        b=SEDGraph(ax,fdatapoint,col=col,xfac=xfac,yfac=yfac,SED=SED,kwdata=kwdata,kwul=kwul)
        return a,b
    else:
        print 'DrawSED error: arguments are invalid.'

if __name__ == '__main__':

    dir='GammaCygni'
    # xml=dir+'/GammaCygni_output_model.xml'
    dat=dir+'/results.dat'
    src='gamma Cygni'
    fdatapoint=dir+'/DataPoints/SpectrumData_'+src.replace(' ','')+'.txt'
    E=[5e3,5e5]

    ax=pyplot.subplot(111)
    ax.loglog()

    # DrawSED_XML(src,ax,E,xml,dat,fdatapoint)
    DrawSED(ax,E,dat,src,fdatapoint)

    ax.set_xlabel('E [MeV]')
    ax.set_ylabel('E2 dN/dE [MeV cm-2 s-1]')
    pyplot.show()
