import matplotlib
from matplotlib import pyplot
import numpy as np
from numpy import array

from SpectralModel import SpectralModelXML, SpectralModel, SpectralModelSEDdict

ver=matplotlib.__version__.split('.')
ver=float('.'.join(ver[:2]))
lims={'uplims':True} if ver>=1.4 else {'lolims':True}

def DrawSED_XML(src,ax,E,xml,dat,fdatapoint,col='b',SED=True):

    Ep=np.logspace( np.log10(E[0]), np.log10(E[-1]), num=100 )

    f,ferr,_=SpectralModelXML(xml,src,lamb=True,butt=True,resultsdat=dat)
    F=f(Ep)
    Ferr=ferr(Ep,F)

    if SED:
        F*=Ep**2
        Ferr*=Ep**2
    ax.fill_between(Ep,F-Ferr,F+Ferr,facecolor=col,alpha=0.2)
    ax.plot(Ep,F,color=col)

    pf=np.loadtxt(fdatapoint).T
    ul =np.where(pf[3]==-1)
    nul=np.where(pf[3]!=-1)

    if SED:
        pf[1]*=pf[0]*pf[0]
        pf[3]*=pf[0]*pf[0]

    ax.errorbar(pf[0][nul],pf[1][nul],yerr=pf[3][nul], fmt='o',linestyle='None',capsize=4,color=col)
    ax.errorbar(pf[0][ul],pf[1][ul],yerr=pf[1][ul]*0.3, fmt='',linestyle='None',elinewidth=0.75,capsize=4,color=col,**lims)

def Spectrum(ax,E,dat,src=None,col='b',xfac=1.,yfac=1.,SED=True):
    Ep=np.logspace( np.log10(E[0]), np.log10(E[-1]), num=100 )

    if src is None:
        f,_=SpectralModelSEDdict(dat,lamb=True)
    else:
        f,_=SpectralModel(dat,src,lamb=True)
    F=f(Ep)
    if SED:
        F*=Ep**2
    ax.plot(Ep*xfac,F*yfac,color=col)


def Butterfly(ax,E,dat,src,col='b',xfac=1.,yfac=1.,SED=True):
    Ep=np.logspace( np.log10(E[0]), np.log10(E[-1]), num=100 )

    f,ferr,_=SpectralModel(dat,src,lamb=True,butt=True)
    F=f(Ep)
    Ferr=ferr(Ep,F)
    if SED:
        F*=Ep**2
        Ferr*=Ep**2

    ax.fill_between(Ep*xfac,(F-Ferr)*yfac,(F+Ferr)*yfac,
                    facecolor=col,alpha=0.2)
    ax.plot(Ep*xfac,F*yfac,color=col)

def SEDGraph(ax,fdatapoint,col='b',xfac=1.,yfac=1.,SED=True):
    pf=np.loadtxt(fdatapoint).T
    ul =np.where(pf[3]==-1)
    nul=np.where(pf[3]!=-1)

    if SED:
        pf[1]*=pf[0]**2
        pf[3]*=pf[0]**2

    ax.errorbar(pf[0][nul]*xfac,pf[1][nul]*yfac,yerr=pf[3][nul]*yfac, 
                fmt='o',linestyle='None',capsize=4,color=col)
    ax.errorbar(pf[0][ul]*xfac,pf[1][ul]*yfac,yerr=pf[1][ul]*0.3*yfac, 
                fmt='',linestyle='None',elinewidth=0.75,
                capsize=4,color=col,**lims)


def SEDGraphFromSEDdict(ax,dat,col='b',xfac=1.,yfac=1.,tsmin=4.,SED=True):
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

    ax.errorbar(E[nul]*xfac,F[nul]*yfac,yerr=ErrF[nul]*yfac, 
                fmt='o',linestyle='None',capsize=4,color=col)
    # ax.plot(E[nul]*xfac,F[nul]*yfac,'o',linestyle='None',color=col)
    ax.errorbar(E[ul]*xfac,UL[ul]*yfac,yerr=UL[ul]*0.3*yfac, 
                fmt='',linestyle='None',elinewidth=0.75,
                capsize=4,color=col,**lims)

def DrawSED(ax,E,resdat=None,src=None,seddat=None,fdatapoint=None,col='b',xfac=1.,yfac=1.,SED=True):
    if   seddat is not None  and  resdat is src is fdatapoint is None:
        Spectrum(ax,E,seddat,col=col,xfac=xfac,yfac=yfac,SED=SED)
        SEDGraphFromSEDdict(ax,seddat,col=col,xfac=xfac,yfac=yfac,SED=SED)
    elif not None in (resdat,src,seddat)  and  fdatapoint is None:
        Butterfly(ax,E,resdat,src,col=col,xfac=xfac,yfac=yfac,SED=SED)
        SEDGraphFromSEDdict(ax,seddat,col=col,xfac=xfac,yfac=yfac,SED=SED)
    elif not None in (resdat,src,fdatapoint)  and  seddat is None:
        Butterfly(ax,E,resdat,src,col=col,xfac=xfac,yfac=yfac,SED=SED)
        SEDGraph(ax,fdatapoint,col=col,xfac=xfac,yfac=yfac,SED=SED)
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
