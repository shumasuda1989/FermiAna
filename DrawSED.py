from matplotlib import pyplot
import numpy as np

from SpectralModel import SpectralModel

def DrawSED(E,xml,dat,fdatapoint):

    Ep=np.logspace( np.log10(E[0]), np.log10(E[-1]), num=100 )

    f,ferr,_=SpectralModel(xml,src,lamb=True,butt=True,resultsdat=dat)
    F=f(Ep)
    Ferr=ferr(Ep,F)

    pyplot.fill_between(Ep,Ep**2*(F-Ferr),Ep**2*(F+Ferr),facecolor='r',alpha=0.2)
    pyplot.plot(Ep,Ep**2*F,color='r')

    pf=np.loadtxt(fdatapoint).T
    ul =np.where(pf[3]==-1)
    nul=np.where(pf[3]!=-1)

    pf[1]*=pf[0]*pf[0]
    pf[3]*=pf[0]*pf[0]

    pyplot.errorbar(pf[0][nul],pf[1][nul],yerr=pf[3][nul], fmt='o',linestyle='None',capsize=4)
    pyplot.errorbar(pf[0][ul],pf[1][ul],yerr=pf[1][ul]*0.3, fmt='',linestyle='None',elinewidth=0.75,lolims=True,capsize=4,color='b') #uplims=True



if __name__ == '__main__':

    dir='GammaCygni'
    xml=dir+'/GammaCygni_output_model.xml'
    dat=dir+'/results.dat'
    src="gamma Cygni"
    fdatapoint=dir+'/DataPoints/SpectrumData_'+src.replace(' ','')+'.txt'
    E=[5e3,5e5]

    pyplot.subplot(111)
    pyplot.loglog()

    DrawSED(E,xml,dat,fdatapoint)

    pyplot.xlabel('E [MeV]')
    pyplot.ylabel('E2 dN/dE [MeV cm-2 s-1]')
    pyplot.show()
