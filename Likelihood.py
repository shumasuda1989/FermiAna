#!/usr/bin/env python

import os
import time
# import pickle
from BinnedAnalysis import BinnedObs, BinnedAnalysis, pyLike
from UpperLimits import UpperLimits
from IntegralUpperLimit import calc_int

like = []
likeobj = []
ullist = []

def Fit(likelist,likeobjlist,obs,modelin,optimizer,out=None,like1=None,obj=None,tol=1e-3):

    if like1 is None:
        like1 = BinnedAnalysis(obs,modelin,optimizer=optimizer)

    like1.tol = tol
    print 'tolerance = ',like1.tol
    if obj is not None:
        like1obj = obj
    elif optimizer=='NEWMINUIT':
        like1obj = pyLike.NewMinuit(like1.logLike)
    elif optimizer=='MINUIT':
        like1obj = pyLike.Minuit(like1.logLike)
    else:
        like1obj = pyLike.Optimizer(like1.logLike)

    logL=like1.fit(verbosity=3,covar=True,optObject=like1obj)

    if len(likelist)==0 or  not like1 in likelist:
        likelist.append(like1)
    if len(likeobjlist)==0 or  not like1obj in likeobjlist:
        likeobjlist.append(like1obj)

    if out is not None:
        like1.logLike.writeXml(out)

    print '\n'
    if optimizer=='MINUIT':
        print 'Quality=', like1obj.getQuality()
    print 'ReturnCode=', like1obj.getRetCode()

    if optimizer=='MINUIT' and like1obj.getQuality() == 3:
        return True

    if optimizer=='NEWMINUIT' and like1obj.getRetCode() == 0:
        return True

    return False



def Likelihood(srcMaps,expCube,binnedExpMap,modelin,modelout,optimizer,statistic,specfile,results,slist,Bayes=False):

    obs = BinnedObs(srcMaps=srcMaps,expCube=expCube,binnedExpMap=binnedExpMap,irfs='CALDB')

    nloop=0
    fitxml_pre=modelin
    fitxml=modelout.replace('output_model','output_fit%d'%(nloop+1))
    while nloop < 3 and not Fit(like,likeobj,obs,fitxml_pre,optimizer,fitxml):
        nloop+=1
        fitxml_pre=fitxml
        fitxml=modelout.replace('output_model','output_fit%d'%(nloop+1))

    if nloop == 3:
        print 'could not converge'
        fitxml=modelout.replace('output_model','output_fit3')

    # pkl=results+'.pkl'
    # pkl=pkl.replace('.dat.pkl','.pkl')
    # with open(pkl, mode='wb') as f:
    #     pickle.dump(like, f)
    #     pickle.dump(likeobj, f)

    like1=like[-1]
    like1obj=likeobj[-1]
    E=like1.energies
    emin=E[0]
    emax=E[-1]

    file = open(results,'w')

    print '\nComputing TS values for each extended source\n'
    print 'Photon fluxes are computed for the energy range '+repr(emin)+' to '+repr(emax)+' MeV\n\n'
    nsts=0
    for sname in like1.model.srcNames:
        flag=False
        for param in like1.model[sname].funcs['Spectrum'].paramNames:
            flag=flag or like1.model[sname].funcs['Spectrum'].params[param].isFree()
        file.write("'%s': {\n" % sname)
        flux=like1.flux(sname,emin=emin,emax=emax)

        if flag and sname.find('gll')<0 and sname.find('iso')<0:
            nsts+=1
            eflx=like1.fluxError(sname,emin=emin,emax=emax)
            file.write("'Flux': '"+repr(flux)+" +- "+repr(eflx)+"',\n")

            ts=like1.Ts(sname)
            file.write("'TS value': '"+repr(ts)+"',\n")
        else:
            file.write("'Flux': '"+repr(flux)+"',\n")

        if like1.model[sname].getType() == 'Diffuse':
            file.write("'Diff Flux':\n")
            nbin=len(E)-1
            for ie in range(nbin):
                flux=like1.energyFlux(sname,emin=E[ie],emax=E[ie+1])
                file.write(repr(flux))
                if ie != nbin-1:
                    file.write(' ')
            file.write('\n')

        file.write('},\n')

    if len(slist) >0:
        file.write('\nFlux UL:\n')
    if Bayes:
        for sname in slist:
            print sname
            file.write("'%s':\n"%sname)
            try:
                flux_ul,results = calc_int(like1, sname, emin=emin, emax=emax)
                s='%lg ph/cm^2/s for emin=%.1f, emax=%.1f (Bayesian UL)'%(flux_ul,results['flux_emin'],results['flux_emax'])
                file.write(s+'\n')
            except RuntimeError as e:
                import traceback
                print(traceback.format_exc())
                print e
                file.write('could not compute upper limit')
    else:
        ul = UpperLimits(like1)
        ullist.append(ul)

        for sname in slist:
            print sname
            file.write("'%s':\n"%sname)
            try:
                flux_ul,_=ul[sname].compute(emin=emin,emax=emax)
                s=repr(ul[sname].results[-1])
                stmp=s.split(); stmp.pop(0); stmp.insert(0,str(flux_ul))
                s=' '.join(stmp)
                file.write(s+'\n')
            except RuntimeError as e:
                import traceback
                print(traceback.format_exc())
                print e
                file.write('could not compute upper limit')
    file.close()
    # with open(pkl, mode='ab') as f:
    #     pickle.dump(ul, f)

    try:
        if specfile!='None':
            like1.writeCountsSpectra(outfile=specfile,nee=len(E)-1)
    except RuntimeError as e:
        print e
        pass


    cnt=0.
    for name in like1.model.srcNames:
        cnt+=like1.NpredValue(name)
    print '\nTotal number of observed counts:', int(like1.total_nobs())
    print 'Total number of model events:', cnt
    print '\n-log(Likelihood):', like1()
    print '\n'

    print 'Writing fitted model to', modelout
    # like1.logLike.writeXml(modelout)
    os.rename(fitxml,modelout)

    # like2 = BinnedAnalysis(obs,modelout1,optimizer=optimizer2)
    # like2.tol = tolerance2
    # like2obj = pyLike.NewMinuit(like2.logLike)
    # like2.fit(verbosity=0,covar=True,optObject=like2obj)
    # like2.logLike.writeXml(modelout2)
    # print like2obj.getRetCode()
    # like2.Ts(sname)
    # like2.model[sname]
    # like2.flux(sname,emin=emin,emax=emax)
    # like2.fluxError(sname,emin=emin,emax=emax)

def GetEnv():

    env=os.environ

    srcMaps=env['srcmap']
    expCube=env['lvtime']
    binnedExpMap=env['bexpcube']
    modelin  =env['srcmdlin']
    modelout =env['srcmdlout']
    optimizer=env['optimizer']
    statistic=env['METHOD']
    specfile=env.get('specfile','counts_spectra.fits')
    results=env.get('results','results.dat')
    USE_BL_EDISP=env['USE_BL_EDISP']
    refit=env.get('refit')
    plot=env.get('plot')

    snamelist=env.get('slist','')
    if snamelist=='':
        slist=[]
    else:
        slist=snamelist.split(',')

    return srcMaps,expCube,binnedExpMap,modelin,modelout,optimizer,statistic,specfile,results,USE_BL_EDISP,refit,plot,slist






if __name__ == '__main__':

    srcMaps,expCube,binnedExpMap,modelin,modelout,optimizer,statistic,specfile,results,USE_BL_EDISP,refit,plot,slist = GetEnv()
    bayes=True if os.environ.get('BAYES','')!='' else False

    start=time.time()
    Likelihood(srcMaps,expCube,binnedExpMap,modelin,modelout,optimizer,statistic,specfile,results,slist,bayes)
    etime=time.time()-start
    print 'Elapsed time:',etime,'sec'


'''
import sys
sys.path.append('..')
from Likelihood import *
srcMaps,expCube,binnedExpMap,modelin,modelout,optimizer,statistic,specfile,results,USE_BL_EDISP,refit,plot,slist = GetEnv()
obs = BinnedObs(srcMaps=srcMaps,expCube=expCube,binnedExpMap=binnedExpMap,irfs='CALDB')

#optimizer='MINUIT'
Fit(like,likeobj,obs,modelin,optimizer,out=,like1=like1,obj=like1obj)

like1=like[-1]
like1obj=likeobj[-1]
E=like1.energies
emin=E[0]
emax=E[-1]

sname='

ul = UpperLimits(like1)
ul[sname].compute(emin=emin,emax=emax)
ul[sname].results[-1]

flux_ul,results = calc_int(like1, sname, emin=emin, emax=emax)



'''
