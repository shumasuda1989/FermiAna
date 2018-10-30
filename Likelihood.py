#!/usr/bin/env python

import os
import time
# import pickle
from BinnedAnalysis import BinnedObs, BinnedAnalysis, pyLike
from UpperLimits import UpperLimits
from IntegralUpperLimit import calc_int
from collections import OrderedDict
from SED import SED 

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


def Formatter(dic,results):
    s='{'
    for k, v in dic.items():
        if k=='Energies':
            tmps=repr(v).replace('array','').replace('[','').replace(']','')
            s+="'%s': %s,\n"%(k,tmps)
        elif k=='Free Parameters':
            s+="'%s': %s,\n"%(k,repr(v))
        elif k=='Covariance':
            s+="'%s': (\n" % k
            for row in v:
                s+='('
                for val in row:
                    s+='%6g, '%val
                s+='),\n'
            s+='),\n'
        else:
            s+="'%s': {" % k
            for kk, vv in v.items():
                if isinstance(vv,tuple) and len(vv)==2:
                    s+="'%s': '%6g +/- %6g',\n" % (kk,vv[0],vv[1])
                elif isinstance(vv,list):
                    tmps=''
                    for val in vv:
                        tmps+='%6g, '%val
                    s+="'%s': (%s),\n" % (kk,tmps)
                elif isinstance(vv,str):
                    s+="'%s': '%s',\n" % (kk,vv)
                else:
                    s+="'%s': '%6g',\n" % (kk,vv)
            s+='},\n'
    s+='}\n'
    with open(results,mode='w') as f:
        f.write(s)


def GetSED(like1,sname, min_ts=9, ul_alg='bayesian',be=None, fbasename=None):
    ### ul_choices = ['frequentist', 'bayesian']
    try:
        sed = SED(like1,sname, min_ts=min_ts, bin_edges=be ,ul_algorithm=ul_alg,do_minos=False)
        if  fbasename is None:
            fbasename='sed_%s_%c'%(sname.replace(' ',''),ul_alg[0])
        sed.save(fbasename+'.dat')
        if os.environ.get('DISPLAY','')!='':
            sed.plot(fbasename+'.png')
        return sed.todict()

    except Exception as e:
        print e
        pass



def Likelihood(srcMaps,expCube,binnedExpMap,modelin,modelout,optimizer,statistic,specfile,results,slist,SkipUL=False,Bayes=False,binedge=None):

    obs = BinnedObs(srcMaps=srcMaps,expCube=expCube,binnedExpMap=binnedExpMap,irfs='CALDB')

    nloop=0
    fitxml_pre=modelin
    fitxml=modelout.replace('.xml','_fit%d.xml'%(nloop+1))
    while nloop < 3 and not Fit(like,likeobj,obs,fitxml_pre,optimizer,fitxml):
        nloop+=1
        fitxml_pre=fitxml
        fitxml=modelout.replace('.xml','_fit%d.xml'%(nloop+1))

    if nloop == 3:
        print 'could not converge'
        fitxml=modelout.replace('.xml','_fit3.xml')

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

    print '\nComputing TS values for each extended source\n'
    print 'Photon fluxes are computed for the energy range '+repr(emin)+' to '+repr(emax)+' MeV\n\n'

    dic=OrderedDict()
    free=[]
    for sname in like1.model.srcNames:
        flag=False
        src=OrderedDict()
        func=like1.model[sname].funcs['Spectrum']
        for param in func.paramNames:
            flag=(flag or func.params[param].isFree())
            if func.params[param].isFree():
                free.append('%s:%s'%(sname,param))
                src[param]=(func[param],func.params[param].error())
            else:
                src[param]=func[param]

        flux=like1.flux(sname,emin=emin,emax=emax)
        if flag:
            if sname.find('gll')<0 and sname.find('iso')<0:
                src['TS value']=like1.Ts(sname)
            eflx=like1.fluxError(sname,emin=emin,emax=emax)
            src['Flux']=(flux,eflx)
        else:
            src['Flux']=flux

        if like1.model[sname].getType() == 'Diffuse':
            tmplist=[]
            nbin=len(E)-1
            for ie in range(nbin):
                tmplist.append(like1.flux(sname,emin=E[ie],emax=E[ie+1]))
            src['Diff Flux']=tmplist
        dic[sname]=src
    dic['Energies']=E
    if len(free)>1:
        dic['Free Parameters']=free
        dic['Covariance']=like1.covariance

    Formatter(dic,results)


    '''
      Calculate Upper Limit. You can choose bayesian or frequentist algorithm.

    '''
    if Bayes and not SkipUL:
        ul = {}
        ullist.append(ul)
        
        for sname in slist:
            print sname
            try:
                flux_ul,ulresults = calc_int(like1, sname, emin=emin, emax=emax)
                print '%lg ph/cm^2/s for emin=%.1f, emax=%.1f (Bayesian UL)'%(flux_ul,ulresults['flux_emin'],ulresults['flux_emax'])
                ul[sname]=ulresults
                dic[sname]['Flux UL']=flux_ul
                dic[sname]['UL algo']='bayesian'
            except RuntimeError as e:
                import traceback
                print(traceback.format_exc())
                print e
                print('could not compute upper limit')
    elif not SkipUL:
        ul = UpperLimits(like1)
        ullist.append(ul)

        for sname in slist:
            print sname
            try:
                flux_ul,_=ul[sname].compute(emin=emin,emax=emax)
                print ul[sname].results[-1]
                dic[sname]['Flux UL']=flux_ul
                dic[sname]['UL algo']='frequentist'
                dic[sname]['UL dlogL']=ul[sname].results[-1].delta
            except RuntimeError as e:
                import traceback
                print(traceback.format_exc())
                print e
                print('could not compute upper limit')

    Formatter(dic,results)

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


    '''
      Calculate Fermi SED data points.

    '''
    if binedge is not None:
        print 'calculating SED...'
        ul_alg='bayesian' if Bayes else 'frequentist'
        for sname in slist:
            GetSED(like1,sname, min_ts=9, ul_alg=ul_alg,be=binedge)
        print 'Done!'


    return dic

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

    skipul=True if os.environ.get('SKIP_UL','')!='' else False
    bayes=True if os.environ.get('BAYES','')!='' else False
    be=os.environ.get('E_bin')
    if be is not None:
        import numpy as np
        be=np.fromstring(be.replace(',',' '),sep=' ')
        print 'E_bin:',be

    start=time.time()
    Likelihood(srcMaps,expCube,binnedExpMap,modelin,modelout,optimizer,statistic,specfile,results,slist,skipul,bayes,be)
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