#!/usr/bin/env python

import os
import time
import traceback
import gc
# import pickle
from BinnedAnalysis import BinnedObs, BinnedAnalysis, pyLike
from SummedLikelihood import SummedLikelihood
from LikelihoodState import LikelihoodState
from UpperLimits import UpperLimits
from IntegralUpperLimit import calc_int
from collections import OrderedDict
from SED import SED 

ullist = []


def OptimizeModel(like1,out=None,slist=[],TSmax=4.):

    print 'All sources with TS < %g will be deleted.'%TSmax
    like1.optimize(optObject=like1.optObject)
    for source in like1.sourceNames():
        TS=like1.Ts(source)
        if TS < TSmax and source not in slist:
            print "Deleting",source,'(TS = %g)'%TS
            like1.deleteSource(source)
    like1.syncSrcParams()
    print
    if out is not None:
        like1.logLike.writeXml(out)


def GetObsList(lobssum):
    obs=[]
    lobs=lobssum.split('|')
    for l0 in lobs:
        l=[ s.strip() for s in l0.split(',') ]
        print l
        obs.append(BinnedObs(srcMaps=l[0],expCube=l[1],
                             binnedExpMap=l[2],irfs='CALDB'))
        print "IRFs:", obs[-1].irfs

    return obs

def GetLikeObj(obs,modelin,optimizer,tol=1e-3):

    if isinstance(obs,list) or isinstance(obs,tuple):
        like1 = SummedLikelihood(optimizer=optimizer)
        for il,obs1 in enumerate(obs):
            print 'Creating a Likelihood object for obs%d...'%il
            like2 = BinnedAnalysis(obs1,modelin,optimizer=optimizer)
            like1.addComponent(like2)
    elif isinstance(obs,BinnedObs):
        like1 = BinnedAnalysis(obs,modelin,optimizer=optimizer)
    else:
        print 'error: invalid ObsObject'
        exit(1)

    return like1

def GetOptimizer(like1,optimizer):

    if   optimizer=='NEWMINUIT':
        like1obj = pyLike.NewMinuit(like1.logLike)
    elif optimizer=='MINUIT':
        like1obj = pyLike.Minuit(like1.logLike)
    else:
        optFactory = pyLike.OptimizerFactory_instance()
        like1obj = optFactory.create(optimizer, like1.logLike)

    return like1obj


def Fit(like1,out=None,covar=True):

    logL=like1.fit(verbosity=3,covar=covar,optObject=like1.optObject)

    if out is not None:
        # like1.logLike.writeXml(out)
        like1.writeXml(out)

    like1obj=like1.optObject
    optimizer=like1.optimizer

    if optimizer=='MINUIT':
        print 'Quality=', like1obj.getQuality()
    print 'ReturnCode=', like1obj.getRetCode()

    if optimizer=='MINUIT' and like1obj.getQuality() == 3:
        return True

    if optimizer=='NEWMINUIT' and like1obj.getRetCode() == 0:
        return True

    return False


def WriteResult(dic,results):
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
        elif k=='-LogLike':
            s+="'%s': %f,\n"%(k,v)
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

def PrintResult(dic):
    for k, v in dic.items():
        if k in ['Energies','Free Parameters','Covariance','-LogLike'] \
           or not isinstance(v,dict):
            continue
        else:
            print k+':'
            for kk, vv in v.items():
                if kk in ['Spectrum','EnFlux','Diff Flux'] \
                   or kk.startswith('scale ')  or 'UL' in kk:
                    continue
                elif kk=='Flux':
                    if isinstance(vv,tuple) and len(vv)==2:
                        print '%s: %6g +/- %6g photons/cm^2/s'% (kk,vv[0],vv[1])
                    else:
                        print '%s: %6g photons/cm^2/s'% (kk,vv)
                else:
                    if isinstance(vv,tuple) and len(vv)==2:
                        print '%s: %6g +/- %6g'% (kk,vv[0],vv[1])
                    else:
                        print '%s: %6g'% (kk,vv)
            print ''

def GetSED(like1,sname, index=-2., be=None, min_ts=4., ul_alg='bayesian', fbasename=None):
    ### ul_choices = ['frequentist', 'bayesian']
    try:
        sed = SED(like1,sname, min_ts=min_ts, bin_edges=be, ul_algorithm=ul_alg,do_minos=False, always_upper_limit=True, powerlaw_index=index)
        if  fbasename is None:
            fbasename='sed_%s_%c'%(sname.replace(' ',''),ul_alg[0])
        sed.save(fbasename+'.dat')
        if os.environ.get('DISPLAY') is not None:
            sed.plot(fbasename+'.png')
        return sed.todict()

    except Exception:
        print(traceback.format_exc())
        pass



def Likelihood(like1,modelout,optimizer,statistic,specfile,results,plot,slist,optmodel=(False,),SkipUL=False,Bayes=False,binedge=None,dellist=None):

    if statistic != 'BINNED':
        print '%s method is not implemented in this script'%statistic
        return None

    if dellist is not None:
        for sdel in dellist:
            like1.deleteSource(sdel)
            print '%s is deleted'%sdel
        like1.syncSrcParams()
    if optmodel[0]:
        OptimizeModel(like1,slist=slist,TSmax=optmodel[1])

    like1.optObject = GetOptimizer(like1,optimizer)
    if like1.optimizer != optimizer:
        like1.optimizer = optimizer
    if isinstance(like1,SummedLikelihood): # not needed, results aren't changed
        for l in like1.components:
            l.optimizer=optimizer

    nloop=0
    fitxml=modelout+'_fit%d'%(nloop+1)
    while nloop < 3 and not Fit(like1,fitxml):
        print; print
        nloop+=1
        fitxml=modelout+'_fit%d'%(nloop+1)
        like1.optObject=GetOptimizer(like1,optimizer)
        gc.collect()

    if nloop == 3:
        print 'could not converge'
        fitxml=modelout+'_fit3'
    # pkl=results+'.pkl'
    # pkl=pkl.replace('.dat.pkl','.pkl')
    # with open(pkl, mode='wb') as f:
    #     pickle.dump(like, f)
    #     pickle.dump(likeobj, f)

    if isinstance(like1,SummedLikelihood):
        is_esame=True
        for l in like1.components[1:]:
            is_esame=(is_esame and 
                      np.allclose(like1.components[0].energies, l.energies) )
        if is_esame:
            E=like1.components[0].energies
        else:
            E=(min([l.energies[0]  for l in like1.components]),
               max([l.energies[-1] for l in like1.components]))
    else:
        E=like1.energies
    emin=E[0]
    emax=E[-1]

    print '\nComputing TS values for each extended source\n'
    print 'Photon fluxes are computed for the energy range '+repr(emin)+' to '+repr(emax)+' MeV\n'

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
        enflux=like1.energyFlux(sname,emin=emin,emax=emax)
        if flag:
            if sname.find('gll')<0 and sname.find('iso')<0:
                for param in func.paramNames:
                    src['scale '+param]=func.params[param].getScale()
                src['Spectrum']=func.genericName()
                src['TS value']=like1.Ts(sname)
            eflx=like1.fluxError(sname,emin=emin,emax=emax)
            eenflx=like1.energyFluxError(sname,emin=emin,emax=emax)
            src['Flux']=(flux,eflx)
            src['EnFlux']=(enflux,eenflx)
        else:
            src['Flux']=flux
            src['EnFlux']=enflux

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
    dic['-LogLike']=like1()

    WriteResult(dic,results)
    PrintResult(dic)

    '''
      Calculate Upper Limit. You can choose bayesian or frequentist algorithm.

    '''
    if SkipUL or len(slist)==0:
        print 'skip UL calculation'
    elif Bayes:
        ul = {}
        ullist.append(ul)
        
        for sname in slist:
            print sname
            flux_ul,ulresults = calc_int(like1, sname, emin=emin, emax=emax)
            try_ul_calc=1
            par=like1.normPar(sname)
            while flux_ul==-1 and try_ul_calc<=10:
                par.setBounds(par.getBounds()[0],par.getBounds()[1]*10)
                like1.syncSrcParams(sname)
                flux_ul,ulresults = calc_int(like1, sname, emin=emin, emax=emax)
                try_ul_calc+=1
            if flux_ul==-1:
                print 'could not compute upper limit'
            else:
                print '%lg ph/cm^2/s for emin=%.1f, emax=%.1f (Bayesian UL)'%(flux_ul,ulresults['flux_emin'],ulresults['flux_emax'])
                ul[sname]=ulresults
                dic[sname]['dNdE UL']=ulresults['ul_value']*par.getScale()
                ## dN/dE UL at a reference enrergy (ph/cm^2/s)
                dic[sname]['Flux UL']=flux_ul
                saved_state = LikelihoodState(like1)
                par.setValue(ulresults['ul_value'])
                dic[sname]['EnFlux UL']=like1.energyFlux(sname,emin,emax)
                saved_state.restore()
                dic[sname]['UL algo']='bayesian'
                WriteResult(dic,results)

    else:
        ul = UpperLimits(like1)
        ullist.append(ul)

        for sname in slist:
            print sname
            try:
                flux_ul,pref_ul=ul[sname].compute(emin=emin,emax=emax)
                print ul[sname].results[-1]
                dic[sname]['dNdE UL']=pref_ul*like1.normPar(sname).getScale()
                dic[sname]['Flux UL']=flux_ul
                saved_state = LikelihoodState(like1)
                like1.normPar(sname).setValue(pref_ul)
                dic[sname]['EnFlux UL']=like1.energyFlux(sname,emin,emax)
                saved_state.restore()
                dic[sname]['UL algo']='frequentist'
                dic[sname]['UL dlogL']=ul[sname].results[-1].delta
                WriteResult(dic,results)
            except RuntimeError:
                print(traceback.format_exc())
                print('could not compute upper limit')


    # with open(pkl, mode='ab') as f:
    #     pickle.dump(ul, f)

    try:
        if specfile!='None' and specfile!='' and specfile is not None:
            like1.writeCountsSpectra(specfile,len(E)-1)
    except RuntimeError:
        print(traceback.format_exc())
        pass
    except NotImplementedError as e:
        print 'NotImplementedError:',e
        pass

    cnt=0.
    for name in like1.model.srcNames:
        cnt+=like1.NpredValue(name)
    print '\nTotal number of observed counts:', int(like1.total_nobs())
    print 'Total number of model events:', cnt
    print '\n-log(Likelihood):', like1()

    print '\nWriting fitted model to', modelout
    # like1.logLike.writeXml(modelout)
    os.rename(fitxml,modelout)


    '''
      Calculate Fermi SED data points.

    '''
    if binedge is not None:
        ul_alg='bayesian' if Bayes else 'frequentist'
        print '\n%s Upper Limit will be calculated' % ul_alg
        print 'calculating SED...'
        for sname in slist:
            func=like1.model[sname].funcs['Spectrum']
            if func.genericName()=='PowerLaw':
                idx=func['Index']*func.params['Index'].getScale()
            else:
                idx=-2.
            GetSED(like1,sname,index=idx,be=binedge,ul_alg=ul_alg)
        print 'Done!'


    '''
      Plot count graph.

    '''
    if plot:
        try:
            print ''
            like1.setPlotter('mpl')
            like1.plot()
            print ''
        except:
            print 'could not plot'
            pass


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


def mybool(Input,default=False):
    env=os.environ
    if env.get(Input) is None:
        return default
    elif env[Input].lower() in ['t','true','y','yes','1']:
        return True
    elif env[Input].lower() in ['f','false','n','no','0']:
        return False
    else:
        print 'Invalid value in %s (= "%s")'% (Input,env[Input])
        exit(1)

def GetEnv():

    env=os.environ

    srcMaps=env['srcmap']
    expCube=env['lvtime']
    binnedExpMap=env['bexpcube']
    IRFs=env.get('irfs','CALDB')
    modelin  =env['srcmdlin']
    modelout =env['srcmdlout']
    optimizer=env['optimizer']
    statistic=env['METHOD']
    specfile=env.get('specfile','counts_spectra.fits')
    results=env.get('results','results.dat')
    # USE_BL_EDISP=env['USE_BL_EDISP']
    refit=mybool('refit')
    plot=mybool('plot')

    snamelist=env.get('slist','')
    if snamelist=='':
        slist=[]
    else:
        slist=snamelist.split(',')

    return srcMaps,expCube,binnedExpMap,IRFs,modelin,modelout,optimizer,statistic,specfile,results,refit,plot,slist






if __name__ == '__main__':

    srcMaps,expCube,binnedExpMap,IRFs,modelin,modelout,optimizer,statistic,specfile,results,refit,plot,slist = GetEnv()

    env=os.environ

    TSmax=float(env.get('TS_Max','4'))
    optmdl=(True,TSmax) if mybool('OPT_MODEL') else (False,)
    skipul=mybool('SKIP_UL')
    bayes=mybool('BAYES')
    be=env.get('E_bin')
    if be is not None:
        import numpy as np
        be=np.fromstring(be.replace(',',' '),sep=' ')
        print 'E_bin:',be
    tol=float(env.get('Tolerance','1e-3'))
    dellist=env.get('dellist','')
    if dellist=='':
        dellist=[]
    else:
        dellist=dellist.split(',')

    start=time.time()

    if statistic != 'BINNED':
        print '%s method is not implemented in this script'%statistic
        exit(1)

    lobssum=env.get('SUM_LIST') 
    # srcmap1, ltcube1, expcube1 | srcmap2, ltcube2, expcube2 | ...
    if lobssum is not None:
        obs = GetObsList(lobssum)
    else:
        obs = BinnedObs(srcMaps=srcMaps,expCube=expCube,
                        binnedExpMap=binnedExpMap,irfs=IRFs)

    like1=GetLikeObj(obs,modelin,optimizer)

    if refit:
        optim_refit=env.get('optimizer_prerefit',optimizer)
        tol_refit=env.get('Tolerance_prerefit',tol)
        like1.tol = tol_refit; print 'tolerance = ',like1.tol

        Likelihood(like1,modelout+'_refit',optim_refit,statistic,None,results,plot,slist,optmdl,True,bayes,None,dellist)

        print '\nRefit\n'
        modelin=modelout+'_refit'
        optmdl=(False,)
        dellist=None
        gc.collect()

    like1.tol = tol; print 'tolerance = ',like1.tol
    Likelihood(like1,modelout,optimizer,statistic,specfile,results,plot,slist,optmdl,skipul,bayes,be,dellist)

    etime=time.time()-start
    print 'Elapsed CPU time:',etime,'sec'

    if plot:
        try:
            import code
            code.InteractiveConsole(globals()).interact()
        except:
            pass



'''
import sys
sys.path.append('..')
from Likelihood import *
srcMaps,expCube,binnedExpMap,IRFs,modelin,modelout,optimizer,statistic,specfile,results,refit,plot,slist = GetEnv()
obs = BinnedObs(srcMaps=srcMaps,expCube=expCube,binnedExpMap=binnedExpMap,irfs=IRFs)

#optimizer='MINUIT'
like1=GetLikeObj(obs,modelin,optimizer)
like1.tol=
like1.optObject = GetOptimizer(like1,optimizer)
Fit(like1)

E=like1.energies
emin=E[0]
emax=E[-1]

sname='

ul = UpperLimits(like1)
flux_ul,pref_ul=ul[sname].compute(emin=emin,emax=emax)
ul[sname].results[-1]

flux_ul,ulresults = calc_int(like1, sname, emin=emin, emax=emax)



'''
