import sys
import os
import xml.etree.ElementTree as ET
import numpy as np
from numpy import array
import glob

sqrt=np.sqrt
log=np.log
exp=np.exp
fs=np.fromstring

def GetVal(node,scale=1.):
    return  float(node.get('value'))*float(node.get('scale'))/scale
def GetErr(node,scale=1.):
    return  float(node.get('error','0'))*float(node.get('scale'))/scale

def GetIndex(xml=None,root=None): #obsolete

    if xml is None and root is None or xml is not None and root is not None:
        exit(1)
    elif xml is not None and root is None:
        tree = ET.parse(xml)
        root = tree.getroot()
    free=[]
    ind=0
    for source in root:
        spec=source.find('spectrum')
        for param in spec:
            if param.get('free') == '1':
                # print ind,source.get('name'),'\t',param.get('name')
                free.append('%s:%s'%(source.get('name'),param.get('name')))
                ind+=1

    return free

def GetCov(fitdatatxt): #obsolete

    d=open(fitdatatxt)
    lines=[s.strip() for s in d.readlines()]
    d.close()

    indf,indl=0,0
    for i,line in enumerate(lines):
        if 'MnUserCovariance:' in line:
            indf=i
        elif 'MnUserCovariance Parameter correlations:' in line:
            indl=i

    cov=[]
    scov=lines[indf+2:indl-1]

    for s in scov:
        row=fs(s,sep=' ')
        cov.append(row)

    cov=np.array(cov)

    return cov


def PowerLawXML(spectrum,scale=1.):
    NodePrefactor = spectrum.find("parameter[@name='Prefactor']")
    NodeIndex     = spectrum.find("parameter[@name='Index']")
    NodeScale     = spectrum.find("parameter[@name='Scale']")
    if NodePrefactor is None or NodeScale is None:
        print 'ERROR!!!!!!!! Cannot find prefactor or scale nodes'
        exit(1)
    Prefactor = GetVal(NodePrefactor,scale)
    scale_n   = float(NodePrefactor.get('scale'))/scale
    EPre      = GetErr(NodePrefactor,scale)
    Index     = GetVal(NodeIndex)
    scale_g   = float(NodeIndex.get('scale'))
    EIndex    = GetErr(NodeIndex)
    Scale     = GetVal(NodeScale)

    return ("[0]*pow(x/[2],[1])",(Prefactor,Index,Scale,scale_n,scale_g),
            (EPre,EIndex))

def PowerLaw(src,scale=1.):
    if 'name' in src.keys() and 'method' in src.keys():
        scale_p   = 1/scale
        Prefactor = src['Prefactor']*scale_p
        scale_i   = 1
        Index     = src['Index']
        Scale     = src['Scale']
        return ("[0]*pow(x/[2],[1])",
                (Prefactor,Index,Scale),(0,0))
    else:
        scale_p   = float(src['scale Prefactor'])/scale
        Prefactor = fs(src['Prefactor'],sep='+/-')*scale_p
        scale_i   = float(src['scale Index'])
        Index     = fs(src['Index'],sep='+/-')*scale_i
        Scale     = float(src['Scale'])*float(src['scale Scale'])
        return ("[0]*pow(x/[2],[1])",
                (Prefactor[0],Index[0],Scale,scale_p,scale_i),
                (Prefactor[1],Index[1]))


def ButterflyPL(N0,N0err,E0,gerr,cov_ng):
    dn=1/N0
    dg=lambda E: log(E/E0)
    return lambda E,F: F*sqrt((dn*N0err)**2+(dg(E)*gerr)**2 +2*dn*dg(E)*cov_ng)


def PLSuperExpCutoffXML(spectrum,scale=1.):
    NodePrefactor = spectrum.find("parameter[@name='Prefactor']")
    NodeIndex1    = spectrum.find("parameter[@name='Index1']")
    NodeScale     = spectrum.find("parameter[@name='Scale']")
    NodeCutoff    = spectrum.find("parameter[@name='Cutoff']")
    NodeIndex2    = spectrum.find("parameter[@name='Index2']")
    if NodePrefactor is None or NodeScale is None:
        print 'ERROR!!!!!!!! Cannot find prefactor or scale nodes'
        exit(1)
    Prefactor = GetVal(NodePrefactor,scale)
    scale_n   = float(NodePrefactor.get('scale'))/scale
    EPre      = GetErr(NodePrefactor,scale)
    Index1    = GetVal(NodeIndex1)
    scale_g   = float(NodeIndex1.get('scale'))
    EInd1     = GetErr(NodeIndex1)
    Scale     = GetVal(NodeScale)
    Cutoff    = GetVal(NodeCutoff)
    scale_c   = float(NodeCutoff.get('scale'))
    ECoff     = GetErr(NodeCutoff)
    Index2    = GetVal(NodeIndex2)

    return ("[0]*pow(x/[2],[1])*exp(-pow(x/[3],[4]))",
            (Prefactor,Index1,Scale,Cutoff,Index2,scale_n,scale_g,scale_c), 
            (EPre,EInd1,ECoff))

def PLSuperExpCutoff(src,scale=1.):
    if 'name' in src.keys() and 'method' in src.keys():
        scale_p   = 1/scale
        Prefactor = src['Prefactor']*scale_p
        scale_i1  = 1
        Index1    = src['Index1']
        Scale     = src['Scale']
        scale_c   = 1
        Cutoff    = src['Cutoff']
        Index2    = src['Index2']
        return ("[0]*pow(x/[2],[1])*exp(-pow(x/[3],[4]))",
                (Prefactor,Index1,Scale,Cutoff,Index2), (0,0,0))
    else:
        scale_p   = float(src['scale Prefactor'])/scale
        Prefactor = fs(src['Prefactor'],sep='+/-')*scale_p
        scale_i1  = float(src['scale Index1'])
        Index1    = fs(src['Index1'],sep='+/-')*scale_i1
        Scale     = float(src['Scale'])*float(src['scale Scale'])
        scale_c   = float(src['scale Cutoff'])
        Cutoff    = fs(src['Cutoff'],sep='+/-')*scale_c
        Index2    = float(src['Index2'])*float(src['scale Index2'])
        return ("[0]*pow(x/[2],[1])*exp(-pow(x/[3],[4]))",
                (Prefactor[0],Index1[0],Scale,Cutoff[0],Index2,
                 scale_p,scale_i1,scale_c), 
                (Prefactor[1],Index1[1],Cutoff[1])              )

def ButterflyPLC(Norm,E0,Ec,Nerr,gerr,Ecerr,cov_ng,cov_gc,cov_cn):
    dn=1/Norm
    dg=lambda E: log(E/E0)
    dc=lambda E: E/Ec**2

    return lambda E,F: F*sqrt((dn*Nerr)**2 +(dg(E)*gerr)**2+(dc(E)*Ecerr)**2
                                 +2*dn*dg(E)*cov_ng
                                 +2*dg(E)*dc(E)*cov_gc
                                 +2*dc(E)*dn*cov_cn)

def LogParabolaXML(spectrum,scale=1.):
    NodeNorm  = spectrum.find("parameter[@name='norm']")
    NodeAlpha = spectrum.find("parameter[@name='alpha']")  
    NodeBeta  = spectrum.find("parameter[@name='beta']")
    NodeEb    = spectrum.find("parameter[@name='Eb']")
    if NodeNorm is None or NodeAlpha is None or NodeBeta is None or NodeEb is None:
        print 'ERROR!!!!!!!! Cannot find prefactor or scale nodes'
        exit(1)
    Norm   = GetVal(NodeNorm,scale)
    scale_n= float(NodeNorm.get('scale'))/scale
    ENorm  = GetErr(NodeNorm,scale)
    Alpha  = GetVal(NodeAlpha)
    EAlpha = GetErr(NodeAlpha)
    Beta   = GetVal(NodeBeta)
    EBeta  = GetErr(NodeBeta)
    Eb     = GetVal(NodeEb)

    return  ("[0]*pow(x/[3],-[1]-[2]*log(x/[3]))",
             (Norm,Alpha,Beta,Eb,scale_n), (ENorm,EAlpha,EBeta))

def LogParabola(src,scale=1.):
    if 'name' in src.keys() and 'method' in src.keys():
        scale_n= 1/scale
        Norm   = src['norm']*scale_n
        Alpha  = src['alpha']
        Beta   = src['beta']
        Eb     = src['Eb']
        return  ( "[0]*pow(x/[3],-[1]-[2]*log(x/[3]))",
                  (Norm,Alpha,Beta,Eb), (0,0,0)        )
    else:
        scale_n= float(src['scale norm'])/scale
        Norm   = fs(src['norm'],sep='+/-')*scale_n
        Alpha  = fs(src['alpha'],sep='+/-')*float(src['scale alpha'])
        Beta   = fs(src['beta'], sep='+/-')*float(src['scale beta'])
        Eb     = float(src['Eb'])*float(src['scale Eb'])
        return  ("[0]*pow(x/[3],-[1]-[2]*log(x/[3]))",
                 (Norm[0],Alpha[0],Beta[0],Eb,scale_n), 
                 (Norm[1],Alpha[1],Beta[1])            )

def ButterflyLP(Norm,E0,nerr,aerr,berr,cov_na,cov_ab,cov_bn):
    dn=1/Norm
    da=lambda E: -log(E/E0)
    db=lambda E: -(log(E/E0))**2

    return lambda E,F: F*sqrt((dn*nerr)**2 +(da(E)*aerr)**2+(db(E)*berr)**2
                                 +2*dn*da(E)*cov_na
                                 +2*da(E)*db(E)*cov_ab
                                 +2*db(E)*dn*cov_bn)
    # return (lambda E:  (dn*nerr)**2 +(da(E)*aerr)**2+(db(E)*berr)**2,
    #         lambda E: ((dn*nerr)**2 +(da(E)*aerr)**2+(db(E)*berr)**2
    #                     +2*dn*da(E)*cov_na
    #                     +2*da(E)*db(E)*cov_ab
    #                     +2*db(E)*dn*cov_bn))

def SpectralModelXML(xml,sname,xmin=100.,xmax=5.e5,scale=1,lamb=False,butt=False,resultsdat=None):
    if not hasattr(SpectralModelXML,'num'):
        SpectralModelXML.num=0
    SpectralModelXML.num+=1

    print xml

    tree = ET.parse(xml)
    root = tree.getroot()

    if butt:
        if not os.path.isfile(resultsdat):
            print 'cannot find %s' % resultsdat
            exit(1)
        # free=GetIndex(root=root)
        # cov=GetCov(fitdatatxt)
        with open(resultsdat) as f:
            dic=eval(f.read())
        free=dic['Free Parameters']
        cov=np.array(dic['Covariance'])
        

    spectrum = root.find("./source[@name='%s']/spectrum" % sname)
    if spectrum is not None:
        ftype=spectrum.get('type')
        if ftype == 'PowerLaw':
            print 'PowerLaw'
            formula,lpar,lerr=PowerLawXML(spectrum,scale)
            if butt:
                indn=free.index('%s:Prefactor'%sname)
                indg=free.index('%s:Index'%sname)
                cov_ng=cov[indn,indg]*lpar[3]*lpar[4]
                elamb=ButterflyPL(lpar[0],lerr[0],lpar[2],lerr[1],cov_ng)
            lpar=lpar[:-2]

        elif ftype == 'PLSuperExpCutoff':
            print 'PLSuperExpCutoff'
            formula,lpar,lerr=PLSuperExpCutoffXML(spectrum,scale)
            if butt:
                ind_n=free.index('%s:Prefactor'%sname)
                ind_g=free.index('%s:Index1'%sname)
                ind_c=free.index('%s:Cutoff'%sname)
                cov_ng=cov[ind_n,ind_g]*lpar[5]*lpar[6]
                cov_gc=cov[ind_g,ind_c]*lpar[6]*lpar[7]
                cov_cn=cov[ind_c,ind_n]*lpar[7]*lpar[5]
                elamb=ButterflyPLC(lpar[0],lpar[2],lpar[3],
                                   lerr[0],lerr[1],lerr[2],
                                   cov_ng,cov_gc,cov_cn)
            lpar=lpar[:-3]

        elif ftype == 'LogParabola':
            print 'LogParabola'
            formula,lpar,lerr=LogParabolaXML(spectrum,scale)
            if butt:
                ind_n=free.index('%s:norm'%sname)
                ind_a=free.index('%s:alpha'%sname)
                ind_b=free.index('%s:beta'%sname)
                cov_na=cov[ind_n,ind_a]*lpar[4]
                cov_ab=cov[ind_a,ind_b]
                cov_bn=cov[ind_b,ind_n]*lpar[4]
                elamb=ButterflyLP(lpar[0],lpar[3],lerr[0],lerr[1],lerr[2],
                                  cov_na,cov_ab,cov_bn)
            lpar=lpar[:-1]

        if lamb:
            ftmp=eval('lambda x,lpar: '+ formula.replace('[','lpar['))
            f=lambda x: ftmp(x,lpar)
        else:
            from ROOT import TF1
            f=TF1("f%d" % SpectralModelXML.num,formula,xmin,xmax)
            f.SetParameters(*lpar)
            # f.SetParErrors(*lerr)
            f.SetTitle("%s;E [MeV];dN/dE [cm^{-2} s^{-1} MeV^{-1}]" % sname)

        if butt:
            return f, elamb, ftype
        else:
            return f, ftype

    else:
        print "    ERROR!!!!!!! Cannot find %s" % sname
        exit(1)


def SpectralModel(resdat,sname,xmin=100.,xmax=5.e5,scale=1,lamb=False,butt=False):
    if not hasattr(SpectralModel,'num'):
        SpectralModel.num=0
    SpectralModel.num+=1

    with open(resdat) as f:
        dic=eval(f.read())

    free=dic['Free Parameters']
    cov=np.array(dic['Covariance'])

    src=dic[sname]

    lpar=[]
    if src['Spectrum'] == 'PowerLaw':
        print 'PowerLaw'
        formula,lpar,lerr=PowerLaw(src,scale)
        if butt:
            indn=free.index('%s:Prefactor'%sname)
            indg=free.index('%s:Index'%sname)
            cov_ng=cov[indn,indg]*lpar[3]*lpar[4]
            elamb=ButterflyPL(lpar[0],lerr[0],lpar[2],lerr[1],cov_ng)
        lpar=lpar[:-2]

    elif src['Spectrum'] == 'PLSuperExpCutoff':
        print 'PLSuperExpCutoff'
        formula,lpar,lerr=PLSuperExpCutoff(src,scale)
        if butt:
            ind_n=free.index('%s:Prefactor'%sname)
            ind_g=free.index('%s:Index1'%sname)
            ind_c=free.index('%s:Cutoff'%sname)
            cov_ng=cov[ind_n,ind_g]*lpar[5]*lpar[6]
            cov_gc=cov[ind_g,ind_c]*lpar[6]*lpar[7]
            cov_cn=cov[ind_c,ind_n]*lpar[7]*lpar[5]
            elamb=ButterflyPLC(lpar[0],lpar[2],lpar[3],
                               lerr[0],lerr[1],lerr[2],
                               cov_ng,cov_gc,cov_cn)
        lpar=lpar[:-3]

    elif src['Spectrum'] == 'LogParabola':
        print 'LogParabola'
        formula,lpar,lerr=LogParabola(src,scale)
        if butt:
            ind_n=free.index('%s:norm'%sname)
            ind_a=free.index('%s:alpha'%sname)
            ind_b=free.index('%s:beta'%sname)
            cov_na=cov[ind_n,ind_a]*lpar[4]
            cov_ab=cov[ind_a,ind_b]
            cov_bn=cov[ind_b,ind_n]*lpar[4]
            elamb=ButterflyLP(lpar[0],lpar[3],lerr[0],lerr[1],lerr[2],
                              cov_na,cov_ab,cov_bn)
        lpar=lpar[:-1]

    if lamb:
        ftmp=eval('lambda x,lpar: '+ formula.replace('[','lpar['))
        f=lambda x: ftmp(x,lpar)
    else:
        from ROOT import TF1
        f=TF1("f%d" % SpectralModel.num,formula,xmin,xmax)
        f.SetParameters(*lpar)
        # f.SetParErrors(*lerr)
        f.SetTitle("%s;E [MeV];dN/dE [cm^{-2} s^{-1} MeV^{-1}]" % sname)

    if butt:
        return f, elamb, src['Spectrum']
    else:
        return f, src['Spectrum']

def SpectralModelSEDdict(seddat,xmin=100.,xmax=5.e5,scale=1,lamb=False):
    if not hasattr(SpectralModelSEDdict,'num'):
        SpectralModelSEDdict.num=0
    SpectralModelSEDdict.num+=1

    with open(seddat) as f:
        dic=eval(f.read())

    spec=dic['Spectrum']
    if spec['name'] == 'PowerLaw':
        print 'PowerLaw'
        formula,lpar,lerr=PowerLaw(spec,scale)

    elif spec['name'] == 'PLSuperExpCutoff':
        print 'PLSuperExpCutoff'
        formula,lpar,lerr=PLSuperExpCutoff(spec,scale)

    elif spec['name'] == 'LogParabola':
        print 'LogParabola'
        formula,lpar,lerr=LogParabola(spec,scale)

    if lamb:
        ftmp=eval('lambda x,lpar: '+ formula.replace('[','lpar['))
        f=lambda x: ftmp(x,lpar)
    else:
        from ROOT import TF1
        f=TF1("f%d" % SpectralModelSEDdict.num,formula,xmin,xmax)
        f.SetParameters(*lpar)
        # f.SetParErrors(*lerr)
        f.SetTitle("%s;E [MeV];dN/dE [cm^{-2} s^{-1} MeV^{-1}]" % sname)

    return f, spec['name']
