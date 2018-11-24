#!/usr/bin/env python

import sys
import os
import xml.etree.ElementTree as ET
import numpy as np
import glob

def GetVal(node,scale=1.):
    return  float(node.get('value'))*float(node.get('scale'))/scale
def GetErr(node,scale=1.):
    return  float(node.get('error'))*float(node.get('scale'))/scale

def oldGetUL(fitdata,source):
    f = open(fitdata,'r')

    E=-1.
    UL=-1.

    flag=False
    for line in f:
        if line.find("'%s':" % source)>=0 or flag:
            flag=True
            if line.find('TS value')>=0:
                tmp=line.split("'")
                TS=float(tmp[3])
                flag=False
            elif line.find('ph/cm^2/s')>=0:
                tmp=line.split()
                semin=tmp[3]
                emin=float(semin.split("=")[1].replace(',',''))
                semax=tmp[4]
                emax=float(semax.split("=")[1].replace(',',''))
                E=np.sqrt(emin*emax)
                UL=float(tmp[0])/(emax-emin) ## WRONG calc. Do Not use this.
                flag=False
    f.close()
    return E,TS,UL

def GetUL(dic,source):
    e=dic['Energies']
    E=np.sqrt(e[0]*e[-1])
    TS=float(dic[source]['TS value'])
    # UL=float(dic[source]['Flux UL'])/(-e[0]+e[-1]) ## wrong calc.
    UL=float(dic[source]['dNdE UL'])
    return E,TS,UL

def WriteSpecFromXML(xml,srcname,TS,f):
    tree = ET.parse(xml)
    root = tree.getroot()

    spectrum = root.find("./source[@name='%s']/spectrum" % srcname)
    if spectrum is not None and spectrum.get('type') == 'PowerLaw':

        NodePrefactor = spectrum.find("parameter[@name='Prefactor']")
        NodeScale     = spectrum.find("parameter[@name='Scale']")
        if NodePrefactor is None or NodeScale is None:
            print 'ERROR!!!!!!!! Cannot find prefactor or scale nodes'
            exit(1)
        Prefactor = GetVal(NodePrefactor)
        EPre      = GetErr(NodePrefactor)
        Scale     = GetVal(NodeScale)

        f.write('%s %s 0 %s # TS= %s\n' % (str(Scale), str(Prefactor), str(EPre), repr(TS)) )

    elif spectrum is not None and spectrum.get('type') == 'PLSuperExpCutoff':

        NodePrefactor = spectrum.find("parameter[@name='Prefactor']")
        NodeScale     = spectrum.find("parameter[@name='Scale']")
        NodeCutoff    = spectrum.find("parameter[@name='Cutoff']")
        NodeIndex2    = spectrum.find("parameter[@name='Index2']")

        if NodePrefactor is None or NodeScale is None:
            print 'ERROR!!!!!!!! Cannot find prefactor or scale nodes'
            exit(1)
        Prefactor = GetVal(NodePrefactor)
        Cutoff = GetVal(NodeCutoff)
        EPre      = GetErr(NodePrefactor)
        ECoff     = GetErr(NodeCutoff)
        Scale     = GetVal(NodeScale)

        Value = Prefactor * np.exp(-Scale/Cutoff)
        Error = EPre * np.exp(-Scale/Cutoff)
        # Error = np.exp(Scale/Cutoff) \
        #     *sqrt(EPre**2 + (ECoff*Prefactor*Scale/(Cutoff**2))**2 \
        #               -2*Cov*Prefactor*Scale/(Cutoff**2))

        f.write('%s %s 0 %s # TS= %s\n' % (str(Scale), str(Value), str(Error), repr(TS)) )

    else:
        print "    ERROR!!!!!!! Cannot find %s or Not power law" % srcname
        exit(1)

    
def WriteSpec(dic,srcname,f):
    src=dic[srcname]
    if src.get('Spectrum','') != 'PowerLaw':
        print ' ERROR! Not Power Law'
        exit(1)
    else:
        E=float(src['Scale'])*float(src['scale Scale'])
        flux=[float(x) for x in src['Prefactor'].split('+/-')]
        v=flux[0]*float(src['scale Prefactor'])
        e=flux[1]*float(src['scale Prefactor'])
        TS=float(src['TS value'])
        f.write('%s %s 0 %s # TS= %s\n' % (str(E), str(v), str(e), repr(TS)) )


if __name__ == '__main__':

    minTSdef=4.
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", help="Name of directory which includes E_* directories")
    parser.add_argument("srcnamelist",help="List of names of sources written in xml model file which should be sepalated by comma ','.")
    parser.add_argument("-t","--minTS",metavar="minTS",type=float,default=minTSdef,help="Minimum TS value under which data points are treated as upper limits (default: %(default)s)")
    parser.add_argument("-i","--infile",default="results.dat",help="Likelihood data file (default:  %(default)s)")
    parser.add_argument("-o","--outfile",default='Spec_%s.txt',help="Output file name. '%%s' is replaced by srcname (default: %(default)s)")
    args = parser.parse_args()

    dir=args.dir
    srcnamelist=args.srcnamelist.split(',')

    for srcname in srcnamelist:

        print srcname
        prefix=srcname.replace(' ','')
        if '%s' in args.outfile:
            outputfile=dir+'/'+(args.outfile % prefix)
        else:
            outputfile=dir+'/'+args.outfile

        fline= "## MeV, cm^-2 s^-1 MeV^-1\n"
        f=open(outputfile,'w')
        f.write(fline)

        reslist=glob.glob('%s/E_*_*/%s' % (dir,args.infile))

        if len(reslist) <= 0:
            print 'cannot find any result files'
            exit(1)

        for res in reslist:

            print res

            with open(res) as fin:
                dic=eval(fin.read())
            E, TS, UL = GetUL(dic,srcname)

            if TS < args.minTS:
                f.write('%s %s 0 -1 # UL!!! TS= %s\n'%(repr(E), repr(UL), repr(TS)))
                continue
            else:
                # WriteSpecFromXML(xml,srcname,TS,f)
                WriteSpec(dic,srcname,f)

        f.close()

