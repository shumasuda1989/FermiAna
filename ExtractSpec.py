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
                UL=float(tmp[0])/(emax-emin)
                flag=False
    f.close()
    return E,TS,UL

def GetUL(fitdata,source):
    with open(fitdata) as f:
        dic=eval(f.read())
    e=dic['Energies']
    E=np.sqrt(e[0]*e[-1])
    TS=float(dic[source]['TS value'])
    UL=float(dic[source]['Flux UL'])/(e[-1]-e[0])
    return E,TS,UL


if len(sys.argv) != 3 and len(sys.argv) != 4:
    print 'usage: %s dir srcname [minTS(default:9)]' % sys.argv[0]
    exit(1)

dir=sys.argv[1]
prefix="_%s" % sys.argv[2]
if len(sys.argv) == 4:
    minTS=float(sys.argv[3])
else:
    minTS=9.
prefix=prefix.replace(' ','')
outputfile=dir+'/SpectrumData'+prefix+'.txt'

fline= "## MeV, cm^-2 s^-1 MeV^-1\n"
f=open(outputfile,'w')
f.write(fline)

xmllist=glob.glob('%s/E_*_*/*_output_model*.xml' % dir)

if len(xmllist) <= 0:
    print 'cannot find any xml files'
    exit(1)

for xml in xmllist:

    print xml
    srcname=sys.argv[2]
    fitfile=xml
    tmp=fitfile.split('/')
    tmp[-1]='results_mod.dat'
    fitfile='/'.join([ str(item) for item in tmp ])
    if not os.path.exists(fitfile):
        tmp[-1]='results.dat'
        fitfile='/'.join([ str(item) for item in tmp ])

    E, TS, UL = GetUL(fitfile,srcname)

    if TS < minTS:
        f.write('%s %s 0 -1 # UL!!! TS= %s\n' % (repr(E), repr(UL), repr(TS)) )
        continue

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

f.close()
# os.system('rootl SED.C\(\"%s\"\)' % outputfile)
