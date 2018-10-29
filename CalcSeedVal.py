import sys
import os
import xml.etree.ElementTree as ET
import numpy as np
import re
sys.path.append('/home/smasuda/storage/Fermi')
from BaseFuncForFermiXML import Formatter

infile=sys.argv[1]
srclist=sys.argv[2].split(',')
valScale=float(sys.argv[3])
outfile=sys.argv[4]

if not os.path.exists(infile):
    print infile, 'does not exist.'
    exit(1)
tree = ET.parse(infile)
root = tree.getroot()

for param in root.findall('.//parameter'):
    param.set('free','0')

for srcname in srclist:
    srcname=srcname.strip()
    spectrum = root.find("./source[@name='%s']/spectrum" % srcname)
    if spectrum is not None and spectrum.get('type') == 'PowerLaw':

        NodePrefactor = spectrum.find("parameter[@name='Prefactor']")
        NodeIndex     = spectrum.find("parameter[@name='Index']")
        NodeScale     = spectrum.find("parameter[@name='Scale']")
        if NodePrefactor is None or NodeIndex is None or NodeScale is None:
            print 'ERROR!!!!!!!! Cannot find prefactor or index'
            exit(1)
        Prefactor = float(NodePrefactor.get('value'))
        Index     = float(NodeIndex.get('value'))*float(NodeIndex.get('scale'))
        Scale     = float(NodeScale.get('value'))*float(NodeScale.get('scale'))

        tmp1=re.sub('e.*','',NodePrefactor.get('scale'))
        tmp2=int(re.sub('.*e','',NodePrefactor.get('scale')))

        NewPrefactor=Prefactor*pow(valScale/Scale,Index)
        ShiftIndex=int(np.floor(np.log10(NewPrefactor)))
        NewPrefactor=NewPrefactor/pow(10,ShiftIndex)
        NodePrefactor.set('value',str(NewPrefactor))
        NodePrefactor.set('scale','%se%d' % (tmp1,tmp2+ShiftIndex))
        NodePrefactor.set('free' ,'1')
        NodeScale    .set('value',str(valScale))
        NodeScale    .set('scale','1.0')

    elif spectrum is not None and spectrum.get('type') == 'PLSuperExpCutoff':

        NodePrefactor = spectrum.find("parameter[@name='Prefactor']")
        NodeIndex1    = spectrum.find("parameter[@name='Index1']")
        NodeScale     = spectrum.find("parameter[@name='Scale']")
        if NodePrefactor is None or NodeIndex1 is None or NodeScale is None:
            print 'ERROR!!!!!!!! Cannot find prefactor or index'
            exit(1)
        Prefactor=float(NodePrefactor.get('value'))
        Index    =float(NodeIndex1.get('value'))*float(NodeIndex1.get('scale'))
        Scale    =float(NodeScale.get('value'))*float(NodeScale.get('scale'))

        tmp1=re.sub('e.*','',NodePrefactor.get('scale'))
        tmp2=int(re.sub('.*e','',NodePrefactor.get('scale')))

        NewPrefactor=Prefactor*pow(valScale/Scale,Index)
        ShiftIndex=int(np.floor(np.log10(NewPrefactor)))
        NewPrefactor=NewPrefactor/pow(10,ShiftIndex)
        NodePrefactor.set('value',str(NewPrefactor))
        NodePrefactor.set('scale','%se%d' % (tmp1,tmp2+ShiftIndex))
        NodePrefactor.set('free' ,'1')
        NodeScale    .set('value',str(valScale))
        NodeScale    .set('scale','1.0')

    elif spectrum is not None and spectrum.get('type') == 'LogParabola':

        spectrum.set('type','PowerLaw')

        NodeNorm  = spectrum.find("parameter[@name='norm']")
        NodeAlpha = spectrum.find("parameter[@name='alpha']")
        NodeEb    = spectrum.find("parameter[@name='Eb']")
        NodeBeta  = spectrum.find("parameter[@name='beta']")
        if NodeNorm is None or NodeAlpha is None or NodeEb is None or NodeBeta is None:
            print 'ERROR!!!!!!!! Cannot find prefactor or index'
            exit(1)
        Norm  =float(NodeNorm.get('value'))
        Alpha =float(NodeAlpha.get('value'))*float(NodeAlpha.get('scale'))
        Eb    =float(NodeEb.get('value'))*float(NodeEb.get('scale'))
        Beta  =float(NodeBeta.get('value'))*float(NodeBeta.get('scale'))

        tmp1=re.sub('e.*','',NodeNorm.get('scale'))
        tmp2=int(re.sub('.*e','',NodeNorm.get('scale')))

        Index=Alpha+Beta*np.log(valScale/Eb)
        NewPrefactor=Norm*pow(valScale/Eb,-Index)
        NewIndex=Alpha+2*Beta*np.log(valScale/Eb)

        ShiftIndex=int(np.floor(np.log10(NewPrefactor)))
        NewPrefactor=NewPrefactor/pow(10,ShiftIndex)
        NodeNorm .set('value',str(NewPrefactor))
        NodeNorm .set('scale','%se%d' % (tmp1,tmp2+ShiftIndex))
        NodeNorm .set('free' ,'1')
        NodeNorm .set('name' ,'Prefactor')
        NodeAlpha.set('value',str(NewIndex))
        NodeAlpha.set('name' ,'Index')
        NodeAlpha.set('scale' ,'-1.0')
        NodeEb   .set('value',str(valScale))
        NodeEb   .set('scale','1.0')
        NodeEb   .set('name' ,'Scale')

        spectrum.remove(NodeBeta)

    else:
        print "    ERROR!!!!!!! Cannot find %s or Not PL or PLexpcoff or LP" % srcname
        exit(1)

tree.write(outfile)
Formatter(infile,outfile)
