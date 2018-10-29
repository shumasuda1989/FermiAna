import sys
import os
import xml.etree.ElementTree as ET
import re

infile=sys.argv[1]
outfile="out.txt"

if not os.path.exists(infile):
    print infile, 'does not exist.'
    exit(1)
tree = ET.parse(infile)
root = tree.getroot()

srcname='3FGL J2021.5+4026'
spectrum = root.find("./source[@name='%s']/spectrum" % srcname)

NodePrefactor = spectrum.find("parameter[@name='Prefactor']")
NodeIndex1    = spectrum.find("parameter[@name='Index1']")
NodeCutoff    = spectrum.find("parameter[@name='Cutoff']")
NodeScale     = spectrum.find("parameter[@name='Scale']")
NodeIndex2    = spectrum.find("parameter[@name='Index2']")
if NodePrefactor is None or NodeIndex1 is None or NodeScale is None \
        or NodeCutoff is None or NodeIndex2 is None:
    print 'ERROR!!!!!!!! Cannot find prefactor or index'
    exit(1)
Prefactor = NodePrefactor.get('value')
Index1    = NodeIndex1.get('value')
Cutoff    = NodeCutoff.get('value')
Scale     = NodeScale.get('value')
Index2    = NodeIndex2.get('value')
EPrefactor = NodePrefactor.get('error')
EIndex1    = NodeIndex1.get('error')
ECutoff    = NodeCutoff.get('error')
EScale     = NodeScale.get('error')
EIndex2    = NodeIndex2.get('error')

tmp=re.sub(r'.*e','e', NodePrefactor.get('scale'))

sep=' '
f=open(outfile,"a")
f.write(Prefactor+tmp+sep+EPrefactor+tmp+sep+Index1+sep+EIndex1+sep+Cutoff+sep+ECutoff+'\n')
# f.write(Prefactor+tmp+sep+'None'+sep+Index1+sep+'None'+sep+Cutoff+sep+'None\n')
f.close()

print Prefactor+tmp+',-'+Index1+','+Scale+','+Cutoff+','+Index2


