import sys
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

lcdir=sys.argv[1]
src=sys.argv[2]
reslist=glob('%s/MET_*/results.dat'%lcdir)
tbin=np.loadtxt('%s/metbin.txt'%lcdir)

print tbin

tref=(tbin[:-1]+tbin[1:])/2.

flux=[]
eflx=[]

for res in reslist:
    print res
    with open(res) as f:
        dic=eval(f.read())
    fltmp=[float(x) for x in dic[src]['Flux'].split('+/-')]
    flux.append(fltmp[0])
    eflx.append(fltmp[1])

ax=plt.subplot(111)
ax.errorbar(tref,flux,yerr=eflx, fmt='o',linestyle='None',capsize=4)

ax.set_xlabel('MET')
ax.set_ylabel('Flux [ph cm-2 s-1]')

plt.show()
