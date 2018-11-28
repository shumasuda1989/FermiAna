import sys
from glob import glob
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

if len(sys.argv) not in [3,4]:
    print 'Usage: %s LCdir srcname [minTS:25.]'%sys.argv[0]
    exit()

ver=matplotlib.__version__.split('.')
ver=float('.'.join(ver[:2]))
lims={'uplims':True} if ver>=1.4 else {'lolims':True}

lcdir=sys.argv[1]
src=sys.argv[2]
minTS=25. if len(sys.argv)==3 else float(sys.argv[3])
reslist=glob('%s/MET_*/results.dat'%lcdir)
tbin=np.loadtxt('%s/metbin.txt'%lcdir)

print tbin

tref=(tbin[:-1]+tbin[1:])/2.

flux=[]
eflx=[]
ul  =[]
ts  =[]

for res in reslist:
    print res
    with open(res) as f:
        dic=eval(f.read())
    fltmp=[float(x) for x in dic[src]['Flux'].split('+/-')]
    flux.append(fltmp[0])
    eflx.append(fltmp[1])
    ul.append(float(dic[src]['Flux UL']))
    ts.append(float(dic[src]['TS value']))

flux=np.array(flux)
eflx=np.array(eflx)
ul  =np.array(ul)
ts  =np.array(ts)

whul=np.where(ts<minTS)
whnu=np.where(ts>=minTS)

ax=plt.subplot(111)
col='b'
ax.errorbar(tref[whnu],flux[whnu],yerr=eflx[whnu], 
            fmt='o',linestyle='None',capsize=4,color=col)
ax.errorbar(tref[whul],ul[whul]  ,yerr=0.3*ul[whul], 
            fmt='', linestyle='None',capsize=4,color=col,**lims)

ax.set_xlabel('MET')
ax.set_ylabel('Flux [ph cm-2 s-1]')

plt.show()
