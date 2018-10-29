#!/usr/bin/env python
import sys
import numpy as np

def ReadResult(fname):

    with open(fname) as f:
        lines = f.readlines()

    result={}
    flagdf=False
    for i,line in enumerate(lines):
        
        if flagdf:
            flagdf=False
            result[src]['Diff Flux']=np.fromstring(line,sep=' ')
            continue
        elif line == 'Flux UL:\n':
            indul=i
            break
        elif line == '},\n' or line == '\n':
            continue
        else:
            s=line.split("'")
            if s[2] == ': {\n':
                src=s[1]
                result[src]={}
            else:
                if s[1] == 'Diff Flux':
                    flagdf=True
                    continue
                    
                result[src][s[1]]=np.fromstring(s[3],sep='+-')

    ul={}
    for line in lines[indul+1:]:
        if line == '\n':
            continue
        elif line.endswith("':\n"):
            src=line.split("'")[1]
        else:
            ul[src]=float(line.split()[0])

    return result, ul



if __name__ == '__main__':

    print ReadResult(sys.argv[1])
