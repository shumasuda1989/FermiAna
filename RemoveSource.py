#!/usr/bin/env python

import sys
import xml.etree.ElementTree as ET

def RemoveSource(fname,srcnamelist,outname):

    srclist=srcnamelist.split(',')

    xml=fname
    print xml
    print srclist

    tree = ET.parse(xml)
    root = tree.getroot()

    if srclist[0] == 'Diffuse' or srclist[0] == 'diffuse':
        srcellist = root.findall('./source[@type="DiffuseSource"]')
        for srcel in srcellist:
            root.remove(srcel)
    else:
        for srcname in srclist:
            srcel = root.find("./source[@name='%s']" % srcname)
            if srcel is not None:
                root.remove(srcel)
            else:
                print "    ERROR!!!!!!! Cannot find %s" % srcname
                exit(1)

    tree.write(outname)

if __name__ == '__main__':

    if len(sys.argv) != 4 :
        print 'usage: python %s infile srcnamelist outname' % sys.argv[0]
        exit(1)

    fname=sys.argv[1]
    srcnamelist=sys.argv[2]
    outname=sys.argv[3]
    RemoveSource(fname,srcnamelist,outname)
