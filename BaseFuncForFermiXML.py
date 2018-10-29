import os

def Formatter(infile,outfile,flag=False):

    if not os.path.exists(infile):
        print infile, 'does not exist.'
        exit(1)

    if not os.path.exists(outfile):
        print outfile, 'does not exist.'
        exit(1)

    f = open(outfile,'r')
    xmlmaintext = f.read()
    if flag:
        xmlmaintext = xmlmaintext.replace(' />','/>')
    f.close()

    f = open(infile,'r')
    firstline = f.readline()
    f.close()

    f = open(outfile,'w')
    f.write(firstline+xmlmaintext)
    f.close()

