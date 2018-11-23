#!/bin/bash

srcname=GammaCygni
date
echo

export srcname=${srcname} method=binned npix1=450 npix2=600 #mandatory
## Add below parameters if you need
# export binsz=(default:0.2) Ra= Dec= npix3=(for tsmap) emin= emax= enumbins= proj= ptsrc=

export srcmdlout=${srcname}_output_model.xml

## Uncomment below line and Add the names of sources you are interested in if you want upper limit values (and SED) of them
# export slist="gamma Cygni,3FGL J2021.5+4021" 

export docut=1 skipmkmodel=1 skipgttsmap=1
## For initial run. If you have a event file which should be cut by gtselect in advance, specify 'filetobecut', else i.e. you have only L*PH??.fits photon files, you don't need 'filetobecut'. Put or link them into work dir or Make list.txt

## Uncomment below line if you want to run only pyLikelihood
# export skip=7 skipmkmodel=1 skipgttsmap=1 skipres=3
## Uncomment below line if you want to run only gttsmap
# export skip=8 skipmkmodel=1 skipres=3

if [ -z "$export" ]; then
    ../execFermiTools.sh
else
    source ../execFermiTools.sh #only for getting environment variables
fi

echo
date

