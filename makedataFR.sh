#!/bin/bash

srcname=GammaCygni
date
echo

export srcname=${srcname} method=binned npix1=450 npix2=600 #mandatory
## Uncomment below line if you need
# export binsz= Ra= Dec= npix3= (for tsmap) emin= emax= enumbins= proj= ptsrc=

export srcmdlout=${srcname}_output_model.xml

## Uncomment below line if you want upper limit value (and SED) of interested sources
# export slist="gamma Cygni,3FGL J2021.5+4021" 

## For initial run. If you have a event file which should be cut by gtselect, specify 'filetobecut', else i.e. you have only L*PH??.fits photon files, you don't need 'filetobecut'.
# export docut=1 filetobecut=$srcname.fits skipmkmodel=1 skipgttsmap=1

## Uncomment below line if you want to run only pyLikelihood
# export skip=7 skipmkmodel=1 skipgttsmap=1 skipres=3
## Uncomment below line if you want to run only gttsmap
# export skip=8 skipmkmodel=1 skipgttsmap=1 skipres=3

if [ -z "$export" ]; then
    ../execFermiTools.sh
else
    source ../execFermiTools.sh #for only getting environment variables
fi

echo
date

