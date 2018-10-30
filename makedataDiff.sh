#!/bin/bash

curdir=$PWD
srcname=GammaCygni
export srcname=${srcname} method=binned npix1=450 npix2=600 binsz=0.06666666667 ## MUST be same as the settings in makedataFR.sh

E=(5000 10772 23208 50000 107722 232079 500000)
nE=${#E[*]}
digE=${#E[$((nE-1))]}
nEbin=$((nE-1))

Diffdir=DataPoints
relFRdir=../.. ##relative path of $curdir relative to $workdir

targetxml="$relFRdir/${srcname}_output_model.xml"

slist="gamma Cygni,3FGL J2021.5+4026"
export slist="$slist"
export docut=1 filetobecut=$srcname.fits skip=1 skipgttsmap=1 skipres=3
# export skip=7 skipgttsmap=1 skipres=3

for ((i0=0 ; i0<nEbin ; i0++))
do
    cd $curdir
    i1=$((i0+1))
    erange=$(printf "%0${digE}d" ${E[$i0]})_$(printf "%0${digE}d" ${E[$i1]})
    echo $erange

    workdir=$Diffdir/E_$erange
    [ ! -d $workdir ] && mkdir -p $workdir
    cd $workdir || (echo error cd && exit 1)

    ln -fs $relFRdir/L*SC*.fits . || (echo error ln; exit 1)
    ln -fs $relFRdir/${srcname}.fits || (echo error ln; exit 1)
    ln -fs $relFRdir/${srcname}_binned_ltcube.fits || (echo error ln; exit 1)

    ecenter=$(echo "scale=9;sqrt(${E[$i0]}*${E[$i1]})" | bc)
    if [ -e $targetxml ] && [ ! -e makemodel.sh ]; then
	echo "python $relFRdir/../CalcSeedVal.py $targetxml \"$slist\" $ecenter \$srcmdl" >makemodel.sh
    elif [ ! -e $targetxml ]; then
	echo $targetxml does not exist; exit 1
    fi

    export emin=${E[$i0]} emax=${E[$i1]} enumbins=4
    $relFRdir/../execFermiTools.sh &>Log.log &
done
