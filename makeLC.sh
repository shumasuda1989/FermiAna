#!/bin/bash
### Usage
# 1. Make a working directory and Go to the directory
# 2. Copy this script into the dir and add execute permission
# 3. Create input model xml file
# 4. Change parameters for execFermiTools.sh
# 5. Change T, DataPath, LCdir, reldir
# 6. Run this script
# 7. After this run, run DrawLC.py

curdir=$PWD
srcname=PSRJ2021
export srcname=${srcname} method=binned npix1=450 npix2=600 binsz=0.06666666667

T=(239557417 240162217 240767017 241371817 241976617 242581417 243186217 243791017 244395817 245000617) # MET
nT=${#T[*]}
digT=${#T[$((nT-1))]}
nTbin=$((nT-1))

DataPath=/home/user/Fermi/Data #absolute path to dir which contains PH&SC files
LCdir=.
reldir=.. ##relative path of $curdir relative to $workdir

[ ! -d $LCdir ] && mkdir -p $LCdir
echo ${T[@]} >$LCdir/metbin.txt
ls $DataPath/L*PH*.fits >$LCdir/list.txt || (echo error make list; exit 1)

slist="3FGL J2021.5+4026"
export slist="$slist"
export docut=1 skipmkmodel=1 skipgttsmap=1 skipres=3
# export skip=7 skipgttsmap=1 skipres=3

for ((i0=0 ; i0<nTbin ; i0++))
do
    cd $curdir
    i1=$((i0+1))
    trange=$(printf "%0${digT}d" ${T[$i0]})_$(printf "%0${digT}d" ${T[$i1]})
    echo $trange

    workdir=$LCdir/MET_$trange
    [ ! -d $workdir ] && mkdir -p $workdir
    cd $workdir || (echo error cd && exit 1)

    ln -fs $DataPath/L*SC*.fits . || (echo error ln; exit 1)
    ln -fs $reldir/${srcname}_input_model.xml || (echo error ln; exit 1)
    ln -fs ../list.txt || (echo error ln; exit 1)

    export tmin=${T[$i0]} tmax=${T[$i1]}
    $reldir/../execFermiTools.sh &>Log.log &
done
