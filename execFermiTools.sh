#!/bin/bash

if [ -z "$FERMI_DIR" ]; then
    echo Error: FERMI_DIR is not set >&2 ; exit 1
else
    source $FERMI_DIR/fermi-init.sh || exit 1
    export PYTHONNOUSERSITE=yes
fi

if [ -z "$export" ]; then
    VERDATE=$(stat -c %Y $0) #epoch time
    echo -e "\e[32m-- execFermiTools.sh version "$(date +"%Y/%m/%d %H:%M:%S" -d @$VERDATE)" ($VERDATE)\e[m"
    olddir=$(dirname $0)/.oldexecFermiTools
    if [ -d $olddir ]; then mkdir -p $olddir; fi
    CPFILE=$olddir/execFermiTools.sh.$VERDATE
    if ! [ -e "$CPFILE" ]; then cat $0 >$CPFILE; fi
fi
echo

SECONDS=0
echo PATH=$PATH; echo

: ${skip:=0} ${skipres:=0} ${skipgtlike:=0} ${skipgttsmap:=0}
: ${docut:=0} ${cutonly:=0} # default value(=0) is set if no value is specified

: ${srcname:?} ${method:?} ${npix1:?} ${npix2:?}
# check if values in these vars are set
if [ -z "$srcrad" ]; then
    echo  srcrad is set to 20 deg 
    srcrad=20
fi
METHOD=${method^^}
if [ $skip -lt 8 ] && [ $skipgtlike -eq 0 ] && [ -z "$slist" ]; then
    echo "\e[31mslist is not set. ULs will not be calculated.\e[m"
fi
if [ $skip -lt 9 ] && [ $skipgttsmap -eq 0 ]; then
    : ${npix3:?}
fi
if [ -z "$binsz" ]; then 
    echo pixel size binsz is set to 0.2 deg/pix 
    binsz=0.2
fi
if [ -z "$binszts" ]; then 
    echo pixel size for gttsmap '"binszts"' is set to the same as binsz 
    binszts=$binsz
fi
if [ -z "$inv_pixsize" ]; then 
    echo pixel size 1/inv_pixsize is set to 0.2 deg/pix 
    inv_pixsize=5 # =1/0.2
fi
: ${proj:=AIT} ${ptsrc:=yes} ${coordsys:=CEL}
: ${USE_BL_EDISP:=true} ${refit:=no} ${plot:=no} ${optimizer:=NEWMINUIT}

if [ "${USE_BL_EDISP}" == "true" ]; then
    export USE_BL_EDISP
else
    unset USE_BL_EDISP
fi

DSSkey(){
    if [ -n "$DSSkeyfile" ]; then
	str=$(gtvcut table=EVENTS $DSSkeyfile | grep $1)
	dstyp=${str/: $1/}
	str=$(gtvcut table=EVENTS $DSSkeyfile | grep ${dstyp/TYP/VAL})
	str=$(echo $str | sed 's/DSVAL.: \(.*\)/\1/')
	echo $str
    fi
}

: ${prefix:=${srcname}_${method}}
: ${evfile:=${prefix}.fits}
echo evfile is set to $evfile; echo
if [ ! -e "$evfile" ]; then
    echo Warning: currently evfile does not exist
    DSSkeyfile=$filetobecut
    if [ ! -e "$filetobecut" ]; then
	DSSkeyfile=$(ls L*_PH00.fits)
    fi
else
    DSSkeyfile=$evfile
fi

if [ -z "$Ra" ] || [ -z "$Dec" ]; then
    VAL=$(DSSkey "POS(RA,DEC)")
    Ra=$(echo $VAL | sed 's/CIRCLE(\(.*\))/\1/' | awk -F, '{print $1;}')
    Dec=$(echo $VAL | sed 's/CIRCLE(\(.*\))/\1/' | awk -F, '{print $2;}')
    if [ -z "$Ra" ] || [ -z "$Dec" ]; then
	echo "specify (Ra,Dec) or DSSkey file" >&2
	exit 1
    fi
fi
if [ -z "$emin" ]; then
    VAL=$(DSSkey ENERGY)
    emin=$(echo $VAL | awk -F: '{print $1;}')
    if [ -z "$emin" ]; then
	echo "specify emin or DSSkey file" >&2
	exit 1
    fi
fi
if [ -z "$emax" ]; then
    VAL=$(DSSkey ENERGY)
    emax=$(echo $VAL | awk -F: '{print $2;}')
    if [ -z "$emax" ]; then
	echo "specify emax or DSSkey file" >&2
	exit 1
    fi
fi
if [ -z "$enumbins" ]; then
    enumbins=$(echo "(l($emax)-l($emin))*10/l(10)+0.5" | bc -l)
    enumbins=$(echo "$enumbins/1" | bc)
fi
echo srcname=$srcname method=$method Ra=$Ra Dec=$Dec srcrad=$srcrad npix1=$npix1 npix2=$npix2 npix3=$npix3 binsz=$binsz binszts=$binszts emin=$emin emax=$emax enumbins=$enumbins inv_pixsize=$inv_pixsize proj=$proj ptsrc=$ptsrc coordsys=$coordsys slist="'$slist'"
echo

if [ "$nosrc" == "yes" ] || [ "$nosrc" == "1" ]; then
# For calculating TS value of extended source
    : ${srcname2:=${srcname}_nsrc}
    : ${prefix2:=${srcname2}_${method}}
else
    : ${srcname2:=$srcname}
    : ${prefix2:=$prefix}
fi

: ${scfile:=$(ls L*_SC00.fits)}
if [ -z "$scfile" ]; then echo scfile not found; exit 1; fi
: ${cmap:=${prefix}_cmap.fits}
: ${ccube:=${prefix}_ccube.fits}
: ${lvtime:=${prefix}_ltcube.fits}
: ${expmap:=${prefix}_expmap.fits}
: ${bexpcube:=${prefix}_expcube.fits}
: ${srcmdlin:=${srcname2}_input_model.xml}
: ${srcmdlout:=${srcname2}_output_model_${emin}-${emax}.xml}
: ${srcmap:=${prefix2}_srcmaps.fits}
if [ "$chngnm" == "yes" ] || [ "$chngnm" == "1" ]; then
    opres="results=${results:=${prefix2}_results.dat}"
    opplo="specfile=${specfile:=${prefix2}_counts_spectra.fits}"
    : ${fitso:=${prefix2}_fit_data.txt}
else
    ## 'chngnm', 'opres' and 'opplo' are no longer obsolete. 
    ## Simply use 'results', 'specfile' and 'fitso' to overwrite them.
    : ${fitso:=fit_data.txt}
fi
: ${srcmdlfix:=${srcmdlout/.xml/_fix.xml}}
: ${tsmap:=${prefix2}_tsmap.fits}

echo evfile=$evfile scfile=$scfile cmap=$cmap ccube=$ccube lvtime=$lvtime expmap=$expmap bexpcube=$bexpcube srcmdlin=$srcmdlin srcmdlout=$srcmdlout srcmap=$srcmap fitso=$fitso results=$results specfile=$specfile tsmap=$tsmap 
echo


export method METHOD Ra Dec srcrad binsz binszts emin emax enumbins inv_pixsize proj ptsrc coordsys optimizer refit plot slist
export evfile scfile bexpcube lvtime srcmdlin srcmdlout srcmap cmap ccube tsmap fitso
if [ -n "$results" ]; then export results; fi
if [ -n "$specfile" ]; then export specfile; fi

if [ -n "$export" ]; then
    echo Export done!;  echo
    return
fi

error(){
    echo -e "\e[31merror occured in $1!\e[m" >&2
    exit 1
}

if [ $docut -ne 0 ] || [ $cutonly -ne 0 ] ; then

    if [ -n "$filetobecut" ]; then 
	echo filetobecut is set to $filetobecut
    else
	filetobecut="@list.txt"
	if [ ! -e "list.txt" ]; then ls L*_PH*.fits >list.txt; fi
    fi
    gtselect infile=${filetobecut} outfile=${prefix}_filtered.fits \
	ra=INDEF dec=INDEF rad=INDEF evclass=128 evtype=3 \
	tmin=INDEF tmax=INDEF emin=${emin} emax=${emax} zmax=90 \
	|| error gtselect
    echo "gtselect done. elapsed time: $SECONDS s"
    sleep 3; echo; echo

    gtmktime scfile=${scfile} \
	filter="(DATA_QUAL>0)&&(LAT_CONFIG==1)" roicut=no \
	evfile=${prefix}_filtered.fits outfile=${evfile} || error gtmktime
    echo "gtmktime done. elapsed time: $SECONDS s"

    if [ $cutonly -ne 0 ]; then 
	echo
	echo only cut done
	exit 0
    fi
    sleep 3; echo; echo
fi


if [ $skip -lt 1 ]; then
    gtltcube zmax=90 evfile=${evfile} scfile=${scfile} outfile=${lvtime} \
	dcostheta=0.025 binsz=1 || error gtltcube
    echo "gtltcube done. elapsed time: $SECONDS s"
else echo gtltcube was skipped; fi 

sleep 3; echo; echo

if [ $skip -lt 2 ] && [ $METHOD = "UNBINNED" ]; then
    gtexpmap evfile=${evfile} scfile=${scfile} expcube=${lvtime} \
	outfile=${expmap} irfs=P8R2_SOURCE_V6 srcrad=${srcrad} \
	nlong=$((srcrad*inv_pixsize)) nlat=$((srcrad*inv_pixsize)) \
	nenergies=${enumbins} || error gtexpmap
    echo "gtexpmap done. elapsed time: $SECONDS s"
else echo gtexpmap was skipped; fi

sleep 3; echo; echo

if [ $skip -lt 3 ]; then
    gtbin algorithm=CMAP evfile=${evfile} scfile=NONE outfile=${cmap} \
	nxpix=${npix2} nypix=${npix2} binsz=${binsz} coordsys=${coordsys} \
	xref=${Ra} yref=${Dec} axisrot=0 proj=${proj} || error "gtbin CMAP"
    echo "gtbin CMAP done. elapsed time: $SECONDS s"
else echo gtbin CMAP was skipped; fi 

sleep 3; echo; echo

if [ $skip -lt 4 ]; then
    gtbin algorithm=CCUBE evfile=${evfile} outfile=${ccube} scfile=NONE \
	nxpix=${npix1} nypix=${npix1} binsz=${binsz} coordsys=${coordsys} \
	xref=${Ra} yref=${Dec} axisrot=0 proj=${proj} ebinalg=LOG \
	emin=${emin} emax=${emax} enumbins=${enumbins} || error "gtbin CCUBE"
    echo "gtbin CCUBE done. elapsed time: $SECONDS s"
else echo gtbin CCUBE was skipped; fi 

sleep 3; echo; echo

if [ $skip -lt 5 ]; then
    gtexpcube2 infile=${lvtime} cmap=none outfile=${bexpcube} \
	irfs=P8R2_SOURCE_V6 nxpix=${npix2} nypix=${npix2} binsz=${binsz} \
	coordsys=${coordsys} xref=${Ra} yref=${Dec} \
	axisrot=0 proj=${proj} ebinalg=LOG \
	emin=${emin} emax=${emax} enumbins=${enumbins} || error gtexpcube2
    echo "gtexpcube2 done. elapsed time: $SECONDS s"
else echo gtexpcube2 was skipped; fi 

sleep 3; echo; echo

if [ $skip -lt 6 ] && [ -z "$skipmkmodel" ]; then
    evfile=${evfile} srcmdl=${srcmdlin} bash ${mkmdl:=makemodel.sh} \
	|| error makemodel
    echo "makemodel done. elapsed time: $SECONDS s"
else echo makemodel was skipped; fi 

sleep 3; echo; echo

if [ $skip -lt 7 ]; then
    if [ $ptsrc == "no" ]; then
	gtsrcmaps ptsrc=${ptsrc} emapbnds=no irfs=CALDB scfile=${scfile} \
	    expcube=${lvtime} cmap=${ccube} srcmdl=${srcmdlin} \
	    bexpmap=${bexpcube} outfile=${srcmap/.fits/_dif.fits} \
	    || error gtsrcmaps
	$(dirname $0)/RemoveSource.py ${srcmdlin} diffuse tmptmp.xml
	gtsrcmaps emapbnds=no irfs=CALDB scfile=${scfile} \
	    expcube=${lvtime} cmap=${ccube} srcmdl=tmptmp.xml \
	    bexpmap=${bexpcube} outfile=${srcmap} || error gtsrcmaps
	$(dirname $0)/MergeSrcmap.py ${srcmap} ${srcmap/.fits/_dif.fits}
	rm tmptmp.xml
    else
	gtsrcmaps emapbnds=no irfs=CALDB scfile=${scfile} \
	    expcube=${lvtime} cmap=${ccube} srcmdl=${srcmdlin} \
	    bexpmap=${bexpcube} outfile=${srcmap} || error gtsrcmaps
    fi	
    echo "gtsrcmaps done. elapsed time: $SECONDS s"
else echo gtsrcmaps was skipped; fi 

sleep 3; echo; echo

if [ $skip -lt 8 ] && [ $skipgtlike -eq 0 ] ; then
    #gtlike refit=${refit} plot=${plot} sfile=${srcmdlout} statistic=${METHOD} \
    # cmap=${srcmap} bexpmap=${bexpcube} expcube=${lvtime} srcmdl=${srcmdlin} \
    # irfs=CALDB optimizer=${optimizer} \
    # ${opres} ${opplo} | tee ${fitso}

    $(dirname $0)/Likelihood.py 2>&1 | tee ${fitso}
    if [ ${PIPESTATUS[0]} -ne 0 ]; then 
	error gtlike
    else
	echo "gtlike done. elapsed time: $SECONDS s"
    fi
else echo gtlike was skipped; fi 

sleep 3; echo; echo

if [ $skip -lt 9 ] && [ $skipgttsmap -eq 0 ] ; then
    $(dirname $0)/FixParamInModel.sh ${srcmdlout} > ${srcmdlfix}
    gttsmap statistic=${METHOD} evfile=${evfile} scfile=${scfile} \
	bexpmap=${bexpcube} expcube=${lvtime} srcmdl=${srcmdlfix} \
	cmap=${ccube} outfile=${tsmap} irfs=CALDB optimizer=${optimizer} \
	nxpix=${npix3} nypix=${npix3} binsz=${binszts} xref=${Ra} yref=${Dec} \
	coordsys=${coordsys} proj=${proj} || error gttsmap
    echo "gttsmap done. elapsed time: $SECONDS s"
else echo gttsmap was skipped; fi



sleep 3; echo; echo

if [ $skipres -lt 1 ]; then
    gtbin evfile=${evfile} \
	scfile=NONE outfile=${prefix}_cmap_small.fits \
	algorithm=CMAP nxpix=${npix1} nypix=${npix1} binsz=${binsz} \
	coordsys=${coordsys} \
	xref=${Ra} yref=${Dec} axisrot=0 proj=${proj}  || error "gtbin CMAP"
    echo "gtbin CMAP small done. elapsed time: $SECONDS s"
else echo gtbin CMAP small was skipped; fi 

sleep 3; echo; echo

if [ $skipres -lt 2 ]; then
    gtmodel srcmaps=${srcmap} srcmdl=${srcmdlout} \
	outfile=${srcname2}_model_map.fits irfs=CALDB expcube=${lvtime} \
	bexpmap=${bexpcube} || error gtmodel
    echo "gtmodel done. elapsed time: $SECONDS s"
else echo gtmodel was skipped; fi 

sleep 3; echo; echo

if [ $skipres -lt 3 ]; then
    farith ${prefix}_cmap_small.fits ${srcname2}_model_map.fits \
	${srcname2}_residual.fits SUB || error farith
    echo "farith done. elapsed time: $SECONDS s"
else echo farith was skipped; fi


echo
echo -e "\e[34;1m DONE!!  End time: $SECONDS s\e[m"