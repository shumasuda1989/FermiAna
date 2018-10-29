#!/bin/bash

xml=$1
srcname=$2
valScale=$3
if [ -z "$energyL" ]; then energyL=0; fi
if [ -z "$scalePre" ]; then scalePre=1e-17; fi

IFS=$'\n'
flag=0
flagspec=0
while read line || [ -n "${line}" ]
do

    test1=$(echo "$line" | grep "<source " | grep "name=\"$srcname\"")
    if [ -n "$test1" ]; then 
	flag=1
    fi

    if [ $flag -eq 1 ]; then

	test1=$(echo "$line" | grep "<spectrum ")
	if [ -n "$test1" ]; then 
	    flagspec=1
	fi

	if [ $flagspec -eq 1 ]; then
	    test2=$(echo "$line" | grep "<parameter " | grep 'name="Scale"')
	    test3=$(echo "$line" | grep "<parameter " | grep 'name="Index"')
	    test4=$(echo "$line" | grep "<parameter " | grep 'name="Prefactor"')

	    testfs=$(echo "$line" | grep "</spectrum>")

	    if [ -n "$test2" ]; then
		echo "$line" | sed "s/value=\".*\"/value=\"$valScale\"/"
	    elif [ -n "$test3" ] && [ $energyL -ge 9999 ]; then
		echo "$line" | sed 's/free="1"/free="0"/' 
	    elif [ -n "$test4" ] && [ $energyL -ge 9999 ]; then
		echo "$line" | sed "s/scale=\"[^ ]\+\"/scale=\"$scalePre\"/"
	    elif [ -n "$testfs" ]; then
		flagspec=0
		echo "$line"
	    else
		echo "$line"
	    fi
	else
	    testfs=$(echo "$line" | grep "</source>")
	    if [ -n "$testfs" ]; then
		flag=0
	    fi
	    echo "$line"
	fi

    else 
	test5=$(echo "$line" | grep "<parameter " | grep 'free="1"')
	if [ -n "$test5" ] ; then
	    echo "$line" | sed 's/free="1"/free="0"/' 
	else
	    echo "$line"
	fi
    fi

done <$xml