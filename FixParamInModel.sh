#!/bin/bash

xml=$1
freeBG=$2
IFS=$'\n'

flagBG=0

while read line || [ -n "${line}" ]
do

    if [ -n "$freeBG" ]; then
	test1=$(echo "$line" | grep "<source " | grep 'name="gll_iem_v')
	test2=$(echo "$line" | grep "<source " | grep 'name="iso_P')
	if [ -n "$test1" ] || [ -n "$test2" ]; then 
	    flagBG=1
	else

	    testEND=$(echo "$line" | grep "</source>")
	    if [ -n "$testEND" ]; then 
		flagBG=0
	    fi
	fi
    fi

    test=$(echo "$line" | grep "<parameter " | grep 'free="1"')
    if [ -n "$test" ] && [ $flagBG -eq 0 ]; then
	echo "$line" | sed 's/free="1"/free="0"/' 
    else
	echo "$line"
    fi


done <$xml
