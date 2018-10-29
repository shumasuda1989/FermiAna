#!/bin/bash

xml=$1
Scale=$2

IFS=$'\n'

flagspec=0
while read line || [ -n "${line}" ]
do

    test1=$(echo "$line" | grep "<spectrum ")
    if [ -n "$test1" ]; then 
    	flagspec=1
    fi

    if [ $flagspec -eq 1 ]; then
	test2=$(echo "$line" | grep "<parameter " | grep 'name="Prefactor"')
	test3=$(echo "$line" | grep "<parameter " | grep 'name="norm"')
	test4=$(echo "$line" | grep "<parameter " | grep 'name="Integral"')
	test5=$(echo "$line" | grep "<parameter " | grep 'name="Normalization"')

	if [ -n "$test2" ] || [ -n "$test3" ] || [ -n "$test4" ] || [ -n "$test5" ]; then
	    str=$(echo "$line" | sed "s/.*scale=\(......\).*/\1/")
	    # echo $str
	    teststr1=$(echo "$str" | grep "\"1\"")
	    teststr2=$(echo "$str" | grep "\"1e")
	    if [ -z "$teststr1" ] && [ -z "$teststr2" ]; then
		echo error occured in replaceing; exit 1
	    fi
	    echo "$line" | sed "s/scale=\"1/scale=\"${Scale}/"
	    flagspec=0
	else
	    echo "$line"
	fi
    else
	echo "$line"
    fi

done <$xml