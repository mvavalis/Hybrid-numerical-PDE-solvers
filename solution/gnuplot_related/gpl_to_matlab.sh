#!/bin/bash

if [ -n "$1"  -a  -f "$1" ]
then
	outfile=`echo "$1" | sed 's/\./_/g'`
	
	echo "x = [" > "$outfile.m"
	cat $1 | grep -v '#' | grep -v '^$' | sed 's/$/;/' >> "$outfile.m"
	echo "];" >> "$outfile.m"
else
	echo "need input file :P"
fi