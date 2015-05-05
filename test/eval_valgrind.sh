#!/bin/sh

file=$1
outfile=tmpvalgrind

valgrind --tool=memcheck --leak-check=full $1 > /dev/null 2> $outfile

deflost=`fgrep "definitely lost:" $outfile | xargs | cut -d' ' -f4 | sed "s/,//g"`
indlost=`fgrep "indirectly lost:" $outfile | xargs | cut -d' ' -f4 | sed "s/,//g"`
reachable=`fgrep "still reachable:" $outfile | xargs | cut -d' ' -f4 | sed "s/,//g"`

if [ -z "${deflost}" ]; then deflost=0; fi
if [ -z "${indlost}" ]; then indlost=0; fi
if [ -z "${reachable}" ]; then reachable=0; fi

echo $deflost $indlost $reachable

rm $outfile
