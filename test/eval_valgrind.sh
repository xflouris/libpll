#!/bin/sh
#
#    Copyright (C) 2015 Diego Darriba
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as
#    published by the Free Software Foundation, either version 3 of the
#    License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    Contact: Diego Darriba <Diego.Darriba@h-its.org>,
#    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
#    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
#

file=$1
timestamp=$2
shift
shift
args=$@

outfile=tmpvalgrind

valgrind --tool=memcheck --leak-check=full ${file} ${args} > /dev/null 2> ${outfile}

deflost=`fgrep "definitely lost:" ${outfile} | xargs | cut -d' ' -f4 | sed "s/,//g"`
indlost=`fgrep "indirectly lost:" ${outfile} | xargs | cut -d' ' -f4 | sed "s/,//g"`
reachable=`fgrep "still reachable:" ${outfile} | xargs | cut -d' ' -f4 | sed "s/,//g"`
errors=`fgrep "Invalid" ${outfile}`

if [ -z "${deflost}" ]; then deflost=0; fi
if [ -z "${indlost}" ]; then indlost=0; fi
if [ -z "${reachable}" ]; then reachable=0; fi
if [ -z "${errors}" ]; then errors=0; else errors=1; fi

echo ${deflost} ${indlost} ${reachable} ${errors}

if [ "$((errors+deflost+indlost+reachable))" -gt "0" ]; then
  filename="${file##*/}"
  mv ${outfile} result/valgrind_${filename}_${timestamp}
else
  rm ${outfile}
fi
