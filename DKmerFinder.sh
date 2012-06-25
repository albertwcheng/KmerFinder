#!/bin/bash

if [ $# -lt 7 ]; then
	echo $0 fg bg k topN outDir chunkSize jobsubcommand
	exit 1;
fi

fg=$1
bg=$2
k=$3
topN=$4
outDir=$5
chunkSize=$6
jobsubcommand=$7

numOfLinesPerSplit=`expr $chunkSize "*" 4`

rmrie.sh -R ${outDir}

mkdir.py ${outDir}
mkdir.py ${outDir}/fg
mkdir.py ${outDir}/bg
mkdir.py ${outDir}/jobscripts

split -l ${numOfLinesPerSplit} ${fg} ${outDir}/fg/fg
split -l ${numOfLinesPerSplit} ${bg} ${outDir}/bg/bg

fgfiles=(`ls ${outDir}/fg/*`)
for((i=0;i<${#fgfiles[@]};i++)); do
fgSlaves[$i]=`basename ${fgfiles[$i]}`
done
fgSlavesCS=`echo ${fgSlaves[@]} | tr " " ","`

bgfiles=(`ls ${outDir}/bg/*`)
for((i=0;i<${#bgfiles[@]};i++)); do
bgSlaves[$i]=`basename ${bgfiles[$i]}`
done
bgSlavesCS=`echo ${bgSlaves[@]} | tr " " ","`

#echo $fgfiles $fgSlavesCS $bgfiles $bgSlavesCS

for((i=0;i<${#fgfiles[@]};i++)); do
	echo "./DKmerFinderCounter $outDir ${fgSlaves[$i]} ${fgfiles[$i]} $k" > ${outDir}/jobscripts/${fgSlaves[$i]}.sh
done

for((i=0;i<${#bgfiles[@]};i++)); do
	echo "./DKmerFinderCounter $outDir ${bgSlaves[$i]} ${bgfiles[$i]} $k"  > ${outDir}/jobscripts/${bgSlaves[$i]}.sh
done

echo "./DKmerFinderMaster $outDir $fgSlavesCS $bgSlavesCS $k $topN" > ${outDir}/jobscripts/master.sh

for i in ${outDir}/jobscripts/*.sh; do
	eval "$jobsubcommand $i"
done