#!/bin/bash

if [ $# -lt 9 ]; then
	echo $0 fg bg k topN outDir chunkSize jobsubcommand jobstatcommand workingDir
	exit 1;
fi

date

DEADLINE=300 #mins
fg=$1
bg=$2
k=$3
topN=$4
outDir=$5
chunkSize=$6
jobsubcommand=$7
jobstatcommand=$8
currentDir=$9

numOfLinesPerSplit=`expr $chunkSize "*" 4`

rmrie.sh -R ${outDir}

mkdir.py ${outDir}
mkdir.py ${outDir}/fg
mkdir.py ${outDir}/bg
mkdir.py ${outDir}/jobscripts

echo "spliting foreground files"
split -l ${numOfLinesPerSplit} ${fg} ${outDir}/fg/fg
echo "splitting background files"
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
	echo '#!/bin/bash' > ${outDir}/jobscripts/${fgSlaves[$i]}.sh
	echo "cd $currentDir" >> ${outDir}/jobscripts/${fgSlaves[$i]}.sh
	echo "./DKmerFinderCounter $outDir ${fgSlaves[$i]} ${fgfiles[$i]} $k" >> ${outDir}/jobscripts/${fgSlaves[$i]}.sh
done

for((i=0;i<${#bgfiles[@]};i++)); do
	echo '#!/bin/bash' > ${outDir}/jobscripts/${bgSlaves[$i]}.sh
	echo "cd $currentDir" >> ${outDir}/jobscripts/${bgSlaves[$i]}.sh
	echo "./DKmerFinderCounter $outDir ${bgSlaves[$i]} ${bgfiles[$i]} $k"  >> ${outDir}/jobscripts/${bgSlaves[$i]}.sh
done

echo '#!/bin/bash' > ${outDir}/jobscripts/master.sh
echo "cd $currentDir" >> ${outDir}/jobscripts/master.sh
echo "./DKmerFinderMaster $outDir $fgSlavesCS $bgSlavesCS $k $topN" >> ${outDir}/jobscripts/master.sh

echo "submitting jobs"
idx=0
for i in ${outDir}/jobscripts/*.sh; do
	chmod 777 $i;
	i=`abspath.py $i`
	echo "submit job command $jobsubcommand  -o ${i}.stdout -e ${i}.stderr $i | extractNumbers.py --numOfNumsPerLine 1 | head -n 1"
	jobnum[$idx]=`eval "$jobsubcommand -o ${i}.stdout -e ${i}.stderr $i | extractNumbers.py --numOfNumsPerLine 1 | head -n 1"` 
	idx=`expr $idx + 1`
done

jsubmitted=`tempfile`
echo ${jobnum[@]} | tr " " "\n" > $jsubmitted

echo "job sumitted:"
cat $jsubmitted | tr "\n" ","
echo ""

#for((i=0;i<$DEADLINE;i++)); do
while [ 1 ]; do
	sleep 1m

	jstemp=`tempfile`
	eval "$jobstatcommand" | extractNumbers.py --numOfNumsPerLine 1 > $jstemp
	numOfJobsStillRunning=`joinu.py $jsubmitted $jstemp 2> /dev/null | wc -l`
	
	if [[ $numOfJobsStillRunning == 0 ]]; then
		break
	fi
	
	now=`date`
	echo $numOfJobsStillRunning of jobs still running at $now
	cat $jstemp  | tr "\n" ","
	echo ""
	
done

#clean up

rm -R $outDir/ipm
rm -R $outDir/fg
rm -R $outDir/bg
rm -R $outDir/kmerUpdates

date
