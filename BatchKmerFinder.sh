#!/bin/bash

:<<'COMMENT'

concentration file (inputSettingFile)

conc<tab>foreground<tab>background<tab>
COMMENT

if [ $# -lt 7 ]; then
	
	echo $0 "inputSettingFile k topN chunkJobs outDir bsubcommand bjobcommand"
	exit 1
fi

currentDir=`pwd`

echo "currentDir=$currentDir"

inputSettingFile=$1
k=$2
topN=$3
chunkJobs=$4
outDir=$5
bsubcommand=$6
bjobcommand=$7

mkdir.py $outDir
mkdir.py $outDir/tmp/

concs=(`cuta.py -f1 $inputSettingFile`) #first column
fgfiles=(`cuta.py -f2 $inputSettingFile`) #second
bgfiles=(`cuta.py -f3 $inputSettingFile`)

for((i=0;i<${#concs[@]};i++)); do
	conc=${concs[$i]}
	fgfile=${fgfiles[$i]}
	bgfile=${bgfiles[$i]}
	
	echo running fgfile=$fgfile bgfile=$bgfile for conc=$conc
	bash DKmerFinder.sh $fgfile $bgfile $k $topN $outDir/$conc/ $chunkJobs "$bsubcommand" "$bjobcommand" $currentDir
	
done

jfiles=""
for((i=0;i<${#concs[@]};i++)); do
	conc=${concs[$i]}
	jfiles="$jfiles $outDir/$conc/kmers.txt"
	
done

#kmer    enrichment      control_count   experim_count

echo $jfiles

echo "kmer" > $outDir/tmp/concheaderColumnar.txt
echo ${concs[@]} | tr " " "\n" >> $outDir/tmp/concheaderColumnar.txt

rmrie.sh $outDir/kmers.merged.txt
eval "multijoinu.sh -fNA $outDir/kmers.merged.txt $jfiles"  #joint

columnSelector="enrichment"
cuta.py -f.kmer,@$columnSelector $outDir/kmers.merged.txt > $outDir/tmp/t1.txt
matrixTranspose.py $outDir/tmp/t1.txt > $outDir/tmp/t1.txt.t
cuta.py -f2-_1 $outDir/tmp/t1.txt.t > $outDir/tmp/t1.txt.t.f2
paste  $outDir/tmp/concheaderColumnar.txt $outDir/tmp/t1.txt.t.f2 > $outDir/tmp/t1.txt.t.swapped
matrixTranspose.py $outDir/tmp/t1.txt.t.swapped > $outDir/kmers.merged.$columnSelector.txt

columnSelector="control_count"
cuta.py -f.kmer,@$columnSelector $outDir/kmers.merged.txt > $outDir/tmp/t1.txt
matrixTranspose.py $outDir/tmp/t1.txt > $outDir/tmp/t1.txt.t
cuta.py -f2-_1 $outDir/tmp/t1.txt.t > $outDir/tmp/t1.txt.t.f2
paste  $outDir/tmp/concheaderColumnar.txt $outDir/tmp/t1.txt.t.f2 > $outDir/tmp/t1.txt.t.swapped
matrixTranspose.py $outDir/tmp/t1.txt.t.swapped > $outDir/kmers.merged.$columnSelector.txt

columnSelector="experim_count"
cuta.py -f.kmer,@$columnSelector $outDir/kmers.merged.txt > $outDir/tmp/t1.txt
matrixTranspose.py $outDir/tmp/t1.txt > $outDir/tmp/t1.txt.t
cuta.py -f2-_1 $outDir/tmp/t1.txt.t > $outDir/tmp/t1.txt.t.f2
paste  $outDir/tmp/concheaderColumnar.txt $outDir/tmp/t1.txt.t.f2 > $outDir/tmp/t1.txt.t.swapped
matrixTranspose.py $outDir/tmp/t1.txt.t.swapped > $outDir/kmers.merged.$columnSelector.txt

rm $outDir/kmers.merged.txt
rm -R $outDir/tmp
