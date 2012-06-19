#!/bin/bash

if [ $# -ne 2 ]; then
	echo Usage: $0 RowMatrixFile outTinyRayPWMFile
	echo Description: Convert the row matrix file from constructSimplePWM.py to the PWM format used by http://demo.tinyray.com/weblogo
	exit 1
fi

RowMatrixFile=$1
outTinyRayPWMFile=$2

tfile1=`tempfile`

echo "A: 0.25" > $tfile1
echo "G: 0.25" >> $tfile1
echo "C: 0.25" >> $tfile1
echo "T: 0.25" >> $tfile1

tfile2=`tempfile`
tfile3=`tempfile`
cuta.py --ignore-truncated-rows -f1,3,2,4 $RowMatrixFile > $tfile2
matrixTranspose.py $tfile2 > $tfile3

paste $tfile1 $tfile3 > $outTinyRayPWMFile