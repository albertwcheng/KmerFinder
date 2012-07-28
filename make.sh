#!/bin/bash

if [ $# -lt 1 ];then 
	echo $0 CppUtilClassesPath
	exit 1
fi

CppUtilClassesPath=$1

g++ -O3 -I${CppUtilClassesPath} -o findKmers findKmers.cpp
g++ -O3 -I${CppUtilClassesPath} -o DKmerFinderCounter DKmerFinderCounter.cpp
g++ -O3 -I${CppUtilClassesPath} -o DKmerFinderMaster DKmerFinderMaster.cpp

g++ -O3 -o DNAPrefixTreeTest DNAPrefixTreeTest.cpp