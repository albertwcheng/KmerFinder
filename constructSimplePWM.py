#!/usr/bin/env python

from sys import *

def printUsageAndExit(programName):
	print >> stderr,"Usage:",programName,"kmerFile colKmer colScore seed"
	exit(1)


def toStrList(L):
	sL=[]
	for x in L:
		sL.append(str(x))
	return sL

def replaceSeedAtPosWith(seed,p,x):
	return seed[:p]+x+seed[p+1:]

def isSeq(s):
	for x in s:
		if x.upper() not in ['A','C','G','T','U']:
			return False
		
	return True

def normalizeTo1(L):
	#print >> stderr,L
	sumL=sum(L)
	P=[]
	for s in L:
		P.append(float(s)/sumL)
	return P

if __name__=='__main__':
	programName=argv[0]
	args=argv[1:]
	try:	
		kmerFile,colKmer,colScore,seed=args
	except:
		printUsageAndExit(programName)
	
	
	colKmer=int(colKmer)-1
	colScore=int(colScore)-1
	
	kmerScores=dict()
	
	
	
	fil=open(kmerFile)
	lino=0
	for lin in fil:
		if len(lin)<1 or lin[0]=="#":
			continue
	
		lino+=1
		if lino<2:
			continue
		
		fields=lin.rstrip("\r\n").split("\t")
		#print >> stderr,fields,colScore
		thisKmer=fields[colKmer]
		
		if not isSeq(thisKmer):
			continue
		
		thisScore=float(fields[colScore])
		kmerScores[thisKmer]=thisScore
		k=len(thisKmer)
		
		
	fil.close()
	
	#now construct PWM
	#output row matrix
	# f(A,1) f(C,1) f(G,1) f(T,1)
	# f(A,2) f(C,2) f(G,2) f(T,2)
	#  ::      ::      ::    ::
	# f(A,k) f(C,k) f(G,k) f(T,k)
	
	letters=['A','C','G','T']
	bgfreq=[0.25,0.25,0.25,0.25]
	
	#PWM=[]
	
	for row in range(0,k):
		rowout=[]
		#PWM.append(rowout)
		for letter in letters:
			neighbor=replaceSeedAtPosWith(seed,row,letter)
			#print >> stderr,neighbor
			try:
				rowout.append(kmerScores[neighbor])
			except:
				rowout.append(0)
		rowout=normalizeTo1(rowout)		
		print >> stdout,"\t".join(toStrList(rowout))
		

	
			