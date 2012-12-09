from sys import *
from getopt import getopt
from operator import itemgetter
from os.path import basename

class CountTree:
	def __init__(self):
		self.count=0
		self.D=dict()
	
	def addCount(self,x,keyPath):
	
		self.count+=x
		
		if len(keyPath)>0:
		
			try:
				DSub=self.D[keyPath[0]]
			except KeyError:
				DSub=CountTree()
				self.D[keyPath[0]]=DSub
			
			
			DSub.addCount(x,keyPath[1:])
	def increment(self,keyPath):
		self.addCount(1,keyPath)
		
	def addCount2(self,x,*args):
		self.addCount(x,args)
	
	def increment2(self,*args):
		self.addCount(1,args)
		
	def printTree(self,stream,label="root",prefix=""):
		print >> stream,"%s%s:%d" %(prefix,label,self.count)
		for k,v in self.D.items():
			v.printTree(stream,k,prefix+" ")
		
	def getCount(self,keyPath=[]):
		if len(keyPath)==0:
			return self.count
		
		cur=self
		for k in keyPath:
			try:
				cur=cur.D[k]
			except KeyError:
				return 0
		
		return cur.getCount()
	
	def getCount2(self,*args):
		return self.getCount(args)
	
def printUsageAndExit(programName):
	print >> stderr,"Usage:",basename(programName),"[wordListSeparatedByComma or +wordListFilename] fastqFileList outDir"
	print >> stderr,"--num-seq maxNumSeqToRead"
	print >> stderr,"--control-name name.default is the first file in the fastqFileList"
	print >> stderr,"[=====wordListFile======]"
	print >> stderr,"word1"
	print >> stderr,"word2"
	print >> stderr,"..."
	print >> stderr,"wordN"

	print >> stderr,"[=====fastqFileList======]"
	print >> stderr,"sampleName1<tab>fastqFilePath1"
	print >> stderr,"sampleNameN<tab>fastqFilePathN"
	exit(1)


COMPL={'A':'T',
	   'T':'A',
	   'C':'G',
	   'G':'C',
	   'U':'A',
	   'a':'t',
	   't':'a',
	   'u':'a',
	   'c':'g',
	   'g':'c'}

def complement(b):
	global COMPL
	try:
		return COMPL[b]
	except:
		return b

def reverseComplement(s):
	rc=""
	for b in s:
		rc=complement(b)+rc
		
	return rc

def String_findAll(s,sub,start=0,end=-1):
	if end!=-1:
		endSearch=end
	else:
		endSearch=len(s)
		
	curStart=start
	
	locations=[]
	
	while True:
		idx=s.find(sub,curStart,endSearch)
		if idx<0:
			break
			
		locations.append(idx)
		curStart=idx+1
	
	return locations	
	
def isPalindromic(s):
	return reverseComplement(s.upper())==s.upper()

def flipMatchRange(m,lenSeq):
	lenMatch=m[1]-m[0]
	newStart0=lenSeq-m[0]-lenMatch+1
	return (newStart0,newStart0+lenMatch)

def removeOverlappingMatches(matches):
	for i in range(len(matches)-1,0,-1): #from end to second
		if matches[i][0]<matches[i-1][1]: #overlapping
			del matches[i]

def findMatchesToWordListStruct(wordListStruct,seq):
	matches=[]
	#wordFreq=dict()
	#return [(start0,end1,strand(+,-,.)] sorted by start0
	rseq=reverseComplement(seq)
	lenSeq=len(seq)
	#forward:
	for word,palindromic in wordListStruct:
		lenMatch=len(word)
		founds=String_findAll(seq,word)
		freq=0
		
		if palindromic:
			strand="."
		else:
			strand="+"
		
		freq+=len(founds)
		
		for f in founds:
			matches.append((f,f+lenMatch,strand,word))

		
		if not palindromic:
			#now for the reversecomplement
			strand="-"
			founds=String_findAll(rseq,word)
			freq+=len(founds)
			for f in founds:
				matches.append(flipMatchRange((f,f+lenMatch),lenSeq)+(strand,word))
		
		#wordFreq[word]=freq
		
	matches.sort(key=itemgetter(0))
	#prevStart0=-1
	#for mat in matches:
	#	if mat[0]==prevStart0:
	#		#has duplicate start position!! no!!
	#		print >> stderr,"has duplicate start position"
	#		raise Exception()
	#	prevStart0=mat[0]		
	removeOverlappingMatches(matches)
	
	return matches #,wordFreq
	
def buildCompoundKmerStructs(wordListStruct,fastqFile,numSeqToRead=0): #wordListStruct=[(word,isPalindromic(word))], word in uppers
	spacingTree=CountTree() #3(spacing)->"++"(orientations)->("ACT","AGT")(words)
	copiesTree=CountTree() #4(numCopies)->"+-++"(orientations)->("ACT","ATT","AGT","ACT")(words)
	wordTree=CountTree() #ACT(word)->"+"(orientation)
	fil=open(fastqFile)
	lino=0
	Seqnum=0
	for lin in fil:
		lino+=1
		lin=lin.rstrip("\r\n")
		if lino%4==2: #line 4n+2 contains the sequence
			Seqnum+=1
			if numSeqToRead>0 and Seqnum>numSeqToRead:
				break
			if Seqnum%1000==1:
				print >> stderr,"reading seq",Seqnum
			seq=lin.upper() #ensure everything is uppercase
			
			matches=findMatchesToWordListStruct(wordListStruct,seq)
			#first add the wordTree
			orientationString=""
			wordTuple=[]
			for start0,end1,strand,word in matches:
				wordTree.increment2(word,strand)
				orientationString+=strand
				wordTuple.append(word)
			
			wordTuple=tuple(wordTuple)
			
			copiesTree.increment2(len(matches),orientationString,wordTuple)
			
			#now for each two
			for i in range(0,len(matches)-1):
				matchL=matches[i]
				matchR=matches[i+1]
				
				spacing=matchR[0]-matchL[1]
				orientationString=matchL[2]+matchR[2]
				words=(matchL[3],matchR[3])
				
				spacingTree.increment2(spacing,orientationString,words)
		
		
				
					
			
			
		
	fil.close()
	
	return (spacingTree,copiesTree,wordTree)


def testCountTree():
	tree=CountTree()
	tree.increment2(2,"+-")
	tree.increment2(3,"++",("ACT","GAT"))
	tree.increment2(3,"++",("ACT","GAT"),"leaf")
	tree.increment2(3,"++",("ACT","GAC"))
	tree.increment2(3,"++",("ACT","GAT"))
	tree.printTree(stderr)
	
	print >> stderr,tree.getCount2(2,"+-")
	print >> stderr,tree.getCount2(2)
	print >> stderr,tree.getCount2(3,"++")
	print >> stderr,tree.getCount2(3,"++",("ACT","GAT"))
	
	wordListStruct=[]
	wordListStruct.append(("ACTAGT",True))
	wordListStruct.append(("ACCAGC",False))
	m=findMatchesToWordListStruct(wordListStruct,"ACTAGTCCACCAGCGGTGCTGGT")
	print >> stderr,m
	#print >> stderr,wordFreq
	exit(0)
if __name__=='__main__':
	
	#testCountTree()
	
	
		
	programName=argv[0]
	opts,args=getopt(argv[1:],'',['num-seq=','control-name='])
	numSeq=0
	controlName=None
	for o,v in opts:
		if o=='--num-seq':
			numSeq=int(v)
		elif o=='--control-name':
			controlName=v
	
	
	try:
		wordList,fastqFiles,outDir=args
	except:
		printUsageAndExit(programName)
		
	
	if wordList[0]=='+':
		wordList=[]
		fil=open(wordList[1:])
		for lin in fil:
			lin=lin.rstrip("\r\n")
			wordList.append(lin)
		fil.close()
	else:
		wordList=wordList.split(",")
	
	
	#now process wordList to [(word,palindromic?)]
	wordListStruct=[]
	for word in wordList:
		wordListStruct.append((word.upper(),isPalindromic(word)))
	
	spacingTree,copiesTree,wordTree=buildCompoundKmerStructs(wordListStruct,fastqFiles,numSeq)
	print >> stderr,"spacingTree:"
	spacingTree.printTree(stderr)
	print >> stderr,"-----------"
	print >> stderr,"copiesTree:"
	copiesTree.printTree(stderr)
	print >> stderr,"-----------"	
	print >> stderr,"wordTree:"
	wordTree.printTree(stderr)	
	print >> stderr,"-----------"