from sys import *
from getopt import getopt
from operator import itemgetter
from os.path import basename,exists
from os import makedirs
import os.path
import pickle
import fisher #http://pypi.python.org/pypi/fisher/

class CountTree:
	def __init__(self):
		self.count=0
		self.D=dict()
	
	def items():
		return self.D.items()
	
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
		
	def getCount2(self,*args):
		self.getCount(args)
	
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
	print >> stderr,"Usage 1:",basename(programName),"compileTree [wordListSeparatedByComma or +wordListFilename] fastqFileList outDir"
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
	
	print >> stderr,""
	print >> stderr,"Usage 2:",basename(programName),"printTree outDir"
	
	print >> stderr,""
	print >> stderr,"Usage 3:",basename(programName),"test1 outDir"	
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

def findMatchesToWordListStruct(wordListStruct,seq,onlyFoward):
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

		
		if not palindromic and not onlyFoward:
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
	
def buildCompoundKmerStructs(wordListStruct,fastqFile,onlyFoward,numSeqToRead=0): #wordListStruct=[(word,isPalindromic(word))], word in uppers
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
			if Seqnum%100000==1:
				print >> stderr,fastqFile,"seq",Seqnum
			seq=lin.upper() #ensure everything is uppercase
			
			matches=findMatchesToWordListStruct(wordListStruct,seq,onlyFoward)
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
	
class TreeCollection:
	def __init__(self,wordListStruct,fastqFiles,onlyFoward,numSeq):
		self.spacingTree,self.copiesTree,self.wordTree=buildCompoundKmerStructs(wordListStruct,fastqFiles,onlyFoward,numSeq)
	

class KmerSpacingAnalysisSession:
	def __init__(self,wordList,fastqFileList,outDir,sampleForest,onlyForward,controlName):
		self.wordList=wordList
		self.fastqFileList=fastqFileList
		self.outDir=outDir
		self.sampleForest=sampleForest
		self.onlyForward=onlyForward
		self.controlName=controlName

def printTree(programName,inargs):
	opts,args=getopt(inargs,'',[])
	
	try:
		outdir,=args
	except:
		printUsageAndExit(programName)
	
	pickIn=open(outdir+os.path.sep+"session.pickle")
	session=pickle.load(pickIn)
	pickIn.close()
	print >> stdout,"wordList:",session.wordList
	print >> stdout,"fastqFileList:",session.fastqFileList
	print >> stdout,"outDir:",session.outDir
	print >> stdout,"onlyForward:",session.onlyForward
	print >> stdout,"controlName:",session.controlName
	
	for sampleName,thisSampleTrees in session.sampleForest.items():	
		print >> stdout,"sample",sampleName
		print >> stdout,"spacingTree:"
		thisSampleTrees.spacingTree.printTree(stderr)
		print >> stdout,"-----------"
		print >> stdout,"copiesTree:"
		thisSampleTrees.copiesTree.printTree(stderr)
		print >> stdout,"-----------"	
		print >> stdout,"wordTree:"
		thisSampleTrees.wordTree.printTree(stderr)	
		print >> stdout,"-----------"

class PrintableBoxes:
	def __init__(self,sep="\t"):
		self._boxes=dict()
		self.maxX=-1
		self.maxY=-1
		self.sep=sep
	
	def put(self,x,y,item):
		try:
			xarr=self._boxes[y]
		except:
			xarr=dict()
			self._boxes[y]=xarr
		
		xarr[x]=item
		self.maxX=max(self.maxX,x)
		self.maxY=max(self.maxY,y)
	
	def printToStream(self,stream):
		for y in range(0,self.maxY+1):
			for x in range(0,self.maxX+1):
				if x>0:
					stream.write(self.sep)
				
				try:
					stream.write(str(self._boxes[y][x]))
				except:
					pass
					
			print >> stream,""
	

def test1(programName,inargs):
	opts,args=getopt(inargs,'',['from=','to='])
	
	_from=0
	_to=10
	
	for o,v in opts:
		if o=='--from':
			_from=int(v)
		elif o=='--to':
			_to=int(v)
	
	try:
		outdir,=args
	except:
		printUsageAndExit(programName)
	
	pickIn=open(outdir+os.path.sep+"session.pickle")
	session=pickle.load(pickIn)
	pickIn.close()

	boxes=PrintableBoxes()

	controlTrees=session.sampleForest[session.controlName]
	
	curX=0
	curY=0
	
	for spacing in range(_from,_to+1):
		boxes.put(curX+1,curY,"spacing")
		boxes.put(curX+2,curY,spacing)	
		
		curX+=6
		
	curX=0	
	curY+=1
	
	
	
	for sampleName,thisSampleTrees in session.sampleForest.items():	
		if sampleName==session.controlName:
			continue #skip
			
		boxes.put(curX,curY+1,"sample")
		boxes.put(curX,curY+2,sampleName)
				
		curX+=1

		for spacing in range(_from,_to+1):
			allSpacingInBg=controlTrees.spacingTree.getCount()
			thisSpacingInBg=controlTrees.spacingTree.getCount2(spacing)
			allSpacingInFg=thisSampleTrees.spacingTree.getCount()
			thisSpacingInFg=thisSampleTrees.spacingTree.getCount2(spacing)
			notSpacingInBg=allSpacingInBg-thisSpacingInBg
			notSpacingInFg=allSpacingInFg-thisSpacingInFg
			#now construct 2x2
			
			#now fisher
			fisher_pvalues=fisher.pvalue(thisSpacingInFg,notSpacingInFg,thisSpacingInBg,notSpacingInBg)
			boxes.put(curX+1,curY,"number with this spacing")
			boxes.put(curX+2,curY,"number with other spacing")
			boxes.put(curX+3,curY,"p-values")
			boxes.put(curX,curY+1,"test")
			boxes.put(curX,curY+2,"control")
			boxes.put(curX+1,curY+1,thisSpacingInFg)
			boxes.put(curX+2,curY+1,notSpacingInFg)
			boxes.put(curX+1,curY+2,thisSpacingInBg)
			boxes.put(curX+2,curY+2,notSpacingInBg)
			boxes.put(curX+3,curY+1,fisher_pvalues.right_tail)
			boxes.put(curX+3,curY+2,fisher_pvalues.left_tail)
			
			curX+=6
		
		curX=0
		curY+=6
		
	fout=open(outdir+os.path.sep+"test1.xls","w")
	boxes.printToStream(fout)
	fout.close()

def compileTree(programName,inargs):
	
	opts,args=getopt(inargs,'',['num-seq=','control-name=','only-forward'])
	numSeq=0
	controlName=None
	onlyForward=False
	for o,v in opts:
		if o=='--num-seq':
			numSeq=int(v)
		elif o=='--control-name':
			controlName=v
		elif o=='--only-forward':
			onlyForward=True
	
	try:
		wordList,fastqFileList,outDir=args
	except:
		printUsageAndExit(programName)
	
	makedirs(outDir)
		
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
	
	print >> stderr,args,wordListStruct,wordList
	
	sampleForest=dict()
	controlName=""
	
	fil=open(fastqFileList)
	for lin in fil:
		
		lin=lin.rstrip("\r\n")
		sampleName,filename=lin.split("\t")
		if controlName=="": #if controlName hasn't been specified in optional args, first one is assumed to be control
			controlName=sampleName
			
		thisSampleTrees=TreeCollection(wordListStruct,filename,onlyForward,numSeq)
		sampleForest[sampleName]=thisSampleTrees
		print >> stderr,"spacingTree:"
		thisSampleTrees.spacingTree.printTree(stderr)
		print >> stderr,"-----------"
		print >> stderr,"copiesTree:"
		thisSampleTrees.copiesTree.printTree(stderr)
		print >> stderr,"-----------"	
		print >> stderr,"wordTree:"
		thisSampleTrees.wordTree.printTree(stderr)	
		print >> stderr,"-----------"
	fil.close()	
	
	
	controlSample=sampleForest[controlName]
	
	#now try test1
	session=KmerSpacingAnalysisSession(wordListStruct,fastqFileList,outDir,sampleForest,onlyForward,controlName)
	pickOut=open(outDir+os.path.sep+"session.pickle","w")
	pickle.dump(session,pickOut)
	pickOut.close()


funcTables={ "compileTree" : compileTree, "printTree" : printTree, "test1": test1}

def testPrintableBoxes():
	boxes=PrintableBoxes()
	boxes.put(5,3,"Albert")
	boxes.put(2,2,"Amy")
	boxes.put(0,7,"Jason")
	boxes.put(9,9,12)
	boxes.printToStream(stdout)
	exit(1)

if __name__=='__main__':
	
	#testCountTree()
	#testPrintableBoxes()
	
	programName=argv[0]
	
	try:
		commandName=argv[1]
		if commandName not in funcTables:
			print >> stderr,commandName,"command not defined. Abort"
			raise KeyError
		
		
	except:
		printUsageAndExit(programName)	
		

	funcTables[commandName](programName,argv[2:])