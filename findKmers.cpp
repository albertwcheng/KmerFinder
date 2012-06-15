/***************************************************************************
 Copyright 2012 Wu Albert Cheng <albertwcheng@gmail.com>
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 *******************************************************************************/ 

#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <inttypes.h>
#include <stdlib.h>
#include <time.h>
using namespace std;


#define REQUIRE_BACKGROUND_PRESENCE
//#define USE_MULTISET

#ifdef USE_MULTISET
#define setType multiset
#else
#define setType set
#endif 

#define uint64 uint64_t

#define VERSION "1.0"

template <class T>
class SmartPtr
{
private:
	T* ptr;
public:
	inline SmartPtr(T* _ptr=NULL):ptr(_ptr){}
	
	inline void operator=( T *_ptr)
	{
		ptr=_ptr;
	}
	
	inline operator T*() const
	{
		return ptr;
	}
	
	inline T& operator *() const
	{
		return *ptr;
	}
	
	inline T* operator ->() const
	{
		return ptr;
	}
	
	inline bool operator<(const SmartPtr<T>& right) const
	{
		return (**this)<(*right);
	}
	inline bool operator>(const SmartPtr<T>& right) const
	{
		return (**this)>(*right);
	}
	inline bool operator==(const SmartPtr<T>& right) const
	{
		return (**this)==(*right);
	}
	inline bool operator<=(const SmartPtr<T>& right) const
	{
		return !(**this>*right);
	}
	
	inline bool operator>=(const SmartPtr<T>& right) const
	{
		return !(**this<*right);
	}
	
	inline bool operator!=(const SmartPtr<T>& right) const
	{
		return !((**this)==(*right));
	}
	
	inline bool ptrEquals(const SmartPtr<T>& right) const
	{
		return ptr==right.ptr;
	}
	
};



class SimpleFastqReader
{
	
public:
	virtual bool onNewSeq(const string& name,const string& seq,const string& quals)=0;
	
	inline void readFile(const string& filename){
		ifstream fin(filename.c_str());
		string name;
		string seq;
		string quals;

		uint64 reallino=0;
		uint64 lino=0;
		
		bool continuetoread=true;
		
		while(continuetoread && fin.good()){
			string line;
			getline(fin,line); 
			
			lino++;
			
			
			if(line.length()==0)
				continue;
			
			reallino++;
			
			switch(reallino%4){
				case 1: //a start of new seq
					
					if(line[0]!='@'){
						cerr<<"Warning: improper Fastq formatting at "<<filename<<":line "<<lino<<"="<<line<<".Expecting @"<<endl;
					}
					
					if(seq.length()>0){
						continuetoread=this->onNewSeq(name,seq,quals);
					}
					
					
					name=line.substr(1); //to the end
					seq="";
					quals="";
					

					
					
					break;
				case 2:
					seq=line;
					break;
				case 3:
					if(line[0]!='+'){
						cerr<<"Warning: improper Fastq formatting at "<<filename<<":line "<<lino<<"="<<line<<".Expecting +"<<endl;
					}
					break;
				case 0:
					quals=line;
					break;
					
				default:
					break;
					
			}
		}
		
		if(seq.length()>0){
			this->onNewSeq(name,seq,quals);
		}
		
		fin.close();
	}
	
};


class SeqRecord;

class kmerFinder;

class kmerRecord
{
public:

	string kmerSeq;
	setType<SeqRecord*> seqs; //this is the pointer to the seqs that have this kmer 
	setType<SeqRecord*> bgseqs;
	
	uint64 unsubtractedBgCount;
	uint64 unsubtractedFgCount;
	
	kmerRecord(const string& _kmerSeq):kmerSeq(_kmerSeq),unsubtractedBgCount(0),unsubtractedFgCount(0){
		
	}
	
	inline uint64 fgInstances() const{
		return this->seqs.size();
	}
	
	inline uint64 bgInstances() const{
		return this->bgseqs.size();	
	}
	
	inline uint64 unsubtractedFgInstances() const{
		return unsubtractedFgCount;
	}
	
	inline uint64 unsubtractedBgInstances() const{
		return unsubtractedBgCount;	
	}
	
	inline double unsubtractedEnrichment() const{
#ifdef REQUIRE_BACKGROUND_PRESENCE
		return double(this->unsubtractedFgInstances())/unsubtractedBgInstances();
#else		
		return double(this->unsubtractedFgInstances())/(unsubtractedBgInstances()+1); //background plus 1 to avoid division by zero
		//no need to care about the totalkmers in foreground and background becoz they are the same for all kmers. argmax not important to know these numbers
#endif
	}
	
	inline double unsubtractedNormalizedEnrichment(double totalFgKmers,double totalBgKmers) const{
#ifdef REQUIRE_BACKGROUND_PRESENCE
		return (double(this->unsubtractedFgInstances())/totalFgKmers)/(double(unsubtractedBgInstances())/totalBgKmers);
#else
		return (double(this->unsubtractedFgInstances())/totalFgKmers)/(double(unsubtractedBgInstances()+1)/totalBgKmers);
#endif
		
		
	}	
	inline double enrichment() const{
		
	#ifdef REQUIRE_BACKGROUND_PRESENCE
		return double(this->seqs.size())/this->bgseqs.size();
	#else		
		return double(this->seqs.size())/(this->bgseqs.size()+1); //background plus 1 to avoid division by zero
		//no need to care about the totalkmers in foreground and background becoz they are the same for all kmers. argmax not important to know these numbers
	#endif
		
	}
	
	inline double normalizedEnrichment(double totalFgKmers,double totalBgKmers) const{
	#ifdef REQUIRE_BACKGROUND_PRESENCE
		return (double(this->seqs.size())/totalFgKmers)/(double(this->bgseqs.size())/totalBgKmers);
	#else
		return (double(this->seqs.size())/totalFgKmers)/(double(this->bgseqs.size()+1)/totalBgKmers);
	#endif
		
	}
	
	inline double normalizedEnrichmentByNumSeqs(double totalFgSeq,double totalBgSeq) const{
	#ifdef REQUIRE_BACKGROUND_PRESENCE
		return (double(this->seqs.size())/totalFgSeq)/(double(this->bgseqs.size())/totalBgSeq);
	#else
		return (double(this->seqs.size())/totalFgSeq)/(double(this->bgseqs.size()+1)/totalBgSeq);
	#endif		
	}
	
	inline bool operator < (const kmerRecord& right) const{
		return this->enrichment()<right.enrichment();
	
	}
	
	void kmerDisappearFromSeqs();
	
	void removeSequencesContainingThisKmer();
};

class SeqRecord
{
public:
	set<kmerRecord*> kmersFromThisSeq;
	inline void disappearFromKmersExcept(kmerRecord* except){
		for(set<kmerRecord*>::iterator i=kmersFromThisSeq.begin();i!=kmersFromThisSeq.end();i++){
			kmerRecord* thisKmer=(*i);
			if(thisKmer!=except){
				(*i)->seqs.erase(this);
			}
		}
	}
	inline void bgdisappearFromKmersExcept(kmerRecord* except){
		for(set<kmerRecord*>::iterator i=kmersFromThisSeq.begin();i!=kmersFromThisSeq.end();i++){
			kmerRecord* thisKmer=(*i);
			if(thisKmer!=except){
				(*i)->bgseqs.erase(this);
			}
		}
	}	
	void registerKmer(kmerRecord* _kmer){
		_kmer->seqs.insert(this);
		kmersFromThisSeq.insert(_kmer);
	}
	
	void registerBgKmer(kmerRecord* _kmer){
		_kmer->bgseqs.insert(this);
		kmersFromThisSeq.insert(_kmer);	
	}
	
	
	void unregisterKmer(kmerRecord* _kmer){
		kmersFromThisSeq.erase(_kmer);
	}
};


void kmerRecord::removeSequencesContainingThisKmer(){
	for(setType<SeqRecord*>::iterator i=seqs.begin();i!=seqs.end();i++){
		SeqRecord* thisSeq=(*i);
		thisSeq->disappearFromKmersExcept(this);
	}
	
	for(setType<SeqRecord*>::iterator i=bgseqs.begin();i!=bgseqs.end();i++){
		SeqRecord* thisSeq=(*i);
		thisSeq->bgdisappearFromKmersExcept(this);	
	}
	
}

void kmerRecord::kmerDisappearFromSeqs(){
	for(setType<SeqRecord*>::iterator i=seqs.begin();i!=seqs.end();i++){
		SeqRecord* thisSeq=(*i);
		thisSeq->unregisterKmer(this);
	}
	
	for(setType<SeqRecord*>::iterator i=bgseqs.begin();i!=bgseqs.end();i++){
		SeqRecord* thisSeq=(*i);
		thisSeq->unregisterKmer(this);
	}
}

#define READ_FOREGROUND 1
#define READ_BACKGROUND 2




class kmerFinder:public SimpleFastqReader
{
public:
	list<SmartPtr<kmerRecord> > kmers;
	vector<kmerRecord*> bgOnlyKmers; //for keeping track the memory for those background only kmer (for destructor)
	
	vector<SeqRecord*> seqs;
	vector<SeqRecord*> bgseqs;
	
	int mode;
	string foregroundFilename;
	string backgroundFilename;
	int k; //wordsize
	
	uint64 fgseqtoread;
	uint64 bgseqtoread;
	
	uint64 foregroundTotalKmerCount;
	uint64 backgroundTotalKmerCount;
	
	uint64 numSeqForeground;
	uint64 numSeqBackground;
	
	map<string,SmartPtr<kmerRecord> > kmerInitStruct; //for storing kmers while reading file initially
	
	
	kmerFinder(int _k,uint64 _fgseqtoread,uint64 _bgseqtoread):mode(READ_BACKGROUND),k(_k),foregroundTotalKmerCount(0),backgroundTotalKmerCount(0),numSeqForeground(0),numSeqBackground(0),fgseqtoread(_fgseqtoread),bgseqtoread(_bgseqtoread){
		
	}
	
	inline uint64 numUniqKmers(){
		return kmers.size();
	}
	
	bool isValidSeq(const string& seq){
		for(int i=0;i<seq.length();i++){
			switch(seq[i]){
				case 'A':
				case 'T':
				case 'C':
				case 'G':
				case 'a':
				case 't':
				case 'c':
				case 'g':
				case 'U':
				case 'u':
					break;
				default:
					return false;
			}
			
		}
		return true;
	}
	
	inline void readForeground(const string&filename){
		mode=READ_FOREGROUND;
		foregroundFilename=filename;
		this->readFile(filename);
		
		//now put the kmer into list
		for(map<string,SmartPtr<kmerRecord> >::iterator i=kmerInitStruct.begin();i!=kmerInitStruct.end();i++)
		{
			kmers.push_back(i->second);
		}
		

	}
	
	inline void readBackground(const string&filename){
		mode=READ_BACKGROUND;
		backgroundFilename=filename;
		this->readFile(filename);
		
		
		
		#ifdef REQUIRE_BACKGROUND_PRESENCE
		
		cerr<<"removing backgroundcount=0 kmers..."<<endl;
		
		vector<list<SmartPtr<kmerRecord> >::iterator> toRemove;
		
		for(list<SmartPtr<kmerRecord> >::iterator i=kmers.begin();i!=kmers.end();i++){
			if((*i)->bgInstances()<1){
				toRemove.push_back(i);
			}else if((*i)->fgInstances()<1){
				bgOnlyKmers.push_back(*i);
			}
		}
		
		for(vector<list<SmartPtr<kmerRecord> >::iterator>::iterator i=toRemove.begin();i!=toRemove.end();i++){
			list<SmartPtr<kmerRecord> >::iterator removeeI=*i;
			kmerRecord* removee=*removeeI;
			removee->kmerDisappearFromSeqs();
			
			kmers.erase(removeeI);
			delete removee;//now can remove this from memory
		}
		
		cerr<<toRemove.size()<<" unique backgroundcount=0 kmers removed"<<endl;
		
		#endif
		
		//now we can remove the map kmerInitStruct to save space
		kmerInitStruct.clear();	
		
		//now calculate foregroundTotalKmerCount and backgroundTotalKmerCount
		for(list<SmartPtr<kmerRecord> >::iterator i=kmers.begin();i!=kmers.end();i++){
			
			(*i)->unsubtractedFgCount=(*i)->fgInstances();
			foregroundTotalKmerCount+=(*i)->fgInstances();
			(*i)->unsubtractedBgCount=(*i)->bgInstances();
			backgroundTotalKmerCount+=(*i)->bgInstances();
		}

			
	}
	
	kmerRecord* getNextEnrichedKmers(bool printCPUTimeUsed=false){
		if(kmers.size()==0)
			return NULL;
		
		clock_t t1=clock();
		kmers.sort(); //sort the linked list by enrichment score
		clock_t t2=clock();
		kmerRecord*result=kmers.back();
		result->removeSequencesContainingThisKmer();
		clock_t t3=clock();
		kmers.pop_back();
		
		if(printCPUTimeUsed){
			cerr<<"\tfound kmer "<<result->kmerSeq<<" with occurences:"<<result->fgInstances()<<"; time invested:"<<((t2-t1)/(double)CLOCKS_PER_SEC)<<"s(Sorting) "<<((t3-t2)/(double)CLOCKS_PER_SEC)<<"s(Removing sequences) "<<((t3-t1)/(double)CLOCKS_PER_SEC)<<"s(Total)"<<endl;
		}
		return result;
	}
	
	bool onNewSeq(const string&name,const string&seq,const string&quals){
		if(mode==READ_FOREGROUND){
			
			numSeqForeground++;
			
			if(numSeqForeground%100000==1){
				cerr<<"\treading seq\t"<<numSeqForeground<<endl;
			}
			
			SeqRecord *thisSeqRecord=new SeqRecord;
			
			seqs.push_back(thisSeqRecord);
			
			//for foreground fill the kmers;
			for(int i=0;i<=seq.length()-k;i++){
				string thisKmerSeq(seq.substr(i,k));
				
				if(!isValidSeq(thisKmerSeq)){
					continue;
				}
				
				//this->foregroundTotalKmerCount++;
				
				//now insert this kmer
				SmartPtr<kmerRecord> thisKmer;
				
				map<string,SmartPtr<kmerRecord> >::iterator kmerI=kmerInitStruct.find(thisKmerSeq);
				if(kmerI==kmerInitStruct.end()){
					thisKmer=new kmerRecord(thisKmerSeq);
					kmerInitStruct.insert(map<string,SmartPtr<kmerRecord> >::value_type(thisKmerSeq,thisKmer));
				}else {
					thisKmer=kmerI->second;
				}
				
				//register thisSeqRecord in thisKmer
				
				thisSeqRecord->registerKmer(thisKmer);
				
			}
			
			return (fgseqtoread==0 || numSeqForeground<fgseqtoread);
			
			
		}else {
			
			numSeqBackground++;
			
			SeqRecord *thisSeqRecord=new SeqRecord;
			
			bgseqs.push_back(thisSeqRecord);
			
			if(numSeqBackground%100000==1){
				cerr<<"\treading seq\t"<<numSeqBackground<<endl;
			}
			
			for(int i=0;i<=seq.length()-k;i++){
				string thisKmerSeq(seq.substr(i,k));
				
				if(!isValidSeq(thisKmerSeq)){
					continue;
				}
				
				//this->backgroundTotalKmerCount++;
				
				SmartPtr<kmerRecord> thisKmer;
				map<string,SmartPtr<kmerRecord> >::iterator kmerI=kmerInitStruct.find(thisKmerSeq);
				if(kmerI==kmerInitStruct.end()){
					thisKmer=new kmerRecord(thisKmerSeq);
					kmerInitStruct.insert(map<string,SmartPtr<kmerRecord> >::value_type(thisKmerSeq,thisKmer));
				}else {
					thisKmer=kmerI->second;
				}
				
				thisSeqRecord->registerBgKmer(thisKmer);
				
				
			}		
			
			
			return (bgseqtoread==0 || numSeqBackground<bgseqtoread);
			
		
		}

	}
	
	~kmerFinder(){
		for(list<SmartPtr<kmerRecord> >::iterator i=kmers.begin();i!=kmers.end();i++){
			delete (*i);
		}
		for(vector<SeqRecord*>::iterator i=seqs.begin();i!=seqs.end();i++){
			delete (*i);
		}
		for(vector<SeqRecord*>::iterator i=bgseqs.begin();i!=bgseqs.end();i++){
			delete (*i);
		}
		for(vector<kmerRecord*>::iterator i=bgOnlyKmers.begin();i!=bgOnlyKmers.end();i++){
			delete (*i);	
		}
	}
	
	void printStat(ostream& os){
		os<<"#ForegroundSeqCount="<<numSeqForeground<<endl;
		os<<"#ForegroundKmerCount="<<foregroundTotalKmerCount<<endl;
		
		os<<"#BackgroundSeqCount="<<numSeqBackground<<endl;
		os<<"#BackgroundKmerCount="<<backgroundTotalKmerCount<<endl;
		
		os<<"#UniqKmerCount="<<kmers.size()<<endl;
	}
	
};

class TestSimpleFastqReader:public SimpleFastqReader{
public:
	int seqcount;
	TestSimpleFastqReader():seqcount(0){}
	bool onNewSeq(const string& name,const string& seq,const string& quals)
	{
		seqcount++;
		cerr<<seqcount<<":got new seq! name="<<name<<" seq="<<seq<<" qual="<<quals<<endl;
		
		return seqcount<10; //read only ten sequences
	}
};

int main(int argc,char**argv){
	
	
	/*TestSimpleFastqReader f;
	f.readFile(argv[1]);
	return 0;*/
	
	cerr<<"Find Kmers "<<VERSION<<endl;
	cerr<<"(";
#ifdef REQUIRE_BACKGROUND_PRESENCE
	cerr<<"REQUIRE_BACKGROUND_PRESENCE=ON)"<<endl;
#else
	cerr<<"REQUIRE_BACKGROUND_PRESENCE=OFF)"<<endl;
#endif
	
	if(argc<7){ 
		cerr<<"Usage:" <<argv[0]<<" <fgfilename> <howmanyfgseqtoread> <bgfilename> <howmanybgseqtoread> <k> <howmanyToFind>"<<endl;
		cerr<<"Description: find <howmanyToFind> top enriched <k>-mers from HT-SELEX experiments using <fgfilename> and <bgfilename> fastq files assuming there's only one motif per sequence"<<endl;
		cerr<<"Specify <howmanyfgseqtoread> or <howmanybgseqtoread>=0 to read all"<<endl;
		cerr<<"Specify <howmanyToFind>=0 to print all kmers"<<endl;
		return 1;
	}
	
	
	string fgfilename=argv[1];
	int fgseqtoread=atoi(argv[2]);
	string bgfilename=argv[3];
	int bgseqtoread=atoi(argv[4]);
	int k=atoi(argv[5]);
	int howmanyToFind=atoi(argv[6]);
	
	kmerFinder theFinder(k,fgseqtoread,bgseqtoread);
	
	cerr<<"command:"<<argv[0];
	for(int i=1;i<argc;i++){
		cerr<<" "<<argv[i];
	}
	
	cerr<<endl;
	
	cerr<<"reading foreground..."<<fgfilename<<endl;
	clock_t t1=clock();
	theFinder.readForeground(fgfilename);
	clock_t t2=clock();
	cerr<<"done reading and processing foreground CPU time used="<<((t2-t1)/(double)CLOCKS_PER_SEC)<<" seconds"<<endl;
	cerr<<"reading background..."<<bgfilename<<endl;
	theFinder.readBackground(bgfilename);
	clock_t t3=clock();
	cerr<<"done reading and processing background CPU time used="<<((t3-t2)/(double)CLOCKS_PER_SEC)<<" seconds"<<endl;
	
	cout<<"#command=\""<<argv[0];
	for(int i=1;i<argc;i++){
		cout<<" "<<argv[i];
	}
	
	cout<<"\""<<endl;
	cout<<"#foreground="<<fgfilename<<endl;
	cout<<"#background="<<bgfilename<<endl;
	cout<<"#k="<<k<<endl;
	cout<<"#howmanyToFind="<<howmanyToFind<<endl;

	
	theFinder.printStat(cout);
	cerr<<"finding kmers..."<<endl;
	
	cout<<"kmer\tenrichment\tcontrol_count\texperim_count"<<endl;
	//cout<<"kmer\tenrichment\tnormalizedEnrichment\tfgInstances\tbgInstances\tunsubtractedFgInstances\tunsubtractedBgInstances\tunsubtractedEnrichment\tunsubtractedNormalizedEnrichment"<<endl;
	
	if(howmanyToFind<1){
		howmanyToFind=theFinder.numUniqKmers();
	}
	
	for(int i=0;i<howmanyToFind;i++){
		kmerRecord* nextKmer=theFinder.getNextEnrichedKmers(true);
		if(!nextKmer)
			break;
		//cout<<nextKmer->kmerSeq<<"\t"<<nextKmer->enrichment()<<"\t"<<nextKmer->normalizedEnrichment(theFinder.foregroundTotalKmerCount,theFinder.backgroundTotalKmerCount)<<"\t"<<nextKmer->fgInstances()<<"\t"<<nextKmer->bgInstances()<<"\t"<<nextKmer->unsubtractedFgInstances()<<"\t"<<nextKmer->unsubtractedBgInstances()<<"\t"<<nextKmer->unsubtractedEnrichment()<<"\t"<<nextKmer->unsubtractedNormalizedEnrichment(theFinder.foregroundTotalKmerCount,theFinder.backgroundTotalKmerCount)<<endl;
		cout<<nextKmer->kmerSeq<<"\t"<<nextKmer->enrichment()<<"\t"<<nextKmer->bgInstances()<<"\t"<<nextKmer->fgInstances()<<endl;
	}
	clock_t t4=clock();
	cout<<"#CPUTimeReadForeground="<<((t2-t1)/(double)CLOCKS_PER_SEC)<<"s"<<endl;
	cout<<"#CPUTimeReadBackground="<<((t3-t2)/(double)CLOCKS_PER_SEC)<<"s"<<endl;
	cout<<"#CPUTimeFindKmers="<<((t4-t3)/(double)CLOCKS_PER_SEC)<<"s"<<endl;
	cout<<"#CPUTimeTotal="<<((t4-t1)/(double)CLOCKS_PER_SEC)<<"s"<<endl;
	
	
	cerr<<"done finding kmers CPU time used="<<((t4-t3)/(double)CLOCKS_PER_SEC)<<" seconds"<<endl;
	cerr<<"done total CPU time used="<<((t4-t1)/(double)CLOCKS_PER_SEC)<<" seconds"<<endl;
	return 0;
	
}
