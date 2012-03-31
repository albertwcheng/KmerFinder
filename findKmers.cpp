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



class SimpleFastxReader
{
	
public:
	virtual void onNewSeq(const string& name,const string& seq,const string& quals)=0;
	
	inline void readFile(const string& filename){
		ifstream fin(filename.c_str());
		string name;
		string seq;
		string quals;
		string *seqFeedDest=NULL; //point to either seq or quals for inserting into either one of them
		char seqDirective='>';
		bool firstSeq=true;
		
		while(fin.good()){
			string line;
			getline(fin,line); 
			//cerr<<line<<endl;
			if(line.length()==0)
				continue;
			switch(line[0]){
				case '@': //a start of new seq
				case '>':
					
					if(firstSeq){
						seqDirective=line[0];
						firstSeq=false;
					}else if(line[0]==seqDirective){
						if(seq.length()>0){
							this->onNewSeq(name,seq,quals);
						}
						name=line.substr(1); //to the end
						seq="";
						quals="";
						seqFeedDest=&seq;
					}else{
						//this is potentially ">" and is within quals
						(*seqFeedDest)+=line;
					}
					
					
					break;
				
				case '+':
					quals="";
					seqFeedDest=&quals;
					break;
				
				default:
					if(seqFeedDest){
						(*seqFeedDest)+=line;
					}
					break;
				
			}
		}
		
		if(seq.length()>0){
			this->onNewSeq(name,seq,quals);
		}
		
		fin.close();
	}

};



class SimpleFastqReader
{
	
public:
	virtual void onNewSeq(const string& name,const string& seq,const string& quals)=0;
	
	inline void readFile(const string& filename){
		ifstream fin(filename.c_str());
		string name;
		string seq;
		string quals;

		uint64 reallino=0;
		uint64 lino=0;
		
		while(fin.good()){
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
						this->onNewSeq(name,seq,quals);
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
	//kmerFinder* parent;
	string kmerSeq;
	multiset<SeqRecord*> seqs; //this is the pointer to the seqs that have this kmer
	uint64 backgroundCount; //this will store the number of times bacskground has this kmer
	
	kmerRecord(const string& _kmerSeq):kmerSeq(_kmerSeq),backgroundCount(0){
		
	}
	
	inline uint64 fgInstances() const{
		return this->seqs.size();
	}
	
	
	inline double enrichment() const{
		
	#ifdef REQUIRE_BACKGROUND_PRESENCE
		return double(this->seqs.size())/backgroundCount;
	#else		
		return double(this->seqs.size())/(backgroundCount+1); //background plus 1 to avoid division by zero
		//no need to care about the totalkmers in foreground and background becoz they are the same for all kmers. argmax not important to know these numbers
	#endif
		
	}
	
	inline double normalizedEnrichment(double totalFgKmers,double totalBgKmers) const{
	#ifdef REQUIRE_BACKGROUND_PRESENCE
		return (double(this->seqs.size())/totalFgKmers)/(double(backgroundCount)/totalBgKmers);
	#else
		return (double(this->seqs.size())/totalFgKmers)/(double(backgroundCount+1)/totalBgKmers);
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
	inline void disappearFromKmers(){
		for(set<kmerRecord*>::iterator i=kmersFromThisSeq.begin();i!=kmersFromThisSeq.end();i++){
			(*i)->seqs.erase(this);
		}
	}
	
	void registerKmer(kmerRecord* _kmer){
		_kmer->seqs.insert(this);
		kmersFromThisSeq.insert(_kmer);
	}
	
	void unregisterKmer(kmerRecord* _kmer){
		kmersFromThisSeq.erase(_kmer);
	}
};


void kmerRecord::removeSequencesContainingThisKmer(){
	for(multiset<SeqRecord*>::iterator i=seqs.begin();i!=seqs.end();i++){
		SeqRecord* thisSeq=(*i);
		thisSeq->disappearFromKmers();
	}
	
}

void kmerRecord::kmerDisappearFromSeqs(){
	for(multiset<SeqRecord*>::iterator i=seqs.begin();i!=seqs.end();i++){
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
	vector<SeqRecord*> seqs;
	int mode;
	string foregroundFilename;
	string backgroundFilename;
	int k; //wordsize
	
	uint64 foregroundTotalKmerCount;
	uint64 backgroundTotalKmerCount;
	
	
	uint64 numSeqForeground;
	uint64 numSeqBackground;
	
	map<string,SmartPtr<kmerRecord> > kmerInitStruct; //for storing kmers while reading file initially
	
	
	kmerFinder(int _k):mode(READ_BACKGROUND),k(_k),foregroundTotalKmerCount(0),backgroundTotalKmerCount(0),numSeqForeground(0),numSeqBackground(0){
		
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
			if((*i)->backgroundCount<1){
				toRemove.push_back(i);
			//}else{
				//cerr<<"bigger"<<endl;
			}
		}
		
		for(vector<list<SmartPtr<kmerRecord> >::iterator>::iterator i=toRemove.begin();i!=toRemove.end();i++){
			list<SmartPtr<kmerRecord> >::iterator removeeI=*i;
			kmerRecord* removee=*removeeI;
			removee->kmerDisappearFromSeqs();
			foregroundTotalKmerCount-=removee->fgInstances(); //update foreground total kmer count
			kmers.erase(removeeI);
			delete removee;//now can remove this from member
		}
		
		cerr<<toRemove.size()<<"unique backgroundcount=0 kmers removed"<<endl;
		
		#endif
		
		//now we can remove the map kmerInitStruct to save space
		kmerInitStruct.clear();		
	}
	
	kmerRecord* getNextEnrichedKmers(){
		if(kmers.size()==0)
			return NULL;
		
		kmers.sort(); //sort the linked list by enrichment score
		kmerRecord*result=kmers.back();
		kmers.pop_back();
		return result;
	}
	
	void onNewSeq(const string&name,const string&seq,const string&quals){
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
				
				this->foregroundTotalKmerCount++;
				
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
			
			
		}else {
			
			numSeqBackground++;
			
			if(numSeqBackground%100000==1){
				cerr<<"\treading seq\t"<<numSeqBackground<<endl;
			}
			
			for(int i=0;i<=seq.length()-k;i++){
				string thisKmerSeq(seq.substr(i,k));
				
				if(!isValidSeq(thisKmerSeq)){
					continue;
				}
				
				this->backgroundTotalKmerCount++;
				
				SmartPtr<kmerRecord> thisKmer;
				map<string,SmartPtr<kmerRecord> >::iterator kmerI=kmerInitStruct.find(thisKmerSeq);
				if(kmerI!=kmerInitStruct.end()){
					//only when this kmer is present in foreground is it counted
					
					thisKmer=kmerI->second;
					thisKmer->backgroundCount++;
				}
				
				
			}		
			
			
			
			
		
		}

	}
	
	~kmerFinder(){
		for(list<SmartPtr<kmerRecord> >::iterator i=kmers.begin();i!=kmers.end();i++){
			delete (*i);
		}
		for(vector<SeqRecord*>::iterator i=seqs.begin();i!=seqs.end();i++){
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

class TestSimpleFastxReader:public SimpleFastqReader{
public:
	void onNewSeq(const string& name,const string& seq,const string& quals)
	{
		cerr<<"got new seq! name="<<name<<" seq="<<seq<<" qual="<<quals<<endl;
	}
};

int main(int argc,char**argv){
	
	
	/*TestSimpleFastxReader f;
	f.readFile(argv[1]);
	return 0;*/
	
	cerr<<"Find Kmers "<<VERSION<<endl;
	cerr<<"(";
#ifdef REQUIRE_BACKGROUND_PRESENCE
	cerr<<"REQUIRE_BACKGROUND_PRESENCE=ON)"<<endl;
#else
	cerr<<"REQUIRE_BACKGROUND_PRESENCE=OFF)"<<endl;
#endif
	
	if(argc<5){ 
		cerr<<"Usage:" <<argv[0]<<" <fgfilename> <bgfilename> <k> <howmanyToFind>"<<endl;
		cerr<<"Description: find <howmanyToFind> top enriched <k>-mers from HT-SELEX experiments using <fgfilename> and <bgfilename> fastq files assuming there's only one motif per sequence"<<endl;
		cerr<<"Specify <howmanyToFind>=0 to print all kmers"<<endl;
		return 1;
	}
	
	int k=atoi(argv[3]);
	int howmanyToFind=atoi(argv[4]);
	string fgfilename=argv[1];
	string bgfilename=argv[2];
	
	kmerFinder theFinder(k);
	
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
	cout<<"kmer\tenrichment\tnormalizedEnrichment\tfgInstances\tbgInstances"<<endl;
	
	if(howmanyToFind<1){
		howmanyToFind=theFinder.numUniqKmers();
	}
	
	for(int i=0;i<howmanyToFind;i++){
		kmerRecord* nextKmer=theFinder.getNextEnrichedKmers();
		if(!nextKmer)
			break;
		cout<<nextKmer->kmerSeq<<"\t"<<nextKmer->enrichment()<<"\t"<<nextKmer->normalizedEnrichment(theFinder.foregroundTotalKmerCount,theFinder.backgroundTotalKmerCount)<<"\t"<<nextKmer->fgInstances()<<"\t"<<nextKmer->backgroundCount<<endl;
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
