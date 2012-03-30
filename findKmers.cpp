#include <fstream>
#include <string>
#include <iostream>
#include <map>
#include <list>
#include <set>
#include <vector>
#include <inttypes.h>
#include <stdlib.h>
using namespace std;


#define uint64 uint64_t

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
			fin>>line; //how to read whole line?
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
	
	inline double enrichment() const{
		return double(this->seqs.size())/(backgroundCount+1); //background plus 1 to avoid division by zero
		//no need to care about the totalkmers in foreground and background becoz they are the same for all kmers. argmax not important to know these numbers
	}
	
	inline bool operator < (const kmerRecord& right) const{
		return this->enrichment()<right.enrichment();
	
	}
	
	
	
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
};


void kmerRecord::removeSequencesContainingThisKmer(){
	for(multiset<SeqRecord*>::iterator i=seqs.begin();i!=seqs.end();i++){
		SeqRecord* thisSeq=(*i);
		thisSeq->disappearFromKmers();
	}
	
}

#define READ_FOREGROUND 1
#define READ_BACKGROUND 2




class kmerFinder:public SimpleFastxReader
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
	}
	
	kmerRecord* getNextEnrichedKmers(){
		kmers.sort(); //sort the linked list by enrichment score
		kmerRecord*result=kmers.back();
		kmers.pop_back();
		return result;
	}
	
	void onNewSeq(const string&name,const string&seq,const string&quals){
		if(mode==READ_FOREGROUND){
			
			numSeqForeground++;
			
			if(numSeqForeground%1000==1){
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
			
			if(numSeqBackground%1000==1){
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
		os<<"Foreground Seq Count="<<numSeqForeground<<endl;
		os<<"Foreground Kmer Count="<<foregroundTotalKmerCount<<endl;
		os<<"Background Seq Count="<<numSeqBackground<<endl;
		os<<"Background Kmer Count="<<backgroundTotalKmerCount<<endl;
	}
	
};

class TestSimpleFastxReader:public SimpleFastxReader{
public:
	void onNewSeq(const string& name,const string& seq,const string& quals)
	{
		cerr<<"got new seq! name="<<name<<" seq="<<seq<<" qual="<<quals<<endl;
	}
};

int main(int argc,char**argv){
	
	if(argc<5){ 
		cerr<<"Usage:" <<argv[0]<<" fgfilename bgfilename k howmanyToFind"<<endl;
		return 1;
	}
	//TestSimpleFastxReader f;
	//f.readFile(argv[1]);
	
	int k=atoi(argv[3]);
	int howmanyToFind=atoi(argv[4]);
	string fgfilename=argv[1];
	string bgfilename=argv[2];
	
	kmerFinder theFinder(k);
	cerr<<"reading foreground "<<fgfilename<<endl;
	theFinder.readForeground(fgfilename);
	cerr<<"reading background "<<bgfilename<<endl;
	theFinder.readBackground(bgfilename);
	cerr<<"done reading"<<endl;
	theFinder.printStat(cerr);
	cerr<<"processing..."<<endl;
	for(int i=0;i<howmanyToFind;i++){
		kmerRecord* nextKmer=theFinder.getNextEnrichedKmers();
		cout<<"The first kmer="<<nextKmer->kmerSeq<<" enrichment="<<nextKmer->enrichment()<<endl;
	}
	cerr<<"done"<<endl;
	
	return 0;
	
}
