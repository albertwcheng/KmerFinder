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
#include "SimpleFastqReader.h"
#include "CppUtilClasses/FileIPMLoop.h"
#include "CppUtilClasses/StringUtil.h"
#include "CppUtilClasses/SmartPtr.h"
#include "defaultValues.h"



#ifdef USE_MULTISET
#define setType multiset
#else
#define setType set
#endif 

#define int64 int64_t

#define VERSION "D/1.0"

class SeqRecord;


class kmerRecord
{
public:

	string kmerSeq;
	setType<SeqRecord*> seqs; //this is the pointer to the seqs that have this kmer 
	
	int64 prevCount;

	
	kmerRecord(const string& _kmerSeq):kmerSeq(_kmerSeq),prevCount(0){
		
	}
	
	inline int64 NumOfInstances() const{
		return this->seqs.size();
	}
	
	//do we need sorting?
	inline bool operator < (const kmerRecord& right) const{
		return this->NumOfInstances()<right.NumOfInstances();
	
	}
	
	void printUpdateAndUpdateCountHistory(ostream& os){
		int64 numI=NumOfInstances();
		if(prevCount==numI) //no need to update nor print
			return;
		
		int64 diff=numI-prevCount;
		os<<kmerSeq<<"\t"<<diff<<endl;
		
		prevCount=numI; //record this new point
		
	}
	
	//void kmerDisappearFromSeqs();
	
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

	void registerKmer(kmerRecord* _kmer){
		_kmer->seqs.insert(this);
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
	
	
}

/*void kmerRecord::kmerDisappearFromSeqs(){
	for(setType<SeqRecord*>::iterator i=seqs.begin();i!=seqs.end();i++){
		SeqRecord* thisSeq=(*i);
		thisSeq->unregisterKmer(this);
	}
	

}*/



class DKmerCounter:public SimpleFastqReader,public FileIPMLoop
{
public:
	
	vector<SeqRecord*> seqs;

	int k; //wordsize
	int cycle;
	
	int64 totalKmerInstances;
	int64 totalNumSeqs;

	string outDir;
	
	map<string,SmartPtr<kmerRecord> > kmers; //for storing kmers while reading file initially
	
	
	DKmerCounter(string _outDir,string _name,int _delayTimeInSeconds,int _k):
		FileIPMLoop(_outDir+IPMSubFolder,_name,_delayTimeInSeconds),outDir(_outDir),
		k(_k),cycle(1),totalKmerInstances(0),totalNumSeqs(0){
	}
	
	inline int64 numUniqKmers(){
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
	
	void postUpdate(){
		
		//write Update to file
		if(!SystemUtil::fexists(outDir+DS+KMERUPDATEPATH)){
			cerr<<"makedir "<<(outDir+DS+KMERUPDATEPATH)<<endl;
			SystemUtil::mkdirs(outDir+DS+KMERUPDATEPATH);
		}
		
		ofstream fil((outDir+DS+KMERUPDATEPATH+DS+name+".txt").c_str());
		for(map<string,SmartPtr<kmerRecord> >::iterator i=kmers.begin();i!=kmers.end();i++){
			i->second->printUpdateAndUpdateCountHistory(fil);
		}
		fil<<TERMINATOR<<endl;
		fil.close();
		
		//send message
		sendMessage(MASTERID,StringUtil::str(cycle)+"\t"+"KmerCountUpdate");
		
		//update cycle number
		cycle++;
	}
	
	void removeKmer(const string& kmerToRemove){
		map<string,SmartPtr<kmerRecord> >::iterator i=kmers.find(kmerToRemove);
		if(i==kmers.end()){
			//nothing to remove
			return;	
		}
		
		i->second->removeSequencesContainingThisKmer();
	}
	
    void loadSeqFile(const string&filename){
		
		
		this->readFile(filename);		
		
		//done loading stuff
		postUpdate();
	}
	
	
	
	void onReceivingMessage(string sender,vector<string>& msgs){
		for(vector<string>::iterator i=msgs.begin();i!=msgs.end();i++)
			cerr<<sender<<":"<<(*i)<<endl;
			
		if(sender==MASTERID){
			if(msgs.size()>1){
				cerr<<"msg size > 1"<<endl;
				sendMessage(MASTERID,"Error\tmsg size > 1");
				terminate(); //fatal error	
			}else{
				string message=msgs[0];
				vector<string> splits;
				StringUtil::split(message,"\t",splits);
				if(StringUtil::atoi(splits[0])!=cycle){

					cerr<<"cycle number not equal from master"<<endl;
					sendMessage(MASTERID,"Error\tcycle number not equal from master");
					terminate();
					return;
				}
				if(splits[1]=="RemoveKCSAndUpdate"){
					if(splits.size()!=3){
						
						cerr<<"RemoveKCSAndUpdate <kmerToRemove> not given"<<endl;
						sendMessage(MASTERID,"Error\tRemoveKCSAndUpdate <kmerToRemove> not given");
						terminate();	
						return;
					}else{
						//ok
						cerr<<"remove Kmer "<<splits[2]<<endl;
						removeKmer(splits[2]);
						cerr<<"post update"<<endl;
						postUpdate();	
					}
				}
					
			}	
		}
	}
	
	
	~DKmerCounter(){
		for(map<string,SmartPtr<kmerRecord> >::iterator i=kmers.begin();i!=kmers.end();i++){
			delete i->second;
		}
		for(vector<SeqRecord*>::iterator i=seqs.begin();i!=seqs.end();i++){
			delete (*i);
		}

	}
	
	void printStat(ostream& os){
		os<<"#Seqs="<<totalNumSeqs<<endl;
		os<<"#k="<<k<<endl;
		os<<"#KmerInstances="<<totalKmerInstances<<endl;
		os<<"#UniqKmerCount="<<numUniqKmers()<<endl;
	}
	
	
	bool onNewSeq(const string&name,const string&seq,const string&quals)
	{
			
		totalNumSeqs++;
			
		if(totalNumSeqs%100000==1){
			cerr<<"\treading seq\t"<<totalNumSeqs<<endl;
		}
			
		SeqRecord *thisSeqRecord=new SeqRecord;
			
		seqs.push_back(thisSeqRecord);
			
			//for foreground fill the kmers;
		for(int i=0;i<=seq.length()-k;i++)
		{
			string thisKmerSeq(seq.substr(i,k));
				
			if(!isValidSeq(thisKmerSeq))
			{
				continue;
			}
				
				//this->foregroundTotalKmerCount++;
				
				//now insert this kmer
			SmartPtr<kmerRecord> thisKmer;
				
			map<string,SmartPtr<kmerRecord> >::iterator kmerI=kmers.find(thisKmerSeq);
			if(kmerI==kmers.end()){
				thisKmer=new kmerRecord(thisKmerSeq);
				kmers.insert(map<string,SmartPtr<kmerRecord> >::value_type(thisKmerSeq,thisKmer));
			}
			else 
			{
				thisKmer=kmerI->second;
			}
				
				//register thisSeqRecord in thisKmer
				
			thisSeqRecord->registerKmer(thisKmer);
				
		}
			
		return true;

	}

	
	
};



int main(int argc,char**argv){
	

	cerr<<"DKmerFinderCounter (Slave) "<<VERSION<<endl;
	
	if(argc<5){ 
		cerr<<"Usage:" <<argv[0]<<" <outDir> <name> <filename> <k>"<<endl;
		return 1;
	}
	
	string outDir=argv[1];
	string name=argv[2];
	string filename=argv[3];
	int k=atoi(argv[4]);

	DKmerCounter counterLoop(outDir,name,SLAVEDELAYTIME,k);
	counterLoop.sendMessage("master","1\tHi");
	counterLoop.loadSeqFile(filename);
	counterLoop.loop();
	
	return 0;
}
