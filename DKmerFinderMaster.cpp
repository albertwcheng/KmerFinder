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



class kmerRecord
{
public:

	string kmerSeq;
		
	int64 prevCount;
	
	int64 fgInstances;
	int64 bgInstances;
	
	kmerRecord(const string& _kmerSeq):kmerSeq(_kmerSeq),fgInstances(0),bgInstances(0){
		
	}
	
	inline float enrichment() const{
		if(bgInstances==0){
			return -1; //unknown	
		}
		return float(fgInstances)/bgInstances;
	}
	
	//do we need sorting?
	inline bool operator < (const kmerRecord& right) const{
		return this->enrichment()<right.enrichment();
	
	}
	
	
};


class DKmerFinder:public FileIPMLoop
{
public:
	
	
	int k; //wordsize
	int cycle;
	int respondents;
	int topN;
	
	uint64 tick;
	
	int64 totalKmerInstances;
	int64 totalNumSeqs;

	string outDir;
	
	map<string,SmartPtr<kmerRecord> > kmersMap; //for storing kmers searchable by key
	list< SmartPtr<kmerRecord> > kmers;
	
	vector<string> fgSlaves;
	vector<string> bgSlaves;
	
	int totalSlaves;
	
	int hirespondents;
	
	ofstream *fout;
	
	
	DKmerFinder(string _outDir,string _name,int _delayTimeInSeconds,int _k,string _fgSlaves,string _bgSlaves,int _topN):
		FileIPMLoop(_outDir+IPMSubFolder,_name,_delayTimeInSeconds),outDir(_outDir),
		k(_k),cycle(1),totalKmerInstances(0),totalNumSeqs(0),respondents(0),topN(_topN),hirespondents(0),tick(0)
	{
		StringUtil::split(_fgSlaves,",",fgSlaves);
		StringUtil::split(_bgSlaves,",",bgSlaves);
		totalSlaves=fgSlaves.size()+bgSlaves.size();
		fout=new ofstream((outDir+DS+"kmers.txt").c_str());
	}
	
	inline void updateCycleAndResetRespondents(){
		cycle++;
		respondents=0;	
	}	
	
	//inline void keepOnlyKmersOccuringInBackground(){
	//	//TODO:	
	//}
	
	inline void updateKmerWithDiff(const string&kmerSeq,int64 diff,bool mode,bool allowNewKmerAddition){
		map<string,SmartPtr<kmerRecord> >::iterator i=kmersMap.find(kmerSeq);
		SmartPtr<kmerRecord> thisKmer;
		if(i==kmersMap.end()){
			if(!allowNewKmerAddition)
				return;
			thisKmer=new kmerRecord(kmerSeq);
			kmersMap.insert(map<string,SmartPtr<kmerRecord> >::value_type(kmerSeq,thisKmer));
			kmers.push_back(thisKmer);	
		}else{
			thisKmer=i->second;	
		}
		
		if(mode==UPDATE_FG){
			thisKmer->fgInstances+=diff;
		}else{
			thisKmer->bgInstances+=diff;	
		}
	}
	
	void loadKmerUpdateFile(string filename,bool mode,bool allowNewKmerAddition){
		ifstream fil(filename.c_str());
		
		while(fil.good()){
			string kmerSeq;
			fil>>kmerSeq;
			if(kmerSeq==TERMINATOR){
				break;	
			}
			int64 kmerDiff;
			fil>>kmerDiff;
			
			updateKmerWithDiff(kmerSeq, kmerDiff, mode ,allowNewKmerAddition);
		}
		
		fil.close();	
	}
	
	void terminateAllSlaves(){
		for(vector<string>::iterator i=fgSlaves.begin();i!=fgSlaves.end();i++){
			sendMessage(*i,"!TERMINATE");	
		}
		for(vector<string>::iterator i=bgSlaves.begin();i!=bgSlaves.end();i++){
			sendMessage(*i,"!TERMINATE");	
		}		
	}
	
	kmerRecord* sortAndFindTopKmer(){
		kmers.sort(); //sort the linked list by enrichment score
		kmerRecord*result=kmers.back();
		kmers.pop_back(); //remove the best kmer.
		kmersMap.erase(result->kmerSeq); //remove that also from the kmer map
		return result;
	}
	
	void askSlavesToRemoveAndUpdate(string kmerSeq){
		for(vector<string>::iterator i=fgSlaves.begin();i!=fgSlaves.end();i++){
			sendMessage(*i,StringUtil::str(cycle)+"\t"+"RemoveKCSAndUpdate"+"\t"+kmerSeq);	
		}
		for(vector<string>::iterator i=bgSlaves.begin();i!=bgSlaves.end();i++){
			sendMessage(*i,StringUtil::str(cycle)+"\t"+"RemoveKCSAndUpdate"+"\t"+kmerSeq);	
		}				
	}
	
	void onLoop(){
		tick++;
		if(tick>HIDEADLINE && hirespondents<totalSlaves){
			cerr<<"Not all slaves said hi before a deadline of "<<HIDEADLINE<<" seconds"<<endl;
			cerr<<"Abort"<<endl;
			cerr<<"Sending all slaves !TERMINATE"<<endl;
			terminateAllSlaves();
			terminate();
		}
	}
	
	
	void onReceivingMessage(string sender,vector<string>& msgs){
		
		for(vector<string>::iterator i=msgs.begin();i!=msgs.end();i++)
			cerr<<sender<<":"<<(*i)<<endl;
			
		int mode;
		if(sender.substr(0,2)=="fg"){
			mode=UPDATE_FG;
		}else if(sender.substr(0,2)=="bg"){
			mode=UPDATE_BG;	
		}else{
			//unknown sender. Do nothing
			return;	
		}
		
		vector<string> splits;
		
		for(vector<string>::iterator i=msgs.begin();i!=msgs.end();i++)
		{
			string& msg=*i;
			StringUtil::split(msg,"\t",splits);
			int msgCycleNumber=StringUtil::atoi(splits[0]);
			if(splits[1]=="Error"){
				cerr<<"Error occured at "<<sender<<" : "<<splits[2]<<endl;	
				cerr<<"Abort"<<endl;
				cerr<<"Sending all slaves !TERMINATE"<<endl;
				terminateAllSlaves();
				terminate(); //terminate myself too
			}else if(splits[1]=="Hi"){
				hirespondents++;
			}else if(splits[1]=="KmerCountUpdate"){
				
				if(msgCycleNumber!=cycle){
					cerr<<"Inconsistent cycle number from slave "<<sender<<" msgCycleNumber="<<msgCycleNumber<<" cycleHere="<<cycle<<endl;
					cerr<<"Abort"<<endl;
					cerr<<"Sending all slaves !TERMINATE"<<endl;
					terminateAllSlaves();
					terminate();	
				}else{
					cerr<<"loadKmerUpdate from "<<sender<<" at cycle "<<cycle<<"…";
					loadKmerUpdateFile(outDir+DS+KMERUPDATEPATH+DS+sender+".txt",mode,cycle==1); //only at cycle one allowing new kmer to be added.	
					cerr<<"done"<<endl;
					respondents++;
					if(respondents==totalSlaves){
						cerr<<"all slaves updated kmer counts at cycle "<<cycle<<endl;
						cerr<<"now find top enriched kmer at cycle "<<cycle<<endl;
						kmerRecord* topKmer=sortAndFindTopKmer();
						
						//output this topKmer
						string kmerSeq = topKmer->kmerSeq;
						
						if(cycle==1){
							(*fout)<<"kmer\tenrichment\tcontrol_count\texperim_count"<<endl;	
						}
						(*fout)<<kmerSeq<<"\t"<<topKmer->enrichment()<<"\t"<<topKmer->bgInstances<<"\t"<<topKmer->fgInstances<<endl;
						
						delete topKmer; //free it from memory
						
						if(topN > 0 && cycle==topN){
							//we are done
							cerr<<topN<<" top kmer(s) have been outputted. End slaves and master"<<endl;
							cerr<<"Sending all slaves !TERMINATE"<<endl;
							terminateAllSlaves();
							terminate();
						}else{
							
							//advance
							updateCycleAndResetRespondents();
							askSlavesToRemoveAndUpdate(kmerSeq);
						}
						
						
					}
				}
			}
		}

		
	}
	
	
	~DKmerFinder(){
		for(map<string,SmartPtr<kmerRecord> >::iterator i=kmersMap.begin();i!=kmersMap.end();i++){
			delete i->second;
		}

		
		if(fout){
			fout->close();
			delete fout;	
			
		}

	}
	


	
	
};



int main(int argc,char**argv){
	

	cerr<<"DKmerFinder (Master) "<<VERSION<<endl;
	
	if(argc<5){ 
		cerr<<"Usage:" <<argv[0]<<" <outDir> <fgNames> <bgNames> <k> <topN>"<<endl;
		return 1;
	}
	
	string outDir=argv[1];
	string fgNames=argv[2];
	string bgNames=argv[3];
	int k=atoi(argv[4]);
	int topN=atoi(argv[5]);

	DKmerFinder masterLoop(outDir,MASTERID,SLAVEDELAYTIME,k,fgNames,bgNames,topN); 
	masterLoop.loop();
	
	return 0;
}
