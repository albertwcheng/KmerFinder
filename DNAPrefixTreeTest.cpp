#include "DNAPrefixTree.h"
#include <iostream>
#include <fstream>
//#include <stringstream>
#include <stdlib.h>
using namespace std;

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


class DNATreeKmerFinder:public SimpleFastqReader,public DNATree{
	public:
	uint64 numFgSeq;
	uint64 numBgSeq;
	DNATreeKmerFinder(int _maxK):DNATree(_maxK),numFgSeq(0),numBgSeq(0){}

	void readForegroundFile(string filename){
		DNATree::countMode=MODE_FOREGROUND;
		readFile(filename);
	}
	
	void readBackgroundFile(string filename){
		DNATree::countMode=MODE_BACKGROUND;
		readFile(filename);	
	}
	
	bool onNewSeq(const string& name,const string& seq,const string& quals){
		feedSeq(seq);
		
		if(DNATree::countMode==MODE_FOREGROUND){
			numFgSeq++;
			if(numFgSeq%100000==1){
				cerr<<"processed "<<numFgSeq<<" fg sequences. numNodes="<<numNodes<<endl;
			}	
		}else{
			numBgSeq++;	
			if(numBgSeq%100000==1){
				cerr<<"processed "<<numBgSeq<<" bg sequences. numNodes="<<numNodes<<endl;
			}	
		}
		return true;
	}
};

int main(int argc,const char **argv){
	//DNATree tree(5);
	//tree.feedSeq("ACGTAANC");
	
	//tree.feedSeq("ACGTCCCATGCCAT");
	//tree.printTree();
	
	if(argc<5){
		cerr<<"Usage:"<<argv[0]<<" fgfile bgfile maxK outDir"<<endl;
		return 1;	
	}
	
	

	int maxK=atoi(argv[3]);
	string outDir=argv[4];
	DNATreeKmerFinder finder(maxK);
	finder.readForegroundFile(argv[1]);
	finder.readBackgroundFile(argv[2]);
	
	//finder.printTree();
	cerr<<"sorting lists"<<endl;
	finder.sortLists();
	cerr<<"output"<<endl;
	char ist[10];
	for(int i=1;i<=maxK;i++){
		sprintf(ist,"%d",i);
		string filename=outDir+"/"+ist+".txt";
		ofstream fil(filename.c_str());	
		finder.printLevelNodes(fil,i,double(finder.numFgSeq)/finder.numBgSeq);
		fil.close();
	}
	//cerr<<finder.maxK<<endl;
	
	return 0;	
}