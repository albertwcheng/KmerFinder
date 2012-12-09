#include "DNAPrefixTree.h"
#include <iostream>
#include <fstream>
//#include <stringstream>
#include <stdlib.h>
#include <math.h>
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

void fillSubRow_inner( DNATree& finder,string& topKmerSeqSub,int j,char sub,vector<double>& row){
	topKmerSeqSub[j]=sub;
	DNATreeNode* subNode=finder.getNodeByPath(topKmerSeqSub);
	if(subNode)
		row.push_back(subNode->enrichment());
	else
		row.push_back(0.0); //subnode not present

}


void fillSubRow(DNATree& finder,const string& topKmerSeq,vector<double>& row,int j)
{
	string topKmerSeqSub=topKmerSeq;
			
	fillSubRow_inner(finder,topKmerSeqSub,j,'A',row);	
	fillSubRow_inner(finder,topKmerSeqSub,j,'C',row);
	fillSubRow_inner(finder,topKmerSeqSub,j,'G',row);
	fillSubRow_inner(finder,topKmerSeqSub,j,'T',row);
}

void normalizePWM(vector< vector<double> > &PWM, double pseudoweight){
	for(vector< vector<double > >::iterator i=PWM.begin();i!=PWM.end();i++){
		double sum=0.0;
		for(vector<double>::iterator j=i->begin();j!=i->end();j++){
			sum+=(*j);	
		}
		//now put new values
		for(vector<double>::iterator j=i->begin();j!=i->end();j++){
			(*j)=double(*j)/sum*pseudoweight+0.25*(1-pseudoweight);	
		}
		
	}
}

double InformationContentOfPWM(vector< vector<double> >& PWM,vector<double> backgroundFreq,double logb){
	double I=0.0;
	for(int i=0;i<PWM.size();i++){
		for(int j=0;j<4;j++){
			I+=PWM[i][j]*(-log(PWM[i][j]))/logb; //log(backgroundFreq[j])
		}	
	}
	
	return I;
	
}

void printPWMRow(ostream& os,vector<double>& row,bool withHeader=true){
	if(withHeader)
		os<<"A\tC\tG\tT"<<endl;
	os<<row[0];
	for(int j=1;j<4;j++){
		os<<"\t";
		os<<row[j];
				
	}	
	os<<endl;

}

void printPWM(ostream& os,vector<vector<double> > & PWM){
	os<<"A\tC\tG\tT"<<endl;
	for(int i=0;i<PWM.size();i++){
		printPWMRow(os,PWM[i],false);
	}	
}


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
	cerr<<"sort kmer lists by enrichment"<<endl;
	finder.sortLists();
	
	double pseudoweight=0.99;
	
	vector<double> backgroundFreq;
	double logb=log(2);
	
	
	//backgroundFreq.push_back(double(finder.getNodeByPath("A")->bgcount));
	//backgroundFreq.push_back(double(finder.getNodeByPath("C")->bgcount));
	//backgroundFreq.push_back(double(finder.getNodeByPath("G")->bgcount));
	//backgroundFreq.push_back(double(finder.getNodeByPath("T")->bgcount));
	
	for(int i=0;i<4;i++){
		backgroundFreq.push_back(0.25);	
	}
	
	
	//normalize backgroundFreq
	double bgsum=0.0;
	
	for(int i=0;i<4;i++){
		bgsum+=backgroundFreq[i];
	}
	
	for(int i=0;i<4;i++){
		backgroundFreq[i]/=bgsum;	
	}
	
	cerr<<"backgroundFreq"<<endl;
	printPWMRow(cerr,backgroundFreq);
	
	
	cerr<<"construct k-mer per k"<<endl;
	for(int i=0;i<maxK;i++){
		DNATreeNode* topKmerNode=finder.listsPerK[i].back();
		string topKmerSeq=topKmerNode->getPathFromRoot();
		vector< vector<double> > PWM;
		for(int j=0;j<topKmerSeq.length();j++){
			//substitute position j
			
			vector<double> row;
			
			fillSubRow(finder,topKmerSeq,row,j);
			
			PWM.push_back(row);
		} 	
		
		int k=PWM.size();
		
		//now normalize
		normalizePWM(PWM,pseudoweight);
		
		printPWM(cerr,PWM);
		
		//now calculate information content
		
		
		double I=InformationContentOfPWM(PWM,backgroundFreq,logb);
		double IpB=I/k;
		
		
		cerr<<"consensus "<<topKmerSeq<<" at k="<<k<<"; information="<<I<<" bits; "<<IpB<<" bits/nt"<<endl;
		
	}
	
	return 0;
	
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