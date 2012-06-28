
#define uint64 uint64_t

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