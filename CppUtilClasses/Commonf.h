#ifndef COMMONF_H_
#define COMMONF_H_


#include <limits.h>
#include <limits>
#include "CppUtilClasses/StringUtil.h"
#include "Nucleic.h"
#include <sstream>
#include <vector>
using namespace std;

int getActualPosGappedSeq(const char*seq, const int* sorted_pos, int *out_pos, int num_pos, char gap='-');

template<class P>
class less_ptr
{
public:
	inline bool operator()(const P left, const P right)
	{
		return *left<*right;
	}
};

inline static string toUpperAndNoSpecialDNA(const string& str)
{
	int len=str.length();
	char *tmp=new char[len+1];
	register const char*cur=str.c_str();
	register const char*terminal=cur+len;
	register char *tmp2=tmp;

	for(;cur<terminal;cur++,tmp2++)
	{
		*tmp2=toupper(*cur);
		switch(*tmp2)
		{
		case A_SYM:case C_SYM:case T_SYM:case G_SYM:
			break;
		case U_SYM:
			*tmp2=T_SYM;
			break;
		default:
			*tmp2=N_SYM;
			
		}
	}

	tmp[len]='\0';
	string tstr=string(tmp);
	delete[] tmp;
	return tstr;
}



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



void die_clean_up();
void die_exit_clean_up();

#define die(msg) {if(strlen(msg)>0){ cerr<<"die:"<<(msg)<<endl;} die_clean_up(); return; }/*added die_clean_up defined  to clean up memory*/
#define dier(msg,r) {if(strlen(msg)>0) { cerr<<"die:"<<(msg)<<endl;} die_clean_up(); return (r); }/*added reset everything to clean up memory*/
#define die_exit(msg) {if(strlen(msg)>0) { cerr<<"die:"<<(msg)<<endl;} die_exit_clean_up(); exit(1); }

#define MAX(a,b) (((a)>(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))


inline KeyPair<int,int> b01tob00(const KeyPair<int,int>& b01)
{
	return KeyPair<int,int>(b01.k1,b01.k2-1);
}

inline KeyPair<int,int> b00tob01(const KeyPair<int,int>& b00)
{
	return KeyPair<int,int>(b00.k1,b00.k2+1);
}


inline int offset00(int a0,int len)
{
	return a0+len-1;
}

inline int offset01(int a0,int len)
{
	return a0+len;
}

inline int offset10(int a1,int len)
{
	return a1+len-2;
}

inline int offset11(int a1,int len)
{
	return offset00(a1,len);
}



inline bool isBound00Valid(const KeyPair<int,int>& b00)
{
	return b00.k1>=0 && b00.k2<INT_MAX && b00.k1<=b00.k2;
}

inline bool isBound11Valid(const KeyPair<int,int>& b11)
{
	return b11.k1>=1 && b11.k2<INT_MAX && b11.k1<=b11.k2;
}

inline bool isBound01Valid(const KeyPair<int,int>& b01)
{
	return b01.k1>=0 && b01.k2<INT_MAX && b01.k1<b01.k2;
}

inline int len01(int start0,int end1)
{
	
	int l=end1-start0;
	return MAX(l,0);
}

inline int len00(int start0,int end0)
{
	int l=end0-start0+1;
	return MAX(l,0);
}

inline int len10(int start1,int end0)
{
	
	int l=end0-start1+2;
	return MAX(l,0);
}

inline int len11(int start1,int end1)
{
	return len00( start1, end1);
}

inline int len01(const KeyPair<int,int>& bound01)
{
	return len01(bound01.k1,bound01.k2);
}

inline int len11(const KeyPair<int,int>& bound11)
{
	return len11(bound11.k1,bound11.k2);
}


inline int len10(const KeyPair<int,int>& bound10)
{
	return len10(bound10.k1,bound10.k2);
}
inline int len00(const KeyPair<int,int>& bound00)
{
	return len00(bound00.k1,bound00.k2);
}


inline int overlaplen1111(const KeyPair<int,int>& a11,const KeyPair<int,int> &b11)
{
	if(!isBound11Valid(a11) || !isBound11Valid(b11))
		return 0;

	return len11(MAX(a11.k1,b11.k1),MIN(a11.k2,b11.k2));
	//return MAX(lraw,0);
}


inline int overlaplen0101(const KeyPair<int,int>& a01,const KeyPair<int,int> &b01)
{
	if(!isBound01Valid(a01) || !isBound01Valid(b01))
		return 0;

	return len01(MAX(a01.k1,b01.k1),MIN(a01.k2,b01.k2));
	//return MAX(lraw,0);
}


template<class T>
class SerializableVector
{
public:
	vector<T> _vector;
	inline SerializableVector(string filename)
	{
		load(filename);
	}
	inline int load(string filename,bool clear=true)
	{
		if(clear)
		{
			_vector.clear();
		}
		int numRead=0;
		
		ifstream fin(filename.c_str());
		
		if(!fin.good())
		{
			fin.close();
			return numRead;
		}
		
		while(!fin.eof())
		{
			T t;
			fin>>t;
			_vector.push_back(t);
			numRead++;
		}
		
		fin.close();
		return numRead;
	}
	inline int write(string filename,std::_Ios_Openmode nMode=ios::out|ios::trunc,string newLine="\n")
	{
		int numWrite=0;
		ofstream fout(filename.c_str(),nMode);
		if(!fout.good())
			return numWrite;
		
		for(typename vector<T>::iterator i=_vector.begin();i!=_vector.end();i++)
		{
			fout<<(*i);
			if(newLine.length()>0)
				fout<<newLine;
				
		}
		
		numWrite=_vector.size();
		return numWrite;
	}
};

template< class K, class V>
class One2OneMapping
{
public:
	map<K,V> _map;
	

	
	inline int load(string filename,bool clearMap=true)
	{
		if(clearMap)
		{
			_map.clear();
		}
		int numRead=0;
		
		ifstream fin(filename.c_str());
		
		
		if(!fin.good())
		{	
			fin.close();
			return numRead;
		}
		
		
		while(!fin.eof())
		{
			K k;
			V v;
			fin>>k;
			if(fin.eof())
			{
				break;
			}
			fin>>v;
			_map.insert(typename map<K,V>::value_type(k,v));
			numRead++;
		}
		
		fin.close();
		return numRead;
	}
	
	inline One2OneMapping(string filename)
	{
		load(filename);
	}
	
};

typedef One2OneMapping<int,string> ChrMap;
typedef SerializableVector<string> FileList;


template< class T, class K, class I>
class buffered_ofstream
{
public:
	string filename;
	string endline;
	//ofstream* fout;
	K store;
	int threshold;
	int storedItemSize;
	int lineOutputed;
	
	typedef I ItemIterator;
	
	inline buffered_ofstream(string _filename,bool clearFil=true,int _threshold=500,string _endline="\n"): filename(_filename), threshold(_threshold), store(_threshold), endline(_endline),storedItemSize(0),lineOutputed(0)
	{
		if(clearFil)
			clearFile();
	}
	
	inline void clearFile()
	{
		//just to clear the content of the file
		lineOutputed=0;
		ofstream fout(filename.c_str());
		fout.close();
	}
	
	inline void flush(std::_Ios_Openmode nMode=ios::out|ios::app,bool printEndLine=true)
	{
		ofstream fout(filename.c_str(),nMode);
		if(!fout.good())
		{	
			cerr<<"file error"<<endl;
			fout.close();
			return;
		}
		for(int i=0;i<storedItemSize;i++)
		{
			if(!fout.good())
			{	
				cerr<<"file error"<<endl;
				fout.close();
				return;
			}
			
			if(printEndLine)
				fout<<store[i]<<endline;
			else
				fout<<store[i];
			//cerr<<(*i)<<endline;
		}
		lineOutputed+=storedItemSize;
		//store.clear();
		storedItemSize=0;
		fout.close();
	}
	
	inline void push(const T& element,std::_Ios_Openmode nMode=ios::out|ios::app,bool printEndLine=true)
	{
		store[storedItemSize]=element;
		storedItemSize++;
		if(storedItemSize==threshold)
		{
			flush(nMode,printEndLine);
		}
	}
	
	inline int getStoredItemSize()
	{
		return this->storedItemSize;
	}
	inline int getLineOutputed()
	{
		return this->lineOutputed;
	}
	
	
	
};

#define READLINEBUFFERLENGTH 3

 inline  string consumeStreamUntilNewLine(istream &is)
{
	/*char newLineBuff[READLINEBUFFERLENGTH];
	stringstream strin;
	bool cont=true;
	while(cont)
	{

		is.getline(newLineBuff,READLINEBUFFERLENGTH);
		strin<<newLineBuff;
		cerr<<"read "<<newLineBuff<<" gcount="<<is.gcount()<<" flag="<<is.rdstate()<<endl;
		if( (is.rdstate() & ios::failbit)!=0 && is.gcount()==READLINEBUFFERLENGTH )
		{
			cont=true;
		}
		else cont=false;
	}


	return strin.str();*/
	 string s;

	 std::getline(is,s);
	 return s;
}


#endif /*COMMONF_H_*/
