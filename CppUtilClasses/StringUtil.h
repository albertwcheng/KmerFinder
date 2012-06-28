/*****
 
 Copyright (c) 2011, Albert W. Cheng
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 1. Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 3. All advertising materials mentioning features or use of this software
 must display the following acknowledgement:
 This product includes software developed by Albert W. Cheng.
 4. Neither the name of Albert W. Cheng nor the
 names of its contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY Albert W. Cheng ''AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL Albert W. Cheng BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 ******/


#ifndef STRINGUTIL_H_
#define STRINGUTIL_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
#include <map>
#include <ctype.h>
using namespace std;

#define BUFFER_SIZE 1025

#define THRESHOLD_TRIM_BUFFER_TIME 10
#define THRESHOLD_TRIM_BUFFER_COMP 10240


template<class T>
class buffer
{
private:
	T*data;
	int size;
	int smallerRequested;
public:
	inline T* buf()
	{
		return data;
	}
	int capacity()
	{
		return size;
	}
	buffer(int _size=BUFFER_SIZE):data(NULL),size(0)
	{
		buf(_size);
	}
	T* realloc(int _size)
	{
		if(data)
		{
			delete[] data;
		}

		data=new T[_size];
		size=_size;
		smallerRequested=0;

		return data;
	}
	T& operator[] (int i){
		return data[i];
	}

	T* buf(int _size)
	{
		//cerr<<"buffer for "<<_size<<endl;
		if(_size>size || data==NULL)
		{

			realloc(_size);

		}
		else if (size>_size+THRESHOLD_TRIM_BUFFER_COMP)
		{
			smallerRequested++;
			if(smallerRequested>=THRESHOLD_TRIM_BUFFER_TIME)
			{
				realloc(_size);

			}
		}
		return data;

	}
	~buffer()
	{
		delete[] data;
	}
};

template<class T1,class T2>
class KeyPair
{
public:
	T1 k1;
	T2 k2;
	inline KeyPair():k1(0),k2(0){}
	KeyPair(T1 _k1,T2 _k2): k1(_k1), k2(_k2) {
	}

	inline bool operator <(const KeyPair<T1,T2>& obj) const
	{
		if(this->k1==obj.k1)
		{
			return (this->k2<obj.k2);
		}

		return this->k1<obj.k1;

	}

	inline bool operator !=(const KeyPair<T1,T2>& obj) const
	{
		return !((*this)==obj);
	}

	inline bool operator >=(const KeyPair<T1,T2>& obj) const
	{
		return !((*this)<obj);
		//return (*this)>obj || (*this)==obj;
	}

	inline bool operator <=(const KeyPair<T1,T2>& obj) const
	{
		return !((*this)>obj);
		//return (*this)<obj || (*this)==obj;
	}

	bool operator ==(const KeyPair<T1,T2>& obj) const
	{
		return this->k1==obj.k1 && this->k2==obj.k2;
	}

	bool operator >(const KeyPair<T1,T2>& obj) const
	{
		if(this->k1==obj.k1)
		{
			return (this->k2>obj.k2);
		}

		return this->k1>obj.k1;

	}
};

template<class T1,class T2>
ostream& operator << ( ostream&os,const KeyPair<T1,T2>&k)
{
	os<<k.k1<<"\t"<<k.k2;
	return os;
}

class StringUtil
{
public:
	
	inline static string str(int d)
	{
		return str(d,"%d");
	}

	inline static string str(unsigned int d)
	{
		return str(int(d));
	}

	inline static string str(double d)
	{
		return str(d,"%f");
	}

	inline static string str(double d,int precision)
	{
		return str(d,string("%.")+str(precision)+"f");
	}

	template<class T>
	inline static string str(T d,string format)
	{
		return str(d,format.c_str());
	}

	template<class T>
	inline static string str(T d,const char* format)
	{
		char buffer[100];
		sprintf(buffer,format,d);
		return string(buffer);
	}

	inline static const char* isPrefix(const char*pprefix,const char*fullstring)
	{
		//pprefix : putative prefix
		if(*pprefix=='\0') 
			return fullstring; //emptystring is prefix of everything
		
		char prefixChar;
		char fullStringChar;
		while((prefixChar=*(pprefix++))!='\0')
		{	
			//prefix string not terminated yet
			fullStringChar=*(fullstring++);
			if(fullStringChar=='\0') //fullstring terminated before pprefix, so pprefix is not real prefix
				return NULL;
			
			if(prefixChar!=fullStringChar)
				return NULL;
			
		}
		
		//pprefix exhausted at the same time as fullstring or before fullstring means it is prefix
		
		return fullstring;
		
	}

	inline static string getCommonPrefix(const string& a, const string& b)
	{
		unsigned int i=0;
		const char * astr=a.c_str();
		const char * bstr=b.c_str();
		while(i<a.length() && i<b.length() && astr[i]==bstr[i])
		{
			i++;
		}

		return a.substr(0,i);
	}

	inline static bool isSuffix(const char*psuffix,const char*fullstring)
	{
		int lpsuffix=strlen(psuffix);
		int lfullstring=strlen(fullstring);
		if(lpsuffix>lfullstring)
			return false; //putative suffix is longer than full string! it is of coz not suffix
		
		return !strcmp(fullstring+lfullstring-lpsuffix,psuffix);

	}

	inline static int atoi(const char* s)
	{
		return ::atoi(s);
	}


	inline static double atof(const char *s)
	{
		return ::atof(s);
	}

	template<typename Container,typename Iterator>
	inline static string join(const Container& arr,string sep)
	{
		Iterator i=arr.begin();

		if(i==arr.end()) return "";

		string result=(*i);
		i++;
		for(;i!=arr.end();i++)
		{
			result+=sep+(*i);
		}

		return result;
	}

	inline static string formArrayString( vector<string>& arr,string sep="|")
	{
		vector<string>::iterator i=arr.begin();
		if(i==arr.end()) return "";
		string result=(*i);
		i++;
		for(;i!=arr.end();i++)
		{
			result+=sep+(*i);
		}
		return result;
	}



	inline static string formKeyValueString( map<string,string>& _map,string sep="|",string equ="=")
	{

		map<string,string>::iterator i=_map.begin();
		if(i==_map.end())
			return "";
		string result=(*i).first+equ+(*i).second;
		i++;
		for(;i!=_map.end();i++)
		{
			result+=sep+(*i).first+equ+(*i).second;
		}

		return result;
	}

	 inline static int atoi(string s)
	{
		return atoi(s.c_str());
	}


	inline static double atof(string s)
	{
		return atof(s.c_str());
	}
	inline static string toUpper(const string& str)
	{
		int len=str.length();
		char *tmp=new char[len+1];
		register const char*cur=str.c_str();
		register const char*terminal=cur+len;
		register char *tmp2=tmp;

		for(;cur<terminal;cur++,tmp2++)
		{
			*tmp2=toupper(*cur);
		}

		tmp[len]='\0';
		string tstr=string(tmp);
		delete[] tmp;
		return tstr;
	}
	inline static string escape(const string&str,const char* escapeSet="nt\0",const char* escapeTo="\n\t\0",char escapePrefix='\\')
	{
		int nOriStr=str.length();
		char *buffer=new char[nOriStr+1];
		//strcpy(buffer,str.c_str());
		const char*oriStr=str.c_str();		
		char *bufferWriter=buffer;
		int nEscapes=strlen(escapeSet);		

		for(int i=0;i<nOriStr;i++)
		{
			char oriStrCurChar=oriStr[i];
			if(oriStrCurChar==escapePrefix)
			{
				char escapeSpecifier=oriStr[i+1];
				int j;

				for( j=0;j<nEscapes;j++)
				{
					if(escapeSpecifier==escapeSet[j])
					{
						break;
					}
				}

				if(j==nEscapes)
				{

					(*bufferWriter++)=escapeSpecifier;
				}
				else
				{
					(*bufferWriter++)=escapeTo[j];
				}

				i++;
				
				
			}
			else
			{
				(*bufferWriter++)=oriStrCurChar;
			}

		}
				
		*bufferWriter='\0';

		string returnString(buffer);		
		delete[] buffer;
		return returnString;
	}
	inline static string toLower(const string&str)
	{
		int len=str.length();
		char *tmp=new char[len+1];
		register const char*cur=str.c_str();
		register const char*terminal=cur+len;
		register char *tmp2=tmp;

		for(;cur<terminal;cur++,tmp2++)
		{
			*tmp2=tolower(*cur);
		}

		tmp[len]='\0';
		string tstr=string(tmp);
		delete[] tmp;
		return tstr;
	}
	inline static string stripAll(string haystack,string remove)
	{
		const char * str=haystack.c_str();
		const char * removestr=remove.c_str();
		int rlen=remove.length();
		int len=haystack.length();

		int k=0;

		buffer<char> buff(len+1);
		char *bstr=buff.buf();
		for(int i=0;i<len;i++)
		{
			bool exclude=false;
			for(int j=0;j<rlen;j++)
			{
				if(str[i]==removestr[j])
				{
					exclude=true;

				}
			}

			if(!exclude)
				bstr[k++]=str[i];
		}

		bstr[k]='\0';
		return string(bstr);
	}


	
	static void splitInt2(string haystack,string sep,vector<int>& v,bool clear=true)
	{
		if(clear)
			v.clear();
		
		
		vector<string> spliton;
		split(haystack,sep,spliton);
		//v.resize(spliton.size());
		
		//int j=0;
		for(vector<string>::iterator i=spliton.begin();i!=spliton.end();i++)
			v.push_back(StringUtil::atoi(*i));
	
		
	}

	static vector<int> splitInt(string haystack,string sep="\t")
	{
		/*vector<string> spliton;
		split(haystack,sep,spliton);
		vector<int> v(spliton.size());
		int j=0;
		for(vector<string>::iterator i=spliton.begin();i!=spliton.end();i++)
			v[j++]=StringUtil::atoi(*i);*/
		
		vector<int> v;
		
		StringUtil::splitInt2(haystack,sep,v,true);

		return v;
	}
	
	inline static vector<string>& split(const string haystack,const string& needle,vector<string>& vec,bool clear=true)
	{

		int hslen=haystack.length();
		int needlen=needle.length();
		const char *shaystack=haystack.c_str();
		char *hay=new char[hslen+1];
		strcpy(hay,shaystack);
		const char *sneedle=needle.c_str();

		if(clear)
			vec.clear();

		int start=0;


		for(int i=0;i<hslen;i++)
		{
			char c=hay[i];
			bool found=false;
			for(int j=0;j<needlen;j++)
			{
				if (sneedle[j]==c)
					found=true;
			}
			if(found)
			{
				hay[i]='\0';
				vec.push_back((hay+start));
				hay[i]=' ';
				start=i+1;
			}
		}

		if(start<hslen)
		{
			vec.push_back( (shaystack+start));
		}else if(start==hslen) //new sep 14 2009
		{
			vec.push_back("");
		}


		delete[] hay;
		return vec;
	}

	inline static vector<string>& splitNoEmpty(const string haystack,const string& needle,vector<string>& vec,bool clear=true)
	{

		int hslen=haystack.length();
		int needlen=needle.length();
		const char *shaystack=haystack.c_str();
		char *hay=new char[hslen+1];
		strcpy(hay,shaystack);
		const char *sneedle=needle.c_str();

		if(clear)
			vec.clear();

		int start=0;


		for(int i=0;i<hslen;i++)
		{
			char c=hay[i];
			bool found=false;
			for(int j=0;j<needlen;j++)
			{
				if (sneedle[j]==c)
					found=true;
			}
			if(found)
			{
				hay[i]='\0';
				//const char*toPush=hay+start;
				if (i>start)
					vec.push_back((hay+start));

				hay[i]=' ';
				start=i+1;
			}
		}

		if(start<hslen)
		{
			vec.push_back( (shaystack+start));
		}else if(start==hslen) //new sep 14 2009
		{
			vec.push_back("");
		}


		delete[] hay;
		return vec;
	}

	inline static string replace(const string& haystack,const string& needle,const string& replacer,int numReplace=0)
	{
		const char*phay=haystack.c_str();
		int haysize=haystack.length();
		const char*pneedle=needle.c_str();
		int needlesize=needle.length();
		if(needlesize>haysize)
			return haystack;

		int replacersize=replacer.length();

		string returnvalue=haystack;

		//const char*phayoutstart=phay;
		int i=0;
		int numReplaced=0;
		int ri=0;

		for(i=0;i<haysize-needlesize+1;i++,ri++)
		{
			if(!strncmp(phay+i,pneedle,needlesize)) //match!
			{
				//output until phay

				//cerr<<"match at "<<i;

				//returnvalue.erase(ri,needlesize);
				//cerr<<returnvalue<<endl;
				//returnvalue.insert(ri,replacer);
				returnvalue.replace(ri,needlesize,replacer);
				//cerr<<returnvalue<<endl;
				i+=needlesize-1;
				ri+=replacersize-1;

				//cerr<<" replaced as "<<returnvalue<<endl;
				numReplaced++;
				if(numReplace>0 && numReplaced==numReplace)
					return returnvalue;
			}
		}

		return returnvalue;
	}
	inline static string replaceFromEnd(const string& haystack,const string& needle,const string& replacer,int numReplace)
	{
		const char*phay=haystack.c_str();
		int haysize=haystack.length();
		const char*pneedle=needle.c_str();
		int needlesize=needle.length();
		if(needlesize>haysize)
			return haystack;

		string returnvalue=haystack;

		//const char*phayoutstart=phay;
		int i=0;
		int numReplaced=0;

		for(i=haysize-needlesize+1;i>=0;i--)
		{
			if(!strncmp(phay+i,pneedle,needlesize)) //match!
			{
				//output until phay

				//cerr<<"match at "<<i;
				//returnvalue.erase(i,needlesize);
				//returnvalue.insert(i,replacer);
				returnvalue.replace(i,needlesize,replacer);
				//cerr<<" replaced as "<<returnvalue<<endl;
				numReplaced++;
				if(numReplace>0 && numReplaced==numReplace)
					return returnvalue;
			}
		}

		return returnvalue;
	}

	inline static void shuffleStringInPlace( string& str)
	{
		char tmp;
		int len=str.length();
		for(int i=0;i<len;i++)
		{
			int toSwap=rand()%len;
			tmp=str[i];
			str[i]=str[toSwap];
			str[toSwap]=tmp;
		}
	}
	inline static string shuffleString(const string&str)
	{
		string inter(str);
		shuffleStringInPlace(inter);
		return inter;
	}

};

#endif /*STRINGUTIL_H_*/
