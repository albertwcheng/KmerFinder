#ifndef __MATRIX_H
#define __MATRIX_H

#include <string.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <stdlib.h>

using namespace std;

template<typename T>
class Matrix{
	private:
		vector< vector<T> > dmat;
		T** smat;
		int smatRows;
		int smatCols;
	public:
		inline Matrix():smat(NULL),smatRows(0),smatCols(0){
			
		}
		inline Matrix(const Matrix<T> &right):smat(NULL),smatRows(right.smatRows),smatCols(right.smatCols){
			if(right.smat){
				this->smat=new T*[right.smatRows];
				for(int i=0;i<right.smatRows;i++){
					T* Trow=new T[right.smatCols];
					this->smat[i]=Trow;
					T* TrowRight=right.smat[i];
					for(int j=0;j<right.smatCols;j++){
						Trow[j]=TrowRight[i];	
					}
				}	
			}else{
				this->dmat=right.dmat;
			}
		}
		
		inline bool addRow(const vector<T>& v){
			if(this->smat){
				return false; //finalized already, can't add	
			}
			
			if(v.size()<1){
				return false; //don't add empty row	
			}
			if(this->dmat.size()>1 && this->dmat[0].size()!=v.size()){
				return false; //can't add this inconsistent		
			}	
		
			this->dmat.push_back(v);
			
			
			return true;
		}
		
		inline ~Matrix(){
			if(smat){
				for(int i=0;i<this->smatRows;i++){
					delete[] this->smat[i];
				}
				delete[] this->smat;
			}
		}
		inline void finalize(){
			if(this->smat || this->dmat.size()<1){
				//already finalized or empty dmat
				return;
			}
			
			this->smatRows=this->dmat.size();
			this->smatCols=this->dmat[0].size();
			
			this->smat=new T*[this->smatRows];
			for(int i=0;i<this->smatRows;i++){
				T* Trow=new T[this->smatCols];
				this->smat[i]=Trow;
				for(int j=0;j<this->smatCols;j++){
					Trow[j]=this->dmat[i][j];	
				}
			}
			
			this->dmat.clear();		
		}
		
		inline T& itemAt(int i,int j){
			if(smat){
				return this->smat[i][j];	
			}
			else{
				return this->dmat[i][j];
			}
		}
		
		inline const T& itemAt(int i,int j) const{
			if(smat){
				return this->smat[i][j];	
			}
			else{
				return this->dmat[i][j];
			}
		}		
		
		inline int getNumberOfRows() const{
			if(this->smat){
				return this->smatRows; 	
			}else{
				return this->dmat.size();	
			}
		
		}
		
		inline int getNumberOfCols() const{
			if(this->smat){
				return this->smatCols;	
			}else{
				if(this->dmat.size()<1){
					return 0;	
				}
				else{
					return this->dmat[0].size();	
				}
			}
		}
		
		inline vector<string>& split(const string haystack,const string& needle,vector<string>& vec,bool clear=true) const
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
		
		inline void print(ostream& os){
			int nrows=this->getNumberOfRows();
			int ncols=this->getNumberOfCols();
			
			if(nrows<1 || ncols<1){
				return;	
			}
			for(int i=0;i<nrows;i++){
				os<<this->itemAt(i,0);
				for(int j=1;j<ncols;j++){
					os<<"\t"<<this->itemAt(i,j);	
				}
				os<<endl;
			}
		}
		
		//a function class to convert string to class T
		class StringConverter
		{
			public:
			virtual T operator() (const string& str) const=0;//{return T();}
		};
		
		inline void readFile(const string&filename,const StringConverter& converter,const string&sep="\t")
		{
			
			if(this->smat){
				return; //already finalized, can't do anything	
			}
			
			string buffer;
			fstream fil(filename.c_str());
			while(fil.good()){
				buffer="";
				getline(fil,buffer);
				if(buffer.length()<1){
					continue;
				}
				//cerr<<"read buffer "<<buffer<<endl;
				if(buffer.length()<1 || buffer[0]=='#'){
					continue;
				}
					
				vector<string> splits;
				vector<T> Ts;
				split(buffer,sep,splits);
				//cerr<<"splits size="<<splits.size()<<endl;
				for(vector<string>::iterator i=splits.begin();i!=splits.end();i++){
						
					Ts.push_back(converter(*i));	
				}
				this->addRow(Ts);
				
			}
			fil.close();
		}
		
};


class StringToFloat:public Matrix<double>::StringConverter
{
	double operator ()(const string& str) const
	{
		//cerr<<str<<endl;
		return atof(str.c_str());
	}
};

class StringToString:public Matrix<string>::StringConverter
{
	string operator ()(const string& str) const
	{
		//cerr<<str<<endl;
		return str;
	}
};

class StringToInt:public Matrix<int>::StringConverter
{
	int operator ()(const string& str) const
	{
		//cerr<<str<<endl;
		return atoi(str.c_str());
	}
};


#endif