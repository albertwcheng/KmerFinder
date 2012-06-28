#ifndef ARRAY_H_
#define ARRAY_H_

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

 This is a class providing array arithmetics, the default container is the vector class from std C++

 ******/

#include <iostream>
#include <vector>
#include <stdarg.h>
#include <typeinfo>
using namespace std;

template<typename _S>
class NDimensionVector:protected vector<_S>
{
	private:
		int dimension;
		int* sizes;
	public:
		NDimensionVector(int _dimension, ... ):dimension(_dimension){
			va_list Numbers;
			va_start(Numbers,_dimension);
			sizes=new int[this->dimension];
			cerr<<"dimension:"<<this->dimension<<endl;
			int capacityNeeded=1;
			for(int i=0;i<this->dimension;i++){
				sizes[i]=va_arg(Numbers,int);
				cerr<<"getting size:"<<sizes[i]<<endl;
				capacityNeeded*=sizes[i];
			}
			
			this->reserve(capacityNeeded);
		}
		_S& operator[](const int* idx ){
			int realIdx=0;
			int mul=1;
			for(int i=dimension-1;i>=0;i--){
				realIdx+=mul*idx[i];
				mul*=sizes[i];
			}	
			return vector<_S>::operator[](realIdx);
		}
		~NDimensionVector(){
			delete[] sizes;	
		}

		
};

template<typename _S=double>
class Array;

template<typename _S>
ostream& operator << (ostream& os,const Array<_S>& _array){
	os<<"[";
	typename Array<_S>::const_iterator i=_array.begin();
	if(i!=_array.end()){
		os<<*i;	
	}
	i++;
	for(;i!=_array.end();i++){
		os<<","<<*i;	
	}
	os<<"]";
	return os;
}



template<typename _S>
class Array: public vector<_S>
{	
	
	public:
		class InconsistentSizesException{
			public:
				int size1;
				int size2;	
				InconsistentSizesException(int _size1,int _size2):size1(_size1),size2(_size2){
				}
					
		};
		
		//_S& operator[](int i){
		//	return vector<_S>::operator[](i);
		//}
		
		Array operator+(const Array& right) const{
			if(this->size()!=right.size()){
				throw InconsistentSizesException(this->size(),right.size());
			}
			
			Array<_S> newArray;
			
			typename vector<_S>::const_iterator right_i=right.begin();
			for(typename vector<_S>::const_iterator this_i=this->begin();this_i!=this->end();this_i++,right_i++){
		
				newArray.push_back(*this_i+*right_i);
			}
			
			return newArray;
		}
		
		Array operator-(const Array& right) const{
			if(this->size()!=right.size()){
				throw InconsistentSizesException(this->size(),right.size());
			}
			
			Array<_S> newArray;
			
			typename vector<_S>::const_iterator right_i=right.begin();
			for(typename vector<_S>::const_iterator this_i=this->begin();this_i!=this->end();this_i++,right_i++){
				newArray.push_back(*this_i-*right_i);
			}
			
			return newArray;
		}	
		
		Array operator/(const Array& right) const{
			if(this->size()!=right.size()){
				throw InconsistentSizesException(this->size(),right.size());
			}
			
			Array<_S> newArray;
			
			typename vector<_S>::const_iterator right_i=right.begin();
			for(typename vector<_S>::const_iterator this_i=this->begin();this_i!=this->end();this_i++,right_i++){
				newArray.push_back(*this_i/(*right_i));
			}
			
			return newArray;
		}			
			
		Array operator*(const Array& right) const{
			if(this->size()!=right.size()){
				throw InconsistentSizesException(this->size(),right.size());
			}
			
			Array<_S> newArray;
			
			typename vector<_S>::const_iterator right_i=right.begin();
			for(typename vector<_S>::const_iterator this_i=this->begin();this_i!=this->end();this_i++,right_i++){
				newArray.push_back((*this_i)*(*right_i));
			}
			
			return newArray;
		}	
		Array& operator+=(const Array& right) {
			if(this->size()!=right.size()){
				throw InconsistentSizesException(this->size(),right.size());
			}

			
			typename vector<_S>::const_iterator right_i=right.begin();
			for(typename vector<_S>::iterator this_i=this->begin();this_i!=this->end();this_i++,right_i++){
		
				(*this_i)+=*right_i;
			}
			
			return (*this);
		}
		
		Array& operator-=(const Array& right) {
			if(this->size()!=right.size()){
				throw InconsistentSizesException(this->size(),right.size());
			}
			
			typename vector<_S>::const_iterator right_i=right.begin();
			for(typename vector<_S>::iterator this_i=this->begin();this_i!=this->end();this_i++,right_i++){
				(*this_i)-=*right_i;
			}
			
			return (*this);
		}	
		
		Array& operator/=(const Array& right) {
			if(this->size()!=right.size()){
				throw InconsistentSizesException(this->size(),right.size());
			}

			
			typename vector<_S>::const_iterator right_i=right.begin();
			for(typename vector<_S>::iterator this_i=this->begin();this_i!=this->end();this_i++,right_i++){
				(*this_i)/=*right_i;
			}
			
			return (*this);
		}			
			
		Array& operator*=(const Array& right) {
			if(this->size()!=right.size()){
				throw InconsistentSizesException(this->size(),right.size());
			}
			
			
			typename vector<_S>::const_iterator right_i=right.begin();
			for(typename vector<_S>::iterator this_i=this->begin();this_i!=this->end();this_i++,right_i++){
				(*this_i)*=(*right_i);
			}
			
			return (*this);
		}	
		
		
		Array& operator+=(const _S& _s) {
			for(typename vector<_S>::iterator this_i=this->begin();this_i!=this->end();this_i++){
				(*this_i)+=_s;
			}
			return (*this);
		}
		
		Array& operator-=(const _S& _s) {
			for(typename vector<_S>::iterator this_i=this->begin();this_i!=this->end();this_i++){
				(*this_i)-=_s;
			}
			return (*this);
		}
		
		Array& operator*=(const _S& _s) {
			for(typename vector<_S>::iterator this_i=this->begin();this_i!=this->end();this_i++){
				(*this_i)*=_s;
			}
			return (*this);
		}		
		Array& operator/=(const _S& _s) {
			for(typename vector<_S>::iterator this_i=this->begin();this_i!=this->end();this_i++){
				(*this_i)/=_s;
			}
			return (*this);
		}			
	
		Array operator+(const _S& _s) const{
			Array<_S> newArray;
			for(typename vector<_S>::const_iterator this_i=this->begin();this_i!=this->end();this_i++){
				newArray.push_back((*this_i)+_s);	
			}
			return newArray;
		}
		Array operator-(const _S& _s) const{
			Array<_S> newArray;
			for(typename vector<_S>::const_iterator this_i=this->begin();this_i!=this->end();this_i++){
				newArray.push_back((*this_i)-_s);	
			}
			return newArray;
		}		
		Array operator*(const _S& _s) const{
			Array<_S> newArray;
			for(typename vector<_S>::const_iterator this_i=this->begin();this_i!=this->end();this_i++){
				newArray.push_back((*this_i)*_s);	
			}
			return newArray;
		}
		Array operator/(const _S& _s) const{
			Array<_S> newArray;
			for(typename vector<_S>::const_iterator this_i=this->begin();this_i!=this->end();this_i++){
				newArray.push_back((*this_i)/_s);	
			}
			return newArray;
		}
		
	    Array& operator++(){
			for(typename vector<_S>::iterator this_i=this->begin();this_i!=this->end();this_i++){
				++(*this_i);
			}
			return *this;
		}
	    Array& operator--(){
			for(typename vector<_S>::iterator this_i=this->begin();this_i!=this->end();this_i++){
				--(*this_i);
			}
			return *this;
		}
		
		Array operator++(int){
			//clone this first
			Array<_S> snapShot(*this);
			++(*this);
			return snapShot;
		}	
		
		Array operator--(int){
			Array<_S> snapShot(*this);
			--(*this);
			return snapShot;	
		}
		

		
		Array& append(const _S& _s){
			this->push_back(_s);
			return *this;	
		}
		

		
		friend ostream& operator << <_S>(ostream& os,const Array& _array);
	
};


template<typename C> 
Array<size_t> getDimensions(const Array< Array< Array<C> > > &_arr )
{
	Array<size_t> dimensions;
	size_t firstDim=_arr.size();
	dimensions.append(firstDim);
	if(firstDim>0){
		int secondDim=_arr[0].size();
		dimensions.append(secondDim);
		if(secondDim>0){
			dimensions.append(_arr[0][0].size());	
		}
	}
	else
		dimensions.append(0);
		
	return dimensions;
}

template<typename C> 
Array<size_t> getDimensions(const Array< Array<C> > &_arr )
{
	Array<size_t> dimensions;
	size_t firstDim=_arr.size();
	dimensions.append(firstDim);
	if(firstDim>0)
		dimensions.append(_arr[0].size());
	else
		dimensions.append(0);
		
	return dimensions;
}

template<typename C> 
Array<size_t> getDimensions(const Array< C > &_arr )
{
	Array<size_t> dimensions;
	size_t firstDim=_arr.size();
	dimensions.append(firstDim);
		
	return dimensions;
}

#endif /*ARRAY_H_*/