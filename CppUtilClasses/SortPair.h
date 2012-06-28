#ifndef _SORT_PAIR_H
#define _SORT_PAIR_H

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



template<class T1,class T2>
class SortPair
{
public:
	T1 k1;
	T2 k2;
	inline SortPair(){}
	inline SortPair(T1 _k1,T2 _k2): k1(_k1), k2(_k2) {
	}
	
	inline bool operator <(const SortPair<T1,T2>& obj) const
	{
		if(this->k1==obj.k1)
		{
			return (this->k2<obj.k2);
		}
		
		return this->k1<obj.k1;
		
	}
	
	inline bool operator !=(const SortPair<T1,T2>& obj) const
	{
		return !((*this)==obj);
	}
	
	inline bool operator >=(const SortPair<T1,T2>& obj) const
	{
		return !((*this)<obj);
		
	}
	
	inline bool operator <=(const SortPair<T1,T2>& obj) const
	{
		return !((*this)>obj);
		
	}
	
	bool operator ==(const SortPair<T1,T2>& obj) const
	{
		return this->k1==obj.k1 && this->k2==obj.k2;
	}
	
	bool operator >(const SortPair<T1,T2>& obj) const
	{
		if(this->k1==obj.k1)
		{
			return (this->k2>obj.k2);
		}
		
		return this->k1>obj.k1;
		
	}
	
	inline SortPair<T1,T2>& operator+=(const SortPair<T1,T2>& right){
		this->k1+=right.k1;
		this->k2+=right.k2;
		return *this;
	}
};

#define syn_rwaccessors(NAME,VAR,TYPE) inline void set_ ## NAME(const TYPE& _ ## NAME) \
{ VAR= _ ## NAME; } \
inline TYPE NAME() const \
{ return VAR; }


/* example in SortPairExample.cpp */


#endif