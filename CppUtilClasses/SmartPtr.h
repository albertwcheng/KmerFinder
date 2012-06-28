#ifndef SMARTPTR_H_
#define SMARTPTR_H_

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
 
 This is a class which wraps around pointer and allow comparison of dereferenced object pointed by the pointed instead of pointer address comparisons
 e.g.,
 
 int a[]={5,4,3,2};
 SmartPtr<int> a0=&a[0];
 SmartPtr<int> a1=&a[1];
 
 a0>a1  //which is *a0=5, *a1=4 would returns true
 
 
 ******/

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


#endif /*SMARTPTR_H_*/
