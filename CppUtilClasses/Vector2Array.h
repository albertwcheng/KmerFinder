/*
 * Vector2Array.h
 *
 *  Created on: Sep 14, 2009
 *      Author: awcheng
 */

#ifndef VECTOR2ARRAY_H_
#define VECTOR2ARRAY_H_

#include <vector>

template<class T>
class Vector2Array
{
private:
	T* arr;
	int _size;
public:

	int size() const
	{
		return _size;
	}
	Vector2Array(const vector<T>& V)
	{
		this->_size=V.size();
		arr=new T[this->_size];
		for(int i=0;i<this->_size;i++)
		{
			arr[i]=V[i];
		}
	}
	~Vector2Array()
	{
		delete[] arr;
	}
	operator T*()
	{
		return arr;
	}

	operator const T*() const
	{
		return arr;
	}
};

#endif /* VECTOR2ARRAY_H_ */
