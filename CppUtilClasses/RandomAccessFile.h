/***************************************************************************
 Copyright 2010 Wu Albert Cheng <albertwcheng@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 *******************************************************************************/ 

#ifndef RANDOMACCESSFILE_H_
#define RANDOMACCESSFILE_H_
#include <ostream>
#include <string>
#include <fstream>

#include <string>
using namespace std;

#define READLENGTH 1024
#define READLENGTHP	1025
#ifndef MIN
#define MIN(x,y) (((x)<(y))?(x):(y))
#endif
class RandomAccessFile: public ifstream
{
public:
	char *rafbuf;
	int rafbuflen;
	bool externBuffer;
	//string get(int start,int endEx);
	void transfer(ostream&os,int start,int endEx);
	string get(int start,int endEx);
	RandomAccessFile(string filename,int bufsize=READLENGTHP, char *_rafbuf=NULL);
	~RandomAccessFile();
};


#endif /*RANDOMACCESSFILE_H_*/
