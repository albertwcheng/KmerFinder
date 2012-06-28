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

#include <iostream>
#include <fstream>
#include <string>
#include <map>
using namespace std;
#include "StringUtil.h"
#include <sstream>
#include "RandomAccessFile.h"

#include <limits.h>



/*string RandomAccessFile::get(int start, int endEx) {

	if(start<0)
		start=0;

	if(endEx<0)
		endEx=INT_MAX;

	int lenRead=0;
	int lenToRead;

	seekg(start, ios::beg);



	string returna;

	while(start+lenRead<endEx && !eof())
	{

		lenToRead=MIN(buflen,endEx-start-lenRead);

		read(buf, lenToRead);

		buf[lenToRead+1]='\0';
		returna+=buf;

		lenRead+=lenToRead;
	}

	return returna;
}*/

RandomAccessFile::RandomAccessFile(string filename, int bufsize, char *_rafbuf) :
	ifstream(filename.c_str()), rafbuflen(bufsize) {
	if(_rafbuf)
	{
		rafbuf=_rafbuf;
		externBuffer=true;
	}
	else
	{
		externBuffer=false;
		rafbuf=new char[rafbuflen];
	}
	
	fileSize=getFileSize();

}

int RandomAccessFile::getFileSize()
{
seekg (0, ios::end);

int length = tellg();

seekg (0, ios::beg);

return length;

}

RandomAccessFile::~RandomAccessFile()
{
	close();
	if(!externBuffer)
		delete[] rafbuf;
}

string RandomAccessFile::get(int start,int endEx)
{
	ostringstream ss;
	this->transfer(ss,start,endEx);
	return ss.str();
}

void RandomAccessFile::transfer(ostream&os,int start,int endEx)
{
	if(start<0)
		start=0;

	if(endEx<0)
		endEx=fileSize;//INT_MAX; //was INT_MAX
		
	endEx=MIN(endEx,fileSize); //to protect the stream
	

	int lenRead=0;
	int lenToRead;

	seekg(start, ios::beg);


	while(start+lenRead<endEx && !eof())
	{

		lenToRead=MIN(rafbuflen-1,endEx-start-lenRead);

		read(rafbuf, lenToRead);


		int thisRead=gcount();

		rafbuf[thisRead]='\0';

		os<<rafbuf;
		lenRead+=thisRead;
	}

}


