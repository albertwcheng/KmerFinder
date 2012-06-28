/*
 * SuperConfigFile.h
 *
 *  Created on: Sep 14, 2009
 *      Author: awcheng
 */

/*
 * superString=$@varname$@another var$ followed by string and escape\nthis is another line
 * money=$$315.00 and email address after variable $@username$\@gmail.com
 * array=albert|grace|amy|can also be a superString with $@varname$
 * escapes \n \t \@ \z=#
 *
 *
 */




#ifndef SUPERCONFIGFILE_H_
#define SUPERCONFIGFILE_H_

#include "ConfigFile.h"
#include "CppUtilClasses/StringUtil.h"
#include <list>
#include <set>
using namespace std;

class SuperConfigFile;

istream& operator>>(istream& is,SuperConfigFile& scf);
ostream& operator<<(ostream& os,const SuperConfigFile& scf);


class SuperConfigFile: public ConfigFile
{
private:
	string varDelim;
	char varMarker;
public:
	void setVarDelim(string _varDelim);
	string getVarDelim() const;
	void setVarMarker(char _varMarker);
	char getVarMarker() const;
	list<string> externQueue;
	SuperConfigFile(string filename,
			string delimiter = "=",
            string comment = "#",
			string sentry = "EndConfigFile");
	SuperConfigFile(int argc,const char**argv,string arr_delim="|");

	SuperConfigFile();

	string readSuperString(string _name) const;
	void readArray(string _name,vector<string>& into,bool clear=true,string _delim="|") const;
	vector<string> readArray(string _name,string delim="|") const;
	string parseSuperString(string _preString) const;
	string escape(const string& str) const;
	string readSuperString(string _name,string _default) const;
	bool readBool(string _name) const;
	//void load(string filename);
	void justLoad(string filename);
	istream& justLoad(istream& is);
	void loadExterns();

	set<string> keys() const;

	vector<string> externRecord;

	ostream& listAllItemsAsSuperString(ostream& os,bool includeExterns=true) const;

	friend istream& operator>>(istream& is,SuperConfigFile& scf);
	friend ostream& operator<<(ostream& os,const SuperConfigFile& scf);



};

#endif /* SUPERCONFIGFILE_H_ */
