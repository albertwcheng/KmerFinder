/*
 * SuperConfigFile.cpp
 *
 *  Created on: Sep 14, 2009
 *      Author: awcheng
 */

#include "SuperConfigFile.h"





set<string> SuperConfigFile::keys() const
{
	set<string> keyset;
	for(ConfigFile::mapci i=this->myContents.begin();i!=this->myContents.end();i++)
		keyset.insert(i->first);

	return keyset;
}

SuperConfigFile::SuperConfigFile(string filename,string delimiter,string comment,string sentry):ConfigFile(filename,delimiter,comment,sentry), varDelim("$"), varMarker('@')
{
	loadExterns();
}

void SuperConfigFile::justLoad(string filename)
{
	std::ifstream fin(filename.c_str());
	if(!fin.good())
		throw ConfigFile::file_not_found(filename);
	this->justLoad(fin);
	fin.close();
}

SuperConfigFile::SuperConfigFile():varDelim("$"), varMarker('@')
{

}

SuperConfigFile::SuperConfigFile(int argc,const char**argv,string arr_delim): varDelim("$"), varMarker('@')
{
	bool acceptValue=false;
	string key="";
	vector<string> general_argvs;
	for(int i=0;i<argc;i++)
	{
		if(!acceptValue)
		{
			if(argv[i][0]=='-')
			{
				acceptValue=true;
				const char*curarg=argv[i];
				key=string(curarg+1);
			}
			else
			{
				general_argvs.push_back(argv[i]);
			}
		}
		else
		{
			this->add(key,argv[i]);
			acceptValue=false;
		}

	}

	this->add("@argv",StringUtil::formArrayString(general_argvs,arr_delim));
	this->loadExterns();
}

istream& SuperConfigFile::justLoad(istream& is)
{
	ConfigFile& baseObj=*this;
	is>>baseObj;
	return is;
}


istream& operator>>(istream& is,SuperConfigFile& scf)
{
	scf.justLoad(is);
	scf.loadExterns();
	return is;
}
ostream& operator<<(ostream& os,const SuperConfigFile& scf)
{
	const ConfigFile& baseObj=scf;
	os<<baseObj;
	return os;
}

ostream& SuperConfigFile::listAllItemsAsSuperString(ostream& os,bool includeExterns) const
{
	if(includeExterns)
	{
		os<<"externs:"<<endl;
		for(vector<string>::const_iterator i=this->externRecord.begin();i!=this->externRecord.end();i++)
		{
			os<<"\t"<<(*i)<<endl;
		}


	}
	os<<"items:"<<endl;
	for(ConfigFile::mapci i=this->myContents.begin();i!=this->myContents.end();i++)
	{

		os<<"\t"<<i->first<<"="<<this->readSuperString(i->first)<<endl;
	}

	return os;
}

void SuperConfigFile::loadExterns()
{
	if(this->keyExists("@externs"))
	{
		//whether to initiate the tree
		//load external stuff
		vector<string> curExtern=this->readArray("@externs");
		this->remove("@externs");
		for(vector<string>::iterator i=curExtern.begin();i!=curExtern.end();i++)
			this->externQueue.push_back(*i);

		while(this->externQueue.size()>0)
		{

			string u=this->externQueue.front();

			this->externRecord.push_back(u);

			this->externQueue.pop_front();
			this->justLoad(u);
			if(this->keyExists("@externs"))
			{
				curExtern=this->readArray("@externs");
				for(vector<string>::iterator i=curExtern.begin();i!=curExtern.end();i++)
					this->externQueue.push_back(*i);
				this->remove("@externs");
			}

		}
	}
}

void SuperConfigFile::readArray(string name,vector<string>& resultStrings,bool clear,string _delim) const
{

	string preString=this->read<string>(name);
	vector<string> splits;
	StringUtil::split(preString,_delim,splits);

	if(clear)
	{
		resultStrings.clear();
	}

	for(vector<string>::iterator i=splits.begin();i!=splits.end();i++)
	{
		resultStrings.push_back(this->parseSuperString(*i));
	}


}

bool SuperConfigFile::readBool(string _name) const
{
	string preString=StringUtil::toUpper(this->readSuperString(_name));
	if(preString=="1" || preString=="TRUE" || preString=="YES" || preString=="Y" || preString=="T")
	{
		return true;
	}
	else
	{
		return false;
	}

}

vector<string> SuperConfigFile::readArray(string _name, string delim) const
{
	vector<string> vstring;
	this->readArray(_name,vstring,true,delim);
	return vstring;
}

void SuperConfigFile::setVarDelim(string _varDelim)
{
	this->varDelim=_varDelim;
}

string SuperConfigFile::getVarDelim() const
{
	return this->varDelim;
}

void SuperConfigFile::setVarMarker(char _varMarker)
{
	this->varMarker=_varMarker;
}

char SuperConfigFile::getVarMarker() const
{
	return this->varMarker;
}

string SuperConfigFile::escape(const string& str) const
{
	return StringUtil::escape(str,"ntz@","\n\t#@",'\\');
}

string SuperConfigFile::readSuperString(string _name,string _default) const
{
	if(!this->keyExists(_name))
	{
		return _default;
	}
	else
	{
		return this->readSuperString(_name);
	}
}

string SuperConfigFile::parseSuperString(string preString) const
{
	vector<string> splits;
	StringUtil::split(preString,this->varDelim,splits);

	int nSplit=splits.size();
	string resultString=this->escape(splits[0]);
	for(int i=1;i<nSplit-1;i++)
	{
		//cerr<<"split #"<<(i+1)<<splits[i]<<endl;
		string &curSplit=splits[i];
		if(curSplit.length()==0)
		{
			resultString+=this->varDelim;
		}
		else if(curSplit.length()>1 && curSplit[0]==this->varMarker)
		{
			string ref=curSplit.substr(1);
			resultString+=this->readSuperString(ref);
		}
		else
		{
			resultString+=this->escape(curSplit);
		}
	}

	if(nSplit>1)
	{
		resultString+=splits[nSplit-1];
	}

	return resultString;
}

string SuperConfigFile::readSuperString(string _name) const
{
	string preString=this->read<string>(_name);

	return this->parseSuperString(preString);
}
