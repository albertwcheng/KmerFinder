#include <iostream>
#include <set>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "SortPair.h"
using namespace std;




class TrySortableKeyPairExt:public SortPair<SortPair<int,int>, int>
{
public:
	
	syn_rwaccessors(start0,k1.k1,int)
	syn_rwaccessors(end1,k1.k2,int)
	syn_rwaccessors(iid,k2,int)
	
	TrySortableKeyPairExt(int _start0,int _end1,int _iid){
		set_start0(_start0);
		set_end1(_end1);
		set_iid(_iid);
	}
};

	 
int main (int argc,char **argv)
{
	set<TrySortableKeyPairExt> tryset;
	tryset.insert(TrySortableKeyPairExt(15,23,24));
	tryset.insert(TrySortableKeyPairExt(6,12,34));
	tryset.insert(TrySortableKeyPairExt(15,23,22));
	
	for(set<TrySortableKeyPairExt>::iterator i=tryset.begin();i!=tryset.end();i++){
		cout<<i->start0()<<","<<i->end1()<<":"<<i->iid()<<endl;
	}
}