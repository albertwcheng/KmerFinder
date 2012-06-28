#include "Array.h"

int main(int argc,char**argv){
	/*NDimensionVector<int> M(4,15,3,4,5);
	int acc[]={2,0,3,2};
	M[acc]=5;
	cerr<<M[acc]<<endl;*/
	
	Array<int> a;
	a.push_back(12);
	a.push_back(3);
	Array<int> b;
	b.push_back(15);
	b.push_back(16);
	Array<int> c=a+b*a;
	cerr<<c<<endl;
	Array< Array<int> > _2D;
	_2D.push_back(a);
	_2D.push_back(b);
	cerr<<_2D<<endl;
	Array< Array<int> > _2D2;
	_2D2.push_back(b);
	_2D2.push_back(a);
	cerr<<_2D2<<endl;
	cerr<<(_2D+_2D2)<<endl;
	cerr<<(_2D+a)<<endl;
	
	//test independence
	Array< Array<int> > save2D=_2D;
	cerr<<(--_2D)<<endl;
	cerr<<(save2D)<<endl;
	save2D.append(a);
	cerr<<getDimensions(save2D)<<endl;
	cerr<<getDimensions(a)<<endl;
	Array<Array<Array<int > > > threeD;
	threeD.append(_2D);
	threeD.append(_2D2);
	cerr<<threeD<<endl;
	cerr<<getDimensions(threeD)<<endl;
}