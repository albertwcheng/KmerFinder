#include <queue>
#include <stack>
#include <string>
#include <iostream>
#include <inttypes.h>
#include <list>
using namespace std;

typedef uint64_t uint64;

#define MODE_FOREGROUND 'F'
#define MODE_BACKGROUND 'B'

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



class DNATreeNode{
	public:
		DNATreeNode* parent;
		
		char thisLetter;
		DNATreeNode* A_branch;
		DNATreeNode* C_branch;
		DNATreeNode* G_branch;
		DNATreeNode* T_branch;
		
		
		
		
		uint64 fgcount;
		uint64 bgcount;
		
		double enrichment(double normalization=1.0) const{
			if(bgcount<1){
				return -1;
			}	
			return double(fgcount)/bgcount*normalization;	
		}
		
		inline bool operator < (const DNATreeNode& right) const{
			return this->enrichment()<right.enrichment();
		}
		
		DNATreeNode(DNATreeNode* _parent=NULL,char _thisLetter='.'):parent(_parent),thisLetter(_thisLetter),A_branch(NULL),C_branch(NULL),G_branch(NULL),T_branch(NULL),fgcount(0),bgcount(0){}
		
		DNATreeNode* selectBranch(char c,bool create=false){
			switch( c ){
				case 'A':case 'a':
					if(!A_branch && create){
						A_branch=new DNATreeNode(this,'A');	
					}
					return A_branch;
				case 'C':case 'c':
					if(!C_branch  && create ){
						C_branch=new DNATreeNode(this,'C');	
					}
					return C_branch;	
				case 'G':case 'g':
					if(!G_branch && create){
						G_branch=new DNATreeNode(this,'G');	
					}
					return G_branch;
				case 'T':case 't':
					if(!T_branch && create){
						T_branch=new DNATreeNode(this,'T');
					}
					return T_branch;
				default:
					return NULL;	
			}	
			
			
	
		}
		
		
		DNATreeNode* visitInc(char c,char mode){
			
			DNATreeNode* selectedBranch=selectBranch(c,true);
				
			if(selectedBranch){
				if(mode==MODE_FOREGROUND){
					selectedBranch->fgcount++;	
				}else if(mode==MODE_BACKGROUND){
					selectedBranch->bgcount++;
				}	
			}
			
			return selectedBranch;
		}
		
		
		
		
		DNATreeNode* getNodeByPath(const char*path,int length){
			DNATreeNode* curNode=this;
			for(int i=0;i<length;i++){
				
				curNode=curNode->selectBranch( path[i] );
				if(curNode==NULL){
					break;	
				}
			}	
			
			return curNode;
		}
		
		string getPathFromRoot(){
			DNATreeNode* curNode=this;
			string path;
			while(curNode){
				if(curNode->thisLetter!='.')
					path=curNode->thisLetter+path;
				curNode=curNode->parent;	
			}
			
			return path;	
		}
		
		void putChildrenOntoQ(queue<DNATreeNode*> &Q){
			if (A_branch) Q.push(A_branch);
			if (C_branch) Q.push(C_branch);
			if (G_branch) Q.push(G_branch);
			if (T_branch) Q.push(T_branch);
		}
		
		void putChildrenOntoS(stack<DNATreeNode*> &S){
						
			
			if (T_branch) S.push(T_branch);
			if (G_branch) S.push(G_branch);
			if (C_branch) S.push(C_branch);
			if (A_branch) S.push(A_branch);

		}		
		void freeSubtree(){
			//DFS delete tree
			queue<DNATreeNode*> Q;
			putChildrenOntoQ(Q);
			while(!Q.empty()){
				DNATreeNode* v=Q.front();
				Q.pop();
				v->putChildrenOntoQ(Q);
				//now delete v
				delete v;	
			}
		}
		
		void printSubtree(){
			stack<DNATreeNode*> S;
			putChildrenOntoS(S);
			while(!S.empty()){
				DNATreeNode* v=S.top();
				S.pop();
				cout<<v->getPathFromRoot()<<"\t"<<v->fgcount<<"\t"<<v->bgcount<<endl;
				v->putChildrenOntoS(S);	
			}
		}
	
};


class DNATree{
	public:
		DNATreeNode *root;
		DNATreeNode **tails; //keep the tails
		int maxK;
		
		char countMode;
		uint64 numNodes;
		
		list< SmartPtr<DNATreeNode> > *listsPerK;
		
		void setCountMode(char _countMode){
			countMode=_countMode;	
		}
		
		void addNodeToList(int depth,DNATreeNode* node){
			listsPerK[depth-1].push_back(node);
			numNodes++;
		}
		
		
		void sortLists(){
			for(int i=0;i<maxK;i++){
				listsPerK[i].sort();
			}	
		}
		
		void printLevelNodes(ostream& os,int level,double normalization){
			for(list<SmartPtr<DNATreeNode> >::reverse_iterator i=listsPerK[level-1].rbegin();i!=listsPerK[level-1].rend();i++){
				os<<(*i)->getPathFromRoot()<<"\t"<<(*i)->fgcount<<"\t"<<(*i)->bgcount<<"\t"<<(*i)->enrichment(normalization)<<endl;	
			}	
		}
		
		
		DNATree(int _maxK):maxK(_maxK),numNodes(0){
			root=new DNATreeNode;
			tails=new DNATreeNode*[maxK];
			//init tails
			for(int i=0;i<maxK;i++){
				tails[i]=NULL;	
			}
			countMode=MODE_FOREGROUND;
			listsPerK=new list<SmartPtr<DNATreeNode> >[maxK];
		}
		~DNATree(){
			root->freeSubtree();
			delete root;
			delete[] tails;
			delete[] listsPerK;
			
		}
		
		void feedSafe(char c){
			int depth;
			for(int i=1;i<maxK;i++){
				//shift
				if(!tails[i])
					continue;
				
				tails[i-1]=tails[i]->visitInc( c, countMode);
				depth=maxK-i+1;
				if(tails[i-1]->fgcount+tails[i-1]->bgcount==1){
					//new node
					addNodeToList(depth,tails[i-1]);
					
				}
				
			}
			
			tails[maxK-1]=root->visitInc( c,countMode);
			depth=1;
			if(tails[maxK-1]->fgcount+tails[maxK-1]->bgcount==1){
				addNodeToList(depth,tails[maxK-1]);
				
			}
			
		}
		
		void feed(char c){
			switch( c){
				case 'A':case 'a':case 'C':case 'c':case 'G':case 'g':case 'T':case 't':
					//ok!
					feedSafe( c);
				break;
				
				default: //anything not matching, ignore and set all tails NULL;
					for(int i=0;i<maxK;i++){
						tails[i]=NULL;	
					}
			}
		}
		void feedSeq(const string& S){
			const char* s=S.c_str();
			int len=S.length();
			for(int i=0;i<len;i++){
				feed(s[i]);	
			}
		}
		
		DNATreeNode* getNodeByPath(const string& path){
			return root->getNodeByPath(path.c_str(),path.length());	
		}
		
		
		void printTree(){
			root->printSubtree();
		}
};