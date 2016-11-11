#pragma once
#include "model.h"
#include "factor.h"
using namespace std;

struct Graph{
	// G[i][j] = bitmask(3) -> [b2 b1 b0]
	// b0 = transition edge
	// b1 = skip edge
	// b2 = pair skip edge
	vector<vector<int> > G;
	int n;
	// key = bitmask(n) -> [bn-1 bn-2 ... b1 b0]
	// (bi==1) => ith node is present
	map<int, Factor> factors;

	vector<int> ordering;
	void init(const vector<int>& img1, const vector<int>& img2, const Model& model, int type){
		init_graph(img1, img2, type);
		vector<int> img;
		img.insert(img.end(), img1.begin(), img1.end());
		img.insert(img.end(), img2.begin(), img2.end());
		init_factors(img, model, type);
		init_ordering();
	}
	void init_graph(const vector<int>& img1, const vector<int>& img2, int type){
		int n1 = img1.size(), n2 = img2.size();
		n = img1.size() + img2.size();
		// Generating graph
		G.assign(n, vector<int>(n,0));
		if(type>=2){
			// Transition Factors
			for(int i=0;i<n1-1;i++)SET(G[i][i+1],0),SET(G[i+1][i],0);
			for(int i=0;i<n2-1;i++)SET(G[n1+i][n1+i+1],0),SET(G[n1+i+1][n1+i],0);
		}
		if(type>=3){
			// Skip Factors
			for(int i=0;i<n1;i++)
				for(int j=i+1;j<n1;j++)
					if(img1[i]==img1[j])SET(G[i][j],1), SET(G[j][i],1);
			for(int i=0;i<n2;i++)
				for(int j=i+1;j<n2;j++)
					if(img2[i]==img2[j])SET(G[n1+i][n1+j],1), SET(G[n1+j][n1+i],1);
		}
		if(type>=4){
			for(int i=0;i<n1;i++)
				for(int j=0;j<n2;j++){
					if(img1[i]==img2[j])SET(G[i][n1+j],2), SET(G[n1+j][i],2);
				}
		}		
	}
	void init_factors(const vector<int>& img, const Model& model, int type){
		// phi(i) -> corresponding to OCR factors
		Factor phi;
		for(int i=0;i<n;i++){
			phi.init((1<<i));
			for(int c=0;c<10;c++){
				phi.v[c] *= model.ocr[img[i]][c];
			}
			factors[1<<i] = phi;
		}
		// phi(i,j)
		for(int i=0;i<n;i++){
			for(int j=i+1;j<n;j++){
				if(G[i][j]==0)continue;
				phi.init((1<<i)+(1<<j));
				if(BIT(G[i][j], 0)){														// Transition Factor (note importance of ordering)
					for(int c1=0;c1<10;c1++)											// c1 -> assigned to i
						for(int c2=0;c2<10;c2++)										// c2 -> assigned to j = i+1
							phi.v[10*c2+c1] *= model.trans[c1][c2];
				}
				if(BIT(G[i][j], 1)){													// Skip Factor
					for(int c=0;c<10;c++)
						phi.v[10*c + c] *= model.skip;
				}
				if(BIT(G[i][j], 2)){													// Pair-Skip Factor
					for(int c=0;c<10;c++)
						phi.v[10*c + c] *= model.pair_skip;
				}
				factors[(1<<i) + (1<<j)] = phi;
			}
		}
	}
	vector<int> get_nbrs(const vector<vector<int> >& G, int k){
		vector<int> ret;
		for(int i=0;i<G.size();i++){
			if(G[k][i]!=0)
				ret.push_back(i);
		}
		return ret;
	}
	int min_fill_heuristics(const vector<vector<int> >& G,const vector<int>& flag){
		int min_val = n*n, idx = -1;
		for(int i=0;i<n;i++){
			if(flag[i]==1)continue;
			vector<int> nbrs = get_nbrs(G, i);
			int val=0;
			for(int a=0;a<nbrs.size();a++){
				for(int b=a+1;b<nbrs.size();b++){
					if(G[nbrs[a]][nbrs[b]]==0)val++;
				}
			}
			if(val<min_val)min_val = val, idx=i;
		}
		return idx;
	}
	void init_ordering(){
		vector<vector<int> > induced_G = G;
		vector<int> flag(n, 0);
		ordering.resize(n);
		for(int k=0;k<n;k++){
			int i = min_fill_heuristics(induced_G, flag);
			flag[i]=1;
			vector<int> nbrs = get_nbrs(induced_G, i);
			for(int a=0;a<nbrs.size();a++){
				induced_G[nbrs[a]][i]=0;
				induced_G[i][nbrs[a]]=0;
				for(int b=a+1;b<nbrs.size();b++){
					induced_G[nbrs[a]][nbrs[b]]=1;
					induced_G[nbrs[b]][nbrs[a]]=1;					
				}
			}
			ordering[k] = i;
		}
	}
	void print(){
		printf("Graph:\n");
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++)
				cout<<G[i][j]<<" ";
			cout<<endl;
		}
		printf("Elimination Ordering:\n");
		for(auto o:ordering)cout<<o<<" ";
		cout<<endl;
	}
};