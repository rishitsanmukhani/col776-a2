#pragma once
#include "clique.h"
#include "bethe.h"

class Word{
public:
	char w1[MAX_WORD_LEN+1], w2[MAX_WORD_LEN+1];
	vector<int> w1_i, w2_i;
	vector<int> pw1_i, pw2_i;
	vector<int> img1, img2;
	Graph graph;
	CliqueTree ctree;
	BetheGraph bgraph;
	int n1, n2, n;
	int nc1, nc2;
	double ll, T;
	void init(){
		n1 = strlen(w1);img1.resize(n1);w1_i.resize(n1);pw1_i.resize(n1);
		n2 = strlen(w2);img2.resize(n2);w2_i.resize(n2);pw2_i.resize(n1);
		n = n1 + n2;
	}
	void convert_to_int(map<char,int>& char_to_id){
		for(int i=0;i<n1;i++)
			w1_i[i] = char_to_id[w1[i]];
		for(int i=0;i<n2;i++)
			w2_i[i] = char_to_id[w2[i]];
	}
	void init_graph(const Model& model, int type){
		graph.init(img1, img2, model, type);
		if(MP)
			init_ctree();
		else
			init_bgraph();
	}
	void init_accuracy(){
		nc1=0;
		for(int i=0;i<n1;i++)
			if(pw1_i[i]==w1_i[i])
				nc1++;
		nc2=0;
		for(int i=0;i<n2;i++)
			if(pw2_i[i]==w2_i[i])
				nc2++;
	}
	int max_marginal(const vector<double>& v){
		int ret=0;
		for(int i=1;i<v.size();i++)
			if(v[i]>v[ret])
				ret = i;
		return ret;
	}
	void init_ctree(){
		ctree.init(graph);
		// log-likelihood
		ll=0;
		for(int i=0;i<n1;i++){
			vector<double>& marginal = ctree.nodes[i].marginal.v;
			ll += log(marginal[w1_i[i]]);
			pw1_i[i] = max_marginal(marginal);
		}
		for(int i=0;i<n2;i++){
			vector<double>& marginal = ctree.nodes[i+n1].marginal.v;
			ll += log(marginal[w2_i[i]]);
			pw2_i[i] = max_marginal(marginal);
		}
		T = ctree.time_taken;
		init_accuracy();
	}
	void init_bgraph(){
		bgraph.init(graph);
		// log-likelihood
		ll=0;
		for(int i=0;i<n1;i++){
			vector<double>& marginal = bgraph.var_nodes[i].marginal.v;
			ll += log(marginal[w1_i[i]]);
			pw1_i[i] = max_marginal(marginal);
		}
		for(int i=0;i<n2;i++){
			vector<double>& marginal = bgraph.var_nodes[i+n1].marginal.v;
			ll += log(marginal[w2_i[i]]);
			pw2_i[i] = max_marginal(marginal);
		}
		T = bgraph.time_taken;
		init_accuracy();
	}
	void print(){
		printf("%s %s\n", w1, w2);
		graph.print();
		if(MP) ctree.print();
		else bgraph.print();
		printf("Word-1: %f (%d/%d)\n", double(nc1)/n1, nc1, n1);
		printf("Word-2: %f (%d/%d)\n", double(nc2)/n2, nc2, n2);
		printf("LL    : %f\n", ll/2);
		printf("Time  : %f ms\n", T);
	}
};