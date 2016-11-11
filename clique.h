#pragma once
#include "graph.h"
using namespace std;

struct CliqueTree{
	struct Node{
		int id;
		Factor psi, beta;
		Factor delta_up, delta_down;
		Factor marginal;
		int up_nbr;
		vector<int> down_nbrs;
		Node(){
			id = -1;
		}
		void init(int _id){
			id = _id;
			up_nbr = -1;
			psi.init(1<<id);
		}
		void init_delta(){
			int m = beta.mask;
			CLR(m, id);
			delta_up.init(m);
			delta_down.init(m);
		}
		void init_psi(map<int, Factor>& factors){
			vector<int> to_erase;
			for(auto& p:factors)
				if(BIT(p.first, id)){
					psi = psi * p.second;
					to_erase.push_back(p.first);
				}
			for(auto& p:to_erase)factors.erase(p);
		}
		void init_marginal(){
			marginal = beta.get_factor(1<<id);
			double sum=0;
			for(auto val:marginal.v)sum += val;
			for(auto& val:marginal.v)val /= sum;			
		}
		void print(){
			cout<<"Node - "<<id<<endl;
			printf("(%d %d %d)\n", up_nbr, psi.mask, delta_up.mask);
			marginal.print();
		}
	};
	double time_taken;
	vector<Node> nodes;
	int n;
	void init(const Graph& G){
		n = G.n;
		nodes.resize(n);
		map<int, Factor> factors = G.factors;
		// Start Variable Elimination
		for(auto x:G.ordering){		// x = variable to be eliminated
			nodes[x].init(x);
			nodes[x].init_psi(factors);
			nodes[x].beta = nodes[x].psi;
			for(int i=0;i<n;i++){
				if(i==x || nodes[i].id==-1)continue;
				if(nodes[i].up_nbr==-1 && BIT(nodes[i].delta_up.mask, x)){
					nodes[i].up_nbr = x;
					nodes[x].down_nbrs.push_back(i);
					nodes[x].beta = nodes[x].beta * nodes[i].delta_up;
				}
			}
			nodes[x].init_delta();
		}
		time_taken = 0;
		Timer t;t.reset();
		UP_pass();
		DOWN_pass();
		time_taken = t.elapsed();
	}
	Factor send_message(int i, int j){
		Factor ret;
		ret = nodes[i].psi;
		if(j != nodes[i].up_nbr)
			ret = ret * nodes[i].delta_down;
		for(auto down_nbr:nodes[i].down_nbrs)
			if(j != down_nbr)
				ret = ret * nodes[down_nbr].delta_up;
		return ret;
	}
	void post_order(int k){
		for(int down_nbr:nodes[k].down_nbrs){
			post_order(down_nbr);
		}
		nodes[k].delta_up = send_message(k, nodes[k].up_nbr).get_factor(nodes[k].delta_up.mask);
	}
	void UP_pass(){
		for(int root=0;root<n;root++){
			if(nodes[root].up_nbr!=-1)continue;
			for(auto down_nbr:nodes[root].down_nbrs){
				post_order(down_nbr);
			}
		}
	}
	void pre_order(int k){
		nodes[k].beta = nodes[k].beta * nodes[k].delta_down;
		for(int down_nbr:nodes[k].down_nbrs){
			nodes[k].beta = nodes[k].beta * nodes[down_nbr].delta_up;
			nodes[down_nbr].delta_down = send_message(k, down_nbr).get_factor(nodes[down_nbr].delta_down.mask);
			pre_order(down_nbr);
		}
		nodes[k].init_marginal();
	}
	void DOWN_pass(){
		for(int root=0;root<n;root++){
			if(nodes[root].up_nbr!=-1)continue;
			for(auto down_nbr:nodes[root].down_nbrs){
				nodes[root].beta = nodes[root].beta * nodes[down_nbr].delta_up;
				nodes[down_nbr].delta_down = send_message(root, down_nbr).get_factor(nodes[down_nbr].delta_down.mask);
				pre_order(down_nbr);
			}
			nodes[root].init_marginal();
		}
	}
	void print(){
		printf("Clique Tree:\n");
		for(auto& node:nodes)
			node.print();
	}
};