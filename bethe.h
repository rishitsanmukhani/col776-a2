#pragma once
#include "graph.h"
using namespace std;

struct BetheGraph{
	struct FactorNode{
		int id;
		Factor psi;
		vector<int> nbrs;
		void init(int _id, const Factor& factor){
			id = _id;
			psi = factor;
			int i=0, mask=psi.mask;
			while(mask){
				if(mask&1){
					nbrs.push_back(i);
				}
				mask>>=1;
				i++;
			}
		}
		void print(){
			printf("FactorNode: %d\n", id);
			for(auto n:nbrs){
				cout<<n<<" ";
			}
			cout<<endl;
		}
	};
	struct VarNode{
		int id;
		Factor beta;
		Factor marginal;
		vector<int> nbrs;
		void init(int i){
			id = i;
			beta.init(1<<id);
		}
		void init_marginal(){
			marginal = beta;
			double sum=0;
			for(auto val:marginal.v)sum += val;
			for(auto& val:marginal.v)val /= sum;			
		}
		void print(){
			printf("VarNode: %d\n", id);
			for(auto n:nbrs){
				cout<<n<<" ";
			}
			cout<<endl;
		}
	};
	vector<FactorNode> factor_nodes;
	vector<VarNode> var_nodes;
	map<pair<int, int>, Factor> delta;
	double time_taken;
	void print(){
		for(auto& node: factor_nodes){
			node.print();
		}
		for(auto& node:var_nodes){
			node.print();
		}
	}
	void init(const Graph& G){
		int i=0;
		var_nodes.resize(G.n);
		for(auto& node:var_nodes){
			node.init(i);
			i++;
		}
		factor_nodes.resize(G.factors.size());
		i=0;
		for(auto& p:G.factors){
			factor_nodes[i].init(i, p.second);
			for(auto j:factor_nodes[i].nbrs){
				delta[{i,j}] = Factor();
				delta[{i,j}].init(1<<j);
				delta[{j,i}] = Factor();
				delta[{j,i}].init(1<<j);
				var_nodes[j].nbrs.push_back(i);
			}
			i++;
		}
		calibrate();
	}
	void calibrate(){
		int itr=0;
		time_taken = 0;
		Timer t;t.reset();
		while(1){
			for(auto& node:factor_nodes){
				for(auto nbr:node.nbrs){
					delta[{node.id, nbr}] = node.psi;
					for(auto x:node.nbrs){
						if(x==nbr)continue;
						delta[{node.id, nbr}] = delta[{node.id, nbr}] * delta[{x, node.id}];
					}
					delta[{node.id, nbr}] = delta[{node.id, nbr}].get_factor(1<<nbr);
				}
			}
			double diff=0;
			for(auto& node:var_nodes){
				Factor beta;
				beta.init(1<<node.id);
				for(auto nbr:node.nbrs){
					beta = beta * delta[{nbr, node.id}];
					delta[{node.id, nbr}].init(1<<node.id);
					for(auto x:node.nbrs){
						if(x==nbr)continue;
						delta[{node.id, nbr}] = delta[{node.id, nbr}] * delta[{x, node.id}];
					}
				}
				for(int i=0;i<beta.v.size();i++){
					diff = max(diff, abs(beta.v[i] - node.beta.v[i]));
					node.beta.v[i] = beta.v[i];
				}
			}
			if(itr>100 || diff < (1e-12)){
				break;
			}
			itr++;
		}
		for(auto&node:var_nodes){
			node.init_marginal();
		}
		time_taken = t.elapsed();
	}
};