#pragma once
#include "common.h"

const char chars[NUM_CHARS] = {'a', 'd', 'e', 'h', 'i', 'n', 'o', 'r', 's', 't'};
	
class Model{
public:
	double trans[NUM_CHARS][NUM_CHARS], log_trans[NUM_CHARS][NUM_CHARS];
	double ocr[MAX_IMG_ID+1][NUM_CHARS], log_ocr[MAX_IMG_ID+1][NUM_CHARS];
	double skip, log_skip;
	double pair_skip, log_pair_skip;
	map<char, int> char_to_id;
	
	void init_char_to_id(){
		for(int i=0;i<NUM_CHARS;i++)
			char_to_id[chars[i]]=i;
	}
	int init(string trans_filename, string ocr_filename){
		init_char_to_id();
		if(init_transition_factors(trans_filename)){
			printf("Error parsing %s", trans_filename.c_str());
			return EXIT_FAILURE;
		}
		if(init_ocr_factors(ocr_filename)){
			printf("Error parsing %s", ocr_filename.c_str());
			return EXIT_FAILURE;
		}
		skip = 5.0;
		log_skip = log(skip);
		pair_skip = 5.0;
		log_pair_skip = log(pair_skip);
		return EXIT_SUCCESS;
	}
	int init_transition_factors(string filename){
		FILE* file = fopen(filename.c_str(), "r");
		if(!file){
			return EXIT_FAILURE;
		}
		int num_lines = NUM_CHARS * NUM_CHARS;
		char c1, c2;
		double p;
		for(int i=0;i<num_lines;i++){
			fscanf(file, " %c %c %lf", &c1, &c2, &p);
			trans[char_to_id[c1]][char_to_id[c2]] = p;
			log_trans[char_to_id[c1]][char_to_id[c2]] = log(p);
		}
		fclose(file);
		return EXIT_SUCCESS;
	}
	int init_ocr_factors(string filename){
		FILE* file = fopen(filename.c_str(), "r");
		if(!file){
			return EXIT_FAILURE;
		}
		int num_lines = MAX_IMG_ID * NUM_CHARS, id;
		char c;
		double p;
		for(int i=0;i<num_lines;i++){
			fscanf(file, "%d %c %lf", &id, &c, &p);
			ocr[id][char_to_id[c]] = p;
			log_ocr[id][char_to_id[c]] = log(p);
		}
		fclose(file);
		return EXIT_SUCCESS;
	}
};