#pragma once
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cstdio>
#include <cassert>
#include <vector>
#include <map>
#include <set>
#include <chrono>
using namespace std;
using namespace std::chrono;

#define NUM_CHARS 		(10)
#define MAX_WORD_LEN 	(6)
#define MAX_IMG_ID 		(1000)

#define SET(x, i) 	(x|=(1<<i))
#define CLR(x, i) 	(x&=~(1<<i))
#define TOG(x, i)		(x^=(1<<i))
#define BIT(x, i) 	((x>>i)&1)
#define CNT(x)			(__builtin_popcount(x))

#define UNION(x, y) (x|y)
#define INTER(x, y) (x&y)
#define IS_SUBSET(x, y) (x==(x&y))

#define MP 1
#define MAP 0

class Timer{
private:
	high_resolution_clock::time_point st;
public:
	void reset() {
		st =  high_resolution_clock::now();
	}
	float elapsed() {
		return (duration_cast<nanoseconds>(high_resolution_clock::now() - st)).count()/1000000.0;
	}
};