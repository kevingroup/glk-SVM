//
// Created by cyhong on 10/18/2017.
//

#include <queue>
#include "Sequence.h"
using namespace std;
typedef long long LL;

#ifndef PATTERN_PATTERNEXTRACTOR_H
#define PATTERN_PATTERNEXTRACTOR_H

struct Node{
    int* pos;
    LL hash_value;
    int npos, nneg;
    int diff;
    int dis[25];
    bool operator < (const Node &x) const
    {
        return diff > x.diff;
    }
};


class PatternExtractor {
private:
    LL* pos_cnt;
    LL* neg_cnt;
    LL* Lmer_pos;
    LL* Lmer_neg;
    int* comb_pos;
    int num_pattern = 100;
    priority_queue<Node> heap;

public:
    int nposLmer = 0;
    int nnegLmer = 0;
    int max_diff = 0;
    int sum_diff = 0;

    PatternExtractor(Sequence** pos_seqs, Sequence** neg_seqs, int npos, int nneg, int LEN, int L, int maxmm, int num_pattern);
    PatternExtractor(Sequence** pos_seqs, Sequence** neg_seqs, int npos, int nneg, int LEN, int L, int indel_len, int maxmm, int num_pattern);
    ~PatternExtractor();

    void glkPattern(int L, int maxmm, int indel);
    void generateGapPosition(int x, int cur_mm, int* cur_pos, int L, int maxmm);
    void glkExtractor(int* position, int L, int maxmm);
    void printPattern(FILE* fout, int L, int maxmm, int indel);
};


#endif //PATTERN_PATTERNEXTRACTOR_H
