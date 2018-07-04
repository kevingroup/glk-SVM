//
// Created by cyhong on 4/10/2017.
//
#include "global.h"
#include <cstring>
#include "TreeLeafData.h"
#include "Sequence.h"
#ifndef GLK_KTREE_H
#define GLK_KTREE_H
#include <string>

class TreeNode {
public:
    bool isLeaf;
    TreeLeafData* leafData;
    TreeNode* child_head[MAX_ALPHABET_SIZE];
    TreeNode* child_tail[MAX_ALPHABET_SIZE];
    TreeNode();
};

class KTree {
private:
    TreeNode* root;
public:
    KTree();
    ~KTree();
    void deleteTree(TreeNode* cur, int alphabetSize);
    void printTree(TreeNode* cur, std::string s);
    void addGKmer(int seq_id, int *GKmer, int len_head, int len_tail, int len);
    void addSequence(int seq_id, int* bid, int len, int g, int l, int k, bool onlyMiddle);
    void addAllSeqs(Sequence **seqs, int nseqs, int g, int l, int k, bool addRC, bool onlyMiddle);
    void dfsTree_test(TreeNode* cur, int seq_id, int* GKmer, int x, int len_head, int len, int mmcnt, int g, int k, double** Kernel, std::string s);
    void dfsTree(TreeNode* cur, int seq_id, int* GKmer, int x, int len_head, int len, int mmcnt, int g, int k, double** Kernel);
    void dfsTree_norm(TreeNode* cur, int* GKmer, int x, int len_head, int len, int mmcnt, int g, int k, double* SUM);
    void dfsTree_InnerProd(TreeNode* cur, int* GKmer, int x, int len_head, int len, int mmcnt, int g, int k, double* InnerProd);
    void calKernel(Sequence **seqs, int nseqs, int g, int l, int k, double** Kernel, bool onlyMiddle);
    void calInnerProd(Sequence *sgi, int g, int l, int k, double* InnerProd, bool onlyMiddle);
    double calcNorm(Sequence *sgi, int g, int l, int k, bool onlyMiddle);
};


#endif //GLK_KTREE_H
