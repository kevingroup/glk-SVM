//
// Created by cyhong on 4/10/2017.
//

#include "KTree.h"
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <iostream>
#include <string>

using  namespace std;

int Comb[35][35];

void initComb()
{
    memset(Comb, 0, sizeof(Comb));
    Comb[1][0] = Comb[1][1] = 1;
    for(int i = 2; i <= 30; i++)
        Comb[i][0] = 1;
    for(int i = 2 ; i <= 30; i++)
        for(int j = 1; j <= i; j++)
            Comb[i][j] = Comb[i - 1][j] + Comb[i - 1][j - 1];

}
TreeNode::TreeNode() {
    this->isLeaf = false;
    this->leafData = NULL;
    memset(child_head, NULL, sizeof(TreeNode*) * MAX_ALPHABET_SIZE);
    memset(child_tail, NULL, sizeof(TreeNode*) * MAX_ALPHABET_SIZE);
}

KTree::KTree()
{
    root = new TreeNode();
    initComb();
}

KTree::~KTree()
{
    deleteTree(root, MAX_ALPHABET_SIZE);
}

void KTree::deleteTree(TreeNode* cur, int alphabetSize)
{
    for (int i = 0; i < alphabetSize; i++)
    {
        if (cur->child_head[i] != NULL)
        {
            deleteTree(cur->child_head[i], alphabetSize);
        }
        if (cur->child_tail[i] != NULL)
        {
            deleteTree(cur->child_tail[i], alphabetSize);
        }
    }
    if (cur->leafData != NULL)
        delete cur->leafData;
    if (cur != NULL)
        delete cur;
}

void KTree::addGKmer(int seq_id, int *GKmer, int len_head, int len_tail, int len)
{
    TreeNode *cur = root;
    for (int i = 0; i < len_head; i++)
    {
        int ch = GKmer[i];
        if (cur->child_head[ch] == NULL)
            cur->child_head[ch] = new TreeNode();
        cur = cur->child_head[ch];
    }
    for (int i = len - len_tail; i < len; i++) {
        int ch = GKmer[i];
        if (cur->child_tail[ch] == NULL)
            cur->child_tail[ch] = new TreeNode();
        cur = cur->child_tail[ch];
    }
    cur->isLeaf = true;
    if (cur->leafData == NULL)
        cur->leafData = new TreeLeafData();
    cur->leafData->add(seq_id);
}

void KTree::addSequence(int seq_id, int* bid, int len, int g, int l, int k, bool onlyMiddle)
{
    for (int i = 0; i <= len - (g + k); i++)
    {
        for (int j = g + k; j <= g + k + l; j++)
        {
            if (i + j > len)
                continue;
            if (onlyMiddle)
            {
                int len_head = (g + k) / 2;
                addGKmer(seq_id, bid + i, len_head, g + k - len_head, j);
            }
            else{
                for (int len_head = 0; len_head <= g + k; len_head++)
                    addGKmer(seq_id, bid + i, len_head, g + k - len_head, j);
            }
        }
    }
}

void KTree::printTree(TreeNode* cur, string s)
{
    if (cur->isLeaf)
    {
        cout<<s<<endl;
        TreeLeafData* x = cur->leafData;
        for (int i = 0; i < x->n; i++)
            printf("id: %d times: %d  ", x->seqInfo[2 * i], x->seqInfo[2 * i + 1]);
        printf("\n");
        return;
    }
    for (int i = 0; i < MAX_ALPHABET_SIZE; i++)
    {
        if (cur->child_head[i] != NULL)
        {
            string temp;
            if (i == 0)
                temp = "headA ";
            else if (i == 1)
                temp = "headC ";
            else if (i == 2)
                temp = "headG ";
            else temp = "headT ";
            printTree(cur->child_head[i], s + temp);
        }
        if (cur->child_tail[i] != NULL)
        {
            string temp;
            if (i == 0)
                temp = "tailA ";
            else if (i == 1)
                temp = "tailC ";
            else if (i == 2)
                temp = "tailG ";
            else temp = "tailT ";
            printTree(cur->child_tail[i], s + temp);
        }
    }
}

void KTree::addAllSeqs(Sequence **seqs, int nseqs, int g, int l, int k, bool addRC, bool onlyMiddle)
{
    for (int i = 0; i < nseqs; i++)
    {
        if (addRC)
            addSequence(i, seqs[i]->getSeqBaseIdRC(), seqs[i]->getLength(), g, l, k, onlyMiddle);
        addSequence(i, seqs[i]->getSeqBaseId(), seqs[i]->getLength(), g, l, k, onlyMiddle);
    }
    //printTree(root, "");
}
int cnt[1500];
void KTree::dfsTree_test(TreeNode* cur, int seq_id, int* GKmer, int x, int len_head, int len, int mmcnt, int g, int k, double** Kernel ,string s)
{
    if (x > len || mmcnt > g)
        return;
    if (cur->isLeaf)
    {
        int* seq_list = cur->leafData->seqInfo;
        cout<<s<<endl;
        printf("%d ---\n", cur->leafData->n);
        for (int i = 0; i < cur->leafData->n; i++)
        {
            cnt[seq_list[2 * i]] += seq_list[2 * i + 1];
            if (seq_id > seq_list[2 * i])
                continue;
            Kernel[seq_id][seq_list[2 * i]] += seq_list[2 * i + 1] * Comb[g + k - mmcnt][k];
            Kernel[seq_list[2 * i]][seq_id] += seq_list[2 * i + 1] * Comb[g + k - mmcnt][k];
            //printf("K[%d, %d] + %d * %d\n", seq_id, seq_list[2 * i], seq_list[2 * i + 1],  Comb[g + k - mmcnt][k]);
        }
        return;
    }
    if (x < len_head)                             //head part
    {
        if (cur->child_head[GKmer[x]] != NULL)
        {
            if (x == len_head - 1)
                dfsTree_test(cur->child_head[GKmer[x]], seq_id, GKmer, len - (g + k - len_head), len_head, len, mmcnt, g, k, Kernel, s + char(GKmer[x] + '0'));
            else
                dfsTree_test(cur->child_head[GKmer[x]], seq_id, GKmer, x + 1, len_head, len, mmcnt, g, k, Kernel, s + char(GKmer[x] + '0'));
        }
        for (int i = 0; i < MAX_ALPHABET_SIZE; i++)
        {
            if (cur->child_head[i] != NULL && i != GKmer[x])
            {
                if (x == len_head - 1)
                    dfsTree_test(cur->child_head[i], seq_id, GKmer, len - (g + k - len_head), len_head, len, mmcnt + 1, g, k, Kernel, s + char(i + '0'));
                else
                    dfsTree_test(cur->child_head[i], seq_id, GKmer, x + 1, len_head, len, mmcnt + 1, g, k, Kernel, s + char(i + '0'));
            }
        }
    }
    else{                                       //tail part
        if (x == 0)
            x = len - (g + k - len_head);
        if (cur->child_tail[GKmer[x]] != NULL)
            dfsTree_test(cur->child_tail[GKmer[x]], seq_id, GKmer, x + 1, len_head, len, mmcnt, g, k, Kernel, s + char(GKmer[x] + '0'));
        for (int i = 0; i < MAX_ALPHABET_SIZE; i++)
        {
            if (cur->child_tail[i] != NULL && i != GKmer[x])
                dfsTree_test(cur->child_tail[i], seq_id, GKmer, x + 1, len_head, len, mmcnt + 1, g, k, Kernel, s + char(i + '0'));
        }
    }
}
void KTree::dfsTree(TreeNode* cur, int seq_id, int* GKmer, int x, int len_head, int len, int mmcnt, int g, int k, double** Kernel)
{
    if (x > len || mmcnt > g)
        return;
    if (cur->isLeaf)
    {
        int* seq_list = cur->leafData->seqInfo;
        for (int i = 0; i < cur->leafData->n; i++)
        {
            if (seq_id > seq_list[2 * i])
                continue;
            Kernel[seq_list[2 * i]][seq_id] += seq_list[2 * i + 1] * Comb[g + k - mmcnt][k];
            //printf("K[%d, %d] + %d * %d\n", seq_id, seq_list[2 * i], seq_list[2 * i + 1],  Comb[g + k - mmcnt][k]);
        }
        return;
    }
    if (x < len_head)                             //head part
    {
        if (cur->child_head[GKmer[x]] != NULL)
        {
            if (x == len_head - 1)
                dfsTree(cur->child_head[GKmer[x]], seq_id, GKmer, len - (g + k - len_head), len_head, len, mmcnt, g, k, Kernel);
            else
                dfsTree(cur->child_head[GKmer[x]], seq_id, GKmer, x + 1, len_head, len, mmcnt, g, k, Kernel);
        }
        for (int i = 0; i < MAX_ALPHABET_SIZE; i++)
        {
            if (cur->child_head[i] != NULL && i != GKmer[x])
            {
                if (x == len_head - 1)
                    dfsTree(cur->child_head[i], seq_id, GKmer, len - (g + k - len_head), len_head, len, mmcnt + 1, g, k, Kernel);
                else
                    dfsTree(cur->child_head[i], seq_id, GKmer, x + 1, len_head, len, mmcnt + 1, g, k, Kernel);
            }
        }
    }
    else{                                       //tail part
        if (x == 0)
            x = len - (g + k - len_head);
        if (cur->child_tail[GKmer[x]] != NULL)
            dfsTree(cur->child_tail[GKmer[x]], seq_id, GKmer, x + 1, len_head, len, mmcnt, g, k, Kernel);
        for (int i = 0; i < MAX_ALPHABET_SIZE; i++)
        {
            if (cur->child_tail[i] != NULL && i != GKmer[x])
                dfsTree(cur->child_tail[i], seq_id, GKmer, x + 1, len_head, len, mmcnt + 1, g, k, Kernel);
        }
    }
}

void KTree::dfsTree_norm(TreeNode* cur, int* GKmer, int x, int len_head, int len, int mmcnt, int g, int k, double* norm)
{
    if (x > len || mmcnt > g)
        return;
    if (cur->isLeaf)
    {
        int* seq_list = cur->leafData->seqInfo;
        for (int i = 0; i < cur->leafData->n; i++)
        {
            *norm += seq_list[2 * i + 1] * Comb[g + k - mmcnt][k];
            //printf("K[%d, %d] + %d * %d\n", seq_id, seq_list[2 * i], seq_list[2 * i + 1],  Comb[g + k - mmcnt][k]);
        }
        return;
    }
    if (x < len_head)                             //head part
    {
        if (cur->child_head[GKmer[x]] != NULL)
        {
            if (x == len_head - 1)
                dfsTree_norm(cur->child_head[GKmer[x]], GKmer, len - (g + k - len_head), len_head, len, mmcnt, g, k, norm);
            else
                dfsTree_norm(cur->child_head[GKmer[x]], GKmer, x + 1, len_head, len, mmcnt, g, k, norm);
        }
        for (int i = 0; i < MAX_ALPHABET_SIZE; i++)
        {
            if (cur->child_head[i] != NULL && i != GKmer[x])
            {
                if (x == len_head - 1)
                    dfsTree_norm(cur->child_head[i], GKmer, len - (g + k - len_head), len_head, len, mmcnt + 1, g, k, norm);
                else
                    dfsTree_norm(cur->child_head[i], GKmer, x + 1, len_head, len, mmcnt + 1, g, k, norm);
            }
        }
    }
    else{                                       //tail part
        if (x == 0)
            x = len - (g + k - len_head);
        if (cur->child_tail[GKmer[x]] != NULL)
            dfsTree_norm(cur->child_tail[GKmer[x]], GKmer, x + 1, len_head, len, mmcnt, g, k, norm);
        for (int i = 0; i < MAX_ALPHABET_SIZE; i++)
        {
            if (cur->child_tail[i] != NULL && i != GKmer[x])
                dfsTree_norm(cur->child_tail[i], GKmer, x + 1, len_head, len, mmcnt + 1, g, k, norm);
        }
    }
}

void KTree::dfsTree_InnerProd(TreeNode* cur, int* GKmer, int x, int len_head, int len, int mmcnt, int g, int k, double* InnerProd)
{
    if (x > len || mmcnt > g)
        return;
    if (cur->isLeaf)
    {
        int* seq_list = cur->leafData->seqInfo;
        for (int i = 0; i < cur->leafData->n; i++)
        {
            InnerProd[seq_list[2 * i]] += seq_list[2 * i + 1] * Comb[g + k - mmcnt][k];
            //printf("K[%d, %d] + %d * %d\n", seq_id, seq_list[2 * i], seq_list[2 * i + 1],  Comb[g + k - mmcnt][k]);
        }
        return;
    }
    if (x < len_head)                             //head part
    {
        if (cur->child_head[GKmer[x]] != NULL)
        {
            if (x == len_head - 1)
                dfsTree_InnerProd(cur->child_head[GKmer[x]], GKmer, len - (g + k - len_head), len_head, len, mmcnt, g, k, InnerProd);
            else
                dfsTree_InnerProd(cur->child_head[GKmer[x]], GKmer, x + 1, len_head, len, mmcnt, g, k, InnerProd);
        }
        for (int i = 0; i < MAX_ALPHABET_SIZE; i++)
        {
            if (cur->child_head[i] != NULL && i != GKmer[x])
            {
                if (x == len_head - 1)
                    dfsTree_InnerProd(cur->child_head[i], GKmer, len - (g + k - len_head), len_head, len, mmcnt + 1, g, k, InnerProd);
                else
                    dfsTree_InnerProd(cur->child_head[i], GKmer, x + 1, len_head, len, mmcnt + 1, g, k, InnerProd);
            }
        }
    }
    else{                                       //tail part
        if (x == 0)
            x = len - (g + k - len_head);
        if (cur->child_tail[GKmer[x]] != NULL)
            dfsTree_InnerProd(cur->child_tail[GKmer[x]], GKmer, x + 1, len_head, len, mmcnt, g, k, InnerProd);
        for (int i = 0; i < MAX_ALPHABET_SIZE; i++)
        {
            if (cur->child_tail[i] != NULL && i != GKmer[x])
                dfsTree_InnerProd(cur->child_tail[i], GKmer, x + 1, len_head, len, mmcnt + 1, g, k, InnerProd);
        }
    }
}

void KTree::calKernel(Sequence **seqs, int nseqs, int g, int l, int k, double** Kernel, bool onlyMiddle, int num_of_thread)
{
    omp_set_num_threads(num_of_thread);
    #pragma omp parallel for
    for (int id = 0; id < nseqs; id++)
    {
        //printf("%d\n", id);
        int seq_len = seqs[id]->getLength();
        int* bid = seqs[id]->getSeqBaseId();
        for (int i = 0; i <= seq_len - (g + k); i++)
        {
            for (int j = g + k; j <= g + k + l; j++)
            {
                if (i + j > seq_len)
                    continue;
                if (onlyMiddle)
                {
                    int len_head = (g + k) / 2;
                    dfsTree(root, id, bid + i, 0, len_head, j, 0, g, k, Kernel);
                }
                else{
                    for (int len_head = 0; len_head <= g + k; len_head++)
                        dfsTree(root, id, bid + i, 0, len_head, j, 0, g, k, Kernel);
                }
            }
        }
    }
}


void KTree::calInnerProd(Sequence *sgi, int g, int l, int k, double* InnerProd, bool onlyMiddle)
{
    int seq_len = sgi->getLength();
    int* bid = sgi->getSeqBaseId();
    for (int i = 0; i <= seq_len - (g + k); i++)
    {
        for (int j = g + k; j <= g + k + l; j++)
        {
            if (i + j > seq_len)
                continue;
            if (onlyMiddle)
            {
                int len_head = (g + k) / 2;
                dfsTree_InnerProd(root, bid + i, 0, len_head, j, 0, g, k, InnerProd);
            }
            else{
                for (int len_head = 0; len_head <= g + k; len_head++)
                    dfsTree_InnerProd(root, bid + i, 0, len_head, j, 0, g, k, InnerProd);
            }
        }
    }
}

double KTree::calcNorm(Sequence *sgi, int g, int l, int k, bool onlyMiddle)
{
    double norm = 0;
    int seq_len = sgi->getLength();
    int* bid = sgi->getSeqBaseId();
    for (int i = 0; i <= seq_len - (g + k); i++)
    {
        for (int j = g + k; j <= g + k + l; j++)
        {
            if (i + j > seq_len)
                continue;
            if (onlyMiddle)
            {
                int len_head = (g + k) / 2;
                dfsTree_norm(root, bid + i, 0, len_head, j, 0, g, k, &norm);
            }
            else{
                for (int len_head = 0; len_head <= g + k; len_head++)
                    dfsTree_norm(root, bid + i, 0, len_head, j, 0, g, k, &norm);
            }
        }
    }
    return sqrt(norm);
}

