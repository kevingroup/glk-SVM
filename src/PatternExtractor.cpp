//
// Created by cyhong on 10/18/2017.
//

#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <bitset>
#include "PatternExtractor.h"
#include "global.h"

PatternExtractor::PatternExtractor(Sequence** pos_seqs, Sequence** neg_seqs, int npos, int nneg, int LEN, int L, int maxmm, int npattern)
{
    pos_cnt = new LL[LEN];
    neg_cnt = new LL[LEN];
    num_pattern = npattern;
    comb_pos = new int[NchooseK(L, maxmm) * maxmm];

    Lmer_pos = new LL[LEN];
    nposLmer = 0;
    for (int i = 0; i < npos; i++)
    {
        int seqlen = pos_seqs[i]->getLength();
        int* seq = pos_seqs[i]->getSeqBaseId();
        LL MASK = (1 << (2 * L)) - 1;
        LL hash_value = 0;
        for (int j = 0; j < L; j++)
            hash_value = LL((hash_value << 2) + seq[j]);
        Lmer_pos[nposLmer++] = hash_value;
        for (int j = L; j < seqlen; j++)
        {
            hash_value = LL((hash_value << 2) + seq[j]);
            hash_value &= MASK;
            Lmer_pos[nposLmer++] = hash_value;
        }
    }

    Lmer_neg = new LL[LEN];
    nnegLmer = 0;
    for (int i = 0; i < nneg; i++)
    {
        int seqlen = neg_seqs[i]->getLength();
        int* seq = neg_seqs[i]->getSeqBaseId();
        LL MASK = (1 << (2 * L)) - 1;
        LL hash_value = 0;
        for (int j = 0; j < L; j++)
            hash_value = (hash_value << 2) + seq[j];
        Lmer_neg[nnegLmer++] = hash_value;
        for (int j = L; j < seqlen; j++)
        {
            hash_value = (hash_value << 2) + seq[j];
            hash_value &= MASK;
            Lmer_neg[nnegLmer++] = hash_value;
        }
    }

}

int comb_id;
PatternExtractor::PatternExtractor(Sequence** pos_seqs, Sequence** neg_seqs, int npos, int nneg, int LEN, int L, int indel_len, int maxmm, int npattern)
{
    max_diff = 0;
    sum_diff = 0;
    comb_id = 0;
    pos_cnt = new LL[LEN * (indel_len + 1)];
    neg_cnt = new LL[LEN * (indel_len + 1)];
    num_pattern = npattern;
    comb_pos = new int[NchooseK(L, maxmm) * maxmm];

    Lmer_pos = new LL[LEN * (indel_len + 1)];
    nposLmer = 0;
    for (int i = 0; i < npos; i++)
    {
        int seqlen = pos_seqs[i]->getLength();
        int* seq = pos_seqs[i]->getSeqBaseId();
        LL MASK = (1 << (2 * L)) - 1;
        for (int k = 0; k <= seqlen - L; k++)
            for (int len = 0; len <= indel_len; len++)
            {
                if (k + len + L > seqlen) break;
                LL hash_value = 0;
                for (int j = 0; j < L / 2; j++)
                    hash_value = LL((hash_value << 2) + seq[k + j]);
                for (int j = L / 2 + len; j < L + len; j++)
                    hash_value = LL((hash_value << 2) + seq[k + j]);

                Lmer_pos[nposLmer++] = LL(hash_value << 5) + len;

                //std::cout<<bitset<64>(Lmer_pos[nposLmer - 1])<<endl;
            }
    }

    Lmer_neg = new LL[LEN * (indel_len + 1)];
    nnegLmer = 0;
    for (int i = 0; i < nneg; i++)
    {
        int seqlen = neg_seqs[i]->getLength();
        int* seq = neg_seqs[i]->getSeqBaseId();
        LL MASK = (1 << (2 * L)) - 1;
        for (int k = 0; k <= seqlen - L; k++)
            for (int len = 0; len <= indel_len; len++)
            {
                if (k + len + L > seqlen) break;
                LL hash_value = 0;
                for (int j = 0; j < L / 2; j++)
                    hash_value = LL((hash_value << 2) + seq[k + j]);
                for (int j = L / 2 + len; j < L + len; j++)
                    hash_value = LL((hash_value << 2) + seq[k + j]);

                Lmer_neg[nnegLmer++] = LL(hash_value << 5) + len;
                //std::cout<<bitset<64>(Lmer_neg[nnegLmer - 1])<<endl;
            }
    }

}

PatternExtractor::~PatternExtractor()
{
    while (!heap.empty())
        heap.pop();
    delete []pos_cnt;
    delete []neg_cnt;
    delete []comb_pos;
    delete []Lmer_pos;
    delete []Lmer_neg;
}

void PatternExtractor::generateGapPosition(int x, int cur_mm, int* cur_pos, int L, int maxmm)
{
    if (cur_mm == maxmm)
    {
        for (int i = 0; i < maxmm; i++)
            comb_pos[comb_id * maxmm + i] = cur_pos[i];
        comb_id++;
        return;
    }
    if (x == L)
        return;
    generateGapPosition(x + 1, cur_mm, cur_pos, L, maxmm);
    cur_pos[cur_mm] = x;
    generateGapPosition(x + 1, cur_mm + 1, cur_pos, L, maxmm);
}

void PatternExtractor::glkPattern(int L, int maxmm, int indel) //add gap to Lmer and store in pos_cnt and neg_cnt
{
    int num_of_case = NchooseK(L, maxmm);

    for (int i = 0; i < num_of_case; i++)
    {
        //printf("Case = %d\n", i);
        /*
        for (int j = 0; j < maxmm; j++)
            printf("%d ", comb_pos[i * maxmm + j]);
        printf("\n");
        */

        LL MASK = (1 << (2 * L)) - 1;
        for (int j = 0; j < maxmm; j++)
        {
            //cout<<comb_pos[i * maxmm + j]<<" ";
            int p = L - comb_pos[i * maxmm + j];
            LL MASK_p = ((1 << (2 * p)) - 1) >> (2 * p - 2) << (2 * p - 2);
            //cout<<bitset<64>(MASK)<<"\n"<<bitset<64>(MASK_p)<<endl;
            MASK -= MASK_p;

        }
        //cout<<bitset<64>(MASK)<<endl;

        MASK <<= 5;
        MASK += 31;

        //std::cout<<bitset<32>(MASK)<<endl;

        for (int j = 0; j < nposLmer; j++)
        {
            pos_cnt[j] = Lmer_pos[j] & MASK;
        }

        for (int j = 0; j < nnegLmer; j++)
        {
            neg_cnt[j] = Lmer_neg[j] & MASK;
            //std::cout<<bitset<64>(neg_cnt[j])<<"  "<<bitset<64>(Lmer_neg[j])<<endl;
        }

        //std::cout<<bitset<64>(MASK)<<endl;

        glkExtractor(comb_pos + i * maxmm, L, maxmm);
    }
    /*
    while (!heap.empty())
    {
        sum_diff += heap.top().diff;
        heap.pop();
    }
     */
}

void PatternExtractor::glkExtractor(int* position, int L, int maxmm)
{
    std::sort(pos_cnt, pos_cnt + nposLmer);
    std::sort(neg_cnt, neg_cnt + nnegLmer);

    int p = 0, q = 0;
    int p_cnt, q_cnt;

    int dis[25];
    while (p < nposLmer)
    {
        memset(dis, 0, sizeof(dis));
        dis[pos_cnt[p] % 32]++;
        p_cnt = 1;
        q_cnt = 1;
        while ((pos_cnt[p + 1] >> 5) == (pos_cnt[p] >> 5))
        {
            dis[pos_cnt[p + 1] % 32]++;
            p++;
            p_cnt++;
        }
        int raw_q = q;
        while ((neg_cnt[q + 1] >> 5) == (neg_cnt[q] >> 5))
        {
            q++;
            q_cnt++;
        }

        if ((pos_cnt[p] >> 5) == (neg_cnt[q] >> 5)) // common pattern
        {
            if (p_cnt > q_cnt)
            {
                int cur_diff = p_cnt - q_cnt;

                if (heap.size() < num_pattern || cur_diff > heap.top().diff)
                {
                    Node node;
                    node.pos = new int[maxmm];
                    for (int i = 0; i < 25; i++)
                        node.dis[i] = dis[i];
                    for (int i = 0; i < maxmm; i++)
                        node.pos[i] = position[i];
                    node.hash_value = (pos_cnt[p] >> 5);
                    node.npos = p_cnt;
                    node.nneg = q_cnt;
                    node.diff = cur_diff;

                    if (heap.size() == num_pattern)
                        heap.pop();
                    heap.push(node);
                    if (cur_diff > max_diff)
                        max_diff = cur_diff;
                }
            }
            p++;
            q++;
            continue;
        }

        if ((pos_cnt[p] >> 5) < (neg_cnt[q] >> 5))   // pattern only occur in pos
        {
            int cur_diff = p_cnt;
            if (heap.size() < num_pattern || cur_diff > heap.top().diff)
            {
                Node node;
                node.pos = new int[maxmm];
                for (int i = 0; i < maxmm; i++)
                    node.pos[i] = position[i];
                node.hash_value = (pos_cnt[p] >> 5);
                node.npos = p_cnt;
                node.nneg = 0;
                node.diff = cur_diff;

                if (heap.size() == num_pattern)
                    heap.pop();
                heap.push(node);
            }
            p++;
            q = raw_q; // jump back to original q
            continue;
        }
        while (pos_cnt[p] > neg_cnt[q]) //skip pattern only occur in neg
        {
            q++;
        }
    }
    //printf("%d %d\n", heap.size(), heap.top().diff);
}

void PatternExtractor::printPattern(FILE* fout, int L, int maxmm, int indel)
{
    char converter[4];
    converter[0] = 'A';
    converter[1] = 'C';
    converter[2] = 'G';
    converter[3] = 'T';
    char pattern[L];
    char pattern_indel[L + 1];

    while (!heap.empty())
    {
        LL x = heap.top().hash_value;
        int* position = heap.top().pos;
        for (int i = L - 1; i >= 0; i--)
        {
            pattern[i] = converter[x % 4];
            x >>= 2;
        }
        for (int i = 0; i < maxmm; i++)
            pattern[position[i]] = '*';
        pattern[L] = 0;

        /*
        printf("pos ");
        for (int i = 0; i < maxmm; i++)
            printf("%d ", position[i]);
        printf("\n");

        printf("%s\n", pattern);
        */
        for (int i = 0; i < L / 2; i++)
            pattern_indel[i] = pattern[i];
        pattern_indel[L / 2] = '-';
        for (int i = L / 2; i < L; i++)
            pattern_indel[i + 1] = pattern[i];

        pattern_indel[L+1] = 0;
        fprintf(fout, "%s  pos_count = %d  neg_count = %d  diff = %d\n", pattern_indel, heap.top().npos, heap.top().nneg, heap.top().diff);
        fprintf(fout, "Positive spacing count distribution (from 0 to %d):\n", indel);
        for (int i = 0; i <= indel; i++)
            fprintf(fout, "%d ", heap.top().dis[i]);
        fprintf(fout, "\n");
        //cout<<heap.top().diff<<endl;
        heap.pop();
    }

}