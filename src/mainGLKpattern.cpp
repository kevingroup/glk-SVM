//
// Created by cyhong on 10/18/2017.
//
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <algorithm>
#include <string>
#include <unistd.h>
#include "Sequence.h"
#include "global.h"
#include "PatternExtractor.h"

using namespace std;

typedef long long LL;

clock_t start;
typedef struct {
    int L; // tot_length
    int K;
    int maxmm;
    int indel;
    int from;
    int to;
    int num_pattern;
    int maxseqlen;
    int maxnumseq;
    bool addRC;
    char *posfile;
    char *negfile;
    char *outfile;
} OptsGLK;

int glkParameter(OptsGLK &opt)
{
    int L = opt.L;
    int maxmm = opt.maxmm;
    int from = opt.from;
    int to = opt.to;
    int num_pattern = opt.num_pattern;
    int maxseqlen =	opt.maxseqlen;
    int nMAXSEQUENCES = opt.maxnumseq;
    char* posSeqFN = opt.posfile;
    char* negSeqFN = opt.negfile;
    char* outFN = opt.outfile;
    //bool addRC = opt.addRC;

    Sequence **pos_seqs = new Sequence *[nMAXSEQUENCES];
    Sequence **neg_seqs = new Sequence *[nMAXSEQUENCES];
    int npos = 0;
    int nneg = 0;

    int tot_len = 0;
    FILE *sfi = fopen(posSeqFN, "r");
    if (sfi == NULL)
    {
        perror("open pos_seq file error");
        return 0;
    }
    int seq_cnt = 0;
    while (!feof(sfi))
    {
        Sequence *sgi = new Sequence(maxseqlen + 3);
        sgi->readFasta(sfi);
        if (sgi->getLength() > 0)
        {
            //if (seq_cnt++ % 2 == 0)
            //    continue;
            pos_seqs[npos] = sgi;
            npos++;
            tot_len += sgi->getLength();
        }
    }
    fclose(sfi);

    sfi = fopen(negSeqFN, "r");
    if (sfi == NULL)
    {
        perror("open neg_seq file error");
        return 0;
    }
    seq_cnt = 0;
    while (!feof(sfi))
    {
        Sequence *sgi = new Sequence(maxseqlen + 3);
        sgi->readFasta(sfi);
        if (sgi->getLength() > 0)
        {
            //if (seq_cnt++ % 2 == 0)
            //    continue;
            neg_seqs[nneg] = sgi;
            nneg++;
        }
    }
    fclose(sfi);

    //cout<<npos<<" "<<nneg<<endl;

    FILE* fout = fopen(outFN, "w");

    for (int indel = from; indel <= to; indel += 1)
    {
        PatternExtractor* PE = new PatternExtractor(pos_seqs, neg_seqs, npos, nneg, tot_len, L, indel, maxmm, num_pattern);

        int pos[maxmm];

        PE->generateGapPosition(0, 0, pos, L, maxmm);

        PE->glkPattern(L, maxmm, indel);

        //printf("g = %d, l = %d, k = %d : diff = %d\n", opt.maxmm, indel, opt.L-opt.maxmm, PE->max_diff);
        //printf("g = %d, l = %d, k = %d : diff = %d\n", opt.maxmm, indel, opt.L-opt.maxmm, PE->sum_diff);
        fprintf(fout, "g = %d, l = %d, k = %d : maximal difference = %d\n", opt.maxmm, indel, opt.L-opt.maxmm, PE->max_diff);

        PE->printPattern(fout, L, maxmm, indel);

        fprintf(fout, "\n");

        delete PE;
    }

    fclose(fout);

    for(int i = 0; i < npos; i++)
    {
        delete pos_seqs[i];
    }
    delete []pos_seqs;
    for(int i = 0; i < nneg; i++)
    {
        delete neg_seqs[i];
    }
    delete []neg_seqs;

    return 0;
}

void print_usage_and_exit(char *prog)
{
    cout << endl;
    cout << " Usage: " << prog << " [options] <pos_seqfile> <neg_seqfile> <outfile>" << endl;
    cout << endl;
    cout << "  generating the most differentiate glk patterns between the positive" << endl;
    cout << "  and negative sequences." << endl;
    cout << endl;
    cout << " Arguments:" << endl;
    cout << "  pos_seqfile: positive sequence file name (fasta format)" << endl;
    cout << "  neg_seqfile: negative sequence file name (fasta format)" << endl;
    cout << "  outfile: output file name" << endl;
    cout << endl;
    cout << " Options:" << endl;
    cout << "  -g G           set gap length, default=" << DEF_G << endl;
    cout << "  -f from        set minimal indel length, default=" << DEF_L << endl;
    cout << "  -t to          set maximal indel length, default=" << DEF_L << endl;
    cout << "  -k K           set number of exact match character, default=" << DEF_K << endl;
    cout << "  -m maxSeqLen   set maximum sequence length in the sequence files," << endl;
    cout << "                 default=" << DEF_MAXSEQLEN << endl;
    cout << "  -n numGLK      set the number of top glk pattern with largest occurrence difference," << endl;
    cout << "                 default=" << 10 << endl;
    //cout << "  -R             if set, reverse complement sequences will be considered" << endl;
    cout << endl;

    exit(0);
}

int main(int argc, char** argv)
{

    OptsGLK opt;

    ::optind = 1; // reset getopt()

    int c;

    opt.maxmm = DEF_G;
    opt.K = DEF_K;

    opt.indel = DEF_L;
    opt.from = DEF_L;
    opt.to = DEF_L;

    opt.maxseqlen = DEF_MAXSEQLEN;
    opt.maxnumseq = DEF_MAXNUMSEQ;
    opt.num_pattern = 10;
    //opt.addRC = false;

    while ((c = getopt (argc, argv, "g:f:t:k:d:m:n:R")) != -1)
    {
        switch (c)
        {
            case 'g':
                //printf("\ngetopt  l = %s = %d\n",optarg,atoi(optarg)  );
                opt.maxmm = atoi(optarg);
                break;
            case 'f':
                opt.from = atoi(optarg);
                break;
            case 't':
                opt.to = atoi(optarg);
                break;
            case 'k':
                opt.K = atoi(optarg);
                break;
            case 'm':
                opt.maxseqlen = atoi(optarg);
                break;
            case 'n':
                opt.num_pattern = atoi(optarg);
                break;
            case 'R':
                opt.addRC = true;
                break;
            default:
                return 0;
        }
    }

    if (argc-optind != 3)
    {
        print_usage_and_exit(argv[0]);
        return 0;
    }

    int index = optind;
    opt.posfile = argv[index++];
    opt.negfile = argv[index++];
    opt.outfile = argv[index++];

    opt.L = opt.maxmm + opt.K;

    glkParameter(opt);

    return 0;
}