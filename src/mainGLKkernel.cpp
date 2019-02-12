/*
This code is partially adopted from gkm-SVM(Mahmoud Ghandi, http://www.beerlab.org/gkmsvm/)
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include <algorithm>
#include <set>
#include <unistd.h>
#include "Sequence.h"
#include "global.h"
#include "KTree.h"

using namespace std;
clock_t start;
typedef struct {
    int G;
    int L;
    int K;
    int maxseqlen;
    int maxnumseq;
    int num_of_thread;
    bool addRC;
    bool onlyMiddle;
    char *posfile;
    char *negfile;
    char *outfile;
} OptsGLK;
int C[50][50] = {0};
char getMatch[100] = {0};
void initCombination()
{
    getMatch['A'] = 'T';
    getMatch['T'] = 'A';
    getMatch['C'] = 'G';
    getMatch['G'] = 'C';
    C[1][0] = C[1][1] = 1;
    for(int i = 2; i <= 20; i++)
        C[i][0] = 1;
    for(int i = 2 ; i <= 20; i++)
        for(int j = 1; j <= i; j++)
            C[i][j] = C[i - 1][j] + C[i - 1][j - 1];
}

void print_usage_and_exit(char *prog)
{
    cout << endl;
    cout << " Usage: " << prog << " [options] <pos_seqfile> <neg_seqfile> <outfile>" << endl;
    cout << endl;
    cout << "  generates a lower triangle of kernel matrix (i.e. pairwise similarities)" << endl;
    cout << "  between the sequences." << endl;
    cout << endl;
    cout << " Arguments:" << endl;
    cout << "  pos_seqfile: positive sequence file name (fasta format)" << endl;
    cout << "  neg_seqfile: negative sequence file name (fasta format)" << endl;
    cout << "  outfile: output file name" << endl;
    cout << endl;
    cout << " Options:" << endl;
    cout << "  -g G           set gap length, default=" << DEF_G << endl;
    cout << "  -l L           set indel length, default=" << DEF_L << endl;
    cout << "  -k K           set number of exact match character, default=" << DEF_K << endl;
    cout << "  -m maxSeqLen   set maximum sequence length in the sequence files," << endl;
    cout << "                 default=" << DEF_MAXSEQLEN << endl;
    cout << "  -n maxNumSeq   set maximum number of sequences in the sequence files," << endl;
    cout << "                 default=" << DEF_MAXNUMSEQ << endl;
    cout << "  -t num_thread  set the number of thread used by the program," << endl;
    cout << "                 default=" << DEF_NUM_OF_THREAD << endl;
    cout << "  -R             if set, reverse complement sequences will be considered" << endl;
    cout << "  -A             if set, indel could be anywhere" << endl;
    cout << endl;

    exit(0);
}

int glkKernel(OptsGLK &opt)
{
    int G = opt.G;
    int L = opt.L;
    int K = opt.K;
    int maxseqlen =	opt.maxseqlen;
    int nMAXSEQUENCES = opt.maxnumseq;
    int num_of_thread = opt.num_of_thread;
    char* posSeqFN = opt.posfile;
    char* negSeqFN = opt.negfile;
    char* outFN = opt.outfile;
    bool addRC = opt.addRC;
    bool onlyMiddle = opt.onlyMiddle;

    Sequence **seqs = new Sequence *[nMAXSEQUENCES];
    Sequence *sgi = new Sequence(maxseqlen + 3);
    int npos = 0;
    int nneg = 0;

    int nseqs = 0;

    FILE *sfi = fopen(posSeqFN, "r");
    if (sfi == NULL)
    {
        perror("open pos_seq file error");
        return 0;
    }
    while (!feof(sfi))
    {
        sgi->readFasta(sfi);
        if (sgi->getLength() > 0)
        {
            seqs[nseqs] = new Sequence(maxseqlen + 3, sgi);
            nseqs++;
        }
    }
    fclose(sfi);
    npos = nseqs;

    sfi = fopen(negSeqFN, "r");
    if (sfi == NULL)
    {
        perror("open neg_seq file error");
        return 0;
    }
    while (!feof(sfi))
    {
        sgi->readFasta(sfi);
        if (sgi->getLength() > 0)
        {
            seqs[nseqs] = new Sequence(maxseqlen + 3, sgi);
            nseqs++;
        }
    }
    fclose(sfi);

    nneg = nseqs - npos;
    cout<<npos<<" "<<nneg<<endl;

    double **Kernel = new double*[nseqs];
    for (int i = 0; i < nseqs; i++)
    {
        Kernel[i] = new double[nseqs];
        for (int j = 0; j < nseqs; j++)
            Kernel[i][j] = 0;
    }


    KTree *Tree = new KTree();

    Tree->addAllSeqs(seqs, nseqs, G, L, K, addRC, onlyMiddle);
    Tree->calKernel(seqs, nseqs, G, L, K, Kernel, onlyMiddle, num_of_thread);

    double *norm = new double [nMAXSEQUENCES];
    int i = 0;
    for(i = 0; i < nseqs; i++)
    {
        norm[i] = sqrt(Kernel[i][i]);
        //cout<<"norm: "<<norm[i]<<endl;
    }

    FILE *fo = fopen(outFN, "w");

    for(i = 0; i < nseqs; i++)
    {
        //printf("FILE: %s  i = %d\n", posSeqFN, i);
        for(int j = 0; j < nseqs; j++)
        {
            if(i > j)
            {
                fprintf(fo, "%e\t",(norm[i] * norm[j] < 1E-30) ? 0.0 : Kernel[i][j] / (norm[i] * norm[j]));
            }
            else if (i==j)
            {
                fprintf(fo, "1.0\t");
            }
        }
        fprintf(fo, "\n");
    }

    fclose(fo);
    for(i = 0; i < nseqs; i++)
    {
        delete seqs[i];
        delete []Kernel[i];
    }
    delete []Kernel;
    delete []seqs;
    delete []norm;
    delete Tree;
    delete sgi;
    
    return 0;
}

int main(int argc, char** argv)
{

    OptsGLK opt;

    ::optind = 1; // reset getopt()

    int c;

    opt.G = DEF_G;
    opt.L = DEF_L;
    opt.K = DEF_K;
    opt.maxseqlen = DEF_MAXSEQLEN;
    opt.maxnumseq = DEF_MAXNUMSEQ;
    opt.num_of_thread = DEF_NUM_OF_THREAD;
    opt.addRC = false;
    opt.onlyMiddle = true;

    while ((c = getopt (argc, argv, "g:l:k:d:m:n:t:RA")) != -1)
    {
        switch (c)
        {
            case 'g':
                opt.G = atoi(optarg);
                break;
            case 'l':
                opt.L = atoi(optarg);
                break;
            case 'k':
                opt.K = atoi(optarg);
                break;
            case 'm':
                opt.maxseqlen = atoi(optarg);
                break;
            case 'n':
                opt.maxnumseq = atoi(optarg);
                break;
            case 't':
                opt.num_of_thread = atoi(optarg);
                break;
            case 'R':
                opt.addRC = true;
                break;
            case 'A':
                opt.onlyMiddle = false;
                break;
            default:
                print_usage_and_exit(argv[0]);
        }
    }

    if (argc-optind != 3)
    {
        print_usage_and_exit(argv[0]);
    }

    int index = optind;
    opt.posfile = argv[index++];
    opt.negfile = argv[index++];
    opt.outfile = argv[index++];

    initCombination();

    //clock_t end, tot_Start; // typedef long clock_t
    //tot_Start = clock();
    //start = clock();

    glkKernel(opt);

    //end = clock();
    //double duration =(double)(end - tot_Start)/CLOCKS_PER_SEC;
    //cout<<"Total time: "<<duration / 60<<" minutes"<<endl;

    return 0;
}
