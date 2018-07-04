//
// Created by cyhong on 10/16/2017.
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <iostream>

#include "global.h"
#include "Sequence.h"
#include "KTree.h"
#include "SVMtrain.h"
#include "SequenceNames.h"

using namespace std;

typedef struct {
    int G;
    int L;
    int K;
    int maxseqlen;
    int maxnumseq;
    bool addRC;
    bool onlyMiddle;
    char *seqfile;
    char *svseqfile;
    char *alphafile;
    char *outfile;

} OptsSVMClassify;

double calcnorm(Sequence *sgi, int g, int l, int k, bool onlyMiddle);
double svmScoreunorm(int i, double *c);

int svmClassifySuffixTree(OptsSVMClassify &opt);
int svmClassifySimple(OptsSVMClassify &opt);

void print_usage_and_exit(char *prog)
{
    cout << endl;
    cout << " Usage: " << prog << " [options] <test_seqfile> <sv_seqfile> <sv_alphafile> <outfile>" << endl;
    cout << endl;
    cout << "  given support vectors SVs and corresponding coefficients alphas and a set of " << endl;
    cout << "  sequences, calculates the SVM scores for the sequences." << endl;
    cout << endl;
    cout << " Arguments:" << endl;
    cout << "  test_seqfile: sequence file name to test/score (fasta format)" << endl;
    cout << "  sv_seqfile: sequence file name containing support vectors (fasta format)" << endl;
    cout << "  sv_alphafile: coefficient file name containing alphas for support vectors" << endl;
    cout << "  outfile: output file name" << endl;
    cout << endl;
    cout << " Options:" << endl;
    cout << "  -g G           set number of gap, default=" << DEF_G << endl;
    cout << "  -l L           set word length, default=" << DEF_L << endl;
    cout << "  -k K           set number of exact match characters, default=" << DEF_K << endl;
    cout << "  -m maxSeqLen   set maximum sequence length in the sequence files," << endl;
    cout << "                 default=" << DEF_MAXSEQLEN << endl;
    cout << "  -n maxNumSeq   set maximum number of sequences in the sequence files," << endl;
    cout << "                 default=" << DEF_MAXNUMSEQ << endl;
    cout << "  -R             if set, reverse complement sequences will NOT be considered" << endl;
    cout << "  -A             if set, indel could be anywhere" << endl;
    cout << endl;

    exit(0);
}

//given fasta file for SVs and the corresponding weights, outputs and another file for the test sequences, gives the SVM score
int svmClassifySuffixTree(OptsSVMClassify &opt)
{
    int i;
    int G = opt.G;
    int L = opt.L;
    int K = opt.K;
    int maxseqlen =	opt.maxseqlen;
    int nMAXSEQUENCES = opt.maxnumseq;
    bool addRC = opt.addRC;
    bool onlyMiddle = opt.onlyMiddle;

    char *SeqsFN = opt.seqfile;
    char *SVSeqsFN = opt.svseqfile;
    char *SVSeqIDsFN = opt.alphafile;
    char *outFN = opt.outfile;

    char **seqname = new char *[nMAXSEQUENCES];

    //CSequenceNames *svsn= new CSequenceNames();

    char *seqNames[MAXNSVSeqs];
    double weight[MAXNSVSeqs];
    char stmp[10000];

    int NSVseqs = 0;
    FILE *f = fopen(SVSeqIDsFN, "r") ;
    while (!feof(f))
    {
        if (fgets(stmp, 10000-5, f))
        {
            if (stmp[0]!=0) {

                seqNames[NSVseqs] = new char[MAXSeqnameLENGTH];
                sscanf(stmp, "%s%lf", seqNames[NSVseqs], weight + NSVseqs);
                //printf("%s\n%e\n",seqNames[Nseqs],weight[Nseqs]);
                NSVseqs++;
            }
        }
    }

    printf("\n  %d SV ids read. \n", NSVseqs);

    KTree *SVTree= new KTree(); //keeps all the sequences of (g+k)mer in support vectors

    double *alpha = new double[NSVseqs + 5];
    double *SV_norm = new double[NSVseqs + 5];

    FILE *sfi = fopen(SVSeqsFN, "r"); //read sv sequences

    if (sfi == NULL) {
        perror("error occurred while opening SV sequences file");
        return 0;
    }

    int nsvseqs = 0;
    while (!feof(sfi))
    {
        Sequence *sgi = new Sequence(maxseqlen + 3);
        sgi->readFasta(sfi);
        if(sgi->getLength()>0)
        {
            alpha[nsvseqs] = weight[nsvseqs];
            SVTree->addSequence(nsvseqs, sgi->getSeqBaseId(), sgi->getLength(), G, L, K, onlyMiddle);
            if(addRC)
            {
                SVTree->addSequence(nsvseqs, sgi->getSeqBaseIdRC(), sgi->getLength(), G, L, K, onlyMiddle);
            }

            SV_norm[nsvseqs] = calcnorm(sgi, G, L, K, onlyMiddle);
            printf("%d = %.2lf %.2lf\n", nsvseqs, SV_norm[nsvseqs], alpha[nsvseqs]);
            nsvseqs++;
        }
    }

    printf("  %d SV seqs read \n", NSVseqs);

    sfi = fopen(SeqsFN, "r"); // read sequences which need to classify
    if (sfi == NULL)
    {
        perror ("error occurred while opening a file");
        return 0;
    }
    FILE *fo = fopen(outFN, "w");
    if (fo == NULL)
    {
        perror ("error occurred while opening a file");
        return 0;
    }
    int nseqs = 0;

    double *InnerProd = new double[NSVseqs + 5];

    while (!feof(sfi))
    {
        printf("seq id = %d ", nseqs);

        Sequence *sgi = new Sequence(maxseqlen+3);
        sgi->readFasta(sfi);

        if(sgi->getLength()>0)
        {
            seqname[nseqs] = new char [strlen(sgi->getName()) + 1];
            sprintf(seqname[nseqs], sgi->getName());

            double cur_norm = 0;
            cur_norm = calcnorm(sgi, G, L, K, onlyMiddle);
            printf(" %.2lf\n", cur_norm);

            for (int i = 0; i < NSVseqs; i++)
                InnerProd[i] = 0;
            SVTree->calInnerProd(sgi, G, L, K, InnerProd, onlyMiddle);

            for (int i = 0; i < 10; i++)
                printf("%.2lf ", InnerProd[i] / (cur_norm * SV_norm[i]));
            printf("\n");

            double svmScore = 0;
            for (int i = 0; i < NSVseqs; i++)
            {
                if (SV_norm[i] * cur_norm < 1e-5)
                    continue;
                svmScore += alpha[i] * InnerProd[i] / (SV_norm[i] * cur_norm);
            }
            fprintf(fo, "%s\t%f\n",seqname[nseqs], svmScore);
            nseqs++;
        }

    }

    fclose(fo);

    delete []SV_norm;
    delete []alpha;
    delete []InnerProd;

    delete []seqname;

//	delete sgi;
    delete SVTree;
    return 0;
}


double calcnorm(Sequence *sgi, int g, int l, int k, bool onlyMiddle)
{
    KTree *TreeS = new KTree();
    TreeS->addSequence(0, sgi->getSeqBaseId(), sgi->getLength(), g, l ,k, onlyMiddle);
    return TreeS->calcNorm(sgi, g, l, k, onlyMiddle);
}
/*
int main()
{
    OptsSVMClassify opt;

    opt.G = 2;
    opt.L = 3;
    opt.K = 6;
    opt.maxseqlen = DEF_MAXSEQLEN;
    opt.maxnumseq = DEF_MAXNUMSEQ;
    opt.addRC = false;

    opt.seqfile = "neg.fa";
    opt.svseqfile = "svm_train_svseq.fa";
    opt.alphafile = "svm_train_svalpha.out";
    opt.outfile = "res.out";

    svmClassifySuffixTree(opt);


}
*/
int main(int argc, char** argv) // mainSVMclassify
{
    OptsSVMClassify opt;
    int c;

    opt.G = DEF_G;
    opt.L = DEF_L;
    opt.K = DEF_K;
    opt.maxseqlen = DEF_MAXSEQLEN;
    opt.maxnumseq = DEF_MAXNUMSEQ;
    opt.addRC = false;
    opt.onlyMiddle = true;

    if (argc == 1) { print_usage_and_exit(argv[0]); }

    while ((c = getopt (argc, argv, "g:l:k:m:n:R:A")) != -1)
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

    if (argc-optind != 4) { print_usage_and_exit(argv[0]); }

    int index = optind;
    opt.seqfile = argv[index++];
    opt.svseqfile = argv[index++];
    opt.alphafile = argv[index++];
    opt.outfile = argv[index++];

    svmClassifySuffixTree(opt);

    return 0;
}
