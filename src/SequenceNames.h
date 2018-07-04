//
// Created by cyhong on 10/16/2017.
//

#ifndef GLK_SEQUENCENAMES_H
#define GLK_SEQUENCENAMES_H


#include "global.h"
#include "Sequence.h"

#define MAXNSeqs 1000000
#define MAXSeqnameLENGTH 100

class CSequenceNames
{
public:
    CSequenceNames(void);
    ~CSequenceNames(void);
    int Nseqs;
    char *seqNames[MAXNSeqs];

    double weight[MAXNSeqs];

    int readSeqNames(char *seqNamesFN);
    int readSeqNamesandWeights(char *seqNamesFN);

    void openSeqFile( char *seqFN, int maxSeqLength);
    Sequence *nextSeq();

    Sequence *curSeq;
private:
    FILE *seqf;
    int nSeqsRead;  // number of times a sequence was returned by nextSeq()
    int nextSeqtoRead;  // index to the next potential sequence
};


#endif //GLK_SEQUENCENAMES_H
