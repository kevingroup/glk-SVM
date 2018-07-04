//
// Created by cyhong on 3/8/2017.
//
#ifndef GLK_SEQUENCE_H
#define GLK_SEQUENCE_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "global.h"

class Sequence{
private:
    char *seq;
    char *seqName;
    int *seqBaseId; // 0123 instead of acgt
    int *seqBaseIdRC;
    int length;
    int maxLength;
    double weight;
    char *NameLink;

public:

    Sequence(int maxLength, Sequence *sCopyFrom = NULL );  // makes a replicate of s;
    char *getSeq();
    int getLength();

    //void setName(char *newname);
    char *getName();

    int *getSeqBaseId();
    int *getSeqBaseIdRC();
    double getWeight();
    void setWeight(double w);
    char *getNameLink();
    void setNameLink(char *sp);

    int readFasta(FILE *f);

    virtual ~Sequence();

};

#endif