//
// Created by cyhong on 4/22/2017.
//


//#include "global.h"

#ifndef GLK_CCONVERTER_H
#define GLK_CCONVERTER_H

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

class CConverter
{
public:

    CConverter();
//	CConverter(char *alphabet=NULL,int alphabetSize=0);
    void init();

    virtual ~CConverter();

    int cidx[256];
    char *icidx; // 'ACGT'= 0123
    char *icidxL;  // 'acgt'= 0123
    char bcompl[256];

    char alphabet[256];
    int b; //alphabetSize;

};

#endif