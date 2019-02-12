//
// Created by cyhong on 3/10/2017.
//

#include "Converter.h"

#ifndef GLK_GLOBAL_H
#define GLK_GLOBAL_H

#define MAX_ALPHABET_SIZE 4
#define MAX_LINE_WIDTH 15000
#define DEF_G 2
#define DEF_L 5
#define DEF_K 8
#define DEF_MAXSEQLEN 15000
#define DEF_MAXNUMSEQ 15000
#define DEF_NUM_OF_THREAD 8

#define MAXNSVSeqs 100000
#define DEF_GLKTYPE 1

#define DEBUG 0 /*0 1*/
#define freeMem(x) if(x!=NULL) delete []x

extern CConverter globalConverter;

int stringcompare(char *s1, char*s2, int maxlength) ;
int myrandom(int M);
void randomPermute(double *x, int N);
void randomPermute(int *x, int N);
int NchooseK(int n, int k);

#endif //GLK_GLOBAL_H

