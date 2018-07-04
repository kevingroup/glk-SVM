//
// Created by cyhong on 4/28/2017.
//

#include "Converter.h"

CConverter globalConverter;

int stringcompare(char *s1, char*s2, int maxlength)
{
    for (int i=0; i<maxlength; i++)
    {
        if (s2[i]!=s1[i]) return 0;
        if (((s1[i]==0)||(s1[i]==13))&&((s2[i]==0)||(s2[i]==13))) return 1;
    }
    return 1;
}

int myrandom(int M) // generates uniform random integer between 0..M-1  (M<<2^31)
{
    int b30r;
    b30r = rand()%1024;
    b30r = (b30r<<10)+(rand()%1024);
    b30r = (b30r<<11)+(rand()%2048);
    return b30r%M;
}

void randomPermute(double *x, int N)
{
    int i,j;
    double h;
    for(i=1; i<N; i++)
    {
        j = myrandom(i+1);
        h = x[i];
        x[i]=x[j];
        x[j]=h;
    }
}

void randomPermute(int *x, int N)
{
    int i,j;
    int h;
    for(i=1; i<N; i++)
    {
        j = myrandom(i+1);
        h = x[i];
        x[i]=x[j];
        x[j]=h;
    }
}

int NchooseK(int n, int k)
{
    if (k > n - k)
        k = n - k;
    int p = 1;
    int q = 1;
    for (int i = 1; i <= k; i++)
    {
        p *= n - i + 1;
        q *= i;
    }
    return p / q;
}