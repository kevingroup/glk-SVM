/*
	Created by cyhong on 4/23/2017.
	This code is partially adopted from gkm-SVM(Mahmoud Ghandi, http://www.beerlab.org/gkmsvm/)
*/

#include "TreeLeafData.h"

TreeLeafData::TreeLeafData(void)
{
    this->n = 0;
    this->seqInfo = NULL;
}

TreeLeafData::~TreeLeafData(void) {
    delete []seqInfo;
}

void TreeLeafData::add(int seqID)
{
    if(n == 0)
    {
        n = 1;
        this->seqInfo = new int[2 * 2];
        this->seqInfo[0]= seqID;
        this->seqInfo[1]= 1;
    }else{
        if (n == 1)
        {
            if (this->seqInfo[0] == seqID)   //seqId already existed
            {
                this->seqInfo[1]++;
            }
            else
            {
                this->seqInfo[2] = seqID;
                this->seqInfo[3] = 1;
                n = 2;
            }
        }
        else // n>1
        {
            if (this->seqInfo[2 * n - 2] == seqID) //seqId already existed
            {
                this->seqInfo[2 * n - 1]++;
            }
            else {
                if ((n & (n - 1)) == 0) // n is power of 2
                {
                    // expand memory
                    int *newseqInfo;
                    newseqInfo = new int[n << 2];
                    for (int i = 0; i < (2 * n); i++) {
                        newseqInfo[i] = this->seqInfo[i];
                    }
                    delete[] this->seqInfo;
                    this->seqInfo = newseqInfo;
                }
                this->seqInfo[2 * n] = seqID;
                this->seqInfo[2 * n + 1] = 1;
                n++;
            }
        }
    }
}

