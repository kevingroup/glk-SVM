//
// Created by cyhong on 10/16/2017.
//

#include "SequenceNames.h"
CSequenceNames::CSequenceNames(void)
{
    Nseqs=0;

    seqf = NULL;

    nSeqsRead=0;
    nextSeqtoRead=0;
    curSeq=NULL;
}


CSequenceNames::~CSequenceNames(void)
{
    int i;
    if (Nseqs!=0)
    {
        for(i=0;i<Nseqs;i++)
        {
            delete seqNames[i];
        }
        Nseqs = 0;
    }

    if (seqf!=NULL)
    {
        fclose(seqf);
        seqf = NULL;
    }

    if (this->curSeq!=NULL) delete curSeq;

}

int CSequenceNames::readSeqNames(char *seqNamesFN)
{
    int i,j;
    if (Nseqs!=0)
    {
        for(i=0;i<Nseqs;i++)
        {
            delete seqNames[i];
        }
        Nseqs = 0;
    }

    char stmp[10000];

    FILE *f = fopen(seqNamesFN, "r") ;
    while (!feof(f))
    {
        if (fgets(stmp, 10000-5, f))
        {
            if (stmp[0]!=0) {
                seqNames[Nseqs] = new char[MAXSeqnameLENGTH];
                sscanf(stmp, "%s", seqNames[Nseqs]);
                Nseqs++;
            }
        }
    }
    fclose(f);
    return Nseqs;
}


int CSequenceNames::readSeqNamesandWeights(char *seqNamesFN)
{
    int i,j;
    if (Nseqs!=0)
    {
        for(i=0;i<Nseqs;i++)
        {
            delete seqNames[i];
        }
        Nseqs = 0;
    }

    char stmp[10000];

    FILE *f = fopen(seqNamesFN, "r") ;
    while (!feof(f))
    {
        if (fgets(stmp, 10000-5, f))
        {
            if (stmp[0]!=0) {

                seqNames[Nseqs] = new char[MAXSeqnameLENGTH];
                sscanf(stmp, "%s%lf", seqNames[Nseqs], weight+Nseqs);
                //printf("%s\n%e\n",seqNames[Nseqs],weight[Nseqs]);
                Nseqs++;
            }
        }

    }
    fclose(f);
    return Nseqs;
}

void CSequenceNames::openSeqFile( char *seqFN,  int maxSeqLength)
{
    this->seqf = fopen(seqFN, "r");

    if (this->curSeq!=NULL) delete curSeq;

    this->curSeq = new Sequence(maxSeqLength);
}

Sequence *CSequenceNames::nextSeq()
{
    while (!feof(seqf))
    {
        if (this->nextSeqtoRead ==0)
        {
            curSeq->readFasta(this->seqf);  // read a new sequence
        }

        while (nextSeqtoRead<this->Nseqs)
        {
            if (stringcompare(this->seqNames[nextSeqtoRead], curSeq->getName(), MAXSeqnameLENGTH)) 			//if (strcmp(this->seqNames[nextSeqtoRead], curSeq->getName())==0)
            {
                curSeq->setWeight(weight[nextSeqtoRead]);
                curSeq->setNameLink(seqNames[nextSeqtoRead]);

                nextSeqtoRead++;
                nSeqsRead++;

                if (nSeqsRead==Nseqs)
                {
                    fclose(seqf);
                    seqf = NULL;
                }

                return curSeq;
            }
            nextSeqtoRead++;
        }
        nextSeqtoRead = 0;
    }

    fclose(seqf);
    seqf = NULL;
    return NULL;
}
