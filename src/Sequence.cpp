//
// Created by cyhong on 3/8/2017.
//

#include <cstring>
#include "Sequence.h"

Sequence::Sequence(int maxLength, Sequence *sCopyFrom)
{
    static int serialnumber = 0;

    this->maxLength = maxLength;
    seqName = new char[MAX_LINE_WIDTH];
    seq = new char[maxLength];
    seqBaseId = new int[maxLength];
    seqBaseIdRC = new int[maxLength];
    length = 0;

    sprintf(seqName,"seq_%d",serialnumber++);
    sprintf(seq,"");

    if (sCopyFrom != NULL)
    {
        Sequence *s = sCopyFrom;
        length = s->getLength();
        sprintf(seqName, "%s", s->getName());
        int i;
        char* sseq = s->getSeq();
        int* sseqBaseId = s->getSeqBaseId();
        int* sseqBaseIdRC = s->getSeqBaseIdRC();
        for (i = 0; i < length; i++)
        {
            seq[i] = sseq[i];
            seqBaseId[i] = sseqBaseId[i];
            seqBaseIdRC[i] = sseqBaseIdRC[i];
        }

    }
}
Sequence::~Sequence() {
    freeMem(seq);
    freeMem(seqName);
    freeMem(seqBaseId);
    freeMem(seqBaseIdRC);
    length = 0;
}

char *Sequence::getSeq()
{
    return this->seq;
}

char *Sequence::getName()
{
    return this->seqName;
}

int *Sequence::getSeqBaseId()
{
    return this->seqBaseId;
}

int *Sequence::getSeqBaseIdRC()
{
    return this->seqBaseIdRC;
}

int Sequence::getLength()
{
    return this->length;
}
double Sequence::getWeight()
{
    return this->weight;
}

void Sequence::setWeight(double w)
{
    this->weight=w;
}

char* Sequence::getNameLink()
{
    return this->NameLink;
}

void Sequence::setNameLink(char *sp)
{
    this->NameLink=sp;
}
int Sequence::readFasta(FILE *f)
{
    length = 0;
    static char nextName[MAX_LINE_WIDTH];
    static int hasNextName = 0;

    static char sline[MAX_LINE_WIDTH+3];

    if (f == NULL)
    {
        return 0;
    }

    fgets(sline, MAX_LINE_WIDTH, f);
    if (sline[0] == '>') //seq_name line
    {
        sscanf(sline+1, "%s", nextName);
        fgets(sline, MAX_LINE_WIDTH, f);
        hasNextName = 1;
    }
    sprintf(seqName,"%s",(hasNextName ? nextName : "NA"));
    hasNextName = 0;

    while(true)
    {
        if (feof(f) || (sline[0] =='>'))
        {
            break;
        }

        int i=0;
        if (sline[0] != ';') // ignore comment lines starting with ";"
        {
            while (sline[i] != 0)
            {
				if ((toupper(sline[i]) == 'A') ||
					(toupper(sline[i]) == 'C') ||
					(toupper(sline[i]) == 'G') ||
					(toupper(sline[i]) == 'T'))
                {
                    this->seq[this->length] = sline[i];
                    length++;
                }
                i++;
            }
            //printf("\n---%s\n",this->seq);
        }
        fgets(sline, MAX_LINE_WIDTH, f);
    }
    if (sline[0]=='>') //first line
    {
        sscanf(sline + 1, "%s", nextName);
        hasNextName = 1;
    }

    seq[length] = 0;

    int i;
    int *cidx = globalConverter.cidx;
    for(i = 0; i < length; i++)
    {
        this->seqBaseId[i] = cidx[seq[i]];
    }

    for (int i = 0; i < length; i++)
    {
        this->seqBaseIdRC[i] = 3 - this->seqBaseId[length - i - 1];
    }
    return length;
}