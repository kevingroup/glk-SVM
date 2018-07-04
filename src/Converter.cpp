//
// Created by cyhong on 4/22/2017.
//

#include "Converter.h"
#include "global.h"
#include <cstring>

CConverter::CConverter()
{
    b = 4;
    alphabet[0] = 'A';
    alphabet[1] = 'C';
    alphabet[2] = 'G';
    alphabet[3] = 'T';

    init();

}

void CConverter::init(){

    int ici;
    memset(cidx, 0, sizeof(cidx));

    //	cidx['a'] = 0; 		cidx['c'] = 1; 		cidx['g'] = 2; 		cidx['t'] = 3;
    //	cidx['A'] = 0; 		cidx['C'] = 1; 		cidx['G'] = 2; 		cidx['T'] = 3;
    for(ici = 0; ici < b; ici++)
    {
        cidx[toupper(alphabet[ici])] = ici;
        cidx[tolower(alphabet[ici])] = ici;
    }

    icidx = new char[b];
    icidxL = new char[b];

    //	icidx[0]  = 'A';		icidx[1]  = 'C';		icidx[2]  = 'G';		icidx[3]  = 'T';
    //	icidxL[0] = 'a';		icidxL[1] = 'c';		icidxL[2] = 'g';		icidxL[3] = 't';
    for(ici = 0; ici < b; ici++)
    {
        icidx[ici]=  toupper(alphabet[ici]);
        icidxL[ici] = tolower(alphabet[ici]);
    }


}
CConverter::~CConverter()
{
    delete []icidx;
    delete []icidxL;
}
