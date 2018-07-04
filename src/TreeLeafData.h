//
// Created by cyhong on 4/23/2017.
//

#ifndef GLK_TREELEAFDATA_H
#define GLK_TREELEAFDATA_H

#include <stdio.h>

class TreeLeafData {
public:

    int n; // number of GK-mers in the list
    int* seqInfo;

    TreeLeafData(void);
    ~TreeLeafData(void);

    void add(int seqID);

};


#endif //GLK_TREELEAFDATA_H
