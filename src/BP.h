#ifndef BP_H
#define BP_H

#include <iostream>
#include "Data.h"
#include "Master.h"
#include "Pricing.h"



class BP
{

public:
    Data *data;
    int capacity;
    int n;

    BP(Data *data);
    ~BP();
    void columnGeneration();
};




#endif