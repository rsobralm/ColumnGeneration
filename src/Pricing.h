#ifndef PRICING_H
#define PRICING_H

#include <ilcplex/ilocplex.h>
#include "Data.h"

class Pricing
{

public:
    IloEnv env;
    IloNumArray pi; // variáveis duais
    Data *data;
    IloCplex pricing_problem;
    IloBoolVarArray x;

    Pricing(Data *data, IloEnv env, IloNumArray pi);
    void buildPricingProblem();
    void solvePricingProblem();
    
    ~Pricing();
};




#endif