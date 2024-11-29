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
    IloModel pricing_model;

    Pricing(Data *data, IloEnv env, IloNumArray pi);
    void buildPricingProblem();
    void solvePricingProblem();
    void addBranchingConstraints(std::vector<std::pair<int, int>> &together, std::vector<std::pair<int, int>> &separated);
    
    ~Pricing();
};




#endif