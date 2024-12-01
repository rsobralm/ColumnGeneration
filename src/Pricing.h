#ifndef PRICING_H
#define PRICING_H

#include <ilcplex/ilocplex.h>
#include "Data.h"
#include <set>

class Pricing
{

public:
    IloEnv env;
    IloNumArray pi; // variáveis duais
    Data *data;
    IloCplex pricing_problem;
    IloBoolVarArray x;
    IloModel pricing_model;
    IloObjective objective_function;

    Pricing(Data *data, IloEnv env, IloNumArray pi);
    void buildPricingProblem();
    void setObjectiveFunction(IloNumArray pi);
    void solvePricingProblem();
    void addBranchingConstraints(std::set<std::pair<int, int>> &together, std::set<std::pair<int, int>> &separated);
    
    ~Pricing();
};




#endif