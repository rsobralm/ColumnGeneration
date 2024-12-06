#ifndef PRICING_H
#define PRICING_H

#include <ilcplex/ilocplex.h>
#include <set>
#include "Data.h"


#define BIG_M 1e6
#define EPS 1e-6


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
    void setObjectiveFunction(IloNumArray &pi);
    void solvePricingProblem();   
    double solveCombo(IloNumArray &pi, IloNumArray &entering_col);
    ~Pricing();
};




#endif