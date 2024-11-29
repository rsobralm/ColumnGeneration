#ifndef MASTER_H
#define MASTER_H

#include <ilcplex/ilocplex.h>
#include <iostream>
#include <vector>
#include "Data.h"
#include "Node.h"

#define BIG_M 1e6

class Master
{
public:
    Data *data;
    IloModel master_model;
    IloEnv env;
	IloObjective master_objective;
	IloRangeArray partition_constraint;
	IloNumVarArray lambda;
	IloCplex rmp; // restricted master problem
    std::vector<std::vector<bool>> lambda_items; // 1 if the item is in the lambda, 0 otherwise
    Master(Data *data, IloEnv env);
    void buildMasterProblem();
    void solveMasterProblem();
    void setBounds(Node *node, std::vector<std::vector<bool>> &lambdaItens);
    ~Master();
};



#endif