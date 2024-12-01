#ifndef BP_H
#define BP_H

#include <iostream>
#include <ilcplex/ilocplex.h>
#include <cmath>
#include "Data.h"
#include "Master.h"
#include "Pricing.h"
#include "Node.h"

#define EPS 1e-6



typedef struct
{
    IloAlgorithm::Status status;
    double cost;
} ResultsCG;

class BP
{

public:
    Data *data;
    int capacity;
    int n;

    double UB = std::numeric_limits<double>::infinity();

    IloEnv env;
    Master master;

    std::vector<std::vector<bool>> lambdaItens;

    BP(Data *data);
    ~BP();
    std::pair<int, int> columnGeneration(Node *node);
    void setMatrixZ(std::vector<std::vector<double>> &z, IloNumArray &lambda_values, std::vector<std::vector<bool>> &lambda_items);
    void addConstraintItemsTogether(Master *master, Pricing *pricing, std::vector<std::pair<int, int>> &together, std::vector<std::vector<bool>> &lambdaItens);
    void addConstraintItemsSeparated(Master *master, Pricing *pricing, std::vector<std::pair<int, int>> &together, std::vector<std::vector<bool>> &lambdaItens);
    void BranchAndPrice();
    void prune(IloNumVarArray &lambda);
    bool checkIfIntegerSolution(IloNumArray lambda_values);
};




#endif