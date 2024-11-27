#ifndef BP_H
#define BP_H

#include <iostream>
#include <ilcplex/ilocplex.h>
#include "Data.h"
#include "Master.h"
#include "Pricing.h"
#include "Node.h"



class BP
{

public:
    Data *data;
    int capacity;
    int n;

    BP(Data *data);
    ~BP();
    void columnGeneration();
    std::vector<std::vector<double>> getMatrixZ(IloNumArray lambda_values, std::vector<std::vector<bool>> lambda_items);
    void addConstraintItemsTogether(Master *master, Pricing *pricing, std::vector<std::pair<int, int>> &together, std::vector<std::vector<bool>> &lambdaItens);
    void addConstraintItemsSeparated(Master *master, Pricing *pricing, std::vector<std::pair<int, int>> &together, std::vector<std::vector<bool>> &lambdaItens);
    void BranchAndPrice();
};




#endif