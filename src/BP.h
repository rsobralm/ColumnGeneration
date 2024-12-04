#ifndef BP_H
#define BP_H

#include <iostream>
#include <ilcplex/ilocplex.h>
#include <cmath>
#include "Data.h"
#include "Master.h"
#include "Pricing.h"
#include "Node.h"
#include <list>
#include <sys/resource.h>

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

    double best_integer = std::numeric_limits<double>::infinity();

    IloEnv env;
    Master master;

    BP(Data *data);
    ~BP();
    std::pair<int, int> columnGeneration(Node *node);
    std::pair<int, int> getFractionalPair(double &most_fractional, std::vector<std::vector<double>> &z, IloNumArray &lambda_values, std::vector<std::vector<bool>> &lambda_items);
    void setBoundsAndAddConstraints(Master *master, Pricing *pricing, Node *node);
    void BranchAndPrice();
    void prune(IloNumVarArray &lambda);
};




#endif