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
    Master(Data *data, IloEnv env);
    void buildMasterProblem();
    void solveMasterProblem(IloCplex master_model);
    ~Master();
};

Master::Master(Data *data, IloEnv env) : data(data)
{   
    this->env = env;
    this->master_model = IloModel(env);
}


Master::~Master()
{
}

void Master::buildMasterProblem(){

	IloNumVarArray lambda(env, data->n_items, 0, IloInfinity);

	IloExpr sum_obj(env);
	IloRangeArray partition_constraint(env);

	for (int i = 0; i < data->n_items; i++)
	{
		char var_name[50];
		sprintf(var_name, "y%d", i);

		lambda[i].setName(var_name);
		sum_obj += BIG_M * lambda[i];

		partition_constraint.add(lambda[i] == 1);
	}

	master_model.add(partition_constraint);

	IloObjective master_objective = IloMinimize(env, sum_obj);
	master_model.add(master_objective);
}

void Master::solveMasterProblem(IloCplex master_model){
    IloCplex rmp(master_model);

    rmp.setOut(env.getNullStream()); // disables CPLEX log

    rmp.solve();

    std::cout << "Initial lower bound: " << rmp.getObjValue() << std::endl;
    
}

#endif