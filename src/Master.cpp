#include "Master.h"

Master::Master(Data *data, IloEnv env) : data(data)
{   
    this->env = env;
    this->master_model = IloModel(env);
	this->lambda_items = std::vector<std::vector<bool>>(data->n_items, std::vector<bool>(data->n_items, false));

	for (int i = 0; i < data->n_items; i++){
		lambda_items[i][i] = true;
	}
}

Master::Master(){}


Master::~Master()
{
}

void Master::buildMasterProblem(){

	//IloNumVarArray lambda(env, data->n_items, 0, IloInfinity);

    lambda = IloNumVarArray(env, data->n_items, 0, IloInfinity);


	IloExpr sum_obj(env);
	partition_constraint = IloRangeArray(env);

	for (int i = 0; i < data->n_items; i++)
	{
		char var_name[50];
		sprintf(var_name, "y%d", i);

		lambda[i].setName(var_name);
		sum_obj += BIG_M * lambda[i];

		partition_constraint.add(lambda[i] == 1);
	}

	master_model.add(partition_constraint);

	//IloObjective master_objective = IloMinimize(env, sum_obj);
    master_objective = IloMinimize(env, sum_obj);

	master_model.add(master_objective);
}

void Master::solveMasterProblem(){
    //IloCplex rmp(master_model);

    rmp = IloCplex(master_model);

    rmp.setOut(env.getNullStream()); // disables CPLEX log

    rmp.solve();

    //std::cout << "Initial lower bound: " << rmp.getObjValue() << std::endl;

	// std::cout << "Initial solution: " << std::endl;
	// for (size_t j = 0; j < lambda.getSize(); j++)
	// {
	// 	std::cout << rmp.getValue(lambda[j]) << " ";
	// }
	// std::cout << std::endl;
    
}

void Master::setBounds(Node *node, std::vector<std::vector<bool>> &lambdaItens){
	
	// Set bounds for the lambdas when the items must be together
	for (auto pair : node->merged){
		//std::cout << "Setting bounds for items " << p.first << " and " << p.second << std::endl;

        for (int i = data->n_items; i < lambdaItens.size(); i++){
           
           // None of the items are in the lambda
            if (!lambdaItens[i][pair.first] and !lambdaItens[i][pair.second]){
                continue;
            }
            // Both items are in the lambda
            if (lambdaItens[i][pair.first] and lambdaItens[i][pair.second]){
                continue;
            }
            // Only one of the items is in the lambda, so we set lambda[i] = 0
            lambda[i].setUB(0.0);
        }
    }

	// Set bounds for the lambdas when the items must be separated
	for (auto pair : node->separated){

        for (int i = data->n_items; i < lambdaItens.size(); i++){
            
			// Both items are in the lambda
            if (lambdaItens[i][pair.first] and lambdaItens[i][pair.second]){
				lambda[i].setUB(0.0);
            }

        }
    }
}