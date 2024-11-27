#include "BP.h"


BP::BP(Data *data) : data(data)
{
    
}

BP::~BP()
{
}

std::vector<std::vector<double>> BP::getMatrixZ(IloNumArray lambda_values, std::vector<std::vector<bool>> lambda_items){

	std::vector<std::vector<double>> z = std::vector<std::vector<double>> (data->n_items, std::vector<double>(data->n_items, 0));

	for (int i = 0; i < data->n_items; i++){
        for (int j = i + 1; j < data->n_items; j++){
            for (int k = data->n_items; k < lambda_items.size(); k++){
                if (lambda_items[k][i] && lambda_items[k][j]){
                    z[i][j] += lambda_values[k];
                    z[j][i] += lambda_values[k];
                }
            }
		}
	}

	return z;
}

std::pair<int, int> getFractionalPair(std::vector<std::vector<double>> z, std::vector<std::vector<bool>> lambda_items){
	// Find the most fractional pair of items

	double most_fractional = std::numeric_limits<double>::infinity();
	std::pair<int, int> fractional_pair = std::make_pair(-1, -1);

	int z_size = z.size();
	double dist;

	for (int i = 0; i < z_size; i++){
		for (int j = i+1; j < z_size; j++){
			dist = fabs(z[i][j] - 0.5);
			if (dist < most_fractional){
				most_fractional = dist;
				fractional_pair = std::make_pair(i, j);
			}
		}
	}
	return fractional_pair;
}

void BP::columnGeneration(Node &node){

    IloEnv env;

    Master master(data, env);
    master.buildMasterProblem();
    master.solveMasterProblem();
    //Pricing pricing(data, env);

    int lambda_counter = data->n_items;

    while (true){
        IloNumArray pi(env, data->n_items);
        master.rmp.getDuals(pi, master.partition_constraint);
        Pricing pricing(data, env, pi);
        pricing.buildPricingProblem();
        pricing.solvePricingProblem();

        if (pricing.pricing_problem.getObjValue() < -1e-5){

			std::cout << "Reduced cost is equal to " << pricing.pricing_problem.getObjValue() << ", which is less than 0..." << std::endl;

			IloNumArray entering_col(env, data->n_items);

			pricing.pricing_problem.getValues(pricing.x, entering_col);

			std::cout << std::endl << "Entering column:" << std::endl;
			for (size_t i = 0; i < data->n_items; i++)
			{
				std::cout << (entering_col[i] < 0.5 ? 0 : 1) << std::endl;
			}
			std::cout << std::endl;

			std::vector<bool> items(data->n_items, false);
            for (int i = 0; i < pricing.x.getSize(); i++)
            {
                if (entering_col[i] > 0.5)
                {
                    items[i] = true;
                    // std::cout << i << " ";
                }
            }

			master.lambda_items.push_back(items);


			// Add the column to the master problem
			// (the cost of the new variable is always 1)
			char var_name[50];
			sprintf(var_name, "y%d", lambda_counter++);
			IloNumVar new_lambda(master.master_objective(1) + master.partition_constraint(entering_col), 0, IloInfinity);
			new_lambda.setName(var_name);

			master.lambda.add(new_lambda);

			std::cout << "Solving the RMP again..." << std::endl;

			// ...
			master.rmp.solve();
		}
		else
		{
			std::cout << "No column with negative reduced costs found. The current basis is optimal" << std::endl;
			// cout << "Final master problem: " << endl;
			// system("cat model.lp");
			break;
		}
    }

    master.rmp.solve();
	int n_lambdas = 0;

	for (size_t j = 0; j < master.lambda.getSize(); j++)
	{	
		if (master.rmp.getValue(master.lambda[j]) > 0.5)
			n_lambdas++;
			//cout << <<rmp.getValue(lambda[j]) << " ";
	}
	std::cout << std::endl;
	env.end();

	std::cout << n_lambdas << std::endl;
    
}

void BP::addConstraintItemsTogether(Master *master, Pricing *pricing, std::vector<std::pair<int, int>> &together, std::vector<std::vector<bool>> &lambdaItens){
    
    for (auto &p : together){

        pricing->pricing_model.add(pricing->x[p.first] == pricing->x[p.second]);

        for (int i = data->n_items; i < lambdaItens.size(); i++){
           
           // None of the items are in the lambda
            if (lambdaItens[i][p.first] == false && lambdaItens[i][p.second] == false){
                continue;
            }
            // Both items are in the lambda
            if (lambdaItens[i][p.first] == true && lambdaItens[i][p.second] == true){
                continue;
            }
            // Only one of the items is in the lambda, so we set lambda[i] = 0
            master->lambda[i].setUB(0.0);
        }
    }
}

void BP::addConstraintItemsSeparated(Master *master, Pricing *pricing, std::vector<std::pair<int, int>> &separated, std::vector<std::vector<bool>> &lambdaItens){
    
    for (auto &p : separated){

        pricing->pricing_model.add(pricing->x[p.first] + pricing->x[p.second] <= 1);

        for (int i = data->n_items; i < lambdaItens.size(); i++){
            
			// Both items are in the lambda
            if (lambdaItens[i][p.first] == true && lambdaItens[i][p.second] == true){
				master->lambda[i].setUB(0.0);
            }

        }
    }
}


void BP::BranchAndPrice(){

	Node root(data);
	std::vector<Node> tree;

	tree.push_back(root);
	tree[0].LB = 0;
	tree[0].UB = data->n_items;
	tree[0].separated = std::vector<std::pair<int, int>>(data->n_items, std::make_pair(-1, -1));
	tree[0].merged = std::vector<std::pair<int, int>>(data->n_items, std::make_pair(-1, -1));

	columnGeneration();




	while(!tree.empty()){

		// roda a geração de colunas pro nó atual
		// se a solução for inteira, atualiza o UB, se for ótima, encerra.
		// se for fracionário e o LB for pior que a melhor solução inteira podar
		// se for fracionário e o LB for melhor que a melhor solução inteira, ramifica


	}

}

