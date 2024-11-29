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


void BP::prune(IloNumVarArray &lambda){
	 for (int i = 0; i < lambda.getSize(); i++){
        lambda[i].setUB(IloInfinity);
    }
}

std::pair<int, int> BP::columnGeneration(Node *node){


	// for (auto &p : node->separated){
	// 	std::cout << "Separated: " << p.first << " " << p.second << std::endl;
	// }

	// std::cout << node->merged.size() << std::endl;
	// for (auto &p : node->merged){
	// 	std::cout << "Together: " << p.first << " " << p.second << std::endl;
	// }
	// for (int i = 0; i < node->separated.size(); i++){
	// 	std::cout << node->separated[i].first << " " << node->separated[i].second << std::endl;
	// }
	// for (int i = 0; i < node->merged.size(); i++){
	// 	std::cout << node->merged[i].first << " " << node->merged[i].second << std::endl;
	// }

    IloEnv env;

    Master master(data, env);
    master.buildMasterProblem();

	master.setBounds(node, master.lambda_items);

    master.solveMasterProblem();


    int lambda_counter = data->n_items;

    while (true){
		
		if (master.rmp.getCplexStatus() == IloCplex::Infeasible){
            break;
        }


        IloNumArray pi(env, data->n_items);
        master.rmp.getDuals(pi, master.partition_constraint);
        Pricing pricing(data, env, pi);
        pricing.buildPricingProblem();

		pricing.addBranchingConstraints(node->separated, node->merged);
        pricing.solvePricingProblem();




		// If pricing is infeasible prune the node

		if (pricing.pricing_problem.getCplexStatus() == IloCplex::Infeasible){
			prune(master.lambda);
			return std::make_pair(-1, -1);
		}


		// If reduced cost is negative, add the column to the master problem
        if (pricing.pricing_problem.getObjValue() < -1e-5){

			//std::cout << "Reduced cost is equal to " << pricing.pricing_problem.getObjValue() << ", which is less than 0..." << std::endl;

			IloNumArray entering_col(env, data->n_items);

			pricing.pricing_problem.getValues(pricing.x, entering_col);

			// std::cout << std::endl << "Entering column:" << std::endl;
			// for (size_t i = 0; i < data->n_items; i++)
			// {
			// 	std::cout << (entering_col[i] < 0.5 ? 0 : 1) << std::endl;
			// }
			// std::cout << std::endl;

			std::vector<bool> items(data->n_items, false);
            for (int i = 0; i < pricing.x.getSize(); i++)
            {
                if (entering_col[i] >= 0.5)
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

			//std::cout << "Solving the RMP again..." << std::endl;

			// ...
			master.rmp.solve();
		}
		else
		{
			//std::cout << "No column with negative reduced costs found. The current basis is optimal" << std::endl;
			// cout << "Final master problem: " << endl;
			// system("cat model.lp");
			//pricing.pricing_model.end();
			break;
		}

		//pricing.pricing_model.end();
    }


	if (master.rmp.getCplexStatus() == IloCplex::Infeasible){
		prune(master.lambda);
		return std::make_pair(-1, -1);
	}

	IloNumArray lambda_values(env, master.lambda.getSize());
	master.rmp.getValues(lambda_values, master.lambda);

	// Prune if there's a artificial variable with a positive value
	for (int i = 0; i < data->n_items; i++){
		if (lambda_values[i] > EPS){
			prune(master.lambda);
			//std::cout << "artificial variable" << std::endl;

			return  std::make_pair(-1, -1);
		}
    }

	// Prune if the LB is greater than the best integer solution
	if (std::ceil(master.rmp.getObjValue() - EPS) - UB >= 0){
		prune(master.lambda);
		//std::cout << "LB maior que UB" << std::endl;
		return std::make_pair(-1, -1);
	}


	std::vector<std::vector<double>> z = getMatrixZ(lambda_values, master.lambda_items);

	std::pair<int, int> fractional_pair = getFractionalPair(z, master.lambda_items);

	std::cout << "Objective value: " << master.rmp.getObjValue() << std::endl;
	std::cout << "Z Value: " << z[fractional_pair.first][fractional_pair.second] << std::endl;


	// If the solution is integer, update the UB and prune

	if (!fabs(z[fractional_pair.first][fractional_pair.second] - std::round(z[fractional_pair.first][fractional_pair.second])) > EPS){
		std::cout << "Integer solution found: " << master.rmp.getObjValue() << std::endl;
		if (master.rmp.getObjValue() < UB){
			UB = master.rmp.getObjValue();
		}

		prune(master.lambda);

		return std::make_pair(-1, -1);
	}

	
	// if (std::abs(0.5 - z[fractional_pair.first][fractional_pair.second]) < EPS){

	// 	std::cout << "Integer solution found: " << master.rmp.getObjValue() << std::endl;
    //     if (master.rmp.getObjValue() < UB){
    //         UB = master.rmp.getObjValue();
    //     }

    //     prune(master.lambda);

    //     return std::make_pair(-1, -1);
    // }


	std::cout << "Most fractional pair: " << fractional_pair.first << " " << fractional_pair.second << std::endl;
	return fractional_pair;

    // master.rmp.solve();
	// int n_lambdas = 0;

	// for (size_t j = 0; j < master.lambda.getSize(); j++)
	// {	
	// 	if (master.rmp.getValue(master.lambda[j]) > 0.5)
	// 		n_lambdas++;
	// 		//cout << <<rmp.getValue(lambda[j]) << " ";
	// }
	// std::cout << std::endl;
	// env.end();

	// std::cout << n_lambdas << std::endl;
    
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

bool BP::checkIfIntegerSolution(IloNumArray lambda_values){
	for (int i = 0; i < lambda_values.getSize(); i++){
		if(lambda_values[i] > EPS && fabs(lambda_values[i] - (int)lambda_values[i]) > EPS){
			return false;
		}
	}
	return true;
}


void BP::BranchAndPrice(){

	Node root(data);
	std::vector<Node> tree;

	tree.push_back(root);
	tree[0].LB = 0;
	tree[0].UB = data->n_items;
	//tree[0].separated = std::vector<std::pair<int, int>>(data->n_items, std::make_pair(-1, -1));
	//tree[0].merged = std::vector<std::pair<int, int>>(data->n_items, std::make_pair(-1, -1));


	while(!tree.empty()){

		Node current = tree.back();
		std::pair<int, int> p = columnGeneration(&current);

		//std::cin.get();

		if (p.first == -1 && p.second == -1){
			tree.pop_back();
			continue;
		}

		Node left(data);
		Node right(data);

		left.LB = current.LB;
		left.UB = current.UB;
		left.separated = current.separated;
		left.merged = current.merged;

		right.LB = current.LB;
		right.UB = current.UB;
		right.separated = current.separated;
		right.merged = current.merged;

		left.merged.insert(p);
		right.separated.insert(p);

		tree.push_back(left);
		tree.push_back(right);
		//tree.push_back(left);

		// roda a geração de colunas pro nó atual
		// se a solução for inteira, atualiza o UB, se for ótima, encerra.
		// se for fracionário e o LB for pior que a melhor solução inteira podar
		// se for fracionário e o LB for melhor que a melhor solução inteira, ramifica

	}

	std::cout << "Best solution found: " << UB << std::endl;

}

