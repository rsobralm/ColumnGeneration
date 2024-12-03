#include "BP.h"



BP::BP(Data *data) : data(data)
{	

	this->env = IloEnv();
	this->master = Master(data, env);
    master.buildMasterProblem();

}

BP::~BP()
{
}

void BP::setMatrixZ(std::vector<std::vector<double>> &z, IloNumArray &lambda_values, std::vector<std::vector<bool>> &lambda_items){

	for (int i = 0; i < data->n_items; i++){
        for (int j = i + 1; j < data->n_items; j++){
            for (int k = data->n_items; k < lambda_items.size(); k++){
                if (lambda_items[k][i] and lambda_items[k][j]){
                    z[i][j] += lambda_values[k];
                    z[j][i] += lambda_values[k];
                }
            }
		}
	}
}

std::pair<int, int> getFractionalPair(std::vector<std::vector<double>> &z, std::vector<std::vector<bool>> &lambda_items){
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


// void BP::prune(IloNumVarArray &lambda){
// 	 for (int i = 0; i < lambda.getSize(); i++){
//         lambda[i].setUB(IloInfinity);
//     }
// }

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

	//master.buildMasterProblem();

    
	IloEnv env2;

	master.setBounds(node, master.lambda_items);

	// std::cout << "Building master problem " << master.lambda_items.size() << std::endl;
	// for (int i = 0; i < master.lambda_items.size(); i++){
	// 	for (int j = 0; j < master.lambda_items[i].size(); j++){
	// 		std::cout << master.lambda_items[i][j] << " ";
	// 	}
	// 	std::cout << master.lambda[i] <<std::endl;
	// }
	// std::cin.get();

	

	// for (int i = 0; i < master.lambda.getSize(); i++){
	// 	std::cout << i << " " << master.lambda[i] << "\n";
	// }

    master.solveMasterProblem();
	//master.rmp.solve();

	// std::cout << "n_cols: " << master.rmp.getNcols()  << std::endl;
	// std::cout << "n_rows: " << master.rmp.getNrows()  << std::endl;
	// std::cout << "Obj value: " << master.rmp.getObjValue() << std::endl;


	// IloNumArray pi(env2, data->n_items);
	// Pricing pricing(data, env2, pi);
    // pricing.buildPricingProblem();

	// pricing.addBranchingConstraints(node->separated, node->merged);
	IloNumArray pi(env2, data->n_items);
	Pricing pricing(data, env2, pi);
    pricing.buildPricingProblem();
	pricing.addBranchingConstraints(node->merged, node->separated);


	std::vector<std::vector<double>> z = std::vector<std::vector<double>> (data->n_items, std::vector<double>(data->n_items, 0));



    while (true){
		
		if (master.rmp.getCplexStatus() == IloCplex::Infeasible){
            break;
        }

		IloNumArray pi(env2, data->n_items);
        master.rmp.getDuals(pi, master.partition_constraint);
		pricing.setObjectiveFunction(pi);

		
        pricing.solvePricingProblem();

		//pricing.pricing_problem.exportModel("pricing.lp");


		pi.clear();
        pi.end();

		// If pricing is infeasible prune the node

		if (pricing.pricing_problem.getCplexStatus() == IloCplex::Infeasible){
			prune(master.lambda);

            pricing.pricing_model.end();

            env2.end();

            master.rmp.clear();
            master.rmp.end();
			return std::make_pair(-1, -1);
		}

		//std::cout << "Reduced cost is equal to " << pricing.pricing_problem.getObjValue()  << std::endl;

		// If reduced cost is negative, add the column to the master problem
        if (pricing.pricing_problem.getObjValue() < -EPS){

			//std::cout << "Reduced cost is equal to " << pricing.pricing_problem.getObjValue() << ", which is less than 0..." << std::endl;

			IloNumArray entering_col(env2, data->n_items);

			pricing.pricing_problem.getValues(entering_col, pricing.x);
			

			// std::cout << std::endl << "Entering column:" << std::endl;
			// for (size_t i = 0; i < data->n_items; i++)
			// {
			// 	std::cout << (entering_col[i] < 0.5 ? 0 : 1) << std::endl;
			// }
			// std::cout << std::endl;

			


			// Add the column to the master problem
			// (the cost of the new variable is always 1)
			IloNumVar new_lambda(master.master_objective(1) + master.partition_constraint(entering_col), 0, IloInfinity);
			char var_name[50];
			sprintf(var_name, "y%d", (int)master.lambda.getSize());
			new_lambda.setName(var_name);

			master.lambda.add(new_lambda);

			std::vector<bool> items(data->n_items, false);
            for (int i = 0; i < entering_col.getSize(); i++)
            {
                if (entering_col[i] > 1 - EPS)
                {
                    items[i] = true;
                    // std::cout << i << " ";
                }
            }

			master.lambda_items.push_back(items);

			//std::cout << "Solving the RMP again..." << std::endl;

			// ...
			master.rmp.solve();

			entering_col.clear();
            entering_col.end();
            pricing.pricing_problem.clear();
            pricing.pricing_problem.end();
		}
		else
		{
			//std::cout << "No column with negative reduced costs found. The current basis is optimal" << std::endl;
			// cout << "Final master problem: " << endl;
			// system("cat model.lp");
			//pricing.pricing_model.end();

			// IloBoolVarArray x_values(env2, pricing.x.getSize());
			// pricing.pricing_problem.getValues(x_values, pricing.x);

			// for (int i = 0; i < pricing.x.getSize(); i++)
			// {
			// 	std::cout << "x[" << i << "] " << pricing.pricing_problem.getValue(pricing.x[i]) << "\n";
			// }


			pricing.pricing_problem.clear();
            pricing.pricing_problem.end();
			break;
		}

		//pricing.pricing_model.end();
    }


	if (master.rmp.getCplexStatus() == IloCplex::Infeasible){
		prune(master.lambda);
		return std::make_pair(-1, -1);
	}

	std::cout << "LB: " << master.rmp.getObjValue() << std::endl;
	std::cout << "Ncols: " << master.rmp.getNcols() << std::endl;

	IloNumArray lambda_values(env2, master.lambda.getSize());
	master.rmp.getValues(lambda_values, master.lambda);

	// Prune if there's a artificial variable with a positive value
	for (int i = 0; i < data->n_items; i++){
		if (lambda_values[i] > EPS){
			prune(master.lambda);

			//pricing.pricing_model.end();
			//env2.end();
			//lambda_values.end();
			master.rmp.clear();
			master.rmp.end();
			//std::cout << "artificial variable" << std::endl;

			return  std::make_pair(-1, -1);
		}
    }

	// Prune if the LB is greater than the best integer solution
	if (std::ceil(master.rmp.getObjValue() - EPS) - UB >= 0){
		prune(master.lambda);
		//std::cout << "LB maior que UB" << std::endl;

		//pricing.pricing_model.end();
		//env2.end();
		//lambda_values.end();
		master.rmp.clear();
		master.rmp.end();
		return std::make_pair(-1, -1);
	}


	setMatrixZ(z, lambda_values, master.lambda_items);

	std::pair<int, int> fractional_pair = getFractionalPair(z, master.lambda_items);

	// std::cout << fractional_pair.first << " " << fractional_pair.second << std::endl;

	// std::cout << "Objective value: " << master.rmp.getObjValue() << std::endl;
	// std::cout << "Z Value: " << z[fractional_pair.first][fractional_pair.second] << std::endl;


	//z[fractional_pair.first][fractional_pair.second] = 1;
	// If the solution is integer, update the UB and prune

	//std::cout << fabs(z[fractional_pair.first][fractional_pair.second] - std::round(z[fractional_pair.first][fractional_pair.second])) << std::endl;

	//std::cout << "Fractional value: " << fmod(z[fractional_pair.first][fractional_pair.second], 1) << std::endl;

	if (fmod(z[fractional_pair.first][fractional_pair.second], 1) < EPS){
		//std::cout << "Integer solution found: " << master.rmp.getObjValue() << std::endl;
		if (master.rmp.getObjValue() < UB){
			UB = master.rmp.getObjValue();
		}

		prune(master.lambda);

		//pricing.pricing_model.end();
		//env2.end();
		//lambda_values.end();
		master.rmp.clear();
		master.rmp.end();

		return std::make_pair(-1, -1);
	}

	
	// if (std::abs(0.5 - z[fractional_pair.first][fractional_pair.second]) < EPS){

	// 	std::cout << "Integer solution found: " << master.rmp.getObjValue() << std::endl;
    //     if (master.rmp.getObjValue() < UB){
    //         UB = master.rmp.getObjValue();
    //     }

    //     prune(master.lambda);

    //     return std::make_pair(-1, -1);

	lambda_values.end();

	for (int i = 0; i < master.lambda.getSize(); i++)
    {
        master.lambda[i].setUB(IloInfinity);
    }

	master.rmp.clear();
    master.rmp.end();
    // }


	//std::cout << "Most fractional pair: " << fractional_pair.first << " " << fractional_pair.second << std::endl;
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

	Node current = tree.back();
	std::pair<int, int> p = columnGeneration(&current);
	// while(!tree.empty()){

	// 	Node current = tree.back();
	// 	std::pair<int, int> p = columnGeneration(&current);

	// 	//std::cin.get();

	// 	if (p.first == -1 && p.second == -1){
	// 		tree.pop_back();
	// 		continue;
	// 	}

	// 	Node left(data);
	// 	Node right(data);

	// 	left.LB = current.LB;
	// 	left.UB = current.UB;
	// 	left.separated = current.separated;
	// 	left.merged = current.merged;

	// 	right.LB = current.LB;
	// 	right.UB = current.UB;
	// 	right.separated = current.separated;
	// 	right.merged = current.merged;

	// 	left.merged.insert(p);
	// 	right.separated.insert(p);

	// 	tree.push_back(left);
	// 	tree.push_back(right);
	// 	//tree.push_back(left);

	// 	tree.pop_back();
	// }

	

	std::cout << "Best solution found: " << UB << std::endl;

}

