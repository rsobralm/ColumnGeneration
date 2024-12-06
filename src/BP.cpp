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

void BP::setBoundsAndAddConstraints(Master *master, Pricing *pricing, Node *node){
	// Set bounds for the lambdas in the master problem and add the constraints in the pricing model

	for (auto pair : node->merged){

		// Add the constraint to the pricing model
		pricing->pricing_model.add(pricing->x[pair.first] == pricing->x[pair.second]);

        for (int i = data->n_items; i < master->lambda_items.size(); i++){
           
           // None of the items are in the lambda
            if (!master->lambda_items[i][pair.first] and !master->lambda_items[i][pair.second]){
                continue;
            }
            // Both items are in the lambda
            if (master->lambda_items[i][pair.first] and master->lambda_items[i][pair.second]){
                continue;
            }
            // Only one of the items is in the lambda, so we set lambda[i] = 0
            master->lambda[i].setUB(0.0);
        }
    }

	// Set bounds for the lambdas when the items must be separated
	for (auto pair : node->separated){

		// Add the constraint to the pricing model
		pricing->pricing_model.add(pricing->x[pair.first] + pricing->x[pair.second] <= 1);

        for (int i = data->n_items; i < master->lambda_items.size(); i++){
            
			// Both items are in the lambda
            if (master->lambda_items[i][pair.first] and master->lambda_items[i][pair.second]){
				master->lambda[i].setUB(0.0);
            }

        }
    }

}


std::pair<int, int> BP::getFractionalPair(double &most_fractional, std::vector<std::vector<double>> &z, IloNumArray &lambda_values, std::vector<std::vector<bool>> &lambda_items){
	// Find the most fractional pair of items

	//double most_fractional = std::numeric_limits<double>::infinity();
	double dist;
	std::pair<int, int> fractional_pair = std::make_pair(-1, -1);

	for (int i = 0; i < data->n_items; i++){
        for (int j = i + 1; j < data->n_items; j++){
            for (int k = data->n_items; k < lambda_items.size(); k++){
                if (lambda_items[k][i] and lambda_items[k][j]){
                    z[i][j] += lambda_values[k];
                    z[j][i] += lambda_values[k];
                }
				
            }
			dist = std::abs(0.5 - z[i][j]);
            if (dist < most_fractional)
            {
                fractional_pair = std::make_pair(i, j);
                most_fractional = dist;
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

    IloEnv env2;

		
	IloNumArray pi(env2, data->n_items);
	Pricing pricing(data, env2, pi);

	pricing.buildPricingProblem();
	setBoundsAndAddConstraints(&master, &pricing, node);

	

	std::vector<std::vector<double>> z = std::vector<std::vector<double>>(data->getNItems(), std::vector<double>(data->getNItems(), 0));


    master.solveMasterProblem();

	if (node->is_root){

		while (true){
			if (master.rmp.getCplexStatus() == IloCplex::Infeasible){
				break;
			}

			IloNumArray pi(env2, data->n_items);
        	master.rmp.getDuals(pi, master.partition_constraint);

			// pricing.setObjectiveFunction(pi);

			
			// pricing.solvePricingProblem();

			// std::cout << pi << " ";
			
			// std::cout << std::endl;
			
	
			// IloNumArray entering_col(env2, data->n_items);
			// pricing.pricing_problem.getValues(entering_col, pricing.x);

			// //  for (int i = 0; i < data->n_items; i++)
			// // {
			// // 	std::cout << entering_col[i] << " ";
			// // }

			// // std::cin.get();

			// double pricing_obj = pricing.pricing_problem.getObjValue();
			
			

			IloNumArray entering_col(env2, data->n_items);
			//std::vector<double> column;
			double pricing_obj = pricing.solveCombo(pi, entering_col);

			// for (int i = 0; i < data->n_items; i++)
			// {
			// 	entering_col[i] = column[i];
			// 	//std::cout << entering_col[i] << " ";
			// }
			// std::cin.get();

			


			if (1 - pricing_obj < -EPS){

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
			}

			else
			{
				//std::cout << "No column with negative reduced costs found. The current basis is optimal" << std::endl;

				break;
			}
		
		}
	}
	else{
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

				for (int i = 0; i < entering_col.getSize(); i++)
				{
				    if (entering_col[i] > 0.9)
				    {
				        entering_col[i] = 1;
				    }
				    else
				    {
				        entering_col[i] = 0;
				    }
				}
				

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

				pricing.pricing_problem.clear();
				pricing.pricing_problem.end();
				break;
			}

			//pricing.pricing_model.end();
		}

	}

	
	if (master.rmp.getCplexStatus() == IloCplex::Infeasible){
		prune(master.lambda);
		return std::make_pair(-1, -1);
	}


	node->LB = master.rmp.getObjValue();

	//std::cout << "LB: " << node->LB << std::endl;

	
	IloNumArray lambda_values(env, master.lambda.getSize());
	master.rmp.getValues(lambda_values, master.lambda);

	// Prune if there's a artificial variable with a positive value
	// for (int i = 0; i < data->n_items; i++){
	// 	if (lambda_values[i] > EPS){
	// 		prune(master.lambda);

	// 		//pricing.pricing_model.end();
	// 		//env2.end();
	// 		//lambda_values.end();
	// 		master.rmp.clear();
	// 		master.rmp.end();
	// 		//std::cout << "artificial variable" << std::endl;

	// 		return  std::make_pair(-1, -1);
	// 	}
    // }

	// Prune if the LB is greater than the best integer solution
	if (std::ceil(master.rmp.getObjValue() - EPS) - best_integer >= 0){
		prune(master.lambda);
		//std::cout << "LB maior que UB" << std::endl;

		//pricing.pricing_model.end();
		//env2.end();
		//lambda_values.end();
		master.rmp.clear();
		master.rmp.end();
		return std::make_pair(-1, -1);
	}

	double most_fractional = std::numeric_limits<double>::infinity();
	std::pair<int, int> fractional_pair = getFractionalPair(most_fractional, z, lambda_values, master.lambda_items);

	 //= std::abs(0.5 - z[fractional_pair.first][fractional_pair.second]);

	


	if (std::abs(0.5 - most_fractional) < EPS){
		//std::cout << "Integer solution found: " << master.rmp.getObjValue() << std::endl;
		if (master.rmp.getObjValue() < best_integer){
			best_integer = master.rmp.getObjValue();
		}

		prune(master.lambda);

		//pricing.pricing_model.end();
		//env2.end();
		//lambda_values.end();
		master.rmp.clear();
		master.rmp.end();

		return std::make_pair(-1, -1);
	}


	lambda_values.end();

	for (int i = 0; i < master.lambda.getSize(); i++)
    {
        master.lambda[i].setUB(IloInfinity);
    }

	master.rmp.clear();
    master.rmp.end();


	std::cout << "Most fractional pair: " << fractional_pair.first << " " << fractional_pair.second << std::endl;
	return fractional_pair;


}

void BP::BranchAndPrice(bool use_combo){

	Node root;
	std::vector<Node> tree;

	tree.push_back(root);
	tree[0].LB = 0;
	tree[0].UB = data->n_items;

	if (use_combo)
		tree[0].is_root = true;
	//tree[0].separated = std::vector<std::pair<int, int>>(data->n_items, std::make_pair(-1, -1));
	//tree[0].merged = std::vector<std::pair<int, int>>(data->n_items, std::make_pair(-1, -1));

	while(!tree.empty()){

		Node current = tree.back();
		std::pair<int, int> p = columnGeneration(&current);
		//std::cout << p.first <<" " <<p.second << std::endl;	

		//std::cin.get();

		if (p.first == -1 && p.second == -1){
			tree.pop_back();
			continue;
		}

		Node left;
		Node right;

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

		tree.pop_back();
	}

	

	//std::cout << "Best solution found: " << best_integer << std::endl;

}

