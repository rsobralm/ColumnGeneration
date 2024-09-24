#include <ilcplex/ilocplex.h>
#include <vector>
#include <cstdlib>
#include "Data.h"

using namespace std;

int main(int argc, char **argv) 
{
	const double M = 1e6;
	// vector<int> weight = {2, 1, 3, 3, 5};
	// int capacity = 7;
	// int n = weight.size();

	Data data;
	data.readData(argv[1]);

	cout << data.getBinCapacity() << endl;
	//exit(1);

	int capacity = data.getBinCapacity();
	int n = data.getNItems();
	vector<int> weight;
	
	for (unsigned  i = 0; i < n; i++)
	{
		weight.push_back(data.getItemWeight(i));
	}

	IloEnv env;
	IloModel master_model(env);

	IloNumVarArray lambda(env, n, 0, IloInfinity);

	IloExpr sum_obj(env);
	IloRangeArray partition_constraint(env);

	for (int i = 0; i < n; i++)
	{
		char var_name[50];
		sprintf(var_name, "y%d", i);

		lambda[i].setName(var_name);
		sum_obj += M * lambda[i];

		partition_constraint.add(lambda[i] == 1);
	}

	master_model.add(partition_constraint);

	IloObjective master_objective = IloMinimize(env, sum_obj);
	master_model.add(master_objective);

	IloCplex rmp(master_model);

	rmp.setOut(env.getNullStream()); // disables CPLEX log

	rmp.solve();

	cout << "Initial lower bound: " << rmp.getObjValue() << endl;

	cout << "Initial solution: ";
	for (size_t j = 0; j < lambda.getSize(); j++)
	{
		cout << rmp.getValue(lambda[j]) << " ";
	}
	cout << endl;

	int lambda_counter = n;
	while(true)
	{
		// Get the dual variables
		IloNumArray pi(env, n);

		rmp.getDuals(pi, partition_constraint);

		// for (size_t i = 0; i < n; i++)
		// {
		// 	cout << "Dual variable of constraint " << i << " = " << pi[i] << endl;
		// }

		// Build and solve the pricing problem
		// ...
        IloModel pricing_model(env);
        IloBoolVarArray x(env, n);

        //FO: MINIMIZAR  O CUSTO REDUZIDO
        IloExpr reduced_cost(env);
		reduced_cost += 1;
        for (int i = 0; i < n; i++)
        {
            reduced_cost -= pi[i] * x[i];
        }
        pricing_model.add(IloMinimize(env, reduced_cost));

        //RESTRIÇÃO: SOMA DOS PESOS DOS ITENS SELECIONADOS NÃO PODE SER MAIOR QUE A CAPACIDADE
        IloExpr sum_weights(env);
        for (int i = 0; i < n; i++)
        {
            sum_weights += weight[i] * x[i];
        }
        pricing_model.add(sum_weights <= capacity);

        
        IloCplex pricing_problem(pricing_model);

		pricing_problem.setOut(env.getNullStream()); // disables CPLEX log
		pricing_problem.solve();
        

		if (pricing_problem.getObjValue() < -1e-5)
		{

			cout << "Reduced cost is equal to " << pricing_problem.getObjValue() << ", which is less than 0..." << endl;

			IloNumArray entering_col(env, n);

			pricing_problem.getValues(x, entering_col);

			cout << endl << "Entering column:" << endl;
			// for (size_t i = 0; i < n; i++)
			// {
			// 	cout << (entering_col[i] < 0.5 ? 0 : 1) << endl;
			// }
			// cout << endl;

			// Add the column to the master problem
			// (the cost of the new variable is always 1)
			char var_name[50];
			sprintf(var_name, "y%d", lambda_counter++);
			IloNumVar new_lambda(master_objective(1) + partition_constraint(entering_col), 0, IloInfinity);
			new_lambda.setName(var_name);

			lambda.add(new_lambda);

			cout << "Solving the RMP again..." << endl;

			// ...
			rmp.solve();
		}
		else
		{
			cout << "No column with negative reduced costs found. The current basis is optimal" << endl;
			// cout << "Final master problem: " << endl;
			// system("cat model.lp");
			break;
		}
	}

	cout << endl;
	// cout << "Forcing items 1 and 2 to be separated in the master (for branch-and-price only): " << endl;
	// // 0 1 2 3 4 5 6 7 8 9 10 11
	// //                           
	// // 1 0 0 0 0 1 1 1 0 1  0  0
	// // 0 1 0 0 0 1 1 0 0 0  1  1
	// // 0 0 1 0 0 1 0 1 1 0  0  1
	// // 0 0 0 1 0 0 1 0 1 0  0  1
	// // 0 0 0 0 1 0 0 0 0 1  1  0

	// // itens 1 and 2 are together only on columns 5 and 11
	// lambda[11].setUB(0.0);
	// lambda[5].setUB(0.0);

	// // to allow them again:
	// // lambda[5].setUB(IloInfinity);
	// // lambda[11].setUB(IloInfinity);

	rmp.solve();
	int n_lambdas = 0;

	for (size_t j = 0; j < lambda.getSize(); j++)
	{	
		if (rmp.getValue(lambda[j]) > 0.5)
			n_lambdas++;
			//cout << <<rmp.getValue(lambda[j]) << " ";
	}
	cout << endl;
	env.end();

	cout << n_lambdas << endl;

	return 0;
}