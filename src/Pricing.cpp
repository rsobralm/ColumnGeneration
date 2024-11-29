#include "Pricing.h"

Pricing::Pricing(Data *data, IloEnv env, IloNumArray pi) : data(data), env(env)
{
    this->pi = pi;
}

Pricing::~Pricing()
{
}

void Pricing::buildPricingProblem()
{
    //IloModel pricing_model(env);
    pricing_model = IloModel(env);
    //IloBoolVarArray x(env, data->n_items);
    x = IloBoolVarArray(env, data->n_items);

    //FO: MINIMIZAR  O CUSTO REDUZIDO
    IloExpr reduced_cost(env);
    reduced_cost += 1;
    for (int i = 0; i < data->n_items; i++)
    {
        reduced_cost -= pi[i] * x[i];
    }
    pricing_model.add(IloMinimize(env, reduced_cost));

    //RESTRIÇÃO: SOMA DOS PESOS DOS ITENS SELECIONADOS NÃO PODE SER MAIOR QUE A CAPACIDADE
    IloExpr sum_weights(env);
    for (int i = 0; i < data->n_items; i++)
    {
        sum_weights += data->weights[i] * x[i];
    }
    pricing_model.add(sum_weights <= data->bin_capacity);



}

void Pricing::solvePricingProblem(){
    pricing_problem = IloCplex(pricing_model);
    pricing_problem.setOut(env.getNullStream()); // disables CPLEX log
    pricing_problem.solve();

}

void Pricing::addBranchingConstraints(std::set<std::pair<int, int>> &together, std::set<std::pair<int, int>> &separated){
    for (auto &p : together){

        pricing_model.add(x[p.first] == x[p.second]);
    }

    for (auto &p : separated){

        pricing_model.add(x[p.first] + x[p.second] <= 1);

    }
}


