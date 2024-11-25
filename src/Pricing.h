#ifndef PRICING_H
#define PRICING_H

#include <ilcplex/ilocplex.h>

class Pricing
{

public:
    IloEnv env;
    IloNumArray pi; // variáveis duais
    Data *data;

    Pricing(Data *data, IloEnv env, IloNumArray pi);
    void buildPricingProblem();
    void solvePricingProblem();
    ~Pricing();
};

Pricing::Pricing(Data *data, IloEnv env, IloNumArray pi) : data(data), env(env)
{
    this->pi = pi;
}

Pricing::~Pricing()
{
}

void Pricing::buildPricingProblem()
{
    IloModel pricing_model(env);
    IloBoolVarArray x(env, data->n_items);

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

    
    IloCplex pricing_problem(pricing_model);

    pricing_problem.setOut(env.getNullStream()); // disables CPLEX log
    pricing_problem.solve();

}



#endif