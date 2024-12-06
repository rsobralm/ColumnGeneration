#include "Pricing.h"
#include "combo.c"

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
    // IloExpr reduced_cost(env);
    // reduced_cost += 1;
    // for (int i = 0; i < data->n_items; i++)
    // {
    //     reduced_cost -= pi[i] * x[i];
    // }
    // pricing_model.add(IloMinimize(env, reduced_cost));

    //RESTRIÇÃO: SOMA DOS PESOS DOS ITENS SELECIONADOS NÃO PODE SER MAIOR QUE A CAPACIDADE
    IloExpr sum_weights(env);
    for (int i = 0; i < data->n_items; i++)
    {
        sum_weights += data->weights[i] * x[i];
    }
    pricing_model.add(sum_weights <= data->bin_capacity);

    objective_function = IloMinimize(env);


    pricing_model.add(objective_function);



}

void Pricing::setObjectiveFunction(IloNumArray &pi){
    IloExpr reduced_cost(env);
    reduced_cost += 1;
    for (int i = 0; i < data->n_items; i++)
    {
        reduced_cost -= pi[i] * x[i];
    }

    objective_function.setExpr(reduced_cost);
    //pricing_model.add(IloMinimize(env, reduced_cost));
}

void Pricing::solvePricingProblem(){
    

    pricing_problem = IloCplex(pricing_model);
    pricing_problem.setParam(IloCplex::Param::Threads, 1);
    pricing_problem.setOut(env.getNullStream()); // disables CPLEX log
    pricing_problem.solve();

}

double Pricing::solveCombo(IloNumArray &pi, IloNumArray &entering_col){
    int n_items = data->getNItems();
    item items[n_items];

    for (int i = 0; i < n_items; i++)
    {
        items[i].p = pi[i] * BIG_M;
        items[i].w = data->weights[i];
        items[i].x = 0;
        items[i].id = i;
        //std::cout << "Item " << i << " - p: " << pi[i] << " w: " << items[i].w << std::endl;
    }

    long capacity = data->getBinCapacity();

    //std::vector<double> entering_col;

    // Define bounds
    long lower_bound = 0; 
    long upper_bound = BIG_M; 

    bool define_solution = true; // Update solution vector
    bool relaxed = true;        // No relaxed problem

    item* last_item = &items[n_items - 1];
    item* first_item = &items[0];

    // Call the combo function to solve the knapsack problem
    long cost = combo(first_item, last_item, capacity, lower_bound, upper_bound, define_solution, relaxed);

    

    // std::cout << "Cost: " << cost << std::endl;
    cost *= EPS;

    //entering_col = std::vector<double>(n_items, 0.0);

    for (int i = 0; i < n_items; i++)
    {   
        int item_id = items[i].id;
        //std::cout << "Item " << item_id << " - x: " << items[i].x << " " << entering_col[i] <<std::endl;
        entering_col[item_id] = (items[i].x);
    }

    return cost;

}



