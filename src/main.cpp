#include <iostream>
#include <ilcplex/ilocplex.h>
#include "Data.h"
#include "BP.h"


int main(int argc, char **argv) 
{
    Data data;
    data.readData(argv[1]);

    BP bp(&data);
    bp.columnGeneration();

    return 0;
}