#include <iostream>
#include <ilcplex/ilocplex.h>
#include "Data.h"
#include "BP.h"
#include <chrono>

int main(int argc, char **argv) 
{
    Data data;
    data.readData(argv[1]);
    bool use_combo = false;

    if (argc > 2 and std::string(argv[2]) == "-p")
    {
        use_combo = true;
    }


    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    BP bp(&data);
    bp.BranchAndPrice(use_combo);

    end = std::chrono::system_clock::now();

    int elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << argv[1] << ";" << bp.best_integer << ";" << elapsed_seconds / 1000.0 << "\n";

    return 0;
}