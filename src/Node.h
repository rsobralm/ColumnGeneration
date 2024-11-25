#ifndef NODE_H
#define NODE_H

#include <iostream>
#include "Data.h"

class Node
{
    public:
        Data *data;
        std::vector<std::pair<int, int>> separated, merged;
        int UB;
        int LB;
};

#endif
