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
        Node(Data *data);

    Node::Node(Data *data)
    {
        this->data = data;
        this->UB = 0;
        this->LB = 0;
    }
};

#endif
