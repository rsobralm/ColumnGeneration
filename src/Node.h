#ifndef NODE_H
#define NODE_H

#include <iostream>
#include "Data.h"
#include <set>

class Node
{
    public:
        Data *data;
        std::set<std::pair<int, int>> separated, merged;
        int UB;
        int LB;
        Node(Data *data){
            this->data = data;
        };

        Node(){
            
        };

        ~Node(){
            
        };


};

#endif
