#ifndef NODE_H
#define NODE_H

#include <iostream>
#include <set>

class Node
{
    public:
        std::set<std::pair<int, int>> separated, merged;
        int UB;
        int LB;

        Node(){
            
        };

        ~Node(){
            
        };


};

#endif
