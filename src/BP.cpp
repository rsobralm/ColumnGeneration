#include "BP.h"


BP::BP(Data *data) : data(data)
{
    this->capacity = data->getBinCapacity();
    this->n = data->getNItems();
}

BP::~BP()
{
}

void BP::columnGeneration(){
    
}
