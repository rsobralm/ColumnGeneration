#ifndef _DATA_H_
#define _DATA_H_

#include <vector>
#include <stdio.h>

class Data
{
   public:

      int bin_capacity;
      int n_items;
      std::vector<int> weights;

      void readData(char* filePath);

      int getNItems();

      int getBinCapacity();

      int getItemWeight(unsigned int item);
};

#endif