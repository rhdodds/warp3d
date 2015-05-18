#ifndef MODELCLASSES_H
#define MODELCLASSES_H

#include <vector>
#include <map>
#include <set>

#include "Helpers.h"

using namespace std;

// Short struct to hold the essential data for an element block
class ElemBlock {
      public:
            vector<int>             elNumbers;
            vector< vector<int> >   connectivity;
            ElementType             eType;
};

// Short struct to hold all the model data before it's written out
class ModelStruct {
      public:
            int                           nodes,elements,blocks;
            vector<double>                x,y,z;
            map<int,ElemBlock>            blockData;

            void initialize(int nodesIn, int elementsIn, set<int> blocks);
};

// Short struct to represent a results header (but not data)
class ResultHeader {
      public:
            string            name;
            FieldLocation     location;
            int               index;

};

// Helper to check nodal coincidence
bool checkCoincidence(ModelStruct& theModel, vector<int> connectivity);


#endif
