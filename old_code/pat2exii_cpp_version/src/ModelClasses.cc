#include "ModelClasses.h"

// Helper to set up the model structure
void ModelStruct::initialize(int nodesIn, int elementsIn, set<int> blocks_ns) {
      nodes = nodesIn;
      elements = elementsIn;
      blocks = blocks_ns.size();

      x.resize(nodes);
      y.resize(nodes);
      z.resize(nodes);
      
      // Make our blank blockdata
      for (set<int>::iterator it=blocks_ns.begin(); it!=blocks_ns.end(); ++it) {
            ElemBlock A;
            blockData[*it] = A;
      }
}

// Check nodal coincidence for l3disop versus inter_8 differentiation
bool checkCoincidence(ModelStruct& theModel, vector<int> connectivity) {
      cout << theModel.y.size() << endl;
      return            realComp(theModel.x[connectivity[0]],theModel.x[connectivity[4]])
            &&          realComp(theModel.x[connectivity[1]],theModel.x[connectivity[5]])
            &&          realComp(theModel.x[connectivity[2]],theModel.x[connectivity[6]])
            &&          realComp(theModel.x[connectivity[3]],theModel.x[connectivity[7]])
            &&          realComp(theModel.y[connectivity[0]],theModel.y[connectivity[4]])
            &&          realComp(theModel.y[connectivity[1]],theModel.y[connectivity[5]])
            &&          realComp(theModel.y[connectivity[2]],theModel.y[connectivity[6]])
            &&          realComp(theModel.y[connectivity[3]],theModel.y[connectivity[7]])
            &&          realComp(theModel.z[connectivity[0]],theModel.z[connectivity[4]])
            &&          realComp(theModel.z[connectivity[1]],theModel.z[connectivity[5]])
            &&          realComp(theModel.z[connectivity[2]],theModel.z[connectivity[6]])
            &&          realComp(theModel.z[connectivity[3]],theModel.z[connectivity[7]]);

}
