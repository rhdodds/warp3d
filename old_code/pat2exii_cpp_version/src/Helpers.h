#ifndef HELPERS_H
#define HELPERS_H

#include <fstream>
#include <string>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <sstream>
#include <iterator>
#include <math.h>

using namespace std;

// Definition of element types
enum ElementType { l3disop, ts9isop, ts12isop, ts15isop, q3disop, 
      tet4, tet10, inter_8, trint6, trint12 };

// Helper function to return a human readable element name
std::string readableElementType(ElementType type);
std::string exodusElementType(ElementType type);

// Definition of field locations
enum FieldLocation { element, nodal };

// Various helpers
std::vector<int> splitToInt(std::string in);
std::vector<double> splitToDouble(std::string in);
std::vector<std::string> splitToString(std::string in);
std::string strip(std::string input);
void checkExError(int code);

bool realComp(double a, double b);

// Helper to input a time map file
void readTimeMap(int steps, vector<double>& timeMap);

#endif
