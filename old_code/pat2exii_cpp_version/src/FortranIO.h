#ifndef FORTRANIO_H
#define FORTRANIO_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>

using namespace std;


char* readAndCheck(ifstream& file, int bytes);
int findPacket(ifstream& file);

#endif
