#ifndef TRANSLATOR_H
#define TRANSLATOR_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <set>
#include <map>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

#include "FortranIO.h"
#include "Helpers.h"
#include "exodusII.h"
#include "ModelClasses.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

using namespace std;


// Main class
class Translator {
      public:
            Translator(std::string exfile, std::string patfile, std::string resultsdir);
            void readNeutralFile();
            void processResultsDirectory();

      private:
            void processNeutralFile(bool lookAhead);
            void processNeutralPacket(char* buffer, int length, bool lookAhead, ModelStruct& theModel);
            void readBinaryResults(int step, std::string filename, std::string id);
            void readASCIIResults(int step, std::string filename, std::string id);
            void parseTemplate(std::string file, std::string ID);
            void defaultTemplate(std::string ID);

      private:
            // The three user inputs
            ifstream          patranFile_;
            int               exFileID_;
            std::string       resultsDir_;
            // Also store some common data
            std::string       modelName_;
            int               modelNodes_;
            int               modelElements_;
            int               modelDimensionality_;
            set<int>          modelBlocks_;
            vector< vector<int> > blockNodeLists_;
            // Result set data
            int               numNodeResults_;
            int               numElementResults_;
            map<pair<string,int>, ResultHeader> headerMap_;
            map<int,int>      actualToExodusTimeMap_; // Map actual time step to exodus time step

};

#endif
