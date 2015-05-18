#include "Helpers.h"

std::string readableElementType(ElementType type) {
      std::string name;
      switch (type) {
            case l3disop: 
                  name = "l3disop";
                  break;
            case ts9isop:
                  name = "ts9isop";
                  break;
            case ts12isop:
                  name = "ts12isop";
                  break;
            case ts15isop:
                  name = "ts15isop";
                  break;
            case q3disop:
                  name = "q3disop";
                  break;
            case tet4:
                  name = "tet4";
                  break;
            case tet10:
                  name = "tet10";
                  break;
            case inter_8:
                  name = "inter_8";
                  break;
            case trint6:
                  name = "trint6";
                  break;
            case trint12:
                  name = "trint12";
                  break;
      }
      return name;
}

std::string exodusElementType(ElementType type) {
      std::string name;
      switch (type) {
            case l3disop: 
                  name = "HEX";
                  break;
            case ts9isop:
                  name = "HEX";
                  break;
            case ts12isop:
                  name = "HEX";
                  break;
            case ts15isop:
                  name = "HEX";
                  break;
            case q3disop:
                  name = "HEX";
                  break;
            case tet4:
                  name = "TETRA";
                  break;
            case tet10:
                  name = "TETRA";
                  break;
            case inter_8:
                  name = "HEX";
                  break;
            case trint6:
                  name = "WEDGE";
                  break;
            case trint12:
                  name = "WEDGE";
                  break;
      }
      return name;
}

// Helpers which split up std::strings
std::vector<int> splitToInt(std::string in) {
      std::istringstream iss(in);
      vector<int> tokens;
      copy(istream_iterator<int>(iss),istream_iterator<int>(),
                  back_inserter<vector<int> >(tokens));
      return tokens;

}

std::vector<double> splitToDouble(std::string in) {
      std::istringstream iss(in);
      vector<double> tokens;
      copy(istream_iterator<double>(iss),istream_iterator<double>(),
                  back_inserter<vector<double> >(tokens));
      return tokens;
}

std::vector<std::string> splitToString(std::string in) {
      std::istringstream iss(in);
      vector<std::string> tokens;
      copy(istream_iterator<std::string>(iss),istream_iterator<std::string>(),
                  back_inserter<vector<std::string> >(tokens));
      return tokens;
}

std::string strip(std::string input) {
      std::string::iterator it;
      it = remove(input.begin(),input.end(),' ');
      input.resize(it - input.begin());
      return input;
}

// Check an Exodus error code
void checkExError(int code) {
      if (code < 0) {
            cerr << "Error: Exodus error." << endl;
            exit(1);

      }
}


bool realComp(double a, double b) {
      if ( fabs(a-b) < 0.0000001) return true;
      else return false;

}

void readTimeMap(int steps, vector<double>& timeMap) {
      // Ask for the file location
      cout << "Input path to time file." << endl;
      cout << "File: ";
      string path;
      cin >> path;
      cout << endl;
      
      // Open
      ifstream ifs(path.c_str(), ifstream::in);
      // check
      if (!ifs) {
            cerr << "Error: file could not be opened." << endl;
            exit(1);
      }

      // Read steps lines and insert to timeMap
      timeMap.resize(steps);
      for (int i=0; i<steps; i++) {
            char buffer[256];
            ifs.getline(buffer, 256);
            string cppstring(buffer);
            istringstream iss(cppstring);
            double time;
            iss >> time;
            timeMap[i] = time;
      }
      ifs.close();

}
