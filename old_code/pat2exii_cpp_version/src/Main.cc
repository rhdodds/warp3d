#include "Main.h"

int main(int argc, char** argv) {
      // Ask for the new filename
      cout << "Input filename for new Exodus II file (without extension)." << endl;
      cout << "Filename: ";
      std::string exfilename;
      cin >> exfilename;
      cout << endl;

      // Append ".exo"
      exfilename += ".exo";

      // Ask for the location of the patran neutral
      cout << "Input path to Patran neutral file (with extension)." << endl;
      cout << "Path: ";
      std::string patfilename;
      cin >> patfilename;
      cout << endl;

      // Ask for the location of the results files
      cout << "Input path to DIRECTORY containing Patran results files." << endl;
      cout << "Path: ";
      std::string resultspath;
      cin >> resultspath;
      cout << endl;

      // Create object
      Translator theTranslator(exfilename,patfilename,resultspath);
      // Read/write patran neutral
      theTranslator.readNeutralFile();

      // Read/write results
      theTranslator.processResultsDirectory();

      cout << "Done..." << endl;

      return 0;
}
