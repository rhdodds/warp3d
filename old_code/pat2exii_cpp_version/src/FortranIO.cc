#include "FortranIO.h"



char* readAndCheck(ifstream& file, int bytes) {
      // Check to make sure the file is open.
      if (!file.is_open()) {
            cerr << "Error: File not open." << endl;
            exit(1);
      }

      // Allocate buffer
      char* buffer = new char[bytes+2*sizeof(int)];

      // Get the full buffer
      file.read(buffer,bytes+2*sizeof(int));

      // Pull out our headers
      int* header1p = (int*) buffer; // Address of buffer[0]
      int* header2p = (int*) &(buffer[bytes+sizeof(int)]);
      int header1 = *header1p;
      int header2 = *header2p;

      // Check data
      if (header1!=bytes || header2!=bytes) {
            cerr << "Error: Packet headers do not match requested data" << endl;
            delete[] buffer;
            exit(1);
      }

      // Crop off the first and last header data
      buffer = (char*) memmove(buffer, &(buffer[0+sizeof(int)]), bytes);
      buffer = (char*) realloc(buffer,bytes);

      // Return
      return buffer;
}

// Try to find a packet size
int findPacket(ifstream& file) {
      // Check to make sure the file is open.
      if (!file.is_open()) {
            cerr << "Error: File not open." << endl;
            exit(1);
      }
      char* buffer = (char*) malloc(1); 
      int s;
      for (s=1; ;s++) {
            // Allocate buffer
            buffer = (char*) realloc(buffer,(s+2*sizeof(int))*sizeof(char));

            // Get the full buffer
            file.read(buffer,s+2*sizeof(int));

            // Pull out our headers
            int* header1p = (int*) buffer; // Address of buffer[0]
            int* header2p = (int*) &(buffer[s+sizeof(int)]);
            int header1 = *header1p;
            int header2 = *header2p;

            // Rewind
            // Get the actual count of operations performed
            int count = file.gcount();
            // Position of point
            int pos = file.tellg();
            // Rewind
            file.seekg(pos-count);
            
            // See if we found a valid pointer
            if (header1==header2) break;
      }

      free(buffer);

      return s;

}
