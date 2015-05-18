#include "Translator.h"


Translator::Translator(std::string exfile, std::string patfile, std::string resultsdir) {
      // Open up the patran file
      cout << "Opening file " << patfile << "." << endl;
      patranFile_.open(patfile.c_str(),ifstream::in);
      if (!patranFile_) {
            cerr << "Error: file could not be opened." << endl;
            exit(1);
      }

      // Store the directory of results
      resultsDir_ = resultsdir;

      // Look ahead in the neutral file to get some key data
      processNeutralFile(true);

      // Open up the new exodus file
      int errorCode;
      int comp_ws = sizeof(double); // Size we store our reals in
      int io_ws = sizeof(double); // Size we want to store the reals in
      exFileID_ = ex_create(exfile.c_str(), EX_CLOBBER, &comp_ws, &io_ws);
      checkExError(exFileID_);

      // Just hack the dimensionality to 3
      modelDimensionality_ = 3;
      // Write the look-ahead data
      errorCode = ex_put_init(exFileID_, modelName_.c_str(), modelDimensionality_,
                        modelNodes_, modelElements_, modelBlocks_.size(),
                        1, 0);
      checkExError(errorCode);

}

// Wrapper to "processNeutralFile(false)"
void Translator::readNeutralFile() {
      processNeutralFile(false);
      // Now we're done with the file
      patranFile_.close();
}

// Process an entire neutral file either in look ahead or full mode
void Translator::processNeutralFile(bool lookAhead) {
      // Setup the model and initialize it if this is an actual write pass
      ModelStruct theModel;
      if (!lookAhead) theModel.initialize(modelNodes_,modelElements_,modelBlocks_);
      
      // Loop on input, reading until EOF
      char linebuffer[256];
      while (!patranFile_.eof()) {
            patranFile_.getline(linebuffer,256);
            processNeutralPacket(linebuffer,256,lookAhead, theModel);
      }
      // Rewind the file - issue
      patranFile_.clear();
      patranFile_.seekg(0, ios::beg);

      // If this isn't just the look ahead run, then actually write
      // to the database
      if (!lookAhead) {
            int errorCode;
            // Write the coordinates
            errorCode = ex_put_coord(exFileID_, (void*) &(theModel.x[0]),
                              (void*) &(theModel.y[0]), (void*) &(theModel.z[0]));
            checkExError(errorCode);
            
            // Write the coordinate system ids
            char* coord_names[3];
            coord_names[0] = "X";
            coord_names[1] = "Y";
            coord_names[2] = "Z";
            errorCode = ex_put_coord_names(exFileID_, coord_names);
            checkExError(errorCode);

            // Don't need to map node numbers, defaults to 1...nodes

            // For user convenience, write a node set containing all the nodes
            // Give it ID 1.
            errorCode = ex_put_node_set_param(exFileID_, 1, theModel.nodes, 0);
            checkExError(errorCode);
            // Make a vector of all the nodes
            vector<int> nodeList(theModel.nodes,0);
            for (int i=1; i<=theModel.nodes; i++) {
                  nodeList[i-1]=i;
            }
            errorCode = ex_put_node_set(exFileID_, 1, &(nodeList[0]));
            checkExError(errorCode);
            // Name it All Nodes
            std::string node_name = "All Nodes";
            errorCode = ex_put_name(exFileID_,EX_NODE_SET,1,node_name.c_str());
            checkExError(errorCode);

            // On to elements
            // Make a vector of length elements to hold the global id map
            vector<int> elMap;
            int blockID = 1;
            // Iterate over each type of element in the model
            for (map<int,ElemBlock>::iterator it = theModel.blockData.begin();
                        it != theModel.blockData.end(); ++it) {
                  // Insert the block
                  errorCode = ex_put_elem_block(exFileID_, blockID, 
                              exodusElementType((it->second).eType).c_str(),
                              (it->second).elNumbers.size(),
                              (it->second).connectivity[0].size(),
                              0);
                  checkExError(errorCode);

                  // Having issues with vector of vectors...
                  // Need to copy I guess
                  int* connect = (int*) malloc((it->second).elNumbers.size()*
                               (it->second).connectivity[0].size()*sizeof(int));

                  for (int k=0; k<(it->second).elNumbers.size(); k++) {
                        for (int j=0; j<(it->second).connectivity[0].size(); j++) {
                              connect[k*(it->second).connectivity[0].size()+j]=(it->second).connectivity[k][j];
                        }
                  }


                  // Insert the connectivities
                  errorCode = ex_put_elem_conn(exFileID_, blockID, 
                              connect);
                  checkExError(errorCode);

                  free(connect);

                  // Append the element map to the main map
                  elMap.insert(elMap.end(), (it->second).elNumbers.begin(),
                              (it->second).elNumbers.end());

                  // Store the element numbers for later
                  blockNodeLists_.push_back((it->second).elNumbers);

                  // Give the set a name
                  errorCode = ex_put_name(exFileID_, EX_ELEM_BLOCK, blockID,
                              SSTR(it->first).c_str());
                  checkExError(errorCode);

                  // increment id
                  blockID++;

            }

            errorCode = ex_put_elem_num_map(exFileID_,&elMap[0]);
            checkExError(errorCode);

      }

}

// Read an entire neutral file packet and do the required operation
// on the model.  This function must also advance the state of the
// input stream to the next packet header.
void Translator::processNeutralPacket(char* buffer, int length, bool lookAhead, ModelStruct& theModel) {
      int errorCode;
      // Convert to c++ std::string for ease of use
      std::string header_line(buffer);
      // Split to tokens
      vector<int> tokens = splitToInt(header_line);
      // Grab the data cards into a vector of strings
      vector<std::string> data;
      for (int i=1; i <= tokens[3]; i++) {
            patranFile_.getline(buffer,length);
            data.push_back(std::string(buffer));
      }
      // Switch on the type of packet and call appropriate function to process
      switch(tokens[0]) {
            case 25: // Title
            {
                  if (lookAhead) {
                        modelName_=data[0];
                  }
                  else {
                        // This was done in look ahead
                  }
                  break;
            }

            case 26: // Summary
            {
                  if (lookAhead) {
                        modelNodes_ = tokens[4];
                        modelElements_ = tokens[5];
                  }
                  else { 
                        // Again, done in look ahead
                  }
                  break;
            }

            case 1: // Node Data
            {
                  if (lookAhead) {
                        // Do nothing
                  }
                  else {
                        // Split the data packet to doubles
                        vector<double> coordinates = splitToDouble(data[0]);
                        // Insert to model
                        // tokens[1] is the node number, we have the coordinates
                        // in coordinates.
                        theModel.x[tokens[1]-1]=coordinates[0];
                        theModel.y[tokens[1]-1]=coordinates[1];
                        theModel.z[tokens[1]-1]=coordinates[2];
                  }
                  break;
            }

            case 2: // Element Data
            {
                  int               enumb;
                  vector<std::string>    splits;
                  vector<int>       spliti;
                  ElementType       enumType;
                  int               elnodes;
                  int               etype;
                  int               cards;
                  int               block;
                  vector<int>       connectivity;
                  // The element number is the 2nd position in the header
                  enumb=tokens[1]; 
                  // The first data card, in the first position, has the
                  // number of nodes for the element
                  splits = splitToString(data[0]);
                  istringstream(splits[0]) >> elnodes;
                  istringstream(splits[1]) >> block;
                  // The shape lives in the third position in the header
                  etype=tokens[2];
                  // Determine the number of cards
                  cards = (elnodes + 9)/10;
                  // Append the connectivity onto the connectivity vector
                  for (int i=1; i<=cards; i++) {
                        spliti = splitToInt(data[i]);
                        connectivity.insert(connectivity.end(),spliti.begin(),spliti.end());
                  }
                  // Now we need to figure out what type of element each one is
                  // and insert it.
                  // While we're at it, rearrange the connectivity for EXODUS
                  if (etype==5 && elnodes==4) {
                        enumType = tet4;
                        // Connectivity ok
                  }
                  else if (etype==5 && elnodes==10) {
                        enumType = tet10;
                        // Connectivity ok
                  }
                  // Need to fix for inter_8
                  else if (etype==8 && elnodes==8) {
                        enumType = l3disop;
                        // Connectivity ok
                  }
                  else if (etype==8 && elnodes==9) {
                        enumType = ts9isop;
                        // Connectivity ok
                  }
                  else if (etype==8 && elnodes==12) {
                        enumType = ts12isop;
                        // Connectivity ok
                  }
                  else if (etype==8 && elnodes==15) {
                        enumType = ts15isop;
                  }
                  else if (etype==8 && elnodes==20) {
                        enumType = q3disop;
                  }
                  else if (etype==7 && elnodes==6) {
                        enumType = trint6;
                  }
                  else if (etype==7 && elnodes==15) {
                        enumType = trint12;
                        // There is a hack here -- exodus can't read a wedge
                        // with the midside nodes on the rectangular part
                        // left out.  Instead just make the nodes the 
                        // same on the rectangular faces.
                        connectivity[9] = connectivity[0];
                        connectivity[10] = connectivity[1];
                        connectivity[11] = connectivity[2];
                  }
                  else {
                        cerr << "Error: unknown element type." << endl;
                        exit(1);
                  }
                  if (lookAhead) {
                        // Just insert the block number into our set
                        modelBlocks_.insert(block);
                  }
                  else {
                        //model_.insertelement(enumb,enumType,connectivity);
                        // Insert the element number
                        theModel.blockData[block].elNumbers.push_back(enumb);
                        // Insert the connectivity
                        theModel.blockData[block].connectivity.push_back(connectivity);
                        // Insert type
                        theModel.blockData[block].eType = enumType;
                  }

                  break;
            }
            case 3: // Material properties
            {
                  // Skip
                  break;
            }

            case 4: // Element properties
            {
                  // Skip
                  break;
            }

            case 5: // Coordinate frames
            {
                  // Skip
                  break;
            }

            case 6: // Distributed load constraint
            {
                  // Skip constraints for now
                  break;
            }

            case 7: // Nodal force constraint
            {
                  // Above
                  break;
            }

            case 8: // Nodal displacement constraint
            {
                  // Above
                  break;
            }

            case 10: // Nodal temperatures
            {
                  // Above
                  break;
            }

            case 11: // Element temperatures
            {
                  // Above
                  break;
            }

            case 99: // All done
            {
                  // That's all...
                  break;
            }

            default:
            {
                  cerr << "Warning: unknown packet type " << tokens[0] << " will not be written to Exodus file." << endl;
            }
      }

}

void Translator::processResultsDirectory() {
      int errorCode;

      // Open the directory and read files
      DIR* dp;
      struct dirent *dirp;
      vector<string> files;
      if ((dp = opendir(resultsDir_.c_str())) == NULL) {
            cerr << "Error opening directory." << endl;
            exit(1);
      }
      // Append filenames to vector
      while ((dirp = readdir(dp)) != NULL) {
            files.push_back(string(dirp->d_name));
      }
      closedir(dp);
      // Find all the result files by a ludicrous if statement (used to be a regex)
      for (vector<string>::iterator it = files.begin(); it!=files.end(); ) {
            if (((*it)[0]=='w') && ((*it)[1]=='e' || (*it)[1]=='n') && ((*it)[2]=='b'
                              || (*it)[2]=='f') && ((*it)[3]=='d' || (*it)[3]=='v'
                                    || (*it)[3]=='r' || (*it)[3]=='e' || (*it)[3]=='s' || (*it)[3]=='f'
                                    || (*it)[3]=='m' || (*it)[3]=='t' )) {
                  ++it;
            }
            else {
                  it = files.erase(it);
            }
      }

      // Exit if no files
      if (files.size()==0) {
            cerr << "Warning: No result files found." << endl;
            return;
      }
      
      // 1) Get a list of unique results types
      vector<string> uniques;
      for (vector<string>::iterator it = files.begin(); it!=files.end(); ++it) {
            uniques.push_back(it->substr(1,1)+it->substr(3,1));
      }
      sort(uniques.begin(),uniques.end());
      vector<string>::iterator it;
      it = unique(uniques.begin(),uniques.end());
      uniques.resize(it - uniques .begin() );

      // 2) Get a list of fully qualified file names
      if (resultsDir_[resultsDir_.length()-1]!='/') resultsDir_+="/";
      for (vector<string>::iterator it = files.begin(); it!=files.end(); ++it) {
            *it = resultsDir_ + *it;
      }

      // 2.5) See if the user wants to use default templates or go through and
      //      prompt for each template.
       
      cout << "Would you like to use the default templates provided in WARP3D v17?" << endl;
      cout << "(If you say no you'll be prompted for the path to each template file)." << endl;
      bool defaultTemplates = true;
      while (true) {
            cout << "y or n? ";
            char ans;
            cin >> ans;
            cout << endl;
            if (ans=='y') {
                  defaultTemplates = true;
                  break;
            }
            else if (ans=='n') {
                  defaultTemplates = false;
                  break;
            }
      }

      // 3) Get a list of the template corresponding to each file and parse them.
      map<string,string> template_paths;
      string path;
      numNodeResults_=0;
      numElementResults_=0;
      for (vector<string>::iterator it = uniques.begin(); it!=uniques.end(); ++it) {
            if (defaultTemplates) {
                  defaultTemplate(*it);
            }
            else {
                  cout << "Input path to template file for results files with prefix " << "w"+it->substr(0,1)+"*"+it->substr(1,1) << endl;
                  cout << "Path: ";
                  cin >> path;
                  cout << endl;
                  template_paths.insert(pair<string,string>(*it,path));
                  parseTemplate(path,*it);
            }
      }

      // 3.5) Write the total number of result types
      errorCode = ex_put_var_param(exFileID_, "n", numNodeResults_);
      checkExError(errorCode);
      errorCode = ex_put_var_param(exFileID_, "e", numElementResults_);
      checkExError(errorCode);

      // 3.6) Write the variable names
      vector<char*> nodeResultNames;
      vector<char*> elementResultNames;
      for (map<pair<string,int>, ResultHeader>::iterator it = headerMap_.begin(); 
                  it != headerMap_.end(); ++it) {
            if ((*it).second.location == nodal) {
                  nodeResultNames.push_back( const_cast<char*> ((*it).second.name.c_str()));
            }
            else if ((*it).second.location == element) {
                  elementResultNames.push_back( const_cast<char*> ((*it).second.name.c_str()));
            }
      }

      if (numNodeResults_ != 0) {
            errorCode = ex_put_var_names(exFileID_, "n", numNodeResults_, 
                  &(nodeResultNames[0]));
            checkExError(errorCode);
      }
      if (numElementResults_ != 0) {
            errorCode = ex_put_var_names(exFileID_, "e", numElementResults_,
                  &(elementResultNames[0]));
            checkExError(errorCode);
      }

      // 3.75) Read through all the results files once to get the time step
      //       map.  Alt: ask for a time map file
      int step;
      set<int> realSteps;
      for (vector<string>::iterator it = files.begin(); it!=files.end(); ++it) {
           // Last 5 characters of filename are the step number
           istringstream(it->substr(it->length()-5,5)) >> step;
           realSteps.insert(step);
      }
      step = 1;
      for (set<int>::iterator it = realSteps.begin(); it!=realSteps.end(); ++it) {
            actualToExodusTimeMap_.insert(pair<int,int> (*it, step));
            step++;
      }
      
      // Ask if the user wants to input a time map file.
      cout << "Would you like to input a time map file listing the actual physical" << endl;
      cout << "time corresponding to each output time step?" << endl;
      bool userMap = false;
      vector<double> timeMap;
      while (true) {
            cout << "y or n? ";
            char ans;
            cin >> ans;
            cout << endl;
            if (ans=='y') {
                  userMap = true;
                  readTimeMap(actualToExodusTimeMap_.size(), timeMap);
                  break;
            }
            else if (ans=='n') {
                  break;
            }
      }


      // Insert the time map
      int i=0;
      for (map<int,int>::iterator it = actualToExodusTimeMap_.begin(); 
                  it != actualToExodusTimeMap_.end(); ++it) {
            double time;
            if (userMap) {
                  time = timeMap[i];
            }
            else {
                  time = (double) (*it).first;
            }
            errorCode = ex_put_time(exFileID_, (*it).second, (void*) &time);
            checkExError(errorCode);
            i++;
      }

      // 3.9) Write the variable truth tables (for elements)
      if (numElementResults_ != 0) {
            int num_dim, num_nodes, num_elem, num_elem_blk,
                  num_node_sets, num_side_sets, error, exoid;
            char title[MAX_LINE_LENGTH+1];
            errorCode = ex_get_init (exFileID_, title, &num_dim, &num_nodes,
                  &num_elem, &num_elem_blk, &num_node_sets, &num_side_sets);
            checkExError(errorCode);
            // Each element block has all the variables
            int* truth_tab;
            truth_tab = (int*) calloc((num_elem_blk*numElementResults_),sizeof(int));
            int m,j,k;
            for (m=0, k=0; m<num_elem_blk; m++)
                  for (j=0; j<numElementResults_; j++)
                        truth_tab[k++] = 1;
            errorCode = ex_put_elem_var_tab(exFileID_, num_elem_blk, numElementResults_,truth_tab);
            checkExError(errorCode);
            free(truth_tab);
      }


      // 4) Loop on each file, call the appropriate reader with the 
      //    filename and template name.
      string templatepath;
      for (vector<string>::iterator it = files.begin(); it!=files.end(); ++it) {
            // Last 5 characters of filename are the step number
            istringstream(it->substr(it->length()-5,5)) >> step;

            // 8 + 6 from end are the key to the map
            string id = it->substr(it->length()-8,1)+it->substr(it->length()-6,1);

            // Announce we're opening for debug
            cout << "Opening file " + *it + "." << endl;

            // 7 from the end is whether it's formatted or binary
            if ((*it)[it->length()-7]=='f') {
                  readASCIIResults(step,*it,id);
            }
            else if ((*it)[it->length()-7]=='b') {
                  readBinaryResults(step,*it,id);
            }
            else {
                  cerr << "Error: this is bad - I can't determine if it's a binary or formatted file" << endl;
                  exit(1);
            }
      }



      // All done with exodus file
      errorCode = ex_close(exFileID_);
      checkExError(errorCode);
}

// Read in a binary results file
void Translator::readBinaryResults(int step, std::string filename, std::string id) {
      // The file specification varies between nodal and element results
      FieldLocation loc = headerMap_[pair<string,int> (id, 1)].location;
      
      int errorCode;

      // Get the exodus step number
      int exodusStep = actualToExodusTimeMap_[step];

      ifstream inputFile;
      inputFile.open(filename.c_str(),ifstream::binary);
      char* buffer;
      if (loc == nodal) {
            // NOTE FORTRAN CHARACTERS ARE INTS!!!!!
            // Read 80 character string and then 5 ints
            buffer = readAndCheck(inputFile,85*sizeof(int));
            // The 1st int is the number of nodes (no point in extracting)
            // The 5th int is the number of columns
            int* ncolsp = (int*) &(buffer[80*sizeof(int)+4*sizeof(int)]);
            int ncols = *ncolsp;

            delete [] buffer;

            // Read 2 80 character strings
            buffer = readAndCheck(inputFile,80*sizeof(int));
            delete [] buffer;
            buffer = readAndCheck(inputFile,80*sizeof(int));
            delete [] buffer;

            // Now we're in the actual data
            // For each node read in an int+ncols*float data
            // Extract the node number, then loop over the rest of the data and insert
            
            // Need a temp array...
            vector< vector<double> > tempVars;
            tempVars.resize(ncols);
            for (vector< vector<double> >::iterator it = tempVars.begin();
                        it != tempVars.end(); ++it) {
                  it->resize(modelNodes_);
            }
            
            for (int i=1; i<=modelNodes_; i++) {
                  buffer = readAndCheck(inputFile,sizeof(int)+ncols*sizeof(float));
                  // Extract the node number
                  int* nodenp = (int*) buffer;
                  int noden = *nodenp;
                  // Do a test just in case
                  if (i!=noden) {
                        cerr << "Error: Invalid file, node numbers don't match" << endl;
                        exit(1);
                  }
                  // Now loop over each column (j+1 is column number)
                  for (int j=0; j<ncols; j++) {
                        float* fp = (float*) &(buffer[sizeof(int)+j*sizeof(float)]);
                        float f = *fp;
                        double d = (double) f;
                        // Into temp array
                        tempVars[j][i-1] = d;
                  }
                  delete [] buffer;
            }
            // Insert to actual model
            // First index + 1 is column
            // Second index is nodal index
            for (int j=0; j<ncols; j++) {
                  for (int r=0; r<modelNodes_; r++) {
                  }
            }

            for (int j=0; j<ncols; j++) {
                  ResultHeader resultSet = headerMap_[pair<string,int> (id, j+1) ];
                  errorCode = ex_put_nodal_var(exFileID_, exodusStep, resultSet.index,
                              modelNodes_, (void*) &(tempVars[j][0]));
                  checkExError(errorCode);
            }

      }
      else if (loc == element) {
            // Int + 80 character string
            buffer = readAndCheck(inputFile,324);
            int* ncolsp = (int*) &(buffer[80*sizeof(int)]);
            int ncols = *ncolsp;
            delete [] buffer;

            // 2 80 character strings
            buffer = readAndCheck(inputFile,80*sizeof(int));
            delete [] buffer;
            buffer = readAndCheck(inputFile,80*sizeof(int));
            delete [] buffer;

            // Need a temp array...
            vector< vector<double> > tempVars;
            tempVars.resize(ncols);
            for (vector< vector<double> >::iterator it = tempVars.begin();
                        it != tempVars.end(); ++it) {
                  it->resize(modelElements_);
            }

            // Actual data
            // Loop on elements
            for (int i=1; i<=modelElements_; i++) {
                  buffer = readAndCheck(inputFile,2*sizeof(int)+ncols*sizeof(float));
                  // Extract the element number
                  int* elnp = (int*) buffer;
                  int eln = *elnp;
                  // Do a test just in case
                  if (i!=eln) {
                        cerr << "Error: Invalid file, element numbers don't match" << endl;
                        exit(1);
                  }
                  // This is actually the element type, if I ever need it again
                  int* idkp = (int*) &(buffer[sizeof(int)]);
                  int idk = *idkp;

                  // Now loop over each column (j+1 is column number)
                  for (int j=0; j<ncols; j++) {
                        float* fp = (float*) &(buffer[2*sizeof(int)+j*sizeof(float)]);
                        float f = *fp;
                        double d = (double) f;
                        tempVars[j][i-1]=d;
                  }
                  delete [] buffer;
            }
            // Time to write
            for (int j=0; j<ncols; j++) {
                  ResultHeader resultSet = headerMap_[pair<string,int> (id, j+1) ];
                  // Loop over each element set
                  for (int k=0; k<modelBlocks_.size(); k++) {
                        // k+1 is the set index
                        vector<double> elSetResults(blockNodeLists_[k].size(),0);
                        for (int l=0; l<elSetResults.size(); l++) {
                              elSetResults[l] = tempVars[j][blockNodeLists_[k][l]-1];
                        }
                        errorCode = ex_put_elem_var(exFileID_, exodusStep, resultSet.index,
                              k+1, elSetResults.size(), (void*) &(elSetResults[0]));
                        checkExError(errorCode);
                  }
            }

      }
      else {
            cerr << "Error: Unknown field location type." << endl;
            exit(1);
      }

      inputFile.close();
}

// Read in a ascii results file
void Translator::readASCIIResults(int step, std::string filename, std::string id) {
      // Unfortunately the damn things differ between element and nodal results
      FieldLocation loc = headerMap_[pair<string,int> (id, 1)].location;

      int exodusStep = actualToExodusTimeMap_[step];

      int errorCode;

      ifstream inputFile;
      inputFile.open(filename.c_str(),ifstream::in);
      char linebuffer[256];
      int ss = 256;

      if (loc == nodal) {
            // Skip 1 line
            inputFile.getline(linebuffer,ss);
            // Get 2nd line and convert to number of columns
            inputFile.getline(linebuffer,ss);
            string buf(linebuffer);
            istringstream conv(buf.substr(45,6));
            int total_cols;
            conv >> total_cols;
            // Skip 2 lines
            for (int i=0; i<2; i++) {
                  inputFile.getline(linebuffer,ss);
            }

            // Need a temp array...
            vector< vector<double> > tempVars;
            tempVars.resize(total_cols);
            for (vector< vector<double> >::iterator it = tempVars.begin();
                        it != tempVars.end(); ++it) {
                  it->resize(modelNodes_);
            }

            // Read number of nodes lines
            for (int i=1; i<=modelNodes_; i++) {
                  inputFile.getline(linebuffer,ss);
                  string stringbuffer(linebuffer);
                  // Each line has i8 and number of cols R13.  
                  // However, it won't put more than 5 R13s per line.
                  int node;
                  istringstream iss(stringbuffer.substr(0,8));
                  iss >> node;
                  int nlines = (total_cols-1)/5+1;
                  int line = 1;
                  int start;

                  for (int j=0; j<total_cols; j++) {
                        if (j!=0 && (j%5)==0) {
                              inputFile.getline(linebuffer,ss);
                              stringbuffer = linebuffer;
                              line++;
                        }
                        if (line==1) start = 8;
                        else start = 0;
                        double value=atof(strip(stringbuffer.substr(start+(j%5)*13,13)).c_str());
                        // Insert!
                        tempVars[j][i-1] = value;
                  }
            }

            // Now time to actually insert
            // First index + 1 is column
            // Second index is nodal index
            for (int j=0; j<total_cols; j++) {
                  ResultHeader resultSet = headerMap_[pair<string,int> (id, j+1) ];
                  errorCode = ex_put_nodal_var(exFileID_, exodusStep, resultSet.index,
                              modelNodes_, (void*) &(tempVars[j][0]));
                  checkExError(errorCode);
            }
      }
      else if (loc == element) {
            // Skip 1 line
            inputFile.getline(linebuffer,ss);
            // Get 2nd line and convert to number of columns
            inputFile.getline(linebuffer,ss);
            string buf(linebuffer);
            istringstream con(buf);
            int total_cols;
            con >> total_cols;
            // Skip 2 lines
            inputFile.getline(linebuffer,ss);
            inputFile.getline(linebuffer,ss);

            // Need a temp array...
            vector< vector<double> > tempVars;
            tempVars.resize(total_cols);
            for (vector< vector<double> >::iterator it = tempVars.begin();
                        it != tempVars.end(); ++it) {
                  it->resize(modelElements_);
            }

            // Read number of elements (note: everything very similar now, except the index info is
            // on a separate line.
            // Read number of elements lines
            for (int i=1; i<=modelElements_; i++) {
                  inputFile.getline(linebuffer,ss);
                  string stringbuffer(linebuffer);
                  // First line has the element number and some junk.
                  // Each subsequent has up to 6 R13s
                  int element;
                  istringstream iss(stringbuffer.substr(0,8));
                  iss >> element;

                  int nlines = (total_cols-1)/6+1;
                  int line = 0;
                  int start;
                  for (int j=0; j<total_cols; j++) {
                        if ((j%6)==0) {
                              inputFile.getline(linebuffer,ss);
                              stringbuffer = linebuffer;
                              line++;
                        }
                        start = 0;
                        double value=atof(strip(stringbuffer.substr(start+(j%6)*13,13)).c_str());
                        // Insert!
                        tempVars[j][i-1]=value;
                  }
            }
            // Time to write
            for (int j=0; j<total_cols; j++) {
                  ResultHeader resultSet = headerMap_[pair<string,int> (id, j+1) ];
                  // Loop over each element set
                  for (int k=0; k<modelBlocks_.size(); k++) {
                        // k+1 is the set index
                        vector<double> elSetResults(blockNodeLists_[k].size(),0);
                        for (int l=0; l<elSetResults.size(); l++) {
                              elSetResults[l] = tempVars[j][blockNodeLists_[k][l]-1];
                        }
                        errorCode = ex_put_elem_var(exFileID_, exodusStep, resultSet.index,
                              k+1, elSetResults.size(), (void*) &(elSetResults[0]));
                        checkExError(errorCode);
                  }
            }

      }
      else {
            cerr << "Error: Unknown field location type." << endl;
            exit(1);
      }
      
      inputFile.close();
}

// Use the default warp template for an id
void Translator::defaultTemplate(std::string ID) {
      // Determine type of field from the 1st character of the ID
      FieldLocation location;
      if (ID[0]=='e' ) {
            location = element;
      }
      else if (ID[0]=='n') {
            location = nodal;
      }
      else {
            cerr << "Error: Unknown field type." << endl;
            exit(1);
      }

      int ncols;
      vector<ResultHeader> headerList;

      // Determine which type of result this is
      if (ID.compare("nd")==0) {
            //cout << "Nodal displacements" << endl;
            ncols = 3;
            headerList.resize(ncols);

            numNodeResults_++;
            headerList[0].name="Disps Transx";
            headerList[0].location = location;
            headerList[0].index = numNodeResults_;

            numNodeResults_++;
            headerList[1].name="Disps Transy";
            headerList[1].location = location;
            headerList[1].index = numNodeResults_;

            numNodeResults_++;
            headerList[2].name="Disps Transz";
            headerList[2].location = location;
            headerList[2].index = numNodeResults_;
      }
      else if (ID.compare("na")==0) {
            //cout << "Nodal accelerations" << endl;
            ncols = 3;
            headerList.resize(ncols);

            numNodeResults_++;
            headerList[0].name="Accels Transx";
            headerList[0].location = location;
            headerList[0].index = numNodeResults_;

            numNodeResults_++;
            headerList[1].name="Accels Transy";
            headerList[1].location = location;
            headerList[1].index = numNodeResults_;

            numNodeResults_++;
            headerList[2].name="Accels Transz";
            headerList[2].location = location;
            headerList[2].index = numNodeResults_;

      }
      else if (ID.compare("ee")==0) {
            //cout << "Element strains" << endl;
            ncols = 22;
            headerList.resize(ncols);

            numElementResults_++;
            headerList[0].name="Strain eps-xx";
            headerList[0].location = location;
            headerList[0].index = numElementResults_;

            numElementResults_++;
            headerList[1].name="Strain eps-yy";
            headerList[1].location = location;
            headerList[1].index = numElementResults_;

            numElementResults_++;
            headerList[2].name="Strain eps-zz";
            headerList[2].location = location;
            headerList[2].index = numElementResults_;

            numElementResults_++;
            headerList[3].name="Strain gam-xy";
            headerList[3].location = location;
            headerList[3].index = numElementResults_;

            numElementResults_++;
            headerList[4].name="Strain gam-yz";
            headerList[4].location = location;
            headerList[4].index = numElementResults_;

            numElementResults_++;
            headerList[5].name="Strain gam-xz";
            headerList[5].location = location;
            headerList[5].index = numElementResults_;

            numElementResults_++;
            headerList[6].name="Strain Effective mises strain";
            headerList[6].location = location;
            headerList[6].index = numElementResults_;


            numElementResults_++;
            headerList[7].name="Strain invariant 1";
            headerList[7].location = location;
            headerList[7].index = numElementResults_;

            numElementResults_++;
            headerList[8].name="Strain invariant 2";
            headerList[8].location = location;
            headerList[8].index = numElementResults_;

            numElementResults_++;
            headerList[9].name="Strain invariant 3";
            headerList[9].location = location;
            headerList[9].index = numElementResults_;

            numElementResults_++;
            headerList[10].name="Strain principal (Minimum)";
            headerList[10].location = location;
            headerList[10].index = numElementResults_;

            numElementResults_++;
            headerList[11].name="Strain principal (Inter)";
            headerList[11].location = location;
            headerList[11].index = numElementResults_;

            numElementResults_++;
            headerList[12].name="Strain principal (Maximum)";
            headerList[12].location = location;
            headerList[12].index = numElementResults_;

            numElementResults_++;
            headerList[13].name="Strain dir cosine l1";
            headerList[13].location = location;
            headerList[13].index = numElementResults_;

            numElementResults_++;
            headerList[14].name="Strain dir cosine m1";
            headerList[14].location = location;
            headerList[14].index = numElementResults_;

            numElementResults_++;
            headerList[15].name="Strain dir cosine n1";
            headerList[15].location = location;
            headerList[15].index = numElementResults_;

            numElementResults_++;
            headerList[16].name="Strain dir cosine l2";
            headerList[16].location = location;
            headerList[16].index = numElementResults_;

            numElementResults_++;
            headerList[17].name="Strain dir cosine m2";
            headerList[17].location = location;
            headerList[17].index = numElementResults_;

            numElementResults_++;
            headerList[18].name="Strain dir cosine n2";
            headerList[18].location = location;
            headerList[18].index = numElementResults_;

            numElementResults_++;
            headerList[19].name="Strain dir cosine l3";
            headerList[19].location = location;
            headerList[19].index = numElementResults_;

            numElementResults_++;
            headerList[20].name="Strain dir cosine m3";
            headerList[20].location = location;
            headerList[20].index = numElementResults_;

            numElementResults_++;
            headerList[21].name="Strain dir cosine n3";
            headerList[21].location = location;
            headerList[21].index = numElementResults_;

      }
      else if (ID.compare("es")==0) {
            //cout << "Element stresses" << endl;
            ncols = 26;
            headerList.resize(ncols);

            numElementResults_++;
            headerList[0].name="Stress sigma-xx";
            headerList[0].location = location;
            headerList[0].index = numElementResults_;

            numElementResults_++;
            headerList[1].name="Stress sigma-yy";
            headerList[1].location = location;
            headerList[1].index = numElementResults_;

            numElementResults_++;
            headerList[2].name="Stress sigma-zz";
            headerList[2].location = location;
            headerList[2].index = numElementResults_;

            numElementResults_++;
            headerList[3].name="Stress tau-xy";
            headerList[3].location = location;
            headerList[3].index = numElementResults_;

            numElementResults_++;
            headerList[4].name="Stress tau-yz";
            headerList[4].location = location;
            headerList[4].index = numElementResults_;

            numElementResults_++;
            headerList[5].name="Stress tau-xz";
            headerList[5].location = location;
            headerList[5].index = numElementResults_;

            numElementResults_++;
            headerList[6].name="Stress work density";
            headerList[6].location = location;
            headerList[6].index = numElementResults_;

            numElementResults_++;
            headerList[7].name="Stress von mises stress";
            headerList[7].location = location;
            headerList[7].index = numElementResults_;

            numElementResults_++;
            headerList[8].name="Stress C1-Mat val (m strain)";
            headerList[8].location = location;
            headerList[8].index = numElementResults_;

            numElementResults_++;
            headerList[9].name="Stress C2-Mat val (m stress)";
            headerList[9].location = location;
            headerList[9].index = numElementResults_;

            numElementResults_++;
            headerList[10].name="Stress C3-Mat val (f)";
            headerList[10].location = location;
            headerList[10].index = numElementResults_;

            numElementResults_++;
            headerList[11].name="Stress invariant 1";
            headerList[11].location = location;
            headerList[11].index = numElementResults_;

            numElementResults_++;
            headerList[12].name="Stress invariant 2";
            headerList[12].location = location;
            headerList[12].index = numElementResults_;

            numElementResults_++;
            headerList[13].name="Stress invariant 3";
            headerList[13].location = location;
            headerList[13].index = numElementResults_;

            numElementResults_++;
            headerList[14].name="Stress principal (Minimum)";
            headerList[14].location = location;
            headerList[14].index = numElementResults_;

            numElementResults_++;
            headerList[15].name="Stress principal (Inter)";
            headerList[15].location = location;
            headerList[15].index = numElementResults_;

            numElementResults_++;
            headerList[16].name="Stress principal (Maximum)";
            headerList[16].location = location;
            headerList[16].index = numElementResults_;

            numElementResults_++;
            headerList[17].name="Stress dir cosine l1";
            headerList[17].location = location;
            headerList[17].index = numElementResults_;

            numElementResults_++;
            headerList[18].name="Stress dir cosine m1";
            headerList[18].location = location;
            headerList[18].index = numElementResults_;

            numElementResults_++;
            headerList[19].name="Stress dir cosine n1";
            headerList[19].location = location;
            headerList[19].index = numElementResults_;

            numElementResults_++;
            headerList[20].name="Stress dir cosine l2";
            headerList[20].location = location;
            headerList[20].index = numElementResults_;

            numElementResults_++;
            headerList[21].name="Stress dir cosine m2";
            headerList[21].location = location;
            headerList[21].index = numElementResults_;

            numElementResults_++;
            headerList[22].name="Stress dir cosine n2";
            headerList[22].location = location;
            headerList[22].index = numElementResults_;

            numElementResults_++;
            headerList[23].name="Stress dir cosine l3";
            headerList[23].location = location;
            headerList[23].index = numElementResults_;

            numElementResults_++;
            headerList[24].name="Stress dir cosine m3";
            headerList[24].location = location;
            headerList[24].index = numElementResults_;

            numElementResults_++;
            headerList[25].name="Stress dir cosine n3";
            headerList[25].location = location;
            headerList[25].index = numElementResults_;

      }
      else if (ID.compare("et")==0) {
            ncols = 3;
            headerList.resize(ncols);

            numElementResults_++;
            headerList[0].name="Temp abs";
            headerList[0].location = location;
            headerList[0].index = numElementResults_;

            numElementResults_++;
            headerList[1].name="Temp delta";
            headerList[1].location = location;
            headerList[1].index = numElementResults_;

      }
      else if (ID.compare("ne")==0) {
            //cout << "Nodal strains" << endl;
            ncols = 22;
            headerList.resize(ncols);

            numNodeResults_++;
            headerList[0].name="Strain eps-xx";
            headerList[0].location = location;
            headerList[0].index = numNodeResults_;

            numNodeResults_++;
            headerList[1].name="Strain eps-yy";
            headerList[1].location = location;
            headerList[1].index = numNodeResults_;

            numNodeResults_++;
            headerList[2].name="Strain eps-zz";
            headerList[2].location = location;
            headerList[2].index = numNodeResults_;

            numNodeResults_++;
            headerList[3].name="Strain gam-xy";
            headerList[3].location = location;
            headerList[3].index = numNodeResults_;

            numNodeResults_++;
            headerList[4].name="Strain gam-yz";
            headerList[4].location = location;
            headerList[4].index = numNodeResults_;

            numNodeResults_++;
            headerList[5].name="Strain gam-xz";
            headerList[5].location = location;
            headerList[5].index = numNodeResults_;

            numNodeResults_++;
            headerList[6].name="Strain Effective mises strain";
            headerList[6].location = location;
            headerList[6].index = numNodeResults_;

            numNodeResults_++;
            headerList[7].name="Strain invariant 1";
            headerList[7].location = location;
            headerList[7].index = numNodeResults_;

            numNodeResults_++;
            headerList[8].name="Strain invariant 2";
            headerList[8].location = location;
            headerList[8].index = numNodeResults_;

            numNodeResults_++;
            headerList[9].name="Strain invariant 3";
            headerList[9].location = location;
            headerList[9].index = numNodeResults_;

            numNodeResults_++;
            headerList[10].name="Strain principal (Minimum)";
            headerList[10].location = location;
            headerList[10].index = numNodeResults_;

            numNodeResults_++;
            headerList[11].name="Strain principal (Inter)";
            headerList[11].location = location;
            headerList[11].index = numNodeResults_;

            numNodeResults_++;
            headerList[12].name="Strain principal (Maximum)";
            headerList[12].location = location;
            headerList[12].index = numNodeResults_;

            numNodeResults_++;
            headerList[13].name="Strain dir cosine l1";
            headerList[13].location = location;
            headerList[13].index = numNodeResults_;

            numNodeResults_++;
            headerList[14].name="Strain dir cosine m1";
            headerList[14].location = location;
            headerList[14].index = numNodeResults_;

            numNodeResults_++;
            headerList[15].name="Strain dir cosine n1";
            headerList[15].location = location;
            headerList[15].index = numNodeResults_;

            numNodeResults_++;
            headerList[16].name="Strain dir cosine l2";
            headerList[16].location = location;
            headerList[16].index = numNodeResults_;

            numNodeResults_++;
            headerList[17].name="Strain dir cosine m2";
            headerList[17].location = location;
            headerList[17].index = numNodeResults_;

            numNodeResults_++;
            headerList[18].name="Strain dir cosine n2";
            headerList[18].location = location;
            headerList[18].index = numNodeResults_;

            numNodeResults_++;
            headerList[19].name="Strain dir cosine l3";
            headerList[19].location = location;
            headerList[19].index = numNodeResults_;

            numNodeResults_++;
            headerList[20].name="Strain dir cosine m3";
            headerList[20].location = location;
            headerList[20].index = numNodeResults_;

            numNodeResults_++;
            headerList[21].name="Strain dir cosine n3";
            headerList[21].location = location;
            headerList[21].index = numNodeResults_;


      }
      else if (ID.compare("ns")==0) {
            //cout << "Nodal stresses" << endl;
            ncols = 26;
            headerList.resize(ncols);

            numNodeResults_++;
            headerList[0].name="Stress sigma-xx";
            headerList[0].location = location;
            headerList[0].index = numNodeResults_;

            numNodeResults_++;
            headerList[1].name="Stress sigma-yy";
            headerList[1].location = location;
            headerList[1].index = numNodeResults_;

            numNodeResults_++;
            headerList[2].name="Stress sigma-zz";
            headerList[2].location = location;
            headerList[2].index = numNodeResults_;

            numNodeResults_++;
            headerList[3].name="Stress tau-xy";
            headerList[3].location = location;
            headerList[3].index = numNodeResults_;

            numNodeResults_++;
            headerList[4].name="Stress tau-yz";
            headerList[4].location = location;
            headerList[4].index = numNodeResults_;

            numNodeResults_++;
            headerList[5].name="Stress tau-xz";
            headerList[5].location = location;
            headerList[5].index = numNodeResults_;

            numNodeResults_++;
            headerList[6].name="Stress work density";
            headerList[6].location = location;
            headerList[6].index = numNodeResults_;

            numNodeResults_++;
            headerList[7].name="Stress von mises stress";
            headerList[7].location = location;
            headerList[7].index = numNodeResults_;

            numNodeResults_++;
            headerList[8].name="Stress C1-Mat val (m strain)";
            headerList[8].location = location;
            headerList[8].index = numNodeResults_;

            numNodeResults_++;
            headerList[9].name="Stress C2-Mat val (m stress)";
            headerList[9].location = location;
            headerList[9].index = numNodeResults_;

            numNodeResults_++;
            headerList[10].name="Stress C3-Mat val (f)";
            headerList[10].location = location;
            headerList[10].index = numNodeResults_;

            numNodeResults_++;
            headerList[11].name="Stress invariant 1";
            headerList[11].location = location;
            headerList[11].index = numNodeResults_;

            numNodeResults_++;
            headerList[12].name="Stress invariant 2";
            headerList[12].location = location;
            headerList[12].index = numNodeResults_;

            numNodeResults_++;
            headerList[13].name="Stress invariant 3";
            headerList[13].location = location;
            headerList[13].index = numNodeResults_;

            numNodeResults_++;
            headerList[14].name="Stress principal (Minimum)";
            headerList[14].location = location;
            headerList[14].index = numNodeResults_;

            numNodeResults_++;
            headerList[15].name="Stress principal (Inter)";
            headerList[15].location = location;
            headerList[15].index = numNodeResults_;

            numNodeResults_++;
            headerList[16].name="Stress principal (Maximum)";
            headerList[16].location = location;
            headerList[16].index = numNodeResults_;

            numNodeResults_++;
            headerList[17].name="Stress dir cosine l1";
            headerList[17].location = location;
            headerList[17].index = numNodeResults_;

            numNodeResults_++;
            headerList[18].name="Stress dir cosine m1";
            headerList[18].location = location;
            headerList[18].index = numNodeResults_;

            numNodeResults_++;
            headerList[19].name="Stress dir cosine n1";
            headerList[19].location = location;
            headerList[19].index = numNodeResults_;

            numNodeResults_++;
            headerList[20].name="Stress dir cosine l2";
            headerList[20].location = location;
            headerList[20].index = numNodeResults_;

            numNodeResults_++;
            headerList[21].name="Stress dir cosine m2";
            headerList[21].location = location;
            headerList[21].index = numNodeResults_;

            numNodeResults_++;
            headerList[22].name="Stress dir cosine n2";
            headerList[22].location = location;
            headerList[22].index = numNodeResults_;

            numNodeResults_++;
            headerList[23].name="Stress dir cosine l3";
            headerList[23].location = location;
            headerList[23].index = numNodeResults_;

            numNodeResults_++;
            headerList[24].name="Stress dir cosine m3";
            headerList[24].location = location;
            headerList[24].index = numNodeResults_;

            numNodeResults_++;
            headerList[25].name="Stress dir cosine n3";
            headerList[25].location = location;
            headerList[25].index = numNodeResults_;            

      }
      else if (ID.compare("nr")==0) {
            //cout << "Nodal reactions" << endl;
            ncols = 3;
            headerList.resize(ncols);

            numNodeResults_++;
            headerList[0].name="Int Forces Transx";
            headerList[0].location = location;
            headerList[0].index = numNodeResults_;

            numNodeResults_++;
            headerList[1].name="Int Forces Transy";
            headerList[1].location = location;
            headerList[1].index = numNodeResults_;

            numNodeResults_++;
            headerList[2].name="Int Forces Transz";
            headerList[2].location = location;
            headerList[2].index = numNodeResults_;


      }
      else if (ID.compare("nv")==0) {
            //cout << "Nodal velocities" << endl;
            ncols = 3;
            headerList.resize(ncols);

            numNodeResults_++;
            headerList[0].name="Vels Transx";
            headerList[0].location = location;
            headerList[0].index = numNodeResults_;

            numNodeResults_++;
            headerList[1].name="Vels Transy";
            headerList[1].location = location;
            headerList[1].index = numNodeResults_;

            numNodeResults_++;
            headerList[2].name="Vels Transz";
            headerList[2].location = location;
            headerList[2].index = numNodeResults_;

      }
      else if (ID.compare("nt")==0) {
            ncols = 3;
            headerList.resize(ncols);

            numNodeResults_++;
            headerList[0].name="Temp abs";
            headerList[0].location = location;
            headerList[0].index = numNodeResults_;

            numNodeResults_++;
            headerList[1].name="Temp delta";
            headerList[1].location = location;
            headerList[1].index = numNodeResults_;

      }
      else {
            cerr << "Error: Unknown (non-standard) result file type " << ID << "." << endl;
            exit(1);
      }

      // Actually insert into our map
      int i = 1;
      for (vector<ResultHeader>::iterator it = headerList.begin(); it!=headerList.end(); ++it)
      {
            headerMap_.insert(pair<pair<string,int>,ResultHeader> ( pair<string,int> (ID, i), *it) ) ;
            i++;
      }


}

// Parse a template file.
void Translator::parseTemplate(std::string file,std::string ID) {
      // Determine type of field from the 1st character of the ID
      FieldLocation location;
      if (ID[0]=='e' ) {
            location = element;
      }
      else if (ID[0]=='n') {
            location = nodal;
      }
      else {
            cerr << "Error: Unknown field type." << endl;
            exit(1);
      }
      
      // Open, etc.
      cout << "Opening file " << file << "." << endl;
      ifstream templateFile;
      templateFile.open(file.c_str(),ifstream::in);
      if(!templateFile) {
            cerr << "Error: file could not be opened." << endl;
            exit(1);
      }

      // Read through and parse into our results structure
      char linebuffer[256];
      string stringbuffer;
      while (!templateFile.eof()) {
            templateFile.getline(linebuffer,256);
            // Look for TYPE
            stringbuffer = strip(string(linebuffer));

            if (stringbuffer.substr(0,5).compare("TYPE=")==0 && stringbuffer.compare("TYPE=END")!=0) {
                  // Look for column and cast
                  vector<int> columns;
                  templateFile.getline(linebuffer,256);
                  stringbuffer = strip(string(linebuffer));
                  if (stringbuffer.substr(0,6).compare("COLUMN")==0) {
                        stringbuffer.erase(0,7);
                        // Have CSV, split into columns
                        istringstream ss(stringbuffer);
                        while (ss) {
                              int i;
                              string s;
                              if (!getline(ss,s,',')) break;
                              istringstream coni(s);
                              coni >> i;
                              columns.push_back(i);
                        }

                  }
                  else {
                        cerr << "Error: Expected COLUMN" << endl;
                        exit(1);
                  }
                  
                  // Primary
                  templateFile.getline(linebuffer,256);
                  stringbuffer = string(linebuffer);
                  string stripbuffer = strip(stringbuffer);
                  string primary;
                  if (stripbuffer.substr(0,4).compare("PRI=")==0) {
                        // Find the equals sign
                        size_t loc = stringbuffer.find("=");
                        primary = stringbuffer.erase(0,loc+1);
                        // Delete leading spaces
                        if (primary[0]==' ') primary.erase(0,1);

                  }
                  else {
                        cerr << "Error: Expected PRI" << endl;
                        exit(1);
                  }
                  
                  // Look for secondary
                  templateFile.getline(linebuffer,256);
                  stringbuffer = string(linebuffer);
                  stripbuffer = strip(stringbuffer);
                  string secondary;
                  if (stripbuffer.substr(0,4).compare("SEC=")==0) {
                        // Find the equals sign
                        size_t loc = stringbuffer.find("=");
                        secondary = stringbuffer.erase(0,loc+1);
                        // Delete leading spaces
                        if (secondary[0]==' ') secondary.erase(0,1);
                  }
                  else {
                        cerr << "Error: Expected SEC" << endl;
                        exit(1);
                  }

                  // Concatate
                  string name = primary + " - " + secondary;

                  // Insert to our header map model
                  for (int i=0; i<columns.size(); i++) {
                        ResultHeader header;
                        if (location==nodal) {
                              numNodeResults_++;
                              header.location=nodal;
                              header.index=numNodeResults_;
                        }
                        else if (location==element) {
                              numElementResults_++;
                              header.location=element;
                              header.index=numElementResults_;
                        }
                        string mname = name;
                        if (columns.size()>1) {
                              if (i==0) {
                                    mname = mname+"x"; 
                              }
                              else if (i==1) {
                                    mname = mname+"y";
                              }
                              else if (i==2) {
                                    mname = mname+"z";
                              }
                              else {
                                    cerr << "Error: Tensor data in result file." << endl;
                                    exit(1);
                              }

                        }
                        header.name = mname;

                        // Insert to our map
                        headerMap_.insert(pair<pair<string,int>,ResultHeader> ( pair<string,int> (ID, columns[i]), header) ) ;
                  }


            }
            else {
            }

      }

      templateFile.close();

}

