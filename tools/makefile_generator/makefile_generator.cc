 /***************************************************************************
 
                           Makefile Generator Program
                                                     
    File: makefile_generator.cc             
                                                               
    Programmer: Erick Jenkins
    Date: Sept 2004
                                                               
    Description:
      Creates a functional makefile from a header file, table file,
      and object file.

    Revision History:
      Sept 2004 - Fixed to support Linux header files and
                  supports compiling via g++ for Windows/Linux
      Jan 2006  - Changed DVF to Windows w/ ifort and MKS toolkit
                  Added AMD64 running Linux
                  Added Windows w/ ifort and Cygwin toolkit
                                                                
  ***************************************************************************/


#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <ctype.h>
#include <cassert>
#include <cmath>
#include <iomanip>

// Program-Wide Constants

#define max_files 1000
//#define max_platforms 10
#define max_dependencies 30
#define max_objects 1000
#define max_object_lines 20
#define max_custom_lines 20
#define num_platforms 10
int debug = false;

using namespace std;

// Function Prototypes

void WriteMakeFile(string, string, string[max_files], 
                   string[max_files][max_dependencies],
                   string[max_files][num_platforms], int, int,
                   string[max_objects], int);

void MakeMenu(string &, string &, string &, int &);
void MakeMenu(string, string &);
void ReadObjectFile(string[max_objects], int &);

void ReadTableFile(string[max_files], string[max_files][max_dependencies],
                   string[max_files][num_platforms], int &);

void WriteMakeFile(string, string);
void GenerateOptimizations(int, int, int, string[max_files][num_platforms], bool);
void GenerateFileNames(string &, string &, string &, int);
void GenerateFileNames(string, string);



// ******************************************************************************



int main(int argc, char *argv[])
{
  string filelist[max_files];
  string dependencies[max_files][max_dependencies];
  string objects[max_objects];
  string optimizations[max_files][num_platforms];
  int numfiles=0, num_objects=0, menuselect=0;
  string headerfilename;
  string makefilename;
  string filenameprefix;
  string maketype="";

  cout << "\n\n\n";
  cout << " *************************************************************************************\n";
  cout << "                                WARP3D Makefile Generator                             \n";
  cout << " *************************************************************************************\n";
  cout << "\n\n";
	
  //  Display command-line argument help if necessary
  if(argc > 1)
    {
      if(string(argv[1])=="-help" || string(argv[1])=="--help")
	{
	  cout << "  Usage: makegen.exe [MACHINE] [VERSION] \n\n";
	  cout << "    [MACHINE] - Choose from :" << endl;
	  cout << "     \"hpux11\" \"hpuxit\" \"dec\" \"win_mks\" \"win_cyg\" \"ibm\" " << endl;
          cout << "     \"sgi\" \"sga\" \"linux\" \"lin64\" OR \"all\" " << endl << endl;
          cout << "    Where..." << endl << endl;
          cout << "     \"hpux11\"  Creates a Makefile for HP PA-RISC Machines running HP-UX 11.xx" << endl;
	  cout << "     \"hpuxit\"  Creates a Makefile for HP Itanium-2+ Machines running HP-UX 11.xx" << endl;
	  cout << "     \"dec\"     Creates a Makefile for DEC Alpha Machines running OSF" << endl;
	  cout << "     \"win_mks\" Creates a Makefile for 32-bit PCs running Win2k/WinXP"<< endl;
 	  cout << "                   w/ Intel Fortran Compiler and MKS Toolkit" << endl;
          cout << "     \"win_cyg\" Creates a Makefile for 32-bit PCs running Win2k/WinXP " << endl;
	  cout << "                   w/ Intel Fortran Compiler and Cygwin Toolkit" << endl;
          cout << "     \"ibm\"     Creates a Makefile for IBM RISC6000, SP-2, or SP-3 Machines running AIX" << endl;
	  cout << "     \"sgi\"     Creates a Makefile for SGI Octane/Origin Machines running IRIX" << endl;
	  cout << "     \"sga\"     Creates a Makefile for SGI Altix Machines Running SGI Linux" << endl;
	  cout << "     \"linux\"   Creates a Makefile for IA-32 Machines running Linux" << endl;
  	  cout << "     \"lin64\"   Creates a Makefile for AMD64 or Intel EM64T machines running Linux" << endl;
          cout << "     \"all\"     Creates all of the above Makefiles" << endl << endl;
	  cout << "  For all selections except \"win_mks\", \"win_cyg\", and \"all\", enter a second argument:\n\n";
	  cout << "    [VERSION] - Choose from :" << endl;
	  cout << "     \"mpi\" \"ser\" OR \"both\" " << endl << endl;
          cout << "    Where..." << endl << endl;
	  cout << "     \"mpi\"     Creates a Message Passing Interface version of the Makefile" << endl;
	  cout << "     \"ser\"     Creates a serial (non-parallel) version of the Makefile" << endl;
	  cout << "     \"both\"    Creates both of the above versions of the Makefile" << endl << endl;
	  cout << "  Example: \"makegen.exe dec mpi\" OR \"makegen.exe all\" ...etc..." << endl;
	  goto endprogram;
	}
    }
	
  //  Process objects file
  ReadObjectFile(objects, num_objects);
	
  //  Process table file
  ReadTableFile(filelist, dependencies, optimizations, numfiles);
	
  //  Process command-line arguments if necessary
  if(argc>1)
    {
      menuselect = 0;
      if(string(argv[1])=="hpux11")  menuselect = -1;
      if(string(argv[1])=="hpuxit")  menuselect = -2;
      if(string(argv[1])=="dec")     menuselect = -3;
      if(string(argv[1])=="win_mks") menuselect = -4;
      if(string(argv[1])=="win_cyg") menuselect = -5;
      if(string(argv[1])=="ibm")     menuselect = -6;
      if(string(argv[1])=="sgi")     menuselect = -7;
      if(string(argv[1])=="linux")   menuselect = -8;
      if(string(argv[1])=="sga")     menuselect = -9;
      if(string(argv[1])=="lin64")   menuselect = -10;
      if(string(argv[1])=="all")     menuselect = -11;
		
      if(argc>2)
	{
	  if(string(argv[2])=="ser")  maketype="s";
	  if(string(argv[2])=="mpi")  maketype="m";
	  if(string(argv[2])=="both") maketype="b";
	}
      delete argv;
    }
	
  while(true)
    {
      //  Run processor selection menu if necessary
      if(menuselect>=0)
	{
	  cout << endl << endl;
	  MakeMenu(headerfilename, makefilename, filenameprefix, menuselect);
	}
      if(menuselect==0) goto endprogram;
		
      //  All of the above case:
      else if(abs(menuselect)==11)
	{
	  for(int i=1; i<=num_platforms; i++)
	    {
	      cout << endl << endl;
	      GenerateFileNames(headerfilename, makefilename, filenameprefix, i);
	      WriteMakeFile(headerfilename, makefilename, filelist, dependencies, optimizations,
                            i, numfiles, objects, num_objects); 
	      if(filenameprefix!="") GenerateFileNames(filenameprefix, "b");
	    }
	}
		
      //  Specific processor case:
      else
	{
	  cout << endl << endl;
	  GenerateFileNames(headerfilename, makefilename, filenameprefix, abs(menuselect));
	  WriteMakeFile(headerfilename, makefilename, filelist, dependencies, optimizations,
                        abs(menuselect), numfiles, objects, num_objects); 
	  //  cout << endl;
	  //  Prompt for .ser or .mpi file (if necessary)
	  if(filenameprefix!="")
	    {
	      if(menuselect>=0)
		{
		  cout << endl;
		  MakeMenu(filenameprefix, maketype);
		}
	      if(maketype!="") 
		{
		  GenerateFileNames(filenameprefix, maketype);
		}
	      //  cout << endl;
	    }
	  maketype="";
	}
      // Avoid menu if necessary
      if(menuselect<0) goto endprogram;
    }

  endprogram :

  cout << "\n\n";
  cout << " *************************************************************************************\n";
  cout << "\n\n\n";



  // Debug - output optimization result
  if (debug)
    { 
      cout << "          DEBUG:  Optimization Matrix: \n";
      for ( int i = 0 ; i < max_files ; i++)
	{ 
          cout << "            Row " << i << endl;
	  for ( int j = 0 ; j < num_platforms ; j++ ) { cout << optimizations[i][j] << endl; }
	  cout << endl;
	}
      cout << endl << endl;
    }


  return(0);
}



// ******************************************************************************



void GenerateOptimizations(int row, int col, int opt, 
                           string optimizations[max_files][num_platforms], 
                           bool MPISpecial)

//   Create scenarios in an ordered switch based upon the COLUMN of each
//   machine type in the 'tablefile.txt' file.

{
  switch(col)
    {
    case 0 :  // HPUX11 optimizations
      if (debug) cout << "          DEBUG: HPUX11=" << opt << "\n";
      switch(opt)
	{
	case -1 : optimizations[row][col] = ""; break;
	case 0 : optimizations[row][col] = "+O0"; break;
	case 1 : optimizations[row][col] = "+O1"; break;
	case 2 : optimizations[row][col] = "+O2"; break;
	case 3 : optimizations[row][col] = "+O3"; break;
	case 4 : optimizations[row][col] = "+Oall"; break;
	case 5 : optimizations[row][col] = "+Oall"; break;
	}
      break;
    case 1 :  // HPUXIT optimizations
      if (debug) cout << "          DEBUG: HPUXIT=" << opt << "\n";
      switch(opt)
	{
	case -1 : optimizations[row][col] = ""; break;
	case 0 : optimizations[row][col] = "+O0"; break;
	case 1 : optimizations[row][col] = "+O1"; break;
	case 2 : optimizations[row][col] = "+O2"; break;
	case 3 : optimizations[row][col] = "+O2 +Oloop_block +Ovectorize +Onolimit"; break;
	case 4 : optimizations[row][col] = "+O3 +Oloop_block +Ovectorize +Onolimit"; break;
	case 5 : optimizations[row][col] = "+O3 +Oloop_block +Ovectorize +Onolimit"; break;
	}
      break;
    case 2 :  // DEC optimizations
      if (debug) cout << "          DEBUG: DEC=" << opt << "\n";
      switch(opt)
	{
	case -1 : optimizations[row][col] = ""; break;
	case 0 : optimizations[row][col] = "-O0"; break;
	case 1 : optimizations[row][col] = "-O1"; break;
	case 2 : optimizations[row][col] = "-O2 -pipeline -transform_loops"; break;
	case 3 : optimizations[row][col] = "-03 -pipeline -transform_loops"; break;
	case 4 : optimizations[row][col] = "-04 -pipeline -transform_loops"; break;
	case 5 : optimizations[row][col] = "-04 -pipeline -transform_loops"; break;
	}
      break;
    case 3 :  // Windows w/ MKS optimizations
      if (debug) cout << "          DEBUG: WIN_MKS=" << opt << "\n";
      switch(opt)
	{
	case -1 : optimizations[row][col] = ""; break;
	case 0 : optimizations[row][col] = "-optimize:0"; break;
	case 1 : optimizations[row][col] = "-optimize:1"; break;
	case 2 : optimizations[row][col] = "-optimize:2"; break;
	case 3 : optimizations[row][col] = "-optimize:3"; break;
	case 4 : optimizations[row][col] = "-optimize:3"; break;
	case 5 : optimizations[row][col] = "-optimize:3"; break;
	}
      break;
    case 4 :  // Windows w/ Cygwin optimizations
      if (debug) cout << "          DEBUG: WIN_CYGWIN=" << opt << "\n";
      switch(opt)
	{
	case -1 : optimizations[row][col] = ""; break;
	case 0 : optimizations[row][col] = "-O0 -Qip"; break;
	case 1 : optimizations[row][col] = "-O1 -Qip"; break;
	case 2 : optimizations[row][col] = "-O2 -Qip"; break;
	case 3 : optimizations[row][col] = "-O3 -Qip"; break;
	case 4 : optimizations[row][col] = "-O3 -Qip"; break;
	case 5 : optimizations[row][col] = "-O3 -Qip"; break;
	}
      break;
    case 5 :  // IBM optimizations
      if (debug) cout << "          DEBUG: IBM=" << opt << "\n";
      switch(opt)
	{
	case -1 : optimizations[row][col] = ""; break;
	case 0 : optimizations[row][col] = "-g"; break;
	case 1 : optimizations[row][col] = "-O3 -qMAXMEM=-1 -Q"; break;
	case 2 : optimizations[row][col] = "-O3 -qMAXMEM=-1 -Q"; break;
	case 3 : optimizations[row][col] = "-O3 -qMAXMEM=-1 -Q"; break;
	case 4 : optimizations[row][col] = "-O3 -qMAXMEM=-1 -Q"; break;
	case 5 : optimizations[row][col] = "-O3 -qMAXMEM=-1 -Q"; break;
	}
      break;
    case 6 :  // SGI optimizations
      if (debug) cout << "          DEBUG: SGI=" << opt << "\n";
      switch(opt)
	{
	case -1 : optimizations[row][col] = ""; break;
	case 0 : optimizations[row][col] = "-O0"; break;
	case 1 : optimizations[row][col] = "-O1"; break;
	case 2 : optimizations[row][col] = "-O2"; break;
	case 3 : optimizations[row][col] = "-O3"; break;
	case 4 : optimizations[row][col] = "-O3"; break;
	case 5 : optimizations[row][col] = "-O3"; break;
	}
      break;
    case 7 :  // IA-32 LINUX optimizations
      if (debug) cout << "          DEBUG:  LINUX=" << opt << "\n";
      switch(opt)
	{
	case -1 : optimizations[row][col] = ""; break;
	case 0 : optimizations[row][col] = "-O0 -ip"; break;
	case 1 : optimizations[row][col] = "-O1 -ip"; break;
	case 2 : optimizations[row][col] = "-O2 -ip"; break;
	case 3 : optimizations[row][col] = "-O3 -ip"; break;
	case 4 : optimizations[row][col] = "-O3 -ip"; break;
	case 5 : optimizations[row][col] = "-O3 -ip"; break;
	}
      break;
    case 8 :  // SGA optimizations
      if (debug) cout << "          DEBUG:  SGA=" << opt << "\n";
      switch(opt)
	{
	case -1 : optimizations[row][col] = ""; break;
	case 0 : optimizations[row][col] = "-O0 -ip"; break;
	case 1 : optimizations[row][col] = "-O1 -ip"; break;
	case 2 : optimizations[row][col] = "-O2 -ip"; break;
	case 3 : optimizations[row][col] = "-O3 -ip"; break;
	case 4 : optimizations[row][col] = "-O3 -ip"; break;
	case 5 : optimizations[row][col] = "-O3 -ip"; break;
	}
      break;
    case 9 :  // AMD64 Linux optimizations
      if (debug) cout << "          DEBUG:  L64=" << opt << "\n";
      switch(opt)
	{
	case -1 : optimizations[row][col] = ""; break;
	case 0 : optimizations[row][col] = "-O0 -ip"; break;
	case 1 : optimizations[row][col] = "-O1 -ip"; break;
	case 2 : optimizations[row][col] = "-O2 -ip"; break;
	case 3 : optimizations[row][col] = "-O3 -ip"; break;
	case 4 : optimizations[row][col] = "-O3 -ip"; break;
	case 5 : optimizations[row][col] = "-O3 -ip"; break;
	}
      break;
    }
 
  if(MPISpecial==1)  // Insert Marker for special MPI variable cases!
    {
      optimizations[row][col].insert(0, "[MPI]");
    }

  // End of function
}



// ******************************************************************************



void MakeMenu(string & headerfilename, string & makefilename, string & filenameprefix,
              int & menuselect)
{
  string selection;
	
  cout << " Select target machine for makefile creation: \n\n";
  cout << "  0.  Exit the Program" << endl;
  cout << "  1.  HP PA-RISC Machines running HP-UX 11.xx" << endl;
  cout << "  2.  HP Itanium-2+ Machines running HP-UX 11.xx" << endl;
  cout << "  3.  DEC Alpha Machines running OSF" << endl;
  cout << "  4.  IA-32 Machines running Windows 2000/XP w/ Intel Fortran Compiler and MKS Toolkit" << endl;
  cout << "  5.  IA-32 Machines running Windows 2000/XP w/ Intel Fortran Compiler and Cygwin Toolkit" << endl;
  cout << "  6.  IBM RISC6000, SP-2, or SP-3 Machines running AIX" << endl;
  cout << "  7.  SGI Octane/Origin Machines running IRIX" << endl;
  cout << "  8.  IA-32 or AMD Machines running Linux" << endl;
  cout << "  9.  SGI Altix Machines Running SGI Linux" << endl;
  cout << "  10. AMD64 or Intel EM64T Machines running Linux" << endl;
  cout << "  11. Create all of the above Makefiles" << endl << endl;

  while(true)
    {
      cout << " Make a selection from above: ";
      cin >> selection;
      menuselect=atoi(selection.c_str());
      if((menuselect<=10 && menuselect>=0) || menuselect==11) break;
      // cout << "Choose again!\n";
    }
  // cout << endl;
}



// ******************************************************************************



void GenerateFileNames(string & headerfilename, string & makefilename,
                       string & filenameprefix, int menuselect)
{
  switch(menuselect)
    {
    case 0 : 
      break;
    case 1 : 
      headerfilename = "Header.hpux11.inc";
      makefilename = "Makefile.hpux11.inc";
      // Use next line ONLY when .ser or .mpi files are desired - otherwise leave blank.
      filenameprefix = "hpux11";
      break;
    case 2 : 
      headerfilename = "Header.hpuxit.inc";
      makefilename = "Makefile.hpuxit.inc";
      filenameprefix = "hpuxit";
      break;
    case 3 : 
      headerfilename = "Header.dec.inc";
      makefilename = "Makefile.dec.inc";
      filenameprefix = "dec";
      break;
    case 4 : 
      headerfilename = "Header.windows.mks";
      makefilename = "Makefile.windows.mks";
      filenameprefix = "";
      break;
    case 5 : 
      headerfilename = "Header.windows.cygwin";
      makefilename = "Makefile.windows.cygwin";
      filenameprefix = "";
      break;
    case 6 : 
      headerfilename = "Header.ibm.inc";
      makefilename = "Makefile.ibm.inc";
      filenameprefix = "ibm";
      break;
    case 7 : 
      headerfilename = "Header.sgi.inc";
      makefilename = "Makefile.sgi.inc";
      filenameprefix = "sgi";
      break;
    case 8 : 
      headerfilename = "Header.linux_ia32.inc";
      makefilename = "Makefile.linux_ia32.inc";
      filenameprefix = "linux_ia32";
      break;
    case 9 :
      headerfilename = "Header.sga_linux.inc";
      makefilename = "Makefile.sga_linux.inc";
      filenameprefix = "sga_linux";
      break;
    case 10 :
      headerfilename = "Header.linux_amd64.inc";
      makefilename = "Makefile.linux_amd64.inc";
      filenameprefix = "linux_amd64";
      break;
    }
}



// ******************************************************************************



void MakeMenu(string filenameprefix, string & maketype)
{
  do
    {
      cout << " Would you like a (s)erial makefile, an (m)pi makefile, or (b)oth? ";
      cin >> maketype;
    } while(maketype!="s" && maketype!="m" && maketype!="b");
  // cout << endl;
}



// ******************************************************************************




void GenerateFileNames(string filenameprefix, string maketype)
{
  string header_file="Header.";
  string output_file="Makefile.";
  header_file.append(filenameprefix);
  output_file.append(filenameprefix);
  if(maketype=="s" || maketype=="b")
    {
      cout << endl;
      WriteMakeFile(header_file+".ser", output_file+".ser");
    }
  if(maketype=="m" || maketype=="b")
    {
      cout << endl;
      WriteMakeFile(header_file+".mpi", output_file+".mpi");
    }
}



// ******************************************************************************



void ReadObjectFile(string objects[max_objects], int & num_objects)
{
  string objectfilename;
  ifstream object_file;
  char inputchar = ' ';
  char input[100];
  int char_position=0;
  string inputstring=" ";
	
  // Open the objects file, but report error if open fails
  while(true)
    {
      // cout << " Enter object file name : ";
      // cin >> objectfilename;
      objectfilename="objects.txt";
      object_file.open(objectfilename.c_str());
      if(object_file.fail()) { cout << "  Error: Could not open object file. Exiting Program...\n"; exit(0); }
      else break;
    }
	
  // Read in object file list...
  cout << " Reading " << objectfilename << "..." << flush;
  while(object_file.eof()!=1)  // Loop for each object file entry (one per line)
    {
      object_file.get(inputchar);
		
      if(object_file.eof()==1) break;
		
      if(inputchar=='#')
	{
	  object_file.getline(input, 80, '\n');
	}
      else if(inputchar=='\n') 
	{
	  num_objects++;
	  char_position=0;
	}
      else if(int(inputchar) != 10 && int(inputchar) != 13)
	{
	  inputstring = inputchar;
	  objects[num_objects].insert(char_position, inputstring);
	  char_position++;
	}

    }
	
  cout << num_objects << " objects total... " << flush;  // Output total number of objects
  object_file.close();
  cout << " Done!" << endl;
}



// ******************************************************************************



void ReadTableFile(string filelist[max_files], string dependencies[max_files][max_dependencies],
                   string optimizations[max_files][num_platforms], int & numfiles)
{
  string tablefilename;
  ifstream table_file;
  int numdeps=0;
  char input[200];
  char junkstuff[80];
	
  // Open the table file, but report error if open fails
  while(true)
    {
      tablefilename="tablefile.txt";
      table_file.open(tablefilename.c_str());
      if(table_file.fail()) { cout << "\n  Error: Could not open table file. Exiting Program...\n";  exit(0); }
      else break;
    }
	
  // Read in .f file list, parameters, and optimizations
  cout << " Reading " << tablefilename << "... " << flush;
  while(table_file.eof()!=1)  // Loop for each .f file entry (every two lines)
    {
      numdeps=0;
      do
	{
          // Read in .f filename (again using char array input, since string is failing...)
	  table_file.getline(input, 200, '\n');
	} 
      while(input[0]=='#');  // Ignore lines with # at the start of the file
      filelist[numfiles]=input;
      // Output the .f file to the screen
      if (debug) cout << "\n          DEBUG: Current .f file is: " << input <<  endl;
      if (debug) cout << "          DEBUG: Variable 'numfiles' = " << numfiles << endl;
  
      filelist[numfiles].erase(filelist[numfiles].find("."));
		
      // Read optimizations		
      int opt;
      bool MPISpecial;

      for (int g=0; g<num_platforms; g++)  // Loop to read each optimization entry
	{
          // Read each of the optimization entries, one at a time
	  if(g<(num_platforms-1)) table_file.getline(input, 200, ' '); 
	  else table_file.getline(input, 200, '\n');
	  opt=atoi(input);
          // Output the optimization number
	  if (debug) cout << "          DEBUG: Optimization number is: " << opt << endl;
	  if(input[1]=='*') MPISpecial=1;
	  else MPISpecial=0;
          // Send optimization integer to be converted to string value
	  GenerateOptimizations(numfiles, g, opt, optimizations, MPISpecial);  
	}

      // Read Dependencies
      do
	{
	  do
	    {
              // Read in dependencies (again using char array input, since string is failing...)
	      table_file.getline(input, 80, '\n');  
	    } while(input[0]=='#');
	  if(string(input).length()==0 || string(input).length()==1) break;
	  numdeps++;
	  dependencies[numfiles][numdeps]=string(input).substr(1, string(input).length());
	  int badpos = dependencies[numfiles][numdeps].find(char(13));
	  if(badpos>1) dependencies[numfiles][numdeps].replace(badpos, 1, "");
          
          // Output each dependency	
	  if (debug)
          {
	    cout << "  DEBUG: Dependency matrix is: \n";
	    cout << "   " << dependencies[numfiles][numdeps] << endl;
	  }	
	}
      while(true);


      if (debug) cout << "          DEBUG: POINT C... \n";

      if(numdeps>=10) dependencies[numfiles][0]=48+numdeps/10;
      dependencies[numfiles][0]+=(48+numdeps%10);
		
      numfiles++;  // Increment the number of files
    }
  table_file.close();
  cout << "                   Done!" << endl;
}



// ******************************************************************************



void WriteMakeFile(string header_id, string output_id, string filelist[max_files],
                   string dependencies[max_files][max_dependencies],
                   string optimizations[max_files][num_platforms],
                   int menuselect, int numfiles, string objects[max_objects], int num_objects)
{
  char text_line[100];
  string string_line="";
  string object_code[max_object_lines];
  int object_lines=0;
  string custom_code[max_custom_lines];
  int custom_lines=0;
  char input;
	
  ofstream output_file;
	
  ifstream header_file;
	
  // Open and begin writing the Makefile
  cout << " Reading " << header_id << endl;
  header_file.open(header_id.c_str());
  if(header_file.fail())
    {
      cout << "Error opening header file!\n";
      return;
    }
	
  cout << " Writing " << output_id << "... " << flush;
  output_file.open(output_id.c_str());
  while(header_file.eof() !=1)  // Begin reading header file
    {
      header_file.get(input);  // Read in a character, including '\n' characters
      // Test the input of characters
      if (debug) cout << "          DEBUG: Variable input = " << input << endl;
      if(input=='[')  // Process "[xxx]" headings
	{
	  header_file.getline(text_line, 80, ']');
			
	  // Fix for HP run-time error
	  // if(string(text_line)=="OBJECTS")
	  if(text_line[0]=='O' && text_line[1]=='B' && text_line[2]=='J' && text_line[6]=='S')
	    {
	      object_lines=0;
	      header_file.getline(text_line, 100, '\n');
              // Debug: test the getline
	      if (debug) cout << "          DEBUG: Variable 'text_line' = " << text_line << endl;
	      do       // Copy the OBJECTS code lines to an array of strings
		{
		  header_file.getline(text_line, 100, '\n');
					
		  // Fix for HP run-time error
		  // if(string(text_line)=="[END_OBJECTS]") break;
		  if(text_line[0]=='[' && text_line[1]=='E' && text_line[2]=='N' && text_line[12]==']') break;
		  object_code[object_lines]=text_line;
		  object_lines++;
		} while(true);
              // Process and write object code to the output file
	      for(int i=0; i<num_objects; i++)  
		{
		  for(int j=0; j<object_lines; j++)
		    {
		      string_line = object_code[j];
		      int position = string_line.find("[object]");
		      string fixedstring=objects[i];
						
		      // Fix for issue with outputting one too many characters
		      if(position>=0) string_line.replace(position, 8, objects[i]);
				
                      // Last object listing in file -- omit the final "\" character		
		      if(j==object_lines-1 && i==num_objects-1)  
			{
			  int pos2 = string_line.rfind('\\');
                          // Erase the \ character and the preceding space character
			  if(pos2>=0) string_line.erase(pos2-1); 
			}
		      output_file << string_line << endl;
		    }
		}
	    }
			
	  // Fix for HP run-time error
	  // else if(string(text_line)=="CUSTOM")
	  else if(text_line[0]=='C' && text_line[1]=='U' && text_line[2]=='S' && text_line[5]=='M')
	    {
	      custom_lines=0;
	      header_file.getline(text_line, 100, '\n');
	      do
		{
		  header_file.getline(text_line, 100, '\n');
					
		  // Fix for HP run-time error
		  // if(string(text_line)=="[END_CUSTOM]") break;
		  if(text_line[0]=='[' && text_line[1]=='E' && text_line[2]=='N' && text_line[11]==']') break;
		  custom_code[custom_lines]=text_line;
		  custom_lines++;
		} while(true);
	      for(int i=0; i<numfiles; i++)
		{
		  for(int j=0; j<custom_lines; j++)
		    {
		      string_line = custom_code[j];
		      string opts;
		      int position=0;
		      do
			{
			  position = string_line.find("[f-file]");  // [f-file] occurences
			  if(position>=0) string_line.replace(position, 8, filelist[i]);
			} while(position>=0);
						
		      position = string_line.find("[optimizations]");  // [optimizations] occurences
		      if(position>=0)
			{
			  if(optimizations[i][menuselect-1].find("[MPI]")==0)
			    {
			      // IF NECESSARY, replace $(F90) optimization with $(MPIF90)
			      optimizations[i][menuselect-1].replace(0, 5, "");
			      int mpipos=string_line.find("$(F90)");
			      string_line.replace(mpipos, 6, "$(MPIF90)");
			      position+=3;
			    }
			  string_line.replace(position, 15, optimizations[i][menuselect-1]);
			}
						
		      position = string_line.find("[dependencies]");  // [dependencies] occurences
		      if(position>=0)
			{
			  string deps = "";
			  for(int k=1; k<=atoi(dependencies[i][0].c_str()); k++)
			    {
                              // Start a new line for dependencies every 3 entries
			      if(k%2==0) deps.append("\\\n                   ");  
			      deps.append(dependencies[i][k]);
			      deps.append(" ");
			    }
			  string_line.replace(position, 14, deps);
			}
		      // Write the processed line to the output file		
		      output_file << string_line << endl;  
		    }
		}
	    }
	}
      else 
	{
	  output_file.put(input);
	}
    }
	
  header_file.close();
  output_file.close();
  cout << "             Done!" << endl;
}



// ******************************************************************************



void WriteMakeFile(string input_id, string output_id)
{
  ifstream input_file;
  ofstream output_file;
  char input;
	
  cout << " Reading " << input_id << endl;
  input_file.open(input_id.c_str());
  if(input_file.fail())
    {
      cout << "  Error: could not open " << input_id << endl;
      return;
    }
	
  output_file.open(output_id.c_str());
	
  cout << " Writing " << output_id << "... " << flush;
  while(true)
    {
      input_file.get(input);
      if(input_file.eof()==1) break;
      output_file.put(input);
    }
	
  input_file.close();
  output_file.close();
  cout << "             Done!" << endl;
}


// ******************************************************************************
