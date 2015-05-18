#include <math.h>
#include <stdlib.h>
#include <string>
#include <iomanip.h>
#include <fstream.h>

#include <iostream.h>

int results_type, fringe_type, num_records, num_cycles, btwn_cycles, fringe_style, results_format, load_step;
int records[10000]={0};
char tif_name[200]; 
char stress_results[50];
char template_path[200];
char structure[9];
int num_colors, format, frames;
double max_value;
double min_value;
char edge_option[30];
double scale_factor;
char scale[10];
char movie_name[100];

int debug=false;

int import( char file[], char import_flag[], char template_name[] );
int input( void );
void deformation( void );
int create_extension( int step, char extension[] );
void stress( void );
void strain( void );
void stress_and_deformation(void);
void strain_and_deformation( void );
int get_records( void );
void calculate_fringe( void );
void mpeg_stress( void );
void mpeg_strain( void );
void mpeg_deformation( int record );
void stress_options( void );
void strain_options( void );


// ---------------------------------------------------------------
// ANIMATION.C
//
// Purpose: This program calls supplementary functions to create
//          a PATRAN session file that creates .tif files or .mpg 
//          movies.
//
// Written by: Breanna Bailey
//
// Last modified: 8/28/01
//
// Global Variables Used:
//   1) results_type = switch variable to select which type of
//                     image file to create (deformation,etc.)
//   2) format = flag to indicate mpeg or tif format
//
// Functions Called:
//   1) input.c = separate file
//   2) deformation.c = separate file
//   3) stress.c = separate file
//   4) strain.c = separate file
//   5) stress_and_deformation.c = separate file
//   6) strain_and_deformation.c = separate file
//   7) mpeg_deformation.c = separate file
//   8) mpeg_stress.c = separate file
//   9) mpeg_strain.c = separate file
// ---------------------------------------------------------------

int main()
{

  //     Call user inputs.  Terminate program if unsuccessful.
  if ( input( ) )
    {
      cout << endl << flush;
      cout << "ERROR during input.  Program will terminate." << endl << flush;
      return(0);
    }

 
  //     Open output file.
  ofstream out(reinterpret_cast< const char *>(movie_name), ios::out );


//  ofstream out( "movie.ses.01", ios::out );
  if (!out)
    { 
      cout << "Output file could not be opened." << endl;
      cout << "Program will terminate." << endl << flush;
      exit(1);
    }

  
  //     Call major functions.
  if ( format == 1 )
    {
      switch( results_type )
	{
	case 1:
	  deformation();
	  break;
	case 2:
	  stress();
	  break;
	case 3:
	  strain();
	  break;
	case 4:
	  stress_and_deformation();
	  break;
	case 5:	 
	  strain_and_deformation();
	}
    }

  if ( format == 2 )
   {
     switch( results_type )
       {
       case 1:
	 mpeg_deformation(1);
	 break;
       case 2:
	 mpeg_stress();
      	 break;
       case 3:
	 mpeg_strain();
	 break;
       case 4:
	 mpeg_deformation(0);
         mpeg_stress();
         break;

       case 5:
	 mpeg_deformation(0);
	 mpeg_strain();
	}
    }

  //     Closing comments.
  out.close();
  cout << endl << "Program Completed." << endl << flush;
  cout << "If successful, session file stored in " << movie_name << "." << endl << endl << flush;
  return(0);
}


// ---------------------------------------------------------------
// CALCULATE_FRINGE.C
//
// Purpose: Calculates spectrum range and posts to output file.
//
// Written by: Breanna Bailey
//
// Last modified: 8/30/01
//
// Global Variables Used:
//   1) out = output file
//   2) max_value = maximum value of stress/strain
//   3) min_value = minimum value of stress/strain
//   4) num_colors = number of colors in spectrum; number of 
//                   subranges is one fewer
//
// Local Variables Used:
//   1) i = integer counter
//   2) increment = step size between range max/mins/avgs
//   3) max_fringe = holds maximum values of subranges
//   4) min_fringe = holds minimum subrange values
//   5) avg_fringe = holds average subrange values
//   6) current_max_value  = holds current maximum value for loop
//   7) current_min_value  = holds current minimum value for loop
//   8) current_avg_value  = holds current average value for loop
// 
// ---------------------------------------------------------------

void calculate_fringe( void )
{
  int i;
  double increment, current_max_value,current_min_value, current_avg_value ;
  double max_fringe[15];
  double min_fringe[15];
  double avg_fringe[15];
  ofstream out(reinterpret_cast< const char *>(movie_name), ios::app );

  
  //     Initialize variables.  
  current_max_value = max_value;
  increment = ( max_value - min_value ) / ( num_colors - 1 );
  current_min_value = max_value - increment;
  current_avg_value = max_value - increment/2.0;
  for (i=0; i<15; i++ )
    {
      max_fringe[i] = 0.1;
      min_fringe[i] = 0.1;
      avg_fringe[i] = 0.1;
    }

  //     Loop to set range values.
  for ( i=0; i<(num_colors-1); i++ )
    {
      max_fringe[i] = current_max_value;
      min_fringe[i] = current_min_value;
      avg_fringe[i] = current_avg_value;
      current_max_value -= increment;
      current_min_value -= increment;
      current_avg_value -= increment;
    }

  //     Write to file.
  out << "sys_poll_option( 2 )" << endl << flush;
  out << "ga_spectrum_current_set( \"standard_spectrum\" )" << endl << flush;
  out << "ga_spectrum_delete( \"custom\" )" << endl << flush;
  out << "sys_poll_option( 0 )" << endl << flush;
  out << "sys_poll_option( 2 )" << endl << flush;
  out << "ga_viewport_range_set( \"\", \"standard_range\" )" << endl << flush;
  out << "ga_range_delete( \"custom\" )"  << endl << flush;
  out << "sys_poll_option( 0 )" << endl << flush;
  out << "ga_spectrum_create( \"custom\", " << num_colors << ", [0, 1, 8, 9, 3, 5, 14, 15, 10, 11, 2, 12, 4,  @" << endl << flush;
  out << "13, 6, 7] )" << endl << flush;
  out << "ga_range_create( \"custom\", " << num_colors-1 << " )" << endl << flush;
  out << "ga_range_values_set( \"custom\", [@" << endl << flush;
  //     Write maximums.
  for ( i=0; i<14; i++ )
    { 
      out << max_fringe[i] << ", @" << endl << flush;
    }
  out << max_fringe[14] << "], [@" << endl << flush;

  //     Write minimums.
    for ( i=0; i<14; i++ )
    { 
      out << min_fringe[i] << ", @" << endl << flush;
    }
  out << min_fringe[14] << "], [@" << endl << flush;

  //     Write averages.
    for ( i=0; i<14; i++ )
    { 
      out << avg_fringe[i] << ", @" << endl << flush;
    }
  out << avg_fringe[14] <<"] )" << endl << flush;
  out.close();
  return;
}
      
// ---------------------------------------------------------------
// DEFORMATION.C
//
// Purpose: This program writes instructions to create tif images
//          of deformed shape.
//
// Written by: Breanna Bailey
//
// Last modified: 8/31/01
//   
// Global Variables Used:
//   1) num_cycles = number of times loading pattern repeats
//   2) num_records = number of results files 
//   3) scale = "True" or "Model" for deformation scale
//   4) scale_factor = deformation scale factor
//   5) tif_name = user input name for tif files   
//   6) structure = name of structure analyzed
//   7) btwn_cycles = step size between successive cycles
//   8) records[] = vector of load steps to use
//  
// Local Variables Used:
//   1) file_1 = wnbd; extension for WARP3d results files
//   2) import_flag = indicates type of results to be imported
//   3) template_name = name of WARP3d results template
//   4) i, j, k = integer counters
//
// Functions Called:
//   1) import.c = separate file
// ---------------------------------------------------------------

void deformation( void )
{
  int i,j,k;
  char file_1[80];
  char import_flag[2];
  char template_name[80];
  ofstream out(reinterpret_cast< const char *>(movie_name), ios::app );

  //     Set prefix of results file name and import_flag.
  //     Call function to import results files.
  strcpy( file_1, "wnbd" );
  strcpy( import_flag, "D" );
  strcpy( template_name, "warp_displ.res_tmpl" );

  if (debug) cout << "In deformation..." << endl << flush;
  if ( import( file_1, import_flag, template_name ) )
    return;
  
  //     Select results to use as images.  
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )

	{
	i = records[j] + k*btwn_cycles;    
	out << "res_data_load_dbresult( 0, \"Nodal\", \"Vector\",  @" << endl << flush;
	out.fill(' ');
	out.setf(ios::left);

	out << "\"nodal displacement results for structure " << setw(8) << structure << ", step";
	out.setf(0,ios::left);
	out << setw(7) << i << "\", \"" << file_1;
	out.fill('0');
	out <<  setw(5) << i << "\", @" << endl << flush;
	out << "\"Displacements\", \"Translational\", \"(NON-LAYERED)\", \"\", \"AsIs\", \"\", \"\", \"\", \"\" @" << endl << flush;
	out << ", 0. )" << endl << flush;
	out << "res_data_title( 0, \"Nodal\", \"Vector\", 1, [ @" << endl << flush;
	out.fill(' ');	
	out.setf(ios::left);
	out << "\"nodal displacement results for structure " << setw(8) <<structure << ", step";
	out.setf(0,ios::left);
	out << setw(7) << i << ", " << file_1;
	out.fill('0');
	out << setw(5) << i << "\" // @" << endl << flush;
	out << "\", Displacements, Translational, (NON-LAYERED)\"] )" << endl << flush;	 
       
   	//     Create deformation.
	out << "res_display_deformation_create( \"\", \"Elements\", 0, [\"\"], 9, [ @" << endl << flush;
	out << "\"DeformedStyle:White,Solid,1,Wireframe\", \"DeformedScale:" << scale << scale_factor << "\",  @" << endl << flush;
	out << "\"UndeformedStyle:OFF,Blue,Solid,1,Wireframe\", \"TitleDisplay:OFF\",  @" << endl << flush;
	out << "\"MinMaxDisplay:OFF\", \"ScaleFactor:1.\", \"LabelStyle:Exponential, 12, White, 3\", @" << endl << flush;
	out << "\"DeformDisplay:Resultant\", \"DeformComps:OFF,OFF,OFF\"] )" << endl << flush;
	out << "res_display_deformation_post( \"\", 0 )" << endl << flush;
	out << "gm_write_image( \"TIFF\", \"" << tif_name << ".tif\", \"Increment\", 0., 0., 1., 1., 0 )" << endl << flush;
      }
  out.close();
  return;
}



// ------------------------------------------------------------
// EXTENSION.C
//
// Purpose: Creates extension of warp input file.
//          Called by import.c.
//
// Written by: Breanna Bailey
//
// Last modified: 9/3/01
//
//  Local IMPORT Variables Used:
//    1) step = current warp load step
//    2) extension = five digit extension to be found
//
//  Local Variables Used:
//    1) quotient = variable used to determine significant
//                  figures of extension
//    2) fatal_error = flag to indicate if current load step
//                     is too large                 
//
// -----------------------------------------------------------

int create_extension( int step, char extension[] )
{
  char *digits[10];
  int quotient;

  int fatal_error = false;


  digits[0]= "0";
  digits[1]= "1";
  digits[2]= "2";
  digits[3]= "3";
  digits[4]= "4";  
  digits[5]= "5";
  digits[6]= "6";
  digits[7]= "7"; 
  digits[8]= "8";
  digits[9]= "9";

  quotient = step/10000;


  if ( quotient >= 10 )
    {
      cout << endl << flush;
      cout << "ERROR: Results step is too large." << endl << flush;
      cout << "Program will terminate." << endl << flush;
      fatal_error = true;
    }
  else
    strcat( extension, digits[ quotient ] );  
  step -= quotient*10000;


  quotient = step / 1000;
  strcat( extension, digits[ quotient ] );
  step -= quotient*1000;
  
  quotient = step / 100;
  strcat( extension, digits[ quotient ] );
  step -= quotient*100;

  quotient = step / 10;
  strcat( extension, digits[ quotient ] );
  step -= quotient*10;

  quotient = step / 1;
  strcat( extension, digits[ quotient ] );
  step -= quotient*1;

  

  if ( fatal_error )
    return(1);
  else
    return(0);
}

// ------------------------------------------------------------
// GET_RECORDS.C
//
// Purpose: This program retrieves user input about the loading
//          pattern.  Called by user_input.c.  Three options
//          available: custom, incremental, all.
//
// Written by: Breanna Bailey
//
// Last modified: 9/3/01
//
// Global Variables Used:
//   1) num_cycles = number of times loading pattern repeats
//   2) num_records = number of results files 
//   3) btwn_cycles = step size between successive cycles
//   4) records[] = vector of load steps to use
//   5) results_format = nodal or element results
//   6) results_type = switch variable to select which type of
//                     image file to create (deformation,etc.)
//
//
// Local Variables Used:
//   1) i = integer counter
//   2) current record = counter used if incremental steps
//   3) start = if incremental, first load step
//   4) end = if incremental, last load step
//   5) by = if incremental, distance between load steps
//   6) pattern = flag to indicate which pattern selected
//   7) file1/file2 = file extensions "wXXX"
//   8) file_name_1/file_name_2 = websXXXXX; extensions with
//                                proper numbers
//   9) extension = XXXXX of above; numbers 
//  10) data_file_1/data_file_2 = results files
//
// -----------------------------------------------------------

int get_records( void )
{
  int i, current_record, start, end, by, pattern;
  char file1[200];
  char file2[200];
  char file_name_1[210];
  char file_name_2[210];
  char extension[6];  
  ifstream data_file_1;
  ifstream data_file_2;
  start = 0;
  end = 0;
  by = 0;

  //     Ask about load step pattern.
  cout << "Before continuing, check that the WARP3D results files you wish to use" << endl << flush;
  cout << "are already in the current directory." << endl << endl << flush;
  cout << "This utility has several methods of selecting which WARP3D results files" << endl << flush;
  cout << "will be used.  Which method would you prefer to use?" << endl << flush;
  cout << "1) Custom load step list (increments between load steps are not constant)."  << endl << flush;
  cout << "2) Constant load step list (increments between load steps are constant)." << endl << flush;
  cout << "3) Automatic load step list (program will select all WARP3D results files " << endl << flush;
  cout << "   of specified type in the current directory; not to exceed 10,000 files)." << endl << flush;
  cout << "Load Step List: ";
  cin  >> pattern;
  cout << endl << endl << flush;

  if (debug) cout << "Entering get_records..." << endl << flush;

  //     If automatic, open files.  Enter all found files into records.
  if ( pattern == 3 )
    {
      cout << "Please wait while results files are detected." << endl << flush;
      num_cycles = 1;
      btwn_cycles = 0;
      num_records = 0;
      switch( results_type )
	{
	  case 1:
			strcpy( file1, "wnbd" );
            strcpy( file2, "wnbd" );
            break;
      case 2:
	    if ( results_format == 1 )
	      { strcpy( file1, "wnbs" ); strcpy( file2, "wnbs" ); }
	    if ( results_format == 2 )
	      { strcpy( file1, "webs" ); strcpy( file2, "webs" ); }
	    break;
	  case 3:
	    if ( results_format == 1 )
	      { strcpy( file1, "wnbe" ); strcpy( file2, "wnbe" ); }
	    if ( results_format == 2 )
	      { strcpy( file1, "webe" ); strcpy( file2, "webe" ); } 
	    break;    
      case 4:
	    if ( results_format == 1 )
	      { strcpy( file1, "wnbs" ); strcpy( file2, "wnbd" ); }
	    if ( results_format == 2 )
	      { strcpy( file1, "webs" ); strcpy( file2, "wnbd" ); }
	    break;
      case 5:
   	    if ( results_format == 1 )
	      { strcpy( file1, "wnbe" ); strcpy( file2, "wnbd" ); }
	    if ( results_format == 2 )
	      { strcpy( file1, "webe" ); strcpy( file2, "wnbd" ); }
	}

      for ( i=1; i<=40; i++ )
	{
	  strcpy( file_name_1, file1 );
	  strcpy( file_name_2, file2 );
	  strcpy( extension, "" );
	  if ( create_extension( i, extension ) )
	  { 
		  cout << "An error occured while reading files.  Program will exit." << endl << flush;
		  return(1);
	  }
	  strcat( file_name_1, extension );
	  strcat( file_name_2, extension );
	  if (debug) cout << "File 1: " << file_name_1 << "   File 2: " << file_name_2 << endl << flush;
	  ifstream data_file_1(file_name_1, ios::nocreate);
      ifstream data_file_2(file_name_2, ios::nocreate);
	  if ( data_file_1.ipfx(0) && data_file_2.ipfx(0) )
	    {
	      records[num_records] = i;
		  if (debug) cout << "Filename found: " << file_name_1 << endl << flush;
	      num_records++;
	    }
	  data_file_1.close();
	  data_file_2.close();
	}
      cout << num_records << " records of the specified type were found in this directory." << endl << flush;
      cout << endl << endl << flush;
    }


  //     If constant, find number of records, set other variables.
  if ( pattern == 2 )
    {
      cout << "Please provide the beginning load step, the ending load step, and the increment" << endl << flush;
      cout << "between steps now.  (Enter only integer values.)" << endl << flush;
      cout << "Starting step:  ";
      cin >> start;  
      cout << "Final step:  ";
      cin >> end;
      cout << "Step increment:  ";
      cin >> by;
      cout << endl << endl << flush; 
      num_records = (end-start)/by + 1;
      btwn_cycles = 0;
      num_cycles = 1;
      current_record = start;
      for (i=0; i<num_records; i++ )
	  {
		  records[i]= current_record;
		  current_record += by;
	  }	
    }


  //     If custom, find out if cyclic, distance between cycles, etc.  Have user input.
  if ( pattern == 1 )
    {
      cout << "Custom load step lists require you to enter the load steps corresponding" << endl << flush;
      cout << "to the WARP3D results files to be displayed.  A cyclic ability has been" << endl << flush;
      cout << "added to simplify this process where possible.  For instance, load steps" << endl << flush;
      cout << "1, 3, 8, 10 might form the first load cycle, and steps 21, 23, 28, 30 might" << endl << flush;
      cout << "form the second.  You have the ability to enter only the first cycle and" << endl << flush;
      cout << "the program will calculate the second, and any subsequent, cycle(s)." << endl << endl << flush;
      cout << "How many cycles?" << endl << flush;
      cout << "Cycles: ";
      cin >> num_cycles;
      cout << endl << flush;
      cout << "How many files per cycle?" << endl << flush;
      cout << "Files Per Cycle: ";
      cin >> num_records;
      cout << endl << flush;
      cout << "How many load steps between cycles? (20 in above example)" << endl << flush;
      cout << "Cycle Increment: ";
      cin >> btwn_cycles;
      cout << endl << flush;
      cout << "Please enter the load steps in the first cycle now." << endl << flush;
      for (i = 0; i< num_records;i++ )
	  {
	  cout << "Record " << i+1 << ": ";
	  cin >> records[i];
	  }
      cout << endl << endl << flush;
  }
  return(0);
}


// ------------------------------------------------------------
// IMPORT.C
//
// Purpose: Imports results files into Patran.  Returns a
//          value of one if import fails.
//
// Written by: Breanna Bailey
//
// Last modified: 7/23/01
//
// Global Variables Used:
//   1) template_path = path to WARP3D templates
//   2) num_cycles = number of times loading 
//   3) num_records = number of results files 
//   4) btwn_cycles = step size between successive cycles
//   5) records[] = vector of load steps to use
//
// Local Variables Used:
//   1) file_name = location of file + file name; used to
//                  check if a results file exists
//   2) extension = five digit extension to follow webs, etc.
//   3) i, j, k  = integer counters
//   4) data_file = WARP3d results filename
//
// Functions Called:
//   1) extension.c = separate file
// -----------------------------------------------------------
int import( char file[], char import_flag[], char template_name[] )
{
  int i, j, k;
  char file_name[250];
  char extension[6];
  ifstream data_file;
  ofstream out(reinterpret_cast< const char *>(movie_name), ios::app  );

  out.fill('0');

  if (debug) cout << "Entering import..." << endl << flush;

  for ( k=0; k<num_cycles; k++ )
  {
	  for ( j=0; j<num_records; j++ )
      {
		  i = records[j] + k*btwn_cycles; 
		  if (debug)
		  {
			  cout << "Current i: " << i << endl << flush;	  
			  cout << "Current j: " << j << endl << flush;
			  cout << "Current k: " << k << endl << flush;
		  }

		  //     Generate character representation of file name. 
		  strcpy( file_name, file );
		  strcpy( extension, "" );
		  cout << extension << endl << flush;
		  if ( create_extension( i, extension ) ) {cout << "Error creating file extension." << endl << flush; return(1);}
		  strcat( file_name, extension );
		  if (debug) cout << "File name is: " << file_name << endl << flush;
		  //     Check that results file exists.  If not, exit.
		  data_file.open( file_name, ios::nocreate);
		  if (!data_file.ipfx(0))
		  {
			  cout << "ERROR: Results file " << file_name << " could not be found." << endl << flush;
			  cout << "Program will terminate." << endl << flush;
			  out << "PROGRAM TERMINATED BEFORE SESSION FILE COMPLETE." << endl << flush;
			  return(1);
		  }

		  //     Write commands to import results files. 
		  out << "resold_import_results( \"" << file;
		  out << setw(5) << i << "\",  @ " << endl << flush;
		  out << " \"" << import_flag << "\", 1E-06, \"";
		  out << template_path << "\"// @" << endl << flush;
		  out << "\"" << template_name << "\" ) " << endl << flush;
		  out << "uil_app_results2.set_update_display(  )"<< endl <<flush;
		  data_file.close();
	  }
  }
  out.close();
  return(0);
}


// ---------------------------------------------------------------
// MPEG_DEFORMATION.C
//
// Purpose: This program writes instructions to animate a
//          deformation file.  File is "previewed" or recorded.
//
// Written by: Breanna Bailey
//
// Last modified: 9/6/01
//   
// Global Variables Used:
//   1) num_cycles = number of times loading pattern repeats
//   2) num_records = number of results files 
//   3) scale = "True" or "Model" for deformation scale
//   4) scale_factor = deformation scale factor
//   5) structure = name of structure analyzed
//   6) btwn_cycles = step size between successive cycles
//   7) records[] = vector of load steps to use
//   8) frames = number of frames to record
//
// Local Variables Used:
//   1) file_1 = wnbd; extension for WARP3d results files
//   2) import_flag = indicates type of results to be imported
//   3) template_name = name of WARP3d results template
//   4) i, j, k = integer counters
//  
// Functions Called:
//   1) import.c = separate file
// ---------------------------------------------------------------
void mpeg_deformation( int record )
{
  int i,j,k;
  char file_1[80];  
  char import_flag[2];
  char template_name[80];
  ofstream out(reinterpret_cast< const char *>(movie_name), ios::app  );

  //     Set prefix of results file name and import_flag.
  //     Call function to import results files.
  strcpy( file_1, "wnbd" );
  strcpy( import_flag, "D" );
  strcpy( template_name, "warp_displ.res_tmpl" );
  if ( import( file_1, import_flag, template_name ) )
    return;

  //     Post results for first step.
	    out << "res_data_load_dbresult( 0, \"Nodal\", \"Vector\",  @" << endl << flush;
	    out.fill(' ');
	    out.setf(ios::left );
	    out << "\"nodal displacement results for structure " << setw(8) << structure << ", step";
	    out.setf( 0, ios:: left );
	    out << setw(7) << records[0] << "\", \"" << file_1;
	    out.fill('0');
	    out <<  setw(5) << records[0] << "\", @" << endl << flush;
	    out << "\"Displacements\", \"Translational\", \"(NON-LAYERED)\", \"\", \"AsIs\", \"\", \"\", \"\", \"\" @" << endl << flush;
	    out << ", 0. )" << endl << flush;
	    out << "res_data_title( 0, \"Nodal\", \"Vector\", 1, [ @" << endl << flush;
	    out.fill(' ');	
	    out.setf(ios::left);
	    out << "\"nodal displacement results for structure " << setw(8) <<structure << ", step";
	    out.setf(0,ios::left);
	    out << setw(7) <<  records[0] << ", " << file_1;
	    out.fill('0');
	    out << setw(5) <<  records[0] << "\" // @" << endl << flush;
	    out << "\", Displacements, Translational, (NON-LAYERED)\"] )" << endl << flush;
	    out << "res_data_dbres_list( 0, \"Nodal\", \"Vector\", " << num_records-1 << ", [ @" << endl << flush; 

  //     Write names of files.  
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles; 
	if ( j == 0 && k == 0 )
	  {
	    //  Don't do anything!
	  }
	else
	  {
        out.fill(' ');
	    out.setf(ios::left);
	    out << "\"nodal displacement results for structure " << setw(8) << structure <<", step";
	    out.setf( 0, ios::left );
	    out << setw(7) << i << "\"";
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << "], [ @" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }
      }


  //     Write filenames.  
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles; 
	if ( j == 0 && k == 0 )
	  {
	    //  Don't do anything!
	  }
	else
	  {
	    out.fill('0');
	    out << "\"" << file_1 << setw(5) << i << "\"";
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << "], [ @" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }
      }


  //     Write 'Displacements'.  
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles; 
	if ( j == 0 && k == 0 )
	  {
	    //  Don't do anything!
	  }
	else
	  {
	    out << "\"Displacements\"";
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << "], [ @" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }
      }


  //     Write 'Translational'.  
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles; 
	if ( j == 0 && k == 0 )
	  {
	    //  Don't do anything!
	  }
	else
	  {
	    out << "\"Translational\"";
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << "], [ @" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }
      }


  //     Write 'NON-LAYERED'.  
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles; 
	if ( j == 0 && k == 0 )
	  {
	    //  Don't do anything!
	  }
	else
	  {
	    out << "\"(NON-LAYERED)\"";
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << "] )" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }	      
      }



  //     Create preview/recording of deformation.
  out << "res_display_deformation_create( \"\", \"Elements\", 0, [\"\"], 9, [ @" << endl << flush;
  out << "\"DeformedStyle:White,Solid,1,Wireframe\", \"DeformedScale:"<< scale << scale_factor << "\",  @" << endl << flush;
  out << "\"UndeformedStyle:OFF,Blue,Solid,1,Wireframe\", \"TitleDisplay:OFF\",  @" << endl << flush;
  out << "\"MinMaxDisplay:OFF\", \"ScaleFactor:1.\", \"LabelStyle:Exponential, 12, White, 3\", @" << endl << flush;
  out << "\"DeformDisplay:Resultant\", \"DeformComps:OFF,OFF,OFF\"] )" << endl << flush;
  out << "res_display_deformation_post( \"\", 0 )" << endl << flush;
  out << "res_display_tool_animate_index( \"Deformation\", \"\", ";
  out << "1, " << num_records << " )" << endl << flush;

  //     Finish animation dependent on record variable.
  if ( record )
  {
    out << "res_display_anim_setup_3d( " << frames << ", \"Linear\" )" << endl << flush;
    out << "res_display_anim_save( \"MPEG\", \"Deformation\",";
    out << frames << ", \"Linear\" )"  << endl << flush;
  }
  else
  {
    out << "res_display_anim_setup_preview( ";
    out << frames << ", \"Linear\" )" << endl << flush;
  }
out.close();
return;
}



// ---------------------------------------------------------------
// MPEG_STRAIN.C
//
// Purpose: This program writes instructions to animate and
//          record strain fringes.
//
// Written by: Breanna Bailey
//
// Last modified: 9/6/01
//
// Global Variables Used:
//   1) num_cycles = number of times loading pattern repeats
//   2) num_records = number of results files 
//   3) structure = name of structure analyzed
//   4) btwn_cycles = step size between successive cycles
//   5) records[] = vector of load steps to use
//   6) frames = number of frames to record
//   7) results_format = nodal or element results
//   8) fringe_type = von mises, etc.
//   9) fringe_style = continuous or discrete
//  10) edge_options = flag to indicate if elements drawn
//
// Local Variables Used:
//   1) i, j, k = integer counters
//   2) style = character string for fringe_style
//   3) strain_results = character string for fringe_type
//   4) parameters = strings of Patran options that changes
//                   with nodal vs. element results
//   5) title = element or nodal
//   6) file_1 = wnbe or webe; extension for WARP3d results files
//   7) import_flag = indicates type of results to be imported
//   9) template_name = name of WARP3d results template
//
// Functions Called:
//   1) import.c = separate file
//   2) calculate_fringe.s = separate file
// ---------------------------------------------------------------
void mpeg_strain( void )
{
  int i,j, k;
  char style[25];
  char strain_results[150];
  char parameters[300];
  char title[30];
  char file_1[80];  
  char import_flag[2];
  char template_name[80];
  ofstream out(reinterpret_cast< const char *>(movie_name),  ios::app );

  //     Set user inputs.
  switch( results_format )
    { 
    case 1:
      strcpy( file_1, "wnbe" );
      strcpy( import_flag, "N" );  
      strcpy( template_name, "warp_node_strain.res_tmpl" ); 
      strcpy( parameters, "\"AsIs\", \"\", \"\", \"\", \"\"," );
      strcpy( title, "nodal" );
      break;
    case 2:
      strcpy( file_1, "webe" );
      strcpy( import_flag, "E" );  
      strcpy( template_name, "warp_elem_strain.res_tmpl" );
      strcpy( parameters, "\"AsIs\", \"DeriveAverage\", \"All\", \"ShapeFunc\", \"\"," );
     strcpy( title, "element" );
    }

  switch( fringe_type )
    {
    case 1:
      strcpy( strain_results, "eps-xx" );
      break;
    case 2:
      strcpy( strain_results, "eps-yy" );
      break;
    case 3:
      strcpy( strain_results, "eps-zz" );
      break;
    case 4:
      strcpy( strain_results, "gam-xy" );
      break;
    case 5:
      strcpy( stress_results, "gam-xz" );
      break;
    case 6:
      strcpy( strain_results, "gam-yz" );
      break;
    case 7:
      strcpy( strain_results, "Effective mises strain" );
    }

  switch( fringe_style )
    {
    case 1:
      strcpy( style, "Continuous" );
      break;
    case 2:
      strcpy( style, "Discrete/Smooth" );
    }
  
  if ( import( file_1, import_flag, template_name ) )
    return;
  calculate_fringe(); 

  //     Post first file.
  out << "res_data_load_dbresult( 0, \"Nodal\", \"Scalar\",  @" << endl << flush;
  out.fill(' ');
  out.setf(ios::left);
  out << "\"" << title << " strain results for structure "<< setw(8) << structure << ", step";
  out.setf( 0, ios:: left );
  out << setw(7) << records[0] <<"\", \"" << file_1;
  out.fill('0');
  out << setw(5) << records[0] << "\",  @" << endl << flush;
  out << "\"Strain\", \"" << strain_results << "\", \"(NON-LAYERED)\", \"\", @" << endl << flush;
  out << parameters << " 0. )"  << endl << flush;
  out << "res_data_title( 0, \"Nodal\", \"Scalar\", 1, [ @"  << endl << flush;
  out.fill(' ');
  out.setf(ios::left);
  out << "\"" << title << " strain results for structure " << setw(8) << structure << ", step";
  out.setf(0,ios::left);
  out << setw(7) << records[0] << ", " << file_1;
  out.fill('0');  
  out << setw(5) <<  records[0]<< ", St\" // @"  << endl << flush;
  out << "\"rain, " << strain_results << ", (NON-LAYERED)\"] )"  << endl << flush;
  out << "res_data_dbres_list( 0, \"Nodal\", \"Scalar\", " << num_records-1 << ", [ @" << endl << flush;

  //     Loop to write titles.
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles;
	//     Post first results.
	if ( j==0 && k==0 )
	  {
	    //     Don't do anything!
	  }
	else
	  {
	    out.fill(' ');
	    out.setf(ios::left);
	    out << "\"" << title << " strain results for structure " << setw(8) << structure << ", step";
	    out.setf(0,ios::left);
	    out << setw(7) << i << "\"";
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << "], [ @" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }
      }

  //     Loop to write file names.
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles;
	//     Post first results.
	if ( j==0 && k==0 )
	  {
	    //     Don't do anything!
	  }
	else
	  {
	    out << "\"" << file_1;
	    out.fill('0');
	    out << setw(5) << i << "\""; 
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << "], [ @" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }
      }

 //     Loop to write 'Strain'.
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles;
	//     Post first results.
	if ( j==0 && k==0 )
	  {
	    //     Don't do anything!
	  }
	else
	  {
	    out << "\"Strain\"";
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << "], [ @" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }
      }

 //     Loop to write strain result.
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles;
	//     Post first results.
	if ( j==0 && k==0 )
	  {
	    //     Don't do anything!
	  }
	else
	  {
	    out << "\"" << strain_results << "\"";
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << "], [ @" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }
      }


//     Loop to write 'NON-LAYERED'.
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles;
	//     Post first results.
	if ( j==0 && k==0 )
	  {
	    //     Don't do anything!
	  }
	else
	  {
	    out << "\"(NON-LAYERED)\"";
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << ", \"(NON-LAYERED)\"] )" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }
      }

  //    Create animation.
  out << "res_data_list_max( 0, \"Nodal\", \"Scalar\", \"Algebraic\", \"\" )" << endl << flush;
  out << "res_display_fringe_create( \"\", \"FreeFaces\", 0, [\"\"], 12, [ @" << endl << flush;
  out << "\"Range:custom\", \"RangeOverwrite:OFF\", \"FringeStyle:";
  out << style << "\",  @" << endl << flush;
  out << "\"Shade:None\", \"ElemEdge:" << edge_option << ",Blue,Solid,1\", \"Shrink:0\", \"TitleDisplay:OFF\" @ " << endl << flush;
  out << ", \"MinMaxDisplay:OFF\", \"ValueDisplay:OFF\", \"Filter:None\", \"ScaleFactor:1.\",  @" << endl << flush;
  out << "\"LabelStyle:Exponential, 12, White, 3\"], TRUE )" << endl << flush;
  out << "res_display_fringe_post( \"\", 0, \"Nodal\", FALSE, TRUE )" << endl << flush;
  out << "res_display_tool_animate_index( \"Fringe\", \"\", 1, ";
  out << num_records << " )" << endl << flush;
  out << "res_display_anim_setup_3d( " << frames << ", \"Linear\" )" << endl << flush;
  out << "res_display_anim_save( \"MPEG\", \"Fringe\",";
  out << frames << ", \"Linear\" )"  << endl << flush;
  out.close();
return;
}


// ---------------------------------------------------------------
// MPEG_STRESS.C
//
// Purpose: This program writes instructions to animate and
//          record stress fringes.
//
// Written by: Breanna Bailey
//
// Last modified: 9/6/01
//
// Global Variables Used:
//   1) num_cycles = number of times loading pattern repeats
//   2) num_records = number of results files 
//   3) structure = name of structure analyzed
//   4) btwn_cycles = step size between successive cycles
//   5) records[] = vector of load steps to use
//   6) frames = number of frames to record
//   7) results_format = nodal or element results
//   8) fringe_type = von mises, etc.
//   9) fringe_style = continuous or discrete
//  10) edge_options = flag to indicate if elements drawn
//
// Local Variables Used:
//   1) i, j, k = integer counters
//   2) style = character string for fringe_style
//   3) stress_results = character string for fringe_type
//   4) parameters = strings of Patran options that changes
//                   with nodal vs. element results
//   5) title = element or nodal
//   6) file_1 = wnbs or webs; extension for WARP3d results files
//   7) import_flag = indicates type of results to be imported
//   8) template_name = name of WARP3d results template
//
// Functions Called:
//   1) import.c = separate file
//   2) calculate_fringe.s = separate file
// ---------------------------------------------------------------
void mpeg_stress( void )
{
  int i, j, k;
  char style[25];
  char stress_results[150];
  char parameters[300];
  char title[30];
  char file_1[80];
  char import_flag[2];
  char template_name[80];
  ofstream out(reinterpret_cast< const char *>(movie_name), ios::app  );
  

  //     Set user inputs.
  switch( results_format )
    { 
    case 1:
      strcpy( file_1, "wnbs" );
      strcpy( import_flag, "N" );  
      strcpy( template_name, "warp_node_stress.res_tmpl" ); 
      strcpy( parameters, "\"AsIs\", \"\", \"\", \"\", \"\"," );
      strcpy( title, "nodal" );
      break;
    case 2:
      strcpy( file_1, "webs" );
      strcpy( import_flag, "E" );  
      strcpy( template_name, "warp_elem_stress.res_tmpl" );
      strcpy( parameters, "\"AsIs\", \"DeriveAverage\", \"All\", \"ShapeFunc\", \"\"," );
     strcpy( title, "element" );
    }

  switch( fringe_type )
    {
    case 1:
      strcpy( stress_results, "sigma-x" );
      break;
    case 2:
      strcpy( stress_results, "sigma-y" );
      break;
    case 3:
      strcpy( stress_results, "sigma-z" );
      break;
    case 4:
      strcpy( stress_results, "tau-xy" );
      break;
    case 5:
      strcpy( stress_results, "tau-xz" );
      break;
    case 6:
      strcpy( stress_results, "tau-yz" );
      break;
    case 7:
      strcpy( stress_results, "von mises stress" );
      break;
    case 8:
      strcpy( stress_results, "work density" );
      break;
    case 9:
      strcpy( stress_results, "C1-Material model dependent value (GT: matrix strain)" );
      break;
    case 10:
      strcpy( stress_results, "C2-Material model dependent value (GT: matrix stress)" );
      break;
    case 11:
      strcpy( stress_results, "C3-Material model dependent value (GT: f)" );
    }

  switch( fringe_style )
    {
    case 1:
      strcpy( style, "Continuous" );
      break;
    case 2:
      strcpy( style, "Discrete/Smooth" );
    }
  
  if ( import( file_1, import_flag, template_name ) )
    return;
  calculate_fringe(); 

  //     First record.
  out << "res_data_load_dbresult( 0, \"Nodal\", \"Scalar\",  @" << endl << flush;
  out.fill(' ');
  out.setf(ios::left);
  out << "\"" << title << " stress results for structure "<< setw(8) << structure << ", step";
  out.setf(0,ios::left);
  out << setw(7) << records[0] <<"\", \"" << file_1;
  out.fill('0');
  out << setw(5) << records[0] << "\",  @" << endl << flush;
  out << "\"Stress\", \"" << stress_results << "\", \"(NON-LAYERED)\", \"\", @" << endl << flush;
  out << parameters << " 0. )"  << endl << flush;
  out << "res_data_title( 0, \"Nodal\", \"Scalar\", 1, [ @"  << endl << flush;
  out.fill(' ');
  out.setf(ios::left);
  out << "\"" << title << " stress results for structure " << setw(8) << structure << ", step";
  out.setf(0,ios::left);
  out << setw(7) <<  records[0]<< ", " << file_1;
  out.fill('0');  
  out << setw(5) << records[0] << ", St\" // @"  << endl << flush;
  out << "\"ress, " << stress_results << ", (NON-LAYERED)\"] )"  << endl << flush;
  out << "res_data_dbres_list( 0, \"Nodal\", \"Scalar\", " << num_records-1 << ", [ @" << endl << flush;

  //     Loop to write titles.
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles;
	//     Post first results.
	if ( j==0 && k==0 )
	  {
	    //     Do nothing.
	  }
	else
	  {
	    out.fill(' ');
	    out.setf(ios::left);
	    out << "\"" << title << " stress results for structure " << setw(8) << structure << ", step";
	    out.setf(0,ios::left);
	    out << setw(7) << i << "\"";
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << "], [ @" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }
      }

  //     Loop to write filenames.
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles;
	//     Post first results.
	if ( j==0 && k==0 )
	  {
	    //     Do nothing.
	  }
	else
	  {
	    out << "\"" << file_1;
	    out.fill('0');
	    out << setw(5) << i << "\"";
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << "], [ @" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }
      }


  //     Loop to write 'Stress'.
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles;
	//     Post first results.
	if ( j==0 && k==0 )
	  {
	    //     Do nothing.
	  }
	else
	  {
	    out << "\"Stress\"";
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << "], [ @" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }
      }


  //     Loop to write stress_results.
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles;
	//     Post first results.
	if ( j==0 && k==0 )
	  {
	    //     Do nothing.
	  }
	else
	  {
	    out << "\"" << stress_results << "\"";
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << "], [ @" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }
      }

  //     Loop to write 'NON-LAYERED'.
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles;
	//     Post first results.
	if ( j==0 && k==0 )
	  {
	    //     Don't do anything!
	  }
	else
	  {
	    out << "\"(NON-LAYERED)\"";
	    if ( j==(num_records-1) && k==(num_cycles-1) )
	      out << ", \"(NON-LAYERED)\"] )" << endl << flush;
	    else
	      out << ", @" << endl << flush;
	  }
      }


  //    Create animation.
  out << "res_data_list_max( 0, \"Nodal\", \"Scalar\", \"Algebraic\", \"\" )" << endl << flush;
  out << "res_display_fringe_create( \"\", \"FreeFaces\", 0, [\"\"], 12, [ @" << endl << flush;
  out << "\"Range:custom\", \"RangeOverwrite:OFF\", \"FringeStyle:";
  out << style << "\",  @" << endl << flush;
  out << "\"Shade:None\", \"ElemEdge:" << edge_option << ",Blue,Solid,1\", \"Shrink:0\", \"TitleDisplay:OFF\" @ " << endl << flush;
  out << ", \"MinMaxDisplay:OFF\", \"ValueDisplay:OFF\", \"Filter:None\", \"ScaleFactor:1.\",  @" << endl << flush;
  out << "\"LabelStyle:Exponential, 12, White, 3\"], TRUE )" << endl << flush;
  out << "res_display_fringe_post( \"\", 0, \"Nodal\", FALSE, TRUE )" << endl << flush;
  out << "res_display_tool_animate_index( \"Fringe\", \"\", 1, ";
  out << num_records << " )" << endl << flush;
  out << "res_display_anim_setup_3d( " << frames << ", \"Linear\" )" << endl << flush;
  out << "res_display_anim_save( \"MPEG\", \"Fringe\",";
  out << frames << ", \"Linear\" )"  << endl << flush;
  out.close();
return;
}


// ---------------------------------------------------------------
// STRAIN.C
//
// Purpose: This program writes instructions to create tif files
//          of strain fringes.
//
// Written by: Breanna Bailey
//
// Last modified: 8/31/01
//
// Global Variables Used:
//   1) num_cycles = number of times loading pattern repeats
//   2) num_records = number of results files 
//   3) tif_name = user input name for tif files   
//   4) structure = name of structure analyzed
//   5) btwn_cycles = step size between successive cycles
//   6) records[] = vector of load steps to use
//   7) results_format = nodal or element results
//   8) fringe_type = von mises, etc.
//   9) fringe_style = continuous or discrete
//
// Local Variables Used:
//   1) i, j, k = integer counters
//   2) style = character string for fringe_style
//   3) strain_results = character string for fringe_type
//   4) parameters = strings of Patran options that changes
//                   with nodal vs. element results
//   5) title = element or nodal 
//   6) file_1 = wnbe or webe; extension for WARP3d results files
//   7) import_flag = indicates type of results to be imported
//   8) template_name = name of WARP3d results template
//
// Functions Called:
//   1) import.c = separate file
//   2) calculate_fringe.c = separate file
// ---------------------------------------------------------------
void strain( void )
{
  int i,j, k;
  char style[25];
  char strain_results[150];
  char parameters[300];
  char title[30];
  char file_1[80];  
  char import_flag[2];
  char template_name[80];
  ofstream out(reinterpret_cast< const char *>(movie_name), ios::app  );

  //     Set user inputs.
  switch( results_format )
    { 
    case 1:
      strcpy( file_1, "wnbe" );
      strcpy( import_flag, "N" );  
      strcpy( template_name, "warp_node_strain.res_tmpl" ); 
      strcpy( parameters, "\"AsIs\", \"\", \"\", \"\", \"\"," );
      strcpy( title, "nodal" );
      break;
    case 2:
      strcpy( file_1, "webe" );
      strcpy( import_flag, "E" );  
      strcpy( template_name, "warp_elem_strain.res_tmpl" );
      strcpy( parameters, "\"AsIs\", \"DeriveAverage\", \"All\", \"ShapeFunc\", \"\"," );
     strcpy( title, "element" );
    }

  switch( fringe_type )
    {
    case 1:
      strcpy( strain_results, "eps-xx" );
      break;
    case 2:
      strcpy( strain_results, "eps-yy" );
      break;
    case 3:
      strcpy( strain_results, "eps-zz" );
      break;
    case 4:
      strcpy( strain_results, "gam-xy" );
      break;
    case 5:
      strcpy( stress_results, "gam-xz" );
      break;
    case 6:
      strcpy( strain_results, "gam-yz" );
      break;
    case 7:
      strcpy( strain_results, "Effective mises strain" );
    }

  switch( fringe_style )
    {
    case 1:
      strcpy( style, "Continuous" );
      break;
    case 2:
      strcpy( style, "Discrete/Smooth" );
    }

  if ( import( file_1, import_flag, template_name ) )
    return;
    
   //     FIRST, create fringe for max step.  This sets range.
   //     Write commands to post stress results.
  calculate_fringe();

  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles; 
	//     Write commands to post stress results.
	out << "res_data_load_dbresult( 0, \"Nodal\", \"Scalar\",  @" << endl << flush;
	out.fill(' ');  
	out.setf(ios::left );
	out << "\"" << title << " strain results for structure" << setw(8) << structure << ", step";
	out.setf(0,ios::left );
	out << setw(7) << i <<"\", \"" << file_1;
	out.fill('0');
	out << setw(5) << i << "\",  @" << endl << flush;
	out << "\"Strain\", \"" << strain_results << "\", \"(NON-LAYERED)\", \"\", @" << endl << flush;
	out << parameters << " 0. )"  << endl << flush;
	out << "res_data_title( 0, \"Nodal\", \"Scalar\", 1, [ @"  << endl << flush;
	out.fill(' ');	
	out.setf(ios::left );
	out << "\"" << title << " strain results for structure" << setw(8) << structure << ", step";
	out.setf(0,ios::left );
	out << i << ", " << file_1;
	out.fill('0');
	out << setw(5) << i << ", St\" // @"  << endl << flush;
	out << "\"rain, " << strain_results << ", (NON-LAYERED)\"] )"  << endl << flush;
	
	//    Create fringe plot.
	out << "res_display_fringe_create( \"\", \"FreeFaces\", 0, [\"\"], 12, [ @" << endl << flush;
	out << "\"Range:custom\", \"RangeOverwrite:OFF\", \"FringeStyle:";
	out << style << "\",  @" << endl << flush;
	out << "\"Shade:None\", \"ElemEdge:" << edge_option << ",Blue,Solid,1\", \"Shrink:0\", \"TitleDisplay:OFF\" @ " << endl << flush;
	out << ", \"MinMaxDisplay:OFF\", \"ValueDisplay:OFF\", \"Filter:None\", \"ScaleFactor:1.\",  @" << endl << flush;
	out << "\"LabelStyle:Exponential, 12, White, 3\"], TRUE )" << endl << flush;
	out << "res_display_fringe_post( \"\", 0, \"Nodal\", TRUE, TRUE )" << endl << flush;
	
	//     Create TIF FILE.
	out << "gm_write_image( \"TIFF\", \"" << tif_name << ".tif\", \"Increment\", 0., 0., 1., 1., 0 )" << endl << flush; 
      }
  out.close();
  return;
}



// ---------------------------------------------------------------
// STRAIN_AND_DEFORMATION.C
//
// Purpose: This program creates a deformation and then posts a
//          strain fringe on top of the deformation.
//
// Written by: Breanna Bailey
//
// Last modified: 8/31/01
//   
// Global Variables Used:
//   1) num_cycles = number of times loading pattern repeats
//   2) num_records = number of results files 
//   3) tif_name = user input name for tif files   
//   4) structure = name of structure analyzed
//   5) btwn_cycles = step size between successive cycles
//   6) records[] = vector of load steps to use
//   7) results_format = nodal or element results
//   8) fringe_type = von mises, etc.
//   9) fringe_style = continuous or discrete
//  10) scale = "True" or "Model" for deformation scale
//  11) scale_factor = deformation scale factor
//
// Local Variables Used:
//   1) i, j, k = integer counters
//   2) style = character string for fringe_style
//   3) strain_results = character string for fringe_type
//   4) parameters = strings of Patran options that changes
//                   with nodal vs. element results
//   5) title = element or nodal
//   6) file_1 = file prefix for displacements
//   7) file_2 = file prefix for strains
//   8) import_flag = indicates type of results to be imported
//   9) template_name = name of WARP3d results template
//
// Functions Called:
//   1) import.c = separate file
//   2) calculate_fringe.c = separate file
// ---------------------------------------------------------------
void strain_and_deformation( void )
{
  int i, j, k;
  char style[25];
  char strain_results[150];
  char parameters[300];
  char title[30];
  char file_1[80];
  char file_2[80];  
  char import_flag[2];
  char template_name[80];
  ofstream out(reinterpret_cast< const char *>(movie_name), ios::app  );

  //     FOR DEFORMATION.
  //     Set prefix of results file name and import_flag.
  //     Call function to import results files.
  strcpy( file_1, "wnbd" );
  strcpy( import_flag, "D" );
  strcpy( template_name, "warp_displ.res_tmpl" );
  if ( import( file_1, import_flag, template_name ) )
    return;

  //     FOR STRAINS.
  //     Set variables according to user inputs.
  switch( results_format )
    { 
    case 1:
      strcpy( file_2, "wnbe" );
      strcpy( import_flag, "N" );  
      strcpy( template_name, "warp_node_strain.res_tmpl" ); 
      strcpy( parameters, "\"AsIs\", \"\", \"\", \"\", \"\"," );
      strcpy( title, "nodal" );
      break;
    case 2:
      strcpy( file_2, "webe" );
      strcpy( import_flag, "E" );  
      strcpy( template_name, "warp_elem_strain.res_tmpl" );
      strcpy( parameters, "\"AsIs\", \"DeriveAverage\", \"All\", \"ShapeFunc\", \"\"," );
     strcpy( title, "element" );
    }

  switch( fringe_type )
    {
    case 1:
      strcpy( strain_results, "eps-xx" );
      break;
    case 2:
      strcpy( strain_results, "eps-yy" );
      break;
    case 3:
      strcpy( strain_results, "eps-zz" );
      break;
    case 4:
      strcpy( strain_results, "gam-xy" );
      break;
    case 5:
      strcpy( stress_results, "gam-xz" );
      break;
    case 6:
      strcpy( strain_results, "gam-yz" );
      break;
    case 7:
      strcpy( strain_results, "Effective mises strain" );
    }

  switch( fringe_style )
    {
    case 1:
      strcpy( style, "Continuous" );
      break;
    case 2:
      strcpy( style, "Discrete/Smooth" );
    }

  if ( import( file_2, import_flag, template_name ) )
    return;

  calculate_fringe();
  //     Begin main loop to create deformation, then stress fringe
  //     for each time step.
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles; 
	
	//     DEFORMATION.
	//     Load displacement results.
	out << "res_data_load_dbresult( 0, \"Nodal\", \"Vector\",  @" << endl << flush;
	out.fill(' ');
	out.setf(ios::left );
	out << "\"nodal displacement results for structure " << setw(8) << structure << ", step";
	out.setf(0,ios::left);
	out << setw(7) << i << "\", \"" << file_1;
	out.fill('0');
	out <<  setw(5) << i << "\", @" << endl << flush;
	out << "\"Displacements\", \"Translational\", \"(NON-LAYERED)\", \"\", \"AsIs\", \"\", \"\", \"\", \"\" @" << endl << flush;
	out << ", 0. )" << endl << flush;
	out << "res_data_title( 0, \"Nodal\", \"Vector\", 1, [ @" << endl << flush;
	out.fill(' ');
	out.setf(ios::left );
	out << "\"nodal displacement results for structure " << setw(8) << structure << ", step";
	out.setf(0,ios::left );
	out << setw(7) << i << ", " << file_1;
	out.fill('0');
	out << setw(5) << i << "\" // @" << endl << flush;
	out << "\", Displacements, Translational, (NON-LAYERED)\"] )" << endl << flush;
	//     Create deformation.
	out << "res_display_deformation_create( \"\", \"Elements\", 0, [\"\"], 9, [ @" << endl << flush;
	out << "\"DeformedStyle:White,Solid,1,Wireframe\", \"DeformedScale:" << scale << scale_factor << "\",  @" << endl << flush;
	out << "\"UndeformedStyle:OFF,Blue,Solid,1,Wireframe\", \"TitleDisplay:OFF\",  @" << endl << flush;
	out << "\"MinMaxDisplay:OFF\", \"ScaleFactor:1.\", \"LabelStyle:Exponential, 12, White, 3\", @" << endl << flush;
	out << "\"DeformDisplay:Resultant\", \"DeformComps:OFF,OFF,OFF\"] )" << endl << flush;
	out << "res_display_deformation_post( \"\", 0 )" << endl << flush;

	//     STRAIN.
	//     Write commands to post stress results.
	out << "res_data_load_dbresult( 0, \"Nodal\", \"Scalar\",  @" << endl << flush;
	out.fill(' ');
	out.setf(ios::left );
	out << "\"" << title << " strain results for structure " << setw(8) << structure << ", step";
	out.setf(0,ios::left );
	out << setw(7) << i <<"\", \"" << file_2;
	out.fill('0');
	out << setw(5) << i << "\",  @" << endl << flush;
	out << "\"Strain\", \"" << strain_results << "\", \"(NON-LAYERED)\", \"\", @" << endl << flush;
	out << parameters << " 0. )"  << endl << flush;
	out << "res_data_title( 0, \"Nodal\", \"Scalar\", 1, [ @"  << endl << flush;
	out.fill(' ');	
	out.setf(ios::left );
	out << "\"" << title << " strain results for structure " << setw(8) << structure << ", step";
	out.setf(0,ios::left );
	out << i << ", " << file_2;
	out.fill('0');
	out << setw(5) << i << ", St\" // @"  << endl << flush;
	out << "\"rain, " << strain_results << ", (NON-LAYERED)\"] )"  << endl << flush;
	//    Create fringe plot.
	out << "res_display_fringe_create( \"\", \"FreeFaces\", 0, [\"\"], 12, [ @" << endl << flush;
	out << "\"Range:custom\", \"RangeOverwrite:OFF\", \"FringeStyle:";
	out << style << "\",  @" << endl << flush;
	out << "\"Shade:None\", \"ElemEdge:" << edge_option << ",Blue,Solid,1\", \"Shrink:0\", \"TitleDisplay:OFF\" @ " << endl << flush;
	out << ", \"MinMaxDisplay:OFF\", \"ValueDisplay:OFF\", \"Filter:None\", \"ScaleFactor:1.\",  @" << endl << flush;
	out << "\"LabelStyle:Exponential, 12, White, 3\"], TRUE )" << endl << flush;
	out << "res_display_fringe_post( \"\", 0, \"Nodal\", TRUE, TRUE )" << endl << flush;
	
	//     Create TIF FILE.
	out << "gm_write_image( \"TIFF\", \"" << tif_name << ".tif\", \"Increment\", 0., 0., 1., 1., 0 )" << endl << flush;  

      }
  out.close();
  return;
}



// ---------------------------------------------------------------
// STRAIN_OPTIONS.C
//
// Purpose: To input user options for stress fringes.
//
// Written by: Breanna Bailey
//
// Last modified: 9/4/01
//
// Global Variables Used:
//   1) results_type = stress, strain, or deformation, or combo
//   2) fringe_type = what component of stress, strain, to plot
//   3) template_path = location or warp templates
//   5) fringe_style = continuous or discrete/smooth
//   6) results_location = nodal or element stresses/strains
//   7) num_records = number of results files per cycle
//   8) num_cycles = the number of results cycles
//   9) btwn_cycles = increment between results cycles
//  10) records[i] = list of records to be processed (per cycle)
//  11) scale = model or true scale
//  12) scale_factor = real number used to scale deformation
//  13) format = flag to indicate mpeg or tif output
//  14) edge_option = character string based on edges
//
// Local Variables Used:
//   1) pattern = constant or nonconstant increment option
//   2) deformation_scale = model or true scale
//   3) edges = flag to indicate edge display 
//
//   NO ERROR CHECKING
// ---------------------------------------------------------------
void strain_options( void )
{    
  int edges;

 cout << "Strain Fringe Options: " << endl << endl << flush;
 cout << "What strain component do you wish to display?" << endl << flush;
 cout << "Please select an option from the following menu: " << endl << flush;
 cout << "   1) Epsilon XX." << endl << flush;
 cout << "   2) Epsilon YY." << endl << flush;
 cout << "   3) Epsilon ZZ." << endl << flush;
 cout << "   4) Gamma XY." << endl << flush;
 cout << "   5) Gamma XZ." << endl << flush;
 cout << "   6) Gamma YZ." << endl << flush;
 cout << "   7) Effective Mises Strain." << endl << flush;
 cout << "Option:  ";
 cin >> fringe_type;
 cout << endl << endl << flush;
 cout << "What fringe style do you prefer?" << endl << flush;
 cout << "   1) Continuous." << endl << flush;
 cout << "   2) Discrete/Smooth." << endl << flush;
 cout << "Option:  ";
 cin >> fringe_style;
 cout << endl << endl << flush;
 cout << "Select an edge display option for the strain fringe." << endl << flush;
 cout << "   1) Element Edges." << endl << flush;
 cout << "   2) Free Edges." << endl << flush;
 cout << "   3) No edges." << endl << flush;
 cout << "Edge Display: ";
 cin >> edges;
 cout << endl << endl << flush;
 cout << "Will you be using nodal or element results?" << endl << flush;
 cout << "   1) Nodal." << endl << flush;
 cout << "   2) Element." << endl << flush;
 cout << "Option:  ";
 cin >> results_format;
 cout << endl << flush;
 cout << "How many colors do you wish to appear in the spectrum?" << endl << flush;
 cout << "(No more than sixteen.  No less than three.)" << endl << flush;
 cout << "Number of Colors: ";
 cin >> num_colors;
 cout << endl << flush;
 cout << "To create an appropriate range, please enter the maximum" << endl << flush;
 cout << "and minumum strain values now." << endl << flush;
 cout << "Maximum Strain: ";
 cin >> max_value;
 cout << "Minumum Strain: ";
 cin >> min_value;
 cout << endl  << endl << flush;
 
  //  Set edge_option according to edges input.
  switch (edges)
    {
    case 1:
      strcpy( edge_option, "ElemEdge" );
      break;
    case 2:
      strcpy( edge_option, "FreeEdge" );
      break;
    case 3:
      strcpy( edge_option, "None" );
    } 
}



// ---------------------------------------------------------------
// STRESS.C
//
// Purpose: This program writes instructions to create tif images
//          of stress fringes.
//
// Written by: Breanna Bailey
//
// Last modified: 8/31/01
//
// Global Variables Used:
//   1) num_cycles = number of times loading pattern repeats
//   2) num_records = number of results files 
//   3) tif_name = user input name for tif files   
//   4) structure = name of structure analyzed
//   5) btwn_cycles = step size between successive cycles
//   6) records[] = vector of load steps to use
//   7) results_format = nodal or element results
//   8) fringe_type = von mises, etc.
//   9) fringe_style = continuous or discrete
//
// Local Variables Used:
//   1) i, j, k = integer counters
//   2) style = character string for fringe_style
//   3) stress_results = character string for fringe_type
//   4) parameters = strings of Patran options that changes
//                   with nodal vs. element results
//   5) title = element or nodal
//   6) file_1 = prefix webs or wnbs
//   7) import_flag = indicates type of results to be imported
//   8) template_name = name of WARP3d results template
//
// Functions Called:
//   1) import.c = separate file
//   2) calculate_fringe.c = separate file
// ---------------------------------------------------------------
void stress( void )
{
  int i, j, k;
  char style[25];
  char stress_results[150];
  char parameters[300];
  char title[30];
  char file_1[80]; 
  char import_flag[2];
  char template_name[80];
  ofstream out(reinterpret_cast< const char *>(movie_name),ios::app   );

  //     Set user inputs.
  switch( results_format )
    { 
    case 1:
      strcpy( file_1, "wnbs" );
      strcpy( import_flag, "N" );  
      strcpy( template_name, "warp_node_stress.res_tmpl" ); 
      strcpy( parameters, "\"AsIs\", \"\", \"\", \"\", \"\"," );
      strcpy( title, "nodal" );
      break;
    case 2:
      strcpy( file_1, "webs" );
      strcpy( import_flag, "E" );  
      strcpy( template_name, "warp_elem_stress.res_tmpl" );
      strcpy( parameters, "\"AsIs\", \"DeriveAverage\", \"All\", \"ShapeFunc\", \"\"," );
     strcpy( title, "element" );
    }

  switch( fringe_type )
    {
    case 1:
      strcpy( stress_results, "sigma-x" );
      break;
    case 2:
      strcpy( stress_results, "sigma-y" );
      break;
    case 3:
      strcpy( stress_results, "sigma-z" );
      break;
    case 4:
      strcpy( stress_results, "tau-xy" );
      break;
    case 5:
      strcpy( stress_results, "tau-xz" );
      break;
    case 6:
      strcpy( stress_results, "tau-yz" );
      break;
    case 7:
      strcpy( stress_results, "von mises stress" );
      break;
    case 8:
      strcpy( stress_results, "work density" );
      break;
    case 9:
      strcpy( stress_results, "C1-Material model dependent value (GT: matrix strain)" );
      break;
    case 10:
      strcpy( stress_results, "C2-Material model dependent value (GT: matrix stress)" );
      break;
    case 11:
      strcpy( stress_results, "C3-Material model dependent value (GT: f)" );
    }

  switch( fringe_style )
    {
    case 1:
      strcpy( style, "Continuous" );
      break;
    case 2:
      strcpy( style, "Discrete/Smooth" );
    }


  if ( import( file_1, import_flag, template_name ) )  
    return;
  calculate_fringe(); 
 

  //     Start loop to either list or post additional results files.
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
      i = records[j] + k*btwn_cycles;
      //     Write commands to post or list stress results.  
      out << "res_data_load_dbresult( 0, \"Nodal\", \"Scalar\",  @" << endl << flush;
      out.fill(' ');
      out.setf(ios::left);
      out << "\"" << title << " stress results for structure "<< setw(8) << structure << ", step";
      out.setf(0,ios::left);
      out << setw(7) << i <<"\", \"" << file_1;
      out.fill('0');
      out << setw(5) << i << "\",  @" << endl << flush;
      out << "\"Stress\", \"" << stress_results << "\", \"(NON-LAYERED)\", \"\", @" << endl << flush;
      out << parameters << " 0. )"  << endl << flush;
      out << "res_data_title( 0, \"Nodal\", \"Scalar\", 1, [ @"  << endl << flush;
      out.fill(' ');
      out.setf(ios::left);
      out << "\"" << title << " stress results for structure " << setw(8) << structure << ", step";
      out.setf(0,ios::left);
      out << setw(7) << i << ", " << file_1;
      out.fill('0');  
      out << setw(5) << i << ", St\" // @"  << endl << flush;
      out << "\"ress, " << stress_results << ", (NON-LAYERED)\"] )"  << endl << flush;
      //     Create fringe and record image.
      out << "res_display_fringe_create( \"\", \"FreeFaces\", 0, [\"\"], 12, [ @" << endl << flush;
      out << "\"Range:custom\", \"RangeOverwrite:OFF\", \"FringeStyle:";
      out << style << "\",  @" << endl << flush;
      out << "\"Shade:None\", \"ElemEdge:" << edge_option << ",Blue,Solid,1\", \"Shrink:0\", \"TitleDisplay:OFF\" @ " << endl << flush;
      out << ", \"MinMaxDisplay:OFF\", \"ValueDisplay:OFF\", \"Filter:None\", \"ScaleFactor:1.\",  @" << endl << flush;
      out << "\"LabelStyle:Exponential, 12, White, 3\"], TRUE )" << endl << flush;
      out << "res_display_fringe_post( \"\", 0, \"Nodal\", TRUE, TRUE )" << endl << flush;
      out << "gm_write_image( \"TIFF\", \"" << tif_name << ".tif\", \"Increment\", 0., 0., 1., 1., 0 )" << endl << flush;
      }
  out.close();
  return;
}



// ---------------------------------------------------------------
// STRESS_AND_DEFORMATION.C
//
// Purpose: This program creates a deformation and then posts a
//          stress fringe on top of the deformation.
//
// Written by: Breanna Bailey
//
// Last modified: 8/31/01
//   
// Global Variables Used:
//   1) num_cycles = number of times loading pattern repeats
//   2) num_records = number of results files 
//   3) tif_name = user input name for tif files   
//   4) structure = name of structure analyzed
//   5) btwn_cycles = step size between successive cycles
//   6) records[] = vector of load steps to use
//   7) results_format = nodal or element results
//   8) fringe_type = von mises, etc.
//   9) fringe_style = continuous or discrete
//  10) scale = "True" or "Model" for deformation scale
//  11) scale_factor = deformation scale factor
//
// Local Variables Used:
//   1) i = integer counter
//   2) style = character string for fringe_style
//   3) stress_results = character string for fringe_type
//   4) parameters = strings of Patran options that changes
//                  with nodal vs. element results
//   5) title = element or nodal
//   6) file_1 = file prefix for displacements
//   7) file_2 = file prefix for stresses
//   8) import_flag = indicates type of results to be imported
//   9) template_name = name of WARP3d results template
//
// Functions Called:
//   1) import.c = separate file
//   2) calculate_fringe.c = separate file
// ---------------------------------------------------------------
void stress_and_deformation( void )
{
  int i, j, k;
  char style[25];
  char stress_results[150];
  char parameters[300];
  char title[30];
  char file_1[80];
  char file_2[80];
  char import_flag[2];
  char template_name[80];
  ofstream out(reinterpret_cast< const char *>(movie_name), ios::app );

  //     FOR DEFORMATION.
  //     Set prefix of results file name and import_flag.
  //     Call function to import results files.
  strcpy( file_1, "wnbd" );
  strcpy( import_flag, "D" );
  strcpy( template_name, "warp_displ.res_tmpl" );
  if ( import( file_1, import_flag, template_name ) )
    return;



  //     FOR STRESSES.  
  //     Set variables according to user inputs.
  switch( results_format )
    { 
    case 1:
      strcpy( file_2, "wnbs" );
      strcpy( import_flag, "N" );  
      strcpy( template_name, "warp_node_stress.res_tmpl" ); 
      strcpy( parameters, "\"AsIs\", \"\", \"\", \"\", \"\"," );
      strcpy( title, "nodal" );
      break;
    case 2:
      strcpy( file_2, "webs" );
      strcpy( import_flag, "E" );  
      strcpy( template_name, "warp_elem_stress.res_tmpl" );
      strcpy( parameters, "\"AsIs\", \"DeriveAverage\", \"All\", \"ShapeFunc\", \"\"," );
     strcpy( title, "element" );
    }

  switch( fringe_type )
    {
    case 1:
      strcpy( stress_results, "sigma-x" );
      break;
    case 2:
      strcpy( stress_results, "sigma-y" );
      break;
    case 3:
      strcpy( stress_results, "sigma-z" );
      break;
    case 4:
      strcpy( stress_results, "tau-xy" );
      break;
    case 5:
      strcpy( stress_results, "tau-xz" );
      break;
    case 6:
      strcpy( stress_results, "tau-yz" );
      break;
    case 7:
      strcpy( stress_results, "von mises stress" );
      break;
    case 8:
      strcpy( stress_results, "work density" );
      break;
    case 9:
      strcpy( stress_results, "C1-Material model dependent value (GT: matrix strain)" );
      break;
    case 10:
      strcpy( stress_results, "C2-Material model dependent value (GT: matrix stress)" );
      break;
    case 11:
      strcpy( stress_results, "C3-Material model dependent value (GT: f)" );
    }

  switch( fringe_style )
    {
    case 1:
      strcpy( style, "Continuous" );
      break;
    case 2:
      strcpy( style, "Discrete/Smooth" );
    }

  if ( import( file_2, import_flag, template_name ) )
    return;   


  //     Begin main loop to create deformation, then stress fringe
  //     for each time step.
  calculate_fringe();
  
  for ( k=0; k<num_cycles; k++ )
    for ( j=0; j<num_records; j++ )
      {
	i = records[j] + k*btwn_cycles;

	//     DEFORMATION.
	//     Load displacement results.
	out << "res_data_load_dbresult( 0, \"Nodal\", \"Vector\",  @" << endl << flush;
	out.fill(' ');
	out.setf(ios::left );
	out << "\"nodal displacement results for structure " << setw(8) << structure << ", step";
	out.setf(0,ios::left );
	out << setw(7) << i << "\", \"" << file_1;
	out.fill('0');
	out <<  setw(5) << i << "\", @" << endl << flush;
	out << "\"Displacements\", \"Translational\", \"(NON-LAYERED)\", \"\", \"AsIs\", \"\", \"\", \"\", \"\" @" << endl << flush;
	out << ", 0. )" << endl << flush;
	out << "res_data_title( 0, \"Nodal\", \"Vector\", 1, [ @" << endl << flush;
	out.fill(' ');
	out.setf(ios::left );
	out << "\"nodal displacement results for structure " << setw(8) << structure << ", step";
	out.setf(0,ios::left );
	out << setw(7) << i << ", " << file_1;
	out.fill('0');
	out << setw(5) << i << "\" // @" << endl << flush;
	out << "\", Displacements, Translational, (NON-LAYERED)\"] )" << endl << flush;
	//     Create deformation.
	out << "res_display_deformation_create( \"\", \"Elements\", 0, [\"\"], 9, [ @" << endl << flush;
	out << "\"DeformedStyle:White,Solid,1,Wireframe\", \"DeformedScale:" << scale << scale_factor << "\",  @" << endl << flush;
	out << "\"UndeformedStyle:OFF,Blue,Solid,1,Wireframe\", \"TitleDisplay:OFF\",  @" << endl << flush;
	out << "\"MinMaxDisplay:OFF\", \"ScaleFactor:1.\", \"LabelStyle:Exponential, 12, White, 3\", @" << endl << flush;
	out << "\"DeformDisplay:Resultant\", \"DeformComps:OFF,OFF,OFF\"] )" << endl << flush;
	out << "res_display_deformation_post( \"\", 0 )" << endl << flush;
	
	//     STRESS.
	//     Load stress results.
	out << "res_data_load_dbresult( 0, \"Nodal\", \"Scalar\",  @" << endl << flush;
	out.fill(' ');
	out.setf(ios::left );
	out << "\"" << title << " stress results for structure " << setw(8) << structure << ", step";
	out.setf(0,ios::left );
	out << setw(7) << i <<"\", \"" << file_2;
	out.fill('0');
	out << setw(5) << i << "\",  @" << endl << flush;
	out << "\"Stress\", \"" << stress_results << "\", \"(NON-LAYERED)\", \"\", @" << endl << flush;
	out << parameters << " 0. )"  << endl << flush;
	out << "res_data_title( 0, \"Nodal\", \"Scalar\", 1, [ @"  << endl << flush;
	out.fill(' ');
	out.setf(ios::left );
	out << "\"" << title << " stress results for structure " << setw(8) << structure << ", step";
	out.setf(0,ios::left );
	out << setw(7) << i << ", " << file_2;
	out.fill('0');  
	out << setw(5) << i << ", St\" // @"  << endl << flush;
	out << "\"ress, " << stress_results << ", (NON-LAYERED)\"] )"  << endl << flush;
	//    Create fringe plot.
	out << "res_display_fringe_create( \"\", \"FreeFaces\", 0, [\"\"], 12, [ @" << endl << flush;
	out << "\"Range:custom\", \"RangeOverwrite:OFF\", \"FringeStyle:";
	out << style << "\",  @" << endl << flush;
	out << "\"Shade:None\", \"ElemEdge:"<< edge_option << ",Blue,Solid,1\", \"Shrink:0\", \"TitleDisplay:OFF\" @ " << endl << flush;
	out << ", \"MinMaxDisplay:OFF\", \"ValueDisplay:OFF\", \"Filter:None\", \"ScaleFactor:1.\",  @" << endl << flush;
	out << "\"LabelStyle:Exponential, 12, White, 3\"], TRUE )" << endl << flush;
      out << "res_display_fringe_post( \"\", 0, \"Nodal\", TRUE, TRUE )" << endl << flush;
      
      //     Create TIF FILE.
      out << "gm_write_image( \"TIFF\", \"" << tif_name << ".tif\", \"Increment\", 0., 0., 1., 1., 0 )" << endl << flush;  
      
      }
  out.close();
  return;
}


// ---------------------------------------------------------------
// STRESS_OPTIONS.C
//
// Purpose: To input user options for stress fringes.
//
// Written by: Breanna Bailey
//
// Last modified: 9/4/01
//
// Global Variables Used:
//   1) results_type = stress, strain, or deformation, or combo
//   2) fringe_type = what component of stress, strain, to plot
//   3) template_path = location or warp templates
//   5) fringe_style = continuous or discrete/smooth
//   6) results_location = nodal or element stresses/strains
//   7) num_records = number of results files per cycle
//   8) num_cycles = the number of results cycles
//   9) btwn_cycles = increment between results cycles
//  10) records[i] = list of records to be processed (per cycle)
//  11) scale = model or true scale
//  12) scale_factor = real number used to scale deformation
//  13) format = flag to indicate mpeg or tif output
//  14) edge_option = character string based on edges
//
// Local Variables Used:
//   1) pattern = constant or nonconstant increment option
//   2) deformation_scale = model or true scale
//   3) edges = flag to indicate edge display 
//
//   NO ERROR CHECKING
// ---------------------------------------------------------------
void stress_options( void )
{    
  int edges;


  cout << "Stress Fringe Options: " << endl << endl << flush;
  cout << "What stress component do you wish to display?" << endl << flush;
  cout << "Please select an option from the following menu: " << endl << flush;
  cout << "   1) Sigma X." << endl << flush;
  cout << "   2) Sigma Y." << endl << flush;
  cout << "   3) Sigma Z." << endl << flush;
  cout << "   4) Tau XY." << endl << flush;
  cout << "   5) Tau XZ." << endl << flush;
  cout << "   6) Tau YZ." << endl << flush;
  cout << "   7) Von Mises Stress." << endl << flush;
  cout << "   8) Work Density." << endl << flush;
  cout << "   9) C1-Material model dependent value." << endl << flush;
  cout << "  10) C2-Material model dependent value." << endl << flush;
  cout << "  11) C3-Material model dependent value." << endl << flush;
  cout << "Option:  ";
  cin >> fringe_type;
  cout << endl << endl << flush;
  cout << "What fringe style do you prefer?" << endl << flush;
  cout << "   1) Continuous." << endl << flush;
  cout << "   2) Discrete/Smooth." << endl << flush;
  cout << "Option:  ";
  cin >> fringe_style;      
  cout << endl << endl << flush;
  cout << "Select an edge display option for the stress fringe." << endl << flush;
  cout << "   1) Element Edges." << endl << flush;
  cout << "   2) Free Edges." << endl << flush;
  cout << "   3) No edges." << endl << flush;
  cout << "Edge Display: ";
  cin >> edges;
  cout << endl << endl << flush;
  cout << "Will you be using nodal or element results?" << endl << flush;
  cout << "   1) Nodal." << endl << flush;
  cout << "   2) Element." << endl << flush;
  cout << "Option:  ";
  cin >> results_format;
  cout << endl << endl << flush;
  cout << "How many colors do you wish to appear in the spectrum?" << endl << flush;
  cout << "(No more than sixteen.  No less than three.)" << endl << flush;
  cout << "Number of Colors: ";
  cin >> num_colors;
  cout << endl << flush;
  cout << "To create an appropriate range, please enter the maximum" << endl << flush;
  cout << "and minumum stress values now." << endl << flush;
  cout << "Maximum Stress: ";
  cin >> max_value;
  cout << "Minumum Stress: ";
  cin >> min_value;
  cout << endl << endl << flush;

  //  Set edge_option according to edges input.
  switch (edges)
    {
    case 1:
      strcpy( edge_option, "ElemEdge" );
      break;
    case 2:
      strcpy( edge_option, "FreeEdge" );
      break;
    case 3:
      strcpy( edge_option, "None" );
    }

    
}



// ---------------------------------------------------------------
// STRESS_OPTIONS.C
//
// Purpose: To input user options for stress fringes.
//
// Written by: Breanna Bailey
//
// Last modified: 9/4/01
//
// Global Variables Used:
//   1) results_type = stress, strain, or deformation, or combo
//   2) fringe_type = what component of stress, strain, to plot
//   3) template_path = location or warp templates
//   5) fringe_style = continuous or discrete/smooth
//   6) results_location = nodal or element stresses/strains
//   7) num_records = number of results files per cycle
//   8) num_cycles = the number of results cycles
//   9) btwn_cycles = increment between results cycles
//  10) records[i] = list of records to be processed (per cycle)
//  11) scale = model or true scale
//  12) scale_factor = real number used to scale deformation
//  13) format = flag to indicate mpeg or tif output
//  14) edge_option = character string based on edges
//
// Local Variables Used:
//   1) pattern = constant or nonconstant increment option
//   2) deformation_scale = model or true scale
//   3) edges = flag to indicate edge display 
//
//   NO ERROR CHECKING
// ---------------------------------------------------------------
void formatting( void );
void deformation_options(void);

int input( void )
{  

  char structure_temp[25];

  //     WELCOME Comments.
  cout << endl << endl << endl << flush;
  cout << "Welcome to the WARP3D To PATRAN Animation Utility." << endl << endl << flush;
  cout << "This program creates a PATRAN session file.  When played, this session" << endl  << flush;
  cout << "file drives the execution of PATRAN to read WARP3D results files and to create" << endl  << flush;
  cout << "either a series of image (.tif) files or an animated (.mpg) movie.  If" << endl << flush;
  cout << "desired, TIFF files can be subsequently converted to an animated (.gif) movie" << endl << flush;
  cout << "using the C-shell program, tif_to_gif, which accompanies this utility program." << endl << endl << flush;

  //     Display formatting options.
  formatting();

  //     Display results options.  
  cout << "What type of results are you post-processing?" << endl << flush;
  cout << "Please select an option from the following menu: " << endl << flush;
  cout << "   1) Deformed shapes." << endl << flush;
  cout << "   2) Stress fringes." << endl << flush;
  cout << "   3) Strain fringes." << endl << flush;
  cout << "   4) Stress fringes on deformed shapes." << endl << flush;
  cout << "   5) Strain fringes on deformed shapes." << endl << flush;
  cout << "Option:  ";
  cin >> results_type;
  cout << endl << endl << flush;

  //     Display deformation scale options.
  if ( results_type == 1 || results_type == 4 || results_type == 5 )
    deformation_options();

  //     Display stress fringe options.
  if ( results_type ==2 || results_type == 4 )
    stress_options();

  // Display strain fringe options.
  if ( results_type ==3 || results_type == 5 )
    strain_options();

  //     Determine file pathway.
  cout << "Please enter the location of the WARP3D template files." << endl << flush;
  cout << "The WARP3D template files are most often located in your" << endl << flush;
  cout << "/<PATRAN>/res_templates/ directory.  (Do not include template name.)" << endl << flush;
  cout << "Template Location:  ";
  cin >> template_path;
  cout << endl << endl << flush;

  //     Ask for structure's name.
  cout << "Please enter the name of the structure as defined in the WARP3D" << endl << flush;
  cout << "input file.  (No more than eight characters.)" << endl << flush;
  cout << "Structure: ";
  cin >> structure_temp;
  strncpy( structure, structure_temp, 8 );
  cout << endl << endl << flush;

  //     Loading pattern.
  if ( get_records( ) )
     return(1);
  
  //     Name file.
  cout << "Please enter a name for the PATRAN session file to be created." << endl << flush;
  cout << "(The extension *.ses.01 will be added automatically.)" << endl << flush;
  cout << "Session Name: ";
  cin >> movie_name;
  strcat( movie_name, ".ses.01" );
  cout << endl << endl << flush;

  //  Ending.
  cout << "This completes the user inputs.  Please wait a few moments" << endl;
  cout << "while your session file is generated." << endl << flush;
  
  return(0);
}




// ---------------------------------------------------------------
// FORMATTING FUNCTION
//
// Purpose: Reads formatting options.
//
// Written by: Breanna Bailey
//
// Last modified: 9/4/01
//
// Global Variables Used:
//   1) format = flag to indicate mpeg or tif output
//   2) frames = if mpeg format chosen, number of frames to 
//               include in movie
//   3) tif_name = if tif format chosen, base name for tif images
//
// ---------------------------------------------------------------

void formatting(void)
{
  //     MPEG OR TIF?
  cout << "Do you wish to create an MPEG movie or a series of TIFF image files?" << endl << flush;
  cout << "Please select an output format from the options below." << endl << flush;
  cout << "   1) A series of .tif image files. (One image file per result file.)" << endl << flush;
  cout << "   2) An animated .mpg movie." << endl << flush;
  cout << "Output Format: ";
  cin >> format;
  cout << endl << endl << flush;

  if ( format == 2 )
    {
      cout << "You have chosen to create a .mpg movie." << endl << endl << flush;  
      cout << "WARP3D results files become key frames of the movie.  If you" << endl << flush;
      cout << "select a number of frames higher than the number of results" << endl << flush;
      cout << "files, PATRAN interpolates linearly between results files" << endl << flush;
      cout << "to create the extra frames.  For this reason, animations" << endl << flush;
      cout << "should contain one frame lass than an integer multiple of the" << endl << flush;
      cout << "number of WARP3D results files being animated.  The movie will" << endl << flush;
      cout << "play at 30 frames per second." << endl << flush;
      cout << "How many frames do you want in the movie?" << endl << flush;
      cout << "Frames: ";
      cin >> frames;
      cout << endl << endl << flush;
    }

  if ( format == 1 )
    {
      cout << "You have selected to create a series of .tif files." << endl << flush;
      cout << "One image will be created for each specified WARP3D result file." << endl << flush;
      cout << "These images may be bound afterward into an animated (.gif) file using" << endl << flush;
      cout << "the tif_to_gif conversion program." << endl << endl << flush;
      cout << "You will be asked to select a default name for the series of tif" << endl << flush;
      cout << "files.  The .tif extension will be added automatically.  This" << endl << flush;
      cout << "extension will be preceded by a number indicating the placement of" << endl << flush;
      cout << "an image file in the results sequence.  That is, files will be named" << endl << flush;
      cout << "*_1.tif, *_2.tif, etc., where * indicates the default name." << endl << flush;
      cout << "Please enter a name for the tif files now." << endl << flush;
      cout << "File name for *.tifs: ";
      cin >> tif_name;
      cout << endl << endl << flush;
    }
return;
}



// ---------------------------------------------------------------
// DEFORMATION OPTIONS FUNCTION
//
// Purpose: Reads deformation scale.
//
// Written by: Breanna Bailey
//
// Last modified: 8/30/01
//
// Global Variables Used:
//   1) scale = character string "True" or "Model"
//   2) scale_factor = real number defining deformation scale
//
// Local Variables Used:
//   1) deformation_scale = flag to indicate whether deformation
//                          scale is true or model 
// ---------------------------------------------------------------

void deformation_options(void)
{
  int deformation_scale;

  cout << "Deformation Options: " << endl << endl << flush;
  cout << "Please choose a deformation scale:" << endl << flush;
  cout << "   1) Model Scale." << endl << flush;
  cout << "   2) True Scale." << endl << flush;
  cout << "Deformation Scale: ";
  cin >> deformation_scale;
  if ( deformation_scale == 1 )
    strcpy( scale, "Model=");
  else
    strcpy( scale, "True=");
  cout << endl << flush;
  cout << "What scale factor do you wish to use?" << endl << flush;
  cout << "Scale Factor: ";
  cin >> scale_factor;
  cout << endl << endl << flush;
  return;
}




