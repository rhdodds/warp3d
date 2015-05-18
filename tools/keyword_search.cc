#include <cmath>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
 int main ()
{
  fstream keyword_file, file_to_search;
  string keyword_file_name, file_to_search_name;
  
  cout << endl << endl << endl 
       << "Search a file for a list of keywords..." 
       << endl << endl ;
  cout << "File with list of keywords: "; cin >> keyword_file_name;
  cout << "File to search: "; cin >> file_to_search_name;

  keyword_file.open( keyword_file_name.c_str(),ios::in );
  file_to_search.open( file_to_search_name.c_str(), ios::in );

  if( !keyword_file.is_open() )
    { cout <<".. keyword file not opened.."; exit(0); } 
  if( !file_to_search.is_open() )
    { cout <<".. file to search not opened.."; exit(0); }
  cout << endl << "..files opened ok.." << endl;
//
  string  search_word(50,' '), line(80,' ');
//
//         loop over all keywords. for each keyword, loop over the
//         file to be searched
//
  while( true )
    {
      if( !getline( keyword_file, search_word ) )
	{ cout << endl << "..All Done.." << endl << endl; exit(0); }
      cout << endl << "Next search word: " << search_word << endl;
      file_to_search.clear(); file_to_search.seekg(ios::beg);
      int line_count = 0;
      while( true )
	{
          if( !getline( file_to_search, line ) )
           { cout << endl << " .. done with keyword .." << endl; break; }
          line_count++;
	  int pos = line.find( search_word );
	  if( pos == string::npos ) continue;
          cout << " .. keyword match @ line: "
               << line_count << " " << line << endl;
       }
    }

  exit(0);

}
 
