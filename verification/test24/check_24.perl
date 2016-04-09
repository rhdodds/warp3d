#
#      WARP3D verification system
#      ==========================
#
#      check results for test_24
#
$inputfile = 'test_24_out';
$line = ' ';
print "\n... Check results: $inputfile\n";
open(infile, "$inputfile") or die
"  >> Fatal Error. could not open: $inputfile\n  >> Aborting this verification segment\n\n";
print "   ... output file opened ...\n";
#
find_line( 1, " output stresses 1 8" );
find_line( 1, " c2                  c3" );
#
#    line with element 1 stresses
#
$line = <infile>;$line = <infile>;
@parts = split( / +/, $line);
$answer = "-205.047360";
$partno = 3;
#
$message = " ";
if ( $answer ne $parts[$partno] ) {
 $message = "\t\t  **** difference in solution"; 
}  
#
print "   ... comparison value:       $answer","\n";
print "   ... value from output file: ", "$parts[$partno]$message\n";


$line = <infile>;$line = <infile>;$line = <infile>;
@parts = split( / +/, $line);
$answer = "151.957661";
$partno = 1;
#
$message = " ";
if ( $answer ne $parts[$partno] ) {
 $message = "\t\t  **** difference in solution"; 
}  
#
print "   ... comparison value:       $answer","\n";
print "   ... value from output file: ", "$parts[$partno]$message\n";

print "   ... done\n";
close infile;
exit;


#**********************************************************
#*                                                        *
#*  find_line    ( $line is global )                      *
#*                                                        *
#**********************************************************

sub find_line {
      my ( $type, $string ) = @_;
      my ( $debug );
      $debug = 0;
#
      if( $debug == 1 )
	{
         print " type: ", $type, "\n";
         print " string: ", $string, "\n";
	}
#
      while ( !eof(infile) )
        {
           $line = <infile>; 
	   if( $line =~ /$string/ ) {return};
        }
#
	 print "\n>>> Fatal Error. string search type: ",$type;
         print "\n    Searching for string: ", "\"",$string,"\" failed";
	 print "\n    Aborting this verification segment\n\n";
         exit;
}


#**********************************************************
#*                                                        *
#*  skip_lines                                            *
#*                                                        *
#**********************************************************

sub skip_lines {
      my ( $type, $nlines ) = @_;
      my ( $line, $debug, $count );
      $debug = 0;
#
      if( $debug == 1 )
	{
         print " type: ", $type, "\n";
         print " nlines: ", $nlines, "\n";
	}
#
      $count = 0;
      while ( !eof(infile) )
        {
         $line = <infile>; $count++;
         if( $count == $nlines ) {return};
        }
#
	 print "\n>>> Fatal Error. EOF reached beofre skip lines type: ",$type;
	 print "\n    Aborting this verification segment\n\n";
         exit;
}
