#
#      WARP3D verification system
#      ==========================
#
#      check results for test_61
#
$inputfile = 'test_61_out';
$line = ' ';
print "\n... Check results: $inputfile\n";
open(infile, "$inputfile") or die
"  >> Fatal Error. could not open: $inputfile\n  >> Aborting this verification segment\n\n";
print "   ... output file opened ...\n";
#
find_line( 1, " domain id: right    " );
find_line( 1, "  average      minimum      maximum" );
find_line( 1, "  average      minimum      maximum" );
#
#    line with KI value at right tip
#
$line = <infile>;
@parts = split( / +/, $line);
#
$answer = "0.1034E+01";
$partno = 1;
#
$message = " ";
if ( $answer ne $parts[$partno] ) {
 $message = "\t\t  **** difference in solution"; 
}  
#
print "   ... comparison value:       $answer","\n";
print "   ... value from output file: ", "$parts[$partno]$message\n";
#
#
#    line with KII value at right tip
#
find_line( 1, "  average      minimum      maximum" );
$line = <infile>;
@parts = split( / +/, $line);
#
$answer = "-0.6174E+00";
$partno = 1;
#
$message = " ";
if ( $answer ne $parts[$partno] ) {
 $message = "\t\t  **** difference in solution"; 
}  
#
print "   ... comparison value:       $answer","\n";
print "   ... value from output file: ", "$parts[$partno]$message\n";
#
#   T-stress at right tip
#
find_line( 1, "  average      minimum      maximum" );
find_line( 1, "  average      minimum      maximum" );
$line = <infile>;
@parts = split( / +/, $line);
#
$answer = "-0.2635E+00";
$partno = 1;
#
$message = " ";
if ( $answer ne $parts[$partno] ) {
 $message = "\t\t  **** difference in solution"; 
}  
#
print "   ... comparison value:       $answer","\n";
print "   ... value from output file: ", "$parts[$partno]$message\n";
#
#   Left tip
#
find_line( 1, " domain id: left    " );
find_line( 1, "  average      minimum      maximum" );
find_line( 1, "  average      minimum      maximum" );
#
#    line with KI value at left tip
#
$line = <infile>;
@parts = split( / +/, $line);
#
$answer = "0.5117E+00";
$partno = 1;
#
$message = " ";
if ( $answer ne $parts[$partno] ) {
 $message = "\t\t  **** difference in solution"; 
}  
#
print "   ... comparison value:       $answer","\n";
print "   ... value from output file: ", "$parts[$partno]$message\n";
#
#    line with KII value at left  tip
#
find_line( 1, "  average      minimum      maximum" );
$line = <infile>;
@parts = split( / +/, $line);
#
$answer = "-0.4057E+00";
$partno = 1;
#
$message = " ";
if ( $answer ne $parts[$partno] ) {
 $message = "\t\t  **** difference in solution"; 
}  
#
print "   ... comparison value:       $answer","\n";
print "   ... value from output file: ", "$parts[$partno]$message\n";
#
#   T-stress at left tip
#
find_line( 1, "  average      minimum      maximum" );
find_line( 1, "  average      minimum      maximum" );
$line = <infile>;
@parts = split( / +/, $line);
#
$answer = "-0.2747E+00";
$partno = 1;
#
$message = " ";
if ( $answer ne $parts[$partno] ) {
 $message = "\t\t  **** difference in solution"; 
}  
#
print "   ... comparison value:       $answer","\n";
print "   ... value from output file: ", "$parts[$partno]$message\n";
#
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