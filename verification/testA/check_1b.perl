#
#      WARP3D verification system
#      ==========================
#
#      check results for test_1b
#
$inputfile = 'out_1b';
print "\n\t... Check results: $inputfile\n";
open(infile, "$inputfile") or die
"  >> Fatal Error. could not open: $inputfile\n  >> Aborting this verification segment\n\n";
print "\t   ... output file opened ...\n";
#
find_line( 1, "average      minimum" );
#
#    get line with domain component values. print value
#
$line = <infile>;  @parts = split( / +/, $line);
#
$answer = "0.2574E-01";
$partno = 1;
#
$message = " ";
if ( $answer ne $parts[$partno] ) {
 $message = "\t\t  **** difference in solution"; 
}  
#
print "\t   ... comparison value:       $answer","\n";
print "\t   ... value from output file: ", "$parts[$partno]$message\n";

print "\t   ... done\n";
close infile;
exit;


#**********************************************************
#*                                                        *
#*  find_line                                             *
#*                                                        *
#**********************************************************

sub find_line {
      my ( $type, $string ) = @_;
      my ( $line, $debug );
      $debug = 0;
#
      if( $debug == 1 )
	{
         print " type: ", $type, "\n";
         print " string: ", $string, "\n";
         print " file: ", $file, "\n";
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