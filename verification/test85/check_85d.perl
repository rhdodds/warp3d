#
#      WARP3D verification system
#      ==========================
#
#      check results for Voche
#
$inputfile = 'test_85d_out';
$line = ' ';
open(infile, "$inputfile") or die
"  >> Fatal Error. could not open: $inputfile\n  >> Aborting this verification segment\n\n";
print "   ... output file opened ...\n";
#
#
find_line( 1, " step no.      2 " );
find_line( 1, "1       1" );
#
@parts = split( / +/, $line);
#
$answer = "3.333333";
$partno = 3;
#
$parts[$partno] =~ s/\x0d{0,1}\x0a\Z//s;
$out_value = $parts[$partno];
#
$message = " ";
if ( $answer ne $parts[$partno] ) {
 $message = "\t\t  **** difference in solution";
}
#
print "   ... comparison value:       $answer","\n";
print "   ... value from output file: ", "$out_value$message\n";
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
	 print "\n>>> Fatal Error. EOF reached before skip lines type: ",$type;
	 print "\n    Aborting this verification segment\n\n";
         exit;
}
