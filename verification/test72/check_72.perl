#
#      WARP3D verification system
#      ==========================
#
#      check results for test_72
#
$inputfile = 'test72_out';
$line = ' ';
print "\n... Check results: $inputfile\n";
open(infile, "$inputfile") or die
"  >> Fatal Error. could not open: $inputfile\n  >> Aborting this verification segment\n\n";
print "   ... output file opened ...\n";
#
#         -------------    test  1    -------------------
#
#find_line( 1, "  1 u = 0." );
find_line( 2, " >>>>> number of nodes" );
@parts = split( / +/, $line);
#
$answer = "853";
$partno = 7;
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
print "\n";
#
#         -------------    test 2    -------------------
#
#find_line( 2, "for proximity");
find_line( 2, " number of nodes" );
@parts = split( / +/, $line);
#
$answer = "21";
$partno = 8;
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
print "\n";
#
#         -------------    test 3    -------------------
#
find_line( 2, " >>>>> number of nodes" );
@parts = split( / +/, $line);
#
$answer = "821";
$partno = 8;
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
print "\n";
#
#         -------------    test 4    -------------------
#
find_line( 2, " >>>>> number of nodes" );
@parts = split( / +/, $line);
#
$answer = "1";
$partno = 7;
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
print "\n";
#
#         -------------    test 5    -------------------
#
find_line( 2, " >>>>> number of nodes" );
@parts = split( / +/, $line);
#
$answer = "821";
$partno = 7;
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
print "\n";
#
#         -------------    test 6    -------------------
#
find_line( 2, " >>>>> number of nodes" );
@parts = split( / +/, $line);
#
$answer = "2053";
$partno = 7;
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
print "\n";
#
#         -------------    test 7    -------------------
#
find_line( 2, " >>>>> number of nodes" );
@parts = split( / +/, $line);
#
$answer = "51";
$partno = 7;
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
print "\n";
#
#         -------------    test 8    -------------------
#
find_line( 2, " >>>>> number of nodes" );
@parts = split( / +/, $line);
#
$answer = "1";
$partno = 7;
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
#
print "   ... done\n";
close infile;
#
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
#           print "..... line:  ",$line, "\n";
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