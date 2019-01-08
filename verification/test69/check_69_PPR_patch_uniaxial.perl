#
#      WARP3D verification system
#      ==========================
#
#      check results for test_69_PPR_patch_uniaxial
#
$inputfile = 'test_69_PPR_patch_uniaxial_out';
print "\n... Check results: $inputfile\n";
open(infile, "$inputfile") or die
"  >> Fatal Error. could not open: $inputfile\n  >> Aborting this verification segment\n\n";
print "   ... output file opened ...\n";
#
find_line( 1, "   newton convergence tests step:     140 iteration:   3 " );
find_line( 1, " -------   --------  " );

#skip_lines( 1, 1 );
#
#  T-normal after step 140
#
$line = <infile>;  @parts = split( / +/, $line);
print "   ... comparison value:       0.405393E-01\n";
print "   ... value from output file: ", "$parts[2]\n";
print "   ... done\n";
close infile;
exit;
#
print "   ... done\n";
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