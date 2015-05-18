#!/bin/ksh 


#  ********************************************************************
#                                                                    
#                        Filename: markup_inven.ksh    
#       
#  ********************************************************************
# 
#    Programmer: Adam Carlyle            Date: July 17,2005 
#   
#    Description:
#      Korn Shell script that searches every file in a directory for
#      the number of WARP3D markup lines contained in each file, and
#      prints a table of the inventory. This is helpful to locate the
#      the frequency and location of various machine dependencies.
#      See the WARP3D Programmer's Manual for for information on
#      markup lines.
# 
#  ********************************************************************


# Print Table Header

 echo " ------------------------------------------------------------------------------------------------------------------ "
 echo "|           Filename           |                               Number of Markup Lines                              | "
 echo "|                              | add | dec | hpu | h11 | hpi | lnx | l64 | r60 | sgi | sga | win | sgl | dbl | ptr | "
 echo " ------------------------------------------------------------------------------------------------------------------ "


# Loop over all files in current directory
 
for INPUT in *
do 
        # Print filename
        echo -n "  "; echo $INPUT
        echo -n "                                  "
        
        # Store dependency list
        num_add=`grep -c '$add' $INPUT`
        num_dec=`grep -c '#dec\|!dec' $INPUT`
        num_hpu=`grep -c '#hpu\|!hpu' $INPUT`
        num_h11=`grep -c '#h11\|!h11' $INPUT`
        num_hpi=`grep -c '#hpi\|!hpi' $INPUT`
        num_lnx=`grep -c '#lnx\|!lnx' $INPUT`
        num_l64=`grep -c '#l64\|!l64' $INPUT`
        num_r60=`grep -c '#r60\|!r60' $INPUT`
        num_sgi=`grep -c '#sgi\|!sgi' $INPUT`
        num_sga=`grep -c '#sga\|!sga' $INPUT`
        num_win=`grep -c '#win\|!win' $INPUT`
        num_sgl=`grep -c '#sgl\|!sgl' $INPUT`
        num_dbl=`grep -c '#dbl\|!dbl' $INPUT`
        num_ptr=`grep -c '#ptr\|!ptr' $INPUT`
 
        # Print dependency list
        echo -n $num_add; echo -n "     "
        echo -n $num_dec; echo -n "     "
        echo -n $num_hpu; echo -n "     "
        echo -n $num_h11; echo -n "     "
        echo -n $num_hpi; echo -n "     "
        echo -n $num_lnx; echo -n "     "
        echo -n $num_l64; echo -n "     "
        echo -n $num_r60; echo -n "     "
        echo -n $num_sgi; echo -n "     "
        echo -n $num_sga; echo -n "     "
        echo -n $num_win; echo -n "     "
        echo -n $num_sgl; echo -n "     "
        echo -n $num_dbl; echo -n "     "
        echo -n $num_ptr; echo ""
        echo " ------------------------------------------------------------------------------------------------------------------"

done 
