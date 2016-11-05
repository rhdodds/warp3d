#!/bin/bash
#
#
#
#    SEE MAIN PROGRAM AT BOTTOM TO SET SOME LOCATIONS
#
#
function cleanup {
\rm wnd* wns* wnt* wnt* wem* state* *batch*mess* save* energy  >& /dev/null
\rm states* >& /dev/null
\rm *text   >& /dev/null
\rm woutput >& /dev/null
}


function test_voche {
echo -e "\n>>> Test 1 (Voche)"
echo      "    =============="
cleanup
echo "   ... running analysis"
$warp  < warp3d_voche.inp > woutput
ok=`grep "total job wall time (secs): " woutput | wc -l`
if [ "$ok" -le "0"  ];then
   echo "    ... analysis failed to complete"
   cd ../..
   return
fi
perl check.perl
cleanup
echo " "
}


function test_mts {
echo -e "\n>>> Test 2 (MTS)"
echo      "    ==========="
cleanup
echo "   ... running analysis"
$warp  < warp3d_mts.inp > woutput
ok=`grep "total job wall time (secs): " woutput | wc -l`
if [ "$ok" -le "0"  ];then
   echo "    ... analysis failed to complete"
   cd ../..
   return
fi
perl check.perl
cleanup
echo " "
}
  

function test_mts_multi_crystal {
echo -e "\n>>> Test 3 (MTS multicrystal)"
echo      "    ========================"
cleanup
echo "   ... running analysis"
$warp  < warp3d_mts_multi.inp > woutput
ok=`grep "total job wall time (secs): " woutput | wc -l`
if [ "$ok" -le "0"  ];then
   echo "    ... analysis failed to complete"
   cd ../..
   return
fi
perl check.perl
cleanup
echo " "
}

  

function test_ornl {
echo -e "\n>>> Test 4 (ORNL)"
echo      "    ============"
cleanup
echo "   ... running analysis"
$warp  < warp3d_ornl48.inp > woutput
ok=`grep "total job wall time (secs): " woutput | wc -l`
if [ "$ok" -le "0"  ];then
   echo "    ... analysis failed to complete"
   cd ../..
   return
fi
perl check.perl
cleanup
echo " "
}

  

function test_mrr {
echo -e "\n>>> Test 5 (mrr)"
echo      "    ==========="
cleanup
echo "   ... running analysis"
$warp  < warp3d_mrr.inp > woutput
ok=`grep "total job wall time (secs): " woutput | wc -l`
if [ "$ok" -le "0"  ];then
   echo "    ... analysis failed to complete"
   cd ../..
   return
fi
perl check.perl
cleanup
echo " "
}

function test_mrr_diff1B {
echo -e "\n>>> Test 6 (mrr - diff1B)"
echo      "    ===================="
cleanup
echo "   ... running analysis"
$warp  <  warp3d_mrr1B.inp > woutput
ok=`grep "total job wall time (secs): " woutput | wc -l`
if [ "$ok" -le "0"  ];then
   echo "    ... analysis failed to complete"
   cd ../..
   return
fi
perl check.perl
cleanup
echo " "
}


function test_mrr_diff2A {
echo -e "\n>>> Test 7 (mrr - diff2A)"
echo      "    ===================="
cleanup
echo "   ... running analysis"
$warp  <  warp3d_mrr2A.inp > woutput
ok=`grep "total job wall time (secs): " woutput | wc -l`
if [ "$ok" -le "0"  ];then
   echo "    ... analysis failed to complete"
   cd ../..
   return
fi
perl check.perl
cleanup
echo " "
}


function test_mrr_diff2B {
echo -e "\n>>> Test 8 (mrr - diff2B)"
echo      "    ===================="
cleanup
echo "   ... running analysis"
$warp  <  warp3d_mrr2B.inp > woutput
ok=`grep "total job wall time (secs): " woutput | wc -l`
if [ "$ok" -le "0"  ];then
   echo "    ... analysis failed to complete"
   cd ../..
   return
fi
perl check.perl
cleanup
echo " "
}


function test_mrr_diff3B {
echo -e "\n>>> Test 9 (mrr - diff3B)"
echo      "    ===================="
cleanup
echo "   ... running analysis"
$warp  <  warp3d_mrr3B.inp > woutput
ok=`grep "total job wall time (secs): " woutput | wc -l`
if [ "$ok" -le "0"  ];then
   echo "    ... analysis failed to complete"
   cd ../..
   return
fi
perl check.perl
cleanup
echo " "
}


function test_mrr_diff4B {
echo -e "\n>>> Test 10 (mrr - diff4B)"
echo      "    ===================="
cleanup
echo "   ... running analysis"
$warp  <  warp3d_mrr4B.inp > woutput
ok=`grep "total job wall time (secs): " woutput | wc -l`
if [ "$ok" -le "0"  ];then
   echo "    ... analysis failed to complete"
   cd ../..
   return
fi
perl check.perl
cleanup
echo " "
}


function test_djgm {
echo -e "\n>>> Test 11 (djgm)"
echo      "    ============="
cleanup
echo "   ... running analysis"
$warp  < warp3d_djgm.inp > woutput
ok=`grep "total job wall time (secs): " woutput | wc -l`
if [ "$ok" -le "0"  ];then
   echo "    ... analysis failed to complete"
   cd ../..
   return
fi
perl check.perl
cleanup
echo " "
}


#**************************************************************
#*                                                            *
#*      main programs                                         *
#*                                                            *
#**************************************************************
#

printf "\n\n**** Run verification tests for CP Model ****\n"
printf     "      (Last Updated: October 5, 2016)\n\n"
#
printf ">> making local copy of warp3d executable\n\n"
#
#  *** fix next line as needed ***
#
where_is_warp=~/warp3d_project/run_linux_em64t/warp3d.omp
if [ ! -e $where_is_warp ]
then
 printf "\n\n... set the variable: where_is_warp in the\n"
 printf     "    script and try again ...\n\n"     
 exit
fi  
cp $where_is_warp warp3d.exe
warp=`pwd`/warp3d.exe

export OMP_NUM_THREADS="10"
export MKL_NUM_THREADS="4"
export KMP_AFFINITY=scatter
printf ">> OMP and MKL threads: $OMP_NUM_THREADS $MKL_NUM_THREADS\n" 
cd voche_model
test_voche
cd ..

export OMP_NUM_THREADS="10"
export MKL_NUM_THREADS="4"
export KMP_AFFINITY=scatter
printf ">> OMP and MKL threads: $OMP_NUM_THREADS $MKL_NUM_THREADS\n" 
cd mts_model
test_mts
cd ..

export OMP_NUM_THREADS="30"
export MKL_NUM_THREADS="4"
export KMP_AFFINITY=scatter
printf ">> OMP and MKL threads: $OMP_NUM_THREADS $MKL_NUM_THREADS\n" 
cd mts_model_multiple_crystals
test_mts_multi_crystal
cd ..

export OMP_NUM_THREADS="20"
export MKL_NUM_THREADS="4"
export KMP_AFFINITY=scatter
printf ">> OMP and MKL threads: $OMP_NUM_THREADS $MKL_NUM_THREADS\n" 
cd ornl_model
test_ornl
cd ..

export OMP_NUM_THREADS="20"
export MKL_NUM_THREADS="4"
export KMP_AFFINITY=scatter
printf ">> OMP and MKL threads: $OMP_NUM_THREADS $MKL_NUM_THREADS\n" 
cd mrr_model
test_mrr
cd .. 



export OMP_NUM_THREADS="20"
export MKL_NUM_THREADS="4"
export KMP_AFFINITY=scatter
printf ">> OMP and MKL threads: $OMP_NUM_THREADS $MKL_NUM_THREADS\n" 
cd mrr_model_diff1B
test_mrr_diff1B
cd .. 


export OMP_NUM_THREADS="20"
export MKL_NUM_THREADS="4"
export KMP_AFFINITY=scatter
printf ">> OMP and MKL threads: $OMP_NUM_THREADS $MKL_NUM_THREADS\n" 
cd mrr_model_diff2A
test_mrr_diff2A
cd .. 


export OMP_NUM_THREADS="20"
export MKL_NUM_THREADS="4"
export KMP_AFFINITY=scatter
printf ">> OMP and MKL threads: $OMP_NUM_THREADS $MKL_NUM_THREADS\n" 
cd mrr_model_diff2B
test_mrr_diff2B
cd .. 


export OMP_NUM_THREADS="20"
export MKL_NUM_THREADS="4"
export KMP_AFFINITY=scatter
printf ">> OMP and MKL threads: $OMP_NUM_THREADS $MKL_NUM_THREADS\n" 
cd mrr_model_diff3B
test_mrr_diff3B
cd .. 



export OMP_NUM_THREADS="20"
export MKL_NUM_THREADS="4"
export KMP_AFFINITY=scatter
printf ">> OMP and MKL threads: $OMP_NUM_THREADS $MKL_NUM_THREADS\n" 
cd mrr_model_diff4B
test_mrr_diff4B
cd .. 



export OMP_NUM_THREADS="20"
export MKL_NUM_THREADS="4"
export KMP_AFFINITY=scatter
printf ">> OMP and MKL threads: $OMP_NUM_THREADS $MKL_NUM_THREADS\n" 
cd djgm_model
test_djgm
cd ..

printf "\n** All Done **\n" 

#
#
#
exit


