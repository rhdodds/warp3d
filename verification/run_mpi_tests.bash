#/bin/bash

function run_all {
test1
test2
test3
test4
test5
test6
}

function test1 {
echo -e "\n>>> Hypre solver. Crk growth w/ cohesive. C(T): ~1.5 mins"
echo      "    ====================================================="
cd hybrid_test_1
bash ./run_tests_and_check  2>/dev/null # fixes lower shell issue
cd ..
}

function test2 {
echo -e "\n>>> Hypre sovler. Crk growth w/ cohesive. C(T) ~8 mins"
echo      "    =================================================="
cd hybrid_test_2
bash ./run_tests_and_check 2>/dev/null
cd ..
}

function test3 {
echo -e "\n>>> Pardiso iterative. Assemble/solve on root (Pressure vessel)"
echo      "    =========================================================="
cd hybrid_test_3
bash ./run_tests_and_check 2>/dev/null
cd ..
}

function test4 {
echo -e "\n>>> MPI (short)soluton + combine partial results files"
echo      "    =================================================="
cd hybrid_test_4
bash ./run_tests_and_check 2>/dev/null
cd ..
}

function test5 {
echo -e "\n>>> Cluster Pardiso asymmetric (Pressure vessel)"
echo      "    ==========================================="
cd hybrid_test_5
bash ./run_tests_and_check 2>/dev/null
cd ..
}

function test6 {
echo -e "\n>>> Cluster Pardiso symmetric (Pressure vessel)"
echo      "    ==========================================="
cd hybrid_test_6
bash ./run_tests_and_check 2>/dev/null
cd ..
}

function test7 {
echo -e "\n>>> Crack growth with cells, restart and J-values)"
echo      "    ============================================="
cd hybrid_test_7
bash ./run_tests_and_check 2>/dev/null
cd ..
}
function test8 {
echo -e "\n>>> CP MTS model with restart (Taylor approximation)"
echo      "    ==============================================="
cd hybrid_test_8
bash ./run_tests_and_check 2>/dev/null
cd ..
}


#*********************************************************
#*                                                       *
#*            Main Program                               *
#*                                                       *
#*********************************************************
#
# Make sure WARP3D_HOME is set
[ -z "$WARP3D_HOME" ] && echo "Need to set WARP3D_HOME before proceeding." && exit 1

# Make sure this really is Linux
machine=`uname`
if [ "${machine:0:5}" != 'Linux' ]; then
      printf "\n>> Host does not appear to be Linux.\n"
      printf   "   The MPI version of WARP3D only runs on Linux.\n"
      printf   "   Exiting...\n"
      exit 1
fi

# Make sure our run script exists
name="$WARP3D_HOME/warp3d_script_linux_hybrid"
if [ ! -f "$name" ]; then
      printf "\n>> Expecting to find a script warp3d_script_linux_hybrid that runs the\n"
      printf   "   MPI+OpenMP version of WARP3D.\n"
      printf   "   Exiting...\n"
      exit 1
fi

# Determine number of threads etc

# Export everything
export warp3d=$name
export ranks_warp3d=4
export threads_warp3d=4
echo " "
echo "    ... Driver to run example MPI+OpenMP version of WARP3D ..."
echo "        =================================================="
echo " "
echo " These are relatively small problems than run in at most a few minutes."
echo " They exercise some key features of the code for MPI."
echo " "
echo " **** The solutions run MPI on the current node with **** "
echo " "
echo "       " $ranks_warp3d "ranks and"  $threads_warp3d" threads per rank."
echo " "
echo " **** If your computer has #cores < #ranks x # threads,"
echo "      reduce the number of threads set inside this script"
echo " "
echo " The script to run WARP3D directs MPI to run on a single node and that"
echo " MPI comminucation is via shared-memory."
echo " See file: " $name
#
# Print the menu
printf "\n> Select a problem...\n\n"
m_test1="Hypre solver (Crk growth w/ coheive in  C(T)): Short (2 mins)"
m_test2="Hypre solver (Crk growth w/ coheive in  C(T)): Long (8-10 mins)"
m_test3="Pardiso (threaded) iterative - pressure vessel: (1-2 mins)"
m_test4="Test MPI combine partial result files: Short (0.2 min)"
m_test5="Cluster Pardiso asymmetric - pressure vessel: (1-2 mins)"
m_test6="Cluster Pardiso symmetric - pressure vessel: (1-2 mins)"
m_test7="Cluster Pardiso - crack growth with cells, restart, J-values: (1-2 mins)"
m_test8="CP model with MTS (restart, Taylor approximation)"
#
all="All problems"
quit="Quit"
PS3="Enter your choice (<return> to repeat menu): "
select menu_list in "$all" "$m_test1" "$m_test2" \
    "$m_test3" "$m_test4" "$m_test5"  "$m_test6" "$m_test7"  \
    "$m_test8" "$quit"
do
      case $menu_list in
            $all)
                  run_all;;
            $m_test1)
                  test1;;
            $m_test2)
                  test2;;
            $m_test3)
                  test3;;
            $m_test4)
                  test4;;
            $m_test5)
                  test5;;
            $m_test6)
                  test6;;
            $m_test7)
                  test7;;
            $m_test8)
                  test8;;
            $quit)
                  break;;
            *) printf "You can enter only 1, .....\n";;
      esac
done

echo -e "\n>>> All done with tests ...\n\n"
exit


