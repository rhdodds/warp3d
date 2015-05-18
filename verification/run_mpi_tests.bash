#/bin/bash

function run_all {
hybrid-test-short
hybrid-test-long
threads-hybrid-test
}

function hybrid-test-short {
echo -e "\n>>> Hypre solver. Crk growth w/ cohesive. C(T): ~1.5 mins"
echo      "    ====================================================="
cd hybrid-test-short
./run_tests_and_check
cd ..
}

function hybrid-test-long {
echo -e "\n>>> Hypre sovler. Crk growth w/ cohesive. C(T) ~8 mins"
echo      "    =================================================="
cd hybrid-test-long
./run_tests_and_check
cd ..
}

function threads-hybrid-test {
echo -e "\n>>> MPI + MKL sparse iterative (Pressure vessel)"
echo      "    ============================================"
cd threads-hybrid-test
./run_tests_and_check
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
export threads_warp3d=2

# Print the menu
printf "\n> Select a problem...\n\n"
m_testhybrid_short="Hypre solver (Crk growth w/ coheive in  C(T)): Short (2 mins)"
m_testhybrid_long="Hypre solver (Crk growth w/ coheive in  C(T)): Long (8-10 mins)"
m_testthreadshybrid="MPI + MKL sparse iterative (Pressure vessel)"
#
all="All problems"
quit="Quit"
PS3="Enter your choice (<return> to repeat menu): "
select menu_list in "$all" "$m_testhybrid_short" "$m_testhybrid_long" "$m_testthreadshybrid" "$quit"
do
      case $menu_list in
            $all)
                  run_all
                  break;;
            $m_testhybrid_short)
                  hybrid-test-short;;
            $m_testhybrid_long)
                  hybrid-test-long;;
            $m_testthreadshybrid)
                  threads-hybrid-test;;
            $quit)
                  break;;
            *) printf "You can enter only 1, .....\n";;
      esac
done

echo -e "\n>>> All done with tests ...\n\n"
exit
                  

