#!/bin/bash
#

function run_all {
testA
testB
testC
testD
testE
testF
testG
testH
testI
testJ
testK
testL
testM
testN
testO
testP
testQ
testR
}
function testA {
echo -e "\n>>> Test A (SET, LEFM, Thermal, 20-node, FGMs, Face Loading, P-strain)"
echo      "    ===================================================================="
cd testA
./run_tests_and_check
cd ..
}

function testB {
echo -e "\n>>> Test B (Same as A with globally rotated model)"
echo      "    ============================================="
cd testB
./run_tests_and_check
cd ..
}
function testC {
echo -e "\n>>>  Test C (SET, Small Eps, Thermal, 20-node, Rotated Model, FGMs, P-strain)"
echo      "    ========================================================================="
cd testC
./run_tests_and_check
cd ..
#
}

function testD {
echo -e "\n>>> Test D (SET, NLGEOM, Thermal, 20-node, Rotated Model, FGMs, SS-curves, Focused Mesh)"
echo      "    ===================================================================================="
cd testD
./run_tests_and_check
cd ..
}

function testE {     # MPI Linux only
echo -e "\n>>> Test E (MPI-SET, NLGEOM, Thermal, 20-node, Rotated Model)"
echo      "    ========================================================"

if [ "${machine:0:6}" = 'Darwin' ]; then
       printf "\n\t>> Cannot run MPI jobs on OSX\n"
       return
fi
if [ "$MACHINE_TYPE" = '2' ]; then
       printf "\n\t>> Cannot run MPI jobs on Windows\n"
fi

cd testE
./run_tests_and_check
cd ..
}

function testF {
echo -e "\n>>> Test F (SET, LEFM, Thermal, 8-node, P-strain, FGMs, blunt front)"
echo      "    ================================================================"
cd testF
./run_tests_and_check
cd ..
}
function testG {
echo -e "\n>>> Test G (SET,  Small Eps, Thermal, 8-node, P-strain, FGMs, blunt front)"
echo      "    ======================================================================"
cd testG
./run_tests_and_check
cd ..
}
function testH {
echo -e "\n>>> Test H (SEB, Small Eps, Bending, 8-node, deformation/mises, P. Strains, FGMs, blunt front)"
echo      "    =========================================================================================="
cd testH
./run_tests_and_check
cd ..
}
function testI {
echo -e "\n>>> Test I (SEB,  Small Eps, Bending, 20-node, deformation/mises, P-strain, FGMs, blunt front)"
echo      "    =========================================================================================="
cd testI
./run_tests_and_check
cd ..
}

function testJ {
echo -e "\n>>> Test J (SEB, Small/NLGEOM, Bending, 8-node, P-strain, FGMs, blunt front)"
echo      "    ========================================================================"
cd testJ
./run_tests_and_check
cd ..
}
function testK {
echo -e "\n>>> Test K (SET, Small/NLGEOM, Thermal, 8-node, P-strain, FGMs, blunt front)"
echo      "    ========================================================================"
cd testK
./run_tests_and_check
cd ..
}
function testL {
echo -e "\n>>> Test L (SEB, Small/NLGEOM, Bending, 20-node, P-strain, FGMs, blunt front, refined mesh)"
echo      "    ======================================================================================="
cd testL
./run_tests_and_check
cd ..
}
function testM {
echo -e "\n>>> Test M (SET, Small/NLGEOM, Thermal, 20-node, P-strain, FGMs, blunt front, refined mesh)"
echo      "    ======================================================================================="
cd testM
./run_tests_and_check
cd ..
}
function testN {
echo -e "\n>>> Test N (SEB, Small/NLGEOM, Bending, 8-node, 3D[3-layers], FGMs, blunt front)"
echo      "    ============================================================================"
cd testN
./run_tests_and_check
cd ..
}
function testO {
echo -e "\n>>> Test O (SEB, residual, initial state, 8-node P Strain)"
echo      "    ======================================================"
cd testO
./run_tests_and_check
cd ..
}
function testP {
echo -e "\n>>> Test P (Pipe, weld eigenstrain method, 8-node P Strain)"
echo      "    ======================================================="
cd testP
./run_tests_and_check
cd ..
}
function testQ {
echo -e "\n>>> Test Q (Simulated weld, release nodes, load as M(T))"
echo      "    ======================================================="
cd testQ
./run_tests_and_check
cd ..
}
function testR {
echo -e "\n>>> Test R (Crystal plasticity, SSY- P. strain)"
echo      "    ==========================================="
cd testR
./run_tests_and_check
cd ..
}

function get_num_threads {
#
 printf ">> number of threads to use? "
 read i
 export NUM_THREADS="$i"
}


#**************************************************************
#*                                                            *
#*      set_platform environment variables for passing down   *
#*                                                            *
#**************************************************************

function set_platform {
#
if [ "${machine:0:6}" = 'Darwin' ]; then
    printf " ---->>> this machine is OSX...\n"
    name="$WARP3D_HOME/warp3d_script_mac_os_x"
    if [ ! -f "$name" ]
     then
       printf "\n>> Expecting to find a warp3d_script_mac_os_x file that runs WARP3D\n"
       printf   "   Enter full path name for a Bash shell script that runs the\n"
       printf   "   Mac OSX version of WARP3D, and that takes 1 parameter to set\n"
       printf   "   the number of threads for execution...\n"
       printf   ">> File name: "
       read name; printf "\n"
       if [ ! -f "$name" ]
         then
            echo -e "\n> Cannot find WARP3D executable. Exiting.. Sorry\n"
            exit
     fi
    fi
    export WARP3D_EXE="$name"
    get_num_threads  # sets NUM_THREADS
    export MACHINE_TYPE="0"
    echo ">> Running jobs with threads = "$NUM_THREADS
    return
fi
#
#
if [ "${machine:0:5}" = 'Linux' ]; then
    printf " --->>> this machine is Linux...\n"
    name="$WARP3D_HOME/warp3d_script_linux_openmp"
    if [ ! -f "$name" ]
     then
       printf "\n>> Expecting to find a warp3d_script_linux_openmp file that runs WARP3D\n"
       printf   "   Linux version of WARP3D, and that takes 1 parameter to set\n"
       printf   "   the number of threads for execution...\n"
       printf   ">> File name: "
       read name
       if [ ! -f "$name" ]
         then
          echo -e "\n> Cannot find WARP3D executable. Exiting.. Sorry\n"
          exit
       fi
    fi
    export WARP3D_EXE="$name"
    get_num_threads  # sets NUM_THREADS
    export MACHINE_TYPE="0"
    echo ">> Running jobs with threads = "$NUM_THREADS
    return
fi
#
#          for Windows, we need the "" enclosing names since
#          it is likely that path names will have spaces in them
#
if [ "${machine:0:6}" = 'CYGWIN' ]; then
    mtype=2
    printf ">> Windows 64-bit selected for testing...\n"
    name="$WARP3D_HOME/warp3d_script_windows"
    if [ ! -f "$name" ]
      then
       printf "\n>> Expecting to find a warp3d_script_windows file\n"
       printf   "   (Bash shell) that runs WARP3D, and that takes 1 parameter to set\n"
       printf   "   the number of threads for execution...\n"
       printf   ">> File name: "
       read name
       if [ ! -f "$name" ]
         then
          echo -e "\n> Cannot find WARP3D executable. Exiting.. Sorry\n"
          exit
       fi
    fi
    export WARP3D_EXE_64="$name"
    get_num_threads  # sets NUM_THREADS
    export MACHINE_TYPE=$mtype
    echo ">> Running jobs with threads = "$NUM_THREADS
    return
fi
#
printf "\n\n **** unrecognized machine type. Quitting ****\n"
printf " **** Sorry ....\n\n"

exit
#
}
#**************************************************************
#*                                                            *
#*      main programs                                         *
#*                                                            *
#**************************************************************
#

# Check to make sure the WARP3D_HOME variable is set
[ -z "$WARP3D_HOME" ] && echo "Need to set WARP3D_HOME before proceeding." && exit 1

#
machine=`uname`
printf ">> machine id from uname: $machine"
#
set_platform
#
if [ "$MACHINE_TYPE" = '0' ]; then
   echo -e "\n>> OSX and Linux verification..."
fi
#
printf "\n>> Note: comparison and output values may be different in the last 1 or 2"
printf "\n>> ----- signficant digits with various number of threads"
#
if [ "$MACHINE_TYPE" = '2' ]; then
   echo -e "\n>> Windows 64 verification..."
fi
#
 printf "\n\n> Select a problem (models w/ FGM subcases force J7, J8 computation)...\n\n"

m_testA="Test A: SE(T), LEFM, Thermal, 20-node, FGMs, Face Loading, P-strain"
m_testB="Test B: Same as A with globally rotated model"
m_testC="Test C: SE(T), Small Eps, Thermal, 20-node, Rotated Model, FGMs, P-strain"
m_testD="Test D: SE(T), GEONL, Thermal, 20-node, Rotated Model, P-strain, FGMs, Focused Mesh"
m_testE="Test E: MPI - SE(T), GEONL, Thermal, 20-node, Rotated Model, P-strain"
m_testF="Test F: SE(T), LEFM, Thermal, 8-node, P-strain, FGMs, blunt front"
m_testG="Test G: SE(T), Small Eps, Thermal, 8-node, P-strain, FGMs, blunt front"
m_testH="Test H: SE(B), Small Eps, Bending, 8-node, deformation/mises, P. Strains, FGMs, blunt front"
m_testI="Test I: SE(B), Small Eps, Bending, 20-node, deformation/mises, P-strain, FGMs, blunt front"
m_testJ="Test J: SE(B), Small/NLGEOM, Bending, 8-node, P-strain, FGMs, blunt front, refined mesh"
m_testK="Test K: SE(T), Small/NLGEOM, Thermal, 8-node, P-strain, FGMs, blunt front, refined mesh"
m_testL="Test L: SE(B), Small/NLGEOM, Bending, 20-node, P-strain, FGMs, blunt front, refined mesh"
m_testM="Test M: SE(T), Small/NLGEOM, Thermal, 20-node, P-strain, FGMs, blunt front, refined mesh"
m_testN="Test N: SE(B), Small/NLGEOM, Bending, 8-node, 3D:3-layers, FGMs, blunt front, refined mesh"
m_testO="Test O: SE(B), Residual Stresses, Initial-State, 8-node, P. Strain"
m_testP="Test P: Pipe, Weld, Eigenstrain, 8-node, P. Strain"
m_testQ="Test Q: Simulated weld, release nodes, load as M(T), 8-node P. Strain"
m_testR="Test R: Crystal plasticity, SSY- P. strain"
all="All problems"
quit="Quit"
PS3="Enter your choice (<return> to repeat menu): "

  select menu_list in  "$all" "$m_testA" "$m_testB"  \
     "$m_testC"  \
     "$m_testD" \
     "$m_testE" \
     "$m_testF" "$m_testG" "$m_testH" "$m_testI" "$m_testJ" "$m_testK" "$m_testL" \
     "$m_testM"  "$m_testN"  "$m_testO"  "$m_testP"   "$m_testQ" "$m_testR" "$quit"
#
  do
   case $menu_list in
          $all)
             run_all
             break;;
          $m_testA)
             testA;;
          $m_testB)
             testB;;
          $m_testC)
             testC;;
          $m_testD)
             testD;;
          $m_testE)
             testE;;
          $m_testF)
             testF;;
          $m_testG)
             testG;;
          $m_testH)
             testH;;
          $m_testI)
             testI;;
          $m_testJ)
             testJ;;
          $m_testK)
             testK;;
          $m_testL)
             testL;;
          $m_testM)
             testM;;
          $m_testN)
             testN;;
          $m_testO)
             testO;;
          $m_testP)
             testP;;
          $m_testQ)
             testQ;;
          $m_testR)
             testR;;
          $quit)
           break;;
          *) printf "You can enter only 1, .....\n";;
   esac
  done
#
#
#
echo -e "\n>>> All done with tests ...\n\n"
exit


