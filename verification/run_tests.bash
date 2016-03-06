#!/bin/bash
#

function run_all {
test14
test18
test24
test39
test41
test44
test47
test48
test50
test51
test54
test57
test60
test61
test63
test67
test69
test70
test71
test72
test73
test74
test75
test76
test77
test78
test80
test81
}
function test14 {
echo -e "\n>>> Test 14 (Linear elastic impact, sparse iterative solver)"
echo      "    ========================================================"
cd test14
./run_tests_and_check
cd ..
}

function test18 {
echo -e "\n>>> Test 18 (Gurson model, impact loading, restart)"
echo      "    ==============================================="
cd test18
./run_tests_and_check
cd ..
}
function test24 {
echo -e "\n>>> Test 24 (skew constraints)"
echo      "    =========================="
cd test24
./run_tests_and_check
cd ..
#
}

function test39 {
echo -e "\n>>> Test 39 (CTOA crack growth)"
echo      "    ==========================="
cd test39
./run_tests_and_check
cd ..
}

function test41 {
echo -e "\n>>> Test 41 (rigid contact, large displacements)"
echo      "    ============================================"
cd test41
./run_tests_and_check
cd ..
}

function test44 {
echo -e "\n>>> Test 44 (Cohesive-interface growth, restart)"
echo      "    ============================================"
cd test44
./run_tests_and_check
cd ..
}

function test47 {
echo -e "\n>>> Test 47 (Bilinear, cyclic, thermo-plasticity)"
echo      "    ============================================="
cd test47
./run_tests_and_check
cd ..
}


function test48 {
echo -e "\n>>> Test 48 (Functionally graded materials)"
echo      "    ======================================="
cd test48
./run_tests_and_check
cd ..
}

function test50 {
echo -e "\n>>> Test 50 (Fatigue crack growth using node release)"
echo      "    ================================================="
cd test50
./run_tests_and_check
cd ..
}

function test51 {
echo -e "\n>>> Test 51 (Nonlinear patch tests: tet elements)"
echo      "    ============================================="
echo " "
cd test51
./run_tests_and_check
cd ..
}

function test54 {
echo -e "\n>>> Test 54 (FGM, non-global constraints)"
echo      "    ====================================="
cd test54
./run_tests_and_check
cd ..
}

function test57 {
echo -e "\n>>> Test 57 (SSY, mises + hydrogen efects)"
echo      "    ======================================"
cd test57
./run_tests_and_check
cd ..
}

function test60 {
echo -e "\n>>> Test 60 (Advanced cyclic plasticty model)"
echo      "    ========================================="
cd test60
./run_tests_and_check
cd ..
}

function test61 {
echo -e "\n>>> Test 61 (FGMs, interaction integrals)"
echo      "    ====================================="
cd test61
./run_tests_and_check
cd ..
}

function test63 {
echo -e "\n>>> Test 63 (360-degree SSY model, KI,KII,KIII)"
echo      "    ==========================================="
cd test63
./run_tests_and_check
cd ..
}

function test67 {
echo -e "\n>>> Test 67 (Mesh tieing of two cylinders)"
echo      "    ======================================"
cd test67
./run_tests_and_check
cd ..
}

function test69 {
echo -e "\n>>> Test 69 (Interface-cohesive elements)"
echo      "    ====================================="
cd test69
./run_tests_and_check
cd ..
}
function test70 {
echo -e "\n>>> Test 70 (cube with temperature dependent generalized plasticity)"
echo      "    ================================================================"
cd test70
./run_tests_and_check
cd ..
}

function test71 {
echo -e "\n>>> Test 71 (piston pressure loading on faces on ramp panel)"
echo      "    ========================================================"
cd test71
./run_tests_and_check
cd ..
}

function test72 {
echo -e "\n>>> Test 72 (user-defined named lists of nodes elements)"
echo      "    ===================================================="
cd test72
./run_tests_and_check
cd ..
}

function test73 {
echo -e "\n>>> Test 73 (UMAT included in WARP3D - bilinear plasticity)"
echo      "    ======================================================="
cd test73
./run_tests_and_check
cd ..
}

function test74 {
echo -e "\n>>> Test 74 (Crystal plasticity)"
echo      "    ============================"
cd test74
./run_tests_and_check
cd ..
}

function test75 {
echo -e "\n>>> Test 75 (Release nodal constraints)"
echo      "    ==================================="
cd test75
./run_tests_and_check
cd ..
}

function test76 {
echo -e "\n>>> Test 76 (User-defined multi-point constraints)"
echo      "    =============================================="
cd test76
./run_tests_and_check
cd ..
}

function test77 {
echo -e "\n>>> Test 77 (Finite strain transformations/output)"
echo      "    =============================================="
cd test77
./run_tests_and_check
cd ..
}

function test78 {
echo -e "\n>>> Test 78 (High-rate loading 20-node panel w/ hole)"
echo      "    ================================================="
cd test78
./run_tests_and_check
cd ..
}

function test80 {
echo -e "\n>>> Test 80 (Norton creep model. Non-global constraints)"
echo      "    ===================================================="
cd test80
./run_tests_and_check
cd ..
}


function test81 {
echo -e "\n>>> Test 81 (Flat and Patran result output files)"
echo      "    ============================================="
cd test81
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
 printf "\n\n> Select a problem...\n\n"
m_test14="Test 14: (Linear elastic impact, sparse iterative solver)"
m_test18="Test 18: (Gurson model growth, impact loading, restart)"
m_test24="Test 24: (hollow sphere, skew constraints, int press, large displ) "
m_test39="Test 39: (CTOA crack growth)"
m_test41="Test 41: (rigid contact, large displacements)"
m_test44="Test 44: (Cohesive-interface growth, restart)"
m_test47="Test 47: (Bilinear, cyclic, thermo-plasticity)"
m_test48="Test 48: (Functionally graded materials - nonlinear)"
m_test50="Test 50: (Fatigue crack growth using node release)"
m_test51="Test 51: (Nonlinear patch tests: tet elements)"
m_test54="Test 54: (FGM, linear, SE(T), non-global constraints)"
m_test57="Test 57: (SSY, mises + hydrogen effects)"
m_test60="Test 60: (Advanced cyclic plasticty model)"
m_test61="Test 61: (FGMs, interaction integrals)"
m_test63="Test 63: (360-degree SSY model, KI,KII,KIII)"
m_test67="Test 67: (Mesh tieing of two cylinders)"
m_test69="Test 69: (Interface-cohesive elements)"
m_test70="Test 70: (generalized_plasticity w/ temperature dependent properties)"
m_test71="Test 71: (piston loading for unsteady aerodynamic pressures)"
m_test72="Test 72: (user-defined integer lists of nodes-elements)"
m_test73="Test 73: (UMAT included in WARP3D - bilinear mises, kinematic hardening)"
m_test74="Test 74: (Crystal plasticity)"
m_test75="Test 75: (Release constraints)"
m_test76="Test 76: (User-defined multi-point constraints)"
m_test77="Test 77: (Finite strain transformations/output)"
m_test78="Test 78: (High-rate loading panel with hole. 20-node, rate plasticity)"
m_test80="Test 80: (Norton creep. Non-global constraints)"
m_test81="Test 81: (Flat and Patran result output files"
all="All problems"
quit="Quit"
PS3="Enter your choice (<return> to repeat menu): "

  select menu_list in  "$all" "$m_test14" "$m_test18" "$m_test24" "$m_test39"  \
     "$m_test41"  \
     "$m_test44" \
     "$m_test47" \
     "$m_test48" \
     "$m_test50" \
     "$m_test51" \
     "$m_test54" \
     "$m_test57" \
     "$m_test60" \
     "$m_test61" \
     "$m_test63" \
     "$m_test67" \
     "$m_test69" \
     "$m_test70" \
     "$m_test71" \
     "$m_test72" \
     "$m_test73" \
     "$m_test74" \
     "$m_test75" \
     "$m_test76" \
     "$m_test77" \
     "$m_test78" \
     "$m_test80" \
     "$m_test81" \
     "$quit"
#
  do
   case $menu_list in
          $all)
             run_all
             break;;
          $m_test14)
             test14;;
          $m_test18)
             test18;;
          $m_test24)
             test24;;
          $m_test39)
             test39;;
          $m_test41)
             test41;;
          $m_test44)
             test44;;
          $m_test47)
             test47;;
          $m_test48)
             test48;;
          $m_test50)
             test50;;
          $m_test51)
             test51;;
          $m_test54)
             test54;;
          $m_test57)
             test57;;
          $m_test60)
             test60;;
          $m_test61)
             test61;;
          $m_test63)
             test63;;
          $m_test67)
             test67;;
          $m_test69)
             test69;;
          $m_test70)
             test70;;
          $m_test71)
             test71;;
          $m_test72)
             test72;;
          $m_test73)
             test73;;
          $m_test74)
             test74;;
          $m_test75)
             test75;;
          $m_test76)
             test76;;
          $m_test77)
             test77;;
          $m_test78)
             test78;;
          $m_test80)
             test80;;
          $m_test81)
             test81;;
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


