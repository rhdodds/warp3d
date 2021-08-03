#!/bin/bash
#
#     Makewarp.bash (5th version)
#
#     modified: July 2021
#
#     Description:
#
#           Bash script to interactively drive compilation of Linux and Mac
#           versions of WARP3D.  For Windows, we print message and quit.
#           Run in a Bash shell:  Makewarp.bash
#
#           Mac: only user selectable option is to use gfortran or Intel Fortran.
#                script hunts down the MKL files, checks for versions
#                of gfortran or Intel compilers, and runs the makefile
#                to build the threads (OpenMP) executable.
#
#           Linux: .....  For
#           Linux can either do that (simple mode) or prompt interactively for
#           the libraries, compiler, etc.
#
#           Jan 2019. Now supports gfortran and Intel Fortran.
#                     Uses (free) MKL library from Intel
#
#
#      Main program (function) at bottom of this script
#
# ****************************************************************************
#
#   Function: select Fortran compile to build WARP3D
#
# ****************************************************************************

function select_Fortran_compiler
{
printf "\nSelect build compiler to use:\n"
PS3="Select choice: "
select opt in 'gfortran' 'Intel Fortran' 'Exit'
do
      case $REPLY in
            1 )   printf "\n"
                  GFORTRAN=yes
                  INTEL_FORTRAN=no
                  break
                  ;;
            2 )   printf "\n"
                  GFORTRAN=no
                  INTEL_FORTRAN=yes
                  break
                  ;;
            3 )   exit 0
                  ;;
      esac
done
#
return
}

# ****************************************************************************
#
#   Function: Check Intel Fortran compiler exists & version
#
# ****************************************************************************
#
function check_Intel_Fortran_setup
{
if [ "$INTEL_FORTRAN" = "no" ]; then
 return
fi
#
hash ifort 2>&- || {
printf "[ERROR]\n"
printf "... Cannot find the Intel Fortran compiler (ifort) in your PATH.\n"
printf "... See Intel Fortran install documentation. Most often a line of\n"
printf "... the form: source /opt/intel/... is placed\n"
printf "... in the /etc/bashrc our your ~/.bashrc file.\n"
printf "Quitting...\n\n"
exit 1
}
/bin/rm zqq03 >& /dev/null
ifort --version >zqq03
sed -i -e '2,10d' zqq03
echo -e "\n ... Intel Fortran detected:" `cat zqq03`
count2=`grep "2021.3" zqq03 |wc -l`
count3=`grep "2021.4" zqq03 |wc -l`
/bin/rm zqq03*
ok=0
if [ $count2 -eq "1" ]; then
   ok=1
fi
if [ $count3 -eq "1" ]; then
   ok=1
fi
if [ $ok -eq "0" ]; then
    printf "\n... ERROR: ifort must be one of these versions:"
    printf "\n... 2021.3 or newer"
    printf "\n... other versions have known bugs that affect WARP3D"
    printf "\n... Quitting...\n\n"
    exit 1
fi
#
return
}

# ****************************************************************************
#
#   Function: Check gfortran exists and version 
#
# ****************************************************************************
#
function check_gfortran_setup
{
if [ "$GFORTRAN" = "no" ]; then
    return
fi
#
hash gfortran 2>&- || {
  printf "[ERROR]\n"
  printf "... Cannot find the gfortran compiler in your PATH.\n"
  printf "Quitting...\n\n"
  exit 1
}
#
version=$(echo `gfortran --version` | cut -c 18-20)
#
if [ "$version" -lt "7" ]; then
      printf "\n... ERROR: gfortran must be version 7 or newer"
      echo -e  "\n...        you have version: " $version
      printf "...        from:  gfortran --version"
      printf "\n... Quitting...\n\n"
      exit 1
fi
#
return
#
}
# ****************************************************************************
#
#   Function: Install MPI code
#
# ****************************************************************************

function install_mpi
{

#    Check to see if MPI code is installed.  If not,
#    copy true MPI versions of MPI code to current
#    directory for parallel execution using MPI.

     mpi_match1=`cmp ./mpi_code.f $mpi_dir/mpi_code_real.f | wc -l`
     mpi_match2=`cmp ./mpi_handle_slaves.f $mpi_dir/mpi_handle_slaves_real.f | wc -l`    
#     
     if [ ! $mpi_match1 = "0"  ]; then
        cd $mpi_dir
        ./install_mpi
        cd $WARP3D_HOME/src
     elif [ ! $mpi_match2 = "0"  ]; then
        cd $mpi_dir
        ./install_mpi
        cd $WARP3D_HOME/src
     else
        printf " > MPI source code already installed...\n"
     fi
}


# ****************************************************************************
#
#   Function: Uninstall MPI code
#
# ****************************************************************************

function uninstall_mpi
{

#    Check to see if serial (dummy MPI) code is installed.
#    If not, copy dummy versions of MPI code to current
#    directory for serial execution.

     mpi_match=`cmp ./mpi_code.f $mpi_dir/mpi_code_dummy.f | wc -l`
     if [ ! $mpi_match = "0" ]; then
        cd $mpi_dir
        ./uninstall_mpi
        cd $WARP3D_HOME/src
     else
        printf " > MPI source code already uninstalled...\n"
     fi
}

# ****************************************************************************
#
#   Function: Install hypre code
#
# ****************************************************************************

function install_hypre
{

#    Check to see if hypre code is installed.  If not,
#    copy in the true versions of the code

     hypre_match=`cmp ./iterative_sparse_hypre.f $hypre_dir/src_sparse_hypre.f | wc -l`
     if [ ! $hypre_match = "0" ]; then
        cd $hypre_dir
        ./install_hypre
        cd $WARP3D_HOME/src
     else
        printf " > hypre source code already installed...\n"
     fi
}

# ****************************************************************************
#
#   Function: Uninstall hypre code
#
# ****************************************************************************

function uninstall_hypre
{

#    Check to see if dummy code is already installed.  If not, install it.
     hypre_match=`cmp ./iterative_sparse_hypre.f $hypre_dir/dummy_sparse_hypre.f | wc -l`
     if [ ! $hypre_match = "0" ]; then
        cd $hypre_dir
        ./uninstall_hypre
        cd $WARP3D_HOME/src
     else
        printf " > hypre source code already uninstalled...\n"
     fi
}

# ****************************************************************************
#
#   Function: Issue message to use Makewarp.bat to compile on Windows
#
# ****************************************************************************

function print_windows_message
{
  printf "\n>> To compile on Windows, exit this \n"
  printf "   script and run the 'Makewarp.bat' batch file from within a \n"
  printf "   Windows command prompt shell. The shell must be setup \n"
  printf "   64-bit building with the Intel Fortran Compser suite.\n\n"

  printf "   The Intel Fortran Compiler creates a link to the proper \n"
  printf "   command prompt build environment in the Windows 'start' menu. \n\n"
}

# ****************************************************************************
#
#     Function:   Global defaults and branch point for Linux{openmp,hybrid}
#                 {simple,custom}
#
# ****************************************************************************

function linux_main
{
#
# These are common defaults for hybrid and openmp
#
HYPRE_ROOT=$WARP3D_HOME/linux_packages
INTEL_FORTRAN=no
GFORTRAN=no
MKLQ=yes
#
# Branch on mpi/openmp
#
printf "Compile OpenMP-only version or hybrid (MPI/OpenMP) version?\n"
PS3="Select choice: "
LIST="OpenMP Hybrid Exit"
select OPT in $LIST
do
      case $OPT in
            "Hybrid")
                  printf "\n"
                  COMPILER=mpiifort
                  ALTCOMPILER=mpiifort
                  MAKEFILE=Makefile.linux_Intel.mpi_omp
                  INTEL_FORTRAN=yes
                  MPIQ=yes
                  break
                  ;;
            "OpenMP")
                  printf "\n"
                  COMPILER=ifort    # default
                  ALTCOMPILER=ifort
                  MAKEFILE=Makefile.linux_Intel.omp
                  MPIQ=no
                  break
                  ;;
            "Exit")
                  printf "\n"
                  exit 0
                  ;;
          esac
done
#
if [ $OPT = "OpenMP" ]; then
  select_Fortran_compiler
  if [ $GFORTRAN = "yes" ]; then
     COMPILER=gfortran
     ALTCOMPILER=gfortran
     MAKEFILE=Makefile.linux_gfortran.omp
  fi
fi      
check_Intel_Fortran_setup
check_gfortran_setup
#
# handle gfortran separately since only OpenMP available.
#
if [ $GFORTRAN = "yes" ]; then
    compile_linux_gfortran_OpenMP
    return
fi      
#
# Intel Fortran for build OpenMP or MPI + OpenMP.
#
# Branch on whether we want simple or interactive mode
#
      printf "Simple (defaults, no prompts) or advanced (prompt) mode?\n"
      PS3="Select choice: "
      LIST="Simple Advanced Help Exit"
      select OPT in $LIST
      do
            case $OPT in
                  "Simple")
                        printf "\n"
                        linux_simple
                        break
                        ;;
                  "Advanced")
                        printf "\n"
                        linux_advanced
                        break
                        ;;
                  "Help")
                        printf "\n"
                        linux_help
                        ;;
                  "Exit")
                        printf "\n"
                        exit 0
                        ;;
            esac
      done
return      
}
# ****************************************************************************
#
#     Function:   Prints a help message for the linux compiles
#
# ****************************************************************************
function linux_help
{

      printf "Simple mode will take built-in defaultsems),\n"
      printf "and checks to ensure they will work on your system.\n"
      printf "Advanced mode will prompt you for all the configuration variables, but does\n"
      printf "provide default values.\n\n"

      printf "Simple (defaults, no prompts) or advanced (prompt) mode?\n"
      printf "1) Simple\n"
      printf "2) Advanced\n"
      printf "3) Help\n"
      printf "4) Exit\n"
return
}

# ****************************************************************************
#
#     Function:   Simple (defaults with checks) mode for Linux
#
# ****************************************************************************
function linux_simple
{
#
# Tell the user what this does
#
      printf "... Setting default options & performing checks to \n"
      printf "... ensure your system and the WARP3D directories \n"
      printf "... are configured correctly.\n"
#
#   Is this really a Linux system?
#
      match=`uname | grep Linux | wc -l`
      if [ $match = "0" ]; then
            printf "[ERROR]\n"
            printf "This is not a Linux system.\n Quitting...\n\n"
            exit 1
      fi
##
# Just make hypre by default
#
      if [ "$MPIQ" = "yes" ]; then
            HYPQ=yes
      else
            HYPQ=no
      fi
#
# Test existence of mpi compiler if required
#
      if [ "$MPIQ" = "yes" ]; then
            hash mpiifort 2>&- || {
                  printf "[ERROR]\n"
                  printf "Cannot find the Intel MPI Fortran compiler (mpiifort) in your PATH.\n"
                  printf "Have you sourced the Intel compiler variables?\n"
                  printf "Do you have Intel MPI installed?\n"
                  printf "Quitting...\n"
                  printf "\n"
                  exit 1
            }
      fi
#
# Test libHYPRE (if we need it)
#
      if [ "$HYPQ" = "yes" ]; then
            if [ ! -e $WARP3D_HOME/linux_packages/lib/libHYPRE.a ]; then
                  printf "[ERROR]\n"
                  printf "Cannot find hypre library in $WARP3D_HOME/linux_packages/lib.\n"
                  printf "Is it compiled?\n"
                  printf "\n"
                  exit 1
            fi
      fi

}

# ****************************************************************************
#
#     Function:   Advanced (prompt but don't check) mode for Linux
#
# ****************************************************************************
function linux_advanced
{
      # Just make hypre by default still (we would have linker errors)
      if [ "$MPIQ" = "yes" ]; then
            HYPQ=yes
      else
            HYPQ=no
      fi

      printf "\nIn the following prompts enter the values you want for each setting.\n"
      printf "If you leave a line blank the value will be set to the default.\n"
      printf "Note: these values will not be checked for correctness as in the simple mode.\n\n"

      printf "Fortran compiler\n"
      printf "Default: $COMPILER\n"
      read -p ": " OPT
      [ -n "$OPT" ] && COMPILER=$OPT
      printf "\n"
#
      if [ $HYPQ = "yes" ]; then
            printf "hypre root directory\n"
            printf "Default: $HYPRE_ROOT\n"
            read -p ": " OPT
            [ -n "$OPT" ] && HYPRE_ROOT=$OPT
            printf "\n"
      fi
}

# ****************************************************************************
#
#     Function:   Compile WARP3D for linux using gfortran (OpenMP)
#
# ****************************************************************************

function compile_linux_gfortran_OpenMP
{
#
# Start by going through the packages and installing them if required
#
printf "... Building OpenMP Linux w/ gfortran ...\n"
uninstall_mpi
uninstall_hypre
#
# Setup the directory structure if required
#
if [ ! -d ../run_linux ]; then
    mkdir ../run_linux
    printf ">> Making run_linux directory...\n"
fi
if [ ! -d ../obj_linux_gfortran_omp ]; then
    mkdir ../obj_linux_gfortran_omp
    printf ">> Making obj_linux_gfortran_omp directory...\n"
fi
#
#   prompt the user for the number of concurrent compile processes to use
#
printf " \n"
read -p " > Number of concurrent compile processes allowed? (default 1): " JCOMP
   [ -z "$JCOMP" ] && JCOMP=1
#
#   touch main program so it will always be re-compiled (will include compile
#   date & time in warp3d hearder block)
#
touch main_program.f  
#
#   run the makefile for Linux. we now pass more parameters to the makefile
#
printf "... Starting make program for Linux  OpenMP w/ gfortran.... \n\n"
echo "makefile: " $MAKEFILE
make -j $JCOMP -f $MAKEFILE 
#
return
}


# ****************************************************************************
#
#     Function:   Compile WARP3D for linux
#
# ****************************************************************************

function compile_linux_Intel {
	
#
# Start by going through the packages and installing them if required
#
      printf "\n"

      if [ $MPIQ = "yes" ] ; then
            install_mpi
      else
            uninstall_mpi
      fi

      if [ $HYPQ = "yes" ] ; then
            install_hypre
      else
            uninstall_hypre
      fi
      printf "\n"
#
# Setup the directory structure if required
#
      if [ ! -d ../run_linux ]; then
            mkdir ../run_linux
            printf ">> Making run_linux directory...\n"
      fi
      if [ ! -d ../obj_linux_Intel_omp ]; then
            mkdir ../obj_linux_Intel_omp
            printf ">> Making obj_linux_Intel_omp directory...\n"
      fi
      if [ ! -d ../obj_linux_Intel_mpi ]; then
            mkdir ../obj_linux_Intel_mpi
            printf ">> Making obj_linux_Intel_mpi directory...\n"
      fi
#
#   prompt the user for the number of concurrent compile processes to use
#
      printf " \n"
      read -p "Number of concurrent compile processes allowed? (default 1): " JCOMP
      [ -z "$JCOMP" ] && JCOMP=1
#
#   touch main program so it will always be re-compiled (will include compile
#   date & time in warp3d hearder block)
#
touch main_program.f  
#
#   run the makefile for Linux. we now pass more parameters to the makefile
#
      printf "... Starting make program for Linux .... \n\n"
      make -j $JCOMP -f $MAKEFILE COMPILER="$COMPILER" HYPRE_ROOT="$HYPRE_ROOT"

}

#****************************************************************************
#
#     Function:   Global defaults and tests for Mac OS X  - OpenMP only
#
# ****************************************************************************
function mac_main
{
#
MAKEFILE=Makefile.osx
MKLQ=yes
MPIQ=no
HYPQ=no
GFORTRAN=no
INTEL_FORTRAN=no
#
printf "\nSelect build compiler to use:\n"
PS3="Select choice: "
select opt in 'gfortran' 'Intel Fortran' 'Exit'
do
      case $REPLY in
            1 )   printf "\n"
                  GFORTRAN=yes
                  INTEL_FORTRAN=no
                  break
                  ;;
            2 )   printf "\n"
                  GFORTRAN=no
                  INTEL_FORTRAN=yes
                  break
                  ;;
            3 )   exit 0
                  ;;
      esac
done
#
printf ".... Running a series of tests to ensure\n"
printf ".... your system is correctly configured to build WARP3D.\n"
#
# Is this really an OS X system?
#
match=`uname | grep Darwin | wc -l`
if [ $match = "0" ]; then
printf "[ERROR]\n"
printf "This is not a Mac OS X system.\n Quitting...\n\n"
exit 1
fi

check_Intel_Fortran_setup
check_gfortran_setup

#
# The MKL libraries must be present in WARP3D distribution directory
#
if [ ! -d "$WARP3D_HOME/OSX_MKL_files" ]; then
printf "\n[ERROR]\n"
printf "... Directory OSX_MKL_files does not exist in WARP3D\n"
printf "    distribution directory. Run this shell command to\n" 
printf "    download:  install_OSX_libs_from_remote\n  "
printf "Quitting...\n\n"
exit 1
fi
#
# Done with all checks. Looks good for a build process.
}

# ****************************************************************************
#
# Function: Compile WARP3D for MAC OS X
#
# ****************************************************************************
function compile_mac
{
#
printf " \n"
printf ".... This Mac appears configured properly to build WARP3D\n"
printf ".... Compiling WARP3D for Mac OS\n"
printf ".... Installing WARP3D packages for Mac OS..\n"
#
# modify source code to install or unistall WARP3D packages for Mac OS X.
#
printf " \n"
uninstall_mpi
uninstall_hypre
printf " \n"
#
#
# setup the directory structure object and executable if required
#
mkdir ../run_mac_os_x 2> /dev/null
if [ "$INTEL_FORTRAN" = "yes" ]; then
   mkdir ../obj_mac_os_x_intel 2> /dev/null
fi
if [ "$GFORTRAN" = "yes" ]; then
   mkdir ../obj_mac_os_x_gfortran 2> /dev/null
fi
#
# prompt the user for the number of concurrent compile processes to use
#
printf " \n"
read -p "... Number of concurrent compile processes allowed? (default 1): " JCOMP
[ -z "$JCOMP" ] && JCOMP=1
#
touch main_program.f   # so the compile date is always current
#
#
# run the makefile for Mac OS
#
printf "... Starting make program for Mac OS.... \n"
#
if [ "$INTEL_FORTRAN" = "yes" ]; then
  make  BUILD_SYS=Intel -j $JCOMP -f Makefile.osx
fi
if [ "$GFORTRAN" = "yes" ]; then
   make  BUILD_SYS=gfortran -j $JCOMP -f Makefile.osx
fi
#
}
# ****************************************************************************
#
#   main
#
# ****************************************************************************
#
#
#  Check that WARP3D_HOME evvironment variable is set. Change to
#  src deirctory of distribution. Set local shell variables with directory
#  names
#
printf "\n"
printf "** Driver shell script to build WARP3D on Linux and Mac OSX **\n"
#
if [ -z "$WARP3D_HOME" ]; then
   printf "\n\n[ERROR]\n"
   printf "... Environment variable WARP3D_HOME is not set.\n"
   printf "... Usually set in ~/.bashrc (Linux) or ~/.bash_profile (OS X)\n"
   printf "... Quitting...\n\n"
   exit 1
fi
#
cd $WARP3D_HOME/src
mpi_dir=$WARP3D_HOME/linux_packages/source/mpi_code_dir
hypre_dir=$WARP3D_HOME/linux_packages/source/hypre_code_dir
osx_mkl_dir=$WARP3D_HOME/OSX_MKL_files
#
#  Prompt user for platform choice
#
printf "\nSelect supported platform:\n"
PS3="Select choice: "
select opt in 'Linux (64-bit)' 'Mac OS' 'Windows (7-10)' 'Exit'
do
      case $REPLY in
            1 )   printf "\n"
                  linux_main
                  if [ "$INTEL_FORTRAN" = "yes" ]; then
                     compile_linux_Intel
                  fi 
                  break
                  ;;
            2 )   printf "\n"
                  mac_main
                  compile_mac
                  break
                  ;;
            3 )   printf "\n"
                  print_windows_message
                  exit 0
                  ;;
            4 )   exit 0
                  ;;
      esac
done
exit



