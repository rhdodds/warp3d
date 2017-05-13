#!/bin/bash

# *****************************************************************************
#
#     Makewarp.bash (third version)
#
#     modified: Dec 2016 RHD
#
#     Description:
#           Bash script to interactively drive compilation of Linux and Mac
#           versions of WARP3D.  For Windows, prints message and quits,
#           for Mac, does not allow customization, just hunts down the MKL
#           files, checks intel compilers, and compiles openmp version.  For
#           Linux can either do that (simple mode) or prompt interactively for
#           the libraries, compiler, etc.
#
# *****************************************************************************


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

  printf "   The Intel Visual Fortran Compiler creates a link to the proper \n"
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
                        MAKEFILE=Makefile.linux_em64t.mpi_omp
                        MPIQ=yes
                        break
                        ;;
                  "OpenMP")
                        printf "\n"
                        COMPILER=ifort
                        ALTCOMPILER=ifort
                        MAKEFILE=Makefile.linux_em64t.omp
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
#   Intel Fortran Composer XE system must be installed  
#
      hash ifort 2>&- || {
            printf "[ERROR]\n"
            printf "[ERROR]\n"
            printf "... Cannot find the Intel Fortran compiler (ifort) in your PATH.\n"
            printf "... See Intel Fortran install documentation. Most often a line of\n"
            printf "... the form: source /opt/intel/composerxe/... is placed\n"
            printf "... in the /etc/bash.bashrc file on your system or in your \n"
            printf "... ~/.bashrc file\n"
            printf "Quitting...\n\n"
            exit 1
      }

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
#     Function:   Compile WARP3D for linux
#
# ****************************************************************************

function compile_linux {
	
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
      if [ ! -d ../run_linux_em64t ]; then
            mkdir ../run_linux_em64t
            printf ">> Making run_linux_em64t directory...\n"
      fi
      if [ ! -d ../obj_linux_em64t ]; then
            mkdir ../obj_linux_em64t
            printf ">> Making obj_linux_em64t directory...\n"
      fi
      if [ ! -d ../obj_linux_em64t_mpi ]; then
            mkdir ../obj_linux_em64t_mpi
            printf ">> Making obj_linux_em64t_mpi directory...\n"
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
#
# Intel Fortran Composer XE system must be installed 
#
hash ifort 2>&- || {
printf "[ERROR]\n"
printf "... Cannot find the Intel Fortran compiler (ifort) in your PATH.\n"
printf "... See Intel Fortran install documentation. Most often a line of\n"
printf "... the form: source /opt/intel/composerxe/... is placed\n"
printf "... in the /etc/bashrc file on your Mac.\n"
printf "Quitting...\n\n"
exit 1
}
#
#
# The MKLROOT environment variable must be set. This is normally
# done by the user's login script or system shell initialization
# script, e.g., source /opt/intel/composerxe/bin/compilervars.sh intel64
# Often best to put this in /etc/bashrc
#
if [ -z "$MKLROOT" ]; then
printf "[ERROR]\n"
printf "... Environment variable MKLROOT is not set.\n"
printf "... See Intel Fortran install documentation. Most often a line of\n"
printf "... the form: source /opt/intel/composerxe/... is placed\n"
printf "... in the /etc/bashrc file on your Mac.\n"
printf "Quitting...\n"
exit 1
fi
#
# The MKL libraries must be accesible.
#
if [ ! -d "$MKLROOT/lib" ]; then
printf "[ERROR]\n"
printf "... Directory MKLROOT/lib does not exist.\n"
printf "Quitting...\n\n"
exit 1
fi
#
# Fortran libraries must be accessible for our makefile to function properly.
# We have to use a small kludge by locating them from the MKL directory. Intel
# does not provide an environment variable to locate them.
#
if [ ! -d $MKLROOT/../compiler/lib/intel64 ]; then
printf "[ERROR]\n"
printf "... Could not find the Intel Fortran library directory.\n"
printf "... by examining $MKLROOT/../compiler/lib/intel64\n"
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
if [ ! -d ../run_mac_os_x ]; then
mkdir ../run_mac_os_x
printf ">> Making run_mac_os_x directory...\n"
fi
if [ ! -d ../obj_mac_os_x ]; then
mkdir ../obj_mac_os_x
printf ">> Making obj_mac_os_x directory...\n"
fi
#
# prompt the user for the number of concurrent compile processes to use
#
printf " \n"
read -p "... Number of concurrent compile processes allowed? (default 1): " JCOMP
[ -z "$JCOMP" ] && JCOMP=1
#
touch main_program.f   # so thecompile date is always current
#
#
# run the makefile for Mac OS
#
printf "... Starting make program for Mac OS.... \n"
make  -j $JCOMP -f Makefile.osx
}

# ****************************************************************************
#
#     Function:   main
#
# ****************************************************************************
#
#
#  Check that WARP3D_HOME evvironment variable is set. Change to
#  src deirctory of distribution. Set local shell variables with directory
#  names
#
printf "\n"
printf "** Driver shell script to build WARP3D on Linux and Mac OS **\n"
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
mkl_dir=$WARP3D_HOME/mkl_solver_dir
mpi_dir=$WARP3D_HOME/linux_packages/source/mpi_code_dir
hypre_dir=$WARP3D_HOME/linux_packages/source/hypre_code_dir
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
                  compile_linux
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


