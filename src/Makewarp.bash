#!/bin/bash
#
#     Makewarp.bash (9 th version)
#
#     modified: April 2, 2026 rhd
#               ** remove MPI support **
#               ** use ifx compiler on Windows, Linux, Apple silicon ** 
#               ** add gfortran option for Intel mac, gfortran mac silicon **
#
#     Description:
#
#           Bash script to interactively drive compilation of Linux and macOS
#           versions of WARP3D.  For Windows, we print message to use another
#           script and quit.
#
#           Run this script in a Bash shell:  Makewarp.bash
#
#           Mac and Linux: no user selectable options.
#                script hunts down the MKL files, checks for Intel (ifort) compiler, 
#                and runs the makefile
#                to build the threads (OpenMP) executable.
#
#      Main program (function) at bottom of this script
#
# ****************************************************************************
#
#   Function: Check macOS Intel ifort Fortran compiler exists & version
#
# ****************************************************************************
#
function check_macOS_ifort
{
hash ifort 2>&- || {
printf "[ERROR]\n"
printf "... Cannot find the Intel Fortran compiler (ifort) in your PATH.\n"
printf "Quitting...\n\n"
exit 1
}
#
/bin/rm zqq03 >& /dev/null
ifort --version -diag-disable=10448  >zqq03
sed -i -e '2,10d' zqq03
echo -e "\n ... Intel Fortran (ifort) detected:" `cat zqq03`
count2=`grep "2021." zqq03 |wc -l`
count3=`grep "2022." zqq03 |wc -l`
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
    printf "\n... 2021..1 or newer"
    printf "\n... other versions have known bugs that affect WARP3D"
    printf "\n... Quitting...\n\n"
    exit 1
fi
#
return
}
# ****************************************************************************
#
#   Function: Check Linux Intel ifx Fortran compiler exists & version
#
# ****************************************************************************
#
function check_Linux
{
#
hash ifx 2>&- || {
printf "[ERROR]\n"
printf "... Cannot find the Intel Fortran compiler (ifx) in your PATH.\n"
printf "... See Intel Fortran install documentation. Most often a line of\n"
printf "... the form: source /opt/intel/... is placed\n"
printf "... in the /etc/bashrc our your ~/.bashrc file.\n"
printf "Quitting...\n\n"
exit 1
}
/bin/rm zqq03 >& /dev/null
ifx --version  >zqq03
sed -i -e '2,$d' zqq03
echo -e "\n... Intel Fortran (ifx) detected:" `cat zqq03`
count2=`grep "2025" zqq03 |wc -l`
/bin/rm zqq03*
ok=0
if [ $count2 -eq "1" ]; then
   ok=1
fi
if [ $ok -eq "0" ]; then
    printf "\n... ERROR: ifx must be one of these versions:"
    printf "\n... 2025 or newer"
    printf "\n... Quitting...\n\n"
    exit 1
fi
#
return
}
# ****************************************************************************
#
#   Function: Issue message to use Makewarp.bat to compile on Windows
#
# ****************************************************************************
#
function print_windows_message
{
  printf "\n>> To compile on Windows, exit this \n"
  printf "   script and run the 'Makewarp.bat' batch file from within a \n"
  printf "   Windows command prompt shell. The shell must be setup \n"
  printf "   64-bit building with the Intel ifx compiler.\n\n"

  printf "   The Intel ifx Fortran Compiler creates a link to the proper \n"
  printf "   command prompt build environment in the Windows 'start' menu. \n\n"
}

# ****************************************************************************
#
#     Function:   Global defaults and branch point for Linux
#
# ****************************************************************************
#
function linux_main
{
#
# These are common defaults
#
INTEL_FORTRAN=no
GFORTRAN=no
MKLQ=yes
#
 COMPILER=ifx    # default
 ALTCOMPILER=ifx
 MAKEFILE=Makefile.linux
 MPIQ=no
 #
check_Linux
linux_simple
# 
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
      if [ "$(uname -s)" != "Linux" ]; then 
            printf "[ERROR]\n"
            printf "This is not a Linux system.\n Quitting...\n\n"
            exit 1
      fi
      HYPQ=no
#
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

# Setup the directory structure if required
#
      if [ ! -d ../run_linux ]; then
            mkdir ../run_linux
            printf ">> Making run_linux directory...\n"
      fi
      if [ ! -d ../obj_linux ]; then
            mkdir ../obj_linux
            printf ">> Making obj_linux directory...\n"
      fi
#
#   prompt the user for the number of concurrent compile processes to use
#
      printf " \n"
      read -p "Number of concurrent compile processes allowed? (default 1): " JCOMP
      [ -z "$JCOMP" ] && JCOMP=1
case "$JCOMP" in
    ''|*[!0-9]*)
        printf "\n[ERROR]\n"
        printf "... JCOMP must be an integer >= 1.\n\n" >&2
        exit 1
        ;;
esac

if [ "$JCOMP" -lt 1 ]; then
    printf "\n[ERROR]\n"
    printf "... JCOMP must be >= 1.\n\n" >&2
    exit 1
fi
##
#   touch main program so it will always be re-compiled (will include compile
#   date & time in warp3d hearder block)
#
touch main_program.f  
#
#   run the makefile for Linux. we now pass more parameters to the makefile
#
      printf "... Starting make program for Linux .... \n\n"
      make  -j $JCOMP -f $MAKEFILE 

}

#****************************************************************************
#
#     Function:   Build for macOS (Intel,ifort)
#
# ****************************************************************************
function compile_mac_Intel_ifort
{
#
printf ".... Running a series of tests to ensure\n"
printf ".... your system is correctly configured to build WARP3D.\n"
printf ".... Intel processors, ifort comoiler, MKL + Pardsio solver.\n"
#
# Is this really a macOS system?
#
if [ "$(uname -s)" != "Darwin" ]; then
   printf "[ERROR]\n"
   printf "This is not a macOS system.\n Quitting...\n\n"
   exit 1
fi
#
check_macOS_ifort
#
if [ ! -f "$WARP3D_HOME/src/Makefile.osx" ]; then
   printf "\n[ERROR]\n"
   printf "... File Makefile.osx does not exist in WARP3D/src\n"
   printf "    distribution directory. "
   printf "Quitting...\n\n"
   exit 1
fi
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
# Are the Xcode Command Line Tools installed
#
if pkgutil --pkg-info=com.apple.pkg.CLTools_Executables >/dev/null 2>&1; then
    echo ""    #  Xcode Command Line Tools installed"
else
    printf "\n.... The Xcode Command Line Tools are not installed\n"
    printf "Quitting...\n\n"
    exit 1
fi
#
# Done with all checks. Looks good for a build process.
#
printf " \n"
printf ".... This Mac appears configured properly to build WARP3D\n"
printf ".... Compiling WARP3D for macOS Intel with ifort\n"
#
# setup the directory structure object and executable if required
#
mkdir ../run_macOS 2> /dev/null
mkdir ../obj_macOS_ifort 2> /dev/null
#
# prompt the user for the number of concurrent compile processes to use
#
printf " \n"
read -p "... Number of concurrent compile processes allowed? (default 1): " JCOMP
[ -z "$JCOMP" ] && JCOMP=1
case "$JCOMP" in
    ''|*[!0-9]*)
        printf "\n[ERROR]\n"
        printf "... JCOMP must be an integer >= 1.\n\n" >&2
        exit 1
        ;;
esac

if [ "$JCOMP" -lt 1 ]; then
    printf "\n[ERROR]\n"
    printf "... JCOMP must be >= 1.\n\n" >&2
    exit 1
fi
#

touch main_program.f   # so the compile date is always current
#
#
# run the makefile for macOS
#
printf "... Starting make program for macOS.... \n"
#
printf "... Note: ignore Linker messages: ipo: warning #11109: unable to ..."
printf "\n\n"
make -j $JCOMP -f Makefile.osx PLATFORM=Intel_ifort
}
#
#****************************************************************************
#
#     Function:   Build for macOS  (Intel,gfortran)
#
# ****************************************************************************
function compile_mac_Intel_gfortran
{
#
# Is this really an maxOS system?
#
match=`uname | grep Darwin | wc -l`
if [ $match = "0" ]; then
   printf "[ERROR]\n"
   printf "This is not a Mac OS X system.\n Quitting...\n\n"
   exit 1
fi
#
MAKEFILE=Makefile.osx
#
if [ ! -f "$WARP3D_HOME/src/Makefile.osx" ]; then
   printf "\n[ERROR]\n"
   printf "... File Makefile.osx does not exist in WARP3D/src\n"
   printf "    distribution directory. "
   printf "Quitting...\n\n"
   exit 1
fi
#
# Is gfortran installed?
#
hash gfortran 2>&- || {
printf "[ERROR]\n"
printf "... Cannot find the gfortran compiler in your PATH.\n"
printf "Quitting...\n\n"
exit 1
}
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
# The MUMPS library for Intel mac must be present in WARP3D distribution directory
#
if [ ! -f "$WARP3D_HOME/MUMPS_libs/libmumps_all_gfort_osx_Intel_metis.a" ]; then
printf "\n[ERROR]\n"
printf "... The file MUMPS_libs/libmumps_all_gfort_osx_Intel_metis.a\n"
printf "     does not exist in the WARP3D distribution directory. Run this shell command to\n" 
printf "    download:  install_MUMPS_libs_from_remote\n  "
printf "Quitting...\n\n"
exit 1
fi
#
# Are the Xcode Command Line Tools installed
#
if pkgutil --pkg-info=com.apple.pkg.CLTools_Executables >/dev/null 2>&1; then
    echo "" #Xcode Command Line Tools installed"
else
    printf "\n.... The Xcode Command Line Tools are not installed\n"
    printf "Quitting...\n\n"
    exit 1
fi
#
# setup the directory structure object and executable if required
#
mkdir ../run_macOS 2> /dev/null
mkdir ../obj_macOS_Intel_gfortran 2> /dev/null
#
# prompt the user for the number of concurrent compile processes to use
#
printf " \n"
read -p "... Number of concurrent compile processes allowed? (default 1): " JCOMP
[ -z "$JCOMP" ] && JCOMP=1
case "$JCOMP" in
    ''|*[!0-9]*)
        printf "\n[ERROR]\n"
        printf "... JCOMP must be an integer >= 1.\n\n" >&2
        exit 1
        ;;
esac

if [ "$JCOMP" -lt 1 ]; then
    printf "\n[ERROR]\n"
    printf "... JCOMP must be >= 1.\n\n" >&2
    exit 1
fi
##
touch main_program.f   # so the compile date is always current
#
# run the makefile for macOS Intel gfortran
#
printf "... Starting make program for macOS.... \n"
#
printf "\n\n"
make -j $JCOMP -f Makefile.osx PLATFORM=Intel_gfortran
#
# Done with all checks. Looks good for a build process.
}
#
#****************************************************************************
#
#     Function:   Build for macOS  (Apple Silicon,gfortran)
#
# ****************************************************************************
function compile_mac_silicon
{
#
# Is this really an maxOS system?
#
if [ "$(uname -s)" != "Darwin" ]; then
   printf "[ERROR]\n"
   printf "This is not a Mac OS X system.\n Quitting...\n\n"
   exit 1
fi
#
arch=$(uname -m)
[[ "$arch" == "arm64" ]] || error_exit "This option requires Apple Silicon (arm64)."
#
MAKEFILE=Makefile.osx
#
if [ ! -f "$WARP3D_HOME/src/Makefile.osx" ]; then
   printf "\n[ERROR]\n"
   printf "... File Makefile.osx does not exist in WARP3D/src\n"
   printf "    distribution directory. "
   printf "Quitting...\n\n"
   exit 1
fi
#
# Is gfortran installed?
#
hash gfortran 2>&- || {
printf "[ERROR]\n"
printf "... Cannot find the gfortran compiler in your PATH.\n"
printf "Quitting...\n\n"
exit 1
}
#
# The MUMPS library for Intel mac must be present in WARP3D distribution directory
#
if [ ! -f "$WARP3D_HOME/MUMPS_libs/libmumps_all_gfort_osx_arm64_metis.a" ]; then
printf "\n[ERROR]\n"
printf "... The file MUMPS_libs/libmumps_all_gfort_osx_arm64_metis.a\n"
printf "     does not exist in the WARP3D distribution directory. Run this shell command to\n" 
printf "    download:  install_MUMPS_libs_from_remote\n  "
printf "Quitting...\n\n"
exit 1
fi
#
# Are the Xcode Command Line Tools installed
#
if pkgutil --pkg-info=com.apple.pkg.CLTools_Executables >/dev/null 2>&1; then
    echo "" #Xcode Command Line Tools installed"
else
    printf "\n.... The Xcode Command Line Tools are not installed\n"
    printf "Quitting...\n\n"
    exit 1
fi
#
# setup the directory structure object and executable if required
#
mkdir ../run_macOS 2> /dev/null
mkdir ../obj_macOS_silicon 2> /dev/null
#
# prompt the user for the number of concurrent compile processes to use
#
printf " \n"
read -p "... Number of concurrent compile processes allowed? (default 1): " JCOMP
[ -z "$JCOMP" ] && JCOMP=1
case "$JCOMP" in
    ''|*[!0-9]*)
        printf "\n[ERROR]\n"
        printf "... JCOMP must be an integer >= 1.\n\n" >&2
        exit 1
        ;;
esac

if [ "$JCOMP" -lt 1 ]; then
    printf "\n[ERROR]\n"
    printf "... JCOMP must be >= 1.\n\n" >&2
    exit 1
fi
##
touch main_program.f   # so the compile date is always current
#
# run the makefile for macOS gfortran
#
printf "... Starting make program for macOS Silicon.... \n"
#
printf "\n\n"
make -j $JCOMP -f Makefile.osx PLATFORM=silicon
#
# Done with all checks. Looks good for a build process.
}
# ****************************************************************************
#
#   main
#
# ****************************************************************************
#
#
#  Check that WARP3D_HOME environment variable is set. Change to
#  src deirctory of distribution. Set local shell variables with directory
#  names
#
printf "\n"
printf "** Driver shell script to build WARP3D on Linux and macOS **\n"
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
osx_mkl_dir=$WARP3D_HOME/OSX_MKL_files
MUMPS_libs=$WARP3D_HOME/MUMPS_libs
#
#  Prompt user for platform choice
#
printf "\nSelect supported platform:\n"
PS3="Select choice: "
select opt in 'Linux (Intel,ifx) ' 'macOS (Intel,ifort)'  \
           'macOS (Silicon)' 'Windows (10,11)' 'macOS (Intel,gfortran)' 'Exit'
do
      case $REPLY in
            1 )   printf "\n"
                  linux_main
                  compile_linux_Intel
                  break
                  ;;
            2 )   printf "\n"
                  compile_mac_Intel_ifort
                  break
                  ;;
            3 )   printf "\n"
                  compile_mac_silicon
                  break
                  ;;
            4 )   printf "\n"
                  print_windows_message
                  exit 0
                  ;;
            5 )   printf "\n"
                  compile_mac_Intel_gfortran
                  break
                  ;;
            6 )   exit 0
                  ;;
      esac
done
exit



