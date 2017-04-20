@echo off
cls
::
::
:: ****************************************************************************
::
::                             Makewarp.bat  cmd
::
::   A command prompt batch file to build WARP3D on Windows.
::   Uses the windows nmake utility and the Makefile.windows.nmake makefile.
::
::   nmake compiles only 1 .f at a time. the alternate build method is in
::   file Makewarp_jom.bat. It uses the free "jom" program as a replacement
::   for nmake. Support for concurrent compiles is available. Greatly speeds
::   up the time to re-compile WARP3D.
::   See comments in Makefile_jom.bat to obtain jom.
::
::   the same *.nmake file is used with the namek and jom programs.
::
::   This script builds a 64-bit version of WARP3D.
::   It must be executed from within a command shell that has environment
::   variables initialed to use the 64-bit version of the Intel F-90
::   compiler suite which now inclues the MKL math library.
::
::   usage:     Makewarp.bat
::
::   This set up works for Intel ifort compiler.
::
:: ****************************************************************************
::
::
::       Location of the "nmake" tool must be set here. The 64-bit command
::       shell does not automatically set up to see Microsoft Visual Studio 8
::       located under the 32-bit c:\Program Files (x86)
::
set nmake_exe="c:\Program Files (x86)\Microsoft Visual Studio 12.0\Intel Fortran\Microsoft Files\VC\Bin\nmake.exe"
::
if not exist %nmake_exe% (
  echo.
  echo .... Fatal Error: the Windows nmake program must be available
  echo ....              set the shell variable inside this script to
  echo ....              its location.
  goto done
  echo.
 )
set build_mode=64
::
::
  echo.
  echo    ******************************************************
  echo    *                                                    *
  echo    *       Makewarp batch file for Windows              *
  echo    *             (64-bit architectures)                 *
  echo    *                                                    *
  echo    ******************************************************
  echo.
::
::   ==================================================================
::
::	Create directories for object code and executable files as needed
::
  if not exist ..\run_windows_%build_mode% (
     md ..\run_windows_%build_mode%
     echo -- Making windows executable file directory...
     echo.
     )
  if not exist ..\obj_windows_%build_mode% (
     md ..\obj_windows_%build_mode%
     echo -- Making windows object file directory...
     echo.
     )
::
::   ==================================================================
::
::	Install dummy MPI code for Windows -- we do not yet
::      support MPI and WARP3D on Windows.
::
  echo -- Uninstalling MPI code...
  del .\mpi_code.f 
  del .\mpi_handle_slaves.f
  del .\mod_local_stiffness.f
  del .\distributed_assembly.f
::
  set xdir=..\linux_packages\source\mpi_code_dir
::
  copy %xdir%\mpi_code_dummy.f .\mpi_code.f /y
  echo        - mpi_code_dummy.f
  copy %xdir%\mpi_handle_slaves_dummy.f .\mpi_handle_slaves.f /y
  echo        - mpi_handle_slaves.f
  copy %xdir%\mod_local_stiffness_dummy.f .\mod_local_stiffness.f /y
  echo        - mod_local_stiffness.f
  copy %xdir%\distributed_assembly_dummy.f .\distributed_assembly.f /y
  echo        - distributed_assembly.f
: 
:
::
  echo -- MPI code uninstalled...
  echo.
::
::   ==================================================================
::
::	Install dummy hypre code for Windows -- we do not yet
::      support hypre and WARP3D on Windows.
::
  echo -- Uninstalling hypre code...
  if  exist .\hypre_driver.f (
      del .\hypre_driver.f
     )
  if  exist  .\hypre_warp_only.f (
      del  .\hypre_warp_only.f
     )
  if  exist  .\iterative_sparse_hypre.f (
      del  .\iterative_sparse_hypre.f
     )
::
  set xdir=..\linux_packages\source\hypre_code_dir
::
  copy %xdir%\dummy_sparse_hypre.f .\iterative_sparse_hypre.f /y
  echo        - dummy_sparse_hypre.f
::
  echo -- hypre code uninstalled...
  echo.
::
::   ==================================================================
::
::	Run the makefile.
::
::  touch main so it always gets compiled. internally gets date/time
::  of this compile and includes in warp3d header block
::
copy /b main_program.f +,,
::
if "%build_mode%" == "64" (
  echo -- Compiling WARP3D for 64-bit execution on
  echo -- Windows with nmake utility...
  echo.
  %nmake_exe% /f Makefile.windows.nmake ARCH64=64
 )
  del .\*.exe
::
  echo.
  echo.
  echo -- WARP3D updated. You may need to increase the
  echo    virtual memory limits on Windows implementations.
  echo    Refer to the README file in the main WARP3D directory...
  echo.
  echo -- All done...
  echo.
:done

