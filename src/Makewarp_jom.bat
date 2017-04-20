@echo off
cls
::
::
:: ****************************************************************************
::
::                             Makewarp_jom.bat cmd file
::   
::   Usage:   Makewarp_jom.bat exectured in a Command shell
::
::
::   This version requires jom.exe to support concurrent compiles using
::   the standard .nmake file for building WARP3D.
::
::   jom is a replacement for nmake
::
::   The WARP3D src/Makefile.windows.nmake works with both the usual 
::   "nmake" command and the "jom" replacement for nmake.
:: 
::   http://wiki.qt.io/Jom
::
::   This script builds a 64-bit version of WARP3D.
::   It must be executed from within a command shell that has environment
::   variables initialed to use the 64-bit version of the Intel F-90
::   compiler suite which now inclues the MKL math library.
::
::   This set up works for Intel ifort compiler.
::
:: ****************************************************************************
::
::
::       Location of the "jom" tool must be set here.
::       Number of conccurrent compiles to be allowed
::
set jom_exe="c:\Users\rdodds\bin\jom_1_1_2\jom.exe"
set jobs="7"
::
if not exist %jom_exe% (
  echo.
  echo .... Fatal Error: the jom.exe tool must be available
  echo ....              set the shell variable inside this script to
  echo ....              its location.
  echo ....              http://wiki.qt.io/Jom
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
  echo -- Windows with jom utility...
  echo.
  %jom_exe% /f Makefile.windows.nmake ARCH64=64 /J %jobs%
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

