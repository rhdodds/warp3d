@echo off
cls
::
::
:: ****************************************************************************
::
::                             Makewarp.windows.cmd
::
::   A command prompt batch file to build WARP3D on Windows XP, Vista, Windows 7.
::   Uses the windows nMake utility and the Makefile.windows.nmake makefile.
::
::   This script builds a 64-bit version of WARP3D.
::   It must be executed from within a command shell that has environment
::   variables initialed to use the 64-bit version of the Intel F-90
::   compiler suite which now inclues the MKL math library.
::
::   usage:     Makewarp.bat
::
::   This set up works for Intel ifort compiler version 13.0 and later.
::
:: ****************************************************************************
::
::
::
::
::       Location of the "nmake" tool must be set here. The 64-bit command
::       shell does not automatically set up to see Microsoft Visual Studio 8
::       located under the 32-bit c:\Program Files (x86)
::
set nmake_exe="c:\Program Files (x86)\Microsoft Visual Studio 10.0\Intel Fortran\Microsoft Files\VC\Bin\nmake.exe"
#
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
  echo    *   Makewarp batch file for WinXP/Vista/Windows 7    *
  echo    *             (64-bit architectures)                 *
  echo    *                                                    *
  echo    ******************************************************
  echo.
::
::   ==================================================================
::
::	 Build the usual filter program if required. WARP3D source code
::       files pass thru the filter to make the OS specific version.
::
  if not exist .\filter_intel_windows.exe (
     echo -- Building filter program...
     echo.
     ifort /O2 /nologo /nowarn filter_intel_windows.f /o filter_intel_windows.exe
     if exist .\filter_intel_windows.obj (del .\filter_intel_windows.obj)
     )
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
::
  set xdir=..\linux_packages\source\mpi_code_dir
::
  copy %xdir%\mpi_code_dummy.f .\mpi_code.f /y
  echo        - mpi_code_dummy.f
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
if "%build_mode%" == "64" (
  echo -- Compiling WARP3D for 64-bit execution on
  echo -- Windows XP, Vista, Windows 7 with nmake utility...
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

