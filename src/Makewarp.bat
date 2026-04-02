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
::   The same *.nmake file is used with the Makewarp.bat and Makewarp_jom.bat
::   scripts.
::
::   This script builds a 64-bit version of WARP3D.
::   It must be executed from within a command shell that has environment
::   variables initialed to use the 64-bit version of the Intel ifx
::   compiler suite which includes the MKL math library.
::
::   usage:     Makewarp.bat
::
::
:: ****************************************************************************
::
::
::       Location of the "nmake" tool must be set here. 
::
set nmake_exe="c:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.29.30037\bin\Hostx64\x64\nmake.exe"
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
  echo    *             (x86-64 processors)                    *
  echo    *                                                    *
  echo    ******************************************************
  echo.
::
::   ==================================================================
::
::	Create directories for object code and executable files as needed
::
  if not exist ..\run_windows (
     md ..\run_windows
     echo -- Making windows executable file directory...
     echo.
     )
  if not exist ..\obj_windows (
     md ..\obj_windows
     echo -- Making windows object file directory...
     echo.
     )

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
  del /F /Q .\*pdb .\*.exe
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

