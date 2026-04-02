@echo off
cls
::
::
:: ****************************************************************************
::
::                compile_patwarp_windows.bat
::
::   Build executable for patwarp utility program on Windows (64-bit mode)
::   and put it into the run_windows_64 directory
::
::   This set up works for Intel ifx compiler
::
:: ****************************************************************************
::
::
::
  echo.
  echo    ******************************************************
  echo    *                                                    *
  echo    *   Build executable for patwarp utility on Windows  *
  echo    *                                                    *
  echo    ******************************************************
  echo.   
::
::   ==================================================================
::
::
set fopts1=/O2 /fpconstant /fixed /traceback /QaxAVX /static
set fopts2=/F4000000
set files=patwarp.f
ifx %fopts1% %fopts2% %files% /o patwarp.exe
::
  echo.
  echo ... Build done
  echo ... Cleanup and move .exe to run_windows directory
  echo.
::
del .\patwarp_data.mod
del .\patwarp.obj
del .\patwarp.pdb
copy patwarp.exe ..\run_windows\patwarp_windows.exe
del patwarp.exe
  echo.
  echo ... All Done ....
  echo.
:done

