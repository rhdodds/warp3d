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
::   This set up works for Intel F-90 version 11.x.
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
set fopts1=/O2 /fpconstant /fixed /traceback /QaxSSE4.2 /heap-arrays /static
set fopts2=/F4000000
set files=patwarp_module.f patwarp_windows_mac.f
ifort %fopts1% %fopts2% %files% /o patwarp.exe
::
  echo.
  echo ... Build done
  echo ... Cleanup and move .exe to run_windows_64 directory
  echo.
::
del .\patwarp_data.mod
del .\patwarp_module.obj
del .\patwarp_windows_mac.obj
copy patwarp.exe ..\run_windows_64\patwarp_windows.exe
del patwarp.exe
  echo.
  echo ... All Done ....
  echo.
:done

