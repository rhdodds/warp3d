@echo off
:: Batch file to compile all the makefiles for windows
echo Compiling mesh generators on Windows...
echo.

:: Compiler and nmake
set compiler=ifort
set make="c:\Program Files (x86)\Microsoft Visual Studio 9.0\VC\Bin\nmake.exe"

:: Suffix for windows executables
set suffix=.exe

:: Start with the easy ones
for %%F in (8to20 add_elm mesh_plot mesh_ssy mesh_ssy2_ts) do (
	echo Compiling %%F
	cd %%F
	%compiler% *.f /o %%F%suffix%
	cd ..
	echo.
)

:: Now for the ones with makefiles
echo Compiling pipe_mesh_gen
cd pipe_mesh_gen
%make% /f make_windows f77=%compiler% f90=%compiler%
cd ..
echo.

echo Compiling mesh_scp
cd mesh_scp
%make% /f makefile_windows F90=%compiler%
cd ..
echo.

echo Compiling mesh_cell
cd mesh_cell
%make% /f makefile_windows F77=%compiler%
cd ..
echo.

:: clean up
call clean_windows.bat
echo.
echo Done...
echo.