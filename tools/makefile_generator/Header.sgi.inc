#
#
#               Makefile for WARP-3D on SGI
#               ---------------------------
#
#
#
#         this makefile works with the MIPSpro 7.2 compilers
#         for earlier versions of the f90 compiler, make
#         all optimization levels = 0 for the f90 compiler.
#
#       This is the remainder of the Makefile which is the same for the MPI
#       version as for the serial version.  Do not call this makefile
#       directly from make; call instead Makefile.sgi.mpi or Makefile.sgi.ser
#       to compile the MPI or serial versions of the code.
#
F90 = f90 -64 
#
FILTER = ./filter_sgi.exe
#
ODIR = ././../obj_sgi
O    = .o
RMC = /bin/rm -f
MVC = /bin/mv -f
#
LOPT = -mp 
#
LIBS = -lscs_mp
#
OBJ =	 \
[OBJECTS]
	$(ODIR)/[object]$O	 \
[END_OBJECTS]


$(EXE_NAME) : $(OBJ)
	$(F90) $(LOPT) $(OBJ) $(LIBS) $(MPI) -o $@
	$(MVC) $@ ./../run_sgi/$@
	chmod ugo+rx ./../run_sgi/$@

[CUSTOM]
$(ODIR)/[f-file]$O :  [f-file].f [dependencies]
	$(RMC) $@
	$(RMC) filt_[f-file].f
	$(FILTER) < [f-file].f > filt_[f-file].f
	$(F90) [optimizations]   -c filt_[f-file].f       
	$(MVC) filt_[f-file]$O   [f-file]$O            
	$(MVC) [f-file]$O $@

[END_CUSTOM]
