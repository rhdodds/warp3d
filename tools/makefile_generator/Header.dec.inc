# 
#
#               Makefile for WARP3D on DEC-Alpha
#               --------------------------------
#
#                 requires only f90
#
F90 = f90 -V -fast -fpe3 -I../obj_dec -module ../obj_dec
FILTER = filter_dec.exe
#
ODIR = ././../obj_dec
O = .o
RMC = /bin/rm -f
MVC = /bin/mv -f
#
#                linker options
#
LOPT = 
#
#                libraries to include at link time
#
LIBS = $(ODIR)/libsolver.a $(MPILIB)
#
#                all object files for warp3d
#
OBJ =	 \
[OBJECTS]
	$(ODIR)/[object]$O	 \
[END_OBJECTS]

$(EXE_NAME): $(OBJ)
	$(F90) $(LOPT) $(OBJ) $(LIBS) -ncr -omp -o $@
	$(MVC) $@ ./../run_dec/$@
	chmod ugo+rx  ./../run_dec/$@

[CUSTOM]
$(ODIR)/[f-file]$O :  [f-file].f  [dependencies]
	$(RMC) $@
	$(RMC) filt_[f-file].f
	$(FILTER) < [f-file].f > filt_[f-file].f
	$(F90) [optimizations]   -c filt_[f-file].f       
	$(MVC) filt_[f-file]$O       [f-file]$O            
	$(MVC) [f-file]$O $@

[END_CUSTOM]
