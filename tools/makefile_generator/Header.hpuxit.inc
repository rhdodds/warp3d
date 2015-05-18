#
#
#    Makefile for WARP-3D on HP-UX Itanium 2+ Processors
#    ---------------------------------------------------
#
#    Requires: HP-UX 11.x  f90 (release 2.3 or later)
#
#       This is the remainder of the Makefile which is the same for the MPI
#       version as for the serial version.  Do not call this makefile
#       directly from make. Use the driver shell script Makewarp.ksh !
#
WARPOPTS = +U77 +nofp_exception -v +FPD +Oinfo -R8 +fltconst_strict \
           +DD64 +DSitanium2 +moddir=../obj_hpit -I ../obj_hpit
F90 = /opt/fortran90/bin/f90 $(WARPOPTS)
#
#
FILTER =  filter_hpuxit.exe
#
ODIR = ././../obj_hpit
O    = .o
RMC = /bin/rm -f
MVC = /bin/mv -f
#
#                linker options (must include mlib)
#
LIBS = $(ODIR)/veclib_itan_64.a $(MPILIB)
#
#                all object files for warp3d
#
OBJ =	 \
[OBJECTS]
	$(ODIR)/[object]$O	 \
[END_OBJECTS]

$(EXE_NAME): $(OBJ)
	$(MPIF90)   \
		+O3 +Oparallel  +U77 +nofp_exception -v +FPD +Oinfo -R8 \
                             +DSitanium2 \
	$(OBJ)  $(LOPT) $(LIBS) -o $@
	$(MVC) $@ ./../run_hpit/$@
	chmod ugo+rx  ./../run_hpit/* 

[CUSTOM]
$(ODIR)/[f-file]$O : [f-file].f [dependencies]
	$(RMC)  $@
	$(RMC)  filt_[f-file].f
	$(FILTER) < [f-file].f > filt_[f-file].f
	$(F90) [optimizations]   -c filt_[f-file].f       
	$(MVC) filt_[f-file]$O  [f-file]$O            
	$(MVC) [f-file]$O $@

[END_CUSTOM]
