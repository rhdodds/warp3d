#
#
#    Makefile for WARP-3D on HP-UX PA8xxx Processors
#    -----------------------------------------------
#
#    Requires: HP-UX 11.x  f90 (release 2.3 or later)
#
#       This is the remainder of the Makefile which is the same for the MPI
#       version as for the serial version.  Do not call this makefile
#       directly from make; call instead Makefile.hpux.mpi or Makefile.hpux.ser
#       to compile the MPI or serial versions of the code.
#
WARPOPTS = +U77 +nofp_exception -v +FPD +Oinfo -R8 +fltconst_strict \
                             +DA2.0W +DS2.0 +DO11.0 +Odataprefetch \
                              +moddir=../obj_hp11 \
                             -I ../obj_hp11 
F90 = /opt/fortran90/bin/f90 $(WARPOPTS)
#
#
FILTER =  filter_hpux11.exe
#
ODIR = ././../obj_hp11
O    = .o
RMC = /bin/rm -f
MVC = /bin/mv -f
#
#                linker options
#
#LOPT = -Wl,-a,archive   <-- no longer use archive
LIBS = $(ODIR)/veclib_pa2_0_64.a $(MPILIB)
#
#                all object files for warp3d
#
OBJ =	 \
[OBJECTS]
	$(ODIR)/[object]$O	 \
[END_OBJECTS]

$(EXE_NAME): $(OBJ)
	$(MPIF90) +O3 +Oparallel  +U77 +nofp_exception -v +FPD +Oinfo -R8 \
                             +DA2.0W +DS2.0 +DO11.0 +Odataprefetch \
	$(OBJ)  $(LOPT) $(LIBS) -o $@
	$(MVC) $@ ./../run_hp11/$@
	chmod ugo+rx  ./../run_hp11/* 

[CUSTOM]
$(ODIR)/[f-file]$O : [f-file].f [dependencies]
	$(RMC)  $@
	$(RMC)  filt_[f-file].f
	$(FILTER) < [f-file].f > filt_[f-file].f
	$(F90) [optimizations]  -c filt_[f-file].f       
	$(MVC) filt_[f-file]$O  [f-file]$O            
	$(MVC) [f-file]$O $@

[END_CUSTOM]