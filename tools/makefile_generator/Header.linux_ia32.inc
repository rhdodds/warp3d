#
#               Makefile for WARP3D running on LINUX
#               x86 Processors (32 bit mode)
#               Uses Intel 8.x F90 compiler & MKL 7.x library  
#
#               COMPILER macro value set on command line to run make
#
EXE_NAME = warp3d.ser
ODIR = ././../obj_linux_ia32
MKLIB= $(ODIR)
F90OPTS1 =  -static -fpconstant -fixed  -traceback -module $(ODIR)
F90OPTS2 =  -unroll -ftz -fpstkchk -axN
F90 = $(COMPILER) $(F90OPTS1) $(F90OPTS2)
LINK90 = $(F90)
FILTER = ./filter_linux_ia32.exe
O    = .o
RMC = /bin/rm -f
MVC = /bin/mv -f
#
#                libraries to include at link time. point to static
#                libraries we distribute for MKL. use system pthread 
#
LIBS = $(MKLIB)/libmkl_solver.a \
       $(MKLIB)/libmkl_lapack.a \
       $(MKLIB)/libmkl_ia32.a \
       $(MKLIB)/libguide.a -lpthread
#
#                all object files for warp3d. put all 
#                module files first.
#
OBJ =	 \
[OBJECTS]
	$(ODIR)/[object]$O	 \
[END_OBJECTS]

$(EXE_NAME) : $(OBJ)
	$(LINK90) $(OBJ) $(LOPT) $(LIBS) -o $@
	$(MVC) $@ ./../run_linux_ia32/$@
	chmod ugo+rx  ./../run_linux_ia32/$@

[CUSTOM]
$(ODIR)/[f-file]$O : [f-file].f [dependencies]
	$(RMC)  $@
	$(RMC)  filt_[f-file].f
	$(FILTER) < [f-file].f > filt_[f-file].f
	$(F90) [optimizations]   -c filt_[f-file].f       
	$(MVC) filt_[f-file]$O  [f-file]$O            
	$(MVC) [f-file]$O $@

[END_CUSTOM]
