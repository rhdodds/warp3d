#
# ****************************************************************************
#                                                               
#   File: Header.sga_linux.inc            
#                                                               
#   Programmer: agcarlyle
#   Date: Aug 2005                                                           
#
#   Description:
#     Makefile Header for use with the WARP3D Makefile Generator Program.
#     This header will make a Makefile for SGI Altix Supercomputers with
#     SGI ProPack Linux and Intel Itanium2+ Processors (64 bit mode).
#     Uses Intel 8.x F90 compiler & MKL 7.x library. The "COMPILER" macro
#     value is set in "Makewarp.ksh"
#
#   Revision History:
#
# *****************************************************************************
#
#
#
ODIR     = ././../obj_sga_linux
MKLIB    = $(ODIR)
F90OPTS1 =  -tpp2 -fpconstant -fixed -72 -traceback
F90OPTS2 =  -module $(ODIR) -ftz
F90      = $(COMPILER) $(F90OPTS1) $(F90OPTS2) $(MPIOPT)
LINK90   = $(COMPILER) $(F90OPTS1) $(F90OPTS2) $(MPIOPT)
FILTER   = ./filter_sga_linux.exe
O        = .o
RMC      = /bin/rm -f
MVC      = /bin/mv -f
#
#                Libraries to include at link time. Point to static
#                libraries we distribute for MKL. Use system pthread. 
#
LIBS = $(MKLIB)/libmkl_solver.a \
       $(MKLIB)/libmkl_lapack.a \
       $(MKLIB)/libmkl_ipf.a \
       $(MKLIB)/libguide.a -lpthread
# 
#                Program will insert all object files for warp3d. All 
#                module files are listed first.
#
OBJ =	 \
[OBJECTS]
	$(ODIR)/[object]$O	 \
[END_OBJECTS]

$(EXE_NAME) : $(OBJ)
	$(LINK90) $(OBJ) $(LIBS) -o $@
	$(MVC) $@ ./../run_sga_linux/$@
	chmod ugo+rx  ./../run_sga_linux/$@

[CUSTOM]
$(ODIR)/[f-file]$O : [f-file].f [dependencies]
	$(RMC)  $@
	$(RMC)  filt_[f-file].f
	$(FILTER) < [f-file].f > filt_[f-file].f
	$(F90) [optimizations]  -c filt_[f-file].f       
	$(MVC) filt_[f-file]$O  [f-file]$O            
	$(MVC) [f-file]$O $@

[END_CUSTOM]
