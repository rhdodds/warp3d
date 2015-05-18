#
#
#               Makefile for WARP-3D on RISC-6000, SP-x
#               ---------------------------------------
#
#       requires f90 (xlf_r)
#
#	This is the remainder of the Makefile which is the same for the MPI
#	version as for the serial version.  Do not call this makefile
#	directly from make; call instead Makefile.ibm.mpi or Makefile.ibm.ser
#	to compile the MPI or serial versions of the code.
#
#
#       use compiler which links in re-entrant versions of runtime
#       libraries: only f90 requried now
#
#       we build in 64-bit mode only
#
F90 = xlf_r -q64 -qwarn64 -b64 -qarch=auto -qtune=auto -qcache=auto -qnoipa -qreport
FILTER = ./filter_ibm.exe
C = xlc
#
ODIR = ././../obj_ibm
O = .o
RMC = /bin/rm -f
MVC = /bin/mv -f
#
#                linker options - large virtual addresses (>2GB) and include
#                                 threads for parallel and the math library
#                                 maxdata value here is 32768 MB
#
LOPT = -bmaxdata:0x800000000 -qstrict -lpthread -lm_r
#
#                point at the 64 bit version of the wsmp library (-l)
#                and include the directory (-L) where it is
#                located. two examples listed here. use the generic
#                RS-6000 version (wsmp) or SP-3 version (wsmpp3)
#
LIBS = -L/usr/local/apps/csm/wsmp/aix5.1 -lwsmp64
#
#                all object files for warp3d
#
OBJ =	 \
[OBJECTS]
	$(ODIR)/[object]$O	 \
[END_OBJECTS]


$(EXE_NAME) : $(OBJ)
	$(MPIF90) $(LOPT) $(OBJ) $(LIBS) -bloadmap:PARM -o $@
	$(MVC) $@ ./../run_ibm/$@
	chmod ugo+rx  ./../run_ibm/$@

[CUSTOM]
$(ODIR)/[f-file]$O :  [f-file].f [dependencies]
	$(RMC) $@
	$(RMC) filt_[f-file].f
	$(FILTER) < [f-file].f > filt_[f-file].f
	$(F90) [optimizations]   -c filt_[f-file].f       
	$(MVC) filt_[f-file]$O       [f-file]$O            
	$(MVC) [f-file]$O $@

[END_CUSTOM]
