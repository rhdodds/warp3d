#
# ****************************************************************************
#
#   File: Makefile.linux.inc
#
#   Developer: rdodds
#
#   Last modifications: 12/24/2019 rhd
#
#   This is invoked by a higher-level makefile:
#         Makefile.linux_Intel.omp        for threads-only executable
#         Makefile.linux_Intel.mpi_omp    for MPI + threads executable
#
#         Set before this makefile: COMPILER, HYINC, HYLIB, MPILIB, ODIR
#
#   Description:
#
#     Build OMP version of WARP3D (64-bit) running on
#     Intel/AMD 64-bit processors.
#   
#     Use Intel Fortran Compiler and possibly Intel MPI libraries
#
#     Includes the (free) Intel MKL library to 
#     access the PARDISO solver and LAPACK/BLAS routines. Those
#     routines are included in the WARP3D distribution directory
#     linux_packages/lib
#
#
#     Notes
#
#       1. Key compiler options:
#
#           Intel: by default we use the compiler feature to
#                  build instructions in the executable for all
#                  processor types from ivybridge to current.
#                  During execution the code dispatches the best
#                  instructions for your specific computer.
#
#                  A smaller executable file can be generated using
#                  -xHOST option as offered below
#
#       2. The Pardiso solver package in Intel MKL automatically
#          uses the most efficient instruction set available on
#          your specific computer.
#
#       3. The distributed executables for WARP3D ** require ** 
#          that your processor support the AVX instruction 
#          set from Intel ivybridge and newer processors
#          (as a minimum).
#
#          These processors became available in early 2011 and include:
#             2nd Generation Intel® Core™ i7 Processors
#             2nd Generation Intel® Core™ i5 Processors
#             2nd Generation Intel® Core™ i3 Processors
#             Intel® Xeon® Processor E5 Family
#             Intel® Xeon® Processor E3 Family
#             + all newer processors
#
#           During the next year, we plan to require processors
#           that support the AVX2 instruction set.
#
#       4. We build an executable using shared libraries for 
#          Intel MKL (and MPI if needed)
#
OSNAME   := Linux
BUILDNUM := 4100
#
#                Extension for relocatable object code files and
#                remove/move shell commands.
#
#                ODIR is set before make invoked. for simplicity
#                we use separate object code directories for
#                threaded only and threaded + MPI
#
O        = .o
OD       = $(ODIR)
RMC      = /bin/rm -f
MVC      = /bin/mv -f
#
#              common options other than used here:
#                 FTRAP = -ftrapuv -check uninit
#                 INLINE = -inline-level=0
#                 PROFILE = -pg for use with gprof, -g for Intel Adviser
#
FTRAP      =
INLINE     = -no-inline-factor # -fno-inline
PROFILE    = -g
PROFILE =
OPT_REPORT = -qopt-report=5 -qopt-report-annotate=html -qopt-report-phase=vec
OPT_REPORT = 
OPENMP     = -qopenmp -reentrancy
CHECKFLOAT = -fpe0 -init=snan,arrays -fp-speculation=safe
CHECKFLOAT =
#
ALL_PROCESSORS = -axAVX,CORE-AVX2,COMMON-AVX512 -xIVYBRIDGE \
                 -align array64byte  -qopt-zmm-usage=high
HARDWARE_SPECIFIC = -xHOST  -align array64byte  -qopt-zmm-usage=high
MACHINE_TYPE = $(ALL_PROCESSORS) # $(HARDWARE_SPECIFIC) 
#
FCOPTS = -fpp -O3 -ip -fpconstant -fixed -72 -traceback -ftz $(CHECKFLOAT) \
         $(MACHINE_TYPE) -module $(ODIR) $(HYINC) $(FTRAP) \
         $(INLINE) $(PROFILE) $(MPIOPT) $(OPENMP)
#
FOR      = $(COMPILER) $(OPT_REPORT) $(FCOPTS)
LINK     = $(COMPILER) $(FCOPTS)
#
#                all system level module files that appear as separate *.f in WARP3D
#                => does not include the user routine module .f for Abaqus
#                    compatible routines
#
MODS     := $(OD)/mod_adaptive_steps.o \
			$(OD)/mod_contact.o       \
			$(OD)/mod_crack_growth.o \
			$(OD)/mod_crystals.o \
			$(OD)/mod_damage.o \
			$(OD)/mod_elblk.o \
			$(OD)/mod_eleblocks.o \
			$(OD)/mod_elem_load.o \
			$(OD)/mod_file_info.o \
			$(OD)/mod_hypre.o \
			$(OD)/mod_jintegral.o \
			$(OD)/mod_local_stiffness.o \
			$(OD)/mod_main.o \
			$(OD)/mod_mpc.o \
			$(OD)/mod_mpi_lnpcg.o \
			$(OD)/mod_pconvert.o \
			$(OD)/mod_performance.o \
			$(OD)/mod_scan.o       \
			$(OD)/mod_segmental_curves.o \
			$(OD)/mod_stiffness.o \
			$(OD)/mod_trlist.o \
			$(OD)/distributed_assembly.o  # contain a module definition
#
INCLUDES := param_def include_mm10 include_sig_up include_tan_ek \
            include_usr_parm
#
#                MKL libraries to include at link time. Needed MPI
#                library for CPardiso set by invoking makefile (MPILIB).
#                The MKL files are part of WARP3D distribution
#
MKL_LIBS = ../linux_packages/lib
LIBMKL =  -L$(MKL_LIBS) -lmkl_intel_lp64 -lmkl_core \
          -lmkl_intel_thread -lpthread -lm
LIBS = $(LIBMKL) $(LIBMPI)  $(HYLIB)
#
#                all *.o files for linking.  module files are listed
#                first for readability. try to keep ordering alphabetic
#                after modules
#
OBJ =	 \
	$(OD)/mod_adaptive_steps.o	 \
	$(OD)/mod_contact.o	 \
	$(OD)/mod_crack_growth.o	 \
	$(OD)/mod_crystals.o   \
	$(OD)/mod_damage.o	 \
	$(OD)/mod_elblk.o	 \
	$(OD)/mod_elem_load.o	 \
	$(OD)/mod_eleblocks.o	 \
	$(OD)/mod_file_info.o	 \
	$(OD)/mod_hypre.o     \
	$(OD)/mod_jintegral.o	 \
	$(OD)/mod_local_stiffness.o \
	$(OD)/mod_main.o	 \
	$(OD)/mod_mpc.o	 \
	$(OD)/mod_mpi_lnpcg.o	 \
	$(OD)/mod_performance.o \
	$(OD)/mod_pconvert.o	 \
	$(OD)/mod_scan.o	 \
	$(OD)/mod_segmental_curves.o	 \
	$(OD)/mod_stiffness.o	 \
	$(OD)/mod_trlist.o	 \
	$(OD)/user_routines_module.o \
	$(OD)/addifv.o	 \
	$(OD)/addmas.o	 \
	$(OD)/adapt_check.o	 \
	$(OD)/allocate_damage.o	 \
	$(OD)/allocate_load.o	 \
	$(OD)/assemble_code.o	 \
	$(OD)/blcmp1.o	 \
	$(OD)/bmod.o	 \
	$(OD)/chkcrack.o	 \
	$(OD)/chkcrack_ctoa.o	 \
	$(OD)/chkcrack_gurson.o	 \
	$(OD)/chkcrack_smcs.o	 \
	$(OD)/chkcrack_cohesive.o	 \
	$(OD)/chk_stp_reduce.o	 \
	$(OD)/chk_stp_reduce_ctoa.o	 \
	$(OD)/chk_stp_reduce_gurson.o	 \
	$(OD)/chk_stp_reduce_smcs.o	 \
	$(OD)/chk_stp_reduce_cohes.o	 \
	$(OD)/cmpmas.o	 \
	$(OD)/cohes_elem.o	 \
	$(OD)/con_dump.o	 \
	$(OD)/contact.o	 \
	$(OD)/comput.o	 \
	$(OD)/cpu_time.o	 \
	$(OD)/cvtest.o	 \
	$(OD)/derivs.o	 \
	$(OD)/didriv.o	 \
	$(OD)/distup.o	 \
	$(OD)/distributed_assembly.o  \
	$(OD)/difrar.o	 \
	$(OD)/di1dsf.o	 \
	$(OD)/dimrot.o	 \
	$(OD)/dicmj.o	 \
	$(OD)/dielem_a.o	 \
	$(OD)/dielem_b.o	 \
	$(OD)/dielem_c.o	 \
	$(OD)/diexp4.o	 \
	$(OD)/dptstf.o	 \
	$(OD)/drive_eps_sig_internal_forces.o	 \
	$(OD)/drive_assemble_solve.o	 \
	$(OD)/drive_pardiso.o	 \
	$(OD)/drive_pardiso_asym.o	 \
	$(OD)/dupmas.o	 \
	$(OD)/dupstr.o	 \
	$(OD)/elmas1.o	 \
	$(OD)/elem_load_a.o	 \
	$(OD)/elem_load_b.o	 \
	$(OD)/elem_load_c.o	 \
	$(OD)/elprp.o	 \
	$(OD)/energy.o	 \
	$(OD)/eqivld.o	 \
	$(OD)/errchk.o	 \
	$(OD)/errmsg.o	 \
	$(OD)/find_neighbors.o	 \
	$(OD)/gatemps.o	 \
	$(OD)/get_edest.o	 \
	$(OD)/growth_loads.o	 \
	$(OD)/gp_temp_eps.o	 \
	$(OD)/gpifv1.o	 \
	$(OD)/gpmas1.o	 \
	$(OD)/gptns1.o	 \
	$(OD)/gtlsn1.o	 \
	$(OD)/gtmat1.o	 \
	$(OD)/incomp.o	 \
	$(OD)/inclmass.o	 \
	$(OD)/incontact.o	 \
	$(OD)/incon.o	 \
	$(OD)/incoor.o	 \
	$(OD)/incrystal.o \
	$(OD)/insurf.o	 \
	$(OD)/intied.o	 \
	$(OD)/indypm.o	 \
	$(OD)/inelbk.o	 \
	$(OD)/inelem.o	 \
	$(OD)/infile.o	 \
	$(OD)/incurv.o	 \
	$(OD)/indom.o	 \
	$(OD)/inalpha.o	 \
	$(OD)/incrack.o	 \
	$(OD)/infgm.o	 \
	$(OD)/inicon.o	 \
	$(OD)/ininc.o	 \
	$(OD)/initdm.o	 \
	$(OD)/initst.o	 \
	$(OD)/inlod.o	 \
	$(OD)/inmat.o	 \
	$(OD)/innum.o	 \
	$(OD)/intab.o	 \
	$(OD)/inpatch.o	 \
	$(OD)/iodevn.o	 \
	$(OD)/iterative_sparse_hypre.o \
	$(OD)/kg1.o	 \
	$(OD)/kill_element.o	 \
	$(OD)/license.o	 \
	$(OD)/linked_list.o	 \
	$(OD)/main_program.o	 \
	$(OD)/mnralg.o	 \
	$(OD)/mem_allocate.o	 \
	$(OD)/mm01.o	 \
	$(OD)/mm02.o	 \
	$(OD)/mm03.o	 \
	$(OD)/mm04.o	 \
	$(OD)/mm05.o	 \
	$(OD)/mm06.o	 \
	$(OD)/mm07.o	 \
	$(OD)/mm08.o	 \
	$(OD)/mm09.o	 \
	$(OD)/mm10_a.o	 \
	$(OD)/mm10_b.o	 \
	$(OD)/mm10_c.o	 \
	$(OD)/mm10_d.o	 \
	$(OD)/mm10_e.o	 \
	$(OD)/mm10_f.o   \
	$(OD)/mm11.o       \
	$(OD)/mmnoopt.o	\
	$(OD)/mpi_code.o	 \
	$(OD)/mpi_handle_slaves.o \
	$(OD)/mpc_modify.o	 \
	$(OD)/name_strip.o	 \
	$(OD)/ndpts1.o	 \
	$(OD)/ouddpa.o	 \
	$(OD)/oudriv.o	 \
	$(OD)/oudups.o	 \
	$(OD)/oudva.o	 \
	$(OD)/ouext1.o	 \
	$(OD)/ouext2.o	 \
	$(OD)/ouflnm.o	 \
	$(OD)/ougts1.o	 \
	$(OD)/ouhdva.o	 \
	$(OD)/ouhel.o	 \
	$(OD)/ouhnod.o	 \
	$(OD)/ouhprn.o	 \
	$(OD)/ouhrks.o	 \
	$(OD)/ouhstr.o	 \
	$(OD)/oulbdd.o	 \
	$(OD)/oulbir.o	 \
	$(OD)/oulbst.o	 \
	$(OD)/oulg1.o	 \
	$(OD)/oumkcv.o	 \
	$(OD)/oundhd.o	 \
	$(OD)/ounds1.o	 \
	$(OD)/ouocdd.o	 \
	$(OD)/ouocst.o	 \
	$(OD)/ouneut.o	 \
	$(OD)/ouocst_elem.o	 \
	$(OD)/oupdva.o	 \
	$(OD)/oupele.o	 \
	$(OD)/ouprks.o	 \
	$(OD)/oupstr_elem.o	 \
	$(OD)/oupstr_node.o	 \
	$(OD)/oures.o	 \
	$(OD)/ourlhd.o	 \
	$(OD)/ousthd.o	 \
	$(OD)/oust_elem.o	 \
	$(OD)/oustates.o	 \
	$(OD)/outfil.o	 \
	$(OD)/oustr.o	 \
	$(OD)/outime.o	 \
	$(OD)/packets.o	 \
	$(OD)/pconvert.o	 \
	$(OD)/prcsel.o	 \
	$(OD)/qmply1.o	 \
	$(OD)/getgpts.o	 \
	$(OD)/read_damage.o	 \
	$(OD)/read_write.o	 \
	$(OD)/reopen.o	 \
	$(OD)/rknifv.o	 \
	$(OD)/rknmas.o	 \
	$(OD)/rknstr.o	 \
	$(OD)/rktstf.o	 \
	$(OD)/rplstr.o	 \
	$(OD)/rscmp1.o	 \
	$(OD)/rstgp1.o	 \
	$(OD)/scan.o	 \
	$(OD)/scstr.o	 \
	$(OD)/setelb.o	 \
	$(OD)/setup.o	 \
	$(OD)/shapef.o	 \
	$(OD)/setup_mat_props.o	 \
	$(OD)/splunj.o	 \
	$(OD)/statistics.o	 \
	$(OD)/steptime.o   \
	$(OD)/store.o	 \
	$(OD)/store_ifv.o	 \
	$(OD)/stpdrv.o	 \
	$(OD)/star_com.o	 \
	$(OD)/tanstf.o	 \
	$(OD)/temsto.o	 \
	$(OD)/thyme.o	 \
	$(OD)/trnvec.o	 \
	$(OD)/trnmtx.o	 \
	$(OD)/tied_mesh_a.o	 \
	$(OD)/tied_mesh_b.o	 \
	$(OD)/tied_mesh_c.o	 \
	$(OD)/tilde.o	 \
	$(OD)/update.o	 \
	$(OD)/update_node_elecnt.o	 \
	$(OD)/uppbar.o	 \
	$(OD)/upres.o	 \
	$(OD)/user_list.o	 \
	$(OD)/user_routines_umat.o	 \
	$(OD)/user_routines_other_Abaqus.o \
	$(OD)/user_routines_warp3d.o	 \
	$(OD)/user_routines_support.o	 \
	$(OD)/vol_terms.o	 \
	$(OD)/vol_avg.o	 \
	$(OD)/zero_vector.o	 \
	$(OD)/zero_vol.o
#
$(EXE_NAME) : $(OBJ)
	$(LINK) $(OBJ) $(LIBS) -o $@
	$(MVC) $@ ./../run_linux/$@
	chmod ugo+rx  ./../run_linux/$@
#
#       ---  default rule to make .o from .f ------
#            ===============================
#
#            => for simplicity we assume all the *.f depend on all
#               modules, includes, param_def
#
$(OD)/%.o : %.f $(MODS) $(INCLUDES)
	$(RMC)  $@
	$(FOR) $< -c -o $@
#        $(FOR) -qopt-report-file=$<.html $< -c -o $@
#
#
#       ---  main_program.f gets special ifort options ------
#            =========================================
#
$(OD)/main_program$O : main_program.f $(MODS) $(INCLUDES)
	$(FOR) -D$(OSNAME) -DBuild=$(BUILDNUM) -fpp  $< -c -o $@
#
#       ---  any .f files that require non-std compile options ------
#            =================================================
#
$(OD)/mm10_b$O : mm10_b.f $(OD)/mod_crystals$O param_def include_mm10
	$(RMC)  $@
	$(FOR) -inline-level=1  $< -c -o $@
#
#       ---  Abaqus compatible user routines. Listed explicitly for clarity ------
#            ==============================================================
#
$(ODIR)/user_routines_module$O : user_routines_module.f
	$(RMC)  $@
	$(FOR)  $< -c -o $@

$(OD)/user_routines_umat$O : user_routines_umat.f $(OD)/user_routines_module$O
	$(RMC)  $@
	$(FOR)  -inline-level=1  $< -c -o $@

$(OD)/user_routines_other_Abaqus$O : user_routines_other_Abaqus.f \
                                     $(OD)/user_routines_module$O
	$(RMC)  $@
	$(FOR) -heap-arrays $< -c -o $@

$(OD)/user_routines_warp3d$O : user_routines_warp3d.f  include_usr_parm \
                               $(OD)/user_routines_module$O
	$(RMC)  $@
	$(FOR) -heap-arrays  $< -c -o $@

$(OD)/user_routines_support$O : user_routines_support.f  include_usr_parm \
                                $(OD)/user_routines_module$O
	$(RMC)  $@
	$(FOR)  $< -c -o $@
#
#
#
#    ------   system level module files ----------------
#             =========================
#
$(ODIR)/mod_adaptive_steps$O : mod_adaptive_steps.f
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_contact$O : mod_contact.f
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_crack_growth$O : mod_crack_growth.f
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_crystals$O : mod_crystals.f param_def  $(OD)/mod_main$O
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_damage$O : mod_damage.f
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_elblk$O : mod_elblk.f param_def
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_eleblocks$O : mod_eleblocks.f param_def
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_elem_load$O : mod_elem_load.f param_def
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_file_info$O : mod_file_info.f
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_hypre$O : mod_hypre.f  param_def
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_jintegral$O : mod_jintegral.f param_def
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_local_stiffness$O : mod_local_stiffness.f  param_def
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_main$O : mod_main.f
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_mpc$O : mod_mpc.f param_def
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_mpi_lnpcg$O : mod_mpi_lnpcg.f param_def
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_pconvert$O : mod_pconvert.f
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_performance$O : mod_performance.f
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_scan$O : mod_scan.f
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_segmental_curves$O : mod_segmental_curves.f param_def
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_stiffness$O : mod_stiffness.f param_def
	$(RMC)  $@
	$(FOR) $< -c -o $@

$(ODIR)/mod_trlist$O : mod_trlist.f param_def $(OD)/mod_main$O $(OD)/mod_scan$O
	$(RMC)  $@
	$(FOR) $< -c -o $@

