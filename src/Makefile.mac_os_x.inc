#
# ********************************************************************
#
#   File: mac_os_x.inc
#
#   Developer: rdodds
#
#   Last modifications: December 19, 2016
#
#   Description:
#
#     Build OMP version of WARP3D (64-bit) running on
#     Mac OS X, Intel 64-bit processors.
#
#     This file must be invoked by Makefile.mac_os.omp
#     Several macro values are set in "Makewarp.bash"
#     Users: only build the OS X executable thru Makewarp.bash
#
#     See README for Mac OS X for description of req'd compilers that
#     must be installed.
#
#   Revision History:
#
#   mcm - Changed the location of the MKL and iomp5 libraries to be
#         configuration by passing in variables MKLLIB and FORLIB
#         Updated target list for hypre (dummies and module).
#         Added a third pass (really, Intel, really) over the libraries
#         to pick up the last few symbols.
#
#   rhd - added mcm code for crystal plasticity and support for hypre
#         running with distributed stiffness (MPI only).
#
#   rhd - major cleanup and dependency fixes. make consistent with updated
#         Linux version.
#
#   rhd - remove the nasa and vss solvers, the old l-d-u profile solver
#
#   rhd - remove EBE solver code
#
#   rhd - change driver.f -> main_program.f
#   rhd - add dependencies for umats.f and stpdrv.f
#   rhd - change umats.f to user_routines.f
#   rhd - add compiles for material model 11
#   rhd - add compile for oustates.f
#   rhd - add compile for mm08.f
#   rhd - change to "simple" linking model for mkl (still static)
#   rhd - delete mod_pcg, addstf, oulres
#   rhd - user_routines divided into 4 files
#   rhd - remove newmrk.f upcnld.f
#   rhd - remove outrac.f
#   rhd - remove ouconv.f
#   rhd - remove oultrc.f
#   rhd - add -heap-arrays on some *.f
#   rhd - delete stifup.f
#   rhd - add user_routines_module.f, update dependencies for
#         user_routines_...f, reorder some rules 
#         for convenience
#   rhd   add -xHost to compiler options
#   rhd - replace mm10 with mm10a, b, c
#         add rstgp1 dependency on mod_crystals
#         use -inline-level=1 on some very long .f files
#         to reduce unnecessarily long compile times w/o 
#         improving runtime performance
#   rhd   1/4/2016  remove gplns1.f, lnstff.f, rklstf.f, lcnst1.f,
#         cnst1.f, cnst2.f, cnst3.f, cnst4.f
#   rhd   1/7/2016  mm10_a.f at -O3 for ifort V16.x causes seg fault.
#                   isolated the culprit routine into a new mm10_d.f
#                   and compiled at -O2 which runs fine.
#   rhd   4/6/2016  add mm10_e.f
#   rhd   5/6/2016  add align64 to some .f
#   rhd   6/7/2016  add mm10_f.f and make oustates use -O3
#   rhd   6/18/2016 renove f_clmloc, mod_crystals.f dependency
#   rhd   8/17/2016 added new macro F90TRAP
#   rhd   8/20/2016 combined mm03p.f & mm03.f into mm03.f for more inlining
#   rhd  12/19/2016 more macros to generalize building
#                   removal of filter program and much simplifications
#
#
#
#    NOTE:   we specify the -xHost Fortran option. the compiler will
#            generate code for the highest performance instruction set of the 
#            CPU on your computer.
#
#            Remove this option is you wish the compiler to generate
#            instructions for the broadest range (old !) Intel
#            processors
#
#
# ***********************************************************************
#
#
#        COMPILER parameter set by the Makewarp.bash
#        script before the Makefile.mac_os_x.omp is invoked using make.
#
#        ODIR, OPENMPOPT parameters are set by
#        the "pre" makefile: Makefile.mac_os_x.omp
#
#        Linking stategy. We build a statically linked executable for
#        WARP3D with the simplest linking model for MKL
#
#        Intel F90 (fort) compiler options. see also Makefile.mac_os_x.omp
#
OSNAME   := OSX
BUILDNUM := 3250
PRE      := fpp

#F90TRAP = -ftrapuv -check uninit
F90TRAP = 
#
INLINEOPT =  -no-inline-factor
#INLINEOPT =  -inline-level=0   # no inlining
#
PROFILEOPT =
#PROFILEOPT = -pg  # for use with gprof
#PROFILEOPT = -g   # for use with Intel Adviser
#
F90OPTS1 =  -fpconstant -fixed -72 -traceback $(PROFILEOPT)
F90OPTS2 =  -module $(ODIR) -ftz -xHost  $(F90TRAP) $(INLINEOPT)
#
F90      = $(COMPILER) $(F90OPTS1) $(F90OPTS2) $(OPENMPOPT)
#
LINK90   = $(COMPILER) $(F90OPTS1) $(F90OPTS2) $(OPENMPOPT) 
LIBS     = -mkl
#
#                Define extension for relocatable
#                object code files, macros for the remove command and the
#                move command.
#
O        = .o
OD       = $(ODIR)
RMC      = rm -f
MVC      = mv -f
#
#                Program will insert all object files for warp3d. All
#                module files are listed first.
#
OBJ =	 \
	$(OD)/mod_adaptive_steps$O	 \
	$(OD)/mod_contact$O	 \
	$(OD)/mod_crack_growth$O	 \
	$(OD)/mod_crystals$O   \
	$(OD)/mod_damage$O	 \
	$(OD)/mod_elblk$O	 \
	$(OD)/mod_elem_load$O	 \
	$(OD)/mod_eleblocks$O	 \
	$(OD)/mod_file_info$O	 \
	$(OD)/mod_hypre$O     \
	$(OD)/mod_jintegral$O	 \
	$(OD)/mod_local_stiffness$O \
	$(OD)/mod_main$O	 \
	$(OD)/mod_mpc$O	 \
	$(OD)/mod_mpi_lnpcg$O	 \
	$(OD)/mod_performance$O \
	$(OD)/mod_pconvert$O	 \
	$(OD)/mod_segmental_curves$O	 \
	$(OD)/mod_stiffness$O	 \
	$(OD)/user_routines_module$O \
	$(OD)/addifv$O	 \
	$(OD)/addmas$O	 \
	$(OD)/adapt_check$O	 \
	$(OD)/aplacc$O	 \
	$(OD)/allocate_damage$O	 \
	$(OD)/allocate_load$O	 \
	$(OD)/assemble_code$O	 \
	$(OD)/bdbt$O	 \
	$(OD)/blcmp1$O	 \
	$(OD)/bmod$O	 \
	$(OD)/chkcrack$O	 \
	$(OD)/chkcrack_ctoa$O	 \
	$(OD)/chkcrack_gurson$O	 \
	$(OD)/chkcrack_smcs$O	 \
	$(OD)/chkcrack_cohesive$O	 \
	$(OD)/chk_stp_reduce$O	 \
	$(OD)/chk_stp_reduce_ctoa$O	 \
	$(OD)/chk_stp_reduce_gurson$O	 \
	$(OD)/chk_stp_reduce_smcs$O	 \
	$(OD)/chk_stp_reduce_cohes$O	 \
	$(OD)/cmpmas$O	 \
	$(OD)/cohes_elem$O	 \
	$(OD)/con_dump$O	 \
	$(OD)/contact$O	 \
	$(OD)/comput$O	 \
	$(OD)/cpu_time$O	 \
	$(OD)/cvtest$O	 \
	$(OD)/derivs$O	 \
	$(OD)/didriv$O	 \
	$(OD)/distup$O	 \
	$(OD)/distributed_assembly$O  \
	$(OD)/difrar$O	 \
	$(OD)/difrtn$O	 \
	$(OD)/di1dsf$O	 \
	$(OD)/dimrot$O	 \
	$(OD)/dicmj$O	 \
	$(OD)/dielem_a$O	 \
	$(OD)/dielem_b$O	 \
	$(OD)/dielem_c$O	 \
	$(OD)/direct_driver_sparse$O	 \
	$(OD)/direct_sparse_mkl$O	 \
	$(OD)/diexp4$O	 \
	$(OD)/dptstf$O	 \
	$(OD)/drive_eps_sig_internal_forces$O	 \
	$(OD)/dupmas$O	 \
	$(OD)/dupstr$O	 \
	$(OD)/elmas1$O	 \
	$(OD)/elem_load_a$O	 \
	$(OD)/elem_load_b$O	 \
	$(OD)/elem_load_c$O	 \
	$(OD)/elprp$O	 \
	$(OD)/energy$O	 \
	$(OD)/eqivld$O	 \
	$(OD)/equiv_strain$O	 \
	$(OD)/errchk$O	 \
	$(OD)/errmsg$O	 \
	$(OD)/find_neighbors$O	 \
	$(OD)/gastr$O	 \
	$(OD)/gatemps$O	 \
	$(OD)/get_edest$O	 \
	$(OD)/growth_loads$O	 \
	$(OD)/gp_temp_eps$O	 \
	$(OD)/gpifv1$O	 \
	$(OD)/gpmas1$O	 \
	$(OD)/gptns1$O	 \
	$(OD)/gtlsn1$O	 \
	$(OD)/gtmat1$O	 \
	$(OD)/hypre$O		 \
	$(OD)/incomp$O	 \
	$(OD)/inclmass$O	 \
	$(OD)/incontact$O	 \
	$(OD)/incon$O	 \
	$(OD)/incoor$O	 \
	$(OD)/incrystal$O \
	$(OD)/insurf$O	 \
	$(OD)/intied$O	 \
	$(OD)/indypm$O	 \
	$(OD)/inelbk$O	 \
	$(OD)/inelem$O	 \
	$(OD)/infile$O	 \
	$(OD)/incurv$O	 \
	$(OD)/indom$O	 \
	$(OD)/inalpha$O	 \
	$(OD)/incrack$O	 \
	$(OD)/infgm$O	 \
	$(OD)/inicon$O	 \
	$(OD)/ininc$O	 \
	$(OD)/initdm$O	 \
	$(OD)/initst$O	 \
	$(OD)/inlod$O	 \
	$(OD)/inmat$O	 \
	$(OD)/innum$O	 \
	$(OD)/intab$O	 \
	$(OD)/iodevn$O	 \
	$(OD)/kg1$O	 \
	$(OD)/kill_element$O	 \
	$(OD)/license$O	 \
	$(OD)/linked_list$O	 \
	$(OD)/main_program$O	 \
	$(OD)/mnralg$O	 \
	$(OD)/mem_allocate$O	 \
	$(OD)/mm01$O	 \
	$(OD)/mm02$O	 \
	$(OD)/mm03$O	 \
	$(OD)/mm04$O	 \
	$(OD)/mm05$O	 \
	$(OD)/mm06$O	 \
	$(OD)/mm07$O	 \
	$(OD)/mm08$O	 \
	$(OD)/mm09$O	 \
	$(OD)/mm10_a$O	 \
	$(OD)/mm10_b$O	 \
	$(OD)/mm10_c$O	 \
	$(OD)/mm10_d$O	 \
	$(OD)/mm10_e$O	 \
	$(OD)/mm10_f$O	 \
	$(OD)/mm11$O	 \
	$(OD)/mmnoopt$O	 \
	$(OD)/mpi_code$O	 \
	$(OD)/mpi_handle_slaves$O \
	$(OD)/mpc_modify$O	 \
	$(OD)/name_strip$O	 \
	$(OD)/ndpts1$O	 \
	$(OD)/ou3dpr$O	 \
	$(OD)/oualph$O	 \
	$(OD)/ouddpa$O	 \
	$(OD)/oudriv$O	 \
	$(OD)/oudups$O	 \
	$(OD)/oudva$O	 \
	$(OD)/oueff1$O	 \
	$(OD)/ouext1$O	 \
	$(OD)/ouext2$O	 \
	$(OD)/ouflnm$O	 \
	$(OD)/ougts1$O	 \
	$(OD)/ouhdva$O	 \
	$(OD)/ouhel$O	 \
	$(OD)/ouhnod$O	 \
	$(OD)/ouhprn$O	 \
	$(OD)/ouhrks$O	 \
	$(OD)/ouhstr$O	 \
	$(OD)/ouinv1$O	 \
	$(OD)/oulbdd$O	 \
	$(OD)/oulbir$O	 \
	$(OD)/oulbst$O	 \
	$(OD)/oulbta$O	 \
	$(OD)/oulg1$O	 \
	$(OD)/oulrlh$O	 \
	$(OD)/oumkcv$O	 \
	$(OD)/oundhd$O	 \
	$(OD)/ounds1$O	 \
	$(OD)/ouocdd$O	 \
	$(OD)/ouocst$O	 \
	$(OD)/ouneut$O	 \
	$(OD)/ouocst_elem$O	 \
	$(OD)/oupads$O	 \
	$(OD)/ouparm$O	 \
	$(OD)/oupdva$O	 \
	$(OD)/oupele$O	 \
	$(OD)/oupri1$O	 \
	$(OD)/ouprks$O	 \
	$(OD)/oupstr_elem$O	 \
	$(OD)/oupstr_node$O	 \
	$(OD)/oures$O	 \
	$(OD)/ourlhd$O	 \
	$(OD)/ousthd$O	 \
	$(OD)/oustpa$O	 \
	$(OD)/oust_elem$O	 \
	$(OD)/oustates$O	 \
	$(OD)/outfil$O	 \
	$(OD)/oustr$O	 \
	$(OD)/outime$O	 \
	$(OD)/ouyld1$O	 \
	$(OD)/packets$O	 \
	$(OD)/pconvert$O	 \
	$(OD)/princ_inv_strain$O	 \
	$(OD)/princ_inv_stress$O	 \
	$(OD)/princ_strain$O	 \
	$(OD)/princ_stress$O	 \
	$(OD)/prcsel$O	 \
	$(OD)/qmply1$O	 \
	$(OD)/getgpts$O	 \
	$(OD)/read_damage$O	 \
	$(OD)/read_write$O	 \
	$(OD)/reopen$O	 \
	$(OD)/rknifv$O	 \
	$(OD)/rknmas$O	 \
	$(OD)/rknstr$O	 \
	$(OD)/rktstf$O	 \
	$(OD)/rplstr$O	 \
	$(OD)/rscmp1$O	 \
	$(OD)/rstgp1$O	 \
	$(OD)/scan$O	 \
	$(OD)/scstr$O	 \
	$(OD)/setelb$O	 \
	$(OD)/setup$O	 \
	$(OD)/shapef$O	 \
	$(OD)/setup_mat_props$O	 \
	$(OD)/splunj$O	 \
	$(OD)/statistics$O	 \
	$(OD)/steptime$O   \
	$(OD)/store$O	 \
	$(OD)/store_ifv$O	 \
	$(OD)/stpdrv$O	 \
	$(OD)/star_com$O	 \
	$(OD)/tanstf$O	 \
	$(OD)/temsto$O	 \
	$(OD)/thyme$O	 \
	$(OD)/trnvec$O	 \
	$(OD)/trnmtx$O	 \
	$(OD)/tied_mesh$O	 \
	$(OD)/tilde$O	 \
	$(OD)/update$O	 \
	$(OD)/update_node_elecnt$O	 \
	$(OD)/uppbar$O	 \
	$(OD)/upres$O	 \
	$(OD)/user_list$O	 \
	$(OD)/user_routines_umat$O	 \
	$(OD)/user_routines_other_Abaqus$O \
	$(OD)/user_routines_warp3d$O	 \
	$(OD)/user_routines_support$O	 \
	$(OD)/vol_avg$O	 \
	$(OD)/vol_terms$O	 \
	$(OD)/yield_function$O	 \
	$(OD)/zero_vector$O	 \
	$(OD)/zero_vol$O	 \
#
#
$(EXE_NAME) : $(OBJ)
	$(LINK90) $(OBJ) $(LIBS) -o $@
	$(MVC) $@ ./../run_mac_os_x/$@
	chmod ugo+rx  ./../run_mac_os_x/$@
#
$(OD)/addifv$O : addifv.f param_def \
                   $(OD)/mod_damage$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c addifv.f
	$(MVC) addifv$O $@

$(OD)/adapt_check$O : adapt_check.f $(OD)/mod_adaptive_steps$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c adapt_check.f
	$(MVC) adapt_check$O $@

$(OD)/addmas$O : addmas.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c addmas.f
	$(MVC) addmas$O $@

$(OD)/allocate_damage$O : allocate_damage.f $(OD)/mod_crack_growth$O \
                   common.main param_def \
                   $(OD)/mod_damage$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c allocate_damage.f
	$(MVC) allocate_damage$O $@

$(OD)/allocate_load$O : allocate_load.f $(OD)/mod_elem_load$O \
                   common.main param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c allocate_load.f
	$(MVC) allocate_load$O $@

$(OD)/aplacc$O : aplacc.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c aplacc.f
	$(MVC) aplacc$O $@

$(OD)/assemble_code$O : assemble_code.f common.main \
                   param_def $(OD)/mod_eleblocks$O \
                   $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O3 -ip  -c assemble_code.f
	$(MVC) assemble_code$O $@
#
#       --- b --
#
$(OD)/bdbt$O : bdbt.f
	$(RMC)  $@
	$(F90) -O3 -ip  -c bdbt.f
	$(MVC) bdbt$O $@

$(OD)/blcmp1$O : blcmp1.f param_def
	$(RMC)  $@
	$(F90) -O3 -ip  -c blcmp1.f
	$(MVC) blcmp1$O $@

$(OD)/bmod$O : bmod.f
	$(RMC)  $@
	$(F90) -O3 -ip  -c bmod.f
	$(MVC) bmod$O $@

$(OD)/bscmp1$O : bscmp1.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c bscmp1.f
	$(MVC) bscmp1$O $@

$(OD)/btcmp1$O : btcmp1.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c btcmp1.f
	$(MVC) btcmp1$O $@
#
#       --- c --
#
$(OD)/chkcrack$O : chkcrack.f $(OD)/mod_crack_growth$O \
                   common.main $(OD)/mod_damage$O \
                   param_def $(OD)/mod_eleblocks$O \
                   $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c chkcrack.f
	$(MVC) chkcrack$O $@

$(OD)/chkcrack_ctoa$O : chkcrack_ctoa.f $(OD)/mod_crack_growth$O \
                   common.main $(OD)/mod_damage$O \
                   param_def  $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c chkcrack_ctoa.f
	$(MVC) chkcrack_ctoa$O $@

$(OD)/chkcrack_gurson$O : chkcrack_gurson.f $(OD)/mod_crack_growth$O \
                   common.main $(OD)/mod_damage$O \
                   param_def $(OD)/mod_eleblocks$O \
                   $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c chkcrack_gurson.f
	$(MVC) chkcrack_gurson$O $@

$(OD)/chkcrack_smcs$O : chkcrack_smcs.f $(OD)/mod_crack_growth$O \
                   common.main $(OD)/mod_damage$O \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c chkcrack_smcs.f
	$(MVC) chkcrack_smcs$O $@

$(OD)/chkcrack_cohesive$O : chkcrack_cohesive.f $(OD)/mod_crack_growth$O \
                   common.main $(OD)/mod_damage$O \
                   param_def $(OD)/mod_eleblocks$O \
                   $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c chkcrack_cohesive.f
	$(MVC) chkcrack_cohesive$O $@

$(OD)/chk_stp_reduce$O : chk_stp_reduce.f  \
                   common.main $(OD)/mod_damage$O \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c chk_stp_reduce.f
	$(MVC) chk_stp_reduce$O $@

$(OD)/chk_stp_reduce_ctoa$O : chk_stp_reduce_ctoa.f $(OD)/mod_crack_growth$O \
                   common.main $(OD)/mod_damage$O \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c chk_stp_reduce_ctoa.f
	$(MVC) chk_stp_reduce_ctoa$O $@

$(OD)/chk_stp_reduce_gurson$O : chk_stp_reduce_gurson.f $(OD)/mod_crack_growth$O \
                   common.main $(OD)/mod_damage$O \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c chk_stp_reduce_gurson.f
	$(MVC) chk_stp_reduce_gurson$O $@

$(OD)/chk_stp_reduce_smcs$O : chk_stp_reduce_smcs.f $(OD)/mod_crack_growth$O \
                   common.main $(OD)/mod_damage$O \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c chk_stp_reduce_smcs.f
	$(MVC) chk_stp_reduce_smcs$O $@

$(OD)/chk_stp_reduce_cohes$O : chk_stp_reduce_cohes.f $(OD)/mod_crack_growth$O \
                   common.main $(OD)/mod_damage$O \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c chk_stp_reduce_cohes.f
	$(MVC) chk_stp_reduce_cohes$O $@

$(OD)/cmpmas$O : cmpmas.f $(OD)/mod_eleblocks$O \
                   common.main param_def \
                   $(OD)/mod_main$O $(OD)/mod_mpi_lnpcg$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c cmpmas.f
	$(MVC) cmpmas$O $@

$(OD)/cohes_elem$O : cohes_elem.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c cohes_elem.f
	$(MVC) cohes_elem$O $@

$(OD)/comput$O : comput.f common.main \
                   param_def $(OD)/mod_jintegral$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c comput.f
	$(MVC) comput$O $@

$(OD)/con_dump$O : con_dump.f common.main \
                   param_def $(OD)/mod_main$O \
                   $(OD)/mod_mpc$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c con_dump.f
	$(MVC) con_dump$O $@

$(OD)/contact$O : contact.f common.main \
                   $(OD)/mod_contact$O param_def \
                   $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c contact.f
	$(MVC) contact$O $@

$(OD)/cpu_time$O : cpu_time.f  $(OD)/mod_performance$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c cpu_time.f
	$(MVC) cpu_time$O $@

$(OD)/cvtest$O : cvtest.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c cvtest.f
	$(MVC) cvtest$O $@
#
#       --- d --
#
$(OD)/derivs$O : derivs.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c derivs.f
	$(MVC) derivs$O $@

$(OD)/direct_driver_sparse$O : direct_driver_sparse.f param_def \
                   common.main $(OD)/mod_eleblocks$O \
                   $(OD)/mod_main$O $(OD)/mod_mpc$O \
                   $(OD)/mod_stiffness$O  $(OD)/mod_hypre$O  \
                   $(OD)/mod_local_stiffness$O $(OD)/distributed_assembly$O \
                   $(OD)/mod_performance$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c direct_driver_sparse.f
	$(MVC) direct_driver_sparse$O $@

$(OD)/direct_sparse_mkl$O : direct_sparse_mkl.f
	$(RMC)  $@
	$(F90) -O3 -ip  -c direct_sparse_mkl.f
	$(MVC) direct_sparse_mkl$O $@

$(OD)/distributed_assembly$O : distributed_assembly.f \
                                 $(OD)/mod_local_stiffness$O  \
                                 param_def
	$(RMC)  $@
	$(F90) -O3 -ip -c distributed_assembly.f
	$(MVC) distributed_assembly$O $@

$(OD)/didriv$O : didriv.f common.main  $(OD)/mod_segmental_curves$O \
                   param_def $(OD)/mod_jintegral$O \
                   $(OD)/mod_main$O $(OD)/mod_eleblocks$O
	$(RMC)  $@
	$(F90) -O1 -c didriv.f
	$(MVC) didriv$O $@

$(OD)/diexp4$O : diexp4.f param_def \
                   $(OD)/mod_jintegral$O $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c diexp4.f
	$(MVC) diexp4$O $@

$(OD)/distup$O : distup.f $(OD)/mod_jintegral$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c distup.f
	$(MVC) distup$O $@

$(OD)/difrar$O : difrar.f $(OD)/mod_jintegral$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c difrar.f
	$(MVC) difrar$O $@

$(OD)/difrtn$O : difrtn.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c difrtn.f
	$(MVC) difrtn$O $@

$(OD)/di1dsf$O : di1dsf.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c di1dsf.f
	$(MVC) di1dsf$O $@

$(OD)/dimrot$O : dimrot.f common.main \
                   $(OD)/mod_main$O $(OD)/mod_jintegral$O \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c dimrot.f
	$(MVC) dimrot$O $@
#
#         note lower opt level
#
$(OD)/dicmj$O : dicmj.f param_def \
                   $(OD)/mod_jintegral$O common.main \
                   $(OD)/mod_main$O $(OD)/mod_eleblocks$O
	$(RMC)  $@
	$(F90) -O0 -ip -c dicmj.f
	$(MVC) dicmj$O $@

$(OD)/dielem_a$O : dielem_a.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c dielem_a.f
	$(MVC) dielem_a$O $@

$(OD)/dielem_b$O : dielem_b.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c dielem_b.f
	$(MVC) dielem_b$O $@

$(OD)/dielem_c$O : dielem_c.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c dielem_c.f
	$(MVC) dielem_c$O $@

$(OD)/dptstf$O : dptstf.f common.main \
                   param_def  $(OD)/mod_damage$O \
                   $(OD)/mod_crack_growth$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c dptstf.f
	$(MVC) dptstf$O $@

$(OD)/drive_eps_sig_internal_forces$O : drive_eps_sig_internal_forces.f \
                   $(OD)/mod_crack_growth$O $(OD)/mod_damage$O common.main \
                   param_def include_sig_up \
                   $(OD)/mod_eleblocks$O $(OD)/mod_main$O \
                   $(OD)/mod_segmental_curves$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c drive_eps_sig_internal_forces.f
	$(MVC) drive_eps_sig_internal_forces$O $@

$(OD)/dupmas$O : dupmas.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c dupmas.f
	$(MVC) dupmas$O $@

$(OD)/dupstr$O : dupstr.f common.main \
                   param_def $(OD)/mod_damage$O \
                   $(OD)/mod_crack_growth$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c dupstr.f
	$(MVC) dupstr$O $@
#
#       --- e --
#
$(OD)/elem_load_a$O : elem_load_a.f common.main \
                   $(OD)/mod_elem_load$O param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c elem_load_a.f
	$(MVC) elem_load_a$O $@

$(OD)/elem_load_b$O : elem_load_b.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c elem_load_b.f
	$(MVC) elem_load_b$O $@

$(OD)/elem_load_c$O : elem_load_c.f  common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c elem_load_c.f
	$(MVC) elem_load_c$O $@

$(OD)/elmas1$O : elmas1.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c elmas1.f
	$(MVC) elmas1$O $@

$(OD)/elprp$O : elprp.f common.main \
                   param_def $(OD)/mod_main$O \
                   $(OD)/mod_segmental_curves$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c elprp.f
	$(MVC) elprp$O $@

$(OD)/energy$O : energy.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c energy.f
	$(MVC) energy$O $@

$(OD)/eqivld$O : eqivld.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c eqivld.f
	$(MVC) eqivld$O $@

$(OD)/equiv_strain$O : equiv_strain.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c equiv_strain.f
	$(MVC) equiv_strain$O $@

$(OD)/errchk$O : errchk.f param_def common.main \
                   $(OD)/mod_crystals$O $(OD)/mod_damage$O \
                   $(OD)/mod_contact$O $(OD)/mod_main$O \
                   $(OD)/mod_segmental_curves$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c errchk.f
	$(MVC) errchk$O $@

$(OD)/errmsg$O : errmsg.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c errmsg.f
	$(MVC) errmsg$O $@

$(OD)/evcmp1$O : evcmp1.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c evcmp1.f
	$(MVC) evcmp1$O $@
#
#       --- f --
#
$(OD)/find_neighbors$O : find_neighbors.f param_def \
                   common.main $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c find_neighbors.f
	$(MVC) find_neighbors$O $@
#
#       --- g --
#
$(OD)/gastr$O : gastr.f param_def
	$(RMC)  $@
	$(F90) -O3 -ip  -c gastr.f
	$(MVC) gastr$O $@

$(OD)/gatemps$O : gatemps.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c gatemps.f
	$(MVC) gatemps$O $@

$(OD)/get_edest$O : get_edest.f common.main \
                   param_def $(OD)/mod_main$O \
                   $(OD)/mod_eleblocks$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c get_edest.f
	$(MVC) get_edest$O $@

$(OD)/getgpts$O : getgpts.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c getgpts.f
	$(MVC) getgpts$O $@

$(OD)/gp_temp_eps$O : gp_temp_eps.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c gp_temp_eps.f
	$(MVC) gp_temp_eps$O $@

$(OD)/gpifv1$O : gpifv1.f param_def \
                   include_sig_up $(OD)/mod_segmental_curves$O
	$(RMC)  $@
	$(F90) -O3 -ip  -c gpifv1.f
	$(MVC) gpifv1$O $@

$(OD)/gpmas1$O : gpmas1.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c gpmas1.f
	$(MVC) gpmas1$O $@

$(OD)/gptns1$O : gptns1.f include_tan_ek $(OD)/mod_crystals$O\
                   param_def $(OD)/mod_main$O $(OD)/mod_eleblocks$O
	$(RMC)  $@
	$(F90) -O3 -ip -c gptns1.f 
	$(MVC) gptns1$O $@

$(OD)/growth_loads$O : growth_loads.f $(OD)/mod_crack_growth$O \
                   param_def common.main \
                   $(OD)/mod_damage$O $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c growth_loads.f
	$(MVC) growth_loads$O $@

$(OD)/gtlsn1$O : gtlsn1.f param_def
	$(RMC)  $@
	$(F90) -O3 -ip  -c gtlsn1.f
	$(MVC) gtlsn1$O $@

$(OD)/gtmat1$O : gtmat1.f param_def
	$(RMC)  $@
	$(F90) -O3 -ip  -c gtmat1.f
	$(MVC) gtmat1$O $@
#
#       --- h --
#
$(OD)/hypre$O : iterative_sparse_hypre.f  common.main param_def \
                $(OD)/mod_hypre$O $(OD)/distributed_assembly$O \
                $(OD)/mod_local_stiffness$O
	$(RMC) $@
	$(F90) -O2 -ip -c iterative_sparse_hypre.f -o hypre$O
	$(MVC) hypre$O $@
#
#       --- i --
#
$(OD)/inclmass$O : inclmass.f $(OD)/mod_eleblocks$O \
                   common.main param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c inclmass.f
	$(MVC) inclmass$O $@

$(OD)/inalpha$O : inalpha.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c inalpha.f
	$(MVC) inalpha$O $@

$(OD)/incomp$O : incomp.f common.main \
                   param_def $(OD)/mod_main$O \
                   $(OD)/mod_mpc$O $(OD)/mod_contact$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c incomp.f
	$(MVC) incomp$O $@

$(OD)/incon$O : incon.f common.main \
                   param_def $(OD)/mod_damage$O \
                   $(OD)/mod_main$O $(OD)/mod_mpc$O
	$(RMC)  $@
	$(F90) -O1 -ip  -c incon.f
	$(MVC) incon$O $@

$(OD)/insurf$O : insurf.f $(OD)/mod_mpc$O \
                   common.main param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c insurf.f
	$(MVC) insurf$O $@

$(OD)/intied$O : intied.f $(OD)/mod_mpc$O \
                   common.main param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c intied.f
	$(MVC) intied$O $@

$(OD)/incontact$O : incontact.f common.main \
                   param_def $(OD)/mod_contact$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c incontact.f
	$(MVC) incontact$O $@

$(OD)/incoor$O : incoor.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c incoor.f
	$(MVC) incoor$O $@

$(OD)/incrystal$O : incrystal.f param_def common.main \
                    $(OD)/mod_crystals$O $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c incrystal.f
	$(MVC) incrystal$O $@

$(OD)/incrack$O : incrack.f $(OD)/mod_crack_growth$O \
                   $(OD)/mod_damage$O common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c incrack.f
	$(MVC) incrack$O $@

$(OD)/infgm$O : infgm.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c infgm.f
	$(MVC) infgm$O $@

$(OD)/incurv$O : incurv.f param_def \
                   $(OD)/mod_segmental_curves$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c incurv.f
	$(MVC) incurv$O $@

$(OD)/indom$O : indom.f common.main \
                   param_def $(OD)/mod_jintegral$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c indom.f
	$(MVC) indom$O $@

$(OD)/indypm$O : indypm.f param_def common.main \
                 $(OD)/mod_performance$O \
                 $(OD)/mod_main$O $(OD)/mod_damage$O \
                 $(OD)/mod_mpc$O  $(OD)/mod_hypre$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c indypm.f
	$(MVC) indypm$O $@

$(OD)/inelbk$O : inelbk.f param_def common.main
	$(RMC)  $@
	$(F90) -O2 -ip  -c inelbk.f
	$(MVC) inelbk$O $@

$(OD)/inelem$O : inelem.f param_def common.main
	$(RMC)  $@
	$(F90) -O2 -ip  -c inelem.f
	$(MVC) inelem$O $@

$(OD)/infile$O : infile.f $(OD)/mod_file_info$O \
                   common.main param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c infile.f
	$(MVC) infile$O $@

$(OD)/inicon$O : inicon.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c inicon.f
	$(MVC) inicon$O $@

$(OD)/ininc$O : ininc.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c ininc.f
	$(MVC) ininc$O $@

$(OD)/initst$O : initst.f param_def common.main \
                 $(OD)/mod_damage$O $(OD)/mod_file_info$O \
                 $(OD)/mod_contact$O $(OD)/mod_crystals$O \
                 $(OD)/mod_elem_load$O $(OD)/mod_stiffness$O \
                 $(OD)/mod_main$O \
                 $(OD)/mod_segmental_curves$O $(OD)/mod_hypre$O \
                 $(OD)/mod_performance$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c initst.f
	$(MVC) initst$O $@

$(OD)/initdm$O : initdm.f $(OD)/mod_jintegral$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c initdm.f
	$(MVC) initdm$O $@

$(OD)/inlod$O : inlod.f $(OD)/mod_elem_load$O \
                   common.main param_def \
                   $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c inlod.f
	$(MVC) inlod$O $@

$(OD)/inmat$O : inmat.f common.main \
                   param_def $(OD)/mod_segmental_curves$O \
                   $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c inmat.f
	$(MVC) inmat$O $@

$(OD)/innum$O : innum.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c innum.f
	$(MVC) innum$O $@

$(OD)/intab$O : intab.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c intab.f
	$(MVC) intab$O $@

$(OD)/iodevn$O : iodevn.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c iodevn.f
	$(MVC) iodevn$O $@
#
#       --- k --
#
$(OD)/kg1$O : kg1.f param_def
	$(RMC)  $@
	$(F90) -O3 -ip  -c kg1.f
	$(MVC) kg1$O $@

$(OD)/kill_element$O : kill_element.f common.main \
                   param_def $(OD)/mod_main$O \
                   $(OD)/mod_eleblocks$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c kill_element.f
	$(MVC) kill_element$O $@
#
#       --- l --
#
$(OD)/license$O : license.f param_def \
                   common.main
	$(RMC)  $@
	$(F90) -O2 -ip  -c license.f
	$(MVC) license$O $@

$(OD)/linked_list$O : linked_list.f $(OD)/mod_crack_growth$O \
                   param_def common.main \
                   $(OD)/mod_damage$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c linked_list.f
	$(MVC) linked_list$O $@

#
#       --- m --
#

$(OD)/main_program$O : main_program.f param_def common.main $(OD)/mod_main$O \
                 $(OD)/mod_file_info$O \
                 $(OD)/mod_performance$O
	$(RMC)  $@
	$(F90) -D$(OSNAME) -DBuild=$(BUILDNUM) -fpp -O2 -ip  -c main_program.f
	$(MVC) main_program$O $@

$(OD)/mem_allocate$O : mem_allocate.f $(OD)/mod_main$O \
                   $(OD)/mod_contact$O \
                   common.main param_def \
                   $(OD)/mod_mpi_lnpcg$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c mem_allocate.f
	$(MVC) mem_allocate$O $@

$(OD)/mm01$O : mm01.f param_def
	$(RMC)  $@
	$(F90) -O3 -ip   -c mm01.f
	$(MVC) mm01$O $@

$(OD)/mm02$O : mm02.f param_def
	$(RMC)  $@
	$(F90) -O3 -ip  -c mm02.f
	$(MVC) mm02$O $@

$(OD)/mm03$O : mm03.f
	$(RMC)  $@
	$(F90) -O3 -ip -c mm03.f
	$(MVC) mm03$O $@

$(OD)/mm04$O : mm04.f param_def 
	$(RMC)  $@
	$(F90) -O3 -ip -c mm04.f
	$(MVC) mm04$O $@

$(OD)/mm05$O : mm05.f
	$(RMC)  $@
	$(F90) -O3 -ip  -c mm05.f
	$(MVC) mm05$O $@

$(OD)/mm06$O : mm06.f
	$(RMC)  $@
	$(F90) -O3 -ip   -c mm06.f 
	$(MVC) mm06$O $@

$(OD)/mm07$O : mm07.f
	$(RMC)  $@
	$(F90) -O3 -ip  -c mm07.f
	$(MVC) mm07$O $@

$(OD)/mm08$O : mm08.f
	$(RMC)  $@
	$(F90) -O3 -ip  -c mm08.f
	$(MVC) mm08$O $@

$(OD)/mm09$O : mm09.f
	$(RMC)  $@
	$(F90) -O3 -ip  -c mm09.f
	$(MVC) mm09$O $@

$(OD)/mm10_a$O : mm10_a.f $(OD)/mod_crystals$O param_def include_mm10
	$(RMC)  $@
	$(F90) -O3 -ip  -c mm10_a.f
	$(MVC) mm10_a$O $@

$(OD)/mm10_b$O : mm10_b.f $(OD)/mod_crystals$O param_def include_mm10
	$(RMC)  $@
	$(F90) -O3  -inline-level=1 -ip  -c mm10_b.f
	$(MVC) mm10_b$O $@

$(OD)/mm10_c$O : mm10_c.f $(OD)/mod_crystals$O param_def include_mm10
	$(RMC)  $@
	$(F90) -O3 -ip   -c mm10_c.f
	$(MVC) mm10_c$O $@

$(OD)/mm10_d$O : mm10_d.f $(OD)/mod_crystals$O param_def include_mm10
	$(RMC)  $@
	$(F90) -O3 -ip   -c mm10_d.f
	$(MVC) mm10_d$O $@
	
$(OD)/mm10_e$O : mm10_e.f $(OD)/mod_crystals$O param_def include_mm10
	$(RMC)  $@
	$(F90) -O3 -ip   -c mm10_e.f
	$(MVC) mm10_e$O $@

$(OD)/mm10_f$O : mm10_f.f $(OD)/mod_crystals$O param_def include_mm10
	$(RMC)  $@
	$(F90) -O3 -ip -c mm10_f.f
	$(MVC) mm10_f$O $@
	
$(OD)/mm11$O : mm11.f $(OD)/mod_crystals$O param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c mm11.f
	$(MVC) mm11$O $@

$(OD)/mmnoopt$O : mmnoopt.f $(OD)/mod_crystals$O param_def
	$(RMC)  $@
	$(F90) -O1 -ip  -c mmnoopt.f
	$(MVC) mmnoopt$O $@

$(OD)/mnralg$O : mnralg.f param_def common.main \
                 $(OD)/mod_adaptive_steps$O \
                 $(OD)/mod_stiffness$O $(OD)/mod_mpc$O $(OD)/mod_performance$O \
                 $(OD)/mod_main$O $(OD)/mod_eleblocks$O \
                 $(OD)/mod_hypre$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c mnralg.f
	$(MVC) mnralg$O $@

$(OD)/mpc_modify$O : mpc_modify.f $(OD)/mod_stiffness$O \
                   $(OD)/mod_mpc$O
	$(RMC)  $@
	$(F90) -O3 -ip  -c mpc_modify.f
	$(MVC) mpc_modify$O $@

$(OD)/mpi_code$O : mpi_code.f param_def common.main $(OD)/mod_damage$O \
                   $(OD)/mod_adaptive_steps$O $(OD)/mod_mpi_lnpcg$O \
                   $(OD)/mod_main$O $(OD)/mod_eleblocks$O \
                   $(OD)/mod_segmental_curves$O \
                   $(OD)/mod_crack_growth$O $(OD)/mod_jintegral$O \
                   $(OD)/mod_contact$O $(OD)/mod_crystals$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c mpi_code.f
	$(MVC) mpi_code$O $@

$(OD)/mpi_handle_slaves$O : mpi_handle_slaves.f $(OD)/mod_damage$O \
                   $(OD)/mod_adaptive_steps$O $(OD)/mod_mpi_lnpcg$O \
                   $(OD)/mod_main$O $(OD)/mod_eleblocks$O \
                   $(OD)/mod_segmental_curves$O \
                   $(OD)/mod_crack_growth$O $(OD)/mod_jintegral$O \
                   common.main param_def \
                   $(OD)/mod_contact$O mpi_code.f
	$(RMC)  $@
	$(F90) -O0 -ip  -c mpi_handle_slaves.f
	$(MVC) mpi_handle_slaves$O $@

#
#       --- n --
#
$(OD)/name_strip$O : name_strip.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c name_strip.f
	$(MVC) name_strip$O $@

$(OD)/ndpts1$O : ndpts1.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c ndpts1.f
	$(MVC) ndpts1$O $@

#
#       --- o --
#
$(OD)/ou3dpr$O : ou3dpr.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c ou3dpr.f
	$(MVC) ou3dpr$O $@

$(OD)/oualph$O : oualph.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c oualph.f
	$(MVC) oualph$O $@

$(OD)/ouddpa$O : ouddpa.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouddpa.f
	$(MVC) ouddpa$O $@

$(OD)/oudriv$O : oudriv.f $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c oudriv.f
	$(MVC) oudriv$O $@

$(OD)/oudups$O : oudups.f common.main \
                   param_def $(OD)/mod_main$O \
                   $(OD)/mod_eleblocks$O $(OD)/mod_elblk$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c oudups.f
	$(MVC) oudups$O $@

$(OD)/oudva$O : oudva.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c oudva.f
	$(MVC) oudva$O $@

$(OD)/oueff1$O : oueff1.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c oueff1.f
	$(MVC) oueff1$O $@

$(OD)/ouext1$O : ouext1.f $(OD)/mod_elblk$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouext1.f
	$(MVC) ouext1$O $@

$(OD)/ouext2$O : ouext2.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouext2.f
	$(MVC) ouext2$O $@

$(OD)/ouflnm$O : ouflnm.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouflnm.f
	$(MVC) ouflnm$O $@

$(OD)/ougts1$O : ougts1.f param_def \
                   $(OD)/mod_elblk$O
	$(RMC)  $@
	$(F90) -O3 -ip  -c ougts1.f
	$(MVC) ougts1$O $@

$(OD)/ouhdva$O : ouhdva.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouhdva.f
	$(MVC) ouhdva$O $@

$(OD)/ouhel$O : ouhel.f common.main \
                   param_def $(OD)/mod_main$O \
                   $(OD)/mod_elblk$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouhel.f
	$(MVC) ouhel$O $@

$(OD)/ouhnod$O : ouhnod.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouhnod.f
	$(MVC) ouhnod$O $@

$(OD)/ouhprn$O : ouhprn.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouhprn.f
	$(MVC) ouhprn$O $@

$(OD)/ouhrks$O : ouhrks.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouhrks.f
	$(MVC) ouhrks$O $@

$(OD)/ouhstr$O : ouhstr.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouhstr.f
	$(MVC) ouhstr$O $@

$(OD)/ouinv1$O : ouinv1.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouinv1.f
	$(MVC) ouinv1$O $@

$(OD)/oulbdd$O : oulbdd.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c oulbdd.f
	$(MVC) oulbdd$O $@

$(OD)/oulbir$O : oulbir.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c oulbir.f
	$(MVC) oulbir$O $@

$(OD)/oulbst$O : oulbst.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c oulbst.f
	$(MVC) oulbst$O $@

$(OD)/oulbta$O : oulbta.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c oulbta.f
	$(MVC) oulbta$O $@

$(OD)/oulg1$O : oulg1.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c oulg1.f
	$(MVC) oulg1$O $@

$(OD)/oulrlh$O : oulrlh.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c oulrlh.f
	$(MVC) oulrlh$O $@

$(OD)/oumkcv$O : oumkcv.f param_def \
                   $(OD)/mod_elblk$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c oumkcv.f
	$(MVC) oumkcv$O $@

$(OD)/oundhd$O : oundhd.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c oundhd.f
	$(MVC) oundhd$O $@

$(OD)/ounds1$O : ounds1.f param_def \
                   $(OD)/mod_elblk$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c ounds1.f
	$(MVC) ounds1$O $@

$(OD)/ouneut$O : ouneut.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouneut.f
	$(MVC) ouneut$O $@

$(OD)/ouocdd$O : ouocdd.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouocdd.f
	$(MVC) ouocdd$O $@

$(OD)/ouocst$O : ouocst.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouocst.f
	$(MVC) ouocst$O $@

$(OD)/ouocst_elem$O : ouocst_elem.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouocst_elem.f
	$(MVC) ouocst_elem$O $@

$(OD)/oupads$O : oupads.f param_def \
                   $(OD)/mod_elblk$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c oupads.f
	$(MVC) oupads$O $@

$(OD)/oupele$O : oupele.f param_def \
                   $(OD)/mod_elblk$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c oupele.f
	$(MVC) oupele$O $@

$(OD)/ouparm$O : ouparm.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouparm.f
	$(MVC) ouparm$O $@

$(OD)/oupdva$O : oupdva.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c oupdva.f
	$(MVC) oupdva$O $@

$(OD)/oupri1$O : oupri1.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c oupri1.f
	$(MVC) oupri1$O $@

$(OD)/ouprks$O : ouprks.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouprks.f
	$(MVC) ouprks$O $@

$(OD)/oupstr_elem$O : oupstr_elem.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c oupstr_elem.f
	$(MVC) oupstr_elem$O $@

$(OD)/oupstr_node$O : oupstr_node.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c oupstr_node.f
	$(MVC) oupstr_node$O $@

$(OD)/oures$O : oures.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c oures.f
	$(MVC) oures$O $@

$(OD)/ourlhd$O : ourlhd.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c ourlhd.f
	$(MVC) ourlhd$O $@

$(OD)/ousthd$O : ousthd.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c ousthd.f
	$(MVC) ousthd$O $@

$(OD)/oustpa$O : oustpa.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c oustpa.f
	$(MVC) oustpa$O $@

$(OD)/oust_elem$O : oust_elem.f param_def common.main \
                   $(OD)/mod_crystals$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c oust_elem.f
	$(MVC) oust_elem$O $@

$(OD)/oustates$O : oustates.f common.main \
                   param_def $(OD)/mod_main$O \
                   $(OD)/mod_eleblocks$O $(OD)/mod_elblk$O \
                   $(OD)/mod_crystals$O
	$(RMC)  $@
	$(F90) -O3 -ip  -c oustates.f
	$(MVC) oustates$O $@

$(OD)/oustr$O : oustr.f common.main \
                   param_def $(OD)/mod_elblk$O \
                   $(OD)/mod_eleblocks$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c oustr.f
	$(MVC) oustr$O $@

$(OD)/outime$O : outime.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c outime.f
	$(MVC) outime$O $@

$(OD)/outfil$O : outfil.f $(OD)/mod_file_info$O \
                   common.main param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c outfil.f
	$(MVC) outfil$O $@

$(OD)/ouyld1$O : ouyld1.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c ouyld1.f
	$(MVC) ouyld1$O $@
#
#       --- p --
#
$(OD)/packets$O : packets.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c packets.f
	$(MVC) packets$O $@

$(OD)/pconvert$O : pconvert.f common.main \
                   param_def $(OD)/mod_main$O \
                          $(OD)/mod_pconvert$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c pconvert.f
	$(MVC) pconvert$O $@

$(OD)/princ_inv_strain$O : princ_inv_strain.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c princ_inv_strain.f
	$(MVC) princ_inv_strain$O $@

$(OD)/princ_inv_stress$O : princ_inv_stress.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c princ_inv_stress.f
	$(MVC) princ_inv_stress$O $@

$(OD)/princ_strain$O : princ_strain.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c princ_strain.f
	$(MVC) princ_strain$O $@

$(OD)/princ_stress$O : princ_stress.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c princ_stress.f
	$(MVC) princ_stress$O $@

$(OD)/prcsel$O : prcsel.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c prcsel.f
	$(MVC) prcsel$O $@
#
#       --- q --
#
$(OD)/qmply1$O : qmply1.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c qmply1.f
	$(MVC) qmply1$O $@
#
#       --- r --
#
$(OD)/read_damage$O : read_damage.f $(OD)/mod_crack_growth$O \
                   $(OD)/mod_damage$O common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c read_damage.f
	$(MVC) read_damage$O $@

$(OD)/read_write$O : read_write.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c read_write.f
	$(MVC) read_write$O $@

$(OD)/reopen$O : reopen.f param_def common.main $(OD)/mod_eleblocks$O \
                 $(OD)/mod_damage$O $(OD)/mod_performance$O \
                 $(OD)/mod_elem_load$O param_def $(OD)/mod_hypre$O \
                 $(OD)/mod_main$O $(OD)/mod_mpc$O \
                 $(OD)/mod_contact$O  \
                 $(OD)/mod_segmental_curves$O $(OD)/mod_crystals$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c reopen.f
	$(MVC) reopen$O $@

$(OD)/rknifv$O : rknifv.f common.main \
                   param_def include_sig_up \
                   $(OD)/mod_crack_growth$O $(OD)/mod_segmental_curves$O \
                   $(OD)/mod_damage$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c rknifv.f
	$(MVC) rknifv$O $@

$(OD)/rknmas$O : rknmas.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c rknmas.f
	$(MVC) rknmas$O $@

$(OD)/rknstr$O : rknstr.f param_def \
                   include_sig_up $(OD)/mod_segmental_curves$O \
                   $(OD)/mod_main$O $(OD)/mod_crystals$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c rknstr.f
	$(MVC) rknstr$O $@

$(OD)/rktstf$O : rktstf.f param_def include_tan_ek \
                 $(OD)/mod_main$O $(OD)/mod_crystals$O
	$(RMC)  $@
	$(F90) -O3  -ip -c rktstf.f
	$(MVC) rktstf$O $@

$(OD)/rplstr$O : rplstr.f common.main \
                   param_def include_sig_up \
                   $(OD)/mod_eleblocks$O $(OD)/mod_segmental_curves$O
	$(RMC)  $@
	$(F90) -O3 -ip  -c rplstr.f 
	$(MVC) rplstr$O $@

$(OD)/rscmp1$O : rscmp1.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c rscmp1.f
	$(MVC) rscmp1$O $@

$(OD)/rstgp1$O : rstgp1.f param_def $(OD)/mod_crystals$O \
                   include_sig_up $(OD)/mod_segmental_curves$O \
                   $(OD)/mod_main$O common.main  $(OD)/mod_eleblocks$O
	$(RMC)  $@
	$(F90) -O3 -ip  -c rstgp1.f
	$(MVC) rstgp1$O $@
#
#       --- s --
#
$(OD)/scan$O : scan.f $(OD)/mod_file_info$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c scan.f
	$(MVC) scan$O $@

$(OD)/scstr$O : scstr.f param_def
	$(RMC)  $@
	$(F90) -O3 -ip  -c scstr.f
	$(MVC) scstr$O $@

$(OD)/setelb$O : setelb.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c setelb.f
	$(MVC) setelb$O $@

$(OD)/setup$O : setup.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c setup.f
	$(MVC) setup$O $@

$(OD)/setup_mat_props$O : setup_mat_props.f param_def \
                   $(OD)/mod_segmental_curves$O common.main
	$(RMC)  $@
	$(F90) -O2 -ip  -c setup_mat_props.f
	$(MVC) setup_mat_props$O $@

$(OD)/shapef$O : shapef.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c shapef.f
	$(MVC) shapef$O $@

$(OD)/splunj$O : splunj.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c splunj.f
	$(MVC) splunj$O $@

$(OD)/star_com$O : star_com.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c star_com.f
	$(MVC) star_com$O $@

$(OD)/statistics$O : statistics.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c statistics.f
	$(MVC) statistics$O $@

$(OD)/steptime$O : steptime.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c steptime.f
	$(MVC) steptime$O $@

$(OD)/store$O : store.f param_def common.main $(OD)/mod_eleblocks$O \
                   $(OD)/mod_crack_growth$O $(OD)/mod_performance$O  \
                   $(OD)/mod_damage$O $(OD)/mod_contact$O \
                   $(OD)/mod_elem_load$O $(OD)/mod_crystals$O \
                   $(OD)/mod_main$O \
                   $(OD)/mod_hypre$O \
                   $(OD)/mod_segmental_curves$O $(OD)/mod_mpc$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c store.f
	$(MVC) store$O $@

$(OD)/store_ifv$O : store_ifv.f $(OD)/mod_crack_growth$O \
                   common.main $(OD)/mod_damage$O \
                   param_def $(OD)/mod_main$O \
                   $(OD)/mod_eleblocks$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c store_ifv.f
	$(MVC) store_ifv$O $@

$(OD)/stpdrv$O : stpdrv.f common.main include_usr_parm \
                   param_def $(OD)/mod_damage$O \
                   $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c stpdrv.f
	$(MVC) stpdrv$O $@
#
#       --- t --
#
$(OD)/tanstf$O : tanstf.f $(OD)/mod_eleblocks$O $(OD)/mod_segmental_curves$O \
                   $(OD)/mod_crack_growth$O include_tan_ek \
                   $(OD)/mod_damage$O common.main \
                   param_def $(OD)/mod_main$O \
                   $(OD)/mod_mpi_lnpcg$O
	$(RMC)  $@
	$(F90) -O3 -ip -c tanstf.f
	$(MVC) tanstf$O $@

$(OD)/temsto$O : temsto.f common.main \
                   param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c temsto.f
	$(MVC) temsto$O $@

$(OD)/thyme$O : thyme.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c thyme.f
	$(MVC) thyme$O $@

$(OD)/tied_mesh$O : tied_mesh.f param_def \
                   $(OD)/mod_mpc$O
	$(RMC)  $@
	$(F90) -O1 -ip  -c tied_mesh.f
	$(MVC) tied_mesh$O $@

$(OD)/tilde$O : tilde.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c tilde.f
	$(MVC) tilde$O $@

$(OD)/trnmtx$O : trnmtx.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c trnmtx.f
	$(MVC) trnmtx$O $@

$(OD)/trnvec$O : trnvec.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c trnvec.f
	$(MVC) trnvec$O $@

#
#       --- u --
#

$(OD)/update$O : update.f common.main \
                   param_def $(OD)/mod_main$O \
                   $(OD)/mod_eleblocks$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c update.f
	$(MVC) update$O $@

$(OD)/update_node_elecnt$O : update_node_elecnt.f $(OD)/mod_crack_growth$O \
                   common.main param_def \
                   $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c update_node_elecnt.f
	$(MVC) update_node_elecnt$O $@

$(OD)/user_list$O : user_list.f common.main \
                  param_def $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O3 -ip  -c user_list.f
	$(MVC) user_list$O $@

$(OD)/uppbar$O : uppbar.f common.main \
                   param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c uppbar.f
	$(MVC) uppbar$O $@

$(OD)/upres$O : upres.f common.main \
                   $(OD)/mod_contact$O param_def \
                   $(OD)/mod_mpc$O $(OD)/mod_stiffness$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c upres.f
	$(MVC) upres$O $@

#
#       ---  Abaqus compatible user routines ------
#
$(OD)/user_routines_umat$O : user_routines_umat.f $(OD)/user_routines_module$O 
	$(RMC)  $@
	$(F90) -O3 -ip  -inline-level=1  -c user_routines_umat.f
	$(MVC) user_routines_umat$O $@

$(OD)/user_routines_other_Abaqus$O : user_routines_other_Abaqus.f \
                                     $(OD)/user_routines_module$O 
	$(RMC)  $@
	$(F90) -O3 -ip -heap-arrays -c user_routines_other_Abaqus.f
	$(MVC) user_routines_other_Abaqus$O $@

$(OD)/user_routines_warp3d$O : user_routines_warp3d.f  include_usr_parm \
                               $(OD)/user_routines_module$O 
	$(RMC)  $@
	$(F90) -O3 -ip -heap-arrays  -c user_routines_warp3d.f
	$(MVC) user_routines_warp3d$O $@

$(OD)/user_routines_support$O : user_routines_support.f  include_usr_parm \
                                $(OD)/user_routines_module$O 
	$(RMC)  $@
	$(F90) -O3 -ip -c user_routines_support.f
	$(MVC) user_routines_support$O $@

#
#       --- v --
#
$(OD)/vol_avg$O : vol_avg.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c vol_avg.f
	$(MVC) vol_avg$O $@

$(OD)/vol_terms$O : vol_terms.f
	$(RMC)  $@
	$(F90) -O3 -ip  -c vol_terms.f
	$(MVC) vol_terms$O $@
#
#       --- w --
#
#       --- y --
#
$(OD)/yield_function$O : yield_function.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c yield_function.f
	$(MVC) yield_function$O $@
#
#       --- z --
#
$(OD)/zero_vector$O : zero_vector.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c zero_vector.f
	$(MVC) zero_vector$O $@

$(OD)/zero_vol$O : zero_vol.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c zero_vol.f
	$(MVC) zero_vol$O $@
#
#    ------   all module files ----------------
#
$(ODIR)/mod_adaptive_steps$O : mod_adaptive_steps.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_adaptive_steps.f
	$(MVC) mod_adaptive_steps$O $@

$(ODIR)/mod_contact$O : mod_contact.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_contact.f
	$(MVC) mod_contact$O $@

$(ODIR)/mod_crack_growth$O : mod_crack_growth.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_crack_growth.f
	$(MVC) mod_crack_growth$O $@

$(ODIR)/mod_crystals$O : mod_crystals.f param_def \
                         common.main $(OD)/mod_main$O
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_crystals.f
	$(MVC) mod_crystals$O $@

$(ODIR)/mod_damage$O : mod_damage.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_damage.f
	$(MVC) mod_damage$O $@

$(ODIR)/mod_elblk$O : mod_elblk.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_elblk.f
	$(MVC) mod_elblk$O $@

$(ODIR)/mod_eleblocks$O : mod_eleblocks.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_eleblocks.f
	$(MVC) mod_eleblocks$O $@

$(ODIR)/mod_elem_load$O : mod_elem_load.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_elem_load.f
	$(MVC) mod_elem_load$O $@

$(ODIR)/mod_file_info$O : mod_file_info.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_file_info.f
	$(MVC) mod_file_info$O $@

$(ODIR)/mod_hypre$O : mod_hypre.f \
		                common.main param_def
	$(RMC) $@
	$(F90) -O2 -ip -c mod_hypre.f
	$(MVC) mod_hypre$O $@

$(ODIR)/mod_jintegral$O : mod_jintegral.f param_def
	$(RMC)  $@
	$(F90) -O0 -ip  -c mod_jintegral.f
	$(MVC) mod_jintegral$O $@

$(ODIR)/mod_local_stiffness$O : mod_local_stiffness.f \
		                              common.main param_def
	$(RMC) $@
	$(F90) -O2 -ip -c mod_local_stiffness.f
	$(MVC) mod_local_stiffness$O $@
	
$(ODIR)/mod_main$O : mod_main.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_main.f
	$(MVC) mod_main$O $@

$(ODIR)/mod_mpc$O : mod_mpc.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_mpc.f
	$(MVC) mod_mpc$O $@

$(ODIR)/mod_mpi_lnpcg$O : mod_mpi_lnpcg.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_mpi_lnpcg.f
	$(MVC) mod_mpi_lnpcg$O $@

$(ODIR)/mod_pconvert$O : mod_pconvert.f
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_pconvert.f
	$(MVC) mod_pconvert$O $@

$(ODIR)/mod_performance$O : mod_performance.f
	$(RMC) $@
	$(F90) -O2 -ip -c mod_performance.f
	$(MVC) mod_performance$O $@

$(ODIR)/mod_segmental_curves$O : mod_segmental_curves.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_segmental_curves.f
	$(MVC) mod_segmental_curves$O $@

$(ODIR)/mod_stiffness$O : mod_stiffness.f param_def
	$(RMC)  $@
	$(F90) -O2 -ip  -c mod_stiffness.f
	$(MVC) mod_stiffness$O $@

$(ODIR)/user_routines_module$O : user_routines_module.f 
	$(RMC)  $@
	$(F90) -O2 -ip  -c user_routines_module.f
	$(MVC) user_routines_module$O $@

