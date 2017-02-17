c *******************************************************************
c *                                                                 *
c *        material model # 9 -- <available>                        *  
c *                                                                 *
c *******************************************************************
c
c
      subroutine mm12( 
     &  step, iter, felem, gpn, mxvl, hist_size, nstrs, nstrn, span,
     &  iout, signal_flag, adaptive_possible, cut_step_size_now,
     &  model_ptr,
     &  time_n, time_np1, temp_n, temp_np1,
     &  strain_n, strain_np1, stress_n, stress_np1,
     &  hist_n, hist_np1)
      use iso_c_binding
      implicit none
      include "neml_interface.f"
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &  step, iter, felem, gpn, mxvl, hist_size, span,
     &  iout, nstrs, nstrn
c
      logical
     &   signal_flag, adaptive_possible, cut_step_size_now,
     &   killed_list_vec(mxvl)
c     
c     
      double precision
     & stress_n(mxvl,nstrs), 
     & stress_np1(mxvl,nstrs), deps(mxvl,nstrn),
     & time_n, time_np1, temp_n(mxvl), temp_np1(mxvl),
     & strain_n(mxvl,nstrs), strain_np1(mxvl,nstrs),
     & hist_n(span,hist_size),
     & hist_np1(span,hist_size)
c
c
      type(c_ptr) :: model_ptr


c
c               description of parameters
c               -------------------------
c
c     step              : current load step number
c     iter              : current newton iteration number
c     felem             : first element of the current block
c     gpn               : gauss point number being processed for block
c     mxvl              : maximum no. elements per block
c     hist_size         : number of history words per gauss point
c     nstrs             : number of stress terms (6 + extras)
c     nstrn             : number of incremental strain terms (6)
c     span              : number of elements in current block
c     iout              : write messates to this device number
c     signal_flag       : user wants notification messages for key
c                         events in material response
c     adaptive_possible : .true. if the material model may request
c                         immediate reduction of global load step size.
c                         no stress-histroy update required
c (*) cut_step_size_now : set .true. if material model wants immediate
c                         reduction of global load step size.
c                         no stress-history update required
c     model_ptr         : pointer to external model
c     time_n            : previous time
c     time_np1          : next time
c     temp_n            : previous temperatures
c     temp_np1          : next temperatures
c     strain_n          : previous mechanical strains
c     strain_np1        : next mechanical strains
c (**)stress_n          : stresses at start of load step (n) for all
c                         elements in block for this gauss point
c (*) stress_np1        : stresses at end of load step (n+1) for all
c                         elements in block for this gauss point
c (**)hist_n            : history values at start of load step (n) for all
c                         elements in block for this gauss point
c (*) hist_np1          : history values at end of load step (n+1) for all
c                         elements in block for this gauss point

c     
c    (*)  values to be updated by this material model
c    (**) needs to be initialized on step 1
c     
c   Input ordering, standard Voigt notation:
c
c   strain ordering:
c     deps-xx, deps-yy, deps-zz, gamma-xy, gamma-yz, gamma-xz
c
c   stress ordering (at n and n+1):
c     (1) sig-xx
c     (2) sig-yy
c     (3) sig-zz
c     (4) tau-xy
c     (5) tau-yz
c     (6) tau-xz
c     (7) total work density
c     (8) total plastic work density
c     (9) total plastic strain 
c
c
c   NEML ordering, Mandel notation:
c     xx yy zz sqrt(2) yz sqrt(2) xz sqrt(2) xy
c
c
c   Also need to:
c     1) store tangent in the first 36 slots in hist_np1
c     2) initialize history on step 1
c     3) transpose the row major tangent
c
c   Remember we're doing this on a block
c
c     Locals
      integer :: i, j, k, ier, nhist, nstore
      double precision :: l_stress_n(6), l_stress_np1(6), l_strain_n(6),
     &                    l_strain_np1(6), l_tangent(6,6),
     &                    l_full_tangent(6,6)
      double precision, allocatable, dimension(:) :: l_hist_n,
     &                                               l_hist_np1
      double precision :: l_temp_n, l_temp_np1, l_time_n, l_time_np1,
     &                        l_u_n, l_u_np1, l_p_n, l_p_np1
      double precision :: vm_mult_s(6), vm_mult_e(6)
      integer :: vm_map(6)

c     Setup history sizes
      nhist = nstore_nemlmodel(model_ptr)
      nstore = nhist + 36

      allocate(l_hist_n(nhist))
      allocate(l_hist_np1(nhist))

c     Setup the Voigt -> Mandel map
      do i=1,3
            vm_map(i) = i
            vm_mult_s(i) = 1.0
            vm_mult_e(i) = 1.0
      end do
      vm_map(4) = 5
      vm_mult_s(4) = sqrt(2.0)
      vm_mult_e(4) = sqrt(2.0)/2.0

      vm_map(5) = 6
      vm_mult_s(5) = sqrt(2.0)
      vm_mult_e(5) = sqrt(2.0)/2.0

      vm_map(6) = 4
      vm_mult_s(6) = sqrt(2.0)
      vm_mult_e(6) = sqrt(2.0)/2.0

c      For each entry in the block
      do i=1,span
            ! Read into our local structures
            l_temp_n   = temp_n(i)
            l_temp_np1 = temp_np1(i)
            do j=1,6
                  l_stress_n(j) = stress_n(i,vm_map(j)) * vm_mult_s(j)
                  l_strain_n(j) = strain_n(i,vm_map(j)) * vm_mult_e(j)
                  l_strain_np1(j) = strain_np1(i,vm_map(j)) 
     &                  * vm_mult_e(j)
            end do

            l_hist_n = hist_n(i,37:nstore)

            l_u_n = stress_n(i,7)
            l_p_n = stress_n(i,8)
            
            ! If the first step, initialize history and stress
            if (step .eq. 1) then
                  l_stress_n(1:6) = 0.0
                  call init_store_nemlmodel(model_ptr, l_hist_n, ier)
                  ! Check for error
                  if (ier .ne. 0) then
                        write(*,*) "Error setting up NEML history"
                        call destroy_nemlmodel(model_ptr, ier)
                        deallocate(l_hist_n)
                        deallocate(l_hist_np1)
                        call die_abort
                  end if
                  
                  ! Not sure if this is actually necessary
                  stress_n(i,1:6) = l_stress_n ! zeros
                  hist_n(i,37:nstore) = l_hist_n
            end if

            ! Call the update
            call update_sd_nemlmodel(model_ptr, l_strain_np1, 
     &            l_strain_n,
     &            l_temp_np1, l_temp_n, time_np1, time_n, l_stress_np1,
     &            l_stress_n, l_hist_np1, l_hist_n, l_tangent,
     &            l_u_np1, l_u_n, l_p_np1, l_p_n, ier)
            if (ier .ne. 0) then
                  write(*,*) "Error updating NEML material"
                  call destroy_nemlmodel(model_ptr, ier)
                  deallocate(l_hist_n)
                  deallocate(l_hist_np1)
                  call die_abort
            end if

            ! Store the updated quantities
            do j=1,6
                  stress_np1(i,vm_map(j)) = l_stress_np1(j) /
     &                  vm_mult_s(j)
                  do k=1,6
                        l_full_tangent(vm_map(j),vm_map(k)) = 
     &                        l_tangent(j,k)  * vm_mult_e(k) / 
     &                        vm_mult_s(j)
                  end do
            end do
            l_full_tangent = transpose(l_full_tangent)
            hist_np1(i,1:36) = reshape(l_full_tangent, (/36/))
            hist_np1(i,37:nstore) = l_hist_np1
            stress_np1(i,7) = l_u_np1
            stress_np1(i,8) = l_p_np1

            ! This is suppose to be the equivalent plastic strain, but
            ! I refuse to provide an interface for it
            stress_np1(i,9) = 0.0
      end do


      ! Deallocate
      deallocate(l_hist_n)
      deallocate(l_hist_np1)
       
      end


c *******************************************************************
c *                                                                 *
c *        material model # 9 -- <available>                        *  
c *                                                                 *
c *******************************************************************
c
c
      subroutine cnst12( 
     &  span, felem, gpn, iter, iout, mxvl, nstrn, 
     &  weight, history_n,
     &  history_np1, stress_np1, dmat, det_jac_block)
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &  span, felem, gpn, iter, iout, mxvl, nstrn
c
      logical
     &   first, killed_list_vec(mxvl)
c     
      double precision
     & weight, history_n(span,*),
     & history_np1(span,*), dmat(mxvl,nstrn,nstrn),
     & det_jac_block(mxvl), stress_np1(mxvl,nstrn)
c
c
c
c               description of parameters
c               -------------------------
c
c     step              : current load step number
c     iter              : current newton iteration number. iter 1 is
c                         for application of the "real" load increment.
c     felem             : first element of the current block
c     gpn               : gauss point number being processed for block
c     first             : logical flag indicating if this
c                         is the very first call to the routine
c                         the linear elastic [d] is returned
c                         whenever first = .true.
c     mxvl              : maximum no. elements per block
c     nstrn             : number of strain-stress components (=6)
c     span              : number of elements in current block
c     iout              : write messages to this device number
c     weight            : integration point weight factor
c     history_n         : history values at start of load step (n) for all
c                         elements in block for this gauss point
c     history_np1       : history values at end of load step (n+1) for all
c                         elements in block for this gauss point
c     det_jac_block     : |J| at this gauss point for each element in
c                         block
c (*) dmat              : 6x6 (symmetric) tangent (consistent) for
c                         this gauss point for each element of block
c                         (see stress ordering below) 
c     
c    (*)  values to be updated by this material model
c
c      
c     All we need to do is move tangent from history_np1 -> right spots
c   
c
c   Note: warp3d expects all dmat[] values to be multiplied by
c         weight * det_jac_block(i), where i = relative
c         element number of block
c
c
c                   local variables
c                   ---------------
c     
      integer i    
c
      do i = 1, span
            dmat(i,1:6,1:6) = reshape(history_np1(i,1:36), (/6,6/)) * 
     &            weight * det_jac_block(i)
      end do

c
      return
      end subroutine
c
c *******************************************************************
c *                                                                 *
c *        material model # 9 -- <available>                        *  
c *                                                                 *
c *           set 3 material model dependent output values          *
c *                                                                 *
c *******************************************************************
c
c
      subroutine oumm12( gpn, mxvl, span, iout, elestr,
     &                   stress, history )
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &  gpn, mxvl, span, iout
c
c     
      double precision
     & stress(mxvl,*), elestr(mxvl,*), history(mxvl,*)
c
c               description of parameters
c               -------------------------
c
c     gpn               : gauss point number being processed for block
c     mxvl              : maximum no. elements per block
c     span              : number of elements in current block
c     iout              : write messages to this device number
c     stress            : current stresses for all
c                         elements in block for this gauss point
c (*) elestr            : stresses to be output for elements
c     history           : current history values for all
c                         elements in block for this gauss point
c     
c    (*)  values to be updated by this material model
c
c
c   stress ordering                elestr ordering
c     (1) sig-xx                   (1) sig-xx
c     (2) sig-yy                   (2) sig-yy
c     (3) sig-zz                   (3) sig-zz
c     (4) tau-xy                   (4) tau-xy
c     (5) tau-yz                   (5) tau-yz        
c     (6) tau-xz                   (6) tau-xz
c     (7) total work density       (7) total work density   
c                                  (8) mises equiv. stress
c                             (*)  (9) mat_val1
c                             (*) (10) mat_val2
c                             (*) (11) mat_val3
c
c  NOTE:  do not modify "stress" array or columns 1-8 of the "elestr"
c         array. only modify columns 9-11 of "elestr". These are the
c         3 "material model" dependent values output in the stress
c         values for a gauss point. See Section 2.12 of manual and
c         description of each material model. The output labels for 
c         columns 9-11 are "c1", "c2", "c3".
c
c
c     
c
       return
       end


c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm12_set_sizes_special            *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 7/1/12  rhd                 *
c     *                                                              *
c     *    called by warp3d for each material model to obtain        *
c     *    various sizes of data for the model                       *
c     *                                                              *
c     *     Used to avoid model call                                 *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm12_set_sizes( info_vector )
      dimension info_vector(*)
c
c        set infor_data
c
c         1        number of history values per integration 
c                  point. Abaqus calls these "statev". Values
c                  double or single precsion based on hardware.
c
c                  we require an extra 36 for the tangent 
c
c    
c         2        number of values in the symmetric part of the 
c                  [D] for each integration point. for solid
c                  elements this is 21, for cohesive elements this 6.
c
c         3        = 0, the material model returns "unrotated"
c                       Cauchy stresses at n+1
c                  = 1, the material model returns the standard
c                       Cauchy stresses at n+1
c
c         4        number of state variables per point to be output
c                  when user requests this type of results
c
      
      
      info_vector(1) = -1
      info_vector(2) = 21
      info_vector(3) = 0
      info_vector(4) = 0
c
      return
      end



c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm12_set_sizes_special            *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 7/1/12  rhd                 *
c     *                                                              *
c     *    called by warp3d for each material model to obtain        *
c     *    various sizes of data for the model                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm12_set_sizes_special( info_vector, model_ptr)
      use iso_c_binding
      include "neml_interface.f"
      integer :: ier, actual_hist
      type(c_ptr) :: model, model_ptr ! Actual compiler error, pad 
      dimension info_vector(*)
c
c        set infor_data
c
c         1        number of history values per integration 
c                  point. Abaqus calls these "statev". Values
c                  double or single precsion based on hardware.
c
c                  we require an extra 36 for the tangent 
c
c    
c         2        number of values in the symmetric part of the 
c                  [D] for each integration point. for solid
c                  elements this is 21, for cohesive elements this 6.
c
c         3        = 0, the material model returns "unrotated"
c                       Cauchy stresses at n+1
c                  = 1, the material model returns the standard
c                       Cauchy stresses at n+1
c
c         4        number of state variables per point to be output
c                  when user requests this type of results
c
      actual_hist = nstore_nemlmodel(model_ptr)

      info_vector(1) = actual_hist + 36
      info_vector(2) = 21
      info_vector(3) = 0
      info_vector(4) = 0

c
      return
      end
c            dummy routines for model not yet supporting 
c            states output
c
c     ****************************************************************
c     *                                                              *
c     *             subroutine mm12_states_values                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 1/3/2015 (rhd))                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm12_states_values( itype, elem_states_output,
     &                                 nrow_states, num_states  )
c
c                       access some global data structures
c
      use elem_block_data, only: history_blocks, history_blk_list
      use main_data, only: elems_to_blocks, cohesive_ele_types
c      
      implicit integer (a-z)
      include 'common.main'
c
c                       parameters
c
      integer :: nrow_states, itype, num_states
      double precision :: elem_states_output(nrow_states,*)
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm12_states_labels                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 1/11/2015 (rhd)                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm12_states_labels( size_state,
     &      num_states, state_labels, state_descriptors, out,
     &      comment_lines, max_comment_lines, num_comment_lines )
      implicit none
c
c                       parameters
c
      integer :: size_state, num_states, out, max_comment_lines,
     &           num_comment_lines
      character(len=8)  :: state_labels(size_state)
      character(len=60) :: state_descriptors(size_state)
      character(len=80) :: comment_lines(max_comment_lines)
c
c                       locals
c
      integer :: i
c
      num_states = 0
      num_comment_lines = 0     
      state_labels(1) = "..."
      state_descriptors(1) = "...."
c      
      return
      end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm12_setup_model                  *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *               last modified : 1/06/2017 (mcm)                *
c
c           Setup an mm12 material model by storing a pointer in
c           the correct data structure.
c     *                                                              *
c     ****************************************************************
c
      subroutine mm12_setup_model(input_file, model_name,
     &            itype, success, ptr)
            use iso_c_binding
            implicit none
            include "neml_interface.f"
            character(len=24), intent(in) :: input_file
            character(len=24), intent(in) :: model_name

            integer, intent(out) :: itype
            logical, intent(out) :: success
            type(c_ptr), intent(out) :: ptr
c
            character(len=25,kind=c_char) :: fname
            character(len=25,kind=c_char) :: mname
            integer :: ier
c
            ! Trim strings and add null-termination
            fname = trim(input_file)//C_NULL_CHAR
            mname = trim(model_name)//C_NULL_CHAR
            
            ! Get the model
            ptr = create_nemlmodel(fname, mname, ier)
            if (ier .ne. 0) then
                  write(*,*) "Error setting up NEML material model"
                  itype = -1
                  success = .false.
                  call die_abort
            end if

            ! Set the information values
            itype = 0
            success = .true.

            return

      end subroutine

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm12_cleanup_model                *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *               last modified : 1/06/2017 (mcm)                *
c
c           clean up an allocated mm12 model
c     *                                                              *
c     ****************************************************************
c
      subroutine mm12_cleanup_model(itype, success, ptr)
            use iso_c_binding
            implicit none
            include "neml_interface.f"
c
            integer, intent(inout) :: itype
            logical, intent(inout) :: success
            type(c_ptr), intent(in) :: ptr
c
            integer :: ier
c
            if (itype .ne. 0) then
                  write(*,*) "Invalid external model in mm12!"
                  call die_abort
            end if
c
            if (.not. success) then
                  write(*,*) "Model already dealloced in mm12!"
                  return
            end if
c           
            call destroy_nemlmodel(ptr, ier)
            if (ier .ne. 0) then
                  write(*,*) "Error destroying NEML material model"
                  call die_abort
            end if
c
            success = .false.
            itype = -1
c
      end subroutine

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm12_setup_alpha                  *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *               last modified : 2/12/2017 (mcm)                *
c
c           get alpha from the model
c     *                                                              *
c     ****************************************************************
c
      subroutine mm12_setup_alpha(span, model_ptr, gp_temps, alpha_vec)
            use iso_c_binding
            implicit none
            include "neml_interface.f"
c
            integer, intent(in) :: span
            type(c_ptr), intent(in) :: model_ptr
            double precision, intent(in) :: gp_temps(span)
            double precision, intent(out) :: alpha_vec(span,6)
c
            integer :: i
            double precision :: alpha
c
            do i=1,span
                  alpha = alpha_nemlmodel(model_ptr, gp_temps(i))
                  alpha_vec(i,1:3) = alpha
                  alpha_vec(i,4:6) = 0.0
            end do
c
      end subroutine
