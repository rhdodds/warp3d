C *******************************************************************
c *                                                                 *
c *        material model # 6 -- adv. gurson model                  *  
c *                                                                 *
c *******************************************************************
c
c
      subroutine mm06( 
     &  step, iter, felem, gpn, mxvl, hist_size, nstrs, nstrn, span,
     &  iout, signal_flag, adaptive_possible, cut_step_size_now,
     &  mm_props, f0_vec, q1_vec, q2_vec, q3_vec, nuc_vec,
     &  nuc_s_n_vec, nuc_e_n_vec, nuc_f_n_vec, 
     &  e_vec, tan_e_vec, nu_vec, sigyld_vec, 
     &  n_power_vec, trial_elas_stress_np1, stress_n, stress_np1,
     &  deps, history_n, history_np1, killed_status_vec)
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &  step, iter, felem, gpn, mxvl, hist_size, span, iout, nstrs,
     &  nstrn
c
      logical
     &   signal_flag, adaptive_possible, cut_step_size_now,
     &   nuc_vec(mxvl), killed_status_vec(mxvl)
c     
#dbl      double precision
#sgl      real
     & mm_props(mxvl,5), e_vec(mxvl), tan_e_vec(mxvl), nu_vec(mxvl),
     & sigyld_vec(mxvl), n_power_vec(mxvl), stress_n(mxvl,nstrs), 
     & stress_np1(mxvl,nstrs), deps(mxvl,nstrn),
     & trial_elas_stress_np1(mxvl,nstrn), history_n(span,hist_size),
     & history_np1(span,hist_size),
     & f0_vec(mxvl), q1_vec(mxvl), q2_vec(mxvl), q3_vec(mxvl), 
     & nuc_s_n_vec(mxvl), nuc_e_n_vec(mxvl), nuc_f_n_vec(mxvl),
     & w0_vec(mxvl),chi0_vec(mxvl),emp_vec(mxvl),
     & fcoal_vec(mxvl),delta(3,3),stress_dev(mxvl,3,3),
     & xcons(mxvl,3,3),plasticity(mxvl),logi_vec(mxvl),
     & coalescence(mxvl)
c
      return
      end




c *******************************************************************
c *                                                                 *
c *        material model # 6 -- adv. gurson model                  *  
c *                                                                 *
c *******************************************************************
c
c
      subroutine cnst6( 
     &  span, felem, gpn, first, iter, iout, mxvl, nstrn,
     &  weight, e_vec, nu_vec, sigyld_vec, n_power_vec, f0_vec,
     &  q1_vec, q2_vec, q3_vec, nuc_vec,nuc_s_n_vec, nuc_e_n_vec, 
     &  nuc_f_n_vec, trial_elas_stress, mm_props, history_n,
     &  history_np1, stress_np1, dmat, det_jac_block,
     &  killed_status_vec)
      implicit none
c
c                   parameter declarations
c                   ----------------------
c
      integer
     &  span, felem, gpn, iter, iout, mxvl, nstrn, hist_size
c
      logical
     &   first, nuc_vec(mxvl),killed_status_vec(mxvl)
c     
#dbl      double precision
#sgl      real
     & weight, mm_props(mxvl,5), e_vec(mxvl), nu_vec(mxvl),
     & sigyld_vec(mxvl),n_power_vec(mxvl),f0_vec(mxvl),
     & trial_elas_stress(mxvl,nstrn), history_n(span,*),
     & history_np1(span,*), dmat(mxvl,nstrn,nstrn),
     & det_jac_block(mxvl), q1_vec(mxvl), q2_vec(mxvl), q3_vec(mxvl),
     & nuc_s_n_vec(mxvl), nuc_e_n_vec(mxvl), nuc_f_n_vec(mxvl),
     & stress_np1(mxvl,nstrn) 


      return
      end

c *******************************************************************
c *                                                                 *
c *        material model # 6 -- adv. gurson model                  *
c *                                                                 *
c *           set 3 material model dependent output values          *
c *                                                                 *
c *******************************************************************
c
c
c
       subroutine oumm06( gpn, mxvl, span, iout, elestr,
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
#dbl      double precision
#sgl      real
     & stress(mxvl,*), elestr(mxvl,*), history(mxvl,*)
c
       return
       end

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm06_set_sizes                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified: 12/14/14  rhd               *
c     *                                                              *
c     *    called by warp3d for each material model to obtain        *
c     *    various sizes of data for the model                       *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm06_set_sizes( info_vector )
      dimension info_vector(*)
c
c        set infor_data
c
c         1        number of history values per integration 
c                  point. Abaqus calles these "statev". Values
c                  double or single precsion based on hardware.
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
      info_vector(1) = 12
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
c     *             subroutine mm06_states_values                    *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 1/3/2015 (rhd))                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm06_states_values( itype, elem_states_output,
     &                                 nrow_states, num_states  )
c
c                       access some global data structures
c
      use elem_block_data, only: history_blocks, history_blk_list
      use main_data, only: elems_to_blocks, cohesive_ele_types
c      
      implicit integer (a-z)
$add common.main
c
c                       parameters
c
      integer :: nrow_states, itype, num_states
#dbl      double precision :: elem_states_output(nrow_states,*)
#sgl      real  :: elem_states_output(nrow_states,*)
c
      return
      end
c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm06_states_labels                *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *               last modified : 1/11/2015 (rhd)                *
c     *                                                              *
c     ****************************************************************
c
      subroutine mm06_states_labels( size_state,
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




