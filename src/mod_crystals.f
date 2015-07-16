c ****************************************************************************
c *                                                                          *
c *    mod_crystals - contains all crystal data, struct for each crystal as  *
c *                   well as the array of crystals and some helper methods  *
c *                                                                          *
c *         written by : mcm                                                 *
c *         last modified : 3/21/12 mcm
c *                                                                          *
c ****************************************************************************
c
      module crystal_data
            implicit integer (a-z)
$add param_def
c                 Crystal array data structures
c
            type :: crystal
                  integer ::  slip_type
c                             1) fcc
c                             2) bcc
c                             3) single
                  integer :: elastic_type
c                             1) isotropic
c                             2) cubic
                  integer :: nslip
                  integer :: h_type
c                             1) voche or voce
c                             2) mts
c                             3) user
                  double precision :: e, nu, mu, harden_n, tau_a,
     &                                tau_hat_y, g_o_y, b, p_v, q_v,
     &                                p_y, q_y, boltz, 
     &                                eps_dot_o_y, t_o,
     &                                theta_o, tau_bar_o,
     &                                tau_hat_v, g_o_v,
     &                                eps_dot_o_v, k_o,
     &                                mu_o, D_o, tau_y, tau_v,
     &                                voche_m !  yes it is spelled wrong
                  double precision :: u1, u2, u3, u4, u5, u6
                  double precision, dimension(6,6) :: elast_stiff,
     &                                                elast_flex
                  double precision, dimension(12,3) :: ni,bi
c
                  logical :: valid

            end type crystal
c
            type(crystal), dimension(max_crystals) :: c_array
c
c           Data structures for mm11 Monte Carlo angles
            integer, allocatable :: mc_array(:,:)

            double precision, allocatable :: angle_input(:,:,:)
            double precision, allocatable :: simple_angles(:,:)
            integer, allocatable :: crystal_input(:,:)
            integer, allocatable :: data_offset(:)
            double precision :: cry_multiplier
            logical :: defined_crystal, srequired
            integer :: nangles

      contains
c                 Initializes crystal num with defaults
            subroutine initialize_new_crystal(num, out)
                  integer, intent(in) :: num, out
c                  
                  if (num .gt. max_crystals) then
                        write (out,*)
     &                   "Error: crystal number exceeds max crystals."
                        call die_gracefully
                  end if                  
                  if (num .lt. 1) then
                        write (out,*)
     &                   "Error: invalid crystal number."
                        call die_gracefully
                  end if
c
                  c_array(num)%slip_type = 1
                  c_array(num)%elastic_type = 1
                  c_array(num)%e = 30000
                  c_array(num)%nu = 0.3
                  c_array(num)%mu = 11538.5
                  c_array(num)%harden_n = 20
                  c_array(num)%tau_a = 0
                  c_array(num)%tau_hat_y = -1.0
                  c_array(num)%g_o_y = -1.0
                  c_array(num)%tau_hat_v = -1.0
                  c_array(num)%g_o_v = -1.0
                  c_array(num)%b = 3.5e-10
                  c_array(num)%p_v = 0.5
                  c_array(num)%q_v = 2
                  c_array(num)%p_y = 0.5
                  c_array(num)%q_y = 2
                  c_array(num)%boltz = 1.3806E-29
                  c_array(num)%eps_dot_o_y = 1.0E10
                  c_array(num)%eps_dot_o_v = 1.0E10
                  c_array(num)%t_o = 273.15
                  c_array(num)%theta_o = 57.7
                  c_array(num)%k_o = 0.0
                  c_array(num)%mu_o = 11538.5
                  c_array(num)%D_o = 0.0
                  c_array(num)%tau_y = 0.0
                  c_array(num)%tau_v = 0.0
                  c_array(num)%voche_m = 1.0

                  c_array(num)%u1 = 0.0
                  c_array(num)%u2 = 0.0
                  c_array(num)%u3 = 0.0
                  c_array(num)%u4 = 0.0
                  c_array(num)%u5 = 0.0
                  c_array(num)%u6 = 0.0

                  c_array(num)%h_type = 1

            end subroutine
c
c                 Inserts "derived" properties like slip system geometry
c                 and elasticity tensors
            subroutine finalize_new_crystal(num, out)
                  integer, intent(in) :: num, out
                  double precision :: f,e,v,u
                  integer :: i,j,info
                  double precision, dimension(600) :: work
                  integer, dimension(600) :: lwork
                  integer, dimension(6) :: ipiv
c
                  if (num .gt. max_crystals) then
                        write (out,*)
     &                   "Error: crystal number exceeds max crystals."
                        call die_gracefully
                  end if                  
                  if (num .lt. 1) then
                        write (out,*)
     &                   "Error: invalid crystal number."
                        call die_gracefully
                  end if
c
c
                  if (c_array(num)%slip_type .eq. 1) then
                        c_array(num)%nslip = 12
                        f = 1/sqrt(2.0D0)
                        
                        c_array(num)%bi(1,1)=0
                        c_array(num)%bi(1,2)=-f
                        c_array(num)%bi(1,3)=f
                        c_array(num)%bi(2,1)=f
                        c_array(num)%bi(2,2)=0
                        c_array(num)%bi(2,3)=-f
                        c_array(num)%bi(3,1)=-f
                        c_array(num)%bi(3,2)=f
                        c_array(num)%bi(3,3)=0
                        c_array(num)%bi(4,1)=0
                        c_array(num)%bi(4,2)=f
                        c_array(num)%bi(4,3)=f
                        c_array(num)%bi(5,1)=f
                        c_array(num)%bi(5,2)=0
                        c_array(num)%bi(5,3)=f
                        c_array(num)%bi(6,1)=f
                        c_array(num)%bi(6,2)=-f
                        c_array(num)%bi(6,3)=0
                        c_array(num)%bi(7,1)=0
                        c_array(num)%bi(7,2)=-f
                        c_array(num)%bi(7,3)=f
                        c_array(num)%bi(8,1)=-f
                        c_array(num)%bi(8,2)=0
                        c_array(num)%bi(8,3)=-f
                        c_array(num)%bi(9,1)=f
                        c_array(num)%bi(9,2)=f
                        c_array(num)%bi(9,3)=0
                        c_array(num)%bi(10,1)=0
                        c_array(num)%bi(10,2)=f
                        c_array(num)%bi(10,3)=f
                        c_array(num)%bi(11,1)=f
                        c_array(num)%bi(11,2)=0
                        c_array(num)%bi(11,3)=-f
                        c_array(num)%bi(12,1)=-f
                        c_array(num)%bi(12,2)=-f
                        c_array(num)%bi(12,3)=0

                        f = 1/sqrt(3.0D0)
                        c_array(num)%ni(1,1)=f
                        c_array(num)%ni(1,2)=f
                        c_array(num)%ni(1,3)=f
                        c_array(num)%ni(2,1)=f
                        c_array(num)%ni(2,2)=f
                        c_array(num)%ni(2,3)=f
                        c_array(num)%ni(3,1)=f
                        c_array(num)%ni(3,2)=f
                        c_array(num)%ni(3,3)=f
                        c_array(num)%ni(4,1)=-f
                        c_array(num)%ni(4,2)=-f
                        c_array(num)%ni(4,3)=f
                        c_array(num)%ni(5,1)=-f
                        c_array(num)%ni(5,2)=-f
                        c_array(num)%ni(5,3)=f
                        c_array(num)%ni(6,1)=-f
                        c_array(num)%ni(6,2)=-f
                        c_array(num)%ni(6,3)=f
                        c_array(num)%ni(7,1)=-f
                        c_array(num)%ni(7,2)=f
                        c_array(num)%ni(7,3)=f
                        c_array(num)%ni(8,1)=-f
                        c_array(num)%ni(8,2)=f
                        c_array(num)%ni(8,3)=f
                        c_array(num)%ni(9,1)=-f
                        c_array(num)%ni(9,2)=f
                        c_array(num)%ni(9,3)=f
                        c_array(num)%ni(10,1)=f
                        c_array(num)%ni(10,2)=-f
                        c_array(num)%ni(10,3)=f
                        c_array(num)%ni(11,1)=f
                        c_array(num)%ni(11,2)=-f
                        c_array(num)%ni(11,3)=f
                        c_array(num)%ni(12,1)=f
                        c_array(num)%ni(12,2)=-f
                        c_array(num)%ni(12,3)=f
                  elseif (c_array(num)%slip_type .eq. 2) then
                        c_array(num)%nslip = 12
                        f = 1/sqrt(2.0D0)
                        
                        c_array(num)%ni(1,1)=0
                        c_array(num)%ni(1,2)=-f
                        c_array(num)%ni(1,3)=f
                        c_array(num)%ni(2,1)=f
                        c_array(num)%ni(2,2)=0
                        c_array(num)%ni(2,3)=-f
                        c_array(num)%ni(3,1)=-f
                        c_array(num)%ni(3,2)=f
                        c_array(num)%ni(3,3)=0
                        c_array(num)%ni(4,1)=0
                        c_array(num)%ni(4,2)=f
                        c_array(num)%ni(4,3)=f
                        c_array(num)%ni(5,1)=f
                        c_array(num)%ni(5,2)=0
                        c_array(num)%ni(5,3)=f
                        c_array(num)%ni(6,1)=f
                        c_array(num)%ni(6,2)=-f
                        c_array(num)%ni(6,3)=0
                        c_array(num)%ni(7,1)=0
                        c_array(num)%ni(7,2)=-f
                        c_array(num)%ni(7,3)=f
                        c_array(num)%ni(8,1)=-f
                        c_array(num)%ni(8,2)=0
                        c_array(num)%ni(8,3)=-f
                        c_array(num)%ni(9,1)=f
                        c_array(num)%ni(9,2)=f
                        c_array(num)%ni(9,3)=0
                        c_array(num)%ni(10,1)=0
                        c_array(num)%ni(10,2)=f
                        c_array(num)%ni(10,3)=f
                        c_array(num)%ni(11,1)=f
                        c_array(num)%ni(11,2)=0
                        c_array(num)%ni(11,3)=-f
                        c_array(num)%ni(12,1)=-f
                        c_array(num)%ni(12,2)=-f
                        c_array(num)%ni(12,3)=0

                        f = 1/sqrt(3.0D0)
                        c_array(num)%bi(1,1)=f
                        c_array(num)%bi(1,2)=f
                        c_array(num)%bi(1,3)=f
                        c_array(num)%bi(2,1)=f
                        c_array(num)%bi(2,2)=f
                        c_array(num)%bi(2,3)=f
                        c_array(num)%bi(3,1)=f
                        c_array(num)%bi(3,2)=f
                        c_array(num)%bi(3,3)=f
                        c_array(num)%bi(4,1)=-f
                        c_array(num)%bi(4,2)=-f
                        c_array(num)%bi(4,3)=f
                        c_array(num)%bi(5,1)=-f
                        c_array(num)%bi(5,2)=-f
                        c_array(num)%bi(5,3)=f
                        c_array(num)%bi(6,1)=-f
                        c_array(num)%bi(6,2)=-f
                        c_array(num)%bi(6,3)=f
                        c_array(num)%bi(7,1)=-f
                        c_array(num)%bi(7,2)=f
                        c_array(num)%bi(7,3)=f
                        c_array(num)%bi(8,1)=-f
                        c_array(num)%bi(8,2)=f
                        c_array(num)%bi(8,3)=f
                        c_array(num)%bi(9,1)=-f
                        c_array(num)%bi(9,2)=f
                        c_array(num)%bi(9,3)=f
                        c_array(num)%bi(10,1)=f
                        c_array(num)%bi(10,2)=-f
                        c_array(num)%bi(10,3)=f
                        c_array(num)%bi(11,1)=f
                        c_array(num)%bi(11,2)=-f
                        c_array(num)%bi(11,3)=f
                        c_array(num)%bi(12,1)=f
                        c_array(num)%bi(12,2)=-f
                        c_array(num)%bi(12,3)=f
                  elseif (c_array(num)%slip_type .eq. 3) then
                        c_array(num)%nslip = 1

                        c_array(num)%bi(1,1) = 1.0
                        c_array(num)%bi(1,2) = 0.0
                        c_array(num)%bi(1,3) = 0.0

                        c_array(num)%ni(1,1) = 0.0
                        c_array(num)%ni(1,2) = 1.0
                        c_array(num)%ni(1,3) = 0.0
                  else
                        write (out,*) "Error: invalid slip type."
                        call die_gracefully
                  end if
c
c                
                  e = c_array(num)%e
                  v = c_array(num)%nu
                  u = c_array(num)%mu
                  if (c_array(num)%elastic_type .eq. 1) then
                        c_array(num)%elast_flex(1,1)=1/e
                        c_array(num)%elast_flex(1,2)=-v/e
                        c_array(num)%elast_flex(1,3)=-v/e
                        c_array(num)%elast_flex(1,4)=0
                        c_array(num)%elast_flex(1,5)=0
                        c_array(num)%elast_flex(1,6)=0
                        c_array(num)%elast_flex(2,1)=-v/e
                        c_array(num)%elast_flex(2,2)=1/e
                        c_array(num)%elast_flex(2,3)=-v/e
                        c_array(num)%elast_flex(2,4)=0
                        c_array(num)%elast_flex(2,5)=0
                        c_array(num)%elast_flex(2,6)=0
                        c_array(num)%elast_flex(3,1)=-v/e
                        c_array(num)%elast_flex(3,2)=-v/e
                        c_array(num)%elast_flex(3,3)=1/e
                        c_array(num)%elast_flex(3,4)=0
                        c_array(num)%elast_flex(3,5)=0
                        c_array(num)%elast_flex(3,6)=0
                        c_array(num)%elast_flex(4,1)=0
                        c_array(num)%elast_flex(4,2)=0
                        c_array(num)%elast_flex(4,3)=0
                        c_array(num)%elast_flex(4,4)=2*(1+v)/e
                        c_array(num)%elast_flex(4,5)=0
                        c_array(num)%elast_flex(4,6)=0
                        c_array(num)%elast_flex(5,1)=0
                        c_array(num)%elast_flex(5,2)=0
                        c_array(num)%elast_flex(5,3)=0
                        c_array(num)%elast_flex(5,4)=0
                        c_array(num)%elast_flex(5,5)=2*(1+v)/e
                        c_array(num)%elast_flex(5,6)=0
                        c_array(num)%elast_flex(6,1)=0
                        c_array(num)%elast_flex(6,2)=0
                        c_array(num)%elast_flex(6,3)=0
                        c_array(num)%elast_flex(6,4)=0
                        c_array(num)%elast_flex(6,5)=0
                        c_array(num)%elast_flex(6,6)=2*(1+v)/e
                  elseif (c_array(num)%elastic_type .eq. 2) then
                        c_array(num)%elast_flex(1,1)=1/e
                        c_array(num)%elast_flex(1,2)=-v/e
                        c_array(num)%elast_flex(1,3)=-v/e
                        c_array(num)%elast_flex(1,4)=0
                        c_array(num)%elast_flex(1,5)=0
                        c_array(num)%elast_flex(1,6)=0
                        c_array(num)%elast_flex(2,1)=-v/e
                        c_array(num)%elast_flex(2,2)=1/e
                        c_array(num)%elast_flex(2,3)=-v/e
                        c_array(num)%elast_flex(2,4)=0
                        c_array(num)%elast_flex(2,5)=0
                        c_array(num)%elast_flex(2,6)=0
                        c_array(num)%elast_flex(3,1)=-v/e
                        c_array(num)%elast_flex(3,2)=-v/e
                        c_array(num)%elast_flex(3,3)=1/e
                        c_array(num)%elast_flex(3,4)=0
                        c_array(num)%elast_flex(3,5)=0
                        c_array(num)%elast_flex(3,6)=0
                        c_array(num)%elast_flex(4,1)=0
                        c_array(num)%elast_flex(4,2)=0
                        c_array(num)%elast_flex(4,3)=0
                        c_array(num)%elast_flex(4,4)=1/u
                        c_array(num)%elast_flex(4,5)=0
                        c_array(num)%elast_flex(4,6)=0
                        c_array(num)%elast_flex(5,1)=0
                        c_array(num)%elast_flex(5,2)=0
                        c_array(num)%elast_flex(5,3)=0
                        c_array(num)%elast_flex(5,4)=0
                        c_array(num)%elast_flex(5,5)=1/u
                        c_array(num)%elast_flex(5,6)=0
                        c_array(num)%elast_flex(6,1)=0
                        c_array(num)%elast_flex(6,2)=0
                        c_array(num)%elast_flex(6,3)=0
                        c_array(num)%elast_flex(6,4)=0
                        c_array(num)%elast_flex(6,5)=0
                        c_array(num)%elast_flex(6,6)=1/u
                  else
                        write (out,*) "Error: invalid elasticity type."
                        call die_gracefully
                  end if
c                       Generate (one time) the inverse of the flexibility
                  c_array(num)%elast_stiff = c_array(num)%elast_flex
                  call mm10_invsym(c_array(num)%elast_stiff,6)
c
                  c_array(num)%valid = .true.
            end subroutine
c
c                 Debug routine, dump the definition to STDOUT
            subroutine print_crystal(num)
                  integer, intent(in) :: num
                  if (num .gt. max_crystals) then
                        write (*,*)
     &                   "Error: crystal number exceeds max crystals."
                        call die_gracefully
                  end if                  
                  if (num .lt. 1) then
                        write (*,*)
     &                   "Error: invalid crystal number."
                        call die_gracefully
                  end if
                  write (*,*) c_array(num)%slip_type,
     &                        c_array(num)%elastic_type
                  write (*,*) ""
                  write (*,*) c_array(num)%e,c_array(num)%nu,
     &                        c_array(num)%mu
                  write (*,*) ""
                  write (*,*) c_array(num)%harden_n,c_array(num)%tau_a,
     &                        c_array(num)%tau_hat_v,c_array(num)%g_o_v,
     &                        c_array(num)%tau_hat_y,c_array(num)%g_o_y,
     &                        c_array(num)%b,c_array(num)%p_v,
     &                        c_array(num)%q_v,c_array(num)%p_y,
     &                        c_array(num)%q_y,c_array(num)%boltz,
     &                        c_array(num)%eps_dot_o_y,
     &                        c_array(num)%eps_dot_o_v,
     &                        c_array(num)%mu_o,c_array(num)%t_o,
     &                        c_array(num)%D_o,
     &                        c_array(num)%theta_o,
     &                        c_array(num)%k_o
                  write (*,*) ""
                  write (*,*) c_array(num)%elast_stiff
                  write (*,*) ""
                  write (*,*) c_array(num)%elast_flex
                  write (*,*) ""
                  write (*,*) c_array(num)%nslip
                  write (*,*) c_array(num)%ni
                  write (*,*) c_array(num)%bi
                  write (*,*)
                  write (*,*) c_array(num)%h_type
            end subroutine
      end module crystal_data
c
c     ****************************************************************
c     *                                                              *
c     *                       subroutine read_simple_angles          *
c     *                                                              *
c     *                       written by: mcm                        *
c     *                       last modified: 3/10/14                 *
c     *                                                              *
c     *     Read angles to file for the damaged interface material   *
c     *                                                              *
c     ****************************************************************
c
      subroutine read_simple_angles
      use main_data
      use crystal_data, only: simple_angles, nangles, srequired,
     &      mc_array
      implicit integer (a-z)
$add common.main
c
      integer :: inter_mat, nlines, avail_device
      integer, external :: warp3d_get_device_number
      character (len=24) :: filen
      double precision :: t1, t2, t3
      real, dimension(max_crystals) :: rand_reals
c
c           Check to see if the simple angles are required
c
      srequired = .false.
      do el=1,noelem
            inter_mat = elstor(11,el)
            if (inter_mat .ne. -1) then
                  srequired = .true.
                  filen = smatprp(112,inter_mat) 
                  exit
            end if
      end do
c
      if (.not. srequired) return
c
      avail_device = warp3d_get_device_number()
      if (avail_device .eq. -1) then
         write(*,*) "No device available."
         call die_gracefully
      end if
c
c     Count the number of lines
      nlines = 0
      open(avail_device, file=filen)
      do
        read(avail_device,*,END=100)
        nlines = nlines + 1
      end do
 100  close(avail_device)
c
c     Read in all the data
      nangles = nlines
      allocate(simple_angles(nangles,3))
c
      i = 1
      open(avail_device, file=filen)
      do
        read(avail_device, *, end=200) t1, t2, t3
        simple_angles(i,1) = t1
        simple_angles(i,2) = t2
        simple_angles(i,3) = t3
        i = i + 1
      end do
 200  close(avail_device)
c
c
c           This is an overly-expensive way of doing this, but oh well
c           Setup mc_array
      allocate(mc_array(noelem, max_crystals))
      call init_random_seed()
      do i=1,noelem
        call random_number(rand_reals)
        mc_array(i,1:max_crystals) = int(rand_reals*nangles)+1
      end do
c
      call wmpi_send_simple_angles
c
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine read_crystal_data            *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 4/17/12 mcm                *
c     *                                                              *
c     *      read in crystal data from file to in memory structures  *
c     *                                                              *
c     ****************************************************************
c
      subroutine read_crystal_data
      use crystal_data, only : angle_input, crystal_input,
     &            data_offset, defined_crystal
      use main_data
      implicit integer (a-z)
$add common.main
c
      integer :: mxcry, nelem, ncry matnum, el, i, iodev, p_matnum
      integer :: avail_device
      integer, external :: warp3d_get_device_number
      logical :: countme, crystalsl, anglesl, open_file, send_mess
      character (len=24) :: filen
c     
c     Just skip if we don't have a crystal
c
      if (.not. defined_crystal) then
            return
      end if
c
c
c     Forward pass: get numbers to allocate and offsets
c
      allocate(data_offset(noelem))
      data_offset = 0
      mxcry = 0
      nelem = 0
      do el=1,noelem
            matnum = elstor(2,el)
            mattype = matprp(9,matnum)
            if (mattype .eq. 10) then
                  countme = (imatprp(104,matnum) .eq. 2) .or.
     &                  (imatprp(107,matnum) .eq. 2)
                  if (countme) then
                        nelem = nelem +1
                        data_offset(el) = nelem
                        ncry = imatprp(101,matnum)
                        if (ncry .gt. mxcry) mxcry = ncry
                  end if
            end if
      end do
c      
c     Backward pass: actual data
c
      send_mess = .false.
      if( any(data_offset .gt. 0) ) then
            write(out,*) '' 
            write(out,'(A)') '>> Reading crystal definitions ...'
            send_mess = .true.
            allocate(angle_input(nelem,mxcry,3))
            allocate(crystal_input(nelem,mxcry))
            p_matnum = -1
            iodev = warp3d_get_device_number()
            if (iodev .eq. -1) then
              write(*,*) "No available device."
              call die_gracefully
            end if
            open_file = .false.
            do el=1,noelem
                  if (data_offset(el) .ne. 0) then
                        matnum = elstor(2,el)
                        if (matnum .ne. p_matnum) then
                          if (el .ne. 1) then
                              close(iodev)
                          end if
                          filen = smatprp(112,matnum)
                          open(iodev,FILE=filen,READONLY)
                          open_file = .true.
                        end if
                        ncry = imatprp(101,matnum)
c                        filen = smatprp(112,matnum)
c                        open(iodev,FILE=filen,READONLY)
                        if (imatprp(104,matnum) .eq. 2) then
                              crystalsl = .true.
                        else
                              crystalsl = .false.
                        end if
                        if (imatprp(107,matnum) .eq. 2) then
                              anglesl = .true.
                        else
                              anglesl = .false.
                        end if
                        call read_defs(iodev,el,ncry,anglesl,
     &                        crystalsl, angle_input(el,1:ncry,1:3),
     &                        crystal_input(el,1:ncry),out)
c                        close(iodev)
                        p_matnum = matnum
                  end if
            end do
            if (open_file) then
                  close(iodev)
            end if
      end if
c
c           If we're using MPI, send out the arrays (including the
c           general crystal struct).  If not it's a dummy routine
c
      call wmpi_send_crystals

      if( send_mess) then
         write(out,*) '    ... Done'
         write(out,*) 
      end if
c     
      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine read_defs                       *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 03/23/2012                 *
c     *                                                              *
c     *     Helper which reads crystal numbers and orientations from *
c     *     a flat file                                              *
c     *                                                              *
c     ****************************************************************
c
      subroutine read_defs(ionum,elnum,ncry,angles,crystals,results_ang,
     &            results_cry,out)
            integer, intent(in) :: ionum, elnum, ncry, out
            logical, intent(in) :: angles, crystals
            double precision, dimension(ncry,3), intent(out) :: 
     &            results_ang
            double precision :: a,b,c
            integer :: n, d, i
            integer, dimension(ncry), intent(out) :: results_cry
c
            character (len=24) :: fmat

c           Scan the file until we find the first entry which matches el
            do
c                 Read into dummy variables just in case
                  if (angles .and. crystals) then
                        read(ionum, *,end=10) n,a,b,c,d
                  elseif (angles) then
                        read(ionum,*,end=10) n,a,b,c
                  elseif (crystals) then
                        read(ionum,*,end=10) n,d
                  else
                        exit
                  endif
c                 Test and input
                  i = 1
                  if (n .eq. elnum) then
                        if (angles .and. crystals) then
                              results_ang(i,1) = a
                              results_ang(i,2) = b
                              results_ang(i,3) = c
                              results_cry(i) = d
                        elseif (angles) then
                              results_ang(i,1) = a
                              results_ang(i,2) = b
                              results_ang(i,3) = c
                        elseif (crystals) then
                              results_cry(i) = d
                        else
                              exit
                        endif
c                             Actual loop to read data
                        do i=2,ncry
                              if (angles .and. crystals) then
                                    read(ionum, *,end=12) n,
     &                                    results_ang(i,1),
     &                                    results_ang(i,2),
     &                                    results_ang(i,3),
     &                                    results_cry(i)
                              elseif (angles) then
                                    read(ionum, *,end=12) n,
     &                                    results_ang(i,1),
     &                                    results_ang(i,2),
     &                                    results_ang(i,3)
                              elseif (crystals) then
                                    read(ionum,*,end=12) n,
     &                                    results_cry(i)
                              else
                                    exit
                              endif
                              if (n .ne. elnum) then
c                                   Terrible error
                                    goto 12
                              end if
                        end do
                        exit
                  end if


            end do


            return
10    continue
            write (out,14) elnum
            call die_gracefully

12    continue
            write (out,13) elnum
            call die_gracefully

13    format(/1x,'>>>>> Parse Error: insufficient data for element ',
     &            i4/)
14    format(/1x,'>>>>> Parse Error: element ', i4, ' not found.' /)

      end subroutine
c
c     ****************************************************************
c     *                                                              *
c     *                   subroutine cleanup_crystal                 *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 04/18/2012                 *
c     *                                                              *
c     *     Call at the end of everything to clean up allocatable    *
c     *     arrays.                                                  *
c     *                                                              *
c     ****************************************************************
c
      subroutine cleanup_crystal
      use crystal_data, only : angle_input, crystal_input,
     &            data_offset, simple_angles, mc_array
      implicit integer (a-z)
$add common.main
c
      if (allocated(angle_input)) deallocate(angle_input)
      if (allocated(crystal_input)) deallocate(crystal_input)
      if (allocated(data_offset)) deallocate(data_offset)
      if (allocated(simple_angles)) deallocate(simple_angles)
      if (allocated(mc_array)) deallocate(mc_array)

      call wmpi_dealloc_crystals

      return

      end subroutine


c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine avg_cry_elast_props          *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 4/5/12 mcm                 *
c     *                                                              *
c     *      Set average elastic properties.  For elastic materials  *
c     *           this will be exact, for anisotropic materials we   *
c     *           need to think about it more                        *
c     *                                                              *
c     ****************************************************************
c
      subroutine avg_cry_elast_props
      use crystal_data, only : angle_input, crystal_input,
     &            data_offset, c_array, defined_crystal
      use main_data
      implicit integer (a-z)
$add common.main
c
      integer :: i, j, k, ncrystals, cnum, osn, num, ecount
      real :: e_avg, nu_avg
c
c
c     Just skip if we don't have a crystal
c
      if (.not. defined_crystal) then
            return
      end if
c
c     Run through our materials, find the CP materials, average their elastic
c     constants
      do i=1,nummat
      if (matprp(9,i) .eq. 10) then
        matprp(1,i) = 0.0
        matprp(2,i) = 0.0
        ecount = 0
        do j = 1, noelem
          if (iprops(38,j)  .eq. i) then
            ecount = ecount + 1
            e_avg = 0.0
            nu_avg = 0.0
            ncrystals = imatprp(101, i)
            do k = 1, ncrystals
c             Get the local crystal number
              if (imatprp(104,i) .eq. 1) then
                cnum = imatprp(105,i)
              elseif (imatprp(104,i) .eq. 2) then
                osn = data_offset(elnum)
                cnum = crystal_input(osn,k)
c               Couldn't do this earlier, so check here
                if ((cnum .gt. max_crystals) .or. 
     &                (cnum .lt. 0)) then
                  write (out,'("Crystal ", i3, " not valid")')
     &                 cnum
                  call die_gracefully
                 end if
              else
                write(out,9502) 
                call die_gracefully
              end if
c              
c             INSERT AVERAGING
c   
              e_avg = e_avg + SNGL(c_array(cnum)%e)
              nu_avg = nu_avg + SNGL(c_array(cnum)%nu)
            end do
            e_avg = e_avg / SNGL(ncrystals)
            nu_avg = nu_avg / SNGL(ncrystals)
            props(7,j) = e_avg
            props(8,j) = nu_avg
            matprp(1,i) = matprp(1,i) + e_avg
            matprp(2,i) = matprp(2,i) + nu_avg
          end if
        end do
        matprp(1,i) = matprp(1,i) / SNGL(ecount)
        matprp(2,i) = matprp(2,i) / SNGL(ecount)
      end if
      end do
c
c     I maintain this is not my fault, we now need to go fix the element props
c
      do i=1,noelem
            if (iprops(25, i) .eq. 10) then
                  num = iprops(38,i)
                  props(7,i) = matprp(1,num)
                  props(8,i) = matprp(2,num)
            end if
      end do

      return

 9502 format(/,1x,
     & '>>>> System error: unexpected input type in avg_elast_props!',
     &       ' Aborting.'/)

      end subroutine
c
c
c
c     ****************************************************************
c     *                                                              *
c     *                      subroutine init_random_see              *
c     *                                                              *
c     *                       written by : mcm                       *
c     *                                                              *
c     *                   last modified : 3/24/14                    *
c     *                                                              *
c     *      Setup the random seed for number generation             *
c     *      This function copied from the GNU documentation         *
c     *                                                              *
c     ****************************************************************
c
      subroutine init_random_seed()
      use ifport, only: getpid
      implicit none
      integer, allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8), pid, t(2), s
      integer(8) :: count, tms
          
      call random_seed(size = n)
      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream", 
     &   form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
         read(un) seed
         close(un)
      else
         ! Fallback to XOR:ing the current time and pid. The PID is
         ! useful in case one launches multiple instances of the same
         ! program in parallel.
         call system_clock(count)
         if (count /= 0) then
            t = transfer(count, t)
         else
            call date_and_time(values=dt)
            tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 
     &           + dt(2) * 31_8 * 24 * 60 * 60 * 1000
     &           + dt(3) * 24 * 60 * 60 * 60 * 1000
     &           + dt(5) * 60 * 60 * 1000 
     &           + dt(6) * 60 * 1000 + dt(7) * 1000
     &           + dt(8)
            t = transfer(tms, t)
         end if
         s = ieor(t(1), t(2))
         pid = getpid() + 1099279 ! Add a prime
         s = ieor(s, pid)
         if (n >= 3) then
            seed(1) = t(1) + 36269
            seed(2) = t(2) + 72551
            seed(3) = pid
            if (n > 3) then
               seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
            end if
         else
            seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
         end if
      end if
      call random_seed(put=seed)
      end subroutine init_random_seed
