

c     ****************************************************************
c     *                                                              *
c     *                 subroutine mm_return_values                  *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 4/21/23 rhd                *
c     *                                                              *
c     *  return requested results for an element if available        *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine mm_return_values( value_type, elem, vec, npts )
c 
      use global_data, only : out, iprops, nstr, nstrs, elblks
      use main_data,       only : elems_to_blocks                               
      use elem_block_data, only : history_blocks, eps_n_blocks,                 
     &                            urcs_n_blocks, history_blk_list    
      use constants            
c
      implicit none       
c
c              parameters
c
      character(len=*), intent(in) :: value_type
      integer, intent(in) :: elem, npts
      double precision, intent(out) :: vec(npts)
c
c              local variables
c
      integer :: mat_model, blk, rel_elem, hist_size, hoffset,
     &           epsoffset, sigoffset, ngp, span
      double precision, dimension(:), pointer :: history, urcs_n,
     &          eps_n, urcs_n1                       
c
      mat_model = iprops(25,elem) 
c
      ngp         = iprops(6,elem)                                           
      blk         = elems_to_blocks(elem,1)                                  
      rel_elem    = elems_to_blocks(elem,2)  
      span        = elblks(0,blk)
      hist_size   = history_blk_list(blk)                                    
      hoffset     = (rel_elem-1)*hist_size*ngp + 1                           
      epsoffset   = (rel_elem-1)*nstr*ngp                                    
      sigoffset   = (rel_elem-1)*nstrs*ngp                                   
      urcs_n      => urcs_n_blocks(blk)%ptr                              
      eps_n       => eps_n_blocks(blk)%ptr                                 
      history     => history_blocks(blk)%ptr 
c      
      select case( value_type )  ! what type data is available
c
       case( "plastic_strain" )
         call mm_return_values_eps_plas
       case( "avg_mises" )
         call mm_return_values_avg_mises
       case( "avg_princ_stress" ) ! in decreasing order
         call mm_return_values_avg_princ_stress
       case( "avg_porosity" ) ! in decreasing order
         call mm_return_values_avg_porosity
       case default
         write(out,9000) elem, value_type, 1
         call die_gracefully
c
      end select
c
      return
c
 9000 format(/,'FATAL ERROR: mm_return_values. Contact WARP3D group',             
     &       /,'             element, type: ', i8,2x,a,
     &       /,'             Job terminated at ',i1,//)                         
c
      contains
c     ========
      subroutine mm_return_values_avg_princ_stress
      implicit none
c 
      integer :: gp, j, k
      double precision :: sig(6), evec(3,3), evals(3)
c
      if( npts /= ngp ) then
        write(out,9020) elem, npts, ngp, 3
        call die_gracefully
      end if
c
c              average of integration point stresses
c
      do gp = 1, ngp    
        do j = 1, 6                                                     
         sig(j) = sig(j) + urcs_n(sigoffset+j)  
        end do                                      
        sigoffset = sigoffset + nstrs                                        
      end do 
c
c              principal stresses in decreasing order
c              gfortran has troubles with principal_values code
c              evcmp1_new returns principal values in
c              increasing order.
c
      sig = sig / dble( ngp )
      call evcmp1_new( 1, 1, sig, evals )
      if( npts <= 2 ) then ! see size of vec on entry
        write(out,9020) elem, npts, ngp, 4
        call die_gracefully 
      end if
      vec(1) = evals(3)
      vec(2) = evals(2)
      vec(3) = evals(1)
c      call principal_values( sig, evec, evals )
      return 
c
 9020 format(/,'FATAL ERROR: mm_return_values. Contact WARP3D group',             
     &       /,'             mismatched number of points.',
     &       /,'             element, npts, ngp: ',i8,2i3,
     &       /,'             Job terminated at ',i1,//)                         
c
      end subroutine mm_return_values_avg_princ_stress
c
      subroutine mm_return_values_avg_mises
      implicit none
c 
      integer :: gp, j
      double precision :: sig(6), sig_mean, s11, s22, s33, s12,
     &                    s13, s23, j2
c
      if( npts /= ngp ) then
        write(out,9020) elem, npts, ngp, 3
        call die_gracefully
      end if
c
      do gp = 1, ngp    
        do j = 1, 6                                                     
         sig(j) = sig(j) + urcs_n(sigoffset+j)  
        end do                                      
        sigoffset = sigoffset + nstrs                                        
      end do 
c
      sig = sig / dble( ngp )
      sig_mean  = (sig(1) + sig(2) + sig(3)) * third                    
      s11 = sig(1) - sig_mean
      s22 = sig(2) - sig_mean
      s33 = sig(3) - sig_mean
      s12 = sig(4)
      s13 = sig(5)
      s23 = sig(6)
      j2 = half * ( s11**2 + s22**2 + s33**2 + 
     &              two*(s12**2 + s23**2 + s13**2))
      vec(1) = sqrt( three * j2 )
c
      return 
c
 9020 format(/,'FATAL ERROR: mm_return_values. Contact WARP3D group',             
     &       /,'             mismatched number of points.',
     &       /,'             element, npts, ngp: ',i8,2i3,
     &       /,'             Job terminated at ',i1,//)                         
c
      end subroutine mm_return_values_avg_mises
c
      subroutine mm_return_values_eps_plas ! uniaxial plastic strain
      implicit none
c 
      integer :: eps_plas_loc, gp, j, eps_plas_col
      double precision :: convert_factor
c
      if( npts /= ngp ) then
        write(out,9020) elem, npts, ngp, 3
        call die_gracefully
      end if
c
      select case( mat_model )
        case( 1 ) ! bilinear   verified
           eps_plas_loc = 3
           do gp = 1, ngp
              j = (gp-1)*hist_size
              vec(gp) = history(hoffset+(eps_plas_loc-1)+j)
           end do 
        case( 2 ) ! deformation verified
          eps_plas_col = 9
          do gp = 1, ngp    
            sigoffset = (rel_elem-1)*nstrs*ngp + (gp-1)*nstrs                     
            vec(gp) = urcs_n(sigoffset + eps_plas_col)                                          
          end do
        case( 3 ) ! mises/gurson  ! verified
           eps_plas_loc = 1
           do gp = 1, ngp
              j = (gp-1)*hist_size
              vec(gp) = history(hoffset+(eps_plas_loc-1)+j)
           end do 
        case( 5 ) ! cyclic   verified
           eps_plas_loc = 3
           convert_factor = root23
           do gp = 1, ngp
             j = (gp-1)*hist_size
             vec(gp) = history(hoffset+(eps_plas_loc-1)+j) * 
     &                 convert_factor 
           end do 
        case( 6 ) ! norton creep    verified 
           eps_plas_loc = 1    ! accumulated creep strain
           do gp = 1, ngp
             j = (gp-1)*hist_size
             vec(gp) = history(hoffset+(eps_plas_loc-1)+j) 
           end do 
        case( 7 ) ! mises-hydrogen  ??? run uniaxial to find factor
           eps_plas_loc = 1
           do gp = 1, ngp
              j = (gp-1)*hist_size
              vec(gp) = history(hoffset+(eps_plas_loc-1)+j) 
           end do 
        case( 10 ) ! crystal plasticity   verified
           eps_plas_col = 9
           do gp = 1, ngp    
            sigoffset = (rel_elem-1)*nstrs*ngp + (gp-1)*nstrs                     
            vec(gp) = urcs_n(sigoffset + eps_plas_col)                                          
           end do
        case default
           write(out,9010) elem, mat_model, 2
           call die_gracefully
      end select
c
      return 
c
 9010 format(/,'FATAL ERROR: mm_return_values. Contact WARP3D group',             
     &       /,'             cannot return plastic strains.',
     &       /,'             element, material model: ',i8,i3,
     &       /,'             Job terminated at ',i1,//)                         
 9020 format(/,'FATAL ERROR: mm_return_values. Contact WARP3D group',             
     &       /,'             mismatched number of points.',
     &       /,'             element, npts, ngp: ',i8,2i3,
     &       /,'             Job terminated at ',i1,//)                         
c
      end subroutine mm_return_values_eps_plas
      
      subroutine mm_return_values_avg_porosity
      implicit none
c 
      integer :: gp, offset
      double precision :: porosity
c
      porosity = zero
      do gp = 1, ngp    
        offset = hoffset+4+(gp-1)*hist_size                                 
        porosity = porosity + history(offset)     
      end do                                                      
      porosity = porosity / dble( ngp )
      vec(1) = porosity
c
      return 
c
 9020 format(/,'FATAL ERROR: mm_return_values. Contact WARP3D group',             
     &       /,'             mismatched number of points.',
     &       /,'             element, npts, ngp: ',i8,2i3,
     &       /,'             Job terminated at ',i1,//)                         
c
      end subroutine mm_return_values_avg_porosity

c
      end subroutine mm_return_values
