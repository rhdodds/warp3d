c     ****************************************************************
c     *                                                              *
c     *                      subroutine oupstr_elem                  *
c     *                                                              *
c     *                       written by : kck                       *
c     *                                                              *
c     *                   last modified : 4/11/2014 rhd              *
c     *                                                              *
c     *  drive output of stress or strain element results to         *
c     *  (1) a Patran file in either binary or formatted forms       *
c     *  or (2) a flat file with text or stream format               *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine oupstr_elem( stress, oubin, ouasc, flat_file,
     &                      stream_file, text_file, compressed )
      use main_data, only : cohesive_ele_types
      implicit integer (a-z)
      logical :: stress, oubin, ouasc, flat_file,
     &           stream_file, text_file, compressed  
      include 'common.main'
c
c                       locally allocated. the big one is
c                       elem_results.
c
      logical ::  bbar_flg, geo_non_flg, long_out_flg, nodpts_flg,
     &            center_output, first_block, last_block, cohesive_elem
c
c
      double precision ::
     &     zero, elem_results, small_tol
      allocatable elem_results(:,:)
      integer :: elem_out_map(mxelmp)
      data zero, small_tol / 0.0d0, 1.0d-80 /
c
c
      allocate( elem_results(mxvl,mxstmp) )
      first_block = .true.
      last_block  = .false.
c
c                       compute the stress/strain element output
c                       vectors and assemble them into the global
c                       element results for each block of similar, 
c                       non-conflicting elements.
c
      do blk = 1, nelblk
         if ( elblks(2,blk) .ne. myid ) cycle
         elem_results(1:mxvl,1:mxstmp) = zero
         span              = elblks(0,blk)
         felem             = elblks(1,blk)
         elem_type         = iprops(1,felem)
         int_order         = iprops(5,felem)
         mat_type          = iprops(25,felem)
         num_enodes        = iprops(2,felem)
         num_enode_dof     = iprops(4,felem)
         totdof            = num_enodes * num_enode_dof
         geo_non_flg       = lprops(18,felem)
         bbar_flg          = lprops(19,felem)
         num_int_points    = iprops(6,felem)
         long_out_flg      = lprops(16,felem)
         output_loc        = iprops(12,felem)
         nodpts_flg        = .true.
         center_output     = .true.
         num_short_stress  = 11
         num_short_strain  = 7
         num_vals          = num_short_strain + 15
         if ( stress ) num_vals = num_short_stress + 15
         cohesive_elem     = cohesive_ele_types(elem_type)
c
c                       for cohesive elements, output zeroes for 
c                       element result values.
c
         if ( cohesive_elem ) then
           call oust_elem( stress, oubin, ouasc, num_vals, 
     &                     elem_results, mxvl, span, first_block,
     &                     felem, last_block, flat_file,
     &                     stream_file, text_file, compressed )
           first_block = .false.
           cycle
         end if
c
c                       process all solid element types
c           
         do i = 1, mxelmp
           elem_out_map(i) = outmap(elem_type,i)
         end do
c
c                       duplicate necessary element block data.
c
         call oudups( span, felem, num_int_points, geo_non_flg, stress,
     &                  .false. )
c
c                       compute the element block center stress/strain
c                       data.
         ibkl = blk
         ifelem = felem
         call ouprks( span, iblk, ifelem, elem_type, int_order,
     &                num_int_points,
     &                num_enodes, geo_non_flg, stress, mat_type,
     &                center_output, num_short_stress, num_short_strain,
     &                .true. )
c
c                       add the element block element stress/strain data
c                       into the global element results array.
c
         call oupele( span, num_short_strain, num_short_stress,
     &                stress, elem_results(1,1), mxvl )
c
c                       calculate the extended values of stress and strain
c                       for each element
c
         call ouext2 ( elem_results, mxvl, span, stress )
c
c                       output the element results
c                       to a file compatable with patran for post
c                       processing.
c
c                       zero small values to prevent 3-digit exponents
c                       in formated output
c
         where( abs(elem_results) .lt. small_tol ) elem_results = zero
         call oust_elem( stress, oubin, ouasc, num_vals, elem_results,
     &                   mxvl, span, first_block, felem, last_block,
     &                   flat_file,stream_file, text_file, 
     &                   compressed )
         first_block = .false.
c
      end do
c
c
      where( abs(elem_results) .lt. small_tol ) elem_results = zero
      last_block = .true.
      call oust_elem( stress, oubin, ouasc, num_vals, elem_results,
     &                mxvl, span, first_block, felem, last_block,
     &                flat_file, stream_file, text_file, compressed )
c
      deallocate( elem_results )
c
      return
      end






