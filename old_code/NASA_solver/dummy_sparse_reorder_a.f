c
c     ****************************************************************
c     *                                                              *
c     *           dummy routines for sparse solver                   *
c     ****************************************************************
c
       subroutine sparse_reorder( neqns, ncoeff, k_diag, p_vec,
     &                           k_coeffs,
     &                           k_ptrs, k_indexes, perm, xadj, adjncy,
     &                           reorder_type )
      implicit integer (a-z)
      write(*,*) ' '
      write(*,*) '>> the sparse solver is not available in this'
      write(*,*) '>> version of warp3d. please use the direct'
      write(*,*) '>> solver or the pcg solver....'
      write(*,*) ' '
      write(*,*) '>> please read the file nasa_solver_notes'
      write(*,*) '>> in the main WARP3D directory for information'
      write(*,*) '>> about installing the sparse solver...'      
      write(*,*) ' '
      write(*,*) '>> job terminated...'
      write(*,*) ' '
      call die_gracefully
      stop
      end
