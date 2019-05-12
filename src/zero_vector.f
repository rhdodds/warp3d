c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine zero_vector                  *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 07/28/94                   *          
c     *                                                              *          
c     *     zero a vector of specified length w/ floating zero       *          
c     *     signle or double based on this port                      *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine zero_vector( vec, n )                                          
      double precision                                                          
     &  vec(*), zero                                                            
      data zero / 0.0d00 /                                                      
c                                                                               
!DIR$ IVDEP                                                                     
      vec(1:n) = zero                                                           
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine vec_ops                      *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 02/27/96                   *          
c     *                                                              *          
c     *     perform simple vector-vector operations on single or     *          
c     *     double precision vectors (based on the port)             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine vec_ops( veca, vecb, vecc, n, opcode )                         
      double precision                                                          
     &  veca(*), vecb(*), vecc(*), zero, const                                  
      integer opcode                                                            
      data zero / 0.0d0/                                                        
c                                                                               
      go to ( 100, 200, 300, 400, 500, 600, 700 ), opcode                       
c                                                                               
c            opcode 1:   c = a * b                                              
c                                                                               
 100  continue                                                                  
!DIR$ IVDEP                                                                     
      vecc(1:n) = veca(1:n) * vecb(1:n)                                         
      return                                                                    
c                                                                               
c            opcode 2:   b = b * a                                              
c                                                                               
 200  continue                                                                  
!DIR$ IVDEP                                                                     
      vecb(1:n) = vecb(1:n) * veca(1:n)                                         
      return                                                                    
c                                                                               
c            opcode 3:   c = a / b                                              
c                                                                               
 300  continue                                                                  
!DIR$ IVDEP                                                                     
      vecc(1:n) = veca(1:n) / vecb(1:n)                                         
      return                                                                    
c                                                                               
c            opcode v:   c = zero                                               
c                                                                               
 400  continue                                                                  
!DIR$ IVDEP                                                                     
      vecc(1:n) = zero                                                          
      return                                                                    
c                                                                               
c            opcode v:   a = b                                                  
c                                                                               
 500  continue                                                                  
!DIR$ IVDEP                                                                     
      veca(1:n) = vecb(1:n)                                                     
      return                                                                    
c                                                                               
c            opcode v:   a = a + b                                              
c                                                                               
 600  continue                                                                  
!DIR$ IVDEP                                                                     
      veca(1:n) = veca(1:n) + vecb(1:n)                                         
      return                                                                    
c                                                                               
c            opcode v:   a = const * b ; const = vecc(1)                        
c                                                                               
 700  continue                                                                  
      const = vecc(1)                                                           
!DIR$ IVDEP                                                                     
      veca(1:n) = const * vecb(1:n)                                             
      return                                                                    
c                                                                               
      end                                                                       
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine warp3d_sort_float            *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 02/28/00                   *          
c     *                                                              *          
c     *    sort a vector of floating point numbers and an associated *          
c     *    index vector into ascending order                         *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine warp3d_sort_float( n, dvec, index_vec )                        
      double precision                                                          
     &  dvec(*), a                                                              
      integer index_vec(*), b                                                   
c                                                                               
      do j = 2, n                                                               
        a = dvec(j)                                                             
        b = index_vec(j)                                                        
        do i = j-1, 1, -1                                                       
          if ( dvec(i) .le. a ) go to 10                                        
          dvec(i+1) = dvec(i)                                                   
          index_vec(i+1) = index_vec(i)                                         
        end do                                                                  
        i = 0                                                                   
 10     continue                                                                
        dvec(i+1) = a                                                           
        index_vec(i+1) = b                                                      
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       

!                                                                            
! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, 
! Carla Verdi, Feliciano Giustino 
! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, 
! Feliciano Giustino  
!                                                                            
! This file is distributed under the terms of the GNU General Public         
! License. See the file `LICENSE' in the root directory of the               
! present distribution, or http://www.gnu.org/copyleft.gpl.txt . 
!
! http://epw.org.uk/uploads/Ford2/sourcefile/sort.f90.html            
!                                                                            
! Adapted from flib/hpsort_eps
!---------------------------------------------------------------------
!
       subroutine hpsort_eps_epw (n, ra, ind, eps)
!---------------------------------------------------------------------
! sort an array ra(1:n) into ascending order using heapsort algorithm,
! and considering two elements being equal if their values differ
! for less than "eps".
! n is input, ra is replaced on output by its sorted rearrangement.
! create an index table (ind) by making an exchange in the index array
! whenever an exchange is made on the sorted data array (ra).
! in case of equal values in the data array (ra) the values in the
! index array (ind) are used to order the entries.
! if on input ind(1)  = 0 then indices are initialized in the routine,
! if on input ind(1) != 0 then indices are assumed to have been
!                initialized before entering the routine and these
!                indices are carried around during the sorting process
!
! adapted from Numerical Recipes pg. 329 (new edition)
! see web site listed above for code w/ all comments
!
      implicit none  
      integer, intent(in)   :: n  
      double precision, intent(in)  :: eps
      integer :: ind (n)  
      double precision :: ra (n)
      integer :: i, ir, j, l, iind  
      double precision :: rra  
!
      if( ind (1) .eq. 0 ) then  
        do i = 1, n  
          ind (i) = i  
        enddo
      endif
! 
      if( n .lt. 2 ) return  
      l = n / 2 + 1  
      ir = n  
      do ! sorting
       if( l .gt. 1 ) then  
         l    = l - 1  
         rra  = ra (l)  
         iind = ind (l)  
       else  
         rra  = ra (ir)  
         iind = ind (ir)  
         ra (ir) = ra (1)  
         ind (ir) = ind (1)  
         ir = ir - 1  
         if( ir .eq. 1 ) then  
          ra (1)  = rra  
          ind (1) = iind  
          exit  
         endif
       endif
       i = l  
       j = l + l  
       do while( j .le. ir )  
         if( j .lt. ir ) then  
          if( hslt( ra(j), ra(j+1) ) ) then  
             j = j + 1  
          endif
         endif
         if( hslt( rra, ra(j) ) ) then  
          ra (i) = ra (j)  
          ind (i) = ind (j)  
          i = j  
          j = j + j  
         else
          j = ir + 1  
         endif
       enddo  ! on while
       ra (i) = rra  
       ind (i) = iind  
      end do ! sorting  
!  
      contains 
!     ========
!
      logical function hslt( a, b )
      double precision :: a, b
      if( abs(a-b) <  eps ) then
          hslt = .false.
      else
          hslt = ( a >  b )
      end if
      end function hslt
  !
      end subroutine hpsort_eps_epw
                                                                                
                                                                                
                                                                                
                                                                                
