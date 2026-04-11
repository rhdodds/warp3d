c ********************************************************************          
c *                                                                  *          
c *    All updated 4/8/26 rhd                                        *          
c *                                                                  *          
c ********************************************************************          
c ********************************************************************          
c *                                                                  *          
c *    write a long 2d array on unformatted file using               *          
c *    multiple physical records. each column of array is written    *          
c *    as a logical record                                           *          
c *                                                                  *          
c ********************************************************************          
c                                                                                                                                                             
      subroutine wrt2d_i4( fileno, array, dimrow, numrows, numcols )  
c                   
      implicit none
      integer ::  fileno, dimrow, numrows,
     &            numcols, array(dimrow,numcols)                               
c
      integer :: col
c                                                                                     
      do col = 1, numcols                                                       
          call wrtbk_i4( fileno, array(1,col), numrows )                           
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c ********************************************************************          
c *                                                                  *          
c *    write a long 2d array on unformatted file using               *          
c *    multiple physical records. each column of array is written    *          
c *    as a logical record                                           *
c *                                                                  *          
c ********************************************************************          
c                                                                                                                                                           
      subroutine wrt2d_r4( fileno, array, dimrow, numrows, numcols )  
c                   
      implicit none
      integer :: fileno, dimrow, numrows, numcols                                            
      real :: array(dimrow,numcols)
c
      integer :: col
c                                                                                     
      do col = 1, numcols                                                       
          call wrtbk_r4( fileno, array(1,col), numrows )                           
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                               
c ********************************************************************          
c *                                                                  *          
c *    write a long 2d array on unformatted file using               *          
c *    multiple physical records. each column of array is written    *          
c *    as a logical record                                           *
c *                                                                  *          
c ********************************************************************          
c                                                                                                                                                           
      subroutine wrt2d_r8( fileno, array, dimrow, numrows, numcols )  
c                   
      implicit none
      integer :: fileno, dimrow, numrows, numcols                                            
      double precision :: array(dimrow,numcols)
c
      integer :: col
c                                                                                     
      do col = 1, numcols                                                       
          call wrtbk_r8( fileno, array(1,col), numrows )                           
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c ********************************************************************          
c *                                                                  *          
c *    write the used portion of a 3d array on an unformmated        *          
c *    sequential file.                                              *          
c *                                                                  *          
c ********************************************************************          
c                                                                                                                                                           
      subroutine wrt3d_r8( fileno, array, dimrows, dimcols,                        
     &                     numrows, numcols, numplanes )                           
      implicit none
      integer :: fileno, dimrows, dimcols, numrows, numcols, numplanes
c     
      double precision :: array(dimrows,dimcols,numplanes)      
c
      integer :: plane, col                                        
c                                                                               
      do plane = 1, numplanes                                                   
         do col = 1, numcols                                                    
          call wrtbk_r8( fileno, array(1,col,plane), numrows )                     
         end do                                                                 
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
                                                             
                                                                                
c ********************************************************************          
c *                                                                  *          
c *    write a long vector on unformatted file using                 *          
c *    multiple physical records. last record may not have full      *          
c *    length                                                        *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
      subroutine wrtbk_i4( fileno, vector, nvalues )                                
      implicit none
c                                                          
      integer :: fileno, nvalues, vector(nvalues)
c   
      integer :: nrecs, blkfm1, uaddr, recno, laddr                                                      
      integer, parameter :: max_rec_size = 10000                                                               
c                                                                               
c        fileno       -- unformatted sequential file no.                        
c        vector       -- data vectors of length nvalues to be                    
c                        written                              
c        max_rec_size -- maximum number of single precision words               
c                        allocated per logical record. generally                
c                        hardware dependent. some machines allow                
c                        any length and sub-divide the logical record           
c                        as required. others have maximum size.                 
c                                                                               
      nrecs  = (nvalues-1) / max_rec_size + 1                                    
      blkfm1 = max_rec_size - 1                                                 
      uaddr  = 0                                                                
      do recno = 1, nrecs                                                       
          laddr = uaddr + 1                                                     
          uaddr = min( nvalues, laddr+blkfm1 )                                   
          write(fileno) vector(laddr:uaddr)                       
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                          
                                                                                  
c ********************************************************************          
c *                                                                  *          
c *    write a long vector on unformatted file using                 *          
c *    multiple physical records. last record may not have full      *          
c *    length                                                        *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
      subroutine wrtbk_l4( fileno, vector, nvalues )                                
      implicit none
c                                                          
      integer :: fileno, nvalues
      logical :: vector(nvalues)
c   
      integer :: nrecs, blkfm1, uaddr, recno, laddr                                                      
      integer, parameter :: max_rec_size = 10000                                                               
c                                                                               
c        fileno       -- unformatted sequential file no.                        
c        vector       -- data vectors of length nvalues to be                    
c                        written                              
c        max_rec_size -- maximum number of single precision words               
c                        allocated per logical record. generally                
c                        hardware dependent. some machines allow                
c                        any length and sub-divide the logical record           
c                        as required. others have maximum size.                 
c                                                                               
      nrecs  = (nvalues-1) / max_rec_size + 1                                    
      blkfm1 = max_rec_size - 1                                                 
      uaddr  = 0                                                                
      do recno = 1, nrecs                                                       
          laddr = uaddr + 1                                                     
          uaddr = min( nvalues, laddr+blkfm1 )                                   
          write(fileno) vector(laddr:uaddr)                       
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       

                                                                                  
c ********************************************************************          
c *                                                                  *          
c *    write a long vector on unformatted file using                 *          
c *    multiple physical records. last record may not have full      *          
c *    length                                                        *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
      subroutine wrtbk_r4( fileno, vector, nvalues )                                
      implicit none
c                                                          
      integer :: fileno, nvalues
      real :: vector(nvalues)
c   
      integer :: nrecs, blkfm1, uaddr, recno, laddr                                                      
      integer, parameter :: max_rec_size = 10000                                                               
c                                                                               
c        fileno       -- unformatted sequential file no.                        
c        vector       -- data vectors of length nvalues to be                    
c                        written                              
c        max_rec_size -- maximum number of single precision words               
c                        allocated per logical record. generally                
c                        hardware dependent. some machines allow                
c                        any length and sub-divide the logical record           
c                        as required. others have maximum size.                 
c                                                                               
      nrecs  = (nvalues-1) / max_rec_size + 1                                    
      blkfm1 = max_rec_size - 1                                                 
      uaddr  = 0                                                                
      do recno = 1, nrecs                                                       
          laddr = uaddr + 1                                                     
          uaddr = min( nvalues, laddr+blkfm1 )                                   
          write(fileno) vector(laddr:uaddr)                       
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c                                                                               
c ********************************************************************          
c *                                                                  *          
c *    write a long vector on unformatted file using                 *          
c *    multiple physical records. last record may not have full      *          
c *    length                                                        *          
c *                                                                  *          
c ********************************************************************          
c                                                                                                                                                           
      subroutine wrtbk_r8( fileno, vector, nwords )                                
      implicit none
c                                                          
      integer :: fileno, nwords  
      double precision :: vector(*)
c   
      integer :: nrecs, blkfm1, uaddr, recno, laddr                                                      
      integer, parameter :: max_rec_size = 5000
c                                                                               
c        fileno       -- unformatted sequential file no.                        
c        vector       -- data vectors of length nwords to be                    
c                        written (single precision equivalent                   
c                        words. multiply by 2 before calling for                
c                        double precision)                                      
c        max_rec_size -- maximum number of single precision words               
c                        allocated per logical record. generally                
c                        hardware dependent. some machines allow                
c                        any length and sub-divide the logical record           
c                        as required. others have maximum size.                 
c                                                                               
      nrecs  = (nwords-1) / max_rec_size + 1                                    
      blkfm1 = max_rec_size - 1                                                 
      uaddr  = 0                                                                
      do recno = 1, nrecs                                                       
          laddr = uaddr + 1                                                     
          uaddr = min( nwords, laddr+blkfm1 )                                   
          write(fileno) vector(laddr:uaddr)                       
      end do                                                                    
c                                                                               
      return     
      end                                                               
                                                                            
c======================================================================
c
c
c                                                                                
c ********************************************************************          
c *                                                                  *          
c *    read a long vector on unformatted file using                 *          
c *    multiple physical records. last record may not have full      *          
c *    length                                                        *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
      subroutine rdbk_i4( fileno, vector, nvalues )                                
      implicit none
c                                                          
      integer :: fileno, nvalues, vector(nvalues)
c   
      integer :: nrecs, blkfm1, uaddr, recno, laddr                                                      
      integer, parameter :: max_rec_size = 10000                                                               
c                                                                               
c        fileno       -- unformatted sequential file no.                        
c        vector       -- data vectors of length nvalues to be                    
c                        written                              
c        max_rec_size -- maximum number of single precision words               
c                        allocated per logical record. generally                
c                        hardware dependent. some machines allow                
c                        any length and sub-divide the logical record           
c                        as required. others have maximum size.                 
c                                                                               
      nrecs  = (nvalues-1) / max_rec_size + 1                                    
      blkfm1 = max_rec_size - 1                                                 
      uaddr  = 0                                                                
      do recno = 1, nrecs                                                       
          laddr = uaddr + 1                                                     
          uaddr = min( nvalues, laddr+blkfm1 )                                   
          read(fileno) vector(laddr:uaddr)                       
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                          
                                                                                  
c ********************************************************************          
c *                                                                  *          
c *    read a long vector on unformatted file using                 *          
c *    multiple physical records. last record may not have full      *          
c *    length                                                        *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
      subroutine rdbk_l4( fileno, vector, nvalues )                                
      implicit none
c                                                          
      integer :: fileno, nvalues
      logical :: vector(nvalues)
c   
      integer :: nrecs, blkfm1, uaddr, recno, laddr                                                      
      integer, parameter :: max_rec_size = 10000                                                               
c                                                                               
c        fileno       -- unformatted sequential file no.                        
c        vector       -- data vectors of length nvalues to be                    
c                        written                              
c        max_rec_size -- maximum number of single precision words               
c                        allocated per logical record. generally                
c                        hardware dependent. some machines allow                
c                        any length and sub-divide the logical record           
c                        as required. others have maximum size.                 
c                                                                               
      nrecs  = (nvalues-1) / max_rec_size + 1                                    
      blkfm1 = max_rec_size - 1                                                 
      uaddr  = 0                                                                
      do recno = 1, nrecs                                                       
          laddr = uaddr + 1                                                     
          uaddr = min( nvalues, laddr+blkfm1 )                                   
          read(fileno) vector(laddr:uaddr)                       
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       

                                                                                  
c ********************************************************************          
c *                                                                  *          
c *    read a long vector on unformatted file using                  *          
c *    multiple physical records. last record may not have full      *          
c *    length                                                        *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
      subroutine rdbk_r4( fileno, vector, nvalues )                                
      implicit none
c                                                          
      integer :: fileno, nvalues
      real :: vector(nvalues)
c   
      integer :: nrecs, blkfm1, uaddr, recno, laddr                                                      
      integer, parameter :: max_rec_size = 10000                                                               
c                                                                               
c        fileno       -- unformatted sequential file no.                        
c        vector       -- data vectors of length nvalues to be                    
c                        written                              
c        max_rec_size -- maximum number of single precision words               
c                        allocated per logical record. generally                
c                        hardware dependent. some machines allow                
c                        any length and sub-divide the logical record           
c                        as required. others have maximum size.                 
c                                                                               
      nrecs  = (nvalues-1) / max_rec_size + 1                                    
      blkfm1 = max_rec_size - 1                                                 
      uaddr  = 0                                                                
      do recno = 1, nrecs                                                       
          laddr = uaddr + 1                                                     
          uaddr = min( nvalues, laddr+blkfm1 )                                   
          read(fileno) vector(laddr:uaddr)                       
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
c                                                                               
c ********************************************************************          
c *                                                                  *          
c *    read a long vector on unformatted file using                  *          
c *    multiple physical records. last record may not have full      *          
c *    length                                                        *          
c *                                                                  *          
c ********************************************************************          
c                                                                                                                                                           
      subroutine rdbk_r8( fileno, vector, nwords )                                
      implicit none
c                                                          
      integer :: fileno, nwords  
      double precision :: vector(*)
c   
      integer :: nrecs, blkfm1, uaddr, recno, laddr                                                      
      integer, parameter :: max_rec_size = 5000
c                                                                               
c        fileno       -- unformatted sequential file no.                        
c        vector       -- data vectors of length nwords to be                    
c                        written (single precision equivalent                   
c                        words. multiply by 2 before calling for                
c                        double precision)                                      
c        max_rec_size -- maximum number of single precision words               
c                        allocated per logical record. generally                
c                        hardware dependent. some machines allow                
c                        any length and sub-divide the logical record           
c                        as required. others have maximum size.                 
c                                                                               
      nrecs  = (nwords-1) / max_rec_size + 1                                    
      blkfm1 = max_rec_size - 1                                                 
      uaddr  = 0                                                                
      do recno = 1, nrecs                                                       
          laddr = uaddr + 1                                                     
          uaddr = min( nwords, laddr+blkfm1 )                                   
          read(fileno) vector(laddr:uaddr)                       
      end do                                                                    
c                                                                               
      return     
      end                                                               
                                                                            
c ********************************************************************          
c *                                                                  *          
c *    read along 2d array on unformatted file using                 *          
c *    multiple physical records. each column of array is written    *          
c *    as a logical record                                           *          
c *                                                                  *          
c ********************************************************************          
c                                                                                                                                                             
      subroutine rd2d_i4( fileno, array, dimrow, numrows, numcols )  
c                   
      implicit none
      integer :: fileno, dimrow, numrows,
     &           numcols,  array(dimrow,numcols)                                           
c
      integer :: col
c                                                                                     
      do col = 1, numcols                                                       
          call rdbk_i4( fileno, array(1,col), numrows )                           
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c ********************************************************************          
c *                                                                  *          
c *    write a long 2d array on unformatted file using               *          
c *    multiple physical records. each column of array is written    *          
c *    as a logical record                                           *
c *                                                                  *          
c ********************************************************************          
c                                                                                                                                                           
      subroutine rd2d_r4( fileno, array, dimrow, numrows, numcols )  
c                   
      implicit none
      integer :: fileno, dimrow, numrows, numcols                                            
      real :: array(dimrow,numcols)
c
      integer :: col
c                                                                                     
      do col = 1, numcols                                                       
          call rdbk_r4( fileno, array(1,col), numrows )                           
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                               
c ********************************************************************          
c *                                                                  *          
c *    write a long 2d array on unformatted file using               *          
c *    multiple physical records. each column of array is written    *          
c *    as a logical record                                           *
c *                                                                  *          
c ********************************************************************          
c                                                                                                                                                           
      subroutine rd2d_r8( fileno, array, dimrow, numrows, numcols )  
c                   
      implicit none
      integer :: fileno, dimrow, numrows, numcols                                            
      double precision :: array(dimrow,numcols)
c
      integer :: col
c                                                                                     
      do col = 1, numcols                                                       
          call rdbk_r8( fileno, array(1,col), numrows )                           
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c ********************************************************************          
c *                                                                  *          
c *    write the used portion of a 3d array on an unformmated        *          
c *    sequential file.                                              *          
c *                                                                  *          
c ********************************************************************          
c                                                                                                                                                           
      subroutine rd3d_r8( fileno, array, dimrows, dimcols,                        
     &                     numrows, numcols, numplanes )                           
      implicit none
      integer :: fileno, dimrows, dimcols, numrows, numcols, numplanes
c     
      double precision :: array(dimrows,dimcols,numplanes)      
c
      integer :: plane, col                                        
c                                                                               
      do plane = 1, numplanes                                                   
         do col = 1, numcols                                                    
          call rdbk_r8( fileno, array(1,col,plane), numrows )                     
         end do                                                                 
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                                      
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
