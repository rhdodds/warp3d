c                                                                               
c ********************************************************************          
c *                                                                  *          
c *    write a long vector on unformatted file using                 *          
c *    multiple physical records. last record may not have full      *          
c *    length -- single precision word sizes                         *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c                                                                               
      subroutine wrtbk( fileno, vector, nwords )                                
      implicit integer (a-z)                                                    
      integer vector(*)                                                         
      data max_rec_size                                                         
     &  / 10000 /                                                               
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
          write(fileno) ( vector(ii), ii = laddr, uaddr )                       
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c                                                                               
c ********************************************************************          
c *                                                                  *          
c *    read a long vector on unformatted file using                  *          
c *    multiple physical records. last record may not have full      *          
c *    length - single precision word sizes                          *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c                                                                               
      subroutine rdbk( fileno, vector, nwords )                                 
      implicit integer (a-z)                                                    
      integer vector(*)                                                         
      data max_rec_size                                                         
     &  / 10000 /                                                               
c                                                                               
c        fileno       -- unformatted sequential file no.                        
c        vector       -- data vectors of length nwords to be                    
c                        read (single precision equivalent                      
c                        words. multiply by 2 before calling for                
c                        double precision)                                      
c        max_rec_size -- maximum number of single precision words               
c                        allocated per logical record. generally                
c                        hardware dependent. some machines allow                
c                        any length and sub-divide the logical record           
c                        as required. others have maximum size.                 
c                                                                               
c                                                                               
      nrecs  = (nwords-1) / max_rec_size + 1                                    
      blkfm1 = max_rec_size - 1                                                 
      uaddr  = 0                                                                
      do recno = 1, nrecs                                                       
        laddr = uaddr + 1                                                       
        uaddr = min( nwords, laddr+blkfm1 )                                     
        read(fileno) ( vector(ii), ii = laddr, uaddr )                          
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c ********************************************************************          
c *                                                                  *          
c *    write a long 2d array on unformatted file using               *          
c *    multiple physical records. each column of array is written    *          
c *    as a logical record - single precision word sizes used        *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c                                                                               
      subroutine wrt2d( fileno, array, dimrow, numrows, numcols )               
      implicit integer (a-z)                                                    
      integer array(dimrow,*)                                                   
c                                                                               
      do col = 1, numcols                                                       
          call wrtbk( fileno, array(1,col), numrows )                           
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c ********************************************************************          
c *                                                                  *          
c *    write the used portion of a 3d array on an unformmated        *          
c *    sequential file. - equivalent single precision data is to     *          
c *    be passed in.                                                 *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c                                                                               
      subroutine wrt3d( fileno, array, dimrows, dimcols,                        
     &                  numrows, numcols, numplanes )                           
      implicit integer (a-z)                                                    
      dimension array(dimrows,dimcols,*)                                        
c                                                                               
      do plane = 1, numplanes                                                   
         do col = 1, numcols                                                    
          call wrtbk( fileno, array(1,col,plane), numrows )                     
         end do                                                                 
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c ********************************************************************          
c *                                                                  *          
c *    read the used portion of a 2d array on unformatted file using *          
c *    multiple physical records. each column of array is written    *          
c *    as a logical record - single precision word sizes used        *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c                                                                               
      subroutine rd2d( fileno, array, dimrow, numrows, numcols )                
      implicit integer (a-z)                                                    
      integer array(dimrow,*)                                                   
c                                                                               
      do col = 1, numcols                                                       
          call rdbk( fileno, array(1,col), numrows )                            
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
c                                                                               
c ********************************************************************          
c *                                                                  *          
c *    read the used portion of a 3d array on unformatted file using *          
c *    multiple physical records. each column of array is written    *          
c *    as a logical record - single precision word sizes used        *          
c *                                                                  *          
c ********************************************************************          
c                                                                               
c                                                                               
      subroutine rd3d( fileno, array, dimrows, dimcols,                         
     &                 numrows, numcols, numplanes )                            
      implicit integer (a-z)                                                    
      dimension array(dimrows,dimcols,*)                                        
c                                                                               
      do plane = 1, numplanes                                                   
         do col = 1, numcols                                                    
          call rdbk( fileno, array(1,col,plane), numrows )                      
         end do                                                                 
      end do                                                                    
c                                                                               
      return                                                                    
      end                                                                       
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
                                                                                
