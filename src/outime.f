c     ****************************************************************          
c     *                                                              *          
c     *                subroutine outime                             *          
c     *  output cpu times for each major part of solution at         *          
c     *  job termination                                             *          
c     *                                                              *          
c     *                                                              *          
c     *                       written by : rhd                       *          
c     *                                                              *          
c     *                   last modified : 01/20/2017 rhd             *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
      subroutine outime                                                         
      use global_data ! old common.main
      implicit none                                                             
c                                                                               
      integer :: calc                                                           
      character(len=30) :: clctyp                                               
      real :: t1                                                                
      real, external :: wcputime                                                
c                                                                               
      write(out,9000)                                                           
c                                                                               
      t1 = wcputime( 1 )                                                        
c                                                                               
      do calc = 1, mxtim                                                        
         if( times(calc,2) .lt. 0.01 ) cycle                                    
         if(calc.eq.1) then                                                     
            clctyp= 'internal force vector:        '                            
         else if(calc.eq.2) then                                                
            clctyp= 'tangent stiffness:            '                            
         else if(calc.eq.3) then                                                
            clctyp= 'linear stiffness:             '                            
         else if(calc.eq.4) then                                                
            clctyp= 'tangent mass:                 '                            
         else if(calc.eq.5) then                                                
            clctyp= 'linear mass:                  '                            
         else if(calc.eq.6) then                                                
            clctyp= 'sig-eps & internal force:     '                            
         else if(calc.eq.7) then                                                
            clctyp= 'patran output:                '                            
         else if(calc.eq.9) then                                                
            clctyp= 'step length:                  '                            
         else if(calc.eq.10) then                                               
            clctyp= 'solution matrix factorization:'                            
         else if(calc.eq.11) then                                               
            clctyp= 'system solution:              '                            
         else if(calc.eq.12) then                                               
            clctyp= 'pseudo residual:              '                            
         else if(calc.eq.13) then                                               
            clctyp= 'beta scaling factor:          '                            
         else if(calc.eq.14) then                                               
            clctyp= 'linear residual update:       '                            
         else if(calc.eq.15) then                                               
            clctyp= 'search direction update:      '                            
         else if(calc.eq.16) then                                               
            clctyp= 'pcg solution vector update:   '                            
         else if (calc.eq.17) then                                              
            clctyp= 'elprd2:                       '                            
         else if (calc.eq.18) then                                              
            clctyp= 'assembly of struct. K:        '                            
         else if (calc.eq.19) then                                              
            clctyp= 'triangularize K:              '                            
         else if (calc.eq.20) then                                              
            clctyp= 'backsubstitue:                '                            
         else if (calc.eq.21) then                                              
            clctyp= 'set-up equations:             '                            
         else if (calc.eq.22) then                                              
            clctyp= 'sparse indexes set-up:        '                            
         else if (calc.eq.23) then                                              
            clctyp= 'reordering, symbolic factor:  '                            
         else if (calc.eq.24) then                                              
            cycle                                                               
         else if (calc.eq.25) then                                              
            clctyp= 'sparse fill, factorization:   '                            
         else if (calc.eq.26) then                                              
            clctyp= 'sparse loadpass:              '                            
         end if                                                                 
         write(out,9010) clctyp                                                 
         write(out,9020) times(calc,1),                                         
     &        100.0*times(calc,1)/t1, int(times(calc,2))                        
c                                                                               
      end do                                                                    
c                                                                               
 9000 format(////1x,'>>>>>  solution timings   <<<<<')                          
c                                                                               
 9010 format(//2x,'calculations for ',a30)                                      
c                                                                               
 9020 format(/3x,'wall time (secs): ',f10.4,                                    
     &       f5.1,' (%) ', 'no. calls: ',i8)                                    
c                                                                               
      return                                                                    
      end                                                                       
