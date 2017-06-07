c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine intab                        *          
c     *                                                              *          
c     *                       written by : jcs                       *          
c     *                                                              *          
c     *                   last modified : 12/01/11  jcs              *          
c     *                                                              *          
c     *     this subroutine supervises and conducts the input of the *          
c     *     tables defined for the current problem                   *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
c                                                                               
c                                                                               
      subroutine intab( sbflg1, sbflg2, tabnum, path )                          
      use global_data ! old common.main
      use main_data, only : tables                                              
      implicit integer (a-z)                                                    
c                                                                               
c                       parameter declarations                                  
c                                                                               
      logical sbflg1, sbflg2                                                    
c                                                                               
c                       local declarations                                      
c                                                                               
      double precision :: forval, dumd                                          
      character :: name*80, lname*24, stlnam*8, dums*1, curtyp*4                
      logical, external :: matchs, endcrd, label, scanms, numi                  
      logical :: debug, complete                                                
      real :: dumr                                                              
      data debug / .false. /                                                    
c                                                                               
      if ( debug ) write (*,*) ' >>> inside intab'                              
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *               get the name of the table being input.               *        
c *               each table must have a name for id                   *        
c *               purposes. if there is no name, or if the vector      *        
c *               containing the table names is full, print an         *        
c *               error message and read another high level            *        
c *               command.                                             *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
      if( label(dummy) ) then                                                   
c                                                                               
         name  = ' '                                                            
         lname = ' '                                                            
         call entits(name,nc)                                                   
         if( nc .gt. 24 ) nc = 24                                               
         lname(1:nc) = name(1:nc)                                               
         if ( debug ) write(*,*) 'Table name = ', lname                         
c                                                                               
c                 determine the location of the specified table                 
c                 in the tables object.                                         
c                                                                               
         do tabn = 1,max_tables                                                 
            if ( debug ) write(*,*) 'numtab', tabn,                             
     &           tables(tabn)%num_rows,                                         
     &           tables(tabn)%table_name                                        
c                                                                               
c                 find the first open slot in the table library                 
c                 and assign it to the specified table.                         
c                 intab only reaches this point after checking                  
c                 the previously defined tables.                                
c                                                                               
            if ( tables(tabn)%num_rows .eq. 0 ) then                            
               tabnum = tabn                                                    
               tables(tabnum)%table_name = lname                                
               go to 810                                                        
            end if                                                              
c                                                                               
c                 check if table at tabn has been defined                       
c                 for this table name previously. if so,                        
c                 use tabn as the table number and move to                      
c                 table definitions.                                            
c                                                                               
            if ( scanms(tables(tabn)%table_name,lname,24) ) then                
               call errmsg(335,param,dums,dumr,dumd)                            
               tabnum = tabn                                                    
               go to 810                                                        
            end if                                                              
c                                                                               
c                 table name at tabn does not have the same name                
c                 as the inputted table. cycle.                                 
c                                                                               
         end do                                                                 
c                                                                               
c                 maximum number of tables has been defined, error.             
c                                                                               
         call errmsg(336,param,dums,dumr,dumd)                                  
         go to 9996                                                             
c                                                                               
      else                                                                      
c                                                                               
         call errmsg(54,dum,dums,dumr,dumd)                                     
         go to 9996                                                             
c                                                                               
      end if                                                                    
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *               look for basic table commands.                       *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
 810  continue                                                                  
      if ( debug ) write(*,*) 'table name established.'                         
      nrows   = 0                                                               
      ncols   = 0                                                               
      tabmark = -1                                                              
 850  continue                                                                  
c                                                                               
c                 scan for high level command                                   
c                                                                               
      if (matchs('rows',4)) go to 851                                           
      if (matchs('type',4)) go to 852                                           
      if (endcrd(dum))      go to 900                                           
c                                                                               
c                 there is no match with existing commands.                     
c                                                                               
      call errmsg(343,dum,dums,dumr,dumd)                                       
      call scan                                                                 
      go to 850                                                                 
c                                                                               
c                 get the number of rows for this table, must be                
c                 greater than 1                                                
c                                                                               
 851  continue                                                                  
c                                                                               
      if ( debug ) write(*,*) 'determining nrows'                               
      if ( .not. numi(ival) ) then                                              
         call errmsg(103,dum,dums,dumr,dumd)                                    
         call scan                                                              
         go to 850                                                              
      else if ( ival .lt. 2 ) then                                              
         if ( debug ) write(*,*) ival                                           
         call errmsg(339,dum,dums,dumr,dumd)                                    
         go to 850                                                              
      else                                                                      
         nrows = ival                                                           
         go to 850                                                              
      end if                                                                    
                                                                                
c                                                                               
c                 determine the table type                                      
c                                                                               
 852  continue                                                                  
c                                                                               
      if ( debug ) write(*,*) 'determing table type'                            
      if ( matchs('piston',4) ) then                                            
         tabmark = 1                                                            
         ncols   = 8                                                            
         go to 850                                                              
      end if                                                                    
c                                                                               
      call errmsg(338,dum,dums,dumr,dumd)                                       
      go to 9996                                                                
c                                                                               
c                 check that table type and the number of rows for this         
c                 table have been stored                                        
c                                                                               
 900  continue                                                                  
c                                                                               
      if ( debug ) write(*,*) 'nrows, tabmark', nrows, tabmark                  
      if ( ( nrows .eq. 0 ) .or. ( tabmark .eq. -1 ) ) then                     
         call errmsg(337,dum,dums,dumr,dumd)                                    
         go to 9996                                                             
      end if                                                                    
c                                                                               
c                 allocate the data structure                                   
c                                                                               
      call do_table_allo( tabnum, nrows, ncols, .true. )                        
c                                                                               
c                 branch list for loading types                                 
c                                                                               
      if ( tabmark .eq. 1 ) then                                                
         tables(tabnum)%table_type = 'PISTON  '                                 
         call table_piston( tabnum, nrows, complete )                           
         if ( complete ) go to 9997                                             
         call errmsg(349,dum,dums,dumr,dumd)                                    
         call die_gracefully                                                    
         stop                                                                   
      end if                                                                    
c                                                                               
c                                                                               
c **********************************************************************        
c **********************************************************************        
c                                                                               
c                                                                               
 9995 sbflg1= .true.                                                            
      sbflg2= .true.                                                            
      if (debug) write (*,*) '9995'                                             
      go to 9999                                                                
c                                                                               
 9996 sbflg1= .false.                                                           
      sbflg2= .true.                                                            
      if (debug) write (*,*) '9996'                                             
      go to 9999                                                                
c                                                                               
 9997 sbflg1= .false.                                                           
      sbflg2= .false.                                                           
      if (debug) write (*,*) '9997'                                             
      go to 9999                                                                
c                                                                               
 9998 sbflg1= .true.                                                            
      sbflg2= .false.                                                           
      if (debug) write (*,*) '9998'                                             
c                                                                               
c                                                                               
 9999 continue                                                                  
      if ( debug ) write (*,*) ' >>> leaving intab'                             
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine do_table_allo                *          
c     *                                                              *          
c     *                       written by : jcs                       *          
c     *                                                              *          
c     *     this subroutine does the allocation and deallocation     *          
c     *     of the permanent table input structures                  *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine do_table_allo( tabnum, nrows, ncols, fill )                    
      use global_data ! old common.main
      use main_data, only : tables                                              
      implicit integer (a-z)                                                    
c                                                                               
      real rzero                                                                
      double precision dzero                                                    
      logical fill, debug                                                       
c                                                                               
      data zero, debug, dzero / 0.0, .false., 0.0d00 /                          
c                                                                               
      if ( debug ) write (*,*) ' >>> inside do_table_allo'                      
      if ( debug ) write (*,*) tabnum,  tables(tabnum)%num_rows,                
     &     nrows, ncols                                                         
c                                                                               
c                 Allocation of variables.  First check if the general          
c                 tab_defn structure has been allocated; allocate it            
c                 if needed.  Then set the size of the variables                
c                 from nrow and allocate the loading structures                 
c                 inside tables for the current table definition                
c                 number.                                                       
c                                                                               
      if (.not. allocated(tables)) then                                         
         write (*,*) '     tables vector not allocated. stop'                   
         call die_gracefully                                                    
         stop                                                                   
      endif                                                                     
c                                                                               
      if (nrows.eq.0) then                                                      
         tables(tabnum)%num_rows = 0                                            
         tables(tabnum)%num_cols = 0                                            
         if (debug) write (*,*) '     skipping allocation:'                     
         goto 9999                                                              
      endif                                                                     
c                                                                               
      if (tables(tabnum)%num_rows .ne. 0) then                                  
         if ( debug ) write(*,*) '  deallocating data'                          
         deallocate (tables(tabnum)%table_values_sgl)                           
         deallocate (tables(tabnum)%table_values_dbl)                           
      endif                                                                     
c                                                                               
      if ( debug ) write(*,*) '     allocating data'                            
      allocate (tables(tabnum)%table_values_sgl(nrows,ncols))                   
      allocate (tables(tabnum)%table_values_dbl(nrows,ncols))                   
      tables(tabnum)%num_rows = nrows                                           
      tables(tabnum)%num_cols = ncols                                           
c                                                                               
c      if ( debug ) write(*,*) tables(tabnum)%num_rows,                         
c     &     tables(tabnum)%num_cols                                             
c      if ( debug ) then                                                        
c         do row = 1,nrows                                                      
c            do col = 1,ncols                                                   
c               write(*,*),                                                     
c     &              row, col, tables(tabnum)%table_values_sgl(row,col)         
c            end do                                                             
c         end do                                                                
c      end if                                                                   
c                                                                               
c                 initialize tables vectors                                     
c                                                                               
      if ( fill ) then                                                          
         do row = 1,nrows                                                       
            do col = 1,ncols                                                    
               tables(tabnum)%table_values_sgl(row,col) = zero                  
               tables(tabnum)%table_values_dbl(row,col) = dzero                 
            end do                                                              
         end do                                                                 
      end if                                                                    
c                                                                               
      if ( debug) write (*,*) ' >>> leaving do_table_allo'                      
c                                                                               
 9999 continue                                                                  
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine table_piston                 *          
c     *                                                              *          
c     *                       written by : jcs                       *          
c     *                                                              *          
c     *     this subroutine reads the definition of a piston type    *          
c     *     table given defined by the user                          *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine table_piston( tabnum, nrows, complete )                        
      use global_data ! old common.main
      use main_data, only : tables                                              
      implicit integer (a-z)                                                    
c                                                                               
c                       parameter declarations                                  
c                                                                               
      logical complete                                                          
c                                                                               
c                       local declarations                                      
c                                                                               
      double precision                                                          
     &  dumd                                                                    
      character name*80, lname*8, stlnam*8, dums*1, curtyp*4                    
      logical matchs, endcrd, true, numr, scanms, debug,                        
     &     pist_def(8), consistent                                              
      dimension pist_order(8)                                                   
      real dumr, rwtime, rwtime_prv, zero, pist_load(8), pist_val,              
     &  tnm1, tn, mx, my, mz, norm                                              
      data debug, zero / .false., 0.0 /                                         
c                                                                               
      if ( debug ) write (*,*) ' >>> inside table_piston'                       
c                                                                               
c                 initialize variables                                          
c                                                                               
      complete = .false.                                                        
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *               look for column ordering.                            *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
      call readsc                                                               
 1000 continue                                                                  
c                                                                               
      if( matchs('columns',3) ) call splunj                                     
      pnow = 1                                                                  
      pist_def(1:8) = .false.                                                   
c                                                                               
 1100 continue                                                                  
c                                                                               
c                 branch on column name                                         
c                                                                               
      if ( matchs('time',4) ) then                                              
         pcol = 1                                                               
         go to 1200                                                             
      else if ( matchs('flowpressure',4) ) then                                 
         pcol = 2                                                               
         go to 1200                                                             
      else if ( matchs('velocity',3) ) then                                     
         pcol = 3                                                               
         go to 1200                                                             
      else if ( matchs('mach',4) ) then                                         
         pcol = 4                                                               
         go to 1200                                                             
      else if ( matchs('isentropic',5) ) then                                   
         pcol = 5                                                               
         go to 1200                                                             
      else if ( matchs('xdirection',4) ) then                                   
         pcol = 6                                                               
         go to 1200                                                             
      else if ( matchs('ydirection',4) ) then                                   
         pcol = 7                                                               
         go to 1200                                                             
      else if ( matchs('zdirection',4) ) then                                   
         pcol = 8                                                               
         go to 1200                                                             
      else if ( endcrd(dum) ) then                                              
         go to 1300                                                             
      else                                                                      
         go to 1400                                                             
      end if                                                                    
c                                                                               
c                 column type known                                             
c                                                                               
 1200 continue                                                                  
c                                                                               
c                 check if number of inputted columns greater                   
c                 than max number allowed for a piston loading                  
c                                                                               
      if ( debug ) write(*,*) '(pnow,pcol)', pnow, pcol                         
      if ( debug ) write (*,*) pist_def(1:8), all( pist_def )                   
c                                                                               
      if ( pnow .gt. 8 ) then                                                   
         call errmsg(345,dum,dums,dumr,dumd)                                    
         if ( all( pist_def ) ) go to 2000                                      
         go to 1900                                                             
      end if                                                                    
c                                                                               
c                 set order in pist_order                                       
c                                                                               
      pist_order(pnow) = pcol                                                   
      pist_def(pcol) = .true.                                                   
      pnow = pnow + 1                                                           
      go to 1100                                                                
c                                                                               
c                 reached end of line, check that all columns                   
c                 have been defined                                             
c                                                                               
 1300 continue                                                                  
c                                                                               
      if ( .not. all( pist_def ) ) go to 1900                                   
      go to 2000                                                                
c                                                                               
c                                                                               
c                 column type unknown, error                                    
c                                                                               
 1400 continue                                                                  
c                                                                               
      call errmsg(344,dum,dums,dumr,dumd)                                       
      call scan                                                                 
      go to 1100                                                                
c                                                                               
c                 not all columns have been input at the end of the             
c                 column definition line, output error and return               
c                 to calling function as an incomplete table.                   
c                                                                               
 1900 continue                                                                  
c                                                                               
      call errmsg(346,dum,dums,dumr,dumd)                                       
      go to 9999                                                                
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *               read in piston table data.                           *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
 2000 continue                                                                  
c                                                                               
      if ( debug ) write(*,*) pist_order(1:8)                                   
      row = 0                                                                   
      consistent = .true.                                                       
c                                                                               
c                 start loop to read in each row of table, check                
c                 inputted data while reading.                                  
c                                                                               
 2100 continue                                                                  
c                                                                               
      row = row + 1                                                             
      if ( row .gt. nrows ) go to 3000                                          
      if ( debug ) write(*,*) 'Reading in row #:', row                          
c                                                                               
      call readsc                                                               
      pist_def(1:8) = .false.                                                   
      col = 1                                                                   
c                                                                               
 2200 continue                                                                  
c                                                                               
c                 read in unformatted data in row/column form.                  
c                 store column values in appropriate locations                  
c                 defined for piston tables in table module.                    
c                                                                               
      if ( .not. numr(pist_val) ) then                                          
c                                                                               
         if (endcrd(dum)) then                                                  
            if ( debug ) write(*,*) pist_def, all( pist_def )                   
            if ( debug ) write(*,*)                                             
     &           row, (tables(tabnum)%table_values_sgl(row,k),k=1,8)            
            if ( all(pist_def) ) go to 2100                                     
            call errmsg(348,row,dums,dumr,dumd)                                 
            consistent = .false.                                                
            go to 2100                                                          
         else                                                                   
            call errmsg(18,dum,dums,dumr,dumd)                                  
            if(true(dum)) go to 2200                                            
         end if                                                                 
c                                                                               
      else                                                                      
c                                                                               
         if ( col .le. 8 ) then                                                 
            pcol = pist_order(col)                                              
            if ( debug ) write(*,*) 'col, pcol, pist_val',                      
     &           col, pcol, pist_val                                            
            tables(tabnum)%table_values_sgl(row,pcol) = pist_val                
            pist_def(pcol) = .true.                                             
            if ( debug ) write(*,*) 'col, pcol, pist_val',                      
     &           col, pcol, pist_val                                            
            col = col + 1                                                       
            go to 2200                                                          
         else                                                                   
            call errmsg(347,dum,dums,dumr,dumd)                                 
            go to 2100                                                          
         end if                                                                 
c                                                                               
      end if                                                                    
c                                                                               
c                                                                               
c **********************************************************************        
c *                                                                    *        
c *               consistency checks for piston data:                  *        
c *               1. time must start at zero                           *        
c *               2. time must be monotonically increasing             *        
c *               3. flow direction vector must have a magnitude       *        
c *                  equal to 1. scale vector if needed.               *        
c *                                                                    *        
c **********************************************************************        
c                                                                               
c                                                                               
 3000 continue                                                                  
c                                                                               
      tnm1 = tables(tabnum)%table_values_sgl(1,1)                               
      if ( tnm1 .ne. zero ) then                                                
         call errmsg(342,dum,dums,dumr,dumd)                                    
         consistent = .false.                                                   
      end if                                                                    
c                                                                               
      do i = 2, nrows                                                           
c                                                                               
         tn = tables(tabnum)%table_values_sgl(i,1)                              
         if ( tnm1 .ge. tn ) then                                               
            call errmsg(341,dum,dums,dumr,dumd)                                 
            consistent = .false.                                                
         end if                                                                 
         tnm1 = tn                                                              
c                                                                               
      end do                                                                    
c                                                                               
      do i = 1, nrows                                                           
c                                                                               
         mx = tables(tabnum)%table_values_sgl(i,6)                              
         my = tables(tabnum)%table_values_sgl(i,7)                              
         mz = tables(tabnum)%table_values_sgl(i,8)                              
c                                                                               
         norm = sqrt( mx*mx + my*my + mz*mz )                                   
         if ( norm .eq. zero ) then                                             
            call errmsg(340,dum,dums,dumr,dumd)                                 
            consistent = .false.                                                
         end if                                                                 
c                                                                               
         tables(tabnum)%table_values_sgl(i,6) = mx/norm                         
         tables(tabnum)%table_values_sgl(i,7) = my/norm                         
         tables(tabnum)%table_values_sgl(i,8) = mz/norm                         
c                                                                               
      end do                                                                    
c                                                                               
      if ( .not. consistent ) go to 9999                                        
c                                                                               
c                 table passes all internal consistentcy checks.                
c                 mark table definition as complete and return                  
c                 to calling function.                                          
c                                                                               
      complete = .true.                                                         
      if ( debug ) call dump_table_piston( tabnum, nrows )                      
c                                                                               
 9999 continue                                                                  
      if ( debug ) write (*,*) ' >>> leaving table_piston'                      
      return                                                                    
      end                                                                       
c                                                                               
c     ****************************************************************          
c     *                                                              *          
c     *                      subroutine dump_table                   *          
c     *                                                              *          
c     *                       written by : jcs                       *          
c     *                                                              *          
c     *          this subroutine prints out table entries            *          
c     *                                                              *          
c     ****************************************************************          
c                                                                               
      subroutine dump_table_piston( tabnum, nrows )                             
      use global_data ! old common.main
      use main_data, only : tables                                              
      implicit integer (a-z)                                                    
c                                                                               
      write (*,*) ' ==================='                                        
      write (*,*) '    DUMPING TABLE   '                                        
      write (*,*) ' ==================='                                        
c                                                                               
c                                                                               
      write(*,1000)                                                             
      do i = 1,nrows                                                            
         write(*,1010) i, (tables(tabnum)%table_values_sgl(i,j),j=1,8)          
      end do                                                                    
c                                                                               
      return                                                                    
 1000 format( '  row|         time|   flow pres.|    flow vel.',                
     &             '|  mach number| isent. expo.| flw drc. (x)',                
     &             '| flw drc. (y)| flw drc. (z) ',/,                           
     &        '===============================================',                
     &             '==========================================',                
     &             '=============================')                             
 1010 format( i5, '|', e13.6, '|', e13.6, '|', e13.6, '|', e13.6,               
     &            '|', e13.6, '|', e13.6, '|', e13.6, '|', e13.6 )              
      end                                                                       
