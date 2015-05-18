c
c     Program gen_lmprd
c
c          This program generates fully unrolled versions of the core
c          code of the lmprd subroutine. The basic operation of the 
c          generated code is to first gather the trial displacement vector
c          for an element, to do a matrix multiplication of the element
c          stiffness, and then to scatter the resulting tial load vector to 
c          the global vector. Fully unrolling this subroutine can 
c          significantly improve performance on some platforms.
c          
c          The number of columns in the element stiffness matrix is set 
c          by the num_cols parameter.  Output is to standard output.
c
c
c     = ================================================================
c
      program gen_lmprd
      implicit integer (a-z)
c
      parameter (num_cols = 27)
      parameter (trianum_cols = num_cols * num_cols / 2 + num_cols / 2)
c      
      dimension tops (num_cols)
c
      tops(1) = 1
      do i=2, num_cols
         tops(i) = tops(i-1) + i-1
      enddo
c
c	    output gather
c
      write (*,'(7x,"do i = 1, span")')
      do i=1, num_cols
	 if (i .lt. 10) then
	    write (*,'(10x,"xe",i1,"=vec(edest(",i2,",i))")')i,i
         else
	    write (*,'(10x,"xe",i2,"=vec(edest(",i2,",i))")')i,i
         endif
      enddo
      write (*,'("c ")')
      write (*,'("c ")')
      write (*,'("c ")')
c
c           output matrix vector multiplication
c
      do i=1, num_cols
         if (i .lt. 10) then
            write (*,'("         ye",i1," =")'), i
         else
            write (*,'("         ye",i2," =")'), i
         endif
         write (*,1000)
         do j=1, i-1
            if (j .lt. 10) then
               write (*,1011) tops(i) + j-1, j
            else
               write (*,1010) tops(i) + j-1, j
            endif
            if (mod(j,3).eq.0) then
               write (*,1030)
               write (*,1000)
            endif
         enddo
c     
         do j=i, num_cols
            if (j .eq. num_cols) then
               if (j .lt. 10) then
                  write (*,1021) tops(j) + i-1, j
               else
                  write (*,1020) tops(j) + i-1, j
               endif
               write (*,1040)
            else
               if (j .lt. 10) then
                  write (*,1011) tops(j) + i-1, j
               else
                  write (*,1010) tops(j) + i-1, j
               endif
c     
               if (mod(j,3).eq.0) then
                  write (*,1030)
                  write (*,1000)
               endif
            endif
         enddo
      enddo
c
c	    output scatter
c
      do i=1, num_cols
         if (i .lt. 10) then
            write (*,1051) i,i,i
         else
            write (*,1050) i,i,i
         endif
      enddo
      write (*,'(7x,"enddo")')
      write (*,'("c ")')
      write (*,'("c ")')
      write (*,'("c ")')
c
 1000 format ("     &     ",$)
 1010 format (" cpcm(",i4,",i)*xe",i2," +",$)
 1011 format (" cpcm(",i4,",i)*xe",i1," +",$)
 1020 format (" cpcm(",i4,",i)*xe",i2)
 1021 format (" cpcm(",i4,",i)*xe",i1)
 1030 format (" ")
 1040 format ("c ")
 1050 format (10x,'prd(edest(',i2,',j))=prd(edest(',i2,',j)) ',
     &        '+ ye',i2)
 1051 format (10x,'prd(edest(',i2,',j))=prd(edest(',i2,',j)) ',
     &        '+ ye',i1)
      end
