c
c     Program gen_frwd
c
c          This program generates fully unrolled versions of the core
c          code of the efrwrd subroutine. The basic operation of the 
c          generated code is to first gather the right hand side vector
c          for an element, to do forward solve with a triangular matrix,
c          and then to scatter the resulting left hand vector to the global
c          structure. Fully unrolling this subroutine can significantly 
c          improve performance on some platforms.
c          
c          The number of columns in the element stiffness matrix is set 
c          by the num_cols parameter.  Output is to standard output.
c
c
c     = ================================================================
c
      program gen_frwd
      implicit integer (a-z)
c
      parameter (num_cols = 60)
      parameter (trianum_cols = num_cols * (num_cols+1) / 2  )
c      
      dimension tops (num_cols)
c
c            generate structure to help resolve the triangular storage
c            scheme for the matrix
c
      tops(1) = 1
      do i=2, num_cols
         tops(i) = tops(i-1) + i-1
      enddo
c
c            output gather
c
      write (*,'("      do i = 1, span")')
      do i=1, num_cols
         if (i .lt. 10) then
            write (*,'("         ze",i1," = z(edest(",i1,",i))")') i,i
         else
            write (*,'("         ze",i2," = z(edest(",i2,",i))")') i,i
         endif
      enddo         
c
c            output forward solve
c
      do i = 2, num_cols 
         write (*,1040)
         if (i .lt. 10) then
            write (*,'("         ze",i1," = ze",i1," -")') i,i
         else
            write (*,'("         ze",i2," = ze",i2," -")') i,i
         endif
         write (*,1000)
         do j = 1, i-1
c
            if (j .eq. i-1) then
               if (j .lt. 10) then
                  write (*,1021) tops(i)+j-1,j
               else
                  write (*,1020) tops(i)+j-1,j
               endif
            else
               if (j .lt. 10) then
                  write (*,1011) tops(i)+j-1,j
               else
                  write (*,1010) tops(i)+j-1,j
               endif
               if (mod (j - num_cols-1, 3).eq.0) write (*,1030)
            endif
         enddo
      enddo
c
c             output scatter
c
      do i=1, num_cols
         if (i .lt. 10) then
            write (*,'("         z(edest(",i1,",i)) = ze",i1)') i,i
         else
            write (*,'("         z(edest(",i2,",i)) = ze",i2)') i,i
         endif
      enddo         
      write (*,'("      enddo")')
c
c
c
 1000 format ("     &     ",$)
 1010 format (" pcm(",i4,",i)*ze",i2," -",$)
 1011 format (" pcm(",i4,",i)*ze",i1," -",$)
 1020 format (" pcm(",i4,",i)*ze",i2)
 1021 format (" pcm(",i4,",i)*ze",i1)
 1030 format (" ",/"     &     ",$)
 1040 format ("c ")
      end
