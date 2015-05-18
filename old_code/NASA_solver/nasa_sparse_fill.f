c****************************************************************
c
      subroutine  sparse_fill( neqns, adjncy, xlnz, maxlnz,
     &                         amat, s_matrix, irow, iloc, inum )
c
c****************************************************************
c
#dbl      implicit double precision (a-h,o-z)
#sgl      implicit real (a-h,o-z)
      integer adjncy(*), xlnz(*), irow(*), iloc(*), inum(*)  
      dimension  amat(*), s_matrix(*)
      data zero / 0.0 /
c
      icont = 1
      do i = 1, maxlnz
       s_matrix(i) = zero
      end do
c
      do i = 1, neqns-1
        if( irow(i) .eq. 0 ) go to 101
        ixx = irow(i)
        kcont = 0
        jcont = inum(i)
        do j = xlnz(i), xlnz(i+1) - 1
              if( iloc(jcont) .eq. adjncy(icont) ) then 
                s_matrix(j) = amat(icont)
                icont = icont + 1
                kcont = kcont + 1
                if( kcont .eq. ixx ) go to 101
              end if
              jcont = jcont + 1
         end do
101      continue
         jcont = jcont + 1
      end do
c
      return
      end
