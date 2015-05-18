      subroutine sparse_check( neq, irow, icoln, maxan, ncof, maxn )
      implicit integer (a-z)
c
      dimension irow(*), icoln(*), maxan(*)
      allocatable jtemp(:)
c 
c              nt number of terms in adjncy
c
      allocate( jtemp(neq+1) )
      jtemp(1:neq) = 0
      icont = 1
c
      do i = 1, neq
        do j = 1, irow(i)
           iix = icoln(icont)
           jtemp(iix) = jtemp(iix) + 1
           jtemp(i) = jtemp(i) + 1
           icont = icont + 1
        end do
      end do
c
      jmax = 0
      maxan(1) = 1
      do i = 1, neq
       jmax = max(jmax,jtemp(i))
       maxan(i+1) = maxan(i) + jtemp(i)
      end do
      maxn = jmax + 5
      ncof = maxan(neq+1)
c
      deallocate( jtemp )
c
      return
      end
