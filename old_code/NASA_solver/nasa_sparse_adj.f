      subroutine sparse_adj( neq, irow, icoln, ncof,
     &                     maxn, nt, xadj, adjncy, jmax, maxan )
#dbl      implicit double precision (a-h,o-z)
#sgl      implicit real (a-h,o-z)
      dimension irow(*), icoln(*), maxan(*)
      integer adjncy(*), xadj(*)
      allocatable itemp(:), jtemp(:)
c
c             nt number of terms in adjncy
c
      allocate( itemp(ncof) )
      allocate( jtemp(neq) )
c
      do i = 1, neq
        jtemp(i) = 0
      end do
      icont = 1
      ncont = 1
      lcont = 1
      jmax = 0
c
      xadj(1) =  lcont
      do i = 1, neq
      kcont = 0
        do j = 1, jtemp(i)
          adjncy(ncont) = itemp(maxan(i)+j-1)
          ncont = ncont + 1
          kcont = kcont + 1
        end do
        do j = 1, irow(i)
          iix = icoln(icont)
          adjncy(ncont) = icoln(icont)
          ncont = ncont + 1
          itemp(maxan(i)+jtemp(i)) = iix
          itemp(maxan(iix)+jtemp(iix)) = i
          if( jtemp(iix) .ge. maxn ) then
            write(*,*)' increase the dimension of maxn to ....',maxn
            call die_gracefully
            stop
          endif
          jtemp(iix) = jtemp(iix) + 1
          jtemp(i) = jtemp(i) + 1
          icont = icont + 1
          kcont = kcont + 1
        end do
        lcont = lcont + kcont
        xadj(i+1) = lcont
       jmax = max(jmax,jtemp(i))
      end do
c
      nt = lcont - 1
      xadj(neq+1) = lcont
c
      deallocate( itemp )
      deallocate( jtemp )
c
      return
      end
