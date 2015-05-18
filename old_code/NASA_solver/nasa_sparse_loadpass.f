      subroutine sparse_loadpass( neqns, xlnz, maxlnz, coefs,
     &                            diag, rhs, iloc, inum  )
#dbl      implicit double precision (a-h,o-z)
#sgl      implicit real (a-h,o-z)
      integer    xlnz(*)
      dimension   coefs(*), rhs(*), diag(*), inum(*), iloc(*)
      data zero / 0.0 /
c
      neq = neqns
c
      do j = 1, neq
        rhsj = rhs(j)/diag(j)
        rhsj = rhs(j)
        rhs(j) = rhsj
        istrt = xlnz(j)
        istop = xlnz(j+1) - 1
        if ( istop .lt. istrt )  go to 202
        i = inum(j)
c$dir no_recurrence
        do ii = istrt, istop
          isub = iloc(i)
          rhs(isub) = rhs(isub) - coefs(ii)*rhsj
          i = i + 1
        end do
202     continue
      end do 
c 
c ---- backsolve
c
      rhs(neq) = rhs(neq) / diag(neq)
      icont = maxlnz
      do i = neq-1, 1, -1
        sum = zero
        mm = inum(i) + xlnz(i+1) - xlnz(i)
        jcont = mm-1
        do j = xlnz(i), xlnz(i+1)-1
         xx = coefs(icont)
         nn = iloc(jcont)
         icont = icont - 1
         jcont = jcont - 1
         sum = sum + xx * rhs(nn) * diag(i)
        end do
        rhs(i) = (rhs(i) - sum) / diag(i)
      end do
c
      return
      end
 







