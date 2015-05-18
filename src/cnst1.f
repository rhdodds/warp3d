c     **************************************************************** 
c     *                                                              *
c     *                      subroutine cnst1                        *
c     *                                                              *
c     *                       written by : bh, rhd                   *
c     *                                                              *
c     *                   last modified : 06/29/91, 10/8/97          *
c     *                                                              *
c     *      this subroutine computes the consistent                 *
c     *      tangent operator matrix for simple rate independent     *
c     *      mises plasticity with constant (mixed) hardening        *
c     *      all elements in the block are processed in full         *
c     *      vector mode. this is our fastest material model.        *
c     *                                                              *
c     ****************************************************************
c
c           
      subroutine cnst1( span, cep, rtsg, nu, e, kn1, hprime,
     &                  beta, ldt, dj, w, dstates, felem )
      implicit integer (a-z)
$add param_def        
c
c                       parameter declarations
c
#dbl      double precision
#sgl      real
     &   cep(mxvl,nstr,*), rtsg(mxvl,*), e(*), nu(*), kn1(*),
     &   hprime(*), beta(*), ldt(*), dstates(*), dj(*), w
c
c                      locally allocated arrays
#dbl      double precision
#sgl      real
     &   bb(mxvl), albar(mxvl),
     &   thbar(mxvl), c1(mxvl), c2(mxvl), c3(mxvl), c4(mxvl), g(mxvl),
     &   l(mxvl), k(mxvl), mrtsq(mxvl), gambar(mxvl), gamma(mxvl),
     &   gbar(mxvl), root2, zero, one, two, three, dword
c
      logical yield(mxvl), local_debug
      integer iword(2)
      equivalence ( dword, iword )
c
      data zero, one, two, three, root2 / 0.0, 1.0, 2.0, 3.0,
     &      1.414213562373095d0/
c
c
      local_debug = .false.
c
      do i = 1, span
       dword   = dstates(i)
       iestate = iword(1)
       if( iestate .eq. 1 ) then
          yield(i) = .true.
       else
          yield(i) = .false.
       end if
      end do
c
      if ( local_debug ) then
        write(*,*) '>>> yield flags...'
        write(*,*) i, (yield(j),j=1,span)
        do i = 1, span
         write(*,fmt='(2x,i4,6f10.3)') i+felem-1,(rtsg(i,j),j=1,6)
         write(*,fmt='(2x,i4,3f10.3)') i+felem-1, e(i),
     &                                 kn1(i), hprime(i) 
        end do
      end if
c
      do i = 1, span
c           
c                       branch on whether or not the material
c                       has yielded.
c
         if( .not. yield(i) ) then
c
c                       linear isotropic elastic matrix. there
c                       is no yield.
c
            cep(i,1,4)= zero
            cep(i,1,5)= zero
            cep(i,1,6)= zero
            cep(i,2,4)= zero
            cep(i,2,5)= zero
            cep(i,2,6)= zero
            cep(i,3,4)= zero
            cep(i,3,5)= zero
            cep(i,3,6)= zero
            cep(i,4,1)= zero
            cep(i,4,2)= zero
            cep(i,4,3)= zero
            cep(i,4,5)= zero
            cep(i,4,6)= zero
            cep(i,5,1)= zero
            cep(i,5,2)= zero
            cep(i,5,3)= zero
            cep(i,5,4)= zero
            cep(i,5,6)= zero
            cep(i,6,1)= zero
            cep(i,6,2)= zero
            cep(i,6,3)= zero
            cep(i,6,4)= zero
            cep(i,6,5)= zero
            c1(i)= (e(i)/((one+nu(i))*(one-two*nu(i))))*dj(i)*w
            c2(i)= (one-nu(i))*c1(i)   
            c3(i)= ((one-two*nu(i))/two)*c1(i)
            c4(i)= nu(i)*c1(i)
            cep(i,1,1)= c2(i)
            cep(i,2,2)= c2(i)
            cep(i,3,3)= c2(i)
            cep(i,4,4)= c3(i)
            cep(i,5,5)= c3(i)
            cep(i,6,6)= c3(i)
            cep(i,1,2)= c4(i)
            cep(i,1,3)= c4(i)
            cep(i,2,1)= c4(i)
            cep(i,3,1)= c4(i)
            cep(i,2,3)= c4(i)
            cep(i,3,2)= c4(i)
         else
c
c                       the material point is in the plastic
c                       range. compute the consistent constituitive
c                       matrix.
c
            g(i) = e(i)/(two*(one+nu(i)))
            l(i) = (e(i)*nu(i))/((one+nu(i))*(one-two*nu(i)))
            k(i) = (three*l(i)+two*g(i))/three                         
            mrtsq(i) = rtsg(i,1)**2+rtsg(i,2)**2+rtsg(i,3)**2+two*
     &                (rtsg(i,4)**2+rtsg(i,5)**2+rtsg(i,6)**2)
            bb(i) = (root2*kn1(i)+(two/three)*(one-beta(i))*hprime(i)
     &              *ldt(i))/sqrt(mrtsq(i))
            gamma(i)= one/(one+hprime(i)/(three*g(i)))
            gambar(i)=  gamma(i)-one+bb(i)
            gbar(i)= g(i)*bb(i)
            albar(i)= k(i)- two*gbar(i)/three
            thbar(i)= two*g(i)*gambar(i)
            cep(i,1,1)=  (albar(i)+two*gbar(i)-thbar(i)*(rtsg(i,1)**2)/
     &                    mrtsq(i))*dj(i)*w
            cep(i,2,2)=  (albar(i)+two*gbar(i)-thbar(i)*(rtsg(i,2)**2)/
     &                    mrtsq(i))*dj(i)*w
            cep(i,3,3)=  (albar(i)+two*gbar(i)-thbar(i)*(rtsg(i,3)**2)/
     &                    mrtsq(i))*dj(i)*w
            cep(i,4,4)=  (gbar(i)-thbar(i)*(rtsg(i,4)**2)/mrtsq(i))*
     &                    dj(i)*w
            cep(i,5,5)=  (gbar(i)-thbar(i)*(rtsg(i,5)**2)/mrtsq(i))*
     &                    dj(i)*w
            cep(i,6,6)=  (gbar(i)-thbar(i)*(rtsg(i,6)**2)/mrtsq(i))*
     &                    dj(i)*w
            cep(i,2,1)=  (albar(i)-thbar(i)*rtsg(i,1)*rtsg(i,2)/
     &                    mrtsq(i))*dj(i)*w
            cep(i,3,1)=  (albar(i)-thbar(i)*rtsg(i,1)*rtsg(i,3)/
     &                    mrtsq(i))*dj(i)*w
            cep(i,4,1)= -(thbar(i)*rtsg(i,1)*rtsg(i,4)/mrtsq(i))*dj(i)*w
            cep(i,5,1)= -(thbar(i)*rtsg(i,1)*rtsg(i,5)/mrtsq(i))*dj(i)*w
            cep(i,6,1)= -(thbar(i)*rtsg(i,1)*rtsg(i,6)/mrtsq(i))*dj(i)*w
            cep(i,3,2)=  (albar(i)-thbar(i)*rtsg(i,3)*rtsg(i,2)/
     &                    mrtsq(i))*dj(i)*w
            cep(i,4,2)= -(thbar(i)*rtsg(i,2)*rtsg(i,4)/mrtsq(i))*dj(i)*w
            cep(i,5,2)= -(thbar(i)*rtsg(i,2)*rtsg(i,5)/mrtsq(i))*dj(i)*w
            cep(i,6,2)= -(thbar(i)*rtsg(i,2)*rtsg(i,6)/mrtsq(i))*dj(i)*w
            cep(i,4,3)= -(thbar(i)*rtsg(i,3)*rtsg(i,4)/mrtsq(i))*dj(i)*w
            cep(i,5,3)= -(thbar(i)*rtsg(i,3)*rtsg(i,5)/mrtsq(i))*dj(i)*w
            cep(i,6,3)= -(thbar(i)*rtsg(i,3)*rtsg(i,6)/mrtsq(i))*dj(i)*w
            cep(i,5,4)= -(thbar(i)*rtsg(i,4)*rtsg(i,5)/mrtsq(i))*dj(i)*w
            cep(i,6,4)= -(thbar(i)*rtsg(i,4)*rtsg(i,6)/mrtsq(i))*dj(i)*w
            cep(i,6,5)= -(thbar(i)*rtsg(i,5)*rtsg(i,6)/mrtsq(i))*dj(i)*w
c                  
            cep(i,1,2)=cep(i,2,1)
            cep(i,1,3)=cep(i,3,1)
            cep(i,1,4)=cep(i,4,1)
            cep(i,1,5)=cep(i,5,1)
            cep(i,1,6)=cep(i,6,1)
            cep(i,2,3)=cep(i,3,2)
            cep(i,2,4)=cep(i,4,2)
            cep(i,2,5)=cep(i,5,2)
            cep(i,2,6)=cep(i,6,2)
            cep(i,3,4)=cep(i,4,3)
            cep(i,3,5)=cep(i,5,3)
            cep(i,3,6)=cep(i,6,3)
            cep(i,4,5)=cep(i,5,4)
            cep(i,4,6)=cep(i,6,4)
            cep(i,5,6)=cep(i,6,5)
         end if
      end do
      return
      end
