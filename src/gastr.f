
c     ****************************************************************
c     *                                                              *
c     *                      subroutine gastr                        *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 09/22/2015 rhd             *
c     *                                                              *
c     *     gathers element stresses from the global                 *
c     *     stress data structure to a block of similar,             *
c     *     elements for all gauss points.                           *
c     *                                                              *
c     ****************************************************************
c
c           
      subroutine gastr( ml, mg, ngp, nprm, span )
      implicit integer (a-z)
$add param_def
c
c               parameter declarations
c
#dbl      double precision
#sgl      real
     & ml(mxvl,nprm,*), mg(nprm,ngp,*)
c    
      if ( ngp .ne. 8 ) then                            
@!DIR$ LOOP COUNT MAX=27
        do k = 1, ngp
         do  j = 1, nprm
@!DIR$ LOOP COUNT MAX=###  
            do  i = 1, span
               ml(i,j,k) = mg(j,k,i)
            end do
         end do
        end do
        return
      end if
c
c                number of gauss points = 8, unroll.
      do  j = 1, nprm
@!DIR$ LOOP COUNT MAX=###  
        do  i = 1, span
            ml(i,j,1) = mg(j,1,i)
            ml(i,j,2) = mg(j,2,i)
            ml(i,j,3) = mg(j,3,i)
            ml(i,j,4) = mg(j,4,i)
            ml(i,j,5) = mg(j,5,i)
            ml(i,j,6) = mg(j,6,i)
            ml(i,j,7) = mg(j,7,i)
            ml(i,j,8) = mg(j,8,i)
        end do
      end do
c
      return
      end
