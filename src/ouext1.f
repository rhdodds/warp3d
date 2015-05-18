c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouext1                       *
c     *                                                              *
c     *                       written by : rhd                       *
c     *                                                              *
c     *                   last modified : 3/11/04 (rhd)              *
c     *                                                              *
c     *     this subroutine computes derived stress/strain values    *
c     *     at gauss points or node points for later output.         *
c     *                                                              *
c     ****************************************************************
c
c
      subroutine ouext1( span, stress, nnode, ngp, nodpts,
     &                   num_short_stress, num_short_strain,
     &                   op_code, mxvl )
      use elblk_data, only : elestr
      implicit integer (a-z)
      logical stress, nodpts
c
c                       find the number of strain points.                     
c                                                                  
      nstrpt = ngp
      if ( nodpts ) nstrpt = nnode
c
c                       compute additional data for output at
c                       each of the element nodes or gauss points 
c                       of the element block.
c
c                       op_code = 1 : compute mises equiv. stress or
c                                     equiv. strain.
c                       op_code = 2 : compute invariants, principal
c                                     values.
c
      if ( op_code .eq. 1 ) then
        do strpt = 1, nstrpt
         if( stress) then
            call ouyld1( span, elestr(1,1,strpt), elestr(1,8,strpt),
     &                   mxvl )
         else
            call oueff1( span, elestr(1,1,strpt), elestr(1,7,strpt),
     &                   mxvl )
         end if
       end do
       return
      end if     
c
c                       compute invariants, principal values.
c
      loc = num_short_strain + 1
      if ( stress ) loc = num_short_stress + 1
c
      do strpt = 1, nstrpt
c
c                       compute the principal invariants of desired 
c                       stress or strain.
c             
         call ouinv1( span, elestr(1,1,strpt), elestr(1,loc,strpt),
     &                stress )
c
c                       compute the desired principal stresses or 
c                       principal strains and the direction
c                       cosines of their corresponding normals.
c                       
         call oupri1( span, elestr(1,1,strpt), elestr(1,loc+3,strpt),
     &                elestr(1,loc+6,strpt), stress )
c
      end do
c
c
      return
      end


