c     ****************************************************************
c     *                                                              *
c     *                      subroutine ouhnod                       *
c     *                                                              *
c     *                       written by : bh                        *
c     *                                                              *
c     *                   modified :    12/1/2017 rhd                *
c     *                                                              *
c     *     hardcopy or packet output for the                        *
c     *     displacements, velocities, accelerations, reaction       *
c     *     forces, temperatures for the given node.                 *
c     *                                                              *
c     ****************************************************************
c
c
c
      subroutine ouhnod( dva, node, lnum, pgnum, doflbl, wide, eform,
     &                   prec, ndof, hedtyp, elem, newel,
     &                   react_sums, noheader, react_totals_only,
     &                   write_to_packet, lsttyp )
      use global_data ! old common.main
c
      use main_data, only : trn, trnmat, mdiag, pbar, rload,
     &                      packet_file_no, temper_nodes,
     &                      temper_nodes_ref
      use stiffness_data, only : total_lagrange_forces
c
      implicit none
c
      integer :: node, pgnum, ndof, elem, dva, lsttyp
      double precision :: react_sums(*)
      logical :: wide, eform, prec, newel, noheader,
     &           react_totals_only, write_to_packet
      character(len=8) :: doflbl(*)
      character(len=20) :: hedtyp
c
c                       local declarations
c
      integer :: i, dof, sdof, fnsh, nitm, hfmtyp, fmtyp, nl,
     &           nvalues_out, lnum, cl, strt
      double precision ::
     &     edva(mxvl,mxndof), nfac, one, zero, force_lag,
     &     trnmte(mxvl,mxedof,mxndof)
      logical :: newhed, trne(mxvl,mxndel), hardcopy
      data one, zero / 1.0d00, 0.0d00 /
c
c                       dva = 1 -> 5 for displacements,
c                       velocities, accelerations, reactions,
c                       temperature
c
c                       set the number of lines to be output per node.
c
      hardcopy = .not. write_to_packet
      nitm = 0
c
      if( hardcopy ) then
        if( wide .and. prec ) then
           nitm    = 4
           hfmtyp  = 1
           if(eform) then
              fmtyp = 1
           else
              fmtyp = 2
           end if
        else if( .not. wide .and. prec ) then
           nitm   = 2
           hfmtyp = 2
           if(eform) then
              fmtyp = 3
           else
              fmtyp = 4
           end if
        else if( wide .and. .not. prec ) then
           nitm   = 6
           hfmtyp = 3
           if( eform ) then
            fmtyp = 5
           else
              fmtyp = 6
           end if
        else if( .not. wide .and. .not. prec ) then
           nitm   = 3
           hfmtyp = 4
           if( eform ) then
              fmtyp = 7
           else
              fmtyp = 8
           end if
        end if
c
c                       compute number of output lines needed
c                       for above printing requirements
        nl = ndof / nitm
        if( dva .eq. 5 ) nl = 1 ! only print temperature at node
        if( (ndof-nl*nitm) .gt. 0 ) nl = nl + 1
c
c                       print page headers if needed
c
        if ( .not. noheader ) then
          call oundhd( pgnum, lnum, hedtyp, doflbl, ndof, newhed,
     &                 hfmtyp, nitm, nl )
        end if
c
c                       check to see if node is from a new element,
c                       and if it is, write a message. applicable
c                       only for element directed output (i.e. nodes
c                       connected to element in user-specified list).
c                       write only to main output file, not to
c                       packet file
c
        if( newel ) then
           write(out,910) elem
           newel = .false.
        end if
c
      end if  ! hardcopy if
c
c                       set the vector of dis/vel/acc/reactions/temps
c                       to be printed. for reactions, it is a bit
c                       complicated since we have to remove internal
c                       forces and inertia effects. nfac is from
c                       newmark beta integration scheme.
c
      if( dva .le. 4 ) then ! not temps
         nfac = one / (nbeta*dt*dt)
         do dof = 1, ndof
            sdof = dstmap(node)+dof-1
            if( dva .eq. 1 ) then
               edva(1,dof) = u(sdof)
            else if( dva .eq. 2 ) then
               edva(1,dof) = v(sdof)
            else if( dva .eq. 3 ) then
               edva(1,dof) = a(sdof)
            else if( dva .eq. 4 ) then
               force_lag = zero
               if( allocated( total_lagrange_forces ) )
     &            force_lag =  total_lagrange_forces(sdof)
               edva(1,dof) = pbar(sdof) - ifv(sdof) -
     &                    mdiag(sdof)*du(sdof)*nfac +
     &               force_lag
               if( cstmap(sdof) .ne. 0 ) then
                edva(1,dof) = -rload(sdof)+ifv(sdof)+
     &                        mdiag(sdof)*a(sdof)
               end if
            end if
         end do
c
c                       transform dis/vel/acc/reaction at the
c                       current node to uniform global coordinates.
c
         trne(1,1) = trn(node)
         if( trne(1,1) ) then
           trnmte(1,1,1) = trnmat(node)%mat(1,1)
           trnmte(1,1,2) = trnmat(node)%mat(1,2)
           trnmte(1,1,3) = trnmat(node)%mat(1,3)
           trnmte(1,2,1) = trnmat(node)%mat(2,1)
           trnmte(1,2,2) = trnmat(node)%mat(2,2)
           trnmte(1,2,3) = trnmat(node)%mat(2,3)
           trnmte(1,3,1) = trnmat(node)%mat(3,1)
           trnmte(1,3,2) = trnmat(node)%mat(3,2)
           trnmte(1,3,3) = trnmat(node)%mat(3,3)
           call trnvec( edva, trnmte, trne, ndof, 1, 1, 2 )
         end if
c
c                       for reaction force output, we accumulate
c                       a sum for each direction at the requested
c                       nodes. totals are output at completion
c                       of request.
c
         if ( dva .eq. 4 ) then
           react_sums(1) = react_sums(1) + edva(1,1)
           react_sums(2) = react_sums(2) + edva(1,2)
           react_sums(3) = react_sums(3) + edva(1,3)
         end if
      end if
c
      nvalues_out = ndof
      if( dva .eq. 5 ) then   ! temperature only
         edva(1,1) = temper_nodes(node)
         edva(1,2) = temper_nodes(node) - temper_nodes_ref(node)
         nvalues_out = 2
      end if
c
c                       print the first line of nodal d/v/a/r/T output.
c
      if( nitm .gt. ndof ) then
         fnsh = ndof
      else
         fnsh = nitm
      end if
      if( dva .eq. 5 ) fnsh = 2   ! temperature only
c
c                       if printing reactions, user can request only
c                       printing of just totals -- then skip node
c                       values. output is directed to the packet
c                       file when needed
c
      if ( dva .eq. 4 .and. react_totals_only ) go to 100
c
      if ( hardcopy ) then
c
         if( fmtyp .eq. 1 ) then
            write(out,915) node,(edva(1,i),i=1,fnsh)
         else if( fmtyp .eq. 2 ) then
            write(out,917) node,(edva(1,i),i=1,fnsh)
         else if( fmtyp .eq. 3 ) then
            write(out,920) node,(edva(1,i),i=1,fnsh)
         else if( fmtyp .eq. 4 ) then
            write(out,922) node,(edva(1,i),i=1,fnsh)
         else if( fmtyp .eq. 5 ) then
            write(out,925) node,(edva(1,i),i=1,fnsh)
         else if( fmtyp .eq. 6 ) then
            write(out,927) node,(edva(1,i),i=1,fnsh)
         else if( fmtyp .eq. 7 ) then
            write(out,930) node,(edva(1,i),i=1,fnsh)
         else if( fmtyp .eq. 8 ) then
            write(out,932) node,(edva(1,i),i=1,fnsh)
         end if
c
      else
c
c                    write to binary packet file (nodal packet
c                    or element packet).
c

         if( lsttyp .eq. 1 ) then
           write( packet_file_no ) node,(edva(1,i),i=1,nvalues_out)
         else
           write( packet_file_no ) elem, node,
     &                             (edva(1,i),i=1,nvalues_out)
         end if
         go to 9999 !  return
c
      end if
c
 100  continue
c
      lnum = lnum + 1
c
c                    print subsequent lines of nodal d/v/a/r/T output
c                    if necessary.
c
      cl = 2
 30   if( cl .gt. nl ) go to 9999 ! return
c
      lnum = lnum + 1
      if ( .not. noheader ) then
        call oundhd( pgnum, lnum, hedtyp, doflbl, ndof, newhed, hfmtyp,
     &               nitm, nl )
        if( newhed ) write(out,945)
      end if
c
      strt = (cl-1)*nitm + 1
      if( cl*nitm .gt. ndof ) then
         fnsh = ndof
      else
         fnsh = cl*nitm
      end if
c
      if ( dva .eq. 4 .and. react_totals_only ) go to 200
c
      if( fmtyp .eq. 1 ) then
            write(out,916) (edva(1,i),i=strt,fnsh)
      else if( fmtyp .eq. 2 ) then
            write(out,918) (edva(1,i),i=strt,fnsh)
      else if( fmtyp .eq. 3 ) then
            write(out,921) (edva(1,i),i=strt,fnsh)
      else if( fmtyp .eq. 4 ) then
            write(out,923) (edva(1,i),i=strt,fnsh)
      else if( fmtyp .eq. 5 ) then
            write(out,926) (edva(1,i),i=strt,fnsh)
      else if( fmtyp .eq. 6 ) then
            write(out,928) (edva(1,i),i=strt,fnsh)
      else if( fmtyp .eq. 7 ) then
            write(out,931) (edva(1,i),i=strt,fnsh)
      else if( fmtyp .eq. 8 ) then
            write(out,933) (edva(1,i),i=strt,fnsh)
      end if
c
 200  continue
      cl = cl + 1
      go to 30
c
 910  format(//1x,'>>>>> element ',i7//)
 915  format(6x,i7,4(3x,e26.19))
 916  format(13x,4(3x,e26.19))
 917  format(6x,i7,4(3x,f26.16))
 918  format(13x,4(3x,f26.16))
c
 920  format(6x,i7,2(3x,e26.19))
 921  format(13x,2(3x,e26.19))
 922  format(6x,i7,2(3x,f26.16))
 923  format(13x,2(3x,f26.16))
c
 925  format(6x,i7,6(8x,e12.5))
 926  format(13x,6(8x,e12.5))
 927  format(6x,i7,6(8x,f12.6))
 928  format(13x,6(8x,f12.6))
c
 930  format(6x,i7,3(8x,e12.5))
 931  format(13x,3(8x,e12.5))
 932  format(6x,i7,3(8x,f12.6))
 933  format(13x,3(8x,f12.6))
c
 945  format(1x,' ')
c
 9999 return
      end










