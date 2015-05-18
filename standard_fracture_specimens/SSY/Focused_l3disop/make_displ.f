      implicit double precision (a-h,o-z)
c
c             make SSY boundary displacements for K and T loading.
c             also scale existing coords to a different radius
c             if desired.
c
c
      dimension coord(3,5000), node_list(5000)
      double precision new_radius
      integer elem, outfile
      character * 80 fname, text * 1
c
      write(*,*) ' '
      write(*,*) '>> Compute applied displacements matching K & T'
      write(*,*) ' '
      write(*,*) ' > Coordinate file name: '
      read(*,fmt='(a80)') fname
      infile = 10
      open(unit=infile,file=fname,status='old')
      write(*,*) 
     &   ' > Final K(I), T/sig-o, sig-o, number of load steps: '
      read(*,*) final_k, t_ratio, sig_o, num_steps
      applied_k = final_k / real(num_steps)
      applied_t = sig_o * t_ratio / real(num_steps)
      write(*,*) ' > Start radius, new scaled radius: '
      read(*,*) start_radius, new_radius
c
      read(infile,*) nnode
      do node = 1, nnode
       read(infile,*) ii, x, y, z
       x = x * (new_radius/start_radius)
       y = y * (new_radius/start_radius)
       coord(1,node) = x
       coord(2,node) = y
       coord(3,node) = z
      end do
      write(*,*) ' >> data file read'
c
c             loop to find all nodes on the outer most surface, i.e.,
c             those with radius.
c                                                     
      do node = 1, 5000
       node_list(node) = 0                 
      end do    
c
      tolerance = 0.01
      node_count = 0
      do node = 1, nnode
       x = coord(1,node)
       y = coord(2,node)
       r = abs( sqrt( x*x + y*y ) - new_radius )
       if ( r .le. tolerance ) then
          node_count = node_count + 1
          node_list(node) = 1
       end if
      end do
      write(*,*) '>> Number of nodes on boundary: ',
     &          node_count
c
c             nodes to process have been marked in list.
c
      rnu = 0.3
      e = 30000.0
      pi = 3.14159
      open(unit=10,file='constraints',status='unknown')
      open(unit=11,file='new_coords',status='unknown')
      outfile = 10
      do node = 1, nnode  
       x = coord(1,node) 
       y = coord(2,node) 
       if ( node_list(node) .eq. 1 ) then
            x = coord(1,node) 
            y = coord(2,node) 
            r = sqrt(x*x+y*y)
            theta = atan2( y, x ) 
            ct2 = cos(theta*0.5)
            st2 = sin(theta*0.5)
            ct  = cos(theta)
            st  = sin(theta)
            term1 = (3.0-4.0*rnu-ct) * ct2
            term2 = (3.0-4.0*rnu-ct) * st2
            term3 = sqrt( r * 0.5 / pi )
            term4 = (1.0 + rnu ) * applied_k / e
            u_k = term4 * term3 * term1
            v_k = term4 * term3 * term2
            u_t = (1.0-rnu*rnu) * applied_t * ct * r / e
            v_t = -rnu *(1.0+rnu) * applied_t * st * r / e 
            u_total = u_k + u_t
            v_total = v_k + v_t
            write(outfile,1000) node, u_total, v_total
        end if
       write(11,1100) node, x, y, coord(3,node)
      end do
c
      close(unit=10)
      close(unit=11)
      write(*,*) '.. new coordinates is file: coords'
      write(*,*) '.. K and T constraints in file: constraints'
c
 1000 format(4x,i5,' u  ',e16.8,' v ',e16.8)
 1100 format(2x,i5,3e16.8)
      stop
      end











