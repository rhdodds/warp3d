c
c----67--1---------2---------3---------4---------5---------6---------712
c
	subroutine reduce_coord_under(nod,n1,n2,n3,
     &                 j1,j2,jst,k1,k2,kst,left1)
C--- The routine creates the nodes in the domain under the crack
C--- front for i =< 5 and 1 =< j =< jma-4.
	implicit none
c
        include 'common_nod.f'
c
	integer  n1,n2,n3,nod(n1,n2,n3),j1,j2,jst,k1,k2,kst,j,k
	logical  left1,left
c
	do k=k1, k2, kst
	   left=left1
	   do j=j1, j2, jst
	      if (left) then
C . . . . . . Mid- and surface-nodes
	         npos(nod(3,j+2,k),1) = ( npos(nod(1,j+2,k),1) +
     &                 npos(nod(5,j+2,k),1) ) / 2.d0
	         npos(nod(3,j+2,k),2) = ( npos(nod(1,j+2,k),2) +
     &                 npos(nod(5,j+2,k),2) ) / 2.d0
	         npos(nod(3,j+2,k),3) = ( npos(nod(1,j+2,k),3) +
     &                 npos(nod(5,j+2,k),3) ) / 2.d0
 
	         npos(nod(2,j+1,k),1) = ( npos(nod(1,j+2,k),1) +
     &                 npos(nod(3,j,k),1) ) / 2.d0
	         npos(nod(2,j+1,k),2) = ( npos(nod(1,j+2,k),2) +
     &                 npos(nod(3,j,k),2) ) / 2.d0
	         npos(nod(2,j+1,k),3) = ( npos(nod(1,j+2,k),3) +
     &                 npos(nod(3,j,k),3) ) / 2.d0
 
	         npos(nod(1,j,k),1) = ( npos(nod(1,j+2,k),1) +
     &                 npos(nod(1,j-2,k),1) ) / 2.d0
	         npos(nod(1,j,k),2) = ( npos(nod(1,j+2,k),2) +
     &                 npos(nod(1,j-2,k),2) ) / 2.d0
	         npos(nod(1,j,k),3) = ( npos(nod(1,j+2,k),3) +
     &                 npos(nod(1,j-2,k),3) ) / 2.d0
C . . . . . . Surface- and centroid-nodes
	         npos(nod(3,j+1,k),1) = 0.25d0 *
     &               ( npos(nod(3,j+2,k),1) + npos(nod(4,j,k),1) +
     &                 npos(nod(2,j+1,k),1) + npos(nod(5,j+1,k),1) )
	         npos(nod(3,j+1,k),2) = 0.25d0 *
     &               ( npos(nod(3,j+2,k),2) + npos(nod(4,j,k),2) +
     &                 npos(nod(2,j+1,k),2) + npos(nod(5,j+1,k),2) )
	         npos(nod(3,j+1,k),3) = 0.25d0 *
     &               ( npos(nod(3,j+2,k),3) + npos(nod(4,j,k),3) +
     &                 npos(nod(2,j+1,k),3) + npos(nod(5,j+1,k),3) )
 
	         npos(nod(2,j,k),1) = 0.25d0 *
     &               ( npos(nod(2,j+1,k),1) + npos(nod(2,j-2,k),1) +
     &                 npos(nod(1,j,k),1) + npos(nod(3,j-1,k),1) )
	         npos(nod(2,j,k),2) = 0.25d0 *
     &               ( npos(nod(2,j+1,k),2) + npos(nod(2,j-2,k),2) +
     &                 npos(nod(1,j,k),2) + npos(nod(3,j-1,k),2) )
	         npos(nod(2,j,k),3) = 0.25d0 *
     &               ( npos(nod(2,j+1,k),3) + npos(nod(2,j-2,k),3) +
     &                 npos(nod(1,j,k),3) + npos(nod(3,j-1,k),3) )
	         left=.false.
	      else
C . . . . . . Mid- and surface-nodes
	         npos(nod(3,j-2,k),1) = ( npos(nod(1,j-2,k),1) +
     &                 npos(nod(5,j-2,k),1) ) / 2.d0
	         npos(nod(3,j-2,k),2) = ( npos(nod(1,j-2,k),2) +
     &                 npos(nod(5,j-2,k),2) ) / 2.d0
	         npos(nod(3,j-2,k),3) = ( npos(nod(1,j-2,k),3) +
     &                 npos(nod(5,j-2,k),3) ) / 2.d0
 
	         npos(nod(2,j-1,k),1) = ( npos(nod(1,j-2,k),1) +
     &                 npos(nod(3,j,k),1) ) / 2.d0
	         npos(nod(2,j-1,k),2) = ( npos(nod(1,j-2,k),2) +
     &                 npos(nod(3,j,k),2) ) / 2.d0
	         npos(nod(2,j-1,k),3) = ( npos(nod(1,j-2,k),3) +
     &                 npos(nod(3,j,k),3) ) / 2.d0
 
	         npos(nod(1,j,k),1) = ( npos(nod(1,j-2,k),1) +
     &                 npos(nod(1,j+2,k),1) ) / 2.d0
	         npos(nod(1,j,k),2) = ( npos(nod(1,j-2,k),2) +
     &                 npos(nod(1,j+2,k),2) ) / 2.d0
	         npos(nod(1,j,k),3) = ( npos(nod(1,j-2,k),3) +
     &                 npos(nod(1,j+2,k),3) ) / 2.d0
C . . . . . . Surface- and centroid-nodes
	         npos(nod(3,j-1,k),1) = 0.25d0 *
     &               ( npos(nod(3,j-2,k),1) + npos(nod(4,j,k),1) +
     &                 npos(nod(2,j-1,k),1) + npos(nod(5,j-1,k),1) )
	         npos(nod(3,j-1,k),2) = 0.25d0 *
     &               ( npos(nod(3,j-2,k),2) + npos(nod(4,j,k),2) +
     &                 npos(nod(2,j-1,k),2) + npos(nod(5,j-1,k),2) )
	         npos(nod(3,j-1,k),3) = 0.25d0 *
     &               ( npos(nod(3,j-2,k),3) + npos(nod(4,j,k),3) +
     &                 npos(nod(2,j-1,k),3) + npos(nod(5,j-1,k),3) )
 
	         npos(nod(2,j,k),1) = 0.25d0 *
     &               ( npos(nod(2,j-1,k),1) + npos(nod(2,j+2,k),1) +
     &                 npos(nod(1,j,k),1) + npos(nod(3,j+1,k),1) )
	         npos(nod(2,j,k),2) = 0.25d0 *
     &               ( npos(nod(2,j-1,k),2) + npos(nod(2,j+2,k),2) +
     &                 npos(nod(1,j,k),2) + npos(nod(3,j+1,k),2) )
	         npos(nod(2,j,k),3) = 0.25d0 *
     &               ( npos(nod(2,j-1,k),3) + npos(nod(2,j+2,k),3) +
     &                 npos(nod(1,j,k),3) + npos(nod(3,j+1,k),3) )
	         left=.true.
	      endif
	   enddo
	enddo
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine reduce_coord_2_to_1(nod,n1,n2,n3,ilocal,i1,i2,ist,
     &                        j1,j2,jst,k,left)
C--- The routine creates nodes in reduced mesh regions.
      implicit none
c
      include 'common_nod.f'
c
      integer  n1,n2,n3,nod(n1,n2,n3),i1,i2,ist,j1,j2,jst,k,ilocal,i,j
      logical  left
c
	if (ilocal.eq.1) then
	   do i=i1, i2, ist
	      if (left) then
	         do j=j1, j2, jst
c . . . . . . Mid- and surface-nodes
	            npos(nod(i+1,j,k+1),1) = (npos(nod(i,j,k),1) +
     &                     npos(nod(i+2,j,k+2),1) ) / 2.d0
	            npos(nod(i+1,j,k+1),2) = (npos(nod(i,j,k),2) +
     &                     npos(nod(i+2,j,k+2),2) ) / 2.d0
	            npos(nod(i+1,j,k+1),3) = (npos(nod(i,j,k),3) +
     &                     npos(nod(i+2,j,k+2),3) ) / 2.d0
c 
	            npos(nod(i+2,j,k),1) = (npos(nod(i+2,j,k-2),1) +
     &                     npos(nod(i+2,j,k+2),1) ) / 2.d0
	            npos(nod(i+2,j,k),2) = (npos(nod(i+2,j,k-2),2) +
     &                     npos(nod(i+2,j,k+2),2) ) / 2.d0
	            npos(nod(i+2,j,k),3) = (npos(nod(i+2,j,k-2),3) +
     &                     npos(nod(i+2,j,k+2),3) ) / 2.d0
c
	            npos(nod(i,j,k+2),1) = 0.5 *
     &                 (npos(nod(i-2,j,k+2),1)+npos(nod(i+2,j,k+2),1))
	            npos(nod(i,j,k+2),2) = 0.5 *
     &                 (npos(nod(i-2,j,k+2),2)+npos(nod(i+2,j,k+2),2))
	            npos(nod(i,j,k+2),3) = 0.5 *
     &                 (npos(nod(i-2,j,k+2),3)+npos(nod(i+2,j,k+2),3))
C . . . . . . Surface- and centroid-nodes
	            npos(nod(i,j,k+1),1) = 0.25d0 *
     &               (npos(nod(i-1,j,k),1) + npos(nod(i,j,k+2),1)+
     &                npos(nod(i-2,j,k+1),1) + npos(nod(i+1,j,k+1),1))
	            npos(nod(i,j,k+1),2) = 0.25d0 *
     &               (npos(nod(i-1,j,k),2) + npos(nod(i,j,k+2),2)+
     &                npos(nod(i-2,j,k+1),2) + npos(nod(i+1,j,k+1),2))
	            npos(nod(i,j,k+1),3) = 0.25d0 *
     &               (npos(nod(i-1,j,k),3) + npos(nod(i,j,k+2),3)+
     &                npos(nod(i-2,j,k+1),3) + npos(nod(i+1,j,k+1),3))
 
	            npos(nod(i+1,j,k),1) = 0.25d0 *
     &               (npos(nod(i+1,j,k-2),1) + npos(nod(i+1,j,k+1),1)+
     &                npos(nod(i,j,k-1),1) + npos(nod(i+2,j,k),1))
	            npos(nod(i+1,j,k),2) = 0.25d0 *
     &               (npos(nod(i+1,j,k-2),2) + npos(nod(i+1,j,k+1),2)+
     &                npos(nod(i,j,k-1),2) + npos(nod(i+2,j,k),2))
	            npos(nod(i+1,j,k),3) = 0.25d0 *
     &               (npos(nod(i+1,j,k-2),3) + npos(nod(i+1,j,k+1),3)+
     &                npos(nod(i,j,k-1),3) + npos(nod(i+2,j,k),3))
	         enddo
	         left=.false.
	      else
	         do j=j1, j2, jst
C . . . . . . Mid- and surface-nodes
	            npos(nod(i-1,j,k+1),1) = (npos(nod(i,j,k),1) +
     &                     npos(nod(i-2,j,k+2),1) ) / 2.d0
	            npos(nod(i-1,j,k+1),2) = (npos(nod(i,j,k),2) +
     &                     npos(nod(i-2,j,k+2),2) ) / 2.d0
	            npos(nod(i-1,j,k+1),3) = (npos(nod(i,j,k),3) +
     &                     npos(nod(i-2,j,k+2),3) ) / 2.d0
c
	            npos(nod(i-2,j,k),1) = (npos(nod(i-2,j,k-2),1) +
     &                     npos(nod(i-2,j,k+2),1) ) / 2.d0
	            npos(nod(i-2,j,k),2) = (npos(nod(i-2,j,k-2),2) +
     &                     npos(nod(i-2,j,k+2),2) ) / 2.d0
	            npos(nod(i-2,j,k),3) = (npos(nod(i-2,j,k-2),3) +
     &                     npos(nod(i-2,j,k+2),3) ) / 2.d0
c
	            npos(nod(i,j,k+2),1) = 0.5 *
     &                 (npos(nod(i-2,j,k+2),1)+npos(nod(i+2,j,k+2),1))
	            npos(nod(i,j,k+2),2) = 0.5 *
     &                 (npos(nod(i-2,j,k+2),2)+npos(nod(i+2,j,k+2),2))
	            npos(nod(i,j,k+2),3) = 0.5 *
     &                 (npos(nod(i-2,j,k+2),3)+npos(nod(i+2,j,k+2),3))
c
c . . . . . . Surface- and centroid-nodes
	            npos(nod(i,j,k+1),1) = 0.25d0 *
     &               (npos(nod(i+1,j,k),1) + npos(nod(i,j,k+2),1)+
     &                npos(nod(i+2,j,k+1),1) + npos(nod(i-1,j,k+1),1))
	            npos(nod(i,j,k+1),2) = 0.25d0 *
     &               (npos(nod(i+1,j,k),2) + npos(nod(i,j,k+2),2)+
     &                npos(nod(i+2,j,k+1),2) + npos(nod(i-1,j,k+1),2))
	            npos(nod(i,j,k+1),3) = 0.25d0 *
     &               (npos(nod(i+1,j,k),3) + npos(nod(i,j,k+2),3)+
     &                npos(nod(i+2,j,k+1),3) + npos(nod(i-1,j,k+1),3))
c
	            npos(nod(i-1,j,k),1) = 0.25d0 *
     &               (npos(nod(i-1,j,k-2),1) + npos(nod(i-1,j,k+1),1)+
     &                npos(nod(i,j,k-1),1) + npos(nod(i-2,j,k),1))
	            npos(nod(i-1,j,k),2) = 0.25d0 *
     &               (npos(nod(i-1,j,k-2),2) + npos(nod(i-1,j,k+1),2)+
     &                npos(nod(i,j,k-1),2) + npos(nod(i-2,j,k),2))
	            npos(nod(i-1,j,k),3) = 0.25d0 *
     &               (npos(nod(i-1,j,k-2),3) + npos(nod(i-1,j,k+1),3)+
     &                npos(nod(i,j,k-1),3) + npos(nod(i-2,j,k),3))
	         enddo
	         left=.true.
	      endif
	   enddo
	else
	   do j=j1, j2, jst
	      if (left) then
	         do i=i1, i2, ist
C . . . . . . Mid- and surface-nodes
 	            npos(nod(i,j+1,k+1),1) = (npos(nod(i,j,k),1) +
     &                     npos(nod(i,j+2,k+2),1) ) / 2.d0
	            npos(nod(i,j+1,k+1),2) = (npos(nod(i,j,k),2) +
     &                     npos(nod(i,j+2,k+2),2) ) / 2.d0
	            npos(nod(i,j+1,k+1),3) = (npos(nod(i,j,k),3) +
     &                     npos(nod(i,j+2,k+2),3) ) / 2.d0
c
	            npos(nod(i,j+2,k),1) = (npos(nod(i,j+2,k-2),1) +
     &                     npos(nod(i,j+2,k+2),1) ) / 2.d0
	            npos(nod(i,j+2,k),2) = (npos(nod(i,j+2,k-2),2) +
     &                     npos(nod(i,j+2,k+2),2) ) / 2.d0
	            npos(nod(i,j+2,k),3) = (npos(nod(i,j+2,k-2),3) +
     &                     npos(nod(i,j+2,k+2),3) ) / 2.d0
c
                    npos(nod(i,j,k+2),1) = 0.5 *
     &                 (npos(nod(i,j-2,k+2),1)+npos(nod(i,j+2,k+2),1))
                    npos(nod(i,j,k+2),2) = 0.5 *
     &                 (npos(nod(i,j-2,k+2),2)+npos(nod(i,j+2,k+2),2))
                    npos(nod(i,j,k+2),3) = 0.5 *
     &                 (npos(nod(i,j-2,k+2),3)+npos(nod(i,j+2,k+2),3))
C . . . . . . Surface- and centroid-nodes
	            npos(nod(i,j,k+1),1) = 0.25d0 *
     &               (npos(nod(i,j-1,k),1) + npos(nod(i,j,k+2),1)+
     &                npos(nod(i,j-2,k+1),1) + npos(nod(i,j+1,k+1),1))
	            npos(nod(i,j,k+1),2) = 0.25d0 *
     &               (npos(nod(i,j-1,k),2) + npos(nod(i,j,k+2),2)+
     &                npos(nod(i,j-2,k+1),2) + npos(nod(i,j+1,k+1),2))
	            npos(nod(i,j,k+1),3) = 0.25d0 *
     &               (npos(nod(i,j-1,k),3) + npos(nod(i,j,k+2),3)+
     &                npos(nod(i,j-2,k+1),3) + npos(nod(i,j+1,k+1),3))
c
	            npos(nod(i,j+1,k),1) = 0.25d0 *
     &               (npos(nod(i,j+1,k-2),1) + npos(nod(i,j+1,k+1),1)+
     &                npos(nod(i,j,k-1),1) + npos(nod(i,j+2,k),1))
	            npos(nod(i,j+1,k),2) = 0.25d0 *
     &               (npos(nod(i,j+1,k-2),2) + npos(nod(i,j+1,k+1),2)+
     &                npos(nod(i,j,k-1),2) + npos(nod(i,j+2,k),2))
	            npos(nod(i,j+1,k),3) = 0.25d0 *
     &               (npos(nod(i,j+1,k-2),3) + npos(nod(i,j+1,k+1),3)+
     &                npos(nod(i,j,k-1),3) + npos(nod(i,j+2,k),3))
	         enddo
	         left=.false.
	      else
	         do i=i1, i2, ist
C . . . . . . Mid- and surface-nodes
	            npos(nod(i,j-1,k+1),1) = (npos(nod(i,j,k),1) +
     &                     npos(nod(i,j-2,k+2),1) ) / 2.d0
	            npos(nod(i,j-1,k+1),2) = (npos(nod(i,j,k),2) +
     &                     npos(nod(i,j-2,k+2),2) ) / 2.d0
	            npos(nod(i,j-1,k+1),3) = (npos(nod(i,j,k),3) +
     &                     npos(nod(i,j-2,k+2),3) ) / 2.d0
c
	            npos(nod(i,j-2,k),1) = (npos(nod(i,j-2,k-2),1) +
     &                     npos(nod(i,j-2,k+2),1) ) / 2.d0
	            npos(nod(i,j-2,k),2) = (npos(nod(i,j-2,k-2),2) +
     &                     npos(nod(i,j-2,k+2),2) ) / 2.d0
	            npos(nod(i,j-2,k),3) = (npos(nod(i,j-2,k-2),3) +
     &                     npos(nod(i,j-2,k+2),3) ) / 2.d0
c
                    npos(nod(i,j,k+2),1) = 0.5 *
     &                 (npos(nod(i,j-2,k+2),1)+npos(nod(i,j+2,k+2),1))
                    npos(nod(i,j,k+2),2) = 0.5 *
     &                 (npos(nod(i,j-2,k+2),2)+npos(nod(i,j+2,k+2),2))
                    npos(nod(i,j,k+2),3) = 0.5 *
     &                 (npos(nod(i,j-2,k+2),3)+npos(nod(i,j+2,k+2),3))
c
c . . . . . . Surface- and centroid-nodes
	            npos(nod(i,j,k+1),1) = 0.25d0 *
     &               (npos(nod(i,j+1,k),1) + npos(nod(i,j,k+2),1)+
     &                npos(nod(i,j+2,k+1),1) + npos(nod(i,j-1,k+1),1))
	            npos(nod(i,j,k+1),2) = 0.25d0 *
     &               (npos(nod(i,j+1,k),2) + npos(nod(i,j,k+2),2)+
     &                npos(nod(i,j+2,k+1),2) + npos(nod(i,j-1,k+1),2))
	            npos(nod(i,j,k+1),3) = 0.25d0 *
     &               (npos(nod(i,j+1,k),3) + npos(nod(i,j,k+2),3)+
     &                npos(nod(i,j+2,k+1),3) + npos(nod(i,j-1,k+1),3))
c
	            npos(nod(i,j-1,k),1) = 0.25d0 *
     &               (npos(nod(i,j-1,k-2),1) + npos(nod(i,j-1,k+1),1)+
     &                npos(nod(i,j,k-1),1) + npos(nod(i,j-2,k),1))
	            npos(nod(i,j-1,k),2) = 0.25d0 *
     &               (npos(nod(i,j-1,k-2),2) + npos(nod(i,j-1,k+1),2)+
     &                npos(nod(i,j,k-1),2) + npos(nod(i,j-2,k),2))
	            npos(nod(i,j-1,k),3) = 0.25d0 *
     &               (npos(nod(i,j-1,k-2),3) + npos(nod(i,j-1,k+1),3)+
     &                npos(nod(i,j,k-1),3) + npos(nod(i,j-2,k),3))
	         enddo
	         left=.true.
	      endif
	   enddo
	endif
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine reduce_coord_plane_3_to_1(nod,n1,n2,n3,ilocal,
     &                  i1,i2,ist,j1,j2,jst,k)
C--- The routine creates nodes in reduced mesh regions.
      implicit none
c
      include 'common_nod.f'
c
      integer  n1,n2,n3,nod(n1,n2,n3),i1,i2,ist,j1,j2,jst,k,ilocal,i,j
c
      if (ilocal.eq.1) then
         do i=i1, i2, ist
            do j=j1, j2, jst
c . . . . . Mid- and surface-nodes
c pt.1
               npos(nod(i-3,j,k),1) = 0.5*( npos(nod(i-3,j,k-2),1) +
     &                                      npos(nod(i-3,j,k+2),1) )
               npos(nod(i-3,j,k),2) = 0.5*( npos(nod(i-3,j,k-2),2) +
     &                                      npos(nod(i-3,j,k+2),2) )
               npos(nod(i-3,j,k),3) = 0.5*( npos(nod(i-3,j,k-2),3) +
     &                                      npos(nod(i-3,j,k+2),3) )
c pt.2
               npos(nod(i-2,j,k+1),1) = 0.5*( npos(nod(i-1,j,k),1) +
     &                                        npos(nod(i-3,j,k+2),1) )
               npos(nod(i-2,j,k+1),2) = 0.5*( npos(nod(i-1,j,k),2) +
     &                                        npos(nod(i-3,j,k+2),2) )
               npos(nod(i-2,j,k+1),3) = 0.5*( npos(nod(i-1,j,k),3) +
     &                                        npos(nod(i-3,j,k+2),3) )
c pt.3
               npos(nod(i,j,k+2),1) = 0.5*( npos(nod(i-3,j,k+2),1) +
     &                                      npos(nod(i+3,j,k+2),1) )
               npos(nod(i,j,k+2),2) = 0.5*( npos(nod(i-3,j,k+2),2) +
     &                                      npos(nod(i+3,j,k+2),2) )
               npos(nod(i,j,k+2),3) = 0.5*( npos(nod(i-3,j,k+2),3) +
     &                                      npos(nod(i+3,j,k+2),3) )
c pt.4
               npos(nod(i+2,j,k+1),1) = 0.5*( npos(nod(i+1,j,k),1) +
     &                                        npos(nod(i+3,j,k+2),1) )
               npos(nod(i+2,j,k+1),2) = 0.5*( npos(nod(i+1,j,k),2) +
     &                                        npos(nod(i+3,j,k+2),2) )
               npos(nod(i+2,j,k+1),3) = 0.5*( npos(nod(i+1,j,k),3) +
     &                                        npos(nod(i+3,j,k+2),3) )
c pt.5
               npos(nod(i+3,j,k),1) = 0.5*( npos(nod(i+3,j,k-2),1) +
     &                                      npos(nod(i+3,j,k+2),1) )
               npos(nod(i+3,j,k),2) = 0.5*( npos(nod(i+3,j,k-2),2) +
     &                                      npos(nod(i+3,j,k+2),2) )
               npos(nod(i+3,j,k),3) = 0.5*( npos(nod(i+3,j,k-2),3) +
     &                                      npos(nod(i+3,j,k+2),3) )
c . . . . . Surface- and centroid-nodes
c pt.6
               npos(nod(i-2,j,k),1) =  0.25 *
     &            ( npos(nod(i-2,j,k-2),1) + npos(nod(i-2,j,k+1),1) +
     &              npos(nod(i-3,j,k),1) + npos(nod(i-1,j,k),1) )
               npos(nod(i-2,j,k),2) =  0.25 *
     &            ( npos(nod(i-2,j,k-2),2) + npos(nod(i-2,j,k+1),2) +
     &              npos(nod(i-3,j,k),2) + npos(nod(i-1,j,k),2) )
               npos(nod(i-2,j,k),3) =  0.25 *
     &            ( npos(nod(i-2,j,k-2),3) + npos(nod(i-2,j,k+1),3) +
     &              npos(nod(i-3,j,k),3) + npos(nod(i-1,j,k),3) )
c pt.7
               npos(nod(i+2,j,k),1) =  0.25 *
     &            ( npos(nod(i+2,j,k-2),1) + npos(nod(i+2,j,k+1),1) +
     &              npos(nod(i+1,j,k),1) + npos(nod(i+3,j,k),1) )
               npos(nod(i+2,j,k),2) =  0.25 *
     &            ( npos(nod(i+2,j,k-2),2) + npos(nod(i+2,j,k+1),2) +
     &              npos(nod(i+1,j,k),2) + npos(nod(i+3,j,k),2) )
               npos(nod(i+2,j,k),3) =  0.25 *
     &            ( npos(nod(i+2,j,k-2),3) + npos(nod(i+2,j,k+1),3) +
     &              npos(nod(i+1,j,k),3) + npos(nod(i+3,j,k),3) )
            enddo
         enddo
      else
         do j=j1, j2, jst
            do i=i1, i2, ist
c . . . . . Mid- and surface-nodes
c pt.1
               npos(nod(i,j-3,k),1) = 0.5*( npos(nod(i,j-3,k-2),1) +
     &                                      npos(nod(i,j-3,k+2),1) )
               npos(nod(i,j-3,k),2) = 0.5*( npos(nod(i,j-3,k-2),2) +
     &                                      npos(nod(i,j-3,k+2),2) )
               npos(nod(i,j-3,k),3) = 0.5*( npos(nod(i,j-3,k-2),3) +
     &                                      npos(nod(i,j-3,k+2),3) )
c pt.2
               npos(nod(i,j-2,k+1),1) = 0.5*( npos(nod(i,j-1,k),1) +
     &                                        npos(nod(i,j-3,k+2),1) )
               npos(nod(i,j-2,k+1),2) = 0.5*( npos(nod(i,j-1,k),2) +
     &                                        npos(nod(i,j-3,k+2),2) )
               npos(nod(i,j-2,k+1),3) = 0.5*( npos(nod(i,j-1,k),3) +
     &                                        npos(nod(i,j-3,k+2),3) )
c pt.3
               npos(nod(i,j,k+2),1) = 0.5*( npos(nod(i,j-3,k+2),1) +
     &                                      npos(nod(i,j+3,k+2),1) )
               npos(nod(i,j,k+2),2) = 0.5*( npos(nod(i,j-3,k+2),2) +
     &                                      npos(nod(i,j+3,k+2),2) )
               npos(nod(i,j,k+2),3) = 0.5*( npos(nod(i,j-3,k+2),3) +
     &                                      npos(nod(i,j+3,k+2),3) )
c pt.4
               npos(nod(i,j+2,k+1),1) = 0.5*( npos(nod(i,j+1,k),1) +
     &                                        npos(nod(i,j+3,k+2),1) )
               npos(nod(i,j+2,k+1),2) = 0.5*( npos(nod(i,j+1,k),2) +
     &                                        npos(nod(i,j+3,k+2),2) )
               npos(nod(i,j+2,k+1),3) = 0.5*( npos(nod(i,j+1,k),3) +
     &                                        npos(nod(i,j+3,k+2),3) )
c pt.5
               npos(nod(i,j+3,k),1) = 0.5*( npos(nod(i,j+3,k-2),1) +
     &                                      npos(nod(i,j+3,k+2),1) )
               npos(nod(i,j+3,k),2) = 0.5*( npos(nod(i,j+3,k-2),2) +
     &                                      npos(nod(i,j+3,k+2),2) )
               npos(nod(i,j+3,k),3) = 0.5*( npos(nod(i,j+3,k-2),3) +
     &                                      npos(nod(i,j+3,k+2),3) )
c . . . . . Surface- and centroid-nodes
c pt.6
               npos(nod(i,j-2,k),1) =  0.25 *
     &            ( npos(nod(i,j-2,k-2),1) + npos(nod(i,j-2,k+1),1) +
     &              npos(nod(i,j-3,k),1) + npos(nod(i,j-1,k),1) )
               npos(nod(i,j-2,k),2) =  0.25 *
     &            ( npos(nod(i,j-2,k-2),2) + npos(nod(i,j-2,k+1),2) +
     &              npos(nod(i,j-3,k),2) + npos(nod(i,j-1,k),2) )
               npos(nod(i,j-2,k),3) =  0.25 *
     &            ( npos(nod(i,j-2,k-2),3) + npos(nod(i,j-2,k+1),3) +
     &              npos(nod(i,j-3,k),3) + npos(nod(i,j-1,k),3) )
c pt.7
               npos(nod(i,j+2,k),1) =  0.25 *
     &            ( npos(nod(i,j+2,k-2),1) + npos(nod(i,j+2,k+1),1) +
     &              npos(nod(i,j+1,k),1) + npos(nod(i,j+3,k),1) )
               npos(nod(i,j+2,k),2) =  0.25 *
     &            ( npos(nod(i,j+2,k-2),2) + npos(nod(i,j+2,k+1),2) +
     &              npos(nod(i,j+1,k),2) + npos(nod(i,j+3,k),2) )
               npos(nod(i,j+2,k),3) =  0.25 *
     &            ( npos(nod(i,j+2,k-2),3) + npos(nod(i,j+2,k+1),3) +
     &              npos(nod(i,j+1,k),3) + npos(nod(i,j+3,k),3) )
            enddo
         enddo
      endif
      return
      end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine reduce_coord_3_to_1(nod,n1,n2,n3,j1,j2,jst,k,left)
C--- The routine creates nodes in reduced mesh region under the crack
C--- front.
      implicit none
c
      include 'common_nod.f'
c
      integer  n1,n2,n3,nod(n1,n2,n3),j1,j2,jst,k,j
      logical  left
c
	do j=j1, j2, jst
	   if (left) then
C . . . . . Mid nodes local element no. 1
	      npos(nod(2,j+1,k+1),1) = ( npos(nod(1,j+2,k+2),1) +
     &             npos(nod(3,j,k),1) ) / 2.d0
	      npos(nod(2,j+1,k+1),2) = ( npos(nod(1,j+2,k+2),2) +
     &             npos(nod(3,j,k),2) ) / 2.d0
	      npos(nod(2,j+1,k+1),3) = ( npos(nod(1,j+2,k+2),3) +
     &             npos(nod(3,j,k),3) ) / 2.d0
 
	      npos(nod(2,j-2,k+1),1) = ( npos(nod(1,j-2,k+2),1) +
     &             npos(nod(3,j-2,k),1) ) / 2.d0
	      npos(nod(2,j-2,k+1),2) = ( npos(nod(1,j-2,k+2),2) +
     &             npos(nod(3,j-2,k),2) ) / 2.d0
	      npos(nod(2,j-2,k+1),3) = ( npos(nod(1,j-2,k+2),3) +
     &             npos(nod(3,j-2,k),3) ) / 2.d0
C . . . . . Surface nodes local element no. 1
	      npos(nod(2,j+1,k),1) = 0.25d0 *
     &         ( npos(nod(1,j+2,k),1) + npos(nod(3,j,k-1),1) +
     &           npos(nod(2,j+1,k-2),1) + npos(nod(2,j+1,k+1),1) )
	      npos(nod(2,j+1,k),2) = 0.25d0 *
     &         ( npos(nod(1,j+2,k),2) + npos(nod(3,j,k-1),2) +
     &           npos(nod(2,j+1,k-2),2) + npos(nod(2,j+1,k+1),2) )
	      npos(nod(2,j+1,k),3) = 0.25d0 *
     &         ( npos(nod(1,j+2,k),3) + npos(nod(3,j,k-1),3) +
     &           npos(nod(2,j+1,k-2),3) + npos(nod(2,j+1,k+1),3) )
 
	      npos(nod(2,j,k+1),1) = 0.25d0 *
     &         ( npos(nod(1,j,k+2),1) + npos(nod(3,j-1,k),1) +
     &           npos(nod(2,j+1,k+1),1) + npos(nod(2,j-2,k+1),1) )
	      npos(nod(2,j,k+1),2) = 0.25d0 *
     &         ( npos(nod(1,j,k+2),2) + npos(nod(3,j-1,k),2) +
     &           npos(nod(2,j+1,k+1),2) + npos(nod(2,j-2,k+1),2) )
	      npos(nod(2,j,k+1),3) = 0.25d0 *
     &         ( npos(nod(1,j,k+2),3) + npos(nod(3,j-1,k),3) +
     &           npos(nod(2,j+1,k+1),3) + npos(nod(2,j-2,k+1),3) )
C . . . . . Centroid node local element no. 1
	      npos(nod(2,j,k),1) = 0.25d0 *
     &         ( npos(nod(2,j+1,k),1) + npos(nod(2,j-2,k),1) +
     &           npos(nod(2,j,k-2),1) + npos(nod(2,j,k+1),1) )
	      npos(nod(2,j,k),2) = 0.25d0 *
     &         ( npos(nod(2,j+1,k),2) + npos(nod(2,j-2,k),2) +
     &           npos(nod(2,j,k-2),2) + npos(nod(2,j,k+1),2) )
	      npos(nod(2,j,k),3) = 0.25d0 *
     &         ( npos(nod(2,j+1,k),3) + npos(nod(2,j-2,k),3) +
     &           npos(nod(2,j,k-2),3) + npos(nod(2,j,k+1),3) )
C . . . . . Surface node local element no. 2
	      npos(nod(3,j+1,k+1),1) = 0.25d0 *
     &         ( npos(nod(2,j+1,k+1),1) + npos(nod(5,j+1,k+1),1) +
     &           npos(nod(3,j+2,k+2),1) + npos(nod(4,j,k),1) )
	      npos(nod(3,j+1,k+1),2) = 0.25d0 *
     &         ( npos(nod(2,j+1,k+1),2) + npos(nod(5,j+1,k+1),2) +
     &           npos(nod(3,j+2,k+2),2) + npos(nod(4,j,k),2) )
	      npos(nod(3,j+1,k+1),3) = 0.25d0 *
     &         ( npos(nod(2,j+1,k+1),3) + npos(nod(5,j+1,k+1),3) +
     &           npos(nod(3,j+2,k+2),3) + npos(nod(4,j,k),3) )
C . . . . . Centroid node local element no. 2
	      npos(nod(3,j+1,k),1) = 0.25d0 *
     &         ( npos(nod(2,j+1,k),1) + npos(nod(5,j+1,k),1) +
     &           npos(nod(3,j+1,k-2),1) + npos(nod(3,j+1,k+1),1) )
	      npos(nod(3,j+1,k),2) = 0.25d0 *
     &         ( npos(nod(2,j+1,k),2) + npos(nod(5,j+1,k),2) +
     &           npos(nod(3,j+1,k-2),2) + npos(nod(3,j+1,k+1),2) )
	      npos(nod(3,j+1,k),3) = 0.25d0 *
     &         ( npos(nod(2,j+1,k),3) + npos(nod(5,j+1,k),3) +
     &           npos(nod(3,j+1,k-2),3) + npos(nod(3,j+1,k+1),3) )
C . . . . . Surface node local element no. 3
	      npos(nod(3,j-2,k+1),1) = 0.25d0 *
     &         ( npos(nod(2,j-2,k+1),1) + npos(nod(5,j-2,k+1),1) +
     &           npos(nod(4,j-2,k),1) + npos(nod(3,j-2,k+2),1) )
	      npos(nod(3,j-2,k+1),2) = 0.25d0 *
     &         ( npos(nod(2,j-2,k+1),2) + npos(nod(5,j-2,k+1),2) +
     &           npos(nod(4,j-2,k),2) + npos(nod(3,j-2,k+2),2) )
	      npos(nod(3,j-2,k+1),3) = 0.25d0 *
     &         ( npos(nod(2,j-2,k+1),3) + npos(nod(5,j-2,k+1),3) +
     &           npos(nod(4,j-2,k),3) + npos(nod(3,j-2,k+2),3) )
C . . . . . Centroid node local element no. 3
	      npos(nod(3,j,k+1),1) = 0.25d0 *
     &         ( npos(nod(2,j,k+1),1) + npos(nod(5,j,k+1),1) +
     &           npos(nod(3,j+1,k+1),1) + npos(nod(3,j-2,k+1),1) )
	      npos(nod(3,j,k+1),2) = 0.25d0 *
     &         ( npos(nod(2,j,k+1),2) + npos(nod(5,j,k+1),2) +
     &           npos(nod(3,j+1,k+1),2) + npos(nod(3,j-2,k+1),2) )
	      npos(nod(3,j,k+1),3) = 0.25d0 *
     &         ( npos(nod(2,j,k+1),3) + npos(nod(5,j,k+1),3) +
     &           npos(nod(3,j+1,k+1),3) + npos(nod(3,j-2,k+1),3) )
	      left=.false.
	   else
C . . . . . Mid nodes local element no. 1
	      npos(nod(2,j-1,k+1),1) = ( npos(nod(1,j-2,k+2),1) +
     &             npos(nod(3,j,k),1) ) / 2.d0
	      npos(nod(2,j-1,k+1),2) = ( npos(nod(1,j-2,k+2),2) +
     &             npos(nod(3,j,k),2) ) / 2.d0
	      npos(nod(2,j-1,k+1),3) = ( npos(nod(1,j-2,k+2),3) +
     &             npos(nod(3,j,k),3) ) / 2.d0
 
	      npos(nod(2,j+2,k+1),1) = ( npos(nod(1,j+2,k+2),1) +
     &             npos(nod(3,j+2,k),1) ) / 2.d0
	      npos(nod(2,j+2,k+1),2) = ( npos(nod(1,j+2,k+2),2) +
     &             npos(nod(3,j+2,k),2) ) / 2.d0
	      npos(nod(2,j+2,k+1),3) = ( npos(nod(1,j+2,k+2),3) +
     &             npos(nod(3,j+2,k),3) ) / 2.d0
C . . . . . Surface nodes local element no. 1
	      npos(nod(2,j-1,k),1) = 0.25d0 *
     &         ( npos(nod(1,j-2,k),1) + npos(nod(3,j,k-1),1) +
     &           npos(nod(2,j-1,k-2),1) + npos(nod(2,j-1,k+1),1) )
	      npos(nod(2,j-1,k),2) = 0.25d0 *
     &         ( npos(nod(1,j-2,k),2) + npos(nod(3,j,k-1),2) +
     &           npos(nod(2,j-1,k-2),2) + npos(nod(2,j-1,k+1),2) )
	      npos(nod(2,j-1,k),3) = 0.25d0 *
     &         ( npos(nod(1,j-2,k),3) + npos(nod(3,j,k-1),3) +
     &           npos(nod(2,j-1,k-2),3) + npos(nod(2,j-1,k+1),3) )
 
	      npos(nod(2,j,k+1),1) = 0.25d0 *
     &         ( npos(nod(1,j,k+2),1) + npos(nod(3,j+1,k),1) +
     &           npos(nod(2,j-1,k+1),1) + npos(nod(2,j+2,k+1),1) )
	      npos(nod(2,j,k+1),2) = 0.25d0 *
     &         ( npos(nod(1,j,k+2),2) + npos(nod(3,j+1,k),2) +
     &           npos(nod(2,j-1,k+1),2) + npos(nod(2,j+2,k+1),2) )
	      npos(nod(2,j,k+1),3) = 0.25d0 *
     &         ( npos(nod(1,j,k+2),3) + npos(nod(3,j+1,k),3) +
     &           npos(nod(2,j-1,k+1),3) + npos(nod(2,j+2,k+1),3) )
C . . . . . Centroid node local element no. 1
	      npos(nod(2,j,k),1) = 0.25d0 *
     &         ( npos(nod(2,j-1,k),1) + npos(nod(2,j+2,k),1) +
     &           npos(nod(2,j,k-2),1) + npos(nod(2,j,k+1),1) )
	      npos(nod(2,j,k),2) = 0.25d0 *
     &         ( npos(nod(2,j-1,k),2) + npos(nod(2,j+2,k),2) +
     &           npos(nod(2,j,k-2),2) + npos(nod(2,j,k+1),2) )
	      npos(nod(2,j,k),3) = 0.25d0 *
     &         ( npos(nod(2,j-1,k),3) + npos(nod(2,j+2,k),3) +
     &           npos(nod(2,j,k-2),3) + npos(nod(2,j,k+1),3) )
C . . . . . Surface node local element no. 2
	      npos(nod(3,j-1,k+1),1) = 0.25d0 *
     &         ( npos(nod(2,j-1,k+1),1) + npos(nod(5,j-1,k+1),1) +
     &           npos(nod(3,j-2,k+2),1) + npos(nod(4,j,k),1) )
	      npos(nod(3,j-1,k+1),2) = 0.25d0 *
     &         ( npos(nod(2,j-1,k+1),2) + npos(nod(5,j-1,k+1),2) +
     &           npos(nod(3,j-2,k+2),2) + npos(nod(4,j,k),2) )
	      npos(nod(3,j-1,k+1),3) = 0.25d0 *
     &         ( npos(nod(2,j-1,k+1),3) + npos(nod(5,j-1,k+1),3) +
     &           npos(nod(3,j-2,k+2),3) + npos(nod(4,j,k),3) )
C . . . . . Centroid node local element no. 2
	      npos(nod(3,j-1,k),1) = 0.25d0 *
     &         ( npos(nod(2,j-1,k),1) + npos(nod(5,j-1,k),1) +
     &           npos(nod(3,j-1,k-2),1) + npos(nod(3,j-1,k+1),1) )
	      npos(nod(3,j-1,k),2) = 0.25d0 *
     &         ( npos(nod(2,j-1,k),2) + npos(nod(5,j-1,k),2) +
     &           npos(nod(3,j-1,k-2),2) + npos(nod(3,j-1,k+1),2) )
	      npos(nod(3,j-1,k),3) = 0.25d0 *
     &         ( npos(nod(2,j-1,k),3) + npos(nod(5,j-1,k),3) +
     &           npos(nod(3,j-1,k-2),3) + npos(nod(3,j-1,k+1),3) )
C . . . . . Surface node local element no. 3
	      npos(nod(3,j+2,k+1),1) = 0.25d0 *
     &         ( npos(nod(2,j+2,k+1),1) + npos(nod(5,j+2,k+1),1) +
     &           npos(nod(4,j+2,k),1) + npos(nod(3,j+2,k+2),1) )
	      npos(nod(3,j+2,k+1),2) = 0.25d0 *
     &         ( npos(nod(2,j+2,k+1),2) + npos(nod(5,j+2,k+1),2) +
     &           npos(nod(4,j+2,k),2) + npos(nod(3,j+2,k+2),2) )
	      npos(nod(3,j+2,k+1),3) = 0.25d0 *
     &         ( npos(nod(2,j+2,k+1),3) + npos(nod(5,j+2,k+1),3) +
     &           npos(nod(4,j+2,k),3) + npos(nod(3,j+2,k+2),3) )
C . . . . . Centroid node local element no. 3
	      npos(nod(3,j,k+1),1) = 0.25d0 *
     &         ( npos(nod(2,j,k+1),1) + npos(nod(5,j,k+1),1) +
     &           npos(nod(3,j-1,k+1),1) + npos(nod(3,j+2,k+1),1) )
	      npos(nod(3,j,k+1),2) = 0.25d0 *
     &         ( npos(nod(2,j,k+1),2) + npos(nod(5,j,k+1),2) +
     &           npos(nod(3,j-1,k+1),2) + npos(nod(3,j+2,k+1),2) )
	      npos(nod(3,j,k+1),3) = 0.25d0 *
     &         ( npos(nod(2,j,k+1),3) + npos(nod(5,j,k+1),3) +
     &           npos(nod(3,j-1,k+1),3) + npos(nod(3,j+2,k+1),3) )
	      left=.true.
	   endif
	enddo
	return
	end
c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
