c
c----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine nodposition(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     1                       nodb,nb1,nb2,nb3,
     2   natyp,y,z,omega,job,jobh,pl,rs,etyp,zb_xw,rate_b,rzero)
C-----------------------------------------------------------------C
C   Rutinen skapar koordinater f|r noderna i alla zoner .         C
C                                                                 C
C-----------------------------------------------------------------C
	implicit none
c
        include 'mesh3d_scp_common_nod.f'
c
        integer  ns1,ns2,ns3,nods(ns1,ns2,ns3),
     1           na1,na2,na3,noda(na1,na2,na3),
     2           nb1,nb2,nb3,nodb(nb1,nb2,nb3)
c
        integer  mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     1           natyp,ims,jms,kms,is,ima,jma,imb,jmb,kma,
     2           etyp, kstart,kstep,kmsi,ia1,ia2,ia3,
     3           ja1,ja2,ka1,ka2,nodnum,pl,i,j,k,it,ids,ic,
     4           ksr1,kar1,kar2,imt,jobh,ilocal,zb_xw, rzero,  ppp
c
	real*8   t,w,c,a,kappa,alfa,r1,r2,rn,eta,my,y(0:200),
     1           z(0:200,2),omega,fi(50),r(50,50,50),fi0,
     2           dr(50),zik(50,50),lambda,beta,r2y,ro,ru,dz,dy,ny,dny,
     3           ai,rot,psi,xe,ze,rs(200),teta,red1,red2,alfayz,rate_b,
     4           drho(50),phi(50,50),rho(50),ztmp,ron,lh1,lh2,trc(50,5)
c
c                     Rami : Add/Mod real local vars
c
	real*8   fact1,fact2
 
C    REAL FUNCTIONS !
	real*8   xaxel,zaxel,rhox,rhoz,rho_ellips,hyperbel_length,
     &           calc_beta,calc_lambda
 
	character  job*40
	logical    left
 
C test variables:
C	INTEGER   UNIT
C	CHARACTER SVAR*1
 
 
	common /nod/ mr,mf,mv,mh,m1,m2,ma,na,mb,nb,lt,lred,rtype,sfred,
     1         /max/ ims,jms,kms,is,ima,jma,imb,jmb,kma,
     2         /geom/t,w,c,a,kappa,alfa,r1,r2,rn,eta,my,
     3         /nblock/ kstart,kstep,
     4         /reduce/ ksr1,kar1,kar2
 
C------- MATRISER OCH VEKTORER ----------------------------
C   INTEGER: NODS(i,j,k) , NODA(i,j,k) , NODB(i,j,k)
C   REAL   : R(i,k,j) i zonS & R(i,j,k) i zonA , ZIK(I,K)
C----------------------------------------------------------
	write(*,'(t5,a)') '* generating the nodes . . . .'
 
C.......1) Best{mm vinklarna FI :
	if (natyp.eq.1) then
C.......A) Enligt typ 1 :
	      do j=1, jms
	         fi(j)=90.d0*dble(jms-j)/dble(jms-1)
 	      enddo
	elseif (natyp.eq.2) then
C	B) Enligt typ 2 :
	   fi(1)=90.d0
	   fi(jms)=0.d0
	   do j=5, jms-2, 2
	      i=j-2
	      psi = ( dble(jms-2-i)/dble(jms-3) )**omega
	      fi(j)= dasind( psi )
	   enddo
	   fi(3)=(fi(1)+fi(5))/2.d0
	   fi(jms-2)=dasind( dsind(fi(jms-4))/2.d0 )
	   do j=2, jms, 2
	      fi(j)=(fi(j-1)+fi(j+1))/2.d0
	   enddo
	elseif (natyp.eq.3) then
C	C) Enligt typ 3 :
	   fi(1)=90.d0
	   fi(jms)=0.d0
	   fi0=90.d0*(1.d0-omega)/(1.d0-omega**dble(na) )
	   do j=1, na-1
	      fi(2*j+1)=fi(2*j-1)-fi0*omega**dble(j-1)
	   enddo
	   do j=2, jms, 2
	      fi(j)=(fi(j-1)+fi(j+1))/2.d0
	   enddo
	elseif (natyp.ge.4) then
C	C) Enligt typ 4 :
C.........:N]GON TYP >= 4 EXISTERAR EJ [NNU !!
           write(*,'(t10,a)')
     &         'note that natyp must be 1,2 or 3, fulfilled !'
 	   stop
	endif
 
C---------------------------------------------------
C.......2) Geometri i Zon S :
C---------------------------------------------------
C . . . a) Den inre "cirkeln" !
c 
C	kmsi=kms-2
C	lambda=log(my)
C	dr(1)=rn
C	if (sfred.eq.0) then
C	   beta=calc_beta(my,eta,mr-1)
C	   ids=1
C	else
C	   beta=calc_beta(my,eta,mr-2)
C	   ids=2
C	endif
C	do k=3, kmsi, 2
C	   dr(k)=dr(k-2)+((r1-rn)/eta)*
C     &           dexp(lambda*dble((k-3)/2)**beta)
C	   rs((k-1)/2)=dr(k)
C	enddo
C	if (sfred.gt.0) then
C	   do k=kmsi, ksr1+2, -2
C	      dr(k)=dr(k-2)
C	      rs((k-1)/2)=dr(k-2)
C	   enddo
C           dr(ksr1+2)=0.5d0*(dr(ksr1+4)+dr(ksr1-2))
C	   dr(ksr1)=dr(ksr1+2)*dsind( 45.d0*(1.d0-2.d0/dble(mf)) ) /
C     &             dsind( 45.d0*(1.d0+2.d0/dble(mf)) )
C	   rs((ksr1+1)/2)=dr(ksr1+2)
C	   rs((ksr1-1)/2)=dr(ksr1)
C	   red1=dr(ksr1)-dr(ksr1-2)
C	   red2=dr(ksr1+2)-dr(ksr1)
C	   if ( ((red1/red2).gt.3.d0).or.((red1/red2).lt.0.25d0) ) then
C	      write(*,'(t10,a)')
C     1        '=> the model does not work with sfred > 0 with the',
C     2        '   currant choice of parameters, remedy change mr!'
C	      stop
C	   endif
C	endif
c----------------------   Modified by Rami Haj-Ali 9/27/1996 ------------------
c
c      1) Set r for the first two rings to be equal to 0.12*R11 ( 0.06*R11 each )
c      2) ther rest (0.88*R11) is devided between MR-3:
c           ( 3 because of : first two rings (1,2) + last ring (MR) )
c           Where:  R11 =  R1 - Rn
c 
c      Note:  The element distances from the notch radius are stored in the
c             odd numbers of the array: dr(k), k=1,3,5,7,....
c
        kmsi=kms-2
	ids=1
	fact1=0.06
	fact2=1.0-fact1
c                 Notch Radius
	dr(1)=rn
c                 1st Ring
	dr(3)=dr(1) + fact1*(R1-rn)
	rs(1)=dr(3)
c
c                 2nd Ring
	dr(5)=dr(3) + fact1*(R1-rn)
	rs(2)=dr(5)
c
c       Reset Eta and My according to modification
c
	eta=R1/(dr(3)-dr(1))
	my =(dr(5)-dr(3))/(dr(3)-dr(1))
c
c       Rest of the Rings
c
	do k=7, kmsi, 2
	   dr(k)=dr(k-2)+(fact2/float(mr-3))*(r1-Rn)
	   rs((k-1)/2)=dr(k)
	enddo
c
	if(sfred.gt.0) then
	write(*,*) 
	write(*,*)  ' >>>>>>>   ERROR '
	write(*,*)  '  This is a modified version that does not include'
	write(*,*)  '  SFRED=0  --  USE ORIGENAL VERSION OF MESH3D_SC '
	stop
	end if
c
c-------------------   end of modifications
C...... Radier f|r fi=90 och I=IMS ber{knas.
 
	do k=1, kmsi, 2
	   drho(k)=rhoz(a-rn,90.d0,alfa)-rhoz(a-rn-dr(k),90.d0,alfa)
	enddo
C . . . Best{m vinklarna i Y-X planet.
	imt=(ims-1)/2+1
	do i=1, imt
	   phi(i,1)=180.d0*dble(i-1)/dble(ims-1)
	   phi(ims-i+1,1)=phi(i,1)
	   do k=1, kms, 2
	      phi(i,k)=phi(i,1)
	      phi(ims-i+1,k)=phi(i,1)
	   enddo
	enddo
C . . Radier framf|r sprickan ber{knas ( I =< IMT ).
	do i=1, imt
	   do k=1, kmsi,2
	      r(i,k,1)=rhoz(a-rn,90.d0,alfa)+drho(k)*dcosd(phi(i,k))
	   enddo
	enddo
C . . Radier bakom sprickan ber{knas ( I > IMT ).
C       H{r {r ZAXEL -> YAXEL & XAXEL -> ZAXEL
C               RHOZ -> RHOY  &  RHOX -> RHOZ
	alfayz=0.5d0*dsqrt( dr(5)**2.d0-(dr(5)-rn)**2.d0)
	do k=3, kmsi, 2
	   rho(k)=rhoz(dr(k)-rn,90.d0,alfayz)
	enddo
        ztmp=xaxel(rho(3),phi(ims-2,1),alfayz)/dble(mf/2-1)
	do i=imt+2, ims-4, 2
           phi(i,3)=dacosd(ztmp*0.5d0*dble(i-imt) /
     &             (alfayz*(rho(3)+1.d0/rho(3))))
	   phi(i,5)=0.67d0*phi(i,3)+0.33d0*phi(i,9)
	   phi(i,7)=0.33d0*phi(i,3)+0.67d0*phi(i,9)
	enddo
	do i=imt+2, ims, 2
	   do k=3, kmsi, 2
	      ztmp=(a-rn)-xaxel(rho(k),phi(i,k),alfayz)
	      r(i,k,1)=rhoz(ztmp,90.d0,alfa)
	   enddo
	   r(i,1,1)=r(i,3,1)
	enddo
	r(ims,3,1)=r(ims-2,3,1)
	r(ims,kms,1)=rhoz(a-r2,90.d0,alfa)
C
C.......Koordinater i X-Y-Z -planet skapas.
C . . . H|rnnoder och noder utmed sprickfronten. !!!
	do j=1, jms
	   lh2=hyperbel_length(r(ims,kms,1),r(imt,1,1),fi(j),alfa)
	   do k=1, kmsi,2
C  . . . . . X-Z Planet
	      do i=1,ims,2
	         npos(nods(i,j,k),1)=xaxel(r(i,k,1),fi(j),alfa)
	         npos(nods(i,j,k),3)=zaxel(r(i,k,1),fi(j),alfa)
	      enddo
	      if (k.eq.1) then
	         do i=1,imt
                    npos(nods(i,j,k),1)=xaxel(r(i,k,1),fi(j),alfa)
	            npos(nods(i,j,k),3)=zaxel(r(i,k,1),fi(j),alfa)
	         enddo
	      endif
C  . . . . . Y - Riktningen
	      if (j.lt.2) then
	         ro=npos(nods(1,j,k),3)-npos(nods(imt,j,1),3)
	         ron=npos(nods(1,j,k),3)-npos(nods(1,j,1),3)
	         ru=npos(nods(imt,j,1),3)-npos(nods(ims,j,k),3)
	      elseif (j.lt.jms) then
	         ro=hyperbel_length(r(imt,1,1),r(1,k,1),fi(j),alfa)
	         ron=hyperbel_length(r(1,1,1),r(1,k,1),fi(j),alfa)
	         ru=hyperbel_length(r(ims,k,1),r(imt,1,1),fi(j),alfa)
	      else
	         ro=npos(nods(1,j,k),1)-npos(nods(imt,j,1),1)
	         ron=npos(nods(1,j,k),1)-npos(nods(1,j,1),1)
	         ru=npos(nods(imt,j,1),1)-npos(nods(ims,j,k),1)
	      endif
C . . . . . Framf|r sprickan ( I =< IMT ).
C . . . . . . Mittnoderna i notchen kommer nu att ligga en cirkelb}ge.
	      npos(nods(1,j,k),2)=0.d0
	      do i=2, imt
	         npos(nods(i,j,k),2)=ro*dsind(phi(i,k))
	      enddo
C . . . . . Bakom  sprickan ( I > IMT ). En ellips utnyttjas d{r
C                RON {r lillaxeln och RU {r storaxeln.
	      if (k.gt.1) then
	         do i=imt+2, ims-2, 2
	            npos(nods(i,j,k),2) = ron*dsind(phi(i,k)) +
     &                         npos(nods(imt,j,1),2)
	         enddo
	         lh1=hyperbel_length(r(ims,k,1),r(imt,1,1),fi(j),alfa)
	         npos(nods(ims,j,k),2)= ( rn*lh1 +
     &                npos(nods(imt,j,1),2)*(lh2-lh1) ) / lh2
	      else
	        do i=imt+2, ims-2, 2
	          lh1=hyperbel_length(r(i,k,1),r(imt,1,1),fi(j),alfa)
	          npos(nods(i,j,k),2)= ( rn*lh1 +
     &                npos(nods(imt,j,1),2)*(lh2-lh1) ) / lh2
	        enddo
	      endif
	   enddo
	enddo
C . . . Mittnoder !!!
	do j=1, jms
C . . . .  Mid-nodes at K=constant:
	   do k=2, kmsi,2
	      do i=1,ims,2
	         npos(nods(i,j,k),1)=(npos(nods(i,j,k-1),1)+
     &                                npos(nods(i,j,k+1),1))/2.d0
	         npos(nods(i,j,k),2)=(npos(nods(i,j,k-1),2)+
     &                                npos(nods(i,j,k+1),2))/2.d0
	         npos(nods(i,j,k),3)=(npos(nods(i,j,k-1),3)+
     &                                npos(nods(i,j,k+1),3))/2.d0
	      enddo
	   enddo
C . . . .  Mid-nodes at I=constant:
	   do i=imt+1, ims, 2
	      npos(nods(i,j,1),1)=(npos(nods(i-1,j,1),1)+
     &                                npos(nods(i+1,j,1),1))/2.d0
	      npos(nods(i,j,1),2)=(npos(nods(i-1,j,1),2)+
     &                                npos(nods(i+1,j,1),2))/2.d0
	      npos(nods(i,j,1),3)=(npos(nods(i-1,j,1),3)+
     &                                npos(nods(i+1,j,1),3))/2.d0
	   enddo
	   do k=3, kmsi,2
	      do i=2,ims,2
	         npos(nods(i,j,k),1)=(npos(nods(i-1,j,k),1)+
     &                                npos(nods(i+1,j,k),1))/2.d0
	         npos(nods(i,j,k),2)=(npos(nods(i-1,j,k),2)+
     &                                npos(nods(i+1,j,k),2))/2.d0
	         npos(nods(i,j,k),3)=(npos(nods(i-1,j,k),3)+
     &                                npos(nods(i+1,j,k),3))/2.d0
	      enddo
	   enddo
C . . . .  Centroid or Mid-surface nodes:
	   do k=2, kmsi,2
	      do i=2,ims,2
	         npos(nods(i,j,k),1)= 0.25d0*
     &               ( npos(nods(i-1,j,k),1)+npos(nods(i+1,j,k),1) +
     &                 npos(nods(i,j,k-1),1)+npos(nods(i,j,k+1),1) )
	         npos(nods(i,j,k),2)= 0.25d0*
     &               ( npos(nods(i-1,j,k),2)+npos(nods(i+1,j,k),2) +
     &                 npos(nods(i,j,k-1),2)+npos(nods(i,j,k+1),2) )
	         npos(nods(i,j,k),3)= 0.25d0*
     &               ( npos(nods(i-1,j,k),3)+npos(nods(i+1,j,k),3) +
     &                 npos(nods(i,j,k-1),3)+npos(nods(i,j,k+1),3) )
	      enddo
	   enddo
	enddo
C
C... In the case of a sharp crack tip:
	if (rzero.eq.1) then
           rn=0.d0
C . . . . .Let all the nodes adjacent to the crack front be located
C          at the crack front
	   do j=1, jms
 	      do i=2, imt
	         npos(nods(i,j,1),1)=npos(nods(1,j,1),1)
	         npos(nods(i,j,1),2)=npos(nods(1,j,1),2)
	         npos(nods(i,j,1),3)=npos(nods(1,j,1),3)
	      enddo
	   enddo
	endif
C
C.......1 b)  Justering av noder pg av elementreduktioner i ZONS !
	if (sfred.gt.0) then
C . . . . The area where ksr1+2 =< k =< kmsi :
	   do j=1, jms
	      do k=ksr1+2, kmsi, 2
	         do i=3, ims, 4
	            npos(nods(i,j,k),1)= 0.5d0*
     &                (npos(nods(i-2,j,k),1)+npos(nods(i+2,j,k),1))
	            npos(nods(i,j,k),2)= 0.5d0*
     &                (npos(nods(i-2,j,k),2)+npos(nods(i+2,j,k),2))
	            npos(nods(i,j,k),3)= 0.5d0*
     &                (npos(nods(i-2,j,k),3)+npos(nods(i+2,j,k),3))
		 enddo
	      enddo
	      do k=ksr1+3, kmsi, 2
	         do i=3, ims, 4
	            npos(nods(i,j,k),1)= 0.25d0*
     &               ( npos(nods(i-2,j,k),1)+npos(nods(i+2,j,k),1) +
     &                 npos(nods(i,j,k-1),1)+npos(nods(i,j,k+1),1) )
	            npos(nods(i,j,k),2)= 0.25d0*
     &               ( npos(nods(i-2,j,k),2)+npos(nods(i+2,j,k),2) +
     &                 npos(nods(i,j,k-1),2)+npos(nods(i,j,k+1),2) )
	            npos(nods(i,j,k),3)= 0.25d0*
     &               ( npos(nods(i-2,j,k),3)+npos(nods(i+2,j,k),3) +
     &                 npos(nods(i,j,k-1),3)+npos(nods(i,j,k+1),3) )
		 enddo
	      enddo
	   enddo
C . .  The area where element reductions occur ksr1-2 =< k =< ksr1+2 :
	   left=.true.
	   ilocal=1
	   call reduce_coord_2_to_1(nods,ns1,ns2,ns3,ilocal,3,ims,4,
     &                        1,jms,1,ksr1,left)
	endif
C.......2 b) Lagertjocklekarna (Y-led) i Zon A dao L <= MV  :
	r2y = r2 * (kappa+1.d0) / (2.d0*kappa)
	if (r2y.lt.(0.7d0*r1+0.3d0*r2)) r2y=(0.7d0*r1+0.3d0*r2)
	y(0)=0.d0
	y(mv)=r2y*dsind(180.d0*dble(mv)/dble(mf/ids))
	if (mv.gt.1) then
	   y(1)=npos(nods(ims-2*ids,1,kmsi),2)
	   dy=0.01d0*y(1)
10	   continue
	   lambda=calc_lambda(y(1),y(mv),mv)
	   if (dexp(lambda).gt.(1.5d0)) then
	      y(1)=y(1)+dy
	      goto 10
	   endif
	   if (dexp(lambda).gt.(1.0d0)) then
	      do i=2, mv-1
	        y(i)=y(i-1)+y(1)*dexp(lambda*dble(i-1))
	      enddo
	   else
	      dy=y(mv)/dble(mv)
	      do i=2, mv-1
	         y(i)=y(i-1)+dy
	      enddo
	   endif
	endif
C.......2 c) Den yttre "cirkeln" !
C . . . Radier f|r fi=90 ber{knas.
C . . . I =< 2*MV+1, elementen forts{tter i Z-led
	do i=1, 2*mv+1, 2
           ic=ids*(i-1)+1
	   r(ims+1-ic,kms,1)=rhoz((a-rn)-r2*dsqrt(1.d0-(y((i-1)/2)/
     &                       r2y)**2.d0),90.d0,alfa)
	   r(ic,kms,1)=2.d0*r(imt,1,1)-r(ims+1-ic,kms,1)
	   npos(nods(ic,1,kms),2)=y((i-1)/2)
	   npos(nods(ims+1-ic,1,kms),2)=y((i-1)/2)
	enddo
C . . . I >= 2*MV+1, elementen forts{tter i Y-led
	do i=2*mv+3, mf/ids-1
           ic=ids*(i-1)+1
	   teta=180.d0*dble(ic-1)/dble(ims-1)
	   rho(1)=1.d0/dsqrt( (dcosd(teta)/r2)**2.d0
     &              +(dsind(teta)/r2y)**2.d0)
	   r(ims+1-ic,kms,1)=rhoz((a-rn)-rho(1)*dcosd(teta),90.d0,alfa)
	   r(ic,kms,1)=2.d0*r(imt,1,1)-r(ims+1-ic,kms,1)
	   npos(nods(ic,1,kms),2)=rho(1)*dsind(teta)
	   npos(nods(ims+1-ic,1,kms),2)=npos(nods(ic,1,kms),2)
	enddo
	r(mf+1,kms,1)=r(imt,1,1)
	npos(nods(mf+1,1,kms),2)=r2y
C . . . H|rnnoder i Zon S K=KMS !
	do j=1, jms, 2
	   do i=1, ims, 2*ids
	      npos(nods(i,j,kms),1)=xaxel(r(i,kms,1),fi(j),alfa)
	      npos(nods(i,j,kms),2)=npos(nods(i,1,kms),2)
	      npos(nods(i,j,kms),3)=zaxel(r(i,kms,1),fi(j),alfa)
	   enddo
	   npos(nods(1,j,kms),2)=0.
	   npos(nods(ims,j,kms),2)=rn
	enddo
C.... Justering koordinatpos. da M1-MH =< 2 , beror 2(M1+MH)+1< I < IMS
	if ((m1-mh).lt.3) then
	   do i=2*(mf/ids-mv)+1, ims, 2
              ic=ids*(i-1)+1
	      do j=jms-2, jms, 2
	         npos(nods(ic,j,kms),1)=npos(nods(ic,jms-4,kms),1)
	      enddo
	      npos(nods(ic,jms-2,kms),3)=0.5*npos(nods(ic,jms-4,kms),3)
	   enddo
	endif
 
C . . . Mittnoder i Zon S K=KMS !
	do j=2, jms, 2
	   do i=1, ims, ids
	      npos(nods(i,j,kms),1)=(npos(nods(i,j-1,kms),1)+
     &                               npos(nods(i,j+1,kms),1))/2.d0
	      npos(nods(i,j,kms),2)=(npos(nods(i,j-1,kms),2)+
     &                               npos(nods(i,j+1,kms),2))/2.d0
	      npos(nods(i,j,kms),3)=(npos(nods(i,j-1,kms),3)+
     &                               npos(nods(i,j+1,kms),3))/2.d0
	   enddo
	enddo
	do j=1, jms
	   do i=1, ims, 2*ids
	      npos(nods(i,j,kms-1),1)=(npos(nods(i,j,kmsi),1)+
     &                                 npos(nods(i,j,kms),1))/2.d0
	      npos(nods(i,j,kms-1),2)=(npos(nods(i,j,kmsi),2)+
     &                                 npos(nods(i,j,kms),2))/2.d0
	      npos(nods(i,j,kms-1),3)=(npos(nods(i,j,kmsi),3)+
     &                                 npos(nods(i,j,kms),3))/2.d0
	   enddo
	   do i=1+ids, ims, 2*ids
	      npos(nods(i,j,kms),1)=(npos(nods(i-ids,j,kms),1)+
     &                               npos(nods(i+ids,j,kms),1))/2.d0
	      npos(nods(i,j,kms),2)=(npos(nods(i-ids,j,kms),2)+
     &                               npos(nods(i+ids,j,kms),2))/2.d0
	      npos(nods(i,j,kms),3)=(npos(nods(i-ids,j,kms),3)+
     &                               npos(nods(i+ids,j,kms),3))/2.d0
	   enddo
C . . . . Surface nodes and the centriod node.
	   do i=1+ids, ims, 2*ids
	      npos(nods(i,j,kms-1),1)= 0.25d0*
     &        (npos(nods(i-ids,j,kms-1),1)+npos(nods(i+ids,j,kms-1),1)+
     &            npos(nods(i,j,kmsi),1)+npos(nods(i,j,kms),1) )
	      npos(nods(i,j,kms-1),2)= 0.25d0*
     &        (npos(nods(i-ids,j,kms-1),2)+npos(nods(i+ids,j,kms-1),2)+
     &            npos(nods(i,j,kmsi),2)+npos(nods(i,j,kms),2) )
	      npos(nods(i,j,kms-1),3)= 0.25d0*
     &        (npos(nods(i-ids,j,kms-1),3)+npos(nods(i+ids,j,kms-1),3)+
     &            npos(nods(i,j,kmsi),3)+npos(nods(i,j,kms),3) )
	   enddo
	enddo
C
C.......2 d) Justera randerna i Zon S !
C . . . Symmetrisnittet (x=0) !
	do k=1, kms
	   do i=1, ims
	      npos(nods(i,1,k),1)=0.d0
	   enddo
	enddo
C . . . Undersidan (z=0) !
	do k=1, kms
	   do i=1, ims
	      npos(nods(i,jms,k),3)=0.d0
	   enddo
	enddo
C . . . Sprickplanet justeras (y=0) !
	do j=1, jms
	   do k=1, kms
	      npos(nods(1,j,k),2)=0.d0
	   enddo
	enddo
	if (rzero.eq.1) then
	   do j=1, jms
	      do i=2, ims-2
	         npos(nods(i,j,1),2)=0.d0
	      enddo
	      do k=3, kms
	         npos(nods(ims,j,k),2)=0.d0
	      enddo
	   enddo
	endif
C
C--------- For Testing the Program ! -------------------------
C.......2 f) => SKALL noderna i ZON S skrivas ut ?
C	WRITE(*,'(T10,A,$)') '* Skall noderna i Zon S plottas (J/N) : '
C	READ(*,'(A)') SVAR
C	IF ((SVAR.EQ.'J').OR.(SVAR.EQ.'j')) THEN
C	   CALL PLOT_NOD_ZONS(NODS,NS1,NS2,NS3,FI)
C	   STOP
C	ENDIF
C
C---------------------------------------------------
C.......3) NODKOORDINATER I ZON A
C---------------------------------------------------
C . . . a) Best{m Z-koordinat f|r fi=90 !
	ia1=2*(m1-mh)+1
	ia2=2*(m1+mh)+1
	ia3=ia2+2
C . . . K=1 !!
	zik(1,1)=0.d0
	zik(ia1,1)=npos(noda(ia1,1,1),3)
	zik(2*m1+1,1)=a-rn
	zik(ia2,1)=npos(noda(ia2,1,1),3)
	zik(ima,1)=t
C       WRITE(*,'(T3,A,I2,A,G12.6)') ' ZIK(i=',1, ')=' ,ZIK(1,1)
C       WRITE(*,'(T3,A,I2,A,G12.6)') ' ZIK(ia1=',IA1, ')=' ,ZIK(IA1,1)
C       WRITE(*,'(T3,A,I2,A,G12.6)') ' ZIK(is=',IS, ')=' ,ZIK(IS,1)
C       WRITE(*,'(T3,A,I2,A,G12.6)') ' ZIK(ia2=',IA2, ')=' ,ZIK(IA2,1)
C       WRITE(*,'(T3,A,I2,A,G12.6)') ' ZIK(ima=',IMA, ')=' ,ZIK(IMA,1)
	dz=npos(nods(ims,1,kmsi),3)-zik(ia1,1)
	lambda=calc_lambda( dz, npos(nods(ims,1,kmsi),3),(m1-mh+1))
	if (dexp(lambda).gt.(1.d0)) then
	   it=0
	   do i=ia1-2, 3, -2
	      it=it+1
	      zik(i,1)=zik(i+2,1)-dz*dexp(lambda*dble(it))
C       WRITE(*,'(T3,A,I2,A,G12.6)') ' ZIK(i=',I, ')=' ,ZIK(I,1)
	   enddo
	else
           dz=(zik(ia1,1)-zik(1,1))/dble(m1-mh)
	   do i=3, ia1-2, 2
	      zik(i,1)=zik(i-2,1)+dz
	   enddo
	endif
	dz=zik(ia2,1)-npos(nods(1,1,kmsi),3)
	lambda=calc_lambda( dz,(zik(ima,1)-npos(nods(1,1,kmsi),3)),
     &                     (m2-mh+1) )
	if (dexp(lambda).gt.(1.d0)) then
	   it=0
	   do i=ia3, ima-2, 2
	      it=it+1
	      zik(i,1)=zik(i-2,1)+dz*dexp(lambda*dble(it))
C       WRITE(*,'(T3,A,I2,A,G12.6)') ' ZIK(i=',I, ')=' ,ZIK(I,1)
	   enddo
 	else
           dz=(zik(ima,1)-zik(ia2,1))/dble(m2-mh)
	   do i=ia2+2, ima-2, 2
	      zik(i,1)=zik(i-2,1)+dz
	   enddo
	endif
C ..... 1 < K < 2*(MV+1)+1 !!!!!!
C . . . . .Parablar skapas kring sprickfronten. :
C . . . . . .Under sprickfronten
	ka1=2*(mv+1)+1
	ny=1.d0
	dny=1.d0
	if ((m1-mh).gt.2) dny=dny/dble(m1-mh-2)
	do i=ia1-2, 3, -2
	   do k=3, ka1, 2
	      ai=(a-rn)-zik(i,1)
	      rot=1.d0-(y((k-1)/2)/(ny*ai))**2.d0
	      if (rot.lt.0) then
	         write(*,'(t5,a)')
     &          '=> the y-coordinate at l=mv+1 is to large, decrease!'
	         stop
	      endif
	      zik(i,k)=(a-rn)-ai*dsqrt(rot)
	   enddo
	   ny=ny+dny
	enddo
C . . . Over sprickfronten
	ny=1.d0
	dny=1.d0
	if ((m2-mh).gt.2) dny=dny/dble(m2-mh-2)
	do i=ia2+2, ima-2, 2
	   do k=3, ka1, 2
	      ai=zik(i,1)-(a-rn)
	      rot=1.d0-(y((k-1)/2)/(ny*ai))**2.d0
	      if (rot.lt.0.d0) then
	         write(*,'(t5,a)')
     &          '=> the y-coordinate at l=mv+1 is to large, decrease!'
	         stop
	      endif
	      zik(i,k)=(a-rn)+ai*dsqrt(rot)
	   enddo
	   ny=ny+dny
	enddo
C . . . Mellan omradet ! ! !
	do k=3, 2*mv+1, 2
	   zik(ia1,k)=npos(noda(ia1,1,k),3)
	   zik(ia2,k)=npos(noda(ia2,1,k),3)
	enddo
	do i=ia1+2, ia2-2, 2
	   zik(i,2*mv+1)=npos(noda(i,1,2*mv+1),3)
	enddo
	dz=(zik(ia2+2,ka1)-zik(ia1-2,ka1))/dble(2*mh+2)
	do i=ia1, ia2, 2
	   zik(i,ka1)=zik(i-2,ka1)+dz
	enddo
C...... K = 2*LRED+1 !!!!!!
C . . Elementen skall ha samma langd i Z-led for fi=90. & K >= 2*lred+1
	ka2=2*(lred-1)+1
	do k=1, ka2, 2
	   zik(1,k)=0.d0
	   zik(ima,k)=t
	enddo
	dz=t/dble(ma)
	do i=3, ima-2, 2
	   zik(i,(ka2))=zik(i-2,(ka2))+dz
	enddo
C.... 2*(MV+1)+1  <  K  >  2*LRED+1  !!!!!!
	do k=ka1+2, ka2-2, 2
	   psi=(y((k-1)/2)-y((ka1-1)/2))/(y((ka2-1)/2)-y((ka1-1)/2))
	   do i=3, ima, 2
	      zik(i,k)=(1.d0-psi)*zik(i,ka1)+psi*zik(i,ka2)
	   enddo
	enddo
 
C-------Vid Test av Programmet ! -------------------------
C	WRITE(*,'(T5,A,$)') '* Vill du plotta ZIK(i,k) (J/N) : '
C	READ(*,'(A)') SVAR
C	IF ((SVAR.EQ.'J').OR.(SVAR.EQ.'j')) THEN
C	   J=0
C	   OPEN(UNIT=33,FILE='ZIK.DAT',STATUS='UNKNOWN')
C	   DO K=1, KA2, 2
C	      DO I=1, IMA, 2
C	         IF (NODA(I,1,K).GT.0) THEN
C	            WRITE(33,*) Y((K-1)/2), ZIK(I,K)
C	            J=J+1
C	         ENDIF
C	      ENDDO
C	   ENDDO
C	   WRITE(*,'(T5,A,I3,A)') '=> ',J,' st koord. finns pa ZIK.DAT'
C	   CLOSE(33)
C	   STOP
C	ENDIF
C-----------------------------------------------
C
C.......3 b) Bestam Radierna r(i,j,k) for fi=90.
	do k=1, ka2, 2
C	   WRITE(*,*) '  * K =',K
	   do i=3, ima, 2
	      if (noda(i,1,k).gt.0) then
                 r(i,1,k)=rhoz(zik(i,k),90.d0,alfa)
C	    WRITE(*,'(A,I2,A,G12.6)')'  * R(i=',I,',1,k)',R(I,1,K)
	      endif
	   enddo
	enddo
C . . . Radien ar oberoende av fi for i =< 2*(M1+MH+1)+1
	do j=3, jma, 2
	   do k=1, ka2, 2
	      do i=3, ima, 2
	         if (noda(i,j,k).gt.0) r(i,j,k)=r(i,1,k)
	      enddo
	   enddo
	enddo
C . . . Radierna ( i > 2*(M1+MH+1)+1 ) ar beroende av fi !
C . . . Justera Radierna for i > 2*(M1+MH+1)+1
C . . . Borja med I=IMA
	ja1=2*(na-nb)+1
	do j=1, ja1, 2
	   do k=1, ka2, 2
	      r(ima,j,k)=rhoz(t,fi(j),alfa)
	   enddo
	enddo
	if ( xaxel(r(ima,ja1,1),fi(ja1),alfa).gt.
     &     xaxel(r(ima,jma,1),0.d0,alfa) ) then
	   do k=1, ka2, 2
	      do j=ja1+2, jma, 2
	         r(ima,j,k)=
     &             rhox( xaxel(r(ima,ja1,k),fi(ja1),alfa),fi(j),alfa)
	      enddo
	   enddo
	else
	   xe=xaxel(r(ima,jma,1),0.d0,alfa)
	   ze=t/dsqrt(1.d0-(xaxel(r(ima,ja1,1),fi(ja1),alfa)/xe)**2.d0)
	   do k=1, ka2, 2
	      do j=ja1+2, jma-2, 2
	         r(ima,j,k)=rho_ellips(xe,ze,fi(j),alfa)
	      enddo
	   enddo
	endif
C . . . 2*(M1+MH+2)+1  =<  I  =<  IMA-2
	do k=1, ka2, 2
	   dr(1)=r(ia3,ja1,k)-r(ia2,ja1,k)
	   lambda=calc_lambda(dr(1),r(ima,ja1,k)-r(ia2,ja1,k),m2-mh)
	   it=0
	   do i=ia3+2, ima-2, 2
	      it=it+1
	      r(i,ja1,k)=r(i-2,ja1,k)+dr(1)*dexp(lambda*dble(it))
	      ze=zaxel(r(i,1,k),90.d0,alfa)
	      if (zaxel(r(i,ja1,k),fi(ja1),alfa).gt.ze) then
                 do j=3, ja1, 2
	  	    r(i,j,k)=rhoz(ze, fi(j), alfa)
	         enddo
	      else
	         xe=xaxel(r(i,ja1,k),fi(ja1),alfa)/
     &            dsqrt(1.d0-(zaxel(r(i,ja1,k),fi(ja1),alfa)/ze)**2.d0)
	         do j=3, ja1-2, 2
	            r(i,j,k)=rho_ellips(xe,ze,fi(j),alfa)
	         enddo
	      endif
	      if ( xaxel(r(i,ja1,k),fi(ja1),alfa).gt.
     &           xaxel(r(i,jma,k),0.d0,alfa) ) then
	         do j=ja1+2, jma, 2
	            r(i,j,k)=
     &               rhox( xaxel(r(i,ja1,k),fi(ja1),alfa), fi(j), alfa)
	         enddo
	      else
	         xe=xaxel(r(i,jma,k),0.d0,alfa)
	         ze=zaxel(r(i,ja1,k),fi(ja1),alfa)/
     &            dsqrt(1.d0-(xaxel(r(i,ja1,k),fi(ja1),alfa)/xe)**2.d0)
	         do j=ja1+2, jma-2, 2
	            r(i,j,k)=rho_ellips(xe,ze,fi(j),alfa)
	         enddo
	      endif
	   enddo
	enddo
C
C.......3 c) X-Y-Z  Koordinaterna bestams !!!
C . . . Hornnoder !!!
C . . . K =< 2*(LRED-1) +1
C . . . i >= 5 . . . . . . . . . .
	ja2=jma-4
	do k=1, ka2, 2
	   do i=5, ima, 2
	      do j=1, jma, 2
	         if ( noda(i,j,k).gt.kstart ) then
	            npos(noda(i,j,k),1)=xaxel(r(i,j,k),fi(j),alfa)
	            npos(noda(i,j,k),2)=y((k-1)/2)
	            npos(noda(i,j,k),3)=zaxel(r(i,j,k),fi(j),alfa)
	         endif
	      enddo
	   enddo
	enddo
 
C.......Låt Yka1=(psi*Y(ka1-2) + (1-psi)*Y(ka1+2) ) för ia1=< i =< ia2
C......där psi=(yka1 - y(ka1-2)(i=ima) ) / ( y(ka1+2) - y(ka1-2) )
	psi= (y((ka1-1)/2)-y((ka1-3)/2)) / (y((ka1+1)/2)-y((ka1-3)/2))
	do j=1, jma, 2
	   do i=ia1, ia2, 2
	      npos(noda(i,j,ka1),2)=(1.d0-psi)*npos(noda(i,j,ka1-2),2) +
     &                psi*npos(noda(i,j,ka1+2),2)
	   enddo
	enddo
C . . . i < 5 . . . . . . . . . .
	do k=1, ka2, 2
	   psi=zik(3,k)/zik(5,k)
	   do j=1, jma-6, 2
	      npos(noda(1,j,k),1)=npos(noda(5,j,k),1)
              npos(noda(1,j,k),2)=y((k-1)/2)
	      npos(noda(1,j,k),3)=0.d0
	      npos(noda(3,j,k),1)=npos(noda(5,j,k),1)
              npos(noda(3,j,k),2)=y((k-1)/2)
	      npos(noda(3,j,k),3)=npos(noda(5,j,k),3)*
     &             (psi*dble(ja2-j)+0.5d0*dble(j))/dble(ja2)
	   enddo
	   if (noda(5,ja2,k).gt.kstart) then
	      npos(noda(1,ja2,k),1)=npos(noda(5,ja2,k),1)
              npos(noda(1,ja2,k),2)=y((k-1)/2)
	      npos(noda(1,ja2,k),3)=0.d0
	      npos(noda(3,ja2,k),1)=npos(noda(5,ja2,k),1)
              npos(noda(3,ja2,k),2)=y((k-1)/2)
	      npos(noda(3,ja2,k),3)=npos(noda(5,j,k),3)*0.5d0
	   endif
	enddo
C . . . Korrigera nodpositioner f|r K=1 & I < IA1.
        do j=1, jma
	   do i=1, ia1
	      if ( noda(i,j,k).gt.kstart ) then
	         npos(noda(i,j,1),2)=rn
	      endif
	   enddo
	enddo
C
C . . Justera hörnnoder vid elementreducering, dvs da K > 2*(lred-1) !
	if (rtype.eq.1) then
	   do i=5, ima, 4
	      npos(noda(i,jma-2,ka2),1) =(npos(noda(i,jma-4,ka2),1)+
     &                   npos(noda(i,jma,ka2),1) ) / 2.d0
	      npos(noda(i,jma-2,ka2),3) =(npos(noda(i,jma-4,ka2),3)+
     &                   npos(noda(i,jma,ka2),3) ) / 2.d0
	   enddo
	   do j=1, jma, 2
	      do i=3, ima, 4
	         npos(noda(i,j,ka2),1) = ( npos(noda(i-2,j,ka2),1) +
     &                   npos(noda(i+2,j,ka2),1) ) / 2.d0
	         npos(noda(i,j,ka2),3) = ( npos(noda(i-2,j,ka2),3) +
     &                   npos(noda(i+2,j,ka2),3) ) / 2.d0
	      enddo
	   enddo
	endif
	if (rtype.eq.2) then
	   do j=3, jma, 4
	      do i=1, ima ,2
	         npos(noda(i,j,ka2),1) = ( npos(noda(i,j-2,ka2),1) +
     &                   npos(noda(i,j+2,ka2),1) ) / 2.d0
	         npos(noda(i,j,ka2),3) = ( npos(noda(i,j-2,ka2),3) +
     &                   npos(noda(i,j+2,ka2),3) ) / 2.d0
	      enddo
	   enddo
	   do j=1, jma, 2
	      do i=3, ima, 4
	         npos(noda(i,j,ka2),1) = ( npos(noda(i-2,j,ka2),1) +
     &                   npos(noda(i+2,j,ka2),1) ) / 2.d0
	         npos(noda(i,j,ka2),3) = ( npos(noda(i-2,j,ka2),3) +
     &                   npos(noda(i+2,j,ka2),3) ) / 2.d0
	      enddo
	   enddo
	endif
C... MID-, SURFACE- and CETROID-NODES !!!
C . . . Midnodes for even k:s
	do k=2, ka2, 2
	   do i=1, ima, 2
	      do j=1, jma, 2
	         npos(noda(i,j,k),1) = ( npos(noda(i,j,k-1),1) +
     &                   npos(noda(i,j,k+1),1) ) / 2.d0
	         npos(noda(i,j,k),2) = ( npos(noda(i,j,k-1),2) +
     &                   npos(noda(i,j,k+1),2) ) / 2.d0
	         npos(noda(i,j,k),3) = ( npos(noda(i,j,k-1),3) +
     &                   npos(noda(i,j,k+1),3) ) / 2.d0
	      enddo
	   enddo
	enddo
C . . . Midnodes for even i:s
	do k=1, ka2
	   do i=2, ima, 2
	      do j=1, jma, 2
	         npos(noda(i,j,k),1) = ( npos(noda(i-1,j,k),1) +
     &                   npos(noda(i+1,j,k),1) ) / 2.d0
	         npos(noda(i,j,k),2) = ( npos(noda(i-1,j,k),2) +
     &                   npos(noda(i+1,j,k),2) ) / 2.d0
	         npos(noda(i,j,k),3) = ( npos(noda(i-1,j,k),3) +
     &                   npos(noda(i+1,j,k),3) ) / 2.d0
	      enddo
	   enddo
	enddo
C . . . Midnodes for even j:s
	do k=1, ka2
	   do i=1, ima, 2
	      do j=2, jma,2
	         npos(noda(i,j,k),1) = ( npos(noda(i,j-1,k),1) +
     &                   npos(noda(i,j+1,k),1) ) / 2.d0
	         npos(noda(i,j,k),2) = ( npos(noda(i,j-1,k),2) +
     &                   npos(noda(i,j+1,k),2) ) / 2.d0
	         npos(noda(i,j,k),3) = ( npos(noda(i,j-1,k),3) +
     &                   npos(noda(i,j+1,k),3) ) / 2.d0
	      enddo
	   enddo
	enddo
C . . . Surface and centroid nodes for even i:s and j:s
	do k=1, ka2
	   do i=2, ima, 2
	      do j=2, jma, 2
	         npos(noda(i,j,k),1) = 0.25d0 *
     &               ( npos(noda(i-1,j,k),1) + npos(noda(i+1,j,k),1) +
     &                 npos(noda(i,j-1,k),1) + npos(noda(i,j+1,k),1) )
	         npos(noda(i,j,k),2) = 0.25d0 *
     &               ( npos(noda(i-1,j,k),2) + npos(noda(i+1,j,k),2) +
     &                 npos(noda(i,j-1,k),2) + npos(noda(i,j+1,k),2) )
	         npos(noda(i,j,k),3) = 0.25d0 *
     &               ( npos(noda(i-1,j,k),3) + npos(noda(i+1,j,k),3) +
     &                 npos(noda(i,j-1,k),3) + npos(noda(i,j+1,k),3) )
	      enddo
	   enddo
	enddo
C . .  The domain under the crack front: i =< 5 and 1 =< j =< jma-4
        if (mod(na,4).ne.0) then
	   left=.true.
	else
	   left=.false.
	endif
	call reduce_coord_under(noda,na1,na2,na3,3,jma-4,4,1,ka2,1,left)
 
C . . The domain in the possibly reduced region for k > ka2
	do k=ka2+2, kma, 2
	   do j=1, jma
	      do i=1, ima
	         npos(noda(i,j,k),1)=npos(noda(i,j,ka2),1)
	         npos(noda(i,j,k),2)=y((k-1)/2)
	         npos(noda(i,j,k),3)=npos(noda(i,j,ka2),3)
 
	         npos(noda(i,j,k-1),1)=npos(noda(i,j,ka2),1)
	         npos(noda(i,j,k-1),2)=(y((k-1)/2)+y((k-3)/2))/2.d0
	         npos(noda(i,j,k-1),3)=npos(noda(i,j,ka2),3)
	      enddo
	   enddo
	enddo
 
	if (rtype.ge.1) then
	   if (mod(ma,4).eq.0) then
	      left=.true.
	   else
	      left=.false.
	   endif
	   ilocal=1
	   call reduce_coord_2_to_1(noda,na1,na2,na3,ilocal,7,ima,4,
     &                        1,jma,1,kar1,left)
	endif
	if (rtype.eq.2) then
	   if (mod(na,4).ne.0) then
	      left=.true.
	   else
	      left=.false.
	   endif
	   ilocal=2
	   call reduce_coord_2_to_1(noda,na1,na2,na3,ilocal,5,ima,1,
     &                        3,jma,4,kar2,left)
C . . . . .The small limited area under the crack front
	   if (mod(na,4).ne.0) then
              left=.true.
	   else
              left=.false.
	   endif
	   call reduce_coord_3_to_1(noda,na1,na2,na3,3,jma,4,kar2,left)
	endif
 
C.......3 d) Justera randerna i Zon A !
C . . . Symmetri snittet justeras (x=0) !
	do k=1, kma
	   do i=1, ima
	      npos(noda(i,1,k),1)=0.d0
	   enddo
	enddo
C . . . Over- (z=T) och under- (z=0) sidan justeras !
	do k=1, kma
	   do j=1, jma-jmb+1
	      npos(noda(ima,j,k),3)=t
	   enddo
	   do j=1, jma-4
	      npos(noda(1,j,k),3)=0.d0
	   enddo
	   do i=5, ima
	      npos(noda(i,jma,k),3)=0.d0
	   enddo
	enddo
C--------- For Test av Programmet ! -------------------------
C.......3 e) => SKALL noderna i ZON A skrivas ut ?
C	WRITE(*,'(T10,A,$)') '* Skall noderna i Zon A plottas (J/N) : '
C	READ(*,'(A)') SVAR
C	IF ((SVAR.EQ.'J').OR.(SVAR.EQ.'j')) THEN
C	   CALL PLOT_NOD_ZONA(NODA,NA1,NA2,NA3)
C           STOP
C	ENDIF
 
C---------------------------------------------------
C.......4) NODKOORDINATER I ZON B
C---------------------------------------------------
 
C.......4 a) Bestam nodpos for y=z=0 (hornnoder).
C            Anvand ZIK(i,j) & DR(i) som hjalpvariabler !
C . . . j=jmb
	if (zb_xw.eq.0) then
	   dr(1)=npos(noda(ima,jma,1),1)-npos(noda(ima-2,jma,1),1)
	   lambda=calc_lambda(dr(1),w-npos(noda(ima-2,jma,1),1),mb+1)
	else
	   dr(1)=(w-npos(noda(ima,jma,1),1))*
     &         (1.d0-rate_b)/(rate_b*(1.d0-rate_b**dble(mb)))
	   lambda=dlog(rate_b)
	endif
	it=0
	do i=3, imb, 2
	   it=it+1
	   npos(nodb(i,jmb,1),1)=npos(nodb(i-2,jmb,1),1)+
     &              dr(1)*dexp(lambda*dble(it))
	   npos(nodb(i,jmb,1),2)=y(0)
	   npos(nodb(i,jmb,1),3)=0.d0
	enddo
	npos(nodb(imb,jmb,1),1)=w
C . . . j=1
	do i=5, imb, 2
	   npos(nodb(i,1,1),1)=npos(nodb(i,jmb,1),1)
	   npos(nodb(i,1,1),2)=y(0)
	   npos(nodb(i,1,1),3)=t
	enddo
	npos(nodb(3,1,1),1)=( npos(nodb(5,1,1),1)+
     &    dexp(lambda)*npos(nodb(1,1,1),1) )/(1.d0+dexp(lambda))
	npos(nodb(3,1,1),2)=y(0)
	npos(nodb(3,1,1),3)=t
C.......j >= 5  &  1 < i < imb  ! !
	do i=5, imb, 2
	   do j=3, jmb-2, 2
	      npos(nodb(i,j,1),1)=npos(nodb(i,jmb,1),1)
	      npos(nodb(i,j,1),2)=y(0)
	      npos(nodb(i,j,1),3)=t*dble(jmb-j)/dble(2*nb)
	   enddo
	enddo
C.......j = 3  &  1 < i < imb  ! !
	do j=3, jmb-2, 2
           npos(nodb(3,j,1),3)=
     &      ( npos(nodb(1,j,1),3)+npos(nodb(5,j,1),3) ) / 2.d0
           npos(nodb(3,j,1),2)=y(0)
	enddo
	xe=npos(nodb(3,jmb,1),1)
	rot=1.d0-(npos(nodb(3,1,1),1)/xe)**2.d0
	if (rot.gt.0.d0) then
	  ze=t/dsqrt(rot)
	  do j=3, jmb-2, 2
	   npos(nodb(3,j,1),1)=
     &           xe*dsqrt(1.d0-(npos(nodb(3,j,1),3)/ze)**2.d0)
	  enddo
	else
	  do j=3, jmb-2, 2
	   npos(nodb(3,j,1),1)=xe
	  enddo
	endif
C.......4 b) Skapa alla Hornnoder for k : (3, 5, ... , ka2)
	do k=3, ka2, 2
	   do i=3, imb, 2
	      do j=1, jmb, 2
	         npos(nodb(i,j,k),1)=npos(nodb(i,j,1),1)
	         npos(nodb(i,j,k),2)=y((k-1)/2)
	         npos(nodb(i,j,k),3)=npos(nodb(i,j,1),3)
	      enddo
	   enddo
	enddo
C . . . Justera k=ka2 om element reducering i zon B ! ! !
	if (rtype.eq.2) then
	   do i=3, imb, 2
	      do j=3, jmb-2, 4
	         npos(nodb(i,j,ka2),1) = ( npos(nodb(i,j-2,ka2),1) +
     &                     npos(nodb(i,j+2,ka2),1) ) / 2.d0
	         npos(nodb(i,j,ka2),3) = ( npos(nodb(i,j-2,ka2),3) +
     &                     npos(nodb(i,j+2,ka2),3) ) / 2.d0
	      enddo
	   enddo
	endif
C
C... MID-, SURFACE- and CETROID-NODES !!!
C . . . Midnodes for even k:s
	do k=2, ka2, 2
	   do i=3, imb, 2
	      do j=1, jmb, 2
	         npos(nodb(i,j,k),1) = ( npos(nodb(i,j,k-1),1) +
     &                   npos(nodb(i,j,k+1),1) ) / 2.d0
	         npos(nodb(i,j,k),2) = ( npos(nodb(i,j,k-1),2) +
     &                   npos(nodb(i,j,k+1),2) ) / 2.d0
	         npos(nodb(i,j,k),3) = ( npos(nodb(i,j,k-1),3) +
     &                   npos(nodb(i,j,k+1),3) ) / 2.d0
	      enddo
	   enddo
	enddo
C . . . Midnodes for even i:s
	do k=1, ka2
	   do i=2, imb, 2
	      do j=1, jmb, 2
	         npos(nodb(i,j,k),1) = ( npos(nodb(i-1,j,k),1) +
     &                   npos(nodb(i+1,j,k),1) ) / 2.d0
	         npos(nodb(i,j,k),2) = ( npos(nodb(i-1,j,k),2) +
     &                   npos(nodb(i+1,j,k),2) ) / 2.d0
	         npos(nodb(i,j,k),3) = ( npos(nodb(i-1,j,k),3) +
     &                   npos(nodb(i+1,j,k),3) ) / 2.d0
	      enddo
	   enddo
	enddo
C . . . Midnodes for even j:s
	do k=1, ka2
	   do i=3, imb, 2
	      do j=2, jmb,2
	         npos(nodb(i,j,k),1) = ( npos(nodb(i,j-1,k),1) +
     &                   npos(nodb(i,j+1,k),1) ) / 2.d0
	         npos(nodb(i,j,k),2) = ( npos(nodb(i,j-1,k),2) +
     &                   npos(nodb(i,j+1,k),2) ) / 2.d0
	         npos(nodb(i,j,k),3) = ( npos(nodb(i,j-1,k),3) +
     &                   npos(nodb(i,j+1,k),3) ) / 2.d0
	      enddo
	   enddo
	enddo
C . . . Surface and centroid nodes for even i:s and j:s
	do k=1, ka2
	   do i=2, imb, 2
	      do j=2, jmb, 2
	         npos(nodb(i,j,k),1) = 0.25d0 *
     &               ( npos(nodb(i-1,j,k),1) + npos(nodb(i+1,j,k),1) +
     &                 npos(nodb(i,j-1,k),1) + npos(nodb(i,j+1,k),1) )
	         npos(nodb(i,j,k),2) = 0.25d0 *
     &               ( npos(nodb(i-1,j,k),2) + npos(nodb(i+1,j,k),2) +
     &                 npos(nodb(i,j-1,k),2) + npos(nodb(i,j+1,k),2) )
	         npos(nodb(i,j,k),3) = 0.25d0 *
     &               ( npos(nodb(i-1,j,k),3) + npos(nodb(i+1,j,k),3) +
     &                 npos(nodb(i,j-1,k),3) + npos(nodb(i,j+1,k),3) )
	      enddo
	   enddo
	enddo
C . . The domain in the possibly reduced region for k > ka2
	do k=ka2+2, kma, 2
	   do j=1, jmb
	      do i=2, imb
	         npos(nodb(i,j,k),1)=npos(nodb(i,j,ka2),1)
	         npos(nodb(i,j,k),2)=y((k-1)/2)
	         npos(nodb(i,j,k),3)=npos(nodb(i,j,ka2),3)
 
	         npos(nodb(i,j,k-1),1)=npos(nodb(i,j,ka2),1)
	         npos(nodb(i,j,k-1),2)=(y((k-1)/2)+y((k-3)/2))/2.d0
	         npos(nodb(i,j,k-1),3)=npos(nodb(i,j,ka2),3)
	      enddo
	   enddo
	enddo
C . . If reduction of elements in Zon B:
	if (rtype.eq.2) then
	   if (mod(nb,4).ne.0) then
	      left=.true.
	   else
	      left=.false.
	   endif
	   ilocal=2
	   call reduce_coord_2_to_1(nodb,nb1,nb2,nb3,ilocal,2,imb,1,
     &                    3,jmb,4,kar2,left)
	endif
C
C--------- For Test av Programmet ! -------------------------
C.......4 d) => SKALL noderna i ZON B skrivas ut ?
C	WRITE(*,'(T10,A,$)') '* Skall noderna i Zon B plottas (J/N) : '
C	READ(*,'(A)') SVAR
C	IF ((SVAR.EQ.'J').OR.(SVAR.EQ.'j')) THEN
C	   CALL PLOT_NOD_ZONB(NODB,NB1,NB2,NB3)
C	ENDIF
 
C-----------------------------------------------------------!
C   KOORDINATTRANSFORMATION
C-----------------------------------------------------------!
 
C.... Determine the transformation coefficient. Note tha the trans-
C     formation only will occur in the Z-direction (3-direction).
	call transformation_coefficient(trc,50,5,t,y,z,lt)
 
C....Transformation av noderna i ZON S :
	do i=0, inm
	   nnr(i)=0
	enddo
	do j=1, jms
	   do k=1, kms
	      do i=1, ims
	         if ( (nods(i,j,k).gt.0).and.
     &                (nnr(nods(i,j,k)).eq.0)) then
	            call coord_trans(nods,ns1,ns2,ns3,i,j,k,trc,50,5)
	            nnr(nods(i,j,k))=1
		 endif
	      enddo
	   enddo
	enddo
 
C Zon A & Zon B:
	do k=1, kma
	   do j=1, jma
	      do i=1, ima
	       if ((noda(i,j,k).gt.0).and.(nnr(noda(i,j,k)).eq.0)) then
	          call coord_trans(noda,na1,na2,na3,i,j,k,trc,50,5)
	          nnr(noda(i,j,k))=1
	       endif
	      enddo
	   enddo
	   do j=1, jmb
	      do i=1, imb
	       if ((nodb(i,j,k).gt.0).and.(nnr(nodb(i,j,k)).eq.0)) then
	          call coord_trans(nodb,nb1,nb2,nb3,i,j,k,trc,50,5)
	          nnr(nodb(i,j,k))=1
	       endif
	      enddo
	   enddo
	enddo
	nodnum=0
	do i=1, inm
	   if(nnr(i).gt.0) nodnum=nodnum+1
	enddo
	write(*,'(t15,a,i6,a/)')     '=> ', nodnum,
     &                ' number of nodes has been defined !'
 
C.......Eventuella plotfiler skapas om PL=1
	if (pl.eq.1) then
           call create_plotfil(nods,ns1,ns2,ns3, noda,na1,na2,na3,
     &                         nodb,nb1,nb2,nb3, job,jobh)
	endif
 
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	real*8 function xaxel(r,fi,alfa)
C.......Ber{knar X-axelns v{rde d} rho(i,j) , fi(j) & alfa {r k{nda !
	   real*8 r,fi,alfa
	   xaxel=(r+1.d0/r)*alfa*dcosd(fi)
	   return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	real*8 function zaxel(r,fi,alfa)
C.......Ber{knar Z-axelns v{rde d} rho(i,j) , fi(j) & alfa {r k{nda !
	   real*8 r,fi,alfa
	   zaxel=(r-1.d0/r)*alfa*dsind(fi)
	   return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	real*8 function rhox(x,fi,alfa)
C.......Ber{knar rho(i,j) d} Xij, fi(j) & alfa {r k{nda  !
	   real*8 x,fi,alfa
	   rhox=(x+dsqrt(x**2.d0-4.d0*(alfa*dcosd(fi))**2.d0))/
     &          (2.d0*alfa*dcosd(fi))
	   return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	real*8 function rhoz(z,fi,alfa)
C.......Ber{knar rho(i,j) d} Zij, fi(j) & alfa {r k{nda  !
	   real*8 z,fi,alfa
	   rhoz=(z+dsqrt(z**2.d0+4.d0*(alfa*dsind(fi))**2.d0))/
     &          (2.d0*alfa*dsind(fi))
	   return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	real*8 function rho_ellips(a,b,fi,alfa)
	   real*8 fi,a,b,alfa,s1,c1,k
	   c1=dcosd(fi)
	   s1=dsind(fi)
	   k=( 2.d0*( (b*c1)**2.d0-(a*s1)**2.d0 )-(a*b/alfa)**2.d0 )/
     &         ((b*c1)**2.d0+(a*s1)**2.d0)
	   rho_ellips=dsqrt( (-k+dsqrt(k**2.d0-4.d0))/2.d0 )
	   return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	real*8 function calc_beta(my,eta,n)
C----------------------------------------------------------------C
C  Rutinen ber{knar beta. Parametern beta anv{nds vid gradering  C
C  av n{tet i zon S.                                             C
C----------------------------------------------------------------C
	implicit none
	integer  n,i
	real*8   my,eta,lambda,eps,beta,dbeta,sum,sumold
	eps=1.d-06
	lambda=dlog(my)
	beta=0.7d0
	dbeta=0.1d0
	sum=0.d0
10	continue
	   beta=beta+dbeta
	   sumold=sum
	   sum=0
	   do i=1, n
	      sum=sum+dexp(lambda*dble(i-1)**beta)
	   enddo
	   if (dabs(sum/eta-1).gt.eps) then
	      if (((eta-sum)*(eta-sumold)).lt.0) dbeta=-dbeta/2.d0
	      goto 10
	   endif
	calc_beta=beta
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	real*8 function calc_lambda(l1,ln,n)
C----------------------------------------------------------------C
C  Rutinen ber{knar lambda. Parametern lambda anv{nds vid        C
C  gradering av n{tet i olika omr}den .                          C
C----------------------------------------------------------------C
	implicit none
	integer  n,i
	real*8     l1,ln,lambda,eps,eta,dlambda,sum,sumold
	eta=ln/l1
	eps=1.d-06
	lambda=0.015d0
	dlambda=0.005d0
	sum=0
10	continue
	   lambda=lambda+dlambda
	   sumold=sum
	   sum=0
	   do i=1, n
	      sum=sum+dexp(lambda*dble(i-1))
	   enddo
	   if (dabs(sum/eta-1).gt.eps) then
	      if (((eta-sum)*(eta-sumold)).lt.0) dlambda=-dlambda/2.d0
	      goto 10
	   endif
	calc_lambda=lambda
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	real*8 function hyperbel_length(r1,r2,fi,alfa)
C----------------------------------------------------------------C
C  Rutinen ber{knar en str{cka utmed en hyperbel.                C
C   => Hyperbelns ekv (x/a)**2 - (z/b)**2 = 1                    C
C----------------------------------------------------------------C
	real*8 r1,r2,fi,alfa,a,b,dr
	a=alfa*dcosd(fi)
	b=alfa*dsind(fi)
	dr=r2-r1
	sum=dsqrt((a*(1.d0-1.d0/r1**2.d0))**2.d0+
     &             (b*(1.d0+1.d0/r1**2.d0))**2.d0)
C . . . Dela in intervallet i 50 st lika delar. Anvand Simpsons rule.
	do i=2, 50, 2
	   r=r1+dr*dble(i-1)/50.d0
	   sum=sum+4.d0*dsqrt((a*(1.d0-1.d0/r**2.d0))**2.d0+
     &            (b*(1.d0+1.d0/r**2.d0))**2.d0)
	enddo
	do i=3, 49, 2
	   r=r1+dr*dble(i-1)/50.d0
	   sum=sum+2.d0*dsqrt((a*(1.d0-1.d0/r**2.d0))**2.d0+
     &               (b*(1.d0+1.d0/r**2.d0))**2.d0)
	enddo
	sum=sum+dsqrt((a*(1.d0-1.d0/r2**2.d0))**2.d0+
     &           (b*(1.d0+1.d0/r2**2.d0))**2.d0)
	hyperbel_length=sum*dr/(150.d0)
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine transformation_coefficient(trc,n1,n2,t,y,z,lt)
C--- The routine determines the transformation coefficients
	implicit none
	integer  n1,n2,lt,i
	real*8   trc(n1,n2),t,y(0:200),z(0:200,2), a(0:1),b(0:1)
	do i=1, lt
	   a(1) = ( z(i,1)-z(i-1,1) ) / ( y(i)-y(i-1) )
	   a(0) = z(i,1) - a(1)*y(i)
	   b(1) = ( z(i,2)-z(i-1,2) ) / ( y(i)-y(i-1) )
	   b(0) = z(i,2) - b(1)*y(i)
	   trc(i,1) = a(0)
	   trc(i,2) = a(1)
	   trc(i,3) = ( b(0)-a(0) ) / t
	   trc(i,4) = ( b(1)-a(1) ) / t
	   trc(i,5) = y(i) + ( y(i)-y(i-1) ) * 1.d-5
	enddo
	trc(lt+1,5) = 2.d0*y(lt)
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine coord_trans(nod,n1,n2,n3,i,j,k,trc,n4,n5)
C--- The routine performs cordinate transformation
	implicit none
c
        include 'mesh3d_scp_common_nod.f'
c
	integer  n1,n2,n3,n4,n5,nod(n1,n2,n3),i,j,k,m
	real*8   trc(n4,n5)
	m=0
10	m=m+1
	if (npos(nod(i,j,k),2).gt.trc(m,5)) goto 10
	npos(nod(i,j,k),3) = trc(m,1) + trc(m,2)*npos(nod(i,j,k),2) +
     &                trc(m,3)*npos(nod(i,j,k),3) +
     &               trc(m,4)*npos(nod(i,j,k),2)*npos(nod(i,j,k),3)
	return
	end
c
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine reduce_coord_under(nod,n1,n2,n3,
     &                 j1,j2,jst,k1,k2,kst,left1)
C--- The routine creates the nodes in the domain under the crack
C--- front for i =< 5 and 1 =< j =< jma-4.
	implicit none
c
        include 'mesh3d_scp_common_nod.f'
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
C----67--1---------2---------3---------4---------5---------6---------7-!
c
      subroutine reduce_coord_2_to_1(nod,n1,n2,n3,ilocal,i1,i2,ist,
     &                        j1,j2,jst,k,left)
C--- The routine creates nodes in reduced mesh regions.
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
c
      integer  n1,n2,n3,nod(n1,n2,n3),i1,i2,ist,j1,j2,jst,k,ilocal,i,j
      logical  left
c
	if (ilocal.eq.1) then
	   do i=i1, i2, ist
	      if (left) then
	         do j=j1, j2, jst
C . . . . . . Mid- and surface-nodes
C	            NPOS(NOD(I,J,K+2),1) = (NPOS(NOD(I-2,J,K+2),1) +
C     &                     NPOS(NOD(I+2,J,K+2),1) ) / 2.D0
C	            NPOS(NOD(I,J,K+2),2) = (NPOS(NOD(I-2,J,K+2),2) +
C     &                     NPOS(NOD(I+2,J,K+2),2) ) / 2.D0
C	            NPOS(NOD(I,J,K+2),3) = (NPOS(NOD(I-2,J,K+2),3) +
C     &                     NPOS(NOD(I+2,J,K+2),3) ) / 2.D0
 
	            npos(nod(i+1,j,k+1),1) = (npos(nod(i,j,k),1) +
     &                     npos(nod(i+2,j,k+2),1) ) / 2.d0
	            npos(nod(i+1,j,k+1),2) = (npos(nod(i,j,k),2) +
     &                     npos(nod(i+2,j,k+2),2) ) / 2.d0
	            npos(nod(i+1,j,k+1),3) = (npos(nod(i,j,k),3) +
     &                     npos(nod(i+2,j,k+2),3) ) / 2.d0
 
	            npos(nod(i+2,j,k),1) = (npos(nod(i+2,j,k-2),1) +
     &                     npos(nod(i+2,j,k+2),1) ) / 2.d0
	            npos(nod(i+2,j,k),2) = (npos(nod(i+2,j,k-2),2) +
     &                     npos(nod(i+2,j,k+2),2) ) / 2.d0
	            npos(nod(i+2,j,k),3) = (npos(nod(i+2,j,k-2),3) +
     &                     npos(nod(i+2,j,k+2),3) ) / 2.d0
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
C	            NPOS(NOD(I,J,K+2),1) = (NPOS(NOD(I+2,J,K+2),1) +
C     &                     NPOS(NOD(I-2,J,K+2),1) ) / 2.D0
C	            NPOS(NOD(I,J,K+2),2) = (NPOS(NOD(I+2,J,K+2),2) +
C     &                     NPOS(NOD(I-2,J,K+2),2) ) / 2.D0
C	            NPOS(NOD(I,J,K+2),3) = (NPOS(NOD(I+2,J,K+2),3) +
C     &                     NPOS(NOD(I-2,J,K+2),3) ) / 2.D0
 
	            npos(nod(i-1,j,k+1),1) = (npos(nod(i,j,k),1) +
     &                     npos(nod(i-2,j,k+2),1) ) / 2.d0
	            npos(nod(i-1,j,k+1),2) = (npos(nod(i,j,k),2) +
     &                     npos(nod(i-2,j,k+2),2) ) / 2.d0
	            npos(nod(i-1,j,k+1),3) = (npos(nod(i,j,k),3) +
     &                     npos(nod(i-2,j,k+2),3) ) / 2.d0
 
	            npos(nod(i-2,j,k),1) = (npos(nod(i-2,j,k-2),1) +
     &                     npos(nod(i-2,j,k+2),1) ) / 2.d0
	            npos(nod(i-2,j,k),2) = (npos(nod(i-2,j,k-2),2) +
     &                     npos(nod(i-2,j,k+2),2) ) / 2.d0
	            npos(nod(i-2,j,k),3) = (npos(nod(i-2,j,k-2),3) +
     &                     npos(nod(i-2,j,k+2),3) ) / 2.d0
C . . . . . . Surface- and centroid-nodes
	            npos(nod(i,j,k+1),1) = 0.25d0 *
     &               (npos(nod(i+1,j,k),1) + npos(nod(i,j,k+2),1)+
     &                npos(nod(i+2,j,k+1),1) + npos(nod(i-1,j,k+1),1))
	            npos(nod(i,j,k+1),2) = 0.25d0 *
     &               (npos(nod(i+1,j,k),2) + npos(nod(i,j,k+2),2)+
     &                npos(nod(i+2,j,k+1),2) + npos(nod(i-1,j,k+1),2))
	            npos(nod(i,j,k+1),3) = 0.25d0 *
     &               (npos(nod(i+1,j,k),3) + npos(nod(i,j,k+2),3)+
     &                npos(nod(i+2,j,k+1),3) + npos(nod(i-1,j,k+1),3))
 
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
C	            NPOS(NOD(I,J,K+2),1) = (NPOS(NOD(I,J-2,K+2),1) +
C     &                     NPOS(NOD(I,J+2,K+2),1) ) / 2.D0
C	            NPOS(NOD(I,J,K+2),2) = (NPOS(NOD(I,J-2,K+2),2) +
C     &                     NPOS(NOD(I,J+2,K+2),2) ) / 2.D0
C	            NPOS(NOD(I,J,K+2),3) = (NPOS(NOD(I,J-2,K+2),3) +
C     &                     NPOS(NOD(I,J+2,K+2),3) ) / 2.D0
 
	            npos(nod(i,j+1,k+1),1) = (npos(nod(i,j,k),1) +
     &                     npos(nod(i,j+2,k+2),1) ) / 2.d0
	            npos(nod(i,j+1,k+1),2) = (npos(nod(i,j,k),2) +
     &                     npos(nod(i,j+2,k+2),2) ) / 2.d0
	            npos(nod(i,j+1,k+1),3) = (npos(nod(i,j,k),3) +
     &                     npos(nod(i,j+2,k+2),3) ) / 2.d0
 
	            npos(nod(i,j+2,k),1) = (npos(nod(i,j+2,k-2),1) +
     &                     npos(nod(i,j+2,k+2),1) ) / 2.d0
	            npos(nod(i,j+2,k),2) = (npos(nod(i,j+2,k-2),2) +
     &                     npos(nod(i,j+2,k+2),2) ) / 2.d0
	            npos(nod(i,j+2,k),3) = (npos(nod(i,j+2,k-2),3) +
     &                     npos(nod(i,j+2,k+2),3) ) / 2.d0
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
C	            NPOS(NOD(I,J,K+2),1) = (NPOS(NOD(I,J+2,K+2),1) +
C     &                     NPOS(NOD(I,J-2,K+2),1) ) / 2.D0
C	            NPOS(NOD(I,J,K+2),2) = (NPOS(NOD(I,J+2,K+2),2) +
C     &                     NPOS(NOD(I,J-2,K+2),2) ) / 2.D0
C	            NPOS(NOD(I,J,K+2),3) = (NPOS(NOD(I,J+2,K+2),3) +
C     &                     NPOS(NOD(I,J-2,K+2),3) ) / 2.D0
 
	            npos(nod(i,j-1,k+1),1) = (npos(nod(i,j,k),1) +
     &                     npos(nod(i,j-2,k+2),1) ) / 2.d0
	            npos(nod(i,j-1,k+1),2) = (npos(nod(i,j,k),2) +
     &                     npos(nod(i,j-2,k+2),2) ) / 2.d0
	            npos(nod(i,j-1,k+1),3) = (npos(nod(i,j,k),3) +
     &                     npos(nod(i,j-2,k+2),3) ) / 2.d0
 
	            npos(nod(i,j-2,k),1) = (npos(nod(i,j-2,k-2),1) +
     &                     npos(nod(i,j-2,k+2),1) ) / 2.d0
	            npos(nod(i,j-2,k),2) = (npos(nod(i,j-2,k-2),2) +
     &                     npos(nod(i,j-2,k+2),2) ) / 2.d0
	            npos(nod(i,j-2,k),3) = (npos(nod(i,j-2,k-2),3) +
     &                     npos(nod(i,j-2,k+2),3) ) / 2.d0
C . . . . . . Surface- and centroid-nodes
	            npos(nod(i,j,k+1),1) = 0.25d0 *
     &               (npos(nod(i,j+1,k),1) + npos(nod(i,j,k+2),1)+
     &                npos(nod(i,j+2,k+1),1) + npos(nod(i,j-1,k+1),1))
	            npos(nod(i,j,k+1),2) = 0.25d0 *
     &               (npos(nod(i,j+1,k),2) + npos(nod(i,j,k+2),2)+
     &                npos(nod(i,j+2,k+1),2) + npos(nod(i,j-1,k+1),2))
	            npos(nod(i,j,k+1),3) = 0.25d0 *
     &               (npos(nod(i,j+1,k),3) + npos(nod(i,j,k+2),3)+
     &                npos(nod(i,j+2,k+1),3) + npos(nod(i,j-1,k+1),3))
 
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
C----67--1---------2---------3---------4---------5---------6---------7-!
c
	subroutine reduce_coord_3_to_1(nod,n1,n2,n3,j1,j2,jst,k,left)
C--- The routine creates nodes in reduced mesh region under the crack
C--- front.
      implicit none
c
      include 'mesh3d_scp_common_nod.f'
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
C----67--1---------2---------3---------4---------5---------6---------7-!
c
