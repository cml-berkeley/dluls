c**********************************************************************
c     Subroutines in this file :
c     1. reyneq
c     2. acmdn
c     3. acmup
c     4. flow
c     5. create_dbase
c     6. calc_bearing_number
c**********************************************************************

c======================================================================
	SUBROUTINE reyneq(idirect)
c======================================================================

c----------------------------------------------------------------------
c     This is the reynolds equation solver
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE siml_pram
      USE syst_cnfg
      USE sldr_grid
      USE sldr_arys
      USE rynl_data
      USE rynl_arys
      USE sldr_dynm
      USE aspr_data

	IMPLICIT REAL*8(a-h,o-z)
      
	nx1=2+(nx-2)/2
	ny1=2+(ny-2)/2
	nx2=2+(nx1-2)/2
	ny2=2+(ny1-2)/2
	nx3=2+(nx2-2)/2
	ny3=2+(ny2-2)/2
	nx4=2+(nx3-2)/2
	ny4=2+(ny3-2)/2

c	Calculate:
c	 1) the flow rarefaction factors qn
c	 2) the average flow factor of rough bearing surfaces
c	=====================================================
	nreyneq=0
 5	nreyneq=nreyneq+1

	DO 300 j=2,nym1
		jm1=j-1
		jp1=j+1
		delyjp1=yref(jp1)-yref(j  )
		delyjm1=yref(j	)-yref(jm1)
		IF(j.NE.nym1) delyjp2=yref(j+2)-yref(jp1)
		IF(j.NE.2) delyjm2=yref(jm1)-yref(j-2)
		ymyp   =(yref(jp1)-yref(jm1))/2.d0

		DO 250 i=2,nxm1
			im1=i-1
			ip1=i+1
			delxip1=xref(ip1)-xref(i	)
			delxim1=xref(i  )-xref(im1)
			IF(i.NE.nxm1) delxip2=xref(i+2)-xref(ip1)
			IF(i.NE.2) delxim2=xref(im1)-xref(i-2)
			xmxp	 =(xref(ip1)-xref(im1))/2.d0

			h2im = hydnew(i  ,j  )
			h2ip = hydnew(i+1,j  )
			h2jm = hxlnew(i  ,j  )
			h2jp = hxlnew(i  ,j+1)

			h2im_b = hyunew(i  ,j  )
			h2ip_b = hyunew(i+1,j  )
			h2jm_b = hxrnew(i  ,j  )
			h2jp_b = hxrnew(i  ,j+1)

			ztam=zta(i  ,j  )
			ztap=zta(i  ,j+1)
			etam=eta(i  ,j  )
			etap=eta(i+1,j  )

			p2im = (p(im1,j)+p(i,j))/2.d0
			p2ip = (p(ip1,j)+p(i,j))/2.d0
			p2jm = (p(i,jm1)+p(i,j))/2.d0
			p2jp = (p(i,jp1)+p(i,j))/2.d0

			pn=p2im*h2im
			CALL flow(pn,qnim)
			qnh2im = qnim*h2im*h2im
			fmi = h2im*(bearx(im1,j)+bearx(i,j))/2.d0
			dmi=qnh2im/delxim1
			pemi=fmi/dmi

			IF(etam.EQ.1.d0) THEN
				ami_b=0.d0
				fmi_b=0.d0
			ELSE
				pn_b=p2im*h2im_b
				CALL flow(pn_b,qnim_b)
				qnh2im_b = qnim_b*h2im_b*h2im_b
				fmi_b = h2im_b*(bearx(im1,j)+bearx(i,j))/2.d0
				dmi_b=qnh2im_b/delxim1
				pemi_b=fmi_b/dmi_b
			ENDIF

			pn=p2ip*h2ip
			CALL flow(pn,qnip)
			qnh2ip = qnip*h2ip*h2ip
			fpi = h2ip*(bearx(ip1,j)+bearx(i,j))/2.d0
			dpi=qnh2ip/delxip1
			pepi=fpi/dpi

			IF(etap.EQ.1.d0) THEN
				api_b=0.d0
				fpi_b=0.d0
			ELSE
				pn_b=p2ip*h2ip_b
				CALL flow(pn_b,qnip_b)
				qnh2ip_b = qnip_b*h2ip_b*h2ip_b
				fpi_b = h2ip_b*(bearx(ip1,j)+bearx(i,j))/2.d0
				dpi_b=qnh2ip_b/delxip1
				pepi_b=fpi_b/dpi_b
			ENDIF

			pn=p2jm*h2jm
			CALL flow(pn,qnjm)
			qnh2jm = qnjm*h2jm*h2jm
			fmj = h2jm*(beary(i,jm1)+beary(i,j))/2.d0
			dmj=qnh2jm/delyjm1
			pemj=fmj/dmj

			IF(ztam.EQ.1.d0) THEN
				amj_b=0.d0
				fmj_b=0.d0
			ELSE
				pn_b=p2jm*h2jm_b
				CALL flow(pn_b,qnjm_b)
				qnh2jm_b = qnjm_b*h2jm_b*h2jm_b
				fmj_b = h2jm_b*(beary(i,jm1)+beary(i,j))/2.d0
				dmj_b=qnh2jm_b/delyjm1
				pemj_b=fmj_b/dmj_b
			ENDIF

			pn=p2jp*h2jp
			CALL flow(pn,qnjp)
			qnh2jp = qnjp*h2jp*h2jp
			fpj = h2jp*(beary(i,jp1)+beary(i,j))/2.d0
			dpj=qnh2jp/delyjp1
			pepj=fpj/dpj

			IF(ztap.EQ.1.d0) THEN
				apj_b=0.d0
				fpj_b=0.d0
			ELSE
				pn_b=p2jp*h2jp_b
				CALL flow(pn_b,qnjp_b)
				qnh2jp_b = qnjp_b*h2jp_b*h2jp_b
				fpj_b = h2jp_b*(beary(i,jp1)+beary(i,j))/2.d0
				dpj_b=qnh2jp_b/delyjp1
				pepj_b=fpj_b/dpj_b
			ENDIF

			GOTO (11,12,13,14,15),idisc

c			Power law
c			=========
 11			ami=DMAX1(0.d0,dmi*(1.d0-0.1d0*DABS(pemi))**5)
			api=DMAX1(0.d0,dpi*(1.d0-0.1d0*DABS(pepi))**5)
			amj=DMAX1(0.d0,dmj*(1.d0-0.1d0*DABS(pemj))**5)
			apj=DMAX1(0.d0,dpj*(1.d0-0.1d0*DABS(pepj))**5)

			IF(etam.NE.1.d0) THEN
				ami_b=DMAX1(0.d0,dmi_b*(1.d0-0.1d0*DABS(pemi_b))**5)
			ENDIF
			IF(etap.NE.1.d0) THEN
				api_b=DMAX1(0.d0,dpi_b*(1.d0-0.1d0*DABS(pepi_b))**5)
			ENDIF
			IF(ztam.NE.1.d0) THEN
				amj_b=DMAX1(0.d0,dmj_b*(1.d0-0.1d0*DABS(pemj_b))**5)
			ENDIF
			IF(ztap.NE.1.d0) THEN
				apj_b=DMAX1(0.d0,dpj_b*(1.d0-0.1d0*DABS(pepj_b))**5)
			ENDIF
			GOTO 20

c			Central difference
c			==================
 12			ami=dmi*(1.d0-0.5d0*DABS(pemi))
			api=dpi*(1.d0-0.5d0*DABS(pepi))
			amj=dmj*(1.d0-0.5d0*DABS(pemj))
			apj=dpj*(1.d0-0.5d0*DABS(pepj))

			IF(etam.NE.1.d0) THEN
				ami_b=dmi_b*(1.d0-0.5d0*DABS(pemi_b))
			ENDIF
			IF(etap.NE.1.d0) THEN
				api_b=dpi_b*(1.d0-0.5d0*DABS(pepi_b))
			ENDIF
			IF(ztam.NE.1.d0) THEN
				amj_b=dmj_b*(1.d0-0.5d0*DABS(pemj_b))
			ENDIF
			IF(ztap.NE.1.d0) THEN
				apj_b=dpj_b*(1.d0-0.5d0*DABS(pepj_b))
			ENDIF

			GOTO 20

c			1st-order upwind: idisc = 2
c			===========================
 13			ami=dmi
			api=dpi
			amj=dmj
			apj=dpj

			IF(etam.NE.1.d0) THEN
				ami_b=dmi_b
			ENDIF
			IF(etap.NE.1.d0) THEN
				api_b=dpi_b
			ENDIF
			IF(ztam.NE.1.d0) THEN
				amj_b=dmj_b
			ENDIF
			IF(ztap.NE.1.d0) THEN
				apj_b=dpj_b
			ENDIF
			GOTO 20

c			Hybrid
c			======
 14			ami=DMAX1(0.d0,dmi*(1.0d0-0.5d0*DABS(pemi)))
			api=DMAX1(0.d0,dpi*(1.0d0-0.5d0*DABS(pepi)))
			amj=DMAX1(0.d0,dmj*(1.0d0-0.5d0*DABS(pemj)))
			apj=DMAX1(0.d0,dpj*(1.0d0-0.5d0*DABS(pepj)))

			IF(etam.NE.1.d0) THEN
				ami_b=DMAX1(0.d0,dmi_b*(1.0d0-0.5d0*DABS(pemi_b)))
			ENDIF
			IF(etap.NE.1.d0) THEN
				api_b=DMAX1(0.d0,dpi_b*(1.0d0-0.5d0*DABS(pepi_b)))
			ENDIF
			IF(ztam.NE.1.d0) THEN
				amj_b=DMAX1(0.d0,dmj_b*(1.0d0-0.5d0*DABS(pemj_b)))
			ENDIF
			IF(ztap.NE.1.d0) THEN
				apj_b=DMAX1(0.d0,dpj_b*(1.0d0-0.5d0*DABS(pepj_b)))
			ENDIF
			GOTO 20

c			Exponential
c			===========
 15			ami=dmi*DABS(pemi)/(dexp(DABS(pemi))-1.d0)
			api=dpi*DABS(pepi)/(dexp(DABS(pepi))-1.d0)
			amj=dmj*DABS(pemj)/(dexp(DABS(pemj))-1.d0)
			apj=dpj*DABS(pepj)/(dexp(DABS(pepj))-1.d0)

			IF(etam.NE.1.d0) THEN
				ami_b=dmi_b*DABS(pemi_b)/(dexp(DABS(pemi_b))-1.d0)
			ENDIF
			IF(etap.NE.1.d0) THEN
				api_b=dpi_b*DABS(pepi_b)/(dexp(DABS(pepi_b))-1.d0)
			ENDIF
			IF(ztam.NE.1.d0) THEN
				amj_b=dmj_b*DABS(pemj_b)/(dexp(DABS(pemj_b))-1.d0)
			ENDIF
			IF(ztap.NE.1.d0) THEN
				apj_b=dpj_b*DABS(pepj_b)/(dexp(DABS(pepj_b))-1.d0)
			ENDIF

 20			CONTINUE

			aw(im1,jm1)=-ymyp*(etam*(ami+DMAX1( fmi,0.d0))+
     &			   (1.d0-etam)*(ami_b+DMAX1( fmi_b,0.d0)))
			ae(im1,jm1)=-ymyp*(etap*(api+DMAX1(-fpi,0.d0))+
     &			   (1.d0-etap)*(api_b+DMAX1(-fpi_b,0.d0)))
			as(im1,jm1)=-xmxp*(ztam*(amj+DMAX1( fmj,0.d0))+
     &			   (1.d0-ztam)*(amj_b+DMAX1( fmj_b,0.d0)))
			an(im1,jm1)=-xmxp*(ztap*(apj+DMAX1(-fpj,0.d0))+
     &			   (1.d0-ztap)*(apj_b+DMAX1(-fpj_b,0.d0)))
			ap(im1,jm1)=-aw(im1,jm1)-ae(im1,jm1)
     &				  -as(im1,jm1)-an(im1,jm1)
     &			  -ymyp*(etam*fmi+(1.d0-etam)*fmi_b)
     &			  +ymyp*(etap*fpi+(1.d0-etap)*fpi_b)
     &			  -xmxp*(ztam*fmj+(1.d0-ztam)*fmj_b)
     &			  +xmxp*(ztap*fpj+(1.d0-ztap)*fpj_b)

			IF(idirect.NE.1) THEN
				amult=xmxp*ymyp*si/dt/omega0
				ar(im1,jm1)=amult*h(i,j)*pold(i,j)
				ap(im1,jm1)=ap(im1,jm1)+amult*hnew(i,j)
			ELSE
				ar(im1,jm1)=0.d0
			ENDIF

250		CONTINUE
300	CONTINUE

	akd=0.d0
	akn=0.d0
	DO j=2,nym1
		jm1=j-1
		jp1=j+1
		DO i=2,nxm1
			im1=i-1
			ip1=i+1
			aj(im1)=aw(im1,jm1)
			bj(im1)=ap(im1,jm1)
			cj(im1)=ae(im1,jm1)
			dj(im1)=ar(im1,jm1)-as(im1,jm1)*p(i,jm1)
     &			  -an(im1,jm1)*p(i,jp1)
			IF(i.EQ.2) THEN
				dj(im1)=dj(im1)-aw(im1,jm1)*p(im1,j)
				aj(im1)=0.d0
			ELSEIF (i.EQ.nxm1) THEN
				dj(im1)=dj(im1)-ae(im1,jm1)*p(ip1,j)
				cj(im1)=0.d0
			ENDIF
			akd=akd+DABS(aj(im1)*p(im1,j)+bj(im1)*p(i,j)+
     &			cj(im1)*p(ip1,j)-dj(im1))
			akn=akn+DABS(bj(im1)*p(i,j))
		ENDDO
	ENDDO

	ak=akd/akn
	IF((mod(nreyneq,5).EQ.0)) THEN
c		WRITE(*,990)t,nreyneq,ak
990		FORMAT('	T(s)=',E16.9,',Iteration=',I3,',Residual=',E16.9)
	ENDIF

c     Convergence
c     ===========
	IF(nreyneq.GT.1.AND.((ak.LT.akmax).OR.(nreyneq.GT.nitmax)))RETURN

	rakpre=ak
	kint=0
 50	kint=kint+1

c	First sweep in x-direction
c	==========================
	DO 500 j=2,nym1
		jm1=j-1
		jp1=j+1
		ymyp=(yref(jp1)-yref(jm1))/2.d0

		DO 450 i=2,nxm1
			im1=i-1
			ip1=i+1
			xmxp=(xref(ip1)-xref(im1))/2.d0

			aj(im1)=aw(im1,jm1)
			bj(im1)=ap(im1,jm1)
			cj(im1)=ae(im1,jm1)
			dj(im1)=ar(im1,jm1)-as(im1,jm1)*p(i,jm1)
     &				-an(im1,jm1)*p(i,jp1)

			IF(i.EQ.2) THEN
				dj(im1)=dj(im1)-aw(im1,jm1)*p(im1,j)
				aj(im1)=0.d0
			ELSEIF (i.EQ.nxm1) THEN
				dj(im1)=dj(im1)-ae(im1,jm1)*p(ip1,j)
				cj(im1)=0.d0
			ENDIF

450		CONTINUE
		CALL tridag(nxm2,aj,bj,cj,dj)
		  
		DO 475 i=1,nxm2
			p(i+1,j)=dj(i)
475		CONTINUE

500	CONTINUE

c	Block correction along x direction
c	==================================
	DO i=2,nxm1
		im1=i-1
		ip1=i+1
		xmxp=(xref(ip1)-xref(im1))/2.d0

		aj(im1)=0.d0
		bj(im1)=0.d0
		cj(im1)=0.d0
		dj(im1)=0.d0

		DO j=2,nym1
			jm1=j-1
			jp1=j+1
			ymyp=(yref(jp1)-yref(jm1))/2.d0

			aj(im1)=aj(im1)+aw(im1,jm1)
c			bj(im1)=bj(im1)+ap(im1,jm1)+as(im1,jm1)+an(im1,jm1)
			bj(im1)=bj(im1)+ap(im1,jm1)
			IF(j.NE.2) bj(im1)=bj(im1)+as(im1,jm1)
			IF(j.NE.nym1) bj(im1)=bj(im1)+an(im1,jm1)
			cj(im1)=cj(im1)+ae(im1,jm1)
			dj(im1)=dj(im1)+ar(im1,jm1)-aw(im1,jm1)*p(im1,j)
     &			  -ae(im1,jm1)*p(ip1,j)-as(im1,jm1)*p(i,jm1)
     &			  -an(im1,jm1)*p(i,jp1)-ap(im1,jm1)*p(i  ,j)

		ENDDO
	ENDDO
	  
	aj(1)=0.d0
	cj(nxm2)=0.d0
	
	CALL tridag(nxm2,aj,bj,cj,dj)
	DO i=1,nxm2
		DO j=1,nym2
			p(i+1,j+1)=p(i+1,j+1)+dj(i)
		ENDDO
	ENDDO

c	Then sweep in y-direction
c	=========================
	DO 700 i=2,nxm1
		im1=i-1
		ip1=i+1
		xmxp=(xref(ip1)-xref(im1))/2.d0

		DO 650 j=2,nym1
			jm1=j-1
			jp1=j+1
			ymyp=(yref(jp1)-yref(jm1))/2.d0

			aj(jm1)=as(im1,jm1)
			bj(jm1)=ap(im1,jm1)
			cj(jm1)=an(im1,jm1)
			dj(jm1)=ar(im1,jm1)-aw(im1,jm1)*p(im1,j)
     &				-ae(im1,jm1)*p(ip1,j)

			IF(j.EQ.2) THEN
				dj(jm1)=dj(jm1)-as(im1,jm1)*p(i,jm1)
				aj(jm1)=0.d0
			ELSEIF(j.EQ.nym1) THEN
				dj(jm1)=dj(jm1)-an(im1,jm1)*p(i,jp1)
				cj(jm1)=0.d0
			ENDIF

650		CONTINUE
		CALL tridag(nym2,aj,bj,cj,dj)
		  
		DO 675 j=1,nym2
			p(i,j+1)=dj(j)
675		CONTINUE
700	CONTINUE

c	Block correction along y direction
c	==================================
	DO j=2,nym1
		jm1=j-1
		jp1=j+1
		ymyp=(yref(jp1)-yref(jm1))/2.d0

		aj(jm1)=0.d0
		bj(jm1)=0.d0
		cj(jm1)=0.d0
		dj(jm1)=0.d0

		DO i=2,nxm1
			im1=i-1
			ip1=i+1
			xmxp=(xref(ip1)-xref(im1))/2.d0

			aj(jm1)=aj(jm1)+as(im1,jm1)
c			bj(jm1)=bj(jm1)+ap(im1,jm1)+aw(im1,jm1)+ae(im1,jm1)
			bj(jm1)=bj(jm1)+ap(im1,jm1)
			IF(i.NE.2) bj(jm1)=bj(jm1)+aw(im1,jm1)
			IF(i.NE.nxm1) bj(jm1)=bj(jm1)+ae(im1,jm1)
			cj(jm1)=cj(jm1)+an(im1,jm1)
			dj(jm1)=dj(jm1)+ar(im1,jm1)-aw(im1,jm1)*p(im1,j)
     &			  -ae(im1,jm1)*p(ip1,j)-as(im1,jm1)*p(i,jm1)
     &			  -an(im1,jm1)*p(i,jp1)-ap(im1,jm1)*p(i  ,j)

		ENDDO
	ENDDO
	  
	aj(1)=0.d0
	cj(nym2)=0.d0
	CALL tridag(nym2,aj,bj,cj,dj)
	DO j=1,nym2
		DO i=1,nxm2
			p(i+1,j+1)=p(i+1,j+1)+dj(j)
		ENDDO
	ENDDO

	akd=0.d0
	akn=0.d0
	DO j=2,nym1
		jm1=j-1
		jp1=j+1
		DO i=2,nxm1
			im1=i-1
			ip1=i+1
			aj(im1)=aw(im1,jm1)
			bj(im1)=ap(im1,jm1)
			cj(im1)=ae(im1,jm1)
			dj(im1)=ar(im1,jm1)-as(im1,jm1)*p(i,jm1)
     &			-an(im1,jm1)*p(i,jp1)
			IF(i.EQ.2) THEN
				dj(im1)=dj(im1)-aw(im1,jm1)*p(im1,j)
				aj(im1)=0.d0
			ELSEIF (i.EQ.nxm1) THEN
				dj(im1)=dj(im1)-ae(im1,jm1)*p(ip1,j)
				cj(im1)=0.d0
			ENDIF
			akd=akd+DABS(aj(im1)*p(im1,j)+bj(im1)*p(i,j)+
     &			cj(im1)*p(ip1,j)-dj(im1))
			akn=akn+DABS(bj(im1)*p(i,j))
		ENDDO
	ENDDO

	rak=akd/akn
	ratio=rak/rakpre
	rakpre=rak

	IF(ratio.LE.0.5d0) GOTO 50

c	Multigrid iteration
c	===================

c	Move down to level 1
c	====================
	CALL acmdn(aw,ae,as,an,ap,ar,p,aw1,ae1,as1,an1,ap1,ar1,dp1,
     &	nx,ny,nx1,ny1,nxmax,nymax,nxmax1,nymax1,1)

c	Move down to level 2
c	====================
	CALL acmdn(aw1,ae1,as1,an1,ap1,ar1,dp1,aw2,ae2,as2,an2,ap2,
     &	ar2,dp2,nx1,ny1,nx2,ny2,nxmax1,nymax1,nxmax2,nymax2,2)

c	Move down to level 3
c	====================
	CALL acmdn(aw2,ae2,as2,an2,ap2,ar2,dp2,aw3,ae3,as3,an3,ap3,
     &	ar3,dp3,nx2,ny2,nx3,ny3,nxmax2,nymax2,nxmax3,nymax3,4)

c	Move down to level 4
c	====================
	CALL acmdn(aw3,ae3,as3,an3,ap3,ar3,dp3,aw4,ae4,as4,an4,ap4,
     &	ar4,dp4,nx3,ny3,nx4,ny4,nxmax3,nymax3,nxmax4,nymax4,16)

c	Move up to level 3
c	====================
	CALL acmup(aw3,ae3,as3,an3,ap3,ar3,dp3,dp4,nx3,ny3,
     &	nx4,ny4,nxmax3,nymax3,nxmax4,nymax4,16)

c	Move up to level 2
c	==================
	CALL acmup(aw2,ae2,as2,an2,ap2,ar2,dp2,dp3,nx2,ny2,
     &	nx3,ny3,nxmax2,nymax2,nxmax3,nymax3,16)

c	Move up to level 1
c	==================
	CALL acmup(aw1,ae1,as1,an1,ap1,ar1,dp1,dp2,nx1,ny1,
     &	nx2,ny2,nxmax1,nymax1,nxmax2,nymax2,2)

c	Move up to level 0
c	==================
	nx1m1=nx1-1
	ny1m1=ny1-1

	DO k=2,nx1m1
		DO l=2,ny1m1
			i1=2*(k-1)
			i2=2*(k-1)+1
			j1=2*(l-1)
			j2=2*(l-1)+1
			p(i1,j1)=p(i1,j1)+dp1(k,l)
			p(i2,j1)=p(i2,j1)+dp1(k,l)
			p(i1,j2)=p(i1,j2)+dp1(k,l)
			p(i2,j2)=p(i2,j2)+dp1(k,l)
		ENDDO
	ENDDO

c	First sweep in x-direction
c	==========================
	DO j=2,nym1
		jm1=j-1
		jp1=j+1
		ymyp=(yref(jp1)-yref(jm1))/2.d0

		DO i=2,nxm1
			im1=i-1
			ip1=i+1
			xmxp=(xref(ip1)-xref(im1))/2.d0

			aj(im1)=aw(im1,jm1)
			bj(im1)=ap(im1,jm1)
			cj(im1)=ae(im1,jm1)
			dj(im1)=ar(im1,jm1)-as(im1,jm1)*p(i,jm1)
     &			-an(im1,jm1)*p(i,jp1)

			IF(i.EQ.2) THEN
				dj(im1)=dj(im1)-aw(im1,jm1)*p(im1,j)
				aj(im1)=0.d0
			ELSEIF (i.EQ.nxm1) THEN
				dj(im1)=dj(im1)-ae(im1,jm1)*p(ip1,j)
				cj(im1)=0.d0
			ENDIF

		ENDDO
		CALL tridag(nxm2,aj,bj,cj,dj)
		
		DO i=1,nxm2
			p(i+1,j)=dj(i)
		ENDDO
	ENDDO

c	Block correction along x direction
c	==================================
	DO i=2,nxm1
		im1=i-1
		ip1=i+1
		xmxp=(xref(ip1)-xref(im1))/2.d0

		aj(im1)=0.d0
		bj(im1)=0.d0
		cj(im1)=0.d0
		dj(im1)=0.d0

		DO j=2,nym1
			jm1=j-1
			jp1=j+1
			ymyp=(yref(jp1)-yref(jm1))/2.d0

			aj(im1)=aj(im1)+aw(im1,jm1)
c			bj(im1)=bj(im1)+ap(im1,jm1)+as(im1,jm1)+an(im1,jm1)
			bj(im1)=bj(im1)+ap(im1,jm1)
		  
			IF(j.NE.2) bj(im1)=bj(im1)+as(im1,jm1)
			IF(j.NE.nym1) bj(im1)=bj(im1)+an(im1,jm1)
			cj(im1)=cj(im1)+ae(im1,jm1)
			dj(im1)=dj(im1)+ar(im1,jm1)-aw(im1,jm1)*p(im1,j)
     &			-ae(im1,jm1)*p(ip1,j)-as(im1,jm1)*p(i,jm1)
     &			-an(im1,jm1)*p(i,jp1)-ap(im1,jm1)*p(i  ,j)

		ENDDO
	ENDDO
	
	aj(1)=0.d0
	cj(nxm2)=0.d0
	CALL tridag(nxm2,aj,bj,cj,dj)
	
	DO i=1,nxm2
		DO j=1,nym2
			p(i+1,j+1)=p(i+1,j+1)+dj(i)
		ENDDO
	ENDDO

c	Then sweep in y-direction
c	=========================
	DO i=2,nxm1
		im1=i-1
		ip1=i+1
		xmxp=(xref(ip1)-xref(im1))/2.d0

		DO j=2,nym1
			jm1=j-1
			jp1=j+1
			ymyp=(yref(jp1)-yref(jm1))/2.d0

			aj(jm1)=as(im1,jm1)
			bj(jm1)=ap(im1,jm1)
			cj(jm1)=an(im1,jm1)
			dj(jm1)=ar(im1,jm1)-aw(im1,jm1)*p(im1,j)
     &			-ae(im1,jm1)*p(ip1,j)

			IF(j.EQ.2) THEN
				dj(jm1)=dj(jm1)-as(im1,jm1)*p(i,jm1)
				aj(jm1)=0.d0
			ELSEIF(j.EQ.nym1) THEN
				dj(jm1)=dj(jm1)-an(im1,jm1)*p(i,jp1)
				cj(jm1)=0.d0
			ENDIF
		ENDDO

		CALL tridag(nym2,aj,bj,cj,dj)
		DO j=1,nym2
		  p(i,j+1)=dj(j)
		ENDDO
	ENDDO

c	Block correction along y direction
c	==================================
	DO j=2,nym1
		jm1=j-1
		jp1=j+1
		ymyp=(yref(jp1)-yref(jm1))/2.d0

		aj(jm1)=0.d0
		bj(jm1)=0.d0
		cj(jm1)=0.d0
		dj(jm1)=0.d0

		DO i=2,nxm1
			im1=i-1
			ip1=i+1
			xmxp=(xref(ip1)-xref(im1))/2.d0

			aj(jm1)=aj(jm1)+as(im1,jm1)
c			bj(jm1)=bj(jm1)+ap(im1,jm1)+aw(im1,jm1)+ae(im1,jm1)
			bj(jm1)=bj(jm1)+ap(im1,jm1)
			IF(i.NE.2) bj(jm1)=bj(jm1)+aw(im1,jm1)
			IF(i.NE.nxm1) bj(jm1)=bj(jm1)+ae(im1,jm1)
			cj(jm1)=cj(jm1)+an(im1,jm1)
			dj(jm1)=dj(jm1)+ar(im1,jm1)-aw(im1,jm1)*p(im1,j)
     &			-ae(im1,jm1)*p(ip1,j)-as(im1,jm1)*p(i,jm1)
     &			-an(im1,jm1)*p(i,jp1)-ap(im1,jm1)*p(i  ,j)

		ENDDO
	ENDDO
	
	aj(1)=0.d0
	cj(nym2)=0.d0
	CALL tridag(nym2,aj,bj,cj,dj)
	DO j=1,nym2
		DO i=1,nxm2
			p(i+1,j+1)=p(i+1,j+1)+dj(j)
			IF(p(i+1,j+1).LT.0.d0) THEN
				p(i+1,j+1)=5.d-2
			ENDIF
		ENDDO
	ENDDO

	  akd=0.d0
	  akn=0.d0
	  DO j=2,nym1
		jm1=j-1
		jp1=j+1
		DO i=2,nxm1
			im1=i-1
			ip1=i+1
			aj(im1)=aw(im1,jm1)
			bj(im1)=ap(im1,jm1)
			cj(im1)=ae(im1,jm1)
			dj(im1)=ar(im1,jm1)-as(im1,jm1)*p(i,jm1)
     &			-an(im1,jm1)*p(i,jp1)
			IF(i.EQ.2) THEN
				dj(im1)=dj(im1)-aw(im1,jm1)*p(im1,j)
				aj(im1)=0.d0
			ELSEIF (i.EQ.nxm1) THEN
				dj(im1)=dj(im1)-ae(im1,jm1)*p(ip1,j)
				cj(im1)=0.d0
			ENDIF
			akd=akd+DABS(aj(im1)*p(im1,j)+bj(im1)*p(i,j)+
     &			cj(im1)*p(ip1,j)-dj(im1))
			akn=akn+DABS(bj(im1)*p(i,j))
		ENDDO
	ENDDO

	rak=akd/akn
	rakpre=rak

	IF(rak.GT.akmax.AND.kint.LT.nitmax) GOTO 50
	GOTO 5

	END

c======================================================================
	SUBROUTINE acmdn(aw1,ae1,as1,an1,ap1,ar1,p1,aw2,ae2,as2,
     &	an2,ap2,ar2,dp2,nx1,ny1,nx2,ny2,nx1max,ny1max,nx2max,
     &	ny2max,nlimit)
c======================================================================

	IMPLICIT REAL*8(a-h,o-z)
	
c     Arguments
c     =========	
	DOUBLE PRECISION  aw1(nx1max,ny1max),ae1(nx1max,ny1max),
     &			        as1(nx1max,ny1max),an1(nx1max,ny1max),
     &			        ap1(nx1max,ny1max),ar1(nx1max,ny1max),
     &			        aw2(nx2max,ny2max),ae2(nx2max,ny2max),
     &			        as2(nx2max,ny2max),an2(nx2max,ny2max),
     &			        ap2(nx2max,ny2max),ar2(nx2max,ny2max),
     &			        p1(nx1max,ny1max), dp2(nx2max,ny2max)

c     Local Variables
c     ===============     
      DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: aj,bj,cj,dj

	nx2m1=nx2-1
	ny2m1=ny2-1
	nx2m2=nx2-2
	ny2m2=ny2-2
	
	ALLOCATE(aj(nx2m2))
      ALLOCATE(bj(nx2m2))
      ALLOCATE(cj(nx2m2))
      ALLOCATE(dj(nx2m2))
	
	aj = 0.0d0
	bj = 0.0d0
	cj = 0.0d0
	dj = 0.0d0	

	DO k=2,nx2m1
		k1=k-1
		DO l=2,ny2m1
			l1=l-1
			i1=2*(k-1)-1
			i2=2*(k-1)+1-1
			j1=2*(l-1)-1
			j2=2*(l-1)+1-1
			ap2(k1,l1)=ap1(i1,j1)+ap1(i2,j1)+ap1(i1,j2)+ap1(i2,j2)+
     &		ae1(i1,j1)+ae1(i1,j2)+aw1(i2,j1)+aw1(i2,j2)+
     &		as1(i1,j2)+as1(i2,j2)+an1(i1,j1)+an1(i2,j1)
			aw2(k1,l1)=aw1(i1,j1)+aw1(i1,j2)
			ae2(k1,l1)=ae1(i2,j1)+ae1(i2,j2)
			as2(k1,l1)=as1(i1,j1)+as1(i2,j1)
			an2(k1,l1)=an1(i1,j2)+an1(i2,j2)
			ar2(k1,l1)=ar1(i1,j1)+ar1(i2,j1)+ar1(i1,j2)+ar1(i2,j2)-
     &		ap1(i1,j1)*p1(i1+1,j1+1)-ap1(i2,j1)*p1(i2+1,j1+1)-
     &		ap1(i1,j2)*p1(i1+1,j2+1)-ap1(i2,j2)*p1(i2+1,j2+1)-
     &		aw1(i1,j1)*p1(i1  ,j1+1)-aw1(i1,j2)*p1(i1	,j2+1)-
     &		ae1(i1,j1)*p1(i2+1,j1+1)-ae1(i1,j2)*p1(i2+1,j2+1)-
     &		as1(i1,j1)*p1(i1+1,j1	)-as1(i1,j2)*p1(i1+1,j1+1)-
     &		an1(i1,j1)*p1(i1+1,j2+1)-an1(i1,j2)*p1(i1+1,j2+2)-
     &		aw1(i2,j1)*p1(i1+1,j1+1)-aw1(i2,j2)*p1(i1+1,j2+1)-
     &		ae1(i2,j1)*p1(i2+2,j1+1)-ae1(i2,j2)*p1(i2+2,j2+1)-
     &		as1(i2,j1)*p1(i2+1,j1	)-as1(i2,j2)*p1(i2+1,j1+1)-
     &		an1(i2,j1)*p1(i2+1,j2+1)-an1(i2,j2)*p1(i2+1,j2+2)
		ENDDO
	ENDDO

	DO k=1,nx2
		DO l=1,ny2
			dp2(k,l)=0.d0
		ENDDO
	ENDDO

	kint=0
50	kint=kint+1

c	Sweep in x direction
c	====================
	DO l=2,ny2m1
		l1=l-1
		DO k=2,nx2m1
			k1=k-1
			aj(k1)=aw2(k1,l1)
			bj(k1)=ap2(k1,l1)
			cj(k1)=ae2(k1,l1)
			dj(k1)=-as2(k1,l1)*dp2(k,l1)-an2(k1,l1)*dp2(k,l+1)
     &			 +ar2(k1,l1)
		ENDDO

		aj(1)=0.d0
		cj(nx2m2)=0.d0
		CALL tridag(nx2m2,aj,bj,cj,dj)
		DO k=1,nx2m2
			dp2(k+1,l)=dj(k)
		ENDDO
	ENDDO

c	Block correction along x direction
c	==================================
	DO k=2,nx2m1
		k1=k-1

		aj(k1)=0.d0
		bj(k1)=0.d0
		cj(k1)=0.d0
		dj(k1)=0.d0

		DO l=2,ny2m1
			l1=l-1
			aj(k1)=aj(k1)+aw2(k1,l1)
			bj(k1)=bj(k1)+ap2(k1,l1)
			IF(l.NE.2) bj(k1)=bj(k1)+as2(k1,l1)
			IF(l.NE.ny2m1) bj(k1)=bj(k1)+an2(k1,l1)
			cj(k1)=cj(k1)+ae2(k1,l1)
			dj(k1)=dj(k1)+ar2(k1,l1)-aw2(k1,l1)*dp2(k1,l)
     &			-ae2(k1,l1)*dp2(k+1,l)-as2(k1,l1)*dp2(k,l1)
     &			-an2(k1,l1)*dp2(k,l+1)-ap2(k1,l1)*dp2(k,l)
		ENDDO
	ENDDO
	
	aj(1)=0.d0
	cj(nx2m2)=0.d0
	CALL tridag(nx2m2,aj,bj,cj,dj)
	DO k=1,nx2m2
		DO l=1,ny2m2
			dp2(k+1,l+1)=dp2(k+1,l+1)+dj(k)
		ENDDO
	ENDDO

c	Sweep in y direction
c	====================
	DO k=2,nx2m1
		k1=k-1
		DO l=2,ny2m1
			l1=l-1
			aj(l1)=as2(k1,l1)
			bj(l1)=ap2(k1,l1)
			cj(l1)=an2(k1,l1)
			dj(l1)=-aw2(k1,l1)*dp2(k1,l)-ae2(k1,l1)*dp2(k+1,l)
     &			+ar2(k1,l1)
		ENDDO
		aj(1)=0.d0
		cj(ny2m2)=0.d0
		CALL tridag(ny2m2,aj,bj,cj,dj)
		DO l=1,ny2m2
		  dp2(k,l+1)=dj(l)
		ENDDO
	ENDDO

c	Block correction along y direction
c	==================================
	DO l=2,ny2m1
		l1=l-1

		aj(l1)=0.d0
		bj(l1)=0.d0
		cj(l1)=0.d0
		dj(l1)=0.d0

		DO k=2,nx2m1
			k1=k-1
			aj(l1)=aj(l1)+as2(k1,l1)
			bj(l1)=bj(l1)+ap2(k1,l1)
			IF(k.NE.2) bj(l1)=bj(l1)+aw2(k1,l1)
			IF(k.NE.nx2m1) bj(l1)=bj(l1)+ae2(k1,l1)
			cj(l1)=cj(l1)+an2(k1,l1)
			dj(l1)=dj(l1)+ar2(k1,l1)-aw2(k1,l1)*dp2(k1,l)
     &			-ae2(k1,l1)*dp2(k+1,l)-as2(k1,l1)*dp2(k,l1)
     &			-an2(k1,l1)*dp2(k,l+1)-ap2(k1,l1)*dp2(k,l)
		ENDDO
	ENDDO
	
	aj(1)=0.d0
	cj(ny2m2)=0.d0
	CALL tridag(ny2m2,aj,bj,cj,dj)
	DO l=1,ny2m2
		DO k=1,nx2m2
			dp2(k+1,l+1)=dp2(k+1,l+1)+dj(l)
		ENDDO
	ENDDO

	akd=0.d0
	akn=0.d0
	DO k=2,nx2m1
		k1=k-1
		DO l=2,ny2m1
			l1=l-1
			aj(l1)=as2(k1,l1)
			bj(l1)=ap2(k1,l1)
			cj(l1)=an2(k1,l1)
			dj(l1)=-aw2(k1,l1)*dp2(k1,l)-ae2(k1,l1)*dp2(k+1,l)
     &			 +ar2(k1,l1)
			akd=akd+DABS(aj(l1)*dp2(k,l1)+bj(l1)*dp2(k,l)+
     &			cj(l1)*dp2(k,l+1)-dj(l1))
			akn=akn+DABS(bj(l1)*dp2(k,l))
		ENDDO
	ENDDO
	ak=akd/akn
	
	IF(kint.EQ.1) THEN
		akpre=ak
		GOTO 50
	ENDIF

	ratio=ak/akpre
	akpre=ak
	
	DEALLOCATE(aj)
      DEALLOCATE(bj)
      DEALLOCATE(cj)
      DEALLOCATE(dj)

	RETURN
	END

c======================================================================
	SUBROUTINE acmup(aw2,ae2,as2,an2,ap2,ar2,dp2,dp3,
     &	nx2,ny2,nx3,ny3,nx2max,ny2max,nx3max,ny3max,nsweep)
c======================================================================

	IMPLICIT REAL*8(a-h,o-z)

c     Arguments
c     =========	
	DOUBLE PRECISION  aw2(nx2max,ny2max),ae2(nx2max,ny2max),
     &			        as2(nx2max,ny2max),an2(nx2max,ny2max),
     &			        ap2(nx2max,ny2max),ar2(nx2max,ny2max),
     &			        dp2(nx2max,ny2max),dp3(nx3max,ny3max)
     
c     Local Variables
c     ===============     
      DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE :: aj,bj,cj,dj

      
	nx2m1=nx2-1
	ny2m1=ny2-1
	nx2m2=nx2-2
	ny2m2=ny2-2
	nx3m1=nx3-1
	ny3m1=ny3-1
	
	ALLOCATE(aj(nx2m2))
      ALLOCATE(bj(nx2m2))
      ALLOCATE(cj(nx2m2))
      ALLOCATE(dj(nx2m2))
	
	aj = 0.0d0
	bj = 0.0d0
	cj = 0.0d0
	dj = 0.0d0	

	DO k=2,nx3m1
		DO l=2,ny3m1
			i1=2*(k-1)
			i2=2*(k-1)+1
			j1=2*(l-1)
			j2=2*(l-1)+1
			dp2(i1,j1)=dp2(i1,j1)+dp3(k,l)
			dp2(i2,j1)=dp2(i2,j1)+dp3(k,l)
			dp2(i1,j2)=dp2(i1,j2)+dp3(k,l)
			dp2(i2,j2)=dp2(i2,j2)+dp3(k,l)
		ENDDO
	ENDDO

	DO ksweep=1,nsweep

c		Sweep in x direction
c		====================
		DO l=2,ny2m1
			l1=l-1
			DO k=2,nx2m1
				k1=k-1
				aj(k1)=aw2(k1,l1)
				bj(k1)=ap2(k1,l1)
				cj(k1)=ae2(k1,l1)
				dj(k1)=-as2(k1,l1)*dp2(k,l1)-an2(k1,l1)*dp2(k,l+1)
     &				+ar2(k1,l1)
			ENDDO
			aj(1)=0.d0
			cj(nx2m2)=0.d0
			CALL tridag(nx2m2,aj,bj,cj,dj)
			DO k=1,nx2m2
				dp2(k+1,l)=dj(k)
			ENDDO
		ENDDO

c		Block correction along x direction
c		==================================
		DO k=2,nx2m1
			k1=k-1

			aj(k1)=0.d0
			bj(k1)=0.d0
			cj(k1)=0.d0
			dj(k1)=0.d0

			DO l=2,ny2m1
				l1=l-1
				aj(k1)=aj(k1)+aw2(k1,l1)
				bj(k1)=bj(k1)+ap2(k1,l1)
				IF(l.NE.2) bj(k1)=bj(k1)+as2(k1,l1)
				IF(l.NE.ny2m1) bj(k1)=bj(k1)+an2(k1,l1)
				cj(k1)=cj(k1)+ae2(k1,l1)
				dj(k1)=dj(k1)+ar2(k1,l1)-aw2(k1,l1)*dp2(k1,l)
     &				-ae2(k1,l1)*dp2(k+1,l)-as2(k1,l1)*dp2(k,l1)
     &				-an2(k1,l1)*dp2(k,l+1)-ap2(k1,l1)*dp2(k,l)
			ENDDO
		ENDDO
	
		aj(1)=0.d0
		cj(nx2m2)=0.d0
		CALL tridag(nx2m2,aj,bj,cj,dj)
	
		DO k=1,nx2m2
			DO l=1,ny2m2
				dp2(k+1,l+1)=dp2(k+1,l+1)+dj(k)
			ENDDO
		ENDDO

c		Sweep in y direction
c		====================
		DO k=2,nx2m1
			k1=k-1
			DO l=2,ny2m1
				l1=l-1
				aj(l1)=as2(k1,l1)
				bj(l1)=ap2(k1,l1)
				cj(l1)=an2(k1,l1)
				dj(l1)=-aw2(k1,l1)*dp2(k1,l)-ae2(k1,l1)*dp2(k+1,l)
     &				+ar2(k1,l1)
			ENDDO

			aj(1)=0.d0
			cj(ny2m2)=0.d0
			CALL tridag(ny2m2,aj,bj,cj,dj)
			DO l=1,ny2m2
				dp2(k,l+1)=dj(l)
			ENDDO
		ENDDO

c		Block correction along y direction
c		==================================
		DO l=2,ny2m1
			l1=l-1

			aj(l1)=0.d0
			bj(l1)=0.d0
			cj(l1)=0.d0
			dj(l1)=0.d0

			DO k=2,nx2m1
				k1=k-1
				aj(l1)=aj(l1)+as2(k1,l1)
				bj(l1)=bj(l1)+ap2(k1,l1)
				IF(k.NE.2) bj(l1)=bj(l1)+aw2(k1,l1)
				IF(k.NE.nx2m1) bj(l1)=bj(l1)+ae2(k1,l1)
				cj(l1)=cj(l1)+an2(k1,l1)
				dj(l1)=dj(l1)+ar2(k1,l1)-aw2(k1,l1)*dp2(k1,l)
     &				-ae2(k1,l1)*dp2(k+1,l)-as2(k1,l1)*dp2(k,l1)
     &				-an2(k1,l1)*dp2(k,l+1)-ap2(k1,l1)*dp2(k,l)
			ENDDO
		ENDDO
		aj(1)=0.d0
		cj(ny2m2)=0.d0
		CALL tridag(ny2m2,aj,bj,cj,dj)
		DO l=1,ny2m2
			DO k=1,nx2m2
				dp2(k+1,l+1)=dp2(k+1,l+1)+dj(l)
			ENDDO
		ENDDO

	ENDDO
	
	DEALLOCATE(aj)
      DEALLOCATE(bj)
      DEALLOCATE(cj)
      DEALLOCATE(dj)

	RETURN
	END

c======================================================================
	SUBROUTINE flow(pn,qn)									
c======================================================================

c----------------------------------------------------------------------
c	This subroutine calculates the poiseuille flow						 
c	under the slider based on the boltzman equation						 
c	The first and second derivatives are also 						 
c	computed. an asymptotic form is used for large values.				 
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE siml_pram	
      USE syst_cnfg
      USE rynl_data
      USE sldr_dynm															 

	IMPLICIT REAL*8(a-h,o-z) 

	DOUBLE PRECISION a(20),b(20)

c	Use the continuum model
c	=======================
	IF (iqpo.EQ.0) THEN
		qn	= pn
		RETURN

c	Use the first order slip model
c	==============================
	ELSEIF (iqpo.EQ.1) THEN
		qn	= pn + t1
		RETURN
c	Use the second order slip model
c	===============================
	ELSEIF (iqpo.EQ.2) THEN
		r1	= t2/pn
		r2	= r1/pn
		qn	= pn + t1 + r1
		RETURN

c	Use the asymtotic fk model
c	==========================
	ELSEIF (iqpo.EQ.3) THEN
		r1	= t2/pn
		s1	= t3/(pn*pn)
		u1	= t4/(pn*pn*pn)
		qn	= pn + t1 + r1 - s1 - u1
		RETURN

c	Use the full fk model
c	=====================
	ELSEIF (iqpo.EQ.4) THEN
		z  = d0*pn
		IF (icoe.EQ.1) GOTO 15
		a(1) = 0.d0
		a(2) = -1.d0
		b(1) = -pir
		b(2) = 1.5d0*(1.d0-gama)
		DO i = 3,nter
			xi	= dfloat(i)
			aij = xi*(xi-1.d0)*(xi-2.d0)
			a(i) = (-2.d0*a(i-2))/aij
			b(i) = (-2.d0*b(i-2)-(3.d0*xi*xi-6.d0*xi+2.d0)*a(i))/aij
		ENDDO

		icoe = 1
   15		IF (z.GE.1.1d0)	THEN

			v = 1.889881575d0*z**0.66666666666666666d0
			v1 = 12.d0*v
			v2 = 24.d0*v*v1
			v3 = 36.d0*v*v2
			v4 = 48.d0*v*v3
			v5 = 60.d0*v*v4
			ev = pit*dexp(-v)
			ff0 = ev*(1.d0-1.d0/v1+25.d0/v2-1225.d0/v3
     &			+89425.d0/v4-7263025.d0/v5)
			f1 = ev*dsqrt(v/3.d0)*(1.d0+5.d0/v1-35.d0/v2
     &			+665.d0/v3+9625.d0/v4-9284275.d0/v5)
			f2 = ev*(v/3.d0)*(1.d0+17.d0/v1-35.d0/v2
     &			+1925.d0/v3-175175.d0/v4+22247225.d0/v5)
		ELSE

		ff0	 = 0.d0
		f1	= 0.d0
		f2	= 0.d0
		alz	= dlog(z)
		DO i = nter,1,-1
			xi	= dfloat(i)
			aij = 1.d0/(1.d0+xi)
			IF (i.EQ.1) THEN
				 ff0 = ((alz*xi+1.d0)*a(i) + b(i)*xi + ff0)
				 GOTO 12
			ENDIF

			ff0 = ((alz*xi+1.d0)*a(i) + b(i)*xi + ff0)*z
   12 		f1 = (a(i)*alz + b(i) + f1)*z

			f2 = (aij*((alz-aij)*a(i)+b(i)) + f2)*z

		ENDDO
	   
		ff0 = -0.5d0*ff0
		f1 = 0.5d0*(f1 + 1.d0)
		f2 = pir/4.d0 - 0.5d0*(f2 + 1.d0)*z
		ENDIF
		c11 = 8.d0-z*z*z*(pir/12.d0-z/16.d0)-2.d0*z*(4.d0+z*z)*ff0
     *		-(16.d0+z*z*(8.d0+z*z/8.d0))*f1-z*(16.d0+z*z)*f2
		c22 = 1.d0- 2.d0*f1
		c12 = 2.d0 - (pir/2.d0)*z + z*z/4.d0-2.d0*z*ff0
     *		-(4.d0 +z*z/2.d0)*f1 - 2.d0*z*f2
		del = c11*c22 - c12*c12

		q2 = (pir/del)*(c11 - z*z*(c12/6.d0-z*z*c22/144.d0))
		qn = (-1.d0/z+q2)*(6.d0/z)*pn

	   RETURN


c	Use the data base of precalculated qp based on fk model
c	=======================================================
	ELSEIF (iqpo.EQ.5) THEN
		z0=d0*pn
		loc=dint(200.d0*dlog10(DABS(z0))+401)
		IF (loc.GE.801) GOTO 55
		IF(loc.LT.1)GOTO 56
		z1=zdat(loc)

		qn1=qndat(loc)
		z2=zdat(loc+1)
		qn2=qndat(loc+1)

c		Interpolate
c		===========   
		qn=(qn2-qn1)*(z0-z1)/(z2-z1)+qn1

		RETURN

c		IF inverse knudsen no. is out of range of database
c		use first order slip
c		==================================================
 55		CONTINUE
		r1	= t2/pn
		s1	= t3/(pn*pn)
		u1	= t4/(pn*pn*pn)
		qn	= pn + t1 + r1 - s1 - u1
		RETURN
 56		CONTINUE
		qn=-pn*6.d0/pir*dlog(z0)/z0
		RETURN

	ENDIF

	END
	 

c======================================================================
	SUBROUTINE create_dbase
c======================================================================

c----------------------------------------------------------------------
c	This subroutine creates the data base of precalculated
c	values of poiseuille flow rate as a function of the 
c	inverse knudsen no.  the data base is modelled after 
c	table 2 of fukui and kaneko (j. of tribology vol 112
c	jan 1990, p78-83)
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE syst_cnfg
      USE rynl_data
      USE sldr_dynm

	IMPLICIT REAL*8(a-h,o-z)

	z =0.01d0
	dz=10.d0**(0.005d0)

c	Temporarily change iqpo to 4 to create data base
c	================================================
	iqpo=4

	DO i=1,801
		pn=z/d0
		CALL flow(pn,qn)

		zdat(i)	= z
		qndat(i)= qn
		z=z*dz
	ENDDO

c	Change iqpo back to 5
c	=====================
	iqpo=5

	RETURN
	END

c======================================================================
	SUBROUTINE calc_bearing_number
c======================================================================

c     Shared Data
c     ===========
      USE siml_pram	
      USE syst_cnfg
      USE sldr_grid
      USE rynl_data
      USE sldr_dynm
      USE actr_data

	IMPLICIT REAL*8(a-h,o-z)

c	Actuator parameters
c	===================
	IF (iact.EQ.1) THEN

c	  Inline actuator
c	  ===============
		dex=xact*dsin(dact+ske_init)
		dey=xact*(dcos(ske_init)-dcos(dact+ske_init))
		raref=ra_init+xact*dsin(ske_init)
		alfaa=atan(dey/(raref-dex))
		ske=alfaa+dact+ske_init
		ra=sqrt(dey*dey+(raref-dex)*(raref-dex))
		IF(ra.LT.ra_if) GOTO 310
		IF(ra.GT.ra_of) GOTO 320
		ssk = dsin(ske)
		csk = dcos(ske)
		uax = vact*xact*dsin(ske_init)
		uay = vact*xact*dcos(ske_init)
		DO i=1,nx
			xs0=(xg-xref(i))*xl
			DO j=1,ny
				yc0=(0.5d0*yl-yref(j))*xl
				ux=omega*(ra*csk-yc0)-uax
				uy=omega*(ra*ssk+xs0)+uay
				bearx(i,j)=6.d0*amu*ux*xl/(p0*hm*hm)
				beary(i,j)=6.d0*amu*uy*xl/(p0*hm*hm)
			ENDDO
		ENDDO

		skeeft=datan2(omega*ra*ssk+uay,omega*ra*csk-uax)

	ELSE

c	  No actuator
c	  ===========
		ssk = dsin(ske)
		csk = dcos(ske)
		DO i=1,nx
			xs0=(xg-xref(i))*xl
			DO j=1,ny
				yc0=(0.5d0*yl-yref(j))*xl
				ux=omega*(ra*csk-yc0)
				uy=omega*(ra*ssk+xs0)
				bearx(i,j)=6.d0*amu*ux*xl/(p0*hm*hm)
				beary(i,j)=6.d0*amu*uy*xl/(p0*hm*hm)
			ENDDO
		ENDDO

	ENDIF
	
	RETURN
  310 WRITE(*,311)
  311 FORMAT ('Radial Location < Inner flyable radius')
	STOP
  320 WRITE(*,321)
  321 FORMAT ('Radial Location > Outer flyable radius')
	STOP
	END
c======================================================================