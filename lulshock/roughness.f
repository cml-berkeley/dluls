c**********************************************************************
c     Subroutines in this file :
c     1. get_slider
c     2. fl_ms
c     3. fl_an
c     4. disktop
c     5. dtopog
c     6. wasp
c     7. aasp
c**********************************************************************

c======================================================================
	SUBROUTINE get_slider
c======================================================================

c----------------------------------------------------------------------
c     Computes slider roughness and asperities : hsad(i,j)
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE siml_pram	
      USE syst_cnfg
      USE sldr_grid
      USE sldr_arys
      USE sldr_rghn
      
	IMPLICIT REAL*8(a-h,o-z)

c	Roughness wave
c	==============

	DO k=1,ns_wave
		DO i=1,nx
			DO j=1,ny

c				Type 0.0: Sine Wave 
c				===================
				IF(stype(k).EQ.0.d0) THEN
		 
					IF(sxwlth(k).EQ.0.d0) THEN
						zx=1.d0
					ELSE
						xtemp=xref(i)*dcos(sang(k))+yref(j)*
     &						dsin(sang(k))
						zx=dsin(twopi*xtemp/sxwlth(k))
					ENDIF

					IF(sywlth(k).EQ.0.d0) THEN
						zy=1.d0
					ELSE
						ytemp=yref(j)*dcos(sang(k))-xref(i)*
     &						dsin(sang(k))
						zy=dsin(twopi*ytemp/sywlth(k))
					ENDIF

					IF(sxwlth(k).NE.0.d0.OR.sywlth(k).NE.0.d0) THEN
						hsad(i,j)=hsad(i,j)+smag(k)*zx*zy
					ENDIF

c				Type 1.0: Elliptical Asperity Wave				 
c				==================================

				ELSE IF(stype(k).EQ.1.d0) THEN
					xtemp=xref(i)*dcos(sang(k))+yref(j)*dsin(sang(k))
					ytemp=yref(j)*dcos(sang(k))-xref(i)*dsin(sang(k))
					xtemp=dmod(DABS(xtemp),sxwlth(k))
					ytemp=dmod(DABS(ytemp),sywlth(k))
					xtemp=DMIN1(DABS(xtemp),DABS(sxwlth(k)-xtemp))
					ytemp=DMIN1(DABS(ytemp),DABS(sywlth(k)-ytemp))
					IF(sxdm(k).NE.0.d0.AND.sydm(k).NE.0.d0) THEN
						temp=xtemp*xtemp/sxdm(k)/sxdm(k)+			  
     &						ytemp*ytemp/sydm(k)/sydm(k) - 1.d0
					ELSEIF(sxdm(k).EQ.0.d0.AND.sydm(k).NE.0.d0) THEN
						temp=ytemp*ytemp/sydm(k)/sydm(k) - 1.d0
					ELSEIF(sxdm(k).NE.0.d0.AND.sydm(k).EQ.0.d0) THEN
						temp=xtemp*xtemp/sxdm(k)/sxdm(k) - 1.d0
					ELSE
						temp=0.d0
					ENDIF

					IF(temp.LE.0.d0) THEN
						hsad(i,j)=hsad(i,j)-smag(k)*temp
					ENDIF


c				Type 3.0: Cylindrical Asperity Wave
c				===================================
				ELSE IF(stype(k).EQ.2.d0) THEN

					xtemp=xref(i)*dcos(sang(k))+yref(j)*dsin(sang(k))
					ytemp=yref(j)*dcos(sang(k))-xref(i)*dsin(sang(k))
					xtemp=dmod(DABS(xtemp),sxwlth(k))
					ytemp=dmod(DABS(ytemp),sywlth(k))
					xtemp=DMIN1(DABS(xtemp),DABS(sxwlth(k)-xtemp))
					ytemp=DMIN1(DABS(ytemp),DABS(sywlth(k)-ytemp))
					
					IF(sxdm(k).NE.0.d0.AND.sydm(k).NE.0.d0) THEN
						temp=xtemp*xtemp/sxdm(k)/sxdm(k)+			  
     &						ytemp*ytemp/sydm(k)/sydm(k) - 1.d0
					ELSEIF(sxdm(k).EQ.0.d0.AND.sydm(k).NE.0.d0) THEN
						temp=ytemp*ytemp/sydm(k)/sydm(k) - 1.d0
					ELSEIF(sxdm(k).NE.0.d0.AND.sydm(k).EQ.0.d0) THEN
						temp=xtemp*xtemp/sxdm(k)/sxdm(k) - 1.d0
					ELSE
						temp=0.d0
					ENDIF
			  
					IF(temp.LT.0.d0) THEN
						hsad(i,j)=hsad(i,j)+smag(k)
					ENDIF


c				Type 3.0: Rectangular Asperity Wave
c				===================================
				ELSE IF(stype(k).EQ.3.d0) THEN

					xtemp=xref(i)*dcos(sang(k))+yref(j)*dsin(sang(k))
					ytemp=yref(j)*dcos(sang(k))-xref(i)*dsin(sang(k))
					xtemp=dmod(xtemp,sxwlth(k))
					ytemp=dmod(ytemp,sywlth(k))
			
					IF(xtemp.LT.0.d0) xtemp=xtemp+sxwlth(k)
						
					IF(ytemp.LT.0.d0) ytemp=ytemp+sywlth(k)
							
					IF(sxdm(k).NE.0.d0.AND.sydm(k).NE.0.d0) THEN
						IF(xtemp.LE.sxdm(k).AND.ytemp.LE.sydm(k))
     &						hsad(i,j)=hsad(i,j)+smag(k)
					ELSEIF(sxdm(k).NE.0.d0.AND.sydm(k).EQ.0.d0) THEN
						IF(xtemp.LE.sxdm(k)) 
     &						hsad(i,j)=	hsad(i,j)+smag(k)
					ELSEIF(sxdm(k).EQ.0.d0.AND.sydm(k).NE.0.d0) THEN
						IF(ytemp.LE.sydm(k)) 
     &						hsad(i,j)=hsad(i,j)+smag(k)
					ENDIF

				ENDIF
			ENDDO
		ENDDO
	ENDDO

c	Individual asperities
c	=====================
	DO k=1,ns_asper
		DO i=1,nx
			DO j=1,ny

c				Type 0.0: Sinusoidal Bump
c				=========================
				IF(satype(k).EQ.0.d0) THEN

					temp=0.d0
					xtemp=xref(i)*dcos(saang(k))+yref(j)*
     &					dsin(saang(k))
					ytemp=yref(j)*dcos(saang(k))-xref(i)*
     &					dsin(saang(k))
					xloc=sxloc(k)*dcos(saang(k))+syloc(k)*
     &					dsin(saang(k))
					yloc=syloc(k)*dcos(saang(k))-sxloc(k)*
     &					dsin(saang(k))
					xtemp=DABS(xloc-xtemp)
					ytemp=DABS(yloc-ytemp)
			  
					IF(saxdm(k).NE.0.d0.AND.saydm(k).NE.0.d0) THEN
						IF(xtemp.LE.saxdm(k).AND.ytemp.LE.saydm(k))
     &					   temp=dsin((1.d0+xtemp/saxdm(k))*twopi/4.d0)
					ELSEIF(saxdm(k).NE.0.d0.AND.saydm(k).EQ.0.d0) THEN
						IF(xtemp.LE.saxdm(k))
     &						temp=dsin((1.d0+xtemp/saxdm(k))*twopi/4.d0)
					ELSE
						temp=0.d0
					ENDIF
					hsad(i,j)=hsad(i,j)+temp*samag(k)

c				Type 1.0: Elliptical Asperity
c				=============================
				ELSE IF(satype(k).EQ.1.d0) THEN

					xtemp=xref(i)*dcos(saang(k))+yref(j)*
     &					dsin(saang(k))
					ytemp=yref(j)*dcos(saang(k))-xref(i)*
     &					dsin(saang(k))
					xloc=sxloc(k)*dcos(saang(k))+syloc(k)*
     &					dsin(saang(k))
					yloc=syloc(k)*dcos(saang(k))-sxloc(k)*
     &					dsin(saang(k))
					xtemp=DABS(xloc-xtemp)
					ytemp=DABS(yloc-ytemp)
				
					IF(saxdm(k).NE.0.d0.AND.saydm(k).NE.0.d0) THEN
						temp=xtemp*xtemp/saxdm(k)/saxdm(k)+ 			
     &						ytemp*ytemp/saydm(k)/saydm(k) - 1.d0
					ELSEIF(saxdm(k).EQ.0.d0.AND.saydm(k).NE.0.d0) THEN
						temp=ytemp*ytemp/saydm(k)/saydm(k) - 1.d0
					ELSEIF(saxdm(k).NE.0.d0.AND.saydm(k).EQ.0.d0) THEN
						temp=xtemp*xtemp/saxdm(k)/saxdm(k) - 1.d0
					ELSE
						temp=0.d0
					ENDIF
				
					IF(temp.LE.0.d0) THEN
						hsad(i,j)=hsad(i,j)-samag(k)*temp
					ENDIF

c				Type 2.0: Cylindrical Asperity
c				==============================
				ELSE IF(satype(k).EQ.2.d0) THEN

					xtemp=xref(i)*dcos(saang(k))+yref(j)*
     &					dsin(saang(k))
					ytemp=yref(j)*dcos(saang(k))-xref(i)*
     &					dsin(saang(k))
					xloc=sxloc(k)*dcos(saang(k))+syloc(k)*
     &					dsin(saang(k))
					yloc=syloc(k)*dcos(saang(k))-sxloc(k)*
     &					dsin(saang(k))
					xtemp=DABS(xloc-xtemp)
					ytemp=DABS(yloc-ytemp)
				
					IF(saxdm(k).NE.0.d0.AND.saydm(k).NE.0.d0) THEN
						temp=xtemp*xtemp/saxdm(k)/saxdm(k)+ 			
     &					ytemp*ytemp/saydm(k)/saydm(k) - 1.d0
					ELSEIF(saxdm(k).EQ.0.d0.AND.saydm(k).NE.0.d0) THEN
						temp=ytemp*ytemp/saydm(k)/saydm(k) - 1.d0
					ELSEIF(saxdm(k).NE.0.d0.AND.saydm(k).EQ.0.d0) THEN
						temp=xtemp*xtemp/saxdm(k)/saxdm(k) - 1.d0
					ELSE
						temp=0.d0
					ENDIF
				
					IF(temp.LT.0.d0) THEN
						hsad(i,j)=hsad(i,j)+samag(k)
					ENDIF

c				Type 3.0: Rectangular Asperity
c				==============================
				ELSE IF(satype(k).EQ.3.d0) THEN

					temp=0.d0
					xtemp=xref(i)*dcos(saang(k))+yref(j)*
     &					dsin(saang(k))
					ytemp=yref(j)*dcos(saang(k))-xref(i)*
     &					dsin(saang(k))
					xloc=sxloc(k)*dcos(saang(k))+syloc(k)*
     &					dsin(saang(k))
					yloc=syloc(k)*dcos(saang(k))-sxloc(k)*
     &					dsin(saang(k))
					xtemp=DABS(xloc-xtemp)
					ytemp=DABS(yloc-ytemp)
			  
					IF(saxdm(k).NE.0.d0.AND.saydm(k).NE.0.d0) THEN
						IF(xtemp.LE.saxdm(k).AND.ytemp.LE.saydm(k))
     &					temp=1.d0
					ELSEIF(saxdm(k).NE.0.d0.AND.saydm(k).EQ.0.d0) THEN
						IF(xtemp.LE.saxdm(k)) temp=1.d0
					ELSEIF(saxdm(k).EQ.0.d0.AND.saydm(k).NE.0.d0) THEN
						IF(ytemp.LE.saydm(k)) temp=1.d0
					ELSE
						temp=0.d0
					ENDIF
					hsad(i,j)=hsad(i,j)+temp*samag(k)

				ENDIF
			ENDDO
		ENDDO
	ENDDO

	RETURN
	END

c===========================================================================
	SUBROUTINE fl_ms
c===========================================================================

c---------------------------------------------------------------------------
c	 Read measured track profile
c---------------------------------------------------------------------------

c     Shared Data
c     ===========      
      USE siml_pram
      USE syst_cnfg
      USE disk_rghn
      USE trck_prfl
      USE sldr_dynm
      
	IMPLICIT REAL*8(a-h,o-z)
		
	IF(inwave.EQ.0) RETURN
		
	OPEN(5,ERR=999,FILE='wave.def',STATUS='OLD')
	DO i=1,nfx
		READ(5,*)xfref(i),hfnew(i)
		xfref(i)=xfref(i)/xl
		hfnew(i)=hfnew(i)/hm
	ENDDO

	xfl=xfref(nfx)

	RETURN
999	WRITE(*,995)
995	FORMAT ('Trouble in opening files')
	STOP
	END
	  
	  
c===========================================================================
	SUBROUTINE fl_an(xv,yvi,thh)
c===========================================================================

c---------------------------------------------------------------------------
c     Compute the thh at specified coordinates (xv,yvi) on slider
c---------------------------------------------------------------------------

c     Shared Data
c     ===========      
      USE siml_pram	
      USE syst_cnfg
      USE disk_rghn
      USE trck_prfl
      USE sldr_dynm
      
	IMPLICIT REAL*8(a-h,o-z)

      DOUBLE PRECISION flhgtn(nwtp4),falhgtn(natp4)
      DOUBLE PRECISION fllocn(nwtp4),fallocn(natp4)

	thh=0.d0

c	Roughness wave
c	==============

	DO k=1,nf_wave
		
		IF((yvi.GT.fras(k).AND.yvi.LT.frae(k)).OR.
     &		(yvi.LT.fras(k).AND.yvi.GT.frae(k))) THEN

		yv=yvi-fras(k)

c		Type 0.0: Sine Wave
c		===================			
		IF(ftype(k).EQ.0.d0) THEN

			
			IF(fxwlth(k).EQ.0.d0) THEN
				zx=1.d0
			ELSE
				xtemp=xv*dcos(fang(k))+yv*dsin(fang(k))
				zx=dsin(twopi*xtemp/fxwlth(k))
			ENDIF

			IF(fywlth(k).EQ.0.d0) THEN
				zy=1.d0
			ELSE
				ytemp=yv*dcos(fang(k))-xv*dsin(fang(k))
				zy=dsin(twopi*ytemp/fywlth(k))
			ENDIF
	
			IF(fxwlth(k).NE.0.d0.OR.fywlth(k).NE.0.d0) THEN
				thh=thh+fmag(k)*zx*zy
			ENDIF

c		Type 1.0: Elliptical Asperity Wave
c		==================================
		ELSE IF(ftype(k).EQ.1.d0) THEN

			xtemp=xv*dcos(fang(k))+yv*dsin(fang(k))
			ytemp=yv*dcos(fang(k))-xv*dsin(fang(k))
			xtemp=dmod(DABS(xtemp),fxwlth(k))
			ytemp=dmod(DABS(ytemp),fywlth(k))
			xtemp=DMIN1(DABS(xtemp),DABS(fxwlth(k)-xtemp))
			ytemp=DMIN1(DABS(ytemp),DABS(fywlth(k)-ytemp))
			
			IF(fxdm(k).NE.0.d0.AND.fydm(k).NE.0.d0) THEN
				temp=xtemp*xtemp/fxdm(k)/fxdm(k)+
     &				ytemp*ytemp/fydm(k)/fydm(k)-1.d0
			ELSEIF(fxdm(k).NE.0.d0.AND.fydm(k).EQ.0.d0) THEN
				temp=xtemp*xtemp/fxdm(k)/fxdm(k)-1.d0
			ELSEIF(fxdm(k).EQ.0.d0.AND.fydm(k).NE.0.d0) THEN
				temp=ytemp*ytemp/fydm(k)/fydm(k)-1.d0
			ELSE
				temp=1.d0
			ENDIF
			
			IF(temp.LE.0.d0) THEN
				thh=thh-fmag(k)*temp
			ENDIF


c		Type 2.0: Cylindrical Aperity Wave
c		==================================
		ELSE IF(ftype(k).EQ.2.d0) THEN

			xtemp=xv*dcos(fang(k))+yv*dsin(fang(k))
			ytemp=yv*dcos(fang(k))-xv*dsin(fang(k))
			xtemp=dmod(DABS(xtemp),fxwlth(k))
			ytemp=dmod(DABS(ytemp),fywlth(k))
			xtemp=DMIN1(DABS(xtemp),DABS(fxwlth(k)-xtemp))
			ytemp=DMIN1(DABS(ytemp),DABS(fywlth(k)-ytemp))
			
			IF(fxdm(k).NE.0.d0.AND.fydm(k).NE.0.d0) THEN
				temp=xtemp*xtemp/fxdm(k)/fxdm(k)+
     &				ytemp*ytemp/fydm(k)/fydm(k)-1.d0
			ELSEIF(fxdm(k).NE.0.d0.AND.fydm(k).EQ.0.d0) THEN
				temp=xtemp*xtemp/fxdm(k)/fxdm(k)-1.d0
			ELSEIF(fxdm(k).EQ.0.d0.AND.fydm(k).NE.0.d0) THEN
				temp=ytemp*ytemp/fydm(k)/fydm(k)-1.d0
			ELSE
				temp=1.d0
			ENDIF
			
			IF(temp.LE.0.d0) THEN
				thh=thh+fmag(k)
			ENDIF


c		Type 3.0: Rectangular Asperity Wave
c		===================================
		ELSE IF(ftype(k).EQ.3.d0) THEN
			
			xtemp=xv*dcos(fang(k))+yv*dsin(fang(k))
			ytemp=yv*dcos(fang(k))-xv*dsin(fang(k))
			xtemp=dmod(xtemp,fxwlth(k))
			ytemp=dmod(ytemp,fywlth(k))
		
			IF(xtemp.LT.0.d0) xtemp=xtemp+fxwlth(k)
			IF(ytemp.LT.0.d0) ytemp=ytemp+fywlth(k)
			IF(fxdm(k).NE.0.d0.AND.fydm(k).NE.0.d0) THEN
				IF(xtemp.LE.fxdm(k).AND.ytemp.LE.fydm(k)) 
     &				thh=thh+fmag(k)
			ELSEIF(fxdm(k).NE.0.d0.AND.fydm(k).EQ.0.d0) THEN
				IF(xtemp.LE.fxdm(k)) thh=thh+fmag(k)
			ELSEIF(fxdm(k).EQ.0.d0.AND.fydm(k).NE.0.d0) THEN
				IF(ytemp.LE.fydm(k)) thh=thh+fmag(k)
			ENDIF


c		Type 4.0:  Asperity wave of arbitrary profile (read additional files)
c		=====================================================================
		ELSE IF(ftype(k).EQ.4.d0) THEN

			xtemp=xv*dcos(fang(k))+yv*dsin(fang(k))
			ytemp=yv*dcos(fang(k))-xv*dsin(fang(k))
			xtemp=dmod(DABS(xtemp),fxwlth(k))
			ytemp=dmod(DABS(ytemp),fywlth(k))
			xtemp=DMIN1(DABS(xtemp),DABS(fxwlth(k)-xtemp))
			ytemp=DMIN1(DABS(ytemp),DABS(fywlth(k)-ytemp))
			raxy=dsqrt(xtemp*xtemp+ytemp*ytemp)
		  
			IF(raxy.GE.0.0d0.AND.raxy.LE.flloc(K,nwtp4n)) THEN
				DO NN=1,nwtp4n
					flhgtn(NN)=FLHGT(k,NN)
					fllocn(NN)=FLLOC(K,NN)
				ENDDO
					CALL finteg1(fllocN,flhgtn,nwtp4n,raxy,htemp)
					thh=thh+htemp
				ENDIF

			ENDIF
		ENDIF
	ENDDO

c	Disk zone profile
c	=================
	DO k=1,nf_zone
		IF(yvi.GT.frp1(k).AND.yvi.LT.frp3(k)) THEN
			thh=thh+fcoe2(k)*yvi*yvi+fcoe1(k)*yvi+fcoe0(k)
		ENDIF
	ENDDO


c	Individal Asperities
c	====================
	DO k=1,nf_asper
		yv=yvi-ra/xl

c		Type 0.0: Sinusoidal Bump
c		=========================
		IF(fatype(k).EQ.0.d0) THEN

			temp=0.d0
			xtemp=xv*dcos(faang(k))+yv*dsin(faang(k))
			ytemp=yv*dcos(faang(k))-xv*dsin(faang(k))
			xloc=fxloc(k)*dcos(faang(k))+fyloc(k)*dsin(faang(k))
			yloc=fyloc(k)*dcos(faang(k))-fxloc(k)*dsin(faang(k))
			xtemp=DABS(xloc-xtemp)
			ytemp=DABS(yloc-ytemp)
			
			IF(faxdm(k).NE.0.d0.AND.faydm(k).NE.0.d0) THEN
				IF(xtemp.LE.faxdm(k).AND.ytemp.LE.faydm(k))
     &				temp=dsin((1.d0+xtemp/faxdm(k))*twopi/4.d0)
			ELSEIF(faxdm(k).NE.0.d0.AND.faydm(k).EQ.0.d0) THEN
				IF(xtemp.LE.faxdm(k))
     &				temp=dsin((1.d0+xtemp/faxdm(k))*twopi/4.d0)
			ELSE
				temp=0.d0
			ENDIF
			  thh=thh+temp*famag(k)

c		Type 1.0: Elliptical Asperity
c		=============================
		ELSE IF(fatype(k).EQ.1.d0) THEN

			xtemp=xv*dcos(faang(k))+yv*dsin(faang(k))
			ytemp=yv*dcos(faang(k))-xv*dsin(faang(k))
			xloc=fxloc(k)*dcos(faang(k))+fyloc(k)*dsin(faang(k))
			yloc=fyloc(k)*dcos(faang(k))-fxloc(k)*dsin(faang(k))
			xtemp=DABS(xloc-xtemp)
			ytemp=DABS(yloc-ytemp)
			
			IF(faxdm(k).NE.0.d0.AND.faydm(k).NE.0.d0) THEN
				temp=xtemp*xtemp/faxdm(k)/faxdm(k)+
     &			 ytemp*ytemp/faydm(k)/faydm(k)-1.d0
		  	ELSEIF(faxdm(k).NE.0.d0.AND.faydm(k).EQ.0.d0) THEN
				temp=xtemp*xtemp/faxdm(k)/faxdm(k)-1.d0
			ELSEIF(faxdm(k).EQ.0.d0.AND.faydm(k).NE.0.d0) THEN
				temp=ytemp*ytemp/faydm(k)/faydm(k)-1.d0
			ELSE
				temp=1.d0
			ENDIF
			
			IF(temp.LE.0.d0) THEN
				thh=thh-famag(k)*temp
			ENDIF

c		Type 2.0: Cylindrical Aperity
c		=============================
		ELSE IF(fatype(k).EQ.2.d0) THEN

			xtemp=xv*dcos(faang(k))+yv*dsin(faang(k))
			ytemp=yv*dcos(faang(k))-xv*dsin(faang(k))
			xloc=fxloc(k)*dcos(faang(k))+fyloc(k)*dsin(faang(k))
			yloc=fyloc(k)*dcos(faang(k))-fxloc(k)*dsin(faang(k))
			xtemp=DABS(xloc-xtemp)
			ytemp=DABS(yloc-ytemp)

			IF(faxdm(k).NE.0.d0.AND.faydm(k).NE.0.d0) THEN
				temp=xtemp*xtemp/faxdm(k)/faxdm(k)+
     &				ytemp*ytemp/faydm(k)/faydm(k)-1.d0
			ELSEIF(faxdm(k).NE.0.d0.AND.faydm(k).EQ.0.d0) THEN
				temp=xtemp*xtemp/faxdm(k)/faxdm(k)-1.d0
			ELSEIF(faxdm(k).EQ.0.d0.AND.faydm(k).NE.0.d0) THEN
				temp=ytemp*ytemp/faydm(k)/faydm(k)-1.d0
			ELSE
				temp=1.d0
			ENDIF
			
			IF(temp.LE.0.d0) THEN
				thh=thh+famag(k)
			ENDIF

c		Type 3.0: Rectangular Asperity
c		==============================
		ELSE IF(fatype(k).EQ.3.d0) THEN

			temp=0.d0
			xtemp=xv*dcos(faang(k))+yv*dsin(faang(k))
			ytemp=yv*dcos(faang(k))-xv*dsin(faang(k))
			xloc=fxloc(k)*dcos(faang(k))+fyloc(k)*dsin(faang(k))
			yloc=fyloc(k)*dcos(faang(k))-fxloc(k)*dsin(faang(k))
			xtemp=DABS(xloc-xtemp)
			ytemp=DABS(yloc-ytemp)
			IF(faxdm(k).NE.0.d0.AND.faydm(k).NE.0.d0) THEN
				IF(xtemp.LE.faxdm(k).AND.ytemp.LE.faydm(k))
     &				temp=1.d0
			ELSEIF(faxdm(k).NE.0.d0.AND.faydm(k).EQ.0.d0) THEN
				IF(xtemp.LE.faxdm(k)) temp=1.d0
			ELSEIF(faxdm(k).EQ.0.d0.AND.faydm(k).NE.0.d0) THEN
				IF(ytemp.LE.faydm(k)) temp=1.d0
			ELSE
				temp=0.d0
			ENDIF

			thh=thh+temp*famag(k)

c		Type 4.0:  Asperity of arbitrary cross-section profile
c		======================================================
		ELSE IF(fatype(k).EQ.4.d0) THEN

			xtemp=xv*dcos(faang(k))+yv*dsin(faang(k))
			ytemp=yv*dcos(faang(k))-xv*dsin(faang(k))
			xloc=fxloc(k)*dcos(faang(k))+fyloc(k)*dsin(faang(k))
			yloc=fyloc(k)*dcos(faang(k))-fxloc(k)*dsin(faang(k))
			xtemp=DABS(xloc-xtemp)
			ytemp=DABS(yloc-ytemp)
			raxy=dsqrt(xtemp*xtemp+ytemp*ytemp)
		  
			IF(raxy.GE.0.0d0.AND.raxy.LE.falloc(k,natp4n)) THEN
			
				DO NN=1,NATP4n
					falhgtn(NN)=FALHGT(K,NN)
					fallocn(NN)=FALLOC(K,NN)
				ENDDO
				CALL finteg1(fallocN,falhgtN,natp4n,raxy,htemp)
					thh=thh+htemp
			ENDIF

		ENDIF
	ENDDO

	RETURN
	END
	  
c===========================================================================
	SUBROUTINE disktop
c===========================================================================

c---------------------------------------------------------------------------
c     Output disk topology at the specific time using hfloor
c---------------------------------------------------------------------------

c     Shared Data
c     ===========      
      USE siml_pram
      USE sldr_grid
      USE sldr_arys
      USE disk_rghn
      USE trck_prfl
      USE rynl_data
      USE sldr_dynm

	IMPLICIT REAL*8(a-h,o-z)
	
	IF(nf_wave.NE.0.OR.nf_asper.NE.0.OR.nf_zone.NE.0) THEN
		
		IF((t.GE.tout1).AND.((t-tout1).LT.dt)) THEN
			OPEN(22,file='disktop1.dat',status='unknown',
     &								  recl=5000)
			IF(igp1.EQ.1) THEN
				DO i=1,nx-1
					ii=ixtr(i)
					WRITE(22,700)(-hm*hfloor(ii,jytr(j)),j=1,ny-1)
				END DO
			ELSE
				DO i=1,nx
					WRITE(22,700)(-hm*hfloor(i,j),j=1,ny)
				ENDDO
			ENDIF
			CLOSE(22)

		ELSE IF((t.GE.tout2).AND.((t-tout2).LT.dt)) THEN
			OPEN(23,file='disktop2.dat',status='unknown',
     &								  recl=5000)
			IF(igp1.EQ.1) THEN
				DO i=1,nx-1
					ii=ixtr(i)
					WRITE(23,700)(-hm*hfloor(ii,jytr(j)),j=1,ny-1)
				END DO
			ELSE
				DO i=1,nx
					WRITE(23,700)(-hm*hfloor(i,j),j=1,ny)
				ENDDO
			ENDIF
			CLOSE(23)

		ELSE IF((t.GE.tout3).AND.((t-tout3).LT.dt)) THEN
			OPEN(24,file='disktop3.dat',status='unknown',
     &								  recl=5000)
			IF(igp1.EQ.1) THEN
				DO i=1,nx-1
					ii=ixtr(i)
					WRITE(24,700)(-hm*hfloor(ii,jytr(j)),j=1,ny-1)
				ENDDO
			ELSE
				DO i=1,nx
					WRITE(24,700)(-hm*hfloor(i,j),j=1,ny)
				ENDDO
			ENDIF
		    CLOSE(24)
		ENDIF
	ENDIF

 700	FORMAT(300(E16.9,1X))
	RETURN
	END

c======================================================================
	SUBROUTINE dtopog
c======================================================================

c     Shared Data
c     ===========
      USE siml_pram	
      USE syst_cnfg
      USE sldr_grid
      USE disk_rghn
      USE sldr_dynm

	IMPLICIT REAL*8(a-h,o-z)

	DOUBLE PRECISION  xfl(300),wave3d(300)

	IF(nf_wave.NE.0) THEN
		
		xfloor=DMAX1(fxwlth(1),fywlth(1))

		DO ii=1,nf_wave
			xfloor=DMAX1(xfloor,fxwlth(ii),fywlth(ii))
		ENDDO
		
		DO i=1,300
			xfl(i)=2.d0*(i-1)*xfloor/300.d0
		ENDDO

c       2D waviness
c       ===========		
        OPEN(17,err=999,file='xfw.dat',status='unknown',
     &		 recl=4800)
        WRITE(17,700) (xl*xfl(i),i=1,300)
		CLOSE(17)

c       3D waviness
c       ===========
		OPEN(18,err=999,file='wave3d.dat',status='unknown',
     &		 recl=4800)		
		DO j=1,300
			DO i=1,300
				xv=0.5d0*xfloor+xfl(i)
				yv=fras(1)+0.5d0*xfloor+xfl(j)
				CALL get_floor(xv,yv,htemp)
				wave3d(i)=htemp
			ENDDO
			WRITE(18,700)(-hm*wave3d(i),i=1,300)
		ENDDO
		CLOSE(18)

	ENDIF

	IF(nf_asper.NE.0) THEN


		IF(fatype(1).NE.4.d0) THEN
			xfloor=DMAX1(faxdm(1),faydm(1))
		ELSE
			xfloor=2.d0*falloc(1,natp4n)
		ENDIF
   
		DO ii=1,nf_asper
			IF(fatype(ii).NE.4.d0) THEN
				xfloor=DMAX1(xfloor,faxdm(ii),faydm(ii))
			ELSE
				xfloor=DMAX1(xfloor,2.d0*falloc(ii,natp4n))
			ENDIF
		ENDDO
		
		DO i=1,300
			xfl(i)=1.5d0*(i-1)*xfloor/300.d0
		ENDDO

c       2D asperity
c       ===========		
	  OPEN(19,err=999,file='xfa.dat',status='unknown',
     &		 recl=4800)		
		WRITE(19,700) (xl*xfl(i),i=1,300)
		CLOSE(19)

c       3D asperity
c       ===========
		OPEN(20,err=999,file='asper3d.dat',status='unknown',
     &		 recl=4800)		
		DO j=1,300
			DO i=1,300
				xv=fxloc(1)-0.75d0*xfloor+xfl(i)
				yv=ra/xl+fyloc(1)-0.75d0*xfloor+xfl(j)
				CALL get_floor(xv,yv,htemp)
				wave3d(i)=htemp
			ENDDO
			WRITE(20,700)(-hm*wave3d(i),i=1,300)
		ENDDO
		CLOSE(20)
		
	ENDIF

	RETURN
700	FORMAT(300(E12.3,1X))
999	WRITE(*,*)'Trouble in opening files'
	END
	
c======================================================================
	SUBROUTINE wasp(ko,fname,k)
c======================================================================

c----------------------------------------------------------------------
c     Read cross-sectional profile for fwtype=4
c----------------------------------------------------------------------

      USE syst_cnfg
      USE disk_rghn
      USE sldr_dynm

	IMPLICIT NONE
	INTEGER i,k,ko
	CHARACTER*50 fname

	OPEN(ko,file=fname,status='OLD',form='FORMATTED')
	REWIND(ko)
	READ(ko,*) nwtp4n
	  
	DO i=1,nwtp4n
		READ(ko,*) flloc(k,i),flhgt(k,i)
		flloc(k,i) = flloc(k,i)/xl  !Normalize
		flhgt(k,i) = flhgt(k,i)/hm  !Normalize
	ENDDO
	
	CLOSE(ko)
	RETURN
	
	END
	
c======================================================================
      SUBROUTINE aasp(ko,fname,k)
c======================================================================

c----------------------------------------------------------------------
c	Read cross-sectional profile for fatype=4
c----------------------------------------------------------------------

      USE syst_cnfg
      USE disk_rghn
      USE sldr_dynm
      
	IMPLICIT NONE
	INTEGER i,k,ko
	CHARACTER*50 fname

	OPEN(ko,file=fname,status='OLD',form='FORMATTED')
	REWIND(ko)
	READ(ko,*) natp4n
	
	DO i=1,natp4n
		READ(ko,*) falloc(k,i),falhgt(k,i)
		falloc(k,i) = falloc(k,i)/xl    !Normalize
		falhgt(k,i) = falhgt(k,i)/hm    !Normalize
	ENDDO
	
	CLOSE(ko)
	RETURN
	
	END
c======================================================================