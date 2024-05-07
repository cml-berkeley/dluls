c**********************************************************************
c     Subroutines in this file :
c     1. calc_ABS_force
c     2. calc_contact_force
c     3. calc_impact_force
c     4. calc_contcoeff
c     5. calc_electrostatic
c**********************************************************************

c======================================================================
      SUBROUTINE calc_ABS_force(iwrite)
c======================================================================

c----------------------------------------------------------------------
c	Net force on slider due to ABS Pressure
c----------------------------------------------------------------------
c----------------------------------------------------------------------
c     This subroutine integrates the pressure profile p(x,y,t)
c     to obtain the resultant force, f, and moment, t                     
c----------------------------------------------------------------------

c     Shared Data
c     ===========      
      USE siml_pram
      USE syst_cnfg
      USE rail_data
      USE sldr_grid
      USE sldr_arys
      USE trck_prfl
      USE rynl_data
      USE sldr_dynm
      USE sldr_stat
      USE aspr_data
      USE vdwf_data
      USE luls_data

      IMPLICIT REAL*8(a-h,o-z)

c     Calculate the air bearing force
c	===============================

	IF (imf_opt .EQ. 0) THEN

		fl=0.d0
		tl=0.d0
		rl=0.d0
		fp =0.d0
		fn =0.d0

		fshear = 0.d0
		tshear = 0.d0
		rshear = 0.d0
		zshear = 0.d0

		ssk=dsin(ske)
		csk=dcos(ske)

		tln=0.0d0
		rln=0.0d0
		tlp=0.0d0
		rlp=0.0d0
	  
		f=0.d0

		IF(omega.EQ.0.d0) RETURN

		IF(igp1.EQ.1) THEN
			nxxm1=nxm1-1
			nyym1=nym1-1
		ELSE
			nxxm1=nxm1
			nyym1=nym1
		ENDIF

		DO jj=2,nyym1
			j=jytr(jj)
			jm1=jytr(jj-1)
			jp1=jytr(jj+1)
			dely=(yref(jp1)-yref(jm1))/2.d0
			DO ii=2,nxxm1
				i=ixtr(ii)
				im1=ixtr(ii-1)
				ip1=ixtr(ii+1)
				delx=(xref(ip1)-xref(im1))/2.d0
				dxdy=delx*dely
				a=(p(i,j)-1.d0)*dxdy
				fl=fl+a
				tl=tl+a*xref(i)
				rl=rl+a*(yref(j)-0.5d0*yl)

				fshear = fshear + 
     &				1.85e-5*omega*ra/hnew(i,j)/hm*dxdy*xl*xl

				zshear = zshear + (1.85e-5*omega*ra/hnew(i,j)/
     &            		 hm*dxdy*xl*xl)*(csk*(yref(j)-0.5)*
     &            		 xl+ssk*(xref(i)-0.5)*xl)

				IF (a.LE.0.d0) THEN
					fn=fn+a
              
					tln=tln+a*xref(i)
					rln=rln+a*(yref(j)-0.5d0*yl)
				ELSE
					fp=fp+a
         
					tlp=tlp+a*xref(i)
					rlp=rlp+a*(yref(j)-0.5d0*yl)
				END IF
			END DO
		END DO

		IF (fp.LE.1.0d-14) THEN
			fp=1d-14
			tlp=.5d-14
		ENDIF
		
		IF (fn.GE.-1.0d-14) THEN
			fn=-1d-14
			tln=-.5d-14
		ENDIF
		
		ABS_data(1) = tlp/fp
		ABS_data(2) = rlp/fp
		ABS_data(3) = tln/fn
		ABS_data(4) = rln/fn
        
		f=fl*p0xl*9.81d0

		IF (DABS(f*1000d0/9.81d0).LE.1.0d-7) THEN
			f=sign(9.81d0/1000d0*1.0d-7,f)
			fl=f/(p0xl*9.81d0)
		ENDIF
      	
      	xf=tl/fl
		yf=rl/fl
		fn=fn*p0xl
		fp=fp*p0xl

c       Output the results
c	  ==================
		IF (iwrite.EQ.1) THEN
			WRITE(*,*)
			WRITE(*,649)
			WRITE(*,655) fp,fn
			WRITE(*,675) f/9.81d0
			WRITE(*,676) xf
			WRITE(*,677) yf
		END IF

		tshear = fshear*csk*zg*xl
		rshear = fshear*ssk*zg*xl

		fvdw = 0.d0
		xfvdw = 0.d0
		yfvdw = 0.d0

	ELSE IF (imf_opt .EQ. 1) THEN

		fl=0.d0
		tl=0.d0
		rl=0.d0
		fp =0.d0
		fn =0.d0


		fshear = 0.d0
		tshear = 0.d0
		rshear = 0.d0
        zshear = 0.d0

        ssk=dsin(ske)
        csk=dcos(ske)

		tln=0.0d0
	  rln=0.0d0
	  tlp=0.0d0
	  rlp=0.0d0

	  f=0.d0

	  IF(omega.EQ.0.d0) RETURN

        IF(igp1.EQ.1) THEN
			nxxm1=nxm1-1
			nyym1=nym1-1
		ELSE
			nxxm1=nxm1
			nyym1=nym1
		ENDIF

		DO jj=2,nyym1
			j=jytr(jj)
			jm1=jytr(jj-1)
			jp1=jytr(jj+1)
			dely=(yref(jp1)-yref(jm1))/2.d0
			DO ii=2,nxxm1
				i=ixtr(ii)
				im1=ixtr(ii-1)
				ip1=ixtr(ii+1)
				delx=(xref(ip1)-xref(im1))/2.d0
				dxdy=delx*dely
				a=(p(i,j)-1.d0)*dxdy
				fl=fl+a
				tl=tl+a*xref(i)
				rl=rl+a*(yref(j)-0.5d0*yl)

				fshear = fshear + 1.85e-5*omega*ra/hnew(i,j)
     &					/hm*dxdy*xl*xl

				zshear = zshear + (1.85e-5*omega*ra/hnew(i,j)
     &					/hm*dxdy*xl*xl)
     &					*(csk*(yref(j)-0.5)*xl+ssk*(xref(i)-0.5)*xl)

				IF (a.LE.0.d0) THEN
					fn=fn+a
            
					tln=tln+a*xref(i)
					rln=rln+a*(yref(j)-0.5d0*yl)
				ELSE
					fp=fp+a
             
					tlp=tlp+a*xref(i)
					rlp=rlp+a*(yref(j)-0.5d0*yl)
				END IF
			END DO
		END DO

		IF (fp.LE.1.0d-14) THEN
			fp=1d-14
			tlp=.5d-14
		ENDIF

		IF (fn.GE.-1.0d-14) THEN
			fn=-1d-14
			tln=-.5d-14
		ENDIF
	   
		ABS_data(1) = tlp/fp
		ABS_data(2) = rlp/fp
		ABS_data(3) = tln/fn
		ABS_data(4) = rln/fn
       
        f=fl*p0xl*9.81d0

	  IF (DABS(f*1000d0/9.81d0).LE.1.0d-7) THEN
			f=sign(9.81d0/1000d0*1.0d-7,f)
			fl=f/(p0xl*9.81d0)
		ENDIF
      	
      	xf=tl/fl
		yf=rl/fl
		fn=fn*p0xl
		fp=fp*p0xl
 
c	  New variables for VWF
c	  =====================
	  fvdw	= 0.d0
	  tvdw1	= 0.d0
	  tvdw2	= 0.d0
	  rvdw1	= 0.d0
	  rvdw2	= 0.d0
	  fposvdw = 0.d0
	  fnegvdw	= 0.d0
	  tlnew	= 0.d0
	  rlnew	= 0.d0
	  flnew	= 0.d0       

		hamcon1 = ahc/6/3.1416d0
	  hamcon2 = bhc/45/3.1416d0

c	  vdw normalization function
c	  ==========================
	  vdwnormf = xl*xl/9.81d0

		IF(omega.EQ.0.d0) RETURN

		IF(igp1.EQ.1) THEN
      		nxx=nx-1
			nyy=ny-1
		ELSE
      		nxx=nx
		    nyy=ny
		ENDIF

		arbit=32.0d-11
		DO i=1,nx
			DO j=1,ny
				IF (spa(i,j).LE.arbit) THEN
					spa(i,j)=100
				END IF
			END DO
		END DO

		DO jj=1,nyy
			j=jj
			jm1=jj-1
			jp1=jj+1
			IF(jj.EQ.1) THEN
      			dely=(yref(jp1)-yref(j))/2.d0
			ELSE IF(jj.EQ.nyy) THEN
				dely=(yref(j)-yref(jm1))/2.d0
			ELSE
				dely=(yref(jp1)-yref(jm1))/2.d0
			ENDIF

			DO ii=1,nxx
	          i=ii
	          im1=ii-1
	          ip1=ii+1
				IF(ii.EQ.1) THEN
					delx=(xref(ip1)-xref(i))/2.d0
				ELSE IF(ii.EQ.nxx) THEN
					delx=(xref(i)-xref(im1))/2.d0
				ELSE
					delx=(xref(ip1)-xref(im1))/2.d0
				ENDIF
				dxdy=delx*dely

c	Calculate depending on where on the slider we are
c	boundary is different
c	=================================================
				IF(i.EQ.1.AND.j.EQ.1) THEN
					tempgap1=(spa(i,j)**3+spa(i,jp1)**3
     &      				+spa(ip1,j)**3+spa(ip1,jp1)**3)/4.d0
					tempgap2=(spa(i,j)**9+spa(i,jp1)**9
     &         			+spa(ip1,j)**9+spa(ip1,jp1)**9)/4.d0

				ELSE IF(i.EQ.1.AND.j.EQ.ny) THEN
					tempgap1=(spa(i,j)**3+spa(i,jm1)**3
     &         			+spa(ip1,j)**3+spa(ip1,jm1)**3)/4.d0
					tempgap2=(spa(i,j)**9+spa(i,jm1)**9
     &         			+spa(ip1,j)**9+spa(ip1,jm1)**9)/4.d0

				ELSE IF(i.EQ.nx.AND.j.EQ.1) THEN
					tempgap1=(spa(i,j)**3+spa(i,jp1)**3
     &         			+spa(im1,j)**3+spa(im1,jp1)**3)/4.d0
					tempgap2=(spa(i,j)**9+spa(i,jp1)**9
     &         			+spa(im1,j)**9+spa(im1,jp1)**9)/4.d0

				ELSE IF(i.EQ.1.AND.j.NE.1.AND.j.NE.ny) THEN
					tempgap1=(spa(i,jm1)**3+spa(ip1,jm1)**3
     &         			+spa(i,j)**3+spa(ip1,j)**3
     &         			+spa(i,jp1)**3+spa(ip1,jp1)**3)/6.d0
					tempgap2=(spa(i,jm1)**9+spa(ip1,jm1)**9
     &         			+spa(i,j)**9+spa(ip1,j)**9
     &         			+spa(i,jp1)**9+spa(ip1,jp1)**9)/6.d0

				ELSE IF(j.EQ.1.AND.i.NE.1.AND.i.NE.nx) THEN
					tempgap1=(spa(im1,j)**3+spa(im1,jp1)**3
     &         			+spa(i,j)**3+spa(i,jp1)**3
     &         			+spa(ip1,j)**3+spa(ip1,jp1)**3)/6.d0
					tempgap2=(spa(im1,j)**9+spa(im1,jp1)**9
     &         			+spa(i,j)**9+spa(i,jp1)**9
     &         			+spa(ip1,j)**9+spa(ip1,jp1)**9)/6.d0

				ELSE IF(i.EQ.nx.AND.j.NE.1.AND.j.NE.ny) THEN
					tempgap1=(spa(i,jm1)**3+spa(im1,jm1)**3
     &         			+spa(i,j)**3+spa(im1,j)**3
     &         			+spa(i,jp1)**3+spa(im1,jp1)**3)/6.d0
					tempgap2=(spa(i,jm1)**9+spa(im1,jm1)**9
     &         			+spa(i,j)**9+spa(im1,j)**9
     &         			+spa(i,jp1)**9+spa(im1,jp1)**9)/6.d0

				ELSE IF(j.EQ.ny.AND.i.NE.1.AND.i.NE.nx) THEN
					tempgap1=(spa(im1,j)**3+spa(im1,jm1)**3
     &         			+spa(i,j)**3+spa(i,jm1)**3
     &         			+spa(ip1,j)**3+spa(ip1,jm1)**3)/6.d0
					tempgap2=(spa(im1,j)**9+spa(im1,jm1)**9
     &         			+spa(i,j)**9+spa(i,jm1)**9
     &         			+spa(ip1,j)**9+spa(ip1,jm1)**9)/6.d0
			
				ELSE
					tempgap1=(spa(im1,jm1)**3+spa(im1,j)**3
     &					+spa(im1,jp1)**3
     &        			+spa(i,jm1)**3+spa(i,j)**3+spa(i,jp1)**3
     &        			+spa(ip1,jm1)**3+spa(ip1,j)**3
     &					+spa(ip1,jp1)**3)/9.d0
					tempgap2=(spa(im1,jm1)**9+spa(im1,j)**9
     &					+spa(im1,jp1)**9
     &	         		+spa(i,jm1)**3+spa(i,j)**9+spa(i,jp1)**9
     &      				+spa(ip1,jm1)**9+spa(ip1,j)**9
     &					+spa(ip1,jp1)**9)/9.d0
				END IF
				
				vdw1=hamcon1*dxdy/(tempgap1)
			    vdw2=hamcon2*dxdy/(tempgap2)

		        fvdw=fvdw-vdw1+vdw2

				tvdw1=tvdw1-vdw1*(xref(i))
				tvdw2=tvdw2+vdw2*(xref(i))

				rvdw1=rvdw1-vdw1*(yref(j)-yl/2.d0)
				rvdw2=rvdw2+vdw2*(yref(j)-yl/2.d0)

         			fnegvdw=fnegvdw-vdw1
				fposvdw=fposvdw+vdw2

				IF(DABS(fvdw).GT.1d10) THEN
					f_yoyo = 0.0
				ENDIF

         		ENDDO
		ENDDO

		fvdw=fposvdw+fnegvdw

		yfvdw=(tvdw1+tvdw2)*vdwnormf
		xfvdw=(rvdw1+rvdw2)*vdwnormf
		fvdw=fvdw*vdwnormf
		fposvdw=fposvdw*vdwnormf
		fnegvdw=fnegvdw*vdwnormf

c		tlnew=tl*p0xl+xfvdw
c		rlnew=rl*p0xl+yfvdw
c		flnew=fl*p0xl+fvdw

		fvdw = fvdw*9.81d0
		xfvdw = xfvdw*9.81d0*xl
		yfvdw = yfvdw*9.81d0*xl

c		f=fl*9.81d0

c		xf=tlnew/flnew
c		yf=rlnew/flnew

c		fp=fp+fposvdw
c		fn=fn+fnegvdw

c		Output the results
c		==================
		IF (iwrite.EQ.1) THEN
			WRITE(*,*)
			WRITE(*,649)
			WRITE(*,655) fp,fn
			WRITE(*,675) f/9.81d0
			WRITE(*,676) xf
			WRITE(*,677) yf
		END IF

		tshear = fshear*csk*zg*xl
		rshear = fshear*ssk*zg*xl

	ENDIF

649   FORMAT('Air-bearing results:')
655   FORMAT(/' positive force = ',G16.9,'; negative force = ',
     &         G16.9/)
675   FORMAT(' resultant force in kgf f  = ',E16.9)
676   FORMAT(' location (norm. by xl) xf = ',E16.9)
677   FORMAT(' location (norm. by xl) yf = ',E16.9)
678   FORMAT(3(E12.5,1X))
        
	RETURN
      END
 
c======================================================================
	SUBROUTINE calc_contact_force(iwrite)
c======================================================================

c----------------------------------------------------------------------
c	Net force on slider due to asperity contact 
c     i.e. Min Clearance < Glide Height
c----------------------------------------------------------------------

c     Shared Data
c     ===========      
      USE siml_pram
      USE syst_cnfg
      USE rail_data
      USE sldr_grid
      USE sldr_arys
      USE trck_prfl
      USE rynl_data
      USE sldr_dynm
      USE aspr_data
      USE impc_data
      
	IMPLICIT REAL*8(a-h,o-z)
	
	EXTERNAL fit1,fitb1,fit32

	COMMON/var/hint,dfctn
	DOUBLE PRECISION tgs(110),fgs1(110),fgs2(110)
c
c     Create the functions for the greenwood-williamson model
c	=======================================================
	  data tgs/0.00d0,0.05d0,0.10d0,0.15d0,0.20d0,0.25d0,0.30d0,
     &		   0.35d0,0.40d0,0.45d0,0.50d0,0.55d0,0.60d0,0.65d0,
     &		   0.70d0,0.75d0,0.80d0,0.85d0,0.90d0,0.95d0,1.00d0,
     &		   1.05d0,1.10d0,1.15d0,1.20d0,1.25d0,1.30d0,1.35d0,
     &		   1.40d0,1.45d0,1.50d0,1.55d0,1.60d0,1.65d0,1.70d0,
     &		   1.75d0,1.80d0,1.85d0,1.90d0,1.95d0,2.00d0,2.05d0,
     &		   2.10d0,2.15d0,2.20d0,2.25d0,2.30d0,2.35d0,2.40d0,
     &		   2.45d0,2.50d0,2.55d0,2.60d0,2.65d0,2.70d0,2.75d0,
     &		   2.80d0,2.85d0,2.90d0,2.95d0,3.00d0,3.05d0,3.10d0,
     &		   3.15d0,3.20d0,3.25d0,3.30d0,3.35d0,3.40d0,3.45d0,
     &		   3.50d0,3.55d0,3.60d0,3.65d0,3.70d0,3.75d0,3.80d0,
     &		   3.85d0,3.90d0,3.95d0,4.00d0,4.05d0,4.10d0,4.15d0,
     &		   4.20d0,4.25d0,4.30d0,4.35d0,4.40d0,4.45d0,4.50d0,
     &		   4.55d0,4.60d0,4.65d0,4.70d0,4.75d0,4.80d0,4.85d0,
     &		   4.90d0,4.95d0,5.00d0,5.50d0,6.00d0,7.00d0,8.00d0,
     &		   9.00d0,10.0d0,20.0d0,500.d0,1000000000.d0/
	data  fgs1/0.3989d0,0.3744d0,0.3509d0,0.3284d0,0.3069d0,0.2863d0,
     &		   0.2668d0,0.2481d0,0.2304d0,0.2137d0,0.1978d0,0.1828d0,
     &		   0.1687d0,0.1554d0,0.1429d0,0.1312d0,0.1202d0,0.1100d0,
     &		   0.1004d0,0.9156d-1,0.8332d-1,0.7568d-1,0.6862d-1,
     &		   0.6210d-1,0.5610d-1,0.5059d-1,0.4553d-1,0.4090d-1,
     &		   0.3667d-1,0.3281d-1,0.2930d-1,0.2612d-1,0.2324d-1,
     &		   0.2064d-1,0.1829d-1,0.1617d-1,0.1428d-1,0.1257d-1,
     &		   0.1105d-1,0.9698d-2,0.8490d-2,0.7418d-2,0.6468d-2,
     &		   0.5628d-2,0.4887d-2,0.4235d-2,0.3662d-2,0.3159d-2,
     &		   0.2720d-2,0.2337d-2,0.2004d-2,0.1715d-2,0.1464d-2,
     &		   0.1247d-2,0.1060d-2,0.8992d-3,0.7611d-3,0.6428d-3,
     &		   0.5417d-3,0.4555d-3,0.3822d-3,0.3199d-3,0.2673d-3,
     &		   0.2228d-3,0.1852d-3,0.1537d-3,0.1273d-3,0.1051d-3,
     &		   0.8666d-4,0.7127d-4,0.5848d-4,0.4788d-4,0.3911d-4,
     &		   0.3188d-4,0.2592d-4,0.2103d-4,0.1702d-4,0.1375d-4,
     &		   0.1108d-4,0.8908d-5,0.7145d-5,0.5718d-5,0.4566d-5,
     &		   0.3637d-5,0.2891d-5,0.2292d-5,0.1814d-5,0.1432d-5,
     &		   0.1127d-5,0.8857d-6,0.6942d-6,0.5429d-6,0.4236d-6,
     &		   0.3297d-6,0.2560d-6,0.1984d-6,0.1533d-6,0.1182d-6,
     &		   0.9096d-7,0.6982d-7,0.5346d-7,0.3255d-8,0.1564d-9,
     &		   0.1760d-12,0.7550d-16,0.1225d-19,0.7475d-24,
     &		   0.1370d-89,0.0d0,0.0d0/
	data  fgs2/0.4299d0,0.4000d0,0.3715d0,0.3446d0,0.3191d0,0.2952d0,
     &		   0.2725d0,0.2513d0,0.2313d0,0.2127d0,0.1951d0,0.1789d0,
     &		   0.1636d0,0.1495d0,0.1363d0,0.1241d0,0.1127d0,0.1023d0,
     &		   0.9267d-1,0.8382d-1,0.7567d-1,0.6819d-1,0.6132d-1,
     &		   0.5508d-1,0.4935d-1,0.4417d-1,0.3944d-1,0.3517d-1,
     &		   0.3129d-1,0.2779d-1,0.2463d-1,0.2180d-1,0.1925d-1,
     &		   0.1697d-1,0.1493d-1,0.1311d-1,0.1149d-1,0.1005d-1,
     &		   0.8773d-2,0.7646d-2,0.6646d-2,0.5769d-2,0.4995d-2,
     &		   0.4319d-2,0.3724d-2,0.3207d-2,0.2754d-2,0.2362d-2,
     &		   0.2020d-2,0.1725d-2,0.1469d-2,0.1250d-2,0.1060d-2,
     &		   0.8979d-3,0.7587d-3,0.6396d-3,0.5380d-3,0.4518d-3,
     &		   0.3784d-3,0.3164d-3,0.2639d-3,0.2197d-3,0.1825d-3,
     &		   0.1513d-3,0.1251d-3,0.1032d-3,0.8500d-4,0.6984d-4,
     &		   0.5725d-4,0.4684d-4,0.3823d-4,0.3113d-4,0.2529d-4,
     &		   0.2051d-4,0.1660d-4,0.1340d-4,0.1079d-4,0.8670d-5,
     &		   0.6952d-5,0.5561d-5,0.4438d-5,0.3535d-5,0.2809d-5,
     &		   0.2227d-5,0.1762d-5,0.1391d-5,0.1095d-5,0.8603d-6,
     &		   0.6743d-6,0.5274d-6,0.4115d-6,0.3204d-6,0.2488d-6,
     &		   0.1928d-6,0.1491d-6,0.1150d-6,0.8851d-7,0.6796d-7,
     &		   0.5206d-7,0.3979d-7,0.3034d-7,0.1774d-8,0.8204d-10,
     &		   0.8621d-13,0.3478d-16,0.5341d-20,0.3101d-24,
     &		   0.4059d-90,0.0d0,0.0d0/
	
c	Calculate the contact force
c	===========================

	atot=0.d0
	ac  =0.d0
	fcr =0.d0
	f_maxcrp = 0.d0
	txr =0.d0
	tyr =0.d0
      tzr =0.d0

	IF(ncz.EQ.0.OR.eyoung.EQ.0.d0) RETURN

	DO k=1,ncz
		IF((ra.GE.rcts(k).AND.ra.LE.rcte(k)).OR.
     &		(ra.LE.rcts(k).AND.ra.GE.rcte(k))) THEN
			rsik=rsikm(k)
			cta=ceta(k)
			rasp=rasper(k)
			ght=gldht(k)
			GOTO 2
		ENDIF
	ENDDO
	
	RETURN

 2	IF(rsik.EQ.0.d0.OR.cta.EQ.0.d0.OR.rasp.EQ.0.d0) RETURN
	
	IF(omega.EQ.0.d0) THEN
		fmotion=0.d0
	ELSE
		fmotion=1.d0
	ENDIF

	hcg=2.d0*zg*xl
	ght=ght/rsik

	IF(igp1.EQ.1) THEN
		nxx=nx-1
		nyy=ny-1
	ELSE
		nxx=nx
		nyy=ny
	ENDIF

c	Greewood-Williamson contact model
c	=================================
	IF (icmod.EQ.1) THEN

		hcm=hm/rsik
		ssk=dsin(ske)
		csk=dcos(ske)
		hp=-hx0*hm/xl
		hyy=hy*hm/xl
		coe1=4.0d0*dsqrt(rasp*rsik)*cta*eyoung*rsik*xl*xl/3.d0

		DO jj=1,nyy
			j=jj
			IF(jj.EQ.1) THEN
				dely=(yref(j+1)-yref(j))/2.d0
			ELSE IF(jj.EQ.nyy) THEN
				dely=(yref(j)-yref(j-1))/2.d0
			ELSE
				dely=(yref(j+1)-yref(j-1))/2.d0
			ENDIF
		
			DO ii=1,nxx
				i=ii
				
				IF(ii.EQ.1) THEN
					delx=(xref(i+1)-xref(i))/2.d0
				ELSE IF(ii.EQ.nxx) THEN
					delx=(xref(i)-xref(i-1))/2.d0
				ELSE
					delx=(xref(i+1)-xref(i-1))/2.d0
				ENDIF
				
				dxdy=delx*dely
				atot=atot+dxdy
  
				hint=hnew(i,j)
				IF ((hxlnew(i,j).LT.hint).AND.
     &				(hxlnew(i,j).NE.0d0)) THEN
					hint=hxlnew(i,j)
				ENDIF

				IF ((hxrnew(i,j).LT.hint).AND.
     &				(hxrnew(i,j).NE.0d0)) THEN
					hint=hxrnew(i,j)
				ENDIF

				IF ((hyunew(i,j).LT.hint).AND.
     &				(hyunew(i,j).NE.0d0)) THEN
					hint=hyunew(i,j)
				ENDIF

				IF ((hydnew(i,j).LT.hint).AND.
     &				(hydnew(i,j).NE.0d0)) THEN
					hint=hcm*hydnew(i,j)
				ENDIF

				hint=hcm*hint
	
c	          Glide Height Criteria
c	          =====================
				IF(hint.LE.ght) THEN
                    CALL finteg(tgs,fgs1,fgs2,110,hint,aa,af)
				ELSE
					aa=0.d0
					af=0.d0
				ENDIF
				aa=aa*dxdy
				pcontact(i,j)=coe1*af/xl/xl/p0
				af=coe1*af*dxdy
				ac=ac+aa
				IF (f_maxcrp.LT.pcontact(i,j)) THEN
					f_maxcrp = pcontact(i,j)
				ENDIF
				fcr=fcr+af
				txr=txr+af*(xref(i)-xg)*xl+
     &				fmotion*frcoe*af*csk*hcg
				tyr=tyr+af*(yref(j)-(0.5d0*yl+yg))*xl+
     &				fmotion*frcoe*af*ssk*hcg
                tzr=tzr+af*frcoe*(csk*(xref(i)-0.5d0)*
     &                  xl+ssk*(yref(j)-0.5d0)*xl)
			END DO
		END DO

		ac=twopi*cta*rasp*rsik*xl*xl*ac/2.d0
		
c	Elastic-plastic contact model
c	=============================		
	ELSEIF (icmod.EQ.2) THEN

		dfctn=(rasp*(twopi*ydcoe*ydst/eyoung/4.d0)**2)/rsik
		hcm=hm/rsik
		ssk=dsin(ske)
		csk=dcos(ske)
		hp=-hx0*hm/xl
		hyy=hy*hm/xl

		DO jj=1,nyy
			j=jj
			IF(jj.EQ.1) THEN
				dely=(yref(j+1)-yref(j))/2.d0
			ELSE IF(jj.EQ.nyy) THEN
				dely=(yref(j)-yref(j-1))/2.d0
			ELSE
				dely=(yref(j+1)-yref(j-1))/2.d0
			ENDIF
		
			DO ii=1,nxx
				i=ii
			
				IF(ii.EQ.1) THEN
					delx=(xref(i+1)-xref(i))/2.d0
				ELSE IF(ii.EQ.nxx) THEN
					delx=(xref(i)-xref(i-1))/2.d0
				ELSE
					delx=(xref(i+1)-xref(i-1))/2.d0
				ENDIF
		  
				dxdy=delx*dely
				atot=atot+dxdy
				hint=hcm*hnew(i,j)
				hinte=hint+dfctn
		 
				IF(hint.LE.ght.AND.hint.LT.4.d0) THEN
					CALL qromb(fit1,hint,hinte,fc1)
					fc1=fc1/dsqrt(twopi)
					CALL qromb(fitb1,hinte,6.d0,fcb1)
					fcb1=fcb1/dsqrt(twopi)
					CALL qromb(fit32,hint,hinte,fc32)
					fc32=fc32/dsqrt(twopi)
				ELSE
					fc1=0.d0
					fcb1=0.d0
					fc32=0.d0
				ENDIF
		  
				aa1=fc1+fcb1
				aa2=4.d0*dsqrt(rasp*rsik)*fc32/3.d0+
     &				twopi*rasp*ydcoe*ydst*fcb1/eyoung/2.d0
				pcontact(i,j)=cta*eyoung*rsik*aa2/p0
				ac=ac+twopi*cta*rasp*rsik*xl*xl*aa1*dxdy/2.d0
	      
				IF (f_maxcrp.LT.pcontact(i,j)) THEN
					f_maxcrp = pcontact(i,j)
				ENDIF
			
				fcr=fcr+cta*eyoung*rsik*xl*xl*aa2*dxdy
				txr=txr+cta*rsik*xl*xl*eyoung*aa2*((xref(i)-xg)*xl+
     &				fmotion*hcg*frcoe*csk)*dxdy
				tyr=tyr+cta*rsik*xl*xl*eyoung*aa2*((yref(j)-
     &				(0.5d0*yl+yg))*xl+fmotion*hcg*frcoe*ssk)*dxdy
			END DO
		END DO

	ENDIF

	atot=xl*xl*atot
	aratio=ac/atot

c	Output the results
c	==================
	IF(iwrite.EQ.1) THEN
		WRITE(*,*)
		WRITE(*,750)
		WRITE(*,755) aratio
		WRITE(*,756) fcr
		WRITE(*,759) txr
		WRITE(*,760) tyr
	END IF

750	FORMAT(' Contact results:')
755	FORMAT(' Contact area (norm. by apparent area)  = ',E16.9)
756	FORMAT(' Contact force in kgf fcr= ',E16.9)
759	FORMAT(' Contact moment in kgf*m txr= ',E16.9)
760	FORMAT(' Contact moment in kgf*m tyr= ',E16.9)
	

	RETURN
	END

c======================================================================
      SUBROUTINE calc_impact_force
c======================================================================

c----------------------------------------------------------------------
c	Net impact force: Puneet BEM model
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE siml_pram
      USE syst_cnfg
      USE actr_data
      USE rail_data
      USE sldr_grid
      USE sldr_arys
      USE trck_prfl
      USE rynl_data
      USE sldr_dynm
      USE sldr_stat
      USE aspr_data
      USE impc_data
      USE luls_data
      
      IMPLICIT REAL*8(a-h,o-z)
      
	LOGICAL crash
	
	PARAMETER(nxr=10,nyr=10,maxcp=500,nrad=10)
	
	INTEGER i_xmap(nxr+1),i_ymap(nyr+1)
	INTEGER i_map(nxr,nyr),i_bmap(maxcp,2)
	INTEGER i_sten(4,2,nrad*nrad)
	DOUBLE PRECISION f_hred(nxr,nyr),f_xred(nxr),f_yred(nyr)
	DOUBLE PRECISION f_Ccompliance(maxcp,maxcp),f_Cdisp(maxcp)
	DOUBLE PRECISION f_Cstiffness(maxcp,maxcp),f_pres(maxcp)
	DOUBLE PRECISION f_Vcoords(4,2),f_Vcoeffs(4)
	DOUBLE PRECISION f_xfac(nxr),f_yfac(nyr)

c	Elastic Constant
c	================
      f_const = -(1-2*f_nu)*(1+f_nu)/(2*eycrash*1d-3*3.142)

      f_delx = (xref(nx)-xref(1))/(nxr-1)
      f_dely = (yref(ny)-yref(1))/(nyr-1)
      f_xred(1) = xref(1) + f_delx*1d-6
      
	DO i=2,nxr
		f_xred(i) = (f_xred(i-1)+f_delx)
	ENDDO
      
	f_xred(nxr) = f_xred(nxr) - 2*f_delx*1d-6
      f_yred(1) = yref(1) + f_dely*1d-6
         
	DO j=2,nyr
		f_yred(j) = (f_yred(j-1)+f_dely)
	ENDDO
      
	f_yred(nyr) = f_yred(nyr) - 2*f_dely*1d-6

	iold = 1
	
	DO i=1,nxr
		DO WHILE(f_xred(i)>xref(iold))
			iold = iold+1
	    ENDDO
		i_xmap(i) = iold
	ENDDO
	
	i_xmap(nxr+1) = i_xmap(nxr)
	jold = 1
	
	DO j=1,nyr
		DO WHILE(f_yred(j)>yref(jold))
			jold = jold+1
	    ENDDO
			i_ymap(j) = jold
	ENDDO
	
	i_ymap(nyr+1) = i_ymap(nyr)

	DO i=1,nxr
		DO j=1,nyr
			f_temp = 10.d9
			DO iold=i_xmap(i),i_xmap(i+1)
				DO jold=i_ymap(j),i_ymap(j+1)
					IF (hnew(iold,jold).LE.f_temp) THEN
						f_temp = hnew(iold,jold)
					ENDIF
				ENDDO
			ENDDO
			f_hred(i,j) = f_temp*hm*1d3
		ENDDO
	ENDDO

      DO i=1,nxr
		f_xred(i) = f_xred(i)*xl*1d3
	ENDDO
       
	DO j=1,nyr
		f_yred(j) = f_yred(j)*xl*1d3
	ENDDO
      
	f_delx = f_delx*xl*1d3
      f_dely = f_dely*xl*1d3

      nc = 0
	
	DO i=1,nxr
		DO j=1,nyr
			IF (f_hred(i,j).LT.0.d0) THEN
				nc = nc+1
	            i_map(i,j) = nc
	            i_bmap(nc,1) = i
	            i_bmap(nc,2) = j
	         ELSE
	            i_map(i,j) = 0
	         ENDIF
		ENDDO
	ENDDO

	ncold = nc

      IF (nc.GT.maxcp) THEN
		WRITE(*,*) "Maximum contact points exceeded ",nc
	    STOP
	ENDIF

      DO i_r=1,nrad
		DO j_r=1,nrad
      		i_sten(1,1,(i_r-1)*nrad+j_r) = i_r - nrad
	        i_sten(2,1,(i_r-1)*nrad+j_r) = i_r - nrad + 1
	        i_sten(3,1,(i_r-1)*nrad+j_r) = i_r - nrad + 1
	        i_sten(4,1,(i_r-1)*nrad+j_r) = i_r - nrad
      	    i_sten(1,2,(i_r-1)*nrad+j_r) = j_r - nrad
	        i_sten(2,2,(i_r-1)*nrad+j_r) = j_r - nrad
	        i_sten(3,2,(i_r-1)*nrad+j_r) = j_r - nrad + 1
	        i_sten(4,2,(i_r-1)*nrad+j_r) = j_r - nrad + 1
	    ENDDO
	ENDDO

      DO i=1,nxr
		DO j=1,nyr
			IF (i_map(i,j).GT.0) THEN
                  f_Cdisp(i_map(i,j)) = f_hred(i,j)
				DO k=1,nrad*nrad
					DO m=1,4
						ic = i+i_sten(m,1,k)
						f_Vcoords(m,1) = f_xred(1) + (ic-1)*f_delx
						jc = j+i_sten(m,2,k)
						f_Vcoords(m,2) = f_yred(1) + (jc-1)*f_dely
					ENDDO
					CALL calc_contcoeff(f_Vcoords,f_xred(i),
     &                        f_yred(j),f_Vcoeffs)
					DO m=1,4
						IF ((ic.GT.0).AND.(ic.LE.nxr).AND.
     &						(jc.GT.0).AND.(jc.LE.nyr)) THEN
							f_Ccompliance(i_map(i,j),i_map(ic,jc)) = 
     &                        f_Ccompliance(i_map(i,j),i_map(ic,jc)) +
     &                        f_Vcoeffs(m)*f_const
						ENDIF
					ENDDO
				ENDDO
			ENDIF
		ENDDO
	ENDDO

	i_flag = 1
	
	DO WHILE(i_flag.EQ.1)
	  CALL inverse(f_Ccompliance,nc,maxcp,f_Cstiffness)
        CALL mult(f_Cstiffness,f_Cdisp,nc,maxcp,1,1,f_pres)
	  k = 0
	  i_flag = 0
	  DO WHILE (k<nc)
	      
	      k = k+1
	      IF (f_pres(k).LT.0.d0) THEN
			    i_flag = 1
                i_map(i_bmap(k,1),i_bmap(k,2)) = 0
	      
	          DO i=1,nc
			        DO j=k+1,nc
				        f_Ccompliance(i,j-1) = f_Ccompliance(i,j)
				    ENDDO
			    ENDDO
            
                DO i=k+1,nc
			        f_Cdisp(i-1) = f_Cdisp(i)
				    i_bmap(i-1,1) = i_bmap(i,1)
				    i_bmap(i-1,2) = i_bmap(i,2)
				    DO j=1,nc
					    f_Ccompliance(i-1,j) = f_Ccompliance(i,j)
				    ENDDO
	          ENDDO
            
                nc = nc-1
		    ENDIF
		ENDDO
	ENDDO

	f_minfh = 1.d6
	DO i=1,nxr
		DO j=1,nyr
			IF (f_hred(i,j).LE.f_minfh) THEN
				f_minfh = f_hred(i,j)
			ENDIF
		ENDDO
	ENDDO

	fcont=0.d0
	f_maximp = 0.d0
	fsctx=0.d0
      fscty=0.d0
         
      DO k=1,nc
		IF (f_pres(k).GT.f_maximp) THEN
			f_maximp = f_pres(k)
		ENDIF
	    fct = f_delx*f_dely*f_pres(k)*1d-3
	    fcont = fcont+fct
          fsctx = fsctx-fct*(f_xred(k)-xg*xl*1d3)*1d-3
          fscty = fscty-fct*(f_yred(k)-(0.5d0*yl+yg)*xl)*1d-3
	ENDDO

	RETURN
	END

c======================================================================
      SUBROUTINE calc_contcoeff(f_Vcoords,f_x,f_y,f_Vcoeffs)
c======================================================================
      IMPLICIT REAL*8(a-h,o-z) 
      DOUBLE PRECISION f_Vcoords(4,2), f_Vcoeffs(4)

      f_x1 = f_Vcoords(1,1) - f_x
      f_y1 = f_Vcoords(1,2) - f_y
      f_x2 = f_Vcoords(3,1) - f_x
      f_y2 = f_Vcoords(3,2) - f_y
	   
      f_t111 = f_x1*flog(f_y1+sqrt(f_y1*f_y1+f_x1*f_x1)) +
     &         f_y1*flog(f_x1+sqrt(f_y1*f_y1+f_x1*f_x1))
      f_t112 = f_x1*flog(f_y2+sqrt(f_y2*f_y2+f_x1*f_x1)) +
     &         f_y2*flog(f_x1+sqrt(f_y2*f_y2+f_x1*f_x1))
      f_t122 = f_x2*flog(f_y2+sqrt(f_y2*f_y2+f_x2*f_x2)) +
     &         f_y2*flog(f_x2+sqrt(f_y2*f_y2+f_x2*f_x2))
      f_t121 = f_x2*flog(f_y1+sqrt(f_y1*f_y1+f_x2*f_x2)) +
     &         f_y1*flog(f_x2+sqrt(f_y1*f_y1+f_x2*f_x2))
      f_t211 =(2.d0*f_y1*sqrt(f_y1*f_y1+f_x1*f_x1) + 
     &         2.d0*f_x1*f_x1*flog(f_y1+sqrt(f_y1*f_y1+f_x1*f_x1)))/4.d0
      f_t212 =(2.d0*f_y2*sqrt(f_y2*f_y2+f_x1*f_x1) + 
     &         2.d0*f_x1*f_x1*flog(f_y2+sqrt(f_y2*f_y2+f_x1*f_x1)))/4.d0
      f_t222 =(2.d0*f_y2*sqrt(f_y2*f_y2+f_x2*f_x2) + 
     &         2.d0*f_x2*f_x2*flog(f_y2+sqrt(f_y2*f_y2+f_x2*f_x2)))/4.d0
      f_t221 =(2.d0*f_y1*sqrt(f_y1*f_y1+f_x2*f_x2) + 
     &         2.d0*f_x2*f_x2*flog(f_y1+sqrt(f_y1*f_y1+f_x2*f_x2)))/4.d0
      f_t311 =(2.d0*f_x1*sqrt(f_x1*f_x1+f_y1*f_y1) + 
     &         2.d0*f_y1*f_y1*flog(f_x1+sqrt(f_x1*f_x1+f_y1*f_y1)))/4.d0
      f_t312 =(2.d0*f_x2*sqrt(f_x2*f_x2+f_y1*f_y1) + 
     &         2.d0*f_y1*f_y1*flog(f_x2+sqrt(f_x2*f_x2+f_y1*f_y1)))/4.d0
      f_t322 =(2.d0*f_x2*sqrt(f_x2*f_x2+f_y2*f_y2) + 
     &         2.d0*f_y2*f_y2*flog(f_x2+sqrt(f_x2*f_x2+f_y2*f_y2)))/4.d0
      f_t321 =(2.d0*f_x1*sqrt(f_x1*f_x1+f_y2*f_y2) + 
     &         2.d0*f_y2*f_y2*flog(f_x1+sqrt(f_x1*f_x1+f_y2*f_y2)))/4.d0
      f_t411 = ((f_x1*f_x1+f_y1*f_y1)**(1.5))/3.d0
      f_t412 = ((f_x1*f_x1+f_y2*f_y2)**(1.5))/3.d0
      f_t422 = ((f_x2*f_x2+f_y2*f_y2)**(1.5))/3.d0
      f_t421 = ((f_x2*f_x2+f_y1*f_y1)**(1.5))/3.d0

      f_t1 = f_t111-f_t112-f_t121+f_t122
      f_t2 = f_t211-f_t212-f_t221+f_t222
      f_t3 = f_t311-f_t312-f_t321+f_t322
      f_t4 = f_t411-f_t412-f_t421+f_t422

      f_Vcoeffs(1) = (f_t1*f_x2*f_y2 + f_t2*-f_y2 + 
     &                f_t3*-f_x2 + f_t4)/((f_x2-f_x1)*(f_y2-f_y1))
      f_Vcoeffs(2) = (f_t1*-f_x1*f_y2 + f_t2*f_y2 + 
     &                f_t3*f_x1 - f_t4)/((f_x2-f_x1)*(f_y2-f_y1))
      f_Vcoeffs(3) = (f_t1*f_x1*f_y1 + f_t2*-f_y1 + 
     &                f_t3*-f_x1 + f_t4)/((f_x2-f_x1)*(f_y2-f_y1))
      f_Vcoeffs(4) = (f_t1*-f_x2*f_y1 + f_t2*f_y1 + 
     &                f_t3*f_x2 - f_t4)/((f_x2-f_x1)*(f_y2-f_y1))

      RETURN

	END

c======================================================================
      SUBROUTINE calc_electrostatic
c======================================================================

c----------------------------------------------------------------------
c	Electorostatic force on slider
c----------------------------------------------------------------------

c     Shared data
c     ===========
      USE siml_pram
      USE syst_cnfg
      USE sldr_grid
      USE sldr_arys
      USE sldr_dynm
      USE elec_stat

      IMPLICIT NONE

c     Local variables
c     ===============
      INTEGER i,j,im1,ip1,jm1,jp1
      INTEGER nxxm1,nyym1
      DOUBLE PRECISION dxdy,delx,dely
      DOUBLE PRECISION tempgap,elecst
      DOUBLE PRECISION xelecst,yelecst

c     Disk Stationary
c     ===============      
      IF (omega .EQ. 0.0d0) RETURN

      IF (igp1.eq.1) THEN
        nxxm1=nxm1-1
        nyym1=nym1-1
      ELSE
        nxxm1=nxm1
        nyym1=nym1
      ENDIF
      
      felecst = 0.d0
      telecst = 0.d0
  	relecst = 0.d0
  	
  	IF (elecpot .EQ. 0.0d0) RETURN
  	
  	
  	DO j=2,nyym1
	  jm1=j-1
	  jp1=j+1
	  dely=(yref(jp1)-yref(jm1))/2.d0
	  
	  DO i=2,nxxm1
	      im1=i-1
	      ip1=i+1
	  
	      tempgap = hm*(hnew(i,j))
	  
c           Unbounded forces below 0.2nm
c           ============================
	      IF (tempgap .LT. 0.2d-9) THEN 
	          tempgap = 0.2d-9
	      ENDIF
	  
	      delx=(xref(ip1)-xref(im1))/2.d0
	      dxdy=delx*dely

c           Electrostatic
c           =============
		    elecst=eps0*Ke*(elecpot**2)*dxdy/2.d0/(tempgap**2)
		    
	      felecst = felecst - elecst
	      telecst = telecst - elecst * xref(i)
	      relecst = relecst - elecst * (yref(j) - 0.5d0*yl)
	      
	  ENDDO
      ENDDO
     
      IF (felecst .ne. 0.0) THEN
        xelecst = telecst / felecst !x and y normalized centers of force
        yelecst = relecst / felecst
      ELSE
        xelecst = 0.0
        yelecst = 0.0
      ENDIF
      
c     Electrostatic Force in Z dir and Pitch and Roll Torques due to it
c     =================================================================
      felecst =  felecst * elecnormf
	telecst = (felecst * (xelecst - xg)) * xl !unnormalize torque
      relecst = (felecst * (yelecst - yg)) * xl
      
      fn = fn + felecst

      RETURN
	END
	
c======================================================================