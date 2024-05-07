c**********************************************************************
c     Subroutines in this file :
c     1. adaptive
c     2. adpt_grid
c     3. interp
c     4. interp1
c     5. get_grids
c     6. snapgrid
c     7. get_toth
c     8. get_floor
c     9. gethx1
c     10.gethx
c**********************************************************************

c======================================================================
	SUBROUTINE adaptive
c======================================================================

c     Shared Data
c     ===========
	USE siml_pram
	USE syst_cnfg
	USE rail_data
	USE grid_ctrl
	USE sldr_grid
	USE sldr_arys
	USE sldr_rghn
	USE disk_rghn
	USE trck_prfl
	USE rynl_data
	USE sldr_dynm
	USE aspr_data
	USE POI_cords
	
	IMPLICIT REAL*8(a-h,o-z)
	
      INTEGER, PARAMETER :: max_nx = 1010 
      INTEGER, PARAMETER :: max_ny = 1010
	
c     Local Variable
c     ==============	
	LOGICAL crash
      
      DOUBLE PRECISION dnx(1001),dny(1001)
      DOUBLE PRECISION xtemp(1001),ytemp(1001)
      DOUBLE PRECISION xrefn(max_nx),yrefn(max_ny)
      
      DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: ptemp

            
c	Use the same mesh size
c	======================
	ngrad=1001
	nxn=nx
	nyn=ny
	
c     Allocate
c     ========
	ALLOCATE(ptemp(nx,ny))
	ptemp = 0.0d0

c	Maximum difference between highest and lowest gradient allowed
c		difmax=80.d0
c	Exponential decay factor in smoothing the gradient change
c		decay=60.d0
c	Amplitude for grid density at taper end xt
c	==============================================================
	axt=2.d0

c	Exponential decay rate for xt grid density
c	==========================================	
	dxt=100.d0
c
	dnxx=0.001d0
	dnxmin=1.d20
	dnyx=0.001d0
	dnymin=1.d20

19	CONTINUE

c	calculate average gradient in x direction
c	=========================================

	DO i=1,nx

c	Avoid zero gradient
c	===================
		dnx(i)=0.001d0

		DO j=2,ny-1
			
			IF(ipmax.EQ.1)THEN
				
				IF(i.EQ.nx)THEN
					dnx(i)=DMAX1(dnx(i),+DABS(p(i,j)-p(i-1,j)))
				ELSE IF(i.EQ.1) THEN
					dnx(i)=DMAX1(dnx(i),+DABS(p(i+1,j)-p(i,j)))
				ELSE
					dnx(i)=DMAX1(dnx(i),+DABS(p(i+1,j)-p(i-1,j)))
				ENDIF

			ELSE

				IF(i.EQ.nx)THEN
					dnx(i)=dnx(i)+DABS(p(i,j)-p(i-1,j))
				ELSE IF(i.EQ.1) THEN
					dnx(i)=dnx(i)+DABS(p(i+1,j)-p(i,j))
				ELSE
					dnx(i)=dnx(i)+DABS(p(i+1,j)-p(i-1,j))/2.d0
				ENDIF

			ENDIF

		ENDDO

		IF(i.EQ.nx)THEN
			dnx(i)=dnx(i)/(xref(i)-xref(i-1))
		ELSE IF(i.EQ.1) THEN
			dnx(i)=dnx(i)/(xref(i+1)-xref(i))
		ELSE
			dnx(i)=dnx(i)/(xref(i+1)-xref(i-1))
		ENDIF

c	Maximum and minimum gradient
c	============================
		IF(dnxx.LT.dnx(i))dnxx=dnx(i)
		IF(dnxmin.GT.dnx(i))dnxmin=dnx(i)
	ENDDO
c
c	Calculate average pressure gradient in y direction
c	==================================================
	DO j=1,ny

c	Avoid zero gradient
c	===================
		dny(j)=0.001d0

		DO i=2,nx-1
			IF(ipmax.EQ.1)THEN
				
				IF(j.EQ.ny) THEN
					dny(j)=DMAX1(dny(j),DABS(p(i,j)-p(i,j-1)))
				ELSE IF(j.EQ.1)THEN
					dny(j)=DMAX1(dny(j),DABS(p(i,j+1)-p(i,j)))
				ELSE
					dny(j)=DMAX1(dny(j),DABS(p(i,j+1)-p(i,j-1)))
				ENDIF

			ELSE

				IF(j.EQ.ny) THEN
						dny(j)=dny(j)+DABS(p(i,j)-p(i,j-1))
					ELSE IF(j.EQ.1)THEN
						dny(j)=dny(j)+DABS(p(i,j+1)-p(i,j))
					ELSE
						dny(j)=dny(j)+DABS(p(i,j+1)-p(i,j-1))/2.d0
					ENDIF
			ENDIF

		ENDDO

		IF(j.EQ.ny)THEN
			dny(j)=dny(j)/(yref(j)-yref(j-1))
		ELSE IF(j.EQ.1) THEN
			dny(j)=dny(j)/(yref(j+1)-yref(j))
		ELSE
			dny(j)=dny(j)/(yref(j+1)-yref(j-1))
		ENDIF

c	Maximum and minimum gradient
c	============================
		IF(dnyx.LT.dny(j))dnyx=dny(j)
		IF(dnymin.GT.dny(j))dnymin=dny(j)

	ENDDO

c	Limit the difference between maximum and minimum gradient
c	=========================================================
	DO i=1,nx
		dnx(i)=DMAX1(dnx(i),dnxx/difmax)
	ENDDO

	DO j=1,ny
		dny(j)=DMAX1(dny(j),dnyx/difmax)
	ENDDO

	DO i=1,ngrad
		xtemp(i)=dfloat(i-1)*(1.d0/dfloat(ngrad-1))
		ytemp(i)=xtemp(i)*yl
	ENDDO

c	Interpolate dnx and dny to xtemp and ytemp grids
c	================================================
	CALL interp1(ngrad,nx,ngrad,xref,xtemp,dnx)
	CALL interp1(ngrad,ny,ngrad,yref,ytemp,dny)

c	Create adaptive grids in x and y directions
c	===========================================
	CALL adpt_grid(ngrad,nyn,yl,ytemp,dny,yrefn,
     &			   decay,0.d0,0.d0,0.d0)
     
	IF(xt.GT.1.d-10)THEN
		axt=2.d0
		ienough=0
		iround=1
		DO WHILE(ienough.EQ.0.AND.iround.LE.8)
			CALL adpt_grid(ngrad,nxn,1.d0,xtemp,dnx,
     &				   xrefn,decay,xt,axt,dxt)
     
c	      Interpolating pressure to the new grid
c	      ======================================
            ptemp(1:nx,1:ny) = p(1:nx,1:ny)			
			CALL interp(nx,nxn,ny,nyn,xref,xrefn,yref,yrefn,ptemp)
			
			ind=1
			DO WHILE(.not.(xrefn(ind).LT.xt.AND.xrefn(ind+1).GE.xt))
				ind=ind+1
			ENDDO

			ninc=nx/12
			IF(ind.LT.dfloat(ninc)) GOTO 1

			pmaxind=1.d0
			pmaxi=1.d0
			DO j=1,ny
				pmaxind=DMAX1(pmaxind,ptemp(ind+1,j))
				pmaxi=DMAX1(pmaxi,ptemp(ind-ninc+2,j))
			ENDDO

			IF((pmaxi-1.d0).GE.0.2d0*(pmaxind-1.d0))ienough=1
			iround=iround+1
1			axt=axt+4.d0
		ENDDO
		axt=axt-4.d0
		p(1:nx,1:ny) = ptemp(1:nx,1:ny)
		
	ELSE
	
		CALL adpt_grid(ngrad,nxn,1.d0,xtemp,dnx,
     &					xrefn,decay,0.d0,0.d0,dxt)
		CALL interp(nx,nxn,ny,nyn,xref,xrefn,yref,yrefn,p)
          
	ENDIF

c	Create corresponding adaptive dense mesh
c	use the old names for the grid coordinates
c	==========================================
	DO i=1,nxn
		xref(i)=xrefn(i)
	ENDDO
	DO j=1,nyn
		yref(j)=yrefn(j)
	ENDDO
	nx=nxn
	ny=nyn

	CALL snapgrid
	CALL calc_bearing_number
	CALL gethx1
	CALL gethx

c     Reset the time : why? 
c     ==============
	t=0.d0
	CALL get_toth(crash,ic,jc)
	IF(crash) THEN
		WRITE(*,250)xref(ic)*xl,yref(jc)*xl
		STOP
	ENDIF

c	Resolve the pressure distribution
c	=================================
	CALL reyneq(1)

250	FORMAT(5X,'Initial FH: crash at (',E16.9,'m;',E16.9,'m)')
	
c     Deallocate
c     ==========
      DEALLOCATE(ptemp)
      	
	RETURN
	END

c======================================================================
	SUBROUTINE adpt_grid(nx,nxp,xl,x,dnp,xp,dec,xt,axt,dxt)
c======================================================================
	IMPLICIT REAL*8(a-h,o-z)
	REAL*8  n(1001),dn(1001),dx(1001),x(nx),dnp(nx),xp(nxp)

	DO k=1,nx-1
		dx(k)=x(k+1)-x(k)
	ENDDO
	
	DO k=1,nx
		dn(k)=0.d0
		DO j=1,nx
			dn(k)=dn(k)+dnp(j)*dexp(-DABS(x(k)-x(j))*dec)
		ENDDO
	ENDDO

c	Calculate scaling constant for dn(k)
c	===================================
	const=0.d0
	DO k=1,nx-1
		const=const+dx(k)*(dn(k)+dn(k+1))/2.d0
	ENDDO

	DO k=1,nx
		dn(k)=dn(k)+const*axt*dexp(-DABS(x(k)-xt)*dxt)
	ENDDO

	const=0.d0
	DO k=1,nx-1
		const=const+dx(k)*(dn(k)+dn(k+1))/2.d0
	ENDDO
	const=dfloat(nxp-1)/const

c	Scaling dn(k) with const
c	========================
	DO k=1,nx
		dn(k)=dn(k)*const
	ENDDO

c	Calculate fractional grid number at old grid points
c	===================================================
	n(1)=1.d0
	n(nx)=dfloat(nxp)
	DO k=2,nx
		n(k)=n(k-1)+dx(k-1)*(dn(k-1)+dn(k))/2.d0
	ENDDO

c	Calculate new grid location
c	============================
	DO k=1,nx-1
		a=(dn(k+1)-dn(k))/2.d0/dx(k)
		b=(x(k+1)*dn(k)-x(k)*dn(k+1))/2.d0/dx(k)
	   c=n(k)-x(k)*((2.d0*x(k+1)-x(k))*dn(k)-x(k)*dn(k+1))/2.d0/dx(k)
		DO j=int(n(k)+1),int(n(k+1))
			IF(DABS(a).GT.DABS(b))THEN
				xp(j)=(dsqrt(b**2-a*(c-dfloat(j)))-b)/a
			ELSE
				xp(j)=(dfloat(j)-c)/(dsqrt(b**2-a*(c-dfloat(j)))+b)
			ENDIF
		ENDDO
	ENDDO
	xp(1)=0.d0
	xp(nxp)=xl
	
	RETURN
	END

c======================================================================
	SUBROUTINE interp(nx,nxn,ny,nyn,x,xn,y,yn,p)
c======================================================================

	IMPLICIT NONE

c     Arguments
c     =========
      INTEGER n1,n2
      INTEGER nx,nxn
      INTEGER ny,nyn

      DOUBLE PRECISION p(nx,ny)
	DOUBLE PRECISION x(nx),y(ny)
	DOUBLE PRECISION xn(nxn),yn(nyn)

c     Local Variables
c     ===============
      INTEGER i, j, ix, jy
      INTEGER         , DIMENSION (:)  , ALLOCATABLE :: ixn,jyn
      DOUBLE PRECISION, DIMENSION (:)  , ALLOCATABLE :: cox,coy
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: pn
      
      ALLOCATE(ixn(nx))
      ALLOCATE(jyn(ny))
      ALLOCATE(cox(nx))
      ALLOCATE(coy(ny))
      ALLOCATE(pn(nx,ny))
      
      ixn = 0
      jyn = 0
      cox = 0.0d0
      coy = 0.0d0
      pn  = 0.0d0 

c	Find the neighbouring grid index
c	================================
	ixn(1)=1
	ixn(nxn)=nx-1
	
	DO i=2,nxn-1
		ix=ixn(i-1)
		DO WHILE(.not.(x(ix).LE.xn(i).AND.
     &		x(ix+1).GT.xn(i)))
			ix=ix+1
		ENDDO
		ixn(i)=ix
c	Interpolation coefficient
c	=========================
		cox(i)=(x(ix+1)-xn(i))/(x(ix+1)-x(ix))
	ENDDO
		
	cox(1)=1.d0
	cox(nxn)=0.d0
	jyn(1)=1
	jyn(nyn)=ny-1
	
	DO j=2,nyn-1
		jy=jyn(j-1)
		DO WHILE(.NOT.(y(jy).LE.yn(j).AND.
     &		y(jy+1).GT.yn(j)))
			jy=jy+1
		ENDDO
		jyn(j)=jy
		coy(j)=(y(jy+1)-yn(j))/(y(jy+1)-y(jy))
	ENDDO

	coy(1)  = 1.d0
	coy(nyn)= 0.d0

	DO i=1,nxn
		DO j=1,nyn
			pn(i,j)=cox(i)*(p(ixn(i),jyn(j))*coy(j)+
     &			p(ixn(i),jyn(j)+1)*(1.d0-coy(j)))+(1.d0-cox(i))
     &			*(p(ixn(i)+1,jyn(j))*coy(j)+
     &			p(ixn(i)+1,jyn(j)+1)*(1.d0-coy(j)))
		ENDDO
	ENDDO

	DO i=1,nxn
		DO j=1,nyn
			p(i,j)=pn(i,j)
		ENDDO
	ENDDO
	
	DEALLOCATE(ixn)
      DEALLOCATE(jyn)
      DEALLOCATE(cox)
      DEALLOCATE(coy)
      DEALLOCATE(pn)
	
	RETURN
	
	END

c======================================================================
	SUBROUTINE interp1(n1,nx,nxn,x,xn,p)
c======================================================================


c	Find the neighbouring grid index
c	================================

	IMPLICIT REAL*8(a-h,o-z)
	
	DOUBLE PRECISION cox(1001),pn(1001),p(n1)
	DOUBLE PRECISION x(nx),xn(nxn),ixn(1001)

	ixn(1)=1
	ixn(nxn)=nx-1
	DO i=2,nxn-1
		ix=ixn(i-1)
		DO WHILE(.not.(x(ix).LE.xn(i).AND.
     &		x(ix+1).GT.xn(i)))
		  ix=ix+1
		ENDDO
		ixn(i)=ix

c	Interpolation coefficient
c	=========================
		cox(i)=(x(ix+1)-xn(i))/(x(ix+1)-x(ix))
	ENDDO

	cox(1)=1.d0
	cox(nxn)=0.d0

	DO i=1,nxn
		pn(i)=cox(i)*p(ixn(i))+(1.d0-cox(i))*p(ixn(i)+1)
	ENDDO

	DO i=1,nxn
		p(i)=pn(i)
	ENDDO

	RETURN
	END

c======================================================================
	SUBROUTINE get_grids
c======================================================================

c----------------------------------------------------------------------
c     Generate initial grid : this is not adaptive to pressure
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE siml_pram
      USE syst_cnfg
      USE rail_data
      USE grid_ctrl
      USE sldr_grid
	
	IMPLICIT REAL*8(a-h,o-z)
	
	INTEGER, PARAMETER :: max_nx = 1010 
      INTEGER, PARAMETER :: max_ny = 1010

	DOUBLE PRECISION dxref(max_nx),dyref(max_ny)

c	Boundary grid information
c	=========================
	xnt(1)=0.d0
	nxt(1)=1
	xnt(nsx + 1)= 1.d0
	nxt(nsx + 1)= nx

	ynt(1)=0.d0
	nyt(1)=1
	IF(isymmetry.EQ.1) THEN

c	For symmetric y grid, calculate half first
c	==========================================
		ynt(nsy+1)=yl/2
		nyt(nsy+1)=int(ny/2)+1
	ELSE
		ynt(nsy + 1)= yl
		nyt(nsy + 1)= ny
	ENDIF

c	Start calculate x grid
c	======================
	DO i=1,nsx
		xref(nxt(i))=xnt(i)
		dxref(nxt(i))=1.d0

		DO k=nxt(i)+1,nxt(i+1)
			dxref(k)=dxref(k-1)*dxr(i)
			xref(k)=xref(k-1)+dxref(k-1)
		ENDDO
	
c	Calculate the constant for adjusting the grid step length
c	=========================================================
		const=(xnt(i+1)-xnt(i))/(xref(nxt(i+1))-xnt(i))

		DO k=nxt(i)+1,nxt(i+1)
		dxref(k-1)=dxref(k-1)*const
		xref(k)=xref(k-1)+dxref(k-1)
		ENDDO
	ENDDO

c	Start calculate y grid
c	======================
	DO j=1,nsy
		yref(nyt(j))=ynt(j)
		dyref(nyt(j))=1.d0

		DO k=nyt(j)+1,nyt(j+1)
			dyref(k)=dyref(k-1)*dyr(j)
			yref(k)=yref(k-1)+dyref(k-1)
		ENDDO

		IF(isymmetry.EQ.1.AND.j.EQ.nsy.AND.mod(ny,2).EQ.0)THEN

c	For symmetric y grid and even number of total y grid
c	point, half grid at center
c	=========================================================
			const=(ynt(j+1)-ynt(j))/(yref(nyt(j+1))-ynt(j)
     &			-dyref(nyt(nsy+1)-1)/2.d0)
		ELSE
			const=(ynt(j+1)-ynt(j))/(yref(nyt(j+1))-ynt(j))
		ENDIF

		DO k=nyt(j)+1,nyt(j+1)
			dyref(k-1)=dyref(k-1)*const
			yref(k)=yref(k-1)+dyref(k-1)
		ENDDO
	ENDDO

c	For symmetric y grid, generate the other half
c	=============================================
	IF(isymmetry.EQ.1)THEN
		DO k=ny,int(ny/2)+1,-1
			yref(k)=yl-yref(ny-k+1)
		ENDDO
	ENDIF

	RETURN
	END

c======================================================================
      SUBROUTINE snapgrid
c======================================================================
      
c     Shared Data
c     ===========
      USE syst_cnfg
      USE rail_data
      USE sldr_grid
      
	IMPLICIT REAL*8(a-h,o-z)

	INTEGER, PARAMETER :: max_nx = 1010 
      INTEGER, PARAMETER :: max_ny = 1010
      	
      DOUBLE PRECISION xtemp(max_nx),ytemp(max_ny)
      DOUBLE PRECISION xcontrol(max_nx),ycontrol(max_ny)

	DO i=1,nx
		xtemp(i)=xref(i)
	ENDDO

	DO i=1,ny
		ytemp(i)=yref(i)
	ENDDO

	ipx=0
	ipy=0

	DO k=1,num_rails
		DO i=1, npoints(k)
			IF(xrail(i,k).GT.0.d0.AND.xrail(i,k).LT.1.d0)THEN
				ipx=ipx+1
				xcontrol(ipx)=xrail(i,k)
			ENDIF
			IF(yrail(i,k).GT.0.d0.AND.yrail(i,k).LT.yl)THEN
				ipy=ipy+1
				ycontrol(ipy)=yrail(i,k)
			ENDIF
		ENDDO
	ENDDO

	IF(xt.GT.0.d0)THEN
		ipx=ipx+1
		xcontrol(ipx)=xt
	ENDIF

	IF(ipx.GT.1) CALL sort(1.d-3,xcontrol,ipx)
	IF(ipy.GT.1) CALL sort(1.d-3,ycontrol,ipy)

	IF(ipx.LE.(nx-2))THEN
		index=2
		DO i=1,ipx
			indexold=index
			IF(xref(index).GT.xcontrol(i))THEN
				xref(index)=xcontrol(i)
				index=index+1
			ELSE
				DO ip=index,nx-1
					IF(xcontrol(i).GE.xref(ip).AND.
     &					xcontrol(i).LT.xref(ip+1))THEN
						IF(DABS(xcontrol(i)-xref(ip)).LE
     &						.DABS(xcontrol(i)-xref(ip+1)))THEN
							index=ip
						ELSE
							index=ip+1
						ENDIF
						index=min(nx-1-ipx+i,index)
						xref(index)=xcontrol(i)
						index=index+1
						GOTO 10
					ENDIF
				ENDDO	!DO ip=index,nx-1
10				CONTINUE
			ENDIF
			co=(xref(index-1)-xref(indexold-1))/
     &			(xtemp(index-1)-xtemp(indexold-1))
			DO ik=indexold,index-2
				xref(ik)=xref(indexold-1)+co*(xref(ik)
     &				-xtemp(indexold-1))
			ENDDO

		ENDDO	!DO i=1,ipx

		co=(1.d0-xref(index-1))/
     &		(1.d0-xtemp(index-1))

		DO ik=index,nx-1
			xref(ik)=xref(index-1)+co*(xref(ik)-xtemp(index-1))
		ENDDO

	ENDIF


	IF(ipy.LE.(ny-2))THEN
		index=2
		DO i=1,ipy
			indexold=index
			IF(yref(index).GT.ycontrol(i))THEN
				yref(index)=ycontrol(i)
				index=index+1
			ELSE
				DO ip=index,ny-1
					IF(ycontrol(i).GE.yref(ip).AND.
     &					ycontrol(i).LT.yref(ip+1))THEN
						IF(DABS(ycontrol(i)-yref(ip)).LE
     &						.DABS(ycontrol(i)-yref(ip+1)))THEN
							index=ip
						ELSE
							index=ip+1
						ENDIF

						index=min(ny-1-ipy+i,index)
						yref(index)=ycontrol(i)
						index=index+1
						GOTO 20
					ENDIF
				ENDDO
20			CONTINUE
			ENDIF
			co=(yref(index-1)-yref(indexold-1))/
     &			(ytemp(index-1)-ytemp(indexold-1))
			DO ik=indexold,index-2
				yref(ik)=yref(indexold-1)+co*(yref(ik)-
     &				ytemp(indexold-1))
			ENDDO
		ENDDO
		co=(yl-yref(index-1))/
     &		(yl-ytemp(index-1))
		DO ik=index,ny-1
			yref(ik)=yref(index-1)+co*(yref(ik)-ytemp(index-1))
		ENDDO
	ENDIF

	RETURN
	END

c======================================================================
	SUBROUTINE get_toth(crash,ic,jc)
c======================================================================

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
      
	IMPLICIT REAL*8(a-h,o-z)
	
c     Local Variable
c     ==============	
	LOGICAL crash
	
	crash=.FALSE.
	ssk=dsin(ske)
	csk=dcos(ske)

	DO i=1,nx
		xtemp=xg-xref(i)
		DO j=1,ny
			ytemp=yref(j)-(0.5d0*yl+yg)
			xloc=dfx+(xtemp*csk-ytemp*ssk)
			yloc=dfy+(xtemp*ssk+ytemp*csk)

			CALL get_floor(xloc,yloc,hloc)
			hfloor(i,j)= hloc
		ENDDO
	ENDDO

	DO i=2,nxm1
		DO j=2,nym1
			hxlnew(i,j)=hxlnew(i,j)+0.5d0*(hfloor(i,j)-hsad(i,j)+
     &			hfloor(i,j-1)-hsad(i,j-1))
			hxrnew(i,j)=hxrnew(i,j)+0.5d0*(hfloor(i,j)-hsad(i,j)+
     &			hfloor(i,j-1)-hsad(i,j-1))
			IF(j.EQ.nym1) THEN
				hxlnew(i,j+1)=hxlnew(i,j+1)+0.5d0*(hfloor(i,j)-
     &				hsad(i,j)+hfloor(i,j+1)-hsad(i,j+1))		
				hxrnew(i,j+1)=hxrnew(i,j+1)+0.5d0*(hfloor(i,j)-
     &				hsad(i,j)+hfloor(i,j+1)-hsad(i,j+1))
				IF(hxlnew(i,j+1).LE.0.d0.OR.
     &				hxrnew(i,j+1).LE.0.d0) THEN
						crash=.TRUE.
						ic=i
						jc=j
				ENDIF
			ENDIF
				 
			hydnew(i,j)=hydnew(i,j)+0.5d0*(hfloor(i,j)-hsad(i,j)+
     &			hfloor(i-1,j)-hsad(i-1,j))		
			hyunew(i,j)=hyunew(i,j)+0.5d0*(hfloor(i,j)-hsad(i,j)+
     &			hfloor(i-1,j)-hsad(i-1,j))
		  
			IF(i.EQ.nxm1) THEN
				hydnew(i+1,j)=hydnew(i+1,j)+0.5d0*(hfloor(i,j)-
     &				hsad(i,j)+hfloor(i+1,j)-hsad(i+1,j))		
				hyunew(i+1,j)=hyunew(i+1,j)+0.5d0*(hfloor(i,j)-
     &				hsad(i,j)+hfloor(i+1,j)-hsad(i+1,j))
				IF(hydnew(i+1,j).LE.0.d0.OR.
     &				hyunew(i+1,j).LE.0.d0) THEN
					crash=.TRUE.
					ic=i
					jc=j
				ENDIF
			ENDIF 	 
			hnew(i,j)=hnew(i,j)+hfloor(i,j)-hsad(i,j)
			IF(hxlnew(i,j).LE.0.d0.OR.hxrnew(i,j).LE.0.d0.OR.
     &			hydnew(i,j).LE.0.d0.OR.hyunew(i,j).LE.0.d0.OR.
     &			hnew(i,j).LE.0.d0) THEN
					crash=.TRUE.
					ic=i
					jc=j
			ENDIF
		ENDDO
	ENDDO

	RETURN

38	crash=.TRUE.
	ic=i
	jc=j
	
	DO i=1,nx
		DO j=1,ny
			IF ((i.EQ.1.d0).OR.(j.EQ.1.d0).OR.
     &			(i.EQ.ic.AND.j.GT.jc).OR.(i.GT.ic))
     &			hnew(i,j)=hnew(i,j)+hfloor(i,j)-hsad(i,j)
		ENDDO
	ENDDO
999	RETURN
	END
	  
c======================================================================
	SUBROUTINE get_floor(xloc,yloc,hloc)
c======================================================================

c----------------------------------------------------------------------
c     Compute the floor height (hloc) at a given (xloc,yloc) on slider
c----------------------------------------------------------------------

c     Shared Data
c     ===========	
	USE siml_pram
      USE syst_cnfg
	USE rail_data
	USE sldr_grid
	USE sldr_arys
	USE disk_rghn
	USE disk_fltr
	USE trck_prfl
	USE rynl_data
	USE sldr_dynm
	
	IMPLICIT REAL*8(a-h,o-z)
	
c     Local Variable
c     ==============	
	LOGICAL crash

c	Initialize
c	==========
	hloc = 0.d0

c	thh represent the elevation from the mean disk surface
c	======================================================
	IF((nf_wave.NE.0).OR.(nf_asper.NE.0).OR.(nf_zone.NE.0)) THEN
	  CALL fl_an(xloc,yloc,thh)
		hloc = hloc-thh
	ENDIF

c	Track profile
c	=============
	IF(inwave.NE.0) THEN
		xloc=dfinit+xloc
		xloc=DMOD(xloc,xfl)
		IF(xloc.LE.0.d0) xloc=xloc+xfl
		IF(xloc.GE.xfl) xloc=xloc-xfl
		CALL locate(xfref,nfx,xloc,iloc)
		ai =(xloc-xfref(iloc+1))/(xfref(iloc)-xfref(iloc+1))
		ai1=(xloc-xfref(iloc))/(xfref(iloc+1)-xfref(iloc))
		hloc=hloc-(ai*hfnew(iloc)+ai1*hfnew(iloc+1))
	ENDIF

c	Sinusoidal disk flutter
c	=======================
	IF(kft.EQ.0) RETURN

c     Substracting elevation due to disk flutter
c     ==========================================
	IF(t.GT.ft_ts.AND.t.LT.ft_tf) THEN
		hloc = hloc - ft_mag*dsin(ft_omg*(t-ft_ts))
	ENDIF

	RETURN
	END
	
c======================================================================
	SUBROUTINE gethx1
c======================================================================

c----------------------------------------------------------------------
c     Computes hxy(i,j), hxr(i,j), hxl(i,j), zta(i,j), eta(i,j)
c     Uses : xref(i) and yref(i)
c----------------------------------------------------------------------

c     Shared Data
c     ===========  
      USE siml_pram
      USE syst_cnfg
      USE rail_data
      USE sldr_grid
      USE sldr_arys
      USE sldr_rghn
      USE rynl_data
      USE sldr_dynm

	IMPLICIT REAL*8(a-h,o-z)

c     Local Variable
c     ==============
	DOUBLE PRECISION dsx(30),dsy(30)
	nav=20

	DO i=2,nx-1
		delx=(xref(i+1)-xref(i-1))/2.d0/(nav-1)
		DO ks=1,nav
			dsx(ks)=(xref(i-1)+xref(i))/2.d0+(ks-1)*delx
		ENDDO

		DO j=2,ny-1
			dely=(yref(j+1)-yref(j-1))/2.d0/(nav-1)
			DO ks=1,nav
				dsy(ks)=(yref(j-1)+yref(j))/2.d0+(ks-1)*dely
			ENDDO	
			hav=0.d0

			havx=0.d0
			hxmax=-1.0d+9
			hxmin=1.0d+9

			havy=0.d0
			hymax=-1.0d+9
			hymin=1.0d+9

			havsx=0.d0
			hxsmax=-1.0d+9
			hxsmin=1.0d+9

			havsy=0.d0
			hysmax=-1.0d+9
			hysmin=1.0d+9

			DO ksx=1,nav
				xv=dsx(ksx)
				DO ksy=1,nav
					yv=dsy(ksy)

					IF(ksy.EQ.1) THEN                     
                        CALL pointrecess(xv,yv,frecess)
						thh=frecess
						havx=havx+thh
						IF(thh.GT.hxmax) hxmax=thh
						IF(thh.LT.hxmin) hxmin=thh
					ENDIF

					IF(j.EQ.(ny-1).AND.ksy.EQ.nav) THEN
                        CALL pointrecess(xv,yv,frecess)
						thh=frecess
						havsx=havsx+thh
						IF(thh.GT.hxsmax) hxsmax=thh
						IF(thh.LT.hxsmin) hxsmin=thh
					ENDIF

					IF(ksx.EQ.1) THEN
                        CALL pointrecess(xv,yv,frecess)
						thh=frecess
						havy=havy+thh
				
						IF(thh.GT.hymax) hymax=thh
						IF(thh.LT.hymin) hymin=thh
						
					ENDIF

					IF(i.EQ.(nx-1).AND.ksx.EQ.nav) THEN
                        CALL pointrecess(xv,yv,frecess)
						thh=frecess
						havsy=havsy+thh
					
						IF(thh.GT.hysmax) hysmax=thh
						IF(thh.LT.hysmin) hysmin=thh
						
					ENDIF
				ENDDO
			ENDDO

			havx=havx/nav
			hxl(i,j)=hxmin
			hxr(i,j)=hxmax
			
			IF(DABS(hxmax-hxmin).LT.epsilon) THEN
				zta(i,j)=1.d0
			ELSE
				zta(i,j)=(hxmax-havx)/(hxmax-hxmin)
			ENDIF

			IF(j.EQ.(ny-1)) THEN
				havsx=havsx/nav
				hxl(i,j+1)=hxsmin
				hxr(i,j+1)=hxsmax
				IF(DABS(hxsmax-hxsmin).LT.epsilon) THEN
					zta(i,j+1)=1.d0
				ELSE
					zta(i,j+1)=(hxsmax-havsx)/(hxsmax-hxsmin)
				ENDIF
			ENDIF

			havy=havy/nav
			hyd(i,j)=hymin
			hyu(i,j)=hymax
			
			IF(DABS(hymax-hymin).LT.epsilon) THEN
				eta(i,j)=1.d0
			ELSE
				eta(i,j)=(hymax-havy)/(hymax-hymin)
			ENDIF

			IF(i.EQ.(nx-1)) THEN
				havsy=havsy/nav
				hyd(i+1,j)=hysmin
				hyu(i+1,j)=hysmax
				IF(DABS(hysmax-hysmin).LT.epsilon) THEN
					eta(i+1,j)=1.d0
				ELSE
					eta(i+1,j)=(hysmax-havsy)/(hysmax-hysmin)
				ENDIF
			ENDIF
 		ENDDO
	ENDDO

	DO i=1,nx
		DO j=1,ny
			xv=xref(i)
			yv=yref(j)
            CALL pointrecess(xv,yv,thh)
			hxy(i,j)=thh
		ENDDO
	ENDDO

	RETURN
	END

c======================================================================
	SUBROUTINE gethx
c======================================================================

c----------------------------------------------------------------------
c     This subroutine calculates separation hmin(i,j) given slider zc,
c	hp and hx0 and add base recess hxy(i,j)
c     Output : hnew(i,j)
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE siml_pram
      USE syst_cnfg
      USE rail_data
      USE sldr_grid
      USE sldr_arys
      USE sldr_dynm

	IMPLICIT REAL*8(a-h,o-z)

	DO i=2,nx-1
		delx=(xref(i+1)-xref(i-1))/2.d0
		DO j=2,ny-1
			
			dely=(yref(j+1)-yref(j-1))/2.d0

			xv=(xref(i-1)+xref(i))/2.d0+0.5d0*delx
			yv=(yref(j-1)+yref(j))/2.d0
			zx=hmin-hx0*(1.d0-xv)
			zy=hy*(yv-0.5d0*yl)
			hxlnew(i,j)=hxl(i,j)+zx+zy
			hxrnew(i,j)=hxr(i,j)+zx+zy

			IF(j.EQ.(ny-1)) THEN
				xv=(xref(i-1)+xref(i))/2.d0+0.5d0*delx
				yv=(yref(j)+yref(j+1))/2.d0
				zx=hmin-hx0*(1.d0-xv)
				zy=hy*(yv-0.5d0*yl)
				hxlnew(i,j+1)=hxl(i,j+1)+zx+zy
				hxrnew(i,j+1)=hxr(i,j+1)+zx+zy
			ENDIF

			xv=(xref(i-1)+xref(i))/2.d0
			yv=(yref(j-1)+yref(j))/2.d0+0.5d0*dely
			zx=hmin-hx0*(1.d0-xv)
			zy=hy*(yv-0.5d0*yl)
			hydnew(i,j)=hyd(i,j)+zx+zy
			hyunew(i,j)=hyu(i,j)+zx+zy


			IF(i.EQ.(nx-1)) THEN
				xv=(xref(i)+xref(i+1))/2.d0
				yv=(yref(j-1)+yref(j))/2.d0+0.5d0*dely
				zx=hmin-hx0*(1.d0-xv)
				zy=hy*(yv-0.5d0*yl)
				hydnew(i+1,j)=hyd(i+1,j)+zx+zy
				hyunew(i+1,j)=hyu(i+1,j)+zx+zy
			ENDIF
		  
			xv=xref(i)
			yv=yref(j)
			zx=hmin-hx0*(1.d0-xv)
			zy=hy*(yv-0.5d0*yl)
			hnew(i,j)=hxy(i,j)+zx+zy
		ENDDO
	ENDDO

	yv1=yref(1)
	zy1=hy*(yv1-0.5d0*yl)
	yv2=yref(ny)
	zy2=hy*(yv2-0.5d0*yl)
	
	DO i=1,nx
		xv=xref(i)
		zx=hmin-hx0*(1.d0-xv)
		hnew(i,1)=hxy(i,1)+zx+zy1
		hnew(i,ny)=hxy(i,ny)+zx+zy2
	ENDDO

	xv1=xref(1)
	zx1=hmin-hx0*(1.d0-xv1)
	xv2=xref(nx)
	zx2=hmin-hx0*(1.d0-xv2)

	DO j=1,ny
		yv=yref(j)
		zy=hy*(yv-0.5d0*yl)
		hnew(1,j)=hxy(1,j)+zx1+zy
		hnew(nx,j)=hxy(nx,j)+zx2+zy
	ENDDO

	RETURN
	END
	
c======================================================================