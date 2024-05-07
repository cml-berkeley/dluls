c**********************************************************************
c     Subroutines in this file :
c      1. read_rail
c      2. norm_rail
c      3. wall_profile
c      4. pointrecess
c      5. get_rail_size
c      6. allocate_rail_arrays
c      7. colinear
c      8. intersection
c      9. PNPOLY
c     10. PTONPOLYBOUNDARY
c**********************************************************************

c======================================================================
	SUBROUTINE read_rail
c======================================================================

c     Shared Data
c     ===========      
      USE syst_cnfg
      USE POI_cords
      USE rail_data
      USE wall_data

	IMPLICIT REAL*8(a-h,o-z)
	
	CALL get_rail_size
	
c     Open input file
c     ===============     
	OPEN(1,ERR=999,FILE='rail.dat',STATUS='OLD')
	READ(1,*)
	READ(1,*)
	READ(1,*)

c	Now for each rail
c	=================
	READ(1,*) num_rails, num_walls
	
c     Allocate rail and wall matrices
c     ===============================
	CALL allocate_rail_arrays
	
	DO k=1,num_rails
		READ(1,*) npoints(k),istep(k)

		DO j=1,npoints(k)
			READ(1,*) xrail(j,k),yrail(j,k),indexw(k,j)
		ENDDO

		IF(istep(k).EQ.0) THEN
			READ(1,*)(hramp(i,k),i=1,3)
		ELSE
			READ(1,*)hramp(1,k)
		ENDIF
	ENDDO			 

	IF(num_walls .GT. 0) READ(1,*)(nwpoint(i),i=1,num_walls)

	DO k=1,num_walls
		READ(1,*)(wpoint(k,i),i=1,nwpoint(k))
		READ(1,*)(wrecess(k,i),i=1,nwpoint(k))
	ENDDO

c	Read the different step height
c	=============================
	READ(1,*) xt,ht,rebase
	READ(1,*) crown,camber,twist
	READ(1,*)(xint(i),i=1,4)
	READ(1,*)(yint(i),i=1,4)

	CLOSE(1)
	RETURN

999	WRITE(*,*)'Trouble in opening files.'
	STOP
	END	

c======================================================================
	SUBROUTINE norm_rail
c======================================================================

c     Shared Data
c     ===========      
      USE syst_cnfg
      USE rail_data
      USE wall_data
      USE grid_ctrl
      USE sldr_dynm
      USE POI_cords
      
	IMPLICIT REAL*8(a-h,o-z)
	     
	DOUBLE PRECISION a(3,3),b(3)

	DO k=1,num_rails

		DO j=1,npoints(k)
			xrail(j,k) = xrail(j,k) / xl
			yrail(j,k) = yrail(j,k) / xl
		ENDDO
        
        CALL colinear(k)        
        
		xpt1(k)=1.d0
		xpt2(k)=0.d0
		ypt1(k)=1.d0
		ypt2(k)=0.d0

		DO j=1,npoints(k)
			xpt1(k)=DMIN1(xpt1(k),xrail(j,k))
			xpt2(k)=DMAX1(xpt2(k),xrail(j,k))
			ypt1(k)=DMIN1(ypt1(k),yrail(j,k))
			ypt2(k)=DMAX1(ypt2(k),yrail(j,k))
		ENDDO
	ENDDO

	DO k=1,num_walls
	   DO i=1, nwpoint(k)
		  wpoint(k,i)  = wpoint(k,i)  / xl
		  wrecess(k,i) = wrecess(k,i) / hm
	   ENDDO
	ENDDO

	DO i = 1, 4
	   xint(i) = xint(i) / xl
	   yint(i) = yint(i) / xl
	ENDDO

	ht = xt * ht
	ht = ht / hm
	xt = xt / xl
	tlnginc = tlnginc / xl

	crown  = crown  / hm
	crninc = crninc / hm
	camber = camber / hm
	cbrinc = cbrinc / hm
	twist  = twist  / hm
	twstinc= twstinc/ hm
	rebase = rebase / hm
	rcsinc = rcsinc / hm

	DO k=1,num_rails

		IF(istep(k).EQ.0) THEN
			a(1,1) =xrail(1,k)
			a(1,2) =yrail(1,k)
			a(1,3) =1.d0
			a(2,1) =xrail(2,k)
			a(2,2) =yrail(2,k)
			a(2,3) =1.d0
			a(3,1) =xrail(3,k)
			a(3,2) =yrail(3,k)
			a(3,3) =1.d0
			DO i=1,3
				hramp(i,k)=hramp(i,k)/hm
				b(i)=hramp(i,k)
			ENDDO
			CALL matrix33_inverse(a)
			CALL matrix_multi(3,3,1,a,b,cramp(1,k))
		ELSE
			hramp(1,k)=hramp(1,k)/hm
		ENDIF

c	Calculate the slope of each lines
c	================================
		DO i=1,npoints(k)
			ip1=i+1
			IF (i.EQ.npoints(k)) ip1=1
				xdiff(i,k)=xrail(ip1,k)-xrail(i,k)
				ydiff(i,k)=yrail(ip1,k)-yrail(i,k)
		ENDDO
	ENDDO	

	DO i = 2, nsx
	   xnt(i) = xnt(i) / xl
	ENDDO

	DO i = 2, nsy
	   ynt(i) = ynt(i) / xl
	ENDDO

c	Computes the wallprofile
c	========================
	CALL wall_profile

	RETURN

	END

c======================================================================
	SUBROUTINE wall_profile
c======================================================================

c     Shared Data
c     ===========
      USE syst_cnfg
      USE rail_data
      USE wall_data

	IMPLICIT REAL*8(a-h,o-z)
    
	DO k=1, num_rails
		xrail(npoints(k)+1,k)=xrail(1,k)
		yrail(npoints(k)+1,k)=yrail(1,k)
		xrail(npoints(k)+2,k)=xrail(2,k)
		yrail(npoints(k)+2,k)=yrail(2,k)

		IF(istep(k).NE.0.AND.hramp(1,k).LT.rebase)THEN
			ymax=-1.d0
			ymin=1.d0

c			Obtain the highest or lowest point
c			==================================
			DO i=1, npoints(k)
				IF(yrail(i,k).GT.ymax) THEN
					ymax=yrail(i,k)
					iymax=i
				ENDIF
			 
				IF(yrail(i,k).LT.ymin) THEN
					ymin=yrail(i,k)
					iymin=i
				ENDIF
			ENDDO

			IF(iymax.EQ.1)THEN
				im=iymin
			ELSE
				im=iymax
			ENDIF

			IF((xrail(im,k)-xrail(im-1,k))*(yrail(im+1,k)-
     &			yrail(im,k)).GT. (xrail(im+1,k)-xrail(im,k))*
     &			(yrail(im,k)-yrail(im-1,k))) THEN
					clock=-1.d0
			ELSE
				clock=1.d0
			ENDIF

			np=npoints(k)
			indexw(k,np+1)=indexw(k,1)
			
c			Since subscript starts from 1, we have to shift the loop
c			========================================================
			DO i=2,np+1
				im1=i-1
				ip1=i+1
		
				x1=xrail(im1,k)
				x2=xrail(i,k)
				x3=xrail(ip1,k)
		   
				y1=yrail(im1,k)
				y2=yrail(i,k)
				y3=yrail(ip1,k)
		   
				iwm=indexw(k,im1)
				iw=indexw(k,i)

				IF(iwm.EQ.0)THEN
					di1=0.d0
					do1=0.d0
				ELSE
					di1=wpoint(iwm,1)
					do1=wpoint(iwm,nwpoint(iwm))
				ENDIF

				IF(iw.EQ.0)THEN
					di2=0.d0
					do2=0.d0
				ELSE
					di2=wpoint(iw,1)
					do2=wpoint(iw,nwpoint(iw))
				ENDIF

c				Obtain wall profile vertices (xi,yi) and (xo, yo)
c				=================================================
				CALL intersection(clock,x1,y1,x2,y2,x3,y3,
     &				di1,do1,di2,do2,xwalli(i,k),ywalli(i,k),
     &				xwallo(i,k),ywallo(i,k))
			ENDDO
		
			xwalli(1,k)=xwalli(np+1,k)
			xwallo(1,k)=xwallo(np+1,k)
			ywalli(1,k)=ywalli(np+1,k)
			ywallo(1,k)=ywallo(np+1,k)
			
!......................................................................
!           new wall generation algorithm : 
!           We round the outer perimeter of the 
!           walls at convex corners
!......................................................................

            DO i=1, npoints(k)
                ip1=i+1
                iw=indexw(k, i)
                IF ((iw .LE. 0).OR.((iw-1).GT.num_walls)) THEN
                    flush_lower_xwallo(i,k) = xwalli(i  ,k)
	              flush_upper_xwallo(i,k) = xwalli(ip1,k)
	              flush_lower_ywallo(i,k) = ywalli(i  ,k)
	              flush_upper_ywallo(i,k) = ywalli(ip1,k)
                ELSE
                    x_mag = xwalli(ip1, k) - xwalli(i,k)
	              y_mag = ywalli(ip1, k) - ywalli(i,k)
!                   get length of our line
        	        vec_mag = SQRT(x_mag*x_mag + y_mag*y_mag)
!                   get nominal distance of wall (that's how far out we'll end up going
	              dnom_dist = (wpoint(iw, nwpoint(iw))) - 
     &	              (wpoint(iw, 1))
!                   Make sure we have a nominal distance and don't divide by 0
	              wall_percent = 0.0
	              IF (vec_mag .GT. 0.0) THEN 
	                  wall_percent = dnom_dist/vec_mag
	              ENDIF
!                   rotate 90 degrees
                    tempx =  x_mag
	              x_mag = -y_mag
	              y_mag = tempx
!                   shorten/lengthen our vector so it's the nominal wall length
    	              x_mag = x_mag * wall_percent * clock
	              y_mag = y_mag * wall_percent * clock
!                   translate xwalli line out to x_mag
	              flush_lower_xwallo(i,k) = x_mag + xwalli(i  ,k)
        	        flush_upper_xwallo(i,k) = x_mag + xwalli(ip1,k)
	              flush_lower_ywallo(i,k) = y_mag + ywalli(i  ,k)
        	        flush_upper_ywallo(i,k) = y_mag + ywalli(ip1,k)
                ENDIF
            ENDDO
!......................................................................       


		ENDIF
	ENDDO

	DO k=1,num_rails
		DO i=1,npoints(k)
			xw1(i,k)=DMIN1(xwalli(i,k),xwalli(i+1,k),
     &			xwallo(i,k),xwallo(i+1,k))
			xw2(i,k)=DMAX1(xwalli(i,k),xwalli(i+1,k),
     &			xwallo(i,k),xwallo(i+1,k))
		  
			yw1(i,k)=DMIN1(ywalli(i,k),ywalli(i+1,k),
     &			ywallo(i,k),ywallo(i+1,k))
			yw2(i,k)=DMAX1(ywalli(i,k),ywalli(i+1,k),
     &			ywallo(i,k),ywallo(i+1,k))
		ENDDO
	ENDDO

	RETURN
	END

c======================================================================
	SUBROUTINE pointrecess(xv,yv,recess)
c======================================================================

c----------------------------------------------------------------------
c	Recess from slider reference plane at(xv,yv)
c----------------------------------------------------------------------

c     Shared Data
c     ===========      
      USE syst_cnfg
      USE rail_data
      USE wall_data

      IMPLICIT REAL*8(a-h,o-z)
      
      LOGICAL inside,iPointInFlushBoundary,iPointInWallioBoundary
      DOUBLE PRECISION xw(6),yw(6)
      DOUBLE PRECISION xtemp(max_rail_pts),ytemp(max_rail_pts)
            
c     Assume base recess
c     ==================
      recess = rebase
      
c     Point (xv,yv) belong to rail k?
c     ==============================
      DO k=num_rails,1,-1
      
c       assume outside first
        inside=.FALSE.

c       Skip if outside range of this rail
        IF((yv.LT.ypt1(k)) .OR. (yv.GT.ypt2(k)) .OR. 
     &   (xv.LT.xpt1(k)) .OR. (xv.GT.xpt2(k)))  GOTO 10


c       Copy rail points to temporary arrays
c       ====================================
        DO it=1,npoints(k)
            xtemp(it)=xrail(it,k)
            ytemp(it)=yrail(it,k)
        ENDDO

c       Check whether the pt is inside the boundary
c       ===========================================
        CALL PNPOLY (xv,yv,xtemp,ytemp,npoints(k),i_temp) 
        IF (i_temp .EQ. -1) THEN
            inside=.FALSE.
            CALL PTONPOLYBOUNDARY(xv,yv,xtemp,ytemp,npoints(k),i_temp)
            IF (i_temp .EQ. 1) THEN
                inside=.TRUE.
            ENDIF
        ELSE
            inside=.TRUE.
        ENDIF


	    IF (inside) THEN
	        IF(istep(k).ne.0) THEN
	            recess=hramp(1,k)
	        ELSE
	            recess=xv*cramp(1,k)+yv*cramp(2,k)+cramp(3,k)
	        ENDIF
	        GOTO 20
	    ENDIF
10	    CONTINUE
	ENDDO !num_rails
      
20    CONTINUE

c     Now do walls
c     ============
      recessold = recess
      
c	Point (xv,yv) belong to sloped walls of rail k?	
c     ===============================================
	DO k = num_rails,1,-1
	
c	    Skip if not a raised step
c         =========================
	    IF((istep(k).EQ.0).OR.(hramp(1,k).GE.rebase)) GOTO 18
	    
	    DO ip=1,npoints(k)
	    
c	      Skip if outside range of this rail
c           ==================================	    
            IF((yv .LT. yw1(ip,k)) .OR. (yv .GT. yw2(ip,k)) .OR.
     &  	   (xv .LT. xw1(ip,k)) .OR. (xv .GT. xw2(ip,k)) .OR.
     &         (indexw(k,ip) .EQ. 0))  GOTO 15
     
c           See if we're inside the flush boundary
c           ======================================	    
		    xw(1) = xwalli(ip, k)
		    xw(2) = xwalli(ip+1, k)
		    xw(3) = flush_upper_xwallo(ip,k)
		    xw(4) = flush_lower_xwallo(ip,k)
		    xw(5) = xw(1)
		    xw(6) = xw(2)
		    yw(1) = ywalli(ip,k)
		    yw(2) = ywalli(ip+1,k)
		    yw(3) = flush_upper_ywallo(ip,k)
		    yw(4) = flush_lower_ywallo(ip,k)
		    yw(5) = yw(1)
		    yw(6) = yw(2)

		    iPointInFlushBoundary = .FALSE.
		    IF(indexw(k,ip) .NE. 0) THEN
		  	    iPointInFlushBoundary = IsPointInQuad(xv,yv,xw,yw)
		    ENDIF
		    
c           Now see if we're in the original (non rounded)
c           =============================================
	      xw(1)=xwalli(ip,k)
	      xw(2)=xwalli(ip+1,k)
	      xw(3)=xwallo(ip+1,k)
	      xw(4)=xwallo(ip,k)
	      xw(5)=xw(1)
	      xw(6)=xw(2)
	      yw(1)=ywalli(ip,k)
	      yw(2)=ywalli(ip+1,k)
	      yw(3)=ywallo(ip+1,k)
	      yw(4)=ywallo(ip,k)
	      yw(5)=yw(1)
	      yw(6)=yw(2)

	      iPointInWallioBoundary = .false.
	      IF(indexw(k,ip) .NE. 0) THEN
	          iPointInWallioBoundary = IsPointInQuad(xv,yv,xw,yw)
	      ENDIF
	      
c           Now see what type of corner we're at
c           ====================================
		    IF ((iPointInFlushBoundary  .EQ. .TRUE. ) .AND. 
     &	        (iPointInWallioBoundary .EQ. .FALSE.)) THEN
c		        Convex corner.  We're actually not on a wall
                GOTO 15 !next point on slider
            ELSEIF (iPointInFlushBoundary .EQ. .TRUE.) THEN
            
c               Regular wall boundary calculation
		        iedge=ip		        
	          clambda = ((yw(2)-yw(1))*(xv-xw(1))-(xw(2)-xw(1))*
     &  		          (yv-yw(1)))/((yw(2)-yw(1))*(xw(3)-xw(1))-
     &  	              (xw(2)-xw(1))*(yw(3)-yw(1)))
     
		        iw = indexw(k,ip)
	          ew = (clambda*wpoint(iw,nwpoint(iw)))+
     &	             ((1.d0-clambda)*wpoint(iw,1))
               
 	          DO ipt = 1, nwpoint(iw)-1
		            ew1=wpoint(iw,ipt)
		            ew2=wpoint(iw,ipt+1)
		            
		            IF((ew.GE.ew1) .AND. (ew.LT.ew2)) THEN
		                wr1 = wrecess(iw,ipt)
		                wr2 = wrecess(iw,ipt+1)
		                recess = wr1 + ((ew-ew1)/(ew2-ew1)*(wr2-wr1))
		                GOTO 25
		             ENDIF
	          ENDDO
25	          CONTINUE   

c           We can have negative wall distances.  
c           Check to see if our rail point is inside the rail poly or outside
c           Also include tolerance for floating point errors

                IF (ew .GT. 1.0d-12) recess=DMIN1(recess,recessold)
	          recessold=recess
	          
            ELSEIF (iPointInWallioBoundary .EQ. .TRUE.) THEN
c               New wall calculation for corner region	   
                recess = rebase
                xDistToCorner = xv - xwalli(ip, k)      !corner
                yDistToCorner = yv - ywalli(ip, k)
                fDistLower = SQRT((xDistToCorner*xDistToCorner) + 
     &                            (yDistToCorner*yDistToCorner))
                xDistToCorner = xv - xwalli(ip+1, k)    !corner
                yDistToCorner = yv - ywalli(ip+1, k)
                fDistUpper = SQRT((xDistToCorner*xDistToCorner) + 
     &                            (yDistToCorner*yDistToCorner))
                fDist = -1.0
                
                IF (fDistUpper < fDistLower) THEN 
                    fDist = fDistUpper 
                ELSE
                    fDist = fDistLower
                ENDIF
                
                	          
c               We'll use the existing trick to determine if we
c               should allow this point to be below the top rail           	      
	          clambda=((yw(2)-yw(1))* (xv-xw(1))- (xw(2)-xw(1))*
     & 		             (yv-yw(1)))/((yw(2)-yw(1))*(xw(3)-xw(1))-
     & 		             (xw(2)-xw(1))*(yw(3)-yw(1)))
     
                iw = indexw(k,ip)
	          ew = (clambda*wpoint(iw,nwpoint(iw))) + 
     &               ((1.d0-clambda) * wpoint(iw,1))
                distFromZero = wpoint(iw,1)
                     	          
  	          DO ipt=1,nwpoint(iw)-1
		            ew1 = wpoint(iw,ipt)   - distFromZero
		            ew2 = wpoint(iw,ipt+1) - distFromZero
		            IF(fDist.ge.ew1.and.fDist.lt.ew2) THEN
		                wr1 = wrecess(iw,ipt)
		                wr2 = wrecess(iw,ipt+1)
		                recess = wr1+((fDist-ew1)/(ew2-ew1)*(wr2-wr1))
		                GOTO 26
		            ENDIF
	          ENDDO
26              CONTINUE         
	          
	          IF(ew .GT. 1.0d-12) recess = DMIN1(recess,recessold)
	          recessold=recess         
            
	      ENDIF
	    
15	      CONTINUE 
          ENDDO	  !npoints	
18	    CONTINUE
      ENDDO !num_rails

      crwn=4.d0*crown*xv*(1.d0-xv)
      cmbr=4.d0*camber*yv*(yl-yv)/yl**2
      twst=4.d0*twist*(xv-0.5d0)*(yv-.5d0*yl)/yl
      recess=recess-crwn-cmbr-twst


	RETURN
      END

c======================================================================
      SUBROUTINE get_rail_size
c======================================================================

c     Shared Data
c     ===========
      USE rail_data
      
      IMPLICIT NONE

      INTEGER i,j,iJunk
      INTEGER no_rails, no_walls   

      INTEGER, DIMENSION (:), ALLOCATABLE :: rail_pts_list
      INTEGER, DIMENSION (:), ALLOCATABLE :: wall_pts_list
      
c     Initialize
c     ==========
      OPEN(1,ERR=999,FILE='rail.dat',STATUS='OLD')  
      READ(1, *)
      READ(1, *)
      READ(1, *)

c     Read in the number of rails and walls
c     =====================================
      READ(1,*) no_rails, no_walls
      
      ALLOCATE(rail_pts_list(no_rails))
      ALLOCATE(wall_pts_list(no_walls))
      
      rail_pts_list = 0
      wall_pts_list = 0
      
      DO i=1,no_rails
		READ(1,*) rail_pts_list(i), iJunk

c       Read in all the nodes for this rail and junk them
c       =================================================
	  DO j = 1,rail_pts_list(i)
	      READ(1,*)
	  ENDDO
 
c       Read in the recess/step info
c       ===========================
        READ(1,*)
      ENDDO

c     Read in the number of points in each wall
c     =========================================
      IF(no_walls .GT. 0) THEN
        READ(1,*)(wall_pts_list(i),i=1,no_walls)
      ENDIF	

      max_rail_pts = 1
      DO i = 1, no_rails
        IF (rail_pts_list(i) .GT. max_rail_pts) THEN 
            max_rail_pts = rail_pts_list(i)
        ENDIF
      ENDDO

      max_wall_pts = 1
      DO i = 1, no_walls
        IF (wall_pts_list(i) .GT. max_wall_pts) THEN 
            max_wall_pts = wall_pts_list(i)
        ENDIF
      ENDDO
    
      CLOSE(1)
      
      DEALLOCATE(rail_pts_list)
      DEALLOCATE(wall_pts_list)
      RETURN

999   WRITE(*,*)'Trouble opening rail.dat file.'
      STOP
      
      END
      
c======================================================================
      SUBROUTINE allocate_rail_arrays
c======================================================================

c     Shared Data
c     ===========
      USE rail_data
      USE wall_data

c     Allocate matrices
c     =================      
      ALLOCATE(istep(num_rails))
      ALLOCATE(npoints(num_rails))
      ALLOCATE(indexw(num_rails,max_rail_pts+1))
      
      ALLOCATE(xpt1(num_rails))
      ALLOCATE(xpt2(num_rails))
      ALLOCATE(ypt1(num_rails))
      ALLOCATE(ypt2(num_rails))   
      
      ALLOCATE(cramp(3,num_rails))
      ALLOCATE(hramp(3,num_rails))
      
      ALLOCATE(xrail(max_rail_pts+2,num_rails))
      ALLOCATE(yrail(max_rail_pts+2,num_rails))
      ALLOCATE(xdiff(max_rail_pts+2,num_rails))
      ALLOCATE(ydiff(max_rail_pts+2,num_rails))
      
      ALLOCATE(nwpoint(num_walls))
      
      ALLOCATE(wpoint (num_walls,max_rail_pts))
      ALLOCATE(wrecess(num_walls,max_rail_pts))
      
      ALLOCATE(xwalli(max_rail_pts+1,num_rails))
      ALLOCATE(ywalli(max_rail_pts+1,num_rails))
      ALLOCATE(xwallo(max_rail_pts+1,num_rails))
      ALLOCATE(ywallo(max_rail_pts+1,num_rails))
      
      ALLOCATE(flush_lower_xwallo(max_rail_pts,num_rails))
      ALLOCATE(flush_lower_ywallo(max_rail_pts,num_rails))
      ALLOCATE(flush_upper_xwallo(max_rail_pts,num_rails))
      ALLOCATE(flush_upper_ywallo(max_rail_pts,num_rails))
      
      ALLOCATE(xw1(max_rail_pts,num_rails))
      ALLOCATE(yw1(max_rail_pts,num_rails))
      ALLOCATE(xw2(max_rail_pts,num_rails))
      ALLOCATE(yw2(max_rail_pts,num_rails))
 
c     Initialize Arrays
c     ================
      istep = 0
      indexw = 0
      npoints = 0
      
      xpt1 = 0.0d0
      xpt2 = 0.0d0
      ypt1 = 0.0d0
      ypt2 = 0.0d0
      
      cramp = 0.0d0
      hramp = 0.0d0
      
      xrail = 0.0d0
      yrail = 0.0d0
      xdiff = 0.0d0
      ydiff = 0.0d0
      
      nwpoint = 0
      
      wpoint  = 0.0d0
      wrecess = 0.0d0
      
      xwalli = 0.0d0
      ywalli = 0.0d0
      xwallo = 0.0d0
      ywallo = 0.0d0
      
      flush_lower_xwallo = 0.0d0
      flush_lower_ywallo = 0.0d0
      flush_upper_xwallo = 0.0d0
      flush_upper_ywallo = 0.0d0
      
      xw1 = 0.0d0
      yw1 = 0.0d0
      xw2 = 0.0d0
      yw2 = 0.0d0
      
      RETURN
      END
      
      
c======================================================================
      SUBROUTINE colinear(k)
c======================================================================

c     Shared Data
c     ===========
      USE rail_data
            
      IMPLICIT NONE
 
c     Local variables
c     ===============
      INTEGER n, j, k, nrind1, nrind2, nrind3
      DOUBLE PRECISION a1, a2, b1, b2, determ0


c     Only remove if at least 3 points in rail
c     ========================================
      IF (npoints(k) .GT. 2) THEN
     
        j = 3
        DO WHILE (j.LE.npoints(k))
      
            nrind1 = j
            nrind2 = j+1
            nrind3 = j+2
            
            IF (nrind2 .GT. npoints(k)) THEN
                nrind2 = 1
                nrind3 = 2
            ENDIF
            
            IF(nrind3 .GT. npoints(k)) THEN
                nrind3 = 1
            ENDIF
            
c           Colinearity test
c           =================
            a1 = (yrail(nrind2,k) - yrail(nrind1,k))
            a2 = (yrail(nrind3,k) - yrail(nrind2,k))
            b1 = (xrail(nrind1,k) - xrail(nrind2,k))
            b2 = (xrail(nrind2,k) - xrail(nrind3,k))
            determ0 = a1*b2-a2*b1
            
            IF(dabs(determ0) .LT .1.e-08) THEN
                WRITE(*, '(A27, I3)') ' Removing colinear point = 
     &              in rail ',k
                j = 3
                DO n = nrind2,npoints(k)-1
                    xrail(n,k) = xrail(n+1,k)
                    yrail(n,k) = yrail(n+1,k)
                    indexw(k,n)= indexw(k,n+1)
                ENDDO
                npoints(k) = npoints(k)-1
            ENDIF
            j = j+1
        ENDDO
10    ENDIF

	RETURN

	END	
    
c======================================================================
      SUBROUTINE intersection(clock,x1,y1,x2,y2,x3,y3,
     &						di1,do1,di2,do2,xi,yi,xo,yo)
c======================================================================

	IMPLICIT REAL*8(a-h,o-z)

	a1=y2-y1
	a2=y3-y2
	b1=x1-x2
	b2=x2-x3
	c1=x1*y2-x2*y1
	c2=x2*y3-x3*y2
	coe1=dsqrt(a1*a1+b1*b1)
	coe2=dsqrt(a2*a2+b2*b2)

	c1i=c1-coe1*di1*clock
	c2i=c2-coe2*di2*clock
	c1o=c1-coe1*do1*clock
	c2o=c2-coe2*do2*clock

	determ0=a1*b2-a2*b1

	IF(DABS(determ0) .LT. 1.e-30) THEN
	   WRITE(*,*)'Devide by zero in Subroutine intersection!'
	   WRITE(*,*)'Colinear points are not allowed in the rail!'
	ENDIF

	determxi=c1i*b2-c2i*b1
	determyi=a1*c2i-a2*c1i
	determxo=c1o*b2-c2o*b1
	determyo=a1*c2o-a2*c1o

	xi=determxi/determ0
	yi=determyi/determ0
	xo=determxo/determ0
	yo=determyo/determ0

	RETURN

	END	

c======================================================================
      FUNCTION IsPointInQuad(xv,yv,xw,yw)
c======================================================================
      IMPLICIT NONE

      DOUBLE PRECISION xw(6),yw(6),xv,yv,xcross    
      INTEGER i, IsPointInQuad, ncross
      LOGICAL inside
    
      inside = .FALSE.
      ncross = 0
      DO i = 1,4
	    IF(DABS(yw(i+1)-yw(i)).GT.1.d-10)THEN
	        IF(((yw(i).GT.yv).AND.(yw(i+1).LT.yv))
     &     	    .OR.((yw(i).LT.yv).AND.(yw(i+1).GT.yv)))THEN
		        xcross=xw(i)+(xw(i+1)-xw(i))/(yw(i+1)-yw(i)) *(yv-yw(i))
		        IF(xv.GT.xcross) ncross=ncross+1
		        IF(ncross.GT.1)GOTO 53
		    ELSE IF(yw(i+1).EQ.yv)THEN
		        IF(((yw(i+1)-yw(i+2))*(yw(i+1)-yw(i)).LT.0.d0) 
     &		      .AND. (xv.GT.xw(i+1))) ncross=ncross+1
		        IF(ncross.GT.1) GOTO  53
		    ENDIF
	    ELSE
		    IF(((xv.GE.xw(i).AND.xv.LE.xw(i+1)) 
     &	         .OR.(xv.LE.xw(i).AND.xv.GE.xw(i+1))) 
     &	         .AND.DABS(yv-yw(i)).LT.1.d-10)THEN
		        inside=.TRUE.
		        GOTO 52
		    ENDIF
	    ENDIF
	ENDDO
	
52    CONTINUE
	IF(ncross.EQ.1) inside=.TRUE.
53    CONTINUE
	
	IsPointInQuad = inside
	END
	
c======================================================================
      SUBROUTINE PNPOLY(PX, PY, XX, YY, N, INOUT) 
c======================================================================

!  From BCox 2/17/05: This subroutine doesn't check to see if we're
!  on the boundary of ALL points.  This is bad since we snap the grid
!  to our boundary!  From Prof. Franklin himself:
!  "If you want to know when a point is exactly on the boundary, 
!   you need another program"
!   See PTONPOLYBOUNDARY below

!  Courtesy of Prof. Randolph Franklin, Rensselaer Polytechnic Inst, Troy NY
!     ..................................................................
!                                                                       
!        SUBROUTINE PNPOLY                                              
!                                                                       
!        PURPOSE                                                        
!           TO DETERMINE WHETHER A POINT IS INSIDE A POLYGON            
!                                                                       
!        USAGE                                                          
!           CALL PNPOLY (PX, PY, XX, YY, N, INOUT )                     
!                                                                       
!        DESCRIPTION OF THE PARAMETERS                                  
!           PX      - X-COORDINATE OF POINT IN QUESTION.                
!           PY      - Y-COORDINATE OF POINT IN QUESTION.                
!           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
!                     VERTICES OF POLYGON.                              
!           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
!                     VERTICES OF POLYGON.                              
!           N       - NUMBER OF VERTICES IN THE POLYGON.                
!           INOUT   - THE SIGNAL RETURNED:                              
!                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
!                      0 IF THE POINT IS ON AN EDGE OR AT A VERTEX,     
!                      1 IF THE POINT IS INSIDE OF THE POLYGON.         
!                                                                       
!        REMARKS                                                        
!           THE VERTICES MAY BE LISTED CLOCKWISE OR ANTICLOCKWISE.      
!           THE FIRST MAY OPTIONALLY BE REPEATED, IF SO N MAY           
!           OPTIONALLY BE INCREASED BY 1.                               
!           THE INPUT POLYGON MAY BE A COMPOUND POLYGON CONSISTING      
!           OF SEVERAL SEPARATE SUBPOLYGONS. IF SO, THE FIRST VERTEX    
!           OF EACH SUBPOLYGON MUST BE REPEATED, AND WHEN CALCULATING   
!           N, THESE FIRST VERTICES MUST BE COUNTED TWICE.              
!           INOUT IS THE ONLY PARAMETER WHOSE VALUE IS CHANGED.         
!           THE SIZE OF THE ARRAYS MUST BE INCREASED IF N > MAXDIM      
!           WRITTEN BY RANDOLPH FRANKLIN, UNIVERSITY OF OTTAWA, 7/70.   
!                                                                       
!        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  
!           NONE                                                        
!                                                                       
!        METHOD                                                         
!           A VERTICAL LINE IS DRAWN THRU THE POINT IN QUESTION. IF IT  
!           CROSSES THE POLYGON AN ODD NUMBER OF TIMES, THEN THE        
!           POINT IS INSIDE OF THE POLYGON.                             
!     ..................................................................
                                                                      
      IMPLICIT NONE
	INTEGER maxRailPts
	PARAMETER(maxRailPts=200)
		
	INTEGER i, J, INOUT,N,O
      DOUBLE PRECISION PX, PY, X(maxRailPts),Y(maxRailPts),XX(N),YY(N)
      LOGICAL MX,MY,KX,KY                                              
	
c     OUTPUT UNIT FOR PRINTED MESSAGES                                 
      DATA O/6/                                                         
    
c     MAXDIM=300
c     IF(N.LE.MAXDIM)GO TO 6                                            
c     WRITE(*,*) '0WARNING: TOO MANY POINTS USED TO DEFINE A RAIL : ', N
c     RETURN
                                                            
6     DO 1 I=1,N                                                        
        X(I)=XX(I)-PX
1       Y(I)=YY(I)-PY
        INOUT=-1 
        DO 2 I=1,N
        J=1+MOD(I,N)
        MX=X(I).GE.0.0
        KX=X(J).GE.0.0
        MY=Y(I).GE.0.0
        KY=Y(J).GE.0.0
        IF(.NOT.((MY.OR.KY).AND.(MX.OR.KX)).OR.(MX.AND.KX)) GOTO 2
        IF(.NOT.(MY.AND.KY.AND.(MX.OR.KX).AND..NOT.(MX.AND.KX))) GOTO 3
        INOUT=-INOUT
        GOTO 2  
3       IF((Y(I)*X(J)-X(I)*Y(J))/(X(J)-X(I))) 2,4,5 
4       INOUT=0
        RETURN
5       INOUT=-INOUT
2     CONTINUE                                                          
      RETURN                                                            
      END
      
c======================================================================
      SUBROUTINE PTONPOLYBOUNDARY(PX, PY, XX, YY, N, INOUT) 
c======================================================================

!     ..................................................................
!        DESCRIPTION OF THE PARAMETERS                                  
!           PX      - X-COORDINATE OF POINT IN QUESTION.                
!           PY      - Y-COORDINATE OF POINT IN QUESTION.                
!           XX      - N LONG VECTOR CONTAINING X-COORDINATES OF         
!                     VERTICES OF POLYGON.                              
!           YY      - N LONG VECTOR CONTAING Y-COORDINATES OF           
!                     VERTICES OF POLYGON.                              
!           N       - NUMBER OF VERTICES IN THE POLYGON.                
!           INOUT   - THE SIGNAL RETURNED:                              
!                     -1 IF THE POINT IS OUTSIDE OF THE POLYGON,        
!                      1 IF THE POINT IS ON AN EDGE OR AT A VERTEX
!     ..................................................................

      IMPLICIT NONE
    
      INTEGER INOUT, I, J, N, MAXDIM
      DOUBLE PRECISION PX, PY  
      DOUBLE PRECISION XX(N),YY(N),YMIN,YMAX,XMIN,XMAX
    
      INOUT=-1
      DO I=1,N
        J=1+MOD(I,N)

c       First check if we're on a corner
        IF ((PX .EQ. XX(i)) .AND. (PY .EQ. YY(i))) THEN
            INOUT = 1
            RETURN
        ENDIF

c       Now check if we're sitting on a line 
c       in the X or Y direction only
        IF ((XX(i) - XX(j) .EQ. 0.d0) .AND. (PX .EQ. XX(i))) THEN
        
c       We're inline with a vertical line so see if 
c       we're between the endpoints
            IF (YY(i) .LT. YY(j)) THEN 
                YMIN = YY(i)
                YMAX = YY(j)
            ELSE
                YMIN = YY(j)
                YMAX = YY(i)
            ENDIF
            IF ((PY .GT. YMIN) .AND. (PY .LT. YMAX)) THEN
                INOUT = 1
                RETURN
            ENDIF
        ENDIF
        
        IF ((YY(i) - YY(j) .EQ. 0.d0) .AND. (PY .EQ. YY(i)))THEN
        
c       We're inline with a horizontal line, 
c       see if we're between the endpoints        
            IF (XX(i) .LT. XX(j)) THEN 
                XMIN = XX(i)
                XMAX = XX(j)
            ELSE
                XMIN = XX(j)
                XMAX = XX(i)
            ENDIF
            IF ((PX .GT. XMIN) .AND. (PX .LT. XMAX)) THEN
                INOUT = 1
                RETURN
            ENDIF
        ENDIF
      ENDDO
      
c     Not on a vertical or horizontal line so (for now) we'll assume that 
c     we're not on the boundary (DANGEROUS, PLEASE FIX!!!)

      RETURN                                  
      END
           
c======================================================================    
      