c**********************************************************************
c     Subroutines in this file :
c     1. dynamic_solver
c     2. dynamic_solver_contact
c     3. compute_shock
c     4. adv_actuator
c**********************************************************************

c======================================================================
	SUBROUTINE dynamic_solver
c======================================================================

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
      USE vdwf_data
      USE elec_stat
      USE POI_cords
      USE luls_data
      USE shck_data
      USE cmpt_dfmn
      USE susp_data
      USE disk_data
      USE disk_aray
      USE motr_data
      USE rotr_data
      USE hsng_data
      USE base_data
      USE supp_data


	IMPLICIT REAL*8(a-h,o-z)
	
	DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979
	
	LOGICAL crash
	INTEGER no_stps
	DOUBLE PRECISION disp_cur(3),disp_prv(3)
	DOUBLE PRECISION ang_vel,ang_cur,s_cur,c_cur
      DOUBLE PRECISION xtemp,ytemp,rad_temp,ang_temp 
      DOUBLE PRECISION w_trn_0,Dw_rad_0,Dw_ang_0

c	Local Variables
c	===============
	DOUBLE PRECISION flht(4),hloc(4)

c	POI base recesses from slider surface
c	=====================================
	DO i=1,4
	  CALL pointrecess(xint(i),yint(i),thh(i))
	ENDDO
	
	WRITE(*,*)
	WRITE(*,*) 'START OF TRANSIENT ANALYSIS'
	WRITE(*,*) '==========================='
	
	t=0.d0
	no_stps = INT(f_dtfac)

c	Suspension kinematics variable
c	==============================	
	dst=dact																  
	vst=vact
	aact=0.d0

c	Force (contact)
c	===============
	fcont=0.d0	! Impact force
	fsctx=0.d0	! Roll Torque 
	fscty=0.d0	! Pitch Torque
	aratio1=0.d0
	f0=f0*9.81d0

	IF(rpm.EQ.0.d0) THEN
		omega=0.d0
	ENDIF

	dfx=0.d0
	dfy=ra/xl
	IF(iadapt.NE.1) THEN
		CALL gethx
		CALL get_toth(crash,iic,jjc)
		IF(crash) THEN
			STOP
		ENDIF
	
		IF(omega.NE.0.d0) CALL reyneq(1)
	ENDIF

	CALL calc_ABS_force(0)
	CALL calc_electrostatic
	CALL calc_contact_force(0)
	

c	Initialize disk displacement at t = 0
c	=====================================
	IF ((disk_opt.EQ.1) .AND. (luls_opt.EQ.0)) THEN

c       Time Variables
c       ==============
		t_val = t
		
c       Disk displacement at slider
c       ===========================		
		w_trn		= 0.0d0
		Dw_rad		= 0.0d0
		Dw_ang		= 0.0d0
		Dw_x_susp	= 0.0d0
		Dw_y_susp	= 0.0d0
		
c       Disk center displacement
c       ========================
        x0_disk_r = 0.0d0
        y0_disk_r = 0.0d0
        x0_disk_s = 0.0d0
        y0_disk_s = 0.0d0
        disk_cent_disp = 0.0d0

c       Disk attitude for interpolation
c       ===============================		
		disp_cur = 0.0d0
		disp_prv = 0.0d0
		
c       Disk displacement vector
c       ========================
	  w_disk_cur = 0.0d0
	  
c       Vectors for system dynamics
c       =========================== 
        IF (supp_opt .EQ. 0) THEN
	      v_disk_cur = 0.0d0
	      a_disk_cur = 0.0d0
	      w_disk_prv = 0.0d0
	      v_disk_prv = 0.0d0
	      a_disk_prv = 0.0d0
	      f_disk_cnst= 0.0d0
        ELSEIF (supp_opt .EQ. 1) THEN
	      w_rotr_cur = 0.0d0
	      v_rotr_cur = 0.0d0
	      a_rotr_cur = 0.0d0
	      w_rotr_prv = 0.0d0
	      v_rotr_prv = 0.0d0
	      a_rotr_prv = 0.0d0
	      f_rotr_cnst= 0.0d0
	  ELSEIF (supp_opt .EQ. 2) THEN
	      w_motr_cur = 0.0d0
	      v_motr_cur = 0.0d0
	      a_motr_cur = 0.0d0
	      w_motr_prv = 0.0d0
	      v_motr_prv = 0.0d0
	      a_motr_prv = 0.0d0
	      f_motr_cnst= 0.0d0
	  ELSEIF (supp_opt .EQ. 3) THEN
		    w_supp_cur = 0.0d0
	      v_supp_cur = 0.0d0
	      a_supp_cur = 0.0d0
	      w_supp_prv = 0.0d0
	      v_supp_prv = 0.0d0
	      a_supp_prv = 0.0d0
	      f_supp_cnst= 0.0d0
	  ENDIF	
	  
c       Initialize full disk displacement
c       ==================================
        IF (data_opt_disk .EQ. 1) THEN
            ALLOCATE(w_disk(no_disk_dofs+1))
            w_disk = 0.0d0
        ENDIF
        
c       Initialize disk/motor data vector
c       =================================        
   	  IF ((supp_opt.NE.0) .AND. (data_opt_supp.EQ.1)) THEN

c           Disk - Hub interface displacement
c           =================================
            ALLOCATE(w_disk_intf((3*no_disk_intf_dofs)+1))
            w_disk_intf = 0.0d0
            
c           FDB displacement
c           ================
            IF (supp_opt .EQ. 1) THEN
                no_fdb_dofs_op = no_fdb_dofs_hub
            ELSEIF (supp_opt .EQ. 2) THEN
                no_fdb_dofs_op = no_fdb_dofs_motr
            ELSEIF (supp_opt .EQ. 3) THEN
                no_fdb_dofs_op = no_fdb_dofs_motr    
            ENDIF
            
            ALLOCATE(w_FDB(no_fdb_dofs_op+1))
            w_FDB = 0.0d0
            
		ENDIF
		
	ENDIF

c	Suspension matrices LU decomposition
c	====================================
	CALL susp_block_fact(1)
	
c     Shock vector for actuator frame
c     ===============================

c     Angle substended at VCM center	
	sldr_vcm_ang = DACOS((len_c2c**2 + xact**2 - ra**2)/
     &                     (2*len_c2c*xact))
c     Angle with +x_HDD axis
      sldr_hdd_ang = base_ang - sldr_vcm_ang 

c	Initial suspension displacment
c	==============================
      CALL solve_susp(1)

c	Add disk attitude
c	=================
	zc  =  zc_susp - (w_trn/hm)
	hx0 = hx0_susp + (Dw_x_susp*xl/hm)
	hy  =  hy_susp + (Dw_y_susp*xl/hm)

c     Slider attitude
c     ===============
   	hp   = -hx0*hm/xl
	hyy  =  hy*hm/xl
	hmin =  zc+hx0*(1.d0-xg)

c	Find the point of minimum flying height
c	=======================================
	uuww=1d30
	DO i=1,nx
		DO j=1,ny
			IF (hnew(i,j).LT.uuww) THEN
				uuww=hnew(i,j)
				iuu=i
				juu=j
			ENDIF
		ENDDO
	ENDDO  

c	Output Pressure
c	===============
	IF (outp_opt.EQ.1) THEN
		DO i=1,4
			IF ((t.GE.toutp(i)).AND.(INT(foutp(i)).NE.1)) THEN
				CALL wr_pres(i)     ! Write Pressure data
				foutp(i) = 1.0d0    ! Turn on the flag 
			ENDIF
		ENDDO	
	ENDIF

c     Displacement at POI (4 points) from disk mean surface
c     =====================================================	  
      ssk=dsin(ske)
	csk=dcos(ske)

	DO i=1,4
		xtemp=xg-xint(i)            
		ytemp=yint(i)-(0.5d0*yl+yg)
		xloc=dfx+(xtemp*csk-ytemp*ssk)
		yloc=dfy+(xtemp*ssk+ytemp*csk)
		CALL get_floor(xloc,yloc,hloc(i))
		flht(i)=hmin*hm+hp*(1.d0-xint(i))*xl-
     &			(hy*hm/xl)*(0.5d0*yl-yint(i))*xl+
     &			hm*thh(i)
	ENDDO

c	Calculate min FH and min FH loc
c	===============================
      flhtmin=100.0d0
      xflhtmin=0.0d0
      yflhtmin=0.0d0
      spacemin=100.0d0
      xspacemin=0.0d0
      yspacemin=0.0d0
      DO i=1,nx
		DO j=1,ny
			xloc=xref(i)
			yloc=yref(j)
			CALL pointrecess(xloc,yloc,tempthh)
			tflht=hmin*hm+hp*(1.d0-xloc)*xl-
     &            (hy*hm/xl)*(0.5d0*yl-yloc)*xl+
     &             hm*tempthh
			IF (tflht.LT.flhtmin) THEN
				flhtmin=tflht
				xflhtmin=xloc
				yflhtmin=yloc
			ENDIF

c	      Now calculate minimum spacing
c	      =============================
			xtemp=xg-xref(i)
			ytemp=yref(j)-(0.5d0*yl+yg)
			xfloc=dfx+(xtemp*csk-ytemp*ssk)
			yfloc=dfy+(xtemp*ssk+ytemp*csk)
			CALL get_floor(xfloc,yfloc,temphloc)
			tspace=tflht+(hm*temphloc)
			IF(tspace.LT.spacemin) THEN
				spacemin=tspace
				xspacemin=xloc
				yspacemin=yloc
			ENDIF
			spa(i,j)=tspace*1.d0+0.12d-9
          ENDDO
      ENDDO

c     POI data
c     ========                               
	WRITE(41,101)t,flht(1),flht(2),flht(3),flht(4)                    
	WRITE(42,101)t,-hm*hloc(1),-hm*hloc(2),-hm*hloc(3),              
     &				 -hm*hloc(4)

c     Flying dynamics
c     ===============
      WRITE(51,101)t,hmin*hm,hp,-hyy,zc*hm,
     &	uuww*hm,xref(iuu)*xl,yref(juu)*xl
      WRITE(52,101)t,f,xf,yf,fp*9.81d0,ABS_data(1),ABS_data(2),
     &  fn*9.81d0,ABS_data(3),ABS_data(4)
      WRITE(53,101)t,f_maxcrp,f_maximp
      WRITE(54,101)t,fcr,txr,tyr,aratio
	WRITE(55,101)t,fcont,fsctx,fscty,aratio1
	WRITE(56,101)t,fvdw,yfvdw,xfvdw
	WRITE(57,101)t,felecst,telecst,relecst
	
c     Disk Angular velocity
c     =====================
      ang_vel = disk_RPM*(2.0d0*pi/60.0d0)
      
c     Update actuator position
c     ========================	
	CALL adv_actuator
	IF (luls_opt.EQ.0) THEN
	  WRITE(71,101)t,-acc_val,w_trn,Dw_x_susp,Dw_y_susp,ramp_cntf
	ENDIF

c     Update Parameters
c     =================
      h(1:nx,1:ny) = hnew(1:nx,1:ny)
      pold(1:nx,1:ny) = p(1:nx,1:ny)
	
5	CONTINUE	! Time loop

c     Initialize normalized error
c     ===========================
	zdf1 = 10000.d0
	zdf2 = 10000.d0
	zdf3 = 10000.d0

c	Initialize forces, moments and pressre
c	======================================
	fcont	   = 0.d0
	f_maximp = 0.d0
	fsctx	   = 0.d0
	fscty	   = 0.d0
	
c     Previous time step dynamic variables
c     ====================================
      IF (supp_opt .EQ. 0) THEN
	  w_disk_prv = w_disk_cur
		v_disk_prv = v_disk_cur
		a_disk_prv = a_disk_cur
	ELSEIF (supp_opt .EQ. 1) THEN
		w_rotr_prv = w_rotr_cur
		v_rotr_prv = v_rotr_cur
		a_rotr_prv = a_rotr_cur
	ELSEIF (supp_opt .EQ. 2) THEN
		w_motr_prv = w_motr_cur
		v_motr_prv = v_motr_cur
		a_motr_prv = a_motr_cur
      ELSEIF (supp_opt .EQ. 3) THEN
		w_supp_prv = w_supp_cur
		v_supp_prv = v_supp_cur
		a_supp_prv = a_supp_cur    
      ENDIF
      
      disp_prv(1) = w_trn
      disp_prv(2) = Dw_rad
      disp_prv(3) = Dw_ang
      
c     Shock magnitude (inertia loading : a)
c	====================================
      t_val = t+dt
      CALL compute_shock(t_val)

c     Disk Rotation
	ang_cur = ang_vel*t_val
	s_cur = DSIN(ang_cur)
	c_cur = DCOS(ang_cur)

c	Disk displacement at t = t_coarse
c	=================================
	IF ((luls_opt.EQ.0) .AND. (disk_opt.EQ.1)) THEN

c		Coarse Time Step
c		================
c---------------------------------------------
c		Note: Solving to t+dt step and then 
c		recording t+dt in log file 
c---------------------------------------------
		t_val = t+dt 

c		Disk Response
c		=============
        IF (supp_opt .EQ. 0) THEN      ! Fixed Disk
            CALL solve_disk
        ELSEIF (supp_opt .EQ. 1) THEN  ! Fixed Housing
            CALL solve_rotor
        ELSEIF (supp_opt .EQ. 2) THEN  ! Full Spindle Motor
            CALL solve_motor
        ELSEIF (supp_opt .EQ. 3) THEN  ! Full HDD Modal
            CALL solve_support
        ENDIF
        
c       Slider positioning
c       ==================        
        sldr_rad = ra   ! Radius
        
c       Angle substended at disk center	
	  sldr_disk_ang = DACOS((len_c2c**2 + ra**2 - xact**2 )/
     &                  (2*len_c2c*ra))
        sldr_ang  =  sldr_disk_ang + base_ang - pi

c       Disk center displacement
c       ========================
        disk_cent_disp = 0.0d0
        DO i = 1, no_disk_intf_dofs
            DO j = 1, 2
                disk_cent_disp(j) = disk_cent_disp(j) + 
     &                              disk_ID_disp(i,j)
            ENDDO
        ENDDO
        
        disk_cent_disp(1) = disk_cent_disp(1)/DBLE(no_disk_intf_dofs)
        disk_cent_disp(2) = disk_cent_disp(2)/DBLE(no_disk_intf_dofs)
        
c       Disk Center Displacement (Rotary frame)		
c       =======================================
		x0_disk_r = disk_cent_disp(1)
		y0_disk_r = disk_cent_disp(2)
		
c       Disk Center Displacement (Stationary frame)
c       ===========================================		
		x0_disk_s = (x0_disk_r*c_cur) - (y0_disk_r*s_cur)
		y0_disk_s = (x0_disk_r*s_cur) + (y0_disk_r*c_cur)
		
c       Adjust for disk center movement
c       ===============================
        xtemp = sldr_rad*DCOS(sldr_ang) - x0_disk_s
        ytemp = sldr_rad*DSIN(sldr_ang) - y0_disk_s		        
	
        rad_temp = DSQRT((xtemp**2) + (ytemp**2))
        ang_temp = ATAN(ABS(ytemp/xtemp))
            
        IF ((xtemp .LT. 0) .AND. (ytemp .LT. 0)) THEN
            ang_temp = pi + ang_temp
        ELSEIF ((xtemp .GE. 0) .AND. (ytemp .LT. 0)) THEN
            ang_temp = 2*pi - ang_temp
        ELSEIF ((xtemp .LT. 0) .AND. (ytemp .GE. 0)) THEN
            ang_temp =   pi - ang_temp
        ENDIF

c		Disk apttitude at slider postion
c		================================
		CALL disk_attitude(rad_temp,ang_temp,w_trn_0,Dw_rad_0,Dw_ang_0)

c       Actuator orientation
c       ====================
		ssk = DSIN(ske)
		csk = DCOS(ske)
		
c       Disk displacement for slider
c       ============================
        w_trn  = -w_trn_0
        Dw_rad = -Dw_rad_0
        Dw_ang = Dw_ang_0

c       Disk displacement along suspension
c       ==================================
		Dw_x_susp = -(Dw_rad*ssk) - (Dw_ang*csk)
		Dw_y_susp =  (Dw_rad*csk) - (Dw_ang*ssk)

	ELSE
		w_trn  = 0.d0
		Dw_rad = 0.d0
		Dw_ang = 0.d0
		Dw_x_susp = 0.d0
		Dw_y_susp = 0.d0

	ENDIF
	
c     Current time step disk attitude
c     ===============================		
      disp_cur(1) = w_trn
      disp_cur(2) = Dw_rad
      disp_cur(3) = Dw_ang
	
c	Suspension matrices LU decomposition
c	====================================
	CALL susp_block_fact(1)
	
c     Shock vector for actuator frame
c     ===============================

c     Angle substended at VCM center	
	sldr_vcm_ang = DACOS((len_c2c**2 + xact**2 - ra**2)/
     &                     (2*len_c2c*xact))
c     Angle with +x_HDD axis
      sldr_hdd_ang = base_ang - sldr_vcm_ang

	kint=0        ! Initialize iteration counter

10	CONTINUE	    ! Iteration Loop
	
	kint=kint+1   ! Increase iteration number

c     Check for crash in starting
c     ===========================
	CALL gethx
	CALL get_toth(crash,iic,jjc)
	IF (crash) THEN
	
c       Interpolate disk amplitude for fine time steps
c       ==============================================
	  IF ((disk_opt.EQ.1) .AND. (luls_opt.EQ.0)) THEN
            CALL intrp_disk_atd(disp_cur,disp_prv,no_stps,disk_atd)
	  ELSE
	      disk_atd = 0.0d0
	  ENDIF	
	
		CALL dynamic_solver_contact(0)

		IF(t.LE.tf) GOTO 5
		WRITE(*,*) 'END OF TRANSIENT ANALYSIS'
		RETURN
	ENDIF

c	Previous iteration value of slider attitude 
c	===========================================
	zcpre   = zc
	hx0pre  = hx0
	hypre   = hy

c	Previous iteration value of normalized error
c	============================================
	zdf1pre = zdf1
	zdf2pre = zdf2
	zdf3pre = zdf3

c	Solving reynold equation
c	========================
	IF(omega.NE.0.d0) CALL reyneq(0)
	CALL calc_ABS_force(0)
	CALL calc_electrostatic
	CALL calc_contact_force(0)

c	Suspension displacment
c	======================
      CALL solve_susp(1)

c	Add disk attitude
c	=================
	zc  =  zc_susp - (w_trn/hm)
	hx0 = hx0_susp + (Dw_x_susp*xl/hm)
	hy  =  hy_susp + (Dw_y_susp*xl/hm)

c     Slider attitude
c     ===============
   	hp   = -hx0*hm/xl
	hyy  =  hy*hm/xl
	hmin =  zc+hx0*(1.d0-xg)

c	Normalized error
c	================
	zdf1 = DABS((zcpre-zc)/zc)
	zdf2 = DABS((hx0pre-hx0)/hx0)
	zdf3 = DABS((hypre-hy)/hy)

c	If not converging (no_itrn > 10) use small time step
c	====================================================
	IF(kint.GT.10) THEN
		IF(((zdf1.GT.zdf1pre.AND.zdf1.GT.emax).OR.(zdf2.GT.zdf2pre.
     &		AND.zdf2.GT.emax).OR.(zdf3.GT.zdf3pre.AND.zdf3.GT.emax)).
     &		AND.(t.GT.2d-3)) THEN
			WRITE(*,445)t
			
c           Interpolate disk amplitude for fine time steps
c           ==============================================
	      IF ((disk_opt.EQ.1) .AND. (luls_opt.EQ.0)) THEN
                CALL intrp_disk_atd(disp_cur,disp_prv,no_stps,disk_atd)
	      ELSE
	          disk_atd = 0.0d0
	      ENDIF
			
			CALL dynamic_solver_contact(0)
			IF(t.LE.tf) GOTO 5
			WRITE(*,*) 'END OF TRANSIENT ANALYSIS'
			RETURN
		ENDIF
	ENDIF

c	Check the tolearance limit
c	==========================
	IF((zdf1.GT.emax.OR.zdf2.GT.emax.OR
     &	.zdf3.GT.emax).AND.kint.LT.12) GOTO 10  ! Iterate again

c	Convergence achieved for a particular time step
c	===============================================
	WRITE(*,555)t,kint,ak

c	Output Disk topology at specific time
c	=====================================
	IF(DABS(t-tout1).LT.dt.OR.DABS(t-tout2).LT.dt
     &	.OR.DABS(t-tout3).LT.dt) CALL disktop

c	Advance time step
c	=================
	tprev = t
	t=t+dt

c	Suspension And Slider update
c	============================
	CALL advance_susp
	
c     Disk update
c     ===========
      IF (luls_opt.EQ.0) THEN

c       Disk displacement      
	  WRITE(71,101)t,-acc_val,w_trn,Dw_x_susp,Dw_y_susp,ramp_cntf

c       Full disk displacement		  
	  IF ((disk_opt .EQ. 1) .AND.(data_opt_disk.EQ.1)) THEN
		    w_disk(1) = t_val
            w_disk(2:no_disk_dofs+1) = w_disk_cur(1:no_disk_dofs)
            WRITE(72,'(25000(E16.9,3X))') 
     &          (w_disk(i),i=1,no_disk_dofs+1)
		ENDIF
		
        IF ((supp_opt.NE.0) .AND. (data_opt_supp.EQ.1)) THEN
        
c           Disk Hub interface displacement   	
   	      w_disk_intf(1) = t_val
            DO i = 1 , no_disk_intf_dofs
                DO j = 1 , 3
                    w_disk_intf((3*(i-1))+j+1) = disk_ID_disp(i,j)  
                ENDDO
            ENDDO
            WRITE(73,'(100(E16.9,3X))') 
     &          (w_disk_intf(i),i=1,(3*no_disk_intf_dofs)+1)   	
   	
c           FDB displacement
            WRITE(74,'(25000(E16.9,3X))') 
     &          (w_FDB(i),i=1,no_fdb_dofs_op+1)
   	  
        ENDIF		
		
	ENDIF

c     Update Parameters
c     =================
      h(1:nx,1:ny) = hnew(1:nx,1:ny)
      pold(1:nx,1:ny) = p(1:nx,1:ny)

c     Update actuator position
c     ========================	
	CALL adv_actuator
	
c     Compute Bearing Number
c     ======================	
	IF((omega.NE.0.d0) .AND. ((nap.NE.0).OR.(nsp.NE.0))) THEN
		CALL calc_bearing_number
	ENDIF

c     Displacement at POI (4 points) from disk mean surface
c     =====================================================	  
      ssk=dsin(ske)
	csk=dcos(ske)
	dfx=dfx+ra*dt*omega/xl
	dfy=ra/xl
	
	DO i=1,4
		xtemp=xg-xint(i)
		ytemp=yint(i)-(0.5d0*yl+yg)
		xloc=dfx+(xtemp*csk-ytemp*ssk)
		yloc=dfy+(xtemp*ssk+ytemp*csk)
		CALL get_floor(xloc,yloc,hloc(i))
		flht(i)=hmin*hm+hp*(1.d0-xint(i))*xl-
     &			(hy*hm/xl)*(0.5d0*yl-yint(i))*xl+
     &			hm*thh(i)
	ENDDO

c	Calculate min FH, min spacing, xloc & yloc
c	==========================================
      flhtmin=100.0d0
      xflhtmin=0.0d0
      yflhtmin=0.0d0
	spacemin=100.0d0
	xspacemin=0.0d0
	yspacemin=0.0d0
      DO i=1,nx
		DO j=1,ny
			xloc=xref(i)
			yloc=yref(j)
			CALL pointrecess(xloc,yloc,tempthh)
			tflht=hmin*hm+hp*(1.d0-xloc)*xl-
     &			(hy*hm/xl)*(0.5d0*yl-yloc)*xl+
     &			hm*tempthh
			IF (tflht.LT.flhtmin) THEN
				flhtmin=tflht
				xflhtmin=xloc
				yflhtmin=yloc
			ENDIF

c           Now calculate minimum spacing
c	      =============================
			xtemp=xg-xref(i)
			ytemp=yref(j)-(0.5d0*yl+yg)
			xfloc=dfx+(xtemp*csk-ytemp*ssk)
			yfloc=dfy+(xtemp*ssk+ytemp*csk)
			CALL get_floor(xfloc,yfloc,temphloc)
			tspace=tflht+(hm*temphloc)
			IF(tspace.LT.spacemin) THEN
				spacemin=tspace
				xspacemin=xloc
				yspacemin=yloc
			ENDIF
			spa(i,j)=tspace*1.d0+0.12d-9
          ENDDO
      ENDDO

c	Find the point of minimum flying height
c	=======================================
c     ----------------------------------------------
c     hnew contains the speration of each grid point
c     ----------------------------------------------
	uuww=1d30
	DO i=1,nx
		DO j=1,ny
			IF (hnew(i,j).LT.uuww) THEN
				uuww=hnew(i,j)
				iuu=i
				juu=j
			ENDIF
		ENDDO
	ENDDO

c     POI data
c     ========
	WRITE(41,101)t,flht(1),flht(2),flht(3),flht(4)
	WRITE(42,101)t,-hm*hloc(1),-hm*hloc(2),-hm*hloc(3),
     &				 -hm*hloc(4)

c     Flying data
c     ===========     
      WRITE(51,101)t,hmin*hm,hp,-hyy,zc*hm,
     &	uuww*hm,xref(iuu)*xl,yref(juu)*xl
      WRITE(52,101)t,f,xf,yf,fp*9.81d0,ABS_data(1),ABS_data(2),
     &  fn*9.81d0,ABS_data(3),ABS_data(4)
	WRITE(53,101)t,f_maxcrp,f_maximp
	WRITE(54,101)t,fcr,txr,tyr,aratio
	WRITE(55,101)t,fcont,fsctx,fscty,aratio1
	WRITE(56,101)t,fvdw,yfvdw,xfvdw
	WRITE(57,101)t,felecst,telecst,relecst      

c	Output Pressure
c	===============
	IF (outp_opt.EQ.1) THEN
		DO i=1,4
			IF ((t.GE.toutp(i)).AND.(INT(foutp(i)).NE.1)) THEN
				CALL wr_pres(i)     ! Write Pressure data
				foutp(i) = 1.0d0    ! Turn on the flag 
			ENDIF
		ENDDO	
	ENDIF

	IF(t.LE.tf) GOTO 5
	WRITE(*,*) 'END OF TRANSIENT ANALYSIS'
	      
c	Close disk data files
c	=====================

c     Disk displacement
      CLOSE(71)
      
c     Full disk displacement
	IF ((disk_opt.EQ.1).AND.(luls_opt.EQ.0).AND.
     &	    (data_opt_disk.EQ.1)) THEN
        CLOSE(72)
	ENDIF
	
	IF ((supp_opt.NE.0) .AND. (data_opt_supp.EQ.1)) THEN
        
c       Disk Hub interface displacement
   	  CLOSE(73)
   	
c       FDB displacement
        CLOSE(74)
   	  
      ENDIF

 101  FORMAT(E16.9,9(3X,E16.9))
 445	FORMAT(/' Adaptive time step due to convergence at t = ',
     &		 G12.6,'(s)'/)
 555	FORMAT('	  t(s)= ',E12.6,', itns=',I2,
     &		 ', rsdl= ',E12.6)
	
	RETURN
	END
	  
c======================================================================
	SUBROUTINE dynamic_solver_contact(kwrite)
c======================================================================

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
      USE vdwf_data
      USE elec_stat
      USE POI_cords
      USE luls_data
      USE shck_data
      USE cmpt_dfmn
      USE susp_data
      USE disk_data
      USE disk_aray
      USE base_data
      
	IMPLICIT REAL*8(a-h,o-z) 

c     Local Variables
c     ===============	

	LOGICAL crash
		
	DOUBLE PRECISION flht(4),hloc(4)
	DOUBLE PRECISION ufcnt(8)
	
c     Coarse time stepping data
c     =========================
	orfact = 1.5d0
	dtold  = dt

	kfinal = INT(f_dtfac)	! Number of small time steps
	dt = dtold/f_dtfac

c     Intialize
c     =========
	ufcnt = 0.0d0
	iuu   = 1000

c	Suspension matrices LU decomposition 
c	====================================
c	------------------------------------
c	Only needed when entering fine time 
c	step, new LU fact will be done after 
c	every advance_susp for the new susp
c	configuration
c	------------------------------------
	CALL susp_block_fact(2)
	
	ksub=0	        ! Initialize no. of fine steps

5	CONTINUE			! For time step

	ksub = ksub + 1   ! Fine time step number

c     Intialize normalized error
c     ==========================
	zdf1=10000.d0
	zdf2=10000.d0
	zdf3=10000.d0

c     Check for crash in starting
c     ===========================  
	CALL gethx
	CALL get_toth(crash,iic,jjc)
	IF(crash) THEN
		CALL calc_impact_force
	ELSE
		fcont=0.d0
		fsctx=0.d0
		fscty=0.d0
		aratio1=0.d0
	ENDIF
		
c     Shock magnitude (inertia loading : a)
c	====================================
      t_val = t+dt
      CALL compute_shock(t_val)

c	Spinning Disk Model
c	===================
	IF ((disk_opt.EQ.1) .AND. (luls_opt.EQ.0)) THEN

		ssk=dsin(ske)
		csk=dcos(ske)
        
        w_trn  = disk_atd(ksub,1)
		Dw_x_susp = -(disk_atd(ksub,2)*ssk) - (disk_atd(ksub,3)*csk)
		Dw_y_susp =  (disk_atd(ksub,2)*csk) - (disk_atd(ksub,3)*ssk)
	
	ELSE
		w_trn  = 0.d0
		Dw_x_susp = 0.d0
		Dw_y_susp = 0.d0

	ENDIF
	
c     Shock vector for actuator frame
c     ===============================

c     Angle substended at VCM center	
	sldr_vcm_ang = DACOS((len_c2c**2 + xact**2 - ra**2)/
     &                     (2*len_c2c*xact))
c     Angle with +x_HDD axis
      sldr_hdd_ang = base_ang - sldr_vcm_ang	

	kint=0            ! Initialize

10	CONTINUE			! For Iteration step

	kint=kint+1       ! Iteration number

c     Store previous iteration flying attitude
c     ========================================
	zcpre   = zc
	hx0pre  = hx0
	hypre   = hy

c     Store previous iteration normalized errors
c     ==========================================
	zdf1pre = zdf1
	zdf2pre = zdf2
	zdf3pre = zdf3

	CALL calc_contact_force(0)
	IF(crash) THEN
	  cgap2=1d-10/hm
	  cgap3=2d-6/hm
	  DO i=1, nx
		    DO j=1, ny
				IF (hnew(i,j).LT.cgap2) THEN 
					hnew(i,j)=cgap3
				ENDIF 
		  
				IF (hydnew(i,j).LT.cgap2) THEN
					hydnew(i,j)=cgap3	
				ENDIF
		    
				IF (hyunew(i,j).LT.cgap2) THEN 
					hyunew(i,j)=cgap3	
				ENDIF
		    
				IF (hxlnew(i,j).LT.cgap2) THEN 
					hxlnew(i,j)=cgap3	
				ENDIF
		    
				IF (hxrnew(i,j).LT.cgap2) THEN
					hxrnew(i,j)=cgap3	
				ENDIF
			ENDDO
		ENDDO
	ENDIF
	
c	Solve Reynolds Equation
c	=======================
	IF(omega.NE.0.d0) CALL reyneq(0)
	CALL calc_ABS_force(0)
	CALL calc_electrostatic

c	Suspension displacment
c	======================
      CALL solve_susp(2)

c	Add disk attitude
c	=================
	zc  =  zc_susp - (w_trn/hm)
	hx0 = hx0_susp + (Dw_x_susp*xl/hm)
	hy  =  hy_susp + (Dw_y_susp*xl/hm)

c     Slider attitude
c     ===============
   	hp   = -hx0*hm/xl
	hyy  =  hy*hm/xl
	hmin =  zc+hx0*(1.d0-xg)

c     Check for contact
c     =================
	CALL gethx
	CALL get_toth(crash,iic,jjc)
	IF(crash) THEN
		CALL calc_impact_force
	ELSE
		fcont=0.d0
		f_maximp = 0.d0
		fsctx=0.d0
		fscty=0.d0
		aratio1=0.d0
	ENDIF

c     Normalized error
c     ================
	zdf1=DABS((zcpre-zc)/zc)
	zdf2=DABS((hx0pre-hx0)/hx0)
	zdf3=DABS((hypre-hy)/hy)

c     If solution doesn't converge try with new flying attitude 
c     computed using interpolation between current and previous step
c     ==============================================================
	IF(kint.GT.10) THEN
		IF((zdf1.GT.zdf1pre.AND.zdf1.GT.emax).OR.(zdf2.GT.zdf2pre.
     &		AND.zdf2.GT.emax).OR.(zdf3.GT.zdf3pre.AND.zdf3.GT.emax))
     &		THEN
			zc  = orfact*zc  + (1.d0-orfact)*zcpre 
			hx0 = orfact*hx0 + (1.d0-orfact)*hx0pre 
			hy  = orfact*hy  + (1.d0-orfact)*hypre 
			hmin= zc+hx0*(1.d0-xg)
		ENDIF
	ENDIF

c	Check the tolearance limit
c	==========================
	IF((zdf1.GT.emax.OR.zdf2.GT.emax.OR
     &	.zdf3.GT.emax).AND.kint.LT.12) GOTO 10

	WRITE(*,555)t,kint

c	Output Disk topology at specific time
c	======================================
	IF(DABS(t-tout1).LT.dt.OR.DABS(t-tout2).LT.dt
     &		  .OR.DABS(t-tout3).LT.dt) CALL disktop

c	Convergence achieved for a this time step
c	=========================================
	tprev = t
	t=t+dt

c	Update previous step values
c	===========================
	CALL advance_susp			! Suspension and Slider update

c	Suspension matrices LU decomposition
c	====================================
	CALL susp_block_fact(2)

	CALL calc_contact_force(0)
	IF(crash) THEN
		cgap2=1d-10/hm
		cgap3=2e-6/hm
		DO i=1, nx
			DO j=1, ny
				IF (hnew(i,j).LT.cgap2) THEN 
					hnew(i,j)=cgap3
				ENDIF 
				IF (hydnew(i,j).LT.cgap2) THEN
					hydnew(i,j)=cgap3	
				ENDIF
				IF (hyunew(i,j).LT.cgap2) THEN 
					hyunew(i,j)=cgap3	
				ENDIF
				IF (hxlnew(i,j).LT.cgap2) THEN 
					hxlnew(i,j)=cgap3	
				ENDIF
				IF (hxrnew(i,j).LT.cgap2) THEN 
					hxrnew(i,j)=cgap3	
				ENDIF
			ENDDO
		ENDDO
	ENDIF	

c	Solve Reynolds Equation
c	=======================
	IF(omega.NE.0.d0) CALL reyneq(0)
	CALL calc_ABS_force(0)
	CALL calc_electrostatic

c     Update Parameters
c     =================
      h(1:nx,1:ny) = hnew(1:nx,1:ny)
      pold(1:nx,1:ny) = p(1:nx,1:ny)

c	Suspension displacment
c	======================
      CALL solve_susp(2)

c	Add disk attitude
c	=================
	zc  =  zc_susp - (w_trn/hm)
	hx0 = hx0_susp + (Dw_x_susp*xl/hm)
	hy  =  hy_susp + (Dw_y_susp*xl/hm)

c     Slider attitude
c     ===============
   	hp   = -hx0*hm/xl
	hyy  =  hy*hm/xl
	hmin =  zc+hx0*(1.d0-xg)

c     Asperity contact force
c     ======================
	IF (fcr.GT.ufcnt(1)) THEN
		ufcnt(1)=fcr
		ufcnt(2)=txr
		ufcnt(3)=tyr
		ufcnt(4)=aratio
	ENDIF

c     Impact force
c     ============
	IF (fcont.GT.ufcnt(5)) THEN 
		ufcnt(5)=fcont
		ufcnt(6)=fsctx
		ufcnt(7)=fscty
		ufcnt(8)=aratio1
	ENDIF

	CALL gethx
	CALL get_toth(crash,iic,jjc)
	
c     Mininmum clearance computation
c     ==============================	
	IF (iuu.EQ.1000) THEN
	  uuww=1d30
	ENDIF
	DO i=1,nx
		DO j=1,ny
			IF (hnew(i,j).LT.uuww) THEN
				uuww=hnew(i,j)
				iuu=i
				juu=j
			ENDIF
		ENDDO
	ENDDO

	DO i=2,nx-1
		DO j=2,ny
			IF (hxlnew(i,j).LT.uuww) THEN
				uuww=hxlnew(i,j)
				iuu=i
				juu=j
			ENDIF
		ENDDO
	ENDDO

	DO i=2,nx
		DO j=2,ny-1
			IF (hydnew(i,j).LT.uuww) THEN
				uuww=hydnew(i,j)
				iuu=i
				juu=j
			ENDIF
		ENDDO
	ENDDO

c     Update actuator position
c     ========================	
	CALL adv_actuator

	IF(omega.NE.0.d0.AND.(nap.NE.0.OR.nsp.NE.0)) THEN
	  CALL calc_bearing_number
	ENDIF
	dfx=dfx+ra*dt*omega/xl
	dfy=ra/xl

c	Conditions states that write data in log files
c	only when done with all fine time steps as 
c     kwrite always equal to 0
c	==============================================
	IF((kwrite.EQ.0).AND.(ksub.NE.kfinal)) GOTO 123

	
c     Displacement at POI (4 points) from disk mean surface
c     =====================================================	 
      ssk=dsin(ske)
	csk=dcos(ske) 	
	DO i=1,4
		xtemp=xg-xint(i)
		ytemp=yint(i)-(0.5d0*yl+yg)
		xloc=dfx+(xtemp*csk-ytemp*ssk)
		yloc=dfy+(xtemp*ssk+ytemp*csk)
		CALL get_floor(xloc,yloc,hloc(i))
		flht(i)=hmin*hm+hp*(1.d0-xint(i))*xl-
     &			(hy*hm/xl)*(0.5d0*yl-yint(i))*xl+
     &			hm*thh(i)
	ENDDO

c	Calculate min FH and min FH loc
c	===============================
      flhtmin=100.0d0
      xflhtmin=0.0d0
      yflhtmin=0.0d0
      spacemin=100.0d0
      xspacemin=0.0d0
      yspacemin=0.0d0
      DO i=1,nx
		DO j=1,ny
			xloc=xref(i)
			yloc=yref(j)
			CALL pointrecess(xloc,yloc,tempthh)
			tflht=hmin*hm+hp*(1.d0-xloc)*xl-
     &            (hy*hm/xl)*(0.5d0*yl-yloc)*xl+
     &            hm*tempthh
			IF (tflht.LT.flhtmin) THEN
				flhtmin=tflht
				xflhtmin=xloc
				yflhtmin=yloc
			ENDIF

c     Now calculate minimum spacing
c	=============================
			xtemp=xg-xref(i)
			ytemp=yref(j)-(0.5d0*yl+yg)
			xfloc=dfx+(xtemp*csk-ytemp*ssk)
			yfloc=dfy+(xtemp*ssk+ytemp*csk)
			CALL get_floor(xfloc,yfloc,temphloc)
			tspace=tflht+(hm*temphloc)
			IF(tspace.LT.spacemin) THEN
				spacemin=tspace
				xspacemin=xloc
				yspacemin=yloc
			ENDIF
			spa(i,j)=tspace*1.d0+0.12d-9
		ENDDO
      ENDDO

c     POI data
c     ========	 								    
	WRITE(41,101)t,flht(1),flht(2),flht(3),flht(4)
	WRITE(42,101)t,-hm*hloc(1),-hm*hloc(2),-hm*hloc(3),
     &				 -hm*hloc(4)

c     Flying Dynamics
c     ===============
      WRITE(51,101)t,hmin*hm,hp,-hyy,zc*hm,
     &	uuww*hm,xref(iuu)*xl,yref(juu)*xl
      WRITE(52,101)t,f,xf,yf,fp*9.81d0,ABS_data(1),ABS_data(2),
     &  fn*9.81d0,ABS_data(3),ABS_data(4)
	WRITE(53,101)t,f_maxcrp,f_maximp
	WRITE(54,101)t,ufcnt(1),ufcnt(2),ufcnt(3),ufcnt(4)		   
	WRITE(55,101)t,ufcnt(5),ufcnt(6),ufcnt(7),ufcnt(8)	   
      WRITE(56,101)t,fvdw,yfvdw,xfvdw
      WRITE(57,101)t,felecst,telecst,relecst

c	Disk displacement
c	=================
      IF (luls_opt.EQ.0) THEN
        WRITE(71,101)t,-acc_val,w_trn,Dw_x_susp,Dw_y_susp,ramp_cntf
	  IF ((disk_opt.EQ.1) .AND. (data_opt_disk.EQ.1)) THEN
            w_disk(1) = t_val
            w_disk(2:no_disk_dofs+1) = w_disk_cur(1:no_disk_dofs)
            WRITE(72,'(25000(E16.9,3X))')
     &                (w_disk(i),i=1,no_disk_dofs+1)
        ENDIF    
	ENDIF

 123	IF(ksub.LT.kfinal) GOTO 5	! Complete small time steps
 
c     Updating time step size
c     =======================
      dt=dtold
 
 101  FORMAT(E16.9,9(3X,E16.9)) 
 555	FORMAT('	  t(s)= ',E12.6,', itns=',I2)
 
	RETURN
	END

c======================================================================
	SUBROUTINE compute_shock(t_shk)
c======================================================================

c----------------------------------------------------------------------
c     Computes the shock magnitude at precribed time : t_shk
c----------------------------------------------------------------------

c     Shared Data
c     =========== 
      USE shck_data
      
      IMPLICIT NONE
      
	DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979

c     Argument
c     ========	
	DOUBLE PRECISION t_shk
	
c     Local Variable	
c     ==============
      INTEGER i_shkpvt1,i_shkpvt2
      
c     Half sine wave
c     ==============
      IF (shk_mod.EQ.1) THEN
	    IF (t_shk.GT.T_st) THEN
		      IF (t_shk.LT.(T_st+T_sh)) THEN
			      acc_val = 9.81*acc_amp*
     &              DSIN((t_shk-T_st)*pi/T_sh)
	        ELSE
	            acc_val = 0.0d0
			  ENDIF
	    ELSE
            acc_val = 0.0d0
	    ENDIF

c     Square wave
c     ===========
      ELSEIF (shk_mod.EQ.2) THEN
          IF (t_shk.GT.T_st) THEN
			  IF (t_shk.LT.(T_st+T_sh)) THEN
                  acc_val = 9.81*acc_amp
			  ELSE
	            acc_val = 0.0d0
			  ENDIF
		  ELSE
		      acc_val = 0.0d0
	    ENDIF

c     Constant magnitude sine wave
c     ============================
	ELSEIF (shk_mod.EQ.3) THEN
	    IF (t_shk.GT.T_st) THEN
              acc_val = 9.81*acc_amp*
     &            DSIN((t_shk-T_st)*pi/T_sh)
	    ELSE
              acc_val = 0.0d0
	    ENDIF

c     Decaying sine wave
c     ==================
	ELSEIF (shk_mod.EQ.4) THEN
	  IF (t_shk.GT.T_st) THEN
		    acc_val = 9.81*acc_amp*exp(-t_shk/T_cnst)*
     &			dsin((t_shk-T_st)*pi/T_sh)
	  ELSE
			acc_val = 0.0d0
	  ENDIF
	    
c     User defined shock wave
c     =======================	    
	ELSEIF (shk_mod.EQ.5) THEN
	  i_shkpvt1 = shkpvt
	   
c       Current time out of prescribed list
c       ===================================
	  IF (shock_val(shkpts,1) .LT. t_shk) THEN 
    	      acc_val = 0.0d0

c       Pivot at current time i.e. t_shk
c       ================================	  
	  ELSEIF (shock_val(i_shkpvt1,1) .EQ. t_shk) THEN
	      acc_val = shock_val(i_shkpvt1,2)

c       Interpolate the value
c       =====================
	  ELSE  
	      DO WHILE ((shock_val(i_shkpvt1,1) .LT. t_shk) .AND. 
     &          (i_shkpvt1 .LT. (shkpts-1)))
                i_shkpvt1 = i_shkpvt1 + 1
            ENDDO
            IF (shock_val(i_shkpvt1,1) .EQ. t_shk) THEN
	          acc_val = shock_val(i_shkpvt1,2)
            ELSE
                i_shkpvt2 = i_shkpvt1
                i_shkpvt1 = i_shkpvt1 - 1
                acc_val = ((shock_val(i_shkpvt2,2)-
     &                      shock_val(i_shkpvt1,2))*
     &                      t_shk)+((shock_val(i_shkpvt2,1)*
     &                      shock_val(i_shkpvt1,2))-
     &                      (shock_val(i_shkpvt1,1)*
     &                      shock_val(i_shkpvt2,2)))       
                acc_val = acc_val/(shock_val(i_shkpvt2,1) 
     &                    - shock_val(i_shkpvt1,1))
	      ENDIF
	    ENDIF
	    
          acc_val = -acc_val
		  shkpvt = i_shkpvt1    
	ELSE
	    acc_val = 0.0d0
	ENDIF
      
      RETURN
      END

c======================================================================
	SUBROUTINE adv_actuator
c======================================================================																	  
													  
c----------------------------------------------------------------------															  
c	Finds the motion of the suspention given the acceleration 			  
c	profile of an actuator								  
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE siml_pram
      USE syst_cnfg																			  
      USE actr_data
      USE sldr_dynm
      	
	IMPLICIT REAL*8 (a-h,o-z) 											  

c     Compute disk angular velocity when nsp .NE. 0
c     =============================================
	IF(nsp.NE.0) THEN
		IF(t.LE.someg(1,1)) THEN
			omega=someg(2,1)
		ELSE IF(t.GE.someg(1,nsp)) THEN
			omega=someg(2,nsp)
		ELSE
			DO 10 i=1,nsp-1
				IF(t.GT.someg(1,i).AND.t.LE.someg(1,i+1)) THEN
					u1=someg(1,i)
					u2=someg(1,i+1)
					h1=someg(2,i)
					h2=someg(2,i+1)
					xmi=(h2-h1)/(u2-u1)										   
					ui=t-u1
					omega=xmi*ui+h1
					GOTO 50
				ENDIF
10			CONTINUE
		ENDIF
	ENDIF

50	CONTINUE

	IF(nap.EQ.0.OR.nap.EQ.1) THEN 
	  WRITE(65,102)t,ra,vact*xact,aact*xact,ske*360.d0/twopi,
     &	    skeeft*360.d0/twopi,omega*60.d0/twopi
	  RETURN
	ENDIF
																			  
c	Update actuator parameters : vact and aact											  
c	==========================================														  

c	Before 
c	======																		  
	IF(t.LE.hac(1,1)) THEN												  
	  aact=0.d0															 
		vact=vst															
		dact=vst*t+dst														
		IF(DABS(hac(1,1)-t).LE.dt/2.d0) dst=dact
		WRITE(65,102)t,ra,vact*xact,aact*xact,ske*360.d0/twopi,
     &	    skeeft*360.d0/twopi,omega*60.d0/twopi					   
		RETURN															 
	ENDIF																  
																			  
c	After
c	=====																		  
	IF(t.GT.hac(1,nap)) THEN												  
		aact=0.d0															 
		vact=vst															
		dact=vst*(t-hac(1,nap))+dst
		WRITE(65,102)t,ra,vact*xact,aact*xact,ske*360.d0/twopi,
     &      skeeft*360.d0/twopi,omega*60.d0/twopi									
		RETURN
	ENDIF																  
																			  
c	During
c	======																		  
	DO 176 i=1,nap-1														  
		IF (t.GT.hac(1,i).AND.t.LE.hac(1,i+1)) GO TO 170					
		GO TO  176															
170		u1=hac(1,i) 														
		u2=hac(1,i+1)														
		IF(DABS(u1-u2).LE.1.0d-10) GO TO 176						 
		h1=hac(2,i) 														
		h2=hac(2,i+1)														
		xmi=(h2-h1)/(u2-u1) 										 
		ui=t-u1 															
		aact=xmi*ui+h1													 
		vact=(ui*xmi/2.d0+h1)*ui+vst									   
		dact=((ui*xmi/6.d0+h1/2.d0)*ui+vst)*ui+dst						 
		IF(t+dt.GT.u2) THEN 											   
			ui=u2-u1														 
			dst=((ui*xmi/6.d0+h1/2.d0)*ui+vst)*ui+dst 					
			vst=(ui*xmi/2.d0+h1)*ui+vst							  
		ENDIF
			
		WRITE(65,102)t,ra,vact*xact,aact*xact,ske*360.d0/twopi,
     &      skeeft*360.d0/twopi,omega*60.d0/twopi																
		RETURN															 
176	CONTINUE
     
102	FORMAT(E16.9,6(3X,E16.9))
														  
	RETURN																  
	END															 
c======================================================================