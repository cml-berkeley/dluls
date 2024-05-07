c======================================================================
	PROGRAM response
c======================================================================

CC*********************************************************************
c    lulshock.F 
c    It's based on d510.f 
c    Last modified by Rahul Rai
c    Main changes:
c	 1. Rotating Disk
c      2. Ramp Disk Contact
CC*********************************Log*********************************
c	 09/23/2009 : Can handle 998 X 998 grid & 200 rail points
c	 09/28/2009 : Ramp Contact on both side
c	 10/04/2009 : solve_disk in fine time stepping was corrected
c	 10/08/2009 : suspension displacement corrected
c	 10/28/2009 : New LU factorization scheme for suspension
c	 11/07/2009 : Correction in t_val for disk displacement 
c      03/03/2010 : Modification to wasp,aasp subroutines for IVF_x64
c      03/05/2010 : Disk ramp contact force log in shock.dat
c      03/07/2010 : Uniform output data format
c      05/03/2010 : Contact model update (icmod = 0 => no contact)
c      03/11/2010 : User defined shock profile
c      05/20/2010 : Unused contact block variables elimnated
c      05/24/2010 : Output for grid point recess from slider ref. plane
c      05/28/2010 : Change the input file format dynluls.def
c      05/28/2010 : Change COMMON blocks : suspension,dlul1 to dlul4
c      06/02/2010 : Added 3 subroutines for new rail/wall algorithm
c      06/02/2010 : Array size corrected for maxRailPts
c      06/06/2010 : Mod. subroutine : read_rail,norm_rail,wall_profile
c      06/07/2010 : New subroutine for pointrecess
c      06/07/2010 : This code do not support Quick300 format anymore
c      06/09/2010 : New Mod : luls_size and susp_data
c      06/09/2010 : suspension arrays : variable size
c      06/10/2010 : New Mod : sldr_stat,luls_data,shck_data,cmpt_dfmn
c      06/10/2010 : New Mod : rail_data,wall_data,sldr_dynm,impc_data
c      06/10/2010 : New Mod : siml_pram,syst_cnfg,actr_data,grid_ctrl
c      06/11/2010 : New Mod : aspr_data,vdwf_data,rynl_data,POI_cords  
c      06/12/2010 : New Mod : sldr_grid,sldr_rghn,disk_rghn,disk_fltr
c      06/12/2010 : New Mod : sldr_arys,trck_prfl
c      06/12/2010 : New Mod : disk_data
c      06/13/2010 : New Mod : disk_aray; susp nb & damp vars changed
c      06/13/2010 : New Subroutine : compute_shock 
c      06/13/2010 : New Subroutine : init_disk
c      06/13/2010 : Mod. Sub. : gen_disk_mesh and gen_disk_matrices
c      06/13/2010 : Mod. Sub. : solve_disk,disk_attitude,ramp_cnstrnt
c      06/14/2010 : disk_solve : No disk arrays images
c      06/14/2010 : gen_disk_matrices : omg initialized
c      06/14/2010 : Common disk_map for mapping disk sparse matrices
c      06/14/2010 : Common time stepping for coarse and fine matrices
c      06/14/2010 : Mod. Mod. : syst_cnfg,actr_data
c      06/15/2010 : Suspension variables names changed
c      06/16/2010 : Disk roughness matrices : Fixed to allocatable
c      06/17/2010 : Slider roughness matrices : Fixed to allocatable
c      06/17/2010 : Removed nramax, nrpmax and nfxmax
c      06/17/2010 : New Sub. : get_rail_size
c      06/17/2010 : rail_data and wall_data : Allocatable
c      06/18/2010 : equate : eliminated; tridag, interp : modified
c      06/18/2010 : sldr_grid : Allocatable
c      06/18/2010 : ixtr,jytr removed from calc_ABS/contact_force
c      06/19/2010 : sldr_arys and rynl_data : Allocatable
c      06/19/2010 : max_nx/ny replaced nx/nymax for sub. used once
c      06/19/2010 : acmdn and acmup : Modified with allocatable arrays
c      06/19/2010 : New : allocate_reynolds_arrays
c      06/19/2010 : All arrays allocatable
c      06/20/2010 : Adaptive Routine Fixed (interp)
c      06/20/2010 : Several SUBROUTINES updated 
c      06/30/2010 : Unused SUBROUTINES removed, some renamed
c      06/30/2010 : LU_sparse_fact removed
c      06/30/2010 : Stoped re-intialization of susp_block_fact matrices
c      07/11/2010 : Modified DynLULS.def to include air parameters
c      07/13/2010 : Tempreature dependent air parameters
c      07/13/2010 : Electrostatic Forces added
c      07/20/2010 : Output file rearranged
c      09/01/2010 : Subroutine for removing collinear points from rails
c      09/01/2010 : Disabled shock during Load and Unload mode
c      09/20/2010 : Removed redundant slider attitude computation
c      09/22/2010 : Same value of eycrash,f_nu,ydcrash for LUL & shock.
c      09/26/2010 : Using disk interpolation during fine time step
c      01/11/2011 : New input data files and maths functions
c      01/29/2011 : Spindle motor integration with new disk modules
c      01/21/2011 : Removed output message from spindle motor module
c      02/03/2011 : Options for FDB and Disk Interface data output
c      04/05/2011 : Base plate integration with oblique shock
c      04/05/2011 : Oblique shock for actuator in rotating frame
c      04/06/2011 : Slider angle corrected for HDD frame
c      04/07/2011 : len_ramp has been updated to x_ramp
c      04/12/2011 : w_FDB/disk_intf is fixed for supt_opt
CC*********************************************************************

c     Shared Data
c     ===========
      USE siml_pram
      USE syst_cnfg
      USE actr_data
      USE rail_data
      USE wall_data
      USE grid_ctrl
      USE sldr_grid
      USE sldr_arys
      USE sldr_rghn
      USE disk_rghn
      USE disk_fltr
      USE trck_prfl
      USE rynl_data
      USE rynl_arys
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

	IMPLICIT REAL*8(a-h,o-z)
	 
c     Local Variable
c     ==============
	LOGICAL crash
	INTEGER i,j
	INTEGER t_start(8),t_end(8)
	
c======================================================================
c     Add code for expiration date
c======================================================================
      INTEGER START_TIME(8)
      CHARACTER (LEN = 12) REAL_CLOCK (3)
      CALL DATE_AND_TIME (REAL_CLOCK (1), REAL_CLOCK (2), 
     &                    REAL_CLOCK (3), START_TIME)
	IF (START_TIME(1).GE.2018) then
	  IF (START_TIME(2).GE.9)then
	     WRITE(*,*) 'This solver has expired'
    	     WRITE(*,*) 'Please download the latest copy of CMLAir'
	     STOP
	  ENDIF
	ENDIF
c======================================================================
c     End Code for expiration
c======================================================================

	CALL DATE_AND_TIME (values=t_start)

	WRITE(*,*)
	WRITE(*,*)"++++++++++++++++++++++++++++++++"
	WRITE(*,*)"| CML LU Shock Simulator v 6.4 |"
	WRITE(*,*)"++++++++++++++++++++++++++++++++"
	WRITE(*,*)	

c	Initialize Workspace
c	====================
	CALL initws
	CALL allocate_slider_arrays
	CALL allocate_reynolds_arrays

	WRITE(*,*)		  
	WRITE(*,*) 'PRE - PROCESSING'
	WRITE(*,*) '================'

c	F-K Boltzman model
c	==================
	IF (iqpo.EQ.5) THEN
		CALL create_dbase
		WRITE(*,*) 'FK Database'
	ENDIF

c	No Pre meshed Grid Avialable
c	============================
	IF(ioldgrid.EQ.0) THEN
		CALL get_grids
	ENDIF

c	Snap the grid points to rails
c	=============================
	IF(ioldgrid.EQ.0.OR.iadapt.EQ.1) THEN
		CALL snapgrid
	ENDIF

	WRITE(*,*) 'Dense Meshing'
	CALL gethx1		 
	CALL get_slider 
	CALL fl_ms
	CALL calc_bearing_number

c	Adaptive Meshing
c	================
	IF(iadapt.EQ.1) THEN
		WRITE(*,*) 'Adaptive mesh discretization'
		t=0.0d0
		CALL gethx
		CALL get_toth(crash,iic,jjc)
		IF(crash) THEN
			WRITE(*,250)xref(iic)*xl,yref(jjc)*xl
250			FORMAT(5X,'Initial FH: crash at
     &			  (',E16.9,'m;',E16.9,'m)')
			STOP
		ENDIF
		CALL reyneq(1)
		CALL adaptive
	ENDIF

c     Slider recess
c     =============
c     CALL wr_grid_recess

c     Write the new mesh
c     ==================	
	CALL wr_new_mesh
	
c     Disk Topography
c     ===============
	IF((nf_wave.NE.0) .OR. (nf_asper.NE.0)) THEN
	  CALL dtopog
	ENDIF 

c	Initialize Suspension
c	=====================
	WRITE(*,*) 'Initialize : Suspension'
	CALL init_susp
	
c     Disk Support System Modeling
c     ============================
	IF ((disk_opt.EQ.1) .AND. (luls_opt.EQ.0)) THEN

c	  Initialize Disk 
c	  ===============
        WRITE(*,*) 'Initialize : Rotating Disk'
        CALL init_disk
        
c       Mount the disk on the rotor
c       ===========================
        IF (supp_opt .NE. 0) THEN
            WRITE(*,*) 'Initialize : Rotating Hub'
            CALL mount_disk
        ENDIF
        
c       Initialize the stationary housing
c       =================================
        IF (supp_opt .EQ. 2) THEN
        
            WRITE(*,*) 'Initialize : Housing'
            CALL init_housing
            WRITE(*,*) 'Coupling   : Rotor and Housing'
            CALL couple_rotr_hsng
            
        ENDIF 
        
c       Initialize the stationary housing
c       =================================
        IF (supp_opt .EQ. 3) THEN
        
            WRITE(*,*) 'Initialize : Housing'
            CALL init_housing
            WRITE(*,*) 'Initialize : Base'
            CALL init_base
            WRITE(*,*) 'Coupling   : Housing and Base Plate'
            CALL couple_hsng_base
            WRITE(*,*) 'Coupling   : Rotor and Housing'
            CALL couple_rotr_base
        
        ENDIF                   
        
c       Dynamic set of matrices for disk supporting system
c       ==================================================      
        IF (supp_opt .EQ. 0) THEN         ! Fixed Disk
            CALL newmark_disk   
        ELSEIF (supp_opt .EQ. 1) THEN     ! Fixed Housing
            CALL newmark_rotor
        ELSEIF (supp_opt .EQ. 2) THEN     ! Full Spindle Motor
            CALL newmark_motor
        ELSEIF (supp_opt .EQ. 3) THEN     ! Full HDD Modal
            CALL newmark_support
        ENDIF     
        
	ENDIF

c	Old Dynamic Setup : Mode Superposition
c	======================================
	IF(msusp.EQ.1) THEN
		WRITE(*,*) 'Mode superposition not supported'
		WRITE(*,*) 'Program aborted'
		STOP
	ELSE
		CALL dynamic_solver	                ! Current Method of solution
	ENDIF

c	Write to data files
c	===================
	CALL wr_final_data
    
	CALL DATE_AND_TIME (values=t_end)

c	Execution Time
c	==============
c	OPEN(501,FILE='duration.out',STATUS='UNKNOWN')
c	REWIND(501)
c
c  	WRITE(501,'(A28,(4(I6,3X)))') 
c     &	   't_start (h:m:s:ms) : ' ,
c     &		t_start(5),t_start(6),
c     &		t_start(7),t_start(8)
c  	WRITE(501,'(A28,(4(I6,3X)))') 
c     &	   't_end   (h:m:s:ms) : ' ,
c     &		t_end(5),t_end(6),
c     &		t_end(7),t_end(8)
c
c	CLOSE(501)
	
	END
c======================================================================