c**********************************************************************
c     Subroutines in this file :
c     1. initws
c     2. allocate_slider_arrays
c     3. allocate_reynolds_arrays
c**********************************************************************

c======================================================================
	SUBROUTINE initws
c======================================================================

c----------------------------------------------------------------------
c	Decription
c	==========
c	Read in dynluls.def and initialize
c     Suspension and Disk
c----------------------------------------------------------------------

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
      USE sldr_dynm
      USE sldr_stat
      USE aspr_data
      USE impc_data
      USE vdwf_data
      USE elec_stat
      USE POI_cords
      USE luls_data
      USE shck_data
      USE susp_data
      USE disk_data
      USE disk_aray
      USE motr_data
      USE rotr_data
      USE hsng_data
      USE base_data
      
	IMPLICIT REAL*8(a-h,o-z)
	
	DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979

	CHARACTER*50 fname
	CHARACTER*70 line
	
c	Local Variables
c	===============
      INTEGER i,j,i_temp(10)
	INTEGER i_dtfac,no_shkpts,eof_flag
	DOUBLE PRECISION Ma,Na,MolecularDiameter
	DOUBLE PRECISION x_cent,tht_sldr,gam_sldr 
      DOUBLE PRECISION x_ramp,tht_ramp,gam_ramp
	DOUBLE PRECISION a(3,3),b(3),temp_val,skw
	DOUBLE PRECISION A8(3,3),AINV(3,3),A9(4,4),AINVER(4,4)
	
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xref_temp
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: yref_temp	
	
c**********************************************************************
c***********************Input from : DynLULS.def***********************
c**********************************************************************

      OPEN(02,ERR=999,FILE='dynluls.def',STATUS='OLD')														   
	twopi=8.d0*DATAN(1.d0)
	REWIND(02) 

c************************Problem Definition Menu***********************

	READ(2,*)
	READ(2,*)
	READ(2,*) xl,yl,xg,yg,zg

	READ(2,*)
	READ(2,*) f0,xf0,yf0

	READ(2,*)
	READ(2,*) halt,rpm,dt,tf
	
	READ(2,*)
	READ(2,*) ra,ra_if,ra_of
  
	xf0=xg
	yf0=yg
	dtorig = dt
	
c****************************Suspension Menu***************************

	READ(2,*)
	READ(2,*)
	READ(2,*) iact,xact,dact,vact,skeID
	
c	The IDEMA skew standard issue	
c     =============================
	skw = skeID
	ske = -skeID

c************************Initial Condition Menu************************

	READ(2,*)
	READ(2,*)
	READ(2,*) hm,hp,hsID,xdisp,ydisp,zang
	 
	hs = -hsID
	xdisp = xdisp*1.0d3
      ydisp = ydisp*1.0d3

	hmin=hm/1.0d-08
	hm=1.0d-08
	
c************************Solution Control Menu*************************

	READ(2,*)
	READ(2,*)
	READ(2,*) iqpo,akmax,emax,idisc

c**************************Grid Control Menu***************************

c	Input grid control parameters
c     =============================
	READ(2,*)
	READ(2,*)
	READ(2,*) iadapt,isymmetry,ioldgrid,nx,ny,nsx,nsy
	READ(2,*) !READ end coordinate of each x grid section
	READ(2,*)(xnt(i),i=2,nsx)
	READ(2,*) !READ end grid number of each x grid section
	READ(2,*)(nxt(i),i=2,nsx)
	READ(2,*) !READ expansion ratio of each x grid section
	READ(2,*)(dxr(i),i=1,nsx)
	READ(2,*) !READ end coordinate of each y grid section
	READ(2,*)(ynt(i),i=2,nsy)
	READ(2,*) !READ end grid number of each y grid section
	READ(2,*)(nyt(i),i=2,nsy)
	READ(2,*) !READ expansion ratio of each y grid section
	READ(2,*)(dyr(i),i=1,nsy)
	READ(2,*)
	READ(2,*) difmax,decay,ipmax

c********************Measured Track Profile Menu***********************

	READ(2,*)
	READ(2,*)
	READ(2,*) inwave,nfx,dfinit

c     Allocate and Initialize
c     =======================
	ALLOCATE(xfref(nfx))
	ALLOCATE(hfnew(nfx))
	xfref = 0.0d0
	hfnew = 0.0d0
	
c******************Disk Topography Gneralization Menu******************

	READ(2,*)
      READ(2,*)
	READ(2,*) nf_wave,nf_zone,nf_asper
	
c	(1) Wavy Roughness
c     ==================
	ALLOCATE(ftype(nf_wave))
	ALLOCATE(fmag(nf_wave))
	ALLOCATE(fang(nf_wave))
	ALLOCATE(fxwlth(nf_wave))
	ALLOCATE(fywlth(nf_wave))
	ALLOCATE(fxdm(nf_wave))
	ALLOCATE(fydm(nf_wave))
	ALLOCATE(fras(nf_wave))
	ALLOCATE(frae(nf_wave))
	ALLOCATE(flloc(nf_wave,nwtp4))
	ALLOCATE(flhgt(nf_wave,nwtp4))
	
	ftype  = 0.0d0
	fmag   = 0.0d0
	fang   = 0.0d0
	fxwlth = 0.0d0
	fywlth = 0.0d0
	fxdm   = 0.0d0
	fydm   = 0.0d0
	fras   = 0.0d0
	frae   = 0.0d0
	flloc  = 0.0d0
	flhgt  = 0.0d0

	READ(2,*)
	DO k=1,nf_wave
		READ(2,*)ftype(k),fmag(k),fang(k),fxwlth(k),fywlth(k),
     &			 fxdm(k),fydm(k),fras(k),frae(k)

	  fmagk=fmag(k)
		fmag(k)=fmag(k)/hm
		fang(k)=fang(k)*twopi/360.d0
		fxwlth(k)=fxwlth(k)/xl
		fywlth(k)=fywlth(k)/xl
		fxdm(k)=fxdm(k)/xl
		fydm(k)=fydm(k)/xl
		fras(k)=fras(k)/xl
		frae(k)=frae(k)/xl

		IF(ftype(k).EQ.4.d0) THEN

		  IF(DABS(fmagk-1.0d0).LE.1.0d-5) THEN
		    fname = 'w1.dat'
		    CALL wasp(133,fname,k)
		  ELSE IF(DABS(fmagk-2.0d0).LE.1.0d-5) THEN
		    fname = 'w2.dat'
		    CALL wasp(134,fname,k)
		  ELSE IF(DABS(fmagk-3.0d0).LE.1.0d-5) THEN
		    fname = 'w3.dat'		  
		    CALL wasp(135,fname,k)
		  ELSE IF(DABS(fmagk-4.0d0).LE.1.0d-5) THEN
		    fname = 'w4.dat'		  
		    CALL wasp(136,fname,k)
		  ELSE IF(DABS(fmagk-5.0d0).LE.1.0d-5) THEN
		    fname = 'w5.dat'		  
		    CALL wasp(137,fname,k)
		  ELSE IF(DABS(fmagk-6.0d0).LE.1.0d-5) THEN
		    fname = 'w6.dat'		  
		    CALL wasp(138,fname,k)
		  ELSE IF(DABS(fmagk-7.0d0).LE.1.0d-5) THEN
		    fname = 'w7.dat'
		    CALL wasp(139,fname,k)
		  ELSE IF(DABS(fmagk-8.0d0).LE.1.0d-5) THEN
		    fname = 'w8.dat'		  
		    CALL wasp(140,fname,k)
		  ELSE IF(DABS(fmagk-9.0d0).LE.1.0d-5) THEN
		    fname = 'w9.dat'		  
		    CALL wasp(141,fname,k)
		  ELSE IF(DABS(fmagk-10.0d0).LE.1.0d-5) THEN
		    fname = 'w10.dat'
		    CALL wasp(142,fname,k)
		  ENDIF

		ENDIF
	ENDDO

c     (2) Circumferential Zone Profile
c     ================================
      ALLOCATE(frp1(nf_zone))
      ALLOCATE(frp2(nf_zone))
      ALLOCATE(frp3(nf_zone))
      ALLOCATE(fcoe0(nf_zone))
      ALLOCATE(fcoe1(nf_zone))
      ALLOCATE(fcoe2(nf_zone))
      
      frp1 = 0.0d0
      frp2 = 0.0d0
      frp3 = 0.0d0
      
      fcoe0 = 0.0d0 
      fcoe1 = 0.0d0
      fcoe2 = 0.0d0
      
	READ(2,*)
	DO k=1,nf_zone
	  READ(2,*)frp1(k),fht1,frp2(k),fht2,frp3(k),fht3
	  frp1(k)=frp1(k)/xl
	  frp2(k)=frp2(k)/xl
	  frp3(k)=frp3(k)/xl
	  fht1=fht1/hm
		fht2=fht2/hm
	  fht3=fht3/hm
	  aa1=fht1/(frp1(k)-frp2(k))/(frp1(k)-frp3(k))
		aa2=fht2/(frp2(k)-frp1(k))/(frp2(k)-frp3(k))
		aa3=fht3/(frp3(k)-frp1(k))/(frp3(k)-frp2(k))
	  fcoe2(k)=aa1+aa2+aa3
	  fcoe1(k)=-aa1*(frp2(k)+frp3(k))
     &			 -aa2*(frp1(k)+frp3(k))
     &			 -aa3*(frp1(k)+frp2(k))
	  fcoe0(k)=aa1*frp2(k)*frp3(k)+aa2*frp1(k)*frp3(k)+
     &			 aa3*frp1(k)*frp2(k)
	ENDDO

c     (3) Asperities
c     ==============
      ALLOCATE(fatype(nf_asper))
      ALLOCATE(famag(nf_asper))
      ALLOCATE(faang(nf_asper))
      ALLOCATE(fxloc(nf_asper))
      ALLOCATE(fyloc(nf_asper))
      ALLOCATE(faxdm(nf_asper))
      ALLOCATE(faydm(nf_asper))
      ALLOCATE(falloc(nf_asper,natp4))
      ALLOCATE(falhgt(nf_asper,natp4))
      
      fatype = 0.0d0
      famag  = 0.0d0
      faang  = 0.0d0
      fxloc  = 0.0d0
      fyloc  = 0.0d0
      faxdm  = 0.0d0
      faydm  = 0.0d0
      falloc = 0.0d0
      falhgt = 0.0d0
      
	READ(2,*)
	DO k=1,nf_asper
	  READ(2,*)fatype(k),famag(k),faang(k),fxloc(k),fyloc(k),
     &			 faxdm(k),faydm(k)
     
	  famagk=famag(k)
		famag(k)=famag(k)/hm
		faang(k)=faang(k)*twopi/360.d0
		fxloc(k)=fxloc(k)/xl
		fyloc(k)=fyloc(k)/xl
		faxdm(k)=faxdm(k)/xl/2.d0
		faydm(k)=faydm(k)/xl/2.d0
		
		IF(fatype(k).EQ.4.d0) THEN
			IF(DABS(famagk-1.0d0).LE.1.0d-5) THEN
			    fname = 'a1.dat'	
			    CALL aasp(33,fname,k)
			ELSE IF(DABS(famagk-2.0d0).LE.1.0d-5) THEN
			    fname = 'a2.dat'
			    CALL aasp(34,fname,k)
			ELSE IF(DABS(famagk-3.0d0).LE.1.0d-5) THEN
			    fname = 'a3.dat'
			    CALL aasp(35,fname,k)
			ELSE IF(DABS(famagk-4.0d0).LE.1.0d-5) THEN
			    fname = 'a4.dat'
			    CALL aasp(36,fname,k)
			ELSE IF(DABS(famagk-5.0d0).LE.1.0d-5) THEN
			    fname = 'a5.dat'
			    CALL aasp(37,fname,k)
			ELSE IF(DABS(famagk-6.0d0).LE.1.0d-5) THEN
			    fname = 'a6.dat'
			    CALL aasp(38,fname,k)
			ELSE IF(DABS(famagk-7.0d0).LE.1.0d-5) THEN
			    fname = 'a7.dat'
			    CALL aasp(39,fname,k)
			ELSE IF(DABS(famagk-8.0d0).LE.1.0d-5) THEN
			    fname = 'a8.dat'
			    CALL aasp(40,fname,k)
			ELSE IF(DABS(famagk-9.0d0).LE.1.0d-5) THEN
			    fname = 'a9.dat'
			    CALL aasp(41,fname,k)
			ELSE IF(DABS(famagk-10.0d0).LE.1.0d-5) THEN
			    fname = 'a10.dat'
			    CALL aasp(42,fname,k)
			ENDIF
		ENDIF

	ENDDO
c	Input the 3 particular time for outputting the disk profile under
c	the slider
	IF(nf_wave.NE.0.OR.nf_asper.NE.0) THEN
		READ(2,*)tout1,tout2,tout3
	ENDIF
c
c*********Slider Wavy Roughness/Asperities Generalization Menu*********
c
	READ(2,*)
      READ(2,*)
	READ(2,*) ns_wave,ns_asper
	
	
c	(1) Wavy Roughness
c     ==================
      ALLOCATE(stype(ns_wave))
	ALLOCATE(smag(ns_wave))
	ALLOCATE(sang(ns_wave))
	ALLOCATE(sxwlth(ns_wave))
	ALLOCATE(sywlth(ns_wave))
	ALLOCATE(sxdm(ns_wave))
	ALLOCATE(sydm(ns_wave))
	
	stype  = 0.0d0
	smag   = 0.0d0
	sang   = 0.0d0
	sxwlth = 0.0d0
	sywlth = 0.0d0
	sxdm   = 0.0d0
	sydm   = 0.0d0
		
	READ(2,*)
	DO k=1,ns_wave
		READ(2,*)stype(k),smag(k),sang(k),sxwlth(k),sywlth(k),
     &			 sxdm(k),sydm(k)
		smag(k)=smag(k)/hm
		sang(k)=sang(k)*twopi/360.d0
		sxwlth(k)=sxwlth(k)/xl
		sywlth(k)=sywlth(k)/xl
		sxdm(k)=sxdm(k)/xl
		sydm(k)=sydm(k)/xl
	ENDDO

c	(2) Asperities
c     ==============
      ALLOCATE(satype(ns_asper))
      ALLOCATE(samag(ns_asper))
      ALLOCATE(saang(ns_asper))
      ALLOCATE(sxloc(ns_asper))
      ALLOCATE(syloc(ns_asper))
      ALLOCATE(saxdm(ns_asper))
      ALLOCATE(saydm(ns_asper))
      

      satype = 0.0d0
      samag  = 0.0d0
      saang  = 0.0d0
      sxloc  = 0.0d0
      syloc  = 0.0d0
      saxdm  = 0.0d0
      saydm  = 0.0d0
      
	READ(2,*)
	DO k=1,ns_asper
		READ(2,*)satype(k),samag(k),saang(k),sxloc(k),syloc(k),
     &			 saxdm(k),saydm(k)
		samag(k)=samag(k)/hm
		saang(k)=saang(k)*twopi/360.d0
		sxloc(k)=sxloc(k)/xl
		syloc(k)=syloc(k)/xl
		saxdm(k)=saxdm(k)/xl/2.d0
		saydm(k)=saydm(k)/xl/2.d0
	ENDDO

c*************************Track Accessing Menu*************************

	READ(2,*)
	READ(2,*)
	READ(2,*) nap
	READ(2,*)
	
	ALLOCATE(hac(2,nap))
	hac = 0.0d0
	DO i=1,nap
		READ(2,*) hac(1,i),hac(2,i) 		 
	ENDDO

c************************Slider Start/Stop Menu************************

	READ(2,*)
	READ(2,*)
	READ(2,*) nsp
	READ(2,*)
	
	ALLOCATE(someg(2,nsp))
	someg = 0.0d0
	
	DO i=1,nsp
		READ(2,*) someg(1,i),someg(2,i) 		 
	ENDDO

c*********************Disk Sinusoidal Flutter Menu*********************

	READ(2,*)
	READ(2,*)
	READ(2,*) kft,ft_ts,ft_tf,ft_omg,ft_mag
	ft_omg=ft_omg*twopi
	ft_mag=ft_mag/hm

c****************************Air Parameters****************************
      
      READ(2,*)
      READ(2,*)
      READ(2,*) p00,al0,amu0,tempCel
      
c************************Intermolecular Forces*************************

      READ(2,*)
      READ(2,*)
      READ(2,*) imf_opt,ahc,bhc,elecpot	

c*****************************Contact Menu*****************************

	READ(2,*)
	READ(2,*)
	READ(2,*) icmod,eyoung,ydst,pratio,frcoe
	READ(2,*)
	READ(2,*) ncz
	READ(2,*)
	DO k = 1,ncz	
		READ(2,*)rsikm(k),ceta(k),rasper(k),rcts(k),
     &	    rcte(k),gldht(k)
	ENDDO
	ydcoe = 1.282d0 + 1.158d0*pratio
	
c     Impact force material constants
c     ===============================
	eycrash = eyoung
	ydcrash = ydst
	f_nu = pratio
		
	IF(iact.EQ.0) msusp=0
	IF(inwave.EQ.0) dfinit=0.d0
	dfx=0.d0
	dfy=ra/xl

c***************************Load Unload Menu***************************
	

	READ(2,*) 
	READ(2,*)
	READ(2,*) luls_opt,outp_opt,toutp(1),toutp(2),toutp(3),toutp(4)
	
	DO i = 1,4
	  foutp(i) = 0.0d0    ! Flag for pressure output
	  toutp(i) = toutp(i)/1.0d3   ! Converting data ms -> sec
	ENDDO

c******************************Shock Menu******************************

      READ(2,*)
	READ(2,*)
	READ(2,*) shk_mod,shk_dir(1,1),shk_dir(2,1),shk_dir(3,1),
     &          T_st,T_sh,acc_amp,T_cnst
      acc_amp = -acc_amp    ! LULS code convention

c****************************Suspension Menu***************************

	READ(2,*)
	READ(2,*)
	READ(2,*) i_suspsize,i_dtfac,cp00,cr00
	f_dtfac = DBLE(i_dtfac)
	
c     Change the sign for skew : IDEMA
c     ================================
      cr00 = -cr00
      cz00 =  hm+hp*xl/2.0d0  
		
	READ(2,*)
	READ(2,*) i_dofux,i_dofuy,i_dofuz,i_dofrotx,
     &          i_dofroty,i_dofrotz,i_doftab
	READ(2,*)
	READ(2,*) i_nocele,i_dofflx,i_doffly,i_doflbx,i_doflby
	READ(2,*)
	DO i=1,i_nocele
		READ(2,*) i_dofcu(i),i_dofcd(i),i_constat(i),f_cedist(i),
     &            f_cepl(i)
      ENDDO
      
      READ(2,*)
      READ(2,*) susp_nb_1,susp_nb_2     ! Newmark beta parameters
      READ(2,*) 
      READ(2,*) susp_M_damp,susp_K_damp ! Structural Damping
         
c******************************Ramp Menu*******************************
      
      READ(2,*)
      READ(2,*)
      READ(2,*) ramp_opt,ramp_clrnc
      READ(2,*)
      READ(2,*) x_ramp,tht_ramp
	READ(2,*)
	READ(2,*) nrp
	ALLOCATE(f_rloc(2,nrp))
	READ(2,*)
	DO i=1,nrp
		READ(2,*) f_rloc(1,i),f_rloc(2,i) 		 
      ENDDO

c************************Support System Modeling***********************

      READ(2,*)      
      READ(2,*)
      READ(2,*) supp_opt,data_opt_supp 
      
c******************************Disk Menu*******************************

	READ(2,*)

c	Disk Options
c     ============
	READ(2,*)
	READ(2,*) disk_opt,data_opt_disk

c	Geometry
c     ========
	READ(2,*)
	READ(2,*) ri_disk,ro_disk,t_disk

c	Material Properties
c     ===================
	READ(2,*)
	READ(2,*) E_disk,nu_disk,den_disk

c	Mesh Data
c     =========
	READ(2,*)
	READ(2,*) no_rad_elem,no_ang_elem

c	Newmark Beta Parameters
c     =======================
	READ(2,*)
	READ(2,*) disk_nb_1,disk_nb_2
	
c	Numerical Damping
c     =================
	READ(2,*)
	READ(2,*) disk_M_damp,disk_K_damp

c**************************Spindle Motor Menu**************************

	READ(2,*)

c     Reduced system DOFs
c     ===================	
	READ(2,*)
	READ(2,*) no_hub_dofs,no_hsng_dofs
	
c     Journal Bearing Stiffness (Bottom)
c     ==================================
      READ(2,*)
      READ(2,*) K_fdb(1,1),K_fdb(1,2),K_fdb(1,3)
      READ(2,*) K_fdb(2,1),K_fdb(2,2),K_fdb(2,3)
      READ(2,*) K_fdb(3,1),K_fdb(3,2),K_fdb(3,3)
      
c     Journal Bearing Damping (Bottom)
c     ================================
      READ(2,*)
      READ(2,*) C_fdb(1,1),C_fdb(1,2),C_fdb(1,3)
      READ(2,*) C_fdb(2,1),C_fdb(2,2),C_fdb(2,3)
      READ(2,*) C_fdb(3,1),C_fdb(3,2),C_fdb(3,3)
      
c     Journal Bearing Stiffness (Top)
c     ===============================
      READ(2,*)
      READ(2,*) K_fdb(4,1),K_fdb(4,2),K_fdb(4,3)
      READ(2,*) K_fdb(5,1),K_fdb(5,2),K_fdb(5,3)
      READ(2,*) K_fdb(6,1),K_fdb(6,2),K_fdb(6,3)
      
c     Journal Bearing Damping (Top)
c     =============================
      READ(2,*)
      READ(2,*) C_fdb(4,1),C_fdb(4,2),C_fdb(4,3)
      READ(2,*) C_fdb(5,1),C_fdb(5,2),C_fdb(5,3)
      READ(2,*) C_fdb(6,1),C_fdb(6,2),C_fdb(6,3)
      
c     Thurst Bearing Stiffness
c     ========================
      READ(2,*)
      READ(2,*) K_fdb(7,1),K_fdb(7,2),K_fdb(7,3)
      READ(2,*) K_fdb(8,1),K_fdb(8,2),K_fdb(8,3)
      READ(2,*) K_fdb(9,1),K_fdb(9,2),K_fdb(9,3)
      
c     Thurst Bearing Damping
c     ======================
      READ(2,*)
      READ(2,*) C_fdb(7,1),C_fdb(7,2),C_fdb(7,3)
      READ(2,*) C_fdb(8,1),C_fdb(8,2),C_fdb(8,3)
      READ(2,*) C_fdb(9,1),C_fdb(9,2),C_fdb(9,3)
      
c******************************Rotor Hub*******************************

      READ(2,*)  

c     Disk interface data
c     ===================
	READ(2,*)
	READ(2,*) no_disk_intf_dofs
	
	ALLOCATE(disk_intf_dofs(no_disk_intf_dofs,2))
	ALLOCATE(disk_ID_map(no_disk_intf_dofs,5))
	ALLOCATE(disk_ID_disp(no_disk_intf_dofs,3))
	
      disk_intf_dofs = 0
      disk_ID_map    = 0
      disk_ID_disp   = 0.0d0
      
      READ(2,*)
      DO i = 1, no_disk_intf_dofs
        READ(2,*) temp_val,disk_ID_map(i,1),
     &      disk_ID_map(i,2),disk_intf_dofs(i,1)
      ENDDO
      
c     FDB nodes on bearing
c     ====================
      READ(2,*) 
      READ(2,*) no_crcm_node_hub
      
      no_fdb_node_hub = 3*no_crcm_node_hub  ! 2 Journal and 1 Thrust Bearing
      no_fdb_dofs_hub = 3*no_fdb_node_hub   ! 3 dofs for each node
      
      ALLOCATE(crcm_ang_hub(no_crcm_node_hub))  ! Circumferential angular position
      crcm_ang_hub = 0.0d0

      ALLOCATE(hub_fdb_dofs(no_fdb_dofs_hub))   ! FDB dofs
      hub_fdb_dofs = 0
      
      READ(2,*)
      DO i = 1, no_crcm_node_hub
        READ(2,*) (i_temp(j),j=1,10)
        
        crcm_ang_hub(i) = DBLE(i_temp(1))*(pi/180.0d0)              ! Angular Position
        
        hub_fdb_dofs((3*(i-1))+1) = i_temp(2)                       ! Lower Journal : Ux
        hub_fdb_dofs((3*(i-1))+2) = i_temp(3)                       ! Lower Journal : Uy
        hub_fdb_dofs((3*(i-1))+3) = i_temp(4)                       ! Lower Journal : Uz
        
        hub_fdb_dofs(no_crcm_node_hub*3+(3*(i-1))+1) = i_temp(5)    ! Upper Journal : Ux
        hub_fdb_dofs(no_crcm_node_hub*3+(3*(i-1))+2) = i_temp(6)    ! Upper Journal : Uy
        hub_fdb_dofs(no_crcm_node_hub*3+(3*(i-1))+3) = i_temp(7)    ! Upper Journal : Uz
        
        hub_fdb_dofs(no_crcm_node_hub*6+(3*(i-1))+1) = i_temp(8)    ! Thrust : Ux
        hub_fdb_dofs(no_crcm_node_hub*6+(3*(i-1))+2) = i_temp(9)    ! Thrust : Uy
        hub_fdb_dofs(no_crcm_node_hub*6+(3*(i-1))+3) = i_temp(10)   ! Thrust : Uz
        
      ENDDO
      
c*******************************Housing********************************
      
      READ(2,*)
      
c     FDB nodes on bearing
c     ====================
      READ(2,*) 
      READ(2,*) no_crcm_node_hsng
      
      no_fdb_node_hsng = 3*no_crcm_node_hsng        ! 2 Journal and 1 Thrust Bearing
      no_fdb_dofs_hsng = 3*no_fdb_node_hsng         ! 3 dofs for each node
      
      ALLOCATE(crcm_ang_hsng(no_crcm_node_hsng))    ! Circumferential angular position
      crcm_ang_hsng = 0.0d0      
      
      ALLOCATE(hsng_fdb_dofs(no_fdb_dofs_hsng))     ! FDB dofs
      hsng_fdb_dofs = 0

      READ(2,*)
      DO i = 1, no_crcm_node_hsng
        READ(2,*) (i_temp(j),j=1,10)
        
        crcm_ang_hsng(i) = DBLE(i_temp(1))*(pi/180.0d0)             ! Angular Position
        
        hsng_fdb_dofs((3*(i-1))+1) = i_temp(2)                      ! Lower Journal : Ux
        hsng_fdb_dofs((3*(i-1))+2) = i_temp(3)                      ! Lower Journal : Uy
        hsng_fdb_dofs((3*(i-1))+3) = i_temp(4)                      ! Lower Journal : Uz
        
        hsng_fdb_dofs(no_crcm_node_hsng*3+(3*(i-1))+1) = i_temp(5)  ! Upper Journal : Ux
        hsng_fdb_dofs(no_crcm_node_hsng*3+(3*(i-1))+2) = i_temp(6)  ! Upper Journal : Uy
        hsng_fdb_dofs(no_crcm_node_hsng*3+(3*(i-1))+3) = i_temp(7)  ! Upper Journal : Uz
        
        hsng_fdb_dofs(no_crcm_node_hsng*6+(3*(i-1))+1) = i_temp(8)  ! Thrust : Ux
        hsng_fdb_dofs(no_crcm_node_hsng*6+(3*(i-1))+2) = i_temp(9)  ! Thrust : Uy
        hsng_fdb_dofs(no_crcm_node_hsng*6+(3*(i-1))+3) = i_temp(10) ! Thrust : Uz
        
      ENDDO

c******************************Base Plate******************************
      
      READ(2,*)
      
      READ(2,*) 
      READ(2,*) no_base_dofs,base_ang
      base_ang = (pi/180.0d0)*base_ang  ! Convert to radians
      
c     Base - Housing Interface Nodes
c     ==============================
      READ(2,*) 
      READ(2,*) no_base_intf_node
      
      no_base_intf_dofs = 3*no_base_intf_node           ! Interface DOFs
      
c     Intialize
c     =========
      ALLOCATE(base_intf_dofs(no_base_intf_dofs,2))
      ALLOCATE(crcm_ang_base(no_base_intf_node))
      
      base_intf_dofs = 0
      crcm_ang_base  = 0.0d0
      
      READ(2,*)
      DO i = 1, no_base_intf_node
        
        i_temp = 0      ! Reset
        READ(2,*) (i_temp(j),j=1,7)
        
        crcm_ang_base(i) = DBLE(i_temp(1))*(pi/180.0d0) ! Angular Position

c       Column 1 : Base Plate
        base_intf_dofs(((3*(i-1))+1),1) = i_temp(2)     ! Base Plate : Ux
        base_intf_dofs(((3*(i-1))+2),1) = i_temp(3)     ! Base Plate : Uy
        base_intf_dofs(((3*(i-1))+3),1) = i_temp(4)     ! Base Plate : Uz
        
c       Column 2 : Motor Housing
        base_intf_dofs(((3*(i-1))+1),2) = i_temp(5)     ! Housing : Ux
        base_intf_dofs(((3*(i-1))+2),2) = i_temp(6)     ! Housing : Uy
        base_intf_dofs(((3*(i-1))+3),2) = i_temp(7)     ! Housing : Uz
        
      ENDDO           

c***************************Comment Section****************************

	CLOSE(2)

c**********************************************************************
c********************Finished Reading : DynLULS.def********************
c**********************************************************************

c	Pressure Output Files
c	=====================
	IF (outp_opt .NE. 0) THEN
	    OPEN(31,ERR=999,FILE='press1.dat',STATUS='UNKNOWN')
	    OPEN(32,ERR=999,FILE='press2.dat',STATUS='UNKNOWN')
	    OPEN(33,ERR=999,FILE='press3.dat',STATUS='UNKNOWN')
	    OPEN(34,ERR=999,FILE='press4.dat',STATUS='UNKNOWN')
	ENDIF

c     Point of Interest data
c     ======================
	OPEN(41,ERR=999,FILE='hpoint.dat',STATUS='UNKNOWN')       ! Clearance
	OPEN(42,ERR=999,FILE='hfloor.dat',STATUS='UNKNOWN')       ! Floor heights
	REWIND(41)
	REWIND(42)

c     Flying Dynamics
c     ===============	
      OPEN(51,ERR=999,FILE='fhhist.dat',STATUS='UNKNOWN')       ! Flying attitude		
	OPEN(52,ERR=999,FILE='wload.dat',STATUS='UNKNOWN')        ! ABS force
	OPEN(53,ERR=999,FILE='cpressures.dat',STATUS='UNKNOWN')   ! Maximum contact pressure	
	OPEN(54,ERR=999,FILE='contact.dat',STATUS='UNKNOWN')      ! Asperity contact force
	OPEN(55,ERR=999,FILE='impact.dat',STATUS='UNKNOWN')       ! Impact force
	OPEN(56,ERR=999,FILE='imf.dat',STATUS='UNKNOWN')          ! Intermolecular force
	OPEN(57,ERR=999,FILE='electro.dat',STATUS='UNKNOWN')      ! Electrostatic force
	REWIND(51)
	REWIND(52)
	REWIND(53)
	REWIND(54)
	REWIND(55)
	REWIND(56)
	REWIND(57)

c	Suspension Data
c	===============
      OPEN(61,ERR=999,FILE='slider.dat',STATUS='UNKNOWN')       ! Slider dynamics
	OPEN(62,ERR=999,FILE='cestatii.dat',STATUS='UNKNOWN')     ! Contact elements
	OPEN(63,ERR=999,FILE='lultab.dat',STATUS='UNKNOWN')       ! LUL tab
	OPEN(64,ERR=999,FILE='beam.dat',STATUS='UNKNOWN')         ! Load beam/flexure
	OPEN(65,ERR=999,FILE='actuator.dat',STATUS='UNKNOWN')     ! Actuator data
	REWIND(61)
	REWIND(62)
	REWIND(63)
	REWIND(64)
	REWIND(65)

c     Disk Support System Data
c     ========================	
      IF (luls_opt.EQ.0) THEN
      
c       Disk displacement
	  OPEN(71,ERR=999,FILE='disk_atd.dat',STATUS='UNKNOWN') 
	  REWIND(71)

c       Full disk displacement	  
	  IF ((disk_opt.EQ.1) .AND. (data_opt_disk.EQ.1)) THEN
	      OPEN(72,FILE='w_disk.dat',STATUS='UNKNOWN')
   		    REWIND(72)
            w_disk = 0.0d0
        ENDIF
        
        IF ((supp_opt.NE.0) .AND. (data_opt_supp.EQ.1)) THEN
        
c           Disk Hub interface displacement   	
   	      OPEN(73,FILE='disk_intf.dat',STATUS='UNKNOWN')
   	      REWIND(73)   	
   	
c           FDB displacement
            OPEN(74,FILE='w_FDB.dat',STATUS='UNKNOWN')
   	      REWIND(74)
   	  
        ENDIF
        
      ENDIF
      
c     Time stepping
c     =============
      ALLOCATE(disk_atd(INT(f_dtfac),3))
      disk_atd = 0.0d0
      WRITE(*,'(A16,E9.2)') 'Time step(s) = ',dt
	WRITE(*,'(A31,I3)')   'Time step reduction factor  = ', NINT(f_dtfac)

c     Air parameters
c     ==============
c     --------------------------
c     amu = viscocity
c     p00 = Ambient air pressure
c     al = Mean Free Path (MFP)
c     --------------------------
      g0 = 9.80665d0
      Ma = 0.0289644d0
      Rstar = 8.31432d0
      R = 287.05d0
      stdDensity = 1.2928d0
      S = 110.4d0
      MolecularDiameter = 3.65d-10
      Na  = 6.02213d23
      B0 = 1.458d-6
      tempKel = tempCel + 273.15d0
      
      IF (ABS(halt) .GT. 1) THEN

c       Pressure (barometric formula with temperature lapse)    
        ttt = 288.15d0 - 0.0065d0 * halt
        p0  = 101325.d0 * (ttt/288.15d0)**5.25588d0
        
c       p0 = 101325.0d0 * dexp((-g0 * Ma * calt) / (Rstar * tempKelvins));
        
c       Mean free path (from US Standard atmosphere)
        al = (Rstar * tempKel) / (dsqrt(2.d0) * pi * 
     &        MolecularDiameter*MolecularDiameter * Na * p0);
     
c       Sutherland formula for Viscosity
        amu = (B0 * tempKel**1.5) / (tempKel + 110.4)

c     Zero altitude
c     =============
      ELSE
        p0 = p00
        al = al0
        amu = amu0
      ENDIF

c     Disk Geometry
c     =============
	IF ((disk_opt.EQ.1) .AND. (luls_opt.EQ.0)) THEN

c		Disk Input Parameters
c		=====================
		disk_RPM = rpm	! disk rotates clockwise
        sldr_rad  = ra
        skw      = skw*(pi/180.0d0)       ! Degrees to radian
	  tht_ramp = tht_ramp*(pi/180.0d0)  ! Degrees to radian		
		
c       Triangle : Slider, VCM center and Disk Center
c       =============================================

c       Center to center distance
	  x_cent   = SQRT((ra**2)+(xact**2)-(2*ra*xact*DSIN(skw)))
c       Angle substended at VCM center	
	  tht_sldr = DACOS((x_cent**2 + xact**2 - ra**2)/(2*x_cent*xact))
c       Angle substended at disk ceter	
	  gam_sldr = (0.5*pi) + skw - tht_sldr
c       Slider angle : HDD frame of refernce
        sldr_ang  =  gam_sldr + base_ang - pi

c       Triangle : LUL Ramp, VCM center and Disk Center
c       ===============================================

c       Angle substended at VCM center
        tht_ramp = tht_sldr + tht_ramp      
c       Ramp Radius
        ramp_rad = SQRT((x_ramp**2)+(x_cent**2)-
     &		    (2*x_ramp*x_cent*DCOS(tht_ramp)))
c       Angle substended at disk center     
        gam_ramp = DACOS((x_cent**2 + ramp_rad**2 - x_ramp**2)/
     &		    (2*x_cent*ramp_rad))
c       Ramp angle : HDD frame of refernce 
        ramp_ang  =  gam_ramp + base_ang - pi
	
c       User defined shock 
c       ==================
        IF (shk_mod .LT. 5) THEN
	     shock_val = 0.0d0
	  ELSE
	     no_shkpts = 0
	     eof_flag = 0
           OPEN(501,FILE='shock_data.dat',STATUS='OLD')
	     REWIND(501)
	     DO WHILE (eof_flag .EQ. 0)
    	        no_shkpts = no_shkpts + 1
	        READ(501,*,IOSTAT=eof_flag) shock_val(no_shkpts,1),
     &	         shock_val(no_shkpts,2)
           ENDDO
	     CLOSE(501)
	     shkpts = no_shkpts - 1
	  ENDIF
	ENDIF
	
c     Suspension Positioning
c     ======================
      IF ((disk_opt.EQ.1) .AND. (luls_opt.EQ.0)) THEN
        len_c2c = x_cent	
      ELSE
        len_c2c = 0.0d0	
      ENDIF
c     Initialize      
      sldr_hdd_ang = base_ang - tht_sldr 
      sldr_disk_ang = gam_sldr
	
c     Normalize shock direction vector
c     ================================
      temp_val = sqrt((shk_dir(1,1)**2) + (shk_dir(2,1)**2) + 
     &                (shk_dir(3,1)**2))
      DO i = 1, 3
        shk_dir(i,1) = shk_dir(i,1)/temp_val
      ENDDO	
 
c	Power Term in Convective Law
c	============================
	IF(idisc.LT.1.OR.idisc.GT.5) idisc=1
	epsilon=1.0d-9
	nitmax=100
      
c     Electrostatic Constants
c     =======================
	eps0 = 8.85d-12   ! Permittivity constant
	Ke = 1            ! dielectric constant of medium i.e. air   
	elecnormf = xl*xl ! un-normalization factor
	
	felecst = 0.0d0
	telecst = 0.0d0
	relecst = 0.0d0

      IF(nsp.NE.0) THEN
        IF(someg(1,1).NE.0.d0) THEN
            WRITE(*,*)'Slider start/stop should start at t=0'
		    STOP
		ELSEIF (rpm.NE.someg(2,1)) THEN
		    WRITE(*,*)'RPM does not match initial start/stop speed'
		    STOP
		ENDIF
      ENDIF
      	  
	IF(rpm.EQ.0.d0.AND.nsp.NE.0) THEN
		rpmmax=someg(2,1)
		DO kk=2,nsp
			IF(someg(2,kk).GT.rpmmax) rpmmax=someg(2,kk)
		ENDDO
		omega=twopi*rpmmax/60.d0
	ELSE IF(rpm.EQ.0.d0.AND.nsp.EQ.0) THEN
		omega=twopi*5400.d0/60.d0
	ELSE
		omega=twopi*rpm/60.d0
	ENDIF
		omega0=omega
	DO i=1,nsp
		someg(2,i)=twopi*someg(2,i)/60.d0
	ENDDO

c	Convert to radians
c	==================
	ske = twopi*ske/360.d0
	ske_init=ske
	ra_init=ra

	p0xl=p0*xl*xl/9.81d0
	si =12.d0*amu*omega0*xl*xl/(p0*hm*hm)

	al=al/hm				
		   
	IF (iqpo.EQ.1) al=1.25d0*al
	
	d0 = 0.8862269254527d0/al 	
	icoe = 0					   
	nter = 11 				  
	gama = 0.57721566490153286d0
	pir  = 1.77245385090551603d0
	pit  = 1.023326708d0
	
	IF (iqpo.EQ.1) THEN	   
		t1 = 6.d0*al 
	ELSE IF (iqpo.EQ.2) THEN
		t1 = 6.d0*al
		t2 = 6.d0*al*al
	ELSE IF(iqpo.EQ.3.OR.iqpo.EQ.5) THEN
		t1 = 6.81971863421d0*al
		t2 = 8.31137590768d0*al*al
		t3 = 0.d0
		t4 = -10.9086332458d0*al*al*al*al
	ENDIF
																			 
	IF(xg.LE.0.d0) xg = xf0

c	Slider attitude computation
c	===========================
c	hmin=1.d0
	hx0=-hp*xl/hm
	hy=hs*xl/hm
	zc=-hx0*(1.d0-xg)+hmin
	
c     Slider grid
c     ===========
      no_xpts = max0(nint((nx-2.0)/16.0)*16+2,34)
      no_ypts = max0(nint((ny-2.0)/16.0)*16+2,34)
      
      ALLOCATE(xref_temp(no_xpts))
      ALLOCATE(yref_temp(no_ypts))
      
      xref_temp = 0.0d0
      yref_temp = 0.0d0
      
c     Check grid points for the form of 16n+2 (For New Grid)
c     ==========================================================
	IF(ioldgrid.EQ.0) THEN
		nx=max0(nint((nx-2.0)/16.0)*16+2,34)
		ny=max0(nint((ny-2.0)/16.0)*16+2,34)
		
		igp1=0
		
c	If old grid avialable : Read It
c	==============================
	ELSEIF(ioldgrid.EQ.1) THEN
	
		OPEN(11,file='x.dat',status='old')
		READ(11,*)(xref_temp(i),i=1,nx)
	  CLOSE(11)

	  OPEN(12,file='y.dat',status='old')
		READ(12,*)(yref_temp(j),j=1,ny)
	  CLOSE(12)

	  igp1=0
	
c       Check grid points for the form of 16n+2 (For Existing Grid)
c       ==========================================================
   
	  IF(MOD(nx,16).EQ.1.OR.MOD(ny,16).EQ.1) THEN
		    igp1=1
		    xref_temp(nx+1)=xref_temp(nx)
		    xref_temp(nx)=(xref_temp(nx-1)+xref_temp(nx+1))/2.d0
		    yref_temp(ny+1)=yref_temp(ny)
		    yref_temp(ny)=(yref_temp(ny-1)+yref_temp(ny+1))/2.d0
	  ENDIF

	  IF(MOD(nx,16).EQ.1) nx=nx+1
	  IF(MOD(ny,16).EQ.1) ny=ny+1
	  
	ENDIF

c     Initialize grid
c     ===============	
	ALLOCATE(xref(nx))
	ALLOCATE(yref(ny))
	
	ALLOCATE(ixtr(nx))
      ALLOCATE(jytr(ny))
      
      ixtr  = 0
      jytr  = 0
 	
	xref(1:nx) = xref_temp(1:nx)
	yref(1:ny) = yref_temp(1:ny)
	
	DEALLOCATE(xref_temp)
	DEALLOCATE(yref_temp)
	
	IF(igp1.EQ.1) THEN
		DO i=1,nx-2
			ixtr(i)=i
		ENDDO
		ixtr(nx-1)=nx
		DO j=1,ny-2
			jytr(j)=j
		ENDDO
		jytr(ny-1)=ny
	ELSE
		DO i=1,nx
			ixtr(i)=i
		ENDDO
		DO j=1,ny
			jytr(j)=j
		ENDDO
	ENDIF
	
	nxm1 = nx-1
	nxm2 = nxm1-1
	nym1 = ny-1
	nym2 = ny-2	

c**********************************************************************
c************************Input from : rail.dat*************************
c**********************************************************************

      OPEN(01,ERR=999,FILE='rail.dat',STATUS='OLD')
	REWIND(01)
	READ(01,'(a70)')line
	IF(line(1:3).EQ.'CML') THEN
		CLOSE(01)
		CALL read_rail
		CALL norm_rail

	ELSE
	  WRITE(*,*)'Rail.dat version 3.0 and earlier not supported'
	  CLOSE(01)
        STOP
        	
	ENDIF

c**********************************************************************
c*********************Finished Reading : rail.dat**********************
c**********************************************************************
															 

	WRITE(*,'(A30,I4,A7,I4)')' Actual Grid Size:        nx = '
     &							,nx,'ny = ',ny
	WRITE(*,'(A31,F10.6)')' Suspension Load(grams)      =  ',f0*1000.d0
	WRITE(*,'(A31,F10.6)')' Initial Disk Velocity(m/s)  = ',
     &						rpm*twopi*ra/60.d0
	WRITE(*,'(A31,F10.6)')' Initial Skew Angle(degrees) = ',
     &						 -ske*360.d0/twopi
	WRITE(*,'(A31,F10.6)')' Initial Radial Position(mm) = ',
     &						 ra*1000.d0

	WRITE(*,*) ' '
	WRITE(*,*) 'DYNAMIC INPUTS'
	WRITE(*,*) '=============='

	
	IF(nf_asper.NE.0) THEN
		WRITE(*,'(a20)')'	=> Disk Asperities'
	ENDIF
	 
	IF(nf_wave.NE.0) THEN
		WRITE(*,'(a18)')'	=> Disk Waviness'
	ENDIF
	  
	IF(nf_zone.NE.0) THEN
		WRITE(*,'(a23)')'	=> Disk Zone Profiles'
	ENDIF
	  
	IF(inwave.NE.0) THEN
		WRITE(*,'(a38)')'	=> Point by Point Disk Track Profile'
	ENDIF
	  
	IF(msusp.NE.0) THEN
		WRITE(*,'(a39)')'	=> Integration of Suspension Dynamics'
	ENDIF
	
	IF(nsp.NE.0) THEN
		WRITE(*,'(a33)')'	=> Time-Dependent Disk Velocity'
	ENDIF
	 
	IF(nap.NE.0.AND.nap.NE.1) THEN
		WRITE(*,'(a26)')'  => Track Seeking Motion'
	ENDIF
	 
	IF(nap.EQ.1) THEN
		WRITE(*,'(a23)')'	=> Actuator Slam Stop'
	ENDIF
	  
	IF(ncz.NE.0) THEN
		WRITE(*,'(a22)')'   => Asperity Contact'
		IF(icmod.EQ.1) THEN
			WRITE(*,'(a17)')'	   ( GW model )'
		ELSEIF(icmod.EQ.2) THEN
			WRITE(*,'(a30)')'	   ( Elastic-plastic model )'
	  ENDIF
	ENDIF
	
	IF(kft.NE.0) THEN
		WRITE(*,'(a21)')'	=> Disk Flutter'
	ENDIF
	
	IF (luls_opt.EQ.1) THEN
		WRITE(*,*)"  => Loading Process"
	ELSEIF (luls_opt.EQ.2) THEN
	  WRITE(*,*)"  => Unloading Process"
	ENDIF

      IF ((disk_opt.EQ.1) .AND. (luls_opt.EQ.0)) THEN
		WRITE(*,*)"  => Disk modeling on"
	ENDIF

      IF ((shk_mod.EQ.1).OR.(shk_mod.EQ.2)) THEN
		WRITE(*,*)"  => Shock response"
      ELSE IF ((shk_mod.EQ.3).OR.(shk_mod.EQ.4)) THEN
		WRITE(*,*)"  => Vibration response"
	ELSE
	  WRITE(*,*)"  => User defined shock"
	ENDIF

1000	RETURN

999	WRITE(*,*) 'Trouble in opening files'
	STOP
	
	END
	
c=====================================================================
      SUBROUTINE allocate_slider_arrays
c=====================================================================

c     Shared Data
c     ===========
      USE sldr_grid
      USE sldr_arys
      USE rynl_data
      USE impc_data
            
c     Slider Arrays
c     =============            
      
      ALLOCATE(h(nx,ny))
      ALLOCATE(hnew(nx,ny))
      
      ALLOCATE(hxl(nx,ny))
      ALLOCATE(hxr(nx,ny))
      ALLOCATE(hyd(nx,ny))
      ALLOCATE(hyu(nx,ny))

      ALLOCATE(hxlnew(nx,ny))
      ALLOCATE(hxrnew(nx,ny))
      ALLOCATE(hydnew(nx,ny))
      ALLOCATE(hyunew(nx,ny))
      
      ALLOCATE(hxy(nx,ny))
      ALLOCATE(hsad(nx,ny))
      
      ALLOCATE(spa(nx+1,ny+1))
      ALLOCATE(hfloor(nx,ny))
      
      h = 0.d0
      hnew = 0.d0 

      hxl = 0.d0
      hxr = 0.d0
      hyd = 0.d0
      hyu = 0.d0
      
      hxlnew = 0.d0
      hxrnew = 0.d0
      hydnew = 0.d0
      hyunew = 0.d0
      
      hxy = 0.d0
      hsad = 0.d0
      
      spa = 0.d0           
      hfloor = 0.d0
            
c     Reynold Equation Arrays
c     =======================
      ALLOCATE(p(nx,ny))
      ALLOCATE(pold(nx,ny))
      
      ALLOCATE(zta(nx,ny))
      ALLOCATE(eta(nx,ny))

      ALLOCATE(bearx(nx,ny))
      ALLOCATE(beary(nx,ny))
      
      p     = 1.d0
      pold  = 1.d0
      
      zta   = 0.d0
      eta   = 0.d0
      
      bearx = 0.0d0
      beary = 0.0d0  
      
c     Asperity contact pressure
c     =========================
      ALLOCATE(pcontact(nx,ny)) 
      pcontact = 0.0d0

      RETURN
      
      END
      
c=====================================================================
      SUBROUTINE allocate_reynolds_arrays
c=====================================================================
      
c     Shared Data    
c     ===========
      USE sldr_grid
      USE rynl_arys
      
c     Array Sizes
c     ===========
      nxmax = nx
      nymax = ny
      
      nxmax1 = 2 + (nxmax -2)/2
      nymax1 = 2 + (nymax -2)/2
      
      nxmax2 = 2 + (nxmax1-2)/2
      nymax2 = 2 + (nymax1-2)/2
      
      nxmax3 = 2 + (nxmax2-2)/2
      nymax3 = 2 + (nymax2-2)/2
      
      nxmax4 = 2 + (nxmax3-2)/2
      nymax4 = 2 + (nymax3-2)/2
      

      ALLOCATE(aj(nxmax)) 
      ALLOCATE(bj(nxmax))
      ALLOCATE(cj(nxmax)) 
      ALLOCATE(dj(nxmax)) 
      
      aj = 0.0d0
      bj = 0.0d0
      cj = 0.0d0
      dj = 0.0d0
      
      ALLOCATE(aw (nxmax ,nymax )) 
      ALLOCATE(ap (nxmax ,nymax )) 
      ALLOCATE(ae (nxmax ,nymax )) 
      ALLOCATE(ar (nxmax ,nymax )) 
      ALLOCATE(as (nxmax ,nymax )) 
      ALLOCATE(an (nxmax ,nymax )) 
      
      aw  = 0.0d0
      ap  = 0.0d0
      ae  = 0.0d0
      ar  = 0.0d0
      as  = 0.0d0
      an  = 0.0d0
      
      ALLOCATE(aw1(nxmax1,nymax1))
      ALLOCATE(ap1(nxmax1,nymax1)) 
      ALLOCATE(ae1(nxmax1,nymax1)) 
      ALLOCATE(ar1(nxmax1,nymax1)) 
      ALLOCATE(as1(nxmax1,nymax1)) 
      ALLOCATE(an1(nxmax1,nymax1)) 
      ALLOCATE(dp1(nxmax1,nymax1)) 
      
      aw1 = 0.0d0
      ap1 = 0.0d0
      ae1 = 0.0d0
      ar1 = 0.0d0
      as1 = 0.0d0
      an1 = 0.0d0
      dp1 = 0.0d0
      
      
      ALLOCATE(aw2(nxmax2,nymax2))
      ALLOCATE(ap2(nxmax2,nymax2)) 
      ALLOCATE(ae2(nxmax2,nymax2)) 
      ALLOCATE(ar2(nxmax2,nymax2)) 
      ALLOCATE(as2(nxmax2,nymax2)) 
      ALLOCATE(an2(nxmax2,nymax2)) 
      ALLOCATE(dp2(nxmax2,nymax2)) 
      
      aw2 = 0.0d0
      ap2 = 0.0d0
      ae2 = 0.0d0
      ar2 = 0.0d0
      as2 = 0.0d0
      an2 = 0.0d0
      dp2 = 0.0d0
      
      ALLOCATE(aw3(nxmax3,nymax3))
      ALLOCATE(ap3(nxmax3,nymax3)) 
      ALLOCATE(ae3(nxmax3,nymax3)) 
      ALLOCATE(ar3(nxmax3,nymax3)) 
      ALLOCATE(as3(nxmax3,nymax3)) 
      ALLOCATE(an3(nxmax3,nymax3)) 
      ALLOCATE(dp3(nxmax3,nymax3)) 
      
      aw3 = 0.0d0
      ap3 = 0.0d0
      ae3 = 0.0d0
      ar3 = 0.0d0
      as3 = 0.0d0
      an3 = 0.0d0
      dp3 = 0.0d0
      
      ALLOCATE(aw4(nxmax4,nymax4))
      ALLOCATE(ap4(nxmax4,nymax4)) 
      ALLOCATE(ae4(nxmax4,nymax4)) 
      ALLOCATE(ar4(nxmax4,nymax4)) 
      ALLOCATE(as4(nxmax4,nymax4)) 
      ALLOCATE(an4(nxmax4,nymax4)) 
      ALLOCATE(dp4(nxmax4,nymax4)) 
      
      aw4 = 0.0d0
      ap4 = 0.0d0
      ae4 = 0.0d0
      ar4 = 0.0d0
      as4 = 0.0d0
      an4 = 0.0d0
      dp4 = 0.0d0      
      
      RETURN
      
      END

c=====================================================================