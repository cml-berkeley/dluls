c======================================================================
      MODULE siml_pram
c======================================================================

      INTEGER nitmax,idisc,iadapt,msusp,igp1
      DOUBLE PRECISION akmax,emax,epsilon,p0xl,twopi

      END MODULE
      
c======================================================================
      MODULE syst_cnfg
c======================================================================     

      INTEGER nsp
      DOUBLE PRECISION xl,yl,xg,yg,zg
      DOUBLE PRECISION ra,ra_init,ra_if,ra_of,ske_init
      DOUBLE PRECISION si,amu,p0,al,rebase
      DOUBLE PRECISION rpm,omega0,omega
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: someg
     
      END MODULE
      
c======================================================================
      MODULE actr_data
c======================================================================
           
      INTEGER iact,nap
      DOUBLE PRECISION xact,dact,dst,vact,vst,aact
      DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: hac
      
      END MODULE

c======================================================================
      MODULE rail_data
c======================================================================

      INTEGER num_rails     ! Total number of rails
      INTEGER num_walls     ! Total number of walls
      
      INTEGER max_rail_pts   ! Max # of pts in a rail
	INTEGER max_wall_pts   ! Max # of pts in a wall
      
      DOUBLE PRECISION crown,twist,camber,xt,ht
      
      INTEGER, DIMENSION(:),   ALLOCATABLE :: istep
      INTEGER, DIMENSION(:),   ALLOCATABLE :: npoints
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: indexw
      
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xpt1
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: xpt2
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ypt1
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ypt2
      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: hramp
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: cramp
       
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xrail
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: yrail
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xdiff
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ydiff
     
      END MODULE

c======================================================================
      MODULE wall_data
c======================================================================
     
      INTEGER         , DIMENSION(:),   ALLOCATABLE :: nwpoint
      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: wpoint
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: wrecess
      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xwalli
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ywalli
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xwallo
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ywallo
      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::flush_lower_xwallo
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::flush_lower_ywallo
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::flush_upper_xwallo
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE ::flush_upper_ywallo
      
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xw1
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: yw1
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xw2
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: yw2
      
      END MODULE

c======================================================================
      MODULE grid_ctrl
c======================================================================

      INTEGER, PARAMETER :: nrmax = 50
            
      INTEGER nsx,nsy
      INTEGER nxt(nrmax),nyt(nrmax)
      INTEGER isymmetry,ioldgrid,ipmax
      
      DOUBLE PRECISION difmax,decay
      DOUBLE PRECISION xnt(nrmax),ynt(nrmax)
      DOUBLE PRECISION dxr(nrmax),dyr(nrmax)
      
      END MODULE
      
c======================================================================
      MODULE sldr_grid
c======================================================================

      INTEGER nx,nxm1,nxm2,ny,nym1,nym2

c     Map for grid (Important in grid size of 16n+1)
c     ==============================================
      INTEGER,          DIMENSION (:), ALLOCATABLE :: ixtr
      INTEGER,          DIMENSION (:), ALLOCATABLE :: jytr
      
c     Slider grid points
c     ==================      
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: xref 
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: yref
      
      END MODULE

c======================================================================
      MODULE sldr_arys
c======================================================================

      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: h,hnew
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: hxl,hxr
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: hyd,hyu
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: hxy,hsad
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: hxlnew,hxrnew
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: hydnew,hyunew
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: spa,hfloor

      END MODULE
      
c======================================================================
      MODULE sldr_rghn
c======================================================================
            
      INTEGER ns_wave,ns_asper
      
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: stype
	DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: smag
	DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: sang
	DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: sxwlth
	DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: sywlth
	DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: sxdm
	DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: sydm
	
	DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: satype
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: samag
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: saang
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: sxloc
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: syloc
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: saxdm
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: saydm
            
      END MODULE

c======================================================================
      MODULE disk_rghn
c======================================================================

      INTEGER, PARAMETER :: nwtp4 = 50, natp4 = 50
      
      INTEGER nwtp4n,natp4n
      INTEGER nf_wave,nf_asper,nf_zone
      	
c	Wavy roughness
c     ================
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: ftype
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: fmag
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: fang
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: fxwlth
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: fywlth
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: fxdm
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: fydm
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: fras
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: frae
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: flloc
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: flhgt
      
c     Circumferential zone profile
c     ===============================      
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: frp1
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: frp2
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: frp3
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: fcoe0
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: fcoe1
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: fcoe2

c     Asperities
c     ==========   
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: fatype
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: famag
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: faang
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: fxloc
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: fyloc
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: faxdm
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: faydm
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: falloc
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: falhgt
      
      END MODULE

c======================================================================
      MODULE disk_fltr
c======================================================================

      INTEGER kft
      DOUBLE PRECISION ft_ts,ft_tf,ft_omg,ft_mag
      
      END MODULE

c======================================================================  
      MODULE trck_prfl
c======================================================================  
      
      INTEGER nfx,inwave
      DOUBLE PRECISION dfx,dfy,dfinit,xfl
      
      DOUBLE PRECISION, DIMENSION (:),  ALLOCATABLE :: xfref
      DOUBLE PRECISION, DIMENSION (:),  ALLOCATABLE :: hfnew

      END MODULE

c======================================================================  
      MODULE rynl_data
c======================================================================

      INTEGER icoe,iqpo,nter,nreyneq
      
      DOUBLE PRECISION t1,t2,t3,t4
      DOUBLE PRECISION gama,pir,pit,d0,ak
      DOUBLE PRECISION zdat(1001),qndat(1001)
  
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: p,pold
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: zta,eta
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: bearx,beary
      
      END MODULE
      
c======================================================================
      MODULE rynl_arys
c======================================================================
      
      INTEGER nxmax ,nymax
	INTEGER nxmax1,nymax1
	INTEGER nxmax2,nymax2
	INTEGER nxmax3,nymax3
	INTEGER nxmax4,nymax4
	
	DOUBLE PRECISION,DIMENSION (:),   ALLOCATABLE :: aj,bj,cj,dj
	DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: aw ,ap ,ae 
	DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: ar ,as ,an
	DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: aw1,ap1,ae1
	DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: ar1,as1,an1,dp1
	DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: aw2,ap2,ae2
	DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: ar2,as2,an2,dp2
      DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: aw3,ap3,ae3
	DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: ar3,as3,an3,dp3
      DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: aw4,ap4,ae4
	DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: ar4,as4,an4,dp4
      
      END MODULE
c======================================================================
      MODULE sldr_dynm
c======================================================================

      DOUBLE PRECISION f0,xf0,yf0,f,xf,yf,fp,fn,ABS_data(4)
      DOUBLE PRECISION fshear,tshear,rshear,zshear
      DOUBLE PRECISION zc,hx0,hy,hmin,hm,hs
      DOUBLE PRECISION ske,skeeft,zang,xdisp,ydisp
      DOUBLE PRECISION dt,tf,t,tprev,dtorig,tout1,tout2,tout3
      
      END MODULE 
      
c======================================================================
      MODULE sldr_stat
c======================================================================
      
      DOUBLE PRECISION ztab0
      DOUBLE PRECISION xslider0,yslider0,zslider0
      DOUBLE PRECISION pslider0,rslider0,yawslider0
      
      END MODULE

c======================================================================
      MODULE aspr_data
c======================================================================

      INTEGER icmod,ncz
      DOUBLE PRECISION fcr,txr,tyr,tzr,f_maxcrp,aratio
      DOUBLE PRECISION rsikm(5),ceta(5),rasper(5)
      DOUBLE PRECISION rcts(5),rcte(5),gldht(5)
      
      END MODULE

c======================================================================
      MODULE impc_data
c======================================================================

      DOUBLE PRECISION fcont,fsctx,fscty,f_maximp,aratio1
      DOUBLE PRECISION eyoung,pratio,ydcoe,ydst,frcoe
      DOUBLE PRECISION eycrash,f_nu,ydcrash
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: pcontact

      END MODULE
      
c======================================================================
      MODULE vdwf_data
c======================================================================

      INTEGER imf_opt
      DOUBLE PRECISION ahc,bhc
      DOUBLE PRECISION fvdw,xfvdw,yfvdw

      END MODULE
      
c======================================================================
      MODULE elec_stat
c======================================================================

      DOUBLE PRECISION elecpot
      DOUBLE PRECISION eps0,Ke,elecnormf
      DOUBLE PRECISION felecst,telecst,relecst
      
      END MODULE

c======================================================================
      MODULE POI_cords
c======================================================================

      DOUBLE PRECISION xint(4),yint(4),thh(4)
      
      END MODULE

c======================================================================
      MODULE luls_data
c======================================================================      
      
      INTEGER luls_opt,outp_opt
      DOUBLE PRECISION cz00,cp00,cr00
      DOUBLE PRECISION toutp(4),foutp(4)
          
      END MODULE
c======================================================================
      MODULE shck_data
c======================================================================

c     User Defined Shock
c     ==================      
      INTEGER, PARAMETER :: max_shkpts = 200000 
      
c     Disk Support options
c     ====================      
      INTEGER supp_opt,data_opt_supp

c     Shock Parameters
c     ================
      INTEGER shk_mod,shkpts,shkpvt
      DOUBLE PRECISION T_st,T_sh,T_cnst,acc_amp,shk_dir(3,1)
      DOUBLE PRECISION acc_val,shock_val(max_shkpts,2)

      END MODULE
      
c======================================================================
      MODULE cmpt_dfmn
c======================================================================

c----------------------------------------------------------------------
c     Component deformation at slider position
c----------------------------------------------------------------------

c     Suspension
c     ==========
      DOUBLE PRECISION zc_susp,hx0_susp,hy_susp
      
c     Disk
c     ====      
      DOUBLE PRECISION w_trn,Dw_rad,Dw_ang,Dw_x_susp,Dw_y_susp
            
      END MODULE
            
c======================================================================
      MODULE susp_data
c======================================================================

c----------------------------------------------------------------------
c     Suspension modeling data
c----------------------------------------------------------------------     

      INTEGER max_susp

c     Suspension Reduced Matrices
c     ===========================      
      INTEGER i_suspsize,no_suspdofs,no_cdofs,nrp
      INTEGER i_dofux,i_dofuy,i_dofuz,i_dofrotx,i_dofroty,i_dofrotz
      INTEGER i_dofflx,i_doffly,i_doflbx,i_doflby,i_rampstat,i_doftab
      INTEGER i_celeramp,i_nocele
      INTEGER i_constat(10),i_dofcd(10),i_dofcu(10),i_cele(10)
      DOUBLE PRECISION f_fac,f_dtfac,f_rmass
      DOUBLE PRECISION f_cedist(10),f_cepl(10),f_cmassu(10),f_cmassd(10)
      DOUBLE PRECISION f_ulul,f_ululp1,f_vlul,f_vlulp1,f_alul,f_alulp1
      DOUBLE PRECISION fldz
      
c     HDD Positioning
c     ===============
      DOUBLE PRECISION len_c2c,sldr_hdd_ang,sldr_vcm_ang,sldr_disk_ang
      
c     Numerical Parameters
c     ====================
      DOUBLE PRECISION susp_nb_1,susp_nb_2      ! Newmark Beta
      DOUBLE PRECISION susp_M_damp,susp_K_damp  ! Structural Damping
      
      INTEGER,          DIMENSION (:),   ALLOCATABLE :: i_doftype
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: f_diagmass
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: f_un,f_unp1
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: f_vn,f_vnp1
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: f_an,f_anp1
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: f_rloc
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: f_coords
      
c     Suspension structural dynamics matrices
c     =======================================       
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: susp_M
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: susp_K
      
c     LU Factorization of Newmark Beta matrices
c     =========================================       
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: susp_Lc
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: susp_Lf
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: susp_Uc
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: susp_Uf
      
c     Matrices used for block factorization
c     =====================================
	DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: susp_A12
	DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: susp_A21
	DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: susp_A22
	DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: susp_A22_m

      DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: A_mat
	DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: B_mat
	DOUBLE PRECISION,DIMENSION (:,:), ALLOCATABLE :: C_mat
	      
      END MODULE    

c======================================================================
      MODULE disk_data
c======================================================================

c     Disk Option
c     ============
      INTEGER disk_opt
     
c     Disk Gemetry
c     ============
      DOUBLE PRECISION ri_disk,ro_disk,t_disk,disk_RPM
      
c     Disk Material Properties
c     ========================      
      DOUBLE PRECISION E_disk,nu_disk,den_disk
      
c     Disk Numerical Parameters
c     =========================
      DOUBLE PRECISION disk_nb_1,disk_nb_2      ! Newmark Beta
      DOUBLE PRECISION disk_M_damp,disk_K_damp  ! Structural Damping
      
c     Disk Grid
c     =========
      INTEGER no_disk_elem,no_disk_node,no_disk_dofs
      INTEGER no_rad_node,no_ang_node,no_rad_elem,no_ang_elem
      DOUBLE PRECISION del_rad,del_ang
      
c     Disk Matrix Size
c     ================      
      INTEGER no_fxddof,N_clamp
      
c     Ramp Contact Modeling
c     =======================
      INTEGER ramp_opt,i_rampcnt
      DOUBLE PRECISION ramp_clrnc,ramp_cntf
      
c     HDD Postioning
c     ==============
      DOUBLE PRECISION sldr_rad,sldr_ang
      DOUBLE PRECISION ramp_rad,ramp_ang   
      
c     Output Data Options
c     ===================
      INTEGER data_opt_disk,no_dpts 
      
c     Time Information
c     ================      
      DOUBLE PRECISION t_val      

      END MODULE

c======================================================================
      MODULE disk_aray
c======================================================================

c     Disk Grid
c     =========
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: rad_grid
	DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: ang_grid
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: disk_elem
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: disk_node

c     Disk Dynamic Matrices
c     =====================
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: disk_M
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: disk_G
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: disk_K
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: disk_map

c     Constraint BC
c     =============
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: C_disk
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: C_disk_map

c     LU Factorization Matrices
c     =========================
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_disk
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_disk
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: L_disk_map
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: U_disk_map
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: U_disk_T_map

c     Disk Dynamic Variables
c     ======================
	DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: w_disk_cur
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: v_disk_cur
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: a_disk_cur
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: w_disk_prv
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: v_disk_prv
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: a_disk_prv
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: f_disk_cnst
	
c     Disk deformation
c     ================
	DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: w_disk
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: disk_atd
	      
      END MODULE

c======================================================================
      MODULE disk_soln
c======================================================================   

c----------------------------------------------------------------------
c     Used for solve_disk only    
c----------------------------------------------------------------------
      
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: w_disk_prd
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: v_disk_prd
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: f_disk_val

	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: disk_vec_1
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: disk_vec_2
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: disk_vec_3
	
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: x_disk
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: y_disk
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: z_disk

c     Only for ramp constraint case
c     =============================	
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: ramp_coef_disk
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: C_disk_r
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: L_21_disk
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: U_12_disk
	
      END MODULE
      
c======================================================================
      MODULE motr_data
c======================================================================
       
c     DOF data
c     ========
      INTEGER no_motr_dofs,no_motr_vrbl
      INTEGER no_fdb_dofs_motr,no_fdb_dofs_op

c     Circumferential nodes
c     =====================
      INTEGER no_crcm_node_hub,no_crcm_node_hsng
      
c     FDB dynamic coefficients
c     ========================
      DOUBLE PRECISION K_fdb(9,3),C_fdb(9,3)
      
c     Disk Center Displacement
c     ========================
      DOUBLE PRECISION x0_disk_r,y0_disk_r
      DOUBLE PRECISION x0_disk_s,y0_disk_s
      
c     Inertia loading
c     ===============      
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: motr_inrt_load
      
c     Interface list
c     ==============
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: crcm_ang_hub
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: crcm_ang_hsng
      
c     Circumferential Interpolation
c     =============================
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: crcm_intp_mat
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: hub_cur_ang
      
c     Coupled rotor-housing matrices 
c     ===========================
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: motr_M
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: motr_G
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: motr_K
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: motr_M_map
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: motr_G_map
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: motr_K_map      

c     FDB Interface Matrices
c     ======================
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: intf_M
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: intf_G
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: intf_K

c     Disk Boundary conditions
c     ========================
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: C_motr
      INTEGER,          DIMENSION (:,:), ALLOCATABLE :: C_motr_map	

c     LU Factorization Matrices
c     =========================

c     Newmark Beta Block Matrices
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_motr_22_m
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_motr_23_m
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_motr_32_m
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_motr_33_m

c     Without ramp contact
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_motr
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_motr
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: L_motr_map
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: U_motr_map
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: U_motr_T_map
	
c     Motor Dynamic Variables
c     =========================
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: w_motr_cur
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: v_motr_cur
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: a_motr_cur
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: w_motr_prv
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: v_motr_prv
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: a_motr_prv
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: f_motr_cnst
      
c     Disk ID deformation
c     ===================
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: w_disk_intf
      
c     FDB deformation
c     ================
	DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: w_FDB

      END MODULE

c======================================================================
      MODULE motr_soln
c======================================================================

c----------------------------------------------------------------------
c     Used for solve_motor only    
c----------------------------------------------------------------------

c     Block matrices size
c     ===================
      INTEGER m_motr,n_motr,p_motr
      
c     Temporary Block Matrices
c     ========================
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: temp_G_m
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: temp_K_m

c     Newmark Beta Block Matrices
c     ===========================
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_motr_22

      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_motr_22
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_motr_32
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_motr_33
      
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_motr_22
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_motr_23
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_motr_33
      
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: blck_motr_1
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: blck_motr_2
      
c     Block matrices for disk ramp contact
c     ====================================
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: ramp_coef_motr
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: C_motr_r
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: L_21_motr
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: U_12_motr      
      
c     Newmark Beta Scheme vectors
c     ===========================
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: w_motr_prd
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: v_motr_prd
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: f_motr_val
	
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: motr_vec_1
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: motr_vec_2
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: motr_vec_3
	
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: x_motr
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: y_motr
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: z_motr

      
      END MODULE

c======================================================================
      MODULE rotr_data
c======================================================================

c     Node/Dofs Data
c     ==============
      INTEGER no_disk_bcs,no_disk_intf_dofs
      INTEGER no_fdb_node_hub,no_fdb_dofs_hub
      INTEGER no_hub_dofs,m_rotr
      INTEGER no_rotr_dofs,no_rotr_vrbl
      DOUBLE PRECISION disk_cent_disp(4)

c     Interface list
c     ==============
      INTEGER, DIMENSION (:),   ALLOCATABLE ::  hub_fdb_dofs
      INTEGER, DIMENSION (:,:), ALLOCATABLE ::  hub_intf_map
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: disk_ID_map
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: disk_intf_dofs
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: disk_intf_map

c     Cordinate matrices
c     ==================
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE ::  hub_cords
      
c     Disk ID deformation
c     ===================
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: disk_ID_disp
      
c     Inertia loading
c     ===============
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: hub_inrt_map
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: hub_inrt_acc
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: hub_inrt_vec
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: rotr_inrt_load
      
c     System structural matrices
c     ==========================      
      
c     Rotating Hub matrices      
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: hub_M
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: hub_G
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: hub_K
      
c     Coupled disk-rotor matrices 
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: rotr_M
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: rotr_G
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: rotr_K
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: rotr_M_map
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: rotr_G_map
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: rotr_K_map
	
c     Disk Boundary conditions
c     ========================
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: C_rotr
      INTEGER,          DIMENSION (:,:), ALLOCATABLE :: C_rotr_map

c     LU Factorization Matrices
c     =========================
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_rotr
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_rotr
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: L_rotr_map
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: U_rotr_map
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: U_rotr_T_map
	
c     Rotor Dynamic Variables
c     =======================
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: w_rotr_cur
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: v_rotr_cur
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: a_rotr_cur
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: w_rotr_prv
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: v_rotr_prv
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: a_rotr_prv
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: f_rotr_cnst

      END MODULE

c======================================================================
      MODULE rotr_soln
c======================================================================   

c----------------------------------------------------------------------
c     Used for solve_rotor only    
c----------------------------------------------------------------------

      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: w_rotr_prd
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: v_rotr_prd
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: f_rotr_val
	
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: rotr_vec_1
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: rotr_vec_2
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: rotr_vec_3
	
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: x_rotr
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: y_rotr
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: z_rotr
	
c     Only for ramp constraint case
c     =============================	
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: ramp_coef_rotr
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: C_rotr_r
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: L_21_rotr
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: U_12_rotr

      END MODULE
      
c======================================================================
      MODULE hsng_data
c======================================================================

c     Node/Dofs Data
c     ==============
      INTEGER no_hsng_dofs
      INTEGER no_fdb_node_hsng,no_fdb_dofs_hsng
      
c     Interface list
c     ==============
      INTEGER, DIMENSION (:),   ALLOCATABLE :: hsng_fdb_dofs
      
c     Structural matrices
c     ===================
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: hsng_M
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: hsng_G
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: hsng_K      
      
c     Cordinate matrices
c     ==================
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: hsng_cords
      
c     Inertia loading
c     ===============
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: hsng_inrt_acc
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: hsng_inrt_load
	
      END MODULE
      
c======================================================================
      MODULE base_data
c======================================================================
      
c     Base Geometry
c     =============
      DOUBLE PRECISION base_ang
      
c     DOF Data
c     ========
      INTEGER no_base_dofs,no_asmb_dofs
      INTEGER no_base_intf_node,no_base_intf_dofs
      
c     Housing - Base Interface
c     ========================
      INTEGER,          DIMENSION (:,:), ALLOCATABLE :: base_intf_dofs
      DOUBLE PRECISION, DIMENSION (:)  , ALLOCATABLE :: crcm_ang_base
      
c     Base plate structural matrices
c     ==============================

c     Base Plate
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: base_M
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: base_G
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: base_K

c     Base plate - motor housing assembly
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: asmb_M
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: asmb_G
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: asmb_K
	

c     Cordinate matrices
c     ==================
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: base_cords
      
c     Inertia loading
c     ===============
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: base_inrt_acc
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: asmb_inrt_load

      END MODULE

c======================================================================
      MODULE supp_data
c======================================================================

c     DOF Data
c     ========
      INTEGER no_supp_dofs,no_supp_vrbl,no_fdb_dofs_supp

c     Inertia loading
c     ===============      
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: supp_inrt_load
      
c     Disk Support system matrices 
c     ============================
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: supp_M
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: supp_G
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: supp_K
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: supp_M_map
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: supp_G_map
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: supp_K_map
	
c     Disk Boundary conditions
c     ========================
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: C_supp
      INTEGER,          DIMENSION (:,:), ALLOCATABLE :: C_supp_map
      
c     LU Factorization Matrices
c     =========================

c     Newmark Beta Block Matrices
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_supp_22_m
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_supp_23_m
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_supp_32_m
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_supp_33_m

c     Without ramp contact
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_supp
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_supp
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: L_supp_map
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: U_supp_map
	INTEGER,          DIMENSION (:,:), ALLOCATABLE :: U_supp_T_map      
      
      
c     Motor Dynamic Variables
c     =========================
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: w_supp_cur
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: v_supp_cur
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: a_supp_cur
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: w_supp_prv
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: v_supp_prv
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: a_supp_prv
      DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: f_supp_cnst

      END MODULE

c======================================================================
      MODULE supp_soln
c======================================================================

c----------------------------------------------------------------------
c     Used for solve_support only    
c----------------------------------------------------------------------

c     Block matrices size
c     ===================
      INTEGER m_supp,n_supp,p_supp
      
c     Temporary Block Matrices
c     ========================
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: temp_G_s
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: temp_K_s
      
c     Newmark Beta Block Matrices
c     ===========================
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_supp_22

      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_supp_22
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_supp_32
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_supp_33
      
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_supp_22
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_supp_23
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_supp_33
      
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: blck_supp_1
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: blck_supp_2
      
c     Block matrices for disk ramp contact
c     ====================================
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: ramp_coef_supp
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: C_supp_r
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: L_21_supp
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: U_12_supp      
      
c     Newmark Beta Scheme vectors
c     ===========================
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: w_supp_prd
      DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: v_supp_prd
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: f_supp_val
	
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: supp_vec_1
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: supp_vec_2
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: supp_vec_3
	
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: x_supp
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: y_supp
	DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: z_supp      
      
      END MODULE
c======================================================================