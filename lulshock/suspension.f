c**********************************************************************
c     Subroutines in this file :
c     1. init_susp
c     2. stat_susp
c     3. susp_0
c     4. solve_susp
c     5. advance_susp
c     6. susp_block_fact
c     7. insert_elem
c     8. remove_elem
c     9. apply_BC
c**********************************************************************

c======================================================================
      SUBROUTINE init_susp
c======================================================================

c----------------------------------------------------------------------
c     Initialize the suspension before the dynamic module
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE syst_cnfg
      USE actr_data
      USE sldr_arys
      USE trck_prfl
      USE sldr_dynm
      USE sldr_stat
      USE aspr_data
      USE impc_data
      USE luls_data
      USE susp_data     

	IMPLICIT REAL*8(a-h,o-z)

	max_susp = i_suspsize + i_nocele + 1
	
c     Allocate suspension matrices	
c     ============================
	ALLOCATE(i_doftype(max_susp))
      ALLOCATE(f_diagmass(max_susp))
      ALLOCATE(f_un(max_susp))
      ALLOCATE(f_vn(max_susp))
      ALLOCATE(f_an(max_susp))
      ALLOCATE(f_unp1(max_susp))
      ALLOCATE(f_vnp1(max_susp))
      ALLOCATE(f_anp1(max_susp))
      ALLOCATE(f_coords(max_susp,3))
      ALLOCATE(susp_M(max_susp,max_susp))
      ALLOCATE(susp_K(max_susp,max_susp))
      ALLOCATE(susp_Lc(max_susp,max_susp))
      ALLOCATE(susp_Uc(max_susp,max_susp))
      ALLOCATE(susp_Lf(max_susp,max_susp))
      ALLOCATE(susp_Uf(max_susp,max_susp))
      
c     Initialize Matrices
c     ===================
      i_doftype = 0
      f_diagmass = 0.0d0
      f_un = 0.0d0
      f_vn = 0.0d0
      f_an = 0.0d0
      f_unp1 = 0.0d0
      f_vnp1 = 0.0d0
      f_anp1 = 0.0d0
      f_coords = 0.0d0
      susp_M = 0.0d0
      susp_K = 0.0d0
      susp_Lc = 0.0d0
      susp_Uc = 0.0d0
      susp_Lf = 0.0d0
      susp_Uf = 0.0d0
     
      OPEN(90,ERR=999,FILE='susp_K.dat',STATUS='OLD')
      OPEN(91,ERR=999,FILE='susp_M.dat',STATUS='OLD')
	OPEN(92,ERR=999,FILE='susp_cords.dat',STATUS='OLD')

      DO i=1,i_suspsize
		READ(90,*) (susp_K(i,j),j=1,i_suspsize)
        READ(91,*) (susp_M(i,j),j=1,i_suspsize)
	  READ(92,*) (f_coords(i,j),j=1,3),i_doftype(i)
	ENDDO
      CLOSE(90)
      CLOSE(91)
	CLOSE(92)

c	Record the original size
c	========================
	no_suspdofs = i_suspsize
	no_cdofs = i_nocele + 1
     
c     The contact force scaling factor
c	================================ 
      f_fac = 10.d0/dt/dt

c     Tab DOF and the ramp contact element
c	====================================
	i_celeramp = i_suspsize + i_nocele + 1

c     Contact element definitions
c	===========================
	DO i=1,i_nocele
		i_cele(i) = i_suspsize + i
	ENDDO

c     Calculate the lumped mass parameters
c	====================================
	f_rmass = 0.d0
      DO k=1,i_nocele
		f_cmassu(k) = 0.d0
	  f_cmassd(k) = 0.d0
      ENDDO
      
	DO i=1,i_suspsize
	  f_diagmass(i) = 0.d0
	ENDDO

	DO i=1,i_suspsize
		f_rmass = f_rmass + susp_M(i_doftab,i)
	  DO k=1,i_nocele
	      f_cmassu(k) = f_cmassu(k) + susp_M(i_dofcu(k),i)
	      f_cmassd(k) = f_cmassd(k) + susp_M(i_dofcd(k),i)
	  ENDDO
	    
		DO j=1,i_suspsize
			f_diagmass(i) = f_diagmass(i) + susp_M(j,i)
	  ENDDO
	ENDDO

c     Calculate the contact element separations
c	=========================================
      DO k=1,i_nocele
		IF (f_cedist(k).LT.-10.d0)
     &	    f_cedist(k) = -f_coords(i_dofcu(k),i_doftype(i_dofcu(k)))
     &      +f_coords(i_dofcd(k),i_doftype(i_dofcd(k)))
	ENDDO

c	Get the static suspension displacements
c	=======================================
	CALL stat_susp

c     Assemble the global stffnesses with the contact elements
c	======================================================
      DO k = 1,i_nocele
		IF (i_constat(k).EQ.1) THEN
			susp_K(i_cele(k),i_dofcu(k)) = -1.d0/f_fac
            susp_K(i_cele(k),i_dofcd(k)) =  1.d0/f_fac
            susp_K(i_dofcu(k),i_cele(k)) = -1.d0/f_fac
            susp_K(i_dofcd(k),i_cele(k)) =  1.d0/f_fac
          ELSE
            susp_K(i_cele(k),i_cele(k)) = 1.d0/f_fac
          ENDIF
	ENDDO

      i_suspsize = i_suspsize+i_nocele

c	Initial state
c	=============
      CALL susp_0

c	The ramp contact element
c	========================
      i_suspsize = i_suspsize+1

c	DOF Value for next time step
c	============================

c	Unload
c	======
      IF (luls_opt.EQ.2) THEN
		i_rampstat = 0
        susp_K(i_celeramp,i_celeramp) = 1.d0/f_fac

        f_alulp1 = 0.d0
	  f_vlulp1 = 0.d0
	  f_ululp1 = ztab0 + f_rloc(2,1)

c	Load
c	====
      ELSEIF (luls_opt.EQ.1) THEN
		i_rampstat = 1

		susp_K(i_celeramp,i_celeramp) = 0.d0
        susp_K(i_celeramp,i_doftab) = -1.d0/f_fac
        susp_K(i_doftab,i_celeramp) = -1.d0/f_fac

      	f_alulp1 = 0.d0
	  f_vlulp1 = 0.d0
        f_ululp1 = f_un(i_doftab)

c	Shock
c	====
      ELSE
		i_rampstat = 0
        susp_K(i_celeramp,i_celeramp) = 1.d0/f_fac

        f_alulp1 = 0.d0
	  f_vlulp1 = 0.d0
	  f_ululp1 = 0.d0

      ENDIF	      

c	Update Values
c	=============
      f_alul = f_alulp1
	f_vlul = f_vlulp1
	f_ulul = f_ululp1

      hp = -hx0*hm/xl
	hyy = hy*hm/xl

	ALLOCATE(A_mat(no_suspdofs,no_suspdofs))

c	Initialize
c	==========
	
	A_mat = 0.0d0
      susp_Lc = 0.0d0
      susp_Uc = 0.0d0

c	Corase time step
c	================	
	DO i=1,no_suspdofs
		DO j=1,no_suspdofs
			A_mat(i,j)=susp_M(i,j)*(1.d0+susp_M_damp*
     &			(2.0d0*susp_nb_2)*dtorig)+susp_K(i,j)*
     &			(susp_nb_1*dtorig*dtorig+susp_K_damp*
     &          (2.0d0*susp_nb_2)*dtorig)
		ENDDO
	ENDDO
	CALL LU_fact(A_mat,
     &	susp_Lc(1:no_suspdofs,1:no_suspdofs),
     &	susp_Uc(1:no_suspdofs,1:no_suspdofs),
     &	no_suspdofs)

c	Reset
c	=====
	A_mat = 0.0d0
	susp_Lf = 0.0d0
      susp_Uf = 0.0d0

c	Fine time step
c	==============	
	DO i=1,no_suspdofs
		DO j=1,no_suspdofs
			A_Mat(i,j) = susp_M(i,j)*(1.d0+susp_M_damp*
     &			(2.0d0*susp_nb_2)*(dtorig/f_dtfac))+
     &			susp_K(i,j)*(susp_nb_1*(dtorig/f_dtfac)*
     &			(dtorig/f_dtfac)
     &			+susp_K_damp*(2.0d0*susp_nb_2)*(dtorig/f_dtfac))
		ENDDO
	ENDDO
	CALL LU_fact(A_mat,
     &	susp_Lf(1:no_suspdofs,1:no_suspdofs),
     &	susp_Uf(1:no_suspdofs,1:no_suspdofs),
     &	no_suspdofs)

      
c     Suspension Newmark Beta Block Factorization
c     ===========================================
	ALLOCATE(B_mat(no_suspdofs,no_cdofs))
	ALLOCATE(C_mat(no_cdofs,no_cdofs))

	ALLOCATE(susp_A12(no_suspdofs,no_cdofs))
	ALLOCATE(susp_A21(no_suspdofs,no_cdofs))
	ALLOCATE(susp_A22(no_cdofs,no_cdofs))
	ALLOCATE(susp_A22_m(no_cdofs,no_cdofs))
	      
      RETURN  

999	WRITE(*,*) 'Trouble in opening suspension files'
	STOP

	END

c======================================================================
	SUBROUTINE stat_susp
c======================================================================

c----------------------------------------------------------------------
c	Static deformation of suspension due to gram load on slider
c	Store this deformations as dof_0
c----------------------------------------------------------------------

c     Shared Data
c     ===========     
      USE syst_cnfg
      USE actr_data
      USE sldr_arys
      USE trck_prfl
      USE sldr_dynm
      USE sldr_stat
      USE aspr_data
      USE impc_data
      USE luls_data
      USE susp_data

	IMPLICIT REAL*8(a-h,o-z)
     
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: f_Kinv
      DOUBLE PRECISION, DIMENSION (:),  ALLOCATABLE :: f_u
      DOUBLE PRECISION, DIMENSION (:),  ALLOCATABLE :: f_F

      ALLOCATE(f_Kinv(max_susp,max_susp))
      ALLOCATE(f_u(max_susp))
      ALLOCATE(f_F(max_susp))
      
      f_Kinv = 0.0d0
      f_u    = 0.0d0
      f_F    = 0.0d0
		
c	Get the static deflections
c	==========================

c     Dimple
c     ======
      susp_K(i_cele(1),i_dofcu(1)) = -1.d0/f_fac
      susp_K(i_cele(1),i_dofcd(1)) =  1.d0/f_fac
      susp_K(i_dofcu(1),i_cele(1)) = -1.d0/f_fac
      susp_K(i_dofcd(1),i_cele(1)) =  1.d0/f_fac
      
c     Limiters
c     ========
      DO k = 2,i_nocele
		susp_K(i_cele(k),i_cele(k)) = 1.d0/f_fac
	ENDDO

      i_suspsize = i_suspsize + i_nocele

c	Gram load reaction on suspension
c	================================

	f_F(i_dofuz) = f0*9.81*1d3

      CALL inverse(susp_K,i_suspsize,max_susp,f_Kinv)      
	CALL mult(f_Kinv,f_F,i_suspsize,max_susp,1,1,f_u)

c	Slider Initial Conditions
c	=========================
      ztab0      = f_u(i_doftab)
	zslider0   = f_u(i_dofuz)
	pslider0   = f_u(i_dofroty)
      rslider0   = f_u(i_dofrotx)
      xslider0   = f_u(i_dofux)
      yslider0   = f_u(i_dofuy)
      yawslider0 = f_u(i_dofrotz)


      susp_K(i_cele(1),i_dofcu(1)) = 0.d0
      susp_K(i_cele(1),i_dofcd(1)) = 0.d0
      susp_K(i_dofcu(1),i_cele(1)) = 0.d0
      susp_K(i_dofcd(1),i_cele(1)) = 0.d0

      DO k = 2,i_nocele
 		susp_K(i_cele(k),i_cele(k)) = 0.d0
 	ENDDO

      i_suspsize = i_suspsize-i_nocele


      DEALLOCATE(f_Kinv)
	DEALLOCATE(f_u)

	END

c======================================================================
	SUBROUTINE susp_0
c======================================================================

c----------------------------------------------------------------------
c	The suspension initial deformation: Initial flying attitude + 
c	deformation due to gram load (suspension deformation after loading)
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE syst_cnfg
      USE actr_data
      USE sldr_arys
      USE trck_prfl
      USE sldr_dynm
      USE sldr_stat
      USE aspr_data
      USE impc_data
      USE luls_data    
      USE susp_data

	IMPLICIT REAL*8(a-h,o-z)

      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: f_F
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: f_u
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: f_coeff
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: f_temp3
      DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: f_temp4
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: f_Kinv
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: f_temp1
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: f_temp2

	ALLOCATE(f_u(max_susp))
      ALLOCATE(f_F(max_susp))
	ALLOCATE(f_coeff(max_susp))
      ALLOCATE(f_temp3(max_susp))
      ALLOCATE(f_temp4(max_susp))
      ALLOCATE(f_Kinv(max_susp,max_susp))
      ALLOCATE(f_temp1(max_susp,max_susp))
      ALLOCATE(f_temp2(max_susp,max_susp))
      
      f_u = 0.0d0
      f_F = 0.0d0
      f_Kinv = 0.0d0
      f_coeff = 0.0d0
      f_temp1 = 0.0d0
      f_temp2 = 0.0d0
      f_temp3 = 0.0d0
      f_temp4 = 0.0d0

c	Initialize the suspension (flying slider)
c	=========================================
c	Unload or shock : Flying Slider

      IF ((luls_opt.EQ.2).OR.(luls_opt.EQ.0)) THEN 
        f_un(i_dofux) = xdisp+xslider0
	  f_un(i_dofuy) = ydisp+yslider0
	  f_un(i_dofuz) = (zc*1d3*hm+zslider0)
	  f_un(i_dofrotx) = (-hy*hm/xl+rslider0-cr00)
	  f_un(i_dofroty) = (-hx0*hm/xl+pslider0-cp00)
	  f_un(i_dofrotz) = zang+yawslider0

        DO i=1,i_suspsize
		    f_F(i) = 0.d0
	  ENDDO

c       Contact element separations
c	  ===========================
		DO i=1,i_nocele
			IF (i_constat(i).EQ.1) THEN
			    f_F(i_cele(i)) = f_cedist(i)/f_fac
	      ENDIF
	      f_F(i_dofcu(i)) = f_F(i_dofcu(i))-f_cepl(i)
            f_F(i_dofcd(i)) = f_F(i_dofcd(i))+f_cepl(i)
	  ENDDO

c       Inertia forces
c	  ==============
		DO i=1,i_suspsize
			IF (i_doftype(i).EQ.1) THEN
				f_temp3(i) = ((vact)**2)*f_coords(i,1)
      	    ELSEIF(i_doftype(i).EQ.2) THEN
	            f_temp3(i) = (aact)*f_coords(i,1)
      	    ELSEIF(i_doftype(i).EQ.6) THEN
	            f_temp3(i) = (aact)
	        ELSE
	            f_temp3(i) = 0.d0
      	    ENDIF
		ENDDO

c	  f_rotation for Suspension
c	  =========================      
	  CALL mult(susp_M,f_temp3,i_suspsize,max_susp,1,1,f_temp4)

        DO i=1,i_suspsize
		    f_F(i) = f_F(i) + f_temp4(i)
	  ENDDO

c       Apply the initial displacements
c	  ===============================
        i_temp = i_suspsize

c       Introduced the prescribed u for slider dof
c	  ==========================================

        CALL remove_elem(f_F,i_temp,max_susp,i_dofux)
        CALL apply_BC(susp_K,f_temp1,i_temp,max_susp,
     &		f_coeff,i_dofux)
	  i_temp = i_temp-1
        DO i=1,i_temp
			f_F(i) = f_F(i) - f_coeff(i)*f_un(i_dofux)
	  ENDDO

		CALL remove_elem(f_F,i_temp,max_susp,i_dofux)
        CALL apply_BC(f_temp1,f_temp2,i_temp,max_susp,
     &		f_coeff,i_dofux)
		i_temp = i_temp-1
        DO i=1,i_temp
			f_F(i) = f_F(i) - f_coeff(i)*f_un(i_dofuy)
	  ENDDO

        CALL remove_elem(f_F,i_temp,max_susp,i_dofux)
		CALL apply_BC(f_temp2,f_temp1,i_temp,max_susp,
     &        f_coeff,i_dofux)
	  i_temp = i_temp-1
        DO i=1,i_temp
			f_F(i) = f_F(i) - f_coeff(i)*f_un(i_dofuz)
	  ENDDO

        CALL remove_elem(f_F,i_temp,max_susp,i_dofux)
		CALL apply_BC(f_temp1,f_temp2,i_temp,max_susp,
     &		f_coeff,i_dofux)
	  i_temp = i_temp-1
        DO i=1,i_temp
			f_F(i) = f_F(i) - f_coeff(i)*f_un(i_dofrotx)
	  ENDDO

        CALL remove_elem(f_F,i_temp,max_susp,i_dofux)
        CALL apply_BC(f_temp2,f_temp1,i_temp,max_susp,
     &        f_coeff,i_dofux)
	  i_temp = i_temp-1
        DO i=1,i_temp
	      f_F(i) = f_F(i) - f_coeff(i)*f_un(i_dofroty)
	  ENDDO

        CALL remove_elem(f_F,i_temp,max_susp,i_dofux)
        CALL apply_BC(f_temp1,f_temp2,i_temp,max_susp,
     &         f_coeff,i_dofux)
	  i_temp = i_temp-1
        DO i=1,i_temp
	      f_F(i) = f_F(i) - f_coeff(i)*f_un(i_dofrotz)
	  ENDDO

        CALL inverse(f_temp2,i_temp,max_susp,f_Kinv)
        CALL mult(f_Kinv,f_F,i_temp,max_susp,1,1,f_u)

c	  Reconstruct dof vector after solving u = K^-1*f
c	  ===============================================

		CALL insert_elem(f_u,i_temp,max_susp,i_dofux,f_un(i_dofux))
	  i_temp = i_temp+1
        CALL insert_elem(f_u,i_temp,max_susp,i_dofuy,f_un(i_dofuy))
	  i_temp = i_temp+1
        CALL insert_elem(f_u,i_temp,max_susp,i_dofuz,f_un(i_dofuz))
	  i_temp = i_temp+1
        CALL insert_elem(f_u,i_temp,max_susp,i_dofrotx,f_un(i_dofrotx))
	  i_temp = i_temp+1
        CALL insert_elem(f_u,i_temp,max_susp,i_dofroty,f_un(i_dofroty))
	  i_temp = i_temp+1
        CALL insert_elem(f_u,i_temp,max_susp,i_dofrotz,f_un(i_dofrotz))
		i_temp = i_temp+1

c       Non flying initial condition i.e. suspension is loaded on the ramp
c	  ==================================================================
	  ELSEIF (luls_opt.EQ.1) THEN

		i_rseg = 1	! Ramp segment

c	  Searching for ramp segment
c	  ==========================          
		DO WHILE ((-dact).GT.f_rloc(1,i_rseg+1)) 
			i_rseg = i_rseg+1
	  ENDDO

        IF (i_rseg.LT.nrp) THEN
		    f_ulul = ztab0+(f_rloc(2,i_rseg)
     &            *(f_rloc(1,i_rseg+1)+dact)
     &            +f_rloc(2,i_rseg+1)*(-dact-f_rloc(1,i_rseg)))
     &            /(f_rloc(1,i_rseg+1)-f_rloc(1,i_rseg))
        ELSE 
            WRITE(*,*) 'Tab is off the ramp initially for 
     &		    loading process'
	      STOP
        ENDIF

        f_un(i_doftab) = f_ulul

        DO i=1,i_suspsize
			f_F(i) = 0.d0
	  ENDDO

c       Contact element separations
c	  ===========================
        DO i=1,i_nocele
			IF (i_constat(i).EQ.1) THEN
	          f_F(i_cele(i)) = f_cedist(i)/f_fac
	      ENDIF
	      f_F(i_dofcu(i)) = f_F(i_dofcu(i))-f_cepl(i)
            f_F(i_dofcd(i)) = f_F(i_dofcd(i))+f_cepl(i)
	  ENDDO


c       Inertia forces
c	  ==============
        DO i=1,i_suspsize
			IF (i_doftype(i).EQ.1) THEN
	            f_temp3(i) = ((vact)**2)*f_coords(i,1)
      	    ELSEIF(i_doftype(i).EQ.2) THEN
	            f_temp3(i) = (aact)*f_coords(i,1)
      	    ELSEIF(i_doftype(i).EQ.6) THEN
	            f_temp3(i) = (aact)
	      ELSE
	            f_temp3(i) = 0.d0
      	    ENDIF
	  ENDDO
        CALL mult(susp_M,f_temp3,i_suspsize,max_susp,1,1,f_temp4)

        DO i=1,i_suspsize
			f_F(i) = f_F(i) + f_temp4(i)
	  ENDDO

c       Initial displcaments
c	  ====================      
        i_temp = i_suspsize

        CALL remove_elem(f_F,i_temp,max_susp,i_doftab)
        CALL apply_BC(susp_K,f_temp1,i_temp,max_susp,
     &         f_coeff,i_doftab)
	  i_temp = i_temp-1
        DO i=1,i_temp
			f_F(i) = f_F(i) - f_coeff(i)*f_un(i_doftab)
	  ENDDO

        CALL inverse(f_temp1,i_temp,max_susp,f_Kinv)
        CALL mult(f_Kinv,f_F,i_temp,max_susp,1,1,f_u)

        CALL insert_elem(f_u,i_temp,max_susp,i_doftab,f_un(i_doftab))
		i_temp = i_temp+1


	ENDIF

	zc  =  (f_un(i_dofuz)-zslider0)*1d-3/hm
	hx0 = -(f_un(i_dofroty)-pslider0+cp00)*xl/hm
	hy  = -(f_un(i_dofrotx)-rslider0+cr00)*xl/hm

      DO i=1,i_suspsize
		f_un(i) = f_u(i)
	ENDDO

      DEALLOCATE(f_F)
      DEALLOCATE(f_Kinv)
      DEALLOCATE(f_u)
      DEALLOCATE(f_temp1)
      DEALLOCATE(f_temp2)
      DEALLOCATE(f_temp3)
      DEALLOCATE(f_temp4)
      DEALLOCATE(f_coeff)

      RETURN
         
      END   

c======================================================================
      SUBROUTINE solve_susp(i_dtmode)
c======================================================================

c----------------------------------------------------------------------
c	Computes deformation of suspension and disk at time t
c----------------------------------------------------------------------
      
c     Shared Data      
c     ===========
      USE syst_cnfg
      USE actr_data
      USE sldr_arys
      USE trck_prfl
      USE sldr_dynm
      USE sldr_stat
      USE aspr_data
      USE impc_data
      USE vdwf_data
      USE elec_stat
      USE luls_data
      USE shck_data
      USE susp_data
      USE cmpt_dfmn
      
	IMPLICIT REAL*8(a-h,o-z) 
	
      DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: f_F
      DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: f_temp1
      DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: f_temp2
      DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: f_temp3
      DOUBLE PRECISION,DIMENSION(:),  ALLOCATABLE :: f_temp4

	DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: susp_L
	DOUBLE PRECISION,DIMENSION(:,:),ALLOCATABLE :: susp_U

      ssk=dsin(ske)
	csk=dcos(ske)

c	Time step size 
c	==============
	dts = t-tprev
	
	IF (ABS(t-tprev-dtorig).LE.(dtorig/100.d0)) THEN
		i_dt = 1
	ELSE
		i_dt = 2
	ENDIF
	
	IF (ABS(dt-dtorig).LE.(dtorig/100.d0)) THEN
		i_dt = 1
		dts = dt
	ELSE
		i_dt = 2
		dts = dt
	ENDIF
		
	IF (t.EQ.0.d0) THEN
		dts = dtorig
		i_dt = 1
	ENDIF

c	Suspension loading on slider
c	============================
      sfor = (f+fcr+fcont+fvdw+felecst)*1d3 ! Net z      force on slider
      sptor=-((f*(xf-xg))*xl+txr+fsctx      ! Net pitch torque on slider
     &      +tshear+yfvdw+telecst)*1d6
      srtor=-((f*(yf-yg))*xl+tyr+fscty      ! Net roll  torque on slider
     &		-rshear-xfvdw+relecst)*1d6

c	Shear force
c	===========
      sxfor=(fshear+(fcr+fcont)*frcoe)*csk*1d3
      syfor=(fshear+(fcr+fcont)*frcoe)*ssk*1d3
      sztor=(zshear+tzr)*1d6

c     Shock magnitude (inertia loading : a)
c	=====================================
      IF (luls_opt .EQ. 0) THEN
        f_accmag = -acc_val*1.0d3
      ELSE ! Load or Unload Mode
        f_accmag = 0.0d0
      ENDIF
      
      ALLOCATE(f_F(max_susp))
      ALLOCATE(f_temp1(max_susp))
      ALLOCATE(f_temp2(max_susp))
      ALLOCATE(f_temp3(max_susp))
      ALLOCATE(f_temp4(max_susp))
      
      f_F = 0.0d0
      f_temp1 = 0.0d0
      f_temp2 = 0.0d0
      f_temp3 = 0.0d0
      f_temp4 = 0.0d0
	
c	L/UL Process
c	============
      IF ((luls_opt.EQ.1).OR.(luls_opt.EQ.2)) THEN
		i_rseg=1
		DO WHILE ((-dact).GT.f_rloc(1,i_rseg+1)) 
			i_rseg = i_rseg+1
	    ENDDO

         	IF (i_rseg.LT.nrp) THEN
              f_ululp1= ztab0+(f_rloc(2,i_rseg)
     &			*(f_rloc(1,i_rseg+1)+dact)
     &            +f_rloc(2,i_rseg+1)*(-dact-f_rloc(1,i_rseg)))
     &            /(f_rloc(1,i_rseg+1)-f_rloc(1,i_rseg))
	     ELSE
			f_ululp1 = f_ulul
	     ENDIF

           f_vlulp1 = (f_ululp1-f_ulul)/dts	! Velocity
           f_alulp1 = (f_vlulp1-f_vlul)/dts	! Accelaration

	ENDIF

      DO i=1,i_suspsize
	      f_temp1(i) = f_un(i) + f_vn(i)*dts +
     &         0.5d0*((1.d0-2.d0*susp_nb_1)*f_an(i))*dts*dts + 
     &         susp_K_damp*(f_vn(i) + (1.d0-(2.0d0*susp_nb_2))*
     &         f_an(i)*dts)
            f_temp2(i) = susp_M_damp*(f_vn(i) + 
     &         (1.d0-(2.0d0*susp_nb_2))*f_an(i)*dts)
	ENDDO

      CALL mult(susp_K,f_temp1,i_suspsize,max_susp,1,1,f_temp3)
      CALL mult(susp_M,f_temp2,i_suspsize,max_susp,1,1,f_temp4)

      DO i=1,i_suspsize
		f_F(i) = -f_temp3(i) - f_temp4(i)
	ENDDO

c     Contact element forces
c	======================
	DO i=1,i_nocele
		IF (i_constat(i).EQ.1) THEN
				f_F(i_cele(i)) = f_F(i_cele(i))+f_cedist(i)/f_fac
	        ENDIF
			f_F(i_dofcu(i)) = f_F(i_dofcu(i))-f_cepl(i)
			f_F(i_dofcd(i)) = f_F(i_dofcd(i))+f_cepl(i)
		ENDDO

c     Shock forces in actuator frame
c     ==============================
      c_phi = DCOS(sldr_hdd_ang)
      s_phi = DSIN(sldr_hdd_ang)
      f_accmag_x = f_accmag*((shk_dir(1,1)*c_phi)+
     &                       (shk_dir(2,1)*s_phi))
      f_accmag_y = f_accmag*((shk_dir(1,1)*s_phi)-
     &                       (shk_dir(2,1)*c_phi))
      f_accmag_z =-f_accmag*shk_dir(3,1) 

c     Shock and actuator rotation forces
c	==================================
      DO i = 1 , i_suspsize
		IF (i_doftype(i).EQ.1) THEN
			f_temp1(i) = (((vact)**2)*f_coords(i,1)) + f_accmag_x
		ELSEIF(i_doftype(i).EQ.2) THEN
			f_temp1(i) = ((aact)*f_coords(i,1)) + f_accmag_y
		ELSEIF(i_doftype(i).EQ.3) THEN
			f_temp1(i) = f_accmag_z
		ELSEIF(i_doftype(i).EQ.6) THEN
			f_temp1(i) = (aact)
		ELSE
			f_temp1(i) = 0.d0
		ENDIF
	ENDDO

      CALL mult(susp_M,f_temp1,i_suspsize,max_susp,1,1,f_temp2)

      DO i=1,i_suspsize
		f_F(i) = f_F(i) + f_temp2(i)
	ENDDO


c     Imposing the ramp displacement
c	==============================
      IF (i_rampstat.EQ.1) THEN
		f_F(i_celeramp) = f_F(i_celeramp) -
     &         f_ululp1/f_fac - susp_K_damp*f_vlulp1/f_fac
	ENDIF

c     Air bearing forces
c	==================
      f_F(i_dofux) = f_F(i_dofux) + sxfor
      f_F(i_dofuy) = f_F(i_dofuy) + syfor
      f_F(i_dofuz) = f_F(i_dofuz) + sfor
      f_F(i_dofrotx) = f_F(i_dofrotx) + srtor
      f_F(i_dofroty) = f_F(i_dofroty) + sptor
      f_F(i_dofrotz) = f_F(i_dofrotz) + sztor


	ALLOCATE(susp_L(i_suspsize,i_suspsize))
	ALLOCATE(susp_U(i_suspsize,i_suspsize))
	
	susp_L = 0.0d0
	susp_U = 0.0d0

	IF (i_dtmode .EQ. 1) THEN
		susp_L(1:i_suspsize,1:i_suspsize) = 
     &		susp_Lc(1:i_suspsize,1:i_suspsize)
		susp_U(1:i_suspsize,1:i_suspsize) = 
     &		susp_Uc(1:i_suspsize,1:i_suspsize)
	
	ELSE
		susp_L(1:i_suspsize,1:i_suspsize) = 
     &		susp_Lf(1:i_suspsize,1:i_suspsize)
		susp_U(1:i_suspsize,1:i_suspsize) = 
     &		susp_Uf(1:i_suspsize,1:i_suspsize)
	ENDIF

	CALL LU_solve(susp_L,susp_U,f_F(1:i_suspsize),	
     &	f_anp1(1:i_suspsize),i_suspsize)

      DO i=1,i_suspsize
            f_unp1(i) = f_un(i) + f_vn(i)*dts +
     &         1/2*((1-2*susp_nb_1)*f_an(i)+
     &         2*susp_nb_1*f_anp1(i))*dts*dts
            f_vnp1(i) = f_vn(i) + 
     &         ((1-(2.0d0*susp_nb_2))*f_an(i)+
     &         (2.0d0*susp_nb_2)*f_anp1(i))*dts
	ENDDO

	fldz = f_F(i_dofuz)						! Suspension Load
      zang = f_unp1(i_dofrotz)                  ! Yaw

c	Reynolds equation input
c	=======================
	zc_susp  = (f_unp1(i_dofuz)-zslider0)*1d-3/hm
	hx0_susp = -(f_unp1(i_dofroty)-pslider0+cp00)*xl/hm
	hy_susp  = -(f_unp1(i_dofrotx)-rslider0+cr00)*xl/hm
			
      DEALLOCATE(f_F)
      DEALLOCATE(f_temp1)
      DEALLOCATE(f_temp2)
      DEALLOCATE(f_temp3)
      DEALLOCATE(f_temp4)
	
	DEALLOCATE(susp_L)
	DEALLOCATE(susp_U)

      RETURN
	END

c======================================================================
      SUBROUTINE advance_susp
c======================================================================

c----------------------------------------------------------------------
c	Update values for next time step
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE syst_cnfg
      USE actr_data
      USE sldr_arys
      USE trck_prfl
      USE sldr_dynm
      USE sldr_stat
      USE aspr_data
      USE impc_data
      USE luls_data
      USE shck_data
      USE susp_data
      USE disk_data
      USE cmpt_dfmn
      
	IMPLICIT REAL*8(a-h,o-z)

      DOUBLE PRECISION  f_celefor(10)

      CHARACTER(len=60) :: forspec

c	Update stiffness matrix for LUL tab and ramp contact
c	====================================================

c	For LUL Mode : Check ramp contact status
c	========================================
      IF ((luls_opt.EQ.1).OR.(luls_opt.EQ.2))THEN 


		IF(i_rampstat.EQ.0)	THEN
			f_d = f_unp1(i_doftab)-f_ululp1
	        IF (f_d.LE.-1.d-6)	THEN

				f_vnp1(i_doftab)	= f_vlul
                  f_unp1(i_celeramp)	= f_unp1(i_celeramp)-
     &               f_rmass*(f_anp1(i_doftab)-f_alul)
	            f_anp1(i_doftab) = f_alul

	            i_rampstat = 1

	            susp_K(i_celeramp,i_celeramp) = 0.d0
	            susp_K(i_celeramp,i_doftab) = -1.d0/f_fac
	            susp_K(i_doftab,i_celeramp) = -1.d0/f_fac

			ENDIF

		ELSEIF(i_rampstat.EQ.1) THEN
			IF (f_unp1(i_celeramp).LT.-f_fac/1.d3) THEN

				f_vnp1(i_doftab) = f_vnp1(i_doftab) + 
     &				dt*f_un(i_celeramp)/(2.d0*f_rmass)/f_fac

      	        f_anp1(i_doftab) = f_anp1(i_doftab) - 
     &				f_unp1(i_celeramp)/f_rmass/f_fac

                  f_unp1(i_celeramp) = 0.d0

      	        i_rampstat = 0

				susp_K(i_celeramp,i_celeramp) = 1.d0/f_fac
	            susp_K(i_celeramp,i_doftab) = 0.d0
	            susp_K(i_doftab,i_celeramp) = 0.d0


			ENDIF
		ENDIF

	ENDIF


c	Update stiffness matrix for dimple and limiters contact
c	=======================================================

      DO k=1,i_nocele
		IF(i_constat(k).EQ.0) THEN

			f_d = f_unp1(i_dofcu(k))-f_unp1(i_dofcd(k))

			IF (f_d.LE.(-f_cedist(k)-1.d-8)) THEN
				f_unp1(i_cele(k)) = f_unp1(i_cele(k))-
     &				f_cmassd(k)*f_cmassu(k)/(f_cmassd(k)+
     &				f_cmassu(k))*(f_anp1(i_dofcu(k))-
     &				f_anp1(i_dofcd(k)))


      			i_constat(k) = 1

				susp_K(i_cele(k),i_cele(k))  =  0.d0/f_fac
				susp_K(i_cele(k),i_dofcu(k)) = -1.d0/f_fac
				susp_K(i_cele(k),i_dofcd(k)) =  1.d0/f_fac
				susp_K(i_dofcu(k),i_cele(k)) = -1.d0/f_fac
				susp_K(i_dofcd(k),i_cele(k)) =  1.d0/f_fac


			ENDIF

		ELSEIF(i_constat(k).EQ.1)THEN
			IF (f_unp1(i_cele(k)).LT.-f_fac/1.d3) THEN

				f_unp1(i_cele(k)) = 0.d0

     				i_constat(k) = 0

				susp_K(i_cele(k),i_cele(k)) = 1.d0/f_fac
				susp_K(i_cele(k),i_dofcu(k)) = 0.d0
				susp_K(i_cele(k),i_dofcd(k)) = 0.d0
				susp_K(i_dofcu(k),i_cele(k)) = 0.d0
				susp_K(i_dofcd(k),i_cele(k)) = 0.d0

			ENDIF

		ENDIF
	ENDDO

c	Update Previous step value
c	==========================

c	Suspension
c	==========
      DO i=1,i_suspsize
		f_un(i) = f_unp1(i)
	    f_vn(i) = f_vnp1(i)
	    f_an(i) = f_anp1(i)
	ENDDO

c	Slider
c	======
      f_alul = f_alulp1
      f_vlul = f_vlulp1
      f_ulul = f_ululp1

	IF (f_unp1(i_celeramp).LT.0) THEN
		f_rampfor = 0.d0
	ELSE
			f_rampfor = f_unp1(i_celeramp)/f_fac
	ENDIF
		
	DO i = 1,i_nocele
		IF (f_unp1(i_cele(i)).LT.0.d0) THEN
			f_celefor(i) = 0.d0
		ELSE
			f_celefor(i) = f_unp1(i_cele(i))/f_fac
		ENDIF
	ENDDO

c	Slider Dynamics
c	===============
      WRITE(61,201)t,fldz,
     &   f_unp1(i_dofux)-xslider0,f_unp1(i_dofuy)-yslider0,
     &   f_unp1(i_dofuz)-zslider0,f_unp1(i_dofrotx)-rslider0,
     &   f_unp1(i_dofroty)-pslider0,f_unp1(i_dofrotz)-yawslider0

c     Contact Element
c     ===============
      WRITE(forspec,'(A10,I1,A38)') 
     &   '(E15.9,3X,',i_nocele,'(I1,3X,E16.9,3X,E16.9,3X),2(E16.9,3X))'

      WRITE(62,forspec)t,(i_constat(j),
     &   f_unp1(i_dofcu(j))-f_unp1(i_dofcd(j))+f_cedist(j),
     &   f_celefor(j),
     &   j= 1,i_nocele)

c     LUL Tab
c     =======
      WRITE(63,202)t,i_rampstat,
     &   f_unp1(i_doftab)-ztab0,f_rampfor,
     &   dact,f_ulul-ztab0,f_vlul,f_alul

c     Load/Flexure Beam
c     =================
      WRITE(64,203)t,f_unp1(i_dofflx),f_unp1(i_doffly),
     &   f_unp1(i_doflbx),f_unp1(i_doflby)
     
201	FORMAT(E16.9,7(3x,E16.9)) 
202	FORMAT(E16.9,3x,I1,6(3x,E16.9))
203	FORMAT(E16.9,4(3x,E16.9)) 
      
      RETURN
	END

c======================================================================
      SUBROUTINE susp_block_fact(i_dtmode)
c======================================================================

c----------------------------------------------------------------------
c	Compute the LU factorization for suspension newmark beta matrices
c	i_dtmode = 1 : for corse time stepping
c	i_dtmode = 2 : for fine  time stepping
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE sldr_dynm
      USE susp_data

	IMPLICIT NONE 

	INTEGER i_dtmode
	DOUBLE PRECISION d_t
	
c	Coarse time step
c	=================
	IF (i_dtmode .EQ. 1) THEN
	
	  d_t = dtorig

c	  Newmark Beta sub matrices
c	  =========================
	  susp_A12(1:no_suspdofs,1:no_cdofs) = 
     &	    (susp_nb_1*d_t**2 + susp_K_damp*(2.0d0*susp_nb_2)*d_t)*
     &      susp_K(1:no_suspdofs,no_suspdofs+1:i_suspsize)

	  susp_A21(1:no_suspdofs,1:no_cdofs) = 
     &	    (susp_nb_1*d_t**2 + susp_K_damp*(2.0d0*susp_nb_2)*d_t)*
     &	    TRANSPOSE(susp_K(no_suspdofs+1:i_suspsize,1:no_suspdofs))

	  susp_A22(1:no_cdofs,1:no_cdofs) = 
     &	    (susp_nb_1*d_t**2 + susp_K_damp*(2.0d0*susp_nb_2)*d_t)*
     &	    susp_K(no_suspdofs+1:i_suspsize,no_suspdofs+1:i_suspsize)

c	  Block factorization
c	  ===================
	  CALL L_solve_mult(susp_Lc(1:no_suspdofs,1:no_suspdofs),
     &	    susp_A12(1:no_suspdofs,1:no_cdofs),
     &	    susp_Uc(1:no_suspdofs,no_suspdofs+1:i_suspsize),
     &	    no_cdofs,no_suspdofs)

	  A_mat(1:no_suspdofs,1:no_suspdofs) = 
     &	    TRANSPOSE(susp_Uc(1:no_suspdofs,1:no_suspdofs))

	  CALL L_solve_mult(A_mat(1:no_suspdofs,1:no_suspdofs),
     &	    susp_A21(1:no_suspdofs,1:no_cdofs),
     &	    B_mat(1:no_suspdofs,1:no_cdofs),
     &	    no_cdofs,no_suspdofs)

	  susp_Lc(no_suspdofs+1:i_suspsize,1:no_suspdofs) = 
     &	    TRANSPOSE(B_mat(1:no_suspdofs,1:no_cdofs))
     
        C_mat = MATMUL(susp_Lc(no_suspdofs+1:i_suspsize,1:no_suspdofs),
     &                 susp_Uc(1:no_suspdofs,no_suspdofs+1:i_suspsize))

	  susp_A22_m = susp_A22 - C_mat
	  
	  CALL LU_fact(susp_A22_m(1:no_cdofs,1:no_cdofs),
     &	    susp_Lc(no_suspdofs+1:i_suspsize,
     &	    no_suspdofs+1:i_suspsize),
     &	    susp_Uc(no_suspdofs+1:i_suspsize,
     &	    no_suspdofs+1:i_suspsize),no_cdofs)

c	Fine time step
c	==============	
	ELSE
	
	  d_t = (dtorig/f_dtfac)

c	  Newmark Beta sub matrices
c	  =========================
	  susp_A12(1:no_suspdofs,1:no_cdofs) = 
     &	    (susp_nb_1*d_t**2 + susp_K_damp*(2.0d0*susp_nb_2)*d_t)*
     &      susp_K(1:no_suspdofs,no_suspdofs+1:i_suspsize)

	  susp_A21(1:no_suspdofs,1:no_cdofs) = 
     &	    (susp_nb_1*d_t**2 + susp_K_damp*(2.0d0*susp_nb_2)*d_t)*
     &	    TRANSPOSE(susp_K(no_suspdofs+1:i_suspsize,1:no_suspdofs))

	  susp_A22(1:no_cdofs,1:no_cdofs) = 
     &	    (susp_nb_1*d_t**2 + susp_K_damp*(2.0d0*susp_nb_2)*d_t)*
     &	    susp_K(no_suspdofs+1:i_suspsize,no_suspdofs+1:i_suspsize)

c	  Block factorization
c	  ===================
	  CALL L_solve_mult(susp_Lf(1:no_suspdofs,1:no_suspdofs),
     &	    susp_A12(1:no_suspdofs,1:no_cdofs),
     &	    susp_Uf(1:no_suspdofs,no_suspdofs+1:i_suspsize),
     &	    no_cdofs,no_suspdofs)

	  A_mat(1:no_suspdofs,1:no_suspdofs) = 
     &	    TRANSPOSE(susp_Uf(1:no_suspdofs,1:no_suspdofs))

	  CALL L_solve_mult(A_mat(1:no_suspdofs,1:no_suspdofs),
     &	    susp_A21(1:no_suspdofs,1:no_cdofs),
     &	    B_mat(1:no_suspdofs,1:no_cdofs),
     &	    no_cdofs,no_suspdofs)

	  susp_Lf(no_suspdofs+1:i_suspsize,1:no_suspdofs) = 
     &	    TRANSPOSE(B_mat(1:no_suspdofs,1:no_cdofs))
        
        C_mat = MATMUL(susp_Lf(no_suspdofs+1:i_suspsize,1:no_suspdofs),
     &                 susp_Uf(1:no_suspdofs,no_suspdofs+1:i_suspsize))
     
	  susp_A22_m = susp_A22 - C_mat
	  
	  CALL LU_fact(susp_A22_m(1:no_cdofs,1:no_cdofs),
     &	    susp_Lf(no_suspdofs+1:i_suspsize,
     &	    no_suspdofs+1:i_suspsize),
     &	    susp_Uf(no_suspdofs+1:i_suspsize,
     &	    no_suspdofs+1:i_suspsize),no_cdofs)


	ENDIF
	      
      RETURN
	END

c======================================================================
      SUBROUTINE insert_elem(f_vec,N,NMAXSIZE,i_loc,f_val)
c======================================================================

	IMPLICIT NONE
	 
      INTEGER i,i_loc,N,NMAXSIZE
	DOUBLE PRECISION f_val,f_vec(NMAXSIZE)

	DO i=N,i_loc,-1
		f_vec(i+1) = f_vec(i)
	ENDDO

	f_vec(i_loc) = f_val
      
      RETURN
	END

c======================================================================
      SUBROUTINE remove_elem(f_vec,N,NMAXSIZE,i_loc)
c======================================================================
	
	IMPLICIT NONE 
	
      INTEGER i,i_loc,N,NMAXSIZE
	DOUBLE PRECISION  f_vec(NMAXSIZE)

	DO i=i_loc,N
		f_vec(i) = f_vec(i+1)
	ENDDO
      
      RETURN
	END

c======================================================================
      SUBROUTINE apply_BC(susp_M,f_Mres,N,NMAXSIZE,f_coeff,i_dof)
c======================================================================
	
	IMPLICIT NONE  

      INTEGER i,j,i_dof,N,NMAXSIZE
      DOUBLE PRECISION susp_M(NMAXSIZE,NMAXSIZE)
      DOUBLE PRECISION f_Mres(NMAXSIZE,NMAXSIZE)
      DOUBLE PRECISION f_coeff(NMAXSIZE)

	   DO i=1,N
	      IF (i.LT.i_dof) THEN
      	      f_coeff(i) = susp_M(i,i_dof)
	         DO j=1,N
	            IF (j.LT.i_dof) THEN
	               f_Mres(i,j) = susp_M(i,j)
	            ELSEIF (j.GT.i_dof) THEN
	               f_Mres(i,j-1) = susp_M(i,j)
	            ENDIF
	         ENDDO
	      ELSEIF (i.GT.i_dof) THEN
      	      f_coeff(i-1) = susp_M(i,i_dof)
	         DO j=1,N
	            IF (j.LT.i_dof) THEN
	               f_Mres(i-1,j) = susp_M(i,j)
	            ELSEIF (j.GT.i_dof) THEN
	               f_Mres(i-1,j-1) = susp_M(i,j)
	            ENDIF
	         ENDDO
	      ENDIF
	   ENDDO
      
      RETURN
	END

c======================================================================