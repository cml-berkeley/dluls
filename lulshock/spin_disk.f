c**********************************************************************
c     Subroutines in this file :
c     01. init_disk
c     02. gen_disk_mesh
c     03. gen_disk_matrices
c     04. newmark_disk
c     05. solve_disk
c     06. disk_attitude
c     07. rad_loc
c     08. ang_loc
c     09. ramp_cnstrnt
c     10. intrp_disk_atd
c**********************************************************************

c======================================================================
      SUBROUTINE init_disk
c======================================================================

c----------------------------------------------------------------------
c     Initialize disk
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE shck_data
      USE disk_data
      USE disk_aray
      USE motr_data
            
      IMPLICIT NONE
            
      DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979
      
c     Local Parameters
c     ================
      INTEGER i,j,info
	INTEGER row_cnt,dof_no,node_no
	
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: G_damp
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: G_net

      
c     Mesh Data
c     =========	
      no_rad_node  = no_rad_elem + 1
	no_ang_node  = no_ang_elem
	no_disk_elem = no_rad_elem*no_ang_elem
	no_disk_node = no_rad_node*no_ang_node
	no_disk_dofs = 4*no_disk_node
	no_fxddof    = 4*no_ang_node
 	N_clamp = no_disk_dofs + no_fxddof
 	del_rad = (DBLE(ro_disk) - DBLE(ri_disk))/DBLE(no_rad_elem)
 	del_ang = (2.0d0*pi)/DBLE(no_ang_elem)
 	
c	Mesh Generation
c	===============
	ALLOCATE(disk_node(no_disk_node,2))
	ALLOCATE(disk_elem(no_disk_elem,4))
	ALLOCATE( rad_grid(no_rad_node))
	ALLOCATE( ang_grid(no_ang_node))
	CALL gen_disk_mesh
	
c	System Matrices
c	===============
	ALLOCATE(disk_M(no_disk_dofs,no_disk_dofs))
	ALLOCATE(disk_G(no_disk_dofs,no_disk_dofs))
	ALLOCATE(disk_K(no_disk_dofs,no_disk_dofs))
	CALL gen_disk_matrices
	
c	Structural Damping 
c	==================
      ALLOCATE(G_damp(no_disk_dofs,no_disk_dofs))
	ALLOCATE( G_net(no_disk_dofs,no_disk_dofs))
	
	G_damp = (disk_M_damp*disk_M + disk_K_damp*disk_K)
	G_net  = G_damp + disk_G
	disk_G = G_net

	DEALLOCATE (G_damp)
	DEALLOCATE (G_net)
	
c     Constraint at r = ri (all dof fixed)
c	====================================
	ALLOCATE(C_disk(no_fxddof,no_disk_dofs))
	C_disk = 0.0d0
	
	row_cnt = 0
	DO i = 1 , no_ang_node
	  node_no = ((no_rad_node)*(i-1)) + 1
		DO j = 1 , 4
		    row_cnt = row_cnt + 1
			dof_no  = 4*(node_no - 1) + j
			C_disk(row_cnt,dof_no) = 1.0d0
		ENDDO
	ENDDO	

	ALLOCATE(C_disk_map(no_fxddof,no_disk_dofs+1))
	CALL sparse_mat_map(C_disk,C_disk_map,
     &	no_fxddof,no_disk_dofs)

      RETURN 
      END

c======================================================================
      SUBROUTINE gen_disk_mesh
c======================================================================

c----------------------------------------------------------------------
c     Generate annular disk mesh
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE disk_data
      USE disk_aray
      
      IMPLICIT NONE

c	Local Variables
c	================
	INTEGER i,j,k,node(4)
	INTEGER node_no,elem_no

c     Radial Divison
c	==============
 	DO i = 1,no_rad_node
		rad_grid(i) = ri_disk + (i-1)*del_rad
 	ENDDO

c	Circumfrential Divison
c	======================
 	DO i = 1,no_ang_node
 		ang_grid(i) = 0 + (i-1)*del_ang
 	ENDDO

c	Nodal Matrix
c	============
c	a) Polar Cordinates
c	b) Cartesian Cordinates
c     Note that last node is not included (Repeated in Circle)

	DO i = 1 , no_ang_node
		DO j = 1 , no_rad_node
			node_no = ((no_rad_node)*(i-1))+j
			disk_node(node_no,1) = rad_grid(j)
			disk_node(node_no,2) = ang_grid(i)
		ENDDO
	ENDDO

c	Element Connectivity
c	====================

c	All angular sectors except last
	DO i = 1 , no_ang_elem-1
		DO j = 1 , no_rad_elem
			node(1) = ((no_rad_node)*(i-1))+(j)
			node(2) = ((no_rad_node)*(i-1))+(j+1)
			node(3) = ((no_rad_node)*(i))+(j+1)
			node(4) = ((no_rad_node)*(i))+(j)
			elem_no = (no_rad_elem*(i-1))+j
			DO k = 1 , 4
				disk_elem(elem_no,k) = node(k)
			ENDDO
		ENDDO
	ENDDO

c	Last angular sector
	DO j = 1 , no_rad_elem
		node(1) = ((no_rad_node)*(no_ang_elem-1))+(j);
		node(2) = ((no_rad_node)*(no_ang_elem-1))+(j+1);
		node(3) = (j+1);
		node(4) = (j);
		elem_no = (no_rad_elem*(no_ang_elem-1))+j;
		DO k = 1 , 4
			disk_elem(elem_no,k) = node(k)
		ENDDO
	ENDDO

c     Output grid data
c     ================
      IF (data_opt_disk .EQ. 1) THEN
      
        OPEN(01,FILE='disk_node.dat',STATUS='UNKNOWN')
	  REWIND(01)
	  WRITE(01,'(E12.6E2,3X,E12.6E2)') 
     &	    ((disk_node(i,j), j=1,2), i=1,no_disk_node)
	  CLOSE(01)

	  OPEN(02,FILE='disk_elem.dat',STATUS='UNKNOWN')
	  REWIND(02)
	  WRITE(02,'(4(I6,3X))') 
     &      ((disk_elem(i,j), j=1,4), i=1,no_disk_elem)
	  CLOSE(02)
      
      ENDIF
      
      RETURN
      END
      
c======================================================================
	SUBROUTINE gen_disk_matrices
c======================================================================

c----------------------------------------------------------------------
c     Generate following disk structural matrices : 
c     - M : Intertia
c     - G : Gyroscopic
c     - K : Stiffness
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE disk_data
      USE disk_aray
      
      IMPLICIT NONE
      
      DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979
      
c	Local Variables
c	================
	INTEGER maplg(4),mapvec(16)
	INTEGER gauss_od,i_elem,i,j,m,n,i_vec,j_vec
	DOUBLE PRECISION c1,c2,c3,denm,qrr,qtt,D0,D_coef,Jacb,omg 
	DOUBLE PRECISION D(3,3),Q(2,2),B(3,16),F(2,16)
	DOUBLE PRECISION B_T(16,3),F_T(16,2),DN0_T(16,1)
	DOUBLE PRECISION M_elem(16,16),G_elem(16,16),K_elem(16,16)
	DOUBLE PRECISION re(4,1),te(4,1),r0,t0,rg,tg
	DOUBLE PRECISION DN0(16),DNr(16),DNt(16),
     &				   DNrr(16),DNtt(16),DNrt(16)
	DOUBLE PRECISION temp01(16,16),temp02(16,16),temp03(16,16),
     &				   temp04(16,16),temp05(16,16)
	DOUBLE PRECISION temp01_sum(16,16),temp02_sum(16,16),
     &				   temp03_sum(16,16),temp04_sum(16,16),
     &				   temp05_sum(16,16)
	DOUBLE PRECISION temp01_mult(3,16),temp02_mult(2,16)
	
 	DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: gauss_wt
	DOUBLE PRECISION, DIMENSION (:),   ALLOCATABLE :: gauss_pt
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: mat_pop
	
c	Angular Speed (rad/sec)
c	=======================	
	omg = disk_RPM*(2*pi/60)

c	Gauss Quadrature Order
c	======================
	gauss_od = 3

c	In plane stress coefficient(qr,qt)
c	==================================
	D0   = E_disk*t_disk/(1.0d0-(nu_disk**2))
	c3   = (den_disk*t_disk*(omg**2))/(8.0d0*D0)
	denm = ((ro_disk**2)*(1.0d0 + nu_disk)) - 
     &	   ((ri_disk**2)*(1.0d0 - nu_disk))	
	c1   = (c3/denm)*(((ro_disk**4)*(3.0d0 + nu_disk)) + 
     &	   ((ri_disk**4)*(1.0d0 - nu_disk)))
	c2   =-((c3*(ri_disk**2)*(ro_disk**2))/denm)*
     &	   (((ro_disk**2)*(3.0d0 + nu_disk)) - 
     &	   ((ri_disk**2)*(1.0d0 + nu_disk)))

	D_coef = (E_disk*t_disk**3)/(12.0d0*(1.0d0-nu_disk**2))

	D(1,1) = D_coef*1.0d0
	D(1,2) = D_coef*nu_disk
	D(1,3) = 0.0d0
	D(2,1) = D_coef*nu_disk
	D(2,2) = D_coef*1.0d0
	D(2,3) = 0.0d0
	D(3,1) = 0.0d0
	D(3,2) = 0.0d0
	D(3,3) = D_coef*0.5d0*(1.0d0-nu_disk)

 	ALLOCATE(gauss_wt(gauss_od))
	ALLOCATE(gauss_pt(gauss_od))
	
      gauss_wt = 0.0d0
      gauss_pt = 0.0d0	

  	CALL gauss_quad(gauss_od,gauss_wt,gauss_pt)
  	
c	Initialize global level matrices
c	================================
      disk_M = 0.0d0
      disk_G = 0.0d0
      disk_K = 0.0d0

c     Mark the locations that are populated
c     =====================================    
      ALLOCATE(mat_pop(no_disk_dofs,no_disk_dofs))
      ALLOCATE(disk_map(no_disk_dofs,no_disk_dofs+1))
      mat_pop = 0.0d0

c	System Matrices Computation
c	===========================
	DO i_elem = 1 , no_disk_elem

c	  Initialize element level matrices
c	  ================================
		DO i = 1 , 16
			DO j = 1 , 16
				M_elem(i,j) = 0
				G_elem(i,j) = 0
				K_elem(i,j) = 0
			ENDDO
		ENDDO

c	  Mapping vector : Relating local node number 
c	                 to global node number 
c	  ===========================================
 		DO i = 1 , 4
 			maplg(i) = disk_elem(i_elem,i)
 		ENDDO

c	  Global co-ordiantes for local nodes
c	  ====================================
		DO i = 1 , 4
			re(i,1) = disk_node(maplg(i),1)
			te(i,1) = disk_node(maplg(i),2)
		ENDDO

c       Centroid Polar Cordinates
c	  =========================
	  r0 = re(1,1) + 0.5d0 * del_rad
	  t0 = te(2,1) + 0.5d0 * del_ang

c	  Variables for Numerical Integration
c	  ===================================
	  temp01 = 0.0d0
	  temp02 = 0.0d0
	  temp03 = 0.0d0
	  temp04 = 0.0d0
	  temp05 = 0.0d0
	  temp01_sum = 0.0d0
	  temp02_sum = 0.0d0
	  temp03_sum = 0.0d0
	  temp04_sum = 0.0d0
	  temp05_sum = 0.0d0
	  temp01_mult = 0.0d0
	  temp02_mult = 0.0d0

c	  Numerical Integration : Gauss Quadrature
c	  ========================================
		DO m = 1 , gauss_od
			DO n = 1 , gauss_od

c	      Polar Cordinates
c	      ================
			rg = r0 + ((del_rad/2.0d0)*(gauss_pt(m)))
			tg = t0 + ((del_ang/2.0d0)*(gauss_pt(n)))

			CALL    N_annular(gauss_pt(m),gauss_pt(n),
     &	                      del_rad,del_ang,DN0)
			CALL  DNr_annular(gauss_pt(m),gauss_pt(n),
     &		                  del_rad,del_ang,DNr)
			CALL  DNt_annular(gauss_pt(m),gauss_pt(n),
     &		                  del_rad,del_ang,DNt)	
			CALL DNrr_annular(gauss_pt(m),gauss_pt(n),
     &		                  del_rad,del_ang,DNrr)
     			CALL DNtt_annular(gauss_pt(m),gauss_pt(n),
     &		                  del_rad,del_ang,DNtt)
			CALL DNrt_annular(gauss_pt(m),gauss_pt(n),
     &		                  del_rad,del_ang,DNrt)		

			Jacb = (rg*del_rad*del_ang/4.0d0)

	      qrr = D0*((c1*(1.0d0+nu_disk))-((c2/(rg**2))*
     &			     (1.0d0-nu_disk))-((c3*(rg**2))*(3.0d0+nu_disk)))
            qtt = D0*((c1*(1.0d0+nu_disk))+((c2/(rg**2))*
     &				 (1.0d0-nu_disk))-((c3*(rg**2))*
     &                 (1.0d0+(3.0d0*nu_disk))))

			DO i = 1 , 16
				B(1,i) = (-4.0d0/(del_rad**2))*DNrr(i)
				B(2,i) = -((2.0d0/(del_rad*rg))*DNr(i)) - 
     &                      ((4.0d0/((del_ang**2)*(rg**2)))*DNtt(i))
				B(3,i) = ((-8.0d0/(rg*del_rad*del_ang))*DNrt(i)) + 
     &					 ((4.0d0/((rg**2)*(del_ang)))*DNt(i))

	          F(1,i) = (2.0d0/del_rad)*DNr(i)
				F(2,i) = (2.0d0/(del_ang*rg))*(DNt(i))
			ENDDO
            
            B_T = TRANSPOSE(B)
            F_T = TRANSPOSE(F)			
			DO i_vec = 1,16
			    DN0_T(i_vec,1) = DN0(i_vec)
			    DO j_vec = 1,16
			        temp01(i_vec,j_vec) = DN0(i_vec)*DN0(j_vec)
			        temp02(i_vec,j_vec) = DN0(i_vec)*DNt(j_vec)
			        temp03(i_vec,j_vec) = DNt(i_vec)*DNt(j_vec)
			    ENDDO
			ENDDO
			
			Q(1,1) = qrr
			Q(1,2) = 0.0d0
			Q(2,1) = 0.0d0
			Q(2,2) = qtt
			
            temp01_mult = MATMUL(D,B)
            temp02_mult = MATMUL(Q,F)
            temp04 = MATMUL(B_T,temp01_mult)
            temp05 = MATMUL(F_T,temp02_mult)			

			temp01_sum = temp01_sum + (gauss_wt(m)*gauss_wt(n)*
     &					 Jacb*den_disk*t_disk*temp01)
			temp02_sum = temp02_sum + (gauss_wt(m)*gauss_wt(n)*
     &					 Jacb*den_disk*t_disk*
     &					 omg*(4.0d0/del_ang)*temp02)
			temp03_sum = temp03_sum - (gauss_wt(m)*gauss_wt(n)*
     &					 Jacb*den_disk*t_disk*
     &					 (omg**2)*(4.0d0/(del_ang**2))*temp03)
			temp04_sum = temp04_sum + (gauss_wt(m)*gauss_wt(n)*
     &					 Jacb*temp04)
			temp05_sum = temp05_sum + (gauss_wt(m)*gauss_wt(n)*
     &					 Jacb*temp05)
			
			ENDDO
		ENDDO

c	  Element level system matrices
c	  =============================
	  M_elem = temp01_sum
		G_elem = temp02_sum
		K_elem = temp03_sum + temp04_sum + temp05_sum

c       Mapping for rows/columns A_elem -> A_glob
c	  =========================================
		DO i = 1 , 4
			DO j = 1 , 4
				mapvec(((4*(j-1))+i)) = 4*(maplg(i)-1)+j
			ENDDO
		ENDDO

c	  Assembly in global matrix
c	  =========================
		DO m = 1 , 16
			DO n = 1 , 16
			  disk_M(mapvec(m),mapvec(n))= disk_M(mapvec(m),mapvec(n))
     &									    + M_elem(m,n)
		      disk_G(mapvec(m),mapvec(n))= disk_G(mapvec(m),mapvec(n))
     &									    + G_elem(m,n)
			  disk_K(mapvec(m),mapvec(n))= disk_K(mapvec(m),mapvec(n))
     &									    + K_elem(m,n)
              mat_pop(mapvec(m),mapvec(n)) = 1.0d0
     
			ENDDO
		ENDDO
	ENDDO
	
c     Common sparse map for disk Matrices
c     ===================================
	CALL sparse_mat_map(mat_pop,disk_map,no_disk_dofs,no_disk_dofs)

	DEALLOCATE(gauss_wt)
	DEALLOCATE(gauss_pt)
	DEALLOCATE(mat_pop)
	
      RETURN
      END

c======================================================================
      SUBROUTINE newmark_disk
c======================================================================
 
c     Shared Data
c     ===========
      USE sldr_dynm
      USE disk_data
      USE disk_aray
      USE disk_soln 
      
      IMPLICIT NONE
      
c     Local variables
c     ===============
      INTEGER i
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_mat
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_disk
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: C_disk_T
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_disk_T
	
c     Newmark Beta Formulation
c	========================
	ALLOCATE(A_mat(no_disk_dofs,no_disk_dofs))
	A_mat = 0.0d0   ! Intialize
 	A_mat = disk_M + (disk_nb_1*dt*disk_G) +
     &	(0.5d0*disk_nb_2*(dt**2)*disk_K)

c     C_disk_T for constrained formulation
c     ====================================	
	ALLOCATE(C_disk_T(no_disk_dofs,no_fxddof))
	C_disk_T = 0.0d0
	C_disk_T(1:no_disk_dofs,1:no_fxddof) = 
     &	TRANSPOSE(C_disk(1:no_fxddof,1:no_disk_dofs))

c	Constraint BC at ID (Lagrange Multiplier Method)
c	================================================
	ALLOCATE(A_disk(N_clamp,N_clamp))
	A_disk = 0.0d0

	A_disk(1:no_disk_dofs,1:no_disk_dofs) = 
     &	A_mat(1:no_disk_dofs,1:no_disk_dofs)
	
 	A_disk(no_disk_dofs+1:N_clamp,1:no_disk_dofs) = 
     &	C_disk(1:no_fxddof,1:no_disk_dofs)

 	A_disk(1:no_disk_dofs,no_disk_dofs+1:N_clamp) = 
     &	C_disk_T(1:no_disk_dofs,1:no_fxddof)
 	
c     LU Factorization for the constant A_disk	
c     =========================================	
	ALLOCATE(L_disk(N_clamp,N_clamp))
	ALLOCATE(U_disk(N_clamp,N_clamp))

	L_disk = 0.0d0
	U_disk = 0.0d0

	ALLOCATE(L_disk_map(N_clamp,N_clamp+1))
	ALLOCATE(U_disk_map(N_clamp,N_clamp+1))
	
	L_disk_map = 0
	U_disk_map = 0
	
 	CALL LU_sparse_fact(A_disk,L_disk,U_disk,
     &  L_disk_map,U_disk_map,N_clamp)

c	U_T for Block Factorization
c	===========================
      ALLOCATE(U_disk_T(N_clamp,N_clamp))
	U_disk_T(1:N_clamp,1:N_clamp) = 
     &	TRANSPOSE(U_disk(1:N_clamp,1:N_clamp))

      ALLOCATE(U_disk_T_map(N_clamp,N_clamp+1))
	CALL sparse_mat_map(U_disk_T,U_disk_T_map,
     &	N_clamp,N_clamp)
     
      DEALLOCATE(U_disk_T)
     
c	Initialize Dynamic Variables
c	============================
	ALLOCATE(w_disk_cur(no_disk_dofs))
	ALLOCATE(v_disk_cur(no_disk_dofs))
	ALLOCATE(a_disk_cur(no_disk_dofs))
	ALLOCATE(w_disk_prv(no_disk_dofs))
	ALLOCATE(v_disk_prv(no_disk_dofs))
	ALLOCATE(a_disk_prv(no_disk_dofs))
	ALLOCATE(f_disk_cnst(no_fxddof)) 
	
c     Disk solution variables
c     =======================
      ALLOCATE(w_disk_prd(no_disk_dofs))
	ALLOCATE(v_disk_prd(no_disk_dofs))
	ALLOCATE(f_disk_val(no_disk_dofs))

	ALLOCATE(disk_vec_1(no_disk_dofs))
	ALLOCATE(disk_vec_2(no_disk_dofs))
	ALLOCATE(disk_vec_3(no_fxddof))
	
	ALLOCATE(x_disk(N_clamp))
	ALLOCATE(y_disk(N_clamp))
	ALLOCATE(z_disk(N_clamp))
	
	ALLOCATE(ramp_coef_disk(no_disk_dofs))
	ALLOCATE(C_disk_r(N_clamp))
	ALLOCATE(L_21_disk(N_clamp))
	ALLOCATE(U_12_disk(N_clamp))
		
	C_disk_r = 0.0d0
	
c     Deallocate temporary matrices
c     =============================	
	DEALLOCATE (A_mat)
	DEALLOCATE (A_disk)
	DEALLOCATE (C_disk_T)
	
      RETURN
      END 

c======================================================================
	SUBROUTINE solve_disk
c======================================================================

c----------------------------------------------------------------------
c	Compute the disk respone at next time step
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE sldr_dynm
      USE shck_data
      USE disk_data
      USE disk_aray
      USE disk_soln
      
      IMPLICIT NONE
      
      DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979

c	Local Variables
c	===============
      INTEGER i,j
      INTEGER rad_gridpt(2),ang_gridpt(2)
      DOUBLE PRECISION time_fac
      DOUBLE PRECISION w_cnst,w_ramp,U_22,lmbd,nu_val
	DOUBLE PRECISION zeta,eta,zeta_temp
	DOUBLE PRECISION ang_cur,ang_vel,eps_ang
 
c	Predicted displacement
c	======================
	w_disk_prd = w_disk_prv + (dt*v_disk_prv) + 
     &	(0.5d0*(1.0d0-disk_nb_2)*(dt**2)*a_disk_prv)

c	Predicted velocity
c	==================
	v_disk_prd = v_disk_prv + ((1.0d0-disk_nb_1)*dt*a_disk_prv)
		
c	Inertia load Vector
c	===================
      f_disk_val  = 0.0d0
	DO i = 1 , no_disk_dofs
		DO j = 1 , no_disk_node
 			f_disk_val(i) = f_disk_val(i) - disk_M(i,4*j-3)*acc_val
		ENDDO
	ENDDO

c	Newmark Beta without Ramp Contact
c	=================================
	CALL d_mult_sparse_mat_vec(disk_G,disk_map,v_disk_prd,
     &	disk_vec_1,no_disk_dofs,no_disk_dofs,no_disk_dofs)
	CALL d_mult_sparse_mat_vec(disk_K,disk_map,w_disk_prd,
     &	disk_vec_2,no_disk_dofs,no_disk_dofs,no_disk_dofs)
      CALL d_mult_sparse_mat_vec(C_disk,C_disk_map,w_disk_prd,
     &	disk_vec_3,no_fxddof,no_disk_dofs,no_disk_dofs)
      time_fac = -(2.0d0/((dt**2)*disk_nb_2))
      disk_vec_3 = time_fac*disk_vec_3
	
	y_disk(1:no_disk_dofs) = f_disk_val - disk_vec_1 - disk_vec_2
	y_disk(no_disk_dofs+1:N_clamp) = disk_vec_3
	
c     Solution of system of equations
c     ===============================
      CALL LU_solve_sparse(L_disk,U_disk,L_disk_map,U_disk_map,
     &	    y_disk(1:N_clamp),x_disk(1:N_clamp),N_clamp)

c	Compute Current Kinematic Value
c	===============================
      a_disk_cur = x_disk(1:no_disk_dofs)
 	v_disk_cur = v_disk_prd + (disk_nb_1*dt*a_disk_cur)
	w_disk_cur = w_disk_prd + (0.5d0*disk_nb_2*(dt**2)*a_disk_cur)

c	Constraint Force
c	================
      f_disk_cnst(1:no_fxddof) = x_disk(no_disk_dofs+1:N_clamp)
	f_disk_cnst(no_fxddof) = 0.0d0

c	Ramp Contact
c	============
 	IF (ramp_opt .NE. 0) THEN

c		Flag for ramp contact
c		=====================
		i_rampcnt = 0
	  ramp_cntf = 0.0d0
	  
c		Ramp Constraint Value
c		=====================
		w_cnst = 0.0d0
		
c	  Tolerance
c	  =========
	  eps_ang = 1.0d-4*del_ang		
		
c	  Angular Speed
c	  =============
	  ang_vel = disk_RPM*(2.0d0*pi/60.0d0)

c		Ramp Location
c		=============
		ang_cur = ramp_ang  - (ang_vel*t_val)
		
c       0 <= Current angle <= 2*pi
c       ==========================
        IF (ABS(ang_cur) .LT. eps_ang ) THEN 
            ang_cur = 0.0d0
        ELSEIF (ang_cur .GT. 0.0d0) THEN
            DO WHILE (ang_cur .GT. (2.0d0*pi)) 
                ang_cur = ang_cur - 2.0d0*pi
            ENDDO
        ELSEIF (ang_cur .LT. 0.0d0) THEN 
            DO WHILE (ang_cur .LT. 0.0d0) 
                ang_cur = ang_cur + 2.0d0*pi
            ENDDO
        ENDIF

c		Position of ramp lowest point in disk (rotating) frame
c		======================================================
		CALL rad_loc(ramp_rad,rad_grid,rad_gridpt,zeta)
		CALL ang_loc(ang_cur,ang_grid,ang_gridpt,eta)
		zeta_temp = zeta

c		DOF Coefficient for w_ramp displacement
c		=======================================
  		CALL ramp_cnstrnt(rad_gridpt,ang_gridpt,
     &		zeta_temp,eta,ramp_coef_disk)

c		Disk Displacement under the ramp
c		================================
		w_ramp = DOT_PRODUCT(ramp_coef_disk(1:no_disk_dofs),
     &		w_disk_cur(1:no_disk_dofs))

		IF (ABS(w_ramp) .GT. ramp_clrnc) THEN

c			Case 1 : Only top ramp
c			======================
			IF (w_ramp .GT. 0.0d0) THEN
				i_rampcnt = 1
				w_cnst = ramp_clrnc
				
c			Case 2 : Both ramp present
c			==========================
			ELSEIF (ramp_opt .EQ. 2)  THEN
				i_rampcnt = 1
				w_cnst = -ramp_clrnc

			ENDIF
		ENDIF
			
		IF (i_rampcnt .EQ. 1) THEN
			
			w_ramp = DOT_PRODUCT(ramp_coef_disk(1:no_disk_dofs),
     &		    w_disk_prd(1:no_disk_dofs))

			lmbd = (2.0d0/((dt**2)*disk_nb_2))* (w_cnst - w_ramp)
			C_disk_r(1:no_disk_dofs) = ramp_coef_disk(1:no_disk_dofs)

c           Solution of system of equations
c           ===============================
            CALL L_solve_sparse(L_disk,L_disk_map,
     &			    C_disk_r,U_12_disk,N_clamp)
			CALL U_T_solve_sparse(U_disk,U_disk_T_map,
     &			    C_disk_r,L_21_disk,N_clamp)
			U_22 = DOT_PRODUCT(L_21_disk(1:N_clamp),
     &		        U_12_disk(1:N_clamp))

c           Solution of system of equations
c           ===============================

c           Forward : L solution
            CALL L_solve_sparse(L_disk,L_disk_map,y_disk,
     &           z_disk,N_clamp)
            nu_val = lmbd - DOT_PRODUCT(L_21_disk,z_disk)
            
c           Ramp Contact Force
            ramp_cntf = -nu_val/U_22

c           Backward : U solution
            DO i = 1 , N_clamp
                z_disk(i) = z_disk(i) - ramp_cntf*U_12_disk(i)
            ENDDO
            CALL U_solve_sparse(U_disk,U_disk_map,z_disk,
     &           x_disk,N_clamp)
            			
c			Compute Current Kinematic Value
c			===============================
			a_disk_cur = x_disk(1:no_disk_dofs)
 			v_disk_cur = v_disk_prd + (disk_nb_1*dt*a_disk_cur)
			w_disk_cur = w_disk_prd + 
     &		          (0.5d0*disk_nb_2*(dt**2)*a_disk_cur)

c			Constraint Force
c			================
            f_disk_cnst(1:no_fxddof) = x_disk(no_disk_dofs+1:N_clamp)

		ENDIF
 	ENDIF

      RETURN
	END

c======================================================================
	SUBROUTINE disk_attitude(rad_pos,ang_pos,w_trn,Dw_rad,Dw_ang)
c======================================================================

c----------------------------------------------------------------------
c	Compute Disk Displacement at specified point
c	a)w   b)w_rad   c)w_ang
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE disk_data
      USE disk_aray
      
      IMPLICIT NONE

      DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979

c	Arguments
c	=========
	DOUBLE PRECISION rad_pos,ang_pos
	DOUBLE PRECISION w_trn,Dw_rad,Dw_ang

c	Local Variable
c	==============
	INTEGER i_node,i_dof
	INTEGER rad_no,ang_no
	INTEGER node_no,elem_no
	INTEGER rad_gridpt(2),ang_gridpt(2)
	DOUBLE PRECISION ang_cur,ang_vel
	DOUBLE PRECISION eps,zeta,eta
	DOUBLE PRECISION eps_ang
	DOUBLE PRECISION w_elem(16),N0(16)
	DOUBLE PRECISION Nz(16),Nr(16),Ne(16),Na(16)
      
c	Tolerance
c	=========
	eps_ang = 1.0d-4*del_ang

c	Tolrance for natural coordinates
c	================================
	eps = 1.0d-2

c	Angular Speed
c	=============
	ang_vel = disk_RPM*(2.0d0*pi/60.0d0)

c	Ramp Location
c	=============
	ang_cur = ang_pos  - (ang_vel*t_val)

c	CCW Rotation
c	============
	IF  (disk_RPM .GT. 0.0d0) THEN
		IF (ABS(ang_cur) .GT. eps_ang) THEN
			DO WHILE (ang_cur .LT. 0.0d0)
				ang_cur = ang_cur + (2.0d0*pi)					
			ENDDO
		ENDIF
c	CW Rotation
c	===========
	ELSE
		IF (ABS(ang_cur - (2.0d0*pi)) .GT. eps_ang) THEN
			DO WHILE (ang_cur .GT. 2.0d0*pi)
				ang_cur = ang_cur - (2.0d0*pi)
			ENDDO
		ENDIF
	ENDIF

c	Point of Interest Position in disk (rotating) frame
c	===================================================
	CALL rad_loc(rad_pos,rad_grid,rad_gridpt,zeta)
	CALL ang_loc(ang_cur,ang_grid,ang_gridpt,eta)

c	Initialize
c	==========
	N0 = 0.0d0
	Nz = 0.0d0
	Nr = 0.0d0
	Ne = 0.0d0
	Na = 0.0d0
	w_elem(16) = 0.0d0

c	Interpolation case
c	==================

	IF (ABS(zeta + 1.0d0) .GT. eps) THEN 
		
c	  Case 1: Lying inside the element
c	  ================================

c	  Case 3: Lying on radial edge
c	  ============================
	  rad_no  = rad_gridpt(1)
		ang_no  = ang_gridpt(1)
		elem_no = (no_rad_elem*(ang_no - 1)) + rad_no
		
c		Shape Function Coefficients
c		===========================
		CALL   N_annular(zeta,eta,del_rad,del_ang,N0)
		CALL DNr_annular(zeta,eta,del_rad,del_ang,Nz)
		CALL DNt_annular(zeta,eta,del_rad,del_ang,Ne)

		Nr = (2/del_rad)*Nz				! J = (2/del_rad)
		Na = (2/(rad_pos*del_ang))*Ne	! J = (2/del_ang)

		DO i_node = 1 , 4
			node_no = disk_elem(elem_no,i_node)
			DO i_dof = 1 , 4
				w_elem(i_node + (4*(i_dof-1))) = 
     &				w_disk_cur((4*(node_no-1)) + i_dof)
			ENDDO
		ENDDO

		w_trn  = DOT_PRODUCT(N0(1:16),w_elem(1:16))
		Dw_rad = DOT_PRODUCT(Nr(1:16),w_elem(1:16))
		Dw_ang = DOT_PRODUCT(Na(1:16),w_elem(1:16))

	ELSE IF (ABS( eta + 1.0d0) .GT. eps) THEN

c	Case 2: Lying on circumfrential edge
c	====================================
c
c		Lying at OD	
c		===========	 
		IF (rad_gridpt(1) .EQ. no_rad_node) THEN
			zeta    = 1.0d0
			rad_no  = rad_gridpt(1)
			ang_no  = ang_gridpt(1)
			elem_no = (no_rad_elem*(ang_no - 1)) + rad_no - 1
		ELSE
			rad_no  = rad_gridpt(1)
			ang_no  = ang_gridpt(1)
			elem_no = (no_rad_elem*(ang_no - 1)) + rad_no
		ENDIF
		
c		Shape Function Coefficients
c		===========================
		CALL   N_annular(zeta,eta,del_rad,del_ang,N0)
		CALL DNr_annular(zeta,eta,del_rad,del_ang,Nz)
		CALL DNt_annular(zeta,eta,del_rad,del_ang,Ne)

		Nr = (2/del_rad)*Nz				! J = (2/del_rad)
		Na = (2/(rad_pos*del_ang))*Ne	! J = (2/del_ang)


		DO i_node = 1 , 4
			node_no = disk_elem(elem_no,i_node)
			DO i_dof = 1 , 4
				w_elem(i_node + (4*(i_dof-1))) = 
     &				w_disk_cur((4*(node_no-1)) + i_dof)
			ENDDO
		ENDDO

		w_trn  = DOT_PRODUCT(N0(1:16),w_elem(1:16))
		Dw_rad = DOT_PRODUCT(Nr(1:16),w_elem(1:16))
		Dw_ang = DOT_PRODUCT(Na(1:16),w_elem(1:16))			
				
	ELSE

c	Case 4: Overlapping with grid node
c	==================================
	  rad_no  = rad_gridpt(1)
		ang_no  = ang_gridpt(1)

		IF (rad_no .EQ. 1) THEN	! Node on ID
			 w_trn = 0.0d0 
			Dw_rad = 0.0d0
			Dw_ang = 0.0d0

		ELSE IF (rad_no .EQ. no_rad_node) THEN	! Node on 0D
			elem_no = (no_rad_elem*(ang_no - 1)) + rad_no - 1
			zeta    = 1.0d0

c			Shape Function Coefficients
c			===========================
			CALL   N_annular(zeta,eta,del_rad,del_ang,N0)
			CALL DNr_annular(zeta,eta,del_rad,del_ang,Nz)
			CALL DNt_annular(zeta,eta,del_rad,del_ang,Ne)

			Nr = (2/del_rad)*Nz				! J = (2/del_rad)
			Na = (2/(rad_pos*del_ang))*Ne	! J = (2/del_ang)

			DO i_node = 1 , 4
				node_no = disk_elem(elem_no,i_node)
				DO i_dof = 1 , 4
					w_elem(i_node + (4*(i_dof-1))) = 
     &				w_disk_cur((4*(node_no-1)) + i_dof)
				ENDDO
			ENDDO
			
			w_trn  = DOT_PRODUCT(N0(1:16),w_elem(1:16))
			Dw_rad = DOT_PRODUCT(Nr(1:16),w_elem(1:16))
			Dw_ang = DOT_PRODUCT(Na(1:16),w_elem(1:16))

		ELSE 

			elem_no = (no_rad_elem*(ang_no - 1)) + rad_no

c			Shape Function Coefficients
c			===========================
			CALL   N_annular(zeta,eta,del_rad,del_ang,N0)
			CALL DNr_annular(zeta,eta,del_rad,del_ang,Nz)
			CALL DNt_annular(zeta,eta,del_rad,del_ang,Ne)

			Nr = (2/del_rad)*Nz				! J = (2/del_rad)
			Na = (2/(rad_pos*del_ang))*Ne	! J = (2/del_ang)

			DO i_node = 1 , 4
				node_no = disk_elem(elem_no,i_node)
				DO i_dof = 1 , 4
					w_elem(i_node + (4*(i_dof-1))) = 
     &				w_disk_cur((4*(node_no-1)) + i_dof)
				ENDDO
			ENDDO

			w_trn  = DOT_PRODUCT(N0(1:16),w_elem(1:16))
			Dw_rad = DOT_PRODUCT(Nr(1:16),w_elem(1:16))
			Dw_ang = DOT_PRODUCT(Na(1:16),w_elem(1:16))

		ENDIF

	ENDIF

      RETURN
	END

c=====================================================================
 	SUBROUTINE rad_loc(rad_pos,rad_grid,rad_gridpt,zeta)
c=====================================================================

c---------------------------------------------------------------------
c	Radial position of a point in mesh
c	==================================
c	rad_gridpt: rad_gridpt(1) <= rad_pos <= rad_gridpt(2)
c	zeta: Natural Coordinate (-1 to 1)
c---------------------------------------------------------------------

c	Global Variables
c	================
      USE disk_data
      
      IMPLICIT NONE

c     Arguments
c     =========     
	INTEGER rad_gridpt(2)
	DOUBLE PRECISION rad_pos,zeta
	DOUBLE PRECISION rad_grid(no_rad_node)

c	Local Variables
c	===============
	INTEGER cnt_rnode,rflag
	DOUBLE PRECISION eps_rad,r0,r1,r2

c     Tolerance
c     =========
	eps_rad = 1.0d-4*del_rad

	IF (ABS(rad_pos - rad_grid(1)) < eps_rad) THEN
	    rad_gridpt(1) = 1
		rad_gridpt(2) = 1
		zeta = -1.0d0
	ELSE
	    cnt_rnode = 1 ! First node already accounted
		rflag = 0       ! Status Flag: 0 == Interval not found
		DO WHILE ((rflag == 0) .AND. (cnt_rnode < no_rad_node))
			cnt_rnode = cnt_rnode + 1
			IF ((rad_grid(cnt_rnode) - rad_pos) > eps_rad) THEN
				rad_gridpt(1) = cnt_rnode-1
				rad_gridpt(2) = cnt_rnode
				rflag = 1
				r0 = rad_pos
				r1 = rad_grid(cnt_rnode-1)
				r2 = rad_grid(cnt_rnode)
				zeta = (2.0d0*r0 - (r1 + r2))/(r2-r1)
            
			ELSE IF (ABS(rad_pos - rad_grid(cnt_rnode))
     &		        < eps_rad) THEN
					rad_gridpt(1) = cnt_rnode
					rad_gridpt(2) = cnt_rnode
					rflag = 1
					zeta = -1.0d0
			ENDIF
		ENDDO
	ENDIF
      
      RETURN
	END

c=====================================================================
 	SUBROUTINE ang_loc(ang_pos,ang_grid,ang_gridpt,eta)
c=====================================================================

c---------------------------------------------------------------------
c	Angular position of a point in mesh
c	===================================
c	ang_gridpt: ang_gridpt(1) <= ang_pos <= ang_gridpt(2)
c	eta: Natural Coordinate (-1 to 1)
c---------------------------------------------------------------------

c	Global Variables
c	================
      USE disk_data
      
      IMPLICIT NONE
	
	DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979

c     Arguments
c     =========
	INTEGER ang_gridpt(2)
	DOUBLE PRECISION ang_pos,eta
	DOUBLE PRECISION ang_grid(no_ang_node)

c	Local Variables
c	===============
	INTEGER cnt_anode,aflag
	DOUBLE PRECISION eps_ang,a0,a1,a2

	eps_ang = 1.0d-4*del_ang	! Tolerance

	cnt_anode = 0
	aflag = 0

c	Exclude the last sector
	DO WHILE ((aflag == 0) .AND. (cnt_anode < no_ang_node-1))  
		cnt_anode = cnt_anode + 1
	    IF (ABS(ang_pos - ang_grid(cnt_anode)) < eps_ang) THEN
			ang_gridpt(1) = cnt_anode
			ang_gridpt(2) = cnt_anode
			eta = -1.0d0
			aflag = 1
		ELSE IF ((ang_grid(cnt_anode+1) - ang_pos) > eps_ang) THEN
			ang_gridpt(1) = cnt_anode
			ang_gridpt(2) = cnt_anode + 1
			a0 = ang_pos
			a1 = ang_grid(cnt_anode)
			a2 = ang_grid(cnt_anode+1)
			eta = (2.0d0*a0 - (a1 + a2))/(a2-a1)
			aflag = 1
		ENDIF
	ENDDO

c	Node present in last sector
	IF (aflag == 0) THEN
		cnt_anode = cnt_anode + 1

		IF (ABS(ang_pos - ang_grid(cnt_anode)) < eps_ang) THEN
			ang_gridpt(1) = cnt_anode
			ang_gridpt(2) = cnt_anode
			eta = -1.0d0
			aflag = 1
		ELSE IF (ABS(ang_pos - 2.0d0*pi) < eps_ang) THEN
			ang_gridpt(1) = 1
			ang_gridpt(2) = 1
			eta = -1.0d0
			aflag = 1
		ELSE IF (ang_pos < 2.0d0*pi) THEN
			ang_gridpt(1) = cnt_anode
			ang_gridpt(2) = 1.0d0
			a0 = ang_pos
			a1 = ang_grid(cnt_anode)
			a2 = 2.0d0*pi
			eta = (2.0d0*a0 - (a1 + a2))/(a2-a1)
			aflag = 1
		ENDIF
	ENDIF

      RETURN
	END

c=====================================================================
 	SUBROUTINE ramp_cnstrnt(rad_gridpt,ang_gridpt,zeta,eta,ramp_coef)
c=====================================================================

c---------------------------------------------------------------------
c	Constraint Vector for ramp contact
c	==================================.
c	w_ramp = DOT(ramp_coef.w_dof)
c	rad_gridpt: rad_gridpt(1) <= rad_pos <= rad_gridpt(2)
c	ang_gridpt: ang_gridpt(1) <= ang_pos <= ang_gridpt(2)
c	zeta: Natural Coordinate (-1 to 1)
c	eta : Natural Coordinate (-1 to 1)
c---------------------------------------------------------------------

c	Shared Data
c	===========
      USE disk_data
      USE disk_aray
      
      IMPLICIT NONE

c     Arguments
c     =========
	INTEGER rad_gridpt(2),ang_gridpt(2)
	DOUBLE PRECISION zeta,eta
	DOUBLE PRECISION ramp_coef(no_disk_dofs)
     
c	Local Variables
c	===============
	INTEGER i_node,i_dof
	INTEGER rad_no,ang_no
	INTEGER node_no,elem_no
	DOUBLE PRECISION eps
	DOUBLE PRECISION elem_coef(16)

c	Initialize
c	==========
	ramp_coef = 0.0d0
	elem_coef = 0.0d0

c	Tolrance for natural coordinates
c	================================	
	eps = 1.0d-2

c	Interpolation case
c	==================

	IF (ABS(zeta + 1.0d0) .GT. eps) THEN 
		
c	Case 1: Lying inside the element
c	================================

c	Case 3: Lying on radial edge
c	============================
	  rad_no  = rad_gridpt(1)
		ang_no  = ang_gridpt(1)
		elem_no = (no_rad_elem*(ang_no - 1)) + rad_no
		CALL N_annular(zeta,eta,del_rad,del_ang,elem_coef)
		DO i_node = 1 , 4
			node_no = disk_elem(elem_no,i_node)
			DO i_dof = 1 , 4
				ramp_coef((4*(node_no - 1)) + i_dof ) =
     &            elem_coef((4*(i_dof   - 1)) + i_node)
			ENDDO
		ENDDO

	ELSE IF (ABS( eta + 1.0d0) .GT. eps) THEN

c	Case 2: Lying on circumfrential edge
c	====================================
c
c		Lying at OD	
c		===========	 
		IF (rad_gridpt(1) .EQ. no_rad_node) THEN
			zeta    = 1.0d0
			rad_no  = rad_gridpt(1)
			ang_no  = ang_gridpt(1)
			elem_no = (no_rad_elem*(ang_no - 1)) + rad_no - 1
		ELSE
			rad_no  = rad_gridpt(1)
			ang_no  = ang_gridpt(1)
			elem_no = (no_rad_elem*(ang_no - 1)) + rad_no
		ENDIF
		
		CALL N_annular(zeta,eta,del_rad,del_ang,elem_coef)
		DO i_node = 1 , 4
			node_no = disk_elem(elem_no,i_node)
			DO i_dof = 1 , 4
				ramp_coef((4*(node_no - 1)) + i_dof ) =
     &            elem_coef((4*(i_dof   - 1)) + i_node)
			ENDDO
		ENDDO		
				
	ELSE

c	Case 4: Overlapping with grid node
c	==================================
	  rad_no  = rad_gridpt(1)
		ang_no  = ang_gridpt(1)
		node_no = (no_rad_node*(ang_no - 1)) + rad_no
		ramp_coef(4*node_no-3) = 1.0d0

	ENDIF
      
      RETURN
	END
	
c======================================================================
      SUBROUTINE intrp_disk_atd(disp_cur,disp_prv,no_stps,disk_atd)
c======================================================================

      IMPLICIT NONE
      
c     Arguments
c     =========      
      DOUBLE PRECISION disp_cur(3),disp_prv(3)
      INTEGER no_stps
      DOUBLE PRECISION disk_atd(no_stps,3)
      
c     Local Variable
c     ==============
      INTEGER i, i_fine
      DOUBLE PRECISION del_disp(3)
      DOUBLE PRECISION f_dtfac
      
      f_dtfac = DBLE(no_stps)
      
      DO i = 1, 3
        del_disp(i) = (disp_cur(i) - disp_prv(i))/f_dtfac
      ENDDO

      DO i_fine = 1, no_stps
        DO i = 1, 3
            disk_atd(i_fine,i) =  disp_prv(i) + 
     &                           (del_disp(i)*DBLE(i_fine))
      	ENDDO
      ENDDO	    
    
      RETURN
	END
	
c====================================================================== 