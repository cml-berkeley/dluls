c**********************************************************************
c     Subroutines in this file :
c     01. mount_disk
c     02. newmark_rotor
c     03. solve_rotor
c     04. init_housing
c     05. couple_rotr_hsng
c     06. newmark_motor
c     07. solve_motor
c     08. init_base
c     09. couple_hsng_base
c     10. couple_rotr_base
c     11. newmark_support
c     12. solve_support
c     13. calc_hub_inrt
c     14. crcm_fdb_intp 
c     15. nodal_FDB_dyn_coef
c**********************************************************************

c======================================================================
      SUBROUTINE mount_disk
c======================================================================

c     Shared data
c     ===========
      USE shck_data
      USE disk_data
      USE disk_aray
      USE motr_data
      USE rotr_data
      
      IMPLICIT NONE
      
c     Local variables
c     ===============
      INTEGER i,j,i_row,j_col,col_no,row_no
      INTEGER i_pvt,i_cnt,no_nodes,row_cnt
      INTEGER i_del_node,node_no,dof_no,rotr_dof
      
c     Temporary vectors/matrices
c     ==========================
      INTEGER, DIMENSION (:), ALLOCATABLE :: rearrange_index
      
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: temp_mat
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: disk_inrt_acc
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: disk_inrt_vec
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: rotr_inrt_map

      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: M_row
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: G_row
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: K_row
      
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: M_col
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: G_col
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: K_col
 
c     Dofs calculation
c     ================
c     Disk slope at ID is zero
      no_disk_bcs = 3*no_ang_node
      
c     DOFs for rotating part
      no_rotr_dofs = no_hub_dofs + no_disk_dofs - no_disk_intf_dofs
      
c     Size of Newmark matrices
      no_rotr_vrbl = no_disk_bcs + no_rotr_dofs 
      
c     Initialize inertia load map
c     ===========================
      ALLOCATE(rotr_inrt_map(no_rotr_dofs,3))
      rotr_inrt_map = 0.0d0 
            
c     Checking node compatibility      
c     ===========================      
      IF (MOD(no_ang_node,no_disk_intf_dofs) .NE. 0) THEN
        WRITE(*,*)
        WRITE(*,*) '      No. of circumferential nodes are not '
        WRITE(*,*) '      intergral multiple of no. of interface '
        WRITE(*,*) '      nodes hence simulation terminated....  '
        STOP
      ENDIF

c     Disk-Rotor interface list
c     =========================    
      i_del_node = no_ang_node/no_disk_intf_dofs
      DO i = 1 , no_disk_intf_dofs
        node_no = (i_del_node*no_rad_node)*(i-1) + 1
        disk_intf_dofs(i,2) = 4*node_no - 3
      ENDDO
      
c     Map for rearranging disk dofs
c     =============================
      ALLOCATE(disk_intf_map(no_disk_dofs,2))
      ALLOCATE(rearrange_index(no_disk_dofs))
      
      disk_intf_map = 0
      rearrange_index = 0
      
      DO i = 1 , no_disk_intf_dofs
        rearrange_index(disk_intf_dofs(i,2)) = 1
      ENDDO

c     Rearranging DOFs
c     ================
      i_cnt = 0
      
c     Non-interface dofs
      DO i = 1 , no_disk_dofs
        IF (rearrange_index(i) .EQ. 0) THEN
            i_cnt = i_cnt + 1
            disk_intf_map(i_cnt,1) = i
        ENDIF 
      ENDDO
      
c     Interface dofs
      DO i = 1 , no_disk_intf_dofs
        disk_intf_map(i_cnt+i,1) = disk_intf_dofs(i,2)
      ENDDO
      
c     Inverse disk map for implementing boundary condtions
c     ====================================================
      DO i = 1 , no_disk_dofs
        disk_intf_map(disk_intf_map(i,1),2) = i
      ENDDO
     
c     Disk inertia load map
c     =====================
      DO i = 1 , no_disk_node
        dof_no = 4*i-3
        rotr_inrt_map(disk_intf_map(dof_no,2),3) = 1.0d0
      ENDDO

c     Deallocate temprory matrix
c     ==========================
      DEALLOCATE(rearrange_index)
      
c     Rearrange disk matrices
c     =======================
      ALLOCATE(M_row(no_disk_dofs,no_disk_dofs))
      ALLOCATE(G_row(no_disk_dofs,no_disk_dofs))
      ALLOCATE(K_row(no_disk_dofs,no_disk_dofs))
      
      M_row = 0.0d0
      G_row = 0.0d0
      K_row = 0.0d0
      
      ALLOCATE(M_col(no_disk_dofs,no_disk_dofs))
      ALLOCATE(G_col(no_disk_dofs,no_disk_dofs))
      ALLOCATE(K_col(no_disk_dofs,no_disk_dofs))
      
      M_col = 0.0d0
      G_col = 0.0d0
      K_col = 0.0d0
      
c     Rearrange rows
c     ==============
      DO i = 1, no_disk_dofs
        M_row(i,1:no_disk_dofs) = 
     &      disk_M(disk_intf_map(i,1),1:no_disk_dofs)
        G_row(i,1:no_disk_dofs) = 
     &      disk_G(disk_intf_map(i,1),1:no_disk_dofs)
        K_row(i,1:no_disk_dofs) = 
     &      disk_K(disk_intf_map(i,1),1:no_disk_dofs)     
      ENDDO

c     Rearrange columns
c     =================
      DO i = 1, no_disk_dofs
        M_col(1:no_disk_dofs,i) = 
     &      M_row(1:no_disk_dofs,disk_intf_map(i,1))
        G_col(1:no_disk_dofs,i) = 
     &      G_row(1:no_disk_dofs,disk_intf_map(i,1))
        K_col(1:no_disk_dofs,i) = 
     &      K_row(1:no_disk_dofs,disk_intf_map(i,1))
      ENDDO
 
c     Deallocate disk temporary matrices
c     ==================================
      DEALLOCATE(M_row)
      DEALLOCATE(G_row)
      DEALLOCATE(K_row)
     
c     Rearranged disk matrices
c     ========================  
      disk_M = M_col
      disk_G = G_col
      disk_K = K_col
      
c     Deallocate disk temporary matrices
c     ==================================      
      DEALLOCATE(M_col)
      DEALLOCATE(G_col)
      DEALLOCATE(K_col)      

c     Hub structural matrices
c     =======================
      ALLOCATE(hub_M(no_hub_dofs,no_hub_dofs))
      ALLOCATE(hub_G(no_hub_dofs,no_hub_dofs))
      ALLOCATE(hub_K(no_hub_dofs,no_hub_dofs))
      
      hub_M = 0.0d0
      hub_G = 0.0d0
      hub_K = 0.0d0
      
c     Hub coordinate list      
c     ===================
      ALLOCATE(hub_cords(no_hub_dofs,4))
      hub_cords = 0.0d0
            
c     Read input ANSYS data
c     =====================      
      OPEN(81,ERR=999,FILE='hub_M.dat',STATUS='OLD')
      OPEN(82,ERR=999,FILE='hub_K.dat',STATUS='OLD')
      OPEN(83,ERR=999,FILE='hub_cords.dat',STATUS='OLD')

      DO i=1,no_hub_dofs
		READ(81,*) (hub_M(i,j),j=1,no_hub_dofs)
        READ(82,*) (hub_K(i,j),j=1,no_hub_dofs)
        READ(83,*) (hub_cords(i,j),j=1,4)
	ENDDO
	
      CLOSE(81)
      CLOSE(82)
      CLOSE(83)
      
c     Proportional damping matrix
c     ===========================      
      hub_G = (disk_M_damp*hub_M) + (disk_K_damp*hub_K)
      
c     Rearrange hub dofs for coupling
c     ===============================
      ALLOCATE(hub_intf_map(no_hub_dofs,2))
      ALLOCATE(rearrange_index(no_hub_dofs))
      
      hub_intf_map = 0
      rearrange_index = 0
      
c     Disk interface dofs
c     ===================      
      DO i = 1 , no_disk_intf_dofs
        rearrange_index(disk_intf_dofs(i,1)) = 1
      ENDDO
      
c     FDB interface dofs
c     ==================      
      DO i = 1 , no_fdb_dofs_hub
        rearrange_index(hub_fdb_dofs(i)) = 2
      ENDDO

c     Rearranging DOFs
c     ================
      i_cnt = 0

c     Disk interface dofs
      DO i = 1 , no_disk_intf_dofs
        i_cnt = i_cnt + 1
        hub_intf_map(i_cnt,1) = disk_intf_dofs(i,1)
      ENDDO
      
c     Non-interface dofs
      DO i = 1 , no_hub_dofs
        IF (rearrange_index(i) .EQ. 0) THEN
            i_cnt = i_cnt + 1
            hub_intf_map(i_cnt,1) = i
        ENDIF 
      ENDDO
      
c     FDB interface dofs
      DO i = 1 , no_fdb_dofs_hub
        i_cnt = i_cnt + 1
        hub_intf_map(i_cnt,1) = hub_fdb_dofs(i)
      ENDDO

c     Inverse hub map for implementing boundary condtions
c     ===================================================
      DO i = 1 , no_hub_dofs
        hub_intf_map(hub_intf_map(i,1),2) = i
      ENDDO
      
c     Disk interface dofs map
c     =======================
      DO i = 1 , no_disk_intf_dofs
        disk_ID_map(i,3) = no_disk_dofs - no_disk_intf_dofs + 
     &      hub_intf_map(disk_ID_map(i,1),2)
        disk_ID_map(i,4) = no_disk_dofs - no_disk_intf_dofs + 
     &      hub_intf_map(disk_ID_map(i,2),2)
        disk_ID_map(i,5) = no_disk_dofs - no_disk_intf_dofs + i
      ENDDO

c     Deallocate temprory matrix
c     ==========================
      DEALLOCATE(rearrange_index)      
      
c     Rearrange hub matrices
c     ======================
      ALLOCATE(M_row(no_hub_dofs,no_hub_dofs))
      ALLOCATE(G_row(no_hub_dofs,no_hub_dofs))
      ALLOCATE(K_row(no_hub_dofs,no_hub_dofs))
      
      M_row = 0.0d0
      G_row = 0.0d0
      K_row = 0.0d0
      
      ALLOCATE(M_col(no_hub_dofs,no_hub_dofs))
      ALLOCATE(G_col(no_hub_dofs,no_hub_dofs))
      ALLOCATE(K_col(no_hub_dofs,no_hub_dofs))
      
      M_col = 0.0d0
      G_col = 0.0d0
      K_col = 0.0d0

c     Rearrange rows
c     ==============
      DO i = 1, no_hub_dofs
        M_row(i,1:no_hub_dofs) = 
     &      hub_M(hub_intf_map(i,1),1:no_hub_dofs)
        G_row(i,1:no_hub_dofs) = 
     &      hub_G(hub_intf_map(i,1),1:no_hub_dofs)
        K_row(i,1:no_hub_dofs) = 
     &      hub_K(hub_intf_map(i,1),1:no_hub_dofs)
      ENDDO
      
c     Rearrange columns
c     =================
      DO i = 1, no_hub_dofs
        M_col(1:no_hub_dofs,i) =
     &    M_row(1:no_hub_dofs,hub_intf_map(i,1))
        G_col(1:no_hub_dofs,i) =
     &    G_row(1:no_hub_dofs,hub_intf_map(i,1))
        K_col(1:no_hub_dofs,i) =
     &    K_row(1:no_hub_dofs,hub_intf_map(i,1))
      ENDDO
      
c     Rearranged hub matrices
c     =======================  
      hub_M = M_col
      hub_G = G_col
      hub_K = K_col
      
c     Deallocate hub temporary matrices
c     =================================
      DEALLOCATE(M_row)
      DEALLOCATE(G_row)
      DEALLOCATE(K_row)
     
      DEALLOCATE(M_col)
      DEALLOCATE(G_col)
      DEALLOCATE(K_col)         
      
c     Rotor structural matrices
c     =========================
      ALLOCATE(rotr_M(no_rotr_dofs,no_rotr_dofs))
      ALLOCATE(rotr_G(no_rotr_dofs,no_rotr_dofs))
      ALLOCATE(rotr_K(no_rotr_dofs,no_rotr_dofs))
      
      rotr_M = 0.0d0
      rotr_G = 0.0d0
      rotr_K = 0.0d0  

c     Integrate the disk matrices with rotor
c     ======================================
      DO i = 1 , no_disk_dofs
        DO j = 1 , no_disk_dofs
            rotr_M(i,j) = rotr_M(i,j) + disk_M(i,j)
            rotr_G(i,j) = rotr_G(i,j) + disk_G(i,j)
            rotr_K(i,j) = rotr_K(i,j) + disk_K(i,j)
        ENDDO
      ENDDO

c     Deallocate disk matrices
c     ========================
      DEALLOCATE(disk_M)
      DEALLOCATE(disk_G)
      DEALLOCATE(disk_K)
      
c     Integrate the hub matrices with rotor
c     =====================================
      DO i = 1 , no_hub_dofs
        i_row = no_disk_dofs - no_disk_intf_dofs + i
        DO j = 1 , no_hub_dofs
            j_col = no_disk_dofs - no_disk_intf_dofs + j
            rotr_M(i_row,j_col) = rotr_M(i_row,j_col) + hub_M(i,j)
            rotr_G(i_row,j_col) = rotr_G(i_row,j_col) + hub_G(i,j)
            rotr_K(i_row,j_col) = rotr_K(i_row,j_col) + hub_K(i,j)
        ENDDO
      ENDDO   

c     Deallocate hub matrices
c     =======================
c     DEALLOCATE(hub_M) : Used for calculating hub inertia forces
      DEALLOCATE(hub_G)
      DEALLOCATE(hub_K)
           
c     Sparse Map for rotor matrices
c     =============================

c     Inertia matrix (does not get affected by FDB)
      ALLOCATE(rotr_M_map(no_rotr_dofs,no_rotr_dofs+1))
	CALL sparse_mat_map(rotr_M,rotr_M_map,
     &	no_rotr_dofs,no_rotr_dofs)

      ALLOCATE(temp_mat(no_rotr_dofs,no_rotr_dofs))

c     Damping matrix get affected by FDB
      temp_mat = rotr_G
      DO i = 1 , no_fdb_dofs_hub
        row_no = no_rotr_dofs - no_fdb_dofs_hub + i
        DO j = 1 , no_fdb_dofs_hub
            col_no = no_rotr_dofs - no_fdb_dofs_hub + j
            temp_mat(row_no,col_no) = 1.0d0
        ENDDO
      ENDDO
      ALLOCATE(rotr_G_map(no_rotr_dofs,no_rotr_dofs+1))
      CALL sparse_mat_map(temp_mat,rotr_G_map,
     &	no_rotr_dofs,no_rotr_dofs) 
      
c     Stiffness matrix get affected by FDB
      temp_mat = rotr_K
      DO i = 1 , no_fdb_dofs_hub
        row_no = no_rotr_dofs - no_fdb_dofs_hub + i
        DO j = 1 , no_fdb_dofs_hub
            col_no = no_rotr_dofs - no_fdb_dofs_hub + j
            temp_mat(row_no,col_no) = 1.0d0
        ENDDO
      ENDDO      
	ALLOCATE(rotr_K_map(no_rotr_dofs,no_rotr_dofs+1))
	CALL sparse_mat_map(temp_mat,rotr_K_map,
     &	no_rotr_dofs,no_rotr_dofs)
     
      DEALLOCATE(temp_mat)       

c     Rotor inertia loading : Hub
c     ===========================      
      DO i = 1 , no_hub_dofs
        dof_no = INT(hub_cords(i,4))
        rotr_dof =  no_disk_dofs - no_disk_intf_dofs +
     &              hub_intf_map(i,2)
        rotr_inrt_map(rotr_dof,dof_no) = 1.0d0
      ENDDO
      
c     Intialize vector for load calculation
c     =====================================

c     --------------------------------------
c     Disk part of rotor displacement has  
c     only Uz DOF which does not change  
c     with disk rotation
c     --------------------------------------
      m_rotr = no_disk_dofs - no_disk_intf_dofs
      ALLOCATE(temp_mat(m_rotr,m_rotr))
      temp_mat = 0.0d0  ! Intiailize
      
      temp_mat(1:m_rotr,1:m_rotr) = 
     &  rotr_M(1:m_rotr,1:m_rotr)
      
      ALLOCATE(disk_inrt_acc(m_rotr,1))
      disk_inrt_acc = 0.0d0  ! Intiailize
      
      disk_inrt_acc = MATMUL(rotr_inrt_map(1:m_rotr,1:3),shk_dir)
      
      ALLOCATE(disk_inrt_vec(m_rotr,1))
      disk_inrt_vec = 0.0d0
      
      disk_inrt_vec = MATMUL(temp_mat,disk_inrt_acc)
      DEALLOCATE(temp_mat)
      DEALLOCATE(disk_inrt_acc)
      
      ALLOCATE(rotr_inrt_load(no_rotr_dofs))
      rotr_inrt_load = 0.0d0
      
      DO i = 1, m_rotr
        rotr_inrt_load(i) = disk_inrt_vec(i,1)
      ENDDO
      
c     --------------------------------------
c     Hub part of rotor displacement changes
c     with disk rotation
c     --------------------------------------
      
      ALLOCATE(hub_inrt_map(no_hub_dofs,3))
      hub_inrt_map = 0.0d0
      
      hub_inrt_map(1:no_hub_dofs,1:3) = 
     &  rotr_inrt_map(m_rotr+1:no_rotr_dofs,1:3)
      
      DEALLOCATE(rotr_inrt_map)
      
      ALLOCATE(hub_inrt_acc(no_hub_dofs,1))
      hub_inrt_acc = 0.0d0
      
      ALLOCATE(hub_inrt_vec(no_hub_dofs,1))
      hub_inrt_vec = 0.0d0
      
c     Boundary Condition : Clamped disk center
c     ========================================
      ALLOCATE(C_rotr(no_disk_bcs,no_rotr_dofs))
      C_rotr = 0.0d0
      
c     Slope at disk ID is zero
      row_cnt = 0
	DO i = 1 , no_ang_node
	  node_no = ((no_rad_node)*(i-1)) + 1
		DO j = 2 , 4    ! Displacement coupled so not constrained
		    row_cnt = row_cnt + 1
			dof_no  = 4*(node_no - 1) + j
	      col_no  = disk_intf_map(dof_no,2)
			C_rotr(row_cnt,col_no) = 1.0d0
		ENDDO
	ENDDO
	      
	ALLOCATE(C_rotr_map(no_disk_bcs,no_rotr_dofs+1))
	CALL sparse_mat_map(C_rotr,C_rotr_map,
     &	no_disk_bcs,no_rotr_dofs)
      
      RETURN
     
999	WRITE(*,*) 'Trouble in opening hub files'
	STOP
	      
      END

c======================================================================
      SUBROUTINE newmark_rotor
c======================================================================

c     Shared data
c     ===========
      USE sldr_dynm
      USE disk_data
      USE disk_aray
      USE motr_data
      USE rotr_data
      USE rotr_soln
      
      IMPLICIT NONE
      
c     Local variables
c     ===============
      INTEGER i,j,i_brng,i_node,brng_row,pvt
      INTEGER FDB_node_map(3)
      DOUBLE PRECISION K_nodal(3,3),C_nodal(3,3)

	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_mat
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_rotr
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_rotr_T
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: C_rotr_T
	
c     Add FDB Stiffness and damping
c     =============================
      DO i_brng = 1 , 3 ! Bearing type : 3

        brng_row  = 3*(i_brng-1)
        
c       Reset    
        C_nodal = 0.0d0    
        K_nodal = 0.0d0
        
c       Nodal Dynanic Coefficients
        DO i = 1 , 3
            DO j = 1 , 3
                C_nodal(i,j) = C_fdb(brng_row+i,j)/
     &                         DBLE(no_crcm_node_hub)
                K_nodal(i,j) = K_fdb(brng_row+i,j)/
     &                         DBLE(no_crcm_node_hub)
            ENDDO
        ENDDO

c       Add nodal coefficents
        DO i_node = 1 , no_crcm_node_hub
            
c           Find location of inteface FDB node of rotr_dofs list           
            pvt = 3*no_crcm_node_hub*(i_brng-1) + 3*(i_node-1)
            DO i = 1 , 3
                FDB_node_map(i) = no_disk_dofs - no_disk_intf_dofs + 
     &                            hub_intf_map(hub_fdb_dofs(pvt+i),2)
            ENDDO
            
            DO i = 1 , 3
                DO j = 1 , 3
                
                    rotr_G(FDB_node_map(i),FDB_node_map(j)) = 
     &                  rotr_G(FDB_node_map(i),FDB_node_map(j)) +
     &                  C_nodal(i,j)
     
                    rotr_K(FDB_node_map(i),FDB_node_map(j)) = 
     &                  rotr_K(FDB_node_map(i),FDB_node_map(j)) +
     &                  K_nodal(i,j)
     
                ENDDO
            ENDDO
            
        ENDDO
      ENDDO	
	
c     Newmark Beta Formulation
c	========================
	ALLOCATE(A_mat(no_rotr_dofs,no_rotr_dofs))
	A_mat = 0.0d0   ! Intialize
 	A_mat = rotr_M + (disk_nb_1*dt*rotr_G) +
     &	(0.5d0*disk_nb_2*(dt**2)*rotr_K) 
     
c     C_rotr_T for constrained formulation
c     =====================================	
	ALLOCATE(C_rotr_T(no_rotr_dofs,no_disk_bcs))
	C_rotr_T = 0.0d0
	C_rotr_T(1:no_rotr_dofs,1:no_disk_bcs) = 
     &	TRANSPOSE(C_rotr(1:no_disk_bcs,1:no_rotr_dofs))
     
c	Constraint BC at ID (Lagrange Multiplier Method)
c	================================================
	ALLOCATE(A_rotr(no_rotr_vrbl,no_rotr_vrbl))
	A_rotr = 0.0d0

	A_rotr(1:no_rotr_dofs,1:no_rotr_dofs) = 
     &	A_mat(1:no_rotr_dofs,1:no_rotr_dofs)
	
 	A_rotr(no_rotr_dofs+1:no_rotr_vrbl,1:no_rotr_dofs) = 
     &	C_rotr(1:no_disk_bcs,1:no_rotr_dofs)

 	A_rotr(1:no_rotr_dofs,no_rotr_dofs+1:no_rotr_vrbl) = 
     &	C_rotr_T(1:no_rotr_dofs,1:no_disk_bcs)     


c     LU Factorization for the constant A_rotr
c     ========================================	
	ALLOCATE(L_rotr(no_rotr_vrbl,no_rotr_vrbl))
	ALLOCATE(U_rotr(no_rotr_vrbl,no_rotr_vrbl))
	
      L_rotr = 0.0d0
	U_rotr = 0.0d0

	ALLOCATE(L_rotr_map(no_rotr_vrbl,no_rotr_vrbl+1))
	ALLOCATE(U_rotr_map(no_rotr_vrbl,no_rotr_vrbl+1))
	
	L_rotr_map = 0
	U_rotr_map = 0
	
 	CALL LU_sparse_fact(A_rotr,L_rotr,U_rotr,
     &  L_rotr_map,U_rotr_map,no_rotr_vrbl)
	
c	U_T for Block Factorization
c	===========================
      ALLOCATE(U_rotr_T(no_rotr_vrbl,no_rotr_vrbl))
	U_rotr_T(1:no_rotr_vrbl,1:no_rotr_vrbl) = 
     &	TRANSPOSE(U_rotr(1:no_rotr_vrbl,1:no_rotr_vrbl))
     
      ALLOCATE(U_rotr_T_map(no_rotr_vrbl,no_rotr_vrbl+1))
	CALL sparse_mat_map(U_rotr_T,U_rotr_T_map,
     &	no_rotr_vrbl,no_rotr_vrbl)
     
      DEALLOCATE(U_rotr_T)

c	Initialize Dynamic Variables
c	============================
      ALLOCATE(w_disk_cur(no_disk_dofs))
	ALLOCATE(w_rotr_cur(no_rotr_dofs))
	ALLOCATE(v_rotr_cur(no_rotr_dofs))
	ALLOCATE(a_rotr_cur(no_rotr_dofs))
	ALLOCATE(w_rotr_prv(no_rotr_dofs))
	ALLOCATE(v_rotr_prv(no_rotr_dofs))
	ALLOCATE(a_rotr_prv(no_rotr_dofs))
	ALLOCATE(f_rotr_cnst(no_disk_bcs))
	
c     Disk-Hub solution variables
c     ===========================
      ALLOCATE(w_rotr_prd(no_rotr_dofs))
	ALLOCATE(v_rotr_prd(no_rotr_dofs))
	ALLOCATE(f_rotr_val(no_rotr_dofs))

      ALLOCATE(rotr_vec_1(no_rotr_dofs))
	ALLOCATE(rotr_vec_2(no_rotr_dofs))
	ALLOCATE(rotr_vec_3(no_disk_bcs))
	
	ALLOCATE(x_rotr(no_rotr_vrbl))
	ALLOCATE(y_rotr(no_rotr_vrbl))
	ALLOCATE(z_rotr(no_rotr_vrbl))
	
	ALLOCATE(ramp_coef_rotr(no_disk_dofs))
	ALLOCATE(C_rotr_r(no_rotr_vrbl))
	ALLOCATE(L_21_rotr(no_rotr_vrbl))
	ALLOCATE(U_12_rotr(no_rotr_vrbl))
	
	C_rotr_r = 0.0d0
	     
c     Deallocate temporary matrices
c     =============================
      DEALLOCATE(A_mat)
      DEALLOCATE(A_rotr)
      DEALLOCATE(C_rotr_T)
      
      RETURN
      END
      
c======================================================================
      SUBROUTINE solve_rotor
c======================================================================

c     Shared data
c     ===========
      USE sldr_dynm
      USE shck_data
      USE disk_data
      USE disk_aray
      USE motr_data
      USE rotr_data
      USE rotr_soln
      
      IMPLICIT NONE
      
      DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979

c     Local variables
c     ===============
      INTEGER i,j,pvt
      INTEGER rad_gridpt(2),ang_gridpt(2)
      DOUBLE PRECISION time_fac
      DOUBLE PRECISION zeta,eta,zeta_temp
      DOUBLE PRECISION ang_cur,ang_vel,eps_ang
      DOUBLE PRECISION w_cnst,w_ramp,U_22,lmbd,nu_val
	DOUBLE PRECISION xtemp,ytemp,rad_temp,ang_temp
	
c	Angular Speed
c	=============
	ang_vel = disk_RPM*(2.0d0*pi/60.0d0)	

c	Predicted displacement
c	======================
	w_rotr_prd = w_rotr_prv + (dt*v_rotr_prv) + 
     &	(0.5d0*(1.0d0-disk_nb_2)*(dt**2)*a_rotr_prv)
     
c	Predicted velocity
c	==================
	v_rotr_prd = v_rotr_prv + ((1.0d0-disk_nb_1)*dt*a_rotr_prv)

c	Inertia load vector
c	===================
      CALL calc_hub_inrt(ang_vel,t_val)
      
c     Add hub part of rotor inertia vector
      DO i = 1 , no_hub_dofs
        rotr_inrt_load(m_rotr+i) = hub_inrt_vec(i,1)
      ENDDO
      
      f_rotr_val = -acc_val*rotr_inrt_load
      
c	Newmark Beta without Ramp Contact
c	=================================
	CALL d_mult_sparse_mat_vec(rotr_G,rotr_G_map,v_rotr_prd,
     &	rotr_vec_1,no_rotr_dofs,no_rotr_dofs,no_rotr_dofs)
	CALL d_mult_sparse_mat_vec(rotr_K,rotr_K_map,w_rotr_prd,
     &	rotr_vec_2,no_rotr_dofs,no_rotr_dofs,no_rotr_dofs)
      CALL d_mult_sparse_mat_vec(C_rotr,C_rotr_map,w_rotr_prd,
     &	rotr_vec_3,no_disk_bcs,no_rotr_dofs,no_rotr_dofs)
      time_fac = -(2.0d0/((dt**2)*disk_nb_2))
      rotr_vec_3 = time_fac*rotr_vec_3

	y_rotr(1:no_rotr_dofs) = f_rotr_val - rotr_vec_1 - rotr_vec_2
	y_rotr(no_rotr_dofs+1:no_rotr_vrbl) = rotr_vec_3
	
c     Solution of system of equations
c     ===============================
	CALL LU_solve_sparse(L_rotr,U_rotr,L_rotr_map,U_rotr_map,
     &	    y_rotr(1:no_rotr_vrbl),x_rotr(1:no_rotr_vrbl),
     &      no_rotr_vrbl)

c	Compute current kinematic Value
c	===============================     
	a_rotr_cur = x_rotr(1:no_rotr_dofs)
 	v_rotr_cur = v_rotr_prd + (disk_nb_1*dt*a_rotr_cur)
	w_rotr_cur = w_rotr_prd + (0.5d0*disk_nb_2*(dt**2)*a_rotr_cur)

c	Constraint Force
c	================
	f_rotr_cnst(1:no_disk_bcs) = x_rotr(no_rotr_dofs+1:no_rotr_vrbl)
	
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
		
c		Ramp Location
c		=============
		ang_cur = ramp_ang  - (ang_vel*t_val)
		
c       Adjust for disk center movement
c       ===============================
        xtemp = ramp_rad*DCOS(ang_cur) - x0_disk_s
        ytemp = ramp_rad*DSIN(ang_cur) - y0_disk_s
        
        rad_temp = DSQRT((xtemp**2) + (ytemp**2))
        ang_temp = ATAN(ABS(ytemp/xtemp))
        
        IF ((xtemp .LT. 0) .AND. (ytemp .LT. 0)) THEN
            ang_temp = pi + ang_temp
        ELSEIF ((xtemp .GE. 0) .AND. (ytemp .LT. 0)) THEN
            ang_temp = 2*pi - ang_temp
        ELSEIF ((xtemp .LT. 0) .AND. (ytemp .GE. 0)) THEN
            ang_temp =   pi - ang_temp
        ENDIF
		
c		Position of ramp lowest point in disk (rotating) frame
c		======================================================
		CALL rad_loc(rad_temp,rad_grid,rad_gridpt,zeta)
		CALL ang_loc(ang_temp,ang_grid,ang_gridpt,eta)
		zeta_temp = zeta		
		    
c		DOF Coefficient for w_ramp displacement
c		=======================================
  		CALL ramp_cnstrnt(rad_gridpt,ang_gridpt,
     &		zeta_temp,eta,ramp_coef_rotr)
     
c       Rearrange using mapping
c       =======================
        C_rotr_r(1:no_disk_dofs) = 0.0d0
        DO i = 1, no_disk_dofs
            C_rotr_r(disk_intf_map(i,2)) = ramp_coef_rotr(i)
        ENDDO
        
c		Disk Displacement under the ramp
c		================================
		w_ramp = DOT_PRODUCT(C_rotr_r(1:no_disk_dofs),
     &		w_rotr_cur(1:no_disk_dofs))  
     
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
			w_ramp = DOT_PRODUCT(C_rotr_r(1:no_rotr_dofs),
     &		    w_rotr_prd(1:no_rotr_dofs))

			lmbd = (2.0d0/((dt**2)*disk_nb_2))*(w_cnst - w_ramp)  
     
c           Block Factorization
c           ===================           
            CALL L_solve_sparse(L_rotr,L_rotr_map,
     &			    C_rotr_r,U_12_rotr,no_rotr_vrbl)
			CALL U_T_solve_sparse(U_rotr,U_rotr_T_map,
     &			    C_rotr_r,L_21_rotr,no_rotr_vrbl)
			U_22 = DOT_PRODUCT(L_21_rotr(1:no_rotr_vrbl),
     &		        U_12_rotr(1:no_rotr_vrbl))
     
c           Solution of system of equations
c           ===============================

c           Forward : L solution
            CALL L_solve_sparse(L_rotr,L_rotr_map,y_rotr,
     &           z_rotr,no_rotr_vrbl)
            nu_val = lmbd - DOT_PRODUCT(L_21_rotr,z_rotr)
            
c           Ramp Contact Force
            ramp_cntf = -nu_val/U_22

c           Backward : U solution
            DO i = 1 , no_rotr_vrbl
                z_rotr(i) = z_rotr(i) - ramp_cntf*U_12_rotr(i)
            ENDDO
            CALL U_solve_sparse(U_rotr,U_rotr_map,z_rotr,
     &           x_rotr,no_rotr_vrbl)
           
c			Compute Current Kinematic Value
c			===============================
			a_rotr_cur = x_rotr(1:no_rotr_dofs)
 			v_rotr_cur = v_rotr_prd + (disk_nb_1*dt*a_rotr_cur)
			w_rotr_cur = w_rotr_prd + 
     &		          (0.5d0*disk_nb_2*(dt**2)*a_rotr_cur)
     
c	      Disk ID Constraint Forces
c	      =========================
	      f_rotr_cnst(1:no_disk_bcs) = 
     &	        x_rotr(no_rotr_dofs+1:no_rotr_vrbl)     

		ENDIF
           
      ENDIF
	
c     Compute disk displacement
c     =========================
      DO i = 1 , no_disk_dofs
        w_disk_cur(disk_intf_map(i,1)) = w_rotr_cur(i)
      ENDDO
      
c     Disk ID displacement
c     ====================
      DO i = 1, no_disk_intf_dofs
        DO j = 1 , 3
            disk_ID_disp(i,j) = w_rotr_cur(disk_ID_map(i,j+2))
        ENDDO
      ENDDO
 
c     FDB dispalcement
c     ================        
      IF (data_opt_supp.EQ.1) THEN
      
        w_FDB(1) = t_val
        pvt = no_rotr_dofs - no_fdb_dofs_hub

        DO i = 1 , no_fdb_dofs_hub 
            w_FDB(i+1) = w_rotr_cur(pvt+i)
        ENDDO
        
      ENDIF  

      RETURN
      END

c======================================================================
      SUBROUTINE init_housing
c======================================================================

c     Shared data
c     ===========
      USE shck_data
      USE disk_data
      USE motr_data
      USE rotr_data
      USE hsng_data
      USE base_data
      
      IMPLICIT NONE
      
c     Local variables
c     ===============
      INTEGER i,j,i_cnt
      INTEGER dof_no,hsng_dof
      
c     Temporary matrices/vectors
c     ==========================
      INTEGER, DIMENSION (:),   ALLOCATABLE :: rearrange_index
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: hsng_intf_map

      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: M_row
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: G_row
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: K_row
      
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: M_col
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: G_col
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: K_col

	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: hsng_inrt_map
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: hsng_inrt_vec
      
c     Stationary Housing structural matrices
c     ======================================
      ALLOCATE(hsng_M(no_hsng_dofs,no_hsng_dofs))
      ALLOCATE(hsng_G(no_hsng_dofs,no_hsng_dofs))
      ALLOCATE(hsng_K(no_hsng_dofs,no_hsng_dofs))

c     Initialize
c     ==========      
      hsng_M = 0.0d0
      hsng_G = 0.0d0
      hsng_K = 0.0d0
      
c     Housing coordinate list      
c     =======================
      ALLOCATE(hsng_cords(no_hsng_dofs,4))
      hsng_cords = 0.0d0
            
c     Read input ANSYS data
c     =====================      
      OPEN(84,ERR=999,FILE='hsng_M.dat',STATUS='OLD')
      OPEN(85,ERR=999,FILE='hsng_K.dat',STATUS='OLD')
      OPEN(86,ERR=999,FILE='hsng_cords.dat',STATUS='OLD')

      DO i=1,no_hsng_dofs
		READ(84,*) (hsng_M(i,j),j=1,no_hsng_dofs)
        READ(85,*) (hsng_K(i,j),j=1,no_hsng_dofs)
        READ(86,*) (hsng_cords(i,j),j=1,4)
	ENDDO
	
      CLOSE(84)
      CLOSE(85)
      CLOSE(86)
      
c     Proportional damping matrix
c     ===========================      
      hsng_G = (disk_M_damp*hsng_M) + (disk_K_damp*hsng_K)

c     Rearrange housing dofs for coupling
c     ===================================
      ALLOCATE(hsng_intf_map(no_hsng_dofs,2))
      ALLOCATE(rearrange_index(no_hsng_dofs))
      
c     Initialize
c     ==========      
      hsng_intf_map = 0
      rearrange_index = 0
      
c     FDB interface dofs
c     ==================      
      DO i = 1 , no_fdb_dofs_hsng
        rearrange_index(hsng_fdb_dofs(i)) = 1
      ENDDO
            
c     Base Interface
c     ==============
      IF (supp_opt .EQ. 3) THEN
        DO i = 1 , no_base_intf_dofs
            rearrange_index(base_intf_dofs(i,2)) = 2
        ENDDO
      ENDIF

c     Rearranging DOFs
c     ================
      
c     FDB interface dofs
      DO i = 1 , no_fdb_dofs_hsng
        hsng_intf_map(i,1) = hsng_fdb_dofs(i)
      ENDDO
      
      i_cnt = no_fdb_dofs_hsng
      
c     Non interface dofs
      DO i = 1, no_hsng_dofs
        IF (rearrange_index(i) .EQ. 0) THEN
            i_cnt = i_cnt + 1
            hsng_intf_map(i_cnt,1) = i
        ENDIF
      ENDDO
      
c     Base interface dofs      
      IF (supp_opt .EQ. 3) THEN
        DO i = 1, no_hsng_dofs
            IF (rearrange_index(i) .EQ. 2) THEN
                i_cnt = i_cnt + 1
                hsng_intf_map(i_cnt,1) = i
            ENDIF
        ENDDO
      ENDIF
      
c     Inverse map
c     ===========
      DO i = 1 , no_hsng_dofs
        hsng_intf_map(hsng_intf_map(i,1),2) = i
      ENDDO
      
c     Deallocate temprory matrix
c     ==========================      
      DEALLOCATE(rearrange_index)
      
c     Rearrange housing matrices
c     ==========================
      ALLOCATE(M_row(no_hsng_dofs,no_hsng_dofs))
      ALLOCATE(G_row(no_hsng_dofs,no_hsng_dofs))
      ALLOCATE(K_row(no_hsng_dofs,no_hsng_dofs))      
      
      M_row = 0.0d0
      G_row = 0.0d0
      K_row = 0.0d0
      
      ALLOCATE(M_col(no_hsng_dofs,no_hsng_dofs))
      ALLOCATE(G_col(no_hsng_dofs,no_hsng_dofs))
      ALLOCATE(K_col(no_hsng_dofs,no_hsng_dofs))
      
      M_col = 0.0d0
      G_col = 0.0d0
      K_col = 0.0d0
      
c     Rearrange rows
c     ==============
      DO i = 1, no_hsng_dofs
        M_row(i,1:no_hsng_dofs) = 
     &      hsng_M(hsng_intf_map(i,1),1:no_hsng_dofs)
        G_row(i,1:no_hsng_dofs) = 
     &      hsng_G(hsng_intf_map(i,1),1:no_hsng_dofs)
        K_row(i,1:no_hsng_dofs) = 
     &      hsng_K(hsng_intf_map(i,1),1:no_hsng_dofs)
      ENDDO
      
c     Rearrange columns
c     =================
      DO i = 1, no_hsng_dofs
        M_col(1:no_hsng_dofs,i) =
     &    M_row(1:no_hsng_dofs,hsng_intf_map(i,1))
        G_col(1:no_hsng_dofs,i) =
     &    G_row(1:no_hsng_dofs,hsng_intf_map(i,1))
        K_col(1:no_hsng_dofs,i) =
     &    K_row(1:no_hsng_dofs,hsng_intf_map(i,1))
      ENDDO      
      
c     Rearranged housing matrices
c     ===========================  
      hsng_M = M_col
      hsng_G = G_col
      hsng_K = K_col
      
c     Deallocate temporary matrices
c     =============================
      DEALLOCATE(M_row)
      DEALLOCATE(G_row)
      DEALLOCATE(K_row)
     
      DEALLOCATE(M_col)
      DEALLOCATE(G_col)
      DEALLOCATE(K_col)
      
c     Housing inertia load map
c     ========================
      ALLOCATE(hsng_inrt_map(no_hsng_dofs,3))
      hsng_inrt_map = 0.0d0
      
      DO i = 1 , no_hsng_dofs
        dof_no = INT(hsng_cords(i,4))
        hsng_dof =  hsng_intf_map(i,2)   
        hsng_inrt_map(hsng_dof,dof_no) = 1.0d0
      ENDDO

c     Housing DOFs accelaration
c     =========================
      ALLOCATE(hsng_inrt_acc(no_hsng_dofs,1))
      hsng_inrt_acc = 0.0d0
      hsng_inrt_acc = MATMUL(hsng_inrt_map,shk_dir)
      DEALLOCATE(hsng_inrt_map)

c     Housing inertia load vector
c     ===========================      
      IF (supp_opt .EQ. 2) THEN
        ALLOCATE(hsng_inrt_vec(no_hsng_dofs,1))
        hsng_inrt_vec = 0.0d0
        hsng_inrt_vec = MATMUL(hsng_M,hsng_inrt_acc)
        DEALLOCATE(hsng_inrt_acc)
      
        ALLOCATE(hsng_inrt_load(no_hsng_dofs))
        hsng_inrt_load = 0.0d0
        DO i = 1, no_hsng_dofs
            hsng_inrt_load(i) = hsng_inrt_vec(i,1)
        ENDDO
        DEALLOCATE(hsng_inrt_vec)
      ENDIF
            
      RETURN
      
999	WRITE(*,*) 'Trouble in opening housing files'
	STOP      
      
      END

c======================================================================
      SUBROUTINE couple_rotr_hsng
c======================================================================

c     Shared data
c     ===========
      USE motr_data
      USE rotr_data
      USE hsng_data
      
      IMPLICIT NONE
      
c     Local variables
c     ===============
      INTEGER i,j,row_no,col_no
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: temp_mat
      
c     DOF data
c     ========     
      no_motr_dofs = no_rotr_dofs + no_hsng_dofs
      no_motr_vrbl = no_motr_dofs + no_disk_bcs
      no_fdb_dofs_motr = no_fdb_dofs_hub+no_fdb_dofs_hsng

c     Initialize coupled matrices
c     ===========================
      ALLOCATE(motr_M(no_motr_dofs,no_motr_dofs))
      ALLOCATE(motr_G(no_motr_dofs,no_motr_dofs))
      ALLOCATE(motr_K(no_motr_dofs,no_motr_dofs))
      motr_M = 0.0d0
      motr_G = 0.0d0
      motr_K = 0.0d0
      
c     Coupling
c     ========

c     Add Rotor matrices
      motr_M(1:no_rotr_dofs,1:no_rotr_dofs) = 
     &  rotr_M(1:no_rotr_dofs,1:no_rotr_dofs)
      motr_G(1:no_rotr_dofs,1:no_rotr_dofs) = 
     &  rotr_G(1:no_rotr_dofs,1:no_rotr_dofs)
      motr_K(1:no_rotr_dofs,1:no_rotr_dofs) = 
     &  rotr_K(1:no_rotr_dofs,1:no_rotr_dofs)
     
c     Add Housing matrices
      motr_M(no_rotr_dofs+1:no_rotr_dofs+no_hsng_dofs,
     &  no_rotr_dofs+1:no_rotr_dofs+no_hsng_dofs) = 
     &  hsng_M(1:no_hsng_dofs,1:no_hsng_dofs) 
      motr_G(no_rotr_dofs+1:no_rotr_dofs+no_hsng_dofs,
     &  no_rotr_dofs+1:no_rotr_dofs+no_hsng_dofs) = 
     &  hsng_G(1:no_hsng_dofs,1:no_hsng_dofs)
      motr_K(no_rotr_dofs+1:no_rotr_dofs+no_hsng_dofs,
     &  no_rotr_dofs+1:no_rotr_dofs+no_hsng_dofs) = 
     &  hsng_K(1:no_hsng_dofs,1:no_hsng_dofs)

c     Deallocate system matrices
c     ==========================
      DEALLOCATE(rotr_M)
      DEALLOCATE(rotr_G)
      DEALLOCATE(rotr_K)
      
      DEALLOCATE(hsng_M)
      DEALLOCATE(hsng_G)
      DEALLOCATE(hsng_K)
      
c     Sparse Map for stationary housing matrices
c     ==========================================
      ALLOCATE(motr_M_map(no_motr_dofs,no_motr_dofs+1))
	CALL sparse_mat_map(motr_M,motr_M_map,
     &	no_motr_dofs,no_motr_dofs)
     
      ALLOCATE(temp_mat(no_motr_dofs,no_motr_dofs))

c     Damping matrix get affected by FDB coupling
      temp_mat = motr_G
      DO i = 1 , no_fdb_dofs_motr
        row_no = no_rotr_dofs - no_fdb_dofs_hub + i
        DO j = 1 , no_fdb_dofs_motr
            col_no = no_rotr_dofs - no_fdb_dofs_hub + j
            temp_mat(row_no,col_no) = 1.0d0
        ENDDO
      ENDDO    
      ALLOCATE(motr_G_map(no_motr_dofs,no_motr_dofs+1))
	CALL sparse_mat_map(temp_mat,motr_G_map,
     &	no_motr_dofs,no_motr_dofs)

c     Stiffness matrix get affected by FDB coupling
      temp_mat = motr_K
      DO i = 1 , no_fdb_dofs_motr
        row_no = no_rotr_dofs - no_fdb_dofs_hub + i
        DO j = 1 , no_fdb_dofs_motr
            col_no = no_rotr_dofs - no_fdb_dofs_hub + j
            temp_mat(row_no,col_no) = 1.0d0
        ENDDO
      ENDDO 
	ALLOCATE(motr_K_map(no_motr_dofs,no_motr_dofs+1))
	CALL sparse_mat_map(temp_mat,motr_K_map,
     &	no_motr_dofs,no_motr_dofs)
     
      DEALLOCATE(temp_mat)    
      
c     Boundary Condition : Clamped disk center
c     ========================================

c     Initialize
      ALLOCATE(C_motr(no_disk_bcs,no_motr_dofs))
      C_motr = 0.0d0
      
c     C_motr is same as C_rotr
      C_motr(1:no_disk_bcs,1:no_rotr_dofs) = 
     &  C_rotr(1:no_disk_bcs,1:no_rotr_dofs)
     
c     Deallocate
      DEALLOCATE(C_rotr)       
      
c     Sparse map     
	ALLOCATE(C_motr_map(no_disk_bcs,no_motr_dofs+1))
	CALL sparse_mat_map(C_motr,C_motr_map,
     &	no_disk_bcs,no_motr_dofs)
     
c     Motor inertia load vector
c     =========================
      ALLOCATE(motr_inrt_load(no_motr_dofs))
      motr_inrt_load = 0.0d0

c     -------------------------------------------------------------      
c     Note : Inertia matrix are uncoupled since FDB has no inertia
c     -------------------------------------------------------------

c     Add disk part of rotor inertia load vector
c     ==========================================
      motr_inrt_load(1:m_rotr) = rotr_inrt_load(1:m_rotr)
     
c     Add Housing vector
c     ==================
      motr_inrt_load(no_rotr_dofs+1:no_rotr_dofs+no_hsng_dofs) = 
     &  hsng_inrt_load(1:no_hsng_dofs)
     
c     Deallocate
c     ==========
      DEALLOCATE(rotr_inrt_load)
      DEALLOCATE(hsng_inrt_load)

      RETURN
      END
c======================================================================
      SUBROUTINE newmark_motor
c======================================================================

c     Shared data
c     ===========
      USE sldr_dynm
      USE disk_data
      USE disk_aray
      USE motr_data
      USE rotr_data
      USE hsng_data
      USE motr_soln
      
      IMPLICIT NONE
      
c     Local Variables
c     ===============
      INTEGER i,j,m,n,p 
                  
c     Block Newmark Beta Matrices
c     ===========================
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_mat
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: C_motr_T
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_motr_T
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: temp_mat
	
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_motr
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_motr_11
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_motr_12  
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_motr_13
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_motr_21
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_motr_23
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_motr_31
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_motr_32
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_motr_33
      
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_motr_11
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_motr_21
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_motr_31
      
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_motr_11
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_motr_12
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_motr_13

c     Newmark Beta Formulation
c	========================
	ALLOCATE(A_mat(no_motr_dofs,no_motr_dofs))
	A_mat = 0.0d0   ! Intialize
 	A_mat = motr_M + (disk_nb_1*dt*motr_G) +
     &	(0.5d0*disk_nb_2*(dt**2)*motr_K)
     
c     C_motr_T for constrained formulation
c     =====================================	
	ALLOCATE(C_motr_T(no_motr_dofs,no_disk_bcs))
	C_motr_T = 0.0d0
	C_motr_T(1:no_motr_dofs,1:no_disk_bcs) = 
     &	TRANSPOSE(C_motr(1:no_disk_bcs,1:no_motr_dofs))
     
c	Constraint BC at ID (Lagrange Multiplier Method)
c	================================================
	ALLOCATE(A_motr(no_motr_vrbl,no_motr_vrbl))
	A_motr = 0.0d0

	A_motr(1:no_motr_dofs,1:no_motr_dofs) = 
     &	A_mat(1:no_motr_dofs,1:no_motr_dofs)
	
 	A_motr(no_motr_dofs+1:no_motr_vrbl,1:no_motr_dofs) = 
     &	C_motr(1:no_disk_bcs,1:no_motr_dofs)

 	A_motr(1:no_motr_dofs,no_motr_dofs+1:no_motr_vrbl) = 
     &	C_motr_T(1:no_motr_dofs,1:no_disk_bcs)
     
      DEALLOCATE(A_mat)
      DEALLOCATE(C_motr_T)
      
c     Block factorization : Time Independent
c     ======================================
      m = no_rotr_dofs - no_fdb_dofs_hub
      n = no_fdb_dofs_motr
      p = no_motr_vrbl - m - n
      
c     A11 = L11*U11
      ALLOCATE(A_motr_11(m,m))
      ALLOCATE(L_motr_11(m,m))
      ALLOCATE(U_motr_11(m,m))
      
      A_motr_11 = 0.0d0
      L_motr_11 = 0.0d0
      U_motr_11 = 0.0d0
      
      A_motr_11(1:m,1:m) = A_motr(1:m,1:m)
      CALL LU_fact(A_motr_11,L_motr_11,U_motr_11,m)
      DEALLOCATE(A_motr_11)
      
c     A12 = L11*U12
      ALLOCATE(A_motr_12(m,n))
      ALLOCATE(U_motr_12(m,n))
      
      A_motr_12 = 0.0d0
      U_motr_12 = 0.0d0
      
      A_motr_12(1:m,1:n) = A_motr(1:m,m+1:m+n)
      CALL L_solve_mult(L_motr_11,A_motr_12,U_motr_12,n,m)
      DEALLOCATE(A_motr_12)
      
c     A13 = L11*U13
      ALLOCATE(A_motr_13(m,p))
      ALLOCATE(U_motr_13(m,p))
      
      A_motr_13 = 0.0d0
      U_motr_13 = 0.0d0
      
      A_motr_13(1:m,1:p) = A_motr(1:m,m+n+1:m+n+p)
      CALL L_solve_mult(L_motr_11,A_motr_13,U_motr_13,p,m)
      DEALLOCATE(A_motr_13)
      
c     A21 = L21*U11 => A21_T = U11_T*L_21*T
      ALLOCATE(A_motr_21(n,m))
      ALLOCATE(L_motr_21(n,m))
      
      A_motr_21 = 0.0d0
      L_motr_21 = 0.0d0

      A_motr_21(1:n,1:m) = A_motr(m+1:m+n,1:m)
      CALL U_T_solve_mult(U_motr_11,A_motr_21,L_motr_21,n,m)
      DEALLOCATE(A_motr_21)
      
c     A31 = L31*U11 => A31_T = U11_T*L_31*T
      ALLOCATE(A_motr_31(p,m))
      ALLOCATE(L_motr_31(p,m))
      
      A_motr_31 = 0.0d0
      L_motr_31 = 0.0d0

      A_motr_31(1:p,1:m) = A_motr(m+n+1:m+n+p,1:m)
      CALL U_T_solve_mult(U_motr_11,A_motr_31,L_motr_31,p,m)
      DEALLOCATE(A_motr_31)      
      
c     A22 = L21*U12 + L22U22 : L21*U12 = A_motr_22_m      
      ALLOCATE(A_motr_22_m(n,n))
      A_motr_22_m = MATMUL(L_motr_21,U_motr_12)
      
c     A23 = L21*U13 + L22*U23 => A23 - L21*U13 = A_motr_23_m
      ALLOCATE(A_motr_23(n,p))
      ALLOCATE(A_motr_23_m(n,p))
      
      A_motr_23   = 0.0d0
      A_motr_23_m = 0.0d0
      
      A_motr_23(1:n,1:p) = A_motr(m+1:m+n,m+n+1:m+n+p)
      A_motr_23_m = A_motr_23 - MATMUL(L_motr_21,U_motr_13)
      DEALLOCATE(A_motr_23)
      
c     A32 = L31*U12 + L32*U22 => A32 - L31*U12 = A_motr_32_m
      ALLOCATE(A_motr_32(p,n))
      ALLOCATE(A_motr_32_m(p,n))
      
      A_motr_32   = 0.0d0
      A_motr_32_m = 0.0d0
      
      A_motr_32(1:p,1:n) = A_motr(m+n+1:m+n+p,m+1:m+n)
      A_motr_32_m = A_motr_32 - MATMUL(L_motr_31,U_motr_12)
      DEALLOCATE(A_motr_32)

c     A33 = L31*U13 + L32*U23 + L33*U33 => A33 - L31*U13 = A_motr_33_m
      ALLOCATE(A_motr_33(p,p))
      ALLOCATE(A_motr_33_m(p,p))
      
      A_motr_33   = 0.0d0
      A_motr_33_m = 0.0d0
      
      A_motr_33(1:p,1:p) = A_motr(m+n+1:m+n+p,m+n+1:m+n+p)
      A_motr_33_m = A_motr_33 - MATMUL(L_motr_31,U_motr_13)
      DEALLOCATE(A_motr_33)
      
c     Factorization done      
      DEALLOCATE(A_motr)
      
c     Assemble L-U Blocks
c     ===================
      ALLOCATE(L_motr(no_motr_vrbl,no_motr_vrbl))
      ALLOCATE(U_motr(no_motr_vrbl,no_motr_vrbl))
      
      L_motr = 0.0d0
      U_motr = 0.0d0
     
      L_motr(1:m,1:m) = L_motr_11(1:m,1:m)
      L_motr(m+1:m+n,1:m) = L_motr_21(1:n,1:m)
      L_motr(m+n+1:m+n+p,1:m) = L_motr_31(1:p,1:m)
      
      U_motr(1:m,1:m) = U_motr_11(1:m,1:m)
      U_motr(1:m,m+1:m+n) = U_motr_12(1:m,1:n)
      U_motr(1:m,m+n+1:m+n+p) = U_motr_13(1:m,1:p)
      
c     Sparse Map
c     ==========
      ALLOCATE(temp_mat(no_motr_vrbl,no_motr_vrbl))
      
c     L map      
      temp_mat = 0.0d0
      temp_mat = L_motr
      
c     L22
      DO i = m+1 , m+n
        DO j = m+1 , i
            temp_mat(i,j) = 1.0d0
        ENDDO
      ENDDO

c     L32
      DO i = m+n+1 , m+n+p
        DO j = m+1 , m+n
            temp_mat(i,j) = 1.0d0
        ENDDO
      ENDDO
      
c     L33 
      DO i = m+n+1 , m+n+p
        DO j = m+n+1 , i
            temp_mat(i,j) = 1.0d0
        ENDDO
      ENDDO
      
      ALLOCATE(L_motr_map(no_motr_vrbl,no_motr_vrbl+1))
      CALL sparse_mat_map(temp_mat,L_motr_map,
     &	no_motr_vrbl,no_motr_vrbl)

c     U map
      temp_mat = 0.0d0
      temp_mat = U_motr
      
c     U22
      DO i = m+1 , m+n
        DO j = i , m+n
            temp_mat(i,j) = 1.0d0
        ENDDO
      ENDDO  
      
c     U23
      DO i = m+1 , m+n
        DO j = m+n+1 , m+n+p
            temp_mat(i,j) = 1.0d0
        ENDDO
      ENDDO
      
c     U33
      DO i = m+n+1 , m+n+p
        DO j = i , m+n+p
            temp_mat(i,j) = 1.0d0
        ENDDO
      ENDDO
      
      ALLOCATE(U_motr_map(no_motr_vrbl,no_motr_vrbl+1))
      CALL sparse_mat_map(temp_mat,U_motr_map,
     &	no_motr_vrbl,no_motr_vrbl)
     
c     U_T map     
      ALLOCATE(U_motr_T(no_motr_vrbl,no_motr_vrbl))
      U_motr_T = TRANSPOSE(temp_mat)
            
      ALLOCATE(U_motr_T_map(no_motr_vrbl,no_motr_vrbl+1))
      CALL sparse_mat_map(U_motr_T,U_motr_T_map,
     &	no_motr_vrbl,no_motr_vrbl)
      DEALLOCATE(U_motr_T)

      DEALLOCATE(temp_mat)
      
c     FDB Interface Block Matrices
c     ============================
      ALLOCATE(intf_M(n,n))
      ALLOCATE(intf_G(n,n))
      ALLOCATE(intf_K(n,n))
      
      intf_M = 0.0d0
      intf_G = 0.0d0
      intf_K = 0.0d0
      
      intf_M(1:n,1:n) = motr_M(m+1:m+n,m+1:m+n)
      intf_G(1:n,1:n) = motr_G(m+1:m+n,m+1:m+n)
      intf_K(1:n,1:n) = motr_K(m+1:m+n,m+1:m+n)   
      
c     Initialize matrices used in motor solver
c     ========================================

c     Initialize circumferential interpolation
      ALLOCATE(crcm_intp_mat(no_crcm_node_hub,4))
      crcm_intp_mat = 0.0d0
      
c     Hub current angle list
      ALLOCATE(hub_cur_ang(no_crcm_node_hub))
      hub_cur_ang = 0.0d0
      
c     Temporary Newmark Beta Block Matrices
      ALLOCATE(temp_G_m(n,n))
      ALLOCATE(temp_K_m(n,n))
      temp_G_m = 0.0d0
      temp_K_m = 0.0d0
      
c     Newmark Beta Block Matrices
      ALLOCATE(A_motr_22(n,n))
      ALLOCATE(L_motr_22(n,n))
      ALLOCATE(U_motr_22(n,n))
      
      ALLOCATE(U_motr_23(n,p))
      ALLOCATE(L_motr_32(p,n))
      
      ALLOCATE(L_motr_33(p,p))
      ALLOCATE(U_motr_33(p,p))
      
      A_motr_22 = 0.0d0
      L_motr_22 = 0.0d0
      U_motr_22 = 0.0d0
      
      L_motr_32 = 0.0d0
      U_motr_23 = 0.0d0
      
      L_motr_33 = 0.0d0
      U_motr_33 = 0.0d0
      
      ALLOCATE(blck_motr_1(n,n))
      ALLOCATE(blck_motr_2(p,p))
      blck_motr_1 = 0.0d0
      blck_motr_2 = 0.0d0

c     Block Matrices size      
      m_motr = m
      n_motr = n
      p_motr = p
      
c	Initialize Dynamic Variables
      ALLOCATE(w_disk_cur(no_disk_dofs))
	ALLOCATE(w_motr_cur(no_motr_dofs))
	ALLOCATE(v_motr_cur(no_motr_dofs))
	ALLOCATE(a_motr_cur(no_motr_dofs))
	ALLOCATE(w_motr_prv(no_motr_dofs))
	ALLOCATE(v_motr_prv(no_motr_dofs))
	ALLOCATE(a_motr_prv(no_motr_dofs))
	ALLOCATE(f_motr_cnst(no_disk_bcs))
	
c     Motor solution variables
      ALLOCATE(w_motr_prd(no_motr_dofs))
	ALLOCATE(v_motr_prd(no_motr_dofs))
	ALLOCATE(f_motr_val(no_motr_dofs))

      ALLOCATE(motr_vec_1(no_motr_dofs))
	ALLOCATE(motr_vec_2(no_motr_dofs))
	ALLOCATE(motr_vec_3(no_disk_bcs))
	
	ALLOCATE(x_motr(no_motr_vrbl))
	ALLOCATE(y_motr(no_motr_vrbl))
	ALLOCATE(z_motr(no_motr_vrbl))
	
	ALLOCATE(ramp_coef_motr(no_disk_dofs))
	ALLOCATE(C_motr_r(no_motr_vrbl))
	ALLOCATE(L_21_motr(no_motr_vrbl))
	ALLOCATE(U_12_motr(no_motr_vrbl))	
	
	C_motr_r = 0.0d0  ! Initialize

      RETURN    
      END      

c======================================================================
      SUBROUTINE solve_motor
c======================================================================

c     Shared data
c     ===========
      USE sldr_dynm
      USE shck_data
      USE disk_data
      USE disk_aray
      USE hsng_data
      USE rotr_data
      USE motr_data
      USE motr_soln
      
      IMPLICIT NONE
      
c     Parameters
c     ==========
      DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979

c     Local variables
c     ===============
      INTEGER i,j,i_brng,i_node,brng_row,pvt
      INTEGER FDB_node_map(9),hsng_node_no
      DOUBLE PRECISION ang_pos,hub_ang
      DOUBLE PRECISION K_local(3,3),C_local(3,3)
      DOUBLE PRECISION K_nodal(9,9),C_nodal(9,9)
      
      INTEGER rad_gridpt(2),ang_gridpt(2)
      DOUBLE PRECISION time_fac
      DOUBLE PRECISION zeta,eta,zeta_temp
      DOUBLE PRECISION w_cnst,w_ramp,U_22,lmbd,nu_val
	DOUBLE PRECISION ang_cur,ang_vel,eps_ang
	DOUBLE PRECISION xtemp,ytemp,rad_temp,ang_temp

c	Angular Speed
	ang_vel = disk_RPM*(2.0d0*pi/60.0d0)
	
c     Current angular postion in HDD frame
      ang_pos = ang_vel*t_val

c     Circumferential interpolation matrix
      CALL crcm_fdb_intp(ang_pos)

c     Hub current angle positions
      DO i = 1 , no_crcm_node_hub
        hub_cur_ang(i) = ang_pos + crcm_ang_hub(i)
      ENDDO
      
c     Initialize FDB interface dynamic coefficient matrices
      temp_G_m = intf_G
      temp_K_m = intf_K
      
c     Add FDB Stiffness and damping
      DO i_brng = 1 , 3 ! Bearing type : 3
      
        brng_row  = 3*(i_brng-1)
        
c       Reset    
        C_local = 0.0d0    
        K_local = 0.0d0
        
c       Local Dynanic Coefficients
        DO i = 1 , 3
            DO j = 1 , 3
                C_local(i,j) = C_fdb(brng_row+i,j)/
     &                         DBLE(no_crcm_node_hub)
                K_local(i,j) = K_fdb(brng_row+i,j)/
     &                         DBLE(no_crcm_node_hub)
            ENDDO
        ENDDO
        
c       Add nodal coefficents
        DO i_node = 1 , no_crcm_node_hub
        
c           Nodal Dynamic Coefficients
            CALL nodal_FDB_dyn_coef(C_local,K_local,C_nodal,K_nodal,
     &           ang_pos,crcm_intp_mat(i_node,1),
     &           crcm_intp_mat(i_node,2))
        
c           Find location of inteface FDB nodes on hub
            pvt = 3*no_crcm_node_hub*(i_brng-1) + 3*(i_node-1)
            DO j = 1 , 3
                FDB_node_map(j) = pvt + j
            ENDDO
            
c           Find location of interface FDB nodes on housing
            DO i = 1 , 2
                hsng_node_no = INT(crcm_intp_mat(i_node,i+2))
                pvt = (3*3*no_crcm_node_hub) + (3*no_crcm_node_hsng*
     &                (i_brng-1)) + 3*(hsng_node_no-1)
                DO j = 1 , 3
                    FDB_node_map(3*(i-1)+j+3) = pvt + j
                ENDDO
            ENDDO
           
            DO i = 1 , 9
                DO j = 1 , 9
                
                    temp_G_m(FDB_node_map(i),FDB_node_map(j)) = 
     &                  temp_G_m(FDB_node_map(i),FDB_node_map(j)) +
     &                  C_nodal(i,j)
     
                    temp_K_m(FDB_node_map(i),FDB_node_map(j)) = 
     &                  temp_K_m(FDB_node_map(i),FDB_node_map(j)) +
     &                  K_nodal(i,j)
     
                ENDDO
            ENDDO
            
        ENDDO

      ENDDO

c     Block Factorization
c     ===================
      
c     Newmark Beta Block Matrix : A22
      A_motr_22 = intf_M + (disk_nb_1*dt*temp_G_m) +
     &	(0.5d0*disk_nb_2*(dt**2)*temp_K_m)
     
c     A22 - L21*U12 = A_motr_22_m = L22*U22
      blck_motr_1 = A_motr_22 - A_motr_22_m
      CALL LU_fact(blck_motr_1,L_motr_22,U_motr_22,n_motr)
      
c     A23 - L21*U13 = A_motr_23_m = L22*U23
      CALL L_solve_mult(L_motr_22,A_motr_23_m,U_motr_23,
     &  p_motr,n_motr)

c     A32 - L31*U12 = A_motr_32_m = L32*U22 => U22'*L32' = A_motr_32_m'
      CALL U_T_solve_mult(U_motr_22,A_motr_32_m,L_motr_32,
     &  p_motr,n_motr)

c     A33 - L31*U13 = A_motr_33_m = L32*U23 + L33*U33
      blck_motr_2 = A_motr_33_m - MATMUL(L_motr_32,U_motr_23)
      CALL LU_fact(blck_motr_2,L_motr_33,U_motr_33,p_motr)
      
c     Assemble L-U Matrices
c     =====================

      L_motr(m_motr+1:m_motr+n_motr,
     &  m_motr+1:m_motr+n_motr) = 
     &  L_motr_22(1:n_motr,1:n_motr)                    ! L22

      L_motr(m_motr+n_motr+1:m_motr+n_motr+p_motr,
     &  m_motr+1:m_motr+n_motr) = 
     &  L_motr_32(1:p_motr,1:n_motr)                    ! L32
     
      L_motr(m_motr+n_motr+1:m_motr+n_motr+p_motr,
     &  m_motr+n_motr+1:m_motr+n_motr+p_motr) = 
     &  L_motr_33(1:p_motr,1:p_motr)                    ! L33
     
      U_motr(m_motr+1:m_motr+n_motr,
     &  m_motr+1:m_motr+n_motr) = 
     &  U_motr_22(1:n_motr,1:n_motr)                    ! U22
     
      U_motr(m_motr+1:m_motr+n_motr,
     &  m_motr+n_motr+1:m_motr+n_motr+p_motr) = 
     &  U_motr_23(1:n_motr,1:p_motr)                    ! U23

      U_motr(m_motr+n_motr+1:m_motr+n_motr+p_motr,
     &  m_motr+n_motr+1:m_motr+n_motr+p_motr) = 
     &  U_motr_33(1:p_motr,1:p_motr)                    ! U33
     
c     Update Motor damping and stiffness matrices
c     ===========================================
      motr_G(m_motr+1:m_motr+n_motr,
     &  m_motr+1:m_motr+n_motr) = 
     &  temp_G_m(1:n_motr,1:n_motr)
     
      motr_K(m_motr+1:m_motr+n_motr,
     &  m_motr+1:m_motr+n_motr) = 
     &  temp_K_m(1:n_motr,1:n_motr)

c	Predicted displacement
c	======================
	w_motr_prd = w_motr_prv + (dt*v_motr_prv) + 
     &	(0.5d0*(1.0d0-disk_nb_2)*(dt**2)*a_motr_prv)     
     
c	Predicted velocity
c	==================
	v_motr_prd = v_motr_prv + ((1.0d0-disk_nb_1)*dt*a_motr_prv)

c	Inertia load vector
c	===================

      CALL calc_hub_inrt(ang_vel,t_val)

c     Add hub part of rotor inertia vector      
      DO i = 1 , no_hub_dofs
        motr_inrt_load(m_rotr+i) = hub_inrt_vec(i,1)
      ENDDO
     
      f_motr_val = -acc_val*motr_inrt_load
      
c	Newmark Beta without Ramp Contact
c	=================================
	CALL d_mult_sparse_mat_vec(motr_G,motr_G_map,v_motr_prd,
     &	motr_vec_1,no_motr_dofs,no_motr_dofs,no_motr_dofs)
	CALL d_mult_sparse_mat_vec(motr_K,motr_K_map,w_motr_prd,
     &	motr_vec_2,no_motr_dofs,no_motr_dofs,no_motr_dofs)
      CALL d_mult_sparse_mat_vec(C_motr,C_motr_map,w_motr_prd,
     &	motr_vec_3,no_disk_bcs,no_motr_dofs,no_motr_dofs)
      time_fac = -(2.0d0/((dt**2)*disk_nb_2))
      motr_vec_3 = time_fac*motr_vec_3

	y_motr(1:no_motr_dofs) = f_motr_val - motr_vec_1 - motr_vec_2
	y_motr(no_motr_dofs+1:no_motr_vrbl) = motr_vec_3
	
c     Solution of system of equations
c     ===============================
	CALL LU_solve_sparse(L_motr,U_motr,L_motr_map,U_motr_map,
     &	    y_motr(1:no_motr_vrbl),x_motr(1:no_motr_vrbl),
     &      no_motr_vrbl)
     
c	Compute current kinematic Value
c	===============================     
	a_motr_cur = x_motr(1:no_motr_dofs)
 	v_motr_cur = v_motr_prd + (disk_nb_1*dt*a_motr_cur)
	w_motr_cur = w_motr_prd + (0.5d0*disk_nb_2*(dt**2)*a_motr_cur)
	
c	Constraint Force
c	================
	f_motr_cnst(1:no_disk_bcs) = x_motr(no_motr_dofs+1:no_motr_vrbl)

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
	  
c		Ramp Location
c		=============
		ang_cur = ramp_ang  - (ang_vel*t_val)
		
c       Adjust for disk center movement
c       ===============================
        xtemp = ramp_rad*DCOS(ang_cur) - x0_disk_s
        ytemp = ramp_rad*DSIN(ang_cur) - y0_disk_s
        
        rad_temp = DSQRT((xtemp**2) + (ytemp**2))
        ang_temp = ATAN(ABS(ytemp/xtemp))
        
        IF ((xtemp .LT. 0) .AND. (ytemp .LT. 0)) THEN
            ang_temp = pi + ang_temp
        ELSEIF ((xtemp .GE. 0) .AND. (ytemp .LT. 0)) THEN
            ang_temp = 2*pi - ang_temp
        ELSEIF ((xtemp .LT. 0) .AND. (ytemp .GE. 0)) THEN
            ang_temp =   pi - ang_temp
        ENDIF
        
c		Position of ramp lowest point in disk (rotating) frame
c		======================================================
		CALL rad_loc(rad_temp,rad_grid,rad_gridpt,zeta)
		CALL ang_loc(ang_temp,ang_grid,ang_gridpt,eta)
		zeta_temp = zeta
		
c		DOF Coefficient for w_ramp displacement
c		=======================================
  		CALL ramp_cnstrnt(rad_gridpt,ang_gridpt,
     &		zeta_temp,eta,ramp_coef_motr)
     
c       Rearrange using mapping
c       =======================
        C_motr_r(1:no_disk_dofs) = 0.0d0
        DO i = 1, no_disk_dofs
            C_motr_r(disk_intf_map(i,2)) = ramp_coef_motr(i)
        ENDDO
        
c		Disk Displacement under the ramp
c		================================
		w_ramp = DOT_PRODUCT(C_motr_r(1:no_disk_dofs),
     &		w_motr_cur(1:no_disk_dofs))  

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
			w_ramp = DOT_PRODUCT(C_motr_r(1:no_disk_dofs),
     &		    w_motr_prd(1:no_disk_dofs))

			lmbd = (2.0d0/((dt**2)*disk_nb_2))*(w_cnst - w_ramp) 

c           Block Factorization
c           ===================           
            CALL L_solve_sparse(L_motr,L_motr_map,
     &			    C_motr_r,U_12_motr,no_motr_vrbl)
			CALL U_T_solve_sparse(U_motr,U_motr_T_map,
     &			    C_motr_r,L_21_motr,no_motr_vrbl)
			U_22 = DOT_PRODUCT(L_21_motr(1:no_motr_vrbl),
     &		        U_12_motr(1:no_motr_vrbl))
     
c           Solution of system of equations
c           ===============================

c           Forward : L solution
            CALL L_solve_sparse(L_motr,L_motr_map,y_motr,
     &           z_motr,no_motr_vrbl)
            nu_val = lmbd - DOT_PRODUCT(L_21_motr,z_motr)
            
c           Ramp Contact Force
            ramp_cntf = -nu_val/U_22

c           Backward : U solution
            DO i = 1 , no_motr_vrbl
                z_motr(i) = z_motr(i) - ramp_cntf*U_12_motr(i)
            ENDDO
            CALL U_solve_sparse(U_motr,U_motr_map,z_motr,
     &           x_motr,no_motr_vrbl)
     
c			Compute Current Kinematic Value
c			===============================
			a_motr_cur = x_motr(1:no_motr_dofs)
 			v_motr_cur = v_motr_prd + (disk_nb_1*dt*a_motr_cur)
			w_motr_cur = w_motr_prd + 
     &		          (0.5d0*disk_nb_2*(dt**2)*a_motr_cur)   
     
c	      Disk ID Constraint Forces
c	      =========================
	      f_motr_cnst(1:no_disk_bcs) = 
     &	        x_motr(no_motr_dofs+1:no_motr_vrbl)        

        ENDIF
 	
 	ENDIF
 	
c     Compute disk displacement
c     =========================
      DO i = 1 , no_disk_dofs
        w_disk_cur(disk_intf_map(i,1)) = w_motr_cur(i)
      ENDDO 
      
c     Disk ID displacement
c     ====================
      DO i = 1, no_disk_intf_dofs
        DO j = 1 , 3
            disk_ID_disp(i,j) = w_motr_cur(disk_ID_map(i,j+2))
        ENDDO
      ENDDO
      
c     FDB dispalcement
c     ================        
      IF (data_opt_supp.EQ.1) THEN
          
        w_FDB(1) = t_val
          
        pvt = no_rotr_dofs - no_fdb_dofs_hub

        DO i = 1 , no_fdb_dofs_hub 
            w_FDB(i+1) = w_motr_cur(pvt+i)
        ENDDO  
      
        DO i = 1 , no_fdb_dofs_hsng
            w_FDB(i+no_fdb_dofs_hub+1) = 
     &          w_motr_cur(no_rotr_dofs+i)
        ENDDO 
      
      ENDIF  
      
      RETURN    
      END

c======================================================================
      SUBROUTINE init_base  
c======================================================================

c     Shared Data
c     ===========
      USE shck_data
      USE disk_data
      USE base_data

      IMPLICIT NONE
      
c     Local Variable
c     ==============
      INTEGER i,j,i_cnt
      INTEGER dof_no,base_dof

c     Temporary matrices/vectors
c     ==========================
      INTEGER, DIMENSION (:),   ALLOCATABLE :: rearrange_index
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: base_intf_map
      
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: M_row
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: G_row
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: K_row
      
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: M_col
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: G_col
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: K_col

	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: base_inrt_map
      
c     Base plate structural matrices
c     ==============================
      ALLOCATE(base_M(no_base_dofs,no_base_dofs))
      ALLOCATE(base_G(no_base_dofs,no_base_dofs))
      ALLOCATE(base_K(no_base_dofs,no_base_dofs))
      
c     Intiailize
c     ==========
      base_M = 0.0d0
      base_G = 0.0d0
      base_K = 0.0d0      
      
c     Base coordinate list      
c     ===================
      ALLOCATE(base_cords(no_base_dofs,4))
      base_cords = 0.0d0
      
c     Read input ANSYS data
c     =====================      
      OPEN(87,ERR=999,FILE='base_M.dat',STATUS='OLD')
      OPEN(88,ERR=999,FILE='base_K.dat',STATUS='OLD')
      OPEN(89,ERR=999,FILE='base_cords.dat',STATUS='OLD')

      DO i=1,no_base_dofs
		READ(87,*) (base_M(i,j),j=1,no_base_dofs)
        READ(88,*) (base_K(i,j),j=1,no_base_dofs)
        READ(89,*) (base_cords(i,j),j=1,4)
	ENDDO
	
      CLOSE(87)
      CLOSE(88)
      CLOSE(89)
      
c     Proportional damping matrix
c     ===========================      
      base_G = (disk_M_damp*base_M) + (disk_K_damp*base_K)
      
c     Rearrange base dofs for coupling
c     ================================
      ALLOCATE(base_intf_map(no_base_dofs,2))
      ALLOCATE(rearrange_index(no_base_dofs))      

c     Intiailize
c     ==========
      base_intf_map = 0
      rearrange_index = 0
      
c     Base Housing interface dofs
c     ===========================
      DO i = 1 , no_base_intf_dofs
        rearrange_index(base_intf_dofs(i,1)) = 1
      ENDDO
      
c     Rearranging DOFs
c     ================
      DO i = 1 , no_base_intf_dofs
        base_intf_map(i,1) = base_intf_dofs(i,1)
      ENDDO
      
      i_cnt = no_base_intf_dofs
      
      DO i = 1 , no_base_dofs
        IF (rearrange_index(i) .EQ. 0) THEN
            i_cnt = i_cnt + 1
            base_intf_map(i_cnt,1) = i
        ENDIF
      ENDDO
      
c     Inverse map
c     ===========
      DO i = 1 , no_base_dofs
        base_intf_map(base_intf_map(i,1),2) = i
      ENDDO

c     Deallocate temprory matrix
c     ==========================
      DEALLOCATE(rearrange_index)

c     Rearrange base plate matrices
c     =============================
      ALLOCATE(M_row(no_base_dofs,no_base_dofs))
      ALLOCATE(G_row(no_base_dofs,no_base_dofs))
      ALLOCATE(K_row(no_base_dofs,no_base_dofs))
      
      M_row = 0.0d0
      G_row = 0.0d0
      K_row = 0.0d0
      
      ALLOCATE(M_col(no_base_dofs,no_base_dofs))
      ALLOCATE(G_col(no_base_dofs,no_base_dofs))
      ALLOCATE(K_col(no_base_dofs,no_base_dofs))
      
      M_col = 0.0d0
      G_col = 0.0d0
      K_col = 0.0d0
      
c     Rearrange rows
c     ==============
      DO i = 1, no_base_dofs
        M_row(i,1:no_base_dofs) = 
     &      base_M(base_intf_map(i,1),1:no_base_dofs)
        G_row(i,1:no_base_dofs) = 
     &      base_G(base_intf_map(i,1),1:no_base_dofs)
        K_row(i,1:no_base_dofs) = 
     &      base_K(base_intf_map(i,1),1:no_base_dofs)
      ENDDO      

c     Rearrange columns
c     =================
      DO i = 1, no_base_dofs
        M_col(1:no_base_dofs,i) =
     &    M_row(1:no_base_dofs,base_intf_map(i,1))
        G_col(1:no_base_dofs,i) =
     &    G_row(1:no_base_dofs,base_intf_map(i,1))
        K_col(1:no_base_dofs,i) =
     &    K_row(1:no_base_dofs,base_intf_map(i,1))
      ENDDO
      
c     Rearranged housing matrices
c     ===========================  
      base_M = M_col
      base_G = G_col
      base_K = K_col
      
c     Deallocate temporary matrices
c     =============================
      DEALLOCATE(M_row)
      DEALLOCATE(G_row)
      DEALLOCATE(K_row)
     
      DEALLOCATE(M_col)
      DEALLOCATE(G_col)
      DEALLOCATE(K_col)

c     Base plate inertia load map
c     ===========================
      ALLOCATE(base_inrt_map(no_base_dofs,3))
      base_inrt_map = 0.0d0

      DO i = 1 , no_base_dofs
        dof_no = INT(base_cords(i,4))
        base_dof =  base_intf_map(i,2)   
        base_inrt_map(base_dof,dof_no) = 1.0d0
      ENDDO

c     Base plate DOFs accelaration vector
c     ===================================
      ALLOCATE(base_inrt_acc(no_base_dofs,1))
      base_inrt_acc = 0.0d0
      base_inrt_acc = MATMUL(base_inrt_map,shk_dir)
      DEALLOCATE(base_inrt_map)
      
      RETURN
      
999	WRITE(*,*) 'Trouble in opening base files'
	STOP
	           
      END

c======================================================================
      SUBROUTINE couple_hsng_base
c======================================================================

c     Shared Data
c     ===========
      USE hsng_data
      USE base_data
      
      IMPLICIT NONE
      
c     Local Variables
c     ===============
      INTEGER i,j,i_row,j_col

c     Temporary Matrices/Vectors
c     ==========================
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: asmb_inrt_acc
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: asmb_inrt_vec

c     DOF Data
c     ========
      no_asmb_dofs = no_hsng_dofs + no_base_dofs - no_base_intf_dofs

c     Intialize assembly matrics
c     ==========================
      ALLOCATE(asmb_M(no_asmb_dofs,no_asmb_dofs))
      ALLOCATE(asmb_G(no_asmb_dofs,no_asmb_dofs))
      ALLOCATE(asmb_K(no_asmb_dofs,no_asmb_dofs))
      
      asmb_M = 0.0d0
      asmb_G = 0.0d0
      asmb_K = 0.0d0
      
c     Housing Matrices
c     ================
      DO i = 1 , no_hsng_dofs
        DO j = 1 , no_hsng_dofs
            asmb_M(i,j) = hsng_M(i,j)
            asmb_G(i,j) = hsng_G(i,j)
            asmb_K(i,j) = hsng_K(i,j)
        ENDDO
      ENDDO

c     Base Plate Matrices
c     ===================
      DO i = 1 , no_base_dofs
        i_row = no_hsng_dofs - no_base_intf_dofs + i
        DO j = 1 , no_base_dofs
            j_col = no_hsng_dofs - no_base_intf_dofs + j
            asmb_M(i_row,j_col) = asmb_M(i_row,j_col) + base_M(i,j)
            asmb_G(i_row,j_col) = asmb_G(i_row,j_col) + base_G(i,j)
            asmb_K(i_row,j_col) = asmb_K(i_row,j_col) + base_K(i,j)
        ENDDO
      ENDDO
      
c     Deallocate system matrices
c     ==========================
      DEALLOCATE(hsng_M)
      DEALLOCATE(hsng_G)
      DEALLOCATE(hsng_K)
      
      DEALLOCATE(base_M)
      DEALLOCATE(base_G)
      DEALLOCATE(base_K)
      
c     Base plate housing assembly DOFs accelaration
c     =============================================
      ALLOCATE(asmb_inrt_acc(no_asmb_dofs,1))
      asmb_inrt_acc = 0.0d0
      
      DO i = 1 , no_hsng_dofs
        asmb_inrt_acc(i,1) = hsng_inrt_acc(i,1)
      ENDDO
      
      DO i = no_base_intf_dofs + 1, no_base_dofs
        i_row = no_hsng_dofs - no_base_intf_dofs + i
        asmb_inrt_acc(i_row,1) = base_inrt_acc(i,1)
      ENDDO
      
c     Deallocate inertia matrices
c     ===========================
      DEALLOCATE(hsng_inrt_acc)
      DEALLOCATE(base_inrt_acc)
      
      ALLOCATE(asmb_inrt_vec(no_asmb_dofs,1))
      asmb_inrt_vec = 0.0d0
      
      asmb_inrt_vec = MATMUL(asmb_M,asmb_inrt_acc)
      DEALLOCATE(asmb_inrt_acc)

      ALLOCATE(asmb_inrt_load(no_asmb_dofs))
      asmb_inrt_load = 0.0d0
      DO i = 1, no_asmb_dofs
        asmb_inrt_load(i) = asmb_inrt_vec(i,1)
      ENDDO
      DEALLOCATE(asmb_inrt_vec)

      RETURN
      END

c======================================================================
      SUBROUTINE couple_rotr_base
c======================================================================

c     Shared Data
c     ===========
      USE motr_data
      USE rotr_data
      USE hsng_data
      USE base_data
      USE supp_data
      
      IMPLICIT NONE
      
c     Local Variables
c     ===============
      INTEGER i,j,row_no,col_no
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: temp_mat

c     DOF data
c     ========     
      no_supp_dofs = no_rotr_dofs + no_asmb_dofs
      no_supp_vrbl = no_supp_dofs + no_disk_bcs
      no_fdb_dofs_motr = no_fdb_dofs_hub+no_fdb_dofs_hsng

c     Initialize support system matrices
c     ===================================
      ALLOCATE(supp_M(no_supp_dofs,no_supp_dofs))
      ALLOCATE(supp_G(no_supp_dofs,no_supp_dofs))
      ALLOCATE(supp_K(no_supp_dofs,no_supp_dofs))
      supp_M = 0.0d0
      supp_G = 0.0d0
      supp_K = 0.0d0
      
c     Coupling
c     ========

c     Add rotor matrices
      supp_M(1:no_rotr_dofs,1:no_rotr_dofs) = 
     &  rotr_M(1:no_rotr_dofs,1:no_rotr_dofs)
      supp_G(1:no_rotr_dofs,1:no_rotr_dofs) = 
     &  rotr_G(1:no_rotr_dofs,1:no_rotr_dofs)
      supp_K(1:no_rotr_dofs,1:no_rotr_dofs) = 
     &  rotr_K(1:no_rotr_dofs,1:no_rotr_dofs)
     
c     Add housing - base plate assembly matrices
      supp_M(no_rotr_dofs+1:no_rotr_dofs+no_asmb_dofs,
     &       no_rotr_dofs+1:no_rotr_dofs+no_asmb_dofs) = 
     &  asmb_M(1:no_asmb_dofs,1:no_asmb_dofs)
      supp_G(no_rotr_dofs+1:no_rotr_dofs+no_asmb_dofs,
     &       no_rotr_dofs+1:no_rotr_dofs+no_asmb_dofs) = 
     &  asmb_G(1:no_asmb_dofs,1:no_asmb_dofs)
      supp_K(no_rotr_dofs+1:no_rotr_dofs+no_asmb_dofs,
     &       no_rotr_dofs+1:no_rotr_dofs+no_asmb_dofs) = 
     &  asmb_K(1:no_asmb_dofs,1:no_asmb_dofs)
     
c     Deallocate system matrices
c     ==========================
      DEALLOCATE(rotr_M)
      DEALLOCATE(rotr_G)
      DEALLOCATE(rotr_K)
      
      DEALLOCATE(asmb_M)
      DEALLOCATE(asmb_G)
      DEALLOCATE(asmb_K)
      
c     Sparse Map for stationary hosuing - base plate assembly
c     =======================================================
      ALLOCATE(supp_M_map(no_supp_dofs,no_supp_dofs+1))
      CALL sparse_mat_map(supp_M,supp_M_map,
     &	no_supp_dofs,no_supp_dofs)

      ALLOCATE(temp_mat(no_supp_dofs,no_supp_dofs))
      
c     Damping matrix get affected by FDB coupling      
      temp_mat = supp_G
      DO i = 1 , no_fdb_dofs_motr
        row_no = no_rotr_dofs - no_fdb_dofs_hub + i
        DO j = 1 , no_fdb_dofs_motr
            col_no = no_rotr_dofs - no_fdb_dofs_hub + j
            temp_mat(row_no,col_no) = 1.0d0
        ENDDO
      ENDDO
      ALLOCATE(supp_G_map(no_supp_dofs,no_supp_dofs+1))
	CALL sparse_mat_map(temp_mat,supp_G_map,
     &	no_supp_dofs,no_supp_dofs)
     
c     Stiffness matrix get affected by FDB coupling      
      temp_mat = supp_K
      DO i = 1 , no_fdb_dofs_motr
        row_no = no_rotr_dofs - no_fdb_dofs_hub + i
        DO j = 1 , no_fdb_dofs_motr
            col_no = no_rotr_dofs - no_fdb_dofs_hub + j
            temp_mat(row_no,col_no) = 1.0d0
        ENDDO
      ENDDO
      ALLOCATE(supp_K_map(no_supp_dofs,no_supp_dofs+1))
	CALL sparse_mat_map(temp_mat,supp_K_map,
     &	no_supp_dofs,no_supp_dofs)
     
      DEALLOCATE(temp_mat)
     
c     Boundary Condition : Clamped disk center
c     ========================================     

c     Initialize
      ALLOCATE(C_supp(no_disk_bcs,no_supp_dofs))
      C_supp = 0.0d0
     
c     C_supp is same as C_rotr
      C_supp(1:no_disk_bcs,1:no_rotr_dofs) = 
     &  C_rotr(1:no_disk_bcs,1:no_rotr_dofs)
         
c     Deallocate
      DEALLOCATE(C_rotr)
      
c     Sparse map     
	ALLOCATE(C_supp_map(no_disk_bcs,no_supp_dofs+1))
	CALL sparse_mat_map(C_supp,C_supp_map,
     &	no_disk_bcs,no_supp_dofs)
      
c     Disk support system inertia load vector
c     =======================================
      ALLOCATE(supp_inrt_load(no_supp_dofs))
      supp_inrt_load = 0.0d0

c     -------------------------------------------------------------      
c     Note : Inertia matrix are uncoupled since FDB has no inertia
c     -------------------------------------------------------------

c     Add disk part of rotor inertia load vector
c     ==========================================
      supp_inrt_load(1:m_rotr) = rotr_inrt_load(1:m_rotr)

c     Add hosuing - base plate assembly inertia load vector
c     =====================================================
      supp_inrt_load(no_rotr_dofs+1:no_rotr_dofs+no_asmb_dofs) = 
     &  asmb_inrt_load(1:no_asmb_dofs)

c     Deallocate
c     ==========
      DEALLOCATE(rotr_inrt_load)
      DEALLOCATE(asmb_inrt_load)
      
      RETURN
      END
      
c======================================================================
      SUBROUTINE newmark_support
c======================================================================

c     Shared data
c     ===========
      USE sldr_dynm
      USE disk_data
      USE disk_aray
      USE motr_data
      USE rotr_data
      USE supp_data
      USE supp_soln

      IMPLICIT NONE

c     Local Variables
c     ===============
      INTEGER i,j,m,n,p 

c     Block Newmark Beta Matrices
c     ===========================
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_mat
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: C_supp_T
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_supp_T
	DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: temp_mat

      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_supp
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_supp_11
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_supp_12  
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_supp_13
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_supp_21
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_supp_23
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_supp_31
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_supp_32
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: A_supp_33

      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_supp_11
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_supp_21
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: L_supp_31
      
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_supp_11
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_supp_12
      DOUBLE PRECISION, DIMENSION (:,:), ALLOCATABLE :: U_supp_13
      
c     Newmark Beta Formulation
c	========================
	ALLOCATE(A_mat(no_supp_dofs,no_supp_dofs))
      A_mat = 0.0d0   ! Intialize
      A_mat = supp_M + (disk_nb_1*dt*supp_G) +
     &	(0.5d0*disk_nb_2*(dt**2)*supp_K)

c     C_supp_T for constrained formulation
c     =====================================
      ALLOCATE(C_supp_T(no_supp_dofs,no_disk_bcs))
	C_supp_T = 0.0d0
	C_supp_T(1:no_supp_dofs,1:no_disk_bcs) = 
     &	TRANSPOSE(C_supp(1:no_disk_bcs,1:no_supp_dofs))

c	Constraint BC at ID (Lagrange Multiplier Method)
c	================================================
      ALLOCATE(A_supp(no_supp_vrbl,no_supp_vrbl))
	A_supp = 0.0d0

      A_supp(1:no_supp_dofs,1:no_supp_dofs) = 
     &	A_mat(1:no_supp_dofs,1:no_supp_dofs)
     
      A_supp(no_supp_dofs+1:no_supp_vrbl,1:no_supp_dofs) = 
     &	C_supp(1:no_disk_bcs,1:no_supp_dofs)

      A_supp(1:no_supp_dofs,no_supp_dofs+1:no_supp_vrbl) = 
     &	C_supp_T(1:no_supp_dofs,1:no_disk_bcs)

      DEALLOCATE(A_mat)
      DEALLOCATE(C_supp_T)

c     Block factorization : Time Independent
c     ======================================      
      m = no_rotr_dofs - no_fdb_dofs_hub
      n = no_fdb_dofs_motr
      p = no_supp_vrbl - m - n
      
c     A11 = L11*U11
      ALLOCATE(A_supp_11(m,m))
      ALLOCATE(L_supp_11(m,m))
      ALLOCATE(U_supp_11(m,m))
      
      A_supp_11 = 0.0d0
      L_supp_11 = 0.0d0
      U_supp_11 = 0.0d0
      
      A_supp_11(1:m,1:m) = A_supp(1:m,1:m)
      CALL LU_fact(A_supp_11,L_supp_11,U_supp_11,m)
      DEALLOCATE(A_supp_11)

c     A12 = L11*U12
      ALLOCATE(A_supp_12(m,n))
      ALLOCATE(U_supp_12(m,n))
      
      A_supp_12 = 0.0d0
      U_supp_12 = 0.0d0
      
      A_supp_12(1:m,1:n) = A_supp(1:m,m+1:m+n)
      CALL L_solve_mult(L_supp_11,A_supp_12,U_supp_12,n,m)
      DEALLOCATE(A_supp_12)

c     A13 = L11*U13
      ALLOCATE(A_supp_13(m,p))
      ALLOCATE(U_supp_13(m,p))
      
      A_supp_13 = 0.0d0
      U_supp_13 = 0.0d0
      
      A_supp_13(1:m,1:p) = A_supp(1:m,m+n+1:m+n+p)
      CALL L_solve_mult(L_supp_11,A_supp_13,U_supp_13,p,m)
      DEALLOCATE(A_supp_13)
      
c     A21 = L21*U11 => A21_T = U11_T*L_21*T
      ALLOCATE(A_supp_21(n,m))
      ALLOCATE(L_supp_21(n,m))
      
      A_supp_21 = 0.0d0
      L_supp_21 = 0.0d0
      
      A_supp_21(1:n,1:m) = A_supp(m+1:m+n,1:m)
      CALL U_T_solve_mult(U_supp_11,A_supp_21,L_supp_21,n,m)
      DEALLOCATE(A_supp_21)
      
c     A31 = L31*U11 => A31_T = U11_T*L_31*T
      ALLOCATE(A_supp_31(p,m))
      ALLOCATE(L_supp_31(p,m))
      
      A_supp_31 = 0.0d0
      L_supp_31 = 0.0d0
      
      A_supp_31(1:p,1:m) = A_supp(m+n+1:m+n+p,1:m)
      CALL U_T_solve_mult(U_supp_11,A_supp_31,L_supp_31,p,m)
      DEALLOCATE(A_supp_31)
      
c     A22 = L21*U12 + L22U22 : L21*U12 = A_supp_22_m      
      ALLOCATE(A_supp_22_m(n,n))
      A_supp_22_m = MATMUL(L_supp_21,U_supp_12)
      
c     A23 = L21*U13 + L22*U23 => A23 - L21*U13 = A_supp_23_m
      ALLOCATE(A_supp_23(n,p))
      ALLOCATE(A_supp_23_m(n,p))
      
      A_supp_23   = 0.0d0
      A_supp_23_m = 0.0d0
      
      A_supp_23(1:n,1:p) = A_supp(m+1:m+n,m+n+1:m+n+p)
      A_supp_23_m = A_supp_23 - MATMUL(L_supp_21,U_supp_13)
      DEALLOCATE(A_supp_23)
      
c     A32 = L31*U12 + L32*U22 => A32 - L31*U12 = A_supp_32_m
      ALLOCATE(A_supp_32(p,n))
      ALLOCATE(A_supp_32_m(p,n))      
      
      A_supp_32   = 0.0d0
      A_supp_32_m = 0.0d0
      
      A_supp_32(1:p,1:n) = A_supp(m+n+1:m+n+p,m+1:m+n)
      A_supp_32_m = A_supp_32 - MATMUL(L_supp_31,U_supp_12)
      DEALLOCATE(A_supp_32)

c     A33 = L31*U13 + L32*U23 + L33*U33 => A33 - L31*U13 = A_supp_33_m
      ALLOCATE(A_supp_33(p,p))
      ALLOCATE(A_supp_33_m(p,p))
      
      A_supp_33   = 0.0d0
      A_supp_33_m = 0.0d0
      
      A_supp_33(1:p,1:p) = A_supp(m+n+1:m+n+p,m+n+1:m+n+p)
      A_supp_33_m = A_supp_33 - MATMUL(L_supp_31,U_supp_13)
      DEALLOCATE(A_supp_33)
      
c     Factorization done      
      DEALLOCATE(A_supp)
      
c     Assemble L-U Blocks
c     ===================
      ALLOCATE(L_supp(no_supp_vrbl,no_supp_vrbl))
      ALLOCATE(U_supp(no_supp_vrbl,no_supp_vrbl))
      
      L_supp = 0.0d0
      U_supp = 0.0d0
      
      L_supp(1:m,1:m) = L_supp_11(1:m,1:m)
      L_supp(m+1:m+n,1:m) = L_supp_21(1:n,1:m)
      L_supp(m+n+1:m+n+p,1:m) = L_supp_31(1:p,1:m)
      
      U_supp(1:m,1:m) = U_supp_11(1:m,1:m)
      U_supp(1:m,m+1:m+n) = U_supp_12(1:m,1:n)
      U_supp(1:m,m+n+1:m+n+p) = U_supp_13(1:m,1:p)

c     Sparse Map
c     ==========
      ALLOCATE(temp_mat(no_supp_vrbl,no_supp_vrbl))
      
c     L map      
      temp_mat = 0.0d0
      temp_mat = L_supp
      
c     L22
      DO i = m+1 , m+n
        DO j = m+1 , i
            temp_mat(i,j) = 1.0d0
        ENDDO
      ENDDO
      
c     L32
      DO i = m+n+1 , m+n+p
        DO j = m+1 , m+n
            temp_mat(i,j) = 1.0d0
        ENDDO
      ENDDO
      
c     L33 
      DO i = m+n+1 , m+n+p
        DO j = m+n+1 , i
            temp_mat(i,j) = 1.0d0
        ENDDO
      ENDDO

      ALLOCATE(L_supp_map(no_supp_vrbl,no_supp_vrbl+1))
      CALL sparse_mat_map(temp_mat,L_supp_map,
     &	no_supp_vrbl,no_supp_vrbl)
     
c     U map
      temp_mat = 0.0d0
      temp_mat = U_supp
      
c     U22
      DO i = m+1 , m+n
        DO j = i , m+n
            temp_mat(i,j) = 1.0d0
        ENDDO
      ENDDO
      
c     U23
      DO i = m+1 , m+n
        DO j = m+n+1 , m+n+p
            temp_mat(i,j) = 1.0d0
        ENDDO
      ENDDO
      
c     U33
      DO i = m+n+1 , m+n+p
        DO j = i , m+n+p
            temp_mat(i,j) = 1.0d0
        ENDDO
      ENDDO
      
      ALLOCATE(U_supp_map(no_supp_vrbl,no_supp_vrbl+1))
      CALL sparse_mat_map(temp_mat,U_supp_map,
     &	no_supp_vrbl,no_supp_vrbl)
     
c     U_T map     
      ALLOCATE(U_supp_T(no_supp_vrbl,no_supp_vrbl))
      U_supp_T = TRANSPOSE(temp_mat)
      
      ALLOCATE(U_supp_T_map(no_supp_vrbl,no_supp_vrbl+1))
      CALL sparse_mat_map(U_supp_T,U_supp_T_map,
     &	no_supp_vrbl,no_supp_vrbl)
      DEALLOCATE(U_supp_T)

      DEALLOCATE(temp_mat)
      
c     FDB Interface Block Matrices
c     ============================
      ALLOCATE(intf_M(n,n))
      ALLOCATE(intf_G(n,n))
      ALLOCATE(intf_K(n,n))
      
      intf_M = 0.0d0
      intf_G = 0.0d0
      intf_K = 0.0d0
      
      intf_M(1:n,1:n) = supp_M(m+1:m+n,m+1:m+n)
      intf_G(1:n,1:n) = supp_G(m+1:m+n,m+1:m+n)
      intf_K(1:n,1:n) = supp_K(m+1:m+n,m+1:m+n)   
      
c     Initialize matrices used in motor solver
c     ========================================
      
c     Initialize circumferential interpolation
      ALLOCATE(crcm_intp_mat(no_crcm_node_hub,4))
      crcm_intp_mat = 0.0d0
      
c     Hub current angle list
      ALLOCATE(hub_cur_ang(no_crcm_node_hub))
      hub_cur_ang = 0.0d0
      
c     Temporary Newmark Beta Block Matrices
      ALLOCATE(temp_G_s(n,n))
      ALLOCATE(temp_K_s(n,n))
      temp_G_s = 0.0d0
      temp_K_s = 0.0d0
      
c     Newmark Beta Block Matrices
      ALLOCATE(A_supp_22(n,n))
      ALLOCATE(L_supp_22(n,n))
      ALLOCATE(U_supp_22(n,n))
            
      ALLOCATE(U_supp_23(n,p))
      ALLOCATE(L_supp_32(p,n))
      
      ALLOCATE(L_supp_33(p,p))
      ALLOCATE(U_supp_33(p,p))
      
      A_supp_22 = 0.0d0
      L_supp_22 = 0.0d0
      U_supp_22 = 0.0d0
      
      L_supp_32 = 0.0d0
      U_supp_23 = 0.0d0
      
      L_supp_33 = 0.0d0
      U_supp_33 = 0.0d0
      
      ALLOCATE(blck_supp_1(n,n))
      ALLOCATE(blck_supp_2(p,p))
      blck_supp_1 = 0.0d0
      blck_supp_2 = 0.0d0
      
c     Block Matrices size      
      m_supp = m
      n_supp = n
      p_supp = p
      
c	Initialize Dynamic Variables
      ALLOCATE(w_disk_cur(no_disk_dofs))
	ALLOCATE(w_supp_cur(no_supp_dofs))
	ALLOCATE(v_supp_cur(no_supp_dofs))
	ALLOCATE(a_supp_cur(no_supp_dofs))
	ALLOCATE(w_supp_prv(no_supp_dofs))
	ALLOCATE(v_supp_prv(no_supp_dofs))
	ALLOCATE(a_supp_prv(no_supp_dofs))
	ALLOCATE(f_supp_cnst(no_disk_bcs))
	
c     Disk support system solution variables
      ALLOCATE(w_supp_prd(no_supp_dofs))
	ALLOCATE(v_supp_prd(no_supp_dofs))
	ALLOCATE(f_supp_val(no_supp_dofs))

	ALLOCATE(supp_vec_1(no_supp_dofs))
	ALLOCATE(supp_vec_2(no_supp_dofs))
	ALLOCATE(supp_vec_3(no_disk_bcs))
	
      ALLOCATE(x_supp(no_supp_vrbl))
	ALLOCATE(y_supp(no_supp_vrbl))
	ALLOCATE(z_supp(no_supp_vrbl))
	
      ALLOCATE(ramp_coef_supp(no_disk_dofs))
	ALLOCATE(C_supp_r(no_supp_vrbl))
	ALLOCATE(L_21_supp(no_supp_vrbl))
	ALLOCATE(U_12_supp(no_supp_vrbl))
	
	C_supp_r = 0.0d0  ! Initialize
      
      RETURN
      END
      
c======================================================================
      SUBROUTINE solve_support
c======================================================================

c     Shared data
c     ===========
      USE sldr_dynm
      USE shck_data
      USE disk_data
      USE disk_aray
      USE motr_data
      USE rotr_data
      USE hsng_data
      USE supp_data    
      USE supp_soln
      
      IMPLICIT NONE
      
c     Parameters
c     ==========
      DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979      
      
c     Local variables
c     ===============
      INTEGER i,j,i_brng,i_node,brng_row,pvt
      INTEGER FDB_node_map(9),hsng_node_no
      DOUBLE PRECISION ang_pos,hub_ang
      DOUBLE PRECISION K_local(3,3),C_local(3,3)
      DOUBLE PRECISION K_nodal(9,9),C_nodal(9,9)
      
      INTEGER rad_gridpt(2),ang_gridpt(2)
      DOUBLE PRECISION time_fac
      DOUBLE PRECISION zeta,eta,zeta_temp
      DOUBLE PRECISION w_cnst,w_ramp,U_22,lmbd,nu_val
	DOUBLE PRECISION ang_cur,ang_vel,eps_ang
	DOUBLE PRECISION xtemp,ytemp,rad_temp,ang_temp
	
c	Angular Speed
	ang_vel = disk_RPM*(2.0d0*pi/60.0d0)
	
c     Current angular postion in HDD frame
      ang_pos = ang_vel*t_val
      
c     Circumferential interpolation matrix
      CALL crcm_fdb_intp(ang_pos)
      
c     Hub current angle positions
      DO i = 1 , no_crcm_node_hub
        hub_cur_ang(i) = ang_pos + crcm_ang_hub(i)
      ENDDO
      
c     Initialize FDB interface dynamic coefficient matrices
      temp_G_s = intf_G
      temp_K_s = intf_K
      
c     Add FDB Stiffness and damping
      DO i_brng = 1 , 3 ! Bearing type : 3
      
        brng_row  = 3*(i_brng-1)
        
c       Reset    
        C_local = 0.0d0    
        K_local = 0.0d0
        
c       Local Dynanic Coefficients
        DO i = 1 , 3
            DO j = 1 , 3
                C_local(i,j) = C_fdb(brng_row+i,j)/
     &                         DBLE(no_crcm_node_hub)
                K_local(i,j) = K_fdb(brng_row+i,j)/
     &                         DBLE(no_crcm_node_hub)
            ENDDO
        ENDDO
        
c       Add nodal coefficents
        DO i_node = 1 , no_crcm_node_hub
        
c           Nodal Dynamic Coefficients
            CALL nodal_FDB_dyn_coef(C_local,K_local,C_nodal,K_nodal,
     &           ang_pos,crcm_intp_mat(i_node,1),
     &           crcm_intp_mat(i_node,2))

c           Find location of inteface FDB nodes on hub
            pvt = 3*no_crcm_node_hub*(i_brng-1) + 3*(i_node-1)
            DO j = 1 , 3
                FDB_node_map(j) = pvt + j
            ENDDO

c           Find location of interface FDB nodes on housing
            DO i = 1 , 2
                hsng_node_no = INT(crcm_intp_mat(i_node,i+2))
                pvt = (3*3*no_crcm_node_hub) + (3*no_crcm_node_hsng*
     &                (i_brng-1)) + 3*(hsng_node_no-1)
                DO j = 1 , 3
                    FDB_node_map(3*(i-1)+j+3) = pvt + j
                ENDDO
            ENDDO
            
            DO i = 1 , 9
                DO j = 1 , 9
                
                    temp_G_s(FDB_node_map(i),FDB_node_map(j)) = 
     &                  temp_G_s(FDB_node_map(i),FDB_node_map(j)) +
     &                  C_nodal(i,j)
     
                    temp_K_s(FDB_node_map(i),FDB_node_map(j)) = 
     &                  temp_K_s(FDB_node_map(i),FDB_node_map(j)) +
     &                  K_nodal(i,j)
     
                ENDDO
            ENDDO            
        
        ENDDO
        
      ENDDO  

c     Block Factorization
c     ===================

c     Newmark Beta Block Matrix : A22
      A_supp_22 = intf_M + (disk_nb_1*dt*temp_G_s) +
     &	(0.5d0*disk_nb_2*(dt**2)*temp_K_s)
     
c     A22 - L21*U12 = A_supp_22_m = L22*U22
      blck_supp_1 = A_supp_22 - A_supp_22_m
      CALL LU_fact(blck_supp_1,L_supp_22,U_supp_22,n_supp)

c     A23 - L21*U13 = A_supp_23_m = L22*U23
      CALL L_solve_mult(L_supp_22,A_supp_23_m,U_supp_23,
     &  p_supp,n_supp)
     
c     A32 - L31*U12 = A_supp_32_m = L32*U22 => U22'*L32' = A_supp_32_m'
      CALL U_T_solve_mult(U_supp_22,A_supp_32_m,L_supp_32,
     &  p_supp,n_supp)
     
c     A33 - L31*U13 = A_supp_33_m = L32*U23 + L33*U33
      blck_supp_2 = A_supp_33_m - MATMUL(L_supp_32,U_supp_23)
      CALL LU_fact(blck_supp_2,L_supp_33,U_supp_33,p_supp)
      
c     Assemble L-U Matrices
c     =====================

      L_supp(m_supp+1:m_supp+n_supp,
     &  m_supp+1:m_supp+n_supp) = 
     &  L_supp_22(1:n_supp,1:n_supp)                    ! L22
     
      L_supp(m_supp+n_supp+1:m_supp+n_supp+p_supp,
     &  m_supp+1:m_supp+n_supp) = 
     &  L_supp_32(1:p_supp,1:n_supp)                    ! L32

      L_supp(m_supp+n_supp+1:m_supp+n_supp+p_supp,
     &  m_supp+n_supp+1:m_supp+n_supp+p_supp) = 
     &  L_supp_33(1:p_supp,1:p_supp)                    ! L33     
     
      U_supp(m_supp+1:m_supp+n_supp,
     &  m_supp+1:m_supp+n_supp) = 
     &  U_supp_22(1:n_supp,1:n_supp)                    ! U22
     
      U_supp(m_supp+1:m_supp+n_supp,
     &  m_supp+n_supp+1:m_supp+n_supp+p_supp) = 
     &  U_supp_23(1:n_supp,1:p_supp)                    ! U23
     
      U_supp(m_supp+n_supp+1:m_supp+n_supp+p_supp,
     &  m_supp+n_supp+1:m_supp+n_supp+p_supp) = 
     &  U_supp_33(1:p_supp,1:p_supp)                    ! U33
     
c     Update Motor damping and stiffness matrices
c     ===========================================
      supp_G(m_supp+1:m_supp+n_supp,
     &  m_supp+1:m_supp+n_supp) = 
     &  temp_G_s(1:n_supp,1:n_supp)

      supp_K(m_supp+1:m_supp+n_supp,
     &  m_supp+1:m_supp+n_supp) = 
     &  temp_K_s(1:n_supp,1:n_supp)
     
c	Predicted displacement
c	======================
	w_supp_prd = w_supp_prv + (dt*v_supp_prv) + 
     &	(0.5d0*(1.0d0-disk_nb_2)*(dt**2)*a_supp_prv)
     
c	Predicted velocity
c	==================
	v_supp_prd = v_supp_prv + ((1.0d0-disk_nb_1)*dt*a_supp_prv)

c	Inertia load vector
c	===================	
	
	CALL calc_hub_inrt(ang_vel,t_val)
	
c     Add hub part of rotor inertia vector
      DO i = 1 , no_hub_dofs
        supp_inrt_load(m_rotr+i) = hub_inrt_vec(i,1)
      ENDDO	
     
      f_supp_val = -acc_val*supp_inrt_load
      
c	Newmark Beta without Ramp Contact
c	================================= 
      CALL d_mult_sparse_mat_vec(supp_G,supp_G_map,v_supp_prd,
     &	supp_vec_1,no_supp_dofs,no_supp_dofs,no_supp_dofs)
	CALL d_mult_sparse_mat_vec(supp_K,supp_K_map,w_supp_prd,
     &	supp_vec_2,no_supp_dofs,no_supp_dofs,no_supp_dofs)
      CALL d_mult_sparse_mat_vec(C_supp,C_supp_map,w_supp_prd,
     &	supp_vec_3,no_disk_bcs,no_supp_dofs,no_supp_dofs)
      time_fac = -(2.0d0/((dt**2)*disk_nb_2))
      supp_vec_3 = time_fac*supp_vec_3

	y_supp(1:no_supp_dofs) = f_supp_val - supp_vec_1 - supp_vec_2
	y_supp(no_supp_dofs+1:no_supp_vrbl) = supp_vec_3
	
c     Solution of system of equations
c     ===============================
	CALL LU_solve_sparse(L_supp,U_supp,L_supp_map,U_supp_map,
     &	    y_supp(1:no_supp_vrbl),x_supp(1:no_supp_vrbl),
     &      no_supp_vrbl)

c	Compute current kinematic Value
c	===============================     
	a_supp_cur = x_supp(1:no_supp_dofs)
 	v_supp_cur = v_supp_prd + (disk_nb_1*dt*a_supp_cur)
	w_supp_cur = w_supp_prd + (0.5d0*disk_nb_2*(dt**2)*a_supp_cur)
	
c	Constraint Force
c	================
	f_supp_cnst(1:no_disk_bcs) = x_supp(no_supp_dofs+1:no_supp_vrbl)
	
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
	  
c		Ramp Location
c		=============
		ang_cur = ramp_ang  - (ang_vel*t_val)
		
c       Adjust for disk center movement
c       ===============================
        xtemp = ramp_rad*DCOS(ang_cur) - x0_disk_s
        ytemp = ramp_rad*DSIN(ang_cur) - y0_disk_s

        rad_temp = DSQRT((xtemp**2) + (ytemp**2))
        ang_temp = ATAN(ABS(ytemp/xtemp))
        
        IF ((xtemp .LT. 0) .AND. (ytemp .LT. 0)) THEN
            ang_temp = pi + ang_temp
        ELSEIF ((xtemp .GE. 0) .AND. (ytemp .LT. 0)) THEN
            ang_temp = 2*pi - ang_temp
        ELSEIF ((xtemp .LT. 0) .AND. (ytemp .GE. 0)) THEN
            ang_temp =   pi - ang_temp
        ENDIF
        
c		Position of ramp lowest point in disk (rotating) frame
c		======================================================
		CALL rad_loc(rad_temp,rad_grid,rad_gridpt,zeta)
		CALL ang_loc(ang_temp,ang_grid,ang_gridpt,eta)
		zeta_temp = zeta
		
c		DOF Coefficient for w_ramp displacement
c		=======================================
  		CALL ramp_cnstrnt(rad_gridpt,ang_gridpt,
     &		zeta_temp,eta,ramp_coef_supp)
     
c       Rearrange using mapping
c       =======================
        C_supp_r(1:no_disk_dofs) = 0.0d0
        DO i = 1, no_disk_dofs
            C_supp_r(disk_intf_map(i,2)) = ramp_coef_supp(i)
        ENDDO     
     
c		Disk Displacement under the ramp
c		================================
		w_ramp = DOT_PRODUCT(C_supp_r(1:no_disk_dofs),
     &		w_supp_cur(1:no_disk_dofs)) 
     
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
        
            w_ramp = DOT_PRODUCT(C_supp_r(1:no_disk_dofs),
     &		    w_supp_prd(1:no_disk_dofs))

			lmbd = (2.0d0/((dt**2)*disk_nb_2))*(w_cnst - w_ramp)

c           Block Factorization
c           ===================
            CALL L_solve_sparse(L_supp,L_supp_map,
     &			    C_supp_r,U_12_supp,no_supp_vrbl)
			CALL U_T_solve_sparse(U_supp,U_supp_T_map,
     &			    C_supp_r,L_21_supp,no_supp_vrbl)
			U_22 = DOT_PRODUCT(L_21_supp(1:no_supp_vrbl),
     &		        U_12_supp(1:no_supp_vrbl))
     
c           Solution of system of equations
c           ===============================

c           Forward : L solution
            CALL L_solve_sparse(L_supp,L_supp_map,y_supp,
     &           z_supp,no_supp_vrbl)
            nu_val = lmbd - DOT_PRODUCT(L_21_supp,z_supp)
            
c           Ramp Contact Force
            ramp_cntf = -nu_val/U_22
            
c           Backward : U solution
            DO i = 1 , no_supp_vrbl
                z_supp(i) = z_supp(i) - ramp_cntf*U_12_supp(i)
            ENDDO
            
            CALL U_solve_sparse(U_supp,U_supp_map,z_supp,
     &           x_supp,no_supp_vrbl)
     
c			Compute Current Kinematic Value
c			===============================
			a_supp_cur = x_supp(1:no_supp_dofs)
 			v_supp_cur = v_supp_prd + (disk_nb_1*dt*a_supp_cur)
			w_supp_cur = w_supp_prd + 
     &		          (0.5d0*disk_nb_2*(dt**2)*a_supp_cur)
     
c	      Disk ID Constraint Forces
c	      =========================
	      f_supp_cnst(1:no_disk_bcs) = 
     &	        x_supp(no_supp_dofs+1:no_supp_vrbl)       
        
 	  ENDIF
 	  
 	ENDIF

c     Compute disk displacement
c     =========================
      DO i = 1 , no_disk_dofs
        w_disk_cur(disk_intf_map(i,1)) = w_supp_cur(i)
      ENDDO
      
c     Disk ID displacement
c     ====================
      DO i = 1, no_disk_intf_dofs
        DO j = 1 , 3
            disk_ID_disp(i,j) = w_supp_cur(disk_ID_map(i,j+2))
        ENDDO
      ENDDO
      
c     FDB dispalcement
c     ================        
      IF (data_opt_supp.EQ.1) THEN
        w_FDB(1) = t_val

        pvt = no_rotr_dofs - no_fdb_dofs_hub
        
        DO i = 1 , no_fdb_dofs_hub 
            w_FDB(i+1) = w_supp_cur(pvt+i)
        ENDDO  

        DO i = 1 , no_fdb_dofs_hsng
            w_FDB(i+no_fdb_dofs_hub+1) = 
     &          w_supp_cur(no_rotr_dofs+i)
        ENDDO
  
      ENDIF                

	
      RETURN
      END
            
c======================================================================
      SUBROUTINE calc_hub_inrt(ang_vel,t_val)
c======================================================================

c     Shared Data
c     ===========
      USE shck_data
      USE motr_data
      USE rotr_data

      IMPLICIT NONE

c     Arguments
c     =========
      DOUBLE PRECISION ang_vel,t_val
      
c     Local Variables
c     ===============
      INTEGER i
      DOUBLE PRECISION ang_frame,csk,ssk,shk_dir_rot(3,1)
      
c     Inertia load direction in rotational frame
c     ==========================================
      ang_frame = ang_vel*t_val
      csk = DCOS(ang_frame)
      ssk = DSIN(ang_frame)
      
      shk_dir_rot(1,1) =  csk*shk_dir(1,1) + ssk*shk_dir(2,1)
      shk_dir_rot(2,1) = -ssk*shk_dir(1,1) + csk*shk_dir(2,1)
      shk_dir_rot(3,1) = shk_dir(3,1)
      
c     Hub DOFs accelaration vector
c     ============================
      hub_inrt_acc = MATMUL(hub_inrt_map,shk_dir_rot)
      
c     Hub inertia load
c     ================
      hub_inrt_vec = MATMUL(hub_M,hub_inrt_acc)



      RETURN    
      END
      
c======================================================================
      SUBROUTINE crcm_fdb_intp(ang_pos)
c======================================================================

c     Shared data
c     ===========
      USE motr_data
      USE motr_soln
      
      IMPLICIT NONE
      
c     Parameters
c     ==========
      DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979      
      
c     Argument
c     ========
      DOUBLE PRECISION ang_pos
      
c     Local variables
c     ===============
      INTEGER i,j,cnt
      LOGICAL flag
      DOUBLE PRECISION den,num_1,num_2
      DOUBLE PRECISION ang_cur,ang_1,ang_2
      
c     Reset matrix
      crcm_intp_mat = 0.0d0

c     Single row represent a hub node information      
      DO i = 1, no_crcm_node_hub
      
c       Current position      
        ang_cur = crcm_ang_hub(i) + ang_pos 
        
c       Input angle in 0 to 2*pi range
        IF (ang_cur .GT. 2*pi) THEN
            DO WHILE (ang_cur .GT. 2.0d0*pi)
                ang_cur = ang_cur - 2.0d0*pi
            ENDDO
        ELSE IF (ang_cur .LT. 0.0d0) THEN
            DO WHILE (ang_cur .LT. 0.0d0)
                ang_cur = ang_cur + 2.0d0*pi
            ENDDO
        ENDIF 
        
        flag = .TRUE.
        cnt  = 1
c       First no_crcm_node_hsng-1 sectors
        DO WHILE ((flag) .AND. (cnt .LE. no_crcm_node_hsng-1))
            ang_1 = crcm_ang_hsng(cnt)
            ang_2 = crcm_ang_hsng(cnt+1)
c           Coincide with the first node
            IF (DABS(ang_cur-ang_1) .LT. 1.0d-3) THEN
                crcm_intp_mat(i,1) = 1.0d0
                crcm_intp_mat(i,2) = 0.0d0
                crcm_intp_mat(i,3) = DBLE(cnt)
                crcm_intp_mat(i,4) = DBLE(cnt+1)
                flag = .FALSE.
c           Coincide with the second node
            ELSEIF (DABS(ang_cur-ang_2) .LT. 1.0d-3) THEN
                crcm_intp_mat(i,1) = 0.0d0
                crcm_intp_mat(i,2) = 1.0d0
                crcm_intp_mat(i,3) = DBLE(cnt)
                crcm_intp_mat(i,4) = DBLE(cnt+1)
                flag = .FALSE.
c           Lies between first and second node
            ELSEIF ((ang_1 .LT. ang_cur).AND.(ang_cur .LT. ang_2))THEN
                den = ang_2 - ang_1
                num_1 = ang_2 - ang_cur
                num_2 = ang_cur - ang_1
                crcm_intp_mat(i,1) = num_1/den
                crcm_intp_mat(i,2) = num_2/den
                crcm_intp_mat(i,3) = DBLE(cnt)
                crcm_intp_mat(i,4) = DBLE(cnt+1)
                
                flag = .FALSE.         
c           Not in the current interval
            ELSE
                cnt = cnt + 1
            ENDIF
        ENDDO
        
c       Last interval
        IF (flag) THEN
            ang_1 = crcm_ang_hsng(no_crcm_node_hsng)
            ang_2 = crcm_ang_hsng(1) + (2*pi)
c           Coincide with the first node
            IF (DABS(ang_cur-ang_1) .LT. 1.0d-3) THEN
                crcm_intp_mat(i,1) = 1.0d0
                crcm_intp_mat(i,2) = 0.0d0
                crcm_intp_mat(i,3) = DBLE(no_crcm_node_hsng)
                crcm_intp_mat(i,4) = DBLE(1)
c           Lies between last and first node
            ELSEIF ((ang_1 .LT. ang_cur).AND.(ang_cur .LT. ang_2))THEN
                den = ang_2 - ang_1
                num_1 = ang_2 - ang_cur
                num_2 = ang_cur - ang_1
                crcm_intp_mat(i,1) = num_1/den
                crcm_intp_mat(i,2) = num_2/den
                crcm_intp_mat(i,3) = DBLE(no_crcm_node_hsng)
                crcm_intp_mat(i,4) = DBLE(1)
c           Node cant be located
            ELSE
                WRITE(*,*) ' Node cant be located'
                RETURN
            ENDIF
        ENDIF        
      ENDDO
      

      RETURN    
      END 

c======================================================================
      SUBROUTINE nodal_FDB_dyn_coef(C,K,C_nodal,K_nodal,theta,c1,c2)
c====================================================================== 

      IMPLICIT NONE 
      
c     Argument
c     ========
      DOUBLE PRECISION c1,c2,theta
      DOUBLE PRECISION C(3,3),C_nodal(9,9)
      DOUBLE PRECISION K(3,3),K_nodal(9,9)
      
c     Local Variable
      DOUBLE PRECISION c_tht,s_tht
      
c     Components      
      c_tht = COS(theta)
      s_tht = SIN(theta)
 
c     Damping
c     =======

c     Row : 1      
      C_nodal(1,1) = C(1,1)
      C_nodal(1,2) = C(1,2)
      C_nodal(1,3) = 0
      
      C_nodal(1,4) = -c1*((C(1,1)*c_tht)-(C(1,2)*s_tht))
      C_nodal(1,5) = -c1*((C(1,2)*c_tht)+(C(1,1)*s_tht))
      C_nodal(1,6) = 0
      
      C_nodal(1,7) = -c2*((C(1,1)*c_tht)-(C(1,2)*s_tht))
      C_nodal(1,8) = -c2*((C(1,2)*c_tht)+(C(1,1)*s_tht))
      C_nodal(1,9) = 0

c     Row : 2      
      C_nodal(2,1) = C(2,1)
      C_nodal(2,2) = C(2,2)
      C_nodal(2,3) = 0
      
      C_nodal(2,4) = -c1*((C(2,1)*c_tht)-(C(2,2)*s_tht))
      C_nodal(2,5) = -c1*((C(2,2)*c_tht)+(C(2,1)*s_tht))
      C_nodal(2,6) = 0
      
      C_nodal(2,7) = -c2*((C(2,1)*c_tht)-(C(2,2)*s_tht))
      C_nodal(2,8) = -c2*((C(2,2)*c_tht)+(C(2,1)*s_tht))
      C_nodal(2,9) = 0 

c     Row : 3      
      C_nodal(3,1) = 0
      C_nodal(3,2) = 0
      C_nodal(3,3) = C(3,3)
      
      C_nodal(3,4) = 0
      C_nodal(3,5) = 0
      C_nodal(3,6) = -c1*C(3,3)
      
      C_nodal(3,7) = 0
      C_nodal(3,8) = 0
      C_nodal(3,9) = -c2*C(3,3)

c     Row : 4 
      C_nodal(4,1) = -c1*((C(1,1)*c_tht)-(C(2,1)*s_tht))
      C_nodal(4,2) = -c1*((C(1,2)*c_tht)-(C(2,2)*s_tht))
      C_nodal(4,3) = 0
      
      C_nodal(4,4) = (c1**2)*((C(1,1)*(c_tht**2))-((C(1,2)+C(2,1))*
     &                c_tht*s_tht)+(C(2,2)*(s_tht**2)))
      C_nodal(4,5) = (c1**2)*((C(1,2)*(c_tht**2))+((C(1,1)-C(2,2))*
     &                c_tht*s_tht)-(C(2,1)*(s_tht**2)))
      C_nodal(4,6) = 0
      
      C_nodal(4,7) = (c1*c2)*((C(1,1)*(c_tht**2))-((C(1,2)+C(2,1))*
     &                c_tht*s_tht)+(C(2,2)*(s_tht**2)))
      C_nodal(4,8) = (c1*c2)*((C(1,2)*(c_tht**2))+((C(1,1)-C(2,2))*
     &                c_tht*s_tht)-(C(2,1)*(s_tht**2)))
      C_nodal(4,9) = 0
      
c     Row : 5
      C_nodal(5,1) = -c1*((C(2,1)*c_tht)+(C(1,1)*s_tht))
      C_nodal(5,2) = -c1*((C(2,2)*c_tht)+(C(1,2)*s_tht))
      C_nodal(5,3) = 0
      
      C_nodal(5,4) = (c1**2)*((C(2,1)*(c_tht**2))+((C(1,1)-C(2,2))*
     &                c_tht*s_tht)-(C(1,2)*(s_tht**2)))
      C_nodal(5,5) = (c1**2)*((C(2,2)*(c_tht**2))+((C(1,2)+C(2,1))*
     &                c_tht*s_tht)+(C(1,1)*(s_tht**2)))
      C_nodal(5,6) = 0
      
      C_nodal(5,7) = (c1*c2)*((C(2,1)*(c_tht**2))+((C(1,1)-C(2,2))*
     &                c_tht*s_tht)-(C(1,2)*(s_tht**2)))
      C_nodal(5,8) = (c1*c2)*((C(2,2)*(c_tht**2))+((C(1,2)+C(2,1))*
     &                c_tht*s_tht)+(C(1,1)*(s_tht**2)))
      C_nodal(5,9) = 0

c     Row : 6      
      C_nodal(6,1) = 0
      C_nodal(6,2) = 0
      C_nodal(6,3) = -c1*C(3,3)
      
      C_nodal(6,4) = 0
      C_nodal(6,5) = 0
      C_nodal(6,6) = (c1**2)*C(3,3)
      
      C_nodal(6,7) = 0
      C_nodal(6,8) = 0
      C_nodal(6,9) = (c1*c2)*C(3,3)

c     Row : 7 
      C_nodal(7,1) = -c2*((C(1,1)*c_tht)-(C(2,1)*s_tht))
      C_nodal(7,2) = -c2*((C(1,2)*c_tht)-(C(2,2)*s_tht))
      C_nodal(7,3) = 0
      
      C_nodal(7,4) = (c1*c2)*((C(1,1)*(c_tht**2))-((C(1,2)+C(2,1))*
     &                c_tht*s_tht)+(C(2,2)*(s_tht**2)))
      C_nodal(7,5) = (c1*c2)*((C(1,2)*(c_tht**2))+((C(1,1)-C(2,2))*
     &                c_tht*s_tht)-(C(2,1)*(s_tht**2)))
      C_nodal(7,6) = 0
      
      C_nodal(7,7) = (c2**2)*((C(1,1)*(c_tht**2))-((C(1,2)+C(2,1))*
     &                c_tht*s_tht)+(C(2,2)*(s_tht**2)))
      C_nodal(7,8) = (c2**2)*((C(1,2)*(c_tht**2))+((C(1,1)-C(2,2))*
     &                c_tht*s_tht)-(C(2,1)*(s_tht**2)))
      C_nodal(7,9) = 0

c     Row : 8
      C_nodal(8,1) = -c2*((C(2,1)*c_tht)+(C(1,1)*s_tht))
      C_nodal(8,2) = -c2*((C(2,2)*c_tht)+(C(1,2)*s_tht))
      C_nodal(8,3) = 0
      
      C_nodal(8,4) = (c1*c2)*((C(2,1)*(c_tht**2))+((C(1,1)-C(2,2))*
     &                c_tht*s_tht)-(C(1,2)*(s_tht**2)))
      C_nodal(8,5) = (c1*c2)*((C(2,2)*(c_tht**2))+((C(1,2)+C(2,1))*
     &                c_tht*s_tht)+(C(1,1)*(s_tht**2)))
      C_nodal(8,6) = 0
      
      C_nodal(8,7) = (c2**2)*((C(2,1)*(c_tht**2))+((C(1,1)-C(2,2))*
     &                c_tht*s_tht)-(C(1,2)*(s_tht**2)))
      C_nodal(8,8) = (c2**2)*((C(2,2)*(c_tht**2))+((C(1,2)+C(2,1))*
     &                c_tht*s_tht)+(C(1,1)*(s_tht**2)))
      C_nodal(8,9) = 0

c     Row : 9      
      C_nodal(9,1) = 0
      C_nodal(9,2) = 0
      C_nodal(9,3) = -c2*C(3,3)
      
      C_nodal(9,4) = 0
      C_nodal(9,5) = 0
      C_nodal(9,6) = (c1*c2)*C(3,3)
      
      C_nodal(9,7) = 0
      C_nodal(9,8) = 0
      C_nodal(9,9) = (c2**2)*C(3,3)

c     Stiffness
c     =========
      
c     Row : 1      
      K_nodal(1,1) = K(1,1)
      K_nodal(1,2) = K(1,2)
      K_nodal(1,3) = 0
      
      K_nodal(1,4) = -c1*((K(1,1)*c_tht)-(K(1,2)*s_tht))
      K_nodal(1,5) = -c1*((K(1,2)*c_tht)+(K(1,1)*s_tht))
      K_nodal(1,6) = 0
      
      K_nodal(1,7) = -c2*((K(1,1)*c_tht)-(K(1,2)*s_tht))
      K_nodal(1,8) = -c2*((K(1,2)*c_tht)+(K(1,1)*s_tht))
      K_nodal(1,9) = 0     
      
c     Row : 2      
      K_nodal(2,1) = K(2,1)
      K_nodal(2,2) = K(2,2)
      K_nodal(2,3) = 0
      
      K_nodal(2,4) = -c1*((K(2,1)*c_tht)-(K(2,2)*s_tht))
      K_nodal(2,5) = -c1*((K(2,2)*c_tht)+(K(2,1)*s_tht))
      K_nodal(2,6) = 0
      
      K_nodal(2,7) = -c2*((K(2,1)*c_tht)-(K(2,2)*s_tht))
      K_nodal(2,8) = -c2*((K(2,2)*c_tht)+(K(2,1)*s_tht))
      K_nodal(2,9) = 0 
      
c     Row : 3      
      K_nodal(3,1) = 0
      K_nodal(3,2) = 0
      K_nodal(3,3) = K(3,3)
      
      K_nodal(3,4) = 0
      K_nodal(3,5) = 0
      K_nodal(3,6) = -c1*K(3,3)
      
      K_nodal(3,7) = 0
      K_nodal(3,8) = 0
      K_nodal(3,9) = -c2*K(3,3)

c     Row : 4 
      K_nodal(4,1) = -c1*((K(1,1)*c_tht)-(K(2,1)*s_tht))
      K_nodal(4,2) = -c1*((K(1,2)*c_tht)-(K(2,2)*s_tht))
      K_nodal(4,3) = 0
      
      K_nodal(4,4) = (c1**2)*((K(1,1)*(c_tht**2))-((K(1,2)+K(2,1))*
     &                c_tht*s_tht)+(K(2,2)*(s_tht**2)))
      K_nodal(4,5) = (c1**2)*((K(1,2)*(c_tht**2))+((K(1,1)-K(2,2))*
     &                c_tht*s_tht)-(K(2,1)*(s_tht**2)))
      K_nodal(4,6) = 0
      
      K_nodal(4,7) = (c1*c2)*((K(1,1)*(c_tht**2))-((K(1,2)+K(2,1))*
     &                c_tht*s_tht)+(K(2,2)*(s_tht**2)))
      K_nodal(4,8) = (c1*c2)*((K(1,2)*(c_tht**2))+((K(1,1)-K(2,2))*
     &                c_tht*s_tht)-(K(2,1)*(s_tht**2)))
      K_nodal(4,9) = 0
      
c     Row : 5
      K_nodal(5,1) = -c1*((K(2,1)*c_tht)+(K(1,1)*s_tht))
      K_nodal(5,2) = -c1*((K(2,2)*c_tht)+(K(1,2)*s_tht))
      K_nodal(5,3) = 0
      
      K_nodal(5,4) = (c1**2)*((K(2,1)*(c_tht**2))+((K(1,1)-K(2,2))*
     &                c_tht*s_tht)-(K(1,2)*(s_tht**2)))
      K_nodal(5,5) = (c1**2)*((K(2,2)*(c_tht**2))+((K(1,2)+K(2,1))*
     &                c_tht*s_tht)+(K(1,1)*(s_tht**2)))
      K_nodal(5,6) = 0
      
      K_nodal(5,7) = (c1*c2)*((K(2,1)*(c_tht**2))+((K(1,1)-K(2,2))*
     &                c_tht*s_tht)-(K(1,2)*(s_tht**2)))
      K_nodal(5,8) = (c1*c2)*((K(2,2)*(c_tht**2))+((K(1,2)+K(2,1))*
     &                c_tht*s_tht)+(K(1,1)*(s_tht**2)))
      K_nodal(5,9) = 0
      
c     Row : 6      
      K_nodal(6,1) = 0
      K_nodal(6,2) = 0
      K_nodal(6,3) = -c1*K(3,3)
      
      K_nodal(6,4) = 0
      K_nodal(6,5) = 0
      K_nodal(6,6) = (c1**2)*K(3,3)
      
      K_nodal(6,7) = 0
      K_nodal(6,8) = 0
      K_nodal(6,9) = (c1*c2)*K(3,3)
      
c     Row : 7 
      K_nodal(7,1) = -c2*((K(1,1)*c_tht)-(K(2,1)*s_tht))
      K_nodal(7,2) = -c2*((K(1,2)*c_tht)-(K(2,2)*s_tht))
      K_nodal(7,3) = 0
      
      K_nodal(7,4) = (c1*c2)*((K(1,1)*(c_tht**2))-((K(1,2)+K(2,1))*
     &                c_tht*s_tht)+(K(2,2)*(s_tht**2)))
      K_nodal(7,5) = (c1*c2)*((K(1,2)*(c_tht**2))+((K(1,1)-K(2,2))*
     &                c_tht*s_tht)-(K(2,1)*(s_tht**2)))
      K_nodal(7,6) = 0
      
      K_nodal(7,7) = (c2**2)*((K(1,1)*(c_tht**2))-((K(1,2)+K(2,1))*
     &                c_tht*s_tht)+(K(2,2)*(s_tht**2)))
      K_nodal(7,8) = (c2**2)*((K(1,2)*(c_tht**2))+((K(1,1)-K(2,2))*
     &                c_tht*s_tht)-(K(2,1)*(s_tht**2)))
      K_nodal(7,9) = 0
      
c     Row : 8
      K_nodal(8,1) = -c2*((K(2,1)*c_tht)+(K(1,1)*s_tht))
      K_nodal(8,2) = -c2*((K(2,2)*c_tht)+(K(1,2)*s_tht))
      K_nodal(8,3) = 0
      
      K_nodal(8,4) = (c1*c2)*((K(2,1)*(c_tht**2))+((K(1,1)-K(2,2))*
     &                c_tht*s_tht)-(K(1,2)*(s_tht**2)))
      K_nodal(8,5) = (c1*c2)*((K(2,2)*(c_tht**2))+((K(1,2)+K(2,1))*
     &                c_tht*s_tht)+(K(1,1)*(s_tht**2)))
      K_nodal(8,6) = 0
      
      K_nodal(8,7) = (c2**2)*((K(2,1)*(c_tht**2))+((K(1,1)-K(2,2))*
     &                c_tht*s_tht)-(K(1,2)*(s_tht**2)))
      K_nodal(8,8) = (c2**2)*((K(2,2)*(c_tht**2))+((K(1,2)+K(2,1))*
     &                c_tht*s_tht)+(K(1,1)*(s_tht**2)))
      K_nodal(8,9) = 0
      
c     Row : 9      
      K_nodal(9,1) = 0
      K_nodal(9,2) = 0
      K_nodal(9,3) = -c2*K(3,3)
      
      K_nodal(9,4) = 0
      K_nodal(9,5) = 0
      K_nodal(9,6) = (c1*c2)*K(3,3)
      
      K_nodal(9,7) = 0
      K_nodal(9,8) = 0
      K_nodal(9,9) = (c2**2)*K(3,3)
      

      RETURN    
      END    

c======================================================================