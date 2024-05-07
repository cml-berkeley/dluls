c**********************************************************************
c     Subroutines in this file :
c     1. wr_pres
c     2. wr_new_mesh
c     3. wr_grid_recess
c     4. wr_final_data
c     5. write_vec_i
c     6. write_vec_f
c     7. write_mat_i
c     8. write_mat_f
c**********************************************************************

c======================================================================
	SUBROUTINE wr_pres(file_id)
c======================================================================

c     Shared Data
c     ===========      
      USE siml_pram
      USE sldr_grid
      USE rynl_data
      USE sldr_dynm
      USE sldr_stat
      USE luls_data
      
      IMPLICIT NONE
	
	INTEGER file_id
	INTEGER i,j,jj

	IF(igp1.EQ.1) THEN
		DO j=1,ny-1
			jj=jytr(j)
			WRITE(30+file_id,731)((p(ixtr(i),jj)-1.d0),i=1,nx)
		ENDDO
	ELSE
		DO j=1,ny
			WRITE(30+file_id,731)((p(i,j)-1.d0),i=1,nx) 
		ENDDO
	ENDIF
	
	CLOSE(30+file_id)
	
	RETURN

731	FORMAT(1010(E16.9,1X))

	END
	
c======================================================================
	SUBROUTINE wr_new_mesh
c======================================================================

c     Shared Data
c     ===========
      USE siml_pram
      USE grid_ctrl
      USE sldr_grid
      USE sldr_arys
      USE sldr_dynm
      
      IMPLICIT NONE
      
c     Local Variables
c     ===============
      INTEGER i,j,ii
      
c	Write new mesh
c	==============
      OPEN(11,ERR=998,FILE='x.dat',STATUS='UNKNOWN')
	OPEN(12,ERR=998,FILE='y.dat',STATUS='UNKNOWN')

	IF(ioldgrid.EQ.0) THEN
		WRITE(11,'(1500(E16.9,3X))') (xref(i),i=1,nx)
		WRITE(12,'(1500(E16.9,3X))') (yref(j),j=1,ny)
	ENDIF
	
	CLOSE(11)
	CLOSE(12)
      
      OPEN(13,ERR=998,FILE='topog.dat',STATUS='UNKNOWN')
	IF(igp1.EQ.1) THEN
		DO i=1,nx-1
	      ii=ixtr(i)
			WRITE(13,700) (-hm*hxy(ii,jytr(j)),j=1,ny-1)
		ENDDO
	ELSE
		DO i=1,nx
			WRITE(13,700) (-hm*hxy(i,j),j=1,ny)
		ENDDO
	ENDIF
	
700	FORMAT(1010(E16.9,1X))
      CLOSE(13)

      RETURN
      
998	WRITE(*,*) 'Trouble in opening files'
	STOP
      END

c======================================================================
	SUBROUTINE wr_grid_recess
c======================================================================

c----------------------------------------------------------------------
c     Computes the recess h(i,j) at the final grid used in simulation
c----------------------------------------------------------------------

c     Shared Data
c     ===========
      USE sldr_grid
      USE sldr_dynm
      
      IMPLICIT NONE
      
c     Local Variables
c     ===============
      INTEGER i, j
      DOUBLE PRECISION x_loc,y_loc
      DOUBLE PRECISION, DIMENSION (:,:),ALLOCATABLE:: sldr_recess
      
c     Grid point recess from slider refrence surface
c     ==============================================
      ALLOCATE(sldr_recess(nx,ny))
      sldr_recess = 0.d0    ! Initialize      
      DO i = 1, nx
        DO j = 1, ny
            x_loc = xref(i)
            y_loc = yref(j)
            CALL pointrecess(xref(i),yref(j),sldr_recess(i,j))
            sldr_recess(i,j) = sldr_recess(i,j)*hm
        ENDDO
      ENDDO

c     Write Recess Data
c     =================      
      OPEN(14,FILE='x_adapt.dat',STATUS='UNKNOWN')
      OPEN(15,FILE='y_adapt.dat',STATUS='UNKNOWN')
      OPEN(16,FILE='sldr_recess.dat',STATUS='UNKNOWN')
      REWIND(14)
      REWIND(15)
      REWIND(16)
      WRITE(14,'(1500(E16.9,3X))') (xref(i),i=1,nx)
 	WRITE(15,'(1500(E16.9,3X))') (yref(j),j=1,ny)
 	DO j = 1, ny
 	  WRITE(16,'(1500(E16.9,3X))') (sldr_recess(i,j),i=1,nx)
 	ENDDO
 
c     Close files
c     ===========
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      
      DEALLOCATE(sldr_recess)
      
      RETURN
      END

c======================================================================
	SUBROUTINE wr_final_data
c======================================================================

c     Shared Data
c     ===========	
	USE siml_pram
	USE syst_cnfg
	USE rail_data
	USE sldr_grid
	USE sldr_arys
	USE disk_rghn
	USE trck_prfl
	USE rynl_data
	USE sldr_dynm
	USE aspr_data
	USE impc_data
	
	IMPLICIT NONE
	
	INTEGER i,j,ii,jj
	
	OPEN(35,ERR=999,FILE='presf_ABS.dat',STATUS='UNKNOWN')
	OPEN(36,ERR=999,FILE='presf_cont.dat',STATUS='UNKNOWN')
	
	REWIND(35)
	REWIND(36)

	IF(igp1.EQ.1) THEN
		DO j=1,ny-1
		    jj=jytr(j)
			WRITE(35,700)((p(ixtr(i),jj)-1.d0),i=1,nx-1)
			WRITE(36,700)(pcontact(ixtr(i),jj),i=1,nx-1)
		ENDDO
	ELSE
		DO j=1,ny
			WRITE(35,700)((p(i,j)-1.d0),i=1,nx)
			WRITE(36,700)(pcontact(i,j),i=1,nx)
		ENDDO
	ENDIF

	CLOSE(35)
	CLOSE(36)

c	Output Disk Topology under the Current Slider Position
c	======================================================
	IF(nf_wave.NE.0.OR.nf_asper.NE.0.OR.nf_zone.NE.0) THEN
		OPEN(21,file='disktop.dat',status='unknown',recl=5000)
		IF(igp1.EQ.1) THEN
			DO i=1,nx-1
				ii=ixtr(i)
				WRITE(21,700)(-hm*hfloor(ii,jytr(j)),j=1,ny-1)
			ENDDO
		ELSE
		    DO i=1,nx
			    WRITE(21,700)(-hm*hfloor(i,j),j=1,ny)
	      ENDDO
		ENDIF
		CLOSE(21)
	ENDIF


700	FORMAT(1010(E16.9,1X))
	RETURN   
	
999	WRITE(*,*) 'Trouble in opening files'
	STOP	
	END 

c======================================================================
	SUBROUTINE write_vec_i(vec_A,filename,m)
c======================================================================

c----------------------------------------------------------------------
c     Write an integer column vector of length m in a file : <filename>
c----------------------------------------------------------------------

      IMPLICIT NONE
	
	INTEGER i,m
	CHARACTER*70 filename
	INTEGER vec_A(m)

	OPEN(501,FILE=filename,STATUS='UNKNOWN')
	REWIND(501)
	
	DO i = 1, m
	  WRITE(501,'(I6)') vec_A(i)
 	ENDDO

	CLOSE(501)
      
      RETURN
	END
	
c======================================================================
	SUBROUTINE write_vec_f(vec_A,filename,m)
c======================================================================

c----------------------------------------------------------------------
c     Write a real column vector of length m in a file : <filename>
c----------------------------------------------------------------------

      IMPLICIT NONE
	
	INTEGER i,m
	CHARACTER*70 filename
	DOUBLE PRECISION vec_A(m)

	OPEN(502,FILE=filename,STATUS='UNKNOWN')
	REWIND(502)
	
	DO i = 1, m
	  WRITE(502,'(E16.9)') vec_A(i)
 	ENDDO

	CLOSE(502)
      
      RETURN
	END

c======================================================================
	SUBROUTINE write_mat_i(mat_A,filename,m,n)
c======================================================================

c----------------------------------------------------------------------
c     Write an integer matrix of size m x n in a file : <filename>
c----------------------------------------------------------------------

	IMPLICIT NONE
	
	INTEGER i,j,m,n
	CHARACTER*70 filename
	INTEGER mat_A(m,n)

	OPEN(503,FILE=filename,STATUS='UNKNOWN')
	REWIND(503)
	
	DO i = 1, m
	  WRITE(503,'(25000(I6,3X))') 
     &		(mat_A(i,j), j=1,n)
	ENDDO
	
	CLOSE(503)
      
      RETURN
	END
	
c======================================================================	
	SUBROUTINE write_mat_f(mat_A,filename,m,n)
c======================================================================

c----------------------------------------------------------------------
c     Write a real matrix of size m x n in a file : <filename>
c----------------------------------------------------------------------

	IMPLICIT NONE
	
	INTEGER i,j,m,n
	CHARACTER*70 filename
	DOUBLE PRECISION mat_A(m,n)

	OPEN(504,FILE=filename,STATUS='UNKNOWN')
	REWIND(504)
	
	DO i = 1, m
	  WRITE(504,'(25000(E16.9,3X))') 
     &		(mat_A(i,j), j=1,n)
	ENDDO
	
	CLOSE(504)
      
      RETURN
	END
	
c======================================================================