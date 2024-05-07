c**********************************************************************
c     Subroutines in this file :
c      1. qromb
c      2. trapzd
c      3. polint
c      4. finteg
c      5. finteg1
c      6. locate
c      7. matrix_multi
c      8. matrix33_inverse
c      9. tridag
c     10. sort
c     11. mult
c     12. inverse
c     13. gauss_quad
c     14. d_mult_sparse_mat_vec
c     15. sparse_mat_map
c     16. LU_fact
c     17. LU_sparse_fact
c     18. LU_solve
c     19. L_solve
c     20. U_solve
c     21. LU_solve_sparse
c     22. L_solve_sparse
c     23. U_solve_sparse
c     24. U_T_solve_sparse
c     25. L_solve_mult
c     26. U_T_solve_mult
c**********************************************************************

c======================================================================
	FUNCTION fit1(x)
c======================================================================

	IMPLICIT NONE
	
	DOUBLE PRECISION s,sc,x,fit1
	COMMON/var/s,sc
	
	fit1=(x-s)*EXP(-0.50d0*x*x)
	
	RETURN
	END

c======================================================================
	FUNCTION fitb1(x)
c======================================================================

	IMPLICIT NONE
	
	DOUBLE PRECISION s,sc,x,fitb1
	COMMON/var/s,sc
	
	fitb1=(2.d0*(x-s)-sc)*EXP(-0.50d0*x*x)
	
	RETURN
	END

c======================================================================
	FUNCTION fit32(x)
c======================================================================

	IMPLICIT NONE
	
	DOUBLE PRECISION s,sc,x,fit32
	COMMON/var/s,sc
	
	fit32=(x-s)*DSQRT(x-s)*EXP(-0.50d0*x*x)
	
	RETURN
	END

c======================================================================
	SUBROUTINE qromb(func,a,b,ss)
c====================================================================== 
	IMPLICIT REAL*8(a-h,o-z)
	EXTERNAL func
	PARAMETER (eps=1.d-6,jmax=1,jmaxp=jmax+1,k=5,km=k-1)
	DIMENSION s(jmaxp),h(jmaxp)
	h(1)=1.d0
	DO j=1,jmax
	   CALL trapzd(func,a,b,s(j),j)
	   IF(j.GE.k) THEN
			CALL polint(h(j-km),s(j-km),k,0.d0,ss,dss)
			IF(DABS(dss).LT.eps*DABS(ss)) RETURN
	   ENDIF
	   s(j+1)=s(j)
	   h(j+1)=0.25d0*h(j)
	ENDDO
	RETURN
	
	WRITE(*,*) 'Too many steps in qromb'
	STOP
	
	END

c======================================================================
	SUBROUTINE trapzd(func,a,b,s,n)
c======================================================================

	IMPLICIT NONE
	
	INTEGER j,it,n
	DOUBLE PRECISION func,a,b,s,tnm,del,x,sum
	
	EXTERNAL func
	
	IF(n.EQ.1) THEN
		s=0.5d0*(b-a)*(func(a)+func(b))
		it=1
	ELSE
		tnm = DBLE(it)
		del = (b-a)/tnm
		x = a+0.5d0*del
		sum = 0.0d0
		DO j = 1,it
			sum=sum+func(x)
			x=x+del
		ENDDO
		s=0.5d0*(s+(b-a)*sum/tnm)
		it=2*it
	ENDIF

	RETURN
	END
	 
c======================================================================	
	SUBROUTINE polint(xa,ya,n,x,y,dy)
c======================================================================

	IMPLICIT REAL*8(a-h,o-z)
	PARAMETER (nmax=10)
	DIMENSION xa(n),ya(n),c(nmax),d(nmax)
	
	ns=1
	dif=DABS(x-xa(1))
	
	DO i=1,n
		dift=DABS(x-xa(i))
		IF(dift.LT.dif) THEN
			ns=i
			dif=dift
		ENDIF
		c(i)=ya(i)
		d(i)=ya(i)
  	ENDDO

	y=ya(ns)
	ns=ns-1
	
	DO m=1,n-1
		DO i=1,n-m
			ho=xa(i)-x
			hp=xa(i+m)-x
			w=c(i+1)-d(i)
			den=ho-hp
		 
			IF(den.EQ.0.d0) PAUSE
			den=w/den
			d(i)=hp*den
			c(i)=ho*den
 		ENDDO
	   
		IF(2*ns.LT.n-m) THEN
			dy=c(ns+1)
		ELSE
			dy=d(ns)
			ns=ns-1
		ENDIF
		y=y+dy
	ENDDO
	
	RETURN
	END
	 
c======================================================================
	SUBROUTINE finteg(tgs,fgs1,fgs2,n,x,y1,y2)
c======================================================================

c----------------------------------------------------------------------
c	Integral calculation for the Contact models
c----------------------------------------------------------------------

	IMPLICIT REAL*8(a-h,o-z)												 
	 DIMENSION tgs(n),fgs1(n),fgs2(n)

	CALL locate(tgs,n,x,j)
	
	a=(x-tgs(j+1))/(tgs(j)-tgs(j+1))
	b=(x-tgs(j))/(tgs(j+1)-tgs(j))
	y1=a*fgs1(j)+b*fgs1(j+1)
	y2=a*fgs2(j)+b*fgs2(j+1)
	
	RETURN
	END

c======================================================================
	SUBROUTINE finteg1(tgs,fgs1,n,x,y)
c======================================================================	  
	
	IMPLICIT REAL*8(a-h,o-z)
	DIMENSION tgs(n),fgs1(n)

	IF(x.LE.tgs(1)) THEN
		y=fgs1(1)
	ELSE IF(x.GE.tgs(n)) THEN
		y=fgs1(n)
	ELSE
		CALL locate(tgs,n,x,j)
		a=(x-tgs(j+1))/(tgs(j)-tgs(j+1))
		b=(x-tgs(j))/(tgs(j+1)-tgs(j))
		y=a*fgs1(j)+b*fgs1(j+1)
	ENDIF
	
	RETURN
	END

c======================================================================	
	SUBROUTINE locate(xx,n,x,j)
c======================================================================	  
	 
	IMPLICIT REAL*8(a-h,o-z)
	DIMENSION xx(n)												 

	jl=1
	ju=n

10	IF(ju-jl.GT.1) THEN
		jm=(ju+jl)/2
		IF((xx(n).GT.xx(1)).eqv.(x.GT.xx(jm)))THEN
			jl=jm
		ELSE
			ju=jm
		ENDIF
		GO TO 10
	ENDIF
      j=jl

	RETURN
	END

c======================================================================
      REAL*8 FUNCTION flog(fxvar)
c======================================================================

	IMPLICIT REAL*8(a-h,o-z) 

      f_tmp1 = 0.37252903d-8/fxvar
         
      f_an = (f_tmp1+1.d0)/2.d0
      f_gn = SQRT(f_tmp1*1.d0)
         
	i_cn = 0
    
      DO WHILE ((ABS(f_an-f_gn).GT.1d-8).AND.(i_cn.LT.25))
		f_anp1 = (f_an+f_gn)/2
          f_gnp1 = SQRT(f_an*f_gn)
          f_an = f_anp1
          f_gn = f_gnp1
	    i_cn = i_cn + 1
      ENDDO
         
      flog = 1.57079633/f_an - 20.79441542

      RETURN
	END
	
c======================================================================
	SUBROUTINE matrix_multi(n1,n2,n3,a,b,c)
c======================================================================
	IMPLICIT REAL*8(a-h,o-z)
	DIMENSION a(n1,n2),b(n2,n3),c(n1,n3)
	
	DO i=1,n1
		DO j=1,n3
			c(i,j)=0.d0
			DO k=1,n2
				c(i,j)=c(i,j)+ a(i,k)*b(k,j)
			ENDDO
		ENDDO
	ENDDO
	
	RETURN
	END

c======================================================================
	SUBROUTINE matrix33_inverse(a)
c======================================================================

c----------------------------------------------------------------------
c	Computes the inverse for 3 x 3 matrix 
c----------------------------------------------------------------------

	IMPLICIT REAL*8(a-h,o-z)
	DIMENSION a(3,3)
	  
	deter=a(1,1)*a(2,2)*a(3,3) + a(1,2)*a(2,3)*a(3,1) +
     &		a(2,1)*a(3,2)*a(1,3) - a(1,3)*a(2,2)*a(3,1) -
     &		a(1,2)*a(2,1)*a(3,3) - a(3,2)*a(2,3)*a(1,1)
	  
	ai11 =  a(2,2)*a(3,3)-a(2,3)*a(3,2)
	ai12 =-(a(2,1)*a(3,3)-a(2,3)*a(3,1))
	ai13 =  a(2,1)*a(3,2)-a(3,1)*a(2,2)
	ai21 =-(a(1,2)*a(3,3)-a(1,3)*a(3,2))
	ai22 =  a(1,1)*a(3,3)-a(1,3)*a(3,1)
	ai23 =-(a(1,1)*a(3,2)-a(1,2)*a(3,1))
	ai31 =  a(1,2)*a(2,3)-a(1,3)*a(2,2)
	ai32 =-(a(1,1)*a(2,3)-a(1,3)*a(2,1))
	ai33 =  a(1,1)*a(2,2)-a(1,2)*a(2,1)
	  
	a(1,1) = ai11/deter
	a(1,2) = ai21/deter
	a(1,3) = ai31/deter
	a(2,1) = ai12/deter
	a(2,2) = ai22/deter
	a(2,3) = ai32/deter
	a(3,1) = ai13/deter
	a(3,2) = ai23/deter
	a(3,3) = ai33/deter

	RETURN
	END

c======================================================================
	 SUBROUTINE tridag(n,a,b,c,d)											 
c======================================================================	 

c----------------------------------------------------------------------
c     This subroutine inverts a tri-diagonal matrix
c----------------------------------------------------------------------
	
	IMPLICIT REAL*8(a-h,o-z)												 
	DIMENSION a(1),b(1),c(1),d(1) 										 
	DIMENSION beta(n),gamma(n)
																			 
	beta(1)=b(1)															 
	gamma(1)=d(1)/beta(1) 												 
	
	DO i = 2 , n															 
		im1=i-1														 
		beta(i)=b(i)-a(i)*c(im1)/beta(im1)							 
		gamma(i)=(d(i)-a(i)*gamma(im1))/beta(i)
	ENDDO
																		 
	d(n)=gamma(n) 														 
	
	DO i=n-1,1,-1														 
		d(i)=gamma(i)-c(i)*d(i+1)/beta(i) 							 
	ENDDO																 
																			 
	RETURN																 
	END
	
c======================================================================
	SUBROUTINE sort(eps,x,ix)
c======================================================================
	
	IMPLICIT REAL*8(a-h,o-z)
	DIMENSION x(ix)

	DO k=1,ix-1
		DO i=1,ix-k
			IF(x(i).GT.x(i+1)) THEN
				xtp=x(i)
				x(i)=x(i+1)
				x(i+1)=xtp
			ENDIF
		ENDDO
	ENDDO

	icount=1

	DO i=2,ix
		IF(DABS(x(i)-x(i-1)).GT.eps) THEN
			icount=icount+1
			x(icount)=x(i)
		ENDIF
	ENDDO

	ix=icount

	RETURN
	END
		
c======================================================================
      SUBROUTINE mult(f_mat1,f_mat2,N,NMAXSIZE,M,MMAXSIZE,f_res)
c======================================================================
	
	IMPLICIT REAL*8(a-h,o-z)

	DIMENSION f_mat1(NMAXSIZE,NMAXSIZE), f_mat2(NMAXSIZE,MMAXSIZE),
     &	f_res(NMAXSIZE,MMAXSIZE)

	DO i=1,N
		DO j=1,M
			f_temp = 0.d0
	        DO k=1,N
				f_temp = f_temp+f_mat1(i,k)*f_mat2(k,j)
	        ENDDO
	        f_res(i,j) = f_temp
	    ENDDO
	ENDDO
	END

c======================================================================
      SUBROUTINE inverse(f_mata,N,NMAXSIZE,f_inv)
c======================================================================
	
	IMPLICIT REAL*8(a-h,o-z)
	DIMENSION f_mata(NMAXSIZE,NMAXSIZE),f_inv(NMAXSIZE,NMAXSIZE)

      INTEGER,DIMENSION(:),  ALLOCATABLE::i_bclist
      REAL*8, DIMENSION(:,:),ALLOCATABLE::f_mat

      ALLOCATE(i_bclist(NMAXSIZE))
      ALLOCATE(f_mat(NMAXSIZE,NMAXSIZE))

      DO i=1,N
		DO j=1,N
			f_mat(i,j) = f_mata(i,j)
              IF (i.EQ.j) THEN
				f_inv(i,j) = 1.d0
	        ELSE
				f_inv(i,j) = 0.d0
	        ENDIF
		ENDDO
	ENDDO

      DO i=1,N
          f_big = 0.d0
		DO j=i,N
              DO k=1,N
				IF (ABS(f_mat(j,k)).GT.f_big) THEN
					f_big = f_mat(j,k)
					i_br = j
					i_bc = k
				ENDIF
			ENDDO
          ENDDO
    
          i_bclist(i) = i_bc
    
		IF (i_br.NE.i) THEN
              DO j=1,N
                  f_temp = f_mat(i,j)
                  f_mat(i,j) = f_mat(i_br,j)
                  f_mat(i_br,j) = f_temp

                  f_temp = f_inv(i,j)
                  f_inv(i,j) = f_inv(i_br,j)
				f_inv(i_br,j) = f_temp
			ENDDO
          ENDIF
          
		IF (i_bc.NE.i) THEN
              DO j=1,N
                  f_temp = f_mat(j,i)
                  f_mat(j,i) = f_mat(j,i_bc)
				f_mat(j,i_bc) = f_temp
			ENDDO
          ENDIF
    
          f_fac = f_mat(i,i)
          
		DO j=1,N
              f_mat(i,j) = f_mat(i,j)/f_fac

			f_inv(i,j) = f_inv(i,j)/f_fac
          ENDDO
    
		DO j=1,N
			IF (j.NE.i) THEN
                  f_fac = f_mat(j,i)
                  DO k=1,N
					f_mat(j,k) = f_mat(j,k) - f_fac*f_mat(i,k)
					f_inv(j,k) = f_inv(j,k) - f_fac*f_inv(i,k)
				ENDDO
			ENDIF
		ENDDO
    
      ENDDO

      DO i=N,1,-1
		IF (i_bclist(i).NE.i) THEN
              DO j=1,N
                  f_temp = f_inv(i,j)
                  f_inv(i,j) = f_inv(i_bclist(i),j)
				f_inv(i_bclist(i),j) = f_temp
			ENDDO
		ENDIF
      ENDDO

      DEALLOCATE(i_bclist)
      DEALLOCATE(f_mat)
      
      RETURN
      END  


c=====================================================================
	SUBROUTINE gauss_quad(gauss_od,gauss_wt,gauss_pt)
c=====================================================================

	INTEGER gauss_od
	DOUBLE PRECISION gauss_wt(gauss_od), gauss_pt(gauss_od)

	IF (gauss_od.EQ.2) THEN
		gauss_wt = (/  1.0d0,             1.0d0/)
		gauss_pt = (/ -0.57735026918963d0,0.57735026918963d0/)
	ELSE IF (gauss_od.EQ.3) THEN
		gauss_wt = (/0.55555555d0, 0.88888889d0,0.55555555d0/)
		gauss_pt = (/-0.77459667d0,0.0d0,       0.77459667d0/)
	ELSE IF (gauss_od.EQ.4) THEN
		gauss_wt = (/ 0.34785485d0, 0.65214515d0,
     &	              0.65214515d0, 0.34785485d0/)
		gauss_pt = (/-0.86113631d0,-0.33998104d0,
     &				  0.33998104d0,0.86113631d0/)
	ELSE IF (gauss_od.EQ.5) THEN
		gauss_wt = (/ 0.23692689d0, 0.47862867d0,
     &				  0.56888889d0, 0.47862867d0, 0.23692689d0/)
		gauss_pt = (/-0.90617985d0,-0.53846931d0,
     &	              0.0d0       , 0.53846931d0, 0.90617985d0/)
	ELSE
		WRITE(*,*) 'Gauss order value out of range'
		STOP
	
	ENDIF
      
      RETURN
	END

c======================================================================
	SUBROUTINE d_mult_sparse_mat_vec(mat,sparse_map,vec_in,vec_out,
     &	nx_mat,ny_mat,n_vec)
c======================================================================

c----------------------------------------------------------------------
c	Multiply Double Precision Sparse Matrices with vector
c	=====================================================
c
c	Input: 
c	Matrix     : (nx_mat X ny_mat)
c	Sparse Map :  (nx_mat X ny_mat + 1)
c	Vector     : n_vec
c
c	Output: 
c	Vector : (nx_mat)
c
c	Error:
c	ny_mat ~= n_vec
c----------------------------------------------------------------------

      IMPLICIT NONE

	INTEGER i,j,k,col_no
	INTEGER nx_mat,ny_mat,n_vec
	INTEGER sparse_map(nx_mat,ny_mat+1)
	DOUBLE PRECISION sum_val
	DOUBLE PRECISION mat(nx_mat,ny_mat)
	DOUBLE PRECISION vec_in(n_vec), vec_out(nx_mat)

	IF (ny_mat.NE.n_vec) THEN 
		WRITE(*,*) 'Input matrix size mistmatch'
		STOP
	ELSE
		DO i = 1 , nx_mat
			sum_val = 0.0d0
			DO j = 1 , sparse_map(i,1)
				col_no = sparse_map(i,j+1)
				sum_val = sum_val + (mat(i,col_no)*vec_in(col_no))
			ENDDO
			vec_out(i) = sum_val
		ENDDO
	ENDIF
	
	RETURN
	END

c======================================================================
	SUBROUTINE sparse_mat_map(A,B,m,n)
c======================================================================

c---------------------------------------------------------------------
c	Find Non zero entries of sparse matrix
c	======================================
c
c	Input: 
c	Matrix A : m X n
c	
c
c	Output: 
c	Matrix B : m X (n + 1)
c	Non zero entries : cnt	

c	Note: Work with DOUBLE PRECISION data type only
c		  : Convert into DOUBLE PRECISION A = DBLE(B)
c---------------------------------------------------------------------
      
      IMPLICIT NONE
      
	INTEGER i,j
	INTEGER m,n,cnt,col_cnt
	INTEGER B(m,n+1)
	DOUBLE PRECISION A(m,n)

	B = 0
	DO i = 1 , m
		col_cnt = 0
		DO j = 1 , n
			IF (A(i,j) .NE. 0.0d0) THEN
				col_cnt = col_cnt + 1
				B(i,col_cnt+1) = j
			ENDIF
		ENDDO
		B(i,1) = col_cnt
	ENDDO
	
	RETURN
	END

c======================================================================
	SUBROUTINE LU_fact(A,L,U,n)
c======================================================================

c	LU Factorization : A = LU
c	-------------------------

	INTEGER i,j,k,n
	DOUBLE PRECISION sum
	DOUBLE PRECISION A(n,n), L(n,n), U(n,n)

c	Initialization
c`	--------------
	L = 0.0d0
	U = 0.0d0

	DO i = 1 , n
		L(i,i) = 1.0d0
	ENDDO

c	Doolittle Algorithm
c	-------------------

c	**********************************************
c	First and Second row require special treatment
c	to avoid error DO loop
c	**********************************************

c	First Row
c	---------
	DO j = 1 , n
		U(1,j) = A(1,j)
	ENDDO
	
c	Second Row
c	----------
	L(2,1) = a(2,1)/U(1,1)
	DO j = 1 , n
		U(2,j) = A(2,j) - (L(2,1)*U(1,j))
	ENDDO

c	Third Row Onward
c	----------------			
	DO i = 3, n
		L(i,1) = A(i,1)/U(1,1)
		
		DO j = 2 , i - 1
			sum = 0.0d0
			DO k = 1 , j - 1
				sum = sum + (L(i,k)*U(k,j))
			ENDDO
			L(i,j) = (A(i,j) - sum)/U(j,j)
		ENDDO
		
		DO j = i , n
			sum = 0.0d0
			DO k = 1 , i - 1
				sum = sum + (L(i,k)*U(k,j))
			ENDDO
			U(i,j) = A(i,j) - sum
		ENDDO

	ENDDO

	RETURN
	END

c======================================================================
	SUBROUTINE LU_sparse_fact(A,L,U,L_map,U_map,n)
c======================================================================
      
      IMPLICIT NONE

c	LU Factorization : A = LU
c	-------------------------

	INTEGER i,j,k,n
	INTEGER col_cnt_L,col_cnt_U
	DOUBLE PRECISION sum
	INTEGER	L_map(n,n+1),U_map(n,n+1)
	DOUBLE PRECISION A(n,n), L(n,n), U(n,n)

c	Initialization
c`	--------------
	L = 0.0d0
	U = 0.0d0

	DO i = 1 , n
		L(i,i) = 1.0d0
	ENDDO

c	Doolittle Algorithm
c	-------------------

c	**********************************************
c	First and Second row require special treatment
c	to avoid error DO loop
c	**********************************************

c	First Row
c	---------
	DO j = 1 , n
		U(1,j) = A(1,j)
	ENDDO
	
c	Second Row
c	----------
	L(2,1) = a(2,1)/U(1,1)
	DO j = 1 , n
		U(2,j) = A(2,j) - (L(2,1)*U(1,j))
	ENDDO

c	Third Row Onward
c	----------------
	DO i = 3, n
		L(i,1) = A(i,1)/U(1,1)
		
		DO j = 2 , i - 1
			sum = 0.0d0
			DO k = 1 , j - 1
				sum = sum + (L(i,k)*U(k,j))
			ENDDO
			L(i,j) = (A(i,j) - sum)/U(j,j)
		ENDDO
		
		DO j = i , n
			sum = 0.0d0
			DO k = 1 , i - 1
				sum = sum + (L(i,k)*U(k,j))
			ENDDO
			U(i,j) = A(i,j) - sum
		ENDDO

	ENDDO

c	Sparse Map
c	----------
	L_map = 0
	U_map = 0
	DO i = 1 , n
		col_cnt_L = 0
		col_cnt_U = 0
		DO j = i , n
			IF (U(i,j) .NE. 0.0d0) THEN
				col_cnt_U = col_cnt_U + 1
				U_map(i,col_cnt_U+1) = j
			ENDIF
			IF (L(n-i+1,j-i+1) .NE. 0.0d0) THEN
				col_cnt_L = col_cnt_L + 1
				L_map(n-i+1,col_cnt_L+1) = j-i+1
			ENDIF

		ENDDO
		L_map(n-i+1,1) = col_cnt_L
		U_map(i,1)     = col_cnt_U
	ENDDO
       
      RETURN
	END

c======================================================================
	SUBROUTINE LU_solve(L,U,b,x,n)
c======================================================================

c	LU Factorization : A = LU
c	-------------------------

	INTEGER i,j,n
	DOUBLE PRECISION b(n),x(n),y(n)
	DOUBLE PRECISION L(n,n), U(n,n)

c	Forward Substitution Ly = b
c	===========================
	y = b
	y(1) = y(1)/L(1,1)
	DO i = 2, n
		DO j = 1, i-1
			y(i) = y(i) - (L(i,j)*y(j))
		ENDDO
		y(i) = y(i)/L(i,i)
	ENDDO

c	Backward Substitution Ux = y
c	============================
	x = y
	x(n) = x(n)/U(n,n)
	DO i = n - 1 , 1, -1
		DO j = i+1 , n
			x(i) = x(i) - (U(i,j)*x(j))
		ENDDO
		x(i) = x(i)/U(i,i)
	ENDDO
      
      RETURN
	END

c======================================================================
	SUBROUTINE L_solve(L,b,y,n)
c======================================================================

c	Lower triangle matrix: L
c	------------------------

	INTEGER i,j,n
	DOUBLE PRECISION b(n), y(n)
	DOUBLE PRECISION L(n,n)

c	Forward Substitution Ly = b
c	===========================
	y = b
	y(1) = y(1)/L(1,1)
	DO i = 2, n
		DO j = 1, i-1
			y(i) = y(i) - (L(i,j)*y(j))
		ENDDO
		y(i) = y(i)/L(i,i)
	ENDDO

      RETURN
	END

c======================================================================
	SUBROUTINE U_solve(U,b,x,n)
c======================================================================

c	Upper triangle matrix: U
c	------------------------

	INTEGER i,j,n
	DOUBLE PRECISION b(n), x(n)
	DOUBLE PRECISION U(n,n)

c	Backward Substitution Ux = b
c	============================
	x = b
	x(n) = x(n)/U(n,n)
	DO i = n - 1 , 1, -1
		DO j = i+1 , n
			x(i) = x(i) - (U(i,j)*x(j))
		ENDDO
		x(i) = x(i)/U(i,i)
	ENDDO

      RETURN
	END

c======================================================================
	SUBROUTINE L_solve_sparse(L,L_map,b,y,n)
c======================================================================

      IMPLICIT NONE
      
c	Lower triangle matrix: L
c	------------------------
	INTEGER i,j,n,col_no
	INTEGER L_map(n,n+1)
	DOUBLE PRECISION b(n), y(n)
	DOUBLE PRECISION L(n,n)

c	Forward Substitution Ly = b
c	===========================
	y = b
	y(1) = y(1)/L(1,1)
	DO i = 2, n
		DO j = 1, (L_map(i,1)-1)
			col_no = L_map(i,j+1)
			y(i) = y(i) - (L(i,col_no )*y(col_no))
		ENDDO
		y(i) = y(i)/L(i,i)
	ENDDO
      
      RETURN
	END

c======================================================================
	SUBROUTINE U_solve_sparse(U,U_map,b,x,n)
c======================================================================


c	Upper triangle matrix: U
c	------------------------
	INTEGER i,j,n,col_no
	INTEGER U_map(n,n+1)
	DOUBLE PRECISION b(n), x(n)
	DOUBLE PRECISION U(n,n)

c	Backward Substitution Ux = b
c	============================
	x = b
	x(n) = x(n)/U(n,n)
	DO i = n - 1 , 1, -1
		DO j = 2 , (U_map(i,1))
			col_no = U_map(i,j+1)
			x(i) = x(i) - (U(i,col_no)*x(col_no))
		ENDDO
		x(i) = x(i)/U(i,i)
	ENDDO
      
      RETURN
	END
	
c======================================================================
	SUBROUTINE U_T_solve_sparse(U,U_T_map,b,y,n)
c======================================================================

      IMPLICIT NONE
      
c	Local Variables
c	===============
	INTEGER i,j,n,col_no
	INTEGER U_T_map(n,n+1)
	DOUBLE PRECISION b(n), y(n)
	DOUBLE PRECISION U(n,n)

c	Forward Substitution U'y' = b'
c	==============================
	y = b
	y(1) = y(1)/U(1,1)
	DO i = 2, n
		DO j = 1, (U_T_map(i,1)-1)
			col_no = U_T_map(i,j+1)
			y(i) = y(i) - (U(col_no,i)*y(col_no))
		ENDDO
		y(i) = y(i)/U(i,i)
	ENDDO
      
      RETURN
	END

c======================================================================
	SUBROUTINE LU_solve_sparse(L,U,L_map,U_map,b,x,n)
c======================================================================

      IMPLICIT NONE

c	LU Factorization : A = LU
c	-------------------------
	INTEGER i,j,n,col_no
	INTEGER L_map(n,n+1),U_map(n,n+1)
	DOUBLE PRECISION b(n),x(n),y(n)
	DOUBLE PRECISION L(n,n), U(n,n)

	x = 0.0d0

c	Forward Substitution Ly = b
c	===========================
	y = b
	y(1) = y(1)/L(1,1)
	DO i = 2, n
		DO j = 1, (L_map(i,1)-1)
			col_no = L_map(i,j+1)
			y(i) = y(i) - (L(i,col_no )*y(col_no))
		ENDDO
		y(i) = y(i)/L(i,i)
	ENDDO

c	Backward Substitution Ux = y
c	============================
	x = y
	x(n) = x(n)/U(n,n)
	DO i = n - 1 , 1, -1
		DO j = 2 , (U_map(i,1))
			col_no = U_map(i,j+1)
			x(i) = x(i) - (U(i,col_no)*x(col_no))
		ENDDO
		x(i) = x(i)/U(i,i)
	ENDDO
	
      RETURN
	END


c======================================================================
	SUBROUTINE L_solve_mult(L,b,x,m,n)
c======================================================================

c----------------------------------------------------------------------
c	Solve multiple (m) forward substitution problem at a time
c	L: Lower triangular matrix
c	b: m known column vectors (size: n)
c	x: m unknown column vectors (size: n)
c----------------------------------------------------------------------

c	Lower triangle matrix: L
c	------------------------

	INTEGER i,j,m,n,i_case
	DOUBLE PRECISION x(n,m), b(n,m)
	DOUBLE PRECISION L(n,n)

	DO i_case = 1 , m
c		Forward Substitution L*b(i) = x(i)
c		==================================
		x(1:n,i_case) = b(1:n,i_case)
		x(1,i_case) = x(1,i_case)/L(1,1)
		DO i = 2, n
			DO j = 1, i-1
				x(i,i_case) = x(i,i_case) - (L(i,j)*x(j,i_case))
			ENDDO
			x(i,i_case) = x(i,i_case)/L(i,i)
		ENDDO
	ENDDO
      
      RETURN
	END

c======================================================================
	SUBROUTINE U_T_solve_mult(U,b,x,m,n)
c======================================================================

c----------------------------------------------------------------------
c	Solve multiple (m) forward substitution problem : xU = b => b'=U'x'
c	U: Upper triangular matrix => U' = Lower triangular matrix
c	b: m known   row vectors (size: n) => b' : m known   column vectors
c	x: m unknown row vectors (size: n) => x' : m unknown column vectors
c----------------------------------------------------------------------

c     Local Variables
c     ================
	INTEGER i,j,m,n,i_case
	DOUBLE PRECISION x(m,n), b(m,n)
	DOUBLE PRECISION U(n,n)

	DO i_case = 1 , m
	
c		Forward Substitution U'*b'(i_case) = x'(i_case)

		x(i_case,1:n) = b(i_case,1:n)
		x(i_case,1) = x(i_case,1)/U(1,1)
		DO i = 2, n
			DO j = 1, i-1
				x(i_case,i) = x(i_case,i) - (U(j,i)*x(i_case,j))
			ENDDO
			x(i_case,i) = x(i_case,i)/U(i,i)
		ENDDO
	ENDDO
      
      RETURN
	END	
	
c======================================================================	
