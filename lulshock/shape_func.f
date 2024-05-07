c**********************************************************************
c     Subroutines in this file :
c     1. N_annular
c     2. DNr_annular
c     3. DNt_annular
c     4. DNrr_annular
c     5. DNtt_annular
c     6. DNrt_annular
c**********************************************************************

c======================================================================
	SUBROUTINE N_annular(a,b,drad,dtht,y)
c======================================================================

c----------------------------------------------------------------------

c	Shape function for the annular element
c	--------------------------------------
c	a = zeta
c	b = eta
c	drad = Delta radius
c	dtht = Delta theta
c----------------------------------------------------------------------

      IMPLICIT NONE

	DOUBLE PRECISION a,b,drad,dtht
	DOUBLE PRECISION y(16)

	y(1)  =  (1.0d0/16.0d0)*((1.0d0-a)**2)*
     &((1.0d0-b)**2)*(2.0d0+a)*(2.0d0+b)
	y(2)  =  (1.0d0/16.0d0)*((1.0d0+a)**2)*
     &((1.0d0-b)**2)*(2.0d0-a)*(2.0d0+b)
	y(3)  =  (1.0d0/16.0d0)*((1.0d0+a)**2)*
     &((1.0d0+b)**2)*(2.0d0-a)*(2.0d0-b)
	y(4)  =  (1.0d0/16.0d0)*((1.0d0-a)**2)*
     &((1.0d0+b)**2)*(2.0d0+a)*(2.0d0-b)
	y(5)  =  (drad/2.0d0)*(1.0d0/16)*((1.0d0-a)**2)*
     &((1.0d0-b)**2)*(1.0d0+a)*(2.0d0+b)
	y(6)  = -(drad/2.0d0)*(1.0d0/16)*((1.0d0+a)**2)*
     &((1.0d0-b)**2)*(1.0d0-a)*(2.0d0+b)
	y(7)  = -(drad/2.0d0)*(1.0d0/16)*((1.0d0+a)**2)*
     &((1.0d0+b)**2)*(1.0d0-a)*(2.0d0-b)
	y(8)  =  (drad/2.0d0)*(1.0d0/16)*((1.0d0-a)**2)*
     &((1.0d0+b)**2)*(1.0d0+a)*(2.0d0-b)
	y(9)  =  (dtht/2.0d0)*(1.0d0/16)*((1.0d0-a)**2)*
     &((1.0d0-b)**2)*(2.0d0+a)*(1.0d0+b)
	y(10) =  (dtht/2.0d0)*(1.0d0/16)*((1.0d0+a)**2)*
     &((1.0d0-b)**2)*(2.0d0-a)*(1.0d0+b)
	y(11) = -(dtht/2.0d0)*(1.0d0/16)*((1.0d0+a)**2)*
     &((1.0d0+b)**2)*(2.0d0-a)*(1.0d0-b)
	y(12) = -(dtht/2.0d0)*(1.0d0/16)*((1.0d0-a)**2)*
     &((1.0d0+b)**2)*(2.0d0+a)*(1.0d0-b)
	y(13) =  (drad*dtht/4.0d0)*(1.0d0/16)*((1.0d0-a)**2)*
     &((1.0d0-b)**2)*(1.0d0+a)*(1.0d0+b)
	y(14) = -(drad*dtht/4.0d0)*(1.0d0/16)*((1.0d0+a)**2)*
     &((1.0d0-b)**2)*(1.0d0-a)*(1.0d0+b)
	y(15) =  (drad*dtht/4.0d0)*(1.0d0/16)*((1.0d0+a)**2)*
     &((1.0d0+b)**2)*(1.0d0-a)*(1.0d0-b)
	y(16) = -(drad*dtht/4.0d0)*(1.0d0/16)*((1.0d0-a)**2)*
     &((1.0d0+b)**2)*(1.0d0+a)*(1.0d0-b)

      RETURN
	END

c======================================================================
	SUBROUTINE DNr_annular(a,b,drad,dtht,y)
c======================================================================

c----------------------------------------------------------------------
c	Derivative of Shape function w.r.t radius
c	-----------------------------------------
c	a = zeta
c	b = eta
c	drad = Delta radius
c	dtht = Delta theta
c----------------------------------------------------------------------
      
      IMPLICIT NONE

	DOUBLE PRECISION a,b,drad,dtht
	DOUBLE PRECISION y(16)

	y(1)  = -(3.0d0/16.0d0)*(1.0d0-(a**2))*
     &((1.0d0-b)**2)*(2.0d0+b)
	y(2)  =  (3.0d0/16.0d0)*(1.0d0-(a**2))*
     &((1.0d0-b)**2)*(2.0d0+b)
	y(3)  =  (3.0d0/16.0d0)*(1.0d0-(a**2))*
     &((1.0d0+b)**2)*(2.0d0-b)
	y(4)  = -(3.0d0/16.0d0)*(1.0d0-(a**2))*
     &((1.0d0+b)**2)*(2.0d0-b)
	y(5)  = -(drad/32.0d0)*(1.0d0-a)*(1.0d0+(3.0d0*a))*
     &((1.0d0-b)**2)*(2.0d0+b)
	y(6)  = -(drad/32.0d0)*(1.0d0+a)*(1.0d0-(3.0d0*a))*
     &((1.0d0-b)**2)*(2.0d0+b)
	y(7)  = -(drad/32.0d0)*(1.0d0+a)*(1.0d0-(3.0d0*a))*
     &((1.0d0+b)**2)*(2.0d0-b)
	y(8)  = -(drad/32.0d0)*(1.0d0-a)*(1.0d0+(3.0d0*a))*
     &((1.0d0+b)**2)*(2.0d0-b)
	y(9)  = -(3.0d0*dtht/32.0d0)*(1.0d0-(a**2))*
     &((1.0d0-b)**2)*(1.0d0+b)
	y(10) =  (3.0d0*dtht/32.0d0)*(1.0d0-(a**2))*
     &((1.0d0-b)**2)*(1.0d0+b)
	y(11) = -(3.0d0*dtht/32.0d0)*(1.0d0-(a**2))*
     &((1.0d0+b)**2)*(1.0d0-b)
	y(12) =  (3.0d0*dtht/32.0d0)*(1.0d0-(a**2))*
     &((1.0d0+b)**2)*(1.0d0-b)
	y(13) = -(drad*dtht/64.0d0)*(1.0d0-a)*(1.0d0+(3.0d0*a))*
     &((1.0d0-b)**2)*(1.0d0+b)
	y(14) = -(drad*dtht/64.0d0)*(1.0d0+a)*(1.0d0-(3.0d0*a))*
     &((1.0d0-b)**2)*(1.0d0+b)
	y(15) =  (drad*dtht/64.0d0)*(1.0d0+a)*(1.0d0-(3.0d0*a))*
     &((1.0d0+b)**2)*(1.0d0-b)
	y(16) =  (drad*dtht/64.0d0)*(1.0d0-a)*(1.0d0+(3.0d0*a))*
     &((1.0d0+b)**2)*(1.0d0-b)
      
      RETURN
	END

c======================================================================
	SUBROUTINE DNt_annular(a,b,drad,dtht,y)
c======================================================================

c----------------------------------------------------------------------
c	Derivative of Shape function w.r.t angle
c	----------------------------------------
c	a = zeta
c	b = eta
c	drad = Delta radius
c	dtht = Delta theta
c----------------------------------------------------------------------
      
      IMPLICIT NONE
      
	DOUBLE PRECISION a,b,drad,dtht
	DOUBLE PRECISION y(16)

	y(1)  = -(3.0d0/16.0d0)*(1.0d0-(b**2))*
     &((1.0d0-a)**2)*(2.0d0+a)
	y(2)  = -(3.0d0/16.0d0)*(1.0d0-(b**2))*
     &((1.0d0+a)**2)*(2.0d0-a)
	y(3)  =  (3.0d0/16.0d0)*(1.0d0-(b**2))*
     &((1.0d0+a)**2)*(2.0d0-a)
	y(4)  =  (3.0d0/16.0d0)*(1.0d0-(b**2))*
     &((1.0d0-a)**2)*(2.0d0+a)
	y(5)  = -(3.0d0*drad/32.0d0)*(1.0d0-(b**2))*
     &((1.0d0-a)**2)*(1.0d0+a)
	y(6)  =  (3.0d0*drad/32.0d0)*(1.0d0-(b**2))*
     &((1.0d0+a)**2)*(1.0d0-a)
	y(7)  = -(3.0d0*drad/32.0d0)*(1.0d0-(b**2))*
     &((1.0d0+a)**2)*(1.0d0-a)
	y(8)  =  (3.0d0*drad/32.0d0)*(1.0d0-(b**2))*
     &((1.0d0-a)**2)*(1.0d0+a)
	y(9)  = -(dtht/32.0d0)*(1.0d0-b)*(1.0d0+(3.0d0*b))*
     &((1.0d0-a)**2)*(2.0d0+a)
	y(10) = -(dtht/32.0d0)*(1.0d0-b)*(1.0d0+(3.0d0*b))*
     &((1.0d0+a)**2)*(2.0d0-a)
	y(11) = -(dtht/32.0d0)*(1.0d0+b)*(1.0d0-(3.0d0*b))*
     &((1.0d0+a)**2)*(2.0d0-a)
	y(12) = -(dtht/32.0d0)*(1.0d0+b)*(1.0d0-(3.0d0*b))*
     &((1.0d0-a)**2)*(2.0d0+a)
	y(13) = -(drad*dtht/64.0d0)*(1.0d0-b)*(1.0d0+(3.0d0*b))*
     &((1.0d0-a)**2)*(1.0d0+a)
	y(14) =  (drad*dtht/64.0d0)*(1.0d0-b)*(1.0d0+(3.0d0*b))*
     &((1.0d0+a)**2)*(1.0d0-a)
	y(15) =  (drad*dtht/64.0d0)*(1.0d0+b)*(1.0d0-(3.0d0*b))*
     &((1.0d0+a)**2)*(1.0d0-a)
	y(16) = -(drad*dtht/64.0d0)*(1.0d0+b)*(1.0d0-(3.0d0*b))*
     &((1.0d0-a)**2)*(1.0d0+a)
      
      RETURN
	END

c======================================================================
	SUBROUTINE DNrr_annular(a,b,drad,dtht,y)
c======================================================================

c----------------------------------------------------------------------
c	Second Derivative of Shape function w.r.t Radius
c----------------------------------------------------------------------
c	a = zeta
c	b = eta
c	drad = Delta radius
c	dtht = Delta theta
c---------------------------------------------------------------------
      
      IMPLICIT NONE

	DOUBLE PRECISION a,b,drad,dtht
	DOUBLE PRECISION y(16)
	
	y(1)  =  (6.0d0/16.0d0)*a*((1.0d0-b)**2)*(2.0d0+b)
	y(2)  = -(6.0d0/16.0d0)*a*((1.0d0-b)**2)*(2.0d0+b)
	y(3)  = -(6.0d0/16.0d0)*a*((1.0d0+b)**2)*(2.0d0-b)
	y(4)  =  (6.0d0/16.0d0)*a*((1.0d0+b)**2)*(2.0d0-b)
	y(5)  = -(drad/16.0d0)*(1.0d0-(3.0d0*a))*((1.0d0-b)**2)*(2.0d0+b) 
	y(6)  =  (drad/16.0d0)*(1.0d0+(3.0d0*a))*((1.0d0-b)**2)*(2.0d0+b)
	y(7)  =  (drad/16.0d0)*(1.0d0+(3.0d0*a))*((1.0d0+b)**2)*(2.0d0-b) 
	y(8)  = -(drad/16.0d0)*(1.0d0-(3.0d0*a))*((1.0d0+b)**2)*(2.0d0-b) 
	y(9)  =  (3.0d0*dtht/16.0d0)*(a)*((1.0d0-b)**2)*(1.0d0+b)
	y(10) = -(3.0d0*dtht/16.0d0)*(a)*((1.0d0-b)**2)*(1.0d0+b)
	y(11) =  (3.0d0*dtht/16.0d0)*(a)*((1.0d0+b)**2)*(1.0d0-b)
	y(12) = -(3.0d0*dtht/16.0d0)*(a)*((1.0d0+b)**2)*(1.0d0-b)
	y(13) = -(drad*dtht/32.0d0)*(1.0d0-(3.0d0*a))*((1.0d0-b)**2)*
     &		 (1.0d0+b)
	y(14) =  (drad*dtht/32.0d0)*(1.0d0+(3.0d0*a))*((1.0d0-b)**2)*
     &		 (1.0d0+b)
	y(15) = -(drad*dtht/32.0d0)*(1.0d0+(3.0d0*a))*((1.0d0+b)**2)*
     &		 (1.0d0-b)
	y(16) =  (drad*dtht/32.0d0)*(1.0d0-(3.0d0*a))*((1.0d0+b)**2)*
     &		 (1.0d0-b)

      RETURN
	END

c======================================================================
	SUBROUTINE DNtt_annular(a,b,drad,dtht,y)
c======================================================================

c----------------------------------------------------------------------
c	Second Derivative of Shape function w.r.t Angle
c	-----------------------------------------------
c	a = zeta
c	b = eta
c	drad = Delta radius
c	dtht = Delta theta
c----------------------------------------------------------------------
      
      IMPLICIT NONE

	DOUBLE PRECISION a,b,drad,dtht
	DOUBLE PRECISION y(16)
	
	y(1)  =  (6.0d0/16.0d0)*b*((1.0d0-a)**2)*(2.0d0+a)
	y(2)  =  (6.0d0/16.0d0)*b*((1.0d0+a)**2)*(2.0d0-a)
	y(3)  = -(6.0d0/16.0d0)*b*((1.0d0+a)**2)*(2.0d0-a)
	y(4)  = -(6.0d0/16.0d0)*b*((1.0d0-a)**2)*(2.0d0+a)
	y(5)  =  (3.0d0*drad/16.0d0)*(b)*((1.0d0-a)**2)*(1.0d0+a)
	y(6)  = -(3.0d0*drad/16.0d0)*(b)*((1.0d0+a)**2)*(1.0d0-a)
	y(7)  =  (3.0d0*drad/16.0d0)*(b)*((1.0d0+a)**2)*(1.0d0-a)
	y(8)  = -(3.0d0*drad/16.0d0)*(b)*((1.0d0-a)**2)*(1.0d0+a)
	y(9)  = -(dtht/16.0d0)*(1.0d0-(3.0d0*b))*((1.0d0-a)**2)*(2.0d0+a)
	y(10) = -(dtht/16.0d0)*(1.0d0-(3.0d0*b))*((1.0d0+a)**2)*(2.0d0-a)
	y(11) =  (dtht/16.0d0)*(1.0d0+(3.0d0*b))*((1.0d0+a)**2)*(2.0d0-a)
	y(12) =  (dtht/16.0d0)*(1.0d0+(3.0d0*b))*((1.0d0-a)**2)*(2.0d0+a)
	y(13) = -(drad*dtht/32.0d0)*(1.0d0-(3.0d0*b))*((1.0d0-a)**2)*
     &		 (1.0d0+a)
	y(14) =  (drad*dtht/32.0d0)*(1.0d0-(3.0d0*b))*((1.0d0+a)**2)*
     &		 (1.0d0-a)
	y(15) = -(drad*dtht/32.0d0)*(1.0d0+(3.0d0*b))*((1.0d0+a)**2)*
     &		 (1.0d0-a)
	y(16) =  (drad*dtht/32.0d0)*(1.0d0+(3.0d0*b))*((1.0d0-a)**2)*
     &		 (1.0d0+a)

	RETURN
	END

c======================================================================
	SUBROUTINE DNrt_annular(a,b,drad,dtht,y)
c======================================================================

c----------------------------------------------------------------------
c	Second Derivative of Shape function w.r.t Radius and Angle
c	----------------------------------------------------------
c	a = zeta
c	b = eta
c	drad = Delta radius
c	dtht = Delta theta
c----------------------------------------------------------------------

      IMPLICIT NONE

	DOUBLE PRECISION a,b,drad,dtht
	DOUBLE PRECISION y(16)

	y(1)  =  (9.0d0/16.0d0)*(1.0d0-(a**2))*
     &(1.0d0-(b**2))
	y(2)  = -(9.0d0/16.0d0)*(1.0d0-(a**2))*
     &(1.0d0-(b**2))
	y(3)  =  (9.0d0/16.0d0)*(1.0d0-(a**2))*
     &(1.0d0-(b**2))
	y(4)  = -(9.0d0/16.0d0)*(1.0d0-(a**2))*
     &(1.0d0-(b**2))
	y(5)  =  (3.0d0*drad/32.0d0)*(1.0d0-(b)**2)*
     &(1.0d0-a)*(1.0d0+(3.0d0*a))
	y(6)  =  (3.0d0*drad/32.0d0)*(1.0d0-(b)**2)*
     &(1.0d0+a)*(1.0d0-(3.0d0*a))
	y(7)  = -(3.0d0*drad/32.0d0)*(1.0d0-(b)**2)*
     &(1.0d0+a)*(1.0d0-(3.0d0*a))
	y(8)  = -(3.0d0*drad/32.0d0)*(1.0d0-(b)**2)*
     &(1.0d0-a)*(1.0d0+(3.0d0*a))
	y(9)  =  (3.0d0*dtht/32.0d0)*(1.0d0-(a)**2)*
     &(1.0d0-b)*(1.0d0+(3.0d0*b))
	y(10) = -(3.0d0*dtht/32.0d0)*(1.0d0-(a)**2)*
     &(1.0d0-b)*(1.0d0+(3.0d0*b))
	y(11) = -(3.0d0*dtht/32.0d0)*(1.0d0-(a)**2)*
     &(1.0d0+b)*(1.0d0-(3.0d0*b))
	y(12) =  (3.0d0*dtht/32.0d0)*(1.0d0-(a)**2)*
     &(1.0d0+b)*(1.0d0-(3.0d0*b))
	y(13) =  (drad*dtht/64.0d0)*(1.0d0-a)*(1.0d0-b)*
     &(1.0d0+(3.0d0*(a)))*(1.0d0+(3.0d0*(b)))
	y(14) =  (drad*dtht/64.0d0)*(1.0d0+a)*(1.0d0-b)*
     &(1.0d0-(3.0d0*(a)))*(1.0d0+(3.0d0*(b)))
	y(15) =  (drad*dtht/64.0d0)*(1.0d0+a)*(1.0d0+b)*
     &(1.0d0-(3.0d0*(a)))*(1.0d0-(3.0d0*(b)))
	y(16) =  (drad*dtht/64.0d0)*(1.0d0-a)*(1.0d0+b)*
     &(1.0d0+(3.0d0*(a)))*(1.0d0-(3.0d0*(b)))

      RETURN
	END

c======================================================================