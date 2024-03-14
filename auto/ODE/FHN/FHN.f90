!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!  FHN : array of diffusionally coupled fitzhugh-nagumo oscillators with PBC
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE FUNC(NDIM,UV,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: UV(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION A, B, K, D, E
      DOUBLE PRECISION U(NDIM/2), V(NDIM/2), D2V(NDIM/2)
      INTEGER N, I
       N=INT(NDIM/2)
	   U=UV(1:N)
	   V=UV(N+1:NDIM)
	   

		IF (N > 1) THEN

	   D2V(1) = (2.0*V(1)-V(2)-V(N))/2.0
		   D2V(1) = (V(1)-V(2))
		   DO I = 2, N-1
		   D2V(I) = (2.0*V(I)-V(I+1)-V(I-1))
		   END DO
		   !D2V(N) = (2.0*V(N)-V(1)-V(N-1))/2.0
		   D2V(N) = (V(N)-V(N-1))
       ELSE
			D2V(1) = 0D0
       END If
       
       B=PAR(1)
       D=PAR(2)
       A=PAR(3)
       K=PAR(4)
       E=PAR(5)



       F(1:N)=U*(1d0-U)*(U-A)-V
       F(N+1:NDIM)=E*(K*U-V-B)-(D*D2V)
       
       IF (IJAC .NE. 0) THEN
       
		   DO I = 1, N
			   DFDU(I,I)=(1d0-U(I))*(U(I)-A)  +  U(I)*(-(U(I)-A)+(1d0-U(I)))
			   
			   DFDU(I,I+N) = -1d0
			   
			   DFDU(I+N,I) = E*K
			   !!!! THIS NEEDS FIXING TO WORK WITH ZERO FLUX
			   IF (N > 1) THEN
				   IF (I .EQ. 1) THEN
					   DFDU(I+N,I+N) = -1d0*D-E
					   DFDU(N+1, N+2) = 1d0*D
				   ELSEIF (I .EQ. N) THEN
					   DFDU(NDIM, NDIM-1) = 1d0*D
					   DFDU(NDIM,NDIM) = -1d0*D-E
				   ELSE 
					   DFDU(I+N,I+N) = -2d0*D-E
					   DFDU(N+I, N+I+1) = 1d0*D
					   DFDU(N+I, N+I-1) = 1d0*D
				   END IF
			   END IF
		   END DO
		   
		   IF(IJAC.EQ.1)RETURN
		   !print *, IJAC
		   DO I = 1, N
		   
		   DFDP(I,3) = -U(I)*(1d0-U(I)) !A
		   
		   DFDP(I+N,1) = -E !B
		   IF (N > 1) THEN
				DFDP(I+N,2) = -(D2V(I)/2.0) !D
		   END IF
		   DFDP(I+N,4) = E*U(I) ! K
		   DFDP(I+N,5) = (K*U(I)-V(I)-B) ! E
		   END DO
	   END IF

      END SUBROUTINE FUNC
!---------------------------------------------------------------------- 

      SUBROUTINE STPNT(NDIM,U,PAR,T) 
!     ---------- ----- 

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

       PAR(1)=0.0   !b, simplest approach is to set to zero and continue from trivial soln
       PAR(2)=0.0003 !D
       PAR(3)=0.139  !A
       PAR(4)=0.6    !k
       PAR(5)=1d-2   !eps

       U(:)=0d0

      END SUBROUTINE STPNT
!---------------------------------------------------------------------- 

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT

           SUBROUTINE PVLS(NDIM,U,PAR)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)

      DOUBLE PRECISION, EXTERNAL :: GETP,GETU2
      INTEGER NDX,NCOL,NTST, I, N
      DOUBLE PRECISION  T1, PHI
      
      ! The first argument of GETP may be one of the following:
!        'NRM' (L2-norm),     'MAX' (maximum),
!        'INT' (integral),    'BV0 (left boundary value),
!        'MIN' (minimum),     'BV1' (right boundary value).
!        'MNT' (t value for minimum)
!        'MXT' (t value for maximum)
!        'NDIM', 'NDX' (effective (active) number of dimensions)
!        'NTST' (NTST from constant file)
!        'NCOL' (NCOL from constant file)
!        'NBC'  (active NBC)
!        'NINT' (active NINT)
!        'DTM'  (delta t for all t values, I=1...NTST)
!        'WINT' (integration weights used for interpolation, I=0...NCOL)
       T1 =  GETP('MXT',1,U)
       N = INT(NDIM/2)
       PHI=0d0
       IF (N>1) THEN
       DO I = 2, N
       PHI = PHI + ABS(GETP('MXT',I,U)-T1)
	   END DO
		PAR(6)=PHI/(N-1)
		END IF
      END SUBROUTINE PVLS
