	  
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 


      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION, DIMENSION(NDIM/2) :: R, THETA, ROOT, THETADOT, X, Y, T, SECRETE, SENSE
  
      DOUBLE PRECISION, SAVE :: PI = 4d0*ATAN(1d0)
      DOUBLE PRECISION MU, EPS, T0, X0
      INTEGER PDE, J, LEFT, RIGHT, N
	

        
        
        
        N = NDIM/2
		X = U(1:N)
		Y = U(N+1:NDIM)

		
		
		MU=PAR(1)
		EPS=PAR(2)


		
		
		DO 	J=1,N
			IF (N.GT.1) THEN
			X0 = DBLE(J-1)/DBLE(N-1)
			T(J) = (1D0+2d0*PAR(3)*ABS(X0-0.5d0))
			ELSE
			T(1) = 1D0
			 ENDIF
			
		END DO
		
		R=SQRT(X**2+Y**2)
		THETA= ATAN2(Y, X)
		
		ROOT = (MU-R**2)

		SENSE = (1d0-SIN(THETA))/2d0
		
		
		SECRETE(1:N-1) =  (1d0+COS(THETA(2:N)))/2d0
		SECRETE(N) = 0d0 
		SECRETE(2:N) = SECRETE(2:N) + (1d0+COS(THETA(1:N-1)))/2d0
		SECRETE(2:N-1) = SECRETE(2:N-1)/2d0

!~ 		PRINT *, SENSE
		
		THETADOT = 1d0/T  - EPS*SECRETE*SENSE


	  F(1:N) = ROOT*X - Y*THETADOT
	  F(N+1:NDIM) = ROOT*Y + X*THETADOT
	  

	  
      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 
! Define the starting stationary solution on the spatial mesh


      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      DOUBLE PRECISION A, B, K, D, E, BETA, P, L

! Set the parameter values
       PAR(1)=1d0   !START THIS NEGATIVE, TO OBSERVE THE (DEGNERATE) HOPF
       PAR(2)=1d0 !EPS
	   PAR(3)=0d0 !GAMMA

		U(:)=1d-200

      END SUBROUTINE STPNT
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!                Problem-independent subroutines
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
	 
	  
	  
	  

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      
	  SUBROUTINE PHASEDIFF(PI,PJ,DP)
	  DOUBLE PRECISION, INTENT(IN) :: PI, PJ
	  DOUBLE PRECISION, INTENT(INOUT) :: DP
	  DOUBLE PRECISION :: PJP

	  IF (PI.LT.PJ) THEN
		PJP = PJ-1D0
		ELSE
		PJP = PJ
	  END IF
	  
	  DP = PI-PJP
	  
	  END SUBROUTINE
      

      SUBROUTINE PVLS(NDIM,U,PAR)
!     ---------- ----
      use, intrinsic :: ieee_arithmetic
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      DOUBLE PRECISION, EXTERNAL :: GETP
      DOUBLE PRECISION PHI1, PHI2, PHI3, PHIR, PHIL, CHI, PSIL, PSIR, PSI, CHI_BARR, CHI_BARL, CHI_BARM, CHI_PREV, CHI_PREV_PREV
      DOUBLE PRECISION CROSS, CROSS_POS, ABS_CROSS_POS
      LOGICAL LASTCROSS
      LOGICAL, SAVE :: ifrst = .TRUE.
      INTEGER :: N, I, IHAT, MID

      DOUBLE PRECISION, SAVE :: PI = 4d0*ATAN(1d0)


      

      IHAT=0
      PAR(14)=0d0
      PAR(13)=0d0
      PAR(18)=0d0
      PAR(19)=0d0
      PAR(20)=0d0
      

	 PSI=0d0
	 
	 N=NINT(NDIM/2D0)
	 
	 IF (N.GE.3) THEN 
	 CHI_PREV = IEEE_VALUE(CHI_PREV, IEEE_QUIET_NAN)
	 CROSS=0d0
	 CROSS_POS=0d0
	 MID = CEILING(DBLE(N)/2d0)
	 PHI1=GETP('MXT',1,U)
	 IF (.NOT.ISNAN(PHI1)) THEN !for this to work properly, AUTO source must be modified to output NaN for GETP('MXT',1,U) when IPS<=1
	 PSIR=GETP('MXT',3,U)
	 PSIL=GETP('MXT',1,U)
	 CALL PHASEDIFF(PSIR,PSIL, PSI)
	 PSI =  COS(PI*PSI)
	 PAR(17) = NDIM - GETP('STA',0,U)
	 DO  I=1,N-2
		 PHI1=GETP('MXT',I,U)
		 PHI2=GETP('MXT',I+1,U)
		 PHI3=GETP('MXT',I+2,U)
		 
		 
		 CALL PHASEDIFF(PHI2,PHI1, PHIL)
		 CALL PHASEDIFF(PHI2,PHI3, PHIR)
		 
		 CHI = SIN(PI*(PHIR-PHIL))/PI
		 
		 
		 PAR(13)=PAR(13)+(ABS(SIN(PI*PHIR)*SIN(PI*PHIL)))**3
		 PAR(14)=PAR(14)+CHI
		 IF ((I+1).GT.MID) THEN
		 PAR(19)=PAR(19) + CHI/DBLE(MID-2)

		 IHAT=IHAT+1
		 ELSEIF ((I+1).LT.MID) THEN
		 PAR(18)=PAR(18) + CHI/DBLE(MID-2)

		 IHAT=IHAT+1
		 ELSE
		 PAR(21)=CHI
		 END IF
		 
		 IF (ieee_is_finite(CHI_PREV)) THEN
			 IF ((ABS(SIGN(1d0,CHI_PREV)-SIGN(1d0,CHI)).GE.1d0).AND.((ABS(CHI).GT.1d-3).OR.(ABS(CHI_PREV).GT.1d-3))) THEN
			    IF (CROSS.EQ.0d0) THEN
					CROSS_POS=0d0
			    END IF
				CROSS = CROSS + 1d0
				CROSS_POS = CROSS_POS +  dble(I-MID)-CHI_PREV/(CHI-CHI_PREV)
				ABS_CROSS_POS = ABS_CROSS_POS + ABS(dble(I-MID)-CHI_PREV/(CHI-CHI_PREV))
				!PRINT *, I+1, SIGN(1d0,CHI_PREV), PAR(21)
			 END IF 
		 
		 END IF
		 CHI_PREV = CHI
     END DO
     PAR(13)=PAR(13)/DBLE(N-2) 
     PAR(14)=PAR(14)/DBLE(N-2)
     PAR(16)=PAR(19)-PAR(18)
     PAR(20)=PAR(16)/4d0+PAR(14)
     PAR(22)=CROSS
     IF (CROSS.GT.0) THEN
		PAR(23)=(CROSS_POS/CROSS)/DBLE(MID-1)
		PAR(24)=(ABS_CROSS_POS/CROSS)/DBLE(MID-1)
	 ELSE
		PAR(23)=0d0
		PAR(24)=0d0
     END IF
     END IF
		
     PAR(15) =(GETP('MAX',2,U)+(GETP('MAX',1,U)-GETP('MAX',3,U))/15d0 + PSI/2d0)/8d0
     END IF

      END SUBROUTINE PVLS



