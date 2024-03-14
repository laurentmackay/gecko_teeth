!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!     brf  :  A parabolic PDE (the Brusselator)
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!     (Discretized in space by fourth order finite differences)
!---------------------------------------------------------------------- 
!----------------------------------------------------------------------
! NOTE: The value of the constant NE is defined below in a module.
!
!      NE  :  the dimension of the PDE system
!
!      NX  :  the number of space intervals for the discretization is
!             derived from the AUTO-constant NDIM:
!             NX = NDIM/NE + 1
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      MODULE wtv
        SAVE
        INTEGER, PARAMETER :: NCOMP=5
        INTEGER, PARAMETER :: NODE=2
        INTEGER, PARAMETER :: NPDE=1
		INTEGER :: NX=100
		DOUBLE PRECISION :: DX
		DOUBLE PRECISION :: DX2
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DI
        INTEGER, ALLOCATABLE, DIMENSION(:) :: IJUMP ! WE HAVE JUMPS AT THE INTERFACE BETWEEN CELLS I AND I+1
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: STEN
        
      END MODULE wtv

	  
      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      USE wtv
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION, DIMENSION(NCOMP) :: R, THETA, ROOT, THETADOT, X, Y, T, SECRETE, SENSE, COMPFLUX(NCOMP), OMEGA0
      DOUBLE PRECISION INH(NX), DIVFLUX(NX)

      DOUBLE PRECISION MU, T0, E, KAPPA2, P, X0 , N, OMEGA
      DOUBLE PRECISION, SAVE :: PI = 4d0*ATAN(1d0)
      INTEGER PDE, J, LEFT, RIGHT

! Problem-independent initialization :
        
        IF (.NOT.ALLOCATED(IJUMP)) THEN
            CALL INITGRID(NDIM, PAR)
        END IF
        
        PDE = NODE*NCOMP+1
        
		X = U(1:NCOMP)
		Y = U(NCOMP+1:2*NCOMP)
		INH = U(PDE:NDIM)
		
		MU=PAR(1)
		N=DBLE(NCOMP)
!~ 	    L=PAR(7)
		OMEGA=PAR(6)
		T0 = PAR(3)
		
		P=2d0*PI/(OMEGA*T0)
		
		F = PAR(2)
		
		DO 	J=1,NCOMP
			IF (NCOMP.GT.1) THEN
			X0 = DBLE(J-1)/(N-1d0)
			T(J) = T0*(1D0+2D0*PAR(4)*ABS(X0-0.5d0))
			ELSE
			T(1) = T0
			 ENDIF
			
		END DO
		
		KAPPA2 = PAR(5)
		
!~ 		D=P*(par(5)*L/DBLE(NCOMP))**2
		
		R=SQRT(X**2+Y**2)
		THETA= ATAN2(Y, X)
		
		ROOT = (MU-R**2)

		SENSE = (1d0-SIN(THETA))/2d0
		
		
		SECRETE =  (1d0+COS(THETA))/2d0

		
		OMEGA0 = 2d0*PI/T0

		THETADOT = OMEGA0*(T0/T  - 0.5d0*F*SENSE*(INH(IJUMP)+INH(IJUMP+1)))

!~ 		THETADOT = (2d0*PI/T0)*(T0/T - F*L*p*(INH(IJUMP)+INH(IJUMP+1))*SENSE/N)
		F(1:NCOMP) = ROOT*X - Y*THETADOT
		F(NCOMP+1:2*NCOMP) = ROOT*Y + X*THETADOT
		
!~ 		print *, p*SECRETE
		
		CALL DIFFUSE(INH, -p*2d0*SECRETE, KAPPA2*p, DIVFLUX)
		
		
		
		F(PDE:NDIM) = -DIVFLUX - p*INH
		
!~ 		print *, "U(:)"
!~ 		print *, U
		
!~ 		print *, "F(:)"
!~ 		PRINT *, F
		
!~ 		print *, D
!~ 		print *, p
		
!~ 		call exit(1)
		
!~        IF (IJAC .NE. 0) THEN
!~         call exit(1) !THIS NEEDS TO BE UPDATED TO MATCH THE PHASE OSCILLATOR DYANMICS

!~ 		   DFDU(PDE, PDE) = -KAPPA/DX2 - p
!~ 		   DFDU(PDE, PDE+1) = KAPPA/DX2
!~ 		   DO J = 1, NX-2
!~ 			   DFDU(PDE+J, PDE+J-1) = D/DX2
!~ 			   DFDU(PDE+J, PDE+J) = -2*D/DX2 - p
!~ 			   DFDU(PDE+J, PDE+J+1) = D/DX2
!~ 		   END DO
!~ 		   DFDU(PDE+NX-1, PDE+NX-1) = -D/DX2 - p
!~ 		   DFDU(PDE+NX-1, PDE+NX-2) = D/DX2
		 
!~ 		   DO J = 1, NCOMP
		   

		   
!~ 			DFDU(J,J)=1d0 - 3d0*X(J)**2 - (E*Cos(THETA(J))*(Inh(IJUMP(J)) + Inh(1 + IJUMP(J)))*X(J)*Y(J))/(4d0*R(J)**2) - Y(J)**2
!~ 			DFDU(J,J + NCOMP)=-OMEGA0(J) - (E*(Inh(IJUMP(J)) + Inh(1 + IJUMP(J)))*(-1d0 + Sin(THETA(J))))/4d0 - 2*X(J)*Y(J) + (E*Cos(THETA(J))*(Inh(IJUMP(J)) + Inh(1 + IJUMP(J)))*Y(J)**2)/(4d0*R(J)**2)
!~ 			DFDU(J,-1 + PDE + IJUMP(J))=-(E*(-1d0 + Sin(THETA(J)))*Y(J))/4d0
!~ 			DFDU(J,PDE + IJUMP(J))=-(E*(-1d0 + Sin(THETA(J)))*Y(J))/4d0
!~ 			DFDU(J + NCOMP,J)=OMEGA0(J) + (E*(Inh(IJUMP(J)) + Inh(1 + IJUMP(J)))*(-1d0 + Sin(THETA(J))))/4d0 + (E*Cos(THETA(J))*(Inh(IJUMP(J)) + Inh(1 + IJUMP(J)))*X(J)**2)/(4d0*R(J)**2) - 2d0*X(J)*Y(J)
!~ 			DFDU(J + NCOMP,J + NCOMP)=1d0 - X(J)**2 - (E*Cos(THETA(J))*(Inh(IJUMP(J)) + Inh(1 + IJUMP(J)))*X(J)*Y(J))/(4d0*R(J)**2) - 3d0*Y(J)**2
!~ 			DFDU(J + NCOMP,-1 + PDE + IJUMP(J))=(E*(-1d0 + Sin(THETA(J)))*X(J))/4d0
!~ 			DFDU(J + NCOMP,PDE + IJUMP(J))=(E*(-1d0 + Sin(THETA(J)))*X(J))/4d0
!~ 			DFDU(-1 + PDE + IJUMP(J),J)=-(Sin(THETA(J))*X(J))/(4d0*DX*R(J)**2)
!~ 			DFDU(-1 + PDE + IJUMP(J),J + NCOMP)=(Sin(THETA(J))*Y(J))/(4d0*DX*R(J)**2)

!~ 			DFDU(PDE + IJUMP(J),J)=-(Sin(THETA(J))*X(J))/(4d0*DX*R(J)**2)
!~ 			DFDU(PDE + IJUMP(J),J + NCOMP)=(Sin(THETA(J))*Y(J))/(4d0*DX*R(J)**2)

!~ 			DFDU(J,J) = 1 - 3*X(J)**2 - (E*Cos(THETA(J))*(Inh(IJUMP(J)) + Inh(1 + IJUMP(J)))*X(J)*Y(J))/(4.*R(J)**2) - Y(J)**2
!~ 			DFDU(J,J + NCOMP) = -OMEGA0(J) - (E*(Inh(IJUMP(J)) + Inh(1 + IJUMP(J)))*(-1 + Sin(THETA(J))))/4. - 2*X(J)*Y(J) + (E*Cos(THETA(J))*(Inh(IJUMP(J)) + Inh(1 + IJUMP(J)))*Y(J)**2)/(4.*R(J)**2)
!~ 			DFDU(J,-1 + PDE + IJUMP(J)) = -(E*(-1 + Sin(THETA(J)))*Y(J))/4.
!~ 			DFDU(J,PDE + IJUMP(J)) = -(E*(-1 + Sin(THETA(J)))*Y(J))/4.
!~ 			DFDU(J + NCOMP,J) = OMEGA0(J) + (E*(Inh(IJUMP(J)) + Inh(1 + IJUMP(J)))*(-1 + Sin(THETA(J))))/4. + (E*Cos(THETA(J))*(Inh(IJUMP(J)) + Inh(1 + IJUMP(J)))*X(J)**2)/(4.*R(J)**2) - 2*X(J)*Y(J)
!~ 			DFDU(J + NCOMP,J + NCOMP) = 1 - X(J)**2 - (E*Cos(THETA(J))*(Inh(IJUMP(J)) + Inh(1 + IJUMP(J)))*X(J)*Y(J))/(4.*R(J)**2) - 3*Y(J)**2
!~ 			DFDU(J + NCOMP,-1 + PDE + IJUMP(J)) = (E*(-1 + Sin(THETA(J)))*X(J))/4.
!~ 			DFDU(J + NCOMP,PDE + IJUMP(J)) = (E*(-1 + Sin(THETA(J)))*X(J))/4.
!~ 			DFDU(-1 + PDE + IJUMP(J),J) = -(Sin(THETA(J))*X(J))/(4.*DX*R(J)**2)
!~ 			DFDU(-1 + PDE + IJUMP(J),J + NCOMP) = (Sin(THETA(J))*Y(J))/(4.*DX*R(J)**2)

!~ 			DFDU(PDE + IJUMP(J),J) = -(Sin(THETA(J))*X(J))/(4.*DX*R(J)**2)
!~ 			DFDU(PDE + IJUMP(J),J + NCOMP) = (Sin(THETA(J))*Y(J))/(4.*DX*R(J)**2)

		   
!~ 		   END DO
		   
		   
!~ 		   IF(IJAC.EQ.1)RETURN
!~ 		   call exit(1)
		   !print *, IJAC
!~ 		   DO I = 1, N
		   
!~ 		   DFDP(I,3) = -U(I)*(1d0-U(I)) !A
		   
!~ 		   DFDP(I+N,1) = -E !B
!~ 		   DFDP(I+N,2) = -(D2V(I)/2.0) !D
!~ 		   DFDP(I+N,4) = E*U(I) ! K
!~ 		   DFDP(I+N,5) = (K*U(I)-V(I)-B) ! E
!~ 		   END DO
!~ 	   END IF
		

		
      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 
! Define the starting stationary solution on the spatial mesh

      USE wtv
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T
      DOUBLE PRECISION, SAVE :: PI = 4d0*ATAN(1d0)
      DOUBLE PRECISION A, B, K, D, E, BETA, P, L

! Set the parameter values
       PAR(1)=1d0   !START THIS NEGATIVE, TO OBSERVE THE (DEGNERATE) HOPF
       PAR(2)=0.1d0 !f
	   PAR(3)=30d0 !T0
	   PAR(4)=0d0 !GRADIENT
	   PAR(5)=2d0 !kappa
	   PAR(6)=2d0*PI !omega


	   U(:)=1d-200

      END SUBROUTINE STPNT
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!                Problem-independent subroutines
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
	  SUBROUTINE DIFFUSE(U, FJUMP, D, DIVFLUX)
	  USE wtv
	  IMPLICIT NONE
	  DOUBLE PRECISION, INTENT(IN) :: U(NX),  FJUMP(NCOMP), D
	  DOUBLE PRECISION, INTENT(OUT) :: DIVFLUX(NX)
	  DOUBLE PRECISION  FLUXR(NX), FLUXL(NX), UINT(NCOMP)
	  INTEGER I

     ! APPROXIMATION FOR U AT THE INTERFACES,
     ! SOLVES D ((U(IJUMP+1) - UJUMP)/(DX/2) - (UJUMP - U(IJUMP))/(DX/2)) == FJUMP
	  UINT =  (U(IJUMP) + U(IJUMP+1))/2d0 - DX*FJUMP/(4d0*D) 
	

	  
	  
	  ! IMPLEMENTATION SPECIFIC CODE: HIGHER-ORDER FINITE DIFFERENCING COULD BE USED (JACOBIAN WILL NEED ADJUSTING)
	  FLUXR(1:NX) = -D*(U(2:NX)-U(1:NX-1))/DX !FIRST ORDER FORWARD DIFFERENCE	
	  FLUXR(NX) = 0d0 !NEUMANN BOUNDARY CONDITIONS
	  
	  FLUXL(1) = 0d0 !NEUMANN BOUNDARY CONDITIONS
	  FLUXL(2:NX) = FLUXR(1:NX-1)
	  
	  FLUXR(IJUMP) = -2d0*D*(UINT-U(IJUMP))/DX ! FLUX TO THE LEFT OF THE INTERFACE
	  FLUXL(IJUMP+1) = -2d0*D*(U(IJUMP+1)-UINT)/DX ! FLUX TO THE RIGHT OF THE INTERFACE
	  
	  

	  DIVFLUX = (FLUXR-FLUXL)/DX

	  
	  END SUBROUTINE DIFFUSE
	  
	  
	  

      SUBROUTINE BCND 
      END SUBROUTINE BCND

      SUBROUTINE ICND 
      END SUBROUTINE ICND

      SUBROUTINE FOPT 
      END SUBROUTINE FOPT
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 

      SUBROUTINE INITGRID(NDIM, PAR)
      use wtv
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: PAR(*)
      INTEGER I, J
      DOUBLE PRECISION TEMP, L, XCOMP, H, SPACE
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: XINT !POSITION OF CELL INTERFACES
      
      TEMP = DBLE(NDIM - (NCOMP*NODE)) / DBLE(NPDE)
      
      IF (MOD(TEMP, 1d0) .NE. 0D0) THEN
		  PRINT *, "error: NDIM cannot be split into an integer number of grid points for the PDEs"
		  call EXIT(1)
      ELSE
		NX = INT(TEMP)
      END IF
      
      ALLOCATE(IJUMP(NCOMP))
      ALLOCATE(XINT(NX-1))
      
      
      L=DBLE(NCOMP)
      DX = L/DBLE(NX)
      DX2 = DX*DX
      H = L/DBLE(NCOMP)
      
      DO I=1,NX-1
		XINT(I) = (I-0.5)*DX 
      END DO
      
      SPACE = DBLE(NX)/DBLE(NCOMP)
      
      DO I=1, NCOMP
	      IJUMP(I) = NINT((I-0.5)*SPACE)
!~ 		  XCOMP =(I-0.5)*H 
!~ 		  DO J=1,NX-2
!~ 			IF ( XCOMP.GE.XINT(J)  ) THEN
!~ 				IJUMP(I)=J
!~ 			END IF
!~ 		  END DO
      END DO
      
      
      PRINT *,"SPACING: ", SPACE
      
      IF (MOD(SPACE,1.0).NE.0d0) THEN
		PRINT *, "WARNING: COMPARTMENTS CANNOT BE EVENLY SPREAD THROUGH THE GRID"
		PRINT *, "SUGGESTED NDIM:"
		PRINT *, NDIM + CEILING(SPACE)*NCOMP-NX
      END IF
      
      PRINT *, NX
      PRINT *, IJUMP

      
      END SUBROUTINE INITGRID
      
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
      use wtv
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      DOUBLE PRECISION, EXTERNAL :: GETP
      DOUBLE PRECISION PHI1, PHI2, PHI3, PHIR, PHIL, CHI
      LOGICAL, SAVE :: ifrst = .TRUE.
      INTEGER :: N, I, IHAT
      DOUBLE PRECISION, SAVE :: PI = 4d0*ATAN(1d0)

! Problem-independent initialization :
      IF(ifrst.AND.(.NOT.ALLOCATED(IJUMP)))THEN
         PRINT *, "first time"
         CALL INITGRID(NDIM, PAR)
         ifrst=.FALSE.
      ENDIF
      

      IHAT=0
      PAR(14)=0d0
      PAR(13)=0d0
      PAR(16)=0d0

	 N=NCOMP-2
	 IF (N.GE.1) THEN 
	 DO  I=1,N
		 PHI1=GETP('MXT',I,U)
		 PHI2=GETP('MXT',I+1,U)
		 PHI3=GETP('MXT',I+2,U)
		 
		 
		 CALL PHASEDIFF(PHI2,PHI1, PHIL)
		 CALL PHASEDIFF(PHI2,PHI3, PHIR)
		 
		 CHI = SIN(PI*(PHIR-PHIL))
		 
		 PAR(13)=PAR(13)+(1d0-ABS(COS(PI*(PHIR+PHIL)/2d0)))
		 PAR(14)=PAR(14)+CHI
		 IF ((I.GT.1).AND.(I.LT.NCOMP)) THEN
			 IF (I.GE.CEILING(DBLE(NCOMP)/2D0)) THEN
				 PAR(16)=PAR(16)+CHI
				 IHAT=IHAT+1		 
			 ELSEIF (I.LE.FLOOR(DBLE(NCOMP)/2D0)) THEN
				 PAR(16)=PAR(16)-CHI
				 IHAT=IHAT+1
			 END IF
		 END IF
		 
     END DO
     par(19)=0d0
     DO  I=1,NDIM
     par(19)=par(19)+GETP('NRM',I,U)
     END DO
     
     PAR(13)=PAR(13)/DBLE(N) 
     PAR(14)=PAR(14)/DBLE(N)
     
     PAR(15)=SQRT(PAR(5))
     PAR(16)=PAR(16)/DBLE(IHAT)
     PAR(17) = NDIM - GETP('STA',0,U)
     PAR(18) = PAR(16) 
     PAR(19) = PAR(19)/(2.5d0*DBLE(NDIM))-1.75d0
     PAR(20) = LOG(PAR(6))
     END IF

      END SUBROUTINE PVLS
