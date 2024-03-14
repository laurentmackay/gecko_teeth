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
      MODULE sensing
		DOUBLE PRECISION, DIMENSION(8, 2) :: X0, NORMAL
      END MODULE sensing
      
      MODULE RD
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
        
      END MODULE RD

	  
      SUBROUTINE FUNC(NDIM,SOL,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      USE RD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: SOL(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION, DIMENSION(NCOMP) :: T, SENSE, GAMMA, IAVE, U,V
      DOUBLE PRECISION I(NX), DIVFLUX(NX), N

      DOUBLE PRECISION A, B, K, GRAD, E, KAPPA2, P, F0, SENSE_PHASE, EPS,  X
      DOUBLE PRECISION, DIMENSION(NCOMP, 2) :: UV, DIFF
      DOUBLE PRECISION, DIMENSION(2) :: X0, NORMAL
      DOUBLE PRECISION, SAVE :: PI = 4d0*ATAN(1d0), T0=2d-1
      INTEGER PDE, J, LEFT, RIGHT

! Problem-independent initialization :
        
        IF (.NOT.ALLOCATED(IJUMP)) THEN
            CALL INITGRID(NDIM, PAR)
        END IF
        
        PDE = NODE*NCOMP+1
        

		U = SOL(1:NCOMP)
		V = SOL(NCOMP+1:2*NCOMP)
		UV = RESHAPE(SOL(1:2*NCOMP), (/ NCOMP, 2/))
		I = SOL(NODE*NCOMP+1:NDIM)
		
		A=PAR(1)
		B=PAR(2)
		E=PAR(3)
		K=PAR(4)
        KAPPA2 = PAR(5)
		P=PAR(6)
		GRAD = PAR(7)
		SENSE_PHASE = PAR(8)
		EPS = PAR(9)
		F0 = PAR(10)
		
		
		N=DBLE(NCOMP)

		
		CALL GET_SENSING_VECTORS(X0, NORMAL, SENSE_PHASE, F0)
		
		DO 	J=1,NCOMP
			IF (NCOMP.GT.1) THEN
			X = DBLE(J-1)/(N-1d0)
			T(J) = T0*(1D0+PAR(7)*(X-0.5d0))
			ELSE
			T(1) = T0
			 ENDIF
			SENSE(J) = 1d0 - (1d0 + TANH(SUM( (UV(J,:)-X0)*NORMAL )/1D-2 ))/2D0
		END DO
		
		
	
		GAMMA = V
        IAVE = (I(IJUMP)+I(IJUMP+1))/2d0


		F(1:NCOMP) = (U*(1d0-U)*(U-A) - V + F0)/T - EPS*SENSE*IAVE
		F(NCOMP+1:2*NCOMP) = E*(K*U-V-B)/T
		
		
		CALL DIFFUSE(I, -p*GAMMA, KAPPA2*p, DIVFLUX)
		
		
		
		F(PDE:NDIM) = -DIVFLUX - p*I
		
       IF (IJAC .NE. 0) THEN
        call exit(1)
	   END IF
		
      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 
! Define the starting stationary solution on the spatial mesh

      USE RD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      DOUBLE PRECISION A, B, K, GRAD, E, KAPPA2, P, F0, SENSE_PHASE, EPS
      
      A=0.169
      B=0.16d0
      E=1d-2
      K=0.6d0
      KAPPA2 =1d0
      EPS = 1d0
      GRAD=0d0
      P=1d-1
      F0=0.06d0
      SENSE_PHASE=6d0
      
! Set the parameter values
		PAR(1) = A
		PAR(2) = B
		PAR(3) = E
		PAR(4) = K
        PAR(5) = KAPPA2
		PAR(6) = P
	    PAR(7) = GRAD
		PAR(8) = SENSE_PHASE
		PAR(9) = EPS
		PAR(10) = F0

	   U(:)=1d-200

      END SUBROUTINE STPNT
      
      SUBROUTINE GET_SENSING_VECTORS(ORIGIN,DIR,PHASE, F0)
      USE sensing
      DOUBLE PRECISION, INTENT(INOUT) :: ORIGIN(2), DIR(2)
      DOUBLE PRECISION, INTENT(IN) :: PHASE, F0
      DOUBLE PRECISION :: ALPHA    
      I=NINT(MOD(FLOOR(PHASE)-1D0, 8D0)+1D0)
      J=NINT(MOD(CEILING(PHASE)-1D0,8D0)+1D0)
      ALPHA = MOD(PHASE, 1D0)
      ORIGIN = (1d0-ALPHA)*X0(I,:) + ALPHA*X0(J,:)
      DIR = (1d0-ALPHA)*NORMAL(I,:) + ALPHA*NORMAL(J,:)
      
      ORIGIN(2)=ORIGIN(2)+F0
      END SUBROUTINE
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
!                Problem-independent subroutines
!---------------------------------------------------------------------- 
!---------------------------------------------------------------------- 
	  SUBROUTINE DIFFUSE(U, FJUMP, D, DIVFLUX)
	  USE RD
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
      use RD
      use sensing
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
      
		X0(1,:)=(/0.75000d0, 0.00000d0/)
		NORMAL(1,:)=(/-0.55470d0, -0.83205d0/)
		 
		X0(2,:)=(/0.80000d0, 0.08000d0/)
		NORMAL(2,:)=(/-0.12403d0, -0.99228d0/)
		 
		X0(3,:)=(/0.00000d0, 0.11000d0/)
		NORMAL(3,:)=(/0.01000d0, -0.99995d0/)
		 
		X0(4,:)=(/0.00000d0, 0.08000d0/)
		NORMAL(4,:)=(/0.16440d0, -0.98639d0/)
		 
		X0(5,:)=(/0.00000d0, 0.04000d0/)
		NORMAL(5,:)=(/0.44721d0, 0.89443d0/)
		 
		X0(6,:)=(/0.00000d0, 0.02000d0/)
		NORMAL(6,:)=(/0.11043d0, 0.99388d0/)
		 
		X0(7,:)=(/0.00000d0, 0.00500d0/)
		NORMAL(7,:)=(/-0.03331d0, 0.99944d0/)
		 
		X0(8,:)=(/0.80000d0, 0.04000d0/)
		NORMAL(8,:)=(/-0.24254d0, 0.97014d0/)          
            
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
      use RD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      DOUBLE PRECISION, EXTERNAL :: GETP
      DOUBLE PRECISION PHI1, PHI2, PHI3, PHIR, PHIL, CHI, CHIR(NCOMP)
      LOGICAL, SAVE :: ifrst = .TRUE.
      INTEGER :: N, I, IR
      DOUBLE PRECISION, SAVE :: PI = 4d0*ATAN(1d0)

! Problem-independent initialization :
      IF(ifrst.AND.(.NOT.ALLOCATED(IJUMP)))THEN
         PRINT *, "first time"
         CALL INITGRID(NDIM, PAR)
         ifrst=.FALSE.
      ENDIF
      
      CHIR=2d0
      IR=1
      PAR(14)=0d0
      PAR(13)=0d0

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
		 IF (I.GT.CEILING(DBLE(N)/2d0)) THEN
		 CHIR(IR)=CHI
		 IR=IR+1
		 END IF
     END DO
     PAR(13)=PAR(13)/DBLE(N) + (GETP('MAX',2,U)-1d0)/4d0
     PAR(14)=PAR(14)/DBLE(N)
     
     PAR(15)=SQRT(PAR(5))
     PAR(16)=MINVAL(CHIR(1:IR-1))
     PAR(17) = NDIM - GETP('STA',0,U)
     END IF

      END SUBROUTINE PVLS
