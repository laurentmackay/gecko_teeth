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
        INTEGER, PARAMETER :: NCOMP=3
        INTEGER, PARAMETER :: NODE=2
        INTEGER, PARAMETER :: NPDE=1
		INTEGER :: NX
		DOUBLE PRECISION :: DX
		DOUBLE PRECISION :: DX2
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DI
        INTEGER, ALLOCATABLE, DIMENSION(:) :: IJUMP ! WE HAVE JUMPS AT THE INTERFACE BETWEEN CELLS I AND I+1
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: STEN
        
      END MODULE wtv

	  
      SUBROUTINE FUNC(NDIM,SOL,ICP,PAR,IJAC,F,DFDU,DFDP) 
!     ---------- ---- 

      USE wtv
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: SOL(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)

      DOUBLE PRECISION U(NCOMP), V(NCOMP), COMPFLUX(NCOMP)
      DOUBLE PRECISION I(NX), DIVFLUX(NX)

      DOUBLE PRECISION A, B, K, D, E, BETA, P
      INTEGER PDE, J, LEFT, RIGHT

! Problem-independent initialization :

		U = SOL(1:NCOMP)
		V = SOL(NCOMP+1:2*NCOMP)
		I = SOL(NODE*NCOMP+1:NDIM)
		
		B=PAR(1)
		D=PAR(2)
		A=PAR(3)
		K=PAR(4)
		E=PAR(5)
		BETA=PAR(6)
		P=PAR(7)
		
		
		COMPFLUX = BETA*((I(IJUMP)+I(IJUMP+1))/2d0 -V)

		
		F(1:NCOMP) = U*(1d0-U)*(U-A)-V
		F(NCOMP+1:2*NCOMP) = E*(K*U-V-B) + COMPFLUX
		
		
		
		CALL DIFFUSE(I, COMPFLUX, D, DIVFLUX)
		
		PDE = NODE*NCOMP+1
		
		F(PDE:NDIM) = - DIVFLUX - p*I 
		
		
       IF (IJAC .NE. 0) THEN
       
		   DO J = 1, NCOMP
			   DFDU(J,J)=(1d0-U(J))*(U(J)-A)  +  U(J)*(-(U(J)-A)+(1d0-U(J)))
			   
			   DFDU(J,J+NCOMP) = -1d0
			   
			   DFDU(J+NCOMP,J) = E*K
			   DFDU(J+NCOMP,J+NCOMP) = - E - BETA
			   
			   DFDU(J+NCOMP, PDE+IJUMP(J) ) = BETA/2d0
			   DFDU(J+NCOMP, PDE+IJUMP(J)-1 ) = BETA/2d0
			   
			   DFDU(PDE+IJUMP(J), J+NCOMP) = BETA/(2d0*DX)
			   DFDU(PDE+IJUMP(J)-1, J+NCOMP) = BETA/(2d0*DX)
			   
		   END DO
		   
		   DFDU(PDE, PDE) = -D/DX2 - p
		   DFDU(PDE, PDE+1) = D/DX2
		   DO J = 1, NX-2
			   DFDU(PDE+J, PDE+J-1) = D/DX2
			   DFDU(PDE+J, PDE+J) = -2*D/DX2 - p
			   DFDU(PDE+J, PDE+J+1) = D/DX2
		   END DO
		   DFDU(PDE+NX-1, PDE+NX-1) = -D/DX2 - p
		   DFDU(PDE+NX-1, PDE+NX-2) = D/DX2
		 
		   DO J = 1, NCOMP
		   LEFT = PDE+IJUMP(J)-1
		   RIGHT = PDE+IJUMP(J)
		   DFDU(LEFT, LEFT) =  DFDU(LEFT, LEFT) - BETA/(4d0*DX)
		    DFDU(LEFT, RIGHT) = DFDU(PDE+IJUMP(J)-1, PDE+IJUMP(J)) - BETA/(4d0*DX)
		   
		   
		   DFDU(PDE+IJUMP(J), PDE+IJUMP(J)-1) = DFDU(PDE+IJUMP(J), PDE+IJUMP(J)-1) - BETA/(4d0*DX)
		   DFDU(PDE+IJUMP(J), PDE+IJUMP(J)) = DFDU(PDE+IJUMP(J), PDE+IJUMP(J)) - BETA/(4d0*DX)
		   END DO
		   IF(IJAC.EQ.1)RETURN
		   !print *, IJAC
!~ 		   DO I = 1, N
		   
!~ 		   DFDP(I,3) = -U(I)*(1d0-U(I)) !A
		   
!~ 		   DFDP(I+N,1) = -E !B
!~ 		   DFDP(I+N,2) = -(D2V(I)/2.0) !D
!~ 		   DFDP(I+N,4) = E*U(I) ! K
!~ 		   DFDP(I+N,5) = (K*U(I)-V(I)-B) ! E
!~ 		   END DO
	   END IF
		

		
      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- ----- 
! Define the starting stationary solution on the spatial mesh

      USE wtv
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      DOUBLE PRECISION A, B, K, D, E, BETA, P, L

! Set the parameter values
       PAR(1)=0.0   !b, for simplicity continue from b=0 (also ensure that other parameters are chosen to place the steady state on the "lower" branch)
       PAR(2)=0.001 !D
       PAR(3)=0.169  !A
       PAR(4)=0.6    !k
       PAR(5)=1d-2   !eps
       PAR(6)=100d0   !beta
	   PAR(7)=0.01   !P
	   PAR(8)=1d0 !L
	   
		B=PAR(1)
		D=PAR(2)
		A=PAR(3)
		K=PAR(4)
		E=PAR(5)
		BETA=PAR(6)
		P=PAR(7)
		L=PAR(8)

		U(:)=0d0

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
	  
	  
!~ 	  FLUXL(2) = FLUXR(1)
!~ 	  FLUXL(3:NX-1) = -D*(1.5D0*U(3:NX-1) - 2d0*  U(2:NX-2) + 0.5D0* U(1:NX-3))/DX!SECOND ORDER BACKWARD DIFFERENCE	
!~ 	  FLUXL(NX) = FLUXR(NX-1)
	  
	  
	  ! END IMPLEMENTATION SPECIFIC CODE
	  
	  ! FLUXJUMP CONDITIONS AT THE COMPARTMNETS
	   ! FLUXL(IJUMP+1) = FLUXR(IJUMP) + FJUMP
!~ 		PRINT *, (FLUXL(IJUMP+1) - (FLUXR(IJUMP) + FJUMP))/FLUXL(IJUMP+1)
	  

	  DIVFLUX = (FLUXR-FLUXL)/DX
	  DIVFLUX(1) = 2*DIVFLUX(1)
	  DIVFLUX(NX) = 2*DIVFLUX(NX)
	  
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
      
      
      L=PAR(8)
      DX = L/DBLE(NX-1)
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
      
      SUBROUTINE PVLS(NDIM,U,PAR)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      LOGICAL, SAVE :: ifrst = .TRUE.

! Problem-independent initialization :
      IF(ifrst)THEN
         CALL INITGRID(NDIM, PAR)
         ifrst=.FALSE.
      ENDIF

      END SUBROUTINE PVLS
