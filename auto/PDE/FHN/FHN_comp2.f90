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
        INTEGER, PARAMETER :: ORDER = 1 ! ORDER FOR FIRST DIFFERENCE SCHEME (BETWEEN 1 AND 3). FOR DIFFUSION, WORST-CASE ACCURACY OF THE NUMERICS IS ORDER+1 (NEAR COMPARTMENTS) BUT AS HIGH AS 2*ORDER AWAY FROM COMPARTMENTS.
		INTEGER :: NX
		INTEGER :: ORDERINT
		DOUBLE PRECISION :: DX
		DOUBLE PRECISION :: DX2
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DI
        INTEGER, ALLOCATABLE, DIMENSION(:) :: IJUMP ! WE HAVE JUMPS AT THE INTERFACE BETWEEN CELLS I AND I+1
        INTEGER, ALLOCATABLE, DIMENSION(:) :: BP !POSITION OF CELL INTERFACES WITH BOUNDARY CONDITIONS (IJUMP + ENDPOINTS)
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: STENFOR, STENBACK, STEN
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::  STENINT
        DOUBLE PRECISION  :: UP1(2) = (/ -1.0, 1.0 /) !FIRST ORDER UPWIND SCHEME
        DOUBLE PRECISION  :: UP2(3) = (/  -1.5d0,  2d0,   	-0.5d0 /) !SECOND ORDER UPWIND SCHEME
        DOUBLE PRECISION  :: UP3(4) = (/ -11d0/6d0,	 3d0,	-1.5d0,	1d0/3d0 /) !THIRD ORDER UPWIND SCHEME
        DOUBLE PRECISION  :: UP4(5) = (/ -25d0/12d0,	4d0,	-3d0,	4d0/3d0,	-1d0/4d0 /) !FOURTH ORDER UPWIND SCHEME
        DOUBLE PRECISION  :: UP5(6) = (/ -137d0/60d0,	5d0,	-5d0,	10d0/3d0,	-5d0/4d0,	1d0/5d0 /) !FIFTH ORDER UPWIND SCHEME
		DOUBLE PRECISION  :: LAP2(3) = (/ 1d0,-2d0,1d0 /) 
        DOUBLE PRECISION  :: LAP4(5) = (/ -1d0/12d0,	4d0/3d0,	-5d0/2d0,	4d0/3d0,	-1d0/12d0 /) 

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
		
		
		COMPFLUX = BETA*(V-(I(IJUMP)+I(IJUMP+1))/2d0)

		
		F(1:NCOMP) = U*(1d0-U)*(U-A)-V
		F(NCOMP+1:2*NCOMP) = E*(K*U-V-B) - COMPFLUX
		
		
		
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
		   !!!!! THIS NEEDS TO BE REDONE COMPLETELY, BUT SHOULD BE "SIMPLER" SOMEHOW WITH THE STENCILS
		   DFDU(PDE, PDE) = -D/DX2 - p
		   DFDU(PDE, PDE+1) = D/DX2
		   DO J = 1, NX-2
			   DFDU(PDE+J, PDE+J-1) = D/DX2
			   DFDU(PDE+J, PDE+J) = -2*D/DX2 - p
			   DFDU(PDE+J, PDE+J+1) = D/DX2
		   END DO
		   DFDU(PDE+NX-1, PDE+NX-1) = -D/DX2 - p
		   DFDU(PDE+NX-1, PDE+NX-2) = D/DX2
		 !!!!! THIS WILL NEED TWEAKING...BUT IS NOT COMPLETELY WRONG...PROBABLY CHANGE TO A FACTOR OF BETA/2d0
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
       PAR(2)=0.0012 !D
       PAR(3)=0.139  !A
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
	  DOUBLE PRECISION  FLUXR(NX), FLUXL(NX), UINT(NCOMP), D2U(NX)
	  INTEGER I, J, J0

     ! J0=BP(I-1)+1-------------------BP(I)
	
	  DO I=2,NCOMP+2
	  J0=BP(I-1)+1
!~ 	  DO J=1,ORDER-1
!~ 		FLUXL(J0+J) = D*SUM(STEN(:,J)*U(J0+J:J0:-1))/DX
!~ 	  END DO
	  
!~ 	  DO J=J0+ORDER,BP(I)
!~ 		FLUXL(J) = D*SUM(STEN(:,ORDER)*U(J:J-ORDER:-1))/DX
!~ 	  END DO


	  
	  DO J=J0,BP(I)-ORDER
	   FLUXR(J) = -D*SUM(STEN(:,ORDER)*U(J:J+ORDER))/DX
	  END DO
	  
	  DO J=ORDER-1,1,-1
	  FLUXR(BP(I)-J) = -D*SUM(STEN(:,J)*U(BP(I)-J:BP(I)))/DX
	  END DO
	  
	  IF (I.GT.2) THEN
      FLUXR(J0-1) = -D*SUM(STENINT*U(J0-1:J0-1-ORDERINT:-1))/DX + D*SUM(STENINT*U(J0:J0+ORDERINT))/DX

	  END IF
	  
	  END DO
	  

	  	  
	    !NEUMANN BOUNDARY CONDITIONS
	  FLUXL(1) = 0d0 
	  FLUXR(NX) = 0d0
	  
	  FLUXL(2:NX)=FLUXR(1:NX-1)
	  
	  
	  D2U(2)=SUM(LAP2*U(1:3))/DX2
	  D2U(NX-1)=SUM(LAP2*U(NX-2:NX))/DX2
	  DO I=3,Nx-2
	  D2U(I)=SUM(LAP4*U(I-2:I+2))/DX2
	  END DO
	  
	  
!~ 	  PRINT *, FLUXR

	  
	  
	  !!!!! NEED SOMETHING HERE TO DEAL WITH THE BIDIRECTIONALITY OF THE FLUX DEPENDENCE AT THE COMPARTMENT BOUNDARIES......
	  DO I=1,NCOMP
	  J0=IJUMP(I)+1
	  
	  END DO
	  
      FLUXL(IJUMP+1) = FLUXL(IJUMP+1) + FJUMP/2d0 
	  FLUXR(IJUMP) = FLUXR(IJUMP) - FJUMP/2d0

	  
!~ 	  DIVFLUX = -D * D2U
!~       DIVFLUX(IJUMP) = (FLUXR(IJUMP)-FLUXL(IJUMP))/DX
!~       DIVFLUX(NX) = (FLUXR(NX)-FLUXL(NX))/DX
!~       DIVFLUX(1) = (FLUXR(1)-FLUXL(1))/DX
      
      
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
          
      PRINT *, NDIM, NCOMP, NODE, NPDE
      TEMP = DBLE(NDIM - (NCOMP*NODE)) / DBLE(NPDE)
      
      IF (MOD(TEMP, 1d0) .NE. 0D0) THEN
		  PRINT *, "error: NDIM cannot be split into an integer number of grid points for the PDEs"
		  call EXIT(1)
      ELSE
		NX = INT(TEMP)
      END IF
      
      ALLOCATE(IJUMP(NCOMP))
      ALLOCATE(BP(NCOMP+2))
      
      
      L=PAR(8)
      DX = L/DBLE(NX-1)
      DX2 = DX*DX
      H = L/DBLE(NCOMP)
      

      SPACE = DBLE(NX)/DBLE(NCOMP)
      
      DO I=1, NCOMP
	      IJUMP(I) = NINT((I-0.5)*SPACE)
      END DO
      
      BP(1)=0
      BP(2:NCOMP+1)=IJUMP
      BP(NCOMP+2)=NX
      
      PRINT *, NX
      PRINT *,"SPACING: ", SPACE
      
      IF (MOD(SPACE,1.0).NE.0d0) THEN
		PRINT *, "WARNING: COMPARTMENTS CANNOT BE EVENLY SPREAD THROUGH THE GRID"
		PRINT *, "SUGGESTED NDIM:"
		PRINT *, NDIM + CEILING(SPACE)*NCOMP-NX
      END IF
      
      PRINT *, NX
      PRINT *, IJUMP
      
      
	  ALLOCATE(STEN(ORDER+1, ORDER))
	  
	  ! THIS CRAP SHOULD BE MOVED INTO AN INITSTEN() SUBROUTINE
	  
	  STEN(:,:)=0d0
	  STEN(1:2,1)=UP1
	  IF (ORDER.GT.1) THEN
		STEN(1:3,2)=UP2
	  END IF
	  IF (ORDER.GT.2) THEN
		STEN(1:4,3)=UP3
	  END IF
	  IF (ORDER.GT.3) THEN
		STEN(1:5,4)=UP4
	  END IF
	  IF (ORDER.GT.4) THEN
		STEN(1:6,5)=UP5
	  END IF
	  
	  DO I=1,ORDER
	  PRINT *,"STEN(:",I,")", STEN(:,I)
	  END DO
	  
	  !STENCIL FOR INTERFACE BOUNDARIES DERIVATIVES
	  ORDERINT = CEILING(DBLE(ORDER)/2D0)+1-MOD(ORDER,2)
!~ 	  ORDERINT=2
	  ALLOCATE(STENINT(ORDERINT))
	  STENINT(:)=0d0
!~ 	  STENINT(1)=1d0
!~ 	  STENINT(2)=-1d0

	  STENINT(1)=-STEN(2, ORDER)
	  DO I=3,ORDER+1
	  J=CEILING(DBLE(I-1)/2D0)
	  IF (MOD(I,2).NE.0) THEN !IF I IS ODD
	      STENINT(J+1)=STENINT(J+1)-STEN(I, ORDER)/2d0
	      STENINT(J)=STENINT(J)-STEN(I, ORDER)/2d0
	  ELSE
	     STENINT(J)=STENINT(J)-STEN(I, ORDER)
	  END IF
	  END DO
	  PRINT *, "STENINT", STENINT
	  
	   
	   
	   
	   
	   
!~ 	  ALLOCATE(STENFOR(ORDER+1, NX))
      
	   
!~ 	   PRINT *, BP
	  
!~ 	  DO I=2,NCOMP+2
!~ 	  DO J=BP(I-1)+1,BP(I)-ORDER
!~ 		STENFOR(1:ORDER+1,J)=STEN(:,ORDER)
!~ 	  END DO
!~ 	  DO J=ORDER-1,1,-1
!~ 		STENFOR(1:J+1,BP(I)-J)=STEN(:,J)
!~ 	  END DO
!~ 	  STENFOR(:,BP(I)) = STENINT
!~ 	  END DO
!~ 	  STENFOR(:,NX)=0d0 !NEUMANN BOUNDARIES
	  
!~ 	  ALLOCATE(STENBACK(ORDER+1, NX))
	  
	  
!~ 	  DO I=2,NCOMP+2
	  
!~ 	  DO J=1,ORDER-1
!~ 		STENBACK(1:J+1,BP(I-1)+1+J)=-STEN(:,J)
!~ 	  END DO
	  
!~ 	  DO J=BP(I-1)+ORDER,BP(I)
!~ 		STENBACK(1:ORDER+1,J)=-STEN(:,ORDER)
!~ 	  END DO

!~ 	  STENBACK(:,BP(I-1)+1) = -STENINT
!~ 	  END DO
!~ 	  STENBACK(:,1)=0d0 !NEUMANN BOUNDARIES
!~ 	  DO I=1,NX
!~ 	    PRINT *, I, "::", STENBACK(:,I)
!~ 	  END DO
	  
!~      call EXIT(1)
     
     
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

subroutine weights (z,x,n,nd,m,c)
!~  c---
!~  c--- Input Parameters
!~  c--- z location where approximations are to be accurate,
!~  c--- x(O:nd) grid point locations, found in x(O:n)
!~  c--- n one less than total number of grid points; n must
!~  c--- not exceed the parameter nd below,
!~  c--- nd dimension of x- and c-arrays in calling program
!~  c--- x(O:nd) and c(O:nd,O:m), respectively,
!~  c--- m highest derivative for which weights are sought,
!~  c--- Output Parameter
!~  c--- c(O:nd,O:m) weights at grid locations x(O:n) for derivatives
!~  c--- of order O:m, found in c(O:n,O:m)
!~  c---
 implicit real*8 (a-h,o-z)
 dimension x(nd+1),c(nd+1,m+1)
 INTEGER I,J,K, N, c1, c2, c3, c4
	 cl = 1d0
	 c4 = x(1)-z
	 DO k=1,m+1
	 DO j=1,n+1
		c(j,k) = 0d0
	 END DO
	 END DO
	 
	 c(1,1) = 1.0d0
	 do i=2,n+1
	 mn = min(i,m)
	 c2 = 1.0d0
	 c5 = c4
	  c4 = x(i)-z
	 do j=1,i-l
		 c3 = x(i)-x(j)
		 c2 = c2*c3
		 if (j.eq.i-1) then
		 do  k=mn,I,-I
		  c(i,k) = cl*(k*c(i-l,k-1)-c5*c(i-l,k))/c2
		 END DO
		 c(i,1) = -cl*c5*c(i-1,1)/c2
		 endif
			 do  k=mn,I,-I
			  c(j,k) = (c4*c(j,k)-k*c(j,k-1))/c3
			 END DO
		  c(j,1) = c4*c(j,1)/c3
	 END DO
	  cl = c2
	 END DO
	 return
 end
