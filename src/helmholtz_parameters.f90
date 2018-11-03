module helmholtz_parameters
 
 implicit none
 
 INTEGER,        PARAMETER :: dpp   = KIND(1.d0)
 REAL(KIND=dpp), PARAMETER :: OMEGA = 10.d0  !wavenumber
 REAL(KIND=dpp), PARAMETER :: PI    = ACOS(-1.d0)
 REAL(KIND=dpp), PARAMETER :: SIGMA = 37.d0
 contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in a boundary node and specifies
 !! the Dirichlet boundary condition at time t. 
 !! Inputs:
 !! 	- x       : boundary node
 !!     - t       : current time
 !! Outputs:
 !!		- u       : boundary value
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine dirichlet_bc(x,t,u)

  implicit none 

  REAL(KIND=dpp), INTENT(IN)    :: x
  REAL(KIND=dpp), INTENT(IN)    :: t
  REAL(KIND=dpp), INTENT(INOUT) :: u

  u = 0.d0
 end subroutine dirichlet_bc

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in a boundary node and specifies
 !! the Neumann boundary condition at time t. 
 !! Inputs:
 !! 	- x       : boundary node
 !!     - t       : current time
 !! Outputs:
 !!		- u       : boundary value
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine neumann_bc(x,t,u)

  implicit none 

  REAL(KIND=dpp), INTENT(IN)    :: x
  REAL(KIND=dpp), INTENT(IN)    :: t
  REAL(KIND=dpp), INTENT(INOUT) :: u

  u = 0.d0
 end subroutine neumann_bc

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in a vector of spatial values and 
 !! returns the desired initial condition for the wave eq.
 !! (We are implicitly assuming we do not need an initial 
 !! velocity).
 !! Inputs:
 !! 	- x       : spatial grid
 !!     - N       : number of grid points
 !! Outputs:
 !!		- u       : initial condition to a wave equation
 !!     - ut      : initial velocity to a wave equation
 !!				    FMG remove this for helmholtz solve?
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine initial_condition(x,N,u,ut)
  
  implicit none 
  
  REAL(KIND=dpp), INTENT(IN)    :: x(:)
  INTEGER,        INTENT(IN)    :: N
  REAL(KIND=dpp), INTENT(INOUT) :: u(:)
  REAL(KIND=dpp), INTENT(INOUT) :: ut(:)

  INTEGER :: i
  !! here we choose a constant
  ! u = 0.d0
  do i =1,N
   ! u(i) = COS(0.5d0*x(i)*PI)
   ! u(i) = EXP(-SIGMA*x(i)**2.d0)
   u(i) = COS(0.25d0*x(i)*PI) - SIN(0.25d0*x(i)*PI)
   ut(i) = -1.d0*u(i)*PI*0.25d0
  end do 

 end subroutine initial_condition

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in a vector of spatial values and 
 !! returns the speed of the medium c
 !! Inputs:
 !! 	- x       : spatial grid
 !!     - N       : number of grid points
 !! Outputs:
 !!		- c       : speed of the medium
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine medium_speed(x,N,c)
  
  implicit none 
  
  REAL(KIND=dpp), INTENT(IN)    :: x(:)
  INTEGER,        INTENT(IN)    :: N
  REAL(KIND=dpp), INTENT(INOUT) :: c(:)

  !! here we choose a constant
  c = 1.d0
 end subroutine medium_speed

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in a vector of spatial values and 
 !! returns the density at that point
 !! Inputs:
 !! 	- x       : spatial grid
 !!     - N       : number of grid points
 !! Outputs:
 !!		- rho     : density
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine density(x,N,rho)
  
  implicit none 
  
  REAL(KIND=dpp), INTENT(IN)    :: x(:)
  INTEGER,        INTENT(IN)    :: N
  REAL(KIND=dpp), INTENT(INOUT) :: rho(:)

  !! here we choose a constant
  rho = 1.d0
 end subroutine density

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in a vector of spatial values and 
 !! returns the spatial part of the forcing at that point
 !! Inputs:
 !! 	- x       : spatial grid
 !!     - N       : number of grid points
 !! Outputs:
 !!		- f       : forcing term
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine spatial_forcing(x,N,f)
  
  implicit none 
  
  REAL(KIND=dpp), INTENT(IN)    :: x(:)
  INTEGER,        INTENT(IN)    :: N
  REAL(KIND=dpp), INTENT(INOUT) :: f(:)

  INTEGER        :: i
  REAL(KIND=dpp) :: coeff

  !! here we choose a constant
  f = 0.d0

  ! do i = 1,N
  !  coeff = 4.d0*(SIGMA*x(i))**2.d0 + OMEGA**2.d0 - 2.d0*SIGMA
  !  f(i) = coeff*EXP(-SIGMA*x(i)**2.d0)
  ! end do

 end subroutine spatial_forcing

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in a vector of spatial values and 
 !! returns the temporal part of the forcing at that point
 !! Inputs:
 !! 	- x       : spatial grid
 !!     - N       : number of grid points
 !! Outputs:
 !!		- g       : forcing term
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine temporal_forcing(x,t,N,g)
  
  implicit none 
  
  REAL(KIND=dpp), INTENT(IN)    :: x(:)
  REAL(KIND=dpp), INTENT(IN)    :: t
  INTEGER,        INTENT(IN)    :: N
  REAL(KIND=dpp), INTENT(INOUT) :: g

  !! here we choose a constant
  ! g = COS(omega*t)
  g = 1.d0
 end subroutine temporal_forcing

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in the current solution vector u
 !! and returns the Laplacian portion of the PDE.
 !! Inputs:
 !! 	- u       : current solution
 !! 	- c       : speed of sound in medium
 !!     - dx      : spatial step size
 !!     - N       : number of grid points
 !! Outputs:
 !!		- lap     : the laplacian 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine compute_laplacian(u,c,dx,N,lap)

  implicit none
  
  REAL(KIND=dpp),   INTENT(IN)    :: u(:)
  REAL(KIND=dpp),   INTENT(IN)    :: c(:)
  REAL(KIND=dpp),   INTENT(IN)    :: dx
  INTEGER,          INTENT(IN)    :: N
  REAL(KIND=dpp),   INTENT(INOUT) :: lap(:)

  INTEGER        :: i
  REAL(KIND=dpp) :: denom

  denom = 0.5d0/(dx**2.0d0)

  lap(1) = 0.d0
  do i=2,N+1
   lap(i) = SUM(c(i:i+1))*u(i+1) - SUM(c(i-1:i+1))*u(i)
   lap(i) = lap(i) - c(i)*u(i) + SUM(c(i-1:i))*u(i-1)
   lap(i) = denom*lap(i)
  end do
  lap(N+2) = 0.d0 

 end subroutine compute_laplacian

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in the current solution vector u
 !! and enforces the desired boundary conditions (this
 !! is to be modified if different conditions are desired).
 !! Inputs:
 !!     - dx      : spatial step size
 !!     - N       : number of interior grid points
 !!     - BC      : 'D' for Dirichlet, 'N' for Neumann BCs
 !!     - t       : currrent time
 !! Outputs:
 !! 	- u       : current solution
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine enforce_bcs(x,dx,N,BC,t,u)
  
  implicit none 
  
  REAL(KIND=dpp),   INTENT(IN)    :: x(:)
  REAL(KIND=dpp),   INTENT(IN)    :: dx
  INTEGER,          INTENT(IN)    :: N
  CHARACTER(len=*), INTENT(IN)    :: BC
  REAL(KIND=dpp),   INTENT(IN)    :: t
  REAL(KIND=dpp),   INTENT(INOUT) :: u(:)

  INTEGER        :: i
  REAL(KIND=dpp) :: val

  SELECT CASE (BC)
  	CASE('DD')
  		!! Purely Dirichlet data
  		CALL dirichlet_bc(x(2),t,val)
  		u(2)   = val
  		CALL dirichlet_bc(x(N+1),t,val)
  		u(N+1) = val
  	CASE('NL')
  		!! Neumann condition on left
  		CALL neumann_bc(x(2),t,val)
  		u(1)   = u(3) - 2.d0*dx*val
   		CALL dirichlet_bc(x(N+1),t,val)
  		u(N+1) = val
  	CASE('NR')
  		!! Neumann condition on right
   		CALL dirichlet_bc(x(2),t,val)
  		u(2)   = val
  		CALL neumann_bc(x(N+1),t,val)
  		u(N+2) = u(N) + 2.d0*dx*val
  	END SELECT
 end subroutine enforce_bcs

end module helmholtz_parameters