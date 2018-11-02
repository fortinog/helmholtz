module helmholtz_parameters
 
 implicit none
 
 INTEGER,        PARAMETER :: dpp   = KIND(1.d0)
 REAL(KIND=dpp), PARAMETER :: omega = 10.d0  !wavenumber
 REAL(KIND=dpp), PARAMETER :: PI    = ACOS(-1.d0)

 contains


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
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine initial_condition(x,N,u)
  
  implicit none 
  
  REAL(KIND=dpp), intent(in)    :: x(:)
  INTEGER,        intent(in)    :: N
  REAL(KIND=dpp), intent(inout) :: u(:)

  INTEGER :: i
  !! here we choose a constant
  ! u = 0.d0
  do i =1,N
   u(i) = COS(0.5d0*x(i)*PI)
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
  
  REAL(KIND=dpp), intent(in)    :: x(:)
  INTEGER,        intent(in)    :: N
  REAL(KIND=dpp), intent(inout) :: c(:)

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
  
  REAL(KIND=dpp), intent(in)    :: x(:)
  INTEGER,        intent(in)    :: N
  REAL(KIND=dpp), intent(inout) :: rho(:)

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
  
  REAL(KIND=dpp), intent(in)    :: x(:)
  INTEGER,        intent(in)    :: N
  REAL(KIND=dpp), intent(inout) :: f(:)

  !! here we choose a constant
  f = 0.d0
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
  
  REAL(KIND=dpp), intent(in)    :: x(:)
  REAL(KIND=dpp), intent(in)    :: t
  INTEGER,        intent(in)    :: N
  REAL(KIND=dpp), intent(inout) :: g

  !! here we choose a constant
  ! g = COS(omega*t)
  g = 0.d0
 end subroutine temporal_forcing

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in the current solution vector u
 !! and returns the Laplacian portion of the PDE.
 !! Inputs:
 !! 	- u       : current solution
 !! 	- c       : speed of sound in medium
 !!     - dx      : spatial step size
 !!     - N       : number of grid points
 !!     - BC      : 'D' for Dirichlet, 'N' for Neumann BCs
 !! Outputs:
 !!		- lap     : the laplacian 
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine compute_laplacian(u,c,dx,N,BC,lap)

  implicit none
  
  REAL(KIND=dpp),   intent(in)    :: u(:)
  REAL(KIND=dpp),   intent(in)    :: c(:)
  REAL(KIND=dpp),   intent(in)    :: dx
  INTEGER,          intent(in)    :: N
  CHARACTER(len=*), intent(in)    :: BC
  REAL(KIND=dpp),   intent(inout) :: lap(:)

  INTEGER        :: i
  REAL(KIND=dpp) :: denom

  denom = 0.5d0/(dx**2.0d0)

  if(BC .eq. 'D') then 

   !! compute laplacian at left boundary
   i = 2
   lap(1) = SUM(c(i:i+1))*u(i+1) - SUM(c(i-1:i+1))*u(i)
   lap(1) = lap(1) - c(i)*u(i)
   lap(1) = denom*lap(1)

   !! Account for Dirichlet boundary conditions
   lap(i) = SUM(c(i:i+1))*u(i+1) - SUM(c(i-1:i+1))*u(i)
   lap(i) = lap(i) - c(i)*u(i)
   lap(i) = denom*lap(i)

   !! compute laplacian on nodes with support 
   !! contained in the interior of the domain
   do i=3,N-2
    lap(i) = SUM(c(i:i+1))*u(i+1) - SUM(c(i-1:i+1))*u(i)
    lap(i) = lap(i) - c(i)*u(i) + SUM(c(i-1:i))*u(i-1)
    lap(i) = denom*lap(i)
   end do 

   i = N-1
   !! Account for Dirichlet boundary conditions
   lap(i) = -1.d0*SUM(c(i-1:i+1))*u(i)
   lap(i) = lap(i) - c(i)*u(i) + SUM(c(i-1:i))*u(i-1)
   lap(i) = denom*lap(i)

   !! compute laplacian at right boundary
   lap(N) = -1.d0*SUM(c(i-1:i+1))*u(i)
   lap(N) = lap(N) - c(i)*u(i) + SUM(c(i-1:i))*u(i-1)
   lap(N) = denom*lap(N)

  elseif(BC .eq. 'N') then 
   
   write(*,*) "This boundary condition is not supported."
   stop 234
  
  else 
   
   write(*,*) "This boundary condition is not supported."
   stop 123
  
  end if

 end subroutine compute_laplacian


end module helmholtz_parameters