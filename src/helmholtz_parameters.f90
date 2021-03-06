module helmholtz_parameters
 
 implicit none
 
 !! Define parameters of the boundary (i.e. geometry
 !! and discretization related values)

 INTEGER,          PARAMETER :: dpp   = KIND(1.d0)
 INTEGER,          PARAMETER :: N     = 130000
 INTEGER,          PARAMETER :: MAXIT = 120
 REAL(KIND=dpp),   PARAMETER :: PI    = ACOS(-1.d0)
 REAL(KIND=dpp),   PARAMETER :: OMEGA = 62777.3d0 !!wavenumber
 REAL(KIND=dpp),   PARAMETER :: Tend  = 2.0d0*PI/OMEGA
 REAL(KIND=dpp),   PARAMETER :: SIGMA = 37.d0
 REAL(KIND=dpp),   PARAMETER :: A     = -1.d0 !!left endpt
 REAL(KIND=dpp),   PARAMETER :: B     = 1.d0  !!right endpt
 REAL(KIND=dpp),   PARAMETER :: CFL   = 0.1d0
 REAL(KIND=dpp),   PARAMETER :: TOL   = 1.0d-15
 CHARACTER(len=2), PARAMETER :: BC    = 'DD'

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
  u  = 0.d0
  ut = 0.d0
  ! do i =1,N
  !  ! u(i) = COS(0.5d0*x(i)*PI)
  !  ! u(i) = EXP(-1.d0*SIGMA*x(i)**2.d0)
  !  ! u(i) = COS(0.25d0*x(i)*PI) - SIN(0.25d0*x(i)*PI)
  !  ! ut(i) = -1.d0*u(i)*PI*0.25d0
  ! end do 

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

  do i = 1,N
   ! coeff = 4.d0*(SIGMA*x(i))**2.d0 + OMEGA**2.d0 - 2.d0*SIGMA
   ! f(i) = coeff*EXP(-SIGMA*(x(i)**2.d0))
   f(i) = (OMEGA**2.d0 - (0.5d0*PI)**2.d0)*COS(0.5d0*PI*x(i))
   ! f(i) = (x(i)**2.d0 - 1.d0)*OMEGA**2.d0 + 2.d0
  end do

 end subroutine spatial_forcing

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in a vector of spatial values and 
 !! returns the temporal part of the forcing at that point
 !! Inputs:
 !! 	- x       : spatial grid
 !!     - t       : the current time
 !!     - N       : number of grid points
 !! Outputs:
 !!		- g       : forcing term
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine temporal_forcing(x,t,N,num_freq,g)
  
  implicit none 
  
  REAL(KIND=dpp), INTENT(IN)    :: x(:)
  REAL(KIND=dpp), INTENT(IN)    :: t
  INTEGER,        INTENT(IN)    :: N
  INTEGER,        INTENT(IN)    :: num_freq
  REAL(KIND=dpp), INTENT(INOUT) :: g

  !! here we choose a constant
  if(num_freq .eq. 1) then 
  	g = COS(OMEGA*t)
  else
  	g = COS(OMEGA*t) + COS(2*OMEGA*t)
  end if

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
 !!     - x       : spatial nodes
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

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in the current solution vector v
 !! and solves a wave equation (assuming 0 initial vel.)
 !! and projects it to the frequency OMEGA to solve 
 !! Helmholtz.
 !! Inputs:
 !!     - N        : number of interior grid points
 !!     - dx       : spatial step size
 !!     - end_time : time to solve wave equation to
 !!     - freq     : frequencies we wish to project to
 !!     - num_freq : number of frequencies we project to
 !! 	- c        : speed of sound in medium
 !!		- f        : forcing term
 !!		- rho      : density
 !!     - x        : spatial nodes
 !!     - v        : intial condition for wave solve
 !! Outputs:
 !! 	- v       : projected wave solution
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine evolve_and_project(N,dx,end_time,freq, &
 							   num_freq,c,f,rho,x,v)
  
  implicit none

  INTEGER,        INTENT(IN)    :: N
  REAL(KIND=dpp), INTENT(IN)    :: dx
  REAL(KIND=dpp), INTENT(IN)    :: end_time
  REAL(KIND=dpp), INTENT(IN)    :: freq(:)
  INTEGER,        INTENT(IN)    :: num_freq
  REAL(KIND=dpp), INTENT(IN)    :: c(:)
  REAL(KIND=dpp), INTENT(IN)    :: f(:)
  REAL(KIND=dpp), INTENT(IN)    :: rho(:)
  REAL(KIND=dpp), INTENT(IN)    :: x(:)
  REAL(KIND=dpp), INTENT(INOUT) :: v(:)

  REAL(KIND=dpp) :: t, g, dt, dt2
  INTEGER        :: i, j, nsteps

  REAL(KIND=dpp), ALLOCATABLE :: utt(:), up(:), u(:), un(:)

  t      = 0.d0
  dt     = CFL*dx
  nsteps = floor(end_time/dt) + 1
  dt     = end_time/dble(nsteps)
  dt2    = dt**2.d0

  !! Allocate working arrays
  ALLOCATE(u(N+2),utt(N+2),up(N+2),un(N+2))

  !! Calculate RHS of equation at t=0
  CALL compute_laplacian(v,c,dx,N,utt)

  !! Calculate u_{tt} at t=0 
  CALL temporal_forcing(x,t,N+2,num_freq,g)
  do i = 1,N
   utt(i) = (utt(i) - f(i)*g)/rho(i)
  end do

  !! Via a Taylor expansion, we note that 
  !! u(-dt) = v - dt*ut + 0.5*dt^2*u_{tt}
  !! FMG: We may remove the dt*ut term 
  !! for Helmholtz solve
  ! up = v + 0.5d0*utt*dt2 - dt*ut
  ! up = v + 0.5d0*utt*dt2 + v*(dt2*OMEGA**2.d0)**2.d0
  up = v + 0.5d0*utt*dt2

  up = v
  u  = v

  !! March forward in time. Here we use 
  !! a central difference in time 
  do i = 1,nsteps

   !! Solution at next time level
   un = 2.0d0*u - up + utt*dt2
   t = t + dt

   !! compute quantities at next time level
   CALL temporal_forcing(x,t,N+2,num_freq,g)
   CALL enforce_bcs(x,dx,N,BC,t,un)
   CALL compute_laplacian(un,c,dx,N,utt)

   !! calculate u_{tt} at current time
   do j = 1,N
    utt(j) = (utt(j) - f(j)*g)/rho(j)
   end do

   !! overwrite solutions
   up = u
   u  = un

   !! collect t = 0
   if(i .eq. 1) then 
    ! v = 0.375d0*dt*v*num_freq
    v = 0.375d0*dt*v*num_freq
   end if

   !! collect projection for each frequency
   do j = 1,num_freq
    v = v + dt*(COS(freq(j)*t)-0.25d0/num_freq)*u
   end do

  end do

  !! do a correction for trap rule
  do j = 1,num_freq
    v = v - 0.5d0*dt*(COS(freq(j)*t)-0.25d0/num_freq)*u
  end do 
  v = 2.d0*v/end_time

  !! Deallocate working arrays
  DEALLOCATE(un)
  DEALLOCATE(u)
  DEALLOCATE(up)
  DEALLOCATE(utt)
 end subroutine evolve_and_project


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in an initial condition
 !! and solves a wave equation (assuming 0 initial vel.)
 !! Inputs:
 !!     - N        : number of interior grid points
 !!     - dx       : spatial step size
 !!     - end_time : time to solve wave equation to
 !! 	- c        : speed of sound in medium
 !!		- f        : forcing term
 !!		- rho      : density
 !!     - x        : spatial nodes
 !!     - v        : intial condition for wave solve
 !! Outputs:
 !! 	- v       : wave equation solution at time end_time
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine wave_solve(N,num_freq,dx,end_time,c,f,rho,x,v)
  implicit none

  INTEGER,        INTENT(IN)    :: N
  INTEGER,        INTENT(IN)    :: num_freq
  REAL(KIND=dpp), INTENT(IN)    :: dx
  REAL(KIND=dpp), INTENT(IN)    :: end_time
  REAL(KIND=dpp), INTENT(IN)    :: c(:)
  REAL(KIND=dpp), INTENT(IN)    :: f(:)
  REAL(KIND=dpp), INTENT(IN)    :: rho(:)
  REAL(KIND=dpp), INTENT(IN)    :: x(:)
  REAL(KIND=dpp), INTENT(INOUT) :: v(:)

  REAL(KIND=dpp) :: t, g, dt, dt2
  INTEGER        :: i, j, nsteps

  REAL(KIND=dpp), ALLOCATABLE :: utt(:), up(:), u(:), un(:)

  t        = 0.d0
  dt       = CFL*dx
  nsteps   = floor(end_time/dt) + 1
  dt       = end_time/dble(nsteps)
  dt2      = dt**2.d0

  !! Allocate working arrays
  ALLOCATE(u(N+2),utt(N+2),up(N+2),un(N+2))

  !! Calculate RHS of equation at t=0
  CALL compute_laplacian(v,c,dx,N,utt)

  !! Calculate u_{tt} at t=0 
  CALL temporal_forcing(x,t,N+2,num_freq,g)
  do i = 1,N
   utt(i) = (utt(i) - f(i)*g)/rho(i)
  end do

  !! Via a Taylor expansion, we note that 
  !! u(-dt) = v - dt*ut + 0.5*dt^2*u_{tt}
  !! FMG: We may remove the dt*ut term 
  !! for Helmholtz solve
  up = v + 0.5d0*utt*dt2
  ! up = v
  u  = v

  !! March forward in time. Here we use 
  !! a central difference in time 
  do i = 1,nsteps

   !! Solution at next time level
   un = 2.0d0*u - up + utt*dt2
   t = t + dt

   !! compute quantities at next time level
   CALL temporal_forcing(x,t,N+2,num_freq,g)
   CALL enforce_bcs(x,dx,N,BC,t,un)
   CALL compute_laplacian(un,c,dx,N,utt)

   !! calculate u_{tt} at current time
   do j = 1,N
    utt(j) = (utt(j) - f(j)*g)/rho(j)
   end do

   !! overwrite solutions
   up = u
   u  = un
  end do
  v = u
  !! Deallocate working arrays
  DEALLOCATE(un)
  DEALLOCATE(u)
  DEALLOCATE(up)
  DEALLOCATE(utt)
 end subroutine wave_solve

end module helmholtz_parameters