program helmholtz_wave_fd
 use helmholtz_parameters
 implicit none

 INTEGER,        PARAMETER :: N = 100
 REAL(KIND=dpp), PARAMETER :: Tend = 2.0d0/omega
 
 REAL(KIND=dpp)     :: x(N), u(N), rho(N), ut(N), utt(N)
 REAL(KIND=dpp)     :: v(N), c(N), f(N), up(N), un(N)
 REAL(KIND=dpp)     :: dx, dt, dt2, t
 REAL(KIND=dpp)     :: a, b, CFL, g
 INTEGER            :: i, j, nsteps
 CHARACTER(len=1)   :: BC

 !! Set up domain
 a      = -1.0d0 !left endpt
 b      = 1.0d0  !rt endpt
 dx     = (b-a)/REAL(N-1,dpp)
 CFL    = 0.4d0
 dt     = CFL*dx
 dt2    = dt**2.d0
 nsteps = INT(Tend/dt)
 t      = 0.d0
 BC     = 'D'    !Dirichlet boundary condition

 !! spatial grid
 do i = 0,N-1
  x(i+1) = a + REAL(i,dpp)*dx
 end do

 !! get initial conditions and physical constants
 CALL initial_condition(x,N,v)
 CALL medium_speed(x,N,c)
 CALL density(x,N,rho)
 CALL spatial_forcing(x,N,f)
 CALL temporal_forcing(x,t,N,g)

 !! Calculate RHS of equation at t=0
 CALL compute_laplacian(v,c,dx,N,BC,utt)

 !! Calculate u_{tt} at t=0 
 do i = 1,N
  utt(i) = (utt(i) + f(i)*g)/rho(i)
 end do

 !! Via a Taylor expansion, we note that 
 !! u(-dt) = v + 0.5*dt^2*u_{tt}
 up = v + 0.5d0*utt*dt2
 u = v

 !! March forward in time. Here we use 
 !! a central difference in time 
 do i = 1,nsteps

  !! Solution at next time level
  un = 2.0d0*u - up + utt*dt2
  t = t + dt
  
  !! compute quantities at next time level
  CALL temporal_forcing(x,t,N,g)
  CALL compute_laplacian(un,c,dx,N,BC,utt)
  
  !! calculate u_{tt} at current time
  do j = 1,N
   utt(j) = (utt(j) + f(j)*g)/rho(j)
  end do

  !! overwrite solutions
  up = u
  u = un
 end do 

 !! write data to file
 open(2,file='x.txt',status='unknown')
 do i=1,N
  write(2,fmt='((E24.16))') x(i)
 end do
 close(2)

 open(3,file='u.txt',status='unknown')
 do i=1,N
  write(3,fmt='((E24.16))') u(i)
 end do
 close(3)


end program helmholtz_wave_fd