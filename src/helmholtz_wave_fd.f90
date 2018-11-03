program helmholtz_wave_fd
 use helmholtz_parameters
 implicit none

 INTEGER,        PARAMETER      :: N = 100
 REAL(KIND=dpp), PARAMETER      :: Tend = 2.0d0/omega
 
 REAL(KIND=dpp), ALLOCATABLE :: x(:),u(:),rho(:),ut(:)
 REAL(KIND=dpp), ALLOCATABLE :: utt(:),v(:),c(:),f(:),up(:)
 REAL(KIND=dpp), ALLOCATABLE :: un(:), sol(:)
 
 REAL(KIND=dpp)     			:: dx, dt, dt2, t
 REAL(KIND=dpp)     			:: a, b, CFL, g
 INTEGER                        :: i, j, nsteps
 CHARACTER(len=2)               :: BC

 !! Set up domain
 a      = -1.0d0 !left endpt
 b      = 1.0d0  !rt endpt
 dx     = (b-a)/REAL(N-1,dpp)
 CFL    = 0.4d0
 dt     = CFL*dx
 dt2    = dt**2.d0
 nsteps = INT(Tend/dt)
 t      = 0.d0
 BC     = 'NL'    ! Boundary Condition

 !! Allocate working arrays
 ALLOCATE(x(N+2), u(N+2), rho(N+2), ut(N+2), utt(N+2))
 ALLOCATE(v(N+2), c(N+2), f(N), up(N+2), un(N+2), sol(N+2))

 !! spatial grid (include ghost points)
 do i = -1,N
  x(i+2) = a + REAL(i,dpp)*dx
 end do

 !! get initial conditions and physical constants
 CALL initial_condition(x,N+2,v,ut)
 CALL medium_speed(x,N+2,c)
 CALL density(x,N+2,rho)
 CALL spatial_forcing(x,N+2,f)
 CALL temporal_forcing(x,t,N+2,g)

 !! Calculate RHS of equation at t=0
 CALL compute_laplacian(v,c,dx,N,utt)

 !! Calculate u_{tt} at t=0 
 do i = 1,N
  utt(i) = (utt(i) + f(i)*g)/rho(i)
 end do

 !! Via a Taylor expansion, we note that 
 !! u(-dt) = v - dt*ut 0.5*dt^2*u_{tt}
 !! FMG: We may remove the dt*ut term 
 !! for Helmholtz solve
 up = v + 0.5d0*utt*dt2 - dt*ut
 u = v

 !! March forward in time. Here we use 
 !! a central difference in time 
 do i = 1,nsteps

  !! Solution at next time level
  un = 2.0d0*u - up + utt*dt2
  t = t + dt
  
  !! compute quantities at next time level
  CALL temporal_forcing(x,t,N,g)
  CALL enforce_bcs(x,dx,N,BC,t,un)
  CALL compute_laplacian(un,c,dx,N,utt)
  
  !! calculate u_{tt} at current time
  do j = 1,N
   utt(j) = (utt(j) + f(j)*g)/rho(j)
  end do

  !! overwrite solutions
  up = u
  u = un

  do j = 2,N+1 
  	sol(j) = v(j)*(COS(0.25*PI*t) - SIN(0.25*PI*t))
  end do 
  !! Output error at each time step
  write(*,'((E24.16))') MAXVAL(ABS(u(2:N+1) - sol(2:N+1)))
 end do 

 !! write data to file
 open(2,file='x.txt',status='unknown')
 do i=2,N+1
  write(2,fmt='((E24.16))') x(i)
 end do
 close(2)

 open(3,file='u.txt',status='unknown')
 do i=2,N+1
  write(3,fmt='((E24.16))') u(i)
 end do
 close(3)

 !! Deallocate all working arrays
 DEALLOCATE(un)
 DEALLOCATE(up)
 DEALLOCATE(f)
 DEALLOCATE(c)
 DEALLOCATE(v)
 DEALLOCATE(utt)
 DEALLOCATE(ut)
 DEALLOCATE(rho)
 DEALLOCATE(u)
 DEALLOCATE(x)

end program helmholtz_wave_fd