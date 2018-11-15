program helmholtz_wave_fd
 
 use helmholtz_parameters
 implicit none
 
 REAL(KIND=dpp), ALLOCATABLE :: x(:),u(:),rho(:),ut(:)
 REAL(KIND=dpp), ALLOCATABLE :: utt(:),v(:),c(:),f(:),up(:)
 REAL(KIND=dpp), ALLOCATABLE :: un(:),vp(:),b_vec(:),v_tmp(:)
 REAL(KIND=dpp), ALLOCATABLE :: u_exact(:)
 REAL(KIND=dpp)     		 :: dx,dt,dt2,t,g,residual,omega2,Tend2
 INTEGER                     :: i, j, k, nsteps, iter

 !! Set up domain
 dx       = (b-a)/REAL(N-1,dpp)
 dt       = CFL*dx
 dt2      = dt**2.d0
 nsteps   = INT(Tend/dt)
 residual = 1.d0
 iter     = 0
 omega2   = 2.d0*OMEGA 
 Tend2    = 2.d0*PI/omega2

 !! Allocate working arrays
 ALLOCATE(x(N+2),u(N+2),rho(N+2),ut(N+2),utt(N+2))
 ALLOCATE(v(N+2),c(N+2),f(N+2),up(N+2),un(N+2),vp(N+2))
 ALLOCATE(v_tmp(N+2),b_vec(N+2),u_exact(N+2))

 !! spatial grid (include ghost points)
 do i = -1,N
  x(i+2) = a + REAL(i,dpp)*dx
  ! u_exact(i+2) = EXP(-SIGMA*x(i+2)**2.d0)
  u_exact(i+2) = COS(0.5d0*PI*x(i+2))
 end do

 !! get initial conditions and physical constants
 CALL initial_condition(x,N+2,v,ut)
 CALL medium_speed(x,N+2,c)
 CALL density(x,N+2,rho)
 CALL spatial_forcing(x,N+2,f)

 !! Before we begin, apply to 0
 b_vec = 0.0_dpp
 CALL evolve_and_project(N,dx,Tend,omega,c,f,rho,x,b_vec)
 vp    = b_vec
 v_tmp = b_vec
 v     = b_vec

 !! Fixed point iteration to find time-harmonic solution
 do while((residual .gt. tol) .AND. (iter .lt. MAXIT))

  CALL evolve_and_project(N,dx,Tend,omega,c,f,rho,x,vp)
  ! v = v - vp + b_vec

  residual = MAXVAL(ABS(vp(2:N+1) - v_tmp(2:N+1)))
  ! residual = MAXVAL(ABS(vp(2:N+1) - u_exact(2:N+1)))

  write(*,'((E24.16))') residual
  
  ! vp    = v
  v_tmp = vp
  iter  = iter + 1

  !! set initial conditions for next problem
  ut = 0.d0
 end do 
 
 ! v = vp
 ! !! Now that we converged, take another step half a period
 ! CALL evolve_and_project(N,dx,Tend,omega,c,f,rho,x,v)

 !! write data to file
 open(2,file='x.txt',status='unknown')
 do i=2,N+1
  write(2,fmt='((E24.16))') x(i)
 end do
 close(2)

 open(3,file='u.txt',status='unknown')
 do i=2,N+1
  write(3,fmt='((E24.16))') vp(i)
  ! write(3,fmt='((E24.16))') 0.5d0*(vp(i) - v(i))
 end do
 close(3)

 !! Deallocate all working arrays
 DEALLOCATE(b_vec)
 DEALLOCATE(v_tmp)
 DEALLOCATE(vp)
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