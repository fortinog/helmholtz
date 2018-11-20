program helmholtz_wave_fd
 
 use helmholtz_parameters
 implicit none
 
 INTEGER,        PARAMETER   :: num_freq = 5
 ! INTEGER,        PARAMETER   :: num_freq = 1
 REAL(KIND=dpp), ALLOCATABLE :: x(:),u(:),rho(:),ut(:)
 REAL(KIND=dpp), ALLOCATABLE :: utt(:),v(:),c(:),f(:),up(:)
 REAL(KIND=dpp), ALLOCATABLE :: un(:),vp(:),b_vec(:),v_tmp(:)
 REAL(KIND=dpp), ALLOCATABLE :: u_exact1(:), u_exact2(:), u_exact3(:), u_exact4(:),u_exact5(:)
 REAL(KIND=dpp), ALLOCATABLE :: u1(:), u2(:), u3(:), u4(:), u5(:)
 REAL(KIND=dpp), ALLOCATABLE :: w1(:), w2(:), w3(:), w4(:), w5(:)

 REAL(KIND=dpp)     		 :: dx,dt,dt2,t,g,residual, coeff2, coeff3, coeff4
 REAL(KIND=dpp)              :: freq(num_freq), Tend2, Tend3, coeff, Tend4, Tend5

 INTEGER                     :: i, j, k, nsteps, iter

 !! Set up domain
 dx       = (b-a)/REAL(N-1,dpp)
 dt       = CFL*dx
 dt2      = dt**2.d0
 nsteps   = INT(Tend/dt)
 residual = 1.d0
 iter     = 0
 Tend2    = Tend/2.d0
 Tend3    = Tend/4.d0
 Tend4    = Tend/8.d0
 Tend5    = Tend/15.d0

 !! Allocate working arrays
 ALLOCATE(x(N+2),u(N+2),rho(N+2),ut(N+2),utt(N+2))
 ALLOCATE(v(N+2),c(N+2),f(N+2),up(N+2),un(N+2),vp(N+2))
 ALLOCATE(v_tmp(N+2),b_vec(N+2),u_exact1(N+2),u_exact2(N+2),u_exact3(N+2), u_exact4(N+2),u_exact5(N+2))
 ALLOCATE(u1(N+2),u2(N+2),u3(N+2),u4(N+2),u5(N+2))
 ALLOCATE(w1(N+2),w2(N+2),w3(N+2),w4(N+2),w5(N+2))

 do i =0,num_freq-1
  freq(i+1) = (2.d0**REAL(i,dpp))*OMEGA
 end do

 coeff   = (PI**2.d0/4.d0 - OMEGA**2.d0)
 coeff   = coeff/(PI**2.d0/4.d0 - 4.d0*OMEGA**2.d0) 
 coeff2  = (PI**2.d0/4.d0 - OMEGA**2.d0)/(PI**2.d0/4.d0 - 16.d0*OMEGA**2.d0)
 coeff3  = (PI**2.d0/4.d0 - OMEGA**2.d0)/(PI**2.d0/4.d0 - 64.d0*OMEGA**2.d0)
 coeff3  = (PI**2.d0/4.d0 - OMEGA**2.d0)/(PI**2.d0/4.d0 - 256.d0*OMEGA**2.d0)

 !! spatial grid (include ghost points)
 do i = -1,N
  x(i+2) = a + REAL(i,dpp)*dx
  ! u_exact(i+2) = EXP(-SIGMA*x(i+2)**2.d0)
  u_exact1(i+2) = COS(0.5d0*PI*x(i+2))
  u_exact2(i+2) = coeff*u_exact1(i+2)
  u_exact3(i+2) = coeff2*u_exact1(i+2)
  u_exact4(i+2) = coeff3*u_exact1(i+2)
  u_exact5(i+2) = coeff4*u_exact1(i+2)

  ! u_exact1(i+2) = x(i+2)**2.d0 - 1.d0
  ! u_exact2(i+2) = 0.25d0*(x(i+2)**2.d0 - 1.d0) + 3.d0/(8.d0*OMEGA**2.d0)
 end do

 !! get initial conditions and physical constants
 CALL initial_condition(x,N+2,v,ut)
 CALL medium_speed(x,N+2,c)
 CALL density(x,N+2,rho)
 CALL spatial_forcing(x,N+2,f)

 !! Before we begin, apply to 0
 b_vec = 0.0_dpp
 CALL evolve_and_project(N,dx,Tend,freq,num_freq,c,f,rho,x,b_vec)
 vp    = b_vec
 v_tmp = b_vec
 v     = b_vec

 !! Fixed point iteration to find time-harmonic solution
 do while((residual .gt. tol) .AND. (iter .lt. MAXIT))

  CALL evolve_and_project(N,dx,Tend,freq,num_freq,c,f,rho,x,vp)

  ! residual = MAXVAL(ABS(vp(2:N+1) - v_tmp(2:N+1)))
  residual   = NORM2(vp(2:N+1) - v_tmp(2:N+1))

  write(*,'((E24.16))') residual
  
  ! vp    = v
  v_tmp = vp
  iter  = iter + 1

  !! set initial conditions for next problem
  ut = 0.d0
 end do 
 
 if(num_freq .eq. 1) then 
 	write(*,*) MAXVAL(ABS(vp(2:N+1) - u_exact1(2:N+1)))
 	stop 888
 end if 

 b_vec = vp 
 w1 = vp
 w2 = vp 
 w3 = vp 
 w4 = vp
 w5 = vp
 !! Now that we've converged step forward in time 
 !! one full period of largest frequency
 CALL wave_solve(N,num_freq,dx,Tend2,c,f,rho,x,w2)
 CALL wave_solve(N,num_freq,dx,Tend3,c,f,rho,x,w3)
 CALL wave_solve(N,num_freq,dx,Tend4,c,f,rho,x,w4)
 CALL wave_solve(N,num_freq,dx,Tend5,c,f,rho,x,w5)

 write(*,*) "We have evolved half a period."
 !! Calculate highest frequency solution
 u5 = 0.5d0*w5 + 0.25d0*w4 + 0.301776695296637d0*w3 + 0.293469883127822d0*w2 -0.345246578424459d0*w1
 u4 = -0.5d0*w5 + 0.25d0*w4 - 0.051776695296637d0*w3 + 0.008306812168815d0*w2 + 0.293469883127822d0*w1
 u3 =  0.301776695296637d0*w1 - 0.051776695296637d0*w2 + 0.25d0*w3 - 0.5d0*w4
 u2 =  0.25*w1 + 0.25d0*w2 - 0.5d0*w3
 u1 = 0.5d0*w1 - 0.5d0*w2
 ! stop 888
 ! u3 = 0.25*(vp + v_tmp + 2.d0*b_vec)
 ! u2 = u3 - b_vec
 ! u1 = vp - u2 - u3

 write(*,*) "Error in first solution:", MAXVAL(ABS(u1(2:N+1) - u_exact1(2:N+1)))
 write(*,*) "Error in second solution:", MAXVAL(ABS(u2(2:N+1) - u_exact2(2:N+1)))
 write(*,*) "Error in third solution:", MAXVAL(ABS(u3(2:N+1) - u_exact3(2:N+1)))
 write(*,*) "Error in fourth solution:", MAXVAL(ABS(u4(2:N+1) - u_exact4(2:N+1)))
 write(*,*) "Error in highest frequency solution:", MAXVAL(ABS(u5(2:N+1) - u_exact5(2:N+1)))
 stop 777


 write(*,*) "Error in first solution:", MAXVAL(ABS(0.5d0*(vp(2:N+1) - v_tmp(2:N+1)) - u_exact1(2:N+1)))
 ! write(*,*) "Error in second solution:", MAXVAL(ABS(0.5d0*(vp(2:N+1) + v_tmp(2:N+1)) - u_exact2(2:N+1)))

 !! write data to file
 open(2,file='x.txt',status='unknown')
 do i=2,N+1
  write(2,fmt='((E24.16))') x(i)
 end do
 close(2)

 open(3,file='u.txt',status='unknown')
 do i=2,N+1
  ! write(3,fmt='((E24.16))') vp(i)
  write(3,fmt='((E24.16))') 0.5d0*(vp(i) + v_tmp(i))
 end do
 close(3)

 !! Deallocate all working arrays
 DEALLOCATE(u3)
 DEALLOCATE(u2)
 DEALLOCATE(u1) 
 DEALLOCATE(u_exact3)
 DEALLOCATE(u_exact2)
 DEALLOCATE(u_exact1)
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