!!gfortran -g -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all -Wall -c gmres_module.f90
!! gfortran -g -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all -Wall -c test_GMRES.f90
!! gfortran -o gmres.x gmres_module.o test_GMRES.o -g -ffpe-trap=invalid,zero,overflow -fbacktrace -fcheck=all -Wall

program test_GMRES
 
 use gmres_module
 
 implicit none

 INTEGER, PARAMETER :: N = 30, maxiter = 10
 real(kind=dpp), PARAMETER :: tol = 1.0d-13
 real(kind=dpp) :: A(N,N), b(N), x(N), y(N)
 integer :: i,j

 !! set right hand side
 b = 0.d0
 b(1) = 1.d0
 A = 0.d0
 A(1,1) = 2.d0
 A(1,2) = -1.d0
 x = b
 y = b

 ! x = 0.d0
 ! y = x

 do i =2,N-1
 	A(i,i-1) = -1.d0
 	A(i,i) = 2.d0
 	A(i,i+1) = -1.d0
 end do

 A(N,N-1) = -1.d0
 A(N,N) = 2.d0

 A = -1.0d0*A

 b = MATMUL(A,x)

 ! y = 0.d0
 ! write(*,*) b 
 ! stop 829
 CALL gmres(A,b,tol,y)
 ! CALL gmres_func(b,tol,x)


 !! post output
 do i=1,N
 	! write(*,fmt='(2(E24.16))') y(i), x(i)
 	write(*,fmt='(1(E24.16))') y(i)

 end do

end program test_GMRES