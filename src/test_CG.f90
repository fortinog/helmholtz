program test_CG
 
 use conjugate_gradient
 
 implicit none

 INTEGER, PARAMETER :: N = 10, maxiter = 10
 real(kind=dpp), PARAMETER :: tol = 1.0d-13
 real(kind=dpp) :: A(N,N), b(N), x(N)
 integer :: i,j

 !! set right hand side
 b = 0.d0
 b(1) = 1.d0
 A = 0.d0
 A(1,1) = 2.d0
 A(1,2) = -1.d0
 x = 0.5d0

 do i =2,N-1
 	A(i,i-1) = -1.d0
 	A(i,i) = 2.d0
 	A(i,i+1) = -1.d0
 end do

 A(N,N-1) = -1.d0
 A(N,N) = 2.d0

 CALL CG(A,b,maxiter,tol,x)

 !! post output
 do i=1,N
 	write(*,fmt='((E24.16))') x(i)  
 end do


end program test_CG