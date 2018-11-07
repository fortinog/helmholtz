program test_qrdelete

 use anderson_acceleration
 implicit none 

 INTEGER, PARAMETER :: N = 10
 INTEGER        :: i, j, m, maxiter
 REAL(KIND=dpp) :: tol, x(N)
 
 m = 10
 tol = 1d-15
 maxiter = 20

 x = 1.d0

 ! REAL(KIND=dpp) :: A(N,N), b(N)

 ! A(1,1) = 2.d0
 ! A(1,2) = -1.d0
 ! do i =2,N-1
 ! 	A(i,i-1) = -1.d0
 ! 	A(i,i) = 2.d0
 ! 	A(i,i+1) = -1.d0
 ! end do

 ! A(N,N-1) = -1.d0
 ! A(N,N) = 2.d0

 ! A = -1.0d0*A
 ! b = 0.d0
 ! b(1) = 1.d0

 CALL anderson(m,tol,maxiter,N,x)

 write(*,*) x
!  open(3,file='Q.txt',status='unknown')
!  open(4,file='R.txt',status='unknown')
!  do j =1,3
!   do i = 1,3
!    read(3,*) Q(i,j)
!    read(4,*) R(i,j)
!   end do 
! end do 
! close(3)
! close(4)

!  CALL QRdelete(3,3,Q,R)
!  write(*,*) Q
!  write(*,*) "Here is a break"
!  write(*,*) R

end program test_qrdelete