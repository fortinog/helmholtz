module conjugate_gradient
 implicit none
 INTEGER,      PARAMETER :: dpp = KIND(1.d0)
 
 contains

 !! This subroutine performs the standard Conjugate
 !! Gradient method to find a solution to the system
 !! Ax = b.
 !! Inputs:
 !! 	- A       : a pos. def. matrix A
 !!		- b		  : RHS vector
 !!		- maxiter : maximum number of iterations
 !! 	- tol     : tolerance in error/residual
 !! Outputs:
 !!		- x       : On input, initial guess to 
 !!					solution. On ouput, approx. 
 !!					solution to Ax=b.
 subroutine CG(A,b,maxiter,tol,x)
  
  implicit none
  
  real(kind=dpp), intent(in)    :: A(:,:)
  real(kind=dpp), intent(in)    :: b(:)
  integer,        intent(in)    :: maxiter
  real(kind=dpp), intent(in)    :: tol
  real(kind=dpp), intent(inout) :: x(:)

  integer                     :: i, n
  real(kind=dpp)              :: alpha, beta, err, tmp
  real(kind=dpp), allocatable :: res(:), dir(:)

  !! allocate and compute 
  !! initial vectors
  n = SIZE(x)
  allocate(res(1:n), dir(1:n))
  dir = b - MATMUL(A,x)
  res = dir
  err = DOT_PRODUCT(res,res)
  i = 0

  !! travel along conjugate directions
  do while((i .lt. maxiter) .and. (tol .lt. ABS(err)))
   alpha = err/(DOT_PRODUCT(dir,MATMUL(A,dir)))
   x = x + alpha*dir
   res = res - alpha*MATMUL(A,dir)
   tmp = DOT_PRODUCT(res,res)
   beta = tmp/err
   err = tmp
   dir = res + beta*dir
   i = i + 1
  end do

  DEALLOCATE(res)
  DEALLOCATE(dir)

 end subroutine CG


end module conjugate_gradient