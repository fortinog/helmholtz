module gmres_module
 
 implicit none
 
 INTEGER,      PARAMETER :: dpp = KIND(1.d0)
 
 contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine performs GMRES on a matrix A with rhs
 !! b (with restarts) using Givens rotations
 !! Inputs:
 !! 	- A       : a matrix
 !!		- b		    : RHS vector
 !! 	- tol     : desired tolerance in error/residual
 !! Outputs:
 !!		- x       : On input the initial guess to the
 !!					      solution. On output, approx. solution to 
 !!					      Ax=b.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine gmres(A,b,tol,x)
  
  implicit none
  
  real(kind=dpp), intent(in)    :: A(:,:)
  real(kind=dpp), intent(in)    :: b(:)
  real(kind=dpp), intent(in)    :: tol
  real(kind=dpp), intent(inout) :: x(:)

  integer                     :: i, j, n, maxiter, ind
  real(kind=dpp)              :: beta, err, denom
  real(kind=dpp)              :: residual, val
  real(kind=dpp), allocatable :: res(:), h(:,:), tmp(:)
  real(kind=dpp), allocatable :: V(:,:), y(:), h_vec(:)
  real(kind=dpp), allocatable :: c_vec(:), s_vec(:)

  n = SIZE(x)
  maxiter = MIN(n,20)

  ALLOCATE(res(n), tmp(n), V(n,maxiter+1))
  ALLOCATE(h(maxiter+1,maxiter), y(maxiter+1))
  ALLOCATE(h_vec(maxiter), c_vec(maxiter), s_vec(maxiter))
  V = 0.d0
  h = 0.d0

  !! start here on a restart
  777 CONTINUE
  res = b - MATMUL(A,x)

  !! normalize the residual
  beta = SQRT(DOT_PRODUCT(res,res))
  res = res/beta
  V(:,1) = res

  !! preset some values
  y = 0.d0
  y(1) = beta
  ind = maxiter

  !! perform Gram-Schmidt on successive Krylov subspaces
  !! and solve a least squares to find an approximate
  !! solution
  do j = 1,maxiter
    
    !! new vector in Krylov subspace
    tmp = MATMUL(A,V(:,j))

    !! find component along previous basis vectors
    do i = 1,j
      h(i,j) = DOT_PRODUCT(tmp,V(:,i))
    end do 

    !! remove linear dependence
    do i = 1,j
      tmp = tmp - h(i,j)*V(:,i)
    end do

    !! normalize and store new basis vector
    h(j+1,j) = SQRT(DOT_PRODUCT(tmp,tmp))
    tmp = tmp/h(j+1,j)
    V(:,j+1) = tmp

    !! apply all previous Givens rotations
    do i = 1,j-1
      val = h(i,j)
      h(i,j) = c_vec(i)*val + s_vec(i)*h(i+1,j)
      h(i+1,j) = -s_vec(i)*val + c_vec(i)*h(i+1,j)
    end do

    !! Givens rotation factors for update
    denom = SQRT(h(j,j)**2.0d0 + h(j+1,j)**2.0d0)
    s_vec(j) = h(j+1,j)/denom
    c_vec(j) = h(j,j)/denom

    !! update RHS vector
    val = y(j)
    y(j) = c_vec(j)*val + s_vec(j)*y(j+1)
    y(j+1) = -s_vec(j)*val + c_vec(j)*y(j+1)

    !! update Hessenberg matrix
    val = h(j,j)
    h(j,j) = c_vec(j)*val + s_vec(j)*h(j+1,j)
    h(j+1,j) = -s_vec(j)*val + c_vec(j)*h(j+1,j)

    !! Exit loop if we get desired error
    if(ABS(y(j+1)) .lt. tol) then 
      ind = j
      exit
    end if

  end do

  !! Backsolve to find the coefficients (ignore last row
  !! of y and H)
  y(ind) = y(ind)/h(ind,ind)
  do i = ind-1,1,-1
    do j = ind,i+1,-1
      y(i) = y(i) - h(i,j)*y(j)
    end do
    y(i) = y(i)/h(i,i)
  end do 

  !! Build our approximate solution
  do j = 1,ind
    x = x + y(j)*V(:,j)
  end do 

  !! Check for a restart if necessary
  if(ind .eq. maxiter) then 
    residual = y(maxiter+1)
    if(ABS(residual) .gt. tol) then 
      GOTO 777
    end if 
  end if

  !! deallocate all working arrays
  DEALLOCATE(s_vec)
  DEALLOCATE(c_vec)
  DEALLOCATE(h_vec)
  DEALLOCATE(y)
  DEALLOCATE(h)
  DEALLOCATE(V)
  DEALLOCATE(tmp)
  DEALLOCATE(res)

 end subroutine gmres

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine performs GMRES on a matrix A with rhs
 !! b (with restarts) using Givens rotations
 !! Inputs:
 !!   - A       : a matrix
 !!   - b       : RHS vector
 !!   - tol     : desired tolerance in error/residual
 !! Outputs:
 !!   - x       : On input the initial guess to the
 !!               solution. On output, approx. solution to 
 !!               Ax=b.
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 subroutine gmres_func(b,tol,x)
  
  implicit none
  
  real(kind=dpp), intent(in)    :: b(:)
  real(kind=dpp), intent(in)    :: tol
  real(kind=dpp), intent(inout) :: x(:)

  integer                     :: i, j, n, maxiter, ind
  real(kind=dpp)              :: beta, err, denom
  real(kind=dpp)              :: residual, val
  real(kind=dpp), allocatable :: res(:), h(:,:), tmp(:)
  real(kind=dpp), allocatable :: V(:,:), y(:), h_vec(:)
  real(kind=dpp), allocatable :: c_vec(:), s_vec(:)

  n = SIZE(x)
  maxiter = MIN(n,20)
  ALLOCATE(res(n), tmp(n), V(n,maxiter+1))
  ALLOCATE(h(maxiter+1,maxiter), y(maxiter+1))
  ALLOCATE(h_vec(maxiter), c_vec(maxiter), s_vec(maxiter))
  V = 0.d0
  h = 0.d0

  !! start here on a restart
  777 CONTINUE
  res = b - MY_MATMUL(x, n)

  !! normalize the residual
  beta = SQRT(DOT_PRODUCT(res,res))
  res = res/beta
  V(:,1) = res

  !! preset some values
  y = 0.d0
  y(1) = beta
  ind = maxiter

  !! perform Gram-Schmidt on successive Krylov subspaces
  !! and solve a least squares to find an approximate
  !! solution
  do j = 1,maxiter
    
    !! new vector in Krylov subspace
    tmp = MY_MATMUL(V(:,j), n)

    !! find component along previous basis vectors
    do i = 1,j
      h(i,j) = DOT_PRODUCT(tmp,V(:,i))
    end do 

    !! remove linear dependence
    do i = 1,j
      tmp = tmp - h(i,j)*V(:,i)
    end do

    !! normalize and store new basis vector
    h(j+1,j) = SQRT(DOT_PRODUCT(tmp,tmp))
    tmp = tmp/h(j+1,j)
    V(:,j+1) = tmp

    !! apply all previous Givens rotations
    do i = 1,j-1
      val = h(i,j)
      h(i,j) = c_vec(i)*val + s_vec(i)*h(i+1,j)
      h(i+1,j) = -s_vec(i)*val + c_vec(i)*h(i+1,j)
    end do

    !! Givens rotation factors for update
    denom = SQRT(h(j,j)**2.0d0 + h(j+1,j)**2.0d0)
    s_vec(j) = h(j+1,j)/denom
    c_vec(j) = h(j,j)/denom

    !! update RHS vector
    val = y(j)
    y(j) = c_vec(j)*val + s_vec(j)*y(j+1)
    y(j+1) = -s_vec(j)*val + c_vec(j)*y(j+1)

    !! update Hessenberg matrix
    val = h(j,j)
    h(j,j) = c_vec(j)*val + s_vec(j)*h(j+1,j)
    h(j+1,j) = -s_vec(j)*val + c_vec(j)*h(j+1,j)

    !! Exit loop if we get desired error
    if(ABS(y(j+1)) .lt. tol) then 
      ind = j
      exit
    end if

  end do

  !! Backsolve to find the coefficients (ignore last row
  !! of y and H)
  y(ind) = y(ind)/h(ind,ind)
  do i = ind-1,1,-1
    do j = ind,i+1,-1
      y(i) = y(i) - h(i,j)*y(j)
    end do
    y(i) = y(i)/h(i,i)
  end do 

  !! Build our approximate solution
  do j = 1,ind
    x = x + y(j)*V(:,j)
  end do 

  !! Check for a restart if necessary
  if(ind .eq. maxiter) then 
    residual = y(maxiter+1)
    if(ABS(residual) .gt. tol) then 
      GOTO 777
    end if 
  end if

  !! deallocate all working arrays
  DEALLOCATE(s_vec)
  DEALLOCATE(c_vec)
  DEALLOCATE(h_vec)
  DEALLOCATE(y)
  DEALLOCATE(h)
  DEALLOCATE(V)
  DEALLOCATE(tmp)
  DEALLOCATE(res)

 end subroutine gmres_func


function MY_MATMUL(c, N)
  implicit none

  real(kind=dpp), intent(inout) :: c(:)
  integer,        intent(in) :: N

  integer        :: i
  real(kind=dpp) :: MY_MATMUL(N)

  MY_MATMUL(1) = -2.0d0*c(1) + c(2)
  do i = 2, N-1
    MY_MATMUL(i) = c(i-1) - 2.0d0*c(i) + c(i+1)
  end do 
  MY_MATMUL(N) = c(N-1) - 2.0d0*c(N)

 end function MY_MATMUL

end module gmres_module