module anderson_acceleration

 implicit none
 INTEGER, PARAMETER :: dpp   = KIND(1.d0)

contains

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in an initial guess and attempts
 !! to find a stationary point defined in the function 
 !! rhs_eval via Anderson acceleration.
 !! Inputs:
 !!     - m       : number of residual vectors to track
 !!     - tol     : desired level of tolerance in residual
 !!     - maxiter : maximum number of allowable iterations
 !!     - N       : number of spatial nodes
 !! Outputs:
 !! 	- x       : approximate solution to fixed point
 !!  				iteration
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine anderson(m,tol,maxiter,N,x)

  implicit none

  INTEGER,        INTENT(IN)    :: m
  REAL(KIND=dpp), INTENT(IN)    :: tol 
  INTEGER,        INTENT(IN)    :: maxiter 
  INTEGER,        INTENT(IN)    :: N 
  REAL(KIND=dpp), INTENT(INOUT) :: x(:)

  REAL(KIND=dpp), ALLOCATABLE :: Fmat(:,:), Gmat(:,:)
  REAL(KIND=dpp), ALLOCATABLE :: Q(:,:), R(:,:), F(:), G(:)
  REAL(KIND=dpp), ALLOCATABLE :: F_old(:), G_old(:)
  REAL(KIND=dpp), ALLOCATABLE :: df(:), dg(:), gamma_vec(:)
  REAL(KIND=dpp)              :: residual
  INTEGER                     :: i, j, k, num_res

  !! Allocate working arrays
  ALLOCATE(Fmat(N,m),Gmat(N,m),Q(N,m),R(m,m),F(N),G(N))
  ALLOCATE(F_old(N),G_old(N),df(N),dg(N),gamma_vec(m))

  Q       = 0.d0
  R       = 0.d0
  num_res = 0   !! number of residuals stored
  F_old = 0.d0
  G_old = 0.d0

  do k=0,maxiter
   CALL rhs_eval(x,N,G)
   F = G - x
   residual = SQRT(DOT_PRODUCT(F,F))

   !! We have achieved the sought error, exit
   if(residual .le. tol) then 
    exit
   end if

   if(m .eq. 0) then 
    !! No acceleration, do a usual FPI
    CALL rhs_eval(x,N,G)
    x = G
    cycle
   else 
    !! Perform an Anderson acceleration
    df = F - F_old
    if(num_res .lt. m) then 
     !! Append to G
     Gmat(:,num_res+1) = G - G_old
    else
     !! Delete first column of G, and append new column
     do j = 1,m-1
      Gmat(:,j) = Gmat(:,j+1)
     end do 
     Gmat(:,m) = G - G_old
    end if
    num_res = num_res + 1 
    !! this if statement not where it needs to be
   end if
   G_old = G 
   F_old = F

   if(num_res .eq. 1) then 
    !! form initial QR decomposition
    R(1,1) = SQRT(DOT_PRODUCT(df,df))
    Q(:,1) = df/R(1,1)
   else 
    !! Update QR decomposition
    if(num_res .gt. m) then 
     !! Update Q,R. Here we don't delete, we simply
     !! shift columns up one and allow the last column to
     !! be deleted.
     CALL QRdelete(m,N,Q,R) 
     num_res = num_res - 1
    end if

    !! Now we update the QR decomposition to incorporate
    !! the new column
    do j = 1, num_res-1
     R(j,num_res) = DOT_PRODUCT(Q(:,j),df)
     df = df - R(j,num_res)*Q(:,j)
    end do
    R(num_res,num_res) = SQRT(DOT_PRODUCT(df,df))
    Q(:,num_res) = df/R(num_res,num_res)
   end if

   !! Now we need to solve the least-squares problem
   gamma_vec = 0.d0
   do j=1,num_res
    gamma_vec = gamma_vec + DOT_PRODUCT(Q(:,j),F)
   end do 

   !! Do a backsolve of R to get gamma_vec
   gamma_vec(num_res) = gamma_vec(num_res)/R(num_res,num_res)
   do i=num_res-1,1,-1
    do j=num_res,i+1,-1
     gamma_vec(i) = gamma_vec(i) - R(i,j)*gamma_vec(j)
    end do 
    gamma_vec(i) = gamma_vec(i)/R(i,i)
   end do 

   !! update approximate solution
   x = G - MATMUL(Gmat,gamma_vec)
  end do 

  !! Deallocate working arrays
  DEALLOCATE(gamma_vec)
  DEALLOCATE(dg)
  DEALLOCATE(df)
  DEALLOCATE(G_old)
  DEALLOCATE(F_old)
  DEALLOCATE(G)
  DEALLOCATE(F)
  DEALLOCATE(R)
  DEALLOCATE(Q)
  DEALLOCATE(Gmat)
  DEALLOCATE(Fmat)

 end subroutine anderson


 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in spatial nodes and returns 
 !! a vector of function values for which we wish to find
 !! a stationary point of.
 !! Inputs:
 !! 	- x       : spatial nodes
 !!     - N       : number of spatial nodes
 !! Outputs:
 !!		- g       : function values
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine rhs_eval(x,N,out)

  implicit none
  
  REAL(KIND=dpp), INTENT(IN)    :: x(:)
  INTEGER,        INTENT(IN)    :: N
  REAL(KIND=dpp), INTENT(INOUT) :: out(:)

  INTEGER :: i

  ! do i = 1,N 
  !  out(i) = x(i)**2.d0
  ! end do 

  !! Solve Poissons eq.
  out(1) = -2.d0*x(1) + x(2)
  do i = 2,N-1 
   out(i) = x(i-1) - 2.d0*x(i) + x(i+1)
  end do
  out(N) = x(N-1) - 2.d0*x(N)

  out = x + out
  !! subtract b
  out(1) = out(1) - 1.d0
 end subroutine rhs_eval

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !! This subroutine takes in Q,R from the QR decomposition
 !! of the Anderson matrix and deletes the first column 
 !! of A and appropriately updates the factorization.
 !! Inputs:
 !! 	- x       : spatial nodes
 !!     - N       : number of spatial nodes
 !! Outputs:
 !!		- g       : function values
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine QRdelete(m,N,Q,R)
  implicit none 

  INTEGER,        INTENT(IN)    :: m
  INTEGER,        INTENT(IN)    :: N
  REAL(KIND=dpp), INTENT(INOUT) :: Q(:,:)
  REAL(KIND=dpp), INTENT(INOUT) :: R(:,:)

  INTEGER        :: i, j
  REAL(KIND=dpp) :: temp, s, c

  do i = 1,m-1
   temp       = SQRT(R(i,i+1)**2.d0 + R(i+1,i+1)**2.d0)
   c          = R(i,i+1)/temp
   s          = R(i+1,i+1)/temp
   R(i,i+1)   = temp
   R(i+1,i+1) = 0.d0
   if(i .lt. m-1) then 
    do j = i+2,m 
     temp = c*R(i,j) + s*R(i+1,j)
     R(i+1,j) = -s*R(i,j) + c*R(i+1,j)
     R(i,j) = temp
    end do 
   end if

   do j = 1,N
    temp = c*Q(j,i) + s*Q(j,i+1)
    Q(j,i+1) = -s*Q(j,i) + c*Q(j,i+1)
    Q(j,i) = temp
   end do 
  end do

  !! Now we remove the columns of interest
  do i =1,m-1
   R(1:m-1,i) = R(1:m-1,i+1)
  end do  

 end subroutine QRdelete

end module anderson_acceleration