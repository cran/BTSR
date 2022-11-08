module main_mod
  implicit none
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  !      This module contains the user defined variables and parameters.
  !
  !      This new types are used to pass extra information to subroutines.
  !
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! double precision: 
  integer, parameter :: dp = kind(1.d0)

  real(dp), parameter :: pi = acos(-1.d0)

  !-------------------------------------------------------------
  !
  !  Arguments for the Link function
  !
  !-------------------------------------------------------------
  type argslink
     ! here we can include link parameters if necessary
     integer :: link       ! the type of link  
     real(dp) :: lower     ! x.lower
     real(dp) :: upper     ! x.upper
     real(dp) :: a = 1.d0  ! a constant for the link
  end type argslink

  !-------------------------------------------------------------
  !
  !  Arguments for the distribution function (if necessary)
  !
  !-------------------------------------------------------------
  type argsDist
     real(dp) :: lower   ! y.lower
     real(dp) :: upper   ! y.upper
     real(dp) :: arg1    ! extra arguments
  end type argsDist


  !-------------------------------------------------------------
  !
  !  Bounds for parameters (for parameter's transformation)
  !
  !-------------------------------------------------------------  
  type par_bounds
     integer, allocatable :: nbd(:)  ! type of bound: 0, 1, 2, 3, 4
     real(dp), allocatable :: lower(:) 
     real(dp), allocatable :: upper(:)
  end type par_bounds

  !-------------------------------------------------------------
  !
  !  Vectors of parameters 
  !
  !-------------------------------------------------------------
  type vec_parameter
     integer :: length = 0
     integer :: fit = 0                ! number of non-fixed values
     integer, allocatable :: lags(:)   ! non-fixed lags
     integer, allocatable :: flags(:)  ! fixed lags
     real(dp), allocatable :: fv(:)    ! fixed values
     real(dp), allocatable :: par(:)   ! vector of parameters
  end type vec_parameter

  !-------------------------------------------------------------
  !
  !  mu/sigma and corresponding regressors
  !
  !-------------------------------------------------------------
  type conditional_ts
     integer :: n      ! size of ut
     integer :: nreg   ! number of regressors
     integer :: xregar ! Indicates if xreg must be included in the AR equation
     real(dp), allocatable :: xstart(:) ! starting values for xreg
     real(dp), allocatable :: ut(:)     ! conditional mean/variance
     real(dp), allocatable :: eta(:)    ! transformed mean/variance
     real(dp), allocatable :: xreg(:,:) ! regressors
     real(dp), allocatable :: orbit(:)  ! T**t(u0)
  end type conditional_ts

  !-------------------------------------------------------------
  !
  ! Partial derivatives deta/dgamma
  !
  !-------------------------------------------------------------
  type deta_d
     real(dp), allocatable :: dalpha(:,:)  ! constant
     real(dp), allocatable :: dbeta(:,:)   ! covariates 
     real(dp), allocatable :: dphi(:,:)      ! ar   
     real(dp), allocatable :: dtheta(:,:)  ! ma   
     real(dp), allocatable :: dd(:,:)      ! d
     real(dp), allocatable :: dthetaT(:,:) ! chaotic     
  end type deta_d

  !-------------------------------------------------------------
  !
  ! Matrices for the score vectors U(gamma)
  !
  !-------------------------------------------------------------
  type score
     real(dp), allocatable :: Unu(:)     ! precision
     real(dp), allocatable :: Ualpha(:)  ! constant
     real(dp), allocatable :: Ubeta(:)   ! covariates 
     real(dp), allocatable :: Uphi(:)      ! ar   
     real(dp), allocatable :: Utheta(:)  ! ma   
     real(dp), allocatable :: Ud(:)      ! d
     real(dp), allocatable :: UthetaT(:) ! chaotic     
  end type score

  !-------------------------------------------------------------
  !
  ! Generic Vector
  !
  !-------------------------------------------------------------
  type vetor
     real(dp), allocatable :: par(:)
  end type vetor

  !-------------------------------------------------------------
  !
  !  Score vector and related matrices
  !
  !-------------------------------------------------------------
  type argsSI
     !--------------------------------------------
     ! deta(1,1) deta1_drho;     deta(2,1) deta2_drho     
     ! deta(1,2) deta1_dlambda;  deta(2,2) deta2_dlambda
     !--------------------------------------------
     type(deta_d) :: deta(2,2)  
     !----------------------
     ! Score Vector
     !----------------------
     type(score) :: U(2)
     !----------------------
     ! T matrix
     !----------------------     
     type(vetor) :: T(2)
     !----------------------
     ! h vector
     !----------------------     
     type(vetor) :: h(2)
     !----------------------
     ! E matrices
     !----------------------
     real(dp), allocatable :: E(:,:)
  end type argsSI

  !-------------------------------------------------------------
  !
  !  Model specification
  !
  !-------------------------------------------------------------
  type argsModel
     integer :: m = 0       ! starting point for log-likelihood
     integer :: n           ! sample size
     integer :: fixnu = 1   ! type of model: 1 = constant precision, 0 = time varying precision
     integer :: npar(2) = 0 ! number of non-fixed parameters
     integer :: inf(2) = 0  ! number of terms in the infinite sum expansion
     integer :: llk = 0     ! 1 = calculate loglikelihood 
     integer :: sco = 0     ! 1 = calculate score vector
     integer :: info = 0    ! 1 = calculate information matrix
     type(argslink), allocatable, dimension(:) :: argsL  ! arguments to be passed to the link function
     type(argsSI) :: SI     ! score vector and information matrix
     real(dp) :: r20        ! starting value: squared error term

     ! parameters of the model
     type(vec_parameter) :: nu     
     type(vec_parameter) :: alpha(2) 
     type(vec_parameter) :: beta(2) 
     type(vec_parameter) :: ar(2)  
     type(vec_parameter) :: ma(2)
     type(vec_parameter) :: d(2)
     type(vec_parameter) :: thetaT
     type(vec_parameter) :: u0

     ! parameter bounds for optimization
     type(par_bounds) :: bounds

     ! chaotic transformation
     integer :: map = 0

     real(dp), allocatable, dimension(:) :: y      ! original time series
     real(dp), allocatable, dimension(:) :: gy     ! transformed time series
     type(conditional_ts) :: cts(2)                ! conditional series and related regressors
     real(dp), allocatable, dimension(:) :: error  ! error term
     integer :: escale = 1  ! 1: error = g(y) - g(mu); 0: error = y - mu
     real(dp) :: ystart     ! starting value for y

     ! fixed arguments related to the distribution
     type(argsDist) :: argsD
  end type argsModel

  ! Interface - to use the same call for any allocation subroutine
  interface safe_allocate
     module procedure safe_allocateR1  ! vector 1:n
     module procedure safe_allocateR1n ! vector m:n
     module procedure safe_allocateR2  ! matrix 1:n1, 1:n2
     module procedure safe_allocateR2n ! matrix m1:n1, m2:n2
     module procedure safe_allocateI1  ! vector of integers
  end interface safe_allocate

contains

  !****************************************************************
  !
  !  Subroutines used to safe allocate vectors and matrices
  !
  !  Just to make sure that the program will not try to allocate
  !  a vector/matrix that is already allocated
  !
  !****************************************************************
  subroutine safe_allocateR1(x,n)
    implicit none
    integer, intent(in) :: n
    real(dp), allocatable, intent(inout) ::  x(:)
    if(allocated(x)) deallocate(x)
    allocate(x(n))
    return
  end subroutine safe_allocateR1

  subroutine safe_allocateR1n(x,n1,n2)
    implicit none
    integer, intent(in) :: n1, n2
    real(dp), allocatable, intent(inout) ::  x(:)
    if(allocated(x)) deallocate(x)
    allocate(x(n1:n2))
    return
  end subroutine safe_allocateR1n

  subroutine safe_allocateR2(x,m,n)
    implicit none
    integer, intent(in) :: m,n
    real(dp), allocatable, intent(inout) ::  x(:,:)
    if(allocated(x)) deallocate(x)
    allocate(x(m,n))
    return
  end subroutine safe_allocateR2

  subroutine safe_allocateR2n(x,m1,m2,n1,n2)
    implicit none
    integer, intent(in) :: m1, m2, n1, n2
    real(dp), allocatable, intent(inout) ::  x(:,:)
    if(allocated(x)) deallocate(x)
    allocate(x(m1:m2,n1:n2))
    return
  end subroutine safe_allocateR2n

  subroutine safe_allocateI1(x,n)
    implicit none
    integer, intent(in) :: n
    integer, allocatable, intent(inout) ::  x(:)
    if(allocated(x)) deallocate(x)
    allocate(x(n))
    return
  end subroutine safe_allocateI1



  !------------------------------------------------------------------------------------------------------
  !
  !               Subroutines for parameter transformation
  !
  !
  ! Based on the Matlab subroutine fminsearchbnd
  ! Reference:
  ! https://www.mathworks.com/matlabcentral/fileexchange/8277-
  !                 fminsearchbnd-fminsearchcon?focused=5216898&tab=function
  !
  !------------------------------------------------------------------------------------------------------
  subroutine set_bounds(bounds, lower, upper, nbd, npar)
    !
    ! allocates the vectors
    !    lower bound, upper bound and type of bounds (nbd)
    ! in the bounds object and save the values passed by the user.
    !
    implicit none
    type(par_bounds), intent(inout) :: bounds
    integer, intent(in) :: npar
    integer, intent(in) :: nbd(npar)
    real(dp), intent(in) :: lower(npar), upper(npar)

    call safe_allocate(bounds%lower, npar)
    call safe_allocate(bounds%upper, npar)
    call safe_allocate(bounds%nbd, npar)

    bounds%lower = lower
    bounds%upper = upper
    bounds%nbd = nbd   
    return    
  end subroutine set_bounds

  function xtransform(npars, x, bounds) result(xtrans)
    !------------------------------------------------------------------
    ! converts unconstrained variables into their original domains
    ! bounds type:
    !  0 = none
    !  1 = lower
    !  2 = both
    !  3 = upper
    !  4 = constant
    !------------------------------------------------------------------
    implicit none
    type(par_bounds), intent(in) :: bounds
    integer, intent(in) :: npars
    real(dp), intent(in) :: x(npars) 
    real(dp) :: xtrans(npars)
    integer ::  i

    do i = 1, npars
       select case(bounds%nbd(i))
       case(0)  ! none
          !% unconstrained variable.
          xtrans(i) = x(i)
       case(1)  ! lower
          !% lower bound only
          xtrans(i) = bounds%lower(i) + x(i)**2
       case(2)  ! both
          !% lower and upper bounds
          xtrans(i) = (sin(x(i))+1.d0)/2.d0
          xtrans(i) = xtrans(i)*(bounds%upper(i) - bounds%lower(i)) + bounds%lower(i)
          !% just in case of any floating point problems
          xtrans(i) = max(bounds%lower(i),min(bounds%upper(i),xtrans(i)))
       case(3)  ! upper
          ! upper bound only
          xtrans(i) = bounds%upper(i) - x(i)**2
       case(4)  ! constant
          !% fixed variable, bounds are equal, set it at either bound
          xtrans(i) = bounds%lower(i)
       end select
    end do
    return
  end function xtransform

  function xtransformstart(npars, x, bounds) result(xtrans)
    !--------------------------------------------------------
    !% transform starting values into their unconstrained
    !% surrogates. check for infeasible starting guesses.
    !--------------------------------------------------------
    implicit none
    type(par_bounds), intent(in) :: bounds
    integer, intent(in) :: npars
    real(dp), intent(in) :: x(npars) 
    real(dp) :: xtrans(npars)
    integer ::  i

    do i = 1, npars
       select case(bounds%nbd(i))
       case(0) ! none
          !% unconstrained variable.
          xtrans(i) = x(i)
       case(1) ! lower
          !% lower bound only
          if(x(i) <= bounds%lower(i)) then
             !% infeasible starting value. use bound.
             xtrans(i) = 0.d0
          else
             xtrans(i) = sqrt(x(i) - bounds%lower(i))
          end if
       case(2) ! both
          !% lower and upper bounds
          if(x(i) <= bounds%lower(i)) then
             !% infeasible starting value
             xtrans(i) = -pi/2.d0
          elseif(x(i) >= bounds%upper(i)) then
             !% infeasible starting value
             xtrans(i) = pi/2.d0
          else
             xtrans(i) = 2*(x(i) - bounds%lower(i))/(bounds%upper(i)-bounds%lower(i)) - 1.d0
             !% shift by 2*pi to avoid problems at zero in fminsearch
             !% otherwise, the initial simplex is vanishingly small
             xtrans(i) = 2.d0*pi+asin(max(-1.d0,min(1.d0,xtrans(i))))
          end if
       case(3) ! upper
          ! upper bound only
          if(x(i) >= bounds%upper(i)) then
             !% infeasible starting value. use bound.
             xtrans(i) = 0.d0
          else
             xtrans(i) = sqrt(bounds%upper(i) - x(i))
          end if
       case(4) ! constant
          !% fixed variable, bounds are equal, set it at either bound
          xtrans(i) = bounds%lower(i)
       end select
    end do
    return      
  end function xtransformstart


  subroutine transform_par(par, npar, bounds, inverse)
    !---------------------------------------------------------------------------
    ! transforms the y = par (bounded) in x = f(y) (unbounded)
    ! if inverse = .true., transforms x (unbounded) in y = f^{-1}(x) (bounded)
    !---------------------------------------------------------------------------
    implicit none
    integer, intent(in) :: npar
    real(dp), intent(inout) :: par(npar)
    type(par_bounds), intent(in) :: bounds
    logical, intent(in) :: inverse

    if(sum(bounds%nbd) == 0) return    

    if(inverse .eqv. .true.) goto 100
    ! bounded -> unbounded
    par = xtransformstart(npar, par, bounds)
    return

100 continue
    ! unbounded -> bounded
    par = xtransform(npar, par, bounds)
    return
  end subroutine transform_par

end module main_mod
