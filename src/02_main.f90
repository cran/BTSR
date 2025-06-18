module main_mod
  use lib_utils ! parameters and special functions
  use distrib   ! distribution related subrotuines
  implicit none
  !-------------------------------------------------------------------------------------------
  !
  !      This module contains the user defined variables and parameters.
  !      This new types are used to pass extra information to subroutines.
  !
  !-------------------------------------------------------------------------------------------
  ! July 2023
  !  - major revision: made the necessary changes so that the parameter nu can vary on time.
  !
  ! Last revision: February 2025
  !  - converted some loops to vectorized operations
  !  - merged subroutines to allocate vectors and matrices
  !  - moved argsDist to distrib module
  !-------------------------------------------------------------------------------------------

  !-------------------------------------------------------------
  !
  !  Arguments for the Link function
  !
  !-------------------------------------------------------------
  type argslink
     ! here we can inclUde link parameters if necessary
     ! July, 2023:
     !  - changed "a" to "ctt"
     !
     ! Last revision: February, 2025
     !   - renamed ctt to par
     !   - par now has size 2 (a and b)
     integer  :: link = 0              ! the type of link
     real(dp) :: lower = 0.d0          ! x.lower
     real(dp) :: upper = 1.d0          ! x.upper
     real(dp) :: par(2) = [1.d0, 1.d0] ! link parameters
     logical  :: update = .true.       ! to compare two links
  end type argslink

  !-------------------------------------------------------------
  !
  !  Bounds for parameters (for parameter's transformation)
  !
  !-------------------------------------------------------------
  type par_bounds
     ! Last revision: July, 2023
     integer,  allocatable :: nbd(:)  ! type of bound: 0, 1, 2, 3, 4
     real(dp), allocatable :: lower(:)
     real(dp), allocatable :: upper(:)
  end type par_bounds

  !-------------------------------------------------------------
  !
  !  Vectors of parameters
  !
  !-------------------------------------------------------------
  type vec_parameter
     ! Last revision: February, 2025
     ! - renamed par to cf
     !
     integer :: length = 0
     integer :: fit = 0                 ! number of non-fixed values
     integer,  allocatable :: lags(:)   ! non-fixed lags
     integer,  allocatable :: flags(:)  ! fixed lags
     real(dp), allocatable :: fv(:)     ! fixed values
     real(dp), allocatable :: cf(:)     ! vector of parameters
  end type vec_parameter

  !-------------------------------------------------------------
  !
  ! Partial derivatives deta / dgamma
  !
  !-------------------------------------------------------------
  type deta_d
     ! Last revision: July, 2023
     real(dp), allocatable :: dalpha(:, :)  ! constant
     real(dp), allocatable :: dbeta(:, :)   ! covariates
     real(dp), allocatable :: dphi(:, :)    ! ar
     real(dp), allocatable :: dtheta(:, :)  ! ma
     real(dp), allocatable :: dd(:, :)      ! d
     real(dp), allocatable :: dthetaT(:, :) ! chaotic
  end type deta_d

  !-------------------------------------------------------------
  !
  ! Matrices for the score vectors U(gamma)
  !
  !-------------------------------------------------------------
  type score
     ! Last revision: July, 2023
     !  - Removed Unu(:)
     real(dp), allocatable :: Ualpha(:)  ! constant
     real(dp), allocatable :: Ubeta(:)   ! covariates
     real(dp), allocatable :: Uphi(:)    ! ar
     real(dp), allocatable :: Utheta(:)  ! ma
     real(dp), allocatable :: Ud(:)      ! d
     real(dp), allocatable :: UthetaT(:) ! chaotic
  end type score

  !-------------------------------------------------------------
  !
  !  The time series for each part
  !
  !-------------------------------------------------------------
  type conditional_ts
     ! July, 2023: added vt(:)
     !
     ! Last revision: February, 2025
     !  - renamed vt to w
     !  - renamed ut to gw
     !  - added estart and et
     !  - added gi1 and gi2
     !
     integer :: n      ! size of wt (the conditional time series)
     integer :: nreg   ! number of regressors
     integer :: xregar ! Indicates if xreg must be included in the AR recursion
     real(dp), allocatable :: xstart(:)  ! starting values for xreg
     real(dp) :: g2start                 ! starting values for g12(yt) or g22(gwt)
     real(dp) :: estart                  ! starting values for et
     real(dp), allocatable :: w(:)       ! omega_it, the original parameter, e.g, mu_t
     real(dp), allocatable :: gw(:)      ! g_i(omega_it), the transformed parameter
     real(dp), allocatable :: eta(:)     ! the linear predictor eta_it = g_i1(gwt)
     real(dp), allocatable :: gi1(:)      ! g_{i1}(.), eg. g11(yt)
     real(dp), allocatable :: gi2(:)      ! g_{i2}(.), eg. g12(yt)
     real(dp), allocatable :: xreg(:, :) ! regressors
     real(dp), allocatable :: orbit(:)   ! T**t(u0), for BARC models
     real(dp), allocatable :: et(:)      ! the error term rt or g23(rt)
  end type conditional_ts

  !-------------------------------------------------------------
  !
  !  model parameters and configurations for each part
  !
  !-------------------------------------------------------------
  type model_part
     ! February, 2025: added this types
     !
     ! parameters realated configurations
     integer :: npar = 0           ! number of non-fixed parameters
     integer :: inf = 0            ! number of terms in the infinite sum expansion
     integer :: skip = 0           ! 1: mu is known npar = 0
     integer :: map = 0            ! chaotic transformation
     type(vec_parameter) :: alpha  ! constant
     type(vec_parameter) :: beta   ! regressors
     type(vec_parameter) :: ar     ! AR recursion
     type(vec_parameter) :: ma     ! MA recursion
     type(vec_parameter) :: d      ! long - memory
     type(vec_parameter) :: thetaT ! constant for BARC transformation
     type(vec_parameter) :: u0     ! starting point for BARC iteration
     type(argsLink) :: linkg(4)    ! link functions g, gi1, gi2, gi3
     type(argsLink) :: linkh       ! link for BARC models: h(T**t(U0))
  end type model_part


  !-------------------------------------------------------------
  !
  ! Generic Vector
  !
  !-------------------------------------------------------------
  type vetor
     ! Last revision: July, 2025
     ! - renamed par to z
     !
     real(dp), allocatable :: z(:)
  end type vetor

  !-------------------------------------------------------------
  !
  !  Score vector and related matrices
  !
  !-------------------------------------------------------------
  type argsSI
     ! Last revision: July, 2023
     !--------------------------------------------
     ! deta(1, 1) deta1_drho;     deta(2, 1) deta2_drho
     ! deta(1, 2) deta1_dlambda;  deta(2, 2) deta2_dlambda
     !--------------------------------------------
     type(deta_d) :: deta(2, 2)
     !----------------------
     ! Score Vectors
     !----------------------
     type(score) :: U(2)
     !----------------------
     ! T matrices
     !----------------------
     type(vetor) :: T(2)
     !----------------------
     ! h vectors
     !----------------------
     type(vetor) :: h(2)
     !----------------------
     ! E matrix
     !----------------------
     real(dp), allocatable :: E(:, :)
  end type argsSI

  !-------------------------------------------------------------
  !
  !  Model specification
  !
  !-------------------------------------------------------------
  type :: argsModel
     ! July, 2023
     ! - removed nu and fixnu
     ! - added skipmu and skipnu
     ! - argsL now has size 6 instead of allocatable
     !
     ! Last revision: February, 2025
     ! - gy and error are now part of the conditional_ts object
     ! - link is now part of the conditional_ts object
     ! - added the model_part type and moved the configurations realated
     !   to each part of the model there.
     !-------------------------------------------------------------------------
     ! general configurations
     integer               :: n      ! sample size
     real(dp), allocatable :: y(:)   ! original time series
     real(dp)              :: ystart ! starting value for y

     ! part related configurations
     type(conditional_ts)  :: cts(2) ! conditional series and related regressors
     type(model_part)      :: pt(2)  ! model configurations for each part

     ! distribution related configurations
     type(argsDist)   :: argsD    ! arguments and functions related to the distribution
     integer          :: m = 0    ! starting point for log - likelihood
     integer          :: llk = 0  ! 1 = calculate loglikelihood
     integer          :: sco = 0  ! 1 = calculate score vector
     integer          :: info = 0 ! 1 = calculate information matrix
     type(argsSI)     :: SI       ! score vector and information matrix
     type(par_bounds) :: bounds   ! parameter bounds for optimization
  end type argsModel

  type :: optimFunc
     character(len = 60) :: message
     logical :: dummy = .true. ! to avoid compiler warning
     procedure(llk_generic),        pointer :: loglik => null()
     procedure(llk_nelder_generic), pointer :: functn => null()
  end type optimFunc

  abstract interface
     subroutine llk_generic(loglik, model, nop, p, func, U)
       import :: dp, optimFunc, argsmodel
       implicit none
       class(optimFunc), intent(inout) :: loglik
       type(argsmodel), intent(inout) :: model
       integer, intent(in) :: nop
       real (dp), intent(in)  :: p(nop)
       real (dp), intent(out) :: func, U(nop)
     end subroutine llk_generic
  end interface


  abstract interface
     subroutine llk_nelder_generic(loglik, model, nop, p, func)
       import :: dp, optimFunc, argsmodel
       implicit none
       class(optimFunc), intent(inout) :: loglik
       type(argsmodel), intent(inout) :: model
       integer, intent(in) :: nop
       real (dp), intent(in)  :: p(nop)
       real (dp), intent(out) :: func
     end subroutine llk_nelder_generic
  end interface

contains
  !-------------------------------------------------------------------------------------------
  !
  !               Subroutines for parameter transformation
  !
  !
  ! Based on the Matlab subroutine fminsearchbnd
  ! Reference:
  ! https://www.mathworks.com/matlabcentral/fileexchange/8277-
  !                 fminsearchbnd-fminsearchcon?focused=5216898&tab=function
  !
  !-------------------------------------------------------------------------------------------
  subroutine set_bounds(bounds, bds, nbd, npar)
    !-------------------------------------------------------------------------------
    ! allocates the vectors
    !    lower bound, upper bound and type of bounds (nbd)
    ! in the bounds object and save the values passed by the user.
    !-------------------------------------------------------------------------------
    ! Last revision: July, 2023
    implicit none
    type(par_bounds), intent(inout) :: bounds
    integer, intent(in) :: npar
    integer, intent(in) :: nbd(npar)
    real(dp), intent(in) :: bds(npar, 2)

    call safe_allocate(bounds%lower, npar)
    call safe_allocate(bounds%upper, npar)
    call safe_allocate(bounds%nbd, npar)

    bounds%lower = bds(:,1)
    bounds%upper = bds(:,2)
    bounds%nbd = nbd
    return
  end subroutine set_bounds

  function xtransform(npars, x, bounds) result(xtrans)
    !-------------------------------------------------------------------------------
    ! converts unconstrained variables into their original domains
    ! bounds type:
    !  0 = none
    !  1 = lower
    !  2 = both
    !  3 = upper
    !  4 = constant
    !-------------------------------------------------------------------------------
    ! Last revision: February, 2025
    !  -  replaced the loop by a vectorized version
    implicit none
    type(par_bounds), intent(in) :: bounds
    integer, intent(in) :: npars
    real(dp), intent(in) :: x(npars)
    real(dp) :: xtrans(npars)
    real(dp) :: temp(npars)  ! Temporary array for intermediate calculations

    ! Initialize xtrans with default values (case 0: unconstrained)
    xtrans = x

    ! Case 1: Lower bound only
    where (bounds%nbd == 1)
       xtrans = bounds%lower + x**2
    end where

    ! Case 2: Both lower and upper bounds
    where (bounds%nbd == 2)
       temp = (sin(x) + 1.d0) / 2.d0
       xtrans = temp * (bounds%upper - bounds%lower) + bounds%lower
       ! just in case of any floating point problems
       xtrans = max(bounds%lower, min(bounds%upper, xtrans))
    end where

    ! Case 3: Upper bound only
    where (bounds%nbd == 3)
       xtrans = bounds%upper - x**2
    end where

    ! Case 4: Constant (fixed variable)
    where (bounds%nbd == 4)
       xtrans = bounds%lower
    end where

    return
  end function xtransform

  function xtransformstart(npars, x, bounds) result(xtrans)
    !--------------------------------------------------------
    ! transform starting values into their unconstrained
    ! surrogates. check for infeasible starting guesses.
    !--------------------------------------------------------
    ! Last revision: February, 2025
    !  -  conv
    implicit none
    type(par_bounds), intent(in) :: bounds
    integer, intent(in) :: npars
    real(dp), intent(in) :: x(npars)
    real(dp) :: xtrans(npars)
    real(dp), parameter :: pi_half = pi / 2, two_pi = 2 * pi
    real(dp) :: temp(npars)  ! Temporary array for intermediate calculations

    ! Initialize xtrans with default values (case 0: unconstrained)
    xtrans = x

    ! Case 1: Lower bound only
    where (bounds%nbd == 1)
       where (x <= bounds%lower)
          ! Infeasible starting value, use bound
          xtrans = 0.d0
       elsewhere
          xtrans = sqrt(x - bounds%lower)
       end where
    end where

    ! Case 2: Both lower and upper bounds
    where (bounds%nbd == 2)
       where (x <= bounds%lower)
          ! Infeasible starting value
          xtrans = - pi_half
       elsewhere (x >= bounds%upper)
          ! Infeasible starting value
          xtrans = pi_half
       elsewhere
          temp = 2.d0 * (x - bounds%lower) / (bounds%upper - bounds%lower) - 1.d0
          xtrans = two_pi + asin(max(-1.d0, min(1.d0, temp)))
       end where
    end where

    ! Case 3: Upper bound only
    where (bounds%nbd == 3)
       where (x >= bounds%upper)
          ! Infeasible starting value, use bound
          xtrans = 0.d0
       elsewhere
          xtrans = sqrt(bounds%upper - x)
       end where
    end where

    ! Case 4: Constant (fixed variable)
    where (bounds%nbd == 4)
       xtrans = bounds%lower
    end where

    return
  end function xtransformstart

  subroutine transform_par(par, npar, bounds, inverse)
    !-------------------------------------------------------------------------------
    ! transforms the y = par (bounded) in x = f(y) (unbounded)
    ! if inverse = .true., transforms x (unbounded) in y = f^{-1}(x) (bounded)
    !-------------------------------------------------------------------------------
    ! Last revision: July, 2023
    implicit none
    integer, intent(in) :: npar
    real(dp), intent(inout) :: par(npar)
    type(par_bounds), intent(in) :: bounds
    logical, intent(in) :: inverse

    if (sum(bounds%nbd) == 0) return

    ! unbounded  - > bounded
    if (inverse) then
       par = xtransform(npar, par, bounds)
       return
    end if

    ! bounded  - > unbounded
    par = xtransformstart(npar, par, bounds)
    return
  end subroutine transform_par

end module main_mod
