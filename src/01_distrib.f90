module distrib
  use lib_utils ! parameters and special functions
  !****************************************************************************
  !
  ! This module contains distributions related subroutines
  !
  !****************************************************************************
  ! Major update: July, 2023
  !  - Removed several subroutines and replaced them with ones provided in R.
  !  - mu and nu now enter the subroutines with size n instead of n1 and n2
  !  - changed the code accordingly
  !
  ! Last update: March, 2025
  !  - added the procedures in the argsDist argument
  !  - moved all distribution related subrotutines to the same module
  !****************************************************************************
  implicit none
  private
  public :: rng_uniform
  public :: argsDist

  real(dp), private, parameter :: em = 0.57721566490153286061d0 ! Euler-Mascheroni

  !-------------------------------------------------------------
  !
  !  Arguments for distribution related function
  !
  !-------------------------------------------------------------
  type :: argsDist
     !
     ! Last revision: March, 2025
     !  - added the procedures
     !
     logical  :: dummy = .true. ! dummy to avoid compiler warning
     real(dp) :: lower = 0.d0  ! y.lower
     real(dp) :: upper = 1.d0  ! y.upper
     real(dp) :: par = 0.5d0   ! extra arguments
     character(len = 8) :: model

     procedure(rdist_generic),      pointer :: rdist       => null()
     procedure(llkdist_generic),    pointer :: llk_dist    => null()
     procedure(dllkdist_generic),   pointer :: dllk_dist   => null()
     procedure(Ed2llkdist_generic), pointer :: Ed2llk_dist => null()

   contains
     procedure :: init => init_mydist  ! Initialization subroutine
  end type argsDist


  !-------------------------------------------------------------------------------------
  ! Define the abstract interfaces
  !  - rdist_generic: random number generator
  !  - llkdist_generic: the sum of the loglikelihood
  !  - dllkdist_generic: computes the loglikelihood derivatives
  !  - Ed2llkdist_generic: computes the expected value of the second derivative
  !--------------------------------------------------------------------------------------
  abstract interface
     function rdist_generic(argsD, npar, par) result(fnval)
       import :: dp, argsDist
       implicit none
       class(argsDist), intent(inout)  :: argsD
       integer,  intent(in)         :: npar
       real(dp), target, intent(in) :: par(npar)
       real(dp)  ::  fnval
     end function rdist_generic
  end interface

  abstract interface
     function llkdist_generic(argsD, m, n, y, mu, nu) result(sll)
       !********************************************************************
       ! llk_dis Computes
       !           sll = s'um(log(f(y| mu, nu)))
       !  where f(y| mu, nu) is the corresponding density function
       !********************************************************************
       ! Input
       !   argsD: aditional arguments to be passed to the density function
       !   m    : starting time to calculate the log - likelihood.
       !          For t < m + 1, sll = 0
       !   n    : sample size
       !   y    : the observed time series
       !   mu   : the conditional time series mut
       !   nu   : the conditional time series nut
       !
       ! Output
       !   sll: sum of the log-likelihood from times m + 1 to n.
       !********************************************************************
       import :: dp, argsDist
       implicit none
       class(argsDist), intent(inout) :: argsD
       integer,  intent(in)        :: m, n
       real(dp), intent(in)        :: y(n), mu(n), nu(n)
       real(dp)                    :: sll
     end function llkdist_generic
  end interface

  abstract interface
     subroutine dllkdist_generic(argsD, m, n, y, mu, nu, skip, dllmu, dllnu)
       !**********************************************************************
       !  dllk_dist calculates dl/dmu and dl/dnu
       !**********************************************************************
       ! Input
       !   argsD: aditional arguments to be passed to the density function
       !   m    : starting time to calculate the derivatives.
       !          For t < m + 1, dl/dmu = dl/dnu = 0
       !   n    : sample size
       !   y    : the observed time series
       !   mu   : the conditional time series mut
       !   nu   : the conditional time series nut
       !   skip : indicates if dl/dmu and/or dl/dmu must be calculated
       !           skip(1) = 1 indicates that mu is known in the model
       !           skip(2) = 1 indicates that mu is known in the model
       !
       ! Output:
       !   dllmu and dldnu: the derivatives of the log-likelihood function
       !                    with respect to mu and nu, respectively
       !**********************************************************************
       import :: dp, argsDist
       implicit none
       class(argsDist), intent(inout) :: argsD
       integer,  intent(in)        :: m, n, skip(2)
       real(dp), intent(in)        :: y(n), mu(n), nu(n)
       real(dp), intent(out)       :: dllmu(min(n, n * (1 - skip(1)) + 1))
       real(dp), intent(out)       :: dllnu(min(n, n * (1 - skip(2)) + 1))
     end subroutine dllkdist_generic
  end interface

  abstract interface
     subroutine Ed2llkdist_generic(argsD, m, n, mu, nu, skip, E)
       !**********************************************************************
       !  Ed2llk_dist calculates
       !           E(dl2/dmu2), E(dl2/dmUdnu) and E(dl2/dnu2)
       !  the exepcted values of the second derivatives of the log-likelihood
       !**********************************************************************
       ! Input
       !   argsD: aditional arguments to be passed to the density function
       !   m    : starting time to calculate the log-likelihood.
       !           For t < m + 1, sll = 0
       !   n    : sample size
       !   mu   : the conditional time series mut
       !   nu   : the conditional time series nut
       !   skip : indicates if dl/dmu and/or dl/dmu must be calculated
       !           skip(1) = 1 indicates that mu is known in the model
       !           skip(2) = 1 indicates that mu is known in the model
       !
       ! Output
       !    E: matrix with the expected values
       !**********************************************************************
       import :: dp, argsDist
       implicit none
       class(argsDist), intent(inout) :: argsD
       integer,  intent(in)        :: m, n, skip(2)
       real(dp), intent(in)        :: mu(n), nu(n)
       real(dp), intent(out)       :: E(n, 3)
     end subroutine Ed2llkdist_generic
  end interface

contains
  !***************************************************************************
  ! Initialization subroutine to set the procedure pointer based on the code
  !***************************************************************************
  subroutine init_mydist(argsD, model)
    class(argsDist), intent(inout) :: argsD
    character(len = 8) :: model
    argsD%model = model
    select case (trim(adjustl(argsD%model)))
    case("beta")
       argsD%rdist       => rbeta
       argsD%llk_dist    => llk_beta
       argsD%dllk_dist   => dllk_beta
       argsD%Ed2llk_dist => Ed2llk_beta
    case("gamma")
       argsD%rdist       => rgamma
       argsD%llk_dist    => llk_gamma
       argsD%dllk_dist   => dllk_gamma
       argsD%Ed2llk_dist => Ed2llk_gamma
    case("kuma")
       argsD%rdist       => rkuma
       argsD%llk_dist    => llk_kuma
       argsD%dllk_dist   => dllk_kuma
       argsD%Ed2llk_dist => Ed2llk_kuma
    case("matsu")
       argsD%rdist       => rmatsu
       argsD%llk_dist    => llk_matsu
       argsD%dllk_dist   => dllk_matsu
       argsD%Ed2llk_dist => Ed2llk_matsu
    case("ul")
       argsD%rdist       => rul
       argsD%llk_dist    => llk_ul
       argsD%dllk_dist   => dllk_ul
       argsD%Ed2llk_dist => Ed2llk_ul
    case("uw")
       argsD%rdist       => ruw
       argsD%llk_dist    => llk_uw
       argsD%dllk_dist   => dllk_uw
       argsD%Ed2llk_dist => Ed2llk_uw
    case default
       argsD%rdist       => rbeta
       argsD%llk_dist    => llk_beta
       argsD%dllk_dist   => dllk_beta
       argsD%Ed2llk_dist => Ed2llk_beta
       return
    end select
  end subroutine init_mydist

  !****************************************************************************
  ! Uniform Random Number Generation
  !****************************************************************************
  function rng_uniform() result(y)
    !
    ! generates uniform random variates using unif_rand() from R.
    !
    ! Implemented by Taiane S. Prass - PPGEst/UFRGS
    ! July, 2023 - Porto Alegre
    !
    implicit none
    real (dp) :: y, unifrnd
    call rndstart()
    y = unifrnd()
    call rndend()
    return
  end function rng_uniform

  !*****************************************************************************
  !
  ! Beta Distribution
  !
  !*****************************************************************************
  function rbeta(argsD, npar, par) result(y)
    !**********************************************************************
    ! Generates beta random variates parametrized by the mean and nu.
    ! Uses rbeta() from R which implements the traditional parametrization
    !
    ! Implemented by Taiane S. Prass - PPGEst/UFRGS
    ! July, 2023 - Porto Alegre
    !**********************************************************************
    ! Last revision: February, 2025
    !  - added pointers
    !**********************************************************************
    implicit none
    class(argsDist), intent(inout)  :: argsD
    integer,  intent(in)         :: npar   ! must be 2
    real(dp), target, intent(in) :: par(npar)
    real(dp), pointer :: mu, nu
    real(dp) :: y, betarnd
    argsD%dummy = .true.
    mu => par(1)
    nu => par(2)
    call rndstart()
    y = betarnd(mu * nu, (1.d0 - mu) * nu)
    call rndend()
    return
  end function rbeta

  function dbeta(y, npar, par, give_log) result(fn_val)
    !**********************************************************************
    ! Computes the density parametrized by the mean and nu.
    ! Uses dbeta() from R which implements the traditional parametrization
    !
    ! Implemented by Taiane S. Prass - PPGEst/UFRGS
    ! July, 2023 - Porto Alegre
    !**********************************************************************
    ! Last revision: February, 2025
    !  - added pointers
    !**********************************************************************
    implicit none
    integer,  intent(in) :: npar
    real(dp), intent(in) :: y
    real(dp), target, intent(in) :: par(npar)
    real(dp), pointer    :: mu, nu
    integer,  intent(in) :: give_log
    real(dp) :: fn_val, betadens
    mu => par(1)
    nu => par(2)
    fn_val = betadens(y, mu * nu, (1.d0 - mu) * nu, give_log)
    return
  end function dbeta

  function llk_beta(argsD, m, n, y, mu, nu) result(sll)
    !------------------------------------------------
    ! Computes sll for beta regression models
    !------------------------------------------------
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in) :: m, n
    real(dp), intent(in) :: y(n), mu(n), nu(n)
    real(dp) :: sll
    integer  :: t
    argsD%dummy = .true.
    sll = 0.d0
    do t = (m  + 1), n
       sll = sll + dbeta(y(t), 2, [mu(t), nu(t)], 1)
    end do
    return
  end function llk_beta

  subroutine dllk_beta(argsD, m, n, y, mu, nu, skip, dllmu, dllnu)
    !**********************************************************************
    !  Computes dl/dmu and dl/dnu for beta regression models
    !**********************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in)  :: m, n, skip(2)
    real(dp), intent(in)  :: y(n), mu(n), nu(n)
    real(dp), intent(out) :: dllmu(min(n, n * (1 - skip(1)) + 1))
    real(dp), intent(out) :: dllnu(min(n, n * (1 - skip(2)) + 1))
    real(dp) :: log1y(n), ymustar(n), dig(n)
    integer  :: t

    argsD%dummy = .true.
    dllmu = 0.d0
    dllnu = 0.d0
    if (sum(skip) == 2) return

    !----------------------------------------------------------------------
    ! ymustart = ystar - mustar
    ! dig = digamma((1 - mu) * nu)
    !----------------------------------------------------------------------
    log1y = log(1.d0 - y)
    do t = (m + 1), n
       dig(t) = digamma((1.d0 - mu(t)) * nu(t))
       ymustar(t) = log(y(t)) - log1y(t) - digamma(mu(t) * nu(t)) + dig(t)
    end do

    !----------------------------------------------------------------------
    ! dl/dmu
    !----------------------------------------------------------------------
    if (skip(1) == 0) then
       do t = (m  + 1), n
          dllmu(t) = nu(t) * ymustar(t)
       end do
    end if

    !----------------------------------------------------------------------
    ! dl/dnu
    !----------------------------------------------------------------------
    if (skip(2) == 0) then
       do t = (m  + 1), n
          dllnu(t) = mu(t) * ymustar(t) + log1y(t) - dig(t) + digamma(nu(t))
       end do
    end if
    return
  end subroutine dllk_beta

  subroutine Ed2llk_beta(argsD, m, n, mu, nu, skip, E)
    !**********************************************************************
    !  Computes
    !          E(dl2/dmu2), E(dl2/dmUdnu) and E(dl2/dnu2)
    !  for beta regression models
    !**********************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in)  :: m, n, skip(2)
    real(dp), intent(in)  :: mu(n), nu(n)
    real(dp), intent(out) :: E(n, 3)
    real(dp) :: psi1, psi2
    integer  :: t, ifault

    argsD%dummy = .true.
    E = 0.d0
    if (sum(skip) == 2) return

    do t = (m + 1), n
       psi1 = trigamma(mu(t) * nu(t), ifault)
       psi2 = trigamma((1.d0 - mu(t)) * nu(t), ifault)
       if (skip(1) == 0) E(t, 1) = (psi1 + psi2) * nu(t)**2
       if (sum(skip) == 0) E(t, 2) = ((psi1 + psi2) * mu(t) - psi2) * nu(t)
       if (skip(2) == 0) E(t, 3) =  psi1 * mu(t)**2 + psi2 * (1.d0 - mu(t))**2  - &
            trigamma(nu(t), ifault)
    end do
    return
  end subroutine Ed2llk_beta

  !*****************************************************************************
  !
  ! Kumaraswamy Distribution
  !
  !*****************************************************************************
  function rkuma_default(delta, nu, a, b) result(y)
    !**********************************************************************
    ! Generates Kumaraswamy random variates using the inversion method.
    !
    ! nu/varphi > 0
    ! delta > 0
    ! a  = lower limit
    ! b  = upper limit
    !
    ! Implemented by Taiane S. Prass - PPGEst/UFRGS
    ! April, 2020 - Porto Alegre
    !**********************************************************************
    ! Last revision: July, 2023
    !   - Removed the "rng" argument
    !**********************************************************************
    implicit none
    real(dp), intent(in) :: delta, nu,  a, b
    real(dp) :: u, y
    u = rng_uniform()
    u = (1.d0 / delta) * log(1.d0 - u) ! Avoid overflow / underflow
    u = exp(u)                         ! z = (1 - u)**(1 / delta)
    y = (1.d0 / nu) * log(1.d0 - u)    ! x = (1 - z)**(1 / nu)
    y = a + (b - a) * exp(y)           ! y = a + (b - a) * exp(x)
    return
  end function rkuma_default

  function rkuma(argsD, npar, par) result(y)
    !************************************************************************
    ! Generates Kumaraswamy random variates parametrized by the mean and nu,
    ! par = (mu, nu, a, b, rho)
    !
    ! Implemented by Taiane S. Prass - PPGEst/UFRGS
    ! April, 2020 - Porto Alegre
    !************************************************************************
    ! July, 2023:
    !  - Removed the "rng" argument
    !
    ! Last revision: February, 2025
    !  - added pointers
    !************************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in) :: npar   ! must be 5
    real(dp), target, intent(in) :: par(npar)
    real(dp), pointer    :: mu, nu, a, b, rho
    real(dp) :: y, delta, mu01
    argsD%dummy = .true.
    mu  => par(1)
    nu  => par(2)
    rho => par(3)
    a   => par(4)
    b   => par(5)
    mu01 = (mu - a) / (b - a)
    delta = log(1.d0 - rho) / log(1.d0 - mu01**nu)
    y = rkuma_default(delta, nu, a, b)
    return
  end function rkuma

  function dkuma(y, npar, par, give_log) result(fn)
    !**********************************************************************
    ! Density funtion - Kumaraswamy distribution parametrized
    ! by the mean and nu.
    !
    ! a, b real numbers
    ! mu = rho - quantile, 0 < rho < 1
    ! par = c(mu, nu, rho, a, b)
    !
    ! Implemented by Taiane S. Prass - PPGEst/UFRGS
    ! April, 2020
    !**********************************************************************
    !  April, 2021:
    !    Replaced original dkuma by dkuma and dkuma_default
    !
    ! Last revision: February, 2025
    !  - added pointers
    !**********************************************************************
    implicit none
    integer,  intent(in) :: npar   ! must be 5
    real(dp), target, intent(in) :: par(npar)
    real(dp), pointer    :: mu, nu, a, b, rho
    real(dp), intent(in) :: y
    integer,  intent(in) :: give_log
    real(dp) :: fn, mu01, delta
    mu  => par(1)
    nu  => par(2)
    rho => par(3)
    a   => par(4)
    b   => par(5)
    mu01 = (mu - a) / (b - a)
    delta = log(1.d0 - rho) / log(1.d0 - mu01**nu)
    fn = dkuma_default(y, delta, nu, a, b, give_log)
    return
  end function dkuma

  function dkuma_default(y, delta, nu, a, b, give_log) result(fn)
    !**********************************************************************
    ! Computes the density of the Kumaraswamy distribution
    ! (traditional parametrization).
    !
    ! Implemented by Taiane S. Prass - PPGEst/UFRGS
    ! April, 2020
    !**********************************************************************
    ! Last revision: April, 2021
    !    Replaced original dkuma by dkuma and dkuma_default
    !**********************************************************************
    implicit none
    real(dp), intent(in) :: y, delta, nu, a, b
    integer,  intent(in) :: give_log
    real(dp) :: fn, b_minus_a, y_minus_a, log_b_minus_a
    b_minus_a = b - a
    log_b_minus_a = log(b - a)
    y_minus_a = y - a
    fn = log(nu) + log(delta) - log_b_minus_a
    fn = fn + (nu - 1.d0) * (log(y_minus_a) - log_b_minus_a)
    fn = fn + (delta - 1.d0) * log(1.d0 - ((y_minus_a) / (b_minus_a))**nu)
    if (give_log == 0) fn = exp(fn)
    return
  end function dkuma_default

  function llk_kuma(argsD, m, n, y, mu, nu) result(sll)
    !***********************************************************************************
    ! Computes sll for Kumaraswamy regression models
    ! argsD contains the quantile of interest rho, for which mu is the rho-th quantile
    !***********************************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in) :: m, n
    real(dp), intent(in) :: y(n), mu(n), nu(n)
    real(dp) :: sll
    integer  :: t
    argsD%dummy = .true.
    sll = 0.d0
    do t = (m  + 1), n
       sll = sll + dkuma(y(t), 5, [mu(t), nu(t), argsD%par, argsD%lower, argsD%upper], 1)
    end do
    return
  end function llk_kuma

  subroutine dllk_kuma(argsD, m, n, y, mu, nu, skip, dllmu, dllnu)
    !*************************************************************************
    !  Computes dl/dmu and dl/dnu for Kumaraswamy regression models
    !*************************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in)  :: m, n, skip(2)
    real(dp), intent(in)  :: y(n), mu(n), nu(n)
    real(dp), intent(out) :: dllmu(min(n, n * (1 - skip(1)) + 1))
    real(dp), intent(out) :: dllnu(min(n, n * (1 - skip(2)) + 1))
    real(dp) :: ct(n), mu1, y01(n), delta(n), mu01(n)
    integer  :: t

    argsD%dummy = .true.
    dllmu = 0.d0
    dllnu = 0.d0
    if (sum(skip) == 2) return

    y01 = (y - argsD%lower) / (argsD%upper - argsD%lower)
    mu01 = (mu - argsD%lower) / (argsD%upper - argsD%lower)
    delta = 0.d0
    ct = 0.d0

    !----------------------------------------------------------------------
    ! calculating delta(t) and c(t)
    !----------------------------------------------------------------------
    do t = (m + 1), n
       mu1 =  1.d0 - mu01(t)**nu(t)
       delta(t) = log(1.d0 - argsD%par) / log(mu1)
       ct(t) = mu01(t)**(nu(t) - 1.d0) / (mu1 * log(mu1))
       ct(t) = ct(t) * (delta(t) * log(1.d0 - y01(t)**nu(t)) + 1.d0)
    end do

    !----------------------------------------------------------------------
    ! dl/dmu
    !----------------------------------------------------------------------
    if (skip(1) == 0) then
       dllmu = nu * ct /  (argsD%upper - argsD%lower)
    end if

    !----------------------------------------------------------------------
    ! dl/dnu
    !----------------------------------------------------------------------
    if (skip(2) == 0) then
       do t = (m + 1), n
          dllnu(t) =  1 / nu(t) + log(y01(t)) + ct(t) * mu01(t) * log(mu01(t)) - &
               (delta(t) - 1.0d0) * y01(t)**nu(t) * log(y01(t)) / (1.d0 - y01(t)**nu(t))
       end do
    end if
    return
  end subroutine dllk_kuma

  subroutine Ed2llk_kuma(argsD, m, n, mu, nu, skip, E)
    !**********************************************************************
    !  Computes
    !       E(dl2/dmu2), E(dl2/dmUdnu) and E(dl2/dnu2)
    !  for Kumaraswamy regression models
    !**********************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in)  :: m, n, skip(2)
    real(dp), intent(in)  :: mu(n), nu(n)
    real(dp), intent(out) :: E(n, 3)
    real(dp) :: mu1, tau1(n), tau2(n), aux
    real(dp) :: delta(n), mu01(n), dgd, ddgd, s
    integer  :: t, ifault
    real(dp), parameter :: k0 = 0.82368066085287928d0

    argsD%dummy = .true.
    E = 0.d0
    if (sum(skip) == 2) return

    mu01 = (mu - argsD%lower) / (argsD%upper - argsD%lower)
    delta = 0.d0
    tau1 = 0.d0
    tau2 = 0.d0
    aux = 0.d0
    do t = (m + 1), n
       mu1 =  1.d0 - mu01(t)**nu(t)
       delta(t) = log(1.d0 - argsD%par) / log(mu1)
       tau1(t) = mu01(t)**(nu(t) - 2.d0) / (mu1 * log(mu1))
       tau2(t) = tau1(t)**2 * mu01(t)**2
    end do

    do t = (m + 1), n
       if (skip(1) == 0) E(t, 1) = tau2(t) * (nu(t) / (argsD%upper - argsD%lower))**2
       ! setting the values used by dl2/dmUdnu and dl2/dnu2
       if (skip(2) == 0) then
          dgd = digamma(delta(t))
          ddgd = trigamma(delta(t), ifault)
          mu1 = mu01(t) * log(mu01(t))
          ! digamma(x + 1) = digamma(x) + 1 / x
          aux = (1.d0 - dgd - 1.d0 / delta(t) - em) / (delta(t) - 1.d0)
          aux = delta(t) * tau1(t) * aux
          if (skip(1) == 0) then
             E(t, 2) = E(t, 1) * mu1 * (argsD%upper - argsD%lower) / nu(t)
             E(t, 2) = E(t, 2) + aux * mu01(t) / (argsD%upper - argsD%lower)
          end if
          s = 2 * aux * mu1 * mu01(t) / nu(t)
          E(t, 3) = 1.d0 / nu(t)**2 + tau2(t) * mu1**2  + s
          s = (dgd * (dgd + 2.d0 * (em - 1.d0)) - ddgd + k0) / nu(t)**2
          E(t, 3) = E(t, 3) + s * delta(t) / (delta(t) - 2.d0)
       end if
    end do
    return
  end subroutine Ed2llk_kuma

  !*****************************************************************************
  !
  ! Gamma Distribution
  !
  !*****************************************************************************
  function rgamma(argsD, npar, par) result(y)
    !*********************************************************************
    ! Generates gamma random variates parametrized by the mean and nu.
    !
    ! Uses rgamma() from R to generate random deviates from the gamma
    ! distribution whose density is
    !    1/((scale**shape)*Gamma(shape))*x**(shape-1)*Exp(-x/scale)
    ! Arguments
    !   scale --> Location parameter of Gamma distribution  (scale > 0)
    !   shape --> Shape parameter of Gamma distribution (shape > 0)
    !
    ! Implemented by Taiane S. Prass - PPGEst/UFRGS
    ! July, 2023 - Porto Alegre
    !*********************************************************************
    ! Last revision: February, 2025
    !  - added pointers
    !*********************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in) :: npar ! must be 2
    real(dp), target, intent(in) :: par(npar)
    real(dp), pointer :: mu, nu
    real(dp) :: y, gammarnd
    argsD%dummy = .true.
    mu => par(1)
    nu => par(2)
    call rndstart()
    y = gammarnd(nu, mu / nu)
    call rndend()
    return
  end function rgamma

  function d_gamma(y, npar, par, give_log) result(fn)
    !*********************************************************************
    ! Computes the density of the gamma distribution parametrized by
    ! the mean and nu.
    !
    ! Uses dgamma() from R which implements
    !  1/((scale**shape)*Gamma(shape))*x**(shape - 1)*Exp(-x/scale)
    !
    ! Implemented by Taiane S. Prass - PPGEst/UFRGS
    ! July, 2023 - Porto Alegre
    !*********************************************************************
    ! Last revision: February, 2025
    !  - added pointers
    !*********************************************************************
    implicit none
    integer,  intent(in) :: npar
    real(dp), intent(in) :: y
    real(dp), target, intent(in) :: par(npar)
    real(dp), pointer    :: mu, nu
    integer,  intent(in) :: give_log
    real(dp) :: fn, gammadens
    mu => par(1)
    nu => par(2)
    fn = gammadens(y, nu, mu / nu, give_log)
    return
  end function d_gamma

  function llk_gamma(argsD, m, n, y, mu, nu) result(sll)
    !***********************************************************************
    ! Computes sll for gamma regression models
    !***********************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in) :: m, n
    real(dp), intent(in) :: y(n), mu(n), nu(n)
    real(dp) :: sll
    integer  :: t
    argsD%dummy = .true.
    sll = 0.d0
    do t = (m  + 1), n
       sll = sll + d_gamma(y(t), 2, [mu(t), nu(t)], 1)
    end do
    return
  end function llk_gamma

  subroutine dllk_gamma(argsD, m, n, y, mu, nu, skip, dllmu, dllnu)
    !*************************************************************************
    !  Computes dl/dmu and dl/dnu for gamma regression models
    !*************************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in) :: m, n, skip(2)
    real(dp), intent(in) :: y(n), mu(n), nu(n)
    real(dp), intent(out) :: dllmu(min(n, n * (1 - skip(1)) + 1))
    real(dp), intent(out) :: dllnu(min(n, n * (1 - skip(2)) + 1))
    integer  :: t
    real(dp) :: ym

    argsD%dummy = .true.
    dllnu = 0.d0
    if (sum(skip) == 2) return

    !----------------------------------------------------------------------
    ! dl/dmu
    !----------------------------------------------------------------------
    if (skip(1) == 0) then
       do t = (m  + 1), n
          dllmu(t) = nu(t) / mu(t) * (y(t) / mu(t)  - 1.d0)
       end do
    end if

    !----------------------------------------------------------------------
    ! dl/dnu
    !----------------------------------------------------------------------
    if (skip(2) == 0) then
       do t = (m  + 1), n
          ym = y(t) / mu(t)
          dllnu(t) = 1.d0 - digamma(nu(t)) + log(nu(t) * ym) - ym
       end do
    end if
    return
  end subroutine dllk_gamma

  subroutine Ed2llk_gamma(argsD, m, n, mu, nu, skip, E)
    !**********************************************************************
    !  Computes
    !      E(dl2/dmu2), E(dl2/dmUdnu) and E(dl2/dnu2)
    !  for gamma regression models
    !**********************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in)  :: m, n, skip(2)
    real(dp), intent(in)  :: mu(n), nu(n)
    real(dp), intent(out) :: E(n, 3)
    integer :: t, ifault

    argsD%dummy = .true.
    E = 0.d0
    if (sum(skip) == 2) return

    do t = (m + 1), n
       if (skip(1) == 0) E(t, 1) = nu(t) / mu(t)**2
       if (skip(2) == 0) E(t, 3) = trigamma(nu(t), ifault) - 1 / nu(t)
    end do
    return
  end subroutine Ed2llk_gamma

  !*****************************************************************************
  !
  ! Unit Weibull Distribution
  !
  !*****************************************************************************
  function ruw(argsD, npar, par) result(y)
    !*********************************************************************
    ! Generates Unit Weibull random variates using the inversion method.
    ! par = c(mu, nu, rho)
    !
    ! 0 < rho < 1, 0 < mu < 1, nu > 0
    !
    ! Implemented by Guilherme Pumi - PPGEst/UFRGS
    ! September, 2021
    !*********************************************************************
    ! Revision: Taiane Schaedler Prass - PPGEst/UFRGS
    ! July, 2023 - Porto Alegre
    !  - Removed the "rng" argument
    !
    ! Last revision: February, 2025
    !  - renamed lambda to nu
    !  - added pointers
    !*********************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in) :: npar   ! must be 3
    real(dp), target, intent(in) :: par(npar)
    real(dp), pointer    :: mu, nu, rho
    real(dp) :: y, u
    argsD%dummy = .true.
    mu  => par(1)
    nu  => par(2)
    rho => par(3)
    u = rng_uniform()
    y = log(mu) * (log(u) / log(rho))**(1 / nu)
    y = exp(y)
    return
  end function ruw

  function duw(y, npar, par, give_log) result(fn)
    !**********************************************************************
    ! Computes the density of the Unit Weibull distribution.
    ! par = c(mu, lambda, rho)
    !
    ! 0 < rho < 1,    0 < mu < 1 is a scale parameter
    ! and nu > 0 is a shape parameter.
    ! In this parameterization, mu = rho-th quantile of the
    ! Unit-Weibull distribution with shape parameter lambda
    !
    ! Implemented by Guilherme Pumi - PPGEst/UFRGS
    ! September, 2021
    !**********************************************************************
    ! Last revision: February, 2025
    !  - renamed lambda to nu
    !  - added pointers
    !**********************************************************************
    implicit none
    integer,  intent(in) :: npar
    real(dp), intent(in) :: y
    real(dp), target, intent(in) :: par(npar)
    real(dp), pointer    :: mu, nu, rho
    integer,  intent(in) :: give_log
    real(dp) :: A,  fn
    mu  => par(1)
    nu  => par(2)
    rho => par(3)
    A = log(y) / log(mu)
    fn = log(nu) - log(y) + log(log(rho) / log(mu))
    fn = fn + (nu - 1) * log(A) + log(rho) * (A**nu)
    if (give_log == 0) fn = exp(fn)
    return
  end function duw

  function llk_uw(argsD, m, n, y, mu, nu) result(sll)
    !---------------------------------------------------------------
    ! Computes sll Unit-Weibull regression models
    ! argsD contains the quantile of interest rho, for which mu is
    ! the rho-th quantile
    !---------------------------------------------------------------
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in) :: m, n
    real(dp), intent(in) :: y(n), mu(n), nu(n) ! nu = lambda
    real(dp) :: sll
    integer  :: t
    argsD%dummy = .true.
    sll = 0.d0
    do t = (m + 1), n
       sll = sll + duw(y(t), 3, [mu(t), nu(t), argsD%par], 1)
    end do
    return
  end function llk_uw

  subroutine dllk_uw(argsD, m, n, y, mu, nu, skip, dllmu, dllnu)
    !*************************************************************************
    !  Computes dl/dmu and dl/dnu for Unit-Weibull regression models
    !  argsD contains rho, which is a constant in the model.
    !*************************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in)  :: m, n, skip(2)
    real(dp), intent(in)  :: y(n), mu(n), nu(n)
    real(dp), intent(out) :: dllmu(min(n, n * (1 - skip(1)) + 1))
    real(dp), intent(out) :: dllnu(min(n, n * (1 - skip(2)) + 1))
    real(dp) :: At(n), num(n), denum(n)
    integer  :: t

    argsD%dummy = .true.
    dllmu = 0.d0
    dllnu = 0.d0
    if (sum(skip) == 2) return

    !----------------------------------------------------------------------
    ! At
    !----------------------------------------------------------------------
    do t = (m + 1), n
       At(t) = log(y(t)) / log(mu(t))
    end do

    !----------------------------------------------------------------------
    ! dl/dmu
    !----------------------------------------------------------------------
    if (skip(1) == 0) then
       do t = (m + 1), n
          num(t) = nu(t) * (1 + log(argsD%par) * At(t)**nu(t))
          denum(t) = mu(t) * log(mu(t))
          dllmu(t) = - num(t) / denum(t)
       end do
    end if

    !----------------------------------------------------------------------
    ! dl/dnu
    !----------------------------------------------------------------------
    if (skip(2) == 0) then
       do t = (m + 1), n
          dllnu(t) = 1 / nu(t) + (1 + log(argsD%par) * At(t)**nu(t)) * log(At(t))
       end do
    end if
    return
  end subroutine dllk_uw

  subroutine Ed2llk_uw(argsD, m, n, mu, nu, skip, E)
    !**********************************************************************
    !  Computes
    !      E(dl2/dmu2), E(dl2/dmUdnu) and E(dl2/dnu2)
    !  for Unit-Weibull regression models.
    !**********************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer, intent(in) :: m, n, skip(2)
    real(dp), intent(in) :: mu(n), nu(n)
    real(dp), intent(out) :: E(n, 3)
    real(dp) :: cons0, cons1, cons2, lr
    integer  :: t

    argsD%dummy = .true.
    E = 0.d0
    if (sum(skip) == 2) return

    lr = log(-log(argsD%par))
    cons0 = em + lr
    cons1 = cons0 - 1
    cons2 = lr * (lr + 2 * em - 2) + pi**2 / 6 + (em - 2) * em

    do t = (m + 1), n
       if (skip(1) == 0) E(t, 1) = (nu(t) / (mu(t) * log(mu(t))))**2
       if (sum(skip) == 0) E(t, 2) = cons1 / (mu(t) * log(mu(t)))
       if (skip(2) == 0) E(t, 3) = (1.d0 + cons2) / nu(t)**2
    end do
    return
  end subroutine Ed2llk_uw

  !*****************************************************************************
  !
  ! Matsuoka's Distribution
  !
  !*****************************************************************************
  function rmatsu(argsD, npar, par) result(y)
    !*********************************************************************
    ! Generates Matsuoka's random variates using the relation
    !     X ~ M(lambda) <==>   - ln(X) ~ Gamma(3 / 2, 1 / lambda)
    ! with E(-ln(X)) = 3 / (2 * lambda), so that scale = 1 / lambda
    !
    ! par = mu,  mu = (1 / (1 + lambda))^(2 / 3) so that 0 < mu < 1
    !
    ! Uses rgamma() from R to generate random deviates from the gamma
    ! distribution whose density is
    !      1/((scale**shape)*Gamma(shape))*X**(shape - 1)*Exp(-X/scale)
    ! Arguments
    !   scale --> Location parameter of Gamma distribution  scale > 0 )
    !   shape --> Shape parameter of Gamma distribution ( shape > 0 )
    !
    ! Implemented by Taiane S. Prass - PPGEst/UFRGS
    ! June, 2024
    !*********************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer, intent(in) :: npar   ! dummy argument, must be 1
    real(dp), target, intent(in) :: par(npar)
    real(dp) :: y, lambda, gammarnd
    argsD%dummy = .true.
    lambda = par(1)**(2.d0 / 3.d0)
    lambda = lambda / (1.d0 - lambda)
    call rndstart()
    y = gammarnd(3.d0 / 2.d0, 1.d0 / lambda)
    y = exp(-y)
    call rndend()
    return
  end function rmatsu

  function dmatsu(y, npar, par, give_log) result(fn)
    !*********************************************************************
    ! Computes the density of Matsuoka's distribution.
    !
    ! par = mu,  mu = (1 / (1 + lambda))^(2 / 3) so that 0 < mu < 1
    !
    ! Implemented by Taiane S. Prass - PPGEst/UFRGS
    ! June, 2024
    !*********************************************************************
    implicit none
    integer,  intent(in) :: npar
    real(dp), intent(in) :: y
    real(dp), target, intent(in) :: par(npar)
    integer,  intent(in) :: give_log
    real(dp) :: fn, ly, m
    ly = log(y)
    m = par(1)**(2.d0 / 3.d0)
    fn = log(2.d0) - 0.5d0 * (log(pi) - log(-ly)) + log(par(1))
    fn = fn - 1.5d0 * log(1.d0 - m) + ly * (2 * m - 1) / (1 - m)
    if (give_log == 0) fn = exp(fn)
    return
  end function dmatsu

  function llk_matsu(argsD, m, n, y, mu, nu) result(sll)
    !***********************************************************************
    ! Computes sll Matsuoka regression models
    !***********************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in) :: m, n
    real(dp), intent(in) :: y(n), mu(n), nu(n) ! nu = dummy argument
    real(dp) :: sll
    integer  :: t
    argsD%dummy = .true.
    sll = nu(1)*0.d0 ! to avoid compiler warning
    do t = (m + 1), n
       sll = sll + dmatsu(y(t), 1, [mu(t)], 1)
    end do
    return
  end function llk_matsu

  subroutine dllk_matsu(argsD, m, n, y, mu, nu, skip, dllmu, dllnu)
    !*************************************************************************
    !  Computes dl/dmu and dl/dnu for Matsuoka regression models
    !  argsD contains rho, which is a constant in the model.
    !*************************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in) :: m, n, skip(2)
    real(dp), intent(in) :: y(n), mu(n), nu(n)
    real(dp), intent(out) :: dllmu(min(n, n * (1 - skip(1)) + 1))
    real(dp), intent(out) :: dllnu(min(n, n * (1 - skip(2)) + 1)) ! dummy
    real(dp) :: m23, mm
    integer  :: t

    argsD%dummy = .true.
    dllmu = 0.d0
    dllnu = nu(1)*0.d0 ! to avoid compiler warning
    if (skip(1) == 1) return

    !----------------------------------------------------------------------
    ! dl/dmu
    !----------------------------------------------------------------------
    do t = (m + 1), n
       m23 = mu(t)**(2.d0 / 3.d0)
       mm = 3.d0 * (1.d0 - m23)
       dllmu(t) = (2.d0 * m23 * log(y(t)) + mm) / mm**2
       dllmu(t) = 3.d0 * dllmu(t) / mu(t)
    end do
    return
  end subroutine dllk_matsu

  subroutine Ed2llk_matsu(argsD, m, n, mu, nu, skip, E)
    !**********************************************************************
    !   Computes
    !      E(dl2/dmu2), E(dl2/dmUdnu) and E(dl2/dnu2)
    !   for Matsuoka regression models.
    !**********************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in)  :: m, n, skip(2)
    real(dp), intent(in)  :: mu(n), nu(n)
    real(dp), intent(out) :: E(n, 3)
    real(dp) :: m2, m3
    integer  :: t

    argsD%dummy = .true.
    E = nu(1)*0.d0 ! to avoid compiler warning
    if (skip(1) == 1) return

    do t = (m + 1), n
       m2 = mu(t)**2
       m3 = m2**(1.d0 / 3.d0)
       E(t, 1) = (4.d0 - 10.d0 * m3) / (3.d0 * (1 - m3)**2 * m2)
    end do
    return
  end subroutine Ed2llk_matsu

  !*****************************************************************************
  !
  ! Unit Lindley Distribution
  !
  !*****************************************************************************
  function rul(argsD, npar, par) result(y)
    !**************************************************************************
    ! Generates Unit Lindley random variates using the relation
    !   1) X ~ Lindley(lambda) <=>  X/(1 + X) ~ UL(lambda)
    !   2) Lindley distribution is a special mixture of Exponential(lambda) and
    !      gamma(2, lambda), where lambda = rate and shape = 1 / lambda
    !
    ! par = mu,  mu = 1/(1 + lambda), so that 0 < mu < 1
    !
    ! Uses rgamma() from R to generate random deviates from the gamma
    ! distribution whose density is
    !   1/((scale**shape)*Gamma(shape))*X**(shape - 1)*Exp(-X / scale)
    ! Arguments
    !   scale --> Location parameter of Gamma distribution  scale > 0 )
    !   shape --> Shape parameter of Gamma distribution ( shape > 0 )
    !
    ! Implemented by Taiane S. Prass - PPGEst/UFRGS
    ! June, 2024
    !**************************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer, intent(in) :: npar   ! dummy argument, must be 1
    real(dp), target, intent(in) :: par(npar)
    real(dp) :: y, lambda, gammarnd, p, u
    argsD%dummy = .true.
    lambda = (1 - par(1)) / par(1)
    ! Generating a Lindlay random variate
    u = rng_uniform()
    p = lambda / (lambda + 1)
    call rndstart()
    if (u > p) then
       y = gammarnd(2.d0, 1.d0 / lambda)  ! Generate from Gamma(2, lambda)
    else
       y = gammarnd(1.d0, 1.d0 / lambda)  ! Generate from Exponential(lambda)
    end if
    call rndend()
    y = y / (y + 1)  ! Transform to Unit Lindley
    return
  end function rul

  function dul(y, npar, par, give_log) result(fn)
    !**********************************************************************
    ! Computes the density of the Unit Lindley distribution.
    ! par = mu,  mu = 1 / (1 + lambda) so that 0 < mu < 1
    !
    ! Implemented by Taiane S. Prass - PPGEst/UFRGS
    ! June, 2024
    !**********************************************************************
    implicit none
    integer,  intent(in) :: npar
    real(dp), intent(in) :: y
    real(dp), target, intent(in) :: par(npar)
    integer,  intent(in) :: give_log
    real(dp) :: fn
    fn = 2.d0 * log(1.d0 - par(1)) - log(par(1)) - 3.d0 * log(1.d0 - y)
    fn = fn + y / par(1) * (par(1) - 1.d0) / (1.d0 - y)
    if (give_log == 0) fn = exp(fn)
    return
  end function dul

  function llk_ul(argsD, m, n, y, mu, nu) result(sll)
    !***********************************************************************
    ! Computes sll Unit Lindley regression models
    !***********************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in) :: m, n
    real(dp), intent(in) :: y(n), mu(n), nu(n) ! nu = dummy argument
    real(dp) :: sll
    integer  :: t
    argsD%dummy = .true.
    sll = nu(1)*0.d0 ! to avoid compiler warnings
    do t = (m + 1), n
       sll = sll + dul(y(t), 1, [mu(t)], 1)
    end do
    return
  end function llk_ul

  subroutine dllk_ul(argsD, m, n, y, mu, nu, skip, dllmu, dllnu)
    !*************************************************************************
    !  Computes dl/dmu and dl/dnu for Unit Lindley regression models
    !  argsD contains rho, which is a constant in the model.
    !*************************************************************************
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in)  :: m, n, skip(2)
    real(dp), intent(in)  :: y(n), mu(n), nu(n)
    real(dp), intent(out) :: dllmu(min(n, n * (1 - skip(1)) + 1))
    real(dp), intent(out) :: dllnu(min(n, n * (1 - skip(2)) + 1)) ! dummy
    integer :: t

    argsD%dummy = .true.
    dllmu = 0.d0
    dllnu = nu(1)*0.d0 ! to avoid compiler warning
    if (skip(1) == 1) return

    !----------------------------------------------------------------------
    ! dl/dmu
    !----------------------------------------------------------------------
    do t = (m + 1), n
       dllmu(t) = y(t) / (1.d0 - y(t))
       dllmu(t) = dllmu(t) / mu(t)**2 - 1.d0 / mu(t) - 2.d0 / (1.d0 - mu(t))
    end do
    return
  end subroutine dllk_ul

  subroutine Ed2llk_ul(argsD, m, n, mu, nu, skip, E)
    !--------------------------------------------------------
    !   Computes
    !        E(dl2/dmu2), E(dl2/dmUdnu) and E(dl2/dnu2)
    !   for Unit Lindley regression models.
    !--------------------------------------------------------
    implicit none
    class(argsDist), intent(inout) :: argsD
    integer,  intent(in)  :: m, n, skip(2)
    real(dp), intent(in)  :: mu(n), nu(n)
    real(dp), intent(out) :: E(n, 3)
    real(dp) :: m1
    integer  :: t

    argsD%dummy = .true.
    E = nu(1)*0.d0 ! to avoid compiler warning
    if (skip(1) == 1) return

    do t = (m + 1), n
       m1 = (1.d0 - mu(t))**2
       E(t, 1) = (m1 - 2.d0) / mu(t)**2
       E(t, 1) = E(t, 1) / m1
    end do
    return
  end subroutine Ed2llk_ul

end module distrib
