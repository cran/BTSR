module base
  !-------------------------------------------------------------------------------------------------
  !
  !  This module contains the base subroutines
  !
  !    Y ~ f(y | mu, nu),  gnu = g(nu)
  !    g11(mu) = a1 + X1 * b1 + sum(phi1 * [g12(y) - X1 * b1 * I(X1ar)]) + sum(ck1 * r)      (part1)
  !    g21(gnu) = a2 + X2 * b2 + sum(phi2 * [g22(gnu) - X2 * b2 * I(X2ar)]) + sum(ck2 * r2)  (part2)
  !
  !-------------------------------------------------------------------------------------------------
  use alloc       ! allocation subroutines
  use lbfgsb      ! L-BFGS-B algorithm
  use Nelder_mead ! Nelder-Mead algorithm
  implicit none

  !-------------------------------------------------------------
  !
  !  Arguments to pass to optimization function
  !
  !-------------------------------------------------------------
  interface set_link
     module procedure set_link_to_link
     module procedure set_link_to_model
  end interface set_link

  interface vc_f
     module procedure vc_f1
     module procedure vc_f2
  end interface vc_f

  interface start_par
     module procedure start_par1
     module procedure start_par2
     module procedure start_par12
  end interface start_par

  type :: optimizer
     logical :: dummy = .true. ! to avoid compiler warnings
     procedure(optim_generic), pointer :: optim => null()
  end type optimizer

  abstract interface
     subroutine optim_generic(optim, loglik, model, npar, par, nbd, bounds, sll, score, cf1, nc2, cf2, &
          neval, conv)
       import :: dp, optimFunc, argsModel, optimizer
       implicit none
       class(optimizer), intent(inout) :: optim
       type(optimFunc), target, intent(inout) :: loglik
       type(argsModel), intent(inout) :: model
       integer,  target, intent(in)   :: cf1(3)
       integer,  intent(in)    :: npar, nc2
       integer,  intent(in)    :: nbd(npar)
       integer,  intent(out)   :: neval, conv
       real(dp), intent(in)    :: bounds(npar, 2)
       real(dp), target, intent(in) :: cf2(nc2)
       real(dp), intent(out)   :: sll, score(max(1, npar*cf1(3)))
       real(dp), intent(inout) :: par(npar)
     end subroutine optim_generic
  end interface


contains

  !-------------------------------------------------------------------------------------------------
  !
  !                       Subroutines related to the link function
  !
  !-------------------------------------------------------------------------------------------------
  subroutine make_shift(x, xlower, xupper, part, rev, iprint)
    !---------------------------------------------------------------------------------------
    ! Check if the generated value is out-of-bounds
    ! and replace by default values, if necessary
    !---------------------------------------------------------------------------------------
    ! February, 2025: added this function
    implicit none
    real(dp), intent(inout) :: x
    real(dp), intent(in)    :: xlower, xupper
    integer, intent(in)     :: part
    integer, intent(inout)  :: rev
    logical, intent(in)     :: iprint
    logical :: lfin, ufin
    real(dp), parameter :: hg = Huge(1.d0)

    rev = 0
    if(x > xlower .and. x < xupper) return

    lfin = is_finite(xlower)
    ufin = is_finite(xupper)

    rev = 1
    if(iprint) then
       call labelpr('----------------------------------------------------', -1)
       call labelpr(' Warning:', -1)
       if(part == 1) call labelpr('  - mu(t) out of bounds.', -1)
       if(part == 2) call labelpr('  - nu(t) or g(nu(t)) out of bounds.', -1)
    end if

    ! use huge instead of the infinite bounds
    ! so the code does not break when fitting a model
    if (x <= xlower) then
       if (.not. lfin) then  ! lower bound = -Inf
          if(part == 1) rev = 11
          if(part == 2) rev = 12
          x = -hg + epsmch
       else
          x = xlower + epsmch
       end if
       if(iprint) call labelpr('  - Replacing it by the default upper bound', -1)
    else
       if(.not. ufin) then ! upper bound = Inf
          if(part == 1) rev = 21
          if(part == 2) rev = 22
          x = hg - epsmch
       else
          x = xupper - epsmch
       end if
       if(iprint) call labelpr('  - Replacing it by the default lower bound', -1)
    end if
    if (iprint) call labelpr('----------------------------------------------------', -1)
    return
  end subroutine make_shift

  function g_err1(y, mu, g11y, eta, code) result(rt)
    !---------------------------------------------------------------------------------------
    ! computes the error term rt
    !---------------------------------------------------------------------------------------
    ! Input
    ! y   : the value of yt
    ! mu  : the value of mut
    ! g11y: the value of g11(yt)
    ! eta : the value of g11(mut)
    ! code: the error scale.
    !         0: error = y - mu,
    !         1: error = g12(y) - g12(mu)
    !
    ! Output
    ! rt: the error term
    !---------------------------------------------------------------------------------------
    ! February, 2025: added this function
    implicit none
    real(dp), intent(in) :: y, mu, g11y, eta
    integer,  intent(in) :: code
    real(dp) :: rt

    rt = 0.d0
    select case(code)
    case(0)
       rt = y - mu
    case(1)
       rt = g11y - eta
    end select
    return
  end function g_err1

  subroutine check_update(args1, args2)
    !---------------------------------------------------------------------------------------
    ! compares two link arguments to check if the links
    ! are the same
    !---------------------------------------------------------------------------------------
    ! Input/Output
    !   args1, args2: argslink type variable
    !---------------------------------------------------------------------------------------
    ! February, 2025: added this function
    type(argslink), intent(inout) :: args1, args2
    logical :: update

    args1%update = .false.
    args2%update = .false.
    update = args1%link == args2%link
    update = update .and. (args1%lower == args2%lower)
    update = update .and. (args1%upper == args2%upper)
    update = update .and. (all(args1%par == args2%par))
    if(update) then
       args1%update = .true.
       args2%update = .true.
    end if
    return
  end subroutine check_update

  function linkfun(x, args) result(lk)
    !---------------------------------------------------------------------------------------
    ! link function
    !            g:(a, b)  -> (c, d)
    !
    ! For g11 and g21 we have (c, d) = (-Inf, Inf)
    ! For g12, g22, g and h the interval (c, d) can be different
    !---------------------------------------------------------------------------------------
    ! input
    !   x   : the argument of the function g
    !   args: an argslink type variable with extra arguments for the link
    !
    ! output
    !   lk: the value of g(x)
    !---------------------------------------------------------------------------------------
    ! Last revision: February, 2024
    ! - updated the description
    ! - added the links used for gnut in Beta and Gamma regression
    implicit none
    type(argslink), intent(in) :: args
    real(dp), intent(in) :: x
    real(dp) :: lk
    real(dp) :: a, b, ctt, power

    a = args%lower ! lower bound for x
    b = args%upper ! upper bound for x
    ctt = args%par(1)
    power = args%par(2)
    lk = 0.d0

    select case(args%link)
    case(0)
       ! linear (identity, if ctt = 1 and power = 1)
       ! otherwhise, polynomial
       lk  = x
       if (abs(power - 1.d0) > epsmch) lk = lk**power
       if (abs(ctt - 1.d0) > epsmch) lk  = ctt * lk
    case(1)
       ! logit
       lk = log((x - a) / (b - x))
    case(2)
       ! log
       lk = log(x - a)
    case(3)
       ! loglog
       lk = log(-log((x - a) / (b - a)))
    case(4)
       ! cloglog
       lk = log(-log(1.d0 - (x - a) / (b - a)))
    case(5)
       if (abs(ctt) > epsmch) then
          ! link for Beta regression
          lk = 1.d0 / (1.d0 + x)
       else
          ! link for Gamma regression
          lk = 1.d0 / x
       end if
       ! general formula
       if (abs(power - 1.d0) > epsmch) lk = lk**power
    end select
    return
  end function linkfun

  function linkinv(x, args) result(inv)
    !---------------------------------------------------------------------------------------
    ! Inverse link function
    !              g^{-1}:(c, d)  -> (a, b)
    !
    ! For g11 and g21 we have (c, d) = (-Inf, Inf)
    ! For g12, g22, g and h the interval (c, d) can be different
    !---------------------------------------------------------------------------------------
    ! input
    !   x   : the argument of the function g^{-1}
    !   args: an argslink type variable with extra arguments for the link
    !
    ! output
    !   inv: the value of g^{-1}(x)
    !---------------------------------------------------------------------------------------
    ! Last revision: July, 2023
    implicit none
    type(argslink), target, intent(in) :: args
    real(dp), intent(in) :: x
    real(dp) :: inv
    real(dp), pointer :: a, b, ctt, power

    a => args%lower ! lower bound for x
    b => args%upper ! upper bound for x
    ctt => args%par(1)
    power => args%par(2)
    inv = 0.d0

    select case(args%link)
    case(0)
       ! linear (identity, if ctt = 1 and power = 1)
       inv  = x
       if (abs(ctt - 1.d0) > epsmch) inv  = inv / ctt
       if (abs(power - 1.d0) > epsmch) inv = (inv)**(1 / power)
    case(1)
       ! logit
       inv = b / (1.d0 + exp(-x))
       if (abs(a) > epsmch) inv = inv + a / (exp(x) + 1.d0)
    case(2)
       ! log
       inv = a + exp(x)
    case(3)
       ! loglog
       inv =  a + (b - a) * exp(-exp(x))
    case(4)
       ! cloglog
       inv = b - (b - a) * exp(-exp(x))
    case(5)
       ! general case
       if (abs(power - 1.d0) > epsmch) then
          inv = x**(1.d0 / power)
       else
          inv = x
       end if
       ! Gamma Regression
       inv = 1.d0 / inv
       ! Beta Regression
       if (abs(ctt) > epsmch) inv = inv - 1.d0
    end select
    return
  end function linkinv

  function diflink(x, args) result(dl)
    !---------------------------------------------------------------------------------------
    ! link function derivative
    !                   g'(x) = dg(x) / dx
    !---------------------------------------------------------------------------------------
    ! input
    !   x   : the argument of the function g'
    !   args: an argslink type variable with extra arguments for the link
    !
    ! output
    !   dl: the value of g'(x)
    !---------------------------------------------------------------------------------------
    ! Last revision: July, 2023
    implicit none
    type(argslink), target, intent(in) :: args
    real(dp) :: x
    real(dp) :: dl
    real(dp), pointer :: a, b, ctt, power

    a => args%lower ! lower bound for x
    b => args%upper ! upper bound for x
    ctt => args%par(1)
    power => args%par(2)
    dl = 0.d0

    select case(args%link)
    case(0)
       ! linear (identity, if ctt = 1 and power = 1)
       dl = 1.d0
       if (abs(power - 1.d0) > epsmch) then
          dl = (power - 1) * x**(power - 1)
       end if
       if (abs(ctt - 1.d0) > epsmch) dl  = dl * ctt
    case(1)
       ! logit
       dl = (b - a) / ((a + b) * x - x**2 - a * b)
    case(2)
       ! log
       dl = 1.d0 / (x - a)
    case(3)
       ! loglog
       dl = 1.d0 / ((x - a) * log((x - a) / (b - a)))
    case(4)
       ! cloglog
       dl = 1.d0 / ((x - b) * log((b - x) / (b - a)))
    case(5)
       ! Gamma Regression
       dl = x
       ! Beta regression
       if (abs(ctt) > epsmch) dl = dl + 1.d0
       if (abs(power - 1.d0) > epsmch) then
          dl = -power / (dl)**(power + 1.d0)
       else
          dl =  -1.d0 / dl**2
       end if
    end select
    return
  end function diflink

  !-------------------------------------------------------------------------------------------------
  !
  !                         Subroutines related to polinomials
  !
  !-------------------------------------------------------------------------------------------------
  subroutine pi_f(d, inf, pik)
    !---------------------------------------------------------------------------------------
    !  coeficients of
    !           pik(z) = (1 - z)^{ - d}
    !                  = sum_{k = 0}^infty pi_k * z^k
    !---------------------------------------------------------------------------------------
    ! Input
    !   d  : fractional order
    !   inf: trunction point
    !
    ! Output
    !   pik: the coefficients pi_k, for k = 0, ..., Inf
    !---------------------------------------------------------------------------------------
    ! Last revision: February, 2025
    !  - changed d == 0 for abs(d) < epsilon(1.d0)
    implicit none
    integer,  intent(in)  :: inf
    real(dp), intent(in)  :: d
    real(dp), intent(out) :: pik(0:inf)
    integer :: j
    pik = 0.d0
    pik(0) = 1.d0
    ! check if d is zero
    if (abs(d) < epsmch) return
    do j = 1, inf
       pik(j) = pik(j - 1) * (j - 1 + d) / j
    end do
    return
  end subroutine pi_f

  subroutine vc_f1(d, q, theta, inf, ck)
    !---------------------------------------------------------------------------------------
    !  coeficients of
    !               ck(z) = theta(z) * (1 - z)^{ - d}
    !                     = sum_{j = 0}^infty c_k * z^k
    !  where
    !                c_k = sum_{i=0}^{min{k, q}} theta_i * pik_(k - i)
    !---------------------------------------------------------------------------------------
    ! Input
    !   d    : fractional order
    !   q    : degree of theta(z)
    !   theta: (1, theta_1, ..., theta_q), the coefficients of theta(z)
    !   inf  : truncation point
    !
    ! Output
    !   ck: the coefficients c_k, k = 1, ..., inf
    !---------------------------------------------------------------------------------------
    ! Last revision: February, 2025
    !  - changed d == 0 for abs(d) < epsilon(1.d0)
    implicit none
    integer,   intent(in)  :: inf, q
    real (dp), intent(in)  :: d
    real (dp), intent(in)  :: theta(0:q)
    real (dp), intent(out) :: ck(0:inf)
    real (dp) :: pik(0:inf)
    integer:: i, k

    ck = 0.d0
    ! check if d is zero
    if (abs(d) < epsmch) then
       ! to avoid extra calculation
       ck(0:q) = theta
       return
    end if

    call pi_f(d, inf, pik)
    do k = 0, q
       do i = 0, k
          ck(k) = ck(k) + pik(k - i) * theta(i)
       end do
    end do
    do k = (q + 1), inf
       do i = 0, q
          ck(k) = ck(k) + pik(k - i) * theta(i)
       end do
    end do
    return
  end subroutine vc_f1

  subroutine vc_f2(model, vc, part)
    !---------------------------------------------------------------------------------------
    ! Computes the coefficients of (1 - L)^{-d}
    !---------------------------------------------------------------------------------------
    ! Input
    !   model: the model's configurations
    !   part: which part of the model needs initialization
    !
    ! Output
    !  vc: the coefficients of (1 - L)^{ - d}
    !---------------------------------------------------------------------------------------
    ! March, 2025: added this subroutine
    implicit none
    type(argsModel), intent(in) :: model
    integer,  intent(in)  :: part
    real(dp), intent(out) :: vc(0:model%pt(part)%inf)
    real(dp) :: theta(0:model%pt(part)%ma%length)

    theta(0) = 1.d0
    if (model%pt(part)%ma%length > 0) theta(1:model%pt(part)%ma%length) = model%pt(part)%ma%cf
    call vc_f1(model%pt(part)%d%cf(1), model%pt(part)%ma%length, theta, model%pt(part)%inf, vc)

    return
  end subroutine vc_f2


  !-------------------------------------------------------------------------------------------------
  !
  !                Subroutines used to calculate the conditional time series
  !                associated to the model
  !
  !-------------------------------------------------------------------------------------------------
  subroutine mu_calc(n, yt, g2start, g11y, g12y, nreg, xreg, xstart, mu, eta, error, &
       alpha, beta, p, phi, xregar, inf, vc, m, argsLg)
    !---------------------------------------------------------------------------------------
    ! Calculates conditional time series mu using a ARMA-type structure
    !---------------------------------------------------------------------------------------
    ! Input
    !   n      : sample size
    !   yt     : observed time series
    !   g2start: starting value for g(yt) (will raplace g(yt), for t < 1)
    !   g11y   : g_11(yt). Transformed time series to be used to calculate the residuals
    !   g12y   : g_12(yt). Transformed time series to be used in the AR part of the model
    !   nreg   : number of regressors
    !   xreg   : matrix of regressors
    !   xstart : starting value for xreg (will raplace xreg, for t < 1)
    !   alpha  : intercept
    !   beta   : coefficients for the regressors
    !   p      : order of AR polinomial
    !   phi    : the AR coefficients
    !   xregar : 0 = xreg is included only in the intercept
    !            1 = xreg is also included in the AR part.
    !   inf    : trunction point for infinite sums
    !   vc     : coeficients ck
    !   m      : starting point for the sum of the log-likelihood
    !   argsLg : arguments for link functions g_11, g_12 and g_13 (the error term)
    !
    ! Output
    !   mu   : conditional time series mut
    !   eta  : conditional time series g_11(mut)
    !   error: the error term rt. Depends on escale.
    !---------------------------------------------------------------------------------------
    ! Last revision: February, 2025
    !  - removed escale
    !  - renamed gy to g12y and added g11y
    !  - now errors are computed using a function
    implicit none
    type(argslink), intent(in) :: argsLg(3)
    integer,  intent(in)  :: m, n, nreg, p, inf, xregar
    real(dp), intent(out) :: mu(n), eta(n), error(n)
    real(dp), intent(in)  :: yt(n), g2start, g11y(n), g12y(n)
    real(dp), intent(in)  :: xreg(n, max(1, nreg)), xstart(max(1, nreg))
    real(dp), intent(in)  :: alpha, phi(max(1, p)), beta(max(1, nreg)), vc(0:inf)
    integer  :: t, j, rev
    real(dp) :: a, b
    real(dp) :: xb, g12ytemp

    ! To avoid the bounds:
    a = argsLg(1)%lower
    b = argsLg(1)%upper

    ! initialization
    error = 0.d0
    eta = 0.d0
    g12ytemp = g2start ! g_12(y)
    xb = 0.d0

    ! starting values for xreg
    if (p > 0 .and. xregar == 1 .and. nreg > 0) xb = sum(xstart * beta)

    do t = (m + 1), n
       !------------------------------------
       ! eta(t) = alpha + x * b + AR + MA
       !------------------------------------
       eta(t) = alpha
       if (nreg > 0) eta(t) = eta(t) + sum(xreg(t, :) * beta)

       do j = 1, p
          if (t - j > 0) then
             ! update g12ytemp if necessary
             g12ytemp = g12y(t - j)
             ! check if xreg is used in AR recursion and update xb
             if (xregar == 1 .and. nreg > 0) xb = sum(xreg(t - j, :) * beta)
          end if
          eta(t) = eta(t) + (g12ytemp - xb) * phi(j)
       end do

       do j = 1, min(t - 1, inf)
          ! sum(c(j) * r(t - j)), only exists if d > 0 or q > 0
          ! (if inf < 1 does not enters the loop)
          ! assuming error(t) = 0, for t < 0
          eta(t) =  eta(t) + vc(j) * error(t - j)
       end do

       ! compute mut and check the bounds
       mu(t) = linkinv(eta(t), argsLg(1))

       ! to avoid the bounds:
       call make_shift(mu(t), a, b, 1, rev, .false.)
       if(rev > 0) eta(t) = linkfun(mu(t), argsLg(1))

       ! calculating rt
       error(t) = g_err1(yt(t), mu(t), g11y(t), eta(t), argsLg(3)%link)
    end do

    return
  end subroutine mu_calc

  subroutine nu_calc(n, error, e2, e2start, nreg, xreg, xstart, wt, gwt, g22gw, g2start, eta, &
       alpha, beta, p, phi, xregar, inf, vc, m, argsLg)
    !---------------------------------------------------------------------------------------
    ! Calculates the conditional time series gnut and nut using a ARMA-type structure
    !---------------------------------------------------------------------------------------
    ! Input:
    !   n      : sample size
    !   error  : the error term rt corresponding to part 1 of the model
    !   e2start: starting values for the square of rt
    !   nreg   : number of regressors
    !   xreg   : matrix of regressors
    !   xstart : starting value for xreg (will raplace xreg, for t < 1)
    !   g2start: starting values for g12(gwt)
    !   alpha  : intercept
    !   beta   : coefficients for the regressors
    !   p      : order of AR polinomial
    !   phi    : the AR coefficients
    !   xregar : 0 = xreg is included only in the intercept
    !            1 = xreg is also included in the AR part.
    !   inf    : trunction point for infinite sums
    !   vc     : coeficients ck
    !   m      : starting point for the sum of the log-likelihood
    !   argsLg : arguments for link functions g2, g_21, g_22 and g_23 (error term)
    !
    ! Output
    !   e2   : the error term for g23(e1)
    !   wt   : conditional time series nut = g^{-1}(gnut)
    !   gwt  : conditional time series gnut = g(wt)
    !   eta  : conditional time series g_21(gnut)
    !   g22gw: conditional time series g_22(gnut)
    !
    !---------------------------------------------------------------------------------------
    ! July, 2023
    !   - added vt
    !   - argsL now has size 3 and corresponds to g, g_21 and g_22.
    !     Code was corrected accordingly
    !
    ! Last revision: February, 2025
    !   - replace r^2 by g(r) in part 2.
    !   - renamed: argsL, vt, ut and ustart (argsLg, wt, gwt and g2start)
    implicit none
    type(argslink), intent(in) :: argsLg(4)
    integer,  intent(in)  :: m, n, nreg, p, inf, xregar
    real(dp), intent(in)  :: error(n), e2start, g2start
    real(dp), intent(in)  :: xreg(n, max(1, nreg)), xstart(max(1, nreg))
    real(dp), intent(in)  :: phi(max(1, p)), beta(max(1, nreg))
    real(dp), intent(in)  :: alpha, vc(0:inf)
    real(dp), intent(out) :: wt(n), gwt(n), g22gw(n), eta(n), e2(n)
    real(dp) :: error2(-inf:n), g22, xb
    integer  :: t, j, rev

    ! initializing
    XB = 0.d0
    wt = 0.d0        ! g^{-1}(gnu)
    eta = 0.d0       ! g_21(gnut)
    g22gw = 0.d0     ! g_22(gnu)
    g22 = g2start    ! temporary g_22(gnu)
    e2 = 0.d0        ! g_23(rt)
    error2 = e2start ! g_23(rt) complete vector

    if (p > 0 .and. xregar == 1 .and. nreg > 0) xb = sum(xstart * beta)

    do t = (m + 1), n
       !------------------------------------
       ! eta(t) = alpha + x * b + AR + MA
       !------------------------------------
       eta(t) = alpha
       if (nreg > 0) eta(t) = eta(t) + sum(xreg(t, :) * beta)

       do j = 1, p
          if (t - j > 0) then
             ! update g22(gnu)
             g22 = g22gw(t-j)
             ! check if xreg is used in AR recursion
             if (xregar == 1 .and. nreg > 0) xb = sum(xreg(t - j, :) * beta)
          end if
          eta(t) = eta(t) + (g22 - xb) * phi(j)
       end do

       do j = 1, inf
          ! sum(c(j) * r2(t - j)), only exists if d > 0 or q > 0 (if inf < 1 does not enters the loop)
          eta(t) =  eta(t) + vc(j) * error2(t - j)
       end do

       gwt(t) = linkinv(eta(t), argsLg(2)) ! g(wt) = g21^{-1}(eta)
       wt(t) = linkinv(gwt(t), argsLg(1))  ! wt = g^{-1}(g(wt))
       call make_shift(wt(t), argsLg(1)%lower, argsLg(1)%upper, 2, rev, .false.)
       if(rev > 0) then ! update
          gwt(t) = linkfun(wt(t), argsLg(1))  ! g(wt) = g2(wt)
          eta(t) = linkfun(gwt(t), argsLg(2)) ! eta(t) = g21(gwt)
       end if

       if (argsLg(3)%update) then
          g22gw(t) = linkfun(gwt(t), argsLg(3))
       else
          g22gw(t) = eta(t)
       end if

       error2(t) = linkfun(error(t), argsLg(4)) ! g23(rt)
       e2(t) = error2(t)
    end do
    return
  end subroutine nu_calc

  subroutine mu_forecast(model, nnew, xhat, forecast)
    !---------------------------------------------------------------------------------------
    !  Calculates the out-of-sample predicted values for mut, t > n
    !---------------------------------------------------------------------------------------
    ! Input
    !   model: a argsModel object with the model's configurations
    !   nnew : number of new predictions
    !   xhat : new values for the regressors Xt
    !
    ! Output
    !   forcast: the predicted values for mut and eta1t
    !---------------------------------------------------------------------------------------
    ! July, 2023
    !  - renamed yhat to muhat and gyhat to etahat
    !
    ! Last revision: March 2025
    !  - replaced muhat and g12hat with forecast.
    !  - forecast now returns g12yt, muhat and eta1hat
    !  - vc is now an internal variable
    implicit none
    type(argsModel), intent(in) :: model
    integer,  intent(in)  :: nnew
    real(dp), intent(in)  :: xhat(nnew, max(1, model%pt(1)%beta%length))
    real(dp), target, intent(out) :: forecast(nnew, 2)
    real(dp), pointer     :: muhat(:), etahat(:)
    real(dp) :: XB((model%n + 1 - model%pt(1)%ar%length):(model%n + nnew))
    real(dp) :: g12y((model%n + 1 - model%pt(1)%ar%length):(model%n + nnew))
    real(dp)  :: vc(0:model%pt(1)%inf), g12hat(nnew)
    integer  :: t, i

    forecast = 0.d0
    muhat   => forecast(:,1)
    etahat  => forecast(:,2)
    g12hat = 0.d0

    ! compute the coefficients ck
    call vc_f(model, vc, 1)

    ! xreg * beta, t = n + 1 - p, ..., n + 1, ..., n + nnew
    XB = 0.d0
    if (model%pt(1)%beta%length > 0) then
       do t = (model%n + 1 - model%pt(1)%ar%length), model%n
          XB(t) = sum(model%cts(1)%xreg(t, :) * model%pt(1)%beta%cf)
       end do
       do t = 1, nnew
          XB(model%n + t) = sum(xhat(t, :) * model%pt(1)%beta%cf)
       end do
    end if

    ! g_12(y)
    if (model%pt(1)%ar%length > 0) then
       i = model%n + 1 - model%pt(1)%ar%length
       g12y(i:model%n) = model%cts(1)%gi2(i:model%n)
    end if

    do t = 1, nnew
       ! etahat = a + x * b
       etahat(t) = model%pt(1)%alpha%cf(1) + XB(model%n + t)
       ! etahat = a + x * b + AR
       do i = 1, model%pt(1)%ar%length   ! if p < 1 ignores the loop
          etahat(t) = etahat(t) + model%pt(1)%ar%cf(i) * g12y(model%n + t - i)
          if (model%cts(1)%xregar == 1) then
             etahat(t) = etahat(t) - model%pt(1)%ar%cf(i) * XB(model%n + t - i)
          end if
       end do
       ! etahat = a + x * b + AR + MA
       do i = t, min(model%n + t - 1, model%pt(1)%inf)
          etahat(t) = etahat(t) + vc(i) * model%cts(1)%et(model%n + t - i)
       end do
       ! muhat = g_11^{-1}(etahat)
       muhat(t) = linkinv(etahat(t), model%pt(1)%linkg(2))
       ! g_12(y) to be used in the AR recursion
       if (model%pt(1)%linkg(3)%update) then
          ! get the prediction for g12(yt)
          g12y(model%n + t) = linkfun(muhat(t), model%pt(1)%linkg(3))
       else
          g12y(model%n + t) = etahat(t)
       end if
       g12hat(t) = g12y(model%n + t)
    end do
    return
  end subroutine mu_forecast

  subroutine nu_forecast(model, nnew, xhat, forecast)
    !---------------------------------------------------------------------------------------
    !  Predicted values for nut and gnut, t > n
    !---------------------------------------------------------------------------------------
    ! Input
    !   model: a argsModel object with the model's configurations
    !   nnew : number of new predictions
    !   xhat : new values for the regressors Xt
    !
    ! Output
    !   forecast : the predicted values for
    !              wt, g(wt) and g21(g(wt))
    !---------------------------------------------------------------------------------------
    ! Implemented by Taiane S. Prass - PPGEst / UFRGS
    ! July, 2023 - Porto Alegre
    !
    ! Last revision: February, 2025
    !   - replace r^2 by g(r) in part 2.
    !   - replaced nuhat and vhat with forecast
    !   - vc is now an internal variable
    implicit none
    type(argsModel), intent(in) :: model
    integer,  intent(in)  :: nnew
    real(dp), intent(in)  :: xhat(nnew, max(1, model%pt(2)%beta%length))
    real(dp), target, intent(out) :: forecast(nnew,3)
    real(dp), pointer     :: what(:), gwhat(:), etahat(:)
    real(dp) :: XB((model%n + 1 - model%pt(2)%ar%length):(model%n + nnew))
    real(dp) :: g22((model%n + 1 - model%pt(2)%ar%length):(model%n + nnew))
    real(dp) :: error2(-model%pt(2)%inf:(model%n + nnew))
    real(dp) :: vc(0:model%pt(2)%inf), g22hat(nnew)
    integer  :: t, i

    forecast = 0.d0
    what   => forecast(:,1)
    gwhat  => forecast(:,2)
    etahat => forecast(:,3)
    g22hat = 0.d0

    if (model%pt(2)%inf > 0) then
       ! initializing
       error2 = model%cts(2)%estart ! g_{23}(rt)
       ! 1 <= t < n: use g23(e)
       error2(1:model%n) = model%cts(2)%et
       ! t > n: use the sample mean for g23(e)
       error2(model%n:(model%n + nnew)) = sum(error2(1:model%n)) / dble(model%n)

       ! compute the coefficients ck
       call vc_f(model, vc, 2)
    end if

    ! xreg * beta, t = n + 1 - p, ..., n + 1, ..., n + nnew
    XB = 0.d0
    if (model%pt(2)%beta%length > 0) then
       do t = (model%n + 1 - model%pt(2)%ar%length), model%n
          XB(t) = sum(model%cts(2)%xreg(t, :) * model%pt(2)%beta%cf)
       end do
       do t = 1, nnew
          XB(model%n + t) = sum(xhat(t, :) * model%pt(2)%beta%cf)
       end do
    end if

    ! g_22(sigma^2), t = n - p, ..., n + nnew
    if (model%pt(2)%linkg(3)%update) then
       i = model%n + 1 - model%pt(2)%ar%length
       g22(i:model%n) = model%cts(2)%gi2(i:model%n)
    end if

    do t = 1, nnew
       ! etahat = a + x * b
       etahat(t) = model%pt(2)%alpha%cf(1) + XB(model%n + t)
       ! etahat = a + x * b + AR
       do i = 1, model%pt(2)%ar%length   ! if p < 1 ignores the loop
          etahat(t) = etahat(t) + model%pt(2)%ar%cf(i) * g22(model%n + t - i)
          if (model%cts(2)%xregar == 1) then
             etahat(t) = etahat(t) - model%pt(2)%ar%cf(i) * XB(model%n + t - i)
          end if
       end do
       ! etahat = a + x * b + AR + MA
       do i = 1, model%pt(2)%inf
          etahat(t) = etahat(t) + vc(i) * error2(model%n + t - i)
       end do
       ! gwhat = g_21^{-1}(etahat)
       gwhat(t) = linkinv(etahat(t), model%pt(2)%linkg(2))
       ! nuhat = g2^{-1}(vthat)
       what(t) = linkinv(gwhat(t), model%pt(2)%linkg(1))
       ! g_22(gnu) to be used in the AR recursion
       if (model%pt(2)%linkg(3)%update) then
          g22(model%n + t) = linkfun(gwhat(t), model%pt(2)%linkg(3))
       else
          g22(model%n + t) = etahat(t)
       end if
       g22hat(t) = g22(model%n + t)
    end do
    return
  end subroutine nu_forecast

  subroutine forecast_model(model, nnew, xnew1, xnew2, forecast)
    !---------------------------------------------------------------------------------------
    !  Calculates the out-of-sample predicted values for mut and nut
    !---------------------------------------------------------------------------------------
    ! Input
    !   model: a argsModel object with the model's configurations
    !   nnew : number of new predictions
    !   xnew1 : new values for the regressors Xt1
    !   xnew2 : new values for the regressors Xt2
    !
    ! Output
    !   forcast: the predicted values
    !---------------------------------------------------------------------------------------
    ! March 2025
    !  - added this subroutine
    implicit none
    type(argsModel), intent(in) :: model
    integer,  intent(in)  :: nnew
    real(dp), intent(in)  :: xnew1(nnew, max(1, model%pt(1)%beta%length))
    real(dp), intent(in)  :: xnew2(nnew, max(1, model%pt(2)%beta%length))
    real(dp), intent(out) :: forecast(nnew, 5)
    forecast = 0.d0
    call mu_forecast(model, nnew, xnew1, forecast(:, 1:2))
    call nu_forecast(model, nnew, xnew2, forecast(:, 3:5))
    return
  end subroutine forecast_model

  !-------------------------------------------------------------------------------------------------
  !
  !           Subroutines to calculate the generic log-likelihood
  !
  !-------------------------------------------------------------------------------------------------
  subroutine loglik_generic(model, npar, par, sll, U)
    !---------------------------------------------------------------------------------------
    !   Calculates the log-likelihood and score vector using a user defined density function
    !---------------------------------------------------------------------------------------
    ! Input
    !   npar     : the number of unknown parameters in the model
    !   par      : a vector of size npar with the non-fixed parameters of the model
    !
    ! Input/Output
    !   model: input, the model's configurations.
    !          output, updated values of the conditional time series and the corresponding score
    !          vector (if required) and related matrices and vectors
    !
    ! Output
    !   sll: the sum of the log-likelihood
    !   U  : the score vector corresponding tho par
    !---------------------------------------------------------------------------------------
    ! Last revision: July, 2023
    !   - added vc2 and renamed vc to vc1
    !   - added parameter initialization for part 2 of the model
    !   - llk_dist now has nu(n) as argument instead of nu
    !   - added the call to nu_calc
    implicit none
    type(argsModel), intent(inout) :: model
    integer,  intent(in)  :: npar
    real(dp), intent(in)  :: par(npar)
    real(dp), intent(out) :: sll, U(npar)
    real(dp) :: vc1(0:model%pt(1)%inf)
    real(dp) :: vc2(0:model%pt(2)%inf)

    ! Initializing parameters
    call start_par(par, model, vc1, 1)
    call start_par(par, model, vc2, 2)

    ! Calculating recursions that do not depend on the distribution
    call mu_calc(model%n, model%y, model%cts(1)%g2start, model%cts(1)%gi1, model%cts(1)%gi2, &
         model%cts(1)%nreg, model%cts(1)%xreg, model%cts(1)%xstart, model%cts(1)%w, &
         model%cts(1)%eta, model%cts(1)%et, model%pt(1)%alpha%cf(1), model%pt(1)%beta%cf,  &
         model%pt(1)%ar%length, model%pt(1)%ar%cf, model%cts(1)%xregar, model%pt(1)%inf, &
         vc1, model%m, model%pt(1)%linkg(2:4))

    call nu_calc(model%n, model%cts(1)%et, model%cts(2)%et, model%cts(2)%estart, model%cts(2)%nreg, &
         model%cts(2)%xreg, model%cts(2)%xstart, model%cts(2)%w, model%cts(2)%gw, model%cts(2)%gi2,&
         model%cts(2)%g2start, model%cts(2)%eta, model%pt(2)%alpha%cf(1), &
         model%pt(2)%beta%cf, model%pt(2)%ar%length, model%pt(2)%ar%cf, model%cts(2)%xregar, &
         model%pt(2)%inf, vc2, model%m, model%pt(2)%linkg)

    !-----------------------------------------
    ! log - likelihood for Generic model
    !-----------------------------------------
    sll = model%argsD%llk_dist(model%m, model%n, model%y, model%cts(1)%w, model%cts(2)%w)
    sll = -sll

    U = 0.d0
    if (model%sco == 0) return

    !-----------------------------------------
    !    Score vector
    !-----------------------------------------
    call U_generic(model, vc1, vc2, U)
    U = -U
    return
  end subroutine loglik_generic

  !-------------------------------------------------------------------------------------------------
  !
  !   Subroutines to calculate  and organize the derivatives deta_t/dgamma_j
  !   and the associated matrices Drho and Dlambda
  !
  !-------------------------------------------------------------------------------------------------
  subroutine deta1_drho(model, vc)
    !---------------------------------------------------------------------------------------
    ! Calculates the derivative deta_{1t}/drho, m <= t <= n
    !---------------------------------------------------------------------------------------
    ! Input
    !   vc   : the coefficients in the MA expansion in part1
    !
    ! Input/Output
    !   model: a argsModel object with the model's configurations.
    !          on exit, updates SI, an argsSI object containing the values
    !          of deta_{1t}/drho_j
    !---------------------------------------------------------------------------------------
    ! Last revision: March, 2025
    !  - added pointers
    implicit none
    type(argsModel), target, intent(in) :: model
    real(dp), intent(in)  :: vc(0:model%pt(1)%inf)
    integer  :: i, j, jj,  k, t, n1, n2
    real(dp) :: soma, drdeta, a
    real(dp) :: pik(0:min(model%pt(1)%inf, model%n - 1))
    real(dp) :: difdg(0:min(model%pt(1)%inf, model%n - 1))
    integer, pointer :: es, m, n, inf
    integer, pointer :: fita, fitb, fitar, fitma, fitd
    real(dp), pointer :: dalpha1(:), dd1(:)
    real(dp), pointer :: dbeta1(:,:), dphi1(:,:), dtheta1(:,:)

    !--------------------------------------------------------------
    ! If part 1 is fixed, there is no calculation to be done here
    !--------------------------------------------------------------
    if (model%pt(1)%npar == 0) return

    n   => model%n
    m   => model%m
    inf => model%pt(1)%inf
    fita  => model%pt(1)%alpha%fit
    fitb  => model%pt(1)%beta%fit
    fitar => model%pt(1)%ar%fit
    fitma => model%pt(1)%ma%fit
    fitd  => model%pt(1)%d%fit

    ! error scale
    es => model%pt(1)%linkg(4)%link
    ! dr/deta1 =  -1, if es = 0;
    !            -T1, if es = 1.
    drdeta = -1.d0

    !--------------------------------------------------------------
    ! compute vc and start the recurrency
    !--------------------------------------------------------------
    if (fitma + fitd > 0) then
       call pi_f(model%pt(1)%d%cf(1), min(inf, n - 1), pik)
       if(fitd > 0) then
          a = digamma(model%pt(1)%d%cf(1))
          do i = 0, min(inf, n-1)
             difdg(i) = digamma(model%pt(1)%d%cf(1) + i) - a
          end do
       end if
    end if

    n1 = 1
    n2 = fita
    ! alpha
    if (n2 >= n1) then
       !
       ! diff(t) = 1 + sum(ck * dr/deta1 * diff(t - k))
       !
       dalpha1 => model%SI%deta(1, 1)%dalpha(:, 1)
       dalpha1 = 0.d0
       do t = (m  + 1), n
          dalpha1(t) = 1.d0
          do j = 1, min(t - 1, inf)
             if (es == 0) drdeta = -model%SI%T(1)%z(t - j)
             dalpha1(t) = dalpha1(t) + vc(j) * drdeta * dalpha1(t - j)
          end do
       end do
    end if

    n1 = n2 + 1
    n2 = n2 + fitb
    ! betas
    if (n2 >= n1) then
       !
       ! diff(t) = x(t, j) - sum(phi(k) * x(t - k, j)) + sum(c(k) * dr/deta1 * diff(t - k))
       !
       dbeta1 => model%SI%deta(1, 1)%dbeta
       dbeta1 = 0.d0
       do j = 1, model%pt(1)%beta%fit
          jj = model%pt(1)%beta%lags(j)  ! j-th non-fixed lag
          do t = (m  + 1), n
             dbeta1(t, j) = model%cts(1)%xreg(t, jj)
             ! if p = 0, then ar = 0
             if (model%cts(1)%xregar == 1) then
                do k = 1, min(t - 1, model%pt(1)%ar%length)
                   dbeta1(t, j) = dbeta1(t, j) - model%pt(1)%ar%cf(k) * model%cts(1)%xreg(t - k, jj)
                end do
             end if
             do k = 1, min(t - 1, inf)
                if (es == 0) drdeta = -model%SI%T(1)%z(t - k)
                dbeta1(t, j) = dbeta1(t, j) + vc(k) * drdeta * dbeta1(t - k, jj)
             end do
          end do
       end do
    end if

    n1 = n2 + 1
    n2 = n2 + fitar
    ! phis
    if (n2 >= n1) then
       !
       ! diff(t) = g12(y(t - j)) - I_X*X(t - j) * b -  sum(ck * dr/deta1 * diff(t - k))
       !
       dphi1 => model%SI%deta(1, 1)%dphi
       dphi1 = 0.d0
       do j = 1,  model%pt(1)%ar%fit
          jj = model%pt(1)%ar%lags(j)  ! j-th non-fixed lag
          do t =  max((jj + 1), (m  + 1)), n   ! so that t - jj > 0
             dphi1(t, j) =  model%cts(1)%gi2(t - jj)
             if (model%cts(1)%xregar == 1) then
                do i = 1, model%cts(1)%nreg  ! if nreg = 0 does not enter the loop
                   dphi1(t, j) = dphi1(t, j) - model%cts(1)%xreg(t - jj, i) * model%pt(1)%beta%cf(i)
                end do
             end if
             do k = 1, min(t - 1, inf)
                if (es == 0) drdeta = -model%SI%T(1)%z(t - k)
                dphi1(t, j) = dphi1(t, j) + vc(k) * drdeta * dphi1(t - k, j)
             end do
          end do
       end do
    end if

    n1 = n2 + 1
    n2 = n2 + fitma
    ! thetas
    if (n2 >= n1) then
       !
       ! diff(t) = sum(pik(k - s) * r(t - k), k = s, ..., infty) +  sum(ck * dr/deta1 * diff(t - k))
       !
       dtheta1 => model%SI%deta(1, 1)%dtheta
       dtheta1 = 0.d0
       do j = 1, model%pt(1)%ma%fit  !if q > 0, then vc(k) is different from zero for some k =1, ..., inf
          jj = model%pt(1)%ma%lags(j)  ! j-th non-fixed lag
          do t = (m  + 1), n
             do k = 1, min(t - 1, (jj - 1))
                if (es == 0) drdeta = -model%SI%T(1)%z(t - k)
                dtheta1(t, j) = dtheta1(t, j) + vc(k) * drdeta * dtheta1(t - k, j)
             end do
             do k = jj, min(t - 1, inf)
                if (es == 0) drdeta = -model%SI%T(1)%z(t - k)
                dtheta1(t, j) =  dtheta1(t, j) +  pik(k - jj) * model%cts(1)%et(t - k) &
                     + vc(k) * drdeta * dtheta1(t - k, j)
             end do
          end do
       end do
    end if

    n1 = n2 + 1
    n2 = n2 + fitd
    ! d
    if (n2 >= n1) then
       !
       ! diff(t) = sum(r(t - k) * s(k)) +  sum(ck * dr/deta1 * diff(t - k)),
       !  where s(k) = sum(theta(i) * pi(k - i) * digamas(i))
       !
       dd1 => model%SI%deta(1, 1)%dd(:, 1)
       dd1 = 0.d0
       do t = (m  + 1), n
          do k = 1, min(t - 1, inf)
             ! i = 0
             soma = pik(k) * difdg(k)
             do i = 1, min(k, model%pt(1)%ma%length) ! theta(0) = 1. If q = 0 does not enter the loop
                soma = soma + model%pt(1)%ma%cf(i) * pik(k - i) * difdg(k - i)
             end do
             if (es == 0) drdeta = -model%SI%T(1)%z(t - k)
             dd1(t) = dd1(t) + model%cts(1)%et(t - k) * soma + vc(k) * drdeta * dd1(t - k)
          end do
       end do
    end if
    return
  end subroutine deta1_drho

  subroutine deta2_drho(model, vc)
    !---------------------------------------------------------------------------------------
    ! Calculates the derivative deta_{2t}/drho, m <= t <= n
    !---------------------------------------------------------------------------------------
    ! Input
    !   vc   : the coefficients in the MA expansion in part2
    !
    ! input/Output
    !   model: a argsModel object with the model's configurations.
    !          on exit, updates SI, an argsSI object containing the values
    !          of deta_{2t}/drho_j
    !---------------------------------------------------------------------------------------
    ! Last revision: March, 2025
    !  - added pointers
    !  - replaced 2 * r by g23'(r)
    implicit none
    type(argsModel), target, intent(in) :: model
    real(dp), intent(in)  :: vc(0:model%pt(2)%inf)
    real(dp) :: drdeta
    real(dp) :: dgdr(model%n), dgg(model%n)
    integer  :: j, k, t, n1, n2
    integer, pointer :: es, m, n, inf
    integer, pointer :: fita, fitb, fitar, fitma, fitd
    real(dp), pointer :: dalpha1(:), dd1(:)
    real(dp), pointer :: dbeta1(:,:), dphi1(:,:), dtheta1(:,:)

    ! If any part is fixed, there is no calculation to be done here
    if (product(model%pt%npar) == 0) return

    !
    !  diff2(t) = sum(phi(k) * g'_22(g2(nu(t - k))) / g'_21(g2(nu(t - k))) * diff2(t - k))
    !             + sum(c(k) * g'_23(r(t - k)) * dr/deta1 * diff1(t - k))
    !

    n   => model%n
    m   => model%m
    inf => model%pt(2)%inf
    fita  => model%pt(1)%alpha%fit
    fitb  => model%pt(1)%beta%fit
    fitar => model%pt(1)%ar%fit
    fitma => model%pt(1)%ma%fit
    fitd  => model%pt(1)%d%fit

    ! error scale
    es => model%pt(1)%linkg(4)%link
    ! dr/deta =  -1, if es = 0;
    !           -T1, if es = 1.
    drdeta = -1.d0

    ! g'_23(rt)
    do t = 1, n
       dgdr(t) = diflink(model%cts(1)%et(t), model%pt(2)%linkg(4))
    end do

    ! g22'(g(nu)) / g21'(g(nu))
    if (model%pt(2)%linkg(3)%update) then
       do t = 1,n
          dgg(t) = diflink(model%cts(2)%gw(t), model%pt(2)%linkg(3))
          dgg(t) = dgg(t) / diflink(model%cts(2)%gw(t), model%pt(2)%linkg(2))
       end do
    else
       dgg = 1.d0
    end if

    n1 = 1
    n2 = fita
    ! alpha
    if (n2 >= n1) then
       dalpha1 => model%SI%deta(2, 1)%dalpha(:, 1)
       dalpha1 = 0.d0
       do t = (m  + 1), n
          do k = 1, min(t - 1, model%pt(2)%ar%length)
             dalpha1(t) = dalpha1(t) + model%pt(2)%ar%cf(k) * dgg(t - k) * dalpha1(t - k)
          end do
          do k = 1, min(t - 1, inf)
             if (es == 0) drdeta = -model%SI%T(1)%z(t - k)
             dalpha1(t) = dalpha1(t) + vc(k) * dgdr(t - k) * drdeta * model%SI%deta(1, 1)%dalpha(t - k, 1)
          end do
       end do
    end if

    n1 = n2 + 1
    n2 = n2 + fitb
    ! betas
    if (n2 >= n1) then
       dbeta1 => model%SI%deta(2, 1)%dbeta
       dbeta1 = 0.d0
       do j = 1, model%pt(1)%beta%fit
          do t = (m  + 1), n
             do k = 1, min(t - 1, model%pt(2)%ar%length)
                dbeta1(t, j) = dbeta1(t, j) + model%pt(2)%ar%cf(k) * dgg(t - k) * dbeta1(t - k, j)
             end do
             do k = 1, min(t - 1, inf)
                if (es == 0) drdeta = -model%SI%T(1)%z(t - k)
                dbeta1(t, j) = dbeta1(t, j) + vc(k) * dgdr(t - k) * drdeta * model%SI%deta(1, 1)%dbeta(t - k, j)
             end do
          end do
       end do
    end if

    n1 = n2 + 1
    n2 = n2 + fitar
    ! phis
    if (n2 >= n1) then
       dphi1 => model%SI%deta(2, 1)%dphi
       dphi1 = 0.d0
       do j = 1, model%pt(1)%ar%fit
          do t = (m  + 1), n
             do k = 1, min(t - 1, model%pt(2)%ar%length)
                dphi1(t, j) = dphi1(t, j) + model%pt(2)%ar%cf(k) * dgg(t - k) * dphi1(t - k, j)
             end do
             do k = 1, min(t - 1, inf)
                if (es == 0) drdeta = -model%SI%T(1)%z(t - k)
                dphi1(t, j) = dphi1(t, j) + vc(k) * dgdr(t - k) * drdeta * model%SI%deta(1, 1)%dphi(t - k, j)
             end do
          end do
       end do
    end if

    n1 = n2 + 1
    n2 = n2 + fitma
    ! thetas
    if (n2 >= n1) then
       dtheta1 => model%SI%deta(2, 1)%dtheta
       dtheta1 = 0.d0
       do j = 1, model%pt(1)%ma%fit
          do t = (m  + 1), n
             do k = 1, min(t - 1, model%pt(2)%ar%length)
                dtheta1(t, j) = dtheta1(t, j) + model%pt(2)%ar%cf(k) * dgg(t - k) * dtheta1(t - k, j)
             end do
             do k = 1, min(t - 1, inf)
                if (es == 0) drdeta = -model%SI%T(1)%z(t - k)
                dtheta1(t, j) = dtheta1(t, j) + vc(k) * dgdr(t - k) * drdeta * model%SI%deta(1, 1)%dtheta(t - k, j)
             end do
          end do
       end do
    end if

    n1 = n2 + 1
    n2 = n2 +  fitd
    ! d
    if (n2 >= n1) then
       dd1 => model%SI%deta(2, 1)%dd(:, 1)
       dd1 =  0.d0
       do t = (m  + 1), n
          do k = 1, min(t - 1, model%pt(2)%ar%length)
             dd1(t) = dd1(t) + model%pt(2)%ar%cf(k) * dgg(t - k) * dd1(t - k)
          end do
          do k = 1, min(t - 1, inf)
             if (es == 0) drdeta = -model%SI%T(1)%z(t - k)
             dd1(t) = dd1(t) + vc(k) * dgdr(t - k) * drdeta * model%SI%deta(1, 1)%dd(t - k, 1)
          end do
       end do
    end if
    return
  end subroutine deta2_drho

  subroutine deta2_dlambda(model)
    !---------------------------------------------------------------------------------------
    ! Calculates the derivative deta_{2t}/dlambda, m <= t <= n
    !---------------------------------------------------------------------------------------
    ! Input/output
    !   model: a argsModel object with the model's configurations.
    !          on exit, updates SI, an argsSI object containing the values
    !          of deta_{2t}/dlambda_j
    !---------------------------------------------------------------------------------------
    ! Last revision: March, 2025
    !  - added pointers
    implicit none
    type(argsModel), target, intent(in) :: model
    real(dp) :: dgg(model%n)
    real(dp) :: pik(0:model%pt(2)%inf), difdg(0:model%pt(2)%inf)
    real(dp) :: error2(-model%pt(2)%inf:model%n), soma, a
    integer  :: i, j, jj,  k, t, n1, n2
    integer, pointer :: m, n, inf
    integer, pointer :: fita, fitb, fitar, fitma, fitd
    real(dp), pointer :: dalpha2(:), dd2(:)
    real(dp), pointer :: dbeta2(:,:), dphi2(:,:), dtheta2(:,:)

    !-----------------------------------------------------------------
    ! If part 2 is fixed, there is no calculation to be done here
    if (model%pt(2)%npar == 0) return
    !-----------------------------------------------------------------

    n   => model%n
    m   => model%m
    inf => model%pt(2)%inf
    fita  => model%pt(2)%alpha%fit
    fitb  => model%pt(2)%beta%fit
    fitar => model%pt(2)%ar%fit
    fitma => model%pt(2)%ma%fit
    fitd  => model%pt(2)%d%fit

    ! initializing
    error2 = model%cts(2)%estart ! g_{23}(rt)
    ! 1 <= t < n : use g23(rt)
    if (inf > 0) error2(1:n) = model%cts(2)%et(1:n)

    ! g22'(g(nu)) / g21'(g(nu))
    if (model%pt(2)%linkg(3)%update) then
       do t = 1, n
          dgg(t) = diflink(model%cts(2)%gw(t), model%pt(2)%linkg(3))
          dgg(t) = dgg(t) / diflink(model%cts(2)%gw(t), model%pt(2)%linkg(2))
       end do
    else
       dgg = 1.d0
    end if

    n1 = 1
    n2 = fita
    ! alpha
    if (n2 >= n1) then
       !
       ! diff(t) = 1 + sum(phi(k) * g'_22(g2(nu(t - k))) / g'_21(g2(nu(t - k))) * diff(t - k))
       !
       dalpha2 => model%SI%deta(2, 2)%dalpha(:, 1)
       dalpha2 = 0.d0
       do t = (m  + 1), n
          dalpha2(t) = 1.d0
          do k = 1, min(t - 1, model%pt(2)%ar%length)
             dalpha2(t) = dalpha2(t) + model%pt(2)%ar%cf(k) * dgg(t - k) * dalpha2(t - k)
          end do
       end do
    end if

    n1 = n2 + 1
    n2 = n2 + fitb
    ! betas
    if (n2 >= n1) then
       !
       ! diff(t, j) = x(t, j) + sum(phi(k) * g'_22(g2(nu(t - k))) / g'_21(g2(nu(t - k))) *  [diff(t - k, j) - I_X * x(t - k, j)])
       !
       dbeta2 => model%SI%deta(2, 2)%dbeta
       dbeta2 = 0.d0
       do j = 1, model%pt(2)%beta%fit
          jj = model%pt(2)%beta%lags(j)  ! j-th non-fixed lag
          do t = (m  + 1), n
             dbeta2(t, j) = model%cts(2)%xreg(t, jj)
             ! if p = 0, then ar = 0
             do k = 1, min(t - 1, model%pt(2)%ar%length)
                dbeta2(t, j) = dbeta2(t, j) + model%pt(2)%ar%cf(k) * dgg(t - k) * dbeta2(t - k, j)
                if (model%cts(2)%xregar == 1) dbeta2(t, j) = dbeta2(t, j) - &
                     model%pt(2)%ar%cf(k) * model%cts(2)%xreg(t - k, jj)
             end do
          end do
       end do
    end if

    n1 = n2 + 1
    n2 = n2 + fitar
    ! phis
    if (n2 >= n1) then
       !
       ! diff(t) = g(wt(t - j)) - I_X * X(t - j) * b +  sum(phi(k) * g'_22(g2(nu(t - k))) / g'_21(g2(nu(t - k))) *  diff(t - k, j))
       !
       dphi2 => model%SI%deta(2, 2)%dphi
       dphi2 = 0.d0
       do j = 1, model%pt(2)%ar%fit
          jj = model%pt(2)%ar%lags(j)    ! j-th non-fixed lag
          do t =  max((jj + 1), (m  + 1)), n   ! so that t - jj > 0
             dphi2(t, j) =  model%cts(2)%gi2(t - jj)
             if (model%cts(2)%xregar == 1) then
                do i = 1, model%cts(2)%nreg
                   dphi2(t, j) = dphi2(t, j) - model%cts(2)%xreg(t - jj, i) * model%pt(2)%beta%cf(i)
                end do
             end if
             do k = 1, min(t - 1, model%pt(2)%ar%length)
                dphi2(t, j) = dphi2(t, j) + model%pt(2)%ar%cf(k) * dgg(t - k) * dphi2(t - k, j)
             end do
          end do
       end do
    end if

    n1 = n2 + 1
    n2 = n2 + fitma
    ! thetas
    if (n2 >= n1) then
       !
       ! diff(t, j) = g_23(r(t - j)) +  sum(phi(k) * g'_22(g2(nu(t - k))) / g'_21(g2(nu(t - k))) *  diff(t - k, j))
       !
       dtheta2 => model%SI%deta(2, 2)%dtheta
       dtheta2 = 0.d0
       do j = 1,  model%pt(2)%ma%fit
          jj = model%pt(2)%ma%lags(j)  ! j-th non-fixed lag
          do t = (m  + 1), n
             dtheta2(t, j) =  error2(t - jj)
             do k = 1, min(t - 1, model%pt(2)%ar%length)
                dtheta2(t, j) = dtheta2(t, j) + model%pt(2)%ar%cf(k) * dgg(t - k) * dtheta2(t - k, j)
             end do
          end do
       end do
    end if

    n1 = n2 + 1
    n2 = n2 + fitd
    ! d
    if (n2 >= n1) then
       !
       ! diff(t) =  sum(phi(k) * g'_22(g2(nu(t - k))) / g'_21(g2(nu(t - k))) *  diff(t - k, j)) +  sum(g_23(r(t - k)) * s(k)),
       ! where s(k) = sum(theta(i) * pi(k - i) * digamas(i))
       !
       call pi_f(model%pt(2)%d%cf(1), inf, pik)
       a = digamma(model%pt(2)%d%cf(1))
       do i = 0, inf
          difdg(i) = digamma(model%pt(2)%d%cf(1) + i) - a
       end do

       dd2 => model%SI%deta(2, 2)%dd(:, 1)
       dd2 = 0.d0
       do t = (m  + 1), n
          do k = 1, min(t - 1, model%pt(2)%ar%length)
             dd2(t) = dd2(t) + model%pt(2)%ar%cf(k) * dgg(t - k) * dd2(t - k)
          end do
          do k = 1, inf
             soma = pik(k) * difdg(k)
             do i = 1, model%pt(2)%ma%length ! theta(0) = 1. If q = 0 does not enter the loop
                soma = soma + model%pt(2)%ma%cf(i) * pik(k - i) * difdg(k-i)
             end do
             dd2(t) = dd2(t) + error2(t - k) * soma
          end do
       end do
    end if

    return
  end subroutine deta2_dlambda

  subroutine fill_D(deta, fita, fitb, fitar, fitma, fitd, n, nd, D, zero)
    !---------------------------------------------------------------------------------------
    ! Creates a matrix with the derivatives deta1/drho, deta2/drho and deta2/dlambda
    ! Will be used to return values when calling from R or to calculate the matrix K.
    !---------------------------------------------------------------------------------------
    ! Input
    !   SI   : an argsSI object containing the values of deta_{*t}/dgamma_j
    !   fita : 1 = fit alpha, 0 = alpha is known
    !   fitb : number of unknown parameters in the beta vector
    !   fitar: number of unknown parameters in the phi (AR) vector
    !   fitma: number of unknown parameters in the theta (MA) vector
    !   fitd : 1 = fit d, 0 = d is known
    !   n    : sample size
    !   nd   : number of unknown parameters
    !   i1   : indicates which eta (1 = eta_{1t}, 2 = eta_{2t})
    !   i2   : indicates which parameter vector (1 = rho or 2 = lambda)
    !
    ! Output
    !   D: matriz of derivatives.
    !   zero : logical, indicates if the matrix is zero
    !---------------------------------------------------------------------------------------
    ! Last revision: July, 2023
    !  - added the variables i1 and i2 so that the same subroutine can be used
    !    for Drho, Dlambda and Mrho.
    implicit none
    type(deta_d), intent(in) :: deta
    integer,  intent(in)  :: n, nd
    integer,  intent(in)  :: fita, fitb, fitar, fitma, fitd
    real(dp), intent(out) :: D(n, nd)
    logical,  intent(out) :: zero
    integer :: n1, n2, nfill

    zero = .true.
    nfill = 0

    n1 = 1
    n2 = fita
    ! alpha
    if (n2 >= n1) then
       D(:, n1:n2) = deta%dalpha
       nfill = nfill + 1
    end if

    n1 = n2 + 1
    n2 = n2 + fitb
    ! beta
    if (n2 >= n1) then
       D(:, n1:n2) = deta%dbeta
       nfill = nfill + 1
    end if

    n1 = n2 + 1
    n2 = n2 + fitar
    ! phi
    if (n2 >= n1) then
       D(:, n1:n2) = deta%dphi
       nfill = nfill + 1
    end if

    n1 = n2 + 1
    n2 = n2 + fitma
    ! theta
    if (n2 >= n1) then
       D(:, n1:n2) = deta%dtheta
       nfill = nfill + 1
    end if

    n1 = n2 + 1
    n2 = n2 + fitd
    ! d
    if (n2 >= n1) then
       D(:, n1:n2) = deta%dd
       nfill = nfill + 1
    end if

    if(nfill > 0) zero = .false.
    return
  end subroutine fill_D

  !-------------------------------------------------------------------------------------------------
  !
  !           Subroutines to calculate extra matrices used to calcualte U an K
  !
  !-------------------------------------------------------------------------------------------------
  subroutine calc_Ts(argsL, m, n, mu, nu, gnut, T1, T2, skip)
    !---------------------------------------------------------------------------------------
    ! Calculates T1 = dmu/deta1 and T2 = dnu/deta2
    !---------------------------------------------------------------------------------------
    ! Input
    !   args  : an argslink type variable with extra arguments for the link
    !   m     : starting point for the sum of the log-likelihood
    !   n     : sample size
    !   mu    : conditional time series mut
    !   nu    : conditional time series nut
    !   gnut  : conditional time series g(nut)
    !   skip  : if skip(1) = 1, then mu is known and fixed in the model. 0, otherwise
    !           if skip(2) = 1, then nu is known and fixed in the model. 0, otherwise
    !
    ! Output
    !   T1: the derivatives dmu/deta1
    !   T2: the derivatives dnu/deta2
    !---------------------------------------------------------------------------------------
    ! July, 2023
    !  - renamed to calc_Ts and added calculations for T2
    !  - argsL now has size 3
    !
    ! Last revision: February, 2025
    !  - replaced skipmu and skipnu by the vector skip with size 2
    implicit none
    type(argsLink), intent(in) :: argsL(6)  ! g1, g11, g12, g2, g21, g22
    integer,  intent(in)  ::  n, m, skip(2)
    real(dp), intent(in)  ::  mu(n), nu(n), gnut(n)
    real(dp), intent(out) :: T1(min(n, n * (1 - skip(1)) + 1))
    real(dp), intent(out) :: T2(min(n, n * (1 - skip(2)) + 1))
    integer  :: i
    real(dp) :: a, b

    !-------------------------------------------------------------
    ! T1 = dmu/deta1  (skip if mu is fixed)
    ! T2 = dnu/deta2  (skip if nu is fixed)
    !--------------------------------------------------------------
    T1 = 0.d0
    T2 = 0.d0
    do i = (m  + 1), n
       if (skip(1) == 0) T1(i) = 1.d0 / diflink(mu(i), argsL(2)) ! g_11'(mu)
       if (skip(2) == 0) then
          a = diflink(nu(i), argsL(4))    ! g'(nu)
          b = diflink(gnut(i), argsL(5))  ! g_21'(gnut)
          T2(i) = 1.d0 / (a * b)
       end if
    end do
    return
  end subroutine calc_Ts

  subroutine calc_hs(argsD, m, n, y, mu, nu, skip, h1, h2)
    !---------------------------------------------------------------------------------------
    !  Calculates h1 = dl/dmu and h2 = dl/dnu
    !---------------------------------------------------------------------------------------
    ! Input
    !   argsD    : an argsDist variable with extra information to be passed to the density function
    !   m        : starting point for the sum of the log-likelihood
    !   n        : sample size
    !   y        : the observed time series
    !   mu       : conditional time series mut
    !   nu       : conditional time series nut
    !   skip     : if skip(1) = 1, then mut is known and fixed in the model. 0, otherwise
    !            : if skip(2) = 1, then nu is known and fixed in the model. 0, otherwise
    !
    ! Output
    !   h1: the derivatives dl/dmu
    !   h2: the derivatives dl/dnu
    !---------------------------------------------------------------------------------------
    ! July, 2023
    !  - renamed the surbroutine from Unuh_dist to calc_hs
    !  - renamed the arguments Unu and h to h2 and h1 and changed the order
    !  - removed npar, renamed fitnu to skipnu and added skipmu
    !  - nu now has size n.
    !  - changed the subroutine accordingly to return h1 and h2
    !
    ! Last revision: February, 2025
    !  - replaced skipmu and skipnu by the vector skip with size 2
    implicit none
    type(argsDist), intent(inout) :: argsD
    integer,  intent(in)  :: m, n, skip(2)
    real(dp), intent(out) :: h1(min(n, n * (1 - skip(1)) + 1))
    real(dp), intent(out) :: h2(min(n, n * (1 - skip(2)) + 1))
    real(dp), intent(in)  :: y(n), mu(n), nu(n)
    real(dp), allocatable :: dllmu(:), dllnu(:)

    !-------------------------------------------------------------
    ! h1 = dl/dmu  (skip if mu is fixed)
    ! h2 = dl/dnu  (skip if nu is fixed)
    !--------------------------------------------------------------
    h1 = 0.d0
    h2 = 0.d0
    if (product(skip) == 1) return

    call safe_allocate(dllmu, min(n, n * (1 - skip(1)) + 1))
    call safe_allocate(dllnu, min(n, n * (1 - skip(2)) + 1))

    call argsD%dllk_dist(m, n, y, mu, nu, skip, dllmu, dllnu)

    if (skip(1) == 0) h1 = dllmu
    if (skip(2) == 0) h2 = dllnu

    return
  end subroutine calc_hs

  !-------------------------------------------------------------------------------------------------
  !
  !       Subroutines used to calculate the product and the sum of matrices that
  !       appear in the score U and in the information matrix K
  !
  !-------------------------------------------------------------------------------------------------
  subroutine ATh(nra, nca, A, Th, P)
    !---------------------------------------------------------------------------------------
    ! Calculates the product of a matrix (transposed) and a vector.
    ! Used to calculate A'Th and M'Th
    !---------------------------------------------------------------------------------------
    ! Input
    !   nra: number of rows in A
    !   nca: number of columns in A
    !   A  : a matrix
    !   Th : a vector
    !
    ! Output
    !   P: the matrix A'Th
    !---------------------------------------------------------------------------------------
    ! Last revision: July, 2023
    implicit none
    integer,  intent(in)  :: nra, nca
    real(dp), intent(in)  ::  A(nra, nca), Th(nra)
    real(dp), intent(out) :: P(nca)
    P = matmul(transpose(A), Th)
    return
  end subroutine ATh

  subroutine calc_DTh(model, SI, part, Th)
    !---------------------------------------------------------------------------------------
    ! Calculates the product D'Th, where D is Drho or Dlambda and Th is T1 * h1 or T2 * h2
    !---------------------------------------------------------------------------------------
    ! Input
    !   model: a argsModel object with the model's configurations
    !   part : an integer indicating if the derivatives correspond to part 1 or 2 of the model
    !   Th   : a vector with the product of T by h
    !
    ! Input/Output:
    !   SI: input, values of deta1/drho, deta2/drho, deta2/dlambda
    !       output, the values of DTh saved in the argument U
    !---------------------------------------------------------------------------------------
    ! Last revision: July, 2023
    implicit none
    type(argsModel), intent(in) :: model
    type(argsSI), intent(inout) :: SI
    integer,  intent(in) :: part
    real(dp), intent(in) :: Th(model%n * (1 - model%pt(part)%skip) + model%pt(part)%skip)

    ! If "part" is fixed, there is no calculation to be done here
    if (model%pt(part)%npar == 0) return

    ! alpha
    if (model%pt(part)%alpha%fit == 1) then
       call ATh(model%n, 1, SI%deta(part, part)%dalpha, Th, SI%U(part)%Ualpha)
    end if

    ! betas
    if (model%pt(part)%beta%fit > 0) then
       call ATh(model%n, model%pt(part)%beta%fit, SI%deta(part, part)%dbeta, Th, SI%U(part)%Ubeta)
    end if

    ! phis
    if (model%pt(part)%ar%fit > 0) then
       call ATh(model%n, model%pt(part)%ar%fit, SI%deta(part, part)%dphi, Th, SI%U(part)%Uphi)
    end if

    ! thetas
    if (model%pt(part)%ma%fit > 0) then
       call ATh(model%n, model%pt(part)%ma%fit, SI%deta(part, part)%dtheta, Th, SI%U(part)%Utheta)
    end if

    ! d
    if (model%pt(part)%d%fit == 1) then
       call ATh(model%n, 1, SI%deta(part, part)%dd, Th, SI%U(part)%Ud)
    end if
    return
  end subroutine Calc_DTh

  subroutine AddM(model, SI, Th)
    !---------------------------------------------------------------------------------------
    ! Calculates the product Mrho' * Th (Mlambda = 0) and adds this quantity to D'Th
    ! previously calculated
    !---------------------------------------------------------------------------------------
    ! Input
    !   model: a argsModel object with the model's configurations
    !   Th   : a vector with the product of T by h
    !
    ! Input/Output:
    !   SI: input, values of deta1/drho, deta2/drho, deta2/dlambda
    !       output, the values of D'Th + Mrho' * T * h saved (updated) in the argument U
    !---------------------------------------------------------------------------------------
    ! Last revision: July, 2023
    implicit none
    type(argsModel), intent(in) :: model
    type(argsSI), intent(inout) :: SI
    real(dp), intent(in)  :: Th(model%n)
    real(dp), allocatable :: dummy(:)

    ! If part 1 is fixed, there is no calculation to be done here
    if (model%pt(1)%npar == 0) return

    ! alpha
    if (model%pt(1)%alpha%fit  == 1) then
       call safe_allocate(dummy, 1)
       call ATh(model%n, 1, SI%deta(2, 1)%dalpha, Th, dummy)
       SI%U(1)%Ualpha = SI%U(1)%Ualpha + dummy
    end if

    ! betas
    if (model%pt(1)%beta%fit > 0) then
       call safe_allocate(dummy, model%pt(1)%beta%fit)
       call ATh(model%n, model%pt(1)%beta%fit, SI%deta(2, 1)%dbeta, Th, dummy)
       SI%U(1)%Ubeta = SI%U(1)%Ubeta + dummy
    end if

    ! phis
    if (model%pt(1)%ar%fit > 0) then
       call safe_allocate(dummy, model%pt(1)%ar%fit)
       call ATh(model%n, model%pt(1)%ar%fit, SI%deta(2, 1)%dphi, Th, dummy)
       SI%U(1)%Uphi = SI%U(1)%Uphi  + dummy
    end if

    ! thetas
    if (model%pt(1)%ma%fit > 0) then
       call safe_allocate(dummy, model%pt(1)%ma%fit)
       call ATh(model%n, model%pt(1)%ma%fit, SI%deta(2, 1)%dtheta, Th, dummy)
       SI%U(1)%Utheta = SI%U(1)%Utheta + dummy
    end if

    ! d
    if (model%pt(1)%d%fit == 1) then
       call safe_allocate(dummy, 1)
       call ATh(model%n, 1, SI%deta(2, 1)%dd, Th, dummy)
       SI%U(1)%Ud  = SI%U(1)%Ud + dummy
    end if
    return
  end subroutine AddM

  !-------------------------------------------------------------------------------------------------
  !
  !  Subroutines used to calculate the score vector U(gamma).
  !
  !-------------------------------------------------------------------------------------------------
  subroutine fill_U(SI, fita, fitb, fitar, fitma, fitd, npar, U)
    !---------------------------------------------------------------------------------------
    ! Creates a vector with the values of the score vector U(gamma).
    ! Will be used to return values when calling from R.
    !---------------------------------------------------------------------------------------
    ! Input
    !   SI   : an argsSI object with previously calculated values of U
    !   fita : 1 = fit alpha, 0 = alpha is known
    !   fitb : number of unknown parameters in the beta vector
    !   fitar: number of unknown parameters in the phi (AR) vector
    !   fitma: number of unknown parameters in the theta (MA) vector
    !   fitd : 1 = fit d, 0 = d is known
    !   npar : number of unknown parameters
    !
    ! Output:
    !   U: a vector with the values of dl/dgamma_j
    !---------------------------------------------------------------------------------------
    ! Last revision: July, 2023
    !  - removed any reference to nu
    !  - added loop to cover part 1 and 2
    implicit none
    type(argsSI), intent(in) :: SI
    integer,  intent(in)  :: npar, fita(2), fitb(2), fitar(2)
    integer,  intent(in)  :: fitma(2), fitd(2)
    real(dp), intent(out) :: U(npar)
    integer :: n1, n2, part

    n2 = 0
    do part = 1, 2
       n1 = n2 + 1
       n2 = n2 + fita(part)
       ! alpha
       if (n2 >= n1) U(n1:n2) = SI%U(part)%Ualpha
       n1 = n2 + 1
       n2 = n2 + fitb(part)
       ! beta
       if (n2 >= n1) U(n1:n2) = SI%U(part)%Ubeta
       n1 = n2 + 1
       n2 = n2 + fitar(part)
       ! phi
       if (n2 >= n1) U(n1:n2) = SI%U(part)%Uphi
       n1 = n2 + 1
       n2 = n2 + fitma(part)
       ! theta
       if (n2 >= n1) U(n1:n2) = SI%U(part)%Utheta
       n1 = n2 + 1
       n2 = n2 + fitd(part)
       ! d
       if (n2 >= n1) U(n1:n2) = SI%U(part)%Ud
    end do
    return
  end subroutine fill_U

  subroutine U_generic(model, vc1, vc2, U)
    !---------------------------------------------------------------------------------------
    ! Calculates the score vector using a user defined density function
    !     Urho(gamma) = Drho' * T1 * h1 + Mrho' * T2 * h2 and
    !     Ulambda(gamma) = Dlambda' * T2 * h2
    ! with
    !      Drho = deta1/drho, Mrho = deta2/drho, Dlambda = deta2/dlambda
    !      T1 = dmu/deta1, T2 = dnu/deta2, h1 = dl/dmu, h2 = dl/dnu
    !---------------------------------------------------------------------------------------
    ! Input
    !   vc1, vc2 : the coefficients of the MA recursion for part 1 and 2
    !
    ! Input/Output
    !   model: input, the model's configurations
    !          output, updated values of U and related matrices and vectors.
    !
    ! Output
    !   U: the score vector
    !---------------------------------------------------------------------------------------
    ! Last revision: July, 2023
    !  - merged calc_Us1 and calc_Us2 with U_generic
    !  - added coded needed to calculate the quantities for part 2.
    implicit none
    type(argsModel), intent(inout) :: model
    type(argsLink) :: linkg(6)
    real(dp), intent(in)  :: vc1(0:model%pt(1)%inf), vc2(0:model%pt(2)%inf)
    real(dp), intent(out) :: U(sum(model%pt%npar))
    real(dp) :: T1h1(min(model%n,model% n * (1 - model%pt(1)%skip) + 1))
    real(dp) :: T2h2(min(model%n, model%n * (1 - model%pt(2)%skip) + 1))

    ! calculating h1 and h2
    call calc_hs(model%argsD, model%m, model%n, model%y, model%cts(1)%w, model%cts(2)%w, &
         model%pt(1:2)%skip, model%SI%h(1)%z, model%SI%h(2)%z)

    linkg(1:3) = model%pt(1)%linkg(1:3)
    linkg(4:6) = model%pt(2)%linkg(1:3)
    ! calculating T1 and T2
    call calc_Ts(linkg, model%m, model%n, model%cts(1)%w, model%cts(2)%w, &
         model%cts(2)%gw, model%SI%T(1)%z, model%SI%T(2)%z, model%pt(1:2)%skip)

    ! Drho - no calculations done if npar(1) = 0
    call deta1_drho(model, vc1)
    ! Mrho - no calculations done if npar(1) = 0 or npar(2) = 0
    call deta2_drho(model, vc2)
    ! Dlambda - no calculations done if npar(2) = 0
    call deta2_dlambda(model)

    ! T1 * h1
    T1h1 = 0.d0
    if (model%pt(1)%npar > 0) T1h1 = model%SI%T(1)%z * model%SI%h(1)%z
    ! T2 * h2
    T2h2 = 0.d0
    if (model%pt(2)%npar > 0) T2h2 = model%SI%T(2)%z * model%SI%h(2)%z

    ! if part 1 is known, skip Drho * T1 * h1 + Mrho * T2 * h2
    if (model%pt(1)%npar > 0) then
       ! Drho * T1 * h1
       call calc_DTh(model, model%SI, 1, T1h1)
       ! Drho * T1 * h1 + Mrho * T2 * h2
       if (model%pt(2)%npar > 0) call AddM(model, model%SI, T2h2)
    end if

    ! if part 2 is known, skip Dlambda * T2 * h2
    if (model%pt(2)%npar > 0) then
       ! Dlambda * T2 * h2
       call calc_DTh(model, model%SI, 2, T2h2)
    end if

    ! Fill the score vector
    call fill_U(model%SI, model%pt(1:2)%alpha%fit, model%pt(1:2)%beta%fit, model%pt(1:2)%ar%fit, &
         model%pt(1:2)%ma%fit, model%pt(1:2)%d%fit, sum(model%pt(1:2)%npar), U)
    return
  end subroutine U_generic

  !-------------------------------------------------------------------------------------------------
  !
  !  Subroutines to calculate the information matrix K
  !
  !-------------------------------------------------------------------------------------------------
  subroutine calc_K(n, nr, nl, T1, T2, E, Dr, Mr, Dl, zero, K)
    !---------------------------------------------------------------------------------------
    !  Calculates the information matrix K with blocks
    !      Krr = Dr' * T1 * Em * T1 * Dr + Mr' * T2 * Emn * T1 * Dr +
    !            Dr' * T1 * Emn * T2 * Mr + Mr' * T2 * En * T2 * Mr
    !      Krl = Klr' = Dr' * T1 * Emn * T2 * Dl + Mr' * T2 * En * T2 * Dl
    !      Kll = Dl' * T2 * En * T2 * Dl
    !
    ! Input
    !  n   : sample size
    !  nr  : number of unknown parameters in part 1
    !  nl  : number of unknown parameters in part 2
    !  T1  : vector with the derivatives dmu/deta1
    !  T2  : vector with the derivatives dnu/deta2
    !  E   : matrix with minus the expected values of the log-likelihood's second derivatives
    !        E_mu = -E(d2l/dmu2), E_(mu, nu) = -E(d2l/dmUdnu) and E_nu = -E(d2l/dnu2)
    !  Dr  : matrix with the derivatives deta1/drho
    !  Mr  : matrix with the derivatives deta2/drho
    !  Dl  : matrix with the derivatives deta2/dlambda
    !  zero: indicates if Dr, Mr and Dl are zero to avoid extra calculations
    !
    ! Output
    !   K: the information matrix
    !---------------------------------------------------------------------------------------
    ! Last revision: July, 2023
    !  - removed the subroutine calc_K1 and renamed calc_K2 to calc_K
    implicit none
    integer,  intent(in)  :: n, nr, nl
    real(dp), intent(in)  :: T1(min(n, n * nr + 1)), T2(min(n, n * nl + 1)), E(n, 3)
    real(dp), intent(in)  :: Dr(n, max(1, nr)), Mr(n, max(1, nr)), Dl(n, max(1, nl))
    real(dp), intent(out) :: K(max(1, nr + nl), max(1, nr + nl))
    logical,  intent(in)  :: zero(3)
    integer :: j, i

    ! initialization
    K = 0.d0

    !--------------------------------------------------------------------------
    ! K_{rho, rho}
    !   = Dr' * T1 * Em * T1 * Dr
    !     + Mr' * T2 * Emn * T1 * Dr
    !     + (Dr' * T1 * Emn + Mr' * T2 * En) * T2 * Mr
    !---------------------------------------------------------------------------
    ! if mu is known,
    !  - K_{rho, rho} does not exist
    !  - nr = 0 => Dr = 0 (no need to check if zero(1) = .true.)
    !
    ! mu is not fixed
    !  - nr > 0 => Dr not zero (no need to check if zero(1) = .true.)
    !  - Mr depends on nu => must check if zero(2) = .true.
    !---------------------------------------------------------------------------
    if (nr > 0) then
       do i = 1, nr
          do j = 1, i
             K(i, j) = sum(Dr(:, j) * E(:, 1) * T1**2 * Dr(:, i))
             K(j, i) = K(i, j)
          end do
       end do
       if (.not. zero(2)) then
          do i = 1, nr
             do j = 1, i
                K(i, j) = K(i, j) + sum(Mr(:, j) * T2 * E(:, 2) *  T1 * Dr(:, i) + &
                     (Dr(:, j) *  T1 * E(:, 2) + Mr(:, j) * T2 * E(:, 3)) * T2 * Mr(:, i))
                K(j, i) = K(i, j)
             end do
          end do
       end if
    end if

    !--------------------------------------------------------------------------
    ! K_{rho, lambda} = K_{lambda, rho}'
    !     = Dr' * T1 * Emn * T2 * Dl + Mr' * T2 * En * T2 * Dl
    !---------------------------------------------------------------------------
    ! if mu is known or nu is known
    !  - K_{rho, lambda} does not exist
    !  - nr = 0 => Dr = 0 and Mr = 0
    !      (no need to check if zero(1) = .true. and zero(2) = .true.)
    !  - nl = 0 => Dl = 0 (no need to check if zero(3) = .true.)
    !
    ! if nu is not fixed
    !  - Mr = 0 if nu does not have the MA (p = 0 and d = 0) recursion
    !  - Mr not zero otherwise
    !  - Here we need to check if zero(2) = .true. and avoid extra computation
    !---------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! K_{lambda, lambda}
    !     = Dl' * T2 * En * T2 * Dl
    !---------------------------------------------------------------------------
    ! if nu is known
    !  - K_{lambda, lambda} does not exist
    !  - nl = 0 => Dl = 0 (no need to check if zero(3) = .true.)
    !---------------------------------------------------------------------------
    if (nl == 0) return

    if(nr > 0 .and. (.not. zero(3))) then
       ! K_{lambda, rho}
       !   = Dr' * T1 * Emn * T2 * Dl + Mr' * T2 * En * T2 * Dl
       do i = (nr + 1), (nr + nl)
          do j = 1, nr
             K(i, j) = sum(Dr(:, j) * T1 * E(:, 2) * T2 * Dl(:, i - nr))
             K(j, i) = K(i, j)
          end do
       end do
       ! if nu does not have tha MA part, skip
       if(.not. zero(2)) then
          do i = (nr + 1), (nr + nl)
             do j = 1, nr
                K(i, j) = K(i, j) + sum(Mr(:, j) * E(:, 3) * T2 **2 * Dl(:, i - nr))
                K(j, i) = K(i, j)
             end do
          end do
       end if
    end if

    ! K_{lambda, lambda}
    !  = Dl' * T2 * En * T2 * Dl
    do i = (nr + 1), (nr + nl)
       do j = (nr + 1), i
          K(i, j) = sum(E(:, 3) * T2**2 * Dl(:, j - nr) * Dl(:, i - nr))
          K(j, i) = K(i, j)
       end do
    end do
    return
  end subroutine calc_K

  subroutine K_generic(SI, mu, nu, fita, fitb, fitar, fitma, fitd, m, n, npar, K, argsD)
    !---------------------------------------------------------------------------------------
    !  Calculate the information matrix K for a generic model. The blocks in K are
    !      Krr = Dr' * T1 * Em * T1 * Dr + Mr' * T2 * Emn * T1 * Dr +
    !            Dr' * T1 * Emn * T2 * Mr + Mr' * T2 * En * T2 * Mr
    !      Krl = Klr' = Dr' * T1 * Emn * T2 * Dl + Mr' * T2 * En * T2 * Dl
    !      Kll = Dl' * T2 * En * T2 * Dl
    ! Requires: mu, nu, Drho, Mrho, Dlambda, T1 and T2
    !---------------------------------------------------------------------------------------
    ! Input
    !   SI   : an argsSI object with previously calculated values of
    !              T1, T2, Drho = deta1/drho, Mrho = deta2/drho, and Dlambda = deta2/dlambda
    !          On exit, will also have the values of
    !              E_mu = -E(d2l/dmu2), E_(mu, nu) = -E(d2l/dmUdnu) and E_nu = -E(d2l/dnu2)
    !          saved in the argument E.
    !   mu   : the conditional time series mut
    !   nu   : the conditional time series nut
    !   fita : 1 = fit alpha, 0 = alpha is known
    !   fitb : number of unknown parameters in the beta vector
    !   fitar: number of unknown parameters in the phi (AR) vector
    !   fitma: number of unknown parameters in the theta (MA) vector
    !   fitd : 1 = fit d, 0 = d is known
    !   m    : starting point for the sum of the log-likelihood
    !   n    : sample size
    !   npar : number of unknown parameters
    !   argsD: an argsDist variable with extra information to be passed to the density function
    !
    ! Output
    !   K: the information matrix
    !---------------------------------------------------------------------------------------
    ! Last revision: July, 2023
    ! - removed fitnu argument
    ! - nu now has size n
    ! - fita, ... fitd and npar now have size 2
    ! - changed the code accordingly
    implicit none
    type(argsSI), intent(inout) :: SI
    integer,  intent(in)  :: fita(2), fitb(2), fitar(2), fitma(2), fitd(2)
    integer,  intent(in)  :: npar(2), n, m
    real(dp), intent(in)  :: mu(n), nu(n)
    real(dp), intent(out) :: K(max(1, sum(npar)), max(1, sum(npar)))
    type(argsDist), intent(inout) :: argsD
    real(dp) :: Drho(n, max(1, npar(1)))
    real(dp) :: Mrho(n, max(1, npar(1)))
    real(dp) :: Dlambda(n, max(1, npar(2)))
    integer  :: skip(2)
    logical  :: zero(3)

    ! checking if any part of the model is known
    skip(1) = 0
    skip(2) = 0
    if (npar(1) == 0) skip(1) = 1
    if (npar(2) == 0) skip(2) = 1

    Drho = 0.d0
    Mrho = 0.d0
    Dlambda = 0.d0
    zero(1:3) = .true.

    ! Calculate E = (E_mu, E_(mu, nu), E_nu)
    call safe_allocate(SI%E, n, 3)
    call argsD%Ed2llk_dist(m, n, mu, nu, skip, SI%E)


    ! Drho, Mrho and Dlambda
    if (skip(1) == 0) then
       ! if part 1 is fixed, skip Drho
       call fill_D(SI%deta(1,1), fita(1), fitb(1), fitar(1), fitma(1), fitd(1), n, npar(1), &
            Drho, zero(1))
    end if

    ! if part 1 or part 2 is fixed skip Mrho
    if (sum(skip) == 0) then
       call fill_D(SI%deta(2,1), fita(1), fitb(1), fitar(1), fitma(1), fitd(1), n, npar(1), &
            Mrho, zero(2))
    end if

    ! if part 2 is fixed skip Dlambda
    if (skip(2) == 0) then
       call fill_D(SI%deta(2,2), fita(2), fitb(2), fitar(2), fitma(2), fitd(2), n, npar(2),  &
            Dlambda, zero(3))
    end if

    ! calculate K
    call calc_K(n, npar(1), npar(2), SI%T(1)%z, SI%T(2)%z, SI%E, Drho, Mrho, Dlambda, zero, K)
    return
  end subroutine K_generic

  !-------------------------------------------------------------------------------------------------
  !
  !                 Optimization subroutines
  !
  !  Available methods
  !  - Nelder-Mead : requires loglik(model, npar, par, sll)
  !  - L-BFGS-B : requires loglik(model, npar, par, sll, score)
  !
  !-------------------------------------------------------------------------------------------------
  subroutine loglik_dist_nelder(loglik, model, npar, par, sll)
    !------------------------------------------------------------------
    !
    !   Subroutine to be used in Nelder-Mead optimization subroutine
    !
    !------------------------------------------------------------------
    implicit none
    class(optimFunc), intent(inout) :: loglik
    type(argsModel), intent(inout) :: model
    integer,  intent(in)  :: npar
    real(dp), intent(in)  :: par(npar)
    real(dp), intent(out) :: sll
    real(dp) :: par_aux(npar), U(npar)

    loglik%dummy = .true.
    par_aux = par
    !----------------------------------------------------
    ! Back to original scale
    !----------------------------------------------------
    call  transform_par(par_aux, npar, model%bounds, .true.)
    model%sco = 0

    call loglik_generic(model, npar, par_aux, sll, U)
    if (sll <  -Huge(1.d0)) sll =  -Huge(1.d0)
    if (sll > Huge(1.d0)) sll = Huge(1.d0)
    return
  end subroutine loglik_dist_nelder

  subroutine loglik_dist(loglik, model, npar, par, sll, U)
    !------------------------------------------------------------------
    !
    !   Subroutine to be used in L-BFGS-B optimization subroutine
    !
    !------------------------------------------------------------------
    implicit none
    class(optimFunc), intent(inout) :: loglik
    type(argsModel), intent(inout) :: model
    integer, intent(in) :: npar
    real(dp), intent(in) :: par(npar)
    real(dp), intent(out) :: sll, U(npar)

    loglik%dummy = .true.
    call loglik_generic(model, npar, par, sll, U)
    if (sll <  - Huge(1.d0)) sll =  - Huge(1.d0)
    if (sll > Huge(1.d0)) sll = Huge(1.d0)
    return
  end subroutine loglik_dist

  subroutine optim_nelder(optim, loglik, model, npar, par, nbd, bounds, sll, score, cf1, nc2, cf2, &
       neval, conv)
    !---------------------------------------------------------------------------------------
    !  Subroutine to perform parameter estimation using the Nelder-Mead algorithm
    !---------------------------------------------------------------------------------------
    ! Input
    !   loglik: the subroutine used to calculate the sum of the log-likelihood, with the following
    !           arguments loglik(model, npar, par, sll)
    !   npar  : number of unknown parameters
    !   nbd   : an integer array of dimension npar indicating the type of bounds imposed on the
    !           variables:
    !           nbd(i) = 0 if par(i) is unbounded,
    !                    1 if par(i) has only a lower bound,
    !                    2 if par(i) has both lower and upper bounds,
    !                    3 if par(i) has only an upper bound.
    !   bounds: lower and upper bounds for the parameters.
    !           lower(i)/upper(i) is ignored if par(i) only has an upper/lower bound
    !           or if it is unbounded.
    !   iprint: print control parameter
    !           < 0 No printing
    !           = 0 Printing of parameter values and the function value after initial evidence
    !               of convergence.
    !           > 0 As for iprint = 0 plus progress reports after every Iprint evaluations, plus
    !               printing for the initial simplex.
    !   stopcr: stopping criterion. The criterion is applied to the standard deviation of the values
    !           of the objective function at the points of the simplex.
    !   maxit : the maximum no. of function evaluations allowed.
    !
    ! Input/Output
    !   model: input, the model's configurations
    !          output, updated values of the conditional time series.
    !   par  : input, starting values of parameters
    !          output, final values of parameters
    !
    ! Output
    !   sll  : The function value (sum of the log-likelihood) corresponding to the final
    !          parameter's values.
    !   neval: number of function evaluations performed by the algorithm
    !   conv : convergence control
    !          = 0 for successful termination
    !          = 1 if maximum no. of function evaluations exceeded
    !          = 2 if nop < 1 .or. stopcr <= 0
    !---------------------------------------------------------------------------------------
    ! Last revision: July, 2023
    implicit none
    class(optimizer), intent(inout) :: optim
    type(optimFunc), target, intent(inout) :: loglik
    type(argsModel), intent(inout) :: model
    integer,  target, intent(in)   :: cf1(3)
    integer,  pointer       :: iprint, maxit
    integer,  intent(in)    :: npar, nc2
    integer,  intent(in)    :: nbd(npar)
    integer,  intent(out)   :: neval, conv
    real(dp), intent(in)    :: bounds(npar, 2)
    real(dp), target, intent(in) :: cf2(nc2)
    real(dp), pointer       :: stopcr
    real(dp), intent(out)   :: sll, score(max(1, npar*cf1(3)))
    real(dp), intent(inout) :: par(npar)
    real(dp) :: step(npar)

    optim%dummy = .true.
    conv = 0
    iprint => cf1(1)
    maxit  => cf1(2)
    stopcr => cf2(1)
    score = 0.d0

    ! allocating bounds and transforming the parameters to (-infty, infty)
    call set_bounds(model%bounds, bounds, nbd, max(1, npar))
    call transform_par(par, npar, model%bounds, .false.)

    ! step for Nelder-Mead. Based on the transformed parameters
    step = max(0.10d0 * abs(par), 0.00025d0)
    where(bounds(:,1) == bounds(:,2) .and. nbd == 2) step = 0.d0

    ! Nelder-Mead
    call minim(par, step, npar, sll, maxit, iprint, stopcr, loglik, conv, neval, model)

    ! transforming the parameters for the original scale
    call transform_par(par, npar, model%bounds, .true.)

    return
  end subroutine optim_nelder

  subroutine optim_lbfgsb(optim, loglik, model, npar, par, nbd, bounds, sll, score, cf1, nc2, cf2, &
       neval, conv)
    !---------------------------------------------------------------------------------------
    !  Subroutine to perform parameter estimation using the L-BFGS-B algorithm
    !---------------------------------------------------------------------------------------
    ! Input
    !   loglik: the subroutine used to calculate the sum of the log-likelihood, with the following
    !           arguments loglik(model, npar, par, sll, score)
    !   npar  : number of unknown parameters
    !   nbd   : an integer array of dimension npar indicating the type of bounds imposed on the
    !           variables:
    !           nbd(i) = 0 if par(i) is unbounded,
    !                    1 if par(i) has only a lower bound,
    !                    2 if par(i) has both lower and upper bounds,
    !                    3 if par(i) has only an upper bound.
    !   bounds: lower and upper bounds for the parameters.
    !           lower(i)/upper(i) is ignored if par(i) only has an upper/lower bound
    !           or if it is unbounded.
    !   iprint: print control parameter
    !           < 0 No printing
    !           = 0 Printing of parameter values and the function value after initial evidence
    !               of convergence.
    !           > 0 As for iprint = 0 plus progress reports after every Iprint evaluations, plus
    !               printing for the initial simplex.
    !   factr : stopping criterion. The iteration will stop when
    !               (f^k - f^{k + 1}) / max{|f^k|, |f^{k + 1}|, 1} <= factr * epsmch
    !           where epsmch is the machine precision, which is automatically generated by the code.
    !           Typical values for factr: 1.d+12 for low accuracy;
    !                                     1.d+7 for moderate accuracy;
    !                                     1.d+1 for extremely high accuracy.
    !   pgtol : stopping criterion. The iteration will stop when
    !                   max{|proj g_i | i = 1, ..., n} <= pgtol
    !           where pg_i is the ith component of the projected gradient.
    !   maxit : the maximum no. of function evaluations allowed.
    !
    ! Input/Output
    !   model: input, the model's configurations
    !          output, updated values of the conditional time series, score vector and related
    !          matrices and vectors
    !   par  : input, starting values of parameters
    !          output, final values of parameters
    !   nbd  : input, bounds for the parameters
    !
    ! Output
    !   sll  : The function value (sum of the log-likelihood) corresponding to the final
    !          parameter's values.
    !   score: the score vector corresponding to the final parameter's values
    !   neval: number of function evaluations performed by the algorithm
    !   conv : convergence control
    !          = 0 for successful termination
    !          = 1 fail. See L-BFGS-B for details
    !   convm: convergence message
    !---------------------------------------------------------------------------------------
    ! Last revision: February, 2025
    !  - removed the goto loops
    implicit none
    class(optimizer), intent(inout) :: optim
    type(optimFunc), target, intent(inout) :: loglik
    type(argsModel), intent(inout) :: model
    integer,  target, intent(in)   :: cf1(3)
    integer,  pointer       :: iprint, maxit
    integer,  intent(in)    :: npar, nc2, nbd(npar)
    integer,  intent(out)   :: neval, conv
    real(dp), target, intent(in) :: cf2(nc2)
    real(dp), pointer       :: factr, pgtol
    real(dp), intent(in) :: bounds(npar, 2)
    real(dp), intent(out)   :: sll, score(max(1, npar*cf1(3)))
    real(dp), intent(inout) :: par(npar)
    character(len = 60), pointer :: convm
    integer, parameter ::  mmax = 5 ! default in r
    real(dp) ::  dsave(29)
    ! wa must have size 2 * mn + 11m**2 + 5n + 8m
    real (dp) :: wa(2 * mmax * npar + 12 * mmax * mmax + 5 * npar + 8 * mmax)
    integer   :: iwa(3 * npar), isave(44)
    integer   :: niter
    logical   :: lsave(4)
    character (len=60) :: csave

    optim%dummy = .true.
    convm  => loglik%message
    iprint => cf1(1)
    maxit  => cf1(2)
    factr  => cf2(1)
    pgtol  => cf2(2)

    wa = 0.d0
    iwa = 0
    niter = 0

    ! l-bfgs-b
    !---------------------------------------------------------
    ! setting lower bounds, upper bounds and fixed values
    !---------------------------------------------------------
    convm = 'start'
    model%sco = 0
    call loglik%loglik(model, npar, par, sll, score)
    model%sco = 1
    conv = 0

    if (maxit == 0) return
    conv = 1

    do while (.true.)

       ! Check if we've exceeded max iterations
       if (niter > maxit) then
          convm = "max number of iteration reached"
          exit
       end if

       niter = niter + 1

       call setulb(npar, mmax, par, bounds(:,1), bounds(:,2), nbd, sll, score, factr, pgtol, &
            wa, iwa, convm, iprint, csave, lsave, isave, dsave)

       if (convm(1:2) .eq. 'fg') then
          call loglik%loglik(model, npar, par, sll, score)
          neval = isave(13)
          cycle  ! Continue to the next iteration
       else if (convm(1:5) .eq. 'new_x') then
          cycle  ! Continue to the next iteration
       else
          ! If task is neither 'fg' nor 'new_x', we terminate execution.
          exit
       end if
    end do

    if (convm(1:4) == 'conv') conv = 0

    return
  end subroutine optim_lbfgsb

  !-------------------------------------------------------------------------------------------------
  !
  !  Subroutines used to initialize/update the vector of parameters
  !
  !-------------------------------------------------------------------------------------------------
  subroutine start_par1(par, model, part)
    !---------------------------------------------------------------------------------------
    ! Initializes the non-fixed parameters for each part of the model
    ! This subroutine is called to set the parameter values during the log-likelihood
    ! evaluation (during and after estimation process)
    !---------------------------------------------------------------------------------------
    ! Input
    !   par : the value of the non-fixed parameters
    !   part: which part of the model needs initialization
    !
    ! Input/Output
    !   model: input, the model's configurations partially defined.
    !          output, updated model's configurations
    !---------------------------------------------------------------------------------------
    ! Last revision: July, 2023
    ! - removed nu related settings
    implicit none
    type(argsModel), intent(inout) :: model
    real(dp), intent(in) :: par(sum(model%pt%npar))
    integer,  intent(in) :: part
    integer :: n1, n2

    n1 = (2 - part) + (part - 1) * (model%pt(1)%npar + 1)
    n2 = n1 - 1

    ! alpha
    n1 = n2 + 1
    n2 = n2 + model%pt(part)%alpha%fit
    if (n2 >= n1) model%pt(part)%alpha%cf = par(n1:n2)

    ! beta
    n1 = n2 + 1
    n2 = n2 + model%pt(part)%beta%fit
    if (n2 >= n1) model%pt(part)%beta%cf(model%pt(part)%beta%lags) = par(n1:n2)

    ! phi - ar
    n1 = n2 + 1
    n2 = n2 + model%pt(part)%ar%fit
    if (n2 >= n1) model%pt(part)%ar%cf(model%pt(part)%ar%lags) = par(n1:n2)

    ! theta - ma
    n1 = n2 + 1
    n2 = n2 + model%pt(part)%ma%fit
    if (n2 >= n1) model%pt(part)%ma%cf(model%pt(part)%ma%lags) = par(n1:n2)

    ! d
    n1 = n2 + 1
    n2 = n2 + model%pt(part)%d%fit
    if (n2 >= n1) model%pt(part)%d%cf = par(n1:n2)

    return
  end subroutine start_par1

  subroutine start_par2(par, model, vc, part)
    !---------------------------------------------------------------------------------------
    ! Initializes the non-fixed parameters for each part of the model. Also initializes the
    ! coefficients of (1 - L)^{-d}
    ! This subroutine is called to set the parameter values during the log-likelihood
    ! evaluation (during and after estimation process)
    !---------------------------------------------------------------------------------------
    ! Input
    !   par : the value of the non-fixed parameters
    !   part: which part of the model needs initialization
    !
    ! Input/Output
    !   model: input, the model's configurations partially defined.
    !          output, updated model's configurations
    !
    ! Output
    !  vc: the coefficients of (1 - L)^{ - d}
    !---------------------------------------------------------------------------------------
    ! Last revision: July, 2023
    implicit none
    type(argsModel), intent(inout) :: model
    real(dp), intent(in)  :: par(sum(model%pt%npar))
    integer,  intent(in)  :: part
    real(dp), intent(out) :: vc(0:model%pt(part)%inf)

    ! initializing the parameters
    call start_par1(par, model, part)

    ! initializing the coefficients c_k
    call vc_f(model, vc, part)
    return
  end subroutine start_par2


  subroutine start_par12(par, model)
    !---------------------------------------------------------------------------------------
    ! Initializes the non-fixed parameters the entire model
    !---------------------------------------------------------------------------------------
    ! Input
    !   par : the value of the non-fixed parameters
    !
    ! Input/Output
    !   model: input, the model's configurations partially defined.
    !          output, updated model's configurations
    !---------------------------------------------------------------------------------------
    ! March, 2025: added this surboutine
    implicit none
    type(argsModel), intent(inout) :: model
    real(dp), intent(in) :: par(sum(model%pt%npar))
    integer :: part

    do part = 1, 2
       call start_par1(par, model, part)
    end do
    return
  end subroutine start_par12

  !-------------------------------------------------------------------------------------------------
  !
  !                  Subroutines used to process information in and out the model
  !
  !-------------------------------------------------------------------------------------------------
  subroutine set_link_to_model(model, link)
    implicit none
    type(argsModel), intent(inout) :: model
    type(argsLink), intent(in) :: link(8)
    model%pt(1)%linkg = link(1:4)
    model%pt(2)%linkg = link(5:8)
    return
  end subroutine set_link_to_model

  subroutine set_link_to_link(link, lconfig, argsL)
    !---------------------------------------------------------------------------
    ! Initializes the arguments in the link and updates the limits for gnu
    !---------------------------------------------------------------------------
    ! Input
    !  link: a size 8 vector with the link codes.
    !        g1 (dummy), g11, g12, g13, g2, g21, g22, g23
    !
    ! Input/Output
    !  lconfig: a matrix of size 8 by 4 with configurations for the links
    !
    ! Output
    !   argsL: a size 8 argslink type variable with information on the links
    !---------------------------------------------------------------------------
    ! February, 2025: added this subroutine
    implicit none
    integer,  intent(in)    :: link(8)
    real(dp), intent(inout) :: lconfig(8, 4)
    type(argslink), intent(inout) :: argsL(8)
    real(dp) :: a, b

    argsL(:)%link   = link
    argsL(:)%lower  = lconfig(:, 1)
    argsL(:)%upper  = lconfig(:, 2)
    argsL(:)%par(1) = lconfig(:, 3)
    argsL(:)%par(2) = lconfig(:, 4)

    ! check if g11 and g12 are the same
    call check_update(argsL(2), argsL(3))

    !----------------------------------------
    ! Update for part 2
    !----------------------------------------
    !  - update bounds for  gi1(g(wt))
    !  - update bounds for  gi2(g(wt))
    a = linkfun(argsL(5)%lower, argsL(5))  ! compute g2(wt_lower)
    b = linkfun(argsL(5)%upper, argsL(5))  ! compute g2(wt_upper)
    lconfig(6:7, 1) = min(a, b)
    lconfig(6:7, 2) = max(a, b)
    argsL(6:7)%lower = lconfig(6:7, 1)
    argsL(6:7)%upper = lconfig(6:7, 2)

    ! check if g21 and g22 are the same
    call check_update(argsL(6), argsL(7))

    !----------------------------------------------
    ! error term:
    !  - in part 1: only the code is used
    !  - in part 2: depends on the code of part 1
    !----------------------------------------------
    if (link(4) == 0) then ! a < y, mu < b
       argsL(8)%lower =  lconfig(3, 1) - lconfig(3, 2) ! a - b
       argsL(8)%upper =  lconfig(3, 2) - lconfig(3, 1) ! b - a
    else ! -Inf < g(y), g(mu) < Inf
       argsL(8)%lower =  -Infinity
       argsL(8)%upper =  Infinity
    end if
    lconfig([4,8],1) = argsL(8)%lower
    lconfig([4,8],2) = argsL(8)%upper
    argsL(8)%par = lconfig(8,3:4)

    return
  end subroutine set_link_to_link

  subroutine gy_update(n, y, gy, escale, p, linkg, ierr)
    !------------------------------------------------------
    ! Compute g11 and g12
    !------------------------------------------------------
    ! Input
    !   n      : sample size
    !   y      : the observed time series
    !   escale : controls the scale of the error.
    !           (0 = data scale. 1 = predictive scale)
    !   p      : the order phi (AR) polynomial
    !   linkg  : a size 2 argslink type variable with
    !            information on the links
    !
    ! output
    !   gy  : matrix with two columns corresponding to the
    !         transformed time series g11(y) and g12(y)
    !   ierr: error code.
    !           0    = no error
    !           11   = fail to compute g11
    !           12   = fail to compute g12
    !           1112 = fail to compute g11 and g12
    !-----------------------------------------------------
    ! February, 2025: added this subroutine
    implicit none
    integer,  intent(in)  :: n, escale, p
    integer,  intent(out) :: ierr
    real(dp), intent(in)  :: y(n)
    real(dp), intent(out) :: gy(n,2)
    type(argsLink), intent(inout) :: linkg(2)
    integer :: t, code

    gy = 0.d0
    code = 0
    ierr = 0

    ! g11(y) - only needed if escale == 1
    if(escale == 1) then
       do t = 1, n
          gy(t,1) = linkfun(y(t), linkg(1))
       end do
       code = code + 1
    end if

    ! g_12(y) is necessary only if p > 0
    if (p > 0) then
       ! No need to recompute only if link(1) = link(2) and escale = 1
       if ((.not. linkg(2)%update) .and. (escale == 1)) then
          gy(:,2) = gy(:,1)
       else
          do t = 1, n
             gy(t,2) = linkfun(y(t), linkg(2))
          end do
       end if
       code = code + 1
    end if

    if(code == 0) return

    ! check for errors in computing g11 and g12
    if(.not. all(is_finite(gy))) then
       code = 0
       call labelpr('----------------------------------------------------', -1)
       call labelpr(' Please select another link', -1)
       if (any(isnan(gy(:,1)))) then
          call labelpr(' Fail to evaluate g11(y)', -1)
          code = 1
       else if (linkg(2)%update .and. any(isnan(gy(:,2)))) then
          call labelpr(' Fail to evaluate g12(y)', -1)
          code = code + 2
       end if
       call labelpr('----------------------------------------------------', -1)
       select case (code)
       case (1)
          ierr = 11
       case (2)
          ierr = 12
       case (3)
          ierr = 1112
       end select
    end if

    return
  end subroutine gy_update

  subroutine g_start_update(x, gx, link, ierr, part)
    !------------------------------------------------------
    ! Compute g(x) for a given x.
    ! Used to update g(ystart) and g(nu0)
    !------------------------------------------------------
    ! Input
    !   x     : the starting value
    !   linkg : an argslink type variable with
    !            information on the links
    !
    ! output
    !   gx  : the transformed value g(x)
    !   ierr: error code. 0 = no error, 999 = fail
    !-----------------------------------------------------
    ! February, 2025: added this subroutine
    implicit none
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: gx
    integer,  intent(in)  :: part
    integer,  intent(out) :: ierr
    type(argslink), intent(in) :: link

    gx = 0.d0
    ierr = 0

    if(x < link%lower .or. x > link%upper) then
       call labelpr('----------------------------------------------------', -1)
       call labelpr(' Please select another starting value', -1)
       if(part == 1) then
          ierr = 91
          call labelpr(' Fail to compute g12(y.start)', -1)
          call labelpr(' y.start is out of bounds', -1)
       end if
       if(part == 2) then
          ierr = 92
          call labelpr(' Fail to compute g22(vt.start)', -1)
          call labelpr(' vt.start is out of bounds', -1)
       end if
       call labelpr('----------------------------------------------------', -1)
       return
    end if

    gx = linkfun(x, link)
    if (.not. is_finite(gx)) then
       call labelpr('----------------------------------------------------', -1)
       call labelpr(' Please select another starting value', -1)
       if(part == 1) then
          call labelpr(' Fail to compute g12(y.start)', -1)
          ierr = 9991
       end if
       if(part == 2) then
          call labelpr(' Fail to compute g22(vt.start)', -1)
          ierr = 9992
       end if
       call labelpr('----------------------------------------------------', -1)
    end if
    return

  end subroutine g_start_update

  subroutine g_update(n, y, gy, escale, p, linkg1, ystart, &
       g2ystart, vtstart, g2vstart, linkg2, ierr)
    !------------------------------------------------------
    ! Compute g11(y), g12(y)
    ! Updates g12(ystart) and g12(g(nu0))
    !------------------------------------------------------
    ! Input
    !   n      : sample size
    !   y      : the observed time series
    !   escale : controls the scale of the error.
    !            (0 = data scale. 1 = predictive scale)
    !   p      : the order phi (AR) polynomial
    !   ystart : starting value for yt
    !   vtstart: starting value for g(nut)
    !   linkg2 : an argslink type variable with information
    !            on the links for part 2
    !
    ! Input/Output
    !   linkg1 : a size 2 argslink type variable with
    !            information on the links for part 1
    !
    ! Output
    !   gy  : matrix with two columns corresponding to the
    !         transformed time series g11(y) and g12(y)
    !   ierr: error code.
    !           0    = no error
    !           11   = fail to compute g11
    !           12   = fail to compute g12
    !           1112 = fail to compute g11 and g12
    !           999  = fail to update starting values
    !   g2ystart : starting value for g12(yt)
    !   g2vtstart: starting value for g22(g(nut))
    !----------------------------------- ------------------
    ! February, 2025: added this subroutine
    implicit none
    integer,  intent(in)  :: n, escale, p(2)
    integer,  intent(out) :: ierr
    real(dp), intent(in)  :: y(n), ystart, vtstart
    real(dp), intent(out) :: gy(n,2), g2ystart, g2vstart
    type(argsLink), intent(inout) :: linkg1(2)
    type(argsLink), intent(in) :: linkg2

    ! g11(yt) and g12(yt)
    call gy_update(n, y, gy, escale, p(1), linkg1, ierr)
    if(ierr > 0) return

    if(p(1) > 0) then
       ! g12(ystart)
       call g_start_update(ystart, g2ystart, linkg1(2), ierr, 1)
       if(ierr > 0) return
    end if

    if(p(2) > 0) then
       ! g22(vtstart); vtstart = g(nustart)
       call g_start_update(vtstart, g2vstart, linkg2, ierr, 2)
       return
    end if
  end subroutine g_update

  subroutine get_model(model, n, order, y, gy, xreg1, xreg2, tstart, xstart,  &
       link, lconfig, skippar, npar, par, xregar, nfix, alpha, flagsb, fvbeta, &
       flagsar, fvar, flagsma, fvma, d, extras, full, ierr)
    !---------------------------------------------------------------------------------------
    !  Subrotuine used to pass the values entered in the main program to the user defined variables
    !  that will be passed to the generic subroutines
    !---------------------------------------------------------------------------------------
    ! Input
    !  n      : sample size
    !  order  : a 2 by 4 matrix with (nreg, p, q, inf)
    !  y      : the observed time series yt
    !  xreg1  : matrix of regressors corresponding to part 1
    !  xreg2  : matrix of regressors corresponding to part 2
    !  tstart : starting value for y, g2(nu) and g23(e1). To be used when t < 1
    !  xstart : matrix of maxval(nreg) rows by 2 rows with starting values for the regressors
    !           (to be used when t < 1)
    !  link   : a vector of size 8, indicating the links to be used in the model
    !  lconfig: a matrix of size 8 by 4 with
    !            - the lower and upper bounds related to the link functions
    !            - the constant and power arguments for polynomial/SIP links
    !  skippar: 0 = parameters must be allocated. 1 = skip allocation (for BARC models)
    !  npar   : number of non-fixed parameters
    !  par    : the starting / final value of the non-fixed parameters
    !  xregar : a vector of size 2.
    !           0 = xreg is included only in the intercept
    !           1 = xreg is also included in the AR part.
    !  nfix   : a 2 by 5 matrix with the number of fixed parameters
    !  alpha  : a vector of size 2.
    !            if fixa(i) = 1, the value of the parameter alpha corresponding to part i.
    !            if fixa(i) = 0, a dummy argument.
    !  flagsb : matrix with 2 rows with the lags that must be fixed in each beta
    !  fvbeta : matrix with 2 rows with the value of the fixed parameters in each beta
    !  flagsar: matrix with 2 rows with the lags that must be fixed in each phi
    !  fvar   : matrix with 2 rows with the value of the fixed parameters in each phi
    !  flagsma: matrix with 2 rows with the lags that must be fixed in each theta
    !  fvma   : matrix with 2 rows with the value of the fixed parameters in each theta
    !  d      : a vector of size 2.
    !            if fixd(i) = 1, the value of the parameter d corresponding to part i.
    !            if fixd(i) = 0, a dummy argument.
    !  extras  : a vector containing
    !          - m: an integer indicating where the sum of the log-likelihood starts
    !          - llk: indicates if the sum of the log-likelihood must be returned
    !          - sco: indicates if the score vector must be returned
    !          - info: indicates if the information matrix must be returned
    !  full   : logical, indicates if all calculations must be done.
    !           If set to false, some quatities are skiped and must be previously computed
    !
    ! Input / Output
    !   model: input, the model's configurations partially defined.
    !          output, updated model's configurations
    !
    ! Output
    !   gy :  matrix with two columns corresponding to the transformed time series
    !         g11(y) and g12(y)
    !  ierr: error code
    !
    !---------------------------------------------------------------------------------------
    ! July, 2023
    !  - removed fixnu, nu and related code
    !  - nreg, npar, fixa, ... now have size 2
    !  - the vectors of fixed lags and values are now matrices with two columns
    !  - model%argsL no longer needs allocation
    !  - removed ylower and yupper and added lconfig
    !  - made the necessary changes to work with part 1 and 2
    !
    ! February, 2024
    !  - Updated the description
    !
    ! Last revision: March, 2025
    !  - gy now is a matrix with 2 columns.
    !  - renamed lconfig to linkinfo and added two more columns.
    !    Now the columns correspond to lower, upper, ctt and power
    !  - removed escale and added the code to the link variable
    !  - lconfig now has 8 rows: g1 (dummy), g11, g12, g13, g2, g21, g22, g23
    !  -  replaced nreg, p, q, inf by order
    !  - replaced ystart and gnustart by tstart
    !  - replaced fixa, fixb, fixar, fixma and fixd by the matrix nfix
    implicit none
    type(argsModel), intent(inout) :: model
    integer,  intent(in)    :: order(2,4), nfix(2,5)
    integer,  intent(in)    :: n, link(8), skippar
    integer,  intent(in)    :: npar(2), xregar(2), extras(4)
    integer,  intent(in)    :: flagsb(2, maxval([1, nfix(:, 2)]))
    integer,  intent(in)    :: flagsar(2, maxval([1, nfix(:, 3)]))
    integer,  intent(in)    :: flagsma(2, maxval([1, nfix(:, 4)]))
    integer,  intent(inout) :: ierr
    real(dp), intent(in)    :: y(n), tstart(2)
    real(dp), intent(in)    :: xstart(2, maxval([1, order(:, 1)]))
    real(dp), intent(in)    :: xreg1(n, max(1, order(1, 1)))
    real(dp), intent(in)    :: xreg2(n, max(1, order(2, 1)))
    real(dp), intent(in)    :: par(sum(npar)), alpha(2), d(2)
    real(dp), intent(in)    :: fvbeta(2, maxval([1, nfix(:, 2)]))
    real(dp), intent(in)    :: fvar(2, maxval([1, nfix(:, 3)]))
    real(dp), intent(in)    :: fvma(2, maxval([1, nfix(:, 4)]))
    real(dp), intent(inout) :: lconfig(8, 4)
    real(dp), intent(out)   :: gy(n,2)
    logical,  intent(in)    :: full
    type(argsLink) :: linkg(8)

    !----------------------------------------------------------------
    ! setting the link parameters and calculating
    ! - g11(y) and g12(y)
    ! - g21(g(nu0))
    ! - g23(e10)
    !----------------------------------------------------------------
    call set_link(link, lconfig, linkg)
    call set_link(model, linkg)

    ! skip this step when invoking the subroutine to initialize the model for predictions
    if(full) then
       ! g11(yt), g12(yt), g12(ystart) and g22(vtstart); vtstart = g(nustart)
       call g_update(n, y, gy, link(4), order(:,2), model%pt(1)%linkg(2:3), tstart(1), &
            model%cts(1)%g2start, tstart(2), model%cts(2)%g2start, model%pt(2)%linkg(3), ierr)
       if(ierr > 0) return
    end if

    !-------------------------------------------------------------
    ! allocating the time series and parameters
    ! setting fixed values / lags of parameters
    !-------------------------------------------------------------
    call allocate_model(model, n, order, y, gy, xreg1, xreg2, xstart, xregar, nfix, alpha, &
         flagsb, fvbeta, flagsar, fvar, flagsma, fvma, d)

    !--------------------------------------------------------------
    !  Setting the parameter's non-fixed values:
    !  alpha, beta, ar, ma and d
    !  BARC models must skip this step an call start_par later
    !--------------------------------------------------------------
    if (skippar == 0) then
       call start_par(par, model)  ! part 1 and 2
    end if

    ! the prediction subroutine does not require likelihood related variables
    if(full) then
       !-------------------------------------------------------------
       ! log-likelihood, score vector and information matrix
       !-------------------------------------------------------------
       model%m = extras(1)
       model%llk = extras(2)
       model%sco = extras(3)
       model%info = extras(4)

       if (model%sco + model%info == 0) return
       !---------------------------------------------------------------------
       ! allocating score - vector, Information matrix and related matrices
       !---------------------------------------------------------------------
       call allocate_SI(model, model%SI)
    end if

    return
  end subroutine get_model

  subroutine return_model(model, n, ts, inf, extra, nd, D, T, E, h)
    !---------------------------------------------------------------------------------------
    !  Subrotuine used to pass the values in the model variable to matrices and vectors
    !  to return to R.
    !---------------------------------------------------------------------------------------
    ! Input
    !   model: the model's configurations and estimation results
    !   n    : the sample size
    !   extra: 1 = extra matrices must be returned. 0 = matrices are not required
    !   nd   : number of non-fixed parameters
    !
    ! Output
    !   ts   : the extracted time series (mut, eta1t, e1t, nut, gnut, eta2t)
    !           - the conditional time series mut
    !           - the conditional time series g_11(mut)
    !           - the error term e1t
    !           - the conditional time series nut
    !           - the conditional time series g(nut)
    !           - the conditional time series g_21(gnu)
    !           - the error term e2t
    !   inf  : the efective value used as truncation point in the infinite sums
    !   D    : matrix with the derivatives deta/dlambda
    !   T    : matrix with the derivatives dmu/deta and dnu/deta2
    !   E    : matrix with the expected value of the second derivatives
    !          d2l/dmu2, d2l/dmUdnu, d2l/dnu2
    !   h    : matrix with the derivatives dl/dmu and dl/dnu
    !---------------------------------------------------------------------------------------
    ! July, 2023
    !  - added nu, gnu, eta2
    !  - inf now has size 2
    !  - T and h now can have two columns
    !  - E now can have 3 columns
    !---------------------------------------------------------------------------------------
    ! Last revision: February, 2024
    ! - Updated the description
    !
    implicit none
    type(argsModel), intent(in) :: model
    integer,  intent(in)  :: n, nd, extra
    integer,  intent(out) :: inf(2)
    real(dp), intent(out) :: ts(n, 7)
    real(dp), intent(out) :: T(max(n * extra, 1), max(2 * extra, 1))
    real(dp), intent(out) :: h(max(n * extra, 1), max(2 * extra, 1))
    real(dp), intent(out) :: D(max(n * extra, 1), max(nd * extra, 1))
    real(dp), intent(out) :: E(max(n * extra, 1), max(3 * extra, 1))
    integer :: i1, i2
    logical :: zero

    ! conditional time series
    ts(:,1) = model%cts(1)%w   ! mu
    ts(:,2) = model%cts(1)%eta ! eta1
    ts(:,3) = model%cts(1)%et  ! e1
    ts(:,4) = model%cts(2)%w   ! nu
    ts(:,5) = model%cts(2)%gw  ! g(nu)
    ts(:,6) = model%cts(2)%eta ! eta2
    ts(:,7) = model%cts(2)%et  ! e2

    inf = model%pt%inf

    D = 0.d0
    E = 0.d0
    h = 0.d0
    T = 0.d0
    if (extra == 0 .or. (model%sco + model%info) == 0) return

    i1 = 1
    i2 = 0

    ! if part 1 is fixed, skip
    if (model%pt(1)%npar > 0) then
       i2 = model%pt(1)%npar
       call fill_D(model%SI%deta(1,1), model%pt(1)%alpha%fit, model%pt(1)%beta%fit, &
            model%pt(1)%ar%fit, model%pt(1)%ma%fit, model%pt(1)%d%fit, n, i2, D(:, i1:i2), zero)
       T(:, 1) = model%SI%T(1)%z
       if (model%sco == 1) h(:, 1) = model%SI%h(1)%z

       ! if part1 or part 2 is fixed, skip
       if (model%pt(2)%npar > 0) then
          i1 = i2 + 1
          i2 = i2 + model%pt(1)%npar
          call fill_D(model%SI%deta(2,1), model%pt(1)%alpha%fit, model%pt(1)%beta%fit,  &
               model%pt(1)%ar%fit, model%pt(1)%ma%fit, model%pt(1)%d%fit, n, i2 - i1 + 1, &
               D(:, i1:i2), zero)
       end if
    end if

    ! if part 2 is fixed, skip
    if (model%pt(2)%npar > 0) then
       i1 = i2 + 1
       i2 = i2 + model%pt(2)%npar
       call fill_D(model%SI%deta(2,2), model%pt(2)%alpha%fit, model%pt(2)%beta%fit, &
            model%pt(2)%ar%fit, model%pt(2)%ma%fit, model%pt(2)%d%fit, n, i2 - i1 + 1, &
            D(:, i1:i2), zero)
       T(:, 1 + extra) = model%SI%T(2)%z
       if (model%sco == 1) h(:, 1 + extra) = model%SI%h(2)%z
    end if

    if (model%info == 1) E = model%SI%E(:, 1:(1 + 2 * extra))

    return
  end subroutine return_model

  subroutine final_model(model, npar, par, length, ts, nreg, xnew1, xnew2, forecast, &
       inf, sll, extras, U, K, nd, D, T, E, h)
    !---------------------------------------------------------------------------------------
    !   Reports the final model and related matrices
    !---------------------------------------------------------------------------------------
    ! Input
    !  model : the model's configurations and estimation results
    !  npar  : number of non-fixed parameters
    !  par   : estimated values for non-fixed parameters
    !  length: size 2 vector with n = sample size and nnew =  number of forecasts
    !  nreg  : a vector of size 2 with the number of regressors in each part of the model
    !  xnew1 : matrix with the out-of-sample values for the regressors in part 1
    !  xnew2 : matrix with the out-of-sample values for the regressors in part 2
    !  extras: a size 3 vector with the values of
    !          - sco: 0(1) = score vector is not required (must be returned)
    !          - info: 0(1) = information matrix is not required (must be returned)
    !          - extra: 0(1) = extra matrices are not required (must be returned)
    !
    ! Output
    !  ts    : the extracted time series mut, eta1t, e1t, nut, gnut, eta2t, e2t
    !  forecast: matrix with 6 columns with the predicted values of
    !            g12(y), mu, eta1, nu, gnu, eta2.
    !  inf     : the efective value used as truncation point in the infinite sums
    !  sll     : sum of the log-likelihood values
    !  U       : score vector
    !  K       : the information matrix
    !  D       : matrix with the derivatives deta/dlambda
    !  T       : matrix with the derivatives dmu/deta and dnu/deta2
    !  E       : matrix with the expected value of the second derivatives
    !            d2l/dmu2, d2l/dmUdnu, d2l/dnu2
    !  h       : matrix with the derivatives dl/dmu and dl/dnu
    !---------------------------------------------------------------------------------------
    ! August, 2023
    !  - removed fixnu
    !  - added nu, gnu, eta2
    !  - inf and nreg now have size 2
    !  - renamed xnew and added xnew2
    !  - renamed Drho to D
    !  - T and h now can have 2 columns
    !  - E now can have 3 columns
    !  - renamed ynew to forecast
    !  - added code to forecast nu and gnu
    !
    ! February, 2024
    ! - Updated the description
    !
    ! Last revision: March, 2025
    ! - replaced n and nreg with the vector length
    ! - replaced mu, eta1, error, nu, gnu and eta2 with the matrix ts
    ! - added pointers
    ! - replaced sco, info and extra  with the vector extra
    ! - forecast now has 6 columns
    implicit none
    type(argsModel), intent(inout) :: model
    integer,  intent(in)  :: npar, length(2), nreg(2), nd
    integer,  target, intent(in) :: extras(3)
    integer,  pointer     :: sco, info, extra
    integer,  intent(out) :: inf(2)
    real(dp), intent(in)  :: par(npar)
    real(dp), intent(in)  :: xnew1(max(1, length(2)), max(1, nreg(1)))
    real(dp), intent(in)  :: xnew2(max(1, length(2)), max(1, nreg(2)))
    real(dp), target, intent(out) :: ts(length(1), 7)
    real(dp), intent(out) :: sll, U(max(1, npar * extras(1)))
    real(dp), intent(out) :: forecast(max(1, length(2)), 5)
    real(dp), intent(out) :: K(max(1, npar * extras(2)), max(1, npar * extras(2)))
    real(dp), intent(out) :: T(max(length(1) * extras(3), 1), max(2 * extras(3), 1))
    real(dp), intent(out) :: h(max(length(1) * extras(3), 1), max(2 * extras(3), 1))
    real(dp), intent(out) :: D(max(length(1) * extras(3), 1), max(nd * extras(3), 1))
    real(dp), intent(out) :: E(max(length(1) * extras(3), 1), max(3 * extras(3), 1))
    real(dp) :: Uaux(npar)

    sco   => extras(1)
    info  => extras(2)
    extra => extras(3)

    !------------------------------------------------------------------------------
    ! reporting score vector and information matrix (if required)
    !------------------------------------------------------------------------------
    model%llk = 1
    if (sco + info > 0) then
       ! setting model%sco  = 1 so D and T will be calculated
       model%sco = max(sco, info)
       model%info = info
       call allocate_SI(model, model%SI)
    end if

    sll = 0.d0
    U = 0.d0
    K = 0.d0
    Uaux = 0.d0

    !---------------------------------------------------------------------------------
    ! Reports the final model.
    ! Calculates: mu, eta1, nu, gnu, eta2, rt,  - sll (if llk = 1) and  - U (if sco = 1).
    !---------------------------------------------------------------------------------
    call loglik_generic(model, npar, par, sll, Uaux)

    if (sco == 1) U = Uaux

    if (info == 1) then
       !------------------------------------------------------------------------------
       ! Information matrix
       !------------------------------------------------------------------------------
       call K_generic(model%SI, model%cts(1)%w, model%cts(2)%w, &
            model%pt(1:2)%alpha%fit, model%pt(1:2)%beta%fit, model%pt(1:2)%ar%fit, &
            model%pt(1:2)%ma%fit, model%pt(1:2)%d%fit, model%m, model%n, &
            model%pt(1:2)%npar, K, model%argsD)
    end if

    !---------------------------------------------------------------------------------
    ! calling the subroutine to fill the matrices with the calculated values
    !--------------------------------------------------------------------------------
    call return_model(model, length(1), ts, inf, extra, nd, D, T, E, h)
    !---------------------------------------------------------------------------------
    ! positive log-likelihood and score vector
    !---------------------------------------------------------------------------------
    sll =  -sll
    U =  -U

    !---------------------------------------------------------------------------------
    !  out-of-sample forecast
    !---------------------------------------------------------------------------------
    if (length(2) > 0) call forecast_model(model, length(2), xnew1, xnew2, forecast)
    return
  end subroutine final_model

  !-------------------------------------------------------------------------------------------------
  !
  !                  Subroutines used to simulate and predict
  !
  !-------------------------------------------------------------------------------------------------
  subroutine sim_model(rdist, n, burn, np, pdist, muonly, order, xregar, alpha, beta,  &
       phi, theta, d, link, lconfig, xreg1, xreg2, xstart, tstart, ts, rev)
    !---------------------------------------------------------------------------------------
    !
    !  Simulating a F-ARFIMA model
    !
    !---------------------------------------------------------------------------------------
    ! Input:
    !  rdist  : random number generator for the chosen distribution
    !  n      : final sample size
    !  burn   : the bur-in size
    !  np     : number of parameters in pdist
    !  pdist  : nuisance parameters to the conditional distribution
    !  muonly : indicates if the distribution only has a mu parameters
    !  order  : a 2 by 4 matrix with (nreg, p, q, inf)
    !  xregar : indicates if xreg is to be included in the AR recursion
    !  alpha  : a vector of size 2 with the intercept for each part
    !  beta   : a matrix with 2 rows (regressors' coefficients)
    !  phi    : a matrix with 2 rows (AR coefficients)
    !  theta  : a matrix with 2 rows (MA coefficients)
    !  d      : a vector of size 2 with the long memory parameter
    !  link   : a vector of size 8, indicating the links to be used in the model
    !           g1 (dummy), g11, g12, g13, g2, g21, g22, g23
    !  lconfig: a matrix of size 8 by 4 with
    !            - the lower and upper bounds related to the link functions
    !            - the constant and power arguments for polynomial/SIP links
    !  xreg1  : matrix of regressors corresponding to part 1
    !  xreg2  : matrix of regressors corresponding to part 2
    !  xstart : matrix of maxval(nreg) rows by 2 rows with starting values for the regressors
    !           (to be used when t < 1)
    !  tstart : starting value for y, g(nu) and g23(e1t)
    !  ts     : the generated time series yt, mut, eta1t, e1t, nut, gnut, eta2t, e2t
    !  rev    : flag to return error code
    !---------------------------------------------------------------------------------------
    ! August, 2023
    !  - adapted the code to work with part 1 and 2. (added new variables)
    !  - changed d == 0.0d0 to abs(d) < epsilon(1.d0)
    !  - removed ns, seed and rngtype
    !  - removed ylower and yupper
    !  - added lconfig (limits for the link)
    !  - fixed a typo in the inicialization of g_12(y) and xreg * beta.
    !    The correct check is p > 0 (now) instead of p > 1 (then)
    !
    ! Last revision: February, 2025
    !  - replaced nreg, p, q and inf by the matrix order
    !  - added pointers
    use distrib
    implicit none
    integer,  intent(in)    :: n, burn, np
    integer,  intent(in)    :: link(8), xregar(2)
    integer,  intent(inout) :: rev
    integer,  target, intent(inout) :: order(2,4)
    real(dp), target, intent(inout) :: ts(n + burn, 8)
    real(dp), target, intent(in)    :: tstart(3)
    integer,  pointer :: nreg(:), p(:), q(:), inf(:)
    real(dp), pointer :: yt(:), mu(:), eta1(:), error(:)
    real(dp), pointer :: nu(:), gnu(:), eta2(:), error2(:)
    real(dp), pointer :: ystart, e2start, vtstart
    real(dp), intent(inout) :: lconfig(8, 4)
    real(dp), intent(in)    :: pdist(np), alpha(2), d(2)
    real(dp), intent(in)    :: xstart(2, maxval([1, order(:,1)]))
    real(dp), intent(inout) :: beta(2, maxval([1, order(:,1)]))
    real(dp), intent(inout) :: phi(2, maxval([1, order(:,2)]))
    real(dp), intent(inout) :: theta(2, maxval([1, order(:,3)]))
    real(dp), intent(inout) :: xreg1(n + burn, max(1, order(1,1)))
    real(dp), intent(inout) :: xreg2(n + burn, max(1, order(2,1)))
    type(argsDist), intent(inout) :: rdist
    logical,  intent(in)    :: muonly
    type(argslink) :: argsLg(8)
    integer  ::  t, i, ierr, nlag, i_y, i_mu, i_nu, i_gnu
    real(dp) :: g12y(n + burn), g22(n + burn), e2(-order(2,4):(n + burn))
    real(dp) :: xb1(n + burn), xb2(n + burn)
    real(dp) :: xb1temp, xb2temp, g12ytemp, g22temp, e2temp, g11ytemp
    type(vetor), target :: vc(2)
    real(dp), pointer :: vc1(:), vc2(:)
    logical :: xar1, xar2, upx1, upx2, e20

    !--------------------------------------------------------
    ! at the end, if no revision is required, rev is set to 0
    !--------------------------------------------------------
    rev = 0
    nreg   => order(:, 1)
    p      => order(:, 2)
    q      => order(:, 3)
    inf    => order(:, 4)
    yt     => ts(:, 1)
    mu     => ts(:, 2)
    eta1   => ts(:, 3)
    error  => ts(:, 4)
    nu     => ts(:, 5)
    gnu    => ts(:, 6)
    eta2   => ts(:, 7)
    error2 => ts(:, 8)
    ystart  => tstart(1)
    vtstart => tstart(2)
    e2start => tstart(3)
    i_mu  = 2 ! g11(mut)
    i_y   = 3 ! g12(y)
    i_nu  = 5 ! g2(nut)
    i_gnu = 6 ! g12(gnut)
    upx1 = nreg(1) > 0
    upx2 = nreg(2) > 0

    !--------------------------------------------------------
    ! setting the link parameters
    !--------------------------------------------------------
    call set_link(link, lconfig, argsLg)

    !--------------------------------------------------------
    ! checking the order of the polynomials and initializing
    !--------------------------------------------------------
    do i = 1, 2
       !  If d = 0 uses inf = q
       inf(i) = max(inf(i), q(i))
       if (abs(d(i)) < epsmch) inf(i) = q(i)

       ! check for regressors
       if (nreg(i) == 0) then
          if (i == 1) xreg1 = 0.d0
          if (i == 2) xreg2 = 0.d0
          beta(i, :) = 0.d0
       end if

       ! check for AR component
       if (p(i) == 0) phi(i, :) = 0.d0

       ! check for MA component
       if (q(i) == 0) theta(i, :) = 0.d0
    end do

    do i = 1, 2
       ! vector of coefficients for the infinite sum
       call safe_allocate(vc(i)%z, 0, inf(i))
       if (q(i) == 0) then
          call vc_f(d(i), q(i), [1.d0], inf(i), vc(i)%z)
       else
          call vc_f(d(i), q(i), [1.d0, theta(i, 1:q(i))], inf(i), vc(i)%z)
       end if
    end do
    vc1 => vc(1)%z
    vc2 => vc(2)%z

    !--------------------------------------------------------
    ! initializing variables in the model and setting
    ! starting values
    !--------------------------------------------------------

    ! part 1
    yt    = 0.d0
    g12y  = 0.d0
    eta1  = 0.d0
    mu    = 0.d0
    error = 0.d0
    xb1   = 0.d0

    ! part 2
    eta2 = 0.d0
    nu   = 0.d0
    gnu  = 0.d0
    g22  = 0.d0
    e2   = e2start
    xb2  = 0.d0
    e20  = (abs(e2start) < epsmch) ! e2start = 0?

    ! auxiliar variables
    g12ytemp = 0.d0 ! g12(ystart)
    g22temp  = 0.d0 ! g22(vtstart)
    xar1     = xregar(1) == 1
    xar2     = xregar(2) == 1
    e2temp   = 0.d0 ! g23(e1t)
    xb1temp  = 0.d0 ! x1 * beta1
    xb2temp  = 0.d0 ! x2 * beta2

    ! update starting values
    if (p(1) > 0) then
       call g_start_update(ystart, g12ytemp, argsLg(i_y), ierr, 1)        ! g12(y), t < 1
       if(ierr > 0) then
          rev = ierr
          return
       end if
       if (xar1) xb1temp = sum(xstart(1, 1:nreg(1)) * beta(1, 1:nreg(1))) ! x1'beta1, for t < 1
    end if

    if (p(2) > 0) then
       call g_start_update(vtstart, g22temp, argsLg(6), ierr, 2)          ! g22(g(nu)), t < 1
       if(ierr > 0) then
          rev = ierr
          return
       end if
       if (xar2) xb2temp = sum(xstart(2, 1:nreg(2)) * beta(2, 1:nreg(2))) ! x2'beta2, t < 1
    end if

    if (inf(2) > 0) e2temp = e2start ! g23(e1t), t < 1

    !--------------------------------------------------------
    ! Loop to calculate mu_t, nu_t and y_t
    ! eta1(t) = alpha1 + x1 * b1 + AR1 + MA1
    ! eta1(2) = alpha2 + x2 * b2 + AR2 + MA2
    !--------------------------------------------------------
    do t = 1, (n + burn)

       ! PART 1:
       ! compute: eta1 = alpha + x * beta
       eta1(t) = alpha(1)
       if(upx1) then
          xb1(t) = sum(xreg1(t, :) * beta(1, 1:nreg(1)))
          eta1(t) = eta1(t) + xb1(t)
       end if

       ! update: eta1 = alpha + x * beta + AR
       do i = 1, p(1)
          if (t - i > 0) then
             g12ytemp = g12y(t - i)                            ! update g12(y)
             if (xar1) xb1temp = xb1(t-i)                      ! update x * beta
          end if
          eta1(t) = eta1(t) + (g12ytemp - xb1temp) * phi(1, i) ! add phi* (g(y) - x * beta)
       end do

       ! final update: eta1 = alpha + x * beta + AR + MA
       do i = 1, min(t - 1, inf(1))
          eta1(t) =  eta1(t) + vc1(i) * error(t - i)  ! error(j) = 0, if j < 0
       end do

       ! compute: mu(t) = g11^{-1}(eta1(t)) and check the bounds
       mu(t) = linkinv(eta1(t), argsLg(i_mu))
       call make_shift(mu(t), lconfig(i_mu, 1), lconfig(i_mu, 2), 1, rev, .true.)
       if(rev > 0) then
          if(rev > 1) return
          eta1(t) = linkfun(mu(t), argsLg(i_mu)) ! update eta(t) = g11(mut)
       end if

       if(.not. muonly) then
          ! PART 2
          ! compute eta2 = alpha + x * b
          if(upx2) xb2(t) = sum(xreg2(t, :) * beta(2, 1:nreg(2)))
          eta2(t) = alpha(2) + xb2(t)

          ! update: eta2 = alpha + x * beta + AR
          do i = 1, p(2)
             if (t - i > 0) then
                g22temp = g22(t - i)                             ! update g_22(gnu)
                if (xar2) xb2temp = xb2(t - i)                   ! update x * beta
             end if
             eta2(t) = eta2(t) + (g22temp - xb2temp) * phi(2, i) ! add phi* (g(y) - x * beta)
          end do

          ! final update: eta1 = alpha + x * beta + AR + MA
          if (e20) then
             nlag = min(t-1, inf(2))  ! if e2start = 0, avoid extra calculation
          else
             nlag = inf(2)
          end if
          do i = 1, nlag
             eta2(t) =  eta2(t) + vc2(i) * e2(t - i)
          end do

          ! compute gnu(t) = g21^{-1}(eta2(t))
          ! compute nu(t) = g^{-1}(gnu(t)) and check bounds
          gnu(t) = linkinv(eta2(t), argsLg(i_gnu))
          nu(t) = linkinv(gnu(t), argsLg(i_nu))
          call make_shift(nu(t), lconfig(i_nu, 1), lconfig(i_nu, 2), 2, rev, .true.)
          if(rev > 0) then
             if(rev > 1) return
             gnu(t) = linkfun(nu(t), argsLg(i_nu))    ! update g(nu(t))
             eta2(t) = linkfun(gnu(t), argsLg(i_gnu)) ! update eta(t)
          end if

          ! update g22(gnu(t))
          if (argsLg(i_gnu)%update) then
             g22(t) = linkfun(gnu(t), argsLg(i_gnu+1))
          else
             g22(t) = eta2(t)
          end if
       end if

       ! SAMPLING STEP
       ! sample: y(t) ~ f(y, mu, nu, pdist)
       yt(t) = rdist%rdist(np + 2, [mu(t), nu(t), pdist])
       if (yt(t) <= lconfig(i_y, 1)) then
          yt(t) = lconfig(i_y, 1) + epsmch
       elseif(yt(t) >= lconfig(i_y, 2)) then
          yt(t) = lconfig(i_y, 2) - epsmch
       end if

       g12y(t) =  linkfun(yt(t), argsLg(i_y)) ! g12(y)

       if (argsLg(4)%link == 1 .and. argsLg(i_y)%update) then
          g11ytemp = linkfun(yt(t), argsLg(i_mu)) ! compute g11(y)
       else
          g11ytemp = g12y(t) ! no update need or g11 is not required
       end if

       ! compute the error term
       error(t) =  g_err1(yt(t), mu(t), g11ytemp, eta1(t), argsLg(4)%link)
       e2(t) = linkfun(error(t), argsLg(8))
       error2(t) = e2(t)
    end do

    rev = 0 ! no revision needed
    return
  end subroutine sim_model

end module base
