module barc
  use base
  implicit none

contains

  subroutine set_link_barc(link, lconfig, argsL)
    !---------------------------------------------------------------------------
    ! Initializes the arguments in the link
    !---------------------------------------------------------------------------
    ! Input
    !  link: a size 4 vector with the link codes for g11, g12, g13 and h
    !
    ! Input/Output
    !  lconfig: a matrix of size 4 by 4 with configurations for the links
    !
    ! Output
    !   argsL: a size 4 argslink type variable with information on the links
    !---------------------------------------------------------------------------
    ! February, 2025: added this subroutine
    implicit none
    integer,  intent(in)    :: link(4)
    real(dp), intent(inout) :: lconfig(4, 4)
    type(argslink), intent(inout) :: argsL(4)

    argsL(:)%link   = link
    argsL(:)%lower  = lconfig(:, 1)
    argsL(:)%upper  = lconfig(:, 2)
    argsL(:)%par(1) = lconfig(:, 3)
    argsL(:)%par(2) = lconfig(:, 4)

    ! check if g11 and g12 are the same
    call check_update(argsL(1), argsL(2))

    ! error term:
    if (link(3) == 0) then ! a < y, mu < b
       argsL(3)%lower =  lconfig(2, 1) - lconfig(2, 2) ! a - b
       argsL(3)%upper =  lconfig(2, 2) - lconfig(2, 1) ! b - a
    else ! -Inf < g(y), g(mu) < Inf
       argsL(3)%lower =  -Infinity
       argsL(3)%upper =  Infinity
    end if
    return
  end subroutine set_link_barc

  function map_T(x, r, theta, mtype) result(Tx)
    implicit none
    integer,  intent(in) :: r, mtype
    real(dp), intent(in) :: x
    real(dp), intent(in) :: theta(r)
    real(dp) :: Tx

    Tx = 0.d0
    select case(mtype)
    case (1)
       ! (kx)(mod 1), k must be an integer greater or equal than 2
       Tx = theta(1) * x
       Tx = Tx - int(Tx)
    case (2)
       ! Rafael's map. 0 <= theta <= 1
       if (x < theta(1)) then
          Tx = x / theta(1)
       else
          Tx = theta(1) * (x - theta(1)) / (1.d0 - theta(1))
       end if
    case (3)
       ! logistic map. 0 <= theta  <= 4
       Tx = theta(1) * x * (1.d0 - x)
    case (4)
       ! Manneville-Pomeau. 0 < theta < 1
       Tx = x + x**(1 + theta(1))
       Tx = Tx - int(Tx)
    case (5)
       ! Lasota-Mackey's map. No theta
       if (x <= 0.5d0) then
          Tx = x / (1.d0 - x)
       else
          Tx = 2.d0 * x - 1.d0
       end if
    end select

  end function map_T

  subroutine start_par_barc(par, model)
    !--------------------------------------------------------
    !
    !             Parameter Initialization
    !
    ! Uses the subroutine from base module to allocate
    ! alpha, beta and phi. Allocates thetaT and u0
    !
    !---------------------------------------------------------
    ! Last revision: May, 2024
    ! - nu is replaced by alpha(2)
    implicit none
    type(argsModel), intent(inout) :: model
    real(dp), intent(in) :: par(sum(model%pt%npar))
    integer  :: n1, n2

    ! setting the values of alpha, beta and phi
    call start_par(par, model, 1)

    ! setting the initial values theta_T and u0 (if needed)
    n1 = model%pt(1)%alpha%fit + model%pt(1)%beta%fit + model%pt(1)%ar%fit
    n2 = n1

    ! ThetaT
    n1 = n2 + 1
    n2 = n2 + model%pt(1)%thetaT%fit
    if (n2 >= n1) model%pt(1)%thetaT%cf = par(n1:n2)

    ! u0
    n1 = n2 + 1
    n2 = n2 + model%pt(1)%u0%fit
    if (n2 >= n1) model%pt(1)%u0%cf = par(n1:n2)

    ! setting nu
    n1 = n2 + 1
    n2 = n2 + model%pt(2)%alpha%fit
    if (n2 >= n1) model%pt(2)%alpha%cf = par(n1:n2)

    return
  end subroutine start_par_barc

  subroutine mu_calc_barc(model)
    !------------------------------------------------------------------------
    !
    ! Recurrence for mu(t). Here nu(t) is fixed
    !
    ! Uses the subroutine from base module to calculate the recurrence
    ! assuming a MA(0) and d = 0. Here the chaotic part is added
    !
    !------------------------------------------------------------------------
    ! Last revision: May, 2024
    ! - nu is replaced by alpha(2)
    ! - the subrotine now also fills model%cts(2)%w
    implicit none
    type(argsModel), intent(inout) :: model
    real(dp) :: vc(1)
    integer  :: t, rev

    vc = 0.d0
    ! error scale makes no difference here:
    call mu_calc(model%n, model%y, model%cts(1)%g2start, model%cts(1)%gi1, model%cts(1)%gi2, &
         model%cts(1)%nreg, model%cts(1)%xreg, model%cts(1)%xstart, model%cts(1)%w, &
         model%cts(1)%eta, model%cts(1)%et, model%pt(1)%alpha%cf(1), model%pt(1)%beta%cf, &
         model%pt(1)%ar%length, model%pt(1)%ar%cf, model%cts(1)%xregar, 0, vc, &
         model%m, model%pt(1)%linkg(2:4))

    ! calculating T**{t - 1}(u0) and adding to eta_t
    model%cts(1)%orbit(1) = model%pt(1)%u0%cf(1)
    model%cts(1)%eta(1) = model%cts(1)%eta(1) + linkfun(model%cts(1)%orbit(1), model%pt(1)%linkh)
    model%cts(1)%w(1) = linkinv(model%cts(1)%eta(1), model%pt(1)%linkg(2))
    model%cts(1)%et(1) = g_err1(model%y(1), model%cts(1)%w(1), model%cts(1)%gi1(1), &
         model%cts(1)%eta(1), model%pt(1)%linkg(4)%link)

    ! to avoid the bounds:
    call make_shift(model%cts(1)%w(1), model%pt(1)%linkg(2)%lower, &
         model%pt(1)%linkg(2)%upper, 1, rev, .false.)
    if(rev > 0) model%cts(1)%eta(1) = linkfun(model%cts(1)%w(1), model%pt(1)%linkg(2))

    do t = 2, model%n
       ! T**{t - 1}(u0)
       model%cts(1)%orbit(t) = map_T(model%cts(1)%orbit(t - 1), model%pt(1)%thetaT%length, &
            model%pt(1)%thetaT%cf, model%pt(1)%map)
       ! eta + h(T**{t - 1}(u0))
       model%cts(1)%eta(t) = model%cts(1)%eta(t) + linkfun(model%cts(1)%orbit(t), model%pt(1)%linkh)
       ! mut = g**{-1}(eta_t)
       model%cts(1)%w(t) = linkinv(model%cts(1)%eta(t), model%pt(1)%linkg(2))
       ! to avoid the bounds:
       call make_shift(model%cts(1)%w(t), model%pt(1)%linkg(2)%lower, &
            model%pt(1)%linkg(2)%upper, 1, rev, .false.)
       if(rev > 0) model%cts(1)%eta(t) = linkfun(model%cts(1)%w(t), model%pt(1)%linkg(2))
       ! calculating rt
       model%cts(1)%et(t) = g_err1(model%y(t), model%cts(1)%w(t), model%cts(1)%gi1(t), &
            model%cts(1)%eta(t), model%pt(1)%linkg(4)%link)
    end do

    ! nu is fixed and g(nu) = nu = alpha(2)
    model%cts(2)%w = model%pt(2)%alpha%cf(1)

    return
  end subroutine mu_calc_barc

  subroutine mu_forecast_barc(model, nnew, xhat, forecast)
    !--------------------------------------------------------
    !
    !   Forecast for a BARC model
    !
    !---------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    integer,  intent(in)  :: nnew
    real(dp), intent(in)  :: xhat(nnew, max(1, model%pt(1)%beta%length))
    real(dp), target, intent(out) :: forecast(nnew, 3)
    real(dp), pointer :: muhat(:), etahat(:), That(:)
    real(dp) :: Ts(0:nnew)
    real(dp) :: XB((model%n + 1 - model%pt(1)%ar%length):(model%n + nnew))
    real(dp) :: gy((model%n + 1 - model%pt(1)%ar%length):(model%n + nnew))
    integer :: t, i

    forecast = 0.d0
    muhat  => forecast(:,1)
    etahat => forecast(:,2)
    That   => forecast(:,3)


    XB = 0.d0
    if (model%pt(1)%beta%length > 0) then
       do t = (model%n + 1 - model%pt(1)%ar%length), model%n
          XB(t) = sum(model%cts(1)%xreg(t, :) * model%pt(1)%beta%cf)
       end do
       do t = 1, nnew
          XB(model%n + t) = sum(xhat(t, :) * model%pt(1)%beta%cf)
       end do
    end if

    ! g12(y)
    if (model%pt(1)%ar%length > 0) gy((model%n + 1 - model%pt(1)%ar%length):model%n) = &
         model%cts(1)%gi2((model%n + 1 - model%pt(1)%ar%length):model%n)

    Ts(0) = model%cts(1)%orbit(model%n)
    do t = 1, nnew
       ! etahat = a + x * b + h(T^{t - 1}(u0))
       Ts(t) = map_T(Ts(t - 1), model%pt(1)%thetaT%length, model%pt(1)%thetaT%cf, model%pt(1)%map)
       etahat(t) =  model%pt(1)%alpha%cf(1) + XB(model%n + t) + linkfun(Ts(t), model%pt(1)%linkh)
       ! etahat = a + x * b + h(T^{t - 1}(u0)) + AR
       do i = 1, model%pt(1)%ar%length   ! if p < 1 ignores the loop
          etahat(t) = etahat(t) + model%pt(1)%ar%cf(i) * gy(model%n + t - i)
          if (model%cts(1)%xregar == 1) then
             etahat(t) = etahat(t) - model%pt(1)%ar%cf(i) * XB(model%n + t - i)
          end if
       end do
       ! muhat = g11^{-1}(etahat)
       muhat(t) = linkinv(etahat(t), model%pt(1)%linkg(2))
       ! g12(y) to be used in the AR recursion
       if (model%pt(1)%linkg(3)%update) then
          gy(model%n + t) = linkinv(muhat(t), model%pt(1)%linkg(3))
       else
          gy(model%n + t) = etahat(t)
       end if
    end do
    That = Ts(1:nnew)
    return
  end subroutine mu_forecast_barc

  subroutine get_model_barc(model, n, order, y, gy, xreg, ystart, xstart, &
       link, lconfig, map, npar, par, xregar,  nfix, alpha, flagsb, fvbeta, &
       flagsar, fvar, flagst, fvtheta, u0, extras, full, ierr)
    !-----------------------------------------------------------------------------------------
    !
    !  Subrotuine used to pass the values enter in the main program to the user defined
    !  variable that will be passed to the  generic subroutines
    !
    !-----------------------------------------------------------------------------------------
    implicit none
    type(argsModel), intent(inout) :: model
    integer,  intent(in)    :: order(2, 4), nfix(2, 5)
    integer,  intent(in)    :: n, link(9), map
    integer,  intent(in)    :: npar(2), xregar(2), extras(4)
    integer,  intent(in)    :: flagsb(2, maxval([1, nfix(:, 2)]))
    integer,  intent(in)    :: flagsar(2, maxval([1, nfix(:, 3)]))
    integer,  intent(in)    :: flagst(2, maxval([1, nfix(:, 4)]))
    integer,  intent(inout) :: ierr
    real(dp), intent(in)    :: y(n), ystart, u0
    real(dp), intent(in)    :: xstart(2, maxval([1, order(:, 1)]))
    real(dp), intent(in)    :: xreg(n, max(1, order(1, 1)))
    real(dp), intent(in)    :: par(sum(npar)), alpha(2)
    real(dp), intent(in)    :: fvbeta(2, maxval([1, nfix(:, 2)]))
    real(dp), intent(in)    :: fvar(2, maxval([1, nfix(:, 3)]))
    real(dp), intent(in)    :: fvtheta(2, maxval([1, nfix(:, 4)]))
    real(dp), intent(inout) :: lconfig(9, 4)
    real(dp), intent(out)   :: gy(n, 2)
    logical,  intent(in)    :: full
    ! dummy
    real(dp) :: xreg2(n, max(1, order(2, 1)))
    integer :: nfixx(2, 5), oorder(2, 4)

    if(full) then
       call model%argsD%init(current_models(1))  ! Set the model
    end if

    oorder = order
    oorder(:,3:4) = 0 ! ma order and inf
    nfixx = nfix
    nfixx(:, 5) = 1   ! set d as fixed and equals to zero.

    ! set some dummy arguments and call the main function
    call get_model(model, n, oorder, y, gy, xreg, xreg2, [ystart, 0.d0, 0.d0], xstart, &
         link(1:8), lconfig(1:8,:), 1, npar, par, xregar, nfixx, alpha, flagsb, fvbeta, &
         flagsar, fvar, flagst, [0.d0, 0.d0], [0.d0, 0.d0], extras, full, ierr)

    ! chaotic map
    model%pt(1)%map = map

    ! setting the link parameters for h(T**t(u0))
    model%pt(1)%linkh%link = link(9)
    model%pt(1)%linkh%lower = lconfig(9, 1)
    model%pt(1)%linkh%upper = lconfig(9, 2)
    model%pt(1)%linkh%par = lconfig(9, 3:4)

    ! allocating the orbit
    call safe_allocate(model%cts(1)%orbit, n)

    ! theta_T and u0
    call allocate_parvec(model%pt(1)%thetaT, order(1,3), nfix(1,4), flagst, fvtheta)
    call allocate_parvec(model%pt(1)%u0, 1, nfix(1,5), [1], [u0])
    model%pt(1)%npar = model%pt(1)%npar + order(1,3) + 1 - sum(nfix(1,4:5))
    if (model%pt(1)%npar > 0) model%pt(1)%skip = 0

    !  Setting the parameter's fixed values
    call start_par_barc(par, model)

    ! for now, skipp allocating score-vector, information matrix and related matrices
    return
  end subroutine get_model_barc


  subroutine return_model_barc(model, ts)
    implicit none
    type(argsModel), intent(in) :: model
    real(dp), intent(out) :: ts(model%n, 4)

    ts(:,1) = model%cts(1)%w
    ts(:,2) = model%cts(1)%eta
    ts(:,3) = model%cts(1)%et(1:model%n)
    ts(:,4) = model%cts(1)%orbit
    return
  end subroutine return_model_barc

  subroutine U_barc_numeric(model, npar, par, U)
    implicit none
    type(argsModel), intent(inout) :: model
    integer,  intent(in)  :: npar
    real(dp), intent(in)  :: par(npar)
    real(dp), intent(out) :: U(npar)
    integer :: i
    real(dp), parameter :: eps = 1.d-4
    real(dp) :: par1(npar), par2(npar), f1, f2

    ! Score vector - numerical derivative
    U(:) = 0.d0
    do i = 1, npar
       par1 = par
       par2 = par
       par1(i) = par1(i) + eps
       par2(i) = par2(i) - eps
       call start_par_barc(par1, model)
       call mu_calc_barc(model)
       f1 = model%argsD%llk_dist(model%m, model%n, model%y, model%cts(1)%w, model%cts(2)%w)
       call start_par_barc(par2, model)
       call mu_calc_barc(model)
       f2 = model%argsD%llk_dist(model%m, model%n, model%y, model%cts(1)%w, model%cts(2)%w)
       U(i) = (f1 - f2) / (2. * eps)
    end do

    ! Restoring parameters
    call start_par_barc(par, model)
    return
  end subroutine U_barc_numeric

  subroutine K_barc_numeric(model, npar, par, K)
    implicit none
    type(argsModel), intent(inout) :: model
    integer,  intent(in)  :: npar
    real(dp), intent(in)  :: par(npar)
    real(dp), intent(out) :: K(npar, npar)
    real(dp), parameter   :: eps = 1.d-4
    real(dp) :: np1(npar), np2(npar), np3(npar), np4(npar), h1, h2, h3, h4
    integer  :: i, j


    ! Numerical Information matrix
    do i = 1, npar
       do j = 1, i
          np1 = par; np2 = par; np3 = par; np4 = par
          np1(i) = np1(i) + eps; np1(j) = np1(j) + eps
          np2(i) = np2(i) + eps; np2(j) = np2(j) - eps
          np3(i) = np3(i) - eps; np3(j) = np3(j) + eps
          np4(i) = np4(i) - eps; np4(j) = np4(j) - eps
          call start_par_barc(np1, model)
          call mu_calc_barc(model)
          h1 = model%argsD%llk_dist(model%m, model%n, model%y, model%cts(1)%w, model%cts(2)%w)
          call start_par_barc(np2, model)
          call mu_calc_barc(model)
          h2 = model%argsD%llk_dist(model%m, model%n, model%y, model%cts(1)%w, model%cts(2)%w)
          call start_par_barc(np3, model)
          call mu_calc_barc(model)
          h3 = model%argsD%llk_dist(model%m, model%n, model%y, model%cts(1)%w, model%cts(2)%w)
          call start_par_barc(np4, model)
          call mu_calc_barc(model)
          h4 = model%argsD%llk_dist(model%m, model%n, model%y, model%cts(1)%w, model%cts(2)%w)
          K(i, j) = (h1  - h2 - h3  + h4 ) / (4. * eps * eps)
          K(j, i) = K(i, j)
       end do
    end do
    ! Information is -Hess.
    K =  -K

    ! Restoring parameters
    call start_par_barc(par, model)
    return
  end subroutine K_barc_numeric


  subroutine loglik_barc(loglik, model, npar, par, sll, U)
    !------------------------------------------------------------------
    !
    !   Log - likelihood: BARC model
    !
    !------------------------------------------------------------------
    implicit none
    class(optimFunc), intent(inout) :: loglik
    type(argsModel), intent(inout)  :: model
    integer,  intent(in)  :: npar
    real(dp), intent(in)  :: par(npar)
    real(dp), intent(out) :: sll, U(npar)

    loglik%dummy = .true.

    ! Initializing parameters
    call start_par_barc(par, model)

    U = 0.d0
    if (model%sco == 1) then
       ! Score vector
       !  here we use a different order of calculation because of the
       !  numerical derivative
       call U_barc_numeric(model, npar, par, U)
       U =  - U
    end if

    ! Calculating recursions that do not depend on the distribution
    call mu_calc_barc(model)

    ! log - likelihood for BARC model
    sll = model%argsD%llk_dist(model%m, model%n, model%y, model%cts(1)%w, model%cts(2)%w)
    sll =  -sll
    if (sll <  - Huge(1.d0)) sll =  - Huge(1.d0)
    if (sll > Huge(1.d0)) sll = Huge(1.d0)

    return
  end subroutine loglik_barc

  subroutine loglik_barc_nelder(loglik, model, npar, par, sll)
    !------------------------------------------------------------------
    !
    !   Subroutine to be used in Nelder-Mead optimization subroutine
    !
    !------------------------------------------------------------------
    implicit none
    class(optimFunc), intent(inout)   :: loglik
    type(argsModel), intent(inout) :: model
    integer,  intent(in)  :: npar
    real(dp), intent(in)  :: par(npar)
    real(dp), intent(out) :: sll
    real(dp) :: par_aux(npar), U(npar)

    ! Back to original scale
    par_aux = par
    call  transform_par(par_aux, npar, model%bounds, .true.)

    model%llk = 1
    model%sco = 0
    call loglik_barc(loglik, model, npar, par_aux, sll, U)

    return
  end subroutine loglik_barc_nelder

end module barc

subroutine simbarc(length, order, ts, xreg, ystart, xstart, link, lconfig, map, &
     xregar, alpha, beta, phi, theta, u0,  rev)
  use barc
  !---------------------------------------------------------------------
  !
  !  Simulating a BARC model
  !
  !---------------------------------------------------------------------
  implicit none
  integer, target  :: length(2), order(3)
  integer, pointer :: n, burn
  integer  :: map, link(4), xregar, rev
  real(dp), target  :: ts(sum(length), 5)
  real(dp), pointer :: yt(:), mu(:), eta(:), error(:), Tt(:)
  real(dp) :: xreg(sum(length), max(1, order(1)))
  real(dp) :: ystart, xstart(max(1, order(1)))
  real(dp) :: lconfig(4, 4)
  real(dp) :: alpha(2), beta(max(1, order(1))), u0
  real(dp) :: phi(max(1, order(2))), theta(max(1, order(3)))
  real(dp) :: g12y(sum(length)), xb(sum(length))
  integer  ::  t, i, ierr
  integer, pointer :: nreg, p, r
  type(argslink) :: argsL(4)
  type(argsDist) :: rdist
  real(dp) :: g11ytemp, g12ytemp, xbtemp, nu
  logical  :: xar, upx

  ! Initialize the argsModel object
  call rdist%init(current_models(1))  ! Set the model

  ! revision required
  rev = 1

  nreg => order(1)
  p    => order(2)
  r    => order(3)
  n    => length(1)
  burn => length(2)
  yt    => ts(:,1)
  mu    => ts(:,2)
  eta   => ts(:,3)
  error => ts(:,4)
  Tt    => ts(:,5)
  xar = xregar == 1
  upx = nreg > 0
  nu = alpha(2)

  ! setting the link parameters
  call set_link_barc(link, lconfig, argsL)

  ! check for regressors
  if (nreg == 0) then
     xreg = 0.d0
     beta = 0.d0
  end if

  ! check for AR component
  if (p == 0) phi = 0.d0

  ! initializing variables in the model
  yt    = 0.d0
  g12y  = 0.d0
  xb    = 0.d0
  mu    = 0.d0
  eta   = 0.d0
  error = 0.d0
  Tt    = 0.d0

  ! auxiliar variables
  g12ytemp = 0.d0 ! g12(ystart)
  xbtemp   = 0.d0   ! x * beta

  if (p > 1) then
     ! g12(y), t < 1
     call g_start_update(ystart, g12ytemp, argsL(2), ierr, 1)
     if(ierr > 0) then
        rev = ierr
        return
     end if
     if (xar) xbtemp = sum(xstart * beta) ! x'beta, for t < 1
  end if

  ! t = 1
  Tt(1) = u0
  eta(1) = alpha(1) + sum(xreg(1, :) * beta) + &
       (g12ytemp - xbtemp) * sum(phi) + linkfun(Tt(1), argsL(4))
  ! compute: mu(t) = g11^{-1}(eta(t)) and check the bounds
  mu(1) = linkinv(eta(1), argsL(1))
  call make_shift(mu(1), lconfig(1, 1), lconfig(1, 2), 1, rev, .true.)
  if(rev > 0) then
     if(rev > 1) return
     eta(1) = linkfun(mu(1), argsL(1)) ! update eta(t) = g12(mut)
  end if
  yt(1) = rdist%rdist(2, [mu(1), nu])
  g12y(1) = linkfun(yt(1), argsL(2))
  if (argsL(3)%link == 0) then
     error(1) = yt(1) - mu(1)
  else
     if (argsL(2)%update) then
        g11ytemp = linkfun(yt(1), argsL(1)) ! compute g11(y)
     else
        g11ytemp = g12y(1) ! no update need
     end if
     error(1) = g11ytemp - eta(1)
  end if

  ! t > 1
  do t = 2, (n + burn)
     ! T^{t - 1}(u0)
     Tt(t) =  map_T(Tt(t - 1), r, theta, map)

     ! eta(t) = alpha + x * b + AR + h(T^{t - 1}(u0))
     eta(t) = alpha(1)
     if(upx) then
        xb(t) = sum(xreg(t, :) * beta)
        eta(t) = eta(t) + xb(t)
     end if

     do i = 1, p
        ! updating the auxiliar variables g12ytemp and xb
        if (t - i > 0) then
           ! update g12(y)
           g12ytemp = g12y(t - i)
           ! check if xreg is used in AR recursion
           if (xar) xbtemp = xb(t-i)
        end if
        eta(t) = eta(t) + (g12ytemp - xbtemp) * phi(i)
     end do
     eta(t)  = eta(t) +  linkfun(Tt(t), argsL(4))

     ! mu(t) = g^{-1}(eta(t))
     mu(t) = linkinv(eta(t), argsL(1))
     call make_shift(mu(t), lconfig(1, 1), lconfig(1, 2), 1, rev, .true.)
     if(rev > 0) then
        if(rev > 1) return
        eta(t) = linkfun(mu(t), argsL(1)) ! update eta(t) = g11(mut)
     end if

     ! y(t) ~ beta
     yt(t) = rdist%rdist(2, [mu(t), nu])
     if (yt(t) <= lconfig(2, 1)) then
        yt(t) = lconfig(2, 1) + epsmch
     elseif(yt(t) >= lconfig(2, 2)) then
        yt(t) = lconfig(2, 2) - epsmch
     end if
     g12y(t) =  linkfun(yt(t), argsL(2))

     if (argsL(3)%link == 0) then
        error(t) = yt(t) - mu(t)
     else
        if (argsL(2)%update) then
           g11ytemp = linkfun(yt(t), argsL(1)) ! compute g11(y)
        else
           g11ytemp = g12y(t) ! no update need
        end if
        error(t) = g11ytemp - eta(t)
     end if
  end do

  rev = 0
  return
end subroutine simbarc

subroutine extractbarc(length, order, ts, xreg, ystart, xstart, &
     xnew, forecast, link, lconfig, map, npar, par, xregar, nfix, alpha, &
     flagsb, beta, flagsar, phi, flagst, theta, u0, extras, sll, U, K, ierr)
  use barc
  implicit none
  ! time series
  integer  :: ierr, map
  integer  :: length(2), link(9), order(2, 4)
  real(dp) :: lconfig(9, 4)
  real(dp) :: ts(length(1), 7)
  real(dp) :: xreg(length(1), max(1, order(1, 1)))
  real(dp) :: ystart, xstart(2, maxval([1, order(:, 1)]))
  real(dp) :: xnew(max(1, length(2)), max(1, order(1, 1)))
  real(dp) :: forecast(max(1, length(2)), 3)
  ! parameters
  integer  :: npar(2), xregar(2), nfix(2, 5)
  integer  :: flagsb(2, maxval([1, nfix(:, 2)]))
  integer  :: flagsar(2, maxval([1, nfix(:, 3)]))
  integer  :: flagst(2, maxval([1, nfix(:, 4)]))
  real(dp) :: par(sum(npar)), alpha(2), u0
  real(dp) :: beta(2, maxval([1, nfix(:, 2)]))
  real(dp) :: phi(2, maxval([1, nfix(:, 3)]))
  real(dp) :: theta(2, maxval([1, nfix(:, 4)]))
  ! likelihood
  integer  :: extras(5)
  real(dp) :: sll, U(max(1, sum(npar) * extras(3)))
  real(dp) :: K(max(1, sum(npar) * extras(4)), max(1, sum(npar) * extras(4)))
  ! auxiliar
  type(argsModel) :: model
  type(optimFunc) :: loglik
  real(dp) :: Utemp(sum(npar))

  ! allocating matrices and vectors and setting variables fixed values
  call get_model_barc(model, length(1), order, ts(:,1), ts(:,2:3), xreg, ystart, xstart, &
       link, lconfig, map, npar, par, xregar,  nfix, alpha, flagsb, beta, &
       flagsar, phi, flagst, theta, u0, extras(1:4), .true., ierr)

  ! changed the order because we are using the numerical hessian
  K = 0.d0
  if (extras(4) == 1) then
     call K_barc_numeric(model, sum(npar), par, K)
  end if

  sll = 0.d0
  U = 0.d0
  ! calculates:
  !    mu, eta, Tt
  !     - sll: log - likelihood (if llk = 1)
  !     - U: score vector (if sco = 1)
  call loglik_barc(loglik, model, sum(npar), par, sll, Utemp)

  call return_model_barc(model, ts(:,4:7))
  sll =  -sll
  if (extras(3) == 1) U =  -Utemp

  if (length(2) > 0) then
     call mu_forecast_barc(model, length(2), xnew, forecast)
  end if
  return
end subroutine extractbarc

subroutine optimbarc(method, length, order, ts, xreg, ystart, xstart, xnew, forecast, &
     link, lconfig, map, npar, par, nbd, bounds, xregar, nfix, alpha, flagsb, beta, flagsar, &
     phi, flagst, theta, u0, extras, sll, U, K, cf1, nc2, cf2, neval, conv)
  use barc
  implicit none
  ! time series
  integer  :: method, length(2), link(9), order(2, 4), map
  real(dp) :: lconfig(9, 4)
  real(dp) :: ts(length(1), 7)
  real(dp) :: xreg(length(1), max(1, order(1, 1)))
  real(dp) :: ystart, xstart(2, maxval([1, order(:,1)]))
  real(dp) :: xnew(max(1, length(2)), max(1, order(1, 1)))
  real(dp) :: forecast(max(1, length(2)), 3)
  ! parameters
  integer  :: npar(2), nfix(2,5), xregar(2)
  integer  :: flagsb(2, maxval([1, nfix(:, 2)]))
  integer  :: flagsar(2, maxval([1, nfix(:, 3)]))
  integer  :: flagst(2, maxval([1, nfix(:, 4)]))
  real(dp) :: par(sum(npar)), alpha(2), u0
  real(dp) :: beta(2, maxval([1, nfix(:, 2)]))
  real(dp) :: phi(2, maxval([1, nfix(:, 3)]))
  real(dp) :: theta(2, maxval([1, nfix(:, 4)]))
  ! optimization
  !--------------------------------------------
  ! iprint < 0 = no print
  ! stopcr = stopping critereon  (1.d-4)
  !---------------------------------------------
  integer  :: extras(5), nbd(sum(npar))
  integer  :: cf1(2), neval, conv, nc2
  real(dp) :: bounds(sum(npar), 2), cf2(nc2)
  real(dp) :: sll, U(max(1, sum(npar) * extras(3)))
  real(dp) :: K(max(1, sum(npar) * extras(4)), max(1, sum(npar) * extras(4)))
  !
  type(argsModel) :: model
  type(optimFunc) :: loglik
  type(optimizer) :: optim
  integer  :: extr(4), cf(3)
  real(dp) :: Utemp(sum(npar))

  select case (method)
  case(0)
     ! L-BFGS-B
     extr = [0, 1, 1, 0]  ! m, llk, sco, info
     cf = [cf1, 1]  ! cf(3) is a multiplier for the size of the score vector
     loglik%loglik => loglik_barc
     optim%optim   => optim_lbfgsb
  case(1)
     ! Nelder-Mead
     extr = [0, 1, 0, 0]  ! m, llk, sco, info
     cf = [cf1, 0] ! cf(3) is a multiplier for the size of the score vector
     loglik%functn => loglik_barc_nelder
     optim%optim   => optim_nelder
  end select

  ! allocating matrices and vectors and setting variables fixed values
  call get_model_barc(model, length(1), order, ts(:,1), ts(:,2:3), xreg, ystart, xstart, &
       link, lconfig, map, npar, par, xregar, nfix, alpha, flagsb, beta, &
       flagsar, phi, flagst, theta, u0, extr, .true., conv)

  ! Optimization subroutine
  call optim%optim(loglik, model, sum(npar), par, nbd, bounds, sll, U, cf, nc2, cf2, neval, conv)

  K = 0.d0
  if (extras(4) > 0) then
     call K_barc_numeric(model, sum(npar), par, K)
  end if

  ! Reports the final model
  model%sco = extras(3)
  call loglik_barc(loglik, model, sum(npar), par, sll, Utemp)

  call return_model_barc(model, ts(:,4:7))
  sll =  -sll
  U = 0.d0
  if (extras(3) == 1) U =  -Utemp

  if (length(2) > 0) then
     call mu_forecast_barc(model, length(2), xnew, forecast)
  end if

  return
end subroutine optimbarc

subroutine predictbarc(length, order, ts, xreg, xnew, forecast, link, lconfig,  map, &
     npar, par, xregar, nfix, alpha, flagsb, fvbeta, flagsar, fvar, flagst, fvtheta, u0)
  use barc
  implicit none
  !-----------------------------------------------------------------------------------------
  !
  !  Subrotuine used for prediction.
  !  The values of y, g12y, eta, error,
  !
  !-----------------------------------------------------------------------------------------
  integer  :: length(2), link(9), order(2,4), nfix(2,5)
  integer  :: npar(2), xregar(2), map
  integer  :: flagsb(2, maxval([1, nfix(:, 2)]))
  integer  :: flagsar(2, maxval([1, nfix(:, 3)]))
  integer  :: flagst(2, maxval([1, nfix(:, 4)]))
  real(dp) :: ts(length(1), 5)
  real(dp) :: xreg(length(1), max(1, order(1, 1)))
  real(dp) :: xnew(max(1, length(2)), max(1, order(1, 1)))
  real(dp) :: forecast(max(1, length(2)), 3)
  real(dp) :: lconfig(9, 4)
  real(dp) :: par(sum(npar)), alpha(2), u0
  real(dp) :: fvbeta(2, maxval([1, nfix(:, 2)]))
  real(dp) :: fvar(2, maxval([1, nfix(:, 3)]))
  real(dp) :: fvtheta(2, maxval([1, nfix(:, 4)]))
  type(argsModel) :: model
  real(dp) :: xstart(2, maxval([1, order(:,1)]))
  integer :: ierr

  xstart = 0.d0     ! dummy

  ! allocating matrices and vectors and setting variables fixed values
  call get_model_barc(model, length(1), order, ts(:,1), ts(:,2:3), xreg, 0.d0, xstart, &
       link, lconfig, map, npar, par, xregar, nfix, alpha, flagsb, fvbeta, &
       flagsar, fvar, flagst, fvtheta, u0, [0,0,0,0], .false., ierr)

  ! setting the values of eta, error and orbit (from 1 to n)
  model%cts(1)%et = ts(:, 4)
  model%cts(1)%eta = 0.d0
  model%cts(1)%orbit = ts(:, 5)

  !  predicted values
  call mu_forecast_barc(model, length(2), xnew, forecast)

  return
end subroutine predictbarc
